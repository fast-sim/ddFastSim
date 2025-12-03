// DD4hep includes
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/DetFactoryHelper.h"

// todo: remove gaudi logging and properly capture output
#define endmsg std::endl
#define lLog std::cout
namespace MSG {
const std::string ERROR = " Error: ";
const std::string DEBUG = " Debug: ";
const std::string INFO  = " Info: ";
}


using dd4hep::Volume;
using dd4hep::DetElement;
using dd4hep::xml::Dimension;
using dd4hep::PlacedVolume;
using dd4hep::Position;

namespace det {

static dd4hep::Ref_t createPlanes(dd4hep::Detector& lcdd, xml_h xmlElement,
                                          dd4hep::SensitiveDetector sensDet) {

  xml_det_t xmlDet = xmlElement;
  std::string detName = xmlDet.nameStr();
  // Make DetElement
  DetElement eCal(detName, xmlDet.id());

  // Make volume that envelopes the whole barrel; set material to air
  Dimension dimensions(xmlDet.dimensions());
  dd4hep::Box envelopeShape(dimensions.x()/2., dimensions.y()/2., dimensions.z()/2.);
  Volume envelopeVolume(detName, envelopeShape, lcdd.air());
  // Invisibility seems to be broken in visualisation tags, have to hardcode that
  envelopeVolume.setVisAttributes(lcdd, dimensions.visStr());

  xml_comp_t layer = xmlElement.child(_Unicode(layer));
  int numLayers = layer.repeat();
  int numSlicesPerLayer = layer.numChildren(_Unicode(slice));
  double thicknessLayer = 0;
  double sliceOneThickness = 0;
  int layerOneFlag = 0;
  for (xml_coll_t x_slice(layer, _Unicode(slice)); x_slice; ++x_slice) {
    xml_comp_t xSlice = x_slice;
    thicknessLayer += xSlice.thickness();
	    if (layerOneFlag++ == 0){
		sliceOneThickness += xSlice.thickness();
	    }
    layerOneFlag++;
  }
  lLog << MSG::INFO << " Creating " << numLayers << " layers with " << numSlicesPerLayer << "slices each" <<
    " and thickness of " << thicknessLayer << " cm." << endmsg;

  // set the sensitive detector type to the DD4hep calorimeter
  dd4hep::SensitiveDetector sd = sensDet;
  dd4hep::xml::Dimension sdType = xmlDet.child(_U(sensitive));
  sd.setType(sdType.typeStr());

  // loop over layers
  double layerFacePosition = -dimensions.z()/2.+sliceOneThickness; /// Box is placed at origin
  for (int iLayer = 0; iLayer < numLayers; iLayer++) {
    int sliceNumber = 0;
    for (xml_coll_t x_slice(layer, _Unicode(slice)); x_slice; ++x_slice) {
      xml_comp_t xSlice = x_slice;
      std::string sliceMaterial = xSlice.materialStr();
      double sliceThickness = xSlice.thickness();
      dd4hep::Box layerShape(dimensions.x()/2., dimensions.y()/2., sliceThickness/2.); /// double check this placement is correct!
      lLog << MSG::INFO << " LayerShape.x(): " << layerShape.x() << " " << " LayerShape.y(): " << layerShape.y() << " " << " LayerShape.z(): " <<  layerShape.z() <<   endmsg;
      Volume layerVol(sliceMaterial+dd4hep::_toString(iLayer, "_layer%d")+dd4hep::_toString(sliceNumber, "_slice%d"), layerShape, lcdd.material(sliceMaterial));
      PlacedVolume placedLayer = envelopeVolume.placeVolume(layerVol, 1, Position(0., 0., layerFacePosition + sliceThickness/2.));
      lLog  << MSG::INFO << " Placing slice " << sliceNumber << " made of " << sliceMaterial << " within layer " << iLayer << " at layer face posiion " << layerFacePosition << " of thickness " << sliceThickness << endmsg;
      //placedLayer.addPhysVolID("layer", i);
      DetElement caloLayer(eCal, sliceMaterial+"_layer" +std::to_string(iLayer) + "_slice" + std::to_string(sliceNumber), numSlicesPerLayer * iLayer + sliceNumber);
      placedLayer.addPhysVolID("layer", iLayer);
      placedLayer.addPhysVolID("slice", sliceNumber);
      caloLayer.setPlacement(placedLayer);
      if (xSlice.isSensitive())
        layerVol.setSensitiveDetector(sensDet);
      sliceNumber++;
      layerFacePosition += sliceThickness;
    }
  }

  // Place envelope (or barrel) volume
  Volume motherVol = lcdd.pickMotherVolume(eCal);
  PlacedVolume placedECal = motherVol.placeVolume(envelopeVolume, 1, Position(0., 0., -dimensions.z()/2.));
  placedECal.addPhysVolID("system", eCal.id());
  eCal.setPlacement(placedECal);
  return eCal;
}
}  // namespace det

DECLARE_DETELEMENT(SandwichPlanes, det::createPlanes)
