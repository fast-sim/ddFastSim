#include "DD4hep/DetFactoryHelper.h"

namespace det {
static dd4hep::Ref_t
createSimpleBox(dd4hep::Detector& lcdd, xml_h e, dd4hep::SensitiveDetector sensDet) {
  xml_det_t x_det = e;
  std::string name = x_det.nameStr();
  dd4hep::DetElement boxDet(name, x_det.id());

  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(boxDet);

  xml_comp_t boxDim(x_det.child(_U(dimensions)));

  // define 
  dd4hep::Box caloBox(boxDim.x()/2., boxDim.y()/2., boxDim.z()/2);

      //cylinderDim.rmin(), cylinderDim.rmax(), cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

  dd4hep::Volume boxVol(
      x_det.nameStr() + "_SimpleBox", caloBox, lcdd.material(boxDim.materialStr()));

  if (x_det.isSensitive()) {
    dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
    boxVol.setSensitiveDetector(sensDet);
    sensDet.setType(sdType.typeStr());
  }

    boxVol.setAttributes(lcdd, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  dd4hep::PlacedVolume boxPhys;

  //double zoff = boxDim.z_offset();
  //if (fabs(zoff) > 0.000000000001) {
  //  dd4hep::Position trans(0., 0., zoff);
  //  boxPhys = experimentalHall.placeVolume(boxVol,
  //                                              dd4hep::Transform3D(dd4hep::RotationZ(0.), trans));
  //} else
  boxPhys = experimentalHall.placeVolume(boxVol);

  boxPhys.addPhysVolID("system", x_det.id());

  boxDet.setPlacement(boxPhys);

  boxDet.setVisAttributes(lcdd, x_det.visStr(), boxVol);

  return boxDet;
}
}
DECLARE_DETELEMENT(SingleBox, det::createSimpleBox)
