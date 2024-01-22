#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Geant4EventAction.h"
#include "ddfastsim/ScoreMeshEventInformation.h"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "DDG4/Geant4FastSimSpot.h"


/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim   {

    /// Class to implement the standard sensitive detector for scoring mesh calorimeters
    struct Geant4ShowerMeshCalorimeter : public Geant4Calorimeter{
      /// Number of mesh readout cells in cylindrical coordinates
      int numOfRhoCells = 18;
      int numOfZCells = 45;
      /// Size of mesh readout cells in cylindrical coordinates.
      double sizeOfRhoCells = 2.325;
      double sizeOfPhiCells = 2 * CLHEP::pi / 50.;
      double sizeOfZCells = 5.05;
      /// Entrance data
      G4ThreeVector  entrancePosition = G4ThreeVector(-1,-1,-1);
      G4ThreeVector  entranceDirection = G4ThreeVector(-1,-1,-1);
      Geant4ShowerMeshCalorimeter() : Geant4Calorimeter()
				//	, _firstLayerNumber(1)
      {}
    };
    /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<Geant4ShowerMeshCalorimeter>::defineCollections() {
      m_collectionID = defineCollection<Geant4Calorimeter::Hit>(m_sensitive.readout().name());
    }
    /// template specialization for c'tor in order to define property: FirstLayerNumber
    template <> 
    Geant4SensitiveAction<Geant4ShowerMeshCalorimeter>::Geant4SensitiveAction(Geant4Context* ctxt,
										const std::string& nam,
										DetElement det,
										Detector& lcdd_ref)
      : Geant4Sensitive(ctxt,nam,det,lcdd_ref), m_collectionID(0)
    {
      initialize();
      defineCollections();
      InstanceCount::increment(this);
      declareProperty("numOfRhoCells", m_userData.numOfRhoCells = 18 );
      declareProperty("numOfZCells", m_userData.numOfZCells = 45 );
      declareProperty("sizeOfRhoCells", m_userData.sizeOfRhoCells = 2.325 );
      declareProperty("sizeOfPhiCells", m_userData.sizeOfPhiCells = 2 * CLHEP::pi / 50. );
      declareProperty("sizeOfZCells", m_userData.sizeOfZCells = 5.05 );
    }
    /// Method for generating hit(s) using the information of G4Step object.
    template <> bool Geant4SensitiveAction<Geant4ShowerMeshCalorimeter>::process(const G4Step* aStep,G4TouchableHistory*) {
      G4double edep = aStep->GetTotalEnergyDeposit();
      if(edep == 0.)
        return true;

      if(m_userData.entrancePosition.x() == -1) {
        auto info = dynamic_cast<ScoreMeshEventInformation*>(
              G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation());
        if(info == nullptr)
          return true;
        m_userData.entrancePosition  = info->GetPosition();
        m_userData.entranceDirection = info->GetDirection();
        std::cout << "Size of rho SD cell = " << m_userData.sizeOfRhoCells << std::endl;
        std::cout << "Size of z SD cell = " << m_userData.sizeOfZCells << std::endl;
        std::cout << "Entrance position = " << m_userData.entrancePosition << std::endl;
        std::cout << "Entrance direction = " << m_userData.entranceDirection<< std::endl;
      }

      auto delta = aStep->GetPostStepPoint()->GetPosition() - m_userData.entrancePosition;

      // Calculate rotation matrix along the particle momentum direction
      // It will rotate the shower axes to match the incoming particle direction
      G4RotationMatrix rotMatrix = G4RotationMatrix();
      double particleTheta       = m_userData.entranceDirection.theta();
      double particlePhi         = m_userData.entranceDirection.phi();
      rotMatrix.rotateZ(-particlePhi);
      rotMatrix.rotateY(-particleTheta);
      G4RotationMatrix rotMatrixInv = CLHEP::inverseOf(rotMatrix);

      delta = rotMatrix * delta;

      int rhoNo = std::floor(delta.perp() / m_userData.sizeOfRhoCells);
      int phiNo = std::floor((CLHEP::pi + delta.phi()) / m_userData.sizeOfPhiCells);
      int zNo   = std::floor(delta.z() / m_userData.sizeOfZCells);
      std::size_t hitID =
        m_userData.numOfRhoCells * m_userData.numOfZCells * phiNo + m_userData.numOfZCells * rhoNo + zNo;
      
      if(zNo >= m_userData.numOfZCells || rhoNo >= m_userData.numOfRhoCells || zNo < 0)
      {
        return true;
      }

      dd4hep::sim::Geant4Calorimeter::Hit* hit = collection(m_collectionID)->findByKey<dd4hep::sim::Geant4Calorimeter::Hit>(hitID);
      if (hit==nullptr) {
        hit = new dd4hep::sim::Geant4Calorimeter::Hit();
        hit->cellID = hitID;
        hit->position.SetXYZ(rhoNo, phiNo, zNo);
        collection(m_collectionID)->add(hitID, hit);      
    }

      // Add energy deposit from G4Step
      hit->energyDeposit += edep;
      return true;
    }

    template <> bool
    Geant4SensitiveAction<Geant4ShowerMeshCalorimeter>::processFastSim(const Geant4FastSimSpot* spot, G4TouchableHistory* /* hist */)
    {
      G4double edep = spot->energy();
      if(edep == 0.)
        return true;

      if(m_userData.entrancePosition.x() == -1) {
        auto info = dynamic_cast<ScoreMeshEventInformation*>(
              G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation());
        if(info == nullptr)
          return true;
        m_userData.entrancePosition  = info->GetPosition();
        m_userData.entranceDirection = info->GetDirection();
      }

      auto delta = spot->hitPosition() - m_userData.entrancePosition;

      // Calculate rotation matrix along the particle momentum direction
      // It will rotate the shower axes to match the incoming particle direction
      G4RotationMatrix rotMatrix = G4RotationMatrix();
      double particleTheta       = m_userData.entranceDirection.theta();
      double particlePhi         = m_userData.entranceDirection.phi();
      rotMatrix.rotateZ(-particlePhi);
      rotMatrix.rotateY(-particleTheta);
      G4RotationMatrix rotMatrixInv = CLHEP::inverseOf(rotMatrix);

      delta = rotMatrix * delta;

      int rhoNo = std::floor(delta.perp() / m_userData.sizeOfRhoCells);
      int phiNo = std::floor((CLHEP::pi + delta.phi()) / m_userData.sizeOfPhiCells);
      int zNo   = std::floor(delta.z() / m_userData.sizeOfZCells);

      std::size_t hitID =
        m_userData.numOfRhoCells * m_userData.numOfZCells * phiNo + m_userData.numOfZCells * rhoNo + zNo;
      
      if(zNo >= m_userData.numOfZCells || rhoNo >= m_userData.numOfRhoCells || zNo < 0)
      {
        return true;
      }

      dd4hep::sim::Geant4Calorimeter::Hit* hit = collection(m_collectionID)->findByKey<dd4hep::sim::Geant4Calorimeter::Hit>(hitID);
      if (hit==nullptr) {
        hit = new dd4hep::sim::Geant4Calorimeter::Hit();
        hit->cellID = hitID;
        hit->position.SetXYZ(rhoNo, phiNo, zNo);
        collection(m_collectionID)->add(hitID, hit);      
    }

      // Add energy deposit from G4Step
      hit->energyDeposit += edep;
      return true;
    }

    typedef Geant4SensitiveAction<Geant4ShowerMeshCalorimeter> Geant4ShowerMeshCalorimeterAction;

  }}
#include "DDG4/Factories.h"
DECLARE_GEANT4SENSITIVE(Geant4ShowerMeshCalorimeterAction)
