#include "ddfastsim/CalibrationEventAction.h"
#include "G4AnalysisManager.hh"
#include "DD4hep/InstanceCount.h"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ThreeVector.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SDManager.hh"                // for G4SDManager
#include "DDG4/Geant4Data.h"
#include "DDG4/Geant4HitCollection.h"
#include "DDSegmentation/BitFieldCoder.h"

dd4hep::sim::CalibrationEventAction::CalibrationEventAction(dd4hep::sim::Geant4Context* c, const std::string& n):
    dd4hep::sim::Geant4EventAction(c,n)
  {
    dd4hep::InstanceCount::increment(this);
    declareProperty("sensitiveEMSlice", fSensitiveEMSlice = 0 );
    declareProperty("sensitiveHadSlice", fSensitiveHadSlice = 0 );
  }

        /// Default destructor
        dd4hep::sim::CalibrationEventAction::~CalibrationEventAction() {
  dd4hep::InstanceCount::decrement(this);
        }

void dd4hep::sim::CalibrationEventAction::begin(const G4Event* aEvent)
{}

void dd4hep::sim::CalibrationEventAction::end(const G4Event* aEvent)
{
  // Get hits collection ID (only once)
  if(fHitCollectionEMID == -1)
  {
    fHitCollectionEMID = G4SDManager::GetSDMpointer()->GetCollectionID("ECalBarrelCollection");
  }
  const BitFieldCoder decoderEM("system:8,barrel:3,module:4,layer:6,slice:5,x:32:-16,y:-16" );
  if(fHitCollectionHadID == -1)
  {
    fHitCollectionHadID = G4SDManager::GetSDMpointer()->GetCollectionID("HCalBarrelCollection");
  }
  const BitFieldCoder decoderHad("system:8,barrel:3,module:4,layer:6,slice:5,x:32:-16,y:-16" );
  // Get hits collection
  auto hitsCollectionEM =
    static_cast<dd4hep::sim::Geant4HitCollection*>(aEvent->GetHCofThisEvent()->GetHC(fHitCollectionEMID));
  auto hitsCollectionHad =
    static_cast<dd4hep::sim::Geant4HitCollection*>(aEvent->GetHCofThisEvent()->GetHC(fHitCollectionHadID));
  
  if(hitsCollectionEM == nullptr)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << fHitCollectionEMID;
    G4Exception("CalibrationEventAction::GetHitsCollection()", "MyCode0001", FatalException, msg);
  }
  if(hitsCollectionHad == nullptr)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << fHitCollectionHadID;
    G4Exception("CalibrationEventAction::GetHitsCollection()", "MyCode0002", FatalException, msg);
  }
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Retrieve information from primary vertex and primary particle
  // To calculate shower axis and entry point to the detector
  auto primaryVertex =
    G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetPrimaryVertex();
  auto primaryParticle   = primaryVertex->GetPrimary(0);
  G4double primaryEnergy = primaryParticle->GetTotalEnergy();
  // Estimate from vertex and particle direction the entry point to the detector
  // Calculate entrance point to the detector located at z = 0
  auto primaryDirection = primaryParticle->GetMomentumDirection();
  auto primaryEntrance =
    primaryVertex->GetPosition() - primaryVertex->GetPosition().z() * primaryDirection;

  // Fill histograms
  dd4hep::sim::Geant4Calorimeter::Hit* hit                  = nullptr;
  G4double hitEn                 = 0;
  G4double totalEnergyEM         = 0;
  G4double totalEnergyHad        = 0;
  G4double activeEnergyEM         = 0;
  G4double activeEnergyHad        = 0;
  G4int hitSlice                 = -1;
  G4int hitID                   = 0;
  G4int restEM                   = 0;
  G4int restHad                  = 0;
  std::cout << "Num of hits in EM: " << hitsCollectionEM->GetSize() << " , had: " << hitsCollectionHad->GetSize() << std::endl;
  for(size_t iHit = 0; iHit < hitsCollectionEM->GetSize(); iHit++)
  {
    hit     = static_cast<dd4hep::sim::Geant4Calorimeter::Hit*>(hitsCollectionEM->hit(iHit));
    hitID = hit->cellID;
    hitSlice = decoderEM.get(hitID, "slice");
    hitEn  = hit->energyDeposit;

    if(hitSlice == fSensitiveEMSlice && hitEn > 0) {
      activeEnergyEM += hitEn;
    } else {restEM++;}
    totalEnergyEM += hitEn;
  }
  for(size_t iHit = 0; iHit < hitsCollectionHad->GetSize(); iHit++)
  {
    hit     = static_cast<dd4hep::sim::Geant4Calorimeter::Hit*>(hitsCollectionHad->hit(iHit));
    hitID = hit->cellID;
    hitSlice = decoderHad.get(hitID, "slice");
    hitEn  = hit->energyDeposit;

    if(hitSlice == fSensitiveHadSlice && hitEn > 0) {
      activeEnergyHad += hitEn;
    } else {restHad++;}
    totalEnergyHad += hitEn;
  }
  analysisManager->FillH1(0, primaryEnergy / CLHEP::GeV);
  analysisManager->FillH1(1, totalEnergyEM / CLHEP::GeV);
  analysisManager->FillH1(2, totalEnergyHad / CLHEP::GeV);
  analysisManager->FillH1(3, activeEnergyEM / totalEnergyEM);
  analysisManager->FillH1(4, activeEnergyHad / totalEnergyHad);
  analysisManager->FillH1(5, (activeEnergyEM + activeEnergyHad) / (totalEnergyEM+totalEnergyHad));
  std::cout << "Non-summed hits in EM: " << restEM << " , had: " << restHad << std::endl;
  std::cout << " energy EM: " << totalEnergyEM <<  " -> sf = " << activeEnergyEM/totalEnergyEM << std::endl;
  std::cout << " energy Had: " << totalEnergyHad << " -> sf = " << activeEnergyHad/totalEnergyHad << std::endl;
  std::cout << " energy sum: " << totalEnergyEM << " , " << totalEnergyHad << " -> " << (activeEnergyEM + activeEnergyHad) / (totalEnergyEM+totalEnergyHad) << std::endl;
}

#include "DDG4/Factories.h"
DECLARE_GEANT4ACTION(CalibrationEventAction)
