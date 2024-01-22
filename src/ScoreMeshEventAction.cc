#include "ddfastsim/ScoreMeshEventAction.h"
#include "ddfastsim/ScoreMeshEventInformation.h"
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

dd4hep::sim::ScoreMeshEventAction::ScoreMeshEventAction(dd4hep::sim::Geant4Context* c, const std::string& n):
    dd4hep::sim::Geant4EventAction(c,n)
  {
    dd4hep::InstanceCount::increment(this);
    declareProperty("numOfRhoCells", fCellNbRho = 18 );
    declareProperty("numOfZCells", fCellNbZ = 45 );
    declareProperty("sizeOfRhoCells", fCellSizeRho = 2.325 );
    declareProperty("sizeOfPhiCells", fCellSizePhi = 2 * CLHEP::pi / 50. );
    declareProperty("sizeOfZCells", fCellSizeZ = 5.05 );
  }

        /// Default destructor
        dd4hep::sim::ScoreMeshEventAction::~ScoreMeshEventAction() {
  dd4hep::InstanceCount::decrement(this);
        }

void dd4hep::sim::ScoreMeshEventAction::begin(const G4Event* aEvent)
{
  fTimer.Start();
  fCellNbPhi  = 2 * CLHEP::pi / fCellSizePhi;
  fCalEdep.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalRho.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalPhi.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalZ.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
}

void dd4hep::sim::ScoreMeshEventAction::end(const G4Event* aEvent)
{
  fTimer.Stop();
  // Get hits collection ID (only once)
  if(fHitCollectionID == -1)
  {
    fHitCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("ECalBarrelCollection");
  }
  // Get hits collection
  auto hitsCollection =
    static_cast<dd4hep::sim::Geant4HitCollection*>(aEvent->GetHCofThisEvent()->GetHC(fHitCollectionID));
  
  if(hitsCollection == nullptr)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << fHitCollectionID;
    G4Exception("ScoreMeshEventAction::GetHitsCollection()", "MyCode0001", FatalException, msg);
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

  // Resize back to initial mesh size
  fCalEdep.resize(fCellNbRho * fCellNbPhi * fCellNbZ,0);
  fCalRho.resize(fCellNbRho * fCellNbPhi * fCellNbZ,0);
  fCalPhi.resize(fCellNbRho * fCellNbPhi * fCellNbZ,0);
  fCalZ.resize(fCellNbRho * fCellNbPhi * fCellNbZ,0);

  std::cout << "Event action: Num of phi x rho x z cells: " << fCellNbPhi << " x "
            << fCellNbRho << " x " << fCellNbZ << std::endl;
  std::cout << "Number of hits in a collection: " <<  hitsCollection->GetSize() << std::endl;

  // Fill histograms
  dd4hep::sim::Geant4Calorimeter::Hit* hit                  = nullptr;
  G4double hitEn                 = 0;
  G4double totalEnergy           = 0;
  G4int hitZ                     = -1;
  G4int hitRho                   = -1;
  G4int hitPhi                   = -1;
  G4int hitID                   = 0;
  G4int numNonZeroThresholdCells = 0;
  G4double tDistance = 0., rDistance = 0., phiDistance = 0.;
  G4double tFirstMoment = 0., tSecondMoment = 0.;
  G4double rFirstMoment = 0., rSecondMoment = 0.;
  G4double phiMean = 0.;
  for(size_t iHit = 0; iHit < hitsCollection->GetSize(); iHit++)
  {
    hit     = static_cast<dd4hep::sim::Geant4Calorimeter::Hit*>(hitsCollection->hit(iHit));
    hitZ    = hit->position.z();
    hitRho  = hit->position.x();
    hitPhi  = hit->position.y();
    hitEn   = hit->energyDeposit;
    hitID = hit->cellID;

    if (hitPhi != (hitID / (fCellNbRho *fCellNbZ))) {
    G4cout << "ERROR, problem with  PHI IDs: " << hitPhi << " != " <<  hitID / (fCellNbRho *fCellNbZ) << G4endl;
    //return;
  }
    if (hitRho != (hitID % (fCellNbRho *fCellNbZ)/fCellNbZ)) {
    G4cout << "ERROR, problem with  Rho IDs: " << hitRho << " != " <<  hitID % (fCellNbRho *fCellNbZ) /fCellNbZ << G4endl;
    //return;
  }
  if (hitZ != (hitID % (fCellNbRho *fCellNbZ) %fCellNbZ)) {
    G4cout << "ERROR, problem with  Z IDs: " << hitZ << " != " <<  hitID % (fCellNbRho *fCellNbZ) %fCellNbZ << G4endl;
    //return;
  }
    
    // std::size_t hitID =
    //       fCellNbRho *fCellNbZ * phiNo
    //     +fCellNbZ * rhoNo
    //     + zNo;
    // hitID / (fCellNbRho *fCellNbZ) = phi
    // hitID % (fCellNbRho *fCellNbZ) /fCellNbZ = rho
    // hitID % (fCellNbRho *fCellNbZ) %fCellNbZ = z
    // hitID = 2 * 3 * 1   +    3 * 1   +   1 = 6+3+1 = 10
    //        hitID / 6 = 1. hitID %6 / 3 = 1  hitID %6 %3 = 1
    // hitID = 2 * 3 * 0   +    3 * 1   +   2 = 3+2 = 5
    //        hitID / 6 = 0 hitID %6 / 3 = 1  hitID %6 %3  = 2
    // hitID = 2 * 3 * 1   +    3 * 0   +   2 = 6+2 = 8
    //        hitID / 6 = 1. hitID %6 / 3 = 0  hitID %6 %3 = 2
    

    if(hitEn > 0)
    {
      totalEnergy += hitEn;
      tDistance = hitZ * fCellSizeZ;
      rDistance = hitRho * fCellSizeRho;
      phiDistance = hitPhi * fCellSizePhi;
      tFirstMoment += hitEn * tDistance;
      rFirstMoment += hitEn * rDistance;
      phiMean += hitEn * phiDistance;
      analysisManager->FillH1(4, tDistance, hitEn);
      analysisManager->FillH1(5, rDistance, hitEn);
      if(hitEn > 0.0005)
      {  // e > 0.5 keV
        fCalEdep[numNonZeroThresholdCells] = hitEn;
        fCalRho[numNonZeroThresholdCells]  = hitRho;
        fCalPhi[numNonZeroThresholdCells]  = hitPhi;
        fCalZ[numNonZeroThresholdCells]    = hitZ;
        numNonZeroThresholdCells++;
        analysisManager->FillH1(12, std::log10(hitEn));
      }
    }
  }
  tFirstMoment /= totalEnergy;
  rFirstMoment /= totalEnergy;
  phiMean /= totalEnergy;
  analysisManager->FillH1(0, primaryEnergy / CLHEP::GeV);
  analysisManager->FillH1(1, totalEnergy / CLHEP::GeV);
  analysisManager->FillH1(2, totalEnergy / primaryEnergy);
  analysisManager->FillH1(3, fTimer.GetRealElapsed());
  analysisManager->FillH1(6, tFirstMoment);
  analysisManager->FillH1(7, rFirstMoment);
  analysisManager->FillH1(11, numNonZeroThresholdCells);
  // Resize to store only energy hits above threshold
  fCalEdep.resize(numNonZeroThresholdCells);
  fCalRho.resize(numNonZeroThresholdCells);
  fCalPhi.resize(numNonZeroThresholdCells);
  fCalZ.resize(numNonZeroThresholdCells);
  analysisManager->FillNtupleDColumn(0, 0, primaryEnergy);
  analysisManager->FillNtupleDColumn(0, 1, fTimer.GetRealElapsed());
  // Second loop over hits to calculate second moments
  for(size_t iHit = 0; iHit < hitsCollection->GetSize(); iHit++)
  {
    hit    = static_cast<dd4hep::sim::Geant4Calorimeter::Hit*>(hitsCollection->hit(iHit));
    hitEn  = hit->energyDeposit;
    hitZ   = hit->position.z();
    hitRho = hit->position.x();
    hitPhi = hit->position.y();
    if(hitEn > 0)
    {
      tDistance = hitZ * fCellSizeZ;
      rDistance = hitRho * fCellSizeRho;
      phiDistance = hitPhi * fCellSizePhi;
      tSecondMoment += hitEn * std::pow(tDistance - tFirstMoment, 2);
      rSecondMoment += hitEn * std::pow(rDistance - rFirstMoment, 2);
      analysisManager->FillH1(10, phiDistance - phiMean, hitEn);
    }
  }
  tSecondMoment /= totalEnergy;
  rSecondMoment /= totalEnergy;
  analysisManager->FillH1(8, tSecondMoment);
  analysisManager->FillH1(9, rSecondMoment);

  analysisManager->AddNtupleRow(0);
  //if(numNonZeroThresholdCells > 0) analysisManager->AddNtupleRow(1);
}

#include "DDG4/Factories.h"
DECLARE_GEANT4ACTION(ScoreMeshEventAction)
