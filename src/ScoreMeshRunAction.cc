#include "ddfastsim/ScoreMeshRunAction.h"
#include "ddfastsim/ScoreMeshEventAction.h"
#include "G4AnalysisManager.hh"
#include "G4EventManager.hh"
#include "DD4hep/InstanceCount.h"
#include "G4Run.hh"

/// Standard constructor with initializing arguments
dd4hep::sim::ScoreMeshRunAction::ScoreMeshRunAction(dd4hep::sim::Geant4Context* c, const std::string& n):
        dd4hep::sim::Geant4RunAction(c, n) {
  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");

  // Default filename, can be overriden with /analysis/setFileName
  analysisManager->SetFileName("Output");
  dd4hep::InstanceCount::increment(this);
  
      declareProperty("numOfRhoCells", cellNumRho = 18 );
      declareProperty("numOfZCells", cellNumZ = 45 );
      declareProperty("sizeOfRhoCells", cellSizeRho = 2.325 );
      declareProperty("sizeOfPhiCells", cellSizePhi = 2 * CLHEP::pi / 50. );
      declareProperty("sizeOfZCells", cellSizeZ = 5.05 );
  }
        /// Default destructor
        dd4hep::sim::ScoreMeshRunAction::~ScoreMeshRunAction() {
  dd4hep::InstanceCount::decrement(this);
}
/// begin-of-run callback
void dd4hep::sim::ScoreMeshRunAction::begin(const G4Run*) {
           // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Create directories
  analysisManager->SetVerboseLevel(0);
  // Default max value of energy stored in histogram (in GeV)
  G4double maxEnergy = 1000;
  cellNumPhi = 2 * CLHEP::pi / cellSizePhi;

  // Creating control histograms
  analysisManager->CreateH1("energyParticle", "Primary energy;E_{MC} (GeV);Entries", 1024, 0,
                            1.1 * maxEnergy);
  analysisManager->CreateH1("energyDepositedInVirtual", "Deposited energy;E_{MC} (GeV);Entries",
                            1024, 0, 1.1 * maxEnergy);
  analysisManager->CreateH1(
    "energyRatioInVirtual", "Ratio of energy deposited to primary;E_{dep} /  E_{MC};Entries",
    1024, 0, 1);
  analysisManager->CreateH1("time", "Simulation time; time (s);Entries", 2048, 0, 100);
  analysisManager->CreateH1("longProfile", "Longitudinal profile;t (mm);#LTE#GT (MeV)", cellNumZ,
                            -0.5 * cellSizeZ, (cellNumZ - 0.5) * cellSizeZ);
  analysisManager->CreateH1("transProfile", "Transverse profile;r (mm);#LTE#GT (MeV)", cellNumRho,
                            -0.5 * cellSizeRho, (cellNumRho - 0.5) * cellSizeRho);
  analysisManager->CreateH1("longFirstMoment",
                            "First moment of longitudinal distribution;#LT#lambda#GT (mm);Entries",
                            1024, -0.5 * cellSizeZ,
                            cellNumZ * cellSizeZ / 2);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1("transFirstMoment",
                            "First moment of transverse distribution;#LTr#GT "
                            "(mm);Entries",
                            1024, -0.5 * cellSizeRho,
                            cellNumRho * cellSizeRho /
                              1);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1(
    "longSecondMoment",
    "Second moment of longitudinal distribution;#LT#lambda^{2}#GT "
    "(mm^{2});Entries",
    1024, 0, std::pow(cellNumZ * cellSizeZ, 2) / 25);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1(
    "transSecondMoment", "Second moment of transverse distribution;#LTr^{2}#GT (mm^{2});Entries",
    1024, 0, std::pow(cellNumRho * cellSizeRho, 2) / 5);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1("phiProfile",
                            "Azimuthal angle profile, centred at mean;phi;#LTE#GT (MeV)",
                            cellNumPhi, - (cellNumPhi - 0.5) * cellSizePhi,
                            (cellNumPhi - 0.5) * cellSizePhi);
  analysisManager->CreateH1("numHitsInVirtual", "Number of hits above 0.5 keV", 4048, 0, 20000);
  analysisManager->CreateH1("cellEnergy", "Cell energy distribution;log10(E/MeV);Entries",
                            1024, -4, 2);

  // Creating ntuple
  analysisManager->CreateNtuple("global", "Event data");
  analysisManager->CreateNtupleDColumn("EnergyMC");
  analysisManager->CreateNtupleDColumn("SimTime");
  analysisManager->FinishNtuple();

  ScoreMeshEventAction* fEventAction = dynamic_cast<ScoreMeshEventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());
  analysisManager->CreateNtuple("virtualReadout", "Cylindrical mesh readout");
  analysisManager->CreateNtupleDColumn("EnergyCell", fEventAction->GetCalEdep());
  analysisManager->CreateNtupleIColumn("rhoCell", fEventAction->GetCalRho());
  analysisManager->CreateNtupleIColumn("phiCell", fEventAction->GetCalPhi());
  analysisManager->CreateNtupleIColumn("zCell", fEventAction->GetCalZ());
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("run", "Run data");
  analysisManager->CreateNtupleDColumn("N");
  analysisManager->CreateNtupleDColumn("RunTime");
  analysisManager->FinishNtuple();

  analysisManager->OpenFile();
  fTimer.Start();
        }
/// End-of-run callback
void dd4hep::sim::ScoreMeshRunAction::end(const G4Run* aRun) {
  fTimer.Stop();
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(2, 0, aRun->GetNumberOfEvent());
  analysisManager->FillNtupleDColumn(2, 1, fTimer.GetRealElapsed());
  analysisManager->AddNtupleRow(2);
  analysisManager->Write();
  analysisManager->CloseFile();}

#include "DDG4/Factories.h"
DECLARE_GEANT4ACTION(ScoreMeshRunAction)
