#include "ddfastsim/CalibrationRunAction.h"
#include "ddfastsim/CalibrationEventAction.h"
#include "G4AnalysisManager.hh"
#include "G4EventManager.hh"
#include "DD4hep/InstanceCount.h"
#include "G4Run.hh"

/// Standard constructor with initializing arguments
dd4hep::sim::CalibrationRunAction::CalibrationRunAction(dd4hep::sim::Geant4Context* c, const std::string& n):
        dd4hep::sim::Geant4RunAction(c, n) {
  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");

  // Default filename, can be overriden with /analysis/setFileName
  analysisManager->SetFileName("Output");
  dd4hep::InstanceCount::increment(this);
  }
        /// Default destructor
        dd4hep::sim::CalibrationRunAction::~CalibrationRunAction() {
  dd4hep::InstanceCount::decrement(this);
}
/// begin-of-run callback
void dd4hep::sim::CalibrationRunAction::begin(const G4Run*) {
           // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Create directories
  analysisManager->SetVerboseLevel(0);
  // Default max value of energy stored in histogram (in GeV)
  G4double maxEnergy = 1000;

  // Creating control histograms
  analysisManager->CreateH1("energyParticle", "Primary energy;E_{MC} (GeV);Entries", 1024, 0,
                            1.1 * maxEnergy);
  analysisManager->CreateH1("energyDepositedInECal", "Deposited energy;E_{EM,MC} (GeV);Entries",
                            1024, 0, 1.1 * maxEnergy);
  analysisManager->CreateH1("energyDepositedInHCal", "Deposited energy;E_{Had,MC} (GeV);Entries",
                            1024, 0, 1.1 * maxEnergy);
  analysisManager->CreateH1(
    "energyRatioECal", "Ratio of energy deposited in ECal active/all;E_{dep} /  E_{MC};Entries",
    1024, 0, 1);
  analysisManager->CreateH1(
    "energyRatioHCal", "Ratio of energy deposited in HCal active/all;E_{dep} /  E_{MC};Entries",
    1024, 0, 1);
  analysisManager->CreateH1(
    "energyRatioBothCal", "Ratio of energy deposited in ECal and HCal active/all;E_{dep} /  E_{MC};Entries",
    1024, 0, 1);

  analysisManager->OpenFile();
        }
/// End-of-run callback
void dd4hep::sim::CalibrationRunAction::end(const G4Run* aRun) {
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

#include "DDG4/Factories.h"
DECLARE_GEANT4ACTION(CalibrationRunAction)
