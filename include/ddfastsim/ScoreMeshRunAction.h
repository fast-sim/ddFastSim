#ifndef SCOREMESHRUNACTION_HH
#define SCOREMESHRUNACTION_HH

#include <CLHEP/Units/SystemOfUnits.h>       // for GeV
#include <G4String.hh>                       // for G4String
#include <G4ThreeVector.hh>                  // for G4ThreeVector
#include <G4Types.hh>                        // for G4int
#include "G4Timer.hh"            // for G4Timer
class G4ParticleDefinition;
#include "DDG4/Geant4Handle.h"
#include "DDG4/Geant4Kernel.h"
#include "DDG4/Geant4RunAction.h"
/**
 * @brief Run action
 *
 * Create analysis file and define control histograms.
 */
namespace dd4hep{
  namespace sim {
class ScoreMeshRunAction: public dd4hep::sim::Geant4RunAction {
      public:
        ScoreMeshRunAction() = delete;
        /// Standard constructor with initializing arguments
        ScoreMeshRunAction(dd4hep::sim::Geant4Context* c, const std::string& n) ;
        /// Default destructor
        virtual ~ScoreMeshRunAction() ;
        /// begin-of-run callback
        void begin(const G4Run*) override;
        /// End-of-run callback
        void end(const G4Run*) override;
        /// begin-of-event callback
        void beginEvent(const G4Event*) ;
        /// End-of-event callback
        void endEvent(const G4Event*) ;
        private:
        /// Timer measurement from Geant4
        G4Timer fTimer;
          /// Guard shared analysis output ownership across multiple instances
          bool fOwnsAnalysisOutput = false;

  // Get detector dimensions
  G4int cellNumZ       = 45;
  G4int cellNumRho     = 18;
  G4int cellNumPhi     = 50;
  G4double cellSizeZ   = 5.05;
  G4double cellSizeRho = 2.325;
  G4double cellSizePhi = 2*CLHEP::pi/50;
      };
  }
}
#endif /* RUNACTION_HH */
