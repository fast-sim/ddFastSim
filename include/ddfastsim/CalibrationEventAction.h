#ifndef CALIBRATIONEVENTACTION_HH
#define CALIBRATIONEVENTACTION_HH

#include <G4Types.hh>            // for G4int, G4double
#include <vector>                // for vector
#include "G4UserEventAction.hh"  // for G4UserEventAction
class G4Event;
#include "DDG4/Geant4Handle.h"
#include "DDG4/Geant4Kernel.h"
#include "DDG4/Geant4EventAction.h"

namespace dd4hep{
  namespace sim {
    class CalibrationEventAction: public dd4hep::sim::Geant4EventAction{
      public:
        /// Standard constructor with initializing arguments
        CalibrationEventAction(dd4hep::sim::Geant4Context* c, const std::string& n);
        /// Default destructor
        virtual ~CalibrationEventAction();
        /// begin-of-event callback
        inline virtual void begin(const G4Event*) override;
        /// End-of-event callback
        virtual void end(const G4Event*) override;
        /// begin-of-run callback
        void beginRun(const G4Run*);
        /// End-of-run callback
        void endRun(const G4Run*);
      private:
        /// ID of a hit collection to analyse
        G4int fHitCollectionEMID = -1;
        G4int fHitCollectionHadID = -1;
        /// Slice ID of the sensitive layer
        int fSensitiveEMSlice = 0;
        int fSensitiveHadSlice = 0;
      };
}
}
#endif /* EVENTACTION_HH */
