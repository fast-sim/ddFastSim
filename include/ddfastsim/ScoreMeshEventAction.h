#ifndef SCOREMESHEVENTACTION_HH
#define SCOREMESHEVENTACTION_HH

#include <G4Types.hh>            // for G4int, G4double
#include <vector>                // for vector
#include "G4UserEventAction.hh"  // for G4UserEventAction
#include "G4Timer.hh"            // for G4Timer
class G4Event;
#include "DDG4/Geant4Handle.h"
#include "DDG4/Geant4Kernel.h"
#include "DDG4/Geant4EventAction.h"   
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

/**
       *  \author  M.Frank
       *  \version 1.0
       *  \ingroup DD4HEP_SIMULATION
       */
       namespace dd4hep{
        namespace sim {
      class ScoreMeshEventAction: public dd4hep::sim::Geant4EventAction{
      public:
        /// Standard constructor with initializing arguments
        ScoreMeshEventAction(dd4hep::sim::Geant4Context* c, const std::string& n);
        /// Default destructor
        virtual ~ScoreMeshEventAction();
        /// begin-of-event callback
        inline virtual void begin(const G4Event*) override;
        /// End-of-event callback
        virtual void end(const G4Event*) override;
        /// begin-of-run callback
        void beginRun(const G4Run*);
        /// End-of-run callback
        void endRun(const G4Run*);
        inline std::vector<G4double>& GetCalEdep() { return fCalEdep; }
        inline std::vector<G4int>& GetCalRho() { return fCalRho; }
        inline std::vector<G4int>& GetCalPhi() { return fCalPhi; }
        inline std::vector<G4int>& GetCalZ() { return fCalZ; }
        private:
        /// Timer measurement from Geant4
        G4Timer fTimer;
        /// ID of a hit collection to analyse
        G4int fHitCollectionID = -1;
        /// Size of cell along Z axis
        double fCellSizeZ =  5.05*mm;
        /// Size of cell along radius of cylinder
        double fCellSizeRho =  2.325*mm;
        /// Size of cell in azimuthal angle
        double fCellSizePhi = CLHEP::pi/25.;
        /// Number of readout cells along radius
        int fCellNbRho = 18;
        /// Number of readout cells in azimuthal angle
        int fCellNbPhi = 50;
        /// Number of readout cells along z axis
        int fCellNbZ = 45;
        /// Cell energy deposits to be stored in ntuple
        std::vector<double> fCalEdep;
        /// Cell ID of radius to be stored in ntuple
        std::vector<int> fCalRho;
        /// Cell ID of azimuthal angle to be stored in ntuple
        std::vector<int> fCalPhi;
        /// Cell ID of z axis to be stored in ntuple
        std::vector<int> fCalZ;
      };
}
}
#endif /* EVENTACTION_HH */
