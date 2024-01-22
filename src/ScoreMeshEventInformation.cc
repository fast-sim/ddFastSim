#include "ddfastsim/ScoreMeshEventInformation.h"
#include <CLHEP/Vector/ThreeVector.h>  // for operator<<
#include <G4VUserEventInformation.hh>  // for G4VUserEventInformation
#include <G4ios.hh>                    // for G4cout, G4endl
#include <ostream>                     // for operator<<, basic_ostream, ost...
#include "G4AnalysisManager.hh"          // for G4AnalysisManager
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScoreMeshEventInformation::ScoreMeshEventInformation()
  : G4VUserEventInformation()
  , fDirection()
  , fPosition()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScoreMeshEventInformation::~ScoreMeshEventInformation() = default;


void ScoreMeshEventInformation::Print() const
{
  std::cout << "Event information\nPrimary particle direction = " << fDirection
         << "\nPrimary particle position = " << fPosition << std::endl;
}

