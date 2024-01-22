#ifndef ScoreMeshEVENTINFORMATION_HH
#define ScoreMeshEVENTINFORMATION_HH

#include <G4Types.hh>                    // for G4bool
#include "G4ThreeVector.hh"              // for G4ThreeVector
#include "G4VUserEventInformation.hh"    // for G4VUserEventInformation
#include <array>

/**
 * @brief Event information
 *
 */

class ScoreMeshEventInformation : public G4VUserEventInformation
{
 public:
  ScoreMeshEventInformation();
  virtual ~ScoreMeshEventInformation();
  
  /// Set particle direction
  inline void SetDirection(const G4ThreeVector& aDirection) { fDirection = aDirection; };
  /// Get particle direction
  inline G4ThreeVector GetDirection() const { return fDirection; };
  /// Set particle position
  inline void SetPosition(const G4ThreeVector& aPosition) { fPosition = aPosition; };
  /// Get particle position
  inline G4ThreeVector GetPosition() const { return fPosition; };
  /// Set flag
  inline void SetFlag(G4bool aFlag) { fIfSet = aFlag; };
  /// Get flag
  inline G4bool GetFlag() const { return fIfSet; };
  /// Print
  void Print() const final;

 private:
  /// Particle direction. By default equal to the default particle gun direction.
  G4ThreeVector fDirection = { 0, 1, 0 };
  /// Particle position. By default equal to the default inner radius.
  G4ThreeVector fPosition = { 0, 800, 0 };
  /// Flag
  G4bool fIfSet = false;
};

#endif /* ScoreMeshEVENTINFORMATION_HH */
