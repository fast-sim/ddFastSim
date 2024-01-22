#include "ddfastsim/DefineMeshModel.h"
#include "ddfastsim/ScoreMeshEventInformation.h"
#include "G4EventManager.hh"
#include "G4Event.hh"

dd4hep::sim::DefineMeshModel::DefineMeshModel(dd4hep::sim::Geant4Context* context, const std::string& nam)
      : Geant4FastSimShowerModel(context, nam)
    {}
  bool dd4hep::sim::DefineMeshModel::check_trigger(const G4FastTrack& track) {
   auto  info = dynamic_cast<ScoreMeshEventInformation*>(
    G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation());
    // check if particle direction and position were already set for this event
    if(info != nullptr)
      return !info->GetFlag();
    else
      return true;
  }

     void dd4hep::sim::DefineMeshModel::modelShower(const G4FastTrack& aTrack, G4FastStep& aStep) {
      auto  info = dynamic_cast<ScoreMeshEventInformation*>(
        G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation());
      if(info == nullptr)
      {
        info = new ScoreMeshEventInformation();
        G4EventManager::GetEventManager()->GetNonconstCurrentEvent()->SetUserInformation(info);
      }
      info->SetPosition(aTrack.GetPrimaryTrack()->GetPosition());
      info->SetDirection(aTrack.GetPrimaryTrack()->GetMomentumDirection());
      info->SetFlag(true);
      info->Print();
    }


#include <DDG4/Factories.h>
DECLARE_GEANT4ACTION(DefineMeshModel)
