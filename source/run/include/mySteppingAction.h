#ifndef mySteppingAction_h
#define mySteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "TG4SteppingAction.h"

class DetectorConstruction;
class myEventAction;


class mySteppingAction : public TG4SteppingAction
{
public:
  mySteppingAction(const DetectorConstruction* detectorConstruction,
                    myEventAction* eventAction);
  virtual ~mySteppingAction();

  virtual void UserSteppingAction(const G4Step* step);
    
private:
  const DetectorConstruction* fDetConstruction;
  myEventAction*  fEventAction;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
