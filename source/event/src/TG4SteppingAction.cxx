// $Id: TG4SteppingAction.cxx,v 1.5 2004/11/10 11:39:27 brun Exp $
// Category: event
//
// Class TG4SteppingAction
// -----------------------
// See the class description in the header file.
//
// Author: I.Hrivnacova

#include "TG4SteppingAction.h"
#include "TG4SensitiveDetector.h"
#include "TG4SDServices.h"
#include "TG4Globals.h"

#include <G4Track.hh>
#include <G4SteppingManager.hh>

#include <TVirtualMCApplication.h>

// static data members
TG4SteppingAction* TG4SteppingAction::fgInstance = 0;

//_____________________________________________________________________________
TG4SteppingAction::TG4SteppingAction() 
  : fMessenger(this),
    fMaxNofSteps(kMaxNofSteps),
    fStandardVerboseLevel(-1),
    fLoopVerboseLevel(1),
    fLoopStepCounter(0)
 {
//
  if (fgInstance) { 
    TG4Globals::Exception("TG4SteppingAction constructed twice."); 
  }

  fgInstance = this;
}

//_____________________________________________________________________________
TG4SteppingAction::TG4SteppingAction(const TG4SteppingAction& right)  
  : fMessenger(this)
 {
//
  TG4Globals::Exception("TG4SteppingAction is protected from copying.");
}

//_____________________________________________________________________________
TG4SteppingAction::~TG4SteppingAction() {
//
}

//
// operators
//

//_____________________________________________________________________________
TG4SteppingAction& 
TG4SteppingAction::operator=(const TG4SteppingAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  TG4Globals::Exception("TG4SteppingAction is protected from assigning.");

  return *this;
}

//
// protected methods
//

//_____________________________________________________________________________
void TG4SteppingAction::PrintTrackInfo(const G4Track* track) const
{
/// Print the track info, 
/// taken from private G4TrackingManager::Verbose(), 
/// and the standard header for verbose tracking, 
/// taken from G4SteppingVerbose::TrackingStarted().

  // print track info
  G4cout << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  G4cout << "* G4Track Information: " 
         << "  Particle = " << track->GetDefinition()->GetParticleName() 
         << "," 
	 << "   Track ID = " << track->GetTrackID() 
         << "," 
	 << "   Parent ID = " << track->GetParentID() 
         << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  G4cout << G4endl;
      
  // print header    
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
    G4cout << std::setw( 5) << "Step#"  << " "
           << std::setw( 8) << "X"      << "     "
	   << std::setw( 8) << "Y"      << "     "  
	   << std::setw( 8) << "Z"      << "     "
	   << std::setw( 9) << "KineE"  << "     "
	   << std::setw( 8) << "dE"     << "     "  
	   << std::setw(12) << "StepLeng"   << " "  
	   << std::setw(12) << "TrackLeng"  << " "
	   << std::setw(12) << "NextVolume" << " "
	   << std::setw( 8) << "ProcName"   << G4endl;	     
#else
    G4cout << std::setw( 5) << "Step#"      << " "
	   << std::setw( 8) << "X(mm)"      << " "
	   << std::setw( 8) << "Y(mm)"      << " "  
	   << std::setw( 8) << "Z(mm)"      << " "
	   << std::setw( 9) << "KinE(MeV)"  << " "
	   << std::setw( 8) << "dE(MeV)"    << " "  
	   << std::setw( 8) << "StepLeng"   << " "  
	   << std::setw( 9) << "TrackLeng"  << " "
	   << std::setw(11) << "NextVolume" << " "
	   << std::setw( 8) << "ProcName"   << G4endl;	     
#endif
}

//
// public methods
//

//_____________________________________________________________________________
void TG4SteppingAction::UserSteppingAction(const G4Step* step)
{
/// Called by G4 kernel at the end of each step.
 
  G4Track* track = step->GetTrack();  
  G4int stepNumber = track->GetCurrentStepNumber();

  // stop track if maximum number of steps has been reached
  // 
  if (fLoopStepCounter) {    
     if (stepNumber == 1 ) { 
        // reset parameters at beginning of tracking
        fpSteppingManager->SetVerboseLevel(fStandardVerboseLevel);
        fLoopStepCounter = 0;
  	        // necessary in case particle has reached fMaxNofSteps
	        // but has stopped before processing kMaxNofLoopSteps
     }	
     else {
      // count steps after detecting looping track
      fLoopStepCounter++;
      if (fLoopStepCounter == kMaxNofLoopSteps) {

         // stop the looping track
         track->SetTrackStatus(fStopAndKill); 

         // reset parameters back
         fpSteppingManager->SetVerboseLevel(fStandardVerboseLevel);
         fLoopStepCounter = 0;
      }          
    }
  }    
  else if (stepNumber> fMaxNofSteps) { 
      
    // print looping info
    if (fLoopVerboseLevel > 0) {
       G4cout << "*** Particle reached max step number ("
	      << fMaxNofSteps << "). ***" << G4endl;
	  if (fStandardVerboseLevel == 0) PrintTrackInfo(track);
    }	
   
    // keep standard verbose level
    if (fStandardVerboseLevel<0) 
      fStandardVerboseLevel = fpSteppingManager->GetverboseLevel();

    // set loop verbose level 
    fpSteppingManager->SetVerboseLevel(fLoopVerboseLevel);
      
    // start looping counter
    fLoopStepCounter++;
  }

  // stop track if a user defined tracking region
  // has been reached
  G4ThreeVector position 
    = step->GetPostStepPoint()->GetPosition();

  if (position.mag()    > TVirtualMCApplication::Instance()->TrackingRmax() ||
      std::abs(position.z()) > TVirtualMCApplication::Instance()->TrackingZmax()) {
 
    // print looping info
    if (fLoopVerboseLevel > 0) {
      G4cout << "*** Particle has reached user defined tracking region. ***" 
             << G4endl;
      if (fStandardVerboseLevel == 0) PrintTrackInfo(step->GetTrack());
    }  

    // stop the track
    step->GetTrack()->SetTrackStatus(fStopAndKill);  
  }  
        
  // call stepping action of derived class
  SteppingAction(step);

  // let sensitive detector process boundary step
  // if crossing geometry border
  // (this ensures compatibility with G3 that
  // makes boundary step of zero length)

  if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary &&
      step->GetTrack()->GetTrackStatus() == fAlive &&
      step->GetTrack()->GetNextVolume() != 0) {

#ifdef MCDEBUG
    TG4SensitiveDetector* tsd
      = TG4SDServices::Instance()
           ->GetSensitiveDetector(
	        step->GetPostStepPoint()->GetPhysicalVolume()
		    ->GetLogicalVolume()->GetSensitiveDetector());

    if (tsd) tsd->ProcessHitsOnBoundary((G4Step*)step);
#else
    TG4SensitiveDetector* tsd
      = (TG4SensitiveDetector*) step->GetPostStepPoint()->GetPhysicalVolume()
          ->GetLogicalVolume()->GetSensitiveDetector();

    if (tsd) tsd->ProcessHitsOnBoundary((G4Step*)step);
#endif     
  }  
}

