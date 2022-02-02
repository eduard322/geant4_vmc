//------------------------------------------------
// The Geant4 Virtual Monte Carlo package
// Copyright (C) 2007 - 2014 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file TG4EventAction.cxx
/// \brief Implementation of the TG4EventAction class
///
/// \author I. Hrivnacova; IPN, Orsay

#include "TG4EventAction.h"
#include "TG4Globals.h"
#include "TG4ParticlesManager.h"
#include "TG4SDServices.h"
#include "TG4StateManager.h"
#include "TG4TrackManager.h"
#include "TG4TrackingAction.h"

#include <G4Event.hh>
#include <G4Trajectory.hh>
#include <G4TrajectoryContainer.hh>
#include <G4UImanager.hh>
#include <G4VVisManager.hh>
#include <Randomize.hh>

#include <RVersion.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TVirtualMCApplication.h>
#include <TVirtualMCSensitiveDetector.h>
#include <TVirtualMCStack.h>

#include <math.h>
#include<TG4SteppingAction.h>

G4int Nev;
G4int wiriteFlag;
std::vector<eventInfo> eInf;
//_____________________________________________________________________________
TG4EventAction::TG4EventAction()
  : TG4Verbose("eventAction"),
    fMessenger(this),
    fTimer(),
    fMCApplication(0),
    fMCStack(0),
    fTrackingAction(0),
    fTrackManager(0),
    fStateManager(0),
    fPrintMemory(false),
    fSaveRandomStatus(false),
    fIsInterruptibleEvent(false)
{

  /// Default constructor
}

//_____________________________________________________________________________
TG4EventAction::~TG4EventAction()
{
  fclose(fp);
  /// Destructor
}

//
// public methods
//

//_____________________________________________________________________________
void TG4EventAction::LateInitialize()
{
  /// Cache thread-local pointers
  Nev = 0;
  fp=fopen( "muData.csv","a");
  fprintf(fp,"eventID  trackID  pid  cScat  muBrems  pre_E  pre_px  pre_py  pre_pz  post_E  post_px  post_py  post_pz  x  y  z\n");
  fMCApplication = TVirtualMCApplication::Instance();
  fTrackingAction = TG4TrackingAction::Instance();
  fTrackManager = TG4TrackManager::Instance();
  fStateManager = TG4StateManager::Instance();
}

//_____________________________________________________________________________
void TG4EventAction::BeginOfEventAction(const G4Event* event)
{
  /// Called by G4 kernel at the beginning of event.

  // reset the tracks counters
  fTrackingAction->PrepareNewEvent();

  // fill primary particles in VMC stack if stack is empty
  if (fMCStack->GetNtrack() == 0) {
    if (VerboseLevel() > 0)
      G4cout << "Filling VMC stack with primaries" << G4endl;

    for (G4int iv = 0; iv < event->GetNumberOfPrimaryVertex(); iv++) {
      G4PrimaryVertex* vertex = event->GetPrimaryVertex(iv);

      for (G4int ip = 0; ip < vertex->GetNumberOfParticle(); ip++) {
        G4PrimaryParticle* particle = vertex->GetPrimary(ip);
        fTrackManager->PrimaryToStack(vertex, particle);
      }
    }
  }

  // save the event random number status per event
  if (fSaveRandomStatus) {
    G4UImanager::GetUIpointer()->ApplyCommand("/random/saveThisEvent");
    if (VerboseLevel() > 0) G4cout << "Saving random status: " << G4endl;
    CLHEP::HepRandom::showEngineStatus();
    G4cout << G4endl;
  }

  if (VerboseLevel() > 0) {
    G4cout << ">>> Event " << event->GetEventID() << G4endl;
    fTimer.Start();
  }
  Nev = event->GetEventID();
}

//_____________________________________________________________________________
void TG4EventAction::EndOfEventAction(const G4Event* event)
{
  /// Called by G4 kernel at the end of event.
  if (wiriteFlag){
    for (auto info:eInf)
    {
      fprintf(fp,"%d  %d  %d  %d  %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e\n",
               info.eventID, info.trackID, info.pid, info.cScat, info.muBrems, info.pre_E, info.pre_px, info.pre_py, 
               info.pre_pz, info.post_E, info.post_px, info.post_py, info.post_pz, info.x, info.y, info.z);
    }
  }
  wiriteFlag = 0;
  eInf.clear();
  // finish the last primary track of the current event
  // G4cout << "Finish primary from event action" << G4endl;
  fTrackingAction->FinishPrimaryTrack();

  if (VerboseLevel() > 1) {
    G4cout << G4endl;
    G4cout << ">>> End of Event " << event->GetEventID() << G4endl;
  }

  if (VerboseLevel() > 2) {
    G4int nofPrimaryTracks = fMCStack->GetNprimary();
    G4int nofSavedTracks = fMCStack->GetNtrack();

    G4cout << "    " << nofPrimaryTracks << " primary tracks processed."
           << G4endl;
    G4cout << "    " << nofSavedTracks << " tracks saved." << G4endl;

    G4int nofAllTracks = fTrackManager->GetNofTracks();
    G4cout << "    " << nofAllTracks << " all tracks processed." << G4endl;
  }

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 18, 0)
  // VMC application end of event
  fMCApplication->EndOfEvent();
#endif

  // User SDs finish event
  if (TG4SDServices::Instance()->GetUserSDs()) {
    for (auto& userSD : (*TG4SDServices::Instance()->GetUserSDs())) {
      userSD->EndOfEvent();
    }
  }

  // VMC application finish event
  if (!fIsInterruptibleEvent) {
    fMCApplication->FinishEvent();
  }
  fStateManager->SetNewState(kNotInApplication);

  if (VerboseLevel() > 1) {
    // print time
    fTimer.Stop();
    fTimer.Print();
  }

  if (fPrintMemory) {
    ProcInfo_t procInfo;
    gSystem->GetProcInfo(&procInfo);
    G4cout << "Current memory usage: resident " << procInfo.fMemResident
           << ", virtual " << procInfo.fMemVirtual << G4endl;
  }
}
