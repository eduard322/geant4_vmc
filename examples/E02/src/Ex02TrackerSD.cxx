// $Id: Ex02TrackerSD.cxx,v 1.1.1.1 2002/06/16 15:57:36 hristov Exp $
//
// Geant4 ExampleN02 adapted to Virtual Monte Carlo 
//
/// Id: ExN02TrackerSD.cc,v 1.6 2002/01/09 17:24:10 ranjard Exp 
// GEANT4 tag Name: geant4-04-00-patch-02 
//
// by Ivana Hrivnacova, 21.4.2002

#include <iostream>

#include <TVirtualMC.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include "Ex02TrackerSD.h"
#include "Ex02TrackerHit.h"
#include "Ex02RootManager.h"

ClassImp(Ex02TrackerSD)

using namespace std;

//_____________________________________________________________________________
Ex02TrackerSD::Ex02TrackerSD(const char* name)
  : TNamed(name, ""),
    fTrackerCollection(0),
    fVerboseLevel(1)
{
  fTrackerCollection = new TClonesArray("Ex02TrackerHit");
}

//_____________________________________________________________________________
Ex02TrackerSD::Ex02TrackerSD()
  : TNamed(),
    fTrackerCollection(0),
    fVerboseLevel(1)
{}

//_____________________________________________________________________________
Ex02TrackerSD::~Ex02TrackerSD()
{}

//
// private methods
//

//_____________________________________________________________________________
Ex02TrackerHit* Ex02TrackerSD::AddHit()
{
// Creates a new hit in the TClonesArray.
// ---

  TClonesArray& ref = *fTrackerCollection;
  Int_t size = ref.GetEntriesFast();

  return new(ref[size]) Ex02TrackerHit();
}

//
// public methods
//

//_____________________________________________________________________________
void Ex02TrackerSD::Initialize()
{
// Registers hits collection in Root manager;
// sets sensitive volumes.
// ---
  
  Register();
  
  fSensitiveVolumeID = gMC->VolId("CHMB");
}

//_____________________________________________________________________________
Bool_t Ex02TrackerSD::ProcessHits()
{
// Creates hits (in stepping).
// ---

  Int_t copyNo;
  Int_t id = gMC->CurrentVolID(copyNo);

  if (id != fSensitiveVolumeID) return false;

  Double_t edep = gMC->Edep();

  if (edep==0.) return false;

  Ex02TrackerHit* newHit = AddHit();

  // Track ID
  newHit->SetTrackID  (gMC->GetStack()->CurrentTrack());

  // Chamber no
  newHit->SetChamberNb(copyNo);

  // Energy deposit
  newHit->SetEdep     (edep);

  // Position
  TLorentzVector pos;
  gMC->TrackPosition(pos);
  newHit->SetPos (TVector3(pos.X(), pos.Y(), pos.Z()));
  
  //newHit->Print();
  //newHit->Draw();

  return true;
}

//_____________________________________________________________________________
void Ex02TrackerSD::EndOfEvent()
{
// Prints hits collection (if verbose)
// and deletes hits afterwards.
// ---

  if (fVerboseLevel>0)  Print();
    
  // Reset hits collection
  fTrackerCollection->Delete();  
}

//_____________________________________________________________________________
void Ex02TrackerSD::Register()
{
// Registers the hits collection in Root manager.
// ---
  
  Ex02RootManager::Instance()->Register("hits", &fTrackerCollection);
}

//_____________________________________________________________________________
void Ex02TrackerSD::Print() const
{
// Prints the hits collection.
// ---
  
   Int_t nofHits = fTrackerCollection->GetEntriesFast();
     
   cout << "\n-------->Hits Collection: in this event they are " << nofHits 
        << " hits in the tracker chambers: " << endl;
	    
   for (Int_t i=0; i<nofHits; i++) (*fTrackerCollection)[i]->Print();          
}