// $Id: Ex02MCApplication.h,v 1.1 2002/06/20 11:52:35 hristov Exp $
//
// Geant4 ExampleN01 adapted to Virtual Monte Carlo 
//
// Class Ex02MCApplication
// ----------------------- 
// Implementation of the TVirtualMCApplication
//
// by Ivana Hrivnacova, 21.4.2002

#ifndef EX02_MC_APPLICATION_H
#define EX02_MC_APPLICATION_H

#include <TVirtualMCApplication.h>

#include "Ex02DetectorConstruction.h"
#include "Ex02TrackerSD.h"
#include "Ex02RootManager.h"

class Ex02MCStack;

class Ex02MCApplication : public TVirtualMCApplication
{
  public:
    Ex02MCApplication(const char* name,  const char *title, 
                      FileMode fileMode = kWrite);
    Ex02MCApplication();
    virtual ~Ex02MCApplication();
  
    // static access method
    static Ex02MCApplication* Instance(); 

    // methods
    void InitMC(const char *setup);
    void RunMC(Int_t nofEvents);
    void FinishRun();
 
    virtual void ConstructGeometry();
    virtual void InitGeometry();
    virtual void GeneratePrimaries();
    virtual void BeginEvent();
    virtual void BeginPrimary();
    virtual void PreTrack();
    virtual void Stepping();
    virtual void PostTrack();
    virtual void FinishPrimary();
    virtual void FinishEvent();
    
    virtual Double_t TrackingRmax() const;
    virtual Double_t TrackingZmax() const;
    virtual void Field(const Double_t* x, Double_t* b) const;
    
    void   ReadEvent(Int_t i);

  private:
    // methods
    void RegisterStack();
  
    // data members
    Ex02MCStack*             fStack;
    Ex02DetectorConstruction fDetConstruction;
    Ex02TrackerSD            fTrackerSD;
    Double_t*                fFieldB;
    Ex02RootManager          fRootManager;

  ClassDef(Ex02MCApplication,1)  //Interface to MonteCarlo application
};

// inline functions

inline Ex02MCApplication* Ex02MCApplication::Instance()
{ return (Ex02MCApplication*)(TVirtualMCApplication::Instance()); }

#endif //EX02_MC_APPLICATION_H
