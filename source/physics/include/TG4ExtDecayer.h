// $Id: TG4ExtDecayer.h,v 1.1 2002/06/20 11:57:45 hristov Exp $
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4ExtDecayer
// -------------------
// TG4ExtDecayer class implements the G4VExtDecayer abstract class
// with the TVirtualMCDecayer.
// In case a particle has not defined any decay channel
// and has not pre-assigned decay products,
// the external decayer is called.

#ifndef TG4_EXT_DECAYER_H
#define TG4_EXT_DECAYER_H

#include "TG4Verbose.h"

#include "G4VExtDecayer.hh"
#include "globals.hh"

class TVirtualMCDecayer;
class TG4ParticlesManager;

class G4Track;
class G4DecayProducts;

class TClonesArray;

class TG4ExtDecayer : public G4VExtDecayer, 
                      public TG4Verbose
{
  public:
    TG4ExtDecayer(TVirtualMCDecayer* externalDecayer);
    // --> protected
    //TG4ExtDecayer(const TG4ExtDecayer& right);
    virtual ~TG4ExtDecayer();

    virtual G4DecayProducts* ImportDecayProducts(const G4Track& track);
    
  protected:  
    TG4ExtDecayer(const TG4ExtDecayer& right);

    // operators
    TG4ExtDecayer& operator=(const TG4ExtDecayer& right);

  private:
    TG4ParticlesManager* fParticlesManager;  //particles manager 
    TVirtualMCDecayer*   fExternalDecayer;   //the external decayer
    TClonesArray*        fDecayProductsArray;//array of decay products
};

#endif //TG4_EXT_DECAYER_H