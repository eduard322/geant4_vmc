// @(#)root/vmc:$Name:  $:$Id$
// Author: Ivana Hrivnacova, IPN Orsay 17/02/2012

/*************************************************************************
 * Copyright (C) 2006, Rene Brun and Fons Rademakers.                    *
 * Copyright (C) 2012, Ivana Hrivnacova.                         *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMCRootManagerImpl
#define ROOT_TMCRootManagerImpl

#include "TVirtualMCRootManager.h"

class TParticle;
class TFile;
class TTree;

//
// Class TMCRootManagerImpl
// ------------------------
// The common implementation of the TVirtualMCRootManager interface
// for the Root IO managers for VMC examples.
// It is used in TMCRootManager (for single-threaded applications)
// and TMCRootManagerMT (for multi-threaded applications)

class TMCRootManagerImpl
{
  public:
    TMCRootManagerImpl(const char* projectName, 
                       TVirtualMCRootManager::FileMode fileMode 
                         = TVirtualMCRootManager::kWrite, 
                       Int_t threadRank = -1);
    virtual ~TMCRootManagerImpl();     
  
    // methods
    void  Register(const char* name, const char* className, void* objAddress);
    void  Fill();
    void  WriteAll();
    void  Close();
    void  ReadEvent(Int_t i);
    
  private:
    // not implemented
    TMCRootManagerImpl(const TMCRootManagerImpl& rhs);
    TMCRootManagerImpl& operator=(const TMCRootManagerImpl& rhs);
    
    // data members
    TFile*  fFile;       // Root output file
    TTree*  fTree;       // Root output tree 
    Bool_t  fIsClosed;   // Info whether its file was closed
};

#endif //ROOT_TMCRootManagerImpl