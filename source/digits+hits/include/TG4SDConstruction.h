// $Id: TG4SDConstruction.h,v 1.1 2002/06/20 11:53:25 hristov Exp $
// Category: digits+hits
//
// Author: I.Hrivnacova
//
// Class TG4SDConstruction
// ------------------------
// Abstract class for construction of user sensitive detectors.
// It has one pure virtual method Construct()
// that has to be be implemented by a derived class.
// Constructed sensitive detectors have to inherit from 
// TG4VSensitiveDetector (see TG4VSensitiveDetector.h description);
// all cloned logical volumes (which a single G3 volume correspond to)
// have to share the same sensitive detector instance.

#ifndef TG4_SD_CONSTRUCTION_H
#define TG4_SD_CONSTRUCTION_H

#include "TG4Verbose.h"

class G4LogicalVolume;

class TG4SDConstruction : public TG4Verbose

{
  public:
    TG4SDConstruction();
    virtual ~TG4SDConstruction();

    // methods
    virtual void Construct();

  private:
    // methods
    void CreateSD(G4LogicalVolume* lv) const;
};

#endif //TG4_SD_CONSTRUCTION_H
