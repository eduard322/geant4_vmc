// $Id: TG4RootGeometryConvertor.h,v 1.3 2003/12/19 14:39:43 brun Exp $
//
// Author: I. Hrivnacova, 8.1.2003
//
// Class TG4RootGeometryConvertor
// ------------------------------
// Class for converting Root geometry in G4 geometry.

#ifndef TG4_ROOT_GEOMETRY_CONVERTOR_H
#define TG4_ROOT_GEOMETRY_CONVERTOR_H

#include <map>

#include <globals.hh>
#include <geomdefs.hh>

class TGeoVolume;
class TGeoPatternFinder;

class G4VPhysicalVolume;
class G4LogicalVolume;

class TG4MaterialConvertor;
class TG4ShapeConvertor;

class TG4RootGeometryConvertor
{
  public:
    typedef std::map<const TGeoVolume*, G4LogicalVolume*,  
                     std::less<const TGeoVolume*> > VolumesMap;
    typedef VolumesMap::const_iterator  VolumesMapIterator;	       

  public:
    TG4RootGeometryConvertor();
    // --> protected
    // TG4RootGeometryConvertor(const TG4RootGeometryConvertor& right);
    virtual ~TG4RootGeometryConvertor();

    // methods
    G4VPhysicalVolume* Convert(const TGeoVolume* world);
    
    // set methods
    void SetSeparator(char separator);

    // get methods
    char GetSeparator() const;
    
  protected:  
    TG4RootGeometryConvertor(const TG4RootGeometryConvertor& right);

  private:  
    // methods
    void Append(G4String& string, G4int number) const;
    void SetUniqueName(G4LogicalVolume* lv, G4int number) const;
    EAxis GetAxis(const TGeoPatternFinder* finder) const;
    Int_t GetIAxis(EAxis axis) const;
    Bool_t IsDivided(const TGeoVolume* volume);
    G4LogicalVolume* CreateLV(const TGeoVolume* volume);
    void CreatePlacements(const TGeoVolume* volume, G4LogicalVolume* mLV);
    void CreateDivision(const TGeoVolume* volume, G4LogicalVolume* mLV);
    
    void SetUniqueNames();
    void ProcessDaughters(G4LogicalVolume* motherLV, const TGeoVolume* mother);
    void ProcessPositions();
    
    // static data members
    static const char  fgDefaultSeparator;

    // data members
    TG4MaterialConvertor* fMaterialConvertor;
    TG4ShapeConvertor*    fShapeConvertor;
    VolumesMap            fVolumesMap;
    char                  fSeparator;
};

// inline functions

inline void TG4RootGeometryConvertor::SetSeparator(char separator)
{ fSeparator = separator; }

inline char TG4RootGeometryConvertor::GetSeparator() const
{ return fSeparator; }

#endif //TG4_ROOT_GEOMETRY_CONVERTOR_H

