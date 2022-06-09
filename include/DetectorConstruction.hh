// Class Definition of DetectorConstruction().
// Created by Codey Olson on May 8, 2021.

/// \file DetectorConstruction.hh
/// \file Definition of DetectorConstruction class.

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

// Define detector/world geomteries and materials.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction(G4bool);
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

  private:
    std::vector<G4String> fdetNames;
    G4int fnumDets;
    std::map<std::string, G4Material*> fmats;
    G4bool isHe3;

  public:
    inline std::vector<G4String> GetDetNames() const { return fdetNames; }
    inline G4int GetDetNum() const { return fnumDets; }
    void ConstructMaterials();


};

#endif