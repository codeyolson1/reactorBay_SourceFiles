// Definition of world geometry and detectors.
// Created by Codey Olson on May 10, 2021.

/// \file DetectorConstruction.cc
/// \brief Definition of world geometry and detectors.

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PSCellFluxToDose.hh"
#include "G4PSCellFluxToEffDose.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSCellEnergy.hh"
#include "G4PSCellFlux.hh"
#include "G4PSIncidentAngle.hh"
#include "G4PSIncidentEnergy.hh"
#include "G4PSIncidentPosition.hh"

#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#define _USE_MATH_DEFINES 
#include <math.h>
#include <iomanip>
#include <iostream>
#include <string>

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(), fdetNames(0.), fnumDets(0)
{
  fdetNames = {"ceiling", "office", "lab", "mesanine", "bayFloor", "surface", "30cm"};
  fnumDets = fdetNames.size(); 
  fmats = {};

  ConstructMaterials();
}


//
//

DetectorConstruction::~DetectorConstruction()
{}

//
//

void DetectorConstruction::ConstructMaterials()
{
  // Get instance of nist material manager:
  G4NistManager* nist = G4NistManager::Instance();

  // Create materials and input into dictionary.
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  fmats["air"] = air;
  G4Material* concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
  fmats["concrete"] = concrete;
  // Plaster data from MCNP Materials Compendium:
  G4Material* plaster = new G4Material("Plaster", 2.32*g/cm3, 4);
  G4Material* hydrogen = nist->FindOrBuildMaterial("G4_H");
  G4Material* oxygen = nist->FindOrBuildMaterial("G4_O");
  G4Material* sulfur = nist->FindOrBuildMaterial("G4_S");
  G4Material* calcium = nist->FindOrBuildMaterial("G4_Ca");
  plaster->AddMaterial(hydrogen, 0.023416*perCent);
  plaster->AddMaterial(oxygen, 0.557572*perCent);
  plaster->AddMaterial(sulfur, 0.186215*perCent);
  plaster->AddMaterial(calcium, 0.232797*perCent);
  fmats["plaster"] = plaster;

  G4Material* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  fmats["steel"] = steel;

  // Material characteristics from ShieldWerx.
  G4Material* poly = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Element* boron = new G4Element("Boron", "B", 2);
  G4Isotope* b10 = new G4Isotope("Boron10", 5, 10, 10.012936862*g/mole);
  G4Isotope* b11 = new G4Isotope("Boron11", 5, 11, 11.0093051662*g/mole);
  boron->AddIsotope(b10, 19.6*perCent);
  boron->AddIsotope(b11, 80.4*perCent);
  G4Material* BPoly5 = new G4Material("5% Borated Polyethylene", 0.96*g/cm3, 2);
  BPoly5->AddElement(boron, 5.0*perCent);
  BPoly5->AddMaterial(poly, 95.0*perCent);
  G4Material* polyPellets = new G4Material("Polyethylene Pellets", 0.564*g/cm3, 1);
  polyPellets->AddMaterial(poly, 100*perCent);
  fmats["poly"] = poly;
  fmats["BPoly5"] = BPoly5;
  fmats["polyPellets"] = polyPellets;

  G4Material* Tantalum = nist->FindOrBuildMaterial("G4_Ta");
  G4Material* Berillyum = nist->FindOrBuildMaterial("G4_Be");
  G4Element* Plutonium = new G4Element("Plutonium", "Pu", 4);
  G4Isotope* Pu239 = new G4Isotope("Plutonium239", 94, 239, 239.0521617*g/mole);
  G4Isotope* Pu240 = new G4Isotope("Plutonium240", 94, 240, 240.0538118*g/mole);
  G4Isotope* Pu241 = new G4Isotope("Plutonium241", 94, 241, 241.0568497*g/mole);
  G4Isotope* Pu242 = new G4Isotope("Plutonium242", 94, 242, 242.0587410*g/mole);
  // Decayed percentages:
  Plutonium->AddIsotope(Pu239, 91.66472859*perCent);
  Plutonium->AddIsotope(Pu240, 8.25386571*perCent);
  Plutonium->AddIsotope(Pu241, 0.0498354615*perCent);
  Plutonium->AddIsotope(Pu242, 0.0315702471*perCent);
  G4Material* PuBe = new G4Material("Pu-Be Mixture", 3.764135004*g/cm3, 2); 
  PuBe->AddElement(Plutonium, 66.39667705*perCent);
  PuBe->AddMaterial(Berillyum, 33.60332295*perCent);
  fmats["PuBe"] = PuBe;
  fmats["Tantalum"] = Tantalum;
  
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;

  //
  // World:
  // Params:
  G4double worldX = 594.36*cm; 
  G4double worldY = 546.1*cm; 
  G4double worldZ = 642.62*cm;
  // Construction:
  G4Box* solidWorld = new G4Box("World", 0.5*worldX, 0.5*worldY,0.5*worldZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, fmats["air"], "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

  //
  // Floor
  // Params:
  G4double floorX = 558.8*cm;
  G4double floorY = 520.7*cm;
  G4double floorZ = 30.48*cm;
  G4ThreeVector floorCenter = G4ThreeVector(0, 0, -290.83*cm);
  // Construction
  G4Box* floorSolid = new G4Box("Floor", 0.5*floorX, 0.5*floorY, 0.5*floorZ);
  G4LogicalVolume* floorLogic = new G4LogicalVolume(floorSolid, fmats["concrete"], "Floor");
  new G4PVPlacement(0, floorCenter, floorLogic, "Floor", logicWorld, false, 0, checkOverlaps);

  //
  // South Wall
  // Params:
  G4double sWallX = 558.8*cm;
  G4double sWallY = 167.64*cm;
  G4double sWallZ = 198.12*cm;
  G4ThreeVector sWallCenter = G4ThreeVector(0, -173.99*cm, -176.53*cm);
  // Construction:
  G4Box* sWallSolid = new G4Box("SouthWall", 0.5*sWallX, 0.5*sWallY, 0.5*sWallZ);
  G4LogicalVolume* sWallLogic = new G4LogicalVolume(sWallSolid, fmats["concrete"], "SouthWall");
  new G4PVPlacement(0, sWallCenter, sWallLogic, "SouthWall", logicWorld, false, 0, checkOverlaps);

  //
  // East Wall
  // Params
  G4double eWallX = 172.72*cm;
  G4double eWallY = 347.98*cm;
  G4double eWallZ = 198.12*cm;
  G4ThreeVector eWallCenter = G4ThreeVector(195.58*cm, 83.82*cm, -176.53*cm);
  // Construction
  G4Box* eWallSolid = new G4Box("EastWall", 0.5*eWallX, 0.5*eWallY, 0.5*eWallZ);
  G4LogicalVolume* eWallLogic = new G4LogicalVolume(eWallSolid, fmats["concrete"], "EastWall");
  new G4PVPlacement(0, eWallCenter, eWallLogic, "EastWall", logicWorld, false, 0, checkOverlaps);

  //
  // Concrete Platform
  // Params:
  G4double platX = 284.48*cm;
  G4double platY = 134.62*cm;
  G4double platZ = 76.2*cm;
  G4ThreeVector platCenter = G4ThreeVector(-139.7*cm, -22.86*cm, -237.49*cm);
  // Construction:
  G4Box* platformSolid = new G4Box("Platform", 0.5*platX, 0.5*platY, 0.5*platZ);
  G4LogicalVolume* platformLogic = new G4LogicalVolume(platformSolid, fmats["concrete"], "Platform");
  new G4PVPlacement(0, platCenter, platformLogic, "Platform", logicWorld, false, 0, checkOverlaps);

  //
  // Pillar
  // Params:
  G4double pillarX = 106.68*cm;
  G4double pillarY = 134.62*cm;
  G4double pillarZ = 246.38*cm;
  G4ThreeVector pillarCenter = G4ThreeVector(55.88*cm, -22.86*cm, -152.4*cm);
  // Construction:
  G4Box* pillarSolid = new G4Box("Pillar", 0.5*pillarX, 0.5*pillarY, 0.5*pillarZ);
  G4LogicalVolume* pillarLogic = new G4LogicalVolume(pillarSolid, fmats["concrete"], "Pillar");
  new G4PVPlacement(0, pillarCenter, pillarLogic, "Pillar", logicWorld, false, 0, checkOverlaps);

  //
  // Mesanine
  // Params:
  G4double mesX = 162.56*cm;
  G4double mesY = 213.36*cm;
  G4double mesZ = 198.12*cm;
  G4ThreeVector mesCenter = G4ThreeVector(27.94*cm, 151.13*cm, -176.53*cm);
  // Construction:
  G4Box* mesanineSolid = new G4Box("Mesanine", 0.5*mesX, 0.5*mesY, 0.5*mesZ);
  G4LogicalVolume* mesanineLogic = new G4LogicalVolume(mesanineSolid, fmats["concrete"], "Mesanine");
  new G4PVPlacement(0, mesCenter, mesanineLogic, "Mesanine", logicWorld, false, 0, checkOverlaps);

  //
  // South Plaster Wall
  // Params:
  G4double sWallPX = 558.8*cm;
  G4double sWallPY = 15.24*cm;
  G4double sWallPZ = 373.38*cm;
  G4ThreeVector sWallPCenter = G4ThreeVector(0, -97.79*cm, 109.22*cm);
  // Construction
  G4Box* sWallPSolid = new G4Box("SouthWallPlaster", 0.5*sWallPX, 0.5*sWallPY, 0.5*sWallPZ);
  G4LogicalVolume* sWallPLogic = new G4LogicalVolume(sWallPSolid, fmats["concrete"], "SouthWallPlaster");
  new G4PVPlacement(0, sWallPCenter, sWallPLogic, "SouthWallPlaster", logicWorld, false, 0, checkOverlaps);

  //
  // East Plaster Wall
  // Params:
  G4double eWallPX = 20.32*cm;
  G4double eWallPY = 347.98*cm;
  G4double eWallPZ = 373.38*cm;
  G4ThreeVector eWallPCenter = G4ThreeVector(119.38*cm, 83.82*cm, 109.22*cm);
  // Construction:
  G4Box* eWallPSolid = new G4Box("EastWallPlaster", 0.5*eWallPX, 0.5*eWallPY, 0.5*eWallPZ);
  G4LogicalVolume* eWallPLogic = new G4LogicalVolume(eWallPSolid, fmats["concrete"], "EastWallPlaster");
  new G4PVPlacement(0, eWallPCenter, eWallPLogic, "EastWallPlaster", logicWorld, false, 0, checkOverlaps);

  //
  // Ceiling Steel
  // Params:
  G4double steelCeilX = 558.8*cm;
  G4double steelCeilY = 520.7*cm;
  G4double steelCeilZ = 0.47625*cm;
  G4ThreeVector steelCeilCenter = G4ThreeVector(0, 0, 296.148125*cm);
  // Construction:
  G4Box* steelCeilSolid = new G4Box("SteelCeiling", 0.5*steelCeilX, 0.5*steelCeilY, 0.5*steelCeilZ);
  G4LogicalVolume* steelCeilLogic = new G4LogicalVolume(steelCeilSolid, fmats["steel"], "SteelCeiling");
  new G4PVPlacement(0, steelCeilCenter, steelCeilLogic, "SteelCeiling", logicWorld, false, 0, checkOverlaps);

  //
  // Ceiling Concrete
  // Params:
  G4double concCeilX = 558.58*cm;
  G4double concCeilY = 520.7*cm;
  G4double concCeilZ = 9.68375*cm;
  G4ThreeVector concCeilCenter = G4ThreeVector(0, 0, 301.228125*cm);
  // Construction:
  G4Box* concCeilSolid = new G4Box("ConcreteCeiling", 0.5*concCeilX, 0.5*concCeilY, 0.5*concCeilZ);
  G4LogicalVolume* concCeilLogic = new G4LogicalVolume(concCeilSolid, fmats["concrete"], "ConcreteCeiling");
  new G4PVPlacement(0, concCeilCenter, concCeilLogic, "ConcreteCeiling", logicWorld, false, 0, checkOverlaps);

  //
  // Storage Drum
  // Params
  G4double drumD = 59.69*cm;
  G4double drumH = 88.265*cm;
  G4ThreeVector drumCenter = G4ThreeVector(-48.26*cm, -22.86*cm, -155.2575*cm);
  // Construction
  G4Tubs* storageDrumOuter = new G4Tubs("StorageDrumOuter", 0, 0.5*drumD, 0.5*drumH, 0, 360.*deg); // Outer Shell
  G4Tubs* storageDrumInner = new G4Tubs("StorageDrumInner", 0, 0.5*drumD - 0.1*cm, 0.5*drumH - 0.1*cm, 0, 360.*deg);
  G4VSolid* storageDrumSolid = new G4SubtractionSolid("StorageDrum", storageDrumOuter, storageDrumInner, 0, G4ThreeVector());
  G4LogicalVolume* storageDrumLogic = new G4LogicalVolume(storageDrumSolid, fmats["steel"], "StorageDrum");
  new G4PVPlacement(0, drumCenter, storageDrumLogic, "StorageDrum", logicWorld, false, 0, checkOverlaps);

  //
  // Main Shield
  // Params:
  G4double shieldTopD = 38.1*cm;
  G4double shieldBottomD = 33.02*cm;
  G4double shieldH = 45.72*cm;
  G4ThreeVector shieldCenter = G4ThreeVector(-48.26*cm, -22.86*cm, -148.49*cm);
  // Construction:
  G4Cons* shieldSolid = new G4Cons("Shield", 0, shieldBottomD*0.5, 0, shieldTopD*0.5, shieldH*0.5, 0, 360.*deg);
  G4LogicalVolume* shieldLogic = new G4LogicalVolume(shieldSolid, fmats["BPoly5"], "Shield");
  new G4PVPlacement(0, shieldCenter, shieldLogic, "Shield", logicWorld, false, 0, checkOverlaps);

  // 
  // Poly Pellets in Storage Drum
  // Params:
  G4double pelletsD = 59.49*cm;
  G4double pelletsH = 63.5*cm;
  G4ThreeVector pelletsCenter = G4ThreeVector(-48.26*cm, -22.86*cm, -167.54*cm);
  // Construction:
  G4Tubs* pelletSolidOuter = new G4Tubs("PelletsOuter", 0, 0.5*pelletsD, 0.5*pelletsH, 0, 360.*deg);
  G4VSolid* pelletSolid = new G4SubtractionSolid("Pellets", pelletSolidOuter, shieldSolid, 0, G4ThreeVector(0, 0, 19.05*cm));
  G4LogicalVolume* pelletLogic = new G4LogicalVolume(pelletSolid, fmats["polyPellets"], "Pellets");
  new G4PVPlacement(0, pelletsCenter, pelletLogic, "Pellets", logicWorld, false, 0, checkOverlaps);

  //
  // Source:
  // Params:
  G4double stainlessD = 2.53492*cm;
  G4double stainlessH = 5.461*cm;
  G4double tantalumD = 2.37236*cm;
  G4double tantalumH = 4.572*cm;
  G4double PuBeD = 2.06756*cm;
  G4double PuBeH = 3.814*cm;
  G4ThreeVector PuBeCenter = G4ThreeVector(0, 0, -6.7675*cm); // Local coords for shield:
  // Construction:
  G4Tubs* stainlessSolid = new G4Tubs("StainlessShell", 0, 0.5*stainlessD, 0.5*stainlessH, 0, 360.*deg);
  G4LogicalVolume* stainlessLogic = new G4LogicalVolume(stainlessSolid, fmats["steel"], "StainlessShell");
  G4VisAttributes* sourceAttr =  new G4VisAttributes(G4Colour(1.,0.,0.));
  sourceAttr->SetForceSolid(true);
  stainlessLogic->SetVisAttributes(sourceAttr);
  new G4PVPlacement(0, PuBeCenter, stainlessLogic, "StainlessShell", shieldLogic, false, 0, checkOverlaps);
  G4Tubs* tantalumSolid = new G4Tubs("TantalumShell", 0, 0.5*tantalumD, 0.5*tantalumH, 0, 360.*deg);
  G4LogicalVolume* tantalumLogic = new G4LogicalVolume(tantalumSolid, fmats["Tantalum"], "TantalumShell");
  new G4PVPlacement(0, G4ThreeVector(), tantalumLogic, "TantalumShell", stainlessLogic, false, 0, checkOverlaps);
  G4Tubs* PuBeSolid = new G4Tubs("PuBeSource", 0, 0.5*PuBeD, 0.5*PuBeH, 0, 360.*deg);
  G4LogicalVolume* PuBeLogic = new G4LogicalVolume(PuBeSolid, fmats["PuBe"], "PuBeSource");
  new G4PVPlacement(0, G4ThreeVector(), PuBeLogic, "PuBe Source", tantalumLogic, false, 0, checkOverlaps);

  //
  // Detectors:
  // Params:
  G4double detectorD = 15.*cm;
  std::map<std::string, G4ThreeVector> storageDetectors;
  std::map<std::string, G4ThreeVector> operationDetectors;
  storageDetectors[fdetNames[0]] = G4ThreeVector(-48.26*cm, -22.86*cm, 288.41*cm); // Ceiling
  storageDetectors[fdetNames[1]] = G4ThreeVector(-48.26*cm, -166.37*cm, 44.45*cm); // Office
  storageDetectors[fdetNames[2]] = G4ThreeVector(190.5*cm, -22.86*cm, 44.45*cm); // Lab
  storageDetectors[fdetNames[3]] = G4ThreeVector(27.94*cm, 105.41*cm, 44.45*cm); // Mesanine
  storageDetectors[fdetNames[4]] = G4ThreeVector(-167.64*cm, 151.13*cm, -153.67*cm);
  storageDetectors[fdetNames[5]] = G4ThreeVector(-85.605*cm, -22.86*cm, -155.2575*cm); // Surface
  storageDetectors[fdetNames[6]] = G4ThreeVector(-115.605*cm, -22.86*cm, -155.2575*cm); //  30 cm
  // Construction:
   for (G4int i = 0; i < fnumDets; i++) {
    G4Sphere* sphere_Solid = new G4Sphere(fdetNames[i], 0, 0.5*detectorD, 0, 360.*deg, 0, 180.*deg);
    G4LogicalVolume* sphere_Logic = new G4LogicalVolume(sphere_Solid, fmats["air"], fdetNames[i]);
    new G4PVPlacement(0, storageDetectors[fdetNames[i]], sphere_Logic, fdetNames[i], logicWorld, false, 0, checkOverlaps);
   }
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  auto nFilter = new G4SDParticleFilter("n-Filter", "neutron");
  auto gFilter = new G4SDParticleFilter("g-Filter", "gamma");

  for (G4int i = 0; i < fnumDets; i++) {
    G4MultiFunctionalDetector* nDetector = new G4MultiFunctionalDetector(fdetNames[i]+"N");
    G4SDManager::GetSDMpointer()->AddNewDetector(nDetector);
    G4VPrimitiveScorer* nDoseEq = new G4PSCellFluxToDose("nEquivalentDose");
    nDoseEq->SetFilter(nFilter);
    nDetector->RegisterPrimitive(nDoseEq);
    G4VPrimitiveScorer* nDoseEff = new G4PSCellFluxToEffDose("nEffectiveDose");
    nDoseEff->SetFilter(nFilter);
    nDetector->RegisterPrimitive(nDoseEff);
    SetSensitiveDetector(fdetNames[i], nDetector);
    
    G4MultiFunctionalDetector* gDetector = new G4MultiFunctionalDetector(fdetNames[i]+"G");
    G4SDManager::GetSDMpointer()->AddNewDetector(gDetector);
    G4VPrimitiveScorer* gDoseEq = new G4PSCellFluxToDose("gEquivalentDose");
    gDoseEq->SetFilter(gFilter);
    gDetector->RegisterPrimitive(gDoseEq);
    G4VPrimitiveScorer* gDoseEff = new G4PSCellFluxToEffDose("gEffectiveDose");
    gDoseEff->SetFilter(gFilter);
    gDetector->RegisterPrimitive(gDoseEff);
    SetSensitiveDetector(fdetNames[i], gDetector);
  }

}