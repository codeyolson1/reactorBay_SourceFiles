//
// Codey Olson (8/8/2020)

// G4 Primitive Scorer:
// Incident Angle

#include "G4PSIncidentEnergy.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"
#include <math.h>

G4PSIncidentEnergy::G4PSIncidentEnergy(G4String name, G4int direction, G4int depth) 
  : G4VPrimitiveScorer(name, depth), HCID(-1), EvtMap(0), fDirection(direction)
{
  //DefineUnitAndCategory();
  SetUnit("MeV");
}

G4PSIncidentEnergy::G4PSIncidentEnergy(G4String name, G4int direction, const G4String& unit, G4int depth)
  : G4VPrimitiveScorer(name, depth), HCID(-1), EvtMap(0), fDirection(direction)
  {
    //DefineUnitAndCategory();
    SetUnit(unit);
  }

G4PSIncidentEnergy::~G4PSIncidentEnergy()
{ ; }

G4bool G4PSIncidentEnergy::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4VPhysicalVolume* physVol = preStep->GetPhysicalVolume();
  G4VSolid* solid = 0;

  // Ordinary volume:
  solid = physVol->GetLogicalVolume()->GetSolid();

  G4Box* boxSolid = (G4Box*)(solid);
  G4int dirFlag = IsSelectedSurface(aStep, boxSolid);
  if (dirFlag > 0) {
    if (fDirection == fFlux_InOut || fDirection == dirFlag) {
      G4StepPoint* thisStep = 0;
      if (dirFlag == fFlux_In) {
        thisStep = preStep;
      } else if (dirFlag == fFlux_Out) {
        thisStep = aStep->GetPostStepPoint();
      } else {
        return false;
      }
      G4double KinE = thisStep->GetKineticEnergy();
      KinE *= thisStep->GetWeight();
      G4int index = GetIndex(aStep);
      EvtMap->add(index, KinE);
    }
  }

  return true;
}


G4int G4PSIncidentEnergy::IsSelectedSurface(G4Step* aStep, G4Box* boxSolid)
{
  G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
 
 if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {
   // Entering Geometry:
   G4ThreeVector stppos1 = aStep->GetPreStepPoint()->GetPosition();
   G4ThreeVector localpos1 = theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
   if (std::fabs(localpos1.z() + boxSolid->GetZHalfLength()) < kCarTolerance) {
     return fFlux_In;
   }
 }

 if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
   // Exiting Geometry:
   G4ThreeVector stppos2 = aStep->GetPostStepPoint()->GetPosition();
   G4ThreeVector localpos2 = theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
   if (std::fabs(localpos2.z() + boxSolid->GetZHalfLength()) < kCarTolerance) {
     return fFlux_Out;
   }
 }

  return -1;
}

void G4PSIncidentEnergy::Initialize(G4HCofThisEvent* HCE) 
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(), GetName());
  if (HCID < 0) {
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSIncidentEnergy::EndOfEvent(G4HCofThisEvent*)
{ ; }

void G4PSIncidentEnergy::clear() 
{
  EvtMap->clear();
}

void G4PSIncidentEnergy::DrawAll()
{ ; }

void G4PSIncidentEnergy::PrintAll() 
{
  G4cout << " MultiFunctionalDet " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Numer of entries " << EvtMap->entries() << G4endl;
  std::map<G4int, G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for (; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "   copy no.: " << itr->first
    << "   Energy: " << *(itr->second)/GetUnitValue()
    << " [" << GetUnit() << "] " << G4endl;
  }
}

void G4PSIncidentEnergy::SetUnit(const G4String& unit)
{
  CheckAndSetUnit(unit, "Energy");
}