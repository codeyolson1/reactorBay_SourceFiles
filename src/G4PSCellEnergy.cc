//
// Codey Olson (8/3/2020)

// G4 Primitive Scorer"
// Average cell energy per event

#include "G4PSCellEnergy.hh"
#include "G4UnitsTable.hh"

G4PSCellEnergy::G4PSCellEnergy(G4String name, G4int depth) :
G4VPrimitiveScorer(name, depth), HCID(-1), EvtMap(0)
{
    SetUnit("MeV");
}

G4PSCellEnergy::G4PSCellEnergy(G4String name, const G4String& unit, G4int depth) :
G4VPrimitiveScorer(name, depth), HCID(-1), EvtMap(0)
{
    SetUnit(unit);
}

G4PSCellEnergy::~G4PSCellEnergy()
{;}

G4bool G4PSCellEnergy::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    // Get energy:
    G4double KinE;
    KinE = aStep->GetPreStepPoint()->GetKineticEnergy();
    KinE *= aStep->GetPreStepPoint()->GetWeight();
    // Add to hits map:
    G4int index = GetIndex(aStep);
    EvtMap->add(index, KinE);
    return true;
}

void G4PSCellEnergy::Initialize(G4HCofThisEvent* HCE)
{
    EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(), GetName());
    if (HCID < 0) { 
        HCID = GetCollectionID(0);
    }
    HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSCellEnergy::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSCellEnergy::clear()
{
    EvtMap->clear();
}

void G4PSCellEnergy::DrawAll()
{;}

void G4PSCellEnergy::PrintAll()
{
    G4cout << " MultiFuntionalDet  " << detector->GetName() << G4endl;
    G4cout << " PrimitiveScorer " << GetName() << G4endl;
    G4cout << " Number of entries " << EvtMap->entries() << G4endl;
    std::map<G4int, G4double*>::iterator itr = EvtMap->GetMap()->begin();
    for (; itr != EvtMap->GetMap()->end(); itr++) {
        G4cout << "  copy no,: " << itr->first
         << " energy: " << *(itr->second)/GetUnitValue() << " [" << GetUnit() <<"] " << G4endl;
    }
}

void G4PSCellEnergy::SetUnit(const G4String& unit)
{
    CheckAndSetUnit(unit, "Energy");
}