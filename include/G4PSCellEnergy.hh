//
// Codey Olson (8/3/2020)

// G4 Primitive Scorer: 
// Average cell energy per event.

#ifndef G4PSCellEnergy_h
#define G4PSCellEnergy_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

class G4PSCellEnergy : public G4VPrimitiveScorer
{
    public:
        G4PSCellEnergy(G4String name, G4int depth=0);
        G4PSCellEnergy(G4String name, const G4String& unit, G4int depth=0);
        virtual ~G4PSCellEnergy();

    protected:
        virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    
    public:
        virtual void Initialize(G4HCofThisEvent*);
        virtual void EndOfEvent(G4HCofThisEvent*);
        virtual void clear();
        virtual void DrawAll();
        virtual void PrintAll();
        virtual void SetUnit(const G4String& unit);

    private:
        G4int HCID;
        G4THitsMap<G4double>* EvtMap;

};
#endif