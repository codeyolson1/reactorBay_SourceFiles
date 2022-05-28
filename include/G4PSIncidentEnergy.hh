//
// Codey Olson (8/8/2020)

// G4 Primitive Scorer:
// Incident energy
//
// Surface is defined at the -Z surface.
// Direction                  -Z   +Z
//   0  IN || OUT            ->|<-  |        fFlux_InOut
//   1  IN                   ->|    |        fFlux_In
//   2  OUT                    |<-  |        fFlux_Out

#ifndef G4PSIncidentEnergy_h
#define G4PSIncidenEnergy_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"
#include "G4Box.hh"
#include "G4PSDirectionFlag.hh"

class G4PSIncidentEnergy : public G4VPrimitiveScorer
{
    public :
        G4PSIncidentEnergy(G4String name, G4int direction, G4int depth=0);
        G4PSIncidentEnergy(G4String name, G4int direction, const G4String& unit, G4int depth=0);
        virtual ~G4PSIncidentEnergy();

    protected :
        virtual G4bool ProcessHits(G4Step* , G4TouchableHistory* );
        G4int IsSelectedSurface(G4Step* , G4Box* );

    public :
        virtual void Initialize(G4HCofThisEvent* );
        virtual void EndOfEvent(G4HCofThisEvent* );
        virtual void clear();
        virtual void DrawAll();
        virtual void PrintAll();
        virtual void SetUnit(const G4String& unit);
        //virtual void DefineUnitAndCategory();

    private :
        G4int HCID;
        G4THitsMap<G4double>* EvtMap;
        G4int fDirection;


};
#endif
