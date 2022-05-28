// Class definition for EventAction().
// Created by Codey Olson on May 7, 2021.

/// \file EventAction.hh
/// \brief Definition of EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

// Event Action:
// Define actions during Geant4 events:

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* );
    virtual void EndOfEventAction(const G4Event* );

};

#endif