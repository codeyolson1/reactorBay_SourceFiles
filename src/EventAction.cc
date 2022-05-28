// Source code for EventAction().
// Created by Codey Olson on May 8, 2021.

/// \file EventAction.cc
/// \brief Source code for EventAction class.

#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4PrimaryVertex.hh"

#include <fstream>
#include <iostream>

EventAction::EventAction() : G4UserEventAction()
{}

//
//

EventAction::~EventAction()
{}

//
//

void EventAction::BeginOfEventAction(const G4Event* )
{

}

//
//

void EventAction::EndOfEventAction(const G4Event* )
{
  
}