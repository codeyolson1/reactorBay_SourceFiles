// Source file for ActionInitialization().
// Created by Codey Olson on May 8, 2021.

/// \file ActionInitialization.cc
/// \brief Source code for ActionInitialization class.

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

ActionInitialization::ActionInitialization(G4bool he3)
: G4VUserActionInitialization()
{
  isHe3 = he3;
}

//
//

ActionInitialization::~ActionInitialization()
{}

//
//

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction(isHe3));
}

//
//

void ActionInitialization::Build() const
{
  SetUserAction(new RunAction(isHe3));
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(new EventAction);
}