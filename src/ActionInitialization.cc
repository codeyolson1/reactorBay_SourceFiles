// Source file for ActionInitialization().
// Created by Codey Olson on May 8, 2021.

/// \file ActionInitialization.cc
/// \brief Source code for ActionInitialization class.

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

ActionInitialization::ActionInitialization()
: G4VUserActionInitialization()
{}

//
//

ActionInitialization::~ActionInitialization()
{}

//
//

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
}

//
//

void ActionInitialization::Build() const
{
  SetUserAction(new RunAction);
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(new EventAction);
}