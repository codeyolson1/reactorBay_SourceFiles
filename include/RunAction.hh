// Definition of class RunAction().
// Created by Codey Olson on May 8, 2021.

/// \file RunAction.hh
/// \brief Definition of the RunAction class.

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Accumulable.hh"
#include "G4StatAnalysis.hh"
#include "Analysis.hh"
#include "G4SDManager.hh"
#include <iostream>

class G4Run;

// RunAction class
// Dictate beginning and end of run actions:

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run* );
    virtual void EndOfRunAction(const G4Run* );

  private:
    G4SDManager* fsdPointer;
    G4AnalysisManager* fanalysisManager;
    std::vector<G4double> nEdges;
    G4String fName;
};

#endif