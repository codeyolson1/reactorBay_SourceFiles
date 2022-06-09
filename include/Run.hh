// Class definition of Run().
// Created by Codey Olson on May 7, 2021;

/// \file Run.hh
/// \brief Class definition of Run.

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4StatAnalysis.hh"
#include "G4THitsMap.hh"
#include "G4THitsVector.hh"
#include "Analysis.hh"
#include "G4ThreeVector.hh"

class G4Event;

// User defined run class for run statistics.
// Use for thread local run data.

template <typename _Tp> using G4StatContainer = G4THitsDeque<_Tp>;

class Run : public G4Run
{
  public:
    Run(G4bool);
    virtual ~Run();

    virtual void RecordEvent(const G4Event* );
    virtual void Merge(const G4Run* );

  private:
    G4bool isHe3;
};

#endif