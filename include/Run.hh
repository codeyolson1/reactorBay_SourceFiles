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
    Run();
    virtual ~Run();

    virtual void RecordEvent(const G4Event* );
    virtual void Merge(const G4Run* );
    void ConstructDetector(const G4String& );

    G4THitsMap<G4double>* GetHitsMap(const G4String&) const;
    G4THitsMap<G4ThreeVector>* GetVectorMaps(const G4String&) const;
    G4StatContainer<G4StatAnalysis>* GetStatMap(const G4String&) const;
    std::vector<std::vector<G4double>> GetUnitVals() const {return funitVals;}
    std::vector<std::vector<G4String>> GetUnitStrings() const {return funitStrings;}
    std::vector<std::vector<G4String>> GetCollNames() const {return fcollNames;}
    std::vector<G4String> GetDetNames() const {return fdetNames;}

  private:
    std::vector<std::vector<G4String>> fcollNames;
    std::vector<std::vector<G4int>> fcollIDs;
    std::vector<std::vector<G4THitsMap<G4double>*>> frunMaps;
    std::vector<std::vector<G4StatContainer<G4StatAnalysis>*>> fstatMaps;
    std::vector<std::vector<G4double>> funitVals;
    std::vector<std::vector<G4String>> funitStrings;
    std::vector<G4String> fdetNames;
    G4AnalysisManager* fanalysisManager;
};

#endif