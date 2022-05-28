// Source code for class Run().
// Created by Codey Olson on May 8, 2021.

/// \file Run.cc
/// \brief Source code for Run class.

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Threading.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4THitsMap.hh"
#include "G4HCofThisEvent.hh"
#include "G4PrimaryVertex.hh"
#include "G4PSCellFluxToDose.hh"
#include "G4PSCellEnergy.hh"
#include "G4PSIncidentAngle.hh"
#include "G4PSIncidentEnergy.hh"
#include "G4PSIncidentPosition.hh"
#include "G4PSCellFlux.hh"
#include "G4ThreeVector.hh"
#include "G4WorkerThread.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4StatAnalysis.hh"

Run::Run()
: G4Run()
{
  fdetNames = {"ceilingN", "ceilingG", "officeN", "officeG", "labN", "labG","mesanineN", "mesanineG", "bayFloorN", "bayFloorG", "surfaceN", "surfaceG","30cmN", "30cmG"};

  for (G4int i = 0; i < (G4int)fdetNames.size(); i++) {
    ConstructDetector(fdetNames[i]);
  }
}

//
//

Run::~Run()
{
  for (G4int i = 0; i < (G4int)fdetNames.size(); i++) {
    for (G4int j = 0; j < (G4int)frunMaps[i].size(); j++) {
      delete frunMaps[i].at(j);
    }
  }

}

// 
//

void Run::ConstructDetector(const G4String& detName)
{
  G4cout << "Constructing Detector: " << detName << G4endl;
  // Connect to instance of Sensitive Detector Manager:
  G4SDManager* sdMan = G4SDManager::GetSDMpointer();
  // Find detector from the passed name:
  G4MultiFunctionalDetector* mfd = (G4MultiFunctionalDetector*)(sdMan->FindSensitiveDetector(detName));
  // Make sure the detector was found:
  if (mfd) {
    // Create temporary vectors
    std::vector<G4String> temp_collNames;
    std::vector<G4int> temp_collIDs;
    std::vector<G4THitsMap<G4double>*> temp_runMaps;
    std::vector<G4StatContainer<G4StatAnalysis>*> temp_statMaps;
    std::vector<G4double> temp_unitVals;
    std::vector<G4String> temp_unitStrings;
    // Loop through PS and append to vectors:
    G4int numPrims = mfd->GetNumberOfPrimitives();
    for (G4int i = 0; i < numPrims; i++) {
      // Get ith PS attached to mfd:
      G4VPrimitiveScorer* scorer = mfd->GetPrimitive(i);
      // Collect data from the scorer:
      G4String collName = scorer->GetName();
      G4String fullName = detName +"/"+ collName;
      G4int collID = sdMan->GetCollectionID(fullName);
      // Make sure the collection is found and not empty:
      if (collID >= 0) {
        G4cout << "---> Primitive Scorer: " << fullName << ", ID: " << collID << G4endl;
        temp_collNames.push_back(fullName);
        temp_collIDs.push_back(collID);
        temp_runMaps.push_back(new G4THitsMap<G4double>(detName, collName));
        temp_statMaps.push_back(new G4StatContainer<G4StatAnalysis>(detName, collName));
        temp_unitVals.push_back(scorer->GetUnitValue());
        temp_unitStrings.push_back(scorer->GetUnit());
      }
    }
    fcollNames.push_back(temp_collNames);
    fcollIDs.push_back(temp_collIDs);
    frunMaps.push_back(temp_runMaps);
    fstatMaps.push_back(temp_statMaps);
    funitVals.push_back(temp_unitVals);
    funitStrings.push_back(temp_unitStrings);
  }

}

//
//

void Run::RecordEvent(const G4Event* event)
{
  // Always record event!!
  G4Run::RecordEvent(event);
  // Get the event ID for this event:
  const G4int eID = event->GetEventID();
  // Print event info every 10000 events:
  if (eID % 10000 == 0) {
    std::cout << "---> Event # " << eID << " started on thread: " << G4Threading::G4GetThreadId() << std::endl;
  }
  // Connect to analysis manager:
  fanalysisManager = G4AnalysisManager::Instance();
  // Get Primary Energy information:
  G4PrimaryVertex* pVertex = event->GetPrimaryVertex();
  G4PrimaryParticle* primary = pVertex->GetPrimary();
  G4double energy = primary->GetKineticEnergy();
  fanalysisManager->FillH1(0, energy);

  // Get the hit collections for this event:
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  // Make sure the collection was found:
  if (!hce) return;
  // Loop through all detectors:
  for (G4int i = 0; i< (G4int)fdetNames.size(); i++) {
    G4bool addRow = false;
    for (G4int j = 0; j < (G4int)fcollIDs[i].size(); j++) {
      G4int collID = fcollIDs[i].at(j);
      // Gather hits from the event into the run data:
      G4THitsMap<G4double>* eventMap = 0;
      if (collID >= 0) {
        eventMap = static_cast<G4THitsMap<G4double>*>(hce->GetHC(collID));
      } else {
        G4cout << "** Event Map Not Found **" << G4endl;
      }
      // Make sure there is at least 1 entry in the map:
      if (eventMap && eventMap->entries() >= 1) {
        addRow = true;
        *frunMaps[i][j] += *eventMap;
        *fstatMaps[i][j] += *eventMap;
        // Accumulate all entries in the map:
        G4double val = 0;
        for (auto itr = eventMap->begin(); itr != eventMap->end(); itr++) {
          val += *itr->second;
        }
        fanalysisManager->FillNtupleDColumn(i, j, val/funitVals[i][j]);
      }
    }
    if (addRow) fanalysisManager->AddNtupleRow(i);
  }

}

//
//

void Run::Merge(const G4Run* aRun)
{
  const Run* localRun = static_cast<const Run*>(aRun);
  // Merge all maps for each run (1 run per thread):
  for (G4int i = 0; i < (G4int)fdetNames.size(); i++) {
    for (G4int j = 0; j < (G4int)frunMaps[i].size(); j++) {
      *frunMaps[i][j] += *localRun->frunMaps[i][j];
      *fstatMaps[i][j] += *localRun->fstatMaps[i][j];
    }
  }
  
  G4Run::Merge(aRun);

}

//
//

G4THitsMap<G4double>* Run::GetHitsMap(const G4String& collName) const
{
  // Loop through collections to find the requested collection name:
  for (G4int i = 0; i < (G4int)fcollNames.size(); i++) {
    for (G4int j = 0; j < (G4int)fcollNames[i].size(); j++) {
      if (collName == fcollNames[i][j]) {
        return frunMaps[i][j];
      }
    }
  }

  // None were found -> Raise exception:
  G4Exception("Run", collName.c_str(), JustWarning, "GetHitsMap failed to locate requested HitsMap.");
  return nullptr;

}

//
//

G4StatContainer<G4StatAnalysis>* Run::GetStatMap(const G4String& collName) const
{
  // Loop through collections to find the requested collection name:
  for (G4int i = 0; i< (G4int)fcollNames.size(); i++) {
    for (G4int j = 0; j < (G4int)fcollNames[i].size(); j++) {
      if (collName == fcollNames[i][j]) {
        return fstatMaps[i][j];
      }
    }
  }

  // None were found -> Raise exception:
    G4Exception("Run", collName.c_str(), JustWarning, "GetStatMap failed to locate requested StatMap.");
  return nullptr;
}