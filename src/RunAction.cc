// Source file for RunAction() class.
// Created by Codey Olson on May 10, 2021.

/// \file RunAction.cc
/// \brief Source code for class RunAction.

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "Analysis.hh"
#include <G4WorkerThread.hh>

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSSphereSurfaceFlux.hh"
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4StatAnalysis.hh"
#include "G4RunManager.hh"
#include <iostream>
#include <fstream>

RunAction::RunAction()
{
  // Connect to pointers:
  fsdPointer = G4SDManager::GetSDMpointer();
  fanalysisManager = G4AnalysisManager::Instance();
  // Neutron energy bins:
  nEdges = {1.0E-07, 4.41E-07, 8.76E-07, 1.86E-06, 5.04E-06, 1.07E-05, 3.73E-05, 1.01E-04, 2.14E-04, 4.54E-04, 1.58E-03, 3.35E-03, 7.10E-03, 1.50E-02, 2.19E-02, 2.42E-02, 2.61E-02, 3.18E-02, 4.09E-02, 6.74E-02, 1.11E-01, 1.83E-01, 2.97E-01, 3.69E-01, 4.98E-01, 6.08E-01, 7.43E-01, 8.21E-01, 1, 1.35, 1.65, 1.92, 2.23, 2.35, 2.37, 2.47, 2.73, 3.01, 3.68, 4.97, 6.07, 7.41, 8.61, 10., 12.2, 14.2, 17.3};
  
  fName = "reactorBay_Storage_400M";
}

//
//

RunAction::~RunAction()
{
  delete fanalysisManager;
}

//
//

G4Run* RunAction::GenerateRun()
{
  return new Run;
}

//
//

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " started." << G4endl;
  const Run* myRun = static_cast<const Run*>(aRun);
  std::vector<G4String> detNames = myRun->GetDetNames();
  std::vector<std::vector<G4String>> collNames = myRun->GetCollNames();
  std::vector<std::vector<G4String>> unitStrings = myRun->GetUnitStrings();
  // Save random number seed in the runManager:
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  fanalysisManager->SetNtupleMerging(true);
  fanalysisManager->SetVerboseLevel(1);
  fanalysisManager->SetFileName(fName); // This will be a root file:
  fanalysisManager->SetFirstHistoId(0);
  fanalysisManager->CreateH1("nEnergy", "Primary Neutron Energies", nEdges);
  // Create ntuple for each detector:
  for (G4int i = 0; i < (G4int)detNames.size(); i ++) {
    // Create column for each primitive on the ith detector:
    fanalysisManager->CreateNtuple(detNames[i], "Bay");
    for (G4int j = 0; j <(G4int)collNames[i].size(); j++) {
      // Because the collection name is in the form
      // "detName/primName", I make the ntuple the name of the detecotr
      // with the columns the names of the primitives. The erase command
      // strips the detector portion from this and leaves the primitve.
      std::string collName = (std::string)collNames[i][j];
      std::string delim = "/";
      collName.erase(0, collName.find(delim)+1);
      fanalysisManager->CreateNtupleDColumn(collName);
    }
    fanalysisManager->FinishNtuple();
  }
  fanalysisManager->OpenFile();

}

//
//

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  const Run* myRun = static_cast<const Run*>(aRun);
  G4int nEvents = myRun->GetNumberOfEvent();
  if (nEvents == 0) return;
  // Check for local run (i.e. thread)
  if (!IsMaster()) {
    G4cout << "---> End of Local Run. Total Events processed: " << nEvents << G4endl;
  } else {
    G4cout << "---> End of Global Run. Total Events Processed: " << nEvents << G4endl;
    // Open output file:
    std::ofstream fileStream(fName+".txt");
    // Detector, primitive and unit vector:
    std::vector<G4String> detNames = myRun->GetDetNames();
    std::vector<std::vector<G4String>> collNames = myRun->GetCollNames();
    std::vector<std::vector<G4double>> unitVals = myRun->GetUnitVals();
    std::vector<std::vector<G4String>> unitStrings = myRun->GetUnitStrings();
    // Loop through detectors:
    for (G4int i = 0; i < (G4int)detNames.size(); i++) {
      // Loop through primitives on this detector:
      for (G4int j = 0; j < (G4int)collNames[i].size(); j++) {
        // Get units:
        G4double unitVal = unitVals[i][j];
        G4String unitString = unitStrings[i][j];
        // Get hit maps:
        G4THitsMap<G4double>* hitMap = myRun->GetHitsMap(collNames[i][j]);
        G4StatContainer<G4StatAnalysis>* statMap = myRun->GetStatMap(collNames[i][j]);
        // Iterate through hitMap:
        if (!hitMap || hitMap->size() == 0) continue;
        for (auto itr = hitMap->begin(); itr != hitMap->end(); itr++) {
          // Check for entries:
          if (!hitMap->GetObject(itr)) continue;
          fileStream << collNames[i][j] << ": " << itr->first << ", " << *itr->second/unitVal;
          fileStream << ", " << unitString << std::endl;
        }
        // Iterate through statmap:
        for (auto itr = statMap->begin(); itr != statMap->end(); itr++) {
          // check for enetries:
          if (!statMap->GetObject(itr)) continue;
          G4StatAnalysis* stats = statMap->GetObject(itr);
          auto stat = (*stats);
          fileStream << stat << std::endl;
        }
      }
    }
    fileStream.close();
  }
  fanalysisManager->Write();
  fanalysisManager->CloseFile();
}

