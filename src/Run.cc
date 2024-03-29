// Source code for class Run().
// Created by Codey Olson on June 2, 2021.

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
#include "G4ThreeVector.hh"
#include "G4WorkerThread.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4StatAnalysis.hh"

#include <atomic>

Run::Run(G4bool he3)
{
  isHe3 = he3;
}

//
//

Run::~Run()
{}

//
//

void Run::Merge(const G4Run* aRun)
{
  G4Run::Merge(aRun);


}

//
//

void Run::RecordEvent(const G4Event* anEvent)
{
  G4int eventNum = anEvent->GetEventID();
  if (eventNum % 100000 == 0) {
    std::cout << "Event " << eventNum << " started." << std::endl;
  }

  Analysis* myAnalysis = Analysis::GetAnalysis(isHe3);
  G4SDManager* sdMan = G4SDManager::GetSDMpointer();
  // Get Primary Energy information:
  G4PrimaryVertex* pVertex = anEvent->GetPrimaryVertex();
  G4ThreeVector primPos = pVertex->GetPosition();
  G4PrimaryParticle* primary = pVertex->GetPrimary();
  if (primary->GetG4code()->GetParticleName() == "neutron") {
    G4double primEnergy = primary->GetKineticEnergy();
    myAnalysis->FillPrimaryEne(primEnergy/MeV);
    myAnalysis->FillPrimaryPos(primPos.getX()/cm, primPos.getY()/cm);
  }
  if (isHe3) {
    //G4cout << "Primary Energy is: " << energy/MeV << G4endl;
    G4HCofThisEvent* hce = anEvent->GetHCofThisEvent();
    G4int collID = sdMan->GetCollectionID("Helium-3/EnergyDep");
    if (!hce) return;
    G4THitsMap<G4double>* eventMap = 0;
    eventMap = static_cast<G4THitsMap<G4double>*>(hce->GetHC(collID));
    if (eventMap && eventMap->entries() >= 1) {
      G4double val = 0.;
      for (auto itr = eventMap->begin(); itr != eventMap->end(); itr++) {
        val += *itr->second;
      }
      if (val > 0.) {
        myAnalysis->FillEDep(val/MeV, 0);
      }
    }
  } else {
    //G4cout << "Primary Energy is: " << energy/MeV << G4endl;
    G4HCofThisEvent* hce = anEvent->GetHCofThisEvent();

    G4int collID1 = sdMan->GetCollectionID("BF31/EnergyDep1");
    if (!hce) return;
    G4THitsMap<G4double>* eventMap1 = 0;
    eventMap1 = static_cast<G4THitsMap<G4double>*>(hce->GetHC(collID1));
    if ((eventMap1 && eventMap1->entries() >= 1)) {
      G4double val1 = 0.;
      for (auto itr = eventMap1->begin(); itr != eventMap1->end(); itr++) {
        val1 += *itr->second;
      }
      //G4cout << "Detector 1: " << val1/MeV << G4endl;
      if ((val1) > 0./MeV) {
        myAnalysis->FillEDep((val1)/MeV, 1);
        myAnalysis->FillEDep((val1)/MeV, 3);
      }
    }
    G4int collID2 = sdMan->GetCollectionID("BF32/EnergyDep2");
    if (!hce) return;
    G4THitsMap<G4double>* eventMap2 = 0;
    eventMap2 = static_cast<G4THitsMap<G4double>*>(hce->GetHC(collID2));
    if ((eventMap2 && eventMap2->entries() >= 1)) {
      G4double val2 = 0.;
      for (auto itr = eventMap2->begin(); itr != eventMap2->end(); itr++) {
        val2 += *itr->second;
      }
      //G4cout << "Detector 2: " << val2/MeV << G4endl;
      if ((val2) > 0./MeV) {
        myAnalysis->FillEDep((val2)/MeV, 2);
        myAnalysis->FillEDep((val2)/MeV, 3);
      }
    }
  }
  G4Run::RecordEvent(anEvent);
}