// Source file for Analysis().
// Created by Codey Olson on May 30, 2022.

/// \file Analysis.cc
/// \brief Source code for Analysis class.

#include "G4AutoDelete.hh"
#include "G4SystemOfUnits.hh"
#include "Analysis.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4SDManager.hh"
#include "g4root.hh"
#include "G4RootAnalysisManager.hh"
#include "G4ConvergenceTester.hh"

G4ThreadLocal Analysis* theAnalysis = 0;

namespace {
  G4Mutex aMutex = G4MUTEX_INITIALIZER;
  G4ConvergenceTester* fConvTest = new G4ConvergenceTester("ConvTest");
}

Analysis::Analysis()
{
  eDepHist = 0;
  primEneHist = 0;
  primPosHist = 0;
  convergenceName = "";
}

//
//

Analysis::~Analysis() 
{
}

//
//

Analysis* Analysis::GetAnalysis()
{
  if (!theAnalysis) {
    theAnalysis = new Analysis();
    G4AutoDelete::Register(theAnalysis);
  }
  return theAnalysis;
}

//
//

void Analysis::Book(G4String runName)
{
  convergenceName = runName;
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->SetVerboseLevel(2);
  #ifdef G4MULTITHREADED
    man->SetNtupleMerging(true);
  #endif
  man->SetFirstNtupleId(0);
  man->SetFirstNtupleColumnId(0);

  eDepHist = man->CreateH1("He3EnergyDep", "He3EnergyDep", 512, 0., 5.);
  primEneHist = man->CreateH1("PrimaryEnergy", "PrimaryEnergy", 50, 0., 20.);
  primPosHist = man->CreateH2("PrimaryPosition", "PrimaryPosition", 110, -5.5, 5.5, 90, -4.5, 4.5);

  return; 
}

//
//

void Analysis::OpenFile(const G4String& fileName)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->OpenFile(fileName.c_str());

  return;
}

//
//

void Analysis::Save()
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();

  return;
}

//
//

void Analysis::Close(G4bool reset)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->CloseFile(reset);

  return;
}

//
//

void Analysis::FillEDep(G4double eDep)
{
  //G4cout << "Adding Energy Deposittion. " << G4endl;+
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(eDepHist, eDep);
  G4AutoLock l(&aMutex);
  fConvTest->AddScore(eDep);
  return;
}

//
//

void Analysis::FillPrimaryEne(G4double energy)
{ 
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(primEneHist, energy);
}

//
//

void Analysis::FillPrimaryPos(G4double xPos, G4double yPos)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH2(primPosHist, xPos, yPos);
}

//
//

void Analysis::CheckConvergence()
{
  std::ofstream convOutput;
  convOutput.open(convergenceName+"-conv.txt");
  fConvTest->ShowResult(convOutput);
  fConvTest->ShowHistory(convOutput);
  convOutput.close();

  return;
}