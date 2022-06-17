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

Analysis::Analysis(G4bool he3)
{
  eDepHist = 0;
  eDepHist1 = 0;
  eDepHist2 = 0;
  eDepHistTot = 0;
  primEneHist = 0;
  primPosHist = 0;
  convergenceName = "";
  isHe3 = he3;
}

//
//

Analysis::~Analysis() 
{
}

//
//

Analysis* Analysis::GetAnalysis(G4bool he3)
{
  if (!theAnalysis) {
    theAnalysis = new Analysis(he3);
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

  if (isHe3) {
    eDepHist = man->CreateH1("He3EnergyDep", "He3EnergyDep", 512, 0., 5.);
  } else {
    eDepHist1 = man->CreateH1("BF3EnergyDep1", "BF3EnergyDep1", 512, 0., 5.);
    eDepHist2 = man->CreateH1("BF3EnergyDep2", "BF3EnergyDep2", 512, 0., 5.);
    eDepHistTot = man->CreateH1("BF3EnergyDepTot", "BF3EnergyDepTot", 512, 0., 5.);
  }
  const std::vector<G4double> binEdges = {1.00E-11,1.00E-07,4.14E-07,8.76E-07,1.86E-06,5.04E-06,1.07E-05,3.73E-05,1.01E-04,2.14E-04,4.54E-04,1.58E-03,3.35E-03,7.10E-03,1.50E-02,2.19E-02,2.42E-02,3.18E-02,4.09E-02,6.74E-02,1.11E-01,1.83E-01,2.97E-01,3.69E-01,4.98E-01,6.08E-01,7.43E-01,8.23E-01,1.00E+00,1.35E+00,1.65E+00,1.92E+00,2.23E+00,2.35E+00,2.37E+00,2.47E+00,2.73E+00,3.01E+00,3.68E+00,4.97E+00,6.07E+00,7.41E+00,8.61E+00,1.00E+01,1.22E+01,1.42E+01,1.73E+01};
  primEneHist = man->CreateH1("PrimaryEnergy", "PrimaryEnergy", binEdges);
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

void Analysis::FillEDep(G4double eDep, G4int hist)
{
  //G4cout << "Adding Energy Deposittion. " << G4endl;+
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  if (hist == 0) {
  man->FillH1(eDepHist, eDep);
  } else if (hist == 1) {
    man->FillH1(eDepHist1, eDep);
  } else if (hist == 2) {
    man->FillH1(eDepHist2, eDep);
  } else if (hist == 3) {
    man->FillH1(eDepHistTot, eDep);
  }
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