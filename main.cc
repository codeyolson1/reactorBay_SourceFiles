// Main program for Pu-Be Source in MEB 1205 E (Reactor Bay)
// Created by Codey Olson on May 7, 2021

/// \file
/// \brief // Main program for Pu-Be Source in MEB 1205 E

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else 
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "QGSP_BIC_HP.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"

int main(int argc, char** argv)
{
  G4UIExecutive* ui = 0;
  if ( argc == 1 ){
    ui = new G4UIExecutive(argc, argv);
  }

  //G4Random::setTheEngine(new CLHEP::MixMaxRng);
  G4MTRunManager* runManager = new G4MTRunManager;

  runManager->SetUserInitialization(new DetectorConstruction());

  G4VModularPhysicsList* physicsList = new QGSP_BIC_HP();
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
  runManager->SetVerboseLevel(0);

  runManager->SetUserInitialization(new ActionInitialization());

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( !ui ) {
    // batch mode - Apply macros directly 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  } else {
    // Interactive Mode:
    UImanager->ApplyCommand("/control/macroPath /reactorBay_sourceFiles");
    UImanager->ApplyCommand("/control/execute init.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;

}