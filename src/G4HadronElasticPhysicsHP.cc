//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronElasticPhysicsHP
//
// Author: 3 June 2010 V. Ivanchenko
//
// Modified:
// 03.06.2011 V.Ivanchenko change design - now first default constructor 
//            is called, HP model and cross section are added on top
//
//----------------------------------------------------------------------------
//
// HP model for n with E < 20 MeV

#include "G4HadronElasticPhysicsHP.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPThermalScattering.hh"
#include "G4NeutronHPThermalScatteringData.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronElasticPhysicsHP);

G4HadronElasticPhysicsHP::G4HadronElasticPhysicsHP(G4int ver)
  : G4HadronElasticPhysics(ver, "hElasticWEL_CHIPS_HP")
{
  if(verbose > 1) { 
    G4cout << "### G4HadronElasticPhysicsHP: " << GetPhysicsName() 
	   << G4endl; 
  }
}

G4HadronElasticPhysicsHP::~G4HadronElasticPhysicsHP()
{}

void G4HadronElasticPhysicsHP::ConstructProcess()
{
  G4HadronElasticPhysics::ConstructProcess();

  // From
 // https://indico.cern.ch/event/245281/contributions/1564676/attachments/420136/583408/thermal_physics_validation_argarcia.pdf
  // and  
  // Bulk material interrogation experimental results and validation with
  // Geant4 for replacement of dangerous radiological sources in oil-well
  // logging industries.  Sharma et. al. 

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronElasticProcess* elasticnp = new G4HadronElasticProcess();
  G4NeutronHPElastic* elasticn = new G4NeutronHPElastic();
  elasticn->SetMinEnergy(4.*eV);
  G4NeutronHPThermalScattering* thermaln = new G4NeutronHPThermalScattering;
  thermaln->SetMaxEnergy(4.*eV);
  elasticnp->RegisterMe(thermaln);

  G4NeutronHPThermalScatteringData* hpthermaldata = new G4NeutronHPThermalScatteringData;

  elasticnp->AddDataSet(hpthermaldata);

  if(verbose > 1) {
    G4cout << "### HadronElasticPhysicsHP is constructed " 
	   << G4endl;
  }
}


