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
//
// G4PSCellFlux
#include "G4PSCellFluxToEffDose.hh"

#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring whole body effective dose (AP)
//   using the ICRP 116 conversion factors.
//
//   If you want to score only tracks passing through the geometry volume,
//  please use G4PSPassageCellFlux.
//
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
// 2010-07-22   Add weighted option
// 
///////////////////////////////////////////////////////////////////////////////

G4PSCellFluxToEffDose::G4PSCellFluxToEffDose(G4String name, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0),weighted(true), fEnergies(0.), fFlu2Dose(0.)
{
    DefineUnitAndCategory();
    SetUnit("Rem");
    //verboseLevel = 10;
}

G4PSCellFluxToEffDose::G4PSCellFluxToEffDose(G4String name, const G4String& unit, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0),weighted(true), fEnergies(0.), fFlu2Dose(0.)
{
    DefineUnitAndCategory();
    SetUnit(unit);
}

G4PSCellFluxToEffDose::~G4PSCellFluxToEffDose()
{;}

G4bool G4PSCellFluxToEffDose::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double stepLength = aStep->GetStepLength();
  if ( stepLength == 0. ) return FALSE;

  G4int idx = ((G4TouchableHistory*)
	       (aStep->GetPreStepPoint()->GetTouchable()))
               ->GetReplicaNumber(indexDepth);
  G4double cubicVolume = ComputeVolume(aStep, idx);

  G4double CellFlux = stepLength / cubicVolume;
  if (weighted) CellFlux *= aStep->GetPreStepPoint()->GetWeight(); 

  // Convert to Dose:
  G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
  G4double weight = Flu2Dose(kineticEnergy);
  CellFlux *= weight;
  G4int index = GetIndex(aStep);
  EvtMap->add(index,CellFlux);

  return TRUE;
}

void G4PSCellFluxToEffDose::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),
				    GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,EvtMap);
 
  // Effective dose data from ICRP 116 using AP orientation:
  std::vector<G4double> nEnerg, nFlu2Dose;
  std::vector<G4double> gEnerg, gFlu2Dose;
  nEnerg = {1.00E-09, 1.00E-08, 2.50E-08, 1.00E-07, 2.00E-07, 5.00E-07, 1.00E-06, 2.00E-06, 5.00E-06, 1.00E-05, 2.00E-05, 5.00E-05, 1.00E-04, 2.00E-04, 5.00E-04, 1.00E-03, 2.00E-03, 5.00E-03, 1.00E-02, 2.00E-02, 3.00E-02, 5.00E-02, 7.00E-02, 1.00E-01, 1.50E-01, 2.00E-01, 3.00E-01, 5.00E-01, 7.00E-01, 9.00E-01, 1.00E+00, 1.20E+00, 1.50E+00, 2.00E+00, 3.00E+00, 4.00E+00, 5.00E+00, 6.00E+00, 7.00E+00, 8.00E+00, 9.00E+00, 1.00E+01, 1.20E+01, 1.40E+01, 1.50E+01, 1.60E+01, 1.80E+01, 2.00E+01};
  nFlu2Dose = {3.09E-10, 3.55E-10, 4.E-10, 5.2E-10, 5.87E-10, 6.59E-10, 7.03E-10, 7.39E-10, 7.71E-10, 7.82E-10, 7.84E-10, 7.82E-10, 7.79E-10, 7.73E-10, 7.54E-10, 7.54E-10, 7.61E-10, 7.97E-10, 9.11E-10, 1.22E-09, 1.57E-09, 2.3E-09, 3.06E-09, 4.19E-09, 6.06E-09, 7.88E-09, 1.14E-08, 1.77E-08, 2.32E-08, 2.79E-08, 3.01E-08, 0.000000033, 3.65E-08, 4.07E-08, 4.58E-08, 4.83E-08, 4.94E-08, 4.98E-08, 4.99E-08, 4.99E-08, 0.00000005, 0.00000005, 4.99E-08, 4.95E-08, 4.93E-08, 0.000000049, 4.84E-08 , 4.77E-08};
  gEnerg = {0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.511, 0.6, 0.662, 0.8, 1., 1.117, 1.33, 1.5, 2., 3., 4., 5., 6., 6.129, 8., 10., 15., 20.};
  gFlu2Dose = {6.85E-12, 1.56E-11, 2.25E-11, 3.13E-11, 3.51E-11, 3.7E-11, 3.9E-11, 4.13E-11, 4.44E-11, 5.19E-11, 7.48E-11, 1.E-10, 1.51E-10, 2.E-10, 2.47E-10, 2.52E-10, 2.91E-10, 3.17E-10, 3.73E-10, 4.49E-10, 4.9E-10, 5.59E-10, 6.12E-10, 7.48E-10, 9.75E-10, 1.17E-09, 1.34E-09, 1.5E-09, 1.51E-09, 1.78E-09, 2.05E-09, 2.61E-09, 3.08E-09};

  // Get filter type for conversion factors:
  auto filter = GetFilter();
  if (filter == nullptr) {
      G4ExceptionDescription ED;
      ED << "No Particle Filter Assigned : \n" << "Neutron conversion factors being used." << G4endl;
      G4Exception("G4PSCellFluxToEffDose::Initialize","DetPS0001",JustWarning,ED);
      fEnergies = nEnerg; fFlu2Dose = nFlu2Dose;
      return;
  }
  G4String filterName = filter->GetName();
  if (filterName == "n-Filter") {
    fEnergies = nEnerg; fFlu2Dose = nFlu2Dose;
  } else if (filterName == "g-Filter") {
    fEnergies = gEnerg; fFlu2Dose = gFlu2Dose;
  } else {
      G4ExceptionDescription ED;
      ED << "No Particle Filter Assigned : \n" << "Neutron conversion factors being used." << G4endl;
      G4Exception("G4PSCellFluxToEffDose::Initialize","DetPS0001",JustWarning,ED);
      fEnergies = nEnerg; fFlu2Dose = nFlu2Dose;
  }
}

void G4PSCellFluxToEffDose::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSCellFluxToEffDose::clear(){
  EvtMap->clear();
}

void G4PSCellFluxToEffDose::DrawAll()
{;}

void G4PSCellFluxToEffDose::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  Dose : " << *(itr->second)/GetUnitValue() 
	   << " [" << GetUnit() << "]"
	   << G4endl;
  }
}

void G4PSCellFluxToEffDose::SetUnit(const G4String& unit)
{
    CheckAndSetUnit(unit,"Effective Dose");
}

void G4PSCellFluxToEffDose::DefineUnitAndCategory(){
   // Per Unit Surface
   new G4UnitDefinition("percentimeter2","percm2","Per Unit Surface",(1./cm2));
   new G4UnitDefinition("permillimeter2","permm2","Per Unit Surface",(1./mm2));
   new G4UnitDefinition("permeter2","perm2","Per Unit Surface",(1./m2));

  // Effective Dose Scorer
  new G4UnitDefinition("Sievert","Sv","Equivalent Dose",(gray));
  new G4UnitDefinition("milliSievert","mSv","Equivalent Dose",(1000*gray));
  new G4UnitDefinition("Rem","Rem","Effective Dose",(gray*100));
  new G4UnitDefinition("mRem","mRem","Effective Dose",(gray*100000));
  
}

G4double G4PSCellFluxToEffDose::ComputeVolume(G4Step* aStep, G4int idx){

  G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid = 0;
  if(physParam)
  { // for parameterized volume
    if(idx<0)
    {
      G4ExceptionDescription ED;
      ED << "Incorrect replica number --- GetReplicaNumber : " << idx << G4endl;
      G4Exception("G4PSCellFluxToEffDose::ComputeVolume","DetPS0001",JustWarning,ED);
    }
    solid = physParam->ComputeSolid(idx, physVol);
    solid->ComputeDimensions(physParam,idx,physVol);
  }
  else
  { // for ordinary volume
    solid = physVol->GetLogicalVolume()->GetSolid();
  }
  
  return solid->GetCubicVolume();
}

G4double G4PSCellFluxToEffDose::Flu2Dose(G4double kinE) {

  G4double weight = 0.;
  // Log-log interpolation between points:
  G4double energ1 = 0., energ2 = 0., f2d1 = 0., f2d2 = 0.;
  for (G4int i=0; i< (G4int) fEnergies.size() - 1; i++) {
    if (fEnergies[i] <= kinE && kinE < fEnergies[i+1]){
      energ1 = fEnergies[i]; energ2 = fEnergies[i+1];
      f2d1 = fFlu2Dose[i]; f2d2 = fFlu2Dose[i+1];
      break;
    } else { 
      energ1 = 0.; energ2 = 0.;
      f2d1 = 0.; f2d2 = 0.;
    }
  }

  if (energ1 == 0. || energ2 == 0. || f2d1 == 0. || f2d2 == 0.) {
		/*
      G4ExceptionDescription ED;
      ED << " No conversion data for particle energy : " << kinE <<
       "\n Zero dose computed." << G4endl;
      G4Exception("G4PSCellFluxToEffDose::Flu2Dose","DetPS0001",JustWarning,ED);
		 */
      return weight;
  } 

  // Fitting power to points of form y = bx^m (straight line for log-log plot)
  G4double M;
  M = (log10(f2d2) - log10(f2d1))/(log10(energ1) - log10(energ2));
  G4double b;
  b = f2d1/(pow(energ1, M));

  weight = b*pow(kinE, M)*GetUnitValue()*cm2;
  return weight;
}
