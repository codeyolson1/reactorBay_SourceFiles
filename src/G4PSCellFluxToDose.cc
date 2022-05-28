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
#include "G4PSCellFluxToDose.hh"

#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring cell flux the equivalent
//   dose using the cell flux scorere and ANSI conversion factors.
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

G4PSCellFluxToDose::G4PSCellFluxToDose(G4String name, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0),weighted(true), fEnergies(0.), fFlu2Dose(0.)
{
    DefineUnitAndCategory();
    SetUnit("RemEQ");
    //verboseLevel = 10;
}

G4PSCellFluxToDose::G4PSCellFluxToDose(G4String name, const G4String& unit, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0),weighted(true), fEnergies(0.), fFlu2Dose(0.)
{
    DefineUnitAndCategory();
    SetUnit(unit);
}

G4PSCellFluxToDose::~G4PSCellFluxToDose()
{;}

G4bool G4PSCellFluxToDose::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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

void G4PSCellFluxToDose::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),
				    GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,EvtMap);
 
  // Neutron Data from ANSI 6.1.1:
  std::vector<G4double> nEnerg, nFlu2Dose;
  nEnerg = {2.50E-08, 1.00E-07, 1.00E-06, 1.00E-05, 1.00E-04, 1.00E-03, 1.00E-02, 
  1.00E-01, 5.00E-01, 1.00E+00, 2.00E+00, 2.50E+00, 5.00E+00, 7.00E+00, 1.00E+01, 1.40E+01, 
  2.00E+01};
  nFlu2Dose = {1.02E-09, 1.02E-09, 1.24E-09, 1.26E-09, 1.16E-09, 1.04E-09, 9.89E-10, 
  6.03E-09, 2.57E-08, 3.67E-08, 3.57E-08, 3.47E-08, 4.33E-08, 4.08E-08, 4.08E-08, 5.78E-08, 
  6.31E-08};
  // Gamma Data from ANSI 6.1.1:
  std::vector<G4double> gEnerg, gFlu2Dose;
  gEnerg = {0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 2.8, 3.25, 3.75, 4.25, 4.75, 5.0, 5.25, 5.75, 6.25, 6.75, 7.5, 9.0, 11.0, 13.0, 15.0};
  gFlu2Dose = {1.10E-09, 1.62E-10, 8.06E-11, 7.17E-11, 7.86E-11, 1.05E-10, 1.39E-10, 1.75E-10, 2.11E-10, 2.44E-10, 2.74E-10, 3.00E-10, 3.25E-10, 3.53E-10, 3.78E-10, 4.00E-10, 4.22E-10, 4.67E-10, 5.50E-10, 6.97E-10, 8.31E-10, 9.50E-10, 1.06E-09, 1.11E-09, 1.23E-09, 1.34E-09, 1.45E-09, 1.56E-09, 1.61E-09, 1.67E-09, 1.77E-09, 1.87E-09, 1.98E-09, 2.13E-09, 2.44E-09, 2.86E-09, 3.28E-09, 3.69E-09};

  // Get filter type for conversion factors:
  auto filter = GetFilter();
  if (filter == nullptr) {
      G4ExceptionDescription ED;
      ED << "No Particle Filter Assigned : \n" << "Neutron conversion factors being used." << G4endl;
      G4Exception("G4PSCellFluxToDose::Initialize","DetPS0001",JustWarning,ED);
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
      G4Exception("G4PSCellFluxToDose::Initialize","DetPS0001",JustWarning,ED);
      fEnergies = nEnerg; fFlu2Dose = nFlu2Dose;
  }
}

void G4PSCellFluxToDose::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSCellFluxToDose::clear(){
  EvtMap->clear();
}

void G4PSCellFluxToDose::DrawAll()
{;}

void G4PSCellFluxToDose::PrintAll()
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

void G4PSCellFluxToDose::SetUnit(const G4String& unit)
{
    CheckAndSetUnit(unit,"Equivalent Dose");
}

void G4PSCellFluxToDose::DefineUnitAndCategory(){
   // Per Unit Surface
   new G4UnitDefinition("percentimeter2","percm2","Per Unit Surface",(1./cm2));
   new G4UnitDefinition("permillimeter2","permm2","Per Unit Surface",(1./mm2));
   new G4UnitDefinition("permeter2","perm2","Per Unit Surface",(1./m2));

  // Equivelent Dose Scorer
    new G4UnitDefinition("Sievert","Sv","Equivalent Dose",(gray));
    new G4UnitDefinition("milliSievert","mSv","Equivalent Dose",(1000*gray));
    new G4UnitDefinition("RemEQ","RemEQ","Equivalent Dose",(gray*100));
    new G4UnitDefinition("mRemEQ","mRemEQ","Equivalent Dose",(gray*100000));
}

G4double G4PSCellFluxToDose::ComputeVolume(G4Step* aStep, G4int idx){

  G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid = 0;
  if(physParam)
  { // for parameterized volume
    if(idx<0)
    {
      G4ExceptionDescription ED;
      ED << "Incorrect replica number --- GetReplicaNumber : " << idx << G4endl;
      G4Exception("G4PSCellFluxToDose::ComputeVolume","DetPS0001",JustWarning,ED);
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

G4double G4PSCellFluxToDose::Flu2Dose(G4double kinE) {

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
      G4Exception("G4PSCellFluxToDose::Flu2Dose","DetPS0001",JustWarning,ED);
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
