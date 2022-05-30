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

#ifndef G4PSCellFluxToEffDose_h
#define G4PSCellFluxToEffDose_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

class G4VSolid;

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring cell flux.
//   The Cell Flux is defined by  a sum of track length divided 
//   by the geometry volume, where all of the tracks in the geometry 
//   are taken into account. e.g. the unit of Cell Flux is mm/mm3.
//    
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


class G4PSCellFluxToEffDose : public G4VPrimitiveScorer
{
   public: // with description
      G4PSCellFluxToEffDose(G4String name, G4int depth=0);
      G4PSCellFluxToEffDose(G4String name, const G4String& unit, G4int depth=0);
      virtual ~G4PSCellFluxToEffDose();

      inline void Weighted(G4bool flg=true) { weighted = flg; }
      // Multiply track weight

  protected: // with description
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);

      virtual G4double ComputeVolume(G4Step*, G4int idx);

  public: 
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);
      virtual void clear();
      virtual void DrawAll();
      virtual void PrintAll();
      virtual void SetUnit(const G4String& unit);    
      inline G4double Flu2Dose(G4double kinE);

  protected:
      virtual void DefineUnitAndCategory();

  private:
      G4int HCID;
      G4THitsMap<G4double>* EvtMap;
      G4bool  weighted;
      std::vector<G4double> fEnergies;
      std::vector<G4double> fFlu2Dose;

};
#endif
