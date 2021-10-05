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
/// \file field/field05/include/F05DetectorConstruction.hh
/// \brief Definition of the F05DetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F05DetectorConstruction_h
#define F05DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "DetectorConfiguration.hh"

class G4Material;
 
class G4Box;
class G4Tubs;
class G4Trd;
class G4VSolid;

class G4LogicalVolume;
class G4VPhysicalVolume;

class F05Field;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F05DetectorConstruction : public G4VUserDetectorConstruction{
public:
  F05DetectorConstruction();
  virtual ~F05DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
private: 
  G4Material* fVacuum;
  G4Material* KrGas;
  G4Material* Air;
  G4Material* Silicon;
  G4Material* Scintillator;
  G4Material* Aluminum;
  G4Material* Copper;
  G4Material* Kapton;
  G4Material* lead;
  G4Material* Kr_HeGas;
  G4Material* FR4;
  G4Material* Permalloy;
  G4Material* TargetGas;

  G4double           fWorldSizeXY;
  G4double           fWorldSizeZ;

  G4Box*             fSolidWorld;    //pointer to the solid World
  G4LogicalVolume*   fLogicWorld;    //pointer to the logical World
  G4VPhysicalVolume* fPhysiWorld;    //pointer to the physical World

  G4Tubs* fSolidChamber;
  G4LogicalVolume* fLogicChamber;
  G4VPhysicalVolume* fPhysiChamber;

  G4Tubs* fSolidChamberFlangeU;
  G4LogicalVolume* fLogicChamberFlangeU;
  G4VPhysicalVolume* fPhysiChamberFlangeU;

  G4Tubs* fSolidChamberFlangeD;
  G4LogicalVolume* fLogicChamberFlangeD;
  G4VPhysicalVolume* fPhysiChamberFlangeD;

  G4Tubs* fSolidGas;
  G4LogicalVolume* fLogicGas;
  G4VPhysicalVolume* fPhysiGas;
  G4Tubs* fSolidGasC;
  G4LogicalVolume* fLogicGasC;
  G4VPhysicalVolume* fPhysiGasC;

  G4Tubs* fSolidCavity;
  G4LogicalVolume* fLogicCavity;
  G4VPhysicalVolume* fPhysiCavity;

  G4Tubs* fSolidChamberFoil;
  G4LogicalVolume* fLogicChamberFoil;
  G4VPhysicalVolume* fPhysiChamberFoil;

  G4Tubs* fSolidCavityFlange;
  G4LogicalVolume* fLogicCavityFlange;
  G4VPhysicalVolume* fPhysiCavityFlangeU;
  G4VPhysicalVolume* fPhysiCavityFlangeD;

  G4Tubs* fSolidCavityFoil;
  G4LogicalVolume* fLogicCavityFoil;
  G4VPhysicalVolume* fPhysiCavityFoilU;
  G4VPhysicalVolume* fPhysiCavityFoilD;

  G4Box* fSolidBPM;
  G4LogicalVolume* fLogicBPM;
  G4VPhysicalVolume* fPhysiBPMU;
  G4VPhysicalVolume* fPhysiBPMD;

  G4Tubs* fSolidKapton;
  G4LogicalVolume* fLogicKapton;
  G4VPhysicalVolume* fPhysiKapton;

  G4Tubs* fSolidVacuum;
  G4LogicalVolume* fLogicVacuum;
  G4VPhysicalVolume* fPhysiVacuum;
  G4Tubs* fSolidBlock;
  G4LogicalVolume* fLogicBlock;
  G4VPhysicalVolume* fPhysiBlock;

  G4Box* fSolidDetector;
  G4LogicalVolume* fLogicDetector;
  G4VPhysicalVolume* fPhysiDetectorU;
  G4VPhysicalVolume* fPhysiDetectorD;

  G4Box* soli_si;
  G4LogicalVolume* logi_si_xstrip;
  G4LogicalVolume* logi_si_ystrip;
  G4VPhysicalVolume* phys_siU;
  G4VPhysicalVolume* phys_siD;

  G4Box *SolidAlCouver;
  G4LogicalVolume* LogicAlCouver;
  G4VPhysicalVolume *PhysAlCouverU;
  G4VPhysicalVolume *PhysAlCouverD ;

  G4Box* div_soli_silicon;
  G4LogicalVolume* div_logi_silicon_xstrip;
  G4LogicalVolume* div_logi_silicon_ystrip;
  G4VPhysicalVolume* div_silicon_xstrip;
  G4VPhysicalVolume* div_silicon_ystrip;

  G4Box* SiliconCircuit;
  G4Box* MotherFrame;
  G4Box* BaseSubtraction;

  G4VSolid* SoliCircuit;
  G4LogicalVolume* LogicCircuit;
  G4VPhysicalVolume* PhysCircuitU;
  G4VPhysicalVolume* PhysCircuitD;
  G4VSolid* SoliMotherFrame;
  G4LogicalVolume* LogicMotherFrame;
  G4VPhysicalVolume* PhysMotherFrameU;
  G4VPhysicalVolume* PhysMotherFrameD;

  G4Box* BasePlate;
  G4Box* BasePlateFoot;
  G4VSolid* BasePlateWana;
  G4VSolid* SoliBasePlate;
  G4LogicalVolume* LogicBasePlate;
  G4VPhysicalVolume* PhysBasePlate;
  
  static G4ThreadLocal F05Field* fField;

private: 
  void DefineMaterials();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
