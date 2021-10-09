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
/// \file field/field05/src/F05DetectorConstruction.cc
/// \brief Implementation of the F05DetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

// for silicon detector
#include "G4PVDivision.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"

#include "F05Field.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

//#include "G4RepleteEofM.hh"
#include "G4EqEMFieldWithSpin.hh"

#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"

#include "DetectorConfiguration.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05DetectorConstruction::F05DetectorConstruction()
 : fVacuum(0), fWorldSizeXY(0), fWorldSizeZ(0), 
   fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0)
{
  // materials
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05DetectorConstruction::~F05DetectorConstruction()
{
  if (fField) delete fField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05DetectorConstruction::DefineMaterials()
{
  G4NistManager* nistMan = G4NistManager::Instance();

  fVacuum = nistMan->FindOrBuildMaterial("G4_Galactic");
  Air = nistMan->FindOrBuildMaterial("G4_AIR");
  Silicon = nistMan->FindOrBuildMaterial("G4_Si");
  Kapton = nistMan->FindOrBuildMaterial("G4_KAPTON");
  Aluminum = nistMan->FindOrBuildMaterial("G4_Al");
  Copper = nistMan->FindOrBuildMaterial("G4_Cu");
  lead = nistMan->FindOrBuildMaterial("G4_Pb");
  KrGas = nistMan->ConstructNewGasMaterial("KrGas","G4_Kr",
                                          300*kelvin,//temperature                                       
                                          1*atmosphere);//gas pressure;
  
  std::vector<G4String> v_sci_components={"C","H"};
  std::vector<double> v_sci_ratio={9.*nistMan->GetAtomicMassAmu(6),
                                   10.*nistMan->GetAtomicMassAmu(1)};
  Scintillator=nistMan->ConstructNewMaterial("Scintillator",v_sci_components,
                                             v_sci_ratio,
                                             1.032*g/cm3);
  Scintillator->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  std::vector<G4String> v_FR4_components = {"Si","O","C","H"};
  std::vector<double> v_FR4_ratio = {1.*nistMan->GetAtomicMassAmu(14),
				     2.*nistMan->GetAtomicMassAmu(8),
				     3.*nistMan->GetAtomicMassAmu(6),
				     3.*nistMan->GetAtomicMassAmu(1)};
  FR4 = nistMan->ConstructNewMaterial("FR4",
				      v_FR4_components,
				      v_FR4_ratio,
				      1.700*g/cm3);

  std::ofstream fout_material;
  fout_material.open("material.dat",std::ios::out);
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  fout_material << *(G4Material::GetMaterialTable()) << std::endl; 
  fout_material.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* F05DetectorConstruction::Construct(){
  //----------------------------------
  // World
  //----------------------------------
  fWorldSizeXY = 1.2*m;
  fWorldSizeZ = 3.6*m;
  fSolidWorld = new G4Box("World",//its name
			  fWorldSizeXY/2,fWorldSizeXY/2,fWorldSizeZ/2);//its size
 
  fLogicWorld = new G4LogicalVolume(fSolidWorld,//its solid
                                    fVacuum,//its material, fVacuum,
                                    "World");//its name
 
  fPhysiWorld = new G4PVPlacement(0,//no rotation
                                  G4ThreeVector(),//at (0,0,0)
                                  fLogicWorld,//its logical volume
                                  "World",//its name
                                  0,//its mother  volume
                                  false,//no boolean operation
                                  0);//copy number
  
  //---------------------------------                                                                   
  // Chamber                                                                                            
  //---------------------------------                                                                   
  fSolidChamber
    =new G4Tubs("Chamber",//its name                                                                    
                 chamber_diameter*0.5*mm,//inside radius                                                 
                 (chamber_diameter*0.5+chamber_thickness)*mm,//outside radius                            
                 chamber_length*0.5*mm,//length                                                          
                 0*deg,360*deg);//size                                                                   

  fLogicChamber
    =new G4LogicalVolume(fSolidChamber,//its solid                                                      
                          Aluminum,//its material                                                        
                          "Chamber");//its name                                                          

  fPhysiChamber
    =new G4PVPlacement(0,//no rotation                                                                  
                        G4ThreeVector(0,0,Magnet_center*mm),
                        fLogicChamber,//its fLogical volume                                              
                        "Chamber",//its name                                                             
                        fLogicWorld,//its mother volume                                                  
                        false,//no boolean operation                                                     
                        0);//copy number
  
  //---------------------------------                                                                  
  // ChamberFlange                                                                                      
  //---------------------------------                                                                   
  //--- Up stream ---                                                                                   
  fSolidChamberFlangeU
    = new G4Tubs("ChamberFlangeU",//its name                                                             
                 chamber_window_diameter*0.5*mm,
                 chamber_flange_diameter*0.5*mm,
                 chamber_flange_thickness_u*0.5*mm,
                 0*deg,360*deg);//size                                                                   

  fLogicChamberFlangeU
    = new G4LogicalVolume(fSolidChamberFlangeU,//its solid                                               
                          Aluminum,//its material                                                        
                          "ChamberFlangeU");//its name                                                   
  
  fPhysiChamberFlangeU
    = new G4PVPlacement(0,//no rotation                                                                  
                        G4ThreeVector(0,0,
                                      (Magnet_center-(chamber_length+chamber_flange_thickness_u)*0.5)*mm),
                        fLogicChamberFlangeU,//its fLogical volume                                       
                        "ChamberFlangeU",//its name                                                      
                        fLogicWorld,//its mother volume                                                  
                        false,//no boolean operation                                                     
                        0);//copy number
   
  //--- Down stream ---                                                                                 
  fSolidChamberFlangeD
    = new G4Tubs("ChamberFlangeD",//its name                                                             
                 0,
                 chamber_flange_diameter*0.5*mm,
                 chamber_flange_thickness_d*0.5*mm,
                 0*deg,360*deg);//size                                                                   
  
  // Aditional overhang                                                                                 
  G4VSolid* Overhang
    = new G4Tubs("ChamberFlangeD",
		 0,
		 overhang_diameter*0.5*mm,
		 overhang_thickness*0.5*mm,
		 0*deg,360*deg);//size

  G4RotationMatrix* prot = new G4RotationMatrix();
  prot->rotateZ(0*deg);
  
  G4VSolid* FlangeD
    = new G4UnionSolid("ChamberFlangeD",
		       fSolidChamberFlangeD,Overhang,prot,
		       G4ThreeVector(0,0,-(chamber_flange_thickness_d)*0.5*mm));
  
  fLogicChamberFlangeD = new G4LogicalVolume(FlangeD,Aluminum,"ChamberFlangeD");

  fPhysiChamberFlangeD
    = new G4PVPlacement(0,//no rotation                                                                  
			G4ThreeVector(0,0,
				      (Magnet_center+(chamber_length+chamber_flange_thickness_d)*0.5)*mm),
			fLogicChamberFlangeD,//its fLogical volume                                       
			"ChamberFlangeD",//its name                                                      
			fLogicWorld,//its mother volume                                                  
			false,//no boolean operation                                                     
			0);//copy number
  
  //------------------------                                                                        
  // ChamberFoil                                                                                    
  //------------------------                                                                        
  fSolidChamberFoil
    = new G4Tubs("ChamberFoil",
		 0*mm,
		 chamber_window_diameter*0.5*mm,
		 chamber_foil_thickness*0.5*mm,
		 0*deg,360*deg);

  fLogicChamberFoil = new G4LogicalVolume(fSolidChamberFoil,Aluminum,"ChamberFoil");

  fPhysiChamberFoil
    = new G4PVPlacement(0,//no rotation                                                              
			G4ThreeVector(0,0,
				      (Magnet_center-(chamber_length+chamber_foil_thickness)*0.5
				       -chamber_flange_thickness_u)*mm),
			fLogicChamberFoil,//its fLogical volume                                      
			"ChamberFoil",//its name                                                     
			fLogicWorld,//its mother volume                                              
			false,//no boolean operation                                                 
			0);//copy number
  
  //----------------------                                                                          
  // Gas Volume                                                                                     
  //----------------------                                                                          
  //--- Main volume ---                                                                             
  fSolidGasC
    =new G4Tubs("TargetGas",//its name                                                              
                0.0,
                chamber_diameter*0.5*mm,
                (chamber_length-overhang_thickness)*0.5*mm,
                0*deg,360*deg);//size                                                               
  fLogicGasC
    =new G4LogicalVolume(fSolidGasC,//its solid                                                     
                         KrGas,//its material                                                   
                         "TargetGas");//its name                                                    
  fPhysiGasC
    =new G4PVPlacement(0,//no rotation                                                              
                       G4ThreeVector(0,0,
                                     (Magnet_center-overhang_thickness*0.5)*mm), //Magnet_center-overhang_thickness*0.5
                       fLogicGasC,//its fLogical volume                                             
                       "TargetGas",//its name                                                       
                       fLogicWorld,//its mother volume                                              
                       false,//no boolean operation                                                 
                       0);//copy number
  
  //--- Upstream Flange volume ---                                                                  
  fSolidGas
    = new G4Tubs("TargetGas",//its name                                                              
		 0,
		 chamber_window_diameter*0.5*mm,
		 chamber_flange_thickness_u*0.5*mm,
		 0*deg,360*deg);//size                                                               

  fLogicGas
    = new G4LogicalVolume(fSolidGas,//its solid                                                      
			  KrGas,//its material                                                   
			  "TargetGas");//its name                                                    

  fPhysiGas
    = new G4PVPlacement(0,//no rotation                                                              
			G4ThreeVector(0,0,
				      (Magnet_center-(chamber_length+chamber_flange_thickness_u)*0.5)*mm),
			fLogicGas,//its fLogical volume                                              
			"TargetGas",//its name                                                       
			fLogicWorld,//its mother volume                                              
			false,//no boolean operation                                                 
			0);//copy number

  //-------------------------                                                                           
  // Cavity                                                                                             
  //-------------------------                                                                           
  fSolidCavity
    =new G4Tubs("Cavity",//its name                                                                     
                cavity_diameter*0.5*mm,
                (cavity_diameter*0.5+cavity_thickness)*mm,
                cavity_length*0.5*mm,
                0*deg,360*deg);//size                                                                   
  fLogicCavity
    =new G4LogicalVolume(fSolidCavity,//its solid                                                       
                         Copper,//its material                                                          
                         "Cavity");//its name                                                           
  fPhysiCavity
    =new G4PVPlacement(0,//no rotation                                                                  
                       G4ThreeVector(0,0,cavity_center*mm),
                       fLogicCavity,//its fLogical volume                                               
                       "Cavity",//its name                                                              
                       fLogicGasC,//its mother volume                                                   
                       false,//no boolean operation                                                     
                       0);//copy number
  
  //------------------------------                                                                  
  // Cavity Flange                                                                                  
  //------------------------------                                                                  
  fSolidCavityFlange
    =new G4Tubs("CavityFlange",//its name                                                           
                cavity_diameter*0.5*mm,
                (cavity_diameter*0.5+cavity_flange_thickness)*mm,
                cavity_flange_length*0.5*mm,
                0*deg,360*deg);//size                                                               

  fLogicCavityFlange
    =new G4LogicalVolume(fSolidCavityFlange,//its solid                                             
                         Copper,//its material                                                      
                         "CavityFlange");//its name                                                 

  //--- Up stream ---                                                                               
  fPhysiCavityFlangeU
    =new G4PVPlacement(0,//no rotation                                                              
                       G4ThreeVector(0,0,
                                     (cavity_center-(cavity_flange_length+cavity_length)*0.5)*mm),
                       fLogicCavityFlange,//its fLogical volume                                     
                       "CavityFlange",//its name                                                    
                       fLogicGasC,//its mother volume                                               
                       false,//no boolean operation                                                 
                       0);//copy number
  //--- Down stream ---                                                                             
  fPhysiCavityFlangeD
    =new G4PVPlacement(0,//no rotation                                                              
                       G4ThreeVector(0,0,
                                     (cavity_center+(cavity_flange_length+cavity_length)*0.5)*mm),
                       fLogicCavityFlange,//its fLogical volume                                     
                       "CavityFlange",//its name                                                    
                       fLogicGasC,//its mother volume                                               
                       false,//no boolean operation                                                 
                       1);//copy number
  
  //-------------------------------                                                                 
  // CavityFoil                                                                                     
  //-------------------------------
  fSolidCavityFoil                                                                                                          
    =new G4Tubs("CavityFoil",//its name                                                                              
                0*mm,                                                                                                            
                cavity_diameter*0.5*mm,                                                                              
                cavity_foil_thickness*0.5*mm,                                                                                              
                0*deg,360*deg);//size                                                                                                                                                 
                                                                                                                                                                                      
  fLogicCavityFoil                                                                                                                                                                    
    =new G4LogicalVolume(fSolidCavityFoil,//its solid                                                                                                                                 
                         Copper,//its material                                                                                                                                        
                         "CavityFoil");//its name
  //--- Up stream ---                                                                               
  fPhysiCavityFoilU
    =new G4PVPlacement(0,//no rotation                                                              
                       G4ThreeVector(0,0,
                                     (cavity_center-26.-cavity_length*0.5)*mm),
                       fLogicCavityFoil,//its fLogical volume                                       
                       "CavityFoil",//its name                                                      
                       fLogicGasC,//its mother volume                                               
                       false,//no boolean operation                                                 
                       0);//copy number                                                             
  //--- Down stream ---                                                                             
  fPhysiCavityFoilD
    =new G4PVPlacement(0,//no rotation                                                              
                       G4ThreeVector(0,0,(cavity_center+26.+cavity_length*0.5)*mm),
                       fLogicCavityFoil,//its fLogical volume                                       
                       "CavityFoil",//its name                                                      
                       fLogicGasC,//its mother volume                                               
                       false,//no boolean operation                                                 
                       1);//copy number
  
  
  //--------------------------------------                                                          
  // BeamProfileMonitor                                                                             
  //--------------------------------------                                                          
  fSolidBPM
    = new G4Box("BPM",
		bpm_sizeXY*0.5*mm,
		bpm_sizeXY*0.5*mm,
		bpm_thickness*0.5*mm);
  fLogicBPM
    = new G4LogicalVolume(fSolidBPM,
			  Scintillator,
			  "BPM");
  
  //--- Up stream ---                                                                               
  fPhysiBPMU
    = new G4PVPlacement(0,
			G4ThreeVector(0,0,(Magnet_center+bpm_positionU)*mm),
			fLogicBPM,
			"BPM_U",
			fLogicWorld,
			false,
			0);
  //--- Down stream ---                                                                             
  fPhysiBPMD
    = new G4PVPlacement(0,
			G4ThreeVector(0,0,(Magnet_center+bpm_positionD)*mm),
			fLogicBPM,
			"BPM_D",
			fLogicWorld,
			false,
			1);
  
  //--------------------------------                                                                
  // Kapton foil                                                                                    
  //--------------------------------                                                                
  fSolidKapton
    = new G4Tubs("Kapton",//its name                                                                 
		 0.0,
		 kapton_diameter*0.5*mm,
		 kapton_thickness*0.5*mm,
		 0*deg,360*deg);//size                                                               
  fLogicKapton
    = new G4LogicalVolume(fSolidKapton,//its solid                                                   
			  Kapton,//its material                                                      
			  "Kapton");//its name                                                       

  fPhysiKapton
    = new G4PVPlacement(0,//no rotation                                                              
			G4ThreeVector(0,0,(Magnet_center+kapton_center)*mm), // 700 mm
			fLogicKapton,//its fLogical volume                                           
			"Kapton",//its name                                                          
			fLogicWorld,//its mother volume                                              
			false,//no boolean operation                                                 
			0);//copy number

  //--------------------------------------                                                          
  // Positron Counter                                                                               
  //---------------------------------------                                                         
  //---- main scintillator ----                                                                     
  //const int n_positron_counter=2; 
  fSolidDetector
    =new G4Box("Detector",
               counter_sizeXY*0.5*mm,
               counter_sizeXY*0.5*mm,
               counter_thickness*0.5*mm);

  fLogicDetector=new G4LogicalVolume(fSolidDetector,Scintillator,"Detector");
  
  //---- Al couver ------                                                                           
  SolidAlCouver
    =new G4Box("AlCouv",
               counter_sizeXY*0.5*mm,
               counter_sizeXY*0.5*mm,
               counter_couver_thickness*0.5*mm);

  LogicAlCouver=new G4LogicalVolume(SolidAlCouver,Aluminum,"AlCouv");

  //--- Up stream ---                                                                               
  fPhysiDetectorU
    =new G4PVPlacement(0,//no rotation                                                              
                       G4ThreeVector(0,0,(Magnet_center+counter_centerU)*mm),
                       fLogicDetector,//its fLogical volume                                         
                       "Detector",//its name                                                        
                       fLogicWorld,//its mother volume                                              
                       false,//no boolean operation                                                 
                       0);//copy number                                                             

  PhysAlCouverU
    =new G4PVPlacement(0,//no rotation                                                              
                       G4ThreeVector(0,0,
                                     (Magnet_center+counter_centerU
                                     -(counter_thickness+counter_couver_thickness)*0.5)*mm),
                       LogicAlCouver,//its fLogical volume                                          
                       "AlCouv",//its name                                                          
                       fLogicWorld,//its mother volume                                              
                       false,//no boolean operation                                                 
                       0);

  //--- Down stream ---                                                                             
  fPhysiDetectorD
    =new G4PVPlacement(0,//no rotation                                                              
                       G4ThreeVector(0,0,(Magnet_center+counter_centerD)*mm),
                       fLogicDetector,//its fLogical volume                                         
                       "Detector",//its name                                                        
                       fLogicWorld,//its mother volume                                              
                       false,//no boolean operation                                                 
                       0);//copy number                                                             
  PhysAlCouverD
    =new G4PVPlacement(0,//no rotation                                                              
                       G4ThreeVector(0,0,
                                     (Magnet_center+counter_centerD
                                     -(counter_thickness+counter_couver_thickness)*0.5)*mm),
                       LogicAlCouver,//its fLogical volume                                          
                       "AlCouv",//its name                                                          
                       fLogicWorld,//its mother volume                                              
                       false,//no boolean operation                                                 
                       0);
  
  //--------------------------------
  // Vacuum region
  //--------------------------------
  double vacuum_z = (Magnet_center+kapton_center
                   +(fWorldSizeZ-kapton_thickness)*0.5)*mm; // 1200-500+(3600-0.075)*0.5=700+1800-0.0325=2499.9675
  fSolidVacuum
    = new G4Tubs("Vacuum",//its name 
                 0.0,
                 kapton_diameter*0.5*mm,
                 vacuum_z*0.5,
                 0*deg,360*deg);//size
  
  fLogicVacuum
    = new G4LogicalVolume(fSolidVacuum,//its solid
			  fVacuum,//its material
			  "Vacuum");//its name

  fPhysiVacuum
    = new G4PVPlacement(0,//no rotation
			G4ThreeVector(0,0,(vacuum_z-fWorldSizeZ)*0.5),//  1249.98375-1800=-550.01625,
			fLogicVacuum,//its fLogical volume
			"Vacuum",//its name
			fLogicWorld,//its mother volume
			false,//no boolean operation
			0);//copy number
  
  //--------------------------------
  //Block                                
  //--------------------------------                                                              
  double Block_z = vacuum_z;
  fSolidBlock
    = new G4Tubs("Block",//its name                                      
                 kapton_diameter*0.5*mm,
                 fWorldSizeXY*0.5,
                 Block_z*0.5,
                 0*deg,360*deg);//size
  
  fLogicBlock
    = new G4LogicalVolume(fSolidBlock,//its solid
                          lead,//its material
                          "Block");//its name
  
  fPhysiBlock
    = new G4PVPlacement(0,//no rotation
                        G4ThreeVector(0,0,(Block_z-fWorldSizeZ)*0.5),
                        fLogicBlock,//its fLogical volume
                        "Block",//its name
                        fLogicWorld,//its mother volume
                        false,//no boolean operation
                        0);//copy number

  //-----------------------------------
  // Silicon detector
  //  added by Shoichiro Nishimura
  //-----------------------------------                                                                                                                                                      
  soli_si
    = new G4Box("silicon",
                silicon_sizeXY*0.5*mm,
                silicon_sizeXY*0.5*mm,
                silicon_thickness*0.5*mm);
  
  logi_si_xstrip
    = new G4LogicalVolume(soli_si,
			  Silicon,
			  "silicon");
  
  logi_si_ystrip
    = new G4LogicalVolume(soli_si,
			  Silicon,
			  "silicon");

  div_soli_silicon
    = new G4Box("silicon",
                silicon_sizeXY*0.5*mm,
                silicon_sizeXY*0.5*mm,
                silicon_thickness*0.5*mm);
  
  div_logi_silicon_xstrip
    = new G4LogicalVolume(div_soli_silicon,
			  Silicon,
			  "divsilicon");
  
  div_logi_silicon_ystrip
    = new G4LogicalVolume(div_soli_silicon,
			  Silicon,
			  "divsilicon");

  // --- Up stream --- 
  phys_siU
    =new G4PVPlacement(0,
                       G4ThreeVector(0,0,(Magnet_center+silicon_centerU)*mm),
                       logi_si_xstrip,
                       "silicon",
                       fLogicWorld,
                       false,
                       2000);
  
  // --- divide to strips ---                                                                           
  div_silicon_xstrip
    =new G4PVDivision("silicon",div_logi_silicon_xstrip,logi_si_xstrip,
                      kXAxis,silicon_strip_num,
                      silicon_strip_pitch*mm,silicon_edge*mm);

  div_silicon_ystrip
    =new G4PVDivision("silicon",div_logi_silicon_ystrip,logi_si_ystrip,
                      kYAxis,silicon_strip_num,
                     silicon_strip_pitch*mm,silicon_edge*mm);

  //------------------------------------------  
  // Silicon circuit board and Mother Frame
  //------------------------------------------                                      
  SiliconCircuit
    = new G4Box("circuit",
                circuit_sizeX*0.5*mm,
                circuit_sizeY*0.5*mm,
                circuit_thickness*0.5*mm);
  MotherFrame
    = new G4Box("motherframe",
                circuit_sizeX*0.5*mm,
                circuit_sizeY*0.5*mm,
                mother_frame_thickness*0.5*mm);

  BaseSubtraction//subtraction module (center sensor hole)                   
    = new G4Box("sub",
                silicon_sizeXY*0.5*mm,
                silicon_sizeXY*0.5*mm,
                mother_frame_thickness*0.5*mm);

  //--- Circuit ---
  SoliCircuit
    = new G4SubtractionSolid("circuit",
                             SiliconCircuit,
                             BaseSubtraction,
                             0,
                             G4ThreeVector(0,0,0));
 
  LogicCircuit
    = new G4LogicalVolume(SoliCircuit,
			  FR4,
			  "circuit");
  PhysCircuitU
    = new G4PVPlacement(0,
                        G4ThreeVector(0,0,(Magnet_center+silicon_centerU+silicon_thickness*0.5)*mm),
                        LogicCircuit,
                        "circuit",
                        fLogicWorld,
                        false,
                        0);

  PhysCircuitD
    = new G4PVPlacement(0,
                        G4ThreeVector(0,0,(Magnet_center+silicon_centerD+silicon_thickness*0.5)*mm),
                        LogicCircuit,
                        "circuit",
                        fLogicWorld,
                        false,
                        1);

  //--- Mother Frame ---  
  SoliMotherFrame
    = new G4SubtractionSolid("MF",
                             MotherFrame,
                             BaseSubtraction,
                             0,
                             G4ThreeVector(0,0,0));

  LogicMotherFrame
    = new G4LogicalVolume(SoliMotherFrame,
			  Aluminum,
			  "MF");

  PhysMotherFrameU
    = new G4PVPlacement(0,
                        G4ThreeVector(0,0,
                                      (Magnet_center+silicon_centerU+circuit_thickness
                                       +(silicon_thickness+mother_frame_thickness)*0.5)*mm),
                        LogicMotherFrame,
                        "MF",
                        fLogicWorld,
                        false,
                        0);

  PhysMotherFrameD
    = new G4PVPlacement(0,
                        G4ThreeVector(0,0,
                                      (Magnet_center+silicon_centerD+circuit_thickness
                                      +(silicon_thickness+mother_frame_thickness)*0.5)*mm),
                        LogicMotherFrame,
			"MF",
			fLogicWorld,
			false,
			1);

  //----- Base Plate ------
  BasePlate
    = new G4Box("BasePlate",
                circuit_sizeX*0.5*mm,
                circuit_sizeY*0.5*mm,
                baseplate_thickness*0.5*mm);

  BasePlateFoot
    = new G4Box("BasePlate",
                baseplate_footwidth*0.5*mm,
                baseplate_footheight*0.5*mm,
                baseplate_thickness*0.5*mm);

  BasePlateWana
    = new G4SubtractionSolid("BasePlateWana",
                             BasePlate,
                             BaseSubtraction,
                             0,
                             G4ThreeVector(0,0,0));

  SoliBasePlate
    = new G4UnionSolid("BP",
                       BasePlateWana,
                       BasePlateFoot,
                       0,
                       G4ThreeVector(0,
                                     -(circuit_sizeY+baseplate_footheight)*0.5*mm,
                                     0));
  LogicBasePlate
    = new G4LogicalVolume(SoliBasePlate,
			  Aluminum,
			  "BP");

  PhysBasePlate
    = new G4PVPlacement(0,
                        G4ThreeVector(0,0,
                                      (Magnet_center+silicon_centerD
                                       +circuit_thickness+mother_frame_thickness
                                       +(silicon_thickness+baseplate_thickness)*0.5)*mm),
                        LogicBasePlate,
                        "BP",
                        fLogicWorld,
                        false,
                        0);

  //----- Step Limiter ------
  G4UserLimits* stepLimit;
  G4UserLimits* stepLimit2;
  stepLimit = new G4UserLimits(10*mm);
  stepLimit2 = new G4UserLimits(0.5*mm);
  fLogicWorld->SetUserLimits(stepLimit); // set the steplimiter of the world

  std::ofstream fout_volume;
  fout_volume.open("volume.dat",std::ios::out);

  // the steplimiter in the world's daughter
  G4VPhysicalVolume *myPVolume, *myPVolumeDaughter;
  G4LogicalVolume *myLVolume, *myLVolumeDaughter;
  G4int NofWorldDaughter = fLogicWorld->GetNoDaughters(); // get the number of daughter in the world

  fout_volume << "World's Daughters list:" << NofWorldDaughter << "\n" << std::endl;  
  
  for(G4int i=0;i<NofWorldDaughter;i++){
    myPVolume = fLogicWorld->GetDaughter(i);
    myLVolume = myPVolume->GetLogicalVolume();
    if(myLVolume->GetName()=="Vacuum"||myLVolume->GetName()=="Block") myLVolume->SetUserLimits(stepLimit); // set the steplimiter of the world's daughter
    else myLVolume->SetUserLimits(stepLimit2);
    fout_volume << myLVolume->GetName() << ":" << myLVolume->GetMaterial() << std::endl;
    G4int NofDaughter = myLVolume->GetNoDaughters();
    fout_volume << myLVolume->GetName() <<"'s Daughters list:" << NofDaughter << "\n" <<std::endl;
    if(NofDaughter==0) continue;
    else{
      for(G4int j=0; j<NofDaughter; j++){
        myPVolumeDaughter = myLVolume->GetDaughter(j);
        myLVolumeDaughter = myPVolumeDaughter->GetLogicalVolume();
        myLVolumeDaughter->SetUserLimits(stepLimit2);
        int NofDaughter2 = myLVolumeDaughter->GetNoDaughters();
	fout_volume << myLVolumeDaughter->GetName() << "'s Daughters:" << NofDaughter2 << "\n" <<std::endl;
	fout_volume << myLVolumeDaughter->GetName() << ":" << myLVolumeDaughter->GetMaterial() << std::endl;
      }
    }
  }
  fout_volume.close();
    
  //
  // Visualization attributes
  //
  // fLogicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());

  //
  //always return the physical World
  //
  G4VisAttributes* ChamberFoilAtt = new G4VisAttributes(TRUE,G4Colour(1.0,0.25,0.0));
  fLogicChamberFoil->SetVisAttributes(ChamberFoilAtt);
  G4VisAttributes* CavityAtt = new G4VisAttributes(TRUE,G4Colour(1,0.2,0.0));
  fLogicCavity->SetVisAttributes(CavityAtt);
  G4VisAttributes* WorldAtt = new G4VisAttributes(FALSE,G4Colour(0.9,0.9,0.9));
  fLogicWorld->SetVisAttributes(WorldAtt);
  G4VisAttributes* CavityFlangeAtt = new G4VisAttributes(TRUE,G4Colour(1,0.2,0.0));
  fLogicCavityFlange->SetVisAttributes(CavityFlangeAtt);
  G4VisAttributes* CavityFoilAtt = new G4VisAttributes(TRUE,G4Colour(1.0,0.25,0.0));
  fLogicCavityFoil->SetVisAttributes(CavityFoilAtt);
  G4VisAttributes* BPMAtt = new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.2));
  fLogicBPM->SetVisAttributes(BPMAtt);
  G4VisAttributes* KaptonAtt = new G4VisAttributes(TRUE,G4Colour(0.3,0.2,0.0));
  fLogicKapton->SetVisAttributes(KaptonAtt);

  G4VisAttributes* VaccumAtt = new G4VisAttributes(FALSE);
  fLogicVacuum->SetVisAttributes(VaccumAtt);
  G4VisAttributes* BlockAtt = new G4VisAttributes(FALSE);
  fLogicBlock->SetVisAttributes(BlockAtt);

  G4VisAttributes* SiliconAtt = new G4VisAttributes(TRUE,G4Colour(0.1,0.5,0.8));
  logi_si_xstrip->SetVisAttributes(SiliconAtt);
  logi_si_ystrip->SetVisAttributes(SiliconAtt);
  G4VisAttributes* CircuitAtt = new G4VisAttributes(TRUE,G4Colour(0.2,0.7,0.3));
  LogicCircuit->SetVisAttributes(CircuitAtt);
  G4VisAttributes* MFAtt = new G4VisAttributes(TRUE,G4Colour(0.6,0.6,0.6));
  LogicMotherFrame->SetVisAttributes(MFAtt);
  G4VisAttributes* BPAtt = new G4VisAttributes(TRUE,G4Colour(0.6,0.6,0.6));
  LogicBasePlate->SetVisAttributes(BPAtt);
    
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal F05Field* F05DetectorConstruction::fField = 0;

void F05DetectorConstruction::ConstructSDandField(){
  if (!fField) {
    fField = new F05Field();

//     G4RepleteEofM* equation = new G4RepleteEofM(fField);
     G4EqEMFieldWithSpin* equation = new G4EqEMFieldWithSpin(fField);
//     equation->SetBField();
//     equation->SetEField();
//     equation->SetSpin();

     G4FieldManager* fieldManager
       = G4TransportationManager::GetTransportationManager()->GetFieldManager();
     fieldManager->SetDetectorField(fField);

     G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation,12);
     G4double minStep = 0.01*mm;
     G4ChordFinder* chordFinder =
       new G4ChordFinder((G4MagneticField*)fField,minStep,stepper);

     // Set accuracy parameters
     G4double deltaChord = 3.0*mm;
     chordFinder->SetDeltaChord( deltaChord );

     G4double deltaOneStep = 0.01*mm;
     fieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

     G4double deltaIntersection = 0.1*mm;
     fieldManager->SetDeltaIntersection(deltaIntersection);
     
     G4TransportationManager* transportManager =
                           G4TransportationManager::GetTransportationManager();

     G4PropagatorInField* fieldPropagator =
                                      transportManager->GetPropagatorInField();

     G4double epsMin = 2.5e-7*mm;
     G4double epsMax = 0.05*mm;

     fieldPropagator->SetMinimumEpsilonStep(epsMin);
     fieldPropagator->SetMaximumEpsilonStep(epsMax);
     
     fieldManager->SetChordFinder(chordFinder);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
