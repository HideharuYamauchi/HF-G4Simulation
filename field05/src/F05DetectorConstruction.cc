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
  Silicon=nistMan->FindOrBuildMaterial("G4_Si");
  Kapton=nistMan->FindOrBuildMaterial("G4_KAPTON");
  Aluminum=nistMan->FindOrBuildMaterial("G4_Al");
  Copper=nistMan->FindOrBuildMaterial("G4_Cu");
  lead=nistMan->FindOrBuildMaterial("G4_Pb");
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
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* F05DetectorConstruction::Construct()
{
  //----------------------------------
  // World
  //----------------------------------
  fWorldSizeXY = 1.2*m;
  fWorldSizeZ = 3.6*m;

  fSolidWorld = new G4Box("World",                               //its name
                   fWorldSizeXY/2,fWorldSizeXY/2,fWorldSizeZ/2); //its size
 
  fLogicWorld = new G4LogicalVolume(fSolidWorld,        //its solid
                                    fVacuum,            //its material, fVacuum,
                                    "World");           //its name
 
  fPhysiWorld = new G4PVPlacement(0,                    //no rotation
                                  G4ThreeVector(),      //at (0,0,0)
                                  fLogicWorld,          //its logical volume
                                  "World",              //its name
                                  0,                    //its mother  volume
                                  false,                //no boolean operation
                                  0);                   //copy number
  
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
                       G4ThreeVector(0,0,(Magnet_center+kapton_center)*mm),
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
  
  
  G4UserLimits* stepLimit;
  stepLimit = new G4UserLimits(5*mm);
  fLogicWorld->SetUserLimits(stepLimit);
  
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
  G4VisAttributes* KaptonAtt=new G4VisAttributes(TRUE,G4Colour(0.3,0.2,0.0));
  fLogicKapton->SetVisAttributes(KaptonAtt);
  
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal F05Field* F05DetectorConstruction::fField = 0;

void F05DetectorConstruction::ConstructSDandField()

{
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

     G4double minStep           = 0.01*mm;

     G4ChordFinder* chordFinder =
                    new G4ChordFinder((G4MagneticField*)fField,minStep,stepper);

     // Set accuracy parameters
     G4double deltaChord        = 3.0*mm;
     chordFinder->SetDeltaChord( deltaChord );

     G4double deltaOneStep      = 0.01*mm;
     fieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

     G4double deltaIntersection = 0.1*mm;
     fieldManager->SetDeltaIntersection(deltaIntersection);

     G4TransportationManager* transportManager =
                           G4TransportationManager::GetTransportationManager();

     G4PropagatorInField* fieldPropagator =
                                      transportManager->GetPropagatorInField();

     G4double epsMin            = 2.5e-7*mm;
     G4double epsMax            = 0.05*mm;

     fieldPropagator->SetMinimumEpsilonStep(epsMin);
     fieldPropagator->SetMaximumEpsilonStep(epsMax);

     fieldManager->SetChordFinder(chordFinder);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
