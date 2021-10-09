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
/// \file field/field05/src/F05SteppingAction.cc
/// \brief Implementation of the F05SteppingAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05SteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ProcessTable.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05SteppingAction::F05SteppingAction(void){;}
F05SteppingAction::~F05SteppingAction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05SteppingAction::UserSteppingAction(const G4Step* aStep){
  G4String processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  G4Track* aTrack = aStep->GetTrack();
  G4String particleName = aTrack->GetDefinition()->GetParticleName();  
  const G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  const G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  G4String volume = aTrack->GetTouchableHandle()->GetVolume()->GetName(); //volume = postvolume
  G4String prevolume = preStepPoint->GetTouchableHandle()->GetVolume()->GetName();
  //G4String postvolume = postStepPoint->GetTouchableHandle()->GetVolume()->GetName(); //error 
  
  std::ofstream fout;
  fout.open("output.dat",std::ios::app);
  /*
  if((processName=="DecayWithSpin"||volume=="silicon"||volume=="Detector"||volume=="Chamber"
      ||volume=="ChamberFlangeU"||volume=="ChamberFlangeD"||volume=="ChamberFoil"
      ||volume=="Cavity"||volume=="CavityFlange"||volume=="CavityFoil"||volume=="TargetGas")
     &&(particleName == "e+"||particleName == "mu+")){
  */
  if(processName=="DecayWithSpin"){
    //G4ThreeVector pos  = aTrack->GetPosition()/mm; // pos = postpos
    G4ThreeVector prepos  = preStepPoint->GetPosition()/mm;
    G4ThreeVector postpos  = postStepPoint->GetPosition()/mm;
    //G4ThreeVector pol  = aTrack->GetPolarization(); // pol = postpol
    G4ThreeVector prepol = preStepPoint->GetPolarization();
    G4ThreeVector postpol = postStepPoint->GetPolarization();
    
    G4int track_id = aTrack->GetTrackID();
    G4int step_num = aTrack->GetCurrentStepNumber();
    
    //G4ThreeVector momDir = aTrack->GetMomentumDirection(); // momDir = postmomDir
    G4ThreeVector premomDir = preStepPoint->GetMomentumDirection();
    G4ThreeVector postmomDir = postStepPoint->GetMomentumDirection();
    //G4double gtime = aTrack->GetGlobalTime()/ns; //gtime = postgtime
    G4double pregtime = preStepPoint->GetGlobalTime()/ns;
    G4double postgtime = postStepPoint->GetGlobalTime()/ns;
    //G4double kE = aTrack->GetKineticEnergy()/keV; // kE = postkE
    G4double prekE = preStepPoint->GetKineticEnergy()/keV;
    G4double postkE = postStepPoint->GetKineticEnergy()/keV;
   
    fout<<pregtime<<"\t"<<postgtime<<"\t" 
	<<prevolume<<"\t"<<volume<<"\t"<<processName<<"\t"
	<<prekE<<"\t"<<postkE<<"\t"
	<<particleName<<"\t"<<track_id<<"\t"<<step_num<<"\t"
	<<prepos.x()<<"\t"<<prepos.y()<<"\t"<<prepos.z()<<"\t"
	<<postpos.x()<<"\t"<<postpos.y()<<"\t"<<postpos.z()<<"\t"
	<<premomDir.x()<<"\t"<<premomDir.y()<<"\t"<<premomDir.z()<<"\t"
	<<postmomDir.x()<<"\t"<<postmomDir.y()<<"\t"<<postmomDir.z()<<"\t"
	<<std::endl;
  }
  fout.close(); 
}
