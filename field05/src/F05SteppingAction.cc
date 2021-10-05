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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05SteppingAction::F05SteppingAction(void){;}
F05SteppingAction::~F05SteppingAction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4String processName = aStep->GetPostStepPoint()->
                               GetProcessDefinedStep()->GetProcessName();
  G4Track* aTrack= aStep->GetTrack();
  G4String particleName = aStep->GetTrack()->
                                    GetDefinition()->GetParticleName();
  //const G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  //const G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
 
  std::ofstream fout;
  fout.open("output.dat",std::ios::app);
  /*
  fout<<"time"<<"\t"<<"track_id"<<"\t"<<"step_num"<<"\t"<<"X"<<"\t"<<"Y"<<"\t"<<"Z"<<"\t"
      <<"Pxdir"<<"\t"<<"Pydir"<<"\t"<<"Pzdir"<<std::endl;
  */
  
  if(processName != "DecayWithSpin" && particleName == "mu+"){
    G4ThreeVector pos  = aTrack->GetPosition()/mm;
    G4ThreeVector pol  = aTrack->GetPolarization();
    G4int track_id = aTrack->GetTrackID();
    G4int step_num = aTrack->GetCurrentStepNumber();
    G4ThreeVector momDir  = aTrack->GetMomentumDirection();
    G4double gtime = aTrack->GetGlobalTime()/ns;

    fout<<gtime<<"\t"<<particleName<<"\t"
	<<track_id<<"\t"<<step_num<<"\t"<<pos.x()<<"\t"<<pos.y()<<"\t"<<pos.z()<<"\t"
	<<momDir.x()<<"\t"<<momDir.y()<<"\t"<<momDir.z()
	<<std::endl;
    /*
      if (momDir * polDir < (1.-1.E-7)) {
	G4double cos_theta = momDir * polDir;
	G4cout << " *** ERROR - WARNING *** " << G4endl;
	G4cout << "processName: " << processName << G4endl;
	G4cout << "particleName " << particleName << G4endl;
	G4cout << "Global Time: " << gTime/ns << "nsec" << G4endl;
	G4cout << "Angle between spin and momentum:" << cos_theta << G4endl;
	G4Exception("SteppingAction::UserSteppingAction","Error",
		    FatalException,
		    "Angle between spin and momentum too large");
      }
    */
  }
}

