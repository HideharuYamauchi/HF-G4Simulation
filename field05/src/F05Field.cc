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
/// \file field/field05/src/F05Field.cc
/// \brief Implementation of the F05Field class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05Field.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05Field::F05Field() : G4ElectroMagneticField(){
  file = TFile::Open("/home/hideharu/g4sample/field05/src/at_hline-20130429-1.root","read");
  if(file->IsOpen()){
    tree = (TTree*)file->Get("FNT/FNTuple");
    entries=tree->GetEntries();
    tree->SetBranchAddress("x",&x);
    tree->SetBranchAddress("z",&z);
    tree->SetBranchAddress("Bx",&Bx_data);
    tree->SetBranchAddress("By",&By_data);
    tree->SetBranchAddress("Bz",&Bz_data);
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("x",1);
    tree->SetBranchStatus("z",1);
    tree->SetBranchStatus("Bx",1);
    tree->SetBranchStatus("By",1);
    tree->SetBranchStatus("Bz",1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05Field::~F05Field(){
  file->Close();
  delete tree;
  delete file;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05Field::GetFieldValue( const G4double Point[4], G4double* Bfield ) const{
  // map_center=5355, z position
  // Point[0],Point[1],Point[2] are x-, y-, z-cordinates, Point[3] is time
  //std::ofstream fout_magfield;
  //fout_magfield.open("magfield.dat",std::ios::app); // additional
  G4double Bx = 0.*tesla;
  G4double By = 0.*tesla;
  G4double Bz = 0.*tesla;
  G4double Ex = 0.*volt/m;
  G4double Ey = 0.*volt/m;
  G4double Ez = 0.*volt/m;
  /*
  double r = std::sqrt(Point[0]*Point[0]+Point[1]*Point[1]);
  double posz = Point[2]-Magnet_center;
  int cut_r = r*0.1;
  int cut_z = (posz+map_center+100)*0.1;
  if(cut_z>=660) cut_z-=1;
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      int entry = (cut_r+i)*661+cut_z+j;
      //std::cout << "entry:" << entry << std::endl;
      if(entry>=0&&entry<entries){
	tree->GetEntry(entry);
	double mapz = z-map_center;
	Bx+=Bx_data*(1.-std::abs(x-r)*0.1)*(1.-std::abs(mapz-posz)*0.1);
	By+=By_data*(1.-std::abs(x-r)*0.1)*(1.-std::abs(mapz-posz)*0.1);
	Bz+=Bz_data*(1.-std::abs(x-r)*0.1)*(1.-std::abs(mapz-posz)*0.1);
      }
    }
  }
  double angle=0;
  double Bxy=std::sqrt(Bx*Bx+By*By);
  if(Bxy>0){
    if(r>0) angle=acos(Point[0]/r);
    if(Point[1]<0) angle*=-1;
    if(By>0) angle+=acos(Bx/Bxy);
    else angle-=acos(Bx/Bxy);
  }
  Bx=Bxy*std::cos(angle)*tesla;
  By=Bxy*std::sin(angle)*tesla;
  Bz*=tesla;
  */
  Bfield[0] = Bx;
  Bfield[1] = By;
  Bfield[2] = Bz;
  Bfield[3] = Ex;
  Bfield[4] = Ey;
  Bfield[5] = Ez;

  /*
  fout_magfield <<Point[0]<<"\t"<<Point[1]<<"\t"<<Point[2]<<"\t"<<Point[3]<<"\t"
		<<Bx<<"\t"<<By<<"\t"<<Bz
		<<std::endl;
  fout_magfield.close();
  */
  return;
}
