//----------------------------------------------------------------------
//This program defines the configuration of position and size of MuSEUM materials.
//All units are written in mm.
//If you want to change detector size, change appropriate parameters.
//2017/04/08
//written by Shoichiro Nishimura
//----------------------------------------------------------------------

#ifndef _DetectorConfiguration_hh_
#define _DetectorConfiguration_hh_

const double Magnet_center=1200.;
//--- Gas Chamber ---
const double chamber_length=390.0;//z direction length
const double chamber_diameter=400.0;//inside diameter
const double chamber_thickness=30.;//outside - inside diameter

//--- Gas Chamber Flange ---
const double chamber_flange_diameter=430.;
const double chamber_flange_thickness_u=30.;//z direction length
const double chamber_flange_thickness_d=30.;//z direction length
const double chamber_window_diameter=100.;

const double overhang_diameter=390.;//alternative absorber
const double overhang_thickness=25.;

//--- Gas Chamber Foil ---
const double chamber_foil_thickness=0.1;

//--- RF Cavity ---
const double cavity_length=244.;
const double cavity_diameter=187.;
const double cavity_thickness=15.;
const double cavity_center=0;//the origin is the center of gas chamber

//--- RF Cavity Flange ---
const double cavity_flange_length=44.;
const double cavity_flange_thickness=25.;//z direction length

//--- RF Cavity Foil ---
const double cavity_foil_thickness=0.025;//z direction length

//--- Beam Window Kapton Foil ---
const double kapton_diameter=300.;
const double kapton_thickness=0.075;
const double kapton_center=-500.0;

//--- Beam Profile Monitor ---
const double bpm_sizeXY=100.;//BPM size
const double bpm_thickness=0.15;
//const double bpm_positionU=-335;//up stream bpm center position
//const double bpm_positionD=-330;//down stream bpm center position
//const double bpm_positionU=-241;//up stream bpm center position
//const double bpm_positionD=-240;//down stream bpm center position
const double bpm_positionU=kapton_center+5;
const double bpm_positionD=kapton_center+10;

//--- Al Absorber ---
const double al_sizeXY=240.;
const double al_thickness=0.;//absorber thickness
const double al_position=
  chamber_length*0.5+chamber_flange_thickness_d+al_thickness*0.5;

//--- Positron Counter ---
const double counter_sizeXY=240.0;
const double counter_thickness=3;
const double counter_centerU=300.;
const double counter_centerD=340.;
const double counter_couver_thickness=2;

//--- Silicon Detector ---
const int silicon_strip_num=512;
const double silicon_sizeXY=98.77;
const double silicon_thickness=0.32;
const double silicon_strip_length=48.595;
const double silicon_strip_pitch=0.19;
const double silicon_edge=0.745;
const double silicon_strip_gap=0.5;

const double silicon_centerU=245;
const double silicon_centerD=265;
/*
const double silicon_centerU=380;
const double silicon_centerD=400;
*/
const double circuit_sizeX=400.;//circuit board whose color is green
const double circuit_sizeY=160.;
const double circuit_thickness=1.6;

const double mother_frame_thickness=3;

const double baseplate_thickness=3;
const double baseplate_footheight=115;
const double baseplate_footwidth=260;

//--- Magnetic Shield ---
const double shield_thickness=1.5;
const double b_shield_posz=375+shield_thickness*0.5;
const double b_shield_posy=80;
const double b_shield_width=500;
const double b_shield_height=550;
const double r_shieldhole=40;

#endif
