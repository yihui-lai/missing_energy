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
// * institutes, nor the agencies providing financial support for this *
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
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#include "DetectorConstruction.hh"
#include "DR_PMTSD.hh"
#include "SurfaceProperty.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4UserLimits.hh"

#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include "CreateTree.hh"
#include <algorithm>
#include <string>
#include <sstream>
#include "G4MultiUnion.hh"

using namespace CLHEP;

//validation of DREAM, with one SiPMs at each ends 

DetectorConstruction::DetectorConstruction(const string &configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------

  ConfigFile config(configFileName);

  config.readInto(checkOverlaps, "checkOverlaps");

  config.readInto(world_material, "world_material");
  config.readInto(bar_length, "bar_length");

  config.readInto(core_radius_x, "core_radius_x");
  config.readInto(core_radius_y, "core_radius_y");
  config.readInto(core_material, "core_material");
  config.readInto(core_rIndex, "core_rIndex");
  config.readInto(core_absLength, "core_absLength");

  config.readInto(gap_l, "gap_l");
  config.readInto(gap_size_x, "gap_size_x");
  config.readInto(gap_size_y, "gap_size_y");
  config.readInto(gap_material, "gap_material");

  config.readInto(det_l, "det_l");
  config.readInto(det_size_x, "det_size_x");
  config.readInto(det_size_y, "det_size_y");
  config.readInto(det_material, "det_material");
  config.readInto(ecal_det_size, "ecal_det_size");

  config.readInto(depth, "depth");
  config.readInto(cryst_dist, "cryst_dist");
  config.readInto(trackerX0, "trackerX0");
  config.readInto(services_thick, "services_thick");

  config.readInto(ecal_incline, "ecal_incline");
  config.readInto(ecal_xy_gap, "ecal_xy_gap");
  config.readInto(fiber_type, "fiber_type");
  config.readInto(ecal_material, "ecal_material");
  config.readInto(scinti_material, "scinti_material");
  config.readInto(Cherenc_material, "Cherenc_material");
  config.readInto(Cherenp_material, "Cherenp_material");
  config.readInto(wrap_material, "wrap_material");
  config.readInto(wrap_ref, "wrap_ref");
  config.readInto(rear_filter, "rear_filter");
  config.readInto(front_filter, "front_filter");
  config.readInto(ecal_front_length, "ecal_front_length");
  config.readInto(ecal_rear_length, "ecal_rear_length");
  config.readInto(ecal_front_face, "ecal_front_face");
  config.readInto(ecal_rear_face, "ecal_rear_face");
  config.readInto(hole_diameter, "hole_diameter");
  config.readInto(fiber_diameter, "fiber_diameter");
  config.readInto(wrapping_thick, "wrapping_thick");

  config.readInto(ecal_timing_distance, "ecal_timing_distance");

  config.readInto(hcal_width, "hcal_width");
  config.readInto(hcalTile_width, "hcalTile_width");
  config.readInto(hcalAbs_1_thick, "hcalAbs_1_thick");
  config.readInto(hcalAbs_2_thick, "hcalAbs_2_thick");
  config.readInto(solenoid_thick, "solenoid_thick");
  config.readInto(hcalTile_thick, "hcalTile_thick");

  B_field_intensity = config.read<double>("B_field_intensity") * tesla;

  expHall_x = 300. * cm;
  expHall_y = 300. * cm;
  expHall_z = 300. * cm;

  B_field_IsInitialized = false;

  initializeMaterials();
  initializeSurface();

  CreateTree::Instance()->inputTrackerX0 = trackerX0;
  CreateTree::Instance()->inputServiceAlmm = services_thick;
  CreateTree::Instance()->inputTimingThick = core_radius_x * 2;
  CreateTree::Instance()->inputE1Thick = ecal_front_length;
  CreateTree::Instance()->inputE2Thick = ecal_rear_length;
  CreateTree::Instance()->inputE1Width = ecal_front_face;
  CreateTree::Instance()->inputTimingECAL_dist = ecal_timing_distance;
}

//---- ---- ---- ---- ---- ---- ---- ---- ----  ---- ---- ---- ---- ---- ----

DetectorConstruction::~DetectorConstruction()
{
  delete stepLimit;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl;

  //------------------------------------
  //------------- Geometry -------------
  //------------------------------------

  // The experimental Hall
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid *worldS = new G4Box("worldS", 5 * expHall_x, 5 * expHall_y, 5 * expHall_z);
  G4LogicalVolume *worldLV = new G4LogicalVolume(worldS, WoMaterial, "worldLV", 0, 0, 0);
  G4VPhysicalVolume *worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "worldPV", 0, false, 0, checkOverlaps);

  // the Tracker
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid *trackerS;
  G4LogicalVolume *trackerLV = NULL;
  float gapTRK_Timing = 20 * mm;
  if (trackerX0 > 0)
  {

    int NLAYERS = 6;
    float X0 = 93.7 * mm;

    float layerSpacing = 1800 / NLAYERS * mm;
    float layerThick = trackerX0 * X0 / NLAYERS;

    trackerS = new G4Box("trackerS", 200 * mm, 200 * mm, 0.5 * layerThick);
    trackerLV = new G4LogicalVolume(trackerS, MyMaterials::Silicon(), "trackerLV");

    for (int iLayer = 0; iLayer < NLAYERS; iLayer++)
    {
//      new G4PVPlacement(0, G4ThreeVector(0., 0., -layerSpacing * (NLAYERS - iLayer) - gapTRK_Timing), trackerLV, Form("trackerPV_Layer%d", iLayer), worldLV, false, 0, checkOverlaps);
    }
  }

  // the cooling-services Dead material layer: 5 cm
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid *servicesS;
  G4LogicalVolume *servicesLV = NULL;
  float pos_Z_DT1 = -20 * mm;
  float pos_Z_DT2 = 30 * mm;
  float pos_Z_DE1 = 70 * mm;
  float pos_Z_DE2 = ecal_timing_distance + ecal_front_length + ecal_rear_length + 20 * mm;

  if (services_thick > 0)
  {
    servicesS = new G4Box("servicesS", 100 * mm, 100 * mm, 0.5 * services_thick);
    servicesLV = new G4LogicalVolume(servicesS, MyMaterials::Aluminium(), "servicesLV");
/*    new G4PVPlacement(0, G4ThreeVector(0., 0., pos_Z_DT1), servicesLV, "servicesPV_T1", worldLV, false, 0, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0., 0., pos_Z_DT2), servicesLV, "servicesPV_T2", worldLV, false, 0, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0., 0., pos_Z_DE1), servicesLV, "servicesPV_E1", worldLV, false, 0, checkOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0., 0., pos_Z_DE2), servicesLV, "servicesPV_E2", worldLV, false, 0, checkOverlaps);
*/  }

  //*********************
  //  TIMING LAYER
  //*********************

  float pos_Z_T1 = 0 * mm;
  float pos_Z_T2 = 6 * mm;

  G4RotationMatrix *piRotTiming_front = new G4RotationMatrix;
  G4RotationMatrix *piRotTiming_rear = new G4RotationMatrix;

  //  piRotTiming->rotateY(-pointingAngle*rad*(iX + 0.5));
  piRotTiming_front->rotateX(90 * deg);
  piRotTiming_rear->rotateY(90 * deg);
  //  piRotTiming->rotateY(180*rad);

  int NBARS = 18;
  float bar_spacing = core_radius_x * 2;

  //TIMING SOLID
  // the crystals
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid *coreS;
  coreS = new G4Box("coreS", core_radius_x, core_radius_y, 0.5 * bar_length);
  G4LogicalVolume *coreLV = new G4LogicalVolume(coreS, CoMaterial, "coreLV");

  // bar end gaps for photon counting
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  //end gaps
  G4VSolid *gapLayerS;
  G4VSolid *gapS;
  if (gap_size_x == -1)
  {
    //	gapLayerS = new G4Box  ("gapLayerS", core_radius_x, core_radius_y, 0.5*depth) ;
    gapS = new G4Box("gapS", core_radius_x, core_radius_y, 0.5 * (gap_l - depth));
  }
  else
  {
    //  	gapLayerS = new G4Box  ("gapLayerS", gap_size_x*0.5, gap_size_y*0.5, 0.5*depth) ;
    gapS = new G4Box("gapS", gap_size_x * 0.5, gap_size_y * 0.5, 0.5 * (gap_l - depth));
  }

  // Si detector for photon counting at T1/T2 bar end
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid *detLayerS = NULL;
  G4VSolid *detS = NULL;

  if (det_size_x == -1) //set the SiPM active area equal to the crystal surface
  {
    detLayerS = new G4Box("detLayerS", core_radius_x * mm, core_radius_y * mm, 0.5 * depth);
    detS = new G4Box("detS", core_radius_x * mm, core_radius_y * mm, 0.5 * (det_l - depth));
  }
  else if (det_size_x > 0)
  {
    detLayerS = new G4Box("detLayerS", det_size_x * 0.5, det_size_y * 0.5, 0.5 * depth);
    detS = new G4Box("detS", det_size_x * 0.5, det_size_y * 0.5, 0.5 * (det_l - depth));
  }

  // TIMING LOGIC
  G4LogicalVolume *gapLV = new G4LogicalVolume(gapS, GaMaterial, "gapLV");
  //  G4LogicalVolume* gapLayerLV = new G4LogicalVolume (gapLayerS, GaMaterial, "gapLayerLV") ;

  G4LogicalVolume *detLayerLV = NULL;
  G4LogicalVolume *detLV = NULL;
  detLayerLV = new G4LogicalVolume(detLayerS, DeMaterial, "detLayerLV");
  detLV = new G4LogicalVolume(detS, DeMaterial, "detLV");

  // TIMING PHYSICAL PLACEMENT

  G4PVPlacement *corePV;

  for (int iBar = 0; iBar < NBARS; iBar++)
  {
    float x = iBar * bar_spacing - bar_spacing * NBARS / 2; //+ bar_spacing/2;
    float y = -bar_spacing / 2.;
/*
    //first crystal layer
    corePV = new G4PVPlacement(piRotTiming_front, G4ThreeVector(x, y, pos_Z_T1), coreLV, Form("corePV_front_%d", iBar), worldLV, false, 0, checkOverlaps);
    //front, left and right gaps + sipms
    new G4PVPlacement(piRotTiming_front, G4ThreeVector(x, y - bar_length * 0.5 - gap_l * 0.5, pos_Z_T1), gapLV, Form("timingGap_front_1_PV_%d", iBar), worldLV, false, 0, checkOverlaps);
    new G4PVPlacement(piRotTiming_front, G4ThreeVector(x, y + bar_length * 0.5 + gap_l * 0.5, pos_Z_T1), gapLV, Form("timingGap_front_2_PV_%d", iBar), worldLV, false, 0, checkOverlaps);
    new G4PVPlacement(piRotTiming_front, G4ThreeVector(x, y - bar_length * 0.5 - gap_l - det_l * 0.5, pos_Z_T1), detLV, Form("timingDet_front_1_PV_%d", iBar), worldLV, false, 0, checkOverlaps);
    new G4PVPlacement(piRotTiming_front, G4ThreeVector(x, y + bar_length * 0.5 + gap_l + det_l * 0.5, pos_Z_T1), detLV, Form("timingDet_front_2_PV_%d", iBar), worldLV, false, 0, checkOverlaps);

    //second crystal layer
    corePV = new G4PVPlacement(piRotTiming_rear, G4ThreeVector(y, x, pos_Z_T2), coreLV, Form("corePV_rear_%d", iBar), worldLV, false, 0, checkOverlaps);
    //front, left and right gaps + sipms
    new G4PVPlacement(piRotTiming_rear, G4ThreeVector(y - bar_length * 0.5 - gap_l * 0.5, x, pos_Z_T2), gapLV, Form("timingGap_rear_1_PV_%d", iBar), worldLV, false, 0, checkOverlaps);
    new G4PVPlacement(piRotTiming_rear, G4ThreeVector(y + bar_length * 0.5 + gap_l * 0.5, x, pos_Z_T2), gapLV, Form("timingGap_rear_2_PV_%d", iBar), worldLV, false, 0, checkOverlaps);
    new G4PVPlacement(piRotTiming_rear, G4ThreeVector(y - bar_length * 0.5 - gap_l - det_l * 0.5, x, pos_Z_T2), detLV, Form("timingDet_rear_1_PV_%d", iBar), worldLV, false, 0, checkOverlaps);
    new G4PVPlacement(piRotTiming_rear, G4ThreeVector(y + bar_length * 0.5 + gap_l + det_l * 0.5, x, pos_Z_T2), detLV, Form("timingDet_rear_2_PV_%d", iBar), worldLV, false, 0, checkOverlaps);
*/
  }

  //SURFACE STATE IMPLEMENTATION
  /*
   //front face surface state (opposite to photodetector)
   G4Box* Front_skin 			= new G4Box("Front_skin", core_radius_x, core_radius_y, 0.5*depth);
   G4LogicalVolume* Front_skin_log 	= new G4LogicalVolume(Front_skin, MyMaterials::Air(), "Front_skin_log", 0, 0, 0);
//   G4VPhysicalVolume* Front_skin_phys 	= new G4PVPlacement(0, G4ThreeVector(0., 0., -(0.5*bar_length+0.5*depth)) , Front_skin_log, "Front_skin_phys", worldLV, false, 0);

   //side surface state
   G4VSolid* dummySk = new G4Box ("dummySk", core_radius_x+depth, core_radius_y+depth, 0.5*bar_length) ;
   G4VSolid* subSk   = new G4Box ("subSk", core_radius_x, core_radius_y, 0.51*bar_length);
   
   G4SubtractionSolid* Side_skin 	= new G4SubtractionSolid("Side_skin", dummySk, subSk);				
   G4LogicalVolume* Side_skin_log 	= new G4LogicalVolume(Side_skin, MyMaterials::Air(), "Side_skin_log", 0, 0, 0);
  // G4VPhysicalVolume* Side_skin_phys 	= new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Side_skin_log, "Side_skin_phys", worldLV, false, 0);

   G4LogicalBorderSurface *CrystalFrontSkin   	= NULL;
   G4LogicalBorderSurface *CrystalSideSkin   	= NULL;
   G4OpticalSurface *OpCrystalSurface		= NULL;
	
   ///-------CRYSTAL SURFACE-------
   OpCrystalSurface = new G4OpticalSurface("crystal");
   OpCrystalSurface->SetType(dielectric_dielectric);
   OpCrystalSurface->SetModel(unified);

*/

  //********************************************
  //  ELECTROMAGNETIC CALORIMETER
  //********************************************

  //  float pos_Z_E1 =  0*mm;
  //  float pos_Z_E2 =  pos_Z_T26*mm;
  //double pointingAngle = 0.001976598*M_PI;	//~0.36° -- arctan (29.31/4720)
  double pointingAngle = ecal_incline; //~0.36° -- arctan (29.31/4720)
                                       //  double alveola_thickness = 0.2 * mm;

  //  double pointingAngle = ecal_incline; //~0.36° -- arctan (29.31/4720)
  double alveola_thickness = ecal_xy_gap * mm;

  //  float total_ecal_length = ecal_front_length + ecal_rear_length;

  // ECAL solid
  G4Trd *ecalCrystalS_f = new G4Trd("ecalCrystalS_f", 0.5 * ecal_front_face, 0.5 * ecal_rear_face, 0.5 * ecal_front_face, 0.5 * ecal_rear_face, 0.5 * ecal_front_length);
  G4Trd *ecalCrystalS_r = new G4Trd("ecalCrystalS_r", 0.5 * ecal_front_face, 0.5 * ecal_rear_face, 0.5 * ecal_front_face, 0.5 * ecal_rear_face, 0.5 * ecal_rear_length);
  double window_l=0.6*mm;
  double Case_l=0.1*mm;
  //construct dream detector
  //
  double dchange_rotation = 0.01*mm;//sqrt((ecal_front_length * 0.5 + det_l + gap_l+window_l) * (ecal_front_length * 0.5 + det_l + gap_l+window_l) + ecal_front_face * ecal_front_face / 4.) + sqrt((ecal_rear_length * 0.5 + det_l + gap_l+window_l) * (ecal_rear_length * 0.5 + det_l + gap_l+window_l) + ecal_rear_face * ecal_rear_face / 4.) - ecal_front_length * 0.5 - ecal_rear_length * 0.5;

//  double dchange_rotation =  sqrt(ecal_front_face* ecal_front_face + (ecal_front_length+2*det_l+2*gap_l+2*window_l) * (ecal_front_length+2*det_l+2*gap_l+2*window_l)) * abs(cos(pointingAngle * deg - atan(ecal_front_face / (ecal_front_length+2*det_l+2*gap_l+2*window_l)))) / 2 - (ecal_front_length+2*det_l+2*gap_l+2*window_l) * 0.5 + sqrt(ecal_front_face * ecal_front_face + (ecal_rear_length+2*det_l+2*gap_l+2*window_l) * (ecal_rear_length+2*det_l+2*gap_l+2*window_l)) * abs(cos(pointingAngle * deg - atan(ecal_front_face / (ecal_rear_length+2*det_l+2*gap_l+2*window_l)))) / 2 - (ecal_rear_length+2*det_l+2*gap_l+2*window_l) * 0.5;


  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  //--------------------- build a dream detector here ------------------
  G4Box *ecalCrystalS_f_hole_s;
  G4Box *ecalCrystalS_r_hole_s;
  G4Box *ecalCrystalS_f_fiber_s;
  G4Box *ecalCrystalS_r_fiber_s;
  G4Box *ecalCrystalS_f_outwrap_s;
  G4Box *ecalCrystalS_r_outwrap_s;
  G4Box *ecalCrystalS_f_innerwrap_s;
  G4Box *ecalCrystalS_r_innerwrap_s;

  G4SubtractionSolid *ecalCrystalS_f_hollowwrap;
  G4SubtractionSolid *ecalCrystalS_r_hollowwrap;
  G4SubtractionSolid *ecalCrystalS_f_absorb;
  G4SubtractionSolid *ecalCrystalS_r_absorb;

  G4SubtractionSolid *ecalCrystalS_f_fiberGap;
  G4SubtractionSolid *ecalCrystalS_r_fiberGap;

  G4LogicalVolume *ecalCrystalL_f_absorb;
  G4LogicalVolume *ecalCrystalL_r_absorb;
  G4LogicalVolume *ecalCrystalL_f_fiber_wrapping;
  G4LogicalVolume *ecalCrystalL_r_fiber_wrapping;
  G4LogicalVolume *ecalCrystalL_f_fiberGap;
  G4LogicalVolume *ecalCrystalL_r_fiberGap;

  G4LogicalVolume *ecalCrystalL_f_fiber;
  G4LogicalVolume *ecalCrystalL_r_fiber;

  G4Box *ecalGapS;
  G4Box *ecalDetS;
  G4Box *ecalDetWindowS;

  G4Box *ecalDetCaseOuterS;
  G4Box *ecalDetCaseInnerS;
  G4SubtractionSolid *ecalDetCase_fS;
  G4SubtractionSolid *ecalDetCase_rS;

  G4LogicalVolume *ecalGapL;
  G4LogicalVolume *ecalDetL;
  G4LogicalVolume *ecalDetWindowL;
  G4LogicalVolume *ecalDetCase_fL;
  G4LogicalVolume *ecalDetCase_rL;

  G4Box *ecalStop = new G4Box("ecalStop", ecal_rear_face*10, ecal_rear_face*10, 0.3*dchange_rotation);
  G4Box *ecalStop_gap = new G4Box("ecalStop_gap", ecal_rear_face*10, ecal_rear_face*10, 0.1*dchange_rotation);
  G4LogicalVolume *ecalDetStopL = new G4LogicalVolume(ecalStop, WrapMaterial, "ecalStopL", 0, 0, 0);
  G4LogicalVolume *ecalDetStopL_gap = new G4LogicalVolume(ecalStop_gap, GaMaterial, "ecalStopL_gap", 0, 0, 0);

  if (fiber_type == 0)
  {
    //circular fiber
  }
  else
  {
    //sqaure fiber

//---------------------- wrapping  --------------------------
    ecalCrystalS_f_outwrap_s = new G4Box("ecalCrystalS_f_outwrap_s", ecal_front_face * 0.5 + wrapping_thick + depth, ecal_front_face * 0.5 + wrapping_thick + depth, ecal_front_length * 0.5);
    ecalCrystalS_f_innerwrap_s = new G4Box("ecalCrystalS_f_innerwrap_s", ecal_front_face * 0.5 + depth, ecal_front_face * 0.5 + depth, ecal_front_length);
    ecalCrystalS_r_outwrap_s = new G4Box("ecalCrystalS_r_outwrap_s", ecal_front_face * 0.5 + wrapping_thick + depth, ecal_front_face * 0.5 + wrapping_thick + depth, ecal_rear_length * 0.5);
    ecalCrystalS_r_innerwrap_s = new G4Box("ecalCrystalS_r_innerwrap_s", ecal_front_face * 0.5 + depth, ecal_front_face * 0.5 + depth, ecal_rear_length);
    ecalCrystalS_f_hollowwrap = new G4SubtractionSolid("ecalCrystalS_f_hollowwrap", ecalCrystalS_f_outwrap_s, ecalCrystalS_f_innerwrap_s);
    ecalCrystalS_r_hollowwrap = new G4SubtractionSolid("ecalCrystalS_r_hollowwrap", ecalCrystalS_r_outwrap_s, ecalCrystalS_r_innerwrap_s);
    ecalCrystalL_f_fiber_wrapping = new G4LogicalVolume(ecalCrystalS_f_hollowwrap, WrapMaterial, "ecalCrystalL_f_fiber_wrapping", 0, 0, 0);
    ecalCrystalL_r_fiber_wrapping = new G4LogicalVolume(ecalCrystalS_r_hollowwrap, WrapMaterial, "ecalCrystalL_r_fiber_wrapping", 0, 0, 0);
//---------------------- fiber  --------------------------
    ecalCrystalL_f_fiber = new G4LogicalVolume(ecalCrystalS_f, CherencMaterial, "ecalCrystalL_f_fiber", 0, 0, 0);
    ecalCrystalL_r_fiber = new G4LogicalVolume(ecalCrystalS_r, CherencMaterial, "ecalCrystalL_r_fiber", 0, 0, 0);
//---------------------- fiberGap  --------------------------
    ecalCrystalS_f_fiberGap  = new G4SubtractionSolid("ecalCrystalS_f_fiberGap", ecalCrystalS_f_outwrap_s, ecalCrystalS_f_hollowwrap);
    ecalCrystalS_f_fiberGap  = new G4SubtractionSolid("ecalCrystalS_f_fiberGap", ecalCrystalS_f_fiberGap, ecalCrystalS_f);
    ecalCrystalS_r_fiberGap  = new G4SubtractionSolid("ecalCrystalS_r_fiberGap", ecalCrystalS_r_outwrap_s, ecalCrystalS_r_hollowwrap);
    ecalCrystalS_r_fiberGap  = new G4SubtractionSolid("ecalCrystalS_r_fiberGap", ecalCrystalS_r_fiberGap, ecalCrystalS_r);
    ecalCrystalL_f_fiberGap = new G4LogicalVolume(ecalCrystalS_f_fiberGap, GaMaterial, "ecalCrystalL_f_fiberGap", 0, 0, 0);
    ecalCrystalL_r_fiberGap = new G4LogicalVolume(ecalCrystalS_r_fiberGap, GaMaterial, "ecalCrystalL_r_fiberGap", 0, 0, 0);


  ecalGapS = new G4Box("ecalGapS", ecal_front_face*0.5, ecal_front_face*0.5 , 0.5*gap_l);
//  ecalGapS = new G4Box("ecalGapS", ecal_det_size, ecal_det_size , 0.5 * (gap_l));
  ecalDetS = new G4Box("ecalDetS", ecal_det_size*0.5, ecal_det_size*0.5 , 0.5 * (det_l));
  ecalDetWindowS = new G4Box("ecalDetWindowS", ecal_det_size*0.5, ecal_det_size*0.5 , 0.5 * (window_l));
//  ecalDetWindowS = new G4Box("ecalDetWindowS",ecal_front_face*0.5, ecal_front_face*0.5 , 0.5 * (window_l));
//  ecalDetCaseOuterS=new G4Box("ecalDetCaseOuterS", ecal_det_size+Case_l, ecal_det_size+Case_l, 0.5 * (det_l+window_l)+Case_l);
//  ecalDetCaseInnerS=new G4Box("ecalDetCaseInnerS", ecal_det_size, ecal_det_size , 0.5 * (det_l+window_l)+Case_l);
//  ecalDetCase_fS= new G4SubtractionSolid( "ecalDetCase_fS",ecalDetCaseOuterS,ecalDetCaseInnerS, 0, G4ThreeVector(0,0,Case_l) );
//  ecalDetCase_rS= new G4SubtractionSolid( "ecalDetCase_rS",ecalDetCaseOuterS,ecalDetCaseInnerS, 0, G4ThreeVector(0,0,-1*Case_l) );

  ecalDetCaseOuterS=new G4Box("ecalDetCaseOuterS", ecal_front_face*0.5, ecal_front_face*0.5, 0.5 * (det_l+window_l)+Case_l);
  ecalDetCaseInnerS=new G4Box("ecalDetCaseInnerS", ecal_front_face*0.5, ecal_front_face*0.5, 0.5 * (det_l+window_l)+Case_l);
  ecalDetCase_fS= new G4SubtractionSolid( "ecalDetCase_fS",ecalDetCaseOuterS,ecalDetCaseInnerS, 0, G4ThreeVector(0,0,Case_l) );
  ecalDetCase_rS= new G4SubtractionSolid( "ecalDetCase_rS",ecalDetCaseOuterS,ecalDetCaseInnerS, 0, G4ThreeVector(0,0,-1*Case_l) );



  ecalGapL = new G4LogicalVolume(ecalGapS, PMTGapMaterial, "ecalGapL", 0, 0, 0);
  ecalDetL = new G4LogicalVolume(ecalDetS, DeMaterial, "ecalDetL", 0, 0, 0);
  ecalDetWindowL = new G4LogicalVolume(ecalDetWindowS, WindowMaterial, "ecalDetWindowL", 0, 0, 0);
  ecalDetCase_fL = new G4LogicalVolume(ecalDetCase_fS, WindowMaterial, "ecalDetCase_fL", 0, 0, 0);
  ecalDetCase_rL = new G4LogicalVolume(ecalDetCase_rS, WindowMaterial, "ecalDetCase_rL", 0, 0, 0);
  
  }

  // ECAL physical placement
  const int NECAL_CRYST = 9; //6400; //2500;
  G4VPhysicalVolume *ecalCrystalP_f_wrap[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_r_wrap[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_f_fiber[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_r_fiber[NECAL_CRYST];
  G4VPhysicalVolume *ecalGapP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalGapP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalGapP_rf[NECAL_CRYST];
  G4VPhysicalVolume *ecalGapP_rr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetP_rf[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetP_rr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetWindowP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetWindowP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetWindowP_rf[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetWindowP_rr[NECAL_CRYST];

  G4VPhysicalVolume *ecalDetCaseP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetCaseP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetCaseP_rf[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetCaseP_rr[NECAL_CRYST];

  G4VPhysicalVolume *ecalCrystalP_f_fiberGap[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_r_fiberGap[NECAL_CRYST];


  G4LogicalSkinSurface* FiberWrap_f_Surface = new G4LogicalSkinSurface( "FiberWrap_f_Surface", ecalCrystalL_f_fiber_wrapping, fFiberWrapSurface );
  G4LogicalSkinSurface* FiberWrap_r_Surface = new G4LogicalSkinSurface( "FiberWrap_r_Surface", ecalCrystalL_r_fiber_wrapping, fFiberWrapSurface );
  G4LogicalSkinSurface* PMTSurface = new G4LogicalSkinSurface( "PMTSurface",ecalDetL, fPMTSurface );
  G4LogicalSkinSurface* PMTCase_fSurface = new G4LogicalSkinSurface( "PMTCase_fSurface",ecalDetCase_fL, fPMTCaseSurface );
  G4LogicalSkinSurface* PMTCase_rSurface = new G4LogicalSkinSurface( "PMTCase_rSurface",ecalDetCase_rL, fPMTCaseSurface );
  G4LogicalBorderSurface* FilterSurface_ff[NECAL_CRYST];
  G4LogicalBorderSurface* FilterSurface_fr[NECAL_CRYST];
  G4LogicalBorderSurface* Fiber_scinti_Surface[NECAL_CRYST];
  G4LogicalBorderSurface* Fiber_cherenp_Surface[NECAL_CRYST];
  G4LogicalBorderSurface* Fiber_cherenc_Surface[NECAL_CRYST];
  G4LogicalSkinSurface* crystalSurfacef = new G4LogicalSkinSurface( "crystalSurface",ecalCrystalL_f_fiber, fIdealPolishedOpSurface );
  G4LogicalSkinSurface* crystalSurfacer = new G4LogicalSkinSurface( "crystalSurface",ecalCrystalL_r_fiber, fIdealPolishedOpSurface );

  char name[60];
  G4double x_pos[NECAL_CRYST];
  G4double y_pos[NECAL_CRYST];
  int nArrayECAL = (int)sqrt(NECAL_CRYST);

  int iCrystal;
  for (int iX = 0; iX < nArrayECAL; iX++)
  {
    for (int iY = 0; iY < nArrayECAL; iY++)
    {

      G4RotationMatrix *piRotEcal = new G4RotationMatrix;
      piRotEcal->rotateX(pointingAngle * deg);
      G4RotationMatrix *piRotEcal_op = new G4RotationMatrix;
      piRotEcal_op->rotateX((pointingAngle+180) * deg);

      iCrystal = nArrayECAL * iX + iY;
      x_pos[iCrystal] = (iX - nArrayECAL / 2) * (ecal_front_face + alveola_thickness); // position the baricenter of crystals and then rotating them by

      y_pos[iCrystal] = (iY - nArrayECAL / 2) * (ecal_front_face / cos(pointingAngle * deg) + alveola_thickness);
      cout << " x_pos [" << iCrystal << "] = " << x_pos[iCrystal] << " :: y_pos[" << iCrystal << "] = " << y_pos[iCrystal] << " :: angle = [" << pointingAngle * iX << ", " << pointingAngle * iY << "] " << endl;

      if(y_pos[iCrystal]>ecal_front_face) x_pos[iCrystal]+= (ecal_front_face + alveola_thickness)*0.5;
      if(y_pos[iCrystal]<-ecal_front_face) x_pos[iCrystal]-= (ecal_front_face + alveola_thickness)*0.5;
if( iCrystal==0 || iCrystal==8  ) continue;

//---------------------- absorb and wrap --------------------------
      //add dream detector
      sprintf(name, "ecalCrystalP_f_wrap_%d", iCrystal);
      ecalCrystalP_f_wrap[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + ecal_front_length * 0.5), ecalCrystalL_f_fiber_wrapping, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalCrystalP_r_wrap_%d", iCrystal);
//      ecalCrystalP_r_wrap[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + dchange_rotation + ecal_front_length + ecal_rear_length * 0.5), ecalCrystalL_r_fiber_wrapping, name, worldLV, false, 0, checkOverlaps);

      sprintf(name, "ecalCrystalP_f_fiberGap_%d", iCrystal);
      ecalCrystalP_f_fiberGap[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + ecal_front_length * 0.5), ecalCrystalL_f_fiberGap, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalCrystalP_r_fiberGap_%d", iCrystal);
//      ecalCrystalP_r_fiberGap[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + dchange_rotation + ecal_front_length + ecal_rear_length * 0.5), ecalCrystalL_r_fiberGap, name, worldLV, false, 0, checkOverlaps);

//---------------------- fibers --------------------------
      sprintf(name, "ecalCrystalP_f_fiber_cherenc_%d", iCrystal);
      ecalCrystalP_f_fiber[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + ecal_front_length * 0.5), ecalCrystalL_f_fiber, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalCrystalP_r_fiber_cherenc_%d", iCrystal);
//      ecalCrystalP_r_fiber[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + dchange_rotation + ecal_front_length + ecal_rear_length * 0.5), ecalCrystalL_r_fiber, name, worldLV, false, 0, checkOverlaps);

//---------------------- detectors at both ends --------------------------
      sprintf(name, "ecalDetP_ff_%d", iCrystal);
      ecalDetP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l+window_l + det_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l+window_l + det_l * 0.5) * cos(pointingAngle * deg)), ecalDetL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetP_rr_%d", iCrystal);
      ecalDetP_rr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_rear_length + gap_l+window_l + det_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l+window_l + det_l * 0.5) * cos(pointingAngle * deg)), ecalDetL, name, worldLV, false, 0, checkOverlaps);
//      cout << " z_pos [" << iCrystal << "] = " << ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l+window_l + det_l * 0.5) * cos(pointingAngle * deg) << " :: z_pos[" << iCrystal << "] = " << ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l+window_l + det_l * 0.5) * cos(pointingAngle * deg) << " " << endl;
//---------------------- PMT windows --------------------------
      sprintf(name, "ecalDetWindowP_ff_%d", iCrystal);
      ecalDetWindowP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l+window_l* 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l+window_l * 0.5) * cos(pointingAngle * deg)), ecalDetWindowL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetWindowP_rr_%d", iCrystal);
      ecalDetWindowP_rr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_rear_length + gap_l+window_l* 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l+window_l* 0.5) * cos(pointingAngle * deg)), ecalDetWindowL, name, worldLV, false, 0, checkOverlaps);
//---------------------- PMT Case --------------------------
      sprintf(name, "ecalDetCaseP_ff_%d", iCrystal);
      ecalDetCaseP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l+0.5 * (det_l+window_l)+Case_l) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l+0.5 * (det_l+window_l)+Case_l) * cos(pointingAngle * deg)), ecalDetCase_fL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetCaseP_rr_%d", iCrystal);
      ecalDetCaseP_rr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_rear_length + gap_l+0.5 * (det_l+window_l)+Case_l) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length* 0.5 + (ecal_front_length* 0.5 + gap_l+0.5 * (det_l+window_l)+Case_l) * cos(pointingAngle * deg)), ecalDetCase_rL, name, worldLV, false, 0, checkOverlaps);
//---------------------- gap between detectors and fibers --------------------------
      sprintf(name, "ecalGapP_ff_%d", iCrystal);
      ecalGapP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l* 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l* 0.5) * cos(pointingAngle * deg)), ecalGapL, name, worldLV, false, iCrystal, checkOverlaps);
      sprintf(name, "ecalGapP_rr_%d", iCrystal);
      ecalGapP_rr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_rear_length + gap_l*0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l* 0.5) * cos(pointingAngle * deg)), ecalGapL, name, worldLV, false, iCrystal, checkOverlaps);
//---------------------- filters --------------------------
  FilterSurface_ff[iCrystal] = new G4LogicalBorderSurface( "FilterSurface_ff", ecalGapP_ff[iCrystal], ecalDetWindowP_ff[iCrystal], fFilterSurface_ff );
  FilterSurface_fr[iCrystal] = new G4LogicalBorderSurface( "FilterSurface_fr", ecalGapP_rr[iCrystal], ecalDetWindowP_rr[iCrystal], fFilterSurface_fr );

    }
  }

  //********************************************
  //  HADRONIC CALORIMETER
  //********************************************

  // HCAL solid
  G4VSolid *hcalAbs_t1_S;
  G4VSolid *hcalAbs_t2_S;

  G4VSolid *solenoidS;
  G4VSolid *hcalTileS;

  hcalAbs_t1_S = new G4Box("hcalAbs_t1_S", hcal_width * 0.5, hcal_width * 0.5, hcalAbs_1_thick * 0.5);
  hcalAbs_t2_S = new G4Box("hcalAbs_t2_S", hcal_width * 0.5, hcal_width * 0.5, hcalAbs_2_thick * 0.5);

  solenoidS = new G4Box("solenoidS", hcal_width * 0.5, hcal_width * 0.5, solenoid_thick * 0.5);
  hcalTileS = new G4Box("hcalTileS", hcalTile_width * 0.5, hcalTile_width * 0.5, hcalTile_thick * 0.5);

  // HCAL logical
  G4LogicalVolume *hcalAbs_t1_LV = new G4LogicalVolume(hcalAbs_t1_S, MyMaterials::Aluminium(), " hcalAbs_t1_LV");
  G4LogicalVolume *hcalAbs_t2_LV = new G4LogicalVolume(hcalAbs_t2_S, MyMaterials::Iron(), " hcalAbs_t2_LV");

  G4LogicalVolume *solenoidLV = new G4LogicalVolume(solenoidS, MyMaterials::Iron(), " solenoidLV");
  G4LogicalVolume *hcalTileLV = new G4LogicalVolume(hcalTileS, MyMaterials::PlasticBC418(), " hcalTileLV");

  // HCAL physical
  G4PVPlacement *hcalAbs_t1_PV;
  G4PVPlacement *hcalAbs_t2_PV;
  G4PVPlacement *solenoidPV;
  G4PVPlacement *hcalTilePV;

  float pos_Z_H1 = pos_Z_DE2 + 25 * mm + dchange_rotation;
  float hLayer_gap = 10 * mm;

  //HCAL 1 (front)
  const int NLAYERS_1 = 3;
  const int NTILES = 25;

  G4double x_pos_tile[625];
  G4double y_pos_tile[625];

  int iTile;
  for (int iX = 0; iX < NTILES; iX++)
  {
    for (int iY = 0; iY < NTILES; iY++)
    {
      iTile = NTILES * iX + iY;
      x_pos_tile[iTile] = iX * hcalTile_width - hcalTile_width * NTILES / 2. + hcalTile_width * 0.5;
      y_pos_tile[iTile] = iY * hcalTile_width - hcalTile_width * NTILES / 2. + hcalTile_width * 0.5;
      //      cout << " x_pos_tile [" << iTile << "] = " << x_pos_tile[iTile] << " :: y_pos_tile[" << iTile << "] = " << y_pos_tile[iTile] << endl;
    }
  }

  if (hcal_width > 0.1)
  {
    for (int iLayer = 0; iLayer < NLAYERS_1; iLayer++)
    {

      float z = pos_Z_H1 + hcalAbs_1_thick * 0.5 + iLayer * (hcalAbs_1_thick + hLayer_gap);

      //passive layer of absorber
//      hcalAbs_t1_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, z), hcalAbs_t1_LV, Form("hcalAbs_t1_%d", iLayer), worldLV, false, 0, checkOverlaps);
      iTile = 0;
      //active layer of scintillator
      for (int iX = 0; iX < NTILES; iX++)
      {
        for (int iY = 0; iY < NTILES; iY++)
        {
          iTile = NTILES * iX + iY;
          sprintf(name, "hcalTile_layer_%d_tile_%d", iLayer, iTile);
//          hcalTilePV = new G4PVPlacement(0, G4ThreeVector(x_pos_tile[iTile], y_pos_tile[iTile], z + (hcalAbs_1_thick + hLayer_gap) * 0.5), hcalTileLV, name, worldLV, false, 0);
        }
      }
    }
  }

  //Solenoid
  if (hcal_width > 0.1)
  {
    for (int iLayer = 0; iLayer < 1; iLayer++)
    {
      float z = pos_Z_H1 + NLAYERS_1 * (hcalAbs_1_thick + hLayer_gap) + solenoid_thick * 0.5;
      //passive layer of solenoid
      //solenoidPV = new G4PVPlacement(0, G4ThreeVector(0, 0, z), solenoidLV, "solenoid_PV", worldLV, false, 0, checkOverlaps);
      //active layer of scintillator
      iTile = 0;
      for (int iX = 0; iX < NTILES; iX++)
      {
        for (int iY = 0; iY < NTILES; iY++)
        {
          iTile = NTILES * iX + iY;
          sprintf(name, "hcalTile_layer_%d_tile_%d", iLayer + NLAYERS_1, iTile);
  //        hcalTilePV = new G4PVPlacement(0, G4ThreeVector(x_pos_tile[iTile], y_pos_tile[iTile], z + (solenoid_thick + hLayer_gap) * 0.5), hcalTileLV, name, worldLV, false, 0);
        }
      }
    }
  }

  //HCAL 2 (rear)
  const int NLAYERS_2 = 11;

  if (hcal_width > 0.1)
  {
    for (int iLayer = 0; iLayer < NLAYERS_2; iLayer++)
    {

      float z = pos_Z_H1 + NLAYERS_1 * (hcalAbs_1_thick + hLayer_gap) + (solenoid_thick + hLayer_gap) + hcalAbs_2_thick * 0.5 + iLayer * (hcalAbs_2_thick + hLayer_gap);

      //passive layer of absorber
//      hcalAbs_t2_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, z), hcalAbs_t2_LV, Form("hcalAbs_t2_%d", iLayer), worldLV, false, 0, checkOverlaps);

      //active layer of scintillator
      iTile = 0;
      for (int iX = 0; iX < NTILES; iX++)
      {
        for (int iY = 0; iY < NTILES; iY++)
        {
          iTile = NTILES * iX + iY;
          sprintf(name, "hcalTile_layer_%d_tile_%d", iLayer + 1 + NLAYERS_1, iTile);
//          hcalTilePV = new G4PVPlacement(0, G4ThreeVector(x_pos_tile[iTile], y_pos_tile[iTile], z + (hcalAbs_2_thick + hLayer_gap) * 0.5), hcalTileLV, name, worldLV, false, 0);
        }
      }
    }
  }

  //-----------------------------------------------------
  //------------- sensitive detector & surface --------------
  //-----------------------------------------------------

  auto sdman = G4SDManager::GetSDMpointer();
  auto fDRDetSD = new DR_PMTSD("/DR_Det");
  sdman->AddNewDetector(fDRDetSD);
  ecalDetL->SetSensitiveDetector(fDRDetSD);


  //-----------------------------------------------------
  //------------- Visualization attributes --------------
  //-----------------------------------------------------

  G4Colour white(1.00, 1.00, 1.00); // white
  G4Colour gray(0.50, 0.50, 0.50);  // gray
  G4Colour black(0.00, 0.00, 0.00); // black
  G4Colour red(1.00, 0.00, 0.00);   // red
  G4Colour green(0.00, 1.00, 0.00); // green
  G4Colour blue(0.00, 0.00, 1.00);  // blue
  G4Colour lightblue(0.678, 0.847, 0.902);
  G4Colour cyan(0.00, 1.00, 1.00);    // cyan
  G4Colour air(0.90, 0.94, 1.00);     // cyan
  G4Colour magenta(1.00, 0.00, 1.00); // magenta
  G4Colour yellow(1.00, 1.00, 0.00);  // yellow
  G4Colour brass(0.80, 0.60, 0.40);   // brass
  G4Colour brown(0.70, 0.40, 0.10);   // brown

  G4VisAttributes *VisAttWorld = new G4VisAttributes(white);
  VisAttWorld->SetVisibility(true);
  VisAttWorld->SetForceWireframe(true);
  worldLV->SetVisAttributes(VisAttWorld);

  G4VisAttributes *VisAttCore = new G4VisAttributes(green);
  VisAttCore->SetVisibility(false);
  VisAttCore->SetForceWireframe(true);
  coreLV->SetVisAttributes(VisAttCore);

  G4VisAttributes *VisCrystalCore = new G4VisAttributes(brass);
  VisCrystalCore->SetVisibility(false);
  VisCrystalCore->SetForceWireframe(true);
  //  ecalCrystalL_f -> SetVisAttributes(VisCrystalCore);


  G4VisAttributes *VisAttGap = new G4VisAttributes(blue);
  VisAttGap->SetVisibility(false);
  VisAttGap->SetForceWireframe(true);
  gapLV->SetVisAttributes(VisAttGap);

  G4VisAttributes *VisAttecalDet = new G4VisAttributes(gray);
  VisAttecalDet->SetVisibility(true);
  VisAttecalDet->SetForceWireframe(true);
  ecalDetL->SetVisAttributes(VisAttecalDet);
//  ecalDetL_f->SetVisAttributes(VisAttecalDet);
//  ecalDetL_r->SetVisAttributes(VisAttecalDet);
//  ecalDetL->SetVisAttributes(VisAttecalDet);
  ecalGapL->SetVisAttributes(VisAttecalDet);
  //dream detector

  G4VisAttributes *VisAttecalDetW = new G4VisAttributes(red);
  VisAttecalDetW->SetVisibility(true);
  VisAttecalDetW->SetForceWireframe(true);
  ecalDetWindowL->SetVisAttributes(VisAttecalDetW);

  G4VisAttributes *Visabs = new G4VisAttributes(brass);
  Visabs->SetVisibility(true);
  Visabs->SetForceWireframe(true);
//  ecalCrystalL_f_absorb->SetVisAttributes(Visabs);
//  ecalCrystalL_r_absorb->SetVisAttributes(Visabs);

  G4VisAttributes *Viswrap = new G4VisAttributes(cyan);
  Viswrap->SetVisibility(true);
  //Viswrap->SetForceWireframe(true);
  Viswrap->SetForceSolid(true);
  ecalCrystalL_f_fiber_wrapping->SetVisAttributes(Viswrap);
  ecalCrystalL_r_fiber_wrapping->SetVisAttributes(Viswrap);

  G4VisAttributes *VisFiber = new G4VisAttributes(lightblue);
  VisFiber->SetVisibility(true);
  VisFiber->SetForceWireframe(false);
  ecalCrystalL_f_fiber->SetVisAttributes(VisFiber);
  ecalCrystalL_r_fiber->SetVisAttributes(VisFiber);
  //
  //
  G4VisAttributes *VisHCALTile = new G4VisAttributes(brass);
  VisHCALTile->SetVisibility(true);
  VisHCALTile->SetForceWireframe(true);
  hcalTileLV->SetVisAttributes(VisHCALTile);

  G4VisAttributes *VisHCALAbs1 = new G4VisAttributes(gray);
  VisHCALAbs1->SetVisibility(true);
  hcalAbs_t1_LV->SetVisAttributes(VisHCALAbs1);

  G4VisAttributes *VisHCALAbs2 = new G4VisAttributes(red);
  VisHCALAbs2->SetVisibility(true);
  hcalAbs_t2_LV->SetVisAttributes(VisHCALAbs2);

  G4VisAttributes *VisAttTracker = new G4VisAttributes(gray);
  VisAttTracker->SetVisibility(false);
  VisAttTracker->SetForceWireframe(false);
  trackerLV->SetVisAttributes(VisAttTracker);
  servicesLV->SetVisAttributes(VisAttTracker);

  G4VisAttributes *VisAttDet = new G4VisAttributes(gray);
  VisAttDet->SetVisibility(false);
  VisAttDet->SetForceWireframe(false);
  detLV->SetVisAttributes(VisAttDet);

  if (B_field_intensity > 0.1 * tesla)
    ConstructField();

  // //-----------------------------------------------
  // //------------- Fast photon timing --------------
  // //-----------------------------------------------

  // std::vector<std::pair<double,double> > rIndVecCore;
  // std::vector<std::pair<double,double> > rIndVecClad;
  // std::vector<std::pair<double,double> > rIndVecAir;
  // std::vector<std::pair<double,double> > rIndVecGap;

  // G4MaterialPropertyVector* mpVec;

  // mpVec = ClMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  // for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
  // {
  //   std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  //   rIndVecCore.push_back(dummy);
  // }

  // mpVec = WoMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  // for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
  // {
  //   std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  //   std::pair<double,double> dummy2(mpVec->GetLowEdgeEnergy(it)/eV,fibre_cladRIndex);
  //   rIndVecAir.push_back(dummy);
  //   rIndVecClad.push_back(dummy2);
  // }

  // mpVec = GaMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
  // for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
  // {
  //   std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  //   rIndVecGap.push_back(dummy);
  // }

  // fib = FiberInit(bar_length,fibre_radius,CreateTree::Instance()->attLengths,rIndVecCore,rIndVecClad,rIndVecAir,rIndVecGap) ;

  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl;
  return worldPV;
}
/*
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
G4VSolid *DetectorConstruction::ConstructHollowWrap_f() const
{

  G4Box *ecalCrystalS_f_outwrap_s;
  G4Box *ecalCrystalS_f_innerwrap_s;

  G4Tubs *ecalCrystalS_f_outwrap_c;
  G4Tubs *ecalCrystalS_f_innerwrap_c;

  if (fiber_type == 0)
  {
    ecalCrystalS_f_outwrap_c = new G4Tubs("ecalCrystalS_f_outwrap_c", 0, fiber_diameter * 0.5 + wrapping_thick + depth, ecal_front_length * 0.5 + gap_l, 0, 360 * deg);
    ecalCrystalS_f_innerwrap_c = new G4Tubs("ecalCrystalS_f_innerwrap_c", 0, fiber_diameter * 0.5 + depth, ecal_front_length, 0, 360 * deg);
    return new G4SubtractionSolid("ecalCrystalS_f_hollowwrap", ecalCrystalS_f_outwrap_c, ecalCrystalS_f_innerwrap_c);
  }
  else
  {

    ecalCrystalS_f_outwrap_s = new G4Box("ecalCrystalS_f_outwrap_s", fiber_diameter * 0.5 + wrapping_thick + depth, fiber_diameter * 0.5 + wrapping_thick + depth, ecal_front_length * 0.5 + gap_l);
    ecalCrystalS_f_innerwrap_s = new G4Box("ecalCrystalS_f_innerwrap_s", fiber_diameter * 0.5 + depth, fiber_diameter * 0.5 + depth, ecal_front_length);
    return new G4SubtractionSolid("ecalCrystalS_f_hollowwrap", ecalCrystalS_f_outwrap_s, ecalCrystalS_f_innerwrap_s);
  }
}

G4VSolid *DetectorConstruction::ConstructHollowWrap_r() const
{

  G4Box *ecalCrystalS_r_outwrap_s;
  G4Box *ecalCrystalS_r_innerwrap_s;

  G4Tubs *ecalCrystalS_r_outwrap_c;
  G4Tubs *ecalCrystalS_r_innerwrap_c;

  if (fiber_type == 0)
  {
    ecalCrystalS_r_outwrap_c = new G4Tubs("ecalCrystalS_r_outwrap_c", 0, fiber_diameter * 0.5 + wrapping_thick + depth, ecal_rear_length * 0.5 + gap_l, 0, 360 * deg);
    ecalCrystalS_r_innerwrap_c = new G4Tubs("ecalCrystalS_r_innerwrap_c", 0, fiber_diameter * 0.5 + depth, ecal_rear_length, 0, 360 * deg);
    return new G4SubtractionSolid("ecalCrystalS_r_hollowwrap", ecalCrystalS_r_outwrap_c, ecalCrystalS_r_innerwrap_c);
  }
  else
  {

    ecalCrystalS_r_outwrap_s = new G4Box("ecalCrystalS_r_outwrap_s", fiber_diameter * 0.5 + wrapping_thick + depth, fiber_diameter * 0.5 + wrapping_thick + depth, ecal_rear_length * 0.5 + gap_l);
    ecalCrystalS_r_innerwrap_s = new G4Box("ecalCrystalS_r_innerwrap_s", fiber_diameter * 0.5 + depth, fiber_diameter * 0.5 + depth, ecal_rear_length);
    return new G4SubtractionSolid("ecalCrystalS_r_hollowwrap", ecalCrystalS_r_outwrap_s, ecalCrystalS_r_innerwrap_s);
  }

}
*/
void DetectorConstruction::initializeMaterials()
{
  //-----------------
  // define materials

  WoMaterial = NULL;
  if (world_material == 0)
    WoMaterial = MyMaterials::Vacuum();
  else if (world_material == 1)
    WoMaterial = MyMaterials::Air();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre world material specifier " << world_material << G4endl;
    exit(-1);
  }
  G4cout << "Wo. material: " << WoMaterial << G4endl;

  CoMaterial = NULL;
  if (core_material == 1)
    CoMaterial = MyMaterials::Quartz();
  else if (core_material == 2)
    CoMaterial = MyMaterials::SiO2();
  else if (core_material == 3)
    CoMaterial = MyMaterials::SiO2_Ce();
  else if (core_material == 4)
    CoMaterial = MyMaterials::LuAG_Ce();
  else if (core_material == 5)
    CoMaterial = MyMaterials::YAG_Ce();
  else if (core_material == 6)
    CoMaterial = MyMaterials::LSO();
  else if (core_material == 7)
    CoMaterial = MyMaterials::LYSO();
  else if (core_material == 8)
    CoMaterial = MyMaterials::LuAG_undoped();
  else if (core_material == 9)
    CoMaterial = MyMaterials::GAGG_Ce();
  else if (core_material == 11)
    CoMaterial = MyMaterials::LuAG_Pr();
  else if (core_material == 12)
    CoMaterial = MyMaterials::PbF2();
  else if (core_material == 13)
    CoMaterial = MyMaterials::PlasticBC408();
  else if (core_material == 14)
    CoMaterial = MyMaterials::PlasticBC418();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << core_material << G4endl;
    exit(-1);
  }
  G4cout << "Co. material: " << CoMaterial << G4endl;

  EcalMaterial = NULL;
  if (ecal_material == 1)
    EcalMaterial = MyMaterials::Quartz();
  else if (ecal_material == 2)
    EcalMaterial = MyMaterials::SiO2();
  else if (ecal_material == 3)
    EcalMaterial = MyMaterials::SiO2_Ce();
  else if (ecal_material == 4)
    EcalMaterial = MyMaterials::LuAG_Ce();
  else if (ecal_material == 5)
    EcalMaterial = MyMaterials::YAG_Ce();
  else if (ecal_material == 6)
    EcalMaterial = MyMaterials::LSO();
  else if (ecal_material == 7)
    EcalMaterial = MyMaterials::LYSO();
  else if (ecal_material == 8)
    EcalMaterial = MyMaterials::LuAG_undoped();
  else if (ecal_material == 9)
    EcalMaterial = MyMaterials::GAGG_Ce();
  else if (ecal_material == 10)
    EcalMaterial = MyMaterials::LuAG_Pr();
  else if (ecal_material == 11)
    EcalMaterial = MyMaterials::PbF2();
  else if (ecal_material == 12)
    EcalMaterial = MyMaterials::PlasticBC408();
  else if (ecal_material == 13)
    EcalMaterial = MyMaterials::PlasticBC418();
  else if (ecal_material == 14)
    EcalMaterial = MyMaterials::PWO();
  else if (ecal_material == 15)
    EcalMaterial = MyMaterials::Acrylic();
  else if (ecal_material == 16)
    EcalMaterial = MyMaterials::copper();
  else if (ecal_material == 17)
    EcalMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << ecal_material << G4endl;
    exit(-1);
  }
  G4cout << "ECAL material: " << EcalMaterial << G4endl;

  /************************************************************************************/
  ScintiMaterial = NULL;
  if (scinti_material == 1)
    ScintiMaterial = MyMaterials::Quartz();
  else if (scinti_material == 2)
    ScintiMaterial = MyMaterials::SiO2();
  else if (scinti_material == 3)
    ScintiMaterial = MyMaterials::SiO2_Ce();
  else if (scinti_material == 4)
    ScintiMaterial = MyMaterials::LuAG_Ce();
  else if (scinti_material == 5)
    ScintiMaterial = MyMaterials::YAG_Ce();
  else if (scinti_material == 6)
    ScintiMaterial = MyMaterials::LSO();
  else if (scinti_material == 7)
    ScintiMaterial = MyMaterials::LYSO();
  else if (scinti_material == 8)
    ScintiMaterial = MyMaterials::LuAG_undoped();
  else if (scinti_material == 9)
    ScintiMaterial = MyMaterials::GAGG_Ce();
  else if (scinti_material == 10)
    ScintiMaterial = MyMaterials::LuAG_Pr();
  else if (scinti_material == 11)
    ScintiMaterial = MyMaterials::PbF2();
  else if (scinti_material == 12)
    ScintiMaterial = MyMaterials::PlasticBC408();
  else if (scinti_material == 13)
    ScintiMaterial = MyMaterials::PlasticBC418();
  else if (scinti_material == 14)
    ScintiMaterial = MyMaterials::PWO();
  else if (scinti_material == 15)
    ScintiMaterial = MyMaterials::Acrylic();
  else if (scinti_material == 16)
    ScintiMaterial = MyMaterials::copper();
  else if (scinti_material == 17)
    ScintiMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << scinti_material << G4endl;
    exit(-1);
  }

  CherencMaterial = NULL;
  if (Cherenc_material == 1)
    CherencMaterial = MyMaterials::Quartz();
  else if (Cherenc_material == 2)
    CherencMaterial = MyMaterials::SiO2();
  else if (Cherenc_material == 3)
    CherencMaterial = MyMaterials::SiO2_Ce();
  else if (Cherenc_material == 4)
    CherencMaterial = MyMaterials::LuAG_Ce();
  else if (Cherenc_material == 5)
    CherencMaterial = MyMaterials::YAG_Ce();
  else if (Cherenc_material == 6)
    CherencMaterial = MyMaterials::LSO();
  else if (Cherenc_material == 7)
    CherencMaterial = MyMaterials::LYSO();
  else if (Cherenc_material == 8)
    CherencMaterial = MyMaterials::LuAG_undoped();
  else if (Cherenc_material == 9)
    CherencMaterial = MyMaterials::GAGG_Ce();
  else if (Cherenc_material == 10)
    CherencMaterial = MyMaterials::LuAG_Pr();
  else if (Cherenc_material == 11)
    CherencMaterial = MyMaterials::PbF2();
  else if (Cherenc_material == 12)
    CherencMaterial = MyMaterials::PlasticBC408();
  else if (Cherenc_material == 13)
    CherencMaterial = MyMaterials::PlasticBC418();
  else if (Cherenc_material == 14)
    CherencMaterial = MyMaterials::PWO();
  else if (Cherenc_material == 15)
    CherencMaterial = MyMaterials::Acrylic();
  else if (Cherenc_material == 16)
    CherencMaterial = MyMaterials::copper();
  else if (Cherenc_material == 17)
    CherencMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Cherenc_material << G4endl;
    exit(-1);
  }

  CherenpMaterial = NULL;
  if (Cherenp_material == 1)
    CherenpMaterial = MyMaterials::Quartz();
  else if (Cherenp_material == 2)
    CherenpMaterial = MyMaterials::SiO2();
  else if (Cherenp_material == 3)
    CherenpMaterial = MyMaterials::SiO2_Ce();
  else if (Cherenp_material == 4)
    CherenpMaterial = MyMaterials::LuAG_Ce();
  else if (Cherenp_material == 5)
    CherenpMaterial = MyMaterials::YAG_Ce();
  else if (Cherenp_material == 6)
    CherenpMaterial = MyMaterials::LSO();
  else if (Cherenp_material == 7)
    CherenpMaterial = MyMaterials::LYSO();
  else if (Cherenp_material == 8)
    CherenpMaterial = MyMaterials::LuAG_undoped();
  else if (Cherenp_material == 9)
    CherenpMaterial = MyMaterials::GAGG_Ce();
  else if (Cherenp_material == 10)
    CherenpMaterial = MyMaterials::LuAG_Pr();
  else if (Cherenp_material == 11)
    CherenpMaterial = MyMaterials::PbF2();
  else if (Cherenp_material == 12)
    CherenpMaterial = MyMaterials::PlasticBC408();
  else if (Cherenp_material == 13)
    CherenpMaterial = MyMaterials::PlasticBC418();
  else if (Cherenp_material == 14)
    CherenpMaterial = MyMaterials::PWO();
  else if (Cherenp_material == 15)
    CherenpMaterial = MyMaterials::Acrylic();
  else if (Cherenp_material == 16)
    CherenpMaterial = MyMaterials::copper();
  else if (Cherenp_material == 17)
    CherenpMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Cherenp_material << G4endl;
    exit(-1);
  }

  WrapMaterial = NULL;
  if (wrap_material == 1)
    WrapMaterial = MyMaterials::Quartz();
  else if (wrap_material == 2)
    WrapMaterial = MyMaterials::SiO2();
  else if (wrap_material == 3)
    WrapMaterial = MyMaterials::SiO2_Ce();
  else if (wrap_material == 4)
    WrapMaterial = MyMaterials::LuAG_Ce();
  else if (wrap_material == 5)
    WrapMaterial = MyMaterials::YAG_Ce();
  else if (wrap_material == 6)
    WrapMaterial = MyMaterials::LSO();
  else if (wrap_material == 7)
    WrapMaterial = MyMaterials::LYSO();
  else if (wrap_material == 8)
    WrapMaterial = MyMaterials::LuAG_undoped();
  else if (wrap_material == 9)
    WrapMaterial = MyMaterials::GAGG_Ce();
  else if (wrap_material == 10)
    WrapMaterial = MyMaterials::LuAG_Pr();
  else if (wrap_material == 11)
    WrapMaterial = MyMaterials::PbF2();
  else if (wrap_material == 12)
    WrapMaterial = MyMaterials::PlasticBC408();
  else if (wrap_material == 13)
    WrapMaterial = MyMaterials::PlasticBC418();
  else if (wrap_material == 14)
    WrapMaterial = MyMaterials::PWO();
  else if (wrap_material == 15)
    WrapMaterial = MyMaterials::Acrylic();
  else if (wrap_material == 16)
    WrapMaterial = MyMaterials::copper();
  else if (wrap_material == 17)
    WrapMaterial = MyMaterials::Epoxy();
  else if (wrap_material == 18)
    WrapMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Al" );
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << wrap_material << G4endl;
    exit(-1);
  }

  /************************************************************************************/
  WindowMaterial = MyMaterials::PyrexGlass();
  PMTGapMaterial = MyMaterials::silicone();
  GaMaterial = NULL;
  if (gap_material == 1)
    GaMaterial = MyMaterials::Air();
  else if (gap_material == 2)
    GaMaterial = MyMaterials::OpticalGrease();
  else if (gap_material == 3)
    GaMaterial = MyMaterials::MeltMount168();
  else if (gap_material == 4)
    GaMaterial = MyMaterials::OpticalGrease155();
  else if (gap_material == 5)
    GaMaterial = MyMaterials::silicone();
  else if (gap_material == 6)
    GaMaterial = MyMaterials::PyrexGlass();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid gap material specifier " << gap_material << G4endl;
    exit(-1);
  }
  G4cout << "Gap material: " << gap_material << G4endl;

  DeMaterial = NULL;
  if (det_material == 1)
    DeMaterial = MyMaterials::Silicon();
  else if (det_material == 2)
    DeMaterial = MyMaterials::Quartz();
  else if (det_material == 3)
    DeMaterial = MyMaterials::Air();
  else if (det_material == 4)
    DeMaterial = MyMaterials::Bialkali();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid detector material specifier " << det_material << G4endl;
    exit(-1);
  }
  G4cout << "Detector material: " << det_material << G4endl;

  //------------------
  // change properties

  if (core_absLength > 0)
  {
    const G4int nEntries_ABS = 2;
    G4double PhotonEnergy_ABS[nEntries_ABS] = {1. * eV, 10. * eV};
    G4double Absorption[nEntries_ABS] = {core_absLength * mm, core_absLength * mm};

    CoMaterial->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
    CoMaterial->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH", PhotonEnergy_ABS, Absorption, nEntries_ABS);
  }
  if (core_rIndex > 0)
  {
    const G4int nEntries_RI = 2;
    G4double PhotonEnergy_RI[nEntries_RI] = {1. * eV, 10. * eV};
    G4double RefractiveIndex[nEntries_RI] = {core_rIndex, core_rIndex};

    CoMaterial->GetMaterialPropertiesTable()->RemoveProperty("RINDEX");
    CoMaterial->GetMaterialPropertiesTable()->AddProperty("RINDEX", PhotonEnergy_RI, RefractiveIndex, nEntries_RI);
  }
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void DetectorConstruction::ConstructField()
{
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::begin <<<<<<" << G4endl;

  static G4TransportationManager *trMgr = G4TransportationManager::GetTransportationManager();

  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager *globalFieldMgr = trMgr->GetFieldManager();

  if (!B_field_IsInitialized)
  {
    // magnetic field parallel to the beam direction (w/ tilt)
    G4ThreeVector fieldVector(0.0522 * B_field_intensity, 0.0522 * B_field_intensity, 0.9973 * B_field_intensity);

    B_field = new G4UniformMagField(fieldVector);
    globalFieldMgr->SetDetectorField(B_field);
    globalFieldMgr->CreateChordFinder(B_field);
    globalFieldMgr->GetChordFinder()->SetDeltaChord(0.005 * mm);
    B_field_IsInitialized = true;
  }

  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::end <<< " << G4endl;
  return;
}

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit) && (maxStep > 0.))
    stepLimit->SetMaxAllowedStep(maxStep);
}

// //-----------------------------------------------
// //------------- Surface properties --------------
// //-----------------------------------------------

void DetectorConstruction::initializeSurface()
{
  fFiberWrapSurface = MakeS_wrap();
  fPMTSurface = MakeS_PMT();
  fPMTCaseSurface = MakeS_IdealWhiteSurface();

  static const unsigned nentries = 2;
  static double phoE[nentries]   = {0.01*eV, 10.0*eV};
  double reflectivity[nentries]  = {wrap_ref, wrap_ref};
  G4MaterialPropertiesTable* table = fFiberWrapSurface->GetMaterialPropertiesTable();
  if( table ){
    table->RemoveProperty( "REFLECTIVITY" );
    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  } else {
    table = new G4MaterialPropertiesTable();
    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
    fFiberWrapSurface->SetMaterialPropertiesTable( table );
  }
  fIdealPolishedOpSurface = MakeS_IdealPolished();
/*
//perfect reflection
  static const unsigned nentriesc = 2;
  static double phoEc[nentries]   = {0.01*eV, 10.0*eV};
  double reflectivityc[nentries]  = {1,1};
  G4MaterialPropertiesTable* tablecrys = fIdealPolishedOpSurface->GetMaterialPropertiesTable();
  if( tablecrys ){
    tablecrys->RemoveProperty( "REFLECTIVITY" );
    tablecrys->AddProperty( "REFLECTIVITY", phoEc, reflectivityc, nentriesc );
  } else {
    tablecrys = new G4MaterialPropertiesTable();
    table->AddProperty( "REFLECTIVITY", phoEc, reflectivityc, nentriesc );
    fIdealPolishedOpSurface->SetMaterialPropertiesTable( tablecrys );
  }
*/


  fFilterSurface_ff = MakeS_IdealPolished();
  fFilterSurface_fr = MakeS_IdealPolished();

  const G4int nEntries_tran_u330 = 18;
  G4double PhotonEnergy_tran_u330[nEntries_tran_u330] = {4.103235582*eV, 4.048912761*eV,3.944470896     *eV,3.81331407      *eV,3.750947295     *eV,3.690598697     *eV,3.603625797     *eV,3.534216059     *eV,3.456978466     *eV,3.380554073     *eV,3.309819592     *eV,3.230572067     *eV,3.185706353     *eV,3.131340814     *eV,3.087086539     *eV,3.050146549     *eV,2.992445212     *eV,2.933127681     *eV};
  G4double transIndex_u330[nEntries_tran_u330] = {0.201372        ,0.202705        ,0.211043        ,0.227125        ,0.234102        ,0.233987        ,0.235942        ,0.235798        ,0.229958        ,0.219856        ,0.199831        ,0.155664        ,0.115833        ,0.068881        ,0.037554        ,0.017592        ,0.00466 ,0.000935        };

  const G4int nEntries_tran_ug5 = 25;
  G4double PhotonEnergy_tran_ug5[nEntries_tran_ug5] = {4.092143963	*eV,4.034407045	*eV,3.943981544	*eV,3.825267479	*eV,3.743879546	*eV,3.62234533	*eV,3.530100421	*eV,3.414187953	*eV,3.300875562	*eV,3.233225815	*eV,3.168293273	*eV,3.137871012	*eV,3.099604675	*eV,3.074608111	*eV,3.041899835	*eV,3.001980276	*eV,2.947821354	*eV,2.873756177	*eV,2.764363413	*eV,2.697530944	*eV,2.642988727	*eV,2.564470256	*eV,2.529030177	*eV,2.498643446	*eV,2.482374634	*eV};
  G4double transIndex_ug5[nEntries_tran_ug5] = {0.197435	,0.199462	,0.209962	,0.224666	,0.230184	,0.23489	,0.238226	,0.232947	,0.219825	,0.204682	,0.176017	,0.153858	,0.125979	,0.106681	,0.089491	,0.06658	,0.04575	,0.031245	,0.029397	,0.029832	,0.02319	,0.016426	,0.012699	,0.009703	,0.008199	};


    G4MaterialPropertiesTable* Filtertable_u330 = new G4MaterialPropertiesTable();
    Filtertable_u330->AddProperty("REFLECTIVITY", PhotonEnergy_tran_u330, transIndex_u330, nEntries_tran_u330);
    G4MaterialPropertiesTable* Filtertable_ug5 = new G4MaterialPropertiesTable();
    Filtertable_ug5->AddProperty("REFLECTIVITY", PhotonEnergy_tran_ug5, transIndex_ug5, nEntries_tran_ug5);

    if(front_filter== -1){
      std::cout<<"No front Filter!"<<std::endl;
    }else if(front_filter==0){
      fFilterSurface_ff->SetMaterialPropertiesTable( Filtertable_u330 );
    }else if(front_filter==1){
      fFilterSurface_ff->SetMaterialPropertiesTable( Filtertable_ug5 );
    }
    if(rear_filter== -1){
      std::cout<<"No rear Filter!"<<std::endl;
    }else if(rear_filter==0){
      fFilterSurface_fr->SetMaterialPropertiesTable( Filtertable_u330 );
    }else if(rear_filter==1){
      fFilterSurface_fr->SetMaterialPropertiesTable( Filtertable_ug5 );
    }

  return;
}


/*
// Initialization classes
//
void DetectorConstruction::initializeSurface(G4OpticalSurface *mySurface, string surfaceType)
{
    if(surfaceType == "crystal") 
    {
//         cout << "CRISTALLO " << crystalSurfinish << endl;

        surfinish   	= crystalSurfinish;
        RefFile     	= cReffile;
        reflectivity    = cReflectivity;
        surrefind   	= cSurrefind;
        surtype     	= cSurtype;
        specularspike   = cSpecularspike;
        specularlobe    = cSpecularlobe;
        sigmaalpha  	= cSigmaalpha;
        backscatter     = cBackscatter;
        lambertian  	= cLambertian;
    }



    if(this->surfinish <= 5)
    {
        G4cout << "Using unified model." << G4endl;
        mySurface -> SetModel(unified);
        switch(this->surtype) {
        case 0:
            mySurface -> SetType(dielectric_metal);
            G4cout << "Surface type: dielectric_metal" << G4endl;
            break;
        case 1:
            mySurface -> SetType(dielectric_dielectric);
            G4cout << "Surface type: dielectric_dielectric" << G4endl;
            break;
        }
    }

    else if(this->surfinish > 5 && surfaceType == "wrapping") G4cout << "Value not allowed" << G4endl;

    else {
        G4cout << "Using LUT for surface treatment." << G4endl;
        mySurface -> SetModel(LUT);
        mySurface -> SetType(dielectric_LUT);
    }

    switch(this->surfinish) {
    case 0:
        mySurface -> SetFinish(polished);
        G4cout << "Surface finish: polished" << G4endl;
        break;
    case 1:
        mySurface -> SetFinish(polishedfrontpainted);
        G4cout << "Surface finish: polishedfrontpainted" << G4endl;
        break;
    case 2:
        mySurface -> SetFinish(polishedbackpainted);
        G4cout << "Surface finish: polishedbackpainted" << G4endl;
        break;
    case 3:
        mySurface -> SetFinish(ground);
        G4cout << "Surface finish: ground" << G4endl;
        break;
    case 4:
        mySurface -> SetFinish(groundfrontpainted);
        G4cout << "Surface finish: groundfrontpainted" << G4endl;
        break;
    case 5:
        mySurface -> SetFinish(groundbackpainted);
        G4cout << "Surface finish: groundbackpainted" << G4endl;
        break;
    case 17:
        mySurface -> SetFinish(polishedteflonair);
        G4cout << "Surface finish: polishedteflonair" << G4endl;
        break;
    case 18:
        mySurface -> SetFinish(polishedtioair);
        G4cout << "Surface finish: polishedtioair" << G4endl;
        break;
    case 26:
        mySurface -> SetFinish(etchedtioair);
        G4cout << "Surface finish: etchedtioair" << G4endl;
        break;
    case 34:
        mySurface -> SetFinish(groundtioair);
        G4cout << "Surface finish: groundtioair" << G4endl;
        break;
    case 36:
        mySurface -> SetFinish(polishedtyvekair);
        G4cout << "Surface finish: polishedtyvekair" << G4endl;
        break;
    default:
        G4cout << "Surface finish unkown!" << G4endl;
        exit(0);
    }
}*/

//
// reflectivity
//
/*
void DetectorConstruction::initializeReflectivitySurface(G4OpticalSurface *surface, string surfaceType)
{

    int NumRefl = 0;
    G4double EphotonRefl[1000];
    G4double Refl[1000];
    
    EphotonRefl[0] = 0.0001 * eV;
    EphotonRefl[1] = 1.0 * eV;
    EphotonRefl[2] = 4.08 * eV;
    Refl[0] = 0.0; // suppress photons with energy < 1eV (will not be detected)
    Refl[1] = this->crystal_reflectivity;
    Refl[2] = this->crystal_reflectivity;
    NumRefl = 100;

    G4cout << "Reflectivities as a function of the photon energy:" << G4endl;
    for(int i = 0; i < NumRefl; i++)         G4cout << i << "   " << EphotonRefl[i] << "   " << Refl[i] << G4endl;
    

    Ephoton[0] = 0.0001 * eV;
    Ephoton[1] = 1.0 * eV;
    Ephoton[2] = 4.08 * eV;
    G4double RefractiveIndex[3] = {this->surrefind, this->surrefind, this->surrefind};
    G4double SpecularLobe[3]    = {this->specularlobe, this->specularlobe, this->specularlobe};
    G4double SpecularSpike[3]   = {this->specularspike, this->specularspike, this->specularspike};
    G4double Backscatter[3]     = {this->backscatter, this->backscatter, this->backscatter};
    G4double Lambertian[3]      = {this->lambertian, this->lambertian, this->lambertian};
    G4MaterialPropertiesTable *myST = new G4MaterialPropertiesTable();
    G4cout << "Read from config-file: " << G4endl;
    G4cout << "Read SPECULARLOBECONSTANT  : " << SpecularLobe[0] << G4endl;
    G4cout << "Read SPECULARSPIKECONSTANT : " << SpecularSpike[0] << G4endl;
    G4cout << "Read BACKSCATTERCONSTANT   : " << Backscatter[0] << G4endl;
    G4cout << "Read LAMBERTIAN            : " << Lambertian[0] << G4endl;
    G4cout << "Read ref. index            : " << RefractiveIndex[0] << G4endl;

    myST->AddProperty("RINDEX", Ephoton, RefractiveIndex, 3);
    if(this->specularlobe >= 0) {
        G4cout << "Setting SPECULARLOBECONSTANT to : " << SpecularLobe[0] << G4endl;
        myST->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    3);
    }
    if(this->specularspike >= 0) {
        G4cout << "Setting SPECULARSPIKECONSTANT to : " << SpecularSpike[0] << G4endl;
        myST->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   3);
    }
    if(this->backscatter >= 0) {
        G4cout << "Setting BACKSCATTERCONSTANT to : " << Backscatter[0] << G4endl;
        myST->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     3);
    }
    if(this->lambertian >= 0) {
        G4cout << "Setting LAMBERTIAN to : " << Lambertian[0] << G4endl;
        myST->AddProperty("LAMBERTIAN",            Ephoton, Lambertian,      3);
    }

//     setWrappingIndices(myST,1.82,"teflon");

    surface->SetMaterialPropertiesTable(myST);
    if(this->sigmaalpha >= 0) surface->SetSigmaAlpha(this->sigmaalpha);


}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
