#ifndef CreateTree_H
#define CreateTree_H 1

#include <iostream>
#include <vector>
#include "TString.h"
#include <map>

#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"

class CreateTree
{
private:
  TTree *ftree;
  TString fname;

public:
  CreateTree(TString name);
  ~CreateTree();

  TTree *GetTree() const { return ftree; };
  TString GetName() const { return fname; };
  void AddEnergyDeposit(int index, float deposit);
  void AddScintillationPhoton(int index);
  void AddCerenkovPhoton(int index);
  int Fill();
  bool Write(TFile *);
  void Clear();

  static CreateTree *Instance() { return fInstance; };
  static CreateTree *fInstance;

  int Event;

  int inputTrackerX0;
  int inputServiceAlmm;
  int inputTimingThick;
  int inputE1Thick;
  int inputE2Thick;
  int inputE1Width;
  int inputTimingECAL_dist;
float ParentEK1;
float ParentEK2;
float ParentEtot1;
float ParentEtot2;

float ParentnonionEdep;
float ParenttotEdep;

  std::vector<int> *primary_second_id;
  std::vector<float> *primary_second_Ek;
  std::vector<float> *primary_second_Et;
  std::vector<float> *primary_second_mass;
/*
  std::vector<float> *primary_Ek_Et_deuteron;
  std::vector<float> *primary_Ek_Et_alpha;
  std::vector<float> *primary_Ek_Et_neutron;
  std::vector<float> *primary_Ek_Et_proton;
  std::vector<float> *primary_Ek_Et_electron;
  std::vector<float> *primary_Ek_Et_positron;
  std::vector<float> *primary_Ek_Et_gamma;
  std::vector<float> *primary_Ek_Et_mup;
  std::vector<float> *primary_Ek_Et_mun;
  std::vector<float> *primary_Ek_Et_pip;
  std::vector<float> *primary_Ek_Et_pin;
  std::vector<float> *primary_Ek_Et_pi0;
  std::vector<float> *primary_Ek_Et_kaonp;
  std::vector<float> *primary_Ek_Et_kaonn;
  std::vector<float> *primary_Ek_Et_kaon0L;
  std::vector<float> *primary_Ek_Et_eta;
  std::vector<float> *primary_Ek_Et_other;
*/
  std::vector<float> *inputMomentum;        // Px Py Pz E
  std::vector<float> *inputInitialPosition; // x, y, z

  std::vector<float> *primaryMomT1; // Px Py Pz E
  std::vector<float> *primaryPosT1; // x, y, z

  std::vector<float> *primaryMomE1; // Px Py Pz E
  std::vector<float> *primaryPosE1; // x, y, z

  int nTracksT1;
  int nTracksT2;
  int nTracksE1;
  int nTracksE2;
  int nTracksTRK[6];

  //integrated energy in each longitudinal layer
  float depositedEnergyEscapeWorld;
  float depositedEnergyTotal;
  float depositedEnergyTiming_f;
  float depositedEnergyTiming_r;
  float depositedEnergyECAL_f[3];
  float depositedEnergyECAL_r[3];
  float depositedEnergyHCALAct;
  float depositedEnergyHCALPas;
  float depositedEnergyServices;
  float depositedEnergyTimingGap;
  float depositedEnergyEcalGap;
  float depositedEnergyEcalDet;
  float depositedEnergySolenoid;
  float depositedEnergyWorld;

  float depositedIonEnergyTotal;
  float depositedIonEnergyTiming_f;
  float depositedIonEnergyTiming_r;
  float depositedIonEnergyECAL_f[3];
  float depositedIonEnergyECAL_r[3];
  float depositedIonEnergyHCALAct;
  float depositedIonEnergyHCALPas;
  float depositedIonEnergyServices;
  float depositedIonEnergyTimingGap;
  float depositedIonEnergyEcalGap;
  float depositedIonEnergyEcalDet;
  float depositedIonEnergySolenoid;
  float depositedIonEnergyWorld;

  float depositedElecEnergyTotal;
  float depositedElecEnergyTiming_f;
  float depositedElecEnergyTiming_r;
  float depositedElecEnergyECAL_f[3];
  float depositedElecEnergyECAL_r[3];
  float depositedElecEnergyHCALAct;
  float depositedElecEnergyHCALPas;
  float depositedElecEnergyServices;
  float depositedElecEnergyTimingGap;
  float depositedElecEnergyEcalGap;
  float depositedElecEnergyEcalDet;
  float depositedElecEnergySolenoid;
  float depositedElecEnergyWorld;

  //store the energy deposition by components

  float depositedEnergyECAL_absorb_f_particleID[8];
  float depositedEnergyECAL_absorb_r_particleID[8];
  float depositedIonEnergyECAL_absorb_f_particleID[8];
  float depositedIonEnergyECAL_absorb_r_particleID[8];

  float depositedEnergyECAL_scinti_f_particleID[8];
  float depositedEnergyECAL_scinti_r_particleID[8];
  float depositedIonEnergyECAL_scinti_f_particleID[8];
  float depositedIonEnergyECAL_scinti_r_particleID[8];

  float depositedEnergyECAL_cheren_f_particleID[8];
  float depositedEnergyECAL_cheren_r_particleID[8];
  float depositedIonEnergyECAL_cheren_f_particleID[8];
  float depositedIonEnergyECAL_cheren_r_particleID[8];

  int tot_phot_cer_Timing_f_total;
  int tot_phot_cer_Timing_r_total;
  int ECAL_f_total_S;
  int ECAL_r_total_S;
  int ECAL_f_total_C;
  int ECAL_r_total_C;
  int tot_phot_cer_ECAL_scinti_f_particleID[8];
  int tot_phot_cer_ECAL_scinti_r_particleID[8];
  int tot_phot_cer_ECAL_cheren_f_total;
  int tot_phot_cer_ECAL_cheren_r_total;
  int tot_phot_cer_ECAL_cheren_f_particleID[8];
  int tot_phot_cer_ECAL_cheren_r_particleID[8];
  int tot_phot_cer_HCAL;

//  int tot_phot_cer_SDdetected_ff; 
  int tot_phot_cer_SDdetected_fr;
  int tot_phot_cer_SDdetected_rf;
//  int tot_phot_cer_SDdetected_rr;
  int SDdetected_ff_S;
  int SDdetected_ff_C;
  int SDdetected_rr_S;
  int SDdetected_rr_C;
  int tot_phot_cer_should_det1;
  int tot_phot_cer_should_det2;
  int tot_phot_cer_should_det3;
  int tot_phot_cer_should_det4;
  /***************** begin to seperate energy into different channels    ******************/
  float Edep_Tracker_layer[6];

  //energy deposit in each trasnversally segmented channel
  float Edep_Timing_f_ch[18];
  float Edep_Timing_r_ch[18];

  float Edep_ECAL_f_ch[6400];
  float Edep_ECAL_r_ch[6400];

  float IonEdep_ECAL_f_ch[6400];
  float IonEdep_ECAL_r_ch[6400];

  //  float ElecEdep_ECAL_f_ch[6400];
  //  float ElecEdep_ECAL_r_ch[6400];

  //  float Edep_ECAL_f_scinti_ch[6400];
  //  float Edep_ECAL_r_scinti_ch[6400];
  //  float Edep_ECAL_f_cheren_ch[6400];
  //  float Edep_ECAL_r_cheren_ch[6400];

  //  float IonEdep_ECAL_f_scinti_ch[6400];
  //  float IonEdep_ECAL_r_scinti_ch[6400];
  //  float IonEdep_ECAL_f_cheren_ch[6400];
  //  float IonEdep_ECAL_r_cheren_ch[6400];

  //  float ElecEdep_ECAL_f_scinti_ch[6400];
  //  float ElecEdep_ECAL_r_scinti_ch[6400];
  //  float ElecEdep_ECAL_f_cheren_ch[6400];
  //  float ElecEdep_ECAL_r_cheren_ch[6400];

  float E_Zdep_0to5000mm_total[2500];
  float E_Tdep_0to5ns_total[2500];
  float E_Zdep_0to5000mm_Pion_n[2500];
  float E_Tdep_0to5ns_Pion_n[2500];
  float E_Zdep_0to5000mm_Positron[2500];
  float E_Tdep_0to5ns_Positron[2500];
  float E_Zdep_0to5000mm_Electron[2500];
  float E_Tdep_0to5ns_Electron[2500];
  float E_Zdep_0to5000mm_Photon[2500];
  float E_Tdep_0to5ns_Photon[2500];
  float E_Zdep_0to5000mm_Pion_p[2500];
  float E_Tdep_0to5ns_Pion_p[2500];
  float E_Zdep_0to5000mm_Kion[2500];
  float E_Tdep_0to5ns_Kion[2500];
  float E_Zdep_0to5000mm_Neutron[2500];
  float E_Tdep_0to5ns_Neutron[2500];
  float E_Zdep_0to5000mm_Proton[2500];
  float E_Tdep_0to5ns_Proton[2500];

  /***************** get the energy distribution in time and space   ******************/
  /*
  TH2F *h_time_z_st;
  TH2F *h_time_z_mt;
  TH2F *h_time_z_lt;

  TH2F *h_localtime_z_st;
  TH2F *h_localtime_z_mt;
  TH2F *h_localtime_z_lt;

  TH2F *h_localtime_z_st_en;
  TH2F *h_localtime_z_mt_en;
  TH2F *h_localtime_z_lt_en;

  TH2F *h_localtime_z_st_ep;
  TH2F *h_localtime_z_mt_ep;
  TH2F *h_localtime_z_lt_ep;

  TH2F *h_localtime_z_st_gamma;
  TH2F *h_localtime_z_mt_gamma;
  TH2F *h_localtime_z_lt_gamma;

  TH2F *h_localtime_z_st_n;
  TH2F *h_localtime_z_mt_n;
  TH2F *h_localtime_z_lt_n;

  TH2F *h_localtime_z_st_p;
  TH2F *h_localtime_z_mt_p;
  TH2F *h_localtime_z_lt_p;

  TH2F *h_localtime_z_st_pin;
  TH2F *h_localtime_z_mt_pin;
  TH2F *h_localtime_z_lt_pin;

  TH2F *h_localtime_z_st_pip;
  TH2F *h_localtime_z_mt_pip;
  TH2F *h_localtime_z_lt_pip;

  TH1F *h_parentID_Edep;

  TH1F *h_phot_cer_lambda_Timing_f;
  TH1F *h_phot_cer_lambda_Timing_r;
  TH1F *h_phot_cer_lambda_ECAL_f;
  TH1F *h_phot_cer_lambda_ECAL_r;
  TH1F *h_phot_cer_lambda_HCAL;
  TH1F *h_phot_cer_parentID;
*/
  TH1F *h_phot_cer_lambda_ECAL_f_Ceren;
  TH1F *h_phot_cer_lambda_ECAL_r_Ceren;
  TH1F *h_phot_cer_lambda_ECAL_f_produce_Ceren;
  TH1F *h_phot_cer_lambda_ECAL_r_produce_Ceren;
  TH2F *h_photon_2D_produce_Ceren;
  TH2F *h_photon_2D_receive_Ceren;

  TH1F *h_signal;

  TH1F *h_phot_cer_lambda_ECAL_f_Scin;
  TH1F *h_phot_cer_lambda_ECAL_r_Scin;
  TH1F *h_phot_cer_lambda_ECAL_f_produce_Scin;
  TH1F *h_phot_cer_lambda_ECAL_r_produce_Scin;
  TH2F *h_photon_2D_produce_Scin;
  TH2F *h_photon_2D_receive_Scin;
float EnergyECAL_r_enter;
float EnergyECAL_r_leave;

float EnergyECAL_r_enterK;
float EnergyECAL_r_leaveK;

  std::vector<int> *pdgid_escape;
  std::vector<float> *energy_escape;
  std::vector<float> *positionx_escape; 
  std::vector<float> *positiony_escape;
  std::vector<float> *positionz_escape;

};

#endif
