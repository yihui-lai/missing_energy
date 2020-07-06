//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul  6 14:48:08 2020 by ROOT version 6.16/00
// from TTree tree/tree
// found on file: pion_5GeV.root
//////////////////////////////////////////////////////////

#ifndef missing_h
#define missing_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class missing {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event;
   Float_t         inputTrackerX0;
   Float_t         inputServiceAlmm;
   Float_t         inputTimingThick;
   Float_t         inputE1Thick;
   Float_t         inputE2Thick;
   Float_t         inputE1Width;
   Float_t         inputTimingECAL_dist;
   vector<int>     *primary_second_id;
   vector<float>   *primary_second_Ek;
   vector<float>   *primary_second_Et;
   vector<float>   *primary_second_mass;
   Float_t         ParentEK1;
   Float_t         ParentEK2;
   Float_t         ParentEtot1;
   Float_t         ParentEtot2;
   Float_t         ParentnonionEdep;
   Float_t         ParenttotEdep;
   vector<float>   *inputInitialPosition;
   vector<float>   *inputMomentum;
   vector<float>   *primaryPosT1;
   vector<float>   *primaryMomT1;
   vector<float>   *primaryPosE1;
   vector<float>   *primaryMomE1;
   Int_t           nTracksT1;
   Int_t           nTracksT2;
   Int_t           nTracksE1;
   Int_t           nTracksE2;
   Float_t         nTracksTRK[6];
   Float_t         depositedEnergyTotal;
   Float_t         depositedEnergyEscapeWorld;
   Float_t         depositedEnergyTiming_f;
   Float_t         depositedEnergyTiming_r;
   Float_t         depositedEnergyECAL_f[3];
   Float_t         depositedEnergyECAL_r[3];
   Float_t         depositedEnergyHCALAct;
   Float_t         depositedEnergyHCALPas;
   Float_t         depositedEnergyWorld;
   Float_t         depositedEnergyTimingGap;
   Float_t         depositedEnergyServices;
   Float_t         depositedEnergyEcalGap;
   Float_t         depositedEnergyEcalDet;
   Float_t         depositedEnergySolenoid;
   Float_t         depositedEnergyECAL_absorb_f_particleID[8];
   Float_t         depositedEnergyECAL_absorb_r_particleID[8];
   Float_t         depositedEnergyECAL_scinti_f_particleID[8];
   Float_t         depositedEnergyECAL_scinti_r_particleID[8];
   Float_t         depositedEnergyECAL_cheren_f_particleID[8];
   Float_t         depositedEnergyECAL_cheren_r_particleID[8];
   Float_t         depositedIonEnergyTotal;
   Float_t         depositedIonEnergyTiming_f;
   Float_t         depositedIonEnergyTiming_r;
   Float_t         depositedIonEnergyECAL_f[3];
   Float_t         depositedIonEnergyECAL_r[3];
   Float_t         depositedIonEnergyHCALAct;
   Float_t         depositedIonEnergyHCALPas;
   Float_t         depositedIonEnergyWorld;
   Float_t         depositedIonEnergyTimingGap;
   Float_t         depositedIonEnergyServices;
   Float_t         depositedIonEnergyEcalGap;
   Float_t         depositedIonEnergyEcalDet;
   Float_t         depositedIonEnergySolenoid;
   Float_t         depositedIonEnergyECAL_absorb_f_particleID[8];
   Float_t         depositedIonEnergyECAL_absorb_r_particleID[8];
   Float_t         depositedIonEnergyECAL_scinti_f_particleID[8];
   Float_t         depositedIonEnergyECAL_scinti_r_particleID[8];
   Float_t         depositedIonEnergyECAL_cheren_f_particleID[8];
   Float_t         depositedIonEnergyECAL_cheren_r_particleID[8];
   Float_t         depositedElecEnergyTotal;
   Float_t         depositedElecEnergyTiming_f;
   Float_t         depositedElecEnergyTiming_r;
   Float_t         depositedElecEnergyECAL_f[3];
   Float_t         depositedElecEnergyECAL_r[3];
   Float_t         depositedElecEnergyHCALAct;
   Float_t         depositedElecEnergyHCALPas;
   Float_t         depositedElecEnergyWorld;
   Float_t         depositedElecEnergyTimingGap;
   Float_t         depositedElecEnergyServices;
   Float_t         depositedElecEnergyEcalGap;
   Float_t         depositedElecEnergyEcalDet;
   Float_t         depositedElecEnergySolenoid;
   Float_t         Edep_Tracker_layer[6];
   Float_t         Edep_Timing_f_ch[18];
   Float_t         Edep_Timing_r_ch[18];
   Float_t         Edep_ECAL_f_ch[6400];
   Float_t         Edep_ECAL_r_ch[6400];
   Float_t         IonEdep_ECAL_f_ch[6400];
   Float_t         IonEdep_ECAL_r_ch[6400];
   Int_t           tot_phot_cer_Timing_f_total;
   Int_t           tot_phot_cer_Timing_r_total;
   Int_t           ECAL_f_total_S;
   Int_t           ECAL_f_total_C;
   Int_t           ECAL_r_total_S;
   Int_t           ECAL_r_total_C;
   Int_t           tot_phot_cer_ECAL_scinti_f_particleID[8];
   Int_t           tot_phot_cer_ECAL_scinti_r_particleID[8];
   Int_t           tot_phot_cer_ECAL_cheren_f_particleID[8];
   Int_t           tot_phot_cer_ECAL_cheren_r_particleID[8];
   Int_t           tot_phot_cer_HCAL;
   Int_t           SDdetected_ff_S;
   Int_t           SDdetected_ff_C;
   Int_t           tot_phot_cer_SDdetected_fr;
   Int_t           tot_phot_cer_SDdetected_rf;
   Int_t           SDdetected_rr_S;
   Int_t           SDdetected_rr_C;
   Int_t           tot_phot_cer_should_det1;
   Int_t           tot_phot_cer_should_det2;
   Int_t           tot_phot_cer_should_det3;
   Int_t           tot_phot_cer_should_det4;
   Float_t         E_Zdep_0to5000mm_total[2500];
   Float_t         E_Zdep_0to5000mm_Pion_n[2500];
   Float_t         E_Zdep_0to5000mm_Positron[2500];
   Float_t         E_Zdep_0to5000mm_Electron[2500];
   Float_t         E_Zdep_0to5000mm_Photon[2500];
   Float_t         E_Zdep_0to5000mm_Pion_p[2500];
   Float_t         E_Zdep_0to5000mm_Kion[2500];
   Float_t         E_Zdep_0to5000mm_Neutron[2500];
   Float_t         E_Zdep_0to5000mm_Proton[2500];
   Float_t         E_Tdep_0to5ns_total[2500];
   Float_t         E_Tdep_0to5ns_Pion_n[2500];
   Float_t         E_Tdep_0to5ns_Positron[2500];
   Float_t         E_Tdep_0to5ns_Electron[2500];
   Float_t         E_Tdep_0to5ns_Photon[2500];
   Float_t         E_Tdep_0to5ns_Pion_p[2500];
   Float_t         E_Tdep_0to5ns_Kion[2500];
   Float_t         E_Tdep_0to5ns_Neutron[2500];
   Float_t         E_Tdep_0to5ns_Proton[2500];
   Float_t         EnergyECAL_r_enter;
   Float_t         EnergyECAL_r_leave;
   Float_t         EnergyECAL_r_enterK;
   Float_t         EnergyECAL_r_leaveK;
   vector<int>     *pdgid_escape;
   vector<float>   *energy_escape;
   vector<float>   *positionx_escape;
   vector<float>   *positiony_escape;
   vector<float>   *positionz_escape;

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_inputTrackerX0;   //!
   TBranch        *b_inputServiceAlmm;   //!
   TBranch        *b_inputTimingThick;   //!
   TBranch        *b_inputE1Thick;   //!
   TBranch        *b_inputE2Thick;   //!
   TBranch        *b_inputE1Width;   //!
   TBranch        *b_inputTimingECAL_dist;   //!
   TBranch        *b_primary_second_id;   //!
   TBranch        *b_primary_second_Ek;   //!
   TBranch        *b_primary_second_Et;   //!
   TBranch        *b_primary_second_mass;   //!
   TBranch        *b_ParentEK1;   //!
   TBranch        *b_ParentEK2;   //!
   TBranch        *b_ParentEtot1;   //!
   TBranch        *b_ParentEtot2;   //!
   TBranch        *b_ParentnonionEdep;   //!
   TBranch        *b_ParenttotEdep;   //!
   TBranch        *b_inputInitialPosition;   //!
   TBranch        *b_inputMomentum;   //!
   TBranch        *b_primaryPosT1;   //!
   TBranch        *b_primaryMomT1;   //!
   TBranch        *b_primaryPosE1;   //!
   TBranch        *b_primaryMomE1;   //!
   TBranch        *b_nTracksT1;   //!
   TBranch        *b_nTracksT2;   //!
   TBranch        *b_nTracksE1;   //!
   TBranch        *b_nTracksE2;   //!
   TBranch        *b_nTracksTRK;   //!
   TBranch        *b_depositedEnergyTotal;   //!
   TBranch        *b_depositedEnergyEscapeWorld;   //!
   TBranch        *b_depositedEnergyTiming_f;   //!
   TBranch        *b_depositedEnergyTiming_r;   //!
   TBranch        *b_depositedEnergyECAL_f;   //!
   TBranch        *b_depositedEnergyECAL_r;   //!
   TBranch        *b_depositedEnergyHCALAct;   //!
   TBranch        *b_depositedEnergyHCALPas;   //!
   TBranch        *b_depositedEnergyWorld;   //!
   TBranch        *b_depositedEnergyTimingGap;   //!
   TBranch        *b_depositedEnergyServices;   //!
   TBranch        *b_depositedEnergyEcalGap;   //!
   TBranch        *b_depositedEnergyEcalDet;   //!
   TBranch        *b_depositedEnergySolenoid;   //!
   TBranch        *b_depositedEnergyECAL_absorb_f_particleID;   //!
   TBranch        *b_depositedEnergyECAL_absorb_r_particleID;   //!
   TBranch        *b_depositedEnergyECAL_scinti_f_particleID;   //!
   TBranch        *b_depositedEnergyECAL_scinti_r_particleID;   //!
   TBranch        *b_depositedEnergyECAL_cheren_f_particleID;   //!
   TBranch        *b_depositedEnergyECAL_cheren_r_particleID;   //!
   TBranch        *b_depositedIonEnergyTotal;   //!
   TBranch        *b_depositedIonEnergyTiming_f;   //!
   TBranch        *b_depositedIonEnergyTiming_r;   //!
   TBranch        *b_depositedIonEnergyECAL_f;   //!
   TBranch        *b_depositedIonEnergyECAL_r;   //!
   TBranch        *b_depositedIonEnergyHCALAct;   //!
   TBranch        *b_depositedIonEnergyHCALPas;   //!
   TBranch        *b_depositedIonEnergyWorld;   //!
   TBranch        *b_depositedIonEnergyTimingGap;   //!
   TBranch        *b_depositedIonEnergyServices;   //!
   TBranch        *b_depositedIonEnergyEcalGap;   //!
   TBranch        *b_depositedIonEnergyEcalDet;   //!
   TBranch        *b_depositedIonEnergySolenoid;   //!
   TBranch        *b_depositedIonEnergyECAL_absorb_f_particleID;   //!
   TBranch        *b_depositedIonEnergyECAL_absorb_r_particleID;   //!
   TBranch        *b_depositedIonEnergyECAL_scinti_f_particleID;   //!
   TBranch        *b_depositedIonEnergyECAL_scinti_r_particleID;   //!
   TBranch        *b_depositedIonEnergyECAL_cheren_f_particleID;   //!
   TBranch        *b_depositedIonEnergyECAL_cheren_r_particleID;   //!
   TBranch        *b_depositedElecEnergyTotal;   //!
   TBranch        *b_depositedElecEnergyTiming_f;   //!
   TBranch        *b_depositedElecEnergyTiming_r;   //!
   TBranch        *b_depositedElecEnergyECAL_f;   //!
   TBranch        *b_depositedElecEnergyECAL_r;   //!
   TBranch        *b_depositedElecEnergyHCALAct;   //!
   TBranch        *b_depositedElecEnergyHCALPas;   //!
   TBranch        *b_depositedElecEnergyWorld;   //!
   TBranch        *b_depositedElecEnergyTimingGap;   //!
   TBranch        *b_depositedElecEnergyServices;   //!
   TBranch        *b_depositedElecEnergyEcalGap;   //!
   TBranch        *b_depositedElecEnergyEcalDet;   //!
   TBranch        *b_depositedElecEnergySolenoid;   //!
   TBranch        *b_Edep_Tracker_layer;   //!
   TBranch        *b_Edep_Timing_f_ch;   //!
   TBranch        *b_Edep_Timing_r_ch;   //!
   TBranch        *b_Edep_ECAL_f_ch;   //!
   TBranch        *b_Edep_ECAL_r_ch;   //!
   TBranch        *b_IonEdep_ECAL_f_ch;   //!
   TBranch        *b_IonEdep_ECAL_r_ch;   //!
   TBranch        *b_tot_phot_cer_Timing_f_total;   //!
   TBranch        *b_tot_phot_cer_Timing_r_total;   //!
   TBranch        *b_ECAL_f_total_S;   //!
   TBranch        *b_ECAL_f_total_C;   //!
   TBranch        *b_ECAL_r_total_S;   //!
   TBranch        *b_ECAL_r_total_C;   //!
   TBranch        *b_tot_phot_cer_ECAL_scinti_f_particleID;   //!
   TBranch        *b_tot_phot_cer_ECAL_scinti_r_particleID;   //!
   TBranch        *b_tot_phot_cer_ECAL_cheren_f_particleID;   //!
   TBranch        *b_tot_phot_cer_ECAL_cheren_r_particleID;   //!
   TBranch        *b_tot_phot_cer_HCAL;   //!
   TBranch        *b_SDdetected_ff_S;   //!
   TBranch        *b_SDdetected_ff_C;   //!
   TBranch        *b_tot_phot_cer_SDdetected_fr;   //!
   TBranch        *b_tot_phot_cer_SDdetected_rf;   //!
   TBranch        *b_SDdetected_rr_S;   //!
   TBranch        *b_SDdetected_rr_C;   //!
   TBranch        *b_tot_phot_cer_should_det1;   //!
   TBranch        *b_tot_phot_cer_should_det2;   //!
   TBranch        *b_tot_phot_cer_should_det3;   //!
   TBranch        *b_tot_phot_cer_should_det4;   //!
   TBranch        *b_E_Zdep_0to5000mm_total;   //!
   TBranch        *b_E_Zdep_0to5000mm_Pion_n;   //!
   TBranch        *b_E_Zdep_0to5000mm_Positron;   //!
   TBranch        *b_E_Zdep_0to5000mm_Electron;   //!
   TBranch        *b_E_Zdep_0to5000mm_Photon;   //!
   TBranch        *b_E_Zdep_0to5000mm_Pion_p;   //!
   TBranch        *b_E_Zdep_0to5000mm_Kion;   //!
   TBranch        *b_E_Zdep_0to5000mm_Neutron;   //!
   TBranch        *b_E_Zdep_0to5000mm_Proton;   //!
   TBranch        *b_E_Tdep_0to5ns_total;   //!
   TBranch        *b_E_Tdep_0to5ns_Pion_n;   //!
   TBranch        *b_E_Tdep_0to5ns_Positron;   //!
   TBranch        *b_E_Tdep_0to5ns_Electron;   //!
   TBranch        *b_E_Tdep_0to5ns_Photon;   //!
   TBranch        *b_E_Tdep_0to5ns_Pion_p;   //!
   TBranch        *b_E_Tdep_0to5ns_Kion;   //!
   TBranch        *b_E_Tdep_0to5ns_Neutron;   //!
   TBranch        *b_E_Tdep_0to5ns_Proton;   //!
   TBranch        *b_EnergyECAL_r_enter;   //!
   TBranch        *b_EnergyECAL_r_leave;   //!
   TBranch        *b_EnergyECAL_r_enterK;   //!
   TBranch        *b_EnergyECAL_r_leaveK;   //!
   TBranch        *b_pdgid_escape;   //!
   TBranch        *b_energy_escape;   //!
   TBranch        *b_positionx_escape;   //!
   TBranch        *b_positiony_escape;   //!
   TBranch        *b_positionz_escape;   //!

   missing(TTree *tree=0);
   virtual ~missing();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef missing_cxx
missing::missing(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pion_5GeV.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pion_5GeV.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

missing::~missing()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t missing::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t missing::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void missing::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   primary_second_id = 0;
   primary_second_Ek = 0;
   primary_second_Et = 0;
   primary_second_mass = 0;
   inputInitialPosition = 0;
   inputMomentum = 0;
   primaryPosT1 = 0;
   primaryMomT1 = 0;
   primaryPosE1 = 0;
   primaryMomE1 = 0;
   pdgid_escape = 0;
   energy_escape = 0;
   positionx_escape = 0;
   positiony_escape = 0;
   positionz_escape = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("inputTrackerX0", &inputTrackerX0, &b_inputTrackerX0);
   fChain->SetBranchAddress("inputServiceAlmm", &inputServiceAlmm, &b_inputServiceAlmm);
   fChain->SetBranchAddress("inputTimingThick", &inputTimingThick, &b_inputTimingThick);
   fChain->SetBranchAddress("inputE1Thick", &inputE1Thick, &b_inputE1Thick);
   fChain->SetBranchAddress("inputE2Thick", &inputE2Thick, &b_inputE2Thick);
   fChain->SetBranchAddress("inputE1Width", &inputE1Width, &b_inputE1Width);
   fChain->SetBranchAddress("inputTimingECAL_dist", &inputTimingECAL_dist, &b_inputTimingECAL_dist);
   fChain->SetBranchAddress("primary_second_id", &primary_second_id, &b_primary_second_id);
   fChain->SetBranchAddress("primary_second_Ek", &primary_second_Ek, &b_primary_second_Ek);
   fChain->SetBranchAddress("primary_second_Et", &primary_second_Et, &b_primary_second_Et);
   fChain->SetBranchAddress("primary_second_mass", &primary_second_mass, &b_primary_second_mass);
   fChain->SetBranchAddress("ParentEK1", &ParentEK1, &b_ParentEK1);
   fChain->SetBranchAddress("ParentEK2", &ParentEK2, &b_ParentEK2);
   fChain->SetBranchAddress("ParentEtot1", &ParentEtot1, &b_ParentEtot1);
   fChain->SetBranchAddress("ParentEtot2", &ParentEtot2, &b_ParentEtot2);
   fChain->SetBranchAddress("ParentnonionEdep", &ParentnonionEdep, &b_ParentnonionEdep);
   fChain->SetBranchAddress("ParenttotEdep", &ParenttotEdep, &b_ParenttotEdep);
   fChain->SetBranchAddress("inputInitialPosition", &inputInitialPosition, &b_inputInitialPosition);
   fChain->SetBranchAddress("inputMomentum", &inputMomentum, &b_inputMomentum);
   fChain->SetBranchAddress("primaryPosT1", &primaryPosT1, &b_primaryPosT1);
   fChain->SetBranchAddress("primaryMomT1", &primaryMomT1, &b_primaryMomT1);
   fChain->SetBranchAddress("primaryPosE1", &primaryPosE1, &b_primaryPosE1);
   fChain->SetBranchAddress("primaryMomE1", &primaryMomE1, &b_primaryMomE1);
   fChain->SetBranchAddress("nTracksT1", &nTracksT1, &b_nTracksT1);
   fChain->SetBranchAddress("nTracksT2", &nTracksT2, &b_nTracksT2);
   fChain->SetBranchAddress("nTracksE1", &nTracksE1, &b_nTracksE1);
   fChain->SetBranchAddress("nTracksE2", &nTracksE2, &b_nTracksE2);
   fChain->SetBranchAddress("nTracksTRK", nTracksTRK, &b_nTracksTRK);
   fChain->SetBranchAddress("depositedEnergyTotal", &depositedEnergyTotal, &b_depositedEnergyTotal);
   fChain->SetBranchAddress("depositedEnergyEscapeWorld", &depositedEnergyEscapeWorld, &b_depositedEnergyEscapeWorld);
   fChain->SetBranchAddress("depositedEnergyTiming_f", &depositedEnergyTiming_f, &b_depositedEnergyTiming_f);
   fChain->SetBranchAddress("depositedEnergyTiming_r", &depositedEnergyTiming_r, &b_depositedEnergyTiming_r);
   fChain->SetBranchAddress("depositedEnergyECAL_f", depositedEnergyECAL_f, &b_depositedEnergyECAL_f);
   fChain->SetBranchAddress("depositedEnergyECAL_r", depositedEnergyECAL_r, &b_depositedEnergyECAL_r);
   fChain->SetBranchAddress("depositedEnergyHCALAct", &depositedEnergyHCALAct, &b_depositedEnergyHCALAct);
   fChain->SetBranchAddress("depositedEnergyHCALPas", &depositedEnergyHCALPas, &b_depositedEnergyHCALPas);
   fChain->SetBranchAddress("depositedEnergyWorld", &depositedEnergyWorld, &b_depositedEnergyWorld);
   fChain->SetBranchAddress("depositedEnergyTimingGap", &depositedEnergyTimingGap, &b_depositedEnergyTimingGap);
   fChain->SetBranchAddress("depositedEnergyServices", &depositedEnergyServices, &b_depositedEnergyServices);
   fChain->SetBranchAddress("depositedEnergyEcalGap", &depositedEnergyEcalGap, &b_depositedEnergyEcalGap);
   fChain->SetBranchAddress("depositedEnergyEcalDet", &depositedEnergyEcalDet, &b_depositedEnergyEcalDet);
   fChain->SetBranchAddress("depositedEnergySolenoid", &depositedEnergySolenoid, &b_depositedEnergySolenoid);
   fChain->SetBranchAddress("depositedEnergyECAL_absorb_f_particleID", depositedEnergyECAL_absorb_f_particleID, &b_depositedEnergyECAL_absorb_f_particleID);
   fChain->SetBranchAddress("depositedEnergyECAL_absorb_r_particleID", depositedEnergyECAL_absorb_r_particleID, &b_depositedEnergyECAL_absorb_r_particleID);
   fChain->SetBranchAddress("depositedEnergyECAL_scinti_f_particleID", depositedEnergyECAL_scinti_f_particleID, &b_depositedEnergyECAL_scinti_f_particleID);
   fChain->SetBranchAddress("depositedEnergyECAL_scinti_r_particleID", depositedEnergyECAL_scinti_r_particleID, &b_depositedEnergyECAL_scinti_r_particleID);
   fChain->SetBranchAddress("depositedEnergyECAL_cheren_f_particleID", depositedEnergyECAL_cheren_f_particleID, &b_depositedEnergyECAL_cheren_f_particleID);
   fChain->SetBranchAddress("depositedEnergyECAL_cheren_r_particleID", depositedEnergyECAL_cheren_r_particleID, &b_depositedEnergyECAL_cheren_r_particleID);
   fChain->SetBranchAddress("depositedIonEnergyTotal", &depositedIonEnergyTotal, &b_depositedIonEnergyTotal);
   fChain->SetBranchAddress("depositedIonEnergyTiming_f", &depositedIonEnergyTiming_f, &b_depositedIonEnergyTiming_f);
   fChain->SetBranchAddress("depositedIonEnergyTiming_r", &depositedIonEnergyTiming_r, &b_depositedIonEnergyTiming_r);
   fChain->SetBranchAddress("depositedIonEnergyECAL_f", depositedIonEnergyECAL_f, &b_depositedIonEnergyECAL_f);
   fChain->SetBranchAddress("depositedIonEnergyECAL_r", depositedIonEnergyECAL_r, &b_depositedIonEnergyECAL_r);
   fChain->SetBranchAddress("depositedIonEnergyHCALAct", &depositedIonEnergyHCALAct, &b_depositedIonEnergyHCALAct);
   fChain->SetBranchAddress("depositedIonEnergyHCALPas", &depositedIonEnergyHCALPas, &b_depositedIonEnergyHCALPas);
   fChain->SetBranchAddress("depositedIonEnergyWorld", &depositedIonEnergyWorld, &b_depositedIonEnergyWorld);
   fChain->SetBranchAddress("depositedIonEnergyTimingGap", &depositedIonEnergyTimingGap, &b_depositedIonEnergyTimingGap);
   fChain->SetBranchAddress("depositedIonEnergyServices", &depositedIonEnergyServices, &b_depositedIonEnergyServices);
   fChain->SetBranchAddress("depositedIonEnergyEcalGap", &depositedIonEnergyEcalGap, &b_depositedIonEnergyEcalGap);
   fChain->SetBranchAddress("depositedIonEnergyEcalDet", &depositedIonEnergyEcalDet, &b_depositedIonEnergyEcalDet);
   fChain->SetBranchAddress("depositedIonEnergySolenoid", &depositedIonEnergySolenoid, &b_depositedIonEnergySolenoid);
   fChain->SetBranchAddress("depositedIonEnergyECAL_absorb_f_particleID", depositedIonEnergyECAL_absorb_f_particleID, &b_depositedIonEnergyECAL_absorb_f_particleID);
   fChain->SetBranchAddress("depositedIonEnergyECAL_absorb_r_particleID", depositedIonEnergyECAL_absorb_r_particleID, &b_depositedIonEnergyECAL_absorb_r_particleID);
   fChain->SetBranchAddress("depositedIonEnergyECAL_scinti_f_particleID", depositedIonEnergyECAL_scinti_f_particleID, &b_depositedIonEnergyECAL_scinti_f_particleID);
   fChain->SetBranchAddress("depositedIonEnergyECAL_scinti_r_particleID", depositedIonEnergyECAL_scinti_r_particleID, &b_depositedIonEnergyECAL_scinti_r_particleID);
   fChain->SetBranchAddress("depositedIonEnergyECAL_cheren_f_particleID", depositedIonEnergyECAL_cheren_f_particleID, &b_depositedIonEnergyECAL_cheren_f_particleID);
   fChain->SetBranchAddress("depositedIonEnergyECAL_cheren_r_particleID", depositedIonEnergyECAL_cheren_r_particleID, &b_depositedIonEnergyECAL_cheren_r_particleID);
   fChain->SetBranchAddress("depositedElecEnergyTotal", &depositedElecEnergyTotal, &b_depositedElecEnergyTotal);
   fChain->SetBranchAddress("depositedElecEnergyTiming_f", &depositedElecEnergyTiming_f, &b_depositedElecEnergyTiming_f);
   fChain->SetBranchAddress("depositedElecEnergyTiming_r", &depositedElecEnergyTiming_r, &b_depositedElecEnergyTiming_r);
   fChain->SetBranchAddress("depositedElecEnergyECAL_f", depositedElecEnergyECAL_f, &b_depositedElecEnergyECAL_f);
   fChain->SetBranchAddress("depositedElecEnergyECAL_r", depositedElecEnergyECAL_r, &b_depositedElecEnergyECAL_r);
   fChain->SetBranchAddress("depositedElecEnergyHCALAct", &depositedElecEnergyHCALAct, &b_depositedElecEnergyHCALAct);
   fChain->SetBranchAddress("depositedElecEnergyHCALPas", &depositedElecEnergyHCALPas, &b_depositedElecEnergyHCALPas);
   fChain->SetBranchAddress("depositedElecEnergyWorld", &depositedElecEnergyWorld, &b_depositedElecEnergyWorld);
   fChain->SetBranchAddress("depositedElecEnergyTimingGap", &depositedElecEnergyTimingGap, &b_depositedElecEnergyTimingGap);
   fChain->SetBranchAddress("depositedElecEnergyServices", &depositedElecEnergyServices, &b_depositedElecEnergyServices);
   fChain->SetBranchAddress("depositedElecEnergyEcalGap", &depositedElecEnergyEcalGap, &b_depositedElecEnergyEcalGap);
   fChain->SetBranchAddress("depositedElecEnergyEcalDet", &depositedElecEnergyEcalDet, &b_depositedElecEnergyEcalDet);
   fChain->SetBranchAddress("depositedElecEnergySolenoid", &depositedElecEnergySolenoid, &b_depositedElecEnergySolenoid);
   fChain->SetBranchAddress("Edep_Tracker_layer", Edep_Tracker_layer, &b_Edep_Tracker_layer);
   fChain->SetBranchAddress("Edep_Timing_f_ch", Edep_Timing_f_ch, &b_Edep_Timing_f_ch);
   fChain->SetBranchAddress("Edep_Timing_r_ch", Edep_Timing_r_ch, &b_Edep_Timing_r_ch);
   fChain->SetBranchAddress("Edep_ECAL_f_ch", Edep_ECAL_f_ch, &b_Edep_ECAL_f_ch);
   fChain->SetBranchAddress("Edep_ECAL_r_ch", Edep_ECAL_r_ch, &b_Edep_ECAL_r_ch);
   fChain->SetBranchAddress("IonEdep_ECAL_f_ch", IonEdep_ECAL_f_ch, &b_IonEdep_ECAL_f_ch);
   fChain->SetBranchAddress("IonEdep_ECAL_r_ch", IonEdep_ECAL_r_ch, &b_IonEdep_ECAL_r_ch);
   fChain->SetBranchAddress("tot_phot_cer_Timing_f_total", &tot_phot_cer_Timing_f_total, &b_tot_phot_cer_Timing_f_total);
   fChain->SetBranchAddress("tot_phot_cer_Timing_r_total", &tot_phot_cer_Timing_r_total, &b_tot_phot_cer_Timing_r_total);
   fChain->SetBranchAddress("ECAL_f_total_S", &ECAL_f_total_S, &b_ECAL_f_total_S);
   fChain->SetBranchAddress("ECAL_f_total_C", &ECAL_f_total_C, &b_ECAL_f_total_C);
   fChain->SetBranchAddress("ECAL_r_total_S", &ECAL_r_total_S, &b_ECAL_r_total_S);
   fChain->SetBranchAddress("ECAL_r_total_C", &ECAL_r_total_C, &b_ECAL_r_total_C);
   fChain->SetBranchAddress("tot_phot_cer_ECAL_scinti_f_particleID", tot_phot_cer_ECAL_scinti_f_particleID, &b_tot_phot_cer_ECAL_scinti_f_particleID);
   fChain->SetBranchAddress("tot_phot_cer_ECAL_scinti_r_particleID", tot_phot_cer_ECAL_scinti_r_particleID, &b_tot_phot_cer_ECAL_scinti_r_particleID);
   fChain->SetBranchAddress("tot_phot_cer_ECAL_cheren_f_particleID", tot_phot_cer_ECAL_cheren_f_particleID, &b_tot_phot_cer_ECAL_cheren_f_particleID);
   fChain->SetBranchAddress("tot_phot_cer_ECAL_cheren_r_particleID", tot_phot_cer_ECAL_cheren_r_particleID, &b_tot_phot_cer_ECAL_cheren_r_particleID);
   fChain->SetBranchAddress("tot_phot_cer_HCAL", &tot_phot_cer_HCAL, &b_tot_phot_cer_HCAL);
   fChain->SetBranchAddress("SDdetected_ff_S", &SDdetected_ff_S, &b_SDdetected_ff_S);
   fChain->SetBranchAddress("SDdetected_ff_C", &SDdetected_ff_C, &b_SDdetected_ff_C);
   fChain->SetBranchAddress("tot_phot_cer_SDdetected_fr", &tot_phot_cer_SDdetected_fr, &b_tot_phot_cer_SDdetected_fr);
   fChain->SetBranchAddress("tot_phot_cer_SDdetected_rf", &tot_phot_cer_SDdetected_rf, &b_tot_phot_cer_SDdetected_rf);
   fChain->SetBranchAddress("SDdetected_rr_S", &SDdetected_rr_S, &b_SDdetected_rr_S);
   fChain->SetBranchAddress("SDdetected_rr_C", &SDdetected_rr_C, &b_SDdetected_rr_C);
   fChain->SetBranchAddress("tot_phot_cer_should_det1", &tot_phot_cer_should_det1, &b_tot_phot_cer_should_det1);
   fChain->SetBranchAddress("tot_phot_cer_should_det2", &tot_phot_cer_should_det2, &b_tot_phot_cer_should_det2);
   fChain->SetBranchAddress("tot_phot_cer_should_det3", &tot_phot_cer_should_det3, &b_tot_phot_cer_should_det3);
   fChain->SetBranchAddress("tot_phot_cer_should_det4", &tot_phot_cer_should_det4, &b_tot_phot_cer_should_det4);
   fChain->SetBranchAddress("E_Zdep_0to5000mm_total", E_Zdep_0to5000mm_total, &b_E_Zdep_0to5000mm_total);
   fChain->SetBranchAddress("E_Zdep_0to5000mm_Pion_n", E_Zdep_0to5000mm_Pion_n, &b_E_Zdep_0to5000mm_Pion_n);
   fChain->SetBranchAddress("E_Zdep_0to5000mm_Positron", E_Zdep_0to5000mm_Positron, &b_E_Zdep_0to5000mm_Positron);
   fChain->SetBranchAddress("E_Zdep_0to5000mm_Electron", E_Zdep_0to5000mm_Electron, &b_E_Zdep_0to5000mm_Electron);
   fChain->SetBranchAddress("E_Zdep_0to5000mm_Photon", E_Zdep_0to5000mm_Photon, &b_E_Zdep_0to5000mm_Photon);
   fChain->SetBranchAddress("E_Zdep_0to5000mm_Pion_p", E_Zdep_0to5000mm_Pion_p, &b_E_Zdep_0to5000mm_Pion_p);
   fChain->SetBranchAddress("E_Zdep_0to5000mm_Kion", E_Zdep_0to5000mm_Kion, &b_E_Zdep_0to5000mm_Kion);
   fChain->SetBranchAddress("E_Zdep_0to5000mm_Neutron", E_Zdep_0to5000mm_Neutron, &b_E_Zdep_0to5000mm_Neutron);
   fChain->SetBranchAddress("E_Zdep_0to5000mm_Proton", E_Zdep_0to5000mm_Proton, &b_E_Zdep_0to5000mm_Proton);
   fChain->SetBranchAddress("E_Tdep_0to5ns_total", E_Tdep_0to5ns_total, &b_E_Tdep_0to5ns_total);
   fChain->SetBranchAddress("E_Tdep_0to5ns_Pion_n", E_Tdep_0to5ns_Pion_n, &b_E_Tdep_0to5ns_Pion_n);
   fChain->SetBranchAddress("E_Tdep_0to5ns_Positron", E_Tdep_0to5ns_Positron, &b_E_Tdep_0to5ns_Positron);
   fChain->SetBranchAddress("E_Tdep_0to5ns_Electron", E_Tdep_0to5ns_Electron, &b_E_Tdep_0to5ns_Electron);
   fChain->SetBranchAddress("E_Tdep_0to5ns_Photon", E_Tdep_0to5ns_Photon, &b_E_Tdep_0to5ns_Photon);
   fChain->SetBranchAddress("E_Tdep_0to5ns_Pion_p", E_Tdep_0to5ns_Pion_p, &b_E_Tdep_0to5ns_Pion_p);
   fChain->SetBranchAddress("E_Tdep_0to5ns_Kion", E_Tdep_0to5ns_Kion, &b_E_Tdep_0to5ns_Kion);
   fChain->SetBranchAddress("E_Tdep_0to5ns_Neutron", E_Tdep_0to5ns_Neutron, &b_E_Tdep_0to5ns_Neutron);
   fChain->SetBranchAddress("E_Tdep_0to5ns_Proton", E_Tdep_0to5ns_Proton, &b_E_Tdep_0to5ns_Proton);
   fChain->SetBranchAddress("EnergyECAL_r_enter", &EnergyECAL_r_enter, &b_EnergyECAL_r_enter);
   fChain->SetBranchAddress("EnergyECAL_r_leave", &EnergyECAL_r_leave, &b_EnergyECAL_r_leave);
   fChain->SetBranchAddress("EnergyECAL_r_enterK", &EnergyECAL_r_enterK, &b_EnergyECAL_r_enterK);
   fChain->SetBranchAddress("EnergyECAL_r_leaveK", &EnergyECAL_r_leaveK, &b_EnergyECAL_r_leaveK);
   fChain->SetBranchAddress("pdgid_escape", &pdgid_escape, &b_pdgid_escape);
   fChain->SetBranchAddress("energy_escape", &energy_escape, &b_energy_escape);
   fChain->SetBranchAddress("positionx_escape", &positionx_escape, &b_positionx_escape);
   fChain->SetBranchAddress("positiony_escape", &positiony_escape, &b_positiony_escape);
   fChain->SetBranchAddress("positionz_escape", &positionz_escape, &b_positionz_escape);
   Notify();
}

Bool_t missing::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void missing::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t missing::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef missing_cxx
