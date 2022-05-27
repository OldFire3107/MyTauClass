//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 11 16:49:26 2022 by ROOT version 6.22/09
// from TTree tree/tree
// found on file: flat_tuple_1.root
//////////////////////////////////////////////////////////

#ifndef MyTauClass_h
#define MyTauClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TH1.h>
#include <TMath.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <map>
#include <string>

using namespace std;

class MyTauClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxHLT_BPH_isFired = 1;

   // Declaration of leaf types
   Int_t           genParticle_N;
   vector<float>   *genParticle_pt;
   vector<float>   *genParticle_eta;
   vector<float>   *genParticle_phi;
   vector<float>   *genParticle_mass;
   vector<int>     *genParticle_pdgId;
   vector<int>     *genParticle_status;
   vector<int>     *genParticle_isPrompt;
   vector<int>     *genParticle_isDirectPromptTauDecayProduct;
   vector<int>     *genParticle_isDirectHardProcessTauDecayProductFinalState;
   vector<int>     *genParticle_fromHardProcessFinalState;
   vector<vector<int> > *genParticle_mother;
   vector<vector<float> > *genParticle_mother_pt;
   vector<int>     *genParticle_nMoth;
   vector<int>     *genParticle_nDau;
   vector<vector<int> > *genParticle_dau;
   vector<vector<int> > *genParticle_pdgs;
   vector<vector<int> > *genParticle_layers;
   vector<vector<float> > *genParticle_ppt;
   vector<vector<float> > *genParticle_peta;
   vector<vector<float> > *genParticle_pphi;
   vector<vector<int> > *genParticle_isfinal;
   Float_t         lheV_pt;
   Float_t         lheHT;
   Int_t           lheNj;
   Int_t           lheNb;
   Int_t           lheNl;
   Float_t         lheV_mass;
   Float_t         genWeight;
   Float_t         genFacWeightUp;
   Float_t         genFacWeightDown;
   Float_t         genRenWeightUp;
   Float_t         genRenWeightDown;
   Float_t         genFacRenWeightUp;
   Float_t         genFacRenWeightDown;
   Float_t         qScale;
   Float_t         PDF_rms;
   vector<float>   *PDF_x;
   vector<float>   *PDF_xPDF;
   vector<int>     *PDF_id;
   Int_t           EVENT_event;
   Int_t           EVENT_run;
   Int_t           EVENT_lumiBlock;
   vector<float>   *nPuVtxTrue;
   vector<int>     *nPuVtx;
   vector<int>     *bX;
   Int_t           PV_N;
   Bool_t          PV_filter;
   vector<float>   *PV_chi2;
   vector<float>   *PV_ndof;
   vector<float>   *PV_rho;
   vector<float>   *PV_z;
   vector<float>   *BeamSpot_x0;
   vector<float>   *BeamSpot_y0;
   vector<float>   *BeamSpot_z0;
   Int_t           HLT_BPH_isFired_;
   string          HLT_BPH_isFired_first[kMaxHLT_BPH_isFired];
   Bool_t          HLT_BPH_isFired_second[kMaxHLT_BPH_isFired];   //[HLT_BPH_isFired_]
   Int_t           JpsiTau_nCandidates;
   Bool_t          JpsiTau_isJpsiMu;
   Bool_t          JpsiTau_isJpsiTau2Mu;
   vector<double>  *truth_tau_dipion1_mass;
   vector<double>  *truth_tau_dipion1_pt;
   vector<double>  *truth_tau_dipion1_eta;
   vector<double>  *truth_tau_dipion1_phi;
   vector<double>  *truth_tau_dipion2_mass;
   vector<double>  *truth_tau_dipion2_pt;
   vector<double>  *truth_tau_dipion2_eta;
   vector<double>  *truth_tau_dipion2_phi;
   Float_t         JpsiTau_mu1_pt;
   Float_t         JpsiTau_mu1_eta;
   Float_t         JpsiTau_mu1_phi;
   Float_t         JpsiTau_mu1_mass;
   Int_t           JpsiTau_mu1_q;
   Int_t           JpsiTau_mu1_isLoose;
   Int_t           JpsiTau_mu1_isTight;
   Int_t           JpsiTau_mu1_isPF;
   Int_t           JpsiTau_mu1_isGlobal;
   Int_t           JpsiTau_mu1_isTracker;
   Int_t           JpsiTau_mu1_isSoft;
   Float_t         JpsiTau_mu1_vx;
   Float_t         JpsiTau_mu1_vy;
   Float_t         JpsiTau_mu1_vz;
   Float_t         JpsiTau_mu1_dbiso;
   Float_t         JpsiTau_mu2_pt;
   Float_t         JpsiTau_mu2_eta;
   Float_t         JpsiTau_mu2_phi;
   Float_t         JpsiTau_mu2_mass;
   Int_t           JpsiTau_mu2_q;
   Int_t           JpsiTau_mu2_isLoose;
   Int_t           JpsiTau_mu2_isTight;
   Int_t           JpsiTau_mu2_isPF;
   Int_t           JpsiTau_mu2_isGlobal;
   Int_t           JpsiTau_mu2_isTracker;
   Int_t           JpsiTau_mu2_isSoft;
   Float_t         JpsiTau_mu2_vx;
   Float_t         JpsiTau_mu2_vy;
   Float_t         JpsiTau_mu2_vz;
   Float_t         JpsiTau_mu2_dbiso;
   vector<float>   *JpsiTau_tau_pt;
   vector<float>   *JpsiTau_tau_eta;
   vector<float>   *JpsiTau_tau_phi;
   vector<float>   *JpsiTau_tau_mass;
   vector<int>     *JpsiTau_tau_q;
   vector<float>   *JpsiTau_tau_vx;
   vector<float>   *JpsiTau_tau_vy;
   vector<float>   *JpsiTau_tau_vz;
   vector<float>   *JpsiTau_tau_max_dr_3prong;
   vector<float>   *JpsiTau_tau_lip;
   vector<float>   *JpsiTau_tau_lips;
   vector<float>   *JpsiTau_tau_pvip;
   vector<float>   *JpsiTau_tau_pvips;
   vector<float>   *JpsiTau_tau_fl3d;
   vector<float>   *JpsiTau_tau_fls3d;
   vector<float>   *JpsiTau_tau_alpha;
   vector<float>   *JpsiTau_tau_vprob;
   vector<float>   *JpsiTau_tau_fl3d_wjpsi;
   vector<float>   *JpsiTau_tau_fls3d_wjpsi;
   vector<float>   *JpsiTau_tau_sumofdnn;
   vector<float>   *JpsiTau_tau_sumofdnn_1prong;
   vector<float>   *JpsiTau_tau_sumofdnn_otherB;
   vector<float>   *JpsiTau_tau_sumofdnn_pu;
   vector<float>   *JpsiTau_tau_rhomass1;
   vector<float>   *JpsiTau_tau_rhomass2;
   vector<float>   *JpsiTau_tau_rhopt1;
   vector<float>   *JpsiTau_tau_rhopt2;
   vector<float>   *JpsiTau_tau_rhomass_ss;
   vector<float>   *JpsiTau_tau_pi1_pt;
   vector<float>   *JpsiTau_tau_pi1_eta;
   vector<float>   *JpsiTau_tau_pi1_phi;
   vector<float>   *JpsiTau_tau_pi1_mass;
   vector<int>     *JpsiTau_tau_pi1_q;
   vector<float>   *JpsiTau_tau_pi1_doca3d;
   vector<float>   *JpsiTau_tau_pi1_doca3de;
   vector<float>   *JpsiTau_tau_pi1_doca2d;
   vector<float>   *JpsiTau_tau_pi1_doca2de;
   vector<float>   *JpsiTau_tau_pi1_dz;
   vector<float>   *JpsiTau_tau_pi1_near_dz;
   vector<bool>    *JpsiTau_tau_pi1_isAssociate;
   vector<int>     *JpsiTau_tau_pi1_pvAssociationQuality;
   vector<bool>    *JpsiTau_tau_pi1_isBdecay;
   vector<int>     *JpsiTau_tau_pi1_isBdecaypdg;
   vector<int>     *JpsiTau_tau_pi1_isBdecayppdg;
   vector<bool>    *JpsiTau_tau_pi1_isSignal;
   vector<int>     *JpsiTau_tau_pi1_nprong;
   vector<int>     *JpsiTau_tau_pi1_nprong_pi0;
   vector<float>   *JpsiTau_tau_pi1_dnn;
   vector<float>   *JpsiTau_tau_pi1_dnn_1prong;
   vector<float>   *JpsiTau_tau_pi1_dnn_otherB;
   vector<float>   *JpsiTau_tau_pi1_dnn_pu;
   vector<bool>    *JpsiTau_tau_pi1_trigMatch;
   vector<float>   *JpsiTau_tau_pi1_trigMatch_dr;
   vector<float>   *JpsiTau_tau_pi2_pt;
   vector<float>   *JpsiTau_tau_pi2_eta;
   vector<float>   *JpsiTau_tau_pi2_phi;
   vector<float>   *JpsiTau_tau_pi2_mass;
   vector<int>     *JpsiTau_tau_pi2_q;
   vector<float>   *JpsiTau_tau_pi2_doca3d;
   vector<float>   *JpsiTau_tau_pi2_doca3de;
   vector<float>   *JpsiTau_tau_pi2_doca2d;
   vector<float>   *JpsiTau_tau_pi2_doca2de;
   vector<float>   *JpsiTau_tau_pi2_dz;
   vector<float>   *JpsiTau_tau_pi2_near_dz;
   vector<bool>    *JpsiTau_tau_pi2_isAssociate;
   vector<int>     *JpsiTau_tau_pi2_pvAssociationQuality;
   vector<bool>    *JpsiTau_tau_pi2_isBdecay;
   vector<int>     *JpsiTau_tau_pi2_isBdecaypdg;
   vector<int>     *JpsiTau_tau_pi2_isBdecayppdg;
   vector<bool>    *JpsiTau_tau_pi2_isSignal;
   vector<int>     *JpsiTau_tau_pi2_nprong;
   vector<int>     *JpsiTau_tau_pi2_nprong_pi0;
   vector<float>   *JpsiTau_tau_pi2_dnn;
   vector<float>   *JpsiTau_tau_pi2_dnn_1prong;
   vector<float>   *JpsiTau_tau_pi2_dnn_otherB;
   vector<float>   *JpsiTau_tau_pi2_dnn_pu;
   vector<bool>    *JpsiTau_tau_pi2_trigMatch;
   vector<float>   *JpsiTau_tau_pi2_trigMatch_dr;
   vector<float>   *JpsiTau_tau_pi3_pt;
   vector<float>   *JpsiTau_tau_pi3_eta;
   vector<float>   *JpsiTau_tau_pi3_phi;
   vector<float>   *JpsiTau_tau_pi3_mass;
   vector<int>     *JpsiTau_tau_pi3_q;
   vector<float>   *JpsiTau_tau_pi3_doca3d;
   vector<float>   *JpsiTau_tau_pi3_doca3de;
   vector<float>   *JpsiTau_tau_pi3_doca2d;
   vector<float>   *JpsiTau_tau_pi3_doca2de;
   vector<float>   *JpsiTau_tau_pi3_dz;
   vector<float>   *JpsiTau_tau_pi3_near_dz;
   vector<bool>    *JpsiTau_tau_pi3_isAssociate;
   vector<int>     *JpsiTau_tau_pi3_pvAssociationQuality;
   vector<bool>    *JpsiTau_tau_pi3_isBdecay;
   vector<int>     *JpsiTau_tau_pi3_isBdecaypdg;
   vector<int>     *JpsiTau_tau_pi3_isBdecayppdg;
   vector<bool>    *JpsiTau_tau_pi3_isSignal;
   vector<int>     *JpsiTau_tau_pi3_nprong;
   vector<int>     *JpsiTau_tau_pi3_nprong_pi0;
   vector<float>   *JpsiTau_tau_pi3_dnn;
   vector<float>   *JpsiTau_tau_pi3_dnn_1prong;
   vector<float>   *JpsiTau_tau_pi3_dnn_otherB;
   vector<float>   *JpsiTau_tau_pi3_dnn_pu;
   vector<bool>    *JpsiTau_tau_pi3_trigMatch;
   vector<float>   *JpsiTau_tau_pi3_trigMatch_dr;
   vector<float>   *JpsiTau_tau_delta_chi2;
   vector<int>     *JpsiTau_tau_delta_n_ch;
   vector<int>     *JpsiTau_tau_delta_n_mu;
   vector<float>   *JpsiTau_tau_vweight;
   vector<float>   *JpsiTau_tau_refit_vx;
   vector<float>   *JpsiTau_tau_refit_vy;
   vector<float>   *JpsiTau_tau_refit_vz;
   vector<float>   *JpsiTau_tau_refit_chi2;
   vector<float>   *JpsiTau_tau_refit_ndof;
   vector<float>   *JpsiTau_tau_refit_rho;
   vector<vector<float> > *JpsiTau_tau_iso;
   vector<vector<int> > *JpsiTau_tau_iso_ntracks;
   vector<vector<float> > *JpsiTau_tau_iso_mindoca;
   vector<float>   *JpsiTau_ptbal;
   vector<float>   *JpsiTau_jpsi_tau_alpha;
   Float_t         JpsiTau_PV_vx;
   Float_t         JpsiTau_PV_vy;
   Float_t         JpsiTau_PV_vz;
   Float_t         JpsiTau_bbPV_vx;
   Float_t         JpsiTau_bbPV_vy;
   Float_t         JpsiTau_bbPV_vz;
   Float_t         JpsiTau_bbPV_chi2;
   Float_t         JpsiTau_bbPV_ndof;
   Float_t         JpsiTau_bbPV_rho;
   Float_t         JpsiTau_Jpsi_pt;
   Float_t         JpsiTau_Jpsi_eta;
   Float_t         JpsiTau_Jpsi_phi;
   Float_t         JpsiTau_Jpsi_mass;
   Float_t         JpsiTau_Jpsi_vprob;
   Float_t         JpsiTau_Jpsi_lip;
   Float_t         JpsiTau_Jpsi_lips;
   Float_t         JpsiTau_Jpsi_pvip;
   Float_t         JpsiTau_Jpsi_pvips;
   Float_t         JpsiTau_Jpsi_fl3d;
   Float_t         JpsiTau_Jpsi_fls3d;
   Float_t         JpsiTau_Jpsi_alpha;
   Float_t         JpsiTau_Jpsi_maxdoca;
   Float_t         JpsiTau_Jpsi_mindoca;
   Float_t         JpsiTau_Jpsi_vx;
   Float_t         JpsiTau_Jpsi_vy;
   Float_t         JpsiTau_Jpsi_vz;
   vector<float>   *JpsiTau_B_pt;
   vector<float>   *JpsiTau_B_eta;
   vector<float>   *JpsiTau_B_phi;
   vector<float>   *JpsiTau_B_mass;
   vector<float>   *JpsiTau_B_mcorr;
   vector<float>   *JpsiTau_B_vprob;
   vector<float>   *JpsiTau_B_lip;
   vector<float>   *JpsiTau_B_lips;
   vector<float>   *JpsiTau_B_pvip;
   vector<float>   *JpsiTau_B_pvips;
   vector<float>   *JpsiTau_B_fl3d;
   vector<float>   *JpsiTau_B_fls3d;
   vector<float>   *JpsiTau_B_alpha;
   vector<float>   *JpsiTau_B_maxdoca;
   vector<float>   *JpsiTau_B_mindoca;
   vector<float>   *JpsiTau_B_vx;
   vector<float>   *JpsiTau_B_vy;
   vector<float>   *JpsiTau_B_vz;
   vector<float>   *JpsiTau_B_q2;
   vector<float>   *JpsiTau_B_mm2;
   vector<float>   *JpsiTau_B_ptmiss;
   vector<float>   *JpsiTau_B_Es;
   vector<float>   *JpsiTau_B_ptback;
   vector<float>   *JpsiTau_B_pt_simple;
   vector<float>   *JpsiTau_B_eta_simple;
   vector<float>   *JpsiTau_B_phi_simple;
   vector<float>   *JpsiTau_B_mass_simple;
   vector<float>   *JpsiTau_B_q2_simple;
   vector<float>   *JpsiTau_B_mm2_simple;
   vector<float>   *JpsiTau_B_ptmiss_simple;
   vector<float>   *JpsiTau_B_Es_simple;
   vector<float>   *JpsiTau_B_ptback_simple;
   Float_t         genWeightBkgB;
   Float_t         JpsiTau_genPV_vx;
   Float_t         JpsiTau_genPV_vy;
   Float_t         JpsiTau_genPV_vz;
   Float_t         JpsiTau_genSV_vx;
   Float_t         JpsiTau_genSV_vy;
   Float_t         JpsiTau_genSV_vz;
   Int_t           JpsiTau_ngenmuons;
   Int_t           JpsiTau_isgenmatched;
   Float_t         JpsiTau_q2_gen;
   Int_t           JpsiTau_nBc;
   Float_t         JpsiTau_B_pt_gen;
   Float_t         JpsiTau_B_eta_gen;
   Float_t         JpsiTau_B_phi_gen;
   Float_t         JpsiTau_B_mass_gen;
   vector<float>   *JpsiTau_hammer_ebe;
   vector<float>   *JpsiTau_hammer_ebe_e0_up;
   vector<float>   *JpsiTau_hammer_ebe_e0_down;
   vector<float>   *JpsiTau_hammer_ebe_e1_up;
   vector<float>   *JpsiTau_hammer_ebe_e1_down;
   vector<float>   *JpsiTau_hammer_ebe_e2_up;
   vector<float>   *JpsiTau_hammer_ebe_e2_down;
   vector<float>   *JpsiTau_hammer_ebe_e3_up;
   vector<float>   *JpsiTau_hammer_ebe_e3_down;
   vector<float>   *JpsiTau_hammer_ebe_e4_up;
   vector<float>   *JpsiTau_hammer_ebe_e4_down;
   vector<float>   *JpsiTau_hammer_ebe_e5_up;
   vector<float>   *JpsiTau_hammer_ebe_e5_down;
   vector<float>   *JpsiTau_hammer_ebe_e6_up;
   vector<float>   *JpsiTau_hammer_ebe_e6_down;
   vector<float>   *JpsiTau_hammer_ebe_e7_up;
   vector<float>   *JpsiTau_hammer_ebe_e7_down;
   vector<float>   *JpsiTau_hammer_ebe_e8_up;
   vector<float>   *JpsiTau_hammer_ebe_e8_down;
   vector<float>   *JpsiTau_hammer_ebe_e9_up;
   vector<float>   *JpsiTau_hammer_ebe_e9_down;
   vector<float>   *JpsiTau_hammer_ebe_e10_up;
   vector<float>   *JpsiTau_hammer_ebe_e10_down;
   vector<float>   *JpsiTau_hammer_ebe_e11_up;
   vector<float>   *JpsiTau_hammer_ebe_e11_down;
   vector<float>   *JpsiTau_hammer_ebe_e12_up;
   vector<float>   *JpsiTau_hammer_ebe_e12_down;
   vector<float>   *JpsiTau_hammer_ebe_e13_up;
   vector<float>   *JpsiTau_hammer_ebe_e13_down;
   vector<float>   *JpsiTau_hammer_ebe_e14_up;
   vector<float>   *JpsiTau_hammer_ebe_e14_down;
   vector<float>   *JpsiTau_hammer_ebe_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e0_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e0_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e1_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e1_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e2_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e2_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e3_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e3_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e4_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e4_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e5_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e5_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e6_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e6_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e7_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e7_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e8_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e8_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e9_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e9_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e10_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e10_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e11_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e11_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e12_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e12_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e13_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e13_down_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e14_up_lattice;
   vector<float>   *JpsiTau_hammer_ebe_e14_down_lattice;
   Int_t           JpsiTau_nch;
   Int_t           JpsiTau_nch_before;
   vector<float>   *JpsiTau_gen_pion_pt;
   vector<float>   *JpsiTau_gen_pion_eta;
   vector<float>   *JpsiTau_gen_pion_phi;
   vector<bool>    *JpsiTau_gen_pion_matched;
   vector<float>   *JpsiTau_gen_tau_pt;
   vector<float>   *JpsiTau_gen_tau_eta;
   vector<float>   *JpsiTau_gen_tau_phi;
   vector<int>     *JpsiTau_gen_tau_nprong;
   vector<int>     *JpsiTau_gen_tau_nmatched;
   Int_t           JpsiTau_st_nch;
   Int_t           JpsiTau_st_nch_matched;
   Int_t           JpsiTau_st_n_charged_pions;
   Int_t           JpsiTau_st_n_neutral_pions;
   Int_t           JpsiTau_st_n_mu_decay;
   Int_t           JpsiTau_st_n_e_decay;
   Int_t           JpsiTau_st_n_occurance;
   Int_t           JpsiTau_st_decayid;
   Float_t         JpsiTau_st_gentau_pt;
   Float_t         JpsiTau_st_gentau_eta;
   Float_t         JpsiTau_st_gentau_phi;
   Float_t         JpsiTau_st_genjpsi_pt;
   Float_t         JpsiTau_st_genjpsi_eta;
   Float_t         JpsiTau_st_genjpsi_phi;
   Int_t           JpsiTau_general_ntau;
   vector<bool>    *JpsiTau_st_isBdecay;
   vector<int>     *JpsiTau_st_isBdecaypdg;
   vector<int>     *JpsiTau_st_isBdecayppdg;
   vector<bool>    *JpsiTau_st_isSignal;
   vector<int>     *JpsiTau_st_nprong;
   vector<int>     *JpsiTau_st_nprong_pi0;
   vector<int>     *JpsiTau_st_idx;
   vector<float>   *JpsiTau_st_doca3d;
   vector<float>   *JpsiTau_st_doca2d;
   vector<float>   *JpsiTau_st_doca3ds;
   vector<float>   *JpsiTau_st_doca2ds;
   vector<float>   *JpsiTau_st_doca3de;
   vector<float>   *JpsiTau_st_doca2de;
   vector<float>   *JpsiTau_st_dz;
   vector<bool>    *JpsiTau_st_isAssociate;
   vector<float>   *JpsiTau_st_near_dz;
   vector<float>   *JpsiTau_st_dr_jpsi;
   vector<bool>    *JpsiTau_st_trigMatch;
   vector<float>   *JpsiTau_st_trigMatch_dr;
   vector<int>     *JpsiTau_st_pvAssociationQuality;
   vector<float>   *JpsiTau_st_pt;
   vector<float>   *JpsiTau_st_eta;
   vector<float>   *JpsiTau_st_phi;
   vector<int>     *JpsiTau_st_charge;
   vector<float>   *JpsiTau_st_mass;
   vector<float>   *JpsiTau_st_dnn;
   vector<float>   *JpsiTau_st_dnn_1prong;
   vector<float>   *JpsiTau_st_dnn_otherB;
   vector<float>   *JpsiTau_st_dnn_pu;
   vector<int>     *JpsiTau_st_matchidx;
   Float_t         JpsiTau_perEVT_mc;
   Float_t         JpsiTau_perEVT_data;

   // List of branches
   TBranch        *b_genParticle_N;   //!
   TBranch        *b_genParticle_pt;   //!
   TBranch        *b_genParticle_eta;   //!
   TBranch        *b_genParticle_phi;   //!
   TBranch        *b_genParticle_mass;   //!
   TBranch        *b_genParticle_pdgId;   //!
   TBranch        *b_genParticle_status;   //!
   TBranch        *b_genParticle_isPrompt;   //!
   TBranch        *b_genParticle_isDirectPromptTauDecayProduct;   //!
   TBranch        *b_genParticle_isDirectHardProcessTauDecayProductFinalState;   //!
   TBranch        *b_genParticle_fromHardProcessFinalState;   //!
   TBranch        *b_genParticle_mother;   //!
   TBranch        *b_genParticle_mother_pt;   //!
   TBranch        *b_genParticle_nMoth;   //!
   TBranch        *b_genParticle_nDau;   //!
   TBranch        *b_genParticle_dau;   //!
   TBranch        *b_genParticle_pdgs;   //!
   TBranch        *b_genParticle_layers;   //!
   TBranch        *b_genParticle_ppt;   //!
   TBranch        *b_genParticle_peta;   //!
   TBranch        *b_genParticle_pphi;   //!
   TBranch        *b_genParticle_isfinal;   //!
   TBranch        *b_lheV_pt;   //!
   TBranch        *b_lheHT;   //!
   TBranch        *b_lheNj;   //!
   TBranch        *b_lheNb;   //!
   TBranch        *b_lheNl;   //!
   TBranch        *b_lheV_mass;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genFacWeightUp;   //!
   TBranch        *b_genFacWeightDown;   //!
   TBranch        *b_genRenWeightUp;   //!
   TBranch        *b_genRenWeightDown;   //!
   TBranch        *b_genFacRenWeightUp;   //!
   TBranch        *b_genFacRenWeightDown;   //!
   TBranch        *b_qScale;   //!
   TBranch        *b_PDF_rms;   //!
   TBranch        *b_PDF_x;   //!
   TBranch        *b_PDF_xPDF;   //!
   TBranch        *b_PDF_id;   //!
   TBranch        *b_EVENT_event;   //!
   TBranch        *b_EVENT_run;   //!
   TBranch        *b_EVENT_lumiBlock;   //!
   TBranch        *b_nPuVtxTrue;   //!
   TBranch        *b_nPuVtx;   //!
   TBranch        *b_bX;   //!
   TBranch        *b_PV_N;   //!
   TBranch        *b_PV_filter;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_rho;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_BeamSpot_x0;   //!
   TBranch        *b_BeamSpot_y0;   //!
   TBranch        *b_BeamSpot_z0;   //!
   TBranch        *b_HLT_BPH_isFired_;   //!
   TBranch        *b_HLT_BPH_isFired_first;   //!
   TBranch        *b_HLT_BPH_isFired_second;   //!
   TBranch        *b_JpsiTau_nCandidates;   //!
   TBranch        *b_JpsiTau_isJpsiMu;   //!
   TBranch        *b_JpsiTau_isJpsiTau2Mu;   //!
   TBranch        *b_truth_tau_dipion1_mass;   //!
   TBranch        *b_truth_tau_dipion1_pt;   //!
   TBranch        *b_truth_tau_dipion1_eta;   //!
   TBranch        *b_truth_tau_dipion1_phi;   //!
   TBranch        *b_truth_tau_dipion2_mass;   //!
   TBranch        *b_truth_tau_dipion2_pt;   //!
   TBranch        *b_truth_tau_dipion2_eta;   //!
   TBranch        *b_truth_tau_dipion2_phi;   //!
   TBranch        *b_JpsiTau_mu1_pt;   //!
   TBranch        *b_JpsiTau_mu1_eta;   //!
   TBranch        *b_JpsiTau_mu1_phi;   //!
   TBranch        *b_JpsiTau_mu1_mass;   //!
   TBranch        *b_JpsiTau_mu1_q;   //!
   TBranch        *b_JpsiTau_mu1_isLoose;   //!
   TBranch        *b_JpsiTau_mu1_isTight;   //!
   TBranch        *b_JpsiTau_mu1_isPF;   //!
   TBranch        *b_JpsiTau_mu1_isGlobal;   //!
   TBranch        *b_JpsiTau_mu1_isTracker;   //!
   TBranch        *b_JpsiTau_mu1_isSoft;   //!
   TBranch        *b_JpsiTau_mu1_vx;   //!
   TBranch        *b_JpsiTau_mu1_vy;   //!
   TBranch        *b_JpsiTau_mu1_vz;   //!
   TBranch        *b_JpsiTau_mu1_dbiso;   //!
   TBranch        *b_JpsiTau_mu2_pt;   //!
   TBranch        *b_JpsiTau_mu2_eta;   //!
   TBranch        *b_JpsiTau_mu2_phi;   //!
   TBranch        *b_JpsiTau_mu2_mass;   //!
   TBranch        *b_JpsiTau_mu2_q;   //!
   TBranch        *b_JpsiTau_mu2_isLoose;   //!
   TBranch        *b_JpsiTau_mu2_isTight;   //!
   TBranch        *b_JpsiTau_mu2_isPF;   //!
   TBranch        *b_JpsiTau_mu2_isGlobal;   //!
   TBranch        *b_JpsiTau_mu2_isTracker;   //!
   TBranch        *b_JpsiTau_mu2_isSoft;   //!
   TBranch        *b_JpsiTau_mu2_vx;   //!
   TBranch        *b_JpsiTau_mu2_vy;   //!
   TBranch        *b_JpsiTau_mu2_vz;   //!
   TBranch        *b_JpsiTau_mu2_dbiso;   //!
   TBranch        *b_JpsiTau_tau_pt;   //!
   TBranch        *b_JpsiTau_tau_eta;   //!
   TBranch        *b_JpsiTau_tau_phi;   //!
   TBranch        *b_JpsiTau_tau_mass;   //!
   TBranch        *b_JpsiTau_tau_q;   //!
   TBranch        *b_JpsiTau_tau_vx;   //!
   TBranch        *b_JpsiTau_tau_vy;   //!
   TBranch        *b_JpsiTau_tau_vz;   //!
   TBranch        *b_JpsiTau_tau_max_dr_3prong;   //!
   TBranch        *b_JpsiTau_tau_lip;   //!
   TBranch        *b_JpsiTau_tau_lips;   //!
   TBranch        *b_JpsiTau_tau_pvip;   //!
   TBranch        *b_JpsiTau_tau_pvips;   //!
   TBranch        *b_JpsiTau_tau_fl3d;   //!
   TBranch        *b_JpsiTau_tau_fls3d;   //!
   TBranch        *b_JpsiTau_tau_alpha;   //!
   TBranch        *b_JpsiTau_tau_vprob;   //!
   TBranch        *b_JpsiTau_tau_fl3d_wjpsi;   //!
   TBranch        *b_JpsiTau_tau_fls3d_wjpsi;   //!
   TBranch        *b_JpsiTau_tau_sumofdnn;   //!
   TBranch        *b_JpsiTau_tau_sumofdnn_1prong;   //!
   TBranch        *b_JpsiTau_tau_sumofdnn_otherB;   //!
   TBranch        *b_JpsiTau_tau_sumofdnn_pu;   //!
   TBranch        *b_JpsiTau_tau_rhomass1;   //!
   TBranch        *b_JpsiTau_tau_rhomass2;   //!
   TBranch        *b_JpsiTau_tau_rhopt1;   //!
   TBranch        *b_JpsiTau_tau_rhopt2;   //!
   TBranch        *b_JpsiTau_tau_rhomass_ss;   //!
   TBranch        *b_JpsiTau_tau_pi1_pt;   //!
   TBranch        *b_JpsiTau_tau_pi1_eta;   //!
   TBranch        *b_JpsiTau_tau_pi1_phi;   //!
   TBranch        *b_JpsiTau_tau_pi1_mass;   //!
   TBranch        *b_JpsiTau_tau_pi1_q;   //!
   TBranch        *b_JpsiTau_tau_pi1_doca3d;   //!
   TBranch        *b_JpsiTau_tau_pi1_doca3de;   //!
   TBranch        *b_JpsiTau_tau_pi1_doca2d;   //!
   TBranch        *b_JpsiTau_tau_pi1_doca2de;   //!
   TBranch        *b_JpsiTau_tau_pi1_dz;   //!
   TBranch        *b_JpsiTau_tau_pi1_near_dz;   //!
   TBranch        *b_JpsiTau_tau_pi1_isAssociate;   //!
   TBranch        *b_JpsiTau_tau_pi1_pvAssociationQuality;   //!
   TBranch        *b_JpsiTau_tau_pi1_isBdecay;   //!
   TBranch        *b_JpsiTau_tau_pi1_isBdecaypdg;   //!
   TBranch        *b_JpsiTau_tau_pi1_isBdecayppdg;   //!
   TBranch        *b_JpsiTau_tau_pi1_isSignal;   //!
   TBranch        *b_JpsiTau_tau_pi1_nprong;   //!
   TBranch        *b_JpsiTau_tau_pi1_nprong_pi0;   //!
   TBranch        *b_JpsiTau_tau_pi1_dnn;   //!
   TBranch        *b_JpsiTau_tau_pi1_dnn_1prong;   //!
   TBranch        *b_JpsiTau_tau_pi1_dnn_otherB;   //!
   TBranch        *b_JpsiTau_tau_pi1_dnn_pu;   //!
   TBranch        *b_JpsiTau_tau_pi1_trigMatch;   //!
   TBranch        *b_JpsiTau_tau_pi1_trigMatch_dr;   //!
   TBranch        *b_JpsiTau_tau_pi2_pt;   //!
   TBranch        *b_JpsiTau_tau_pi2_eta;   //!
   TBranch        *b_JpsiTau_tau_pi2_phi;   //!
   TBranch        *b_JpsiTau_tau_pi2_mass;   //!
   TBranch        *b_JpsiTau_tau_pi2_q;   //!
   TBranch        *b_JpsiTau_tau_pi2_doca3d;   //!
   TBranch        *b_JpsiTau_tau_pi2_doca3de;   //!
   TBranch        *b_JpsiTau_tau_pi2_doca2d;   //!
   TBranch        *b_JpsiTau_tau_pi2_doca2de;   //!
   TBranch        *b_JpsiTau_tau_pi2_dz;   //!
   TBranch        *b_JpsiTau_tau_pi2_near_dz;   //!
   TBranch        *b_JpsiTau_tau_pi2_isAssociate;   //!
   TBranch        *b_JpsiTau_tau_pi2_pvAssociationQuality;   //!
   TBranch        *b_JpsiTau_tau_pi2_isBdecay;   //!
   TBranch        *b_JpsiTau_tau_pi2_isBdecaypdg;   //!
   TBranch        *b_JpsiTau_tau_pi2_isBdecayppdg;   //!
   TBranch        *b_JpsiTau_tau_pi2_isSignal;   //!
   TBranch        *b_JpsiTau_tau_pi2_nprong;   //!
   TBranch        *b_JpsiTau_tau_pi2_nprong_pi0;   //!
   TBranch        *b_JpsiTau_tau_pi2_dnn;   //!
   TBranch        *b_JpsiTau_tau_pi2_dnn_1prong;   //!
   TBranch        *b_JpsiTau_tau_pi2_dnn_otherB;   //!
   TBranch        *b_JpsiTau_tau_pi2_dnn_pu;   //!
   TBranch        *b_JpsiTau_tau_pi2_trigMatch;   //!
   TBranch        *b_JpsiTau_tau_pi2_trigMatch_dr;   //!
   TBranch        *b_JpsiTau_tau_pi3_pt;   //!
   TBranch        *b_JpsiTau_tau_pi3_eta;   //!
   TBranch        *b_JpsiTau_tau_pi3_phi;   //!
   TBranch        *b_JpsiTau_tau_pi3_mass;   //!
   TBranch        *b_JpsiTau_tau_pi3_q;   //!
   TBranch        *b_JpsiTau_tau_pi3_doca3d;   //!
   TBranch        *b_JpsiTau_tau_pi3_doca3de;   //!
   TBranch        *b_JpsiTau_tau_pi3_doca2d;   //!
   TBranch        *b_JpsiTau_tau_pi3_doca2de;   //!
   TBranch        *b_JpsiTau_tau_pi3_dz;   //!
   TBranch        *b_JpsiTau_tau_pi3_near_dz;   //!
   TBranch        *b_JpsiTau_tau_pi3_isAssociate;   //!
   TBranch        *b_JpsiTau_tau_pi3_pvAssociationQuality;   //!
   TBranch        *b_JpsiTau_tau_pi3_isBdecay;   //!
   TBranch        *b_JpsiTau_tau_pi3_isBdecaypdg;   //!
   TBranch        *b_JpsiTau_tau_pi3_isBdecayppdg;   //!
   TBranch        *b_JpsiTau_tau_pi3_isSignal;   //!
   TBranch        *b_JpsiTau_tau_pi3_nprong;   //!
   TBranch        *b_JpsiTau_tau_pi3_nprong_pi0;   //!
   TBranch        *b_JpsiTau_tau_pi3_dnn;   //!
   TBranch        *b_JpsiTau_tau_pi3_dnn_1prong;   //!
   TBranch        *b_JpsiTau_tau_pi3_dnn_otherB;   //!
   TBranch        *b_JpsiTau_tau_pi3_dnn_pu;   //!
   TBranch        *b_JpsiTau_tau_pi3_trigMatch;   //!
   TBranch        *b_JpsiTau_tau_pi3_trigMatch_dr;   //!
   TBranch        *b_JpsiTau_tau_delta_chi2;   //!
   TBranch        *b_JpsiTau_tau_delta_n_ch;   //!
   TBranch        *b_JpsiTau_tau_delta_n_mu;   //!
   TBranch        *b_JpsiTau_tau_vweight;   //!
   TBranch        *b_JpsiTau_tau_refit_vx;   //!
   TBranch        *b_JpsiTau_tau_refit_vy;   //!
   TBranch        *b_JpsiTau_tau_refit_vz;   //!
   TBranch        *b_JpsiTau_tau_refit_chi2;   //!
   TBranch        *b_JpsiTau_tau_refit_ndof;   //!
   TBranch        *b_JpsiTau_tau_refit_rho;   //!
   TBranch        *b_JpsiTau_tau_iso;   //!
   TBranch        *b_JpsiTau_tau_iso_ntracks;   //!
   TBranch        *b_JpsiTau_tau_iso_mindoca;   //!
   TBranch        *b_JpsiTau_ptbal;   //!
   TBranch        *b_JpsiTau_jpsi_tau_alpha;   //!
   TBranch        *b_JpsiTau_PV_vx;   //!
   TBranch        *b_JpsiTau_PV_vy;   //!
   TBranch        *b_JpsiTau_PV_vz;   //!
   TBranch        *b_JpsiTau_bbPV_vx;   //!
   TBranch        *b_JpsiTau_bbPV_vy;   //!
   TBranch        *b_JpsiTau_bbPV_vz;   //!
   TBranch        *b_JpsiTau_bbPV_chi2;   //!
   TBranch        *b_JpsiTau_bbPV_ndof;   //!
   TBranch        *b_JpsiTau_bbPV_rho;   //!
   TBranch        *b_JpsiTau_Jpsi_pt;   //!
   TBranch        *b_JpsiTau_Jpsi_eta;   //!
   TBranch        *b_JpsiTau_Jpsi_phi;   //!
   TBranch        *b_JpsiTau_Jpsi_mass;   //!
   TBranch        *b_JpsiTau_Jpsi_vprob;   //!
   TBranch        *b_JpsiTau_Jpsi_lip;   //!
   TBranch        *b_JpsiTau_Jpsi_lips;   //!
   TBranch        *b_JpsiTau_Jpsi_pvip;   //!
   TBranch        *b_JpsiTau_Jpsi_pvips;   //!
   TBranch        *b_JpsiTau_Jpsi_fl3d;   //!
   TBranch        *b_JpsiTau_Jpsi_fls3d;   //!
   TBranch        *b_JpsiTau_Jpsi_alpha;   //!
   TBranch        *b_JpsiTau_Jpsi_maxdoca;   //!
   TBranch        *b_JpsiTau_Jpsi_mindoca;   //!
   TBranch        *b_JpsiTau_Jpsi_vx;   //!
   TBranch        *b_JpsiTau_Jpsi_vy;   //!
   TBranch        *b_JpsiTau_Jpsi_vz;   //!
   TBranch        *b_JpsiTau_B_pt;   //!
   TBranch        *b_JpsiTau_B_eta;   //!
   TBranch        *b_JpsiTau_B_phi;   //!
   TBranch        *b_JpsiTau_B_mass;   //!
   TBranch        *b_JpsiTau_B_mcorr;   //!
   TBranch        *b_JpsiTau_B_vprob;   //!
   TBranch        *b_JpsiTau_B_lip;   //!
   TBranch        *b_JpsiTau_B_lips;   //!
   TBranch        *b_JpsiTau_B_pvip;   //!
   TBranch        *b_JpsiTau_B_pvips;   //!
   TBranch        *b_JpsiTau_B_fl3d;   //!
   TBranch        *b_JpsiTau_B_fls3d;   //!
   TBranch        *b_JpsiTau_B_alpha;   //!
   TBranch        *b_JpsiTau_B_maxdoca;   //!
   TBranch        *b_JpsiTau_B_mindoca;   //!
   TBranch        *b_JpsiTau_B_vx;   //!
   TBranch        *b_JpsiTau_B_vy;   //!
   TBranch        *b_JpsiTau_B_vz;   //!
   TBranch        *b_JpsiTau_B_q2;   //!
   TBranch        *b_JpsiTau_B_mm2;   //!
   TBranch        *b_JpsiTau_B_ptmiss;   //!
   TBranch        *b_JpsiTau_B_Es;   //!
   TBranch        *b_JpsiTau_B_ptback;   //!
   TBranch        *b_JpsiTau_B_pt_simple;   //!
   TBranch        *b_JpsiTau_B_eta_simple;   //!
   TBranch        *b_JpsiTau_B_phi_simple;   //!
   TBranch        *b_JpsiTau_B_mass_simple;   //!
   TBranch        *b_JpsiTau_B_q2_simple;   //!
   TBranch        *b_JpsiTau_B_mm2_simple;   //!
   TBranch        *b_JpsiTau_B_ptmiss_simple;   //!
   TBranch        *b_JpsiTau_B_Es_simple;   //!
   TBranch        *b_JpsiTau_B_ptback_simple;   //!
   TBranch        *b_genWeightBkgB;   //!
   TBranch        *b_JpsiTau_genPV_vx;   //!
   TBranch        *b_JpsiTau_genPV_vy;   //!
   TBranch        *b_JpsiTau_genPV_vz;   //!
   TBranch        *b_JpsiTau_genSV_vx;   //!
   TBranch        *b_JpsiTau_genSV_vy;   //!
   TBranch        *b_JpsiTau_genSV_vz;   //!
   TBranch        *b_JpsiTau_ngenmuons;   //!
   TBranch        *b_JpsiTau_isgenmatched;   //!
   TBranch        *b_JpsiTau_q2_gen;   //!
   TBranch        *b_JpsiTau_nBc;   //!
   TBranch        *b_JpsiTau_B_pt_gen;   //!
   TBranch        *b_JpsiTau_B_eta_gen;   //!
   TBranch        *b_JpsiTau_B_phi_gen;   //!
   TBranch        *b_JpsiTau_B_mass_gen;   //!
   TBranch        *b_JpsiTau_hammer_ebe;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e0_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e0_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e1_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e1_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e2_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e2_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e3_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e3_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e4_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e4_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e5_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e5_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e6_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e6_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e7_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e7_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e8_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e8_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e9_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e9_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e10_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e10_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e11_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e11_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e12_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e12_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e13_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e13_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e14_up;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e14_down;   //!
   TBranch        *b_JpsiTau_hammer_ebe_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e0_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e0_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e1_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e1_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e2_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e2_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e3_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e3_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e4_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e4_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e5_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e5_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e6_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e6_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e7_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e7_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e8_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e8_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e9_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e9_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e10_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e10_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e11_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e11_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e12_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e12_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e13_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e13_down_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e14_up_lattice;   //!
   TBranch        *b_JpsiTau_hammer_ebe_e14_down_lattice;   //!
   TBranch        *b_JpsiTau_nch;   //!
   TBranch        *b_JpsiTau_nch_before;   //!
   TBranch        *b_JpsiTau_gen_pion_pt;   //!
   TBranch        *b_JpsiTau_gen_pion_eta;   //!
   TBranch        *b_JpsiTau_gen_pion_phi;   //!
   TBranch        *b_JpsiTau_gen_pion_matched;   //!
   TBranch        *b_JpsiTau_gen_tau_pt;   //!
   TBranch        *b_JpsiTau_gen_tau_eta;   //!
   TBranch        *b_JpsiTau_gen_tau_phi;   //!
   TBranch        *b_JpsiTau_gen_tau_nprong;   //!
   TBranch        *b_JpsiTau_gen_tau_nmatched;   //!
   TBranch        *b_JpsiTau_st_nch;   //!
   TBranch        *b_JpsiTau_st_nch_matched;   //!
   TBranch        *b_JpsiTau_st_n_charged_pions;   //!
   TBranch        *b_JpsiTau_st_n_neutral_pions;   //!
   TBranch        *b_JpsiTau_st_n_mu_decay;   //!
   TBranch        *b_JpsiTau_st_n_e_decay;   //!
   TBranch        *b_JpsiTau_st_n_occurance;   //!
   TBranch        *b_JpsiTau_st_decayid;   //!
   TBranch        *b_JpsiTau_st_gentau_pt;   //!
   TBranch        *b_JpsiTau_st_gentau_eta;   //!
   TBranch        *b_JpsiTau_st_gentau_phi;   //!
   TBranch        *b_JpsiTau_st_genjpsi_pt;   //!
   TBranch        *b_JpsiTau_st_genjpsi_eta;   //!
   TBranch        *b_JpsiTau_st_genjpsi_phi;   //!
   TBranch        *b_JpsiTau_general_ntau;   //!
   TBranch        *b_JpsiTau_st_isBdecay;   //!
   TBranch        *b_JpsiTau_st_isBdecaypdg;   //!
   TBranch        *b_JpsiTau_st_isBdecayppdg;   //!
   TBranch        *b_JpsiTau_st_isSignal;   //!
   TBranch        *b_JpsiTau_st_nprong;   //!
   TBranch        *b_JpsiTau_st_nprong_pi0;   //!
   TBranch        *b_JpsiTau_st_idx;   //!
   TBranch        *b_JpsiTau_st_doca3d;   //!
   TBranch        *b_JpsiTau_st_doca2d;   //!
   TBranch        *b_JpsiTau_st_doca3ds;   //!
   TBranch        *b_JpsiTau_st_doca2ds;   //!
   TBranch        *b_JpsiTau_st_doca3de;   //!
   TBranch        *b_JpsiTau_st_doca2de;   //!
   TBranch        *b_JpsiTau_st_dz;   //!
   TBranch        *b_JpsiTau_st_isAssociate;   //!
   TBranch        *b_JpsiTau_st_near_dz;   //!
   TBranch        *b_JpsiTau_st_dr_jpsi;   //!
   TBranch        *b_JpsiTau_st_trigMatch;   //!
   TBranch        *b_JpsiTau_st_trigMatch_dr;   //!
   TBranch        *b_JpsiTau_st_pvAssociationQuality;   //!
   TBranch        *b_JpsiTau_st_pt;   //!
   TBranch        *b_JpsiTau_st_eta;   //!
   TBranch        *b_JpsiTau_st_phi;   //!
   TBranch        *b_JpsiTau_st_charge;   //!
   TBranch        *b_JpsiTau_st_mass;   //!
   TBranch        *b_JpsiTau_st_dnn;   //!
   TBranch        *b_JpsiTau_st_dnn_1prong;   //!
   TBranch        *b_JpsiTau_st_dnn_otherB;   //!
   TBranch        *b_JpsiTau_st_dnn_pu;   //!
   TBranch        *b_JpsiTau_st_matchidx;   //!
   TBranch        *b_JpsiTau_perEVT_mc;   //!
   TBranch        *b_JpsiTau_perEVT_data;   //!

   TString        FileNameIn;
   TFile          *FileIn=NULL;
   TString        FileNameOut;
   TFile          *FileOut=NULL;
   Float_t        pi_pt[3];
   Float_t        pi_eta[3];
   Float_t        pi_phi[3];
   Float_t        pir_pt[3];
   Float_t        pir_eta[3];
   Float_t        pir_phi[3];
   Int_t          pir_q[3];
   Float_t        pir_dce3d[3];
   Float_t        pir_dz[3];
   Float_t        tau_true_pt;
   Bool_t         pi_flag;
   TTree          *tree1 = NULL;

   // Book here my histograms
   TH1F *h1pttau=NULL;
   TH1F *h1etatau=NULL;
   TH1F *h1phitau=NULL;
   TH1F *h1nprongtau=NULL;
   TH1F *h1nmatchedtau=NULL;

   TH1F *h1ptpi=NULL;
   TH1F *h1etapi=NULL;
   TH1F *h1phipi=NULL;
   

   TH1F *h1pt1=NULL;
   TH1F *h1pt2=NULL;
   TH1F *h1pt3=NULL;

   TH1F *h1ratio=NULL;
   TH1F *h15pt=NULL;
   TH1F *h110pt=NULL;
   TH1F *h120pt=NULL;
   TH1F *h1expt=NULL;

   // PART 2
   TH2F *h2etapt=NULL;
   TH2F *h2phipt=NULL;
   TH2F *h2etaphi=NULL;
   TH2F *h2dRpt=NULL;

   MyTauClass();
   // MyTauClass(TTree *tree);
   MyTauClass(const char *file, TString fname);
   virtual ~MyTauClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     InitOut();
   virtual void     Loop();
   virtual void     show(const char* file="final.root");
   virtual void     MegaLoop(const char *str);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     SwapValue(Float_t &a, Float_t &b);
   virtual void     SwapValue(Int_t &a, Int_t &b);
   virtual void     CombAdd(vector<vector<Int_t>>& comb, Int_t arr[3], int n);
   virtual float    deltaPhi(float phi1, float phi2);
   virtual Float_t  getDeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2);
};

#endif

#ifdef MyTauClass_cxx
MyTauClass::MyTauClass() : fChain(0)
{
   // TTree *tree = NULL;
   // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("flat_tuple_1.root");
   // if (!f || !f->IsOpen()) {
   //    f = new TFile("flat_tuple_1.root");
   // }
   // TDirectory * dir = (TDirectory*)f->Get("flat_tuple_1.root:/ntuplizer");
   // dir->GetObject("tree",tree);
}

// MyTauClass::MyTauClass(TTree *tree) : fChain(0) 
// {
// // if parameter tree is not specified (or zero), connect the file
// // used to generate this class and read the Tree.
//    Init(tree);
// }

MyTauClass::MyTauClass(const char *str, TString fname) : fChain(0) 
{
// Overloads the previous function to open a tree
   const char* fileo="files/histos_";
   TString fileout = fileo;
   fileout += fname;
   fileout += ".root";
   FileNameOut = fileout;

   if(FileOut) delete FileOut;
   FileOut = new TFile(fileout, "RECREATE");
   tree1=NULL;
   InitOut();

   if(FileIn) delete FileIn;
   FileNameIn = str;
   FileIn = new TFile(str);
   TTree *tree=NULL;

   TDirectory *dir = (TDirectory*)FileIn->Get("ntuplizer");
   dir->GetObject("tree",tree);

   Init(tree);
}

MyTauClass::~MyTauClass()
{
   delete h1pttau;
   delete h1etatau;
   delete h1phitau;
   delete h1nprongtau;
   delete h1nmatchedtau;

   delete h1ptpi;
   delete h1etapi;
   delete h1phipi;
   

   delete h1pt1;
   delete h1pt2;
   delete h1pt3;

   delete h1ratio;
   delete h15pt;
   delete h110pt;
   delete h120pt;
   delete h1expt;

   delete tree1;
   delete FileIn;
   delete FileOut;

   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyTauClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyTauClass::LoadTree(Long64_t entry)
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

void MyTauClass::InitOut()
{

// The tree
   if (tree1) delete tree1;

   tree1 = new TTree("triplet","triplet");
   tree1->Branch("pi1_pt", &pi_pt[0], "pi1_pt/F");
   tree1->Branch("pi2_pt", &pi_pt[1], "pi2_pt/F");
   tree1->Branch("pi3_pt", &pi_pt[2], "pi3_pt/F");
   tree1->Branch("pi1_eta", &pi_eta[0], "pi1_eta/F");
   tree1->Branch("pi2_eta", &pi_eta[1], "pi2_eta/F");
   tree1->Branch("pi3_eta", &pi_eta[2], "pi3_eta/F");
   tree1->Branch("pi1_phi", &pi_phi[0], "pi1_phi/F");
   tree1->Branch("pi2_phi", &pi_phi[1], "pi2_phi/F");
   tree1->Branch("pi3_phi", &pi_phi[2], "pi3_phi/F");
   tree1->Branch("tau_pt", &tau_true_pt, "tau_pt/F");
   // TTree *tree1 = new TTree("reco_triplet","reco_triplet");
   tree1->Branch("pi1r_pt", &pir_pt[0], "pi1r_pt/F");
   tree1->Branch("pi2r_pt", &pir_pt[1], "pi2r_pt/F");
   tree1->Branch("pi3r_pt", &pir_pt[2], "pi3r_pt/F");
   tree1->Branch("pi1r_eta", &pir_eta[0], "pi1r_eta/F");
   tree1->Branch("pi2r_eta", &pir_eta[1], "pi2r_eta/F");
   tree1->Branch("pi3r_eta", &pir_eta[2], "pi3r_eta/F");
   tree1->Branch("pi1r_phi", &pir_phi[0], "pi1r_phi/F");
   tree1->Branch("pi2r_phi", &pir_phi[1], "pi2r_phi/F");
   tree1->Branch("pi3r_phi", &pir_phi[2], "pi3r_phi/F");
   tree1->Branch("pi1r_q", &pir_q[0], "pi1r_q/I");
   tree1->Branch("pi2r_q", &pir_q[1], "pi2r_q/I");
   tree1->Branch("pi3r_q", &pir_q[2], "pi3r_q/I");
   tree1->Branch("pi1r_dce3d", &pir_dce3d[0], "pi1r_dce3d/F");
   tree1->Branch("pi2r_dce3d", &pir_dce3d[1], "pi2r_dce3d/F");
   tree1->Branch("pi3r_dce3d", &pir_dce3d[2], "pi3r_dce3d/F");
   tree1->Branch("pi1r_dz", &pir_dz[0], "pi1r_dz/F");
   tree1->Branch("pi2r_dz", &pir_dz[1], "pi2r_dz/F");
   tree1->Branch("pi3r_dz", &pir_dz[2], "pi3r_dz/F");
   tree1->Branch("flag", &pi_flag, "flag/O");

   if (h1pttau)
   {
      delete h1pttau;
      delete h1etatau;
      delete h1phitau;
      delete h1nprongtau;
      delete h1nmatchedtau;

      delete h1ptpi;
      delete h1etapi;
      delete h1phipi;
      

      delete h1pt1;
      delete h1pt2;
      delete h1pt3;

      delete h1ratio;
      delete h15pt;
      delete h110pt;
      delete h120pt;
      delete h1expt;
   }

   // Book here my histograms
   h1pttau = new TH1F("h1pttau","pt of the gen tau",50,1.,40.);
   h1etatau = new TH1F("h1etatau","eta of the gen tau",25,-2.5,2.5);
   h1phitau = new TH1F("h1phitau","phi of the gen tau",25,-1*TMath::Pi(),TMath::Pi());
   h1nprongtau = new TH1F("h1nprongtau","nprong of the gen tau",4,0,4);
   h1nmatchedtau = new TH1F("h1nmatchedtau","nmatched of the gen tau",4,0,4);

   h1ptpi = new TH1F("h1ptpi","pt of the gen pion",50,0.,20.);
   h1etapi = new TH1F("h1etapi","eta of the gen pion",25,-2.5,2.5);
   h1phipi = new TH1F("h1phipi","phi of the gen pion",25,-1*TMath::Pi(),TMath::Pi());
   

   h1pt1 = new TH1F("h1pt1","pt of the gen pion1",50,0.,10.);
   h1pt2 = new TH1F("h1pt2","pt of the gen pion2",50,0.,10.);
   h1pt3 = new TH1F("h1pt3","pt of the gen pion3",50,0.,10.);

   h1ratio = new TH1F("h1ratio","pt of pion1/pt of tau",50,0,1);
   h15pt = new TH1F("h15pt", "dR for pt upto 5 GeV", 25, 0, 2);
   h110pt = new TH1F("h110pt", "dR for pt in 5-10 GeV range", 25, 0, 2);
   h120pt = new TH1F("h120pt", "dR for pt in 10-20 GeV range", 25, 0, 1);
   h1expt = new TH1F("h1expt", "dR for pt above 20 GeV", 25, 0, 1);

   // PART 2
   h2etapt = new TH2F("h2etapt", "eta vs pt", 50, 1, 40, 25, -2.5, 2.5);
   h2phipt = new TH2F("h2phipt", "phi vs pt", 50, 1, 40, 25, -1*TMath::Pi(), TMath::Pi());
   h2etaphi = new TH2F("h2etaphi", "eta vs phi", 25, -1*TMath::Pi(), TMath::Pi(), 25, -2.5, 2.5);
   h2dRpt = new TH2F("h2dRpt", "dR vs pt", 50, 1, 40, 25, 0, 1.8);
}

void MyTauClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genParticle_pt = 0;
   genParticle_eta = 0;
   genParticle_phi = 0;
   genParticle_mass = 0;
   genParticle_pdgId = 0;
   genParticle_status = 0;
   genParticle_isPrompt = 0;
   genParticle_isDirectPromptTauDecayProduct = 0;
   genParticle_isDirectHardProcessTauDecayProductFinalState = 0;
   genParticle_fromHardProcessFinalState = 0;
   genParticle_mother = 0;
   genParticle_mother_pt = 0;
   genParticle_nMoth = 0;
   genParticle_nDau = 0;
   genParticle_dau = 0;
   genParticle_pdgs = 0;
   genParticle_layers = 0;
   genParticle_ppt = 0;
   genParticle_peta = 0;
   genParticle_pphi = 0;
   genParticle_isfinal = 0;
   PDF_x = 0;
   PDF_xPDF = 0;
   PDF_id = 0;
   nPuVtxTrue = 0;
   nPuVtx = 0;
   bX = 0;
   PV_chi2 = 0;
   PV_ndof = 0;
   PV_rho = 0;
   PV_z = 0;
   BeamSpot_x0 = 0;
   BeamSpot_y0 = 0;
   BeamSpot_z0 = 0;
   truth_tau_dipion1_mass = 0;
   truth_tau_dipion1_pt = 0;
   truth_tau_dipion1_eta = 0;
   truth_tau_dipion1_phi = 0;
   truth_tau_dipion2_mass = 0;
   truth_tau_dipion2_pt = 0;
   truth_tau_dipion2_eta = 0;
   truth_tau_dipion2_phi = 0;
   JpsiTau_tau_pt = 0;
   JpsiTau_tau_eta = 0;
   JpsiTau_tau_phi = 0;
   JpsiTau_tau_mass = 0;
   JpsiTau_tau_q = 0;
   JpsiTau_tau_vx = 0;
   JpsiTau_tau_vy = 0;
   JpsiTau_tau_vz = 0;
   JpsiTau_tau_max_dr_3prong = 0;
   JpsiTau_tau_lip = 0;
   JpsiTau_tau_lips = 0;
   JpsiTau_tau_pvip = 0;
   JpsiTau_tau_pvips = 0;
   JpsiTau_tau_fl3d = 0;
   JpsiTau_tau_fls3d = 0;
   JpsiTau_tau_alpha = 0;
   JpsiTau_tau_vprob = 0;
   JpsiTau_tau_fl3d_wjpsi = 0;
   JpsiTau_tau_fls3d_wjpsi = 0;
   JpsiTau_tau_sumofdnn = 0;
   JpsiTau_tau_sumofdnn_1prong = 0;
   JpsiTau_tau_sumofdnn_otherB = 0;
   JpsiTau_tau_sumofdnn_pu = 0;
   JpsiTau_tau_rhomass1 = 0;
   JpsiTau_tau_rhomass2 = 0;
   JpsiTau_tau_rhopt1 = 0;
   JpsiTau_tau_rhopt2 = 0;
   JpsiTau_tau_rhomass_ss = 0;
   JpsiTau_tau_pi1_pt = 0;
   JpsiTau_tau_pi1_eta = 0;
   JpsiTau_tau_pi1_phi = 0;
   JpsiTau_tau_pi1_mass = 0;
   JpsiTau_tau_pi1_q = 0;
   JpsiTau_tau_pi1_doca3d = 0;
   JpsiTau_tau_pi1_doca3de = 0;
   JpsiTau_tau_pi1_doca2d = 0;
   JpsiTau_tau_pi1_doca2de = 0;
   JpsiTau_tau_pi1_dz = 0;
   JpsiTau_tau_pi1_near_dz = 0;
   JpsiTau_tau_pi1_isAssociate = 0;
   JpsiTau_tau_pi1_pvAssociationQuality = 0;
   JpsiTau_tau_pi1_isBdecay = 0;
   JpsiTau_tau_pi1_isBdecaypdg = 0;
   JpsiTau_tau_pi1_isBdecayppdg = 0;
   JpsiTau_tau_pi1_isSignal = 0;
   JpsiTau_tau_pi1_nprong = 0;
   JpsiTau_tau_pi1_nprong_pi0 = 0;
   JpsiTau_tau_pi1_dnn = 0;
   JpsiTau_tau_pi1_dnn_1prong = 0;
   JpsiTau_tau_pi1_dnn_otherB = 0;
   JpsiTau_tau_pi1_dnn_pu = 0;
   JpsiTau_tau_pi1_trigMatch = 0;
   JpsiTau_tau_pi1_trigMatch_dr = 0;
   JpsiTau_tau_pi2_pt = 0;
   JpsiTau_tau_pi2_eta = 0;
   JpsiTau_tau_pi2_phi = 0;
   JpsiTau_tau_pi2_mass = 0;
   JpsiTau_tau_pi2_q = 0;
   JpsiTau_tau_pi2_doca3d = 0;
   JpsiTau_tau_pi2_doca3de = 0;
   JpsiTau_tau_pi2_doca2d = 0;
   JpsiTau_tau_pi2_doca2de = 0;
   JpsiTau_tau_pi2_dz = 0;
   JpsiTau_tau_pi2_near_dz = 0;
   JpsiTau_tau_pi2_isAssociate = 0;
   JpsiTau_tau_pi2_pvAssociationQuality = 0;
   JpsiTau_tau_pi2_isBdecay = 0;
   JpsiTau_tau_pi2_isBdecaypdg = 0;
   JpsiTau_tau_pi2_isBdecayppdg = 0;
   JpsiTau_tau_pi2_isSignal = 0;
   JpsiTau_tau_pi2_nprong = 0;
   JpsiTau_tau_pi2_nprong_pi0 = 0;
   JpsiTau_tau_pi2_dnn = 0;
   JpsiTau_tau_pi2_dnn_1prong = 0;
   JpsiTau_tau_pi2_dnn_otherB = 0;
   JpsiTau_tau_pi2_dnn_pu = 0;
   JpsiTau_tau_pi2_trigMatch = 0;
   JpsiTau_tau_pi2_trigMatch_dr = 0;
   JpsiTau_tau_pi3_pt = 0;
   JpsiTau_tau_pi3_eta = 0;
   JpsiTau_tau_pi3_phi = 0;
   JpsiTau_tau_pi3_mass = 0;
   JpsiTau_tau_pi3_q = 0;
   JpsiTau_tau_pi3_doca3d = 0;
   JpsiTau_tau_pi3_doca3de = 0;
   JpsiTau_tau_pi3_doca2d = 0;
   JpsiTau_tau_pi3_doca2de = 0;
   JpsiTau_tau_pi3_dz = 0;
   JpsiTau_tau_pi3_near_dz = 0;
   JpsiTau_tau_pi3_isAssociate = 0;
   JpsiTau_tau_pi3_pvAssociationQuality = 0;
   JpsiTau_tau_pi3_isBdecay = 0;
   JpsiTau_tau_pi3_isBdecaypdg = 0;
   JpsiTau_tau_pi3_isBdecayppdg = 0;
   JpsiTau_tau_pi3_isSignal = 0;
   JpsiTau_tau_pi3_nprong = 0;
   JpsiTau_tau_pi3_nprong_pi0 = 0;
   JpsiTau_tau_pi3_dnn = 0;
   JpsiTau_tau_pi3_dnn_1prong = 0;
   JpsiTau_tau_pi3_dnn_otherB = 0;
   JpsiTau_tau_pi3_dnn_pu = 0;
   JpsiTau_tau_pi3_trigMatch = 0;
   JpsiTau_tau_pi3_trigMatch_dr = 0;
   JpsiTau_tau_delta_chi2 = 0;
   JpsiTau_tau_delta_n_ch = 0;
   JpsiTau_tau_delta_n_mu = 0;
   JpsiTau_tau_vweight = 0;
   JpsiTau_tau_refit_vx = 0;
   JpsiTau_tau_refit_vy = 0;
   JpsiTau_tau_refit_vz = 0;
   JpsiTau_tau_refit_chi2 = 0;
   JpsiTau_tau_refit_ndof = 0;
   JpsiTau_tau_refit_rho = 0;
   JpsiTau_tau_iso = 0;
   JpsiTau_tau_iso_ntracks = 0;
   JpsiTau_tau_iso_mindoca = 0;
   JpsiTau_ptbal = 0;
   JpsiTau_jpsi_tau_alpha = 0;
   JpsiTau_B_pt = 0;
   JpsiTau_B_eta = 0;
   JpsiTau_B_phi = 0;
   JpsiTau_B_mass = 0;
   JpsiTau_B_mcorr = 0;
   JpsiTau_B_vprob = 0;
   JpsiTau_B_lip = 0;
   JpsiTau_B_lips = 0;
   JpsiTau_B_pvip = 0;
   JpsiTau_B_pvips = 0;
   JpsiTau_B_fl3d = 0;
   JpsiTau_B_fls3d = 0;
   JpsiTau_B_alpha = 0;
   JpsiTau_B_maxdoca = 0;
   JpsiTau_B_mindoca = 0;
   JpsiTau_B_vx = 0;
   JpsiTau_B_vy = 0;
   JpsiTau_B_vz = 0;
   JpsiTau_B_q2 = 0;
   JpsiTau_B_mm2 = 0;
   JpsiTau_B_ptmiss = 0;
   JpsiTau_B_Es = 0;
   JpsiTau_B_ptback = 0;
   JpsiTau_B_pt_simple = 0;
   JpsiTau_B_eta_simple = 0;
   JpsiTau_B_phi_simple = 0;
   JpsiTau_B_mass_simple = 0;
   JpsiTau_B_q2_simple = 0;
   JpsiTau_B_mm2_simple = 0;
   JpsiTau_B_ptmiss_simple = 0;
   JpsiTau_B_Es_simple = 0;
   JpsiTau_B_ptback_simple = 0;
   JpsiTau_hammer_ebe = 0;
   JpsiTau_hammer_ebe_e0_up = 0;
   JpsiTau_hammer_ebe_e0_down = 0;
   JpsiTau_hammer_ebe_e1_up = 0;
   JpsiTau_hammer_ebe_e1_down = 0;
   JpsiTau_hammer_ebe_e2_up = 0;
   JpsiTau_hammer_ebe_e2_down = 0;
   JpsiTau_hammer_ebe_e3_up = 0;
   JpsiTau_hammer_ebe_e3_down = 0;
   JpsiTau_hammer_ebe_e4_up = 0;
   JpsiTau_hammer_ebe_e4_down = 0;
   JpsiTau_hammer_ebe_e5_up = 0;
   JpsiTau_hammer_ebe_e5_down = 0;
   JpsiTau_hammer_ebe_e6_up = 0;
   JpsiTau_hammer_ebe_e6_down = 0;
   JpsiTau_hammer_ebe_e7_up = 0;
   JpsiTau_hammer_ebe_e7_down = 0;
   JpsiTau_hammer_ebe_e8_up = 0;
   JpsiTau_hammer_ebe_e8_down = 0;
   JpsiTau_hammer_ebe_e9_up = 0;
   JpsiTau_hammer_ebe_e9_down = 0;
   JpsiTau_hammer_ebe_e10_up = 0;
   JpsiTau_hammer_ebe_e10_down = 0;
   JpsiTau_hammer_ebe_e11_up = 0;
   JpsiTau_hammer_ebe_e11_down = 0;
   JpsiTau_hammer_ebe_e12_up = 0;
   JpsiTau_hammer_ebe_e12_down = 0;
   JpsiTau_hammer_ebe_e13_up = 0;
   JpsiTau_hammer_ebe_e13_down = 0;
   JpsiTau_hammer_ebe_e14_up = 0;
   JpsiTau_hammer_ebe_e14_down = 0;
   JpsiTau_hammer_ebe_lattice = 0;
   JpsiTau_hammer_ebe_e0_up_lattice = 0;
   JpsiTau_hammer_ebe_e0_down_lattice = 0;
   JpsiTau_hammer_ebe_e1_up_lattice = 0;
   JpsiTau_hammer_ebe_e1_down_lattice = 0;
   JpsiTau_hammer_ebe_e2_up_lattice = 0;
   JpsiTau_hammer_ebe_e2_down_lattice = 0;
   JpsiTau_hammer_ebe_e3_up_lattice = 0;
   JpsiTau_hammer_ebe_e3_down_lattice = 0;
   JpsiTau_hammer_ebe_e4_up_lattice = 0;
   JpsiTau_hammer_ebe_e4_down_lattice = 0;
   JpsiTau_hammer_ebe_e5_up_lattice = 0;
   JpsiTau_hammer_ebe_e5_down_lattice = 0;
   JpsiTau_hammer_ebe_e6_up_lattice = 0;
   JpsiTau_hammer_ebe_e6_down_lattice = 0;
   JpsiTau_hammer_ebe_e7_up_lattice = 0;
   JpsiTau_hammer_ebe_e7_down_lattice = 0;
   JpsiTau_hammer_ebe_e8_up_lattice = 0;
   JpsiTau_hammer_ebe_e8_down_lattice = 0;
   JpsiTau_hammer_ebe_e9_up_lattice = 0;
   JpsiTau_hammer_ebe_e9_down_lattice = 0;
   JpsiTau_hammer_ebe_e10_up_lattice = 0;
   JpsiTau_hammer_ebe_e10_down_lattice = 0;
   JpsiTau_hammer_ebe_e11_up_lattice = 0;
   JpsiTau_hammer_ebe_e11_down_lattice = 0;
   JpsiTau_hammer_ebe_e12_up_lattice = 0;
   JpsiTau_hammer_ebe_e12_down_lattice = 0;
   JpsiTau_hammer_ebe_e13_up_lattice = 0;
   JpsiTau_hammer_ebe_e13_down_lattice = 0;
   JpsiTau_hammer_ebe_e14_up_lattice = 0;
   JpsiTau_hammer_ebe_e14_down_lattice = 0;
   JpsiTau_gen_pion_pt = 0;
   JpsiTau_gen_pion_eta = 0;
   JpsiTau_gen_pion_phi = 0;
   JpsiTau_gen_pion_matched = 0;
   JpsiTau_gen_tau_pt = 0;
   JpsiTau_gen_tau_eta = 0;
   JpsiTau_gen_tau_phi = 0;
   JpsiTau_gen_tau_nprong = 0;
   JpsiTau_gen_tau_nmatched = 0;
   JpsiTau_st_isBdecay = 0;
   JpsiTau_st_isBdecaypdg = 0;
   JpsiTau_st_isBdecayppdg = 0;
   JpsiTau_st_isSignal = 0;
   JpsiTau_st_nprong = 0;
   JpsiTau_st_nprong_pi0 = 0;
   JpsiTau_st_idx = 0;
   JpsiTau_st_doca3d = 0;
   JpsiTau_st_doca2d = 0;
   JpsiTau_st_doca3ds = 0;
   JpsiTau_st_doca2ds = 0;
   JpsiTau_st_doca3de = 0;
   JpsiTau_st_doca2de = 0;
   JpsiTau_st_dz = 0;
   JpsiTau_st_isAssociate = 0;
   JpsiTau_st_near_dz = 0;
   JpsiTau_st_dr_jpsi = 0;
   JpsiTau_st_trigMatch = 0;
   JpsiTau_st_trigMatch_dr = 0;
   JpsiTau_st_pvAssociationQuality = 0;
   JpsiTau_st_pt = 0;
   JpsiTau_st_eta = 0;
   JpsiTau_st_phi = 0;
   JpsiTau_st_charge = 0;
   JpsiTau_st_mass = 0;
   JpsiTau_st_dnn = 0;
   JpsiTau_st_dnn_1prong = 0;
   JpsiTau_st_dnn_otherB = 0;
   JpsiTau_st_dnn_pu = 0;
   JpsiTau_st_matchidx = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("genParticle_N", &genParticle_N, &b_genParticle_N);
   fChain->SetBranchAddress("genParticle_pt", &genParticle_pt, &b_genParticle_pt);
   fChain->SetBranchAddress("genParticle_eta", &genParticle_eta, &b_genParticle_eta);
   fChain->SetBranchAddress("genParticle_phi", &genParticle_phi, &b_genParticle_phi);
   fChain->SetBranchAddress("genParticle_mass", &genParticle_mass, &b_genParticle_mass);
   fChain->SetBranchAddress("genParticle_pdgId", &genParticle_pdgId, &b_genParticle_pdgId);
   fChain->SetBranchAddress("genParticle_status", &genParticle_status, &b_genParticle_status);
   fChain->SetBranchAddress("genParticle_isPrompt", &genParticle_isPrompt, &b_genParticle_isPrompt);
   fChain->SetBranchAddress("genParticle_isDirectPromptTauDecayProduct", &genParticle_isDirectPromptTauDecayProduct, &b_genParticle_isDirectPromptTauDecayProduct);
   fChain->SetBranchAddress("genParticle_isDirectHardProcessTauDecayProductFinalState", &genParticle_isDirectHardProcessTauDecayProductFinalState, &b_genParticle_isDirectHardProcessTauDecayProductFinalState);
   fChain->SetBranchAddress("genParticle_fromHardProcessFinalState", &genParticle_fromHardProcessFinalState, &b_genParticle_fromHardProcessFinalState);
   fChain->SetBranchAddress("genParticle_mother", &genParticle_mother, &b_genParticle_mother);
   fChain->SetBranchAddress("genParticle_mother_pt", &genParticle_mother_pt, &b_genParticle_mother_pt);
   fChain->SetBranchAddress("genParticle_nMoth", &genParticle_nMoth, &b_genParticle_nMoth);
   fChain->SetBranchAddress("genParticle_nDau", &genParticle_nDau, &b_genParticle_nDau);
   fChain->SetBranchAddress("genParticle_dau", &genParticle_dau, &b_genParticle_dau);
   fChain->SetBranchAddress("genParticle_pdgs", &genParticle_pdgs, &b_genParticle_pdgs);
   fChain->SetBranchAddress("genParticle_layers", &genParticle_layers, &b_genParticle_layers);
   fChain->SetBranchAddress("genParticle_ppt", &genParticle_ppt, &b_genParticle_ppt);
   fChain->SetBranchAddress("genParticle_peta", &genParticle_peta, &b_genParticle_peta);
   fChain->SetBranchAddress("genParticle_pphi", &genParticle_pphi, &b_genParticle_pphi);
   fChain->SetBranchAddress("genParticle_isfinal", &genParticle_isfinal, &b_genParticle_isfinal);
   fChain->SetBranchAddress("lheV_pt", &lheV_pt, &b_lheV_pt);
   fChain->SetBranchAddress("lheHT", &lheHT, &b_lheHT);
   fChain->SetBranchAddress("lheNj", &lheNj, &b_lheNj);
   fChain->SetBranchAddress("lheNb", &lheNb, &b_lheNb);
   fChain->SetBranchAddress("lheNl", &lheNl, &b_lheNl);
   fChain->SetBranchAddress("lheV_mass", &lheV_mass, &b_lheV_mass);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genFacWeightUp", &genFacWeightUp, &b_genFacWeightUp);
   fChain->SetBranchAddress("genFacWeightDown", &genFacWeightDown, &b_genFacWeightDown);
   fChain->SetBranchAddress("genRenWeightUp", &genRenWeightUp, &b_genRenWeightUp);
   fChain->SetBranchAddress("genRenWeightDown", &genRenWeightDown, &b_genRenWeightDown);
   fChain->SetBranchAddress("genFacRenWeightUp", &genFacRenWeightUp, &b_genFacRenWeightUp);
   fChain->SetBranchAddress("genFacRenWeightDown", &genFacRenWeightDown, &b_genFacRenWeightDown);
   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
   fChain->SetBranchAddress("PDF_rms", &PDF_rms, &b_PDF_rms);
   fChain->SetBranchAddress("PDF_x", &PDF_x, &b_PDF_x);
   fChain->SetBranchAddress("PDF_xPDF", &PDF_xPDF, &b_PDF_xPDF);
   fChain->SetBranchAddress("PDF_id", &PDF_id, &b_PDF_id);
   fChain->SetBranchAddress("EVENT_event", &EVENT_event, &b_EVENT_event);
   fChain->SetBranchAddress("EVENT_run", &EVENT_run, &b_EVENT_run);
   fChain->SetBranchAddress("EVENT_lumiBlock", &EVENT_lumiBlock, &b_EVENT_lumiBlock);
   fChain->SetBranchAddress("nPuVtxTrue", &nPuVtxTrue, &b_nPuVtxTrue);
   fChain->SetBranchAddress("nPuVtx", &nPuVtx, &b_nPuVtx);
   fChain->SetBranchAddress("bX", &bX, &b_bX);
   fChain->SetBranchAddress("PV_N", &PV_N, &b_PV_N);
   fChain->SetBranchAddress("PV_filter", &PV_filter, &b_PV_filter);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_rho", &PV_rho, &b_PV_rho);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("BeamSpot_x0", &BeamSpot_x0, &b_BeamSpot_x0);
   fChain->SetBranchAddress("BeamSpot_y0", &BeamSpot_y0, &b_BeamSpot_y0);
   fChain->SetBranchAddress("BeamSpot_z0", &BeamSpot_z0, &b_BeamSpot_z0);
   fChain->SetBranchAddress("HLT_BPH_isFired", &HLT_BPH_isFired_, &b_HLT_BPH_isFired_);
   fChain->SetBranchAddress("HLT_BPH_isFired.first", HLT_BPH_isFired_first, &b_HLT_BPH_isFired_first);
   fChain->SetBranchAddress("HLT_BPH_isFired.second", HLT_BPH_isFired_second, &b_HLT_BPH_isFired_second);
   fChain->SetBranchAddress("JpsiTau_nCandidates", &JpsiTau_nCandidates, &b_JpsiTau_nCandidates);
   fChain->SetBranchAddress("JpsiTau_isJpsiMu", &JpsiTau_isJpsiMu, &b_JpsiTau_isJpsiMu);
   fChain->SetBranchAddress("JpsiTau_isJpsiTau2Mu", &JpsiTau_isJpsiTau2Mu, &b_JpsiTau_isJpsiTau2Mu);
   fChain->SetBranchAddress("truth_tau_dipion1_mass", &truth_tau_dipion1_mass, &b_truth_tau_dipion1_mass);
   fChain->SetBranchAddress("truth_tau_dipion1_pt", &truth_tau_dipion1_pt, &b_truth_tau_dipion1_pt);
   fChain->SetBranchAddress("truth_tau_dipion1_eta", &truth_tau_dipion1_eta, &b_truth_tau_dipion1_eta);
   fChain->SetBranchAddress("truth_tau_dipion1_phi", &truth_tau_dipion1_phi, &b_truth_tau_dipion1_phi);
   fChain->SetBranchAddress("truth_tau_dipion2_mass", &truth_tau_dipion2_mass, &b_truth_tau_dipion2_mass);
   fChain->SetBranchAddress("truth_tau_dipion2_pt", &truth_tau_dipion2_pt, &b_truth_tau_dipion2_pt);
   fChain->SetBranchAddress("truth_tau_dipion2_eta", &truth_tau_dipion2_eta, &b_truth_tau_dipion2_eta);
   fChain->SetBranchAddress("truth_tau_dipion2_phi", &truth_tau_dipion2_phi, &b_truth_tau_dipion2_phi);
   fChain->SetBranchAddress("JpsiTau_mu1_pt", &JpsiTau_mu1_pt, &b_JpsiTau_mu1_pt);
   fChain->SetBranchAddress("JpsiTau_mu1_eta", &JpsiTau_mu1_eta, &b_JpsiTau_mu1_eta);
   fChain->SetBranchAddress("JpsiTau_mu1_phi", &JpsiTau_mu1_phi, &b_JpsiTau_mu1_phi);
   fChain->SetBranchAddress("JpsiTau_mu1_mass", &JpsiTau_mu1_mass, &b_JpsiTau_mu1_mass);
   fChain->SetBranchAddress("JpsiTau_mu1_q", &JpsiTau_mu1_q, &b_JpsiTau_mu1_q);
   fChain->SetBranchAddress("JpsiTau_mu1_isLoose", &JpsiTau_mu1_isLoose, &b_JpsiTau_mu1_isLoose);
   fChain->SetBranchAddress("JpsiTau_mu1_isTight", &JpsiTau_mu1_isTight, &b_JpsiTau_mu1_isTight);
   fChain->SetBranchAddress("JpsiTau_mu1_isPF", &JpsiTau_mu1_isPF, &b_JpsiTau_mu1_isPF);
   fChain->SetBranchAddress("JpsiTau_mu1_isGlobal", &JpsiTau_mu1_isGlobal, &b_JpsiTau_mu1_isGlobal);
   fChain->SetBranchAddress("JpsiTau_mu1_isTracker", &JpsiTau_mu1_isTracker, &b_JpsiTau_mu1_isTracker);
   fChain->SetBranchAddress("JpsiTau_mu1_isSoft", &JpsiTau_mu1_isSoft, &b_JpsiTau_mu1_isSoft);
   fChain->SetBranchAddress("JpsiTau_mu1_vx", &JpsiTau_mu1_vx, &b_JpsiTau_mu1_vx);
   fChain->SetBranchAddress("JpsiTau_mu1_vy", &JpsiTau_mu1_vy, &b_JpsiTau_mu1_vy);
   fChain->SetBranchAddress("JpsiTau_mu1_vz", &JpsiTau_mu1_vz, &b_JpsiTau_mu1_vz);
   fChain->SetBranchAddress("JpsiTau_mu1_dbiso", &JpsiTau_mu1_dbiso, &b_JpsiTau_mu1_dbiso);
   fChain->SetBranchAddress("JpsiTau_mu2_pt", &JpsiTau_mu2_pt, &b_JpsiTau_mu2_pt);
   fChain->SetBranchAddress("JpsiTau_mu2_eta", &JpsiTau_mu2_eta, &b_JpsiTau_mu2_eta);
   fChain->SetBranchAddress("JpsiTau_mu2_phi", &JpsiTau_mu2_phi, &b_JpsiTau_mu2_phi);
   fChain->SetBranchAddress("JpsiTau_mu2_mass", &JpsiTau_mu2_mass, &b_JpsiTau_mu2_mass);
   fChain->SetBranchAddress("JpsiTau_mu2_q", &JpsiTau_mu2_q, &b_JpsiTau_mu2_q);
   fChain->SetBranchAddress("JpsiTau_mu2_isLoose", &JpsiTau_mu2_isLoose, &b_JpsiTau_mu2_isLoose);
   fChain->SetBranchAddress("JpsiTau_mu2_isTight", &JpsiTau_mu2_isTight, &b_JpsiTau_mu2_isTight);
   fChain->SetBranchAddress("JpsiTau_mu2_isPF", &JpsiTau_mu2_isPF, &b_JpsiTau_mu2_isPF);
   fChain->SetBranchAddress("JpsiTau_mu2_isGlobal", &JpsiTau_mu2_isGlobal, &b_JpsiTau_mu2_isGlobal);
   fChain->SetBranchAddress("JpsiTau_mu2_isTracker", &JpsiTau_mu2_isTracker, &b_JpsiTau_mu2_isTracker);
   fChain->SetBranchAddress("JpsiTau_mu2_isSoft", &JpsiTau_mu2_isSoft, &b_JpsiTau_mu2_isSoft);
   fChain->SetBranchAddress("JpsiTau_mu2_vx", &JpsiTau_mu2_vx, &b_JpsiTau_mu2_vx);
   fChain->SetBranchAddress("JpsiTau_mu2_vy", &JpsiTau_mu2_vy, &b_JpsiTau_mu2_vy);
   fChain->SetBranchAddress("JpsiTau_mu2_vz", &JpsiTau_mu2_vz, &b_JpsiTau_mu2_vz);
   fChain->SetBranchAddress("JpsiTau_mu2_dbiso", &JpsiTau_mu2_dbiso, &b_JpsiTau_mu2_dbiso);
   fChain->SetBranchAddress("JpsiTau_tau_pt", &JpsiTau_tau_pt, &b_JpsiTau_tau_pt);
   fChain->SetBranchAddress("JpsiTau_tau_eta", &JpsiTau_tau_eta, &b_JpsiTau_tau_eta);
   fChain->SetBranchAddress("JpsiTau_tau_phi", &JpsiTau_tau_phi, &b_JpsiTau_tau_phi);
   fChain->SetBranchAddress("JpsiTau_tau_mass", &JpsiTau_tau_mass, &b_JpsiTau_tau_mass);
   fChain->SetBranchAddress("JpsiTau_tau_q", &JpsiTau_tau_q, &b_JpsiTau_tau_q);
   fChain->SetBranchAddress("JpsiTau_tau_vx", &JpsiTau_tau_vx, &b_JpsiTau_tau_vx);
   fChain->SetBranchAddress("JpsiTau_tau_vy", &JpsiTau_tau_vy, &b_JpsiTau_tau_vy);
   fChain->SetBranchAddress("JpsiTau_tau_vz", &JpsiTau_tau_vz, &b_JpsiTau_tau_vz);
   fChain->SetBranchAddress("JpsiTau_tau_max_dr_3prong", &JpsiTau_tau_max_dr_3prong, &b_JpsiTau_tau_max_dr_3prong);
   fChain->SetBranchAddress("JpsiTau_tau_lip", &JpsiTau_tau_lip, &b_JpsiTau_tau_lip);
   fChain->SetBranchAddress("JpsiTau_tau_lips", &JpsiTau_tau_lips, &b_JpsiTau_tau_lips);
   fChain->SetBranchAddress("JpsiTau_tau_pvip", &JpsiTau_tau_pvip, &b_JpsiTau_tau_pvip);
   fChain->SetBranchAddress("JpsiTau_tau_pvips", &JpsiTau_tau_pvips, &b_JpsiTau_tau_pvips);
   fChain->SetBranchAddress("JpsiTau_tau_fl3d", &JpsiTau_tau_fl3d, &b_JpsiTau_tau_fl3d);
   fChain->SetBranchAddress("JpsiTau_tau_fls3d", &JpsiTau_tau_fls3d, &b_JpsiTau_tau_fls3d);
   fChain->SetBranchAddress("JpsiTau_tau_alpha", &JpsiTau_tau_alpha, &b_JpsiTau_tau_alpha);
   fChain->SetBranchAddress("JpsiTau_tau_vprob", &JpsiTau_tau_vprob, &b_JpsiTau_tau_vprob);
   fChain->SetBranchAddress("JpsiTau_tau_fl3d_wjpsi", &JpsiTau_tau_fl3d_wjpsi, &b_JpsiTau_tau_fl3d_wjpsi);
   fChain->SetBranchAddress("JpsiTau_tau_fls3d_wjpsi", &JpsiTau_tau_fls3d_wjpsi, &b_JpsiTau_tau_fls3d_wjpsi);
   fChain->SetBranchAddress("JpsiTau_tau_sumofdnn", &JpsiTau_tau_sumofdnn, &b_JpsiTau_tau_sumofdnn);
   fChain->SetBranchAddress("JpsiTau_tau_sumofdnn_1prong", &JpsiTau_tau_sumofdnn_1prong, &b_JpsiTau_tau_sumofdnn_1prong);
   fChain->SetBranchAddress("JpsiTau_tau_sumofdnn_otherB", &JpsiTau_tau_sumofdnn_otherB, &b_JpsiTau_tau_sumofdnn_otherB);
   fChain->SetBranchAddress("JpsiTau_tau_sumofdnn_pu", &JpsiTau_tau_sumofdnn_pu, &b_JpsiTau_tau_sumofdnn_pu);
   fChain->SetBranchAddress("JpsiTau_tau_rhomass1", &JpsiTau_tau_rhomass1, &b_JpsiTau_tau_rhomass1);
   fChain->SetBranchAddress("JpsiTau_tau_rhomass2", &JpsiTau_tau_rhomass2, &b_JpsiTau_tau_rhomass2);
   fChain->SetBranchAddress("JpsiTau_tau_rhopt1", &JpsiTau_tau_rhopt1, &b_JpsiTau_tau_rhopt1);
   fChain->SetBranchAddress("JpsiTau_tau_rhopt2", &JpsiTau_tau_rhopt2, &b_JpsiTau_tau_rhopt2);
   fChain->SetBranchAddress("JpsiTau_tau_rhomass_ss", &JpsiTau_tau_rhomass_ss, &b_JpsiTau_tau_rhomass_ss);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_pt", &JpsiTau_tau_pi1_pt, &b_JpsiTau_tau_pi1_pt);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_eta", &JpsiTau_tau_pi1_eta, &b_JpsiTau_tau_pi1_eta);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_phi", &JpsiTau_tau_pi1_phi, &b_JpsiTau_tau_pi1_phi);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_mass", &JpsiTau_tau_pi1_mass, &b_JpsiTau_tau_pi1_mass);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_q", &JpsiTau_tau_pi1_q, &b_JpsiTau_tau_pi1_q);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_doca3d", &JpsiTau_tau_pi1_doca3d, &b_JpsiTau_tau_pi1_doca3d);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_doca3de", &JpsiTau_tau_pi1_doca3de, &b_JpsiTau_tau_pi1_doca3de);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_doca2d", &JpsiTau_tau_pi1_doca2d, &b_JpsiTau_tau_pi1_doca2d);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_doca2de", &JpsiTau_tau_pi1_doca2de, &b_JpsiTau_tau_pi1_doca2de);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_dz", &JpsiTau_tau_pi1_dz, &b_JpsiTau_tau_pi1_dz);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_near_dz", &JpsiTau_tau_pi1_near_dz, &b_JpsiTau_tau_pi1_near_dz);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_isAssociate", &JpsiTau_tau_pi1_isAssociate, &b_JpsiTau_tau_pi1_isAssociate);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_pvAssociationQuality", &JpsiTau_tau_pi1_pvAssociationQuality, &b_JpsiTau_tau_pi1_pvAssociationQuality);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_isBdecay", &JpsiTau_tau_pi1_isBdecay, &b_JpsiTau_tau_pi1_isBdecay);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_isBdecaypdg", &JpsiTau_tau_pi1_isBdecaypdg, &b_JpsiTau_tau_pi1_isBdecaypdg);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_isBdecayppdg", &JpsiTau_tau_pi1_isBdecayppdg, &b_JpsiTau_tau_pi1_isBdecayppdg);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_isSignal", &JpsiTau_tau_pi1_isSignal, &b_JpsiTau_tau_pi1_isSignal);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_nprong", &JpsiTau_tau_pi1_nprong, &b_JpsiTau_tau_pi1_nprong);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_nprong_pi0", &JpsiTau_tau_pi1_nprong_pi0, &b_JpsiTau_tau_pi1_nprong_pi0);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_dnn", &JpsiTau_tau_pi1_dnn, &b_JpsiTau_tau_pi1_dnn);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_dnn_1prong", &JpsiTau_tau_pi1_dnn_1prong, &b_JpsiTau_tau_pi1_dnn_1prong);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_dnn_otherB", &JpsiTau_tau_pi1_dnn_otherB, &b_JpsiTau_tau_pi1_dnn_otherB);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_dnn_pu", &JpsiTau_tau_pi1_dnn_pu, &b_JpsiTau_tau_pi1_dnn_pu);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_trigMatch", &JpsiTau_tau_pi1_trigMatch, &b_JpsiTau_tau_pi1_trigMatch);
   fChain->SetBranchAddress("JpsiTau_tau_pi1_trigMatch_dr", &JpsiTau_tau_pi1_trigMatch_dr, &b_JpsiTau_tau_pi1_trigMatch_dr);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_pt", &JpsiTau_tau_pi2_pt, &b_JpsiTau_tau_pi2_pt);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_eta", &JpsiTau_tau_pi2_eta, &b_JpsiTau_tau_pi2_eta);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_phi", &JpsiTau_tau_pi2_phi, &b_JpsiTau_tau_pi2_phi);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_mass", &JpsiTau_tau_pi2_mass, &b_JpsiTau_tau_pi2_mass);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_q", &JpsiTau_tau_pi2_q, &b_JpsiTau_tau_pi2_q);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_doca3d", &JpsiTau_tau_pi2_doca3d, &b_JpsiTau_tau_pi2_doca3d);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_doca3de", &JpsiTau_tau_pi2_doca3de, &b_JpsiTau_tau_pi2_doca3de);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_doca2d", &JpsiTau_tau_pi2_doca2d, &b_JpsiTau_tau_pi2_doca2d);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_doca2de", &JpsiTau_tau_pi2_doca2de, &b_JpsiTau_tau_pi2_doca2de);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_dz", &JpsiTau_tau_pi2_dz, &b_JpsiTau_tau_pi2_dz);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_near_dz", &JpsiTau_tau_pi2_near_dz, &b_JpsiTau_tau_pi2_near_dz);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_isAssociate", &JpsiTau_tau_pi2_isAssociate, &b_JpsiTau_tau_pi2_isAssociate);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_pvAssociationQuality", &JpsiTau_tau_pi2_pvAssociationQuality, &b_JpsiTau_tau_pi2_pvAssociationQuality);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_isBdecay", &JpsiTau_tau_pi2_isBdecay, &b_JpsiTau_tau_pi2_isBdecay);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_isBdecaypdg", &JpsiTau_tau_pi2_isBdecaypdg, &b_JpsiTau_tau_pi2_isBdecaypdg);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_isBdecayppdg", &JpsiTau_tau_pi2_isBdecayppdg, &b_JpsiTau_tau_pi2_isBdecayppdg);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_isSignal", &JpsiTau_tau_pi2_isSignal, &b_JpsiTau_tau_pi2_isSignal);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_nprong", &JpsiTau_tau_pi2_nprong, &b_JpsiTau_tau_pi2_nprong);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_nprong_pi0", &JpsiTau_tau_pi2_nprong_pi0, &b_JpsiTau_tau_pi2_nprong_pi0);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_dnn", &JpsiTau_tau_pi2_dnn, &b_JpsiTau_tau_pi2_dnn);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_dnn_1prong", &JpsiTau_tau_pi2_dnn_1prong, &b_JpsiTau_tau_pi2_dnn_1prong);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_dnn_otherB", &JpsiTau_tau_pi2_dnn_otherB, &b_JpsiTau_tau_pi2_dnn_otherB);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_dnn_pu", &JpsiTau_tau_pi2_dnn_pu, &b_JpsiTau_tau_pi2_dnn_pu);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_trigMatch", &JpsiTau_tau_pi2_trigMatch, &b_JpsiTau_tau_pi2_trigMatch);
   fChain->SetBranchAddress("JpsiTau_tau_pi2_trigMatch_dr", &JpsiTau_tau_pi2_trigMatch_dr, &b_JpsiTau_tau_pi2_trigMatch_dr);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_pt", &JpsiTau_tau_pi3_pt, &b_JpsiTau_tau_pi3_pt);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_eta", &JpsiTau_tau_pi3_eta, &b_JpsiTau_tau_pi3_eta);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_phi", &JpsiTau_tau_pi3_phi, &b_JpsiTau_tau_pi3_phi);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_mass", &JpsiTau_tau_pi3_mass, &b_JpsiTau_tau_pi3_mass);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_q", &JpsiTau_tau_pi3_q, &b_JpsiTau_tau_pi3_q);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_doca3d", &JpsiTau_tau_pi3_doca3d, &b_JpsiTau_tau_pi3_doca3d);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_doca3de", &JpsiTau_tau_pi3_doca3de, &b_JpsiTau_tau_pi3_doca3de);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_doca2d", &JpsiTau_tau_pi3_doca2d, &b_JpsiTau_tau_pi3_doca2d);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_doca2de", &JpsiTau_tau_pi3_doca2de, &b_JpsiTau_tau_pi3_doca2de);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_dz", &JpsiTau_tau_pi3_dz, &b_JpsiTau_tau_pi3_dz);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_near_dz", &JpsiTau_tau_pi3_near_dz, &b_JpsiTau_tau_pi3_near_dz);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_isAssociate", &JpsiTau_tau_pi3_isAssociate, &b_JpsiTau_tau_pi3_isAssociate);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_pvAssociationQuality", &JpsiTau_tau_pi3_pvAssociationQuality, &b_JpsiTau_tau_pi3_pvAssociationQuality);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_isBdecay", &JpsiTau_tau_pi3_isBdecay, &b_JpsiTau_tau_pi3_isBdecay);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_isBdecaypdg", &JpsiTau_tau_pi3_isBdecaypdg, &b_JpsiTau_tau_pi3_isBdecaypdg);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_isBdecayppdg", &JpsiTau_tau_pi3_isBdecayppdg, &b_JpsiTau_tau_pi3_isBdecayppdg);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_isSignal", &JpsiTau_tau_pi3_isSignal, &b_JpsiTau_tau_pi3_isSignal);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_nprong", &JpsiTau_tau_pi3_nprong, &b_JpsiTau_tau_pi3_nprong);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_nprong_pi0", &JpsiTau_tau_pi3_nprong_pi0, &b_JpsiTau_tau_pi3_nprong_pi0);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_dnn", &JpsiTau_tau_pi3_dnn, &b_JpsiTau_tau_pi3_dnn);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_dnn_1prong", &JpsiTau_tau_pi3_dnn_1prong, &b_JpsiTau_tau_pi3_dnn_1prong);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_dnn_otherB", &JpsiTau_tau_pi3_dnn_otherB, &b_JpsiTau_tau_pi3_dnn_otherB);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_dnn_pu", &JpsiTau_tau_pi3_dnn_pu, &b_JpsiTau_tau_pi3_dnn_pu);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_trigMatch", &JpsiTau_tau_pi3_trigMatch, &b_JpsiTau_tau_pi3_trigMatch);
   fChain->SetBranchAddress("JpsiTau_tau_pi3_trigMatch_dr", &JpsiTau_tau_pi3_trigMatch_dr, &b_JpsiTau_tau_pi3_trigMatch_dr);
   fChain->SetBranchAddress("JpsiTau_tau_delta_chi2", &JpsiTau_tau_delta_chi2, &b_JpsiTau_tau_delta_chi2);
   fChain->SetBranchAddress("JpsiTau_tau_delta_n_ch", &JpsiTau_tau_delta_n_ch, &b_JpsiTau_tau_delta_n_ch);
   fChain->SetBranchAddress("JpsiTau_tau_delta_n_mu", &JpsiTau_tau_delta_n_mu, &b_JpsiTau_tau_delta_n_mu);
   fChain->SetBranchAddress("JpsiTau_tau_vweight", &JpsiTau_tau_vweight, &b_JpsiTau_tau_vweight);
   fChain->SetBranchAddress("JpsiTau_tau_refit_vx", &JpsiTau_tau_refit_vx, &b_JpsiTau_tau_refit_vx);
   fChain->SetBranchAddress("JpsiTau_tau_refit_vy", &JpsiTau_tau_refit_vy, &b_JpsiTau_tau_refit_vy);
   fChain->SetBranchAddress("JpsiTau_tau_refit_vz", &JpsiTau_tau_refit_vz, &b_JpsiTau_tau_refit_vz);
   fChain->SetBranchAddress("JpsiTau_tau_refit_chi2", &JpsiTau_tau_refit_chi2, &b_JpsiTau_tau_refit_chi2);
   fChain->SetBranchAddress("JpsiTau_tau_refit_ndof", &JpsiTau_tau_refit_ndof, &b_JpsiTau_tau_refit_ndof);
   fChain->SetBranchAddress("JpsiTau_tau_refit_rho", &JpsiTau_tau_refit_rho, &b_JpsiTau_tau_refit_rho);
   fChain->SetBranchAddress("JpsiTau_tau_iso", &JpsiTau_tau_iso, &b_JpsiTau_tau_iso);
   fChain->SetBranchAddress("JpsiTau_tau_iso_ntracks", &JpsiTau_tau_iso_ntracks, &b_JpsiTau_tau_iso_ntracks);
   fChain->SetBranchAddress("JpsiTau_tau_iso_mindoca", &JpsiTau_tau_iso_mindoca, &b_JpsiTau_tau_iso_mindoca);
   fChain->SetBranchAddress("JpsiTau_ptbal", &JpsiTau_ptbal, &b_JpsiTau_ptbal);
   fChain->SetBranchAddress("JpsiTau_jpsi_tau_alpha", &JpsiTau_jpsi_tau_alpha, &b_JpsiTau_jpsi_tau_alpha);
   fChain->SetBranchAddress("JpsiTau_PV_vx", &JpsiTau_PV_vx, &b_JpsiTau_PV_vx);
   fChain->SetBranchAddress("JpsiTau_PV_vy", &JpsiTau_PV_vy, &b_JpsiTau_PV_vy);
   fChain->SetBranchAddress("JpsiTau_PV_vz", &JpsiTau_PV_vz, &b_JpsiTau_PV_vz);
   fChain->SetBranchAddress("JpsiTau_bbPV_vx", &JpsiTau_bbPV_vx, &b_JpsiTau_bbPV_vx);
   fChain->SetBranchAddress("JpsiTau_bbPV_vy", &JpsiTau_bbPV_vy, &b_JpsiTau_bbPV_vy);
   fChain->SetBranchAddress("JpsiTau_bbPV_vz", &JpsiTau_bbPV_vz, &b_JpsiTau_bbPV_vz);
   fChain->SetBranchAddress("JpsiTau_bbPV_chi2", &JpsiTau_bbPV_chi2, &b_JpsiTau_bbPV_chi2);
   fChain->SetBranchAddress("JpsiTau_bbPV_ndof", &JpsiTau_bbPV_ndof, &b_JpsiTau_bbPV_ndof);
   fChain->SetBranchAddress("JpsiTau_bbPV_rho", &JpsiTau_bbPV_rho, &b_JpsiTau_bbPV_rho);
   fChain->SetBranchAddress("JpsiTau_Jpsi_pt", &JpsiTau_Jpsi_pt, &b_JpsiTau_Jpsi_pt);
   fChain->SetBranchAddress("JpsiTau_Jpsi_eta", &JpsiTau_Jpsi_eta, &b_JpsiTau_Jpsi_eta);
   fChain->SetBranchAddress("JpsiTau_Jpsi_phi", &JpsiTau_Jpsi_phi, &b_JpsiTau_Jpsi_phi);
   fChain->SetBranchAddress("JpsiTau_Jpsi_mass", &JpsiTau_Jpsi_mass, &b_JpsiTau_Jpsi_mass);
   fChain->SetBranchAddress("JpsiTau_Jpsi_vprob", &JpsiTau_Jpsi_vprob, &b_JpsiTau_Jpsi_vprob);
   fChain->SetBranchAddress("JpsiTau_Jpsi_lip", &JpsiTau_Jpsi_lip, &b_JpsiTau_Jpsi_lip);
   fChain->SetBranchAddress("JpsiTau_Jpsi_lips", &JpsiTau_Jpsi_lips, &b_JpsiTau_Jpsi_lips);
   fChain->SetBranchAddress("JpsiTau_Jpsi_pvip", &JpsiTau_Jpsi_pvip, &b_JpsiTau_Jpsi_pvip);
   fChain->SetBranchAddress("JpsiTau_Jpsi_pvips", &JpsiTau_Jpsi_pvips, &b_JpsiTau_Jpsi_pvips);
   fChain->SetBranchAddress("JpsiTau_Jpsi_fl3d", &JpsiTau_Jpsi_fl3d, &b_JpsiTau_Jpsi_fl3d);
   fChain->SetBranchAddress("JpsiTau_Jpsi_fls3d", &JpsiTau_Jpsi_fls3d, &b_JpsiTau_Jpsi_fls3d);
   fChain->SetBranchAddress("JpsiTau_Jpsi_alpha", &JpsiTau_Jpsi_alpha, &b_JpsiTau_Jpsi_alpha);
   fChain->SetBranchAddress("JpsiTau_Jpsi_maxdoca", &JpsiTau_Jpsi_maxdoca, &b_JpsiTau_Jpsi_maxdoca);
   fChain->SetBranchAddress("JpsiTau_Jpsi_mindoca", &JpsiTau_Jpsi_mindoca, &b_JpsiTau_Jpsi_mindoca);
   fChain->SetBranchAddress("JpsiTau_Jpsi_vx", &JpsiTau_Jpsi_vx, &b_JpsiTau_Jpsi_vx);
   fChain->SetBranchAddress("JpsiTau_Jpsi_vy", &JpsiTau_Jpsi_vy, &b_JpsiTau_Jpsi_vy);
   fChain->SetBranchAddress("JpsiTau_Jpsi_vz", &JpsiTau_Jpsi_vz, &b_JpsiTau_Jpsi_vz);
   fChain->SetBranchAddress("JpsiTau_B_pt", &JpsiTau_B_pt, &b_JpsiTau_B_pt);
   fChain->SetBranchAddress("JpsiTau_B_eta", &JpsiTau_B_eta, &b_JpsiTau_B_eta);
   fChain->SetBranchAddress("JpsiTau_B_phi", &JpsiTau_B_phi, &b_JpsiTau_B_phi);
   fChain->SetBranchAddress("JpsiTau_B_mass", &JpsiTau_B_mass, &b_JpsiTau_B_mass);
   fChain->SetBranchAddress("JpsiTau_B_mcorr", &JpsiTau_B_mcorr, &b_JpsiTau_B_mcorr);
   fChain->SetBranchAddress("JpsiTau_B_vprob", &JpsiTau_B_vprob, &b_JpsiTau_B_vprob);
   fChain->SetBranchAddress("JpsiTau_B_lip", &JpsiTau_B_lip, &b_JpsiTau_B_lip);
   fChain->SetBranchAddress("JpsiTau_B_lips", &JpsiTau_B_lips, &b_JpsiTau_B_lips);
   fChain->SetBranchAddress("JpsiTau_B_pvip", &JpsiTau_B_pvip, &b_JpsiTau_B_pvip);
   fChain->SetBranchAddress("JpsiTau_B_pvips", &JpsiTau_B_pvips, &b_JpsiTau_B_pvips);
   fChain->SetBranchAddress("JpsiTau_B_fl3d", &JpsiTau_B_fl3d, &b_JpsiTau_B_fl3d);
   fChain->SetBranchAddress("JpsiTau_B_fls3d", &JpsiTau_B_fls3d, &b_JpsiTau_B_fls3d);
   fChain->SetBranchAddress("JpsiTau_B_alpha", &JpsiTau_B_alpha, &b_JpsiTau_B_alpha);
   fChain->SetBranchAddress("JpsiTau_B_maxdoca", &JpsiTau_B_maxdoca, &b_JpsiTau_B_maxdoca);
   fChain->SetBranchAddress("JpsiTau_B_mindoca", &JpsiTau_B_mindoca, &b_JpsiTau_B_mindoca);
   fChain->SetBranchAddress("JpsiTau_B_vx", &JpsiTau_B_vx, &b_JpsiTau_B_vx);
   fChain->SetBranchAddress("JpsiTau_B_vy", &JpsiTau_B_vy, &b_JpsiTau_B_vy);
   fChain->SetBranchAddress("JpsiTau_B_vz", &JpsiTau_B_vz, &b_JpsiTau_B_vz);
   fChain->SetBranchAddress("JpsiTau_B_q2", &JpsiTau_B_q2, &b_JpsiTau_B_q2);
   fChain->SetBranchAddress("JpsiTau_B_mm2", &JpsiTau_B_mm2, &b_JpsiTau_B_mm2);
   fChain->SetBranchAddress("JpsiTau_B_ptmiss", &JpsiTau_B_ptmiss, &b_JpsiTau_B_ptmiss);
   fChain->SetBranchAddress("JpsiTau_B_Es", &JpsiTau_B_Es, &b_JpsiTau_B_Es);
   fChain->SetBranchAddress("JpsiTau_B_ptback", &JpsiTau_B_ptback, &b_JpsiTau_B_ptback);
   fChain->SetBranchAddress("JpsiTau_B_pt_simple", &JpsiTau_B_pt_simple, &b_JpsiTau_B_pt_simple);
   fChain->SetBranchAddress("JpsiTau_B_eta_simple", &JpsiTau_B_eta_simple, &b_JpsiTau_B_eta_simple);
   fChain->SetBranchAddress("JpsiTau_B_phi_simple", &JpsiTau_B_phi_simple, &b_JpsiTau_B_phi_simple);
   fChain->SetBranchAddress("JpsiTau_B_mass_simple", &JpsiTau_B_mass_simple, &b_JpsiTau_B_mass_simple);
   fChain->SetBranchAddress("JpsiTau_B_q2_simple", &JpsiTau_B_q2_simple, &b_JpsiTau_B_q2_simple);
   fChain->SetBranchAddress("JpsiTau_B_mm2_simple", &JpsiTau_B_mm2_simple, &b_JpsiTau_B_mm2_simple);
   fChain->SetBranchAddress("JpsiTau_B_ptmiss_simple", &JpsiTau_B_ptmiss_simple, &b_JpsiTau_B_ptmiss_simple);
   fChain->SetBranchAddress("JpsiTau_B_Es_simple", &JpsiTau_B_Es_simple, &b_JpsiTau_B_Es_simple);
   fChain->SetBranchAddress("JpsiTau_B_ptback_simple", &JpsiTau_B_ptback_simple, &b_JpsiTau_B_ptback_simple);
   fChain->SetBranchAddress("genWeightBkgB", &genWeightBkgB, &b_genWeightBkgB);
   fChain->SetBranchAddress("JpsiTau_genPV_vx", &JpsiTau_genPV_vx, &b_JpsiTau_genPV_vx);
   fChain->SetBranchAddress("JpsiTau_genPV_vy", &JpsiTau_genPV_vy, &b_JpsiTau_genPV_vy);
   fChain->SetBranchAddress("JpsiTau_genPV_vz", &JpsiTau_genPV_vz, &b_JpsiTau_genPV_vz);
   fChain->SetBranchAddress("JpsiTau_genSV_vx", &JpsiTau_genSV_vx, &b_JpsiTau_genSV_vx);
   fChain->SetBranchAddress("JpsiTau_genSV_vy", &JpsiTau_genSV_vy, &b_JpsiTau_genSV_vy);
   fChain->SetBranchAddress("JpsiTau_genSV_vz", &JpsiTau_genSV_vz, &b_JpsiTau_genSV_vz);
   fChain->SetBranchAddress("JpsiTau_ngenmuons", &JpsiTau_ngenmuons, &b_JpsiTau_ngenmuons);
   fChain->SetBranchAddress("JpsiTau_isgenmatched", &JpsiTau_isgenmatched, &b_JpsiTau_isgenmatched);
   fChain->SetBranchAddress("JpsiTau_q2_gen", &JpsiTau_q2_gen, &b_JpsiTau_q2_gen);
   fChain->SetBranchAddress("JpsiTau_nBc", &JpsiTau_nBc, &b_JpsiTau_nBc);
   fChain->SetBranchAddress("JpsiTau_B_pt_gen", &JpsiTau_B_pt_gen, &b_JpsiTau_B_pt_gen);
   fChain->SetBranchAddress("JpsiTau_B_eta_gen", &JpsiTau_B_eta_gen, &b_JpsiTau_B_eta_gen);
   fChain->SetBranchAddress("JpsiTau_B_phi_gen", &JpsiTau_B_phi_gen, &b_JpsiTau_B_phi_gen);
   fChain->SetBranchAddress("JpsiTau_B_mass_gen", &JpsiTau_B_mass_gen, &b_JpsiTau_B_mass_gen);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe", &JpsiTau_hammer_ebe, &b_JpsiTau_hammer_ebe);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e0_up", &JpsiTau_hammer_ebe_e0_up, &b_JpsiTau_hammer_ebe_e0_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e0_down", &JpsiTau_hammer_ebe_e0_down, &b_JpsiTau_hammer_ebe_e0_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e1_up", &JpsiTau_hammer_ebe_e1_up, &b_JpsiTau_hammer_ebe_e1_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e1_down", &JpsiTau_hammer_ebe_e1_down, &b_JpsiTau_hammer_ebe_e1_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e2_up", &JpsiTau_hammer_ebe_e2_up, &b_JpsiTau_hammer_ebe_e2_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e2_down", &JpsiTau_hammer_ebe_e2_down, &b_JpsiTau_hammer_ebe_e2_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e3_up", &JpsiTau_hammer_ebe_e3_up, &b_JpsiTau_hammer_ebe_e3_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e3_down", &JpsiTau_hammer_ebe_e3_down, &b_JpsiTau_hammer_ebe_e3_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e4_up", &JpsiTau_hammer_ebe_e4_up, &b_JpsiTau_hammer_ebe_e4_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e4_down", &JpsiTau_hammer_ebe_e4_down, &b_JpsiTau_hammer_ebe_e4_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e5_up", &JpsiTau_hammer_ebe_e5_up, &b_JpsiTau_hammer_ebe_e5_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e5_down", &JpsiTau_hammer_ebe_e5_down, &b_JpsiTau_hammer_ebe_e5_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e6_up", &JpsiTau_hammer_ebe_e6_up, &b_JpsiTau_hammer_ebe_e6_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e6_down", &JpsiTau_hammer_ebe_e6_down, &b_JpsiTau_hammer_ebe_e6_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e7_up", &JpsiTau_hammer_ebe_e7_up, &b_JpsiTau_hammer_ebe_e7_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e7_down", &JpsiTau_hammer_ebe_e7_down, &b_JpsiTau_hammer_ebe_e7_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e8_up", &JpsiTau_hammer_ebe_e8_up, &b_JpsiTau_hammer_ebe_e8_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e8_down", &JpsiTau_hammer_ebe_e8_down, &b_JpsiTau_hammer_ebe_e8_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e9_up", &JpsiTau_hammer_ebe_e9_up, &b_JpsiTau_hammer_ebe_e9_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e9_down", &JpsiTau_hammer_ebe_e9_down, &b_JpsiTau_hammer_ebe_e9_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e10_up", &JpsiTau_hammer_ebe_e10_up, &b_JpsiTau_hammer_ebe_e10_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e10_down", &JpsiTau_hammer_ebe_e10_down, &b_JpsiTau_hammer_ebe_e10_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e11_up", &JpsiTau_hammer_ebe_e11_up, &b_JpsiTau_hammer_ebe_e11_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e11_down", &JpsiTau_hammer_ebe_e11_down, &b_JpsiTau_hammer_ebe_e11_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e12_up", &JpsiTau_hammer_ebe_e12_up, &b_JpsiTau_hammer_ebe_e12_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e12_down", &JpsiTau_hammer_ebe_e12_down, &b_JpsiTau_hammer_ebe_e12_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e13_up", &JpsiTau_hammer_ebe_e13_up, &b_JpsiTau_hammer_ebe_e13_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e13_down", &JpsiTau_hammer_ebe_e13_down, &b_JpsiTau_hammer_ebe_e13_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e14_up", &JpsiTau_hammer_ebe_e14_up, &b_JpsiTau_hammer_ebe_e14_up);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e14_down", &JpsiTau_hammer_ebe_e14_down, &b_JpsiTau_hammer_ebe_e14_down);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_lattice", &JpsiTau_hammer_ebe_lattice, &b_JpsiTau_hammer_ebe_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e0_up_lattice", &JpsiTau_hammer_ebe_e0_up_lattice, &b_JpsiTau_hammer_ebe_e0_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e0_down_lattice", &JpsiTau_hammer_ebe_e0_down_lattice, &b_JpsiTau_hammer_ebe_e0_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e1_up_lattice", &JpsiTau_hammer_ebe_e1_up_lattice, &b_JpsiTau_hammer_ebe_e1_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e1_down_lattice", &JpsiTau_hammer_ebe_e1_down_lattice, &b_JpsiTau_hammer_ebe_e1_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e2_up_lattice", &JpsiTau_hammer_ebe_e2_up_lattice, &b_JpsiTau_hammer_ebe_e2_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e2_down_lattice", &JpsiTau_hammer_ebe_e2_down_lattice, &b_JpsiTau_hammer_ebe_e2_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e3_up_lattice", &JpsiTau_hammer_ebe_e3_up_lattice, &b_JpsiTau_hammer_ebe_e3_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e3_down_lattice", &JpsiTau_hammer_ebe_e3_down_lattice, &b_JpsiTau_hammer_ebe_e3_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e4_up_lattice", &JpsiTau_hammer_ebe_e4_up_lattice, &b_JpsiTau_hammer_ebe_e4_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e4_down_lattice", &JpsiTau_hammer_ebe_e4_down_lattice, &b_JpsiTau_hammer_ebe_e4_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e5_up_lattice", &JpsiTau_hammer_ebe_e5_up_lattice, &b_JpsiTau_hammer_ebe_e5_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e5_down_lattice", &JpsiTau_hammer_ebe_e5_down_lattice, &b_JpsiTau_hammer_ebe_e5_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e6_up_lattice", &JpsiTau_hammer_ebe_e6_up_lattice, &b_JpsiTau_hammer_ebe_e6_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e6_down_lattice", &JpsiTau_hammer_ebe_e6_down_lattice, &b_JpsiTau_hammer_ebe_e6_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e7_up_lattice", &JpsiTau_hammer_ebe_e7_up_lattice, &b_JpsiTau_hammer_ebe_e7_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e7_down_lattice", &JpsiTau_hammer_ebe_e7_down_lattice, &b_JpsiTau_hammer_ebe_e7_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e8_up_lattice", &JpsiTau_hammer_ebe_e8_up_lattice, &b_JpsiTau_hammer_ebe_e8_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e8_down_lattice", &JpsiTau_hammer_ebe_e8_down_lattice, &b_JpsiTau_hammer_ebe_e8_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e9_up_lattice", &JpsiTau_hammer_ebe_e9_up_lattice, &b_JpsiTau_hammer_ebe_e9_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e9_down_lattice", &JpsiTau_hammer_ebe_e9_down_lattice, &b_JpsiTau_hammer_ebe_e9_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e10_up_lattice", &JpsiTau_hammer_ebe_e10_up_lattice, &b_JpsiTau_hammer_ebe_e10_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e10_down_lattice", &JpsiTau_hammer_ebe_e10_down_lattice, &b_JpsiTau_hammer_ebe_e10_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e11_up_lattice", &JpsiTau_hammer_ebe_e11_up_lattice, &b_JpsiTau_hammer_ebe_e11_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e11_down_lattice", &JpsiTau_hammer_ebe_e11_down_lattice, &b_JpsiTau_hammer_ebe_e11_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e12_up_lattice", &JpsiTau_hammer_ebe_e12_up_lattice, &b_JpsiTau_hammer_ebe_e12_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e12_down_lattice", &JpsiTau_hammer_ebe_e12_down_lattice, &b_JpsiTau_hammer_ebe_e12_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e13_up_lattice", &JpsiTau_hammer_ebe_e13_up_lattice, &b_JpsiTau_hammer_ebe_e13_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e13_down_lattice", &JpsiTau_hammer_ebe_e13_down_lattice, &b_JpsiTau_hammer_ebe_e13_down_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e14_up_lattice", &JpsiTau_hammer_ebe_e14_up_lattice, &b_JpsiTau_hammer_ebe_e14_up_lattice);
   fChain->SetBranchAddress("JpsiTau_hammer_ebe_e14_down_lattice", &JpsiTau_hammer_ebe_e14_down_lattice, &b_JpsiTau_hammer_ebe_e14_down_lattice);
   fChain->SetBranchAddress("JpsiTau_nch", &JpsiTau_nch, &b_JpsiTau_nch);
   fChain->SetBranchAddress("JpsiTau_nch_before", &JpsiTau_nch_before, &b_JpsiTau_nch_before);
   fChain->SetBranchAddress("JpsiTau_gen_pion_pt", &JpsiTau_gen_pion_pt, &b_JpsiTau_gen_pion_pt);
   fChain->SetBranchAddress("JpsiTau_gen_pion_eta", &JpsiTau_gen_pion_eta, &b_JpsiTau_gen_pion_eta);
   fChain->SetBranchAddress("JpsiTau_gen_pion_phi", &JpsiTau_gen_pion_phi, &b_JpsiTau_gen_pion_phi);
   fChain->SetBranchAddress("JpsiTau_gen_pion_matched", &JpsiTau_gen_pion_matched, &b_JpsiTau_gen_pion_matched);
   fChain->SetBranchAddress("JpsiTau_gen_tau_pt", &JpsiTau_gen_tau_pt, &b_JpsiTau_gen_tau_pt);
   fChain->SetBranchAddress("JpsiTau_gen_tau_eta", &JpsiTau_gen_tau_eta, &b_JpsiTau_gen_tau_eta);
   fChain->SetBranchAddress("JpsiTau_gen_tau_phi", &JpsiTau_gen_tau_phi, &b_JpsiTau_gen_tau_phi);
   fChain->SetBranchAddress("JpsiTau_gen_tau_nprong", &JpsiTau_gen_tau_nprong, &b_JpsiTau_gen_tau_nprong);
   fChain->SetBranchAddress("JpsiTau_gen_tau_nmatched", &JpsiTau_gen_tau_nmatched, &b_JpsiTau_gen_tau_nmatched);
   fChain->SetBranchAddress("JpsiTau_st_nch", &JpsiTau_st_nch, &b_JpsiTau_st_nch);
   fChain->SetBranchAddress("JpsiTau_st_nch_matched", &JpsiTau_st_nch_matched, &b_JpsiTau_st_nch_matched);
   fChain->SetBranchAddress("JpsiTau_st_n_charged_pions", &JpsiTau_st_n_charged_pions, &b_JpsiTau_st_n_charged_pions);
   fChain->SetBranchAddress("JpsiTau_st_n_neutral_pions", &JpsiTau_st_n_neutral_pions, &b_JpsiTau_st_n_neutral_pions);
   fChain->SetBranchAddress("JpsiTau_st_n_mu_decay", &JpsiTau_st_n_mu_decay, &b_JpsiTau_st_n_mu_decay);
   fChain->SetBranchAddress("JpsiTau_st_n_e_decay", &JpsiTau_st_n_e_decay, &b_JpsiTau_st_n_e_decay);
   fChain->SetBranchAddress("JpsiTau_st_n_occurance", &JpsiTau_st_n_occurance, &b_JpsiTau_st_n_occurance);
   fChain->SetBranchAddress("JpsiTau_st_decayid", &JpsiTau_st_decayid, &b_JpsiTau_st_decayid);
   fChain->SetBranchAddress("JpsiTau_st_gentau_pt", &JpsiTau_st_gentau_pt, &b_JpsiTau_st_gentau_pt);
   fChain->SetBranchAddress("JpsiTau_st_gentau_eta", &JpsiTau_st_gentau_eta, &b_JpsiTau_st_gentau_eta);
   fChain->SetBranchAddress("JpsiTau_st_gentau_phi", &JpsiTau_st_gentau_phi, &b_JpsiTau_st_gentau_phi);
   fChain->SetBranchAddress("JpsiTau_st_genjpsi_pt", &JpsiTau_st_genjpsi_pt, &b_JpsiTau_st_genjpsi_pt);
   fChain->SetBranchAddress("JpsiTau_st_genjpsi_eta", &JpsiTau_st_genjpsi_eta, &b_JpsiTau_st_genjpsi_eta);
   fChain->SetBranchAddress("JpsiTau_st_genjpsi_phi", &JpsiTau_st_genjpsi_phi, &b_JpsiTau_st_genjpsi_phi);
   fChain->SetBranchAddress("JpsiTau_general_ntau", &JpsiTau_general_ntau, &b_JpsiTau_general_ntau);
   fChain->SetBranchAddress("JpsiTau_st_isBdecay", &JpsiTau_st_isBdecay, &b_JpsiTau_st_isBdecay);
   fChain->SetBranchAddress("JpsiTau_st_isBdecaypdg", &JpsiTau_st_isBdecaypdg, &b_JpsiTau_st_isBdecaypdg);
   fChain->SetBranchAddress("JpsiTau_st_isBdecayppdg", &JpsiTau_st_isBdecayppdg, &b_JpsiTau_st_isBdecayppdg);
   fChain->SetBranchAddress("JpsiTau_st_isSignal", &JpsiTau_st_isSignal, &b_JpsiTau_st_isSignal);
   fChain->SetBranchAddress("JpsiTau_st_nprong", &JpsiTau_st_nprong, &b_JpsiTau_st_nprong);
   fChain->SetBranchAddress("JpsiTau_st_nprong_pi0", &JpsiTau_st_nprong_pi0, &b_JpsiTau_st_nprong_pi0);
   fChain->SetBranchAddress("JpsiTau_st_idx", &JpsiTau_st_idx, &b_JpsiTau_st_idx);
   fChain->SetBranchAddress("JpsiTau_st_doca3d", &JpsiTau_st_doca3d, &b_JpsiTau_st_doca3d);
   fChain->SetBranchAddress("JpsiTau_st_doca2d", &JpsiTau_st_doca2d, &b_JpsiTau_st_doca2d);
   fChain->SetBranchAddress("JpsiTau_st_doca3ds", &JpsiTau_st_doca3ds, &b_JpsiTau_st_doca3ds);
   fChain->SetBranchAddress("JpsiTau_st_doca2ds", &JpsiTau_st_doca2ds, &b_JpsiTau_st_doca2ds);
   fChain->SetBranchAddress("JpsiTau_st_doca3de", &JpsiTau_st_doca3de, &b_JpsiTau_st_doca3de);
   fChain->SetBranchAddress("JpsiTau_st_doca2de", &JpsiTau_st_doca2de, &b_JpsiTau_st_doca2de);
   fChain->SetBranchAddress("JpsiTau_st_dz", &JpsiTau_st_dz, &b_JpsiTau_st_dz);
   fChain->SetBranchAddress("JpsiTau_st_isAssociate", &JpsiTau_st_isAssociate, &b_JpsiTau_st_isAssociate);
   fChain->SetBranchAddress("JpsiTau_st_near_dz", &JpsiTau_st_near_dz, &b_JpsiTau_st_near_dz);
   fChain->SetBranchAddress("JpsiTau_st_dr_jpsi", &JpsiTau_st_dr_jpsi, &b_JpsiTau_st_dr_jpsi);
   fChain->SetBranchAddress("JpsiTau_st_trigMatch", &JpsiTau_st_trigMatch, &b_JpsiTau_st_trigMatch);
   fChain->SetBranchAddress("JpsiTau_st_trigMatch_dr", &JpsiTau_st_trigMatch_dr, &b_JpsiTau_st_trigMatch_dr);
   fChain->SetBranchAddress("JpsiTau_st_pvAssociationQuality", &JpsiTau_st_pvAssociationQuality, &b_JpsiTau_st_pvAssociationQuality);
   fChain->SetBranchAddress("JpsiTau_st_pt", &JpsiTau_st_pt, &b_JpsiTau_st_pt);
   fChain->SetBranchAddress("JpsiTau_st_eta", &JpsiTau_st_eta, &b_JpsiTau_st_eta);
   fChain->SetBranchAddress("JpsiTau_st_phi", &JpsiTau_st_phi, &b_JpsiTau_st_phi);
   fChain->SetBranchAddress("JpsiTau_st_charge", &JpsiTau_st_charge, &b_JpsiTau_st_charge);
   fChain->SetBranchAddress("JpsiTau_st_mass", &JpsiTau_st_mass, &b_JpsiTau_st_mass);
   fChain->SetBranchAddress("JpsiTau_st_dnn", &JpsiTau_st_dnn, &b_JpsiTau_st_dnn);
   fChain->SetBranchAddress("JpsiTau_st_dnn_1prong", &JpsiTau_st_dnn_1prong, &b_JpsiTau_st_dnn_1prong);
   fChain->SetBranchAddress("JpsiTau_st_dnn_otherB", &JpsiTau_st_dnn_otherB, &b_JpsiTau_st_dnn_otherB);
   fChain->SetBranchAddress("JpsiTau_st_dnn_pu", &JpsiTau_st_dnn_pu, &b_JpsiTau_st_dnn_pu);
   fChain->SetBranchAddress("JpsiTau_st_matchidx", &JpsiTau_st_matchidx, &b_JpsiTau_st_matchidx);
   fChain->SetBranchAddress("JpsiTau_perEVT_mc", &JpsiTau_perEVT_mc, &b_JpsiTau_perEVT_mc);
   fChain->SetBranchAddress("JpsiTau_perEVT_data", &JpsiTau_perEVT_data, &b_JpsiTau_perEVT_data);
   Notify();
}

Bool_t MyTauClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyTauClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyTauClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyTauClass_cxx
