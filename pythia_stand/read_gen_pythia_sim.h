
//****************************************************************************************
// Analize MC tree.
// New version filling TH2 and saving th1 matrix
// Michele Pennisi
//****************************************************************************************

#if !defined(__CINT__) || defined(__CLING__)
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TGrid.h"
#endif

Double_t weight_single_muon(TH1D *h_weigh[5], Double_t eta_muon, Double_t pt_muon);
Int_t N_HFquarks_gen; // gen c/cbar or b/bar HFquarks in the event
// Int_t N_HFquarks_rec; // gen c/cbar or b/bar HFquarks in the event

Int_t NMuons_gen; // gen muon in the event
Int_t NDimu_gen;  // gen dimuons in the event
Int_t NMuons_rec; // rec muon tracks in the event
Int_t NDimu_rec;  // rec dimuons in the event

const Int_t fMuons_dim = 1000;
const Int_t fDimu_dim = 1000;

// Int_t PDG_HFquark_rec[fMuons_dim];   // single rec c/cbar PDG mum
// Double_t Pt_HFquark_rec[fMuons_dim]; // single rec c/cbar or b/bbar HFquark pT
// Double_t Y_HFquark_rec[fMuons_dim];  // single rec c/cbar or b/bbar HFquark y

Int_t PDGmum_rec[fMuons_dim];           // single rec mu PDG mum
Double_t Pt_rec[fMuons_dim];            // single rec mu pT
Double_t E_rec[fMuons_dim];             // single rec mu E
Double_t Px_rec[fMuons_dim];            // single rec mu px
Double_t Py_rec[fMuons_dim];            // single rec mu py
Double_t Pz_rec[fMuons_dim];            // single rec mu pz
Double_t Y_rec[fMuons_dim];             // single rec mu y
Double_t Eta_rec[fMuons_dim];           // single rec mu eta
Int_t MatchTrig_rec[fMuons_dim];        // single rec mu match trigger
Double_t TrackChi2_rec[fMuons_dim];     // single rec mu chi2 track
Double_t MatchTrigChi2_rec[fMuons_dim]; // single rec mu chi2 of match trigger
Int_t Charge_rec[fMuons_dim];           // single rec mu charge
Double_t RAtAbsEnd_rec[fMuons_dim];     // single rec mu distance from beam center at end abs
Int_t pDCA_rec[fMuons_dim];             // single rec mu charge
Double_t Phi_rec[fMuons_dim];           // single rec mu phi
Double_t Theta_rec[fMuons_dim];         // single rec mu theta

Int_t PDG_HFquark_gen[fMuons_dim]; // single gen c/cbar PDG mum
Int_t PDG_HFquark_gen_daughter1[fMuons_dim];
Int_t PDG_HFquark_gen_daughter2[fMuons_dim];
Double_t Pt_HFquark_gen[fMuons_dim]; // single gen c/cbar or b/bbar HFquark pT
Double_t Y_HFquark_gen[fMuons_dim];  // single gen c/cbar or b/bbar HFquark y

Int_t PDGmum_gen[fMuons_dim];   // single gen mu PDG mum
Double_t Pt_gen[fMuons_dim];    // single gen mu pT
Double_t E_gen[fMuons_dim];     // single gen mu E
Double_t Px_gen[fMuons_dim];    // single gen mu px
Double_t Py_gen[fMuons_dim];    // single gen mu py
Double_t Pz_gen[fMuons_dim];    // single gen mu pz
Double_t Y_gen[fMuons_dim];     // single gen mu y
Double_t Eta_gen[fMuons_dim];   // single gen mu eta
Double_t Phi_gen[fMuons_dim];   // single gen mu phi
Double_t Theta_gen[fMuons_dim]; // single gen mu theta

Int_t NHadron_gen = 0;
Int_t PDGHadron_gen[fMuons_dim];
Int_t PromptHadron_gen[fMuons_dim];
Double_t Hadron_Pt_gen[fMuons_dim];
Double_t Hadron_E_gen[fMuons_dim];
Double_t Hadron_Px_gen[fMuons_dim];
Double_t Hadron_Py_gen[fMuons_dim];
Double_t Hadron_Pz_gen[fMuons_dim];
Double_t Hadron_Y_gen[fMuons_dim];
Double_t Hadron_Eta_gen[fMuons_dim];
Double_t Hadron_Phi_gen[fMuons_dim];
Double_t Hadron_Theta_gen[fMuons_dim];

Int_t DimuMu_gen[fDimu_dim][2];   // reference to single gen mus
Double_t DimuPt_gen[fDimu_dim];   // gen dimuon pT
Double_t DimuPx_gen[fDimu_dim];   // gen dimuon px
Double_t DimuPy_gen[fDimu_dim];   // gen dimuon py
Double_t DimuPz_gen[fDimu_dim];   // gen dimuon pz
Double_t DimuY_gen[fDimu_dim];    // gen dimuon y
Double_t DimuMass_gen[fDimu_dim]; // gen dimuon invariant mass
Int_t DimuCharge_gen[fDimu_dim];  // gen dimuon charge

Int_t DimuMu_rec[fDimu_dim][2];    // reference to single rec mus
Double_t DimuPt_rec[fDimu_dim];    // rec dimuon pT
Double_t DimuPx_rec[fDimu_dim];    // rec dimuon px
Double_t DimuPy_rec[fDimu_dim];    // rec dimuon py
Double_t DimuPz_rec[fDimu_dim];    // rec dimuon pz
Double_t DimuY_rec[fDimu_dim];     // rec dimuon y
Double_t DimuMass_rec[fDimu_dim];  // rec dimuon invariant mass
Int_t DimuCharge_rec[fDimu_dim];   // rec dimuon charge
Int_t DimuMatch_rec[fDimu_dim];    // rec dimuon match
Double_t DimuPhi_rec[fDimu_dim];   // rec dimuon phi
Double_t DimuTheta_rec[fDimu_dim]; // rec dimuon theta

// 0 For Generated, 1 For Generated -4.0<Y<-2.5, 2 For Generated in ALICEacc, 3 For Generated with DQ Cut, 4 For Reconstructed, 5 For Reconstructed with OLD all cut applied, 6 For Reconstructed OLD with all cut, 7 For Reconstructed with DQ CUT, 8 Reconstructed with DQ CUT and match
const Int_t n_MuCut = 9;
const Int_t n_MuCut_rec = 5;
TString name_MuCut[n_MuCut];

const Int_t n_DiMuCut = 9;
const Int_t n_DiMuCut_rec = 5;
TString name_DiMuCut[n_DiMuCut];

//-----------------------------------------------------//
const Int_t n_Mass_Cut = 9;
// 0 For no cut on Dimu Mass, 1 For M>4 GeV/c^2, 2 For 4<M<=9 GeV/c^2, 3 For 9<M<=11 GeV/c^2, 4 For 11<M<=15 GeV/c^2, 5 For M>15 GeV/c^2, 6 For 4<M<9 GeV/c^2 && 11<M<30 GeV/c^2, 7 For 11<M<30 GeV/c^2, 8 For 4<M<9 GeV/c^2 && pT<10 GeV/c
TString name_Mass_Cut[n_Mass_Cut];

Double_t Mass_Edge[n_Mass_Cut][2];

Double_t Mass_Bin[n_Mass_Cut];

//-----------------------------------------------------//
// 0 For All, 1 For HF, 2 For Charm, 3 For Beauty, 4 For LF, 5 for others
const Int_t n_MuSelection = 6;
TString name_MuSelection[n_MuSelection];

// 0 For All, 1 For HF, 2 For Charm, 3 For Beauty, 4 For HF Mixed (one muon from Charm, one muon from Beauty), 5 For LF (two muons from LF), 6 For LF Mixed (one muon from LH, one muon from HF), 7 Others
const Int_t n_DiMuSelection = 8;
TString name_DiMuSelection[n_DiMuSelection];

//-----------------------------------------------------//

// 0 for Plus, 1 for Minus, 2 for total
const Int_t n_Mu_Charge = 3;
TString name_Mu_Charge[n_Mu_Charge];

// 0 for ULS, 1 for LS++, 2 for LS--, 3 LS total
const Int_t n_DiMu_Charge = 4;
TString name_DiMu_Charge[n_DiMu_Charge];

Int_t nMu_xevent[n_Mass_Cut][n_MuCut][n_MuSelection];

Int_t nDiMu_xevent[n_Mass_Cut][n_DiMuCut][n_DiMu_Charge][n_DiMuSelection];

TH2D *h_PtYMu[n_Mass_Cut][n_MuCut][n_MuSelection];
TH1D *h_pdgMu[n_Mass_Cut][n_MuCut][n_MuSelection];
TH1D *h_nMu_xevent[n_Mass_Cut][n_MuCut][n_MuSelection];

TH1D *h_pdgDimuMu[n_Mass_Cut][n_DiMuCut][n_DiMu_Charge][n_DiMuSelection];
TH2D *h_PtYDiMu[n_Mass_Cut][n_DiMuCut][n_DiMu_Charge][n_DiMuSelection];
TH2D *h_PtMDiMu[n_Mass_Cut][n_DiMuCut][n_DiMu_Charge][n_DiMuSelection];
TH1D *h_MDiMu[n_Mass_Cut][n_DiMuCut][n_DiMu_Charge][n_DiMuSelection];
TH1D *h_nDiMu_xevent[n_Mass_Cut][n_DiMuCut][n_DiMu_Charge][n_DiMuSelection];

TH1D *h_nMu_xevent_rec_charm = new TH1D("h_nMu_xevent_rec_charm", "h_nMu_xevent_rec_charm", 10, -0.5, 9.5);
TH1D *h_nMu_xevent_rec_beauty = new TH1D("h_nMu_xevent_rec_beauty", "h_nMu_xevent_rec_beauty", 10, -0.5, 9.5);

TH1D *h_Ncharm_quark = new TH1D("h_Ncharm_quark", "h_Ncharm_quark", 60, 0, 60);
TH1D *h_pt_charm_quark = new TH1D("h_pt_charm_quark", "h_pt_charm_quark", 300, -0.5, 299.5);
TH1D *h_pt_beauty_quark = new TH1D("h_pt_beauty_quark", "h_pt_beauty_quark", 300, -0.5, 299.5);

TH1D *h_Nbeauty_quark = new TH1D("h_Nbeauty_quark", "h_Nbeauty_quark", 60, 0, 60);
TH1D *h_y_charm_quark = new TH1D("h_y_charm_quark", "h_y_charm_quark", 200, -10.0, 10.0);
TH1D *h_y_beauty_quark = new TH1D("h_y_beauty_quark", "h_y_beauty_quark", 200, -10.0, 10.0);

TH1D *h_Ncharm_antiquark = new TH1D("h_Ncharm_antiquark", "h_Ncharm_antiquark", 60, 0, 60);
TH1D *h_pt_charm_antiquark = new TH1D("h_pt_charm_antiquark", "h_pt_charm_antiquark", 300, -0.5, 299.5);
TH1D *h_pt_beauty_antiquark = new TH1D("h_pt_beauty_antiquark", "h_pt_beauty_antiquark", 300, -0.5, 299.5);

TH1D *h_Nbeauty_antiquark = new TH1D("h_Nbeauty_antiquark", "h_Nbeauty_antiquark", 300, -0.5, 60.5);
TH1D *h_y_charm_antiquark = new TH1D("h_y_charm_antiquark", "h_y_charm_antiquark", 200, -10.0, 10.0);
TH1D *h_y_beauty_antiquark = new TH1D("h_y_beauty_antiquark", "h_y_beauty_antiquark", 200, -10.0, 10.0);

TH1D *h_Nbeauty_pairs = new TH1D("h_Nbeauty_pairs", "h_Nbeauty_pairs", 60, 0.0, 60);
TH1D *h_ptbeauty_pairs = new TH1D("h_ptbeauty_pairs", "h_ptbeauty_pairs", 300, 0, 300);
TH1D *h_ybeauty_pairs = new TH1D("h_ybeauty_pairs", "h_ybeauty_pairs", 200, -10, 10);

TH1D *h_Ncharm_pairs = new TH1D("h_Ncharm_pairs", "h_Ncharm_pairs", 60, 0.0, 60);
TH1D *h_ptcharm_pairs = new TH1D("h_ptcharm_pairs", "h_ptcharm_pairs", 300, 0, 300);
TH1D *h_ycharm_pairs = new TH1D("h_ycharm_pairs", "h_ycharm_pairs", 200, -10, 10);

TH1D *h_Nbeauty_pairs_v4 = new TH1D("h_Nbeauty_pairs_v4", "h_Nbeauty_pairs_v4", 60, 0.0, 60);
TH1D *h_ptbeauty_pairs_v4 = new TH1D("h_ptbeauty_pairs_v4", "h_ptbeauty_pairs_v4", 300, 0, 300);
TH1D *h_ybeauty_pairs_v4 = new TH1D("h_ybeauty_pairs_v4", "h_ybeauty_pairs_v4", 200, -10, 10);

TH1D *h_Ncharm_pairs_v4 = new TH1D("h_Ncharm_pairs_v4", "h_Ncharm_pairs_v4", 60, 0.0, 60);
TH1D *h_ptcharm_pairs_v4 = new TH1D("h_ptcharm_pairs_v4", "h_ptcharm_pairs_v4", 300, 0, 300);
TH1D *h_ycharm_pairs_v4 = new TH1D("h_ycharm_pairs_v4", "h_ycharm_pairs_v4", 200, -10, 10);

TH1D *h_Nbeauty_pairs_v5 = new TH1D("h_Nbeauty_pairs_v5", "h_Nbeauty_pairs_v5", 60, 0.0, 60);
TH1D *h_ptbeauty_pairs_v5 = new TH1D("h_ptbeauty_pairs_v5", "h_ptbeauty_pairs_v5", 300, 0, 300);
TH1D *h_ybeauty_pairs_v5 = new TH1D("h_ybeauty_pairs_v5", "h_ybeauty_pairs_v5", 200, -10, 10);

TH1D *h_Ncharm_pairs_v5 = new TH1D("h_Ncharm_pairs_v5", "h_Ncharm_pairs_v5", 60, 0.0, 60);
TH1D *h_ptcharm_pairs_v5 = new TH1D("h_ptcharm_pairs_v5", "h_ptcharm_pairs_v5", 300, 0, 300);
TH1D *h_ycharm_pairs_v5 = new TH1D("h_ycharm_pairs_v5", "h_ycharm_pairs_v5", 200, -10, 10);

TH1D *h_Nbeauty_pairs_v6 = new TH1D("h_Nbeauty_pairs_v6", "h_Nbeauty_pairs_v6", 60, 0.0, 60);
TH1D *h_ptbeauty_pairs_v6 = new TH1D("h_ptbeauty_pairs_v6", "h_ptbeauty_pairs_v6", 300, 0, 300);
TH1D *h_ybeauty_pairs_v6 = new TH1D("h_ybeauty_pairs_v6", "h_ybeauty_pairs_v6", 200, -10, 10);

TH1D *h_Ncharm_pairs_v6 = new TH1D("h_Ncharm_pairs_v6", "h_Ncharm_pairs_v6", 60, 0.0, 60);
TH1D *h_ptcharm_pairs_v6 = new TH1D("h_ptcharm_pairs_v6", "h_ptcharm_pairs_v6", 300, 0, 300);
TH1D *h_ycharm_pairs_v6 = new TH1D("h_ycharm_pairs_v6", "h_ycharm_pairs_v6", 200, -10, 10);

const Int_t n_Hadron_studied = 4;
const Int_t n_Rapidity_studied = 6;
TString name_Hadron_studied[n_Hadron_studied];
TString name_rapidity_studied[n_Rapidity_studied];

TH1D *h_ptHadron_prompt[n_Hadron_studied][n_Rapidity_studied];
TH1D *h_yHadron_prompt[n_Hadron_studied][n_Rapidity_studied];
TH1D *h_pdgHadron_prompt[n_Hadron_studied][n_Rapidity_studied];

TH1D *h_ptHadron_Notprompt[n_Hadron_studied][n_Rapidity_studied];
TH1D *h_yHadron_Notprompt[n_Hadron_studied][n_Rapidity_studied];
TH1D *h_pdgHadron_Notprompt[n_Hadron_studied][n_Rapidity_studied];

TH1D *h_Nevents = new TH1D("h_Nevents", "h_Nevents", 2, 0, 2);

TH2D *h_Dimu_deltaphi_deltaeta_ULS = new TH2D("h_Dimu_deltaphi_deltaeta_ULS", "h_Dimu_deltaphi_deltaeta_ULS", 100, -10, 10, 100, -10, 10);
TH2D *h_Dimu_deltaphi_deltaeta_fromCharm_ULS = new TH2D("h_Dimu_deltaphi_deltaeta_fromCharm_ULS", "h_Dimu_deltaphi_deltaeta_fromCharm_ULS", 100, -10, 10, 100, -10, 10);
TH2D *h_Dimu_deltaphi_deltaeta_fromBeauty_ULS = new TH2D("h_Dimu_deltaphi_deltaeta_fromBeauty_ULS", "h_Dimu_deltaphi_deltaeta_fromBeauty_ULS", 100, -10, 10, 100, -10, 10);
TH2D *h_Dimu_deltaphi_deltaeta_fromMixed_ULS = new TH2D("h_Dimu_deltaphi_deltaeta_fromMixed_ULS", "h_Dimu_deltaphi_deltaeta_fromMixed_ULS", 100, -10, 10, 100, -10, 10);

TH2D *h_Dimu_deltaphi_deltaeta_LS = new TH2D("h_Dimu_deltaphi_deltaeta_LS", "h_Dimu_deltaphi_deltaeta_LS", 100, -10, 10, 100, -10, 10);
TH2D *h_Dimu_deltaphi_deltaeta_fromCharm_LS = new TH2D("h_Dimu_deltaphi_deltaeta_fromCharm_LS", "h_Dimu_deltaphi_deltaeta_fromCharm_LS", 100, -10, 10, 100, -10, 10);
TH2D *h_Dimu_deltaphi_deltaeta_fromBeauty_LS = new TH2D("h_Dimu_deltaphi_deltaeta_fromBeauty_LS", "h_Dimu_deltaphi_deltaeta_fromBeauty_LS", 100, -10, 10, 100, -10, 10);
TH2D *h_Dimu_deltaphi_deltaeta_fromMixed_LS = new TH2D("h_Dimu_deltaphi_deltaeta_fromMixed_LS", "h_Dimu_deltaphi_deltaeta_fromMixed_LS", 100, -10, 10, 100, -10, 10);

const Int_t eta_bins = 5;
Double_t eta_binning[eta_bins + 1] = {-4.0, -3.7, -3.4, -3.1, -2.8, -2.5};
TH1D *h_PtMuonsRec_notweightedfortrigger[eta_bins];
TH1D *h_PtMuonsRec_weightedfortrigger[eta_bins];

TH1D *h_PtDiMuonsRec_weightedfortrigger_ULS = new TH1D("h_PtDiMuonsRec_weightedfortrigger_ULS", "h_PtDiMuonsRec_weightedfortrigger_ULS", 300, 0, 30.0);
TH1D *h_PtDiMuonsRec_weightedfortrigger_fromCharm_ULS = new TH1D("h_PtDiMuonsRec_weightedfortrigger_fromCharm_ULS", "h_PtDiMuonsRec_weightedfortrigger_fromCharm_ULS", 300, 0, 30.0);
TH1D *h_PtDiMuonsRec_weightedfortrigger_fromBeauty_ULS = new TH1D("h_PtDiMuonsRec_weightedfortrigger_fromBeauty_ULS", "h_PtDiMuonsRec_weightedfortrigger_fromBeauty_ULS", 300, 0, 30.0);
TH1D *h_PtDiMuonsRec_weightedfortrigger_fromMixed_ULS = new TH1D("h_PtDiMuonsRec_weightedfortrigger_fromMixed_ULS", "h_PtDiMuonsRec_weightedfortrigger_fromMixed_ULS", 300, 0, 30.0);

TH1D *h_PtDiMuonsRec_weightedfortrigger_LS = new TH1D("h_PtDiMuonsRec_weightedfortrigger_LS", "h_PtDiMuonsRec_weightedfortrigger_LS", 300, 0, 30.0);
TH1D *h_PtDiMuonsRec_weightedfortrigger_fromCharm_LS = new TH1D("h_PtDiMuonsRec_weightedfortrigger_fromCharm_LS", "h_PtDiMuonsRec_weightedfortrigger_fromCharm_LS", 300, 0, 30.0);
TH1D *h_PtDiMuonsRec_weightedfortrigger_fromBeauty_LS = new TH1D("h_PtDiMuonsRec_weightedfortrigger_fromBeauty_LS", "h_PtDiMuonsRec_weightedfortrigger_fromBeauty_LS", 300, 0, 30.0);
TH1D *h_PtDiMuonsRec_weightedfortrigger_fromMixed_LS = new TH1D("h_PtDiMuonsRec_weightedfortrigger_fromMixed_LS", "h_PtDiMuonsRec_weightedfortrigger_fromMixed_LS", 300, 0, 30.0);

TH1D *h_MDiMuonsRec_weightedfortrigger_ULS = new TH1D("h_MDiMuonsRec_weightedfortrigger_ULS", "h_MDiMuonsRec_weightedfortrigger_ULS", 260, 4, 30.0);
TH1D *h_MDiMuonsRec_weightedfortrigger_fromCharm_ULS = new TH1D("h_MDiMuonsRec_weightedfortrigger_fromCharm_ULS", "h_MDiMuonsRec_weightedfortrigger_fromCharm_ULS", 260, 4, 30.0);
TH1D *h_MDiMuonsRec_weightedfortrigger_fromBeauty_ULS = new TH1D("h_MDiMuonsRec_weightedfortrigger_fromBeauty_ULS", "h_MDiMuonsRec_weightedfortrigger_fromBeauty_ULS", 260, 4, 30.0);
TH1D *h_MDiMuonsRec_weightedfortrigger_fromMixed_ULS = new TH1D("h_MDiMuonsRec_weightedfortrigger_fromMixed_ULS", "h_MDiMuonsRec_weightedfortrigger_fromMixed_ULS", 260, 4, 30.0);

TH1D *h_MDiMuonsRec_weightedfortrigger_LS = new TH1D("h_MDiMuonsRec_weightedfortrigger_LS", "h_MDiMuonsRec_weightedfortrigger_LS", 260, 4, 30.0);
TH1D *h_MDiMuonsRec_weightedfortrigger_fromCharm_LS = new TH1D("h_MDiMuonsRec_weightedfortrigger_fromCharm_LS", "h_MDiMuonsRec_weightedfortrigger_fromCharm_LS", 260, 4, 30.0);
TH1D *h_MDiMuonsRec_weightedfortrigger_fromBeauty_LS = new TH1D("h_MDiMuonsRec_weightedfortrigger_fromBeauty_LS", "h_MDiMuonsRec_weightedfortrigger_fromBeauty_LS", 260, 4, 30.0);
TH1D *h_MDiMuonsRec_weightedfortrigger_fromMixed_LS = new TH1D("h_MDiMuonsRec_weightedfortrigger_fromMixed_LS", "h_MDiMuonsRec_weightedfortrigger_fromMixed_LS", 260, 4, 30.0);

void SetHist()
{
  TString name_root_files[n_Hadron_studied + 1];
  name_root_files[0].Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/hf_hadron_cs/root_files/D0-Lc-Xsec-pp13.root");
  name_root_files[1].Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/hf_hadron_cs/root_files/HFPtSpectrum_13TeV_Dplus_withOldSys_RenuBala.root");
  name_root_files[2].Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/hf_hadron_cs/root_files/HFPtSpectrum_13TeV_Ds_JulienHamon.root");
  name_root_files[3].Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/hf_hadron_cs/root_files/D0-Lc-Xsec-pp13.root");
  name_root_files[4].Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/hf_hadron_cs/root_files/HEPData-ins1396331-v2-root.root");

  TString name_histo[n_Hadron_studied + 1];
  name_histo[0].Form("hxsecD0");
  name_histo[1].Form("histoSigmaCorr");
  name_histo[2].Form("histoSigmaCorr");
  name_histo[3].Form("hxsecLc");

  name_Hadron_studied[0].Form("Dzero");
  name_Hadron_studied[1].Form("Dplus");
  name_Hadron_studied[2].Form("Dstrange");
  name_Hadron_studied[3].Form("Lambda");

  name_rapidity_studied[0].Form("MidY");
  name_rapidity_studied[1].Form("FwdY1");
  name_rapidity_studied[2].Form("FwdY2");
  name_rapidity_studied[3].Form("FwdY3");
  name_rapidity_studied[4].Form("FwdY4");
  name_rapidity_studied[5].Form("FwdY5");

  for (size_t i_hadron = 0; i_hadron < n_Hadron_studied; i_hadron++)
  {
    for (size_t i_rapidity = 0; i_rapidity < n_Rapidity_studied; i_rapidity++)
    {
      if (i_rapidity == 0)
      {
        TFile *fIn_data = new TFile(Form("%s", name_root_files[i_hadron].Data()), "READ");
        TH1D *hist_data = (TH1D *)fIn_data->Get(Form("%s", name_histo[i_hadron].Data()));
        hist_data->SetDirectory(0);
        fIn_data->Close();
        h_ptHadron_prompt[i_hadron][i_rapidity] = (TH1D *)hist_data->Clone(Form("h_pt%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
        h_ptHadron_prompt[i_hadron][i_rapidity]->Reset();
        h_ptHadron_prompt[i_hadron][i_rapidity]->SetTitle(Form("h_pt%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
      }
      else
      {
        if (i_hadron < n_Hadron_studied - 1)
        {
          TFile *fIn_data = new TFile(Form("%s", name_root_files[n_Hadron_studied].Data()), "READ");
          TH1D *hist_data = (TH1D *)fIn_data->Get(Form("Table %zu/Hist1D_y%zu", i_hadron+1, i_rapidity));
          hist_data->SetDirectory(0);
          fIn_data->Close();
          h_ptHadron_prompt[i_hadron][i_rapidity] = (TH1D *)hist_data->Clone(Form("h_pt%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
          h_ptHadron_prompt[i_hadron][i_rapidity]->Reset();
          h_ptHadron_prompt[i_hadron][i_rapidity]->SetTitle(Form("h_pt%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
        }
        else
        {
          h_ptHadron_prompt[i_hadron][i_rapidity] = new TH1D(Form("h_pt%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_pt%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 300, 0, 300);
        }
      }

      h_yHadron_prompt[i_hadron][i_rapidity] = new TH1D(Form("h_y%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_y%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 200, -10, 10);
      h_pdgHadron_prompt[i_hadron][i_rapidity] = new TH1D(Form("h_pdg%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_pdg%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 12000, -6000, 6000);

      h_ptHadron_Notprompt[i_hadron][i_rapidity] = new TH1D(Form("h_pt%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_pt%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 300, 0, 300);
      h_yHadron_Notprompt[i_hadron][i_rapidity] = new TH1D(Form("h_y%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_y%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 200, -10, 10);
      h_pdgHadron_Notprompt[i_hadron][i_rapidity] = new TH1D(Form("h_pdg%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_pdg%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 12000, -6000, 6000);
    }
  }

  for (size_t w = 0; w < eta_bins; w++)
  {
    h_PtMuonsRec_notweightedfortrigger[w] = new TH1D(Form("h_PtMuonsRec_notweightedfortrigger_eta%zu", w), Form("notweighted %0.1f < #eta < %0.1f", eta_binning[w], eta_binning[w + 1]), 300, 0, 30.0);
    h_PtMuonsRec_weightedfortrigger[w] = new TH1D(Form("h_PtMuonsRec_weightedfortrigger_eta%zu", w), Form("weighted %0.1f < #eta < %0.1f", eta_binning[w], eta_binning[w + 1]), 300, 0, 30.0);
  }

  name_MuCut[0].Form("Gen"); // All single muon produced
  name_MuCut[1].Form("Gen_ycut");
  name_MuCut[2].Form("Gen_ALICEacc"); //
  name_MuCut[3].Form("Gen_DQ_cut");   //
  name_MuCut[n_MuCut - n_MuCut_rec].Form("Rec_nocut");
  name_MuCut[n_MuCut - n_MuCut_rec + 1].Form("Rec_old_cut");         // OLD cut, over single muon rapidity and matchtrig>=1
  name_MuCut[n_MuCut - n_MuCut_rec + 2].Form("Rec_old_cut_match");   // OLD cut, over single muon rapidity and matchtrig>1
  name_MuCut[n_MuCut - n_MuCut_rec + 3].Form("Rec_DQ_cut_match_AT"); // DQ cut with matchtrig>=1
  name_MuCut[n_MuCut - n_MuCut_rec + 4].Form("Rec_DQ_cut_match_LT"); // DQ cut with matchtrig>1

  name_DiMuCut[0].Form("Gen");
  name_DiMuCut[1].Form("Gen_ycut");
  name_DiMuCut[2].Form("Gen_ALICEacc");
  name_DiMuCut[3].Form("Gen_DQ_cut"); //
  name_DiMuCut[n_DiMuCut - n_DiMuCut_rec].Form("Rec_nocut");
  name_DiMuCut[n_DiMuCut - n_DiMuCut_rec + 1].Form("Rec_old_cut");         // OLD cut, over single muon rapidity, with 1o match required
  name_DiMuCut[n_DiMuCut - n_DiMuCut_rec + 2].Form("Rec_old_cut_match");   // OLD cut, over single muon rapidity, dimum1tch=2
  name_DiMuCut[n_DiMuCut - n_DiMuCut_rec + 3].Form("Rec_DQ_cut_match_AT"); // DQ cut with no match required
  name_DiMuCut[n_DiMuCut - n_DiMuCut_rec + 4].Form("Rec_DQ_cut_match_LT"); // DQ cut with dimumatch=2

  name_Mass_Cut[0].Form("M0");
  name_Mass_Cut[1].Form("M4");
  name_Mass_Cut[2].Form("LowMass");
  name_Mass_Cut[3].Form("Yres");
  name_Mass_Cut[4].Form("InterMass");
  name_Mass_Cut[5].Form("VeryHighMass");
  name_Mass_Cut[6].Form("Cut_Yres");
  name_Mass_Cut[7].Form("HighMass");
  name_Mass_Cut[8].Form("LowMass_LowPt");

  Mass_Edge[0][0] = 0.0;
  Mass_Edge[1][0] = 4.0;
  Mass_Edge[2][0] = 4.0;
  Mass_Edge[3][0] = 9.0;
  Mass_Edge[4][0] = 11.0;
  Mass_Edge[5][0] = 15.0;
  Mass_Edge[6][0] = 4.0;
  Mass_Edge[7][0] = 11.0;
  Mass_Edge[8][0] = 4.0;

  Mass_Edge[0][1] = 30.0;
  Mass_Edge[1][1] = 30.0;
  Mass_Edge[2][1] = 9.0;
  Mass_Edge[3][1] = 11.0;
  Mass_Edge[4][1] = 15.0;
  Mass_Edge[5][1] = 30.0;
  Mass_Edge[6][1] = 30.0;
  Mass_Edge[7][1] = 30.0;
  Mass_Edge[8][1] = 9.0;

  Mass_Bin[0] = 300;
  Mass_Bin[1] = 260;
  Mass_Bin[2] = 50;
  Mass_Bin[3] = 200;
  Mass_Bin[4] = 40;
  Mass_Bin[5] = 150;
  Mass_Bin[6] = 260;
  Mass_Bin[7] = 190;
  Mass_Bin[8] = 50;

  name_MuSelection[0].Form("All");
  name_MuSelection[1].Form("fromHF");
  name_MuSelection[2].Form("fromCharm");
  name_MuSelection[3].Form("fromBeauty");
  name_MuSelection[4].Form("fromLF");
  name_MuSelection[5].Form("fromOthers");

  name_DiMuSelection[0].Form("All");
  name_DiMuSelection[1].Form("fromHF");
  name_DiMuSelection[2].Form("fromCharm");
  name_DiMuSelection[3].Form("fromBeauty");
  name_DiMuSelection[4].Form("fromMixed");
  name_DiMuSelection[5].Form("fromLF");
  name_DiMuSelection[6].Form("fromLFMixed");
  name_DiMuSelection[7].Form("fromOthers");

  name_DiMu_Charge[0].Form("ULS");
  name_DiMu_Charge[1].Form("LSplus");
  name_DiMu_Charge[2].Form("LSminus");
  name_DiMu_Charge[3].Form("LS");

  for (Int_t a = 0; a < n_Mass_Cut; a++)
  {
    for (Int_t Mu_b = 0; Mu_b < n_MuCut; Mu_b++)
    {
      for (Int_t Mu_c = 0; Mu_c < n_MuSelection; Mu_c++)
      {

        h_PtYMu[a][Mu_b][Mu_c] = new TH2D(Form("h_PtYMu_%s_%s_%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data(), name_MuSelection[Mu_c].Data()), Form("h_PtYMu_%s_%s_%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data(), name_MuSelection[Mu_c].Data()), 300, 0.0, 30.0, 200, -10.0, 10.0);
        h_PtYMu[a][Mu_b][Mu_c]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h_PtYMu[a][Mu_b][Mu_c]->GetYaxis()->SetTitle("Y");

        h_pdgMu[a][Mu_b][Mu_c] = new TH1D(Form("h_pdgMu_%s_%s_%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data(), name_MuSelection[Mu_c].Data()), Form("h_pdgMu_%s_%s_%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data(), name_MuSelection[Mu_c].Data()), 12000, -5999.5, 5999.5);
        h_pdgMu[a][Mu_b][Mu_c]->GetXaxis()->SetTitle("PDG Code");
        h_pdgMu[a][Mu_b][Mu_c]->GetYaxis()->SetTitle("Counts");

        h_nMu_xevent[a][Mu_b][Mu_c] = new TH1D(Form("h_nMu_xevent_%s_%s_%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data(), name_MuSelection[Mu_c].Data()), Form("h_nMu_xevent_%s_%s_%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data(), name_MuSelection[Mu_c].Data()), 10, -0.5, 9.5);
        h_nMu_xevent[a][Mu_b][Mu_c]->GetXaxis()->SetTitle("Number of #mu x event");
        h_nMu_xevent[a][Mu_b][Mu_c]->GetYaxis()->SetTitle("Counts");

      } // End definition over Mu charge

    } // End definition over Mu selection
    for (Int_t DiMu_b = 0; DiMu_b < n_DiMuCut; DiMu_b++)
    {
      for (Int_t DiMu_c = 0; DiMu_c < n_DiMu_Charge; DiMu_c++)
      {
        for (Int_t DiMu_d = 0; DiMu_d < n_DiMuSelection; DiMu_d++)
        {

          h_PtYDiMu[a][DiMu_b][DiMu_c][DiMu_d] = new TH2D(Form("h_PtYDiMu_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), Form("h_PtYDiMu_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), 300, 0.0, 30.0, 150, -4.0, -2.5);
          h_PtYDiMu[a][DiMu_b][DiMu_c][DiMu_d]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
          h_PtYDiMu[a][DiMu_b][DiMu_c][DiMu_d]->GetYaxis()->SetTitle("Y");

          h_PtMDiMu[a][DiMu_b][DiMu_c][DiMu_d] = new TH2D(Form("h_PtMDiMu_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), Form("h_PtMDiMu_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), 300, 0.0, 30.0, Mass_Bin[a], Mass_Edge[a][0], Mass_Edge[a][1]);
          h_PtMDiMu[a][DiMu_b][DiMu_c][DiMu_d]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
          h_PtMDiMu[a][DiMu_b][DiMu_c][DiMu_d]->GetYaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");

          h_MDiMu[a][DiMu_b][DiMu_c][DiMu_d] = new TH1D(Form("h_MDiMu_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), Form("h_MDiMu_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), Mass_Bin[a], Mass_Edge[a][0], Mass_Edge[a][1]);
          h_MDiMu[a][DiMu_b][DiMu_c][DiMu_d]->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
          h_MDiMu[a][DiMu_b][DiMu_c][DiMu_d]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

          h_pdgDimuMu[a][DiMu_b][DiMu_c][DiMu_d] = new TH1D(Form("h_pdgDimuMu_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), Form("h_pdgDimuMu_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), 12000, -5999.5, 5999.5);
          h_pdgDimuMu[a][DiMu_b][DiMu_c][DiMu_d]->GetXaxis()->SetTitle("PDG Code");
          h_pdgDimuMu[a][DiMu_b][DiMu_c][DiMu_d]->GetYaxis()->SetTitle("Counts");

          h_nDiMu_xevent[a][DiMu_b][DiMu_c][DiMu_d] = new TH1D(Form("h_nDiMu_xevent_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), Form("h_nDiMu_xevent_%s_%s_%s_%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data(), name_DiMuSelection[DiMu_d].Data()), 10, -0.5, 9.5);
          h_nDiMu_xevent[a][DiMu_b][DiMu_c][DiMu_d]->GetXaxis()->SetTitle("Number of #mu x event");
          h_nDiMu_xevent[a][DiMu_b][DiMu_c][DiMu_d]->GetYaxis()->SetTitle("Counts");

        } // End definition over DiMu selection

      } // End definition over DiMu cut
    }
  }
}
