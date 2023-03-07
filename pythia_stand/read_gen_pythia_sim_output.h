
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

Int_t PDG_HFquark_gen[fMuons_dim]; // single gen c/cbar PDG mum
Int_t PDG_HFquark_gen_daughter1[fMuons_dim];
Int_t PDG_HFquark_gen_daughter2[fMuons_dim];
Double_t Pt_HFquark_gen[fMuons_dim]; // single gen c/cbar or b/bbar HFquark pT
Double_t Y_HFquark_gen[fMuons_dim];  // single gen c/cbar or b/bbar HFquark y

Int_t PDGmum_gen[fMuons_dim];    // single gen mu PDG mum
Double_t Pt_mum_gen[fMuons_dim]; // pt of single mu mum
Double_t Y_mum_gen[fMuons_dim];  // Y of single mu mum
Int_t PDG_gen[fMuons_dim];       // single gen mu PDG mum
Int_t Promptmum_gen[fMuons_dim]; // single gen mu PDG mum
Double_t Pt_gen[fMuons_dim];     // single gen mu pT
Double_t E_gen[fMuons_dim];      // single gen mu E
Double_t Px_gen[fMuons_dim];     // single gen mu px
Double_t Py_gen[fMuons_dim];     // single gen mu py
Double_t Pz_gen[fMuons_dim];     // single gen mu pz
Double_t Y_gen[fMuons_dim];      // single gen mu y
Double_t Eta_gen[fMuons_dim];    // single gen mu eta
Double_t Phi_gen[fMuons_dim];    // single gen mu phi
Double_t Theta_gen[fMuons_dim];  // single gen mu theta
Int_t Charge_gen[fMuons_dim];    // single gen mu theta

Int_t DimuMu_gen[fDimu_dim][2];   // reference to single gen mus
Double_t DimuPt_gen[fDimu_dim];   // gen dimuon pT
Double_t DimuPx_gen[fDimu_dim];   // gen dimuon px
Double_t DimuPy_gen[fDimu_dim];   // gen dimuon py
Double_t DimuPz_gen[fDimu_dim];   // gen dimuon pz
Double_t DimuY_gen[fDimu_dim];    // gen dimuon y
Double_t DimuMass_gen[fDimu_dim]; // gen dimuon invariant mass
Int_t DimuCharge_gen[fDimu_dim];  // gen dimuon charge

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

//-----------------------------------------------------//
// 0 For All, 1 For HF, 2 For Charm, 3 For Beauty, 4 For Charm Mesons, 5 For Charm Barions, 6 For Beauty Mesons, 7 For Beauty Barions
const Int_t n_MuSelection = 8;
TString name_MuSelection[n_MuSelection];

// 0 For All, 1 For HF, 2 For Charm, 3 For Beauty, 4 For HF Mixed (one muon from Charm, one muon from Beauty), 5 For Charm Mesons (two muons from Charm Mesons), 6 For Charm Barions (two muon from Charm Barions), 7 For Charm Mixed (one muon from Charm Mesons, one muon from Charm Barions), 8 For Beauty Mesons (two muons from Beauty Mesons), 9 For Beauty Barions (two muon from Beauty Barions), 10 For Beauty Mixed (one muon from Beauty Mesons, one muon from Beauty Barions).

const Int_t n_DiMuSelection = 11;
TString name_DiMuSelection[n_DiMuSelection];

//-----------------------------------------------------//

// 0 for Plus, 1 for Minus, 2 for total
const Int_t n_Mu_Charge = 3;
TString name_MuCharge[n_Mu_Charge];

// 0 for ULS, 1 for LS++, 2 for LS--, 3 LS total
const Int_t n_DiMu_Charge = 4;
TString name_DiMu_Charge[n_DiMu_Charge];

TH2D *h_PtYMu[n_MuSelection][n_Mu_Charge];
TH2D *h_PtMu_PtMum[n_MuSelection][n_Mu_Charge];
TH2D *h_YMu_YMum[n_MuSelection][n_Mu_Charge];
TH1D *h_PhiMu[n_MuSelection][n_Mu_Charge];
TH1D *h_pdgMu[n_MuSelection][n_Mu_Charge];
TH1D *h_nMu_xevent[n_MuSelection][n_Mu_Charge];

TH1D *h_pdgDimuMu[n_DiMuSelection][n_DiMu_Charge];
TH2D *h_PtYDiMu[n_DiMuSelection][n_DiMu_Charge];
TH2D *h_PtMDiMu[n_DiMuSelection][n_DiMu_Charge];
TH1D *h_MDiMu[n_DiMuSelection][n_DiMu_Charge];
TH1D *h_nDiMu_xevent[n_DiMuSelection][n_DiMu_Charge];

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

TH1D *h_ptHadron_prompt_forCScalc[n_Hadron_studied][n_Rapidity_studied];
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
        h_ptHadron_prompt_forCScalc[i_hadron][i_rapidity] = (TH1D *)hist_data->Clone(Form("h_pt%s_prompt_%s_forCScalc", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
        h_ptHadron_prompt_forCScalc[i_hadron][i_rapidity]->Reset();
        h_ptHadron_prompt_forCScalc[i_hadron][i_rapidity]->SetTitle(Form("h_pt%s_prompt_%s_forCScalc", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
      }
      else
      {
        if (i_hadron < n_Hadron_studied - 1)
        {
          TFile *fIn_data = new TFile(Form("%s", name_root_files[n_Hadron_studied].Data()), "READ");
          TH1D *hist_data = (TH1D *)fIn_data->Get(Form("Table %zu/Hist1D_y%zu", i_hadron + 1, i_rapidity));
          hist_data->SetDirectory(0);
          fIn_data->Close();
          h_ptHadron_prompt_forCScalc[i_hadron][i_rapidity] = (TH1D *)hist_data->Clone(Form("h_pt%s_prompt_%s_forCScalc", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
          h_ptHadron_prompt_forCScalc[i_hadron][i_rapidity]->Reset();
          h_ptHadron_prompt_forCScalc[i_hadron][i_rapidity]->SetTitle(Form("h_pt%s_prompt_%s_forCScalc", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
        }
        else
        {
          h_ptHadron_prompt_forCScalc[i_hadron][i_rapidity] = new TH1D(Form("h_pt%s_prompt_%s_forCScalc", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_pt%s_prompt_%s_forCScalc", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 300, 0, 300);
        }
      }

      h_ptHadron_prompt[i_hadron][i_rapidity] = new TH1D(Form("h_pt%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_pt%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 300, 0, 30);
      h_yHadron_prompt[i_hadron][i_rapidity] = new TH1D(Form("h_y%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_y%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 200, -10, 10);
      h_pdgHadron_prompt[i_hadron][i_rapidity] = new TH1D(Form("h_pdg%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_pdg%s_prompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 12000, -6000, 6000);

      h_ptHadron_Notprompt[i_hadron][i_rapidity] = new TH1D(Form("h_pt%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_pt%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 300, 0, 30);
      h_yHadron_Notprompt[i_hadron][i_rapidity] = new TH1D(Form("h_y%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_y%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 200, -10, 10);
      h_pdgHadron_Notprompt[i_hadron][i_rapidity] = new TH1D(Form("h_pdg%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), Form("h_pdg%s_Notprompt_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()), 12000, -6000, 6000);
    }
  }

  name_MuCharge[0].Form("All");
  name_MuCharge[1].Form("Plus");
  name_MuCharge[2].Form("Minus");

  name_MuSelection[0].Form("All");
  name_MuSelection[1].Form("fromHF");
  name_MuSelection[2].Form("fromCharm");
  name_MuSelection[3].Form("fromBeauty");
  name_MuSelection[4].Form("fromCharm_Mesons");
  name_MuSelection[5].Form("fromCharm_Barions");
  name_MuSelection[6].Form("fromBeauty_Mesons");
  name_MuSelection[7].Form("fromBeauty_Barions");

  name_DiMuSelection[0].Form("All");
  name_DiMuSelection[1].Form("fromHF");
  name_DiMuSelection[2].Form("fromCharm");
  name_DiMuSelection[3].Form("fromBeauty");
  name_DiMuSelection[4].Form("fromMixed");
  name_DiMuSelection[5].Form("fromCharm_Mesons");
  name_DiMuSelection[6].Form("fromCharm_Barions");
  name_DiMuSelection[7].Form("fromCharm_Mixed");
  name_DiMuSelection[8].Form("fromBeauty_Mesons");
  name_DiMuSelection[9].Form("fromBeauty_Barions");
  name_DiMuSelection[10].Form("fromBeauty_Mixed");

  name_DiMu_Charge[0].Form("ULS");
  name_DiMu_Charge[1].Form("LSplus");
  name_DiMu_Charge[2].Form("LSminus");
  name_DiMu_Charge[3].Form("LS");

  for (Int_t Mu_a = 0; Mu_a < n_MuSelection; Mu_a++)
  {
    for (Int_t Mu_b = 0; Mu_b < n_Mu_Charge; Mu_b++)
    {

      h_PtYMu[Mu_a][Mu_b] = new TH2D(Form("h_PtYMu_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), Form("h_PtYMu_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), 300, 0.0, 30.0, 150, -4.0, -2.5);
      h_PtYMu[Mu_a][Mu_b]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      h_PtYMu[Mu_a][Mu_b]->GetYaxis()->SetTitle("Y");

      h_PtMu_PtMum[Mu_a][Mu_b] = new TH2D(Form("h_PtMu_PtMum_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), Form("h_PtMu_PtMum_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), 300, 0.0, 30.0, 300, 0.0, 30.0);
      h_PtMu_PtMum[Mu_a][Mu_b]->GetXaxis()->SetTitle("#it{p}_{T} #mu (GeV/#it{c})");
      h_PtMu_PtMum[Mu_a][Mu_b]->GetYaxis()->SetTitle("#it{p}_{T} mum (GeV/#it{c})");
      
      h_YMu_YMum[Mu_a][Mu_b] = new TH2D(Form("h_YMu_YMum_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), Form("h_YMu_YMum_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), 150, -4.0, -2.5, 150, -4.0, -2.5);
      h_YMu_YMum[Mu_a][Mu_b]->GetXaxis()->SetTitle("Y #mu");
      h_YMu_YMum[Mu_a][Mu_b]->GetYaxis()->SetTitle("Y mum");

      h_PhiMu[Mu_a][Mu_b] = new TH1D(Form("h_PhiMu_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), Form("h_PhiMu_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), 314, -3.14, 3.14);
      h_PhiMu[Mu_a][Mu_b]->GetXaxis()->SetTitle("#phi");
      
      h_pdgMu[Mu_a][Mu_b] = new TH1D(Form("h_pdgMu_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), Form("h_pdgMu_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), 12000, -5999.5, 5999.5);
      h_pdgMu[Mu_a][Mu_b]->GetXaxis()->SetTitle("PDG Code");
      h_pdgMu[Mu_a][Mu_b]->GetYaxis()->SetTitle("Counts");

      h_nMu_xevent[Mu_a][Mu_b] = new TH1D(Form("h_nMu_xevent_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), Form("h_nMu_xevent_%s_%s", name_MuSelection[Mu_a].Data(), name_MuCharge[Mu_b].Data()), 10, -0.5, 9.5);
      h_nMu_xevent[Mu_a][Mu_b]->GetXaxis()->SetTitle("Number of #mu x event");
      h_nMu_xevent[Mu_a][Mu_b]->GetYaxis()->SetTitle("Counts");

    } // End definition over Mu charge

  } // End definition over Mu selection
  for (Int_t DiMu_a = 0; DiMu_a < n_DiMuSelection; DiMu_a++)
  {
    for (Int_t DiMu_b = 0; DiMu_b < n_DiMu_Charge; DiMu_b++)
    {

      h_PtYDiMu[DiMu_a][DiMu_b] = new TH2D(Form("h_PtYDiMu_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), Form("h_PtYDiMu_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), 300, 0.0, 30.0, 150, -4.0, -2.5);
      h_PtYDiMu[DiMu_a][DiMu_b]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      h_PtYDiMu[DiMu_a][DiMu_b]->GetYaxis()->SetTitle("Y");

      h_PtMDiMu[DiMu_a][DiMu_b] = new TH2D(Form("h_PtMDiMu_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), Form("h_PtMDiMu_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), 300, 0.0, 30.0, 300, 0.0, 30.0);
      h_PtMDiMu[DiMu_a][DiMu_b]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      h_PtMDiMu[DiMu_a][DiMu_b]->GetYaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");

      h_MDiMu[DiMu_a][DiMu_b] = new TH1D(Form("h_MDiMu_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), Form("h_MDiMu_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), 300, 0.0, 30.0);
      h_MDiMu[DiMu_a][DiMu_b]->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
      h_MDiMu[DiMu_a][DiMu_b]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

      h_pdgDimuMu[DiMu_a][DiMu_b] = new TH1D(Form("h_pdgDimuMu_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), Form("h_pdgDimuMu_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), 12000, -5999.5, 5999.5);
      h_pdgDimuMu[DiMu_a][DiMu_b]->GetXaxis()->SetTitle("PDG Code");
      h_pdgDimuMu[DiMu_a][DiMu_b]->GetYaxis()->SetTitle("Counts");

      h_nDiMu_xevent[DiMu_a][DiMu_b] = new TH1D(Form("h_nDiMu_xevent_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), Form("h_nDiMu_xevent_%s_%s", name_DiMuSelection[DiMu_a].Data(), name_DiMu_Charge[DiMu_b].Data()), 10, -0.5, 9.5);
      h_nDiMu_xevent[DiMu_a][DiMu_b]->GetXaxis()->SetTitle("Number of #mu x event");
      h_nDiMu_xevent[DiMu_a][DiMu_b]->GetYaxis()->SetTitle("Counts");

    } // End definition over DiMu selection

  } // End definition over DiMu cut
}
