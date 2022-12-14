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

Int_t N_HFquarks_gen; // gen c/cbar or b/bar HFquarks in the event
Int_t N_HFquarks_rec; // gen c/cbar or b/bar HFquarks in the event

Int_t NMuons_gen; // gen muon in the event
Int_t NDimu_gen;  // gen dimuons in the event
Int_t NMuons_rec; // rec muon tracks in the event
Int_t NDimu_rec;  // rec dimuons in the event

const Int_t fMuons_dim = 1000;
const Int_t fDimu_dim = 1000;

Int_t PDG_HFquark_rec[fMuons_dim];   // single rec c/cbar PDG mum
Double_t Pt_HFquark_rec[fMuons_dim]; // single rec c/cbar or b/bbar HFquark pT
Double_t Y_HFquark_rec[fMuons_dim];  // single rec c/cbar or b/bbar HFquark y

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

Int_t PDG_HFquark_gen[fMuons_dim];   // single gen c/cbar PDG mum
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

const Int_t n_MuSelection = 3;
TString name_MuSelection[n_MuSelection];
const Int_t n_Mu_Charge = 3;
TString name_Mu_Charge[n_Mu_Charge];

const Int_t n_DiMuSelection = 3;
TString name_DiMuSelection[n_DiMuSelection];
const Int_t n_DiMu_Charge = 4;
TString name_DiMu_Charge[n_DiMu_Charge];

const Int_t n_Fraction = 11;
TString name_n_Fraction[n_Fraction];

Double_t Fraction[n_Fraction] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};

TH2D *h_PtYMu_ccut[n_Mu_Charge][n_MuSelection][n_Fraction];
TH1D *h_pdgMu_ccut[n_Mu_Charge][n_MuSelection][n_Fraction];
TH1D *h_nMu_xevent_ccut[n_Mu_Charge][n_MuSelection][n_Fraction];

TH2D *h_PtYMu_bcut[n_Mu_Charge][n_MuSelection][n_Fraction];
TH1D *h_pdgMu_bcut[n_Mu_Charge][n_MuSelection][n_Fraction];
TH1D *h_nMu_xevent_bcut[n_Mu_Charge][n_MuSelection][n_Fraction];

TH2D *h_PtYDiMu_ccut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];
TH2D *h_PtMDiMu_ccut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];
TH1D *h_MDiMu_ccut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];
TH1D *h_pdgDimuMu_ccut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];
TH1D *h_nDiMu_xevent_ccut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];

TH2D *h_PtYDiMu_bcut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];
TH2D *h_PtMDiMu_bcut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];
TH1D *h_MDiMu_bcut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];
TH1D *h_pdgDimuMu_bcut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];
TH1D *h_nDiMu_xevent_bcut[n_DiMu_Charge][n_DiMuSelection][n_Fraction];

void SetHist()
{
    name_MuSelection[0].Form("fromCharm");
    name_MuSelection[1].Form("fromBeauty");
    name_MuSelection[2].Form("fromHF");

    name_Mu_Charge[0].Form("plus");
    name_Mu_Charge[1].Form("minus");
    name_Mu_Charge[2].Form("total");

    for (Int_t Mu_a = 0; Mu_a < n_Mu_Charge; Mu_a++)
    {
        for (Int_t Mu_b = 0; Mu_b < n_MuSelection; Mu_b++)
        {
            for (Int_t Mu_c = 0; Mu_c < n_Fraction; Mu_c++)
            {
                h_PtYMu_ccut[Mu_a][Mu_b][Mu_c] = new TH2D(Form("h_PtYMu_ccut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), Form("h_PtYMu_ccut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), 300, 0.0, 30.0, 200, -10.0, 10.0);
                h_PtYMu_ccut[Mu_a][Mu_b][Mu_c]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                h_PtYMu_ccut[Mu_a][Mu_b][Mu_c]->GetYaxis()->SetTitle("Y");

                h_pdgMu_ccut[Mu_a][Mu_b][Mu_c] = new TH1D(Form("h_pdgMu_ccut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), Form("h_pdgMu_ccut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), 12000, -5999.5, 5999.5);
                h_pdgMu_ccut[Mu_a][Mu_b][Mu_c]->GetXaxis()->SetTitle("PDG Code");
                h_pdgMu_ccut[Mu_a][Mu_b][Mu_c]->GetYaxis()->SetTitle("Counts");

                h_nMu_xevent_ccut[Mu_a][Mu_b][Mu_c] = new TH1D(Form("h_nMu_xevent_ccut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), Form("h_nMu_xevent_ccut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), 10, -0.5, 9.5);
                h_nMu_xevent_ccut[Mu_a][Mu_b][Mu_c]->GetXaxis()->SetTitle("Number of #mu x event");
                h_nMu_xevent_ccut[Mu_a][Mu_b][Mu_c]->GetYaxis()->SetTitle("Counts");

                h_PtYMu_bcut[Mu_a][Mu_b][Mu_c] = new TH2D(Form("h_PtYMu_bcut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), Form("h_PtYMu_bcut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), 300, 0.0, 30.0, 200, -10.0, 10.0);
                h_PtYMu_bcut[Mu_a][Mu_b][Mu_c]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                h_PtYMu_bcut[Mu_a][Mu_b][Mu_c]->GetYaxis()->SetTitle("Y");

                h_pdgMu_bcut[Mu_a][Mu_b][Mu_c] = new TH1D(Form("h_pdgMu_bcut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), Form("h_pdgMu_bcut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), 12000, -5999.5, 5999.5);
                h_pdgMu_bcut[Mu_a][Mu_b][Mu_c]->GetXaxis()->SetTitle("PDG Code");
                h_pdgMu_bcut[Mu_a][Mu_b][Mu_c]->GetYaxis()->SetTitle("Counts");

                h_nMu_xevent_bcut[Mu_a][Mu_b][Mu_c] = new TH1D(Form("h_nMu_xevent_bcut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), Form("h_nMu_xevent_bcut_%s_%s_fraction%0.2f", name_Mu_Charge[Mu_a].Data(), name_MuSelection[Mu_b].Data(), Fraction[Mu_c]), 10, -0.5, 9.5);
                h_nMu_xevent_bcut[Mu_a][Mu_b][Mu_c]->GetXaxis()->SetTitle("Number of #mu x event");
                h_nMu_xevent_bcut[Mu_a][Mu_b][Mu_c]->GetYaxis()->SetTitle("Counts");
            }
        }
    }

    name_DiMuSelection[0].Form("fromCharm");
    name_DiMuSelection[1].Form("fromBeauty");
    name_DiMuSelection[2].Form("fromMixed");
    // name_DiMuSelection[3].Form("fromHF");

    name_DiMu_Charge[0].Form("ULS");
    name_DiMu_Charge[1].Form("LSplus");
    name_DiMu_Charge[2].Form("LSminus");
    name_DiMu_Charge[3].Form("LS");

    for (Int_t DiMu_a = 0; DiMu_a < n_DiMu_Charge; DiMu_a++)
    {
        for (Int_t DiMu_b = 0; DiMu_b < n_DiMuSelection; DiMu_b++)
        {
            for (Int_t DiMu_c = 0; DiMu_c < n_Fraction; DiMu_c++)
            {
                h_PtYDiMu_ccut[DiMu_a][DiMu_b][DiMu_c] = new TH2D(Form("h_PtYDiMu_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_PtYDiMu_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 300, 0.0, 30.0, 200, -10.0, 10.0);
                h_PtYDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                h_PtYDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("Y");

                h_PtMDiMu_ccut[DiMu_a][DiMu_b][DiMu_c] = new TH2D(Form("h_PtMDiMu_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_PtMDiMu_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 300, 0.0, 30.0, 260, 4.0, 30.0);
                h_PtMDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                h_PtMDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");

                h_MDiMu_ccut[DiMu_a][DiMu_b][DiMu_c] = new TH1D(Form("h_MDiMu_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_MDiMu_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 260, 4.0, 30.0);
                h_MDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
                h_MDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

                h_pdgDimuMu_ccut[DiMu_a][DiMu_b][DiMu_c] = new TH1D(Form("h_pdgDimuMu_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_pdgDimuMu_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 12000, -5999.5, 5999.5);
                h_pdgDimuMu_ccut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("PDG Code");
                h_pdgDimuMu_ccut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("Counts");

                h_nDiMu_xevent_ccut[DiMu_a][DiMu_b][DiMu_c] = new TH1D(Form("h_nDiMu_xevent_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_nDiMu_xevent_ccut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 10, -0.5, 9.5);
                h_nDiMu_xevent_ccut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("Number of #mu x event");
                h_nDiMu_xevent_ccut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("Counts");

                h_PtYDiMu_bcut[DiMu_a][DiMu_b][DiMu_c] = new TH2D(Form("h_PtYDiMu_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_PtYDiMu_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 300, 0.0, 30.0, 200, -10.0, 10.0);
                h_PtYDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                h_PtYDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("Y");

                h_PtMDiMu_bcut[DiMu_a][DiMu_b][DiMu_c] = new TH2D(Form("h_PtMDiMu_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_PtMDiMu_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 300, 0.0, 30.0, 260, 4.0, 30.0);
                h_PtMDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                h_PtMDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");

                h_MDiMu_bcut[DiMu_a][DiMu_b][DiMu_c] = new TH1D(Form("h_MDiMu_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_MDiMu_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 260, 4.0, 30.0);
                h_MDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
                h_MDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

                h_pdgDimuMu_bcut[DiMu_a][DiMu_b][DiMu_c] = new TH1D(Form("h_pdgDimuMu_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_pdgDimuMu_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 12000, -5999.5, 5999.5);
                h_pdgDimuMu_bcut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("PDG Code");
                h_pdgDimuMu_bcut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("Counts");

                h_nDiMu_xevent_bcut[DiMu_a][DiMu_b][DiMu_c] = new TH1D(Form("h_nDiMu_xevent_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), Form("h_nDiMu_xevent_bcut_%s_%s_fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), name_DiMuSelection[DiMu_b].Data(), Fraction[DiMu_c]), 10, -0.5, 9.5);
                h_nDiMu_xevent_bcut[DiMu_a][DiMu_b][DiMu_c]->GetXaxis()->SetTitle("Number of #mu x event");
                h_nDiMu_xevent_bcut[DiMu_a][DiMu_b][DiMu_c]->GetYaxis()->SetTitle("Counts");
            }
        } // End definition over DiMu selection
    }
}