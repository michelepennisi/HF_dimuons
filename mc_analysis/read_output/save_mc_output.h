#include "/home/michele_pennisi/cernbox/common_include.h"
Int_t N_HFquarks_gen; // gen c/cbar or b/bar HFquarks in the event

Int_t NMuons_gen; // gen muon in the event
Int_t NDimu_gen;  // gen dimuons in the event
Int_t NMuons_rec; // rec muon tracks in the event
Int_t NDimu_rec;  // rec dimuons in the event

const Int_t fMuons_dim = 1000;
const Int_t fDimu_dim = 1000;

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
Int_t Charge_gen[fMuons_dim];   // single gen mu theta

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

TH1D *h_Nevents = new TH1D("h_Nevents", "h_Nevents", 2, 0, 2);

const Int_t n_Muon_origin = 5; // 0 for All sources, 1 for Charm, 2 for Beauty, 3 for LF, 4 for DY

//------Declaration Hist for Generated Muons----------------//
TH2F *h_PtPdg_Muon_Gen[n_Muon_origin];
TH2F *h_YPdg_Muon_Gen[n_Muon_origin];
TH1F *h_Nperevent_Muon_Gen[n_Muon_origin];
//------Declaration Hist for Reconstructed Muons----------------//
TH2F *h_PtPdg_Muon_Rec[n_Muon_origin];
TH2F *h_YPdg_Muon_Rec[n_Muon_origin];
TH1F *h_Nperevent_Muon_Rec[n_Muon_origin];

TString Muon_origin[n_Muon_origin];

const Int_t n_DiMuon_origin = 6; // 0 for Charm, 1 for Beauty, 2 for HF Mixed, 3 for LF, 4 for LF-HF Mixed, 5 for DY

//------Declaration Hist for Generated Dimuons----------------//
TH3F *h_PtPdg1Pdg2_DiMuon_Gen[n_DiMuon_origin];
TH3F *h_YPdg1Pdg2_DiMuon_Gen[n_DiMuon_origin];
TH3F *h_MPdg1Pdg2_DiMuon_Gen[n_DiMuon_origin];
TH2F *h_PtM_DiMuon_Gen[n_DiMuon_origin];
TH1F *h_Nperevent_DiMuon_Gen[n_DiMuon_origin];

//------Declaration Hist for Reconstructed Dimuons----------------//
TH3F *h_PtPdg1Pdg2_DiMuon_Rec[n_DiMuon_origin];
TH3F *h_YPdg1Pdg2_DiMuon_Rec[n_DiMuon_origin];
TH3F *h_MPdg1Pdg2_DiMuon_Rec[n_DiMuon_origin];
TH2F *h_PtM_DiMuon_Rec[n_DiMuon_origin];
TH1F *h_Nperevent_DiMuon_Rec[n_DiMuon_origin];

TString DiMuon_origin[n_DiMuon_origin];

Double_t Charm_Hadron[8] = {400, 410, 420, 430, 440, 4100, 4200, 4300};
void Set_Histograms()
{
    // 0 for All sources, 1 for Charm, 2 for Beauty, 3 for LF, 4 for DY
    Muon_origin[0].Form("All");
    Muon_origin[1].Form("Charm");
    Muon_origin[2].Form("Beauty");
    Muon_origin[3].Form("LF");
    Muon_origin[4].Form("DY");

    for (size_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
    {
        h_PtPdg_Muon_Gen[i_Muon_origin] = new TH2F(Form("h_PtPdg_Muon_Gen_%s", Muon_origin[i_Muon_origin].Data()), "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 300, 0, 30.0, 12000, -6000, 6000);
        h_YPdg_Muon_Gen[i_Muon_origin] = new TH2F(Form("h_YPdg_Muon_Gen_%s", Muon_origin[i_Muon_origin].Data()), "; #it{y} ; PDG code", 120, -6, 6, 12000, -6000, 6000);
        h_Nperevent_Muon_Gen[i_Muon_origin] = new TH1F(Form("h_Nperevent_Muon_Gen_%s", Muon_origin[i_Muon_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);

        h_PtPdg_Muon_Rec[i_Muon_origin] = new TH2F(Form("h_PtPdg_Muon_Rec_%s", Muon_origin[i_Muon_origin].Data()), "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 300, 0, 30.0, 12000, -6000, 6000);
        h_YPdg_Muon_Rec[i_Muon_origin] = new TH2F(Form("h_YPdg_Muon_Rec_%s", Muon_origin[i_Muon_origin].Data()), "; #it{y} ; PDG code", 150, -4.0, -2.5, 12000, -6000, 6000);
        h_Nperevent_Muon_Rec[i_Muon_origin] = new TH1F(Form("h_Nperevent_Muon_Rec_%s", Muon_origin[i_Muon_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);
    }
    Int_t n_bin_pt = h_PtPdg_Muon_Gen[0]->GetXaxis()->GetNbins();
    Double_t low_bin_pt[n_bin_pt+1];
    for (size_t i_lowpt = 0; i_lowpt <= (size_t)n_bin_pt; i_lowpt++)
        low_bin_pt[i_lowpt] = h_PtPdg_Muon_Gen[0]->GetXaxis()->GetBinLowEdge(i_lowpt);
    
    Int_t n_bin_Y_Gen = h_YPdg_Muon_Gen[0]->GetXaxis()->GetNbins();
    Double_t low_bin_Y_Gen[n_bin_Y_Gen+1];
    for (size_t i_lowY_Gen = 0; i_lowY_Gen <= (size_t)n_bin_Y_Gen; i_lowY_Gen++)
        low_bin_Y_Gen[i_lowY_Gen] = h_YPdg_Muon_Gen[0]->GetXaxis()->GetBinLowEdge(i_lowY_Gen);

    Int_t n_bin_Y_Rec = h_YPdg_Muon_Rec[0]->GetXaxis()->GetNbins();
    Double_t low_bin_Y_Rec[n_bin_Y_Rec+1];
    for (size_t i_lowY_Rec = 0; i_lowY_Rec <= (size_t)n_bin_Y_Rec; i_lowY_Rec++)
        low_bin_Y_Rec[i_lowY_Rec] = h_YPdg_Muon_Rec[0]->GetXaxis()->GetBinLowEdge(i_lowY_Rec);

    // 0 for Charm, 1 for Beauty, 2 for HF Mixed, 3 for LF, 4 for LF-HF Mixed, 5 for DY
    DiMuon_origin[0].Form("Charm");
    DiMuon_origin[1].Form("Beauty");
    DiMuon_origin[2].Form("HF_Mixed");
    DiMuon_origin[3].Form("LF");
    DiMuon_origin[4].Form("LF_HF_Mixed");
    DiMuon_origin[5].Form("DY");
    for (size_t i_DiMuon_origin = 0; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
    {
        h_PtPdg1Pdg2_DiMuon_Gen[i_DiMuon_origin] = new TH3F(Form("h_PtPdg1Pdg2_DiMuon_Gen_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; PDG code mum 1; PDG code mum 2; #it{p}_{T} (GeV/#it{c}) ", 7, Charm_Hadron, 7, Charm_Hadron, n_bin_pt, low_bin_pt);
        h_YPdg1Pdg2_DiMuon_Gen[i_DiMuon_origin] = new TH3F(Form("h_YPdg1Pdg2_DiMuon_Gen_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; PDG code mum 1; PDG code mum 2; #it{y}) ", 7, Charm_Hadron, 7, Charm_Hadron, n_bin_Y_Gen, low_bin_Y_Gen);
        h_MPdg1Pdg2_DiMuon_Gen[i_DiMuon_origin] = new TH3F(Form("h_MPdg1Pdg2_DiMuon_Gen_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; PDG code mum 1; PDG code mum 2; #it{m}_{#mu#mu} (GeV/#it{c}^{2}) ", 7, Charm_Hadron, 7, Charm_Hadron, n_bin_pt, low_bin_pt);
        h_PtM_DiMuon_Gen[i_DiMuon_origin] = new TH2F(Form("h_PtM_DiMuon_Gen_%s", DiMuon_origin[i_DiMuon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 300, 0, 30.0, 260, 4, 300);
        h_Nperevent_DiMuon_Gen[i_DiMuon_origin] = new TH1F(Form("h_Nperevent_DiMuon_Gen_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);

        h_PtPdg1Pdg2_DiMuon_Rec[i_DiMuon_origin] = new TH3F(Form("h_PtPdg1Pdg2_DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; PDG code mum 1; PDG code mum 2; #it{p}_{T} (GeV/#it{c}) ", 7, Charm_Hadron, 7, Charm_Hadron, n_bin_pt, low_bin_pt);
        h_YPdg1Pdg2_DiMuon_Rec[i_DiMuon_origin] = new TH3F(Form("h_YPdg1Pdg2_DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; PDG code mum 1; PDG code mum 2; #it{y}) ", 7, Charm_Hadron, 7, Charm_Hadron, n_bin_Y_Rec, low_bin_Y_Rec);
        h_MPdg1Pdg2_DiMuon_Rec[i_DiMuon_origin] = new TH3F(Form("h_MPdg1Pdg2_DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; PDG code mum 1; PDG code mum 2; #it{m}_{#mu#mu} (GeV/#it{c}^{2}) ", 7, Charm_Hadron, 7, Charm_Hadron, n_bin_pt, low_bin_pt);
        h_PtM_DiMuon_Rec[i_DiMuon_origin] = new TH2F(Form("h_PtM_DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 300, 0, 30.0, 260, 4, 300);
        h_Nperevent_DiMuon_Rec[i_DiMuon_origin] = new TH1F(Form("h_Nperevent_DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);
    }
}