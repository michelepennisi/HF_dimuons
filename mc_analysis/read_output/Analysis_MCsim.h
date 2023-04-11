#include "/home/michele_pennisi/cernbox/common_include.h"

void axis_name(TH3F *histo)
{
    histo->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    histo->GetYaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");

    TString name_axis[3];
    name_axis[0].Form("Charm");
    name_axis[1].Form("Beauty");
    name_axis[2].Form("Mixed");

    for (Int_t i_bin = 0; i_bin < 3; i_bin++)
    {
        histo->GetZaxis()->SetBinLabel(i_bin + 1, Form("%s", name_axis[i_bin].Data()));
    }
    return;
}

Int_t N_HFquarks_gen; // gen c/cbar or b/bar HFquarks in the event
// Int_t N_HFquarks_rec; // gen c/cbar or b/bar HFquarks in the event

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

// Pt-PDG for Generated Muons +
TH2F *h_PtPdg_Muon_Gen_Meson = new TH2F("h_PtPdg_Muon_Gen_Meson", "h_PtPdg_Muon_Gen_Meson", 300, 0.0, 30.0, 200, 400, 600);
TH2F *h_PtPdg_Muon_Gen_Barion = new TH2F("h_PtPdg_Muon_Gen_Barion", "h_PtPdg_Muon_Gen_Barion", 300, 0.0, 30.0, 2000, 4000, 6000);

// Pt-Pdg for Reconstructed Muons +
TH2F *h_PtPdg_Muon_Rec_Meson = new TH2F("h_PtPdg_Muon_Rec_Meson", "h_PtPdg_Muon_Rec_Meson", 300, 0.0, 30.0, 200, 400, 600);
TH2F *h_PtPdg_Muon_Rec_Barion = new TH2F("h_PtPdg_Muon_Rec_Barion", "h_PtPdg_Muon_Rec_Barion", 300, 0.0, 30.0, 2000, 4000, 6000);

// Pt-Pdg for Corrected Muons +
TH2F *h_PtPdg_Muon_Corr_Meson = new TH2F("h_PtPdg_Muon_Corr_Meson", "h_PtPdg_Muon_Corr_Meson", 300, 0.0, 30.0, 200, 400, 600);
TH2F *h_PtPdg_Muon_Corr_Barion = new TH2F("h_PtPdg_Muon_Corr_Barion", "h_PtPdg_Muon_Corr_Barion", 300, 0.0, 30.0, 2000, 4000, 6000);

const Int_t n_Muon_origin = 3;

TH3F *h_PtYPdg_Muon_Gen[n_Muon_origin];
TH3F *h_PtYPdg_Muon_Rec[n_Muon_origin];

TString Muon_origin[n_Muon_origin];

const Int_t n_DiMu_origin = 5;
// Generated Dimuons
TH3F *h_PtMPdg_DiMu_Gen_ULS[n_DiMu_origin];
TH3F *h_PtYPdg_DiMu_Gen_ULS[n_DiMu_origin];
TH3F *h_PtYPdg_DiMu_Gen_ULS_M49_Pt010[n_DiMu_origin];
TH2F *h_PtPt_DiMuULS_Muon_Gen[n_DiMu_origin];
TH2F *h_PtPt_DiMuULS_M49_Pt010_Muon_Gen[n_DiMu_origin];
// Reconstructed Dimuons
TH3F *h_PtMPdg_DiMu_Rec_ULS[n_DiMu_origin];
TH3F *h_PtYPdg_DiMu_Rec_ULS[n_DiMu_origin];
TH3F *h_PtYPdg_DiMu_Rec_ULS_M49_Pt010[n_DiMu_origin];
TH2F *h_PtPt_DiMuULS_Muon_Rec[n_DiMu_origin];
TH2F *h_PtPt_DiMuULS_M49_Pt010_Muon_Rec[n_DiMu_origin];
// Corrected Dimuons
TH3F *h_PtMPdg_DiMu_DiMu_Corr_ULS[n_DiMu_origin];
TH3F *h_PtYPdg_DiMu_DiMu_Corr_ULS[n_DiMu_origin];
TH3F *h_PtYPdg_DiMu_DiMu_Corr_ULS_M49_Pt010[n_DiMu_origin];

TH3F *h_PtMPdg_DiMu_Muon_Corr_ULS[n_DiMu_origin];
TH3F *h_PtYPdg_DiMu_Muon_Corr_ULS[n_DiMu_origin];
TH3F *h_PtYPdg_DiMu_Muon_Corr_ULS_M49_Pt010[n_DiMu_origin];

TString DiMu_origin[n_DiMu_origin];

void Set_Hist()
{
    Muon_origin[0].Form("Charm");
    Muon_origin[1].Form("Beauty");
    Muon_origin[2].Form("LF");

    for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
    {
        h_PtYPdg_Muon_Gen[i_Muon_origin] = new TH3F(Form("h_PtYPdg_Muon_Gen%s", Muon_origin[i_Muon_origin].Data()), Form("h_PtYPdg_Muon_Gen%s", Muon_origin[i_Muon_origin].Data()), 300, 0.0, 30.0, 150, -4.0, -2.5, 3, 0, 3);
        h_PtYPdg_Muon_Rec[i_Muon_origin] = new TH3F(Form("h_PtYPdg_Muon_Rec%s", Muon_origin[i_Muon_origin].Data()), Form("h_PtYPdg_Muon_Rec%s", Muon_origin[i_Muon_origin].Data()), 300, 0.0, 30.0, 150, -4.0, -2.5, 3, 0, 3);
        ;
    }

    DiMu_origin[0].Form("Charm");
    DiMu_origin[1].Form("Beauty");
    DiMu_origin[2].Form("HF_Mixed");
    DiMu_origin[3].Form("LF");
    DiMu_origin[4].Form("LF_HF_Mixed");

    for (Int_t i_DiMu_origin = 0; i_DiMu_origin < n_DiMu_origin; i_DiMu_origin++)
    {
        h_PtMPdg_DiMu_Gen_ULS[i_DiMu_origin] = new TH3F(Form("h_PtMPdg_DiMu_Gen_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtMPdg_DiMu_Gen_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 0, 3);
        h_PtYPdg_DiMu_Gen_ULS[i_DiMu_origin] = new TH3F(Form("h_PtYPdg_DiMu_Gen_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtYPdg_DiMu_Gen_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0.0, 30.0, 150, -4.0, -2.5, 3, 0, 3);
        h_PtYPdg_DiMu_Gen_ULS_M49_Pt010[i_DiMu_origin] = new TH3F(Form("h_PtYPdg_DiMu_Gen_ULS_M49_Pt010_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtYPdg_DiMu_Gen_ULS_M49_Pt010_%s", DiMu_origin[i_DiMu_origin].Data()), 100, 0.0, 10.0, 150, -4.0, -2.5, 3, 0, 3);
        h_PtPt_DiMuULS_Muon_Gen[i_DiMu_origin] = new TH2F(Form("h_PtPt_DiMuULS_Muon_Gen_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtPt_DiMuULS_Muon_Gen_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0, 30.0, 300, 0, 30.0);
        h_PtPt_DiMuULS_M49_Pt010_Muon_Gen[i_DiMu_origin] = new TH2F(Form("h_PtPt_DiMuULS_M49_Pt010_Gen_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtPt_DiMuULS_M49_Pt010_Gen_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0, 30.0, 300, 0, 30.0);

        h_PtMPdg_DiMu_Rec_ULS[i_DiMu_origin] = new TH3F(Form("h_PtMPdg_DiMu_Rec_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtMPdg_DiMu_Rec_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 0, 3);
        h_PtYPdg_DiMu_Rec_ULS[i_DiMu_origin] = new TH3F(Form("h_PtYPdg_DiMu_Rec_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtYPdg_DiMu_Rec_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0.0, 30.0, 150, -4.0, -2.5, 3, 0, 3);
        h_PtYPdg_DiMu_Rec_ULS_M49_Pt010[i_DiMu_origin] = new TH3F(Form("h_PtYPdg_DiMu_Rec_ULS_M49_Pt010_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtYPdg_DiMu_Rec_ULS_M49_Pt010_%s", DiMu_origin[i_DiMu_origin].Data()), 100, 0.0, 10.0, 150, -4.0, -2.5, 3, 0, 3);
        h_PtPt_DiMuULS_Muon_Rec[i_DiMu_origin] = new TH2F(Form("h_PtPt_DiMuULS_Muon_Rec_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtPt_DiMuULS_Muon_Rec_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0, 30.0, 300, 0, 30.0);
        h_PtPt_DiMuULS_M49_Pt010_Muon_Rec[i_DiMu_origin] = new TH2F(Form("h_PtPt_DiMuULS_M49_Pt010_Rec_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtPt_DiMuULS_M49_Pt010_Rec_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0, 30.0, 300, 0, 30.0);

        h_PtMPdg_DiMu_DiMu_Corr_ULS[i_DiMu_origin] = new TH3F(Form("h_PtMPdg_DiMu_Corr_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtMPdg_DiMu_Corr_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 0, 3);
        h_PtYPdg_DiMu_DiMu_Corr_ULS[i_DiMu_origin] = new TH3F(Form("h_PtYPdg_DiMu_Corr_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtYPdg_DiMu_Corr_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0.0, 30.0, 150, -4.0, -2.5, 3, 0, 3);
        h_PtYPdg_DiMu_DiMu_Corr_ULS_M49_Pt010[i_DiMu_origin] = new TH3F(Form("h_PtYPdg_DiMu_Corr_ULS_M49_Pt010_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtYPdg_DiMu_Corr_ULS_M49_Pt010_%s", DiMu_origin[i_DiMu_origin].Data()), 100, 0.0, 10.0, 150, -4.0, -2.5, 3, 0, 3);

        h_PtMPdg_DiMu_Muon_Corr_ULS[i_DiMu_origin] = new TH3F(Form("h_PtMPdg_Muon_Corr_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtMPdg_Muon_Corr_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 0, 3);
        h_PtYPdg_DiMu_Muon_Corr_ULS[i_DiMu_origin] = new TH3F(Form("h_PtYPdg_Muon_Corr_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtYPdg_Muon_Corr_ULS_%s", DiMu_origin[i_DiMu_origin].Data()), 300, 0.0, 30.0, 150, -4.0, -2.5, 3, 0, 3);
        h_PtYPdg_DiMu_Muon_Corr_ULS_M49_Pt010[i_DiMu_origin] = new TH3F(Form("h_PtYPdg_Muon_Corr_ULS_M49_Pt010_%s", DiMu_origin[i_DiMu_origin].Data()), Form("h_PtYPdg_Muon_Corr_ULS_M49_Pt010_%s", DiMu_origin[i_DiMu_origin].Data()), 100, 0.0, 10.0, 150, -4.0, -2.5, 3, 0, 3);
    }
}

// TH3F *h_PtMPdg_DiMu_Gen_Meson_LSsign = new TH3F("h_PtMPdg_DiMu_Gen_Meson_LSsign", "h_PtMPdg_DiMu_Gen_Meson_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);
// TH3F *h_PtMPdg_DiMu_Rec_Meson_LSsign = new TH3F("h_PtMPdg_DiMu_Rec_Meson_LSsign", "h_PtMPdg_DiMu_Rec_Meson_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);
// TH3F *h_PtMPdg_DiMu_Gen_Meson_LSsign = new TH3F("h_PtMPdg_DiMu_Gen_Meson_LSsign", "h_PtMPdg_DiMu_Gen_Meson_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);
// TH3F *h_PtMPdg_DiMu_Rec_Meson_LSsign = new TH3F("h_PtMPdg_DiMu_Rec_Meson_LSsign", "h_PtMPdg_DiMu_Rec_Meson_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);

// TH3F *h_PtMPdg_DiMu_Gen_Barion_LSsign = new TH3F("h_PtMPdg_DiMu_Gen_Barion_LSsign", "h_PtMPdg_DiMu_Gen_Barion_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);
// TH3F *h_PtMPdg_DiMu_Rec_Barion_LSsign = new TH3F("h_PtMPdg_DiMu_Rec_Barion_LSsign", "h_PtMPdg_DiMu_Rec_Barion_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);
// TH3F *h_PtMPdg_DiMu_Gen_Barion_LSsign = new TH3F("h_PtMPdg_DiMu_Gen_Barion_LSsign", "h_PtMPdg_DiMu_Gen_Barion_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);
// TH3F *h_PtMPdg_DiMu_Rec_Barion_LSsign = new TH3F("h_PtMPdg_DiMu_Rec_Barion_LSsign", "h_PtMPdg_DiMu_Rec_Barion_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);

TChain *Importing_Tree(TString dir_filename, TString filename)
{
    TChain *tree = nullptr;

    tree = new TChain("MCTree");
    //    printf("%s \n",Form("%s/%s",dir_filename.Data(),filename.Data()));
    tree->AddFile(Form("%s/%s", dir_filename.Data(), filename.Data()));

    tree->SetBranchAddress("N_HFquarks_gen", &N_HFquarks_gen);

    tree->SetBranchAddress("PDG_HFquark_gen", PDG_HFquark_gen);
    tree->SetBranchAddress("Pt_HFquark_gen", Pt_HFquark_gen);
    tree->SetBranchAddress("Y_HFquark_gen", Y_HFquark_gen);

    tree->SetBranchAddress("NDimu_gen", &NDimu_gen);
    tree->SetBranchAddress("DimuMu_gen", DimuMu_gen);
    tree->SetBranchAddress("DimuPt_gen", DimuPt_gen);
    tree->SetBranchAddress("DimuPx_gen", DimuPx_gen);
    tree->SetBranchAddress("DimuPy_gen", DimuPy_gen);
    tree->SetBranchAddress("DimuPz_gen", DimuPz_gen);
    tree->SetBranchAddress("DimuY_gen", DimuY_gen);
    tree->SetBranchAddress("DimuMass_gen", DimuMass_gen);
    tree->SetBranchAddress("DimuCharge_gen", DimuCharge_gen);

    tree->SetBranchAddress("NMuons_gen", &NMuons_gen);
    tree->SetBranchAddress("PDGmum_gen", PDGmum_gen);
    tree->SetBranchAddress("Pt_gen", Pt_gen);
    tree->SetBranchAddress("E_gen", E_gen);
    tree->SetBranchAddress("Px_gen", Px_gen);
    tree->SetBranchAddress("Py_gen", Py_gen);
    tree->SetBranchAddress("Pz_gen", Pz_gen);
    tree->SetBranchAddress("Y_gen", Y_gen);
    tree->SetBranchAddress("Eta_gen", Eta_gen);
    tree->SetBranchAddress("Phi_gen", Phi_gen);
    tree->SetBranchAddress("Theta_gen", Theta_gen);
    tree->SetBranchAddress("Charge_gen", Charge_gen);

    tree->SetBranchAddress("NDimu_rec", &NDimu_rec);
    tree->SetBranchAddress("DimuMu_rec", DimuMu_rec);
    tree->SetBranchAddress("DimuPt_rec", DimuPt_rec);
    tree->SetBranchAddress("DimuPx_rec", DimuPx_rec);
    tree->SetBranchAddress("DimuPy_rec", DimuPy_rec);
    tree->SetBranchAddress("DimuPz_rec", DimuPz_rec);
    tree->SetBranchAddress("DimuY_rec", DimuY_rec);
    tree->SetBranchAddress("DimuMass_rec", DimuMass_rec);
    tree->SetBranchAddress("DimuCharge_rec", DimuCharge_rec);
    tree->SetBranchAddress("DimuMatch_rec", DimuMatch_rec);
    tree->SetBranchAddress("DimuPhi_rec", DimuPhi_rec);
    tree->SetBranchAddress("DimuTheta_rec", DimuTheta_rec);

    tree->SetBranchAddress("NMuons_rec", &NMuons_rec);
    tree->SetBranchAddress("PDGmum_rec", PDGmum_rec);
    tree->SetBranchAddress("E_rec", E_rec);
    tree->SetBranchAddress("Px_rec", Px_rec);
    tree->SetBranchAddress("Pt_rec", Pt_rec);
    tree->SetBranchAddress("Py_rec", Py_rec);
    tree->SetBranchAddress("Pz_rec", Pz_rec);
    tree->SetBranchAddress("Y_rec", Y_rec);
    tree->SetBranchAddress("Eta_rec", Eta_rec);
    tree->SetBranchAddress("MatchTrig_rec", MatchTrig_rec);
    tree->SetBranchAddress("TrackChi2_rec", TrackChi2_rec);
    tree->SetBranchAddress("MatchTrigChi2_rec", MatchTrigChi2_rec);
    tree->SetBranchAddress("Charge_rec", Charge_rec);
    tree->SetBranchAddress("RAtAbsEnd_rec", RAtAbsEnd_rec);
    tree->SetBranchAddress("pDCA_rec", pDCA_rec);
    tree->SetBranchAddress("Phi_rec", Phi_rec);
    tree->SetBranchAddress("Theta_rec", Theta_rec);

    return tree;
}

void progress_status(Int_t i_Event, Int_t total_entries)
{
    Double_t progress = (Double_t)i_Event / total_entries;
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " % (" << i_Event << "/ " << total_entries << ")\r" << std::scientific;
    std::cout.flush();

    std::cout << std::endl;
}