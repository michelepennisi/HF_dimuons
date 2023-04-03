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
Int_t Charge_gen[fMuons_dim]; // single gen mu theta

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

//Pt-Y for Generated Muons +
TH3F *h_PtYPdg_Muon_Gen_Meson = new TH3F("h_PtYPdg_Muon_Gen_Meson", "h_PtYPdg_Muon_Gen_Meson", 300, 0.0, 30.0, 150, -4.0, -2.5, 2, 400, 600);
TH3F *h_PtYPdg_Muon_Gen_Barion = new TH3F("h_PtYPdg_Muon_Gen_Barion", "h_PtYPdg_Muon_Gen_Barion", 300, 0.0, 30.0, 150, -4.0, -2.5, 2, 4000, 6000);

//Pt-PDG for Generated Muons +
TH2F *h_PtPdg_Muon_Gen_Meson = new TH2F("h_PtPdg_Muon_Gen_Meson", "h_PtPdg_Muon_Gen_Meson", 300, 0.0, 30.0, 200, 400, 600);
TH2F *h_PtPdg_Muon_Gen_Barion = new TH2F("h_PtPdg_Muon_Gen_Barion", "h_PtPdg_Muon_Gen_Barion", 300, 0.0, 30.0, 200, 4000, 6000);

//Pt-Y for Reconstructed Muons +
TH3F *h_PtYPdg_Muon_Rec_Meson = new TH3F("h_PtYPdg_Muon_Rec_Meson", "h_PtYPdg_Muon_Rec_Meson", 300, 0.0, 30.0, 150, -4.0, -2.5, 2, 400, 600);
TH3F *h_PtYPdg_Muon_Rec_Barion = new TH3F("h_PtYPdg_Muon_Rec_Barion", "h_PtYPdg_Muon_Rec_Barion", 300, 0.0, 30.0, 150, -4.0, -2.5, 2, 4000, 6000);

//Pt-Pdg for Reconstructed Muons +
TH2F *h_PtPdg_Muon_Rec_Meson = new TH2F("h_PtPdg_Muon_Rec_Meson", "h_PtPdg_Muon_Rec_Meson", 300, 0.0, 30.0, 200, 400, 600);
TH2F *h_PtPdg_Muon_Rec_Barion = new TH2F("h_PtPdg_Muon_Rec_Barion", "h_PtPdg_Muon_Rec_Barion", 300, 0.0, 30.0, 200, 4000, 6000);


//Pt-Mass for Generated Dimuons
TH3F *h_PtMPdg_DiMu_Gen_Meson_ULS = new TH3F("h_PtMPdg_DiMu_Gen_Meson_ULS", "h_PtMPdg_DiMu_Gen_Meson_ULS", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);
TH3F *h_PtMPdg_DiMu_Gen_Barion_ULS = new TH3F("h_PtMPdg_DiMu_Gen_Barion_ULS", "h_PtMPdg_DiMu_Gen_Barion_ULS", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);
//Pt-Y for Generated Dimuons
TH3F *h_PtYPdg_DiMu_Gen_Meson_ULS = new TH3F("h_PtYPdg_DiMu_Gen_Meson_ULS", "h_PtYPdg_DiMu_Gen_Meson_ULS", 300, 0.0, 30.0, 150, -4.0, -2.5, 3, 400, 700);
TH3F *h_PtYPdg_DiMu_Gen_Barion_ULS = new TH3F("h_PtYPdg_DiMu_Gen_Barion_ULS", "h_PtYPdg_DiMu_Gen_Barion_ULS", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);

//Pt-Mass for Reconstructed Dimuons
TH3F *h_PtMPdg_DiMu_Rec_Meson_ULS = new TH3F("h_PtMPdg_DiMu_Rec_Meson_ULS", "h_PtMPdg_DiMu_Rec_Meson_ULS", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);
TH3F *h_PtMPdg_DiMu_Rec_Barion_ULS = new TH3F("h_PtMPdg_DiMu_Rec_Barion_ULS", "h_PtMPdg_DiMu_Rec_Barion_ULS", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);
//Pt-Y for Reconstructed Dimuons
TH3F *h_PtYPdg_DiMu_Rec_Meson_ULS = new TH3F("h_PtYPdg_DiMu_Rec_Meson_ULS", "h_PtYPdg_DiMu_Rec_Meson_ULS", 300, 0.0, 30.0, 150, -4.0, -2.5, 3, 400, 700);
TH3F *h_PtYPdg_DiMu_Rec_Barion_ULS = new TH3F("h_PtYPdg_DiMu_Rec_Barion_ULS", "h_PtYPdg_DiMu_Rec_Barion_ULS", 300, 0.0, 30.0, 150, -4.0, -2.5, 3, 4000, 7000);

// TH3F *h_PtMPdg_DiMu_Gen_Meson_LSsign = new TH3F("h_PtMPdg_DiMu_Gen_Meson_LSsign", "h_PtMPdg_DiMu_Gen_Meson_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);
// TH3F *h_PtMPdg_DiMu_Rec_Meson_LSsign = new TH3F("h_PtMPdg_DiMu_Rec_Meson_LSsign", "h_PtMPdg_DiMu_Rec_Meson_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);
// TH3F *h_PtMPdg_DiMu_Gen_Meson_LSsign = new TH3F("h_PtMPdg_DiMu_Gen_Meson_LSsign", "h_PtMPdg_DiMu_Gen_Meson_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);
// TH3F *h_PtMPdg_DiMu_Rec_Meson_LSsign = new TH3F("h_PtMPdg_DiMu_Rec_Meson_LSsign", "h_PtMPdg_DiMu_Rec_Meson_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 400, 700);

// TH3F *h_PtMPdg_DiMu_Gen_Barion_LSsign = new TH3F("h_PtMPdg_DiMu_Gen_Barion_LSsign", "h_PtMPdg_DiMu_Gen_Barion_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);
// TH3F *h_PtMPdg_DiMu_Rec_Barion_LSsign = new TH3F("h_PtMPdg_DiMu_Rec_Barion_LSsign", "h_PtMPdg_DiMu_Rec_Barion_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);
// TH3F *h_PtMPdg_DiMu_Gen_Barion_LSsign = new TH3F("h_PtMPdg_DiMu_Gen_Barion_LSsign", "h_PtMPdg_DiMu_Gen_Barion_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);
// TH3F *h_PtMPdg_DiMu_Rec_Barion_LSsign = new TH3F("h_PtMPdg_DiMu_Rec_Barion_LSsign", "h_PtMPdg_DiMu_Rec_Barion_LSsign", 300, 0.0, 30.0, 300, 0.0, 30.0, 3, 4000, 7000);
