#include "/home/michele_pennisi/cernbox/common_include.h"
// #include "/alidata/mpennisi/HF_dimuons/common_include.h"
Int_t N_HFquarks_gen; // gen c/cbar or b/bar HFquarks in the event

Int_t NMuons_gen;   // gen muon in the event
Int_t NHadrons_gen; // gen muon in the event
Int_t NDimu_gen;    // gen dimuons in the event
Int_t NMuons_rec;   // rec muon tracks in the event
Int_t NDimu_rec;    // rec dimuons in the event
Int_t fN_gamma;     // gen gamma* in the event

const Int_t fMuons_dim = 5000;
const Int_t fDimu_dim = 100000;

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
Int_t fFrom_Powheg_rec[fMuons_dim];     // check muon gen origin
Int_t fInitial_Parton_rec[fMuons_dim];  // check muon gen orecin
Double_t fVzmother_rec[fMuons_dim];     // check muon gen origin
Int_t fFrom_Geant_rec[fMuons_dim];      // check muon gen origin

Int_t PDG_HFquark_gen[fMuons_dim];   // single gen c/cbar PDG mum
Double_t Pt_HFquark_gen[fMuons_dim]; // single gen c/cbar or b/bbar HFquark pT
Double_t Y_HFquark_gen[fMuons_dim];  // single gen c/cbar or b/bbar HFquark y
Int_t Mother_index[fMuons_dim];      // single gen c/cbar or b/bbar HFquark y

Double_t fPt_gamma[fMuons_dim]; // Gamma star pt
Double_t fM_gamma[fMuons_dim];  // Gamma star M
Double_t fY_gamma[fMuons_dim];  // Gamma star Y

Int_t PDGmum_gen[fMuons_dim];       // single gen mu PDG mum
Double_t Pt_gen[fMuons_dim];        // single gen mu pT
Double_t E_gen[fMuons_dim];         // single gen mu E
Double_t Px_gen[fMuons_dim];        // single gen mu px
Double_t Py_gen[fMuons_dim];        // single gen mu py
Double_t Pz_gen[fMuons_dim];        // single gen mu pz
Double_t Y_gen[fMuons_dim];         // single gen mu y
Double_t Eta_gen[fMuons_dim];       // single gen mu eta
Double_t Phi_gen[fMuons_dim];       // single gen mu phi
Double_t Theta_gen[fMuons_dim];     // single gen mu theta
Int_t Charge_gen[fMuons_dim];       // single gen mu theta
Int_t fFrom_Powheg_gen[fMuons_dim]; // check muon gen origin
Int_t fInitial_Parton_gen[fMuons_dim];
Double_t fRadius_gen[fMuons_dim];
Double_t fVz_gen[fMuons_dim];
Int_t fFrom_Geant_gen[fMuons_dim];
Double_t fVzmother_gen[fMuons_dim]; // check muon gen origin

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

Int_t PDGmum_Hadron_gen[fMuons_dim];      // gen Hadron PDG mum
Int_t PDG_Hadron_gen[fMuons_dim];         // gen Hadron PDG
Double_t Pt_Hadron_gen[fMuons_dim];       // gen Hadron pT
Double_t E_Hadron_gen[fMuons_dim];        // gen Hadron E
Double_t Px_Hadron_gen[fMuons_dim];       // gen Hadron px
Double_t Py_Hadron_gen[fMuons_dim];       // gen Hadron py
Double_t Pz_Hadron_gen[fMuons_dim];       // gen Hadron pz
Double_t Y_Hadron_gen[fMuons_dim];        // gen Hadron y
Double_t Eta_Hadron_gen[fMuons_dim];      // gen Hadron eta
Int_t fHadronFrom_Powheg_gen[fMuons_dim]; // check muon gen origin
Double_t fVzHadron_gen[fMuons_dim];
Int_t fHadronFrom_Geant_gen[fMuons_dim];

TH1D *h_Nevents = new TH1D("h_Nevents", "h_Nevents", 2, 0, 2);

const Int_t n_Muon_origin = 5; // 0 for All sources, 1 for Charm, 2 for Beauty, 3 for LF, 4 for DY

TH2F *h_PtY_Charm_quark;
TH1F *h_NCharm_event;
TH1F *h_NCharm_event_fwd;
TH2F *h_PtY_Charm_quark_PowhegOnly;
TH1F *h_NCharm_event_PowhegOnly;
TH1F *h_NCharm_event_fwd_PowhegOnly;
TH2F *h_PtY_Beauty_quark;
TH1F *h_NBeauty_event;
TH1F *h_NBeauty_event_fwd;
TH2F *h_PtY_Beauty_quark_PowhegOnly;
TH1F *h_NBeauty_event_PowhegOnly;
TH1F *h_NBeauty_event_fwd_PowhegOnly;

//------Declaration Hist for Generated Muons----------------//
TH2F *h_RadiusEta_Muon_Gen_LF;
TH2F *h_VzEta_Muon_Gen_LF;

TH2F *h_VzmotherEta_Muon_Gen_PYTHIA[n_Muon_origin];
TH2F *h_VzmotherEta_Muon_Rec_PYTHIA[n_Muon_origin];

TH2F *h_VzmotherEta_Muon_Gen_Geant[n_Muon_origin];
TH2F *h_VzmotherEta_Muon_Rec_Geant[n_Muon_origin];

TH2F *h_VzHadronEta_Hadron_Gen_Charm_PYTHIA;
TH2F *h_VzHadronEta_Hadron_Gen_Beauty_PYTHIA;
TH2F *h_VzHadronEta_Hadron_Gen_LF_PYTHIA;

TH2F *h_VzHadronEta_Hadron_Gen_Charm_Geant;
TH2F *h_VzHadronEta_Hadron_Gen_Beauty_Geant;
TH2F *h_VzHadronEta_Hadron_Gen_LF_Geant;

TH2F *h_PtPdg_Muon_Gen;
TH2F *h_YPdg_Muon_Gen;
TH2F *h_EtaPdg_Muon_Gen;

TH2F *h_PtY_Muon_Gen[n_Muon_origin];
TH2F *h_PtEta_Muon_Gen[n_Muon_origin];
TH1F *h_Nperevent_Muon_Gen[n_Muon_origin];

TH2F *h_PtQ_Muon_Gen_LF_Powheg;
TH2F *h_PtQ_Muon_Gen_LF_Pythia;

TH2F *h_YQ_Muon_Gen_LF_Powheg;
TH2F *h_YQ_Muon_Gen_LF_Pythia;

TH2F *h_PtPdg_Muon_Gen_DQcut;
TH2F *h_YPdg_Muon_Gen_DQcut;
TH2F *h_EtaPdg_Muon_Gen_DQcut;

TH2F *h_PtY_Muon_Gen_DQcut[n_Muon_origin];
TH2F *h_PtEta_Muon_Gen_DQcut[n_Muon_origin];
TH1F *h_Nperevent_Muon_Gen_DQcut[n_Muon_origin];

//------ "" from PYTHIA only -----------//

TH2F *h_PtPdg_Muon_Gen_PYTHIAOnly;
TH2F *h_YPdg_Muon_Gen_PYTHIAOnly;
TH2F *h_EtaPdg_Muon_Gen_PYTHIAOnly;

TH2F *h_PtPdg_Muon_Gen_DQcut_PYTHIAOnly;
TH2F *h_YPdg_Muon_Gen_DQcut_PYTHIAOnly;
TH2F *h_EtaPdg_Muon_Gen_DQcut_PYTHIAOnly;

TH2F *h_PtY_Muon_Gen_PYTHIAOnly[n_Muon_origin];
TH2F *h_PtEta_Muon_Gen_PYTHIAOnly[n_Muon_origin];
TH1F *h_Nperevent_Muon_Gen_PYTHIAOnly[n_Muon_origin];

TH2F *h_PtY_Muon_Gen_DQcut_PYTHIAOnly[n_Muon_origin];
TH2F *h_PtEta_Muon_Gen_DQcut_PYTHIAOnly[n_Muon_origin];
TH1F *h_Nperevent_Muon_Gen_DQcut_PYTHIAOnly[n_Muon_origin];

//------ "" from GEANT only -----------//

TH2F *h_PtPdg_Muon_Gen_GeantOnly;
TH2F *h_YPdg_Muon_Gen_GeantOnly;
TH2F *h_EtaPdg_Muon_Gen_GeantOnly;

TH2F *h_PtY_Muon_Gen_GeantOnly[n_Muon_origin];
TH2F *h_PtEta_Muon_Gen_GeantOnly[n_Muon_origin];
TH1F *h_Nperevent_Muon_Gen_GeantOnly[n_Muon_origin];

TH2F *h_PtPdg_Muon_Gen_DQcut_GeantOnly;
TH2F *h_YPdg_Muon_Gen_DQcut_GeantOnly;
TH2F *h_EtaPdg_Muon_Gen_DQcut_GeantOnly;

TH2F *h_PtY_Muon_Gen_DQcut_GeantOnly[n_Muon_origin];
TH2F *h_PtEta_Muon_Gen_DQcut_GeantOnly[n_Muon_origin];
TH1F *h_Nperevent_Muon_Gen_DQcut_GeantOnly[n_Muon_origin];

//------ "" from POWHEG only -----------//

TH2F *h_PtPdg_Muon_Gen_PowhegOnly;
TH2F *h_YPdg_Muon_Gen_PowhegOnly;
TH2F *h_EtaPdg_Muon_Gen_PowhegOnly;

TH2F *h_PtY_Muon_Gen_PowhegOnly[n_Muon_origin];
TH2F *h_PtEta_Muon_Gen_PowhegOnly[n_Muon_origin];
TH1F *h_Nperevent_Muon_Gen_PowhegOnly[n_Muon_origin];

TH2F *h_PtPdg_Muon_Gen_DQcut_PowhegOnly;
TH2F *h_YPdg_Muon_Gen_DQcut_PowhegOnly;
TH2F *h_EtaPdg_Muon_Gen_DQcut_PowhegOnly;

TH2F *h_PtY_Muon_Gen_DQcut_PowhegOnly[n_Muon_origin];
TH2F *h_PtEta_Muon_Gen_DQcut_PowhegOnly[n_Muon_origin];
TH1F *h_Nperevent_Muon_Gen_DQcut_PowhegOnly[n_Muon_origin];

//------Declaration Hist for Reconstructed Muons----------------//
TH2F *h_PtPdg_Muon_Rec;
TH2F *h_YPdg_Muon_Rec;
TH2F *h_EtaPdg_Muon_Rec;

TH2F *h_PtY_Muon_Rec[n_Muon_origin];
TH2F *h_PtEta_Muon_Rec[n_Muon_origin];
TH1F *h_Nperevent_Muon_Rec[n_Muon_origin];

//------ "" from GEANT only -----------//
TH2F *h_PtPdg_Muon_Rec_GeantOnly;
TH2F *h_YPdg_Muon_Rec_GeantOnly;
TH2F *h_EtaPdg_Muon_Rec_GeantOnly;

TH2F *h_PtY_Muon_Rec_GeantOnly[n_Muon_origin];
TH2F *h_PtEta_Muon_Rec_GeantOnly[n_Muon_origin];
TH1F *h_Nperevent_Muon_Rec_GeantOnly[n_Muon_origin];

//------ "" from POWHEG only -----------//
TH2F *h_PtPdg_Muon_Rec_PowhegOnly;
TH2F *h_YPdg_Muon_Rec_PowhegOnly;
TH2F *h_EtaPdg_Muon_Rec_PowhegOnly;

TH2F *h_PtY_Muon_Rec_PowhegOnly[n_Muon_origin];
TH2F *h_PtEta_Muon_Rec_PowhegOnly[n_Muon_origin];
TH1F *h_Nperevent_Muon_Rec_PowhegOnly[n_Muon_origin];

//------ "" from PYTHIA only -----------//
TH2F *h_PtPdg_Muon_Rec_PYTHIAOnly;
TH2F *h_YPdg_Muon_Rec_PYTHIAOnly;
TH2F *h_EtaPdg_Muon_Rec_PYTHIAOnly;

TH2F *h_PtY_Muon_Rec_PYTHIAOnly[n_Muon_origin];
TH2F *h_PtEta_Muon_Rec_PYTHIAOnly[n_Muon_origin];
TH1F *h_Nperevent_Muon_Rec_PYTHIAOnly[n_Muon_origin];

TString Muon_origin[n_Muon_origin];

const Int_t n_DiMuon_origin = 6; // 0 for Charm, 1 for Beauty, 2 for HF Mixed, 3 for LF, 4 for LF-HF Mixed, 5 for DY

//------Declaration Hist for Generated Dimuons----------------//

TH3F *h_Pdg1Pdg2Pt_DiMuon_Gen;
TH3F *h_Pdg1Pdg2Y_DiMuon_Gen;
TH3F *h_Pdg1Pdg2M_DiMuon_Gen;

TH2F *h_PtM_DiMuon_Gen[n_DiMuon_origin];
TH2F *h_PtY_DiMuon_Gen[n_DiMuon_origin];
TH1F *h_Nperevent_DiMuon_Gen[n_DiMuon_origin];

TH3F *h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut;
TH3F *h_Pdg1Pdg2Y_DiMuon_Gen_DQcut;
TH3F *h_Pdg1Pdg2M_DiMuon_Gen_DQcut;

TH2F *h_PtM_DiMuon_Gen_DQcut[n_DiMuon_origin];
TH2F *h_PtY_DiMuon_Gen_DQcut[n_DiMuon_origin];
TH1F *h_Nperevent_DiMuon_Gen_DQcut[n_DiMuon_origin];

TH2F *h_PtM_DiMuon_Gen_DQcut_Charm_corrected = new TH2F("h_PtM_DiMuon_Gen_DQcut_Charm_corrected", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 300, 0, 30.0, 300, 0, 30);
TH2F *h_PtY_DiMuon_Gen_DQcut_Charm_corrected = new TH2F("h_PtY_DiMuon_Gen_DQcut_Charm_corrected", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 300, 0, 30.0, 150, -4.0, -2.5);
TH3F *h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_Charm_corrected;
TH3F *h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_Charm_corrected;
TH3F *h_Pdg1Pdg2M_DiMuon_Gen_DQcut_Charm_corrected;

//------ "" from powheg only -----------//

TH3F *h_Pdg1Pdg2Pt_DiMuon_Gen_PowhegOnly;
TH3F *h_Pdg1Pdg2Y_DiMuon_Gen_PowhegOnly;
TH3F *h_Pdg1Pdg2M_DiMuon_Gen_PowhegOnly;

TH2F *h_PtM_DiMuon_Gen_PowhegOnly[n_DiMuon_origin];
TH2F *h_PtY_DiMuon_Gen_PowhegOnly[n_DiMuon_origin];
TH1F *h_Nperevent_DiMuon_Gen_PowhegOnly[n_DiMuon_origin];

TH3F *h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_PowhegOnly;
TH3F *h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_PowhegOnly;
TH3F *h_Pdg1Pdg2M_DiMuon_Gen_DQcut_PowhegOnly;

TH2F *h_PtM_DiMuon_Gen_DQcut_PowhegOnly[n_DiMuon_origin];
TH2F *h_PtY_DiMuon_Gen_DQcut_PowhegOnly[n_DiMuon_origin];
TH1F *h_Nperevent_DiMuon_Gen_DQcut_PowhegOnly[n_DiMuon_origin];

//------Declaration Hist for Reconstructed Dimuons----------------//

TH3F *h_Pdg1Pdg2Pt_DiMuon_Rec;
TH3F *h_Pdg1Pdg2Y_DiMuon_Rec;
TH3F *h_Pdg1Pdg2M_DiMuon_Rec;

TH2F *h_PtM_DiMuon_Rec[n_DiMuon_origin];
TH2F *h_PtY_DiMuon_Rec[n_DiMuon_origin];
TH1F *h_Nperevent_DiMuon_Rec[n_DiMuon_origin];

TH2F **h_PtM_DiMuon_Rec_fromLF;
TH2F **h_PtY_DiMuon_Rec_fromLF;
TH1F **h_Nperevent_DiMuon_Rec_fromLF;

TH2F **h_PtM_DiMuon_Rec_fromLF_HF_Mixed;
TH2F **h_PtY_DiMuon_Rec_fromLF_HF_Mixed;
TH1F **h_Nperevent_DiMuon_Rec_fromLF_HF_Mixed;

TH3F *h_Pdg1Pdg2Pt_DiMuon_Rec_PowhegOnly;
TH3F *h_Pdg1Pdg2Y_DiMuon_Rec_PowhegOnly;
TH3F *h_Pdg1Pdg2M_DiMuon_Rec_PowhegOnly;

TH2F *h_PtM_DiMuon_Rec_PowhegOnly[n_DiMuon_origin];
TH2F *h_PtY_DiMuon_Rec_PowhegOnly[n_DiMuon_origin];
TH1F *h_Nperevent_DiMuon_Rec_PowhegOnly[n_DiMuon_origin];

//------ old staff for Z analysis on the fly -----------//
TH2F *h_PtM_DiMuon_Gen_Z_ptmucut09;
TH2F *h_PtY_DiMuon_Gen_Z_ptmucut09;
TH1F *h_Nperevent_DiMuon_Gen_Z_ptmucut09;

TH2F *h_PtM_DiMuon_Gen_Z_ptmucut10;
TH2F *h_PtY_DiMuon_Gen_Z_ptmucut10;
TH1F *h_Nperevent_DiMuon_Gen_Z_ptmucut10;

TH2F *h_PtM_DiMuon_Gen_Z_ptmucut20;
TH2F *h_PtY_DiMuon_Gen_Z_ptmucut20;
TH1F *h_Nperevent_DiMuon_Gen_Z_ptmucut20;

TH2F *h_PtM_DiMuon_Gen_Z_DQcut_ptmucut09;
TH2F *h_PtY_DiMuon_Gen_Z_DQcut_ptmucut09;
TH1F *h_Nperevent_DiMuon_Gen_Z_DQcut_ptmucut09;

TH2F *h_PtM_DiMuon_Gen_Z_DQcut_ptmucut10;
TH2F *h_PtY_DiMuon_Gen_Z_DQcut_ptmucut10;
TH1F *h_Nperevent_DiMuon_Gen_Z_DQcut_ptmucut10;

TH2F *h_PtM_DiMuon_Gen_Z_DQcut_ptmucut20;
TH2F *h_PtY_DiMuon_Gen_Z_DQcut_ptmucut20;
TH1F *h_Nperevent_DiMuon_Gen_Z_DQcut_ptmucut20;

TH2F *h_PtM_DiMuon_Rec_Z_ptmucut09;
TH2F *h_PtY_DiMuon_Rec_Z_ptmucut09;
TH1F *h_Nperevent_DiMuon_Rec_Z_ptmucut09;

TH2F *h_PtM_DiMuon_Rec_Z_ptmucut10;
TH2F *h_PtY_DiMuon_Rec_Z_ptmucut10;
TH1F *h_Nperevent_DiMuon_Rec_Z_ptmucut10;

TH2F *h_PtM_DiMuon_Rec_Z_ptmucut20;
TH2F *h_PtY_DiMuon_Rec_Z_ptmucut20;
TH1F *h_Nperevent_DiMuon_Rec_Z_ptmucut20;

//---------------------------------//

TString DiMuon_origin[n_DiMuon_origin];
TString *DiMuon_fromLF_Generator;

const Int_t n_PDG_selection = 46;
Double_t PDG_Selection[n_PDG_selection] = {0, 23, 24, 125, 135, 200, 210, 220, 230, 240, 250, 300, 310, 320, 330, 340, 400, 410, 420, 430, 440, 450, 500, 510, 520, 530, 540, 550, 600, 2000, 3000, 4000, 4000, 4100, 4120, 4130, 4140, 4200, 4230, 4240, 4300, 5000, 5100, 5200, 5300, 6000};
const Int_t n_charm_hadrons = 7;
Int_t PDG_charm_hadrons[n_charm_hadrons] = {411, 421, 431, 443, 4122, 4132, 4232};

Double_t BR_charm_hadrons2mu_PYTHIA8_Monash[n_charm_hadrons] = {16.5, 6.45, 7.5, 5.9, 4.5, 2.5, 3.5};
Double_t BR_charm_hadrons2mu_MEAS[n_charm_hadrons] = {17.6, 6.80, 6.33, 5.96, 3.95, 2.5, 3.5};

Double_t Frag_charm_hadrons_PYTHIA8_Monash[n_charm_hadrons] = {29.3, 56.1, 9.59, 0.4, 3.8, 0.49, 0.49};
Double_t Frag_charm_hadrons_MEAS[n_charm_hadrons] = {19.1, 38.2, 6.1, 0.37, 16.8, 9.9, 9.6};

//---- Declarion hist for Hadrons --------//

TH2F *h_PdgPt_HFHadron_prompt;
TH2F *h_PdgPt_HFHadron_notprompt;
TH2F *h_PdgY_HFHadron_prompt;
TH2F *h_PdgY_HFHadron_notprompt;
TH2F *h_PdgEta_HFHadron_prompt;
TH2F *h_PdgEta_HFHadron_notprompt;
TH3F *h_PdgPtY_HFHadron_prompt;
TH3F *h_PdgPtY_HFHadron_notprompt;

//----- "" for pythia only ------------//

TH2F *h_PdgPt_HFHadron_prompt_PYTHIAOnly;
TH2F *h_PdgPt_HFHadron_notprompt_PYTHIAOnly;
TH2F *h_PdgY_HFHadron_prompt_PYTHIAOnly;
TH2F *h_PdgY_HFHadron_notprompt_PYTHIAOnly;
TH2F *h_PdgEta_HFHadron_prompt_PYTHIAOnly;
TH2F *h_PdgEta_HFHadron_notprompt_PYTHIAOnly;
TH3F *h_PdgPtY_HFHadron_prompt_PYTHIAOnly;
TH3F *h_PdgPtY_HFHadron_notprompt_PYTHIAOnly;

//----- "" for Geant only ------------//

TH2F *h_PdgPt_HFHadron_prompt_GeantOnly;
TH2F *h_PdgPt_HFHadron_notprompt_GeantOnly;
TH2F *h_PdgY_HFHadron_prompt_GeantOnly;
TH2F *h_PdgY_HFHadron_notprompt_GeantOnly;
TH2F *h_PdgEta_HFHadron_prompt_GeantOnly;
TH2F *h_PdgEta_HFHadron_notprompt_GeantOnly;

TH3F *h_PdgPtY_HFHadron_prompt_GeantOnly;
TH3F *h_PdgPtY_HFHadron_notprompt_GeantOnly;

//----- "" for powheg only ------------//

TH2F *h_PdgPt_HFHadron_prompt_PowhegOnly;
TH2F *h_PdgPt_HFHadron_notprompt_PowhegOnly;
TH2F *h_PdgY_HFHadron_prompt_PowhegOnly;
TH2F *h_PdgY_HFHadron_notprompt_PowhegOnly;

TH2F *h_PdgEta_HFHadron_prompt_PowhegOnly;
TH2F *h_PdgEta_HFHadron_notprompt_PowhegOnly;

TH3F *h_PdgPtY_HFHadron_prompt_PowhegOnly;
TH3F *h_PdgPtY_HFHadron_notprompt_PowhegOnly;

TH2F *h_PtM_Gamma;
TH2F *h_PtY_Gamma;

TH2F *h_YGamma_YDimuon;
void Set_Histograms(TString Generator)
{
    h_VzEta_Muon_Gen_LF = new TH2F("h_VzEta_Muon_Gen_LF", ";Vz (cm?);#eta", 6000, -3000, 3000, 200, -10, 10);
    h_RadiusEta_Muon_Gen_LF = new TH2F("h_RadiusEta_Muon_Gen_LF", ";radius (cm?); #eta", 5000, 0, 5000, 200, -10, 10);

    h_VzHadronEta_Hadron_Gen_Charm_PYTHIA = new TH2F("h_VzHadronEta_Hadron_Gen_Charm_PYTHIA", ";Vz(cm?);#eta", 6000, -3000, 3000, 200, -10, 10);
    h_VzHadronEta_Hadron_Gen_Beauty_PYTHIA = new TH2F("h_VzHadroEtaY_Hadron_Gen_Beauty_PYTHIA", ";Vz(cm?);#eta", 6000, -3000, 3000, 200, -10, 10);
    h_VzHadronEta_Hadron_Gen_LF_PYTHIA = new TH2F("h_VzHadronEta_Hadron_Gen_LF_PYTHIA", ";Vz(cm?);#eta", 6000, -3000, 3000, 200, -10, 10);

    h_VzHadronEta_Hadron_Gen_Charm_Geant = new TH2F("h_VzHadronEta_Hadron_Gen_Charm_Geant", ";Vz(cm?);#eta", 6000, -3000, 3000, 200, -10, 10);
    h_VzHadronEta_Hadron_Gen_Beauty_Geant = new TH2F("h_VzHadroEtaY_Hadron_Gen_Beauty_Geant", ";Vz(cm?);#eta", 6000, -3000, 3000, 200, -10, 10);
    h_VzHadronEta_Hadron_Gen_LF_Geant = new TH2F("h_VzHadronEta_Hadron_Gen_LF_Geant", ";Vz(cm?);#eta", 6000, -3000, 3000, 200, -10, 10);

    h_PtY_Charm_quark = new TH2F("h_PtY_Charm_quark", "#it{p}_{T} (GeV/#it{c}); #it{y}", 400, 0.0, 40.0, 200, -10.0, 10.0);
    h_NCharm_event = new TH1F("h_NCharm_event", "c quark x ev", 20, -0.5, 19.5);
    h_NCharm_event_fwd = new TH1F("h_NCharm_event_fwd", "c quark x ev fwd", 20, -0.5, 19.5);

    h_PtY_Charm_quark_PowhegOnly = new TH2F("h_PtY_Charm_quark_PowhegOnly", "#it{p}_{T} (GeV/#it{c}); #it{y}", 400, 0.0, 40.0, 200, -10.0, 10.0);
    h_NCharm_event_PowhegOnly = new TH1F("h_NCharm_event_PowhegOnly", "c quark x ev", 20, -0.5, 19.5);
    h_NCharm_event_fwd_PowhegOnly = new TH1F("h_NCharm_event_fwd_PowhegOnly", "c quark x ev fwd", 20, -0.5, 19.5);

    h_PtY_Beauty_quark = new TH2F("h_PtY_Beauty_quark", "#it{p}_{T} (GeV/#it{c}); #it{y}", 400, 0.0, 40.0, 200, -10.0, 10.0);
    h_NBeauty_event = new TH1F("h_NBeauty_event", "b quark x ev", 20, -0.5, 19.5);
    h_NBeauty_event_fwd = new TH1F("h_NBeauty_event_fwd", "b quark x ev fwd", 20, -0.5, 19.5);

    h_PtY_Beauty_quark_PowhegOnly = new TH2F("h_PtY_Beauty_quark_PowhegOnly", "#it{p}_{T} (GeV/#it{c}); #it{y}", 400, 0.0, 40.0, 200, -10.0, 10.0);
    h_NBeauty_event_PowhegOnly = new TH1F("h_NBeauty_event_PowhegOnly", "b quark x ev", 20, -0.5, 19.5);
    h_NBeauty_event_fwd_PowhegOnly = new TH1F("h_NBeauty_event_fwd_PowhegOnly", "b quark x ev fwd", 20, -0.5, 19.5);

    h_PtQ_Muon_Gen_LF_Powheg = new TH2F("h_PtQ_Muon_Gen_LF_Powheg", "; #it{p}_{T} (GeV/#it{c}) ; initial parton", 400, 0.0, 40.0, 4000, 0, 4000);
    h_PtQ_Muon_Gen_LF_Pythia = new TH2F("h_PtQ_Muon_Gen_LF_Pythia", "; #it{p}_{T} (GeV/#it{c}) ; initial parton", 400, 0.0, 40.0, 4000, 0, 4000);

    h_YQ_Muon_Gen_LF_Powheg = new TH2F("h_YQ_Muon_Gen_LF_Powheg", "; #it{y}; initial parton", 160, -8, 8, 4000, 0, 4000);
    h_YQ_Muon_Gen_LF_Pythia = new TH2F("h_YQ_Muon_Gen_LF_Pythia", "; #it{y}; initial parton", 160, -8, 8, 4000, 0, 4000);
    // Histograms for saving Gamma star kinematics
    h_PtM_Gamma = new TH2F("h_PtM_Gamma", "; #it{p}_{T,#gamma^{*}} (GeV/#it{c}); #it{m}_{#gamma^{*}} (GeV/#it{c}^{2})", 400, 0, 40.0, 400, 0, 40.0);
    h_PtY_Gamma = new TH2F("h_PtY_Gamma", "; #it{p}_{T,#gamma^{*}} (GeV/#it{c}); #it{y}_{#gamma^{*}}", 400, 0, 40.0, 160, -8, 8);
    h_YGamma_YDimuon = new TH2F("h_YGamma_YDimuon", "; #it{y}_{#gamma^{*}}; #it{y}_{#mu#mu}", 160, -8, 8, 160, -8, 8);

    // Saving Generated Muons kinematic/PDG ancestor correlation

    h_PtPdg_Muon_Gen = new TH2F("h_PtPdg_Muon_Gen", "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 400, 0, 40.0, 6001, 0, 6000);
    h_YPdg_Muon_Gen = new TH2F("h_YPdg_Muon_Gen", "; #it{y} ; PDG code", 160, -8.0, 8.0, 6001, 0, 6000);
    h_EtaPdg_Muon_Gen = new TH2F("h_EtaPdg_Muon_Gen", "; #eta ; PDG code", 160, -8.0, 8.0, 6001, 0, 6000);

    if (Generator.Contains("Geant"))
    {

        h_PtPdg_Muon_Gen_PYTHIAOnly = new TH2F("h_PtPdg_Muon_Gen_PYTHIAOnly", "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 400, 0, 40.0, 6001, 0, 6000);
        h_YPdg_Muon_Gen_PYTHIAOnly = new TH2F("h_YPdg_Muon_Gen_PYTHIAOnly", "; #it{y} ; PDG code", 160, -8.0, 8.0, 6001, 0, 6000);
        h_EtaPdg_Muon_Gen_PYTHIAOnly = new TH2F("h_EtaPdg_Muon_Gen_PYTHIAOnly", "; #eta ; PDG code", 160, -8.0, 8.0, 6001, 0, 6000);

        h_PtPdg_Muon_Gen_GeantOnly = new TH2F("h_PtPdg_Muon_Gen_GeantOnly", "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 400, 0, 40.0, 6001, 0, 6000);
        h_YPdg_Muon_Gen_GeantOnly = new TH2F("h_YPdg_Muon_Gen_GeantOnly", "; #it{y} ; PDG code", 160, -8.0, 8.0, 6001, 0, 6000);
        h_EtaPdg_Muon_Gen_GeantOnly = new TH2F("h_EtaPdg_Muon_Gen_GeantOnly", "; #eta ; PDG code", 160, -8.0, 8.0, 6001, 0, 6000);

        if (Generator.Contains("Powheg"))
        {
            h_PtPdg_Muon_Gen_PowhegOnly = new TH2F("h_PtPdg_Muon_Gen_PowhegOnly", "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 400, 0, 40.0, 6001, 0, 6000);
            h_YPdg_Muon_Gen_PowhegOnly = new TH2F("h_YPdg_Muon_Gen_PowhegOnly", "; #it{y} ; PDG code", 160, -8.0, 8.0, 6001, 0, 6000);
            h_EtaPdg_Muon_Gen_PowhegOnly = new TH2F("h_EtaPdg_Muon_Gen_PowhegOnly", "; #eta ; PDG code", 160, -8.0, 8.0, 6001, 0, 6000);
        }
    }

    h_PtPdg_Muon_Gen_DQcut = new TH2F("h_PtPdg_Muon_Gen_DQcut", "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 400, 0, 40.0, 6001, 0, 6000);
    h_YPdg_Muon_Gen_DQcut = new TH2F("h_YPdg_Muon_Gen_DQcut", "; #it{y} ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);
    h_EtaPdg_Muon_Gen_DQcut = new TH2F("h_EtaPdg_Muon_Gen_DQcut", "; #eta ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);

    if (Generator.Contains("Geant"))
    {
        h_PtPdg_Muon_Gen_DQcut_PYTHIAOnly = new TH2F("h_PtPdg_Muon_Gen_DQcut_PYTHIAOnly", "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 400, 0, 40.0, 6001, 0, 6000);
        h_YPdg_Muon_Gen_DQcut_PYTHIAOnly = new TH2F("h_YPdg_Muon_Gen_DQcut_PYTHIAOnly", "; #it{y} ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);
        h_EtaPdg_Muon_Gen_DQcut_PYTHIAOnly = new TH2F("h_EtaPdg_Muon_Gen_DQcut_PYTHIAOnly", "; #eta ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);

        h_PtPdg_Muon_Gen_DQcut_GeantOnly = new TH2F("h_PtPdg_Muon_Gen_DQcut_GeantOnly", "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 400, 0, 40.0, 6001, 0, 6000);
        h_YPdg_Muon_Gen_DQcut_GeantOnly = new TH2F("h_YPdg_Muon_Gen_DQcut_GeantOnly", "; #it{y} ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);
        h_EtaPdg_Muon_Gen_DQcut_GeantOnly = new TH2F("h_EtaPdg_Muon_Gen_DQcut_GeantOnly", "; #eta ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);

        if (Generator.Contains("Powheg"))
        {
            h_PtPdg_Muon_Gen_DQcut_PowhegOnly = new TH2F("h_PtPdg_Muon_Gen_DQcut_PowhegOnly", "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 400, 0, 40.0, 6001, 0, 6000);
            h_YPdg_Muon_Gen_DQcut_PowhegOnly = new TH2F("h_YPdg_Muon_Gen_DQcut_PowhegOnly", "; #it{y} ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);
            h_EtaPdg_Muon_Gen_DQcut_PowhegOnly = new TH2F("h_EtaPdg_Muon_Gen_DQcut_PowhegOnly", "; #eta ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);
        }
    }

    h_PtPdg_Muon_Rec = new TH2F("h_PtPdg_Muon_Rec", "; #it{p}_{T} (GeV/#it{c}) ; PDG code", 400, 0, 40.0, 6001, 0, 6000);
    h_YPdg_Muon_Rec = new TH2F("h_YPdg_Muon_Rec", "; #it{y} ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);
    h_EtaPdg_Muon_Rec = new TH2F("h_EtaPdg_Muon_Rec", "; #eta ; PDG code", 150, -4.0, -2.5, 6001, 0, 6000);

    if (Generator.Contains("Geant"))
    {

        h_PtPdg_Muon_Rec_PYTHIAOnly = (TH2F *)h_PtPdg_Muon_Rec->Clone("h_PtPdg_Muon_Rec_PYTHIAOnly");
        h_YPdg_Muon_Rec_PYTHIAOnly = (TH2F *)h_YPdg_Muon_Rec->Clone("h_YPdg_Muon_Rec_PYTHIAOnly");
        h_EtaPdg_Muon_Rec_PYTHIAOnly = (TH2F *)h_YPdg_Muon_Rec->Clone("h_EtaPdg_Muon_Rec_PYTHIAOnly");

        h_PtPdg_Muon_Rec_GeantOnly = (TH2F *)h_PtPdg_Muon_Rec->Clone("h_PtPdg_Muon_Rec_GeantOnly");
        h_YPdg_Muon_Rec_GeantOnly = (TH2F *)h_YPdg_Muon_Rec->Clone("h_YPdg_Muon_Rec_GeantOnly");
        h_EtaPdg_Muon_Rec_GeantOnly = (TH2F *)h_YPdg_Muon_Rec->Clone("h_EtaPdg_Muon_Rec_GeantOnly");

        if (Generator.Contains("Powheg"))
        {
            h_PtPdg_Muon_Rec_PowhegOnly = (TH2F *)h_PtPdg_Muon_Rec->Clone("h_PtPdg_Muon_Rec_PowhegOnly");
            h_YPdg_Muon_Rec_PowhegOnly = (TH2F *)h_YPdg_Muon_Rec->Clone("h_YPdg_Muon_Rec_PowhegOnly");
            h_EtaPdg_Muon_Rec_PowhegOnly = (TH2F *)h_YPdg_Muon_Rec->Clone("h_EtaPdg_Muon_Rec_PowhegOnly");
        }
    }

    // 0 for All sources, 1 for Charm, 2 for Beauty, 3 for LF, 4 for DY

    Muon_origin[0].Form("All");
    Muon_origin[1].Form("Charm");
    Muon_origin[2].Form("Beauty");
    Muon_origin[3].Form("LF");
    Muon_origin[4].Form("DY");

    for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
    {
        // Gen Muons
        h_PtY_Muon_Gen[i_Muon_origin] = new TH2F(Form("h_PtY_Muon_Gen_%s", Muon_origin[i_Muon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ;#it{y}", 500, 0, 100, 900, -8.0, 8.0);
        h_PtEta_Muon_Gen[i_Muon_origin] = new TH2F(Form("h_PtEta_Muon_Gen_%s", Muon_origin[i_Muon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ;#eta", 500, 0, 100, 900, -8.0, 8.0);
        h_Nperevent_Muon_Gen[i_Muon_origin] = new TH1F(Form("h_Nperevent_Muon_Gen_%s", Muon_origin[i_Muon_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);
        // Gen Muons DQ Cut
        h_PtY_Muon_Gen_DQcut[i_Muon_origin] = new TH2F(Form("h_PtY_Muon_Gen_DQcut_%s", Muon_origin[i_Muon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ;#it{y}", 500, 0, 100, 150, -4.0, -2.5);
        h_PtEta_Muon_Gen_DQcut[i_Muon_origin] = new TH2F(Form("h_PtEta_Muon_Gen_DQcut_%s", Muon_origin[i_Muon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ;#it{y}", 500, 0, 100, 150, -4.0, -2.5);
        h_Nperevent_Muon_Gen_DQcut[i_Muon_origin] = new TH1F(Form("h_Nperevent_Muon_Gen_DQcut_%s", Muon_origin[i_Muon_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);
        // Rec Muons
        h_PtY_Muon_Rec[i_Muon_origin] = new TH2F(Form("h_PtY_Muon_Rec_%s", Muon_origin[i_Muon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ;#it{y}", 500, 0, 100, 150, -4.0, -2.5);
        h_PtEta_Muon_Rec[i_Muon_origin] = new TH2F(Form("h_PtEta_Muon_Rec%s", Muon_origin[i_Muon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ;#eta", 500, 0, 100, 150, -4.0, -2.5);
        h_Nperevent_Muon_Rec[i_Muon_origin] = new TH1F(Form("h_Nperevent_Muon_Rec_%s", Muon_origin[i_Muon_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);

        if (Generator.Contains("Geant"))
        {
            // Gen Muons PYTHIA Only
            h_VzmotherEta_Muon_Gen_PYTHIA[i_Muon_origin] = new TH2F(Form("h_VzmotherEta_Muon_Gen_PYTHIA_%s", Muon_origin[i_Muon_origin].Data()), ";v_{z} mother (cm?);#eta", 60000, -3000, 3000, 200, -10, 10);
            h_PtY_Muon_Gen_PYTHIAOnly[i_Muon_origin] = (TH2F *)h_PtY_Muon_Gen[i_Muon_origin]->Clone(Form("%s_PYTHIAOnly", h_PtY_Muon_Gen[i_Muon_origin]->GetName()));
            h_PtEta_Muon_Gen_PYTHIAOnly[i_Muon_origin] = (TH2F *)h_PtEta_Muon_Gen[i_Muon_origin]->Clone(Form("%s_PYTHIAOnly", h_PtEta_Muon_Gen[i_Muon_origin]->GetName()));
            h_Nperevent_Muon_Gen_PYTHIAOnly[i_Muon_origin] = (TH1F *)h_Nperevent_Muon_Gen[i_Muon_origin]->Clone(Form("%s_PYTHIAOnly", h_Nperevent_Muon_Gen[i_Muon_origin]->GetName()));
            // Gen Muons DQ Cut PYTHIA Only
            h_PtY_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_origin] = (TH2F *)h_PtY_Muon_Gen_DQcut[i_Muon_origin]->Clone(Form("%s_PYTHIAOnly", h_PtY_Muon_Gen_DQcut[i_Muon_origin]->GetName()));
            h_PtEta_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_origin] = (TH2F *)h_PtEta_Muon_Gen_DQcut[i_Muon_origin]->Clone(Form("%s_PYTHIAOnly", h_PtEta_Muon_Gen_DQcut[i_Muon_origin]->GetName()));
            h_Nperevent_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_origin] = (TH1F *)h_Nperevent_Muon_Gen_DQcut[i_Muon_origin]->Clone(Form("%s_PYTHIAOnly", h_Nperevent_Muon_Gen_DQcut[i_Muon_origin]->GetName()));
            // Rec Muons PYTHIA Only
            h_VzmotherEta_Muon_Rec_PYTHIA[i_Muon_origin] = new TH2F(Form("h_VzmotherEta_Muon_Rec_PYTHIA_%s", Muon_origin[i_Muon_origin].Data()), ";v_{z} mother (cm?);#eta", 6000, -3000, 3000, 150, -4.0, -2.5);
            h_PtY_Muon_Rec_PYTHIAOnly[i_Muon_origin] = (TH2F *)h_PtY_Muon_Rec[i_Muon_origin]->Clone(Form("%s_PYTHIAOnly", h_PtY_Muon_Rec[i_Muon_origin]->GetName()));
            h_PtEta_Muon_Rec_PYTHIAOnly[i_Muon_origin] = (TH2F *)h_PtEta_Muon_Rec[i_Muon_origin]->Clone(Form("%s_PYTHIAOnly", h_PtEta_Muon_Rec[i_Muon_origin]->GetName()));
            h_Nperevent_Muon_Rec_PYTHIAOnly[i_Muon_origin] = (TH1F *)h_Nperevent_Muon_Rec[i_Muon_origin]->Clone(Form("%s_PYTHIAOnly", h_Nperevent_Muon_Rec[i_Muon_origin]->GetName()));

            // Gen Muons Geant Only
            h_VzmotherEta_Muon_Gen_Geant[i_Muon_origin] = new TH2F(Form("h_VzmotherEta_Muon_Gen_Geant_%s", Muon_origin[i_Muon_origin].Data()), ";v_{z} mother (cm?);#eta", 60000, -3000, 3000, 200, -10, 10);
            h_PtY_Muon_Gen_GeantOnly[i_Muon_origin] = (TH2F *)h_PtY_Muon_Gen[i_Muon_origin]->Clone(Form("%s_GeantOnly", h_PtY_Muon_Gen[i_Muon_origin]->GetName()));
            h_PtEta_Muon_Gen_GeantOnly[i_Muon_origin] = (TH2F *)h_PtEta_Muon_Gen[i_Muon_origin]->Clone(Form("%s_GeantOnly", h_PtEta_Muon_Gen[i_Muon_origin]->GetName()));
            h_Nperevent_Muon_Gen_GeantOnly[i_Muon_origin] = (TH1F *)h_Nperevent_Muon_Gen[i_Muon_origin]->Clone(Form("%s_GeantOnly", h_Nperevent_Muon_Gen[i_Muon_origin]->GetName()));
            // Gen Muons DQ Cut Geant Only
            h_PtY_Muon_Gen_DQcut_GeantOnly[i_Muon_origin] = (TH2F *)h_PtY_Muon_Gen_DQcut[i_Muon_origin]->Clone(Form("%s_GeantOnly", h_PtY_Muon_Gen_DQcut[i_Muon_origin]->GetName()));
            h_PtEta_Muon_Gen_DQcut_GeantOnly[i_Muon_origin] = (TH2F *)h_PtEta_Muon_Gen_DQcut[i_Muon_origin]->Clone(Form("%s_GeantOnly", h_PtEta_Muon_Gen_DQcut[i_Muon_origin]->GetName()));
            h_Nperevent_Muon_Gen_DQcut_GeantOnly[i_Muon_origin] = (TH1F *)h_Nperevent_Muon_Gen_DQcut[i_Muon_origin]->Clone(Form("%s_GeantOnly", h_Nperevent_Muon_Gen_DQcut[i_Muon_origin]->GetName()));
            // Rec Muons Geant Only
            h_VzmotherEta_Muon_Rec_Geant[i_Muon_origin] = new TH2F(Form("h_VzmotherEta_Muon_Rec_Geant_%s", Muon_origin[i_Muon_origin].Data()), ";v_{z} mother (cm?);#eta", 6000, -3000, 3000, 150, -4.0, -2.5);
            h_PtY_Muon_Rec_GeantOnly[i_Muon_origin] = (TH2F *)h_PtY_Muon_Rec[i_Muon_origin]->Clone(Form("%s_GeantOnly", h_PtY_Muon_Rec[i_Muon_origin]->GetName()));
            h_PtEta_Muon_Rec_GeantOnly[i_Muon_origin] = (TH2F *)h_PtEta_Muon_Rec[i_Muon_origin]->Clone(Form("%s_GeantOnly", h_PtEta_Muon_Rec[i_Muon_origin]->GetName()));
            h_Nperevent_Muon_Rec_GeantOnly[i_Muon_origin] = (TH1F *)h_Nperevent_Muon_Rec[i_Muon_origin]->Clone(Form("%s_GeantOnly", h_Nperevent_Muon_Rec[i_Muon_origin]->GetName()));

            if (Generator.Contains("Powheg"))
            {
                // Gen Muons POWHEG Only
                h_PtY_Muon_Gen_PowhegOnly[i_Muon_origin] = (TH2F *)h_PtY_Muon_Gen[i_Muon_origin]->Clone(Form("%s_PowhegOnly", h_PtY_Muon_Gen[i_Muon_origin]->GetName()));
                h_PtEta_Muon_Gen_PowhegOnly[i_Muon_origin] = (TH2F *)h_PtEta_Muon_Gen[i_Muon_origin]->Clone(Form("%s_PowhegOnly", h_PtEta_Muon_Gen[i_Muon_origin]->GetName()));
                h_Nperevent_Muon_Gen_PowhegOnly[i_Muon_origin] = (TH1F *)h_Nperevent_Muon_Gen[i_Muon_origin]->Clone(Form("%s_PowhegOnly", h_Nperevent_Muon_Gen[i_Muon_origin]->GetName()));
                // Gen Muons DQ Cut POWHEG Only
                h_PtY_Muon_Gen_DQcut_PowhegOnly[i_Muon_origin] = (TH2F *)h_PtY_Muon_Gen_DQcut[i_Muon_origin]->Clone(Form("%s_PowhegOnly", h_PtY_Muon_Gen_DQcut[i_Muon_origin]->GetName()));
                h_PtEta_Muon_Gen_DQcut_PowhegOnly[i_Muon_origin] = (TH2F *)h_PtEta_Muon_Gen_DQcut[i_Muon_origin]->Clone(Form("%s_PowhegOnly", h_PtEta_Muon_Gen_DQcut[i_Muon_origin]->GetName()));
                h_Nperevent_Muon_Gen_DQcut_PowhegOnly[i_Muon_origin] = (TH1F *)h_Nperevent_Muon_Gen_DQcut[i_Muon_origin]->Clone(Form("%s_PowhegOnly", h_Nperevent_Muon_Gen_DQcut[i_Muon_origin]->GetName()));
                // Rec Muons POWHEG Only
                h_PtY_Muon_Rec_PowhegOnly[i_Muon_origin] = (TH2F *)h_PtY_Muon_Rec[i_Muon_origin]->Clone(Form("%s_PowhegOnly", h_PtY_Muon_Rec[i_Muon_origin]->GetName()));
                h_PtEta_Muon_Rec_PowhegOnly[i_Muon_origin] = (TH2F *)h_PtEta_Muon_Rec[i_Muon_origin]->Clone(Form("%s_PowhegOnly", h_PtEta_Muon_Rec[i_Muon_origin]->GetName()));
                h_Nperevent_Muon_Rec_PowhegOnly[i_Muon_origin] = (TH1F *)h_Nperevent_Muon_Rec[i_Muon_origin]->Clone(Form("%s_PowhegOnly", h_Nperevent_Muon_Rec[i_Muon_origin]->GetName()));
            }
        }
    }

    Int_t n_bin_pt = h_PtPdg_Muon_Gen_DQcut->GetXaxis()->GetNbins();
    Double_t low_bin_pt[n_bin_pt + 1];
    for (Int_t i_lowpt = 0; i_lowpt <= (Int_t)n_bin_pt; i_lowpt++)
        low_bin_pt[i_lowpt] = h_PtPdg_Muon_Gen_DQcut->GetXaxis()->GetBinLowEdge(i_lowpt);

    Int_t n_bin_Y_Gen = h_YPdg_Muon_Gen->GetXaxis()->GetNbins();
    Double_t low_bin_Y_Gen[n_bin_Y_Gen + 1];
    for (Int_t i_lowY_Gen = 0; i_lowY_Gen <= (Int_t)n_bin_Y_Gen; i_lowY_Gen++)
        low_bin_Y_Gen[i_lowY_Gen] = h_YPdg_Muon_Gen->GetXaxis()->GetBinLowEdge(i_lowY_Gen);

    Int_t n_bin_Y_Rec = h_YPdg_Muon_Rec->GetXaxis()->GetNbins();
    Double_t low_bin_Y_Rec[n_bin_Y_Rec + 1];
    for (Int_t i_lowY_Rec = 0; i_lowY_Rec <= (Int_t)n_bin_Y_Rec; i_lowY_Rec++)
        low_bin_Y_Rec[i_lowY_Rec] = h_YPdg_Muon_Rec->GetXaxis()->GetBinLowEdge(i_lowY_Rec);

    // 0 for Charm, 1 for Beauty, 2 for HF Mixed, 3 for LF, 4 for LF-HF Mixed, 5 for DY
    DiMuon_origin[0].Form("Charm");
    DiMuon_origin[1].Form("Beauty");
    DiMuon_origin[2].Form("HF_Mixed");
    DiMuon_origin[3].Form("LF");
    DiMuon_origin[4].Form("LF_HF_Mixed");
    DiMuon_origin[5].Form("DY");

    h_Pdg1Pdg2Pt_DiMuon_Gen = new TH3F("h_Pdg1Pdg2Pt_DiMuon_Gen", "; PDG code mum 1; PDG code mum 2; #it{p}_{T} (GeV/#it{c}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);
    h_Pdg1Pdg2Y_DiMuon_Gen = new TH3F("h_Pdg1Pdg2Y_DiMuon_Gen", "; PDG code mum 1; PDG code mum 2; #it{y}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_Y_Gen, low_bin_Y_Gen);
    h_Pdg1Pdg2M_DiMuon_Gen = new TH3F("h_Pdg1Pdg2M_DiMuon_Gen", "; PDG code mum 1; PDG code mum 2; #it{m}_{#mu#mu} (GeV/#it{c}^{2}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);

    h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut = new TH3F("h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut", "; PDG code mum 1; PDG code mum 2; #it{p}_{T} (GeV/#it{c}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);
    h_Pdg1Pdg2Y_DiMuon_Gen_DQcut = new TH3F("h_Pdg1Pdg2Y_DiMuon_Gen_DQcut", "; PDG code mum 1; PDG code mum 2; #it{y}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_Y_Rec, low_bin_Y_Rec);
    h_Pdg1Pdg2M_DiMuon_Gen_DQcut = new TH3F("h_Pdg1Pdg2M_DiMuon_Gen_DQcut", "; PDG code mum 1; PDG code mum 2; #it{m}_{#mu#mu} (GeV/#it{c}^{2}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);

    h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_Charm_corrected = new TH3F("h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_Charm_corrected", "; PDG code mum 1; PDG code mum 2; #it{p}_{T} (GeV/#it{c}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);
    h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_Charm_corrected = new TH3F("h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_Charm_corrected", "; PDG code mum 1; PDG code mum 2; #it{y}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_Y_Rec, low_bin_Y_Rec);
    h_Pdg1Pdg2M_DiMuon_Gen_DQcut_Charm_corrected = new TH3F("h_Pdg1Pdg2M_DiMuon_Gen_DQcut_Charm_corrected", "; PDG code mum 1; PDG code mum 2; #it{m}_{#mu#mu} (GeV/#it{c}^{2}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);

    h_Pdg1Pdg2Pt_DiMuon_Rec = new TH3F("h_Pdg1Pdg2Pt_DiMuon_Rec", "; PDG code mum 1; PDG code mum 2; #it{p}_{T} (GeV/#it{c}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);
    h_Pdg1Pdg2Y_DiMuon_Rec = new TH3F("h_Pdg1Pdg2Y_DiMuon_Rec", "; PDG code mum 1; PDG code mum 2; #it{y}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_Y_Rec, low_bin_Y_Rec);
    h_Pdg1Pdg2M_DiMuon_Rec = new TH3F("h_Pdg1Pdg2M_DiMuon_Rec", "; PDG code mum 1; PDG code mum 2; #it{m}_{#mu#mu} (GeV/#it{c}^{2}) ", n_PDG_selection - 1, PDG_Selection, n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);

    if (Generator.Contains("Powheg"))
    {
        h_Pdg1Pdg2Pt_DiMuon_Gen_PowhegOnly = (TH3F *)h_Pdg1Pdg2Pt_DiMuon_Gen->Clone(Form("%s_PowhegOnly", h_Pdg1Pdg2Pt_DiMuon_Gen->GetName()));
        h_Pdg1Pdg2Y_DiMuon_Gen_PowhegOnly = (TH3F *)h_Pdg1Pdg2Y_DiMuon_Gen->Clone(Form("%s_PowhegOnly", h_Pdg1Pdg2Y_DiMuon_Gen->GetName()));
        h_Pdg1Pdg2M_DiMuon_Gen_PowhegOnly = (TH3F *)h_Pdg1Pdg2M_DiMuon_Gen->Clone(Form("%s_PowhegOnly", h_Pdg1Pdg2M_DiMuon_Gen->GetName()));
        h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_PowhegOnly = (TH3F *)h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut->Clone(Form("%s_PowhegOnly", h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut->GetName()));
        h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_PowhegOnly = (TH3F *)h_Pdg1Pdg2Y_DiMuon_Gen_DQcut->Clone(Form("%s_PowhegOnly", h_Pdg1Pdg2Y_DiMuon_Gen_DQcut->GetName()));
        h_Pdg1Pdg2M_DiMuon_Gen_DQcut_PowhegOnly = (TH3F *)h_Pdg1Pdg2M_DiMuon_Gen_DQcut->Clone(Form("%s_PowhegOnly", h_Pdg1Pdg2M_DiMuon_Gen_DQcut->GetName()));
        h_Pdg1Pdg2Pt_DiMuon_Rec_PowhegOnly = (TH3F *)h_Pdg1Pdg2Pt_DiMuon_Rec->Clone(Form("%s_PowhegOnly", h_Pdg1Pdg2Pt_DiMuon_Rec->GetName()));
        h_Pdg1Pdg2Y_DiMuon_Rec_PowhegOnly = (TH3F *)h_Pdg1Pdg2Y_DiMuon_Rec->Clone(Form("%s_PowhegOnly", h_Pdg1Pdg2Y_DiMuon_Rec->GetName()));
        h_Pdg1Pdg2M_DiMuon_Rec_PowhegOnly = (TH3F *)h_Pdg1Pdg2M_DiMuon_Rec->Clone(Form("%s_PowhegOnly", h_Pdg1Pdg2M_DiMuon_Rec->GetName()));
    }
    // ---------------------- Inizialization hist for Dimuons from LF study, separating PYTHIA, GEANT and PYTHIA-GEANT components ---------------------------//

    Int_t n_LF_DiMuon_Generator = 999;

    if (Generator.Contains("Powheg"))
    {
        n_LF_DiMuon_Generator = 6;
        DiMuon_fromLF_Generator = new TString[n_LF_DiMuon_Generator]{"GeantOnly", "PYTHIAGeant", "PYTHIAOnly", "PowhegOnly", "PowhegGeant", "PowhegPYTHIA"};
    }
    else if (Generator.Contains("Geant"))
    {
        n_LF_DiMuon_Generator = 3;
        DiMuon_fromLF_Generator = new TString[n_LF_DiMuon_Generator]{"GeantOnly", "PYTHIAGeant", "PYTHIAOnly"};
    }
    else
        n_LF_DiMuon_Generator = 0;

    h_PtM_DiMuon_Rec_fromLF = new TH2F *[n_LF_DiMuon_Generator];
    h_PtY_DiMuon_Rec_fromLF = new TH2F *[n_LF_DiMuon_Generator];
    h_Nperevent_DiMuon_Rec_fromLF = new TH1F *[n_LF_DiMuon_Generator];

    h_PtM_DiMuon_Rec_fromLF_HF_Mixed = new TH2F *[n_LF_DiMuon_Generator];
    h_PtY_DiMuon_Rec_fromLF_HF_Mixed = new TH2F *[n_LF_DiMuon_Generator];
    h_Nperevent_DiMuon_Rec_fromLF_HF_Mixed = new TH1F *[n_LF_DiMuon_Generator];

    for (Int_t i_LF_DiMuon_Generator = 0; i_LF_DiMuon_Generator < n_LF_DiMuon_Generator; i_LF_DiMuon_Generator++)
    {
        h_PtM_DiMuon_Rec_fromLF[i_LF_DiMuon_Generator] = new TH2F(Form("h_PtM_DiMuon_Rec_fromLF_%s", DiMuon_fromLF_Generator[i_LF_DiMuon_Generator].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 400, 0, 40.0, 400, 0, 40);
        h_PtY_DiMuon_Rec_fromLF[i_LF_DiMuon_Generator] = new TH2F(Form("h_PtY_DiMuon_Rec_fromLF_%s", DiMuon_fromLF_Generator[i_LF_DiMuon_Generator].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{y}_{#mu#mu}", 400, 0, 40.0, 150, -4.0, -2.5);
        h_Nperevent_DiMuon_Rec_fromLF[i_LF_DiMuon_Generator] = new TH1F(Form("h_Nperevent_DiMuon_Rec_fromLF_%s", DiMuon_fromLF_Generator[i_LF_DiMuon_Generator].Data()), "; #mu#mu x ev", 10, -0.5, 9.5);

        h_PtM_DiMuon_Rec_fromLF_HF_Mixed[i_LF_DiMuon_Generator] = new TH2F(Form("h_PtM_DiMuon_Rec_fromLF_HF_Mixed_%s", DiMuon_fromLF_Generator[i_LF_DiMuon_Generator].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 400, 0, 40.0, 400, 0, 40);
        h_PtY_DiMuon_Rec_fromLF_HF_Mixed[i_LF_DiMuon_Generator] = new TH2F(Form("h_PtY_DiMuon_Rec_fromLF_HF_Mixed_%s", DiMuon_fromLF_Generator[i_LF_DiMuon_Generator].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{y}_{#mu#mu}", 400, 0, 40.0, 150, -4.0, -2.5);
        h_Nperevent_DiMuon_Rec_fromLF_HF_Mixed[i_LF_DiMuon_Generator] = new TH1F(Form("h_Nperevent_DiMuon_Rec_fromLF_HF_Mixed_%s", DiMuon_fromLF_Generator[i_LF_DiMuon_Generator].Data()), "; #mu#mu x ev", 10, -0.5, 9.5);
    }

    // ---------------------- Inizialization hist for Dimuons ---------------------------//

    for (Int_t i_DiMuon_origin = 0; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
    {
        // Inizialization hist for Generated Dimuons
        h_PtM_DiMuon_Gen[i_DiMuon_origin] = new TH2F(Form("h_PtM_DiMuon_Gen_%s", DiMuon_origin[i_DiMuon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 400, 0, 40.0, 400, 0, 40);
        h_PtY_DiMuon_Gen[i_DiMuon_origin] = new TH2F(Form("h_PtY_DiMuon_Gen_%s", DiMuon_origin[i_DiMuon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{y}_{#mu#mu}", 400, 0, 40.0, 160, -8.0, 8.0);
        h_Nperevent_DiMuon_Gen[i_DiMuon_origin] = new TH1F(Form("h_Nperevent_DiMuon_Gen_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; #mu#mu x ev", 10, -0.5, 9.5);
        // Inizialization hist for Generated Dimuons with DQ cuts
        h_PtM_DiMuon_Gen_DQcut[i_DiMuon_origin] = new TH2F(Form("h_PtM_DiMuon_Gen_DQcut_%s", DiMuon_origin[i_DiMuon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 400, 0, 40.0, 400, 0, 40);
        h_PtY_DiMuon_Gen_DQcut[i_DiMuon_origin] = new TH2F(Form("h_PtY_DiMuon_Gen_DQcut_%s", DiMuon_origin[i_DiMuon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 400, 0, 40.0, 150, -4.0, -2.5);
        h_Nperevent_DiMuon_Gen_DQcut[i_DiMuon_origin] = new TH1F(Form("h_Nperevent_DiMuon_Gen_DQcut_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; #mu#mu x ev", 10, -0.5, 9.5);
        // Inizialization hist for Reconstructed Dimuons
        h_PtM_DiMuon_Rec[i_DiMuon_origin] = new TH2F(Form("h_PtM_DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 400, 0, 40.0, 400, 0, 40);
        h_PtY_DiMuon_Rec[i_DiMuon_origin] = new TH2F(Form("h_PtY_DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 400, 0, 40.0, 150, -4.0, -2.5);
        h_Nperevent_DiMuon_Rec[i_DiMuon_origin] = new TH1F(Form("h_Nperevent_DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), "; #mu#mu x ev", 10, -0.5, 9.5);

        if (Generator.Contains("Powheg"))
        {
            h_PtM_DiMuon_Gen_PowhegOnly[i_DiMuon_origin] = (TH2F *)h_PtM_DiMuon_Gen[i_DiMuon_origin]->Clone(Form("%s_PowhegOnly", h_PtM_DiMuon_Gen[i_DiMuon_origin]->GetName()));
            h_PtY_DiMuon_Gen_PowhegOnly[i_DiMuon_origin] = (TH2F *)h_PtY_DiMuon_Gen[i_DiMuon_origin]->Clone(Form("%s_PowhegOnly", h_PtY_DiMuon_Gen[i_DiMuon_origin]->GetName()));
            h_Nperevent_DiMuon_Gen_PowhegOnly[i_DiMuon_origin] = (TH1F *)h_Nperevent_DiMuon_Gen[i_DiMuon_origin]->Clone(Form("%s_PowhegOnly", h_Nperevent_DiMuon_Gen[i_DiMuon_origin]->GetName()));
            h_PtM_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin] = (TH2F *)h_PtM_DiMuon_Gen_DQcut[i_DiMuon_origin]->Clone(Form("%s_PowhegOnly", h_PtM_DiMuon_Gen_DQcut[i_DiMuon_origin]->GetName()));
            h_PtY_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin] = (TH2F *)h_PtY_DiMuon_Gen_DQcut[i_DiMuon_origin]->Clone(Form("%s_PowhegOnly", h_PtY_DiMuon_Gen_DQcut[i_DiMuon_origin]->GetName()));
            h_Nperevent_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin] = (TH1F *)h_Nperevent_DiMuon_Gen_DQcut[i_DiMuon_origin]->Clone(Form("%s_PowhegOnly", h_Nperevent_DiMuon_Gen_DQcut[i_DiMuon_origin]->GetName()));
            h_PtM_DiMuon_Rec_PowhegOnly[i_DiMuon_origin] = (TH2F *)h_PtM_DiMuon_Rec[i_DiMuon_origin]->Clone(Form("%s_PowhegOnly", h_PtM_DiMuon_Rec[i_DiMuon_origin]->GetName()));
            h_PtY_DiMuon_Rec_PowhegOnly[i_DiMuon_origin] = (TH2F *)h_PtY_DiMuon_Rec[i_DiMuon_origin]->Clone(Form("%s_PowhegOnly", h_PtY_DiMuon_Rec[i_DiMuon_origin]->GetName()));
            h_Nperevent_DiMuon_Rec_PowhegOnly[i_DiMuon_origin] = (TH1F *)h_Nperevent_DiMuon_Rec[i_DiMuon_origin]->Clone(Form("%s_PowhegOnly", h_Nperevent_DiMuon_Rec[i_DiMuon_origin]->GetName()));
        }
    }

    //------ Inizialization hist for HF hadrons----------------//

    h_PdgPtY_HFHadron_prompt = new TH3F("h_PdgPtY_HFHadron_prompt", "; PDG hadron; #it{p}_{T} (GeV/#it{c}) ; #it{y}", n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt, n_bin_Y_Gen, low_bin_Y_Gen);
    h_PdgPtY_HFHadron_notprompt = new TH3F("h_PdgPtY_HFHadron_notprompt", "; PDG hadron; #it{p}_{T} (GeV/#it{c}) ; #it{y}", n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt, n_bin_Y_Gen, low_bin_Y_Gen);

    h_PdgY_HFHadron_notprompt = new TH2F("h_PdgY_HFHadron_notprompt", "; PDG hadron; #it{y}", n_PDG_selection - 1, PDG_Selection, n_bin_Y_Gen, low_bin_Y_Gen);
    h_PdgY_HFHadron_prompt = new TH2F("h_PdgY_HFHadron_prompt", "; PDG hadron; #it{y}", n_PDG_selection - 1, PDG_Selection, n_bin_Y_Gen, low_bin_Y_Gen);

    h_PdgEta_HFHadron_notprompt = new TH2F("h_PdgEta_HFHadron_notprompt", "; PDG hadron; #eta", n_PDG_selection - 1, PDG_Selection, n_bin_Y_Gen, low_bin_Y_Gen);
    h_PdgEta_HFHadron_prompt = new TH2F("h_PdgEta_HFHadron_prompt", "; PDG hadron; #eta", n_PDG_selection - 1, PDG_Selection, n_bin_Y_Gen, low_bin_Y_Gen);

    h_PdgPt_HFHadron_prompt = new TH2F("h_PdgPt_HFHadron_prompt", "; PDG hadron; #it{p}_{T} (GeV/#it{c})", n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);
    h_PdgPt_HFHadron_notprompt = new TH2F("h_PdgPt_HFHadron_notprompt", "; PDG hadron; #it{p}_{T} (GeV/#it{c})", n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);

    if (Generator.Contains("Geant"))
    {

        h_PdgY_HFHadron_notprompt_PYTHIAOnly = (TH2F *)h_PdgY_HFHadron_notprompt->Clone(Form("%s_PYTHIAOnly", h_PdgY_HFHadron_notprompt->GetName()));
        h_PdgY_HFHadron_prompt_PYTHIAOnly = (TH2F *)h_PdgY_HFHadron_prompt->Clone(Form("%s_PYTHIAOnly", h_PdgY_HFHadron_prompt->GetName()));

        h_PdgEta_HFHadron_notprompt_PYTHIAOnly = (TH2F *)h_PdgEta_HFHadron_notprompt->Clone(Form("%s_PYTHIAOnly", h_PdgEta_HFHadron_notprompt->GetName()));
        h_PdgEta_HFHadron_prompt_PYTHIAOnly = (TH2F *)h_PdgEta_HFHadron_prompt->Clone(Form("%s_PYTHIAOnly", h_PdgEta_HFHadron_prompt->GetName()));

        h_PdgPt_HFHadron_notprompt_PYTHIAOnly = (TH2F *)h_PdgPt_HFHadron_notprompt->Clone(Form("%s_PYTHIAOnly", h_PdgPt_HFHadron_notprompt->GetName()));
        h_PdgPt_HFHadron_prompt_PYTHIAOnly = (TH2F *)h_PdgPt_HFHadron_prompt->Clone(Form("%s_PYTHIAOnly", h_PdgPt_HFHadron_prompt->GetName()));

        h_PdgY_HFHadron_notprompt_GeantOnly = (TH2F *)h_PdgY_HFHadron_notprompt->Clone(Form("%s_GeantOnly", h_PdgY_HFHadron_notprompt->GetName()));
        h_PdgY_HFHadron_prompt_GeantOnly = (TH2F *)h_PdgY_HFHadron_prompt->Clone(Form("%s_GeantOnly", h_PdgY_HFHadron_prompt->GetName()));

        h_PdgEta_HFHadron_notprompt_GeantOnly = (TH2F *)h_PdgEta_HFHadron_notprompt->Clone(Form("%s_GeantOnly", h_PdgEta_HFHadron_notprompt->GetName()));
        h_PdgEta_HFHadron_prompt_GeantOnly = (TH2F *)h_PdgEta_HFHadron_prompt->Clone(Form("%s_GeantOnly", h_PdgEta_HFHadron_prompt->GetName()));

        h_PdgPt_HFHadron_notprompt_GeantOnly = (TH2F *)h_PdgPt_HFHadron_notprompt->Clone(Form("%s_GeantOnly", h_PdgPt_HFHadron_notprompt->GetName()));
        h_PdgPt_HFHadron_prompt_GeantOnly = (TH2F *)h_PdgPt_HFHadron_prompt->Clone(Form("%s_GeantOnly", h_PdgPt_HFHadron_prompt->GetName()));

        if (Generator.Contains("Powheg"))
        {
            h_PdgY_HFHadron_notprompt_PowhegOnly = (TH2F *)h_PdgY_HFHadron_notprompt->Clone(Form("%s_PowhegOnly", h_PdgY_HFHadron_notprompt->GetName()));
            h_PdgY_HFHadron_prompt_PowhegOnly = (TH2F *)h_PdgY_HFHadron_prompt->Clone(Form("%s_PowhegOnly", h_PdgY_HFHadron_prompt->GetName()));

            h_PdgEta_HFHadron_notprompt_PowhegOnly = (TH2F *)h_PdgEta_HFHadron_notprompt->Clone(Form("%s_PowhegOnly", h_PdgEta_HFHadron_notprompt->GetName()));
            h_PdgEta_HFHadron_prompt_PowhegOnly = (TH2F *)h_PdgEta_HFHadron_prompt->Clone(Form("%s_PowhegOnly", h_PdgEta_HFHadron_prompt->GetName()));

            h_PdgPt_HFHadron_notprompt_PowhegOnly = (TH2F *)h_PdgPt_HFHadron_notprompt->Clone(Form("%s_PowhegOnly", h_PdgPt_HFHadron_notprompt->GetName()));
            h_PdgPt_HFHadron_prompt_PowhegOnly = (TH2F *)h_PdgPt_HFHadron_prompt->Clone(Form("%s_PowhegOnly", h_PdgPt_HFHadron_prompt->GetName()));
        }

        // h_PdgEta_HFHadron_notprompt = new TH2F("h_PdgEta_HFHadron_notprompt", "; PDG hadron; #eta", n_PDG_selection - 1, PDG_Selection, n_bin_Y_Gen, low_bin_Y_Gen);
        // h_PdgEta_HFHadron_prompt = new TH2F("h_PdgEta_HFHadron_prompt", "; PDG hadron; #eta", n_PDG_selection - 1, PDG_Selection, n_bin_Y_Gen, low_bin_Y_Gen);

        // h_PdgPt_HFHadron_prompt = new TH2F("h_PdgPt_HFHadron_prompt", "; PDG hadron; #it{p}_{T} (GeV/#it{c})", n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);
        // h_PdgPt_HFHadron_notprompt = new TH2F("h_PdgPt_HFHadron_notprompt", "; PDG hadron; #it{p}_{T} (GeV/#it{c})", n_PDG_selection - 1, PDG_Selection, n_bin_pt, low_bin_pt);
    }
    // if (forZ_sim)
    // {
    //     Double_t low_mass = 0;
    //     Double_t high_mass = 120;
    //     Int_t N_Bin = 1200;
    //     h_PtM_DiMuon_Gen[limit_DiMu_origin] = new TH2F(Form("h_PtM_DiMuon_Gen_%s", DiMuon_origin[limit_DiMu_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Gen[limit_DiMu_origin] = new TH2F(Form("h_PtY_DiMuon_Gen_%s", DiMuon_origin[limit_DiMu_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{y}_{#mu#mu}", 1200, 0, 120.0, 160, -8.0, 8.0);
    //     h_Nperevent_DiMuon_Gen[limit_DiMu_origin] = new TH1F(Form("h_Nperevent_DiMuon_Gen_%s", DiMuon_origin[limit_DiMu_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);
    //     // Inizialization hist for Generated Dimuons with DQ cuts
    //     h_PtM_DiMuon_Gen_DQcut[limit_DiMu_origin] = new TH2F(Form("h_PtM_DiMuon_Gen_DQcut_%s", DiMuon_origin[limit_DiMu_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Gen_DQcut[limit_DiMu_origin] = new TH2F(Form("h_PtY_DiMuon_Gen_DQcut_%s", DiMuon_origin[limit_DiMu_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 150, -4.0, -2.5);
    //     h_Nperevent_DiMuon_Gen_DQcut[limit_DiMu_origin] = new TH1F(Form("h_Nperevent_DiMuon_Gen_DQcut_%s", DiMuon_origin[limit_DiMu_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);
    //     // Inizialization hist for Reconstructed Dimuons
    //     h_PtM_DiMuon_Rec[limit_DiMu_origin] = new TH2F(Form("h_PtM_DiMuon_Rec_%s", DiMuon_origin[limit_DiMu_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Rec[limit_DiMu_origin] = new TH2F(Form("h_PtY_DiMuon_Rec_%s", DiMuon_origin[limit_DiMu_origin].Data()), ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 150, -4.0, -2.5);
    //     h_Nperevent_DiMuon_Rec[limit_DiMu_origin] = new TH1F(Form("h_Nperevent_DiMuon_Rec_%s", DiMuon_origin[limit_DiMu_origin].Data()), "; #mu x ev", 10, -0.5, 9.5);

    //     h_PtM_DiMuon_Gen_Z_ptmucut09 = new TH2F("h_PtM_DiMuon_Gen_Z_ptmucut09", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Gen_Z_ptmucut09 = new TH2F("h_PtY_DiMuon_Gen_Z_ptmucut09", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 160, -8.0, 8.0);
    //     h_Nperevent_DiMuon_Gen_Z_ptmucut09 = new TH1F("h_Nperevent_DiMuon_Gen_Z_ptmucut09", "; #mu x ev", 10, -0.5, 9.5);

    //     h_PtM_DiMuon_Gen_Z_ptmucut10 = new TH2F("h_PtM_DiMuon_Gen_Z_ptmucut10", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Gen_Z_ptmucut10 = new TH2F("h_PtY_DiMuon_Gen_Z_ptmucut10", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 160, -8.0, 8.0);
    //     h_Nperevent_DiMuon_Gen_Z_ptmucut10 = new TH1F("h_Nperevent_DiMuon_Gen_Z_ptmucut10", "; #mu x ev", 10, -0.5, 9.5);

    //     h_PtM_DiMuon_Gen_Z_ptmucut20 = new TH2F("h_PtM_DiMuon_Gen_Z_ptmucut20", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Gen_Z_ptmucut20 = new TH2F("h_PtY_DiMuon_Gen_Z_ptmucut20", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 160, -8.0, 8.0);
    //     h_Nperevent_DiMuon_Gen_Z_ptmucut20 = new TH1F("h_Nperevent_DiMuon_Gen_Z_ptmucut20", "; #mu x ev", 10, -0.5, 9.5);

    //     h_PtM_DiMuon_Gen_Z_DQcut_ptmucut09 = new TH2F("h_PtM_DiMuon_Gen_Z_DQcut_ptmucut09", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Gen_Z_DQcut_ptmucut09 = new TH2F("h_PtY_DiMuon_Gen_Z_DQcut_ptmucut09", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 150, -4.0, -2.5);
    //     h_Nperevent_DiMuon_Gen_Z_DQcut_ptmucut09 = new TH1F("h_Nperevent_DiMuon_Gen_Z_DQcut_ptmucut09", "; #mu x ev", 10, -0.5, 9.5);

    //     h_PtM_DiMuon_Gen_Z_DQcut_ptmucut10 = new TH2F("h_PtM_DiMuon_Gen_Z_DQcut_ptmucut10", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Gen_Z_DQcut_ptmucut10 = new TH2F("h_PtY_DiMuon_Gen_Z_DQcut_ptmucut10", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 150, -4.0, -2.5);
    //     h_Nperevent_DiMuon_Gen_Z_DQcut_ptmucut10 = new TH1F("h_Nperevent_DiMuon_Gen_Z_DQcut_ptmucut10", "; #mu x ev", 10, -0.5, 9.5);

    //     h_PtM_DiMuon_Gen_Z_DQcut_ptmucut20 = new TH2F("h_PtM_DiMuon_Gen_Z_DQcut_ptmucut20", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Gen_Z_DQcut_ptmucut20 = new TH2F("h_PtY_DiMuon_Gen_Z_DQcut_ptmucut20", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 150, -4.0, -2.5);
    //     h_Nperevent_DiMuon_Gen_Z_DQcut_ptmucut20 = new TH1F("h_Nperevent_DiMuon_Gen_Z_DQcut_ptmucut20", "; #mu x ev", 10, -0.5, 9.5);

    //     h_PtM_DiMuon_Rec_Z_ptmucut09 = new TH2F("h_PtM_DiMuon_Rec_Z_ptmucut09", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Rec_Z_ptmucut09 = new TH2F("h_PtY_DiMuon_Rec_Z_ptmucut09", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 150, -4.0, -2.5);
    //     h_Nperevent_DiMuon_Rec_Z_ptmucut09 = new TH1F("h_Nperevent_DiMuon_Rec_Z_ptmucut09", "; #mu x ev", 10, -0.5, 9.5);

    //     h_PtM_DiMuon_Rec_Z_ptmucut10 = new TH2F("h_PtM_DiMuon_Rec_Z_ptmucut10", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Rec_Z_ptmucut10 = new TH2F("h_PtY_DiMuon_Rec_Z_ptmucut10", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 150, -4.0, -2.5);
    //     h_Nperevent_DiMuon_Rec_Z_ptmucut10 = new TH1F("h_Nperevent_DiMuon_Rec_Z_ptmucut10", "; #mu x ev", 10, -0.5, 9.5);

    //     h_PtM_DiMuon_Rec_Z_ptmucut20 = new TH2F("h_PtM_DiMuon_Rec_Z_ptmucut20", ";#it{p}_{T} (GeV/#it{c}) ; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", 1200, 0, 120.0, N_Bin, low_mass, high_mass);
    //     h_PtY_DiMuon_Rec_Z_ptmucut20 = new TH2F("h_PtY_DiMuon_Rec_Z_ptmucut20", ";#it{p}_{T} (GeV/#it{c}) ; #it{y}", 1200, 0, 120.0, 150, -4.0, -2.5);
    //     h_Nperevent_DiMuon_Rec_Z_ptmucut20 = new TH1F("h_Nperevent_DiMuon_Rec_Z_ptmucut20", "; #mu x ev", 10, -0.5, 9.5);
    // }
}

TChain *Importing_Tree(TString dir_filename, TString filename, TString Generator)
{
    TChain *tree = nullptr;

    tree = new TChain("MCTree");
    //    printf("%s \n",Form("%s/%s",dir_filename.Data(),filename.Data()));
    tree->AddFile(Form("%s/%s", dir_filename.Data(), filename.Data()));
    if (Generator.Contains("HF"))
    {
        tree->SetBranchAddress("N_HFquarks_gen", &N_HFquarks_gen);
        tree->SetBranchAddress("PDG_HFquark_gen", PDG_HFquark_gen);
        tree->SetBranchAddress("Pt_HFquark_gen", Pt_HFquark_gen);
        tree->SetBranchAddress("Y_HFquark_gen", Y_HFquark_gen);
        if (Generator.Contains("Powheg"))
            tree->SetBranchAddress("Mother_index", Mother_index);
        tree->SetBranchAddress("NHadrons_gen", &NHadrons_gen);
        tree->SetBranchAddress("PDGmum_Hadron_gen", PDGmum_Hadron_gen);
        tree->SetBranchAddress("PDG_Hadron_gen", PDG_Hadron_gen);
        tree->SetBranchAddress("Pt_Hadron_gen", Pt_Hadron_gen);
        tree->SetBranchAddress("E_Hadron_gen", E_Hadron_gen);
        tree->SetBranchAddress("Px_Hadron_gen", Px_Hadron_gen);
        tree->SetBranchAddress("Py_Hadron_gen", Py_Hadron_gen);
        tree->SetBranchAddress("Pz_Hadron_gen", Pz_Hadron_gen);
        tree->SetBranchAddress("Y_Hadron_gen", Y_Hadron_gen);
        tree->SetBranchAddress("Eta_Hadron_gen", Eta_Hadron_gen);
        if (Generator.Contains("Powheg"))
            tree->SetBranchAddress("HadronFrom_Powheg_gen", fHadronFrom_Powheg_gen);
    }
    if (Generator.Contains("DY"))
    {
        tree->SetBranchAddress("N_gamma_gen", &fN_gamma);
        tree->SetBranchAddress("Pt_gamma_gen", fPt_gamma);
        tree->SetBranchAddress("M_gamma_gen", fM_gamma);
        tree->SetBranchAddress("Y_gamma_gen", fY_gamma);
    }
    if (Generator.Contains("Geant"))
    {
        tree->SetBranchAddress("HadronFrom_Geant_gen", fHadronFrom_Geant_gen);
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
        tree->SetBranchAddress("From_Powheg_rec", fFrom_Powheg_rec);
        tree->SetBranchAddress("Initial_Parton_rec", fInitial_Parton_rec);
        tree->SetBranchAddress("Vzmother_rec", fVzmother_rec);
        tree->SetBranchAddress("From_Geant_rec", fFrom_Geant_rec);

        tree->SetBranchAddress("VzHadron_gen", fVzHadron_gen);
        tree->SetBranchAddress("From_Geant_gen", fFrom_Geant_gen);
        tree->SetBranchAddress("Vzmother_gen", fVzmother_gen);
        tree->SetBranchAddress("Radius_gen", fRadius_gen);
        tree->SetBranchAddress("Vz_gen", fVz_gen);

        if (Generator.Contains("Powheg"))
        {
            tree->SetBranchAddress("From_Powheg_gen", fFrom_Powheg_gen);
            tree->SetBranchAddress("Initial_Parton_gen", fInitial_Parton_gen);
        }
    }

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