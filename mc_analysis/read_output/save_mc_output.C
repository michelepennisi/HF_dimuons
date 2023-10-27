#include "save_mc_output.h"

void save_mc_output(
    // TString RunMode = "test_new_prompt_LHC22b3",
    TString RunMode = "SoftQCD_Def",
    TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/sim",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LF_test",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/test_beauty_sim_2",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/test_charm_sim",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC18p_DY_100k_Version2_AOD",
    TString dir_fileOut = "test_new_sim",
    Int_t RunNumber = 200,
    TString Generator = "HF_Pythia",
    TString prefix_filename = "pythia_sim_12345_DefaultBR")
{

    Set_Histograms(Generator);

    TString filename;
    filename.Form("%s_%s_%d.root", RunMode.Data(), prefix_filename.Data(), RunNumber);

    TString file_out;
    file_out.Form("%s_MC_output_Hist_%d.root", RunMode.Data(), RunNumber);
    // file_out.Form("old_%s_MC_output_Hist_%d.root", RunMode.Data(), RunNumber);

    TString file_out_tree;
    file_out_tree.Form("%s_MC_output_Tree_%d.root", RunMode.Data(), RunNumber);

    printf("%s/%s\n", dir_fileOut.Data(), file_out.Data());
    printf("%s/%s\n", dir_fileOut.Data(), file_out_tree.Data());

    Double_t *Pt_Dimu_Rec = new Double_t;
    Double_t *M_Dimu_Rec = new Double_t;

    Double_t *Pt_Dimu_Rec_PowhegOnly = new Double_t;
    Double_t *M_Dimu_Rec_PowhegOnly = new Double_t;

    TTree *Tree_DiMuon_Rec[n_DiMuon_origin];
    TTree *Tree_DiMuon_Rec_PowhegOnly[n_DiMuon_origin];
    for (Int_t i_DiMuon_origin = 0; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
    {
        Tree_DiMuon_Rec[i_DiMuon_origin] = new TTree(Form("Tree_DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
        Tree_DiMuon_Rec[i_DiMuon_origin]->Branch("m", M_Dimu_Rec, "m/D");
        Tree_DiMuon_Rec[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec, "pt/D");

        if (Generator.Contains("Powheg"))
        {
            Tree_DiMuon_Rec_PowhegOnly[i_DiMuon_origin] = new TTree(Form("Tree_DiMuon_Rec_PowhegOnly_%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
            Tree_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_PowhegOnly, "m/D");
            Tree_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_PowhegOnly, "pt/D");
        }
    }

    printf("Input File: %s\n", filename.Data());
    printf("Saving in dir: %s \nFile: %s\n", dir_fileOut.Data(), file_out.Data());

    TChain *input_tree = Importing_Tree(dir_fileIn, filename, Generator);
    input_tree->ls();

    Int_t total_entries = input_tree->GetEntries();

    if (total_entries == 0)
        return;

    Bool_t Verbose = kFALSE;
    Int_t counter_test = 0;

    for (Int_t i_Event = 0; i_Event < total_entries; i_Event++)
    {
        h_Nevents->Fill(1);
        if (i_Event % 250000 == 0)
        {
            progress_status(i_Event, total_entries);
        }
        input_tree->GetEntry(i_Event);
        Int_t n_Muon_Gen[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_PowhegOnly[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_DQcut[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_DQcut_PowhegOnly[n_Muon_origin] = {0};
        Int_t n_Muon_Rec[n_Muon_origin] = {0};
        Int_t n_Muon_Rec_PowhegOnly[n_Muon_origin] = {0};

        Int_t n_DiMuon_Gen[n_DiMuon_origin] = {0};
        Int_t n_DiMuon_Gen_DQcut[n_DiMuon_origin] = {0};
        Int_t n_DiMuon_Rec[n_DiMuon_origin] = {0};

        Int_t n_DiMuon_Gen_PowhegOnly[n_DiMuon_origin] = {0};
        Int_t n_DiMuon_Gen_DQcut_PowhegOnly[n_DiMuon_origin] = {0};
        Int_t n_DiMuon_Rec_PowhegOnly[n_DiMuon_origin] = {0};

        Double_t Pt_Gamma = fPt_gamma[0];
        Double_t Y_Gamma = fY_gamma[0];
        Double_t M_Gamma = fM_gamma[0];

        h_PtM_Gamma->Fill(Pt_Gamma, M_Gamma);
        h_PtY_Gamma->Fill(Pt_Gamma, Y_Gamma);

        for (Int_t i_NHadron_gen = 0; i_NHadron_gen < NHadrons_gen; i_NHadron_gen++)
        {
            Int_t PDGmum_Hadron = PDGmum_Hadron_gen[i_NHadron_gen];           // gen Hadron PDG mum
            Int_t PDG_Hadron = PDG_Hadron_gen[i_NHadron_gen];                 // gen Hadron PDG
            Double_t Pt_Hadron = Pt_Hadron_gen[i_NHadron_gen];                // gen Hadron pT
            Double_t E_Hadron = E_Hadron_gen[i_NHadron_gen];                  // gen Hadron E
            Double_t Px_Hadron = Px_Hadron_gen[i_NHadron_gen];                // gen Hadron px
            Double_t Py_Hadron = Py_Hadron_gen[i_NHadron_gen];                // gen Hadron py
            Double_t Pz_Hadron = Pz_Hadron_gen[i_NHadron_gen];                // gen Hadron pz
            Double_t Y_Hadron = Y_Hadron_gen[i_NHadron_gen];                  // gen Hadron y
            Double_t Eta_Hadron = Eta_Hadron_gen[i_NHadron_gen];              // gen Hadron eta
            Int_t IsHadronfromPowheg = fHadronFrom_Powheg_gen[i_NHadron_gen]; // check muon gen origin

            if ((TMath::Abs(PDG_Hadron) > 400 && TMath::Abs(PDG_Hadron) < 500) || (TMath::Abs(PDG_Hadron) > 4000 && TMath::Abs(PDG_Hadron) < 5000))
            {
                if ((PDGmum_Hadron == 0) || (PDGmum_Hadron == 4) || (TMath::Abs(PDGmum_Hadron) > 400 && TMath::Abs(PDGmum_Hadron) < 500) || (TMath::Abs(PDGmum_Hadron) > 4000 && TMath::Abs(PDGmum_Hadron) < 5000))
                {
                    h_PdgPt_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                    h_PdgY_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);

                    h_PdgPtY_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron, Y_Hadron);
                    if (Generator.Contains("Powheg") && IsHadronfromPowheg == -1)
                        h_PdgPtY_HFHadron_prompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron, Y_Hadron);
                }
                else
                {
                    h_PdgPtY_HFHadron_notprompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron, Y_Hadron);
                    h_PdgPt_HFHadron_notprompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                    h_PdgY_HFHadron_notprompt->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                    if (Generator.Contains("Powheg") && IsHadronfromPowheg == -1)
                        h_PdgPtY_HFHadron_notprompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron, Y_Hadron);
                }
            }
            else
            {
                h_PdgPtY_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron, Y_Hadron);
                h_PdgPt_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                h_PdgY_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                if (Generator.Contains("Powheg") && IsHadronfromPowheg == -1)
                    h_PdgPtY_HFHadron_prompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron, Y_Hadron);
            }

            // printf("i) %d | PDG %d | PDG mum %d | Pt %f | Y %f\n", i_NHadron_gen, PDG_Hadron,PDGmum_Hadron, Pt_Hadron, Y_Hadron);
        }
        for (Int_t i_NMuons_gen = 0; i_NMuons_gen < NMuons_gen; i_NMuons_gen++)
        {

            Int_t PDG_Mu = PDGmum_gen[i_NMuons_gen];                  // single gen mu PDG mum
            Double_t Pt_Mu = Pt_gen[i_NMuons_gen];                    // single gen mu pT
            Double_t E_Mu = E_gen[i_NMuons_gen];                      // single gen mu E
            Double_t Px_Mu = Px_gen[i_NMuons_gen];                    // single gen mu px
            Double_t Py_Mu = Py_gen[i_NMuons_gen];                    // single gen mu py
            Double_t Pz_Mu = Pz_gen[i_NMuons_gen];                    // single gen mu pz
            Double_t Y_Mu = Y_gen[i_NMuons_gen];                      // single gen mu y
            Double_t Eta_Mu = Eta_gen[i_NMuons_gen];                  // single gen mu eta
            Double_t Phi_Mu = Phi_gen[i_NMuons_gen];                  // single gen mu phi
            Double_t Theta_Mu = Theta_gen[i_NMuons_gen];              // single gen mu theta
            Double_t Charge_Mu = Charge_gen[i_NMuons_gen];            // single gen mu theta
            Int_t IsFrom_Powheg_gen = fFrom_Powheg_gen[i_NMuons_gen]; // single gen mu theta
            Int_t Initial_parton = fInitial_Parton_gen[i_NMuons_gen]; // single gen mu theta

            Bool_t DQ_Muon = kFALSE;
            Bool_t Selection_Muon[n_Muon_origin] = {kFALSE};

            Selection_Muon[0] = kTRUE; // Selection All Muons
            if ((TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500) || (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000))
            {
                Selection_Muon[1] = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600) || (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000))
            {
                Selection_Muon[2] = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu) > 100 && TMath::Abs(PDG_Mu) < 400) || (TMath::Abs(PDG_Mu) > 1000 && TMath::Abs(PDG_Mu) < 4000))
            {
                Selection_Muon[3] = kTRUE; // Selection Muons from LF
            }
            else if (TMath::Abs(PDG_Mu) == 23)
            {
                Selection_Muon[4] = kTRUE;
            }

            // printf("Initial_parton %d\n",Initial_parton);
            if (Selection_Muon[3] && Generator.Contains("Powheg"))
            {
                if (IsFrom_Powheg_gen == -1)
                {
                    h_PtQ_Muon_Gen_LF_Powheg->Fill(Pt_Mu, Initial_parton);
                    h_YQ_Muon_Gen_LF_Powheg->Fill(Y_Mu, Initial_parton);
                }
                else
                {
                    h_PtQ_Muon_Gen_LF_Pythia->Fill(Pt_Mu, Initial_parton);
                    h_YQ_Muon_Gen_LF_Pythia->Fill(Y_Mu, Initial_parton);
                }
            }

            for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
            {
                if (Selection_Muon[i_Muon_origin])
                {
                    h_PtY_Muon_Gen[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    n_Muon_Gen[i_Muon_origin]++;
                    if (Generator.Contains("Powheg") && IsFrom_Powheg_gen == -1)
                    {
                        h_PtY_Muon_Gen_PowhegOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                        n_Muon_Gen_PowhegOnly[i_Muon_origin]++;
                    }
                }
            }

            h_PtPdg_Muon_Gen->Fill(Pt_Mu, PDG_Mu);
            h_YPdg_Muon_Gen->Fill(Y_Mu, PDG_Mu);
            if (Generator.Contains("Powheg") && IsFrom_Powheg_gen == -1)
            {
                h_PtPdg_Muon_Gen_PowhegOnly->Fill(Pt_Mu, PDG_Mu);
                h_YPdg_Muon_Gen_PowhegOnly->Fill(Y_Mu, PDG_Mu);
            }

            if (Eta_Mu > -4.0 && Eta_Mu < -2.5)
                DQ_Muon = kTRUE;
            if (!DQ_Muon)
                continue;

            h_PtPdg_Muon_Gen_DQcut->Fill(Pt_Mu, PDG_Mu);
            h_YPdg_Muon_Gen_DQcut->Fill(Y_Mu, PDG_Mu);

            if (Generator.Contains("Powheg") && IsFrom_Powheg_gen == -1)
            {
                h_PtPdg_Muon_Gen_DQcut_PowhegOnly->Fill(Pt_Mu, PDG_Mu);
                h_YPdg_Muon_Gen_DQcut_PowhegOnly->Fill(Y_Mu, PDG_Mu);
            }

            for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
            {
                if (Selection_Muon[i_Muon_origin])
                {
                    h_PtY_Muon_Gen_DQcut[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    n_Muon_Gen_DQcut[i_Muon_origin]++;
                    if (Generator.Contains("Powheg") && IsFrom_Powheg_gen == -1)
                    {
                        n_Muon_Gen_DQcut_PowhegOnly[i_Muon_origin]++;
                        h_PtY_Muon_Gen_DQcut_PowhegOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    }
                }
            }
        }

        for (Int_t i_Muon_Origin = 0; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
        {
            h_Nperevent_Muon_Gen[i_Muon_Origin]->Fill(n_Muon_Gen[i_Muon_Origin]);
            h_Nperevent_Muon_Gen_DQcut[i_Muon_Origin]->Fill(n_Muon_Gen_DQcut[i_Muon_Origin]);
            if (Generator.Contains("Powheg"))
            {
                h_Nperevent_Muon_Gen_PowhegOnly[i_Muon_Origin]->Fill(n_Muon_Gen_PowhegOnly[i_Muon_Origin]);
                h_Nperevent_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]->Fill(n_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]);
            }
        }

        for (Int_t i_NMuons_rec = 0; i_NMuons_rec < NMuons_rec; i_NMuons_rec++)
        {

            Int_t PDG_Mu = PDGmum_rec[i_NMuons_rec];       // single rec mu PDG mum
            Double_t Pt_Mu = Pt_rec[i_NMuons_rec];         // single rec mu pT
            Double_t E_Mu = E_rec[i_NMuons_rec];           // single rec mu E
            Double_t Px_Mu = Px_rec[i_NMuons_rec];         // single rec mu px
            Double_t Py_Mu = Py_rec[i_NMuons_rec];         // single rec mu py
            Double_t Pz_Mu = Pz_rec[i_NMuons_rec];         // single rec mu pz
            Double_t Y_Mu = Y_rec[i_NMuons_rec];           // single rec mu y
            Double_t Eta_Mu = Eta_rec[i_NMuons_rec];       // single rec mu eta
            Double_t Phi_Mu = Phi_rec[i_NMuons_rec];       // single rec mu phi
            Double_t Theta_Mu = Theta_rec[i_NMuons_rec];   // single rec mu theta
            Double_t Charge_Mu = Charge_rec[i_NMuons_rec]; // single rec mu theta
            Double_t RAbs_Mu = RAtAbsEnd_rec[i_NMuons_rec];
            Int_t MatchTrig_Murec = MatchTrig_rec[i_NMuons_rec];
            Double_t pDCA_Mu = pDCA_rec[i_NMuons_rec];
            Int_t IsFromPowheg = fFrom_Powheg_rec[i_NMuons_rec];
            Bool_t DQ_Muon = kFALSE;
            Bool_t Selection_Muon[n_Muon_origin] = {kFALSE};

            Selection_Muon[0] = kTRUE; // Selection All Muons
            if ((TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500) || (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000))
            {
                Selection_Muon[1] = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600) || (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000))
            {
                Selection_Muon[2] = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu) > 100 && TMath::Abs(PDG_Mu) < 400) || (TMath::Abs(PDG_Mu) > 1000 && TMath::Abs(PDG_Mu) < 4000))
            {
                Selection_Muon[3] = kTRUE; // Selection Muons from LF
            }
            else if (TMath::Abs(PDG_Mu) == 23)
            {
                Selection_Muon[4] = kTRUE;
            }

            if (pDCA_Mu == 1)
                if (RAbs_Mu > 17.6 && RAbs_Mu < 89.5)
                    if (Eta_Mu > -4.0 && Eta_Mu < -2.5)
                        if (MatchTrig_Murec > 1)
                            DQ_Muon = kTRUE;

            if (!DQ_Muon)
                continue;
            // printf("Pt_Mu %f Y_Mu %f PDG_Mu %d",Pt_Mu,Y_Mu,PDG_Mu);
            h_PtPdg_Muon_Rec->Fill(Pt_Mu, PDG_Mu);
            h_YPdg_Muon_Rec->Fill(Y_Mu, PDG_Mu);
            if (Generator.Contains("Powheg") && IsFromPowheg == -1)
            {
                h_PtPdg_Muon_Rec_PowhegOnly->Fill(Pt_Mu, PDG_Mu);
                h_YPdg_Muon_Rec_PowhegOnly->Fill(Y_Mu, PDG_Mu);
            }

            for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
            {
                if (Selection_Muon[i_Muon_origin])
                {
                    h_PtY_Muon_Rec[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    n_Muon_Rec[i_Muon_origin]++;
                    if (Generator.Contains("Powheg") && IsFromPowheg == -1)
                    {
                        h_PtY_Muon_Rec_PowhegOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                        n_Muon_Rec_PowhegOnly[i_Muon_origin]++;
                    }
                }
            }
        }

        for (Int_t i_Muon_Origin = 0; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
        {
            h_Nperevent_Muon_Rec[i_Muon_Origin]->Fill(n_Muon_Rec[i_Muon_Origin]);
            if (Generator.Contains("Powheg"))
                h_Nperevent_Muon_Rec_PowhegOnly[i_Muon_Origin]->Fill(n_Muon_Rec_PowhegOnly[i_Muon_Origin]);
        }

        for (Int_t i_NDimu_gen = 0; i_NDimu_gen < NDimu_gen; i_NDimu_gen++)
        {
            
            Double_t Pt_DiMu = DimuPt_gen[i_NDimu_gen];      // gen dimuon pT
            Double_t Px_DiMu = DimuPx_gen[i_NDimu_gen];      // gen dimuon px
            Double_t Py_DiMu = DimuPy_gen[i_NDimu_gen];      // gen dimuon py
            Double_t Pz_DiMu = DimuPz_gen[i_NDimu_gen];      // gen dimuon pz
            Double_t Y_DiMu = DimuY_gen[i_NDimu_gen];        // gen dimuon y
            Double_t M_DiMu = DimuMass_gen[i_NDimu_gen];     // gen dimuon invariant mass
            Int_t Charge_DiMu = DimuCharge_gen[i_NDimu_gen]; // gen dimuon charge
            
            Double_t Pt_Mu0 = Pt_gen[DimuMu_gen[i_NDimu_gen][0]];
            Double_t Y_Mu0 = Y_gen[DimuMu_gen[i_NDimu_gen][0]];
            Double_t Eta_Mu0 = Eta_gen[DimuMu_gen[i_NDimu_gen][0]];

            Double_t Pt_Mu1 = Pt_gen[DimuMu_gen[i_NDimu_gen][1]];
            Double_t Y_Mu1 = Y_gen[DimuMu_gen[i_NDimu_gen][1]];
            Double_t Eta_Mu1 = Eta_gen[DimuMu_gen[i_NDimu_gen][1]];

            Int_t PDG_Mu0 = PDGmum_gen[DimuMu_gen[i_NDimu_gen][0]];
            Int_t PDG_Mu1 = PDGmum_gen[DimuMu_gen[i_NDimu_gen][1]];
            Int_t IsFromPowheg_Mu0 = 999;
            Int_t IsFromPowheg_Mu1 = 999;
            if (Generator.Contains("Powheg"))
            {
                IsFromPowheg_Mu0 = fFrom_Powheg_gen[DimuMu_gen[i_NDimu_gen][0]];
                IsFromPowheg_Mu1 = fFrom_Powheg_gen[DimuMu_gen[i_NDimu_gen][1]];
            }

            Bool_t HF_mu0 = kFALSE;
            Bool_t Charm_mu0 = kFALSE;
            Bool_t Beauty_mu0 = kFALSE;
            Bool_t LF_mu0 = kFALSE;

            Bool_t HF_mu1 = kFALSE;
            Bool_t Charm_mu1 = kFALSE;
            Bool_t Beauty_mu1 = kFALSE;
            Bool_t LF_mu1 = kFALSE;
            // cout<<"pdg code 1° mu: "<<PDG_Mu0<<"pdg code 2° mu: "<<PDG_Mu1<<endl;

            Double_t BR_times_frag_PYTHIA_Mu0 = 1.;
            Double_t BR_times_frag_MEAS_Mu0 = 1.;

            Double_t BR_times_frag_PYTHIA_Mu1 = 1.;
            Double_t BR_times_frag_MEAS_Mu1 = 1.;

            if ((TMath::Abs(PDG_Mu0) == 4) || (TMath::Abs(PDG_Mu0) > 400 && TMath::Abs(PDG_Mu0) < 500) || (TMath::Abs(PDG_Mu0) > 4000 && TMath::Abs(PDG_Mu0) < 5000))
            {
                HF_mu0 = kTRUE;
                Charm_mu0 = kTRUE;

                for (Int_t i_charm_hadron = 0; i_charm_hadron < n_charm_hadrons; i_charm_hadron++)
                {
                    if (TMath::Abs(PDG_Mu0) == PDG_charm_hadrons[i_charm_hadron])
                    {

                        BR_times_frag_PYTHIA_Mu0 = BR_charm_hadrons2mu_PYTHIA8_Monash[i_charm_hadron] * Frag_charm_hadrons_PYTHIA8_Monash[i_charm_hadron];
                        BR_times_frag_MEAS_Mu0 = BR_charm_hadrons2mu_MEAS[i_charm_hadron] * Frag_charm_hadrons_MEAS[i_charm_hadron];
                        // printf("PDG Mu0: %d --- BR_charm_hadrons2mu_PYTHIA8_Monash[i_charm_hadron] %0.2f --- Frag_charm_hadrons_PYTHIA8_Monash[i_charm_hadron] %0.2f\n",PDG_Mu0,BR_charm_hadrons2mu_PYTHIA8_Monash[i_charm_hadron],Frag_charm_hadrons_PYTHIA8_Monash[i_charm_hadron]);
                        // printf("PDG Mu0: %d --- BR_charm_hadrons2mu_MEAS[i_charm_hadron] %0.2f --- Frag_charm_hadrons_MEAS[i_charm_hadron] %0.2f\n",PDG_Mu0,BR_charm_hadrons2mu_MEAS[i_charm_hadron],Frag_charm_hadrons_MEAS[i_charm_hadron]);
                    }
                }
            }
            else if ((TMath::Abs(PDG_Mu0) == 5) || (TMath::Abs(PDG_Mu0) > 500 && TMath::Abs(PDG_Mu0) < 600) || (TMath::Abs(PDG_Mu0) > 5000 && TMath::Abs(PDG_Mu0) < 6000))
            {
                HF_mu0 = kTRUE;
                Beauty_mu0 = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu0) > 0 && TMath::Abs(PDG_Mu0) < 4) || (TMath::Abs(PDG_Mu0) > 100 && TMath::Abs(PDG_Mu0) < 400) || (TMath::Abs(PDG_Mu0) > 1000 && TMath::Abs(PDG_Mu0) < 4000))
            {
                LF_mu0 = kTRUE;
                Charm_mu0 = kFALSE;
                Beauty_mu0 = kFALSE;
                HF_mu0 = kFALSE;
            }

            if ((TMath::Abs(PDG_Mu1) == 4) || (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500) || (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000))
            {
                HF_mu1 = kTRUE;
                Charm_mu1 = kTRUE;
                for (Int_t i_charm_hadron = 0; i_charm_hadron < n_charm_hadrons; i_charm_hadron++)
                {
                    if (TMath::Abs(PDG_Mu1) == PDG_charm_hadrons[i_charm_hadron])
                    {
                        BR_times_frag_PYTHIA_Mu1 = BR_charm_hadrons2mu_PYTHIA8_Monash[i_charm_hadron] * Frag_charm_hadrons_PYTHIA8_Monash[i_charm_hadron];
                        BR_times_frag_MEAS_Mu1 = BR_charm_hadrons2mu_MEAS[i_charm_hadron] * Frag_charm_hadrons_MEAS[i_charm_hadron];
                    }
                }
            }
            else if ((TMath::Abs(PDG_Mu1) == 5) || (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600) || (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000))
            {
                HF_mu1 = kTRUE;
                Beauty_mu1 = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu1) > 0 && TMath::Abs(PDG_Mu1) < 4) || (TMath::Abs(PDG_Mu1) > 100 && TMath::Abs(PDG_Mu1) < 400) || (TMath::Abs(PDG_Mu1) > 1000 && TMath::Abs(PDG_Mu1) < 4000))
            {
                LF_mu1 = kTRUE;
                Charm_mu1 = kFALSE;
                Beauty_mu1 = kFALSE;
                HF_mu1 = kFALSE;
            }
            Bool_t DiMu_origin_Selection[n_DiMuon_origin] = {kFALSE};
            // 0 for Charm, 1 for Beauty, 2 for HF Mixed, 3 for LF, 4 for LF-HF Mixed, 5 for DY

            h_Pdg1Pdg2Pt_DiMuon_Gen->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Pt_DiMu);
            h_Pdg1Pdg2Y_DiMuon_Gen->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Y_DiMu);
            h_Pdg1Pdg2M_DiMuon_Gen->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), M_DiMu);

            if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
            {
                h_Pdg1Pdg2Pt_DiMuon_Gen_PowhegOnly->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Pt_DiMu);
                h_Pdg1Pdg2Y_DiMuon_Gen_PowhegOnly->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Y_DiMu);
                h_Pdg1Pdg2M_DiMuon_Gen_PowhegOnly->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), M_DiMu);
            }

            if (Charm_mu0 && Charm_mu1)
                DiMu_origin_Selection[0] = kTRUE;

            else if (Beauty_mu0 && Beauty_mu1)
                DiMu_origin_Selection[1] = kTRUE;

            else if ((Charm_mu0 && Beauty_mu1) || (Charm_mu1 && Beauty_mu0))
                DiMu_origin_Selection[2] = kTRUE;

            else if (LF_mu0 && LF_mu1)
                DiMu_origin_Selection[3] = kTRUE;

            else if ((LF_mu0 && HF_mu1) || (LF_mu1 && HF_mu0))
                DiMu_origin_Selection[4] = kTRUE;

            else if (PDG_Mu0 == 23 && PDG_Mu1 == 23)
            {
                DiMu_origin_Selection[5] = kTRUE;
                // if ((Sim_for_Z) && Pt_Mu0 > 0.9 && Pt_Mu1 > 0.9)
                // {
                //     h_PtM_DiMuon_Gen_Z_ptmucut09->Fill(Pt_DiMu, M_DiMu);
                //     h_PtY_DiMuon_Gen_Z_ptmucut09->Fill(Pt_DiMu, Y_DiMu);
                // }

                // if ((Sim_for_Z) && Pt_Mu0 > 10 && Pt_Mu1 > 10)
                // {
                //     h_PtM_DiMuon_Gen_Z_ptmucut10->Fill(Pt_DiMu, M_DiMu);
                //     h_PtY_DiMuon_Gen_Z_ptmucut10->Fill(Pt_DiMu, Y_DiMu);
                // }

                // if ((Sim_for_Z) && Pt_Mu0 > 20 && Pt_Mu1 > 20)
                // {
                //     h_PtM_DiMuon_Gen_Z_ptmucut20->Fill(Pt_DiMu, M_DiMu);
                //     h_PtY_DiMuon_Gen_Z_ptmucut20->Fill(Pt_DiMu, Y_DiMu);
                // }

                h_YGamma_YDimuon->Fill(Y_DiMu, Y_Gamma);
            }

            Bool_t DQ_Dimuon = kFALSE;
            if ((Y_DiMu > -4.0 && Y_DiMu < -2.5) && (Eta_Mu0 > -4.0 && Eta_Mu0 < -2.5) && (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5) && (Charge_DiMu == 0))
                DQ_Dimuon = kTRUE;

            Double_t charm_correction = (BR_times_frag_MEAS_Mu0 * BR_times_frag_MEAS_Mu1) / (BR_times_frag_PYTHIA_Mu0 * BR_times_frag_PYTHIA_Mu1);
            if (DQ_Dimuon)
            {
                h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Pt_DiMu);
                h_Pdg1Pdg2Y_DiMuon_Gen_DQcut->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Y_DiMu);
                h_Pdg1Pdg2M_DiMuon_Gen_DQcut->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), M_DiMu);
                // if ((TMath::Abs(PDG_Mu0) == 411 || TMath::Abs(PDG_Mu0) == 421) && TMath::Abs(PDG_Mu1) == 4122)
                //     printf("charm correction %0.4f \n", charm_correction);

                h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_Charm_corrected->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Pt_DiMu, charm_correction);
                h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_Charm_corrected->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Y_DiMu, charm_correction);
                h_Pdg1Pdg2M_DiMuon_Gen_DQcut_Charm_corrected->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), M_DiMu, charm_correction);

                if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                {
                    h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_PowhegOnly->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Pt_DiMu);
                    h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_PowhegOnly->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Y_DiMu);
                    h_Pdg1Pdg2M_DiMuon_Gen_DQcut_PowhegOnly->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), M_DiMu);
                }

                // if ((DiMu_origin_Selection[5]) && (Sim_for_Z) && Pt_Mu0 > 0.9 && Pt_Mu1 > 0.9)
                // {
                //     h_PtM_DiMuon_Gen_Z_DQcut_ptmucut09->Fill(Pt_DiMu, M_DiMu);
                //     h_PtY_DiMuon_Gen_Z_DQcut_ptmucut09->Fill(Pt_DiMu, Y_DiMu);
                // }

                // if ((DiMu_origin_Selection[5]) && (Sim_for_Z) && Pt_Mu0 > 10 && Pt_Mu1 > 10)
                // {
                //     h_PtM_DiMuon_Gen_Z_DQcut_ptmucut10->Fill(Pt_DiMu, M_DiMu);
                //     h_PtY_DiMuon_Gen_Z_DQcut_ptmucut10->Fill(Pt_DiMu, Y_DiMu);
                // }

                // if ((DiMu_origin_Selection[5]) && (Sim_for_Z) && Pt_Mu0 > 20 && Pt_Mu1 > 20)
                // {
                //     h_PtM_DiMuon_Gen_Z_DQcut_ptmucut20->Fill(Pt_DiMu, M_DiMu);
                //     h_PtY_DiMuon_Gen_Z_DQcut_ptmucut20->Fill(Pt_DiMu, Y_DiMu);
                // }
            }

            for (Int_t i_DiMuon_origin = 0; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
            {
                if (DiMu_origin_Selection[i_DiMuon_origin])
                {
                    // 0 for Charm, 1 for Beauty, 2 for HF Mixed, 3 for LF, 4 for LF-HF Mixed, 5 for DY
                    h_PtM_DiMuon_Gen[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                    h_PtY_DiMuon_Gen[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                    n_DiMuon_Gen[i_DiMuon_origin]++;

                    if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                    {
                        h_PtM_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                        h_PtY_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                        n_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]++;
                    }

                    if (DQ_Dimuon)
                    {

                        h_PtM_DiMuon_Gen_DQcut[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                        h_PtY_DiMuon_Gen_DQcut[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                        n_DiMuon_Gen_DQcut[i_DiMuon_origin]++;

                        if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                        {
                            h_PtM_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                            h_PtY_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                            n_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]++;
                        }

                        if (i_DiMuon_origin == 0)
                        {
                            h_PtM_DiMuon_Gen_DQcut_Charm_corrected->Fill(Pt_DiMu, M_DiMu, charm_correction);
                            h_PtY_DiMuon_Gen_DQcut_Charm_corrected->Fill(Pt_DiMu, Y_DiMu, charm_correction);
                        }
                    }
                }
            }
        }

        for (Int_t i_Dimuon_Origin = 0; i_Dimuon_Origin < n_DiMuon_origin; i_Dimuon_Origin++)
        {

            h_Nperevent_DiMuon_Gen[i_Dimuon_Origin]->Fill(n_DiMuon_Gen[i_Dimuon_Origin]);
            h_Nperevent_DiMuon_Gen_DQcut[i_Dimuon_Origin]->Fill(n_DiMuon_Gen_DQcut[i_Dimuon_Origin]);
            if (Generator.Contains("Powheg"))
            {
                h_Nperevent_DiMuon_Gen_PowhegOnly[i_Dimuon_Origin]->Fill(n_DiMuon_Gen_PowhegOnly[i_Dimuon_Origin]);
                h_Nperevent_DiMuon_Gen_DQcut_PowhegOnly[i_Dimuon_Origin]->Fill(n_DiMuon_Gen_DQcut_PowhegOnly[i_Dimuon_Origin]);
            }
        }

        for (Int_t i_NDimu_rec = 0; i_NDimu_rec < NDimu_rec; i_NDimu_rec++)
        {

            Double_t Pt_DiMu = DimuPt_rec[i_NDimu_rec];      // rec dimuon pT
            Double_t Px_DiMu = DimuPx_rec[i_NDimu_rec];      // rec dimuon px
            Double_t Py_DiMu = DimuPy_rec[i_NDimu_rec];      // rec dimuon py
            Double_t Pz_DiMu = DimuPz_rec[i_NDimu_rec];      // rec dimuon pz
            Double_t Y_DiMu = DimuY_rec[i_NDimu_rec];        // rec dimuon y
            Double_t M_DiMu = DimuMass_rec[i_NDimu_rec];     // rec dimuon invariant mass
            Int_t Charge_DiMu = DimuCharge_rec[i_NDimu_rec]; // rec dimuon charge

            Double_t Pt_Mu0 = Pt_rec[DimuMu_rec[i_NDimu_rec][0]];
            Double_t Y_Mu0 = Y_rec[DimuMu_rec[i_NDimu_rec][0]];
            Int_t PDG_Mu0 = PDGmum_rec[DimuMu_rec[i_NDimu_rec][0]];
            Double_t Charge_Mu0 = Charge_rec[DimuMu_rec[i_NDimu_rec][0]];
            Double_t RAbs_Mu0 = RAtAbsEnd_rec[DimuMu_rec[i_NDimu_rec][0]];
            Double_t pDCA_Mu0 = pDCA_rec[DimuMu_rec[i_NDimu_rec][0]];
            Double_t Eta_Mu0 = Eta_rec[DimuMu_rec[i_NDimu_rec][0]];
            Double_t Phi_Mu0 = Phi_rec[DimuMu_rec[i_NDimu_rec][0]];
            Int_t IsFromPowheg_Mu0 = fFrom_Powheg_rec[DimuMu_rec[i_NDimu_rec][0]];

            Double_t Y_Mu1 = Y_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Pt_Mu1 = Pt_rec[DimuMu_rec[i_NDimu_rec][1]];
            Int_t PDG_Mu1 = PDGmum_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Charge_Mu1 = Charge_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t RAbs_Mu1 = RAtAbsEnd_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t pDCA_Mu1 = pDCA_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Eta_Mu1 = Eta_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Phi_Mu1 = Phi_rec[DimuMu_rec[i_NDimu_rec][1]];
            Int_t IsFromPowheg_Mu1 = fFrom_Powheg_rec[DimuMu_rec[i_NDimu_rec][1]];

            /*
            printf("PDG_Mu0 %0.2f\n", Eta_Mu0);
            printf("DimuMu_rec[i_NDimu_rec][0] %i\n", DimuMu_rec[i_NDimu_rec][0]);

            if ((Charge_DiMu == 0) && (Y_DiMu > -4.0 && Y_DiMu < -2.5))
            {

                printf("PDG_Mu0 %i (DimuMu_rec[i_NDimu_rec][0] %i) || PDG_Mu1 %i || Y_DiMu %0.2f || Eta_Mu0 %0.2f\n", PDG_Mu0, DimuMu_rec[i_NDimu_rec][0], PDG_Mu1, Y_DiMu, Eta_Mu0);
                counter_test++;
            }
            */
            Bool_t DQ_Dimuon = kFALSE;
            if (pDCA_Mu0 == 1 && pDCA_Mu1 == 1 && (Charge_DiMu == 0))

                if ((RAbs_Mu0 > 17.6 && RAbs_Mu0 < 89.5) && (RAbs_Mu1 > 17.6 && RAbs_Mu1 < 89.5))

                    if ((Y_DiMu > -4.0 && Y_DiMu < -2.5) && (Eta_Mu0 > -4.0 && Eta_Mu0 < -2.5) && (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5))

                        if (DimuMatch_rec[i_NDimu_rec] == 2)

                            DQ_Dimuon = kTRUE;

            if (!DQ_Dimuon)
                continue;

            Bool_t HF_mu0 = kFALSE;
            Bool_t Charm_mu0 = kFALSE;
            Bool_t Beauty_mu0 = kFALSE;
            Bool_t LF_mu0 = kFALSE;

            Bool_t HF_mu1 = kFALSE;
            Bool_t Charm_mu1 = kFALSE;
            Bool_t Beauty_mu1 = kFALSE;
            Bool_t LF_mu1 = kFALSE;

            if ((TMath::Abs(PDG_Mu0) == 4) || (TMath::Abs(PDG_Mu0) > 400 && TMath::Abs(PDG_Mu0) < 500) || (TMath::Abs(PDG_Mu0) > 4000 && TMath::Abs(PDG_Mu0) < 5000))
            {
                HF_mu0 = kTRUE;
                Charm_mu0 = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu0) == 5) || (TMath::Abs(PDG_Mu0) > 500 && TMath::Abs(PDG_Mu0) < 600) || (TMath::Abs(PDG_Mu0) > 5000 && TMath::Abs(PDG_Mu0) < 6000))
            {
                HF_mu0 = kTRUE;
                Beauty_mu0 = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu0) > 0 && TMath::Abs(PDG_Mu0) < 4) || (TMath::Abs(PDG_Mu0) > 100 && TMath::Abs(PDG_Mu0) < 400) || (TMath::Abs(PDG_Mu0) > 1000 && TMath::Abs(PDG_Mu0) < 4000))
            {
                LF_mu0 = kTRUE;
                Charm_mu0 = kFALSE;
                Beauty_mu0 = kFALSE;
                HF_mu0 = kFALSE;
            }

            if ((TMath::Abs(PDG_Mu1) == 4) || (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500) || (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000))
            {
                HF_mu1 = kTRUE;
                Charm_mu1 = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu1) == 5) || (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600) || (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000))
            {
                HF_mu1 = kTRUE;
                Beauty_mu1 = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu1) > 0 && TMath::Abs(PDG_Mu1) < 4) || (TMath::Abs(PDG_Mu1) > 100 && TMath::Abs(PDG_Mu1) < 400) || (TMath::Abs(PDG_Mu1) > 1000 && TMath::Abs(PDG_Mu1) < 4000))
            {
                LF_mu1 = kTRUE;
                Charm_mu1 = kFALSE;
                Beauty_mu1 = kFALSE;
                HF_mu1 = kFALSE;
            }
            Bool_t DiMu_origin_Selection[n_DiMuon_origin] = {kFALSE};

            // 0 for Charm, 1 for Beauty, 2 for HF Mixed, 3 for LF, 4 for LF-HF Mixed, 5 for DY

            h_Pdg1Pdg2Pt_DiMuon_Rec->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Pt_DiMu);
            h_Pdg1Pdg2Y_DiMuon_Rec->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Y_DiMu);
            h_Pdg1Pdg2M_DiMuon_Rec->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), M_DiMu);

            if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
            {
                h_Pdg1Pdg2Pt_DiMuon_Rec_PowhegOnly->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Pt_DiMu);
                h_Pdg1Pdg2Y_DiMuon_Rec_PowhegOnly->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Y_DiMu);
                h_Pdg1Pdg2M_DiMuon_Rec_PowhegOnly->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), M_DiMu);
            }

            if (Charm_mu0 && Charm_mu1)
                DiMu_origin_Selection[0] = kTRUE;

            else if (Beauty_mu0 && Beauty_mu1)
                DiMu_origin_Selection[1] = kTRUE;

            else if ((Charm_mu0 && Beauty_mu1) || (Charm_mu1 && Beauty_mu0))
                DiMu_origin_Selection[2] = kTRUE;

            else if (LF_mu0 && LF_mu1)
                DiMu_origin_Selection[3] = kTRUE;

            else if ((LF_mu0 && HF_mu1) || (LF_mu1 && HF_mu0))
                DiMu_origin_Selection[4] = kTRUE;

            else if (PDG_Mu0 == 23 && PDG_Mu1 == 23)
            {

                DiMu_origin_Selection[5] = kTRUE;

                // if ((Sim_for_Z) && Pt_Mu0 > 0.9 && Pt_Mu1 > 0.9)
                // {
                //     h_PtM_DiMuon_Rec_Z_ptmucut09->Fill(Pt_DiMu, M_DiMu);
                //     h_PtY_DiMuon_Rec_Z_ptmucut09->Fill(Pt_DiMu, Y_DiMu);
                // }

                // if ((Sim_for_Z) && Pt_Mu0 > 10 && Pt_Mu1 > 10)
                // {
                //     h_PtM_DiMuon_Rec_Z_ptmucut10->Fill(Pt_DiMu, M_DiMu);
                //     h_PtY_DiMuon_Rec_Z_ptmucut10->Fill(Pt_DiMu, Y_DiMu);
                // }

                // if ((Sim_for_Z) && Pt_Mu0 > 20 && Pt_Mu1 > 20)
                // {
                //     h_PtM_DiMuon_Rec_Z_ptmucut20->Fill(Pt_DiMu, M_DiMu);
                //     h_PtY_DiMuon_Rec_Z_ptmucut20->Fill(Pt_DiMu, Y_DiMu);
                // }
            }

            for (Int_t i_DiMuon_origin = 0; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
            {
                *Pt_Dimu_Rec = 999;
                *M_Dimu_Rec = 999;
                if (DiMu_origin_Selection[i_DiMuon_origin])
                {
                    h_PtM_DiMuon_Rec[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                    h_PtY_DiMuon_Rec[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                    n_DiMuon_Rec[i_DiMuon_origin]++;

                    if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                    {
                        h_PtM_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                        h_PtY_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                        n_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]++;
                    }
                    if (M_DiMu > 4 && M_DiMu < 40)
                    {
                        *Pt_Dimu_Rec = Pt_DiMu;
                        *M_Dimu_Rec = M_DiMu;
                        Tree_DiMuon_Rec[i_DiMuon_origin]->Fill();

                        if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                        {
                            *Pt_Dimu_Rec_PowhegOnly = Pt_DiMu;
                            *M_Dimu_Rec_PowhegOnly = M_DiMu;
                            Tree_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Fill();
                        }
                    }
                }
            }
        }

        for (Int_t i_Dimuon_Origin = 0; i_Dimuon_Origin < n_DiMuon_origin; i_Dimuon_Origin++)
        {
            h_Nperevent_DiMuon_Rec[i_Dimuon_Origin]->Fill(n_DiMuon_Rec[i_Dimuon_Origin]);
            if (Generator.Contains("Powheg"))
                h_Nperevent_DiMuon_Rec_PowhegOnly[i_Dimuon_Origin]->Fill(n_DiMuon_Rec_PowhegOnly[i_Dimuon_Origin]);
        }
    }
    printf("counter_test %i\n", counter_test);

    TFile fOut_Tree(Form("%s/%s", dir_fileOut.Data(), file_out_tree.Data()), "RECREATE");
    fOut_Tree.cd();
    h_Nevents->Write();
    for (Int_t i_Dimuon_origin = 0; i_Dimuon_origin < n_DiMuon_origin; i_Dimuon_origin++)
    {
        if (Tree_DiMuon_Rec[i_Dimuon_origin]->GetEntries() > 0)
            Tree_DiMuon_Rec[i_Dimuon_origin]->Write(0, 2, 0);
        if (Generator.Contains("Powheg"))
        {
            if (Tree_DiMuon_Rec_PowhegOnly[i_Dimuon_origin]->GetEntries() > 0)
                Tree_DiMuon_Rec_PowhegOnly[i_Dimuon_origin]->Write(0, 2, 0);
        }
    }
    // Tree_DiMuon_Rec->Write(0, 2, 0);
    fOut_Tree.Close();

    TFile fOut(Form("%s/%s", dir_fileOut.Data(), file_out.Data()), "RECREATE");
    fOut.cd();
    h_Nevents->Write();
    if (!fOut.GetDirectory("Muon_Gen"))
        fOut.mkdir("Muon_Gen");

    if (!fOut.GetDirectory("Muon_Rec"))
        fOut.mkdir("Muon_Rec");

    if (!fOut.GetDirectory("DiMuon_Gen"))
        fOut.mkdir("DiMuon_Gen");

    if (!fOut.GetDirectory("DiMuon_Rec"))
        fOut.mkdir("DiMuon_Rec");

    if (Generator.Contains("Powheg"))
    {
        if (!fOut.GetDirectory("Muon_Gen/PowhegOnly"))
            fOut.mkdir("Muon_Gen/PowhegOnly");
        if (!fOut.GetDirectory("Muon_Rec/PowhegOnly"))
            fOut.mkdir("Muon_Rec/PowhegOnly");

        if (!fOut.GetDirectory("DiMuon_Gen/PowhegOnly"))
            fOut.mkdir("DiMuon_Gen/PowhegOnly");
        if (!fOut.GetDirectory("DiMuon_Rec/PowhegOnly"))
            fOut.mkdir("DiMuon_Rec/PowhegOnly");
    }

    fOut.cd("Muon_Gen");
    if (Generator.Contains("Powheg"))
    {
        h_PtQ_Muon_Gen_LF_Powheg->Write();
        h_PtQ_Muon_Gen_LF_Pythia->Write();
        h_YQ_Muon_Gen_LF_Powheg->Write();
        h_YQ_Muon_Gen_LF_Pythia->Write();
    }

    h_PtPdg_Muon_Gen->Write(0, 2, 0);
    h_YPdg_Muon_Gen->Write(0, 2, 0);

    h_PtPdg_Muon_Gen_DQcut->Write(0, 2, 0);
    h_YPdg_Muon_Gen_DQcut->Write(0, 2, 0);

    if (Generator.Contains("Powheg"))
    {
        fOut.cd("Muon_Gen/PowhegOnly");
        h_PtPdg_Muon_Gen_PowhegOnly->Write(0, 2, 0);
        h_YPdg_Muon_Gen_PowhegOnly->Write(0, 2, 0);
        h_PtPdg_Muon_Gen_DQcut_PowhegOnly->Write(0, 2, 0);
        h_YPdg_Muon_Gen_DQcut_PowhegOnly->Write(0, 2, 0);
    }

    fOut.cd();
    fOut.cd("Muon_Rec");
    h_PtPdg_Muon_Rec->Write(0, 2, 0);
    h_YPdg_Muon_Rec->Write(0, 2, 0);
    if (Generator.Contains("Powheg"))
    {
        fOut.cd("Muon_Rec/PowhegOnly");
        h_PtPdg_Muon_Rec_PowhegOnly->Write(0, 2, 0);
        h_YPdg_Muon_Rec_PowhegOnly->Write(0, 2, 0);
    }

    for (Int_t i_Muon_Origin = 0; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
    {
        // Gen Muon saving
        fOut.cd();
        fOut.cd("Muon_Gen");
        if (h_Nperevent_Muon_Gen[i_Muon_Origin]->GetMean() > 0)
            h_Nperevent_Muon_Gen[i_Muon_Origin]->Write(0, 2, 0);

        if (h_PtY_Muon_Gen[i_Muon_Origin]->GetMean() > 0)
            h_PtY_Muon_Gen[i_Muon_Origin]->Write(0, 2, 0);

        if (h_Nperevent_Muon_Gen_DQcut[i_Muon_Origin]->GetMean() > 0)
            h_Nperevent_Muon_Gen_DQcut[i_Muon_Origin]->Write(0, 2, 0);

        if (h_PtY_Muon_Gen_DQcut[i_Muon_Origin]->GetMean() > 0)
            h_PtY_Muon_Gen_DQcut[i_Muon_Origin]->Write(0, 2, 0);

        // Gen Muon from Powheg saving

        if (Generator.Contains("Powheg"))
        {
            fOut.cd("Muon_Gen/PowhegOnly");
            if (h_Nperevent_Muon_Gen_PowhegOnly[i_Muon_Origin]->GetMean() > 0)
                h_Nperevent_Muon_Gen_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);

            if (h_PtY_Muon_Gen_PowhegOnly[i_Muon_Origin]->GetMean() > 0)
                h_PtY_Muon_Gen_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);

            if (h_Nperevent_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]->GetMean() > 0)
                h_Nperevent_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);

            if (h_PtY_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]->GetMean() > 0)
                h_PtY_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
        }
        // Rec Muon saving
        fOut.cd();
        fOut.cd("Muon_Rec");
        if (h_PtY_Muon_Rec[i_Muon_Origin]->GetMean() > 0)
            h_PtY_Muon_Rec[i_Muon_Origin]->Write(0, 2, 0);

        if (h_Nperevent_Muon_Rec[i_Muon_Origin]->GetMean() > 0)
            h_Nperevent_Muon_Rec[i_Muon_Origin]->Write(0, 2, 0);
        // Rec Muon from Powheg saving
        if (Generator.Contains("Powheg"))
        {
            fOut.cd("Muon_Rec/PowhegOnly");
            if (h_PtY_Muon_Rec_PowhegOnly[i_Muon_Origin]->GetMean() > 0)
                h_PtY_Muon_Rec_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);

            if (h_Nperevent_Muon_Rec_PowhegOnly[i_Muon_Origin]->GetMean() > 0)
                h_Nperevent_Muon_Rec_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
        }
    }

    fOut.cd();
    fOut.cd("DiMuon_Gen");
    h_Pdg1Pdg2Pt_DiMuon_Gen->Write(0, 2, 0);
    h_Pdg1Pdg2Y_DiMuon_Gen->Write(0, 2, 0);
    h_Pdg1Pdg2M_DiMuon_Gen->Write(0, 2, 0);

    h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut->Write(0, 2, 0);
    h_Pdg1Pdg2Y_DiMuon_Gen_DQcut->Write(0, 2, 0);
    h_Pdg1Pdg2M_DiMuon_Gen_DQcut->Write(0, 2, 0);

    fOut.cd();
    fOut.cd("DiMuon_Rec");
    h_Pdg1Pdg2Pt_DiMuon_Rec->Write(0, 2, 0);
    h_Pdg1Pdg2Y_DiMuon_Rec->Write(0, 2, 0);
    h_Pdg1Pdg2M_DiMuon_Rec->Write(0, 2, 0);

    fOut.cd();
    for (Int_t i_DiMuon_origin = 0; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
    {
        fOut.cd("DiMuon_Gen");
        if (h_Nperevent_DiMuon_Gen[i_DiMuon_origin]->GetMean() > 0)
            h_Nperevent_DiMuon_Gen[i_DiMuon_origin]->Write(0, 2, 0);
        if (h_PtM_DiMuon_Gen[i_DiMuon_origin]->GetMean() > 0)
            h_PtM_DiMuon_Gen[i_DiMuon_origin]->Write(0, 2, 0);
        if (h_PtY_DiMuon_Gen[i_DiMuon_origin]->GetMean() > 0)
            h_PtY_DiMuon_Gen[i_DiMuon_origin]->Write(0, 2, 0);

        if (h_Nperevent_DiMuon_Gen_DQcut[i_DiMuon_origin]->GetMean() > 0)
            h_Nperevent_DiMuon_Gen_DQcut[i_DiMuon_origin]->Write(0, 2, 0);
        if (h_PtM_DiMuon_Gen_DQcut[i_DiMuon_origin]->GetMean() > 0)
            h_PtM_DiMuon_Gen_DQcut[i_DiMuon_origin]->Write(0, 2, 0);
        if (h_PtY_DiMuon_Gen_DQcut[i_DiMuon_origin]->GetMean() > 0)
            h_PtY_DiMuon_Gen_DQcut[i_DiMuon_origin]->Write(0, 2, 0);

        if (Generator.Contains("Powheg"))
        {
            fOut.cd("DiMuon_Gen/PowhegOnly");
            if (h_Nperevent_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->GetMean() > 0)
                h_Nperevent_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
            if (h_PtM_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->GetMean() > 0)
                h_PtM_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
            if (h_PtY_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->GetMean() > 0)
                h_PtY_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);

            if (h_Nperevent_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->GetMean() > 0)
                h_Nperevent_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
            if (h_PtM_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->GetMean() > 0)
                h_PtM_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
            if (h_PtY_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->GetMean() > 0)
                h_PtY_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
        }

        fOut.cd();
        fOut.cd("DiMuon_Rec");
        if (h_PtM_DiMuon_Rec[i_DiMuon_origin]->GetMean() > 0)
            h_PtM_DiMuon_Rec[i_DiMuon_origin]->Write(0, 2, 0);
        if (h_PtY_DiMuon_Rec[i_DiMuon_origin]->GetMean() > 0)
            h_PtY_DiMuon_Rec[i_DiMuon_origin]->Write(0, 2, 0);
        if (h_Nperevent_DiMuon_Rec[i_DiMuon_origin]->GetMean() > 0)
            h_Nperevent_DiMuon_Rec[i_DiMuon_origin]->Write(0, 2, 0);

        if (Generator.Contains("Powheg"))
        {
            fOut.cd("DiMuon_Rec/PowhegOnly");
            if (h_Nperevent_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->GetMean() > 0)
                h_Nperevent_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
            if (h_PtM_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->GetMean() > 0)
                h_PtM_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
            if (h_PtY_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->GetMean() > 0)
                h_PtY_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
        }
    }
    // if (Sim_for_Z)
    // {
    //     fOut.cd();
    //     h_PtM_Gamma->Write(0, 2, 0);
    //     h_PtY_Gamma->Write(0, 2, 0);
    //     h_YGamma_YDimuon->Write(0, 2, 0);
    //     if (!fOut.GetDirectory("DiMuon_Gen_Z"))
    //         fOut.mkdir("DiMuon_Gen_Z");
    //     fOut.cd("DiMuon_Gen_Z");
    //     h_PtM_DiMuon_Gen_Z_ptmucut09->Write(0, 2, 0);
    //     h_PtY_DiMuon_Gen_Z_ptmucut09->Write(0, 2, 0);
    //     h_PtM_DiMuon_Gen_Z_ptmucut10->Write(0, 2, 0);
    //     h_PtY_DiMuon_Gen_Z_ptmucut10->Write(0, 2, 0);
    //     h_PtM_DiMuon_Gen_Z_ptmucut20->Write(0, 2, 0);
    //     h_PtY_DiMuon_Gen_Z_ptmucut20->Write(0, 2, 0);

    //     h_PtM_DiMuon_Gen_Z_DQcut_ptmucut09->Write(0, 2, 0);
    //     h_PtY_DiMuon_Gen_Z_DQcut_ptmucut09->Write(0, 2, 0);
    //     h_PtM_DiMuon_Gen_Z_DQcut_ptmucut10->Write(0, 2, 0);
    //     h_PtY_DiMuon_Gen_Z_DQcut_ptmucut10->Write(0, 2, 0);
    //     h_PtM_DiMuon_Gen_Z_DQcut_ptmucut20->Write(0, 2, 0);
    //     h_PtY_DiMuon_Gen_Z_DQcut_ptmucut20->Write(0, 2, 0);

    //     fOut.cd();
    //     if (!fOut.GetDirectory("DiMuon_Rec_Z"))
    //         fOut.mkdir("DiMuon_Rec_Z");
    //     fOut.cd("DiMuon_Rec_Z");
    //     h_PtM_DiMuon_Rec_Z_ptmucut09->Write(0, 2, 0);
    //     h_PtY_DiMuon_Rec_Z_ptmucut09->Write(0, 2, 0);
    //     h_PtM_DiMuon_Rec_Z_ptmucut10->Write(0, 2, 0);
    //     h_PtY_DiMuon_Rec_Z_ptmucut10->Write(0, 2, 0);
    //     h_PtM_DiMuon_Rec_Z_ptmucut20->Write(0, 2, 0);
    //     h_PtY_DiMuon_Rec_Z_ptmucut20->Write(0, 2, 0);
    // }
    if (Generator.Contains("HF"))
    {
        fOut.cd();
        if (!fOut.GetDirectory("DiMuon_Charm_corrected"))
            fOut.mkdir("DiMuon_Charm_corrected");
        fOut.cd("DiMuon_Charm_corrected");

        h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);
        h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);
        h_Pdg1Pdg2M_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);
        h_PtM_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);
        h_PtY_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);

        fOut.cd();
        if (!fOut.GetDirectory("HF_Hadron"))
            fOut.mkdir("HF_Hadron");
        fOut.cd("HF_Hadron");
        h_PdgPtY_HFHadron_prompt->Write(0, 2, 0);
        h_PdgPtY_HFHadron_notprompt->Write(0, 2, 0);
        if (Generator.Contains("Powheg"))
        {
            h_PdgPtY_HFHadron_prompt_PowhegOnly->Write(0, 2, 0);
            h_PdgPtY_HFHadron_notprompt_PowhegOnly->Write(0, 2, 0);
        }
        h_PdgPt_HFHadron_prompt->Write(0, 2, 0);
        h_PdgY_HFHadron_prompt->Write(0, 2, 0);
        h_PdgPt_HFHadron_notprompt->Write(0, 2, 0);
        h_PdgY_HFHadron_notprompt->Write(0, 2, 0);
    }

    fOut.Close();
    return;
}
