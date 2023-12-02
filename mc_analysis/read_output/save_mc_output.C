#include "save_mc_output.h"
// sistema bene powheg muoni dentro not geant nel caso rec
// fai lo stesso per dimuoni
// salvo low mass muons

void save_mc_output(
    // TString RunMode = "test_new_prompt_LHC22b3",
    // TString RunMode = "pythia8_purifykineoff_test",
    // TString RunMode = "SoftQCD_inel_LFoff_Def",
    TString RunMode = "Merged_LHC22c1",
    // TString RunMode = "LHC23i1",
    // TString RunMode = "SoftQCD_inel_Def",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC23i1_Version5_AliAOD_HF_LF",
    TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC22c1/294009/output",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/pythia8_purifykineoff_test",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LF_test",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/test_beauty_sim_2",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/test_charm_sim",
    // TString dir_fileIn = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC18p_DY_100k_Version2_AOD",
    // TString dir_fileOut = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC23i1_Version5_AliAOD_HF_LF",
    TString dir_fileOut = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC22c1/294009/output",
    // TString dir_fileOut = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/pythia8_purifykineoff_test",
    // TString dir_fileOut = "/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim",
    // TString dir_fileOut = "/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/sim",
    // Int_t RunNumber = 100000,
    Int_t RunNumber = 294009,
    TString Generator = "Geant_HF",
    TString prefix_filename = "MCDimuHFTree"
    // TString prefix_filename = "pythia_sim_74_DefaultBR"
    )
{

    Set_Histograms(Generator);

    TString filename;
    filename.Form("%s_%s_%d.root", RunMode.Data(), prefix_filename.Data(), RunNumber);

    TString file_out;
    file_out.Form("%s_MC_output_Hist_%d.root", RunMode.Data(), RunNumber);
    if (Generator.Contains("Pythia"))
        file_out.Form("%s_%s_output_Hist_%d.root", RunMode.Data(), prefix_filename.Data(), RunNumber);
    // file_out.Form("old_%s_MC_output_Hist_%d.root", RunMode.Data(), RunNumber);

    TString file_out_tree;
    file_out_tree.Form("%s_MC_output_Tree_%d.root", RunMode.Data(), RunNumber);

    printf("%s/%s\n", dir_fileOut.Data(), file_out.Data());
    printf("%s/%s\n", dir_fileOut.Data(), file_out_tree.Data());

    Double_t *Pt_Dimu_Rec = new Double_t;
    Double_t *Pt_Dimu_Rec_FullMass = new Double_t;
    Double_t *Pt_Dimu_Rec_LowMass = new Double_t;
    Double_t *Pt_Dimu_Rec_LowMass_LowPt = new Double_t;
    Double_t *Pt_Dimu_Rec_HighMass = new Double_t;

    Double_t *M_Dimu_Rec = new Double_t;
    Double_t *M_Dimu_Rec_FullMass = new Double_t;
    Double_t *M_Dimu_Rec_LowMass = new Double_t;
    Double_t *M_Dimu_Rec_LowMass_LowPt = new Double_t;
    Double_t *M_Dimu_Rec_HighMass = new Double_t;

    Double_t *Pt_Dimu_Rec_PowhegOnly = new Double_t;
    Double_t *Pt_Dimu_Rec_FullMass_PowhegOnly = new Double_t;
    Double_t *Pt_Dimu_Rec_LowMass_PowhegOnly = new Double_t;
    Double_t *Pt_Dimu_Rec_LowMass_LowPt_PowhegOnly = new Double_t;
    Double_t *Pt_Dimu_Rec_HighMass_PowhegOnly = new Double_t;

    Double_t *M_Dimu_Rec_PowhegOnly = new Double_t;
    Double_t *M_Dimu_Rec_FullMass_PowhegOnly = new Double_t;
    Double_t *M_Dimu_Rec_LowMass_PowhegOnly = new Double_t;
    Double_t *M_Dimu_Rec_LowMass_LowPt_PowhegOnly = new Double_t;
    Double_t *M_Dimu_Rec_HighMass_PowhegOnly = new Double_t;

    TTree *Tree_DiMuon_Rec[n_DiMuon_origin];
    TTree *Tree_DiMuon_Rec_PowhegOnly[n_DiMuon_origin];

    TTree *Tree_DiMuon_Rec_FullMass[n_DiMuon_origin];
    TTree *Tree_DiMuon_Rec_FullMass_PowhegOnly[n_DiMuon_origin];

    TTree *Tree_DiMuon_Rec_LowMass[n_DiMuon_origin];
    TTree *Tree_DiMuon_Rec_LowMass_PowhegOnly[n_DiMuon_origin];

    TTree *Tree_DiMuon_Rec_LowMass_LowPt[n_DiMuon_origin];
    TTree *Tree_DiMuon_Rec_LowMass_LowPt_PowhegOnly[n_DiMuon_origin];

    TTree *Tree_DiMuon_Rec_HighMass[n_DiMuon_origin];
    TTree *Tree_DiMuon_Rec_HighMass_PowhegOnly[n_DiMuon_origin];

    for (Int_t i_DiMuon_origin = 0; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
    {
        Tree_DiMuon_Rec[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with 4 < m <30 Gev");
        Tree_DiMuon_Rec[i_DiMuon_origin]->Branch("m", M_Dimu_Rec, "m/D");
        Tree_DiMuon_Rec[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec, "pt/D");

        Tree_DiMuon_Rec_FullMass[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_FullMass_%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
        Tree_DiMuon_Rec_FullMass[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_FullMass, "m/D");
        Tree_DiMuon_Rec_FullMass[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_FullMass, "pt/D");

        Tree_DiMuon_Rec_LowMass[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_LowMass_%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
        Tree_DiMuon_Rec_LowMass[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_LowMass, "m/D");
        Tree_DiMuon_Rec_LowMass[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_LowMass, "pt/D");

        Tree_DiMuon_Rec_LowMass_LowPt[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_LowMass_LowPt_%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
        Tree_DiMuon_Rec_LowMass_LowPt[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_LowMass_LowPt, "m/D");
        Tree_DiMuon_Rec_LowMass_LowPt[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_LowMass_LowPt, "pt/D");

        Tree_DiMuon_Rec_HighMass[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_HighMass_%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
        Tree_DiMuon_Rec_HighMass[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_HighMass, "m/D");
        Tree_DiMuon_Rec_HighMass[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_HighMass, "pt/D");

        if (Generator.Contains("Powheg"))
        {
            Tree_DiMuon_Rec_PowhegOnly[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_PowhegOnly_%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev Powheg Onlys");
            Tree_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_PowhegOnly, "m/D");
            Tree_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_PowhegOnly, "pt/D");

            Tree_DiMuon_Rec_FullMass_PowhegOnly[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_FullMass_PowhegOnly%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
            Tree_DiMuon_Rec_FullMass_PowhegOnly[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_FullMass_PowhegOnly, "m/D");
            Tree_DiMuon_Rec_FullMass_PowhegOnly[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_FullMass_PowhegOnly, "pt/D");

            Tree_DiMuon_Rec_LowMass_PowhegOnly[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_LowMass_PowhegOnly%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
            Tree_DiMuon_Rec_LowMass_PowhegOnly[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_LowMass_PowhegOnly, "m/D");
            Tree_DiMuon_Rec_LowMass_PowhegOnly[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_LowMass_PowhegOnly, "pt/D");

            Tree_DiMuon_Rec_LowMass_LowPt_PowhegOnly[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_LowMass_LowPt_PowhegOnly%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
            Tree_DiMuon_Rec_LowMass_LowPt_PowhegOnly[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_LowMass_LowPt_PowhegOnly, "m/D");
            Tree_DiMuon_Rec_LowMass_LowPt_PowhegOnly[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_LowMass_LowPt_PowhegOnly, "pt/D");

            Tree_DiMuon_Rec_HighMass_PowhegOnly[i_DiMuon_origin] = new TTree(Form("DiMuon_Rec_HighMass_PowhegOnly%s", DiMuon_origin[i_DiMuon_origin].Data()), "Dimuons with mass > 4 Gev");
            Tree_DiMuon_Rec_HighMass_PowhegOnly[i_DiMuon_origin]->Branch("m", M_Dimu_Rec_HighMass_PowhegOnly, "m/D");
            Tree_DiMuon_Rec_HighMass_PowhegOnly[i_DiMuon_origin]->Branch("pt", Pt_Dimu_Rec_HighMass_PowhegOnly, "pt/D");
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
        if (i_Event % (Int_t)(total_entries * 0.125) == 0)
        {
            progress_status(i_Event, total_entries);
        }
        input_tree->GetEntry(i_Event);
        Int_t n_Muon_Gen[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_PowhegOnly[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_PYTHIAOnly[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_GeantOnly[n_Muon_origin] = {0};

        Int_t n_Muon_Gen_DQcut[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_DQcut_PowhegOnly[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_DQcut_PYTHIAOnly[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_DQcut_GeantOnly[n_Muon_origin] = {0};

        Int_t n_Muon_Rec[n_Muon_origin] = {0};
        Int_t n_Muon_Rec_PowhegOnly[n_Muon_origin] = {0};
        Int_t n_Muon_Rec_PYTHIAOnly[n_Muon_origin] = {0};
        Int_t n_Muon_Rec_GeantOnly[n_Muon_origin] = {0};

        Int_t n_DiMuon_Gen[n_DiMuon_origin] = {0};
        Int_t n_DiMuon_Gen_DQcut[n_DiMuon_origin] = {0};
        Int_t n_DiMuon_Rec[n_DiMuon_origin] = {0};

        Int_t n_DiMuon_Gen_PowhegOnly[n_DiMuon_origin] = {0};
        Int_t n_DiMuon_Gen_DQcut_PowhegOnly[n_DiMuon_origin] = {0};
        Int_t n_DiMuon_Rec_PowhegOnly[n_DiMuon_origin] = {0};

        Double_t Pt_Gamma = fPt_gamma[0];
        Double_t Y_Gamma = fY_gamma[0];
        Double_t M_Gamma = fM_gamma[0];

        Int_t N_Charm_event = 0;
        Int_t N_Charm_event_fwd = 0;
        Int_t N_Charm_event_PowhegOnly = 0;
        Int_t N_Charm_event_fwd_PowhegOnly = 0;

        Int_t N_Beauty_event = 0;
        Int_t N_Beauty_event_fwd = 0;
        Int_t N_Beauty_event_PowhegOnly = 0;
        Int_t N_Beauty_event_fwd_PowhegOnly = 0;

        h_PtM_Gamma->Fill(Pt_Gamma, M_Gamma);
        h_PtY_Gamma->Fill(Pt_Gamma, Y_Gamma);

        for (Int_t i_N_HFquarks_gen = 0; i_N_HFquarks_gen < N_HFquarks_gen; i_N_HFquarks_gen++)
        {
            Int_t PDG_HFquark = PDG_HFquark_gen[i_N_HFquarks_gen];
            Double_t Pt_HFquark = Pt_HFquark_gen[i_N_HFquarks_gen];
            Double_t Y_HFquark = Y_HFquark_gen[i_N_HFquarks_gen];
            Int_t HF_quark_origin = Mother_index[i_N_HFquarks_gen];

            if (TMath::Abs(PDG_HFquark) == 4)
            {
                h_PtY_Charm_quark->Fill(Pt_HFquark, Y_HFquark);
                N_Charm_event++;
                if (Generator.Contains("Powheg") && HF_quark_origin == -1)
                {
                    h_PtY_Charm_quark_PowhegOnly->Fill(Pt_HFquark, Y_HFquark);
                    N_Charm_event_PowhegOnly++;
                }

                if (Y_HFquark > -4.0 && Y_HFquark < -2.5)
                {
                    N_Charm_event_fwd++;
                    if (Generator.Contains("Powheg") && HF_quark_origin == -1)
                        N_Charm_event_fwd_PowhegOnly++;
                }
            }
            else if (TMath::Abs(PDG_HFquark) == 5)
            {
                h_PtY_Beauty_quark->Fill(Pt_HFquark, Y_HFquark);
                N_Beauty_event++;
                if (Generator.Contains("Powheg") && HF_quark_origin == -1)
                {
                    h_PtY_Beauty_quark_PowhegOnly->Fill(Pt_HFquark, Y_HFquark);
                    N_Beauty_event_PowhegOnly++;
                }

                if (Y_HFquark > -4.0 && Y_HFquark < -2.5)
                {
                    N_Beauty_event_fwd++;
                    if (Generator.Contains("Powheg") && HF_quark_origin == -1)
                        N_Beauty_event_fwd_PowhegOnly++;
                }
            }
        }
        h_NCharm_event->Fill(N_Charm_event);
        h_NCharm_event_fwd->Fill(N_Charm_event_fwd);
        h_NCharm_event_PowhegOnly->Fill(N_Charm_event_PowhegOnly);
        h_NCharm_event_fwd_PowhegOnly->Fill(N_Charm_event_fwd_PowhegOnly);

        h_NBeauty_event->Fill(N_Beauty_event);
        h_NBeauty_event_fwd->Fill(N_Beauty_event_fwd);
        h_NBeauty_event_PowhegOnly->Fill(N_Beauty_event_PowhegOnly);
        h_NBeauty_event_fwd_PowhegOnly->Fill(N_Beauty_event_fwd_PowhegOnly);

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
            Double_t VzHadron = fVzHadron_gen[i_NHadron_gen];
            Int_t IsHadronfromGeant = fHadronFrom_Geant_gen[i_NHadron_gen];

            // if ((TMath::Abs(PDG_Hadron) > 100 && TMath::Abs(PDG_Hadron) < 400) || (TMath::Abs(PDG_Hadron) > 1000 && TMath::Abs(PDG_Hadron) < 4000))
            // {
            //     if (IsHadronfromGeant == -1)
            //         h_VzHadronEta_Hadron_Gen_LF_Geant->Fill(VzHadron, Eta_Hadron);
            //     else
            //         h_VzHadronEta_Hadron_Gen_LF_PYTHIA->Fill(VzHadron, Eta_Hadron);
            // }
            // else if ((TMath::Abs(PDG_Hadron) > 400 && TMath::Abs(PDG_Hadron) < 500) || (TMath::Abs(PDG_Hadron) > 4000 && TMath::Abs(PDG_Hadron) < 5000))
            // {
            //     if ((PDGmum_Hadron == 0) || (PDGmum_Hadron == 4) || (TMath::Abs(PDGmum_Hadron) > 400 && TMath::Abs(PDGmum_Hadron) < 500) || (TMath::Abs(PDGmum_Hadron) > 4000 && TMath::Abs(PDGmum_Hadron) < 5000))
            //     {
            //         if (IsHadronfromGeant == -1)
            //             h_VzHadronEta_Hadron_Gen_Charm_Geant->Fill(VzHadron, Eta_Hadron);
            //         else
            //             h_VzHadronEta_Hadron_Gen_Charm_PYTHIA->Fill(VzHadron, Eta_Hadron);
            //     }
            // }
            // else if ((TMath::Abs(PDG_Hadron) > 500 && TMath::Abs(PDG_Hadron) < 600) || (TMath::Abs(PDG_Hadron) > 5000 && TMath::Abs(PDG_Hadron) < 6000))
            // {
            //     if (IsHadronfromGeant == -1)
            //         h_VzHadronEta_Hadron_Gen_Beauty_Geant->Fill(VzHadron, Eta_Hadron);
            //     else
            //         h_VzHadronEta_Hadron_Gen_Beauty_PYTHIA->Fill(VzHadron, Eta_Hadron);
            // }

            if ((TMath::Abs(PDG_Hadron) > 400 && TMath::Abs(PDG_Hadron) < 600) || (TMath::Abs(PDG_Hadron) > 4000 && TMath::Abs(PDG_Hadron) < 6000))
            {
                if ((PDGmum_Hadron == 0) || (PDGmum_Hadron == 4) || (TMath::Abs(PDGmum_Hadron) > 400 && TMath::Abs(PDGmum_Hadron) < 500) || (TMath::Abs(PDGmum_Hadron) > 4000 && TMath::Abs(PDGmum_Hadron) < 5000))
                {
                    h_PdgPtY_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron, Y_Hadron);

                    h_PdgPt_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                    h_PdgY_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                    h_PdgEta_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);

                    if (Generator.Contains("Geant"))
                    {
                        if (IsHadronfromGeant == -1)
                        {
                            h_PdgY_HFHadron_prompt_GeantOnly->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                            h_PdgEta_HFHadron_prompt_GeantOnly->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                            h_PdgPt_HFHadron_prompt_GeantOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                        }
                        else
                        {
                            if (Generator.Contains("Powheg") && IsHadronfromPowheg == -1)
                            {

                                h_PdgY_HFHadron_prompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                                h_PdgEta_HFHadron_prompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                                h_PdgPt_HFHadron_prompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                            }
                            else
                            {
                                h_PdgY_HFHadron_prompt_PYTHIAOnly->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                                h_PdgEta_HFHadron_prompt_PYTHIAOnly->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                                h_PdgPt_HFHadron_prompt_PYTHIAOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                            }
                        }
                    }
                }
                else
                {
                    h_PdgPtY_HFHadron_notprompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron, Y_Hadron);

                    h_PdgPt_HFHadron_notprompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                    h_PdgY_HFHadron_notprompt->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                    h_PdgEta_HFHadron_notprompt->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                    if (Generator.Contains("Geant"))
                    {
                        if (IsHadronfromGeant == -1)
                        {
                            h_PdgY_HFHadron_notprompt_GeantOnly->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                            h_PdgEta_HFHadron_notprompt_GeantOnly->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                            h_PdgPt_HFHadron_notprompt_GeantOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                        }
                        else
                        {
                            if (Generator.Contains("Powheg") && IsHadronfromPowheg == -1)
                            {

                                h_PdgY_HFHadron_notprompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                                h_PdgEta_HFHadron_notprompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                                h_PdgPt_HFHadron_notprompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                            }
                            else
                            {
                                h_PdgY_HFHadron_notprompt_PYTHIAOnly->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                                h_PdgEta_HFHadron_notprompt_PYTHIAOnly->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                                h_PdgPt_HFHadron_notprompt_PYTHIAOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                            }
                        }
                    }
                }
            }
            else
            {
                h_PdgPtY_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron, Y_Hadron);

                h_PdgPt_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                h_PdgY_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                h_PdgEta_HFHadron_prompt->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                if (Generator.Contains("Geant"))
                {
                    if (IsHadronfromGeant == -1)
                    {
                        h_PdgY_HFHadron_prompt_GeantOnly->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                        h_PdgEta_HFHadron_prompt_GeantOnly->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                        h_PdgPt_HFHadron_prompt_GeantOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                    }
                    else
                    {
                        if (Generator.Contains("Powheg") && IsHadronfromPowheg == -1)
                        {

                            h_PdgY_HFHadron_prompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                            h_PdgEta_HFHadron_prompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                            h_PdgPt_HFHadron_prompt_PowhegOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                        }
                        else
                        {
                            h_PdgY_HFHadron_prompt_PYTHIAOnly->Fill(TMath::Abs(PDG_Hadron), Y_Hadron);
                            h_PdgEta_HFHadron_prompt_PYTHIAOnly->Fill(TMath::Abs(PDG_Hadron), Eta_Hadron);
                            h_PdgPt_HFHadron_prompt_PYTHIAOnly->Fill(TMath::Abs(PDG_Hadron), Pt_Hadron);
                        }
                    }
                }
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
            Double_t Radius_gen = 999;
            Double_t Vz_gen = 999;
            Double_t Vzmother_gen = 999;
            Int_t IsFrom_Geant_gen = 999;

            if (Generator.Contains("Geant"))
            {
                Radius_gen = fRadius_gen[i_NMuons_gen] / 1e+03;
                Vz_gen = fVz_gen[i_NMuons_gen];
                Vzmother_gen = fVzmother_gen[i_NMuons_gen];
                IsFrom_Geant_gen = fFrom_Geant_gen[i_NMuons_gen];
            }

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

            h_PtPdg_Muon_Gen->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
            h_YPdg_Muon_Gen->Fill(Y_Mu, TMath::Abs(PDG_Mu));
            h_EtaPdg_Muon_Gen->Fill(Eta_Mu, TMath::Abs(PDG_Mu));

            // Separating PYTHIA, GEANT and POWHEG (if necessary) components for generated muons in centralized simulations
            if (Generator.Contains("Geant"))
            {
                if (IsFrom_Geant_gen == -1)
                {
                    h_PtPdg_Muon_Gen_GeantOnly->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                    h_YPdg_Muon_Gen_GeantOnly->Fill(Y_Mu, TMath::Abs(PDG_Mu));
                    h_EtaPdg_Muon_Gen_GeantOnly->Fill(Eta_Mu, TMath::Abs(PDG_Mu));
                }
                else
                {
                    if (Generator.Contains("Powheg") && (IsFrom_Powheg_gen == -1))
                    {

                        h_PtPdg_Muon_Gen_PowhegOnly->Fill(Pt_Mu, PDG_Mu);
                        h_YPdg_Muon_Gen_PowhegOnly->Fill(Y_Mu, PDG_Mu);
                        // h_EtaPdg_Muon_Gen_PowhegOnly->Fill(Y_Mu, PDG_Mu);
                    }
                    else
                    {
                        h_PtPdg_Muon_Gen_PYTHIAOnly->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                        h_YPdg_Muon_Gen_PYTHIAOnly->Fill(Y_Mu, TMath::Abs(PDG_Mu));
                        h_EtaPdg_Muon_Gen_PYTHIAOnly->Fill(Eta_Mu, TMath::Abs(PDG_Mu));
                    }
                }
            }

            if (Selection_Muon[3] && (Generator.Contains("Geant")))
            {
                h_RadiusEta_Muon_Gen_LF->Fill(Radius_gen, Eta_Mu);
                h_VzEta_Muon_Gen_LF->Fill(Vz_gen, Eta_Mu);
                if (Generator.Contains("Powheg"))
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
            }

            for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
            {
                if (Selection_Muon[i_Muon_origin])
                {
                    h_PtY_Muon_Gen[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    h_PtEta_Muon_Gen[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                    n_Muon_Gen[i_Muon_origin]++;

                    if (!(Generator.Contains("Geant")))
                        continue;

                    if (IsFrom_Geant_gen == -1)
                    {
                        h_VzmotherEta_Muon_Gen_Geant[i_Muon_origin]->Fill(Vzmother_gen, Eta_Mu);
                        h_PtY_Muon_Gen_GeantOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                        h_PtEta_Muon_Gen_GeantOnly[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                        n_Muon_Gen_GeantOnly[i_Muon_origin]++;
                    }
                    else
                    {
                        if (Generator.Contains("Powheg") && (IsFrom_Powheg_gen == -1))
                        {
                            h_PtY_Muon_Gen_PowhegOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                            h_PtEta_Muon_Gen_PowhegOnly[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                            n_Muon_Gen_PowhegOnly[i_Muon_origin]++;
                        }
                        else
                        {
                            h_VzmotherEta_Muon_Gen_PYTHIA[i_Muon_origin]->Fill(Vzmother_gen, Eta_Mu);
                            h_PtY_Muon_Gen_PYTHIAOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                            h_PtEta_Muon_Gen_PYTHIAOnly[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                            n_Muon_Gen_PYTHIAOnly[i_Muon_origin]++;
                        }
                    }
                }
            }

            if ((Eta_Mu > -4.0 && Eta_Mu < -2.5) && Pt_Mu > 0.5)
                DQ_Muon = kTRUE;
            if (!DQ_Muon)
                continue;

            // Separating PYTHIA, GEANT and POWHEG (if necessary) components for generated muons with DQ cuts in centralized simulations

            h_PtPdg_Muon_Gen_DQcut->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
            h_YPdg_Muon_Gen_DQcut->Fill(Y_Mu, TMath::Abs(PDG_Mu));
            h_EtaPdg_Muon_Gen_DQcut->Fill(Eta_Mu, TMath::Abs(PDG_Mu));

            if (Generator.Contains("Geant"))
            {
                if (IsFrom_Geant_gen == -1)
                {
                    h_PtPdg_Muon_Gen_DQcut_GeantOnly->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                    h_YPdg_Muon_Gen_DQcut_GeantOnly->Fill(Y_Mu, TMath::Abs(PDG_Mu));
                    h_EtaPdg_Muon_Gen_DQcut_GeantOnly->Fill(Eta_Mu, TMath::Abs(PDG_Mu));
                }
                else
                {
                    if (Generator.Contains("Powheg") && (IsFrom_Powheg_gen == -1))
                    {

                        h_PtPdg_Muon_Gen_DQcut_PowhegOnly->Fill(Pt_Mu, PDG_Mu);
                        h_YPdg_Muon_Gen_DQcut_PowhegOnly->Fill(Y_Mu, PDG_Mu);
                    }
                    else
                    {
                        h_PtPdg_Muon_Gen_DQcut_PYTHIAOnly->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                        h_YPdg_Muon_Gen_DQcut_PYTHIAOnly->Fill(Y_Mu, TMath::Abs(PDG_Mu));
                        h_EtaPdg_Muon_Gen_DQcut_PYTHIAOnly->Fill(Eta_Mu, TMath::Abs(PDG_Mu));
                    }
                }
            }

            for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
            {
                if (Selection_Muon[i_Muon_origin])
                {
                    h_PtY_Muon_Gen_DQcut[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    h_PtEta_Muon_Gen_DQcut[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                    n_Muon_Gen_DQcut[i_Muon_origin]++;

                    if (Generator.Contains("Geant"))
                    {
                        if (IsFrom_Geant_gen == -1)
                        {
                            n_Muon_Gen_DQcut_GeantOnly[i_Muon_origin]++;
                            h_PtY_Muon_Gen_DQcut_GeantOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                            h_PtEta_Muon_Gen_DQcut_GeantOnly[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                        }
                        else
                        {
                            if (Generator.Contains("Powheg") && (IsFrom_Powheg_gen == -1))
                            {
                                n_Muon_Gen_DQcut_PowhegOnly[i_Muon_origin]++;
                                h_PtY_Muon_Gen_DQcut_PowhegOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                                h_PtEta_Muon_Gen_DQcut_PowhegOnly[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                            }
                            else
                            {
                                n_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_origin]++;
                                h_PtY_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                                h_PtEta_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                            }
                        }
                    }
                }
            }
        }

        for (Int_t i_Muon_Origin = 0; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
        {
            h_Nperevent_Muon_Gen[i_Muon_Origin]->Fill(n_Muon_Gen[i_Muon_Origin]);
            h_Nperevent_Muon_Gen_DQcut[i_Muon_Origin]->Fill(n_Muon_Gen_DQcut[i_Muon_Origin]);
            if (Generator.Contains("Geant"))
            {
                h_Nperevent_Muon_Gen_PYTHIAOnly[i_Muon_Origin]->Fill(n_Muon_Gen_PYTHIAOnly[i_Muon_Origin]);
                h_Nperevent_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_Origin]->Fill(n_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_Origin]);

                h_Nperevent_Muon_Gen_GeantOnly[i_Muon_Origin]->Fill(n_Muon_Gen_GeantOnly[i_Muon_Origin]);
                h_Nperevent_Muon_Gen_DQcut_GeantOnly[i_Muon_Origin]->Fill(n_Muon_Gen_DQcut_GeantOnly[i_Muon_Origin]);

                if (Generator.Contains("Powheg"))
                {
                    h_Nperevent_Muon_Gen_PowhegOnly[i_Muon_Origin]->Fill(n_Muon_Gen_PowhegOnly[i_Muon_Origin]);
                    h_Nperevent_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]->Fill(n_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]);
                }
            }
        }

        for (Int_t i_NMuons_rec = 0; i_NMuons_rec < NMuons_rec; i_NMuons_rec++)
        {
            if (!(Generator.Contains("Geant")))
                continue;

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
            Double_t Vzmother_rec = fVzmother_rec[i_NMuons_rec];
            Int_t IsFrom_Geant_rec = fFrom_Geant_rec[i_NMuons_rec];

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
            h_PtPdg_Muon_Rec->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
            h_YPdg_Muon_Rec->Fill(Y_Mu, TMath::Abs(PDG_Mu));
            h_EtaPdg_Muon_Rec->Fill(Eta_Mu, TMath::Abs(PDG_Mu));

            if (Generator.Contains("Geant"))
            {
                if (IsFrom_Geant_rec == -1)
                {
                    h_PtPdg_Muon_Rec_GeantOnly->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                    h_YPdg_Muon_Rec_GeantOnly->Fill(Y_Mu, TMath::Abs(PDG_Mu));
                    h_EtaPdg_Muon_Rec_GeantOnly->Fill(Eta_Mu, TMath::Abs(PDG_Mu));
                }
                else
                {
                    if (Generator.Contains("Powheg") && (IsFromPowheg == -1))
                    {
                        h_PtPdg_Muon_Rec_PowhegOnly->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                        h_YPdg_Muon_Rec_PowhegOnly->Fill(Y_Mu, TMath::Abs(PDG_Mu));
                        h_EtaPdg_Muon_Rec_PowhegOnly->Fill(Eta_Mu, TMath::Abs(PDG_Mu));
                    }
                    else
                    {
                        h_PtPdg_Muon_Rec_PYTHIAOnly->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                        h_YPdg_Muon_Rec_PYTHIAOnly->Fill(Y_Mu, TMath::Abs(PDG_Mu));
                        h_EtaPdg_Muon_Rec_PYTHIAOnly->Fill(Eta_Mu, TMath::Abs(PDG_Mu));
                    }
                }
            }

            for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
            {
                if (Selection_Muon[i_Muon_origin])
                {
                    h_PtY_Muon_Rec[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    h_PtEta_Muon_Rec[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                    n_Muon_Rec[i_Muon_origin]++;
                    if (!(Generator.Contains("Geant")))
                        continue;

                    if (IsFrom_Geant_rec == -1)
                    {
                        h_VzmotherEta_Muon_Rec_Geant[i_Muon_origin]->Fill(Vzmother_rec, Eta_Mu);
                        h_PtY_Muon_Rec_GeantOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                        h_PtEta_Muon_Rec_GeantOnly[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                        n_Muon_Rec_GeantOnly[i_Muon_origin]++;
                    }
                    else
                    {

                        if (Generator.Contains("Powheg") && (IsFromPowheg == -1))
                        {

                            h_PtY_Muon_Rec_PowhegOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                            h_PtEta_Muon_Rec_PowhegOnly[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                            n_Muon_Rec_PowhegOnly[i_Muon_origin]++;
                        }
                        else
                        {
                            h_VzmotherEta_Muon_Rec_PYTHIA[i_Muon_origin]->Fill(Vzmother_rec, Eta_Mu);
                            h_PtY_Muon_Rec_PYTHIAOnly[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                            h_PtEta_Muon_Rec_PYTHIAOnly[i_Muon_origin]->Fill(Pt_Mu, Eta_Mu);
                            n_Muon_Rec_PYTHIAOnly[i_Muon_origin]++;
                        }
                    }
                }
            }
        }

        for (Int_t i_Muon_Origin = 0; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
        {
            h_Nperevent_Muon_Rec[i_Muon_Origin]->Fill(n_Muon_Rec[i_Muon_Origin]);
            if (!(Generator.Contains("Geant")))
                continue;
            h_Nperevent_Muon_Rec_GeantOnly[i_Muon_Origin]->Fill(n_Muon_Rec_GeantOnly[i_Muon_Origin]);
            h_Nperevent_Muon_Rec_PYTHIAOnly[i_Muon_Origin]->Fill(n_Muon_Rec_PYTHIAOnly[i_Muon_Origin]);
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
            Int_t IsFromPowheg_Mu0 = 999;

            Double_t Y_Mu1 = Y_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Pt_Mu1 = Pt_rec[DimuMu_rec[i_NDimu_rec][1]];
            Int_t PDG_Mu1 = PDGmum_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Charge_Mu1 = Charge_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t RAbs_Mu1 = RAtAbsEnd_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t pDCA_Mu1 = pDCA_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Eta_Mu1 = Eta_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Phi_Mu1 = Phi_rec[DimuMu_rec[i_NDimu_rec][1]];
            Int_t IsFromPowheg_Mu1 = 999;

            Int_t IsFrom_Geant_gen_Mu0 = 999;
            Int_t IsFrom_Geant_gen_Mu1 = 999;

            Bool_t LF_Generator[n_LF_DiMuon_Generator] = {kFALSE};

            if (Generator.Contains("Geant"))
            {
                IsFrom_Geant_gen_Mu0 = fFrom_Geant_gen[DimuMu_rec[i_NDimu_rec][0]];
                IsFrom_Geant_gen_Mu1 = fFrom_Geant_gen[DimuMu_rec[i_NDimu_rec][1]];
                if (IsFrom_Geant_gen_Mu0 == -1 && IsFrom_Geant_gen_Mu1 == -1)
                    LF_Generator[0] = kTRUE;
                else if ((IsFrom_Geant_gen_Mu0 == -1 && IsFrom_Geant_gen_Mu1 == 0) || (IsFrom_Geant_gen_Mu1 == -1 && IsFrom_Geant_gen_Mu0 == 0))
                    LF_Generator[1] = kTRUE;
                else
                    LF_Generator[2] = kTRUE;
                if (Generator.Contains("Powheg"))
                {
                    IsFromPowheg_Mu0 = fFrom_Powheg_rec[DimuMu_rec[i_NDimu_rec][0]];
                    IsFromPowheg_Mu1 = fFrom_Powheg_rec[DimuMu_rec[i_NDimu_rec][1]];
                }
            }

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

                    if (i_DiMuon_origin == 3) // selecting LF
                    {
                        for (Int_t i_LF_Generator = 0; i_LF_Generator < n_LF_DiMuon_Generator; i_LF_Generator++)
                        {
                            if (LF_Generator[i_LF_Generator])
                            {
                                h_PtM_DiMuon_Rec_fromLF[i_LF_Generator]->Fill(Pt_DiMu, M_DiMu);
                                h_PtY_DiMuon_Rec_fromLF[i_LF_Generator]->Fill(Pt_DiMu, Y_DiMu);
                            }
                        }
                    }

                    if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                    {
                        h_PtM_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                        h_PtY_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                        n_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]++;
                    }
                    if (M_DiMu > 4 && M_DiMu < 30)
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
                    if (M_DiMu > 4)
                    {
                        *Pt_Dimu_Rec_FullMass = Pt_DiMu;
                        *M_Dimu_Rec_FullMass = M_DiMu;
                        Tree_DiMuon_Rec_FullMass[i_DiMuon_origin]->Fill();

                        if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                        {
                            *Pt_Dimu_Rec_FullMass_PowhegOnly = Pt_DiMu;
                            *M_Dimu_Rec_FullMass_PowhegOnly = M_DiMu;
                            Tree_DiMuon_Rec_FullMass_PowhegOnly[i_DiMuon_origin]->Fill();
                        }
                    }
                    if (M_DiMu > 4 && M_DiMu < 9)
                    {
                        *Pt_Dimu_Rec_LowMass = Pt_DiMu;
                        *M_Dimu_Rec_LowMass = M_DiMu;
                        Tree_DiMuon_Rec_LowMass[i_DiMuon_origin]->Fill();

                        if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                        {
                            *Pt_Dimu_Rec_LowMass_PowhegOnly = Pt_DiMu;
                            *M_Dimu_Rec_LowMass_PowhegOnly = M_DiMu;
                            Tree_DiMuon_Rec_LowMass_PowhegOnly[i_DiMuon_origin]->Fill();
                        }

                        if (Pt_DiMu < 10)
                        {
                            *Pt_Dimu_Rec_LowMass_LowPt = Pt_DiMu;
                            *M_Dimu_Rec_LowMass_LowPt = M_DiMu;
                            Tree_DiMuon_Rec_LowMass_LowPt[i_DiMuon_origin]->Fill();

                            if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                            {
                                *Pt_Dimu_Rec_LowMass_LowPt_PowhegOnly = Pt_DiMu;
                                *M_Dimu_Rec_LowMass_LowPt_PowhegOnly = M_DiMu;
                                Tree_DiMuon_Rec_LowMass_LowPt_PowhegOnly[i_DiMuon_origin]->Fill();
                            }
                        }
                    }
                    if (M_DiMu > 11 && M_DiMu < 30)
                    {
                        *Pt_Dimu_Rec_HighMass = Pt_DiMu;
                        *M_Dimu_Rec_HighMass = M_DiMu;
                        Tree_DiMuon_Rec_HighMass[i_DiMuon_origin]->Fill();

                        if (Generator.Contains("Powheg") && IsFromPowheg_Mu0 == -1 && IsFromPowheg_Mu1 == -1)
                        {
                            *Pt_Dimu_Rec_HighMass_PowhegOnly = Pt_DiMu;
                            *M_Dimu_Rec_HighMass_PowhegOnly = M_DiMu;
                            Tree_DiMuon_Rec_HighMass_PowhegOnly[i_DiMuon_origin]->Fill();
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
    // printf("counter_test %i\n", counter_test);
    TFile fOut_Tree(Form("%s/%s", dir_fileOut.Data(), file_out_tree.Data()), "RECREATE");
    fOut_Tree.cd();
    h_Nevents->Write();

    for (Int_t i_Dimuon_origin = 0; i_Dimuon_origin < n_DiMuon_origin; i_Dimuon_origin++)
    {
        if (Tree_DiMuon_Rec[i_Dimuon_origin]->GetEntries() > 0)
            Tree_DiMuon_Rec[i_Dimuon_origin]->Write(0, 2, 0);

        if (Tree_DiMuon_Rec_FullMass[i_Dimuon_origin]->GetEntries() > 0)
            Tree_DiMuon_Rec_FullMass[i_Dimuon_origin]->Write(0, 2, 0);

        if (Tree_DiMuon_Rec_LowMass[i_Dimuon_origin]->GetEntries() > 0)
            Tree_DiMuon_Rec_LowMass[i_Dimuon_origin]->Write(0, 2, 0);

        if (Tree_DiMuon_Rec_LowMass_LowPt[i_Dimuon_origin]->GetEntries() > 0)
            Tree_DiMuon_Rec_LowMass_LowPt[i_Dimuon_origin]->Write(0, 2, 0);

        if (Tree_DiMuon_Rec_HighMass[i_Dimuon_origin]->GetEntries() > 0)
            Tree_DiMuon_Rec_HighMass[i_Dimuon_origin]->Write(0, 2, 0);
        if (Generator.Contains("Powheg"))
        {
            if (Tree_DiMuon_Rec_PowhegOnly[i_Dimuon_origin]->GetEntries() > 0)
                Tree_DiMuon_Rec_PowhegOnly[i_Dimuon_origin]->Write(0, 2, 0);

            if (Tree_DiMuon_Rec_FullMass_PowhegOnly[i_Dimuon_origin]->GetEntries() > 0)
                Tree_DiMuon_Rec_FullMass_PowhegOnly[i_Dimuon_origin]->Write(0, 2, 0);

            if (Tree_DiMuon_Rec_LowMass_PowhegOnly[i_Dimuon_origin]->GetEntries() > 0)
                Tree_DiMuon_Rec_LowMass_PowhegOnly[i_Dimuon_origin]->Write(0, 2, 0);

            if (Tree_DiMuon_Rec_LowMass_LowPt_PowhegOnly[i_Dimuon_origin]->GetEntries() > 0)
                Tree_DiMuon_Rec_LowMass_LowPt_PowhegOnly[i_Dimuon_origin]->Write(0, 2, 0);

            if (Tree_DiMuon_Rec_HighMass_PowhegOnly[i_Dimuon_origin]->GetEntries() > 0)
                Tree_DiMuon_Rec_HighMass_PowhegOnly[i_Dimuon_origin]->Write(0, 2, 0);
        }
    }
    // Tree_DiMuon_Rec->Write(0, 2, 0);
    fOut_Tree.Close();

    TFile fOut(Form("%s/%s", dir_fileOut.Data(), file_out.Data()), "RECREATE");
    fOut.cd();
    h_Nevents->Write();
    if (Generator.Contains("HF"))
    {
        if (!fOut.GetDirectory("HF_quarks"))
            fOut.mkdir("HF_quarks");

        fOut.cd("HF_quarks");

        h_PtY_Charm_quark->Write(0, 2, 0);
        h_NCharm_event->Write(0, 2, 0);
        h_NCharm_event_fwd->Write(0, 2, 0);

        h_PtY_Beauty_quark->Write(0, 2, 0);
        h_NBeauty_event->Write(0, 2, 0);
        h_NBeauty_event_fwd->Write(0, 2, 0);

        if (Generator.Contains("Powheg"))
        {
            if (!fOut.GetDirectory("HF_quarks/Powheg"))
                fOut.mkdir("HF_quarks/Powheg");

            fOut.cd("HF_quarks/Powheg");

            h_PtY_Charm_quark_PowhegOnly->Write(0, 2, 0);
            h_NCharm_event_PowhegOnly->Write(0, 2, 0);
            h_NCharm_event_fwd_PowhegOnly->Write(0, 2, 0);

            h_PtY_Beauty_quark_PowhegOnly->Write(0, 2, 0);
            h_NBeauty_event_PowhegOnly->Write(0, 2, 0);
            h_NBeauty_event_fwd_PowhegOnly->Write(0, 2, 0);
        }
    }

    if (!fOut.GetDirectory("Muon_Gen"))
        fOut.mkdir("Muon_Gen");

    if (!fOut.GetDirectory("Muon_Rec"))
        fOut.mkdir("Muon_Rec");

    if (!fOut.GetDirectory("DiMuon_Gen"))
        fOut.mkdir("DiMuon_Gen");

    if (!fOut.GetDirectory("DiMuon_Rec"))
        fOut.mkdir("DiMuon_Rec");

    if (Generator.Contains("Geant"))
    {
        if (!fOut.GetDirectory("Muon_Gen/Geant"))
            fOut.mkdir("Muon_Gen/Geant");

        if (!fOut.GetDirectory("Muon_Gen/PYTHIA"))
            fOut.mkdir("Muon_Gen/PYTHIA");

        if (!fOut.GetDirectory("Muon_Rec/Geant"))
            fOut.mkdir("Muon_Rec/Geant");

        if (!fOut.GetDirectory("Muon_Rec/PYTHIA"))
            fOut.mkdir("Muon_Rec/PYTHIA");

        if (!fOut.GetDirectory("DiMuon_Rec/LF_origin"))
            fOut.mkdir("DiMuon_Rec/LF_origin");

        if (Generator.Contains("Powheg"))
        {
            if (!fOut.GetDirectory("Muon_Gen/Powheg"))
                fOut.mkdir("Muon_Gen/Powheg");
            if (!fOut.GetDirectory("Muon_Rec/Powheg"))
                fOut.mkdir("Muon_Rec/Powheg");

            if (!fOut.GetDirectory("DiMuon_Gen/Powheg"))
                fOut.mkdir("DiMuon_Gen/Powheg");
            if (!fOut.GetDirectory("DiMuon_Rec/Powheg"))
                fOut.mkdir("DiMuon_Rec/Powheg");
        }
    }

    fOut.cd("Muon_Gen");
    h_PtPdg_Muon_Gen->Write(0, 2, 0);
    h_YPdg_Muon_Gen->Write(0, 2, 0);
    h_EtaPdg_Muon_Gen->Write(0, 2, 0);

    h_PtPdg_Muon_Gen_DQcut->Write(0, 2, 0);
    h_YPdg_Muon_Gen_DQcut->Write(0, 2, 0);
    h_EtaPdg_Muon_Gen_DQcut->Write(0, 2, 0);

    h_RadiusEta_Muon_Gen_LF->Write(0, 2, 0);
    h_VzEta_Muon_Gen_LF->Write(0, 2, 0);

    if (fOut.GetDirectory("Muon_Gen/Powheg"))
    {
        fOut.cd("Muon_Gen/Powheg");
        if (h_PtQ_Muon_Gen_LF_Powheg->GetEntries())
        {
            h_PtQ_Muon_Gen_LF_Powheg->Write();
            h_PtQ_Muon_Gen_LF_Pythia->Write();
            h_YQ_Muon_Gen_LF_Powheg->Write();
            h_YQ_Muon_Gen_LF_Pythia->Write();
        }
        if (h_PtPdg_Muon_Gen_PowhegOnly->GetEntries() > 0.0)
        {
            h_PtPdg_Muon_Gen_PowhegOnly->Write(0, 2, 0);
            h_YPdg_Muon_Gen_PowhegOnly->Write(0, 2, 0);
            h_PtPdg_Muon_Gen_DQcut_PowhegOnly->Write(0, 2, 0);
            h_YPdg_Muon_Gen_DQcut_PowhegOnly->Write(0, 2, 0);
        }
    }

    if (fOut.GetDirectory("Muon_Gen/Geant"))
    {
        fOut.cd("Muon_Gen/Geant");
        if (h_PtPdg_Muon_Gen_GeantOnly->GetEntries() > 0.0)
        {
            h_PtPdg_Muon_Gen_GeantOnly->Write(0, 2, 0);
            h_YPdg_Muon_Gen_GeantOnly->Write(0, 2, 0);
            h_EtaPdg_Muon_Gen_GeantOnly->Write(0, 2, 0);
        }
        if (h_PtPdg_Muon_Gen_DQcut_GeantOnly->GetEntries() > 0.0)
        {
            h_PtPdg_Muon_Gen_DQcut_GeantOnly->Write(0, 2, 0);
            h_YPdg_Muon_Gen_DQcut_GeantOnly->Write(0, 2, 0);
            h_EtaPdg_Muon_Gen_DQcut_GeantOnly->Write(0, 2, 0);
        }
    }

    if (fOut.GetDirectory("Muon_Gen/PYTHIA"))
    {
        fOut.cd("Muon_Gen/PYTHIA");
        if (h_PtPdg_Muon_Gen_PYTHIAOnly->GetEntries() > 0.0)
        {
            h_PtPdg_Muon_Gen_PYTHIAOnly->Write(0, 2, 0);
            h_YPdg_Muon_Gen_PYTHIAOnly->Write(0, 2, 0);
            h_EtaPdg_Muon_Gen_PYTHIAOnly->Write(0, 2, 0);
        }
        if (h_PtPdg_Muon_Gen_DQcut_PYTHIAOnly->GetEntries() > 0.0)
        {
            h_PtPdg_Muon_Gen_DQcut_PYTHIAOnly->Write(0, 2, 0);
            h_YPdg_Muon_Gen_DQcut_PYTHIAOnly->Write(0, 2, 0);
            h_EtaPdg_Muon_Gen_DQcut_PYTHIAOnly->Write(0, 2, 0);
        }
    }

    fOut.cd();
    fOut.cd("Muon_Rec");
    h_PtPdg_Muon_Rec->Write(0, 2, 0);
    h_YPdg_Muon_Rec->Write(0, 2, 0);
    h_EtaPdg_Muon_Rec->Write(0, 2, 0);

    if (fOut.GetDirectory("Muon_Rec/Powheg"))
    {
        fOut.cd("Muon_Rec/Powheg");
        if (h_PtPdg_Muon_Rec_PowhegOnly->GetEntries() > 0.0)
        {
            h_PtPdg_Muon_Rec_PowhegOnly->Write(0, 2, 0);
            h_YPdg_Muon_Rec_PowhegOnly->Write(0, 2, 0);
            h_EtaPdg_Muon_Rec_PowhegOnly->Write(0, 2, 0);
        }
    }

    if (fOut.GetDirectory("Muon_Rec/PYTHIA"))
    {
        fOut.cd("Muon_Rec/PYTHIA");
        if (h_PtPdg_Muon_Rec_PYTHIAOnly->GetEntries() > 0.0)
        {
            h_PtPdg_Muon_Rec_PYTHIAOnly->Write(0, 2, 0);
            h_YPdg_Muon_Rec_PYTHIAOnly->Write(0, 2, 0);
            h_EtaPdg_Muon_Rec_PYTHIAOnly->Write(0, 2, 0);
        }
    }

    if (fOut.GetDirectory("Muon_Rec/Geant"))
    {
        fOut.cd("Muon_Rec/Geant");
        if (h_PtPdg_Muon_Rec_GeantOnly->GetEntries() > 0.0)
        {
            h_PtPdg_Muon_Rec_GeantOnly->Write(0, 2, 0);
            h_YPdg_Muon_Rec_GeantOnly->Write(0, 2, 0);
            h_EtaPdg_Muon_Rec_GeantOnly->Write(0, 2, 0);
        }
    }

    for (Int_t i_Muon_Origin = 0; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
    {
        // Gen Muon saving
        fOut.cd();
        fOut.cd("Muon_Gen");
        if (h_Nperevent_Muon_Gen[i_Muon_Origin]->GetMean() > 0.0)
        {
            h_Nperevent_Muon_Gen[i_Muon_Origin]->Write(0, 2, 0);
            h_PtY_Muon_Gen[i_Muon_Origin]->Write(0, 2, 0);
            h_PtEta_Muon_Gen[i_Muon_Origin]->Write(0, 2, 0);
            h_Nperevent_Muon_Gen_DQcut[i_Muon_Origin]->Write(0, 2, 0);
            h_PtY_Muon_Gen_DQcut[i_Muon_Origin]->Write(0, 2, 0);
        }

        // Gen Muon from Powheg saving

        if (fOut.GetDirectory("Muon_Gen/Geant"))
        {
            fOut.cd("Muon_Gen/Geant");
            if (h_Nperevent_Muon_Gen_GeantOnly[i_Muon_Origin]->GetMean() > 0.0)
            {
                h_VzmotherEta_Muon_Gen_Geant[i_Muon_Origin]->Write(0, 2, 0);

                h_Nperevent_Muon_Gen_GeantOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtY_Muon_Gen_GeantOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtEta_Muon_Gen_GeantOnly[i_Muon_Origin]->Write(0, 2, 0);

                h_VzmotherEta_Muon_Gen_PYTHIA[i_Muon_Origin]->Write(0, 2, 0);

                h_Nperevent_Muon_Gen_DQcut_GeantOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtY_Muon_Gen_DQcut_GeantOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtEta_Muon_Gen_DQcut_GeantOnly[i_Muon_Origin]->Write(0, 2, 0);
            }
        }

        if (fOut.GetDirectory("Muon_Gen/PYTHIA"))
        {
            fOut.cd("Muon_Gen/PYTHIA");
            if (h_Nperevent_Muon_Gen_PYTHIAOnly[i_Muon_Origin]->GetMean() > 0.0)
            {
                h_Nperevent_Muon_Gen_PYTHIAOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtY_Muon_Gen_PYTHIAOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtEta_Muon_Gen_PYTHIAOnly[i_Muon_Origin]->Write(0, 2, 0);

                h_Nperevent_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtY_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtEta_Muon_Gen_DQcut_PYTHIAOnly[i_Muon_Origin]->Write(0, 2, 0);
            }
        }

        if (fOut.GetDirectory("Muon_Gen/Powheg"))
        {
            fOut.cd("Muon_Gen/Powheg");
            if (h_Nperevent_Muon_Gen_PowhegOnly[i_Muon_Origin]->GetMean() > 0.0)
            {
                h_Nperevent_Muon_Gen_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtY_Muon_Gen_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtEta_Muon_Gen_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);

                h_Nperevent_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtY_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtEta_Muon_Gen_DQcut_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
            }
        }
        // Rec Muon saving
        fOut.cd();
        fOut.cd("Muon_Rec");
        if (h_PtY_Muon_Rec[i_Muon_Origin]->GetEntries() > 0.0)
        {
            h_PtY_Muon_Rec[i_Muon_Origin]->Write(0, 2, 0);
            h_PtEta_Muon_Rec[i_Muon_Origin]->Write(0, 2, 0);
            h_Nperevent_Muon_Rec[i_Muon_Origin]->Write(0, 2, 0);
        }

        // Rec Muon from Powheg saving
        if (fOut.GetDirectory("Muon_Rec/Powheg"))
        {
            fOut.cd("Muon_Rec/Powheg");

            if (h_PtY_Muon_Rec_PowhegOnly[i_Muon_Origin]->GetEntries() > 0.0)
            {
                h_PtY_Muon_Rec_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtEta_Muon_Rec_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_Nperevent_Muon_Rec_PowhegOnly[i_Muon_Origin]->Write(0, 2, 0);
            }
        }

        if (fOut.GetDirectory("Muon_Rec/Geant"))
        {
            fOut.cd("Muon_Rec/Geant");
            if (h_PtY_Muon_Rec_GeantOnly[i_Muon_Origin]->GetEntries() > 0.0)
            {
                h_VzmotherEta_Muon_Rec_Geant[i_Muon_Origin]->Write(0, 2, 0);
                h_PtY_Muon_Rec_GeantOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtEta_Muon_Rec_GeantOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_Nperevent_Muon_Rec_GeantOnly[i_Muon_Origin]->Write(0, 2, 0);
            }
        }

        if (fOut.GetDirectory("Muon_Rec/PYTHIA"))
        {
            fOut.cd("Muon_Rec/PYTHIA");
            if (h_PtY_Muon_Rec_PYTHIAOnly[i_Muon_Origin]->GetEntries() > 0.0)
            {
                h_PtY_Muon_Rec_PYTHIAOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_PtEta_Muon_Rec_PYTHIAOnly[i_Muon_Origin]->Write(0, 2, 0);
                h_VzmotherEta_Muon_Rec_PYTHIA[i_Muon_Origin]->Write(0, 2, 0);
                h_Nperevent_Muon_Rec_PYTHIAOnly[i_Muon_Origin]->Write(0, 2, 0);
            }
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
        if (h_Nperevent_DiMuon_Gen[i_DiMuon_origin]->GetEntries() > 0.0)
        {
            h_Nperevent_DiMuon_Gen[i_DiMuon_origin]->Write(0, 2, 0);
            h_PtM_DiMuon_Gen[i_DiMuon_origin]->Write(0, 2, 0);
            h_PtY_DiMuon_Gen[i_DiMuon_origin]->Write(0, 2, 0);
        }

        if (h_Nperevent_DiMuon_Gen_DQcut[i_DiMuon_origin]->GetMean() > 0.0)
        {
            h_Nperevent_DiMuon_Gen_DQcut[i_DiMuon_origin]->Write(0, 2, 0);
            h_PtM_DiMuon_Gen_DQcut[i_DiMuon_origin]->Write(0, 2, 0);
            h_PtY_DiMuon_Gen_DQcut[i_DiMuon_origin]->Write(0, 2, 0);
        }

        if (Generator.Contains("Powheg"))
        {
            fOut.cd("DiMuon_Gen/Powheg");
            if (h_Nperevent_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->GetEntries() > 0.0)
            {
                h_Nperevent_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
                h_PtM_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
                h_PtY_DiMuon_Gen_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
            }
            if (h_Nperevent_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->GetEntries() > 0.0)
            {
                h_Nperevent_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
                h_PtM_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
                h_PtY_DiMuon_Gen_DQcut_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
            }
        }

        fOut.cd();
        fOut.cd("DiMuon_Rec");
        if (h_PtM_DiMuon_Rec[i_DiMuon_origin]->GetEntries() > 0.0)
        {
            h_PtM_DiMuon_Rec[i_DiMuon_origin]->Write(0, 2, 0);
            h_PtY_DiMuon_Rec[i_DiMuon_origin]->Write(0, 2, 0);
            h_Nperevent_DiMuon_Rec[i_DiMuon_origin]->Write(0, 2, 0);
        }

        if (Generator.Contains("Powheg"))
        {
            fOut.cd("DiMuon_Rec/Powheg");
            if (h_Nperevent_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->GetEntries() > 0.0)
            {
                h_Nperevent_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
                h_PtM_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
                h_PtY_DiMuon_Rec_PowhegOnly[i_DiMuon_origin]->Write(0, 2, 0);
            }
        }
    }

    if (fOut.GetDirectory("DiMuon_Rec/LF_origin"))
    {
        fOut.cd("DiMuon_Rec/LF_origin");
        for (Int_t i_LF_Generator = 0; i_LF_Generator < n_LF_DiMuon_Generator; i_LF_Generator++)
        {
            if (h_PtM_DiMuon_Rec_fromLF[i_LF_Generator]->GetEntries() > 0.0)
            {
                h_PtM_DiMuon_Rec_fromLF[i_LF_Generator]->Write(0, 2, 0);
                h_PtY_DiMuon_Rec_fromLF[i_LF_Generator]->Write(0, 2, 0);
            }
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
        if (!fOut.GetDirectory("Hadron"))
            fOut.mkdir("Hadron");
        fOut.cd("Hadron");

        if (Generator.Contains("Geant"))
        {
            if (!fOut.GetDirectory("Hadron/Geant"))
                fOut.mkdir("Hadron/Geant");

            if (!fOut.GetDirectory("Hadron/Pythia"))
                fOut.mkdir("Hadron/Pythia");
            fOut.cd("Hadron/Pythia");

            if (Generator.Contains("Powheg"))
            {
                if (!fOut.GetDirectory("Hadron/Pythia/Powheg"))
                    fOut.mkdir("Hadron/Pythia/Powheg");
            }

            fOut.cd("Hadron/Geant");
            h_PdgPt_HFHadron_prompt_GeantOnly->Write(0, 2, 0);
            h_PdgY_HFHadron_prompt_GeantOnly->Write(0, 2, 0);
            h_PdgEta_HFHadron_prompt_GeantOnly->Write(0, 2, 0);

            h_PdgPt_HFHadron_notprompt_GeantOnly->Write(0, 2, 0);
            h_PdgY_HFHadron_notprompt_GeantOnly->Write(0, 2, 0);
            h_PdgEta_HFHadron_notprompt_GeantOnly->Write(0, 2, 0);

            fOut.cd("Hadron/Pythia");
            h_PdgPt_HFHadron_prompt_PYTHIAOnly->Write(0, 2, 0);
            h_PdgY_HFHadron_prompt_PYTHIAOnly->Write(0, 2, 0);
            h_PdgEta_HFHadron_prompt_PYTHIAOnly->Write(0, 2, 0);

            h_PdgPt_HFHadron_notprompt_PYTHIAOnly->Write(0, 2, 0);
            h_PdgY_HFHadron_notprompt_PYTHIAOnly->Write(0, 2, 0);
            h_PdgEta_HFHadron_notprompt_PYTHIAOnly->Write(0, 2, 0);

            if (Generator.Contains("Powheg"))
            {
                fOut.cd("Hadron/Pythia/Powheg");
                h_PdgPt_HFHadron_prompt_PowhegOnly->Write(0, 2, 0);
                h_PdgY_HFHadron_prompt_PowhegOnly->Write(0, 2, 0);
                h_PdgEta_HFHadron_prompt_PowhegOnly->Write(0, 2, 0);

                h_PdgPt_HFHadron_notprompt_PowhegOnly->Write(0, 2, 0);
                h_PdgY_HFHadron_notprompt_PowhegOnly->Write(0, 2, 0);
                h_PdgEta_HFHadron_notprompt_PowhegOnly->Write(0, 2, 0);
            }
        }
        fOut.cd("Hadron");

        h_PdgPtY_HFHadron_prompt->Write(0, 2, 0);
        h_PdgPtY_HFHadron_notprompt->Write(0, 2, 0);

        h_PdgPt_HFHadron_prompt->Write(0, 2, 0);
        h_PdgY_HFHadron_prompt->Write(0, 2, 0);
        h_PdgEta_HFHadron_prompt->Write(0, 2, 0);

        h_PdgPt_HFHadron_notprompt->Write(0, 2, 0);
        h_PdgY_HFHadron_notprompt->Write(0, 2, 0);
        h_PdgEta_HFHadron_notprompt->Write(0, 2, 0);
    }

    fOut.Close();
    return;
}
