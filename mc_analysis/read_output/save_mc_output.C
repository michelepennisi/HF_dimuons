#include "save_mc_output.h"

void save_mc_output(
    // TString RunMode = "powheg_DY_mass_3_35",
    TString RunMode = "HF",
    Int_t RunNumber = 294009,
    TString Task_Version = "Version1",
    Bool_t test = kTRUE,
    TString prefix_filename = "MCDimuHFTree")
{

    Set_Histograms();

    TString dir_fileOut;
    TString dir_fileIn;
    if (RunMode.Contains("powheg"))
    {
        dir_fileIn.Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Powheg_Sim/%s/%s", RunMode.Data(), Task_Version.Data());                 // official with files saved locally
        dir_fileOut.Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Powheg_Sim/%s/%s/Analysis_MCsim", RunMode.Data(), Task_Version.Data()); // official with files saved locally
    }
    else
    {
        dir_fileIn.Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s", Task_Version.Data(), RunMode.Data());                 // official with files saved locally
        dir_fileOut.Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s/save_mc_output", Task_Version.Data(), RunMode.Data()); // official with files saved locally
    }

    if (test)
        dir_fileOut.Form("test");

    TString filename;
    filename.Form("%s_MCDimuHFTree_%d.root", RunMode.Data(), RunNumber);

    TString file_out;
    file_out.Form("%s_MC_output_Hist_%d.root", RunMode.Data(), RunNumber);

    TString file_out_tree;
    file_out_tree.Form("%s_MC_output_Tree_%d.root", RunMode.Data(), RunNumber);

    printf("%s/%s\n", dir_fileOut.Data(), file_out.Data());
    TTree *rec_tree_mucharm;
    TTree *rec_tree_muDY;

    Double_t *pt_dimu_DY = new Double_t;
    Double_t *m_dimu_DY = new Double_t;
    rec_tree_muDY = new TTree("rec_tree_muDY", "Tree pt dimuon from DY with 11<M<30");
    rec_tree_muDY->Branch("pt", pt_dimu_DY, "pt/D");
    rec_tree_muDY->Branch("m", m_dimu_DY, "m/D");
    if (RunMode.Contains("powheg"))
    {
    }
    else
    {
        rec_tree_mucharm = new TTree("rec_tree_mucharm", "Tree pt dimuon charm with 4<M<30");
    }

    printf("Input File: %s\n", filename.Data());
    printf("Saving in dir: %s \nFile: %s\n", dir_fileOut.Data(), file_out.Data());

    // TChain *input_tree = Importing_Tree(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s", Task_Version.Data(), RunMode.Data()), filename);
    TChain *input_tree = Importing_Tree(dir_fileIn, filename);
    input_tree->ls();

    Int_t total_entries = input_tree->GetEntries();

    if (total_entries == 0)
        return;

    Bool_t Verbose = kFALSE;

    for (Int_t i_Event = 0; i_Event < total_entries; i_Event++)
    {
        h_Nevents->Fill(1);
        if (i_Event % 250000 == 0)
        {
            progress_status(i_Event, total_entries);
        }
        input_tree->GetEntry(i_Event);
        Int_t n_Muon_Gen[n_Muon_origin] = {0};
        Int_t n_Muon_Gen_DQcut[n_Muon_origin] = {0};

        Int_t n_Muon_Rec[n_Muon_origin] = {0};

        Int_t n_DiMuon_Gen[n_DiMuon_origin] = {0};
        Int_t n_DiMuon_Gen_DQcut[n_DiMuon_origin] = {0};

        Int_t n_DiMuon_Rec[n_DiMuon_origin] = {0};

        for (Int_t i_NMuons_gen = 0; i_NMuons_gen < NMuons_gen; i_NMuons_gen++)
        {

            Int_t PDG_Mu = PDGmum_gen[i_NMuons_gen];       // single gen mu PDG mum
            Double_t Pt_Mu = Pt_gen[i_NMuons_gen];         // single gen mu pT
            Double_t E_Mu = E_gen[i_NMuons_gen];           // single gen mu E
            Double_t Px_Mu = Px_gen[i_NMuons_gen];         // single gen mu px
            Double_t Py_Mu = Py_gen[i_NMuons_gen];         // single gen mu py
            Double_t Pz_Mu = Pz_gen[i_NMuons_gen];         // single gen mu pz
            Double_t Y_Mu = Y_gen[i_NMuons_gen];           // single gen mu y
            Double_t Eta_Mu = Eta_gen[i_NMuons_gen];       // single gen mu eta
            Double_t Phi_Mu = Phi_gen[i_NMuons_gen];       // single gen mu phi
            Double_t Theta_Mu = Theta_gen[i_NMuons_gen];   // single gen mu theta
            Double_t Charge_Mu = Charge_gen[i_NMuons_gen]; // single gen mu theta

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

            for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
            {
                if (Selection_Muon[i_Muon_origin])
                {
                    h_PtY_Muon_Gen[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    n_Muon_Gen[i_Muon_origin]++;
                }
            }

            h_PtPdg_Muon_Gen->Fill(Pt_Mu, PDG_Mu);
            h_YPdg_Muon_Gen->Fill(Y_Mu, PDG_Mu);

            if (Eta_Mu > -4.0 && Eta_Mu < -2.5)
                DQ_Muon = kTRUE;
            if (!DQ_Muon)
                continue;

            h_PtPdg_Muon_Gen_DQcut->Fill(Pt_Mu, PDG_Mu);
            h_YPdg_Muon_Gen_DQcut->Fill(Y_Mu, PDG_Mu);

            for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
            {
                if (Selection_Muon[i_Muon_origin])
                {
                    h_PtY_Muon_Gen_DQcut[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    n_Muon_Gen_DQcut[i_Muon_origin]++;
                }
            }
        }

        for (Int_t i_Muon_Origin = 0; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
        {
            h_Nperevent_Muon_Gen[i_Muon_Origin]->Fill(n_Muon_Gen[i_Muon_Origin]);
            h_Nperevent_Muon_Gen_DQcut[i_Muon_Origin]->Fill(n_Muon_Gen_DQcut[i_Muon_Origin]);
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

            h_PtPdg_Muon_Rec->Fill(Pt_Mu, PDG_Mu);
            h_YPdg_Muon_Rec->Fill(Y_Mu, PDG_Mu);

            for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
            {
                if (Selection_Muon[i_Muon_origin])
                {
                    h_PtY_Muon_Rec[i_Muon_origin]->Fill(Pt_Mu, Y_Mu);
                    n_Muon_Rec[i_Muon_origin]++;
                }
            }
        }

        for (Int_t i_Muon_Origin = 0; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
            h_Nperevent_Muon_Rec[i_Muon_Origin]->Fill(n_Muon_Rec[i_Muon_Origin]);

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
                        BR_times_frag_PYTHIA_Mu0 = BR_charm_hadrons2mu_PYTHIA[i_charm_hadron] * Frag_charm_hadrons_PYTHIA[i_charm_hadron];
                        BR_times_frag_MEAS_Mu0 = BR_charm_hadrons2mu_MEAS[i_charm_hadron] * Frag_charm_hadrons_MEAS[i_charm_hadron];
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
                        BR_times_frag_PYTHIA_Mu1 = BR_charm_hadrons2mu_PYTHIA[i_charm_hadron] * Frag_charm_hadrons_PYTHIA[i_charm_hadron];
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
                DiMu_origin_Selection[5] = kTRUE;

            Bool_t DQ_Dimuon = kFALSE;
            if ((Y_DiMu > -4.0 && Y_DiMu < -2.5) && (Eta_Mu0 > -4.0 && Eta_Mu0 < -2.5) && (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5) && (Charge_DiMu == 0))
                DQ_Dimuon = kTRUE;

            Double_t charm_correction = (BR_times_frag_MEAS_Mu0 * BR_times_frag_MEAS_Mu1) / (BR_times_frag_PYTHIA_Mu0 * BR_times_frag_PYTHIA_Mu1);

            if (DQ_Dimuon)
            {
                h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Pt_DiMu);
                h_Pdg1Pdg2Y_DiMuon_Gen_DQcut->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Y_DiMu);
                h_Pdg1Pdg2M_DiMuon_Gen_DQcut->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), M_DiMu);
                std::cout<<"PDG_Mu0: "<<TMath::Abs(PDG_Mu0)<<" PDG_Mu1:"<<TMath::Abs(PDG_Mu1)<<endl;
                std::cout << "BR_times_frag_MEAS_Mu0: " << BR_times_frag_MEAS_Mu0 << " BR_times_frag_MEAS_Mu1: " << BR_times_frag_MEAS_Mu1 << "BR_times_frag_PYTHIA_Mu0: " << BR_times_frag_PYTHIA_Mu0 << " BR_times_frag_PYTHIA_Mu1" << BR_times_frag_PYTHIA_Mu0 << endl;
                h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_Charm_corrected->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Pt_DiMu, charm_correction);
                h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_Charm_corrected->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), Y_DiMu, charm_correction);
                h_Pdg1Pdg2M_DiMuon_Gen_DQcut_Charm_corrected->Fill(TMath::Abs(PDG_Mu0), TMath::Abs(PDG_Mu1), M_DiMu, charm_correction);
            }

            for (Int_t i_DiMuon_origin = 0; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
            {
                if (DiMu_origin_Selection[i_DiMuon_origin])
                {
                    // 0 for Charm, 1 for Beauty, 2 for HF Mixed, 3 for LF, 4 for LF-HF Mixed, 5 for DY
                    h_PtM_DiMuon_Gen[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                    h_PtY_DiMuon_Gen[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                    n_DiMuon_Gen[i_DiMuon_origin]++;
                    if (DQ_Dimuon)
                    {

                        h_PtM_DiMuon_Gen_DQcut[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                        h_PtY_DiMuon_Gen_DQcut[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                        n_DiMuon_Gen_DQcut[i_DiMuon_origin]++;
                        if (i_DiMuon_origin == 0)
                        {
                            h_PtM_DiMuon_Gen_DQcut_Charm_corrected->Fill(Pt_DiMu, M_DiMu);
                            h_PtY_DiMuon_Gen_DQcut_Charm_corrected->Fill(Pt_DiMu, Y_DiMu);
                        }
                    }
                }
            }
        }

        for (Int_t i_Dimuon_Origin = 0; i_Dimuon_Origin < n_DiMuon_origin; i_Dimuon_Origin++)
        {

            h_Nperevent_DiMuon_Gen[i_Dimuon_Origin]->Fill(n_DiMuon_Gen[i_Dimuon_Origin]);
            h_Nperevent_DiMuon_Gen_DQcut[i_Dimuon_Origin]->Fill(n_DiMuon_Gen_DQcut[i_Dimuon_Origin]);
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

            Double_t Y_Mu1 = Y_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Pt_Mu1 = Pt_rec[DimuMu_rec[i_NDimu_rec][1]];
            Int_t PDG_Mu1 = PDGmum_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Charge_Mu1 = Charge_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t RAbs_Mu1 = RAtAbsEnd_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t pDCA_Mu1 = pDCA_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Eta_Mu1 = Eta_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Phi_Mu1 = Phi_rec[DimuMu_rec[i_NDimu_rec][1]];

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
                DiMu_origin_Selection[5] = kTRUE;

            for (Int_t i_DiMuon_origin = 0; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
            {
                if (DiMu_origin_Selection[i_DiMuon_origin])
                {
                    h_PtM_DiMuon_Rec[i_DiMuon_origin]->Fill(Pt_DiMu, M_DiMu);
                    h_PtY_DiMuon_Rec[i_DiMuon_origin]->Fill(Pt_DiMu, Y_DiMu);
                    n_DiMuon_Rec[i_DiMuon_origin]++;
                }
            }
        }

        for (Int_t i_Dimuon_Origin = 0; i_Dimuon_Origin < n_DiMuon_origin; i_Dimuon_Origin++)
            h_Nperevent_DiMuon_Rec[i_Dimuon_Origin]->Fill(n_DiMuon_Rec[i_Dimuon_Origin]);
    }

    Int_t starting_Muon_origin = 999;
    if (RunMode.Contains("DY"))
        starting_Muon_origin = n_Muon_origin - 1;
    else
        starting_Muon_origin = 0;

    Int_t startingDiMuon_origin = 999;
    if (RunMode.Contains("DY"))
        startingDiMuon_origin = n_Muon_origin - 1;
    else
        startingDiMuon_origin = 0;

    TFile fOut(Form("%s/%s", dir_fileOut.Data(), file_out.Data()), "RECREATE");
    fOut.cd();
    h_Nevents->Write();
    if (!fOut.GetDirectory("Muon_Gen"))
        fOut.mkdir("Muon_Gen");
    fOut.cd("Muon_Gen");
    h_PtPdg_Muon_Gen->Write(0, 2, 0);
    h_YPdg_Muon_Gen->Write(0, 2, 0);

    h_PtPdg_Muon_Gen_DQcut->Write(0, 2, 0);
    h_YPdg_Muon_Gen_DQcut->Write(0, 2, 0);
    for (Int_t i_Muon_Origin = starting_Muon_origin; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
    {
        h_PtY_Muon_Gen_DQcut[i_Muon_Origin]->Write(0, 2, 0);
        h_Nperevent_Muon_Gen[i_Muon_Origin]->Write(0, 2, 0);
        h_Nperevent_Muon_Gen_DQcut[i_Muon_Origin]->Write(0, 2, 0);
    }
    if (!fOut.GetDirectory("Muon_Rec"))
        fOut.mkdir("Muon_Rec");
    fOut.cd("Muon_Rec");
    h_PtPdg_Muon_Rec->Write(0, 2, 0);
    h_YPdg_Muon_Rec->Write(0, 2, 0);
    for (Int_t i_Muon_Origin = starting_Muon_origin; i_Muon_Origin < n_Muon_origin; i_Muon_Origin++)
    {
        h_PtY_Muon_Rec[i_Muon_Origin]->Write(0, 2, 0);
        h_Nperevent_Muon_Rec[i_Muon_Origin]->Write(0, 2, 0);
    }

    if (!fOut.GetDirectory("DiMuon_Gen"))
        fOut.mkdir("DiMuon_Gen");
    fOut.cd("DiMuon_Gen");
    h_Pdg1Pdg2Pt_DiMuon_Gen->Write(0, 2, 0);
    h_Pdg1Pdg2Y_DiMuon_Gen->Write(0, 2, 0);
    h_Pdg1Pdg2M_DiMuon_Gen->Write(0, 2, 0);

    h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut->Write(0, 2, 0);
    h_Pdg1Pdg2Y_DiMuon_Gen_DQcut->Write(0, 2, 0);
    h_Pdg1Pdg2M_DiMuon_Gen_DQcut->Write(0, 2, 0);

    h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);
    h_Pdg1Pdg2Y_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);
    h_Pdg1Pdg2M_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);

    if (!fOut.GetDirectory("DiMuon_Rec"))
        fOut.mkdir("DiMuon_Rec");
    fOut.cd("DiMuon_Rec");
    h_Pdg1Pdg2Pt_DiMuon_Rec->Write(0, 2, 0);
    h_Pdg1Pdg2Y_DiMuon_Rec->Write(0, 2, 0);
    h_Pdg1Pdg2M_DiMuon_Rec->Write(0, 2, 0);

    fOut.cd();
    for (Int_t i_DiMuon_origin = startingDiMuon_origin; i_DiMuon_origin < n_DiMuon_origin; i_DiMuon_origin++)
    {
        fOut.cd("DiMuon_Gen");
        h_PtM_DiMuon_Gen[i_DiMuon_origin]->Write(0, 2, 0);
        h_PtY_DiMuon_Gen[i_DiMuon_origin]->Write(0, 2, 0);
        h_Nperevent_DiMuon_Gen[i_DiMuon_origin]->Write(0, 2, 0);

        h_PtM_DiMuon_Gen_DQcut[i_DiMuon_origin]->Write(0, 2, 0);
        h_PtY_DiMuon_Gen_DQcut[i_DiMuon_origin]->Write(0, 2, 0);
        h_Nperevent_DiMuon_Gen_DQcut[i_DiMuon_origin]->Write(0, 2, 0);

        if (i_DiMuon_origin == 0)
        {
            h_PtM_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);
            h_PtY_DiMuon_Gen_DQcut_Charm_corrected->Write(0, 2, 0);
        }

        fOut.cd();
        fOut.cd("DiMuon_Rec");
        h_PtM_DiMuon_Rec[i_DiMuon_origin]->Write(0, 2, 0);
        h_PtY_DiMuon_Rec[i_DiMuon_origin]->Write(0, 2, 0);
        h_Nperevent_DiMuon_Rec[i_DiMuon_origin]->Write(0, 2, 0);
    }

    fOut.Close();
    return;
}
