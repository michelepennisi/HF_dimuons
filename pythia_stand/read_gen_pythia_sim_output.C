#include "read_gen_pythia_sim_output.h"

TChain *Getting_Tree(TString dir_filename, TString filename, bool from_data)
{
    TChain *tree = nullptr;

    // if (!gSystem->OpenDirectory(dir_filename.Data())) {
    //     printf("Dir not found\n");
    //     return tree;
    // }

    tree = new TChain("MCTree");
    //    printf("%s \n",Form("%s/%s",dir_filename.Data(),filename.Data()));
    tree->AddFile(Form("%s/%s", dir_filename.Data(), filename.Data()));

    tree->SetBranchAddress("N_HFquarks_gen", &N_HFquarks_gen);

    tree->SetBranchAddress("PDG_HFquark_gen", PDG_HFquark_gen);
    tree->SetBranchAddress("PDG_HFquark_gen_daughter1", PDG_HFquark_gen_daughter1);
    tree->SetBranchAddress("PDG_HFquark_gen_daughter2", PDG_HFquark_gen_daughter2);
    tree->SetBranchAddress("Pt_HFquark_gen", Pt_HFquark_gen);
    tree->SetBranchAddress("Y_HFquark_gen", Y_HFquark_gen);

    // tree->SetBranchAddress("N_HFquarks_rec", &N_HFquarks_rec);

    // tree->SetBranchAddress("PDG_HFquark_rec", PDG_HFquark_rec);
    // tree->SetBranchAddress("Pt_HFquark_rec", Pt_HFquark_rec);
    // tree->SetBranchAddress("Y_HFquark_rec", Y_HFquark_rec);

    tree->SetBranchAddress("NMuons_gen", &NMuons_gen);
    tree->SetBranchAddress("Pt_mum_gen", Pt_mum_gen);
    tree->SetBranchAddress("Y_mum_gen", Y_mum_gen);
    tree->SetBranchAddress("PDGmum_gen", PDGmum_gen);
    tree->SetBranchAddress("Promptmum_gen", Promptmum_gen);
    tree->SetBranchAddress("PDG_gen", PDG_gen);
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

    tree->SetBranchAddress("NHadron_gen", &NHadron_gen);
    tree->SetBranchAddress("PDGHadron_gen", PDGHadron_gen);
    tree->SetBranchAddress("PromptHadron_gen", PromptHadron_gen);
    tree->SetBranchAddress("Hadron_Pt_gen", Hadron_Pt_gen);
    tree->SetBranchAddress("Hadron_E_gen", Hadron_E_gen);
    tree->SetBranchAddress("Hadron_Px_gen", Hadron_Px_gen);
    tree->SetBranchAddress("Hadron_Py_gen", Hadron_Py_gen);
    tree->SetBranchAddress("Hadron_Pz_gen", Hadron_Pz_gen);
    tree->SetBranchAddress("Hadron_Y_gen", Hadron_Y_gen);
    tree->SetBranchAddress("Hadron_Eta_gen", Hadron_Eta_gen);
    tree->SetBranchAddress("Hadron_Phi_gen", Hadron_Phi_gen);
    tree->SetBranchAddress("Hadron_Theta_gen", Hadron_Theta_gen);

    if (from_data)
        printf("Data option not avaible\n");

    return tree;
}

void read_gen_pythia_sim_output()
{
    SetHist();
    TString dir_filename("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/sim");
    // TString filename("pythia_sim_SoftQCD_Def_100000_3216_DefaultBR.root");
    // TString filename("pythia_sim_SoftQCD_Def_100000_9170_DefaultBR.root");
    // TString filename("monash_100Mev.root");
    // TString filename("pythia_sim_SoftQCD_Mode2_10000_3216_DefaultBR.root");
    TString filename("mode2_sim_fixed_10Mev.root");

    printf("Creating hist for %s/%s\n", dir_filename.Data(), filename.Data());
    TString fileout(Form("test/Hist_%s", filename.Data()));

    TChain *input_tree = Getting_Tree(dir_filename, filename, false);

    if (!input_tree)
        return;

    input_tree->ls();
    Int_t tot_dimuon = 0;
    Int_t tot_cquark = 0;
    Int_t tot_cquark_2 = 0;
    Int_t PDG_hadron[n_Hadron_studied] = {421, 411, 431, 4122};
    Double_t Y_hadron_studied[n_Rapidity_studied] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5};
    for (Int_t w = 0; w < input_tree->GetEntries(); w++)
    // for (Int_t w = 0; w < 1000; w++)
    {
        if (w % 1000000 == 0)
            printf("Evento :i %d\n", w);
        input_tree->GetEntry(w);
        h_Nevents->Fill(1);
        tot_cquark += N_HFquarks_gen;
        Int_t N_charm_xevent = 0;
        Int_t N_anticharm_xevent = 0;
        Int_t N_beauty_xevent = 0;
        Int_t N_antibeauty_xevent = 0;
        Int_t N_charmpair_xevent = 0;
        Int_t N_beautypair_xevent = 0;

        // printf("N_HFquarks_gen %d\n", N_HFquarks_gen);
        Int_t N_charm = 0;
        Int_t N_anticharm = 0;

        Int_t N_beauty = 0;
        Int_t N_antibeauty = 0;
        for (Int_t i = 0; i < N_HFquarks_gen; i++)
        {
            tot_cquark_2++;

            Int_t PDG_quark_gen1 = PDG_HFquark_gen[i];  // single gen mu PDG mum
            Double_t Pt_quark_gen1 = Pt_HFquark_gen[i]; // single gen mu pT
            Double_t Y_quark_gen1 = Y_HFquark_gen[i];   // single gen mu E

            if (PDG_quark_gen1 == 4)
            {
                h_pt_charm_quark->Fill(Pt_quark_gen1);
                h_y_charm_quark->Fill(Y_quark_gen1);
                N_charm_xevent++;
                N_charm++;
            }
            if (PDG_quark_gen1 == -4)
            {
                h_pt_charm_quark->Fill(Pt_quark_gen1);
                h_y_charm_quark->Fill(Y_quark_gen1);
                N_charm_xevent++;
                N_anticharm++;
            }
            else if (PDG_quark_gen1 == 5)
            {
                h_pt_beauty_quark->Fill(Pt_quark_gen1);
                h_y_beauty_quark->Fill(Y_quark_gen1);
                N_beauty_xevent++;
                N_beauty++;
            }

            else if (PDG_quark_gen1 == -5)
            {
                h_pt_beauty_quark->Fill(Pt_quark_gen1);
                h_y_beauty_quark->Fill(Y_quark_gen1);
                N_beauty_xevent++;
                N_antibeauty++;
            }

            for (Int_t j = i + 1; j < N_HFquarks_gen; j++)
            {
                Int_t PDG_quark_gen2 = PDG_HFquark_gen[j];  // sjngle gen mu PDG mum
                Double_t Pt_quark_gen2 = Pt_HFquark_gen[j]; // single gen mu pT
                Double_t Y_quark_gen2 = Y_HFquark_gen[j];   // single gen mu E

                if (TMath::Abs(PDG_quark_gen1) == 4 && TMath::Abs(PDG_quark_gen2) == 4)
                {
                    if ((PDG_quark_gen1 + PDG_quark_gen2) == 0)
                    {
                        N_charmpair_xevent++;
                    }
                }

                if (TMath::Abs(PDG_quark_gen1) == 5 && TMath::Abs(PDG_quark_gen2) == 5)
                {
                    if ((PDG_quark_gen1 + PDG_quark_gen2) == 0)
                    {
                        N_beautypair_xevent++;
                    }
                }
            }
        } // End loop on quark generated
        h_Ncharm_pairs->Fill(N_charmpair_xevent);
        h_Nbeauty_pairs->Fill(N_beautypair_xevent);

        h_Ncharm_quark->Fill(N_charm_xevent);
        h_Nbeauty_quark->Fill(N_beauty_xevent);

        if (N_charm > 0 && N_anticharm > 0)
        {
            if (N_charm >= N_anticharm)
            {
                h_Ncharm_pairs_v4->Fill(N_anticharm);
            }
            if (N_charm < N_anticharm)
            {
                h_Ncharm_pairs_v4->Fill(N_charm);
            }
        }

        if (N_beauty > 0 && N_antibeauty > 0)
        {
            if (N_beauty >= N_antibeauty)
            {
                h_Nbeauty_pairs_v4->Fill(N_antibeauty);
            }
            if (N_beauty < N_antibeauty)
            {
                h_Nbeauty_pairs_v4->Fill(N_beauty);
            }
        }
        if (N_charm > 0 || N_anticharm > 0)
        {
            // printf("N_charm %d || N_anticharm %d\n", N_charm, N_anticharm);
            if (N_charm >= N_anticharm)
            {
                h_Ncharm_pairs_v5->Fill(N_charm);
            }
            else if (N_anticharm > N_charm)
            {
                h_Ncharm_pairs_v5->Fill(N_anticharm);
            }
        }
        if (N_beauty > 0 || N_antibeauty > 0)
        {
            if (N_beauty >= N_antibeauty)
            {
                h_Nbeauty_pairs_v5->Fill(N_beauty);
            }
            else if (N_antibeauty > N_beauty)
            {
                h_Nbeauty_pairs_v5->Fill(N_antibeauty);
            }
        }
        Double_t charm_mean = (N_charm + N_anticharm);

        Double_t beauty_mean = (N_beauty + N_antibeauty);

        h_Ncharm_pairs_v6->Fill(charm_mean);
        h_Nbeauty_pairs_v6->Fill(beauty_mean);

        Int_t nMu_xevent[n_MuSelection][n_Mu_Charge] = {0};
        // Loop on generated single muons
        for (Int_t i = 0; i < NMuons_gen; i++)
        {
            Int_t PDG_Mu = PDGmum_gen[i];       // single gen mu PDG mum
            Double_t Pt_Mu_Mum = Pt_mum_gen[i]; // single gen mu pT
            Double_t Y_Mu_Mum = Y_mum_gen[i];   // single gen mu pT
            Double_t Pt_Mu = Pt_gen[i];         // single gen mu pT
            Double_t E_Mu = E_gen[i];           // single gen mu E
            Double_t Px_Mu = Px_gen[i];         // single gen mu px
            Double_t Py_Mu = Py_gen[i];         // single gen mu py
            Double_t Pz_Mu = Pz_gen[i];         // single gen mu pz
            Double_t Y_Mu = Y_gen[i];           // single gen mu y
            Double_t Eta_Mu = Eta_gen[i];       // single gen mu eta
            Double_t Phi_Mu = Phi_gen[i];       // single gen mu phi
            Double_t Theta_Mu = Theta_gen[i];   // single gen mu theta
            Int_t Charge_Mu = Charge_gen[i];    // single gen mu theta

            Bool_t Selection_Mu_gen[n_MuSelection] = {kFALSE};
            Selection_Mu_gen[0] = kTRUE;

            Bool_t Charge_Mu_gen[n_Mu_Charge] = {kFALSE};
            Charge_Mu_gen[0] = kTRUE;
            // 0 For All, 1 For HF, 2 For Charm, 3 For Beauty, 4 For Charm Mesons, 5 For Charm Barions, 6 For Beauty Mesons, 7 For Beauty Barions

            if (TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500)
            {
                Selection_Mu_gen[1] = kTRUE;
                Selection_Mu_gen[2] = kTRUE;
                Selection_Mu_gen[4] = kTRUE;
            }
            else if (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000)
            {
                Selection_Mu_gen[1] = kTRUE;
                Selection_Mu_gen[2] = kTRUE;
                Selection_Mu_gen[5] = kTRUE;
            }
            else if (TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600)
            {
                Selection_Mu_gen[1] = kTRUE;
                Selection_Mu_gen[3] = kTRUE;
                Selection_Mu_gen[6] = kTRUE;
            }
            else if (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000)
            {
                Selection_Mu_gen[1] = kTRUE;
                Selection_Mu_gen[3] = kTRUE;
                Selection_Mu_gen[7] = kTRUE;
            }

            // printf("Charge_Mu %d\n", Charge_Mu);
            if (Charge_Mu == +1)
                Charge_Mu_gen[1] = kTRUE;
            else if (Charge_Mu == -1)
                Charge_Mu_gen[2] = kTRUE;

            for (Int_t Mu_a = 0; Mu_a < n_MuSelection; Mu_a++)
            {
                for (Int_t Mu_b = 0; Mu_b < n_Mu_Charge; Mu_b++)
                {
                    if ((Eta_Mu > -4.0 && Eta_Mu < -2.5) && Selection_Mu_gen[Mu_a] && Charge_Mu_gen[Mu_b])
                    {
                        h_PtYMu[Mu_a][Mu_b]->Fill(Pt_Mu, Y_Mu);
                        h_PtMu_PtMum[Mu_a][Mu_b]->Fill(Pt_Mu, Pt_Mu_Mum);
                        h_YMu_YMum[Mu_a][Mu_b]->Fill(Y_Mu, Y_Mu_Mum);
                        h_PhiMu[Mu_a][Mu_b]->Fill(Phi_Mu);
                        h_pdgMu[Mu_a][Mu_b]->Fill(PDG_Mu);
                        nMu_xevent[Mu_a][Mu_b]++;
                    }
                }

            } // End loop on muon selection

        } // End loop on generated muons

        // h_nMu_xevent[0][0]->Fill(nMu_xevent[0][0]);

        for (Int_t Mu_a = 0; Mu_a < n_MuSelection; Mu_a++)
        {
            for (Int_t Mu_b = 0; Mu_b < n_Mu_Charge; Mu_b++)
            {
                // printf("%d nMu_xevent[Mu_a][Mu_b]\n", nMu_xevent[Mu_a][Mu_b]);
                h_nMu_xevent[Mu_a][Mu_b]->Fill(nMu_xevent[Mu_a][Mu_b]);
            }

        } // End loop on muon selection

        Int_t nDiMu_xevent[n_DiMuSelection][n_DiMu_Charge] = {0};

        for (Int_t i_firstMu = 0; i_firstMu < NMuons_gen; i_firstMu++)
        {

            Int_t PDG_Mu1 = PDGmum_gen[i_firstMu];     // single gen mu PDG mum
            Double_t Pt_Mu1 = Pt_gen[i_firstMu];       // single gen mu pT
            Double_t E_Mu1 = E_gen[i_firstMu];         // single gen mu E
            Double_t Px_Mu1 = Px_gen[i_firstMu];       // single gen mu px
            Double_t Py_Mu1 = Py_gen[i_firstMu];       // single gen mu py
            Double_t Pz_Mu1 = Pz_gen[i_firstMu];       // single gen mu pz
            Double_t Y_Mu1 = Y_gen[i_firstMu];         // single gen mu y
            Double_t Eta_Mu1 = Eta_gen[i_firstMu];     // single gen mu eta
            Double_t Phi_Mu1 = Phi_gen[i_firstMu];     // single gen mu phi
            Double_t Theta_Mu1 = Theta_gen[i_firstMu]; // single gen mu theta
            Int_t Charge_Mu1 = Charge_gen[i_firstMu];  // single gen mu theta

            Bool_t Charm_Mu1 = kFALSE;
            Bool_t Charm_Mu1_Meson = kFALSE;
            Bool_t Charm_Mu1_Barion = kFALSE;

            Bool_t Beauty_Mu1 = kFALSE;
            Bool_t Beauty_Mu1_Meson = kFALSE;
            Bool_t Beauty_Mu1_Barion = kFALSE;

            Bool_t HF_Mu1 = kFALSE;

            if (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500)
            {
                HF_Mu1 = kTRUE;

                Charm_Mu1 = kTRUE;
                Charm_Mu1_Meson = kTRUE;
                Charm_Mu1_Barion = kFALSE;

                Beauty_Mu1 = kFALSE;
                Beauty_Mu1_Meson = kFALSE;
                Beauty_Mu1_Barion = kFALSE;
            }
            else if (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000)
            {
                HF_Mu1 = kTRUE;
                Charm_Mu1 = kTRUE;
                Charm_Mu1_Meson = kFALSE;
                Charm_Mu1_Barion = kTRUE;

                Beauty_Mu1 = kFALSE;
                Beauty_Mu1_Meson = kFALSE;
                Beauty_Mu1_Barion = kFALSE;
            }
            else if (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600)
            {
                HF_Mu1 = kTRUE;

                Beauty_Mu1 = kTRUE;
                Beauty_Mu1_Meson = kTRUE;
                Beauty_Mu1_Barion = kFALSE;

                Charm_Mu1 = kFALSE;
                Charm_Mu1_Meson = kFALSE;
                Charm_Mu1_Barion = kFALSE;
            }
            else if (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000)
            {
                HF_Mu1 = kTRUE;

                Beauty_Mu1 = kTRUE;
                Beauty_Mu1_Meson = kFALSE;
                Beauty_Mu1_Barion = kTRUE;

                Charm_Mu1 = kFALSE;
                Charm_Mu1_Meson = kFALSE;
                Charm_Mu1_Barion = kFALSE;
            }

            if (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5)
            {
                for (Int_t i_secMu = i_firstMu + 1; i_secMu < NMuons_gen; i_secMu++)
                {
                    Int_t PDG_Mu2 = PDGmum_gen[i_secMu];     // single gen mu PDG mum
                    Double_t Pt_Mu2 = Pt_gen[i_secMu];       // single gen mu pT
                    Double_t E_Mu2 = E_gen[i_secMu];         // single gen mu E
                    Double_t Px_Mu2 = Px_gen[i_secMu];       // single gen mu px
                    Double_t Py_Mu2 = Py_gen[i_secMu];       // single gen mu py
                    Double_t Pz_Mu2 = Pz_gen[i_secMu];       // single gen mu pz
                    Double_t Y_Mu2 = Y_gen[i_secMu];         // single gen mu y
                    Double_t Eta_Mu2 = Eta_gen[i_secMu];     // single gen mu eta
                    Double_t Phi_Mu2 = Phi_gen[i_secMu];     // single gen mu phi
                    Double_t Theta_Mu2 = Theta_gen[i_secMu]; // single gen mu theta
                    Int_t Charge_Mu2 = Charge_gen[i_secMu];  // single gen mu theta

                    Bool_t Charm_Mu2 = kFALSE;
                    Bool_t Charm_Mu2_Meson = kFALSE;
                    Bool_t Charm_Mu2_Barion = kFALSE;

                    Bool_t Beauty_Mu2 = kFALSE;
                    Bool_t Beauty_Mu2_Meson = kFALSE;
                    Bool_t Beauty_Mu2_Barion = kFALSE;

                    Bool_t HF_Mu2 = kFALSE;

                    if (TMath::Abs(PDG_Mu2) > 400 && TMath::Abs(PDG_Mu2) < 500)
                    {
                        HF_Mu2 = kTRUE;

                        Charm_Mu2 = kTRUE;
                        Charm_Mu2_Meson = kTRUE;
                        Charm_Mu2_Barion = kFALSE;

                        Beauty_Mu2 = kFALSE;
                        Beauty_Mu2_Meson = kFALSE;
                        Beauty_Mu2_Barion = kFALSE;
                    }
                    else if (TMath::Abs(PDG_Mu2) > 4000 && TMath::Abs(PDG_Mu2) < 5000)
                    {
                        HF_Mu2 = kTRUE;
                        Charm_Mu2 = kTRUE;
                        Charm_Mu2_Meson = kFALSE;
                        Charm_Mu2_Barion = kTRUE;

                        Beauty_Mu2 = kFALSE;
                        Beauty_Mu2_Meson = kFALSE;
                        Beauty_Mu2_Barion = kFALSE;
                    }
                    else if (TMath::Abs(PDG_Mu2) > 500 && TMath::Abs(PDG_Mu2) < 600)
                    {
                        HF_Mu2 = kTRUE;

                        Beauty_Mu2 = kTRUE;
                        Beauty_Mu2_Meson = kTRUE;
                        Beauty_Mu2_Barion = kFALSE;

                        Charm_Mu2 = kFALSE;
                        Charm_Mu2_Meson = kFALSE;
                        Charm_Mu2_Barion = kFALSE;
                    }
                    else if (TMath::Abs(PDG_Mu2) > 5000 && TMath::Abs(PDG_Mu2) < 6000)
                    {
                        HF_Mu2 = kTRUE;

                        Beauty_Mu2 = kTRUE;
                        Beauty_Mu2_Meson = kFALSE;
                        Beauty_Mu2_Barion = kTRUE;

                        Charm_Mu2 = kFALSE;
                        Charm_Mu2_Meson = kFALSE;
                        Charm_Mu2_Barion = kFALSE;
                    }

                    if (Eta_Mu2 > -4.0 && Eta_Mu2 < -2.5)
                    {
                        TLorentzVector Mu1(Px_Mu1, Py_Mu1, Pz_Mu1, E_Mu1);
                        TLorentzVector Mu2(Px_Mu2, Py_Mu2, Pz_Mu2, E_Mu2);
                        TLorentzVector DiMu = Mu1 + Mu2;

                        Double_t Pt_DiMu = DiMu.Pt();

                        Double_t fHadron_E_gen = DiMu.E();
                        Double_t Px_DiMu = DiMu.Px();
                        Double_t Py_DiMu = DiMu.Py();
                        Double_t Pz_DiMu = DiMu.Pz();
                        Double_t Y_DiMu = DiMu.Rapidity();
                        Double_t M_DiMu = DiMu.M();                  // gen dimuon invariant mass
                        Int_t Charge_DiMu = Charge_Mu1 + Charge_Mu2; // gen dimuon charge

                        if ((Y_DiMu > -4.0 && Y_DiMu < -2.5))
                        {
                            tot_dimuon++;
                            Bool_t Charge_DiMu_gen[n_DiMu_Charge] = {kFALSE};
                            //  0 => For All,
                            //  1 => For HF,
                            //  2 => For Charm
                            //  3 => For Beauty,
                            //  4 => For HF Mixed (one muon from Charm, one muon from Beauty)
                            //  5 => For Charm Mesons (two muons from Charm Mesons)
                            //  6 => For Charm Barions (two muon from Charm Barions)
                            //  7 => For Charm Mixed (one muon from Charm Mesons, one muon from Charm Barions)
                            //  8 => For Beauty Mesons (two muons from Beauty Mesons)
                            //  9 => For Beauty Barions (two muon from Beauty Barions)
                            // 10 => For Beauty Mixed (one muon from Beauty Mesons, one muon from Beauty Barions

                            Bool_t Selection_DiMu_gen[n_DiMuSelection] = {kFALSE};
                            Bool_t isGoodDimuon = kFALSE;

                            Selection_DiMu_gen[0] = kTRUE; // All selection

                            if (HF_Mu1 && HF_Mu2)
                            {
                                Selection_DiMu_gen[1] = kTRUE; // HF selection
                                if (Charm_Mu1 && Charm_Mu2)
                                    Selection_DiMu_gen[2] = kTRUE; // Charm selection
                                else if (Beauty_Mu1 && Beauty_Mu2)
                                    Selection_DiMu_gen[3] = kTRUE; // Beauty selection
                                else if ((Charm_Mu1 && Beauty_Mu2) || (Beauty_Mu1 && Charm_Mu2))
                                    Selection_DiMu_gen[4] = kTRUE; // Mixed selection

                                if (Charm_Mu1_Meson && Charm_Mu2_Meson)
                                    Selection_DiMu_gen[5] = kTRUE; // Charm mesons selection
                                else if (Charm_Mu1_Barion && Charm_Mu2_Barion)
                                    Selection_DiMu_gen[6] = kTRUE; // Charm Barion selection
                                else if ((Charm_Mu1_Meson && Charm_Mu2_Barion) || (Charm_Mu2_Meson && Charm_Mu1_Barion))
                                    Selection_DiMu_gen[7] = kTRUE; // Charm Mixed selection

                                if (Beauty_Mu1_Meson && Beauty_Mu2_Meson)
                                    Selection_DiMu_gen[8] = kTRUE; // Beauty mesons selection
                                else if (Beauty_Mu1_Barion && Beauty_Mu2_Barion)
                                    Selection_DiMu_gen[9] = kTRUE; // Beauty Barion selection
                                else if ((Beauty_Mu1_Meson && Beauty_Mu2_Barion) || (Beauty_Mu2_Meson && Beauty_Mu1_Barion))
                                    Selection_DiMu_gen[10] = kTRUE; // Beauty Mixed selection
                            }
                            // Gen dimuon DQ cuts

                            if (Charge_DiMu == 0)
                            {
                                Charge_DiMu_gen[0] = kTRUE;
                                Charge_DiMu_gen[1] = kFALSE;
                                Charge_DiMu_gen[2] = kFALSE;
                                Charge_DiMu_gen[3] = kFALSE;
                            }
                            else if (Charge_DiMu == 2)
                            {
                                Charge_DiMu_gen[0] = kFALSE;
                                Charge_DiMu_gen[1] = kTRUE;
                                Charge_DiMu_gen[2] = kFALSE;
                                Charge_DiMu_gen[3] = kTRUE;
                            }
                            else if (Charge_DiMu == -2)
                            {
                                Charge_DiMu_gen[0] = kFALSE;
                                Charge_DiMu_gen[1] = kFALSE;
                                Charge_DiMu_gen[2] = kTRUE;
                                Charge_DiMu_gen[3] = kTRUE;
                            }

                            for (Int_t DiMu_a = 0; DiMu_a < n_DiMuSelection; DiMu_a++)
                            {
                                for (Int_t DiMu_b = 0; DiMu_b < n_DiMu_Charge; DiMu_b++)
                                {
                                    if (Selection_DiMu_gen[DiMu_a] && Charge_DiMu_gen[DiMu_b])
                                    {
                                        // printf("Filling %s with Charge %d\n", name_DiMu_Charge[DiMu_a].Data(), Charge_DiMu);
                                        h_PtYDiMu[DiMu_a][DiMu_b]->Fill(Pt_DiMu, Y_DiMu);
                                        h_PtMDiMu[DiMu_a][DiMu_b]->Fill(Pt_DiMu, M_DiMu);
                                        h_MDiMu[DiMu_a][DiMu_b]->Fill(M_DiMu);
                                        h_pdgDimuMu[DiMu_a][DiMu_b]->Fill(PDG_Mu1);
                                        h_pdgDimuMu[DiMu_a][DiMu_b]->Fill(PDG_Mu2);
                                        nDiMu_xevent[DiMu_a][DiMu_b]++;
                                    }

                                } // End definition over DiMu selection

                            } // End definition over DiMu charge cut
                        }
                        // End selection over dimuon rapidity
                    } // End selection over second muon eta
                }     // Loop over second muon
            }         // End selection over first muon eta

        } // End loop on muon entries to form dimuons

        for (Int_t DiMu_a = 0; DiMu_a < n_DiMuSelection; DiMu_a++)
        {
            for (Int_t DiMu_b = 0; DiMu_b < n_DiMu_Charge; DiMu_b++)
            {
                // printf("%d nMu_xevent[Mu_a][Mu_b]\n", nMu_xevent[Mu_a][Mu_b]);
                h_nDiMu_xevent[DiMu_a][DiMu_b]->Fill(nMu_xevent[DiMu_a][DiMu_b]);
            }

        } // End loop on muon selection

        for (Int_t iHadron = 0; iHadron < NHadron_gen; iHadron++)
        {
            Int_t PDGHadron = PDGHadron_gen[iHadron];       // single gen hadron PDG
            Int_t PromptHadron = PromptHadron_gen[iHadron]; // hadron prompt

            Double_t Pt_hadron = Hadron_Pt_gen[iHadron]; // single gen hadron pT
            Double_t Y_hadron = Hadron_Y_gen[iHadron];   // single gen hadron E

            for (size_t i_hadron = 0; i_hadron < n_Hadron_studied; i_hadron++)
            {
                if (TMath::Abs(PDGHadron) == PDG_hadron[i_hadron])
                {
                    if (Y_hadron > -0.5 && Y_hadron < 0.5)
                    {
                        // if (TMath::Abs(PDGHadron) == 421)
                        //     // Dzero_MidY++;
                        if (TMath::Abs(PromptHadron) == 4)
                        {
                            h_ptHadron_prompt[i_hadron][0]->Fill(Pt_hadron);
                            h_ptHadron_prompt_forCScalc[i_hadron][0]->Fill(Pt_hadron);
                            h_yHadron_prompt[i_hadron][0]->Fill(Y_hadron);
                            h_pdgHadron_prompt[i_hadron][0]->Fill(PDGHadron);
                        }
                        else if (TMath::Abs(PromptHadron) == 5)
                        {
                            h_ptHadron_Notprompt[i_hadron][0]->Fill(Pt_hadron);
                            h_yHadron_Notprompt[i_hadron][0]->Fill(Y_hadron);
                            h_pdgHadron_Notprompt[i_hadron][0]->Fill(PDGHadron);
                        }
                    }
                }

                for (size_t i_rapidity = 1; i_rapidity < n_Rapidity_studied; i_rapidity++)
                {
                    if (Y_hadron > Y_hadron_studied[i_rapidity - 1] && Y_hadron < Y_hadron_studied[i_rapidity])
                    {
                        // if (TMath::Abs(PDGHadron) == 421)
                        //     // Dzero_FwdY++;
                        if (TMath::Abs(PDGHadron) == PDG_hadron[i_hadron])
                        {
                            if (TMath::Abs(PromptHadron) == 4)
                            {
                                h_ptHadron_prompt[i_hadron][i_rapidity]->Fill(Pt_hadron);
                                h_ptHadron_prompt_forCScalc[i_hadron][i_rapidity]->Fill(Pt_hadron);
                                h_yHadron_prompt[i_hadron][i_rapidity]->Fill(Y_hadron);
                                h_pdgHadron_prompt[i_hadron][i_rapidity]->Fill(PDGHadron);
                            }
                            else if (TMath::Abs(PromptHadron) == 5)
                            {
                                h_ptHadron_Notprompt[i_hadron][i_rapidity]->Fill(Pt_hadron);
                                h_yHadron_Notprompt[i_hadron][i_rapidity]->Fill(Y_hadron);
                                h_pdgHadron_Notprompt[i_hadron][i_rapidity]->Fill(PDGHadron);
                            }
                        }
                    }
                }
            }
        }
    }

    printf("totcquark %d", tot_cquark);
    printf("tot_cquark_2 %d", tot_cquark_2);
    printf("tot_dimuon %d", tot_dimuon);

    TFile *fOut = new TFile(Form("%s", fileout.Data()), "RECREATE");
    fOut->cd();
    h_Nevents->Write(0, 2, 0);
    TDirectory *dir_fOut = fOut->GetDirectory("CS");
    if (!dir_fOut)
        dir_fOut = fOut->mkdir("CS");
    else
        printf("%s already exists \n", "CS");
    fOut->cd("CS");

    h_Ncharm_quark->Write(0, 2, 0);
    h_pt_charm_quark->Write(0, 2, 0);
    h_y_charm_quark->Write(0, 2, 0);
    h_Nbeauty_quark->Write(0, 2, 0);
    h_pt_beauty_quark->Write(0, 2, 0);
    h_y_beauty_quark->Write(0, 2, 0);

    h_Ncharm_pairs->Write(0, 2, 0);
    h_Nbeauty_pairs->Write(0, 2, 0);

    h_Ncharm_pairs_v4->Write(0, 2, 0);
    h_Nbeauty_pairs_v4->Write(0, 2, 0);

    h_Ncharm_pairs_v5->Write(0, 2, 0);
    h_Nbeauty_pairs_v5->Write(0, 2, 0);

    h_Ncharm_pairs_v6->Write(0, 2, 0);
    h_Nbeauty_pairs_v6->Write(0, 2, 0);

    dir_fOut = fOut->GetDirectory("muon");
    if (!dir_fOut)
        dir_fOut = fOut->mkdir("muon");
    else
        printf("%s already exists \n", "muon");
    fOut->cd("muon");

    for (Int_t Mu_a = 0; Mu_a < n_MuSelection; Mu_a++)
    {
        for (Int_t Mu_b = 0; Mu_b < n_Mu_Charge; Mu_b++)
        {
            fOut->cd();
            dir_fOut = fOut->GetDirectory(Form("muon/%s", name_MuCharge[Mu_b].Data()));
            if (!dir_fOut)
                dir_fOut = fOut->mkdir(Form("muon/%s", name_MuCharge[Mu_b].Data()));
            else
                printf("%s already exists \n", Form("muon/%s", name_MuCharge[Mu_b].Data()));
            fOut->cd(Form("muon/%s", name_MuCharge[Mu_b].Data()));
            h_PtYMu[Mu_a][Mu_b]->Write(0, 2, 0);
            h_PtMu_PtMum[Mu_a][Mu_b]->Write(0, 2, 0);
            h_YMu_YMum[Mu_a][Mu_b]->Write(0, 2, 0);
            h_PhiMu[Mu_a][Mu_b]->Write(0, 2, 0);
            h_pdgMu[Mu_a][Mu_b]->Write(0, 2, 0);
            // nMu_xevent[Mu_a][Mu_b]++;
        }

    } // End loop on muon selection
    fOut->cd();
    dir_fOut = fOut->GetDirectory("dimuon");
    if (!dir_fOut)
        dir_fOut = fOut->mkdir("dimuon");
    else
        printf("%s already exists \n", "dimuon");
    fOut->cd("dimuon");
    for (Int_t DiMu_a = 0; DiMu_a < n_DiMuSelection; DiMu_a++)
    {
        for (Int_t DiMu_b = 0; DiMu_b < n_DiMu_Charge; DiMu_b++)
        {
            fOut->cd();
            dir_fOut = fOut->GetDirectory(Form("dimuon/%s", name_DiMu_Charge[DiMu_b].Data()));
            if (!dir_fOut)
                dir_fOut = fOut->mkdir(Form("dimuon/%s", name_DiMu_Charge[DiMu_b].Data()));
            else
                printf("%s already exists \n", Form("dimuon/%s", name_DiMu_Charge[DiMu_b].Data()));
            fOut->cd(Form("dimuon/%s", name_DiMu_Charge[DiMu_b].Data()));

            h_PtYDiMu[DiMu_a][DiMu_b]->Write(0, 2, 0);
            h_PtMDiMu[DiMu_a][DiMu_b]->Write(0, 2, 0);
            h_MDiMu[DiMu_a][DiMu_b]->Write(0, 2, 0);
            h_pdgDimuMu[DiMu_a][DiMu_b]->Write(0, 2, 0);
            h_pdgDimuMu[DiMu_a][DiMu_b]->Write(0, 2, 0);
            // nDiMu_xevent[DiMu_a][DiMu_b]++;

        } // End definition over DiMu selection

    } // End definition over DiMu charge cut
    fOut->cd();
    dir_fOut = fOut->GetDirectory("hf_hadron");
    if (!dir_fOut)
        dir_fOut = fOut->mkdir("hf_hadron");
    else
        printf("%s already exists \n", "hf_hadron");
    fOut->cd("hf_hadron");
    dir_fOut = fOut->GetDirectory("hf_hadron/prompt");
    if (!dir_fOut)
        dir_fOut = fOut->mkdir("hf_hadron/prompt");
    else
        printf("%s already exists \n", "hf_hadron/prompt");
    dir_fOut = fOut->GetDirectory("hf_hadron/not_prompt");
    if (!dir_fOut)
        dir_fOut = fOut->mkdir("hf_hadron/not_prompt");
    else
        printf("%s already exists \n", "hf_hadron/not_prompt");

    for (size_t i_hadron = 0; i_hadron < n_Hadron_studied; i_hadron++)
    {
        for (size_t i_rapidity = 0; i_rapidity < n_Rapidity_studied; i_rapidity++)
        {
            fOut->cd();
            fOut->cd("hf_hadron/prompt");
            h_ptHadron_prompt[i_hadron][i_rapidity]->Write(0, 2, 0);
            h_ptHadron_prompt_forCScalc[i_hadron][i_rapidity]->Write(0, 2, 0);
            h_yHadron_prompt[i_hadron][i_rapidity]->Write(0, 2, 0);
            h_pdgHadron_prompt[i_hadron][i_rapidity]->Write(0, 2, 0);

            fOut->cd();
            fOut->cd("hf_hadron/not_prompt");
            h_ptHadron_Notprompt[i_hadron][i_rapidity]->Write(0, 2, 0);
            h_yHadron_Notprompt[i_hadron][i_rapidity]->Write(0, 2, 0);
            h_pdgHadron_Notprompt[i_hadron][i_rapidity]->Write(0, 2, 0);
        }
    }

    fOut->Close();
}
