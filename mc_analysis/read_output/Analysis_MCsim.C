#include "Analysis_MCsim.h"

void Analysis_MCsim(
    TString RunMode = "HF",
    Int_t RunNumber = 294241,
    TString Task_Version = "Version1",
    Bool_t test = kTRUE,
    TString prefix_filename = "MCDimuHFTree")
{
    TString fileout;
    fileout.Form("%s_MCDimuHFTree_%d.root", RunMode.Data(), RunNumber);

    TString dir_fileout;
    if (test)
        dir_fileout = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/root_files/test";
    else
        dir_fileout.Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s/Analysis_MCsim", Task_Version.Data(), RunMode.Data()); // official with files saved locally

    printf("%s/%s\n", dir_fileout.Data(), fileout.Data());

    TString filename;
    filename.Form("%s_MCDimuHFTree_%d.root", RunMode.Data(), RunNumber);

    TString file_out;
    file_out.Form("%s_Analysis_MCsim_%d.root", RunMode.Data(), RunNumber);

    printf("Input File: %s\n", filename.Data());
    printf("Saving in dir: %s \nFile: %s\n", dir_fileout.Data(), file_out.Data());

    TChain *input_tree = Importing_Tree(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s", Task_Version.Data(), RunMode.Data()), filename);
    input_tree->ls();
    Int_t total_entries = input_tree->GetEntries();

    Set_Hist();
    for (Int_t i_Event = 0; i_Event < total_entries; i_Event++)
    {
        if (i_Event % 500000 == 0)
        {
            progress_status(i_Event, total_entries);
        }
        input_tree->GetEntry(i_Event);

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

            Bool_t HF_mu = kFALSE;

            Bool_t Charm_mu = kFALSE;
            Bool_t Charm_mu_Meson = kFALSE;
            Bool_t Charm_mu_Barion = kFALSE;

            Bool_t Beauty_mu = kFALSE;
            Bool_t Beauty_mu_Meson = kFALSE;
            Bool_t Beauty_mu_Barion = kFALSE;

            Bool_t DQ_Muon = kFALSE;
            if (Pt_Mu > 0.0 && Pt_Mu < 30.0)
                if (Eta_Mu > -4.0 && Eta_Mu < -2.5)
                    DQ_Muon = kTRUE;

            if (!DQ_Muon)
                continue;

            if ((TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500) || (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000))
            {
                HF_mu = kTRUE;
                Charm_mu = kTRUE;
                h_PtYPdg_Muon_Gen[0]->Fill(Pt_Mu, Y_Mu, 0);
                if (TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500)
                {
                    Charm_mu_Meson = kTRUE;
                    Charm_mu_Barion = kFALSE;
                    h_PtYPdg_Muon_Gen[0]->Fill(Pt_Mu, Y_Mu, 1);
                    h_PtPdg_Muon_Gen_Meson->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                }
                else if (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000)
                {
                    Charm_mu_Meson = kFALSE;
                    Charm_mu_Barion = kTRUE;
                    h_PtYPdg_Muon_Gen[0]->Fill(Pt_Mu, Y_Mu, 2);
                    h_PtPdg_Muon_Gen_Barion->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                }
            }
            else if ((TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600) || (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000))
            {
                HF_mu = kTRUE;
                Beauty_mu = kTRUE;
                h_PtYPdg_Muon_Gen[1]->Fill(Pt_Mu, Y_Mu, 0);
                if (TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600)
                {
                    Beauty_mu_Meson = kTRUE;
                    Beauty_mu_Barion = kFALSE;
                    h_PtYPdg_Muon_Gen[1]->Fill(Pt_Mu, Y_Mu, 1);
                    h_PtPdg_Muon_Gen_Meson->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                }
                else if (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000)
                {
                    Beauty_mu_Meson = kFALSE;
                    Beauty_mu_Barion = kTRUE;
                    h_PtYPdg_Muon_Gen[1]->Fill(Pt_Mu, Y_Mu, 2);
                    h_PtPdg_Muon_Gen_Barion->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                }
            }
            else if ((TMath::Abs(PDG_Mu) > 100 && TMath::Abs(PDG_Mu) < 400) || (TMath::Abs(PDG_Mu) > 1000 && TMath::Abs(PDG_Mu) < 4000))
            {
                HF_mu = kFALSE;
                // LF_mu=kTRUE;
                h_PtYPdg_Muon_Gen[2]->Fill(Pt_Mu, Y_Mu, 0);
                // if (TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600)
                // {
                //     Beauty_mu_Meson = kTRUE;
                //     Beauty_mu_Barion = kFALSE;
                // }
                // else if (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000)
                // {
                //     Beauty_mu_Meson = kFALSE;
                //     Beauty_mu_Barion = kTRUE;
                // }
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
            Double_t pDCA_Mu = pDCA_rec[i_NMuons_rec];

            Bool_t HF_mu = kFALSE;
            Bool_t LF_mu = kFALSE;

            Bool_t Charm_mu = kFALSE;
            Bool_t Charm_mu_Meson = kFALSE;
            Bool_t Charm_mu_Barion = kFALSE;

            Bool_t Beauty_mu = kFALSE;
            Bool_t Beauty_mu_Meson = kFALSE;
            Bool_t Beauty_mu_Barion = kFALSE;

            Bool_t DQ_Muon = kFALSE;
            if (pDCA_Mu == 1)
                if (RAbs_Mu > 17.6 && RAbs_Mu < 89.5)
                    if (Pt_Mu > 0.0 && Pt_Mu < 30.0)
                        if (Eta_Mu > -4.0 && Eta_Mu < -2.5)
                            DQ_Muon = kTRUE;

            if (!DQ_Muon)
                continue;

            // Bool_t Muon_origin_Selection[n_Muon_origin] = {kFALSE};

            if ((TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500) || (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000))
            {
                HF_mu = kTRUE;
                Charm_mu = kTRUE;
                h_PtYPdg_Muon_Rec[0]->Fill(Pt_Mu, Y_Mu, 0);
                if (TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500)
                {
                    Charm_mu_Meson = kTRUE;
                    Charm_mu_Barion = kFALSE;
                    h_PtYPdg_Muon_Rec[0]->Fill(Pt_Mu, Y_Mu, 1);
                    h_PtPdg_Muon_Rec_Meson->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                }
                else if (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000)
                {
                    Charm_mu_Meson = kFALSE;
                    Charm_mu_Barion = kTRUE;
                    h_PtYPdg_Muon_Rec[0]->Fill(Pt_Mu, Y_Mu, 2);
                    h_PtPdg_Muon_Rec_Barion->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                }
            }
            else if ((TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600) || (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000))
            {
                HF_mu = kTRUE;
                Beauty_mu = kTRUE;
                h_PtYPdg_Muon_Rec[1]->Fill(Pt_Mu, Y_Mu, 0);
                if (TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600)
                {
                    Beauty_mu_Meson = kTRUE;
                    Beauty_mu_Barion = kFALSE;
                    h_PtYPdg_Muon_Rec[1]->Fill(Pt_Mu, Y_Mu, 1);
                    h_PtPdg_Muon_Rec_Meson->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                }
                else if (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000)
                {
                    Beauty_mu_Meson = kFALSE;
                    Beauty_mu_Barion = kTRUE;
                    h_PtYPdg_Muon_Rec[1]->Fill(Pt_Mu, Y_Mu, 2);
                    h_PtPdg_Muon_Rec_Barion->Fill(Pt_Mu, TMath::Abs(PDG_Mu));
                }
            }
            else if ((TMath::Abs(PDG_Mu) > 100 && TMath::Abs(PDG_Mu) < 400) || (TMath::Abs(PDG_Mu) > 1000 && TMath::Abs(PDG_Mu) < 4000))
            {
                HF_mu = kFALSE;
                LF_mu = kTRUE;
                h_PtYPdg_Muon_Rec[2]->Fill(Pt_Mu, Y_Mu, 0);
                // if (TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600)
                // {
                //     Beauty_mu_Meson = kTRUE;
                //     Beauty_mu_Barion = kFALSE;
                // }
                // else if (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000)
                // {
                //     Beauty_mu_Meson = kFALSE;
                //     Beauty_mu_Barion = kTRUE;
                // }
            }
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

            Bool_t Charm_mu0 = kFALSE;
            Bool_t Charm_mu0_Meson = kFALSE;
            Bool_t Charm_mu0_Barion = kFALSE;

            Bool_t Charm_mu1 = kFALSE;
            Bool_t Charm_mu1_Meson = kFALSE;
            Bool_t Charm_mu1_Barion = kFALSE;

            Bool_t Beauty_mu0 = kFALSE;
            Bool_t Beauty_mu0_Meson = kFALSE;
            Bool_t Beauty_mu0_Barion = kFALSE;

            Bool_t Beauty_mu1 = kFALSE;
            Bool_t Beauty_mu1_Meson = kFALSE;
            Bool_t Beauty_mu1_Barion = kFALSE;

            Bool_t HF_mu0 = kFALSE;
            Bool_t HF_mu1 = kFALSE;
            Bool_t LF_mu0 = kFALSE;
            Bool_t LF_mu1 = kFALSE;

            Bool_t DQ_Dimuon = kFALSE;
            if ((Y_DiMu > -4.0 && Y_DiMu < -2.5) && (Eta_Mu0 > -4.0 && Eta_Mu0 < -2.5) && (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5))
                DQ_Dimuon = kTRUE;

            if (!DQ_Dimuon)
                continue;

            if ((TMath::Abs(PDG_Mu0) == 4) || (TMath::Abs(PDG_Mu0) > 400 && TMath::Abs(PDG_Mu0) < 500) || (TMath::Abs(PDG_Mu0) > 4000 && TMath::Abs(PDG_Mu0) < 5000))
            {
                HF_mu0 = kTRUE;
                Charm_mu0 = kTRUE;

                if (TMath::Abs(PDG_Mu0) > 400 && TMath::Abs(PDG_Mu0) < 500)
                {
                    Charm_mu0_Meson = kTRUE;
                    Charm_mu0_Barion = kFALSE;
                }
                else if (TMath::Abs(PDG_Mu0) > 4000 && TMath::Abs(PDG_Mu0) < 5000)
                {
                    Charm_mu0_Meson = kFALSE;
                    Charm_mu0_Barion = kTRUE;
                }
            }
            else if ((TMath::Abs(PDG_Mu0) == 5) || (TMath::Abs(PDG_Mu0) > 500 && TMath::Abs(PDG_Mu0) < 600) || (TMath::Abs(PDG_Mu0) > 5000 && TMath::Abs(PDG_Mu0) < 6000))
            {
                HF_mu0 = kTRUE;
                Beauty_mu0 = kTRUE;
                if (TMath::Abs(PDG_Mu0) > 500 && TMath::Abs(PDG_Mu0) < 600)
                {
                    Beauty_mu0_Meson = kTRUE;
                    Beauty_mu0_Barion = kFALSE;
                }
                else if (TMath::Abs(PDG_Mu0) > 5000 && TMath::Abs(PDG_Mu0) < 6000)
                {
                    Beauty_mu0_Meson = kFALSE;
                    Beauty_mu0_Barion = kTRUE;
                }
            }

            if ((TMath::Abs(PDG_Mu1) == 4) || (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500) || (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000))
            {
                HF_mu1 = kTRUE;
                Charm_mu1 = kTRUE;
                if (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500)
                {
                    Charm_mu1_Meson = kTRUE;
                    Charm_mu1_Barion = kFALSE;
                }
                else if (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000)
                {
                    Charm_mu1_Meson = kFALSE;
                    Charm_mu1_Barion = kTRUE;
                }
            }
            else if ((TMath::Abs(PDG_Mu1) == 5) || (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600) || (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000))
            {
                HF_mu1 = kTRUE;
                Beauty_mu1 = kTRUE;
                if (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600)
                {
                    Beauty_mu1_Meson = kTRUE;
                    Beauty_mu1_Barion = kFALSE;
                }
                else if (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000)
                {
                    Beauty_mu1_Meson = kFALSE;
                    Beauty_mu1_Barion = kTRUE;
                }
            }

            Bool_t DiMu_origin_Selection[n_DiMu_origin] = {kFALSE};
            Int_t type_DiMu = 999;
            if (Charm_mu0 && Charm_mu1)
            {
                DiMu_origin_Selection[0] = kTRUE;
                if (Charm_mu0_Meson && Charm_mu1_Meson)
                    type_DiMu = 0;
                if (Charm_mu0_Barion && Charm_mu1_Barion)
                    type_DiMu = 1;
                if ((Charm_mu0_Meson && Charm_mu1_Barion) || (Charm_mu1_Meson && Charm_mu0_Barion))
                    type_DiMu = 2;
            }
            else if (Beauty_mu0 && Beauty_mu1)
            {
                DiMu_origin_Selection[1] = kTRUE;
                if (Beauty_mu0_Meson && Beauty_mu1_Meson)
                    type_DiMu = 0;
                if (Beauty_mu0_Barion && Beauty_mu1_Barion)
                    type_DiMu = 1;
                if ((Beauty_mu0_Meson && Beauty_mu1_Barion) || (Beauty_mu1_Meson && Beauty_mu0_Barion))
                    type_DiMu = 2;
            }
            else if ((Charm_mu0 && Beauty_mu1) || (Charm_mu1 && Beauty_mu0))
            {
                DiMu_origin_Selection[2] = kTRUE;
                type_DiMu = 0;
            }
            else if (LF_mu0 && LF_mu1)
            {
                DiMu_origin_Selection[3] = kTRUE;
                type_DiMu = 0;
            }
            else if ((LF_mu0 && HF_mu1) || (LF_mu1 && HF_mu0))
            {
                DiMu_origin_Selection[4] = kTRUE;
                type_DiMu = 0;
            }

            for (Int_t i_DiMu_origin = 0; i_DiMu_origin < n_DiMu_origin; i_DiMu_origin++)
            {
                if (DiMu_origin_Selection[i_DiMu_origin])
                {
                    if (Charge_DiMu == 0)
                    {
                        h_PtMPdg_DiMu_Gen_ULS[i_DiMu_origin]->Fill(Pt_DiMu, M_DiMu, type_DiMu);
                        h_PtYPdg_DiMu_Gen_ULS[i_DiMu_origin]->Fill(Pt_DiMu, Y_DiMu, type_DiMu);
                        h_PtPt_DiMuULS_Muon_Gen[i_DiMu_origin]->Fill(Pt_DiMu, Pt_Mu0);
                        h_PtPt_DiMuULS_Muon_Gen[i_DiMu_origin]->Fill(Pt_DiMu, Pt_Mu1);
                        if ((M_DiMu > 4 && M_DiMu < 9) && Pt_DiMu < 10)
                        {
                            h_PtYPdg_DiMu_Gen_ULS_M49_Pt010[i_DiMu_origin]->Fill(Pt_DiMu, Y_DiMu, type_DiMu);
                            h_PtPt_DiMuULS_M49_Pt010_Muon_Gen[i_DiMu_origin]->Fill(Pt_DiMu, Pt_Mu0);
                            h_PtPt_DiMuULS_M49_Pt010_Muon_Gen[i_DiMu_origin]->Fill(Pt_DiMu, Pt_Mu1);
                        }
                        
                    }
                }
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

            Double_t Y_Mu1 = Y_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Pt_Mu1 = Pt_rec[DimuMu_rec[i_NDimu_rec][1]];
            Int_t PDG_Mu1 = PDGmum_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Charge_Mu1 = Charge_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t RAbs_Mu1 = RAtAbsEnd_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t pDCA_Mu1 = pDCA_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Eta_Mu1 = Eta_rec[DimuMu_rec[i_NDimu_rec][1]];
            Double_t Phi_Mu1 = Phi_rec[DimuMu_rec[i_NDimu_rec][1]];

            Bool_t DQ_Dimuon = kFALSE;
            if (pDCA_Mu0 == 1 && pDCA_Mu1 == 1)

                if ((RAbs_Mu0 > 17.6 && RAbs_Mu0 < 89.5) && (RAbs_Mu1 > 17.6 && RAbs_Mu1 < 89.5))

                    if ((Pt_Mu0 > 0.0 && Pt_Mu0 < 30.0) && (Pt_Mu1 > 0.0 && Pt_Mu1 < 30.0))

                        if ((Y_DiMu > -4.0 && Y_DiMu < -2.5) && (Eta_Mu0 > -4.0 && Eta_Mu0 < -2.5) && (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5))

                            if (DimuMatch_rec[i_NDimu_rec] == 2)

                                DQ_Dimuon = kTRUE;

            if (!DQ_Dimuon)
                continue;
            Bool_t Charm_mu0 = kFALSE;
            Bool_t Charm_mu0_Meson = kFALSE;
            Bool_t Charm_mu0_Barion = kFALSE;

            Bool_t Charm_mu1 = kFALSE;
            Bool_t Charm_mu1_Meson = kFALSE;
            Bool_t Charm_mu1_Barion = kFALSE;

            Bool_t Beauty_mu0 = kFALSE;
            Bool_t Beauty_mu0_Meson = kFALSE;
            Bool_t Beauty_mu0_Barion = kFALSE;

            Bool_t Beauty_mu1 = kFALSE;
            Bool_t Beauty_mu1_Meson = kFALSE;
            Bool_t Beauty_mu1_Barion = kFALSE;

            Bool_t HF_mu0 = kFALSE;
            Bool_t HF_mu1 = kFALSE;
            Bool_t LF_mu0 = kFALSE;
            Bool_t LF_mu1 = kFALSE;

            if ((TMath::Abs(PDG_Mu0) == 4) || (TMath::Abs(PDG_Mu0) > 400 && TMath::Abs(PDG_Mu0) < 500) || (TMath::Abs(PDG_Mu0) > 4000 && TMath::Abs(PDG_Mu0) < 5000))
            {
                HF_mu0 = kTRUE;
                Charm_mu0 = kTRUE;
                if (TMath::Abs(PDG_Mu0) > 400 && TMath::Abs(PDG_Mu0) < 500)
                {
                    Charm_mu0_Meson = kTRUE;
                    Charm_mu0_Barion = kFALSE;
                }
                else if (TMath::Abs(PDG_Mu0) > 4000 && TMath::Abs(PDG_Mu0) < 5000)
                {
                    Charm_mu0_Meson = kFALSE;
                    Charm_mu0_Barion = kTRUE;
                }
            }
            else if ((TMath::Abs(PDG_Mu0) == 5) || (TMath::Abs(PDG_Mu0) > 500 && TMath::Abs(PDG_Mu0) < 600) || (TMath::Abs(PDG_Mu0) > 5000 && TMath::Abs(PDG_Mu0) < 6000))
            {
                HF_mu0 = kTRUE;
                Beauty_mu0 = kTRUE;
                if (TMath::Abs(PDG_Mu0) > 500 && TMath::Abs(PDG_Mu0) < 600)
                {
                    Beauty_mu0_Meson = kTRUE;
                    Beauty_mu0_Barion = kFALSE;
                }
                else if (TMath::Abs(PDG_Mu0) > 5000 && TMath::Abs(PDG_Mu0) < 6000)
                {
                    Beauty_mu0_Meson = kFALSE;
                    Beauty_mu0_Barion = kTRUE;
                }
            }
            else if ((TMath::Abs(PDG_Mu0) > 0 && TMath::Abs(PDG_Mu0) < 4) || (TMath::Abs(PDG_Mu0) > 100 && TMath::Abs(PDG_Mu0) < 400) || (TMath::Abs(PDG_Mu0) > 1000 && TMath::Abs(PDG_Mu0) < 4000))
            {
                LF_mu0 = kTRUE;
                Charm_mu0 = kFALSE;
                Beauty_mu0 = kFALSE;
                HF_mu0 = kFALSE;
            }
            else
            {
                LF_mu0 = kFALSE;
                Charm_mu0 = kFALSE;
                Beauty_mu0 = kFALSE;
                HF_mu0 = kFALSE;
            }

            if ((TMath::Abs(PDG_Mu1) == 4) || (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500) || (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000))
            {
                HF_mu1 = kTRUE;
                Charm_mu1 = kTRUE;
                if (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500)
                {
                    Charm_mu1_Meson = kTRUE;
                    Charm_mu1_Barion = kFALSE;
                }
                else if (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000)
                {
                    Charm_mu1_Meson = kFALSE;
                    Charm_mu1_Barion = kTRUE;
                }
            }
            else if ((TMath::Abs(PDG_Mu1) == 5) || (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600) || (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000))
            {
                HF_mu1 = kTRUE;
                Beauty_mu1 = kTRUE;
                if (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600)
                {
                    Beauty_mu1_Meson = kTRUE;
                    Beauty_mu1_Barion = kFALSE;
                }
                else if (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000)
                {
                    Beauty_mu1_Meson = kFALSE;
                    Beauty_mu1_Barion = kTRUE;
                }
            }
            else if ((TMath::Abs(PDG_Mu1) > 0 && TMath::Abs(PDG_Mu1) < 4) || (TMath::Abs(PDG_Mu1) > 100 && TMath::Abs(PDG_Mu1) < 400) || (TMath::Abs(PDG_Mu1) > 1000 && TMath::Abs(PDG_Mu1) < 4000))
            {
                LF_mu1 = kTRUE;
                Charm_mu1 = kFALSE;
                Beauty_mu1 = kFALSE;
                HF_mu1 = kFALSE;
            }
            else
            {
                LF_mu1 = kFALSE;
                Charm_mu1 = kFALSE;
                Beauty_mu1 = kFALSE;
                HF_mu1 = kFALSE;
            }

            Bool_t DiMu_origin_Selection[n_DiMu_origin] = {kFALSE};
            Int_t type_DiMu = 999;
            if (Charm_mu0 && Charm_mu1)
            {
                DiMu_origin_Selection[0] = kTRUE;
                if (Charm_mu0_Meson && Charm_mu1_Meson)
                    type_DiMu = 0;
                if (Charm_mu0_Barion && Charm_mu1_Barion)
                    type_DiMu = 1;
                if ((Charm_mu0_Meson && Charm_mu1_Barion) || (Charm_mu1_Meson && Charm_mu0_Barion))
                    type_DiMu = 2;
            }
            else if (Beauty_mu0 && Beauty_mu1)
            {
                DiMu_origin_Selection[1] = kTRUE;
                if (Beauty_mu0_Meson && Beauty_mu1_Meson)
                    type_DiMu = 0;
                if (Beauty_mu0_Barion && Beauty_mu1_Barion)
                    type_DiMu = 1;
                if ((Beauty_mu0_Meson && Beauty_mu1_Barion) || (Beauty_mu1_Meson && Beauty_mu0_Barion))
                    type_DiMu = 2;
            }
            else if ((Charm_mu0 && Beauty_mu1) || (Charm_mu1 && Beauty_mu0))
            {
                DiMu_origin_Selection[2] = kTRUE;
                type_DiMu = 0;
            }
            else if (LF_mu0 && LF_mu1)
            {
                DiMu_origin_Selection[3] = kTRUE;
                type_DiMu = 0;
            }
            else if ((LF_mu0 && HF_mu1) || (LF_mu1 && HF_mu0))
            {
                DiMu_origin_Selection[4] = kTRUE;
                type_DiMu = 0;
            }

            for (Int_t i_DiMu_origin = 0; i_DiMu_origin < n_DiMu_origin; i_DiMu_origin++)
            {

                if (DiMu_origin_Selection[i_DiMu_origin])
                {
                    if (Charge_DiMu == 0)
                    {
                        h_PtMPdg_DiMu_Rec_ULS[i_DiMu_origin]->Fill(Pt_DiMu, M_DiMu, type_DiMu);
                        h_PtYPdg_DiMu_Rec_ULS[i_DiMu_origin]->Fill(Pt_DiMu, Y_DiMu, type_DiMu);
                        h_PtPt_DiMuULS_Muon_Rec[i_DiMu_origin]->Fill(Pt_DiMu, Pt_Mu0);
                        h_PtPt_DiMuULS_Muon_Rec[i_DiMu_origin]->Fill(Pt_DiMu, Pt_Mu1);
                        if ((M_DiMu > 4 && M_DiMu < 9) && Pt_DiMu < 10)
                        {

                            h_PtYPdg_DiMu_Rec_ULS_M49_Pt010[i_DiMu_origin]->Fill(Pt_DiMu, Y_DiMu, type_DiMu);
                            h_PtPt_DiMuULS_M49_Pt010_Muon_Rec[i_DiMu_origin]->Fill(Pt_DiMu, Pt_Mu0);
                            h_PtPt_DiMuULS_M49_Pt010_Muon_Rec[i_DiMu_origin]->Fill(Pt_DiMu, Pt_Mu1);
                        }
                    }
                }
            }
        }
    }
    TFile fOut(Form("%s/%s", dir_fileout.Data(), file_out.Data()), "RECREATE");
    fOut.cd();
    if (!fOut.GetDirectory("Muon_Gen"))
        fOut.mkdir("Muon_Gen");
    fOut.cd("Muon_Gen");
    h_PtPdg_Muon_Gen_Meson->Write();
    h_PtPdg_Muon_Gen_Barion->Write();
    if (!fOut.GetDirectory("Muon_Rec"))
        fOut.mkdir("Muon_Rec");

    fOut.cd("Muon_Rec");
    h_PtPdg_Muon_Rec_Meson->Write();
    h_PtPdg_Muon_Rec_Barion->Write();
    if (!fOut.GetDirectory("DiMu_Gen"))
        fOut.mkdir("DiMu_Gen");

    if (!fOut.GetDirectory("DiMu_Rec"))
        fOut.mkdir("DiMu_Rec");

    for (Int_t i_Muon_origin = 0; i_Muon_origin < n_Muon_origin; i_Muon_origin++)
    {
        fOut.cd("Muon_Gen");
        h_PtYPdg_Muon_Gen[i_Muon_origin]->Write();

        fOut.cd("Muon_Rec");
        h_PtYPdg_Muon_Rec[i_Muon_origin]->Write();
    }
    for (Int_t i_DiMu_origin = 0; i_DiMu_origin < n_DiMu_origin; i_DiMu_origin++)
    {
        fOut.cd("DiMu_Gen");
        h_PtMPdg_DiMu_Gen_ULS[i_DiMu_origin]->Write();
        h_PtYPdg_DiMu_Gen_ULS[i_DiMu_origin]->Write();
        h_PtYPdg_DiMu_Gen_ULS_M49_Pt010[i_DiMu_origin]->Write();
        h_PtPt_DiMuULS_Muon_Gen[i_DiMu_origin]->GetXaxis()->SetTitle("#it{p}_{T} dimu");
        h_PtPt_DiMuULS_Muon_Gen[i_DiMu_origin]->GetYaxis()->SetTitle("#it{p}_{T} mu");
        h_PtPt_DiMuULS_Muon_Gen[i_DiMu_origin]->Write();
        h_PtPt_DiMuULS_M49_Pt010_Muon_Gen[i_DiMu_origin]->GetXaxis()->SetTitle("#it{p}_{T} dimu");
        h_PtPt_DiMuULS_M49_Pt010_Muon_Gen[i_DiMu_origin]->GetYaxis()->SetTitle("#it{p}_{T} mu");
        h_PtPt_DiMuULS_M49_Pt010_Muon_Gen[i_DiMu_origin]->Write();

        fOut.cd("DiMu_Rec");
        h_PtMPdg_DiMu_Rec_ULS[i_DiMu_origin]->Write();
        h_PtYPdg_DiMu_Rec_ULS[i_DiMu_origin]->Write();
        h_PtYPdg_DiMu_Rec_ULS_M49_Pt010[i_DiMu_origin]->Write();
        h_PtPt_DiMuULS_Muon_Rec[i_DiMu_origin]->GetXaxis()->SetTitle("#it{p}_{T} dimu");
        h_PtPt_DiMuULS_Muon_Rec[i_DiMu_origin]->GetYaxis()->SetTitle("#it{p}_{T} mu");
        h_PtPt_DiMuULS_Muon_Rec[i_DiMu_origin]->Write();
        h_PtPt_DiMuULS_M49_Pt010_Muon_Rec[i_DiMu_origin]->GetXaxis()->SetTitle("#it{p}_{T} dimu");
        h_PtPt_DiMuULS_M49_Pt010_Muon_Rec[i_DiMu_origin]->GetYaxis()->SetTitle("#it{p}_{T} mu");
        h_PtPt_DiMuULS_M49_Pt010_Muon_Rec[i_DiMu_origin]->Write();
    }

    fOut.Close();

    return;
}