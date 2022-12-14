#include "mixedHF.h"

// Fallo con array di bool e poi così applica la condizione

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
    tree->SetBranchAddress("Pt_HFquark_gen", Pt_HFquark_gen);
    tree->SetBranchAddress("Y_HFquark_gen", Y_HFquark_gen);

    // tree->SetBranchAddress("N_HFquarks_rec", &N_HFquarks_rec);

    // tree->SetBranchAddress("PDG_HFquark_rec", PDG_HFquark_rec);
    // tree->SetBranchAddress("Pt_HFquark_rec", Pt_HFquark_rec);
    // tree->SetBranchAddress("Y_HFquark_rec", Y_HFquark_rec);

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

    if (from_data)
        printf("Data option not avaible\n");

    return tree;
}

void mixedHF(
    const char *RunMode = "HF",
    Int_t RunNumber = 294009,
    TString prefix_dir_filename = "30_09_2022",
    TString dir_fileout = "/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Mixed_Analysis/12_09_22",
    TString prefix_filename = "MCDimuHFTree")
{

    SetHist();
    TString dir_filename;
    // dir_filename.Form("/home/michele_pennisi/dimuon_HF_pp/Grid_Sim/%s/RunMCDimuHF/%s", RunMode, prefix_dir_filename.Data());//Local test
    dir_filename.Form("/media/michele_pennisi/DataBoi/Grid_Sim/%s/RunMCDimuHF/%s", RunMode, prefix_dir_filename.Data()); // official
    printf("%s", dir_filename.Data());

    TString filename;
    filename.Form("%s_MCDimuHFTree_%d.root", RunMode, RunNumber);

    TString fileout;
    // fileout.Form("%s/%s/%s_HFmixedAnalysis_%d.root", dir_fileout.Data(), RunMode, RunMode, RunNumber);
    fileout.Form("HFmixedAnalysis_%d.root", RunNumber);
    printf("Creating hist for %s/%s \n fileout %s \n", dir_filename.Data(), filename.Data(), fileout.Data());
    TChain *input_tree = Getting_Tree(dir_filename, filename, false);
    if (!input_tree)
        return;

    input_tree->ls();

    Bool_t **accepted_muon_ccut = new Bool_t *[n_Fraction];
    Bool_t **accepted_muon_bcut = new Bool_t *[n_Fraction];
    // vector<Int_t> accepted_muon[n_Fraction];

    for (Int_t i = 0; i < input_tree->GetEntries(); i++)
    {
        input_tree->GetEntry(i);
        if (i % 500000 == 0)
            printf("Evento :i %d\n", i);

        Int_t nMu_xevent_ccut[n_MuSelection][n_Mu_Charge][n_Fraction] = {0};
        Int_t nMu_xevent_bcut[n_MuSelection][n_Mu_Charge][n_Fraction] = {0};

        Int_t nDiMu_xevent_ccut[n_DiMuSelection][n_DiMu_Charge][n_Fraction] = {0};
        Int_t nDiMu_xevent_bcut[n_DiMuSelection][n_DiMu_Charge][n_Fraction] = {0};

        for (Int_t a = 0; a < NMuons_rec; a++)
        {
            Int_t PDG_Mu = PDGmum_rec[a]; // single rec mu PDG mum
            accepted_muon_ccut[a] = new Bool_t[n_Fraction];
            accepted_muon_bcut[a] = new Bool_t[n_Fraction];
            Double_t prob_Rand_rejection_charm = gRandom->Rndm();
            Double_t prob_Rand_rejection_beauty = gRandom->Rndm();

            // if (NMuons_rec > 1)
            // {
            //     if ((TMath::Abs(PDG_Mu) == 4) || (TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500) || (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000))
            //     {
            //         printf("c muon prob %0.3f", prob_Rand_rejection);
            //     }
            //     else if ((TMath::Abs(PDG_Mu) == 5) || (TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600) || (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000))
            //     {
            //         printf("b muon prob %0.3f", prob_Rand_rejection);
            //     }
            // }

            for (Int_t cut = 0; cut < n_Fraction; cut++)
            {
                accepted_muon_ccut[a][cut] = kFALSE;
                accepted_muon_bcut[a][cut] = kFALSE;
                if ((TMath::Abs(PDG_Mu) == 4) || (TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500) || (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000))
                {
                    if (prob_Rand_rejection_charm < Fraction[cut])
                    {
                        accepted_muon_ccut[a][cut] = kTRUE;
                    }
                }
                else if ((TMath::Abs(PDG_Mu) == 5) || (TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600) || (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000))
                {
                    if (prob_Rand_rejection_beauty < Fraction[cut])
                    {
                        accepted_muon_bcut[a][cut] = kTRUE;
                    }
                    

                }
            }
        }

        // Int_t frac = 5;

        for (Int_t q = 0; q < NMuons_rec; q++)
        {

            Int_t PDG_Mu1 = PDGmum_rec[q];             // single rec mu PDG mum
            Double_t Pt_Mu1 = Pt_rec[q];               // single rec mu pT
            Double_t E_Mu1 = E_rec[q];                 // single rec mu E
            Double_t Px_Mu1 = Px_rec[q];               // single rec mu px
            Double_t Py_Mu1 = Py_rec[q];               // single rec mu py
            Double_t Pz_Mu1 = Pz_rec[q];               // single rec mu pz
            Double_t Y_Mu1 = Y_rec[q];                 // single rec mu y
            Double_t Eta_Mu1 = Eta_rec[q];             // single rec mu eta
            Double_t Phi_Mu1 = Phi_rec[q];             // single rec mu phi
            Double_t Theta_Mu1 = Theta_rec[q];         // single rec mu theta
            Double_t MatchTrig_Mu1 = MatchTrig_rec[q]; // single gen mu theta
            Double_t RAbs_Mu1 = RAtAbsEnd_rec[q];
            Double_t pDCA_Mu1 = pDCA_rec[q];
            Int_t Charge_Mu1 = Charge_rec[q]; // single rec mu charge

            Bool_t From_HF1 = kFALSE;
            Bool_t From_Charm1 = kFALSE;
            Bool_t From_Beauty1 = kFALSE;
            Bool_t Good_Mu1 = kFALSE;

            Bool_t Selection_Mu_rec[n_MuSelection] = {kFALSE};
            Bool_t Charge_Mu_rec[n_Mu_Charge] = {kFALSE};

            TLorentzVector vector_mcp1;

            vector_mcp1.SetPxPyPzE(Px_Mu1, Py_Mu1, Pz_Mu1, E_Mu1);

            if ((TMath::Abs(PDG_Mu1) == 4) || (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500) || (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000))
            {
                Selection_Mu_rec[0] = kTRUE;
                Selection_Mu_rec[2] = kTRUE;
                From_HF1 = kTRUE;
                From_Charm1 = kTRUE;
            }
            else if ((TMath::Abs(PDG_Mu1) == 5) || (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600) || (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000))
            {
                Selection_Mu_rec[1] = kTRUE;
                Selection_Mu_rec[2] = kTRUE;
                From_HF1 = kTRUE;
                From_Beauty1 = kTRUE;
            }

            if ((RAbs_Mu1 > 17.6 && RAbs_Mu1 < 89.5) && pDCA_Mu1 == 1)
            {
                if ((Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5) && (MatchTrig_Mu1 >= 2) && (Pt_Mu1 > 0.0 && Pt_Mu1 < 30.0))
                {
                    Good_Mu1 = kTRUE;
                }
            }

            if ((Charge_Mu1) == -1)
            {
                Charge_Mu_rec[0] = kTRUE;
                Charge_Mu_rec[2] = kTRUE;
            }
            else if (Charge_Mu1 == 1)
            {
                // printf("Charge_Mu1 %d + Charge_Mu2 %d == %d || PDG1 %d PDG2 %d \n", Charge_Mu1, Charge_Mu2, Charge_Mu1 + Charge_Mu2, PDG_Mu1, PDG_Mu2);
                Charge_Mu_rec[1] = kTRUE;
                Charge_Mu_rec[2] = kTRUE;
            }

            for (Int_t Mu_a = 0; Mu_a < n_Mu_Charge; Mu_a++)
            {
                for (Int_t Mu_b = 0; Mu_b < n_MuSelection; Mu_b++)
                {
                    for (Int_t Mu_c = 0; Mu_c < n_Fraction; Mu_c++)
                    {
                        if ((Good_Mu1) && Charge_Mu_rec[Mu_a] && Selection_Mu_rec[Mu_b])
                        {
                            if (accepted_muon_ccut[q][Mu_c])
                            {
                                h_PtYMu_ccut[Mu_a][Mu_b][Mu_c]->Fill(Pt_Mu1, Y_Mu1);

                                h_pdgMu_ccut[Mu_a][Mu_b][Mu_c]->Fill(PDG_Mu1);

                                h_nMu_xevent_ccut[Mu_a][Mu_b][Mu_c]++;
                            }

                            if (accepted_muon_bcut[q][Mu_c])
                            {
                                h_PtYMu_bcut[Mu_a][Mu_b][Mu_c]->Fill(Pt_Mu1, Y_Mu1);

                                h_pdgMu_bcut[Mu_a][Mu_b][Mu_c]->Fill(PDG_Mu1);

                                h_nMu_xevent_bcut[Mu_a][Mu_b][Mu_c]++;
                            }
                        }
                    }
                }
            }

            for (Int_t j = q + 1; j < NMuons_rec; j++)
            {

                Int_t PDG_Mu2 = PDGmum_rec[j];             // single rec mu PDG mum
                Double_t Pt_Mu2 = Pt_rec[j];               // single rec mu pT
                Double_t E_Mu2 = E_rec[j];                 // single rec mu E
                Double_t Px_Mu2 = Px_rec[j];               // single rec mu px
                Double_t Py_Mu2 = Py_rec[j];               // single rec mu py
                Double_t Pz_Mu2 = Pz_rec[j];               // single rec mu pz
                Double_t Y_Mu2 = Y_rec[j];                 // single rec mu y
                Double_t Eta_Mu2 = Eta_rec[j];             // single rec mu eta
                Double_t Phi_Mu2 = Phi_rec[j];             // single rec mu phi
                Double_t Theta_Mu2 = Theta_rec[j];         // single rec mu theta
                Double_t MatchTrig_Mu2 = MatchTrig_rec[j]; // single gen mu theta
                Double_t RAbs_Mu2 = RAtAbsEnd_rec[j];
                Double_t pDCA_Mu2 = pDCA_rec[j];
                Int_t Charge_Mu2 = Charge_rec[j]; // single rec mu charge

                Bool_t From_HF2 = kFALSE;
                Bool_t From_Charm2 = kFALSE;
                Bool_t From_Beauty2 = kFALSE;
                Bool_t Good_Mu2 = kFALSE;

                Bool_t Good_Dimuon = kFALSE;
                Bool_t Selection_DiMu_rec[n_DiMuSelection] = {kFALSE};
                Bool_t Charge_DiMu_rec[n_DiMu_Charge] = {kFALSE};

                // printf("PDG %d Charge %d \n",PDG_Mu2,);

                if ((TMath::Abs(PDG_Mu2) == 4) || (TMath::Abs(PDG_Mu2) > 400 && TMath::Abs(PDG_Mu2) < 500) || (TMath::Abs(PDG_Mu2) > 4000 && TMath::Abs(PDG_Mu2) < 5000))
                {
                    From_HF2 = kTRUE;
                    From_Charm2 = kTRUE;
                }
                else if ((TMath::Abs(PDG_Mu2) == 5) || (TMath::Abs(PDG_Mu2) > 500 && TMath::Abs(PDG_Mu2) < 600) || (TMath::Abs(PDG_Mu2) > 5000 && TMath::Abs(PDG_Mu2) < 6000))
                {
                    From_HF2 = kTRUE;
                    From_Beauty2 = kTRUE;
                }

                if ((RAbs_Mu2 > 17.6 && RAbs_Mu2 < 89.5) && pDCA_Mu2 == 1)
                {
                    if ((Eta_Mu2 > -4.0 && Eta_Mu2 < -2.5) && (MatchTrig_Mu2 >= 2) && (Pt_Mu2 > 0.0 && Pt_Mu2 < 30.0))
                    {
                        Good_Mu2 = kTRUE;
                    }
                } // End if to apply the kinematic cut applied in the data

                if (From_HF1 && From_HF2)
                {
                    // Selection_DiMu_rec[3] = kTRUE;
                    if (From_Charm1 && From_Charm2)
                    {
                        Selection_DiMu_rec[0] = kTRUE;
                    }
                    else if (From_Beauty1 && From_Beauty2)
                    {
                        // printf("PDG1 %d PDG2 %d \n", PDG_Mu1, PDG_Mu2);
                        Selection_DiMu_rec[1] = kTRUE;
                    }
                    else if ((From_Charm1 && From_Beauty2) || (From_Charm2 && From_Beauty1))
                    {
                        Selection_DiMu_rec[2] = kTRUE;
                    }
                }

                if ((Charge_Mu1 + Charge_Mu2) == 0)
                {
                    Charge_DiMu_rec[0] = kTRUE;
                }
                else if (Charge_Mu1 + Charge_Mu2 == 2)
                {
                    // printf("Charge_Mu1 %d + Charge_Mu2 %d == %d || PDG1 %d PDG2 %d \n", Charge_Mu1, Charge_Mu2, Charge_Mu1 + Charge_Mu2, PDG_Mu1, PDG_Mu2);
                    Charge_DiMu_rec[3] = kTRUE;
                    Charge_DiMu_rec[1] = kTRUE;
                }
                else if (Charge_Mu1 + Charge_Mu2 == -2)
                {
                    Charge_DiMu_rec[3] = kTRUE;
                    Charge_DiMu_rec[2] = kTRUE;
                }

                TLorentzVector vector_mcp2;
                vector_mcp2.SetPxPyPzE(Px_Mu2, Py_Mu2, Pz_Mu2, E_Mu2);

                TLorentzVector dimu = vector_mcp1 + vector_mcp2;
                Double_t Y_DiMu = dimu.Rapidity();
                if (Y_DiMu > -4.0 && Y_DiMu < -2.5 && Good_Mu1 && Good_Mu2)
                {
                    Good_Dimuon = kTRUE;
                }

                for (Int_t DiMu_a = 0; DiMu_a < n_DiMu_Charge; DiMu_a++)
                {
                    for (Int_t DiMu_b = 0; DiMu_b < n_DiMuSelection; DiMu_b++)
                    {
                        if ((dimu.M() > 4 && dimu.M() <= 30.0) && (Good_Dimuon) && Charge_DiMu_rec[DiMu_a] && Selection_DiMu_rec[DiMu_b])
                        {
                            for (Int_t DiMu_c = 0; DiMu_c < n_Fraction; DiMu_c++)
                            {
                                if (DiMu_b == 2)
                                {
                                    if (accepted_muon_ccut[q][DiMu_c] || accepted_muon_ccut[j][DiMu_c])
                                    {
                                        h_PtYDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.Pt(), dimu.Rapidity());
                                        h_PtMDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.Pt(), dimu.M());
                                        h_MDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.M());
                                        h_pdgDimuMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(PDG_Mu1);
                                        h_pdgDimuMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(PDG_Mu2);
                                        nDiMu_xevent_ccut[DiMu_a][DiMu_b][DiMu_c]++;
                                    }
                                    if (accepted_muon_bcut[q][DiMu_c] || accepted_muon_bcut[j][DiMu_c])
                                    {
                                        h_PtYDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.Pt(), dimu.Rapidity());
                                        h_PtMDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.Pt(), dimu.M());
                                        h_MDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.M());
                                        h_pdgDimuMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(PDG_Mu1);
                                        h_pdgDimuMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(PDG_Mu2);
                                        nDiMu_xevent_bcut[DiMu_a][DiMu_b][DiMu_c]++;
                                    }
                                }
                                else
                                {
                                    if (accepted_muon_ccut[q][DiMu_c] && accepted_muon_ccut[j][DiMu_c])
                                    {
                                        h_PtYDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.Pt(), dimu.Rapidity());
                                        h_PtMDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.Pt(), dimu.M());
                                        h_MDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.M());
                                        h_pdgDimuMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(PDG_Mu1);
                                        h_pdgDimuMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Fill(PDG_Mu2);
                                        nDiMu_xevent_ccut[DiMu_a][DiMu_b][DiMu_c]++;
                                    }
                                    else if (accepted_muon_bcut[q][DiMu_c] && accepted_muon_bcut[j][DiMu_c])
                                    {
                                        h_PtYDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.Pt(), dimu.Rapidity());
                                        h_PtMDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.Pt(), dimu.M());
                                        h_MDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(dimu.M());
                                        h_pdgDimuMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(PDG_Mu1);
                                        h_pdgDimuMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Fill(PDG_Mu2);
                                        nDiMu_xevent_bcut[DiMu_a][DiMu_b][DiMu_c]++;
                                    }
                                }
                            }
                        }
                    } // End definition over n_DiMuSelection
                }     // End loop over n_DiMu_Charge
            }         // End Secondo loop over rec muon

        } // End loop on reconstructed muons
    }

    TFile *outFile = new TFile(fileout.Data(), "UPDATE");
    outFile->cd();
    outFile->ls();

    TDirectory *dir = outFile->GetDirectory("Muon");
    if (!dir)
        dir = outFile->mkdir("Muon");
    else
        printf("Muon already exists \n");

    for (Int_t Mu_a = 0; Mu_a < n_Mu_Charge; Mu_a++)
    {
        for (Int_t Mu_b = 0; Mu_b < n_MuSelection; Mu_b++)
        {
            dir = outFile->GetDirectory(Form("Muon/%s", name_Mu_Charge[Mu_a].Data()));
            if (!dir)
                dir = outFile->mkdir(Form("Muon/%s", name_Mu_Charge[Mu_a].Data()));
            // else
            //     printf("%s already exists \n", Form("Muon/%s", name_Mu_Charge[Mu_a].Data()));

            for (Int_t Mu_c = 0; Mu_c < n_Fraction; Mu_c++)
            {
                dir = outFile->GetDirectory(Form("Muon/%s/Fraction%0.2f", name_Mu_Charge[Mu_a].Data(), Fraction[Mu_c]));
                if (!dir)
                    dir = outFile->mkdir(Form("Muon/%s/Fraction%0.2f", name_Mu_Charge[Mu_a].Data(), Fraction[Mu_c]));
                // else
                //     printf("%s already exists \n", Form("Muon/%s/Fraction%0.2f", name_Mu_Charge[Mu_a].Data(), Fraction[Mu_c]));

                outFile->cd(Form("Muon/%s/Fraction%0.2f", name_Mu_Charge[Mu_a].Data(), Fraction[Mu_c]));

                h_PtYMu_ccut[Mu_a][Mu_b][Mu_c]->Write(0, 2, 0);
                h_pdgMu_ccut[Mu_a][Mu_b][Mu_c]->Write(0, 2, 0);

                h_PtYMu_bcut[Mu_a][Mu_b][Mu_c]->Write(0, 2, 0);
                h_pdgMu_bcut[Mu_a][Mu_b][Mu_c]->Write(0, 2, 0);
            }

        } // End definition over n_DiMuSelection

    } // End loop over n_DiMu_Charge

    outFile->cd();
    dir = outFile->GetDirectory("Dimuon");
    if (!dir)
        dir = outFile->mkdir("Dimuon");
    else
        printf("Dimuon already exists \n");

    for (Int_t DiMu_a = 0; DiMu_a < n_DiMu_Charge; DiMu_a++)
    {
        for (Int_t DiMu_b = 0; DiMu_b < n_DiMuSelection; DiMu_b++)
        {
            dir = outFile->GetDirectory(Form("Dimuon/%s", name_DiMu_Charge[DiMu_a].Data()));
            if (!dir)
                dir = outFile->mkdir(Form("Dimuon/%s", name_DiMu_Charge[DiMu_a].Data()));
            // else
            //     printf("%s already exists \n", Form("Dimuon/%s", name_DiMu_Charge[DiMu_a].Data()));

            for (Int_t DiMu_c = 0; DiMu_c < n_Fraction; DiMu_c++)
            {
                dir = outFile->GetDirectory(Form("Dimuon/%s/Fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), Fraction[DiMu_c]));
                if (!dir)
                    dir = outFile->mkdir(Form("Dimuon/%s/Fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), Fraction[DiMu_c]));
                // else
                //     printf("%s already exists \n", Form("Dimuon/%s/Fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), Fraction[DiMu_c]));

                outFile->cd(Form("Dimuon/%s/Fraction%0.2f", name_DiMu_Charge[DiMu_a].Data(), Fraction[DiMu_c]));

                h_PtYDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);
                h_PtMDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);
                h_MDiMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);
                h_pdgDimuMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);
                h_pdgDimuMu_ccut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);

                h_PtYDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);
                h_PtMDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);
                h_MDiMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);
                h_pdgDimuMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);
                h_pdgDimuMu_bcut[DiMu_a][DiMu_b][DiMu_c]->Write(0, 2, 0);
            }
        } // End definition over n_DiMuSelection

    } // End loop over n_DiMu_Charge
}