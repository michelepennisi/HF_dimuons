#include "Analysis_MCsim.h"

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
        dir_fileout.Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s", Task_Version.Data(), RunMode.Data()); // official with files saved locally

    printf("%s/%s\n", dir_fileout.Data(), fileout.Data());

    TString filename;
    filename.Form("%s_MCDimuHFTree_%d.root", RunMode.Data(), RunNumber);

    TString file_out;
    file_out.Form("%s_Analysis_MCsim_%d.root", RunMode.Data(), RunNumber);

    printf("Input File: %s\n", filename.Data());
    printf("Saving in dir: %s \nFile: %s\n", dir_fileout.Data(), file_out.Data());

    TChain *input_tree = Importing_Tree(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s", Task_Version.Data(), RunMode.Data()), filename);
    input_tree->ls();
    for (Int_t i_Event = 0; i_Event < input_tree->GetEntries(); i_Event++)
    {
        if (i_Event % 500000 == 0)
        {
            Double_t progress = (Double_t)i_Event / input_tree->GetEntries();
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
            std::cout << "] " << int(progress * 100.0) << " % (" << i_Event << "/ " << input_tree->GetEntries() << ")\r";
            std::cout.flush();

            std::cout << std::endl;
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

            Bool_t DQ_Muon=kFALSE;
            if (Pt_Mu > 0.0 && Pt_Mu < 30.0)
                if (Eta_Mu > -4.0 && Eta_Mu < -2.5)
                    DQ_Muon = kTRUE;

            if (!DQ_Muon)
                continue;

            if (TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 600)
            {
                h_PtYPdg_Muon_Gen_Meson->Fill(Pt_Mu, Y_Mu, PDG_Mu);
                h_PtPdg_Muon_Gen_Meson->Fill(Pt_Mu, PDG_Mu);
            }
            else if (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 6000)
            {
                h_PtYPdg_Muon_Gen_Barion->Fill(Pt_Mu, Y_Mu, PDG_Mu);
                h_PtPdg_Muon_Gen_Barion->Fill(Pt_Mu, PDG_Mu);
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

            if (TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 600)
            {
                h_PtYPdg_Muon_Rec_Meson->Fill(Pt_Mu, Y_Mu, PDG_Mu);
                h_PtPdg_Muon_Rec_Meson->Fill(Pt_Mu, PDG_Mu);
            }
            else if (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 6000)
            {
                h_PtYPdg_Muon_Rec_Barion->Fill(Pt_Mu, Y_Mu, PDG_Mu);
                h_PtPdg_Muon_Rec_Barion->Fill(Pt_Mu, PDG_Mu);
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
            Int_t DiMu_PDG = 999;

            if (Beauty_mu0_Meson && Beauty_mu1_Meson)
            {
                DiMu_PDG = 501;
            }
            else if (Charm_mu0_Meson && Charm_mu1_Meson)
            {
                DiMu_PDG = 401;
            }
            else if ((Charm_mu0_Meson && Beauty_mu1_Meson) || (Beauty_mu0_Meson && Charm_mu1_Meson))
            {
                DiMu_PDG = 601;
            }

            if (Beauty_mu0_Barion && Beauty_mu1_Barion)
            {
                DiMu_PDG = 5001;
            }
            else if (Charm_mu0_Barion && Charm_mu1_Barion)
            {
                DiMu_PDG = 4001;
            }
            else if ((Charm_mu0_Barion && Beauty_mu1_Barion) || (Beauty_mu0_Barion && Charm_mu1_Barion))
            {
                DiMu_PDG = 6001;
            }

            Bool_t DQ_Dimuon = kFALSE;
            if ((Y_DiMu > -4.0 && Y_DiMu < -2.5) && (Eta_Mu0 > -4.0 && Eta_Mu0 < -2.5) && (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5))
                DQ_Dimuon = kTRUE;

            if (!DQ_Dimuon)
                continue;

            if (Charge_DiMu == 0 && DiMu_PDG < 700)
            {
                h_PtMPdg_DiMu_Gen_Meson_ULS->Fill(Pt_DiMu, M_DiMu, DiMu_PDG);
                h_PtYPdg_DiMu_Gen_Meson_ULS->Fill(Pt_DiMu, Y_DiMu, DiMu_PDG);
            }

            else if (Charge_DiMu == 0 && (DiMu_PDG > 4000 && DiMu_PDG < 7000))
            {

                h_PtMPdg_DiMu_Gen_Barion_ULS->Fill(Pt_DiMu, M_DiMu, DiMu_PDG);
                h_PtYPdg_DiMu_Gen_Barion_ULS->Fill(Pt_DiMu, Y_DiMu, DiMu_PDG);
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

            Int_t DiMu_PDG = 999;

            if (Beauty_mu0_Meson && Beauty_mu1_Meson)
            {
                DiMu_PDG = 501;
            }
            else if (Charm_mu0_Meson && Charm_mu1_Meson)
            {
                DiMu_PDG = 401;
            }
            else if ((Charm_mu0_Meson && Beauty_mu1_Meson) || (Beauty_mu0_Meson && Charm_mu1_Meson))
            {
                DiMu_PDG = 601;
            }

            if (Beauty_mu0_Barion && Beauty_mu1_Barion)
            {
                DiMu_PDG = 5001;
            }
            else if (Charm_mu0_Barion && Charm_mu1_Barion)
            {
                DiMu_PDG = 4001;
            }
            else if ((Charm_mu0_Barion && Beauty_mu1_Barion) || (Beauty_mu0_Barion && Charm_mu1_Barion))
            {
                DiMu_PDG = 6001;
            }

            if (Charge_DiMu == 0 && DiMu_PDG < 700)
            {
                h_PtMPdg_DiMu_Rec_Meson_ULS->Fill(Pt_DiMu, M_DiMu, DiMu_PDG);
                h_PtYPdg_DiMu_Rec_Meson_ULS->Fill(Pt_DiMu, Y_DiMu, DiMu_PDG);
            }

            else if (Charge_DiMu == 0 && (DiMu_PDG > 4000 && DiMu_PDG < 7000))
            {

                h_PtMPdg_DiMu_Rec_Barion_ULS->Fill(Pt_DiMu, M_DiMu, DiMu_PDG);
                h_PtYPdg_DiMu_Rec_Barion_ULS->Fill(Pt_DiMu, Y_DiMu, DiMu_PDG);
            }
        }
    }
    TFile fOut(Form("%s/%s", dir_fileout.Data(), file_out.Data()), "RECREATE");
    fOut.cd();
    if (!fOut.GetDirectory("Muon_Gen"))
        fOut.mkdir("Muon_Gen");

    fOut.cd("Muon_Gen");
    h_PtYPdg_Muon_Gen_Meson->Write();
    h_PtYPdg_Muon_Gen_Barion->Write();

    h_PtPdg_Muon_Gen_Meson->Write();
    h_PtPdg_Muon_Gen_Barion->Write();

    if (!fOut.GetDirectory("Muon_Rec"))
        fOut.mkdir("Muon_Rec");

    fOut.cd("Muon_Rec");
    h_PtYPdg_Muon_Rec_Meson->Write();
    h_PtYPdg_Muon_Rec_Barion->Write();

    h_PtPdg_Muon_Rec_Meson->Write();
    h_PtPdg_Muon_Rec_Barion->Write();

    if (!fOut.GetDirectory("DiMu_Gen"))
        fOut.mkdir("DiMu_Gen");

    fOut.cd("DiMu_Gen");
    h_PtMPdg_DiMu_Gen_Meson_ULS->Write();
    h_PtMPdg_DiMu_Gen_Barion_ULS->Write();

    h_PtYPdg_DiMu_Gen_Meson_ULS->Write();
    h_PtYPdg_DiMu_Gen_Barion_ULS->Write();

    if (!fOut.GetDirectory("DiMu_Rec"))
        fOut.mkdir("DiMu_Rec");

    fOut.cd("DiMu_Rec");
    h_PtMPdg_DiMu_Rec_Meson_ULS->Write();
    h_PtMPdg_DiMu_Rec_Barion_ULS->Write();

    h_PtYPdg_DiMu_Rec_Meson_ULS->Write();
    h_PtYPdg_DiMu_Rec_Barion_ULS->Write();

    fOut.Close();

    return;
}