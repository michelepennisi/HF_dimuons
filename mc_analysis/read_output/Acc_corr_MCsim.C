#include "Analysis_MCsim.h"

void Acc_corr_MCsim(
    TString RunMode = "HF",
    Int_t RunNumber = 294925,
    // Int_t RunNumber = 294241,
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
        dir_fileout.Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s/Acc_corr_MCsim", Task_Version.Data(), RunMode.Data()); // official with files saved locally

    printf("%s/%s\n", dir_fileout.Data(), fileout.Data());

    TString filename;
    filename.Form("%s_MCDimuHFTree_%d.root", RunMode.Data(), RunNumber);

    TString file_out;
    file_out.Form("%s_Acc_corr_MCsim_%d.root", RunMode.Data(), RunNumber);

    printf("Input File: %s\n", filename.Data());
    printf("Saving in dir: %s \nFile: %s\n", dir_fileout.Data(), file_out.Data());

    TString file_corr;
    TString dir_file_corr;
    // file_corr.Form("%s_Analysis_MCsim_%d_Acc_Table.root", RunMode.Data(), RunNumber);
    dir_file_corr.Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s/Analysis_MCsim", Task_Version.Data(), RunMode.Data()); // official with files saved locally
    file_corr.Form("%s_Analysis_MCsim_Acc_Table_merged.root", RunMode.Data());

    printf("Saving in dir: %s \nFile: %s\n", dir_file_corr.Data(), file_corr.Data());

    TFile fIn_corr(Form("%s/%s", dir_file_corr.Data(), file_corr.Data()), "READ");
    fIn_corr.ls();
    fIn_corr.cd();
    TH2D *h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm = (TH2D *)fIn_corr.Get("DiMu_Corr/h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm");
    h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm->SetDirectory(0);

    TH2D *h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty = (TH2D *)fIn_corr.Get("DiMu_Corr/h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty");
    h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty->SetDirectory(0);

    TH2D *h_PtY_Muon_Charm_Acc = (TH2D *)fIn_corr.Get("Muon_Corr/h_PtY_Muon_Charm_Acc");
    h_PtY_Muon_Charm_Acc->SetDirectory(0);

    TH2D *h_PtY_Muon_Beauty_Acc = (TH2D *)fIn_corr.Get("Muon_Corr/h_PtY_Muon_Beauty_Acc");
    h_PtY_Muon_Beauty_Acc->SetDirectory(0);

    fIn_corr.Close();
    Set_Hist();
    TChain *input_tree = Importing_Tree(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s", Task_Version.Data(), RunMode.Data()), filename);
    input_tree->ls();
    Int_t total_entries = input_tree->GetEntries();
    for (Int_t i_Event = 0; i_Event < input_tree->GetEntries(); i_Event++)
    {
        if (i_Event % 500000 == 0)
        {
            progress_status(i_Event, total_entries);
        }
        input_tree->GetEntry(i_Event);

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

            Bool_t DQ_Dimuon = kFALSE;
            if ((Y_DiMu > -4.0 && Y_DiMu < -2.5) && (Eta_Mu0 > -4.0 && Eta_Mu0 < -2.5) && (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5))
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
            Double_t DiMu_weight_DiMu_Corr = 999;
            Double_t DiMu_weight_Muon_Corr = 999;

            Bool_t Verbose = kFALSE;

            if (Charm_mu0 && Charm_mu1)
            {
                DiMu_origin_Selection[0] = kTRUE;
                if (Charm_mu0_Meson && Charm_mu1_Meson)
                    type_DiMu = 0;
                if (Charm_mu0_Barion && Charm_mu1_Barion)
                    type_DiMu = 1;
                if ((Charm_mu0_Meson && Charm_mu1_Barion) || (Charm_mu1_Meson && Charm_mu0_Barion))
                    type_DiMu = 2;
                Double_t Bin_Pt_Mu0 = h_PtY_Muon_Charm_Acc->GetXaxis()->FindBin(Pt_Mu0);
                Double_t Bin_Y_Mu0 = h_PtY_Muon_Charm_Acc->GetYaxis()->FindBin(Y_Mu0);

                Double_t Bin_Pt_Mu1 = h_PtY_Muon_Charm_Acc->GetXaxis()->FindBin(Pt_Mu1);
                Double_t Bin_Y_Mu1 = h_PtY_Muon_Charm_Acc->GetYaxis()->FindBin(Y_Mu1);

                Double_t Weight_Mu0 = h_PtY_Muon_Charm_Acc->GetBinContent(Bin_Pt_Mu0, Bin_Y_Mu0);

                Double_t Weight_Mu1 = h_PtY_Muon_Charm_Acc->GetBinContent(Bin_Pt_Mu1, Bin_Y_Mu1);

                DiMu_weight_Muon_Corr = Weight_Mu0 * Weight_Mu1;

                // h_PtY_Muon_Charm_Acc->Draw("COLZ");
                // return;

                Double_t Bin_Pt = h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm->GetXaxis()->FindBin(Pt_DiMu);
                Double_t Bin_Y = h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm->GetYaxis()->FindBin(Y_DiMu);

                DiMu_weight_DiMu_Corr = h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm->GetBinContent(Bin_Pt, Bin_Y);

                if (Verbose)
                {
                    printf("Pt_Mu0 %0.5f || Y_Mu0 %0.5f \n", Pt_Mu0, Y_Mu0);
                    printf("Bin_Pt_Mu0 %0.0f || Bin_Y_Mu0 %0.0f \n", Bin_Pt_Mu0, Bin_Y_Mu0);
                    printf("Weight_Mu0 %0.3f \n", Weight_Mu0);
                    printf("Pt_Mu1 %0.5f || Y_Mu1 %0.5f \n", Pt_Mu1, Y_Mu1);
                    printf("Bin_Pt_Mu1 %0.0f || Bin_Y_Mu1 %0.0f \n", Bin_Pt_Mu1, Bin_Y_Mu1);
                    printf("Weight_Mu1 %0.3f \n", Weight_Mu1);

                    printf("DiMu_weight_Muon_Corr %0.2f\n", DiMu_weight_Muon_Corr);

                    printf("Pt_DiMu %0.5f || Y_DiMu %0.5f \n", Pt_DiMu, Y_DiMu);
                    printf("Bin_Pt %005f || Bin_Y %0.0f \n", Bin_Pt, Bin_Y);
                    printf("DiMu_weight_DiMu_Corr %0.2f\n", DiMu_weight_DiMu_Corr);
                }
                // h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm->Draw("COLZ");
                // return;
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

                Double_t Bin_Pt_Mu0 = h_PtY_Muon_Beauty_Acc->GetXaxis()->FindBin(Pt_Mu0);
                Double_t Bin_Y_Mu0 = h_PtY_Muon_Beauty_Acc->GetYaxis()->FindBin(Y_Mu0);

                Double_t Bin_Pt_Mu1 = h_PtY_Muon_Beauty_Acc->GetXaxis()->FindBin(Pt_Mu1);
                Double_t Bin_Y_Mu1 = h_PtY_Muon_Beauty_Acc->GetYaxis()->FindBin(Y_Mu1);

                Double_t Weight_Mu0 = h_PtY_Muon_Beauty_Acc->GetBinContent(Bin_Pt_Mu0, Bin_Y_Mu0);

                Double_t Weight_Mu1 = h_PtY_Muon_Beauty_Acc->GetBinContent(Bin_Pt_Mu1, Bin_Y_Mu1);

                DiMu_weight_Muon_Corr = Weight_Mu0 * Weight_Mu1;

                Double_t Bin_Pt = h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty->GetXaxis()->FindBin(Pt_DiMu);
                Double_t Bin_Y = h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty->GetYaxis()->FindBin(Y_DiMu);

                DiMu_weight_DiMu_Corr = h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty->GetBinContent(Bin_Pt, Bin_Y);

                if (Verbose)
                {
                    printf("Pt_Mu0 %0.5f || Y_Mu0 %0.5f \n", Pt_Mu0, Y_Mu0);
                    printf("Bin_Pt_Mu0 %0.0f || Bin_Y_Mu0 %0.0f \n", Bin_Pt_Mu0, Bin_Y_Mu0);
                    printf("Weight_Mu0 %0.3f \n", Weight_Mu0);
                    printf("Pt_Mu1 %0.5f || Y_Mu1 %0.5f \n", Pt_Mu1, Y_Mu1);
                    printf("Bin_Pt_Mu1 %0.0f || Bin_Y_Mu1 %0.0f \n", Bin_Pt_Mu1, Bin_Y_Mu1);
                    printf("Weight_Mu1 %0.3f \n", Weight_Mu1);

                    printf("DiMu_weight_Muon_Corr %0.2f\n", DiMu_weight_Muon_Corr);

                    h_PtY_Muon_Beauty_Acc->Draw("COLZ");

                    printf("Pt_DiMu %0.5f || Y_DiMu %0.5f \n", Pt_DiMu, Y_DiMu);
                    printf("Bin_Pt %005f || Bin_Y %0.0f \n", Bin_Pt, Bin_Y);
                    printf("DiMu_weight_DiMu_Corr %0.2f\n", DiMu_weight_DiMu_Corr);
                }
                // h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty->Draw("COLZ");

                // return;
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

            for (Int_t i_DiMu_origin = 0; i_DiMu_origin < n_DiMu_origin - 3; i_DiMu_origin++)
            {

                if (DiMu_origin_Selection[i_DiMu_origin])
                {
                    if (Charge_DiMu == 0)
                    {
                        if ((M_DiMu > 4 && M_DiMu < 9) && Pt_DiMu < 10)
                        {
                            h_PtYPdg_DiMu_DiMu_Corr_ULS_M49_Pt010[i_DiMu_origin]->Fill(Pt_DiMu, Y_DiMu,type_DiMu, DiMu_weight_DiMu_Corr);
                            h_PtYPdg_DiMu_Muon_Corr_ULS_M49_Pt010[i_DiMu_origin]->Fill(Pt_DiMu, Y_DiMu,type_DiMu, DiMu_weight_Muon_Corr);
                        }
                    }
                }
            }
        }
    }

    TFile fOut(Form("%s/%s", dir_fileout.Data(), file_out.Data()), "RECREATE");
    fOut.cd();
    if (!fOut.GetDirectory("DiMu_corr"))
        fOut.mkdir("DiMu_corr");

    if (!fOut.GetDirectory("Muon_corr"))
        fOut.mkdir("Muon_corr");
        
    for (Int_t i_DiMu_origin = 0; i_DiMu_origin < n_DiMu_origin - 3; i_DiMu_origin++)
    {
        fOut.cd("DiMu_corr");
        h_PtYPdg_DiMu_DiMu_Corr_ULS_M49_Pt010[i_DiMu_origin]->Write();
        fOut.cd("Muon_corr");
        h_PtYPdg_DiMu_Muon_Corr_ULS_M49_Pt010[i_DiMu_origin]->Write();
    }

    if (!fOut.GetDirectory("Single_Muon"))
        fOut.mkdir("Single_Muon");

    fOut.cd("Single_Muon");
}
