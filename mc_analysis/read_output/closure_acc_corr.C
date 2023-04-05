#include "/home/michele_pennisi/cernbox/common_include.h"

void closure_acc_corr()
{
    TString Origin[2] = {"Charm", "Beauty"};
    TString Type[2] = {"Meson", "Barion"};
    Int_t i_Type = 0;
    Int_t i_Origin = 0;

    TFile *fIn_Sim=new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/root_files/test/HF_Analysis_MCsim_Big.root", "READ");
    fIn_Sim->ls();
    TH3F *h_PtYPdg_DiMu_Rec_ULS = (TH3F *)fIn_Sim->Get(Form("DiMu_Gen/h_PtMPdg_DiMu_Gen_%s_ULS", Type[0].Data()));
    TH3F *h_PtYPdg_DiMu_Rec_ULS_Cloned = (TH3F *)h_PtYPdg_DiMu_Rec_ULS->Clone(Form("h_PtMPdg_DiMu_Gen_%s_%s_ULS_M4cut", Origin[i_Origin].Data(), Type[i_Type].Data()));
    h_PtYPdg_DiMu_Rec_ULS_Cloned->GetZaxis()->SetRange(1, 1);
    h_PtYPdg_DiMu_Rec_ULS_Cloned->GetXaxis()->SetRangeUser(0, 10);
    h_PtYPdg_DiMu_Rec_ULS_Cloned->GetYaxis()->SetRangeUser(4, 9);
    h_PtYPdg_DiMu_Rec_ULS_Cloned->Draw("LEGO");
    // return;
    // if (Origin[i_Origin].Contains("Charm"))
    // else if (Origin[i_Origin].Contains("Beauty"))
    //     h_PtYPdg_DiMu_Rec_ULS_Cloned->GetZaxis()->SetRange(2, 2);

    TH1F *h_Pt_DiMu_Rec_ULS = (TH1F *)h_PtYPdg_DiMu_Rec_ULS_Cloned->Project3D("xe");
    h_Pt_DiMu_Rec_ULS->RebinX(5);
    
    
    // h_Pt_DiMu_Rec_ULS->RebinY(20);

    TFile *fIn_Corr=new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/root_files/test/HF_Acc_corr_MCsim_Big.root", "READ");

    TH3F *h_PtYPdg_DiMu_Rec_ULS_Corr = (TH3F *)fIn_Corr->Get(Form("DiMu_corr/h_PtYPdg_DiMu_Corr_%s_ULS_M4cut", Type[0].Data()));
    
    TH3F *h_PtYPdg_DiMu_Rec_ULS_Cloned_Corr = (TH3F *)h_PtYPdg_DiMu_Rec_ULS_Corr->Clone(Form("h_PtYPdg_DiMu_Rec_%s_%s_ULS_M4cut", Origin[i_Origin].Data(), Type[i_Type].Data()));
    if (Origin[i_Origin].Contains("Charm"))
        h_PtYPdg_DiMu_Rec_ULS_Cloned_Corr->GetZaxis()->SetRange(1, 1);
    else if (Origin[i_Origin].Contains("Beauty"))
        h_PtYPdg_DiMu_Rec_ULS_Cloned_Corr->GetZaxis()->SetRange(2, 2);
    
    TH1F *h_Pt_DiMu_Rec_ULS_Corr = (TH1F *)h_PtYPdg_DiMu_Rec_ULS_Cloned_Corr->Project3D("xe");
    
    h_Pt_DiMu_Rec_ULS_Corr->RebinX(5);
    h_Pt_DiMu_Rec_ULS_Corr->SetLineColor(kRed);
    h_Pt_DiMu_Rec_ULS_Corr->SetMarkerColor(kRed);


    h_Pt_DiMu_Rec_ULS->Draw();
    h_Pt_DiMu_Rec_ULS_Corr->Draw("PESAME");
    TCanvas *canvas = new TCanvas(Form("canvas_%s_%s_corr_test", Origin[i_Origin].Data(), Type[i_Type].Data()), Form("canvas_%s_%s", Origin[i_Origin].Data(), Type[i_Type].Data()), 1000, 1000);
    canvas->cd();
    return;
    canvas->Divide(1, 2);
    canvas->cd(1);

    

    canvas->cd(2);

    TH1F *ratio = (TH1F *)h_Pt_DiMu_Rec_ULS_Corr->Clone(Form("ratio_%s_%s_ULS_M4cut", Origin[i_Origin].Data(), Type[i_Type].Data()));
    ratio->SetName(Form("ratio_%s_%s_ULS_M4cut", Origin[i_Origin].Data(), Type[i_Type].Data()));
    ratio->Divide(h_Pt_DiMu_Rec_ULS);
    ratio->Draw();
}