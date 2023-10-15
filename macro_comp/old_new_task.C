#include "/home/michele_pennisi/cernbox/common_include.h"

void old_new_task(){
    TFile *fIn_oldTask=new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Hist_fromSim/Version1/HF/HistLite_HF_MCDimuHFTree_294009.root","READ");
    // fIn_oldTask->cd("DiMuon/M4/Rec_DQ_cut_match_LT/ULS");
    // fIn_oldTask->ls();

    // TH2F *h_PtYDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask=(TH2F*) fIn_oldTask->Get("DiMuon/M4/Rec_DQ_cut_match_LT/ULS/h_PtYDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm");
    // TH1F *h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask=(TH1F *)h_PtYDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->ProjectionX();
    // TH1F *h_YDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask=(TH1F *)h_PtYDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->ProjectionY();
    TH2F *h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask=(TH2F*) fIn_oldTask->Get("DiMuon/M4/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm");
    TH1F *h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask=(TH1F *)h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->ProjectionX();
    TH1F *h_MDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask=(TH1F *)h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->ProjectionY();

    TFile *fIn_newTask=new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/read_output/test/LHC22b3_MC_output_Hist_294009.root","READ");
    fIn_newTask->cd("DiMuon_Rec");
    fIn_newTask->ls();

    // TH2F *h_PtYDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask=(TH2F*) fIn_newTask->Get("DiMuon_Rec/h_PtY_DiMuon_Rec_Charm");

    // TH1F *h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask=(TH1F *)h_PtYDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask->ProjectionX();
    // TH1F *h_YDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask=(TH1F *)h_PtYDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask->ProjectionY();

    TH2F *h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask=(TH2F*) fIn_newTask->Get("DiMuon_Rec/h_PtM_DiMuon_Rec_Charm");
    h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask->GetYaxis()->SetRangeUser(4,30);
    TH1F *h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask=(TH1F *)h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask->ProjectionX();
    h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask->SetLineColor(kRed);
    TH1F *h_MDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask=(TH1F *)h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask->ProjectionY();
    h_MDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask->SetLineColor(kRed);

    TCanvas *c=new TCanvas("c","c",1800,1200);
    c->Divide(2,2);
    c->cd(1);
    gPad->SetLogy();
    h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->DrawCopy();
    h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask->Draw("SAME");
    c->cd(2);
    gPad->SetLogy();
    h_MDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->DrawCopy();
    h_MDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask->Draw("SAME");
    c->cd(3);
    h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->Divide(h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask);
    h_PtDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->Draw();
    c->cd(4);
    h_MDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->Divide(h_MDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_newTask);
    h_MDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm_oldTask->Draw("PE");
}