#if !defined(__CINT__) || defined(__CLING__)
#include <iostream>
#include <fstream>
// #include "Riostream.h"
#include <sstream>
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TRatioPlot.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#endif
// Macro riferita al file di simulazione DY_sigle_new.cc. Consente di analizzare le sim prodotte con questo + il confronto con il vecchio sistema di simulazione.

using namespace std;

void fast_test(Int_t choice = 1)
{

    TString type("fromBeauty");

    TFile *fIn_data = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/root_files/test/HF/HistLite_HF_MCDimuHFTree_294009.root", "READ");

    TH2D *hDimuPtM_data_old = (TH2D *)fIn_data->Get(Form("DiMuon/M4/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromHF"));
    printf("n dimu tot %0.0f\n", hDimuPtM_data_old->GetEntries());
    TH1D *hDimuPt_data_old = hDimuPtM_data_old->ProjectionX();
    // hDimuPt_data_old->Rebin(10);
    hDimuPt_data_old->SetMarkerSize(1.0);
    hDimuPt_data_old->SetMarkerStyle(20);
    hDimuPt_data_old->SetMarkerColor(kMagenta + 2);
    hDimuPt_data_old->SetLineColor(kMagenta + 2);
    hDimuPt_data_old->SetLineWidth(2);
    TH1D *hDimuM_data_old = hDimuPtM_data_old->ProjectionY();
    hDimuM_data_old->SetMarkerSize(1.0);
    hDimuM_data_old->SetMarkerStyle(20);
    hDimuM_data_old->SetMarkerColor(kMagenta + 2);
    hDimuM_data_old->SetLineColor(kMagenta + 2);
    hDimuM_data_old->SetLineWidth(2);

    TFile *fIn_new = TFile::Open("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/root_files/test/HF/Tree_HF_MCDimuHFTree_294009.root", "UPDATE");
    // h_PtMdiMu_CMUL7_old_cut_LS;1
    TTree *hDimuPtM_data_new = (TTree *)fIn_new->Get("Rec/rec_tree_muhf");

    // printf("n dimu tot %0.0f\n", hDimuPtM_data_new->GetEntries());

    hDimuPtM_data_new->Draw("m>>h_Pt_tree(260,4.,30.)");

    new TCanvas();
    TH1F *h_Pt_tree = (TH1F *)gDirectory->Get("h_Pt_tree");
    // printf("dhdhdh %0.0f\n",h_Pt_tree->GetNbinsX());
    // h_Pt_tree->Divide(hDimuPt_data_old);
    hDimuM_data_old->Draw("PE");
    h_Pt_tree->SetMarkerSize(1.5);
    h_Pt_tree->SetMarkerStyle(24);
    h_Pt_tree->SetMarkerColor(kMagenta + 2);
    h_Pt_tree->SetLineColor(kMagenta + 2);
    h_Pt_tree->SetLineWidth(2);
    h_Pt_tree->Draw("PESAME");

    /*
    TH1D *hDimuPt_data_new = hDimuPtM_data_new->ProjectionX();
    hDimuPt_data_new->Rebin(10);

    hDimuPt_data_new->SetMarkerSize(1.5);
    hDimuPt_data_new->SetMarkerStyle(20);
    hDimuPt_data_new->SetMarkerColor(kOrange + 2);
    hDimuPt_data_new->SetLineColor(kOrange + 2);
    hDimuPt_data_new->SetLineWidth(2);
    TH1D *hDimuM_data_new = hDimuPtM_data_new->ProjectionY();

    TCanvas *pt_test = new TCanvas("pt_test", "pt_test", 1000, 1000);
    pt_test->cd();
    hDimuPt_data_old->Draw("PE");
    hDimuPt_data_new->Draw("PESAME");
    TCanvas *pt_test_ratio = new TCanvas("pt_test_ratio", "pt_test_ratio", 1000, 1000);
    pt_test_ratio->cd();
    TH1D *hDimuPt_data_old_clone = (TH1D *)hDimuPt_data_old->Clone("hDimuPt_data_old_clone");
    hDimuPt_data_old_clone->Divide(hDimuPt_data_new);
    hDimuPt_data_old_clone->Draw("PE");
    */
}