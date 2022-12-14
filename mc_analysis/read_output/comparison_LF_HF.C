#if !defined(__CINT__) || defined(__CLING__)
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TClassTable.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLine.h>
#include <TGrid.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliPDG.h"
#include "AliMC.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TClonesArray.h"
#include "TRandom3.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TH3D.h"
#include "THStack.h"
#include "TLegend.h"

#include "TGraphAsymmErrors.h"
// #include "/Users/michelepennisi/cernbox/Analisi/Graphics.h"
#endif

TH1D *Get1Dhist(const char *filename, const char *hist_name, bool norm);
TH2D *Get2Dhist(const char *filename, const char *hist_name);
TCanvas *allprint2hist(TH1D *h1, TH1D *h2, double h_grid_minx, Bool_t norm);
TCanvas *allprint3hist(TH1D *h1, TH1D *h2, TH1D *h3, double h_grid_minx, double h_grid_maxx, TString Info);

TCanvas *allprint4hist(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4, double h_grid_minx, double h_grid_maxx, TString Info);
TCanvas *allprint5hist(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4, TH1D *h5, double h_grid_minx, double h_grid_maxx, TString Info);
double funcPt(double *x, double *par){
    return par[0] * x[0]/TMath::Power(1 + TMath::Power(x[0]/par[1],par[2]),par[3]);
}
// dimu_gen_masscut0/dimu_gen_ycut/Al

// dimu_gen_masscut0/dumyu/Al

//------------------------------
// Double_t Mass_cut=0.0,
// TString ptHist_name="dimu_gen_masscut0/dimu_gen_ycut/All/h_PtYdiMu_gen_ycut",
// TString massHist_name="dimu_gen_masscut0/dimu_gen_ycut/All/h_MdiMu_gen_ycut",
// TString Info="Gen. #mu#mu, -4.0<y_{#mu}<-2.5, Stat. Unc. Only",
// TString Canvas="Gen_0cut"

// Double_t Mass_cut=4.0,
// TString ptHist_name="dimu_gen_masscut4/dimu_gen_ycut/All/h_PtYdiMu_gen_ycut",
// TString massHist_name="dimu_gen_masscut4/dimu_gen_ycut/All/h_MdiMu_gen_ycut",
// TString Info="Gen. #mu#mu, -4.0<y_{#mu}<-2.5, Stat. Unc. Only",
// TString Canvas="Gen_4cut"
//...............,,,,,,,,,
// Double_t Mass_cut=0.0,
// TString ptHist_name="dimu_gen_masscut0/dimu_rec/All/h_PtYdiMu_rec",
// TString massHist_name="dimu_gen_masscut0/dimu_rec/All/h_MdiMu_rec",
// TString Info="Rec. #mu#mu, -4.0<y_{#mu}<-2.5, Stat. Unc. Only",
// TString Canvas="Rec_0Cut"

// Double_t Mass_cut=4.0,
// TString ptHist_name="dimu_gen_masscut4/dimu_rec/All/h_PtYdiMu_rec",
// TString massHist_name="dimu_gen_masscut4/dimu_rec/All/h_MdiMu_rec",
// TString Info="Rec. #mu#mu, -4.0<y_{#mu}<-2.5, Stat. Unc. Only",
// TString Canvas="Rec_4Cut"

void Check_DimuLF_HF(
    Double_t Mass_cut = 0.0,
    TString ptHist_name = "DiMu_M0/Rec_data_sel/ULS/h_PtYDiMu_M0_Rec_data_sel_ULS",
    TString massHist_name = "DiMu_M0/Rec_data_sel/ULS/h_MDiMu_M0_Rec_data_sel_ULS",
    TString Info = "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu} < 4.0",
    TString Canvas = "Rec_0Cut")
{

    const Int_t N_bin=15;
    Double_t new_bin[N_bin] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0};

    Int_t LF_events = 47967000;
    Double_t LF_sec = 78.05; // PYTHIA integrated cross section in mb
    Double_t LF_eqPythia = 1;
    //-----------------------------------------//
    Int_t HF_events = 2 * TMath::Power(10, 8);
    Double_t HF_sec = 56.42; // PYTHIA integrated cross section in mb
    Double_t HF_eqPythia = 216.1;

    Double_t LF_norm = LF_sec / (HF_sec * LF_events * LF_eqPythia);
    TString LF_filename;
    // LF from local simulation
    LF_filename.Form("~/dimuon_HF_pp/data/MBHistMCDimuHFTree_merged.root");

    TH2D *LFh_PtYDimu_ycut = Get2Dhist(LF_filename.Data(), Form("%s_fromLF", ptHist_name.Data()));
    TH1D *LFh_Pt_ycut = LFh_PtYDimu_ycut->ProjectionX();
    LFh_Pt_ycut->Sumw2();
    LFh_Pt_ycut->SetTitle("#mu^{#plus}#mu^{#minus} #leftarrow Sum LF");
    LFh_Pt_ycut->SetMarkerSize(1.5);
    LFh_Pt_ycut->SetMarkerStyle(20);
    LFh_Pt_ycut->SetMarkerColor(kCyan + 2);
    LFh_Pt_ycut->SetLineColor(kCyan + 2);
    LFh_Pt_ycut->SetLineWidth(2);

    TH1D *LFh_Mass_ycut = Get1Dhist(LF_filename.Data(), Form("%s_fromLF", massHist_name.Data()), false);
    LFh_Mass_ycut->Sumw2();
    LFh_Mass_ycut->SetTitle("#mu^{#plus}#mu^{#minus} #leftarrow Sum LF");
    LFh_Mass_ycut->SetMarkerSize(1.5);
    LFh_Mass_ycut->SetMarkerStyle(20);
    LFh_Mass_ycut->SetMarkerColor(kCyan + 2);
    LFh_Mass_ycut->SetLineColor(kCyan + 2);
    LFh_Mass_ycut->SetLineWidth(2);

    TH1D *LFh_Pt_ycut_rebinned = dynamic_cast<TH1D *>(LFh_Pt_ycut->Rebin(N_bin-1, "LFh_Pt_ycut_rebinned", new_bin));
    LFh_Pt_ycut_rebinned->Scale(1. / LF_events, "width");

    TH1D *LFh_Mass_ycut_rebinned = dynamic_cast<TH1D *>(LFh_Mass_ycut->Rebin(N_bin-1, "LFh_Mass_ycut_rebinned", new_bin));
    LFh_Mass_ycut_rebinned->Scale(1. / LF_events, "width");

    TH2D *LFMixedh_PtYDimu_ycut = Get2Dhist(LF_filename.Data(), Form("%s_fromLFMixed", ptHist_name.Data()));
    TH1D *LFMixedh_Pt_ycut = LFMixedh_PtYDimu_ycut->ProjectionX();
    LFMixedh_Pt_ycut->Sumw2();
    LFMixedh_Pt_ycut->SetTitle("#mu^{#plus}#mu^{#minus} #leftarrow LF, HF (incl. resonances)");
    LFMixedh_Pt_ycut->SetMarkerSize(1.5);
    LFMixedh_Pt_ycut->SetMarkerStyle(20);
    LFMixedh_Pt_ycut->SetMarkerColor(kMagenta + 2);
    LFMixedh_Pt_ycut->SetLineColor(kMagenta + 2);
    LFMixedh_Pt_ycut->SetLineWidth(2);

    TH1D *LFMixedh_Mass_ycut = Get1Dhist(LF_filename.Data(), Form("%s_fromLFMixed", massHist_name.Data()), false);
    LFMixedh_Mass_ycut->Sumw2();
    LFMixedh_Mass_ycut->SetTitle("#mu^{#plus}#mu^{#minus} #leftarrow LF, HF (incl. resonances)");
    LFMixedh_Mass_ycut->SetMarkerSize(1.5);
    LFMixedh_Mass_ycut->SetMarkerStyle(20);
    LFMixedh_Mass_ycut->SetMarkerColor(kMagenta + 2);
    LFMixedh_Mass_ycut->SetLineColor(kMagenta + 2);
    LFMixedh_Mass_ycut->SetLineWidth(2);
    TF1 *func = new TF1("func", funcPt, 0, 30, 4);
    func->SetParameter(0, 1e-6);
    func->SetParameter(1, 2.85);
    func->SetParameter(2, 2.81);
    func->SetParameter(3, 2.43);
    
    TH1D *LFMixedh_Pt_ycut_rebinned = dynamic_cast<TH1D *>(LFMixedh_Pt_ycut->Rebin(N_bin-1, "LFMixedh_Pt_ycut_rebinned", new_bin));
    LFMixedh_Pt_ycut_rebinned->Scale(1. / LF_events, "width");

    TH1D *LFMixedh_Mass_ycut_rebinned = dynamic_cast<TH1D *>(LFMixedh_Mass_ycut->Rebin(N_bin-1, "LFMixedh_Mass_ycut_rebinned", new_bin));
    LFMixedh_Mass_ycut_rebinned->Scale(1. / LF_events, "width");
    LFMixedh_Mass_ycut_rebinned -> Fit(func, "R");

    // TH1D *ALT_LFh_Pt_ycut_rebinned=dynamic_cast<TH1D*>(LFh_Pt_ycut->Rebin(N_bin-1,"LFh_Pt_ycut_rebinned",new_bin));
    // for (Int_t i = 0; i <= ALT_LFh_Pt_ycut_rebinned->GetNbinsX(); i++) {
    //   printf("Bin %d) Content %0.4f | Width %0.2f\n",i,ALT_LFh_Pt_ycut_rebinned->GetBinContent(i),ALT_LFh_Pt_ycut_rebinned->GetBinWidth(i));
    //   Double_t value=LF_norm*ALT_LFh_Pt_ycut_rebinned->GetBinContent(i)/ALT_LFh_Pt_ycut_rebinned->GetBinWidth(i);
    //   ALT_LFh_Pt_ycut_rebinned->SetBinContent(i,value);
    // }
    //
    // for (Int_t i = 0; i < ALT_LFh_Pt_ycut_rebinned->GetNbinsX(); i++) {
    //   printf("Bin %d) Content with Scale %0.4e | Content Manual %0.4e\n",i,LFh_Pt_ycut_rebinned->GetBinContent(i),ALT_LFh_Pt_ycut_rebinned->GetBinContent(i));
    // }
    // ALT_LFh_Pt_ycut_rebinned->SetMarkerColor(kBlack);
    //
    // TH1D *ALT_LFMixedh_Pt_ycut_rebinned=dynamic_cast<TH1D*>(LFMixedh_Pt_ycut->Rebin(N_bin-1,"LFMixedh_Pt_ycut_rebinned",new_bin));
    // for (Int_t i = 0; i <=ALT_LFMixedh_Pt_ycut_rebinned->GetNbinsX(); i++) {
    //   printf("Bin %d) Content %0.4f | Width %0.2f\n",i,ALT_LFMixedh_Pt_ycut_rebinned->GetBinContent(i),ALT_LFMixedh_Pt_ycut_rebinned->GetBinWidth(i));
    //   Double_t value=LF_norm*ALT_LFMixedh_Pt_ycut_rebinned->GetBinContent(i)/ALT_LFMixedh_Pt_ycut_rebinned->GetBinWidth(i);
    //   ALT_LFMixedh_Pt_ycut_rebinned->SetBinContent(i,value);
    // }
    // for (Int_t i = 0; i <=ALT_LFMixedh_Pt_ycut_rebinned->GetNbinsX(); i++) {
    //   printf("Bin %d) Content with Scale %0.4e | Content Manual %0.4e\n",i,LFMixedh_Pt_ycut_rebinned->GetBinContent(i),ALT_LFMixedh_Pt_ycut_rebinned->GetBinContent(i));
    // }
    // ALT_LFMixedh_Pt_ycut_rebinned->SetMarkerColor(kBlack);
    //
    // TCanvas *c_ptLF=allprint2hist(LFh_Pt_ycut_rebinned,ALT_LFh_Pt_ycut_rebinned, 0.0, false);
    // c_ptLF->SetName("c_ptLF");
    // TCanvas *c_ptLFMixed=allprint2hist(LFMixedh_Pt_ycut_rebinned,ALT_LFMixedh_Pt_ycut_rebinned, 0.0, false);
    // c_ptLFMixed->SetName("c_ptLFMixed");

    // HF from grid Simulation

    TString HF_filename_MB("~/dimuon_HF_pp/data/MBHistMCDimuHFTree_merged.root");

    TH2D *HF_MB_h_PtYDimu_ycut = Get2Dhist(HF_filename_MB.Data(), Form("%s_fromHF", ptHist_name.Data()));
    TH1D *HF_MB_h_Pt_ycut = HF_MB_h_PtYDimu_ycut->ProjectionX();
    HF_MB_h_Pt_ycut->Sumw2();
    HF_MB_h_Pt_ycut->SetTitle("#mu^{#plus}#mu^{#minus} #leftarrow Sum HF (incl. resonances)");
    HF_MB_h_Pt_ycut->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    HF_MB_h_Pt_ycut->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    HF_MB_h_Pt_ycut->SetMarkerSize(1.5);
    HF_MB_h_Pt_ycut->SetMarkerStyle(20);
    HF_MB_h_Pt_ycut->SetMarkerColor(kOrange + 8);
    HF_MB_h_Pt_ycut->SetLineColor(kOrange + 8);
    HF_MB_h_Pt_ycut->SetLineWidth(2);

    TH1D *HF_MB_h_Mass_ycut = Get1Dhist(HF_filename_MB.Data(), Form("%s_fromHF", massHist_name.Data()), false);
    HF_MB_h_Mass_ycut->Sumw2();
    HF_MB_h_Mass_ycut->SetTitle("#mu^{#plus}#mu^{#minus} #leftarrow Sum HF (incl. resonances)");
    HF_MB_h_Mass_ycut->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");
    HF_MB_h_Mass_ycut->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
    HF_MB_h_Mass_ycut->SetMarkerSize(1.5);
    HF_MB_h_Mass_ycut->SetMarkerStyle(20);
    HF_MB_h_Mass_ycut->SetMarkerColor(kOrange + 8);
    HF_MB_h_Mass_ycut->SetLineColor(kOrange + 8);
    HF_MB_h_Mass_ycut->SetLineWidth(2);

    TH1D *HF_MB_h_Pt_ycut_rebinned = dynamic_cast<TH1D *>(HF_MB_h_Pt_ycut->Rebin(N_bin-1, "HF_MB_h_Pt_ycut_rebinned", new_bin));
    HF_MB_h_Pt_ycut_rebinned->Scale(1. / LF_events, "width");

    TH1D *HF_MB_h_Mass_ycut_rebinned = dynamic_cast<TH1D *>(HF_MB_h_Mass_ycut->Rebin(N_bin-1, "HF_MB_h_Mass_ycut_rebinned", new_bin));
    HF_MB_h_Mass_ycut_rebinned->Scale(1. / LF_events, "width");

    // printf("Entries %0.1f \n",HF_MB_h_Pt_ycut->GetEntries());
    // printf("norm %0.7e \n",LF_norm);
    // printf("Manual %0.7e\n",HF_MB_h_Pt_ycut->GetEntries()/(LF_events));
    // printf("Intagral %0.7e\n",HF_MB_h_Pt_ycut_rebinned->Integral("width"));
    //
    // TString HF_filename("/media/michele_pennisi/DataBoi/Grid_Sim/HF/HistMCDimuHF/AllRun_HistMCDimuHFTree.root");
    //
    // TH2D *HFh_PtYDimu_ycut=Get2Dhist(HF_filename.Data(),Form("%s_fromHF",ptHist_name.Data()));
    // TH1D *HFh_Pt_ycut=HFh_PtYDimu_ycut->ProjectionX();
    // printf("Integral HF %0.10f\n", HFh_Pt_ycut->GetEntries()/(HF_events*HF_eqPythia));
    // printf("Entries %0.1f \n",HFh_Pt_ycut->GetEntries());
    // printf("norm %0.7e \n",1./(HF_events*HF_eqPythia));
    // printf("tot %0.7e\n",HFh_Pt_ycut->GetEntries()/(HF_events*HF_eqPythia) );
    // HFh_Pt_ycut->Sumw2();
    // HFh_Pt_ycut->SetTitle("#mu#mu #leftarrow Sum HF");
    // HFh_Pt_ycut->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} [GeV/#it{c}]^{-1}");
    // HFh_Pt_ycut->GetXaxis()->SetTitle("#it{p}_{T} [GeV/#it{c}]");
    // HFh_Pt_ycut->SetMarkerSize(1.5);
    // HFh_Pt_ycut->SetMarkerStyle(20);
    // HFh_Pt_ycut->SetMarkerColor(kOrange+8);
    // HFh_Pt_ycut->SetLineColor(kOrange+8);
    // HFh_Pt_ycut->SetLineWidth(2);
    //
    //
    // printf("Int LF MB %0.7e\n",LFh_Pt_ycut_rebinned->Integral("width")/LF_norm);
    // printf("Int LF MB %0.7e\n",LFMixedh_Pt_ycut_rebinned->Integral("width")/LF_norm);
    // printf("------------------------------------------------------------\n" );
    // TH1D *HFh_Mass_ycut=Get1Dhist(HF_filename.Data(),Form("%s_fromHF",massHist_name.Data()),false);
    // printf("Integral HF %0.10f\n", HFh_Mass_ycut->GetEntries()/(HF_events*HF_eqPythia));
    // printf("Entries %0.1f \n",HFh_Mass_ycut->GetEntries());
    // printf("norm %0.7e \n",1./(HF_events*HF_eqPythia));
    // printf("tot %0.7e\n",HFh_Mass_ycut->GetEntries()/(HF_events*HF_eqPythia) );
    // HFh_Mass_ycut->Sumw2();
    // printf("Int HF MB %0.7e\n",HFh_Mass_ycut->Integral("width")/LF_norm);
    // printf("Int LF MB %0.7e\n",LFh_Mass_ycut_rebinned->Integral("width")/LF_norm);
    // printf("Int LF MB %0.7e\n",LFMixedh_Mass_ycut_rebinned->Integral("width")/LF_norm);
    //
    // HFh_Mass_ycut->SetTitle("#mu#mu #leftarrow Sum HF");
    // HFh_Mass_ycut->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} [GeV/#it{c}^{2}]^{-1}");
    // HFh_Mass_ycut->GetXaxis()->SetTitle("#it{m}_{#mu#mu} [GeV/#it{c}^{2}]");
    // HFh_Mass_ycut->SetMarkerSize(1.5);
    // HFh_Mass_ycut->SetMarkerStyle(20);
    // HFh_Mass_ycut->SetMarkerColor(kOrange+8);
    // HFh_Mass_ycut->SetLineColor(kOrange+8);
    // HFh_Mass_ycut->SetLineWidth(2);
    //
    // TH1D *HFh_Pt_ycut_rebinned=dynamic_cast<TH1D*>(HFh_Pt_ycut->Rebin(14,"HFh_Pt_ycut_rebinned",new_bin));
    // HFh_Pt_ycut_rebinned->Scale(1./(HF_events*HF_eqPythia),"width");
    // printf("Int HF MB Pt%0.7e\n",HFh_Pt_ycut_rebinned->Integral("width")*LF_events);
    //
    // TH1D *HFh_Mass_ycut_rebinned=dynamic_cast<TH1D*>(HFh_Mass_ycut->Rebin(14,"HFh_Mass_ycut_rebinned",new_bin));
    // HFh_Mass_ycut_rebinned->Scale(1./(HF_events*HF_eqPythia),"width");

    // printf("Int HF MB Mass %0.7e\n",HFh_Mass_ycut_rebinned->Integral("width")*LF_events);
    TLatex *letexTitle = new TLatex();
    letexTitle->SetTextSize(0.045);
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);
    Double_t h_grid_minx = 0;

    if (Mass_cut == 0)
        h_grid_minx = -0.5;
    else if (Mass_cut == 4.0)
        h_grid_minx = 3.8;

    // TCanvas *cpt_LFHF=allprint4hist(HFh_Pt_ycut_rebinned,HF_MB_h_Pt_ycut_rebinned,LFh_Pt_ycut_rebinned,LFMixedh_Pt_ycut_rebinned,-0.5,20.0,Info);
    // cpt_LFHF->cd();
    // if(Mass_cut==4.0)letexTitle -> DrawLatex(0.225,0.68, Form("M_{#mu#mu}>%0.0f GeV/#it{c}^{2}",Mass_cut) );
    // cpt_LFHF->SetName(Form("%s_cpt_LFHF",Canvas.Data()));
    // cpt_LFHF->SetTitle(Form("%s_cpt_LFHF",Canvas.Data()));
    //
    // TCanvas *cmass_LFHF=allprint4hist(HFh_Mass_ycut_rebinned,HF_MB_h_Mass_ycut_rebinned,LFh_Mass_ycut_rebinned,LFMixedh_Mass_ycut_rebinned,h_grid_minx,20.0,Info);
    // cmass_LFHF->SetName(Form("%s_cmass_LFHF",Canvas.Data()));
    // cmass_LFHF->SetTitle(Form("%s_cmass_LFHF",Canvas.Data()));

    TCanvas *cpt_LFHF = allprint3hist(HF_MB_h_Pt_ycut_rebinned, LFh_Pt_ycut_rebinned, LFMixedh_Pt_ycut_rebinned, -0.5, 20.0, Info);
    cpt_LFHF->cd();
    if (Mass_cut == 4.0)
        letexTitle->DrawLatex(0.225, 0.68, Form("M_{#mu#mu}>%0.0f GeV/#it{c}^{2}", Mass_cut));
    cpt_LFHF->SetName(Form("%s_cpt_LFHF", Canvas.Data()));
    cpt_LFHF->SetTitle(Form("%s_cpt_LFHF", Canvas.Data()));

    TCanvas *cmass_LFHF = allprint3hist(HF_MB_h_Mass_ycut_rebinned, LFh_Mass_ycut_rebinned, LFMixedh_Mass_ycut_rebinned, h_grid_minx, 20.0, Info);
    cmass_LFHF->SetName(Form("%s_cmass_LFHF", Canvas.Data()));
    cmass_LFHF->SetTitle(Form("%s_cmass_LFHF", Canvas.Data()));

    // TCanvas *cpt_CBLF=allprint5hist(Charmh_Pt_ycut,Beautyh_Pt_ycut,HFMixedh_Pt_ycut,LFh_Pt_ycut_rebinned,LFMixedh_Pt_ycut_rebinned,-0.5,20.0,Info);
    // cpt_LFHF->cd();
    // if(Mass_cut==4.0)letexTitle -> DrawLatex(0.225,0.68, Form("M_{#mu#mu}>%0.0f GeV/#it{c}^{2}",Mass_cut) );
    // cpt_CBLF->SetName(Form("%s_cpt_CBLF",Canvas.Data()));
    // cpt_CBLF->SetTitle(Form("%s_cpt_CBLF",Canvas.Data()));
    //
    // TCanvas *cmass_CBLF=allprint5hist(Charmh_Mass_ycut,Beautyh_Mass_ycut,HFMixedh_Mass_ycut,LFh_Mass_ycut_rebinned,LFMixedh_Mass_ycut_rebinned,h_grid_minx,20.0,Info);
    // cmass_CBLF->SetName(Form("%s_cmass_CBLF",Canvas.Data()));
    // cmass_CBLF->SetTitle(Form("%s_cmass_CBLF",Canvas.Data()));

    printf("HF_MB_h_Pt_ycut (%.3f)\n", HF_MB_h_Pt_ycut->GetEntries());
    printf("LFh_Pt_ycut (%.3f)\n", LFh_Pt_ycut->GetEntries());
    printf("LFMixedh_Pt_ycut (%.3f)\n", LFMixedh_Pt_ycut->GetEntries());
    gSystem->cd("Image");
    cpt_LFHF->Print(Form("%s_ULSnormMB.eps", cpt_LFHF->GetName()));
    cmass_LFHF->Print(Form("%s_ULSnormMB.eps", cmass_LFHF->GetName()));
    cpt_LFHF->Print(Form("%s_ULSnormMB.pdf", cpt_LFHF->GetName()));
    cmass_LFHF->Print(Form("%s_ULSnormMB.pdf", cmass_LFHF->GetName()));
    // cpt_CBLF->Print(Form("%s.pdf",cpt_CBLF->GetName()));
    // cmass_CBLF->Print(Form("%s.pdf",cmass_CBLF->GetName()));

    TString fileout(Form("/home/michele_pennisi/cernbox/Read_AOD/Grid_Hist/Comparison%s", LF_filename.Data()));

    printf("Creating hist for %s, with Mass Cut = %0.1f\n", LF_filename.Data(), Mass_cut);
    printf("Output %s\n", fileout.Data());
}

TH1D *Get1Dhist(const char *filename, const char *hist_name, bool norm)
{

    TFile *myFile = TFile::Open(filename, "UPDATE");

    if (!myFile || myFile->IsZombie())
    {
        printf("%s non trovato \n", filename);
        //        cout<<"File non trovato"<<endl;
        return nullptr;
    }
    TH1D *h_hist = (TH1D *)myFile->Get(hist_name);
    if (!h_hist || h_hist->IsZombie())
    {
        printf("%s non trovato \n", hist_name);
        //        cout<<"Istogramma non presente o con nome errato"<<endl;
        return nullptr;
    }
    h_hist->SetDirectory(0);
    h_hist->SetMinimum(1);
    if (norm)
        h_hist->Scale(1, "width");
    myFile->Close();
    //    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptStat(0);
    h_hist->SetMarkerStyle(0);
    return h_hist;
}

TH2D *Get2Dhist(const char *filename, const char *hist_name)
{

    TFile *myFile = TFile::Open(filename, "UPDATE");

    if (!myFile || myFile->IsZombie())
    {
        printf("%s non trovato \n", filename);
        //        cout<<"File non trovato"<<endl;
        return nullptr;
    }
    TH2D *h_hist = (TH2D *)myFile->Get(hist_name);
    if (!h_hist || h_hist->IsZombie())
    {
        printf("%s non trovato \n", hist_name);
        //        cout<<"Istogramma non presente o con nome errato"<<endl;
        return nullptr;
    }
    h_hist->SetDirectory(0);
    h_hist->SetMinimum(1);
    myFile->Close();
    //    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptStat(0);
    h_hist->SetMarkerStyle(0);
    return h_hist;
}

TCanvas *allprint2hist(TH1D *h1, TH1D *h2, double h_grid_minx, Bool_t norm)
{

    TString type(h1->GetName());

    TCanvas *canvasSpectra = new TCanvas("canvasSpectra", "canvasSpectra", 1000, 800);
    canvasSpectra->SetName("canvasSpectra");
    canvasSpectra->cd();
    TLegend *legend = new TLegend(0.45, 0.45, 0.85, 0.7);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);

    legend->AddEntry(h1, Form("%s (%0.0f)", h1->GetTitle(), h1->GetEntries()));
    legend->AddEntry(h2, Form("%s (%0.0f)", h2->GetTitle(), h2->GetEntries()));

    if (norm)
    {

        // h1->Scale(1,"width");
        // h2->Scale(1,"width");
    }
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetLogy();

    h1->Draw("PE SAME");
    h2->Draw("PE SAME");

    legend->Draw();

    TLatex *letexTitle = new TLatex();
    letexTitle->SetTextSize(0.045);
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);

    letexTitle->DrawLatex(0.4, 0.88, "ALICE Performance,pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.4, 0.81, "PYTHIA8 Monash Tune, Normalized to N_{eq. Pythia}");
    letexTitle->DrawLatex(0.4, 0.74, "-4.0<y_{#mu}<-2.5, Stat. Unc. Only");
    // letexTitle -> DrawLatex(0.4,0.66,"M_{#mu#mu}>4 GeV/#it{c}^{2}");

    return canvasSpectra;
}

TCanvas *allprint3hist(TH1D *h1, TH1D *h2, TH1D *h3, double h_grid_minx, double h_grid_maxx, TString Info)
{

    TH2D *h_grid = new TH2D("h_grid", "", 100, h_grid_minx, h_grid_maxx, 10, 0.00000000001, h1->GetMaximum() * 1000);
    h_grid->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.045);
    h_grid->GetXaxis()->SetLabelSize(0.04);
    h_grid->GetXaxis()->SetNdivisions(505);
    h_grid->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
    h_grid->GetYaxis()->SetTitleOffset(1.3);
    h_grid->GetYaxis()->SetTitleSize(0.045);
    h_grid->GetYaxis()->SetLabelSize(0.04);
    h_grid->GetYaxis()->SetNdivisions(505);
    TCanvas *canvasSpectra = new TCanvas("canvasSpectra", "canvasSpectra", 1200, 1000);
    canvasSpectra->SetTicks();
    canvasSpectra->cd();
    TLegend *legend = new TLegend(0.5, 0.525, 0.7, 0.725);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0325);
    legend->SetFillStyle(0);

    // legend->AddEntry(h1,Form("%s (%.3e)",h1->GetTitle(),h1->Integral()));
    // legend->AddEntry(h2,Form("%s (%.3e)",h2->GetTitle(),h2->Integral()));
    // legend->AddEntry(h3,Form("%s (%.3e)",h3->GetTitle(),h3->Integral()));

    legend->AddEntry(h1, Form("%s", h1->GetTitle()));
    legend->AddEntry(h2, Form("%s", h2->GetTitle()));
    legend->AddEntry(h3, Form("%s", h3->GetTitle()));

    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.045);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetLogy();
    TLatex *letexTitle = new TLatex();
    letexTitle->SetTextSize(0.045);
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);

    h_grid->Draw();
    h1->Draw("PE SAME");
    h2->Draw("PE SAME");
    h3->Draw("PE SAME");

    letexTitle->DrawLatex(0.175, 0.885, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.175, 0.815, "PYTHIA8 Monash Tune, Norm. to #it{N}_{ev}= #it{N}_{Inel} #plus #it{N}_{Diff}");
    letexTitle->DrawLatex(0.175, 0.75, Info.Data());

    legend->Draw();

    return canvasSpectra;
}

TCanvas *allprint4hist(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4, double h_grid_minx, double h_grid_maxx, TString Info)
{
    TH2D *h_grid = new TH2D("h_grid", "", 100, h_grid_minx, h_grid_maxx, 10, 0.00000000001, h1->GetMaximum() * 100000);
    h_grid->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.045);
    h_grid->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
    h_grid->GetYaxis()->SetTitleOffset(1.3);
    h_grid->GetYaxis()->SetTitleSize(0.045);
    TCanvas *canvasSpectra = new TCanvas("canvasSpectra", "canvasSpectra", 1200, 1000);
    canvasSpectra->cd();
    TLegend *legend = new TLegend(0.45, 0.5, 0.65, 0.7);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);

    legend->AddEntry(h1, Form("%s (%.3e)", h1->GetTitle(), h1->Integral()));
    legend->AddEntry(h2, Form("%s (%.3e)", h2->GetTitle(), h2->Integral()));
    legend->AddEntry(h3, Form("%s (%.3e)", h3->GetTitle(), h3->Integral()));
    legend->AddEntry(h4, Form("%s (%.3e)", h4->GetTitle(), h4->Integral()));

    // legend->AddEntry(h1,Form("%s",h1->GetTitle()));
    // legend->AddEntry(h2,Form("%s",h2->GetTitle()));
    // legend->AddEntry(h3,Form("%s",h3->GetTitle()));

    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.045);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetLogy();
    TLatex *letexTitle = new TLatex();
    letexTitle->SetTextSize(0.045);
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);

    h_grid->Draw();
    h1->Draw("PE SAME");
    h2->Draw("PE SAME");
    h3->Draw("PE SAME");
    h4->Draw("PE SAME");

    letexTitle->DrawLatex(0.225, 0.885, "ALICE Simulation,pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.225, 0.815, "PYTHIA8 Monash Tune, Normalized to N_{eq. Pythia}");
    letexTitle->DrawLatex(0.225, 0.75, Info.Data());

    legend->Draw();

    return canvasSpectra;
}

TCanvas *allprint5hist(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4, TH1D *h5, double h_grid_minx, double h_grid_maxx, TString Info)
{

    TH2D *h_grid = new TH2D("h_grid", "", 100, h_grid_minx, h_grid_maxx, 10, 1e-11, h1->GetMaximum() * 2500);
    h_grid->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.045);
    h_grid->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
    h_grid->GetYaxis()->SetTitleOffset(1.3);
    h_grid->GetYaxis()->SetTitleSize(0.045);
    TCanvas *canvasSpectra = new TCanvas("canvasSpectra", "canvasSpectra", 1200, 1000);
    canvasSpectra->cd();
    TLegend *legend = new TLegend(0.45, 0.5, 0.7, 0.7);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);

    legend->AddEntry(h1, Form("%s (%.3e)", h1->GetTitle(), h1->Integral()));
    legend->AddEntry(h2, Form("%s (%.3e)", h2->GetTitle(), h2->Integral()));
    legend->AddEntry(h3, Form("%s (%.3e)", h3->GetTitle(), h3->Integral()));
    legend->AddEntry(h4, Form("%s (%.3e)", h4->GetTitle(), h4->Integral()));
    legend->AddEntry(h5, Form("%s (%.3e)", h5->GetTitle(), h5->Integral()));

    // legend->AddEntry(h1,Form("%s",h1->GetTitle()));
    // legend->AddEntry(h2,Form("%s",h2->GetTitle()));
    // legend->AddEntry(h3,Form("%s",h3->GetTitle()));
    // legend->AddEntry(h4,Form("%s",h4->GetTitle()));
    // legend->AddEntry(h5,Form("%s",h5->GetTitle()));

    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.045);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetLogy();
    TLatex *letexTitle = new TLatex();
    letexTitle->SetTextSize(0.04);
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);

    h_grid->Draw();
    h1->Draw("PE SAME");
    h2->Draw("PE SAME");
    h3->Draw("PE SAME");
    h4->Draw("PE SAME");
    h5->Draw("PE SAME");

    letexTitle->DrawLatex(0.225, 0.885, "ALICE Simulation,pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.225, 0.815, "PYTHIA8 Monash Tune, Normalized to N_{eq. Pythia}");
    letexTitle->DrawLatex(0.225, 0.75, Info.Data());
    legend->Draw();

    return canvasSpectra;
}
