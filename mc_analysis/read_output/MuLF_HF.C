#if !defined (__CINT__) || defined (__CLING__)
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

TH1D *Get1Dhist(const char *filename,const char *hist_name,bool norm);
TH2D *Get2Dhist(const char *filename,const char *hist_name);
TCanvas *allprint2hist(TH2D *h_grid,TH1D *h1,TH1D *h2,double h_grid_minx, Bool_t norm);
TCanvas *allprint3hist(TH1D *h1,TH1D *h2,TH1D *h3,double h_grid_minx,double h_grid_maxx,TString Info,Double_t Norm);

TCanvas *allprint4hist(TH1D *h1,TH1D *h2,TH1D *h3,TH1D *h4);
TCanvas *allprint5hist(TH1D *h1,TH1D *h2,TH1D *h3,TH1D *h4,TH1D *h5,double h_grid_minx,double h_grid_maxx,TString Info);
//dimu_gen_masscut0/dimu_gen_ycut/Al

//dimu_gen_masscut0/dumyu/Al

//------------------------------
// TString ptHist_name="dimu_gen_masscut0/dimu_gen_ycut/All/h_PtYdiMu_gen_ycut",
// TString massHist_name="dimu_gen_masscut0/dimu_gen_ycut/All/h_MdiMu_gen_ycut",
// TString Info="#splitline{Gen. #mu#mu, -4.0<y_{#mu}<-2.5, Stat. Unc. Only}{M_{#mu#mu}>0 GeV/c^{2}}",
// TString Canvas="Gen_0cut"
//...............,,,,,,,,,
// TString ptHist_name="dimu_gen_masscut0/dimu_rec/All/h_PtYdiMu_rec",
// TString massHist_name="dimu_gen_masscut0/dimu_rec/All/h_MdiMu_rec",
// TString Info="#splitline{Rec. #mu#mu, -4.0<y_{#mu}<-2.5, Stat. Unc. Only}{M_{#mu#mu}>0 GeV/c^{2}}",
// TString Canvas="Rec_0Cut"

void MuLF_HF(
  Double_t Mass_cut=4.0,
  Double_t h_grid_minx=-0.5,
  Int_t MC_choice=1 ,
  TString ptHist_name="muon/h_PtYmu_gen",
  TString Info="LF MC Gen. #mu, Stat. Unc. Only",
  TString Canvas="Gen_0cut"
){

  Int_t LF_events=2*TMath::Power(10,5);

  Double_t LF_sec=78.05; //PYTHIA integrated cross section in mb

  Double_t LF_eqPythia=201.9;

  Int_t HF_events=2*TMath::Power(10,8);

  Double_t HF_sec=56.42; //PYTHIA integrated cross section in mb

  Double_t HF_eqPythia=216.1;

  Double_t Norm=0.0;
  TString filename;
  if (MC_choice==0) {
    filename.Form("/media/michele_pennisi/DataBoi/Grid_Sim/HF/HistMCDimuHF/AllRun_HistMCDimuHFTree.root");
    Norm=(HF_events*HF_eqPythia);
  }else if (MC_choice==1) {
    filename.Form("/media/michele_pennisi/DataBoi/LF_study/HistLF_study_AnalysisDimuonHF_Runs0_1999_Ev198700.root");
    Norm=(HF_sec*LF_events*LF_eqPythia)/LF_sec;
  }else if (MC_choice==2) {
    filename.Form("/media/michele_pennisi/DataBoi/MB_study/HistMB_study_AnalysisDimuonHF_Runs0_2000_Ev199300.root");
    Norm=(HF_sec*LF_events)/LF_sec;
  }
  printf("Norm %0.5e \n",Norm);
  TH2D *LFh_PtYMu_ycut=Get2Dhist(filename.Data(),Form("%s_fromLF",ptHist_name.Data()));

  TH1D *LFh_Pt_ycut=LFh_PtYMu_ycut->ProjectionX();
  TH1D *LFh_Pt_ycut_norm=(TH1D*)LFh_Pt_ycut->Clone(Form("%s_norm",LFh_Pt_ycut->GetName()));

  LFh_Pt_ycut_norm->Sumw2();
  LFh_Pt_ycut_norm->Scale(1./Norm);
  LFh_Pt_ycut_norm->SetTitle("LF#rightarrow#mu");
  LFh_Pt_ycut_norm->SetMarkerSize(1.5);
  LFh_Pt_ycut_norm->SetMarkerStyle(20);
  LFh_Pt_ycut_norm->SetMarkerColor(kCyan+2);
  LFh_Pt_ycut_norm->SetLineColor(kCyan+2);
  LFh_Pt_ycut_norm->SetLineWidth(2);

  TH1D *LFh_Y_ycut=LFh_PtYMu_ycut->ProjectionY();
  TH1D *LFh_Y_ycut_norm=(TH1D*)LFh_Y_ycut->Clone(Form("%s_norm",LFh_Y_ycut->GetName()));

  LFh_Y_ycut_norm->Sumw2();
  LFh_Y_ycut_norm->Scale(1./Norm);
  LFh_Y_ycut_norm->SetTitle("LF#rightarrow#mu");
  LFh_Y_ycut_norm->SetMarkerSize(1.5);
  LFh_Y_ycut_norm->SetMarkerStyle(20);
  LFh_Y_ycut_norm->SetMarkerColor(kCyan+2);
  LFh_Y_ycut_norm->SetLineColor(kCyan+2);
  LFh_Y_ycut_norm->SetLineWidth(2);

  TH2D *HFh_PtYMu_ycut=Get2Dhist(filename.Data(),Form("%s_fromHF",ptHist_name.Data()));

  TH1D *HFh_Pt_ycut=HFh_PtYMu_ycut->ProjectionX();
  TH1D *HFh_Pt_ycut_norm=(TH1D*)HFh_Pt_ycut->Clone(Form("%s_norm",HFh_Pt_ycut->GetName()));

  HFh_Pt_ycut_norm->Sumw2();
  HFh_Pt_ycut_norm->Scale(1./Norm);
  // printf("Integral HF %0.10f\n", HFh_Mass_ycut->GetEntries()/(HF_events*HF_eqPythia));
  HFh_Pt_ycut_norm->SetTitle("HF#rightarrow#mu");
  HFh_Pt_ycut_norm->GetYaxis()->SetTitle("d#it{N}/dp_{#it{T}} [GeV/c]^{-1}");
  HFh_Pt_ycut_norm->GetXaxis()->SetTitle("p_{#it{T}} [GeV/c]");
  HFh_Pt_ycut_norm->SetMarkerSize(1.5);
  HFh_Pt_ycut_norm->SetMarkerStyle(20);
  HFh_Pt_ycut_norm->SetMarkerColor(kOrange+8);
  HFh_Pt_ycut_norm->SetLineColor(kOrange+8);
  HFh_Pt_ycut_norm->SetLineWidth(2);

  TH1D *HFh_Y_ycut=HFh_PtYMu_ycut->ProjectionY();
  TH1D *HFh_Y_ycut_norm=(TH1D*)HFh_Y_ycut->Clone(Form("%s_norm",HFh_Y_ycut->GetName()));

  HFh_Y_ycut_norm->Sumw2();
  HFh_Y_ycut_norm->Scale(1./Norm);
  // printf("Integral HF %0.10f\n", HFh_Mass_ycut->GetEntries()/(HF_events*HF_eqPythia));
  HFh_Y_ycut_norm->SetTitle("HF#rightarrow#mu");
  HFh_Y_ycut_norm->GetYaxis()->SetTitle("d#it{N}/dY");
  HFh_Y_ycut_norm->GetXaxis()->SetTitle("Y");
  HFh_Y_ycut_norm->SetMarkerSize(1.5);
  HFh_Y_ycut_norm->SetMarkerStyle(20);
  HFh_Y_ycut_norm->SetMarkerColor(kOrange+8);
  HFh_Y_ycut_norm->SetLineColor(kOrange+8);
  HFh_Y_ycut_norm->SetLineWidth(2);

  TH2D *Otherh_PtYMu_ycut=Get2Dhist(filename.Data(),Form("%s_fromOthers",ptHist_name.Data()));

  TH1D *Otherh_Pt_ycut=Otherh_PtYMu_ycut->ProjectionX();
  TH1D *Otherh_Pt_ycut_norm=(TH1D*)Otherh_Pt_ycut->Clone(Form("%s_norm",Otherh_Pt_ycut->GetName()));

  Otherh_Pt_ycut_norm->Sumw2();
  Otherh_Pt_ycut_norm->Scale(1./Norm);
  // printf("Integral Other %0.10f\n", Otherh_Mass_ycut->GetEntries()/(Other_events*Other_eqPythia));
  Otherh_Pt_ycut_norm->SetTitle("Other #rightarrow#mu");
  Otherh_Pt_ycut_norm->GetYaxis()->SetTitle("d#it{N}/dp_{#it{T}} [GeV/c]^{-1}");
  Otherh_Pt_ycut_norm->GetXaxis()->SetTitle("p_{#it{T}} [GeV/c]");
  Otherh_Pt_ycut_norm->SetMarkerSize(1.5);
  Otherh_Pt_ycut_norm->SetMarkerStyle(20);
  Otherh_Pt_ycut_norm->SetMarkerColor(kBlack);
  Otherh_Pt_ycut_norm->SetLineColor(kBlack);
  // Otherh_Pt_ycut_norm->SetLineWidth(2);

  TH1D *Otherh_Y_ycut=Otherh_PtYMu_ycut->ProjectionY();
  TH1D *Otherh_Y_ycut_norm=(TH1D*)Otherh_Y_ycut->Clone(Form("%s_norm",Otherh_Y_ycut->GetName()));

  Otherh_Y_ycut_norm->Sumw2();
  Otherh_Y_ycut_norm->Scale(1./Norm);
  // printf("Integral Other %0.10f\n", Otherh_Mass_ycut->GetEntries()/(Other_events*Other_eqPythia));
  Otherh_Y_ycut_norm->SetTitle("Other #rightarrow#mu");

  Otherh_Y_ycut_norm->SetMarkerSize(1.5);
  Otherh_Y_ycut_norm->SetMarkerStyle(20);
  Otherh_Y_ycut_norm->SetMarkerColor(kBlack);
  Otherh_Y_ycut_norm->SetLineColor(kBlack);
  Otherh_Y_ycut_norm->SetLineWidth(2);
  printf("%0.5e\n",Otherh_Y_ycut->GetEntries() );

// ULSptdimu_Mixed_rebinned->SetMarkerColor(kAzure+9);

  TH2D *h_grid;

  TCanvas *cpt_LFHF=allprint3hist(HFh_Pt_ycut_norm,LFh_Pt_ycut_norm,Otherh_Pt_ycut_norm,-0.2,30.2,Info,Norm);
  cpt_LFHF->SetName(Form("%s_cpt_LFHF",Canvas.Data()));
  cpt_LFHF->SetTitle(Form("%s_cpt_LFHF",Canvas.Data()));

  TCanvas *cy_LFHF=allprint3hist(HFh_Y_ycut_norm,LFh_Y_ycut_norm,Otherh_Y_ycut_norm,-10.0,14.0,Info,Norm);
  cy_LFHF->SetName(Form("%s_cy_LFHF",Canvas.Data()));
  cy_LFHF->SetTitle(Form("%s_cy_LFHF",Canvas.Data()));
  //
  //
  // TString fileout(Form("/home/michele_pennisi/cernbox/Read_AOD/Grid_Hist/Comparison%s",LF_filename.Data()));
  //
  //
  //
  // printf("Creating hist for %s, with Mass Cut = %0.1f\n", LF_filename.Data(),Mass_cut);
  // printf("Output %s\n",fileout.Data());
}


TH1D *Get1Dhist(const char *filename,const char *hist_name,bool norm) {

    TFile *myFile = TFile::Open(filename,"UPDATE");

    if (!myFile || myFile->IsZombie()) {
        printf("%s non trovato \n",filename);
//        cout<<"File non trovato"<<endl;
        return nullptr;
    }
    TH1D *h_hist = (TH1D*)myFile->Get(hist_name);
    if (!h_hist || h_hist->IsZombie()) {
        printf("%s non trovato \n",hist_name);
//        cout<<"Istogramma non presente o con nome errato"<<endl;
        return nullptr;
    }
    h_hist->SetDirectory(0);
    h_hist->SetMinimum(1);
    if (norm) h_hist->Scale(1,"width");
    myFile->Close();
//    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptStat(0);
    h_hist->SetMarkerStyle(0);
    return h_hist;
}

TH2D *Get2Dhist(const char *filename,const char *hist_name) {

    TFile *myFile = TFile::Open(filename,"UPDATE");

    if (!myFile || myFile->IsZombie()) {
        printf("%s non trovato \n",filename);
//        cout<<"File non trovato"<<endl;
        return nullptr;
    }
    TH2D *h_hist = (TH2D*)myFile->Get(hist_name);
    if (!h_hist || h_hist->IsZombie()) {
        printf("%s non trovato \n",hist_name);
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

TCanvas *allprint2hist(TH2D *h_grid,TH1D *h1,TH1D *h2,double h_grid_minx, Bool_t norm){
    h_grid=nullptr;
    h_grid=new TH2D("h_grid","",100,h_grid_minx,30.0,8,10e3,h1->GetMaximum()*100000);
    h_grid->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.045);
    h_grid->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
    h_grid->GetYaxis()->SetTitleOffset(1.3);
    h_grid->GetYaxis()->SetTitleSize(0.045);

    TString type(h1->GetName());

    TCanvas *canvasSpectra=new TCanvas ("canvasSpectra","canvasSpectra",1000,800);
    canvasSpectra->SetName("canvasSpectra");
    canvasSpectra->cd();
    TLegend *legend=new TLegend(0.45,0.45,0.85,0.7);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);

    legend->AddEntry(h1,Form("%s (%0.0f)",h1->GetTitle(),h1->GetEntries()));
    legend->AddEntry(h2,Form("%s (%0.0f)",h2->GetTitle(),h2->GetEntries()));

    if (norm) {

      h1->Scale(1,"width");
      h2->Scale(1,"width");

    }
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetLogy();

    h_grid->Draw();
    h1->Draw("PE SAME");
    h2->Draw("PE SAME");

    legend->Draw();

    TLatex *letexTitle = new TLatex();
    letexTitle -> SetTextSize(0.045);
    letexTitle -> SetNDC();
    letexTitle -> SetTextFont(42);

    letexTitle -> DrawLatex(0.4,0.88,"ALICE Performance,pp #sqrt{#it{s}} = 13 TeV");
    letexTitle -> DrawLatex(0.4,0.81,"PYTHIA8 Monash Tune, Normalized to N_{eq. Pythia}");
    letexTitle -> DrawLatex(0.4,0.74,"-4.0<y_{#mu}<-2.5, Stat. Unc. Only");
    // letexTitle -> DrawLatex(0.4,0.66,"M_{#mu#mu}>4 GeV/#it{c}^{2}");


    return canvasSpectra;
}

TCanvas *allprint3hist(TH1D *h1,TH1D *h2,TH1D *h3,double h_grid_minx,double h_grid_maxx,TString Info,Double_t Norm){

    TH2D* h_grid=new TH2D("h_grid","",100,h_grid_minx,h_grid_maxx,10,0.000000005,h2->GetMaximum()*1500.0);
    h_grid->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.045);
    h_grid->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
    h_grid->GetYaxis()->SetTitleOffset(1.3);
    h_grid->GetYaxis()->SetTitleSize(0.045);
    TCanvas *canvasSpectra=new TCanvas ("canvasSpectra","canvasSpectra",1200,1000);
    canvasSpectra->cd();
    TLegend *legend=new TLegend(0.45,0.49,0.65,0.69);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);

    legend->AddEntry(h1,Form("%s (%.5e)",h1->GetTitle(),(h1->GetEntries())/Norm));
    legend->AddEntry(h2,Form("%s (%.5e)",h2->GetTitle(),(h2->GetEntries())/Norm));
    legend->AddEntry(h3,Form("%s (%.5e)",h3->GetTitle(),(h3->GetEntries())/Norm));

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
    letexTitle -> SetTextSize(0.045);
    letexTitle -> SetNDC();
    letexTitle -> SetTextFont(42);

    h_grid->Draw();
    h1->Draw("PE SAME");
    h2->Draw("PE SAME");
    h3->Draw("PE SAME");

    letexTitle -> DrawLatex(0.225,0.885,"ALICE Simulation,pp #sqrt{#it{s}} = 13 TeV");
    letexTitle -> DrawLatex(0.225,0.815,"PYTHIA8 Monash Tune, Normalized to N_{eq. Pythia}");
    letexTitle -> DrawLatex(0.225,0.745,Info.Data());

    legend->Draw();

    return canvasSpectra;
}

TCanvas *allprint4hist(TH1D *h1,TH1D *h2,TH1D *h3,TH1D *h4){
    TH2D *h_grid=new TH2D("h_grid","",100,0,16.,100,h1->GetMinimum()/2,h1->GetMaximum()*25);
    h_grid->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h_grid->GetXaxis()->SetTitleOffset(1.35);
    h_grid->GetXaxis()->SetTitleSize(0.045);
    h_grid->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
    h_grid->GetYaxis()->SetTitleOffset(1.4);
    h_grid->GetYaxis()->SetTitleSize(0.045);


    TCanvas *canvasSpectra=new TCanvas ("canvasSpectra","canvasSpectra",850,950);
    canvasSpectra->SetName("canvasSpectra");
    canvasSpectra->cd();
    TLegend *legend=new TLegend(0.2,0.15,0.4,0.45);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);

    legend->AddEntry(h1,Form("%s (%0.0f)",h1->GetTitle(),h1->GetEntries()));
    legend->AddEntry(h2,Form("%s (%0.0f)",h2->GetTitle(),h2->GetEntries()));
    legend->AddEntry(h3,Form("%s (%0.0f)",h3->GetTitle(),h3->GetEntries()));
    legend->AddEntry(h4,Form("%s (%0.0f)",h4->GetTitle(),h4->GetEntries()));

    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.04);
    gPad->SetLeftMargin(0.17);
    gPad->SetBottomMargin(0.15);
    // gPad->SetGridx();
    // gPad->SetGridy();
    gPad->SetLogy();

    h_grid->Draw();
    h1->Draw("PE SAME");
    h2->Draw("PE SAME");
    h3->Draw("PE SAME");
    h4->Draw("PE SAME");

    legend->Draw();

    // TLatex *letexTitle = new TLatex();
    // letexTitle -> SetTextSize(0.05);
    // letexTitle -> SetNDC();
    // letexTitle -> SetTextFont(42);
    //
    // letexTitle -> DrawLatex(0.350,0.88,"ALICE Simulation,pp #sqrt{#it{s}} = 13 TeV");
    // letexTitle -> DrawLatex(0.29,0.81,"PYTHIA8 Monash Tune, N_{ev}=297900");
    // letexTitle -> DrawLatex(0.3854,0.74,"Reconstructed #mu#mu, Stat. Unc. Only");
    // letexTitle -> DrawLatex(0.35,0.67,"-4.0<y<-2.5");

    return canvasSpectra;
}

TCanvas *allprint5hist(TH1D *h1,TH1D *h2,TH1D *h3,TH1D *h4,TH1D *h5,double h_grid_minx,double h_grid_maxx,TString Info){

    TH2D* h_grid=new TH2D("h_grid","",100,h_grid_minx,h_grid_maxx,10,1e-11 ,h1->GetMaximum()*1500);
    h_grid->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.045);
    h_grid->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());
    h_grid->GetYaxis()->SetTitleOffset(1.3);
    h_grid->GetYaxis()->SetTitleSize(0.045);
    TCanvas *canvasSpectra=new TCanvas ("canvasSpectra","canvasSpectra",1200,1000);
    canvasSpectra->cd();
    TLegend *legend=new TLegend(0.45,0.49,0.7,0.69);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);

    legend->AddEntry(h1,Form("%s (%.5e)",h1->GetTitle(),h1->Integral(0.,35.0)));
    legend->AddEntry(h2,Form("%s (%.5e)",h2->GetTitle(),h2->Integral(0.,35.0)));
    legend->AddEntry(h3,Form("%s (%.5e)",h3->GetTitle(),h3->Integral(0.,35.0)));
    legend->AddEntry(h4,Form("%s (%.5e)",h4->GetTitle(),h4->Integral(0.,35.0)));
    legend->AddEntry(h5,Form("%s (%.5e)",h5->GetTitle(),h5->Integral(0.,35.0)));

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
    letexTitle -> SetTextSize(0.045);
    letexTitle -> SetNDC();
    letexTitle -> SetTextFont(42);

    h_grid->Draw();
    h1->Draw("PE SAME");
    h2->Draw("PE SAME");
    h3->Draw("PE SAME");
    h4->Draw("PE SAME");
    h5->Draw("PE SAME");

    letexTitle -> DrawLatex(0.225,0.885,"ALICE Simulation,pp #sqrt{#it{s}} = 13 TeV");
    letexTitle -> DrawLatex(0.225,0.815,"PYTHIA8 Monash Tune, Normalized to N_{eq. Pythia}");
    letexTitle -> DrawLatex(0.225,0.715,Info.Data());

    legend->Draw();

    return canvasSpectra;
}
