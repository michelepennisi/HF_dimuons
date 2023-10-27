#if !defined(__CINT__) || defined(__CLING__)
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TGrid.h"
#endif

void cs(TString opt = "LowMass_LowPt")
{
    if (opt.Contains("Cut_Yres"))
        printf("cross-section removing Y region");
    if (opt.Contains("LowMass_LowPt"))
        printf("cross-section removing 4<M<9 and pt<10");
    // TFile *_file0 = TFile::Open("~/dimuon_HF_pp/data/LHC18p/Hist_MC/3_11_22/HF/Tree_HF_MCDimuHFTree_merged.root");
    TFile *new_file = new TFile(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Hist_fromSim/Version1/HF/HistLite_HF_MCDimuHFTree_merged.root"), "READ");
    TH1D *histo_charm = (TH1D *)new_file->Get(Form("DiMuon/%s/Rec_DQ_cut_match_LT/ULS/h_PtYDiMu_%s_Rec_DQ_cut_match_LT_ULS_fromCharm", opt.Data(), opt.Data()));
    TH1D *histo_beauty = (TH1D *)new_file->Get(Form("DiMuon/%s/Rec_DQ_cut_match_LT/ULS/h_PtYDiMu_%s_Rec_DQ_cut_match_LT_ULS_fromBeauty", opt.Data(), opt.Data()));
    histo_charm->SetDirectory(0);
    histo_beauty->SetDirectory(0);

    new_file->Close();
    // TFile *_file0 = TFile::Open("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/Hist_pythia_sim_SoftQCD_Def_1000000_2710_DefaultBR_HFcount3.root");
    TFile *_file0 = TFile::Open("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/test/Hist_monash_sim_fixed_100Mev.root");
    //  TFile *_file0 = TFile::Open("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/test/Hist_pythia_sim_SoftQCD_Def_100000_3216_DefaultBR.root");

    // TFile *_file0 = TFile::Open("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/read_output/HFTreeHF_MCDimuHFTree_294925.root");
    // TH1D *hcharmquark = (TH1D *)_file0->Get("h_Ncharm_pairs");
    // TH1D *hcharmquark = (TH1D *)_file0->Get("h_Ncharm_pairs_v6");
    TH1D *hcharmquark = (TH1D *)_file0->Get("CS/h_Ncharm_pairs_v6");

    printf("Integral %0.4e\n", hcharmquark->Integral(2, 100));
    // TH1D *hbeautyquark = (TH1D *)_file0->Get("h_Nbeauty_pairs");
    // TH1D *hbeautyquark = (TH1D *)_file0->Get("h_Nbeauty_pairs_v6");
    TH1D *hbeautyquark = (TH1D *)_file0->Get("CS/h_Nbeauty_pairs_v6");

    Double_t tot_ccbar = 0;
    Double_t tot_bbbar = 0;

    for (Int_t i = 1; i < hcharmquark->GetNbinsX(); i++)
    {
        // printf("Bin %d || Content %0.0f\n", i - 1, hcharmquark->GetBinContent(i));
        tot_ccbar = tot_ccbar + (i - 1) * hcharmquark->GetBinContent(i);
    }

    for (Int_t i = 1; i < hbeautyquark->GetNbinsX(); i++)
    {
        // printf("Bin %d || Content %0.0f\n", i - 1, hbeautyquark->GetBinContent(i));
        tot_bbbar = tot_bbbar + (i - 1) * hbeautyquark->GetBinContent(i);
    }
    printf("Numero coppie ccbar %0.0f\n", tot_ccbar);
    printf("Numero coppie bbbar %0.0f\n", tot_bbbar);
    Double_t ccbar_cs_PYTHIA = (tot_ccbar / (1.5 * 2*1e+08)) * 56.42;
    // Double_t ccbar_cs_PYTHIA = (tot_ccbar / (1.5 * 2*9.88e+07)) * 56.42; //for mode2
    Double_t bbbar_cs_PYTHIA = (tot_bbbar / (1.5 * 2*1e+08)) * 56.42;
    // Double_t bbbar_cs_PYTHIA = (tot_bbbar / (1.5 * 2*9.88e+07)) * 56.42; //for mode2

    printf("Pythia cc bar cs at fwdy %0.3e in mb\n", ccbar_cs_PYTHIA);
    printf("Pythia bb bar cs at fwdy %0.3e in mb\n", bbbar_cs_PYTHIA);

    TFile *_file1 = TFile::Open("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Hist_fromSim/Version1/HF/Tree_HF_MCDimuHFTree_merged.root");
    TH1D *h_ev_PYTHIA = (TH1D *)_file1->Get("h_Nevents");
    Double_t n_ev_PYTHIA = 216 * h_ev_PYTHIA->GetBinContent(2);

    TFile *fIn_data = new TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/HistResults_merged.root", "READ");
    TH1D *fhNEv = (TH1D *)fIn_data->Get("fhNEv");

    printf("N ev PYTHIA MB from HF enriched sim%0.3e\n", n_ev_PYTHIA);
    printf("N ev DATA MB %0.3e\n", 2384.73 * fhNEv->GetBinContent(3));

    // Double_t c_frac_fit = 52287.2;
    Double_t c_frac_fit = 51153.9; //for low Mass low pt
    // Double_t c_frac_fit = 49051; // using mode2 tuning

    Double_t c_frac_fit_MB = c_frac_fit / (2384.73 * fhNEv->GetBinContent(3));

    printf("\n ------------------------------------- \n");

    printf("c number from fit %0.0f\n", c_frac_fit);
    printf("c number from fit norm MB %0.10f\n", c_frac_fit_MB);
    // Double_t c_frac_th = 26.7;
    printf("\n ------------------------------------- \n");

    // Double_t b_frac_fit = 19281;
    Double_t b_frac_fit = 15756.3; //for low Mass low pt
    // Double_t b_frac_fit = 17595.9; // using mode 2 tuning Low Mass Low pt

    Double_t b_frac_fit_MB = b_frac_fit / (2384.73 * fhNEv->GetBinContent(3));

    printf("b number from fit %0.0f\n", b_frac_fit);
    printf("b number from fit norm MB %0.10f\n", b_frac_fit_MB);

    printf("\n ------------------------------------- \n");

    Double_t c_frac_th = histo_charm->GetEntries();
    Double_t b_frac_th = histo_beauty->GetEntries();

    Double_t c_frac_th_MB = (c_frac_th) / (n_ev_PYTHIA);

    printf("c number from PYTHIA %0.0f\n", c_frac_th);
    printf("c number from PYTHIA norm MB %0.10f\n", c_frac_th_MB);

    Double_t b_frac_th_MB = (b_frac_th) / (n_ev_PYTHIA);
    // Double_t b_frac_th_MB = (b_frac_th) / (n_ev_PYTHIA)*0.94495413; //mode2 tune

    printf("b number from PYTHIA %0.0f\n", b_frac_th);
    printf("b number from PYTHIA norm MB %0.10f\n", b_frac_th_MB);
    printf("\n ------------------------------------- \n");
    // Double_t b_frac_th = 60.5;

    printf("tot_ccbar %0.7f || tot_bbbar %0.7f\n", tot_ccbar, tot_bbbar);
    printf("CC CS %0.7f || BB CS %0.7f\n", ccbar_cs_PYTHIA, bbbar_cs_PYTHIA);
    // printf("CC CS %0.7f || BB CS %0.7f\n", ccbar_cs_meas, bbbar_cs_meas);

    Double_t ccbar_cs_meas = (c_frac_fit_MB / c_frac_th_MB) * ccbar_cs_PYTHIA;

    Double_t bbbar_cs_meas = (b_frac_fit_MB / b_frac_th_MB) * bbbar_cs_PYTHIA;

    printf("\n ------------------------------------- \n");
    printf("CS cc bar measured %0.10f in mub\n", ccbar_cs_meas * 1000);
    printf("CS bb bar measured %0.10f in mub\n", bbbar_cs_meas * 1000);
}