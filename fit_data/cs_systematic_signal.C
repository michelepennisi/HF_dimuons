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

void cs_systematic_signal()
{

    TFile *_file0 = TFile::Open("~/dimuon_HF_pp/data/LHC18p/Hist_MC/3_11_22/HF/Tree_HF_MCDimuHFTree_merged.root");
    TH1D *hcharmquark = (TH1D *)_file0->Get("h_Ncharm_pairs");
    printf("Integral %0.2f\n", hcharmquark->Integral());
    TH1D *hbeautyquark = (TH1D *)_file0->Get("h_Nbeauty_pairs");
    // hbeautyquark->Draw();
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

    Double_t ccbar_cs_PYTHIA = (tot_ccbar / (1.5 * 216 * hcharmquark->GetEntries())) * 56.42;

    Double_t bbbar_cs_PYTHIA = (tot_bbbar / (1.5 * 216 * hbeautyquark->GetEntries())) * 56.42;

    // Double_t ccbar_cs_meas = (c_frac_fit / c_frac_th) * (ccbar_cs_PYTHIA/1.5);

    // Double_t bbbar_cs_meas = (b_frac_fit / b_frac_th) * bbbar_cs_PYTHIA;

    // Scaling_Factor[0] = 2384.73 * fhNEv->GetBinContent(3);

    // _file0->Close();

    TFile *fIn_data = new TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/HistResults_merged.root", "READ");
    TH1D *fhNEv = (TH1D *)fIn_data->Get("fhNEv");

    // Double_t N_ev = 73000;
    // Double_t N_ev_MB = N_ev / (2384.73 * fhNEv->GetBinContent(3));

    // printf("N_ev number from fit %0.0f\n", N_ev);
    // printf("N_ev number from fit norm MB %0.10f\n", N_ev_MB);

    Double_t c_frac_th = 16520;
    Double_t b_frac_th = 27089;

    Double_t c_frac_th_MB = c_frac_th / (216 * hcharmquark->GetEntries());

    printf("c number from PYTHIA %0.0f\n", c_frac_th);
    printf("c number from PYTHIA norm MB %0.10f\n", c_frac_th_MB);

    Double_t b_frac_th_MB = b_frac_th / (216 * hbeautyquark->GetEntries());

    printf("b number from PYTHIA %0.0f\n", b_frac_th);
    printf("b number from PYTHIA norm MB %0.10f\n", b_frac_th_MB);
    printf("\n ------------------------------------- \n");
    // Double_t b_frac_th = 60.5;

    printf("tot_ccbar %0.7f || tot_bbbar %0.7f\n", tot_ccbar, tot_bbbar);
    printf("CC CS %0.7f || BB CS %0.7f\n", ccbar_cs_PYTHIA, bbbar_cs_PYTHIA);
    // printf("CC CS %0.7f || BB CS %0.7f\n", ccbar_cs_meas, bbbar_cs_meas);

    const Int_t n_options = 4;
    Double_t c_frac_fit[n_options] = {51884, 52267.9, 46848.2, 58189.1};
    Double_t b_frac_fit[n_options] = {19478, 19078.8, 24492.8, 13145.2};

    Double_t c_frac_fit_MB[n_options];
    Double_t b_frac_fit_MB[n_options];

    Double_t ccbar_cs_meas[n_options];
    Double_t bbbar_cs_meas[n_options];
    TString name_options[n_options];
    name_options[0].Form("unbinned");
    name_options[1].Form("original");
    name_options[2].Form("scaled_Low2Up");
    name_options[3].Form("scaled_Up2Low");

    printf("N ev PYTHIA MB %0.3e\n", 216 * hcharmquark->GetEntries());
    printf("N ev DATA MB %0.3e\n", 2384.73 * fhNEv->GetBinContent(3));

    printf("\n ------------------------------------- \n");
    for (size_t opt = 0; opt < 4; opt++)
    {
        printf("\n ------------------------------------- \n");
        printf("\n ------------------||%s||------------------- \n",name_options[opt].Data());
        printf("\n ------------------------------------- \n");
        c_frac_fit_MB[opt] = c_frac_fit[opt] / (2384.73 * fhNEv->GetBinContent(3));
        printf("c number from fit %0.1f\n", c_frac_fit[opt]);
        printf("c number from fit norm MB %0.10f\n", c_frac_fit_MB[opt]);

        b_frac_fit_MB[opt] = b_frac_fit[opt] / (2384.73 * fhNEv->GetBinContent(3));
        

        printf("b number from fit %0.1f\n", b_frac_fit[opt]);
        printf("b number from fit norm MB %0.10f\n", b_frac_fit_MB[opt]);

        ccbar_cs_meas[opt] = (c_frac_fit_MB[opt] / c_frac_th_MB) * ccbar_cs_PYTHIA;

        bbbar_cs_meas[opt] = (b_frac_fit_MB[opt] / b_frac_th_MB) * bbbar_cs_PYTHIA;

        printf("\n ------------------------------------- \n");
        printf("CS cc bar measured %0.10f in mub\n", ccbar_cs_meas[opt] * 1000);
        printf("CS bb bar measured %0.10f in mub\n", bbbar_cs_meas[opt] * 1000);
        printf("\n ------------------------------------- \n");
    }

    // Double_t c_frac_th = 26.7;
}