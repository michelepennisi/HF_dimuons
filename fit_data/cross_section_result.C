#include "/home/michele_pennisi/cernbox/HF_dimuons/common_include.h"

struct opt
{
    TString HF = "Charm";
    TString Generator = "powheg";

    //----MY measurament----//
    Double_t PYTHIA_ds_dy_fwd_LowMass_LowPt_noLF_HF[2] = {1582.03, 22.97};
    Double_t PYTHIA_stat_ds_dy_fwd_LowMass_LowPt_noLF_HF[2] = {1.43, 4.5};
    Double_t PYTHIA_syst_ds_dy_fwd_LowMass_LowPt_noLF_HF[2] = {14.89, 48.36};

    Double_t PYTHIA_ds_dy_fwd_LowMass_LowPt_withLF_HF[2] = {1335.77, 19.78};
    Double_t PYTHIA_stat_ds_dy_fwd_LowMass_LowPt_withLF_HF[2] = {1.721, 5.32};
    Double_t PYTHIA_syst_ds_dy_fwd_LowMass_LowPt_withLF_HF[2] = {15.17, 47.87};

    Double_t PYTHIA_ds_dy_fwd_FullMass_withLF_HF[2] = {1177.83, 25.38};
    Double_t PYTHIA_stat_ds_dy_fwd_FullMass_withLF_HF[2] = {1.63, 3.322};
    Double_t PYTHIA_syst_ds_dy_fwd_FullMass_withLF_HF[2] = {14.02, 28.98};

    Double_t POWHEG_ds_dy_fwd_LowMass_LowPt_noLF_HF[2] = {2009.21, 17.68};
    Double_t POWHEG_stat_ds_dy_fwd_LowMass_LowPt_noLF_HF[2] = {2.36, 4.5};
    Double_t POWHEG_syst_ds_dy_fwd_LowMass_LowPt_noLF_HF[2] = {8.44, 16.22};

    Double_t POWHEG_ds_dy_fwd_LowMass_LowPt_withLF_HF[2] = {1532.3, 18.996};
    Double_t POWHEG_stat_ds_dy_fwd_LowMass_LowPt_withLF_HF[2] = {3.24, 4.16};
    Double_t POWHEG_syst_ds_dy_fwd_LowMass_LowPt_withLF_HF[2] = {7.82, 10.06};

    Double_t POWHEG_ds_dy_fwd_FullMass_withLF_HF[2] = {1533.14, 16.13};
    Double_t POWHEG_stat_ds_dy_fwd_FullMass_withLF_HF[2] = {2.98, 4.68};
    Double_t POWHEG_syst_ds_dy_fwd_FullMass_withLF_HF[2] = {5.57, 8.73};

    //----ELECTRONS measurament----//

    Double_t PYTHIA_ds_dy_ALICE_ELECTRON[2] = {974., 79.};
    Double_t PYTHIA_stat_ds_dy_ALICE_ELECTRON[2] = {138., 14.};
    Double_t PYTHIA_syst_ds_dy_ALICE_ELECTRON[2] = {140., 11.};

    Double_t POWHEG_ds_dy_ALICE_ELECTRON[2] = {1417, 48.};
    Double_t POWHEG_stat_ds_dy_ALICE_ELECTRON[2] = {184., 14.};
    Double_t POWHEG_syst_ds_dy_ALICE_ELECTRON[2] = {204., 7.};

    //----HF measurament----//

    Double_t ds_dy_ALICE_HF[2] = {2021., 75.2};
    Double_t stat_ds_dy_ALICE_HF[2] = {140., 3.2};
    Double_t syst_ds_dy_ALICE_HF[2] = {150., 5.2};

    //----LHCb measurament----//

    Double_t ds_dy_Charm_LHCb[2] = {1086.3, 1086.3};
    Double_t stat_ds_dy_Charm_LHCb[2] = {27., 27.};
    Double_t syst_ds_dy_Charm_LHCb[2] = {282., 282.};

    Double_t ds_deta_Beauty_LHCb[12] = {29.2, 43.2, 54.6, 58.4, 57.4, 45.2, 45.2, 57.4, 58.4, 54.6, 43.2, 29.2};
    Double_t stat_ds_deta_Beauty_LHCb[12] = {1.6, 0.8, 1.2, 0.8, 1.0, 1.0, 1.6, 0.8, 1.2, 0.8, 1.0, 1.0};
    Double_t syst_ds_deta_Beauty_LHCb[12] = {4.8, 6.0, 5.8, 5.4, 4.4, 3.0, 4.8, 6.0, 5.8, 5.4, 4.4, 3.0};

    Bool_t Fullerror = kFALSE;
    TString Res_selector = "OnlyPrel";
    // TString Res_selector = "OnlyNew";
};

void cross_section_result()
{
    opt info;
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetStatBorderSize(1);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetDrawBorder(0);
    gStyle->SetTextFont(42);
    gStyle->SetStatFontSize(0.05);
    gStyle->SetStatX(0.97);
    gStyle->SetStatY(0.98);
    gStyle->SetStatH(0.03);
    gStyle->SetStatW(0.3);
    gStyle->SetTickLength(0.02, "y");
    gStyle->SetEndErrorSize(3);
    gStyle->SetLabelSize(0.04, "xyz");
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetLabelOffset(0.01, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetTitleOffset(0.9, "x");
    gStyle->SetTitleOffset(1.1, "y");
    gStyle->SetTitleSize(0.045, "xyz");
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetHatchesSpacing(0.5);

    Double_t y_fwd_ALICE[2] = {-3.25, 3.25};
    Double_t dy_fwd_ALICE[2] = {0.75, 0.75};

    Double_t y_fwd_LHCb[2] = {-2.75, 2.75};
    Double_t dy_fwd_LHCb[2] = {1.75, 1.75};

    Double_t eta_fwd_LHCb[12] = {-4.75, -4.25, -3.75, -3.25, -2.75, -2.25, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75};
    Double_t deta_fwd_LHCb[12] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};

    Double_t ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[2];
    Double_t stat_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[2];
    Double_t syst_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[2];

    Double_t ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[2];
    Double_t stat_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[2];
    Double_t syst_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[2];

    Double_t ds_dy_fwd_ALICE_FullMass_withLF_HF[2];
    Double_t stat_ds_dy_fwd_ALICE_FullMass_withLF_HF[2];
    Double_t syst_ds_dy_fwd_ALICE_FullMass_withLF_HF[2];

    // bb_cs_PYTHIA_fwd->GetAttLine(0)->SetLineColor(kRed);
    // bb_cs_PYTHIA_fwd->GetAttLine(0)->SetLineWidth(1);
    // bb_cs_PYTHIA_fwd->GetAttLine(1)->SetLineColor(kRed);
    // bb_cs_PYTHIA_fwd->GetAttLine(1)->SetLineWidth(1);
    // bb_cs_PYTHIA_fwd->GetAttFill(1)->SetFillStyle(0);

    Double_t y_mid[1] = {0};
    Double_t dy_mid[1] = {0.5};

    Double_t ds_dy_ALICE_ELECTRON[1];
    Double_t stat_ds_dy_ALICE_ELECTRON[1];
    Double_t syst_ds_dy_ALICE_ELECTRON[1];

    Double_t ds_dy_ALICE_HF[1];
    Double_t stat_ds_dy_ALICE_HF[1];
    Double_t syst_ds_dy_ALICE_HF[1];

    Double_t ds_dy_LHCb_CHARM[2] = {info.ds_dy_Charm_LHCb[0], info.ds_dy_Charm_LHCb[1]};
    Double_t stat_ds_dy_LHCb_CHARM[2] = {info.stat_ds_dy_Charm_LHCb[0], info.stat_ds_dy_Charm_LHCb[1]};
    Double_t syst_ds_dy_LHCb_CHARM[2] = {info.syst_ds_dy_Charm_LHCb[0], info.syst_ds_dy_Charm_LHCb[1]};

    Double_t ds_deta_LHCb_BEAUTY[12];
    Double_t stat_ds_deta_LHCb_BEAUTY[12];
    Double_t syst_ds_deta_LHCb_BEAUTY[12];
    for (Int_t i = 0; i < 12; i++)
    {
        ds_deta_LHCb_BEAUTY[i] = info.ds_deta_Beauty_LHCb[i];
        stat_ds_deta_LHCb_BEAUTY[i] = info.stat_ds_deta_Beauty_LHCb[i];
        syst_ds_deta_LHCb_BEAUTY[i] = info.syst_ds_deta_Beauty_LHCb[i];
    }

    TString FONLL_filename;
    Int_t HF_Selector = 999; // 0 for CHARM, 1 for BEAUTY

    if (info.HF.Contains("Charm"))
    {
        FONLL_filename.Form("FONLL_cc_Pt_0_30_CTEQ6.txt");
        HF_Selector = 0;
        ds_dy_ALICE_HF[0] = info.ds_dy_ALICE_HF[HF_Selector];
        stat_ds_dy_ALICE_HF[0] = info.stat_ds_dy_ALICE_HF[HF_Selector];
        syst_ds_dy_ALICE_HF[0] = info.syst_ds_dy_ALICE_HF[HF_Selector];
    }
    else if (info.HF.Contains("Beauty"))
    {
        FONLL_filename.Form("FONLL_bb_Pt_0_30_CTEQ6.txt");
        HF_Selector = 1;
        ds_dy_ALICE_HF[0] = info.ds_dy_ALICE_HF[HF_Selector];
        stat_ds_dy_ALICE_HF[0] = info.stat_ds_dy_ALICE_HF[HF_Selector];
        syst_ds_dy_ALICE_HF[0] = info.syst_ds_dy_ALICE_HF[HF_Selector];
    }

    if (info.Generator.Contains("pythia"))
    {
        ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] = info.PYTHIA_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector];
        ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] = info.PYTHIA_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector];
        stat_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] = ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] * info.PYTHIA_stat_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector] / 100.;
        stat_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] = ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] * info.PYTHIA_stat_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] = ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] * info.PYTHIA_syst_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] = ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] * info.PYTHIA_syst_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector] / 100.;

        ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] = info.PYTHIA_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector];
        ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] = info.PYTHIA_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector];
        stat_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] = ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] * info.PYTHIA_stat_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector] / 100.;
        stat_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] = ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] * info.PYTHIA_stat_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] = ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] * info.PYTHIA_syst_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] = ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] * info.PYTHIA_syst_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector] / 100.;

        ds_dy_fwd_ALICE_FullMass_withLF_HF[0] = info.PYTHIA_ds_dy_fwd_FullMass_withLF_HF[HF_Selector];
        ds_dy_fwd_ALICE_FullMass_withLF_HF[1] = info.PYTHIA_ds_dy_fwd_FullMass_withLF_HF[HF_Selector];
        stat_ds_dy_fwd_ALICE_FullMass_withLF_HF[0] = ds_dy_fwd_ALICE_FullMass_withLF_HF[0] * info.PYTHIA_stat_ds_dy_fwd_FullMass_withLF_HF[HF_Selector] / 100.;
        stat_ds_dy_fwd_ALICE_FullMass_withLF_HF[1] = ds_dy_fwd_ALICE_FullMass_withLF_HF[1] * info.PYTHIA_stat_ds_dy_fwd_FullMass_withLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_FullMass_withLF_HF[0] = ds_dy_fwd_ALICE_FullMass_withLF_HF[0] * info.PYTHIA_syst_ds_dy_fwd_FullMass_withLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_FullMass_withLF_HF[1] = ds_dy_fwd_ALICE_FullMass_withLF_HF[1] * info.PYTHIA_syst_ds_dy_fwd_FullMass_withLF_HF[HF_Selector] / 100.;

        ds_dy_ALICE_ELECTRON[0] = info.PYTHIA_ds_dy_ALICE_ELECTRON[HF_Selector];
        stat_ds_dy_ALICE_ELECTRON[0] = info.PYTHIA_stat_ds_dy_ALICE_ELECTRON[HF_Selector];
        syst_ds_dy_ALICE_ELECTRON[0] = info.PYTHIA_syst_ds_dy_ALICE_ELECTRON[HF_Selector];
    }
    else if (info.Generator.Contains("powheg"))
    {
        ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] = info.POWHEG_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector];
        ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] = info.POWHEG_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector];
        stat_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] = ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] * info.POWHEG_stat_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector] / 100.;
        stat_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] = ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] * info.POWHEG_stat_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] = ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[0] * info.POWHEG_syst_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] = ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF[1] * info.POWHEG_syst_ds_dy_fwd_LowMass_LowPt_noLF_HF[HF_Selector] / 100.;

        ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] = info.POWHEG_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector];
        ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] = info.POWHEG_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector];
        stat_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] = ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] * info.POWHEG_stat_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector] / 100.;
        stat_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] = ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] * info.POWHEG_stat_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] = ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[0] * info.POWHEG_syst_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] = ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF[1] * info.POWHEG_syst_ds_dy_fwd_LowMass_LowPt_withLF_HF[HF_Selector] / 100.;

        ds_dy_fwd_ALICE_FullMass_withLF_HF[0] = info.POWHEG_ds_dy_fwd_FullMass_withLF_HF[HF_Selector];
        ds_dy_fwd_ALICE_FullMass_withLF_HF[1] = info.POWHEG_ds_dy_fwd_FullMass_withLF_HF[HF_Selector];
        stat_ds_dy_fwd_ALICE_FullMass_withLF_HF[0] = ds_dy_fwd_ALICE_FullMass_withLF_HF[0] * info.POWHEG_stat_ds_dy_fwd_FullMass_withLF_HF[HF_Selector] / 100.;
        stat_ds_dy_fwd_ALICE_FullMass_withLF_HF[1] = ds_dy_fwd_ALICE_FullMass_withLF_HF[1] * info.POWHEG_stat_ds_dy_fwd_FullMass_withLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_FullMass_withLF_HF[0] = ds_dy_fwd_ALICE_FullMass_withLF_HF[0] * info.POWHEG_syst_ds_dy_fwd_FullMass_withLF_HF[HF_Selector] / 100.;
        syst_ds_dy_fwd_ALICE_FullMass_withLF_HF[1] = ds_dy_fwd_ALICE_FullMass_withLF_HF[1] * info.POWHEG_syst_ds_dy_fwd_FullMass_withLF_HF[HF_Selector] / 100.;

        ds_dy_ALICE_ELECTRON[0] = info.POWHEG_ds_dy_ALICE_ELECTRON[HF_Selector];
        stat_ds_dy_ALICE_ELECTRON[0] = info.POWHEG_stat_ds_dy_ALICE_ELECTRON[HF_Selector];
        syst_ds_dy_ALICE_ELECTRON[0] = info.POWHEG_syst_ds_dy_ALICE_ELECTRON[HF_Selector];
    }

    TGraphMultiErrors *cs_fwd_ALICE_LowMass_LowPt_noLH_HF = new TGraphMultiErrors("cs_fwd_ALICE_LowMass_LowPt_noLH_HF", "TGraphMultiErrors Example", 2, y_fwd_ALICE, ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF, dy_fwd_ALICE, dy_fwd_ALICE, stat_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF, stat_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->AddYError(2, syst_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF, syst_ds_dy_fwd_ALICE_LowMass_LowPt_noLF_HF);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->SetMarkerStyle(24);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->SetMarkerColor(kRed);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->SetLineColor(kRed);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->SetLineWidth(2);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->SetMarkerSize(2);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->GetAttLine(0)->SetLineColor(kRed);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->GetAttLine(0)->SetLineWidth(2);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->GetAttLine(1)->SetLineColor(kRed);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->GetAttLine(1)->SetLineWidth(2);
    cs_fwd_ALICE_LowMass_LowPt_noLH_HF->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors *cs_fwd_ALICE_LowMass_LowPt_withLH_HF = new TGraphMultiErrors("cs_fwd_ALICE_LowMass_LowPt_withLH_HF", "TGraphMultiErrors Example", 2, y_fwd_ALICE, ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF, dy_fwd_ALICE, dy_fwd_ALICE, stat_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF, stat_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->AddYError(2, syst_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF, syst_ds_dy_fwd_ALICE_LowMass_LowPt_withLF_HF);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->SetMarkerStyle(28);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->SetMarkerColor(kRed);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->SetLineColor(kRed);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->SetLineWidth(2);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->SetMarkerSize(2);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->GetAttLine(0)->SetLineColor(kRed);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->GetAttLine(0)->SetLineWidth(2);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->GetAttLine(1)->SetLineColor(kRed);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->GetAttLine(1)->SetLineWidth(2);
    cs_fwd_ALICE_LowMass_LowPt_withLH_HF->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors *cs_fwd_ALICE_FullMass_withLH_HF = new TGraphMultiErrors("cs_fwd_ALICE_FullMass_withLH_HF", "TGraphMultiErrors Example", 2, y_fwd_ALICE, ds_dy_fwd_ALICE_FullMass_withLF_HF, dy_fwd_ALICE, dy_fwd_ALICE, stat_ds_dy_fwd_ALICE_FullMass_withLF_HF, stat_ds_dy_fwd_ALICE_FullMass_withLF_HF);
    cs_fwd_ALICE_FullMass_withLH_HF->AddYError(2, syst_ds_dy_fwd_ALICE_FullMass_withLF_HF, syst_ds_dy_fwd_ALICE_FullMass_withLF_HF);
    cs_fwd_ALICE_FullMass_withLH_HF->SetMarkerStyle(46);
    cs_fwd_ALICE_FullMass_withLH_HF->SetMarkerColor(kRed);
    cs_fwd_ALICE_FullMass_withLH_HF->SetLineColor(kRed);
    cs_fwd_ALICE_FullMass_withLH_HF->SetLineWidth(2);
    cs_fwd_ALICE_FullMass_withLH_HF->SetMarkerSize(2);
    cs_fwd_ALICE_FullMass_withLH_HF->GetAttLine(0)->SetLineColor(kRed);
    cs_fwd_ALICE_FullMass_withLH_HF->GetAttLine(0)->SetLineWidth(2);
    cs_fwd_ALICE_FullMass_withLH_HF->GetAttLine(1)->SetLineColor(kRed);
    cs_fwd_ALICE_FullMass_withLH_HF->GetAttLine(1)->SetLineWidth(2);
    cs_fwd_ALICE_FullMass_withLH_HF->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors *cs_ALICE_ELECTRON = new TGraphMultiErrors("cs_ALICE_ELECTRON", "TGraphMultiErrors Example", 1, y_mid, ds_dy_ALICE_ELECTRON, dy_mid, dy_mid, stat_ds_dy_ALICE_ELECTRON, stat_ds_dy_ALICE_ELECTRON);
    cs_ALICE_ELECTRON->AddYError(1, syst_ds_dy_ALICE_ELECTRON, syst_ds_dy_ALICE_ELECTRON);
    cs_ALICE_ELECTRON->SetMarkerStyle(21);
    cs_ALICE_ELECTRON->SetMarkerColor(kMagenta + 2);
    cs_ALICE_ELECTRON->SetLineColor(kMagenta + 2);
    cs_ALICE_ELECTRON->SetLineWidth(2);
    cs_ALICE_ELECTRON->SetMarkerSize(2);
    cs_ALICE_ELECTRON->GetAttLine(0)->SetLineColor(kMagenta + 2);
    cs_ALICE_ELECTRON->GetAttLine(0)->SetLineWidth(2);
    cs_ALICE_ELECTRON->GetAttLine(1)->SetLineColor(kMagenta + 2);
    cs_ALICE_ELECTRON->GetAttLine(1)->SetLineWidth(2);
    cs_ALICE_ELECTRON->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors *cs_ALICE_HF = new TGraphMultiErrors("cs_ALICE_HF", "TGraphMultiErrors Example", 1, y_mid, ds_dy_ALICE_HF, dy_mid, dy_mid, stat_ds_dy_ALICE_HF, stat_ds_dy_ALICE_HF);
    cs_ALICE_HF->AddYError(1, syst_ds_dy_ALICE_ELECTRON, syst_ds_dy_ALICE_ELECTRON);
    cs_ALICE_HF->SetMarkerStyle(34);
    cs_ALICE_HF->SetMarkerColor(kAzure + 2);
    cs_ALICE_HF->SetLineColor(kAzure + 2);
    cs_ALICE_HF->SetLineWidth(2);
    cs_ALICE_HF->SetMarkerSize(2);
    cs_ALICE_HF->GetAttLine(0)->SetLineColor(kAzure + 2);
    cs_ALICE_HF->GetAttLine(0)->SetLineWidth(2);
    cs_ALICE_HF->GetAttLine(1)->SetLineColor(kAzure + 2);
    cs_ALICE_HF->GetAttLine(1)->SetLineWidth(2);
    cs_ALICE_HF->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors *cs_LHCb_CHARM = new TGraphMultiErrors("cs_LHCb_CHARM", "TGraphMultiErrors Example", 2, y_fwd_LHCb, ds_dy_LHCb_CHARM, dy_fwd_LHCb, dy_fwd_LHCb, stat_ds_dy_LHCb_CHARM, stat_ds_dy_LHCb_CHARM);
    cs_LHCb_CHARM->AddYError(2, syst_ds_dy_LHCb_CHARM, syst_ds_dy_LHCb_CHARM);
    cs_LHCb_CHARM->SetMarkerStyle(47);
    cs_LHCb_CHARM->SetMarkerColor(kGreen + 2);
    cs_LHCb_CHARM->SetLineColor(kGreen + 2);
    cs_LHCb_CHARM->SetLineWidth(2);
    cs_LHCb_CHARM->SetMarkerSize(2);
    cs_LHCb_CHARM->GetAttLine(0)->SetLineColor(kGreen + 2);
    cs_LHCb_CHARM->GetAttLine(0)->SetLineWidth(2);
    cs_LHCb_CHARM->GetAttLine(1)->SetLineColor(kGreen + 2);
    cs_LHCb_CHARM->GetAttLine(1)->SetLineWidth(2);
    cs_LHCb_CHARM->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors *cs_LHCb_BEAUTY = new TGraphMultiErrors("cs_LHCb_BEAUTY", "TGraphMultiErrors Example", 12, eta_fwd_LHCb, ds_deta_LHCb_BEAUTY, deta_fwd_LHCb, deta_fwd_LHCb, stat_ds_deta_LHCb_BEAUTY, stat_ds_deta_LHCb_BEAUTY);
    cs_LHCb_BEAUTY->AddYError(12, syst_ds_deta_LHCb_BEAUTY, syst_ds_deta_LHCb_BEAUTY);
    cs_LHCb_BEAUTY->SetMarkerStyle(47);
    cs_LHCb_BEAUTY->SetMarkerColor(kGreen + 2);
    cs_LHCb_BEAUTY->SetLineColor(kGreen + 2);
    cs_LHCb_BEAUTY->SetLineWidth(1);
    cs_LHCb_BEAUTY->SetMarkerSize(2);
    cs_LHCb_BEAUTY->GetAttLine(0)->SetLineColor(kGreen + 2);
    cs_LHCb_BEAUTY->GetAttLine(0)->SetLineWidth(1);
    cs_LHCb_BEAUTY->GetAttLine(1)->SetLineColor(kGreen + 2);
    cs_LHCb_BEAUTY->GetAttLine(1)->SetLineWidth(1);
    cs_LHCb_BEAUTY->GetAttFill(1)->SetFillStyle(0);

    vector<double> vecY, vecL_Y, vecH_Y;
    vector<double> vec_central, vec_min_central, vec_max_central;
    vector<double> low_central_error, high_central_error;

    vector<double> vec_min_mass, vec_max_mass;
    vector<double> low_mass_error, high_mass_error;

    vector<double> vec_min_scale, vec_max_scale;
    vector<double> low_scale_error, high_scale_error;

    vector<double> vec_min_pdf, vec_max_pdf;
    vector<double> low_pdf_error, high_pdf_error;

    double Y, L_Y, H_Y, central, min_central, max_central, min_scale, max_scale, min_mass, max_mass, min_pdf, max_pdf;

    ifstream inputFile(FONLL_filename.Data());

    while (inputFile >> Y >> L_Y >> H_Y >> central >> min_central >> max_central >> min_scale >> max_scale >> min_mass >> max_mass >> min_pdf >> max_pdf)
    {
        vecY.push_back(Y);
        vecL_Y.push_back(L_Y);
        vecH_Y.push_back(H_Y);

        vec_central.push_back(central);
        vec_min_central.push_back(min_central);
        vec_max_central.push_back(max_central);

        vec_min_scale.push_back(min_scale);
        vec_max_scale.push_back(max_scale);

        vec_min_mass.push_back(min_mass);
        vec_max_mass.push_back(max_mass);

        vec_min_pdf.push_back(min_pdf);
        vec_max_pdf.push_back(max_pdf);
    }

    for (size_t i(0); i < vec_central.size(); i++)
    {
        low_central_error.push_back(vec_central[i] - vec_min_central[i]);
        high_central_error.push_back(vec_max_central[i] - vec_central[i]);

        low_scale_error.push_back(vec_central[i] - vec_min_scale[i]);
        high_scale_error.push_back(vec_max_scale[i] - vec_central[i]);

        low_mass_error.push_back(vec_central[i] - vec_min_mass[i]);
        high_mass_error.push_back(vec_max_mass[i] - vec_central[i]);

        low_pdf_error.push_back(vec_central[i] - vec_min_pdf[i]);
        high_pdf_error.push_back(vec_max_pdf[i] - vec_central[i]);

        vec_central[i] = vec_central[i] * 1e-6;

        low_central_error[i] = low_central_error[i] * 1e-6;
        high_central_error[i] = high_central_error[i] * 1e-6;

        low_scale_error[i] = low_scale_error[i] * 1e-6;
        high_scale_error[i] = high_scale_error[i] * 1e-6;

        low_mass_error[i] = low_mass_error[i] * 1e-6;
        high_mass_error[i] = high_mass_error[i] * 1e-6;

        low_pdf_error[i] = low_pdf_error[i] * 1e-6;
        high_pdf_error[i] = high_pdf_error[i] * 1e-6;
    }

    TGraphAsymmErrors *FONLL_bb_cs_NNPDF_minmaxerror = new TGraphAsymmErrors(vec_central.size(), &vecY[0], &vec_central[0], &vecL_Y[0], &vecH_Y[0], &low_central_error[0], &high_central_error[0]);
    FONLL_bb_cs_NNPDF_minmaxerror->SetLineWidth(2);

    TGraphAsymmErrors *FONLL_bb_cs_NNPDF_scale_error = new TGraphAsymmErrors(vec_central.size(), &vecY[0], &vec_central[0], &vecL_Y[0], &vecH_Y[0], &low_scale_error[0], &high_scale_error[0]);
    FONLL_bb_cs_NNPDF_scale_error->SetLineWidth(2);

    TGraphAsymmErrors *FONLL_bb_cs_NNPDF_mass_error = new TGraphAsymmErrors(vec_central.size(), &vecY[0], &vec_central[0], &vecL_Y[0], &vecH_Y[0], &low_mass_error[0], &high_mass_error[0]);
    FONLL_bb_cs_NNPDF_mass_error->SetLineWidth(2);

    TGraphAsymmErrors *FONLL_bb_cs_NNPDF_pdf_error = new TGraphAsymmErrors(vec_central.size(), &vecY[0], &vec_central[0], &vecL_Y[0], &vecH_Y[0], &low_pdf_error[0], &high_pdf_error[0]);
    FONLL_bb_cs_NNPDF_pdf_error->SetLineWidth(2);

    if (info.Fullerror)
    {
        FONLL_bb_cs_NNPDF_minmaxerror->SetLineColorAlpha(kGray, 0.9);
        FONLL_bb_cs_NNPDF_minmaxerror->SetFillColorAlpha(kGray, 0.7);
        FONLL_bb_cs_NNPDF_scale_error->SetFillColorAlpha(kCyan - 10, 0.7);
        FONLL_bb_cs_NNPDF_scale_error->SetLineColorAlpha(kCyan - 10, 0.9);
        FONLL_bb_cs_NNPDF_mass_error->SetFillColorAlpha(kGreen, 0.7);
        FONLL_bb_cs_NNPDF_mass_error->SetLineColorAlpha(kGreen, 0.9);
        FONLL_bb_cs_NNPDF_pdf_error->SetFillColorAlpha(kMagenta - 9, 0.7);
        FONLL_bb_cs_NNPDF_pdf_error->SetLineColorAlpha(kMagenta - 9, 0.9);
    }
    else
    {
        FONLL_bb_cs_NNPDF_minmaxerror->SetLineColorAlpha(kOrange + 1, 0.9);
        FONLL_bb_cs_NNPDF_minmaxerror->SetFillColorAlpha(kOrange + 1, 0.7);
        FONLL_bb_cs_NNPDF_minmaxerror->SetLineWidth(4);
        FONLL_bb_cs_NNPDF_minmaxerror->SetFillStyle(3005);
        FONLL_bb_cs_NNPDF_scale_error->SetFillColorAlpha(kOrange + 1, 0.7);
        FONLL_bb_cs_NNPDF_scale_error->SetLineColorAlpha(kOrange + 1, 0.9);
        FONLL_bb_cs_NNPDF_scale_error->SetFillStyle(3005);
        FONLL_bb_cs_NNPDF_mass_error->SetFillColorAlpha(kOrange + 1, 0.7);
        FONLL_bb_cs_NNPDF_mass_error->SetLineColorAlpha(kOrange + 1, 0.9);
        FONLL_bb_cs_NNPDF_mass_error->SetFillStyle(3005);
        FONLL_bb_cs_NNPDF_pdf_error->SetFillColorAlpha(kOrange + 1, 0.7);
        FONLL_bb_cs_NNPDF_pdf_error->SetLineColorAlpha(kOrange + 1, 0.9);
        FONLL_bb_cs_NNPDF_pdf_error->SetFillStyle(3005);
    }

    TMultiGraph *bb_bar_cs_NNPDF = new TMultiGraph();
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_minmaxerror, "CX");
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_scale_error);
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_mass_error);
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_pdf_error);
    if (info.Res_selector.Contains("OnlyPrel"))
        bb_bar_cs_NNPDF->Add(cs_fwd_ALICE_LowMass_LowPt_noLH_HF, "APS; Z ; 5 s=0.5");
    else if (info.Res_selector.Contains("OnlyNew"))
    {
        bb_bar_cs_NNPDF->Add(cs_fwd_ALICE_LowMass_LowPt_withLH_HF, "APS; Z ; 5 s=0.5");
        bb_bar_cs_NNPDF->Add(cs_fwd_ALICE_FullMass_withLH_HF, "APS; Z ; 5 s=0.5");
    }
    else
    {
        bb_bar_cs_NNPDF->Add(cs_fwd_ALICE_LowMass_LowPt_noLH_HF, "APS; Z ; 5 s=0.5");
        bb_bar_cs_NNPDF->Add(cs_fwd_ALICE_LowMass_LowPt_withLH_HF, "APS; Z ; 5 s=0.5");
        bb_bar_cs_NNPDF->Add(cs_fwd_ALICE_FullMass_withLH_HF, "APS; Z ; 5 s=0.5");
    }
    bb_bar_cs_NNPDF->Add(cs_ALICE_ELECTRON, "APS; Z ; 5 s=0.5");
    bb_bar_cs_NNPDF->Add(cs_ALICE_HF, "APS; Z ; 5 s=0.5");
    if (info.HF.Contains("Charm"))
        bb_bar_cs_NNPDF->Add(cs_LHCb_CHARM, "APS; Z ; 5 s=0.5");
    else if (info.HF.Contains("Beauty"))
        bb_bar_cs_NNPDF->Add(cs_LHCb_BEAUTY, "APS; Z ; 5 s=0.5");

    if (info.HF.Contains("Beauty"))
    {
        bb_bar_cs_NNPDF->GetYaxis()->SetRangeUser(6.5e-1, 8e+2);
        bb_bar_cs_NNPDF->GetYaxis()->SetTitle("d#sigma_{b#bar{b}} / d#it{y} or d#sigma_{b#bar{b}} / d#it{eta} (#mub)");
    }
    else if (info.HF.Contains("Charm"))
    {
        bb_bar_cs_NNPDF->GetYaxis()->SetRangeUser(12.2e-0, 1.25e+4);
        bb_bar_cs_NNPDF->GetYaxis()->SetTitle("d#sigma_{c#bar{c}} / d#it{y} (#mub)");
    }

    bb_bar_cs_NNPDF->GetXaxis()->SetTitle("#it{y}");

    TCanvas *canvas = canvas_for_prel();
    canvas->cd();
    bb_bar_cs_NNPDF->Draw("A3");

    TLegend *Legend_bb_cs_NNPDF_FONLL;
    if (info.Fullerror)
    {
        Legend_bb_cs_NNPDF_FONLL = new TLegend(0.175, 0.18, 0.475, 0.525, " ", "brNDC");
        Legend_bb_cs_NNPDF_FONLL->SetHeader("       #it{FONLL}");
        Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_minmaxerror, "unc. tot.", "F");
        Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_scale_error, "unc. scale", "F");
        Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_mass_error, "unc. mass", "F");
        Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_pdf_error, "unc. pdf", "F");
    }
    else
    {
        Legend_bb_cs_NNPDF_FONLL = new TLegend(0.725, 0.18, 0.875, 0.325, " ", "brNDC");
        Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_minmaxerror, "FONLL", "F");
    }

    Legend_bb_cs_NNPDF_FONLL->SetBorderSize(0);
    Legend_bb_cs_NNPDF_FONLL->SetFillColor(10);
    Legend_bb_cs_NNPDF_FONLL->SetFillStyle(1);
    Legend_bb_cs_NNPDF_FONLL->SetLineStyle(0);
    Legend_bb_cs_NNPDF_FONLL->SetLineColor(0);
    Legend_bb_cs_NNPDF_FONLL->SetTextFont(42);
    Legend_bb_cs_NNPDF_FONLL->SetTextSize(0.04);
    Legend_bb_cs_NNPDF_FONLL->Draw("SAME");

    TLegend *Legend_bb_cs_NNPDF_Meas = new TLegend(0.135, 0.175, 0.8, 0.475, " ", "brNDC");
    if (info.Generator.Contains("pythia"))
    {
        if (info.Res_selector.Contains("OnlyPrel"))
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_LowMass_LowPt_noLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) PYTHIA8 fit no LF-HF}{4 < #it{m_{#mu#mu}} < 9 GeV/#it{c}^{2}, #it{p}_{T} < 10 GeV/#it{c}}", "EP");
        else if (info.Res_selector.Contains("OnlyNew"))
        {
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_LowMass_LowPt_withLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) PYTHIA8 fit with LF-HF}{4 < #it{m_{#mu#mu}} < 9 GeV/#it{c}^{2}, #it{p}_{T} < 10 GeV/#it{c}}", "EP");
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_FullMass_withLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) PYTHIA8 fit with LF-HF}{4 < #it{m_{#mu#mu}} < 30 GeV/#it{c}^{2}, #it{p}_{T} < 30 GeV/#it{c}}", "EP");
        }
        else
        {
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_LowMass_LowPt_noLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) PYTHIA8 fit no LF-HF}{4 < #it{m_{#mu#mu}} < 9 GeV/#it{c}^{2}, #it{p}_{T} < 10 GeV/#it{c}}", "EP");
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_LowMass_LowPt_withLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) PYTHIA8 fit with LF-HF}{4 < #it{m_{#mu#mu}} < 9 GeV/#it{c}^{2}, #it{p}_{T} < 10 GeV/#it{c}}", "EP");
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_FullMass_withLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) PYTHIA8 fit with LF-HF}{4 < #it{m_{#mu#mu}} < 30 GeV/#it{c}^{2}, #it{p}_{T} < 30 GeV/#it{c}}", "EP");
        }
        Legend_bb_cs_NNPDF_Meas->AddEntry(cs_ALICE_ELECTRON, "#splitline{(#it{m}_{e^{#plus}e^{#minus}}, #it{p}_{T, e^{#plus}e^{#minus}}) PYTHIA6 fit}{Phys. Lett. B788 (2019) 505}", "EP");
    }
    else if (info.Generator.Contains("powheg"))
    {
        if (info.Res_selector.Contains("OnlyPrel"))
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_LowMass_LowPt_noLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) POWHEG fit}{4 < #it{m_{#mu#mu} no LF-HF} < 9 GeV/#it{c}^{2}, #it{p}_{T} < 10 GeV/#it{c}}", "EP");
        else if (info.Res_selector.Contains("OnlyNew"))
        {
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_LowMass_LowPt_withLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) POWHEG fit}{4 < #it{m_{#mu#mu} with LF-HF} < 9 GeV/#it{c}^{2}, #it{p}_{T} < 10 GeV/#it{c}}", "EP");
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_FullMass_withLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) POWHEG fit}{4 < #it{m_{#mu#mu} with LF-HF} < 30 GeV/#it{c}^{2}, #it{p}_{T} < 30 GeV/#it{c}}", "EP");
        }
        else
        {
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_LowMass_LowPt_noLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) POWHEG fit}{4 < #it{m_{#mu#mu} no LF-HF} < 9 GeV/#it{c}^{2}, #it{p}_{T} < 10 GeV/#it{c}}", "EP");
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_LowMass_LowPt_withLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) POWHEG fit}{4 < #it{m_{#mu#mu} with LF-HF} < 9 GeV/#it{c}^{2}, #it{p}_{T} < 10 GeV/#it{c}}", "EP");
            Legend_bb_cs_NNPDF_Meas->AddEntry(cs_fwd_ALICE_FullMass_withLH_HF, "#splitline{(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) POWHEG fit}{4 < #it{m_{#mu#mu} with LF-HF} < 30 GeV/#it{c}^{2}, #it{p}_{T} < 30 GeV/#it{c}}", "EP");
        }

        Legend_bb_cs_NNPDF_Meas->AddEntry(cs_ALICE_ELECTRON, "#splitline{(#it{m}_{e^{#plus}e^{#minus}}, #it{p}_{T, e^{#plus}e^{#minus}}) POWHEG fit}{Phys. Lett. B788 (2019) 505}", "EP");
    }
    Legend_bb_cs_NNPDF_Meas->AddEntry(cs_ALICE_HF, "ALICE HF |#it{y}|<0.5", "EP");
    if (info.HF.Contains("Charm"))
        Legend_bb_cs_NNPDF_Meas->AddEntry(cs_LHCb_CHARM, "LHCb (est. from D^{0},D^{#plus},D^{#plus}_{s})", "EP");
    else if (info.HF.Contains("Beauty"))
        Legend_bb_cs_NNPDF_Meas->AddEntry(cs_LHCb_BEAUTY, "LHCb from B decays", "EP");
    Legend_bb_cs_NNPDF_Meas->SetBorderSize(0);
    Legend_bb_cs_NNPDF_Meas->SetFillColor(10);
    Legend_bb_cs_NNPDF_Meas->SetFillStyle(1);
    Legend_bb_cs_NNPDF_Meas->SetLineStyle(0);
    Legend_bb_cs_NNPDF_Meas->SetLineColor(0);
    Legend_bb_cs_NNPDF_Meas->SetTextFont(42);
    Legend_bb_cs_NNPDF_Meas->SetTextSize(0.02);
    Legend_bb_cs_NNPDF_Meas->SetNColumns(2);
    Legend_bb_cs_NNPDF_Meas->Draw("SAME");

    TLatex letexTitle;
    letexTitle.SetTextSize(0.055);
    letexTitle.SetNDC();
    letexTitle.SetTextFont(42);
    letexTitle.DrawLatex(0.2, 0.82, "ALICE Preliminary, pp#sqrt{#it{s}} = 13 TeV");
    letexTitle.SetTextSize(0.0375);
    letexTitle.DrawLatex(0.2, 0.74, "FONLL CTEQ6");
    canvas->SetName(Form("cs_%s_%s", info.HF.Data(), info.Generator.Data()));
    canvas->SetTitle(Form("cs_%s_%s", info.HF.Data(), info.Generator.Data()));
    canvas->SaveAs(Form("results/fit_result/cs_%s_%s_%s.pdf", info.HF.Data(), info.Generator.Data(),info.Res_selector.Data()));
}