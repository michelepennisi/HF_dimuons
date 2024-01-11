#include "/home/michele_pennisi/cernbox/common_include.h"

TCanvas *canvas_for_prel()
{
    gStyle->SetImageScaling(3.);
    int font = 42;

    TCanvas *canvas = new TCanvas("c", "c", 1600, 1200);
    canvas->cd();
    canvas->SetLogy();
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.10);
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(0);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.1518219);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameBorderMode(0);

    return canvas;
}

void cross_section_result(TString HF = "Beauty", TString Generator = "pythia")
{

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

    Double_t y_fwd[2] = {-3.25, 3.25};
    Double_t dy_fwd[2] = {0.75, 0.75};

    Double_t ds_bb_dy_PYTHIA_fwd[2];
    Double_t stat_ds_bb_dy_PYTHIA_fwd[2];
    Double_t syst_ds_bb_dy_PYTHIA_fwd[2];

    // bb_cs_PYTHIA_fwd->GetAttLine(0)->SetLineColor(kRed);
    // bb_cs_PYTHIA_fwd->GetAttLine(0)->SetLineWidth(1);
    // bb_cs_PYTHIA_fwd->GetAttLine(1)->SetLineColor(kRed);
    // bb_cs_PYTHIA_fwd->GetAttLine(1)->SetLineWidth(1);
    // bb_cs_PYTHIA_fwd->GetAttFill(1)->SetFillStyle(0);

    Double_t y_mid[1] = {0};
    Double_t dy_mid[1] = {0.5};

    Double_t ds_bb_dy_PYTHIA_mid[1];
    Double_t stat_ds_bb_dy_PYTHIA_mid[1];
    Double_t syst_ds_bb_dy_PYTHIA_mid[1];

    TString FONLL_filename;
    if (HF.Contains("Beauty"))
    {
        FONLL_filename.Form("FONLL_bb_Pt_0_30_CTEQ6.txt");
        if (Generator.Contains("pythia"))
        {
            ds_bb_dy_PYTHIA_fwd[0] = 22.4;
            ds_bb_dy_PYTHIA_fwd[1] = 22.4;
            stat_ds_bb_dy_PYTHIA_fwd[0] = 1.0;
            stat_ds_bb_dy_PYTHIA_fwd[1] = 1.0;
            syst_ds_bb_dy_PYTHIA_fwd[0] = 8.26;
            syst_ds_bb_dy_PYTHIA_fwd[1] = 8.26;

            ds_bb_dy_PYTHIA_mid[0] = 79.;
            stat_ds_bb_dy_PYTHIA_mid[0] = 14.;
            syst_ds_bb_dy_PYTHIA_mid[0] = 11.;
        }
        else if (Generator.Contains("powheg"))
        {
            ds_bb_dy_PYTHIA_fwd[0] = 14.93;
            ds_bb_dy_PYTHIA_fwd[1] = 14.93;
            stat_ds_bb_dy_PYTHIA_fwd[0] = 1.0;
            stat_ds_bb_dy_PYTHIA_fwd[1] = 1.0;
            syst_ds_bb_dy_PYTHIA_fwd[0] = 8.26;
            syst_ds_bb_dy_PYTHIA_fwd[1] = 8.26;

            ds_bb_dy_PYTHIA_mid[0] = 48.;
            stat_ds_bb_dy_PYTHIA_mid[0] = 14.;
            syst_ds_bb_dy_PYTHIA_mid[0] = 7.;
        }
    }
    else if (HF.Contains("Charm"))
    {
        FONLL_filename.Form("FONLL_cc_Pt_0_30_CTEQ6.txt");
        if (Generator.Contains("pythia"))
        {
            ds_bb_dy_PYTHIA_fwd[0] = 1540.;
            ds_bb_dy_PYTHIA_fwd[1] = 1540.;
            stat_ds_bb_dy_PYTHIA_fwd[0] = 21.2;
            stat_ds_bb_dy_PYTHIA_fwd[1] = 21.2;
            syst_ds_bb_dy_PYTHIA_fwd[0] = 203.4;
            syst_ds_bb_dy_PYTHIA_fwd[1] = 203.4;

            ds_bb_dy_PYTHIA_mid[0] = 974.;
            stat_ds_bb_dy_PYTHIA_mid[0] = 138.;
            syst_ds_bb_dy_PYTHIA_mid[0] = 140.;
        }
        else if (Generator.Contains("powheg"))
        {
            ds_bb_dy_PYTHIA_fwd[0] = 2258;
            ds_bb_dy_PYTHIA_fwd[1] = 2258;
            stat_ds_bb_dy_PYTHIA_fwd[0] = 1.0;
            stat_ds_bb_dy_PYTHIA_fwd[1] = 1.0;
            syst_ds_bb_dy_PYTHIA_fwd[0] = 8.26;
            syst_ds_bb_dy_PYTHIA_fwd[1] = 8.26;

            ds_bb_dy_PYTHIA_mid[0] = 1417.;
            stat_ds_bb_dy_PYTHIA_mid[0] = 184.;
            syst_ds_bb_dy_PYTHIA_mid[0] = 204.;
        }
    }
    TGraphMultiErrors *bb_cs_PYTHIA_fwd = new TGraphMultiErrors("bb_cs_PYTHIA_fwd", "TGraphMultiErrors Example", 2, y_fwd, ds_bb_dy_PYTHIA_fwd, dy_fwd, dy_fwd, stat_ds_bb_dy_PYTHIA_fwd, stat_ds_bb_dy_PYTHIA_fwd);
    bb_cs_PYTHIA_fwd->AddYError(2, syst_ds_bb_dy_PYTHIA_fwd, syst_ds_bb_dy_PYTHIA_fwd);
    bb_cs_PYTHIA_fwd->SetMarkerStyle(20);
    bb_cs_PYTHIA_fwd->SetMarkerColor(kRed);
    bb_cs_PYTHIA_fwd->SetLineColor(kRed);
    bb_cs_PYTHIA_fwd->SetLineWidth(1);
    bb_cs_PYTHIA_fwd->SetMarkerSize(2);
    bb_cs_PYTHIA_fwd->GetAttLine(0)->SetLineColor(kRed);
    bb_cs_PYTHIA_fwd->GetAttLine(0)->SetLineWidth(2);
    bb_cs_PYTHIA_fwd->GetAttLine(1)->SetLineColor(kRed);
    bb_cs_PYTHIA_fwd->GetAttLine(1)->SetLineWidth(2);
    bb_cs_PYTHIA_fwd->GetAttFill(1)->SetFillStyle(0);

    TGraphMultiErrors *bb_cs_PYTHIA_mid = new TGraphMultiErrors("bb_cs_PYTHIA_mid", "TGraphMultiErrors Example", 1, y_mid, ds_bb_dy_PYTHIA_mid, dy_mid, dy_mid, stat_ds_bb_dy_PYTHIA_mid, stat_ds_bb_dy_PYTHIA_mid);
    bb_cs_PYTHIA_mid->AddYError(1, syst_ds_bb_dy_PYTHIA_mid, syst_ds_bb_dy_PYTHIA_mid);
    bb_cs_PYTHIA_mid->SetMarkerStyle(21);
    bb_cs_PYTHIA_mid->SetMarkerColor(kMagenta + 2);
    bb_cs_PYTHIA_mid->SetLineColor(kMagenta + 2);
    bb_cs_PYTHIA_mid->SetLineWidth(1);
    bb_cs_PYTHIA_mid->SetMarkerSize(2);
    bb_cs_PYTHIA_mid->GetAttLine(0)->SetLineColor(kMagenta + 2);
    bb_cs_PYTHIA_mid->GetAttLine(0)->SetLineWidth(1);
    bb_cs_PYTHIA_mid->GetAttLine(1)->SetLineColor(kMagenta + 2);
    bb_cs_PYTHIA_mid->GetAttLine(1)->SetLineWidth(1);
    bb_cs_PYTHIA_mid->GetAttFill(1)->SetFillStyle(0);

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
    FONLL_bb_cs_NNPDF_minmaxerror->SetLineColorAlpha(kGray, 0.9);
    FONLL_bb_cs_NNPDF_minmaxerror->SetFillColorAlpha(kGray, 0.7);

    TGraphAsymmErrors *FONLL_bb_cs_NNPDF_scale_error = new TGraphAsymmErrors(vec_central.size(), &vecY[0], &vec_central[0], &vecL_Y[0], &vecH_Y[0], &low_scale_error[0], &high_scale_error[0]);
    FONLL_bb_cs_NNPDF_scale_error->SetLineWidth(2);
    FONLL_bb_cs_NNPDF_scale_error->SetFillColorAlpha(kCyan - 10, 0.7);
    FONLL_bb_cs_NNPDF_scale_error->SetLineColorAlpha(kCyan - 10, 0.9);

    TGraphAsymmErrors *FONLL_bb_cs_NNPDF_mass_error = new TGraphAsymmErrors(vec_central.size(), &vecY[0], &vec_central[0], &vecL_Y[0], &vecH_Y[0], &low_mass_error[0], &high_mass_error[0]);
    FONLL_bb_cs_NNPDF_mass_error->SetLineWidth(2);
    FONLL_bb_cs_NNPDF_mass_error->SetFillColorAlpha(kGreen, 0.7);
    FONLL_bb_cs_NNPDF_mass_error->SetLineColorAlpha(kGreen, 0.9);

    TGraphAsymmErrors *FONLL_bb_cs_NNPDF_pdf_error = new TGraphAsymmErrors(vec_central.size(), &vecY[0], &vec_central[0], &vecL_Y[0], &vecH_Y[0], &low_pdf_error[0], &high_pdf_error[0]);
    FONLL_bb_cs_NNPDF_pdf_error->SetLineWidth(2);
    FONLL_bb_cs_NNPDF_pdf_error->SetFillColorAlpha(kMagenta - 9, 0.7);
    FONLL_bb_cs_NNPDF_pdf_error->SetLineColorAlpha(kMagenta - 9, 0.9);

    TMultiGraph *bb_bar_cs_NNPDF = new TMultiGraph();
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_minmaxerror);
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_scale_error);
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_mass_error);
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_pdf_error);
    bb_bar_cs_NNPDF->Add(bb_cs_PYTHIA_fwd, "APS; Z ; 5 s=0.5");
    bb_bar_cs_NNPDF->Add(bb_cs_PYTHIA_mid, "APS; Z ; 5 s=0.5");
    
    if (HF.Contains("Beauty")){
        bb_bar_cs_NNPDF->GetYaxis()->SetRangeUser(0.08, 8e+2);
        bb_bar_cs_NNPDF->GetYaxis()->SetTitle("d#sigma_{b#bar{b}} / d#it{y} (#mub)");
    }
    else if (HF.Contains("Charm")){
        bb_bar_cs_NNPDF->GetYaxis()->SetRangeUser(1.2e-1, 1.2e+5);
        bb_bar_cs_NNPDF->GetYaxis()->SetTitle("d#sigma_{c#bar{c}} / d#it{y} (#mub)");
    }

    bb_bar_cs_NNPDF->GetXaxis()->SetTitle("#it{y}");

    TCanvas *canvas = canvas_for_prel();
    canvas->cd();
    bb_bar_cs_NNPDF->Draw("A3");

    TLegend *Legend_bb_cs_NNPDF_FONLL = new TLegend(0.175, 0.18, 0.475, 0.525, " ", "brNDC");
    Legend_bb_cs_NNPDF_FONLL->SetHeader("       #it{FONLL}");
    Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_minmaxerror, "unc. tot.", "F");
    Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_scale_error, "unc. scale", "F");
    Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_mass_error, "unc. mass", "F");
    Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_pdf_error, "unc. pdf", "F");

    Legend_bb_cs_NNPDF_FONLL->SetBorderSize(0);
    Legend_bb_cs_NNPDF_FONLL->SetFillColor(10);
    Legend_bb_cs_NNPDF_FONLL->SetFillStyle(1);
    Legend_bb_cs_NNPDF_FONLL->SetLineStyle(0);
    Legend_bb_cs_NNPDF_FONLL->SetLineColor(0);
    Legend_bb_cs_NNPDF_FONLL->SetTextFont(42);
    Legend_bb_cs_NNPDF_FONLL->SetTextSize(0.04);
    Legend_bb_cs_NNPDF_FONLL->Draw("SAME");

    TLegend *Legend_bb_cs_NNPDF_Meas = new TLegend(0.475, 0.18, 0.725, 0.525, " ", "brNDC");
    if (Generator.Contains("pythia"))
        Legend_bb_cs_NNPDF_Meas->AddEntry(bb_cs_PYTHIA_fwd, "(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) PYTHIA8 fit", "EP");
    else if (Generator.Contains("powheg"))
        Legend_bb_cs_NNPDF_Meas->AddEntry(bb_cs_PYTHIA_fwd, "(#it{m}_{#mu^{#plus}#mu^{#minus}}, #it{p}_{T, #mu^{#plus}#mu^{#minus}}) POWHEG fit", "EP");
    if (Generator.Contains("pythia"))
        Legend_bb_cs_NNPDF_Meas->AddEntry(bb_cs_PYTHIA_mid, "#splitline{(#it{m}_{e^{#plus}e^{#minus}}, #it{p}_{T, e^{#plus}e^{#minus}}) PYTHIA6 fit}{Phys. Lett. B788 (2019) 505}", "EP");
    else if (Generator.Contains("powheg"))
        Legend_bb_cs_NNPDF_Meas->AddEntry(bb_cs_PYTHIA_mid, "#splitline{(#it{m}_{e^{#plus}e^{#minus}}, #it{p}_{T, e^{#plus}e^{#minus}}) POWHEG fit}{Phys. Lett. B788 (2019) 505}", "EP");
    Legend_bb_cs_NNPDF_Meas->SetBorderSize(0);
    Legend_bb_cs_NNPDF_Meas->SetFillColor(10);
    Legend_bb_cs_NNPDF_Meas->SetFillStyle(1);
    Legend_bb_cs_NNPDF_Meas->SetLineStyle(0);
    Legend_bb_cs_NNPDF_Meas->SetLineColor(0);
    Legend_bb_cs_NNPDF_Meas->SetTextFont(42);
    Legend_bb_cs_NNPDF_Meas->SetTextSize(0.035);
    Legend_bb_cs_NNPDF_Meas->Draw("SAME");

    TLatex letexTitle;
    letexTitle.SetTextSize(0.055);
    letexTitle.SetNDC();
    letexTitle.SetTextFont(42);
    letexTitle.DrawLatex(0.2, 0.82, "ALICE Preliminary, pp#sqrt{#it{s}} = 13 TeV");
    letexTitle.SetTextSize(0.0375);
    letexTitle.DrawLatex(0.2, 0.74, "FONLL CTEQ6");
    canvas->SetName(Form("cs_%s_%s",HF.Data(),Generator.Data()));
    canvas->SetTitle(Form("cs_%s_%s",HF.Data(),Generator.Data()));
    canvas->SaveAs(Form("cs_%s_%s.pdf",HF.Data(),Generator.Data()));
}