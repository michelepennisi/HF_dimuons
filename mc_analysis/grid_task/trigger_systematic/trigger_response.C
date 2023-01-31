TCanvas *draw_response(TH1D *response_data, TH1D *response_MC, Double_t max_x);
double fit_function(double *x, double *par);
TCanvas *draw_response_withFit(TH1D *response_data, TH1D *response_MC, TF1 *fit_response_data, TF1 *fit_response_MC, Double_t max_x, TString title, Bool_t with_ratio);

void trigger_response()
{
    const Int_t eta_bins = 5;

    Double_t eta_binning[eta_bins + 1] = {-4.0, -3.7, -3.4, -3.1, -2.8, -2.5};

    // retrieve th2 pt-eta histograms from data
    TFile *fIn_data = new TFile("AnalysisResults_data_LHC18q.root", "READ");
    fIn_data->cd("PWG_MTRResponse");
    AliMergeableCollection *col_data = static_cast<AliMergeableCollection *>(fIn_data->FindObjectAny("MTRResponseOut"));

    TH2 *h_AllPt_data = static_cast<TH2 *>(col_data->GetObject("MB", "histoMatchAllPtPerEta"));
    h_AllPt_data->Sumw2();
    h_AllPt_data->Draw("COLZ");

    TH2 *h_LowPt_data = static_cast<TH2 *>(col_data->GetObject("MB", "histoMatchLowPtPerEta"));
    h_LowPt_data->Sumw2();
    h_LowPt_data->Draw("COLZ");

    // retrieve th2 pt-eta histograms from MC
    TFile *fIn_MC = new TFile("AnalysisResults_mc_strarlight_LHC18q.root", "READ");
    fIn_MC->cd("PWG_MTRResponse");

    AliMergeableCollection *col_MC = static_cast<AliMergeableCollection *>(fIn_MC->FindObjectAny("MTRResponseOut"));
    TH2 *h_AllPt_MC = static_cast<TH2 *>(col_MC->GetObject("MB", "histoMatchAllPtPerEta"));
    h_AllPt_MC->Sumw2();
    h_AllPt_MC->Draw("COLZ");

    TH2 *h_LowPt_MC = static_cast<TH2 *>(col_MC->GetObject("MB", "histoMatchLowPtPerEta"));
    h_LowPt_MC->Sumw2();
    h_LowPt_MC->Draw("COLZ");

    TH1D *h_PtAllPt_data[eta_bins];
    TH1D *h_PtLowPt_data[eta_bins];
    TH1D *response_data[eta_bins];
    // TGraphAsymmErrors *response_data[eta_bins];

    TF1 *fit_response_data[eta_bins];

    TH1D *h_PtAllPt_MC[eta_bins];
    TH1D *h_PtLowPt_MC[eta_bins];
    TH1D *response_MC[eta_bins];

    TF1 *fit_response_MC[eta_bins];

    TCanvas *canvas_response_withFit[eta_bins];

    TF1 *fit_response_data_shifted[eta_bins];
    TF1 *fit_response_MC_shifted[eta_bins];

    TH1D *response_data_shifted[eta_bins];
    TH1D *response_MC_shifted[eta_bins];

    TH1D *weight_shifted[eta_bins];

    TCanvas *canvas_check_shift_data[eta_bins];
    TCanvas *canvas_check_shift_MC[eta_bins];
    TCanvas *canvas_weight_shift[eta_bins];

    TFile *fOut = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/root_files/trigger_efficiency_weight.root", "UPDATE");

    for (size_t q = 0; q < eta_bins; q++)
    {
        h_PtAllPt_data[q] = (TH1D *)h_AllPt_data->ProjectionX("h_PtAllPt_data", q * 3, 3 + q * 3, "e ");
        // h_PtAllPt_data[q] = (TH1D *)h_AllPt_data->ProjectionX();
        h_PtAllPt_data[q]->Sumw2();

        h_PtLowPt_data[q] = (TH1D *)h_LowPt_data->ProjectionX("h_PtLowPt_data", q * 3, 3 + q * 3, "e ");
        // h_PtLowPt_data[q] = (TH1D *)h_LowPt_data->ProjectionX();
        h_PtLowPt_data[q]->Sumw2();

        response_data[q] = (TH1D *)h_PtLowPt_data[q]->Clone(Form("response_data_eta%zu", q));
        response_data[q]->SetTitle("Data");
        response_data[q]->Reset();
        response_data[q]->Divide(h_PtLowPt_data[q], h_PtAllPt_data[q], 1, 1, "B");
        response_data[q]->Sumw2();

        // response_data[q]->Sumw2();
        // response_data[q] = new TGraphAsymmErrors(h_PtLowPt_data[q], h_PtAllPt_data[q], " cpe0");
        // response_data[q]->Draw();
        // return;

        fit_response_data[q] = new TF1("fit_response_data", fit_function, 0.5, 30, 3);
        fit_response_data[q]->SetParameter(0, 1.2);
        fit_response_data[q]->SetParameter(1, 4.5);
        fit_response_data[q]->SetParameter(2, 1.0);
        fit_response_data[q]->SetTitle("Fit (Data)");
        response_data[q]->Fit(fit_response_data[q], "R0I");

        h_PtAllPt_MC[q] = (TH1D *)h_AllPt_MC->ProjectionX("h_PtAllPt_MC", q * 3, 3 + q * 3, "e");
        h_PtAllPt_MC[q]->Sumw2();

        h_PtLowPt_MC[q] = (TH1D *)h_LowPt_MC->ProjectionX("h_PtLowPt_MC", q * 3, 3 + q * 3, "e");
        h_PtLowPt_MC[q]->Sumw2();

        response_MC[q] = (TH1D *)h_PtLowPt_MC[q]->Clone(Form("response_MC_eta%zu", q));
        response_MC[q]->SetTitle("MC");
        // response_MC[q]->Sumw2();
        // response_MC[q]->Divide(h_PtAllPt_MC[q]);
        response_MC[q]->Reset();
        response_MC[q]->Divide(h_PtLowPt_MC[q], h_PtAllPt_MC[q], 1, 1, "B");
        response_MC[q]->Sumw2();

        fit_response_MC[q] = new TF1(Form("fit_response_MC_eta%zu", q), fit_function, 0.5, 30, 3);
        fit_response_MC[q]->SetParameter(0, 0.7);
        fit_response_MC[q]->SetParameter(1, 4);
        fit_response_MC[q]->SetParameter(2, 1.0);
        fit_response_MC[q]->SetTitle("Fit (MC)");
        response_MC[q]->Fit(fit_response_MC[q], "R0I");

        //         new TCanvas();
        // //         h_PtAllPt_data[q]->Draw("PE");
        // // h_PtLowPt_data[q]->Draw("SAMEPE");
        //         response_data[q]->Draw("PE");
        //         fit_response_data[q]->Draw("lsame");

        canvas_response_withFit[q] = draw_response_withFit(response_data[q], response_MC[q], fit_response_data[q], fit_response_MC[q], 30, Form("%0.1f < #eta_{#mu} %0.1f", eta_binning[q], eta_binning[q + 1]), kTRUE);
        canvas_response_withFit[q]->SetName(Form("canvas_response_withFit_eta%zu", q));
        canvas_response_withFit[q]->SetTitle();
        canvas_response_withFit[q]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/plot/%s.pdf", canvas_response_withFit[q]->GetTitle()));
        canvas_response_withFit[q]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/plot/%s.png", canvas_response_withFit[q]->GetTitle()));

        fit_response_data_shifted[q] = new TF1(Form("fit_response_data_shifted_eta%zu", q), fit_function, 0.65, 30, 3);
        fit_response_data_shifted[q]->FixParameter(0, fit_response_data[q]->GetParameter(0));
        fit_response_data_shifted[q]->FixParameter(1, fit_response_data[q]->GetParameter(1));
        fit_response_data_shifted[q]->FixParameter(2, fit_response_data[q]->GetParameter(2) - 0.5);
        fit_response_data_shifted[q]->SetTitle("Fit (Data shifted)");

        fit_response_MC_shifted[q] = new TF1(Form("fit_response_MC_shifted_eta%zu", q), fit_function, 0.65, 30, 3);
        fit_response_MC_shifted[q]->FixParameter(0, fit_response_MC[q]->GetParameter(0));
        fit_response_MC_shifted[q]->FixParameter(1, fit_response_MC[q]->GetParameter(1));
        fit_response_MC_shifted[q]->FixParameter(2, fit_response_MC[q]->GetParameter(2) - 0.5);
        fit_response_MC_shifted[q]->SetTitle("Fit (MC shifted)");

        // TAxis *xAxis = response_data[q]->GetXaxis();

        response_data_shifted[q] = (TH1D *)response_data[q]->Clone(Form("response_data_shifted_eta%zu", q));
        response_data_shifted[q]->Reset();
        response_data_shifted[q]->SetTitle("Data (shifted)");

        response_MC_shifted[q] = (TH1D *)response_MC[q]->Clone(Form("response_MC_shifted_eta%zu", q));
        response_MC_shifted[q]->Reset();
        response_MC_shifted[q]->SetTitle("MC (shifted)");

        for (size_t i = 0; i < response_data_shifted[q]->GetXaxis()->GetNbins(); i++)
        {
            response_data_shifted[q]->SetBinContent(i + 1, fit_response_data_shifted[q]->Eval(response_data_shifted[q]->GetXaxis()->GetBinCenter(i + 1), 0, 0));
            response_MC_shifted[q]->SetBinContent(i + 1, fit_response_MC_shifted[q]->Eval(response_data_shifted[q]->GetXaxis()->GetBinCenter(i + 1), 0, 0));
        }

        canvas_check_shift_data[q] = draw_response_withFit(response_data[q], response_data_shifted[q], fit_response_data[q], fit_response_data_shifted[q], 30, Form("%0.1f < #eta_{#mu} %0.1f", eta_binning[q], eta_binning[q + 1]), kFALSE);
        canvas_check_shift_data[q]->SetName(Form("canvas_check_shift_data_eta%zu", q));
        canvas_check_shift_data[q]->SetTitle(Form("canvas_check_shift_data_eta%zu", q));
        canvas_check_shift_data[q]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/plot/%s.pdf", canvas_check_shift_data[q]->GetTitle()));
        canvas_check_shift_data[q]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/plot/%s.png", canvas_check_shift_data[q]->GetTitle()));

        canvas_check_shift_MC[q] = draw_response_withFit(response_MC[q], response_MC_shifted[q], fit_response_MC[q], fit_response_MC_shifted[q], 30, Form("%0.1f < #eta_{#mu} %0.1f", eta_binning[q], eta_binning[q + 1]), kFALSE);
        canvas_check_shift_MC[q]->SetName(Form("canvas_check_shift_MC_eta%zu", q));
        canvas_check_shift_MC[q]->SetTitle(Form("canvas_check_shift_MC_eta%zu", q));
        canvas_check_shift_MC[q]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/plot/%s.pdf", canvas_check_shift_MC[q]->GetTitle()));
        canvas_check_shift_MC[q]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/plot/%s.png", canvas_check_shift_MC[q]->GetTitle()));

        canvas_weight_shift[q] = draw_response_withFit(response_data_shifted[q], response_MC_shifted[q], fit_response_data_shifted[q], fit_response_MC_shifted[q], 30, Form("%0.1f < #eta_{#mu} %0.1f", eta_binning[q], eta_binning[q + 1]), kTRUE);
        canvas_weight_shift[q]->SetName(Form("canvas_weight_shift_eta%zu", q));
        canvas_weight_shift[q]->SetTitle(Form("canvas_weight_shift_eta%zu", q));
        canvas_weight_shift[q]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/plot/%s.pdf", canvas_weight_shift[q]->GetTitle()));
        canvas_weight_shift[q]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/plot/%s.png", canvas_weight_shift[q]->GetTitle()));

        weight_shifted[q] = (TH1D *)response_data_shifted[q]->Clone(Form("weight_shifted_eta%zu", q));
        weight_shifted[q]->Divide(response_MC_shifted[q]);
        weight_shifted[q]->Sumw2();
        fOut->cd();
        weight_shifted[q]->Write(0, 2, 0);
        fOut->mkdir("plot");
        fOut->cd("plot");
        canvas_response_withFit[q]->Write(0, 2, 0);
        canvas_check_shift_data[q]->Write(0, 2, 0);
        canvas_check_shift_MC[q]->Write(0, 2, 0);
        canvas_weight_shift[q]->Write(0, 2, 0);
    }
    fOut->Close();
    /*
        printf("h_AllPt_data eta bin %d \n", h_AllPt_data->GetNbinsY());

        TH1D *eta = (TH1D *)h_AllPt_data->ProjectionY();
        eta->Draw();
        return;

        TH1D *h_PtAllPt_data = (TH1D *)h_AllPt_data->ProjectionX("h_PtAllPt_data", 0, 15, "e");
        h_PtAllPt_data->Sumw2();

        TH1D *h_PtLowPt_data = (TH1D *)h_LowPt_data->ProjectionX("h_PtLowPt_data", 0, 15, "e");
        h_PtLowPt_data->Sumw2();

        TH1D *response_eta1_data = (TH1D *)h_PtLowPt_data->Clone("response_eta1_data");
        response_eta1_data->SetTitle("Data");
        response_eta1_data->Sumw2();
        response_eta1_data->Divide(h_PtAllPt_data);
        //
        TH1D *h_PtAllPt_MC = (TH1D *)h_AllPt_MC->ProjectionX("h_PtAllPt_MC", 0, 15, "e");
        h_PtAllPt_MC->Sumw2();

        TH1D *h_PtLowPt_MC = (TH1D *)h_LowPt_MC->ProjectionX("h_PtLowPt_MC", 0, 15, "e");
        h_PtLowPt_MC->Sumw2();

        TH1D *response_eta1_MC = (TH1D *)h_PtLowPt_MC->Clone("response_eta1_MC");
        response_eta1_MC->Sumw2();
        response_eta1_MC->Divide(h_PtAllPt_MC);
        response_eta1_MC->SetTitle("MC");
        //
        TF1 *fit_response_data = new TF1("fit_response_data", fit_function, 0.5, 30, 4);
        fit_response_data->SetParameter(0, 2);
        fit_response_data->SetParameter(1, 4);
        fit_response_data->SetParameter(2, 1.0);
        fit_response_data->SetTitle("Fit (Data)");
        response_eta1_data->Fit(fit_response_data, "RL0I");

        TF1 *fit_response_MC = new TF1("fit_response_MC", fit_function, 0.5, 30, 4);
        fit_response_MC->SetParameter(0, 2);
        fit_response_MC->SetParameter(1, 4);
        fit_response_MC->SetParameter(2, 1.0);
        fit_response_MC->SetTitle("Fit (MC)");
        response_eta1_MC->Fit(fit_response_MC, "RL0I");

        // response_eta1_MC->Scale(1./response_eta1_MC->Integral());
        // TCanvas *canvas_response = new TCanvas("canvas_response", "canvas_response", 1000, 800);
        // canvas_response->cd();
        // response_eta1_data->SetLineColor(kRed);
        // response_eta1_data->SetMarkerColor(kRed);
        // response_eta1_data->Draw("PE");
        // response_eta1_MC->Draw("SAMEPE");
        new TCanvas();
        response_eta1_data->Draw();
        fit_response_data->Draw("SAME");
        TCanvas *canvas_response = draw_response(response_eta1_data, response_eta1_MC, 30);
        canvas_response->SetName("canvas_response");
        canvas_response->SetTitle("canvas_response");

        TCanvas *canvas_response_withFit = draw_response_withFit(response_eta1_data, response_eta1_MC, fit_response_data, fit_response_MC, 10);
        canvas_response_withFit->SetName("canvas_response_withFit");
        canvas_response_withFit->SetTitle("canvas_response_withFit");

        TH1D *response_eta1_data_zoom = (TH1D *)response_eta1_data->Clone("response_eta1_data_zoom");
        response_eta1_data_zoom->GetXaxis()->SetRangeUser(0, 10);

        TH1D *response_eta1_MC_zoom = (TH1D *)response_eta1_MC->Clone("response_eta1_MC_zoom");
        response_eta1_MC_zoom->GetXaxis()->SetRangeUser(0, 10);

        // TCanvas *canvas_response_zoom = draw_response(response_eta1_data_zoom, response_eta1_MC_zoom, 10);
        // canvas_response_zoom->SetName("canvas_response_zoom");
        // canvas_response_zoom->SetTitle("canvas_response_zoom");

        TF1 *fit_response_data_shifted = new TF1("fit_response_data_shifted", fit_function, 0.5, 30, 4);
        fit_response_data_shifted->FixParameter(0, fit_response_data->GetParameter(0));
        fit_response_data_shifted->FixParameter(1, fit_response_data->GetParameter(1));
        fit_response_data_shifted->FixParameter(2, fit_response_data->GetParameter(2) - 0.5);
        fit_response_data_shifted->SetTitle("Fit (Data shifted)");

        TCanvas *canvas_shift_tf1data = new TCanvas("canvas_shift_tf1data", "canvas_shift_tf1data");
        canvas_shift_tf1data->cd();

        fit_response_data_shifted->Draw("L");
        fit_response_data->Draw("LSAME");

        TF1 *fit_response_MC_shifted = new TF1("fit_response_MC_shifted", fit_function, 0.5, 30, 4);
        fit_response_MC_shifted->FixParameter(0, fit_response_MC->GetParameter(0));
        fit_response_MC_shifted->FixParameter(1, fit_response_MC->GetParameter(1));
        fit_response_MC_shifted->FixParameter(2, fit_response_MC->GetParameter(2) - 0.5);
        fit_response_MC_shifted->SetTitle("Fit (MC shifted)");

        TCanvas *canvas_shift_tf1_MC = new TCanvas("canvas_shift_tf1_MC", "canvas_shift_tf1_MC");
        canvas_shift_tf1_MC->cd();

        fit_response_MC_shifted->Draw("L");
        fit_response_MC->Draw("LSAME");

        TAxis *xAxis = response_eta1_data->GetXaxis();

        TH1D *response_eta1_data_shifted = (TH1D *)response_eta1_data->Clone("response_eta1_data_shifted");
        response_eta1_data_shifted->Reset();
        response_eta1_data_shifted->SetTitle("Data (shifted)");

        TH1D *response_eta1_MC_shifted = (TH1D *)response_eta1_MC->Clone("response_eta1_MC_shifted");
        response_eta1_MC_shifted->Reset();
        response_eta1_MC_shifted->SetTitle("MC (shifted)");

        for (size_t i = 0; i < xAxis->GetNbins(); i++)
        {
            response_eta1_data_shifted->SetBinContent(i + 1, fit_response_data_shifted->Eval(xAxis->GetBinCenter(i + 1), 0, 0));
            response_eta1_MC_shifted->SetBinContent(i + 1, fit_response_MC_shifted->Eval(xAxis->GetBinCenter(i + 1), 0, 0));
        }

        TCanvas *canvas_check_data_shift = draw_response_withFit(response_eta1_data, response_eta1_data_shifted, fit_response_data, fit_response_data_shifted, 30);
        canvas_check_data_shift->SetName("canvas_check_data_shift");
        canvas_check_data_shift->SetTitle("canvas_check_data_shift");
        TCanvas *canvas_check_MC_shift = draw_response_withFit(response_eta1_MC, response_eta1_MC_shifted, fit_response_MC, fit_response_MC_shifted, 30);
        canvas_check_MC_shift->SetName("canvas_check_MC_shift");
        canvas_check_MC_shift->SetTitle("canvas_check_MC_shift");

        TH1D *weight_shifted = (TH1D *)response_eta1_data_shifted->Clone("weight_shifted");
        weight_shifted->Sumw2();
        weight_shifted->Divide(response_eta1_MC_shifted);

        new TCanvas();
        weight_shifted->Draw();
        // return;

        TCanvas *canvas_weight_shifted = draw_response_withFit(response_eta1_data_shifted, response_eta1_MC_shifted, fit_response_data_shifted, fit_response_MC_shifted, 30);
        canvas_weight_shifted->SetName("canvas_weight_shifted");
        canvas_weight_shifted->SetTitle("canvas_weight_shifted");
        */
}

TCanvas *draw_response_withFit(TH1D *h_response_data, TH1D *h_response_MC, TF1 *h_fit_response_data, TF1 *h_fit_response_MC, Double_t max_x, TString title, Bool_t with_ratio)
{
    TH1D *response_data = (TH1D *)h_response_data->Clone("response_data");
    TH1D *response_MC = (TH1D *)h_response_MC->Clone("response_MC");
    TF1 *fit_response_data = (TF1 *)h_fit_response_data->Clone("response_data");
    TF1 *fit_response_MC = (TF1 *)h_fit_response_MC->Clone("response_MC");

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 800);
    canvas->SetTicks();
    canvas->cd();

    Double_t low_pad = 0;
    Double_t bottom_margin_pad = 0.15;
    Double_t factor = 0.8;
    if (with_ratio)
    {
        low_pad = 0.375;
        bottom_margin_pad = 0.0;
        factor = 1.0;
    }

    TPad *pad1 = new TPad("pad1", "pad1", 0, low_pad, 1.0, 1.0);
    pad1->SetTicks();
    // pad1->SetLogy(1);
    pad1->SetTopMargin(0.05);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.14);
    pad1->SetBottomMargin(bottom_margin_pad);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, low_pad);
    pad2->SetTicks();
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.03);
    pad2->SetLeftMargin(0.14);
    pad2->SetBottomMargin(0.25);

    if (with_ratio)
        pad2->Draw();

    TH2D *h_Grid = new TH2D("h_Grid", Form("%s", title.Data()), max_x, 0, max_x, 100, -0.1, 1.4);

    h_Grid->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    h_Grid->GetYaxis()->SetTitle("Low_{#it{p}_{T}}/All_{#it{p}_{T}}");
    h_Grid->GetXaxis()->SetTitleOffset(1.3 * factor);
    h_Grid->GetXaxis()->SetTitleSize(0.0475 * factor);
    h_Grid->GetXaxis()->SetLabelSize(0.045 * factor);
    h_Grid->GetYaxis()->SetNdivisions(505);
    h_Grid->GetYaxis()->SetTitleOffset(0.9 * factor);
    h_Grid->GetYaxis()->SetTitleSize(0.065 * factor);
    h_Grid->GetYaxis()->SetLabelSize(0.055 * factor);

    response_data->SetMarkerStyle(20);
    response_data->SetMarkerSize(1.1);
    response_data->SetMarkerColor(kRed);
    response_data->SetLineColor(kRed);
    response_data->SetLineWidth(2.0);
    response_data->Sumw2();

    fit_response_data->SetLineColor(kRed);
    fit_response_data->SetLineStyle(7);
    fit_response_data->SetLineWidth(3.0);

    response_MC->SetMarkerStyle(20);
    response_MC->SetMarkerSize(1.1);
    response_MC->SetMarkerColor(kBlue);
    response_MC->SetLineColor(kBlue);
    response_MC->SetLineWidth(2.0);
    response_MC->Sumw2();

    fit_response_MC->SetLineColor(kBlue);
    fit_response_MC->SetLineStyle(7);
    fit_response_MC->SetLineWidth(3.0);

    TH1D *ratio = (TH1D *)response_data->Clone("ratio");
    ratio->Divide(response_MC);
    ratio->SetTitle(" ");

    pad1->cd();

    h_Grid->Draw();
    response_data->Draw("PESAME");
    response_MC->Draw("PESAME");

    fit_response_data->Draw("SAME");
    fit_response_MC->Draw("SAME");

    TLegend *legend = new TLegend(0.315, 0.35, 0.85, 0.5);
    legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0525 * factor);
    // legend->SetHeader("MC");
    legend->SetTextAlign(12);
    legend->AddEntry(response_data, Form("%s", response_data->GetTitle()), "PE");
    legend->AddEntry(response_MC, Form("%s", response_MC->GetTitle()), "PE");
    legend->AddEntry(fit_response_data, Form("%s", fit_response_data->GetTitle()), "L");
    legend->AddEntry(fit_response_MC, Form("%s", fit_response_MC->GetTitle()), "L");
    TLatex *letexTitle = new TLatex();
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);
    letexTitle->SetTextSize(0.0525 * factor);
    letexTitle->DrawLatex(0.165, 0.85, Form("LHC18q, %s", title.Data()));
    letexTitle->SetTextSize(0.0425 * factor);
    letexTitle->DrawLatex(0.325, 0.252, Form("#color[2]{b = %0.3f #pm %0.3f}", fit_response_data->GetParameter(1), fit_response_data->GetParError(1)));
    letexTitle->DrawLatex(0.325, 0.185, Form("#color[2]{x_{0} = %0.3f #pm %0.3f}", fit_response_data->GetParameter(2), fit_response_data->GetParError(2)));

    letexTitle->DrawLatex(0.60, 0.252, Form("#color[4]{b = %0.3f #pm %0.3f}", fit_response_MC->GetParameter(1), fit_response_MC->GetParError(1)));
    letexTitle->DrawLatex(0.60, 0.185, Form("#color[4]{x_{0} = %0.3f #pm %0.3f}", fit_response_MC->GetParameter(2), fit_response_MC->GetParError(2)));

    legend->Draw();

    if (with_ratio)
    {
        pad2->cd();

        TH2D *h_grid = new TH2D("h_grid", " ", max_x, 0, 10, max_x, 0.75, 1.35);
        h_grid->GetYaxis()->SetTitle(Form("Ratio"));
        h_grid->GetYaxis()->CenterTitle();
        h_grid->GetYaxis()->SetNdivisions(504);
        h_grid->GetYaxis()->SetTitleSize(0.08);
        h_grid->GetYaxis()->SetTitleOffset(0.8);
        h_grid->GetYaxis()->SetLabelOffset(0.02);
        h_grid->GetYaxis()->SetLabelSize(0.1);

        h_grid->GetXaxis()->SetTitleSize(0.1);
        h_grid->GetXaxis()->SetTitleOffset(1.1);
        h_grid->GetXaxis()->SetLabelSize(0.1);
        h_grid->GetXaxis()->SetTitle(ratio->GetXaxis()->GetTitle());
        h_grid->Draw();

        ratio->SetMarkerColor(kBlack);
        ratio->SetLineColor(kBlack);
        ratio->Draw("PESAME");
    }

    return canvas;
}

double fit_function(double *x, double *par)
{
    return par[0] / (1 + TMath::Exp(-par[1] * (x[0] - par[2])));
}