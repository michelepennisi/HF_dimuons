void chisquare_data_scan()
{
    const Int_t n_scan = 5;

    Double_t bin_size[n_scan] = {0, 0.2, 0.4, 0.5, 1};

    Double_t pt_zero_ups[n_scan] = {494.5, 975.9, 1201.7, 2392};

    Double_t mass_zero_ups[n_scan] = {607.2, 1141.2, 1436.1, 2751};

    Double_t charm_number[n_scan] = {44174.7, 43815.6, 42953.5, 42255.4, 36257.4};

    Double_t beauty_number[n_scan] = {29737.4, 30105.2, 30968.1, 31667, 37676};

    Double_t charm_number_norm[n_scan];
    Double_t beauty_number_norm[n_scan];

    for (Int_t i = 0; i < n_scan; i++)
    {
        charm_number_norm[i] = charm_number[i] / 77094;
        beauty_number_norm[i] = beauty_number[i] / 77094;
    }

    TCanvas *c = new TCanvas("c", " ", 1000, 800);
    c->cd();
    // c->Divide(2, 1);
    c->cd();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    gStyle->SetOptStat(0);

    TGraph *pt_zero_ups_chi = new TGraph(n_scan, bin_size, pt_zero_ups);
    pt_zero_ups_chi->SetName("pt_zero_ups_chi");
    pt_zero_ups_chi->SetTitle(" ");
    pt_zero_ups_chi->SetMaximum(3220);
    pt_zero_ups_chi->SetMarkerStyle(20);
    pt_zero_ups_chi->SetMarkerSize(1.5);
    pt_zero_ups_chi->SetMarkerColor(kOrange - 2);
    pt_zero_ups_chi->SetLineColor(kOrange - 2);
    pt_zero_ups_chi->GetXaxis()->SetTitle("bin size");
    pt_zero_ups_chi->GetYaxis()->SetTitle("#chi^{2}");
    pt_zero_ups_chi->GetXaxis()->SetTitleOffset(1.3);
    pt_zero_ups_chi->GetXaxis()->SetTitleSize(0.0475);
    pt_zero_ups_chi->GetXaxis()->SetLabelSize(0.045);

    pt_zero_ups_chi->GetYaxis()->SetTitleOffset(1.3);
    pt_zero_ups_chi->GetYaxis()->SetNdivisions(505);
    pt_zero_ups_chi->GetYaxis()->SetTitleSize(0.0475);
    pt_zero_ups_chi->GetYaxis()->SetLabelSize(0.045);
    pt_zero_ups_chi->Draw();

    TGraph *mass_zero_ups_chi = new TGraph(n_scan, bin_size, mass_zero_ups);
    mass_zero_ups_chi->SetName("mass_zero_ups_chi");
    mass_zero_ups_chi->SetTitle(" ");

    mass_zero_ups_chi->SetMarkerStyle(47);
    mass_zero_ups_chi->SetMarkerSize(1.5);
    mass_zero_ups_chi->SetMarkerColor(kBlue - 2);
    mass_zero_ups_chi->SetLineColor(kBlue - 2);
    mass_zero_ups_chi->GetXaxis()->SetTitle("#it{m}_{#mu#mu} bin size (GeV/#it{c}^{2})");
    mass_zero_ups_chi->GetYaxis()->SetTitle("#chi^{2}");
    mass_zero_ups_chi->GetXaxis()->SetTitleOffset(1.3);
    mass_zero_ups_chi->GetXaxis()->SetTitleSize(0.0475);
    mass_zero_ups_chi->GetXaxis()->SetLabelSize(0.045);

    mass_zero_ups_chi->GetYaxis()->SetTitleOffset(1.3);
    mass_zero_ups_chi->GetYaxis()->SetNdivisions(505);
    mass_zero_ups_chi->GetYaxis()->SetTitleSize(0.0475);
    mass_zero_ups_chi->GetYaxis()->SetLabelSize(0.045);
    mass_zero_ups_chi->Draw("LPSAME");

    TLegend *legend1 = new TLegend(0.23, 0.515, 0.5, 0.715);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend1->SetFillStyle(0);
    legend1->SetLineColor(kWhite);
    legend1->SetBorderSize(0);
    legend1->SetTextSize(0.0425);
    // legend1->SetHeader("Toy MC");
    // legend1->SetTextAlign(11);
    legend1->AddEntry("pt_zero_ups_chi", "#it{p}_{T}", "LP");
    legend1->AddEntry("mass_zero_ups_chi", "#it{m}_{#mu#mu}", "LP");
    legend1->Draw();

    TCanvas *c1 = new TCanvas("c1", " ", 1000, 800);
    c1->cd();
    // c->Divide(2, 1);
    c1->cd();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    gStyle->SetOptStat(0);

    TH2D *h_grid = new TH2D("h_grid", " ", 5, -0.05, 1.05, 100, 0.231, 0.67);

    h_grid->GetXaxis()->SetTitle("bin size");
    h_grid->GetYaxis()->SetTitle("Fit output");
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.0475);
    h_grid->GetXaxis()->SetLabelSize(0.045);

    h_grid->GetYaxis()->SetTitleOffset(1.3);
    h_grid->GetYaxis()->SetNdivisions(505);
    h_grid->GetYaxis()->SetTitleSize(0.0475);
    h_grid->GetYaxis()->SetLabelSize(0.045);
    h_grid->GetYaxis()->SetMaxDigits(3);
    h_grid->Draw();

    TGraph *charm_output = new TGraph(n_scan, bin_size, charm_number_norm);
    charm_output->SetName("charm_output");
    charm_output->SetTitle(" ");

    charm_output->SetMinimum(0.36);
    charm_output->SetMarkerStyle(20);
    charm_output->SetMarkerSize(1.5);
    charm_output->SetMarkerColor(kMagenta - 2);
    charm_output->SetLineColor(kMagenta - 2);
    charm_output->GetXaxis()->SetTitle("bin size");
    charm_output->GetYaxis()->SetTitle("Fit output");
    charm_output->GetXaxis()->SetTitleOffset(1.3);
    charm_output->GetXaxis()->SetTitleSize(0.0475);
    charm_output->GetXaxis()->SetLabelSize(0.045);

    charm_output->GetYaxis()->SetTitleOffset(1.3);
    charm_output->GetYaxis()->SetNdivisions(505);
    charm_output->GetYaxis()->SetTitleSize(0.0475);
    charm_output->GetYaxis()->SetLabelSize(0.045);
    charm_output->Draw("LPSAME");

    TGraph *beauty_output = new TGraph(n_scan, bin_size, beauty_number_norm);
    beauty_output->SetName("beauty_output");
    beauty_output->SetTitle(" ");

    beauty_output->SetMarkerStyle(20);
    beauty_output->SetMarkerSize(1.5);
    beauty_output->SetMarkerColor(kSpring - 6);
    beauty_output->SetLineColor(kSpring - 6);
    beauty_output->GetXaxis()->SetTitle("bin size");
    beauty_output->GetYaxis()->SetTitle("Fit output");
    beauty_output->GetXaxis()->SetTitleOffset(1.3);
    beauty_output->GetXaxis()->SetTitleSize(0.0475);
    beauty_output->GetXaxis()->SetLabelSize(0.045);

    beauty_output->GetYaxis()->SetTitleOffset(1.3);
    beauty_output->GetYaxis()->SetNdivisions(505);
    beauty_output->GetYaxis()->SetTitleSize(0.0475);
    beauty_output->GetYaxis()->SetLabelSize(0.045);
    beauty_output->Draw("LPSAME");

    TLegend *legend2 = new TLegend(0.23, 0.515, 0.5, 0.715);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend2->SetFillStyle(0);
    legend2->SetLineColor(kWhite);
    legend2->SetBorderSize(0);
    legend2->SetTextSize(0.0425);
    // legend2->SetHeader("Toy MC");
    // legend2->SetTextAlign(11);
    legend2->AddEntry("charm_output", "charm fraction", "LP");
    legend2->AddEntry("beauty_output", "beauty fraction", "LP");
    legend2->Draw();
}