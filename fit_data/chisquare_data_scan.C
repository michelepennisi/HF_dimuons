void chisquare_data_scan()
{
    const Int_t n_scan = 5;

    Double_t bin_size[n_scan] = {0.1, 0.2, 0.5, 1, 1.5};

    Double_t charm_number_ROOFIT[n_scan] = {44066, 43819, 42260, 36256, 24154};
    Double_t beauty_number_ROOFIT[n_scan] = {29854, 30102, 31662, 37673, 49788};

    Double_t charm_number_ROOFIT_norm[n_scan];
    Double_t beauty_number_ROOFIT_norm[n_scan];

    Double_t PT_charm_number_ROOT_TEST[n_scan] = {40599., 40494., 39753., 37332., 32969.};
    Double_t PT_beauty_number_ROOT_TEST[n_scan] = {33360., 33447., 34066., 36050., 39680.};

    Double_t PT_charm_number_ROOT_TEST_IOPTION[n_scan] = {40630., 40620., 40538., 40493., 40118.};
    Double_t PT_beauty_number_ROOT_TEST_IOPTION[n_scan] = {33334., 33345., 33426., 33472., 33848.};

    Double_t PT_charm_number_ROOT_TEST_norm[n_scan];
    Double_t PT_beauty_number_ROOT_TEST_norm[n_scan];

    Double_t PT_charm_number_ROOT_TEST_IOPTION_norm[n_scan];
    Double_t PT_beauty_number_ROOT_TEST_IOPTION_norm[n_scan];

    Double_t MASS_charm_number_ROOT_TEST[n_scan] = {46847., 46855, 47272, 48295, 51037};
    Double_t MASS_beauty_number_ROOT_TEST[n_scan] = {27077., 27128., 27134., 27618., 27433.};

    Double_t MASS_charm_number_ROOT_TEST_IOPTION[n_scan] = {46820., 46751, 46616, 45670, 45043};
    Double_t MASS_beauty_number_ROOT_TEST_IOPTION[n_scan] = {27083., 27153., 27288., 28233., 28861.};

    Double_t MASS_charm_number_ROOT_TEST_norm[n_scan];
    Double_t MASS_beauty_number_ROOT_TEST_norm[n_scan];

    Double_t MASS_charm_number_ROOT_TEST_IOPTION_norm[n_scan];
    Double_t MASS_beauty_number_ROOT_TEST_IOPTION_norm[n_scan];

    for (Int_t i = 0; i < n_scan; i++)
    {
        charm_number_ROOFIT_norm[i] = charm_number_ROOFIT[i] / 77094;
        beauty_number_ROOFIT_norm[i] = beauty_number_ROOFIT[i] / 77094;

        PT_charm_number_ROOT_TEST_norm[i] = PT_charm_number_ROOT_TEST[i] / 77094;
        PT_beauty_number_ROOT_TEST_norm[i] = PT_beauty_number_ROOT_TEST[i] / 77094;

        MASS_charm_number_ROOT_TEST_norm[i] = MASS_charm_number_ROOT_TEST[i] / 77094;
        MASS_beauty_number_ROOT_TEST_norm[i] = MASS_beauty_number_ROOT_TEST[i] / 77094;

        PT_charm_number_ROOT_TEST_IOPTION_norm[i] = PT_charm_number_ROOT_TEST_IOPTION[i] / 77094;
        PT_beauty_number_ROOT_TEST_IOPTION_norm[i] = PT_beauty_number_ROOT_TEST_IOPTION[i] / 77094;

        MASS_charm_number_ROOT_TEST_IOPTION_norm[i] = MASS_charm_number_ROOT_TEST_IOPTION[i] / 77094;
        MASS_beauty_number_ROOT_TEST_IOPTION_norm[i] = MASS_beauty_number_ROOT_TEST_IOPTION[i] / 77094;
    }

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

    TGraph *charm_output_ROOFIT = new TGraph(n_scan, bin_size, charm_number_ROOFIT_norm);
    charm_output_ROOFIT->SetName("charm_output_ROOFIT");
    charm_output_ROOFIT->SetTitle(" ");

    charm_output_ROOFIT->SetMinimum(0.36);
    charm_output_ROOFIT->SetMarkerStyle(24);
    charm_output_ROOFIT->SetMarkerSize(1.5);
    charm_output_ROOFIT->SetMarkerColor(kMagenta - 2);
    charm_output_ROOFIT->SetLineColor(kMagenta - 2);

    TGraph *beauty_output_ROOFIT = new TGraph(n_scan, bin_size, beauty_number_ROOFIT_norm);
    beauty_output_ROOFIT->SetName("beauty_output_ROOFIT");
    beauty_output_ROOFIT->SetTitle(" ");

    beauty_output_ROOFIT->SetMarkerStyle(24);
    beauty_output_ROOFIT->SetMarkerSize(1.5);
    beauty_output_ROOFIT->SetMarkerColor(kSpring - 6);
    beauty_output_ROOFIT->SetLineColor(kSpring - 6);

    TH2D *PT_h_grid = new TH2D("PT_h_grid", " ", 5, 0.05, 1.55, 100, 0.231, 0.67);

    PT_h_grid->GetXaxis()->SetTitle("#it{p}_{T} bin size (GeV/#it{c})");
    PT_h_grid->GetYaxis()->SetTitle("Fit output");
    PT_h_grid->GetXaxis()->SetTitleOffset(1.3);
    PT_h_grid->GetXaxis()->SetTitleSize(0.0475);
    PT_h_grid->GetXaxis()->SetLabelSize(0.045);

    PT_h_grid->GetYaxis()->SetTitleOffset(1.3);
    PT_h_grid->GetYaxis()->SetNdivisions(505);
    PT_h_grid->GetYaxis()->SetTitleSize(0.0475);
    PT_h_grid->GetYaxis()->SetLabelSize(0.045);
    PT_h_grid->GetYaxis()->SetMaxDigits(3);
    PT_h_grid->Draw();

    TGraph *pt_charm_output_ROOT_TEST = new TGraph(n_scan, bin_size, PT_charm_number_ROOT_TEST_norm);
    pt_charm_output_ROOT_TEST->SetName("pt_charm_output_ROOT_TEST");
    pt_charm_output_ROOT_TEST->SetTitle(" ");

    pt_charm_output_ROOT_TEST->SetMinimum(0.36);
    pt_charm_output_ROOT_TEST->SetMarkerStyle(20);
    pt_charm_output_ROOT_TEST->SetMarkerSize(1.5);
    pt_charm_output_ROOT_TEST->SetMarkerColor(kMagenta - 2);
    pt_charm_output_ROOT_TEST->SetLineColor(kMagenta - 2);

    TGraph *pt_charm_output_ROOT_TEST_IOPTION = new TGraph(n_scan, bin_size, PT_charm_number_ROOT_TEST_IOPTION_norm);
    pt_charm_output_ROOT_TEST_IOPTION->SetName("pt_charm_output_ROOT_TEST_IOPTION");
    pt_charm_output_ROOT_TEST_IOPTION->SetTitle(" ");

    pt_charm_output_ROOT_TEST_IOPTION->SetMinimum(0.36);
    pt_charm_output_ROOT_TEST_IOPTION->SetMarkerStyle(24);
    pt_charm_output_ROOT_TEST_IOPTION->SetMarkerSize(2.5);
    pt_charm_output_ROOT_TEST_IOPTION->SetMarkerColor(kMagenta - 2);
    pt_charm_output_ROOT_TEST_IOPTION->SetLineColor(kMagenta - 2);

    pt_charm_output_ROOT_TEST->Draw("LPSAME");
    pt_charm_output_ROOT_TEST_IOPTION->Draw("LPSAME");
    // charm_output_ROOFIT->Draw("LPSAME");

    TGraph *pt_beauty_output_ROOT_TEST = new TGraph(n_scan, bin_size, PT_beauty_number_ROOT_TEST_norm);
    pt_beauty_output_ROOT_TEST->SetName("pt_beauty_output_ROOT_TEST");
    pt_beauty_output_ROOT_TEST->SetTitle(" ");

    pt_beauty_output_ROOT_TEST->SetMinimum(0.36);
    pt_beauty_output_ROOT_TEST->SetMarkerStyle(20);
    pt_beauty_output_ROOT_TEST->SetMarkerSize(1.5);
    pt_beauty_output_ROOT_TEST->SetMarkerColor(kSpring - 6);
    pt_beauty_output_ROOT_TEST->SetLineColor(kSpring - 6);

    TGraph *pt_beauty_output_ROOT_TEST_IOPTION = new TGraph(n_scan, bin_size, PT_beauty_number_ROOT_TEST_IOPTION_norm);
    pt_beauty_output_ROOT_TEST_IOPTION->SetName("pt_beauty_output_ROOT_TEST_IOPTION");
    pt_beauty_output_ROOT_TEST_IOPTION->SetTitle(" ");

    pt_beauty_output_ROOT_TEST_IOPTION->SetMinimum(0.36);
    pt_beauty_output_ROOT_TEST_IOPTION->SetMarkerStyle(24);
    pt_beauty_output_ROOT_TEST_IOPTION->SetMarkerSize(1.5);
    pt_beauty_output_ROOT_TEST_IOPTION->SetMarkerColor(kSpring - 6);
    pt_beauty_output_ROOT_TEST_IOPTION->SetLineColor(kSpring - 6);

    pt_beauty_output_ROOT_TEST->Draw("LPSAME");
    pt_beauty_output_ROOT_TEST_IOPTION->Draw("LPSAME");
    // beauty_output_ROOFIT->Draw("LPSAME");

    TLegend *legend2 = new TLegend(0.23, 0.15, 0.5, 0.35);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend2->SetFillStyle(0);
    legend2->SetLineColor(kWhite);
    legend2->SetBorderSize(0);
    legend2->SetTextSize(0.0425);
    // legend2->SetHeader("Toy MC");
    // legend2->SetTextAlign(11);
    // legend2->AddEntry("charm_output_ROOFIT", "charm fraction from ROOFIT", "LP");
    // legend2->AddEntry("beauty_output_ROOFIT", "beauty fraction from ROOFIT", "LP");

    legend2->AddEntry("pt_charm_output_ROOT_TEST", "charm fraction from ROOT Fit", "LP");
    legend2->AddEntry("pt_beauty_output_ROOT_TEST", "beauty fraction from ROOT Fit", "LP");

    legend2->AddEntry("pt_charm_output_ROOT_TEST_IOPTION", "charm fraction from ROOT Fit with I option", "LP");
    legend2->AddEntry("pt_beauty_output_ROOT_TEST_IOPTION", "beauty fraction from ROOT Fit with I option", "LP");

    legend2->Draw("SAME");

    c1->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/pt_bin_size_test.pdf"));
    c1->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/pt_bin_size_test.png"));

    TCanvas *c2 = new TCanvas("c2", " ", 1000, 800);
    c2->cd();
    // c->Divide(2, 1);
    c2->cd();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    gStyle->SetOptStat(0);

    TH2D *MASS_h_grid = (TH2D *)PT_h_grid->Clone("MASS_h_grid");

    MASS_h_grid->GetXaxis()->SetTitle("#it{m}_{#mu#mu} bin size (GeV/#it{c}^{2})");
    MASS_h_grid->Draw();

    TGraph *MASS_charm_output_ROOT_TEST = new TGraph(n_scan, bin_size, MASS_charm_number_ROOT_TEST_norm);
    MASS_charm_output_ROOT_TEST->SetName("MASS_charm_output_ROOT_TEST");
    MASS_charm_output_ROOT_TEST->SetTitle(" ");

    MASS_charm_output_ROOT_TEST->SetMinimum(0.36);
    MASS_charm_output_ROOT_TEST->SetMarkerStyle(20);
    MASS_charm_output_ROOT_TEST->SetMarkerSize(1.5);
    MASS_charm_output_ROOT_TEST->SetMarkerColor(kMagenta - 2);
    MASS_charm_output_ROOT_TEST->SetLineColor(kMagenta - 2);

    TGraph *MASS_charm_output_ROOT_TEST_IOPTION = new TGraph(n_scan, bin_size, MASS_charm_number_ROOT_TEST_IOPTION_norm);
    MASS_charm_output_ROOT_TEST_IOPTION->SetName("MASS_charm_output_ROOT_TEST_IOPTION");
    MASS_charm_output_ROOT_TEST_IOPTION->SetTitle(" ");

    MASS_charm_output_ROOT_TEST_IOPTION->SetMinimum(0.36);
    MASS_charm_output_ROOT_TEST_IOPTION->SetMarkerStyle(24);
    MASS_charm_output_ROOT_TEST_IOPTION->SetMarkerSize(2.5);
    MASS_charm_output_ROOT_TEST_IOPTION->SetMarkerColor(kMagenta - 2);
    MASS_charm_output_ROOT_TEST_IOPTION->SetLineColor(kMagenta - 2);

    MASS_charm_output_ROOT_TEST->Draw("LPSAME");
    MASS_charm_output_ROOT_TEST_IOPTION->Draw("LPSAME");
    // charm_output_ROOFIT->Draw("LPSAME");

    TGraph *MASS_beauty_output_ROOT_TEST = new TGraph(n_scan, bin_size, MASS_beauty_number_ROOT_TEST_norm);
    MASS_beauty_output_ROOT_TEST->SetName("MASS_beauty_output_ROOT_TEST");
    MASS_beauty_output_ROOT_TEST->SetTitle(" ");

    MASS_beauty_output_ROOT_TEST->SetMinimum(0.36);
    MASS_beauty_output_ROOT_TEST->SetMarkerStyle(20);
    MASS_beauty_output_ROOT_TEST->SetMarkerSize(1.5);
    MASS_beauty_output_ROOT_TEST->SetMarkerColor(kSpring - 6);
    MASS_beauty_output_ROOT_TEST->SetLineColor(kSpring - 6);

    TGraph *MASS_beauty_output_ROOT_TEST_IOPTION = new TGraph(n_scan, bin_size, MASS_beauty_number_ROOT_TEST_IOPTION_norm);
    MASS_beauty_output_ROOT_TEST_IOPTION->SetName("MASS_beauty_output_ROOT_TEST_IOPTION");
    MASS_beauty_output_ROOT_TEST_IOPTION->SetTitle(" ");

    MASS_beauty_output_ROOT_TEST_IOPTION->SetMinimum(0.36);
    MASS_beauty_output_ROOT_TEST_IOPTION->SetMarkerStyle(24);
    MASS_beauty_output_ROOT_TEST_IOPTION->SetMarkerSize(2.5);
    MASS_beauty_output_ROOT_TEST_IOPTION->SetMarkerColor(kSpring - 6);
    MASS_beauty_output_ROOT_TEST_IOPTION->SetLineColor(kSpring - 6);

    MASS_beauty_output_ROOT_TEST->Draw("LPSAME");
    MASS_beauty_output_ROOT_TEST_IOPTION->Draw("LPSAME");

    // beauty_output_ROOFIT->Draw("LPSAME");

    TLegend *legend3 = new TLegend(0.23, 0.15, 0.5, 0.35);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend3->SetFillStyle(0);
    legend3->SetLineColor(kWhite);
    legend3->SetBorderSize(0);
    legend3->SetTextSize(0.0425);
    // legend3->SetHeader("Toy MC");
    // legend3->SetTextAlign(11);
    // legend3->AddEntry("charm_output_ROOFIT", "charm fraction from ROOFIT", "LP");
    // legend3->AddEntry("beauty_output_ROOFIT", "beauty fraction from ROOFIT", "LP");

    legend3->AddEntry("MASS_charm_output_ROOT_TEST", "charm fraction from ROOT Fit", "LP");
    legend3->AddEntry("MASS_beauty_output_ROOT_TEST", "beauty fraction from ROOT Fit", "LP");

    legend3->AddEntry("MASS_charm_output_ROOT_TEST_IOPTION", "charm fraction from ROOT Fit with I Option", "LP");
    legend3->AddEntry("MASS_beauty_output_ROOT_TEST_IOPTION", "beauty fraction from ROOT Fit with I Option", "LP");

    legend3->Draw();

    c2->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/mass_bin_size_test.pdf"));
    c2->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/mass_bin_size_test.png"));
}