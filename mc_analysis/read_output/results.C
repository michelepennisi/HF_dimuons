void LoadStyle();
void SetLegend(TLegend *);

void results(Bool_t saveResults = kFALSE){
    LoadStyle();
    TLatex *letexTitle = new TLatex();
    letexTitle -> SetTextSize(0.05);
    letexTitle -> SetNDC();
    letexTitle -> SetTextFont(42);

    //**************************************************
    // Using FONLL for scaling across different energies
    double y_scal[] = {0.0000, 0.5000, 1.0000, 1.5000, 2.0000, 2.5000, 3.0000, 3.5000, 4.0000, 4.5000, 5.0000};
    double dy_scal[] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double dSigma_bb_dy_FONLL_200GeV[] = {5.3450e+05, 4.9740e+05, 3.9740e+05, 2.6530e+05, 1.4160e+05, 5.5810e+04, 1.2740e+04, 1.8180e+02, 0.0000e+00, 0.0000e+00, 0.0000e+00};
    double stat_dSigma_bb_dy_FONLL_200GeV[11];
    double dSigma_bb_dy_FONLL_13TeV[]  = {3.4200e+07, 3.3950e+07, 3.3190e+07, 3.1930e+07, 3.0160e+07, 2.7950e+07, 2.5320e+07, 2.2160e+07, 1.8490e+07, 1.4420e+07, 1.0180e+07};
    double stat_dSigma_bb_dy_FONLL_13TeV[11];

    for (int i = 0;i < 11;i++) {
        dSigma_bb_dy_FONLL_200GeV[i] = dSigma_bb_dy_FONLL_200GeV[i] * 1e-06;
        dSigma_bb_dy_FONLL_13TeV[i] = dSigma_bb_dy_FONLL_13TeV[i] * 1e-06;
        stat_dSigma_bb_dy_FONLL_200GeV[i] = dSigma_bb_dy_FONLL_200GeV[i] * 0.10;
        stat_dSigma_bb_dy_FONLL_13TeV[i] = dSigma_bb_dy_FONLL_13TeV[i] * 0.10;
    }
    

    TGraphErrors *gra_stat_dSigma_bb_dy_FONLL_200GeV = new TGraphErrors(11, y_scal, dSigma_bb_dy_FONLL_200GeV, dy_scal, stat_dSigma_bb_dy_FONLL_200GeV);
    gra_stat_dSigma_bb_dy_FONLL_200GeV -> SetMarkerStyle(20);
    gra_stat_dSigma_bb_dy_FONLL_200GeV -> SetMarkerColor(kRed);
    gra_stat_dSigma_bb_dy_FONLL_200GeV -> SetLineWidth(2);
    gra_stat_dSigma_bb_dy_FONLL_200GeV -> SetLineColor(kRed);
    TF1 *func_dSigma_bb_dy_FONLL_200GeV = new TF1("func_dSigma_bb_dy_FONLL_200GeV", "gaus", 0, 5);
    gra_stat_dSigma_bb_dy_FONLL_200GeV -> Fit(func_dSigma_bb_dy_FONLL_200GeV);

    TGraphErrors *gra_stat_dSigma_bb_dy_FONLL_13TeV = new TGraphErrors(11, y_scal, dSigma_bb_dy_FONLL_13TeV, dy_scal, stat_dSigma_bb_dy_FONLL_13TeV);
    gra_stat_dSigma_bb_dy_FONLL_13TeV -> SetMarkerStyle(20);
    gra_stat_dSigma_bb_dy_FONLL_13TeV -> SetMarkerColor(kBlue);
    gra_stat_dSigma_bb_dy_FONLL_13TeV -> SetLineWidth(2);
    gra_stat_dSigma_bb_dy_FONLL_13TeV -> SetLineColor(kBlue);
    TF1 *func_dSigma_bb_dy_FONLL_13TeV = new TF1("func_dSigma_bb_dy_FONLL_13TeV", "gaus", 0, 5);
    gra_stat_dSigma_bb_dy_FONLL_13TeV -> Fit(func_dSigma_bb_dy_FONLL_13TeV);

    double integral_dSigma_bb_dy_FONLL_200GeV = func_dSigma_bb_dy_FONLL_200GeV -> Integral(1.2, 2.2);
    double integral_dSigma_bb_dy_FONLL_13TeV = func_dSigma_bb_dy_FONLL_13TeV -> Integral(1.2, 2.2);
    double ratio_FONLL_200GeV_to_13TeV = integral_dSigma_bb_dy_FONLL_13TeV / integral_dSigma_bb_dy_FONLL_200GeV;
    std::cout << "pp 200 GeV integral = " << integral_dSigma_bb_dy_FONLL_200GeV << std::endl;
    std::cout << "pp 13 TeV integral = " << integral_dSigma_bb_dy_FONLL_13TeV << std::endl;
    std::cout << "Ratio 13 TeV / 200 GeV = " << ratio_FONLL_200GeV_to_13TeV << std::endl;

    double y_PHENIX[] = {1.7};
    double dy_PHENIX[] = {0.5};
    double dSigma_bb_dy_PHENIX[] = {3.75};
    double stat_dSigma_bb_dy_PHENIX[] = {0.24};
    double syst_dSigma_bb_dy_PHENIX[] = {0.67};

    //dSigma_bb_dy_PHENIX[0] = dSigma_bb_dy_PHENIX[0] * ratio_FONLL_200GeV_to_13TeV;
    //stat_dSigma_bb_dy_PHENIX[0] = stat_dSigma_bb_dy_PHENIX[0] * ratio_FONLL_200GeV_to_13TeV;
    //syst_dSigma_bb_dy_PHENIX[0] = syst_dSigma_bb_dy_PHENIX[0] * ratio_FONLL_200GeV_to_13TeV;

    dSigma_bb_dy_PHENIX[0] = dSigma_bb_dy_PHENIX[0] * ratio_FONLL_200GeV_to_13TeV / 8.3;
    stat_dSigma_bb_dy_PHENIX[0] = stat_dSigma_bb_dy_PHENIX[0] * ratio_FONLL_200GeV_to_13TeV / 8.3;
    syst_dSigma_bb_dy_PHENIX[0] = syst_dSigma_bb_dy_PHENIX[0] * ratio_FONLL_200GeV_to_13TeV / 8.3;

    TGraphErrors *gra_stat_dSigma_bb_dy_PHENIX = new TGraphErrors(1, y_PHENIX, dSigma_bb_dy_PHENIX, dy_PHENIX, stat_dSigma_bb_dy_PHENIX);
    gra_stat_dSigma_bb_dy_PHENIX -> SetMarkerStyle(20);
    gra_stat_dSigma_bb_dy_PHENIX -> SetMarkerColor(kGreen+1);
    gra_stat_dSigma_bb_dy_PHENIX -> SetLineWidth(2);
    gra_stat_dSigma_bb_dy_PHENIX -> SetLineColor(kGreen+1);

    TGraphErrors *gra_syst_dSigma_bb_dy_PHENIX = new TGraphErrors(1, y_PHENIX, dSigma_bb_dy_PHENIX, dy_PHENIX, syst_dSigma_bb_dy_PHENIX);
    gra_syst_dSigma_bb_dy_PHENIX -> SetMarkerStyle(20);
    gra_syst_dSigma_bb_dy_PHENIX -> SetMarkerColor(kGreen+1);
    gra_syst_dSigma_bb_dy_PHENIX -> SetLineWidth(2);
    gra_syst_dSigma_bb_dy_PHENIX -> SetLineColor(kGreen+1);
    gra_syst_dSigma_bb_dy_PHENIX -> SetFillStyle(0);

    TCanvas *canvas_dSigma_bb_dy_scal = new TCanvas("canvas_dSigma_bb_dy_scal","",800,600);
    canvas_dSigma_bb_dy_scal -> SetFillColor(0);
    canvas_dSigma_bb_dy_scal -> SetBorderMode(0);
    canvas_dSigma_bb_dy_scal -> SetBorderSize(0);
    canvas_dSigma_bb_dy_scal -> SetTickx(1);
    canvas_dSigma_bb_dy_scal -> SetTicky(1);
    canvas_dSigma_bb_dy_scal -> SetLeftMargin(0.15);
    canvas_dSigma_bb_dy_scal -> SetBottomMargin(0.1518219);
    canvas_dSigma_bb_dy_scal -> SetFrameBorderMode(0);
    canvas_dSigma_bb_dy_scal -> SetFrameBorderMode(0);
    gPad -> SetLogy(1);
    
    TH2D *hist_grid_dSigma_bb_dy_scal = new TH2D("hist_grid_dSigma_bb_dy_scal","",100, 0, 5., 100, 0.001, 10000);
    hist_grid_dSigma_bb_dy_scal -> GetXaxis() -> SetTitle("#it{y}");
    hist_grid_dSigma_bb_dy_scal -> GetXaxis() -> SetTitleOffset(1.2);
    hist_grid_dSigma_bb_dy_scal -> GetXaxis() -> SetTitleSize(0.05);
    hist_grid_dSigma_bb_dy_scal -> GetXaxis() -> SetLabelSize(0.05);
    hist_grid_dSigma_bb_dy_scal -> GetYaxis() -> SetTitle("d#sigma / d#it{y}");
    hist_grid_dSigma_bb_dy_scal -> GetYaxis() -> SetTitleOffset(1.2);
    hist_grid_dSigma_bb_dy_scal -> GetYaxis() -> SetTitleSize(0.055);
    hist_grid_dSigma_bb_dy_scal -> GetYaxis() -> SetLabelSize(0.05);
    hist_grid_dSigma_bb_dy_scal -> Draw();

    gra_stat_dSigma_bb_dy_FONLL_200GeV -> Draw("EPsame");
    gra_stat_dSigma_bb_dy_FONLL_13TeV -> Draw("EPsame");
    gra_stat_dSigma_bb_dy_PHENIX -> Draw("EPsame");
    gra_syst_dSigma_bb_dy_PHENIX -> Draw("E2Psame");
    //**************************************************

    // PYTHIA results
    double y[] = {0, 3.25};
    double dy[] = {0.5, 0.75};

    // // Beauty cross section
    double dSigma_bb_dy_PYTHIA[] = {79, 21.8};
    double stat_dSigma_bb_dy_PYTHIA[] = {14, 6};
    double syst_dSigma_bb_dy_PYTHIA[] = {11, 6};

    TGraphErrors *gra_stat_dSigma_bb_dy_PYTHIA = new TGraphErrors(2, y, dSigma_bb_dy_PYTHIA, dy, stat_dSigma_bb_dy_PYTHIA);
    gra_stat_dSigma_bb_dy_PYTHIA -> SetMarkerStyle(20);
    gra_stat_dSigma_bb_dy_PYTHIA -> SetMarkerColor(kRed);
    gra_stat_dSigma_bb_dy_PYTHIA -> SetLineWidth(2);
    gra_stat_dSigma_bb_dy_PYTHIA -> SetLineColor(kRed);

    TGraphErrors *gra_syst_dSigma_bb_dy_PYTHIA = new TGraphErrors(2, y, dSigma_bb_dy_PYTHIA, dy, syst_dSigma_bb_dy_PYTHIA);
    gra_syst_dSigma_bb_dy_PYTHIA -> SetMarkerStyle(20);
    gra_syst_dSigma_bb_dy_PYTHIA -> SetMarkerColor(kRed);
    gra_syst_dSigma_bb_dy_PYTHIA -> SetLineWidth(2);
    gra_syst_dSigma_bb_dy_PYTHIA -> SetLineColor(kRed);
    gra_syst_dSigma_bb_dy_PYTHIA -> SetFillStyle(0);

    // // Charm cross section
    double dSigma_cc_dy_PYTHIA[] = {974, 545.9};
    double stat_dSigma_cc_dy_PYTHIA[] = {138, 56.6};
    double syst_dSigma_cc_dy_PYTHIA[] = {140, 56.6};

    TGraphErrors *gra_stat_dSigma_cc_dy_PYTHIA = new TGraphErrors(2, y, dSigma_cc_dy_PYTHIA, dy, stat_dSigma_cc_dy_PYTHIA);
    gra_stat_dSigma_cc_dy_PYTHIA -> SetMarkerStyle(20);
    gra_stat_dSigma_cc_dy_PYTHIA -> SetMarkerColor(kRed);
    gra_stat_dSigma_cc_dy_PYTHIA -> SetLineWidth(2);
    gra_stat_dSigma_cc_dy_PYTHIA -> SetLineColor(kRed);

    TGraphErrors *gra_syst_dSigma_cc_dy_PYTHIA = new TGraphErrors(2, y, dSigma_cc_dy_PYTHIA, dy, syst_dSigma_cc_dy_PYTHIA);
    gra_syst_dSigma_cc_dy_PYTHIA -> SetMarkerStyle(20);
    gra_syst_dSigma_cc_dy_PYTHIA -> SetMarkerColor(kRed);
    gra_syst_dSigma_cc_dy_PYTHIA -> SetLineWidth(2);
    gra_syst_dSigma_cc_dy_PYTHIA -> SetLineColor(kRed);
    gra_syst_dSigma_cc_dy_PYTHIA -> SetFillStyle(0);

    // FONLL results
    double y_FONLL[] = {-3.7500,-3.2500,-2.7500,-2.2500,-1.7500,-1.2500,-0.7500,-0.2500,0.2500,0.7500,1.2500,1.7500,2.2500,2.7500,3.2500,3.7500};
    const int size_y_FONLL = sizeof(y_FONLL)/sizeof(y_FONLL[0]);
    double dy_min_FONLL[size_y_FONLL];
    double dy_max_FONLL[size_y_FONLL];
    
    double dSigma_bb_dy_FONLL[] = {2.9500e+07,3.5800e+07,4.1430e+07,4.6310e+07,5.0380e+07,5.3500e+07,5.5600e+07,5.6650e+07,5.6650e+07,5.5600e+07,5.3500e+07,5.0380e+07,4.6310e+07,4.1430e+07,3.5800e+07,2.9500e+07};
    double dev_min_dSigma_bb_dy_FONLL[] = {1.9800e+07,2.3440e+07,2.6520e+07,2.9100e+07,3.1230e+07,3.2860e+07,3.3940e+07,3.4480e+07,3.4480e+07,3.3940e+07,3.2860e+07,3.1230e+07,2.9100e+07,2.6520e+07,2.3440e+07,1.9800e+07};
    double dev_max_dSigma_bb_dy_FONLL[] = {4.2250e+07,5.1100e+07,5.9000e+07,6.5810e+07,7.1490e+07,7.5850e+07,7.8790e+07,8.0260e+07,8.0260e+07,7.8790e+07,7.5850e+07,7.1490e+07,6.5810e+07,5.9000e+07,5.1100e+07,4.2250e+07};

    for (int i = 0;i < size_y_FONLL;i++) {
        dy_min_FONLL[i] = 0.25;
        dy_max_FONLL[i] = 0.25;

        dev_min_dSigma_bb_dy_FONLL[i] = (dSigma_bb_dy_FONLL[i] - dev_min_dSigma_bb_dy_FONLL[i]) * 1e-06;
        dev_max_dSigma_bb_dy_FONLL[i] = (dev_max_dSigma_bb_dy_FONLL[i] - dSigma_bb_dy_FONLL[i]) * 1e-06;
        dSigma_bb_dy_FONLL[i] = dSigma_bb_dy_FONLL[i] * 1e-06;
    }

    TGraphAsymmErrors *gra_dSigma_bb_dy_FONLL = new TGraphAsymmErrors(size_y_FONLL, y_FONLL, dSigma_bb_dy_FONLL, dy_min_FONLL, dy_max_FONLL, dev_min_dSigma_bb_dy_FONLL, dev_max_dSigma_bb_dy_FONLL);
    gra_dSigma_bb_dy_FONLL -> SetMarkerStyle(20);
    gra_dSigma_bb_dy_FONLL -> SetMarkerColor(kGray+1);
    gra_dSigma_bb_dy_FONLL -> SetLineWidth(2);
    gra_dSigma_bb_dy_FONLL -> SetLineColorAlpha(kGray+1, 0.5);
    gra_dSigma_bb_dy_FONLL -> SetFillColorAlpha(kGray+1, 0.5);

    double dSigma_cc_dy_FONLL[] = {6.3410e+08,6.5330e+08,6.6090e+08,6.5980e+08,6.5370e+08,6.4710e+08,6.4690e+08,6.5020e+08,6.5020e+08,6.4690e+08,6.4710e+08,6.5370e+08,6.5980e+08,6.6090e+08,6.5330e+08,6.3410e+08};
    double dev_min_dSigma_cc_dy_FONLL[] = {2.2060e+08,2.0890e+08,1.9410e+08,1.7760e+08,1.6010e+08,1.5140e+08,1.5170e+08,1.5250e+08,1.5250e+08,1.5170e+08,1.5140e+08,1.6010e+08,1.7760e+08,1.9410e+08,2.0890e+08,2.2060e+08};
    double dev_max_dSigma_cc_dy_FONLL[] = {1.4990e+09,1.5220e+09,1.5200e+09,1.5000e+09,1.4690e+09,1.4370e+09,1.4230e+09,1.4250e+09,1.4250e+09,1.4230e+09,1.4370e+09,1.4690e+09,1.5000e+09,1.5200e+09,1.5220e+09,1.4990e+09};

    for (int i = 0;i < size_y_FONLL;i++) {
        dev_min_dSigma_cc_dy_FONLL[i] = (dSigma_cc_dy_FONLL[i] - dev_min_dSigma_cc_dy_FONLL[i]) * 1e-06;
        dev_max_dSigma_cc_dy_FONLL[i] = (dev_max_dSigma_cc_dy_FONLL[i] - dSigma_cc_dy_FONLL[i]) * 1e-06;
        dSigma_cc_dy_FONLL[i] = dSigma_cc_dy_FONLL[i] * 1e-06;
    }

    TGraphAsymmErrors *gra_dSigma_cc_dy_FONLL = new TGraphAsymmErrors(size_y_FONLL, y_FONLL, dSigma_cc_dy_FONLL, dy_min_FONLL, dy_max_FONLL, dev_min_dSigma_cc_dy_FONLL, dev_max_dSigma_cc_dy_FONLL);
    gra_dSigma_cc_dy_FONLL -> SetMarkerStyle(20);
    gra_dSigma_cc_dy_FONLL -> SetMarkerColor(kGray+1);
    gra_dSigma_cc_dy_FONLL -> SetLineWidth(2);
    gra_dSigma_cc_dy_FONLL -> SetLineColorAlpha(kGray+1, 0.5);
    gra_dSigma_cc_dy_FONLL -> SetFillColorAlpha(kGray+1, 0.5);






    TCanvas *canvas_dSigma_bb_dy = new TCanvas("canvas_dSigma_bb_dy","",800,600);
    canvas_dSigma_bb_dy -> SetFillColor(0);
    canvas_dSigma_bb_dy -> SetBorderMode(0);
    canvas_dSigma_bb_dy -> SetBorderSize(0);
    canvas_dSigma_bb_dy -> SetTickx(1);
    canvas_dSigma_bb_dy -> SetTicky(1);
    canvas_dSigma_bb_dy -> SetLeftMargin(0.15);
    canvas_dSigma_bb_dy -> SetBottomMargin(0.1518219);
    canvas_dSigma_bb_dy -> SetFrameBorderMode(0);
    canvas_dSigma_bb_dy -> SetFrameBorderMode(0);
    gPad -> SetLogy(1);

    TH2D *hist_grid_dSigma_bb_dy = new TH2D("hist_grid_dSigma_bb_dy","",100, -5., 5., 100, 0.1, 10000);
    hist_grid_dSigma_bb_dy -> GetXaxis() -> SetTitle("#it{y}");
    hist_grid_dSigma_bb_dy -> GetXaxis() -> SetTitleOffset(1.2);
    hist_grid_dSigma_bb_dy -> GetXaxis() -> SetTitleSize(0.05);
    hist_grid_dSigma_bb_dy -> GetXaxis() -> SetLabelSize(0.05);
    hist_grid_dSigma_bb_dy -> GetYaxis() -> SetTitle("d#sigma / d#it{y}");
    hist_grid_dSigma_bb_dy -> GetYaxis() -> SetTitleOffset(1.2);
    hist_grid_dSigma_bb_dy -> GetYaxis() -> SetTitleSize(0.055);
    hist_grid_dSigma_bb_dy -> GetYaxis() -> SetLabelSize(0.05);
    hist_grid_dSigma_bb_dy -> Draw();

    gra_dSigma_bb_dy_FONLL -> Draw("E2same");
    gra_stat_dSigma_bb_dy_PYTHIA -> Draw("EPsame");
    gra_syst_dSigma_bb_dy_PYTHIA -> Draw("E2same");
    gra_stat_dSigma_bb_dy_PHENIX -> Draw("EPsame");
    gra_syst_dSigma_bb_dy_PHENIX -> Draw("E2same");

    TLegend *legend_dSigma_bb_dy = new TLegend(0.18,0.18,0.35,0.35," ","brNDC");
    SetLegend(legend_dSigma_bb_dy);
    legend_dSigma_bb_dy -> AddEntry(gra_stat_dSigma_bb_dy_PYTHIA,"Data","PE");
    legend_dSigma_bb_dy -> AddEntry(gra_dSigma_bb_dy_FONLL,"FONLL (uncertainty from scale only)","F");
    legend_dSigma_bb_dy -> Draw("same");

    letexTitle -> DrawLatex(0.2,0.82,"ALICE Preliminary, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle -> DrawLatex(0.2,0.75,"Inclusive J/#psi #rightarrow #mu^{#plus}#mu^{#minus}");
    letexTitle -> DrawLatex(0.2,0.68,"2 < #it{p}_{T} < 6 GeV/#it{c} , 2.5 < #it{y} < 4");






    TCanvas *canvas_dSigma_cc_dy = new TCanvas("canvas_dSigma_cc_dy","",800,600);
    canvas_dSigma_cc_dy -> SetFillColor(0);
    canvas_dSigma_cc_dy -> SetBorderMode(0);
    canvas_dSigma_cc_dy -> SetBorderSize(0);
    canvas_dSigma_cc_dy -> SetTickx(1);
    canvas_dSigma_cc_dy -> SetTicky(1);
    canvas_dSigma_cc_dy -> SetLeftMargin(0.15);
    canvas_dSigma_cc_dy -> SetBottomMargin(0.1518219);
    canvas_dSigma_cc_dy -> SetFrameBorderMode(0);
    canvas_dSigma_cc_dy -> SetFrameBorderMode(0);
    gPad -> SetLogy(1);

    TH2D *hist_grid_dSigma_cc_dy = new TH2D("hist_grid_dSigma_cc_dy","",100, -5., 5., 100, 10, 100000);
    hist_grid_dSigma_cc_dy -> GetXaxis() -> SetTitle("#it{y}");
    hist_grid_dSigma_cc_dy -> GetXaxis() -> SetTitleOffset(1.2);
    hist_grid_dSigma_cc_dy -> GetXaxis() -> SetTitleSize(0.05);
    hist_grid_dSigma_cc_dy -> GetXaxis() -> SetLabelSize(0.05);
    hist_grid_dSigma_cc_dy -> GetYaxis() -> SetTitle("d#sigma / d#it{y}");
    hist_grid_dSigma_cc_dy -> GetYaxis() -> SetTitleOffset(1.2);
    hist_grid_dSigma_cc_dy -> GetYaxis() -> SetTitleSize(0.055);
    hist_grid_dSigma_cc_dy -> GetYaxis() -> SetLabelSize(0.05);
    hist_grid_dSigma_cc_dy -> Draw();

    gra_dSigma_cc_dy_FONLL -> Draw("E2same");
    gra_stat_dSigma_cc_dy_PYTHIA -> Draw("EPsame");
    gra_syst_dSigma_cc_dy_PYTHIA -> Draw("E2same");

    TLegend *legend_dSigma_cc_dy = new TLegend(0.18,0.18,0.35,0.35," ","brNDC");
    SetLegend(legend_dSigma_cc_dy);
    legend_dSigma_cc_dy -> AddEntry(gra_stat_dSigma_cc_dy_PYTHIA,"Data","PE");
    legend_dSigma_cc_dy -> AddEntry(gra_dSigma_cc_dy_FONLL,"FONLL (uncertainty from scale only)","F");
    legend_dSigma_cc_dy -> Draw("same");

    letexTitle -> DrawLatex(0.2,0.82,"ALICE Preliminary, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle -> DrawLatex(0.2,0.75,"Inclusive J/#psi #rightarrow #mu^{#plus}#mu^{#minus}");
    letexTitle -> DrawLatex(0.2,0.68,"2 < #it{p}_{T} < 6 GeV/#it{c} , 2.5 < #it{y} < 4");
    
    


    

    if (saveResults) {
        canvas_dSigma_bb_dy -> SaveAs("bb_cross_section.pdf");
        canvas_dSigma_cc_dy -> SaveAs("cc_cross_section.pdf");
    }
}
////////////////////////////////////////////////////////////////////////////////
void LoadStyle(){
    int font = 42;
    //TGaxis::SetMaxDigits(2);
    gStyle -> SetFrameBorderMode(0);
    gStyle -> SetFrameFillColor(0);
    gStyle -> SetCanvasBorderMode(0);
    gStyle -> SetPadBorderMode(0);
    gStyle -> SetPadColor(10);
    gStyle -> SetCanvasColor(10);
    gStyle -> SetTitleFillColor(10);
    gStyle -> SetTitleBorderSize(1);
    gStyle -> SetStatColor(10);
    gStyle -> SetStatBorderSize(1);
    gStyle -> SetLegendBorderSize(1);
    gStyle -> SetDrawBorder(0);
    gStyle -> SetTextFont(font);
    gStyle -> SetStatFontSize(0.05);
    gStyle -> SetStatX(0.97);
    gStyle -> SetStatY(0.98);
    gStyle -> SetStatH(0.03);
    gStyle -> SetStatW(0.3);
    gStyle -> SetTickLength(0.02,"y");
    gStyle -> SetEndErrorSize(3);
    gStyle -> SetLabelSize(0.04,"xyz");
    gStyle -> SetLabelFont(font,"xyz");
    gStyle -> SetLabelOffset(0.01,"xyz");
    gStyle -> SetTitleFont(font,"xyz");
    gStyle -> SetTitleOffset(0.9,"x");
    gStyle -> SetTitleOffset(1.02,"y");
    gStyle -> SetTitleSize(0.04,"xyz");
    gStyle -> SetMarkerSize(1.3);
    gStyle -> SetOptStat(0);
    gStyle -> SetEndErrorSize(0);
    gStyle -> SetCanvasPreferGL(kTRUE);
    gStyle -> SetHatchesSpacing(0.5);
}
////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.05);
}
