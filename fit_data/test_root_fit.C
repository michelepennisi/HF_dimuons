double FuncPtMass(double *x, double *par);
double total(double *x, double *par);

TPad *divide2(TH1D *h_data, TF1 *model,int min);

void test_root_fit(Int_t N_rebin = 10,TString opt="LR0I")
{   
    Double_t Binning_m=260./N_rebin;
    Double_t Binning_pt=300./N_rebin;

    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");

    TFile *fIn = new TFile(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/root_files/pdfMC_unbinned.root"));

    RooWorkspace *w = (RooWorkspace *)fIn->Get("w");
    w->Print();

    RooRealVar *B_DimuMassFromCharm = w->var("B_DimuMassFromCharm");
    B_DimuMassFromCharm->setConstant(kTRUE);
    printf("B_DimuMassFromCharm %0.3f", B_DimuMassFromCharm->getVal());
    RooRealVar *n1_DimuMassFromCharm = w->var("n1_DimuMassFromCharm");
    n1_DimuMassFromCharm->setConstant(kTRUE);
    RooRealVar *n2_DimuMassFromCharm = w->var("n2_DimuMassFromCharm");
    n2_DimuMassFromCharm->setConstant(kTRUE);

    RooRealVar *B_DimuMassFromBeauty = w->var("B_DimuMassFromBeauty");
    B_DimuMassFromBeauty->setConstant(kTRUE);
    RooRealVar *n1_DimuMassFromBeauty = w->var("n1_DimuMassFromBeauty");
    n1_DimuMassFromBeauty->setConstant(kTRUE);
    RooRealVar *n2_DimuMassFromBeauty = w->var("n2_DimuMassFromBeauty");
    n2_DimuMassFromBeauty->setConstant(kTRUE);

    RooRealVar *B_DimuMassFromMixed = w->var("B_DimuMassFromMixed");
    B_DimuMassFromMixed->setConstant(kTRUE);
    RooRealVar *n1_DimuMassFromMixed = w->var("n1_DimuMassFromMixed");
    n1_DimuMassFromMixed->setConstant(kTRUE);
    RooRealVar *n2_DimuMassFromMixed = w->var("n2_DimuMassFromMixed");
    n2_DimuMassFromMixed->setConstant(kTRUE);

    RooRealVar *B_DimuPtFromCharm = w->var("B_DimuPtFromCharm");
    B_DimuPtFromCharm->setConstant(kTRUE);
    RooRealVar *n1_DimuPtFromCharm = w->var("n1_DimuPtFromCharm");
    n1_DimuPtFromCharm->setConstant(kTRUE);
    RooRealVar *n2_DimuPtFromCharm = w->var("n2_DimuPtFromCharm");
    n2_DimuPtFromCharm->setConstant(kTRUE);

    RooRealVar *B_DimuPtFromBeauty = w->var("B_DimuPtFromBeauty");
    B_DimuPtFromBeauty->setConstant(kTRUE);
    RooRealVar *n1_DimuPtFromBeauty = w->var("n1_DimuPtFromBeauty");
    n1_DimuPtFromBeauty->setConstant(kTRUE);
    RooRealVar *n2_DimuPtFromBeauty = w->var("n2_DimuPtFromBeauty");
    n2_DimuPtFromBeauty->setConstant(kTRUE);

    RooRealVar *B_DimuPtFromMixed = w->var("B_DimuPtFromMixed");
    B_DimuPtFromMixed->setConstant(kTRUE);
    RooRealVar *n1_DimuPtFromMixed = w->var("n1_DimuPtFromMixed");
    n1_DimuPtFromMixed->setConstant(kTRUE);
    RooRealVar *n2_DimuPtFromMixed = w->var("n2_DimuPtFromMixed");
    n2_DimuPtFromMixed->setConstant(kTRUE);

    TFile *fIn_data = new TFile("/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/HistResults_merged.root", "READ");

    TH2D *histDimuPtM_fromdata = (TH2D *)fIn_data->Get(Form("Dimuon/CMUL7/DQ_cut_match_LT_Ycut_ULS/h_PtMdiMu_CMUL7_DQ_cut_match_LT_Ycut_ULS"));
    histDimuPtM_fromdata->SetName("histDimuPtM_fromdata");

    TH1D *histDimuPt = (TH1D *)histDimuPtM_fromdata->ProjectionX();
    histDimuPt->Sumw2();
    TH1D *histDimuMass = (TH1D *)histDimuPtM_fromdata->ProjectionY();
    histDimuMass->Sumw2();

    for (Int_t i = 0; i < histDimuMass->GetNbinsX(); i++)
    {
        if (histDimuMass->GetBinLowEdge(i)>9.0 && histDimuMass->GetBinLowEdge(i)<11.0)
        {
            printf("Edge %0.5f\n",histDimuMass->GetBinLowEdge(i));
            printf("Error %0.5f\n",histDimuMass->GetBinError(i));
            histDimuMass->SetBinContent(i,0);
            histDimuMass->SetBinError(i,0);
        }
        
    }
    

    histDimuMass->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    histDimuMass->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    histDimuMass->Rebin(N_rebin);
    histDimuMass->Scale(1., "width");

    histDimuMass->SetMarkerStyle(20);
    histDimuMass->SetMarkerSize(1.5);
    histDimuMass->SetMarkerColor(kBlack);
    histDimuMass->SetLineColor(kBlack);
    histDimuMass->SetLineWidth(2.0);

    auto mass_pdf_total = new TF1("mass_pdf_total", total, 4, 30, 12);
    mass_pdf_total->FixParameter(0, B_DimuMassFromCharm->getVal());
    mass_pdf_total->FixParameter(1, n1_DimuMassFromCharm->getVal());
    mass_pdf_total->FixParameter(2, n2_DimuMassFromCharm->getVal());

    mass_pdf_total->FixParameter(3, B_DimuMassFromBeauty->getVal());
    mass_pdf_total->FixParameter(4, n1_DimuMassFromBeauty->getVal());
    mass_pdf_total->FixParameter(5, n2_DimuMassFromBeauty->getVal());

    mass_pdf_total->FixParameter(6, B_DimuMassFromMixed->getVal());
    mass_pdf_total->FixParameter(7, n1_DimuMassFromMixed->getVal());
    mass_pdf_total->FixParameter(8, n2_DimuMassFromMixed->getVal());

    mass_pdf_total->SetLineColor(kRed);
    mass_pdf_total->SetLineWidth(2);
    mass_pdf_total->SetParameter(9, 500000);
    mass_pdf_total->SetParameter(10, 200000);
    mass_pdf_total->FixParameter(11, 225000);

    histDimuMass->Fit(mass_pdf_total, opt.Data());
    auto mass_pdf_Charm = new TF1("mass_pdf_Charm", FuncPtMass, 4, 30, 4);
    mass_pdf_Charm->SetParameter(3, mass_pdf_total->GetParameter(9));
    mass_pdf_Charm->FixParameter(0, B_DimuMassFromCharm->getVal());
    mass_pdf_Charm->FixParameter(1, n1_DimuMassFromCharm->getVal());
    mass_pdf_Charm->FixParameter(2, n2_DimuMassFromCharm->getVal());
    mass_pdf_Charm->SetLineColor(kMagenta - 2);
    mass_pdf_Charm->SetLineWidth(2);
    mass_pdf_Charm->SetLineStyle(9);
    // mass_pdf_Charm->Draw();

    auto mass_pdf_Beauty = new TF1("mass_pdf_Beauty", FuncPtMass, 4, 30, 4);
    mass_pdf_Beauty->FixParameter(3, mass_pdf_total->GetParameter(10));
    mass_pdf_Beauty->FixParameter(0, B_DimuMassFromBeauty->getVal());
    mass_pdf_Beauty->FixParameter(1, n1_DimuMassFromBeauty->getVal());
    mass_pdf_Beauty->FixParameter(2, n2_DimuMassFromBeauty->getVal());
    mass_pdf_Beauty->SetLineColor(kSpring - 6);
    mass_pdf_Beauty->SetLineWidth(2);
    mass_pdf_Beauty->SetLineStyle(9);
    // mass_pdf_Beauty->Draw("SAME");

    auto mass_pdf_Mixed = new TF1("mass_pdf_Mixed", FuncPtMass, 4, 30, 4);
    mass_pdf_Mixed->FixParameter(3, mass_pdf_total->GetParameter(11));
    mass_pdf_Mixed->FixParameter(0, B_DimuMassFromMixed->getVal());
    mass_pdf_Mixed->FixParameter(1, n1_DimuMassFromMixed->getVal());
    mass_pdf_Mixed->FixParameter(2, n2_DimuMassFromMixed->getVal());
    mass_pdf_Mixed->SetLineColor(kCyan + 2);
    mass_pdf_Mixed->SetLineWidth(2);
    mass_pdf_Mixed->SetLineStyle(9);

    printf("--------------------------------------------------- \n");
    printf("From Mass fit \n");
    printf("N Charm: %0.0f \n", mass_pdf_Charm->Integral(4, 30));
    printf("N Beauty: %0.0f \n", mass_pdf_Beauty->Integral(4, 30));
    printf("N Mixed: %0.0f \n", mass_pdf_Mixed->Integral(4, 30));
    printf("--------------------------------------------------- \n");

    gStyle->SetOptStat(0);
    TCanvas *c_mass = new TCanvas("mass_fit", "mass_fit", 1000, 800);

    c_mass->cd();
    histDimuMass->SetMaximum(1.8e+05);
    histDimuMass->SetMinimum(1.8e-02);
    TPad *pads = divide2(histDimuMass, mass_pdf_total,4);

    pads[0].cd();
    TLegend *legend = new TLegend(0.175, 0.15, 0.35, 0.4);

    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0525);
    // legend->SetHeader("Data");
    // legend->SetTextAlign(11);
    legend->AddEntry(histDimuMass, "Data");
    legend->AddEntry(mass_pdf_total, "Cocktail");
    legend->AddEntry(mass_pdf_Charm, "#mu#mu #leftarrow c#bar{c}");
    legend->AddEntry(mass_pdf_Beauty, "#mu#mu #leftarrow b#bar{b}");
    legend->AddEntry(mass_pdf_Mixed, "#mu#mu #leftarrow c,b");

    TLatex *letexTitle = new TLatex();
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);
    letexTitle->SetTextSize(0.055);
    letexTitle->DrawLatex(0.625, 0.825, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f", mass_pdf_Charm->Integral(4, 30)));
    letexTitle->DrawLatex(0.625, 0.725, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f", mass_pdf_Beauty->Integral(4, 30)));
    // letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.2e #pm %0.2e", fit_output[2]->getVal(), fit_output[2]->getError()));
    letexTitle->DrawLatex(0.625, 0.625, Form("#it{N}^{fixed}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f ", mass_pdf_Mixed->Integral(4, 30)));

    mass_pdf_Charm->Draw("SAME");
    mass_pdf_Beauty->Draw("SAME");
    mass_pdf_Mixed->Draw("SAME");
    legend->Draw();
    c_mass->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/mass_root_test_%0.0fmbin_%s_Opt.pdf",Binning_m,opt.Data()));
    c_mass->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/mass_root_test_%0.0fmbin_%s_Opt.png",Binning_m,opt.Data()));
    //----------------------------------------------------------------------//
    histDimuPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    histDimuPt->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})");
    histDimuPt->Rebin(N_rebin);
    histDimuPt->Scale(1., "width");

    histDimuPt->SetMarkerStyle(20);
    histDimuPt->SetMarkerSize(1.5);
    histDimuPt->SetMarkerColor(kBlack);
    histDimuPt->SetLineColor(kBlack);
    histDimuPt->SetLineWidth(2.0);

    auto Pt_pdf_total = new TF1("Pt_pdf_total", total, 0, 30, 12);
    Pt_pdf_total->FixParameter(0, B_DimuPtFromCharm->getVal());
    Pt_pdf_total->FixParameter(1, n1_DimuPtFromCharm->getVal());
    Pt_pdf_total->FixParameter(2, n2_DimuPtFromCharm->getVal());

    Pt_pdf_total->FixParameter(3, B_DimuPtFromBeauty->getVal());
    Pt_pdf_total->FixParameter(4, n1_DimuPtFromBeauty->getVal());
    Pt_pdf_total->FixParameter(5, n2_DimuPtFromBeauty->getVal());

    Pt_pdf_total->FixParameter(6, B_DimuPtFromMixed->getVal());
    Pt_pdf_total->FixParameter(7, n1_DimuPtFromMixed->getVal());
    Pt_pdf_total->FixParameter(8, n2_DimuPtFromMixed->getVal());

    Pt_pdf_total->SetLineColor(kRed);
    Pt_pdf_total->SetLineWidth(2);
    Pt_pdf_total->SetParameter(9, 500000);
    Pt_pdf_total->SetParameter(10, 200000);
    Pt_pdf_total->FixParameter(11, 500);

    histDimuPt->Fit(Pt_pdf_total, opt.Data());
    auto Pt_pdf_Charm = new TF1("Pt_pdf_Charm", FuncPtMass, 0, 30, 4);
    Pt_pdf_Charm->SetParameter(3, Pt_pdf_total->GetParameter(9));
    Pt_pdf_Charm->FixParameter(0, B_DimuPtFromCharm->getVal());
    Pt_pdf_Charm->FixParameter(1, n1_DimuPtFromCharm->getVal());
    Pt_pdf_Charm->FixParameter(2, n2_DimuPtFromCharm->getVal());
    Pt_pdf_Charm->SetLineColor(kMagenta - 2);
    Pt_pdf_Charm->SetLineWidth(2);
    Pt_pdf_Charm->SetLineStyle(9);
    // Pt_pdf_Charm->Draw();

    auto Pt_pdf_Beauty = new TF1("Pt_pdf_Beauty", FuncPtMass, 0, 30, 4);
    Pt_pdf_Beauty->FixParameter(3, Pt_pdf_total->GetParameter(10));
    Pt_pdf_Beauty->FixParameter(0, B_DimuPtFromBeauty->getVal());
    Pt_pdf_Beauty->FixParameter(1, n1_DimuPtFromBeauty->getVal());
    Pt_pdf_Beauty->FixParameter(2, n2_DimuPtFromBeauty->getVal());
    Pt_pdf_Beauty->SetLineColor(kSpring - 6);
    Pt_pdf_Beauty->SetLineWidth(2);
    Pt_pdf_Beauty->SetLineStyle(9);
    // Pt_pdf_Beauty->Draw("SAME");

    auto Pt_pdf_Mixed = new TF1("Pt_pdf_Mixed", FuncPtMass, 0, 30, 4);
    Pt_pdf_Mixed->FixParameter(3, Pt_pdf_total->GetParameter(11));
    Pt_pdf_Mixed->FixParameter(0, B_DimuPtFromMixed->getVal());
    Pt_pdf_Mixed->FixParameter(1, n1_DimuPtFromMixed->getVal());
    Pt_pdf_Mixed->FixParameter(2, n2_DimuPtFromMixed->getVal());
    Pt_pdf_Mixed->SetLineColor(kCyan + 2);
    Pt_pdf_Mixed->SetLineWidth(2);
    Pt_pdf_Mixed->SetLineStyle(9);

    printf("--------------------------------------------------- \n");
    printf("From Pt fit \n");
    printf("N Charm: %0.0f \n", Pt_pdf_Charm->Integral(0, 30));
    printf("N Beauty: %0.0f \n", Pt_pdf_Beauty->Integral(0, 30));
    printf("N Mixed: %0.0f \n", Pt_pdf_Mixed->Integral(0, 30));
    printf("--------------------------------------------------- \n");

    gStyle->SetOptStat(0);
    TCanvas *c_Pt = new TCanvas("Pt_fit", "Pt_fit", 1000, 800);

    c_Pt->cd();
    histDimuPt->SetMaximum(1.8e+05);
    histDimuPt->SetMinimum(1.8e-02);
    TPad *Pt_pads = divide2(histDimuPt, Pt_pdf_total,0);

    Pt_pads[0].cd();
    TLegend *legend1 = new TLegend(0.175, 0.15, 0.35, 0.4);

    legend1->SetFillStyle(0);
    legend1->SetLineColor(kWhite);
    legend1->SetBorderSize(0);
    legend1->SetTextSize(0.0525);
    // legend1->SetHeader("Data");
    // legend1->SetTextAlign(11);
    legend1->AddEntry(histDimuPt, "Data");
    legend1->AddEntry(Pt_pdf_total, "Cocktail");
    legend1->AddEntry(Pt_pdf_Charm, "#mu#mu #leftarrow c#bar{c}");
    legend1->AddEntry(Pt_pdf_Beauty, "#mu#mu #leftarrow b#bar{b}");
    legend1->AddEntry(Pt_pdf_Mixed, "#mu#mu #leftarrow c,b");

    letexTitle->SetTextSize(0.055);
    letexTitle->DrawLatex(0.625, 0.825, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f", Pt_pdf_Charm->Integral(0, 30)));
    letexTitle->DrawLatex(0.625, 0.725, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f", Pt_pdf_Beauty->Integral(0, 30)));
    // letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.2e #pm %0.2e", fit_output[2]->getVal(), fit_output[2]->getError()));
    letexTitle->DrawLatex(0.625, 0.625, Form("#it{N}^{fixed}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f ", Pt_pdf_Mixed->Integral(0, 30)));

    Pt_pdf_Charm->Draw("SAME");
    Pt_pdf_Beauty->Draw("SAME");
    Pt_pdf_Mixed->Draw("SAME");
    legend1->Draw();
    c_Pt->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/pt_root_test_%0.0fptbin_%s_Opt.pdf",Binning_pt,opt.Data()));
    c_Pt->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/pt_root_test_%0.0fptbin_%s_Opt.png",Binning_pt,opt.Data()));
}

double FuncPtMass(double *x, double *par)
{
    return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
}

double total(double *x, double *par)
{
    return par[9] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2])) + par[10] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[3], par[4]), par[5])) + par[11] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[6], par[7]), par[8]));
}

TPad *divide2(TH1D *h_data, TF1 *model,int min)
{
    TH1D *draw_h_data = (TH1D *)h_data->Clone("draw_h_datav2");
    draw_h_data->SetTitle(" ");
    static TPad *pad[2];
    pad[0] = new TPad("pad1", "pad1", 0, 0.375, 1.0, 1.0);
    pad[0]->SetTicks();
    pad[0]->SetLogy(1);
    pad[0]->SetTopMargin(0.05);
    pad[0]->SetRightMargin(0.03);
    pad[0]->SetLeftMargin(0.14);
    pad[0]->SetBottomMargin(0.0);
    pad[0]->Draw();

    pad[1] = new TPad("pad2", "pad2", 0, 0.0, 1, 0.375);
    pad[1]->SetTicks();
    pad[1]->SetTopMargin(0.0);
    pad[1]->SetRightMargin(0.03);
    pad[1]->SetLeftMargin(0.14);
    pad[1]->SetBottomMargin(0.25);
    pad[1]->Draw();

    pad[0]->cd();
    // TH2D *h_grid = new TH2D("hgrid", " ", 30, 0, 30, 100, h_pt_total->GetMinimum() / 2., h_pt_total->GetMaximum() * 5.);
    draw_h_data->GetXaxis()->SetTitleOffset(0.85);
    draw_h_data->GetXaxis()->SetTitleSize(0.065);
    draw_h_data->GetXaxis()->SetLabelSize(0.055);

    draw_h_data->GetYaxis()->SetTitleOffset(1.);
    draw_h_data->GetYaxis()->SetNdivisions(505);
    draw_h_data->GetYaxis()->SetTitleSize(0.0575);
    draw_h_data->GetYaxis()->SetLabelSize(0.055);

    // draw_h_data->GetXaxis()->SetTitle(draw_h_data->GetXaxis()->GetTitle());
    // draw_h_data->GetYaxis()->SetTitle(draw_h_data->GetYaxis()->GetTitle());

    draw_h_data->Draw();
    model->Draw("SAME");

    // TLatex *letexTitle = new TLatex();
    // letexTitle->SetNDC();
    // letexTitle->SetTextFont(42);
    // letexTitle->SetTextSize(0.06);
    // // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    // letexTitle->DrawLatex(0.65, 0.875, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    // letexTitle->DrawLatex(0.65, 0.785, "LHC18p period");
    // letexTitle->DrawLatex(0.65, 0.695, Form("%0.0f < #it{m}_{#mu #mu} < %0.0f GeV/#it{c}^{2}", Low_Mass, High_Mass));

    pad[1]->cd();
    // TH2D *h_grid_ratio = new TH2D("hgrid_rateo", " ", 30, 0, 30, 100, 0.6, 2.2);
    TH1D *clone_h_data = (TH1D *)h_data->Clone("clone_h_datav2");

    // clone_h_data->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
    clone_h_data->GetYaxis()->SetTitle("Data/Cocktail");

    clone_h_data->GetYaxis()->CenterTitle();
    clone_h_data->GetYaxis()->SetNdivisions(504);
    clone_h_data->GetYaxis()->SetTitleSize(0.08);
    clone_h_data->GetYaxis()->SetTitleOffset(0.8);
    clone_h_data->GetYaxis()->SetLabelOffset(0.02);
    clone_h_data->GetYaxis()->SetLabelSize(0.1);

    clone_h_data->GetXaxis()->SetTitleSize(0.1);
    clone_h_data->GetXaxis()->SetTitleOffset(1.15);
    clone_h_data->GetXaxis()->SetLabelSize(0.1);

    clone_h_data->Sumw2();
    clone_h_data->Divide(model);
    clone_h_data->SetMinimum(-0.2);
    clone_h_data->SetMinimum(-0.2);

    // if (clone_h_data->GetMaximum() > 7.5)
    // {
    clone_h_data->SetMaximum(4.2);
    // }

    // if (clone_h_data->GetMinimum() < -7.5)
    // {
    //     clone_h_data->SetMinimum(-0.2);
    // }
    clone_h_data->SetTitle(" ");
    clone_h_data->Draw("EP");

    TLine *l = new TLine(min, 1.0, 30.0, 1.0);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->SetLineColor(kRed);
    l->Draw();

    return *pad;
}