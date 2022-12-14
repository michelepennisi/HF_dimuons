double FuncPtMass(double *x, double *par);
double total(double *x, double *par);

void test_root_fit()
{
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");

    TFile *fIn = new TFile(Form("root_files/pdfMC_unbinned.root"));

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

    TH2D *histDimuPtM_fromdata = (TH2D *)fIn_data->Get(Form("Dimuon/CMUL7/DQ_cut_match_LT_ULS/h_PtMdiMu_CMUL7_DQ_cut_match_LT_ULS"));
    histDimuPtM_fromdata->SetName("histDimuPtM_fromdata");

    TH1D *histDimuPt = (TH1D *)histDimuPtM_fromdata->ProjectionX();
    histDimuPt->Sumw2();
    TH1D *histDimuMass = (TH1D *)histDimuPtM_fromdata->ProjectionY();
    histDimuMass->Sumw2();

    histDimuMass->Rebin(1);
    histDimuMass->Scale(1., "width");

    histDimuMass->SetMarkerStyle(20);
    histDimuMass->SetMarkerColor(kBlack);

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
    mass_pdf_total->SetParameter(9, 500000);
    mass_pdf_total->SetParameter(10, 200000);
    mass_pdf_total->FixParameter(11, 316000);
    
    histDimuMass->Fit(mass_pdf_total, "R");

    auto mass_pdf_Charm = new TF1("mass_pdf_Charm", FuncPtMass, 4, 30, 4);
    mass_pdf_Charm->SetParameter(3, mass_pdf_total->GetParameter(9));
    mass_pdf_Charm->FixParameter(0, B_DimuMassFromCharm->getVal());
    mass_pdf_Charm->FixParameter(1, n1_DimuMassFromCharm->getVal());
    mass_pdf_Charm->FixParameter(2, n2_DimuMassFromCharm->getVal());
    mass_pdf_Charm->SetLineColor(kMagenta - 2);
    // mass_pdf_Charm->Draw();

    auto mass_pdf_Beauty = new TF1("mass_pdf_Beauty", FuncPtMass, 4, 30, 4);
    mass_pdf_Beauty->FixParameter(3, mass_pdf_total->GetParameter(10));
    mass_pdf_Beauty->FixParameter(0, B_DimuMassFromBeauty->getVal());
    mass_pdf_Beauty->FixParameter(1, n1_DimuMassFromBeauty->getVal());
    mass_pdf_Beauty->FixParameter(2, n2_DimuMassFromBeauty->getVal());
    mass_pdf_Beauty->SetLineColor(kSpring - 6);
    // mass_pdf_Beauty->Draw("SAME");

    auto mass_pdf_Mixed = new TF1("mass_pdf_Mixed", FuncPtMass, 4, 30, 4);
    mass_pdf_Mixed->FixParameter(3, mass_pdf_total->GetParameter(11));
    mass_pdf_Mixed->FixParameter(0, B_DimuMassFromMixed->getVal());
    mass_pdf_Mixed->FixParameter(1, n1_DimuMassFromMixed->getVal());
    mass_pdf_Mixed->FixParameter(2, n2_DimuMassFromMixed->getVal());
    mass_pdf_Mixed->SetLineColor(kCyan + 2);

    histDimuMass->Draw();
    mass_pdf_total->Draw("SAME");
    mass_pdf_Charm->Draw("SAME");
    mass_pdf_Beauty->Draw("SAME");
    mass_pdf_Mixed->Draw("SAME");

    printf("N Charm: %0.0f \n", mass_pdf_Charm->Integral(4, 30));
    printf("N Beauty: %0.0f \n", mass_pdf_Beauty->Integral(4, 30));
    printf("N Mixed: %0.0f \n", mass_pdf_Mixed->Integral(4, 30));
}

double FuncPtMass(double *x, double *par)
{
    return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
}

double total(double *x, double *par)
{
    return par[9] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2])) + par[10] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[3], par[4]), par[5])) + par[11] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[6], par[7]), par[8]));
}