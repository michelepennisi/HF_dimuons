void LoadStyle();
TCanvas *compare_mode(TH1D *h1, TH1D *h2, TH1D *h3);

void mode_pythia()
{
    const Int_t n_Mode_studied = 3;
    const Int_t n_Hadron_studied = 4;
    const Int_t n_Rapidity_studied = 6;

    TString name_File_studied[n_Hadron_studied];
    TString name_Hadron_studied[n_Hadron_studied];
    TString name_rapidity_studied[n_Rapidity_studied];

    name_File_studied[0].Form("Def");
    name_File_studied[1].Form("Config+Atlas");
    name_File_studied[2].Form("Mode2");

    name_Hadron_studied[0].Form("Dzero");
    name_Hadron_studied[1].Form("Dplus");
    name_Hadron_studied[2].Form("Dstrange");
    name_Hadron_studied[3].Form("Lambda");

    name_rapidity_studied[0].Form("MidY");
    name_rapidity_studied[1].Form("FwdY1");
    name_rapidity_studied[2].Form("FwdY2");
    name_rapidity_studied[3].Form("FwdY3");
    name_rapidity_studied[4].Form("FwdY4");
    name_rapidity_studied[5].Form("FwdY5");

    TFile *fIn_Mode[n_Mode_studied];

    TH1D *h_ptHadron_prompt[n_Mode_studied][n_Hadron_studied][n_Rapidity_studied];
    TH1D *h_yHadron_prompt[n_Mode_studied][n_Hadron_studied][n_Rapidity_studied];
    TH1D *h_pdgHadron_prompt[n_Mode_studied][n_Hadron_studied][n_Rapidity_studied];

    for (Int_t i_Mode = 0; i_Mode < n_Mode_studied; i_Mode++)
    {
        fIn_Mode[i_Mode] = new TFile(Form("sim/Hist_HFHadronpythia_sim_SoftQCD_%s_1000000_2710_DefaultBR.root", name_File_studied[i_Mode].Data()), "READ");
        fIn_Mode[i_Mode]->ls();
        for (Int_t i_Hadron = 0; i_Hadron < n_Hadron_studied; i_Hadron++)
        {
            for (Int_t i_Rapidity = 0; i_Rapidity < n_Rapidity_studied; i_Rapidity++)
            {
                h_ptHadron_prompt[i_Mode][i_Hadron][i_Rapidity] = (TH1D *)fIn_Mode[i_Mode]->Get(Form("%s/h_pt%s_prompt_%s", name_rapidity_studied[i_Rapidity].Data(), name_Hadron_studied[i_Hadron].Data(), name_rapidity_studied[i_Rapidity].Data()));
                // printf("%s_%s", h_ptHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->GetName(), name_File_studied[i_Mode].Data());
                h_ptHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->SetName(Form("%s_%s", h_ptHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->GetName(), name_File_studied[i_Mode].Data()));
                h_ptHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->SetTitle(Form("%s_%s", h_ptHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->GetTitle(), name_File_studied[i_Mode].Data()));

                h_yHadron_prompt[i_Mode][i_Hadron][i_Rapidity] = (TH1D *)fIn_Mode[i_Mode]->Get(Form("%s/h_y%s_prompt_%s", name_rapidity_studied[i_Rapidity].Data(), name_Hadron_studied[i_Hadron].Data(), name_rapidity_studied[i_Rapidity].Data()));
                h_yHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->SetName(Form("%s_%s", h_yHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->GetName(), name_File_studied[i_Mode].Data()));
                h_yHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->SetTitle(Form("%s_%s", h_yHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->GetTitle(), name_File_studied[i_Mode].Data()));

                h_pdgHadron_prompt[i_Mode][i_Hadron][i_Rapidity] = (TH1D *)fIn_Mode[i_Mode]->Get(Form("%s/h_pdg%s_prompt_%s", name_rapidity_studied[i_Rapidity].Data(), name_Hadron_studied[i_Hadron].Data(), name_rapidity_studied[i_Rapidity].Data()));
                h_pdgHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->SetName(Form("%s_%s", h_pdgHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->GetName(), name_File_studied[i_Mode].Data()));
                h_pdgHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->SetTitle(Form("%s_%s", h_pdgHadron_prompt[i_Mode][i_Hadron][i_Rapidity]->GetTitle(), name_File_studied[i_Mode].Data()));
            }
        }
    }

    // h_ptHadron_prompt[1][0][0]->Draw();

    TCanvas *c = compare_mode(h_ptHadron_prompt[0][0][0], h_ptHadron_prompt[1][0][0], h_ptHadron_prompt[2][0][0]);
}

TCanvas *compare_mode(TH1D *h1, TH1D *h2, TH1D *h3)
{
    LoadStyle();
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1200);
    canvas->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.45, 1.0, 1.0);
    pad1->SetTicks();
    pad1->SetLogy(1);
    pad1->SetTopMargin(0.05);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.14);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.25, 1, 0.45);
    pad2->SetTicks();
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.03);
    pad2->SetLeftMargin(0.14);
    pad2->SetBottomMargin(0.0);
    pad2->Draw();

    TPad *pad3 = new TPad("pad3", "pad3", 0, 0.0, 1, 0.25);
    pad3->SetTicks();
    pad3->SetTopMargin(0.0);
    pad3->SetRightMargin(0.03);
    pad3->SetLeftMargin(0.14);
    pad3->SetBottomMargin(0.25);
    pad3->Draw();

    pad1->cd();
    h1->SetLineColor(kBlack);
    h1->SetLineWidth(2.0);
    h1->SetMarkerColor(kBlack);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1.2);
    h1->Draw("PE");

    h2->SetLineColor(kRed);
    h2->SetMarkerColor(kRed);
    h2->SetMarkerStyle(24);
    h2->SetLineWidth(2.0);
    h2->SetMarkerSize(1.2);

    h2->Draw("PESAME");

    h3->SetLineColor(kBlue);
    h3->SetMarkerColor(kBlue);
    h3->SetMarkerStyle(24);
    h3->SetLineWidth(2.0);
    h3->SetMarkerSize(1.2);

    h3->Draw("PESAME");

    pad2->cd();

    TH1D *h2_clone = (TH1D *)h2->Clone("h2_clone");
    h2_clone->SetTitle("");
    h2_clone->GetYaxis()->CenterTitle();
    h2_clone->GetYaxis()->SetNdivisions(504);
    h2_clone->GetYaxis()->SetTitleSize(0.07);
    // h_grid_ratio->GetYaxis()->SetTitleOffset(0.8);
    h2_clone->GetYaxis()->SetLabelOffset(0.02);
    h2_clone->GetYaxis()->SetLabelSize(0.09);

    h2_clone->GetXaxis()->SetTitleSize(0.09);
    h2_clone->GetXaxis()->SetTitleOffset(1.1);
    h2_clone->GetXaxis()->SetLabelSize(0.09);
    h2_clone->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h2_clone->Divide(h1);
    h2_clone->SetMaximum(1.4);
    h2_clone->SetMinimum(0.2);
    h2_clone->Draw("PE");

    pad3->cd();
    TH1D *h3_clone = (TH1D *)h3->Clone("h3_clone");
    h3_clone->SetTitle("");
    h3_clone->GetYaxis()->CenterTitle();
    h3_clone->GetYaxis()->SetNdivisions(504);
    h3_clone->GetYaxis()->SetTitleSize(0.07);
    // h_grid_ratio->GetYaxis()->SetTitleOffset(0.8);
    h3_clone->GetYaxis()->SetLabelOffset(0.02);
    h3_clone->GetYaxis()->SetLabelSize(0.09);

    h3_clone->GetXaxis()->SetTitleSize(0.09);
    h3_clone->GetXaxis()->SetTitleOffset(1.1);
    h3_clone->GetXaxis()->SetLabelSize(0.09);
    h3_clone->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h3_clone->Divide(h1);

    h3_clone->SetMaximum(1.4);
    h3_clone->SetMinimum(0.2);
    h3_clone->Draw("PE");
    return canvas;
}

void LoadStyle()
{
    gStyle->SetImageScaling(3.);
    Int_t font = 42;
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
    gStyle->SetTextFont(font);
    gStyle->SetStatFontSize(0.05);
    gStyle->SetStatX(0.97);
    gStyle->SetStatY(0.98);
    gStyle->SetStatH(0.03);
    gStyle->SetStatW(0.3);
    gStyle->SetTickLength(0.02, "y");
    gStyle->SetEndErrorSize(3);
    gStyle->SetLabelSize(0.04, "xyz");
    gStyle->SetLabelFont(font, "xyz");
    gStyle->SetLabelOffset(0.01, "xyz");
    gStyle->SetTitleFont(font, "xyz");
    gStyle->SetTitleOffset(0.9, "x");
    gStyle->SetTitleOffset(1.02, "y");
    gStyle->SetTitleSize(0.04, "xyz");
    gStyle->SetMarkerSize(1.3);
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetHatchesSpacing(0.5);
}