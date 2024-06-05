#include "save_mc_output.h"
#include "TLorentzVector.h"
struct opt
{
    TString HF = "Charm";
    TString Pythia_Pros = "NonDiffr";
};

void cross_section()
{
    opt info;
    Double_t PYTHIA_CS;
    TString Pythia_filename;
    if (info.Pythia_Pros.Contains("NonDiffr"))
    {
        PYTHIA_CS = 56.42;
        Pythia_filename = "cs_Generator_SoftQCD_Def_pythia_sim_12345_DefaultBR_1000000.root";
    }
    else
    {
        PYTHIA_CS = 78.05;
        Pythia_filename = "";
    }
    TFile *fIn_hist = new TFile(Pythia_filename, "READ");
    TH1F *h_Nevents = (TH1F *)fIn_hist->Get("h_Nevents");
    fIn_hist->ls();
    TH1F *dsigma_SP_DiCharm = (TH1F *)fIn_hist->Get("DiCharm_rapidity");
    dsigma_SP_DiCharm->GetXaxis()->SetRangeUser(-4, -2.5);
    dsigma_SP_DiCharm->Draw();
    Double_t Charm_FWD_CS = (dsigma_SP_DiCharm->Integral() / h_Nevents->GetBinContent(2)) * (56.42 / 1.5);
    TH1F *dsigma_SP_DiBeauty = (TH1F *)fIn_hist->Get("DiBeauty_rapidity");
    dsigma_SP_DiBeauty->GetXaxis()->SetRangeUser(-4, -2.5);
    Double_t Beauty_FWD_CS = (dsigma_SP_DiBeauty->Integral() / h_Nevents->GetBinContent(2)) * (56.42 / 1.5);

    cout << "DIQUARK" << endl;
    cout << "Entries CHARM cs: " << dsigma_SP_DiCharm->Integral() << endl;
    cout << "Entries Beauty cs: " << dsigma_SP_DiBeauty->Integral() << endl;

    cout << "FWD CHARM cs: " << Charm_FWD_CS << endl;
    cout << "FWD Beauty cs: " << Beauty_FWD_CS << endl;
    cout << "Counting" << endl;
    TH1F *h_NCharm_event_fwd = (TH1F *)fIn_hist->Get("h_NCharm_event_fwd");

    TH1F *h_NBeauty_event_fwd = (TH1F *)fIn_hist->Get("h_NBeauty_event_fwd");
    Double_t NCharm_fwd = 0;
    Double_t NBeauty_fwd = 0;

    for (Int_t bin = 1; bin <= h_NCharm_event_fwd->GetNbinsX(); bin++)
    {
        NCharm_fwd = (NCharm_fwd + (0.5) * h_NCharm_event_fwd->GetBinContent(bin) * (bin - 1));
        NBeauty_fwd = (NBeauty_fwd + (0.5) * h_NBeauty_event_fwd->GetBinContent(bin) * (bin - 1));
    }
    cout << "Entries CHARM cs: " << NCharm_fwd << endl;
    cout << "Entries Beauty cs: " << NBeauty_fwd << endl;
    Double_t OLD_Charm_FWD_CS = (NCharm_fwd / h_Nevents->GetBinContent(2)) * (56.42 / 1.5);
    Double_t OLD_Beauty_FWD_CS = (NBeauty_fwd / h_Nevents->GetBinContent(2)) * (56.42 / 1.5);
    cout << "OLD FWD CHARM cs: " << OLD_Charm_FWD_CS << endl;
    cout << "OLD FWD Beauty cs: " << OLD_Beauty_FWD_CS << endl;

    //     TH1F *dsigma_SP_DiCharm = (TH1F *)fIn_hist_Charm->Get("dsigma_SP_DiCharm");
    // TH1F *dsigma_SP_DiBeauty = (TH1F *)fIn_hist_Beauty->Get("dsigma_SP_DiBeauty");
}

void b_c_ratio()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    opt info;
    // TFile *fIn_hist = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_inel_LFoff_Def_pythia_sim_2411_DefaultBR_output_Hist_100000.root", "READ");
    TString Pythia_filename;
    Double_t PYTHIA_CS = 999;
    if (info.Pythia_Pros.Contains("NonDiffr"))
    {
        PYTHIA_CS = 56.42;
        Pythia_filename = "cs_Generator_SoftQCD_Def_pythia_sim_12345_DefaultBR_100000.root";
    }
    else
    {
        PYTHIA_CS = 78.05;
        Pythia_filename = "";
    }

    TFile *fIn_hist_Charm = new TFile(Form("Charm_CS_%s.root", info.Pythia_Pros.Data()), "READ");
    TH1F *dsigma_Total_Charm = (TH1F *)fIn_hist_Charm->Get("dsigma_Total_Charm");
    TH1F *dsigma_SP_Charm = (TH1F *)fIn_hist_Charm->Get("dsigma_SP_Charm");
    TH1F *dsigma_SP_DiCharm = (TH1F *)fIn_hist_Charm->Get("dsigma_SP_DiCharm");

    TFile *fIn_hist_Beauty = new TFile(Form("Beauty_CS_%s.root", info.Pythia_Pros.Data()), "READ");
    fIn_hist_Beauty->ls();
    TH1F *dsigma_Total_Beauty = (TH1F *)fIn_hist_Beauty->Get("dsigma_Total_Beauty");
    TH1F *dsigma_SP_Beauty = (TH1F *)fIn_hist_Beauty->Get("dsigma_SP_Beauty");
    TH1F *dsigma_SP_DiBeauty = (TH1F *)fIn_hist_Beauty->Get("dsigma_SP_DiBeauty");

    TH1F *Beauty_Charm_Total_ratio = (TH1F *)dsigma_Total_Beauty->Clone("Beauty_Charm_Total_ratio");
    Beauty_Charm_Total_ratio->Divide(dsigma_Total_Charm);

    TH1F *Beauty_Charm_SP_ratio = (TH1F *)dsigma_SP_Beauty->Clone("Beauty_Charm_SP_ratio");
    Beauty_Charm_SP_ratio->Divide(dsigma_SP_Charm);

    TH1F *DiBeauty_DiCharm_SP_ratio = (TH1F *)dsigma_SP_DiBeauty->Clone("DiBeauty_DiCharm_SP_ratio");
    DiBeauty_DiCharm_SP_ratio->Divide(dsigma_SP_DiCharm);

    vector<double> vecY_Charm, vecL_Y_Charm, vecH_Y_Charm;
    vector<double> vec_central_Charm, vec_min_central_Charm, vec_max_central_Charm;
    vector<double> low_central_error_Charm, high_central_error_Charm;

    vector<double> vec_min_mass_Charm, vec_max_mass_Charm;
    vector<double> low_mass_error_Charm, high_mass_error_Charm;

    vector<double> vec_min_scale_Charm, vec_max_scale_Charm;
    vector<double> low_scale_error_Charm, high_scale_error_Charm;

    vector<double> vec_min_pdf_Charm, vec_max_pdf_Charm;
    vector<double> low_pdf_error_Charm, high_pdf_error_Charm;

    double Y_Charm, L_Y_Charm, H_Y_Charm, central_Charm, min_central_Charm, max_central_Charm, min_scale_Charm, max_scale_Charm, min_mass_Charm, max_mass_Charm, min_pdf_Charm, max_pdf_Charm;

    ifstream inputFile_Charm("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/FONLL_cc_Pt_0_30_CTEQ6.txt");

    while (inputFile_Charm >> Y_Charm >> L_Y_Charm >> H_Y_Charm >> central_Charm >> min_central_Charm >> max_central_Charm >> min_scale_Charm >> max_scale_Charm >> min_mass_Charm >> max_mass_Charm >> min_pdf_Charm >> max_pdf_Charm)
    {
        vecY_Charm.push_back(Y_Charm);
        vecL_Y_Charm.push_back(L_Y_Charm);
        vecH_Y_Charm.push_back(H_Y_Charm);

        vec_central_Charm.push_back(central_Charm);
        vec_min_central_Charm.push_back(min_central_Charm);
        vec_max_central_Charm.push_back(max_central_Charm);

        vec_min_scale_Charm.push_back(min_scale_Charm);
        vec_max_scale_Charm.push_back(max_scale_Charm);

        vec_min_mass_Charm.push_back(min_mass_Charm);
        vec_max_mass_Charm.push_back(max_mass_Charm);

        vec_min_pdf_Charm.push_back(min_pdf_Charm);
        vec_max_pdf_Charm.push_back(max_pdf_Charm);
    }

    vector<double> vecY_Beauty, vecL_Y_Beauty, vecH_Y_Beauty;
    vector<double> vec_central_Beauty, vec_min_central_Beauty, vec_max_central_Beauty;
    vector<double> low_central_error_Beauty, high_central_error_Beauty;

    vector<double> vec_min_mass_Beauty, vec_max_mass_Beauty;
    vector<double> low_mass_error_Beauty, high_mass_error_Beauty;

    vector<double> vec_min_scale_Beauty, vec_max_scale_Beauty;
    vector<double> low_scale_error_Beauty, high_scale_error_Beauty;

    vector<double> vec_min_pdf_Beauty, vec_max_pdf_Beauty;
    vector<double> low_pdf_error_Beauty, high_pdf_error_Beauty;

    double Y_Beauty, L_Y_Beauty, H_Y_Beauty, central_Beauty, min_central_Beauty, max_central_Beauty, min_scale_Beauty, max_scale_Beauty, min_mass_Beauty, max_mass_Beauty, min_pdf_Beauty, max_pdf_Beauty;

    ifstream inputFile_Beauty("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/FONLL_bb_Pt_0_30_CTEQ6.txt");

    while (inputFile_Beauty >> Y_Beauty >> L_Y_Beauty >> H_Y_Beauty >> central_Beauty >> min_central_Beauty >> max_central_Beauty >> min_scale_Beauty >> max_scale_Beauty >> min_mass_Beauty >> max_mass_Beauty >> min_pdf_Beauty >> max_pdf_Beauty)
    {
        vecY_Beauty.push_back(Y_Beauty);
        vecL_Y_Beauty.push_back(L_Y_Beauty);
        vecH_Y_Beauty.push_back(H_Y_Beauty);

        vec_central_Beauty.push_back(central_Beauty);
        vec_min_central_Beauty.push_back(min_central_Beauty);
        vec_max_central_Beauty.push_back(max_central_Beauty);

        vec_min_scale_Beauty.push_back(min_scale_Beauty);
        vec_max_scale_Beauty.push_back(max_scale_Beauty);

        vec_min_mass_Beauty.push_back(min_mass_Beauty);
        vec_max_mass_Beauty.push_back(max_mass_Beauty);

        vec_min_pdf_Beauty.push_back(min_pdf_Beauty);
        vec_max_pdf_Beauty.push_back(max_pdf_Beauty);
    }

    double vec_central_Beauty_Charm[vec_central_Charm.size()];
    double central[vec_central_Charm.size()];
    for (size_t i(0); i < vec_central_Charm.size(); i++)
    {
        vec_central_Beauty_Charm[i] = vec_central_Beauty[i] / vec_central_Charm[i];
        central[i] = vecY_Beauty[i];
    }

    TGraph *FONLL_bb_cs_NNPDF_minmaxerror = new TGraph(vec_central_Charm.size(), central, vec_central_Beauty_Charm);
    FONLL_bb_cs_NNPDF_minmaxerror->SetLineWidth(2);
    FONLL_bb_cs_NNPDF_minmaxerror->SetLineColor(kOrange + 2);

    TCanvas *canvas = canvas_for_prel();
    canvas->cd();

    Beauty_Charm_Total_ratio->GetYaxis()->SetRangeUser(8e-3, Beauty_Charm_Total_ratio->GetMaximum() * 12);
    Beauty_Charm_Total_ratio->GetYaxis()->SetTitle("Beauty/Charm");
    Beauty_Charm_Total_ratio->GetXaxis()->SetTitle("#it{y}");

    Beauty_Charm_Total_ratio->Draw();
    Beauty_Charm_SP_ratio->Draw("same");
    DiBeauty_DiCharm_SP_ratio->Draw("same");
    FONLL_bb_cs_NNPDF_minmaxerror->Draw("LSAME");

    TLegend *Legend_bb_cs_NNPDF_FONLL = new TLegend(0.175, 0.6, 0.575, 0.875, " ", "brNDC");
    Legend_bb_cs_NNPDF_FONLL->SetBorderSize(0);
    Legend_bb_cs_NNPDF_FONLL->SetFillColor(10);
    Legend_bb_cs_NNPDF_FONLL->SetFillStyle(1);
    Legend_bb_cs_NNPDF_FONLL->SetLineStyle(0);
    Legend_bb_cs_NNPDF_FONLL->SetLineColor(0);
    Legend_bb_cs_NNPDF_FONLL->SetTextFont(42);
    Legend_bb_cs_NNPDF_FONLL->SetTextSize(0.04);
    Legend_bb_cs_NNPDF_FONLL->AddEntry(Beauty_Charm_Total_ratio, "Total");
    Legend_bb_cs_NNPDF_FONLL->AddEntry(Beauty_Charm_SP_ratio, "Single Pair");
    Legend_bb_cs_NNPDF_FONLL->AddEntry(DiBeauty_DiCharm_SP_ratio, "Diquark");
    Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_minmaxerror, "FONLL", "L");
    Legend_bb_cs_NNPDF_FONLL->Draw();
}

void generator_cs()
{
    opt info;
    // TFile *fIn_hist = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_inel_LFoff_Def_pythia_sim_2411_DefaultBR_output_Hist_100000.root", "READ");
    TString Pythia_filename;
    Double_t PYTHIA_CS = 999;
    if (info.Pythia_Pros.Contains("NonDiffr"))
    {
        PYTHIA_CS = 56.42;
        Pythia_filename = "cs_Generator_SoftQCD_Def_pythia_sim_12345_DefaultBR_1000000.root";
    }
    else
    {
        PYTHIA_CS = 78.05;
        Pythia_filename = "";
    }

    TFile *fIn_hist = new TFile(Pythia_filename, "READ");
    TH1F *h_Nevents = (TH1F *)fIn_hist->Get("h_Nevents");
    TH1F *h_NQuark_event = (TH1F *)fIn_hist->Get(Form("h_N%s_event", info.HF.Data()));

    Double_t SP_cc = (h_NQuark_event->GetBinContent(3) / h_Nevents->GetBinContent(2)) * PYTHIA_CS * 1000; // in mub
    cout << Form("Single %s cs: ", info.HF.Data()) << SP_cc << " (mub)" << endl;

    Double_t Total_cc_Quark = 0;
    Double_t n_pairs_counting = 0;
    for (Int_t bin = 1; bin < h_NQuark_event->GetNbinsX(); bin++)
        n_pairs_counting = n_pairs_counting + ((bin - 1) * h_NQuark_event->GetBinContent(bin) / 2.);

    cout << "n_pairs_counting: " << n_pairs_counting << endl;
    Total_cc_Quark = n_pairs_counting / (h_Nevents->GetBinContent(2)) * PYTHIA_CS * 1000;
    cout << Form("Total_cc_%s: ", info.HF.Data()) << Total_cc_Quark << " (mub)" << endl;

    // TH2F *h_PtY_Charm_quark = (TH2F *)fIn_hist->Get("HF_quarks/h_PtY_Charm_quark");
    TH1F *dsigma_SP_Quark = (TH1F *)fIn_hist->Get(Form("%s_rapidity", info.HF.Data()));
    dsigma_SP_Quark->SetName(Form("dsigma_SP_%s", info.HF.Data()));
    TH1F *dsigma_Total_Quark = (TH1F *)fIn_hist->Get(Form("%s_rapidity", info.HF.Data()));
    dsigma_Total_Quark->SetName(Form("dsigma_Total_%s", info.HF.Data()));

    TH1F *dsigma_SP_DiQuark = (TH1F *)fIn_hist->Get(Form("Di%s_rapidity", info.HF.Data()));
    dsigma_SP_DiQuark->SetName(Form("dsigma_SP_Di%s", info.HF.Data()));

    Double_t n_diquarks = dsigma_SP_DiQuark->GetEntries();
    cout << "n_diquarks: " << n_diquarks << endl;
    cout << Form("Di%s cs: ", info.HF.Data()) << (n_diquarks / h_Nevents->GetBinContent(2)) * PYTHIA_CS * 1000 << " (mub)" << endl;

    hist1D_graphic_opt(dsigma_SP_DiQuark, kTRUE, 1, 24, kGreen, (1 / h_Nevents->GetBinContent(2)) * PYTHIA_CS * 1000);

    hist1D_graphic_opt(dsigma_SP_Quark, kTRUE, 1, 24, kBlack, SP_cc / dsigma_SP_Quark->Integral());

    hist1D_graphic_opt(dsigma_Total_Quark, kTRUE, 1, 20, kRed, Total_cc_Quark / dsigma_Total_Quark->Integral());

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

    TString FONLL_filename;
    Int_t HF_Selector = 999; // 0 for CHARM, 1 for BEAUTY

    if (info.HF.Contains("Charm"))
        FONLL_filename.Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/FONLL_cc_Pt_0_30_CTEQ6.txt");

    else if (info.HF.Contains("Beauty"))
        FONLL_filename.Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/FONLL_bb_Pt_0_30_CTEQ6.txt");

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

    TMultiGraph *bb_bar_cs_NNPDF = new TMultiGraph();
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_minmaxerror, "CX");
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_scale_error);
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_mass_error);
    bb_bar_cs_NNPDF->Add(FONLL_bb_cs_NNPDF_pdf_error);

    if (info.HF.Contains("Beauty"))
    {
        bb_bar_cs_NNPDF->GetYaxis()->SetRangeUser(6.5e-1, 8e+2);
        bb_bar_cs_NNPDF->GetYaxis()->SetTitle("d#sigma_{b#bar{b}} / d#it{y}(#mub)");
    }
    else if (info.HF.Contains("Charm"))
    {
        bb_bar_cs_NNPDF->GetYaxis()->SetRangeUser(7.2e-0, 6.5e+4);
        bb_bar_cs_NNPDF->GetYaxis()->SetTitle("d#sigma_{c#bar{c}} / d#it{y} (#mub)");
    }
    TCanvas *canvas = canvas_for_prel();
    canvas->cd();

    bb_bar_cs_NNPDF->GetXaxis()->SetTitle("#it{y}");
    bb_bar_cs_NNPDF->Draw("A3");
    dsigma_Total_Quark->Draw("same");
    dsigma_SP_Quark->Draw("same");
    dsigma_SP_DiQuark->Draw("same");

    TLegend *Legend_bb_cs_NNPDF_Meas = new TLegend(0.16, 0.175, 0.85, 0.475, " ", "brNDC");
    Legend_bb_cs_NNPDF_Meas->SetBorderSize(0);
    Legend_bb_cs_NNPDF_Meas->SetFillColor(10);
    Legend_bb_cs_NNPDF_Meas->SetFillStyle(1);
    Legend_bb_cs_NNPDF_Meas->SetLineStyle(0);
    Legend_bb_cs_NNPDF_Meas->SetLineColor(0);
    Legend_bb_cs_NNPDF_Meas->SetTextFont(42);
    Legend_bb_cs_NNPDF_Meas->SetTextSize(0.04);
    Legend_bb_cs_NNPDF_Meas->AddEntry(dsigma_Total_Quark, Form("Total %s cs: %0.2e #mub", info.HF.Data(), Total_cc_Quark));
    Legend_bb_cs_NNPDF_Meas->AddEntry(dsigma_SP_Quark, Form("Single pair %s cs: %0.2e #mub", info.HF.Data(), SP_cc));
    Legend_bb_cs_NNPDF_Meas->AddEntry(dsigma_SP_DiQuark, Form("Di%s cs: %0.2e #mub", info.HF.Data(), (dsigma_SP_DiQuark->GetEntries() / h_Nevents->GetBinContent(2)) * PYTHIA_CS * 1000));
    Legend_bb_cs_NNPDF_Meas->Draw("same");

    TLegend *Legend_bb_cs_NNPDF_FONLL = new TLegend(0.725, 0.18, 0.875, 0.325, " ", "brNDC");
    Legend_bb_cs_NNPDF_FONLL->SetBorderSize(0);
    Legend_bb_cs_NNPDF_FONLL->SetFillColor(10);
    Legend_bb_cs_NNPDF_FONLL->SetFillStyle(1);
    Legend_bb_cs_NNPDF_FONLL->SetLineStyle(0);
    Legend_bb_cs_NNPDF_FONLL->SetLineColor(0);
    Legend_bb_cs_NNPDF_FONLL->SetTextFont(42);
    Legend_bb_cs_NNPDF_FONLL->SetTextSize(0.04);
    Legend_bb_cs_NNPDF_FONLL->AddEntry(FONLL_bb_cs_NNPDF_minmaxerror, "FONLL", "F");
    Legend_bb_cs_NNPDF_FONLL->Draw();

    TLatex letexTitle;
    letexTitle.SetTextSize(0.055);
    letexTitle.SetNDC();
    letexTitle.SetTextFont(42);
    letexTitle.DrawLatex(0.2, 0.82, Form("PYTHIA sim %s, pp#sqrt{#it{s}} = 13 TeV", info.Pythia_Pros.Data()));
    letexTitle.SetTextSize(0.0375);
    letexTitle.DrawLatex(0.2, 0.74, "FONLL CTEQ6");
    canvas->SaveAs(Form("images/cs_generator_%s.pdf", info.HF.Data()));

    TFile *fOut = new TFile(Form("%s_CS_%s.root", info.HF.Data(), Pythia_filename.Data()), "RECREATE");
    dsigma_Total_Quark->Write(0, 2, 0);
    dsigma_SP_Quark->Write(0, 2, 0);
    dsigma_SP_DiQuark->Write(0, 2, 0);
}

void rapidity()
{
    opt info;
    TH1F *h_Nevents = new TH1F("h_Nevents", "h_Nevents", 1, 0, 1);
    TH1F *Charm_rapidity = new TH1F("Charm_rapidity", "Charm_rapidity", 200, -10, 10);
    TH2F *NDiCharm_pt = new TH2F("NDiCharm_pt", "NDiCharm_pt; N_{cc} x ev.; #it{p}_{T} (GeV/#it{c})", 10, 0, 10, 300, 0, 30);
    TH2F *NDiCharm_y = new TH2F("NDiCharm_y", "NDiCharm_y; N_{bb} x ev.; #it{y}", 10, 0, 10, 200, -10, 10);
    TH1F *Beauty_rapidity = new TH1F("Beauty_rapidity", "Beauty_rapidity", 200, -10, 10);

    TH1F *DiCharm_rapidity = new TH1F("DiCharm_rapidity", "DiCharm_rapidity", 200, -10, 10);
    TH2F *NDiBeauty_pt = new TH2F("NDiBeauty_pt", "NDiBeauty_pt; N_{bb} x ev.; #it{p}_{T} (GeV/#it{c})", 10, 0, 10, 300, 0, 30);
    TH2F *NDiBeauty_y = new TH2F("NDiBeauty_y", "NDiBeauty_y; N_{bb} x ev.; #it{y}", 10, 0, 10, 200, -10, 10);
    TH1F *DiBeauty_rapidity = new TH1F("DiBeauty_rapidity", "DiBeauty_rapidity", 200, -10, 10);

    TString Generator = "stand_HF_Pythia";
    Set_Histograms(Generator);
    TString Pythia_filename;
    if (info.Pythia_Pros.Contains("NonDiffr"))
        Pythia_filename = "SoftQCD_Def_pythia_sim_12345_DefaultBR_1000000.root";
    else
        Pythia_filename = "";

    TChain *input_tree = Importing_Tree("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/sim", Pythia_filename, Generator);
    input_tree->ls();
    Int_t total_entries = input_tree->GetEntries();
    for (Int_t i_Event = 0; i_Event < total_entries; i_Event++)
    {
        h_Nevents->Fill(1);
        input_tree->GetEntry(i_Event);
        if (i_Event % (Int_t)(total_entries * 0.05) == 0)
            progress_status(i_Event, total_entries);
        Int_t N_Charm_event = 0;
        Int_t N_Charm_event_fwd = 0;

        Int_t N_Beauty_event = 0;
        Int_t N_Beauty_event_fwd = 0;

        // if (N_HFquarks_gen > 0)
        //     printf("N_HFquarks_gen %d\n", N_HFquarks_gen);
        Int_t n_pair = 0;

        vector<double> DiCharm_pt;
        vector<double> DiCharm_y;

        vector<double> DiBeauty_pt;
        vector<double> DiBeauty_y;

        for (Int_t i_N_HFquarks_gen = 0; i_N_HFquarks_gen < N_HFquarks_gen; i_N_HFquarks_gen++)
        {
            Int_t PDG_HFquark = PDG_HFquark_gen[i_N_HFquarks_gen];
            Double_t Pt_HFquark = Pt_HFquark_gen[i_N_HFquarks_gen];
            Double_t Y_HFquark = Y_HFquark_gen[i_N_HFquarks_gen];
            Double_t Px_HFquark = Px_HFquark_gen[i_N_HFquarks_gen];
            Double_t Py_HFquark = Py_HFquark_gen[i_N_HFquarks_gen];
            Double_t Pz_HFquark = Pz_HFquark_gen[i_N_HFquarks_gen];
            Double_t E_HFquark = E_HFquark_gen[i_N_HFquarks_gen];
            Int_t Quark1_Mother1 = Index_HFquark_gen_mother1[i_N_HFquarks_gen];
            Int_t Quark1_Mother2 = Index_HFquark_gen_mother2[i_N_HFquarks_gen];
            // printf("PDG_HFquark %d\n", PDG_HFquark);
            // printf("Px_HFquark %0.1f\n", Px_HFquark);
            // printf("Py_HFquark %0.1f\n", Py_HFquark);
            // printf("Pz_HFquark %0.1f\n", Pz_HFquark);
            // printf("E_HFquark %0.1f\n", E_HFquark);
            // if (Y_HFquark < -5 || Y_HFquark > 5)
            //     continue;
            TLorentzVector Particle1(Px_HFquark, Py_HFquark, Pz_HFquark, E_HFquark);
            if (Particle1.Rapidity() != Y_HFquark)
            {
                printf("Y_HFquark %0.1f\n", Y_HFquark);
                printf("Particle1.Rapidity() %0.1f\n", Particle1.Rapidity());
            }

            if (TMath::Abs(PDG_HFquark) == 4)
            {
                h_PtY_Charm_quark->Fill(Pt_HFquark, Y_HFquark);
                N_Charm_event++;
                Charm_rapidity->Fill(Y_HFquark);

                if (Y_HFquark > -4.0 && Y_HFquark < -2.5)
                    N_Charm_event_fwd++;
            }
            else if (TMath::Abs(PDG_HFquark) == 5)
            {
                h_PtY_Beauty_quark->Fill(Pt_HFquark, Y_HFquark);
                N_Beauty_event++;
                Beauty_rapidity->Fill(Y_HFquark);

                if (Y_HFquark > -4.0 && Y_HFquark < -2.5)
                    N_Beauty_event_fwd++;
            }
            for (Int_t j_N_HFquarks_gen = i_N_HFquarks_gen + 1; j_N_HFquarks_gen < N_HFquarks_gen; j_N_HFquarks_gen++)
            {
                Int_t PDG_HFquark2 = PDG_HFquark_gen[j_N_HFquarks_gen];
                Double_t Pt_HFquark2 = Pt_HFquark_gen[j_N_HFquarks_gen];
                Double_t Y_HFquark2 = Y_HFquark_gen[j_N_HFquarks_gen];
                Double_t Px_HFquark2 = Px_HFquark_gen[j_N_HFquarks_gen];
                Double_t Py_HFquark2 = Py_HFquark_gen[j_N_HFquarks_gen];
                Double_t Pz_HFquark2 = Pz_HFquark_gen[j_N_HFquarks_gen];
                Double_t E_HFquark2 = E_HFquark_gen[j_N_HFquarks_gen];
                Int_t Quark2_Mother1 = Index_HFquark_gen_mother1[j_N_HFquarks_gen];
                Int_t Quark2_Mother2 = Index_HFquark_gen_mother2[j_N_HFquarks_gen];

                TLorentzVector Particle2(Px_HFquark2, Py_HFquark2, Pz_HFquark2, E_HFquark2);

                if (PDG_HFquark + PDG_HFquark2 != 0)
                    continue;

                // printf("Quark1-> Mother1 %d Mother2 %d || Quark2-> Mother1 %d Mother2 %d\n", Quark1_Mother1, Quark1_Mother2, Quark2_Mother1, Quark2_Mother2);
                if ((Quark1_Mother1 - Quark2_Mother1 != 0) && N_HFquarks_gen > 2)
                    continue;
                TLorentzVector DiQuark = Particle1 + Particle2;

                if (TMath::Abs(PDG_HFquark) == 4 && TMath::Abs(PDG_HFquark2) == 4)
                {
                    DiCharm_rapidity->Fill(DiQuark.Rapidity());
                    DiCharm_pt.push_back(DiQuark.Pt());
                    DiCharm_y.push_back(DiQuark.Rapidity());
                }

                if (TMath::Abs(PDG_HFquark) == 5 && TMath::Abs(PDG_HFquark2) == 5){
                    DiBeauty_rapidity->Fill(DiQuark.Rapidity());
                    DiBeauty_pt.push_back(DiQuark.Pt());
                    DiBeauty_y.push_back(DiQuark.Rapidity());
                }
                n_pair++;
            }
        }
        
        for (Int_t i_DiCharm_x_ev = 0; i_DiCharm_x_ev < DiCharm_pt.size(); i_DiCharm_x_ev++)
        {
            NDiCharm_pt->Fill(DiCharm_pt.size(), DiCharm_pt[i_DiCharm_x_ev]);
            NDiCharm_y->Fill(DiCharm_pt.size(), DiCharm_y[i_DiCharm_x_ev]);
        }

        for (Int_t i_DiBeauty_x_ev = 0; i_DiBeauty_x_ev < DiBeauty_pt.size(); i_DiBeauty_x_ev++)
        {
            NDiBeauty_pt->Fill(DiBeauty_pt.size(), DiBeauty_pt[i_DiBeauty_x_ev]);
            NDiBeauty_y->Fill(DiBeauty_pt.size(), DiBeauty_y[i_DiBeauty_x_ev]);
        }
        DiCharm_pt.clear();
        DiCharm_y.clear();

        DiBeauty_pt.clear();
        DiBeauty_y.clear();
        // cout << n_pair << endl;
        if (n_pair == 0 && N_HFquarks_gen != 0)
            return;
        // if (N_Charm_event > 0)
        //     printf("N_Charm_event %d\n", N_Charm_event);
        // if (N_Beauty_event > 0)
        //     printf("N_Beauty_event %d\n", N_Beauty_event);
        h_NCharm_event->Fill(N_Charm_event);
        h_NCharm_event_fwd->Fill(N_Charm_event_fwd);

        h_NBeauty_event->Fill(N_Beauty_event);
        h_NBeauty_event_fwd->Fill(N_Beauty_event_fwd);
    }
    TFile *fOut = new TFile(Form("test_cs_Generator_%s", Pythia_filename.Data()), "RECREATE");
    h_Nevents->Write();
    h_NCharm_event->Write();
    h_NCharm_event_fwd->Write();
    h_NBeauty_event->Write();
    h_NBeauty_event_fwd->Write();
    DiCharm_rapidity->Write();
    NDiCharm_pt->Write();
    NDiCharm_y->Write();
    Charm_rapidity->Write();
    DiBeauty_rapidity->Write();
    Beauty_rapidity->Write();
    NDiBeauty_pt->Write();
    NDiBeauty_y->Write();

    fOut->Close();
    h_NCharm_event->Draw();
    new TCanvas();
    h_NBeauty_event->Draw();
    new TCanvas();
    DiBeauty_rapidity->SetLineColor(kRed);
    DiCharm_rapidity->SetLineColor(kRed);
    DiCharm_rapidity->Draw();
    Charm_rapidity->Draw("sameS");
    gPad->BuildLegend();
    gPad->SetLogy();
    new TCanvas();
    DiBeauty_rapidity->Draw();
    Beauty_rapidity->Draw("sameS");
    gPad->BuildLegend();
    gPad->SetLogy();
    // h_NCharm_event->Draw("same");
}