#include "/home/michele_pennisi/cernbox/common_include.h"

double FuncPtMass(double *x, double *par);

void dimuon_shape(Int_t origin = 4)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TString Origin_powheg;
    TString Origin;
    TString Header;
    if (origin == 4)
    {
        Header = "#splitline{Gen. #mu^{#plus}#mu^{#minus} #leftarrow c}{No cut}";
        Origin_powheg = "charm";
        Origin = "Charm";
    }
    else if (origin == 5)
    {
        Header = "#splitline{Gen. #mu^{#plus}#mu^{#minus} #leftarrow b}{No cut}";

        Origin_powheg = "beauty";
        Origin = "Beauty";
    }

    TFile *fIn_Powheg = new TFile(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Powheg_Sim/powheg_%s/Version1/save_mc_output/powheg_%s_MC_output_Hist_merged.root", Origin_powheg.Data(), Origin_powheg.Data()), "READ");

    TH2F *h_PtM_Muon_Gen_Charm_POWHEG = (TH2F *)fIn_Powheg->Get(Form("DiMuon_Gen/h_PtM_DiMuon_Gen_%s", Origin.Data()));
    h_PtM_Muon_Gen_Charm_POWHEG->SetName("h_PtM_Muon_Gen_POWHEG");
    TH1F *h_Pt_Muon_Gen_Charm_POWHEG = (TH1F *)h_PtM_Muon_Gen_Charm_POWHEG->ProjectionX();
    hist1D_graphic_opt(h_Pt_Muon_Gen_Charm_POWHEG, 5, 20, kCyan + 1, 1. / h_Pt_Muon_Gen_Charm_POWHEG->GetEntries());
    h_Pt_Muon_Gen_Charm_POWHEG->SetName("h_Pt_Muon_Gen_Charm_POWHEG");
    h_Pt_Muon_Gen_Charm_POWHEG->SetTitle("Distr.");

    TH1F *h_M_Muon_Gen_Charm_POWHEG = (TH1F *)h_PtM_Muon_Gen_Charm_POWHEG->ProjectionY();
    hist1D_graphic_opt(h_M_Muon_Gen_Charm_POWHEG, 5, 20, kCyan + 1, 1. / h_Pt_Muon_Gen_Charm_POWHEG->GetEntries());
    h_M_Muon_Gen_Charm_POWHEG->SetName("h_Pt_Muon_Gen_Charm_POWHEG");
    h_M_Muon_Gen_Charm_POWHEG->SetTitle("Distr.");

    TF1 *pdf_Pt_POWHEG = new TF1("pdf_Pt_POWHEG", FuncPtMass, 0, 40, 4);
    pdf_Pt_POWHEG->SetTitle("Fit");
    pdf_Pt_POWHEG->SetParameter(3, 1);
    pdf_Pt_POWHEG->SetParameter(0, 3.6);
    pdf_Pt_POWHEG->SetParameter(1, 2.81);
    pdf_Pt_POWHEG->SetParameter(2, 2.5);
    pdf_Pt_POWHEG->SetNpx(300);
    pdf_Pt_POWHEG->SetLineColor(kCyan + 2);
    pdf_Pt_POWHEG->SetLineWidth(3);
    h_Pt_Muon_Gen_Charm_POWHEG->Fit(pdf_Pt_POWHEG, "LR0I");
    h_Pt_Muon_Gen_Charm_POWHEG->GetXaxis()->SetRangeUser(0, 10);
    h_Pt_Muon_Gen_Charm_POWHEG->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    TH1F *ratio_Pt_fit_POWHEG = (TH1F *)h_Pt_Muon_Gen_Charm_POWHEG->Clone("ratio_Pt_fit_POWHEG");
    ratio_Pt_fit_POWHEG->Divide(pdf_Pt_POWHEG);
    ratio_Pt_fit_POWHEG->GetYaxis()->SetRangeUser(-0.2, 3.2);
    ratio_Pt_fit_POWHEG->GetYaxis()->SetTitle("point/fit");
    TCanvas *pt_POWHEG = histo_fit_ratio(h_Pt_Muon_Gen_Charm_POWHEG, pdf_Pt_POWHEG, ratio_Pt_fit_POWHEG, "canvas_DiMu_Pt_POWHEG", Header, kTRUE);

    TF1 *pdf_M_POWHEG = new TF1("pdf_M_POWHEG", FuncPtMass, 0, 40, 4);
    pdf_M_POWHEG->SetTitle("Fit");

    pdf_M_POWHEG->SetParameter(3, 1);
    pdf_M_POWHEG->SetParameter(0, 3.6);
    pdf_M_POWHEG->SetParameter(1, 2.81);
    pdf_M_POWHEG->SetParameter(2, 2.5);
    pdf_M_POWHEG->SetNpx(300);
    pdf_M_POWHEG->SetLineWidth(3);
    h_M_Muon_Gen_Charm_POWHEG->Fit(pdf_M_POWHEG, "LR0I");
    h_M_Muon_Gen_Charm_POWHEG->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
    TH1F *ratio_M_fit_POWHEG = (TH1F *)h_M_Muon_Gen_Charm_POWHEG->Clone("ratio_M_fit_POWHEG");
    ratio_M_fit_POWHEG->Divide(pdf_M_POWHEG);
    ratio_M_fit_POWHEG->GetYaxis()->SetRangeUser(-0.2, 3.2);
    ratio_M_fit_POWHEG->GetYaxis()->SetTitle("point/fit");

    TCanvas *M_POWHEG = histo_fit_ratio(h_M_Muon_Gen_Charm_POWHEG, pdf_M_POWHEG, ratio_M_fit_POWHEG, "canvas_DiMu_M_POWHEG", Header, kTRUE);

    TFile *fIn_Pythia = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Pythia_Sim/MB/Version3/save_mc_output/MB_MC_output_Hist_merged.root", "READ");

    TH2F *h_PtM_Muon_Gen_Charm_PYTHIA_MONASH = (TH2F *)fIn_Pythia->Get(Form("DiMuon_Gen/h_PtM_DiMuon_Gen_%s", Origin.Data()));
    h_PtM_Muon_Gen_Charm_PYTHIA_MONASH->SetName("h_PtM_Muon_Gen_MONASH");

    TH1F *h_Pt_Muon_Gen_Charm_PYTHIA_MONASH = (TH1F *)h_PtM_Muon_Gen_Charm_PYTHIA_MONASH->ProjectionX();
    h_Pt_Muon_Gen_Charm_PYTHIA_MONASH->SetName("h_Pt_Muon_Gen_Charm_PYTHIA_MONASH");
    h_Pt_Muon_Gen_Charm_PYTHIA_MONASH->SetTitle("Distr.");
    hist1D_graphic_opt(h_Pt_Muon_Gen_Charm_PYTHIA_MONASH, 5, 20, kGreen + 2, 1.0 / h_Pt_Muon_Gen_Charm_PYTHIA_MONASH->GetEntries());
    TH1F *h_M_Muon_Gen_Charm_PYTHIA_MONASH = (TH1F *)h_PtM_Muon_Gen_Charm_PYTHIA_MONASH->ProjectionY();
    h_M_Muon_Gen_Charm_PYTHIA_MONASH->SetName("h_M_Muon_Gen_Charm_PYTHIA_MONASH");
    h_M_Muon_Gen_Charm_PYTHIA_MONASH->SetTitle("Distr.");
    hist1D_graphic_opt(h_M_Muon_Gen_Charm_PYTHIA_MONASH, 5, 20, kGreen + 2, 1.0 / h_Pt_Muon_Gen_Charm_PYTHIA_MONASH->GetEntries());

    TF1 *pdf_Pt_MONASH = new TF1("pdf_Pt_MONASH", FuncPtMass, 0, 40, 4);
    pdf_Pt_MONASH->SetTitle("Fit");
    pdf_Pt_MONASH->SetParameter(3, 1);
    pdf_Pt_MONASH->SetParameter(0, 3.6);
    pdf_Pt_MONASH->SetParameter(1, 2.81);
    pdf_Pt_MONASH->SetParameter(2, 2.5);
    pdf_Pt_MONASH->SetNpx(300);
    pdf_Pt_MONASH->SetLineWidth(3);
    h_Pt_Muon_Gen_Charm_PYTHIA_MONASH->Fit(pdf_Pt_MONASH, "LR0I");
    h_Pt_Muon_Gen_Charm_PYTHIA_MONASH->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    TH1F *ratio_Pt_fit_MONASH = (TH1F *)h_Pt_Muon_Gen_Charm_PYTHIA_MONASH->Clone("ratio_Pt_fit_MONASH");
    ratio_Pt_fit_MONASH->Divide(pdf_Pt_MONASH);
    ratio_Pt_fit_MONASH->GetYaxis()->SetRangeUser(-0.2, 3.2);
    ratio_Pt_fit_MONASH->GetYaxis()->SetTitle("point/fit");

    TCanvas *Pt_MONASH = histo_fit_ratio(h_Pt_Muon_Gen_Charm_PYTHIA_MONASH, pdf_Pt_MONASH, ratio_Pt_fit_MONASH, "canvas_DiMu_Pt_Monash", Header, kTRUE);

    TF1 *pdf_M_MONASH = new TF1("pdf_M_MONASH", FuncPtMass, 0, 40, 4);
    pdf_M_MONASH->SetTitle("Fit");
    pdf_M_MONASH->SetParameter(3, 1);
    pdf_M_MONASH->SetParameter(0, 3.6);
    pdf_M_MONASH->SetParameter(1, 2.81);
    pdf_M_MONASH->SetParameter(2, 2.5);
    pdf_M_MONASH->SetNpx(300);
    pdf_M_MONASH->SetLineWidth(3);
    h_M_Muon_Gen_Charm_PYTHIA_MONASH->Fit(pdf_M_MONASH, "LR0I");
    h_M_Muon_Gen_Charm_PYTHIA_MONASH->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/it{c}^{2})^{-1}");

    TH1F *ratio_M_fit_MONASH = (TH1F *)h_M_Muon_Gen_Charm_PYTHIA_MONASH->Clone("ratio_M_fit_MONASH");
    ratio_M_fit_MONASH->Divide(pdf_M_MONASH);
    ratio_M_fit_MONASH->GetYaxis()->SetRangeUser(-0.2, 3.2);
    ratio_M_fit_MONASH->GetYaxis()->SetTitle("point/fit");

    TCanvas *M_MONASH = histo_fit_ratio(h_M_Muon_Gen_Charm_PYTHIA_MONASH, pdf_M_MONASH, ratio_M_fit_MONASH, "canvas_DiMu_M_Monash", Header, kTRUE);

    TFile *fIn_Pythia_MNR = new TFile("test/HF_MC_output_Hist_294009.root", "READ");

    TH2F *h_PtM_Muon_Gen_Charm_PYTHIA_MNR = (TH2F *)fIn_Pythia_MNR->Get(Form("DiMuon_Gen/h_PtM_DiMuon_Gen_%s", Origin.Data()));
    h_PtM_Muon_Gen_Charm_PYTHIA_MNR->SetName("h_PtM_Muon_Gen_MNR");

    TH1F *h_Pt_Muon_Gen_Charm_PYTHIA_MNR = (TH1F *)h_PtM_Muon_Gen_Charm_PYTHIA_MNR->ProjectionX();
    h_Pt_Muon_Gen_Charm_PYTHIA_MNR->SetName("h_Pt_Muon_Gen_Charm_PYTHIA_MNR");
    h_Pt_Muon_Gen_Charm_PYTHIA_MNR->SetTitle("Distr.");
    hist1D_graphic_opt(h_Pt_Muon_Gen_Charm_PYTHIA_MNR, 5, 20, kMagenta + 2, 1.0 / h_Pt_Muon_Gen_Charm_PYTHIA_MNR->GetEntries());
    TH1F *h_M_Muon_Gen_Charm_PYTHIA_MNR = (TH1F *)h_PtM_Muon_Gen_Charm_PYTHIA_MNR->ProjectionY();
    h_M_Muon_Gen_Charm_PYTHIA_MNR->SetName("h_M_Muon_Gen_Charm_PYTHIA_MNR");
    h_M_Muon_Gen_Charm_PYTHIA_MNR->SetTitle("Distr.");
    hist1D_graphic_opt(h_M_Muon_Gen_Charm_PYTHIA_MNR, 5, 20, kMagenta + 2, 1.0 / h_Pt_Muon_Gen_Charm_PYTHIA_MNR->GetEntries());

    TF1 *pdf_Pt_MNR = new TF1("pdf_Pt_MNR", FuncPtMass, 0, 40, 4);
    pdf_Pt_MNR->SetTitle("Fit");
    pdf_Pt_MNR->SetParameter(3, 1);
    pdf_Pt_MNR->SetParameter(0, 3.6);
    pdf_Pt_MNR->SetParameter(1, 2.81);
    pdf_Pt_MNR->SetParameter(2, 2.5);
    pdf_Pt_MNR->SetNpx(300);
    pdf_Pt_MNR->SetLineWidth(3);
    h_Pt_Muon_Gen_Charm_PYTHIA_MNR->Fit(pdf_Pt_MNR, "LR0I");
    h_Pt_Muon_Gen_Charm_PYTHIA_MNR->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    TH1F *ratio_Pt_fit_MNR = (TH1F *)h_Pt_Muon_Gen_Charm_PYTHIA_MNR->Clone("ratio_Pt_fit_MNR");
    ratio_Pt_fit_MNR->Divide(pdf_Pt_MNR);
    ratio_Pt_fit_MNR->GetYaxis()->SetRangeUser(-0.2, 3.2);
    ratio_Pt_fit_MNR->GetYaxis()->SetTitle("point/fit");

    TCanvas *Pt_MNR = histo_fit_ratio(h_Pt_Muon_Gen_Charm_PYTHIA_MNR, pdf_Pt_MNR, ratio_Pt_fit_MNR, "canvas_DiMu_Pt_MNR", Header, kTRUE);

    TF1 *pdf_M_MNR = new TF1("pdf_M_MNR", FuncPtMass, 0, 40, 4);
    pdf_M_MNR->SetTitle("Fit");
    pdf_M_MNR->SetParameter(3, 1);
    pdf_M_MNR->SetParameter(0, 3.6);
    pdf_M_MNR->SetParameter(1, 2.81);
    pdf_M_MNR->SetParameter(2, 2.5);
    pdf_M_MNR->SetNpx(300);
    pdf_M_MNR->SetLineWidth(3);
    h_M_Muon_Gen_Charm_PYTHIA_MNR->Fit(pdf_M_MNR, "LR0I");
    h_M_Muon_Gen_Charm_PYTHIA_MNR->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

    TH1F *ratio_M_fit_MNR = (TH1F *)h_M_Muon_Gen_Charm_PYTHIA_MNR->Clone("ratio_M_fit_MNR");
    ratio_M_fit_MNR->Divide(pdf_M_MNR);
    ratio_M_fit_MNR->GetYaxis()->SetRangeUser(-0.2, 3.2);
    ratio_M_fit_MNR->GetYaxis()->SetTitle("point/fit");

    TCanvas *M_MNR = histo_fit_ratio(h_M_Muon_Gen_Charm_PYTHIA_MNR, pdf_M_MNR, ratio_M_fit_MNR, "canvas_DiMu_M_MNR", Header, kTRUE);

    TH1F *hist_M_Monash_PDF = new TH1F("hist_M_Monash_PDF", "PYTHIA8 Monash", 400, 0, 40);
    hist_M_Monash_PDF->Add(pdf_M_MONASH);
    hist_M_Monash_PDF->GetYaxis()->SetTitle("PDF");

    hist1D_graphic_opt(hist_M_Monash_PDF, 5, 20, kGreen + 2, 1.0);
    TH1F *hist_M_POWHEG_PDF = new TH1F("hist_M_POWHEG_PDF", "Powheg+PYTHIA6", 400, 0, 40);
    hist_M_POWHEG_PDF->Add(pdf_M_POWHEG);
    hist_M_POWHEG_PDF->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");

    hist1D_graphic_opt(hist_M_POWHEG_PDF, 5, 20, kCyan + 1, 1.0);
    TH1F *hist_M_MNR_PDF = new TH1F("hist_M_MNR_PDF", "PYTHIA8 MNR", 400, 0, 40);
    hist_M_MNR_PDF->Add(pdf_M_MNR);
    hist1D_graphic_opt(hist_M_MNR_PDF, 5, 20, kMagenta + 2, 1.0);

    TH1F *ratio_M_POWHEG_MONASH = (TH1F *)hist_M_POWHEG_PDF->Clone("ratio_M_POWHEG_MONASH");
    ratio_M_POWHEG_MONASH->Divide(hist_M_Monash_PDF);
    ratio_M_POWHEG_MONASH->GetYaxis()->SetRangeUser(0.6, 2.6);
    ratio_M_POWHEG_MONASH->GetYaxis()->SetTitle("Ratio");
    ratio_M_POWHEG_MONASH->SetTitle("MNR/Monash");

    TH1F *ratio_M_MNR_MONASH = (TH1F *)hist_M_MNR_PDF->Clone("ratio_M_MNR_MONASH");
    ratio_M_MNR_MONASH->Divide(hist_M_Monash_PDF);
    ratio_M_MNR_MONASH->GetYaxis()->SetTitle("Ratio");
    ratio_M_MNR_MONASH->SetTitle("Powheg+PYTHIA6/Monash");

    TCanvas *M_shape_comp = three_histo_ratio(hist_M_Monash_PDF, hist_M_POWHEG_PDF, hist_M_MNR_PDF, ratio_M_POWHEG_MONASH, ratio_M_MNR_MONASH, "canvas_M_shape_comp", Header, kTRUE);

    TH1F *hist_Pt_Monash_PDF = new TH1F("hist_Pt_Monash_PDF", "PYTHIA8 Monash", 400, 0, 40);
    hist_Pt_Monash_PDF->Add(pdf_Pt_MONASH);
    hist_Pt_Monash_PDF->GetYaxis()->SetTitle("PDF");
    hist1D_graphic_opt(hist_Pt_Monash_PDF, 5, 20, kGreen + 2, 1.0);
    TH1F *hist_Pt_POWHEG_PDF = new TH1F("hist_Pt_POWHEG_PDF", "Powheg+PYTHIA6", 400, 0, 40);
    hist_Pt_POWHEG_PDF->Add(pdf_Pt_POWHEG);
    hist_Pt_POWHEG_PDF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hist1D_graphic_opt(hist_Pt_POWHEG_PDF, 5, 20, kCyan + 1, 1.0);
    TH1F *hist_Pt_MNR_PDF = new TH1F("hist_Pt_MNR_PDF", "PYTHIA8 MNR", 400, 0, 40);
    hist_Pt_MNR_PDF->Add(pdf_Pt_MNR);

    hist1D_graphic_opt(hist_Pt_MNR_PDF, 5, 20, kMagenta + 2, 1.0);

    TH1F *ratio_Pt_MNR_MONASH = (TH1F *)hist_Pt_MNR_PDF->Clone("ratio_Pt_MNR_MONASH");
    ratio_Pt_MNR_MONASH->Divide(hist_Pt_Monash_PDF);
    ratio_Pt_MNR_MONASH->GetYaxis()->SetRangeUser(0.6, 2.6);
    ratio_Pt_MNR_MONASH->GetYaxis()->SetTitle("Ratio");
    ratio_Pt_MNR_MONASH->SetTitle("MNR/Monash");

    TH1F *ratio_Pt_POWHEG_MONASH = (TH1F *)hist_Pt_POWHEG_PDF->Clone("ratio_Pt_POWHEG_MONASH");
    ratio_Pt_POWHEG_MONASH->Divide(hist_Pt_Monash_PDF);
    ratio_Pt_POWHEG_MONASH->GetYaxis()->SetTitle("Ratio");
    ratio_Pt_POWHEG_MONASH->SetTitle("Powheg+PYTHIA6/Monash");

    TCanvas *Pt_shape_comp = three_histo_ratio(hist_Pt_Monash_PDF, hist_Pt_POWHEG_PDF, hist_Pt_MNR_PDF, ratio_Pt_POWHEG_MONASH, ratio_Pt_MNR_MONASH, "canvas_Pt_shape_comp", Header, kTRUE);
}

double FuncPtMass(double *x, double *par)
{
    return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
}
