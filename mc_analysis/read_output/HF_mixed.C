#if !defined(__CINT__) || defined(__CLING__)
#include "TROOT.h"
#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TGrid.h"
#endif

void HF_mixed()
{

    TFile *file_in = new TFile("/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Mixed_Analysis/12_10_22/HF/HF_HFmixedAnalysis_merged.root", "READ");

    const Int_t n_Fraction = 11;
    Double_t Fraction[n_Fraction] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    Double_t Charm_Fraction_ccut[n_Fraction] = {0};
    Double_t Beauty_Fraction_bcut[n_Fraction] = {0};

    Double_t Mixed_Fraction_ccut[n_Fraction] = {0};
    Double_t Mixed_Fraction_bcut[n_Fraction] = {0};

    TH1D *h_ULSMixed_ccut;
    TH1D *h_ULSMixed_bcut;
    TH1D *h_ULSCharm;
    TH1D *h_ULSBeauty;

    TH1D *h_ULSCharm_total = (TH1D *)file_in->Get(Form("Dimuon/ULS/Fraction1.00/h_MDiMu_ccut_ULS_fromCharm_fraction1.00"));
    TH1D *h_ULSBeauty_total = (TH1D *)file_in->Get(Form("Dimuon/ULS/Fraction1.00/h_MDiMu_bcut_ULS_fromBeauty_fraction1.00"));

    for (Int_t i = 0; i < n_Fraction; i++)
    {
        h_ULSCharm = (TH1D *)file_in->Get(Form("Dimuon/ULS/Fraction%0.2f/h_MDiMu_ccut_ULS_fromCharm_fraction%0.2f", Fraction[i], Fraction[i]));
        h_ULSBeauty = (TH1D *)file_in->Get(Form("Dimuon/ULS/Fraction%0.2f/h_MDiMu_bcut_ULS_fromBeauty_fraction%0.2f", Fraction[i], Fraction[i]));

        h_ULSMixed_ccut = (TH1D *)file_in->Get(Form("Dimuon/ULS/Fraction%0.2f/h_MDiMu_ccut_ULS_fromMixed_fraction%0.2f", Fraction[i], Fraction[i]));

        h_ULSMixed_bcut = (TH1D *)file_in->Get(Form("Dimuon/ULS/Fraction%0.2f/h_MDiMu_bcut_ULS_fromMixed_fraction%0.2f", Fraction[i], Fraction[i]));

        printf("%0.0f  %0.0f  %0.0f \n", h_ULSCharm->GetEntries(), h_ULSBeauty->GetEntries(), h_ULSMixed_ccut->GetEntries());

        printf("%0.0f  %0.0f  %0.0f \n", h_ULSCharm->GetEntries(), h_ULSBeauty->GetEntries(), h_ULSMixed_bcut->GetEntries());

        Double_t total_ccut = h_ULSCharm->GetEntries() + h_ULSBeauty_total->GetEntries() + h_ULSMixed_ccut->GetEntries();

        Double_t total_bcut = h_ULSCharm_total->GetEntries() + h_ULSBeauty->GetEntries() + h_ULSMixed_bcut->GetEntries();

        Mixed_Fraction_ccut[i] = (h_ULSMixed_ccut->GetEntries() / total_ccut) * 100;
        Mixed_Fraction_bcut[i] = (h_ULSMixed_bcut->GetEntries() / total_bcut) * 100;

        Charm_Fraction_ccut[i] = (h_ULSCharm->GetEntries() / total_ccut) * 100;
        Beauty_Fraction_bcut[i] = (h_ULSBeauty->GetEntries() / total_bcut) * 100;

        Fraction[i] = Fraction[i] * 100;
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->cd();
    c1->SetTicks();
    c1->cd();
    gPad->SetTicks();
    // gPad->SetLogy(1);
    gPad->SetTopMargin(0.10);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.15);

    TH2D *h_grid_dimumixed = new TH2D("h_grid_dimumixed", "Mixed HF Dimuon from MC ", 5, 49.5, 99.5, 10, 0, 10);
    h_grid_dimumixed->SetStats(0);
    for (Int_t q = 0; q < n_Fraction; q++)
    {
        // printf("%s\n", Form("Frac %0.2f", Fractions[q]));
        h_grid_dimumixed->GetXaxis()->SetBinLabel(q + 1, Form("%0.0f", Fraction[q]));

        // Fraction[q]=Fraction[q]*100;
    }

    // h_grid_dimumixed->LabelsDeflate();
    // h_grid_dimumixed->GetXaxis()->SetTitle("% Fraction");
    // h_grid_dimumixed->GetYaxis()->SetTitle("% Mixed Dimuon over total");
    // h_grid_dimumixed->GetXaxis()->SetTitleOffset(1.1);
    // h_grid_dimumixed->GetXaxis()->SetTitleSize(0.045);
    // h_grid_dimumixed->GetXaxis()->SetLabelSize(0.045);

    // h_grid_dimumixed->GetYaxis()->SetNdivisions(505);
    // h_grid_dimumixed->GetYaxis()->SetTitleOffset(1.1);
    // h_grid_dimumixed->GetYaxis()->SetTitleSize(0.0475);
    // h_grid_dimumixed->GetYaxis()->SetLabelSize(0.045);
    // h_grid_dimumixed->Draw("");

    Double_t Fractions[n_Fraction] = {50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
    TGraph *mixed_graph_ccut = new TGraph(n_Fraction, Fractions, Mixed_Fraction_ccut);
    // mixed_graph_ccut->LabelsDeflate();
    mixed_graph_ccut->GetXaxis()->SetTitle("% Fraction");
    mixed_graph_ccut->GetYaxis()->SetTitle("% Mixed Dimuon over total");
    mixed_graph_ccut->GetXaxis()->SetTitleOffset(1.1);
    mixed_graph_ccut->GetXaxis()->SetTitleSize(0.045);
    mixed_graph_ccut->GetXaxis()->SetLabelSize(0.045);

    mixed_graph_ccut->GetYaxis()->SetNdivisions(505);
    mixed_graph_ccut->GetYaxis()->SetTitleOffset(1.1);
    mixed_graph_ccut->GetYaxis()->SetTitleSize(0.0475);
    mixed_graph_ccut->GetYaxis()->SetLabelSize(0.045);

    mixed_graph_ccut->SetName("mixed_graph_ccut");
    mixed_graph_ccut->SetTitle(" ");
    mixed_graph_ccut->SetMarkerSize(1.5);
    mixed_graph_ccut->SetMarkerStyle(20);
    mixed_graph_ccut->SetMarkerColor(kViolet);
    mixed_graph_ccut->SetLineColor(kViolet);
    // mixed_graph_ccut->GetXaxis()->SetTitle("Fraction c");
    // mixed_graph_ccut->GetYaxis()->SetTitle("Fraction Mixed cutting on charm");
    mixed_graph_ccut->SetMinimum(0);
    mixed_graph_ccut->SetMaximum(8);
    mixed_graph_ccut->Draw();
    mixed_graph_ccut->Draw("PESAME");

    TGraph *mixed_graph_bcut = new TGraph(n_Fraction, Fractions, Mixed_Fraction_bcut);
    mixed_graph_bcut->SetName("mixed_graph_bcut");
    mixed_graph_bcut->SetTitle(" ");
    mixed_graph_bcut->SetMarkerSize(1.5);
    mixed_graph_bcut->SetMarkerStyle(20);
    mixed_graph_bcut->SetMarkerColor(kGreen);
    mixed_graph_bcut->SetLineColor(kGreen);
    mixed_graph_bcut->GetXaxis()->SetTitle("Fraction b");
    mixed_graph_bcut->GetYaxis()->SetTitle("Fraction Mixed cutting on beauty");
    mixed_graph_bcut->Draw("PLSAME");

    TLegend *legend = new TLegend(0.175, 0.575, 0.65, 0.775);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);
    legend->SetTextAlign(12);
    legend->AddEntry(mixed_graph_ccut, "Cut on charm", "LP");
    legend->AddEntry(mixed_graph_bcut, "Cut on beauty", "LP");
    legend->Draw();

    c1->SaveAs("MC_mixed.pdf");

    new TCanvas();

    TGraph *charm_graph_ccut = new TGraph(n_Fraction, Fraction, Charm_Fraction_ccut);
    charm_graph_ccut->SetName("charm_graph_ccut");
    charm_graph_ccut->SetTitle(" ");
    charm_graph_ccut->SetMarkerSize(1.5);
    charm_graph_ccut->SetMarkerStyle(20);
    charm_graph_ccut->SetMarkerColor(kViolet);
    charm_graph_ccut->SetLineColor(kViolet);
    charm_graph_ccut->GetXaxis()->SetTitle("Fraction c");
    charm_graph_ccut->GetYaxis()->SetTitle("Fraction Charm cutting on charm");
    charm_graph_ccut->Draw("");

    new TCanvas();

    TGraph *beauty_graph_bcut = new TGraph(n_Fraction, Fraction, Beauty_Fraction_bcut);
    beauty_graph_bcut->SetName("beauty_graph_bcut");
    beauty_graph_bcut->SetTitle(" ");
    beauty_graph_bcut->SetMarkerSize(1.5);
    beauty_graph_bcut->SetMarkerStyle(20);
    beauty_graph_bcut->SetMarkerColor(kGreen);
    beauty_graph_bcut->SetLineColor(kGreen);
    beauty_graph_bcut->GetXaxis()->SetTitle("Fraction bs");
    beauty_graph_bcut->GetYaxis()->SetTitle("Fraction Beauty cutting on beauty");
    beauty_graph_bcut->Draw("");
}