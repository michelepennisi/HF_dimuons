#include <RooRealVar.h>
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TRatioPlot.h"
#include "RooHistPdf.h"
#include "TMath.h"
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "RooMinimizer.h"
#include "TH2.h"
#include "iostream"
#include "TROOT.h"
#include "TChain.h"
#include "TH3.h"
#include "TGraphMultiErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TKey.h"

using namespace RooFit;
using namespace std;

void progress_status(Int_t i_Event, Int_t total_entries)
{
    Double_t progress = (Double_t)i_Event / total_entries;
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " % (" << i_Event << "/ " << total_entries << ")\r" << std::scientific;
    std::cout.flush();

    std::cout << std::endl;
}

TCanvas *canvas_ratio(TString name)
{
    TCanvas *canvas = new TCanvas(name, name, 1100, 1200);
    canvas->Divide(1, 2, 0.0001, 0.0001);
    canvas->cd(1)->SetLogy();
    canvas->cd(1)->SetTopMargin(0.05);
    canvas->cd(1)->SetRightMargin(0.03);
    canvas->cd(1)->SetLeftMargin(0.17);
    canvas->cd(1)->SetBottomMargin(0.0);

    canvas->cd(2)->SetTopMargin(0.0);
    canvas->cd(2)->SetRightMargin(0.03);
    canvas->cd(2)->SetLeftMargin(0.17);
    canvas->cd(2)->SetBottomMargin(0.225);
    return canvas;
}

void Legend_settings(TLegend *legend, Double_t X1, Double_t X2, Double_t Y1, Double_t Y2, TString Header)
{

    legend->SetHeader(Header);
    legend->SetX1(X1);
    legend->SetX2(X2);
    legend->SetY1(Y1);
    legend->SetY2(Y2);
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0425);
};
void hist1D_graphic_opt(TH1F *hist, Bool_t Scale, Int_t Rebin, Style_t style, Color_t color, Double_t scale_factor)
{
    if (hist->IsZombie())
    {
        printf("Hist not found\n");
        return;
    }

    if (Rebin > 0)
        hist->Rebin(Rebin);

    if (Scale)
        hist->Scale(scale_factor, "width");
    else
        hist->Scale(scale_factor);

    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(style);
    hist->SetMarkerSize(1.5);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);

    hist->GetXaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetLabelSize(0.035);
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(2);
    hist->GetYaxis()->SetNdivisions(505);
}

TCanvas *canvas_noratio(TString name)
{
    TCanvas *canvas = new TCanvas(name, name, 1200, 1200);
    canvas->cd();
    canvas->GetPad(0)->SetTopMargin(0.05);
    canvas->GetPad(0)->SetRightMargin(0.05);
    canvas->GetPad(0)->SetLeftMargin(0.18);
    canvas->GetPad(0)->SetBottomMargin(0.115);

    // canvas->GetPad(0)->SetLogy();

    return canvas;
}

TCanvas *canvas_noratio_divide2(TString name)
{
    TCanvas *canvas = new TCanvas(name, name, 1600, 1100);
    canvas->Divide(2, 1);
    canvas->cd(1);
    canvas->GetPad(1)->SetTopMargin(0.05);
    canvas->GetPad(1)->SetRightMargin(0.05);
    canvas->GetPad(1)->SetLeftMargin(0.165);
    canvas->GetPad(1)->SetBottomMargin(0.13);
    canvas->GetPad(1)->SetLogy();

    canvas->cd(2);
    canvas->GetPad(2)->SetTopMargin(0.05);
    canvas->GetPad(2)->SetRightMargin(0.05);
    canvas->GetPad(2)->SetLeftMargin(0.165);
    canvas->GetPad(2)->SetBottomMargin(0.13);
    canvas->GetPad(2)->SetLogy();

    return canvas;
}

TCanvas *canvas_2d_print(TH2F *hist, TString name)
{
    TCanvas *canvas = new TCanvas(name, name, 1000, 1000);
    canvas->GetPad(0)->SetTopMargin(0.05);
    canvas->GetPad(0)->SetRightMargin(0.135);
    canvas->GetPad(0)->SetLeftMargin(0.135);
    canvas->GetPad(0)->SetBottomMargin(0.115);
    canvas->cd();
    gPad->SetGridx();
    gPad->SetGridy();
    hist->GetXaxis()->SetLabelSize(0.025);
    hist->GetXaxis()->SetTitleOffset(1.35);

    hist->GetYaxis()->SetLabelSize(0.025);
    hist->GetYaxis()->SetTitleOffset(1.5);

    hist->GetZaxis()->SetMaxDigits(2);
    hist->Draw("textCOLZSAME");

    return canvas;
}

TCanvas *histo_fit_ratio(TH1F *hist, TF1 *fit, TH1F *ratio_histo_fit, TString name, TString Header, Bool_t draw_lagend = kTRUE)
{
    TCanvas *canvas = new TCanvas(name, name, 1100, 1200);
    canvas->SetTicks();
    canvas->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.50, 1.0, 1.00);
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.50);

    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    pad1->SetTicks();
    pad1->SetLogy(1);
    pad1->SetTopMargin(0.05);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.14);
    pad1->SetBottomMargin(0.0);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->Draw("PE");
    fit->DrawCopy("SAME");
    TLegend *legend = new TLegend(0.575, 0.675, 0.825, 0.875);
    legend->SetHeader(Header);
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0425);

    legend->AddEntry(hist);
    legend->AddEntry(fit);
    if (draw_lagend)
        legend->Draw();

    pad2->cd();
    pad2->SetTicks();
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.03);
    pad2->SetLeftMargin(0.14);
    pad2->SetBottomMargin(0.175);

    ratio_histo_fit->GetXaxis()->SetTitleFont(42);
    ratio_histo_fit->GetXaxis()->SetLabelFont(42);
    ratio_histo_fit->GetXaxis()->SetLabelSize(0.05);
    ratio_histo_fit->GetXaxis()->SetTitleSize(0.055);
    ratio_histo_fit->GetXaxis()->SetTitleOffset(1.2);
    ratio_histo_fit->GetYaxis()->SetTitleOffset(0.95);
    ratio_histo_fit->GetYaxis()->SetLabelFont(42);
    ratio_histo_fit->GetYaxis()->SetNdivisions(505);
    ratio_histo_fit->GetYaxis()->SetLabelOffset(0.0175);
    ratio_histo_fit->GetYaxis()->SetLabelSize(0.05);
    ratio_histo_fit->GetYaxis()->SetTitleSize(0.055);
    ratio_histo_fit->GetYaxis()->SetTitleOffset(1.2);
    ratio_histo_fit->Draw("PESAME");

    return canvas;
}

TCanvas *two_histo_ratio(TH1F *hist1, TH1F *hist2, TH1F *ratio_histo_fit, TString name, TString Header, Bool_t draw_lagend = kTRUE, Bool_t setlogy = kTRUE)
{
    TCanvas *canvas = new TCanvas(name, name, 1100, 1200);
    canvas->SetTicks();
    canvas->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.35, 1.0, 1.00);
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.35);

    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    pad1->SetTicks();
    if (setlogy)
        pad1->SetLogy(1);
    pad1->SetTopMargin(0.05);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.17);
    pad1->SetBottomMargin(0.0);
    // if (hist1->GetMaximum() > hist2->GetMaximum())
    //     hist1->SetMaximum(hist1->GetMaximum() * 3.5);
    // else
    //     hist1->SetMaximum(hist2->GetMaximum() * 3.5);
    hist1->GetYaxis()->SetLabelSize(0.05);
    hist1->GetYaxis()->SetTitleSize(0.055);
    hist1->GetYaxis()->SetTitleOffset(1.4);
    hist1->GetXaxis()->SetNdivisions(510);

    TLegend *legend = new TLegend(0.575, 0.525, 0.825, 0.725);
    legend->SetHeader(Header);
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.045);

    auto copy_hist1 = hist1->DrawCopy("P");
    auto copy_hist2 = hist2->DrawCopy("PSAME");
    legend->AddEntry(copy_hist1);
    legend->AddEntry(copy_hist2);

    if (draw_lagend)
        legend->Draw();

    pad2->cd();
    pad2->SetTicks();
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.03);
    pad2->SetLeftMargin(0.17);
    pad2->SetBottomMargin(0.225);

    // ratio_histo_fit->GetYaxis()->SetRangeUser(ratio_histo_fit->GetMinimum()/2., ratio_histo_fit->GetMaximum()*2);
    // ratio_histo_fit->GetYaxis()->SetTitle(" ");
    ratio_histo_fit->GetXaxis()->SetTitleFont(42);
    ratio_histo_fit->GetXaxis()->SetLabelFont(42);
    ratio_histo_fit->GetXaxis()->SetLabelSize(0.095);
    ratio_histo_fit->GetXaxis()->SetTitleSize(0.095);
    ratio_histo_fit->GetXaxis()->SetTitleOffset(1.1);
    ratio_histo_fit->GetYaxis()->SetTitleOffset(0.95);
    ratio_histo_fit->GetYaxis()->SetLabelFont(42);
    ratio_histo_fit->GetYaxis()->SetNdivisions(505);
    ratio_histo_fit->GetYaxis()->SetLabelOffset(0.0175);
    ratio_histo_fit->GetYaxis()->SetLabelSize(0.095);
    ratio_histo_fit->GetYaxis()->SetTitleSize(0.095);
    ratio_histo_fit->GetYaxis()->SetTitleOffset(0.8);
    ratio_histo_fit->DrawCopy("P");

    return canvas;
}

TCanvas *three_histo_one_ratio(TH1F *hist1, TH1F *hist2, TH1F *hist3, TH1F *ratio_histo_fit, TString name, TString Header, Bool_t draw_lagend = kTRUE, Bool_t setlogy = kTRUE)
{
    TCanvas *canvas = new TCanvas(name, name, 900, 1200);
    canvas->SetTicks();
    canvas->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.50, 1.0, 1.00);
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.50);

    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    pad1->SetTicks();
    if (setlogy)
        pad1->SetLogy(1);
    pad1->SetTopMargin(0.05);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.17);
    pad1->SetBottomMargin(0.0);
    // if (hist1->GetMaximum() > hist2->GetMaximum())
    //     hist1->SetMaximum(hist1->GetMaximum() * 3.5);
    // else
    //     hist1->SetMaximum(hist2->GetMaximum() * 3.5);
    hist1->GetYaxis()->SetLabelSize(0.045);
    hist1->GetYaxis()->SetTitleSize(0.055);
    hist1->GetYaxis()->SetTitleOffset(1.4);
    hist1->Draw("PE");
    hist2->Draw("PESAME");
    hist3->Draw("PESAME");
    TLegend *legend = new TLegend(0.175, 0.675, 0.425, 0.95);
    legend->SetName("Legend");
    legend->SetHeader(Header);
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0425);

    legend->AddEntry(hist1);
    legend->AddEntry(hist2);
    legend->AddEntry(hist3);
    if (draw_lagend)
        legend->Draw();

    pad2->cd();
    pad2->SetTicks();
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.03);
    pad2->SetLeftMargin(0.17);
    pad2->SetBottomMargin(0.175);

    // ratio_histo_fit->GetYaxis()->SetRangeUser(ratio_histo_fit->GetMinimum()/2., ratio_histo_fit->GetMaximum()*2);
    // ratio_histo_fit->GetYaxis()->SetTitle(" ");
    ratio_histo_fit->GetXaxis()->SetTitleFont(42);
    ratio_histo_fit->GetXaxis()->SetLabelFont(42);
    ratio_histo_fit->GetXaxis()->SetLabelSize(0.05);
    ratio_histo_fit->GetXaxis()->SetTitleSize(0.06);
    ratio_histo_fit->GetXaxis()->SetTitleOffset(1.2);
    ratio_histo_fit->GetYaxis()->SetTitleOffset(1.8);
    ratio_histo_fit->GetYaxis()->SetLabelFont(42);
    ratio_histo_fit->GetYaxis()->SetNdivisions(505);
    ratio_histo_fit->GetYaxis()->SetLabelOffset(0.0175);
    ratio_histo_fit->GetYaxis()->SetLabelSize(0.05);
    ratio_histo_fit->Draw("PESAME");

    return canvas;
}

TCanvas *three_histo_ratio(TH1F *hist1, TH1F *hist2, TH1F *hist3, TH1F *ratio_histo12, TH1F *ratio_histo13, TString name, TString Header, Bool_t draw_lagend = kTRUE)
{
    TCanvas *canvas = new TCanvas(name, name, 1100, 1200);
    canvas->SetTicks();
    canvas->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.50, 1.0, 1.00);
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.50);

    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    pad1->SetTicks();
    pad1->SetLogy(1);
    pad1->SetTopMargin(0.05);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.14);
    pad1->SetBottomMargin(0.0);
    hist1->GetYaxis()->SetLabelSize(0.045);
    hist1->GetYaxis()->SetTitleSize(0.055);
    hist1->GetYaxis()->SetTitleOffset(1.1);
    hist1->Draw("PE");
    hist2->Draw("PESAME");
    hist3->Draw("PESAME");
    TLegend *legend = new TLegend(0.375, 0.575, 0.625, 0.875);
    legend->SetName("legend1");
    legend->SetHeader(Header);
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0425);

    legend->AddEntry(hist1);
    legend->AddEntry(hist2);
    legend->AddEntry(hist3);
    if (draw_lagend)
        legend->Draw();

    pad2->cd();
    pad2->SetTicks();
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.03);
    pad2->SetLeftMargin(0.14);
    pad2->SetBottomMargin(0.175);

    ratio_histo12->GetYaxis()->SetTitle(" ");
    ratio_histo12->GetXaxis()->SetTitleFont(42);
    ratio_histo12->GetXaxis()->SetLabelFont(42);
    ratio_histo12->GetXaxis()->SetLabelSize(0.05);
    ratio_histo12->GetXaxis()->SetTitleSize(0.06);
    ratio_histo12->GetXaxis()->SetTitleOffset(1.2);
    ratio_histo12->GetYaxis()->SetLabelFont(42);
    // ratio_histo12->GetYaxis()->SetRangeUser(-0.5, 2.5);
    ratio_histo12->GetYaxis()->SetNdivisions(505);
    ratio_histo12->GetYaxis()->SetLabelOffset(0.0175);
    ratio_histo12->GetYaxis()->SetLabelSize(0.05);
    ratio_histo12->Draw("PE");
    ratio_histo13->Draw("PESAME");

    TLegend *legend2 = new TLegend(0.375, 0.675, 0.625, 0.875);
    legend2->SetName("legend2");
    legend2->SetFillStyle(0);
    legend2->SetLineColor(kWhite);
    legend2->SetBorderSize(0);
    legend2->SetTextSize(0.0425);

    legend2->AddEntry(ratio_histo12);
    legend2->AddEntry(ratio_histo13);

    if (draw_lagend)
        legend2->Draw();

    return canvas;
}

TCanvas *four_histo(TH1F *hist1, TH1F *hist2, TH1F *hist3, TH1F *hist4, TString name, TString Header, Bool_t draw_lagend = kTRUE)
{
    TCanvas *canvas = new TCanvas(name, name, 1200, 1000);
    canvas->SetTicks();
    canvas->cd();

    gPad->SetTicks();
    gPad->SetLogy(1);
    canvas->GetPad(0)->SetTopMargin(0.05);
    canvas->GetPad(0)->SetRightMargin(0.05);
    canvas->GetPad(0)->SetLeftMargin(0.14);
    canvas->GetPad(0)->SetBottomMargin(0.14);

    hist1->GetXaxis()->SetTitleOffset(1.3);
    hist1->GetXaxis()->SetTitleSize(0.045);
    hist1->GetXaxis()->SetLabelSize(0.04);
    hist1->GetXaxis()->SetNdivisions(505);
    hist1->GetYaxis()->SetTitleOffset(1.3);
    hist1->GetYaxis()->SetTitleSize(0.045);
    hist1->GetYaxis()->SetLabelSize(0.04);
    hist1->GetYaxis()->SetNdivisions(505);
    hist1->Draw("PE");
    hist2->Draw("PESAME");
    hist3->Draw("PESAME");
    hist4->Draw("PESAME");
    TLegend *legend = new TLegend(0.5, 0.525, 0.7, 0.825);
    legend->SetHeader(Header);
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0375);

    legend->AddEntry(hist1);
    legend->AddEntry(hist2);
    legend->AddEntry(hist3);
    // legend->AddEntry(hist4);
    if (draw_lagend)
        legend->Draw();

    return canvas;
}