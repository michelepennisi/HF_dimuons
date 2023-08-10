#include "/home/michele_pennisi/cernbox/common_include.h"

// void hist1D_graphic_opt(TH1F *hist, Int_t Rebin, Style_t style, Color_t color, Double_t scale_factor);
// TCanvas *four_histo(TH1F *hist1, TH1F *hist2, TH1F *hist3,TH1F *hist4, TString name, TString Header, Bool_t draw_lagend = kTRUE);

void Z_study()
{
    TFile *fIn_MC = new TFile("test/powheg_Z_Lucas_MC_output_Hist_294925.root", "READ");
    fIn_MC->cd();

    TH1F *h_Nevents_Z_powheg = (TH1F *)fIn_MC->Get("h_Nevents");
    h_Nevents_Z_powheg->SetName("h_Nevents_Powheg");

    Double_t powheg_crosssection = 2e-06; //mb
    Double_t MB_crosssection = 56.42; // mb

    Double_t norm_factor=(powheg_crosssection/(MB_crosssection*h_Nevents_Z_powheg->GetBinContent(2)));

    std::cout<<norm_factor<<endl;

    TH2F *h_PtM_DiMuon_Rec_Z_nocut = (TH2F *)fIn_MC->Get("DiMuon_Rec/h_PtM_DiMuon_Rec_DY");
    TH1F *h_M_DiMuon_Rec_Z_nocut = (TH1F *)h_PtM_DiMuon_Rec_Z_nocut->ProjectionY();

    TH2F *h_PtM_DiMuon_Rec_Z_ptmucut09 = (TH2F *)fIn_MC->Get("DiMuon_Rec_Z/h_PtM_DiMuon_Rec_Z_ptmucut09");
    TH1F *h_M_DiMuon_Rec_Z_ptmucut09 = (TH1F *)h_PtM_DiMuon_Rec_Z_ptmucut09->ProjectionY();

    TH2F *h_PtM_DiMuon_Rec_Z_ptmucut10 = (TH2F *)fIn_MC->Get("DiMuon_Rec_Z/h_PtM_DiMuon_Rec_Z_ptmucut10");
    TH1F *h_M_DiMuon_Rec_Z_ptmucut10 = (TH1F *)h_PtM_DiMuon_Rec_Z_ptmucut10->ProjectionY();

    TH2F *h_PtM_DiMuon_Rec_Z_ptmucut20 = (TH2F *)fIn_MC->Get("DiMuon_Rec_Z/h_PtM_DiMuon_Rec_Z_ptmucut20");
    TH1F *h_M_DiMuon_Rec_Z_ptmucut20 = (TH1F *)h_PtM_DiMuon_Rec_Z_ptmucut20->ProjectionY();

    Int_t N_Rebin = 10;
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    h_M_DiMuon_Rec_Z_nocut->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    h_M_DiMuon_Rec_Z_nocut->GetYaxis()->SetRangeUser(2e-15,2);
    h_M_DiMuon_Rec_Z_nocut->SetTitle("#it{p}_{T,#mu} > 0 GeV/#it{c}");
    h_M_DiMuon_Rec_Z_ptmucut09->SetTitle("#it{p}_{T,#mu} > 0.9 GeV/#it{c}");
    h_M_DiMuon_Rec_Z_ptmucut10->SetTitle("#it{p}_{T,#mu} > 10 GeV/#it{c}");
    h_M_DiMuon_Rec_Z_ptmucut20->SetTitle("#it{p}_{T,#mu} > 20 GeV/#it{c}");


    hist1D_graphic_opt(h_M_DiMuon_Rec_Z_nocut, 10, 20, kGreen + 2, norm_factor);
    hist1D_graphic_opt(h_M_DiMuon_Rec_Z_ptmucut09, 10, 20, kMagenta + 2, norm_factor);
    hist1D_graphic_opt(h_M_DiMuon_Rec_Z_ptmucut10, 10, 20, kCyan + 2, norm_factor);
    hist1D_graphic_opt(h_M_DiMuon_Rec_Z_ptmucut20, 10, 20, kBlack, norm_factor);

    TCanvas *Z_rec_comp = four_histo(h_M_DiMuon_Rec_Z_nocut, h_M_DiMuon_Rec_Z_ptmucut09, h_M_DiMuon_Rec_Z_ptmucut10, h_M_DiMuon_Rec_Z_ptmucut20, "Z_rec_comp", " ", kTRUE);
}

// void hist1D_graphic_opt(TH1F *hist, Int_t Rebin, Style_t style, Color_t color, Double_t scale_factor)
// {
//     hist->Rebin(Rebin);
//     hist->Scale(scale_factor, "width");
//     hist->SetMarkerColor(color);
//     hist->SetMarkerStyle(style);
//     hist->SetMarkerSize(1.5);
//     hist->SetLineColor(color);
//     hist->SetLineWidth(2);

//     hist->GetXaxis()->SetLabelSize(0.04);
//     hist->GetYaxis()->SetLabelSize(0.04);
//     hist->GetXaxis()->SetTitleSize(0.045);
//     hist->GetYaxis()->SetTitleSize(0.045);
//     hist->GetXaxis()->SetTitleOffset(1.4);
//     hist->GetYaxis()->SetTitleOffset(1.6);
// }


// TCanvas *four_histo(TH1F *hist1, TH1F *hist2, TH1F *hist3,TH1F *hist4, TString name, TString Header, Bool_t draw_lagend = kTRUE)
// {
//     TCanvas *canvas = new TCanvas(name, name, 1000, 1600);
//     canvas->SetTicks();
//     canvas->cd();
    
//     gPad->SetTicks();
//     // gPad->SetLogy(1);
//     canvas->GetPad(0)->SetTopMargin(0.05);
//     canvas->GetPad(0)->SetRightMargin(0.05);
//     canvas->GetPad(0)->SetLeftMargin(0.165);
//     canvas->GetPad(0)->SetBottomMargin(0.13);
//     hist1->GetYaxis()->SetLabelSize(0.045);
//     hist1->GetYaxis()->SetTitleSize(0.055);
//     hist1->GetYaxis()->SetTitleOffset(1.1);
//     hist1->Draw("PE");
//     hist2->Draw("PESAME");
//     hist3->Draw("PESAME");
//     hist4->Draw("PESAME");
//     TLegend *legend = new TLegend(0.575, 0.575, 0.825, 0.875);
//     legend->SetHeader(Header);
//     legend->SetFillStyle(0);
//     legend->SetLineColor(kWhite);
//     legend->SetBorderSize(0);
//     legend->SetTextSize(0.0425);

//     legend->AddEntry(hist1);
//     legend->AddEntry(hist2);
//     legend->AddEntry(hist3);
//     legend->AddEntry(hist4);
//     if (draw_lagend)
//         legend->Draw();


//     return canvas;
// }