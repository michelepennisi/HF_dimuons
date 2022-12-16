/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN

#include "RooRealVar.h"
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
using namespace RooFit;

Double_t Low_Pt = 0.0;
Double_t High_Pt = 30.0;

Double_t Low_Mass = 4.0;
Double_t High_Mass = 20.0;

Int_t n_charm_input = 51000;
Int_t n_beauty_input = 84000;

// Tot dimuoni 77114
// Da beauty 77114*0.605=46653,97
// Da charm 77114*0.35=26989,9
// 3084,56
//-------------------------------------------------------------
//  Fit MC to extract pT and mass shapes
//-------------------------------------------------------------
//---------------------------------------------------------------------------------------//
TCanvas *printRooPlot_residual(RooPlot *frame, RooPlot *frame2, TString roohist_name[5]);
TCanvas *printRooPlot_ratio(RooPlot *frame, Bool_t norm, RooRealVar *fit_output[3], Int_t choice, TString roohist_name[5], TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x);

void single_unbinned_fit_data_sample(Int_t choice = 0, Double_t Low_Mass = 4.0, Double_t High_Mass = 30.0, Double_t Low_Pt = 0.0, Double_t High_Pt = 30.0, Int_t Binning_m = 26, Int_t Binning_pt = 30)
{
  gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
  gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");

  TFile *fIn = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/root_files/pdfMC_unbinned.root");

  RooWorkspace *w = (RooWorkspace *)fIn->Get("w");
  w->Print();

  RooRealVar *m = w->var("m");
  RooRealVar *pt = w->var("pt");
  m->setBins(Binning_m);
  pt->setBins(Binning_pt);

  RooCategory sample("sample", "sample");
  sample.defineType("mass");
  sample.defineType("transversemomentum");
  TFile *fIn_data = new TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/TreeResults_merged.root", "READ");

  // Taking data saved in tree
  TTree *tree_data = (TTree *)fIn_data->Get("rec_data_tree");
  RooDataSet *unbinned_M_Dimu_data = new RooDataSet("M_Dimu_data", "M_Dimu_data", RooArgSet(*m), Import(*tree_data));
  RooDataSet *unbinned_Pt_Dimu_data = new RooDataSet("Pt_Dimu_data", "Pt_Dimu_data", RooArgSet(*pt), Import(*tree_data));

  // Fix parameters to the MCRooSimultaneous
  RooRealVar *B_DimuMassFromCharm = w->var("B_DimuMassFromCharm");
  B_DimuMassFromCharm->setConstant(kTRUE);
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

  TString pdfMass_name[3];
  TString pdfPt_name[3];
  pdfMass_name[0].Form("pdfDimuMassFromcharm");
  pdfMass_name[1].Form("pdfDimuMassFrombeauty");
  pdfMass_name[2].Form("pdfDimuMassFrommixed");

  pdfPt_name[0].Form("pdfDimuPtFromcharm");
  pdfPt_name[1].Form("pdfDimuPtFrombeauty");
  pdfPt_name[2].Form("pdfDimuPtFrommixed");

  RooAbsPdf *pdfDimuMassFromCharm = w->pdf("pdfDimuMassFromcharm");
  RooAbsPdf *pdfDimuPtFromCharm = w->pdf("pdfDimuPtFromcharm");

  RooAbsPdf *pdfDimuMassFromBeauty = w->pdf("pdfDimuMassFrombeauty");
  RooAbsPdf *pdfDimuPtFromBeauty = w->pdf("pdfDimuPtFrombeauty");

  RooAbsPdf *pdfDimuMassFromMixed = w->pdf("pdfDimuMassFrommixed");
  RooAbsPdf *pdfDimuPtFromMixed = w->pdf("pdfDimuPtFrommixed");

  RooPlot *frame = m->frame(Title("Imported TH1 with Poisson error bars"));
  pdfDimuMassFromCharm->plotOn(frame, LineColor(kMagenta + 2));
  pdfDimuMassFromBeauty->plotOn(frame, LineColor(kSpring - 6));
  pdfDimuMassFromMixed->plotOn(frame, LineColor(kAzure + 9));
  TCanvas *c2 = new TCanvas();
  c2->cd();
  frame->Draw();

  RooRealVar *Mass_normForC = new RooRealVar("Mass_n_charm_output", "number dimuon from c", 28440, 0, 200000);
  RooRealVar *Mass_normForB = new RooRealVar("Mass_n_beauty_output", "number dimuon from b", 48000, 0, 200000);
  RooRealVar *Mass_normForMixed = new RooRealVar("Mass_n_mixed_output", "number dimuon from b,c", 3160);
  Mass_normForMixed->setConstant(kTRUE);

  RooRealVar *Pt_normForC = new RooRealVar("Pt_n_charm_output", "number dimuon from c", 28440, 0, 200000);
  RooRealVar *Pt_normForB = new RooRealVar("Pt_n_beauty_output", "number dimuon from b", 48000, 0, 200000);
  RooRealVar *Pt_normForMixed = new RooRealVar("Pt_n_mixed_output", "number dimuon from b,c", 3160);
  Pt_normForMixed->setConstant(kTRUE);

  RooAddPdf *m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed), RooArgList(*Mass_normForC, *Mass_normForB, *Mass_normForMixed));

  RooAddPdf *pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed), RooArgList(*Pt_normForC, *Pt_normForB, *Pt_normForMixed));

  // m->setRange("mregion1", 4, 9);
  // m->setRange("mregion2", 11, 30);

  // Define the pdf for the simultaneous fit

  RooFitResult *Mass_r = m_model->fitTo(*unbinned_M_Dimu_data, Minimizer("Minuit2"), Save(), SumW2Error(true));

  RooFitResult *Pt_r = pt_model->fitTo(*unbinned_Pt_Dimu_data, Minimizer("Minuit2"), Save(), SumW2Error(true));

  RooRealVar *Mass_fit_output[3];
  Mass_fit_output[0] = (RooRealVar *)Mass_r->floatParsFinal().find("Mass_n_charm_output");
  Mass_fit_output[1] = (RooRealVar *)Mass_r->floatParsFinal().find("Mass_n_beauty_output");
  Mass_fit_output[2] = Mass_normForMixed;

  RooRealVar *Pt_fit_output[3];
  Pt_fit_output[0] = (RooRealVar *)Pt_r->floatParsFinal().find("Pt_n_charm_output");
  Pt_fit_output[1] = (RooRealVar *)Pt_r->floatParsFinal().find("Pt_n_beauty_output");
  Pt_fit_output[2] = Pt_normForMixed;

  RooPlot *m_frame = m->frame(Title("m_frame"));
  m_frame->SetTitle(" ");
  m_frame->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
  m_frame->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");

  unbinned_M_Dimu_data->plotOn(m_frame, Name("Mass_data_set"), DrawOption("PEZ"));
  m_model->plotOn(m_frame, Name("pdfmass"), LineStyle(kSolid), LineColor(kRed));
  m_model->plotOn(m_frame, Name("pdfmasscharm"), Components(Form("%s", pdfMass_name[0].Data())), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
  m_model->plotOn(m_frame, Name("pdfmassbeauty"), Components(Form("%s", pdfMass_name[1].Data())), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
  m_model->plotOn(m_frame, Name("pdfmassmixed"), Components(Form("%s", pdfMass_name[2].Data())), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

  RooPlot *pt_frame = pt->frame(Title("pt_frame"));
  pt_frame->SetTitle(" ");
  pt_frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  pt_frame->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

  unbinned_Pt_Dimu_data->plotOn(pt_frame, Name("Pt_data_set"), DrawOption("PEZ"));
  pt_model->plotOn(pt_frame, Name("pdfpt"), LineStyle(kSolid), LineColor(kRed));
  pt_model->plotOn(pt_frame, Name("pdfptcharm"), Components(Form("%s", pdfPt_name[0].Data())), LineStyle(kDashed), LineColor(kMagenta - 2), LineWidth(5));
  pt_model->plotOn(pt_frame, Name("pdfptbeauty"), Components(Form("%s", pdfPt_name[1].Data())), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
  pt_model->plotOn(pt_frame, Name("pdfptmixed"), Components(Form("%s", pdfPt_name[2].Data())), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

  TString Pt_name[5] = {"Mass_data_set", "pdfpt", "pdfptcharm", "pdfptbeauty", "pdfptmixed"};

  TString Mass_name[5] = {"Pt_data_set", "pdfmass", "pdfmasscharm", "pdfmassbeauty", "pdfmassmixed"};

  // Save the fit as TF1

  RooArgSet *m_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*m_model).snapshot(true)); // True means copy the PDF and everything it depends on
  auto &m_modelcopiedPdf = static_cast<RooAbsPdf &>((*m_modelcopyOfEverything)["m_model"]);          // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
  RooArgSet *m_modelobs = m_modelcopiedPdf.getObservables(*unbinned_M_Dimu_data);
  RooArgSet *m_modelPars = m_modelcopiedPdf.getParameters(*m_modelobs);
  TF1 *m_modelFunc = m_modelcopiedPdf.asTF(*m_modelobs, *m_modelPars, *m);

  RooArgSet *pt_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*pt_model).snapshot(true)); // True means copy the PDF and everything it depends on
  auto &pt_modelcopiedPdf = static_cast<RooAbsPdf &>((*pt_modelcopyOfEverything)["pt_model"]);         // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
  RooArgSet *pt_modelobs = pt_modelcopiedPdf.getObservables(*unbinned_Pt_Dimu_data);
  RooArgSet *pt_modelPars = pt_modelcopiedPdf.getParameters(*pt_modelobs);
  TF1 *pt_modelFunc = pt_modelcopiedPdf.asTF(*pt_modelobs, *pt_modelPars, *pt);

  // Save the roodataset as histogram

  TH1 *hDimuPt_data = unbinned_Pt_Dimu_data->createHistogram("h_ptdata", *pt, Binning(Binning_pt, 0, High_Pt));
  hDimuPt_data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hDimuPt_data->GetYaxis()->SetTitle("d#it{N}/#it{p}_{T} (GeV/#it{c})^{-1}");

  TH1 *hDimuM_data = unbinned_M_Dimu_data->createHistogram("h_mdata", *m, Binning(Binning_m, Low_Mass, High_Mass));
  hDimuM_data->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
  hDimuM_data->GetYaxis()->SetTitle("d#it{N}/#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

  // Print the fit with ratio

  TString info;
  info.Form("unbinned_data_optfit%d_M%0.0f_%0.0f_BinM%d_BinPt%d", choice, Low_Mass, High_Mass, Binning_m, Binning_pt);
  pt_frame->SetMaximum(1.2e+8);
  pt_frame->SetMinimum(1.5e-1);
  TCanvas *pt_canvas = printRooPlot_ratio(pt_frame, false, Pt_fit_output, choice, Pt_name, pt_modelFunc, hDimuPt_data, 0, High_Pt);
  pt_canvas->SetName(Form("pt_canvas_%s", info.Data()));
  pt_canvas->SetTitle(Form("pt_canvas_%s", info.Data()));
  pt_canvas->SaveAs(Form("plot/%s.pdf", pt_canvas->GetName()));

  m_frame->SetMaximum(1.2e+7);
  m_frame->SetMinimum(1.5e-1);
  TCanvas *m_canvas = printRooPlot_ratio(m_frame, false, Mass_fit_output, choice, Mass_name, m_modelFunc, hDimuM_data, Low_Mass, High_Mass);
  m_canvas->SetName(Form("m_canvas_%s", info.Data()));
  m_canvas->SetTitle(Form("m_canvas_%s", info.Data()));
  m_canvas->SaveAs(Form("plot/%s.pdf", m_canvas->GetName()));

  Mass_r->floatParsFinal().Print("s");
  Pt_r->floatParsFinal().Print("s");
}

TCanvas *printRooPlot_residual(RooPlot *frame, RooPlot *frame2, TString roohist_name[5])
{
  TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 800);
  canvas->SetTicks();
  canvas->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.375, 1.0, 1.0);
  pad1->SetTicks();
  pad1->SetLogy(1);
  pad1->SetTopMargin(0.05);
  pad1->SetRightMargin(0.03);
  pad1->SetLeftMargin(0.14);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.375);
  pad2->SetTicks();
  pad2->SetTopMargin(0.0);
  pad2->SetRightMargin(0.03);
  pad2->SetLeftMargin(0.14);
  pad2->SetBottomMargin(0.25);
  pad2->Draw();

  pad1->cd();
  frame->GetXaxis()->SetTitleOffset(0.85);
  frame->GetXaxis()->SetTitleSize(0.065);
  frame->GetXaxis()->SetLabelSize(0.055);

  frame->GetYaxis()->SetTitleOffset(1.);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleSize(0.0575);
  frame->GetYaxis()->SetLabelSize(0.055);

  frame->Draw();

  TLegend *legend = new TLegend(0.315, 0.15, 0.49, 0.4);
  // legend->SetNColumns(2);

  // legend->AddEntry((TObject*)0, "", "");

  legend->SetFillStyle(0);
  legend->SetLineColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.0525);
  legend->SetHeader("Data");
  legend->SetTextAlign(11);
  // legend->AddEntry(roohist_name[0], " ", "PE");

  TLegend *fit_legend = new TLegend(0.365, 0.15, 0.515, 0.4);
  fit_legend->SetTextAlign(11);
  fit_legend->SetFillStyle(0);
  fit_legend->SetBorderSize(0);
  fit_legend->SetTextSize(0.0525);
  fit_legend->SetHeader("Fit");
  fit_legend->SetTextAlign(11);
  // fit_legend->AddEntry(roohist_name[1], " ", "L");

  TLegend *pdf_legend = new TLegend(0.365, 0.15, 0.515, 0.4);
  pdf_legend->SetTextAlign(11);
  pdf_legend->SetFillStyle(0);
  pdf_legend->SetBorderSize(0);
  pdf_legend->SetTextSize(0.0525);
  pdf_legend->SetHeader("PDF");
  pdf_legend->SetTextAlign(11);
  // pdf_legend->AddEntry(roohist_name[2], " ", "L");
  // pdf_legend->AddEntry(roohist_name[3], " ", "L");
  // pdf_legend->AddEntry(roohist_name[4], " ", "L");

  legend->Draw();
  fit_legend->Draw();

  TLatex *letexTitle = new TLatex();
  letexTitle->SetTextSize(0.055);
  letexTitle->SetNDC();
  letexTitle->SetTextFont(42);

  // if (roohist_name[option].Contains("Charm")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow c,c");
  // if (roohist_name[option].Contains("Beauty")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,b");
  // if (roohist_name[option].Contains("Mixed")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,c");

  letexTitle->SetTextSize(0.065);
  letexTitle->DrawLatex(0.46, 0.86, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
  letexTitle->DrawLatex(0.405, 0.78, "PYTHIA8 Monash Tune, N_{ev} = 2 #upoint 10^{8}");
  // printf("WOW %s\n",roohist_name[option].Data() );
  if (roohist_name[0].BeginsWith("rooHistDimuPt"))
  {
    letexTitle->DrawLatex(0.405, 0.70, "Reconstructed #mu^{#plus}#mu^{#minus}, #it{m}_{#mu^{#plus}#mu^{#minus}} > 4 GeV/#it{c}^{2}");
    letexTitle->DrawLatex(0.405, 0.62, "2.5 < #it{y}_{#mu} < 4.0");
  }
  else
  {
    letexTitle->DrawLatex(0.405, 0.70, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu} < 4.0");
  }
  canvas->cd();
  pad2->cd();

  frame2->GetXaxis()->SetTitleOffset(1.15);
  frame2->GetXaxis()->SetTitleSize(0.10);
  frame2->GetXaxis()->SetLabelSize(0.075);

  frame2->GetYaxis()->SetTitleOffset(0.8);
  frame2->GetYaxis()->SetNdivisions(505);
  frame2->GetYaxis()->SetTitleSize(0.085);
  frame2->GetYaxis()->SetLabelSize(0.075);
  frame2->Draw();

  return canvas;
}

TCanvas *printRooPlot_ratio(RooPlot *frame, Bool_t norm, RooRealVar *fit_output[3], Int_t choice, TString roohist_name[5], TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x)
{
  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1100);
  canvas->SetTicks();
  canvas->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.375, 1.0, 1.0);
  pad1->SetTicks();
  pad1->SetLogy(1);
  pad1->SetTopMargin(0.05);
  pad1->SetRightMargin(0.03);
  pad1->SetLeftMargin(0.14);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.375);
  pad2->SetTicks();
  pad2->SetTopMargin(0.0);
  pad2->SetRightMargin(0.03);
  pad2->SetLeftMargin(0.14);
  pad2->SetBottomMargin(0.25);
  pad2->Draw();

  pad1->cd();
  // if (norm)
  // {
  //   frame->SetMaximum(1.5e-0);
  //   frame->SetMinimum(1.5e-12);
  // }
  // else
  // {
  //   frame->SetMaximum(1.5e+7);
  //   frame->SetMinimum(1.5e-2);
  // }

  frame->GetXaxis()->SetTitleOffset(1.3);
  frame->GetXaxis()->SetTitleSize(0.0475);
  frame->GetXaxis()->SetLabelSize(0.045);

  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleOffset(0.9);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetLabelSize(0.05);

  frame->Draw();

  TLegend *legend = new TLegend(0.675, 0.375, 1.0, 0.595);
  // legend->SetNColumns(2);
  // legend->AddEntry((TObject*)0, "", "");
  legend->SetFillStyle(0);
  legend->SetLineColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.0425);
  legend->SetHeader("Data");
  legend->SetTextAlign(11);
  legend->AddEntry("combDatamass", " ", "LP");

  TLegend *fit_legend = new TLegend(0.775, 0.375, 0.9, 0.595);
  fit_legend->SetTextAlign(11);
  fit_legend->SetFillStyle(0);
  fit_legend->SetBorderSize(0);
  fit_legend->SetTextSize(0.0425);
  fit_legend->SetHeader("Fit");
  if (roohist_name[1].Contains("pt"))
    fit_legend->AddEntry("pdfpt", " ", "L");
  else if (roohist_name[1].Contains("mass"))
    fit_legend->AddEntry("pdfmass", " ", "L");

  TLegend *pdf_legend = new TLegend(0.175, 0.03, 0.40, 0.305);
  pdf_legend->SetTextAlign(11);
  pdf_legend->SetFillStyle(0);
  pdf_legend->SetBorderSize(0);
  pdf_legend->SetTextSize(0.0425);
  pdf_legend->SetHeader("PDF");
  pdf_legend->AddEntry(roohist_name[2], " ", "L");
  pdf_legend->AddEntry(roohist_name[3], " ", "L");
  pdf_legend->AddEntry(roohist_name[4], " ", "L");

  legend->Draw();
  fit_legend->Draw();
  pdf_legend->Draw();

  TLatex *letexTitle = new TLatex();
  letexTitle->SetNDC();
  letexTitle->SetTextFont(42);

  letexTitle->DrawLatex(0.25, 0.19, "#mu^{#plus}#mu^{#minus} #leftarrow c,c");
  letexTitle->DrawLatex(0.25, 0.12, "#mu^{#plus}#mu^{#minus} #leftarrow b,b");
  letexTitle->DrawLatex(0.25, 0.05, "#mu^{#plus}#mu^{#minus} #leftarrow c,b");
  // if (roohist_name[option].Contains("Charm")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow c,c");
  // if (roohist_name[option].Contains("Beauty")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,b");
  // if (roohist_name[option].Contains("Mixed")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,c");

  letexTitle->SetTextSize(0.0475);
  // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
  letexTitle->DrawLatex(0.175, 0.875, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
  letexTitle->DrawLatex(0.175, 0.785, "LHC18p period");
  // printf("WOW %s\n",roohist_name[option].Data() );
  letexTitle->SetTextSize(0.0425);
  if (choice == 0)
  {
    letexTitle->DrawLatex(0.625, 0.825, Form("#it{N}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.1f #pm %0.1f", fit_output[0]->getVal(), fit_output[0]->getError()));
    letexTitle->DrawLatex(0.625, 0.725, Form("#it{N}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.1f #pm %0.1f", fit_output[1]->getVal(), fit_output[1]->getError()));
    letexTitle->DrawLatex(0.625, 0.625, Form("#it{N}^{fixed}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.0f", fit_output[2]->getVal()));
    // letexTitle->DrawLatex(0.675, 0.625, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f", fit_output[3]->getVal()));
  }
  else if (choice == 1)
  {
    letexTitle->DrawLatex(0.625, 0.825, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f %% #pm %0.3f %%", 100 * fit_output[0]->getVal(), 100 * fit_output[0]->getError()));
    letexTitle->DrawLatex(0.625, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f %% #pm %0.3f %%", 100 * fit_output[1]->getVal(), 100 * fit_output[1]->getError()));
    // letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.2e #pm %0.2e", fit_output[2]->getVal(), fit_output[2]->getError()));
    letexTitle->DrawLatex(0.625, 0.625, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f %%", 100 - 100 * fit_output[0]->getVal() - 100 * fit_output[1]->getVal()));
  }

  letexTitle->SetTextSize(0.0425);
  if (roohist_name[1].Contains("pt"))
  {
    letexTitle->DrawLatex(0.175, 0.695, "Reconstructed #mu^{#plus}#mu^{#minus}, #it{m}_{#mu^{#plus}#mu^{#minus}} > 4 GeV/#it{c}^{2}");
    letexTitle->DrawLatex(0.175, 0.605, "2.5 < #it{y}_{#mu} < 4.0");
  }
  else
  {
    letexTitle->DrawLatex(0.175, 0.695, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu} < 4.0");
  }

  pad2->cd();
  pad2->SetTicks();
  TLine *l = new TLine(minx, 1.0, 30.0, 1.0);
  l->SetLineWidth(3);
  l->SetLineStyle(2);
  l->SetLineColor(kRed);

  TLine *l1 = new TLine(minx, 0.5, 30.0, 0.5);
  l1->SetLineWidth(2);
  l1->SetLineStyle(9);
  l1->SetLineColor(kGray + 2);

  TLine *l2 = new TLine(minx, 1.5, 30.0, 1.5);
  l2->SetLineWidth(2);
  l2->SetLineStyle(9);
  l2->SetLineColor(kGray + 2);

  TH2D *h_grid_ratio = new TH2D("h_grid", "", 100, minx, max_x, 100, -0.5, 2.5);
  h_grid_ratio->SetTitle("");

  // TH1D *c_data = (TH1D *)data->Clone("c_data");

  h_grid_ratio->GetYaxis()->SetTitle(Form("#frac{Data}{Cocktail}"));
  h_grid_ratio->GetYaxis()->CenterTitle();
  h_grid_ratio->GetYaxis()->SetNdivisions(504);
  h_grid_ratio->GetYaxis()->SetTitleSize(0.08);
  // h_grid_ratio->GetYaxis()->SetTitleOffset(0.8);
  h_grid_ratio->GetYaxis()->SetLabelOffset(0.02);
  h_grid_ratio->GetYaxis()->SetLabelSize(0.1);

  h_grid_ratio->GetXaxis()->SetTitleSize(0.09);
  h_grid_ratio->GetXaxis()->SetTitleOffset(1.1);
  h_grid_ratio->GetXaxis()->SetLabelSize(0.08);
  h_grid_ratio->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());

  data->SetLineColor(kBlack);
  data->SetMarkerColor(kBlack);
  data->SetMarkerStyle(20);
  // data->Rebin(15);
  data->Scale(1. / data->Integral(), "width");
  data->Divide(pdf);

  h_grid_ratio->Draw();
  data->Draw("PESAME");
  l->Draw();
  l1->Draw();
  l2->Draw();

  return canvas;
}