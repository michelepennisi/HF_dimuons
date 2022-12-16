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

//-------------------------------------------------------------
// Fit MC to extract pT and mass shapes
//-------------------------------------------------------------
//---------------------------------------------------------------------------------------//
TCanvas *printRooPlot_residual(RooPlot *frame, RooPlot *frame2, TString roohist_name[5]);
TCanvas *printRooPlot_ratio(RooPlot *frame, RooRealVar *fit_output[3], TString roohist_name[5], TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x);

void single_binned_fit_data_sample(Double_t Low_Mass = 4.0, Double_t High_Mass = 30.0, Int_t N_rebin_m = 15, Int_t N_rebin_pt = 15)
{
  gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
  gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");

  Double_t Binning_m = 260. / N_rebin_m;
  Double_t Binning_pt = 300. / N_rebin_pt;

  TFile *fIn = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/root_files/pdfMC_unbinned.root");

  RooWorkspace *w = (RooWorkspace *)fIn->Get("w");
  w->Print();

  RooRealVar *m = w->var("m");
  RooRealVar *pt = w->var("pt");
  // m->setBins(Binning_m);
  // pt->setBins(Binning_pt);

  TFile *fIn_data = new TFile("/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/HistResults_merged.root", "READ");

  TH2D *histDimuPtM_fromdata = (TH2D *)fIn_data->Get("Dimuon/CMUL7/DQ_cut_match_LT_ULS/h_PtMdiMu_CMUL7_DQ_cut_match_LT_ULS");

  TH1D *histDimuPt = histDimuPtM_fromdata->ProjectionX();
  histDimuPt->Sumw2();
  histDimuPt->Rebin(N_rebin_pt);

  TH1D *histDimuMass = histDimuPtM_fromdata->ProjectionY();
  histDimuMass->Sumw2();
  histDimuMass->Rebin(N_rebin_m);

  RooDataHist *M_Dimu_data = new RooDataHist("rooHistDimuMassFromCharm", "rooHistDimuMassFromCharm", *m, Import(*histDimuMass));
  RooDataHist *Pt_Dimu_data = new RooDataHist("rooHistDimuPtFromCharm", "rooHistDimuPtFromCharm", *pt, Import(*histDimuPt));

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

  RooAbsPdf *pdfDimuMassFromCharm = w->pdf(Form("%s", pdfMass_name[0].Data()));
  RooAbsPdf *pdfDimuPtFromCharm = w->pdf(Form("%s", pdfPt_name[0].Data()));

  RooAbsPdf *pdfDimuMassFromBeauty = w->pdf(Form("%s", pdfMass_name[1].Data()));
  RooAbsPdf *pdfDimuPtFromBeauty = w->pdf(Form("%s", pdfPt_name[1].Data()));

  RooAbsPdf *pdfDimuMassFromMixed = w->pdf(Form("%s", pdfMass_name[2].Data()));
  RooAbsPdf *pdfDimuPtFromMixed = w->pdf(Form("%s", pdfPt_name[2].Data()));

  RooPlot *frame = m->frame(Title("Imported TH1 with Poisson error bars"));
  pdfDimuMassFromCharm->plotOn(frame, LineColor(kMagenta + 2));
  pdfDimuMassFromBeauty->plotOn(frame, LineColor(kSpring - 6));
  pdfDimuMassFromMixed->plotOn(frame, LineColor(kAzure + 9));
  TCanvas *c2 = new TCanvas();
  c2->cd();
  frame->Draw();

  RooAddPdf *m_model;
  RooAddPdf *pt_model;

  RooRealVar *Mass_normForC = new RooRealVar("Mass_n_charm_output", "number dimuon from c", 28440, 0, 200000000);
  RooRealVar *Mass_normForB = new RooRealVar("Mass_n_beauty_output", "number dimuon from b", 48000, 0, 200000000);
  RooRealVar *Mass_normForMixed = new RooRealVar("Mass_n_mixed_output", "number dimuon from b,c", 3160);
  Mass_normForMixed->setConstant(kTRUE);

  m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed), RooArgList(*Mass_normForB, *Mass_normForC, *Mass_normForMixed));

  RooRealVar *Pt_normForC = new RooRealVar("Pt_n_charm_output", "number dimuon from c", 28440, 0, 200000000);
  RooRealVar *Pt_normForB = new RooRealVar("Pt_n_beauty_output", "number dimuon from b", 48000, 0, 200000000);
  RooRealVar *Pt_normForMixed = new RooRealVar("Pt_n_mixed_output", "number dimuon from b,c", 3160);
  Pt_normForMixed->setConstant(kTRUE);

  pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed), RooArgList(*Pt_normForB, *Pt_normForC, *Pt_normForMixed));

  RooFitResult *Mass_r = m_model->fitTo(*M_Dimu_data, Minimizer("Minuit2"), Save(), SumW2Error(true));
  RooFitResult *Pt_r = pt_model->fitTo(*Pt_Dimu_data, Minimizer("Minuit2"), Save(), SumW2Error(true));

  RooRealVar *Mass_fit_output[3];
  Mass_fit_output[0] = (RooRealVar *)Mass_r->floatParsFinal().find("Mass_n_charm_output");
  Mass_fit_output[1] = (RooRealVar *)Mass_r->floatParsFinal().find("Mass_n_beauty_output");
  Mass_fit_output[2] = (RooRealVar *)Mass_r->floatParsFinal().find("Mass_n_mixed_input");

  RooRealVar *Pt_fit_output[3];
  Pt_fit_output[0] = (RooRealVar *)Pt_r->floatParsFinal().find("Pt_n_charm_output");
  Pt_fit_output[1] = (RooRealVar *)Pt_r->floatParsFinal().find("Pt_n_beauty_output");
  Pt_fit_output[2] = (RooRealVar *)Pt_r->floatParsFinal().find("Pt_n_mixed_input");

  RooPlot *m_frame = m->frame(Title("m_frame"));
  m_frame->SetTitle(" ");
  m_frame->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
  m_frame->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");

  M_Dimu_data->plotOn(m_frame, Name("M_Dimu_data"), DrawOption("PEZ"));
  m_model->plotOn(m_frame, Name("pdfmass"), LineStyle(kSolid), LineColor(kRed));
  m_model->plotOn(m_frame, Name("pdfmasscharm"), Components(Form("%s", pdfMass_name[0].Data())), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
  m_model->plotOn(m_frame, Name("pdfmassbeauty"), Components(Form("%s", pdfMass_name[1].Data())), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
  m_model->plotOn(m_frame, Name("pdfmassmixed"), Components(Form("%s", pdfMass_name[2].Data())), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

  RooPlot *pt_frame = pt->frame(Title("pt_frame"));
  pt_frame->SetTitle(" ");
  pt_frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  pt_frame->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

  Pt_Dimu_data->plotOn(pt_frame, Name("Pt_Dimu_data"), DrawOption("PEZ"));
  pt_model->plotOn(pt_frame, Name("pdfpt"), LineStyle(kSolid), LineColor(kRed));
  pt_model->plotOn(pt_frame, Name("pdfptcharm"), Components(Form("%s", pdfPt_name[0].Data())), LineStyle(kDashed), LineColor(kMagenta - 2), LineWidth(5));
  pt_model->plotOn(pt_frame, Name("pdfptbeauty"), Components(Form("%s", pdfPt_name[1].Data())), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
  pt_model->plotOn(pt_frame, Name("pdfptmixed"), Components(Form("%s", pdfPt_name[2].Data())), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

  //
  // m_frame->Draw();
  TString Pt_name[5] = {"combDatapt", "pdfpt", "pdfptcharm", "pdfptbeauty", "pdfptmixed"};
  TString Mass_name[5] = {"combDatamass", "pdfmass", "pdfmasscharm", "pdfmassbeauty", "pdfmassmixed"};

  RooArgSet *m_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*m_model).snapshot(true)); // True means copy the PDF and everything it depends on
  auto &m_modelcopiedPdf = static_cast<RooAbsPdf &>((*m_modelcopyOfEverything)["m_model"]);          // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
  RooArgSet *m_modelobs = m_modelcopiedPdf.getObservables(*M_Dimu_data);
  ;

  RooArgSet *m_modelPars = m_modelcopiedPdf.getParameters(*m_modelobs);

  TF1 *m_modelFunc = m_modelcopiedPdf.asTF(*m_modelobs, *m_modelPars, *m);

  RooArgSet *pt_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*pt_model).snapshot(true)); // True means copy the PDF and everything it depends on
  auto &pt_modelcopiedPdf = static_cast<RooAbsPdf &>((*pt_modelcopyOfEverything)["pt_model"]);         // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
  RooArgSet *pt_modelobs = pt_modelcopiedPdf.getObservables(*Pt_Dimu_data);
  RooArgSet *pt_modelPars = pt_modelcopiedPdf.getParameters(*pt_modelobs);
  TF1 *pt_modelFunc = pt_modelcopiedPdf.asTF(*pt_modelobs, *pt_modelPars, *pt);

  histDimuPt->Scale(1., "width");
  histDimuPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  histDimuPt->GetYaxis()->SetTitle("d#it{N}/#it{p}_{T} (GeV/#it{c})^{-1}");

  TCanvas *pt_canvas = printRooPlot_ratio(pt_frame, Pt_fit_output, Pt_name, pt_modelFunc, histDimuPt, 0, 30);
  pt_canvas->SetName(Form("single_pt_canvas_%0.0fptbin_%0.0fmbin", Binning_pt, Binning_m));
  pt_canvas->SetTitle(Form("single_pt_canvas_%0.0fptbin_%0.0fmbin", Binning_pt, Binning_m));
  pt_canvas->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/%s_unbinned.pdf", pt_canvas->GetName()));
  pt_canvas->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/%s_unbinned.png", pt_canvas->GetName()));

  histDimuMass->Scale(1., "width");
  histDimuMass->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^2)");
  histDimuMass->GetYaxis()->SetTitle("d#it{N}/#it{m}_{#mu#mu} (GeV/#it{c}^2)^{-1}");

  TCanvas *m_canvas = printRooPlot_ratio(m_frame, Mass_fit_output, Mass_name, m_modelFunc, histDimuMass, Low_Mass, High_Mass);
  m_canvas->SetName(Form("single_m_canvas_%0.0fptbin_%0.0fmbin", Binning_pt, Binning_m));
  m_canvas->SetTitle(Form("single_m_canvas_%0.0fptbin_%0.0fmbin", Binning_pt, Binning_m));
  m_canvas->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/%s_unbinned.pdf", m_canvas->GetName()));
  m_canvas->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/%s_unbinned.png", m_canvas->GetName()));

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

TCanvas *printRooPlot_ratio(RooPlot *frame, RooRealVar *fit_output[3], TString roohist_name[5], TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x)
{
  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 1250, 850);
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
  frame->SetMaximum(1.5e+9);
  frame->SetMinimum(1.5e-3);
  frame->GetXaxis()->SetTitleOffset(1.3);
  frame->GetXaxis()->SetTitleSize(0.0475);
  frame->GetXaxis()->SetLabelSize(0.045);

  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleOffset(0.9);
  frame->GetYaxis()->SetTitleSize(0.065);
  frame->GetYaxis()->SetLabelSize(0.055);

  frame->Draw();

  TLegend *legend = new TLegend(0.675, 0.415, 1.0, 0.615);
  // legend->SetNColumns(2);
  // legend->AddEntry((TObject*)0, "", "");
  legend->SetFillStyle(0);
  legend->SetLineColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.0425);
  legend->SetHeader("Data");
  legend->SetTextAlign(11);
  legend->AddEntry("combDatamass", " ", "LP");

  TLegend *fit_legend = new TLegend(0.775, 0.415, 0.9, 0.615);
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

  letexTitle->SetTextSize(0.06);
  // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
  letexTitle->DrawLatex(0.175, 0.875, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
  letexTitle->DrawLatex(0.175, 0.785, "LHC18p period");
  // printf("WOW %s\n",roohist_name[option].Data() );
  letexTitle->SetTextSize(0.055);
  letexTitle->DrawLatex(0.625, 0.825, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f  #pm %0.3f ", fit_output[0]->getVal(), fit_output[0]->getError()));
  letexTitle->DrawLatex(0.625, 0.725, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f  #pm %0.3f ", fit_output[1]->getVal(), fit_output[1]->getError()));
  // letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.2e #pm %0.2e", fit_output[2]->getVal(), fit_output[2]->getError()));
  letexTitle->DrawLatex(0.625, 0.625, Form("#it{N}^{fixed}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %d ", 3160));
  letexTitle->SetTextSize(0.055);
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

  TLine *l1 = new TLine(minx, 0.75, 30.0, 0.75);
  l1->SetLineWidth(2);
  l1->SetLineStyle(9);
  l1->SetLineColor(kGray + 2);

  TLine *l2 = new TLine(minx, 1.25, 30.0, 1.25);
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

  h_grid_ratio->GetXaxis()->SetTitleSize(0.1);
  h_grid_ratio->GetXaxis()->SetTitleOffset(1.1);
  h_grid_ratio->GetXaxis()->SetLabelSize(0.1);
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