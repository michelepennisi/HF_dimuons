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

void binned_fit_data_sample(Int_t choice = 1, Double_t Low_Mass = 4.0, Double_t High_Mass = 30.0)
{
  gROOT->ProcessLineSync(".x ../fit_library/PtMassExpPdf.cxx+");
  gROOT->ProcessLineSync(".x ../fit_library/PtMassPol1ExpPdf.cxx+");

  TString mass_range;
  Double_t Binning_m;

  if (Low_Mass == 4 && High_Mass == 9)
  {
    mass_range.Form("_lowmass");
    Binning_m = 15;
  }

  else if (Low_Mass == 11 && High_Mass == 30)
  {
    mass_range.Form("_highmass");
    Binning_m = 19;
  }
  else if (Low_Mass == 4 && High_Mass == 30)
  {
    mass_range.Form("_y_cut");
    Binning_m = 26;
  }

  else
  {
    mass_range.Form("");
    Binning_m = 26;
  }
  Int_t choice_Mass_Cut = 0;

  if (Low_Mass == 4.0 && High_Mass == 9.0)
    choice_Mass_Cut = 1;
  else if (Low_Mass == 9.0 && High_Mass == 11.0)
    choice_Mass_Cut = 2;
  else if (Low_Mass == 11.0 && High_Mass == 15.0)
    choice_Mass_Cut = 3;
  else if (Low_Mass == 15.0 && High_Mass == 30.0)
    choice_Mass_Cut = 4;

  const Int_t n_Mass_Cut = 6;
  TString name_DataMass_Cut[n_Mass_Cut];
  name_DataMass_Cut[0].Form("DQ_cut_match_LT_ULS");
  name_DataMass_Cut[1].Form("DQ_cut_match_LT_LowMass_ULS");
  name_DataMass_Cut[2].Form("DQ_cut_match_LT_Yres_ULS");
  name_DataMass_Cut[3].Form("DQ_cut_match_LT_InterMass_ULS");
  name_DataMass_Cut[4].Form("DQ_cut_match_LT_HighMass_ULS");

  TFile *fIn = new TFile("rooWorkspace_4_30_Rec_DQ_cut_match_LT.root");
  

  RooWorkspace *w = (RooWorkspace *)fIn->Get("w");
  w->Print();

  RooRealVar *m = w->var("m");
  RooRealVar *pt = w->var("pt");
  m->setBins(20);
  pt->setBins(20);

  TFile *fIn_data = new TFile("/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Hist_AOD/15_10_2022/HistResults_merged.root", "READ");

  TH2D *histDimuPtM_fromdata = (TH2D *)fIn_data->Get("Dimuon/CMUL7/DQ_cut_match_AT_ULS/h_PtMdiMu_CMUL7_DQ_cut_match_AT_ULS_norm");

  TH1D *histDimuPt = histDimuPtM_fromdata->ProjectionX();
  // histDimuPt->Rebin(30);
  histDimuPt->Scale(1., "width");

  TH1D *histDimuMass = histDimuPtM_fromdata->ProjectionY();
  // histDimuMass->Rebin(Binning_m);
  histDimuMass->Scale(1., "width");
  histDimuMass->Sumw2();

  RooDataHist *M_Dimu_data = new RooDataHist("rooHistDimuMassFromCharm", "rooHistDimuMassFromCharm", *m, Import(*histDimuMass));
  RooDataHist *Pt_Dimu_data = new RooDataHist("rooHistDimuPtFromCharm", "rooHistDimuPtFromCharm", *pt, Import(*histDimuPt));

  RooCategory sample("sample", "sample");
  sample.defineType("mass");
  sample.defineType("transversemomentum");

  RooDataHist *combData_set = new RooDataHist("combData", "combined data", RooArgSet(*m, *pt), Index(sample), Import("mass", *M_Dimu_data), Import("transversemomentum", *Pt_Dimu_data));

  // RooPlot *_frame = pt->frame(Title("pt_frame")) ;
  // combData->plotOn(_frame,Name("combDatapt"),Cut("sample==sample::transversemomentum"),DrawOption("PEZ")) ;
  // new TCanvas();
  // _frame->Draw();
  // hDimuPt_data_rebin->SetMarkerColor(kGreen);
  // hDimuPt_data_rebin->SetMarkerSize(1.2);
  // hDimuPt_data_rebin->Scale(1.,"width");
  // hDimuPt_data_rebin->Draw("PESAME");
  // return;
  // Fix parameters to the MCRooSimultaneous
  RooRealVar *B_DimuMassFromCharm = w->var("B_DimuMassFromCharm");
  B_DimuMassFromCharm->setConstant(kTRUE);
  RooRealVar *n1_DimuMassFromCharm = w->var("n1_DimuMassFromCharm");
  n1_DimuMassFromCharm->setConstant(kTRUE);
  RooRealVar *n2_DimuMassFromCharm = w->var("n2_DimuMassFromCharm");
  n2_DimuMassFromCharm->setConstant(kTRUE);
  RooRealVar *C_DimuMassFromCharm = w->var("C_DimuMassFromCharm");
  C_DimuMassFromCharm->setConstant(kTRUE);
  RooRealVar *n3_DimuMassFromCharm = w->var("n3_DimuMassFromCharm");
  n3_DimuMassFromCharm->setConstant(kTRUE);

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
  pdfMass_name[0].Form("pdfDimuMassFromCharm");
  pdfMass_name[1].Form("pdfDimuMassFromBeauty");
  pdfMass_name[2].Form("pdfDimuMassFromMixed");

  pdfPt_name[0].Form("pdfDimuPtFromCharm");
  pdfPt_name[1].Form("pdfDimuPtFromBeauty");
  pdfPt_name[2].Form("pdfDimuPtFromMixed");

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

  if (choice == 0)
  {
    ////////////////////
    // Fit with normalization
    RooRealVar *normForC = new RooRealVar("n_charm_output", "number dimuon from c", 25000, 0, 2000000);
    RooRealVar *normForB = new RooRealVar("n_beauty_output", "number dimuon from b", 51350, 40000, 2000000);
    RooRealVar *normForMixed = new RooRealVar("n_mixed_output", "number dimuon from b,c", 3160);
    normForMixed->setConstant(kTRUE);

    m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
    pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
  }

  else if (choice == 1)
  {
    ////////////////////
    // Fit with fraction
    RooRealVar *normForC = new RooRealVar("fr_charm_output", "fraction dimuon from c", 0.35, 0., 1.);
    RooRealVar *normForB = new RooRealVar("fr_beauty_output", "fraction dimuon from b", 0.6, 0., 1.);
    // RooRealVar *normForMixed = new RooRealVar("fr_mixed_output", "fraction dimuon from c", 0.05);
    // normForMixed->setConstant(kTRUE);

    RooFormulaVar *normForMixed = new RooFormulaVar("fr_mixed_output", "1-@0-@1", RooArgList(*normForC, *normForB));

    m_model = new RooAddPdf("m_model", " dimuMassFromB+dimuMassFromC + dimuMassFromMixed", RooArgList(*pdfDimuMassFromCharm, *pdfDimuMassFromBeauty, *pdfDimuMassFromMixed), RooArgList(*normForC, *normForB));
    pt_model = new RooAddPdf("pt_model", "dimuMassFromB+dimuMassFromC + dimuPtFromMixed", RooArgList(*pdfDimuPtFromCharm, *pdfDimuPtFromBeauty, *pdfDimuPtFromMixed), RooArgList(*normForC, *normForB));

    // m_model = new RooAddPdf("m_model", " dimuMassFromB+dimuMassFromC + dimuMassFromMixed", RooArgList(*pdfDimuMassFromCharm, *pdfDimuMassFromBeauty), RooArgList(*normForC, *normForB));
    // pt_model = new RooAddPdf("pt_model", "dimuMassFromB+dimuMassFromC + dimuPtFromMixed", RooArgList(*pdfDimuPtFromCharm, *pdfDimuPtFromBeauty), RooArgList(*normForC, *normForB));

    // RooAddPdf *m_model_signal=new RooAddPdf("m_model_signal", " (dimuMassFromB+dimuMassFromC)", RooArgList(*pdfDimuMassFromCharm,*pdfDimuMassFromMixed), RooArgList(*normForC));
    // RooAddPdf *pt_model_signal=new RooAddPdf("pt_model_signal", " (dimuPtFromB+dimuPtFromC)", RooArgList(*pdfDimuPtFromCharm,*pdfDimuPtFromMixed), RooArgList(*normForC));
    //
    // m_model = new RooAddPdf("m_model", " (dimuMassFromB+dimuMassFromC) + dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty,*m_model_signal), RooArgList(*normForB));
    // pt_model = new RooAddPdf("pt_model", "(dimuMassFromB+dimuMassFromC) + dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty,*pt_model_signal), RooArgList(*normForB));
  }
  // m->setRange("region1", 4, 9);
  // m->setRange("region2", 11, 30);
  // Define the pdf for the simultaneous fit

  RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
  simPdf.addPdf(*m_model, "mass");
  simPdf.addPdf(*pt_model, "transversemomentum");
  // simPdf.fitTo(*combData);
  //
  RooFitResult *r = simPdf.fitTo(*combData_set, Save(), SumW2Error(true));
  ;
  // r = simPdf.fitTo(*combData, Save());
  r->floatParsFinal().Print("s");

  RooRealVar *fit_output[3];
  if (choice == 0)
  {
    fit_output[0] = (RooRealVar *)r->floatParsFinal().find("n_charm_output");
    fit_output[1] = (RooRealVar *)r->floatParsFinal().find("n_beauty_output");
    fit_output[2] = (RooRealVar *)r->floatParsFinal().find("n_mixed_input");
  }
  else if (choice == 1)
  {
    fit_output[0] = (RooRealVar *)r->floatParsFinal().find("fr_charm_output");
    printf("OOOOOOOOO %0.5f", fit_output[0]->getError());
    fit_output[1] = (RooRealVar *)r->floatParsFinal().find("fr_beauty_output");
    printf("OOOOOOOOO %0.5f", fit_output[1]->getError());
    fit_output[2] = (RooRealVar *)r->floatParsFinal().find("fr_mixed_output");
  }

  // fit_output[2] = 0;
  RooPlot *m_frame = m->frame(Title("m_frame"));
  m_frame->SetTitle(" ");
  m_frame->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
  m_frame->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");

  combData_set->plotOn(m_frame, Name("combDatamass"), Cut("sample==sample::mass"), DrawOption("PEZ"));
  simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *combData_set), LineStyle(kSolid), LineColor(kRed));
  simPdf.plotOn(m_frame, Name("pdfmasscharm"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[0].Data())), ProjWData(sample, *combData_set), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
  simPdf.plotOn(m_frame, Name("pdfmassbeauty"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[1].Data())), ProjWData(sample, *combData_set), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
  simPdf.plotOn(m_frame, Name("pdfmassmixed"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[2].Data())), ProjWData(sample, *combData_set), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

  RooPlot *pt_frame = pt->frame(Title("pt_frame"));
  pt_frame->SetTitle(" ");
  pt_frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  pt_frame->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

  combData_set->plotOn(pt_frame, Name("combDatapt"), Cut("sample==sample::transversemomentum"), DrawOption("PEZ"));
  simPdf.plotOn(pt_frame, Name("pdfpt"), Slice(sample, "transversemomentum"), ProjWData(sample, *combData_set), LineStyle(kSolid), LineColor(kRed));
  simPdf.plotOn(pt_frame, Name("pdfptcharm"), Slice(sample, "transversemomentum"), Components(Form("%s", pdfPt_name[0].Data())), ProjWData(sample, *combData_set), LineStyle(kDashed), LineColor(kMagenta - 2), LineWidth(5));
  simPdf.plotOn(pt_frame, Name("pdfptbeauty"), Slice(sample, "transversemomentum"), Components(Form("%s", pdfPt_name[1].Data())), ProjWData(sample, *combData_set), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
  simPdf.plotOn(pt_frame, Name("pdfptmixed"), Slice(sample, "transversemomentum"), Components(Form("%s", pdfPt_name[2].Data())), ProjWData(sample, *combData_set), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

  //
  // m_frame->Draw();
  TString Pt_name[5] = {"combDatapt", "pdfpt", "pdfptcharm", "pdfptbeauty", "pdfptmixed"};
  TString Mass_name[5] = {"combDatamass", "pdfmass", "pdfmasscharm", "pdfmassbeauty", "pdfmassmixed"};

  // pt_frame->SetMaximum(hDimuPt_data->GetMaximum() * 1000);
  // pt_frame->SetMinimum(2e-14);

  RooArgSet *m_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*m_model).snapshot(true)); // True means copy the PDF and everything it depends on
  auto &m_modelcopiedPdf = static_cast<RooAbsPdf &>((*m_modelcopyOfEverything)["m_model"]);          // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
  RooArgSet *m_modelobs = m_modelcopiedPdf.getObservables(*combData_set);
  ;

  RooArgSet *m_modelPars = m_modelcopiedPdf.getParameters(*m_modelobs);

  TF1 *m_modelFunc = m_modelcopiedPdf.asTF(*m_modelobs, *m_modelPars, *m);

  RooArgSet *pt_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*pt_model).snapshot(true)); // True means copy the PDF and everything it depends on
  auto &pt_modelcopiedPdf = static_cast<RooAbsPdf &>((*pt_modelcopyOfEverything)["pt_model"]);         // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
  RooArgSet *pt_modelobs = pt_modelcopiedPdf.getObservables(*combData_set);
  RooArgSet *pt_modelPars = pt_modelcopiedPdf.getParameters(*pt_modelobs);
  TF1 *pt_modelFunc = pt_modelcopiedPdf.asTF(*pt_modelobs, *pt_modelPars, *pt);

  TH1 *hDimuPt_data = Pt_Dimu_data->createHistogram("h_ptdata", *pt, Binning(30, 0, 30));
  hDimuPt_data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hDimuPt_data->GetYaxis()->SetTitle("d#it{N}/#it{p}_{T} (GeV/#it{c})^{-1}");

  TCanvas *pt_canvas = printRooPlot_ratio(pt_frame, fit_output, Pt_name, pt_modelFunc, hDimuPt_data, 0, 30);
  pt_canvas->SetName(Form("pt_canvas_opt%d_%s", choice, mass_range.Data()));
  pt_canvas->SetTitle(Form("pt_canvas_opt%d_%s", choice, mass_range.Data()));
  pt_canvas->SaveAs(Form("plot/%s_unbinned.pdf", pt_canvas->GetName()));
  // m_frame->SetMaximum(hDimuM_data->GetMaximum() * 1000);
  // m_frame->SetMinimum(2e-14);
  TH1 *hDimuM_data = M_Dimu_data->createHistogram("h_mdata", *m, Binning(Binning_m, Low_Mass, High_Mass));
  hDimuM_data->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^2)");
  hDimuM_data->GetYaxis()->SetTitle("d#it{N}/#it{m}_{#mu#mu} (GeV/#it{c}^2)^{-1}");
  TCanvas *m_canvas = printRooPlot_ratio(m_frame, fit_output, Mass_name, m_modelFunc, hDimuM_data, Low_Mass, High_Mass);
  m_canvas->SetName(Form("m_canvas_opt%d%s", choice, mass_range.Data()));
  m_canvas->SetTitle(Form("m_canvas_opt%d%s", choice, mass_range.Data()));
  m_canvas->SaveAs(Form("plot/%s_unbinned.pdf", m_canvas->GetName()));

  // printf("fraction c %0.3f error fr c %0.10f",fit_output[1]->getVal(), (RooRealVar *)r->floatParsFinal().find("fr_beauty_output").getError());
  // for (auto &p : r->floatParsFinal())
  // {
  //   auto v = (RooRealVar *)p;
  //   std::cout << "parameter " << v->GetName() << " : " << v->getVal() << " +/- " << v->getError() << std::endl;
  // }
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
  frame->SetMaximum(1.5e+8);
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
  letexTitle->DrawLatex(0.675, 0.825, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f %% #pm %0.3f %%", 100 * fit_output[0]->getVal(), 100 * fit_output[0]->getError()));
  letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f %% #pm %0.3f %%", 100 * fit_output[1]->getVal(), 100 * fit_output[1]->getError()));
  // letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.2e #pm %0.2e", fit_output[2]->getVal(), fit_output[2]->getError()));
  letexTitle->DrawLatex(0.675, 0.625, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f %%", 100 - 100 * fit_output[0]->getVal() - 100 * fit_output[1]->getVal()));
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