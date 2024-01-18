#include "/home/michele_pennisi/cernbox/common_include.h"

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

TString Generator = "Powheg+PYTHIA6";
Bool_t With_DY = kTRUE, Plot_Likehood = kFALSE;

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

void unbinned_fit_data_sample(Int_t choice = 0, Double_t Low_Mass = 4.0, Double_t High_Mass = 30.0, Double_t Low_Pt = 0.0, Double_t High_Pt = 30.0, Int_t Binning_m = 26, Int_t Binning_pt = 30)
{

  gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
  gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");

  TString mass_range_data;

  TString mass_range;

  if (Low_Mass == 4 && High_Mass == 9)
  {
    if (High_Pt == 10)
    {
      mass_range_data.Form("_LowMass_LowPt");
      mass_range.Form("_LowMass_LowPt");
    }
    else if (High_Pt == 30)
    {
      mass_range_data.Form("_LowMass");
      mass_range.Form("_LowMass");
    }
  }

  else if (Low_Mass == 11 && High_Mass == 30)
  {
    mass_range_data.Form("_HighMass");
    mass_range.Form("_HighMass");
  }
  else if (Low_Mass == 4 && High_Mass == 30)
  {
    mass_range_data.Form("_withcut");
    mass_range.Form("_PowhegOnly_");
  }
  TFile *fIn;
  if (Generator.Contains("Powheg"))
  {

    fIn = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/Powheg_pdfMC_unbinned%s.root", mass_range.Data()));
    printf("%s\n", Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/Powheg_pdfMC_unbinned%s.root", mass_range.Data()));
  }
  else
  {
    printf("%s\n", Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/pdfMC_unbinned%s.root", mass_range.Data()));
    fIn = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/pdfMC_unbinned%s.root", mass_range.Data()));
  }

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
  TTree *tree_data = (TTree *)fIn_data->Get(Form("rec_data_tree%s", mass_range_data.Data()));
  RooDataSet *unbinned_M_Dimu_data = new RooDataSet("M_Dimu_data", "M_Dimu_data", RooArgSet(*m), Import(*tree_data), Cut("m<9 || (m>11 && m<30)"));
  RooDataSet *unbinned_Pt_Dimu_data = new RooDataSet("Pt_Dimu_data", "Pt_Dimu_data", RooArgSet(*pt), Import(*tree_data)); // 5
  RooDataSet *unbinned_combData_set = new RooDataSet("combData", "combined data", RooArgSet(*m, *pt), Index(sample), Import("mass", *unbinned_M_Dimu_data), Import("transversemomentum", *unbinned_Pt_Dimu_data));

  TFile *fIn_DY = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/root_files/pdfMC_unbinnedDY.root");
  // RooWorkspace *w = (RooWorkspace *)fIn_DY->Get("w_weighted_mode2");
  fIn_DY->ls();
  RooWorkspace *w_DY = (RooWorkspace *)fIn_DY->Get("w");
  w_DY->Print();

  // Fix parameters to the MCRooSimultaneous
  RooRealVar *B_DimuMassFromCharm = w->var("B_DimuMassFromCharm");
  B_DimuMassFromCharm->setConstant(kTRUE);
  RooRealVar *n1_DimuMassFromCharm = w->var("n1_DimuMassFromCharm");
  n1_DimuMassFromCharm->setConstant(kTRUE);
  RooRealVar *n2_DimuMassFromCharm = w->var("n2_DimuMassFromCharm");
  n2_DimuMassFromCharm->setConstant(kTRUE);

  RooRealVar *B_DimuPtFromCharm = w->var("B_DimuPtFromCharm");
  B_DimuPtFromCharm->setConstant(kTRUE);
  RooRealVar *n1_DimuPtFromCharm = w->var("n1_DimuPtFromCharm");
  n1_DimuPtFromCharm->setConstant(kTRUE);
  RooRealVar *n2_DimuPtFromCharm = w->var("n2_DimuPtFromCharm");
  n2_DimuPtFromCharm->setConstant(kTRUE);

  RooRealVar *B_DimuMassFromBeauty = w->var("B_DimuMassFromBeauty");
  B_DimuMassFromBeauty->setConstant(kTRUE);
  RooRealVar *n1_DimuMassFromBeauty = w->var("n1_DimuMassFromBeauty");
  n1_DimuMassFromBeauty->setConstant(kTRUE);
  RooRealVar *n2_DimuMassFromBeauty = w->var("n2_DimuMassFromBeauty");
  n2_DimuMassFromBeauty->setConstant(kTRUE);

  RooRealVar *B_DimuPtFromBeauty = w->var("B_DimuPtFromBeauty");
  B_DimuPtFromBeauty->setConstant(kTRUE);
  RooRealVar *n1_DimuPtFromBeauty = w->var("n1_DimuPtFromBeauty");
  n1_DimuPtFromBeauty->setConstant(kTRUE);
  RooRealVar *n2_DimuPtFromBeauty = w->var("n2_DimuPtFromBeauty");
  n2_DimuPtFromBeauty->setConstant(kTRUE);

  // RooRealVar *C_DimuPtFromCharm = w->var("C_DimuPtFromCharm");
  // C_DimuPtFromCharm->setConstant(kTRUE);
  // RooRealVar *n3_DimuPtFromCharm = w->var("n3_DimuPtFromCharm");
  // n3_DimuPtFromCharm->setConstant(kTRUE);

  TString Mixed_suffix;
  if (Generator.Contains("Powheg"))
    Mixed_suffix.Form("HF_Mixed");
  else
    Mixed_suffix.Form("Mixed");

  RooRealVar *B_DimuMassFromMixed = w->var(Form("B_DimuMassFrom%s", Mixed_suffix.Data()));
  B_DimuMassFromMixed->setConstant(kTRUE);
  RooRealVar *n1_DimuMassFromMixed = w->var(Form("n1_DimuMassFrom%s", Mixed_suffix.Data()));
  n1_DimuMassFromMixed->setConstant(kTRUE);
  RooRealVar *n2_DimuMassFromMixed = w->var(Form("n2_DimuMassFrom%s", Mixed_suffix.Data()));
  n2_DimuMassFromMixed->setConstant(kTRUE);

  RooRealVar *B_DimuPtFromMixed = w->var(Form("B_DimuPtFrom%s", Mixed_suffix.Data()));
  B_DimuPtFromMixed->setConstant(kTRUE);
  RooRealVar *n1_DimuPtFromMixed = w->var(Form("n1_DimuPtFrom%s", Mixed_suffix.Data()));
  n1_DimuPtFromMixed->setConstant(kTRUE);
  RooRealVar *n2_DimuPtFromMixed = w->var(Form("n2_DimuPtFrom%s", Mixed_suffix.Data()));
  n2_DimuPtFromMixed->setConstant(kTRUE);

  TString pdfMass_name[3];
  TString pdfPt_name[3];
  pdfMass_name[0].Form("pdfDimuMassFromcharm");
  pdfMass_name[1].Form("pdfDimuMassFrombeauty");
  pdfMass_name[2].Form("pdfDimuMassFrommixed");

  pdfPt_name[0].Form("pdfDimuPtFromcharm");
  pdfPt_name[1].Form("pdfDimuPtFrombeauty");
  pdfPt_name[2].Form("pdfDimuPtFrommixed");

  RooRealVar *B_DimuMassFromDY = w_DY->var("B_DimuMassFromDY");
  B_DimuMassFromDY->setConstant(kTRUE);
  RooRealVar *n1_DimuMassFromDY = w_DY->var("n1_DimuMassFromDY");
  n1_DimuMassFromDY->setConstant(kTRUE);
  RooRealVar *n2_DimuMassFromDY = w_DY->var("n2_DimuMassFromDY");
  n2_DimuMassFromDY->setConstant(kTRUE);

  RooRealVar *B_DimuPtFromDY = w_DY->var("B_DimuPtFromDY");
  B_DimuPtFromDY->setConstant(kTRUE);
  RooRealVar *n1_DimuPtFromDY = w_DY->var("n1_DimuPtFromDY");
  n1_DimuPtFromDY->setConstant(kTRUE);
  RooRealVar *n2_DimuPtFromDY = w_DY->var("n2_DimuPtFromDY");
  n2_DimuPtFromDY->setConstant(kTRUE);
  RooAbsPdf *pdfDimuMassFromDY = w_DY->pdf("pdfDimuMassFromDY");
  RooAbsPdf *pdfDimuPtFromDY = w_DY->pdf("pdfDimuPtFromDY");

  RooAbsPdf *pdfDimuMassFromCharm = w->pdf("pdfDimuMassFromcharm");
  RooAbsPdf *pdfDimuPtFromCharm = w->pdf("pdfDimuPtFromcharm");

  RooAbsPdf *pdfDimuMassFromBeauty = w->pdf("pdfDimuMassFrombeauty");
  RooAbsPdf *pdfDimuPtFromBeauty = w->pdf("pdfDimuPtFrombeauty");

  RooAbsPdf *pdfDimuMassFromMixed = w->pdf("pdfDimuMassFrommixed");
  RooAbsPdf *pdfDimuPtFromMixed = w->pdf("pdfDimuPtFrommixed");
  // RooPlot *frame = pt->frame(Title("Imported TH1 with Poisson error bars"));
  // pdfDimuPtFromCharm->plotOn(frame, LineColor(kMagenta + 2));
  // pdfDimuPtFromBeauty->plotOn(frame, LineColor(kSpring - 6));
  // pdfDimuPtFromMixed->plotOn(frame, LineColor(kAzure + 9));
  // TCanvas *c2 = new TCanvas();
  // c2->cd();
  // frame->Draw();

  RooAddPdf *m_model;
  RooAddPdf *m_model_HF;
  RooAddPdf *pt_model;
  RooAddPdf *pt_model_HF;

  RooRealVar *normForC;
  RooRealVar *normForB;
  RooRealVar *normForMixed;
  RooRealVar *normForDY;

  if (choice == 0)
  {
    ////////////////////
    // Fit with normalization

    if (Low_Mass == 4 && High_Mass == 9)
    {
      normForC = new RooRealVar("n_charm_output", "number dimuon from c", 28440, 0, 200000);
      normForB = new RooRealVar("n_beauty_output", "number dimuon from b", 48000, 0, 200000);
      normForMixed = new RooRealVar("n_mixed_output", "number dimuon from b,c", (3.6 / 100) * tree_data->GetEntries());
      normForMixed->setConstant(kTRUE);
    }
    else if (Low_Mass == 4 && High_Mass == 30)
    {
      normForC = new RooRealVar("n_charm_output", "number dimuon from c", 28440, 0, 200000);
      normForB = new RooRealVar("n_beauty_output", "number dimuon from b", 48000, 0, 200000);
      normForMixed = new RooRealVar("n_mixed_output", "number dimuon from b,c", (2.6 / 100) * tree_data->GetEntries());
      normForMixed->setConstant(kTRUE);
    }
    else if (Low_Mass == 11 && High_Mass == 30)
    {
      normForC = new RooRealVar("n_charm_output", "number dimuon from c", 1391, 0, 3000);
      normForB = new RooRealVar("n_beauty_output", "number dimuon from b", 345, 0, 3000);
      normForMixed = new RooRealVar("n_mixed_output", "number dimuon from b,c", 30);
      normForMixed->setConstant(kTRUE);
    }
    if (With_DY)
    {
      normForDY = new RooRealVar("n_DY_output", "number dimuon from DY", 5000, 150, 200000);
      m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed, *pdfDimuMassFromDY), RooArgList(*normForB, *normForC, *normForMixed, *normForDY));
      pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed, *pdfDimuPtFromDY), RooArgList(*normForB, *normForC, *normForMixed, *normForDY));
    }
    else
    {
      m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
      pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
    }
  }

  else if (choice == 1)
  {
    ////////////////////
    // Fit with fraction
    normForC = new RooRealVar("fr_charm_output", "fraction dimuon from c", 0.35, 0., 1.);
    normForB = new RooRealVar("fr_beauty_output", "fraction dimuon from b", 0.605, 0., 1.);

    RooRealVar *normForMixed = new RooRealVar("fr_mixed_output", "fraction dimuon from c", 0.036);
    normForMixed->setConstant(kTRUE);

    if (With_DY)
    {
      m_model_HF = new RooAddPdf("m_model_HF", "dimuMassFromMixed +dimuMassFromB + dimuMassFromC", RooArgList(*pdfDimuMassFromMixed, *pdfDimuMassFromCharm, *pdfDimuMassFromBeauty), RooArgList(*normForMixed, *normForC));
      pt_model_HF = new RooAddPdf("pt_model_HF", "dimuPtFromMixed + dimuMassFromB + dimuMassFromC", RooArgList(*pdfDimuPtFromMixed, *pdfDimuPtFromCharm, *pdfDimuPtFromBeauty), RooArgList(*normForMixed, *normForC));
      normForDY = new RooRealVar("fr_beauty_output", "fraction dimuon from b", 0.05, 0., 1.);

      m_model = new RooAddPdf("m_model", " dimuMassFromDY + HF", RooArgList(*pdfDimuMassFromDY, *m_model_HF), RooArgList(*normForDY));
      pt_model = new RooAddPdf("pt_model", "dimuMassFromDY + HF", RooArgList(*pdfDimuPtFromDY, *pt_model_HF), RooArgList(*normForDY));
    }
    else
    {
      m_model_HF = new RooAddPdf("m_model_HF", "dimuMassFromB + dimuMassFromC", RooArgList(*pdfDimuMassFromCharm, *pdfDimuMassFromBeauty), RooArgList(*normForC));
      pt_model_HF = new RooAddPdf("pt_model_HF", "dimuMassFromB + dimuMassFromC", RooArgList(*pdfDimuPtFromCharm, *pdfDimuPtFromBeauty), RooArgList(*normForC));

      m_model = new RooAddPdf("m_model", " dimuMassFromB+dimuMassFromC + dimuMassFromMixed", RooArgList(*pdfDimuMassFromMixed, *m_model_HF), RooArgList(*normForMixed));
      pt_model = new RooAddPdf("pt_model", "dimuMassFromB+dimuMassFromC + dimuPtFromMixed", RooArgList(*pdfDimuPtFromMixed, *pt_model_HF), RooArgList(*normForMixed));
    }
  }

  // Define the pdf for the simultaneous fit

  RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
  simPdf.addPdf(*m_model, "mass");
  simPdf.addPdf(*pt_model, "transversemomentum");
  simPdf.Print("t");
  RooFitResult *r;
  if (Low_Mass == 4 && High_Mass == 30)
  {
    m->setRange("pluto", 4, 9);
    pt->setRange("pluto", 20, 30);
    m->setRange("pippo", 11, 30);
    pt->setRange("pippo", 0, 20);
    r = simPdf.fitTo(*unbinned_combData_set, Range("pluto,pippo"), SumCoefRange("pluto,pippo"), SplitRange(false), Save(), SumW2Error(true));
  }
  else
    r = simPdf.fitTo(*unbinned_combData_set, Minimizer("Minuit2"), Save(), SumW2Error(true));

  RooRealVar *fit_output[4];
  if (choice == 0)
  {
    fit_output[0] = (RooRealVar *)r->floatParsFinal().find("n_charm_output");
    fit_output[1] = (RooRealVar *)r->floatParsFinal().find("n_beauty_output");
    fit_output[2] = normForMixed;
    if (With_DY)
      fit_output[3] = (RooRealVar *)r->floatParsFinal().find("n_DY_output");
  }
  else if (choice == 1)
  {
    // fit_output[0] = (RooRealVar *)r->floatParsFinal().find("fr_charm_output");
    // printf("OOOOOOOOO %0.5f", fit_output[0]->getError());
    // fit_output[1] = (RooRealVar *)r->floatParsFinal().find("fr_beauty_output");
    // printf("OOOOOOOOO %0.5f", fit_output[1]->getError());
    // fit_output[2] = (RooRealVar *)r->floatParsFinal().find("fr_mixed_output");
  }

  // fit_output[2] = 0;
  RooPlot *m_frame = m->frame(Title("m_frame"));
  m_frame->SetTitle(" ");
  m_frame->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
  m_frame->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");

  RooPlot *pt_frame = pt->frame(Title("pt_frame"));
  pt_frame->SetTitle(" ");
  pt_frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  pt_frame->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
  if (Low_Mass == 4 && High_Mass == 30)
  {
    unbinned_combData_set->plotOn(m_frame, Name("combDatamass"), Cut("sample==sample::mass"), Range("pluto,pippo"), NormRange("pluto,pippo"), DrawOption("PEZ"));
    simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *unbinned_combData_set), Range("pluto,pippo"), NormRange("pluto,pippo"), LineStyle(kSolid), LineColor(kRed));
    simPdf.plotOn(m_frame, Name("pdfmasscharm"), Slice(sample, "mass"), Components("pdfDimuMassFromcharm"), ProjWData(sample, *unbinned_combData_set), Range("pluto,pippo"), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
    simPdf.plotOn(m_frame, Name("pdfmassbeauty"), Slice(sample, "mass"), Components("pdfDimuMassFrombeauty"), ProjWData(sample, *unbinned_combData_set), Range("pluto,pippo"), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
    simPdf.plotOn(m_frame, Name("pdfmassmixed"), Slice(sample, "mass"), Components("pdfDimuMassFrommixed"), ProjWData(sample, *unbinned_combData_set), Range("pluto,pippo"), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

    unbinned_combData_set->plotOn(pt_frame, Name("combDatapt"), Cut("sample==sample::transversemomentum"), DrawOption("PEZ"));
    simPdf.plotOn(pt_frame, Name("pdfpt"), Slice(sample, "transversemomentum"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
    simPdf.plotOn(pt_frame, Name("pdfptcharm"), Slice(sample, "transversemomentum"), Components("pdfDimuPtFromcharm"), ProjWData(sample, *unbinned_combData_set), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(kMagenta - 2), LineWidth(5));
    simPdf.plotOn(pt_frame, Name("pdfptbeauty"), Slice(sample, "transversemomentum"), Components("pdfDimuPtFrombeauty"), ProjWData(sample, *unbinned_combData_set), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
    simPdf.plotOn(pt_frame, Name("pdfptmixed"), Slice(sample, "transversemomentum"), Components("pdfDimuPtFrommixed"), ProjWData(sample, *unbinned_combData_set), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));
  }
  else
  {
    unbinned_combData_set->plotOn(m_frame, Name("combDatamass"), Cut("sample==sample::mass"), DrawOption("PEZ"));
    simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
    simPdf.plotOn(m_frame, Name("pdfmasscharm"), Slice(sample, "mass"), Components("pdfDimuMassFromcharm"), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
    simPdf.plotOn(m_frame, Name("pdfmassbeauty"), Slice(sample, "mass"), Components("pdfDimuMassFrombeauty"), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
    simPdf.plotOn(m_frame, Name("pdfmassmixed"), Slice(sample, "mass"), Components("pdfDimuMassFrommixed"), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

    unbinned_combData_set->plotOn(pt_frame, Name("combDatapt"), Cut("sample==sample::transversemomentum"), DrawOption("PEZ"));
    simPdf.plotOn(pt_frame, Name("pdfpt"), Slice(sample, "transversemomentum"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
    simPdf.plotOn(pt_frame, Name("pdfptcharm"), Slice(sample, "transversemomentum"), Components("pdfDimuPtFromcharm"), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(kMagenta - 2), LineWidth(5));
    simPdf.plotOn(pt_frame, Name("pdfptbeauty"), Slice(sample, "transversemomentum"), Components("pdfDimuPtFrombeauty"), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
    simPdf.plotOn(pt_frame, Name("pdfptmixed"), Slice(sample, "transversemomentum"), Components("pdfDimuPtFrommixed"), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));
  }

  if (With_DY)
  {
    simPdf.plotOn(m_frame, Name("pdfmassDY"), Slice(sample, "mass"), Components("pdfDimuMassFromDY"), ProjWData(sample, *unbinned_combData_set), Range("pluto,pippo"), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(kBlack), LineWidth(5));
    simPdf.plotOn(pt_frame, Name("pdfptDY"), Slice(sample, "transversemomentum"), Components("pdfDimuPtFromDY"), ProjWData(sample, *unbinned_combData_set), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(kBlack), LineWidth(5));
  }

  TString Pt_name[6] = {"combDatapt", "pdfpt", "pdfptcharm", "pdfptbeauty", "pdfptmixed", "pdfptDY"};

  TString Mass_name[6] = {"combDatamass", "pdfmass", "pdfmasscharm", "pdfmassbeauty", "pdfmassmixed", "pdfmassDY"};

  // Save the fit as TF1

  RooArgSet *m_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*m_model).snapshot(true)); // True means copy the PDF and everything it depends on
  auto &m_modelcopiedPdf = static_cast<RooAbsPdf &>((*m_modelcopyOfEverything)["m_model"]);          // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
  RooArgSet *m_modelobs = m_modelcopiedPdf.getObservables(*unbinned_combData_set);
  RooArgSet *m_modelPars = m_modelcopiedPdf.getParameters(*m_modelobs);
  TF1 *m_modelFunc = m_modelcopiedPdf.asTF(*m_modelobs, *m_modelPars, *m);

  RooArgSet *pt_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*pt_model).snapshot(true)); // True means copy the PDF and everything it depends on
  auto &pt_modelcopiedPdf = static_cast<RooAbsPdf &>((*pt_modelcopyOfEverything)["pt_model"]);         // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
  RooArgSet *pt_modelobs = pt_modelcopiedPdf.getObservables(*unbinned_combData_set);
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
  if (With_DY)
    info.Form("unbinned_data_optfit%d_M%0.0f_%0.0f_BinM%d_BinPt%d_withDY", choice, Low_Mass, High_Mass, Binning_m, Binning_pt);
  else
    info.Form("unbinned_data_optfit%d_M%0.0f_%0.0f_BinM%d_BinPt%d_noDY", choice, Low_Mass, High_Mass, Binning_m, Binning_pt);
  pt_frame->SetMaximum(1.2e+8);
  pt_frame->SetMinimum(1.5e-1);
  TCanvas *pt_canvas = printRooPlot_ratio(pt_frame, false, fit_output, choice, Pt_name, pt_modelFunc, hDimuPt_data, 0, High_Pt);
  pt_canvas->SetName(Form("pt_canvas_%s", info.Data()));
  pt_canvas->SetTitle(Form("pt_canvas_%s", info.Data()));
  pt_canvas->SaveAs(Form("images/%s_%s.pdf", Generator.Data(), pt_canvas->GetName()));
  pt_canvas->SaveAs(Form("images/%s_%s.png", Generator.Data(), pt_canvas->GetName()));

  m_frame->SetMaximum(1.2e+9);
  m_frame->SetMinimum(1.5e-1);
  TCanvas *m_canvas = printRooPlot_ratio(m_frame, false, fit_output, choice, Mass_name, m_modelFunc, hDimuM_data, Low_Mass, High_Mass);
  m_canvas->SetName(Form("m_canvas_%s", info.Data()));
  m_canvas->SetTitle(Form("m_canvas_%s", info.Data()));
  m_canvas->SaveAs(Form("images/%s_%s.pdf", Generator.Data(), m_canvas->GetName()));
  m_canvas->SaveAs(Form("images/%s_%s.png", Generator.Data(), m_canvas->GetName()));
  r->floatParsFinal().Print("s");
  const RooArgList &fitParams = r->floatParsFinal();

  // Open new ROOT file save save result
  TString fOut;
  if (With_DY)
    fOut = Form("Fitresult_%s_M_%0.0f_%0.0f_withDY.root", Generator.Data(), Low_Mass, High_Mass);
  else
    fOut = Form("Fitresult_%s_M_%0.0f_%0.0f_noDY.root", Generator.Data(), Low_Mass, High_Mass);
  TFile f(fOut.Data(), "RECREATE");
  r->Write("rf607");
  f.Close();

  const Int_t n_DiMuSelection = 3;
  TString name_DiMuSelection[n_DiMuSelection];
  name_DiMuSelection[0].Form("Charm");
  name_DiMuSelection[1].Form("Beauty");
  name_DiMuSelection[2].Form("Mixed");
  RooRealVar *B_DimuMass[n_DiMuSelection];
  RooRealVar *n1_DimuMass[n_DiMuSelection];
  RooRealVar *n2_DimuMass[n_DiMuSelection];
  RooRealVar *B_DimuPt[n_DiMuSelection];
  RooRealVar *n1_DimuPt[n_DiMuSelection];
  RooRealVar *n2_DimuPt[n_DiMuSelection];

  // for (Int_t i = 0; i < n_DiMuSelection - 1; i++)
  // {
  //   B_DimuMass[i] = w->var(Form("B_DimuMassFrom%s", name_DiMuSelection[i].Data()));
  //   B_DimuMass[i]->setConstant(kTRUE);
  //   n1_DimuMass[i] = w->var(Form("n1_DimuMassFrom%s", name_DiMuSelection[i].Data()));
  //   n1_DimuMass[i]->setConstant(kTRUE);
  //   n2_DimuMass[i] = w->var(Form("n2_DimuMassFrom%s", name_DiMuSelection[i].Data()));
  //   n2_DimuMass[i]->setConstant(kTRUE);

  //   B_DimuPt[i] = w->var(Form("B_DimuPtFrom%s", name_DiMuSelection[i].Data()));
  //   B_DimuPt[i]->setConstant(kTRUE);
  //   n1_DimuPt[i] = w->var(Form("n1_DimuPtFrom%s", name_DiMuSelection[i].Data()));
  //   n1_DimuPt[i]->setConstant(kTRUE);
  //   n2_DimuPt[i] = w->var(Form("n2_DimuPtFrom%s", name_DiMuSelection[i].Data()));
  //   n2_DimuPt[i]->setConstant(kTRUE);

  //   printf("%s\n", name_DiMuSelection[i].Data());
  //   printf("B_DimuPt: %0.3f\n n1_DimuPt: %0.3f\n n2_DimuPt: %0.3f\n", B_DimuPt[i]->getError(), n1_DimuPt[i]->getError(), n2_DimuPt[i]->getError());

  //   printf("B_DimuMass: %0.3f\n n1_DimuMass: %0.3f\n n2_DimuMass: %0.3f\n", B_DimuMass[i]->getError(), n1_DimuMass[i]->getError(), n2_DimuMass[i]->getError());
  // }

  // printf("B_DimuMassFromCharm: %0.3e|| n1_DimuMassFromCharm: %0.3e|| n2_DimuMassFromCharm: %0.3e\n", B_DimuMassFromCharm->getError(), n1_DimuMassFromCharm->getError(), n2_DimuMassFromCharm->getError());

  // printf("B_DimuPtFromCharm: %0.3e|| n1_DimuPtFromCharm: %0.3e|| n2_DimuPtFromCharm: %0.3e\n", B_DimuPtFromCharm->getError(), n1_DimuPtFromCharm->getError(), n2_DimuPtFromCharm->getError());

  // printf("B_DimuMassFromBeauty: %0.3e|| n1_DimuMassFromBeauty: %0.3e|| n2_DimuMassFromBeauty: %0.3e\n", B_DimuMassFromBeauty->getError(), n1_DimuMassFromBeauty->getError(), n2_DimuMassFromBeauty->getError());

  // printf("B_DimuPtFromBeauty: %0.3e|| n1_DimuPtFromBeauty: %0.3e|| n2_DimuPtFromBeauty: %0.3e\n", B_DimuPtFromBeauty->getError(), n1_DimuPtFromBeauty->getError(), n2_DimuPtFromBeauty->getError());

  // printf("B_DimuPt: %0.3e|| n1_DimuPt: %0.3e|| n2_DimuPt: %0.3e\n", B_DimuMassFromMixed->getError(), n1_DimuMassFromMixed->getError(), n2_DimuMassFromMixed->getError());

  // printf("B_DimuPt: %0.3e|| n1_DimuPt: %0.3e|| n2_DimuPt: %0.3f\n", B_DimuPtFromMixed->getError(), n1_DimuPtFromMixed->getError(), n2_DimuPtFromMixed->getError());
  if (Plot_Likehood)
  {
    int font = 42;
    gStyle->SetImageScaling(3.);
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

    TCanvas *canvas = new TCanvas("c", "c", 1600, 1200);
    canvas->cd();
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.10);
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(0);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.15);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameBorderMode(0);
    // C o n s t r u c t   p l a i n   l i k e l i h o o d
    // ---------------------------------------------------
    // Construct unbinned likelihood
    RooAbsReal *nll_charm = simPdf.createNLL(*unbinned_combData_set, NumCPU(4));
    // Minimize likelihood w.r.t all parameters before making plots
    RooMinimizer(*nll_charm).lastMinuitFit();
    // Plot likelihood scan frac
    // RooWorkspace *w_Nll = new RooWorkspace("w_Nll", "workspace");
    // w_Nll->import(*nll) ;
    //
    //
    // w_Nll->writeToFile(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/root_files/likehood_%s.root", mass_range_data.Data()));
    TString title;
    if (mass_range_data.Contains("_LowMass_LowPt"))
      title.Form("4<m_{#mu#mu}<9 GeV/#it{c}^{2} && #it{p}_{T}<10 GeV/#it{c}");
    else if (mass_range_data.Contains("_withcut") & With_DY)
    {
      title.Form("excluding #Upsilon res. with DY");
    }
    else
    {
      title.Form("excluding #Upsilon res.");
    }

    RooPlot *frame_charm = normForC->frame(Bins(800), Range(1, 80000), Title(title));
    nll_charm->plotOn(frame_charm, Name("nll_charm"), ShiftToZero(), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
    // nll_charm->paramOn(*frame_charm);
    frame_charm->GetXaxis()->SetTitle("N_{#mu^{#plus}#mu^{#minus}} from fit");
    frame_charm->GetXaxis()->SetTitleOffset(1.5);
    frame_charm->GetYaxis()->SetTitleOffset(1.8);
    // frame_charm->GetYaxis()->SetTitle("Likehood");
    frame_charm->SetMaximum(52500);
    frame_charm->SetMinimum(10e-7);
    frame_charm->Draw();

    RooAbsReal *nll_beauty = simPdf.createNLL(*unbinned_combData_set, NumCPU(4), Name("nll_beauty"));
    // Minimize likelihood w.r.t all parameters before making plots
    RooMinimizer(*nll_beauty).lastMinuitFit();
    RooPlot *frame_beauty = normForB->frame(Bins(800), Range(1, 80000), Title("LL and profileLL in frac"));
    nll_beauty->plotOn(frame_beauty, Name("nll_beauty"), ShiftToZero(), LineStyle(kDashed), LineColor(kGreen + 2), LineWidth(5));
    frame_beauty->Draw("SAME");

    if (With_DY)
    {
      RooAbsReal *nll_DY = simPdf.createNLL(*unbinned_combData_set, NumCPU(4), Name("nll_DY"));
      // Minimize likelihood w.r.t all parameters before making plots
      RooMinimizer(*nll_DY).lastMinuitFit();
      RooPlot *frame_DY = normForDY->frame(Bins(800), Range(1, 80000), Title("LL and profileLL in frac"));
      nll_DY->plotOn(frame_DY, Name("nll_DY"), ShiftToZero(), LineStyle(kDashed), LineColor(kBlack), LineWidth(5));
      frame_DY->Draw("SAME");
    }

    TLegend *legend = new TLegend(0.475, 0.575, 0.75, 0.795);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0425);
    // legend->SetHeader("Data");
    legend->SetTextAlign(12);
    legend->AddEntry("nll_charm", "Likelihood #it{N}_{#mu^{#plus}#mu^{#minus} #leftarrow c}", "L");
    legend->AddEntry("nll_beauty", "Likelihood #it{N}_{#mu^{#plus}#mu^{#minus} #leftarrow b}", "L");
    if (With_DY)
      legend->AddEntry("nll_beauty", "Likelihood #it{N}_{#mu^{#plus}#mu^{#minus} #leftarrow DY}", "L");
    legend->Draw("SAME");
    if (With_DY)
    {

      canvas->SaveAs(Form("images/LLcanvas%s_withDY.pdf", mass_range_data.Data()));
      canvas->SaveAs(Form("images/LLcanvas%s_withDY.png", mass_range_data.Data()));
    }
    else
    {
      canvas->SaveAs(Form("images/LLcanvas%s_noDY.pdf", mass_range_data.Data()));
      canvas->SaveAs(Form("images/LLcanvas%s_noDY.png", mass_range_data.Data()));
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
  }
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
  pdf_legend->SetTextAlign(12);
  pdf_legend->SetHeader(Form("%s PDF", Generator.Data()));
  pdf_legend->AddEntry(roohist_name[2], "#mu^{#plus}#mu^{#minus} #leftarrow c,c", "L");
  pdf_legend->AddEntry(roohist_name[3], "#mu^{#plus}#mu^{#minus} #leftarrow b,b", "L");
  pdf_legend->AddEntry(roohist_name[4], "#mu^{#plus}#mu^{#minus} #leftarrow c,b", "L");
  if (With_DY)
    pdf_legend->AddEntry(roohist_name[5], "#mu^{#plus}#mu^{#minus} #leftarrow DY", "L");

  legend->Draw();
  fit_legend->Draw();
  pdf_legend->Draw();

  TLatex *letexTitle = new TLatex();
  letexTitle->SetNDC();
  letexTitle->SetTextFont(42);

  // letexTitle->DrawLatex(0.25, 0.26, "#mu^{#plus}#mu^{#minus} #leftarrow c,c");
  // letexTitle->DrawLatex(0.25, 0.19, "#mu^{#plus}#mu^{#minus} #leftarrow b,b");
  // letexTitle->DrawLatex(0.25, 0.12, "#mu^{#plus}#mu^{#minus} #leftarrow c,b");

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
    if (With_DY)
      letexTitle->DrawLatex(0.625, 0.525, Form("#it{N}_{#mu^{#plus}#mu^{#minus} #leftarrow DY} = %0.1f #pm %0.1f", fit_output[3]->getVal(), fit_output[3]->getError()));
    // letexTitle->DrawLatex(0.675, 0.625, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f", fit_output[3]->getVal()));
  }
  else if (choice == 1)
  {
    // letexTitle->DrawLatex(0.625, 0.825, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f %% #pm %0.3f %%", 100 * fit_output[0]->getVal(), 100 * fit_output[0]->getError()));
    // letexTitle->DrawLatex(0.625, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f %% #pm %0.3f %%", 100 * fit_output[1]->getVal(), 100 * fit_output[1]->getError()));
    // // letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.2e #pm %0.2e", fit_output[2]->getVal(), fit_output[2]->getError()));
    // letexTitle->DrawLatex(0.625, 0.625, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f %%", 100 - 100 * fit_output[0]->getVal() - 100 * fit_output[1]->getVal()));
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
  TLine *l = new TLine(minx, 1.0, max_x, 1.0);
  l->SetLineWidth(3);
  l->SetLineStyle(2);
  l->SetLineColor(kRed);

  TLine *l1 = new TLine(minx, 0.5, max_x, 0.5);
  l1->SetLineWidth(2);
  l1->SetLineStyle(9);
  l1->SetLineColor(kGray + 2);

  TLine *l2 = new TLine(minx, 1.5, max_x, 1.5);
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