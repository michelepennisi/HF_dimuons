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
using namespace RooFit;

Double_t n_charm_input = 27755.57;
Double_t n_beauty_input = 46720.55;
Double_t n_mixed_input = 2741.88;

Double_t n_total_input = 77218.;
//-------------------------------------------------------------
// Fit MC to extract pT and mass shapes
//-------------------------------------------------------------
//---------------------------------------------------------------------------------------//
TCanvas *printRooPlot(RooPlot *frame, TString pdf_name[3], RooRealVar *fit_output[3], Int_t option);
TCanvas *printRooPlot_ratio(RooPlot *frame, Bool_t norm, RooRealVar *fit_output[3], Int_t choice, TString roohist_name[5], TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x);
// 0 for fit via normalization, 1 for fit via fraction
void unbinned_test_toy_upsilon(Int_t choice = 0, Double_t Low_Mass = 4.0, Double_t High_Mass = 30.0, Double_t Low_Pt = 4.0, Double_t High_Pt = 30.0, Double_t N_Events = 100000, Double_t N_Upsilon = 5000)
{
   TString pdfMass_name[3];
   TString pdfPt_name[3];
   // TString pdfPt_name[3] = {"pdfptcharm", "pdfptbeauty", "pdfptmixed"};
   gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
   gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");
   TString mass_range;

   Double_t Binning_m;
   Double_t Binning_pt;
   if (Low_Mass == 4 && High_Mass == 9)
   {
      if (High_Pt == 10)
      {
         mass_range.Form("_LowMass_LowPt");
         Binning_m = 15;
         Binning_pt = 20;
      }
      else if (High_Pt == 30)
      {
         mass_range.Form("_LowMass");
         Binning_m = 20;
         Binning_pt = 30;
      }
   }

   else if (Low_Mass == 11 && High_Mass == 30)
   {
      mass_range.Form("_HighMass");
      Binning_m = 19;
      Binning_pt = 30;
   }
   else
   {
      mass_range.Form("");
      Binning_m = 26;
      Binning_pt = 30;
   }

   TFile *fIn = new TFile(Form("rooWorkspace_test_upsilon%s_Nev_%0.0f_NY_%0.0f.root", mass_range.Data(), N_Events, N_Upsilon));

   pdfMass_name[0].Form("pdfDimuMassFromcharm");
   pdfMass_name[1].Form("pdfDimuMassFrombeauty");
   pdfMass_name[2].Form("pdfDimuMassFrommixed");

   pdfPt_name[0].Form("pdfDimuPtFromcharm");
   pdfPt_name[1].Form("pdfDimuPtFrombeauty");
   pdfPt_name[2].Form("pdfDimuPtFrommixed");
   RooWorkspace *w = (RooWorkspace *)fIn->Get("w");
   w->Print();

   RooRealVar *m = w->var("m");
   RooRealVar *pt = w->var("pt");
   m->setBins(Binning_m);
   pt->setBins(Binning_pt);
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
   RooRealVar *C_DimuPtFromCharm = w->var("C_DimuPtFromCharm");
   C_DimuPtFromCharm->setConstant(kTRUE);
   RooRealVar *n3_DimuPtFromCharm = w->var("n3_DimuPtFromCharm");
   n3_DimuPtFromCharm->setConstant(kTRUE);

   RooRealVar *B_DimuPtFromBeauty = w->var("B_DimuPtFromBeauty");
   B_DimuPtFromBeauty->setConstant(kTRUE);
   RooRealVar *n1_DimuPtFromBeauty = w->var("n1_DimuPtFromBeauty");
   n1_DimuPtFromBeauty->setConstant(kTRUE);
   RooRealVar *n2_DimuPtFromBeauty = w->var("n2_DimuPtFromBeauty");
   n2_DimuPtFromBeauty->setConstant(kTRUE);
   RooRealVar *C_DimuPtFromBeauty = w->var("C_DimuPtFromBeauty");
   C_DimuPtFromBeauty->setConstant(kTRUE);
   RooRealVar *n3_DimuPtFromBeauty = w->var("n3_DimuPtFromBeauty");
   n3_DimuPtFromBeauty->setConstant(kTRUE);

   RooRealVar *B_DimuPtFromMixed = w->var("B_DimuPtFromMixed");
   B_DimuPtFromMixed->setConstant(kTRUE);
   RooRealVar *n1_DimuPtFromMixed = w->var("n1_DimuPtFromMixed");
   n1_DimuPtFromMixed->setConstant(kTRUE);
   RooRealVar *n2_DimuPtFromMixed = w->var("n2_DimuPtFromMixed");
   n2_DimuPtFromMixed->setConstant(kTRUE);
   RooRealVar *C_DimuPtFromMixed = w->var("C_DimuPtFromMixed");
   C_DimuPtFromMixed->setConstant(kTRUE);
   RooRealVar *n3_DimuPtFromMixed = w->var("n3_DimuPtFromMixed");
   n3_DimuPtFromMixed->setConstant(kTRUE);

   RooAbsPdf *pdfDimuMassFromCharm = w->pdf(Form("%s", pdfMass_name[0].Data()));
   RooAbsPdf *pdfDimuPtFromCharm = w->pdf(Form("%s", pdfPt_name[0].Data()));

   RooAbsPdf *pdfDimuMassFromBeauty = w->pdf(Form("%s", pdfMass_name[1].Data()));
   RooAbsPdf *pdfDimuPtFromBeauty = w->pdf(Form("%s", pdfPt_name[1].Data()));

   RooAbsPdf *pdfDimuMassFromMixed = w->pdf(Form("%s", pdfMass_name[2].Data()));
   RooAbsPdf *pdfDimuPtFromMixed = w->pdf(Form("%s", pdfPt_name[2].Data()));

   // RooRealVar *normForC;
   // RooRealVar *normForB;
   RooAddPdf *m_model;
   RooAddPdf *pt_model;

   if (choice == 0)
   {
      ////////////////////
      // Fit with normalization
      RooRealVar *normForC = new RooRealVar("n_charm_output", "number dimuon from c", 10000, 0, 2000000);
      RooRealVar *normForB = new RooRealVar("n_beauty_output", "number dimuon from b", 10000, 0, 2000000);

      // RooRealVar *normForMixed = new RooRealVar("n_mixed_input", "number dimuon from b,c", 3200);
      RooRealVar Roo_N_Events("N_events", " ", N_Events);
      RooRealVar Roo_N_Upsilon("N_Upsilon", " ", N_Upsilon);
      RooRealVar Roo_fr_mixed("Roo_fr_mixed", " ", 0.032);
      RooFormulaVar *normForMixed = new RooFormulaVar("n_mixed_input", "(@0-@1)*@2", RooArgList(Roo_N_Events, Roo_N_Upsilon, Roo_fr_mixed));
      // normForMixed->setConstant(kTRUE);

      m_model = new RooAddPdf("m_model", "dimuMassFromC + dimuMassFromB + dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
      pt_model = new RooAddPdf("pt_model", "dimuPtFromC + dimuPtFromB + dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
   }

   else if (choice == 1)
   {
      ////////////////////
      // Fit with fraction
      RooRealVar *normForC = new RooRealVar("fr_charm_output", "fraction dimuon from c", 0.2, 0., 1.);
      RooRealVar *normForB = new RooRealVar("fr_beauty_output", "fraction dimuon from b", 0.6, 0., 1.);
      RooFormulaVar *normForMixed = new RooFormulaVar("fr_mixed_output", "1-0.965", RooArgList(*normForC, *normForB));

      m_model = new RooAddPdf("m_model", "dimuMassFromC + dimuMassFromB + dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed), RooArgList(*normForB, *normForC));
      pt_model = new RooAddPdf("pt_model", "dimuPtFromC + dimuPtFromB + dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed), RooArgList(*normForB, *normForC));
   }

   else if (choice == 2)
   {
      ////////////////////
      // Fit with fraction
      RooRealVar *normForC = new RooRealVar("fr_charm_output", "fraction dimuon from c", 0.2, 0., 1.);
      RooRealVar *normForMixed = new RooRealVar("fr_mixed_output", "fraction dimuon from Mixed", 0.04, 0., 1.);
      RooFormulaVar *normForB = new RooFormulaVar("fr_beauty_output", "1-@0-@1", RooArgList(*normForC, *normForMixed));

      m_model = new RooAddPdf("m_model", "dimuMassFromC + dimuMassFromB + dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
      pt_model = new RooAddPdf("pt_model", "dimuPtFromC + dimuPtFromB + dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
   }
   else if (choice == 3)
   {
      ////////////////////
      // Fit with fraction
      RooRealVar *normForB = new RooRealVar("fr_beauty_output", "fraction dimuon from b", 0.6, 0., 1.);
      RooRealVar *normForMixed = new RooRealVar("fr_mixed_output", "fraction dimuon from Mixed", 0.04, 0., 1.);
      RooFormulaVar *normForC = new RooFormulaVar("fr_charm_output", "1-@0-@1", RooArgList(*normForB, *normForMixed));

      m_model = new RooAddPdf("m_model", "dimuMassFromC + dimuMassFromB + dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
      pt_model = new RooAddPdf("pt_model", "dimuPtFromC + dimuPtFromB + dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
   }

   RooDataSet *dataDimuMassFromCharmAndBeauty = (RooDataSet *)w->data("pdfDimuMassFromcharmData");
   RooDataSet *dataDimuPtFromCharmAndBeauty = (RooDataSet *)w->data("pdfDimuPtFromcharmData");

   TFile *fIn_upsilon = new TFile("/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/HistResults_merged.root", "READ");
   TH2D *h_ptm = (TH2D *)fIn_upsilon->Get("Dimuon/CMUL7/DQ_cut_match_LT_Yres_ULS/h_PtMdiMu_CMUL7_DQ_cut_match_LT_Yres_ULS");
   TH1D *h_pt = (TH1D *)h_ptm->ProjectionX();
   h_pt->Sumw2();
   h_pt->Rebin(10);
   h_pt->Scale(1., "width");

   RooWorkspace *q = new RooWorkspace("q", "q");
   q->factory(Form("PtMassExpPdf::pdf_PtUpisilon(pt[%d,%d], B_DimuPtUpisilon[0,0100], n1_DimuPtUpisilon[0,0,100], n2_DimuPtUpisilon[0,0,100])", 0, 25));

   RooAbsPdf *pdf_PtUpisilon = q->pdf("pdf_PtUpisilon");

   RooDataHist *Pt = new RooDataHist("Pt_Dimu", "Pt_Dimu", RooArgSet(*pt), Import(*h_pt));
   // Pt->Draw();
   auto pt_r = pdf_PtUpisilon->fitTo(*Pt, Minimizer("Minuit2"), Save(), SumW2Error(true));
   RooPlot *frame = pt->frame(Title("Imported TH1 with Poisson error bars"));
   Pt->plotOn(frame, LineColor(kMagenta + 2), DrawOption("PE"));
   pdf_PtUpisilon->plotOn(frame, LineColor(kMagenta + 2));
   frame->Draw();
   // Build gaussian pdf in terms of x,mean and sigma
   RooRealVar mean("mean", "mean of gaussian", 9.4);
   RooRealVar sigma("sigma", "width of gaussian", 0.12);

   RooGaussian gauss("gauss", "gaussian PDF", *m, mean, sigma);
   RooDataSet *m_upsilon = gauss.generate(*m, N_Upsilon);
   RooDataSet *pt_upsilon = pdf_PtUpisilon->generate(*pt, N_Upsilon);

   dataDimuMassFromCharmAndBeauty->append(*m_upsilon);
   dataDimuPtFromCharmAndBeauty->append(*pt_upsilon);

   // Define the combined dataset
   RooCategory sample("sample", "sample");
   sample.defineType("mass");
   sample.defineType("transversemomentum");
   RooDataSet *combData = new RooDataSet("combData", "combined data", RooArgSet(*m, *pt), Index(sample), Import("mass", *dataDimuMassFromCharmAndBeauty), Import("transversemomentum", *dataDimuPtFromCharmAndBeauty));

   RooFitResult *r_test = gauss.fitTo(*m_upsilon, Save());
   RooPlot *test = m->frame(Title("m_test"));
   m_upsilon->plotOn(test, Name("upsilon"), DrawOption("PEZ"));
   gauss.plotOn(test, Name("test"), LineStyle(kSolid), LineColor(kRed));
   new TCanvas();
   test->Draw();
   // Define the pdf for the simultaneous fit
   RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
   simPdf.addPdf(*m_model, "mass");
   simPdf.addPdf(*pt_model, "transversemomentum");
   // simPdf.fitTo(*combData) ;
   // m->setRange("region1", 4, 9);
   // m->setRange("region2", 11, 30);
   // RooFitResult *r = simPdf.fitTo(*combData, Range("region1,region2"), Save());
   RooFitResult *r = simPdf.fitTo(*combData, Save());

   // printf("normForMixed %0.1f\n", normForMixed->GetValue());

   RooPlot *m_frame = m->frame(Title("m_frame"));
   m_frame->SetTitle(" ");
   m_frame->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
   m_frame->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");

   combData->plotOn(m_frame, Name("combDatamass"), Cut("sample==sample::mass"), DrawOption("PEZ"));

   // simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *combData), Range("region1"), LineStyle(kSolid), LineColor(kRed));
   // simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *combData), Range("region2"), LineStyle(kSolid), LineColor(kRed));

   simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *combData), LineStyle(kSolid), LineColor(kRed));

   simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *combData), VisualizeError(*r), FillColor(kRed));

   simPdf.plotOn(m_frame, Name("pdfmasscharm"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[0].Data())), ProjWData(sample, *combData), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
   simPdf.plotOn(m_frame, Name("pdfmassbeauty"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[1].Data())), ProjWData(sample, *combData), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
   simPdf.plotOn(m_frame, Name("pdfmassmixed"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[2].Data())), ProjWData(sample, *combData), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

   // m_frame->Draw();

   RooPlot *pt_frame = pt->frame(Title("pt_frame"));
   pt_frame->SetTitle(" ");
   pt_frame->SetMaximum(60000);
   pt_frame->SetMinimum(2);
   pt_frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   pt_frame->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
   combData->plotOn(pt_frame, Name("combDatapt"), Cut("sample==sample::transversemomentum"), DrawOption("PEZ"));

   simPdf.plotOn(pt_frame, Name("pdfpt"), Slice(sample, "transversemomentum"), ProjWData(sample, *combData), LineStyle(kSolid), LineColor(kRed));
   simPdf.plotOn(pt_frame, Name("pdfptcharm"), Slice(sample, "transversemomentum"), Components(Form("%s", pdfPt_name[0].Data())), ProjWData(sample, *combData), LineStyle(kDashed), LineColor(kMagenta - 2), LineWidth(5));
   simPdf.plotOn(pt_frame, Name("pdfptbeauty"), Slice(sample, "transversemomentum"), Components(Form("%s", pdfPt_name[1].Data())), ProjWData(sample, *combData), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
   simPdf.plotOn(pt_frame, Name("pdfptmixed"), Slice(sample, "transversemomentum"), Components(Form("%s", pdfPt_name[2].Data())), ProjWData(sample, *combData), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

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
      fit_output[1] = (RooRealVar *)r->floatParsFinal().find("fr_beauty_output");
      fit_output[2] = 0;
   }
   else if (choice == 2)
   {
      fit_output[0] = (RooRealVar *)r->floatParsFinal().find("fr_charm_output");
      fit_output[1] = (RooRealVar *)r->floatParsFinal().find("fr_mixed_output");
      fit_output[2] = 0;
   }
   else if (choice == 3)
   {
      fit_output[0] = (RooRealVar *)r->floatParsFinal().find("fr_beauty_output");
      fit_output[1] = (RooRealVar *)r->floatParsFinal().find("fr_mixed_output");
      fit_output[2] = 0;
   }

   // TCanvas *c_fim = printRooPlot(pt_frame, pdfPt_name, fit_output, choice);
   // c_fipt->SetName("c_fipt");

   // TCanvas *c_fim = printRooPlot(m_frame, pdfMass_name, fit_output, choice);
   // c_fim->SetName("c_fipt");

   // RooPlot *m_frame2 = m->frame(Title("m_frame"));
   // RooHist *hMassresid_Mixed = m_frame->residHist("combDatamass", "pdfmass");
   // m_frame2->addPlotable(hMassresid_Mixed, "PE");

   // TCanvas *c_fitmass = new TCanvas();
   // c_fitmass->cd();
   // m_frame2->Draw();

   // // c_fitmass->SetName("c_fitmass");
   // // c_fitmass->SetTitle("c_fitmass");
   // //
   // // TCanvas *c_fipt=printRooPlot(pt_frame,pdfPt_name);
   // // c_fipt->SetName("c_fipt");
   // // c_fipt->SetTitle("c_fipt");

   // TH1 *pt_toy = combData->createHistogram("pt_toy", *pt, Binning(300, 0.0, 30.0));
   // TH1 *mass_toy = combData->createHistogram("mass_toy", *m, Binning(300, Low_Mass, High_Mass));

   // RooWorkspace *w_output = new RooWorkspace("w_output", "workspace");
   // w_output->import(*combData);

   // // w_output->import(simPdf.addPdf(m_model,"mass"));

   // w_output->writeToFile(Form("../data/rooWorkspaceOut%0.0f_%0.0f.root", Low_Mass, High_Mass));
   // w_output->Print();
   // gDirectory->Add(w_output);
   TString Pt_name[5] = {"combDatapt", "pdfpt", "pdfptcharm", "pdfptbeauty", "pdfptmixed"};
   TString Mass_name[5] = {"combDatamass", "pdfmass", "pdfmasscharm", "pdfmassbeauty", "pdfmassmixed"};

   RooArgSet *m_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*m_model).snapshot(true)); // True means copy the PDF and everything it depends on
   auto &m_modelcopiedPdf = static_cast<RooAbsPdf &>((*m_modelcopyOfEverything)["m_model"]);          // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
   RooArgSet *m_modelobs = m_modelcopiedPdf.getObservables(*dataDimuMassFromCharmAndBeauty);

   RooArgSet *m_modelPars = m_modelcopiedPdf.getParameters(*m_modelobs);

   TF1 *m_modelFunc = m_modelcopiedPdf.asTF(*m_modelobs, *m_modelPars, *m);

   RooArgSet *pt_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*pt_model).snapshot(true)); // True means copy the PDF and everything it depends on
   auto &pt_modelcopiedPdf = static_cast<RooAbsPdf &>((*pt_modelcopyOfEverything)["pt_model"]);         // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
   RooArgSet *pt_modelobs = pt_modelcopiedPdf.getObservables(*dataDimuPtFromCharmAndBeauty);
   RooArgSet *pt_modelPars = pt_modelcopiedPdf.getParameters(*pt_modelobs);
   TF1 *pt_modelFunc = pt_modelcopiedPdf.asTF(*pt_modelobs, *pt_modelPars, *pt);

   TH1 *hDimuPt_toy = dataDimuPtFromCharmAndBeauty->createHistogram("h_pttoy", *pt, Binning(Binning_pt, 0, High_Pt));

   hDimuPt_toy->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   hDimuPt_toy->GetYaxis()->SetTitle("d#it{N}/#it{p}_{T} (GeV/#it{c})^{-1}");

   TString info;

   info.Form("unbinned_data_optfit%d_M%0.0f_%0.0f", choice, Low_Mass, High_Mass);

   TCanvas *pt_canvas = printRooPlot_ratio(pt_frame, false, fit_output, choice, Pt_name, pt_modelFunc, hDimuPt_toy, 0, High_Pt);
   pt_canvas->SetName(Form("pt_upsilontest_%s_NEv_%0.0f_NY_%0.0f", info.Data(), N_Events, N_Upsilon));
   pt_canvas->SetTitle(Form("pt_upsilontest_%s_NEv_%0.0f_NY_%0.0f", info.Data(), N_Events, N_Upsilon));
   pt_canvas->SaveAs(Form("plot/%s.pdf", pt_canvas->GetName()));
   // m_frame->SetMaximum(hDimuM_data->GetMaximum() * 1000);
   // m_frame->SetMinimum(2e-14);
   TH1 *hDimuM_toy = dataDimuMassFromCharmAndBeauty->createHistogram("h_mdata", *m, Binning(Binning_m, Low_Mass, High_Mass));
   hDimuM_toy->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
   hDimuM_toy->GetYaxis()->SetTitle("d#it{N}/#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

   TCanvas *m_canvas = printRooPlot_ratio(m_frame, false, fit_output, choice, Mass_name, m_modelFunc, hDimuM_toy, Low_Mass, High_Mass);
   m_canvas->SetName(Form("m_upsilontest_%s_NEv_%0.0f_NY_%0.0f", info.Data(), N_Events, N_Upsilon));
   m_canvas->SetTitle(Form("m_upsilontest_%s_NEv_%0.0f_NY_%0.0f", info.Data(), N_Events, N_Upsilon));
   m_canvas->SaveAs(Form("plot/%s.pdf", m_canvas->GetName()));

   TFile *fOut = new TFile(Form("test_upsilon_rooOutput_unbinned%s_NEv_%0.0f_NY_%0.0f.root", mass_range.Data(), N_Events, N_Upsilon), "UPDATE");
   fOut->cd();
   r->Write("fit_result_opt_unbinned");

   fOut->WriteObject(m_frame, Form("m_frame_opt%d_%s", choice, mass_range.Data()));
   fOut->WriteObject(pt_frame, Form("pt_frame_opt%d_%s", choice, mass_range.Data()));

   // pt_toy->Write(0, 2, 0);
   // mass_toy->Write(0, 2, 0);

   fOut->Close();

   // printf("%s\n", Form(, suffix.Data()));
   // // TH2F* prova=combData->createHistogram(*pt,*m,100,100,"","prova");
   // // hh_data_all->Draw("PE");
   //
   // TH2* hh_sigmag2_frac = (TH2*) simPdf.createHistogram("pt,transversemomentum",100,100) ;
   // hh_sigmag2_frac->SetMarkerColor(kRed);
   // hh_sigmag2_frac->Draw("PE");
   // TH2* hh_sigmag2 = (TH2*) pdfDimuPtFromCharm->createHistogram("pt,pdfDimuPtFromCharm",100,100) ;
   // hh_sigmag2->Draw("PESAME");

   cout << "Mass chi^2 = " << m_frame->chiSquare() << endl;
   cout << "Pt chi^2 = " << pt_frame->chiSquare() << endl;
}

TCanvas *printRooPlot(RooPlot *frame, TString pdf_name[3], RooRealVar *fit_output[3], Int_t option)
{
   TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 800);
   canvas->SetTicks();
   canvas->cd();
   gPad->SetLogy();
   gPad->SetTopMargin(0.05);
   gPad->SetRightMargin(0.03);
   gPad->SetLeftMargin(0.14);
   gPad->SetBottomMargin(0.14);

   frame->GetXaxis()->SetTitleOffset(1.3);
   frame->GetXaxis()->SetTitleSize(0.0475);
   frame->GetXaxis()->SetLabelSize(0.045);

   frame->GetYaxis()->SetTitleOffset(1.3);
   frame->GetYaxis()->SetTitleSize(0.0475);
   frame->GetYaxis()->SetLabelSize(0.045);

   frame->Draw();

   TLegend *legend = new TLegend(0.53, 0.415, 1.0, 0.615);
   // legend->SetNColumns(2);
   // legend->AddEntry((TObject*)0, "", "");
   legend->SetFillStyle(0);
   legend->SetLineColor(kWhite);
   legend->SetBorderSize(0);
   legend->SetTextSize(0.0425);
   legend->SetHeader("Toy MC");
   legend->SetTextAlign(11);
   legend->AddEntry("combDatamass", " ", "LP");

   TLegend *fit_legend = new TLegend(0.655, 0.415, 0.815, 0.615);
   fit_legend->SetTextAlign(11);
   fit_legend->SetFillStyle(0);
   fit_legend->SetBorderSize(0);
   fit_legend->SetTextSize(0.0425);
   fit_legend->SetHeader("Fit");
   if (pdf_name[0].Contains("pt"))
      fit_legend->AddEntry("pdfpt", " ", "L");
   else if (pdf_name[0].Contains("mass"))
      fit_legend->AddEntry("pdfmass", " ", "L");

   TLegend *pdf_legend = new TLegend(0.175, 0.16, 0.475, 0.41);
   pdf_legend->SetTextAlign(11);
   pdf_legend->SetFillStyle(0);
   pdf_legend->SetBorderSize(0);
   pdf_legend->SetTextSize(0.0425);
   pdf_legend->SetHeader("PDF");
   pdf_legend->AddEntry(pdf_name[0], " ", "L");
   pdf_legend->AddEntry(pdf_name[1], " ", "L");
   pdf_legend->AddEntry(pdf_name[2], " ", "L");

   legend->Draw();
   fit_legend->Draw();
   pdf_legend->Draw();

   TLatex *letexTitle = new TLatex();
   letexTitle->SetTextSize(0.045);
   letexTitle->SetNDC();
   letexTitle->SetTextFont(42);

   letexTitle->SetTextSize(0.0375);
   letexTitle->DrawLatex(0.715, 0.455, "#mu^{#plus}#mu^{#minus} #leftarrow c+b");
   letexTitle->DrawLatex(0.25, 0.305, "#mu^{#plus}#mu^{#minus} #leftarrow c,c");
   letexTitle->DrawLatex(0.25, 0.245, "#mu^{#plus}#mu^{#minus} #leftarrow b,b");
   letexTitle->DrawLatex(0.25, 0.185, "#mu^{#plus}#mu^{#minus} #leftarrow c,b");
   letexTitle->SetTextSize(0.0425);
   letexTitle->DrawLatex(0.39, 0.88, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
   if (pdf_name[0].Contains("pt"))
      letexTitle->DrawLatex(0.39, 0.81, "Toy MC, N_{ev} = 82643, #it{m}_{#mu^{#plus}#mu^{#minus}} > 4 GeV/#it{c}^{2}");
   else if (pdf_name[0].Contains("mass"))
   {
      letexTitle->DrawLatex(0.39, 0.81, "Toy MC, N_{ev} = 82643");
   }

   if (option == 0)
   {
      letexTitle->DrawLatex(0.39, 0.74, Form("#it{N}^{Input}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.0f, #it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.0f #pm %0.0f", n_charm_input, fit_output[0]->getVal(), fit_output[0]->getError()));
      letexTitle->DrawLatex(0.39, 0.67, Form("#it{N}^{Input}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.0f, #it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.0f #pm %0.0f", n_beauty_input, fit_output[1]->getVal(), fit_output[1]->getError()));
      letexTitle->DrawLatex(0.39, 0.60, Form("#it{N}^{Input}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.0f, #it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.0f #pm %0.0f", n_mixed_input, fit_output[2]->getVal(), fit_output[2]->getError()));
   }
   else if (option == 1)
   {
      Double_t fr_th_charm = n_charm_input / n_total_input;
      Double_t fr_th_beaut = n_beauty_input / n_total_input;
      Double_t fr_th_mixed = n_mixed_input / n_total_input;

      letexTitle->DrawLatex(0.39, 0.74, Form("#it{N}^{Input}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.2f, #it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.2f #pm %0.2f", fr_th_charm, fit_output[0]->getVal(), fit_output[0]->getError()));
      (
          letexTitle->DrawLatex)(0.39, 0.67, Form("#it{N}^{Input}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.2f, #it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.2f #pm %0.2f", fr_th_beaut, fit_output[1]->getVal(), fit_output[1]->getError()));
      letexTitle->DrawLatex(0.39, 0.60, Form("#it{N}^{Input}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.2f, #it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.2f #pm %0.2f", fr_th_mixed, 1 - fit_output[0]->getVal() - fit_output[1]->getVal(), 0.0));
   }

   return canvas;
}

TCanvas *printRooPlot_ratio(RooPlot *frame, Bool_t norm, RooRealVar *fit_output[3], Int_t choice, TString roohist_name[5], TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x)
{
   gStyle->SetOptStat(0);
   TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1100);
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
   if (norm)
   {
      frame->SetMaximum(1.5e-0);
      frame->SetMinimum(1.5e-12);
   }
   else
   {
      frame->SetMaximum(1.5e+8);
      frame->SetMinimum(1.5e-2);
   }

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

   Double_t fr_th_charm = (n_charm_input / n_total_input)*100000;
   Double_t fr_th_beauty = (n_beauty_input / n_total_input)*100000;
   Double_t fr_th_mixed = (n_mixed_input / n_total_input)*100000;

   if (choice == 0)
   {
      letexTitle->DrawLatex(0.575, 0.825, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.1f #pm %0.1f", fit_output[0]->getVal(), fit_output[0]->getError()));
      letexTitle->DrawLatex(0.575, 0.725, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.1f #pm %0.1f", fit_output[1]->getVal(), fit_output[1]->getError()));

      letexTitle->DrawLatex(0.575, 0.625, Form("#it{N}^{Input}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.1f", fr_th_charm));
      letexTitle->DrawLatex(0.775, 0.625, Form("#it{N}^{Input}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.1f", fr_th_beauty));
      // letexTitle->DrawLatex(0.575, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.2e #pm %0.2e", fit_output[2]->getVal(), fit_output[2]->getError()));
      // letexTitle->DrawLatex(0.575, 0.625, Form("#it{N}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f", fit_output[3]->getVal()));
   }
   else if (choice == 1)
   {
      letexTitle->DrawLatex(0.575, 0.825, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f %% #pm %0.3f %%", 100 * fit_output[0]->getVal(), 100 * fit_output[0]->getError()));
      letexTitle->DrawLatex(0.575, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f %% #pm %0.3f %%", 100 * fit_output[1]->getVal(), 100 * fit_output[1]->getError()));
      // letexTitle->DrawLatex(0.575, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.2e #pm %0.2e", fit_output[2]->getVal(), fit_output[2]->getError()));
      letexTitle->DrawLatex(0.575, 0.625, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f %%", 100 - 100 * fit_output[0]->getVal() - 100 * fit_output[1]->getVal()));
   }

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

   TLine *l1 = new TLine(minx, 0.5, 30.0, 0.5);
   l1->SetLineWidth(2);
   l1->SetLineStyle(9);
   l1->SetLineColor(kGray + 2);

   TLine *l2 = new TLine(minx, 1.5, 30.0, 1.5);
   l2->SetLineWidth(2);
   l2->SetLineStyle(9);
   l2->SetLineColor(kGray + 2);

   TH2D *h_grid_ratio = new TH2D("h_grid", "", 100, minx, max_x, 100, -0.5, 4.25);
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