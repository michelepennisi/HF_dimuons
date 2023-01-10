
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN
/// \author Michele Pennisi <michele.pennisi@cern.ch>, INFN

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

// Define Pt and mass value to sample MC distributions

Double_t Maxy = 150000;
Double_t Lowy = 0.0025;

// Charm=17157 (17157/48590)=0,353097345

// Beauty=29586 (29586/48590)=0,608890718

// Mixed=1847 (1847/48590)=0,038011937

// tot dimuon from HF in simulation=48590

// tot dimuon ULS from data=82643

//-------------------------------------------------------------
// Fit MC to extract pT and mass shapes
//-------------------------------------------------------------
TCanvas *printRooPlot_ratio(Int_t choice_Kin_Cut, RooPlot *frame, TString roohist_name, TF1 *pdf, TH1D *data, Double_t Low_Mass, Double_t High_Mass, Double_t Low_Pt, Double_t High_Pt);
TCanvas *printRooPlot(Int_t choice_Kin_Cut, RooPlot *frame, TString roohist_name[3], TString pdf_name[3], Int_t option);
TCanvas *printMC_ratio(RooPlot *frame, TH1 *data, TF1 *pdf, Color_t color, Int_t minx = 0, Int_t max_x = 30);

void extract_pdf(Int_t Low_Mass = 4, Int_t High_Mass = 30, Int_t Low_Pt = 0, Int_t High_Pt = 30)
{
  gROOT->ProcessLineSync(".x /home/michele_pennisi/dimuon_HF_pp/fit_data/fit_library/PtMassExpPdf.cxx+");
  gROOT->ProcessLineSync(".x /home/michele_pennisi/dimuon_HF_pp/fit_data/fit_library/PtMassPol1ExpPdf.cxx+");

  TFile *fIn = new TFile("/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Hist_MC/3_11_22/HF/Tree_HF_MCDimuHFTree_merged.root", "READ");

  const Int_t n_DiMuSelection = 3;
  TString name_DiMuSelection[n_DiMuSelection];
  name_DiMuSelection[0].Form("charm");
  name_DiMuSelection[1].Form("beauty");
  name_DiMuSelection[2].Form("mixed");

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
      Binning_m = 26;
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
    //  mass_range.Form("_withcut");
    Binning_m = 26;
    Binning_pt = 30;
  }

  TTree *pt_tree_MC[n_DiMuSelection];
  TTree *m_tree_MC[n_DiMuSelection];

  RooDataSet *M_Dimu[n_DiMuSelection];
  RooDataSet *Pt_Dimu[n_DiMuSelection];

  TH1 *m_histo[n_DiMuSelection];
  TH1 *pt_histo[n_DiMuSelection];

  RooPlot *frameDimuMass[n_DiMuSelection];
  RooPlot *frameDimuPt[n_DiMuSelection];

  RooRealVar *m = new RooRealVar("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", Low_Mass, High_Mass);
  // m->setRange("region1", 4, 9);
  // m->setRange("region2", 11, 30);
  m->setBins(Binning_m);
  RooRealVar *pt = new RooRealVar("pt", "#it{p}_{T} (GeV/#it{c})", Low_Pt, High_Pt);
  pt->setBins(Binning_pt);
  RooPlot *prova = pt->frame(Title(Form("prova")));

  Color_t color[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kAzure + 9};
  Color_t fillcolor[n_DiMuSelection] = {kMagenta - 10, kGreen - 10, kCyan - 10};

  for (Int_t i = 0; i < n_DiMuSelection; i++)
  {
    m_tree_MC[i] = (TTree *)fIn->Get(Form("rec_tree_mu%s%s", name_DiMuSelection[i].Data(), mass_range.Data()));
    printf("%s", Form("rec_tree_mu%s%s\n", name_DiMuSelection[i].Data(), mass_range.Data()));
    pt_tree_MC[i] = (TTree *)fIn->Get(Form("rec_tree_mu%s%s", name_DiMuSelection[i].Data(), mass_range.Data()));

    M_Dimu[i] = new RooDataSet(Form("M_Dimu_%s", name_DiMuSelection[i].Data()), Form("M_Dimu_%s", name_DiMuSelection[i].Data()), RooArgSet(*m), Import(*m_tree_MC[i]));

    Pt_Dimu[i] = new RooDataSet(Form("Pt_Dimu_%s", name_DiMuSelection[i].Data()), Form("Pt_Dimu_%s", name_DiMuSelection[i].Data()), RooArgSet(*pt), Import(*pt_tree_MC[i]));
    printf("%s", Form("rec_tree_mu%s%s\n", name_DiMuSelection[i].Data(), mass_range.Data()));
    m_histo[i] = M_Dimu[i]->createHistogram(Form("h_M_Dimu_%s", name_DiMuSelection[i].Data()), *m, Binning(Binning_m, Low_Mass, High_Mass));

    pt_histo[i] = Pt_Dimu[i]->createHistogram(Form("h_Pt_Dimu_%s", name_DiMuSelection[i].Data()), *pt, Binning(Binning_pt, Low_Pt, High_Pt));

    frameDimuMass[i] = m->frame(Title(Form("frameDimuMass_%s", name_DiMuSelection[i].Data())));

    frameDimuMass[i]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");

    // M_Dimu[i]->plotOn(frameDimuMass[i], Name(Form("M_%s", name_DiMuSelection[i].Data())),Normalization(0.00001), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i]), LineColor(color[i]), LineWidth(2), DrawOption("PEZ"));

    frameDimuPt[i] = pt->frame(Title(Form("frameDimuPt_%s", name_DiMuSelection[i].Data())));

    frameDimuPt[i]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

    // Pt_Dimu[i]->plotOn(frameDimuPt[i], Name(Form("Pt_%s", name_DiMuSelection[i].Data())),Normalization(0.00001), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i]), LineColor(color[i]), LineWidth(2), DrawOption("PEZ"));
  }

  RooWorkspace *w = new RooWorkspace("w", "workspace");
  w->factory(Form("PtMassExpPdf::pdfDimuMassFromcharm(m[%d,%d], B_DimuMassFromCharm[2.85,1,100], n1_DimuMassFromCharm[2.81,1,100], n2_DimuMassFromCharm[5,1,100])", Low_Mass, High_Mass));

  // w->factory(Form("PtMassPol1ExpPdf::pdfDimuPtFromcharm(pt[%d,%d], B_DimuPtFromCharm[2.85,1,100], n1_DimuPtFromCharm[2.81,1,100], n2_DimuPtFromCharm[2.43,1,100],C_DimuPtFromCharm[1.5,0,3], n3_DimuPtFromCharm[0.9,0.5,1.5])", Low_Pt, High_Pt));

  w->factory(Form("PtMassExpPdf::pdfDimuPtFromcharm(pt[%d,%d], B_DimuPtFromCharm[2.85,1,100], n1_DimuPtFromCharm[2.81,1,100], n2_DimuPtFromCharm[2.43,1,100])", Low_Pt, High_Pt));

  w->factory(Form("PtMassExpPdf::pdfDimuMassFrombeauty(m[%d,%d], B_DimuMassFromBeauty[4.16,1,100], n1_DimuMassFromBeauty[2.256,1,100], n2_DimuMassFromBeauty[2.83,1,100])", Low_Mass, High_Mass));

  // w->factory(Form("PtMassPol1ExpPdf::pdfDimuPtFrombeauty(pt[%d,%d], B_DimuPtFromBeauty[2.85,1,100], n1_DimuPtFromBeauty[2.81,1,100], n2_DimuPtFromBeauty[2.43,1,100],C_DimuPtFromBeauty[1.5,0,3], n3_DimuPtFromBeauty[0.9,0.5,1.5])", Low_Pt, High_Pt));

  w->factory(Form("PtMassExpPdf::pdfDimuPtFrombeauty(pt[%d,%d], B_DimuPtFromBeauty[2.85,1,100], n1_DimuPtFromBeauty[2.81,1,100], n2_DimuPtFromBeauty[2.43,1,100])", Low_Pt, High_Pt));

  w->factory(Form("PtMassExpPdf::pdfDimuMassFrommixed(m[%d,%d], B_DimuMassFromMixed[2.65,1,100], n1_DimuMassFromMixed[2.81,1,100], n2_DimuMassFromMixed[5,1,100])", Low_Mass, High_Mass));

  // w->factory(Form("PtMassPol1ExpPdf::pdfDimuPtFrommixed(pt[%d,%d], B_DimuPtFromMixed[2.85,1,100], n1_DimuPtFromMixed[2.81,1,100], n2_DimuPtFromMixed[2.43,1,100],C_DimuPtFromMixed[1.5,0,3], n3_DimuPtFromMixed[0.9,0.5,1.5])", Low_Pt, High_Pt));

  w->factory(Form("PtMassExpPdf::pdfDimuPtFrommixed(pt[%d,%d], B_DimuPtFromMixed[2.85,1,100], n1_DimuPtFromMixed[2.81,1,100], n2_DimuPtFromMixed[2.43,1,100])", Low_Pt, High_Pt));

  RooAbsPdf *pdfDimuPt[n_DiMuSelection];
  RooAbsPdf *pdfDimuM[n_DiMuSelection];
  TF1 *pt_Func[n_DiMuSelection];
  TF1 *m_Func[n_DiMuSelection];
  TCanvas *pt_canvas[n_DiMuSelection];
  TCanvas *m_canvas[n_DiMuSelection];

  for (Int_t i = 0; i < n_DiMuSelection; i++)
  {

    pdfDimuPt[i] = w->pdf(Form("pdfDimuPtFrom%s", name_DiMuSelection[i].Data()));
    auto result1 = pdfDimuPt[i]->fitTo(*Pt_Dimu[i], Minimizer("Minuit2"), Save(), SumW2Error(true));
    // pdfDimuPt[i]->fitTo(*Pt_Dimu[i]);
    Pt_Dimu[i]->plotOn(frameDimuPt[i], Name(Form("Pt_%s", name_DiMuSelection[i].Data())), Normalization(0.00001), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i]), LineColor(color[i]), LineWidth(2), DrawOption("PEZ"));

    pdfDimuPt[i]->plotOn(frameDimuPt[i], Name(Form("pdfDimuPtFrom%s", name_DiMuSelection[i].Data())), LineStyle(kSolid), LineColor(color[i]));
    // pdfDimuPt[i]->removeStringAttribute("fitrange");

    pdfDimuM[i] = w->pdf(Form("pdfDimuMassFrom%s", name_DiMuSelection[i].Data()));

    auto result2 = pdfDimuM[i]->fitTo(*M_Dimu[i], Minimizer("Minuit2"), Save(), SumW2Error(true));

    M_Dimu[i]->plotOn(frameDimuMass[i], Name(Form("M_%s", name_DiMuSelection[i].Data())), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i]), LineColor(color[i]), LineWidth(2), DrawOption("PEZ"));
    pdfDimuM[i]->plotOn(frameDimuMass[i], Name(Form("pdfDimuMFrom%s", name_DiMuSelection[i].Data())), LineStyle(kSolid), LineColor(color[i]));

    RooArgSet *pt_model = static_cast<RooArgSet *>(RooArgSet(*pdfDimuPt[i]).snapshot(true));                             // True means copy the PDF and everything it depends on
    auto &pt_modelcopied = static_cast<RooAbsPdf &>((*pt_model)[Form("pdfDimuPtFrom%s", name_DiMuSelection[i].Data())]); // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    RooArgSet *pt_modelobs = pt_modelcopied.getObservables(*Pt_Dimu[i]);
    RooArgSet *pt_modelPars = pt_modelcopied.getParameters(*pt_modelobs);
    pt_Func[i] = pt_modelcopied.asTF(*pt_modelobs, *pt_modelPars, *pt);

    RooArgSet *m_model = static_cast<RooArgSet *>(RooArgSet(*pdfDimuM[i]).snapshot(true));                               // True means copy the PDF and everything it depends on
    auto &m_modelcopied = static_cast<RooAbsPdf &>((*m_model)[Form("pdfDimuMassFrom%s", name_DiMuSelection[i].Data())]); // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    RooArgSet *m_modelobs = m_modelcopied.getObservables(*M_Dimu[i]);
    RooArgSet *m_modelPars = m_modelcopied.getParameters(*m_modelobs);
    m_Func[i] = m_modelcopied.asTF(*m_modelobs, *m_modelPars, *m);

    m_canvas[i] = printMC_ratio(frameDimuMass[i], m_histo[i], m_Func[i], color[i], Low_Mass, High_Mass);
    m_canvas[i]->SetTitle(Form("m_canvas_%s_unbinned", name_DiMuSelection[i].Data()));
    m_canvas[i]->SetName(Form("m_canvas_%s_unbinned", name_DiMuSelection[i].Data()));
    m_canvas[i]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/m_canvas_%s_unbinned%s.pdf", name_DiMuSelection[i].Data(), mass_range.Data()));

    pt_canvas[i] = printMC_ratio(frameDimuPt[i], pt_histo[i], pt_Func[i], color[i], Low_Pt, High_Pt);
    pt_canvas[i]->SetTitle(Form("pt_canvas_%s_unbinned", name_DiMuSelection[i].Data()));
    pt_canvas[i]->SetName(Form("pt_canvas_%s_unbinned", name_DiMuSelection[i].Data()));
    pt_canvas[i]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/pt_canvas_%s_unbinned%s.pdf", name_DiMuSelection[i].Data(), mass_range.Data()));
  }
  w->writeToFile(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/root_files/pdfMC_unbinned%s.root", mass_range.Data()));

  w->Print();
  gDirectory->Add(w);


}
//---------------------------------------------------------------------------------------//

TCanvas *printMC_ratio(RooPlot *frame, TH1 *data, TF1 *pdf, Color_t color, Int_t minx = 0, Int_t max_x = 30)
{
  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 1000);
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
  TString str;
  str.Form("%s", data->GetTitle());

  if (str.Contains("h_M"))
  {
    data->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#plus}} (GeV/#it{c}^{2})");
    data->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#plus}} (GeV/#it{c}^{2})^{-1}");
  }
  else
  {
    data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    data->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
  }

  TH2D *h_grid = new TH2D("h_grid", "", 100, minx, max_x, 100, Lowy, Maxy);
  h_grid->SetTitle("");

  h_grid->GetXaxis()->SetTitleOffset(1.3);
  h_grid->GetXaxis()->SetTitleSize(0.0475);
  h_grid->GetXaxis()->SetLabelSize(0.045);

  h_grid->GetYaxis()->SetNdivisions(505);
  h_grid->GetYaxis()->SetTitleOffset(0.9);
  h_grid->GetYaxis()->SetTitleSize(0.065);
  h_grid->GetYaxis()->SetLabelSize(0.055);

  h_grid->Draw();
  frame->Draw("SAME");

  // data->SetMarkerSize(1.75);
  // data->SetMarkerStyle(24);
  // data->SetMarkerColor(color);
  // data->SetLineColor(color);
  // data->SetLineWidth(2);
  data->Scale(1. / data->Integral(), "width");
  // data->Draw("PESAME");

  // pdf->SetLineColor(color);
  // pdf->SetLineWidth(2);

  // pdf->Draw("LSAME");

  TLegend *legend = new TLegend(0.315, 0.15, 0.55, 0.4);
  // legend->SetNColumns(2);
  // legend->AddEntry((TObject*)0, "", "");
  legend->SetFillStyle(0);
  legend->SetLineColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.0525);
  legend->SetHeader("MC");
  legend->SetTextAlign(11);
  legend->AddEntry(data, " ", "PE");

  TLegend *pdf_legend = new TLegend(0.375, 0.15, 0.635, 0.4);
  pdf_legend->SetTextAlign(11);
  pdf_legend->SetFillStyle(0);
  pdf_legend->SetBorderSize(0);
  pdf_legend->SetTextSize(0.0525);
  pdf_legend->SetHeader("PDF");
  pdf_legend->AddEntry(pdf, " ", "L");

  legend->Draw();
  pdf_legend->Draw();

  TLatex *letexTitle = new TLatex();
  letexTitle->SetNDC();
  letexTitle->SetTextFont(42);
  letexTitle->SetTextSize(0.05);
  if (str.Contains("charm"))
    letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow c,c");
  if (str.Contains("beauty"))
    letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow b,b");
  if (str.Contains("mixed"))
    letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow b,c");
  // if (roohist_name.Contains("Charm")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow c,c");
  // if (roohist_name[option].Contains("Beauty")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,b");
  // if (roohist_name[option].Contains("Mixed")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,c");

  letexTitle->SetTextSize(0.06);
  // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
  // printf("WOW %s\n",roohist_name[option].Data() );
  // letexTitle->SetTextSize(0.055);
  // letexTitle->DrawLatex(0.675, 0.825, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f #pm %0.3f", fit_output[0]->getVal(), fit_output[0]->getError()));
  // letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f #pm %0.3f", fit_output[1]->getVal(), fit_output[1]->getError()));
  // letexTitle->DrawLatex(0.675, 0.625, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f", 1 - fit_output[0]->getVal() - fit_output[1]->getVal()));

  letexTitle->SetTextSize(0.055);
  letexTitle->DrawLatex(0.355, 0.88, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
  letexTitle->DrawLatex(0.355, 0.81, "PYTHIA8 Monash Tune, N_{ev} = 2 #upoint 10^{8}");
  if (str.Contains("h_Pt"))
  {
    letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{#eta}_{#mu} < 4.0");
    letexTitle->DrawLatex(0.355, 0.66, Form("2.5 < #it{y}_{#mu^{#plus}#mu^{#minus}} < 4.0, 4 < #it{m}_{#mu^{#plus}#mu^{#minus}} < 30 GeV/#it{c}^{2}"));
  }
  else
  {
    letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{#eta}_{#mu} < 4.0");
    letexTitle->DrawLatex(0.355, 0.66, Form("2.5 < #it{y}_{#mu^{#plus}#mu^{#minus}} < 4.0"));
  }

  pad2->cd();
  pad2->SetTicks();
  TLine *l = new TLine(minx, 1.0, 30.0, 1.0);
  l->SetLineWidth(3);
  l->SetLineStyle(2);
  l->SetLineColor(kRed);

  TLine *l1 = new TLine(minx, 1.5, 30.0, 1.5);
  l1->SetLineWidth(2);
  l1->SetLineStyle(9);
  l1->SetLineColor(kGray + 2);

  TLine *l2 = new TLine(minx, 0.5, 30.0, 0.5);
  l2->SetLineWidth(2);
  l2->SetLineStyle(9);
  l2->SetLineColor(kGray + 2);

  TH2D *h_grid_ratio = new TH2D("h_grid_ratio", "", 100, minx, max_x, 100, -0.55, 2.55);
  h_grid_ratio->SetTitle("");

  TH1D *c_data = (TH1D *)data->Clone("c_data");

  h_grid_ratio->GetYaxis()->SetTitle(Form("#frac{MC}{PDF}"));
  h_grid_ratio->GetYaxis()->CenterTitle();
  h_grid_ratio->GetYaxis()->SetNdivisions(504);
  h_grid_ratio->GetYaxis()->SetTitleSize(0.08);
  // h_grid_ratio->GetYaxis()->SetTitleOffset(0.8);
  h_grid_ratio->GetYaxis()->SetLabelOffset(0.02);
  h_grid_ratio->GetYaxis()->SetLabelSize(0.1);

  h_grid_ratio->GetXaxis()->SetTitleSize(0.1);
  h_grid_ratio->GetXaxis()->SetTitleOffset(1.1);
  h_grid_ratio->GetXaxis()->SetLabelSize(0.1);
  h_grid_ratio->GetXaxis()->SetTitle(c_data->GetXaxis()->GetTitle());

  c_data->SetLineColor(kBlack);
  c_data->SetMarkerColor(kBlack);
  c_data->SetMarkerStyle(20);
  // c_data->Rebin(10);
  // c_data->Scale(1. / c_data->Integral(), "width");
  c_data->Divide(pdf);

  h_grid_ratio->Draw();
  c_data->Draw("PESAME");
  l->Draw();
  l1->Draw();
  l2->Draw();
  return canvas;
}

TCanvas *printRooPlot_ratio(Int_t choice_Kin_Cut, RooPlot *frame, TString roohist_name, TF1 *pdf, TH1D *data, Double_t Low_Mass, Double_t High_Mass, Double_t Low_Pt, Double_t High_Pt)
{
  gStyle->SetOptStat(0);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 1000);
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
  frame->SetTitle("");
  frame->SetMaximum(Maxy);
  frame->SetMinimum(Lowy);

  frame->GetXaxis()->SetTitleOffset(1.3);
  frame->GetXaxis()->SetTitleSize(0.0475);
  frame->GetXaxis()->SetLabelSize(0.045);

  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleOffset(0.9);
  frame->GetYaxis()->SetTitleSize(0.065);
  frame->GetYaxis()->SetLabelSize(0.055);

  frame->Draw();

  TLegend *legend = new TLegend(0.315, 0.15, 0.55, 0.4);
  // legend->SetNColumns(2);
  // legend->AddEntry((TObject*)0, "", "");
  legend->SetFillStyle(0);
  legend->SetLineColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.0525);
  legend->SetHeader("MC");
  legend->SetTextAlign(11);
  legend->AddEntry(roohist_name, " ", "PE");

  TLegend *pdf_legend = new TLegend(0.375, 0.15, 0.635, 0.4);
  pdf_legend->SetTextAlign(11);
  pdf_legend->SetFillStyle(0);
  pdf_legend->SetBorderSize(0);
  pdf_legend->SetTextSize(0.0525);
  pdf_legend->SetHeader("PDF");
  pdf_legend->AddEntry(roohist_name, " ", "L");

  legend->Draw();
  pdf_legend->Draw();

  TLatex *letexTitle = new TLatex();
  letexTitle->SetNDC();
  letexTitle->SetTextFont(42);
  letexTitle->SetTextSize(0.05);
  if (roohist_name.Contains("Charm"))
    letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow c,c");
  if (roohist_name.Contains("Beauty"))
    letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow b,b");
  if (roohist_name.Contains("Mixed"))
    letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow b,c");
  // if (roohist_name.Contains("Charm")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow c,c");
  // if (roohist_name[option].Contains("Beauty")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,b");
  // if (roohist_name[option].Contains("Mixed")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,c");

  letexTitle->SetTextSize(0.06);
  // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
  // printf("WOW %s\n",roohist_name[option].Data() );
  // letexTitle->SetTextSize(0.055);
  // letexTitle->DrawLatex(0.675, 0.825, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f #pm %0.3f", fit_output[0]->getVal(), fit_output[0]->getError()));
  // letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f #pm %0.3f", fit_output[1]->getVal(), fit_output[1]->getError()));
  // letexTitle->DrawLatex(0.675, 0.625, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f", 1 - fit_output[0]->getVal() - fit_output[1]->getVal()));
  Double_t minx = 0;
  Double_t maxx = 30.0;
  letexTitle->SetTextSize(0.055);
  letexTitle->DrawLatex(0.355, 0.88, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
  letexTitle->DrawLatex(0.355, 0.81, "PYTHIA8 Monash Tune, N_{ev} = 2 #upoint 10^{8}");
  if (roohist_name.BeginsWith("rooHistDimuPt"))
  {
    minx = Low_Pt;
    maxx = High_Pt;
    if (choice_Kin_Cut == 0)
    {
      letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{#eta}_{#mu} < 4.0");
      letexTitle->DrawLatex(0.355, 0.66, Form("2.5 < #it{y}_{#mu^{#plus}#mu^{#minus}} < 4.0, %0.0f < #it{m}_{#mu^{#plus}#mu^{#minus}} < %0.0f GeV/#it{c}^{2}", Low_Mass, High_Mass));
    }
    else if (choice_Kin_Cut == 1)
    {
      letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu} < 4.0");
      letexTitle->DrawLatex(0.355, 0.66, Form("%0.0f < #it{m}_{#mu^{#plus}#mu^{#minus}} < %0.0f GeV/#it{c}^{2}", Low_Mass, High_Mass));
    }
  }
  else
  {
    minx = Low_Mass;
    maxx = High_Mass;
    if (choice_Kin_Cut == 0)
    {

      letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu^{#plus}#mu^{#minus}} < 4.0");
      letexTitle->DrawLatex(0.355, 0.67, "2.5 < #it{#eta}_{#mu} < 4.0");
    }
    else if (choice_Kin_Cut == 1)
    {
      letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu} < 4.0");
    }
  }

  pad2->cd();
  pad2->SetTicks();
  TLine *l = new TLine(minx, 1.0, maxx, 1.0);
  l->SetLineWidth(3);
  l->SetLineStyle(2);
  l->SetLineColor(kRed);

  TLine *l1 = new TLine(minx, 1.5, maxx, 1.5);
  l1->SetLineWidth(2);
  l1->SetLineStyle(9);
  l1->SetLineColor(kGray + 2);

  TLine *l2 = new TLine(minx, 0.5, maxx, 0.5);
  l2->SetLineWidth(2);
  l2->SetLineStyle(9);
  l2->SetLineColor(kGray + 2);

  TH2D *h_grid_ratio = new TH2D("h_grid", "", 100, minx, maxx, 100, -0.5, 4.5);
  h_grid_ratio->SetTitle("");

  TH1D *c_data = (TH1D *)data->Clone("c_data");

  h_grid_ratio->GetYaxis()->SetTitle(Form("#frac{MC}{PDF}"));
  h_grid_ratio->GetYaxis()->CenterTitle();
  h_grid_ratio->GetYaxis()->SetNdivisions(504);
  h_grid_ratio->GetYaxis()->SetTitleSize(0.08);
  // h_grid_ratio->GetYaxis()->SetTitleOffset(0.8);
  h_grid_ratio->GetYaxis()->SetLabelOffset(0.02);
  h_grid_ratio->GetYaxis()->SetLabelSize(0.1);

  h_grid_ratio->GetXaxis()->SetTitleSize(0.1);
  h_grid_ratio->GetXaxis()->SetTitleOffset(1.1);
  h_grid_ratio->GetXaxis()->SetLabelSize(0.1);
  h_grid_ratio->GetXaxis()->SetTitle(c_data->GetXaxis()->GetTitle());

  c_data->SetLineColor(kBlack);
  c_data->SetMarkerColor(kBlack);
  c_data->SetMarkerStyle(20);
  // c_data->Rebin(10);
  c_data->Scale(1. / c_data->Integral(), "width");
  c_data->Divide(pdf);

  h_grid_ratio->Draw();
  c_data->Draw("PESAME");
  l->Draw();
  l1->Draw();
  l2->Draw();
  return canvas;
}

TCanvas *printRooPlot(Int_t choice_Kin_Cut, RooPlot *frame, TString roohist_name[3], TString pdf_name[3], Int_t option)
{
  TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 1000);
  canvas->SetTicks();
  canvas->cd();
  gPad->SetLogy(1);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.165);
  gPad->SetBottomMargin(0.095);

  frame->GetXaxis()->SetTitleOffset(1.025);
  frame->GetXaxis()->SetTitleSize(0.0425);
  frame->GetXaxis()->SetLabelSize(0.0375);

  frame->GetYaxis()->SetTitleOffset(1.65);
  frame->GetYaxis()->SetTitleSize(0.0425);
  frame->GetYaxis()->SetLabelSize(0.0375);

  frame->Draw();

  TLegend *legend = new TLegend(0.33, 0.15, 0.605, 0.4);
  // legend->SetNColumns(2);

  // legend->AddEntry((TObject*)0, "", "");

  legend->SetFillStyle(0);
  legend->SetLineColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.0375);

  TLegend *fit_legend = new TLegend(0.405, 0.15, 0.555, 0.4);
  fit_legend->SetTextAlign(11);
  fit_legend->SetFillStyle(0);
  fit_legend->SetBorderSize(0);
  fit_legend->SetTextSize(0.0375);

  if (option == 3)
  {
    legend->SetHeader("MC");
    legend->SetTextAlign(11);
    legend->AddEntry(roohist_name[0], " ", "PE");
    legend->AddEntry(roohist_name[1], " ", "PE");
    legend->AddEntry(roohist_name[2], " ", "PE");

    fit_legend->SetHeader("Fit");
    fit_legend->SetTextAlign(11);
    fit_legend->AddEntry(pdf_name[0], " ", "L");
    fit_legend->AddEntry(pdf_name[1], " ", "L");
    fit_legend->AddEntry(pdf_name[2], " ", "L");
  }
  else if (option < 3)
  {
    legend->SetHeader("MC");
    legend->SetTextAlign(11);
    legend->AddEntry(roohist_name[option], " ", "PE");

    fit_legend->SetHeader("Fit");
    fit_legend->SetTextAlign(11);
    fit_legend->AddEntry(pdf_name[option], " ", "L");
  }

  legend->Draw();
  fit_legend->Draw();

  TLatex *letexTitle = new TLatex();
  letexTitle->SetTextSize(0.045);
  letexTitle->SetNDC();
  letexTitle->SetTextFont(42);

  if (option == 3)
  {
    letexTitle->SetTextSize(0.03);
    letexTitle->DrawLatex(0.185, 0.292, "#mu^{#plus}#mu^{#minus} #leftarrow c,c");
    letexTitle->DrawLatex(0.185, 0.23, "#mu^{#plus}#mu^{#minus} #leftarrow b,b");
    letexTitle->DrawLatex(0.185, 0.175, "#mu^{#plus}#mu^{#minus} #leftarrow b,c");
    letexTitle->SetTextSize(0.0375);
    letexTitle->DrawLatex(0.35, 0.88, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.35, 0.82, "PYTHIA8 Monash Tune, N_{ev} = 2 #upoint 10^{8}");
    if (roohist_name[0].BeginsWith("pdfDimuPt"))
    {
      letexTitle->DrawLatex(0.35, 0.76, "Reconstructed #mu^{#plus}#mu^{#minus}, #it{m}_{#mu^{#plus}#mu^{#minus}} > 4 GeV/#it{c}^{2}");
      letexTitle->DrawLatex(0.7, 0.69, "2.5 < #it{y}_{#mu} < 4.0");
    }
    else
    {
      letexTitle->DrawLatex(0.35, 0.76, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu} < 4.0");
    }
  }
  else if (option < 3)
  {
    letexTitle->SetTextSize(0.0375);
    if (roohist_name[option].Contains("Charm"))
      letexTitle->DrawLatex(0.21, 0.292, "#mu^{#plus}#mu^{#minus} #leftarrow c,c");
    if (roohist_name[option].Contains("Beauty"))
      letexTitle->DrawLatex(0.21, 0.292, "#mu^{#plus}#mu^{#minus} #leftarrow b,b");
    if (roohist_name[option].Contains("Mixed"))
      letexTitle->DrawLatex(0.21, 0.292, "#mu^{#plus}#mu^{#minus} #leftarrow c,b");

    letexTitle->SetTextSize(0.045);
    letexTitle->DrawLatex(0.39, 0.88, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.39, 0.81, "PYTHIA8 Monash Tune, N_{ev} = 2 #upoint 10^{8}");
    letexTitle->DrawLatex(0.39, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu} < 4.0");
  }

  return canvas;
}
