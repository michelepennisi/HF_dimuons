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

TCanvas *histo2(Int_t i, TH1D *hist_scaled_Low2Up, TH1D *hist_scaled_Up2Low, TH1D *hist_original, TF1 *pdf_scaled_Low2Up, TF1 *pdf_scaled_Up2Low, TF1 *pdf_original, TH1D *hint, Color_t color, Color_t fillcolor, Int_t minx = 0, Int_t max_x = 30);
double FuncPtMass(double *x, double *par);
void param_deviation();
void create_workspace();

TCanvas *printRooPlot_ratio(RooPlot *frame, Bool_t norm, RooRealVar *fit_output[3], Int_t choice, TString roohist_name[5], TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x);

const Int_t n_DiMuSelection = 3;
TString name_DiMuSelection[n_DiMuSelection];

void systematic_pdf()
{
    param_deviation();
    // create_workspace();
}

void create_workspace()
{
    Int_t n_options = 4;
    TString name_options[n_options];
    name_options[0].Form("unbinned");
    name_options[1].Form("original");
    name_options[2].Form("scaled_Low2Up");
    name_options[3].Form("scaled_Up2Low");

    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");
    name_DiMuSelection[0].Form("Charm");
    name_DiMuSelection[1].Form("Beauty");
    name_DiMuSelection[2].Form("Mixed");

    TFile *fIn_param = new TFile("signal_extraction_systematic_pdf.root", "READ");

    RooWorkspace *w[n_options];

    for (size_t p = 0; p < n_options; p++)
    {
        TH1D *Param_Pt[n_DiMuSelection];
        TH1D *Param_M[n_DiMuSelection];

        RooRealVar *m = new RooRealVar("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", 4, 30);
        // m->setBins(Binning_m);
        RooRealVar *pt = new RooRealVar("pt", "#it{p}_{T} (GeV/#it{c})", 0, 30);
        // pt->setBins(Binning_pt);

        RooRealVar *B_DimuMass[n_DiMuSelection];
        RooRealVar *n1_DimuMass[n_DiMuSelection];
        RooRealVar *n2_DimuMass[n_DiMuSelection];
        RooRealVar *B_DimuPt[n_DiMuSelection];
        RooRealVar *n1_DimuPt[n_DiMuSelection];
        RooRealVar *n2_DimuPt[n_DiMuSelection];
        for (Int_t i = 0; i < n_DiMuSelection; i++)
        {
            printf("%s", Form("Param_M_%s_from%s", name_options[p].Data(), name_DiMuSelection[i].Data()));
            Param_Pt[i] = (TH1D *)fIn_param->Get(Form("Param_Pt_%s_from%s", name_options[p].Data(), name_DiMuSelection[i].Data()));
            Param_M[i] = (TH1D *)fIn_param->Get(Form("Param_M_%s_from%s", name_options[p].Data(), name_DiMuSelection[i].Data()));
            Param_M[i]->Draw();

            B_DimuMass[i] = new RooRealVar(Form("B_DimuMass_from%s", name_DiMuSelection[i].Data()), Form("B_DimuMass_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_M[i]->GetBinContent(1));
            n1_DimuMass[i] = new RooRealVar(Form("n1_DimuMass_from%s", name_DiMuSelection[i].Data()), Form("n1_DimuMass_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_M[i]->GetBinContent(2));
            n2_DimuMass[i] = new RooRealVar(Form("n2_DimuMass_from%s", name_DiMuSelection[i].Data()), Form("n2_DimuMass_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_M[i]->GetBinContent(3));

            B_DimuMass[i]->setError(Param_M[i]->GetBinError(1));
            n1_DimuMass[i]->setError(Param_M[i]->GetBinError(2));
            n2_DimuMass[i]->setError(Param_M[i]->GetBinError(3));

            B_DimuPt[i] = new RooRealVar(Form("B_DimuPt_from%s", name_DiMuSelection[i].Data()), Form("B_DimuPt_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_Pt[i]->GetBinContent(1));
            n1_DimuPt[i] = new RooRealVar(Form("n1_DimuPt_from%s", name_DiMuSelection[i].Data()), Form("n1_DimuPt_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_Pt[i]->GetBinContent(2));
            n2_DimuPt[i] = new RooRealVar(Form("n2_DimuPt_from%s", name_DiMuSelection[i].Data()), Form("n2_DimuMass_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_Pt[i]->GetBinContent(3));

            B_DimuPt[i]->setError(Param_Pt[i]->GetBinError(1));
            n1_DimuPt[i]->setError(Param_Pt[i]->GetBinError(2));
            n2_DimuPt[i]->setError(Param_Pt[i]->GetBinError(3));

            printf("%s\n", name_DiMuSelection[i].Data());
            printf("B_DimuPt: %0.3f\n n1_DimuPt: %0.3f\n n2_DimuPt: %0.3f\n", B_DimuPt[i]->getVal(), n1_DimuPt[i]->getVal(), n2_DimuPt[i]->getVal());
        }

        // return;

        w[p] = new RooWorkspace(Form("w_%s", name_options[p].Data()), Form("w_%s", name_options[p].Data()));

        w[p]->factory(Form("PtMassExpPdf::pdfDimuMassFromcharm(m[4.0, 30.0], B_DimuMassFromCharm[%0.10f], n1_DimuMassFromCharm[%0.10f], n2_DimuMassFromCharm[%0.10f])", B_DimuMass[0]->getVal(), n1_DimuMass[0]->getVal(), n2_DimuMass[0]->getVal()));

        w[p]->factory(Form("PtMassExpPdf::pdfDimuPtFromcharm(pt[0.0, 30.0], B_DimuPtFromCharm[%0.10f], n1_DimuPtFromCharm[%0.10f], n2_DimuPtFromCharm[%0.10f])", B_DimuPt[0]->getVal(), n1_DimuPt[0]->getVal(), n2_DimuPt[0]->getVal()));

        w[p]->factory(Form("PtMassExpPdf::pdfDimuMassFrombeauty(m[4.0, 30.0], B_DimuMassFromBeauty[%0.10f], n1_DimuMassFromBeauty[%0.10f], n2_DimuMassFromBeauty[%0.10f])", B_DimuMass[1]->getVal(), n1_DimuMass[1]->getVal(), n2_DimuMass[1]->getVal()));

        w[p]->factory(Form("PtMassExpPdf::pdfDimuPtFrombeauty(pt[0.0, 30.0], B_DimuPtFromBeauty[%0.10f], n1_DimuPtFromBeauty[%0.10f], n2_DimuPtFromBeauty[%0.10f])", B_DimuPt[1]->getVal(), n1_DimuPt[1]->getVal(), n2_DimuPt[1]->getVal()));

        w[p]->factory(Form("PtMassExpPdf::pdfDimuMassFrommixed(m[4.0, 30.0], B_DimuMassFromMixed[%0.10f], n1_DimuMassFromMixed[%0.10f], n2_DimuMassFromMixed[%0.10f])", B_DimuMass[2]->getVal(), n1_DimuMass[2]->getVal(), n2_DimuMass[2]->getVal()));

        w[p]->factory(Form("PtMassExpPdf::pdfDimuPtFrommixed(pt[0.0, 30.0], B_DimuPtFromMixed[%0.10f], n1_DimuPtFromMixed[%0.10f], n2_DimuPtFromMixed[%0.10f])", B_DimuPt[2]->getVal(), n1_DimuPt[2]->getVal(), n2_DimuPt[2]->getVal()));

        w[p]->writeToFile("workspaces_systematic_pdf.root", kFALSE);

        w[p]->Print();
        gDirectory->Add(w[p]);
    }

    return;
}

void not_working()
{
    Int_t Binning_m = 26;
    Int_t Binning_pt = 30;

    name_DiMuSelection[0].Form("Charm");
    name_DiMuSelection[1].Form("Beauty");
    name_DiMuSelection[2].Form("Mixed");

    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");

    TString mass_range_data("_withcut");

    TString mass_range("");

    TFile *fIn_param = new TFile("signal_extraction_systematic_pdf.root", "READ");
    fIn_param->ls();

    // TH1D *Param_Pt_Charm = (TH1D *)fIn_param->Get("Param_Pt_unbinned_fromCharm");
    // TH1D *Param_Pt_Beauty = (TH1D *)fIn_param->Get("Param_Pt_unbinned_fromBeauty");
    // TH1D *Param_Pt_Mixed = (TH1D *)fIn_param->Get("Param_Pt_unbinned_fromMixed");

    // TH1D *Param_M_Charm = (TH1D *)fIn_param->Get("Param_M_unbinned_fromCharm");
    // TH1D *Param_M_Beauty = (TH1D *)fIn_param->Get("Param_M_unbinned_fromBeauty");
    // TH1D *Param_M_Mixed = (TH1D *)fIn_param->Get("Param_M_unbinned_fromMixed");

    TH1D *Param_Pt[n_DiMuSelection];
    TH1D *Param_M[n_DiMuSelection];

    // RooRealVar *B_DimuMass[n_DiMuSelection];
    // RooRealVar *n1_DimuMass[n_DiMuSelection];
    // RooRealVar *n2_DimuMass[n_DiMuSelection];
    // RooRealVar *B_DimuPt[n_DiMuSelection];
    // RooRealVar *n1_DimuPt[n_DiMuSelection];
    // RooRealVar *n2_DimuPt[n_DiMuSelection];
    RooRealVar *m = new RooRealVar("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", 4, 30);
    m->setBins(Binning_m);
    RooRealVar *pt = new RooRealVar("pt", "#it{p}_{T} (GeV/#it{c})", 0, 30);
    pt->setBins(Binning_pt);

    RooRealVar *B_DimuMass[n_DiMuSelection];
    RooRealVar *n1_DimuMass[n_DiMuSelection];
    RooRealVar *n2_DimuMass[n_DiMuSelection];
    RooRealVar *B_DimuPt[n_DiMuSelection];
    RooRealVar *n1_DimuPt[n_DiMuSelection];
    RooRealVar *n2_DimuPt[n_DiMuSelection];
    for (Int_t i = 0; i < n_DiMuSelection; i++)
    {
        printf("%s\n", Form("Param_Pt_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()));
        return;
        Param_Pt[i] = (TH1D *)fIn_param->Get(Form("Param_Pt_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()));
        Param_M[i] = (TH1D *)fIn_param->Get(Form("Param_M_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()));
        // printf("%s", Form("Param_M_unbinned_from%s", name_DiMuSelection[i].Data()));
        Param_M[i]->Draw();
        B_DimuMass[i] = new RooRealVar(Form("B_DimuMass_from%s", name_DiMuSelection[i].Data()), Form("B_DimuMass_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_M[i]->GetBinContent(1));

        n1_DimuMass[i] = new RooRealVar(Form("n1_DimuMass_from%s", name_DiMuSelection[i].Data()), Form("n1_DimuMass_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_M[i]->GetBinContent(2));
        n2_DimuMass[i] = new RooRealVar(Form("n2_DimuMass_from%s", name_DiMuSelection[i].Data()), Form("n2_DimuMass_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_M[i]->GetBinContent(3));

        // B_DimuMass[i]->setVal(Param_M[i]->GetBinContent(1));

        // n1_DimuMass[i]->setVal(Param_M[i]->GetBinContent(2));
        // n2_DimuMass[i]->setVal(Param_M[i]->GetBinContent(3));
        B_DimuMass[i]->setError(Param_M[i]->GetBinError(1));
        n1_DimuMass[i]->setError(Param_M[i]->GetBinError(2));
        n2_DimuMass[i]->setError(Param_M[i]->GetBinError(3));

        // return;
        B_DimuPt[i] = new RooRealVar(Form("B_DimuPt_from%s", name_DiMuSelection[i].Data()), Form("B_DimuPt_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_Pt[i]->GetBinContent(1));

        n1_DimuPt[i] = new RooRealVar(Form("n1_DimuPt_from%s", name_DiMuSelection[i].Data()), Form("n1_DimuPt_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_Pt[i]->GetBinContent(2));
        n2_DimuPt[i] = new RooRealVar(Form("n2_DimuPt_from%s", name_DiMuSelection[i].Data()), Form("n2_DimuMass_from%s", name_DiMuSelection[i].Data()), (Double_t)Param_Pt[i]->GetBinContent(3));

        B_DimuPt[i]->setError(Param_Pt[i]->GetBinError(1));
        n1_DimuPt[i]->setError(Param_Pt[i]->GetBinError(2));
        n2_DimuPt[i]->setError(Param_Pt[i]->GetBinError(3));

        printf("%s\n", name_DiMuSelection[i].Data());
        printf("B_DimuPt: %0.3f\n n1_DimuPt: %0.3f\n n2_DimuPt: %0.3f\n", B_DimuPt[i]->getError(), n1_DimuPt[i]->getError(), n2_DimuPt[i]->getError());

        // B_DimuMass[i] = w->var(Form("B_DimuMassFrom%s", name_DiMuSelection[i].Data()));

        // n1_DimuMass[i] = w->var(Form("n1_DimuMassFrom%s", name_DiMuSelection[i].Data()));

        // n2_DimuMass[i] = w->var(Form("n2_DimuMassFrom%s", name_DiMuSelection[i].Data()));

        // B_DimuPt[i] = w->var(Form("B_DimuPtFrom%s", name_DiMuSelection[i].Data()));

        // n1_DimuPt[i] = w->var(Form("n1_DimuPtFrom%s", name_DiMuSelection[i].Data()));

        // n2_DimuPt[i] = w->var(Form("n2_DimuPtFrom%s", name_DiMuSelection[i].Data()));
    }

    // return;

    RooWorkspace *w = new RooWorkspace("w", "workspace");
    w->factory(Form("PtMassExpPdf::pdfDimuMassFromcharm(m[4.0, 30.0], B_DimuMassFromCharm[%0.10f], n1_DimuMassFromCharm[%0.10f], n2_DimuMassFromCharm[%0.10f])", B_DimuMass[0]->getVal(), n1_DimuMass[0]->getVal(), n2_DimuMass[0]->getVal()));

    w->factory(Form("PtMassExpPdf::pdfDimuPtFromcharm(pt[0.0, 30.0], B_DimuPtFromCharm[%0.10f], n1_DimuPtFromCharm[%0.10f], n2_DimuPtFromCharm[%0.10f])", B_DimuPt[0]->getVal(), n1_DimuPt[0]->getVal(), n2_DimuPt[0]->getVal()));

    w->factory(Form("PtMassExpPdf::pdfDimuMassFrombeauty(m[4.0, 30.0], B_DimuMassFromBeauty[%0.10f], n1_DimuMassFromBeauty[%0.10f], n2_DimuMassFromBeauty[%0.10f])", B_DimuMass[1]->getVal(), n1_DimuMass[1]->getVal(), n2_DimuMass[1]->getVal()));

    w->factory(Form("PtMassExpPdf::pdfDimuPtFrombeauty(pt[0.0, 30.0], B_DimuPtFromBeauty[%0.10f], n1_DimuPtFromBeauty[%0.10f], n2_DimuPtFromBeauty[%0.10f])", B_DimuPt[1]->getVal(), n1_DimuPt[1]->getVal(), n2_DimuPt[1]->getVal()));

    w->factory(Form("PtMassExpPdf::pdfDimuMassFrommixed(m[4.0, 30.0], B_DimuMassFromMixed[%0.10f], n1_DimuMassFromMixed[%0.10f], n2_DimuMassFromMixed[%0.10f])", B_DimuMass[2]->getVal(), n1_DimuMass[2]->getVal(), n2_DimuMass[2]->getVal()));

    w->factory(Form("PtMassExpPdf::pdfDimuPtFrommixed(pt[0.0, 30.0], B_DimuPtFromMixed[%0.10f], n1_DimuPtFromMixed[%0.10f], n2_DimuPtFromMixed[%0.10f])", B_DimuPt[2]->getVal(), n1_DimuPt[2]->getVal(), n2_DimuPt[2]->getVal()));

    for (size_t i = 0; i < n_DiMuSelection; i++)
    {
        // B_DimuMass[i]->setError(Param_M[i]->GetBinError(1));
        // n1_DimuMass[i]->setError(Param_M[i]->GetBinError(2));
        // n2_DimuMass[i]->setError(Param_M[i]->GetBinError(3));

        // B_DimuPt[i]->setError(Param_Pt[i]->GetBinError(1));
        // n1_DimuPt[i]->setError(Param_Pt[i]->GetBinError(2));
        // n2_DimuPt[i]->setError(Param_Pt[i]->GetBinError(3));
        B_DimuMass[i]->setConstant(kTRUE);
        n1_DimuMass[i]->setConstant(kTRUE);
        n2_DimuMass[i]->setConstant(kTRUE);
        B_DimuPt[i]->setConstant(kTRUE);
        n1_DimuPt[i]->setConstant(kTRUE);
        n2_DimuPt[i]->setConstant(kTRUE);
    }

    // RooRealVar *B_DimuMassFromCharm = w->var("B_DimuMassFromCharm");
    // B_DimuMassFromCharm->setConstant(kTRUE);
    // RooRealVar *n1_DimuMassFromCharm = w->var("n1_DimuMassFromCharm");
    // n1_DimuMassFromCharm->setConstant(kTRUE);
    // RooRealVar *n2_DimuMassFromCharm = w->var("n2_DimuMassFromCharm");
    // n2_DimuMassFromCharm->setConstant(kTRUE);

    // RooRealVar *B_DimuMassFromBeauty = w->var("B_DimuMassFromBeauty");
    // B_DimuMassFromBeauty->setConstant(kTRUE);
    // RooRealVar *n1_DimuMassFromBeauty = w->var("n1_DimuMassFromBeauty");
    // n1_DimuMassFromBeauty->setConstant(kTRUE);
    // RooRealVar *n2_DimuMassFromBeauty = w->var("n2_DimuMassFromBeauty");
    // n2_DimuMassFromBeauty->setConstant(kTRUE);

    // RooRealVar *B_DimuMassFromMixed = w->var("B_DimuMassFromMixed");
    // B_DimuMassFromMixed->setConstant(kTRUE);
    // RooRealVar *n1_DimuMassFromMixed = w->var("n1_DimuMassFromMixed");
    // n1_DimuMassFromMixed->setConstant(kTRUE);
    // RooRealVar *n2_DimuMassFromMixed = w->var("n2_DimuMassFromMixed");
    // n2_DimuMassFromMixed->setConstant(kTRUE);

    // RooRealVar *B_DimuPtFromCharm = w->var("B_DimuPtFromCharm");
    // B_DimuPtFromCharm->setConstant(kTRUE);
    // RooRealVar *n1_DimuPtFromCharm = w->var("n1_DimuPtFromCharm");
    // n1_DimuPtFromCharm->setConstant(kTRUE);
    // RooRealVar *n2_DimuPtFromCharm = w->var("n2_DimuPtFromCharm");
    // n2_DimuPtFromCharm->setConstant(kTRUE);
    // // RooRealVar *C_DimuPtFromCharm = w->var("C_DimuPtFromCharm");
    // // C_DimuPtFromCharm->setConstant(kTRUE);
    // // RooRealVar *n3_DimuPtFromCharm = w->var("n3_DimuPtFromCharm");
    // // n3_DimuPtFromCharm->setConstant(kTRUE);

    // RooRealVar *B_DimuPtFromBeauty = w->var("B_DimuPtFromBeauty");
    // B_DimuPtFromBeauty->setConstant(kTRUE);
    // RooRealVar *n1_DimuPtFromBeauty = w->var("n1_DimuPtFromBeauty");
    // n1_DimuPtFromBeauty->setConstant(kTRUE);
    // RooRealVar *n2_DimuPtFromBeauty = w->var("n2_DimuPtFromBeauty");
    // n2_DimuPtFromBeauty->setConstant(kTRUE);
    // // RooRealVar *C_DimuPtFromBeauty = w->var("C_DimuPtFromBeauty");
    // // C_DimuPtFromBeauty->setConstant(kTRUE);
    // // RooRealVar *n3_DimuPtFromBeauty = w->var("n3_DimuPtFromBeauty");
    // // n3_DimuPtFromBeauty->setConstant(kTRUE);

    // RooRealVar *B_DimuPtFromMixed = w->var("B_DimuPtFromMixed");
    // B_DimuPtFromMixed->setConstant(kTRUE);
    // RooRealVar *n1_DimuPtFromMixed = w->var("n1_DimuPtFromMixed");
    // n1_DimuPtFromMixed->setConstant(kTRUE);
    // RooRealVar *n2_DimuPtFromMixed = w->var("n2_DimuPtFromMixed");
    // n2_DimuPtFromMixed->setConstant(kTRUE);
    // // RooRealVar *C_DimuPtFromMixed = w->var("C_DimuPtFromMixed");
    // // C_DimuPtFromMixed->setConstant(kTRUE);
    // // RooRealVar *n3_DimuPtFromMixed = w->var("n3_DimuPtFromMixed");
    // // n3_DimuPtFromMixed->setConstant(kTRUE);
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

    w->Print();

    // TFile *fIn = new TFile("/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Hist_MC/3_11_22/HF/Tree_HF_MCDimuHFTree_merged.root", "READ");
    // TTree *pt_tree_MC[n_DiMuSelection];
    // RooDataSet *Pt_Dimu[n_DiMuSelection];
    Int_t i = 0;

    // pt_tree_MC[i] = (TTree *)fIn->Get(Form("rec_tree_mu%s%s", "charm", mass_range.Data()));
    // Pt_Dimu[i] = new RooDataSet(Form("Pt_Dimu_%s", "Charm"), Form("Pt_Dimu_%s", "Charm"), RooArgSet(*pt), Import(*pt_tree_MC[i]));
    // auto result1 = pdfDimuPtFromCharm->fitTo(*Pt_Dimu[i], Minimizer("Minuit2"), Save(), SumW2Error(true));

    // RooPlot *prova = pt->frame(Title(Form("prova")));
    // // Pt_Dimu[i]->plotOn(prova,DrawOption("PE"));
    // pdfDimuPtFromCharm->plotOn(prova);
    // prova->Draw();

    Int_t choice = 0;
    RooCategory sample("sample", "sample");
    sample.defineType("mass");
    sample.defineType("transversemomentum");
    TFile *fIn_data = new TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/TreeResults_merged.root", "READ");

    // Taking data saved in tree
    TTree *tree_data = (TTree *)fIn_data->Get(Form("rec_data_tree%s", mass_range_data.Data()));
    RooDataSet *unbinned_M_Dimu_data = new RooDataSet("M_Dimu_data", "M_Dimu_data", RooArgSet(*m), Import(*tree_data));
    RooDataSet *unbinned_Pt_Dimu_data = new RooDataSet("Pt_Dimu_data", "Pt_Dimu_data", RooArgSet(*pt), Import(*tree_data)); // 5
    RooDataSet *unbinned_combData_set = new RooDataSet("combData", "combined data", RooArgSet(*m, *pt), Index(sample), Import("mass", *unbinned_M_Dimu_data), Import("transversemomentum", *unbinned_Pt_Dimu_data));

    RooPlot *frame = m->frame(Title("Imported TH1 with Poisson error bars"));
    pdfDimuMassFromCharm->plotOn(frame, LineColor(kMagenta + 2));
    pdfDimuMassFromBeauty->plotOn(frame, LineColor(kSpring - 6));
    pdfDimuMassFromMixed->plotOn(frame, LineColor(kAzure + 9));
    TCanvas *c2 = new TCanvas();
    c2->cd();
    frame->Draw();

    RooAddPdf *m_model;
    RooAddPdf *pt_model;

    RooRealVar *normForC;
    RooRealVar *normForB;
    RooRealVar *normForMixed;

    if (choice == 0)
    {
        ////////////////////
        // Fit with normalization

        normForC = new RooRealVar("n_charm_output", "number dimuon from c", 19474.5, 0, 200000);
        normForB = new RooRealVar("n_beauty_output", "number dimuon from b", 51888.9, 0, 200000);
        normForMixed = new RooRealVar("n_mixed_output", "number dimuon from b,c", 3160);
        normForMixed->setConstant(kTRUE);

        m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *pdfDimuMassFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
        pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pdfDimuPtFromMixed), RooArgList(*normForB, *normForC, *normForMixed));
    }

    else if (choice == 1)
    {
        ////////////////////
        // Fit with fraction
        normForC = new RooRealVar("fr_charm_output", "fraction dimuon from c", 0.35, 0., 1.);
        normForB = new RooRealVar("fr_beauty_output", "fraction dimuon from b", 0.605, 0., 1.);
        // RooRealVar *normForMixed = new RooRealVar("fr_mixed_output", "fraction dimuon from c", 0.045);
        // normForMixed->setConstant(kTRUE);

        RooFormulaVar *normForMixed = new RooFormulaVar("fr_mixed_output", "1-@0-@1", RooArgList(*normForC, *normForB));

        m_model = new RooAddPdf("m_model", " dimuMassFromB+dimuMassFromC + dimuMassFromMixed", RooArgList(*pdfDimuMassFromCharm, *pdfDimuMassFromBeauty, *pdfDimuMassFromMixed), RooArgList(*normForC, *normForB));
        pt_model = new RooAddPdf("pt_model", "dimuMassFromB+dimuMassFromC + dimuPtFromMixed", RooArgList(*pdfDimuPtFromCharm, *pdfDimuPtFromBeauty, *pdfDimuPtFromMixed), RooArgList(*normForC, *normForB));

        // RooAddPdf *m_model_signal = new RooAddPdf("m_model_signal", " (dimuMassFromMixed+dimuMassFromC)", RooArgList(*pdfDimuMassFromMixed), RooArgList(*normForMixed));
        // RooAddPdf *pt_model_signal = new RooAddPdf("pt_model_signal", " (dimuPtFromMixed+dimuPtFromC)", RooArgList(*pdfDimuPtFromMixed), RooArgList(*normForMixed));

        // m_model = new RooAddPdf("m_model", " dimuMassFromB+dimuMassFromC + dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty, *pdfDimuMassFromCharm, *m_model_signal), RooArgList(*normForB, *normForC));
        // pt_model = new RooAddPdf("pt_model", "dimuMassFromB+dimuMassFromC + dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty, *pdfDimuPtFromCharm, *pt_model_signal), RooArgList(*normForB, *normForC));

        // m_model = new RooAddPdf("m_model", " (dimuMassFromB+dimuMassFromC) + dimuMassFromMixed", RooArgList(*pdfDimuMassFromBeauty,*m_model_signal), RooArgList(*normForB));
        // pt_model = new RooAddPdf("pt_model", "(dimuMassFromB+dimuMassFromC) + dimuPtFromMixed", RooArgList(*pdfDimuPtFromBeauty,*pt_model_signal), RooArgList(*normForB));
    }
    // Define mass range in which do the fit

    m->setRange("mregion1", 4, 9);
    m->setRange("mregion2", 11, 30);

    // Define the pdf for the simultaneous fit

    RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
    simPdf.addPdf(*m_model, "mass");
    simPdf.addPdf(*pt_model, "transversemomentum");
    // simPdf.fitTo(*combData);
    //
    // RooFitResult *r = simPdf.fitTo(*unbinned_combData_set, Minimizer("Minuit2"), Save(), SumW2Error(true));
    RooFitResult *r = simPdf.fitTo(*unbinned_combData_set, Minimizer("Minuit2"), Range("mregion1", "mregion2"), Save(), SumW2Error(true));

    RooRealVar *fit_output[3];
    if (choice == 0)
    {
        fit_output[0] = (RooRealVar *)r->floatParsFinal().find("n_charm_output");
        fit_output[1] = (RooRealVar *)r->floatParsFinal().find("n_beauty_output");
        fit_output[2] = normForMixed;
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

    unbinned_combData_set->plotOn(m_frame, Name("combDatamass"), Cut("sample==sample::mass"), DrawOption("PEZ"));
    simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *unbinned_combData_set), Range("mregion1"), NormRange("mregion1,mregion2"), LineStyle(kSolid), LineColor(kRed));
    simPdf.plotOn(m_frame, Name("pdfmasscharm"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[0].Data())), ProjWData(sample, *unbinned_combData_set), Range("mregion1"), NormRange("mregion1,mregion2"), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
    simPdf.plotOn(m_frame, Name("pdfmassbeauty"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[1].Data())), ProjWData(sample, *unbinned_combData_set), Range("mregion1"), NormRange("mregion1,mregion2"), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
    simPdf.plotOn(m_frame, Name("pdfmassmixed"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[2].Data())), ProjWData(sample, *unbinned_combData_set), Range("mregion1"), NormRange("mregion1,mregion2"), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

    simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *unbinned_combData_set), Range("mregion2"), NormRange("mregion1,mregion2"), LineStyle(kSolid), LineColor(kRed));
    simPdf.plotOn(m_frame, Name("pdfmasscharm"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[0].Data())), ProjWData(sample, *unbinned_combData_set), Range("mregion2"), NormRange("mregion1,mregion2"), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
    simPdf.plotOn(m_frame, Name("pdfmassbeauty"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[1].Data())), ProjWData(sample, *unbinned_combData_set), Range("mregion2"), NormRange("mregion1,mregion2"), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
    simPdf.plotOn(m_frame, Name("pdfmassmixed"), Slice(sample, "mass"), Components(Form("%s", pdfMass_name[2].Data())), ProjWData(sample, *unbinned_combData_set), Range("mregion2"), NormRange("mregion1,mregion2"), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

    RooPlot *pt_frame = pt->frame(Title("pt_frame"));
    pt_frame->SetTitle(" ");
    pt_frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    pt_frame->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

    unbinned_combData_set->plotOn(pt_frame, Name("combDatapt"), Cut("sample==sample::transversemomentum"), DrawOption("PEZ"));
    simPdf.plotOn(pt_frame, Name("pdfpt"), NormRange(""), Slice(sample, "transversemomentum"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
    simPdf.plotOn(pt_frame, Name("pdfptcharm"), NormRange(""), Slice(sample, "transversemomentum"), Components(Form("%s", pdfPt_name[0].Data())), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(kMagenta - 2), LineWidth(5));
    simPdf.plotOn(pt_frame, Name("pdfptbeauty"), NormRange(""), Slice(sample, "transversemomentum"), Components(Form("%s", pdfPt_name[1].Data())), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(kSpring - 6), LineWidth(5));
    simPdf.plotOn(pt_frame, Name("pdfptmixed"), NormRange(""), Slice(sample, "transversemomentum"), Components(Form("%s", pdfPt_name[2].Data())), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(kAzure + 9), LineWidth(5));

    TString Pt_name[5] = {"combDatapt", "pdfpt", "pdfptcharm", "pdfptbeauty", "pdfptmixed"};

    TString Mass_name[5] = {"combDatamass", "pdfmass", "pdfmasscharm", "pdfmassbeauty", "pdfmassmixed"};

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

    TH1 *hDimuPt_data = unbinned_Pt_Dimu_data->createHistogram("h_ptdata", *pt, Binning(Binning_pt, 0, 30.0));
    hDimuPt_data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hDimuPt_data->GetYaxis()->SetTitle("d#it{N}/#it{p}_{T} (GeV/#it{c})^{-1}");

    TH1 *hDimuM_data = unbinned_M_Dimu_data->createHistogram("h_mdata", *m, Binning(Binning_m, 4.0, 30.0));
    ;
    hDimuM_data->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    hDimuM_data->GetYaxis()->SetTitle("d#it{N}/#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

    // Print the fit with ratio

    TString info;
    info.Form("unbinned_data_optfit%d_M%0.0f_%0.0f_BinM%d_BinPt%d", choice, 4.0, 30.0, Binning_m, Binning_pt);
    pt_frame->SetMaximum(1.2e+8);
    pt_frame->SetMinimum(1.5e-1);
    TCanvas *pt_canvas = printRooPlot_ratio(pt_frame, false, fit_output, choice, Pt_name, pt_modelFunc, hDimuPt_data, 0, 30.0);
    pt_canvas->SetName(Form("pt_canvas_%s", info.Data()));
    pt_canvas->SetTitle(Form("pt_canvas_%s", info.Data()));
    pt_canvas->SaveAs(Form("plot/%s.pdf", pt_canvas->GetName()));

    m_frame->SetMaximum(1.2e+7);
    m_frame->SetMinimum(1.5e-1);
    TCanvas *m_canvas = printRooPlot_ratio(m_frame, false, fit_output, choice, Mass_name, m_modelFunc, hDimuM_data, 4.0, 30.0);
    m_canvas->SetName(Form("m_canvas_%s", info.Data()));
    m_canvas->SetTitle(Form("m_canvas_%s", info.Data()));
    m_canvas->SaveAs(Form("plot/%s.pdf", m_canvas->GetName()));

    r->floatParsFinal().Print("s");

    for (Int_t i = 0; i < n_DiMuSelection; i++)
    {
        // B_DimuMass[i] = w->var(Form("B_DimuMassFrom%s", name_DiMuSelection[i].Data()));
        // B_DimuMass[i]->setConstant(kTRUE);
        // n1_DimuMass[i] = w->var(Form("n1_DimuMassFrom%s", name_DiMuSelection[i].Data()));
        // n1_DimuMass[i]->setConstant(kTRUE);
        // n2_DimuMass[i] = w->var(Form("n2_DimuMassFrom%s", name_DiMuSelection[i].Data()));
        // n2_DimuMass[i]->setConstant(kTRUE);

        // B_DimuPt[i] = w->var(Form("B_DimuPtFrom%s", name_DiMuSelection[i].Data()));
        // B_DimuPt[i]->setConstant(kTRUE);
        // n1_DimuPt[i] = w->var(Form("n1_DimuPtFrom%s", name_DiMuSelection[i].Data()));
        // n1_DimuPt[i]->setConstant(kTRUE);
        // n2_DimuPt[i] = w->var(Form("n2_DimuPtFrom%s", name_DiMuSelection[i].Data()));
        // n2_DimuPt[i]->setConstant(kTRUE);

        printf("%s\n", name_DiMuSelection[i].Data());
        printf("B_DimuPt: %0.3f\n n1_DimuPt: %0.3f\n n2_DimuPt: %0.3f\n", B_DimuPt[i]->getError(), n1_DimuPt[i]->getError(), n2_DimuPt[i]->getError());

        printf("B_DimuMass: %0.3f\n n1_DimuMass: %0.3f\n n2_DimuMass: %0.3f\n", B_DimuMass[i]->getError(), n1_DimuMass[i]->getError(), n2_DimuMass[i]->getError());
    }
}

void param_deviation()
{

    name_DiMuSelection[0].Form("Charm");
    name_DiMuSelection[1].Form("Beauty");
    name_DiMuSelection[2].Form("Mixed");

    TCanvas *mass_scaled_original_canvas[n_DiMuSelection];
    TCanvas *pt_scaled_original_canvas[n_DiMuSelection];

    TCanvas *control_Low2Up[n_DiMuSelection];
    TCanvas *control_Up2Low[n_DiMuSelection];

    TH2D *h_PtM_MC[n_DiMuSelection];
    TH1D *h_Pt_MC[n_DiMuSelection];
    TH1D *h_M_MC[n_DiMuSelection];

    TH1D *h_M_MC_scaled_Low2Up[n_DiMuSelection];
    TH1D *h_M_MC_scaled_Up2Low[n_DiMuSelection];

    TH1D *h_Pt_MC_scaled_Low2Up[n_DiMuSelection];
    TH1D *h_Pt_MC_scaled_Up2Low[n_DiMuSelection];

    Color_t color[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kAzure + 9};
    Color_t fillcolor[n_DiMuSelection] = {kMagenta - 10, kGreen - 10, kCyan - 10};

    TF1 *mass_pdf_linear_var_Low2Up[n_DiMuSelection];
    TF1 *mass_pdf_linear_var_Up2Low[n_DiMuSelection];

    TF1 *pt_pdf_linear_var_Low2Up[n_DiMuSelection];
    TF1 *pt_pdf_linear_var_Up2Low[n_DiMuSelection];

    TH1D *h_M_pdf_linear_var_Low2Up[n_DiMuSelection];
    TH1D *h_M_pdf_linear_var_Up2Low[n_DiMuSelection];

    TH1D *h_Pt_pdf_linear_var_Low2Up[n_DiMuSelection];
    TH1D *h_Pt_pdf_linear_var_Up2Low[n_DiMuSelection];

    Double_t B[n_DiMuSelection] = {2.5, 2.65, 3.43};
    Double_t n1[n_DiMuSelection] = {2.81, 2.81, 1.00};
    Double_t n2[n_DiMuSelection] = {2.5, 2.5, 7.5};

    TF1 *pdf_M[n_DiMuSelection];
    TF1 *pdf_Pt[n_DiMuSelection];

    TF1 *pdf_M_scaled_Low2Up[n_DiMuSelection];
    TF1 *pdf_M_scaled_Up2Low[n_DiMuSelection];

    TF1 *pdf_Pt_scaled_Low2Up[n_DiMuSelection];
    TF1 *pdf_Pt_scaled_Up2Low[n_DiMuSelection];

    TH1D *Param_Pt_original[n_DiMuSelection];
    TH1D *Param_Pt_scaled_Low2Up[n_DiMuSelection];
    TH1D *Param_Pt_scaled_Up2Low[n_DiMuSelection];

    TH1D *Param_M_original[n_DiMuSelection];
    TH1D *Param_M_scaled_Low2Up[n_DiMuSelection];
    TH1D *Param_M_scaled_Up2Low[n_DiMuSelection];

    TH1D *mass_hint_cl95[n_DiMuSelection];

    TH1D *pt_hint_cl95[n_DiMuSelection];

    TFile *fOut = new TFile("cl_95.root", "READ");
    fOut->cd();

    fOut->ls();

    TH1D *up_mass_hint_error_cl95[n_DiMuSelection];
    TH1D *lo_mass_hint_error_cl95[n_DiMuSelection];

    TH1D *up_pt_hint_error_cl95[n_DiMuSelection];
    TH1D *lo_pt_hint_error_cl95[n_DiMuSelection];

    Int_t i = 0;
    TFile *fOut_param = new TFile("signal_extraction_systematic_pdf.root", "UPDATE");

    for (size_t i = 0; i < n_DiMuSelection; i++)
    {
        up_mass_hint_error_cl95[i] = (TH1D *)fOut->Get(Form("up_mass_hint_error_%s_cl95", name_DiMuSelection[i].Data()));
        lo_mass_hint_error_cl95[i] = (TH1D *)fOut->Get(Form("lo_mass_hint_error_%s_cl95", name_DiMuSelection[i].Data()));

        up_pt_hint_error_cl95[i] = (TH1D *)fOut->Get(Form("up_pt_hint_error_%s_cl95", name_DiMuSelection[i].Data()));
        lo_pt_hint_error_cl95[i] = (TH1D *)fOut->Get(Form("lo_pt_hint_error_%s_cl95", name_DiMuSelection[i].Data()));

        // TLine *lo2up_line = new TLine(4, lo_mass_hint_error_cl95[i]->GetBinContent(1), 30, up_mass_hint_error_cl95[i]->GetBinContent(up_mass_hint_error_cl95[i]->GetNbinsX()));

        // TLine *up2lo_line = new TLine(4, up_mass_hint_error_cl95[i]->GetBinContent(1), 30, lo_mass_hint_error_cl95[i]->GetBinContent(lo_mass_hint_error_cl95[i]->GetNbinsX()));

        TFile *fIn = new TFile("/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Hist_MC/3_11_22/HF/HistLite_HF_MCDimuHFTree_merged.root", "READ");

        h_PtM_MC[i] = (TH2D *)fIn->Get(Form("DiMuon/M4/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_from%s", name_DiMuSelection[i].Data()));

        h_Pt_MC[i] = (TH1D *)h_PtM_MC[i]->ProjectionX();

        h_Pt_MC[i]->Rebin(10);
        h_Pt_MC[i]->Scale(1., "width");
        h_Pt_MC[i]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

        h_M_MC[i] = (TH1D *)h_PtM_MC[i]->ProjectionY();
        h_M_MC[i]->Rebin(10);
        h_M_MC[i]->Scale(1., "width");
        h_M_MC[i]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu #mu} (GeV/#it{c}^2)^{-1}");

        //----------------------------------------Pt------------------------------------------//
        //--------------------------Low 2 Up------------------------//
        Double_t pt_x1_Low2Up = 4;

        Double_t pt_y1_Low2Up = lo_pt_hint_error_cl95[i]->GetBinContent(1);

        Double_t pt_x2_Low2Up = 30;

        Double_t pt_y2_Low2Up = up_pt_hint_error_cl95[i]->GetBinContent(up_pt_hint_error_cl95[i]->GetNbinsX());

        cout << "pt_x1_Low2Up: " << pt_x1_Low2Up << "pt_y1_Low2Up: " << pt_y1_Low2Up << "pt_x2_Low2Up: " << pt_x2_Low2Up << "pt_y2_Low2Up: " << pt_y2_Low2Up << endl;

        Double_t pt_q = pt_y1_Low2Up - ((pt_y2_Low2Up - pt_y1_Low2Up) / (pt_x2_Low2Up - pt_x1_Low2Up)) * pt_x1_Low2Up;

        Double_t pt_m = (pt_y2_Low2Up - pt_y1_Low2Up) / (pt_x2_Low2Up - pt_x1_Low2Up);

        pt_pdf_linear_var_Low2Up[i] = new TF1(Form("pt_pdf_linear_var_Low2Up_from%s", name_DiMuSelection[i].Data()), "pol1", 4, 30);
        pt_pdf_linear_var_Low2Up[i]->FixParameter(0, pt_q);
        pt_pdf_linear_var_Low2Up[i]->FixParameter(1, pt_m);
        control_Low2Up[i] = new TCanvas(Form("control_Low2Up_from%s", name_DiMuSelection[i].Data()), Form("control_Low2Up_from%s", name_DiMuSelection[i].Data()), 600, 800);
        control_Low2Up[i]->cd();
        up_pt_hint_error_cl95[i]->Draw();

        pt_pdf_linear_var_Low2Up[i]->Draw("same");
        // lo2up_line->Draw("same");
        // return;

        //--------------------------Up 2 Low------------------------//

        pt_pdf_linear_var_Up2Low[i] = new TF1(Form("pt_pdf_linear_var_Up2Low_from%s", name_DiMuSelection[i].Data()), "pol1", 4, 30);
        pt_pdf_linear_var_Up2Low[i]->FixParameter(0, -pt_q);
        pt_pdf_linear_var_Up2Low[i]->FixParameter(1, -pt_m);
        // return;
        // control_Up2Low[i]= new TCanvas(Form("control_Low2Up_from%s", name_DiMuSelection[i].Data()), Form("control_Low2Up_from%s", name_DiMuSelection[i].Data()), 600, 800);
        // control_Up2Low[i]->cd();
        lo_pt_hint_error_cl95[i]->Draw("same");
        pt_pdf_linear_var_Up2Low[i]->Draw("same");
        // lo2up_line->Draw("same");

        //--------------------------Low 2 Up------------------------//
        h_Pt_pdf_linear_var_Low2Up[i] = (TH1D *)h_Pt_MC[i]->Clone(Form("h_Pt_pdf_linear_var_Low2Up_from%s", name_DiMuSelection[i].Data()));
        h_Pt_pdf_linear_var_Low2Up[i]->Reset();

        h_Pt_MC_scaled_Low2Up[i] = (TH1D *)h_Pt_MC[i]->Clone(Form("h_Pt_MC_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()));
        h_Pt_MC_scaled_Low2Up[i]->Reset();

        //--------------------------Up 2 Low------------------------//
        h_Pt_pdf_linear_var_Up2Low[i] = (TH1D *)h_Pt_MC[i]->Clone(Form("h_Pt_pdf_linear_var_Up2Low_from%s", name_DiMuSelection[i].Data()));
        h_Pt_pdf_linear_var_Up2Low[i]->Reset();

        h_Pt_MC_scaled_Up2Low[i] = (TH1D *)h_Pt_MC[i]->Clone(Form("h_Pt_MC_scaled_Up2Low_from%s", name_DiMuSelection[i].Data()));
        h_Pt_MC_scaled_Up2Low[i]->Reset();

        printf("--------------Low2Up------------\n");
        for (Int_t q = 1; q < h_Pt_MC[i]->GetNbinsX(); q++)
        {
            h_Pt_pdf_linear_var_Low2Up[i]->SetBinContent(q, h_Pt_MC[i]->GetBinContent(q) * pt_pdf_linear_var_Low2Up[i]->Eval(h_Pt_MC[i]->GetXaxis()->GetBinCenter(q)));

            Double_t new_value = h_Pt_MC[i]->GetBinContent(q) * (1 + pt_pdf_linear_var_Low2Up[i]->Eval(h_Pt_MC[i]->GetXaxis()->GetBinCenter(q)));

            h_Pt_MC_scaled_Low2Up[i]->SetBinContent(q, h_Pt_pdf_linear_var_Low2Up[i]->GetBinContent(q) + h_Pt_MC[i]->GetBinContent(q));

            printf("Bin Center %0.1f|| new_value %0.2f || Scaled content %0.2f || Original Content %0.2f || linear modification %0.2f\n", h_Pt_MC[i]->GetXaxis()->GetBinCenter(q), new_value, h_Pt_MC_scaled_Low2Up[i]->GetBinContent(q), h_Pt_MC[i]->GetBinContent(q), h_Pt_pdf_linear_var_Low2Up[i]->GetBinContent(q));
        }
        printf("--------------Up2Low------------\n");

        for (Int_t q = 1; q < h_Pt_MC[i]->GetNbinsX(); q++)
        {

            h_Pt_pdf_linear_var_Up2Low[i]->SetBinContent(q, h_Pt_MC[i]->GetBinContent(q) * pt_pdf_linear_var_Up2Low[i]->Eval(h_Pt_MC[i]->GetXaxis()->GetBinCenter(q)));

            Double_t new_value = h_Pt_MC[i]->GetBinContent(q) * (1 + pt_pdf_linear_var_Up2Low[i]->Eval(h_Pt_MC[i]->GetXaxis()->GetBinCenter(q)));

            h_Pt_MC_scaled_Up2Low[i]->SetBinContent(q, h_Pt_pdf_linear_var_Up2Low[i]->GetBinContent(q) + h_Pt_MC[i]->GetBinContent(q));

            printf("Bin Center %0.1f|| new_value %0.2f || Scaled content %0.2f || Original Content %0.2f || linear modification %0.2f\n", h_Pt_MC[i]->GetXaxis()->GetBinCenter(q), new_value, h_Pt_MC_scaled_Up2Low[i]->GetBinContent(q), h_Pt_MC[i]->GetBinContent(q), h_Pt_pdf_linear_var_Up2Low[i]->GetBinContent(q));
        }

        // pt_scaled_original_canvas[i] = histo2(i, h_Pt_MC_scaled_Low2Up[i], h_Pt_MC_scaled_Up2Low[i], h_Pt_MC[i], color[i], 0, 30);
        // pt_scaled_original_canvas[i]->SetName(Form("pt_scaled_original_canvas_from%s", name_DiMuSelection[i].Data()));
        // pt_scaled_original_canvas[i]->SetTitle(Form("pt_scaled_original_canvas_from%s", name_DiMuSelection[i].Data()));

        pdf_Pt[i] = new TF1(Form("pt_pdf_%s", name_DiMuSelection[i].Data()), FuncPtMass, 0, 30, 4);
        pdf_Pt[i]->SetParameter(3, 20000);
        pdf_Pt[i]->SetParameter(0, B[i]);
        pdf_Pt[i]->SetParameter(1, n1[i]);
        pdf_Pt[i]->SetParameter(2, n2[i]);

        h_Pt_MC[i]->Fit(pdf_Pt[i], "LR0I");

        pt_hint_cl95[i] = (TH1D *)h_Pt_MC[i]->Clone(Form("pt_hint_%s_cl95", name_DiMuSelection[i].Data()));
        pt_hint_cl95[i]->SetTitle(Form("pt_hint_%s_cl95", name_DiMuSelection[i].Data()));

        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(pt_hint_cl95[i], 0.99);

        Param_Pt_original[i] = new TH1D(Form("Param_Pt_original_from%s", name_DiMuSelection[i].Data()), Form("Param_Pt_original_from%s", name_DiMuSelection[i].Data()), 4, 0, 4);
        Param_Pt_original[i]->SetBinContent(1, pdf_Pt[i]->GetParameter(0));
        Param_Pt_original[i]->SetBinContent(2, pdf_Pt[i]->GetParameter(1));
        Param_Pt_original[i]->SetBinContent(3, pdf_Pt[i]->GetParameter(2));
        Param_Pt_original[i]->SetBinContent(4, pdf_Pt[i]->GetParameter(3));

        Param_Pt_original[i]->SetBinError(1, pdf_Pt[i]->GetParError(0));
        Param_Pt_original[i]->SetBinError(2, pdf_Pt[i]->GetParError(1));
        Param_Pt_original[i]->SetBinError(3, pdf_Pt[i]->GetParError(2));
        Param_Pt_original[i]->SetBinError(4, pdf_Pt[i]->GetParError(3));

        TCanvas *pippo = new TCanvas("pippo", "pippo", 500, 750);
        pippo->cd();

        Param_Pt_original[i]->Draw("Bar");

        // pdf_Pt[i]->Draw("SAME");

        pdf_Pt_scaled_Low2Up[i] = new TF1(Form("pdf_Pt_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()), FuncPtMass, 0, 30, 4);
        pdf_Pt_scaled_Low2Up[i]->SetParameter(3, 100000);
        pdf_Pt_scaled_Low2Up[i]->SetParameter(0, B[i]);
        pdf_Pt_scaled_Low2Up[i]->SetParameter(1, n1[i]);
        pdf_Pt_scaled_Low2Up[i]->SetParameter(2, n2[i]);

        h_Pt_MC_scaled_Low2Up[i]->Fit(pdf_Pt_scaled_Low2Up[i], "LR0I");

        Param_Pt_scaled_Low2Up[i] = new TH1D(Form("Param_Pt_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()), Form("Param_Pt_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()), 4, 0, 4);
        Param_Pt_scaled_Low2Up[i]->SetBinContent(1, pdf_Pt_scaled_Low2Up[i]->GetParameter(0));
        Param_Pt_scaled_Low2Up[i]->SetBinContent(2, pdf_Pt_scaled_Low2Up[i]->GetParameter(1));
        Param_Pt_scaled_Low2Up[i]->SetBinContent(3, pdf_Pt_scaled_Low2Up[i]->GetParameter(2));
        Param_Pt_scaled_Low2Up[i]->SetBinContent(4, pdf_Pt_scaled_Low2Up[i]->GetParameter(3));

        Param_Pt_scaled_Low2Up[i]->SetBinError(1, pdf_Pt_scaled_Low2Up[i]->GetParError(0));
        Param_Pt_scaled_Low2Up[i]->SetBinError(2, pdf_Pt_scaled_Low2Up[i]->GetParError(1));
        Param_Pt_scaled_Low2Up[i]->SetBinError(3, pdf_Pt_scaled_Low2Up[i]->GetParError(2));
        Param_Pt_scaled_Low2Up[i]->SetBinError(4, pdf_Pt_scaled_Low2Up[i]->GetParError(3));

        // pdf_Pt_scaled_Low2Up[i]->Draw("SAME");

        pdf_Pt_scaled_Up2Low[i] = new TF1(Form("pdf_Pt_scaled_Up2Low_from%s", name_DiMuSelection[i].Data()), FuncPtMass, 0, 30, 4);
        pdf_Pt_scaled_Up2Low[i]->SetParameter(3, 100000);
        pdf_Pt_scaled_Up2Low[i]->SetParameter(0, B[i]);
        pdf_Pt_scaled_Up2Low[i]->SetParameter(1, n1[i]);
        pdf_Pt_scaled_Up2Low[i]->SetParameter(2, n2[i]);

        h_Pt_MC_scaled_Up2Low[i]->Fit(pdf_Pt_scaled_Up2Low[i], "LR0I");

        Param_Pt_scaled_Up2Low[i] = new TH1D(Form("Param_Pt_scaled_Up2Low_from%s", name_DiMuSelection[i].Data()), Form("Param_Pt_scaled_Up2Low_from%s", name_DiMuSelection[i].Data()), 4, 0, 4);
        Param_Pt_scaled_Up2Low[i]->SetBinContent(1, pdf_Pt_scaled_Up2Low[i]->GetParameter(0));
        Param_Pt_scaled_Up2Low[i]->SetBinContent(2, pdf_Pt_scaled_Up2Low[i]->GetParameter(1));
        Param_Pt_scaled_Up2Low[i]->SetBinContent(3, pdf_Pt_scaled_Up2Low[i]->GetParameter(2));
        Param_Pt_scaled_Up2Low[i]->SetBinContent(4, pdf_Pt_scaled_Up2Low[i]->GetParameter(3));

        Param_Pt_scaled_Up2Low[i]->SetBinError(1, pdf_Pt_scaled_Up2Low[i]->GetParError(0));
        Param_Pt_scaled_Up2Low[i]->SetBinError(2, pdf_Pt_scaled_Up2Low[i]->GetParError(1));
        Param_Pt_scaled_Up2Low[i]->SetBinError(3, pdf_Pt_scaled_Up2Low[i]->GetParError(2));
        Param_Pt_scaled_Up2Low[i]->SetBinError(4, pdf_Pt_scaled_Up2Low[i]->GetParError(3));

        // pdf_Pt_scaled_Up2Low[i]->Draw("SAME");

        pt_scaled_original_canvas[i] = histo2(i, h_Pt_MC_scaled_Low2Up[i], h_Pt_MC_scaled_Up2Low[i], h_Pt_MC[i], pdf_Pt_scaled_Low2Up[i], pdf_Pt_scaled_Up2Low[i], pdf_Pt[i], pt_hint_cl95[i], color[i], fillcolor[i], 0, 30);
        pt_scaled_original_canvas[i]->SetName(Form("pt_scaled_original_canvas_from%s", name_DiMuSelection[i].Data()));
        pt_scaled_original_canvas[i]->SetTitle(Form("pt_scaled_original_canvas_from%s", name_DiMuSelection[i].Data()));
        pt_scaled_original_canvas[i]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/pt_scaled_original_canvas_from%s.png", name_DiMuSelection[i].Data()));

        //----------------------------------------M------------------------------------------//
        //--------------------------Low 2 Up------------------------//
        Double_t x1_Low2Up = 4;

        Double_t y1_Low2Up = lo_mass_hint_error_cl95[i]->GetBinContent(1);

        Double_t x2_Low2Up = 30;

        Double_t y2_Low2Up = up_mass_hint_error_cl95[i]->GetBinContent(up_mass_hint_error_cl95[i]->GetNbinsX());

        cout << "x1_Low2Up: " << x1_Low2Up << "y1_Low2Up: " << y1_Low2Up << "x2_Low2Up: " << x2_Low2Up << "y2_Low2Up: " << y2_Low2Up << endl;

        Double_t q = y1_Low2Up - ((y2_Low2Up - y1_Low2Up) / (x2_Low2Up - x1_Low2Up)) * x1_Low2Up;

        Double_t m = (y2_Low2Up - y1_Low2Up) / (x2_Low2Up - x1_Low2Up);

        mass_pdf_linear_var_Low2Up[i] = new TF1(Form("mass_pdf_linear_var_Low2Up_from%s", name_DiMuSelection[i].Data()), "pol1", 4, 30);
        mass_pdf_linear_var_Low2Up[i]->FixParameter(0, q);
        mass_pdf_linear_var_Low2Up[i]->FixParameter(1, m);
        control_Low2Up[i] = new TCanvas(Form("control_Low2Up_from%s", name_DiMuSelection[i].Data()), Form("control_Low2Up_from%s", name_DiMuSelection[i].Data()), 600, 800);
        control_Low2Up[i]->cd();
        up_mass_hint_error_cl95[i]->Draw();

        mass_pdf_linear_var_Low2Up[i]->Draw("same");
        // lo2up_line->Draw("same");

        //--------------------------Up 2 Low------------------------//

        mass_pdf_linear_var_Up2Low[i] = new TF1(Form("mass_pdf_linear_var_Up2Low_from%s", name_DiMuSelection[i].Data()), "pol1", 4, 30);
        mass_pdf_linear_var_Up2Low[i]->FixParameter(0, -q);
        mass_pdf_linear_var_Up2Low[i]->FixParameter(1, -m);
        // return;
        // control_Up2Low[i]= new TCanvas(Form("control_Low2Up_from%s", name_DiMuSelection[i].Data()), Form("control_Low2Up_from%s", name_DiMuSelection[i].Data()), 600, 800);
        // control_Up2Low[i]->cd();
        lo_mass_hint_error_cl95[i]->Draw("same");
        mass_pdf_linear_var_Up2Low[i]->Draw("same");
        // lo2up_line->Draw("same");

        //--------------------------Low 2 Up------------------------//
        h_M_pdf_linear_var_Low2Up[i] = (TH1D *)h_M_MC[i]->Clone(Form("h_M_pdf_linear_var_Low2Up_from%s", name_DiMuSelection[i].Data()));
        h_M_pdf_linear_var_Low2Up[i]->Reset();

        h_M_MC_scaled_Low2Up[i] = (TH1D *)h_M_MC[i]->Clone(Form("h_M_MC_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()));
        h_M_MC_scaled_Low2Up[i]->Reset();

        //--------------------------Up 2 Low------------------------//
        h_M_pdf_linear_var_Up2Low[i] = (TH1D *)h_M_MC[i]->Clone(Form("h_M_pdf_linear_var_Up2Low_from%s", name_DiMuSelection[i].Data()));
        h_M_pdf_linear_var_Up2Low[i]->Reset();

        h_M_MC_scaled_Up2Low[i] = (TH1D *)h_M_MC[i]->Clone(Form("h_M_MC_scaled_Up2Low_from%s", name_DiMuSelection[i].Data()));
        h_M_MC_scaled_Up2Low[i]->Reset();

        printf("--------------Low2Up------------\n");
        for (Int_t q = 1; q < h_M_MC[i]->GetNbinsX(); q++)
        {
            h_M_pdf_linear_var_Low2Up[i]->SetBinContent(q, h_M_MC[i]->GetBinContent(q) * mass_pdf_linear_var_Low2Up[i]->Eval(h_M_MC[i]->GetXaxis()->GetBinCenter(q)));

            Double_t new_value = h_M_MC[i]->GetBinContent(q) * (1 + mass_pdf_linear_var_Low2Up[i]->Eval(h_M_MC[i]->GetXaxis()->GetBinCenter(q)));

            h_M_MC_scaled_Low2Up[i]->SetBinContent(q, h_M_pdf_linear_var_Low2Up[i]->GetBinContent(q) + h_M_MC[i]->GetBinContent(q));

            printf("Bin Center %0.1f|| new_value %0.2f || Scaled content %0.2f || Original Content %0.2f || linear modification %0.2f\n", h_M_MC[i]->GetXaxis()->GetBinCenter(q), new_value, h_M_MC_scaled_Low2Up[i]->GetBinContent(q), h_M_MC[i]->GetBinContent(q), h_M_pdf_linear_var_Low2Up[i]->GetBinContent(q));
        }
        printf("--------------Up2Low------------\n");

        for (Int_t q = 1; q < h_M_MC[i]->GetNbinsX(); q++)
        {
            h_M_pdf_linear_var_Up2Low[i]->SetBinContent(q, h_M_MC[i]->GetBinContent(q) * mass_pdf_linear_var_Up2Low[i]->Eval(h_M_MC[i]->GetXaxis()->GetBinCenter(q)));

            Double_t new_value = h_M_MC[i]->GetBinContent(q) * (1 + mass_pdf_linear_var_Up2Low[i]->Eval(h_M_MC[i]->GetXaxis()->GetBinCenter(q)));

            h_M_MC_scaled_Up2Low[i]->SetBinContent(q, h_M_pdf_linear_var_Up2Low[i]->GetBinContent(q) + h_M_MC[i]->GetBinContent(q));

            printf("Bin Center %0.1f|| new_value %0.2f || Scaled content %0.2f || Original Content %0.2f || linear modification %0.2f\n", h_M_MC[i]->GetXaxis()->GetBinCenter(q), new_value, h_M_MC_scaled_Up2Low[i]->GetBinContent(q), h_M_MC[i]->GetBinContent(q), h_M_pdf_linear_var_Up2Low[i]->GetBinContent(q));
        }

        pdf_M[i] = new TF1(Form("mass_pdf_%s", name_DiMuSelection[i].Data()), FuncPtMass, 4, 30, 4);
        pdf_M[i]->SetParameter(3, 20000);
        pdf_M[i]->SetParameter(0, B[i]);
        pdf_M[i]->SetParameter(1, n1[i]);
        pdf_M[i]->SetParameter(2, n2[i]);

        h_M_MC[i]->Fit(pdf_M[i], "LR0I");

        mass_hint_cl95[i] = (TH1D *)h_M_MC[i]->Clone(Form("mass_hint_%s_cl95", name_DiMuSelection[i].Data()));
        mass_hint_cl95[i]->SetTitle(Form("mass_hint_%s_cl95", name_DiMuSelection[i].Data()));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mass_hint_cl95[i], 0.99);

        Param_M_original[i] = new TH1D(Form("Param_M_original_from%s", name_DiMuSelection[i].Data()), Form("Param_M_original_from%s", name_DiMuSelection[i].Data()), 4, 0, 4);
        Param_M_original[i]->SetBinContent(1, pdf_M[i]->GetParameter(0));
        Param_M_original[i]->SetBinContent(2, pdf_M[i]->GetParameter(1));
        Param_M_original[i]->SetBinContent(3, pdf_M[i]->GetParameter(2));
        Param_M_original[i]->SetBinContent(4, pdf_M[i]->GetParameter(3));

        Param_M_original[i]->SetBinError(1, pdf_M[i]->GetParError(0));
        Param_M_original[i]->SetBinError(2, pdf_M[i]->GetParError(1));
        Param_M_original[i]->SetBinError(3, pdf_M[i]->GetParError(2));
        Param_M_original[i]->SetBinError(4, pdf_M[i]->GetParError(3));

        // pdf_M[i]->Draw("SAME");

        pdf_M_scaled_Low2Up[i] = new TF1(Form("pdf_M_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()), FuncPtMass, 4, 30, 4);
        pdf_M_scaled_Low2Up[i]->SetParameter(3, 200000);
        pdf_M_scaled_Low2Up[i]->SetParameter(0, B[i]);
        pdf_M_scaled_Low2Up[i]->SetParameter(1, n1[i]);
        pdf_M_scaled_Low2Up[i]->SetParameter(2, n2[i]);

        h_M_MC_scaled_Low2Up[i]->Fit(pdf_M_scaled_Low2Up[i], "LR0I");

        Param_M_scaled_Low2Up[i] = new TH1D(Form("Param_M_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()), Form("Param_M_scaled_Low2Up_from%s", name_DiMuSelection[i].Data()), 4, 0, 4);
        Param_M_scaled_Low2Up[i]->SetBinContent(1, pdf_M_scaled_Low2Up[i]->GetParameter(0));
        Param_M_scaled_Low2Up[i]->SetBinContent(2, pdf_M_scaled_Low2Up[i]->GetParameter(1));
        Param_M_scaled_Low2Up[i]->SetBinContent(3, pdf_M_scaled_Low2Up[i]->GetParameter(2));
        Param_M_scaled_Low2Up[i]->SetBinContent(4, pdf_M_scaled_Low2Up[i]->GetParameter(3));

        Param_M_scaled_Low2Up[i]->SetBinError(1, pdf_M_scaled_Low2Up[i]->GetParError(0));
        Param_M_scaled_Low2Up[i]->SetBinError(2, pdf_M_scaled_Low2Up[i]->GetParError(1));
        Param_M_scaled_Low2Up[i]->SetBinError(3, pdf_M_scaled_Low2Up[i]->GetParError(2));
        Param_M_scaled_Low2Up[i]->SetBinError(4, pdf_M_scaled_Low2Up[i]->GetParError(3));

        // pdf_M_scaled_Low2Up[i]->Draw("SAME");

        pdf_M_scaled_Up2Low[i] = new TF1(Form("pdf_M_scaled_Up2Low_from%s", name_DiMuSelection[i].Data()), FuncPtMass, 4, 30, 4);
        pdf_M_scaled_Up2Low[i]->SetParameter(3, 200000);
        pdf_M_scaled_Up2Low[i]->SetParameter(0, B[i]);
        pdf_M_scaled_Up2Low[i]->SetParameter(1, n1[i]);
        pdf_M_scaled_Up2Low[i]->SetParameter(2, n2[i]);

        h_M_MC_scaled_Up2Low[i]->Fit(pdf_M_scaled_Up2Low[i], "LR0I");

        Param_M_scaled_Up2Low[i] = new TH1D(Form("Param_M_scaled_Up2Low_from%s", name_DiMuSelection[i].Data()), Form("Param_M_scaled_Up2Low_from%s", name_DiMuSelection[i].Data()), 4, 0, 4);
        Param_M_scaled_Up2Low[i]->SetBinContent(1, pdf_M_scaled_Up2Low[i]->GetParameter(0));
        Param_M_scaled_Up2Low[i]->SetBinContent(2, pdf_M_scaled_Up2Low[i]->GetParameter(1));
        Param_M_scaled_Up2Low[i]->SetBinContent(3, pdf_M_scaled_Up2Low[i]->GetParameter(2));
        Param_M_scaled_Up2Low[i]->SetBinContent(4, pdf_M_scaled_Up2Low[i]->GetParameter(3));

        Param_M_scaled_Up2Low[i]->SetBinError(1, pdf_M_scaled_Up2Low[i]->GetParError(0));
        Param_M_scaled_Up2Low[i]->SetBinError(2, pdf_M_scaled_Up2Low[i]->GetParError(1));
        Param_M_scaled_Up2Low[i]->SetBinError(3, pdf_M_scaled_Up2Low[i]->GetParError(2));
        Param_M_scaled_Up2Low[i]->SetBinError(4, pdf_M_scaled_Up2Low[i]->GetParError(3));

        // pdf_M_scaled_Up2Low[i]->Draw("SAME");
        mass_scaled_original_canvas[i] = histo2(i, h_M_MC_scaled_Low2Up[i], h_M_MC_scaled_Up2Low[i], h_M_MC[i], pdf_M_scaled_Low2Up[i], pdf_M_scaled_Up2Low[i], pdf_M[i], mass_hint_cl95[i], color[i], fillcolor[i], 0, 30);
        mass_scaled_original_canvas[i]->SetName(Form("mass_scaled_original_canvas_from%s", name_DiMuSelection[i].Data()));
        mass_scaled_original_canvas[i]->SetTitle(Form("mass_scaled_original_canvas_from%s", name_DiMuSelection[i].Data()));
        mass_scaled_original_canvas[i]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/mass_scaled_original_canvas_from%s.png", name_DiMuSelection[i].Data()));

        fOut_param->cd();

        Param_Pt_original[i]->Write(0, 2, 0);
        Param_Pt_scaled_Low2Up[i]->Write(0, 2, 0);
        Param_Pt_scaled_Up2Low[i]->Write(0, 2, 0);

        Param_M_original[i]->Write(0, 2, 0);
        Param_M_scaled_Low2Up[i]->Write(0, 2, 0);
        Param_M_scaled_Up2Low[i]->Write(0, 2, 0);
    }
}

TCanvas *histo2(Int_t i, TH1D *hist_scaled_Low2Up, TH1D *hist_scaled_Up2Low, TH1D *hist_original, TF1 *pdf_scaled_Low2Up, TF1 *pdf_scaled_Up2Low, TF1 *pdf_original, TH1D *hint, Color_t color, Color_t fillcolor, Int_t minx = 0, Int_t max_x = 30)
{
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 1000);
    canvas->SetTicks();
    canvas->cd();
    gPad->SetLogy(1);
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.11);

    hist_scaled_Low2Up->SetTitle(" ");
    hist_scaled_Low2Up->GetXaxis()->SetTitleOffset(1.2);
    hist_scaled_Low2Up->SetMinimum(0.3);

    // hist_scaled_Low2Up->GetYaxis()->SetRangeUser(hist_scaled_Low2Up->GetMinimum()*0.02,hist_scaled_Low2Up->GetMaximum()*4.);

    hist_scaled_Low2Up->SetMarkerStyle(22);
    hist_scaled_Low2Up->SetMarkerSize(1.5);
    hist_scaled_Low2Up->SetMarkerColor(color);

    hist_scaled_Up2Low->SetMarkerStyle(23);
    hist_scaled_Up2Low->SetMarkerSize(1.5);
    hist_scaled_Up2Low->SetMarkerColor(color);

    hist_original->SetMarkerStyle(46);
    hist_original->SetMarkerSize(1.5);
    hist_original->SetMarkerColor(color);

    canvas->cd();

    hist_scaled_Low2Up->Draw("PE");
    hint->SetFillColor(fillcolor);
    hint->Draw("e3gsame");
    pdf_scaled_Low2Up->SetLineColor(color);
    pdf_scaled_Low2Up->SetLineWidth(3);
    pdf_scaled_Low2Up->SetLineStyle(1);
    pdf_scaled_Low2Up->SetNpx(10000);
    pdf_scaled_Low2Up->Draw("SAME");
    hist_scaled_Low2Up->Draw("PESAME");

    

    pdf_scaled_Up2Low->SetLineColor(color);
    pdf_scaled_Up2Low->SetLineWidth(3);
    pdf_scaled_Up2Low->SetLineStyle(10);
    pdf_scaled_Up2Low->SetNpx(10000);
    pdf_scaled_Up2Low->Draw("SAME");
    hist_scaled_Up2Low->Draw("PESAME");

    pdf_original->SetLineColor(kRed);
    pdf_original->SetLineWidth(2);
    pdf_original->SetLineStyle(1);
    pdf_original->SetNpx(10000);
    pdf_original->Draw("SAME");
    hist_original->Draw("PESAME");

    TLegend *legend = new TLegend(0.515, 0.65, 0.85, 0.85);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0325);

    legend->SetTextAlign(12);
    legend->AddEntry(hist_scaled_Low2Up, "Scaled Low2Up", "PE");
    legend->AddEntry(hist_original, "Original ", "PE");
    legend->AddEntry(hist_scaled_Up2Low, "Scaled Up2Low", "PE");

    legend->AddEntry(pdf_scaled_Low2Up, "PDF Scaled Low2Up", "L");
    legend->AddEntry(pdf_original, "PDF Original ", "L");
    legend->AddEntry(pdf_scaled_Up2Low, "PDF Scaled Up2Low", "L");

    legend->AddEntry(hint, "99% CL original PDF", "F");

    legend->Draw();

    return canvas;
}

double FuncPtMass(double *x, double *par)
{
    return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
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