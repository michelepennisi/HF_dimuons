#include "/home/michele_pennisi/cernbox/HF_dimuons/common_include.h"
using namespace std;
const Int_t N_DiMu_sel = 3;
const Int_t N_variation = 3;

double FuncMass(double *x, double *par);
double FuncPt(double *x, double *par);

struct info
{
    TString old_syst_file = "";
    TString new_syst_file = "";

    TString Name_DiMu_sel[N_DiMu_sel] = {"Charm", "Beauty", "HF_Mixed"};
    TString Name_variation[N_variation] = {"original", "scaled_Low2Up", "scaled_Up2Low"};

    TString Bin_Label[3] = {"B", "n1", "n2"};

    Int_t Mass_Binning = 20;
    Int_t Low_Mass = 4;
    Int_t High_Mass = 30;
    Double_t LowM_cut = 8.;
    Double_t HighM_cut = 11.;
    Int_t Pt_Binning = 40;
    Int_t Low_Pt = 0;
    Int_t High_Pt = 30;

    Double_t HF_Mixed_fraction = 3.0;
    Double_t LF_HF_Mixed_fraction = 14.9;
    // Double_t LF_HF_Mixed_fraction = 0.;

    TString Generator = "PYTHIA";
    TString stat_MC = "full_stat";
    TString stat_Data = "LHC18p";
    TString LF_HF = "noLF_HF_LHC23i2";
    TString DY = "noDY";
    Color_t color[N_DiMu_sel] = {kMagenta + 2, kSpring - 6, kAzure + 9};
    Color_t fillcolor[N_DiMu_sel] = {kMagenta - 10, kGreen - 10, kCyan - 10};
};

void uncertainties()
{
}

void workspace()
{
    info opt;
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");

    TFile *fIn_param = new TFile(Form("systematic/template_modification_%s.root", opt.Generator.Data()), "READ");
    fIn_param->ls();

    RooWorkspace *w[N_variation];
    RooWorkspace *w_debug[N_variation];

    for (size_t p = 0; p < N_variation; p++)
    {
        TH1D *Param_Pt[N_DiMu_sel];
        TH1D *Param_M[N_DiMu_sel];

        RooRealVar *m = new RooRealVar("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", 4, 30);
        // m->setBins(Binning_m);
        RooRealVar *pt = new RooRealVar("pt", "#it{p}_{T} (GeV/#it{c})", 0, 30);
        // pt->setBins(Binning_pt);

        RooRealVar *B_DimuMass[N_DiMu_sel];
        RooRealVar *n1_DimuMass[N_DiMu_sel];
        RooRealVar *n2_DimuMass[N_DiMu_sel];
        RooRealVar *B_DimuPt[N_DiMu_sel];
        RooRealVar *n1_DimuPt[N_DiMu_sel];
        RooRealVar *n2_DimuPt[N_DiMu_sel];
        w[p] = new RooWorkspace(Form("w_%s", opt.Name_variation[p].Data()), Form("w_%s", opt.Name_variation[p].Data()));
        w_debug[p] = new RooWorkspace(Form("w_debug_%s", opt.Name_variation[p].Data()), Form("w_debug_%s", opt.Name_variation[p].Data()));

        for (Int_t i = 0; i < N_DiMu_sel; i++)
        {
            Param_Pt[i] = (TH1D *)fIn_param->Get(Form("param_syst_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/Param_Pt_%s_%s", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut, opt.Name_variation[p].Data(), opt.Name_DiMu_sel[i].Data()));
            Param_Pt[i]->Draw();
            Param_M[i] = (TH1D *)fIn_param->Get(Form("param_syst_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/Param_Mass_%s_%s", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut, opt.Name_variation[p].Data(), opt.Name_DiMu_sel[i].Data()));

            B_DimuMass[i] = new RooRealVar(Form("B_DimuMass_from%s", opt.Name_DiMu_sel[i].Data()), Form("B_DimuMass_from%s", opt.Name_DiMu_sel[i].Data()), (Double_t)Param_M[i]->GetBinContent(1));
            n1_DimuMass[i] = new RooRealVar(Form("n1_DimuMass_from%s", opt.Name_DiMu_sel[i].Data()), Form("n1_DimuMass_from%s", opt.Name_DiMu_sel[i].Data()), (Double_t)Param_M[i]->GetBinContent(2));
            n2_DimuMass[i] = new RooRealVar(Form("n2_DimuMass_from%s", opt.Name_DiMu_sel[i].Data()), Form("n2_DimuMass_from%s", opt.Name_DiMu_sel[i].Data()), (Double_t)Param_M[i]->GetBinContent(3));

            B_DimuMass[i]->setError(Param_M[i]->GetBinError(1));
            n1_DimuMass[i]->setError(Param_M[i]->GetBinError(2));
            n2_DimuMass[i]->setError(Param_M[i]->GetBinError(3));

            B_DimuPt[i] = new RooRealVar(Form("B_DimuPt_from%s", opt.Name_DiMu_sel[i].Data()), Form("B_DimuPt_from%s", opt.Name_DiMu_sel[i].Data()), (Double_t)Param_Pt[i]->GetBinContent(1));
            n1_DimuPt[i] = new RooRealVar(Form("n1_DimuPt_from%s", opt.Name_DiMu_sel[i].Data()), Form("n1_DimuPt_from%s", opt.Name_DiMu_sel[i].Data()), (Double_t)Param_Pt[i]->GetBinContent(2));
            n2_DimuPt[i] = new RooRealVar(Form("n2_DimuPt_from%s", opt.Name_DiMu_sel[i].Data()), Form("n2_DimuMass_from%s", opt.Name_DiMu_sel[i].Data()), (Double_t)Param_Pt[i]->GetBinContent(3));

            B_DimuPt[i]->setError(Param_Pt[i]->GetBinError(1));
            n1_DimuPt[i]->setError(Param_Pt[i]->GetBinError(2));
            n2_DimuPt[i]->setError(Param_Pt[i]->GetBinError(3));

            // B_DimuPt[i]->setConstant(kTRUE);
            // n1_DimuPt[i]->setConstant(kTRUE);
            // n2_DimuPt[i]->setConstant(kTRUE);

            printf("%s\n", opt.Name_DiMu_sel[i].Data());
            printf("B_DimuPt: %0.3f n1_DimuPt: %0.3f n2_DimuPt: %0.3f\n", B_DimuPt[i]->getVal(), n1_DimuPt[i]->getVal(), n2_DimuPt[i]->getVal());

            printf("B_DimuMass: %0.3f n1_DimuMass: %0.3f n2_DimuMass: %0.3f\n", B_DimuMass[i]->getVal(), n1_DimuMass[i]->getVal(), n2_DimuMass[i]->getVal());

            w[p]->factory(Form("pt[%d,%d], #it{p}_{T} (GeV/#it{c})", opt.Low_Pt, opt.High_Pt));
            w[p]->factory(Form("m[%d,%d], #it{m}_{#mu#mu} (GeV/#it{c}^{2})", opt.Low_Mass, opt.High_Mass));
            RooRealVar *m = w[p]->var("m");

            if (opt.Low_Mass == 4 && opt.High_Mass == 30)
            {
                w[p]->var("m")->setRange("low", 4, 8);
                w[p]->var("m")->setRange("high", 11, 30);
            }

            cout << Form("PtMassExpPdf::pdfDimuPtFrom%s(pt, B_DimuPtFrom%s[%0.10f], n1_DimuPtFrom%s[%0.10f], n2_DimuPtFrom%s[%0.10f])", opt.Name_DiMu_sel[i].Data(), opt.Name_DiMu_sel[i].Data(), B_DimuPt[i]->getVal(), opt.Name_DiMu_sel[i].Data(), n1_DimuPt[i]->getVal(), opt.Name_DiMu_sel[i].Data(), n2_DimuPt[i]->getVal()) << endl;

            w[p]->factory(Form("PtMassExpPdf::pdfDimuPtFrom%s(pt, B_DimuPtFrom%s[%0.10f], n1_DimuPtFrom%s[%0.10f], n2_DimuPtFrom%s[%0.10f])", opt.Name_DiMu_sel[i].Data(), opt.Name_DiMu_sel[i].Data(), B_DimuPt[i]->getVal(), opt.Name_DiMu_sel[i].Data(), n1_DimuPt[i]->getVal(), opt.Name_DiMu_sel[i].Data(), n2_DimuPt[i]->getVal()));

            cout << (Form("PtMassExpPdf::pdfDimuMassFrom%s(m, B_DimuMassFrom%s[%0.10f], n1_DimuMassFrom%s[%0.10f], n2_DimuMassFrom%s[%0.10f])", opt.Name_DiMu_sel[i].Data(), opt.Name_DiMu_sel[i].Data(), B_DimuMass[i]->getVal(), opt.Name_DiMu_sel[i].Data(), n1_DimuMass[i]->getVal(), opt.Name_DiMu_sel[i].Data(), n2_DimuMass[i]->getVal())) << endl;

            w[p]->factory(Form("PtMassExpPdf::pdfDimuMassFrom%s(m, B_DimuMassFrom%s[%0.10f], n1_DimuMassFrom%s[%0.10f], n2_DimuMassFrom%s[%0.10f])", opt.Name_DiMu_sel[i].Data(), opt.Name_DiMu_sel[i].Data(), B_DimuMass[i]->getVal(), opt.Name_DiMu_sel[i].Data(), n1_DimuMass[i]->getVal(), opt.Name_DiMu_sel[i].Data(), n2_DimuMass[i]->getVal()));
        }

        // return;

        w[p]->writeToFile(Form("systematic/syst_workspace_%s_M_%d_%d.root", opt.Generator.Data(), opt.Low_Mass, opt.High_Mass), kTRUE);
        w[p]->Print();
        // gDirectory->Add(w[p]);
    }

    return;
}

void converter()
{
    TString Version_ALIAOD = "Version_5_AliAOD_skimmed_fwd_fullstat";
    TString Dir_name = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output";
    info opt;

    TFile *fOut = new TFile(Form("systematic/hist_for_syst_%s_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f_.root", opt.Generator.Data(), opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut), "UPDATE");
    TString applied_cut;
    if (opt.Low_Mass == 4 && opt.High_Mass == 30)
        applied_cut.Form("((m > %d && m< %0.1f) || (m>%0.1f && m<%d)) && pt<%d", opt.Low_Mass, opt.LowM_cut, opt.HighM_cut, opt.High_Mass, opt.High_Pt);
    else
        applied_cut.Form("(m > %d && m< %d) && pt<%d", opt.Low_Mass, opt.High_Mass, opt.High_Pt);

    TFile *fIn;
    TString FileName;
    TString Tree_name;
    TTree *fTree;
    TH1F *Pt_prov_hist;
    TH1F *M_prov_hist;

    for (Int_t i_DiMu = 0; i_DiMu < N_DiMu_sel; i_DiMu++)
    {
        if (opt.Name_DiMu_sel[i_DiMu].Contains("Charm"))
        {
            if (opt.Generator.Contains("POWHEG"))
            {
                FileName.Form("%s/LHC23i1/%s/LHC23i1_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name.Form("DiMuon_Rec_PowhegOnly");
            }
            else if (opt.Generator.Contains("PYTHIA"))
            {
                FileName.Form("%s/LHC22b3/%s/LHC22b3_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name.Form("DiMuon_Rec");
            }
        }

        else if (opt.Name_DiMu_sel[i_DiMu].Contains("Beauty"))
        {
            if (opt.Generator.Contains("POWHEG"))
            {
                FileName.Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name.Form("DiMuon_Rec_PowhegOnly");
            }
            else if (opt.Generator.Contains("PYTHIA"))
            {
                FileName.Form("%s/LHC22b3/%s/LHC22b3_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name.Form("DiMuon_Rec");
            }
        }
        else if (opt.Name_DiMu_sel[i_DiMu].Contains("HF_Mixed"))
        {
            if (opt.Generator.Contains("POWHEG"))
            {
                FileName.Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name.Form("DiMuon_Rec_PythiaOnly");
            }
            else if (opt.Generator.Contains("PYTHIA"))
            {
                FileName.Form("%s/LHC22b3/%s/LHC22b3_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name.Form("DiMuon_Rec");
            }
        }

        fIn = new TFile(FileName, "READ");
        // fIn->ls();
        fTree = (TTree *)fIn->Get(Form("%s_%s", Tree_name.Data(), opt.Name_DiMu_sel[i_DiMu].Data()));
        std::cout << "Component: " << opt.Name_DiMu_sel[i_DiMu].Data() << " ||file: " << FileName.Data() << " || tree name: " << Tree_name.Data() << endl;

        fOut->cd();

        if (!fOut->GetDirectory(TString::Format("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut)))
            fOut->mkdir(TString::Format("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut));

        fOut->cd(TString::Format("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut));
        Pt_prov_hist = new TH1F(Form("h_Pt_from%s", opt.Name_DiMu_sel[i_DiMu].Data()), "; #it{p}_{T} (GeV/#it{c})", (opt.High_Pt - opt.Low_Pt) * 10, opt.Low_Pt, opt.High_Pt);
        M_prov_hist = new TH1F(Form("h_Mass_from%s", opt.Name_DiMu_sel[i_DiMu].Data()), "; #it{m}_{#mu#mu} (GeV/#it{c}^{2})", (opt.High_Mass - opt.Low_Mass) * 10, opt.Low_Mass, opt.High_Mass);

        fTree->Draw(Form("pt>>h_Pt_from%s", opt.Name_DiMu_sel[i_DiMu].Data()), applied_cut, "goff");
        fTree->Draw(Form("m>>h_Mass_from%s", opt.Name_DiMu_sel[i_DiMu].Data()), applied_cut, "goff");

        Pt_prov_hist->Write(0, 2, 0);
        M_prov_hist->Write(0, 2, 0);
    }
    fOut->Close();
    fIn->Close();

    fOut = new TFile(Form("systematic/hist_for_syst_%s_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f_.root", opt.Generator.Data(), opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut), "READ");
    new TBrowser();
}

TH1F *Linear_deviation(TH1F *h_CL)
{
    TH1F *deviation = (TH1F *)h_CL->Clone("deviation");
    deviation->Reset();
    for (Int_t i_bin = 0; i_bin < h_CL->GetNbinsX(); i_bin++)
    {
        if (h_CL->GetBinContent(i_bin + 1) != 0)
        {
            deviation->SetBinContent(i_bin + 1, (Double_t)h_CL->GetBinError(i_bin + 1) / h_CL->GetBinContent(i_bin + 1));
            deviation->SetBinError(i_bin + 1, 0.1 * (Double_t)h_CL->GetBinError(i_bin + 1) / h_CL->GetBinContent(i_bin + 1));
        }
    }
    return deviation;
}

TH1F *MC_deviation(TH1F *Original_MC, TF1 *linear_modification)
{

    TH1F *Modified_MC = (TH1F *)Original_MC->Clone("Modified_MC");
    Modified_MC->Reset();
    for (Int_t bin = 1; bin <= Modified_MC->GetNbinsX(); bin++)
        Modified_MC->SetBinContent(bin, Original_MC->GetBinContent(bin) * (1 + linear_modification->Eval(Original_MC->GetXaxis()->GetBinCenter(bin))));

    return Modified_MC;
}

void param_deviation()
{
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
    info opt;
    TCanvas *mass_scaled_original_canvas[N_DiMu_sel];
    TCanvas *pt_scaled_original_canvas[N_DiMu_sel];

    TCanvas *control_Low2Up[N_DiMu_sel];
    TCanvas *control_Up2Low[N_DiMu_sel];

    TH1F *h_Pt_MC[N_DiMu_sel];
    TH1F *h_Mass_MC[N_DiMu_sel];

    TH1F *h_Mass_MC_scaled_Low2Up[N_DiMu_sel];
    TH1F *h_Mass_MC_scaled_Up2Low[N_DiMu_sel];

    TH1F *h_Pt_MC_scaled_Low2Up[N_DiMu_sel];
    TH1F *h_Pt_MC_scaled_Up2Low[N_DiMu_sel];

    TF1 *Mass_pdf_linear_var_Low2Up[N_DiMu_sel];
    TF1 *Mass_pdf_linear_var_Up2Low[N_DiMu_sel];

    TF1 *Pt_pdf_linear_var_Low2Up[N_DiMu_sel];
    TF1 *Pt_pdf_linear_var_Up2Low[N_DiMu_sel];

    TH1F *h_Mass_pdf_linear_var_Low2Up[N_DiMu_sel];
    TH1F *h_Mass_pdf_linear_var_Up2Low[N_DiMu_sel];

    TH1F *h_Pt_pdf_linear_var_Low2Up[N_DiMu_sel];
    TH1F *h_Pt_pdf_linear_var_Up2Low[N_DiMu_sel];

    Double_t par3_Pt[N_DiMu_sel];
    Double_t par3_Mass[N_DiMu_sel];
    Double_t B_Pt[N_DiMu_sel];
    Double_t n1_Pt[N_DiMu_sel];
    Double_t n2_Pt[N_DiMu_sel];
    Double_t B_Mass[N_DiMu_sel];
    Double_t n1_Mass[N_DiMu_sel];
    Double_t n2_Mass[N_DiMu_sel];

    if (opt.Generator.Contains("PYTHIA"))
    {
        par3_Pt[0] = {20000};
        par3_Pt[1] = {20000};
        par3_Pt[2] = {20000};

        par3_Mass[0] = {10000};
        par3_Mass[1] = {10000};
        par3_Mass[2] = {10000};

        B_Pt[0] = {2.5};
        B_Pt[1] = {2.5};
        B_Pt[2] = {2.5};

        n1_Pt[0] = {2.81};
        n1_Pt[1] = {2.81};
        n1_Pt[2] = {2.81};

        n2_Pt[0] = {2.5};
        n2_Pt[1] = {2.5};
        n2_Pt[2] = {2.5};

        B_Mass[0] = {2.5};
        B_Mass[1] = {2.5};
        B_Mass[2] = {2.5};

        n1_Mass[0] = {2.81};
        n1_Mass[1] = {2.81};
        n1_Mass[2] = {2.81};

        n2_Mass[0] = {2.5};
        n2_Mass[1] = {2.5};
        n2_Mass[2] = {2.5};
    }
    else
    {
        par3_Pt[0] = {20000};
        par3_Pt[1] = {200000};
        par3_Pt[2] = {20000};
        if (opt.Low_Mass == 4 && opt.High_Mass == 30)
            par3_Mass[0] = {650000};
        else
            par3_Mass[0] = {20000};
        par3_Mass[1] = {10000};
        par3_Mass[2] = {10000};

        B_Pt[0] = {2.5};
        B_Pt[1] = {4.2133};
        B_Pt[2] = {5.9703};

        n1_Pt[0] = {2.81};
        n1_Pt[1] = {2.0563};
        n1_Pt[2] = {1.4663e+00};

        n2_Pt[0] = {2.5};
        n2_Pt[1] = {2.5386};
        n2_Pt[2] = {1.6398e+01};

        B_Mass[0] = {2.5};
        B_Mass[1] = {2.5};
        B_Mass[2] = {4.5224e+00};

        n1_Mass[0] = {2.81};
        n1_Mass[1] = {2.81};
        n1_Mass[2] = {1.6422e+00};

        n2_Mass[0] = {2.5};
        n2_Mass[1] = {2.5};
        n2_Mass[2] = {6.0236e+00};

        // par3_Pt[0] = {20000};
        // par3_Pt[1] = {200000};
        // par3_Pt[2] = {20000};

        // par3_Mass[0] = {10000};
        // par3_Mass[1] = {10000};
        // par3_Mass[2] = {10000};

        // B_Pt[0] = {2.5};
        // B_Pt[1] = {4.2133};
        // B_Pt[2] = {5.9703};

        // n1_Pt[0] = {2.81};
        // n1_Pt[1] = {2.0563};
        // n1_Pt[2] = {1.4663e+00};

        // n2_Pt[0] = {2.5};
        // n2_Pt[1] = {2.5386};
        // n2_Pt[2] = {1.6398e+01};

        // B_Mass[0] = {2.5};
        // B_Mass[1] = {2.5};
        // B_Mass[2] = {4.5224e+00};

        // n1_Mass[0] = {2.81};
        // n1_Mass[1] = {2.81};
        // n1_Mass[2] = {1.6422e+00};

        // n2_Mass[0] = {2.5};
        // n2_Mass[1] = {2.5};
        // n2_Mass[2] = {6.0236e+00};
    }

    TF1 *pdf_Mass[N_DiMu_sel];
    TF1 *pdf_Pt[N_DiMu_sel];

    TF1 *pdf_Mass_scaled_Low2Up[N_DiMu_sel];
    TF1 *pdf_Mass_scaled_Up2Low[N_DiMu_sel];

    TF1 *pdf_Pt_scaled_Low2Up[N_DiMu_sel];
    TF1 *pdf_Pt_scaled_Up2Low[N_DiMu_sel];

    TH1F *Param_Pt_original[N_DiMu_sel];
    TH1F *Param_Pt_scaled_Low2Up[N_DiMu_sel];
    TH1F *Param_Pt_scaled_Up2Low[N_DiMu_sel];

    TH1F *Param_Mass_original[N_DiMu_sel];
    TH1F *Param_Mass_scaled_Low2Up[N_DiMu_sel];
    TH1F *Param_Mass_scaled_Up2Low[N_DiMu_sel];

    // TFile *fOut = new TFile("test_old_extraction.root", "READ");
    TH1F *up_Mass_hint_error_cl95[N_DiMu_sel];
    TH1F *lo_Mass_hint_error_cl95[N_DiMu_sel];

    TH1F *up_Pt_hint_error_cl95[N_DiMu_sel];
    TH1F *lo_Pt_hint_error_cl95[N_DiMu_sel];

    TFile *fOut_param = new TFile(Form("systematic/template_modification_%s.root", opt.Generator.Data()), "UPDATE");
    TString saving_dir = Form("param_syst_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut);
    if (!fOut_param->GetDirectory(saving_dir))
        fOut_param->mkdir(saving_dir);

    TFile *fIn = new TFile(Form("systematic/hist_for_syst_%s_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f_.root", opt.Generator.Data(), opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut), "READ");
    fIn->ls();
    TH1F *Pt_Conf_interval;
    TH1F *M_Conf_interval;
    for (size_t i = 0; i < N_DiMu_sel; i++)
    {
        //----------------------------Pt------------------------------------------------//
        h_Pt_MC[i] = (TH1F *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/h_Pt_from%s", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut, opt.Name_DiMu_sel[i].Data()));
        h_Pt_MC[i]->Rebin(10);
        h_Pt_MC[i]->Scale(1., "width");
        hist1D_graphic_opt(h_Pt_MC[i], kFALSE, 1, 24, opt.color[i], 1.);
        h_Pt_MC[i]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
        h_Pt_MC[i]->SetTitle("Original");
        pdf_Pt[i] = new TF1(Form("pt_pdf_%s", opt.Name_DiMu_sel[i].Data()), FuncPt, opt.Low_Pt, opt.High_Pt, 4);
        pdf_Pt[i]->SetTitle(Form("%s original PDF", opt.Name_DiMu_sel[i].Data()));
        pdf_Pt[i]->SetParameter(3, par3_Pt[i]);
        pdf_Pt[i]->SetParameter(0, B_Pt[i]);
        pdf_Pt[i]->SetParameter(1, n1_Pt[i]);
        pdf_Pt[i]->SetParameter(2, n2_Pt[i]);
        pdf_Pt[i]->SetLineWidth(4);
        pdf_Pt[i]->SetNpx(200);
        pdf_Pt[i]->SetLineColor(opt.fillcolor[i]);
        h_Pt_MC[i]->Fit(pdf_Pt[i], "LR0I");
        h_Pt_MC[i]->Draw("PE");
        pdf_Pt[i]->Draw("same");

        Param_Pt_original[i] = new TH1F(Form("Param_Pt_original_%s", opt.Name_DiMu_sel[i].Data()), Form("Param_Pt_original_%s", opt.Name_DiMu_sel[i].Data()), 4, 0, 4);
        for (Int_t bin = 0; bin < Param_Pt_original[i]->GetNbinsX(); bin++)
        {
            Param_Pt_original[i]->SetBinContent(bin + 1, pdf_Pt[i]->GetParameter(bin));
            Param_Pt_original[i]->SetBinError(bin + 1, pdf_Pt[i]->GetParError(bin));
        }

        fOut_param->cd(saving_dir);
        Param_Pt_original[i]->Write(0, 2, 0);

        Pt_Conf_interval = (TH1F *)h_Pt_MC[i]->Clone(Form("h_Pt_CL_95_%s", opt.Name_DiMu_sel[i].Data()));
        Pt_Conf_interval->Reset();
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(Pt_Conf_interval, 0.95);
        Pt_Conf_interval->SetTitle("95% CL");

        up_Pt_hint_error_cl95[i] = Linear_deviation(Pt_Conf_interval);
        up_Pt_hint_error_cl95[i]->SetName(Form("up_pt_hint_error_cl95_%s_cl95", opt.Name_DiMu_sel[i].Data()));

        lo_Pt_hint_error_cl95[i] = (TH1F *)up_Pt_hint_error_cl95[i]->Clone(Form("lo_pt_hint_error_cl95_%s_cl95", opt.Name_DiMu_sel[i].Data()));
        lo_Pt_hint_error_cl95[i]->Reset();
        for (Int_t bin = 0; bin < lo_Pt_hint_error_cl95[i]->GetNbinsX(); bin++)
        {
            lo_Pt_hint_error_cl95[i]->SetBinContent(bin + 1, -(up_Pt_hint_error_cl95[i]->GetBinContent(bin + 1)));
            lo_Pt_hint_error_cl95[i]->SetBinError(bin + 1, (up_Pt_hint_error_cl95[i]->GetBinError(bin + 1)));
        }

        //--------------------------Low 2 Up------------------------//
        Double_t pt_x1_Low2Up = opt.Low_Pt;

        Double_t pt_y1_Low2Up = lo_Pt_hint_error_cl95[i]->GetBinContent(1);

        Double_t pt_x2_Low2Up = opt.High_Pt;

        Double_t pt_y2_Low2Up = up_Pt_hint_error_cl95[i]->GetBinContent(up_Pt_hint_error_cl95[i]->GetNbinsX());

        cout << "pt_x1_Low2Up: " << pt_x1_Low2Up << "pt_y1_Low2Up: " << pt_y1_Low2Up << "pt_x2_Low2Up: " << pt_x2_Low2Up << "pt_y2_Low2Up: " << pt_y2_Low2Up << endl;

        Double_t pt_q = pt_y1_Low2Up - ((pt_y2_Low2Up - pt_y1_Low2Up) / (pt_x2_Low2Up - pt_x1_Low2Up)) * pt_x1_Low2Up;

        Double_t pt_m = (pt_y2_Low2Up - pt_y1_Low2Up) / (pt_x2_Low2Up - pt_x1_Low2Up);

        Pt_pdf_linear_var_Low2Up[i] = new TF1(Form("pt_pdf_linear_var_Low2Up_from%s", opt.Name_DiMu_sel[i].Data()), "pol1", opt.Low_Pt, opt.High_Pt);
        Pt_pdf_linear_var_Low2Up[i]->FixParameter(0, pt_q);
        Pt_pdf_linear_var_Low2Up[i]->FixParameter(1, pt_m);

        h_Pt_MC_scaled_Low2Up[i] = MC_deviation(h_Pt_MC[i], Pt_pdf_linear_var_Low2Up[i]);
        hist1D_graphic_opt(h_Pt_MC_scaled_Low2Up[i], kFALSE, 1, 22, opt.color[i], 1.);
        h_Pt_MC_scaled_Low2Up[i]->SetTitle("Low to Up");
        pdf_Pt_scaled_Low2Up[i] = new TF1(Form("pdf_Pt_scaled_Low2Up_from%s", opt.Name_DiMu_sel[i].Data()), FuncPt, opt.Low_Pt, opt.High_Pt, 4);
        pdf_Pt_scaled_Low2Up[i]->SetTitle(Form("%s Low2Up mod. PDF", opt.Name_DiMu_sel[i].Data()));
        pdf_Pt_scaled_Low2Up[i]->SetParameter(3, pdf_Pt[i]->GetParameter(3));
        pdf_Pt_scaled_Low2Up[i]->SetParameter(0, pdf_Pt[i]->GetParameter(0));
        pdf_Pt_scaled_Low2Up[i]->SetParameter(1, pdf_Pt[i]->GetParameter(1));
        pdf_Pt_scaled_Low2Up[i]->SetParameter(2, pdf_Pt[i]->GetParameter(2));
        pdf_Pt_scaled_Low2Up[i]->SetLineColor(opt.color[i]);
        pdf_Pt_scaled_Low2Up[i]->SetLineWidth(4);
        pdf_Pt_scaled_Low2Up[i]->SetNpx(200);
        pdf_Pt_scaled_Low2Up[i]->SetLineStyle(kDotted);

        h_Pt_MC_scaled_Low2Up[i]->Fit(pdf_Pt_scaled_Low2Up[i], "LR0I");
        Param_Pt_scaled_Low2Up[i] = new TH1F(Form("Param_Pt_scaled_Low2Up_%s", opt.Name_DiMu_sel[i].Data()), Form("Param_Pt_scaled_Low2Up_%s", opt.Name_DiMu_sel[i].Data()), 4, 0, 4);
        for (Int_t bin = 0; bin < Param_Pt_scaled_Low2Up[i]->GetNbinsX(); bin++)
        {
            Param_Pt_scaled_Low2Up[i]->SetBinContent(bin + 1, pdf_Pt_scaled_Low2Up[i]->GetParameter(bin));
            Param_Pt_scaled_Low2Up[i]->SetBinError(bin + 1, pdf_Pt_scaled_Low2Up[i]->GetParError(bin));
        }
        fOut_param->cd(saving_dir);
        Param_Pt_scaled_Low2Up[i]->Write(0, 2, 0);

        //--------------------------Up 2 Low------------------------//

        Pt_pdf_linear_var_Up2Low[i] = new TF1(Form("pt_pdf_linear_var_Up2Low_from%s", opt.Name_DiMu_sel[i].Data()), "pol1", opt.Low_Pt, opt.High_Pt);
        Pt_pdf_linear_var_Up2Low[i]->FixParameter(0, -pt_q);
        Pt_pdf_linear_var_Up2Low[i]->FixParameter(1, -pt_m);

        h_Pt_MC_scaled_Up2Low[i] = MC_deviation(h_Pt_MC[i], Pt_pdf_linear_var_Up2Low[i]);
        hist1D_graphic_opt(h_Pt_MC_scaled_Up2Low[i], kFALSE, 1, 23, opt.color[i], 1.);
        h_Pt_MC_scaled_Up2Low[i]->SetTitle("Up to Low");
        pdf_Pt_scaled_Up2Low[i] = new TF1(Form("pdf_Pt_scaled_Up2Low_from%s", opt.Name_DiMu_sel[i].Data()), FuncPt, opt.Low_Pt, opt.High_Pt, 4);
        pdf_Pt_scaled_Up2Low[i]->SetTitle(Form("%s Up2Low mod. PDF", opt.Name_DiMu_sel[i].Data()));
        pdf_Pt_scaled_Up2Low[i]->SetParameter(3, pdf_Pt[i]->GetParameter(3));
        pdf_Pt_scaled_Up2Low[i]->SetParameter(0, pdf_Pt[i]->GetParameter(0));
        pdf_Pt_scaled_Up2Low[i]->SetParameter(1, pdf_Pt[i]->GetParameter(1));
        pdf_Pt_scaled_Up2Low[i]->SetParameter(2, pdf_Pt[i]->GetParameter(2));
        pdf_Pt_scaled_Up2Low[i]->SetLineColor(opt.color[i]);
        pdf_Pt_scaled_Up2Low[i]->SetLineWidth(4);
        pdf_Pt_scaled_Up2Low[i]->SetNpx(200);
        pdf_Pt_scaled_Up2Low[i]->SetLineStyle(kDashed);

        h_Pt_MC_scaled_Up2Low[i]->Fit(pdf_Pt_scaled_Up2Low[i], "LR0I");
        Param_Pt_scaled_Up2Low[i] = new TH1F(Form("Param_Pt_scaled_Up2Low_%s", opt.Name_DiMu_sel[i].Data()), Form("Param_Pt_scaled_Up2Low_%s", opt.Name_DiMu_sel[i].Data()), 4, 0, 4);
        for (Int_t bin = 0; bin < Param_Pt_scaled_Up2Low[i]->GetNbinsX(); bin++)
        {
            Param_Pt_scaled_Up2Low[i]->SetBinContent(bin + 1, pdf_Pt_scaled_Up2Low[i]->GetParameter(bin));
            Param_Pt_scaled_Up2Low[i]->SetBinError(bin + 1, pdf_Pt_scaled_Up2Low[i]->GetParError(bin));
        }

        fOut_param->cd(saving_dir);
        Param_Pt_scaled_Up2Low[i]->Write(0, 2, 0);

        TCanvas *C_Pt_linear_deviation = canvas_noratio(Form("C_Pt_linear_deviation_%s", opt.Name_DiMu_sel[i].Data()));
        hist1D_graphic_opt(up_Pt_hint_error_cl95[i], kFALSE, 1, 22, opt.color[i], 1.);
        up_Pt_hint_error_cl95[i]->GetYaxis()->SetRangeUser(lo_Pt_hint_error_cl95[i]->GetMinimum() * 1.5, up_Pt_hint_error_cl95[i]->GetMaximum() * 1.5);
        up_Pt_hint_error_cl95[i]->Draw("PE");
        Pt_pdf_linear_var_Low2Up[i]->Draw("same");
        hist1D_graphic_opt(lo_Pt_hint_error_cl95[i], kFALSE, 1, 23, opt.color[i], 1.);
        lo_Pt_hint_error_cl95[i]->Draw("PE same");
        Pt_pdf_linear_var_Up2Low[i]->Draw("same");
        C_Pt_linear_deviation->SaveAs(Form("systematic/images/%s_%s_M_%d_%d.pdf", C_Pt_linear_deviation->GetName(), opt.Generator.Data(), opt.Low_Mass, opt.High_Mass));

        TCanvas *C_Pt_MC_modification = canvas_noratio(Form("C_Pt_MC_modification_%s", opt.Name_DiMu_sel[i].Data()));
        C_Pt_MC_modification->SetLogy();
        h_Pt_MC_scaled_Up2Low[i]->Draw("PE");
        Pt_Conf_interval->Draw("E3same");
        h_Pt_MC_scaled_Up2Low[i]->Draw("PEsame");
        h_Pt_MC[i]->Draw("PE same");
        h_Pt_MC_scaled_Low2Up[i]->Draw("PE same");
        pdf_Pt_scaled_Low2Up[i]->Draw("same");
        pdf_Pt_scaled_Up2Low[i]->Draw("same");
        Pt_Conf_interval->SetFillColor(opt.fillcolor[i]);
        TLegend *legend = new TLegend(0.625, 0.625, 0.925, 0.925);
        legend->AddEntry(h_Pt_MC[i]);
        legend->AddEntry(h_Pt_MC_scaled_Up2Low[i]);
        legend->AddEntry(h_Pt_MC_scaled_Low2Up[i]);
        legend->AddEntry(pdf_Pt[i]);
        legend->AddEntry(pdf_Pt_scaled_Up2Low[i]);
        legend->AddEntry(pdf_Pt_scaled_Low2Up[i]);
        legend->AddEntry(Pt_Conf_interval);
        C_Pt_MC_modification->cd();
        legend->DrawClone("same");

        C_Pt_MC_modification->SaveAs(Form("systematic/images/%s_%s_M_%d_%d.pdf", C_Pt_MC_modification->GetName(), opt.Generator.Data(), opt.Low_Mass, opt.High_Mass));

        //----------------------------Mass------------------------------------------------//

        h_Mass_MC[i] = (TH1F *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/h_Mass_from%s", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut, opt.Name_DiMu_sel[i].Data()));
        if (opt.Low_Mass == 4 && opt.High_Mass == 30)
            h_Mass_MC[i]->Rebin(10);
        hist1D_graphic_opt(h_Mass_MC[i], kFALSE, 1, 24, opt.color[i], 1.);
        h_Mass_MC[i]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
        h_Mass_MC[i]->SetTitle("Original");
        pdf_Mass[i] = new TF1(Form("M_pdf_%s", opt.Name_DiMu_sel[i].Data()), FuncMass, opt.Low_Mass, opt.High_Mass, 4);
        pdf_Mass[i]->SetTitle(Form("%s original PDF", opt.Name_DiMu_sel[i].Data()));
        pdf_Mass[i]->SetParameter(3, par3_Mass[i]);
        pdf_Mass[i]->SetParameter(0, B_Mass[i]);
        pdf_Mass[i]->SetParameter(1, n1_Mass[i]);
        pdf_Mass[i]->SetParameter(2, n2_Mass[i]);
        pdf_Mass[i]->SetLineWidth(4);
        pdf_Mass[i]->SetNpx(200);
        pdf_Mass[i]->SetLineColor(opt.fillcolor[i]);

        h_Mass_MC[i]->Fit(pdf_Mass[i], "LR0I");
        h_Mass_MC[i]->Draw("PE");
        pdf_Mass[i]->Draw("same");

        Param_Mass_original[i] = new TH1F(Form("Param_Mass_original_%s", opt.Name_DiMu_sel[i].Data()), Form("Param_Mass_original_%s", opt.Name_DiMu_sel[i].Data()), 4, 0, 4);
        for (Int_t bin = 0; bin < Param_Mass_original[i]->GetNbinsX(); bin++)
        {
            Param_Mass_original[i]->SetBinContent(bin + 1, pdf_Mass[i]->GetParameter(bin));
            Param_Mass_original[i]->SetBinError(bin + 1, pdf_Mass[i]->GetParError(bin));
        }

        fOut_param->cd(saving_dir);
        Param_Mass_original[i]->Write(0, 2, 0);

        M_Conf_interval = (TH1F *)h_Mass_MC[i]->Clone(Form("h_M_CL_95_%s", opt.Name_DiMu_sel[i].Data()));
        M_Conf_interval->Reset();
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(M_Conf_interval, 0.95);
        M_Conf_interval->SetTitle("95% CL");

        up_Mass_hint_error_cl95[i] = Linear_deviation(M_Conf_interval);
        up_Mass_hint_error_cl95[i]->SetName(Form("up_Mass_hint_error_cl95_%s_cl95", opt.Name_DiMu_sel[i].Data()));

        lo_Mass_hint_error_cl95[i] = (TH1F *)up_Mass_hint_error_cl95[i]->Clone(Form("lo_Mass_hint_error_cl95_%s_cl95", opt.Name_DiMu_sel[i].Data()));
        lo_Mass_hint_error_cl95[i]->Reset();
        for (Int_t bin = 0; bin < lo_Mass_hint_error_cl95[i]->GetNbinsX(); bin++)
        {
            lo_Mass_hint_error_cl95[i]->SetBinContent(bin + 1, -(up_Mass_hint_error_cl95[i]->GetBinContent(bin + 1)));
            lo_Mass_hint_error_cl95[i]->SetBinError(bin + 1, (up_Mass_hint_error_cl95[i]->GetBinError(bin + 1)));
        }

        Double_t Mass_x1_Low2Up = opt.Low_Mass;

        Double_t Mass_y1_Low2Up = lo_Mass_hint_error_cl95[i]->GetBinContent(1);

        Double_t Mass_x2_Low2Up = opt.High_Mass;

        Double_t Mass_y2_Low2Up = up_Mass_hint_error_cl95[i]->GetBinContent(up_Mass_hint_error_cl95[i]->GetNbinsX());

        cout << "Mass_x1_Low2Up: " << Mass_x1_Low2Up << "Mass_y1_Low2Up: " << Mass_y1_Low2Up << "Mass_x2_Low2Up: " << Mass_x2_Low2Up << "Mass_y2_Low2Up: " << Mass_y2_Low2Up << endl;

        Double_t Mass_q = Mass_y1_Low2Up - ((Mass_y2_Low2Up - Mass_y1_Low2Up) / (Mass_x2_Low2Up - Mass_x1_Low2Up)) * Mass_x1_Low2Up;

        Double_t Mass_m = (Mass_y2_Low2Up - Mass_y1_Low2Up) / (Mass_x2_Low2Up - Mass_x1_Low2Up);

        //--------------------------Low 2 Up------------------------//

        Mass_pdf_linear_var_Low2Up[i] = new TF1(Form("Mass_pdf_linear_var_Low2Up_from%s", opt.Name_DiMu_sel[i].Data()), "pol1", opt.Low_Mass, opt.High_Mass);
        Mass_pdf_linear_var_Low2Up[i]->FixParameter(0, Mass_q);
        Mass_pdf_linear_var_Low2Up[i]->FixParameter(1, Mass_m);

        h_Mass_MC_scaled_Low2Up[i] = MC_deviation(h_Mass_MC[i], Mass_pdf_linear_var_Low2Up[i]);
        hist1D_graphic_opt(h_Mass_MC_scaled_Low2Up[i], kFALSE, 1, 22, opt.color[i], 1.);
        h_Mass_MC_scaled_Low2Up[i]->SetTitle("Low to Up");
        pdf_Mass_scaled_Low2Up[i] = new TF1(Form("pdf_M_scaled_Low2Up_from%s", opt.Name_DiMu_sel[i].Data()), FuncMass, opt.Low_Mass, opt.High_Mass, 4);
        pdf_Mass_scaled_Low2Up[i]->SetTitle(Form("%s Low2Up mod. PDF", opt.Name_DiMu_sel[i].Data()));
        pdf_Mass_scaled_Low2Up[i]->SetParameter(3, pdf_Mass[i]->GetParameter(3));
        pdf_Mass_scaled_Low2Up[i]->SetParameter(0, pdf_Mass[i]->GetParameter(0));
        pdf_Mass_scaled_Low2Up[i]->SetParameter(1, pdf_Mass[i]->GetParameter(1));
        pdf_Mass_scaled_Low2Up[i]->SetParameter(2, pdf_Mass[i]->GetParameter(2));
        pdf_Mass_scaled_Low2Up[i]->SetLineColor(opt.color[i]);
        pdf_Mass_scaled_Low2Up[i]->SetLineStyle(kDotted);
        pdf_Mass_scaled_Low2Up[i]->SetLineWidth(4);
        pdf_Mass_scaled_Low2Up[i]->SetNpx(200);

        h_Mass_MC_scaled_Low2Up[i]->Fit(pdf_Mass_scaled_Low2Up[i], "LR0I");

        Param_Mass_scaled_Low2Up[i] = new TH1F(Form("Param_Mass_scaled_Low2Up_%s", opt.Name_DiMu_sel[i].Data()), Form("Param_Mass_scaled_Low2Up_%s", opt.Name_DiMu_sel[i].Data()), 4, 0, 4);
        for (Int_t bin = 0; bin < Param_Mass_scaled_Low2Up[i]->GetNbinsX(); bin++)
        {
            Param_Mass_scaled_Low2Up[i]->SetBinContent(bin + 1, pdf_Mass_scaled_Low2Up[i]->GetParameter(bin));
            Param_Mass_scaled_Low2Up[i]->SetBinError(bin + 1, pdf_Mass_scaled_Low2Up[i]->GetParError(bin));
        }

        fOut_param->cd(saving_dir);
        Param_Mass_scaled_Low2Up[i]->Write(0, 2, 0);

        //--------------------------Up 2 Low------------------------//

        Mass_pdf_linear_var_Up2Low[i] = new TF1(Form("Mass_pdf_linear_var_Up2Low_from%s", opt.Name_DiMu_sel[i].Data()), "pol1", opt.Low_Mass, opt.High_Mass);
        Mass_pdf_linear_var_Up2Low[i]->FixParameter(0, -Mass_q);
        Mass_pdf_linear_var_Up2Low[i]->FixParameter(1, -Mass_m);

        h_Mass_MC_scaled_Up2Low[i] = MC_deviation(h_Mass_MC[i], Mass_pdf_linear_var_Up2Low[i]);
        hist1D_graphic_opt(h_Mass_MC_scaled_Up2Low[i], kFALSE, 1, 23, opt.color[i], 1.);
        h_Mass_MC_scaled_Up2Low[i]->SetTitle("Up to Low");

        pdf_Mass_scaled_Up2Low[i] = new TF1(Form("pdf_M_scaled_Up2Low_from%s", opt.Name_DiMu_sel[i].Data()), FuncMass, opt.Low_Mass, opt.High_Mass, 4);
        pdf_Mass_scaled_Up2Low[i]->SetTitle(Form("%s Up2Low mod. PDF", opt.Name_DiMu_sel[i].Data()));
        pdf_Mass_scaled_Up2Low[i]->SetParameter(3, pdf_Mass[i]->GetParameter(3));
        pdf_Mass_scaled_Up2Low[i]->SetParameter(0, pdf_Mass[i]->GetParameter(0));
        pdf_Mass_scaled_Up2Low[i]->SetParameter(1, pdf_Mass[i]->GetParameter(1));
        pdf_Mass_scaled_Up2Low[i]->SetParameter(2, pdf_Mass[i]->GetParameter(2));
        pdf_Mass_scaled_Up2Low[i]->SetLineColor(opt.color[i]);
        pdf_Mass_scaled_Up2Low[i]->SetLineStyle(kDashed);
        pdf_Mass_scaled_Up2Low[i]->SetLineWidth(4);
        pdf_Mass_scaled_Up2Low[i]->SetNpx(200);

        h_Mass_MC_scaled_Up2Low[i]->Fit(pdf_Mass_scaled_Up2Low[i], "LR0I");

        Param_Mass_scaled_Up2Low[i] = new TH1F(Form("Param_Mass_scaled_Up2Low_%s", opt.Name_DiMu_sel[i].Data()), Form("Param_Mass_scaled_Up2Low_%s", opt.Name_DiMu_sel[i].Data()), 4, 0, 4);
        for (Int_t bin = 0; bin < Param_Mass_scaled_Up2Low[i]->GetNbinsX(); bin++)
        {
            Param_Mass_scaled_Up2Low[i]->SetBinContent(bin + 1, pdf_Mass_scaled_Up2Low[i]->GetParameter(bin));
            Param_Mass_scaled_Up2Low[i]->SetBinError(bin + 1, pdf_Mass_scaled_Up2Low[i]->GetParError(bin));
        }

        fOut_param->cd(saving_dir);
        Param_Mass_scaled_Up2Low[i]->Write(0, 2, 0);

        TCanvas *C_M_linear_deviation = canvas_noratio(Form("C_M_linear_deviation_%s", opt.Name_DiMu_sel[i].Data()));
        hist1D_graphic_opt(up_Mass_hint_error_cl95[i], kFALSE, 1, 22, opt.color[i], 1.);
        up_Mass_hint_error_cl95[i]->GetYaxis()->SetRangeUser(lo_Mass_hint_error_cl95[i]->GetMinimum() * 1.5, up_Mass_hint_error_cl95[i]->GetMaximum() * 1.5);
        up_Mass_hint_error_cl95[i]->Draw("PE");
        Mass_pdf_linear_var_Low2Up[i]->Draw("same");
        hist1D_graphic_opt(lo_Mass_hint_error_cl95[i], kFALSE, 1, 23, opt.color[i], 1.);
        lo_Mass_hint_error_cl95[i]->Draw("PE same");
        Mass_pdf_linear_var_Up2Low[i]->Draw("same");
        C_M_linear_deviation->SaveAs(Form("systematic/images/%s_%s_M_%d_%d.pdf", C_M_linear_deviation->GetName(), opt.Generator.Data(), opt.Low_Mass, opt.High_Mass));

        TCanvas *C_Mass_MC_modification = canvas_noratio(Form("C_Mass_MC_modification_%s", opt.Name_DiMu_sel[i].Data()));
        C_Mass_MC_modification->SetLogy();
        h_Mass_MC_scaled_Up2Low[i]->Draw("PE");
        M_Conf_interval->SetFillColor(opt.fillcolor[i]);
        if (opt.Low_Mass == 4 && opt.High_Mass == 30)
        {
            M_Conf_interval->GetXaxis()->SetRangeUser(4, 8.5);
            M_Conf_interval->Draw("E3same");
            M_Conf_interval->GetXaxis()->SetRangeUser(11.5, 30);
            M_Conf_interval->Draw("E3same");

            TF1 *fleft = new TF1(Form("left_pdf_M_scaled_Up2Low_from%s", opt.Name_DiMu_sel[i].Data()), FuncMass, opt.Low_Mass, opt.LowM_cut, 4);
            fleft->SetParameters(pdf_Mass_scaled_Up2Low[i]->GetParameters());
            fleft->SetLineColor(opt.color[i]);
            fleft->SetLineStyle(pdf_Mass_scaled_Up2Low[i]->GetLineStyle());
            fleft->SetLineWidth(pdf_Mass_scaled_Up2Low[i]->GetLineWidth());
            fleft->DrawCopy("same");

            fleft = new TF1(Form("right_pdf_M_scaled_Up2Low_from%s", opt.Name_DiMu_sel[i].Data()), FuncMass, opt.HighM_cut, opt.High_Mass, 4);
            fleft->SetParameters(pdf_Mass_scaled_Up2Low[i]->GetParameters());
            fleft->SetLineColor(opt.color[i]);
            fleft->SetLineStyle(pdf_Mass_scaled_Up2Low[i]->GetLineStyle());
            fleft->SetLineWidth(pdf_Mass_scaled_Up2Low[i]->GetLineWidth());
            fleft->DrawCopy("same");

            fleft = new TF1(Form("left_pdf_Mass_scaled_Low2Up_from%s", opt.Name_DiMu_sel[i].Data()), FuncMass, opt.Low_Mass, opt.LowM_cut, 4);
            fleft->SetParameters(pdf_Mass_scaled_Low2Up[i]->GetParameters());
            fleft->SetLineColor(opt.color[i]);
            fleft->SetLineStyle(pdf_Mass_scaled_Low2Up[i]->GetLineStyle());
            fleft->SetLineWidth(pdf_Mass_scaled_Low2Up[i]->GetLineWidth());
            fleft->DrawCopy("same");

            fleft = new TF1(Form("right_pdf_Mass_scaled_Low2Up_from%s", opt.Name_DiMu_sel[i].Data()), FuncMass, opt.HighM_cut, opt.High_Mass, 4);
            fleft->SetParameters(pdf_Mass_scaled_Low2Up[i]->GetParameters());
            fleft->SetLineColor(opt.color[i]);
            fleft->SetLineStyle(pdf_Mass_scaled_Low2Up[i]->GetLineStyle());
            fleft->SetLineWidth(pdf_Mass_scaled_Low2Up[i]->GetLineWidth());
            fleft->DrawCopy("same");

            fleft = new TF1(Form("left_pdf_Mass_original_from%s", opt.Name_DiMu_sel[i].Data()), FuncMass, opt.Low_Mass, opt.LowM_cut, 4);
            fleft->SetParameters(pdf_Mass[i]->GetParameters());
            fleft->SetLineWidth(pdf_Mass[i]->GetLineWidth());
            fleft->SetLineColor(opt.color[i]);
            fleft->DrawCopy("same");

            fleft = new TF1(Form("right_pdf_Mass_original_from%s", opt.Name_DiMu_sel[i].Data()), FuncMass, opt.HighM_cut, opt.High_Mass, 4);
            fleft->SetParameters(pdf_Mass[i]->GetParameters());
            fleft->SetLineWidth(pdf_Mass[i]->GetLineWidth());
            fleft->SetLineColor(opt.color[i]);
            fleft->DrawCopy("same");
        }
        else
        {
            M_Conf_interval->Draw("E3same");
            pdf_Mass_scaled_Low2Up[i]->Draw("same");
            pdf_Mass_scaled_Up2Low[i]->Draw("same");
        }

        h_Mass_MC_scaled_Up2Low[i]->Draw("PEsame");
        h_Mass_MC[i]->Draw("PE same");
        h_Mass_MC_scaled_Low2Up[i]->Draw("PE same");

        legend = new TLegend(0.625, 0.625, 0.925, 0.925);
        legend->AddEntry(h_Mass_MC[i]);
        legend->AddEntry(h_Mass_MC_scaled_Up2Low[i]);
        legend->AddEntry(h_Mass_MC_scaled_Low2Up[i]);
        legend->AddEntry(pdf_Mass[i]);
        legend->AddEntry(pdf_Mass_scaled_Up2Low[i]);
        legend->AddEntry(pdf_Mass_scaled_Low2Up[i]);
        legend->AddEntry(M_Conf_interval);
        C_Mass_MC_modification->cd();
        legend->DrawClone("same");
        C_Mass_MC_modification->SaveAs(Form("systematic/images/%s_%s_M_%d_%d.pdf", C_Mass_MC_modification->GetName(), opt.Generator.Data(), opt.Low_Mass, opt.High_Mass));
    }
}

void comparison_new_old()
{
    gStyle->SetOptTitle(0);
    info opt;

    TFile *old_syst = new TFile("test_old_systematic/signal_extraction_systematic_pdf.root", "READ");

    TFile *new_syst = new TFile("systematic/template_modification_PYTHIA.root", "READ");
    new_syst->ls();

    for (Int_t i_var = 0; i_var < N_variation; i_var++)
    {
        for (Int_t i_dimu = 0; i_dimu < N_DiMu_sel-1; i_dimu++)
        {
            TH1F *Pt_old_syst_hist = (TH1F *)old_syst->Get(Form("Param_Pt_%s_from%s", opt.Name_variation[i_var].Data(), opt.Name_DiMu_sel[i_dimu].Data()));
            hist1D_graphic_opt(Pt_old_syst_hist, kFALSE, 1, 20, kAzure + 7, 1.);
            Pt_old_syst_hist->SetTitle("for the prelimary");
            Pt_old_syst_hist->GetXaxis()->SetRangeUser(0, 3);
            Pt_old_syst_hist->GetXaxis()->SetLabelSize(0.045);
            for (Int_t i_bin = 1; i_bin <= Pt_old_syst_hist->GetNbinsX() - 1; i_bin++)
                Pt_old_syst_hist->GetXaxis()->SetBinLabel(i_bin, opt.Bin_Label[i_bin - 1]);
            TH1F *Pt_new_syst_hist = (TH1F *)new_syst->Get(Form("param_syst_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/Param_Pt_%s_%s", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut, opt.Name_variation[i_var].Data(), opt.Name_DiMu_sel[i_dimu].Data()));
            hist1D_graphic_opt(Pt_new_syst_hist, kFALSE, 1, 20, kRed + 2, 1.);
            Pt_new_syst_hist->SetTitle("Now");
            Pt_new_syst_hist->GetXaxis()->SetRangeUser(0, 3);
            Pt_new_syst_hist->GetXaxis()->SetLabelSize(0.045);
            for (Int_t i_bin = 1; i_bin <= Pt_new_syst_hist->GetNbinsX() - 1; i_bin++)
                Pt_new_syst_hist->GetXaxis()->SetBinLabel(i_bin, opt.Bin_Label[i_bin - 1]);

            TCanvas *C_pt_old_new = canvas_noratio(Form("Pt_comparison_%s_%s", opt.Name_DiMu_sel[i_dimu].Data(), opt.Name_variation[i_var].Data()));
            if (Pt_old_syst_hist->GetMaximum() > Pt_new_syst_hist->GetMaximum())
            {
                Pt_old_syst_hist->SetMaximum(Pt_old_syst_hist->GetMaximum() + 2);
                Pt_old_syst_hist->Draw("hist E TEXT0");
                Pt_new_syst_hist->Draw("hist E text0 same");
            }
            else
            {
                Pt_new_syst_hist->SetMaximum(Pt_new_syst_hist->GetMaximum() + 2);
                Pt_new_syst_hist->Draw("hist E text0");
                Pt_old_syst_hist->Draw("hist E TEXT0 same");
            }
            TLegend *Pt_legend = (TLegend *)gPad->BuildLegend();
            Legend_settings(Pt_legend, 0.2, 0.5, 0.725, 0.885, Form("Parameters %s %s", opt.Name_DiMu_sel[i_dimu].Data(), opt.Name_variation[i_var].Data()));
            // C_pt_old_new->SaveAs(Form("test_old_systematic/comparison_old_new/%s.pdf", C_pt_old_new->GetName()));

            TH1F *M_old_syst_hist = (TH1F *)old_syst->Get(Form("Param_M_%s_from%s", opt.Name_variation[i_var].Data(), opt.Name_DiMu_sel[i_dimu].Data()));
            hist1D_graphic_opt(M_old_syst_hist, kFALSE, 1, 20, kAzure + 7, 1.);
            M_old_syst_hist->SetTitle("for the prelimary");
            M_old_syst_hist->GetXaxis()->SetRangeUser(0, 3);
            M_old_syst_hist->GetXaxis()->SetLabelSize(0.045);
            for (Int_t i_bin = 1; i_bin <= Pt_old_syst_hist->GetNbinsX() - 1; i_bin++)
                M_old_syst_hist->GetXaxis()->SetBinLabel(i_bin, opt.Bin_Label[i_bin - 1]);

            TH1F *M_new_syst_hist = (TH1F *)new_syst->Get(Form("param_syst_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/Param_Mass_%s_%s", opt.Low_Mass, opt.High_Mass, opt.Low_Pt, opt.High_Pt, opt.LowM_cut, opt.HighM_cut, opt.Name_variation[i_var].Data(), opt.Name_DiMu_sel[i_dimu].Data()));
            hist1D_graphic_opt(M_new_syst_hist, kFALSE, 1, 20, kRed + 2, 1.);
            M_new_syst_hist->SetTitle("Now");
            M_new_syst_hist->GetXaxis()->SetRangeUser(0, 3);
            M_new_syst_hist->GetXaxis()->SetLabelSize(0.045);
            for (Int_t i_bin = 1; i_bin <= Pt_new_syst_hist->GetNbinsX() - 1; i_bin++)
                M_new_syst_hist->GetXaxis()->SetBinLabel(i_bin, opt.Bin_Label[i_bin - 1]);

            TCanvas *C_M_old_new = canvas_noratio(Form("M_comparison_%s_%s", opt.Name_DiMu_sel[i_dimu].Data(), opt.Name_variation[i_var].Data()));
            if (M_old_syst_hist->GetMaximum() > M_new_syst_hist->GetMaximum())
            {
                Pt_old_syst_hist->SetMaximum(Pt_old_syst_hist->GetMaximum() + 2);
                M_old_syst_hist->Draw("hist E TEXT0");
                M_new_syst_hist->Draw("hist E text0 same");
            }
            else
            {
                Pt_new_syst_hist->SetMaximum(Pt_new_syst_hist->GetMaximum() + 2);
                M_new_syst_hist->Draw("hist E text0");
                M_old_syst_hist->Draw("hist E TEXT0 same");
            }

            TLegend *M_legend = (TLegend *)gPad->BuildLegend();
            Legend_settings(M_legend, 0.2, 0.5, 0.725, 0.885, Form("Parameters %s %s", opt.Name_DiMu_sel[i_dimu].Data(), opt.Name_variation[i_var].Data()));
            // C_M_old_new->SaveAs(Form("test_old_systematic/comparison_old_new/%s.pdf", C_M_old_new->GetName()));
        }
    }
}

double FuncMass(double *x, double *par)
{
    info opt;
    if (opt.Low_Mass == 4 && opt.High_Mass == 30)
    {
        if (x[0] > opt.LowM_cut && x[0] < opt.HighM_cut)
        {
            TF1::RejectPoint();
            return 0;
        }
        else
            return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
    }
    else
        return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
}

double FuncPt(double *x, double *par)
{
    return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
}