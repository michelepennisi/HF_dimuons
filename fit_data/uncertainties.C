#include "/home/michele_pennisi/cernbox/HF_dimuons/common_include.h"

const Int_t N_DiMu_sel = 3;
const Int_t N_variation = 3;

struct info
{
    TString old_syst_file = "";
    TString new_syst_file = "";

    TString Name_DiMu_sel[N_DiMu_sel] = {"Charm", "Beauty", "Mixed"};
    TString Name_variation[N_variation] = {"original", "scaled_Low2Up", "scaled_Up2Low"};

    TString Bin_Label[3]={"B","n1","n2"};
     
};

void uncertainties()
{
}

void comparison_new_old()
{
    gStyle->SetOptTitle(0);
    info opt;

    TFile *old_syst = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/test_old_systematic/signal_extraction_systematic_pdf.root", "READ");

    TFile *new_syst = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/test_old_systematic/test_old_signal_extraction_systematic_pdf.root", "READ");

    for (Int_t i_var = 0; i_var < N_variation; i_var++)
    {
        for (Int_t i_dimu = 0; i_dimu < N_DiMu_sel; i_dimu++)
        {
            TH1F *Pt_old_syst_hist = (TH1F *)old_syst->Get(Form("Param_Pt_%s_from%s", opt.Name_variation[i_var].Data(), opt.Name_DiMu_sel[i_dimu].Data()));
            hist1D_graphic_opt(Pt_old_syst_hist, kFALSE, 1, 20, kAzure + 7, 1.);
            Pt_old_syst_hist->SetTitle("for the prelimary");
            Pt_old_syst_hist->GetXaxis()->SetRangeUser(0, 3);
            Pt_old_syst_hist->GetXaxis()->SetLabelSize(0.045);
            for (Int_t i_bin=1;i_bin<=Pt_old_syst_hist->GetNbinsX()-1;i_bin++) Pt_old_syst_hist->GetXaxis()->SetBinLabel(i_bin,opt.Bin_Label[i_bin-1]);

            TH1F *Pt_new_syst_hist = (TH1F *)new_syst->Get(Form("Param_Pt_%s_from%s", opt.Name_variation[i_var].Data(), opt.Name_DiMu_sel[i_dimu].Data()));
            hist1D_graphic_opt(Pt_new_syst_hist, kFALSE, 1, 20, kRed + 2, 1.);
            Pt_new_syst_hist->SetTitle("Now");
            Pt_new_syst_hist->GetXaxis()->SetRangeUser(0, 3);
            Pt_new_syst_hist->GetXaxis()->SetLabelSize(0.045);
            for (Int_t i_bin=1;i_bin<=Pt_new_syst_hist->GetNbinsX()-1;i_bin++) Pt_new_syst_hist->GetXaxis()->SetBinLabel(i_bin,opt.Bin_Label[i_bin-1]);

            TCanvas *C_pt_old_new = canvas_noratio(Form("Pt_comparison_%s_%s", opt.Name_DiMu_sel[i_dimu].Data(), opt.Name_variation[i_var].Data()));
            if (Pt_old_syst_hist->GetMaximum() > Pt_new_syst_hist->GetMaximum())
            {
                Pt_old_syst_hist->SetMaximum(Pt_old_syst_hist->GetMaximum()+2);
                Pt_old_syst_hist->Draw("hist E TEXT0");
                Pt_new_syst_hist->Draw("hist E text0 same");
            }
            else
            {
                Pt_new_syst_hist->SetMaximum(Pt_new_syst_hist->GetMaximum()+2);
                Pt_new_syst_hist->Draw("hist E text0");
                Pt_old_syst_hist->Draw("hist E TEXT0 same");
            }
            TLegend *Pt_legend = (TLegend *)gPad->BuildLegend();
            Legend_settings(Pt_legend, 0.2, 0.5, 0.725, 0.885, Form("Parameters %s %s", opt.Name_DiMu_sel[i_dimu].Data(), opt.Name_variation[i_var].Data()));
            C_pt_old_new->SaveAs(Form("test_old_systematic/comparison_old_new/%s.pdf",C_pt_old_new->GetName()));

            TH1F *M_old_syst_hist = (TH1F *)old_syst->Get(Form("Param_M_%s_from%s", opt.Name_variation[i_var].Data(), opt.Name_DiMu_sel[i_dimu].Data()));
            hist1D_graphic_opt(M_old_syst_hist, kFALSE, 1, 20, kAzure + 7, 1.);
            M_old_syst_hist->SetTitle("for the prelimary");
            M_old_syst_hist->GetXaxis()->SetRangeUser(0, 3);
            M_old_syst_hist->GetXaxis()->SetLabelSize(0.045);
            for (Int_t i_bin=1;i_bin<=Pt_old_syst_hist->GetNbinsX()-1;i_bin++) M_old_syst_hist->GetXaxis()->SetBinLabel(i_bin,opt.Bin_Label[i_bin-1]);

            TH1F *M_new_syst_hist = (TH1F *)new_syst->Get(Form("Param_M_%s_from%s", opt.Name_variation[i_var].Data(), opt.Name_DiMu_sel[i_dimu].Data()));
            hist1D_graphic_opt(M_new_syst_hist, kFALSE, 1, 20, kRed + 2, 1.);
            M_new_syst_hist->SetTitle("Now");
            M_new_syst_hist->GetXaxis()->SetRangeUser(0, 3);
            M_new_syst_hist->GetXaxis()->SetLabelSize(0.045);
            for (Int_t i_bin=1;i_bin<=Pt_new_syst_hist->GetNbinsX()-1;i_bin++) M_new_syst_hist->GetXaxis()->SetBinLabel(i_bin,opt.Bin_Label[i_bin-1]);

            TCanvas *C_M_old_new = canvas_noratio(Form("M_comparison_%s_%s", opt.Name_DiMu_sel[i_dimu].Data(), opt.Name_variation[i_var].Data()));
            if (M_old_syst_hist->GetMaximum() > M_new_syst_hist->GetMaximum())
            {
                Pt_old_syst_hist->SetMaximum(Pt_old_syst_hist->GetMaximum()+2);
                M_old_syst_hist->Draw("hist E TEXT0");
                M_new_syst_hist->Draw("hist E text0 same");
            }
            else {
                Pt_new_syst_hist->SetMaximum(Pt_new_syst_hist->GetMaximum()+2);
                M_new_syst_hist->Draw("hist E text0");
                M_old_syst_hist->Draw("hist E TEXT0 same");
            }

            TLegend *M_legend = (TLegend *)gPad->BuildLegend();
            Legend_settings(M_legend, 0.2, 0.5, 0.725, 0.885, Form("Parameters %s %s", opt.Name_DiMu_sel[i_dimu].Data(), opt.Name_variation[i_var].Data()));
            C_M_old_new->SaveAs(Form("test_old_systematic/comparison_old_new/%s.pdf",C_M_old_new->GetName()));
        }
    }
}