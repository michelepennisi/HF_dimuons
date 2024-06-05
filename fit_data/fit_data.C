#include "../common_include.h"
#include <vector>
#include "RooChi2Var.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
using namespace RooStats;

const Int_t n_DiMuSelection = 6;

struct opt
{
    TString path_to_file = TString::Format("%s", getenv("universal_PATH"));
    Int_t Mass_Binning = 20;
    Int_t Low_Mass = 4;
    Int_t High_Mass = 30;
    Double_t LowM_cut = 8.;
    Double_t HighM_cut = 11.;
    Int_t Pt_Binning = 20;
    Int_t Low_Pt = 0;
    Int_t High_Pt = 30;

    Double_t HF_Mixed_fraction = 3.0;
    Double_t LF_HF_Mixed_fraction = 14.9;
    Double_t LF_fraction = 2.9;
    // Double_t LF_HF_Mixed_fraction = 0.;

    TString Generator = "POWHEG";
    TString stat_MC = "full_stat";
    TString stat_Data = "LHC18p";
    // TString DY = "withDY";
    TString LF_HF = "withLF_HF_Mixed_Bkg";
    // TString LF_HF = "noLF_HF_LHC23i2";
    // TString LF_HF = "withLF_HF_LHC23i1";
    // TString LF_HF = "withLF_HF_LHC23i2";
    TString DY = "withDY";
    TString LF = "noLF";
    // TString LF_HF = "noLF_HF_LHC23i1";
    Bool_t Plot_Likehood = kFALSE;
    TString systematic = "";
    // TString systematic = "scaled_Low2Up";
    // TString systematic = "scaled_Up2Low";

    TString Name_DimuSel[n_DiMuSelection] = {"Charm", "Beauty", "DY", "HF_Mixed", "LF_HF_Mixed", "LF"};
    Color_t color[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kOrange + 7, kAzure + 9, kRed - 2, kBlack};
    TString info_label[n_DiMuSelection] = {"#leftarrow c", "#leftarrow b", "#leftarrow DY", "#leftarrow c,b", "#leftarrow LF,HF", "#leftarrow LF"};
};
double FuncMass(double *x, double *par)
{
    opt info;
    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        if (x[0] > info.LowM_cut && x[0] < info.HighM_cut)
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

double FuncPtMass(double *x, double *par)
{
    return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
}
TCanvas *printMC_ratio(TString name, TString Title, RooPlot *frame, RooArgSet *param, TH1 *data, TF1 *pdf, Color_t color, Int_t minx = 0, Int_t max_x = 30);
TCanvas *printRooPlot_ratio(RooPlot *frame, Bool_t norm, RooFitResult *r, Int_t choice, TString roohist_name, TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x, Double_t N_HFMixed, Double_t N_LF_HFMixed);
TCanvas *printMC_ratio_BINNED(TString name, TString Title, TH1 *data, TF1 *pdf, Color_t color);
Double_t Lowy = 0.0025;

void HF_mixing_POWHEG();
void Mixing_background();
void unbinned_fit_data_sample_singleregion();
void shape_comparison();
void plotting_fit();
void syst_plot();

void fit_data()
{
    // unbinned_fit_data_sample_singleregion();
    // plotting_fit();
    // Mixing_background();
    HF_mixing_POWHEG();
    // syst_plot();
}

void shape_comparison()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    const Int_t n_DiMuSelection = 4;

    TString name_DiMu_Sel[n_DiMuSelection];
    name_DiMu_Sel[0] = TString("Charm");
    name_DiMu_Sel[1] = TString("Beauty");
    name_DiMu_Sel[2] = TString("HF_Mixed");
    name_DiMu_Sel[3] = TString("DY");

    TFile *fIn_MC[n_DiMuSelection];
    TString MC_file[n_DiMuSelection];
    TTree *m_tree_MC[n_DiMuSelection];
    TH1F *m_shape_4_30[n_DiMuSelection];
    TH1F *pt_shape_4_30[n_DiMuSelection];

    TH1F *m_shape_11_30[n_DiMuSelection];
    TH1F *pt_shape_11_30[n_DiMuSelection];

    MC_file[0] = TString("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i1/Version3_AliAOD/save_output/LHC23i1_MC_output_Tree_merged.root");
    MC_file[1] = TString("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version3_AliAOD/save_output/LHC23i2_MC_output_Tree_merged.root");
    MC_file[2] = TString("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version3_AliAOD/save_output/LHC23i2_MC_output_Tree_merged.root");
    MC_file[3] = TString("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Powheg_Sim/powheg_DY_mass_3_35/Version1/Analysis_MCsim/powheg_DY_mass_3_35_Analysis_MCsim_Tree_merged.root");

    Color_t color[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kAzure + 9, kOrange + 7};

    for (Int_t i_DiMuSel = 0; i_DiMuSel < n_DiMuSelection; i_DiMuSel++)
    {
        fIn_MC[i_DiMuSel] = new TFile(MC_file[i_DiMuSel].Data(), "READ");
        if (i_DiMuSel < n_DiMuSelection - 1)
            m_tree_MC[i_DiMuSel] = (TTree *)fIn_MC[i_DiMuSel]->Get(Form("DiMuon_Rec_PowhegOnly_%s", name_DiMu_Sel[i_DiMuSel].Data()));
        else
            m_tree_MC[i_DiMuSel] = (TTree *)fIn_MC[i_DiMuSel]->Get("rec_tree_muDY");

        m_shape_4_30[i_DiMuSel] = new TH1F(Form("m_shape_4_30_%s", name_DiMu_Sel[i_DiMuSel].Data()), Form("#mu#mu from %s", name_DiMu_Sel[i_DiMuSel].Data()), 26, 4, 30);
        m_shape_11_30[i_DiMuSel] = new TH1F(Form("m_shape_11_30_%s", name_DiMu_Sel[i_DiMuSel].Data()), Form("#mu#mu from %s", name_DiMu_Sel[i_DiMuSel].Data()), 10, 11, 30);

        pt_shape_4_30[i_DiMuSel] = new TH1F(Form("pt_shape_4_30_%s", name_DiMu_Sel[i_DiMuSel].Data()), Form("#mu#mu from %s", name_DiMu_Sel[i_DiMuSel].Data()), 30, 0, 30);
        pt_shape_11_30[i_DiMuSel] = new TH1F(Form("pt_shape_11_30_%s", name_DiMu_Sel[i_DiMuSel].Data()), Form("#mu#mu from %s", name_DiMu_Sel[i_DiMuSel].Data()), 15, 0, 30);

        m_tree_MC[i_DiMuSel]->Draw(TString::Format("m>>%s", m_shape_4_30[i_DiMuSel]->GetName()), "(m > 4 && m< 30) && (pt > 0 && pt< 30)", "goff");
        m_shape_4_30[i_DiMuSel]->Scale(1. / (m_shape_4_30[i_DiMuSel]->GetEntries()));
        m_shape_4_30[i_DiMuSel]->SetMarkerStyle(20);
        m_shape_4_30[i_DiMuSel]->SetMarkerColor(color[i_DiMuSel]);
        m_shape_4_30[i_DiMuSel]->SetLineColor(color[i_DiMuSel]);
        m_shape_4_30[i_DiMuSel]->SetLineWidth(2);

        m_tree_MC[i_DiMuSel]->Draw(TString::Format("m>>%s", m_shape_11_30[i_DiMuSel]->GetName()), "(m > 11 && m< 30) && (pt > 0 && pt< 30)", "goff");
        m_shape_11_30[i_DiMuSel]->Scale(1. / (m_shape_11_30[i_DiMuSel]->GetEntries()));
        m_shape_11_30[i_DiMuSel]->SetMarkerStyle(20);
        m_shape_11_30[i_DiMuSel]->SetMarkerColor(color[i_DiMuSel]);
        m_shape_11_30[i_DiMuSel]->SetLineColor(color[i_DiMuSel]);
        m_shape_11_30[i_DiMuSel]->SetLineWidth(2);

        m_tree_MC[i_DiMuSel]->Draw(TString::Format("pt>>%s", pt_shape_4_30[i_DiMuSel]->GetName()), "(m > 4 && m< 30) && (pt > 0 && pt< 30)", "goff");
        pt_shape_4_30[i_DiMuSel]->Scale(1. / (pt_shape_4_30[i_DiMuSel]->GetEntries()));
        pt_shape_4_30[i_DiMuSel]->SetMarkerStyle(20);
        pt_shape_4_30[i_DiMuSel]->SetMarkerColor(color[i_DiMuSel]);
        pt_shape_4_30[i_DiMuSel]->SetLineColor(color[i_DiMuSel]);
        pt_shape_4_30[i_DiMuSel]->SetLineWidth(2);

        m_tree_MC[i_DiMuSel]->Draw(TString::Format("pt>>%s", pt_shape_11_30[i_DiMuSel]->GetName()), "(m > 11 && m< 30) && (pt > 0 && pt< 30)", "goff");
        pt_shape_11_30[i_DiMuSel]->Scale(1. / (pt_shape_11_30[i_DiMuSel]->GetEntries()));
        pt_shape_11_30[i_DiMuSel]->SetMarkerStyle(20);
        pt_shape_11_30[i_DiMuSel]->SetMarkerColor(color[i_DiMuSel]);
        pt_shape_11_30[i_DiMuSel]->SetLineColor(color[i_DiMuSel]);
        pt_shape_11_30[i_DiMuSel]->SetLineWidth(2);
    }

    TH1F *m_shape_4_30_sum = (TH1F *)m_shape_4_30[0]->Clone("m_shape_4_30_sum");
    m_shape_4_30_sum->Add(m_shape_4_30[1]);
    m_shape_4_30_sum->Add(m_shape_4_30[2]);

    TH1F *m_shape_4_30_DY_HF_ratio = (TH1F *)m_shape_4_30[3]->Clone("m_shape_4_30_DY_HF_ratio");
    m_shape_4_30_DY_HF_ratio->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    m_shape_4_30_DY_HF_ratio->GetYaxis()->SetTitle("DY over HF sum");
    m_shape_4_30_DY_HF_ratio->GetYaxis()->CenterTitle();
    m_shape_4_30_DY_HF_ratio->GetYaxis()->SetTitleSize(m_shape_4_30_DY_HF_ratio->GetXaxis()->GetTitleSize());
    m_shape_4_30_DY_HF_ratio->Divide(m_shape_4_30_sum);

    TH1F *m_shape_11_30_sum = (TH1F *)m_shape_11_30[0]->Clone("m_shape_11_30_sum");
    m_shape_11_30_sum->Add(m_shape_11_30[1]);
    m_shape_11_30_sum->Add(m_shape_11_30[2]);

    TH1F *m_shape_11_30_DY_HF_ratio = (TH1F *)m_shape_11_30[3]->Clone("m_shape_11_30_DY_HF_ratio");
    m_shape_11_30_DY_HF_ratio->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    m_shape_11_30_DY_HF_ratio->GetYaxis()->SetTitle("DY over HF sum");
    m_shape_11_30_DY_HF_ratio->GetYaxis()->CenterTitle();
    m_shape_11_30_DY_HF_ratio->GetYaxis()->SetTitleSize(m_shape_11_30_DY_HF_ratio->GetXaxis()->GetTitleSize());
    m_shape_11_30_DY_HF_ratio->Divide(m_shape_11_30_sum);

    TH1F *pt_shape_4_30_sum = (TH1F *)pt_shape_4_30[0]->Clone("pt_shape_4_30_sum");
    pt_shape_4_30_sum->Add(pt_shape_4_30[1]);
    pt_shape_4_30_sum->Add(pt_shape_4_30[2]);

    TH1F *pt_shape_4_30_DY_HF_ratio = (TH1F *)pt_shape_4_30[3]->Clone("pt_shape_4_30_DY_HF_ratio");
    pt_shape_4_30_DY_HF_ratio->GetXaxis()->SetTitle("#it{p}_{T, #mu} (GeV/#it{c})");
    pt_shape_4_30_DY_HF_ratio->GetYaxis()->SetTitle("DY over HF sum");
    pt_shape_4_30_DY_HF_ratio->GetYaxis()->CenterTitle();
    pt_shape_4_30_DY_HF_ratio->GetYaxis()->SetTitleSize(pt_shape_4_30_DY_HF_ratio->GetXaxis()->GetTitleSize());
    pt_shape_4_30_DY_HF_ratio->Divide(pt_shape_4_30_sum);

    TH1F *pt_shape_11_30_sum = (TH1F *)pt_shape_11_30[0]->Clone("pt_shape_11_30_sum");
    pt_shape_11_30_sum->Add(pt_shape_11_30[1]);
    pt_shape_11_30_sum->Add(pt_shape_11_30[2]);

    TH1F *pt_shape_11_30_DY_HF_ratio = (TH1F *)pt_shape_11_30[3]->Clone("pt_shape_11_30_DY_HF_ratio");
    pt_shape_11_30_DY_HF_ratio->GetXaxis()->SetTitle("#it{p}_{T, #mu} (GeV/#it{c})");
    pt_shape_11_30_DY_HF_ratio->GetYaxis()->SetTitle("DY over HF sum");
    pt_shape_11_30_DY_HF_ratio->GetYaxis()->CenterTitle();
    pt_shape_11_30_DY_HF_ratio->GetYaxis()->SetTitleSize(pt_shape_11_30_DY_HF_ratio->GetXaxis()->GetTitleSize());
    pt_shape_11_30_DY_HF_ratio->Divide(pt_shape_11_30_sum);

    TCanvas *canvas_m_shape_4_30 = three_histo_one_ratio(m_shape_4_30[0], m_shape_4_30[1], m_shape_4_30[3], m_shape_4_30_DY_HF_ratio, "canvas_m_shape_4_30", "#it{p}_{T, #mu} < 30 GeV/#it{c}", kTRUE, kTRUE);
    TPad *prov = (TPad *)canvas_m_shape_4_30->GetListOfPrimitives()->FindObject("pad1");
    TLegend *leg = (TLegend *)prov->FindObject("Legend");
    prov->cd();
    m_shape_4_30[2]->Draw("PESAME");
    leg->AddEntry(m_shape_4_30[2]);
    leg->SetY1(0.125);
    leg->SetY2(0.375);
    prov->Update();
    canvas_m_shape_4_30->SaveAs(Form("images/%s.png", canvas_m_shape_4_30->GetName()));

    TCanvas *canvas_pt_shape_4_30 = three_histo_one_ratio(pt_shape_4_30[0], pt_shape_4_30[1], pt_shape_4_30[3], pt_shape_4_30_DY_HF_ratio, "canvas_pt_shape_4_30", "4 < #it{m}_{#mu#mu} < 30 GeV/#it{c}^{2}", kTRUE, kTRUE);

    prov = (TPad *)canvas_pt_shape_4_30->GetListOfPrimitives()->FindObject("pad1");
    leg = (TLegend *)prov->FindObject("Legend");
    prov->cd();
    pt_shape_4_30[2]->Draw("PESAME");
    leg->AddEntry(m_shape_4_30[2]);
    leg->SetY1(0.125);
    leg->SetY2(0.375);
    canvas_pt_shape_4_30->SaveAs(Form("images/%s.png", canvas_pt_shape_4_30->GetName()));
    prov->Update();

    TCanvas *canvas_m_shape_11_30 = three_histo_one_ratio(m_shape_11_30[0], m_shape_11_30[1], m_shape_11_30[3], m_shape_11_30_DY_HF_ratio, "canvas_m_shape_11_30", "#it{p}_{T, #mu} < 30 GeV/#it{c}", kTRUE, kTRUE);
    prov = (TPad *)canvas_m_shape_11_30->GetListOfPrimitives()->FindObject("pad1");
    leg = (TLegend *)prov->FindObject("Legend");
    prov->cd();
    m_shape_11_30[2]->Draw("PESAME");
    leg->AddEntry(m_shape_4_30[2]);
    leg->SetY1(0.125);
    leg->SetY2(0.375);
    prov->Update();
    canvas_m_shape_11_30->SaveAs(Form("images/%s.png", canvas_m_shape_11_30->GetName()));

    TCanvas *canvas_pt_shape_11_30 = three_histo_one_ratio(pt_shape_11_30[0], pt_shape_11_30[1], pt_shape_11_30[3], pt_shape_11_30_DY_HF_ratio, "canvas_pt_shape_11_30", "11 < #it{m}_{#mu#mu} < 30 GeV/#it{c}^{2}", kTRUE, kTRUE);
    prov = (TPad *)canvas_pt_shape_11_30->GetListOfPrimitives()->FindObject("pad1");
    leg = (TLegend *)prov->FindObject("Legend");
    prov->cd();
    pt_shape_11_30[2]->Draw("PESAME");
    leg->AddEntry(m_shape_4_30[2]);
    leg->SetY1(0.125);
    leg->SetY2(0.375);
    prov->Update();
    canvas_pt_shape_11_30->SaveAs(Form("images/%s.png", canvas_pt_shape_11_30->GetName()));
}

void cross_section()
{
    opt info;
    gROOT->ProcessLineSync(Form(".x %s/HF_dimuons/PtMassExpPdf.cxx+", info.path_to_file.Data()));

    TFile *fIn_data;
    if (info.stat_Data.Contains("Run2"))
        fIn_data = new TFile(Form("%s/HF_dimuons/data_analysis/Tree_MassPt_MassCut4_Run2.root", info.path_to_file.Data()), "READ");
    else
        fIn_data = new TFile(Form("%s/HF_dimuons/data_analysis/3_11_2022/HistResults_merged.root", info.path_to_file.Data()), "READ");
    fIn_data->ls();
    TH1D *fhNEv = (TH1D *)fIn_data->Get("fhNEv");
    fhNEv->Draw();
    fhNEv->SetDirectory(0);
    fIn_data->Close();
    Double_t NEv_MB_Data = (2384.73 * fhNEv->GetBinContent(3));

    TFile *fIn = fIn = new TFile(Form("%s/HF_dimuons/fit_data/results/fit_result/Result_Fit_Data_%s_MC_%s_%s_%s_%s.root", info.path_to_file.Data(), info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data(), info.LF.Data(), info.DY.Data()), "READ");
    cout << Form("%s/HF_dimuons/fit_data/results/fit_result/Result_Fit_Data_%s_MC_%s_%s_%s_%s.root", info.path_to_file.Data(), info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data(), info.LF.Data(), info.DY.Data()) << endl;
    fIn->ls();

    TH1D *r = (TH1D *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/hist_fit_result_%s_%s", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut, info.DY.Data(), info.LF_HF.Data()));
    r->Draw("hist");
    r->SetDirectory(0);
    fIn->Close();
    const Int_t n_obs = 3;
    TString NameDiMu_Sel[n_obs] = {"Charm", "Beauty", "DY"};

    // TString NameMC_file_x_CS_calc[3] = {"/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_charm_nocut_Version_5_AliAOD_withHF_Q/powheg_charm_nocut_Version_5_AliAOD_withHF_Q_MC_output_Hist_294154.root", "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_beauty_nocut_Version_5_AliAOD_withHF_Q/powheg_beauty_nocut_Version_5_AliAOD_withHF_Q_MC_output_Hist_294154.root", "~/cernbox/HF_dimuons/fit_data/ingredient_cs_powheg/DY_cs.root"};
    // TString NameMC_Hist_x_CS_calc[3] = {"HF_quarks/Powheg/h_NCharm_event_PowhegOnly", "HF_quarks/Powheg/h_NBeauty_event_PowhegOnly", "h_NDY_event"};
    // TString NameMC_Hist_x_CS_calc_fwd[3] = {"HF_quarks/Powheg/h_NCharm_event_fwd_PowhegOnly", "HF_quarks/Powheg/h_NBeauty_event_fwd_PowhegOnly", "h_NDY_event_fwd"};
    TString NameMC_file[n_obs];
    TString NameMC_file_x_CS_calc[n_obs];
    TString NameMC_Hist_x_CS_calc[n_obs];
    TString NameMC_Hist_x_CS_calc_fwd[n_obs];
    Double_t fNorm_MC[n_obs];
    Double_t Powheg_CS[n_obs];

    // DY quantities are the same for PYTHIA or POWHEG fits
    NameMC_file[2] = TString::Format("%s/output_HF_dimuons/mc_analysis_output/LHC18p_DY/Version_2_AliAOD/LHC18p_DY_MC_output_Tree_merged.root", info.path_to_file.Data());
    NameMC_file_x_CS_calc[2] = TString::Format("%s/HF_dimuons/fit_data/ingredient_cs_powheg/DY_cs.root", info.path_to_file.Data());
    NameMC_Hist_x_CS_calc[2] = "h_NDY_event";
    NameMC_Hist_x_CS_calc_fwd[2] = "h_NDY_event_fwd";
    fNorm_MC[2] = 1. * (56.42 / 4.6e-05);
    Powheg_CS[2] = (4.6e-05);

    if (info.Generator.Contains("POWHEG"))
    {
        Powheg_CS[0] = (5.0);
        Powheg_CS[1] = (0.5);

        fNorm_MC[0] = {30.5 * (56.42 / 5.0)};
        fNorm_MC[1] = {15.0 * (56.42 / 0.5)};

        NameMC_file[0].Form("%s/output_HF_dimuons/mc_analysis_output/LHC23i1/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i1_MC_output_Tree_merged.root", info.path_to_file.Data());
        NameMC_file[1].Form("%s/output_HF_dimuons/mc_analysis_output/LHC23i2/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i2_MC_output_Tree_merged.root", info.path_to_file.Data());

        NameMC_file_x_CS_calc[0].Form("%s/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_charm_nocut_Version_5_AliAOD_withHF_Q/powheg_charm_nocut_Version_5_AliAOD_withHF_Q_MC_output_Hist_294154.root", info.path_to_file.Data());
        NameMC_file_x_CS_calc[1].Form("%s/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_beauty_nocut_Version_5_AliAOD_withHF_Q/powheg_beauty_nocut_Version_5_AliAOD_withHF_Q_MC_output_Hist_294154.root", info.path_to_file.Data());

        NameMC_Hist_x_CS_calc[0] = "HF_quarks/Powheg/h_NCharm_event_PowhegOnly";
        NameMC_Hist_x_CS_calc[1] = "HF_quarks/Powheg/h_NBeauty_event_PowhegOnly";

        NameMC_Hist_x_CS_calc_fwd[0] = "HF_quarks/Powheg/h_NCharm_event_fwd_PowhegOnly";
        NameMC_Hist_x_CS_calc_fwd[1] = "HF_quarks/Powheg/h_NBeauty_event_fwd_PowhegOnly";
    }
    else
    {
        Powheg_CS[0] = (56.42);
        Powheg_CS[1] = (56.42);

        fNorm_MC[0] = {216.};
        fNorm_MC[1] = {216.};

        NameMC_file[0].Form("%s/output_HF_dimuons/mc_analysis_output/LHC22b3/Version_5_AliAOD_skimmed_fwd_fullstat/LHC22b3_MC_output_Tree_merged.root", info.path_to_file.Data());
        NameMC_file[1].Form("%s/output_HF_dimuons/mc_analysis_output/LHC22b3/Version_5_AliAOD_skimmed_fwd_fullstat/LHC22b3_MC_output_Tree_merged.root", info.path_to_file.Data());

        // NameMC_file_x_CS_calc[0] = "/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/test/Hist_monash_sim_fixed_100Mev.root";
        // NameMC_file_x_CS_calc[1] = "/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/test/Hist_monash_sim_fixed_100Mev.root";

        NameMC_file_x_CS_calc[0].Form("%s/HF_dimuons/mc_analysis/read_output/cs_Generator_SoftQCD_Def_pythia_sim_12345_DefaultBR_1000000.root", info.path_to_file.Data());
        NameMC_file_x_CS_calc[1].Form("%s/HF_dimuons/mc_analysis/read_output/cs_Generator_SoftQCD_Def_pythia_sim_12345_DefaultBR_1000000.root", info.path_to_file.Data());

        // NameMC_file_x_CS_calc[0] = "/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_inel_LFoff_Def_pythia_sim_74_DefaultBR_output_Hist_100000.root";
        // NameMC_file_x_CS_calc[1] = "/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_inel_LFoff_Def_pythia_sim_74_DefaultBR_output_Hist_100000.root";

        // NameMC_Hist_x_CS_calc[0] = "HF_quarks/h_NCharm_event";
        // NameMC_Hist_x_CS_calc[1] = "HF_quarks/h_NBeauty_event";

        // NameMC_Hist_x_CS_calc_fwd[0] = "HF_quarks/h_NCharm_event_fwd";
        // NameMC_Hist_x_CS_calc_fwd[1] = "HF_quarks/h_NBeauty_event_fwd";

        // NameMC_Hist_x_CS_calc[0] = "CS/h_Ncharm_pairs_v6";
        // NameMC_Hist_x_CS_calc[1] = "CS/h_Nbeauty_pairs_v6";

        NameMC_Hist_x_CS_calc[0] = "DiCharm_rapidity";
        NameMC_Hist_x_CS_calc[1] = "DiBeauty_rapidity";

        NameMC_Hist_x_CS_calc_fwd[0] = "DiCharm_rapidity";
        NameMC_Hist_x_CS_calc_fwd[1] = "DiBeauty_rapidity";
    }
    std::ofstream out(TString::Format("results/fit_result/cross_section_%s_MC_%s_%s.txt", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data()), std::ios_base::app);
    TString output_fit;
    output_fit.Form("==============================================\n");
    out << output_fit.Data();
    output_fit.Form("Generator : %s\n", info.Generator.Data());
    out << output_fit.Data();
    output_fit.Form("Mass Range : [%d,%d] || Mass cut : ]%0.1f,%0.1f[\n", info.Low_Mass, info.High_Mass, info.LowM_cut, info.HighM_cut);
    out << output_fit.Data();
    output_fit.Form("Pt Range : [%d,%d] \n", info.Low_Pt, info.High_Pt);
    out << output_fit.Data();

    for (Int_t i = 0; i < n_obs; i++)
    {
        cout << "Result for " << NameDiMu_Sel[i].Data() << endl;
        cout << "===== MonteCarlo =====" << endl;
        out << "Result for " << NameDiMu_Sel[i].Data() << endl;
        out << "===== MonteCarlo =====" << endl;
        TFile *fIn_MC = new TFile(NameMC_file[i], "READ");
        gROOT->cd();
        TTree *m_tree;
        if (i < 2 && info.Generator.Contains("POWHEG"))
            m_tree = ((TTree *)fIn_MC->Get(Form("DiMuon_Rec_PowhegOnly_%s", NameDiMu_Sel[i].Data())));
        else
            m_tree = ((TTree *)fIn_MC->Get(Form("DiMuon_Rec_%s", NameDiMu_Sel[i].Data())));

        TTree *m_tree_MC_cutted;
        if (info.Low_Mass == 4 && info.High_Mass == 30)
            m_tree_MC_cutted = (TTree *)m_tree->CopyTree(Form("((m>%d && m<9.0) ||(m>11.0 && m<%d)) && pt<%d", info.Low_Mass, info.High_Mass, info.High_Pt)); // Number of rec dimuons in the MC
        else
            m_tree_MC_cutted = (TTree *)m_tree->CopyTree(Form("(m>%d && m<%d) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt)); // Number of rec dimuons in the MC
        // m_tree_MC_cutted->SetName(Form("Tree_%s_M_%d_%d_Pt_%d_%d", NameDiMu_Sel[i].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
        TH1F *NEv_hist;
        Double_t NEv_MB_MC = 0;
        NEv_hist = (TH1F *)(fIn_MC)->Get("h_Nevents");
        NEv_MB_MC = fNorm_MC[i] * NEv_hist->GetBinContent(2); // Number of MB event in the MC
        if (i < 2)
        {
        }
        // else
        //     NEv_MB_MC = fNorm_MC[i] * 637848;

        cout << "N events MB: " << NEv_MB_MC << endl;
        out << "N events MB: " << NEv_MB_MC << endl;
        fIn_MC->Close();
        Double_t NDimu_MC_MB = (m_tree_MC_cutted->GetEntries() / NEv_MB_MC); // Number of rec Dimuons x event in MC
        cout << "N Dimu rec: " << m_tree_MC_cutted->GetEntries() << endl;
        cout << "N Dimu rec x event: " << NDimu_MC_MB << endl;
        out << "N Dimu rec: " << m_tree_MC_cutted->GetEntries() << endl;
        out << "N Dimu rec x event: " << NDimu_MC_MB << endl;

        cout << "===== Data =====" << endl;
        out << "===== Data =====" << endl;

        Double_t fit_result = (r->GetBinContent(i + 1));
        Double_t fit_error = (r->GetBinError(i + 1));
        Double_t NDimu_MC_Data = (fit_result / NEv_MB_Data);
        cout << "N events MB: " << NEv_MB_Data << endl;
        cout << "N Dimu from fit: " << fit_result << " +/ " << fit_error << endl;
        cout << "N Dimu from fit x event: " << NDimu_MC_Data << endl;

        cout << "===== Cross section calculation =====" << endl;

        out << "N events MB: " << NEv_MB_Data << endl;
        out << "N Dimu from fit: " << fit_result << endl;
        out << "N Dimu from fit x event: " << NDimu_MC_Data << endl;

        out << "===== Cross section calculation =====" << endl;

        TFile *fIn_MC_x_CS = new TFile(NameMC_file_x_CS_calc[i], "READ");
        Double_t N_Pair = 0;
        Double_t N_Pair_fwd = 0;
        for (Int_t i_bin = 0; i_bin < ((TH1D *)fIn_MC_x_CS->Get(NameMC_Hist_x_CS_calc[i]))->GetNbinsX(); i_bin++)
        {
            N_Pair = (N_Pair + (0.5) * ((TH1D *)fIn_MC_x_CS->Get(NameMC_Hist_x_CS_calc[i]))->GetBinContent(i_bin) * (i_bin - 1));
            // printf("Bin %d || Content %0.0f\n", i_bin, h_NHF_event_PowhegOnly[i]->GetBinContent(i_bin));
            if (info.Generator.Contains("POWHEG"))
                N_Pair_fwd = (N_Pair_fwd + (0.5) * ((TH1D *)fIn_MC_x_CS->Get(NameMC_Hist_x_CS_calc_fwd[i]))->GetBinContent(i_bin) * (i_bin - 1));
            else
                N_Pair_fwd = 9999999;

            // if (h_NHF_event_Pythia[i] == nullptr)
            //     continue;
        }
        if (info.Generator.Contains("PYTHIA"))
        {
            ((TH1D *)fIn_MC_x_CS->Get(NameMC_Hist_x_CS_calc_fwd[i]))->GetXaxis()->SetRangeUser(-4, -2.5);
            N_Pair_fwd = ((TH1D *)fIn_MC_x_CS->Get(NameMC_Hist_x_CS_calc_fwd[i]))->Integral();
        }
        fIn_MC_x_CS->Close();

        Double_t CS = 999;
        if (info.Generator.Contains("POWHEG"))
            CS = ((Powheg_CS[i] * N_Pair_fwd / N_Pair) / (2 * 1.5));
        else
            CS = ((Powheg_CS[i] * N_Pair_fwd / 1e+6) / (1.5)); // PYTHIA sim equals to 1e+6 events

        Double_t Ratio = NDimu_MC_Data / NDimu_MC_MB;
        printf("N_Pair_fwd %0.1f\n", N_Pair_fwd);
        printf("Generator CS %0.3e\n", Powheg_CS[i]);
        printf("FWD generator CS %0.3e\n", CS);
        printf("Ratio %0.3e\n", Ratio);
        Double_t stat_CS_error = (Ratio * CS) * (fit_error / fit_result);
        printf("Measured Powheg CS %0.3e +/n %0.3e\n", (Ratio * CS), stat_CS_error);

        out << "N_Pair_fwd : " << N_Pair_fwd << endl;
        out << "Generator CS : " << Powheg_CS[i] << endl;
        out << "Powheg CS : " << CS << endl;
        out << "Ratio : " << Ratio << endl;
        out << "Measured Powheg CS : " << (Ratio * CS) << " +/- " << stat_CS_error << endl;

        if (i < 2)
        {
            cout << "===== MPI corr =====" << endl;
            out << "===== MPI corr =====" << endl;
            TFile *fIn_Pythia = new TFile(Form("%s/HF_dimuons/mc_analysis/read_output/cs_Generator_SoftQCD_Def_pythia_sim_12345_DefaultBR_1000000.root", info.path_to_file.Data()), "READ");
            Double_t n_HF_total_PYTHIA = 0.0;
            for (Int_t i_bin = 1; i_bin <= ((TH1F *)fIn_Pythia->Get(Form("h_N%s_event", NameDiMu_Sel[i].Data())))->GetNbinsX(); i_bin++)
            {
                n_HF_total_PYTHIA = (n_HF_total_PYTHIA + (0.5) * ((TH1F *)fIn_Pythia->Get(Form("h_N%s_event", NameDiMu_Sel[i].Data())))->GetBinContent(i_bin) * (i_bin - 1));
            }

            Double_t n_HF_single_pair_PYTHIA = ((TH1F *)fIn_Pythia->Get(Form("h_N%s_event", NameDiMu_Sel[i].Data())))->GetBinContent(3);
            printf("MPI correction with PYTHIA %0.3f\n", n_HF_total_PYTHIA / n_HF_single_pair_PYTHIA);
            printf("Measured CS corrected for MPI %0.4e +/- %0.4e\n", 1000 * (Ratio * CS * n_HF_total_PYTHIA / n_HF_single_pair_PYTHIA), 1000 * stat_CS_error * (n_HF_total_PYTHIA / n_HF_single_pair_PYTHIA));
            out << "MPI correction with PYTHIA %0.3f\n"
                << (n_HF_total_PYTHIA / n_HF_single_pair_PYTHIA) << endl;
            out << "Measured CS corrected for MPI %0.3e\n"
                << (Ratio * CS * n_HF_total_PYTHIA / n_HF_single_pair_PYTHIA) << endl;
        }
        cout << "=================================================" << endl;
        out << "=================================================" << endl;
    }
    out.close();
}

void syst_plot()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    opt info;

    TString dir = Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut);
    Double_t fit_result[n_DiMuSelection];
    Double_t fit_Low2Up[n_DiMuSelection];
    Double_t fit_Up2Low[n_DiMuSelection];

    TFile *fIn_result = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/fit_result/Result_Fit_Data_%s_MC_%s_%s.root", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data()), "READ");
    TH1F *hist_result = (TH1F *)fIn_result->Get(Form("%s/hist_fit_result_%s_%s", dir.Data(), info.DY.Data(), info.LF_HF.Data()));
    hist_result->SetName("hist_result");
    hist_result->SetTitle("original");
    hist1D_graphic_opt(hist_result, kFALSE, 1, 20, kMagenta + 2, 1.);
    for (Int_t i = 0; i < hist_result->GetNbinsX(); i++)
        hist_result->GetXaxis()->SetBinLabel(i + 1, "");
    hist_result->GetXaxis()->SetTitle("");

    TFile *fIn_Low2Up = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/systematic/Syst_scaled_Low2Up_Fit_Data_%s_MC_%s_%s.root", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data()), "READ");
    TH1F *hist_Low2Up = (TH1F *)fIn_Low2Up->Get(Form("%s/hist_fit_result_%s_%s", dir.Data(), info.DY.Data(), info.LF_HF.Data()));
    hist_Low2Up->SetName("hist_Low2Up");
    hist_Low2Up->SetTitle("Low2Up variation");
    hist1D_graphic_opt(hist_Low2Up, kFALSE, 1, 22, kAzure + 2, 1.);

    TFile *fIn_Up2Low = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/systematic/Syst_scaled_Up2Low_Fit_Data_%s_MC_%s_%s.root", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data()), "READ");
    TH1F *hist_Up2Low = (TH1F *)fIn_Up2Low->Get(Form("%s/hist_fit_result_%s_%s", dir.Data(), info.DY.Data(), info.LF_HF.Data()));
    hist_Up2Low->SetName("hist_Up2Low");
    hist_Up2Low->SetTitle("Up2Low variation");
    hist1D_graphic_opt(hist_Up2Low, kFALSE, 1, 23, kGreen + 2, 1.);

    TCanvas *Charm = canvas_noratio(Form("Charm_systematic_%s_M_%d_%d_Pt_%d_%d_%s", info.Generator.Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LF_HF.Data()));
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.095);

    hist_result->GetXaxis()->SetRangeUser(0, 1);
    if (hist_Low2Up->GetBinContent(1) < hist_result->GetBinContent(1))
        hist_result->GetYaxis()->SetRangeUser(hist_Low2Up->GetBinContent(1) * 0.9, hist_Up2Low->GetBinContent(1) * 1.2);
    else
        hist_result->GetYaxis()->SetRangeUser(hist_Up2Low->GetBinContent(1) * 0.9, hist_Low2Up->GetBinContent(1) * 1.2);

    hist_result->DrawCopy("PE text0");
    hist_Low2Up->DrawCopy("PE text0same");
    hist_Up2Low->DrawCopy("PE text0same");
    TLegend *legend = (TLegend *)gPad->BuildLegend();
    TString Header = Form("Charm #mu#mu, %s %d < #it{m}_{#mu#mu} < %d, %s", info.Generator.Data(), info.Low_Mass, info.High_Mass, info.LF_HF.Data());
    Legend_settings(legend, 0.175, 0.3, 0.725, 0.925, Header);
    legend->DrawClone("same");
    Charm->SaveAs(Form("systematic/images/%s.pdf", Charm->GetName()));

    TCanvas *Beauty = canvas_noratio(Form("Beauty_systematic_%s_M_%d_%d_Pt_%d_%d_%s", info.Generator.Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LF_HF.Data()));
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.095);

    hist_result->GetXaxis()->SetRangeUser(1, 2);
    if (hist_Low2Up->GetBinContent(2) < hist_result->GetBinContent(2))
        hist_result->GetYaxis()->SetRangeUser(hist_Low2Up->GetBinContent(2) * 0.9, hist_Up2Low->GetBinContent(2) * 1.3);
    else
        hist_result->GetYaxis()->SetRangeUser(hist_Up2Low->GetBinContent(2) * 0.9, hist_Low2Up->GetBinContent(2) * 1.3);

    hist_result->DrawCopy("PE text0");
    hist_Low2Up->DrawCopy("PE text0same");
    hist_Up2Low->DrawCopy("PE text0same");
    delete legend;
    legend = (TLegend *)gPad->BuildLegend();
    Header.Form("Beauty #mu#mu, %s %d < #it{m}_{#mu#mu} < %d, %s", info.Generator.Data(), info.Low_Mass, info.High_Mass, info.LF_HF.Data());
    Legend_settings(legend, 0.175, 0.3, 0.725, 0.925, Header);
    legend->Draw("same");
    Beauty->SaveAs(Form("systematic/images/%s.pdf", Beauty->GetName()));

    // for (Int_t i_Dimu = 0; i_Dimu < n_DiMuSelection; i_Dimu++)
    // {
    //     fit_result[i_Dimu] = hist_result->GetBinContent(i_Dimu + 1);
    //     fit_Low2Up[i_Dimu] = hist_Low2Up->GetBinContent(i_Dimu + 1);
    //     fit_Up2Low[i_Dimu] = hist_Up2Low->GetBinContent(i_Dimu + 1);

    //     cout << info.Name_DimuSel[i_Dimu].Data() << ": " << fit_result[i_Dimu] << endl;
    //     cout << info.Name_DimuSel[i_Dimu].Data() << ": " << fit_Low2Up[i_Dimu] << endl;
    //     cout << info.Name_DimuSel[i_Dimu].Data() << ": " << fit_Up2Low[i_Dimu] << endl;
    // }
}

void MPI_question(TString HF = "Beauty")
{
    opt info;
    Int_t N_bins_Pt_hist = 999;
    Int_t N_bins_Mass_hist = 999;

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        N_bins_Pt_hist = 60;
        N_bins_Mass_hist = 52;
    }
    else
    {
        N_bins_Pt_hist = 60;
        N_bins_Mass_hist = 30;
    }

    TFile *fIn_LHC23i1_Tree = new TFile(Form("%s/output_HF_dimuons/mc_analysis_output/LHC23i1/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i1_MC_output_Tree_merged.root", info.path_to_file.Data()), "READ");
    TTree *fTree_POWHEGOnly = (TTree *)fIn_LHC23i1_Tree->Get(Form("DiMuon_Rec_PowhegOnly_%s", HF.Data()));
    TH1F *Hist_Pt_POWHEGOnly = new TH1F(Form("Hist_Pt_%s_POWHEGOnly", HF.Data()), Form("Hist_Pt_%s_POWHEGOnly; #it{p}_{T} (GeV/#it{c});", HF.Data()), N_bins_Pt_hist, info.Low_Pt, info.High_Pt);

    TTree *fTree_PYTHIAOnly = (TTree *)fIn_LHC23i1_Tree->Get(Form("DiMuon_Rec_PythiaOnly_%s", HF.Data()));
    TH1F *Hist_Pt_PYTHIAOnly = new TH1F(Form("Hist_Pt_%s_PYTHIAOnly", HF.Data()), Form("Hist_Pt_%s_PYTHIAOnly; #it{p}_{T} (GeV/#it{c});", HF.Data()), N_bins_Pt_hist, info.Low_Pt, info.High_Pt);

    TTree *fTree_All = (TTree *)fIn_LHC23i1_Tree->Get(Form("DiMuon_Rec_%s", HF.Data()));
    TH1F *Hist_Pt_All = new TH1F(Form("Hist_Pt_%s_All", HF.Data()), Form("Hist_Pt_%s_All; #it{p}_{T} (GeV/#it{c});", HF.Data()), N_bins_Pt_hist, info.Low_Pt, info.High_Pt);

    fTree_POWHEGOnly->Draw(Form("pt>>Hist_Pt_%s_POWHEGOnly", HF.Data()), Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    fTree_PYTHIAOnly->Draw(Form("pt>>Hist_Pt_%s_PYTHIAOnly", HF.Data()), Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    fTree_All->Draw(Form("pt>>Hist_Pt_%s_All", HF.Data()), Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");

    hist1D_graphic_opt(Hist_Pt_POWHEGOnly, kTRUE, 1, 20, kMagenta + 2, 1. / Hist_Pt_POWHEGOnly->GetEntries());
    hist1D_graphic_opt(Hist_Pt_PYTHIAOnly, kTRUE, 1, 20, kRed + 2, 1. / Hist_Pt_PYTHIAOnly->GetEntries());
    hist1D_graphic_opt(Hist_Pt_All, kTRUE, 1, 20, kAzure + 2, 1. / Hist_Pt_All->GetEntries());

    // hist1D_graphic_opt(Hist_Pt_POWHEGOnly, kTRUE, 1, 20, kMagenta + 2, 1.);
    // hist1D_graphic_opt(Hist_Pt_PYTHIAOnly, kTRUE, 1, 20, kRed + 2, 1.);
    // hist1D_graphic_opt(Hist_Pt_All, kTRUE, 1, 20, kAzure + 2, 1.);

    Hist_Pt_POWHEGOnly->Draw("PESAME");
    Hist_Pt_PYTHIAOnly->Draw("PESAME");
    Hist_Pt_All->Draw("PESAME");

    gPad->BuildLegend();
}

void HF_mixing_POWHEG()
{

    TVirtualFitter::SetDefaultFitter("Minuit2");
    opt info;
    Int_t N_bins_Pt_hist = 999;
    Int_t N_bins_Mass_hist = 999;

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        N_bins_Pt_hist = 60;
        N_bins_Mass_hist = 52;
    }
    else
    {
        N_bins_Pt_hist = 60;
        N_bins_Mass_hist = 30;
    }

    TFile *fIn_LHC23i1_Tree = new TFile(Form("%s/output_HF_dimuons/mc_analysis_output/LHC23i1/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i1_MC_output_Tree_merged.root", info.path_to_file.Data()), "READ");

    TTree *fTree_LHC23i1 = (TTree *)fIn_LHC23i1_Tree->Get("DiMuon_Rec_PowhegOnly_HF_Mixed");
    TH1F *Hist_Pt_HF_LHC23i1 = new TH1F("Hist_Pt_HF_LHC23i1", "Hist_Pt_HF_LHC23i1; #it{p}_{T} (GeV/#it{c});", N_bins_Pt_hist, info.Low_Pt, info.High_Pt);
    TH1F *Hist_Mass_HF_LHC23i1 = new TH1F("Hist_Mass_HF_LHC23i1", "Hist_Mass_HF_LHC23i1 ;#it{m}_{#mu#mu} (GeV/#it{c}^{2});", N_bins_Mass_hist, info.Low_Mass, info.High_Mass);

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        fTree_LHC23i1->Draw("pt>>Hist_Pt_HF_LHC23i1", Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
        fTree_LHC23i1->Draw("m>>Hist_Mass_HF_LHC23i1", Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    }
    else
    {
        fTree_LHC23i1->Draw("pt>>Hist_Pt_HF_LHC23i1", Form("((m>%d && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
        fTree_LHC23i1->Draw("m>>Hist_Mass_HF_LHC23i1", Form("((m>%d && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    }

    Hist_Pt_HF_LHC23i1->Draw();

    TF1 *Pt_HF_LHC23i1_BINNED = new TF1("Pt_HF_LHC23i1_BINNED", FuncPt, info.Low_Pt, info.High_Pt, 4);
    Pt_HF_LHC23i1_BINNED->SetParameter(0, 1.9703e+01);
    Pt_HF_LHC23i1_BINNED->SetParameter(1, 1.4663e+00);
    Pt_HF_LHC23i1_BINNED->SetParameter(2, 1.6398e+01);
    Pt_HF_LHC23i1_BINNED->SetParameter(3, 2e+2);

    Hist_Pt_HF_LHC23i1->Fit(Pt_HF_LHC23i1_BINNED, "LR0I");
    TCanvas *c_PT_test_LHC23i1 = new TCanvas("c_PT_test_LHC23i1", "c_PT_test_LHC23i1", 1200, 1200);
    TH1F *Error_PT_LHC23i1 = (TH1F *)Hist_Pt_HF_LHC23i1->Clone("Error_PT_LHC23i1");
    Error_PT_LHC23i1->SetTitle("Error_PT_LHC23i1");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(Error_PT_LHC23i1, 0.95);
    Error_PT_LHC23i1->SetFillColor(kRed - 10);
    Error_PT_LHC23i1->Draw("e3");
    Hist_Pt_HF_LHC23i1->Draw("PEsame");
    Pt_HF_LHC23i1_BINNED->Draw("same");

    gPad->BuildLegend();

    TF1 *Mass_HF_LHC23i1_BINNED = new TF1("Mass_HF_LHC23i1_BINNED", FuncMass, info.Low_Mass, info.High_Mass, 4);
    Mass_HF_LHC23i1_BINNED->SetParameter(0, 1.9703e+01);
    Mass_HF_LHC23i1_BINNED->SetParameter(1, 1.4663e+00);
    Mass_HF_LHC23i1_BINNED->SetParameter(2, 1.6398e+01);
    Mass_HF_LHC23i1_BINNED->SetParLimits(3, 0, 5e+6);

    Hist_Mass_HF_LHC23i1->Fit(Mass_HF_LHC23i1_BINNED, "LR0I");
    TCanvas *c_Mass_test_LHC23i1 = new TCanvas("c_Mass_test_LHC23i1", "c_Mass_test_LHC23i1", 1200, 1200);
    TH1F *Error_Mass_LHC23i1 = (TH1F *)Hist_Mass_HF_LHC23i1->Clone("Error_Mass_LHC23i1");
    Error_Mass_LHC23i1->SetTitle("Error_Mass_LHC23i1");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(Error_Mass_LHC23i1, 0.95);
    Error_Mass_LHC23i1->SetFillColor(kRed - 10);
    Error_Mass_LHC23i1->Draw("e3");
    Hist_Mass_HF_LHC23i1->Draw("PEsame");
    Mass_HF_LHC23i1_BINNED->Draw("same");

    gPad->BuildLegend();

    // LHC23i2

    TFile *fIn_LHC23i2_Tree = new TFile(Form("%s/output_HF_dimuons/mc_analysis_output/LHC23i2/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i2_MC_output_Tree_merged.root", info.path_to_file.Data()), "READ");

    TTree *fTree_LHC23i2 = (TTree *)fIn_LHC23i2_Tree->Get("DiMuon_Rec_PowhegOnly_HF_Mixed");
    TH1F *Hist_Pt_HF_LHC23i2 = new TH1F("Hist_Pt_HF_LHC23i2", "Hist_Pt_HF_LHC23i2; #it{p}_{T} (GeV/#it{c});", N_bins_Pt_hist, info.Low_Pt, info.High_Pt);
    TH1F *Hist_Mass_HF_LHC23i2 = new TH1F("Hist_Mass_HF_LHC23i2", "Hist_Mass_HF_LHC23i2 ;#it{m}_{#mu#mu} (GeV/#it{c}^{2});", N_bins_Mass_hist, info.Low_Mass, info.High_Mass);

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        fTree_LHC23i2->Draw("pt>>Hist_Pt_HF_LHC23i2", Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
        fTree_LHC23i2->Draw("m>>Hist_Mass_HF_LHC23i2", Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    }
    else
    {
        fTree_LHC23i2->Draw("pt>>Hist_Pt_HF_LHC23i2", Form("((m>%d && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
        fTree_LHC23i2->Draw("m>>Hist_Mass_HF_LHC23i2", Form("((m>%d && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    }

    Hist_Pt_HF_LHC23i2->Draw();

    TF1 *Pt_HF_LHC23i2_BINNED = new TF1("Pt_HF_LHC23i2_BINNED", FuncPt, info.Low_Pt, info.High_Pt, 4);
    Pt_HF_LHC23i2_BINNED->SetParameter(0, 1.9703e+01);
    Pt_HF_LHC23i2_BINNED->SetParameter(1, 1.4663e+00);
    Pt_HF_LHC23i2_BINNED->SetParameter(2, 1.6398e+01);
    Pt_HF_LHC23i2_BINNED->SetParameter(3, 2e+2);

    Hist_Pt_HF_LHC23i2->Fit(Pt_HF_LHC23i2_BINNED, "LR0I");
    TCanvas *c_PT_test_LHC23i2 = new TCanvas("c_PT_test_LHC23i2", "c_PT_test_LHC23i2", 1200, 1200);
    TH1F *Error_PT_LHC23i2 = (TH1F *)Hist_Pt_HF_LHC23i2->Clone("Error_PT_LHC23i2");
    Error_PT_LHC23i2->SetTitle("Error_PT_LHC23i2");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(Error_PT_LHC23i2, 0.95);
    Error_PT_LHC23i2->SetFillColor(kRed - 10);
    Error_PT_LHC23i2->Draw("e3");
    Hist_Pt_HF_LHC23i2->Draw("PEsame");
    Pt_HF_LHC23i2_BINNED->Draw("same");

    gPad->BuildLegend();

    TF1 *Mass_HF_LHC23i2_BINNED = new TF1("Mass_HF_LHC23i2_BINNED", FuncMass, info.Low_Mass, info.High_Mass, 4);
    Mass_HF_LHC23i2_BINNED->SetParameter(0, 1.9703e+01);
    Mass_HF_LHC23i2_BINNED->SetParameter(1, 1.4663e+00);
    Mass_HF_LHC23i2_BINNED->SetParameter(2, 1.6398e+01);
    Mass_HF_LHC23i2_BINNED->SetParLimits(3, 0, 5e+7);

    Hist_Mass_HF_LHC23i2->Fit(Mass_HF_LHC23i2_BINNED, "LR0I");
    TCanvas *c_Mass_test_LHC23i2 = new TCanvas("c_Mass_test_LHC23i2", "c_Mass_test_LHC23i2", 1200, 1200);
    TH1F *Error_Mass_LHC23i2 = (TH1F *)Hist_Mass_HF_LHC23i2->Clone("Error_Mass_LHC23i2");
    Error_Mass_LHC23i2->SetTitle("Error_Mass_LHC23i2");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(Error_Mass_LHC23i2, 0.95);
    Error_Mass_LHC23i2->SetFillColor(kRed - 10);
    Error_Mass_LHC23i2->Draw("e3");
    Hist_Mass_HF_LHC23i2->Draw("PEsame");
    Mass_HF_LHC23i2_BINNED->Draw("same");

    gPad->BuildLegend();

    Int_t N_bins_Pt_shape = 999;
    Int_t N_bins_Mass_shape = 999;

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        N_bins_Pt_shape = 150;
        N_bins_Mass_shape = 104;
    }
    else
    {
        N_bins_Pt_shape = 100;
        N_bins_Mass_shape = 50;
    }

    TH1F *Shape_Pt_HF_LHC23i1 = new TH1F("Shape_Pt_HF_LHC23i1", "Shape_Pt_HF_LHC23i1", N_bins_Pt_shape, info.Low_Pt, info.High_Pt);
    TH1F *Shape_Pt_HF_LHC23i2 = new TH1F("Shape_Pt_HF_LHC23i2", "Shape_Pt_HF_LHC23i2", N_bins_Pt_shape, info.Low_Pt, info.High_Pt);
    TH1F *Shape_Pt_HF_mid = new TH1F("Shape_Pt_HF_mid", "Shape_Pt_HF_mid", N_bins_Pt_shape, info.Low_Pt, info.High_Pt);

    TH1F *Shape_Mass_HF_LHC23i1 = new TH1F("Shape_Mass_HF_LHC23i1", "Shape_Mass_HF_LHC23i1", N_bins_Mass_shape, info.Low_Mass, info.High_Mass);
    TH1F *Shape_Mass_HF_LHC23i2 = new TH1F("Shape_Mass_HF_LHC23i2", "Shape_Mass_HF_LHC23i2", N_bins_Mass_shape, info.Low_Mass, info.High_Mass);
    TH1F *Shape_Mass_HF_mid = new TH1F("Shape_Mass_HF_mid", "Shape_Mass_HF_mid", N_bins_Mass_shape, info.Low_Mass, info.High_Mass);

    for (Int_t i_pt_bin = 0; i_pt_bin < Shape_Pt_HF_LHC23i1->GetNbinsX(); i_pt_bin++)
    {
        Shape_Pt_HF_LHC23i1->SetBinContent(i_pt_bin + 1, Pt_HF_LHC23i1_BINNED->Eval(Shape_Pt_HF_LHC23i1->GetBinCenter(i_pt_bin + 1)));
        Shape_Pt_HF_LHC23i1->SetBinError(i_pt_bin + 1, Error_PT_LHC23i1->GetBinError(i_pt_bin + 1));
        Shape_Pt_HF_LHC23i2->SetBinContent(i_pt_bin + 1, Pt_HF_LHC23i2_BINNED->Eval(Shape_Pt_HF_LHC23i2->GetBinCenter(i_pt_bin + 1)));
        Shape_Pt_HF_LHC23i2->SetBinError(i_pt_bin + 1, Error_PT_LHC23i2->GetBinError(i_pt_bin + 1));
    }
    Shape_Pt_HF_LHC23i1->Scale(1. / Shape_Pt_HF_LHC23i1->Integral(), "width");
    Shape_Pt_HF_LHC23i2->Scale(1. / Shape_Pt_HF_LHC23i2->Integral(), "width");

    for (Int_t i_pt_bin = 0; i_pt_bin < Shape_Pt_HF_LHC23i1->GetNbinsX(); i_pt_bin++)
    {
        Shape_Pt_HF_mid->SetBinContent(i_pt_bin + 1, (Shape_Pt_HF_LHC23i1->GetBinContent(i_pt_bin + 1) + Shape_Pt_HF_LHC23i2->GetBinContent(i_pt_bin + 1)) / 2);
        Shape_Pt_HF_mid->SetBinError(i_pt_bin + 1, TMath::Sqrt(0.5 * TMath::Power(Shape_Pt_HF_LHC23i1->GetBinError(i_pt_bin + 1), 2) + 0.5 * TMath::Power(Shape_Pt_HF_LHC23i2->GetBinError(i_pt_bin + 1), 2)));
        // cout << "error LHC23i1: " << Shape_Pt_HF_LHC23i1->GetBinError(i_pt_bin + 1) << endl;
        // cout << "error LHC23i2: " << Shape_Pt_HF_LHC23i2->GetBinError(i_pt_bin + 1) << endl;

        // cout << "error mid: " << Shape_Pt_HF_mid->GetBinError(i_pt_bin + 1) << endl;
    }

    for (Int_t i_mass_bin = 0; i_mass_bin < Shape_Mass_HF_LHC23i1->GetNbinsX(); i_mass_bin++)
    {
        Shape_Mass_HF_LHC23i1->SetBinContent(i_mass_bin + 1, Mass_HF_LHC23i1_BINNED->Eval(Shape_Mass_HF_LHC23i1->GetBinCenter(i_mass_bin + 1)));
        Shape_Mass_HF_LHC23i1->SetBinError(i_mass_bin + 1, Error_Mass_LHC23i1->GetBinError(i_mass_bin + 1));
        Shape_Mass_HF_LHC23i2->SetBinContent(i_mass_bin + 1, Mass_HF_LHC23i2_BINNED->Eval(Shape_Mass_HF_LHC23i2->GetBinCenter(i_mass_bin + 1)));
        Shape_Mass_HF_LHC23i2->SetBinError(i_mass_bin + 1, Error_Mass_LHC23i2->GetBinError(i_mass_bin + 1));
    }

    Shape_Mass_HF_LHC23i1->Scale(1. / Shape_Mass_HF_LHC23i1->Integral(), "width");
    Shape_Mass_HF_LHC23i2->Scale(1. / Shape_Mass_HF_LHC23i2->Integral(), "width");

    for (Int_t i_mass_bin = 0; i_mass_bin < Shape_Mass_HF_LHC23i1->GetNbinsX(); i_mass_bin++)
    {
        Shape_Mass_HF_mid->SetBinContent(i_mass_bin + 1, (Shape_Mass_HF_LHC23i1->GetBinContent(i_mass_bin + 1) + Shape_Mass_HF_LHC23i2->GetBinContent(i_mass_bin + 1)) / 2);
        Shape_Mass_HF_mid->SetBinError(i_mass_bin + 1, TMath::Sqrt(0.5 * TMath::Power(Shape_Mass_HF_LHC23i1->GetBinError(i_mass_bin + 1), 2) + 0.5 * TMath::Power(Shape_Mass_HF_LHC23i2->GetBinError(i_mass_bin + 1), 2)));
        // cout << "error LHC23i1: " << Shape_Mass_LF_HF_LHC23i1->GetBinError(i_mass_bin + 1) << endl;
        // cout << "error LHC23i2: " << Shape_Mass_LF_HF_LHC23i2->GetBinError(i_mass_bin + 1) << endl;
        // cout << "error mid: " << Shape_Mass_LF_HF_mid->GetBinError(i_mass_bin + 1) << endl;
    }
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    TCanvas *c_pt = canvas_noratio("c_pt");
    gPad->SetLogy();
    // Shape_Pt_LF_HF_mid->Scale(1. / Shape_Pt_LF_HF_mid->Integral(), "width");
    Shape_Pt_HF_LHC23i1->SetTitle("from LHC23i1");
    Shape_Pt_HF_LHC23i2->SetTitle("from LHC23i2");
    Shape_Pt_HF_mid->SetTitle("Mixing");
    Shape_Pt_HF_LHC23i1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    Shape_Pt_HF_LHC23i1->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    hist1D_graphic_opt(Shape_Pt_HF_LHC23i1, kFALSE, 1, 24, kMagenta + 2, 1.);
    hist1D_graphic_opt(Shape_Pt_HF_LHC23i2, kFALSE, 1, 24, kGreen + 2, 1.);
    hist1D_graphic_opt(Shape_Pt_HF_mid, kFALSE, 1, 20, kRed - 2, 1.);
    gPad->BuildLegend();

    Shape_Pt_HF_LHC23i1->Draw("PE");
    Shape_Pt_HF_LHC23i2->Draw("PE SAME");
    Shape_Pt_HF_mid->DrawCopy("PE SAME");
    gPad->BuildLegend();

    // c_pt->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Pt_shapes_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    TCanvas *c_pt_for_pres = canvas_noratio("c_pt_for_pres");
    gPad->SetLogy();
    Shape_Pt_HF_LHC23i1->DrawCopy("PE");
    Shape_Pt_HF_LHC23i2->DrawCopy("PE SAME");
    gPad->BuildLegend();
    // c_pt_for_pres->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Pt_shapes_forPres_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    hist1D_graphic_opt(Shape_Pt_HF_mid, kFALSE, 1, 24, kRed - 2, 1.);
    TF1 *Pt_HF_MID = new TF1("Pt_HF_MID", FuncPt, info.Low_Pt, info.High_Pt, 4);
    Pt_HF_MID->SetParameter(0, Pt_HF_LHC23i2_BINNED->GetParameter(0));
    Pt_HF_MID->SetParameter(1, Pt_HF_LHC23i2_BINNED->GetParameter(1));
    Pt_HF_MID->SetParameter(2, Pt_HF_LHC23i2_BINNED->GetParameter(2));
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        Pt_HF_MID->SetParLimits(3, -10000, 10000);
    else if (info.Low_Mass == 4 && info.High_Mass == 9)
        Pt_HF_MID->SetParameter(3, 1500);

    Shape_Pt_HF_mid->Fit(Pt_HF_MID, "R0");

    // TCanvas *c_PT_mid = canvas_ratio("c_PT_mid");
    // c_PT_mid->cd(1);
    // Shape_Pt_HF_mid->Draw();
    // Pt_HF_MID->Draw("same");

    // TH1F *Shape_Pt_HF_mid_ratio = (TH1F *)Shape_Pt_HF_mid->Clone("ddjs");
    // Shape_Pt_HF_mid_ratio->Divide(Pt_HF_MID);
    // c_PT_mid->cd(2);
    // Shape_Pt_HF_mid_ratio->Draw();
    Shape_Pt_HF_mid->GetYaxis()->SetRangeUser(Shape_Pt_HF_mid->GetMinimum() * 0.2, Shape_Pt_HF_mid->GetMaximum() * 2e+2);
    TString Title = TString::Format("%d < #it{m}_{#mu#mu} < %d (GeV/#it{c}^{2}), PYTHIA only", info.Low_Mass, info.High_Mass);
    TCanvas *c_PT_mid = printMC_ratio_BINNED("c_PT_mid", Title, Shape_Pt_HF_mid, Pt_HF_MID, kBlack);
    // c_PT_mid->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Pt_extr_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));

    TCanvas *c_Pt_diff = canvas_noratio("c_Pt_diff");
    c_Pt_diff->cd();
    TH1F *Pt_Diff1 = (TH1F *)Shape_Pt_HF_LHC23i2->Clone("Pt_Diff1");
    Pt_Diff1->Add(Shape_Pt_HF_mid, -1);
    TH1F *Pt_Diff2 = (TH1F *)Shape_Pt_HF_LHC23i1->Clone("Pt_Diff2");
    Pt_Diff2->Add(Shape_Pt_HF_mid, -1);
    Pt_Diff1->GetYaxis()->SetRangeUser(-0.04, 0.04);
    Pt_Diff1->Draw("P hist ");
    Pt_Diff2->Draw("P hist same");
    // c_Pt_diff->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Pt_test_diff_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
}

void Mixing_background()
{
    // Getting paramater from unbinned fit to test the difference with binned fit and to anchor the parameters
    opt info;
    Int_t N_bins_Pt_hist = 999;
    Int_t N_bins_Mass_hist = 999;

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        N_bins_Pt_hist = 60;
        N_bins_Mass_hist = 52;
    }
    else
    {
        N_bins_Pt_hist = 60;
        N_bins_Mass_hist = 30;
    }
    TFile *fIn_LHC23i1_TEST = new TFile(Form("%s/HF_dimuons/fit_data/results/pdf_extraction/POWHEG_PDF_withLF_HF_LHC23i1_Mcut_%0.1f_%0.1f.root", info.path_to_file.Data(), info.LowM_cut, info.HighM_cut), "READ");
    TF1 *Pt_LF_HF_LHC23i1_TEST = (TF1 *)fIn_LHC23i1_TEST->Get(Form("pdfDimuPtFromLF_HF_Mixed_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    Pt_LF_HF_LHC23i1_TEST->SetName("Pt_LF_HF_LHC23i1_TEST");
    Pt_LF_HF_LHC23i1_TEST->SetTitle("Pt_LF_HF_LHC23i1_TEST");
    Pt_LF_HF_LHC23i1_TEST->SetLineColor(kGreen);

    TF1 *Mass_LF_HF_LHC23i1_TEST = (TF1 *)fIn_LHC23i1_TEST->Get(Form("pdfDimuMassFromLF_HF_Mixed_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    Mass_LF_HF_LHC23i1_TEST->SetName("Mass_LF_HF_LHC23i1_TEST");
    Mass_LF_HF_LHC23i1_TEST->SetTitle("Mass_LF_HF_LHC23i1_TEST");
    Mass_LF_HF_LHC23i1_TEST->SetLineColor(kGreen);

    TFile *fIn_LHC23i2_TEST = new TFile(Form("%s/HF_dimuons/fit_data/results/pdf_extraction/POWHEG_PDF_withLF_HF_LHC23i2_Mcut_%0.1f_%0.1f.root", info.path_to_file.Data(), info.LowM_cut, info.HighM_cut), "READ");

    TF1 *Pt_LF_HF_LHC23i2_TEST = (TF1 *)fIn_LHC23i2_TEST->Get(Form("pdfDimuPtFromLF_HF_Mixed_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    Pt_LF_HF_LHC23i2_TEST->SetName("Pt_LF_HF_LHC23i2_TEST");
    Pt_LF_HF_LHC23i2_TEST->SetTitle("Pt_LF_HF_LHC23i2_TEST");
    Pt_LF_HF_LHC23i2_TEST->SetLineColor(kOrange);

    TF1 *Mass_LF_HF_LHC23i2_TEST = (TF1 *)fIn_LHC23i2_TEST->Get(Form("pdfDimuMassFromLF_HF_Mixed_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    Mass_LF_HF_LHC23i2_TEST->SetName("Mass_LF_HF_LHC23i2_TEST");
    Mass_LF_HF_LHC23i2_TEST->SetTitle("Mass_LF_HF_LHC23i2_TEST");
    Mass_LF_HF_LHC23i2_TEST->SetLineColor(kOrange);

    // Get tree to fit

    TFile *fIn_LHC23i1_Tree = new TFile(Form("%s/output_HF_dimuons/mc_analysis_output/LHC23i1/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i1_MC_output_Tree_merged.root", info.path_to_file.Data()), "READ");

    TTree *fTree_LHC23i1 = (TTree *)fIn_LHC23i1_Tree->Get("DiMuon_Rec_PythiaOnly_LF_HF_Mixed");
    TH1F *Hist_Pt_LF_HF_LHC23i1 = new TH1F("Hist_Pt_LF_HF_LHC23i1", "Hist_Pt_LF_HF_LHC23i1; #it{p}_{T} (GeV/#it{c});", N_bins_Pt_hist, info.Low_Pt, info.High_Pt);
    TH1F *Hist_Mass_LF_HF_LHC23i1 = new TH1F("Hist_Mass_LF_HF_LHC23i1", "Hist_Mass_LF_HF_LHC23i1 ;#it{m}_{#mu#mu} (GeV/#it{c}^{2});", N_bins_Mass_hist, info.Low_Mass, info.High_Mass);

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        fTree_LHC23i1->Draw("pt>>Hist_Pt_LF_HF_LHC23i1", Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
        fTree_LHC23i1->Draw("m>>Hist_Mass_LF_HF_LHC23i1", Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    }
    else
    {
        fTree_LHC23i1->Draw("pt>>Hist_Pt_LF_HF_LHC23i1", Form("((m>%d && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
        fTree_LHC23i1->Draw("m>>Hist_Mass_LF_HF_LHC23i1", Form("((m>%d && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    }

    TF1 *Pt_LF_HF_LHC23i1_BINNED = new TF1("Pt_LF_HF_LHC23i1_BINNED", FuncPt, info.Low_Pt, info.High_Pt, 4);
    Pt_LF_HF_LHC23i1_BINNED->SetParameter(0, Pt_LF_HF_LHC23i1_TEST->GetParameter(0));
    Pt_LF_HF_LHC23i1_BINNED->SetParameter(1, Pt_LF_HF_LHC23i1_TEST->GetParameter(1));
    Pt_LF_HF_LHC23i1_BINNED->SetParameter(2, Pt_LF_HF_LHC23i1_TEST->GetParameter(2));
    Pt_LF_HF_LHC23i1_BINNED->SetParLimits(3, 0., 1000);

    // Hist_Pt_LF_HF_LHC23i1->Scale(1. / Hist_Pt_LF_HF_LHC23i1->GetEntries(), "width");
    Hist_Pt_LF_HF_LHC23i1->Fit(Pt_LF_HF_LHC23i1_BINNED, "LR0");
    TCanvas *c_PT_test_LHC23i1 = new TCanvas("c_PT_test_LHC23i1", "c_PT_test_LHC23i1", 1200, 1200);
    Hist_Pt_LF_HF_LHC23i1->Draw("PE");
    Pt_LF_HF_LHC23i1_BINNED->Draw("same");
    TH1F *Error_PT_LHC23i1 = (TH1F *)Hist_Pt_LF_HF_LHC23i1->Clone("Error_PT_LHC23i1");
    Error_PT_LHC23i1->SetTitle("Error_PT_LHC23i1");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(Error_PT_LHC23i1, 0.95);
    Error_PT_LHC23i1->SetFillColorAlpha(kBlue, 0.3);
    Error_PT_LHC23i1->Draw("e3same");

    gPad->BuildLegend();
    TH1F *Pt_grid = (TH1F *)Hist_Pt_LF_HF_LHC23i1->Clone("Pt_grid");
    Pt_grid->Reset();
    Pt_grid->GetYaxis()->SetRangeUser(1.25e-6, 1.25e+1);
    TCanvas *c_Pt_roofit_ratio_LHC23i1 = canvas_ratio("c_Pt_roofit_ratio_LHC23i1");
    c_Pt_roofit_ratio_LHC23i1->cd(1);

    Pt_grid->Draw();
    TF1 *Pt_LF_HF_norm_binned_LHC23i1 = new TF1("Pt_LF_HF_LHC23i1_BINNED", FuncPt, info.Low_Pt, info.High_Pt, 4);
    Pt_LF_HF_norm_binned_LHC23i1->FixParameter(0, Pt_LF_HF_LHC23i1_BINNED->GetParameter(0));
    Pt_LF_HF_norm_binned_LHC23i1->FixParameter(1, Pt_LF_HF_LHC23i1_BINNED->GetParameter(1));
    Pt_LF_HF_norm_binned_LHC23i1->FixParameter(2, Pt_LF_HF_LHC23i1_BINNED->GetParameter(2));
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        Pt_LF_HF_norm_binned_LHC23i1->FixParameter(3, 0.2);
    else if (info.Low_Mass == 4 && info.High_Mass == 9)
        Pt_LF_HF_norm_binned_LHC23i1->FixParameter(3, 0.216);
    Pt_LF_HF_norm_binned_LHC23i1->Draw("same");
    Pt_LF_HF_LHC23i1_TEST->Draw("same");
    c_Pt_roofit_ratio_LHC23i1->cd(2);

    TH1F *ratio_Pt_LHC23i1 = (TH1F *)(Pt_LF_HF_LHC23i1_TEST->CreateHistogram())->Clone("ratio_Pt_LHC23i1");
    ratio_Pt_LHC23i1->Divide(Pt_LF_HF_norm_binned_LHC23i1->CreateHistogram());
    ratio_Pt_LHC23i1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ratio_Pt_LHC23i1->GetYaxis()->SetRangeUser(0.95, 1.05);
    ratio_Pt_LHC23i1->Draw();

    TF1 *Mass_LF_HF_LHC23i1_BINNED = new TF1("Mass_LF_HF_LHC23i1_BINNED", FuncMass, info.Low_Mass, info.High_Mass, 4);
    Mass_LF_HF_LHC23i1_BINNED->SetParameter(0, Mass_LF_HF_LHC23i1_TEST->GetParameter(0));
    Mass_LF_HF_LHC23i1_BINNED->SetParameter(1, Mass_LF_HF_LHC23i1_TEST->GetParameter(1));
    Mass_LF_HF_LHC23i1_BINNED->SetParameter(2, Mass_LF_HF_LHC23i1_TEST->GetParameter(2));
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        Mass_LF_HF_LHC23i1_BINNED->SetParLimits(3, 0, 650000);
    else if (info.Low_Mass == 4 && info.High_Mass == 9)
        Mass_LF_HF_LHC23i1_BINNED->SetParameter(3, 300000);

    // Hist_Mass_LF_HF_LHC23i1->Scale(1. / Hist_Mass_LF_HF_LHC23i1->GetEntries(), "width");
    Hist_Mass_LF_HF_LHC23i1->Fit(Mass_LF_HF_LHC23i1_BINNED, "LR0I");
    TCanvas *c_Mass_test_LHC23i1 = new TCanvas("c_Mass_test_LHC23i1", "c_Mass_test_LHC23i1", 1200, 1200);
    Hist_Mass_LF_HF_LHC23i1->Draw();
    Mass_LF_HF_LHC23i1_BINNED->Draw("same");
    TH1F *Error_Mass_LHC23i1 = (TH1F *)Hist_Mass_LF_HF_LHC23i1->Clone("Error_Mass_LHC23i1");
    Error_Mass_LHC23i1->SetTitle("Error_Mass_LHC23i1");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(Error_Mass_LHC23i1, 0.95);
    Error_Mass_LHC23i1->SetFillColorAlpha(kBlue, 0.3);
    Error_Mass_LHC23i1->Draw("e3same");
    gPad->BuildLegend();
    TCanvas *c_Mass_roofit_ratio_LCH23i1 = canvas_ratio("c_Mass_roofit_ratio_LCH23i1");
    c_Mass_roofit_ratio_LCH23i1->cd(1);
    TH1F *Mass_grid = (TH1F *)Hist_Mass_LF_HF_LHC23i1->Clone("Mass_grid");
    Mass_grid->Reset();
    Mass_grid->GetYaxis()->SetRangeUser(1.25e-6, 1.25e+1);
    Mass_grid->Draw();
    TF1 *Mass_LF_HF_norm_binned_LHC23i1 = new TF1("Mass_LF_HF_LHC23i1_BINNED", FuncMass, info.Low_Mass, info.High_Mass, 4);
    Mass_LF_HF_norm_binned_LHC23i1->FixParameter(0, Mass_LF_HF_LHC23i1_BINNED->GetParameter(0));
    Mass_LF_HF_norm_binned_LHC23i1->FixParameter(1, Mass_LF_HF_LHC23i1_BINNED->GetParameter(1));
    Mass_LF_HF_norm_binned_LHC23i1->FixParameter(2, Mass_LF_HF_LHC23i1_BINNED->GetParameter(2));
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        Mass_LF_HF_norm_binned_LHC23i1->FixParameter(3, 312);
    else if (info.Low_Mass == 4 && info.High_Mass == 9)
        Mass_LF_HF_norm_binned_LHC23i1->FixParameter(3, 1089);

    Mass_LF_HF_norm_binned_LHC23i1->Draw("same");
    Mass_LF_HF_LHC23i1_TEST->Draw("same");
    c_Mass_roofit_ratio_LCH23i1->cd(2);

    TH1F *ratio_Mass_LHC23i1 = (TH1F *)(Mass_LF_HF_LHC23i1_TEST->CreateHistogram())->Clone("ratio_Mass_LHC23i1");
    ratio_Mass_LHC23i1->Divide(Mass_LF_HF_norm_binned_LHC23i1->CreateHistogram());
    ratio_Mass_LHC23i1->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    ratio_Mass_LHC23i1->GetYaxis()->SetRangeUser(0.85, 1.15);
    ratio_Mass_LHC23i1->Draw();

    // LHC23i2

    TFile *fIn_LHC23i2_Tree = new TFile(Form("%s/output_HF_dimuons/mc_analysis_output/LHC23i2/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i2_MC_output_Tree_merged.root", info.path_to_file.Data()), "READ");

    TTree *fTree_LHC23i2 = (TTree *)fIn_LHC23i2_Tree->Get("DiMuon_Rec_PythiaOnly_LF_HF_Mixed");
    TH1F *Hist_Pt_LF_HF_LHC23i2 = new TH1F("Hist_Pt_LF_HF_LHC23i2", "Hist_Pt_LF_HF_LHC23i2;  #it{p}_{T} (GeV/#it{c});", N_bins_Pt_hist, info.Low_Pt, info.High_Pt);
    TH1F *Hist_Mass_LF_HF_LHC23i2 = new TH1F("Hist_Mass_LF_HF_LHC23i2", "Hist_Mass_LF_HF_LHC23i2;  #it{m}_{#mu#mu} (GeV/#it{c}^{2});", N_bins_Mass_hist, info.Low_Mass, info.High_Mass);

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        fTree_LHC23i2->Draw("pt>>Hist_Pt_LF_HF_LHC23i2", Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
        fTree_LHC23i2->Draw("m>>Hist_Mass_LF_HF_LHC23i2", Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    }
    else
    {
        fTree_LHC23i2->Draw("pt>>Hist_Pt_LF_HF_LHC23i2", Form("((m>%d && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
        fTree_LHC23i2->Draw("m>>Hist_Mass_LF_HF_LHC23i2", Form("((m>%d && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "goff");
    }

    TF1 *Pt_LF_HF_LHC23i2_BINNED = new TF1("Pt_LF_HF_LHC23i2_BINNED", FuncPt, info.Low_Pt, info.High_Pt, 4);
    Pt_LF_HF_LHC23i2_BINNED->SetParameter(0, Pt_LF_HF_LHC23i2_TEST->GetParameter(0));
    Pt_LF_HF_LHC23i2_BINNED->SetParameter(1, Pt_LF_HF_LHC23i2_TEST->GetParameter(1));
    Pt_LF_HF_LHC23i2_BINNED->SetParameter(2, Pt_LF_HF_LHC23i2_TEST->GetParameter(2));
    Pt_LF_HF_LHC23i2_BINNED->SetParLimits(3, 0, 1000);

    Hist_Pt_LF_HF_LHC23i2->Fit(Pt_LF_HF_LHC23i2_BINNED, "LR0");
    TH1F *Error_PT_LHC23i2 = (TH1F *)Hist_Pt_LF_HF_LHC23i2->Clone("Error_PT_LHC23i2");
    Error_PT_LHC23i2->SetTitle("Error_PT_LHC23i2");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(Error_PT_LHC23i2, 0.95);
    // Hist_Pt_LF_HF_LHC23i2->Scale(1. / Hist_Pt_LF_HF_LHC23i2->GetEntries(), "width");
    TCanvas *c_PT_test_LHC23i2 = canvas_noratio("c_PT_test_LHC23i2");
    Hist_Pt_LF_HF_LHC23i2->Draw();
    Pt_LF_HF_LHC23i2_BINNED->Draw("same");
    Error_PT_LHC23i2->SetFillColorAlpha(kBlue, 0.3);
    Error_PT_LHC23i2->Draw("e3same");
    TCanvas *c_Pt_roofit_ratio_LHC23i2 = canvas_ratio("c_Pt_roofit_ratio_LHC23i2");
    c_Pt_roofit_ratio_LHC23i2->cd(1);

    Pt_grid->Draw();
    TF1 *Pt_LF_HF_norm_binned_LHC23i2 = new TF1("Pt_LF_HF_LHC23i2_BINNED", FuncPt, info.Low_Pt, info.High_Pt, 4);
    Pt_LF_HF_norm_binned_LHC23i2->FixParameter(0, Pt_LF_HF_LHC23i2_BINNED->GetParameter(0));
    Pt_LF_HF_norm_binned_LHC23i2->FixParameter(1, Pt_LF_HF_LHC23i2_BINNED->GetParameter(1));
    Pt_LF_HF_norm_binned_LHC23i2->FixParameter(2, Pt_LF_HF_LHC23i2_BINNED->GetParameter(2));
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        Pt_LF_HF_norm_binned_LHC23i2->FixParameter(3, 0.13);
    else if (info.Low_Mass == 4 && info.High_Mass == 9)
        Pt_LF_HF_norm_binned_LHC23i2->FixParameter(3, 0.1352);

    Pt_LF_HF_norm_binned_LHC23i2->Draw("same");
    Pt_LF_HF_LHC23i2_TEST->Draw("same");
    c_Pt_roofit_ratio_LHC23i2->cd(2);

    TH1F *ratio_Pt_LHC23i2 = (TH1F *)(Pt_LF_HF_LHC23i2_TEST->CreateHistogram())->Clone("ratio_Pt_LHC23i2");
    ratio_Pt_LHC23i2->Divide(Pt_LF_HF_norm_binned_LHC23i2->CreateHistogram());
    ratio_Pt_LHC23i2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ratio_Pt_LHC23i2->GetYaxis()->SetRangeUser(0.95, 1.05);
    ratio_Pt_LHC23i2->Draw();

    TF1 *Mass_LF_HF_LHC23i2_BINNED = new TF1("Mass_LF_HF_LHC23i2_BINNED", FuncMass, info.Low_Mass, info.High_Mass, 4);
    Mass_LF_HF_LHC23i2_BINNED->SetParameter(0, Mass_LF_HF_LHC23i2_TEST->GetParameter(0));
    Mass_LF_HF_LHC23i2_BINNED->SetParameter(1, Mass_LF_HF_LHC23i2_TEST->GetParameter(1));
    Mass_LF_HF_LHC23i2_BINNED->SetParameter(2, Mass_LF_HF_LHC23i2_TEST->GetParameter(2));
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        Mass_LF_HF_LHC23i2_BINNED->SetParLimits(3, 0, 650000000);
    else if (info.Low_Mass == 4 && info.High_Mass == 9)
        Mass_LF_HF_LHC23i2_BINNED->SetParameter(3, 300000);
    // Hist_Mass_LF_HF_LHC23i2->Scale(1. / Hist_Mass_LF_HF_LHC23i2->GetEntries(), "width");
    Hist_Mass_LF_HF_LHC23i2->Fit(Mass_LF_HF_LHC23i2_BINNED, "LR0I");

    TCanvas *c_Mass_test_LHC23i2 = canvas_noratio("c_Mass_test_LHC23i2");
    Hist_Mass_LF_HF_LHC23i2->Draw();
    Mass_LF_HF_LHC23i2_BINNED->Draw("same");
    TH1F *Error_Mass_LHC23i2 = (TH1F *)Hist_Mass_LF_HF_LHC23i2->Clone("Error_Mass_LHC23i2");
    Error_Mass_LHC23i2->SetTitle("Error_Mass_LHC23i2");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(Error_Mass_LHC23i2, 0.95);
    Error_Mass_LHC23i2->SetFillColorAlpha(kBlue, 0.3);
    Error_Mass_LHC23i2->Draw("e3same");
    gPad->BuildLegend();
    TCanvas *c_Mass_roofit_ratio_LHC23i2 = canvas_ratio("c_Mass_roofit_ratio_LHC23i2");
    c_Mass_roofit_ratio_LHC23i2->cd(1);
    Mass_grid->Draw();
    TF1 *Mass_LF_HF_norm_binned_LHC23i2 = new TF1("Mass_LF_HF_LHC23i2_BINNED", FuncMass, info.Low_Mass, info.High_Mass, 4);
    Mass_LF_HF_norm_binned_LHC23i2->FixParameter(0, Mass_LF_HF_LHC23i2_BINNED->GetParameter(0));
    Mass_LF_HF_norm_binned_LHC23i2->FixParameter(1, Mass_LF_HF_LHC23i2_BINNED->GetParameter(1));
    Mass_LF_HF_norm_binned_LHC23i2->FixParameter(2, Mass_LF_HF_LHC23i2_BINNED->GetParameter(2));
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        Mass_LF_HF_norm_binned_LHC23i2->FixParameter(3, 2400);
    else if (info.Low_Mass == 4 && info.High_Mass == 9)
        Mass_LF_HF_norm_binned_LHC23i2->FixParameter(3, 300);
    Mass_LF_HF_norm_binned_LHC23i2->Draw("same");
    Mass_LF_HF_LHC23i2_TEST->Draw("same");
    c_Mass_roofit_ratio_LHC23i2->cd(2);

    TH1F *ratio_Mass_LHC23i2 = (TH1F *)(Mass_LF_HF_LHC23i2_TEST->CreateHistogram())->Clone("ratio_Mass_LHC23i2");
    ratio_Mass_LHC23i2->Divide(Mass_LF_HF_norm_binned_LHC23i2->CreateHistogram());
    ratio_Mass_LHC23i2->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    ratio_Mass_LHC23i2->GetYaxis()->SetRangeUser(0.85, 1.15);
    ratio_Mass_LHC23i2->Draw();

    Int_t N_bins_Pt_shape = 999;
    Int_t N_bins_Mass_shape = 999;

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        N_bins_Pt_shape = 150;
        N_bins_Mass_shape = 104;
    }
    else
    {
        N_bins_Pt_shape = 100;
        N_bins_Mass_shape = 50;
    }

    TH1F *Shape_Pt_LF_HF_LHC23i1 = new TH1F("Shape_Pt_LF_HF_LHC23i1", "Shape_Pt_LF_HF_LHC23i1", N_bins_Pt_shape, info.Low_Pt, info.High_Pt);
    TH1F *Shape_Pt_LF_HF_LHC23i2 = new TH1F("Shape_Pt_LF_HF_LHC23i2", "Shape_Pt_LF_HF_LHC23i2", N_bins_Pt_shape, info.Low_Pt, info.High_Pt);
    TH1F *Shape_Pt_LF_HF_mid = new TH1F("Shape_Pt_LF_HF_mid", "Shape_Pt_LF_HF_mid", N_bins_Pt_shape, info.Low_Pt, info.High_Pt);

    TH1F *Shape_Mass_LF_HF_LHC23i1 = new TH1F("Shape_Mass_LF_HF_LHC23i1", "Shape_Mass_LF_HF_LHC23i1", N_bins_Mass_shape, info.Low_Mass, info.High_Mass);
    TH1F *Shape_Mass_LF_HF_LHC23i2 = new TH1F("Shape_Mass_LF_HF_LHC23i2", "Shape_Mass_LF_HF_LHC23i2", N_bins_Mass_shape, info.Low_Mass, info.High_Mass);
    TH1F *Shape_Mass_LF_HF_mid = new TH1F("Shape_Mass_LF_HF_mid", "Shape_Mass_LF_HF_mid", N_bins_Mass_shape, info.Low_Mass, info.High_Mass);

    for (Int_t i_pt_bin = 0; i_pt_bin < Shape_Pt_LF_HF_LHC23i1->GetNbinsX(); i_pt_bin++)
    {
        Shape_Pt_LF_HF_LHC23i1->SetBinContent(i_pt_bin + 1, Pt_LF_HF_LHC23i1_BINNED->Eval(Shape_Pt_LF_HF_LHC23i1->GetBinCenter(i_pt_bin + 1)));
        Shape_Pt_LF_HF_LHC23i1->SetBinError(i_pt_bin + 1, Error_PT_LHC23i1->GetBinError(i_pt_bin + 1));
        Shape_Pt_LF_HF_LHC23i2->SetBinContent(i_pt_bin + 1, Pt_LF_HF_LHC23i2_BINNED->Eval(Shape_Pt_LF_HF_LHC23i2->GetBinCenter(i_pt_bin + 1)));
        Shape_Pt_LF_HF_LHC23i2->SetBinError(i_pt_bin + 1, Error_PT_LHC23i2->GetBinError(i_pt_bin + 1));
    }
    Shape_Pt_LF_HF_LHC23i1->Scale(1. / Shape_Pt_LF_HF_LHC23i1->Integral(), "width");
    Shape_Pt_LF_HF_LHC23i2->Scale(1. / Shape_Pt_LF_HF_LHC23i2->Integral(), "width");

    for (Int_t i_pt_bin = 0; i_pt_bin < Shape_Pt_LF_HF_LHC23i1->GetNbinsX(); i_pt_bin++)
    {
        Shape_Pt_LF_HF_mid->SetBinContent(i_pt_bin + 1, (Shape_Pt_LF_HF_LHC23i1->GetBinContent(i_pt_bin + 1) + Shape_Pt_LF_HF_LHC23i2->GetBinContent(i_pt_bin + 1)) / 2);
        Shape_Pt_LF_HF_mid->SetBinError(i_pt_bin + 1, TMath::Sqrt(0.5 * TMath::Power(Shape_Pt_LF_HF_LHC23i1->GetBinError(i_pt_bin + 1), 2) + 0.5 * TMath::Power(Shape_Pt_LF_HF_LHC23i2->GetBinError(i_pt_bin + 1), 2)));
        // cout << "error LHC23i1: " << Shape_Pt_LF_HF_LHC23i1->GetBinError(i_pt_bin + 1) << endl;
        // cout << "error LHC23i2: " << Shape_Pt_LF_HF_LHC23i2->GetBinError(i_pt_bin + 1) << endl;

        // cout << "error mid: " << Shape_Pt_LF_HF_mid->GetBinError(i_pt_bin + 1) << endl;
    }

    for (Int_t i_mass_bin = 0; i_mass_bin < Shape_Mass_LF_HF_LHC23i1->GetNbinsX(); i_mass_bin++)
    {
        Shape_Mass_LF_HF_LHC23i1->SetBinContent(i_mass_bin + 1, Mass_LF_HF_LHC23i1_BINNED->Eval(Shape_Mass_LF_HF_LHC23i1->GetBinCenter(i_mass_bin + 1)));
        Shape_Mass_LF_HF_LHC23i1->SetBinError(i_mass_bin + 1, Error_Mass_LHC23i1->GetBinError(i_mass_bin + 1));
        Shape_Mass_LF_HF_LHC23i2->SetBinContent(i_mass_bin + 1, Mass_LF_HF_LHC23i2_BINNED->Eval(Shape_Mass_LF_HF_LHC23i2->GetBinCenter(i_mass_bin + 1)));
        Shape_Mass_LF_HF_LHC23i2->SetBinError(i_mass_bin + 1, Error_Mass_LHC23i2->GetBinError(i_mass_bin + 1));
    }

    Shape_Mass_LF_HF_LHC23i1->Scale(1. / Shape_Mass_LF_HF_LHC23i1->Integral(), "width");
    Shape_Mass_LF_HF_LHC23i2->Scale(1. / Shape_Mass_LF_HF_LHC23i2->Integral(), "width");

    for (Int_t i_mass_bin = 0; i_mass_bin < Shape_Mass_LF_HF_LHC23i1->GetNbinsX(); i_mass_bin++)
    {
        Shape_Mass_LF_HF_mid->SetBinContent(i_mass_bin + 1, (Shape_Mass_LF_HF_LHC23i1->GetBinContent(i_mass_bin + 1) + Shape_Mass_LF_HF_LHC23i2->GetBinContent(i_mass_bin + 1)) / 2);
        Shape_Mass_LF_HF_mid->SetBinError(i_mass_bin + 1, TMath::Sqrt(0.5 * TMath::Power(Shape_Mass_LF_HF_LHC23i1->GetBinError(i_mass_bin + 1), 2) + 0.5 * TMath::Power(Shape_Mass_LF_HF_LHC23i2->GetBinError(i_mass_bin + 1), 2)));
        // cout << "error LHC23i1: " << Shape_Mass_LF_HF_LHC23i1->GetBinError(i_mass_bin + 1) << endl;
        // cout << "error LHC23i2: " << Shape_Mass_LF_HF_LHC23i2->GetBinError(i_mass_bin + 1) << endl;
        // cout << "error mid: " << Shape_Mass_LF_HF_mid->GetBinError(i_mass_bin + 1) << endl;
    }
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    TCanvas *c_pt = canvas_noratio("c_pt");
    gPad->SetLogy();
    // Shape_Pt_LF_HF_mid->Scale(1. / Shape_Pt_LF_HF_mid->Integral(), "width");
    Shape_Pt_LF_HF_LHC23i1->SetTitle("from LHC23i1");
    Shape_Pt_LF_HF_LHC23i2->SetTitle("from LHC23i2");
    Shape_Pt_LF_HF_mid->SetTitle("Mixing");
    Shape_Pt_LF_HF_LHC23i1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    Shape_Pt_LF_HF_LHC23i1->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    hist1D_graphic_opt(Shape_Pt_LF_HF_LHC23i1, kFALSE, 1, 24, kMagenta + 2, 1.);
    hist1D_graphic_opt(Shape_Pt_LF_HF_LHC23i2, kFALSE, 1, 24, kGreen + 2, 1.);
    hist1D_graphic_opt(Shape_Pt_LF_HF_mid, kFALSE, 1, 20, kRed - 2, 1.);
    gPad->BuildLegend();

    Shape_Pt_LF_HF_LHC23i1->Draw("PE");
    Shape_Pt_LF_HF_LHC23i2->Draw("PE SAME");
    Shape_Pt_LF_HF_mid->DrawCopy("PE SAME");
    gPad->BuildLegend();
    c_pt->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Pt_shapes_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    TCanvas *c_pt_for_pres = canvas_noratio("c_pt_for_pres");
    gPad->SetLogy();
    Shape_Pt_LF_HF_LHC23i1->DrawCopy("PE");
    Shape_Pt_LF_HF_LHC23i2->DrawCopy("PE SAME");
    gPad->BuildLegend();
    c_pt_for_pres->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Pt_shapes_forPres_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    hist1D_graphic_opt(Shape_Pt_LF_HF_mid, kFALSE, 1, 24, kRed - 2, 1.);
    TF1 *Pt_LF_HF_MID = new TF1("Pt_LF_HF_MID", FuncPt, info.Low_Pt, info.High_Pt, 4);
    Pt_LF_HF_MID->SetParameter(0, Pt_LF_HF_LHC23i2_TEST->GetParameter(0));
    Pt_LF_HF_MID->SetParameter(1, Pt_LF_HF_LHC23i2_TEST->GetParameter(1));
    Pt_LF_HF_MID->SetParameter(2, Pt_LF_HF_LHC23i2_TEST->GetParameter(2));
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        Pt_LF_HF_MID->SetParLimits(3, -10000, 10000);
    else if (info.Low_Mass == 4 && info.High_Mass == 9)
        Pt_LF_HF_MID->SetParameter(3, 1500);

    Shape_Pt_LF_HF_mid->Fit(Pt_LF_HF_MID, "R0");

    // TCanvas *c_PT_mid = canvas_ratio("c_PT_mid");
    // c_PT_mid->cd(1);
    // Shape_Pt_LF_HF_mid->Draw();
    // Pt_LF_HF_MID->Draw("same");

    // TH1F *Shape_Pt_LF_HF_mid_ratio = (TH1F *)Shape_Pt_LF_HF_mid->Clone("ddjs");
    // Shape_Pt_LF_HF_mid_ratio->Divide(Pt_LF_HF_MID);
    // c_PT_mid->cd(2);
    // Shape_Pt_LF_HF_mid_ratio->Draw();
    Shape_Pt_LF_HF_mid->GetYaxis()->SetRangeUser(Shape_Pt_LF_HF_mid->GetMinimum() * 0.2, Shape_Pt_LF_HF_mid->GetMaximum() * 2e+2);
    TString Title = TString::Format("%d < #it{m}_{#mu#mu} < %d (GeV/#it{c}^{2}), PYTHIA only", info.Low_Mass, info.High_Mass);
    TCanvas *c_PT_mid = printMC_ratio_BINNED("c_PT_mid", Title, Shape_Pt_LF_HF_mid, Pt_LF_HF_MID, kBlack);
    c_PT_mid->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Pt_extr_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));

    TCanvas *c_Pt_diff = canvas_noratio("c_Pt_diff");
    c_Pt_diff->cd();
    TH1F *Pt_Diff1 = (TH1F *)Shape_Pt_LF_HF_LHC23i2->Clone("Pt_Diff1");
    Pt_Diff1->Add(Shape_Pt_LF_HF_mid, -1);
    TH1F *Pt_Diff2 = (TH1F *)Shape_Pt_LF_HF_LHC23i1->Clone("Pt_Diff2");
    Pt_Diff2->Add(Shape_Pt_LF_HF_mid, -1);
    Pt_Diff1->GetYaxis()->SetRangeUser(-0.04, 0.04);
    Pt_Diff1->Draw("P hist ");
    Pt_Diff2->Draw("P hist same");
    c_Pt_diff->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Pt_test_diff_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));

    cout << "===MASS====" << endl;
    TCanvas *c_Mass = canvas_noratio("c_Mass");
    gPad->SetLogy();
    // Shape_Pt_LF_HF_mid->Scale(1. / Shape_Pt_LF_HF_mid->Integral(), "width");
    Shape_Mass_LF_HF_LHC23i1->SetTitle("from LHC23i1");
    Shape_Mass_LF_HF_LHC23i2->SetTitle("from LHC23i2");
    Shape_Mass_LF_HF_mid->SetTitle("Mixing");
    Shape_Mass_LF_HF_LHC23i1->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    Shape_Mass_LF_HF_LHC23i1->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
    hist1D_graphic_opt(Shape_Mass_LF_HF_LHC23i1, kFALSE, 1, 24, kMagenta + 2, 1.);
    hist1D_graphic_opt(Shape_Mass_LF_HF_LHC23i2, kFALSE, 1, 24, kGreen + 2, 1.);
    hist1D_graphic_opt(Shape_Mass_LF_HF_mid, kFALSE, 1, 20, kRed - 2, 1.);
    Shape_Mass_LF_HF_LHC23i1->Draw("PE");
    Shape_Mass_LF_HF_LHC23i2->Draw("PE SAME");
    Shape_Mass_LF_HF_mid->DrawCopy("PE SAME");
    gPad->BuildLegend();
    c_Mass->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Mass_shapes_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    TCanvas *c_Mass_for_pres = canvas_noratio("c_Mass_for_pres");
    gPad->SetLogy();
    Shape_Mass_LF_HF_LHC23i1->DrawCopy("PE");
    Shape_Mass_LF_HF_LHC23i2->DrawCopy("PE SAME");
    gPad->BuildLegend();
    c_Mass_for_pres->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Mass_shape_forPres_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));

    hist1D_graphic_opt(Shape_Mass_LF_HF_mid, kFALSE, 1, 24, kRed - 2, 1.);

    // TCanvas *c_Mass = new TCanvas("c_Mass", "c_Mass", 1200, 1200);

    // Shape_Mass_LF_HF_LHC23i1->SetLineColor(kGreen);
    // Shape_Mass_LF_HF_LHC23i2->SetLineColor(kOrange);
    // Shape_Mass_LF_HF_mid->SetLineColor(kBlack);

    // Shape_Mass_LF_HF_LHC23i1->Draw();
    // Shape_Mass_LF_HF_LHC23i2->Draw("SAME");
    // Shape_Mass_LF_HF_mid->Draw("SAME");
    TF1 *Mass_LF_HF_MID = new TF1("Mass_LF_HF_MID", FuncMass, info.Low_Mass, info.High_Mass, 4);
    Mass_LF_HF_MID->SetParameter(0, Mass_LF_HF_LHC23i2_TEST->GetParameter(0));
    Mass_LF_HF_MID->SetParameter(1, Mass_LF_HF_LHC23i2_TEST->GetParameter(1));
    Mass_LF_HF_MID->SetParameter(2, Mass_LF_HF_LHC23i2_TEST->GetParameter(2));

    if (info.Low_Mass == 4 && info.High_Mass == 30)
        Mass_LF_HF_MID->SetParLimits(3, 0, 10000);

    else if (info.Low_Mass == 4 && info.High_Mass == 9)
        Mass_LF_HF_MID->SetParLimits(3, 0, 10000);

    Shape_Mass_LF_HF_mid->Fit(Mass_LF_HF_MID, "LR0");
    // Mass_LF_HF_MID->Draw("same");
    // gPad->BuildLegend();
    Title.Form("#it{p}_{T} < %d (GeV/#it{c}), PYTHIA only", info.High_Pt);
    TCanvas *c_Mass_mid = printMC_ratio_BINNED("c_Mass_mid", Title, Shape_Mass_LF_HF_mid, Mass_LF_HF_MID, kBlack);
    c_Mass_mid->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Mass_extr_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    // TCanvas *c_Mass_mid = canvas_ratio("c_Mass_mid");
    // c_Mass_mid->cd(1);
    // Shape_Mass_LF_HF_mid->Draw();
    // Mass_LF_HF_MID->Draw("same");
    // TH1F *Shape_Mass_LF_HF_mid_ratio = (TH1F *)Shape_Mass_LF_HF_mid->Clone("ddjs");
    // Shape_Mass_LF_HF_mid_ratio->Divide(Mass_LF_HF_MID);
    // c_Mass_mid->cd(2);
    // Shape_Mass_LF_HF_mid_ratio->Draw();

    TCanvas *c_Mass_diff = canvas_noratio("c_Mass_diff");
    c_Mass_diff->cd();
    TH1F *Mass_Diff1 = (TH1F *)Shape_Mass_LF_HF_LHC23i2->Clone("Mass_Diff1");
    Mass_Diff1->Add(Shape_Mass_LF_HF_mid, -1);
    TH1F *Mass_Diff2 = (TH1F *)Shape_Mass_LF_HF_LHC23i1->Clone("Mass_Diff2");
    Mass_Diff2->Add(Shape_Mass_LF_HF_mid, -1);
    Mass_Diff1->GetYaxis()->SetRangeUser(-0.04, 0.04);
    Mass_Diff1->Draw("P hist ");
    Mass_Diff2->Draw("P hist same");

    c_Mass_diff->SaveAs(Form("results/LF_HF_bkg/images/LF_HF_Mass_test_diff_M_%d_%d_Pt_%d_%d.pdf", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));

    TCanvas *test_hist = canvas_noratio("test_hist");
    test_hist->Divide(2, 1);
    test_hist->cd(1);
    TH1F *Param_Pt_LF_HF = new TH1F(Form("Param_Pt_LF_HF_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), Form("Param_Pt_LF_HF_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), 4, 0, 4);
    for (Int_t bin = 0; bin < Param_Pt_LF_HF->GetNbinsX(); bin++)
    {
        Param_Pt_LF_HF->SetBinContent(bin + 1, Pt_LF_HF_MID->GetParameter(bin));
        Param_Pt_LF_HF->SetBinError(bin + 1, Pt_LF_HF_MID->GetParError(bin));
    }

    Param_Pt_LF_HF->Draw();
    test_hist->cd(2);

    TH1F *Param_Mass_LF_HF = new TH1F(Form("Param_Mass_LF_HF_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), Form("Param_Mass_LF_HF_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), 4, 0, 4);
    for (Int_t bin = 0; bin < Param_Mass_LF_HF->GetNbinsX(); bin++)
    {
        Param_Mass_LF_HF->SetBinContent(bin + 1, Mass_LF_HF_MID->GetParameter(bin));
        Param_Mass_LF_HF->SetBinError(bin + 1, Mass_LF_HF_MID->GetParError(bin));
    }

    Param_Mass_LF_HF->Draw();
    // fOut_param->cd(saving_dir);
    // Param_Pt_LF_HF->Write(0, 2, 0);
    cout << "Output file: " << Form("results/LF_HF_bkg/LF_HF_Mixed_bkg_fromPOWHEG_M_%d_%d_Pt_%d_%d.root", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt) << endl;
    TFile *fOut_param = new TFile(Form("results/LF_HF_bkg/LF_HF_Mixed_bkg_fromPOWHEG_M_%d_%d_Pt_%d_%d.root", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "UPDATE");
    Param_Pt_LF_HF->Write(0, 2, 0);
    Param_Mass_LF_HF->Write(0, 2, 0);

    return;
}

void workspace_Mixing_bkg()
{
    opt info;
    gROOT->ProcessLineSync(Form(".x %s/HF_dimuons/PtMassExpPdf.cxx+", info.path_to_file.Data()));

    TFile *fIn_extra = new TFile(Form("results/pdf_extraction/%s_PDF_%s_Mcut_%0.1f_%0.1f.root", info.Generator.Data(), "withLF_HF_LHC23i2", info.LowM_cut, info.HighM_cut), "READ");
    RooWorkspace *w_extra = (RooWorkspace *)fIn_extra->Get(Form("w_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    w_extra->Print();
    TFile *fIn_param = new TFile(Form("results/LF_HF_bkg/LF_HF_Mixed_bkg_fromPOWHEG_M_%d_%d_Pt_%d_%d.root", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "READ");
    fIn_param->ls();

    const Int_t N_Signal = 5;
    const Int_t N_Bkg = 1;

    TString Name_Signal[N_Signal] = {"Beauty", "Charm", "DY", "HF_Mixed", "LF"};
    TString Name_Bkg[N_Bkg] = {"LF_HF_Mixed"};

    TH1D *Param_Pt[N_Signal];
    TH1D *Param_M[N_Signal];

    RooRealVar *m = new RooRealVar("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", info.Low_Mass, info.High_Mass);
    RooRealVar *pt = new RooRealVar("pt", "#it{p}_{T} (GeV/#it{c})", info.Low_Pt, info.High_Pt);

    RooRealVar *B_DimuMass_extra[N_Signal];
    RooRealVar *n1_DimuMass_extra[N_Signal];
    RooRealVar *n2_DimuMass_extra[N_Signal];
    RooRealVar *B_DimuPt_extra[N_Signal];
    RooRealVar *n1_DimuPt_extra[N_Signal];
    RooRealVar *n2_DimuPt_extra[N_Signal];

    RooRealVar *B_DimuMass[N_Bkg];
    RooRealVar *n1_DimuMass[N_Bkg];
    RooRealVar *n2_DimuMass[N_Bkg];
    RooRealVar *B_DimuPt[N_Bkg];
    RooRealVar *n1_DimuPt[N_Bkg];
    RooRealVar *n2_DimuPt[N_Bkg];

    RooWorkspace *w = new RooWorkspace(Form("w_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), Form("w_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    w->factory(Form("pt[%d,%d], #it{p}_{T} (GeV/#it{c})", info.Low_Pt, info.High_Pt));
    w->factory(Form("m[%d,%d], #it{m}_{#mu#mu} (GeV/#it{c}^{2})", info.Low_Mass, info.High_Mass));

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        w->var("m")->setRange("low", 4, 8);
        w->var("m")->setRange("high", 11, 30);
    }

    TString txt_infos = Form("results/pdf_extractio/workspace_info_%s_withLF_HF_MixedBkg.txt", info.Generator.Data());
    std::ofstream out(txt_infos.Data(), std::ios_base::app);
    TString output_fit;
    out << output_fit.Data();
    output_fit.Form("Generator : %s\n", info.Generator.Data());
    out << output_fit.Data();
    output_fit.Form("Mass Range : [%d,%d] || Mass cut : ]%0.1f,%0.1f[\n", info.Low_Mass, info.High_Mass, info.LowM_cut, info.HighM_cut);
    out << output_fit.Data();
    output_fit.Form("Signals\n");
    out << output_fit.Data();

    for (Int_t i_Signal = 0; i_Signal < N_Signal; i_Signal++)
    {
        cout << Name_Signal[i_Signal].Data() << endl;
        B_DimuMass_extra[i_Signal] = w_extra->var(Form("B_DimuMassFrom%s", Name_Signal[i_Signal].Data()));
        B_DimuMass_extra[i_Signal]->setConstant(kTRUE);
        n1_DimuMass_extra[i_Signal] = w_extra->var(Form("n1_DimuMassFrom%s", Name_Signal[i_Signal].Data()));
        n1_DimuMass_extra[i_Signal]->setConstant(kTRUE);
        n2_DimuMass_extra[i_Signal] = w_extra->var(Form("n2_DimuMassFrom%s", Name_Signal[i_Signal].Data()));
        n2_DimuMass_extra[i_Signal]->setConstant(kTRUE);

        B_DimuPt_extra[i_Signal] = w_extra->var(Form("B_DimuPtFrom%s", Name_Signal[i_Signal].Data()));
        B_DimuPt_extra[i_Signal]->setConstant(kTRUE);
        n1_DimuPt_extra[i_Signal] = w_extra->var(Form("n1_DimuPtFrom%s", Name_Signal[i_Signal].Data()));
        n1_DimuPt_extra[i_Signal]->setConstant(kTRUE);
        n2_DimuPt_extra[i_Signal] = w_extra->var(Form("n2_DimuPtFrom%s", Name_Signal[i_Signal].Data()));
        n2_DimuPt_extra[i_Signal]->setConstant(kTRUE);

        printf("B_DimuPt: %0.3e|| n1_DimuPt: %0.3e|| n2_DimuPt: %0.3e\n", B_DimuPt_extra[i_Signal]->getError(), n1_DimuPt_extra[i_Signal]->getError(), n2_DimuPt_extra[i_Signal]->getError());

        cout << Form("PtMassExpPdf::pdfDimuPtFrom%s(pt, B_DimuPtFrom%s[%0.10f], n1_DimuPtFrom%s[%0.10f], n2_DimuPtFrom%s[%0.10f])", Name_Signal[i_Signal].Data(), Name_Signal[i_Signal].Data(), B_DimuPt_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n1_DimuPt_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n2_DimuPt_extra[i_Signal]->getVal()) << endl;

        w->factory(Form("PtMassExpPdf::pdfDimuPtFrom%s(pt, B_DimuPtFrom%s[%0.10f], n1_DimuPtFrom%s[%0.10f], n2_DimuPtFrom%s[%0.10f])", Name_Signal[i_Signal].Data(), Name_Signal[i_Signal].Data(), B_DimuPt_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n1_DimuPt_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n2_DimuPt_extra[i_Signal]->getVal()));

        printf("B_DimuMass: %0.3e|| n1_DimuMass: %0.3e|| n2_DimuMass: %0.3e\n", B_DimuMass_extra[i_Signal]->getError(), n1_DimuMass_extra[i_Signal]->getError(), n2_DimuMass_extra[i_Signal]->getError());
        cout << (Form("PtMassExpPdf::pdfDimuMassFrom%s(m, B_DimuMassFrom%s[%0.10f], n1_DimuMassFrom%s[%0.10f], n2_DimuMassFrom%s[%0.10f])", Name_Signal[i_Signal].Data(), Name_Signal[i_Signal].Data(), B_DimuMass_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n1_DimuMass_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n2_DimuMass_extra[i_Signal]->getVal())) << endl;

        w->factory(Form("PtMassExpPdf::pdfDimuMassFrom%s(m, B_DimuMassFrom%s[%0.10f], n1_DimuMassFrom%s[%0.10f], n2_DimuMassFrom%s[%0.10f])", Name_Signal[i_Signal].Data(), Name_Signal[i_Signal].Data(), B_DimuMass_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n1_DimuMass_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n2_DimuMass_extra[i_Signal]->getVal()));

        output_fit.Form("==============================================\n");
        out << output_fit.Data();
        output_fit.Form("%s\n", Name_Signal[i_Signal].Data());
        out << output_fit.Data();
        output_fit.Form("B_DimuPt: %0.3e|| n1_DimuPt: %0.3e|| n2_DimuPt: %0.3e\n", B_DimuPt_extra[i_Signal]->getVal(), n1_DimuPt_extra[i_Signal]->getVal(), n2_DimuPt_extra[i_Signal]->getVal());
        out << output_fit.Data();
        output_fit.Form("PtMassExpPdf::pdfDimuPtFrom%s(pt, B_DimuPtFrom%s[%0.10f], n1_DimuPtFrom%s[%0.10f], n2_DimuPtFrom%s[%0.10f])\n", Name_Signal[i_Signal].Data(), Name_Signal[i_Signal].Data(), B_DimuPt_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n1_DimuPt_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n2_DimuPt_extra[i_Signal]->getVal());
        out << output_fit.Data();
        output_fit.Form("B_DimuMass: %0.3e|| n1_DimuMass: %0.3e|| n2_DimuMass: %0.3e\n", B_DimuMass_extra[i_Signal]->getVal(), n1_DimuMass_extra[i_Signal]->getVal(), n2_DimuMass_extra[i_Signal]->getVal());
        out << output_fit.Data();
        output_fit.Form("PtMassExpPdf::pdfDimuMassFrom%s(pt, B_DimuMassFrom%s[%0.10f], n1_DimuMassFrom%s[%0.10f], n2_DimuMassFrom%s[%0.10f])\n", Name_Signal[i_Signal].Data(), Name_Signal[i_Signal].Data(), B_DimuMass_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n1_DimuMass_extra[i_Signal]->getVal(), Name_Signal[i_Signal].Data(), n2_DimuMass_extra[i_Signal]->getVal());
        out << output_fit.Data();
    }
    output_fit.Form("Mixed Backgrounds\n");
    out << output_fit.Data();
    for (Int_t i_Bkg = 0; i_Bkg < N_Bkg; i_Bkg++)
    {
        Param_Pt[i_Bkg] = (TH1D *)fIn_param->Get(Form("Param_Pt_LF_HF_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
        Param_Pt[i_Bkg]->Draw();
        Param_M[i_Bkg] = (TH1D *)fIn_param->Get(Form("Param_Mass_LF_HF_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));

        B_DimuMass[i_Bkg] = new RooRealVar(Form("B_DimuMass_from%s", Name_Bkg[i_Bkg].Data()), Form("B_DimuMass_from%s", Name_Bkg[i_Bkg].Data()), (Double_t)Param_M[i_Bkg]->GetBinContent(1));
        n1_DimuMass[i_Bkg] = new RooRealVar(Form("n1_DimuMass_from%s", Name_Bkg[i_Bkg].Data()), Form("n1_DimuMass_from%s", Name_Bkg[i_Bkg].Data()), (Double_t)Param_M[i_Bkg]->GetBinContent(2));
        n2_DimuMass[i_Bkg] = new RooRealVar(Form("n2_DimuMass_from%s", Name_Bkg[i_Bkg].Data()), Form("n2_DimuMass_from%s", Name_Bkg[i_Bkg].Data()), (Double_t)Param_M[i_Bkg]->GetBinContent(3));

        B_DimuMass[i_Bkg]->setError(Param_M[i_Bkg]->GetBinError(1));
        n1_DimuMass[i_Bkg]->setError(Param_M[i_Bkg]->GetBinError(2));
        n2_DimuMass[i_Bkg]->setError(Param_M[i_Bkg]->GetBinError(3));

        B_DimuPt[i_Bkg] = new RooRealVar(Form("B_DimuPt_from%s", Name_Bkg[i_Bkg].Data()), Form("B_DimuPt_from%s", Name_Bkg[i_Bkg].Data()), (Double_t)Param_Pt[i_Bkg]->GetBinContent(1));
        n1_DimuPt[i_Bkg] = new RooRealVar(Form("n1_DimuPt_from%s", Name_Bkg[i_Bkg].Data()), Form("n1_DimuPt_from%s", Name_Bkg[i_Bkg].Data()), (Double_t)Param_Pt[i_Bkg]->GetBinContent(2));
        n2_DimuPt[i_Bkg] = new RooRealVar(Form("n2_DimuPt_from%s", Name_Bkg[i_Bkg].Data()), Form("n2_DimuMass_from%s", Name_Bkg[i_Bkg].Data()), (Double_t)Param_Pt[i_Bkg]->GetBinContent(3));

        B_DimuPt[i_Bkg]->setError(Param_Pt[i_Bkg]->GetBinError(1));
        n1_DimuPt[i_Bkg]->setError(Param_Pt[i_Bkg]->GetBinError(2));
        n2_DimuPt[i_Bkg]->setError(Param_Pt[i_Bkg]->GetBinError(3));

        // B_DimuPt[i]->setConstant(kTRUE);
        // n1_DimuPt[i]->setConstant(kTRUE);
        // n2_DimuPt[i]->setConstant(kTRUE);

        printf("%s\n", Name_Bkg[i_Bkg].Data());
        printf("B_DimuPt: %0.3f n1_DimuPt: %0.3f n2_DimuPt: %0.3f\n", B_DimuPt[i_Bkg]->getVal(), n1_DimuPt[i_Bkg]->getVal(), n2_DimuPt[i_Bkg]->getVal());
        printf("B_DimuMass: %0.3f n1_DimuMass: %0.3f n2_DimuMass: %0.3f\n", B_DimuMass[i_Bkg]->getVal(), n1_DimuMass[i_Bkg]->getVal(), n2_DimuMass[i_Bkg]->getVal());

        cout << Form("PtMassExpPdf::pdfDimuPtFrom%s(pt, B_DimuPtFrom%s[%0.10f], n1_DimuPtFrom%s[%0.10f], n2_DimuPtFrom%s[%0.10f])", Name_Bkg[i_Bkg].Data(), Name_Bkg[i_Bkg].Data(), B_DimuPt[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n1_DimuPt[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n2_DimuPt[i_Bkg]->getVal()) << endl;

        w->factory(Form("PtMassExpPdf::pdfDimuPtFrom%s(pt, B_DimuPtFrom%s[%0.10f], n1_DimuPtFrom%s[%0.10f], n2_DimuPtFrom%s[%0.10f])", Name_Bkg[i_Bkg].Data(), Name_Bkg[i_Bkg].Data(), B_DimuPt[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n1_DimuPt[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n2_DimuPt[i_Bkg]->getVal()));

        cout << (Form("PtMassExpPdf::pdfDimuMassFrom%s(m, B_DimuMassFrom%s[%0.10f], n1_DimuMassFrom%s[%0.10f], n2_DimuMassFrom%s[%0.10f])", Name_Bkg[i_Bkg].Data(), Name_Bkg[i_Bkg].Data(), B_DimuMass[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n1_DimuMass[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n2_DimuMass[i_Bkg]->getVal())) << endl;

        w->factory(Form("PtMassExpPdf::pdfDimuMassFrom%s(m, B_DimuMassFrom%s[%0.10f], n1_DimuMassFrom%s[%0.10f], n2_DimuMassFrom%s[%0.10f])", Name_Bkg[i_Bkg].Data(), Name_Bkg[i_Bkg].Data(), B_DimuMass[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n1_DimuMass[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n2_DimuMass[i_Bkg]->getVal()));

        output_fit.Form("==============================================\n");
        out << output_fit.Data();
        output_fit.Form("%s\n", Name_Bkg[i_Bkg].Data());
        out << output_fit.Data();
        output_fit.Form("B_DimuPt: %0.3e|| n1_DimuPt: %0.3e|| n2_DimuPt: %0.3e\n", B_DimuPt[i_Bkg]->getVal(), n1_DimuPt[i_Bkg]->getVal(), n2_DimuPt[i_Bkg]->getVal());
        out << output_fit.Data();
        output_fit.Form("PtMassExpPdf::pdfDimuPtFrom%s(pt, B_DimuPtFrom%s[%0.10f], n1_DimuPtFrom%s[%0.10f], n2_DimuPtFrom%s[%0.10f])\n", Name_Bkg[i_Bkg].Data(), Name_Bkg[i_Bkg].Data(), B_DimuPt[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n1_DimuPt[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n2_DimuPt[i_Bkg]->getVal());
        out << output_fit.Data();
        output_fit.Form("B_DimuMass: %0.3e|| n1_DimuMass: %0.3e|| n2_DimuMass: %0.3e\n", B_DimuMass[i_Bkg]->getVal(), n1_DimuMass[i_Bkg]->getVal(), n2_DimuMass[i_Bkg]->getVal());
        out << output_fit.Data();
        output_fit.Form("PtMassExpPdf::pdfDimuMassFrom%s(pt, B_DimuMassFrom%s[%0.10f], n1_DimuMassFrom%s[%0.10f], n2_DimuMassFrom%s[%0.10f])\n", Name_Bkg[i_Bkg].Data(), Name_Bkg[i_Bkg].Data(), B_DimuMass[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n1_DimuMass[i_Bkg]->getVal(), Name_Bkg[i_Bkg].Data(), n2_DimuMass[i_Bkg]->getVal());
        out << output_fit.Data();
    }

    output_fit.Form("Backgrounds\n");
    out << output_fit.Data();
    out.close();
    // return;

    w->writeToFile(Form("results/pdf_extraction/%s_PDF_%s_Mcut_%0.1f_%0.1f.root", info.Generator.Data(), "withLF_HF_Mixed_Bkg", info.LowM_cut, info.HighM_cut), kFALSE);
    w->Print();
    // gDirectory->Add(w[p]);

    return;
}

void unbinned_fit_data_sample_singleregion()
{

    opt info;
    gROOT->ProcessLineSync(Form(".x %s/HF_dimuons/PtMassExpPdf.cxx+", info.path_to_file.Data()));

    // TFile *fIn = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/PROVA_LF_HF_mixed_LHC23i2_removing_Y_POWHEG_PDF_%s.root", info.stat_MC.Data()), "READ");
    TFile *fIn;
    RooWorkspace *w;
    cout << info.systematic.Data() << endl;
    if (info.systematic.Contains("scaled_Low2Up") || info.systematic.Contains("scaled_Up2Low"))
    {
        cout << Form("systematic/syst_workspace_%s_M_%d_%d.root", info.Generator.Data(), info.Low_Mass, info.High_Mass) << endl;
        fIn = new TFile(Form("systematic/syst_workspace_%s_M_%d_%d.root", info.Generator.Data(), info.Low_Mass, info.High_Mass), "READ");
        fIn->ls();
        w = (RooWorkspace *)fIn->Get(Form("w_%s", info.systematic.Data()));
    }
    else
    {
        cout << Form("results/pdf_extraction/%s_PDF_%s_Mcut_%0.1f_%0.1f.root", info.Generator.Data(), info.LF_HF.Data(), info.LowM_cut, info.HighM_cut) << endl;
        fIn = new TFile(Form("results/pdf_extraction/%s_PDF_%s_Mcut_%0.1f_%0.1f.root", info.Generator.Data(), info.LF_HF.Data(), info.LowM_cut, info.HighM_cut), "READ");
        w = (RooWorkspace *)fIn->Get(Form("w_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    }

    w->Print("s");
    // RooRealVar *m = new RooRealVar("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", Low_Mass, High_Mass);
    // RooRealVar *pt = new RooRealVar("pt", "#it{p}_{T} (GeV/#it{c})", Low_Pt, High_Pt);

    RooRealVar *m = w->var("m");
    RooRealVar *pt = w->var("pt");
    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        w->var("m")->setRange("low", 4, info.LowM_cut);
        w->var("m")->setRange("high", info.HighM_cut, 30);

        w->var("pt")->setRange("low", 0, 10);
        w->var("pt")->setRange("high", 10, 30);
    }

    m->setBins(info.Mass_Binning);
    pt->setBins(info.Pt_Binning);

    RooCategory sample("sample", "sample");
    sample.defineType("mass");
    sample.defineType("transversemomentum");
    TFile *fIn_data;
    if (info.stat_Data.Contains("Run2"))
        fIn_data = new TFile(Form("%s/HF_dimuons/data_analysis/Tree_MassPt_MassCut4_Run2.root", info.path_to_file.Data()), "READ");
    else
        fIn_data = new TFile(Form("%s/HF_dimuons/data_analysis/3_11_2022/TreeResults_merged.root", info.path_to_file.Data()), "READ");
    fIn_data->ls();
    // Taking data saved in tree
    TTree *tree_data = (TTree *)fIn_data->Get("rec_data_tree");
    // tree_data->Draw("m", Form("((m>%d && m<9) || (m>11 && m<%d)) && (pt > %d && pt <%d)", Low_Mass, High_Mass, Low_Pt, High_Pt));
    gROOT->cd();

    TTree *tree_data_cutted;
    TTree *tree_data_forM_norm = (TTree *)tree_data->CopyTree(Form("pt<%d", info.High_Pt));

    if (info.Low_Mass == 4 && info.High_Mass == 30)
        tree_data_cutted = (TTree *)tree_data->CopyTree(Form("((m>%d && m<%0.1f) ||(m>%0.1f && m<%d)) && pt<%d", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.High_Pt));
    else
        tree_data_cutted = (TTree *)tree_data->CopyTree(Form("(m>%d && m<%d) && (pt <%d)", info.Low_Mass, info.High_Mass, info.High_Pt));

    // tree_data_cutted->Draw("m");
    tree_data_cutted->SetName("tree_data_cutted");
    RooDataSet *unbinned_M_Dimu_data;
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        unbinned_M_Dimu_data = new RooDataSet("M_Dimu_data", "M_Dimu_data", RooArgSet(*m), Import(*tree_data_cutted), Cut(Form("m<%0.1f || (m>%0.1f && m<30)", info.LowM_cut, info.HighM_cut)));
    else
        unbinned_M_Dimu_data = new RooDataSet("M_Dimu_data", "M_Dimu_data", RooArgSet(*m), Import(*tree_data_cutted));
    RooDataSet *unbinned_M_Dimu_data_PLOTTING;
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        unbinned_M_Dimu_data_PLOTTING = new RooDataSet("M_Dimu_data", "M_Dimu_data", RooArgSet(*m), Import(*tree_data_forM_norm));
    else
        unbinned_M_Dimu_data_PLOTTING = new RooDataSet("M_Dimu_data", "M_Dimu_data", RooArgSet(*m), Import(*tree_data_cutted));

    RooDataSet *unbinned_Pt_Dimu_data = new RooDataSet("Pt_Dimu_data", "Pt_Dimu_data", RooArgSet(*pt), Import(*tree_data_cutted), Cut(Form("pt<%d", info.High_Pt)));
    RooDataSet *unbinned_combData_set = new RooDataSet("combData", "combined data", RooArgSet(*m, *pt), Index(sample), Import("mass", *unbinned_M_Dimu_data), Import("transversemomentum", *unbinned_Pt_Dimu_data));

    RooRealVar *B_DimuMass[n_DiMuSelection];
    RooRealVar *n1_DimuMass[n_DiMuSelection];
    RooRealVar *n2_DimuMass[n_DiMuSelection];
    RooRealVar *B_DimuPt[n_DiMuSelection];
    RooRealVar *n1_DimuPt[n_DiMuSelection];
    RooRealVar *n2_DimuPt[n_DiMuSelection];
    RooAbsPdf *pdfDimuMass[n_DiMuSelection];
    RooAbsPdf *pdfDimuPt[n_DiMuSelection];
    for (Int_t i_DiMu_Sel = 0; i_DiMu_Sel < n_DiMuSelection; i_DiMu_Sel++)
    {
        cout << info.Name_DimuSel[i_DiMu_Sel].Data() << endl;
        B_DimuMass[i_DiMu_Sel] = w->var(Form("B_DimuMassFrom%s", info.Name_DimuSel[i_DiMu_Sel].Data()));
        B_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        n1_DimuMass[i_DiMu_Sel] = w->var(Form("n1_DimuMassFrom%s", info.Name_DimuSel[i_DiMu_Sel].Data()));
        n1_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        n2_DimuMass[i_DiMu_Sel] = w->var(Form("n2_DimuMassFrom%s", info.Name_DimuSel[i_DiMu_Sel].Data()));
        n2_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        B_DimuPt[i_DiMu_Sel] = w->var(Form("B_DimuPtFrom%s", info.Name_DimuSel[i_DiMu_Sel].Data()));
        B_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        n1_DimuPt[i_DiMu_Sel] = w->var(Form("n1_DimuPtFrom%s", info.Name_DimuSel[i_DiMu_Sel].Data()));
        n1_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        n2_DimuPt[i_DiMu_Sel] = w->var(Form("n2_DimuPtFrom%s", info.Name_DimuSel[i_DiMu_Sel].Data()));
        n2_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        printf("B_DimuMass: %0.3e|| n1_DimuMass: %0.3e|| n2_DimuMass: %0.3e\n", B_DimuMass[i_DiMu_Sel]->getError(), n1_DimuMass[i_DiMu_Sel]->getError(), n2_DimuMass[i_DiMu_Sel]->getError());

        pdfDimuMass[i_DiMu_Sel] = w->pdf(Form("pdfDimuMassFrom%s", info.Name_DimuSel[i_DiMu_Sel].Data()));
        pdfDimuPt[i_DiMu_Sel] = w->pdf(Form("pdfDimuPtFrom%s", info.Name_DimuSel[i_DiMu_Sel].Data()));
    }

    RooRealVar *normForC = new RooRealVar(Form("n_%s_output", info.Name_DimuSel[0].Data()), "number dimuon from c", 980000, 0, 2000000);
    RooRealVar *normForB = new RooRealVar(Form("n_%s_output", info.Name_DimuSel[1].Data()), "number dimuon from b", 20000, 0, 2000000);
    RooRealVar *normForMixed = new RooRealVar(Form("n_%s_output", info.Name_DimuSel[3].Data()), "number dimuon from b,c", (info.HF_Mixed_fraction / 100) * tree_data_cutted->GetEntries());
    normForMixed->setConstant(kTRUE);
    RooRealVar *normForDY = new RooRealVar(Form("n_%s_output", info.Name_DimuSel[2].Data()), "number dimuon from DY", 380000, 0, 2000000);
    // RooRealVar *normForLF_HF_Mixed = new RooRealVar("n_LF_HF_mixed_output", "number dimuon from LF-HF mixed", 380000, 0, 2000000);
    RooRealVar *normForLF_HF_Mixed = new RooRealVar(Form("n_%s_output", info.Name_DimuSel[4].Data()), "number dimuon from LF-HF mixed", (info.LF_HF_Mixed_fraction / 100) * tree_data_cutted->GetEntries());
    normForLF_HF_Mixed->setConstant(kTRUE);
    RooRealVar *normForLF = new RooRealVar((Form("n_%s_output", info.Name_DimuSel[5].Data())), "number dimuon from b,c", (info.LF_fraction / 100) * tree_data_cutted->GetEntries());
    normForLF->setConstant(kTRUE);

    RooAddPdf *m_model;
    RooAddPdf *pt_model;
    if (info.DY.Contains("noDY") && info.LF_HF.Contains("noLF_HF"))
    {
        m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*pdfDimuMass[0], *pdfDimuMass[1], *pdfDimuMass[3]), RooArgList(*normForC, *normForB, *normForMixed));
        pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*pdfDimuPt[0], *pdfDimuPt[1], *pdfDimuPt[3]), RooArgList(*normForC, *normForB, *normForMixed));
    }
    else if (info.DY.Contains("withDY") && info.LF_HF.Contains("noLF_HF"))
    {
        m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed + n_DY_output*dimuMassFromDY", RooArgList(*pdfDimuMass[0], *pdfDimuMass[1], *pdfDimuMass[3], *pdfDimuMass[2]), RooArgList(*normForC, *normForB, *normForMixed, *normForDY));
        pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed + n_DY_output*dimuPtFromDY", RooArgList(*pdfDimuPt[0], *pdfDimuPt[1], *pdfDimuPt[3], *pdfDimuPt[2]), RooArgList(*normForC, *normForB, *normForMixed, *normForDY));
    }
    else if (info.DY.Contains("noDY") && (info.LF_HF.Contains("withLF_HF_LHC23i2") || info.LF_HF.Contains("withLF_HF_LHC23i1") || info.LF_HF.Contains("withLF_HF_Mixed_Bkg")))
    {
        if (info.LF.Contains("withLF"))
        {
            m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed++ n_DY_output*dimuMassFromDY", RooArgList(*pdfDimuMass[0], *pdfDimuMass[1], *pdfDimuMass[3], *pdfDimuMass[4], *pdfDimuMass[5]), RooArgList(*normForC, *normForB, *normForMixed, *normForLF_HF_Mixed, *normForLF));
            pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed++ DYed_output*dimuPtFromDY", RooArgList(*pdfDimuPt[0], *pdfDimuPt[1], *pdfDimuPt[3], *pdfDimuPt[4], *pdfDimuPt[5]), RooArgList(*normForC, *normForB, *normForMixed, *normForLF_HF_Mixed, *normForLF));
        }
        else if (info.LF.Contains("noLF"))
        {
            m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed++ n_DY_output*dimuMassFromDY", RooArgList(*pdfDimuMass[0], *pdfDimuMass[1], *pdfDimuMass[3], *pdfDimuMass[4]), RooArgList(*normForC, *normForB, *normForMixed, *normForLF_HF_Mixed));
            pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed++ DYed_output*dimuPtFromDY", RooArgList(*pdfDimuPt[0], *pdfDimuPt[1], *pdfDimuPt[3], *pdfDimuPt[4]), RooArgList(*normForC, *normForB, *normForMixed, *normForLF_HF_Mixed));
        }
    }
    else if (info.DY.Contains("withDY") && info.LF_HF.Contains("withLF_HF"))
    {

        // TString Name_DimuSel[n_DiMuSelection] = {"Charm", "Beauty", "DY", "HF_Mixed", "LF_HF_Mixed", "LF"};
        m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed + n_DY_output*dimuMassFromDY + n_LF_HF_mixed_output*dimuMassFromDY", RooArgList(*pdfDimuMass[0], *pdfDimuMass[1], *pdfDimuMass[2], *pdfDimuMass[3], *pdfDimuMass[4], *pdfDimuMass[5]), RooArgList(*normForC, *normForB, *normForDY, *normForMixed, *normForLF_HF_Mixed, *normForLF));
        pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed + n_DY_output*dimuPtFromDY + n_LF_HF_mixed_output*dimuPtFromDY", RooArgList(*pdfDimuPt[0], *pdfDimuPt[1], *pdfDimuPt[2], *pdfDimuPt[3], *pdfDimuPt[4], *pdfDimuPt[5]), RooArgList(*normForC, *normForB, *normForDY, *normForMixed, *normForLF_HF_Mixed, *normForLF));
    }

    RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
    simPdf.addPdf(*m_model, "mass");
    simPdf.addPdf(*pt_model, "transversemomentum");
    simPdf.Print("t");
    RooFitResult *RooFitter;
    if (info.Low_Mass == 4 && info.High_Mass == 30)
        RooFitter = simPdf.fitTo(*unbinned_combData_set, Range("low,high"), SumCoefRange("low,high"), Save(), SumW2Error(true));
    else
        RooFitter = simPdf.fitTo(*unbinned_combData_set, Minimizer("Minuit2"), Save(), SumW2Error(true));

    TH1D *Fit_Result = new TH1D(Form("hist_fit_result_%s_%s", info.DY.Data(), info.LF_HF.Data()), "; coeff x Pt fit", n_DiMuSelection, 0, n_DiMuSelection);
    cout << info.LF_HF.Data() << endl;
    RooFitter->floatParsFinal().Print("s");
    // auto integral_4_30 = pdfDimuMass[2]->createIntegral(*m, RooFit::NormSet(*m), RooFit::Range("low,high"));
    // cout << "integral_4_30->getVal(): " << integral_4_30->getVal() << endl;
    // m->setRange("stupid_range", 15, 30);
    // auto integral = pdfDimuMass[2]->createIntegral(*m, RooFit::NormSet(*m), RooFit::Range("stupid_range"));
    // cout << "integral->getVal(): " << integral->getVal() << endl;
    // cout << "normForC->getVal(): " << normForC->getVal();
    // cout << "normForB->getVal(): " << normForB->getVal();
    // cout << "normForDY->getVal(): " << normForDY->getVal();
    // cout << "normForMixed->getVal(): " << normForMixed->getVal();
    // cout << "normForLF_HF_Mixed->getVal(): " << normForLF_HF_Mixed->getVal();
    // cout << "normForLF->getVal(): " << normForLF->getVal();

    Fit_Result->SetBinContent(1, normForC->getVal());
    Fit_Result->SetBinContent(2, normForB->getVal());
    Fit_Result->SetBinContent(3, normForDY->getVal());
    Fit_Result->SetBinContent(4, normForMixed->getVal());
    Fit_Result->SetBinContent(5, normForLF_HF_Mixed->getVal());
    Fit_Result->SetBinContent(6, normForLF->getVal());

    Fit_Result->SetBinError(1, normForC->getError());
    Fit_Result->SetBinError(2, normForB->getError());
    Fit_Result->SetBinError(3, normForDY->getError());
    Fit_Result->SetBinError(4, TMath::Sqrt(normForMixed->getVal()));
    Fit_Result->SetBinError(5, TMath::Sqrt(normForLF_HF_Mixed->getVal()));
    Fit_Result->SetBinError(6, TMath::Sqrt(normForLF->getVal()));
    // fIn->Close();

    RooPlot *m_frame = m->frame(Title("m_frame"), Name(Form("m_frame_%s_%s", info.DY.Data(), info.LF_HF.Data())));
    m_frame->SetTitle(" ");
    m_frame->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
    m_frame->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");

    RooPlot *pt_frame = pt->frame(Title("pt_frame"), Name(Form("pt_frame_%s_%s", info.DY.Data(), info.LF_HF.Data())));
    pt_frame->SetTitle(" ");
    pt_frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    pt_frame->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

    Double_t Scale_factor = 1.;
    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {

        Scale_factor = (Double_t)tree_data_cutted->GetEntries() / tree_data_forM_norm->GetEntries();
        cout << "Scale_factor: " << Scale_factor << endl;
        unbinned_M_Dimu_data_PLOTTING->plotOn(m_frame, Name("combDatamass"), DrawOption("PEZ"), Rescale((Double_t)tree_data_cutted->GetEntries() / (2 * tree_data_forM_norm->GetEntries())));
        unbinned_combData_set->plotOn(pt_frame, Name("combDatapt"), Cut("sample==sample::transversemomentum"), DrawOption("PEZ"));
    }
    else
    {
        unbinned_combData_set->plotOn(m_frame, Name("combDatamass"), Cut("sample==sample::mass"), DrawOption("PEZ"));
        unbinned_combData_set->plotOn(pt_frame, Name("combDatapt"), Cut("sample==sample::transversemomentum"), DrawOption("PEZ"));
    }

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *unbinned_combData_set), Range("low,high"), NormRange("low,high"), LineStyle(kSolid), LineColor(kRed));
        simPdf.plotOn(pt_frame, Name("pdfpt"), Slice(sample, "transversemomentum"), ProjWData(sample, *unbinned_combData_set), Range("low,high"), NormRange("low,high"), LineStyle(kSolid), LineColor(kRed));
    }
    else
    {
        simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
        simPdf.plotOn(pt_frame, Name("pdfpt"), Slice(sample, "transversemomentum"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
    }

    RooFitter->floatParsFinal().Print("s");

    for (Int_t i_DiMu_Sel = 0; i_DiMu_Sel < n_DiMuSelection; i_DiMu_Sel++)
    {
        if (info.Low_Mass == 4 && info.High_Mass == 30)
        {
            // TString Name_DimuSel[n_DiMuSelection] = {"Charm", "Beauty", "DY", "HF_Mixed", "LF_HF_Mixed", "LF"};
            simPdf.plotOn(m_frame, Name(Form("pdfmass%s", info.Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "mass"), Components(pdfDimuMass[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), Range("low,high"), NormRange("low,high"), LineStyle(kDashed), LineColor(info.color[i_DiMu_Sel]), LineWidth(5));
            simPdf.plotOn(pt_frame, Name(Form("pdfpt%s", info.Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "transversemomentum"), Components(pdfDimuPt[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), Range("low,high"), NormRange("low,high"), LineStyle(kDashed), LineColor(info.color[i_DiMu_Sel]), LineWidth(5));
        }
        else
        {
            simPdf.plotOn(m_frame, Name(Form("pdfmass%s", info.Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "mass"), Components(pdfDimuMass[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(info.color[i_DiMu_Sel]), LineWidth(5));
            simPdf.plotOn(pt_frame, Name(Form("pdfpt%s", info.Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "transversemomentum"), Components(pdfDimuPt[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(info.color[i_DiMu_Sel]), LineWidth(5));
        }
    }
    // RooFitter->floatParsFinal().Print("s");
    // new TCanvas();
    // m_frame->Draw();
    // new TCanvas();
    // pt_frame->Draw();
    // return;
    // convert Roofit to TF1 for ratios
    RooArgSet *m_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*m_model).snapshot(true)); // True means copy the PDF and everything it depends on
    auto &m_modelcopiedPdf = static_cast<RooAbsPdf &>((*m_modelcopyOfEverything)["m_model"]);          // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    RooArgSet *m_modelobs = m_modelcopiedPdf.getObservables(*unbinned_combData_set);
    RooArgSet *m_modelPars = m_modelcopiedPdf.getParameters(*m_modelobs);
    TF1 *m_modelFunc = m_modelcopiedPdf.asTF(*m_modelobs, *m_modelPars, *m);
    m_modelFunc->SetName(Form("m_modelFunc_%s_%s", info.DY.Data(), info.LF_HF.Data()));
    // convert RooDataset to TH1F for ratios
    TH1 *hDimuM_data = unbinned_M_Dimu_data->createHistogram(Form("hDimuM_data_%s_%s", info.DY.Data(), info.LF_HF.Data()), *m, Binning(info.Mass_Binning, info.Low_Mass, info.High_Mass));
    hDimuM_data->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    hDimuM_data->GetYaxis()->SetTitle("d#it{N}/#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

    // convert Roofit to TF1 for ratios
    RooArgSet *pt_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*pt_model).snapshot(true)); // True means copy the PDF and everything it depends on
    auto &pt_modelcopiedPdf = static_cast<RooAbsPdf &>((*pt_modelcopyOfEverything)["pt_model"]);         // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    RooArgSet *pt_modelobs = pt_modelcopiedPdf.getObservables(*unbinned_combData_set);
    RooArgSet *pt_modelPars = pt_modelcopiedPdf.getParameters(*pt_modelobs);
    TF1 *pt_modelFunc = pt_modelcopiedPdf.asTF(*pt_modelobs, *pt_modelPars, *pt);
    pt_modelFunc->SetName(Form("pt_modelFunc_%s_%s", info.DY.Data(), info.LF_HF.Data()));
    // convert RooDataset to TH1F for ratios
    TH1 *hDimuPt_data = unbinned_Pt_Dimu_data->createHistogram(Form("hDimuPt_data_%s_%s", info.DY.Data(), info.LF_HF.Data()), *pt, Binning(info.Pt_Binning, info.Low_Pt, info.High_Pt));
    hDimuPt_data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hDimuPt_data->GetYaxis()->SetTitle("d#it{N}/#it{p}_{T} (GeV/#it{c})^{-1}");
    // Save useful stuff for plotting fit//
    TString txt_infos;
    if (info.systematic.Contains("scaled_Low2Up") || info.systematic.Contains("scaled_Up2Low"))
    {
        txt_infos.Form("systematic/specifics_fit_results_%s_%s_MC_%s_%s.txt", info.systematic.Data(), info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data());
        fIn = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/systematic/Syst_%s_Fit_Data_%s_MC_%s_%s.root", info.systematic.Data(), info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data()), "UPDATE");
    }
    else
    {
        txt_infos.Form("results/fit_result/specifics_fit_results_%s_MC_%s_%s_%s_%s.txt", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data(), info.LF.Data(), info.DY.Data());
        fIn = new TFile(Form("results/fit_result/Result_Fit_Data_%s_MC_%s_%s_%s_%s.root", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data(), info.LF.Data(), info.DY.Data()), "UPDATE");
    }
    fIn->cd();
    if (!fIn->GetDirectory(TString::Format("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut)))
        fIn->mkdir(TString::Format("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut));

    fIn->cd(TString::Format("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut));
    RooFitter->Write(Form("fit_result_%s_%s", info.DY.Data(), info.LF_HF.Data()));
    Fit_Result->Write(0, 2, 0);
    m_frame->Write(0, 2, 0);
    pt_frame->Write(0, 2, 0);
    pt_modelFunc->Write(0, 2, 0);
    hDimuPt_data->Write(0, 2, 0);
    m_modelFunc->Write(0, 2, 0);
    hDimuM_data->Write(0, 2, 0);
    fIn->Close();
    RooFitter->Print("s");

    std::ofstream out(txt_infos.Data(), std::ios_base::app);

    TString output_fit;
    output_fit.Form("==============================================\n");
    out << output_fit.Data();
    output_fit.Form("Generator : %s\n", info.Generator.Data());
    out << output_fit.Data();
    output_fit.Form("Mass Range : [%d,%d] || Mass cut : ]%0.1f,%0.1f[\n", info.Low_Mass, info.High_Mass, info.LowM_cut, info.HighM_cut);
    out << output_fit.Data();
    output_fit.Form("Pt Range : [%d,%d] \n", info.Low_Pt, info.High_Pt);
    out << output_fit.Data();
    output_fit.Form("fOut name : %s\n", fIn->GetName());
    out << output_fit.Data();
    RooFitter->printMultiline(out, 0, true, " ");
    output_fit.Form("OS pair in Data : %lld\n", tree_data_cutted->GetEntries());
    out << output_fit.Data();
    output_fit.Form("Charm number : %0.2f\n", normForC->getVal());
    out << output_fit.Data();
    output_fit.Form("Beauty number : %0.2f\n", normForB->getVal());
    out << output_fit.Data();
    output_fit.Form("HF-mixed number : %0.2f\n", normForMixed->getVal());
    out << output_fit.Data();
    if (info.DY.Contains("withLF_HF"))
    {
        output_fit.Form("LF-HF Mixed fraction : %0.2f\n", normForLF_HF_Mixed->getVal());
        out << output_fit.Data();
    }
    if (info.DY.Contains("withDY"))
    {
        output_fit.Form("LF fraction : %0.2f\n", normForLF->getVal());
        out << output_fit.Data();

        output_fit.Form("DY fraction : %0.2f\n", normForDY->getVal());
        out << output_fit.Data();
    }

    output_fit.Form("Charm fraction : %0.2f\n", 100 * normForC->getVal() / tree_data_cutted->GetEntries());
    out << output_fit.Data();
    output_fit.Form("Beauty fraction : %0.2f\n", 100 * normForB->getVal() / tree_data_cutted->GetEntries());
    out << output_fit.Data();
    output_fit.Form("HF-mixed fraction : %0.2f\n", 100 * normForMixed->getVal() / tree_data_cutted->GetEntries());
    out << output_fit.Data();
    if (info.DY.Contains("withLF_HF"))
    {
        output_fit.Form("LF-HF Mixed fraction : %0.2f\n", 100 * normForLF_HF_Mixed->getVal() / tree_data_cutted->GetEntries());
        out << output_fit.Data();
    }
    if (info.DY.Contains("withDY"))
    {
        output_fit.Form("DY fraction : %0.2f\n", 100 * normForDY->getVal() / tree_data_cutted->GetEntries());
        out << output_fit.Data();
    }
    out.close();
    cout << "Saving in: " << txt_infos.Data() << endl;
    return;

    w->import(simPdf);
    // Save the roodataset as histogram and RooSimultaneous as TF1

    m_frame->SetMaximum(1.2e+6);
    m_frame->SetMinimum(1.2e-2);
    TCanvas *M_Canvas = printRooPlot_ratio(m_frame, kFALSE, RooFitter, 0, "pdfmass", m_modelFunc, hDimuM_data, info.Low_Mass, info.High_Mass, normForMixed->getVal(), normForLF_HF_Mixed->getVal());
    // TCanvas *M_Canvas = printMC_ratio(Form("Mcomp_fit_M_%d_%d_Pt_%d_%d_%s", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.DY.Data()), "Title", m_frame, nullptr, hDimuM_data, m_modelFunc, kBlack, info.Low_Mass, info.High_Mass);
    M_Canvas->SetName(Form("Mcomp_fit_M_%d_%d_Pt_%d_%d_%s", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.DY.Data()));

    pt_frame->SetMaximum(1.2e+6);
    pt_frame->SetMinimum(1.2e-2);
    TCanvas *Pt_Canvas = printRooPlot_ratio(pt_frame, kFALSE, RooFitter, 0, "pdfpt", pt_modelFunc, hDimuPt_data, info.Low_Pt, info.High_Pt, normForMixed->getVal(), normForLF_HF_Mixed->getVal());
    Pt_Canvas->SetName(Form("Ptcomp_fit_M_%d_%d_Pt_%d_%d_%s", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.DY.Data()));

    M_Canvas->SaveAs(Form("prov_im/%s.png", M_Canvas->GetName()));
    Pt_Canvas->SaveAs(Form("prov_im/%s.png", Pt_Canvas->GetName()));
    RooFitter->Print("s");
    if (info.Plot_Likehood)
    {
        RooPlot *m_frame_chisquare[4];
        Int_t m_binnings[4] = {104, 52, 26, 13};
        Double_t chi2[4];
        TCanvas *c = new TCanvas("c", "c", 1200, 1200);
        c->Divide(2, 2);
        for (Int_t i_chi2 = 0; i_chi2 < 4; i_chi2++)
        {
            m_frame_chisquare[i_chi2] = m->frame(Title(TString::Format("m_frame_chisquare_%d", i_chi2)), Binning(m_binnings[i_chi2]));
            unbinned_combData_set->plotOn(m_frame_chisquare[i_chi2], Name(TString::Format("combDatamass_%d", i_chi2)), Cut("sample==sample::mass"), DrawOption("PEZ"));
            simPdf.plotOn(m_frame_chisquare[i_chi2], Name(TString::Format("pdfmass_%d", i_chi2)), Slice(sample, "mass"), ProjWData(sample, *unbinned_combData_set), Range("low,high"), NormRange("low,high"), LineStyle(kSolid), LineColor(kRed));
            chi2[i_chi2] = m_frame_chisquare[i_chi2]->chiSquare(TString::Format("pdfmass_%d", i_chi2), TString::Format("combDatamass_%d", i_chi2), 3);

            cout << "chi2: " << chi2[i_chi2] << endl;
            c->cd(i_chi2 + 1);
            m_frame_chisquare[i_chi2]->Draw();
        }
        TCanvas *LL = new TCanvas("LL", "LL", 1200, 1200);
        LL->Divide(2, 2);
        // ModelConfig model;
        // model.SetWorkspace(*w);
        // model.SetPdf("simPdf");
        // ProfileLikelihoodCalculator pl(*unbinned_combData_set, model);
        // pl.SetConfidenceLevel(0.95); // 95% interval
        // LikelihoodInterval *interval = pl.GetInterval();

        // // print out the interval on the first Parameter of Interest
        // // RooRealVar *firstPOI = (RooRealVar *)model.GetParametersOfInterest()->first();
        // cout << "\n>>>> RESULT : " << 0.95 * 100 << "% interval on " << normForC->GetName() << " is : ["
        //      << interval->LowerLimit(*normForC) << ", " << interval->UpperLimit(*normForC) << "]\n " << endl;

        // // make a plot

        // cout << "making a plot of the profile likelihood function ....(if it is taking a lot of time use less points or the "
        //         "TF1 drawing option)\n";
        // LikelihoodIntervalPlot plot(interval);
        // plot.SetNPoints(10); // do not use too many points, it could become very slow for some models

        // plot.SetRange(0, 80000);
        // TString opt = TString("tf1");
        // LL->cd(1);
        // plot.Draw(opt); // use option TF1 if too slow (plot.Draw("tf1")
        RooAbsReal *nll_charm = simPdf.createNLL(*unbinned_combData_set, NumCPU(6));
        // Minimize likelihood w.r.t all parameters before making plots
        RooMinimizer(*nll_charm).lastMinuitFit();
        // Plot likelihood scan frac
        // RooWorkspace *w_Nll = new RooWorkspace("w_Nll", "workspace");
        // w_Nll->import(*nll) ;
        //
        //
        // w_Nll->writeToFile(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/root_files/likehood_%s.root", mass_range_data.Data()));

        RooPlot *frame_charm = normForC->frame(Bins(1000), Range(55000, 65000), Title(TString::Format("%d < m_{#mu#mu} < %d GeV/#it{c}^{2} && #it{p}_{T} < %d GeV/#it{c}", info.Low_Mass, info.High_Mass, info.High_Pt)));
        nll_charm->plotOn(frame_charm, Name("nll_charm"), ShiftToZero(), LineStyle(kDashed), LineColor(kMagenta + 2), LineWidth(5));
        cout << "nll_charm->getVal()" << nll_charm->getVal() << endl;
        // nll_charm->paramOn(*frame_charm);
        // std::unique_ptr<RooAbsReal> pll_charm{nll_charm->createProfile(*normForC)};
        // pll_charm->plotOn(frame_charm, Name("pll_charm"), ShiftToZero(), LineStyle(kSolid), LineColor(kMagenta + 2), LineWidth(5));
        frame_charm->GetXaxis()->SetTitle("N_{#mu^{#plus}#mu^{#minus}} from fit");
        frame_charm->GetXaxis()->SetTitleOffset(1.5);
        frame_charm->GetYaxis()->SetTitleOffset(1.8);
        // frame_charm->GetYaxis()->SetTitle("Likehood");
        frame_charm->SetMaximum(52500);
        frame_charm->SetMinimum(10e-7);
        LL->cd(1);
        frame_charm->Draw();

        // RooPlot *frame_beauty = normForB->frame(Bins(800), Range(1, 80000), Title(TString::Format("%d < m_{#mu#mu} < %d GeV/#it{c}^{2} && #it{p}_{T} < %d GeV/#it{c}", Low_Mass, High_Mass, High_Pt)));
        // RooAbsReal *nll_beauty = simPdf.createNLL(*unbinned_combData_set, NumCPU(4), Name("nll_beauty"));
        // // Minimize likelihood w.r.t all parameters before making plots
        // RooMinimizer(*nll_beauty).lastMinuitFit();
        // nll_beauty->plotOn(frame_beauty, Name("nll_beauty"), ShiftToZero(), LineStyle(kDashed), LineColor(kGreen + 2), LineWidth(5));
        // std::unique_ptr<RooAbsReal> pll_beauty{nll_beauty->createProfile(*normForB)};
        // pll_beauty->plotOn(frame_beauty, Name("pll_beauty"), ShiftToZero(), LineStyle(kSolid), LineColor(kGreen + 2), LineWidth(5));
        // LL->cd(2);
        // frame_beauty->Draw();

        if (info.DY.Contains("withDY"))
        {
            RooAbsReal *nll_DY = simPdf.createNLL(*unbinned_combData_set, NumCPU(4), Name("nll_DY"));
            // Minimize likelihood w.r.t all parameters before making plots
            RooMinimizer(*nll_DY).lastMinuitFit();
            RooPlot *frame_DY = normForDY->frame(Bins(800), Range(1, 80000), Title(TString::Format("%d < m_{#mu#mu} < %d GeV/#it{c}^{2} && #it{p}_{T} < %d GeV/#it{c}", info.Low_Mass, info.High_Mass, info.High_Pt)));
            nll_DY->plotOn(frame_DY, Name("nll_DY"), ShiftToZero(), LineStyle(kDashed), LineColor(kBlack), LineWidth(5));
            LL->cd(4);
            frame_DY->Draw();
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
        // legend->AddEntry("nll_beauty", "Likelihood #it{N}_{#mu^{#plus}#mu^{#minus} #leftarrow b}", "L");
        if (info.DY.Contains("withDY"))
            legend->AddEntry("nll_beauty", "Likelihood #it{N}_{#mu^{#plus}#mu^{#minus} #leftarrow DY}", "L");
        LL->cd(3);
        legend->Draw();
        // if (With_DY)
        // {

        //     canvas->SaveAs(Form("images/LLcanvas%s_withDY.pdf", mass_range_data.Data()));
        //     canvas->SaveAs(Form("images/LLcanvas%s_withDY.png", mass_range_data.Data()));
        // }
        // else
        // {
        //     canvas->SaveAs(Form("images/LLcanvas%s_noDY.pdf", mass_range_data.Data()));
        //     canvas->SaveAs(Form("images/LLcanvas%s_noDY.png", mass_range_data.Data()));
        // }

        // gStyle->SetOptStat(0);
        // gStyle->SetPalette(1);
    }
}

void plotting_fit()
{
    opt info;
    gROOT->ProcessLineSync(Form(".x %s/HF_dimuons/PtMassExpPdf.cxx+", info.path_to_file.Data()));

    TFile *fIn_data;
    if (info.stat_Data.Contains("Run2"))
        fIn_data = new TFile(Form("%s/HF_dimuons/data_analysis/Tree_MassPt_MassCut4_Run2.root", info.path_to_file.Data()), "READ");
    else
        fIn_data = new TFile(Form("%s/HF_dimuons/data_analysis/3_11_2022/TreeResults_merged.root", info.path_to_file.Data()), "READ");
    fIn_data->ls();

    // Taking data saved in tree
    TTree *tree_data = (TTree *)fIn_data->Get("rec_data_tree");
    gROOT->cd();
    TTree *tree_data_forM_norm = (TTree *)tree_data->CopyTree(Form("pt<%d", info.High_Pt));
    // if (info.Low_Mass == 4 && info.High_Mass == 30)

    TTree *tree_data_cutted;
    RooRealVar *normForMixed;
    RooRealVar *normForLF_HF_Mixed;
    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {
        tree_data_cutted = (TTree *)tree_data->CopyTree(Form("((m>%d && m<%0.1f) ||(m>%0.1f && m<%d)) && pt<%d", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.High_Pt));
        normForMixed = new RooRealVar(Form("n_%s_output", info.Name_DimuSel[2].Data()), "number dimuon from b,c", (info.HF_Mixed_fraction / 100) * tree_data_cutted->GetEntries());
        normForLF_HF_Mixed = new RooRealVar(Form("n_%s_output", info.Name_DimuSel[4].Data()), "number dimuon from LF-HF mixed", (info.LF_HF_Mixed_fraction / 100) * tree_data_cutted->GetEntries());
    }
    else
    {

        normForMixed = new RooRealVar(Form("n_%s_output", info.Name_DimuSel[2].Data()), "number dimuon from b,c", (info.HF_Mixed_fraction / 100) * tree_data->GetEntries());
        normForLF_HF_Mixed = new RooRealVar(Form("n_%s_output", info.Name_DimuSel[4].Data()), "number dimuon from LF-HF mixed", (info.LF_HF_Mixed_fraction / 100) * tree_data->GetEntries());
    }

    normForMixed->setConstant(kTRUE);
    normForLF_HF_Mixed->setConstant(kTRUE);
    TFile *fIn;

    if (info.systematic.Contains("scaled_Low2Up") || info.systematic.Contains("scaled_Up2Low"))
    {
        // txt_infos.Form("systematic/specifics_fit_results_%s_%s_MC_%s_%s.txt", info.systematic.Data(), info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data());
        fIn = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/systematic/Syst_%s_Fit_Data_%s_MC_%s_%s.root", info.systematic.Data(), info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data()), "UPDATE");
    }
    else
    {
        // txt_infos.Form("results/fit_result/specifics_fit_results_%s_MC_%s_%s.txt", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data());
        fIn = new TFile(Form("results/fit_result/Result_Fit_Data_%s_MC_%s_%s_%s_%s.root", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data(), info.LF.Data(), info.DY.Data()), "READ");
        // fIn = new TFile(Form("results/fit_result/Result_Fit_Data_%s_MC_%s_%s.root", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data()), "UPDATE");
    }

    fIn->cd(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut));
    fIn->ls();
    RooFitResult *result = (RooFitResult *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/fit_result_%s_%s", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut, info.DY.Data(), info.LF_HF.Data()));

    RooPlot *m_frame = (RooPlot *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/m_frame_%s_%s", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut, info.DY.Data(), info.LF_HF.Data()));
    TH1F *hDimuM_data = (TH1F *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/hDimuM_data_%s_%s__m", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut, info.DY.Data(), info.LF_HF.Data()));
    TF1 *m_modelFunc = (TF1 *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/m_modelFunc_%s_%s", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut, info.DY.Data(), info.LF_HF.Data()));

    cout << hDimuM_data->GetMinimum() << endl;
    cout << hDimuM_data->GetMaximum() << endl;
    if (info.Low_Mass == 4 && info.High_Mass == 9)
        m_frame->GetYaxis()->SetRangeUser(1.2, hDimuM_data->GetMaximum() * 1.2e+2);
    else
        m_frame->GetYaxis()->SetRangeUser(1.2e-3, hDimuM_data->GetMaximum() * 3.2e+1);
    // m_frame->GetYaxis()->SetRangeUser(1.2, 6000);

    Double_t Integral = 0;

    for (Int_t i = 0; i < hDimuM_data->GetNbinsX(); i++)
    {
        // cout<<"Bin center: "<<hDimuM_data->GetBinCenter(i+1)<<" Content: "<<hDimuM_data->GetBinContent(i+1)/hDimuM_data->Integral(1,hDimuM_data->GetNbinsX()+1,"width")<<" fit value: "<<m_modelFunc->Eval(hDimuM_data->GetBinCenter(i+1))<<" ratio: "<<(hDimuM_data->GetBinContent(i+1)/hDimuM_data->Integral(1,hDimuM_data->GetNbinsX()+1,"width"))/m_modelFunc->Eval(hDimuM_data->GetBinCenter(i+1))<<endl;
        Integral = Integral + hDimuM_data->GetBinContent(i + 1);
    }
    cout << "Integral: " << tree_data_forM_norm->GetEntries() << endl;
    cout << "hDimuM_data->Intagral(): " << hDimuM_data->Integral() << endl;
    if (info.Low_Mass == 4 && info.High_Mass == 9)
        hDimuM_data->Scale(1. / hDimuM_data->Integral(), "width");
    else
        hDimuM_data->Scale(1. / tree_data_forM_norm->GetEntries(), "width");

    TCanvas *M_Canvas = printRooPlot_ratio(m_frame, kFALSE, result, 0, "pdfmass", m_modelFunc, hDimuM_data, info.Low_Mass, info.High_Mass, normForMixed->getVal(), normForLF_HF_Mixed->getVal());
    M_Canvas->SetName("M_Canvas");
    if (info.systematic.Contains("scaled_Low2Up") || info.systematic.Contains("scaled_Up2Low"))
        M_Canvas->SaveAs(Form("systematic/images/M_fit_%s_MC_%s_%s_M_%d_%d_Mcut_%0.1f_%0.1f_%s.pdf", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data(), info.Low_Mass, info.High_Mass, info.LowM_cut, info.HighM_cut, info.systematic.Data()));
    else
        M_Canvas->SaveAs(Form("results/fit_result/images/M_fit_%s_MC_%s_%s_%s_%s_M_%d_%d_Mcut_%0.1f_%0.1f.pdf", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data(), info.LF.Data(), info.DY.Data(), info.Low_Mass, info.High_Mass, info.LowM_cut, info.HighM_cut));

    RooPlot *pt_frame = (RooPlot *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/pt_frame_%s_%s", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut, info.DY.Data(), info.LF_HF.Data()));
    TH1F *hDimuPt_data = (TH1F *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/hDimuPt_data_%s_%s__pt", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut, info.DY.Data(), info.LF_HF.Data()));
    TF1 *pt_modelFunc = (TF1 *)fIn->Get(Form("dir_M_%d_%d_Pt_%d_%d_Mcut_%0.1f_%0.1f/pt_modelFunc_%s_%s", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt, info.LowM_cut, info.HighM_cut, info.DY.Data(), info.LF_HF.Data()));
    if (info.Low_Mass == 4 && info.High_Mass == 9)
        pt_frame->GetYaxis()->SetRangeUser(1.2, hDimuPt_data->GetMaximum() * 1.2e+2);
    else
        pt_frame->GetYaxis()->SetRangeUser(8.5e-2, hDimuPt_data->GetMaximum() * 4.2e+2);
    // pt_frame->GetYaxis()->SetRangeUser(1.2, 4000);

    hDimuPt_data->Scale(1. / hDimuPt_data->Integral(), "width");
    TCanvas *Pt_Canvas = printRooPlot_ratio(pt_frame, kFALSE, result, 0, "pdfpt", pt_modelFunc, hDimuPt_data, info.Low_Pt, info.High_Pt, normForMixed->getVal(), normForLF_HF_Mixed->getVal());
    Pt_Canvas->SetName("Pt_Canvas");
    if (info.systematic.Contains("scaled_Low2Up") || info.systematic.Contains("scaled_Up2Low"))
    {
        // txt_infos.Form("systematic/specifics_fit_results_%s_%s_MC_%s_%s.txt", info.systematic.Data(), info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data());
        Pt_Canvas->SaveAs(Form("systematic/images/Pt_fit_%s_MC_%s_%s_M_%d_%d_Mcut_%0.1f_%0.1f_%s.pdf", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data(), info.Low_Mass, info.High_Mass, info.LowM_cut, info.HighM_cut, info.systematic.Data()));
    }
    else
    {
        // txt_infos.Form("results/fit_result/specifics_fit_results_%s_MC_%s_%s.txt", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data());
        Pt_Canvas->SaveAs(Form("results/fit_result/images/Pt_fit_%s_MC_%s_%s_%s_%s_M_%d_%d_Mcut_%0.1f_%0.1f.pdf", info.stat_Data.Data(), info.Generator.Data(), info.LF_HF.Data(), info.LF.Data(), info.DY.Data(), info.Low_Mass, info.High_Mass, info.LowM_cut, info.HighM_cut));
    }
}

void pdf_extraction()
{
    opt info;
    TString fOut_filename;
    fOut_filename.Form("%s/HF_dimuons/fit_data/results/pdf_extraction/%s_PDF_%s_Mcut_%0.1f_%0.1f.root", info.path_to_file.Data(), info.Generator.Data(), info.LF_HF.Data(), info.LowM_cut, info.HighM_cut);
    std::ofstream out(TString::Format("results/pdf_extraction/specific_extrc_%s_PDF_%s_Mcut_%0.1f_%0.1f.txt", info.Generator.Data(), info.LF_HF.Data(), info.LowM_cut, info.HighM_cut), std::ios_base::app);

    TString output_fit;
    output_fit.Form("==============================================\n");
    out << output_fit.Data();
    output_fit.Form("Generator : %s\n", info.Generator.Data());
    out << output_fit.Data();
    output_fit.Form("Mass Range : [%d,%d] || Mass cut : ]%0.1f,%0.1f[\n", info.Low_Mass, info.High_Mass, info.LowM_cut, info.HighM_cut);
    out << output_fit.Data();
    output_fit.Form("Pt Range : [%d,%d] \n", info.Low_Pt, info.High_Pt);
    out << output_fit.Data();
    output_fit.Form("fOut name : %s\n", fOut_filename.Data());
    out << output_fit.Data();

    gROOT->ProcessLineSync(Form(".x %s/HF_dimuons/PtMassExpPdf.cxx+", info.path_to_file.Data()));

    const Int_t n_DiMuSelection = 7;
    TString Name_DimuSel[n_DiMuSelection] = {"Charm", "Beauty", "DY", "HF_Mixed", "LF_HF_Mixed", "LF", "Data"};
    TString Dir_name = TString::Format("%s/output_HF_dimuons/mc_analysis_output", info.path_to_file.Data());
    TString Version_ALIAOD;
    TString Tree_name[n_DiMuSelection];
    TString FileName[n_DiMuSelection];
    TTree *fTree[n_DiMuSelection];
    TFile *fIn[n_DiMuSelection];
    if (info.stat_MC.Contains("small_stat"))
    {
        Version_ALIAOD = "Version3_AliAOD/save_output";
        // Tree_name = "DiMuon_Rec_FullMass_PowhegOnly";
    }
    else if (info.stat_MC.Contains("full_stat"))
    {
        Version_ALIAOD = "Version_5_AliAOD_skimmed_fwd_fullstat";
        // Tree_name = "DiMuon_Rec_PowhegOnly_";
    }
    for (Int_t i_DiMu = 0; i_DiMu < n_DiMuSelection - 1; i_DiMu++)
    {
        if (Name_DimuSel[i_DiMu].Contains("Charm"))
        {
            if (info.Generator.Contains("POWHEG"))
            {
                FileName[i_DiMu].Form("%s/LHC23i1/%s/LHC23i1_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name[i_DiMu].Form("DiMuon_Rec_PowhegOnly");
            }
            else if (info.Generator.Contains("PYTHIA"))
            {
                FileName[i_DiMu].Form("%s/LHC22b3/%s/LHC22b3_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name[i_DiMu].Form("DiMuon_Rec");
            }
        }

        else if (Name_DimuSel[i_DiMu].Contains("Beauty"))
        {
            if (info.Generator.Contains("POWHEG"))
            {
                FileName[i_DiMu].Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name[i_DiMu].Form("DiMuon_Rec_PowhegOnly");
            }
            else if (info.Generator.Contains("PYTHIA"))
            {
                FileName[i_DiMu].Form("%s/LHC22b3/%s/LHC22b3_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name[i_DiMu].Form("DiMuon_Rec");
            }
        }
        else if (Name_DimuSel[i_DiMu].Contains("LF_HF_Mixed"))
        {
            if (info.LF_HF.Contains("LHC23i1"))
            {
                FileName[i_DiMu].Form("%s/LHC23i1/%s/LHC23i1_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name[i_DiMu].Form("DiMuon_Rec_PythiaOnly");
            }
            else if (info.LF_HF.Contains("LHC23i2"))
            {
                FileName[i_DiMu].Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name[i_DiMu].Form("DiMuon_Rec_PythiaOnly");
            }
        }
        else if (Name_DimuSel[i_DiMu].Contains("HF_Mixed"))
        {
            if (info.Generator.Contains("POWHEG"))
            {
                FileName[i_DiMu].Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name[i_DiMu].Form("DiMuon_Rec_PythiaOnly");
            }
            else if (info.Generator.Contains("PYTHIA"))
            {
                FileName[i_DiMu].Form("%s/LHC22b3/%s/LHC22b3_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
                Tree_name[i_DiMu].Form("DiMuon_Rec");
            }
        }

        else if (Name_DimuSel[i_DiMu].Contains("DY"))
        {
            FileName[i_DiMu].Form("%s/LHC18p_DY/Version_2_AliAOD/LHC18p_DY_MC_output_Tree_merged.root", Dir_name.Data());
            Tree_name[i_DiMu].Form("DiMuon_Rec");
        }

        else if (Name_DimuSel[i_DiMu].Contains("LF"))
        {
            FileName[i_DiMu].Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data());
            Tree_name[i_DiMu].Form("DiMuon_Rec_PythiaOnly");
        }
        fIn[i_DiMu] = new TFile(FileName[i_DiMu], "READ");
        // fIn[i_DiMu]->ls();
        fTree[i_DiMu] = (TTree *)fIn[i_DiMu]->Get(Form("%s_%s", Tree_name[i_DiMu].Data(), Name_DimuSel[i_DiMu].Data()));
        cout << "Component: " << Name_DimuSel[i_DiMu].Data() << " ||file: " << FileName[i_DiMu].Data() << " || tree name: " << Tree_name[i_DiMu].Data() << endl;
        // fTree[i_DiMu]->Print();
        output_fit = Name_DimuSel[i_DiMu];
        out << Form("%s \n", output_fit.Data());
        output_fit = FileName[i_DiMu];
        out << Form("%s \n", output_fit.Data());
        output_fit = Tree_name[i_DiMu];
        out << Form("%s_%s \n", Tree_name[i_DiMu].Data(), Name_DimuSel[i_DiMu].Data());
    }
    gROOT->cd();
    // TFile *fIn[n_DiMuSelection] = {new TFile(Form("%s/LHC23i1/%s/LHC23i1_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())),new TFile(Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())),new TFile(Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())), new TFile(Form("%s/LHC18p_DY/Version_2_AliAOD/LHC18p_DY_MC_output_Tree_merged.root", Dir_name.Data())), new TFile(Form("%s/LHC23i1/%s/LHC23i1_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())), new TFile(Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())), new TFile("/home/michele_pennisi/cernbox/HF_dimuons/data_analysis/Tree_MassPt_MassCut4_Run2.root", "READ")};
    // fIn[0]->ls();
    // TTree *fTree[n_DiMuSelection] = {(TTree *)fIn[0]->Get(Form("%sCharm", Tree_name.Data())), (TTree *)fIn[1]->Get(Form("%sBeauty", Tree_name.Data())), (TTree *)fIn[2]->Get("DiMuon_Rec_PythiaOnly_LF_HF_Mixed"), (TTree *)fIn[3]->Get("DiMuon_Rec_DY"), (TTree *)fIn[4]->Get("DiMuon_Rec_PythiaOnly_LF_HF_Mixed"), (TTree *)fIn[5]->Get("DiMuon_Rec_PythiaOnly_LF"), (TTree *)fIn[6]->Get("rec_data_tree")};

    TTree *fTree_x_pt[n_DiMuSelection];
    TTree *fTree_x_m[n_DiMuSelection];
    TFile *fOut = new TFile(fOut_filename, "UPDATE");
    // TFile *fOut = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/PROVA_LF_HF_mixed_LHC23i1_removing_Y_POWHEG_PDF_%s.root", info.stat_MC.Data()), "UPDATE");
    // TFile *fOut = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/POWHEG_PDF_%s.root", info.stat_MC.Data()), "UPDATE");
    RooDataSet *M_Dimu[n_DiMuSelection];
    RooDataSet *Pt_Dimu[n_DiMuSelection];

    TH1 *m_histo[n_DiMuSelection];
    TH1 *pt_histo[n_DiMuSelection];

    RooPlot *frameDimuMass[n_DiMuSelection];
    RooPlot *frameDimuPt[n_DiMuSelection];

    RooRealVar *B_DimuMass[n_DiMuSelection];
    RooRealVar *n1_DimuMass[n_DiMuSelection];
    RooRealVar *n2_DimuMass[n_DiMuSelection];
    RooRealVar *B_DimuPt[n_DiMuSelection];
    RooRealVar *n1_DimuPt[n_DiMuSelection];
    RooRealVar *n2_DimuPt[n_DiMuSelection];

    TH1D *Param_Pt_unbinned[n_DiMuSelection];
    TH1D *Param_M_unbinned[n_DiMuSelection];

    // RooRealVar *m = new RooRealVar("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", info.Low_Mass, info.High_Mass);
    // RooRealVar *pt = new RooRealVar("pt", "#it{p}_{T} (GeV/#it{c})", info.Low_Pt, info.High_Pt);

    Color_t color[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kOrange + 7, kAzure + 9, kBlue, kRed, kBlack};
    Color_t fillcolor[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kOrange + 7, kAzure + 9, kBlue, kRed, kBlack};

    RooWorkspace *w = new RooWorkspace(Form("w_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "workspace");
    w->factory(Form("pt[%d,%d], #it{p}_{T} (GeV/#it{c})", info.Low_Pt, info.High_Pt));
    RooRealVar *pt = (RooRealVar *)w->var("pt");
    pt->setBins(info.Pt_Binning);
    w->factory(Form("m[%d,%d], #it{m}_{#mu#mu} (GeV/#it{c}^{2})", info.Low_Mass, info.High_Mass));
    RooRealVar *m = w->var("m");
    m->setBins(info.Mass_Binning);

    if (info.Low_Mass == 4 && info.High_Mass == 30)
    {

        w->var("m")->setRange("low", 4, 8);
        w->var("m")->setRange("high", 11, 30);
    }
    m->Print();
    TString Mass_Factory_Info[n_DiMuSelection];
    Mass_Factory_Info[0].Form("PtMassExpPdf::pdfDimuMassFromCharm(m, B_DimuMassFromCharm[2.85,1e-4,1e+4], n1_DimuMassFromCharm[2.81,1e-4,1e+4], n2_DimuMassFromCharm[5,1e-3,100])");
    Mass_Factory_Info[1].Form("PtMassExpPdf::pdfDimuMassFromBeauty(m, B_DimuMassFromBeauty[2.85,1e-4,1e+4], n1_DimuMassFromBeauty[2.256,1e-4,1e+4], n2_DimuMassFromBeauty[2.83,1e-3,100])");
    Mass_Factory_Info[2].Form("PtMassExpPdf::pdfDimuMassFromDY(m, B_DimuMassFromDY[2.65,1e-4,1e+4], n1_DimuMassFromDY[2.81,1e-4,1e+4], n2_DimuMassFromDY[5,1e-4,1e+4])");
    Mass_Factory_Info[3].Form("PtMassExpPdf::pdfDimuMassFromHF_Mixed(m, B_DimuMassFromHF_Mixed[2.65,1e-4,1e+4], n1_DimuMassFromHF_Mixed[2.81,1e-4,1e+4], n2_DimuMassFromHF_Mixed[5,1e-3,100])");
    Mass_Factory_Info[4].Form("PtMassExpPdf::pdfDimuMassFromLF_HF_Mixed(m, B_DimuMassFromLF_HF_Mixed[2.65,1e-4,1e+4], n1_DimuMassFromLF_HF_Mixed[2.81,1e-4,1e+4], n2_DimuMassFromLF_HF_Mixed[5,1e-4,1e+4])");
    Mass_Factory_Info[5].Form("PtMassExpPdf::pdfDimuMassFromLF(m, B_DimuMassFromLF[2.65,1e-4,1e+4], n1_DimuMassFromLF[2.81,1e-4,1e+4], n2_DimuMassFromLF[5,1e-4,1e+4])");
    Mass_Factory_Info[6].Form("PtMassExpPdf::pdfDimuMassFromData(m, B_DimuMassFromData[2.65,1e-4,1e+4], n1_DimuMassFromData[2.81,1e-4,1e+4], n2_DimuMassFromData[5,1e-4,1e+4])");

    TString Pt_Factory_Info[n_DiMuSelection];
    Pt_Factory_Info[0].Form("PtMassExpPdf::pdfDimuPtFromCharm(pt, B_DimuPtFromCharm[2.85,1e-4,1e+4], n1_DimuPtFromCharm[2.81,1e-4,1e+4], n2_DimuPtFromCharm[2.43,1e-4,1e+4])");
    Pt_Factory_Info[1].Form("PtMassExpPdf::pdfDimuPtFromBeauty(pt, B_DimuPtFromBeauty[2.85,1e-4,1e+4], n1_DimuPtFromBeauty[2.81,1e-4,1e+4], n2_DimuPtFromBeauty[2.43,1e-4,1e+4])");
    Pt_Factory_Info[2].Form("PtMassExpPdf::pdfDimuPtFromDY(pt, B_DimuPtFromDY[2.85,1e-4,1e+4], n1_DimuPtFromDY[2.81,1e-4,1e+4], n2_DimuPtFromDY[2.43,1e-4,1e+4])");
    Pt_Factory_Info[3].Form("PtMassExpPdf::pdfDimuPtFromHF_Mixed(pt, B_DimuPtFromHF_Mixed[2.85,1e-4,1e+4], n1_DimuPtFromHF_Mixed[2.81,1e-4,1e+4], n2_DimuPtFromHF_Mixed[2.43,1e-4,1e+4])");
    Pt_Factory_Info[4].Form("PtMassExpPdf::pdfDimuPtFromLF_HF_Mixed(pt, B_DimuPtFromLF_HF_Mixed[2.85,1e-4,1e+4], n1_DimuPtFromLF_HF_Mixed[2.81,1e-4,1e+4], n2_DimuPtFromLF_HF_Mixed[2.43,1e-4,1e+4])");
    Pt_Factory_Info[5].Form("PtMassExpPdf::pdfDimuPtFromLF(pt, B_DimuPtFromLF[2.85,1e-4,1e+4], n1_DimuPtFromLF[2.81,1e-4,1e+4], n2_DimuPtFromLF[2.43,1e-4,1e+4])");
    Pt_Factory_Info[6].Form("PtMassExpPdf::pdfDimuPtFromData(pt, B_DimuPtFromData[2.85,1e-4,1e+4], n1_DimuPtFromData[2.81,1e-4,1e+4], n2_DimuPtFromData[2.43,1e-4,1e+4])");

    RooAbsPdf *pdfDimuPt[n_DiMuSelection];
    RooAbsPdf *pdfDimuM[n_DiMuSelection];
    TF1 *pt_Func[n_DiMuSelection];
    TF1 *m_Func[n_DiMuSelection];
    TCanvas *pt_canvas[n_DiMuSelection];
    TCanvas *m_canvas[n_DiMuSelection];

    for (Int_t i_DimuSel = 0; i_DimuSel < n_DiMuSelection - 1; i_DimuSel++)
    {
        // Select the region of the MC distribution to be extracted and creation of the RooDataSet
        // fTree_cut[i_DimuSel] = (TTree *)fTree[i_DimuSel]->CopyTree(Form("(m>%d && m<%d ) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
        if (info.Low_Mass == 4 && info.High_Mass == 30)
        {
            fTree_x_pt[i_DimuSel] = (TTree *)fTree[i_DimuSel]->CopyTree(Form("((m>%d && m<%0.1f) || (m>%0.1f && m<%d) ) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt));
            fTree_x_pt[i_DimuSel]->SetName(Form("Tree_removeY_%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
            fTree_x_pt[i_DimuSel]->GetBranch("pt")->SetName(pt->GetName());

            fTree_x_m[i_DimuSel] = (TTree *)fTree[i_DimuSel]->CopyTree(Form("(m>%d && m<%d ) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
            fTree_x_m[i_DimuSel]->SetName(Form("Tree_%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
            fTree_x_m[i_DimuSel]->GetBranch("m")->SetName(m->GetName());
        }
        else
        {
            fTree_x_pt[i_DimuSel] = (TTree *)fTree[i_DimuSel]->CopyTree(Form("((m>%d && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
            fTree_x_pt[i_DimuSel]->SetName(Form("Pt_Tree_%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
            fTree_x_pt[i_DimuSel]->GetBranch("pt")->SetName(pt->GetName());

            fTree_x_m[i_DimuSel] = (TTree *)fTree[i_DimuSel]->CopyTree(Form("(m>%d && m<%d ) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
            fTree_x_m[i_DimuSel]->SetName(Form("M_Tree_%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
            fTree_x_m[i_DimuSel]->GetBranch("m")->SetName(m->GetName());
        }

        Pt_Dimu[i_DimuSel] = new RooDataSet(Form("Pt_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), Form("Pt_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), RooArgSet(*pt), Import(*fTree_x_pt[i_DimuSel]));
        M_Dimu[i_DimuSel] = new RooDataSet(Form("M_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), Form("M_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), RooArgSet(*m), Import(*fTree_x_m[i_DimuSel]));

        w->factory(Mass_Factory_Info[i_DimuSel]);
        w->factory(Pt_Factory_Info[i_DimuSel]);

        pdfDimuPt[i_DimuSel] = w->pdf(Form("pdfDimuPtFrom%s", Name_DimuSel[i_DimuSel].Data()));
        auto result1 = pdfDimuPt[i_DimuSel]->fitTo(*Pt_Dimu[i_DimuSel], Minimizer("Minuit2"), Save(), SumW2Error(true));
        result1->printMultiline(out, 0, true, " ");

        pdfDimuM[i_DimuSel] = w->pdf(Form("pdfDimuMassFrom%s", Name_DimuSel[i_DimuSel].Data()));
        RooFitResult *result2;
        if (info.Low_Mass == 4 && info.High_Mass == 30)
            result2 = pdfDimuM[i_DimuSel]->fitTo(*M_Dimu[i_DimuSel], Range("low,high"), SumCoefRange("low,high"), Minimizer("Minuit2"), Save(), SumW2Error(true));
        else
            result2 = pdfDimuM[i_DimuSel]->fitTo(*M_Dimu[i_DimuSel], Minimizer("Minuit2"), Save(), SumW2Error(true));
        result2->printMultiline(out, 0, true, " ");
        result2->Print("v");

        // Drawing the data and pdf on RooPlot

        frameDimuMass[i_DimuSel] = m->frame(Title(Form("frameDimuMass_%s", Name_DimuSel[i_DimuSel].Data())));
        frameDimuMass[i_DimuSel]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");
        M_Dimu[i_DimuSel]->plotOn(frameDimuMass[i_DimuSel], Name(Form("M_%s", Name_DimuSel[i_DimuSel].Data())), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i_DimuSel]), LineColor(color[i_DimuSel]), LineWidth(2), DrawOption("PEZ"), Binning(info.Mass_Binning));
        if (info.Low_Mass == 4 && info.High_Mass == 30)
            pdfDimuM[i_DimuSel]->plotOn(frameDimuMass[i_DimuSel], Name(Form("pdfDimuMFrom%s", Name_DimuSel[i_DimuSel].Data())), Range("low,high"), NormRange("low,high"), LineStyle(kSolid), LineColor(color[i_DimuSel]));
        else
            pdfDimuM[i_DimuSel]->plotOn(frameDimuMass[i_DimuSel], Name(Form("pdfDimuMFrom%s", Name_DimuSel[i_DimuSel].Data())), LineStyle(kSolid), LineColor(color[i_DimuSel]));
        frameDimuMass[i_DimuSel]->SetMaximum(fTree_x_m[i_DimuSel]->GetEntries());
        frameDimuMass[i_DimuSel]->SetMinimum(fTree_x_m[i_DimuSel]->GetEntries() * 1e-5);

        frameDimuPt[i_DimuSel] = pt->frame(Title(Form("frameDimuPt_%s", Name_DimuSel[i_DimuSel].Data())));
        frameDimuPt[i_DimuSel]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
        Pt_Dimu[i_DimuSel]->plotOn(frameDimuPt[i_DimuSel], Name(Form("Pt_%s", Name_DimuSel[i_DimuSel].Data())), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i_DimuSel]), LineColor(color[i_DimuSel]), LineWidth(2), DrawOption("PEZ"), Binning(info.Pt_Binning));
        pdfDimuPt[i_DimuSel]->plotOn(frameDimuPt[i_DimuSel], Name(Form("pdfDimuPtFrom%s", Name_DimuSel[i_DimuSel].Data())), LineStyle(kSolid), LineColor(color[i_DimuSel]));
        frameDimuPt[i_DimuSel]->SetMaximum(fTree_x_pt[i_DimuSel]->GetEntries());
        frameDimuPt[i_DimuSel]->SetMinimum(fTree_x_pt[i_DimuSel]->GetEntries() * 1e-5);

        // Conversion of the RooDataSet in TH1 & RooPdf in a TF1 object to calculate the ratio

        m_histo[i_DimuSel] = M_Dimu[i_DimuSel]->createHistogram(Form("h_M_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), *m, Binning(info.Mass_Binning, info.Low_Mass, info.High_Mass));
        m_histo[i_DimuSel]->SetMarkerStyle(24);
        m_histo[i_DimuSel]->SetMarkerSize(1.75);
        m_histo[i_DimuSel]->SetMarkerColor(color[i_DimuSel]);
        m_histo[i_DimuSel]->SetMarkerColor(color[i_DimuSel]);

        pt_histo[i_DimuSel] = Pt_Dimu[i_DimuSel]->createHistogram(Form("h_Pt_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), *pt, Binning(info.Pt_Binning, info.Low_Pt, info.High_Pt));
        pt_histo[i_DimuSel]->SetMarkerStyle(24);
        pt_histo[i_DimuSel]->SetMarkerSize(1.75);
        pt_histo[i_DimuSel]->SetMarkerColor(color[i_DimuSel]);
        pt_histo[i_DimuSel]->SetMarkerColor(color[i_DimuSel]);

        RooArgSet *pt_model = static_cast<RooArgSet *>(RooArgSet(*pdfDimuPt[i_DimuSel]).snapshot(true));                       // True means copy the PDF and everything it depends on
        auto &pt_modelcopied = static_cast<RooAbsPdf &>((*pt_model)[Form("pdfDimuPtFrom%s", Name_DimuSel[i_DimuSel].Data())]); // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
        RooArgSet *pt_modelobs = pt_modelcopied.getObservables(*Pt_Dimu[i_DimuSel]);
        RooArgSet *pt_modelPars = pt_modelcopied.getParameters(*pt_modelobs);
        pt_modelPars->setName(Name_DimuSel[i_DimuSel].Data());
        pt_Func[i_DimuSel] = pt_modelcopied.asTF(*pt_modelobs, *pt_modelPars, *pt);
        pt_Func[i_DimuSel]->SetName(Form("pdfDimuPtFrom%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
        pt_Func[i_DimuSel]->SetLineWidth(2);
        pt_Func[i_DimuSel]->SetLineColor(color[i_DimuSel]);
        pt_Func[i_DimuSel]->Write(0, 2, 0);

        RooArgSet *m_model = static_cast<RooArgSet *>(RooArgSet(*pdfDimuM[i_DimuSel]).snapshot(true));                         // True means copy the PDF and everything it depends on
        auto &m_modelcopied = static_cast<RooAbsPdf &>((*m_model)[Form("pdfDimuMassFrom%s", Name_DimuSel[i_DimuSel].Data())]); // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
        RooArgSet *m_modelobs = m_modelcopied.getObservables(*M_Dimu[i_DimuSel]);
        RooArgSet *m_modelPars = m_modelcopied.getParameters(*m_modelobs);
        m_modelPars->setName(Name_DimuSel[i_DimuSel].Data());
        m_Func[i_DimuSel] = m_modelcopied.asTF(*m_modelobs, *m_modelPars, *m);
        m_Func[i_DimuSel]->SetName(Form("pdfDimuMassFrom%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
        m_Func[i_DimuSel]->SetLineWidth(2);
        m_Func[i_DimuSel]->SetLineColor(color[i_DimuSel]);
        m_Func[i_DimuSel]->Write(0, 2, 0);
        RooArgSet *param = nullptr;
        m_canvas[i_DimuSel] = printMC_ratio(Form("m_canvas_%s_%s_PDF_%s_Mcut_%0.1f_%0.1f", Name_DimuSel[i_DimuSel].Data(), info.Generator.Data(), info.LF_HF.Data(), info.LowM_cut, info.HighM_cut), TString::Format("#it{p}_{T} < %d GeV/#it{c}", info.High_Pt), frameDimuMass[i_DimuSel], m_modelPars, m_histo[i_DimuSel], m_Func[i_DimuSel], color[i_DimuSel], info.Low_Mass, info.High_Mass);
        m_canvas[i_DimuSel]->SaveAs(Form("images/pdf_extraction/%s.png", m_canvas[i_DimuSel]->GetName()));
        m_canvas[i_DimuSel]->SaveAs(Form("images/pdf_extraction/%s.pdf", m_canvas[i_DimuSel]->GetName()));

        TString Title;
        if (info.Low_Mass == 4 && info.High_Mass == 30)
            Title.Form("%d < #it{m}_{#mu#mu} < %0.1f && %0.1f < #it{m}_{#mu#mu} < %d", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass);
        else
            Title.Form("%d < #it{m}_{#mu#mu} < %d", info.Low_Mass, info.High_Mass);
        pt_canvas[i_DimuSel] = printMC_ratio(Form("pt_canvas_%s_%s_PDF_%s_Mcut_%0.1f_%0.1f", Name_DimuSel[i_DimuSel].Data(), info.Generator.Data(), info.LF_HF.Data(), info.LowM_cut, info.HighM_cut), Title, frameDimuPt[i_DimuSel], pt_modelPars, pt_histo[i_DimuSel], pt_Func[i_DimuSel], color[i_DimuSel], info.Low_Pt, info.High_Pt);
        pt_canvas[i_DimuSel]->SaveAs(Form("images/pdf_extraction/%s.png", pt_canvas[i_DimuSel]->GetName()));
        pt_canvas[i_DimuSel]->SaveAs(Form("images/pdf_extraction/%s.pdf", pt_canvas[i_DimuSel]->GetName()));

        // Saving the paramater of the fit for Systematic errors calculations

        B_DimuMass[i_DimuSel] = w->var(Form("B_DimuMassFrom%s", Name_DimuSel[i_DimuSel].Data()));
        B_DimuMass[i_DimuSel]->setConstant(kTRUE);
        n1_DimuMass[i_DimuSel] = w->var(Form("n1_DimuMassFrom%s", Name_DimuSel[i_DimuSel].Data()));
        n1_DimuMass[i_DimuSel]->setConstant(kTRUE);
        n2_DimuMass[i_DimuSel] = w->var(Form("n2_DimuMassFrom%s", Name_DimuSel[i_DimuSel].Data()));
        n2_DimuMass[i_DimuSel]->setConstant(kTRUE);

        B_DimuPt[i_DimuSel] = w->var(Form("B_DimuPtFrom%s", Name_DimuSel[i_DimuSel].Data()));
        B_DimuPt[i_DimuSel]->setConstant(kTRUE);
        n1_DimuPt[i_DimuSel] = w->var(Form("n1_DimuPtFrom%s", Name_DimuSel[i_DimuSel].Data()));
        n1_DimuPt[i_DimuSel]->setConstant(kTRUE);
        n2_DimuPt[i_DimuSel] = w->var(Form("n2_DimuPtFrom%s", Name_DimuSel[i_DimuSel].Data()));
        n2_DimuPt[i_DimuSel]->setConstant(kTRUE);

        Param_Pt_unbinned[i_DimuSel] = new TH1D(Form("Param_Pt_unbinned_from%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "; coeff x Mass fit", 4, 0, 4);
        Param_Pt_unbinned[i_DimuSel]->SetBinContent(1, B_DimuPt[i_DimuSel]->getVal());
        Param_Pt_unbinned[i_DimuSel]->SetBinContent(2, n1_DimuPt[i_DimuSel]->getVal());
        Param_Pt_unbinned[i_DimuSel]->SetBinContent(3, n2_DimuPt[i_DimuSel]->getVal());
        Param_Pt_unbinned[i_DimuSel]->SetBinContent(4, fTree_x_pt[i_DimuSel]->GetEntries());

        Param_Pt_unbinned[i_DimuSel]->SetBinError(1, B_DimuPt[i_DimuSel]->getError());
        Param_Pt_unbinned[i_DimuSel]->SetBinError(2, n1_DimuPt[i_DimuSel]->getError());
        Param_Pt_unbinned[i_DimuSel]->SetBinError(3, n2_DimuPt[i_DimuSel]->getError());
        Param_Pt_unbinned[i_DimuSel]->SetBinError(4, 1);
        Param_Pt_unbinned[i_DimuSel]->Write(0, 2, 0);

        Param_M_unbinned[i_DimuSel] = new TH1D(Form("Param_M_unbinned_from%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "; coeff x Pt fit", 4, 0, 4);
        Param_M_unbinned[i_DimuSel]->SetBinContent(1, B_DimuMass[i_DimuSel]->getVal());
        Param_M_unbinned[i_DimuSel]->SetBinContent(2, n1_DimuMass[i_DimuSel]->getVal());
        Param_M_unbinned[i_DimuSel]->SetBinContent(3, n2_DimuMass[i_DimuSel]->getVal());
        Param_M_unbinned[i_DimuSel]->SetBinContent(4, fTree_x_pt[i_DimuSel]->GetEntries());

        Param_M_unbinned[i_DimuSel]->SetBinError(1, B_DimuMass[i_DimuSel]->getError());
        Param_M_unbinned[i_DimuSel]->SetBinError(2, n1_DimuMass[i_DimuSel]->getError());
        Param_M_unbinned[i_DimuSel]->SetBinError(3, n2_DimuMass[i_DimuSel]->getError());
        Param_M_unbinned[i_DimuSel]->SetBinError(4, 1);
        Param_M_unbinned[i_DimuSel]->Write(0, 2, 0);
    }

    w->Write();

    fOut->ls();

    for (Int_t i_DiMu = 0; i_DiMu < n_DiMuSelection - 1; i_DiMu++)
    {
        output_fit = info.Name_DimuSel[i_DiMu];
        out << Form("%s \n", output_fit.Data());
        output_fit = FileName[i_DiMu];
        out << Form("%s \n", output_fit.Data());
        output_fit = Mass_Factory_Info[i_DiMu];
        out << Form("%s \n", output_fit.Data());
        output_fit.Form("B_DimuMassFrom%s: %0.3e +/- %0.3e|| n1_DimuMassFrom%s: %0.3e +/- %0.3e || n2_DimuMassFrom%s: %0.3e +/- %0.3e\n", info.Name_DimuSel[i_DiMu].Data(), B_DimuMass[i_DiMu]->getVal(), B_DimuMass[i_DiMu]->getError(), info.Name_DimuSel[i_DiMu].Data(), n1_DimuMass[i_DiMu]->getVal(), n1_DimuMass[i_DiMu]->getError(), info.Name_DimuSel[i_DiMu].Data(), n2_DimuMass[i_DiMu]->getVal(), n2_DimuMass[i_DiMu]->getError());
        out << Form("%s \n", output_fit.Data());
        output_fit = Pt_Factory_Info[i_DiMu];
        out << Form("%s \n", output_fit.Data());
        output_fit.Form("B_DimuPtFrom%s: %0.3e +/- %0.3e|| n1_DimuPtFrom%s: %0.3e +/- %0.3e || n2_DimuPtFrom%s: %0.3e +/- %0.3e\n", info.Name_DimuSel[i_DiMu].Data(), B_DimuPt[i_DiMu]->getVal(), B_DimuPt[i_DiMu]->getError(), info.Name_DimuSel[i_DiMu].Data(), n1_DimuPt[i_DiMu]->getVal(), n1_DimuPt[i_DiMu]->getError(), info.Name_DimuSel[i_DiMu].Data(), n2_DimuPt[i_DiMu]->getVal(), n2_DimuPt[i_DiMu]->getError());
        out << Form("%s \n", output_fit.Data());
    }
    out.close();
}

void MC_Data_shape()
{
    TFile *fY_cut = new TFile("results/removing_Y_POWHEG_PDF_full_stat.root", "READ");
    const Int_t n_DiMuSelection = 7;
    TString Name_DimuSel[n_DiMuSelection] = {"Charm", "Beauty", "HF_Mixed", "DY", "LF_HF_Mixed", "LF", "Data"};
    TH1F *Data_shapes = new TH1F("Data_shapes", ";#it{p}_{T}", 300, 0, 30);
    Data_shapes = (TH1F *)((TF1 *)fY_cut->Get(Form("pdfDimuPtFrom%s_M_4_30_Pt_0_30", Name_DimuSel[6].Data())))->CreateHistogram();
    Data_shapes->SetTitle("Data fit");
    TH1F *MC_shapes_sum = new TH1F("MC_shapes_sum", "Sum MC contr. fit", 300, 0, 30);
    for (Int_t i_DimuSel = 0; i_DimuSel < n_DiMuSelection; i_DimuSel++)
    {

        MC_shapes_sum->Add((TH1F *)((TF1 *)fY_cut->Get(Form("pdfDimuPtFrom%s_M_4_30_Pt_0_30", Name_DimuSel[i_DimuSel].Data())))->CreateHistogram());
    }
    TCanvas *c = new TCanvas("Fit_MC_Data_shape", "Fit_MC_Data_shape", 900, 1000);
    c->Divide(1, 2);
    c->cd(1);
    gPad->SetLogy();
    MC_shapes_sum->Scale(1. / MC_shapes_sum->Integral(), "width");
    Data_shapes->SetLineColor(kBlue);
    MC_shapes_sum->SetLineColor(kRed);
    Data_shapes->DrawCopy();
    MC_shapes_sum->DrawCopy("hist same");
    gPad->BuildLegend();
    c->cd(2);
    TH1F *ratio = (TH1F *)Data_shapes->Clone("ratio");
    ratio->Divide(MC_shapes_sum);
    ratio->DrawCopy("hist");
}

void cut_bias()
{
    gROOT->ProcessLineSync(".x /home/michele_pennisi/dimuon_HF_pp/fit_data/fit_library/PtMassExpPdf.cxx+");
    opt info;
    const Int_t n_DiMuSelection = 7;
    TString Name_DimuSel[n_DiMuSelection] = {"Charm", "Beauty", "HF_Mixed", "DY", "LF_HF_Mixed", "LF", "Data"};
    Color_t color[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kAzure + 9, kOrange + 7, kBlue, kRed, kBlack};

    TString Dir_name = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output";
    TString Version_ALIAOD;
    TString Tree_name;
    if (info.stat_MC.Contains("small_stat"))
    {
        Version_ALIAOD = "Version3_AliAOD/save_output";
        Tree_name = "DiMuon_Rec_FullMass_PowhegOnly";
    }
    else if (info.stat_MC.Contains("full_stat"))
    {
        Version_ALIAOD = "Version_5_AliAOD_skimmed_fwd_fullstat";
        Tree_name = "DiMuon_Rec_PowhegOnly_";
    }
    cout << Version_ALIAOD.Data() << endl;

    TFile *fIn[n_DiMuSelection] = {new TFile(Form("%s/LHC23i1/%s/LHC23i1_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())), new TFile(Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())), new TFile(Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())), new TFile(Form("%s/LHC18p_DY/Version_2_AliAOD/LHC18p_DY_MC_output_Tree_merged.root", Dir_name.Data())), new TFile(Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())), new TFile(Form("%s/LHC23i2/%s/LHC23i2_MC_output_Tree_merged.root", Dir_name.Data(), Version_ALIAOD.Data())), new TFile("/home/michele_pennisi/cernbox/HF_dimuons/data_analysis/Tree_MassPt_MassCut4_Run2.root", "READ")};
    fIn[0]->ls();
    TTree *fTree[n_DiMuSelection] = {(TTree *)fIn[0]->Get(Form("%sCharm", Tree_name.Data())), (TTree *)fIn[1]->Get(Form("%sBeauty", Tree_name.Data())), (TTree *)fIn[2]->Get(Form("%sHF_Mixed", Tree_name.Data())), (TTree *)fIn[3]->Get("DiMuon_Rec_DY"), (TTree *)fIn[4]->Get("DiMuon_Rec_PythiaOnly_LF_HF_Mixed"), (TTree *)fIn[5]->Get("DiMuon_Rec_PythiaOnly_LF"), (TTree *)fIn[6]->Get("rec_data_tree")};
    gROOT->cd();

    Int_t i_DimuSel = 6;

    TTree *fTree_std_cut = (TTree *)fTree[i_DimuSel]->CopyTree(Form("(m>%d && m<%d) && (pt > %d && pt <%d)", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt)); // standard cut for MC, no esclusion of Y region

    fTree_std_cut->SetName(Form("fTree_std_cut_%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    RooWorkspace *w = new RooWorkspace(Form("w_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "workspace");

    w->factory("pt[0,30], #it{p}_{T} (GeV/#it{c})");
    RooRealVar *pt = (RooRealVar *)w->var("pt");
    pt->setBins(info.Pt_Binning);
    w->factory("m[4,30], #it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    RooRealVar *m = w->var("m");
    m->setBins(info.Mass_Binning);

    TString Pt_Factory_Info;
    TString M_Factory_Info;
    // Construction fit without removing Y region
    Pt_Factory_Info.Form("PtMassExpPdf::pdfDimuPtFrom%s(pt, B_DimuPtFrom%s[2.85,1e-4,1e+4], n1_DimuPtFrom%s[2.81,1e-4,1e+4], n2_DimuPtFrom%s[2.43,1e-4,1e+4])", Name_DimuSel[i_DimuSel].Data(), Name_DimuSel[i_DimuSel].Data(), Name_DimuSel[i_DimuSel].Data(), Name_DimuSel[i_DimuSel].Data());
    M_Factory_Info.Form("PtMassExpPdf::pdfDimuMFrom%s(m, B_DimuMassFrom%s[2.85,1e-4,1e+4], n1_DimuMassFrom%s[2.81,1e-4,1e+4], n2_DimuMassFrom%s[5,1e-3,100])", Name_DimuSel[i_DimuSel].Data(), Name_DimuSel[i_DimuSel].Data(), Name_DimuSel[i_DimuSel].Data(), Name_DimuSel[i_DimuSel].Data());

    w->factory(Pt_Factory_Info);
    w->factory(M_Factory_Info);

    RooDataSet *Pt_Dimu_std_cut = new RooDataSet(Form("Pt_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), Form("Pt_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), RooArgSet(*pt), Import(*fTree_std_cut));
    RooDataSet *M_Dimu_std_cut = new RooDataSet(Form("M_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), Form("M_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), RooArgSet(*m), Import(*fTree_std_cut));

    RooAbsPdf *pdfDimuPt = w->pdf(Form("pdfDimuPtFrom%s", Name_DimuSel[i_DimuSel].Data()));
    RooAbsPdf *pdfDimuM = w->pdf(Form("pdfDimuMFrom%s", Name_DimuSel[i_DimuSel].Data()));

    auto Pt_result1 = pdfDimuPt->fitTo(*Pt_Dimu_std_cut, Minimizer("Minuit2"), Save(), SumW2Error(true));
    auto M_result1 = pdfDimuM->fitTo(*M_Dimu_std_cut, Minimizer("Minuit2"), Save(), SumW2Error(true));

    // Rooplot def
    RooPlot *frameDimuPt = pt->frame(Name(Form("std_cut_frameDimuPt_%s", Name_DimuSel[i_DimuSel].Data())), Title(Form("frameDimuPt_%s", Name_DimuSel[i_DimuSel].Data())));
    frameDimuPt->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    Pt_Dimu_std_cut->plotOn(frameDimuPt, Name(Form("Pt_%s", Name_DimuSel[i_DimuSel].Data())), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i_DimuSel]), LineColor(color[i_DimuSel]), LineWidth(2), DrawOption("PEZ"), Binning(info.Pt_Binning));
    pdfDimuPt->plotOn(frameDimuPt, Name(TString::Format("pdfDimuPtFrom%s", Name_DimuSel[i_DimuSel].Data()).Data()), LineStyle(kSolid), LineColor(color[i_DimuSel]));

    RooPlot *frameDimuM = m->frame(Name(Form("std_cut_frameDimuM_%s", Name_DimuSel[i_DimuSel].Data())), Title(Form("frameDimuM_%s", Name_DimuSel[i_DimuSel].Data())));
    frameDimuM->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
    M_Dimu_std_cut->plotOn(frameDimuM, Name(Form("M_%s", Name_DimuSel[i_DimuSel].Data())), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i_DimuSel]), LineColor(color[i_DimuSel]), LineWidth(2), DrawOption("PEZ"), Binning(info.Mass_Binning));
    pdfDimuM->plotOn(frameDimuM, Name(TString::Format("pdfDimuMFrom%s", Name_DimuSel[i_DimuSel].Data()).Data()), LineStyle(kSolid), LineColor(color[i_DimuSel]));

    // TF1 and TH1F Pt extraction for drawing
    RooArgSet *Pt_model = static_cast<RooArgSet *>(RooArgSet(*pdfDimuPt).snapshot(true));                                  // True means copy the PDF and everything it depends on
    auto &Pt_modelcopied = static_cast<RooAbsPdf &>((*Pt_model)[Form("pdfDimuPtFrom%s", Name_DimuSel[i_DimuSel].Data())]); // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    RooArgSet *Pt_modelobs = Pt_modelcopied.getObservables(*Pt_Dimu_std_cut);
    RooArgSet *Pt_modelPars = Pt_modelcopied.getParameters(*Pt_modelobs);
    Pt_modelPars->setName(Name_DimuSel[i_DimuSel].Data());

    TF1 *Pt_Func = Pt_modelcopied.asTF(*Pt_modelobs, *Pt_modelPars, *pt);
    Pt_Func->SetName(Form("std_cut_pdfDimuPtFrom%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    Pt_Func->SetTitle(Form("Fit no Y cut"));

    TH1F *Pt_histo = (TH1F *)Pt_Dimu_std_cut->createHistogram(Form("std_cut_h_Pt_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), *pt, Binning(info.Pt_Binning, info.Low_Pt, info.High_Pt));

    TH1F *ratio_std_Y_cut = (TH1F *)Pt_Func->CreateHistogram();
    ratio_std_Y_cut->SetName("ratio_std_Y_cut");
    ratio_std_Y_cut->SetTitle("fit ratios");
    ratio_std_Y_cut->SetLineColor(kRed);

    // TF1 and TH1F Mass extraction for drawing

    RooArgSet *M_model = static_cast<RooArgSet *>(RooArgSet(*pdfDimuM).snapshot(true));                                 // True means copy the PDF and everything it depends on
    auto &M_modelcopied = static_cast<RooAbsPdf &>((*M_model)[Form("pdfDimuMFrom%s", Name_DimuSel[i_DimuSel].Data())]); // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    RooArgSet *M_modelobs = M_modelcopied.getObservables(*M_Dimu_std_cut);
    RooArgSet *M_modelPars = M_modelcopied.getParameters(*M_modelobs);
    M_modelPars->setName(Name_DimuSel[i_DimuSel].Data());

    TF1 *M_Func = M_modelcopied.asTF(*M_modelobs, *M_modelPars, *m);
    M_Func->SetName(Form("std_cut_pdfDimuMFrom%s_M_%d_%d_M_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Mass, info.High_Mass));
    M_Func->SetTitle(Form("std_cut_pdfDimuMFrom%s_M_%d_%d_M_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Mass, info.High_Mass));

    TH1F *M_histo = (TH1F *)M_Dimu_std_cut->createHistogram(Form("std_cut_h_M_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), *m, Binning(info.Mass_Binning, info.Low_Mass, info.High_Mass));

    // Drawing Canvas

    TString c_Title = "4 < #it{m}_{#mu#mu} < 30 GeV/#it{c}^2";

    frameDimuPt->SetMaximum(1.2e+6);
    frameDimuPt->SetMinimum(1.2e-3);

    TCanvas *Pt_canvas = printMC_ratio(Form("std_cut_Pt_canvas_%s", Name_DimuSel[i_DimuSel].Data()), c_Title, frameDimuPt, Pt_modelPars, Pt_histo, Pt_Func, color[i_DimuSel], info.Low_Pt, info.High_Pt);
    Pt_canvas->SaveAs(Form("images/cut_bias/%s.png", Pt_canvas->GetName()));

    frameDimuM->SetMaximum(1.2e+6);
    frameDimuM->SetMinimum(1.2e-3);

    TCanvas *M_canvas = printMC_ratio(Form("std_cut_M_canvas_%s", Name_DimuSel[i_DimuSel].Data()), c_Title, frameDimuM, M_modelPars, M_histo, M_Func, color[i_DimuSel], info.Low_Mass, info.High_Mass);
    M_canvas->SaveAs(Form("images/cut_bias/%s.png", M_canvas->GetName()));

    TTree *fTree_Y_cut = (TTree *)fTree[i_DimuSel]->CopyTree(Form("((m>%d && m<%0.1f) || (m>%0.0f && m<%d)) && (pt > %d && pt <%d)", info.Low_Mass, info.LowM_cut, info.HighM_cut, info.High_Mass, info.Low_Pt, info.High_Pt)); // cut escluding Y region
    fTree_Y_cut->SetName(Form("fTree_Y_cut_%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));

    TCanvas *Pt_control = new TCanvas(Form("Pt_control_%s", Name_DimuSel[i_DimuSel].Data()), Form("Pt_control_%s", Name_DimuSel[i_DimuSel].Data()), 900, 1000);
    Pt_control->Divide(1, 2);
    Pt_control->cd(1);
    gPad->SetLogy();
    TH1F *htemp = new TH1F("htemp", ";#it{p}_{T} (GeV/#it{c})", info.Pt_Binning, 0, 30);
    fTree_std_cut->Draw("pt>>htemp", "", "goff");
    htemp->Rebin(2);
    htemp->Scale(1. / htemp->GetEntries(), "width");
    htemp->SetTitle(Form("4 < #it{m}< 9 (GeV/#it{c}^{2}), #it{N}_{#mu#mu}: %0.0f", htemp->GetEntries()));
    htemp->DrawCopy();

    TH1F *ratio = (TH1F *)htemp->Clone("ratio");
    ratio->SetTitle("std cut / Y cut");
    fTree_Y_cut->Draw("pt>>htemp", "", "goff");
    htemp->Scale(1. / htemp->GetEntries(), "width");
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerSize(1.5);
    ratio->Divide(htemp);
    ratio->GetYaxis()->SetRangeUser(0.2, 2.3);
    htemp->SetTitle(Form("Y cut, %0.1f < #it{m} < %0.1f  (GeV/#it{c}^{2}), #it{N}_{#mu#mu}: %0.0f", info.LowM_cut, info.HighM_cut, htemp->GetEntries()));
    htemp->SetLineColor(kRed);
    htemp->DrawCopy("same");
    Pt_Func->SetLineWidth(2);
    Pt_Func->SetLineColor(color[i_DimuSel]);
    Pt_Func->SetTitle("fit");
    Pt_Func->DrawCopy("same");
    TFile *fY_cut = new TFile("results/removing_Y_POWHEG_PDF_full_stat.root", "READ");
    TF1 *Pt_Func_Y_removed = (TF1 *)fY_cut->Get(Form("pdfDimuPtFrom%s_M_4_30_Pt_0_30", Name_DimuSel[i_DimuSel].Data()));
    Pt_Func_Y_removed->SetLineStyle(kDashed);
    Pt_Func_Y_removed->SetTitle("fit Y removed");
    Pt_Func_Y_removed->SetLineWidth(2);
    Pt_Func_Y_removed->DrawCopy("same");
    gPad->BuildLegend();
    Pt_control->cd(2);
    ratio->SetTitle("distr ratios");
    ratio->DrawCopy("PE");

    TH1F *ratio_fit = (TH1F *)Pt_Func->CreateHistogram();
    ratio_fit->Divide((TH1F *)Pt_Func_Y_removed->CreateHistogram());
    ratio_fit->SetTitle("fit ratios");
    ratio_fit->DrawCopy("same");
    gPad->BuildLegend();
    Pt_control->SaveAs(Form("images/cut_bias/%s.png", Pt_control->GetName()));

    TCanvas *M_control = new TCanvas(Form("M_control_%s", Name_DimuSel[i_DimuSel].Data()), Form("M_control_%s", Name_DimuSel[i_DimuSel].Data()), 900, 1000);
    M_control->Divide(1, 2);
    M_control->cd(1);
    gPad->SetLogy();
    TH1F *hM_temp = new TH1F("hM_temp", ";#it{m}_{#mu#mu} (GeV/#it{c}^2)", info.Mass_Binning, 4, 30);
    fTree_std_cut->Draw("m>>hM_temp", "", "goff");
    hM_temp->Rebin(2);
    cout << "Entries: " << hM_temp->GetEntries() << endl;
    hM_temp->Scale(1. / hM_temp->GetEntries(), "width");
    hM_temp->SetTitle(Form("4 < #it{m}< 9 (GeV/#it{c}^{2}), #it{N}_{#mu#mu}: %0.0f", hM_temp->GetEntries()));
    hM_temp->DrawCopy();

    ratio = (TH1F *)hM_temp->Clone("ratio");
    ratio->SetTitle("std cut / Y cut");
    Double_t Int = hM_temp->GetEntries();
    fTree_Y_cut->Draw("m>>hM_temp", "", "goff");
    hM_temp->Scale(1. / Int, "width");
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerSize(1.5);
    ratio->Divide(hM_temp);
    ratio->GetYaxis()->SetRangeUser(0.8, ratio->GetMaximum() * 1.1);
    hM_temp->SetTitle(Form("Y cut, %0.1f < #it{m} < %0.1f  (GeV/#it{c}^{2}), #it{N}_{#mu#mu}: %0.0f", info.LowM_cut, info.HighM_cut, hM_temp->GetEntries()));
    hM_temp->SetLineColor(kRed);
    hM_temp->DrawCopy("same");
    M_Func->SetLineWidth(2);
    M_Func->SetLineColor(color[i_DimuSel]);
    M_Func->SetTitle("fit");
    M_Func->DrawCopy("same");
    TF1 *M_Func_Y_removed = (TF1 *)fY_cut->Get(Form("pdfDimuMassFrom%s_M_4_30_Pt_0_30", Name_DimuSel[i_DimuSel].Data()));
    M_Func_Y_removed->SetLineStyle(kDashed);
    M_Func_Y_removed->SetTitle("fit Y removed");
    M_Func_Y_removed->SetLineWidth(2);
    M_Func_Y_removed->DrawCopy("same");
    gPad->BuildLegend();
    M_control->cd(2);
    ratio->SetTitle("distr ratios");
    ratio->DrawCopy("PE");

    ratio_fit = (TH1F *)M_Func->CreateHistogram();
    ratio_fit->Divide((TH1F *)M_Func_Y_removed->CreateHistogram());
    ratio_fit->SetTitle("fit ratios");
    ratio_fit->DrawCopy("same");
    gPad->BuildLegend();
    M_control->SaveAs(Form("images/cut_bias/%s.png", M_control->GetName()));

    // TCanvas *M_control = new TCanvas(Form("M_control_%s", Name_DimuSel[i_DimuSel].Data()), Form("M_control_%s", Name_DimuSel[i_DimuSel].Data()), 900, 1000);
    // M_control->Divide(1, 2);
    // M_control->cd(1);
    // gPad->SetLogy();
    // TH1F *htemp = new TH1F("htemp", ";#it{p}_{T} (GeV/#it{c})", info.Pt_Binning, 0, 30);
    // fTree_std_cut->Draw("pt>>htemp", "", "goff");
    // htemp->Rebin(2);
    // htemp->Scale(1. / htemp->GetEntries(), "width");
    // htemp->SetTitle(Form("4 < #it{m}< 9 (GeV/#it{c}^{2}), #it{N}_{#mu#mu}: %0.0f", htemp->GetEntries()));
    // htemp->DrawCopy();

    // M_control->cd();
    // M_Func->SetLineColor(kRed);
    // M_Func->DrawCopy("same");

    //     Pt_control->cd(1);
    // Pt_Func->SetLineColor(kBlue);
    // Pt_Func->SetLineWidth(2);
    // Pt_Func->DrawCopy("same");

    // Pt_control->cd(1);
    // Pt_Func->SetLineColor(kRed);
    // Pt_Func->SetLineWidth(2);
    // Pt_Func->DrawCopy("same");

    // Pt_control->cd(2);
    // ratio_std_Y_cut->Divide((TH1F *)Pt_Func->CreateHistogram());
    // ratio_std_Y_cut->Draw("same");
    // gPad->BuildLegend();
    // gPad->Modified();
    // gPad->Update();
    // Pt_control->SaveAs(Form("images/cut_bias/%s.png", Pt_control->GetName()));
}

void data_pdf()
{
    opt info;
    gROOT->ProcessLineSync(".x /home/michele_pennisi/dimuon_HF_pp/fit_data/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/dimuon_HF_pp/fit_data/fit_library/PtMassPol1ExpPdf.cxx+");

    TFile *fIn_data;
    if (info.stat_Data.Contains("Run2"))
        fIn_data = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/data_analysis/Tree_MassPt_MassCut4_Run2.root", "READ");
    else
        fIn_data = new TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/TreeResults_merged.root", "READ");
    fIn_data->ls();

    // Taking data saved in tree
    TTree *tree_data = (TTree *)fIn_data->Get("rec_data_tree");
    // tree_data->Draw("m", Form("((m>%d && m<9) || (m>11 && m<%d)) && (pt > %d && pt <%d)", Low_Mass, High_Mass, Low_Pt, High_Pt));
    gROOT->cd();
    RooRealVar m("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", 4., 30.);
    RooRealVar low_m("low_m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", 4., 9.);
    RooRealVar high_m("high_m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", 11., 30.);

    RooCategory sample("sample", "sample");
    sample.defineType("low_mass");
    sample.defineType("high_mass");

    RooRealVar pt("pt", "#it{p}_{T} (GeV/#it{c})", info.Low_Pt, info.High_Pt);
    // pt.setRange("pluto", 4, info.LowM_cut);
    // pt.setRange("pippo", info.HighM_cut, 30);
    // Select the region of the MC distribution to be extracted and creation of the RooDataSet

    TTree *low_m_tree_data_cutted = (TTree *)tree_data->CopyTree(Form("(m>%d && m<%0.1f) && pt<%d", 4, 9., 30));
    low_m_tree_data_cutted->GetBranch("m")->SetName("low_m");
    low_m_tree_data_cutted->Print();

    RooDataSet *low_M_Dimu = new RooDataSet("low_M_Dimu_data", "M_Dimu_data", RooArgSet(low_m), Import(*low_m_tree_data_cutted));
    low_M_Dimu->Print();

    TTree *high_m_tree_data_cutted = (TTree *)tree_data->CopyTree(Form("(m>%d && m<%0.1f) && pt<%d", 11, 30., 30));
    high_m_tree_data_cutted->GetBranch("m")->SetName("high_m");
    high_m_tree_data_cutted->Print();

    RooDataSet *high_M_Dimu = new RooDataSet("high_M_Dimu_data", "M_Dimu_data", RooArgSet(high_m), Import(*high_m_tree_data_cutted));
    high_M_Dimu->Print();

    RooDataSet *unbinned_combData_set = new RooDataSet("combData", "combined data", RooArgSet(low_m, high_m), Index(sample), Import("low_mass", *low_M_Dimu), Import("high_mass", *high_M_Dimu));
    unbinned_combData_set->Print();

    RooWorkspace *w = new RooWorkspace(Form("w_M_%d_%d_Pt_%d_%d", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt), "workspace");

    TString Mass_Factory_Info;
    Mass_Factory_Info.Form("PtMassExpPdf::pdfDimu_Low_MassFromData(%s[%d,%d], B_DimuMassFromData[2.85,1e-4,1e+4], n1_DimuMassFromData[2.81,1e-4,1e+4], n2_DimuMassFromData[5,1e-3,100])", low_m.GetName(), 4, 9);
    w->factory(Mass_Factory_Info);

    RooAbsPdf *low_pdfDimuM = w->pdf("pdfDimu_Low_MassFromData");

    Mass_Factory_Info.Form("PtMassExpPdf::pdfDimu_High_MassFromData(%s[%d,%d], B_DimuMassFromData[2.85,1e-4,1e+4], n1_DimuMassFromData[2.81,1e-4,1e+4], n2_DimuMassFromData[5,1e-3,100])", high_m.GetName(), 11, 30);

    w->factory(Mass_Factory_Info);

    RooAbsPdf *high_pdfDimuM = w->pdf("pdfDimu_High_MassFromData");

    RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
    simPdf.addPdf(*low_pdfDimuM, "low_mass");
    simPdf.addPdf(*high_pdfDimuM, "high_mass");
    simPdf.Print("t");
    // RooFitResult *r = simPdf.fitTo(*unbinned_combData_set, Save(), SumW2Error(true));

    // RooRealVar *B_DimuMass = w->var("B_DimuMassFromData");
    // B_DimuMass->setConstant(kTRUE);
    // RooRealVar *n1_DimuMass = w->var("n1_DimuMassFromData");
    // n1_DimuMass->setConstant(kTRUE);
    // RooRealVar *n2_DimuMass = w->var("n2_DimuMassFromData");
    // n2_DimuMass->setConstant(kTRUE);

    TF1 *pdf_low_m = new TF1("pdf_low_m", FuncPtMass, 4, 9, 4);
    pdf_low_m->SetTitle("Fit");
    // pdf_low_m->FixParameter(0, B_DimuMass->getVal());
    // pdf_low_m->FixParameter(1, n1_DimuMass->getVal());
    // pdf_low_m->FixParameter(2, n2_DimuMass->getVal());

    // pdf_low_m->SetParError(0, B_DimuMass->getError());
    // pdf_low_m->SetParError(1, n1_DimuMass->getError());
    // pdf_low_m->SetParError(2, n2_DimuMass->getError());

    pdf_low_m->FixParameter(0, 1.60856e+01);
    pdf_low_m->FixParameter(1, 2.70776e-01);
    pdf_low_m->FixParameter(2, 4.55889e+01);

    // f_M_PYTHIA6->SetParameter(3, 1);
    pdf_low_m->SetLineColor(kMagenta + 2);
    pdf_low_m->SetLineStyle(kDashed);
    pdf_low_m->SetLineWidth(3);

    TH1F *low_m_histo = new TH1F("low_m_histo", "low_m_histo", 150, 4, 9);
    low_m_tree_data_cutted->Draw("low_m>>low_m_histo", "", "goff");
    low_m_histo->Scale(1., "width");
    low_m_histo->SetMarkerSize(1.5);
    low_m_histo->SetMarkerStyle(20);
    low_m_histo->SetMarkerColor(kBlack);
    low_m_histo->Fit(pdf_low_m, "R0I");

    TH1F *ratio_low_m = (TH1F *)low_m_histo->Clone("ratio_lowm_fit");
    ratio_low_m->Divide(pdf_low_m);
    ratio_low_m->GetYaxis()->SetTitle("MC/Fit");
    ratio_low_m->GetYaxis()->SetRangeUser(-0.2, 2.2);
    TCanvas *c_fit = histo_fit_ratio(low_m_histo, pdf_low_m, ratio_low_m, "low_m_extr", "#splitline{ALICE Simulation, pp #sqrt{s} = 13 TeV}{POWHEG+PYTHIA6 simulation}", kTRUE);

    TF1 *pdf_high_m = new TF1("pdf_high_m", FuncPtMass, info.HighM_cut, 30, 4);
    pdf_high_m->SetTitle("Fit");
    // pdf_high_m->FixParameter(0, B_DimuMass->getVal());
    // pdf_high_m->FixParameter(1, n1_DimuMass->getVal());
    // pdf_high_m->FixParameter(2, n2_DimuMass->getVal());

    // pdf_high_m->SetParError(0, B_DimuMass->getError());
    // pdf_high_m->SetParError(1, n1_DimuMass->getError());
    // pdf_high_m->SetParError(2, n2_DimuMass->getError());

    pdf_high_m->FixParameter(0, 1.60856e+01);
    pdf_high_m->FixParameter(1, 2.70776e-01);
    pdf_high_m->FixParameter(2, 4.55889e+01);

    // f_M_PYTHIA6->SetParameter(3, 1);
    pdf_high_m->SetLineColor(kMagenta + 2);
    pdf_high_m->SetLineStyle(kDashed);
    pdf_high_m->SetLineWidth(3);

    TH1F *high_m_histo = new TH1F("high_m_histo", "high_m_histo", 28, info.HighM_cut, 30);
    high_m_tree_data_cutted->Draw("high_m>>high_m_histo", "", "goff");
    high_m_histo->Scale(1., "width");
    high_m_histo->SetMarkerSize(1.5);
    high_m_histo->SetMarkerStyle(20);
    high_m_histo->SetMarkerColor(kBlack);
    high_m_histo->Fit(pdf_high_m, "R0I");

    TH1F *ratio_high_m = (TH1F *)high_m_histo->Clone("ratio_highm_fit");
    ratio_high_m->Divide(pdf_high_m);
    ratio_high_m->GetYaxis()->SetTitle("MC/Fit");
    ratio_high_m->GetYaxis()->SetRangeUser(-0.2, 2.2);
    c_fit = histo_fit_ratio(high_m_histo, pdf_high_m, ratio_high_m, "high_m_extr", "#splitline{ALICE Simulation, pp #sqrt{s} = 13 TeV}{POWHEG+PYTHIA6 simulation}", kTRUE);

    TCanvas *c = new TCanvas("c", "c", 1200, 1200);
    c->Divide(1, 2);
    c->cd(1);
    gPad->SetLogy();
    TH1F *m_histo = new TH1F("m_histo", "m_histo", 260, 4, 30);
    m_histo->GetYaxis()->SetRangeUser(1.2, 1.2e+6);
    m_histo->DrawCopy();
    low_m_histo->Draw("PEsame");
    pdf_low_m->Draw("same");
    high_m_histo->Draw("PEsame");
    pdf_high_m->Draw("same");

    c->cd(2);
    m_histo->GetYaxis()->SetRangeUser(0.2, 2.2);
    m_histo->DrawCopy();
    ratio_low_m->Draw("PESAME");
    ratio_high_m->Draw("PESAME");

    return;

    gROOT->cd();

    Color_t color = {kBlack};
    Color_t fillcolor = {kBlack};

    return;
    TString Pt_Factory_Info;
    // Pt_Factory_Info.Form("PtMassExpPdf::pdfDimuPtFromData(%s[%d,%d], B_DimuPtFromData[2.85,1e-4,1e+4], n1_DimuPtFromData[2.81,1e-4,1e+4], n2_DimuPtFromData[2.43,1e-4,1e+4])", pt->GetName(), info.Low_Pt, info.High_Pt);

    // w->factory(Pt_Factory_Info);

    // RooAbsPdf *pdfDimuPt = w->pdf("pdfDimuPtFromData");
    // auto result1 = pdfDimuPt->fitTo(*Pt_Dimu, Minimizer("Minuit2"), Save(), SumW2Error(true));
    // RooPlot *frameDimuPt = pt->frame(Title("frameDimuMass_Data"));
    // frameDimuPt->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    // Pt_Dimu->plotOn(frameDimuPt, Name("Pt_Data"), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color), LineColor(color), LineWidth(2), DrawOption("PEZ"), Binning(info.Pt_Binning));
    // pdfDimuPt->plotOn(frameDimuPt, Name("pdfDimuPtFromData"), LineStyle(kSolid), LineColor(color));
    // TH1 *pt_histo = Pt_Dimu->createHistogram("h_Pt_Dimu_Data", *pt, Binning(info.Pt_Binning, info.Low_Pt, info.High_Pt));

    // RooArgSet *pt_model = static_cast<RooArgSet *>(RooArgSet(*pdfDimuPt).snapshot(true)); // True means copy the PDF and everything it depends on
    // auto &pt_modelcopied = static_cast<RooAbsPdf &>((*pt_model)["pdfDimuPtFromData"]);    // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    // RooArgSet *pt_modelobs = pt_modelcopied.getObservables(*Pt_Dimu);
    // RooArgSet *pt_modelPars = pt_modelcopied.getParameters(*pt_modelobs);
    // TF1 *pt_Func = pt_modelcopied.asTF(*pt_modelobs, *pt_modelPars, *pt);
    // pt_Func->SetName(Form("pdfDimuPtFrom%s_M_%d_%d_Pt_%d_%d", "Data", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    // // pt_Func->Write(0, 2, 0);
    // TCanvas *pt_canvas = printMC_ratio("pt_canvas_Data", frameDimuPt, pt_histo, pt_Func, color, info.Low_Pt, info.High_Pt);

    // auto result2 = high_pdfDimuM->fitTo(*high_M_Dimu, Minimizer("Minuit2"), Range("pluto"), Save(), SumW2Error(true));
    // RooPlot *frameDimuMass = m.frame(Title("frameDimuMass_Data"));
    // frameDimuMass->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");
    // high_M_Dimu->plotOn(frameDimuMass, Name("M_Data"), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color), LineColor(color), LineWidth(2), DrawOption("PEZ"), Binning(info.Mass_Binning));
    // high_pdfDimuM->plotOn(frameDimuMass, Name("pdfDimuMFromData"), LineStyle(kSolid), LineColor(color));
    // TH1 *m_histo = high_M_Dimu->createHistogram("h_M_Dimu_Data", m, Binning(info.Mass_Binning, info.Low_Mass, info.High_Mass));

    // RooArgSet *m_model = static_cast<RooArgSet *>(RooArgSet(*high_pdfDimuM).snapshot(true)); // True means copy the PDF and everything it depends on
    // auto &m_modelcopied = static_cast<RooAbsPdf &>((*m_model)["pdfDimuMassFromData"]);       // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    // RooArgSet *m_modelobs = m_modelcopied.getObservables(*high_M_Dimu);
    // RooArgSet *m_modelPars = m_modelcopied.getParameters(*m_modelobs);
    // TF1 *m_Func = m_modelcopied.asTF(*m_modelobs, *m_modelPars, m);
    // m_Func->SetName(Form("pdfDimuMassFrom%s_M_%d_%d_Pt_%d_%d", "Data", info.Low_Mass, info.High_Mass, info.Low_Pt, info.High_Pt));
    // // m_Func->Write(0, 2, 0);

    // TCanvas *m_canvas = printMC_ratio("m_canvas_Data", frameDimuMass, m_histo, m_Func, color, info.Low_Mass, info.High_Mass);
}

TCanvas *printMC_ratio_BINNED(TString name, TString Title, TH1 *data, TF1 *pdf, Color_t color)
{
    opt info;
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas(name, name, 900, 1000);
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
    TString var_str;
    str.Form("%s", data->GetName());

    if (str.Contains("Mass") || str.Contains("h_M") || str.Contains("mass"))
    {
        data->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#plus}} (GeV/#it{c}^{2})");
        data->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#plus}} (GeV/#it{c}^{2})^{-1}");
        var_str.Form("Mass");
    }
    else
    {
        data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        data->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
        var_str.Form("Pt");
    }

    data->SetTitle("");
    data->GetXaxis()->SetTitleOffset(1.3);
    data->GetXaxis()->SetTitleSize(0.0475);
    data->GetXaxis()->SetLabelSize(0.045);

    data->GetYaxis()->SetNdivisions(505);
    data->GetYaxis()->SetTitleOffset(0.9);
    data->GetYaxis()->SetTitleSize(0.065);
    data->GetYaxis()->SetLabelSize(0.055);
    data->Draw();
    pdf->SetLineColor(data->GetLineColor());
    pdf->SetLineWidth(3.);
    pdf->Draw("same");
    TLegend *legend = new TLegend(0.375, 0.15, 0.55, 0.4);

    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0525);
    legend->SetHeader("MC");
    legend->SetTextAlign(11);
    legend->AddEntry(data, " ", "PE");

    TLegend *pdf_legend = new TLegend(0.45, 0.15, 0.635, 0.4);
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
    if (str.Contains("LF_HF"))
        letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow LF,HF");

    letexTitle->SetTextSize(0.055);
    letexTitle->DrawLatex(0.355, 0.88, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    if (info.Generator.Contains("POWHEG"))
        letexTitle->DrawLatex(0.355, 0.81, "POWHEG+PYTHIA6, N_{ev} = 2 #upoint 10^{9}");
    else
        letexTitle->DrawLatex(0.355, 0.81, "PYTHIA8, N_{ev} = 2.1 #upoint 10^{8}");

    letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{#eta}_{#mu} < 4.0");
    letexTitle->DrawLatex(0.355, 0.66, Title);

    letexTitle->DrawLatex(0.5, 0.60, Form("B: %0.3f #pm %0.3f", pdf->GetParameter(0), pdf->GetParError(0)));
    letexTitle->DrawLatex(0.5, 0.55, Form("n1: %0.3f #pm %0.3f", pdf->GetParameter(1), pdf->GetParError(1)));
    letexTitle->DrawLatex(0.5, 0.50, Form("n2: %0.3f #pm %0.3f", pdf->GetParameter(2), pdf->GetParError(2)));

    pad2->cd();
    pad2->SetTicks();
    TLine *l = new TLine(data->GetBinLowEdge(1), 1.0, data->GetBinLowEdge(data->GetNbinsX() + 1), 1.0);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->SetLineColor(kRed);

    TLine *l1 = new TLine(data->GetBinLowEdge(1), 1.5, data->GetBinLowEdge(data->GetNbinsX() + 1), 1.5);
    l1->SetLineWidth(2);
    l1->SetLineStyle(9);
    l1->SetLineColor(kGray + 2);

    TLine *l2 = new TLine(data->GetBinLowEdge(1), 0.5, data->GetBinLowEdge(data->GetNbinsX() + 1), 0.5);
    l2->SetLineWidth(2);
    l2->SetLineStyle(9);
    l2->SetLineColor(kGray + 2);

    TH1D *c_data = (TH1D *)data->Clone(Form("c_data_%s", data->GetName()));

    c_data->SetTitle("");
    c_data->GetYaxis()->SetTitle(Form("#frac{MC}{PDF}"));
    c_data->GetYaxis()->CenterTitle();
    c_data->GetYaxis()->SetNdivisions(504);
    c_data->GetYaxis()->SetTitleSize(0.08);
    c_data->GetYaxis()->SetLabelOffset(0.02);
    c_data->GetYaxis()->SetLabelSize(0.1);

    c_data->GetXaxis()->SetTitleSize(0.1);
    c_data->GetXaxis()->SetTitleOffset(1.1);
    c_data->GetXaxis()->SetLabelSize(0.1);
    c_data->GetXaxis()->SetTitle(c_data->GetXaxis()->GetTitle());

    c_data->SetLineColor(kBlack);
    c_data->SetMarkerColor(kBlack);
    c_data->SetMarkerStyle(20);
    c_data->Divide(pdf);
    c_data->GetYaxis()->SetRangeUser(0.2, 2.2);

    c_data->Draw();
    l->Draw();
    l1->Draw();
    l2->Draw();

    // canvas->ls();
    return canvas;
}

TCanvas *printMC_ratio(TString name, TString Title, RooPlot *frame, RooArgSet *param, TH1 *data, TF1 *pdf, Color_t color, Int_t minx = 0, Int_t max_x = 30)
{
    opt info;
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas(name, name, 900, 1000);
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
    TString var_str;
    str.Form("%s", data->GetTitle());

    if (str.Contains("h_M"))
    {
        data->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#plus}} (GeV/#it{c}^{2})");
        data->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#plus}} (GeV/#it{c}^{2})^{-1}");
        var_str.Form("Mass");
    }
    else
    {
        data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        data->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
        var_str.Form("Pt");
    }

    frame->SetTitle("");
    frame->GetXaxis()->SetTitleOffset(1.3);
    frame->GetXaxis()->SetTitleSize(0.0475);
    frame->GetXaxis()->SetLabelSize(0.045);

    frame->GetYaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetTitleOffset(0.9);
    frame->GetYaxis()->SetTitleSize(0.065);
    frame->GetYaxis()->SetLabelSize(0.055);
    frame->Draw();
    Double_t Integral = 0.;
    for (Int_t i = 0; i < data->GetNbinsX(); i++)
    {
        Integral = Integral + data->GetBinContent(i + 1);
    }
    data->Scale(1. / Integral, "width");

    TLegend *legend = new TLegend(0.315, 0.15, 0.55, 0.4);

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

    letexTitle->SetTextSize(0.055);
    letexTitle->DrawLatex(0.355, 0.88, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    if (info.Generator.Contains("POWHEG"))
        letexTitle->DrawLatex(0.355, 0.81, "POWHEG+PYTHIA6, N_{ev} = 2 #upoint 10^{9}");
    else
        letexTitle->DrawLatex(0.355, 0.81, "PYTHIA8, N_{ev} = 2.1 #upoint 10^{8}");

    letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{#eta}_{#mu} < 4.0");
    letexTitle->DrawLatex(0.355, 0.66, Title);

    if (param != nullptr)
    {
        letexTitle->SetTextSize(0.045);
        param->Print("v");
        double_t var;
        double_t error;

        var = static_cast<RooRealVar &>((*param)[TString::Format("B_Dimu%sFrom%s", var_str.Data(), param->GetName())]).getVal();
        error = static_cast<RooRealVar &>((*param)[TString::Format("B_Dimu%sFrom%s", var_str.Data(), param->GetName())]).getError();
        letexTitle->DrawLatex(0.5, 0.60, Form("B: %0.3f #pm %0.3f", var, error));
        var = static_cast<RooRealVar &>((*param)[TString::Format("n1_Dimu%sFrom%s", var_str.Data(), param->GetName())]).getVal();
        error = static_cast<RooRealVar &>((*param)[TString::Format("n1_Dimu%sFrom%s", var_str.Data(), param->GetName())]).getError();
        letexTitle->DrawLatex(0.5, 0.55, Form("n1: %0.3f #pm %0.3f", var, error));
        var = static_cast<RooRealVar &>((*param)[TString::Format("n2_Dimu%sFrom%s", var_str.Data(), param->GetName())]).getVal();
        error = static_cast<RooRealVar &>((*param)[TString::Format("n2_Dimu%sFrom%s", var_str.Data(), param->GetName())]).getError();
        letexTitle->DrawLatex(0.5, 0.50, Form("n2: %0.3f #pm %0.3f", var, error));
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

    TH1D *c_data = (TH1D *)data->Clone(Form("c_data_%s", data->GetName()));

    c_data->SetTitle("");
    c_data->GetYaxis()->SetTitle(Form("#frac{MC}{PDF}"));
    c_data->GetYaxis()->CenterTitle();
    c_data->GetYaxis()->SetNdivisions(504);
    c_data->GetYaxis()->SetTitleSize(0.08);
    c_data->GetYaxis()->SetLabelOffset(0.02);
    c_data->GetYaxis()->SetLabelSize(0.1);

    c_data->GetXaxis()->SetTitleSize(0.1);
    c_data->GetXaxis()->SetTitleOffset(1.1);
    c_data->GetXaxis()->SetLabelSize(0.1);
    c_data->GetXaxis()->SetTitle(c_data->GetXaxis()->GetTitle());

    c_data->SetLineColor(kBlack);
    c_data->SetMarkerColor(kBlack);
    c_data->SetMarkerStyle(20);
    c_data->Divide(pdf);
    c_data->GetYaxis()->SetRangeUser(0.2, 2.2);
    if (str.Contains("h_M"))
    {
        TH1F *blank = (TH1F *)c_data->Clone("blank");
        blank->Reset();
        blank->Draw();
        TH1F *copy = (TH1F *)c_data->Clone("copy");
        copy->GetXaxis()->SetRangeUser(4, 8);
        copy->DrawCopy("same");
        copy->GetXaxis()->SetRangeUser(11, 30);
        copy->DrawCopy("same");
    }
    else
        c_data->Draw();
    l->Draw();
    l1->Draw();
    l2->Draw();

    // canvas->ls();
    return canvas;
}

TCanvas *printRooPlot_ratio(RooPlot *frame, Bool_t norm, RooFitResult *r, Int_t choice, TString roohist_name, TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x, Double_t N_HFMixed, Double_t N_LF_HFMixed)
{
    Bool_t Debug = kFALSE;
    opt info;
    const RooArgList &fitParams = r->floatParsFinal();

    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 1200);
    canvas->SetTicks();
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(0);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetFrameBorderMode(0);

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
    // if // {
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

    if (roohist_name.Contains("pt"))
        fit_legend->AddEntry("pdfpt", " ", "L");
    else if (roohist_name.Contains("mass"))
        fit_legend->AddEntry("pdfmass", " ", "L");

    TLegend *pdf_legend = new TLegend(0.65, 0.615, 0.9, 0.915);
    pdf_legend->SetTextAlign(11);
    pdf_legend->SetFillStyle(0);
    pdf_legend->SetBorderSize(0);
    pdf_legend->SetTextSize(0.0425);
    pdf_legend->SetTextAlign(12);
    if (info.Generator.Contains("POWHEG"))
        pdf_legend->SetHeader(Form("#splitline{POWHEG+PYTHIA6 PDF}{%s}", info.systematic.Data()));
    else
        pdf_legend->SetHeader(Form("#splitline{PYTHIA8 Monash PDF}{%s}", info.systematic.Data()));

    // pdf_legend->AddEntry(TString::Format("%s%s", roohist_name.Data(), info.Name_DimuSel[0].Data()), Form("#mu^{#plus}#mu^{#minus} %s", info.info_label[0].Data()), "L");
    // pdf_legend->AddEntry(TString::Format("%s%s", roohist_name.Data(), info.Name_DimuSel[1].Data()), Form("#mu^{#plus}#mu^{#minus} %s", info.info_label[1].Data()), "L");
    // pdf_legend->AddEntry(TString::Format("%s%s", roohist_name.Data(), info.Name_DimuSel[3].Data()), Form("#mu^{#plus}#mu^{#minus} %s", info.info_label[3].Data()), "L");

    for (Int_t i = 0; i < fitParams.getSize(); i++)
    {
        auto &fitPar = (RooRealVar &)fitParams[i];
        cout << fitPar.GetName() << endl;
        if (fitPar.getVal() != 0)
            pdf_legend->AddEntry(TString::Format("%s%s", roohist_name.Data(), info.Name_DimuSel[i].Data()), Form("#mu^{#plus}#mu^{#minus} %s", info.info_label[i].Data()), "L");
    }

    pdf_legend->AddEntry(TString::Format("%s%s", roohist_name.Data(), info.Name_DimuSel[3].Data()), Form("#mu^{#plus}#mu^{#minus} %s", info.info_label[3].Data()), "L");
    if (N_LF_HFMixed != 0)
        pdf_legend->AddEntry(TString::Format("%s%s", roohist_name.Data(), info.Name_DimuSel[4].Data()), Form("#mu^{#plus}#mu^{#minus} %s", info.info_label[4].Data()), "L");

    // legend->Draw();
    // fit_legend->Draw();
    pdf_legend->Draw();

    TLatex *letexTitle = new TLatex();
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);
    letexTitle->SetTextSize(0.0475);
    letexTitle->DrawLatex(0.175, 0.875, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    if (info.stat_Data.Contains("LHC18p"))
        letexTitle->DrawLatex(0.175, 0.805, "LHC18p period");
    else if (info.stat_Data.Contains("Run2"))
        letexTitle->DrawLatex(0.175, 0.805, "Run 2");
    letexTitle->SetTextSize(0.0425);
    vector<double> fit_result;
    if (Debug)
    {
        Double_t start_y_latex = 0.875;
        Int_t i_par = 0;
        while (i_par < fitParams.getSize())
        {
            auto &fitPar = (RooRealVar &)fitParams[i_par];
            fit_result.push_back(fitPar.getVal());
            info.info_label[i_par].Append(Form("} = %0.1f #pm %0.1f", fitPar.getVal(), fitPar.getError()));
            cout << Form("#it{N}_{#mu^{#plus}#mu^{#minus}%s", info.info_label[i_par].Data()) << endl;
            letexTitle->DrawLatex(0.625, start_y_latex - 0.08 * i_par, Form("#it{N}_{#mu^{#plus}#mu^{#minus}%s", info.info_label[i_par].Data()));
            i_par++;
        }

        if (N_HFMixed != 0)
        {
            info.info_label[3].Append(Form("} = %0.1f", N_HFMixed));
            cout << Form("#it{N}^{fixed}_{#mu^{#plus}#mu^{#minus}%s", info.info_label[3].Data()) << endl;
            letexTitle->DrawLatex(0.625, start_y_latex - 0.08 * i_par, Form("#it{N}^{fixed}_{#mu^{#plus}#mu^{#minus}%s", info.info_label[3].Data()));
            i_par++;
        }
        if (N_LF_HFMixed != 0)
        {
            info.info_label[4].Append(Form("} = %0.1f", N_LF_HFMixed));
            cout << Form("#it{N}^{fixed}_{#mu^{#plus}#mu^{#minus}%s", info.info_label[4].Data()) << endl;
            letexTitle->DrawLatex(0.625, start_y_latex - 0.08 * i_par, Form("#it{N}^{fixed}_{#mu^{#plus}#mu^{#minus}%s", info.info_label[4].Data()));
            i_par++;
        }
    }

    if (roohist_name.Contains("pt"))
    {
        letexTitle->DrawLatex(0.175, 0.735, "Reconstructed #mu^{#plus}#mu^{#minus}");
        letexTitle->DrawLatex(0.175, 0.665, Form("%d < #it{m}_{#mu^{#plus}#mu^{#minus}} < %d GeV/#it{c}^{2}, 2.5 < #it{y}_{#mu} < 4.0", info.Low_Mass, info.High_Mass));
    }
    else
        letexTitle->DrawLatex(0.175, 0.735, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu} < 4.0");

    pad2->cd();
    pad2->SetTicks();
    TLine *l = new TLine(minx, 1.0, max_x, 1.0);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->SetLineColor(kRed);

    TLine *l1 = new TLine(minx, 0.8, max_x, 0.8);
    l1->SetLineWidth(2);
    l1->SetLineStyle(9);
    l1->SetLineColor(kGray + 2);

    TLine *l2 = new TLine(minx, 1.2, max_x, 1.2);
    l2->SetLineWidth(2);
    l2->SetLineStyle(9);
    l2->SetLineColor(kGray + 2);

    TH2D *h_grid_ratio = new TH2D("h_grid", "", 100, minx, max_x, 100, 0.7, 1.3);
    h_grid_ratio->SetName(Form("%s_ratiodatafit", roohist_name.Data()));
    h_grid_ratio->SetTitle("");
    h_grid_ratio->GetYaxis()->SetTitle(Form("#frac{Data}{Cocktail}"));
    h_grid_ratio->GetYaxis()->CenterTitle();
    h_grid_ratio->GetYaxis()->SetNdivisions(504);
    h_grid_ratio->GetYaxis()->SetTitleSize(0.08);
    h_grid_ratio->GetYaxis()->SetTitleOffset(0.8);
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
    // data->Scale(1. / data->Integral(), "width");
    data->Divide(pdf);

    h_grid_ratio->Draw();
    data->Draw("PESAME");
    l->Draw();
    l1->Draw();
    l2->Draw();

    return canvas;
}

void unbinned_fit_data_multiregion(Int_t Low_Mass = 4, Int_t High_Mass = 9, Int_t Low_Pt = 0, Int_t High_Pt = 10)
{
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");

    Int_t Mass_Binning = (High_Mass - Low_Mass) / 0.5;
    Int_t Pt_Binning = (High_Pt - Low_Pt) / 0.5;
    const Int_t n_DiMuSelection = 3;
    TString Name_DimuSel[n_DiMuSelection] = {"Charm", "Beauty", "HF_Mixed"};
    Color_t color[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kAzure + 9};

    RooRealVar *Low_m = new RooRealVar("Low_m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", 4, 9);
    RooRealVar *Low_pt = new RooRealVar("Low_pt", "#it{p}_{T} (GeV/#it{c})", 0, 10);

    RooRealVar *High_m = new RooRealVar("High_m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", 11, 30);
    RooRealVar *High_pt = new RooRealVar("High_pt", "#it{p}_{T} (GeV/#it{c})", 0, 30);
    // m->setBins(Mass_Binning);
    // pt->setBins(Pt_Binning);

    RooCategory sample("sample", "sample");
    sample.defineType("Low_mass");
    sample.defineType("Low_transversemomentum");
    sample.defineType("High_mass");
    sample.defineType("High_transversemomentum");

    TFile *fIn_data = new TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/TreeResults_merged.root", "READ");
    fIn_data->ls();

    // Taking data saved in tree
    TTree *tree_data = (TTree *)fIn_data->Get("rec_data_tree");
    gROOT->cd();

    TTree *LowM_LowPt_tree_data_cutted = (TTree *)tree_data->CopyTree(Form("(m>%d && m<%d) && (pt > %d && pt <%d)", 4, 9, 0, 10));
    LowM_LowPt_tree_data_cutted->SetName("LowM_LowPt_tree_data_cutted");
    LowM_LowPt_tree_data_cutted->GetBranch("m")->SetName("Low_m");
    LowM_LowPt_tree_data_cutted->GetBranch("pt")->SetName("Low_pt");
    // LowM_LowPt_tree_data_cutted->Draw("Low_m");
    RooDataSet *LowM_LowPt_unbinned_M_Dimu_data = new RooDataSet("LowM_Dimu_data", "LowM_Dimu_data", RooArgSet(*Low_m), Import(*LowM_LowPt_tree_data_cutted));
    RooDataSet *LowM_LowPt_unbinned_Pt_Dimu_data = new RooDataSet("LowPt_Dimu_data", "LowPt_Dimu_data", RooArgSet(*Low_pt), Import(*LowM_LowPt_tree_data_cutted));

    TTree *HighM_HighPt_tree_data_cutted = (TTree *)tree_data->CopyTree(Form("(m>%d && m<%d) && (pt > %d && pt <%d)", 11, 30, 0, 30));
    HighM_HighPt_tree_data_cutted->SetName("HighM_HighPt_tree_data_cutted");
    HighM_HighPt_tree_data_cutted->GetBranch("m")->SetName("High_m");
    HighM_HighPt_tree_data_cutted->GetBranch("pt")->SetName("High_pt");
    RooDataSet *HighM_HighPt_unbinned_M_Dimu_data = new RooDataSet("HighM_Dimu_data", "HighM_Dimu_data", RooArgSet(*High_m), Import(*HighM_HighPt_tree_data_cutted));
    RooDataSet *HighM_HighPt_unbinned_Pt_Dimu_data = new RooDataSet("HighPt_Dimu_data", "HighPt_Dimu_data", RooArgSet(*High_pt), Import(*HighM_HighPt_tree_data_cutted));
    RooDataSet *unbinned_combData_set = new RooDataSet("combData", "combined data", RooArgSet(*Low_m, *Low_pt, *High_m, *High_pt), Index(sample), Import("Low_mass", *LowM_LowPt_unbinned_M_Dimu_data), Import("Low_transversemomentum", *LowM_LowPt_unbinned_Pt_Dimu_data), Import("High_mass", *HighM_HighPt_unbinned_M_Dimu_data), Import("High_transversemomentum", *HighM_HighPt_unbinned_Pt_Dimu_data));

    // RooDataSet *unbinned_combData_set = new RooDataSet("combData", "combined data", RooArgSet(*High_m, *High_pt), Index(sample), Import("High_mass", *HighM_HighPt_unbinned_M_Dimu_data), Import("High_transversemomentum", *HighM_HighPt_unbinned_Pt_Dimu_data));
    //----//
    RooPlot *frame = High_pt->frame(Title("m_frame"));
    unbinned_combData_set->plotOn(frame, Name("combDatapt"), Cut("sample==sample::High_transversemomentum"), DrawOption("PEZ"), MarkerColor(kRed));
    // unbinned_combData_set->plotOn(frame, Name("combDatapt"), Cut("sample==sample::High_mass"), DrawOption("PEZ"));
    frame->Draw();

    TFile *fIn = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/Powheg_pdfMC_unbinnedprova.root", "READ");

    RooWorkspace *LowM_LowPt_w = (RooWorkspace *)fIn->Get(Form("w_M_%d_%d_Pt_%d_%d", 4, 9, 0, 10));
    LowM_LowPt_w->Print("s");

    RooWorkspace *HighM_HighPt_w = (RooWorkspace *)fIn->Get(Form("w_M_%d_%d_Pt_%d_%d", 11, 30, 0, 30));
    HighM_HighPt_w->Print("s");

    RooRealVar *LowM_LowPt_B_DimuMass[n_DiMuSelection];
    RooRealVar *LowM_LowPt_n1_DimuMass[n_DiMuSelection];
    RooRealVar *LowM_LowPt_n2_DimuMass[n_DiMuSelection];
    RooRealVar *LowM_LowPt_B_DimuPt[n_DiMuSelection];
    RooRealVar *LowM_LowPt_n1_DimuPt[n_DiMuSelection];
    RooRealVar *LowM_LowPt_n2_DimuPt[n_DiMuSelection];
    RooAbsPdf *LowM_LowPt_pdfDimuMass[n_DiMuSelection];
    RooAbsPdf *LowM_LowPt_pdfDimuPt[n_DiMuSelection];

    RooRealVar *HighM_HighPt_B_DimuMass[n_DiMuSelection];
    RooRealVar *HighM_HighPt_n1_DimuMass[n_DiMuSelection];
    RooRealVar *HighM_HighPt_n2_DimuMass[n_DiMuSelection];
    RooRealVar *HighM_HighPt_B_DimuPt[n_DiMuSelection];
    RooRealVar *HighM_HighPt_n1_DimuPt[n_DiMuSelection];
    RooRealVar *HighM_HighPt_n2_DimuPt[n_DiMuSelection];
    RooAbsPdf *HighM_HighPt_pdfDimuMass[n_DiMuSelection];
    RooAbsPdf *HighM_HighPt_pdfDimuPt[n_DiMuSelection];

    for (Int_t i_DiMu_Sel = 0; i_DiMu_Sel < n_DiMuSelection; i_DiMu_Sel++)
    {
        LowM_LowPt_B_DimuMass[i_DiMu_Sel] = LowM_LowPt_w->var(Form("B_DimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        LowM_LowPt_B_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        LowM_LowPt_n1_DimuMass[i_DiMu_Sel] = LowM_LowPt_w->var(Form("n1_DimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        LowM_LowPt_n1_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        LowM_LowPt_n2_DimuMass[i_DiMu_Sel] = LowM_LowPt_w->var(Form("n2_DimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        LowM_LowPt_n2_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        LowM_LowPt_B_DimuPt[i_DiMu_Sel] = LowM_LowPt_w->var(Form("B_DimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        LowM_LowPt_B_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        LowM_LowPt_n1_DimuPt[i_DiMu_Sel] = LowM_LowPt_w->var(Form("n1_DimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        LowM_LowPt_n1_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        LowM_LowPt_n2_DimuPt[i_DiMu_Sel] = LowM_LowPt_w->var(Form("n2_DimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        LowM_LowPt_n2_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        LowM_LowPt_pdfDimuMass[i_DiMu_Sel] = LowM_LowPt_w->pdf(Form("pdfDimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        LowM_LowPt_pdfDimuMass[i_DiMu_Sel]->SetName(Form("LowM_LowPt_%s", LowM_LowPt_pdfDimuMass[i_DiMu_Sel]->GetName()));
        LowM_LowPt_pdfDimuPt[i_DiMu_Sel] = LowM_LowPt_w->pdf(Form("pdfDimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        LowM_LowPt_pdfDimuPt[i_DiMu_Sel]->SetName(Form("LowM_LowPt_%s", LowM_LowPt_pdfDimuPt[i_DiMu_Sel]->GetName()));

        HighM_HighPt_B_DimuMass[i_DiMu_Sel] = HighM_HighPt_w->var(Form("B_DimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        HighM_HighPt_B_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        HighM_HighPt_n1_DimuMass[i_DiMu_Sel] = HighM_HighPt_w->var(Form("n1_DimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        HighM_HighPt_n1_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        HighM_HighPt_n2_DimuMass[i_DiMu_Sel] = HighM_HighPt_w->var(Form("n2_DimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        HighM_HighPt_n2_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        HighM_HighPt_B_DimuPt[i_DiMu_Sel] = HighM_HighPt_w->var(Form("B_DimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        HighM_HighPt_B_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        HighM_HighPt_n1_DimuPt[i_DiMu_Sel] = HighM_HighPt_w->var(Form("n1_DimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        HighM_HighPt_n1_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        HighM_HighPt_n2_DimuPt[i_DiMu_Sel] = HighM_HighPt_w->var(Form("n2_DimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        HighM_HighPt_n2_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        HighM_HighPt_pdfDimuMass[i_DiMu_Sel] = HighM_HighPt_w->pdf(Form("pdfDimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        HighM_HighPt_pdfDimuMass[i_DiMu_Sel]->SetName(Form("HighM_HighPt_%s", HighM_HighPt_pdfDimuMass[i_DiMu_Sel]->GetName()));
        HighM_HighPt_pdfDimuPt[i_DiMu_Sel] = HighM_HighPt_w->pdf(Form("pdfDimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        HighM_HighPt_pdfDimuPt[i_DiMu_Sel]->SetName(Form("HighM_HighPt_%s", HighM_HighPt_pdfDimuPt[i_DiMu_Sel]->GetName()));
    }

    // Fit with fraction
    // RooRealVar *normForC = new RooRealVar("fr_charm_output", "fraction dimuon from c", 0.35, 0., 1.);
    // RooRealVar *normForB = new RooRealVar("fr_beauty_output", "fraction dimuon from b", 0.605, 0., 1.);

    // RooRealVar *normForMixed = new RooRealVar("fr_mixed_output", "fraction dimuon from c", 0.036);
    // normForMixed->setConstant(kTRUE);

    RooRealVar *normForC = new RooRealVar("n_charm_output", "number dimuon from c", 28440, 0, 200000);
    RooRealVar *normForB = new RooRealVar("n_beauty_output", "number dimuon from b", 48000, 0, 200000);

    RooRealVar *High_normForC = new RooRealVar("High_M_n_charm_output", "number dimuon from c", 28440, 0, 200000);
    RooRealVar *High_normForB = new RooRealVar("High_M_n_beauty_output", "number dimuon from b", 48000, 0, 200000);

    RooRealVar *LowMnormForC = new RooRealVar("Low_M_n_charm_output", "number dimuon from c", 28440, 0, 200000);
    RooRealVar *LowMnormForB = new RooRealVar("Low_M_n_beauty_output", "number dimuon from b", 48000, 0, 200000);

    RooRealVar *LowM_normForMixed = new RooRealVar("n_mixed_output", "number dimuon from b,c", (2.2 / 100) * LowM_LowPt_tree_data_cutted->GetEntries());
    LowM_normForMixed->setConstant(kTRUE);

    RooRealVar *HighM_normForMixed = new RooRealVar("n_mixed_output", "number dimuon from b,c", (2.2 / 100) * HighM_HighPt_tree_data_cutted->GetEntries());
    HighM_normForMixed->setConstant(kTRUE);

    RooAddPdf *Low_m_model = new RooAddPdf("Low_m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*LowM_LowPt_pdfDimuMass[0], *LowM_LowPt_pdfDimuMass[1], *LowM_LowPt_pdfDimuMass[2]), RooArgList(*LowMnormForC, *LowMnormForB, *LowM_normForMixed));
    RooAddPdf *Low_pt_model = new RooAddPdf("Low_pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*LowM_LowPt_pdfDimuPt[0], *LowM_LowPt_pdfDimuPt[1], *LowM_LowPt_pdfDimuPt[2]), RooArgList(*LowMnormForC, *LowMnormForB, *LowM_normForMixed));

    RooAddPdf *High_m_model = new RooAddPdf("High_m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*HighM_HighPt_pdfDimuMass[0], *HighM_HighPt_pdfDimuMass[1], *HighM_HighPt_pdfDimuMass[2]), RooArgList(*High_normForC, *High_normForB, *HighM_normForMixed));
    RooAddPdf *High_pt_model = new RooAddPdf("High_pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*HighM_HighPt_pdfDimuPt[0], *HighM_HighPt_pdfDimuPt[1], *HighM_HighPt_pdfDimuPt[2]), RooArgList(*High_normForC, *High_normForB, *HighM_normForMixed));

    // m->setRange("ino", 4, 9);
    // pt->setRange("ino", 0, 10);

    // m->setRange("paper", 11, 30);
    // pt->setRange("paper", 12, 30);
    RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
    simPdf.addPdf(*Low_m_model, "Low_mass");
    simPdf.addPdf(*Low_pt_model, "Low_transversemomentum");
    simPdf.addPdf(*High_m_model, "High_mass");
    simPdf.addPdf(*High_pt_model, "High_transversemomentum");
    simPdf.Print("t");

    RooFitResult *r = simPdf.fitTo(*unbinned_combData_set, Minimizer("Minuit2"), Save(), SumW2Error(true));
    // simPdf.fixCoefRange("norm");

    RooPlot *m_frame = High_m->frame(Title("m_frame"));
    m_frame->SetTitle(" ");
    m_frame->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
    m_frame->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");
    unbinned_combData_set->plotOn(m_frame, Name("combDatamass"), Cut("sample==sample::High_mass"), DrawOption("PEZ"), Binning(10));
    // unbinned_combData_set->plotOn(m_frame, Name("combDatamass"), Cut("sample==sample::High_mass"), DrawOption("PEZ"), MarkerColor(kRed));

    RooPlot *pt_frame = High_pt->frame(Title("pt_frame"));
    pt_frame->SetTitle(" ");
    pt_frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    pt_frame->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

    simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "High_mass"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
    unbinned_combData_set->plotOn(pt_frame, Name("combDatapt"), Cut("sample==sample::High_transversemomentum"), DrawOption("PEZ"), Binning(10));
    simPdf.plotOn(pt_frame, Name("pdfpt"), Slice(sample, "High_transversemomentum"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
    for (Int_t i_DiMu_Sel = 0; i_DiMu_Sel < n_DiMuSelection; i_DiMu_Sel++)
    {
        simPdf.plotOn(m_frame, Name(Form("pdfmass%s", Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "High_mass"), Components(HighM_HighPt_pdfDimuMass[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(color[i_DiMu_Sel]), LineWidth(5));

        simPdf.plotOn(pt_frame, Name(Form("pdfpt%s", Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "High_transversemomentum"), Components(HighM_HighPt_pdfDimuPt[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(color[i_DiMu_Sel]), LineWidth(5));
    }
    m_frame->Draw();
    new TCanvas();
    pt_frame->Draw();
}

void conv_DY_cs()
{
    TFile *fIn = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC18p_DY_Version_5/LHC18p_DY_MCDimuHFTree_294009.root", "READ");

    TTree *MCTree = (TTree *)fIn->Get("MCTree");
    TH1F *h_NDY_event = new TH1F("h_NDY_event", "N #gamma^{*} from POWHEG", 10, 0, 10);
    TH1F *h_NDY_event_fwd = new TH1F("h_NDY_event_fwd", "N #gamma^{*} from POWHEG in 2.5 < #it{y} < 4", 10, 0, 10);
    TH1F *h_YDY = new TH1F("h_YDY", "h_YDY", 160, -8, 8);
    TH1F *h_YDY_fwd = new TH1F("h_YDY_fwd", "h_YDY_fwd", 150, -4, -2.5);

    MCTree->Draw("N_gamma_gen>>h_NDY_event", "", "goff");

    MCTree->Draw("N_gamma_gen>>h_NDY_event_fwd", "Y_gamma_gen > -4 && Y_gamma_gen<-2.5", "goff");

    MCTree->Draw("Y_gamma_gen>>h_YDY", "", "goff");
    MCTree->Draw("Y_gamma_gen>>h_YDY_fwd", "Y_gamma_gen > -4 && Y_gamma_gen<-2.5", "goff");

    TCanvas *c = new TCanvas("c", "c", 1200, 1220);
    c->Divide(2, 2);
    c->cd(1);
    gPad->SetLogy();
    h_NDY_event->GetYaxis()->SetRangeUser(0.1, 1.2e+4);
    h_NDY_event->Draw();
    c->cd(2);
    gPad->SetLogy();
    h_NDY_event_fwd->GetYaxis()->SetRangeUser(0.1, 1.2e+4);
    h_NDY_event_fwd->Draw();
    c->cd(3);
    h_YDY->Draw();
    c->cd(4);
    h_YDY_fwd->Draw();

    TFile *fOut = new TFile("~/cernbox/HF_dimuons/fit_data/ingredient_cs_powheg/DY_cs.root", "RECREATE");
    fOut->cd();
    h_NDY_event->Write();
    h_NDY_event_fwd->Write();
}