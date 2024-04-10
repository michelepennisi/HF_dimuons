#include "/home/michele_pennisi/cernbox/HF_dimuons/common_include.h"

using namespace std;
struct opt
{
    TString Hist_PYTHIAMB_fname = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC22c1/Version_5_AliAOD_skimmed_fwd_fullstat/no_HF_mixed_origin/LHC22c1_MC_output_Hist_merged.root";
    TString Tree_PYTHIAMB_fname = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC22c1/Version_5_AliAOD_skimmed_fwd_fullstat/LHC22c1_MC_output_Tree_merged.root";

    TString Hist_PYTHIAHF_enriched_fname = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC22b3/Version_5_AliAOD_skimmed_fwd_fullstat/LHC22b3_MC_output_Hist_merged.root";
    TString Tree_PYTHIAHF_enriched_fname = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC22b3/Version_5_AliAOD_skimmed_fwd_fullstat/LHC22b3_MC_output_Tree_merged.root";

    TString Hist_POWHEGSIM_fname_Charm = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i1/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i1_MC_output_Hist_merged.root";
    TString Tree_POWHEGSIM_fname_Charm = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i1/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i1_MC_output_Tree_merged.root";

    TString Hist_POWHEGSIM_fname_Beauty = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i2_MC_output_Hist_merged.root";
    TString Tree_POWHEGSIM_fname_Beauty = "/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version_5_AliAOD_skimmed_fwd_fullstat/LHC23i2_MC_output_Tree_merged.root";

    TString Generator = "powheg";

    Bool_t Norm_from_fit = kFALSE;

    Int_t Mass_Binning = 52;
    Int_t Low_Mass = 4;
    Int_t High_Mass = 9;
    Double_t LowM_cut = 8.;
    Double_t HighM_cut = 11.;
    Int_t Pt_Binning = 60;
    Int_t Low_Pt = 0;
    Int_t High_Pt = 10;
};

void yield_check()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    opt info;
    // MC stuff;
    TFile *fIn_Charm;
    TFile *fIn_Beauty;
    TH1F *h_Nevents_Charm;
    TH1F *h_Nevents_Beauty;
    Double_t Nev_MC_Charm;
    Double_t Nev_MC_Beauty;

    // MC data;
    TFile *fIn_data = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/data_analysis/Tree_MassPt_MassCut4_Run2.root", "READ");

    // Generator condition
    TString tree_name;
    if (info.Generator.Contains("pythia"))
    {
        fIn_Charm = new TFile(info.Tree_PYTHIAHF_enriched_fname, "READ");
        fIn_Charm->ls();
        fIn_Beauty = new TFile(info.Tree_PYTHIAHF_enriched_fname, "READ");
        h_Nevents_Charm = (TH1F *)fIn_Charm->Get("h_Nevents");
        Nev_MC_Charm = 216 * h_Nevents_Charm->GetBinContent(2);  // Number of MB event in the MC
        Nev_MC_Beauty = 216 * h_Nevents_Charm->GetBinContent(2); // Number of MB event in the MC
        cout << "Nev_MC for PYTHIA HF-enriched: " << Nev_MC_Charm << endl;
        tree_name.Form("DiMuon_Rec");
    }
    else if (info.Generator.Contains("powheg"))
    {
        fIn_Charm = new TFile(info.Tree_POWHEGSIM_fname_Charm, "READ");
        fIn_Charm->ls();
        h_Nevents_Charm = (TH1F *)fIn_Charm->Get("h_Nevents");
        h_Nevents_Charm->SetName("h_Nevents_Charm");
        cout<<h_Nevents_Charm->GetBinContent(2)<<endl;
        fIn_Beauty = new TFile(info.Tree_POWHEGSIM_fname_Beauty, "READ");
        h_Nevents_Beauty = (TH1F *)fIn_Beauty->Get("h_Nevents");
        h_Nevents_Beauty->SetName("h_Nevents_Beauty");
        tree_name.Form("DiMuon_Rec_PowhegOnly");
        cout<<h_Nevents_Beauty->GetBinContent(2)<<endl;
        Nev_MC_Charm = (30.5 * (56.42 / 5.0)) * h_Nevents_Charm->GetBinContent(2);
        Nev_MC_Beauty = (15.0 * (56.42 / 0.5)) * h_Nevents_Beauty->GetBinContent(2);
        cout << "Nev_MC for POWHEG sim for CHARM: " << Nev_MC_Charm << endl;
        cout << "Nev_MC for POWHEG sim for Beauty: " << Nev_MC_Beauty << endl;
    }

    TH1F *h_M_Charm;
    TH1F *h_M_Beauty;
    TH1F *h_M_HF_Mixed;
    TH1F *h_M_Sum;
    TString var[2] = {"m", "pt"};
    TString label_hist[2] = {"#it{m}_{#mu#mu} (GeV/#it{c}^{2})", "#it{p}_{T} (GeV/#it{c})"};
    for (Int_t I_Var = 0; I_Var < 2; I_Var++)
    {

        Double_t low_var;
        Double_t high_var;
        Double_t binning_var;
        if (I_Var == 0)
        {
            low_var = info.Low_Mass;
            high_var = info.High_Mass;
        }
        else
        {
            low_var = info.Low_Pt;
            high_var = info.High_Pt;
        }
        gROOT->cd();
        TTree *tree_data = (TTree *)fIn_data->Get("rec_data_tree");
        TH1F *h_M_Data = new TH1F(Form("h_%s_Data", var[I_Var].Data()), Form("Data ; %s", label_hist[I_Var].Data()), (high_var - low_var) / 0.2, low_var, high_var);
        TH1D *fhNEv = (TH1D *)fIn_data->Get("fhNEv");
        Double_t NEv_MB_Data = (2384.73 * fhNEv->GetBinContent(3));
        cout << "N events for Data: " << NEv_MB_Data << endl;
        tree_data->Draw(Form("%s>>h_%s_Data", var[I_Var].Data(), var[I_Var].Data()), Form("(m>%d && m<%d) && pt<%d", info.Low_Mass, info.High_Mass, info.High_Pt), "goff");
        hist1D_graphic_opt(h_M_Data, kFALSE, 1, 20, kBlack, 1. / NEv_MB_Data);

        h_M_Charm = new TH1F(Form("h_%s_Charm", var[I_Var].Data()), Form("#mu#mu #leftarrow c ; %s", label_hist[I_Var].Data()), (high_var - low_var) / 0.2, low_var, high_var);
        h_M_Beauty = new TH1F(Form("h_%s_Beauty", var[I_Var].Data()), Form("#mu#mu #leftarrow b; %s", label_hist[I_Var].Data()), (high_var - low_var) / 0.2, low_var, high_var);
        h_M_HF_Mixed = new TH1F(Form("h_%s_HF_Mixed", var[I_Var].Data()), Form("#mu#mu #leftarrow c,b; %s", label_hist[I_Var].Data()), (high_var - low_var) / 0.2, low_var, high_var);

        TTree *DiMuon_Rec_Charm = (TTree *)fIn_Charm->Get(Form("%s_Charm", tree_name.Data()));
        DiMuon_Rec_Charm->Draw(Form("%s>>h_%s_Charm", var[I_Var].Data(), var[I_Var].Data()), Form("(m>%d && m<%d) && pt<%d", info.Low_Mass, info.High_Mass, info.High_Pt), "goff");

        TTree *DiMuon_Rec_Beauty = (TTree *)fIn_Beauty->Get(Form("%s_Beauty", tree_name.Data()));
        DiMuon_Rec_Beauty->Draw(Form("%s>>h_%s_Beauty", var[I_Var].Data(), var[I_Var].Data()), Form("(m>%d && m<%d) && pt<%d", info.Low_Mass, info.High_Mass, info.High_Pt), "goff");

        TTree *DiMuon_Rec_HF_Mixed = (TTree *)fIn_Beauty->Get(Form("%s_HF_Mixed", tree_name.Data()));
        DiMuon_Rec_HF_Mixed->Draw(Form("%s>>h_%s_HF_Mixed", var[I_Var].Data(), var[I_Var].Data()), Form("(m>%d && m<%d) && pt<%d", info.Low_Mass, info.High_Mass, info.High_Pt), "goff");

        if (info.Norm_from_fit && info.Generator.Contains("pythia"))
        {
            hist1D_graphic_opt(h_M_Charm, kFALSE, 1, 20, kMagenta + 2, 1.77 / Nev_MC_Charm);
            hist1D_graphic_opt(h_M_Beauty, kFALSE, 1, 20, kSpring - 6, 0.34 / Nev_MC_Beauty);

            hist1D_graphic_opt(h_M_HF_Mixed, kFALSE, 1, 20, kAzure + 9, (h_M_Data->Integral() * 0.036)/h_M_HF_Mixed->GetEntries());
        }
        else if (info.Norm_from_fit && info.Generator.Contains("powheg"))
        {
            hist1D_graphic_opt(h_M_Charm, kFALSE, 1, 20, kMagenta + 2, 6.54 / Nev_MC_Charm);
            hist1D_graphic_opt(h_M_Beauty, kFALSE, 1, 20, kSpring - 6, 0.86 / Nev_MC_Beauty);
            hist1D_graphic_opt(h_M_HF_Mixed, kFALSE, 1, 20, kAzure + 9,  (h_M_Data->Integral() * 0.036)/h_M_HF_Mixed->GetEntries());
        }
        else
        {
            hist1D_graphic_opt(h_M_Charm, kFALSE, 1, 20, kMagenta + 2, 1. / Nev_MC_Charm);
            hist1D_graphic_opt(h_M_Beauty, kFALSE, 1, 20, kSpring - 6, 1. / Nev_MC_Beauty);
            hist1D_graphic_opt(h_M_HF_Mixed, kFALSE, 1, 20, kAzure + 9, 1. / Nev_MC_Beauty);
        }

        h_M_Sum = (TH1F *)h_M_Charm->Clone(Form("h_%s_Sum", var[I_Var].Data()));
        h_M_Sum->SetTitle("#mu#mu #leftarrow HF");
        h_M_Sum->Add(h_M_Beauty);
        h_M_Sum->Add(h_M_HF_Mixed);

        hist1D_graphic_opt(h_M_Sum, kFALSE, 1, 20, kOrange + 7, 1.);

        cout << "Charm fraction: " << h_M_Charm->Integral() / h_M_Sum->Integral() << endl;
        cout << "Beauty fraction: " << h_M_Beauty->Integral() / h_M_Sum->Integral() << endl;
        cout << "HF_Mixed fraction: " << h_M_HF_Mixed->Integral() / h_M_Sum->Integral() << endl;

        TH1F *ratio_Data_Sum = (TH1F *)h_M_Data->Clone(Form("h_%s_ratio_Data_Sum", var[I_Var].Data()));
        ratio_Data_Sum->Divide(h_M_Sum);
        ratio_Data_Sum->GetYaxis()->SetTitle("Data/Cocktail");
        ratio_Data_Sum->GetYaxis()->SetRangeUser(0.65, 1.65);
        TCanvas *canvas = two_histo_ratio(h_M_Data, h_M_Sum, ratio_Data_Sum, Form("%s_yields_check", var[I_Var].Data()), " ", kTRUE, kTRUE, low_var, high_var);
        TPad *pad1 = (TPad *)canvas->GetListOfPrimitives()->FindObject("pad1");
        pad1->SetName(Form("%s_pad1", var[I_Var].Data()));
        TLegend *legend = (TLegend *)pad1->FindObject("Legend");
        legend->SetName(Form("%s_legend", var[I_Var].Data()));
        TString Header;
        if (info.Norm_from_fit)
            Header.Form("with fraction from fit");
        else
            Header.Form("PYTHIA fraction");

        Legend_settings(legend, legend->GetX1() + 0.05, legend->GetX2() + 0.05, legend->GetY1() + 0.075, legend->GetY2() + 0.175, Header.Data());
        legend->AddEntry(h_M_Charm, Form("%s, fr: %0.2f", h_M_Charm->GetTitle(), h_M_Charm->Integral() / h_M_Sum->Integral() * 100));
        legend->AddEntry(h_M_Beauty, Form("%s, fr: %0.2f", h_M_Beauty->GetTitle(), h_M_Beauty->Integral() / h_M_Sum->Integral() * 100));
        legend->AddEntry(h_M_HF_Mixed, Form("%s, fr: %0.2f", h_M_HF_Mixed->GetTitle(), h_M_HF_Mixed->Integral() / h_M_Sum->Integral() * 100));
        pad1->Modified();
        pad1->Update();

        if (info.Norm_from_fit)
            canvas->SaveAs(Form("images/%s_%s_fr_fit.pdf", canvas->GetName(), info.Generator.Data()));
        else
            canvas->SaveAs(Form("images/%s_%s_original.pdf", canvas->GetName(), info.Generator.Data()));
    }

    // kMagenta + 2, kSpring - 6, kOrange + 7, kAzure + 9

    // hist1D_graphic_opt(h_M_Sum, kFALSE, 1, 20, kOrange + 7, 0.3 / Nev_MC);
}