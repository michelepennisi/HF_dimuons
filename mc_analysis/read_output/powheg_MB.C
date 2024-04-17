#include "../../common_include.h"
#include "THStack.h"
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

    TString Hist_POWHEGSIM_fname_Charm_NoCut = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_charm_nocut_Version_5_skimmed_fwd_withHF_Q_local/powheg_charm_nocut_MC_output_Hist_294925.root";
    TString Tree_POWHEGSIM_fname_Charm_NoCut = "/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_charm_nocut_Version_5_skimmed_fwd_withHF_Q_local/powheg_charm_nocut_MC_output_Tree_294925.root";

    Int_t Mass_Binning = 52;
    Int_t Low_Mass = 4;
    Int_t High_Mass = 9;
    Double_t LowM_cut = 8.;
    Double_t HighM_cut = 11.;
    Int_t Pt_Binning = 60;
    Int_t Low_Pt = 0;
    Int_t High_Pt = 10;

    Color_t color[5] = {kMagenta + 2, kSpring - 6, kAzure + 9, kBlack, kRed - 2};
};

void powheg_MB()
{

    opt info;
    TString Name_Dimu[] = {"Charm", "Beauty", "HF_Mixed", "LF", "LF_HF_Mixed"};

    TFile *fIn_POWHEG = new TFile(info.Tree_POWHEGSIM_fname_Beauty, "READ");
    TH1F *h_Nevents_POWHEG = (TH1F *)fIn_POWHEG->Get("h_Nevents");
    h_Nevents_POWHEG->SetName("h_Nevents_POWHEG");

    fIn_POWHEG->ls();

    TIter next_POWHEG(fIn_POWHEG->GetListOfKeys());
    TKey *key_POWHEG = new TKey();
    TObject *Obj_Tree_POWHEG = nullptr;
    TObject *M_Obj_Hist_POWHEG[5]; // Charm, Beauty, HF-Mixed, LF, LF-HF-Mixed
    TCanvas *c_M_POWHEG = new TCanvas("c_M_POWHEG", "c_M_POWHEG", 1200, 1200);

    Int_t i_Hist = 0;
    while ((key_POWHEG = (TKey *)next_POWHEG()))
    {
        if (TString::Format("%s", key_POWHEG->GetClassName()).Contains("TTree") && TString::Format("%s", key_POWHEG->GetName()).EqualTo(Form("DiMuon_Rec_%s", Name_Dimu[i_Hist].Data())))
        {
            Obj_Tree_POWHEG = (TObject *)key_POWHEG->ReadObj();
            Obj_Tree_POWHEG->InheritsFrom(TString::Format("TTree::Class()"));
            ((TTree *)Obj_Tree_POWHEG)->SetDirectory(gROOT);
            Obj_Tree_POWHEG->InheritsFrom(TString::Format("TH1F::Class()"));
            TString Name;
            if (TString::Format("%s", fIn_POWHEG->GetName()).Contains("LHC23i1"))
                Name.Form("POWHEG Charm sim");
            else if (TString::Format("%s", fIn_POWHEG->GetName()).Contains("LHC23i2"))
                Name.Form("POWHEG Beauty sim");
            else if (TString::Format("%s", fIn_POWHEG->GetName()).Contains("LHC22b3"))
                Name.Form("PYTHIA HF-enr.sim");
            else if (TString::Format("%s", fIn_POWHEG->GetName()).Contains("LHC22c1"))
                Name.Form("PYTHIA MB sim");
            // M_Obj_Hist_POWHEG[i_Hist] = new TH1F(Form("h_M_%s", Name_Dimu[i_Hist].Data()), Form("h_M_%s", Name_Dimu[i_Hist].Data()), 130, 4, 30);
            M_Obj_Hist_POWHEG[i_Hist] = new TH1F(Form("h_M_%s", Name_Dimu[i_Hist].Data()), Form("%s from %s", Name_Dimu[i_Hist].Data(), Name.Data()), 30, 0, 30);
            // ((TTree *)Obj_Tree_POWHEG)->Draw(Form("m>>h_M_%s", Name_Dimu[i_Hist].Data()), "((m>4 && m<8)|| (m>11 && m<30)) && pt<30", "goff");
            ((TTree *)Obj_Tree_POWHEG)->Draw(Form("m>>h_M_%s", Name_Dimu[i_Hist].Data()), "", "goff");
            hist1D_graphic_opt((TH1F *)M_Obj_Hist_POWHEG[i_Hist], kFALSE, 1, 20, info.color[i_Hist], 1. / ((TH1F *)M_Obj_Hist_POWHEG[i_Hist])->Integral());
            M_Obj_Hist_POWHEG[i_Hist]->Draw("PESAME");
            i_Hist++;
        }
    }
    gPad->BuildLegend();

    TFile *fIn_PYTHIA = new TFile(info.Tree_POWHEGSIM_fname_Charm, "READ");
    TH1F *h_Nevents_PYTHIA = (TH1F *)fIn_PYTHIA->Get("h_Nevents");
    h_Nevents_PYTHIA->SetName("h_Nevents_PYTHIA");
    fIn_PYTHIA->ls();

    TIter next_PYTHIA(fIn_PYTHIA->GetListOfKeys());
    TKey *key_PYTHIA = new TKey();
    TObject *Obj_Tree_PYTHIA = nullptr;
    TObject *M_Obj_Hist_PYTHIA[5]; // Charm, Beauty, HF-Mixed, LF, LF-HF-Mixed
    TCanvas *c_M_PYTHIA = new TCanvas("c_M_PYTHIA", "c_M_PYTHIA", 1200, 1200);

    i_Hist = 0;
    while ((key_PYTHIA = (TKey *)next_PYTHIA()))
    {
        if (TString::Format("%s", key_PYTHIA->GetClassName()).Contains("TTree") && TString::Format("%s", key_PYTHIA->GetName()).EqualTo(Form("DiMuon_Rec_%s", Name_Dimu[i_Hist].Data())))
        {
            Obj_Tree_PYTHIA = (TObject *)key_PYTHIA->ReadObj();
            Obj_Tree_PYTHIA->InheritsFrom(TString::Format("TTree::Class()"));
            ((TTree *)Obj_Tree_PYTHIA)->SetDirectory(gROOT);
            Obj_Tree_PYTHIA->InheritsFrom(TString::Format("TH1F::Class()"));
            TString Name;
            if (TString::Format("%s", fIn_PYTHIA->GetName()).Contains("LHC23i1"))
                Name.Form("POWHEG Charm sim");
            else if (TString::Format("%s", fIn_PYTHIA->GetName()).Contains("LHC23i2"))
                Name.Form("POWHEG Beauty sim");
            else if (TString::Format("%s", fIn_PYTHIA->GetName()).Contains("LHC22b3"))
                Name.Form("PYTHIA HF-enr.sim");
            else if (TString::Format("%s", fIn_PYTHIA->GetName()).Contains("LHC22c1"))
                Name.Form("PYTHIA MB sim");
            // M_Obj_Hist_PYTHIA[i_Hist] = new TH1F(Form("h_M_%s", Name_Dimu[i_Hist].Data()), Form("h_M_%s", Name_Dimu[i_Hist].Data()), 130, 4, 30);
            M_Obj_Hist_PYTHIA[i_Hist] = new TH1F(Form("h_M_%s", Name_Dimu[i_Hist].Data()), Form("%s from %s", Name_Dimu[i_Hist].Data(), Name.Data()), 30, 0, 30);
            // ((TTree *)Obj_Tree_PYTHIA)->Draw(Form("pt>>h_M_%s", Name_Dimu[i_Hist].Data()), "((m>4 && m<8)|| (m>11 && m<30)) && pt<30", "goff");
            ((TTree *)Obj_Tree_PYTHIA)->Draw(Form("m>>h_M_%s", Name_Dimu[i_Hist].Data()), "", "goff");
            hist1D_graphic_opt((TH1F *)M_Obj_Hist_PYTHIA[i_Hist], kFALSE, 1, 20, info.color[i_Hist], 1. / ((TH1F *)M_Obj_Hist_PYTHIA[i_Hist])->Integral());
            M_Obj_Hist_PYTHIA[i_Hist]->Draw("PESAME");
            i_Hist++;
        }
    }
    gPad->BuildLegend();

    TCanvas *C_M_PYTHIA_POWHEG_ratio[5];

    for (Int_t i_DiMu = 0; i_DiMu < 3; i_DiMu++)
    {
        if (((TH1F *)M_Obj_Hist_PYTHIA[i_DiMu])->GetEntries() == 0)
            continue;
        ((TH1F *)M_Obj_Hist_PYTHIA[i_DiMu])->SetMarkerColor(kRed - 2);
        ((TH1F *)M_Obj_Hist_POWHEG[i_DiMu])->SetMarkerColor(kAzure - 7);
        TH1F *ratio = (TH1F *)M_Obj_Hist_PYTHIA[i_DiMu]->Clone(Form("PYTHIA_POWHEG_ratio_%s", Name_Dimu[i_DiMu].Data()));
        ratio->Divide((TH1F *)M_Obj_Hist_POWHEG[i_DiMu]);

        C_M_PYTHIA_POWHEG_ratio[i_DiMu] = two_histo_ratio((TH1F *)M_Obj_Hist_PYTHIA[i_DiMu], (TH1F *)M_Obj_Hist_POWHEG[i_DiMu], ratio, Form("C_M_PYTHIA_POWHEG_ratio_%s", Name_Dimu[i_DiMu].Data()), "Charm/Beauty", kTRUE, kTRUE, 4, 30);
    }
}

void single_muon()
{
    opt info;
    TString Name_Dimu[] = {"Charm", "Beauty", "LF"};

    TFile *fIn_First = new TFile(info.Hist_POWHEGSIM_fname_Charm, "READ");

    TH1F *h_Nevents_First = (TH1F *)fIn_First->Get("h_Nevents");
    h_Nevents_First->SetName("h_Nevents_First");

    TDirectory *dir_First = (TDirectory *)fIn_First->Get("Muon_Gen/PYTHIA");
    TIter next_First(dir_First->GetListOfKeys());
    dir_First->GetListOfKeys()->Print();
    TKey *key_First = new TKey();
    TObject *Obj_PtPdg_First = nullptr;
    TObject *M_Obj_Hist_First[5]; // Charm, Beauty, HF-Mixed, LF, LF-HF-Mixed
    TCanvas *c_M_First = new TCanvas("c_M_First", "c_M_First", 1200, 1200);

    Int_t i_Hist = 0;
    while ((key_First = (TKey *)next_First()))
    {
        if (TString::Format("%s", key_First->GetClassName()).Contains("TH2F") && TString::Format("%s", key_First->GetName()).EqualTo(Form("h_PtEta_Muon_Gen_%s_PYTHIAOnly", Name_Dimu[i_Hist].Data())))
        {
            Obj_PtPdg_First = (TObject *)key_First->ReadObj();
            Obj_PtPdg_First->InheritsFrom(TString::Format("TH2F::Class()"));
            ((TH2F *)Obj_PtPdg_First)->SetDirectory(gROOT);
            TString Name;
            if (TString::Format("%s", fIn_First->GetName()).Contains("LHC23i1"))
                Name.Form("POWHEG Charm sim");
            else if (TString::Format("%s", fIn_First->GetName()).Contains("powheg_charm"))
                Name.Form("POWHEG Charm sim No Cut fwd");
            else if (TString::Format("%s", fIn_First->GetName()).Contains("LHC23i2"))
                Name.Form("POWHEG Beauty sim");
            else if (TString::Format("%s", fIn_First->GetName()).Contains("powheg_beauty"))
                Name.Form("POWHEG Beauty sim No Cut fwd");
            else if (TString::Format("%s", fIn_First->GetName()).Contains("LHC22b3"))
                Name.Form("First HF-enr.sim");
            else if (TString::Format("%s", fIn_First->GetName()).Contains("LHC22c1"))
                Name.Form("First MB sim");
            // M_Obj_Hist_First[i_Hist] = new TH1F(Form("h_M_%s", Name_Dimu[i_Hist].Data()), Form("h_M_%s", Name_Dimu[i_Hist].Data()), 130, 4, 30);
            M_Obj_Hist_First[i_Hist] = (TH1F *)((TH2F *)Obj_PtPdg_First)->ProjectionY();
            ((TH1F *)M_Obj_Hist_First[i_Hist])->SetTitle(Form("%s from %s", Name_Dimu[i_Hist].Data(), Name.Data()));
            hist1D_graphic_opt((TH1F *)M_Obj_Hist_First[i_Hist], kFALSE, 10, 20, info.color[i_Hist], 1. / ((TH1F *)M_Obj_Hist_First[i_Hist])->Integral());
            M_Obj_Hist_First[i_Hist]->Draw("PESAME");
            i_Hist++;
        }
    }

    TFile *fIn_Second = new TFile(info.Hist_POWHEGSIM_fname_Charm_NoCut, "READ");
    TH1F *h_Nevents_Second = (TH1F *)fIn_Second->Get("h_Nevents");
    h_Nevents_Second->SetName("h_Nevents_Second");

    TDirectory *dir_Second = (TDirectory *)fIn_Second->Get("Muon_Gen/PYTHIA");
    TIter next_Second(dir_Second->GetListOfKeys());
    dir_Second->GetListOfKeys()->Print();
    TKey *key_Second = new TKey();
    TObject *Obj_PtPdg_Second = nullptr;
    TObject *M_Obj_Hist_Second[5]; // Charm, Beauty, HF-Mixed, LF, LF-HF-Mixed
    TCanvas *c_M_Second = new TCanvas("c_M_Second", "c_M_Second", 1200, 1200);

    i_Hist = 0;
    while ((key_Second = (TKey *)next_Second()))
    {
        if (TString::Format("%s", key_Second->GetClassName()).Contains("TH2F") && TString::Format("%s", key_Second->GetName()).EqualTo(Form("h_PtEta_Muon_Gen_%s_PYTHIAOnly", Name_Dimu[i_Hist].Data())))
        {
            Obj_PtPdg_Second = (TObject *)key_Second->ReadObj();
            Obj_PtPdg_Second->InheritsFrom(TString::Format("TH2F::Class()"));
            ((TH2F *)Obj_PtPdg_Second)->SetDirectory(gROOT);
            TString Name;
            if (TString::Format("%s", fIn_Second->GetName()).Contains("LHC23i1"))
                Name.Form("POWHEG Charm sim");
            else if (TString::Format("%s", fIn_First->GetName()).Contains("powheg_charm"))
                Name.Form("POWHEG Charm sim No Cut fwd");
            else if (TString::Format("%s", fIn_Second->GetName()).Contains("LHC23i2"))
                Name.Form("POWHEG Beauty sim");
            else if (TString::Format("%s", fIn_First->GetName()).Contains("powheg_charm"))
                Name.Form("POWHEG Beauty sim No Cut fwd");
            else if (TString::Format("%s", fIn_Second->GetName()).Contains("LHC22b3"))
                Name.Form("PYTHIA HF-enr.sim");
            else if (TString::Format("%s", fIn_Second->GetName()).Contains("LHC22c1"))
                Name.Form("PYTHIA MB sim");
            // M_Obj_Hist_Second[i_Hist] = new TH1F(Form("h_M_%s", Name_Dimu[i_Hist].Data()), Form("h_M_%s", Name_Dimu[i_Hist].Data()), 130, 4, 30);
            M_Obj_Hist_Second[i_Hist] = (TH1F *)((TH2F *)Obj_PtPdg_Second)->ProjectionY();
            ((TH1F *)M_Obj_Hist_Second[i_Hist])->SetTitle(Form("%s from %s", Name_Dimu[i_Hist].Data(), Name.Data()));
            hist1D_graphic_opt((TH1F *)M_Obj_Hist_Second[i_Hist], kFALSE, 10, 20, info.color[i_Hist], 1. / ((TH1F *)M_Obj_Hist_Second[i_Hist])->Integral());
            M_Obj_Hist_Second[i_Hist]->Draw("PESAME");
            i_Hist++;
        }
    }

    TCanvas *C_M_PYTHIA_POWHEG_ratio[2];

    for (Int_t i_DiMu = 0; i_DiMu < 3; i_DiMu++)
    {
        ((TH1F *)M_Obj_Hist_First[i_DiMu])->SetMarkerColor(kRed - 2);
        ((TH1F *)M_Obj_Hist_Second[i_DiMu])->SetMarkerColor(kAzure - 7);
        TH1F *ratio = (TH1F *)M_Obj_Hist_First[i_DiMu]->Clone(Form("PYTHIA_POWHEG_ratio_%s", Name_Dimu[i_DiMu].Data()));
        ratio->Divide((TH1F *)M_Obj_Hist_Second[i_DiMu]);

        C_M_PYTHIA_POWHEG_ratio[i_DiMu] = two_histo_ratio((TH1F *)M_Obj_Hist_First[i_DiMu], (TH1F *)M_Obj_Hist_Second[i_DiMu], ratio, Form("C_M_PYTHIA_POWHEG_ratio_%s", Name_Dimu[i_DiMu].Data()), "Charm/Beauty", kTRUE, kTRUE, 4, 30);
    }
}