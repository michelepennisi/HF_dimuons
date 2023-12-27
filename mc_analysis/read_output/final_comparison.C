#include "/home/michele_pennisi/cernbox/common_include.h"

TFile *fIn_PYTHIA_CENTRALIZED_MB = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC22c1_Version_5_AliAOD/Merged_LHC22c1_MC_output_Hist.root", "READ"); // saving only muons and dimuons
// TFile *fIn_PYTHIA_CENTRALIZED_MB = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC22c1/294925/output/Merged_LHC22c1_MC_output_Hist_294925.root", "READ"); // saving also Hadrons
TH1F *h_Nevents_PYTHIA_MB_CENTRALIZED = (TH1F *)fIn_PYTHIA_CENTRALIZED_MB->Get("h_Nevents");
Int_t Nev_PYTHIA_MB_CENTRALIZED = (Int_t)h_Nevents_PYTHIA_MB_CENTRALIZED->GetBinContent(2);

Double_t Norm_PYTHIA_MB_CENTRALIZED = (78.05 / ((Double_t)Nev_PYTHIA_MB_CENTRALIZED)); // in mb

// TFile *fIn_Pythia_STANDALONE_NoDiffr = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_Def_MC_output_Hist_100000.root", "READ"); // pi and Kaons decay with NonDiffractive SowfQCD
// TH1F *h_Nevents_PYTHIA_STANDALONE_NoDiffr = (TH1F *)fIn_Pythia_STANDALONE_NoDiffr->Get("h_Nevents");
// Int_t Nev_PYTHIA_STANDALONE_NoDiffr = (Int_t)h_Nevents_PYTHIA_STANDALONE_NoDiffr->GetBinContent(2);
// Double_t Norm_PYTHIA_STANDALONE_NoDiffr = (56.42 / ((Double_t)Nev_PYTHIA_STANDALONE_NoDiffr)); // in mb

TFile *fIn_Pythia_STANDALONE_INEL = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_inel_LFoff_Def_pythia_sim_74_DefaultBR_output_Hist_100000.root", "READ"); // pi decay with NonDiffractive SoftQCD
TH1F *h_Nevents_PYTHIA_STANDALONE_INEL = (TH1F *)fIn_Pythia_STANDALONE_INEL->Get("h_Nevents");
Int_t Nev_PYTHIA_STANDALONE_INEL = (Int_t)h_Nevents_PYTHIA_STANDALONE_INEL->GetBinContent(2);
Double_t Norm_PYTHIA_STANDALONE_INEL = (78.05 / ((Double_t)Nev_PYTHIA_STANDALONE_INEL)); // in mb

TFile *fIn_Powheg_Charm = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC23i1_Version_5_AliAOD/Merged_LHC23i1_MC_output_Hist.root", "READ"); // Only muons and dimuons
// TFile *fIn_Powheg_Charm = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC23i1_Version_5_AliAOD/LHC23i1_MC_output_Hist_294009.root", "READ"); // small stat with hadrons
TH1F *h_Nevents_Powheg_Charm = (TH1F *)fIn_Powheg_Charm->Get("h_Nevents");
Int_t Nev_Powheg_Charm = (Int_t)h_Nevents_Powheg_Charm->GetBinContent(2);
Double_t Norm_Powheg_Charm = (56.42 / ((Double_t)Nev_Powheg_Charm)); // in mb

void Powheg_Gen_muons(TString gen);
void Gen_muons();
void LF_Hadron(TString gen);
void Rec_Muons();
void Rec_Dimu(TString gen);
void Pythia6_Pythia8_LF_Hadron();
void Pythia6_Pythia8_LF_Muons(Bool_t cross_section);

void Pythia6_LF_DiMuons();
void Pythia6_Pythia8_LF_DiMuons(Bool_t cross_section, TString DiMu_origin);
void Pythia6_Pythia8_LF_Hadron();
void final_comparison()
{
    // Powheg_Gen_muons("Gen");
    // Powheg_Gen_muons("Rec");

    // Gen_muons();
    // Rec_Muons();
    // Rec_Dimu("Powheg");
    // LF_Hadron("Powheg");
    // Pythia6_Pythia8_LF_Muons(kTRUE);
    // Pythia6_LF_DiMuons();
    Pythia6_Pythia8_LF_DiMuons(kFALSE, "DiMuon_Rec/LF_HF_Mixed_origin");
    // Pythia6_Pythia8_LF_Hadron();
}

void Pythia6_LF_DiMuons()
{
    const Int_t n_bin_pt = 14;
    Double_t pt_array_low_bin[n_bin_pt] = {0.0, 0.5, 1., 1.5, 2.5, 3.5, 5.0, 7.5, 10., 12.5, 15.0, 20.0, 26.0, 30.0};

    const Int_t NCut = 2;
    TString Name_Histogram[NCut] = {"DiMuon_Rec/LF_Origin/h_PtM_DiMuon_Rec_fromLF", "DiMuon_Rec/LF_Origin/h_PtY_DiMuon_Rec_fromLF"};
    TString Title_Canvas[NCut] = {"Reconstructed #mu^{#plus}#mu^{#minus} #leftarrow LF hadrons (no POWHEG particles)", "Reconstructed #mu^{#plus}#mu^{#minus} #leftarrow LF hadrons (no POWHEG particles)"};
    TString Name_Canvas[NCut] = {"entries_dimuonRec_pythia6_8", "entries_dimuonRec_pythia6_8"};
    TLine *p_line = new TLine(0, 1, 10, 1);
    p_line->SetLineColor(kPink - 8);
    p_line->SetLineStyle(9);
    p_line->SetLineWidth(3);
    TH2F *h_PtM_DiMuon_LF_sum_PYTHIA6[6];
    TH2F *h_PtY_DiMuon_LF_sum_PYTHIA6[6];

    TDirectory *dir = fIn_Powheg_Charm->GetDirectory("DiMuon_Rec/LF_origin");
    TIter next(dir->GetListOfKeys());
    TKey *key = new TKey();
    Int_t counter = 0;
    while ((key = (TKey *)next()))
    {
        if (TString::Format("%s", key->GetClassName()).Contains("TH2F"))
        {
            TH2F *hist_PYTHIA6 = (TH2F *)fIn_Powheg_Charm->Get(Form("%s/%s", dir->GetPath(), key->GetName()));
            printf(" key : %s is a %s in %s\n", key->GetName(), key->GetClassName(), dir->GetPath());

            if (TString::Format("%s", key->GetName()).Contains("PtM"))
                h_PtM_DiMuon_LF_sum_PYTHIA6[counter] = (TH2F *)fIn_Powheg_Charm->Get(Form("%s/%s", dir->GetPath(), key->GetName()));

            else if (TString::Format("%s", key->GetName()).Contains("PtY"))
            {
                h_PtY_DiMuon_LF_sum_PYTHIA6[counter] = (TH2F *)fIn_Powheg_Charm->Get(Form("%s/%s", dir->GetPath(), key->GetName()));
                counter++;
            }
        }
    }

    for (Int_t i_gen = 0; i_gen < 6; i_gen++)
    {
        TH1F *h_Pt_LF_DiMuon_Rec_PYTHIA6 = (TH1F *)h_PtM_DiMuon_LF_sum_PYTHIA6[i_gen]->ProjectionX();
        h_Pt_LF_DiMuon_Rec_PYTHIA6->SetName(h_PtM_DiMuon_LF_sum_PYTHIA6[i_gen]->GetName());
        h_Pt_LF_DiMuon_Rec_PYTHIA6->SetTitle(h_PtM_DiMuon_LF_sum_PYTHIA6[i_gen]->GetName());
        TH1F *h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED = (TH1F *)h_Pt_LF_DiMuon_Rec_PYTHIA6->Rebin(n_bin_pt - 1, "h_Pt_LF_DiMuon_Rec_PYTHIA6_rebinned", pt_array_low_bin);
        hist1D_graphic_opt(h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED, kTRUE, 0, 20, kAzure + 2, 1. / h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED->Integral());

        if (i_gen == 0)
            h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED->DrawCopy("PLC PMC");
        else
            h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED->DrawCopy("PLC PMC SAME");
    }
}

void Pythia6_Pythia8_LF_DiMuons(Bool_t cross_section, TString DiMu_origin)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    Int_t Rebin_Eta = 5;
    Int_t Rebin_Pt = 15;

    const Int_t n_bin_pt = 11;
    Double_t pt_array_low_bin[n_bin_pt] = {0.0, 0.5, 1., 1.5, 2.5, 5, 7.5, 10, 12.5, 20., 30.0};

    const Int_t NCut = 2;
    TString Name_Histogram[NCut] = {"DiMuon_Rec/LF_Origin/h_PtM_DiMuon_Rec_fromLF", "DiMuon_Rec/LF_Origin/h_PtY_DiMuon_Rec_fromLF"};

    TString Title_Canvas;
    TString Name_Canvas;
    Double_t correction_factor_PYTHIA6 = 999;

    if (DiMu_origin.Contains("LF_HF"))
    {
        Title_Canvas.Form("Reconstructed #mu^{#plus}#mu^{#minus} #leftarrow LF,HF (no POWHEG particles)");
        correction_factor_PYTHIA6 = 30.5 / 1.96;
        if (cross_section)
            Name_Canvas.Form("LF_HF_mixed_cross_section");
        else
            Name_Canvas.Form("LF_HF_mixed_shape");
    }
    else if (DiMu_origin.Contains("LF"))
    {
        Title_Canvas.Form("Reconstructed #mu^{#plus}#mu^{#minus} #leftarrow LF (no POWHEG particles)");
        correction_factor_PYTHIA6 = 1.;
        if (cross_section)
            Name_Canvas.Form("LF_cross_section");
        else
            Name_Canvas.Form("LF_shape");
    }

    TLine *p_line = new TLine(0, 1, 10, 1);
    p_line->SetLineColor(kPink - 8);
    p_line->SetLineStyle(9);
    p_line->SetLineWidth(3);
    TH2F *h_PtM_DiMuon_LF_sum_PYTHIA6 = nullptr;
    TH2F *h_PtY_DiMuon_LF_sum_PYTHIA6 = nullptr;

    TH2F *h_PtM_DiMuon_LF_sum_PYTHIA8 = nullptr;
    TH2F *h_PtY_DiMuon_LF_sum_PYTHIA8 = nullptr;

    TDirectory *dir = fIn_Powheg_Charm->GetDirectory(Form("%s", DiMu_origin.Data()));
    TIter next(dir->GetListOfKeys());
    TKey *key = new TKey();
    Int_t counter_PYTHIA6 = -1;
    Int_t counter_PYTHIA8 = -1;
    while ((key = (TKey *)next()))
    {
        if (TString::Format("%s", key->GetClassName()).Contains("TH2F") && !TString::Format("%s", key->GetName()).Contains("Powheg"))
        {
            printf(" key : %s is a %s in %s\n", key->GetName(), key->GetClassName(), dir->GetPath());
            TH2F *hist_PYTHIA6 = (TH2F *)fIn_Powheg_Charm->Get(Form("%s/%s", dir->GetPath(), key->GetName()));
            TH2F *hist_PYTHIA8 = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get(Form("%s/%s", DiMu_origin.Data(), key->GetName()));

            if (TString::Format("%s", key->GetName()).Contains("PtM"))
            {
                if (h_PtM_DiMuon_LF_sum_PYTHIA6 == nullptr)
                    h_PtM_DiMuon_LF_sum_PYTHIA6 = (TH2F *)hist_PYTHIA6->Clone("h_PtM_DiMuon_LF_sum_PYTHIA6");
                else
                    h_PtM_DiMuon_LF_sum_PYTHIA6->Add(hist_PYTHIA6);
                if (hist_PYTHIA8 == nullptr)
                    continue;
                if (h_PtM_DiMuon_LF_sum_PYTHIA8 == nullptr)
                    h_PtM_DiMuon_LF_sum_PYTHIA8 = (TH2F *)hist_PYTHIA8->Clone("h_PtM_DiMuon_LF_sum_PYTHIA8");
                else
                    h_PtM_DiMuon_LF_sum_PYTHIA8->Add(hist_PYTHIA8);
            }
            else if (TString::Format("%s", key->GetName()).Contains("PtY"))
            {
                if (h_PtY_DiMuon_LF_sum_PYTHIA6 == nullptr)
                    h_PtY_DiMuon_LF_sum_PYTHIA6 = (TH2F *)hist_PYTHIA6->Clone("h_PtY_DiMuon_LF_sum_PYTHIA6");
                else
                    h_PtY_DiMuon_LF_sum_PYTHIA6->Add(hist_PYTHIA6);

                if (hist_PYTHIA8 == nullptr)
                    continue;

                if (h_PtY_DiMuon_LF_sum_PYTHIA8 == nullptr)
                    h_PtY_DiMuon_LF_sum_PYTHIA8 = (TH2F *)hist_PYTHIA8->Clone("h_PtY_DiMuon_LF_sum_PYTHIA8");
                else
                    h_PtY_DiMuon_LF_sum_PYTHIA8->Add(hist_PYTHIA8);
            }
        }
    }

    TH1F *h_Pt_LF_DiMuon_Rec_PYTHIA6 = (TH1F *)h_PtM_DiMuon_LF_sum_PYTHIA6->ProjectionX();
    h_Pt_LF_DiMuon_Rec_PYTHIA6->SetName("h_Pt_LF_DiMuon_Rec_PYTHIA6");
    h_Pt_LF_DiMuon_Rec_PYTHIA6->SetTitle("from POWHEG+PYTHIA6 sim");
    TH1F *h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED = (TH1F *)h_Pt_LF_DiMuon_Rec_PYTHIA6->Rebin(n_bin_pt - 1, "h_Pt_LF_DiMuon_Rec_PYTHIA6_rebinned", pt_array_low_bin);

    TH1F *h_Pt_LF_DiMuon_Rec_PYTHIA8 = (TH1F *)h_PtM_DiMuon_LF_sum_PYTHIA8->ProjectionX();
    h_Pt_LF_DiMuon_Rec_PYTHIA8->SetName("h_Pt_LF_DiMuon_Rec_PYTHIA8");
    h_Pt_LF_DiMuon_Rec_PYTHIA8->SetTitle("from PYTHIA8 MB sim");
    TH1F *h_Pt_LF_DiMuon_Rec_PYTHIA8_REBINNED = (TH1F *)h_Pt_LF_DiMuon_Rec_PYTHIA8->Rebin(n_bin_pt - 1, "h_Pt_LF_DiMuon_Rec_PYTHIA8_rebinned", pt_array_low_bin);

    if (cross_section)
    {
        h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED->GetYaxis()->SetTitle("d#sigma/d#it{p}_{T} (mb/(GeV/c))");
        hist1D_graphic_opt(h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED, kTRUE, 0, 20, kMagenta + 2, 46.4 / (correction_factor_PYTHIA6 * Nev_Powheg_Charm));

        h_Pt_LF_DiMuon_Rec_PYTHIA8_REBINNED->GetYaxis()->SetTitle("d#sigma/d#it{p}_{T} (mb/(GeV/c))");
        hist1D_graphic_opt(h_Pt_LF_DiMuon_Rec_PYTHIA8_REBINNED, kTRUE, 0, 24, kMagenta + 2, 78.5 / (1. * Nev_PYTHIA_MB_CENTRALIZED));
    }
    else
    {
        h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED->GetYaxis()->SetTitle("shape");
        hist1D_graphic_opt(h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED, kTRUE, 0, 20, kMagenta + 2, 1. / (h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED->Integral()));

        h_Pt_LF_DiMuon_Rec_PYTHIA8_REBINNED->GetYaxis()->SetTitle("shape");
        hist1D_graphic_opt(h_Pt_LF_DiMuon_Rec_PYTHIA8_REBINNED, kTRUE, 0, 24, kMagenta + 2, 1. / (h_Pt_LF_DiMuon_Rec_PYTHIA8->Integral()));
    }

    h_Pt_LF_DiMuon_Rec_PYTHIA8_REBINNED->GetXaxis()->SetRangeUser(0, 20);

    TH1F *Pt_Ratio_PYTHIA8_PYTHIA6 = (TH1F *)h_Pt_LF_DiMuon_Rec_PYTHIA8_REBINNED->Clone("Pt_Ratio_PYTHIA8_PYTHIA6");
    Pt_Ratio_PYTHIA8_PYTHIA6->Divide(h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED);
    Pt_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetTitle("PYTHIA8/PYTHIA6");
    Pt_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetRangeUser(-0.3, 4.3);

    TCanvas *C_Pt_LF_Muons = two_histo_ratio(h_Pt_LF_DiMuon_Rec_PYTHIA8_REBINNED, h_Pt_LF_DiMuon_Rec_PYTHIA6_REBINNED, Pt_Ratio_PYTHIA8_PYTHIA6, Form("c_pt_%s", Name_Canvas.Data()), Title_Canvas.Data(), kTRUE, kTRUE);
    gPad = (TPad *)C_Pt_LF_Muons->GetListOfPrimitives()->FindObject("pad2");
    gPad->cd();
    p_line->SetX2(20);
    p_line->DrawClone("SAME");

    TH1F *h_M_LF_DiMuon_Rec_PYTHIA6 = (TH1F *)h_PtM_DiMuon_LF_sum_PYTHIA6->ProjectionY();
    h_M_LF_DiMuon_Rec_PYTHIA6->SetName("h_M_LF_DiMuon_Rec_PYTHIA6");
    h_M_LF_DiMuon_Rec_PYTHIA6->SetTitle("from POWHEG+PYTHIA6 sim");
    TH1F *h_M_LF_DiMuon_Rec_PYTHIA6_REBINNED = (TH1F *)h_M_LF_DiMuon_Rec_PYTHIA6->Rebin(n_bin_pt - 1, "h_M_LF_DiMuon_Rec_PYTHIA6_rebinned", pt_array_low_bin);

    TH1F *h_M_LF_DiMuon_Rec_PYTHIA8 = (TH1F *)h_PtM_DiMuon_LF_sum_PYTHIA8->ProjectionY();
    h_M_LF_DiMuon_Rec_PYTHIA8->SetName("h_M_LF_DiMuon_Rec_PYTHIA8");
    h_M_LF_DiMuon_Rec_PYTHIA8->SetTitle("from PYTHIA8 MB sim");
    TH1F *h_M_LF_DiMuon_Rec_PYTHIA8_REBINNED = (TH1F *)h_M_LF_DiMuon_Rec_PYTHIA8->Rebin(n_bin_pt - 1, "h_M_LF_DiMuon_Rec_PYTHIA8_rebinned", pt_array_low_bin);

    if (cross_section)
    {
        h_M_LF_DiMuon_Rec_PYTHIA6_REBINNED->GetYaxis()->SetTitle("d#sigma/d#it{m}_{#mu#mu} (mb/(GeV/#it{c}^2))");
        hist1D_graphic_opt(h_M_LF_DiMuon_Rec_PYTHIA6_REBINNED, kTRUE, 0, 20, kMagenta + 2, 46.4 / (correction_factor_PYTHIA6 * Nev_Powheg_Charm));

        h_M_LF_DiMuon_Rec_PYTHIA8_REBINNED->GetYaxis()->SetTitle("d#sigma/d#it{m}_{#mu#mu} (mb/(GeV/#it{c}^2))");
        hist1D_graphic_opt(h_M_LF_DiMuon_Rec_PYTHIA8_REBINNED, kTRUE, 0, 24, kMagenta + 2, 78.5 / (1. * Nev_PYTHIA_MB_CENTRALIZED));
    }
    else
    {
        h_M_LF_DiMuon_Rec_PYTHIA6_REBINNED->GetYaxis()->SetTitle("shape");
        hist1D_graphic_opt(h_M_LF_DiMuon_Rec_PYTHIA6_REBINNED, kTRUE, 0, 20, kMagenta + 2, 1. / (h_M_LF_DiMuon_Rec_PYTHIA6_REBINNED->GetEntries()));

        h_M_LF_DiMuon_Rec_PYTHIA8_REBINNED->GetYaxis()->SetTitle("shape");
        hist1D_graphic_opt(h_M_LF_DiMuon_Rec_PYTHIA8_REBINNED, kTRUE, 0, 24, kMagenta + 2, 1. / (h_M_LF_DiMuon_Rec_PYTHIA8_REBINNED->GetEntries()));
    }

    h_M_LF_DiMuon_Rec_PYTHIA8_REBINNED->GetXaxis()->SetRangeUser(0, 12.5);
    TH1F *M_Ratio_PYTHIA8_PYTHIA6 = (TH1F *)h_M_LF_DiMuon_Rec_PYTHIA8_REBINNED->Clone("M_Ratio_PYTHIA8_PYTHIA6");
    M_Ratio_PYTHIA8_PYTHIA6->Divide(h_M_LF_DiMuon_Rec_PYTHIA6_REBINNED);
    M_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetTitle("PYTHIA8/PYTHIA6");
    M_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetRangeUser(-0.3, 4.3);

    TCanvas *C_M_LF_Muons = two_histo_ratio(h_M_LF_DiMuon_Rec_PYTHIA8_REBINNED, h_M_LF_DiMuon_Rec_PYTHIA6_REBINNED, M_Ratio_PYTHIA8_PYTHIA6, Form("c_M_%s", Name_Canvas.Data()), Title_Canvas.Data(), kTRUE, kTRUE);
    gPad = (TPad *)C_M_LF_Muons->GetListOfPrimitives()->FindObject("pad2");
    gPad->cd();
    p_line->SetX2(12.5);
    p_line->DrawClone("SAME");

    TH1F *h_Y_LF_DiMuon_Rec_PYTHIA6 = (TH1F *)h_PtY_DiMuon_LF_sum_PYTHIA6->ProjectionY();
    h_Y_LF_DiMuon_Rec_PYTHIA6->SetName("h_Y_LF_DiMuon_Rec_PYTHIA6");
    h_Y_LF_DiMuon_Rec_PYTHIA6->SetTitle("from POWHEG+PYTHIA6 sim");

    TH1F *h_Y_LF_DiMuon_Rec_PYTHIA8 = (TH1F *)h_PtY_DiMuon_LF_sum_PYTHIA8->ProjectionY();
    h_Y_LF_DiMuon_Rec_PYTHIA8->SetName("h_Y_LF_DiMuon_Rec_PYTHIA8");
    h_Y_LF_DiMuon_Rec_PYTHIA8->SetTitle("from PYTHIA8 MB sim");

    if (cross_section)
    {
        h_Y_LF_DiMuon_Rec_PYTHIA6->GetYaxis()->SetTitle("d#sigma/d#it{y} (mb)");
        hist1D_graphic_opt(h_Y_LF_DiMuon_Rec_PYTHIA6, kTRUE, Rebin_Pt, 20, kMagenta + 2, 46.4 / (correction_factor_PYTHIA6 * Nev_Powheg_Charm));

        h_Y_LF_DiMuon_Rec_PYTHIA8->GetYaxis()->SetTitle("d#sigma/d#it{y} (mb)");
        hist1D_graphic_opt(h_Y_LF_DiMuon_Rec_PYTHIA8, kTRUE, Rebin_Pt, 24, kMagenta + 2, 78.5 / (1. * Nev_PYTHIA_MB_CENTRALIZED));
    }
    else
    {
        h_Y_LF_DiMuon_Rec_PYTHIA6->GetYaxis()->SetTitle("shape");
        hist1D_graphic_opt(h_Y_LF_DiMuon_Rec_PYTHIA6, kFALSE, Rebin_Pt, 20, kMagenta + 2, 1. / (h_Y_LF_DiMuon_Rec_PYTHIA6->GetEntries()));

        h_Y_LF_DiMuon_Rec_PYTHIA8->GetYaxis()->SetTitle("shape");
        hist1D_graphic_opt(h_Y_LF_DiMuon_Rec_PYTHIA8, kFALSE, Rebin_Pt, 24, kMagenta + 2, 1. / (h_Y_LF_DiMuon_Rec_PYTHIA8->GetEntries()));
    }

    h_Y_LF_DiMuon_Rec_PYTHIA8->GetXaxis()->SetRangeUser(-4.0, -2.5);
    TH1F *Y_Ratio_PYTHIA8_PYTHIA6 = (TH1F *)h_Y_LF_DiMuon_Rec_PYTHIA8->Clone("Y_Ratio_PYTHIA8_PYTHIA6");
    Y_Ratio_PYTHIA8_PYTHIA6->Divide(h_Y_LF_DiMuon_Rec_PYTHIA6);
    Y_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetTitle("PYTHIA8/PYTHIA6");
    Y_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetRangeUser(-0.3, 4.3);

    TCanvas *C_Y_LF_Muons = two_histo_ratio(h_Y_LF_DiMuon_Rec_PYTHIA8, h_Y_LF_DiMuon_Rec_PYTHIA6, Y_Ratio_PYTHIA8_PYTHIA6, Form("c_Y_%s", Name_Canvas.Data()), Title_Canvas.Data(), kTRUE, kTRUE);
    gPad = (TPad *)C_Y_LF_Muons->GetListOfPrimitives()->FindObject("pad2");
    gPad->cd();
    p_line->SetX1(-4.0);
    p_line->SetX2(-2.5);
    p_line->DrawClone("SAME");
}

void Pythia6_Pythia8_LF_Muons(Bool_t cross_section)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    Int_t Rebin_Eta = 5;
    Int_t Rebin_Pt = 1;

    const Int_t NCut = 3;
    TString Name_Histogram[NCut] = {"Muon_Gen/PYTHIA/h_PtEta_Muon_Gen_LF_PYTHIAOnly", "Muon_Gen/PYTHIA/h_PtEta_Muon_Gen_DQcut_LF_PYTHIAOnly", "Muon_Rec/PYTHIA/h_PtEta_Muon_RecLF_PYTHIAOnly"};
    TString Title_Canvas[NCut] = {"#splitline{Generated #mu #leftarrow LF hadrons (PYTHIA Particle)}{-4.5 < #eta < 4.5}", "#splitline{Generated #mu #leftarrow LF hadrons (PYTHIA Particle)}{#it{p}_{T}>0.5, 2.5 < #eta < 4}", "Reconstructed #mu #leftarrow LF hadrons (PYTHIA Particle)"};
    TString Name_Canvas[NCut] = {"cross_section_muongen_pythia6_8", "cross_section_muongen_dqcut_pythia6_8", "cross_section_muonrec_pythia6_8"};
    TLine *p_line = new TLine(0, 1, 20, 1);
    p_line->SetLineColor(kPink - 8);
    p_line->SetLineStyle(9);
    p_line->SetLineWidth(3);

    for (Int_t i_cut = 0; i_cut < NCut; i_cut++)
    {
        TH2F *h_PtEta_Muon_Gen_LF_PYTHIA6 = (TH2F *)fIn_Powheg_Charm->Get(Name_Histogram[i_cut].Data());
        if (i_cut == 0)
        {
            Rebin_Eta = 4;
            h_PtEta_Muon_Gen_LF_PYTHIA6->GetYaxis()->SetRangeUser(-4.5, 4.5);
        }
        else
            Rebin_Eta = 5;

        // h_PtEta_Muon_Gen_LF_PYTHIA6->GetYaxis()->SetRangeUser(-4.1,-2.4);
        TH1F *h_Eta_LF_Muon_PYTHIA6 = (TH1F *)h_PtEta_Muon_Gen_LF_PYTHIA6->ProjectionY();
        h_Eta_LF_Muon_PYTHIA6->SetName("h_Eta_LF_Muon_PYTHIA6");
        h_Eta_LF_Muon_PYTHIA6->SetTitle("from POWHEG+PYTHIA6 sim");
        if (cross_section)
        {
            h_Eta_LF_Muon_PYTHIA6->GetYaxis()->SetTitle("d#sigma/d#eta (mb)");
            hist1D_graphic_opt(h_Eta_LF_Muon_PYTHIA6, kTRUE, Rebin_Eta, 20, kAzure + 2, 46.4 / (1. * Nev_Powheg_Charm));
        }
        else
        {
            h_Eta_LF_Muon_PYTHIA6->GetYaxis()->SetTitle("entries");
            hist1D_graphic_opt(h_Eta_LF_Muon_PYTHIA6, kFALSE, Rebin_Eta, 20, kAzure + 2, 1. / (h_Eta_LF_Muon_PYTHIA6->GetEntries()));
        }

        TH1F *h_Pt_LF_Muon_PYTHIA6 = (TH1F *)h_PtEta_Muon_Gen_LF_PYTHIA6->ProjectionX();
        h_Pt_LF_Muon_PYTHIA6->SetName("h_Pt_LF_Muon_PYTHIA6");
        h_Pt_LF_Muon_PYTHIA6->SetTitle("from POWHEG+PYTHIA6 sim");
        if (cross_section)
        {
            h_Pt_LF_Muon_PYTHIA6->GetYaxis()->SetTitle("d#sigma/d#it{p}_{T} (mb/(GeV/c))");
            hist1D_graphic_opt(h_Pt_LF_Muon_PYTHIA6, kTRUE, Rebin_Pt, 20, kAzure + 2, 46.4 / (1. * Nev_Powheg_Charm));
        }
        else
        {
            h_Pt_LF_Muon_PYTHIA6->GetYaxis()->SetTitle("entries");
            hist1D_graphic_opt(h_Pt_LF_Muon_PYTHIA6, kFALSE, Rebin_Pt, 20, kAzure + 2, 1. / (h_Pt_LF_Muon_PYTHIA6->GetEntries()));
        }

        h_Pt_LF_Muon_PYTHIA6->GetXaxis()->SetRangeUser(0, 20);

        TH2F *h_PtEta_Muon_Gen_LF_PYTHIA8 = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get(Name_Histogram[i_cut].Data());
        if (i_cut == 0)
            h_PtEta_Muon_Gen_LF_PYTHIA8->GetYaxis()->SetRangeUser(-4.5, 4.5);

        // h_PtEta_Muon_Gen_LF_PYTHIA8->GetYaxis()->SetRangeUser(-4.1,-2.4);
        TH1F *h_Eta_LF_Muon_PYTHIA8 = (TH1F *)h_PtEta_Muon_Gen_LF_PYTHIA8->ProjectionY();
        h_Eta_LF_Muon_PYTHIA8->SetTitle("from PYTHIA8 MB sim");
        if (cross_section)
            hist1D_graphic_opt(h_Eta_LF_Muon_PYTHIA8, kTRUE, Rebin_Eta, 24, kAzure + 2, 78.05 / Nev_PYTHIA_MB_CENTRALIZED);
        else
            hist1D_graphic_opt(h_Eta_LF_Muon_PYTHIA8, kFALSE, Rebin_Eta, 24, kAzure + 2, 1.0 / h_Eta_LF_Muon_PYTHIA8->GetEntries());

        TH1F *h_Pt_LF_Muon_PYTHIA8 = (TH1F *)h_PtEta_Muon_Gen_LF_PYTHIA8->ProjectionX();
        h_Pt_LF_Muon_PYTHIA8->SetTitle("from PYTHIA8 MB sim");
        if (cross_section)
            hist1D_graphic_opt(h_Pt_LF_Muon_PYTHIA8, kTRUE, Rebin_Pt, 24, kAzure + 2, 78.05 / Nev_PYTHIA_MB_CENTRALIZED);
        else
            hist1D_graphic_opt(h_Pt_LF_Muon_PYTHIA8, kFALSE, Rebin_Pt, 24, kAzure + 2, 1.0 / h_Pt_LF_Muon_PYTHIA8->GetEntries());
        h_Pt_LF_Muon_PYTHIA8->GetXaxis()->SetRangeUser(0, 20);

        TH1F *Eta_Ratio_PYTHIA8_PYTHIA6 = (TH1F *)h_Eta_LF_Muon_PYTHIA8->Clone("Eta_Ratio_PYTHIA8_PYTHIA6");
        Eta_Ratio_PYTHIA8_PYTHIA6->Divide(h_Eta_LF_Muon_PYTHIA6);
        Eta_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetTitle("PYTHIA8/PYTHIA6");

        TCanvas *C_ETA_LF_Muons = two_histo_ratio(h_Eta_LF_Muon_PYTHIA6, h_Eta_LF_Muon_PYTHIA8, Eta_Ratio_PYTHIA8_PYTHIA6, Form("c_eta_%s", Name_Canvas[i_cut].Data()), Title_Canvas[i_cut], kTRUE, kTRUE);
        // C_ETA_LF_Muons->SaveAs(Form("images/%s.png", C_ETA_LF_Muons->GetName()));

        TH1F *Pt_Ratio_PYTHIA8_PYTHIA6 = (TH1F *)h_Pt_LF_Muon_PYTHIA8->Clone("Pt_Ratio_PYTHIA8_PYTHIA6");
        Pt_Ratio_PYTHIA8_PYTHIA6->Divide(h_Pt_LF_Muon_PYTHIA6);
        Pt_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetTitle("PYTHIA8/PYTHIA6");
        Pt_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetRangeUser(-0.3, 4.3);

        TCanvas *C_Pt_LF_Muons = two_histo_ratio(h_Pt_LF_Muon_PYTHIA6, h_Pt_LF_Muon_PYTHIA8, Pt_Ratio_PYTHIA8_PYTHIA6, Form("c_pt_%s", Name_Canvas[i_cut].Data()), Title_Canvas[i_cut], kTRUE, kTRUE);
        gPad = (TPad *)C_Pt_LF_Muons->GetListOfPrimitives()->FindObject("pad2");
        gPad->cd();
        p_line->Draw("SAME");
        // C_Pt_LF_Muons->SaveAs(Form("images/%s.png", C_Pt_LF_Muons->GetName()));

        new TCanvas();
        TH2F *h_ratio_PtEta_Muon_Gen_LF_PYTHIA6_8 = (TH2F *)h_PtEta_Muon_Gen_LF_PYTHIA8->Clone("h_ratio_PtEta_Muon_Gen_LF_PYTHIA6_8");
        // h_ratio_PtEta_Muon_Gen_LF_PYTHIA6_8->GetZaxis()->SetRangeUser(0.9,1);
        h_ratio_PtEta_Muon_Gen_LF_PYTHIA6_8->Divide(h_PtEta_Muon_Gen_LF_PYTHIA6);
        h_ratio_PtEta_Muon_Gen_LF_PYTHIA6_8->Draw("COLZ");
    }
}

void Pythia6_Pythia8_LF_Hadron()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    Int_t Rebin_Eta = 10;
    Int_t Rebin_Pt = 5;

    TH2F *h_PdgEta_HFHadron_prompt_PYTHIA6 = (TH2F *)fIn_Powheg_Charm->Get("Hadron/PYTHIA/h_PdgEta_HFHadron_prompt_PYTHIAOnly");
    h_PdgEta_HFHadron_prompt_PYTHIA6->GetXaxis()->SetRangeUser(210, 212);
    TH1F *h_Eta_HFHadron_prompt_PYTHIA6 = (TH1F *)h_PdgEta_HFHadron_prompt_PYTHIA6->ProjectionY();
    h_Eta_HFHadron_prompt_PYTHIA6->SetTitle("from POWHEG+PYTHIA6 sim");
    h_Eta_HFHadron_prompt_PYTHIA6->SetName("h_Eta_HFHadron_prompt_PYTHIA6");
    hist1D_graphic_opt(h_Eta_HFHadron_prompt_PYTHIA6, kTRUE, Rebin_Eta, 20, kAzure + 2, 46.4 / Nev_Powheg_Charm);

    TH2F *h_PdgPt_HFHadron_prompt_PYTHIA6 = (TH2F *)fIn_Powheg_Charm->Get("Hadron/PYTHIA/h_PdgPt_HFHadron_prompt_PYTHIAOnly");
    h_PdgPt_HFHadron_prompt_PYTHIA6->GetXaxis()->SetRangeUser(210, 212);
    TH1F *h_Pt_HFHadron_prompt_PYTHIA6 = (TH1F *)h_PdgPt_HFHadron_prompt_PYTHIA6->ProjectionY();
    h_Pt_HFHadron_prompt_PYTHIA6->GetXaxis()->SetRangeUser(0, 20);
    h_Pt_HFHadron_prompt_PYTHIA6->SetTitle("from POWHEG+PYTHIA6 sim");
    h_Pt_HFHadron_prompt_PYTHIA6->SetName("h_Pt_HFHadron_prompt_PYTHIA6");
    hist1D_graphic_opt(h_Pt_HFHadron_prompt_PYTHIA6, kTRUE, Rebin_Pt, 20, kAzure + 2, 46.4 / Nev_Powheg_Charm);

    TH2F *h_PdgEta_HFHadron_prompt_PYTHIA8 = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Hadron/Pythia/h_PdgEta_HFHadron_prompt_PYTHIAOnly");
    h_PdgEta_HFHadron_prompt_PYTHIA8->GetXaxis()->SetRangeUser(210, 212);
    TH1F *h_Eta_HFHadron_prompt_PYTHIA8 = (TH1F *)h_PdgEta_HFHadron_prompt_PYTHIA8->ProjectionY();
    h_Eta_HFHadron_prompt_PYTHIA8->SetTitle("from PYTHIA8 sim");
    hist1D_graphic_opt(h_Eta_HFHadron_prompt_PYTHIA8, kFALSE, Rebin_Eta, 24, kAzure + 2, 78.05 / Nev_PYTHIA_MB_CENTRALIZED);

    TH2F *h_PdgPt_HFHadron_prompt_PYTHIA8 = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Hadron/Pythia/h_PdgPt_HFHadron_prompt_PYTHIAOnly");
    h_PdgPt_HFHadron_prompt_PYTHIA8->GetXaxis()->SetRangeUser(210, 212);
    TH1F *h_Pt_HFHadron_prompt_PYTHIA8 = (TH1F *)h_PdgPt_HFHadron_prompt_PYTHIA8->ProjectionY();
    h_Pt_HFHadron_prompt_PYTHIA8->GetXaxis()->SetRangeUser(0, 20);
    h_Pt_HFHadron_prompt_PYTHIA8->SetTitle("from PYTHIA8 sim");
    hist1D_graphic_opt(h_Pt_HFHadron_prompt_PYTHIA8, kFALSE, Rebin_Pt, 24, kAzure + 2, 78.05 / Nev_PYTHIA_MB_CENTRALIZED);

    TH1F *Eta_Ratio_PYTHIA8_PYTHIA6 = (TH1F *)h_Eta_HFHadron_prompt_PYTHIA8->Clone("Eta_Ratio_PYTHIA8_PYTHIA6");
    Eta_Ratio_PYTHIA8_PYTHIA6->Divide(h_Eta_HFHadron_prompt_PYTHIA6);
    Eta_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetTitle("PYTHIA8/PYTHIA6");

    TCanvas *C_ETA_Pions = two_histo_ratio(h_Eta_HFHadron_prompt_PYTHIA6, h_Eta_HFHadron_prompt_PYTHIA8, Eta_Ratio_PYTHIA8_PYTHIA6, "c_eta_pions_pythia6_8", "Generated #pi (PYTHIA particles)", kTRUE, kTRUE);
    // C_ETA_Pions->SaveAs(Form("images/%s.png", C_ETA_Pions->GetName()));

    TH1F *Pt_Ratio_PYTHIA8_PYTHIA6 = (TH1F *)h_Pt_HFHadron_prompt_PYTHIA8->Clone("Pt_Ratio_PYTHIA8_PYTHIA6");
    Pt_Ratio_PYTHIA8_PYTHIA6->Divide(h_Pt_HFHadron_prompt_PYTHIA6);
    Pt_Ratio_PYTHIA8_PYTHIA6->GetYaxis()->SetTitle("PYTHIA8/PYTHIA6");

    TCanvas *C_Pt_Pions = two_histo_ratio(h_Pt_HFHadron_prompt_PYTHIA6, h_Pt_HFHadron_prompt_PYTHIA8, Pt_Ratio_PYTHIA8_PYTHIA6, "c_pt_pions_pythia6_8", "Generated #pi (PYTHIA particles)", kTRUE, kTRUE);
    // C_Pt_Pions->SaveAs(Form("images/%s.png", C_Pt_Pions->GetName()));
}

void Powheg_Gen_muons(TString gen)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    Int_t Rebin_eta = 15;
    Int_t Rebin_pt = 5;
    TString Muon_fromLF_CENTRALIZED[3];
    TString Muon_fromLF_CENTRALIZED_Title[3] = {"from Geant LF hadrons", "from PYTHIA LF hadrons", "from Powheg LF hadrons"};
    Color_t Muon_fromLF_CENTRALIZED_Color[3] = {kGreen + 2, kMagenta + 2, kOrange + 9};

    Muon_fromLF_CENTRALIZED[0].Form("Geant");
    Muon_fromLF_CENTRALIZED[1].Form("PYTHIA");
    Muon_fromLF_CENTRALIZED[2].Form("Powheg");

    TH2F *h_PtEta_Muon_Gen_CENTRALIZED[3];
    TH1F *h_Pt_Muon_Gen_CENTRALIZED[3];
    TH1F *h_Eta_Muon_Gen_CENTRALIZED[3];
    cout << Nev_Powheg_Charm << endl;
    for (Int_t i_origin = 0; i_origin < 3; i_origin++)
    {
        TString name;
        if (gen.Contains("Gen"))
            name.Form("Muon_%s/%s/h_PtEta_Muon_%s_LF_%sOnly", gen.Data(), Muon_fromLF_CENTRALIZED[i_origin].Data(), gen.Data(), Muon_fromLF_CENTRALIZED[i_origin].Data());
        else
            name.Form("Muon_%s/%s/h_PtEta_Muon_%sLF_%sOnly", gen.Data(), Muon_fromLF_CENTRALIZED[i_origin].Data(), gen.Data(), Muon_fromLF_CENTRALIZED[i_origin].Data());

        h_PtEta_Muon_Gen_CENTRALIZED[i_origin] = (TH2F *)fIn_Powheg_Charm->Get(name.Data());
        printf("%s\n", name.Data());
        h_PtEta_Muon_Gen_CENTRALIZED[i_origin]->Draw("COLZ");
        if (h_PtEta_Muon_Gen_CENTRALIZED[i_origin] == nullptr)
            continue;
        cout << h_PtEta_Muon_Gen_CENTRALIZED[i_origin]->GetName() << " Entries :" << h_PtEta_Muon_Gen_CENTRALIZED[i_origin]->GetEntries() << endl;
        h_Pt_Muon_Gen_CENTRALIZED[i_origin] = (TH1F *)h_PtEta_Muon_Gen_CENTRALIZED[i_origin]->ProjectionX();
        hist1D_graphic_opt(h_Pt_Muon_Gen_CENTRALIZED[i_origin], kFALSE, Rebin_pt, 20, Muon_fromLF_CENTRALIZED_Color[i_origin], 1. / Nev_Powheg_Charm);
        h_Pt_Muon_Gen_CENTRALIZED[i_origin]->GetYaxis()->SetTitle("Entries x ev.");
        h_Pt_Muon_Gen_CENTRALIZED[i_origin]->SetTitle(Muon_fromLF_CENTRALIZED_Title[i_origin]);
        h_Pt_Muon_Gen_CENTRALIZED[i_origin]->GetXaxis()->SetRangeUser(0, 30);

        h_Eta_Muon_Gen_CENTRALIZED[i_origin] = (TH1F *)h_PtEta_Muon_Gen_CENTRALIZED[i_origin]->ProjectionY();
        hist1D_graphic_opt(h_Eta_Muon_Gen_CENTRALIZED[i_origin], kFALSE, Rebin_eta, 20, Muon_fromLF_CENTRALIZED_Color[i_origin], 1. / Nev_Powheg_Charm);
        h_Eta_Muon_Gen_CENTRALIZED[i_origin]->GetYaxis()->SetTitle("Entries x ev.");
        h_Eta_Muon_Gen_CENTRALIZED[i_origin]->SetTitle(Muon_fromLF_CENTRALIZED_Title[i_origin]);
    }

    TH1F *Pt_Ratio_Powheg_PYTHIA = (TH1F *)h_Pt_Muon_Gen_CENTRALIZED[2]->Clone("Pt_Ratio_Powheg_PYTHIA");
    Pt_Ratio_Powheg_PYTHIA->Divide(h_Pt_Muon_Gen_CENTRALIZED[1]);
    Pt_Ratio_Powheg_PYTHIA->SetTitle("Powheg/PYTHIA");

    TH1F *Pt_Ratio_Geant_PYTHIA = (TH1F *)h_Pt_Muon_Gen_CENTRALIZED[0]->Clone("Pt_Ratio_Geant_PYTHIA");
    Pt_Ratio_Geant_PYTHIA->Divide(h_Pt_Muon_Gen_CENTRALIZED[1]);
    Pt_Ratio_Geant_PYTHIA->SetTitle("Geant/PYTHIA");

    TCanvas *C_Pt_MuGen_PowhegSim;
    if (gen.Contains("Gen"))
        C_Pt_MuGen_PowhegSim = three_histo_ratio(h_Pt_Muon_Gen_CENTRALIZED[1], h_Pt_Muon_Gen_CENTRALIZED[2], h_Pt_Muon_Gen_CENTRALIZED[0], Pt_Ratio_Powheg_PYTHIA, Pt_Ratio_Geant_PYTHIA, "C_Pt_MuGen_PowhegSim", "#splitline{POWHEG+PYTHIA6 pp @13 Tev}{Generated #mu #leftarrow LF hadrons}", kTRUE);
    else
        C_Pt_MuGen_PowhegSim = three_histo_ratio(h_Pt_Muon_Gen_CENTRALIZED[1], h_Pt_Muon_Gen_CENTRALIZED[2], h_Pt_Muon_Gen_CENTRALIZED[0], Pt_Ratio_Powheg_PYTHIA, Pt_Ratio_Geant_PYTHIA, "C_Pt_MuRec_PowhegSim", "#splitline{POWHEG+PYTHIA6 pp @13 Tev}{Reconstructed #mu #leftarrow LF hadrons}", kTRUE);

    C_Pt_MuGen_PowhegSim->SaveAs(Form("images/C_Pt_Mu%s_PowhegSim.png", gen.Data()));

    TH1F *Eta_Ratio_Powheg_PYTHIA = (TH1F *)h_Eta_Muon_Gen_CENTRALIZED[2]->Clone("Eta_Ratio_Powheg_PYTHIA");
    Eta_Ratio_Powheg_PYTHIA->Divide(h_Eta_Muon_Gen_CENTRALIZED[1]);
    Eta_Ratio_Powheg_PYTHIA->SetTitle("Powheg/PYTHIA");

    TH1F *Eta_Ratio_Geant_PYTHIA = (TH1F *)h_Eta_Muon_Gen_CENTRALIZED[0]->Clone("Eta_Ratio_Geant_PYTHIA");
    Eta_Ratio_Geant_PYTHIA->Divide(h_Eta_Muon_Gen_CENTRALIZED[1]);
    Eta_Ratio_Geant_PYTHIA->SetTitle("Geant/PYTHIA");

    Eta_Ratio_Powheg_PYTHIA->GetYaxis()->SetRangeUser(2.3e-02, 3.6);
    h_Eta_Muon_Gen_CENTRALIZED[1]->GetYaxis()->SetRangeUser(1e-06, 2e-02);
    TCanvas *C_Eta_MuGen_PowhegSim;

    if (gen.Contains("Gen"))
        C_Eta_MuGen_PowhegSim = three_histo_ratio(h_Eta_Muon_Gen_CENTRALIZED[1], h_Eta_Muon_Gen_CENTRALIZED[2], h_Eta_Muon_Gen_CENTRALIZED[0], Eta_Ratio_Powheg_PYTHIA, Eta_Ratio_Geant_PYTHIA, "C_Eta_MuGen_PowhegSim", "#splitline{POWHEG+PYTHIA6 pp @13 Tev}{Generated #mu #leftarrow LF hadrons}", kTRUE);

    else
        C_Eta_MuGen_PowhegSim = three_histo_ratio(h_Eta_Muon_Gen_CENTRALIZED[1], h_Eta_Muon_Gen_CENTRALIZED[2], h_Eta_Muon_Gen_CENTRALIZED[0], Eta_Ratio_Powheg_PYTHIA, Eta_Ratio_Geant_PYTHIA, "C_Eta_MuRec_PowhegSim", "#splitline{POWHEG+PYTHIA6 pp @13 Tev}{Reconstructed #mu #leftarrow LF hadrons}", kTRUE);

    TPad *g = (TPad *)C_Eta_MuGen_PowhegSim->GetListOfPrimitives()->FindObject("pad2");
    g->SetLogy();

    TLegend *legend1 = (TLegend *)C_Eta_MuGen_PowhegSim->FindObject("legend1");
    legend1->SetY1(0.05);
    legend1->SetY2(0.35);

    TLegend *legend2 = (TLegend *)C_Eta_MuGen_PowhegSim->FindObject("legend2");
    legend2->SetY1(0.3);
    legend2->SetY2(0.5);

    C_Eta_MuGen_PowhegSim->Update();
    C_Eta_MuGen_PowhegSim->SaveAs(Form("images/C_Eta_Mu%s_PowhegSim.png", gen.Data()));
}

void Gen_muons()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    Int_t Rebin = 5;
    fIn_PYTHIA_CENTRALIZED_MB->cd("Muon_Gen/Geant");
    fIn_PYTHIA_CENTRALIZED_MB->ls();
    fIn_PYTHIA_CENTRALIZED_MB->cd("Muon_Gen/Pythia");
    fIn_PYTHIA_CENTRALIZED_MB->ls();

    TString Muon_fromLF_CENTRALIZED[2];
    TString Muon_fromLF_CENTRALIZED_Title[2] = {"from Geant LF hadrons", "from PYTHIA LF hadrons"};
    Color_t Muon_fromLF_CENTRALIZED_Color[2] = {kGreen + 2, kMagenta + 2};

    Muon_fromLF_CENTRALIZED[0].Form("Geant");
    Muon_fromLF_CENTRALIZED[1].Form("PYTHIA");

    TH2F *h_PtEta_Muon_Gen_CENTRALIZED[2];
    TH1F *h_Pt_Muon_Gen_CENTRALIZED[2];
    TH1F *h_Eta_Muon_Gen_CENTRALIZED[2];

    for (Int_t i_origin = 0; i_origin < 2; i_origin++)
    {
        h_PtEta_Muon_Gen_CENTRALIZED[i_origin] = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get(Form("Muon_Gen/%s/h_PtEta_Muon_Gen_LF_%sOnly", Muon_fromLF_CENTRALIZED[i_origin].Data(), Muon_fromLF_CENTRALIZED[i_origin].Data()));
        if (h_PtEta_Muon_Gen_CENTRALIZED[i_origin] == nullptr)
            continue;

        h_Pt_Muon_Gen_CENTRALIZED[i_origin] = (TH1F *)h_PtEta_Muon_Gen_CENTRALIZED[i_origin]->ProjectionX();
        hist1D_graphic_opt(h_Pt_Muon_Gen_CENTRALIZED[i_origin], kFALSE, Rebin, 20, Muon_fromLF_CENTRALIZED_Color[i_origin], 1. / Nev_PYTHIA_MB_CENTRALIZED);
        h_Pt_Muon_Gen_CENTRALIZED[i_origin]->GetYaxis()->SetTitle("Entries x ev.");
        h_Pt_Muon_Gen_CENTRALIZED[i_origin]->SetTitle(Muon_fromLF_CENTRALIZED_Title[i_origin]);
        h_Pt_Muon_Gen_CENTRALIZED[i_origin]->GetXaxis()->SetRangeUser(0, 30);

        h_Eta_Muon_Gen_CENTRALIZED[i_origin] = (TH1F *)h_PtEta_Muon_Gen_CENTRALIZED[i_origin]->ProjectionY();
        hist1D_graphic_opt(h_Eta_Muon_Gen_CENTRALIZED[i_origin], kFALSE, Rebin, 20, Muon_fromLF_CENTRALIZED_Color[i_origin], 1. / Nev_PYTHIA_MB_CENTRALIZED);
        h_Eta_Muon_Gen_CENTRALIZED[i_origin]->GetYaxis()->SetTitle("Entries x ev.");
        h_Eta_Muon_Gen_CENTRALIZED[i_origin]->SetTitle(Muon_fromLF_CENTRALIZED_Title[i_origin]);
    }

    TH2F *h_PtEta_Muon_Gen_STANDALONE_INEL = (TH2F *)fIn_Pythia_STANDALONE_INEL->Get("Muon_Gen/h_PtEta_Muon_Gen_LF");
    TH1F *h_Pt_Muon_Gen_STANDALONE_INEL = (TH1F *)h_PtEta_Muon_Gen_STANDALONE_INEL->ProjectionX();
    hist1D_graphic_opt(h_Pt_Muon_Gen_STANDALONE_INEL, kFALSE, Rebin, 20, kAzure + 4, 1. / Nev_PYTHIA_STANDALONE_INEL);
    h_Pt_Muon_Gen_STANDALONE_INEL->GetXaxis()->SetRangeUser(0, 30);

    h_Pt_Muon_Gen_STANDALONE_INEL->SetTitle("PYTHIA std. Inel");
    TH1F *h_Eta_Muon_Gen_STANDALONE_INEL = (TH1F *)h_PtEta_Muon_Gen_STANDALONE_INEL->ProjectionY();
    hist1D_graphic_opt(h_Eta_Muon_Gen_STANDALONE_INEL, kFALSE, Rebin, 20, kAzure + 4, 1. / Nev_PYTHIA_STANDALONE_INEL);
    h_Eta_Muon_Gen_STANDALONE_INEL->SetTitle("PYTHIA std. Inel");

    h_Pt_Muon_Gen_CENTRALIZED[1]->SetTitle("PYTHIA centr. (from PYTHIA hadrons)");
    TH1F *Pt_Ratio_CENTRALIZED_STANDALONE = (TH1F *)h_Pt_Muon_Gen_CENTRALIZED[1]->Clone("Pt_Ratio_CENTRALIZED_STANDALONE");
    Pt_Ratio_CENTRALIZED_STANDALONE->Divide(h_Pt_Muon_Gen_STANDALONE_INEL);
    Pt_Ratio_CENTRALIZED_STANDALONE->GetYaxis()->SetTitle("CENTRALIZED/STD. INEL");
    Pt_Ratio_CENTRALIZED_STANDALONE->GetYaxis()->CenterTitle();

    TCanvas *c_Pt_CENTR_STD = two_histo_ratio(h_Pt_Muon_Gen_CENTRALIZED[1], h_Pt_Muon_Gen_STANDALONE_INEL, Pt_Ratio_CENTRALIZED_STANDALONE, "C_Pt_Mu_GEN_CENTR_STD", "Generated #mu #leftarrow LF", kTRUE, kTRUE);
    c_Pt_CENTR_STD->SaveAs(Form("%s_LF.png", c_Pt_CENTR_STD->GetName()));

    h_Eta_Muon_Gen_CENTRALIZED[1]->SetTitle("PYTHIA centr. (from PYTHIA hadrons)");
    TH1F *Eta_Ratio_CENTRALIZED_STANDALONE = (TH1F *)h_Eta_Muon_Gen_CENTRALIZED[1]->Clone("Eta_Ratio_CENTRALIZED_STANDALONE");
    Eta_Ratio_CENTRALIZED_STANDALONE->Divide(h_Eta_Muon_Gen_STANDALONE_INEL);
    Eta_Ratio_CENTRALIZED_STANDALONE->GetYaxis()->SetTitle("CENTRALIZED/STD. INEL");
    Eta_Ratio_CENTRALIZED_STANDALONE->GetYaxis()->CenterTitle();

    TCanvas *c_Eta_CENTR_STD = two_histo_ratio(h_Eta_Muon_Gen_CENTRALIZED[1], h_Eta_Muon_Gen_STANDALONE_INEL, Eta_Ratio_CENTRALIZED_STANDALONE, "C_Eta_Mu_GEN_CENTR_STD", "Generated #mu #leftarrow LF", kTRUE, kTRUE);
    // c_Eta_CENTR_STD->SaveAs(Form("%s_LF.png",c_Eta_CENTR_STD->GetName()));

    if (h_PtEta_Muon_Gen_CENTRALIZED[0] != nullptr)
    {
        TH1F *Pt_Ratio_CENTRALIZED_Geant_PYTHIA = (TH1F *)h_Pt_Muon_Gen_CENTRALIZED[0]->Clone("Pt_Ratio_CENTRALIZED_Geant_PYTHIA");
        Pt_Ratio_CENTRALIZED_Geant_PYTHIA->Divide(h_Pt_Muon_Gen_CENTRALIZED[1]);
        Pt_Ratio_CENTRALIZED_Geant_PYTHIA->GetYaxis()->SetTitle("Geant/PYTHIA");
        Pt_Ratio_CENTRALIZED_Geant_PYTHIA->GetYaxis()->CenterTitle();
        h_Pt_Muon_Gen_CENTRALIZED[1]->GetYaxis()->SetRangeUser(1.2e-06, h_Pt_Muon_Gen_CENTRALIZED[1]->GetMaximum() * 4000);
        h_Pt_Muon_Gen_CENTRALIZED[1]->SetTitle(Muon_fromLF_CENTRALIZED_Title[1]);

        TCanvas *c_Pt_Geant_PYTHIA = two_histo_ratio(h_Pt_Muon_Gen_CENTRALIZED[1], h_Pt_Muon_Gen_CENTRALIZED[0], Pt_Ratio_CENTRALIZED_Geant_PYTHIA, "c_Pt_Geant_PYTHIA", "Generated #mu #leftarrow LF centralized sim", kTRUE, kTRUE);
        c_Pt_Geant_PYTHIA->SaveAs(Form("%s_LF.png", c_Pt_Geant_PYTHIA->GetName()));

        TH1F *Eta_Ratio_CENTRALIZED_Geant_PYTHIA = (TH1F *)h_Eta_Muon_Gen_CENTRALIZED[0]->Clone("Eta_Ratio_CENTRALIZED_Geant_PYTHIA");
        Eta_Ratio_CENTRALIZED_Geant_PYTHIA->Divide(h_Eta_Muon_Gen_CENTRALIZED[1]);
        Eta_Ratio_CENTRALIZED_Geant_PYTHIA->GetYaxis()->SetTitle("Geant/PYTHIA");
        Eta_Ratio_CENTRALIZED_Geant_PYTHIA->GetYaxis()->CenterTitle();

        h_Eta_Muon_Gen_CENTRALIZED[1]->GetYaxis()->SetRangeUser(1.2e-06, h_Eta_Muon_Gen_CENTRALIZED[1]->GetMaximum() * 4000);
        h_Eta_Muon_Gen_CENTRALIZED[1]->SetTitle(Muon_fromLF_CENTRALIZED_Title[1]);

        TCanvas *c_Eta_Geant_PYTHIA = two_histo_ratio(h_Eta_Muon_Gen_CENTRALIZED[1], h_Eta_Muon_Gen_CENTRALIZED[0], Eta_Ratio_CENTRALIZED_Geant_PYTHIA, "c_Eta_Geant_PYTHIA", "Generated #mu #leftarrow LF centralized sim", kTRUE, kTRUE);
        c_Eta_Geant_PYTHIA->SaveAs(Form("%s_LF.png", c_Eta_Geant_PYTHIA->GetName()));
    }
}

void Rec_Muons()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    fIn_PYTHIA_CENTRALIZED_MB->cd("Muon_Rec/Geant");
    fIn_PYTHIA_CENTRALIZED_MB->ls();
    std::cout << Nev_PYTHIA_MB_CENTRALIZED << endl;
    // fIn_PYTHIA_CENTRALIZED_MB->cd("Muon_Rec/not_Geant");
    // fIn_PYTHIA_CENTRALIZED_MB->ls();
    TH2F *h_VzmotherEta_Muon_Rec_PYTHIA_LF = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/PYTHIA/h_VzmotherEta_Muon_Rec_PYTHIA_LF");
    TH2F *h_VzmotherEta_Muon_Rec_Geant_LF = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/Geant/h_VzmotherEta_Muon_Rec_Geant_LF");

    TCanvas *C_VzmotherEta = canvas_noratio_divide2("C_VzmotherEta");
    C_VzmotherEta->cd();
    h_VzmotherEta_Muon_Rec_PYTHIA_LF->Draw("COLZ");
    h_VzmotherEta_Muon_Rec_Geant_LF->Draw("COLZSAME0");

    // TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.0, 1.0, 1.0);
    // pad1->Draw();
    // TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 1.0);
    // pad2->SetFillStyle(4000);
    // pad2->Draw();

    // gStyle->SetPalette(kSolar);
    // pad1->cd();
    // gStyle->SetPalette(kBird);
    // pad2->cd();
    // h_VzmotherEta_Muon_Rec_Geant_LF->Draw("COLZ");

    TH1F *h_Nperevent_Muon_Rec_LF_PYTHIAOnly = (TH1F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/PYTHIA/h_Nperevent_Muon_Rec_LF_PYTHIAOnly");
    hist1D_graphic_opt(h_Nperevent_Muon_Rec_LF_PYTHIAOnly, kTRUE, 1, 20, kMagenta + 2, 1.);
    h_Nperevent_Muon_Rec_LF_PYTHIAOnly->SetTitle("from PYTHIA hadrons");
    h_Nperevent_Muon_Rec_LF_PYTHIAOnly->GetYaxis()->SetTitle("Events");

    TH1F *h_Nperevent_Muon_Rec_LF_GeantOnly = (TH1F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/Geant/h_Nperevent_Muon_Rec_LF_GeantOnly");
    hist1D_graphic_opt(h_Nperevent_Muon_Rec_LF_GeantOnly, kTRUE, 1, 20, kGreen + 2, 1.);
    h_Nperevent_Muon_Rec_LF_GeantOnly->SetTitle("from GeantO hadrons");
    TCanvas *C_Nperevent = canvas_noratio_divide2("C_Nperevent");
    C_Nperevent->cd();
    gPad->SetLogy();
    h_Nperevent_Muon_Rec_LF_PYTHIAOnly->Draw("HIST");
    h_Nperevent_Muon_Rec_LF_GeantOnly->Draw("HISTSAME");
    gPad->BuildLegend(0.5, 0.5, 0.8, 0.8, " ", 0);
    C_Nperevent->SaveAs("images/LF_MuRec_xevent.png");

    TH2F *h_PtEta_Muon_Rec_LF_PYTHIAOnly = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/PYTHIA/h_PtEta_Muon_RecLF_PYTHIAOnly");
    TH1F *h_Pt_Muon_Rec_LF_PYTHIAOnly = (TH1F *)h_PtEta_Muon_Rec_LF_PYTHIAOnly->ProjectionX();
    hist1D_graphic_opt(h_Pt_Muon_Rec_LF_PYTHIAOnly, kTRUE, 1, 20, kMagenta + 2, 1.);
    h_Pt_Muon_Rec_LF_PYTHIAOnly->GetXaxis()->SetRangeUser(0, 16);
    h_Pt_Muon_Rec_LF_PYTHIAOnly->SetTitle("from PYTHIA hadron");
    h_Pt_Muon_Rec_LF_PYTHIAOnly->GetYaxis()->SetTitle("Entries");
    cout << "h_Pt_Muon_Rec_LF_PYTHIAOnly Entries() " << h_Pt_Muon_Rec_LF_PYTHIAOnly->GetEntries() << endl;

    TH1F *h_Eta_Muon_Rec_LF_PYTHIAOnly = (TH1F *)h_PtEta_Muon_Rec_LF_PYTHIAOnly->ProjectionY();
    hist1D_graphic_opt(h_Eta_Muon_Rec_LF_PYTHIAOnly, kTRUE, 1, 20, kMagenta + 2, 1.);
    h_Eta_Muon_Rec_LF_PYTHIAOnly->SetTitle("from PYTHIA hadron");
    h_Eta_Muon_Rec_LF_PYTHIAOnly->GetYaxis()->SetTitle("Entries");

    TH2F *h_PtEta_Muon_Rec_LF_GeantOnly = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/Geant/h_PtEta_Muon_RecLF_GeantOnly");
    TH1F *h_Pt_Muon_Rec_LF_GeantOnly = (TH1F *)h_PtEta_Muon_Rec_LF_GeantOnly->ProjectionX();
    hist1D_graphic_opt(h_Pt_Muon_Rec_LF_GeantOnly, kTRUE, 1, 20, kGreen + 2, 1.);
    h_Pt_Muon_Rec_LF_GeantOnly->GetXaxis()->SetRangeUser(0, 16);
    h_Pt_Muon_Rec_LF_GeantOnly->SetTitle("from Geant hadron");
    h_Pt_Muon_Rec_LF_GeantOnly->GetYaxis()->SetTitle("Entries");
    cout << "h_Pt_Muon_Rec_LF_GeantOnly Entries() " << h_Pt_Muon_Rec_LF_GeantOnly->GetEntries() << endl;

    TH1F *h_Eta_Muon_Rec_LF_GeantOnly = (TH1F *)h_PtEta_Muon_Rec_LF_GeantOnly->ProjectionY();
    hist1D_graphic_opt(h_Eta_Muon_Rec_LF_GeantOnly, kTRUE, 1, 20, kGreen + 2, 1.);
    h_Eta_Muon_Rec_LF_GeantOnly->SetTitle("from Geant hadron");
    h_Eta_Muon_Rec_LF_GeantOnly->GetYaxis()->SetTitle("Entries");

    TH1F *ratio_Pt_Geant_PYTHIA_CENTRALIZED = (TH1F *)h_Pt_Muon_Rec_LF_GeantOnly->Clone("ratio_Pt_Geant_PYTHIA_CENTRALIZED");
    ratio_Pt_Geant_PYTHIA_CENTRALIZED->Divide(h_Pt_Muon_Rec_LF_PYTHIAOnly);
    ratio_Pt_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("Geant/PYTHIA");

    h_Pt_Muon_Rec_LF_PYTHIAOnly->GetYaxis()->SetRangeUser(0.8, h_Pt_Muon_Rec_LF_PYTHIAOnly->GetMaximum() * 400.8);
    ratio_Pt_Geant_PYTHIA_CENTRALIZED->GetYaxis()->CenterTitle();
    ratio_Pt_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(ratio_Pt_Geant_PYTHIA_CENTRALIZED->GetMinimum() * 0.3, ratio_Pt_Geant_PYTHIA_CENTRALIZED->GetMaximum() * 1.2);
    TCanvas *c_Pt_Rec = two_histo_ratio(h_Pt_Muon_Rec_LF_PYTHIAOnly, h_Pt_Muon_Rec_LF_GeantOnly, ratio_Pt_Geant_PYTHIA_CENTRALIZED, "c_Pt_Rec", "#splitline{PYTHIA8 pp @13 TeV}{Reconstructed #mu #leftarrow LF hadrons}", kTRUE, kTRUE);
    c_Pt_Rec->SaveAs("images/LF_MuRec_Pt.png");

    TH1F *ratio_Y_Geant_PYTHIA_CENTRALIZED = (TH1F *)h_Eta_Muon_Rec_LF_GeantOnly->Clone("ratio_Y_Geant_PYTHIA_CENTRALIZED");
    ratio_Y_Geant_PYTHIA_CENTRALIZED->Divide(h_Eta_Muon_Rec_LF_PYTHIAOnly);
    ratio_Y_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("Geant/PYTHIA");

    h_Eta_Muon_Rec_LF_PYTHIAOnly->GetYaxis()->SetRangeUser(0.8, h_Eta_Muon_Rec_LF_PYTHIAOnly->GetMaximum() * 400.8);
    ratio_Y_Geant_PYTHIA_CENTRALIZED->GetYaxis()->CenterTitle();
    ratio_Y_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(ratio_Y_Geant_PYTHIA_CENTRALIZED->GetMinimum() * 0.3, ratio_Y_Geant_PYTHIA_CENTRALIZED->GetMaximum() * 1.2);
    TCanvas *c_Y_Rec = two_histo_ratio(h_Eta_Muon_Rec_LF_PYTHIAOnly, h_Eta_Muon_Rec_LF_GeantOnly, ratio_Y_Geant_PYTHIA_CENTRALIZED, "c_Y_Rec", "#splitline{PYTHIA8 pp @13 TeV}{Reconstructed #mu #leftarrow LF hadrons}", kTRUE, kTRUE);
    c_Y_Rec->SaveAs("images/LF_MuRec_Eta.png");
    // c_ETA_Rec->SaveAs(Form("images/%s.png", c_ETA_Rec->GetName()));
}

void Rec_Dimu(TString gen)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    Int_t Rebin_Y = 2;
    Int_t Rebin_Pt = 2;
    Int_t Rebin_M = 2;
    TString *DiMuon_fromLF_Generator;
    TString *DiMuon_fromLF_Generator_Title;
    Color_t *DiMuon_fromLF_Generator_Color;

    TH2F **h_PtM_DiMuon_Rec_fromLF;
    TH2F **h_PtY_DiMuon_Rec_fromLF;
    TH1F **h_Pt_DiMuon_Rec_fromLF;
    TH1F **h_M_DiMuon_Rec_fromLF;
    TH1F **h_Y_DiMuon_Rec_fromLF;

    Int_t N_Gen = 0;
    TCanvas *C_Y_DiMuRec = canvas_noratio_divide2(Form("C_Y_DiMuRec_%s", gen.Data()));
    C_Y_DiMuRec->cd();
    gPad->SetLogy();

    TCanvas *C_Pt_DiMuRec = canvas_noratio_divide2(Form("C_Pt_DiMuRec_%s", gen.Data()));
    C_Pt_DiMuRec->cd();
    gPad->SetLogy();

    TCanvas *C_M_DiMuRec = canvas_noratio_divide2(Form("C_M_DiMuRec_%s", gen.Data()));
    C_M_DiMuRec->cd();
    gPad->SetLogy();

    if (gen.Contains("Pythia"))
    {
        const Int_t N_GenPythia = 3;

        h_PtM_DiMuon_Rec_fromLF = new TH2F *[N_GenPythia];
        h_PtY_DiMuon_Rec_fromLF = new TH2F *[N_GenPythia];
        h_Pt_DiMuon_Rec_fromLF = new TH1F *[N_GenPythia];
        h_M_DiMuon_Rec_fromLF = new TH1F *[N_GenPythia];
        h_Y_DiMuon_Rec_fromLF = new TH1F *[N_GenPythia];

        DiMuon_fromLF_Generator = new TString[N_GenPythia]{"GeantOnly", "PYTHIAGeant", "PYTHIAOnly"};
        DiMuon_fromLF_Generator_Title = new TString[N_GenPythia]{"#mu#mu <- Geant LF", "#mu#mu <- Geant,PYTHIA LF", "#mu#mu <- PYTHIA LF"};
        DiMuon_fromLF_Generator_Color = new Color_t[N_GenPythia]{kGreen + 2, kRed, kMagenta + 2};
        N_Gen = N_GenPythia;
    }
    else if (gen.Contains("Powheg"))
    {
        const Int_t N_GenPowheg = 3;
        h_PtM_DiMuon_Rec_fromLF = new TH2F *[N_GenPowheg];
        h_PtY_DiMuon_Rec_fromLF = new TH2F *[N_GenPowheg];
        h_Pt_DiMuon_Rec_fromLF = new TH1F *[N_GenPowheg];
        h_M_DiMuon_Rec_fromLF = new TH1F *[N_GenPowheg];
        h_Y_DiMuon_Rec_fromLF = new TH1F *[N_GenPowheg];

        DiMuon_fromLF_Generator = new TString[N_GenPowheg]{"PYTHIAOnly", "GeantOnly", "PYTHIAGeant"};
        DiMuon_fromLF_Generator_Title = new TString[N_GenPowheg]{"#mu#mu <- Geant LF", "#mu#mu <- Geant,PYTHIA LF", "#mu#mu <- PYTHIA LF"};
        DiMuon_fromLF_Generator_Color = new Color_t[N_GenPowheg]{kMagenta + 2, kGreen + 2, kRed};
        N_Gen = N_GenPowheg;
    }

    for (Int_t i_LF_Generator = 0; i_LF_Generator < N_Gen; i_LF_Generator++)
    {
        if (gen.Contains("Pythia"))
        {
            h_PtM_DiMuon_Rec_fromLF[i_LF_Generator] = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get(Form("DiMuon_Rec/LF_origin/h_PtM_DiMuon_Rec_fromLF_%s", DiMuon_fromLF_Generator[i_LF_Generator].Data()));
            h_PtY_DiMuon_Rec_fromLF[i_LF_Generator] = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get(Form("DiMuon_Rec/LF_origin/h_PtY_DiMuon_Rec_fromLF_%s", DiMuon_fromLF_Generator[i_LF_Generator].Data()));
        }
        else if (gen.Contains("Powheg"))
        {
            h_PtM_DiMuon_Rec_fromLF[i_LF_Generator] = (TH2F *)fIn_Powheg_Charm->Get(Form("DiMuon_Rec/LF_origin/h_PtM_DiMuon_Rec_fromLF_%s", DiMuon_fromLF_Generator[i_LF_Generator].Data()));
            h_PtY_DiMuon_Rec_fromLF[i_LF_Generator] = (TH2F *)fIn_Powheg_Charm->Get(Form("DiMuon_Rec/LF_origin/h_PtY_DiMuon_Rec_fromLF_%s", DiMuon_fromLF_Generator[i_LF_Generator].Data()));
        }

        cout << h_PtY_DiMuon_Rec_fromLF[i_LF_Generator]->GetName() << " Entries :" << h_PtY_DiMuon_Rec_fromLF[i_LF_Generator]->GetEntries() << endl;
        h_Pt_DiMuon_Rec_fromLF[i_LF_Generator] = (TH1F *)h_PtM_DiMuon_Rec_fromLF[i_LF_Generator]->ProjectionX();
        h_Pt_DiMuon_Rec_fromLF[i_LF_Generator]->SetTitle(DiMuon_fromLF_Generator_Title[i_LF_Generator].Data());
        hist1D_graphic_opt(h_Pt_DiMuon_Rec_fromLF[i_LF_Generator], kFALSE, Rebin_Y, 20, DiMuon_fromLF_Generator_Color[i_LF_Generator], 1.);
        h_Pt_DiMuon_Rec_fromLF[i_LF_Generator]->GetXaxis()->SetRangeUser(0., 10.);
        h_Pt_DiMuon_Rec_fromLF[i_LF_Generator]->GetYaxis()->SetTitle("Entries");
        h_Pt_DiMuon_Rec_fromLF[i_LF_Generator]->GetYaxis()->SetTitleOffset(1.2);
        h_Pt_DiMuon_Rec_fromLF[i_LF_Generator]->GetXaxis()->SetTitleOffset(1.15);

        h_M_DiMuon_Rec_fromLF[i_LF_Generator] = (TH1F *)h_PtM_DiMuon_Rec_fromLF[i_LF_Generator]->ProjectionY();
        h_M_DiMuon_Rec_fromLF[i_LF_Generator]->SetTitle(DiMuon_fromLF_Generator_Title[i_LF_Generator].Data());
        hist1D_graphic_opt(h_M_DiMuon_Rec_fromLF[i_LF_Generator], kFALSE, Rebin_Pt, 20, DiMuon_fromLF_Generator_Color[i_LF_Generator], 1.);
        h_M_DiMuon_Rec_fromLF[i_LF_Generator]->GetXaxis()->SetRangeUser(0., 10.);
        h_M_DiMuon_Rec_fromLF[i_LF_Generator]->GetYaxis()->SetTitle("Entries");
        h_M_DiMuon_Rec_fromLF[i_LF_Generator]->GetYaxis()->SetTitleOffset(1.2);
        h_M_DiMuon_Rec_fromLF[i_LF_Generator]->GetXaxis()->SetTitleOffset(1.15);

        h_Y_DiMuon_Rec_fromLF[i_LF_Generator] = (TH1F *)h_PtY_DiMuon_Rec_fromLF[i_LF_Generator]->ProjectionY();
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->SetTitle(DiMuon_fromLF_Generator_Title[i_LF_Generator].Data());
        hist1D_graphic_opt(h_Y_DiMuon_Rec_fromLF[i_LF_Generator], kFALSE, Rebin_M, 20, DiMuon_fromLF_Generator_Color[i_LF_Generator], 1.);
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->GetXaxis()->SetRangeUser(-4, -2.5);
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->GetYaxis()->SetTitle("Entries");
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->GetYaxis()->SetTitleOffset(1.2);
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->GetXaxis()->SetTitleOffset(1.15);

        if (i_LF_Generator == 0)
        {

            C_Y_DiMuRec->cd();
            h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->Draw("PE");

            C_Pt_DiMuRec->cd();
            h_Pt_DiMuon_Rec_fromLF[i_LF_Generator]->Draw("PE");

            C_M_DiMuRec->cd();
            h_M_DiMuon_Rec_fromLF[i_LF_Generator]->Draw("PE");
        }
        else
        {
            C_Y_DiMuRec->cd();
            h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->Draw("PESAME");

            C_Pt_DiMuRec->cd();
            h_Pt_DiMuon_Rec_fromLF[i_LF_Generator]->Draw("PESAME");

            C_M_DiMuRec->cd();
            h_M_DiMuon_Rec_fromLF[i_LF_Generator]->Draw("PESAME");
        }
    }

    TLegend *leg = C_Y_DiMuRec->BuildLegend(0.7, 0.75, 0.9, 0.95, " ");
    leg->SetFillStyle(0);
    leg->SetLineColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.0225);
    C_Pt_DiMuRec->cd();
    leg->Draw();
    C_M_DiMuRec->cd();
    leg->Draw();

    delete[] DiMuon_fromLF_Generator;
    delete[] DiMuon_fromLF_Generator_Title;
    delete[] DiMuon_fromLF_Generator_Color;
    delete[] h_Y_DiMuon_Rec_fromLF;
    delete[] h_PtM_DiMuon_Rec_fromLF;
    delete[] h_PtY_DiMuon_Rec_fromLF;
    delete[] h_Pt_DiMuon_Rec_fromLF;
    delete[] h_M_DiMuon_Rec_fromLF;
    // C_Y->SaveAs("images/LF_DiMuRec_Y.png");
}

void LF_Hadron(TString gen)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // std::cout<<Nev_PYTHIA_STANDALONE_NoDiffr<<endl;
    // fIn_PYTHIA_CENTRALIZED_MB->cd("Hadron/Geant");
    // fIn_PYTHIA_CENTRALIZED_MB->ls();
    // fIn_PYTHIA_CENTRALIZED_MB->cd("Hadron/Pythia");
    // fIn_PYTHIA_CENTRALIZED_MB->ls();

    Int_t Rebin_Pt = 10;
    Int_t Rebin_Eta = 5;
    TString *DiMuon_fromLF_Generator;
    TString *DiMuon_fromLF_Generator_Title;
    Color_t *DiMuon_fromLF_Generator_Color;
    // return;
    const Int_t n_PDG_tested = 1;

    // TString canvas_header[n_PDG_tested] = {"#pi^{#plus}", "K^{0}", "K^{#plus}"};
    // TString canvas_name[n_PDG_tested] = {"pi_plus", "kzero", "kplus"};
    Int_t low_PDG[n_PDG_tested] = {210};
    Int_t high_PDG[n_PDG_tested] = {220};

    Int_t nGenerator = 999;
    if (gen.Contains("Pythia"))
    {
        nGenerator = 2;
        DiMuon_fromLF_Generator = new TString[nGenerator]{"PYTHIA", "Geant"};
        DiMuon_fromLF_Generator_Color = new Color_t[nGenerator]{kMagenta + 2, kGreen + 2};
    }
    else if (gen.Contains("Powheg"))
    {
        nGenerator = 3;
        DiMuon_fromLF_Generator = new TString[nGenerator]{"PYTHIA", "Geant", "Powheg"};
        DiMuon_fromLF_Generator_Color = new Color_t[nGenerator]{kMagenta + 2, kGreen + 2, kOrange + 2};
    }

    TH2F **h_PdgEta_HFHadron_prompt = new TH2F *[nGenerator];
    TH2F **h_PdgPt_HFHadron_prompt = new TH2F *[nGenerator];
    TH1F **h_Eta_HFHadron_prompt = new TH1F *[nGenerator];
    TH1F **h_Pt_HFHadron_prompt = new TH1F *[nGenerator];

    TH1F **h_Ratio_Eta_HFHadron_prompt = new TH1F *[nGenerator - 1];
    TH1F **h_Ratio_Pt_HFHadron_prompt = new TH1F *[nGenerator - 1];

    for (Int_t i_Generator = 0; i_Generator < nGenerator; i_Generator++)
    {
        if (gen.Contains("Pythia"))
        {
            h_PdgEta_HFHadron_prompt[i_Generator] = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get(Form("Hadron/%s/h_PdgEta_HFHadron_prompt_%sOnly", DiMuon_fromLF_Generator[i_Generator].Data(), DiMuon_fromLF_Generator[i_Generator].Data()));
            h_PdgPt_HFHadron_prompt[i_Generator] = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get(Form("Hadron/%s/h_PdgPt_HFHadron_prompt_%sOnly", DiMuon_fromLF_Generator[i_Generator].Data(), DiMuon_fromLF_Generator[i_Generator].Data()));
        }
        else if (gen.Contains("Powheg"))
        {
            cout << Form("Hadron/%s/h_PdgEta_HFHadron_prompt_%sOnly", DiMuon_fromLF_Generator[i_Generator].Data(), DiMuon_fromLF_Generator[i_Generator].Data()) << endl;
            h_PdgEta_HFHadron_prompt[i_Generator] = (TH2F *)fIn_Powheg_Charm->Get(Form("Hadron/%s/h_PdgEta_HFHadron_prompt_%sOnly", DiMuon_fromLF_Generator[i_Generator].Data(), DiMuon_fromLF_Generator[i_Generator].Data()));
            h_PdgPt_HFHadron_prompt[i_Generator] = (TH2F *)fIn_Powheg_Charm->Get(Form("Hadron/%s/h_PdgPt_HFHadron_prompt_%sOnly", DiMuon_fromLF_Generator[i_Generator].Data(), DiMuon_fromLF_Generator[i_Generator].Data()));
        }
        h_PdgEta_HFHadron_prompt[i_Generator]->GetXaxis()->SetRangeUser(low_PDG[0], high_PDG[0]);
        h_PdgPt_HFHadron_prompt[i_Generator]->GetXaxis()->SetRangeUser(low_PDG[0], high_PDG[0]);

        h_Eta_HFHadron_prompt[i_Generator] = (TH1F *)h_PdgEta_HFHadron_prompt[i_Generator]->ProjectionY();
        hist1D_graphic_opt(h_Eta_HFHadron_prompt[i_Generator], kFALSE, Rebin_Eta, 20, DiMuon_fromLF_Generator_Color[i_Generator], 1.);
        h_Eta_HFHadron_prompt[i_Generator]->SetTitle(Form("from %s", DiMuon_fromLF_Generator[i_Generator].Data()));
        h_Eta_HFHadron_prompt[i_Generator]->GetYaxis()->SetTitle("Entries");

        h_Pt_HFHadron_prompt[i_Generator] = (TH1F *)h_PdgPt_HFHadron_prompt[i_Generator]->ProjectionY();
        hist1D_graphic_opt(h_Pt_HFHadron_prompt[i_Generator], kFALSE, Rebin_Pt, 20, DiMuon_fromLF_Generator_Color[i_Generator], 1.);
        h_Pt_HFHadron_prompt[i_Generator]->SetTitle(Form("from %s", DiMuon_fromLF_Generator[i_Generator].Data()));
        h_Pt_HFHadron_prompt[i_Generator]->GetYaxis()->SetTitle("Entries");

        if (i_Generator > 0)
        {
            h_Ratio_Eta_HFHadron_prompt[i_Generator - 1] = (TH1F *)h_Eta_HFHadron_prompt[i_Generator]->Clone(Form("h_Ratio_Eta_HFHadron_prompt_%s_%s", DiMuon_fromLF_Generator[i_Generator].Data(), DiMuon_fromLF_Generator[0].Data()));
            h_Ratio_Eta_HFHadron_prompt[i_Generator - 1]->SetTitle(Form("%s/PYTHIA", DiMuon_fromLF_Generator[i_Generator].Data()));
            h_Ratio_Eta_HFHadron_prompt[i_Generator - 1]->Divide(h_Eta_HFHadron_prompt[0]);
            h_Ratio_Pt_HFHadron_prompt[i_Generator - 1] = (TH1F *)h_Pt_HFHadron_prompt[i_Generator]->Clone(Form("h_Ratio_Pt_HFHadron_prompt_%s_%s", DiMuon_fromLF_Generator[i_Generator].Data(), DiMuon_fromLF_Generator[0].Data()));
            h_Ratio_Eta_HFHadron_prompt[i_Generator - 1]->SetTitle(Form("%s/PYTHIA", DiMuon_fromLF_Generator[i_Generator].Data()));
            h_Ratio_Pt_HFHadron_prompt[i_Generator - 1]->Divide(h_Pt_HFHadron_prompt[0]);
        }
    }

    TCanvas *C_Eta = three_histo_ratio(h_Eta_HFHadron_prompt[0], h_Eta_HFHadron_prompt[1], h_Eta_HFHadron_prompt[2], h_Ratio_Eta_HFHadron_prompt[0], h_Ratio_Eta_HFHadron_prompt[1], "C_Eta", "#splitline{POWHEG+PYTHIA6 pp @13 Tev}{Generated #pi}", kTRUE);
    // C_Eta->SaveAs("images/Eta_pions_Powheg_Pythia.png");

    TCanvas *C_Pt = three_histo_ratio(h_Pt_HFHadron_prompt[0], h_Pt_HFHadron_prompt[1], h_Pt_HFHadron_prompt[2], h_Ratio_Pt_HFHadron_prompt[0], h_Ratio_Pt_HFHadron_prompt[1], "C_Pt", "#splitline{POWHEG+PYTHIA6 pp @13 Tev}{Generated #pi}", kTRUE);
    // C_Pt->SaveAs("images/Pt_pions_Powheg_Pythia.png");
}