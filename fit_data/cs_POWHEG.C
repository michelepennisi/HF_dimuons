#include "/home/michele_pennisi/cernbox/common_include.h"

void cs_POWHEG()
{
    const Int_t n_DiMuSelection = 2;
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");
    TString mass_range;
    mass_range.Form("_LowMass_LowPt");
    // TFile *fIn = new TFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/Powheg_pdfMC_unbinned%s.root", mass_range.Data()), "READ");
    TFile *fIn = new TFile("rf607_fitresult.root", "READ");

    RooFitResult *r = (RooFitResult *)fIn->Get("rf607");
    const RooArgList &fitParams = r->floatParsFinal();
    Double_t fit_result[n_DiMuSelection];
    for (int i = 0; i < fitParams.getSize(); ++i)
    {
        auto &fitPar = (RooRealVar &)fitParams[i];
        std::cout << fitPar.GetName() << " " << fitPar.getVal() << std::endl;
        fit_result[i] = (Double_t)fitPar.getVal();
    }
    TFile *fIn_MC[2];
    fIn_MC[0] = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version3_AliAOD/save_output/LHC23i2_MC_output_Tree_merged.root", "READ");
    fIn_MC[1] = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i1/Version3_AliAOD/save_output/LHC23i1_MC_output_Tree_merged.root", "READ");

    TString name_2DiMuSelection[n_DiMuSelection];
    name_2DiMuSelection[0].Form("Beauty");
    name_2DiMuSelection[1].Form("Charm");

    TTree *m_tree_MC[n_DiMuSelection];
    Double_t MC_dimu[n_DiMuSelection];
    Double_t n_ev_MC[n_DiMuSelection];
    TH1F *NEv_MC[n_DiMuSelection];
    Double_t feq_MC[n_DiMuSelection] = {15 * (56.42 / 0.5), 30.5 * (56.42 / 5)};

    TFile *fIn_data = new TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/HistResults_merged.root", "READ");
    TH1D *fhNEv = (TH1D *)fIn_data->Get("fhNEv");

    Double_t c_frac_fit_MB[n_DiMuSelection];
    Double_t c_frac_MC_MB[n_DiMuSelection];
    TFile *fIn_Powheg[n_DiMuSelection];
    TH1F *h_NHF_event_PowhegOnly[n_DiMuSelection];
    TH1F *h_NHF_event_PowhegOnly_fwd[n_DiMuSelection];
    Double_t N_Pair[n_DiMuSelection] = {0, 0};
    Double_t N_Pair_fwd[n_DiMuSelection] = {0, 0};
    Double_t Powheg_CS[n_DiMuSelection] = {0.52832510, 5.};
    Double_t CS[n_DiMuSelection] = {999., 999.};
    Double_t MPI_corr[n_DiMuSelection] = {999., 999.};

    fIn_Powheg[0] = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_beauty_nocut_Version_5_AliAOD_withHF_Q/powheg_beauty_nocut_Version_5_AliAOD_withHF_Q_MC_output_Hist_294154.root", "READ");

    h_NHF_event_PowhegOnly[0] = (TH1F *)fIn_Powheg[0]->Get("HF_quarks/Powheg/h_NBeauty_event_PowhegOnly");

    h_NHF_event_PowhegOnly_fwd[0] = (TH1F *)fIn_Powheg[0]->Get("HF_quarks/Powheg/h_NBeauty_event_fwd_PowhegOnly");

    fIn_Powheg[1] = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_charm_nocut_Version_5_AliAOD_withHF_Q/powheg_charm_nocut_Version_5_AliAOD_withHF_Q_MC_output_Hist_294154.root", "READ");

    h_NHF_event_PowhegOnly[1] = (TH1F *)fIn_Powheg[1]->Get("HF_quarks/Powheg/h_NCharm_event_PowhegOnly");

    h_NHF_event_PowhegOnly_fwd[1] = (TH1F *)fIn_Powheg[1]->Get("HF_quarks/Powheg/h_NCharm_event_fwd_PowhegOnly");

    TFile *fIn_Pythia = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_inel_LFoff_Def_pythia_sim_2411_DefaultBR_output_Hist_100000.root", "READ");

    TH1F *h_NHF_event_Pythia[n_DiMuSelection] = {(TH1F *)fIn_Pythia->Get("HF_quarks/h_NBeauty_event"), (TH1F *)fIn_Pythia->Get("HF_quarks/h_NCharm_event")};
    Double_t n_HF_single_pair_PYTHIA[n_DiMuSelection] = {(h_NHF_event_Pythia[0]->GetBinContent(3)), (h_NHF_event_Pythia[1]->GetBinContent(3))};
    Double_t n_HF_total_PYTHIA[n_DiMuSelection] = {0., 0.};

    // for (Int_t i_bin = 1; i_bin < h_NHF_event_Pythia[0]->GetNbinsX(); i_bin++)
    // {
    //     n_HF_total_PYTHIA[0] = n_HF_total_PYTHIA[0] + (0.5) * h_NHF_event_Pythia[0]->GetBinContent(i_bin) * (i_bin - 1);
    // }

    // cout << (double)n_HF_total_PYTHIA[0]/n_HF_single_pair_PYTHIA[0] << endl;
    // return;
    for (Int_t i = 0; i < n_DiMuSelection; i++)
    {
        for (Int_t i_bin = 1; i_bin < h_NHF_event_PowhegOnly[i]->GetNbinsX(); i_bin++)
        {
            // printf("Bin %d || Content %0.0f\n", i_bin, h_NHF_event_PowhegOnly[i]->GetBinContent(i_bin));
            N_Pair[i] = N_Pair[i] + (0.5) * h_NHF_event_PowhegOnly[i]->GetBinContent(i_bin) * (i_bin - 1);
            N_Pair_fwd[i] = N_Pair_fwd[i] + (0.5) * h_NHF_event_PowhegOnly_fwd[i]->GetBinContent(i_bin) * (i_bin - 1);

            n_HF_total_PYTHIA[i] = n_HF_total_PYTHIA[i] + (0.5) * h_NHF_event_Pythia[i]->GetBinContent(i_bin) * (i_bin - 1);
        }

        cout << "Result for " << name_2DiMuSelection[i].Data() << endl;
        NEv_MC[i] = (TH1F *)fIn_MC[i]->Get("h_Nevents");
        n_ev_MC[i] = feq_MC[i] * NEv_MC[i]->GetBinContent(2);
        printf("MC: Nev %0.3e || f norm powheg %0.3e || Nev MB %0.3e\n", NEv_MC[i]->GetBinContent(2), feq_MC[i], n_ev_MC[i]);
        m_tree_MC[i] = (TTree *)fIn_MC[i]->Get(Form("DiMuon_Rec%s_PowhegOnly%s", mass_range.Data(), name_2DiMuSelection[i].Data()));
        printf("Data: Nev = %0.3e || Nev MB (2384.73 * fhNEv->GetBinContent(3)) %0.3e\n", fhNEv->GetBinContent(3), 2384.73 * fhNEv->GetBinContent(3));
        c_frac_fit_MB[i] = (Double_t)fit_result[i] / (2384.73 * fhNEv->GetBinContent(3));
        c_frac_MC_MB[i] = (Double_t)m_tree_MC[i]->GetEntries() / n_ev_MC[i];
        printf("Result\n From fit %0.5f\n norm MB event %0.5e\n",   fit_result[i], c_frac_fit_MB[i]);
        printf("from MC %0.0lld\n norm MB event %0.5e\n", m_tree_MC[i]->GetEntries(), c_frac_MC_MB[i]);
        Double_t Ratio = c_frac_fit_MB[i] / c_frac_MC_MB[i];
        printf("Ratio %0.3e\n", Ratio);
        printf("Fwd fraction %0.3e || (N_Pair_fwd %f  N_Pair_total %f )", (Double_t)N_Pair_fwd[i] / N_Pair[i], N_Pair_fwd[i], N_Pair[i]);
        CS[i] = (Powheg_CS[i] * (Double_t)N_Pair_fwd[i] / N_Pair[i]) / (2 * 1.5);
        printf("CS from MC %0.3e\n", CS[i]);
        printf("Measured CS %0.3f\n", (Ratio * CS[i]));
        printf("MPI correction with PYTHIA %0.3f\n", n_HF_total_PYTHIA[i] / n_HF_single_pair_PYTHIA[i]);
        printf("Measured CS corrected for MPI %0.3e\n", (Ratio * CS[i] * n_HF_total_PYTHIA[i] / n_HF_single_pair_PYTHIA[i]));
        cout << "=================================================" << endl;
    }
}