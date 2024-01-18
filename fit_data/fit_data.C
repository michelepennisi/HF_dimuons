#include "/home/michele_pennisi/cernbox/common_include.h"
#include <vector>
TCanvas *printMC_ratio(TString name, RooPlot *frame, TH1 *data, TF1 *pdf, Color_t color, Int_t minx = 0, Int_t max_x = 30);
TCanvas *printRooPlot_ratio(RooPlot *frame, Bool_t norm, RooFitResult *r, Int_t choice, TString roohist_name, TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x);
Double_t Lowy = 0.0025;

void unbinned_fit_data_sample_singleregion(Int_t Low_Mass = 4, Int_t High_Mass = 9, Int_t Low_Pt = 0, Int_t High_Pt = 10, Bool_t withDY = kTRUE);
void fit_data()
{
    unbinned_fit_data_sample_singleregion(4, 30, 0, 30);
}
void cross_section(Int_t Low_Mass = 11, Int_t High_Mass = 30, Int_t Low_Pt = 0, Int_t High_Pt = 30, Bool_t withDY = kTRUE)
{
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");

    TString Generator = "Powheg+PYTHIA6";
    TFile *fIn = TFile::Open("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/Powheg_pdfMC_unbinnedprova.root", "READ");

    RooFitResult *r;

    if (withDY)
        r = (RooFitResult *)fIn->Get((Form("fit_result_M_%d_%d_Pt_%d_%d_withDY", Low_Mass, High_Mass, Low_Pt, High_Pt)));
    else
        r = (RooFitResult *)fIn->Get((Form("fit_result_M_%d_%d_Pt_%d_%d_noDY", Low_Mass, High_Mass, Low_Pt, High_Pt)));

    const RooArgList &fitParams = r->floatParsFinal();
    int n_DiMuSelection = fitParams.getSize();

    std::vector<Double_t> fit_result;
    TString *MC_file = new TString;
    TFile *fIn_MC[n_DiMuSelection];
    MC_file[0] = TString("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i1/Version3_AliAOD/save_output/LHC23i1_MC_output_Tree_merged.root");
    MC_file[1] = TString("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version3_AliAOD/save_output/LHC23i2_MC_output_Tree_merged.root");
    MC_file[2] = TString("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Powheg_Sim/powheg_DY_mass_3_35/Version1/Analysis_MCsim/powheg_DY_mass_3_35_Analysis_MCsim_Tree_merged.root");

    for (Int_t i = 0; i < fitParams.getSize(); i++)
    {
        fIn_MC[i] = new TFile(MC_file[i].Data(), "READ");
        auto &fitPar = (RooRealVar &)fitParams[fitParams.getSize() - i - 1];
        fit_result.push_back(fitPar.getVal());
        std::cout << fitPar.GetName() << " " << fit_result[i] << std::endl;
    }

    TString *name_DiMu_Sel = new TString;
    name_DiMu_Sel[0] = TString("Charm");
    name_DiMu_Sel[1] = TString("Beauty");
    name_DiMu_Sel[2] = TString("DY");

    TTree **m_tree_MC = new TTree *;
    TTree **m_tree_MC_cutted = new TTree *;
    std::vector<Double_t> MC_dimu;
    std::vector<Double_t> n_ev_MC;
    TH1F **NEv_MC = new TH1F *;
    std::vector<Double_t> feq_MC;
    feq_MC.push_back(30.5 * (56.42 / 5.0));
    feq_MC.push_back(15.0 * (56.42 / 0.5));
    feq_MC.push_back(1. * (56.42 / 4.6e-05));

    TFile *fIn_data = new TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/HistResults_merged.root", "READ");
    TH1D *fhNEv = (TH1D *)fIn_data->Get("fhNEv");

    std::vector<Double_t> c_frac_fit_MB;
    std::vector<Double_t> c_frac_MC_MB;
    TFile **fIn_Powheg = new TFile *;
    fIn_Powheg[0] = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_charm_nocut_Version_5_AliAOD_withHF_Q/powheg_charm_nocut_Version_5_AliAOD_withHF_Q_MC_output_Hist_294154.root", "READ");
    fIn_Powheg[1] = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/powheg_beauty_nocut_Version_5_AliAOD_withHF_Q/powheg_beauty_nocut_Version_5_AliAOD_withHF_Q_MC_output_Hist_294154.root", "READ");
    fIn_Powheg[2] = new TFile("~/cernbox/HF_dimuons/fit_data/ingredient_cs_powheg/DY_cs.root", "READ");

    TH1F **h_NHF_event_PowhegOnly = new TH1F *;
    TH1F **h_NHF_event_PowhegOnly_fwd = new TH1F *;
    h_NHF_event_PowhegOnly[0] = (TH1F *)fIn_Powheg[0]->Get("HF_quarks/Powheg/h_NCharm_event_PowhegOnly");
    h_NHF_event_PowhegOnly_fwd[0] = (TH1F *)fIn_Powheg[0]->Get("HF_quarks/Powheg/h_NCharm_event_fwd_PowhegOnly");

    h_NHF_event_PowhegOnly[1] = (TH1F *)fIn_Powheg[1]->Get("HF_quarks/Powheg/h_NBeauty_event_PowhegOnly");
    h_NHF_event_PowhegOnly_fwd[1] = (TH1F *)fIn_Powheg[1]->Get("HF_quarks/Powheg/h_NBeauty_event_fwd_PowhegOnly");

    h_NHF_event_PowhegOnly[2] = (TH1F *)fIn_Powheg[2]->Get("h_NDY_event");
    h_NHF_event_PowhegOnly_fwd[2] = (TH1F *)fIn_Powheg[2]->Get("h_NDY_event_fwd");
    std::vector<Double_t> Powheg_CS;
    Powheg_CS.push_back(5.0);
    Powheg_CS.push_back(0.50);
    Powheg_CS.push_back(4.6e-05);
    std::vector<Double_t> N_Pair;
    std::vector<Double_t> N_Pair_fwd;
    std::vector<Double_t> CS;
    std::vector<Double_t> MPI_corr;
    TFile *fIn_Pythia = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_inel_LFoff_Def_pythia_sim_2411_DefaultBR_output_Hist_100000.root", "READ");

    TH1F **h_NHF_event_Pythia = new TH1F *;
    h_NHF_event_Pythia[0] = (TH1F *)fIn_Pythia->Get("HF_quarks/h_NCharm_event");
    h_NHF_event_Pythia[1] = (TH1F *)fIn_Pythia->Get("HF_quarks/h_NBeauty_event");

    std::vector<Double_t> n_HF_single_pair_PYTHIA;
    n_HF_single_pair_PYTHIA.push_back(h_NHF_event_Pythia[0]->GetBinContent(3));
    n_HF_single_pair_PYTHIA.push_back(h_NHF_event_Pythia[1]->GetBinContent(3));

    std::vector<Double_t> n_HF_total_PYTHIA;
    N_Pair.push_back(0.);
    N_Pair_fwd.push_back(0.);
    n_HF_total_PYTHIA.push_back(0.);

    for (Int_t i = 0; i < n_DiMuSelection; i++)
    {

        cout << N_Pair[i] << endl;
        cout << "Result for " << name_DiMu_Sel[i].Data() << endl;

        NEv_MC[i] = (TH1F *)fIn_MC[i]->Get("h_Nevents");
        printf("ARRIVO\n");
        for (Int_t i_bin = 1; i_bin < h_NHF_event_PowhegOnly[i]->GetNbinsX(); i_bin++)
        {
            N_Pair[i] = (N_Pair[i] + (0.5) * h_NHF_event_PowhegOnly[i]->GetBinContent(i_bin) * (i_bin - 1));
            // printf("Bin %d || Content %0.0f\n", i_bin, h_NHF_event_PowhegOnly[i]->GetBinContent(i_bin));
            N_Pair_fwd[i] = (N_Pair_fwd[i] + (0.5) * h_NHF_event_PowhegOnly_fwd[i]->GetBinContent(i_bin) * (i_bin - 1));
            if (i < n_DiMuSelection - 1)
                n_HF_total_PYTHIA[i] = (n_HF_total_PYTHIA[i] + (0.5) * h_NHF_event_Pythia[i]->GetBinContent(i_bin) * (i_bin - 1));
            else
                n_HF_total_PYTHIA[i] = 0;
        }

        fIn_MC[i]->cd();

        if (i < n_DiMuSelection - 1)
        {
            n_ev_MC.push_back(feq_MC[i] * NEv_MC[i]->GetBinContent(2));
            m_tree_MC[i] = (TTree *)fIn_MC[i]->Get(Form("DiMuon_Rec_PowhegOnly_%s", name_DiMu_Sel[i].Data()));
            printf("MC: Nev %0.3e || f norm powheg %0.3e || Nev MB %0.3e\n", NEv_MC[i]->GetBinContent(2), feq_MC[i], n_ev_MC[i]);
            CS.push_back((Powheg_CS[i] * N_Pair_fwd[i] / N_Pair[i]) / (2 * 1.5));
        }
        else
        {

            m_tree_MC[i] = (TTree *)fIn_MC[i]->Get("rec_tree_muDY");
            n_ev_MC.push_back(feq_MC[i] * 637848.);
            printf("MC: Nev %0.3e || f norm powheg %0.3e || Nev MB %0.3e\n", 637848., feq_MC[i], n_ev_MC[i]);
            CS.push_back((Powheg_CS[i] * N_Pair_fwd[i] / N_Pair[i]) / (1.5));
        }
        gROOT->cd();
        printf("Fwd fraction %0.3e || (N_Pair_fwd %f  N_Pair_total %f )", N_Pair_fwd[i] / N_Pair[i], N_Pair_fwd[i], N_Pair[i]);
        printf("CS from MC %0.3e\n", CS[i]);

        m_tree_MC_cutted[i] = (TTree *)m_tree_MC[i]->CopyTree(Form("(m>%d && m<%d) && (pt > %d && pt <%d)", Low_Mass, High_Mass, Low_Pt, High_Pt));
        m_tree_MC_cutted[i]->SetName(Form("Tree_%s_M_%d_%d_Pt_%d_%d", name_DiMu_Sel[i].Data(), Low_Mass, High_Mass, Low_Pt, High_Pt));

        printf("Data: Nev = %0.3e || Nev MB (2384.73 * fhNEv->GetBinContent(3)) %0.3e\n", fhNEv->GetBinContent(3), 2384.73 * fhNEv->GetBinContent(3));
        c_frac_fit_MB.push_back(fit_result[i] / (2384.73 * fhNEv->GetBinContent(3)));
        c_frac_MC_MB.push_back(m_tree_MC_cutted[i]->GetEntries() / n_ev_MC[i]);
        printf("Result\n From fit %0.5f\n norm MB event %0.5e\n", fit_result[i], c_frac_fit_MB[i]);
        printf("from MC %0.0lld\n norm MB event %0.5e\n", m_tree_MC_cutted[i]->GetEntries(), c_frac_MC_MB[i]);
        Double_t Ratio = c_frac_fit_MB[i] / c_frac_MC_MB[i];
        printf("Ratio %0.3e\n", Ratio);

        printf("Measured CS %0.3e\n", (Ratio * CS[i]));
        printf("MPI correction with PYTHIA %0.3f\n", n_HF_total_PYTHIA[i] / n_HF_single_pair_PYTHIA[i]);
        printf("Measured CS corrected for MPI %0.3e\n", (Ratio * CS[i] * n_HF_total_PYTHIA[i] / n_HF_single_pair_PYTHIA[i]));
        cout << "=================================================" << endl;
    }
}

void unbinned_fit_data_sample_singleregion(Int_t Low_Mass = 4, Int_t High_Mass = 9, Int_t Low_Pt = 0, Int_t High_Pt = 10, Bool_t withDY = kTRUE)
{

    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/high_mass_dimuons/fit_library/PtMassPol1ExpPdf.cxx+");

    Int_t Mass_Binning = (High_Mass - Low_Mass) / 0.5;
    Int_t Pt_Binning = (High_Pt - Low_Pt) / 0.5;
    const Int_t n_DiMuSelection = 4;
    TString Name_DimuSel[n_DiMuSelection] = {"Charm", "Beauty", "HF_Mixed", "DY"};
    Color_t color[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kAzure + 9, kOrange + 7};
    TFile *fIn = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/Powheg_pdfMC_unbinnedprova.root", "READ");

    RooWorkspace *w = (RooWorkspace *)fIn->Get(Form("w_M_%d_%d_Pt_%d_%d", Low_Mass, High_Mass, Low_Pt, High_Pt));
    w->Print("s");
    fIn->Close();

    // RooRealVar *m = new RooRealVar("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", Low_Mass, High_Mass);
    // RooRealVar *pt = new RooRealVar("pt", "#it{p}_{T} (GeV/#it{c})", Low_Pt, High_Pt);

    RooRealVar *m = w->var("m");
    RooRealVar *pt = w->var("pt");
    m->setBins(Mass_Binning);
    pt->setBins(Pt_Binning);

    RooCategory sample("sample", "sample");
    sample.defineType("mass");
    sample.defineType("transversemomentum");
    TFile *fIn_data = new TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/TreeResults_merged.root", "READ");
    fIn_data->ls();

    // Taking data saved in tree
    TTree *tree_data = (TTree *)fIn_data->Get("rec_data_tree");
    // tree_data->Draw("m", Form("((m>%d && m<9) || (m>11 && m<%d)) && (pt > %d && pt <%d)", Low_Mass, High_Mass, Low_Pt, High_Pt));
    gROOT->cd();

    TTree *tree_data_cutted;
    if (Low_Mass == 4 && High_Mass == 30)
        // tree_data_cutted = (TTree *)tree_data->CopyTree(Form("((m>%d && m<9) || (m>11 && m<%d)) && (pt > %d && pt <%d)", Low_Mass, High_Mass, Low_Pt, High_Pt));
        tree_data_cutted = (TTree *)fIn_data->Get("rec_data_tree_withcut");
    else
        tree_data_cutted = (TTree *)tree_data->CopyTree(Form("(m>%d && m<%d) && (pt > %d && pt <%d)", Low_Mass, High_Mass, Low_Pt, High_Pt));
    tree_data_cutted->SetName("tree_data_cutted");
    RooDataSet *unbinned_M_Dimu_data = new RooDataSet("M_Dimu_data", "M_Dimu_data", RooArgSet(*m), Import(*tree_data_cutted), Cut("m<9 || (m>11 && m<30)"));
    RooDataSet *unbinned_Pt_Dimu_data = new RooDataSet("Pt_Dimu_data", "Pt_Dimu_data", RooArgSet(*pt), Import(*tree_data_cutted));
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
        B_DimuMass[i_DiMu_Sel] = w->var(Form("B_DimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        B_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        n1_DimuMass[i_DiMu_Sel] = w->var(Form("n1_DimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        n1_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        n2_DimuMass[i_DiMu_Sel] = w->var(Form("n2_DimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        n2_DimuMass[i_DiMu_Sel]->setConstant(kTRUE);
        B_DimuPt[i_DiMu_Sel] = w->var(Form("B_DimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        B_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        n1_DimuPt[i_DiMu_Sel] = w->var(Form("n1_DimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        n1_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);
        n2_DimuPt[i_DiMu_Sel] = w->var(Form("n2_DimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        n2_DimuPt[i_DiMu_Sel]->setConstant(kTRUE);

        pdfDimuMass[i_DiMu_Sel] = w->pdf(Form("pdfDimuMassFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
        pdfDimuPt[i_DiMu_Sel] = w->pdf(Form("pdfDimuPtFrom%s", Name_DimuSel[i_DiMu_Sel].Data()));
    }

    RooRealVar *normForC = new RooRealVar("n_charm_output", "number dimuon from c", 28440, 0, 200000);
    RooRealVar *normForB = new RooRealVar("n_beauty_output", "number dimuon from b", 48000, 0, 200000);
    RooRealVar *normForDY = new RooRealVar("n_DY_output", "number dimuon from DY", 48000, 0, 200000);
    RooRealVar *normForMixed = new RooRealVar("n_mixed_output", "number dimuon from b,c", (3.6 / 100) * tree_data_cutted->GetEntries());
    normForMixed->setConstant(kTRUE);

    RooAddPdf *m_model;
    RooAddPdf *pt_model;
    if (withDY)
    {
        m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*pdfDimuMass[0], *pdfDimuMass[1], *pdfDimuMass[2], *pdfDimuMass[3]), RooArgList(*normForC, *normForB, *normForMixed, *normForDY));
        pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*pdfDimuPt[0], *pdfDimuPt[1], *pdfDimuPt[2], *pdfDimuPt[3]), RooArgList(*normForC, *normForB, *normForMixed, *normForDY));
    }
    else
    {
        m_model = new RooAddPdf("m_model", "n_charm_output*dimuMassFromC + n_beauty_output*dimuMassFromB + n_mixed_output*dimuMassFromMixed", RooArgList(*pdfDimuMass[0], *pdfDimuMass[1], *pdfDimuMass[2]), RooArgList(*normForC, *normForB, *normForMixed));
        pt_model = new RooAddPdf("pt_model", "n_charm_output*dimuPtFromC + n_beauty_output*dimuPtFromB + n_mixed_output*dimuPtFromMixed", RooArgList(*pdfDimuPt[0], *pdfDimuPt[1], *pdfDimuPt[2]), RooArgList(*normForC, *normForB, *normForMixed));
    }

    RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
    simPdf.addPdf(*m_model, "mass");
    simPdf.addPdf(*pt_model, "transversemomentum");
    simPdf.Print("t");
    RooFitResult *r;
    if (Low_Mass == 4 && High_Mass == 30)
    {
        m->setRange("pluto", 4, 9);
        pt->setRange("pluto", 20, 30);
        m->setRange("pippo", 11, 30);
        pt->setRange("pippo", 0, 20);
        r = simPdf.fitTo(*unbinned_combData_set, Range("pluto,pippo"), SumCoefRange("pluto,pippo"), SplitRange(false), Save(), SumW2Error(true));
    }
    else
        r = simPdf.fitTo(*unbinned_combData_set, Minimizer("Minuit2"), Save(), SumW2Error(true));
    fIn = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/Powheg_pdfMC_unbinnedprova.root", "UPDATE");
    if (withDY)
    {
        r->Write(Form("fit_result_M_%d_%d_Pt_%d_%d_withDY", Low_Mass, High_Mass, Low_Pt, High_Pt));
    }
    else
    {
        r->Write(Form("fit_result_M_%d_%d_Pt_%d_%d_noDY", Low_Mass, High_Mass, Low_Pt, High_Pt));
    }

    fIn->Close();

    RooPlot *m_frame = m->frame(Title("m_frame"));
    m_frame->SetTitle(" ");
    m_frame->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
    m_frame->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");

    RooPlot *pt_frame = pt->frame(Title("pt_frame"));
    pt_frame->SetTitle(" ");
    pt_frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    pt_frame->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

    unbinned_combData_set->plotOn(m_frame, Name("combDatamass"), Cut("sample==sample::mass"), DrawOption("PEZ"));
    unbinned_combData_set->plotOn(pt_frame, Name("combDatapt"), Cut("sample==sample::transversemomentum"), DrawOption("PEZ"));
    if (Low_Mass == 4 && High_Mass == 30)
    {
        simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *unbinned_combData_set), Range("pluto,pippo"), NormRange("pluto,pippo"), LineStyle(kSolid), LineColor(kRed));
        simPdf.plotOn(pt_frame, Name("pdfpt"), Slice(sample, "transversemomentum"), ProjWData(sample, *unbinned_combData_set), Range("pluto,pippo"), NormRange("pluto,pippo"), LineStyle(kSolid), LineColor(kRed));
    }
    else
    {
        simPdf.plotOn(m_frame, Name("pdfmass"), Slice(sample, "mass"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
        simPdf.plotOn(pt_frame, Name("pdfpt"), Slice(sample, "transversemomentum"), ProjWData(sample, *unbinned_combData_set), LineStyle(kSolid), LineColor(kRed));
    }

    r->floatParsFinal().Print("s");
    for (Int_t i_DiMu_Sel = 0; i_DiMu_Sel < n_DiMuSelection; i_DiMu_Sel++)
    {
        if (Low_Mass == 4 && High_Mass == 30)
        {
            simPdf.plotOn(m_frame, Name(Form("pdfmass%s", Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "mass"), Components(pdfDimuMass[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), Range("pluto,pippo"), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(color[i_DiMu_Sel]), LineWidth(5));
            simPdf.plotOn(pt_frame, Name(Form("pdfpt%s", Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "transversemomentum"), Components(pdfDimuPt[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), Range("pluto,pippo"), NormRange("pluto,pippo"), LineStyle(kDashed), LineColor(color[i_DiMu_Sel]), LineWidth(5));
        }
        else
        {
            simPdf.plotOn(m_frame, Name(Form("pdfmass%s", Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "mass"), Components(pdfDimuMass[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(color[i_DiMu_Sel]), LineWidth(5));
            simPdf.plotOn(pt_frame, Name(Form("pdfpt%s", Name_DimuSel[i_DiMu_Sel].Data())), Slice(sample, "transversemomentum"), Components(pdfDimuPt[i_DiMu_Sel]->GetName()), ProjWData(sample, *unbinned_combData_set), LineStyle(kDashed), LineColor(color[i_DiMu_Sel]), LineWidth(5));
        }
    }
    // Save the roodataset as histogram and RooSimultaneous as TF1

    RooArgSet *m_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*m_model).snapshot(true)); // True means copy the PDF and everything it depends on
    auto &m_modelcopiedPdf = static_cast<RooAbsPdf &>((*m_modelcopyOfEverything)["m_model"]);          // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    RooArgSet *m_modelobs = m_modelcopiedPdf.getObservables(*unbinned_combData_set);
    RooArgSet *m_modelPars = m_modelcopiedPdf.getParameters(*m_modelobs);
    TF1 *m_modelFunc = m_modelcopiedPdf.asTF(*m_modelobs, *m_modelPars, *m);
    TH1 *hDimuM_data = unbinned_M_Dimu_data->createHistogram("h_mdata", *m, Binning(Mass_Binning, Low_Mass, High_Mass));
    hDimuM_data->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    hDimuM_data->GetYaxis()->SetTitle("d#it{N}/#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");

    RooArgSet *pt_modelcopyOfEverything = static_cast<RooArgSet *>(RooArgSet(*pt_model).snapshot(true)); // True means copy the PDF and everything it depends on
    auto &pt_modelcopiedPdf = static_cast<RooAbsPdf &>((*pt_modelcopyOfEverything)["pt_model"]);         // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
    RooArgSet *pt_modelobs = pt_modelcopiedPdf.getObservables(*unbinned_combData_set);
    RooArgSet *pt_modelPars = pt_modelcopiedPdf.getParameters(*pt_modelobs);
    TF1 *pt_modelFunc = pt_modelcopiedPdf.asTF(*pt_modelobs, *pt_modelPars, *pt);
    TH1 *hDimuPt_data = unbinned_Pt_Dimu_data->createHistogram("h_ptdata", *pt, Binning(Pt_Binning, 0, High_Pt));
    hDimuPt_data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hDimuPt_data->GetYaxis()->SetTitle("d#it{N}/#it{p}_{T} (GeV/#it{c})^{-1}");

    m_frame->SetMaximum(1.2e+6);
    m_frame->SetMinimum(1.2e-2);
    TCanvas *M_Canvas = printRooPlot_ratio(m_frame, kFALSE, r, 0, "pdfmass", m_modelFunc, hDimuM_data, Low_Mass, High_Mass);

    if (withDY)
        M_Canvas->SetName(Form("Mcomp_fit_M_%d_%d_Pt_%d_%d_withDY", Low_Mass, High_Mass, Low_Pt, High_Pt));
    else
        M_Canvas->SetName(Form("Mcomp_fit_M_%d_%d_Pt_%d_%d_noDY", Low_Mass, High_Mass, Low_Pt, High_Pt));

    pt_frame->SetMaximum(1.2e+6);
    pt_frame->SetMinimum(1.2e-2);
    TCanvas *Pt_Canvas = printRooPlot_ratio(pt_frame, kFALSE, r, 0, "pdfpt", pt_modelFunc, hDimuPt_data, Low_Pt, High_Pt);
    if (withDY)
        Pt_Canvas->SetName(Form("Ptcomp_fit_M_%d_%d_Pt_%d_%d_withDY", Low_Mass, High_Mass, Low_Pt, High_Pt));
    else
        Pt_Canvas->SetName(Form("Ptcomp_fit_M_%d_%d_Pt_%d_%d_noDY", Low_Mass, High_Mass, Low_Pt, High_Pt));
}

void pdf_extraction(Int_t Mass_Binning = 19, Int_t Low_Mass = 11, Int_t High_Mass = 30, Int_t Pt_Binning = 15, Int_t Low_Pt = 0, Int_t High_Pt = 30)
{
    gROOT->ProcessLineSync(".x /home/michele_pennisi/dimuon_HF_pp/fit_data/fit_library/PtMassExpPdf.cxx+");
    gROOT->ProcessLineSync(".x /home/michele_pennisi/dimuon_HF_pp/fit_data/fit_library/PtMassPol1ExpPdf.cxx+");

    const Int_t n_DiMuSelection = 4;
    TString Name_DimuSel[n_DiMuSelection] = {"Charm", "Beauty", "HF_Mixed", "DY"};

    TFile *fIn[n_DiMuSelection] = {new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i1/Version3_AliAOD/save_output/LHC23i1_MC_output_Tree_merged.root"), new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version3_AliAOD/save_output/LHC23i2_MC_output_Tree_merged.root"), new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version3_AliAOD/save_output/LHC23i2_MC_output_Tree_merged.root"), new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Powheg_Sim/powheg_DY_mass_3_35/Version1/Analysis_MCsim/powheg_DY_mass_3_35_Analysis_MCsim_Tree_merged.root")};

    TTree *fTree[n_DiMuSelection] = {(TTree *)fIn[0]->Get("DiMuon_Rec_FullMass_PowhegOnlyCharm"), (TTree *)fIn[1]->Get("DiMuon_Rec_FullMass_PowhegOnlyBeauty"), (TTree *)fIn[2]->Get("DiMuon_Rec_FullMass_PowhegOnlyHF_Mixed"), (TTree *)fIn[3]->Get("rec_tree_muDY")};

    gROOT->cd();

    TTree *fTree_cut[n_DiMuSelection];
    TFile *fOut = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/Powheg_pdfMC_unbinnedprova.root", "UPDATE");

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

    RooRealVar *m = new RooRealVar("m", "#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})", Low_Mass, High_Mass);
    RooRealVar *pt = new RooRealVar("pt", "#it{p}_{T} (GeV/#it{c})", Low_Pt, High_Pt);

    Color_t color[n_DiMuSelection] = {kMagenta + 2, kSpring - 6, kAzure + 9, kOrange + 7};
    Color_t fillcolor[n_DiMuSelection] = {kMagenta - 10, kGreen - 10, kCyan - 10, kOrange + 7};

    RooWorkspace *w = new RooWorkspace(Form("w_M_%d_%d_Pt_%d_%d", Low_Mass, High_Mass, Low_Pt, High_Pt), "workspace");

    TString Mass_Factory_Info[n_DiMuSelection];
    Mass_Factory_Info[0].Form("PtMassExpPdf::pdfDimuMassFromCharm(%s[%d,%d], B_DimuMassFromCharm[2.85,1e-3,100], n1_DimuMassFromCharm[2.81,1e-3,100], n2_DimuMassFromCharm[5,1e-3,100])", m->GetName(), Low_Mass, High_Mass);
    Mass_Factory_Info[1].Form("PtMassExpPdf::pdfDimuMassFromBeauty(%s[%d,%d], B_DimuMassFromBeauty[4.16,1e-3,100], n1_DimuMassFromBeauty[2.256,1e-3,100], n2_DimuMassFromBeauty[2.83,1e-3,100])", m->GetName(), Low_Mass, High_Mass);
    Mass_Factory_Info[2].Form("PtMassExpPdf::pdfDimuMassFromHF_Mixed(%s[%d,%d], B_DimuMassFromHF_Mixed[2.65,1e-3,100], n1_DimuMassFromHF_Mixed[2.81,1e-3,100], n2_DimuMassFromHF_Mixed[5,1e-3,100])", m->GetName(), Low_Mass, High_Mass);
    Mass_Factory_Info[3].Form("PtMassExpPdf::pdfDimuMassFromDY(%s[%d,%d], B_DimuMassFromDY[2.65,1e-3,100], n1_DimuMassFromDY[2.81,1e-3,100], n2_DimuMassFromDY[5,1e-3,100])", m->GetName(), Low_Mass, High_Mass);

    TString Pt_Factory_Info[n_DiMuSelection];
    Pt_Factory_Info[0].Form("PtMassExpPdf::pdfDimuPtFromCharm(%s[%d,%d], B_DimuPtFromCharm[2.85,1e-3,100], n1_DimuPtFromCharm[2.81,1e-3,100], n2_DimuPtFromCharm[2.43,1e-3,100])", pt->GetName(), Low_Pt, High_Pt);
    Pt_Factory_Info[1].Form("PtMassExpPdf::pdfDimuPtFromBeauty(%s[%d,%d], B_DimuPtFromBeauty[2.85,1e-3,100], n1_DimuPtFromBeauty[2.81,1e-3,100], n2_DimuPtFromBeauty[2.43,1e-3,100])", pt->GetName(), Low_Pt, High_Pt);
    Pt_Factory_Info[2].Form("PtMassExpPdf::pdfDimuPtFromHF_Mixed(%s[%d,%d], B_DimuPtFromHF_Mixed[2.85,1e-3,100], n1_DimuPtFromHF_Mixed[2.81,1e-3,100], n2_DimuPtFromHF_Mixed[2.43,1e-3,100])", pt->GetName(), Low_Pt, High_Pt);
    Pt_Factory_Info[3].Form("PtMassExpPdf::pdfDimuPtFromDY(%s[%d,%d], B_DimuPtFromDY[2.85,1e-3,100], n1_DimuPtFromDY[2.81,1e-3,100], n2_DimuPtFromDY[2.43,1e-3,100])", pt->GetName(), Low_Pt, High_Pt);

    RooAbsPdf *pdfDimuPt[n_DiMuSelection];
    RooAbsPdf *pdfDimuM[n_DiMuSelection];
    TF1 *pt_Func[n_DiMuSelection];
    TF1 *m_Func[n_DiMuSelection];
    TCanvas *pt_canvas[n_DiMuSelection];
    TCanvas *m_canvas[n_DiMuSelection];

    for (Int_t i_DimuSel = 0; i_DimuSel < n_DiMuSelection; i_DimuSel++)
    {
        // Select the region of the MC distribution to be extracted and creation of the RooDataSet

        fTree_cut[i_DimuSel] = (TTree *)fTree[i_DimuSel]->CopyTree(Form("(m>%d && m<%d) && (pt > %d && pt <%d)", Low_Mass, High_Mass, Low_Pt, High_Pt));
        fTree_cut[i_DimuSel]->SetName(Form("Tree_%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), Low_Mass, High_Mass, Low_Pt, High_Pt));
        fTree_cut[i_DimuSel]->GetBranch("m")->SetName(m->GetName());
        fTree_cut[i_DimuSel]->GetBranch("pt")->SetName(pt->GetName());
        M_Dimu[i_DimuSel] = new RooDataSet(Form("M_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), Form("M_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), RooArgSet(*m), Import(*fTree_cut[i_DimuSel]));
        Pt_Dimu[i_DimuSel] = new RooDataSet(Form("Pt_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), Form("Pt_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), RooArgSet(*pt), Import(*fTree_cut[i_DimuSel]));

        w->factory(Mass_Factory_Info[0]);
        w->factory(Mass_Factory_Info[1]);
        w->factory(Mass_Factory_Info[2]);
        w->factory(Mass_Factory_Info[3]);

        w->factory(Pt_Factory_Info[0]);
        w->factory(Pt_Factory_Info[1]);
        w->factory(Pt_Factory_Info[2]);
        w->factory(Pt_Factory_Info[3]);

        pdfDimuPt[i_DimuSel] = w->pdf(Form("pdfDimuPtFrom%s", Name_DimuSel[i_DimuSel].Data()));
        auto result1 = pdfDimuPt[i_DimuSel]->fitTo(*Pt_Dimu[i_DimuSel], Minimizer("Minuit2"), Save(), SumW2Error(true));

        pdfDimuM[i_DimuSel] = w->pdf(Form("pdfDimuMassFrom%s", Name_DimuSel[i_DimuSel].Data()));
        auto result2 = pdfDimuM[i_DimuSel]->fitTo(*M_Dimu[i_DimuSel], Minimizer("Minuit2"), Save(), SumW2Error(true));

        // Drawing the data and pdf on RooPlot

        frameDimuMass[i_DimuSel] = m->frame(Title(Form("frameDimuMass_%s", Name_DimuSel[i_DimuSel].Data())));
        frameDimuMass[i_DimuSel]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})^{-1}");
        M_Dimu[i_DimuSel]->plotOn(frameDimuMass[i_DimuSel], Name(Form("M_%s", Name_DimuSel[i_DimuSel].Data())), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i_DimuSel]), LineColor(color[i_DimuSel]), LineWidth(2), DrawOption("PEZ"), Binning(Mass_Binning));
        pdfDimuM[i_DimuSel]->plotOn(frameDimuMass[i_DimuSel], Name(Form("pdfDimuMFrom%s", Name_DimuSel[i_DimuSel].Data())), LineStyle(kSolid), LineColor(color[i_DimuSel]));

        frameDimuPt[i_DimuSel] = pt->frame(Title(Form("frameDimuPt_%s", Name_DimuSel[i_DimuSel].Data())));
        frameDimuPt[i_DimuSel]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
        Pt_Dimu[i_DimuSel]->plotOn(frameDimuPt[i_DimuSel], Name(Form("Pt_%s", Name_DimuSel[i_DimuSel].Data())), MarkerSize(1.75), MarkerStyle(24), MarkerColor(color[i_DimuSel]), LineColor(color[i_DimuSel]), LineWidth(2), DrawOption("PEZ"), Binning(Pt_Binning));
        pdfDimuPt[i_DimuSel]->plotOn(frameDimuPt[i_DimuSel], Name(Form("pdfDimuPtFrom%s", Name_DimuSel[i_DimuSel].Data())), LineStyle(kSolid), LineColor(color[i_DimuSel]));

        // Conversion of the RooDataSet in TH1 & RooPdf in a TF1 object to calculate the ratio

        m_histo[i_DimuSel] = M_Dimu[i_DimuSel]->createHistogram(Form("h_M_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), *m, Binning(Mass_Binning, Low_Mass, High_Mass));
        pt_histo[i_DimuSel] = Pt_Dimu[i_DimuSel]->createHistogram(Form("h_Pt_Dimu_%s", Name_DimuSel[i_DimuSel].Data()), *pt, Binning(Pt_Binning, Low_Pt, High_Pt));

        RooArgSet *pt_model = static_cast<RooArgSet *>(RooArgSet(*pdfDimuPt[i_DimuSel]).snapshot(true));                       // True means copy the PDF and everything it depends on
        auto &pt_modelcopied = static_cast<RooAbsPdf &>((*pt_model)[Form("pdfDimuPtFrom%s", Name_DimuSel[i_DimuSel].Data())]); // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
        RooArgSet *pt_modelobs = pt_modelcopied.getObservables(*Pt_Dimu[i_DimuSel]);
        RooArgSet *pt_modelPars = pt_modelcopied.getParameters(*pt_modelobs);
        pt_Func[i_DimuSel] = pt_modelcopied.asTF(*pt_modelobs, *pt_modelPars, *pt);

        RooArgSet *m_model = static_cast<RooArgSet *>(RooArgSet(*pdfDimuM[i_DimuSel]).snapshot(true));                         // True means copy the PDF and everything it depends on
        auto &m_modelcopied = static_cast<RooAbsPdf &>((*m_model)[Form("pdfDimuMassFrom%s", Name_DimuSel[i_DimuSel].Data())]); // Get back the copied pdf. It lives in the RooArgSet "copyOfEverything"
        RooArgSet *m_modelobs = m_modelcopied.getObservables(*M_Dimu[i_DimuSel]);
        RooArgSet *m_modelPars = m_modelcopied.getParameters(*m_modelobs);
        m_Func[i_DimuSel] = m_modelcopied.asTF(*m_modelobs, *m_modelPars, *m);

        m_canvas[i_DimuSel] = printMC_ratio(Form("m_canvas_%s", Name_DimuSel[i_DimuSel].Data()), frameDimuMass[i_DimuSel], m_histo[i_DimuSel], m_Func[i_DimuSel], color[i_DimuSel], Low_Mass, High_Mass);
        pt_canvas[i_DimuSel] = printMC_ratio(Form("pt_canvas_%s", Name_DimuSel[i_DimuSel].Data()), frameDimuPt[i_DimuSel], pt_histo[i_DimuSel], pt_Func[i_DimuSel], color[i_DimuSel], Low_Pt, High_Pt);

        fOut->cd();
        fTree_cut[i_DimuSel]->Write();

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

        Param_Pt_unbinned[i_DimuSel] = new TH1D(Form("Param_Pt_unbinned_from%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), Low_Mass, High_Mass, Low_Pt, High_Pt), "; coeff x Mass fit", 4, 0, 4);
        Param_Pt_unbinned[i_DimuSel]->SetBinContent(1, B_DimuPt[i_DimuSel]->getVal());
        Param_Pt_unbinned[i_DimuSel]->SetBinContent(2, n1_DimuPt[i_DimuSel]->getVal());
        Param_Pt_unbinned[i_DimuSel]->SetBinContent(3, n2_DimuPt[i_DimuSel]->getVal());
        Param_Pt_unbinned[i_DimuSel]->SetBinContent(4, 1);

        Param_Pt_unbinned[i_DimuSel]->SetBinError(1, B_DimuPt[i_DimuSel]->getError());
        Param_Pt_unbinned[i_DimuSel]->SetBinError(2, n1_DimuPt[i_DimuSel]->getError());
        Param_Pt_unbinned[i_DimuSel]->SetBinError(3, n2_DimuPt[i_DimuSel]->getError());
        Param_Pt_unbinned[i_DimuSel]->SetBinError(4, 1);

        Param_M_unbinned[i_DimuSel] = new TH1D(Form("Param_M_unbinned_from%s_M_%d_%d_Pt_%d_%d", Name_DimuSel[i_DimuSel].Data(), Low_Mass, High_Mass, Low_Pt, High_Pt), "; coeff x Pt fit", 4, 0, 4);
        Param_M_unbinned[i_DimuSel]->SetBinContent(1, B_DimuMass[i_DimuSel]->getVal());
        Param_M_unbinned[i_DimuSel]->SetBinContent(2, n1_DimuMass[i_DimuSel]->getVal());
        Param_M_unbinned[i_DimuSel]->SetBinContent(3, n2_DimuMass[i_DimuSel]->getVal());
        Param_M_unbinned[i_DimuSel]->SetBinContent(4, 1);

        Param_M_unbinned[i_DimuSel]->SetBinError(1, B_DimuMass[i_DimuSel]->getError());
        Param_M_unbinned[i_DimuSel]->SetBinError(2, n1_DimuMass[i_DimuSel]->getError());
        Param_M_unbinned[i_DimuSel]->SetBinError(3, n2_DimuMass[i_DimuSel]->getError());
        Param_M_unbinned[i_DimuSel]->SetBinError(4, 1);
    }

    w->Write();
    // w->writeToFile(Form("/home/michele_pennisi/cernbox/HF_dimuons/fit_data/results/Powheg_pdfMC_unbinnedprova.root"));

    // fTree_cut[0]->Draw("m");
    // new TCanvas();
    // fTree_cut[0]->Draw("pt");
    // fIn[0]->Close();

    // return;
}

TCanvas *printMC_ratio(TString name, RooPlot *frame, TH1 *data, TF1 *pdf, Color_t color, Int_t minx = 0, Int_t max_x = 30)
{
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

    // TH2D *h_grid = new TH2D("h_grid", "", 100, minx, max_x, 100, Lowy, frame->GetMaximum() * 1000);
    // h_grid->SetTitle("");

    // h_grid->GetXaxis()->SetTitleOffset(1.3);
    // h_grid->GetXaxis()->SetTitleSize(0.0475);
    // h_grid->GetXaxis()->SetLabelSize(0.045);

    // h_grid->GetYaxis()->SetNdivisions(505);
    // h_grid->GetYaxis()->SetTitleOffset(0.9);
    // h_grid->GetYaxis()->SetTitleSize(0.065);
    // h_grid->GetYaxis()->SetLabelSize(0.055);

    // h_grid->Draw();
    frame->SetTitle("");
    frame->GetXaxis()->SetTitleOffset(1.3);
    frame->GetXaxis()->SetTitleSize(0.0475);
    frame->GetXaxis()->SetLabelSize(0.045);

    frame->GetYaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetTitleOffset(0.9);
    frame->GetYaxis()->SetTitleSize(0.065);
    frame->GetYaxis()->SetLabelSize(0.055);
    frame->Draw();

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
    letexTitle->DrawLatex(0.355, 0.81, "POWHEG+PYTHIA6, N_{ev} = 1.4 #upoint 10^{8}");
    if (str.Contains("h_Pt"))
    {
        letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{#eta}_{#mu} < 4.0");
        letexTitle->DrawLatex(0.355, 0.66, Form("2.5 < #it{y}_{#mu^{#plus}#mu^{#minus}} < 4.0, 4 < #it{m}_{#mu^{#plus}#mu^{#minus}} < 9 GeV/#it{c}^{2}"));
    }
    else
    {
        letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{#eta}_{#mu} < 4.0");
        letexTitle->DrawLatex(0.355, 0.66, Form("2.5 < #it{y}_{#mu^{#plus}#mu^{#minus}} < 4.0, #it{p}_{T} < 10 GeV/#it{c}"));
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

    // TH2D *h_grid_ratio = new TH2D("h_grid_ratio", "", 100, minx, max_x, 100, -0.55, 2.55);
    // h_grid_ratio->SetTitle("");

    TH1D *c_data = (TH1D *)data->Clone(Form("c_data_%s", data->GetName()));

    c_data->SetTitle("");
    c_data->GetYaxis()->SetTitle(Form("#frac{MC}{PDF}"));
    c_data->GetYaxis()->CenterTitle();
    c_data->GetYaxis()->SetNdivisions(504);
    c_data->GetYaxis()->SetTitleSize(0.08);
    // c_data->GetYaxis()->SetTitleOffset(0.8);
    c_data->GetYaxis()->SetLabelOffset(0.02);
    c_data->GetYaxis()->SetLabelSize(0.1);

    c_data->GetXaxis()->SetTitleSize(0.1);
    c_data->GetXaxis()->SetTitleOffset(1.1);
    c_data->GetXaxis()->SetLabelSize(0.1);
    c_data->GetXaxis()->SetTitle(c_data->GetXaxis()->GetTitle());

    c_data->SetLineColor(kBlack);
    c_data->SetMarkerColor(kBlack);
    c_data->SetMarkerStyle(20);
    // c_data->Rebin(10);
    // c_data->Scale(1. / c_data->Integral(), "width");
    c_data->Divide(pdf);

    c_data->Draw();
    l->Draw();
    l1->Draw();
    l2->Draw();

    // canvas->ls();
    return canvas;
}

TCanvas *printRooPlot_ratio(RooPlot *frame, Bool_t norm, RooFitResult *r, Int_t choice, TString roohist_name, TF1 *pdf, TH1 *data, Double_t minx, Double_t max_x)
{
    const RooArgList &fitParams = r->floatParsFinal();

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
    TString Name_DimuSel[4] = {"Charm", "Beauty", "DY", "HF_Mixed"};
    TString info_label[4];
    info_label[0].Form("#leftarrow c");
    info_label[1].Form("#leftarrow b");
    info_label[2].Form("#leftarrow DY");
    info_label[3].Form("#leftarrow c,b ");

    if (roohist_name.Contains("pt"))
        fit_legend->AddEntry("pdfpt", " ", "L");
    else if (roohist_name.Contains("mass"))
        fit_legend->AddEntry("pdfmass", " ", "L");

    TLegend *pdf_legend = new TLegend(0.175, 0.03, 0.40, 0.305);
    pdf_legend->SetTextAlign(11);
    pdf_legend->SetFillStyle(0);
    pdf_legend->SetBorderSize(0);
    pdf_legend->SetTextSize(0.0425);
    pdf_legend->SetTextAlign(12);
    pdf_legend->SetHeader("POWHEG+PYTHIA6 PDF");
    for (Int_t i = 0; i < fitParams.getSize() + 1; i++)
    {
        pdf_legend->AddEntry(TString::Format("%s%s", roohist_name.Data(), Name_DimuSel[i].Data()), Form("#mu^{#plus}#mu^{#minus} %s", info_label[i].Data()), "L");
    }

    legend->Draw();
    fit_legend->Draw();
    pdf_legend->Draw();

    TLatex *letexTitle = new TLatex();
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);
    letexTitle->SetTextSize(0.0475);
    letexTitle->DrawLatex(0.175, 0.875, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.175, 0.785, "LHC18p period");
    letexTitle->SetTextSize(0.0425);
    vector<double> fit_result;
    Double_t start_y_latex = 0.825;
    Int_t i_par = 0;
    while (i_par < fitParams.getSize())
    {
        auto &fitPar = (RooRealVar &)fitParams[fitParams.getSize() - i_par - 1];
        fit_result.push_back(fitPar.getVal());
        info_label[i_par].Append(Form("} = %0.1f #pm %0.1f", fitPar.getVal(), fitPar.getError()));
        letexTitle->DrawLatex(0.625, start_y_latex - 0.1 * i_par, Form("#it{N}_{#mu^{#plus}#mu^{#minus}%s", info_label[i_par].Data()));
        i_par++;
    }

    info_label[i_par].Append(Form("} = 0"));
    letexTitle->DrawLatex(0.625, start_y_latex - 0.1 * i_par, Form("#it{N}^{fixed}_{#mu^{#plus}#mu^{#minus}%s", info_label[i_par].Data()));

    if (roohist_name.Contains("pt"))
    {
        letexTitle->DrawLatex(0.175, 0.695, "Reconstructed #mu^{#plus}#mu^{#minus}, #it{m}_{#mu^{#plus}#mu^{#minus}} > 4 GeV/#it{c}^{2}");
        letexTitle->DrawLatex(0.175, 0.605, "2.5 < #it{y}_{#mu} < 4.0");
    }
    else
        letexTitle->DrawLatex(0.175, 0.695, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{y}_{#mu} < 4.0");

    pad2->cd();
    pad2->SetTicks();
    TLine *l = new TLine(minx, 1.0, max_x, 1.0);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->SetLineColor(kRed);

    TLine *l1 = new TLine(minx, 0.5, max_x, 0.5);
    l1->SetLineWidth(2);
    l1->SetLineStyle(9);
    l1->SetLineColor(kGray + 2);

    TLine *l2 = new TLine(minx, 1.5, max_x, 1.5);
    l2->SetLineWidth(2);
    l2->SetLineStyle(9);
    l2->SetLineColor(kGray + 2);

    TH2D *h_grid_ratio = new TH2D("h_grid", "", 100, minx, max_x, 100, -0.5, 2.5);
    h_grid_ratio->SetName(Form("%s_ratiodatafit", roohist_name.Data()));
    h_grid_ratio->SetTitle("");
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
    h_NDY_event->GetYaxis()->SetRangeUser(0.1,1.2e+4);
    h_NDY_event->Draw();
    c->cd(2);
    gPad->SetLogy();
    h_NDY_event_fwd->GetYaxis()->SetRangeUser(0.1,1.2e+4);
    h_NDY_event_fwd->Draw();
    c->cd(3);
    h_YDY->Draw();
    c->cd(4);
    h_YDY_fwd->Draw();


    TFile *fOut=new TFile("~/cernbox/HF_dimuons/fit_data/ingredient_cs_powheg/DY_cs.root","RECREATE");
    fOut->cd();
    h_NDY_event->Write();
    h_NDY_event_fwd->Write();
    
}