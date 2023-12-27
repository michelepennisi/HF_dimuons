#include "/home/michele_pennisi/cernbox/common_include.h"
void mc_comparinson_shape();
void mc_comparinson_dimuon();
void mc_comparinson_hfquark();
void mc_comparinson_Rec_Dimu();

TFile *fIn_PYTHIA_CENTRALIZED_MB = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC22c1/294925/output/Merged_LHC22c1_MC_output_Hist_294925.root", "READ");
TH1F *h_Nevents_PYTHIA_MB_CENTRALIZED = (TH1F *)fIn_PYTHIA_CENTRALIZED_MB->Get("h_Nevents");
Int_t Nev_PYTHIA_MB_CENTRALIZED = (Int_t)h_Nevents_PYTHIA_MB_CENTRALIZED->GetBinContent(2);
// cout<<Nev_PYTHIA_MB_CENTRALIZED<<endl;
Double_t Norm_PYTHIA_MB_CENTRALIZED = (78.05 / ((Double_t)Nev_PYTHIA_MB_CENTRALIZED)); // in mb

TFile *fIn_PYTHIA_HF = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Hist_fromSim/Version1/HF/HistLite_HF_MCDimuHFTree_merged.root", "READ");

TFile *fIn_Pythia_STANDALONE_NoDiffr = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_Def_MC_output_Hist_100000.root", "READ");
TH1F *h_Nevents_PYTHIA_STANDALONE_NoDiffr = (TH1F *)fIn_Pythia_STANDALONE_NoDiffr->Get("h_Nevents");

Int_t Nev_PYTHIA_STANDALONE_NoDiffr = 999;
Double_t Norm_PYTHIA_STANDALONE_NoDiffr = (56.42 / ((Double_t)Nev_PYTHIA_STANDALONE_NoDiffr)); // in mb

TFile *fIn_Pythia_STANDALONE_EL = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/new_pythia_sim/SoftQCD_inel_LFoff_Def_pythia_sim_2411_DefaultBR_output_Hist_100000.root", "READ");
TH1F *h_Nevents_PYTHIA_STANDALONE_EL = (TH1F *)fIn_Pythia_STANDALONE_EL->Get("h_Nevents");
Int_t Nev_PYTHIA_STANDALONE_EL = (Int_t)h_Nevents_PYTHIA_STANDALONE_EL->GetBinContent(2);
Double_t Norm_PYTHIA_STANDALONE_EL = (78.05 / ((Double_t)Nev_PYTHIA_STANDALONE_EL)); // in mb

TFile *fIn_Powheg_Charm = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i1/Version3_AliAOD/save_output/LHC23i1_MC_output_Hist_merged.root", "READ");
// TFile *fIn_Powheg_Charm = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/on_farm_test/LHC23i1_MC_output_Hist_294009.root", "READ");

TH1F *h_Nevents_Powheg_Charm = (TH1F *)fIn_Powheg_Charm->Get("h_Nevents");
Int_t Nev_Powheg_Charm = (Int_t)h_Nevents_Powheg_Charm->GetBinContent(2);
Double_t Norm_Powheg_Charm = (56.42 / ((Double_t)Nev_Powheg_Charm)); // in mb

TFile *fIn_Powheg_Beauty = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version3_AliAOD/save_output/LHC23i2_MC_output_Hist_merged.root", "READ");

void mc_comparison_DiMu_Rec()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TString DiMuon_fromLF_Generator[3];
    TString DiMuon_fromLF_Generator_Title[3] = {"#mu#mu <- Geant LF", "#mu#mu <- Geant,PYTHIA LF", "#mu#mu <- PYTHIA LF"};
    Color_t DiMuon_fromLF_Generator_Color[3] = {kGreen + 2, kCyan + 2, kMagenta + 2};

    DiMuon_fromLF_Generator[0].Form("GeantOnly");
    DiMuon_fromLF_Generator[1].Form("PYTHIAGeant");
    DiMuon_fromLF_Generator[2].Form("PYTHIAOnly");

    TH2F *h_PtM_DiMuon_Rec_fromLF[3];
    TH2F *h_PtY_DiMuon_Rec_fromLF[3];

    TH1F *h_Pt_DiMuon_Rec_fromLF[3];
    TH1F *h_M_DiMuon_Rec_fromLF[3];
    TH1F *h_Y_DiMuon_Rec_fromLF[3];

    for (Int_t i_LF_Generator = 0; i_LF_Generator < 3; i_LF_Generator++)
    {
        h_PtM_DiMuon_Rec_fromLF[i_LF_Generator] = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get(Form("DiMuon_Rec/LF_origin/h_PtM_DiMuon_Rec_fromLF_%s", DiMuon_fromLF_Generator[i_LF_Generator].Data()));
        h_PtY_DiMuon_Rec_fromLF[i_LF_Generator] = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get(Form("DiMuon_Rec/LF_origin/h_PtY_DiMuon_Rec_fromLF_%s", DiMuon_fromLF_Generator[i_LF_Generator].Data()));

        h_Pt_DiMuon_Rec_fromLF[i_LF_Generator] = (TH1F *)h_PtM_DiMuon_Rec_fromLF[i_LF_Generator]->ProjectionX();
        h_Pt_DiMuon_Rec_fromLF[i_LF_Generator]->SetTitle(DiMuon_fromLF_Generator_Title[i_LF_Generator].Data());
        hist1D_graphic_opt(h_Pt_DiMuon_Rec_fromLF[i_LF_Generator], kTRUE, 1, 20, DiMuon_fromLF_Generator_Color[i_LF_Generator], 1.);
        h_M_DiMuon_Rec_fromLF[i_LF_Generator] = (TH1F *)h_PtM_DiMuon_Rec_fromLF[i_LF_Generator]->ProjectionY();
        h_M_DiMuon_Rec_fromLF[i_LF_Generator]->SetTitle(DiMuon_fromLF_Generator_Title[i_LF_Generator].Data());
        hist1D_graphic_opt(h_M_DiMuon_Rec_fromLF[i_LF_Generator], kTRUE, 1, 20, DiMuon_fromLF_Generator_Color[i_LF_Generator], 1.);
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator] = (TH1F *)h_PtY_DiMuon_Rec_fromLF[i_LF_Generator]->ProjectionY();
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->SetTitle(DiMuon_fromLF_Generator_Title[i_LF_Generator].Data());
        hist1D_graphic_opt(h_Y_DiMuon_Rec_fromLF[i_LF_Generator], kTRUE, 1, 20, DiMuon_fromLF_Generator_Color[i_LF_Generator], 1.);
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->GetXaxis()->SetRangeUser(-4, -2.5);
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->GetYaxis()->SetTitle("Entries");
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->GetYaxis()->SetTitleOffset(1.2);
        h_Y_DiMuon_Rec_fromLF[i_LF_Generator]->GetXaxis()->SetTitleOffset(1.15);
    }

    TCanvas *C_Y = canvas_noratio_divide2("C_Y");
    C_Y->cd();
    gPad->SetLogy();
    h_Y_DiMuon_Rec_fromLF[2]->Draw("PE");
    h_Y_DiMuon_Rec_fromLF[0]->Draw("PESAME");
    h_Y_DiMuon_Rec_fromLF[1]->Draw("PESAME");
    TLegend *leg = C_Y->BuildLegend(0.7, 0.75, 0.9, 0.95, " ");
    leg->SetFillStyle(0);
    leg->SetLineColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.0225);
    // C_Y->SaveAs("images/LF_DiMuRec_Y.png");
}

void mc_comparison_Mu_Rec()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    fIn_PYTHIA_CENTRALIZED_MB->cd("Muon_Rec/Geant");
    fIn_PYTHIA_CENTRALIZED_MB->ls();

    // fIn_PYTHIA_CENTRALIZED_MB->cd("Muon_Rec/not_Geant");
    // fIn_PYTHIA_CENTRALIZED_MB->ls();
    TH2F *h_VzmotherEta_Muon_Rec_PYTHIA_LF = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/not_Geant/h_VzmotherEta_Muon_Rec_PYTHIA_LF");
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

    TH1F *h_Nperevent_Muon_Rec_LF_PYTHIAOnly = (TH1F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/not_Geant/h_Nperevent_Muon_Rec_LF_PYTHIAOnly");
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

    TH2F *h_PtEta_Muon_Rec_LF_PYTHIAOnly = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/not_Geant/h_PtEta_Muon_RecLF_PYTHIAOnly");
    TH1F *h_Pt_Muon_Rec_LF_PYTHIAOnly = (TH1F *)h_PtEta_Muon_Rec_LF_PYTHIAOnly->ProjectionX();
    hist1D_graphic_opt(h_Pt_Muon_Rec_LF_PYTHIAOnly, kTRUE, 1, 20, kMagenta + 2, 1.);
    h_Pt_Muon_Rec_LF_PYTHIAOnly->GetXaxis()->SetRangeUser(0, 16);
    h_Pt_Muon_Rec_LF_PYTHIAOnly->SetTitle("from PYTHIA hadron");
    h_Pt_Muon_Rec_LF_PYTHIAOnly->GetYaxis()->SetTitle("Entries");
    cout<<"h_Pt_Muon_Rec_LF_PYTHIAOnly Entries() "<<h_Pt_Muon_Rec_LF_PYTHIAOnly->GetEntries()<<endl;

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
    cout<<"h_Pt_Muon_Rec_LF_GeantOnly Entries() "<<h_Pt_Muon_Rec_LF_GeantOnly->GetEntries()<<endl;

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
    TCanvas *c_Pt_Rec = two_histo_ratio(h_Pt_Muon_Rec_LF_PYTHIAOnly, h_Pt_Muon_Rec_LF_GeantOnly, ratio_Pt_Geant_PYTHIA_CENTRALIZED, "c_Pt_Rec", "#mu #leftarrow LF hadrons", kTRUE, kTRUE);
    c_Pt_Rec->SaveAs("images/LF_MuRec_Pt.png");

    TH1F *ratio_Y_Geant_PYTHIA_CENTRALIZED = (TH1F *)h_Eta_Muon_Rec_LF_GeantOnly->Clone("ratio_Y_Geant_PYTHIA_CENTRALIZED");
    ratio_Y_Geant_PYTHIA_CENTRALIZED->Divide(h_Eta_Muon_Rec_LF_PYTHIAOnly);
    ratio_Y_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("Geant/PYTHIA");

    h_Eta_Muon_Rec_LF_PYTHIAOnly->GetYaxis()->SetRangeUser(0.8, h_Eta_Muon_Rec_LF_PYTHIAOnly->GetMaximum() * 400.8);
    ratio_Y_Geant_PYTHIA_CENTRALIZED->GetYaxis()->CenterTitle();
    ratio_Y_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(ratio_Y_Geant_PYTHIA_CENTRALIZED->GetMinimum() * 0.3, ratio_Y_Geant_PYTHIA_CENTRALIZED->GetMaximum() * 1.2);
    TCanvas *c_Y_Rec = two_histo_ratio(h_Eta_Muon_Rec_LF_PYTHIAOnly, h_Eta_Muon_Rec_LF_GeantOnly, ratio_Y_Geant_PYTHIA_CENTRALIZED, "c_Y_Rec", "#mu #leftarrow LF hadrons", kTRUE, kTRUE);
    c_Y_Rec->SaveAs("images/LF_MuRec_Eta.png");
    // c_ETA_Rec->SaveAs(Form("images/%s.png", c_ETA_Rec->GetName()));
}

void mc_comparison_Mu_Gen()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // std::cout<<Nev_PYTHIA_STANDALONE_NoDiffr<<endl;
    fIn_PYTHIA_CENTRALIZED_MB->cd("Hadron/Geant");
    fIn_PYTHIA_CENTRALIZED_MB->ls();
    fIn_PYTHIA_CENTRALIZED_MB->cd("Hadron/not_Geant");
    fIn_PYTHIA_CENTRALIZED_MB->ls();
    // return;
    const Int_t n_PDG_tested = 3;
    TString canvas_header[n_PDG_tested] = {"#mu #leftarrow #pi^{#plus}", "#mu #leftarrow K^{0}", "#mu #leftarrow K^{#plus}"};
    TString canvas_name[n_PDG_tested] = {"mu_pi_plus", "mu_kzero", "mu_kplus"};
    Int_t low_PDG[n_PDG_tested] = {210, 310, 320};
    Int_t high_PDG[n_PDG_tested] = {220, 320, 330};

    TH2F *h_EtaPdg_Muon_Gen_PYTHIAOnly_CENTRALIZED = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Gen/not_Geant/h_EtaPdg_Muon_Gen_PYTHIAOnly");
    h_EtaPdg_Muon_Gen_PYTHIAOnly_CENTRALIZED->SetName(" h_EtaPdg_Muon_Gen_PYTHIAOnly_CENTRALIZED");
    TH2F *h_EtaPdg_Muon_Gen_GeantOnly_CENTRALIZED = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Gen/Geant/h_EtaPdg_Muon_Gen_GeantOnly");
    h_EtaPdg_Muon_Gen_GeantOnly_CENTRALIZED->SetName("h_EtaPdg_Muon_Gen_GeantOnly_CENTRALIZED");

    TH2F *h_EtaPdg_Muon_Gen_STANDALONE = (TH2F *)fIn_Pythia_STANDALONE_NoDiffr->Get("Muon_Gen/h_EtaPdg_Muon_Gen");
    h_EtaPdg_Muon_Gen_STANDALONE->SetName("h_EtaPdg_Muon_Gen_STANDALONE");
    TH1F *h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED;
    TH1F *h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED;
    TH1F *h_Eta_Hadron_Gen_LF_STANDALONE;

    for (Int_t i_PDG = 0; i_PDG < n_PDG_tested; i_PDG++)
    {
        h_EtaPdg_Muon_Gen_PYTHIAOnly_CENTRALIZED->GetYaxis()->SetRangeUser(low_PDG[i_PDG], high_PDG[i_PDG]);

        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED = (TH1F *)h_EtaPdg_Muon_Gen_PYTHIAOnly_CENTRALIZED->ProjectionX();
        hist1D_graphic_opt(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED, kTRUE, 2, 20, kMagenta + 2, Norm_PYTHIA_MB_CENTRALIZED);
        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->SetName(Form("h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED_%s", canvas_name[i_PDG].Data()));
        std::cout << h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetEntries() << endl;
        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->SetTitle("PYTHIA8 centralized (PYTHIA only)");
        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("d#sigma/d#eta (mb)");

        h_EtaPdg_Muon_Gen_GeantOnly_CENTRALIZED->GetYaxis()->SetRangeUser(low_PDG[i_PDG], high_PDG[i_PDG]);

        h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED = (TH1F *)h_EtaPdg_Muon_Gen_GeantOnly_CENTRALIZED->ProjectionX();
        h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->SetName(Form("h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED_%s", canvas_name[i_PDG].Data()));
        h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->SetTitle("PYTHIA8 centralized (Geant only)");
        std::cout << h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->GetEntries() << endl;
        hist1D_graphic_opt(h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED, kTRUE, 2, 20, kGreen + 2, Norm_PYTHIA_MB_CENTRALIZED);

        //-----------------------------------------//

        h_EtaPdg_Muon_Gen_STANDALONE->SetTitle("PYTHIA8 stand-alone");
        h_EtaPdg_Muon_Gen_STANDALONE->GetYaxis()->SetRangeUser(low_PDG[i_PDG], high_PDG[i_PDG]);

        h_Eta_Hadron_Gen_LF_STANDALONE = (TH1F *)h_EtaPdg_Muon_Gen_STANDALONE->ProjectionX();
        h_Eta_Hadron_Gen_LF_STANDALONE->SetName(Form("h_Eta_Hadron_Gen_LF_STANDALONE_%s", canvas_name[i_PDG].Data()));
        hist1D_graphic_opt(h_Eta_Hadron_Gen_LF_STANDALONE, kTRUE, 2, 20, kCyan + 2, Norm_PYTHIA_STANDALONE_NoDiffr);

        //-----------------------------------------//

        // TH2F *h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm = (TH2F *)fIn_Powheg_Charm->Get("Hadron/not_Geant/h_PdgEta_HFHadron_prompt_PYTHIAOnly");
        // h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm->SetName(" h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm");
        // h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm->GetXaxis()->SetRangeUser(210, 220);

        // TH1F *h_Eta_Hadron_Gen_LF_PYTHIAOnly_Powheg_Charm = (TH1F *)h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm->ProjectionY();
        // hist1D_graphic_opt(h_Eta_Hadron_Gen_LF_PYTHIAOnly_Powheg_Charm, kTRUE, 10, 20, kBlack, Norm_Powheg_Charm);

        TH1F *ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE = (TH1F *)h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->Clone("ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE");
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE->Divide(h_Eta_Hadron_Gen_LF_STANDALONE);
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE->GetYaxis()->SetTitle("CENTRALIZED/STANDALONE");

        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetMinimum() * 0.6, h_Eta_Hadron_Gen_LF_STANDALONE->GetMaximum() * 400.8);
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE->GetYaxis()->CenterTitle();
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE->GetYaxis()->SetRangeUser(ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE->GetMinimum() * 0.6, ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE->GetMaximum() * 1.2);

        TCanvas *c_ETA_GEN = two_histo_ratio(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED, h_Eta_Hadron_Gen_LF_STANDALONE, ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE, canvas_name[i_PDG].Data(), canvas_header[i_PDG].Data(), kTRUE, kTRUE);
        c_ETA_GEN->SaveAs(Form("images/%s.png", c_ETA_GEN->GetName()));

        TH1F *ratio_ETA_Geant_PYTHIA_CENTRALIZED = (TH1F *)h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->Clone("ratio_ETA_Geant_PYTHIA_CENTRALIZED");
        ratio_ETA_Geant_PYTHIA_CENTRALIZED->Divide(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED);
        ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("Geant/PYTHIA");

        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->GetMinimum() * 0.2, h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetMaximum() * 400.8);
        ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetYaxis()->CenterTitle();
        ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetMinimum() * 0.3, ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetMaximum() * 1.2);
        TCanvas *c_ETA_GEANT_PYTHIA = two_histo_ratio(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED, h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED, ratio_ETA_Geant_PYTHIA_CENTRALIZED, Form("Geant_PYTHIA_%s", canvas_name[i_PDG].Data()), canvas_header[i_PDG].Data(), kTRUE, kTRUE);
        c_ETA_GEANT_PYTHIA->SaveAs(Form("images/%s.png", c_ETA_GEANT_PYTHIA->GetName()));
    }
}

void mc_comparison()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // std::cout<<Nev_PYTHIA_STANDALONE_NoDiffr<<endl;
    fIn_PYTHIA_CENTRALIZED_MB->cd("Hadron/Geant");
    fIn_PYTHIA_CENTRALIZED_MB->ls();
    fIn_PYTHIA_CENTRALIZED_MB->cd("Hadron/not_Geant");
    fIn_PYTHIA_CENTRALIZED_MB->ls();
    // return;
    const Int_t n_PDG_tested = 3;
    TString canvas_header[n_PDG_tested] = {"#pi^{#plus}", "K^{0}", "K^{#plus}"};
    TString canvas_name[n_PDG_tested] = {"pi_plus", "kzero", "kplus"};
    Int_t low_PDG[n_PDG_tested] = {211, 310, 320};
    Int_t high_PDG[n_PDG_tested] = {212, 320, 330};

    TH2F *h_PdgEta_HFHadron_prompt_PYTHIAOnly_CENTRALIZED = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Hadron/not_Geant/h_PdgEta_HFHadron_prompt_PYTHIAOnly");
    h_PdgEta_HFHadron_prompt_PYTHIAOnly_CENTRALIZED->SetName(" h_PdgEta_HFHadron_prompt_PYTHIAOnly_CENTRALIZED");
    TH2F *h_PdgEta_HFHadron_prompt_GeantOnly_CENTRALIZED = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Hadron/Geant/h_PdgEta_HFHadron_prompt_GeantOnly");
    h_PdgEta_HFHadron_prompt_GeantOnly_CENTRALIZED->SetName("h_PdgEta_HFHadron_prompt_GeantOnly_CENTRALIZED");
    TH2F *h_PdgEta_HFHadron_prompt_STANDALONE_NoDiffr = (TH2F *)fIn_Pythia_STANDALONE_NoDiffr->Get("Hadron/h_PdgEta_HFHadron_prompt");
    h_PdgEta_HFHadron_prompt_STANDALONE_NoDiffr->SetName("h_PdgEta_HFHadron_prompt_STANDALONE_NoDiffr");

    TH2F *h_PdgEta_HFHadron_prompt_STANDALONE_EL = (TH2F *)fIn_Pythia_STANDALONE_EL->Get("Hadron/h_PdgEta_HFHadron_prompt");
    h_PdgEta_HFHadron_prompt_STANDALONE_EL->SetName("h_PdgEta_HFHadron_prompt_STANDALONE_EL");

    TH1F *h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED;
    TH1F *h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED;
    TH1F *h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr;
    TH1F *h_Eta_Hadron_Gen_LF_STANDALONE_EL;

    for (Int_t i_PDG = 0; i_PDG < n_PDG_tested; i_PDG++)
    {
        h_PdgEta_HFHadron_prompt_PYTHIAOnly_CENTRALIZED->GetXaxis()->SetRangeUser(low_PDG[i_PDG], high_PDG[i_PDG]);

        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED = (TH1F *)h_PdgEta_HFHadron_prompt_PYTHIAOnly_CENTRALIZED->ProjectionY();
        hist1D_graphic_opt(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED, kTRUE, 4, 20, kMagenta + 2, Norm_PYTHIA_MB_CENTRALIZED);
        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->SetName(Form("h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED_%s", canvas_name[i_PDG].Data()));
        std::cout << h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetEntries() << endl;
        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->SetTitle("PYTHIA8 centralized (PYTHIA only)");
        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("d#sigma/d#eta (mb)");

        h_PdgEta_HFHadron_prompt_GeantOnly_CENTRALIZED->GetXaxis()->SetRangeUser(low_PDG[i_PDG], high_PDG[i_PDG]);

        h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED = (TH1F *)h_PdgEta_HFHadron_prompt_GeantOnly_CENTRALIZED->ProjectionY();
        h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->SetName(Form("h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED_%s", canvas_name[i_PDG].Data()));
        h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->SetTitle("PYTHIA8 centralized (Geant only)");
        std::cout << h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->GetEntries() << endl;
        hist1D_graphic_opt(h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED, kTRUE, 4, 20, kGreen + 2, Norm_PYTHIA_MB_CENTRALIZED);

        //-----------------------------------------//

        h_PdgEta_HFHadron_prompt_STANDALONE_NoDiffr->SetTitle("PYTHIA8 stand-alone (SoftQCD ND)");
        h_PdgEta_HFHadron_prompt_STANDALONE_NoDiffr->GetXaxis()->SetRangeUser(low_PDG[i_PDG], high_PDG[i_PDG]);

        h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr = (TH1F *)h_PdgEta_HFHadron_prompt_STANDALONE_NoDiffr->ProjectionY();
        h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr->SetName(Form("h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr_%s", canvas_name[i_PDG].Data()));
        hist1D_graphic_opt(h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr, kTRUE, 4, 20, kCyan + 2, Norm_PYTHIA_STANDALONE_NoDiffr);

        //-----------------------------------------//

        h_PdgEta_HFHadron_prompt_STANDALONE_EL->SetTitle("PYTHIA8 stand-alone (SoftQCD INEL)");
        h_PdgEta_HFHadron_prompt_STANDALONE_EL->GetXaxis()->SetRangeUser(low_PDG[i_PDG], high_PDG[i_PDG]);

        h_Eta_Hadron_Gen_LF_STANDALONE_EL = (TH1F *)h_PdgEta_HFHadron_prompt_STANDALONE_EL->ProjectionY();
        h_Eta_Hadron_Gen_LF_STANDALONE_EL->SetName(Form("h_Eta_Hadron_Gen_LF_STANDALONE_EL_%s", canvas_name[i_PDG].Data()));
        std::cout << h_Eta_Hadron_Gen_LF_STANDALONE_EL->GetEntries() << endl;
        hist1D_graphic_opt(h_Eta_Hadron_Gen_LF_STANDALONE_EL, kTRUE, 4, 20, kOrange + 2, Norm_PYTHIA_STANDALONE_EL);

        //-----------------------------------------//

        // TH2F *h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm = (TH2F *)fIn_Powheg_Charm->Get("Hadron/not_Geant/h_PdgEta_HFHadron_prompt_PYTHIAOnly");
        // h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm->SetName(" h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm");
        // h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm->GetXaxis()->SetRangeUser(210, 220);

        // TH1F *h_Eta_Hadron_Gen_LF_PYTHIAOnly_Powheg_Charm = (TH1F *)h_PdgEta_HFHadron_prompt_PYTHIAOnly_Powheg_Charm->ProjectionY();
        // hist1D_graphic_opt(h_Eta_Hadron_Gen_LF_PYTHIAOnly_Powheg_Charm, kTRUE, 10, 20, kBlack, Norm_Powheg_Charm);

        TH1F *ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr = (TH1F *)h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->Clone("ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr");
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr->Divide(h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr);
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr->SetTitle("CENTRALIZED/STANDALONE ND");
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr->SetMarkerColor(h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr->GetMarkerColor());
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr->SetLineColor(h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr->GetLineColor());

        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetMinimum() * 0.8, h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetMaximum() * 2.8);
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr->GetYaxis()->CenterTitle();
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr->GetYaxis()->SetRangeUser(ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr->GetMinimum() * 0.8, ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr->GetMaximum() * 1.2);

        TH1F *ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL = (TH1F *)h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->Clone("ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL");
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL->Divide(h_Eta_Hadron_Gen_LF_STANDALONE_EL);
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL->SetTitle("CENTRALIZED/STANDALONE INEL");
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL->SetMarkerColor(h_Eta_Hadron_Gen_LF_STANDALONE_EL->GetMarkerColor());
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL->SetLineColor(h_Eta_Hadron_Gen_LF_STANDALONE_EL->GetLineColor());

        h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetMinimum() * 0.8, h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetMaximum() * 2.8);
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL->GetYaxis()->CenterTitle();
        ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL->GetYaxis()->SetRangeUser(ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL->GetMinimum() * 0.8, ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL->GetMaximum() * 1.2);

        TCanvas *c_ETA_GEN = three_histo_ratio(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED, h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr, h_Eta_Hadron_Gen_LF_STANDALONE_EL, ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr, ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_EL, canvas_name[i_PDG].Data(), canvas_header[i_PDG].Data(), kTRUE);

        // TCanvas *c_ETA_GEN = two_histo_ratio(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED, h_Eta_Hadron_Gen_LF_STANDALONE_NoDiffr, ratio_ETA_PYTHIA_CENTRALIZED_PYTHIA_STANDALONE_NoDiffr, canvas_name[i_PDG].Data(), canvas_header[i_PDG].Data(), kTRUE, kTRUE);
        c_ETA_GEN->SaveAs(Form("images/%s_withINEL.png", c_ETA_GEN->GetName()));

        // TH1F *ratio_ETA_Geant_PYTHIA_CENTRALIZED = (TH1F *)h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->Clone("ratio_ETA_Geant_PYTHIA_CENTRALIZED");
        // ratio_ETA_Geant_PYTHIA_CENTRALIZED->Divide(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED);
        // ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("Geant/PYTHIA");

        // h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED->GetMinimum() * 0.2, h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED->GetMaximum() * 400.8);
        // ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetYaxis()->CenterTitle();
        // ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetMinimum() * 0.3, ratio_ETA_Geant_PYTHIA_CENTRALIZED->GetMaximum() * 1.2);
        // TCanvas *c_ETA_GEANT_PYTHIA = two_histo_ratio(h_Eta_Hadron_Gen_LF_PYTHIA_CENTRALIZED, h_Eta_Hadron_Gen_GeantOnly_CENTRALIZED, ratio_ETA_Geant_PYTHIA_CENTRALIZED, Form("Geant_PYTHIA_%s", canvas_name[i_PDG].Data()), canvas_header[i_PDG].Data(), kTRUE, kTRUE);
        // c_ETA_GEANT_PYTHIA->SaveAs(Form("images/%s.png", c_ETA_GEANT_PYTHIA->GetName()));
    }
}

void VZ_mc_comparison_Mu_Gen()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    fIn_PYTHIA_CENTRALIZED_MB->cd("Muon_Gen");
    fIn_PYTHIA_CENTRALIZED_MB->ls();
    //------------CENTRALIZED SIM and comparison pythia only and geant only------------------//
    TH2F *h_VzmotherEta_Muon_Gen_Geant_LF_CENTRALUZED = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Gen/h_VzmotherEta_Muon_Gen_Geant_LF");
    h_VzmotherEta_Muon_Gen_Geant_LF_CENTRALUZED->GetYaxis()->SetRangeUser(-8, 8);
    TH1F *h_Eta_Muon_Gen_Geant_LF_CENTRALIZED = (TH1F *)h_VzmotherEta_Muon_Gen_Geant_LF_CENTRALUZED->ProjectionY();
    h_Eta_Muon_Gen_Geant_LF_CENTRALIZED->SetTitle("from GEANT hadrons");
    h_Eta_Muon_Gen_Geant_LF_CENTRALIZED->GetYaxis()->SetTitle("d#sigma/d#eta (mb)");
    hist1D_graphic_opt(h_Eta_Muon_Gen_Geant_LF_CENTRALIZED, kTRUE, 16, 20, kGreen + 2, Norm_PYTHIA_MB_CENTRALIZED);
    TH1F *h_Vzmother_Muon_Gen_Geant_LF_CENTRALIZED = (TH1F *)h_VzmotherEta_Muon_Gen_Geant_LF_CENTRALUZED->ProjectionX();
    h_Vzmother_Muon_Gen_Geant_LF_CENTRALIZED->SetTitle("from GEANT hadron");
    hist1D_graphic_opt(h_Vzmother_Muon_Gen_Geant_LF_CENTRALIZED, kTRUE, 50, 20, kGreen + 2, 1. / Nev_PYTHIA_MB_CENTRALIZED);

    TH2F *h_VzmotherEta_Muon_Gen_PYTHIA_LF_CENTRALIZED = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Gen/h_VzmotherEta_Muon_Gen_PYTHIA_LF");
    h_VzmotherEta_Muon_Gen_PYTHIA_LF_CENTRALIZED->GetYaxis()->SetRangeUser(-8, 8);

    TH1F *h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED = (TH1F *)h_VzmotherEta_Muon_Gen_PYTHIA_LF_CENTRALIZED->ProjectionY();
    h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED->SetTitle("from PYTHIA hadrons");
    cout << h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED->GetNbinsX() << endl;

    hist1D_graphic_opt(h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED, kTRUE, 16, 20, kMagenta + 2, Norm_PYTHIA_MB_CENTRALIZED);
    cout << h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED->GetNbinsX() << endl;

    TH1F *h_Vzmother_Muon_Gen_PYTHIA_LF = (TH1F *)h_VzmotherEta_Muon_Gen_PYTHIA_LF_CENTRALIZED->ProjectionX();
    h_Vzmother_Muon_Gen_PYTHIA_LF->SetTitle("from PYTHIA hadron");
    hist1D_graphic_opt(h_Vzmother_Muon_Gen_PYTHIA_LF, kTRUE, 5, 20, kMagenta + 2, 1. / Nev_PYTHIA_MB_CENTRALIZED);
    h_Vzmother_Muon_Gen_PYTHIA_LF->GetXaxis()->SetRangeUser(-10, 10);

    TH1F *ratio_ETA_GEANT_PYTHIA_CENTRALIZED = (TH1F *)h_Eta_Muon_Gen_Geant_LF_CENTRALIZED->Clone("ratio_ETA_GEANT_PYTHIA_CENTRALIZED");
    ratio_ETA_GEANT_PYTHIA_CENTRALIZED->Divide(h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED);
    ratio_ETA_GEANT_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("GEANT/PYTHIA");

    h_Eta_Muon_Gen_Geant_LF_CENTRALIZED->GetYaxis()->SetRangeUser(h_Eta_Muon_Gen_Geant_LF_CENTRALIZED->GetMinimum() * 0.6, h_Eta_Muon_Gen_Geant_LF_CENTRALIZED->GetMaximum() * 423.6);
    ratio_ETA_GEANT_PYTHIA_CENTRALIZED->GetYaxis()->CenterTitle();
    ratio_ETA_GEANT_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(ratio_ETA_GEANT_PYTHIA_CENTRALIZED->GetMinimum() * 0.2, ratio_ETA_GEANT_PYTHIA_CENTRALIZED->GetMaximum() * 1.2);
    TCanvas *c_ETA_GEN = two_histo_ratio(h_Eta_Muon_Gen_Geant_LF_CENTRALIZED, h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED, ratio_ETA_GEANT_PYTHIA_CENTRALIZED, "c_ETA_PYTHIA_GEANT", "Gen #mu #leftarrow LF hadrons", kTRUE, kTRUE);
    c_ETA_GEN->SaveAs(Form("images/%s_MuonGen.png", c_ETA_GEN->GetName()));

    TCanvas *c_Vzmother = canvas_noratio_divide2("c_Vzmother");
    c_Vzmother->cd();
    c_Vzmother->GetPad(1)->cd();
    h_Vzmother_Muon_Gen_Geant_LF_CENTRALIZED->Draw();
    c_Vzmother->GetPad(2)->cd();
    h_Vzmother_Muon_Gen_PYTHIA_LF->Draw();

    //------------STANDALONE SIm------------------//

    TH2F *h_PtEta_Muon_Gen_LF_STANDALONE = (TH2F *)fIn_Pythia_STANDALONE_NoDiffr->Get("Muon_Gen/h_PtEta_Muon_Gen_LF");
    TH1F *h_Eta_Muon_Gen_LF_STANDALONE = (TH1F *)h_PtEta_Muon_Gen_LF_STANDALONE->ProjectionY();
    cout << h_Eta_Muon_Gen_LF_STANDALONE->GetNbinsX() << endl;
    hist1D_graphic_opt(h_Eta_Muon_Gen_LF_STANDALONE, kTRUE, 90, 20, kCyan + 2, Norm_PYTHIA_STANDALONE_NoDiffr);

    cout << h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED->GetEntries() << endl;
    cout << h_Eta_Muon_Gen_LF_STANDALONE->GetEntries() << endl;
    TH1F *ratio_ETA_PYTHIA_STANDALONE_CENTRALIZED = (TH1F *)h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED->Clone("ratio_ETA_PYTHIA_STANDALONE_CENTRALIZED");
    ratio_ETA_PYTHIA_STANDALONE_CENTRALIZED->Divide(h_Eta_Muon_Gen_LF_STANDALONE);
    ratio_ETA_PYTHIA_STANDALONE_CENTRALIZED->GetYaxis()->SetTitle("CENTRALIZED/STAND-ALONE");

    h_Eta_Muon_Gen_LF_STANDALONE->SetTitle("stand-alone");
    h_Eta_Muon_Gen_LF_STANDALONE->GetYaxis()->SetTitle("d#sigma/d#eta (mb)");
    h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED->SetTitle("centralized (PYTHIA only)");
    h_Eta_Muon_Gen_LF_STANDALONE->GetYaxis()->SetRangeUser(h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED->GetMinimum() * 0.6, h_Eta_Muon_Gen_LF_STANDALONE->GetMaximum() * 230.6);
    ratio_ETA_PYTHIA_STANDALONE_CENTRALIZED->GetYaxis()->CenterTitle();
    ratio_ETA_PYTHIA_STANDALONE_CENTRALIZED->GetYaxis()->SetRangeUser(ratio_ETA_PYTHIA_STANDALONE_CENTRALIZED->GetMinimum() * 0.2, ratio_ETA_PYTHIA_STANDALONE_CENTRALIZED->GetMaximum() * 1.2);
    TCanvas *c_ETA_CENTRALIZED_STANDALONE = two_histo_ratio(h_Eta_Muon_Gen_LF_STANDALONE, h_Eta_Muon_Gen_PYTHIA_LF_CENTRALIZED, ratio_ETA_PYTHIA_STANDALONE_CENTRALIZED, "c_ETA_CENTRALIZED_STANDALONE", "Gen #mu #leftarrow LF hadrons", kTRUE, kTRUE);
    c_ETA_CENTRALIZED_STANDALONE->SaveAs(Form("images/%s_MuonGen.png", c_ETA_CENTRALIZED_STANDALONE->GetName()));
}

void mc_comparinson_Rec_Mu()

{

    TH1F *h_Nevents_PYTHIA = (TH1F *)fIn_PYTHIA_CENTRALIZED_MB->Get("h_Nevents");
    Int_t Nev_PYTHIA = (Int_t)h_Nevents_PYTHIA->GetBinContent(2);
    Double_t Norm_PYTHIA = (1. / ((Double_t)Nev_PYTHIA * 216)); // in mb

    fIn_PYTHIA_CENTRALIZED_MB->cd("Muon_Rec");
    fIn_PYTHIA_CENTRALIZED_MB->ls();

    TH2F *h_VzmotherY_Muon_Rec_Geant_LF = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/h_VzmotherY_Muon_Rec_Geant_LF");
    TH1F *h_Vzmother_Muon_Rec_Geant_LF = (TH1F *)h_VzmotherY_Muon_Rec_Geant_LF->ProjectionX();
    hist1D_graphic_opt(h_Vzmother_Muon_Rec_Geant_LF, kTRUE, 1, 20, kMagenta + 2, 1.);
    // h_Vzmother_Muon_Rec_Geant_LF->GetXaxis()->SetRangeUser(-500, 200);
    TH1F *h_Y_Muon_Rec_Geant_LF = (TH1F *)h_VzmotherY_Muon_Rec_Geant_LF->ProjectionY();
    hist1D_graphic_opt(h_Y_Muon_Rec_Geant_LF, kTRUE, 1, 20, kMagenta + 2, 1. / h_Y_Muon_Rec_Geant_LF->GetEntries());
    h_Y_Muon_Rec_Geant_LF->SetTitle("from GEANT #pi");

    TH2F *h_VzmotherY_Muon_Rec_PYTHIA_LF = (TH2F *)fIn_PYTHIA_CENTRALIZED_MB->Get("Muon_Rec/h_VzmotherY_Muon_Rec_PYTHIA_LF");
    TH1F *h_Vzmother_Muon_Rec_PYTHIA_LF = (TH1F *)h_VzmotherY_Muon_Rec_PYTHIA_LF->ProjectionX();
    hist1D_graphic_opt(h_Vzmother_Muon_Rec_PYTHIA_LF, kTRUE, 1, 20, kGreen + 2, 1.);
    // h_Vzmother_Muon_Rec_PYTHIA_LF->GetXaxis()->SetRangeUser(-30, 30);
    TH1F *h_Y_Muon_Rec_PYTHIA_LF = (TH1F *)h_VzmotherY_Muon_Rec_PYTHIA_LF->ProjectionY();
    hist1D_graphic_opt(h_Y_Muon_Rec_PYTHIA_LF, kTRUE, 1, 20, kGreen + 2, 1. / h_Y_Muon_Rec_PYTHIA_LF->GetEntries());
    h_Y_Muon_Rec_PYTHIA_LF->SetTitle("from PYTHIA #pi");

    TCanvas *c_Vzmother = canvas_noratio_divide2("c_Vzmother");
    c_Vzmother->cd();
    c_Vzmother->GetPad(1)->cd();
    h_Vzmother_Muon_Rec_Geant_LF->Draw();
    c_Vzmother->GetPad(2)->cd();
    h_Vzmother_Muon_Rec_PYTHIA_LF->Draw();

    TH1F *ratio_histo_fit = (TH1F *)h_Y_Muon_Rec_Geant_LF->Clone("ratio_histo_fit");
    ratio_histo_fit->Divide(h_Y_Muon_Rec_PYTHIA_LF);
    ratio_histo_fit->SetTitle(" ");
    ratio_histo_fit->GetYaxis()->SetTitle("Geant/PYTHIA");

    TCanvas *c_Y = two_histo_ratio(h_Y_Muon_Rec_Geant_LF, h_Y_Muon_Rec_PYTHIA_LF, ratio_histo_fit, "c_Y", " ", kTRUE, kTRUE);
}

void mc_comparinson_Rec_Dimu()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    Int_t rebin_factor = 10;
    //============================POWHEG CHARM=======================================//

    TH1F *h_Nevents_Powheg_Charm = (TH1F *)fIn_Powheg_Charm->Get("h_Nevents");
    h_Nevents_Powheg_Charm->SetName("h_Nevents_Powheg_Charm");
    cout << h_Nevents_Powheg_Charm->GetBinContent(2) << endl;

    TH2F *h_PtM_DiMuon_Rec_Charm = (TH2F *)fIn_Powheg_Charm->Get("DiMuon_Rec/PowhegOnly/h_PtM_DiMuon_Rec_Charm_PowhegOnly");
    h_PtM_DiMuon_Rec_Charm->SetName("h_PtM_DiMuon_Rec_Charm");
    h_PtM_DiMuon_Rec_Charm->GetXaxis()->SetRangeUser(0, 30);
    h_PtM_DiMuon_Rec_Charm->GetYaxis()->SetRangeUser(0, 30);

    // Double_t Powheg_Charm_Norm = (5 / 56.42) / (h_Nevents_Powheg_Charm->GetBinContent(2) * 30.5);

    TH1F *h_Pt_DiMuon_Rec_Charm = (TH1F *)h_PtM_DiMuon_Rec_Charm->ProjectionX();
    hist1D_graphic_opt(h_Pt_DiMuon_Rec_Charm, kTRUE, rebin_factor, 20, kMagenta + 2, 1. / h_Pt_DiMuon_Rec_Charm->Integral());
    // h_Pt_DiMuon_Rec_Charm->SetTitle(Form("POWHEG+PYTHIA6 %0.2e", h_Pt_DiMuon_Rec_Charm->Integral()));
    h_Pt_DiMuon_Rec_Charm->SetTitle(Form("POWHEG+PYTHIA6"));

    TH1F *h_M_DiMuon_Rec_Charm = (TH1F *)h_PtM_DiMuon_Rec_Charm->ProjectionY();
    hist1D_graphic_opt(h_M_DiMuon_Rec_Charm, kTRUE, rebin_factor, 20, kMagenta + 2, 1. / h_M_DiMuon_Rec_Charm->Integral());
    // h_M_DiMuon_Rec_Charm->SetTitle(Form("POWHEG+PYTHIA6"));
    h_M_DiMuon_Rec_Charm->SetTitle(Form("POWHEG+PYTHIA6"));

    //============================POWHEG Beauty=======================================//

    TH1F *h_Nevents_Powheg_Beauty = (TH1F *)fIn_Powheg_Beauty->Get("h_Nevents");
    h_Nevents_Powheg_Beauty->SetName("h_Nevents_Powheg_Beauty");
    cout << h_Nevents_Powheg_Beauty->GetBinContent(2) << endl;

    TH2F *h_PtM_DiMuon_Rec_Beauty = (TH2F *)fIn_Powheg_Beauty->Get("DiMuon_Rec/PowhegOnly/h_PtM_DiMuon_Rec_Beauty_PowhegOnly");
    h_PtM_DiMuon_Rec_Beauty->GetXaxis()->SetRangeUser(0, 30);
    h_PtM_DiMuon_Rec_Beauty->GetYaxis()->SetRangeUser(0, 30);
    h_PtM_DiMuon_Rec_Beauty->SetName("h_PtM_DiMuon_Rec_Beauty");

    // Double_t Powheg_Beauty_Norm = (0.5 / 56.42) / (h_Nevents_Powheg_Beauty->GetBinContent(2) * 15);

    TH1F *h_Pt_DiMuon_Rec_Beauty = (TH1F *)h_PtM_DiMuon_Rec_Beauty->ProjectionX();
    hist1D_graphic_opt(h_Pt_DiMuon_Rec_Beauty, kTRUE, rebin_factor, 20, kGreen + 2, 1. / h_Pt_DiMuon_Rec_Beauty->Integral());
    // h_Pt_DiMuon_Rec_Beauty->SetTitle(Form("POWHEG+PYTHIA6 %0.2e", h_Pt_DiMuon_Rec_Beauty->Integral()));
    h_Pt_DiMuon_Rec_Beauty->SetTitle(Form("POWHEG+PYTHIA6"));

    TH1F *h_M_DiMuon_Rec_Beauty = (TH1F *)h_PtM_DiMuon_Rec_Beauty->ProjectionY();
    hist1D_graphic_opt(h_M_DiMuon_Rec_Beauty, kTRUE, rebin_factor, 20, kGreen + 2, 1. / h_M_DiMuon_Rec_Beauty->Integral());
    // h_M_DiMuon_Rec_Beauty->SetTitle(Form("POWHEG+PYTHIA6"));
    h_M_DiMuon_Rec_Beauty->SetTitle(Form("POWHEG+PYTHIA6"));

    //============================Pythia=======================================//

    TH1F *h_Nevents_PYTHIA = (TH1F *)fIn_PYTHIA_HF->Get("h_Nevents");
    Int_t Nev_PYTHIA = (Int_t)h_Nevents_PYTHIA->GetBinContent(2);
    // Double_t Norm_PYTHIA = (1. / ((Double_t)Nev_PYTHIA * 216)); // in mb

    TH2F *h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm = (TH2F *)fIn_PYTHIA_HF->Get("DiMuon/M0/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_M0_Rec_DQ_cut_match_LT_ULS_fromCharm");

    TH1F *h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm = (TH1F *)h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->ProjectionX();
    hist1D_graphic_opt(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm, kTRUE, rebin_factor, 24, kMagenta + 2, 1. / h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->Integral());
    // h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->SetTitle(Form("PYTHIA8 %0.2e", h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->Integral()));
    h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->SetTitle(Form("PYTHIA8"));

    TH1F *h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm = (TH1F *)h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->ProjectionY();
    hist1D_graphic_opt(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm, kTRUE, rebin_factor, 24, kMagenta + 2, 1. / h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->Integral());
    // h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->SetTitle(Form("PYTHIA8 %0.2e", h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->Integral()));
    h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->SetTitle(Form("PYTHIA8"));

    TH2F *h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty = (TH2F *)fIn_PYTHIA_HF->Get("DiMuon/M0/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_M0_Rec_DQ_cut_match_LT_ULS_fromBeauty");

    TH1F *h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty = (TH1F *)h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->ProjectionX();
    hist1D_graphic_opt(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty, kTRUE, rebin_factor, 24, kGreen + 2, 1. / h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->Integral());
    // h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->SetTitle(Form("PYTHIA8 %0.2e", h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->Integral()));
    h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->SetTitle(Form("PYTHIA8"));

    TH1F *h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty = (TH1F *)h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->ProjectionY();
    hist1D_graphic_opt(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty, kTRUE, rebin_factor, 24, kGreen + 2, 1. / h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->Integral());
    // h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->SetTitle(Form("PYTHIA8 %0.2e", h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->Integral()));
    h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->SetTitle(Form("PYTHIA8"));

    // h_Pt_DiMuon_Rec_Charm->Divide(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm);
    // h_Pt_DiMuon_Rec_Charm->Draw();

    TH1F *h_PtRatio_Powheg_Pythia_Charm = (TH1F *)h_Pt_DiMuon_Rec_Charm->Clone("h_PtRatio_Powheg_Pythia_Charm");
    h_PtRatio_Powheg_Pythia_Charm->Divide(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm);
    h_PtRatio_Powheg_Pythia_Charm->GetYaxis()->SetTitle("Powheg/Pythia");
    // h_Pt_DiMuon_Rec_Charm->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    h_Pt_DiMuon_Rec_Charm->GetYaxis()->SetTitle("a.u.");
    TCanvas *c_Pt_Charm = two_histo_ratio(h_Pt_DiMuon_Rec_Charm, h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm, h_PtRatio_Powheg_Pythia_Charm, "c_Pt_Charm", "ALICE Simulation, Rec. #mu#mu #leftarrow c", kTRUE, kTRUE);
    c_Pt_Charm->SaveAs(Form("images/%s_simfigure.png", c_Pt_Charm->GetName()));

    TH1F *h_MRatio_Powheg_Pythia_Charm = (TH1F *)h_M_DiMuon_Rec_Charm->Clone("h_MRatio_Powheg_Pythia_Charm");
    h_MRatio_Powheg_Pythia_Charm->Divide(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm);
    h_MRatio_Powheg_Pythia_Charm->GetYaxis()->SetTitle("Powheg/Pythia");

    // h_M_DiMuon_Rec_Charm->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
    h_M_DiMuon_Rec_Charm->GetYaxis()->SetTitle("a.u.");
    TCanvas *c_M_Charm = two_histo_ratio(h_M_DiMuon_Rec_Charm, h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm, h_MRatio_Powheg_Pythia_Charm, "c_M_Charm", "ALICE Simulation, Rec. #mu#mu #leftarrow c", kTRUE, kTRUE);
    c_M_Charm->SaveAs(Form("images/%s_simfigure.png", c_M_Charm->GetName()));

    TH1F *h_PtRatio_Powheg_Pythia_Beauty = (TH1F *)h_Pt_DiMuon_Rec_Beauty->Clone("h_PtRatio_Powheg_Pythia_Beauty");
    h_PtRatio_Powheg_Pythia_Beauty->Divide(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty);
    h_PtRatio_Powheg_Pythia_Beauty->GetYaxis()->SetTitle("Powheg/Pythia");
    // h_Pt_DiMuon_Rec_Beauty->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    h_Pt_DiMuon_Rec_Beauty->GetYaxis()->SetTitle("a.u.");
    TCanvas *c_Pt_Beauty = two_histo_ratio(h_Pt_DiMuon_Rec_Beauty, h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty, h_PtRatio_Powheg_Pythia_Beauty, "c_Pt_Beauty", "ALICE Simulation, Rec. #mu#mu #leftarrow b", kTRUE, kTRUE);
    c_Pt_Beauty->SaveAs(Form("images/%s_simfigure.png", c_Pt_Beauty->GetName()));

    TH1F *h_MRatio_Powheg_Pythia_Beauty = (TH1F *)h_M_DiMuon_Rec_Beauty->Clone("h_MRatio_Powheg_Pythia_Beauty");
    h_MRatio_Powheg_Pythia_Beauty->Divide(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty);
    h_MRatio_Powheg_Pythia_Beauty->GetYaxis()->SetTitle("Powheg/Pythia");
    // h_M_DiMuon_Rec_Beauty->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
    h_M_DiMuon_Rec_Beauty->GetYaxis()->SetTitle("a.u.");
    TCanvas *c_M_Beauty = two_histo_ratio(h_M_DiMuon_Rec_Beauty, h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty, h_MRatio_Powheg_Pythia_Beauty, "c_M_Beauty", "ALICE Simulation, Rec. #mu#mu #leftarrow b", kTRUE, kTRUE);
    c_M_Beauty->SaveAs(Form("images/%s_simfigure.png", c_M_Beauty->GetName()));
    // c_Pt_Charm->GetPad(0)->SetGridx();
    // c_Pt_Charm->GetPad(0)->SetGridy();

    printf("Powheg charm: %0.3f\n", h_Pt_DiMuon_Rec_Charm->Integral());
    printf("Pythia charm: %0.3f\n", h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->Integral());
    // h_Pt_DiMuon_Rec_Beauty->Draw("PE");
    // h_Pt_DiMuon_Rec_Charm->Draw("PESAME");
    // h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Charm->Draw("PESAME");
    // h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA_Beauty->Draw("PESAME");
}

void mc_comparinson_LF()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    //============================PYTHIA STANDALONE=======================================//

    fIn_Pythia_STANDALONE_NoDiffr->cd("HF_Hadron");
    fIn_Pythia_STANDALONE_NoDiffr->ls();

    TH1F *h_Nevents_Pythia_STANDALONE = (TH1F *)fIn_Pythia_STANDALONE_NoDiffr->Get("h_Nevents");
    h_Nevents_Pythia_STANDALONE->SetName("h_Nevents_Pythia_STANDALONE");

    TH3F *h_PdgPtY_HFHadron_prompt_PYTHIA_STANDALONE = (TH3F *)fIn_Pythia_STANDALONE_NoDiffr->Get("HF_Hadron/h_PdgPtY_HFHadron_prompt");
    h_PdgPtY_HFHadron_prompt_PYTHIA_STANDALONE->SetName("h_PdgPtY_HFHadron_prompt_PYTHIA_STANDALONE");
    h_PdgPtY_HFHadron_prompt_PYTHIA_STANDALONE->GetXaxis()->SetRangeUser(200, 400);
    // h_PdgPtY_HFHadron_prompt_PYTHIA_STANDALONE->GetZaxis()->SetRangeUser(-4, -2.5);

    TH2F *h_PtY_HFHadron_prompt_PYTHIA_STANDALONE = (TH2F *)h_PdgPtY_HFHadron_prompt_PYTHIA_STANDALONE->Project3D("yze");
    h_PtY_HFHadron_prompt_PYTHIA_STANDALONE->Scale(1. / h_Nevents_Pythia_STANDALONE->GetBinContent(2));

    TH1F *h_Pt_HFHadron_prompt_PYTHIA_STANDALONE = (TH1F *)h_PdgPtY_HFHadron_prompt_PYTHIA_STANDALONE->Project3D("ye");
    h_Pt_HFHadron_prompt_PYTHIA_STANDALONE->SetTitle("PYTHIA8 Monash stand-alone");
    h_Pt_HFHadron_prompt_PYTHIA_STANDALONE->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    hist1D_graphic_opt(h_Pt_HFHadron_prompt_PYTHIA_STANDALONE, kTRUE, 10, 20, kBlack, 1. / h_Nevents_Pythia_STANDALONE->GetBinContent(2));
    h_Pt_HFHadron_prompt_PYTHIA_STANDALONE->GetYaxis()->SetRangeUser(6e-07, 1e+04);

    TH1F *h_Y_HFHadron_prompt_PYTHIA_STANDALONE = (TH1F *)h_PdgPtY_HFHadron_prompt_PYTHIA_STANDALONE->Project3D("ze");
    hist1D_graphic_opt(h_Y_HFHadron_prompt_PYTHIA_STANDALONE, kTRUE, 10, 20, kBlack, 1. / h_Nevents_Pythia_STANDALONE->GetBinContent(2));
    h_Y_HFHadron_prompt_PYTHIA_STANDALONE->SetTitle("PYTHIA8 Monash stand-alone");
    h_Y_HFHadron_prompt_PYTHIA_STANDALONE->GetYaxis()->SetTitle("d#it{N}/d#it{y}");
    h_Y_HFHadron_prompt_PYTHIA_STANDALONE->GetYaxis()->SetRangeUser(8.4e-02, 1e+03);

    //============================PYTHIA CENTRALIZED=======================================//
    TFile *fIn_Pythia_CENTRALIZED = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC22c1/output/Merged_LHC22c1_MC_output_Hist_294009.root", "READ");

    TH1F *h_Nevents_Pythia_CENTRALIZED = (TH1F *)fIn_Pythia_CENTRALIZED->Get("h_Nevents");
    h_Nevents_Pythia_CENTRALIZED->SetName("h_Nevents_Pythia_CENTRALIZED");

    TH3F *h_PdgPtY_HFHadron_prompt_PYTHIA_CENTRALIZED = (TH3F *)fIn_Pythia_CENTRALIZED->Get("HF_Hadron/h_PdgPtY_HFHadron_prompt");
    h_PdgPtY_HFHadron_prompt_PYTHIA_CENTRALIZED->SetName("h_PdgPtY_HFHadron_prompt_PYTHIA_CENTRALIZED");
    h_PdgPtY_HFHadron_prompt_PYTHIA_CENTRALIZED->GetXaxis()->SetRangeUser(200, 400);
    // h_PdgPtY_HFHadron_prompt_PYTHIA_CENTRALIZED->GetZaxis()->SetRangeUser(-4, -2.5);

    TH2F *h_PtY_HFHadron_prompt_PYTHIA_CENTRALIZED = (TH2F *)h_PdgPtY_HFHadron_prompt_PYTHIA_CENTRALIZED->Project3D("yze");
    h_PtY_HFHadron_prompt_PYTHIA_CENTRALIZED->Scale(1. / h_Nevents_Pythia_CENTRALIZED->GetBinContent(2));

    TH1F *h_Pt_HFHadron_prompt_PYTHIA_CENTRALIZED = (TH1F *)h_PdgPtY_HFHadron_prompt_PYTHIA_CENTRALIZED->Project3D("ye");
    h_Pt_HFHadron_prompt_PYTHIA_CENTRALIZED->SetTitle("PYTHIA8 Monash centralized");
    h_Pt_HFHadron_prompt_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    hist1D_graphic_opt(h_Pt_HFHadron_prompt_PYTHIA_CENTRALIZED, kTRUE, 10, 20, kBlue, 1. / h_Nevents_Pythia_CENTRALIZED->GetBinContent(2));
    h_Pt_HFHadron_prompt_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(6e-07, 1e+04);

    TH1F *h_Y_HFHadron_prompt_PYTHIA_CENTRALIZED = (TH1F *)h_PdgPtY_HFHadron_prompt_PYTHIA_CENTRALIZED->Project3D("ze");
    hist1D_graphic_opt(h_Y_HFHadron_prompt_PYTHIA_CENTRALIZED, kTRUE, 10, 20, kBlue, 1. / h_Nevents_Pythia_CENTRALIZED->GetBinContent(2));
    h_Y_HFHadron_prompt_PYTHIA_CENTRALIZED->SetTitle("PYTHIA8 Monash centralized");
    h_Y_HFHadron_prompt_PYTHIA_CENTRALIZED->GetYaxis()->SetTitle("d#it{N}/d#it{y}");
    h_Y_HFHadron_prompt_PYTHIA_CENTRALIZED->GetYaxis()->SetRangeUser(8.4e-02, 1e+03);

    //============================POWHEG CHARM=======================================//
    TFile *fIn_Powheg_Charm = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC23i1_Version5_AliAOD_HF_LF/LHC23i1_MC_output_Hist_294009.root", "READ");

    TH1F *h_Nevents_Powheg_Charm = (TH1F *)fIn_Powheg_Charm->Get("h_Nevents");
    h_Nevents_Powheg_Charm->SetName("h_Nevents_Powheg_Charm");

    TH3F *h_PdgPtY_HFHadron_prompt_Powheg_Charm = (TH3F *)fIn_Powheg_Charm->Get("HF_Hadron/h_PdgPtY_HFHadron_prompt_PYTHIAOnly");
    h_PdgPtY_HFHadron_prompt_Powheg_Charm->SetName("h_PdgPtY_HFHadron_prompt_Powheg_Charm");
    h_PdgPtY_HFHadron_prompt_Powheg_Charm->GetXaxis()->SetRangeUser(200, 400);
    // h_PdgPtY_HFHadron_prompt_Powheg_Charm->GetZaxis()->SetRangeUser(-4, -2.5);

    TH2F *h_PtY_HFHadron_prompt_Powheg_Charm = (TH2F *)h_PdgPtY_HFHadron_prompt_Powheg_Charm->Project3D("yze");
    h_PtY_HFHadron_prompt_Powheg_Charm->Scale(1. / h_Nevents_Powheg_Charm->GetBinContent(2));

    TH1F *h_Pt_HFHadron_prompt_Powheg_Charm = (TH1F *)h_PdgPtY_HFHadron_prompt_Powheg_Charm->Project3D("ye");
    hist1D_graphic_opt(h_Pt_HFHadron_prompt_Powheg_Charm, kTRUE, 10, 20, kRed, 1. / h_Nevents_Powheg_Charm->GetBinContent(2));
    h_Pt_HFHadron_prompt_Powheg_Charm->SetTitle("POWHEG+PYTHIA6");

    TH1F *h_Y_HFHadron_prompt_Powheg_Charm = (TH1F *)h_PdgPtY_HFHadron_prompt_Powheg_Charm->Project3D("ze");
    hist1D_graphic_opt(h_Y_HFHadron_prompt_Powheg_Charm, kTRUE, 10, 20, kRed, 1. / h_Nevents_Powheg_Charm->GetBinContent(2));
    h_Y_HFHadron_prompt_Powheg_Charm->SetTitle("POWHEG+PYTHIA6");

    //============================DRAWING=======================================//

    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    c->Divide(2, 1);
    c->cd(1);
    h_PtY_HFHadron_prompt_PYTHIA_STANDALONE->Draw("COLZ");

    c->cd(2);
    h_PtY_HFHadron_prompt_Powheg_Charm->Draw("COLZ");

    TH1F *h_Pt_ratio_Powheg_Pythia_STANDALONE = (TH1F *)h_Pt_HFHadron_prompt_Powheg_Charm->Clone("h_Pt_ratio_Powheg_Pythia_STANDALONE");
    h_Pt_ratio_Powheg_Pythia_STANDALONE->Divide(h_Pt_HFHadron_prompt_PYTHIA_STANDALONE);
    h_Pt_ratio_Powheg_Pythia_STANDALONE->SetTitle("POWHEG/PYTHIA stand-alone");
    h_Pt_ratio_Powheg_Pythia_STANDALONE->GetYaxis()->SetRangeUser(-0.2, 3.2);

    TH1F *h_Y_ratio_Powheg_Pythia_STANDALONE = (TH1F *)h_Y_HFHadron_prompt_Powheg_Charm->Clone("h_Y_ratio_Powheg_Pythia_STANDALONE");
    h_Y_ratio_Powheg_Pythia_STANDALONE->Divide(h_Y_HFHadron_prompt_PYTHIA_STANDALONE);
    h_Y_ratio_Powheg_Pythia_STANDALONE->SetTitle("POWHEG/PYTHIA standa-alone");
    h_Y_ratio_Powheg_Pythia_STANDALONE->GetYaxis()->SetRangeUser(-0.2, 3.2);

    TH1F *h_Pt_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE = (TH1F *)h_Pt_HFHadron_prompt_PYTHIA_CENTRALIZED->Clone("h_Pt_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE");
    h_Pt_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE->Divide(h_Pt_HFHadron_prompt_PYTHIA_STANDALONE);
    h_Pt_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE->SetTitle("PYTHIA centralized/PYTHIA stand-alone");
    h_Pt_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE->GetYaxis()->SetRangeUser(-0.2, 3.2);

    TH1F *h_Y_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE = (TH1F *)h_Y_HFHadron_prompt_PYTHIA_CENTRALIZED->Clone("h_Y_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE");
    h_Y_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE->Divide(h_Y_HFHadron_prompt_PYTHIA_STANDALONE);
    h_Y_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE->SetTitle("PYTHIA centralized/PYTHIA stand-alone");
    h_Y_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE->GetYaxis()->SetRangeUser(-0.2, 3.2);

    TCanvas *c_Pt_Charm = three_histo_ratio(h_Pt_HFHadron_prompt_PYTHIA_STANDALONE, h_Pt_HFHadron_prompt_Powheg_Charm, h_Pt_HFHadron_prompt_PYTHIA_CENTRALIZED, h_Pt_ratio_Powheg_Pythia_STANDALONE, h_Pt_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE, "c_Pt_Charm", "ALICE Simulation, LF hadron", kTRUE);
    c_Pt_Charm->GetPad(0)->SetGridx();
    c_Pt_Charm->GetPad(0)->SetGridy();
    // c_Pt_Charm->SaveAs("images/mc_comparison/pt_LF_hadrons.png");

    TCanvas *c_Y_Charm = three_histo_ratio(h_Y_HFHadron_prompt_PYTHIA_STANDALONE, h_Y_HFHadron_prompt_Powheg_Charm, h_Y_HFHadron_prompt_PYTHIA_CENTRALIZED, h_Y_ratio_Powheg_Pythia_STANDALONE, h_Y_ratio_Pythia_CENTRALIZED_Pythia_STANDALONE, "c_Y_Charm", "ALICE Simulation, LF hadron", kTRUE);
    c_Y_Charm->GetPad(0)->SetGridx();
    c_Y_Charm->GetPad(0)->SetGridy();
    // c_Y_Charm->SaveAs("images/mc_comparison/pt_LF_hadrons.png");

    // TCanvas *c_Pt_Charm = two_histo_ratio(h_Pt_HFHadron_prompt_PYTHIA_STANDALONE, h_Pt_HFHadron_prompt_Powheg_Charm, h_Pt_ratio_Powheg_Pythia_STANDALONE, "c_Pt_Charm", "ALICE Simulation, LF hadron", kTRUE, kTRUE);
    // c_Pt_Charm->GetPad(0)->SetGridx();
    // c_Pt_Charm->GetPad(0)->SetGridy();
    // c_Pt_Charm->SaveAs("images/mc_comparison/pt_LF_hadrons.png");

    // TCanvas *c_Y_Charm = two_histo_ratio(h_Y_HFHadron_prompt_PYTHIA_STANDALONE, h_Y_HFHadron_prompt_Powheg_Charm, h_Y_ratio_Powheg_Pythia_STANDALONE, "c_Y_Charm", "ALICE Simulation, LF hadron", kTRUE, kTRUE);
    // c_Y_Charm->GetPad(0)->SetGridx();
    // c_Y_Charm->GetPad(0)->SetGridy();
    // c_Y_Charm->SaveAs("images/mc_comparison/y_LF_hadrons.png");
}

void mc_comparinson_hfquark()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TH1F *h_NCharm_event_Pythia = (TH1F *)fIn_Pythia_STANDALONE_NoDiffr->Get("HF_quarks/h_NCharm_event");

    TH2F *h_PtY_Charm_quark_Pythia = (TH2F *)fIn_Pythia_STANDALONE_NoDiffr->Get("HF_quarks/h_PtY_Charm_quark");
    h_PtY_Charm_quark_Pythia->SetName("h_PtY_Charm_quark_Pythia");

    TH1F *h_Y_Charm_quark_Pythia = (TH1F *)h_PtY_Charm_quark_Pythia->ProjectionY();
    h_Y_Charm_quark_Pythia->GetXaxis()->SetTitle("#it{y}");
    h_Y_Charm_quark_Pythia->GetYaxis()->SetTitle("d#sigma/d#it{y}_{c#bar{c}} (#mub)");
    Double_t pythia_cs_charm = (56420 * h_NCharm_event_Pythia->GetBinContent(3)) / (1e+05 * h_Y_Charm_quark_Pythia->Integral());
    hist1D_graphic_opt(h_Y_Charm_quark_Pythia, kTRUE, 10, 20, kRed, pythia_cs_charm);
    h_Y_Charm_quark_Pythia->SetTitle(Form("PYTHIA8 (%0.2e #mub)", h_Y_Charm_quark_Pythia->Integral()));
    h_Y_Charm_quark_Pythia->GetYaxis()->SetRangeUser(3.5e-01, 3.6e+05);

    TFile *fIn_Powheg_Charm_nocut = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/read_output/test_new_sim/powheg_charm_nocut_test_MC_output_Hist_294915.root", "READ");

    TH2F *h_PtY_Charm_quark_Powheg = (TH2F *)fIn_Powheg_Charm_nocut->Get("HF_quarks/h_PtY_Charm_quark");
    h_PtY_Charm_quark_Powheg->SetName("h_PtY_Charm_quark_Powheg");

    TH1F *h_Y_Charm_quark_Powheg = (TH1F *)h_PtY_Charm_quark_Powheg->ProjectionY();
    hist1D_graphic_opt(h_Y_Charm_quark_Powheg, kTRUE, 10, 21, kBlue, 5000. / h_Y_Charm_quark_Powheg->Integral());
    h_Y_Charm_quark_Powheg->SetTitle(Form("POWHEG (%0.0e #mub)", h_Y_Charm_quark_Powheg->Integral()));

    TH1F *h_Y_ratio_Powheg_Pythia = (TH1F *)h_Y_Charm_quark_Pythia->Clone("h_Y_ratio_Powheg_Pythia");
    h_Y_ratio_Powheg_Pythia->Divide(h_Y_Charm_quark_Powheg);
    h_Y_ratio_Powheg_Pythia->GetYaxis()->SetRangeUser(-0.2, 3.2);
    h_Y_ratio_Powheg_Pythia->GetYaxis()->SetTitle("PYTHIA/POWHEG");

    TCanvas *c_Y_Charm = two_histo_ratio(h_Y_Charm_quark_Pythia, h_Y_Charm_quark_Powheg, h_Y_ratio_Powheg_Pythia, "c_Y_Charm", "ALICE Simulation, single c#bar{c} pair cs", kTRUE, kTRUE);
    c_Y_Charm->GetPad(0)->SetGridx();
    c_Y_Charm->GetPad(0)->SetGridy();
    c_Y_Charm->SaveAs("images/charm_cs.png");

    //------------------------------------------//

    TH1F *h_NBeauty_event_Pythia = (TH1F *)fIn_Pythia_STANDALONE_NoDiffr->Get("HF_quarks/h_NBeauty_event");

    TH2F *h_PtY_Beauty_quark_Pythia = (TH2F *)fIn_Pythia_STANDALONE_NoDiffr->Get("HF_quarks/h_PtY_Beauty_quark");
    h_PtY_Beauty_quark_Pythia->SetName("h_PtY_Beauty_quark_Pythia");

    TH1F *h_Y_Beauty_quark_Pythia = (TH1F *)h_PtY_Beauty_quark_Pythia->ProjectionY();
    h_Y_Beauty_quark_Pythia->GetXaxis()->SetTitle("#it{y}");
    h_Y_Beauty_quark_Pythia->GetYaxis()->SetTitle("d#sigma/d#it{y}_{b#bar{b}} (#mub)");
    Double_t pythia_cs_Beauty = (56420 * h_NBeauty_event_Pythia->GetBinContent(3)) / (1e+05 * h_Y_Beauty_quark_Pythia->Integral());
    hist1D_graphic_opt(h_Y_Beauty_quark_Pythia, kTRUE, 10, 20, kRed, pythia_cs_Beauty);
    h_Y_Beauty_quark_Pythia->SetTitle(Form("PYTHIA8 (%0.2e #mub)", h_Y_Beauty_quark_Pythia->Integral()));
    h_Y_Beauty_quark_Pythia->GetYaxis()->SetRangeUser(3.5e-01, 6e+03);

    TFile *fIn_Powheg_Beauty_nocut = new TFile("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/read_output/test_new_sim/powheg_beauty_nocut_test_MC_output_Hist_294915.root", "READ");

    TH2F *h_PtY_Beauty_quark_Powheg = (TH2F *)fIn_Powheg_Beauty_nocut->Get("HF_quarks/h_PtY_Beauty_quark");
    h_PtY_Beauty_quark_Powheg->SetName("h_PtY_Beauty_quark_Powheg");

    TH1F *h_Y_Beauty_quark_Powheg = (TH1F *)h_PtY_Beauty_quark_Powheg->ProjectionY();
    hist1D_graphic_opt(h_Y_Beauty_quark_Powheg, kTRUE, 10, 21, kBlue, 500 / h_Y_Beauty_quark_Powheg->Integral());
    h_Y_Beauty_quark_Powheg->SetTitle(Form("POWHEG (%0.0e #mub)", h_Y_Beauty_quark_Powheg->Integral()));

    TH1F *h_Y_ratio_Powheg_Pythia_Beauty = (TH1F *)h_Y_Beauty_quark_Pythia->Clone("h_Y_ratio_Powheg_Pythia_Beauty");
    h_Y_ratio_Powheg_Pythia_Beauty->Divide(h_Y_Beauty_quark_Powheg);
    h_Y_ratio_Powheg_Pythia_Beauty->GetYaxis()->SetRangeUser(-0.2, 3.2);
    h_Y_ratio_Powheg_Pythia_Beauty->GetYaxis()->SetTitle("PYTHIA/POWHEG");

    TCanvas *c_Y_Beauty = two_histo_ratio(h_Y_Beauty_quark_Pythia, h_Y_Beauty_quark_Powheg, h_Y_ratio_Powheg_Pythia_Beauty, "c_Y_Beauty", "ALICE Simulation, single b#bar{b} pair cs", kTRUE, kTRUE);
    c_Y_Beauty->GetPad(0)->SetGridx();
    c_Y_Beauty->GetPad(0)->SetGridy();
    c_Y_Beauty->SaveAs("images/beauty_cs.png");
}

void mc_comparinson_dimuon()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TH1F *h_Nevents_PYTHIA = (TH1F *)fIn_PYTHIA_HF->Get("h_Nevents");
    Int_t Nev_PYTHIA = (Int_t)h_Nevents_PYTHIA->GetBinContent(2);
    Double_t Norm_PYTHIA = (1. / ((Double_t)Nev_PYTHIA * 216)); // in mb

    TString name_Dimuon_origin[] = {"Charm", "Beauty"};

    TFile *fIn_POWHEG[2];
    TH1F *h_Nevents_POWHEG[2];
    Int_t Nev_POWHEG[2];
    Double_t Norm_POWHEG[2];

    fIn_POWHEG[0] = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i1/Version3_AliAOD/save_output/LHC23i1_MC_output_Hist_merged.root", "READ");
    h_Nevents_POWHEG[0] = (TH1F *)fIn_POWHEG[0]->Get("h_Nevents");
    Nev_POWHEG[0] = (Int_t)h_Nevents_POWHEG[0]->GetBinContent(2);
    Norm_POWHEG[0] = (1. / (Nev_POWHEG[0] * 32.123)); // in mb

    fIn_POWHEG[1] = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version3_AliAOD/save_output/LHC23i2_MC_output_Hist_merged.root", "READ");
    h_Nevents_POWHEG[1] = (TH1F *)fIn_POWHEG[1]->Get("h_Nevents");
    Nev_POWHEG[1] = (Int_t)h_Nevents_POWHEG[1]->GetBinContent(2);
    Norm_POWHEG[1] = (1. / (Nev_POWHEG[1] * 15)); // in mb

    cout << "Nev_PYTHIA    " << Nev_PYTHIA << endl;
    cout << "Nev_POWHEG[0] " << Nev_POWHEG[0] << endl;
    cout << Nev_POWHEG[0] * 32.123 << Nev_POWHEG[0] * 32.123 << endl;
    cout << "Nev_POWHEG[1] " << Nev_POWHEG[1] << endl;

    cout << "Norm_PYTHIA    " << Norm_PYTHIA << endl;
    cout << "Norm_POWHEG[0] " << Norm_POWHEG[0] << endl;
    cout << "Norm_POWHEG[1] " << Norm_POWHEG[1] << endl;

    TString name_MassCuts[] = {"M0", "M4"};

    TH2F *h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[2][2];

    h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0] = (TH2F *)fIn_PYTHIA_HF->Get(Form("DiMuon/%s/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_%s_Rec_DQ_cut_match_LT_ULS_from%s", name_MassCuts[0].Data(), name_MassCuts[0].Data(), name_Dimuon_origin[0].Data()));
    // h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0]->Rebin2D(10, 10);

    TH1F *h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[2][2];
    TH1F *h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[2][2];

    h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0] = (TH1F *)h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0]->ProjectionX();
    h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0]->SetTitle(Form("PYTHIA8 %s", name_Dimuon_origin[0].Data()));
    h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0] = (TH1F *)h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0]->ProjectionY();
    h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0]->SetTitle(Form("PYTHIA8 %s", name_Dimuon_origin[0].Data()));

    Color_t colors_pythia[] = {kGreen, kMagenta};
    Color_t colors_powheg[] = {kRed, kBlack};

    hist1D_graphic_opt(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0], kTRUE, 10, 20, colors_pythia[0], Norm_PYTHIA);
    cout << "N dimu from charm PYTHIA8:   " << h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0]->GetEntries() << endl;
    cout << "Integral from charm PYTHIA8: " << h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0]->Integral() << endl;

    TH2F *h_PtM_DiMuon_Rec_PowhegOnly[2][2];
    TH1F *h_Pt_DiMuon_Rec_PowhegOnly[2][2];
    TH1F *h_M_DiMuon_Rec_PowhegOnly[2][2];

    h_PtM_DiMuon_Rec_PowhegOnly[0][0] = (TH2F *)fIn_POWHEG[0]->Get(Form("DiMuon_Rec/PowhegOnly/h_PtM_DiMuon_Rec_%s_PowhegOnly", name_Dimuon_origin[0].Data()));
    // h_PtM_DiMuon_Rec_PowhegOnly[0][0]->Rebin2D(10, 10);
    h_PtM_DiMuon_Rec_PowhegOnly[0][0]->GetXaxis()->SetRangeUser(0, 30);
    h_Pt_DiMuon_Rec_PowhegOnly[0][0] = (TH1F *)h_PtM_DiMuon_Rec_PowhegOnly[0][0]->ProjectionX();
    h_Pt_DiMuon_Rec_PowhegOnly[0][0]->SetTitle(Form("Powheg %s", name_Dimuon_origin[0].Data()));
    h_M_DiMuon_Rec_PowhegOnly[0][0] = (TH1F *)h_PtM_DiMuon_Rec_PowhegOnly[0][0]->ProjectionY();
    h_M_DiMuon_Rec_PowhegOnly[0][0]->SetTitle(Form("Powheg %s", name_Dimuon_origin[0].Data()));

    hist1D_graphic_opt(h_Pt_DiMuon_Rec_PowhegOnly[0][0], kTRUE, 10, 20, colors_powheg[0], Norm_POWHEG[0]);
    cout << "N dimu from charm POWHEG:    " << h_Pt_DiMuon_Rec_PowhegOnly[0][0]->GetEntries() << endl;
    cout << "Integral from charm POWHEG:  " << h_Pt_DiMuon_Rec_PowhegOnly[0][0]->Integral() << endl;

    TH1F *h_Pt_ratio_Powheg_Pythia[2][2];
    h_Pt_ratio_Powheg_Pythia[0][0] = (TH1F *)h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0]->Clone(Form("h_Pt_ratio_%s_%s", name_MassCuts[0].Data(), name_Dimuon_origin[0].Data()));
    h_Pt_ratio_Powheg_Pythia[0][0]->Divide(h_Pt_DiMuon_Rec_PowhegOnly[0][0]);
    // h_Pt_ratio_Powheg_Pythia[0][0]->GetYaxis()->SetRangeUser(-2,5);

    TCanvas *c_M_Charm[2];

    c_M_Charm[0] = two_histo_ratio(h_Pt_DiMuon_Rec_PowhegOnly[0][0], h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0], h_Pt_ratio_Powheg_Pythia[0][0], Form("c_M_Charm_%s", name_MassCuts[0].Data()), "#splitline{ALICE Simulation, Rec #mu^{#plus}#mu^{#minus}}{#it{p}_{T} < 30 GeV/#it{c} }", kTRUE, kTRUE);

    // new TCanvas();
    // gPad->SetLogy();
    // h_Pt_DiMuon_Rec_PowhegOnly[0][0]->Draw("PE");
    // h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][0]->Draw("PESAME");
}

void mc_comparinson_shape()
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TH1F *h_Nevents_PYTHIA = (TH1F *)fIn_PYTHIA_HF->Get("h_Nevents");
    Int_t Nev_PYTHIA = (Int_t)h_Nevents_PYTHIA->GetBinContent(2);
    Double_t CS_PYTHIA = (1. / ((Double_t)Nev_PYTHIA * 216)); // in mb

    TH2F *h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[2][2];
    TH1F *h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[2][2];
    TH1F *h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[2][2];

    TH2F *h_PtM_DiMuon_Rec_PowhegOnly[2][2];
    TH1F *h_Pt_DiMuon_Rec_PowhegOnly[2][2];
    TH1F *h_M_DiMuon_Rec_PowhegOnly[2][2];

    TString name_Dimuon_origin[] = {"Charm", "Beauty"};

    TFile *fIn_POWHEG[2];
    TH1F *h_Nevents_POWHEG[2];
    Int_t Nev_POWHEG[2];
    Double_t CS_POWHEG[2];

    fIn_POWHEG[0] = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i1/Version3_AliAOD/save_ouput/LHC23i1_MC_output_Hist_merged.root", "READ");
    h_Nevents_POWHEG[0] = (TH1F *)fIn_POWHEG[0]->Get("h_Nevents");
    Nev_POWHEG[0] = (Int_t)h_Nevents_POWHEG[0]->GetBinContent(2);
    CS_POWHEG[0] = (1. / (Nev_POWHEG[0] * 30.5)); // in mb

    fIn_POWHEG[1] = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC23i2/Version3_AliAOD/save_ouput/LHC23i2_MC_output_Hist_merged.root", "READ");
    h_Nevents_POWHEG[1] = (TH1F *)fIn_POWHEG[1]->Get("h_Nevents");
    Nev_POWHEG[1] = (Int_t)h_Nevents_POWHEG[1]->GetBinContent(2);
    CS_POWHEG[1] = (1. / (Nev_POWHEG[1] * 15)); // in mb
    cout << "Nev_POWHEG[0] " << Nev_POWHEG[0] << endl;
    cout << "Nev_POWHEG[1] " << Nev_POWHEG[1] << endl;

    cout << "CS_PYTHIA " << CS_PYTHIA << endl;
    cout << "CS_POWHEG[0] " << CS_POWHEG[0] << endl;
    cout << "CS_POWHEG[1] " << CS_POWHEG[1] << endl;

    Color_t colors_pythia[] = {kGreen, kMagenta};
    Color_t colors_powheg[] = {kRed, kBlack};

    TString name_MassCuts[] = {"M0", "M4"};

    TH1F *h_Pt_ratio_Powheg_Pythia[2][2];
    TH1F *h_M_ratio_Powheg_Pythia[2][2];

    for (Int_t i_Dimuon_origin = 0; i_Dimuon_origin < 2; i_Dimuon_origin++)
    {
        for (Int_t i_MassCuts = 0; i_MassCuts < 2; i_MassCuts++)
        {
            h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts] = (TH2F *)fIn_PYTHIA_HF->Get(Form("DiMuon/%s/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_%s_Rec_DQ_cut_match_LT_ULS_from%s", name_MassCuts[i_MassCuts].Data(), name_MassCuts[i_MassCuts].Data(), name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts]->Rebin2D(10, 10);
            h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts] = (TH1F *)h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts]->ProjectionX();
            h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts]->SetTitle(Form("PYTHIA8 %s", name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts] = (TH1F *)h_PtMDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts]->ProjectionY();
            h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts]->SetTitle(Form("PYTHIA8 %s", name_Dimuon_origin[i_Dimuon_origin].Data()));

            hist1D_graphic_opt(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts], kTRUE, 1, 20, colors_pythia[i_Dimuon_origin], CS_PYTHIA);
            hist1D_graphic_opt(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts], kTRUE, 1, 20, colors_pythia[i_Dimuon_origin], CS_PYTHIA);

            h_PtM_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts] = (TH2F *)fIn_POWHEG[i_Dimuon_origin]->Get(Form("DiMuon_Rec/PowhegOnly/h_PtM_DiMuon_Rec_%s_PowhegOnly", name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_PtM_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]->SetName(Form("h_PtM_DiMuon_Rec_PowhegOnly_%s_%s", name_MassCuts[i_MassCuts].Data(), name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_PtM_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]->Rebin2D(10, 10);
            h_PtM_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]->GetXaxis()->SetRangeUser(0, 30);
            if (i_MassCuts == 0)
                h_PtM_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]->GetYaxis()->SetRangeUser(0, 30.0);
            else
                h_PtM_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]->GetYaxis()->SetRangeUser(4, 30);

            h_Pt_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts] = (TH1F *)h_PtM_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]->ProjectionX();
            h_Pt_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]->SetTitle(Form("Powheg %s", name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_M_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts] = (TH1F *)h_PtM_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]->ProjectionY();
            h_M_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]->SetTitle(Form("Powheg %s", name_Dimuon_origin[i_Dimuon_origin].Data()));

            hist1D_graphic_opt(h_Pt_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts], kTRUE, 1, 20, colors_powheg[i_Dimuon_origin], CS_POWHEG[i_Dimuon_origin]);
            hist1D_graphic_opt(h_M_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts], kTRUE, 1, 20, colors_powheg[i_Dimuon_origin], CS_POWHEG[i_Dimuon_origin]);

            h_Pt_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts] = (TH1F *)h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts]->Clone(Form("h_Pt_ratio_%s_%s", name_MassCuts[i_MassCuts].Data(), name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_Pt_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts]->Divide(h_Pt_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]);
            h_Pt_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts]->SetTitle(Form("%s PYTHIA/Powheg", name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_Pt_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts]->GetYaxis()->SetTitle(Form("%s PYTHIA/Powheg", name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_Pt_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts]->GetYaxis()->SetRangeUser(-0.5, 3.5);

            h_M_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts] = (TH1F *)h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[i_Dimuon_origin][i_MassCuts]->Clone(Form("h_M_ratio_%s_%s", name_MassCuts[i_MassCuts].Data(), name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_M_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts]->Divide(h_M_DiMuon_Rec_PowhegOnly[i_Dimuon_origin][i_MassCuts]);
            h_M_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts]->SetTitle(Form("%s PYTHIA/Powheg", name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_M_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts]->GetYaxis()->SetTitle(Form("%s PYTHIA/Powheg", name_Dimuon_origin[i_Dimuon_origin].Data()));
            h_M_ratio_Powheg_Pythia[i_Dimuon_origin][i_MassCuts]->GetYaxis()->SetRangeUser(-0.5, 3.5);
        }
    }

    TCanvas *c_Pt_Beauty[2];
    TCanvas *c_M_Beauty[2];
    TCanvas *c_Pt_Charm[2];
    TCanvas *c_M_Charm[2];
    cout << h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][0]->Integral() << endl;
    cout << h_Pt_DiMuon_Rec_PowhegOnly[1][0]->Integral() << endl;

    // h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][1]->GetYaxis()->SetRangeUser(2.5e-6, 225);
    // c_Pt_Beauty[1] = two_histo_ratio(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][1], h_Pt_DiMuon_Rec_PowhegOnly[1][1], h_Pt_ratio_Powheg_Pythia[1][1], "c_Pt_Beauty", " #splitline{ALICE Simulation}{Rec #mu^{#plus}#mu^{#minus}, #it{m}_{#mu#mu} > 4 GeV/#it{c}^{2} }", kTRUE, kTRUE);

    // h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][1]->GetYaxis()->SetRangeUser(2.5e-6, 225);
    // c_M_Beauty[1] = two_histo_ratio(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][1], h_M_DiMuon_Rec_PowhegOnly[1][1], h_M_ratio_Powheg_Pythia[1][1], "c_M_Beauty", " #splitline{ALICE Simulation}{Rec #mu^{#plus}#mu^{#minus}, #it{p}_{T} < 30 GeV/#it{c} }", kTRUE, kTRUE);

    // ///

    // h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][1]->GetYaxis()->SetRangeUser(2.5e-6, 225);
    // c_Pt_Charm[1] = two_histo_ratio(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][1], h_Pt_DiMuon_Rec_PowhegOnly[0][1], h_Pt_ratio_Powheg_Pythia[0][1], "c_Pt_Charm", " #splitline{ALICE Simulation}{Rec #mu^{#plus}#mu^{#minus}, #it{m}_{#mu#mu} > 4 GeV/#it{c}^{2} }", kTRUE, kTRUE);

    // h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][1]->GetYaxis()->SetRangeUser(2.5e-6, 225);
    // c_M_Charm[1] = two_histo_ratio(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][1], h_M_DiMuon_Rec_PowhegOnly[0][1], h_M_ratio_Powheg_Pythia[0][1], "c_M_Charm", " #splitline{ALICE Simulation}{Rec #mu^{#plus}#mu^{#minus}, #it{p}_{T} < 30 GeV/#it{c} }", kTRUE, kTRUE);

    TH1F *h_Pt_ratio_Charm_Beauty_Powheg[2];
    TH1F *h_M_ratio_Charm_Beauty_Powheg[2];
    TH1F *h_Pt_ratio_Charm_Beauty_Pythia[2];
    TH1F *h_M_ratio_Charm_Beauty_Pythia[2];

    TCanvas *c_Pt_ratio_Charm_Beauty_Powheg[2];
    TCanvas *c_M_ratio_Charm_Beauty_Powheg[2];
    TCanvas *c_Pt_ratio_Charm_Beauty_Pythia[2];
    TCanvas *c_M_ratio_Charm_Beauty_Pythia[2];

    TString info_m_cut[] = {"", "#it{m}_{#mu#mu} > 4 Gev/#it{c}^{2}"};

    Int_t counter_for_slides = 1;

    for (Int_t i_Masscuts = 0; i_Masscuts < 2; i_Masscuts++)
    {
        //--------Powheg Pythia ratio canvas -----------------//
        h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][i_Masscuts]->GetYaxis()->SetRangeUser(2.5e-6, 225);
        c_Pt_Beauty[i_Masscuts] = two_histo_ratio(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][i_Masscuts], h_Pt_DiMuon_Rec_PowhegOnly[1][i_Masscuts], h_Pt_ratio_Powheg_Pythia[1][i_Masscuts], Form("c_Pt_Beauty_%s", name_MassCuts[i_Masscuts].Data()), Form("#splitline{ALICE Simulation, Rec #mu^{#plus}#mu^{#minus}}{%s}", info_m_cut[i_Masscuts].Data()), kTRUE, kTRUE);
        // c_Pt_Beauty[i_Masscuts]->SaveAs(Form("images/mc_comparison/%d_%s.png",counter_for_slides,c_Pt_Beauty[i_Masscuts]->GetName()));
        counter_for_slides++;

        h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][i_Masscuts]->GetYaxis()->SetRangeUser(2.5e-6, 225);
        c_M_Beauty[i_Masscuts] = two_histo_ratio(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][i_Masscuts], h_M_DiMuon_Rec_PowhegOnly[1][i_Masscuts], h_M_ratio_Powheg_Pythia[1][i_Masscuts], Form("c_M_Beauty_%s", name_MassCuts[i_Masscuts].Data()), "#splitline{ALICE Simulation, Rec #mu^{#plus}#mu^{#minus}}{#it{p}_{T} < 30 GeV/#it{c} }", kTRUE, kTRUE);
        // c_M_Beauty[i_Masscuts]->SaveAs(Form("images/mc_comparison/%d_%s.png",counter_for_slides,c_M_Beauty[i_Masscuts]->GetName()));
        counter_for_slides++;

        h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][i_Masscuts]->GetYaxis()->SetRangeUser(2.5e-6, 225);
        c_Pt_Charm[i_Masscuts] = two_histo_ratio(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][i_Masscuts], h_Pt_DiMuon_Rec_PowhegOnly[0][i_Masscuts], h_Pt_ratio_Powheg_Pythia[0][i_Masscuts], Form("c_Pt_Charm_%s", name_MassCuts[i_Masscuts].Data()), Form("#splitline{ALICE Simulation, Rec #mu^{#plus}#mu^{#minus}}{%s}", info_m_cut[i_Masscuts].Data()), kTRUE, kTRUE);
        // c_Pt_Charm[i_Masscuts]->SaveAs(Form("images/mc_comparison/%d_%s.png",counter_for_slides,c_Pt_Charm[i_Masscuts]->GetName()));
        counter_for_slides++;

        h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][i_Masscuts]->GetYaxis()->SetRangeUser(2.5e-6, 225);
        c_M_Charm[i_Masscuts] = two_histo_ratio(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][i_Masscuts], h_M_DiMuon_Rec_PowhegOnly[0][i_Masscuts], h_M_ratio_Powheg_Pythia[0][i_Masscuts], Form("c_M_Charm_%s", name_MassCuts[i_Masscuts].Data()), "#splitline{ALICE Simulation, Rec #mu^{#plus}#mu^{#minus}}{#it{p}_{T} < 30 GeV/#it{c} }", kTRUE, kTRUE);
        // c_M_Charm[i_Masscuts]->SaveAs(Form("images/mc_comparison/%d_%s.png",counter_for_slides,c_M_Charm[i_Masscuts]->GetName()));
        counter_for_slides++;

        //--------Powheg Beauty ratio canvas -----------------//
        h_Pt_ratio_Charm_Beauty_Powheg[i_Masscuts] = (TH1F *)h_Pt_DiMuon_Rec_PowhegOnly[0][i_Masscuts]->Clone(Form("h_Pt_ratio_Charm_Beauty_Powheg_%s", name_MassCuts[i_Masscuts].Data()));
        h_Pt_ratio_Charm_Beauty_Powheg[i_Masscuts]->Divide(h_Pt_DiMuon_Rec_PowhegOnly[1][i_Masscuts]);
        h_Pt_ratio_Charm_Beauty_Powheg[i_Masscuts]->GetYaxis()->SetRangeUser(-0.5, 3.5);

        h_Pt_DiMuon_Rec_PowhegOnly[0][i_Masscuts]->GetYaxis()->SetRangeUser(2.5e-6, 225);
        c_Pt_ratio_Charm_Beauty_Powheg[i_Masscuts] = two_histo_ratio(h_Pt_DiMuon_Rec_PowhegOnly[0][i_Masscuts], h_Pt_DiMuon_Rec_PowhegOnly[1][i_Masscuts], h_Pt_ratio_Charm_Beauty_Powheg[i_Masscuts], Form("c_Pt_ratio_Charm_Beauty_Powheg_%s", name_MassCuts[i_Masscuts].Data()), Form("#splitline{ALICE Simulation, Rec #mu^{#plus}#mu^{#minus}}{%s}", info_m_cut[i_Masscuts].Data()), kTRUE, kTRUE);
        // c_Pt_ratio_Charm_Beauty_Powheg[i_Masscuts]->SaveAs(Form("images/mc_comparison/%d_%s.png",counter_for_slides,c_Pt_ratio_Charm_Beauty_Powheg[i_Masscuts]->GetName()));
        counter_for_slides++;

        h_M_ratio_Charm_Beauty_Powheg[i_Masscuts] = (TH1F *)h_M_DiMuon_Rec_PowhegOnly[0][i_Masscuts]->Clone(Form("h_M_ratio_Charm_Beauty_Powheg_%s", name_MassCuts[i_Masscuts].Data()));
        h_M_ratio_Charm_Beauty_Powheg[i_Masscuts]->Divide(h_M_DiMuon_Rec_PowhegOnly[1][i_Masscuts]);
        h_M_ratio_Charm_Beauty_Powheg[i_Masscuts]->GetYaxis()->SetRangeUser(-0.5, 3.5);

        h_M_DiMuon_Rec_PowhegOnly[0][i_Masscuts]->GetYaxis()->SetRangeUser(2.5e-6, 225);
        c_M_ratio_Charm_Beauty_Powheg[i_Masscuts] = two_histo_ratio(h_M_DiMuon_Rec_PowhegOnly[0][i_Masscuts], h_M_DiMuon_Rec_PowhegOnly[1][i_Masscuts], h_M_ratio_Charm_Beauty_Powheg[i_Masscuts], Form("c_M_ratio_Charm_Beauty_Powheg_%s", name_MassCuts[i_Masscuts].Data()), "#splitline{ALICE Simulation, Rec #mu^{#plus}#mu^{#minus}}{#it{p}_{T} < 30 GeV/#it{c} }", kTRUE, kTRUE);
        // c_M_ratio_Charm_Beauty_Powheg[i_Masscuts]->SaveAs(Form("images/mc_comparison/%d_%s.png",counter_for_slides,c_M_ratio_Charm_Beauty_Powheg[i_Masscuts]->GetName()));
        counter_for_slides++;

        h_Pt_ratio_Charm_Beauty_Pythia[i_Masscuts] = (TH1F *)h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][i_Masscuts]->Clone(Form("h_Pt_ratio_Charm_Beauty_Pythia_%s", name_MassCuts[i_Masscuts].Data()));
        h_Pt_ratio_Charm_Beauty_Pythia[i_Masscuts]->Divide(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][i_Masscuts]);
        h_Pt_ratio_Charm_Beauty_Pythia[i_Masscuts]->GetYaxis()->SetRangeUser(-0.5, 3.5);

        c_Pt_ratio_Charm_Beauty_Pythia[i_Masscuts] = two_histo_ratio(h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][i_Masscuts], h_PtDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][i_Masscuts], h_Pt_ratio_Charm_Beauty_Pythia[i_Masscuts], Form("c_Pt_ratio_Charm_Beauty_Pythia_%s", name_MassCuts[i_Masscuts].Data()), Form("#splitline{ALICE Simulation, Rec #mu^{#plus}#mu^{#minus}}{%s}", info_m_cut[i_Masscuts].Data()), kTRUE, kTRUE);
        // c_Pt_ratio_Charm_Beauty_Pythia[i_Masscuts]->SaveAs(Form("images/mc_comparison/%d_%s.png",counter_for_slides,c_Pt_ratio_Charm_Beauty_Pythia[i_Masscuts]->GetName()));
        counter_for_slides++;

        h_M_ratio_Charm_Beauty_Pythia[i_Masscuts] = (TH1F *)h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][i_Masscuts]->Clone(Form("h_M_ratio_Charm_Beauty_Pythia_%s", name_MassCuts[i_Masscuts].Data()));
        h_M_ratio_Charm_Beauty_Pythia[i_Masscuts]->Divide(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][i_Masscuts]);
        h_M_ratio_Charm_Beauty_Pythia[i_Masscuts]->GetYaxis()->SetRangeUser(-0.5, 3.5);

        c_M_ratio_Charm_Beauty_Pythia[i_Masscuts] = two_histo_ratio(h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[0][i_Masscuts], h_MDiMu_Rec_DQ_cut_match_LT_ULS_PYTHIA[1][i_Masscuts], h_M_ratio_Charm_Beauty_Pythia[i_Masscuts], Form("c_M_ratio_Charm_Beauty_Pythia_%s", name_MassCuts[i_Masscuts].Data()), "#splitline{ALICE Simulation, Rec #mu^{#plus}#mu^{#minus}}{#it{p}_{T} < 30 GeV/#it{c} }", kTRUE, kTRUE);
        // c_M_ratio_Charm_Beauty_Pythia[i_Masscuts]->SaveAs(Form("images/mc_comparison/%d_%s.png",counter_for_slides,c_Pt_ratio_Charm_Beauty_Pythia[i_Masscuts]->GetName()));
        counter_for_slides++;
    }
}