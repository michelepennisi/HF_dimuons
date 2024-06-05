#include "/home/michele_pennisi/cernbox/HF_dimuons/common_include.h"

void corr()
{

    TFile *fIn = TFile::Open("test_cs_Generator_SoftQCD_Def_pythia_sim_12345_DefaultBR_1000000.root", "READ");

    TH1F *h_NCharm_event = (TH1F *)fIn->Get("h_NCharm_event");
    h_NCharm_event->Draw("text");
    new TCanvas();

    Double_t n_HF_total_PYTHIA = 0;

    for (int i_bin = 1; i_bin <= h_NCharm_event->GetNbinsX(); i_bin++)
        n_HF_total_PYTHIA = (n_HF_total_PYTHIA + (0.5) * (h_NCharm_event->GetBinContent(i_bin) * (i_bin - 1)));

    cout << "Total charm pairs: " << n_HF_total_PYTHIA << endl;
    cout << "N charm single pair: " << h_NCharm_event->GetBinContent(3) << endl;
    cout << "MPI factor: " << n_HF_total_PYTHIA / h_NCharm_event->GetBinContent(3) << endl;

    TH2F *NDiCharm_y = (TH2F *)fIn->Get("NDiCharm_y");
    // NDiCharm_y->GetYaxis()->SetRangeUser(-4.0, -2.5);
    TH1F *NDiCharm = (TH1F *)NDiCharm_y->ProjectionX();
    NDiCharm->Draw("text");

    Double_t n_Diquark_pair_PYTHIA = 0;

    for (int i_bin = 1; i_bin <= NDiCharm->GetNbinsX(); i_bin++){
        // cout<<"i bin: "<<i_bin<<"NDiCharm->GetBinContent(i_bin): "<<NDiCharm->GetBinContent(i_bin)<<endl;
        n_Diquark_pair_PYTHIA = (n_Diquark_pair_PYTHIA + (NDiCharm->GetBinContent(i_bin)));
        // cout<<"n_Diquark_pair_PYTHIA: "<<n_Diquark_pair_PYTHIA<<endl;
    }

    cout << "Total DiCharm at fwd pairs: " << n_Diquark_pair_PYTHIA << endl;
    cout << "N DiCharm at fwd single pair: " << NDiCharm->GetBinContent(2) << endl;
    cout << "MPI factor at fwd: " << n_Diquark_pair_PYTHIA / NDiCharm->GetBinContent(2) << endl;
}

void std_PYTHIA()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TFile *fIn = TFile::Open("test_cs_Generator_SoftQCD_Def_pythia_sim_12345_DefaultBR_1000000.root", "READ");

    TH2F *Pt_test = (TH2F *)fIn->Get("NDiCharm_y");

    Pt_test->GetXaxis()->SetRangeUser(1, 2);

    TH1F *proj_1_2 = (TH1F *)Pt_test->ProjectionY();
    proj_1_2->SetName("proj_1_2");
    proj_1_2->SetTitle("N_{c#bar{c}, ev} = 1");
    Pt_test->GetXaxis()->SetRangeUser(2, 3);

    TH1F *proj_2_10 = (TH1F *)Pt_test->ProjectionY();
    proj_2_10->SetName("proj_2_10");
    proj_2_10->SetTitle("N_{c#bar{c}, ev} = 2");

    TH1F *total = (TH1F *)proj_1_2->Clone("total_dimuons");
    total->Add(proj_2_10,2);
    total->SetTitle("total");

    hist1D_graphic_opt(total, kTRUE, 5, 20, kBlack, 1. / total->Integral());

    hist1D_graphic_opt(proj_1_2, kTRUE, 5, 20, kRed, 1. / proj_1_2->Integral());

    TH1F *total_proj_1_2_ratio = (TH1F *)total->Clone("total_proj_1_2_ratio");
    total_proj_1_2_ratio->Reset();
    total_proj_1_2_ratio->Divide(total, proj_1_2, 1, 1, "B");
    total_proj_1_2_ratio->GetYaxis()->SetRangeUser(0.875, 1.125);

    hist1D_graphic_opt(proj_2_10, kTRUE, 5, 20, kGreen, 1. / proj_2_10->Integral());

    TCanvas *dimuon_rec = two_histo_ratio(total, proj_1_2, total_proj_1_2_ratio, "dimuon_rec", "", kTRUE, kTRUE, 0., 40.0);
    TPad *pad1 = (TPad *)dimuon_rec->GetListOfPrimitives()->FindObject("pad1");
    pad1->cd();
    proj_2_10->Draw("same");
    dimuon_rec->GetListOfPrimitives()->Print();
    TLegend *legend = (TLegend *)dimuon_rec->FindObject("Legend");
    legend->AddEntry(proj_2_10);
    Legend_settings(legend, legend->GetX1(), legend->GetX2(), legend->GetY1() - 0.1, legend->GetY2(), "");
    legend->Draw();
    pad1->Modified();
    pad1->Update();
    return;
}

void MPI_test(TString HF = "Charm", TString mass_cut = "mass4")
{
    Double_t Mass_cut = 999;
    if (mass_cut.Contains("mass0"))
        Mass_cut = 0;
    else if (mass_cut.Contains("mass4"))
        Mass_cut = 4;
    TH2F *Pt_test = new TH2F(Form("Pt_test_%s_%s", HF.Data(), mass_cut.Data()), Form("Pt_test_%s_%s; p_{T} (GeV/#it{c}); N_{#mu#mu}", HF.Data(), mass_cut.Data()), 400, 0, 40, 10, 0, 10);
    TH1F *N_ev_test = new TH1F(Form("N_ev_test_%s_%s", HF.Data(), mass_cut.Data()), Form("N_ev_test_%s_%s; N_{#mu#mu}", HF.Data(), mass_cut.Data()), 10, 0, 10);
    TH2F *Mass_test = new TH2F(Form("Mass_test_%s_%s", HF.Data(), mass_cut.Data()), Form("Mass_test_%s_%s; #it{m}_{#mu#mu} (GeV/#it{c}^{2}); N_{#mu#mu}", HF.Data(), mass_cut.Data()), (40 - Mass_cut) * 10, Mass_cut, 40, 10, 0, 10);
    // TFile *fIn = TFile::Open("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/read_output/test/LHC23i2_MC_output_Tree_294009.root");
    TFile *fIn = TFile::Open("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/LHC22b3/Version_5_AliAOD_skimmed_fwd_fullstat_mod/LHC22b3_MC_output_Tree_merged.root");

    TTree *DiMuon_Rec_Beauty = (TTree *)fIn->Get(Form("DiMuon_Rec_%s", HF.Data()));

    Double_t pt;
    Double_t m;
    Int_t event;
    DiMuon_Rec_Beauty->SetBranchAddress("pt", &pt);
    DiMuon_Rec_Beauty->SetBranchAddress("m", &m);
    DiMuon_Rec_Beauty->SetBranchAddress("event", &event);

    Int_t prev = 999;

    vector<Double_t> Pt;
    vector<Double_t> Mass;

    vector<Double_t> to_skip_event;

    for (Int_t i = 0; i < DiMuon_Rec_Beauty->GetEntries(); i++)
    {
        Int_t N_dimu_x_ev = 1;

        if (count(to_skip_event.begin(), to_skip_event.end(), i))
            continue;
        DiMuon_Rec_Beauty->GetEntry(i);
        if (m < Mass_cut)
            continue;
        Int_t studied_event = event;
        cout << "==============" << endl;
        cout << "studied_event " << studied_event << endl;
        cout << "i " << i << endl;
        cout << "event " << event << endl;
        cout << "pt correct " << pt << endl;
        cout << "m correct " << m << endl;

        Pt.push_back(pt);
        Mass.push_back(m);
        Int_t end = 999;
        if (i < DiMuon_Rec_Beauty->GetEntries() - 10)
            end = i + 9;
        else
            end = DiMuon_Rec_Beauty->GetEntries();
        for (Int_t j = i + 1; j < end; j++)
        {
            DiMuon_Rec_Beauty->GetEntry(j);
            if (m < Mass_cut)
                continue;
            Int_t next_event = event;

            if (studied_event == next_event)
            {
                cout << "second event " << studied_event << endl;
                cout << "j " << j << endl;
                cout << "event " << event << endl;
                cout << "pt correct " << pt << endl;
                cout << "m correct " << m << endl;
                N_dimu_x_ev += 1;
                Pt.push_back(pt);
                Mass.push_back(m);
                to_skip_event.push_back(j);
            }
        }

        cout << "N_dimu_x_ev: " << N_dimu_x_ev << endl;
        cout << "Pt.size(): " << Pt.size() << endl;
        N_ev_test->Fill(Pt.size());
        for (size_t i = 0; i < Pt.size(); i++)
        {
            cout << "Pt element: " << Pt[i] << endl;
            Pt_test->Fill(Pt[i], (Int_t)Pt.size());
            Mass_test->Fill(Mass[i], (Int_t)Mass.size());
        }

        Pt.clear();
        Mass.clear();
    }
    TFile *fout = TFile::Open("MPI_test.root", "UPDATE");
    fout->cd();
    Pt_test->Write(0, 2, 0);
    Mass_test->Write(0, 2, 0);
    // N_ev_test->Draw();
    // new TCanvas();
    // // Pt_test->ProjectionX()->SetLineColor(kRed);
    // Pt_test->Draw("same");
    // new TCanvas();
    // DiMuon_Rec_Beauty->Draw("m>>(400,0,40)");
    // Mass_test->ProjectionX()->Draw("same");
    // gPad->BuildLegend();
}

void test(TString mass_cut = "mass0")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TFile *fIn = TFile::Open("MPI_test.root", "READ");

    TH2F *Pt_test = (TH2F *)fIn->Get(Form("Pt_test_Charm_%s", mass_cut.Data()));

    Pt_test->GetYaxis()->SetRangeUser(1, 2);

    TH1F *proj_1_2 = (TH1F *)Pt_test->ProjectionX();
    proj_1_2->SetName("proj_1_2");
    proj_1_2->SetTitle("N_{#mu#mu, ev} = 1");
    Pt_test->GetYaxis()->SetRangeUser(2, 10);

    TH1F *proj_2_10 = (TH1F *)Pt_test->ProjectionX();
    proj_2_10->SetName("proj_2_10");
    proj_2_10->SetTitle("N_{#mu#mu, ev} > 1");

    TH1F *total = (TH1F *)proj_1_2->Clone("total_dimuons");
    total->Add(proj_2_10);
    total->SetTitle("total");

    hist1D_graphic_opt(total, kTRUE, 5, 20, kBlack, 1. / total->Integral());

    hist1D_graphic_opt(proj_1_2, kTRUE, 5, 20, kRed, 1. / proj_1_2->Integral());

    TH1F *total_proj_1_2_ratio = (TH1F *)total->Clone("total_proj_1_2_ratio");
    total_proj_1_2_ratio->Reset();
    total_proj_1_2_ratio->Divide(total, proj_1_2, 1, 1, "B");
    total_proj_1_2_ratio->GetYaxis()->SetRangeUser(0.875, 1.125);

    hist1D_graphic_opt(proj_2_10, kTRUE, 5, 20, kGreen, 1. / proj_2_10->Integral());

    TCanvas *dimuon_rec = two_histo_ratio(total, proj_1_2, total_proj_1_2_ratio, "dimuon_rec", "", kTRUE, kTRUE, 0., 40.0);
    TPad *pad1 = (TPad *)dimuon_rec->GetListOfPrimitives()->FindObject("pad1");
    pad1->cd();
    proj_2_10->Draw("same");
    dimuon_rec->GetListOfPrimitives()->Print();
    TLegend *legend = (TLegend *)dimuon_rec->FindObject("Legend");
    legend->AddEntry(proj_2_10);
    if (mass_cut.Contains("mass0"))
        Legend_settings(legend, legend->GetX1() - 0.1, legend->GetX2() - 0.1, legend->GetY1() + 0.05, legend->GetY2() + 0.1, "Reconstructed #mu#mu #leftarrow c,c, #it{m} > 0 GeV/#it{c}^{2}");
    else if (mass_cut.Contains("mass4"))
        Legend_settings(legend, legend->GetX1(), legend->GetX2(), legend->GetY1() - 0.1, legend->GetY2(), "Reconstructed #mu#mu #leftarrow c,c, #it{m} > 4 GeV/#it{c}^{2}");
    legend->Draw();
    pad1->Modified();
    pad1->Update();
    return;
}