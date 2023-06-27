#include "/home/michele_pennisi/cernbox/common_include.h"

void dimu_origin()
{
    // gStyle->SetOptStat(0);
    // gStyle->SetOptTitle(0);

    // TFile *fIn = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Version1/HF/save_mc_output/HF_MC_output_Hist_merged.root", "READ");

    TFile *fIn = new TFile("test/HF_MC_output_Hist_294009.root", "READ");

    TH3F *h_Pdg1Pdg2Pt_DiMuon_Gen = (TH3F *)fIn->Get("DiMuon_Gen/h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_Charm_corrected");
    TH3F *h_Pdg1Pdg2M_DiMuon_Gen = (TH3F *)fIn->Get("DiMuon_Gen/h_Pdg1Pdg2M_DiMuon_Gen_DQcut_Charm_corrected");

    h_Pdg1Pdg2Pt_DiMuon_Gen->GetXaxis()->SetRangeUser(400, 500);
    h_Pdg1Pdg2Pt_DiMuon_Gen->GetYaxis()->SetRangeUser(400, 500);
    TH2F *h_Pdg1Pdg2_DiMuon_Gen_mesons = (TH2F *)h_Pdg1Pdg2Pt_DiMuon_Gen->Project3D("xy");
    h_Pdg1Pdg2_DiMuon_Gen_mesons->SetName("charm_mesons_pdg1_pdg2");

    h_Pdg1Pdg2Pt_DiMuon_Gen->GetXaxis()->SetRangeUser(4100, 4300);
    h_Pdg1Pdg2Pt_DiMuon_Gen->GetYaxis()->SetRangeUser(4100, 4300);
    TH2F *h_Pdg1Pdg2_DiMuon_Gen_barions = (TH2F *)h_Pdg1Pdg2Pt_DiMuon_Gen->Project3D("xy");
    h_Pdg1Pdg2_DiMuon_Gen_barions->SetName("charm_barions_pdg1_pdg2");

    h_Pdg1Pdg2Pt_DiMuon_Gen->GetXaxis()->SetRangeUser(400, 500);
    h_Pdg1Pdg2Pt_DiMuon_Gen->GetYaxis()->SetRangeUser(4100, 4300);
    TH2F *h_Pdg1Pdg2_DiMuon_Gen_mesons_barions = (TH2F *)h_Pdg1Pdg2Pt_DiMuon_Gen->Project3D("xy");
    h_Pdg1Pdg2_DiMuon_Gen_mesons_barions->SetName("charm_mesons_barions_pdg1_pdg2");

    h_Pdg1Pdg2Pt_DiMuon_Gen->GetXaxis()->SetRangeUser(4100, 4300);
    h_Pdg1Pdg2Pt_DiMuon_Gen->GetYaxis()->SetRangeUser(400, 500);
    TH2F *h_Pdg1Pdg2_DiMuon_Gen_barions_mesons = (TH2F *)h_Pdg1Pdg2Pt_DiMuon_Gen->Project3D("xy");
    h_Pdg1Pdg2_DiMuon_Gen_barions_mesons->SetName("charm_barions_mesons_pdg1_pdg2");

    TCanvas *c_2d_mesons = canvas_2d_print(h_Pdg1Pdg2_DiMuon_Gen_mesons, "c_2d_mesons");

    TCanvas *c_2d_barion = canvas_2d_print(h_Pdg1Pdg2_DiMuon_Gen_barions, "c_2d_barion");

    TCanvas *c_2d_mesons_barions = canvas_2d_print(h_Pdg1Pdg2_DiMuon_Gen_mesons_barions, "c_2d_mesons_barions");

    TCanvas *c_2d_barions_mesons = canvas_2d_print(h_Pdg1Pdg2_DiMuon_Gen_barions_mesons, "c_2d_barions_mesons");

    const Int_t n_charm_meson_origin = 10;
    Int_t charm_meson_origin[n_charm_meson_origin] = {410, 420, 430, 440, 450, 4120, 4130, 4140, 4230, 4240};
    TString charm_meson_origin_name[n_charm_meson_origin] = {"Dplus", "Dzero", "Dstrange", "Jpsi", " ", "Lambda", "Xi_c0", " ", "Xi_c+"};
    TString charm_meson_origin_title[n_charm_meson_origin] = {"D^{#plus}", "D^{0}", "D^{#plus}_{s}", "J/#psi", " ", "#Lambda^{+}_{c}", "#Xi^{0}_{c}", " ", "#Xi^{+}_{c}"};
    Int_t original_pdg[n_charm_meson_origin] = {411, 421, 431, 443, 99999, 4122, 4132, 99999, 4232};
    TH1F *h_Pt_DiMuon_Gen[49];
    TH1F *h_M_DiMuon_Gen[49];
    Int_t index = 0;
    for (Int_t i_charm_meson_origin_mu1 = 0; i_charm_meson_origin_mu1 < n_charm_meson_origin - 1; i_charm_meson_origin_mu1++)
    {
        if ((charm_meson_origin[i_charm_meson_origin_mu1] == 450 && charm_meson_origin[i_charm_meson_origin_mu1 + 1] == 4120) || (charm_meson_origin[i_charm_meson_origin_mu1] == 4140 && charm_meson_origin[i_charm_meson_origin_mu1 + 1] == 4230))
            continue;
        h_Pdg1Pdg2Pt_DiMuon_Gen->GetXaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu1], charm_meson_origin[i_charm_meson_origin_mu1 + 1]);
        h_Pdg1Pdg2M_DiMuon_Gen->GetXaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu1], charm_meson_origin[i_charm_meson_origin_mu1 + 1]);

        for (Int_t i_charm_meson_origin_mu2 = 0; i_charm_meson_origin_mu2 < n_charm_meson_origin - 1; i_charm_meson_origin_mu2++)
        {
            if ((charm_meson_origin[i_charm_meson_origin_mu2] == 450 && charm_meson_origin[i_charm_meson_origin_mu2 + 1] == 4120) || (charm_meson_origin[i_charm_meson_origin_mu2] == 4140 && charm_meson_origin[i_charm_meson_origin_mu2 + 1] == 4230))
                continue;
            cout << "Bound1 (" << charm_meson_origin[i_charm_meson_origin_mu1] << "," << charm_meson_origin[i_charm_meson_origin_mu1 + 1] << ")  ------> " << charm_meson_origin_name[i_charm_meson_origin_mu1].Data() << " original pdg: " << original_pdg[i_charm_meson_origin_mu1] << endl;
            cout << "Bound2 (" << charm_meson_origin[i_charm_meson_origin_mu2] << "," << charm_meson_origin[i_charm_meson_origin_mu2 + 1] << ")  ------> " << charm_meson_origin_name[i_charm_meson_origin_mu2].Data() << " original pdg: " << original_pdg[i_charm_meson_origin_mu2] << endl;
            //---Pt---//
            h_Pdg1Pdg2Pt_DiMuon_Gen->GetYaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu2], charm_meson_origin[i_charm_meson_origin_mu2 + 1]);
            h_Pt_DiMuon_Gen[index] = (TH1F *)h_Pdg1Pdg2Pt_DiMuon_Gen->Project3D("z");
            h_Pt_DiMuon_Gen[index]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
            h_Pt_DiMuon_Gen[index]->SetName(Form("charm_%s_%s_pt_%d", charm_meson_origin_name[i_charm_meson_origin_mu1].Data(), charm_meson_origin_name[i_charm_meson_origin_mu2].Data(), index));
            h_Pt_DiMuon_Gen[index]->SetTitle(Form("#mu #leftarrow %s, #mu #leftarrow %s", charm_meson_origin_title[i_charm_meson_origin_mu1].Data(), charm_meson_origin_title[i_charm_meson_origin_mu2].Data()));
            //---M---//
            h_Pdg1Pdg2M_DiMuon_Gen->GetYaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu2], charm_meson_origin[i_charm_meson_origin_mu2 + 1]);
            h_M_DiMuon_Gen[index] = (TH1F *)h_Pdg1Pdg2M_DiMuon_Gen->Project3D("z");
            h_M_DiMuon_Gen[index]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
            h_M_DiMuon_Gen[index]->SetName(Form("charm_%s_%s_M_%d", charm_meson_origin_name[i_charm_meson_origin_mu1].Data(), charm_meson_origin_name[i_charm_meson_origin_mu2].Data(), index));
            h_M_DiMuon_Gen[index]->SetTitle(Form("#mu #leftarrow %s, #mu #leftarrow %s", charm_meson_origin_title[i_charm_meson_origin_mu1].Data(), charm_meson_origin_title[i_charm_meson_origin_mu2].Data()));
            cout << "index: " << index << endl;
            index++;
        }
    }
    //------Barions ---------------//

    TH1F *h_Pt_DiMuon_Gen_Barions = (TH1F *)h_Pt_DiMuon_Gen[32]->Clone("h_Pt_DiMuon_Gen_Barions");
    h_Pt_DiMuon_Gen_Barions->Add(h_Pt_DiMuon_Gen[40]);
    h_Pt_DiMuon_Gen_Barions->Add(h_Pt_DiMuon_Gen[48]);

    

    hist1D_graphic_opt(h_Pt_DiMuon_Gen[32], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[32]->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen[40], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[40]->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen[48], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[48]->GetEntries());

    TH1F *h_M_DiMuon_Gen_Barions = (TH1F *)h_M_DiMuon_Gen[32]->Clone("h_M_DiMuon_Gen_Barions");
    h_M_DiMuon_Gen_Barions->Add(h_M_DiMuon_Gen[40]);
    h_M_DiMuon_Gen_Barions->Add(h_M_DiMuon_Gen[48]);

    hist1D_graphic_opt(h_M_DiMuon_Gen[32], 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen[32]->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen[40], 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen[40]->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen[48], 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen[48]->GetEntries());

    TCanvas *c_Barions = canvas_noratio_divide2("c_Barions");
    c_Barions->cd(1);

    h_Pt_DiMuon_Gen[32]->Draw("PE PLC PMC");
    h_Pt_DiMuon_Gen[40]->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen[48]->Draw("PESAME PLC PMC");
    

    TLegend *legend_Barions = new TLegend(0.45, 0.45, 0.85, 0.7);
    legend_Barions->SetBorderSize(0);
    legend_Barions->SetTextSize(0.04);
    legend_Barions->SetFillStyle(0);
    legend_Barions->AddEntry(h_Pt_DiMuon_Gen[32], Form("%s,#Xi^{0}_{c} (%0.0f)", h_Pt_DiMuon_Gen[32]->GetTitle(), h_Pt_DiMuon_Gen[32]->GetEntries()));
    legend_Barions->AddEntry(h_Pt_DiMuon_Gen[40], Form("%s,#Xi^{0}_{c} (%0.0f)", h_Pt_DiMuon_Gen[40]->GetTitle(), h_Pt_DiMuon_Gen[40]->GetEntries()));
    legend_Barions->AddEntry(h_Pt_DiMuon_Gen[48], Form("%s,#Xi^{0}_{c} (%0.0f)", h_Pt_DiMuon_Gen[48]->GetTitle(), h_Pt_DiMuon_Gen[48]->GetEntries()));
    legend_Barions->Draw();

    c_Barions->cd(2);

    h_M_DiMuon_Gen[32]->Draw("PE PLC PMC");
    h_M_DiMuon_Gen[40]->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen[48]->Draw("PESAME PLC PMC");

    

    //---- Dmesons and Xi---------------//

    TH1F *h_Pt_DiMuon_Gen_Dplus_Xizero = (TH1F *)h_Pt_DiMuon_Gen[5]->Clone("h_Pt_DiMuon_Gen_Dplus_Xizero");
    h_Pt_DiMuon_Gen_Dplus_Xizero->Add(h_Pt_DiMuon_Gen[35]);

    TH1F *h_Pt_DiMuon_Gen_Dplus_Xiplus = (TH1F *)h_Pt_DiMuon_Gen[6]->Clone("h_Pt_DiMuon_Gen_Dplus_Xiplus");
    h_Pt_DiMuon_Gen_Dplus_Xiplus->Add(h_Pt_DiMuon_Gen[42]);

    TH1F *h_Pt_DiMuon_Gen_Dplus_Xi = (TH1F *)h_Pt_DiMuon_Gen_Dplus_Xiplus->Clone("h_Pt_DiMuon_Gen_Dplus_Xi");
    h_Pt_DiMuon_Gen_Dplus_Xi->Add(h_Pt_DiMuon_Gen_Dplus_Xizero);

    

    TH1F *h_Pt_DiMuon_Gen_Dzero_Xizero = (TH1F *)h_Pt_DiMuon_Gen[12]->Clone("h_Pt_DiMuon_Gen_Dzero_Xizero");
    h_Pt_DiMuon_Gen_Dzero_Xizero->Add(h_Pt_DiMuon_Gen[36]);

    TH1F *h_Pt_DiMuon_Gen_Dzero_Xiplus = (TH1F *)h_Pt_DiMuon_Gen[13]->Clone("h_Pt_DiMuon_Gen_Dzero_Xiplus");
    h_Pt_DiMuon_Gen_Dzero_Xiplus->Add(h_Pt_DiMuon_Gen[43]);

    TH1F *h_Pt_DiMuon_Gen_Dzero_Xi = (TH1F *)h_Pt_DiMuon_Gen_Dzero_Xiplus->Clone("h_Pt_DiMuon_Gen_Dzero_Xi");
    h_Pt_DiMuon_Gen_Dzero_Xi->Add(h_Pt_DiMuon_Gen_Dzero_Xizero);

    h_Pt_DiMuon_Gen_Dzero_Xi->Draw();
        

    TH1F *h_Pt_DiMuon_Gen_Dstrange_Xizero = (TH1F *)h_Pt_DiMuon_Gen[19]->Clone("h_Pt_DiMuon_Gen_Dstrange_Xizero");
    h_Pt_DiMuon_Gen_Dstrange_Xizero->Add(h_Pt_DiMuon_Gen[37]);

    TH1F *h_Pt_DiMuon_Gen_Dstrange_Xiplus = (TH1F *)h_Pt_DiMuon_Gen[20]->Clone("h_Pt_DiMuon_Gen_Dstrange_Xiplus");
    h_Pt_DiMuon_Gen_Dstrange_Xiplus->Add(h_Pt_DiMuon_Gen[44]);

    TH1F *h_Pt_DiMuon_Gen_Dstrange_Xi = (TH1F *)h_Pt_DiMuon_Gen_Dstrange_Xiplus->Clone("h_Pt_DiMuon_Gen_Dstrange_Xi");
    h_Pt_DiMuon_Gen_Dstrange_Xi->Add(h_Pt_DiMuon_Gen_Dstrange_Xizero);

    TH1F *h_Pt_DiMuon_Gen_Dmesons_Xi = (TH1F *)h_Pt_DiMuon_Gen_Dstrange_Xiplus->Clone("h_Pt_DiMuon_Gen_Dmesons_Xi");
    h_Pt_DiMuon_Gen_Dmesons_Xi->Add(h_Pt_DiMuon_Gen_Dzero_Xi);
    h_Pt_DiMuon_Gen_Dmesons_Xi->Add(h_Pt_DiMuon_Gen_Dplus_Xi);

    


    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dplus_Xi, 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Dplus_Xi->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dzero_Xi, 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Dzero_Xi->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dstrange_Xi, 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Dstrange_Xi->GetEntries());

    TH1F *h_M_DiMuon_Gen_Dplus_Xizero = (TH1F *)h_M_DiMuon_Gen[5]->Clone("h_M_DiMuon_Gen_Dplus_Xizero");
    h_M_DiMuon_Gen_Dplus_Xizero->Add(h_M_DiMuon_Gen[35]);

    TH1F *h_M_DiMuon_Gen_Dplus_Xiplus = (TH1F *)h_M_DiMuon_Gen[6]->Clone("h_M_DiMuon_Gen_Dplus_Xiplus");
    h_M_DiMuon_Gen_Dplus_Xiplus->Add(h_M_DiMuon_Gen[42]);

    TH1F *h_M_DiMuon_Gen_Dplus_Xi = (TH1F *)h_M_DiMuon_Gen_Dplus_Xiplus->Clone("h_M_DiMuon_Gen_Dplus_Xi");
    h_M_DiMuon_Gen_Dplus_Xi->Add(h_M_DiMuon_Gen_Dplus_Xizero);

    TH1F *h_M_DiMuon_Gen_Dzero_Xizero = (TH1F *)h_M_DiMuon_Gen[12]->Clone("h_M_DiMuon_Gen_Dzero_Xizero");
    h_M_DiMuon_Gen_Dzero_Xizero->Add(h_M_DiMuon_Gen[36]);

    TH1F *h_M_DiMuon_Gen_Dzero_Xiplus = (TH1F *)h_M_DiMuon_Gen[13]->Clone("h_M_DiMuon_Gen_Dzero_Xiplus");
    h_M_DiMuon_Gen_Dzero_Xiplus->Add(h_M_DiMuon_Gen[43]);

    TH1F *h_M_DiMuon_Gen_Dzero_Xi = (TH1F *)h_M_DiMuon_Gen_Dzero_Xiplus->Clone("h_M_DiMuon_Gen_Dzero_Xi");
    h_M_DiMuon_Gen_Dzero_Xi->Add(h_M_DiMuon_Gen_Dzero_Xizero);

    TH1F *h_M_DiMuon_Gen_Dstrange_Xizero = (TH1F *)h_M_DiMuon_Gen[19]->Clone("h_M_DiMuon_Gen_Dstrange_Xizero");
    h_M_DiMuon_Gen_Dstrange_Xizero->Add(h_M_DiMuon_Gen[37]);

    TH1F *h_M_DiMuon_Gen_Dstrange_Xiplus = (TH1F *)h_M_DiMuon_Gen[20]->Clone("h_M_DiMuon_Gen_Dstrange_Xiplus");
    h_M_DiMuon_Gen_Dstrange_Xiplus->Add(h_M_DiMuon_Gen[44]);

    TH1F *h_M_DiMuon_Gen_Dstrange_Xi = (TH1F *)h_M_DiMuon_Gen_Dstrange_Xiplus->Clone("h_M_DiMuon_Gen_Dstrange_Xi");
    h_M_DiMuon_Gen_Dstrange_Xi->Add(h_M_DiMuon_Gen_Dstrange_Xizero);

    TH1F *h_M_DiMuon_Gen_Dmesons_Xi = (TH1F *)h_M_DiMuon_Gen_Dstrange_Xiplus->Clone("h_M_DiMuon_Gen_Dmesons_Xi");
    h_M_DiMuon_Gen_Dmesons_Xi->Add(h_M_DiMuon_Gen_Dzero_Xi);
    h_M_DiMuon_Gen_Dmesons_Xi->Add(h_M_DiMuon_Gen_Dplus_Xi);

    hist1D_graphic_opt(h_M_DiMuon_Gen_Dplus_Xi, 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen_Dplus_Xi->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dzero_Xi, 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen_Dzero_Xi->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dstrange_Xi, 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen_Dstrange_Xi->GetEntries());

    TCanvas *c_Dmesons_Xi = canvas_noratio_divide2("c_Dmesons_Xi");
    c_Dmesons_Xi->cd(1);

    h_Pt_DiMuon_Gen_Dplus_Xi->Draw("PE PLC PMC");
    h_Pt_DiMuon_Gen_Dzero_Xi->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Dstrange_Xi->Draw("PESAME PLC PMC");

    TLegend *legend_Dmesons_Xi = new TLegend(0.45, 0.45, 0.85, 0.7);
    legend_Dmesons_Xi->SetBorderSize(0);
    legend_Dmesons_Xi->SetTextSize(0.04);
    legend_Dmesons_Xi->SetFillStyle(0);
    legend_Dmesons_Xi->AddEntry(h_Pt_DiMuon_Gen_Dplus_Xi, Form("%s,#Xi^{0}_{c} (%0.0f)", h_Pt_DiMuon_Gen_Dplus_Xi->GetTitle(), h_Pt_DiMuon_Gen_Dplus_Xi->GetEntries()));
    legend_Dmesons_Xi->AddEntry(h_Pt_DiMuon_Gen_Dzero_Xi, Form("%s,#Xi^{0}_{c} (%0.0f)", h_Pt_DiMuon_Gen_Dzero_Xi->GetTitle(), h_Pt_DiMuon_Gen_Dzero_Xi->GetEntries()));
    legend_Dmesons_Xi->AddEntry(h_Pt_DiMuon_Gen_Dstrange_Xi, Form("%s,#Xi^{0}_{c} (%0.0f)", h_Pt_DiMuon_Gen_Dstrange_Xi->GetTitle(), h_Pt_DiMuon_Gen_Dstrange_Xi->GetEntries()));
    legend_Dmesons_Xi->Draw();

    c_Dmesons_Xi->cd(2);

    h_M_DiMuon_Gen_Dplus_Xi->Draw("PE PLC PMC");
    h_M_DiMuon_Gen_Dzero_Xi->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen_Dstrange_Xi->Draw("PESAME PLC PMC");

    

    //----- Dmesons and J/psi --------------//

    TH1F *h_Pt_DiMuon_Gen_Dplus_Jpsi = (TH1F *)h_Pt_DiMuon_Gen[3]->Clone("h_Pt_DiMuon_Gen_Dplus_Jpsi");
    h_Pt_DiMuon_Gen_Dplus_Jpsi->Add(h_Pt_DiMuon_Gen[21]);

    TH1F *h_Pt_DiMuon_Gen_Dzero_Jpsi = (TH1F *)h_Pt_DiMuon_Gen[10]->Clone("h_Pt_DiMuon_Gen_Dzero_Jpsi");
    h_Pt_DiMuon_Gen_Dzero_Jpsi->Add(h_Pt_DiMuon_Gen[22]);

    TH1F *h_Pt_DiMuon_Gen_Dstrange_Jpsi = (TH1F *)h_Pt_DiMuon_Gen[17]->Clone("h_Pt_DiMuon_Gen_Dstrange_Jpsi");
    h_Pt_DiMuon_Gen_Dstrange_Jpsi->Add(h_Pt_DiMuon_Gen[23]);

    TH1F *h_Pt_DiMuon_Gen_Dmesons_Jpsi = (TH1F *)h_Pt_DiMuon_Gen_Dplus_Jpsi->Clone("h_Pt_DiMuon_Gen_Dmesons_Jpsi");
    h_Pt_DiMuon_Gen_Dmesons_Jpsi->Add(h_Pt_DiMuon_Gen_Dzero_Jpsi);
    h_Pt_DiMuon_Gen_Dmesons_Jpsi->Add(h_Pt_DiMuon_Gen_Dstrange_Jpsi);

    hist1D_graphic_opt(h_Pt_DiMuon_Gen[24], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[24]->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dplus_Jpsi, 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Dplus_Jpsi->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dzero_Jpsi, 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Dzero_Jpsi->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dstrange_Jpsi, 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Dstrange_Jpsi->GetEntries());

    TH1F *h_M_DiMuon_Gen_Dplus_Jpsi = (TH1F *)h_M_DiMuon_Gen[3]->Clone("h_M_DiMuon_Gen_Dplus_Jpsi");
    h_M_DiMuon_Gen_Dplus_Jpsi->Add(h_M_DiMuon_Gen[21]);

    TH1F *h_M_DiMuon_Gen_Dzero_Jpsi = (TH1F *)h_M_DiMuon_Gen[10]->Clone("h_M_DiMuon_Gen_Dzero_Jpsi");
    h_M_DiMuon_Gen_Dzero_Jpsi->Add(h_M_DiMuon_Gen[22]);

    TH1F *h_M_DiMuon_Gen_Dstrange_Jpsi = (TH1F *)h_M_DiMuon_Gen[17]->Clone("h_M_DiMuon_Gen_Dstrange_Jpsi");
    h_M_DiMuon_Gen_Dstrange_Jpsi->Add(h_M_DiMuon_Gen[23]);

    TH1F *h_M_DiMuon_Gen_Dmesons_Jpsi = (TH1F *)h_M_DiMuon_Gen_Dplus_Jpsi->Clone("h_M_DiMuon_Gen_Dmesons_Jpsi");
    h_M_DiMuon_Gen_Dmesons_Jpsi->Add(h_M_DiMuon_Gen_Dzero_Jpsi);
    h_M_DiMuon_Gen_Dmesons_Jpsi->Add(h_M_DiMuon_Gen_Dstrange_Jpsi);

    hist1D_graphic_opt(h_M_DiMuon_Gen[24], 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen[24]->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dplus_Jpsi, 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen_Dplus_Jpsi->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dzero_Jpsi, 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen_Dzero_Jpsi->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dstrange_Jpsi, 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen_Dstrange_Jpsi->GetEntries());

    TCanvas *c_Dmesons_Jpsi = canvas_noratio_divide2("c_Dmesons_Jpsi");
    c_Dmesons_Jpsi->cd(1);
    h_Pt_DiMuon_Gen[24]->Draw("PE PLC PMC");
    h_Pt_DiMuon_Gen_Dplus_Jpsi->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Dzero_Jpsi->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Dstrange_Jpsi->Draw("PESAME PLC PMC");

    TLegend *legend_Dmesons_Jpsi = new TLegend(0.45, 0.45, 0.85, 0.7);
    legend_Dmesons_Jpsi->SetBorderSize(0);
    legend_Dmesons_Jpsi->SetTextSize(0.04);
    legend_Dmesons_Jpsi->SetFillStyle(0);

    legend_Dmesons_Jpsi->AddEntry(h_Pt_DiMuon_Gen[24], Form("%s (%0.0f)", h_Pt_DiMuon_Gen[24]->GetTitle(), h_Pt_DiMuon_Gen[24]->GetEntries()));
    legend_Dmesons_Jpsi->AddEntry(h_Pt_DiMuon_Gen_Dplus_Jpsi, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dplus_Jpsi->GetTitle(), h_Pt_DiMuon_Gen_Dplus_Jpsi->GetEntries()));
    legend_Dmesons_Jpsi->AddEntry(h_Pt_DiMuon_Gen_Dzero_Jpsi, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dzero_Jpsi->GetTitle(), h_Pt_DiMuon_Gen_Dzero_Jpsi->GetEntries()));
    legend_Dmesons_Jpsi->AddEntry(h_Pt_DiMuon_Gen_Dstrange_Jpsi, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dstrange_Jpsi->GetTitle(), h_Pt_DiMuon_Gen_Dstrange_Jpsi->GetEntries()));
    legend_Dmesons_Jpsi->Draw();

    c_Dmesons_Jpsi->cd(2);
    h_M_DiMuon_Gen_Dplus_Jpsi->Draw("PESAME PMC PLC");
    h_M_DiMuon_Gen_Dzero_Jpsi->Draw("PESAME PMC PLC");
    h_M_DiMuon_Gen_Dstrange_Jpsi->Draw("PESAME PMC PLC");

    //-------Light Mesons------------//

    TH1F *h_Pt_DiMuon_Gen_Dplus_Dzero = (TH1F *)h_Pt_DiMuon_Gen[1]->Clone("h_Pt_DiMuon_Gen_Dplus_Dzero");
    h_Pt_DiMuon_Gen_Dplus_Dzero->Add(h_Pt_DiMuon_Gen[7]);

    TH1F *h_Pt_DiMuon_Gen_light_mesons = (TH1F *)h_Pt_DiMuon_Gen[0]->Clone("h_Pt_DiMuon_Gen_light_mesons");
    h_Pt_DiMuon_Gen_light_mesons->Add(h_Pt_DiMuon_Gen[8]);
    h_Pt_DiMuon_Gen_light_mesons->Add(h_Pt_DiMuon_Gen_Dplus_Dzero);

    hist1D_graphic_opt(h_Pt_DiMuon_Gen[0], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[0]->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen[8], 10, 20, kGreen + 2, 1. / h_Pt_DiMuon_Gen[8]->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dplus_Dzero, 10, 20, kAzure + 2, 1. / h_Pt_DiMuon_Gen_Dplus_Dzero->GetEntries());

    TH1F *h_M_DiMuon_Gen_Dplus_Dzero = (TH1F *)h_M_DiMuon_Gen[1]->Clone("h_M_DiMuon_Gen_Dplus_Dzero");
    h_M_DiMuon_Gen_Dplus_Dzero->Add(h_M_DiMuon_Gen[7]);

    TH1F *h_M_DiMuon_Gen_light_mesons = (TH1F *)h_M_DiMuon_Gen[0]->Clone("h_M_DiMuon_Gen_light_mesons");
    h_M_DiMuon_Gen_light_mesons->Add(h_M_DiMuon_Gen[8]);
    h_M_DiMuon_Gen_light_mesons->Add(h_M_DiMuon_Gen_Dplus_Dzero);

    hist1D_graphic_opt(h_M_DiMuon_Gen[0], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[0]->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen[8], 10, 20, kGreen + 2, 1. / h_Pt_DiMuon_Gen[8]->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dplus_Dzero, 10, 20, kAzure + 2, 1. / h_M_DiMuon_Gen_Dplus_Dzero->GetEntries());

    TCanvas *c_light_mesons = canvas_noratio_divide2("c_light_mesons");
    c_light_mesons->cd(1);

    h_Pt_DiMuon_Gen[0]->SetMinimum(1e-05);
    h_Pt_DiMuon_Gen[0]->Draw("PE PLC PMC");
    h_Pt_DiMuon_Gen[8]->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Dplus_Dzero->Draw("PESAME PLC PMC");

    TLegend *legend = new TLegend(0.45, 0.45, 0.85, 0.7);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);

    legend->AddEntry(h_Pt_DiMuon_Gen[0], Form("%s (%0.0f)", h_Pt_DiMuon_Gen[0]->GetTitle(), h_Pt_DiMuon_Gen[0]->GetEntries()));
    legend->AddEntry(h_Pt_DiMuon_Gen[8], Form("%s (%0.0f)", h_Pt_DiMuon_Gen[8]->GetTitle(), h_Pt_DiMuon_Gen[8]->GetEntries()));
    legend->AddEntry(h_Pt_DiMuon_Gen_Dplus_Dzero, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dplus_Dzero->GetTitle(), h_Pt_DiMuon_Gen_Dplus_Dzero->GetEntries()));
    legend->Draw();

    c_light_mesons->cd(2);
    h_M_DiMuon_Gen[0]->SetMinimum(1e-05);
    h_M_DiMuon_Gen[0]->Draw("PE PLC PMC");
    h_M_DiMuon_Gen[8]->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen_Dplus_Dzero->Draw("PESAME PLC PMC");

    //-------Heavy Mesons------------//

    TH1F *h_Pt_DiMuon_Gen_Dplus_Dstrange = (TH1F *)h_Pt_DiMuon_Gen[2]->Clone("h_Pt_DiMuon_Gen_Dplus_Dstrange");
    h_Pt_DiMuon_Gen_Dplus_Dstrange->Add(h_Pt_DiMuon_Gen[14]);

    TH1F *h_Pt_DiMuon_Gen_Dzero_Dstrange = (TH1F *)h_Pt_DiMuon_Gen[9]->Clone("h_Pt_DiMuon_Gen_Dzero_Dstrange");
    h_Pt_DiMuon_Gen_Dzero_Dstrange->Add(h_Pt_DiMuon_Gen[15]);

    TH1F *h_Pt_DiMuon_Gen_heavy_mesons = (TH1F *)h_Pt_DiMuon_Gen[16]->Clone("h_Pt_DiMuon_Gen_heavy_mesons");
    h_Pt_DiMuon_Gen_heavy_mesons->Add(h_Pt_DiMuon_Gen_Dplus_Dstrange);
    h_Pt_DiMuon_Gen_heavy_mesons->Add(h_Pt_DiMuon_Gen_Dzero_Dstrange);

    hist1D_graphic_opt(h_Pt_DiMuon_Gen[16], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[16]->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dzero_Dstrange, 10, 20, kGreen + 2, 1. / h_Pt_DiMuon_Gen_Dzero_Dstrange->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dplus_Dstrange, 10, 20, kAzure + 2, 1. / h_Pt_DiMuon_Gen_Dplus_Dstrange->GetEntries());

    TH1F *h_M_DiMuon_Gen_Dplus_Dstrange = (TH1F *)h_M_DiMuon_Gen[2]->Clone("h_M_DiMuon_Gen_Dplus_Dstrange");
    h_M_DiMuon_Gen_Dplus_Dstrange->Add(h_M_DiMuon_Gen[14]);

    TH1F *h_M_DiMuon_Gen_Dzero_Dstrange = (TH1F *)h_M_DiMuon_Gen[9]->Clone("h_M_DiMuon_Gen_Dzero_Dstrange");
    h_M_DiMuon_Gen_Dzero_Dstrange->Add(h_M_DiMuon_Gen[15]);

    TH1F *h_M_DiMuon_Gen_heavy_mesons = (TH1F *)h_M_DiMuon_Gen[16]->Clone("h_M_DiMuon_Gen_heavy_mesons");
    h_M_DiMuon_Gen_heavy_mesons->Add(h_M_DiMuon_Gen_Dplus_Dstrange);
    h_M_DiMuon_Gen_heavy_mesons->Add(h_M_DiMuon_Gen_Dzero_Dstrange);

    hist1D_graphic_opt(h_M_DiMuon_Gen[16], 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen[16]->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dzero_Dstrange, 10, 20, kGreen + 2, 1. / h_M_DiMuon_Gen_Dzero_Dstrange->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dplus_Dstrange, 10, 20, kAzure + 2, 1. / h_M_DiMuon_Gen_Dplus_Dstrange->GetEntries());

    TCanvas *c_heavy_meson = canvas_noratio_divide2("c_heavy_meson");
    c_heavy_meson->cd(1);

    h_Pt_DiMuon_Gen[16]->SetMinimum(1e-04);
    h_Pt_DiMuon_Gen[16]->Draw("PE PLC PMC");
    h_Pt_DiMuon_Gen_Dzero_Dstrange->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Dplus_Dstrange->Draw("PESAME PLC PMC");

    TLegend *legend_heavy_mesons = new TLegend(0.45, 0.45, 0.85, 0.7);
    legend_heavy_mesons->SetBorderSize(0);
    legend_heavy_mesons->SetTextSize(0.04);
    legend_heavy_mesons->SetFillStyle(0);

    legend_heavy_mesons->AddEntry(h_Pt_DiMuon_Gen[16], Form("%s (%0.0f)", h_Pt_DiMuon_Gen[16]->GetTitle(), h_Pt_DiMuon_Gen[16]->GetEntries()));
    legend_heavy_mesons->AddEntry(h_Pt_DiMuon_Gen_Dzero_Dstrange, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dzero_Dstrange->GetTitle(), h_Pt_DiMuon_Gen_Dzero_Dstrange->GetEntries()));
    legend_heavy_mesons->AddEntry(h_Pt_DiMuon_Gen_Dplus_Dstrange, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dplus_Dstrange->GetTitle(), h_Pt_DiMuon_Gen_Dplus_Dstrange->GetEntries()));
    legend_heavy_mesons->Draw();

    c_heavy_meson->cd(2);

    h_M_DiMuon_Gen[16]->SetMinimum(1e-04);
    h_M_DiMuon_Gen[16]->Draw("PE PLC PMC");
    h_M_DiMuon_Gen_Dzero_Dstrange->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen_Dplus_Dstrange->Draw("PESAME PLC PMC");

    //-------Lambda DMesons------------//

    TH1F *h_Pt_DiMuon_Gen_Dplus_Lambda = (TH1F *)h_Pt_DiMuon_Gen[4]->Clone("h_Pt_DiMuon_Gen_Dplus_Lambda");
    h_Pt_DiMuon_Gen_Dplus_Lambda->Add(h_Pt_DiMuon_Gen[28]);

    TH1F *h_Pt_DiMuon_Gen_Dzero_Lambda = (TH1F *)h_Pt_DiMuon_Gen[11]->Clone("h_Pt_DiMuon_Gen_Dzero_Lambda");
    h_Pt_DiMuon_Gen_Dzero_Lambda->Add(h_Pt_DiMuon_Gen[29]);

    TH1F *h_Pt_DiMuon_Gen_Dstrange_Lambda = (TH1F *)h_Pt_DiMuon_Gen[18]->Clone("h_Pt_DiMuon_Gen_Dstrange_Lambda");
    h_Pt_DiMuon_Gen_Dstrange_Lambda->Add(h_Pt_DiMuon_Gen[30]);

    TH1F *h_Pt_DiMuon_Gen_Dmesons_Lambda = (TH1F *)h_Pt_DiMuon_Gen_Dplus_Lambda->Clone("h_Pt_DiMuon_Gen_Dmesons_Lambda");
    h_Pt_DiMuon_Gen_Dmesons_Lambda->Add(h_Pt_DiMuon_Gen_Dzero_Lambda);
    h_Pt_DiMuon_Gen_Dmesons_Lambda->Add(h_Pt_DiMuon_Gen_Dstrange_Lambda);

    hist1D_graphic_opt(h_Pt_DiMuon_Gen[32], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[32]->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dplus_Lambda, 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Dplus_Lambda->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dzero_Lambda, 10, 20, kGreen + 2, 1. / h_Pt_DiMuon_Gen_Dzero_Lambda->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dstrange_Lambda, 10, 20, kAzure + 2, 1. / h_Pt_DiMuon_Gen_Dstrange_Lambda->GetEntries());

    TH1F *h_M_DiMuon_Gen_Dplus_Lambda = (TH1F *)h_M_DiMuon_Gen[4]->Clone("h_M_DiMuon_Gen_Dplus_Lambda");
    h_M_DiMuon_Gen_Dplus_Lambda->Add(h_M_DiMuon_Gen[28]);

    TH1F *h_M_DiMuon_Gen_Dzero_Lambda = (TH1F *)h_M_DiMuon_Gen[11]->Clone("h_M_DiMuon_Gen_Dzero_Lambda");
    h_M_DiMuon_Gen_Dzero_Lambda->Add(h_M_DiMuon_Gen[29]);

    TH1F *h_M_DiMuon_Gen_Dstrange_Lambda = (TH1F *)h_M_DiMuon_Gen[18]->Clone("h_M_DiMuon_Gen_Dstrange_Lambda");
    h_M_DiMuon_Gen_Dstrange_Lambda->Add(h_M_DiMuon_Gen[30]);

    TH1F *h_M_DiMuon_Gen_Dmesons_Lambda = (TH1F *)h_M_DiMuon_Gen_Dplus_Lambda->Clone("h_M_DiMuon_Gen_Dmesons_Lambda");
    h_M_DiMuon_Gen_Dmesons_Lambda->Add(h_M_DiMuon_Gen_Dzero_Lambda);
    h_M_DiMuon_Gen_Dmesons_Lambda->Add(h_M_DiMuon_Gen_Dstrange_Lambda);

    hist1D_graphic_opt(h_M_DiMuon_Gen[32], 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen[32]->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dplus_Lambda, 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen_Dplus_Lambda->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dzero_Lambda, 10, 20, kGreen + 2, 1. / h_M_DiMuon_Gen_Dzero_Lambda->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dstrange_Lambda, 10, 20, kAzure + 2, 1. / h_M_DiMuon_Gen_Dstrange_Lambda->GetEntries());

    TCanvas *c_lambda_Dmeson = canvas_noratio_divide2("c_lambda_Dmeson");
    c_lambda_Dmeson->cd(1);

    printf("%s\n", h_Pt_DiMuon_Gen[32]->GetName());
    printf("%s\n", h_Pt_DiMuon_Gen[4]->GetName());
    printf("%s\n", h_Pt_DiMuon_Gen[28]->GetName());

    printf("%s\n", h_Pt_DiMuon_Gen[11]->GetName());
    printf("%s\n", h_Pt_DiMuon_Gen[29]->GetName());

    printf("%s\n", h_Pt_DiMuon_Gen[18]->GetName());
    printf("%s\n", h_Pt_DiMuon_Gen[30]->GetName());

    // h_Pt_DiMuon_Gen[32]->Draw("PE");
    h_Pt_DiMuon_Gen_Dplus_Lambda->SetMinimum(1e-04);
    h_Pt_DiMuon_Gen_Dplus_Lambda->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Dzero_Lambda->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Dstrange_Lambda->Draw("PESAME PLC PMC");

    TLegend *legend_Lambda_Dmesons = new TLegend(0.45, 0.45, 0.85, 0.7);
    legend_Lambda_Dmesons->SetBorderSize(0);
    legend_Lambda_Dmesons->SetTextSize(0.04);
    legend_Lambda_Dmesons->SetFillStyle(0);

    // legend_Lambda_Dmesons->AddEntry(h_Pt_DiMuon_Gen[32], Form("%s (%0.0f)", h_Pt_DiMuon_Gen[32]->GetTitle(), h_Pt_DiMuon_Gen[32]->GetEntries()));
    legend_Lambda_Dmesons->AddEntry(h_Pt_DiMuon_Gen_Dplus_Lambda, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dplus_Lambda->GetTitle(), h_Pt_DiMuon_Gen_Dplus_Lambda->GetEntries()));
    legend_Lambda_Dmesons->AddEntry(h_Pt_DiMuon_Gen_Dzero_Lambda, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dzero_Lambda->GetTitle(), h_Pt_DiMuon_Gen_Dzero_Lambda->GetEntries()));
    legend_Lambda_Dmesons->AddEntry(h_Pt_DiMuon_Gen_Dstrange_Lambda, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dstrange_Lambda->GetTitle(), h_Pt_DiMuon_Gen_Dstrange_Lambda->GetEntries()));
    legend_Lambda_Dmesons->Draw();

    c_lambda_Dmeson->cd(2);
    h_M_DiMuon_Gen_Dplus_Lambda->SetMinimum(1e-04);
    h_M_DiMuon_Gen_Dplus_Lambda->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen_Dzero_Lambda->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen_Dstrange_Lambda->Draw("PESAME PLC PMC");

    //---- Comparison of the different contributions-----//
    TCanvas *c_comparison = canvas_noratio_divide2("c_comparison");
    c_comparison->cd(1);

    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dmesons_Jpsi, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_light_mesons, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_heavy_mesons, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dmesons_Lambda, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Barions, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dmesons_Xi, 15, 20, kRed + 2, 1);

    h_Pt_DiMuon_Gen_Dmesons_Jpsi->Draw("PE PLC PMC");
    h_Pt_DiMuon_Gen_light_mesons->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_heavy_mesons->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Dmesons_Lambda->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Barions->Draw("PESAME PLC PMC");
    h_Pt_DiMuon_Gen_Dmesons_Xi->Draw("PESAME PLC PMC");

    TLegend *legend_comparison = new TLegend(0.45, 0.45, 0.85, 0.7);
    legend_comparison->SetBorderSize(0);
    legend_comparison->SetTextSize(0.04);
    legend_comparison->SetFillStyle(0);
    legend_comparison->AddEntry(h_Pt_DiMuon_Gen_Dmesons_Jpsi, Form("#mu #leftarrow D, #mu #leftarrow J/#psi  (%0.0f)", h_Pt_DiMuon_Gen_Dmesons_Jpsi->GetEntries()));
    legend_comparison->AddEntry(h_Pt_DiMuon_Gen_light_mesons, Form("#mu #leftarrow D^{0}/D^{#plus},#mu #leftarrow D^{0}/D^{#plus} (%0.2e)", h_Pt_DiMuon_Gen_light_mesons->GetEntries()));
    legend_comparison->AddEntry(h_Pt_DiMuon_Gen_heavy_mesons, Form("#mu #leftarrow D,#mu #leftarrow D^{#plus}_{s} (%0.0f)", h_Pt_DiMuon_Gen_heavy_mesons->GetEntries()));
    legend_comparison->AddEntry(h_Pt_DiMuon_Gen_Dmesons_Lambda, Form("#mu #leftarrow D/#mu #leftarrow #Lambda^{#plus}_{c} (%0.0f)", h_Pt_DiMuon_Gen_Dmesons_Lambda->GetEntries()));
    legend_comparison->AddEntry(h_Pt_DiMuon_Gen_Barions, Form("#mu #leftarrow Barions (%0.0f)", h_Pt_DiMuon_Gen_Barions->GetEntries()));
    legend_comparison->AddEntry(h_Pt_DiMuon_Gen_Dmesons_Xi, Form("#mu #leftarrow D, #mu #leftarrow #Xi (%0.0f)", h_Pt_DiMuon_Gen_Dmesons_Xi->GetEntries()));
    legend_comparison->Draw();

    c_comparison->cd(2);

    hist1D_graphic_opt(h_M_DiMuon_Gen_Dmesons_Jpsi, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_M_DiMuon_Gen_light_mesons, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_M_DiMuon_Gen_heavy_mesons, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dmesons_Lambda, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_M_DiMuon_Gen_Barions, 15, 20, kRed + 2, 1);
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dmesons_Xi, 15, 20, kRed + 2, 1);

    h_M_DiMuon_Gen_Dmesons_Jpsi->Draw("PE PLC PMC");
    h_M_DiMuon_Gen_light_mesons->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen_heavy_mesons->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen_Dmesons_Lambda->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen_Barions->Draw("PESAME PLC PMC");
    h_M_DiMuon_Gen_Dmesons_Xi->Draw("PESAME PLC PMC");
}
