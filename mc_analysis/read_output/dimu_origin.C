#include "/home/michele_pennisi/cernbox/common_include.h"

void dimu_origin()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile *fIn = new TFile("test/HF_MC_output_Hist_294009.root", "READ");

    TH3F *h_Pdg1Pdg2Pt_DiMuon_Gen = (TH3F *)fIn->Get("DiMuon_Gen/h_Pdg1Pdg2Pt_DiMuon_Gen");
    TH3F *h_Pdg1Pdg2M_DiMuon_Gen = (TH3F *)fIn->Get("DiMuon_Gen/h_Pdg1Pdg2M_DiMuon_Gen");

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

    TCanvas *c_2d_mesons = new TCanvas("c_2d_mesons", "c_2d_mesons", 800, 800);
    c_2d_mesons->cd();
    gPad->SetGridx();
    gPad->SetGridy();
    h_Pdg1Pdg2_DiMuon_Gen_mesons->Draw("textCOLZSAME");

    TCanvas *c_2d_barions = new TCanvas("c_2d_barions", "c_2d_barions", 800, 800);
    c_2d_barions->cd();
    gPad->SetGridx();
    gPad->SetGridy();

    h_Pdg1Pdg2_DiMuon_Gen_barions->Draw("textCOLZSAME");

    TCanvas *c_2d_mesons_barions = new TCanvas("c_2d_mesons_barions", "c_2d_mesons_barions", 800, 800);
    c_2d_mesons_barions->cd();
    gPad->SetGridx();
    gPad->SetGridy();

    h_Pdg1Pdg2_DiMuon_Gen_mesons_barions->Draw("textCOLZSAME");

    TCanvas *c_2d_barions_mesons = new TCanvas("c_2d_barions_mesons", "c_2d_barions_mesons", 800, 800);
    c_2d_barions_mesons->cd();
    gPad->SetGridx();
    gPad->SetGridy();

    h_Pdg1Pdg2_DiMuon_Gen_barions_mesons->Draw("textCOLZSAME");

    const Int_t n_charm_meson_origin = 10;
    Int_t charm_meson_origin[n_charm_meson_origin] = {410, 420, 430, 440, 450, 4120, 4130, 4140, 4230, 4240};
    TString charm_meson_origin_name[n_charm_meson_origin] = {"Dplus", "Dzero", "Dstrange", "Jpsi", " ", "Lambda", "Xi_c0", " ", "Xi_c+"};
    TString charm_meson_origin_title[n_charm_meson_origin] = {"D^{#plus}", "D^{0}", "D^{#plus}_{s}", " ", "J/#psi", "#Lambda^{+}_{c}", "#Xi^{0}_{c}", " ", "#Xi^{+}_{c}"};
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
    //-------Light Mesons------------//
    hist1D_graphic_opt(h_Pt_DiMuon_Gen[0], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[0]->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen[8], 10, 20, kGreen + 2, 1. / h_Pt_DiMuon_Gen[8]->GetEntries());


    TH1F *h_Pt_DiMuon_Gen_Dplus_Dzero = (TH1F *)h_Pt_DiMuon_Gen[1]->Clone("h_Pt_DiMuon_Gen_Dplus_Dzero");
    h_Pt_DiMuon_Gen_Dplus_Dzero->Add(h_Pt_DiMuon_Gen[7]);

    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dplus_Dzero, 10, 20, kAzure + 2, 1. / h_Pt_DiMuon_Gen_Dplus_Dzero->GetEntries());
    
    hist1D_graphic_opt(h_M_DiMuon_Gen[0], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[0]->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen[8], 10, 20, kGreen + 2, 1. / h_Pt_DiMuon_Gen[8]->GetEntries());

    TH1F *h_M_DiMuon_Gen_Dplus_Dzero = (TH1F *)h_M_DiMuon_Gen[1]->Clone("h_M_DiMuon_Gen_Dplus_Dzero");
    h_M_DiMuon_Gen_Dplus_Dzero->Add(h_M_DiMuon_Gen[7]);

    hist1D_graphic_opt(h_M_DiMuon_Gen_Dplus_Dzero, 10, 20, kAzure + 2, 1. / h_M_DiMuon_Gen_Dplus_Dzero->GetEntries());

    TCanvas *c_light_mesons = canvas_noratio_divide2("c_light_mesons");
    c_light_mesons->cd(1);

    h_Pt_DiMuon_Gen[0]->SetMinimum(1e-05);
    h_Pt_DiMuon_Gen[0]->Draw("PE");
    h_Pt_DiMuon_Gen[8]->Draw("PESAME");
    h_Pt_DiMuon_Gen_Dplus_Dzero->Draw("PESAME");

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
    h_M_DiMuon_Gen[0]->Draw("PE");
    h_M_DiMuon_Gen[8]->Draw("PESAME");
    h_M_DiMuon_Gen_Dplus_Dzero->Draw("PESAME");

    //-------Heavy Mesons------------//
    hist1D_graphic_opt(h_Pt_DiMuon_Gen[16], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[16]->GetEntries());

    TH1F *h_Pt_DiMuon_Gen_Dplus_Dstrange = (TH1F *)h_Pt_DiMuon_Gen[2]->Clone("h_Pt_DiMuon_Gen_Dplus_Dstrange");
    h_Pt_DiMuon_Gen_Dplus_Dstrange->Add(h_Pt_DiMuon_Gen[14]);

    TH1F *h_Pt_DiMuon_Gen_Dzero_Dstrange = (TH1F *)h_Pt_DiMuon_Gen[9]->Clone("h_Pt_DiMuon_Gen_Dzero_Dstrange");
    h_Pt_DiMuon_Gen_Dzero_Dstrange->Add(h_Pt_DiMuon_Gen[15]);

    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dzero_Dstrange, 10, 20, kGreen + 2, 1. / h_Pt_DiMuon_Gen_Dzero_Dstrange->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dplus_Dstrange, 10, 20, kAzure + 2, 1. / h_Pt_DiMuon_Gen_Dplus_Dstrange->GetEntries());

    hist1D_graphic_opt(h_M_DiMuon_Gen[16], 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen[16]->GetEntries());

    TH1F *h_M_DiMuon_Gen_Dplus_Dstrange = (TH1F *)h_M_DiMuon_Gen[2]->Clone("h_M_DiMuon_Gen_Dplus_Dstrange");
    h_M_DiMuon_Gen_Dplus_Dstrange->Add(h_M_DiMuon_Gen[14]);

    TH1F *h_M_DiMuon_Gen_Dzero_Dstrange = (TH1F *)h_M_DiMuon_Gen[9]->Clone("h_M_DiMuon_Gen_Dzero_Dstrange");
    h_M_DiMuon_Gen_Dzero_Dstrange->Add(h_M_DiMuon_Gen[15]);

    hist1D_graphic_opt(h_M_DiMuon_Gen_Dzero_Dstrange, 10, 20, kGreen + 2, 1. / h_M_DiMuon_Gen_Dzero_Dstrange->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dplus_Dstrange, 10, 20, kAzure + 2, 1. / h_M_DiMuon_Gen_Dplus_Dstrange->GetEntries());

    TCanvas *c_heavy_meson = canvas_noratio_divide2("c_heavy_meson");
    c_heavy_meson->cd(1);

    h_Pt_DiMuon_Gen[16]->SetMinimum(1e-04);
    h_Pt_DiMuon_Gen[16]->Draw("PE");
    h_Pt_DiMuon_Gen_Dzero_Dstrange->Draw("PESAME");
    h_Pt_DiMuon_Gen_Dplus_Dstrange->Draw("PESAME");

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
    h_M_DiMuon_Gen[16]->Draw("PE");
    h_M_DiMuon_Gen_Dzero_Dstrange->Draw("PESAME");
    h_M_DiMuon_Gen_Dplus_Dstrange->Draw("PESAME");

    //-------Barion Mesons------------//

    hist1D_graphic_opt(h_Pt_DiMuon_Gen[32], 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen[32]->GetEntries());

    TH1F *h_Pt_DiMuon_Gen_Dplus_Lambda = (TH1F *)h_Pt_DiMuon_Gen[4]->Clone("h_Pt_DiMuon_Gen_Dplus_Lambda");
    h_Pt_DiMuon_Gen_Dplus_Lambda->Add(h_Pt_DiMuon_Gen[28]);

    TH1F *h_Pt_DiMuon_Gen_Dzero_Lambda = (TH1F *)h_Pt_DiMuon_Gen[11]->Clone("h_Pt_DiMuon_Gen_Dzero_Lambda");
    h_Pt_DiMuon_Gen_Dzero_Lambda->Add(h_Pt_DiMuon_Gen[29]);

    TH1F *h_Pt_DiMuon_Gen_Dstrange_Lambda = (TH1F *)h_Pt_DiMuon_Gen[18]->Clone("h_Pt_DiMuon_Gen_Dstrange_Lambda");
    h_Pt_DiMuon_Gen_Dstrange_Lambda->Add(h_Pt_DiMuon_Gen[30]);
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dplus_Lambda, 10, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Dplus_Lambda->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dzero_Lambda, 10, 20, kGreen + 2, 1. / h_Pt_DiMuon_Gen_Dzero_Lambda->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Dstrange_Lambda, 10, 20, kAzure + 2, 1. / h_Pt_DiMuon_Gen_Dstrange_Lambda->GetEntries());

    hist1D_graphic_opt(h_M_DiMuon_Gen[32], 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen[32]->GetEntries());

    TH1F *h_M_DiMuon_Gen_Dplus_Lambda = (TH1F *)h_M_DiMuon_Gen[4]->Clone("h_M_DiMuon_Gen_Dplus_Lambda");
    h_M_DiMuon_Gen_Dplus_Lambda->Add(h_M_DiMuon_Gen[28]);

    TH1F *h_M_DiMuon_Gen_Dzero_Lambda = (TH1F *)h_M_DiMuon_Gen[11]->Clone("h_M_DiMuon_Gen_Dzero_Lambda");
    h_M_DiMuon_Gen_Dzero_Lambda->Add(h_M_DiMuon_Gen[29]);

    TH1F *h_M_DiMuon_Gen_Dstrange_Lambda = (TH1F *)h_M_DiMuon_Gen[18]->Clone("h_M_DiMuon_Gen_Dstrange_Lambda");
    h_M_DiMuon_Gen_Dstrange_Lambda->Add(h_M_DiMuon_Gen[30]);

    hist1D_graphic_opt(h_M_DiMuon_Gen_Dplus_Lambda, 10, 20, kRed + 2, 1. / h_M_DiMuon_Gen_Dplus_Lambda->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dzero_Lambda, 10, 20, kGreen + 2, 1. / h_M_DiMuon_Gen_Dzero_Lambda->GetEntries());
    hist1D_graphic_opt(h_M_DiMuon_Gen_Dstrange_Lambda, 10, 20, kAzure + 2, 1. / h_M_DiMuon_Gen_Dstrange_Lambda->GetEntries());

    TCanvas *c_barion_meson = canvas_noratio_divide2("c_barion_meson");
    c_barion_meson->cd(1);
    
    printf("%s\n", h_Pt_DiMuon_Gen[32]->GetName());
    printf("%s\n", h_Pt_DiMuon_Gen[4]->GetName());
    printf("%s\n", h_Pt_DiMuon_Gen[28]->GetName());

    printf("%s\n", h_Pt_DiMuon_Gen[11]->GetName());
    printf("%s\n", h_Pt_DiMuon_Gen[29]->GetName());

    printf("%s\n", h_Pt_DiMuon_Gen[18]->GetName());
    printf("%s\n", h_Pt_DiMuon_Gen[30]->GetName());


    // h_Pt_DiMuon_Gen[32]->Draw("PE");
    h_Pt_DiMuon_Gen_Dplus_Lambda->SetMinimum(1e-04);
    h_Pt_DiMuon_Gen_Dplus_Lambda->Draw("PESAME");
    h_Pt_DiMuon_Gen_Dzero_Lambda->Draw("PESAME");
    h_Pt_DiMuon_Gen_Dstrange_Lambda->Draw("PESAME");

    TLegend *legend_barion_mesons = new TLegend(0.45, 0.45, 0.85, 0.7);
    legend_barion_mesons->SetBorderSize(0);
    legend_barion_mesons->SetTextSize(0.04);
    legend_barion_mesons->SetFillStyle(0);

    // legend_barion_mesons->AddEntry(h_Pt_DiMuon_Gen[32], Form("%s (%0.0f)", h_Pt_DiMuon_Gen[32]->GetTitle(), h_Pt_DiMuon_Gen[32]->GetEntries()));
    legend_barion_mesons->AddEntry(h_Pt_DiMuon_Gen_Dplus_Lambda, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dplus_Lambda->GetTitle(), h_Pt_DiMuon_Gen_Dplus_Lambda->GetEntries()));
    legend_barion_mesons->AddEntry(h_Pt_DiMuon_Gen_Dzero_Lambda, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dzero_Lambda->GetTitle(), h_Pt_DiMuon_Gen_Dzero_Lambda->GetEntries()));
    legend_barion_mesons->AddEntry(h_Pt_DiMuon_Gen_Dstrange_Lambda, Form("%s (%0.0f)", h_Pt_DiMuon_Gen_Dstrange_Lambda->GetTitle(), h_Pt_DiMuon_Gen_Dstrange_Lambda->GetEntries()));
    legend_barion_mesons->Draw();

    c_barion_meson->cd(2);
    h_M_DiMuon_Gen_Dplus_Lambda->SetMinimum(1e-04);
    h_M_DiMuon_Gen_Dplus_Lambda->Draw("PESAME");
    h_M_DiMuon_Gen_Dzero_Lambda->Draw("PESAME");
    h_M_DiMuon_Gen_Dstrange_Lambda->Draw("PESAME");

}