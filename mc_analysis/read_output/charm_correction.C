#include "/home/michele_pennisi/cernbox/common_include.h"

double FuncPtMass(double *x, double *par);

void charm_correction()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TFile *fIn = new TFile("test/HF_MC_output_Hist_294009.root", "UPDATED");

    TH2F *h_PtM_DiMuon_Gen_DQcut = (TH2F *)fIn->Get("DiMuon_Gen/h_PtM_DiMuon_Gen_DQcut_Charm");
    TH1F *h_Pt_DiMuon_Gen_DQcut = (TH1F *)h_PtM_DiMuon_Gen_DQcut->ProjectionX();

    TH3F *h_Pdg1Pdg2Pt_DiMuon_Gen_Original = (TH3F *)fIn->Get("DiMuon_Gen/h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut");
    TH3F *h_Pdg1Pdg2M_DiMuon_Gen_Original = (TH3F *)fIn->Get("DiMuon_Gen/h_Pdg1Pdg2M_DiMuon_Gen_DQcut");
    TH3F *h_Pdg1Pdg2Pt_DiMuon_Gen_Corrected = (TH3F *)fIn->Get("DiMuon_Charm_corrected/h_Pdg1Pdg2Pt_DiMuon_Gen_DQcut_Charm_corrected");
    TH3F *h_Pdg1Pdg2M_DiMuon_Gen_Corrected = (TH3F *)fIn->Get("DiMuon_Charm_corrected/h_Pdg1Pdg2M_DiMuon_Gen_DQcut_Charm_corrected");

    h_Pdg1Pdg2Pt_DiMuon_Gen_Original->GetXaxis()->SetRangeUser(400, 500);
    h_Pdg1Pdg2Pt_DiMuon_Gen_Original->GetYaxis()->SetRangeUser(4000, 5000);
    TH2F *h_Pdg1Pdg2_DiMuon_Gen_mesons_Original = (TH2F *)h_Pdg1Pdg2Pt_DiMuon_Gen_Original->Project3D("xy");
    h_Pdg1Pdg2_DiMuon_Gen_mesons_Original->SetName("charm_mesons_pdg1_pdg2_Original");
    h_Pdg1Pdg2_DiMuon_Gen_mesons_Original->Draw("textCOLZ");

    new TCanvas();

    h_Pdg1Pdg2Pt_DiMuon_Gen_Corrected->GetXaxis()->SetRangeUser(400, 500);
    h_Pdg1Pdg2Pt_DiMuon_Gen_Corrected->GetYaxis()->SetRangeUser(4000, 5000);
    TH2F *h_Pdg1Pdg2_DiMuon_Gen_mesons_Corrected = (TH2F *)h_Pdg1Pdg2Pt_DiMuon_Gen_Corrected->Project3D("xy");
    h_Pdg1Pdg2_DiMuon_Gen_mesons_Corrected->SetName("charm_mesons_pdg1_pdg2_Corrected");
    h_Pdg1Pdg2_DiMuon_Gen_mesons_Corrected->Draw("textCOLZ");

    const Int_t n_charm_hadron_origin = 10;
    Int_t charm_meson_origin[n_charm_hadron_origin] = {410, 420, 430, 440, 450, 4120, 4130, 4140, 4230, 4240};
    TString charm_meson_origin_name[n_charm_hadron_origin] = {"Dplus", "Dzero", "Dstrange", "Jpsi", " ", "Lambda", "Xi_c0", " ", "Xi_c+"};
    TString charm_meson_origin_title[n_charm_hadron_origin] = {"D^{#plus}", "D^{0}", "D^{#plus}_{s}", "J/#psi", " ", "#Lambda^{+}_{c}", "#Xi^{0}_{c}", " ", "#Xi^{+}_{c}"};
    Int_t original_pdg[n_charm_hadron_origin] = {411, 421, 431, 443, 99999, 4122, 4132, 99999, 4232};

    Double_t BR_charm_hadrons2mu_PYTHIA[n_charm_hadron_origin] = {16.5, 6.45, 7.5, 5.9, 99999, 4.5, 2.5, 99999, 3.5};
    Double_t BR_charm_hadrons2mu_MEAS[n_charm_hadron_origin] = {17.6, 6.80, 6.33, 5.96, 99999, 3.95, 2.5, 99999, 3.5};

    Double_t Frag_charm_hadrons_PYTHIA[n_charm_hadron_origin] = {29.3, 56.1, 9.59, 0.4, 99999, 3.8, 0.49, 99999, 0.49};
    Double_t Frag_charm_hadrons_MEAS[n_charm_hadron_origin] = {19.1, 38.2, 6.1, 0.37, 99999, 16.8, 9.9, 999999, 9.6};

    TH1F *h_Pt_DiMuon_Gen_Original[49];
    TH1F *h_M_DiMuon_Gen_Original[49];

    TH1F *h_Pt_DiMuon_Gen_Corrected[49];
    TH1F *h_M_DiMuon_Gen_Corrected[49];

    TH1F *h_Pt_DiMuon_Gen_Sum;
    TH1F *h_Pt_DiMuon_Gen_Sum_Corrected;
    TCanvas *c[49];

    Int_t index = 0;

    for (Int_t i_charm_meson_origin_mu1 = 0; i_charm_meson_origin_mu1 < n_charm_hadron_origin - 1; i_charm_meson_origin_mu1++)
    {
        if ((charm_meson_origin[i_charm_meson_origin_mu1] == 450 && charm_meson_origin[i_charm_meson_origin_mu1 + 1] == 4120) || (charm_meson_origin[i_charm_meson_origin_mu1] == 4140 && charm_meson_origin[i_charm_meson_origin_mu1 + 1] == 4230))
            continue;
        h_Pdg1Pdg2Pt_DiMuon_Gen_Original->GetXaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu1], charm_meson_origin[i_charm_meson_origin_mu1 + 1]);
        h_Pdg1Pdg2M_DiMuon_Gen_Original->GetXaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu1], charm_meson_origin[i_charm_meson_origin_mu1 + 1]);

        h_Pdg1Pdg2Pt_DiMuon_Gen_Corrected->GetXaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu1], charm_meson_origin[i_charm_meson_origin_mu1 + 1]);
        h_Pdg1Pdg2M_DiMuon_Gen_Corrected->GetXaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu1], charm_meson_origin[i_charm_meson_origin_mu1 + 1]);

        for (Int_t i_charm_meson_origin_mu2 = 0; i_charm_meson_origin_mu2 < n_charm_hadron_origin - 1; i_charm_meson_origin_mu2++)
        {
            if ((charm_meson_origin[i_charm_meson_origin_mu2] == 450 && charm_meson_origin[i_charm_meson_origin_mu2 + 1] == 4120) || (charm_meson_origin[i_charm_meson_origin_mu2] == 4140 && charm_meson_origin[i_charm_meson_origin_mu2 + 1] == 4230))
                continue;
            cout << "Bound1 (" << charm_meson_origin[i_charm_meson_origin_mu1] << "," << charm_meson_origin[i_charm_meson_origin_mu1 + 1] << ")  ------> " << charm_meson_origin_name[i_charm_meson_origin_mu1].Data() << " original pdg: " << original_pdg[i_charm_meson_origin_mu1] << endl;
            cout << "Bound2 (" << charm_meson_origin[i_charm_meson_origin_mu2] << "," << charm_meson_origin[i_charm_meson_origin_mu2 + 1] << ")  ------> " << charm_meson_origin_name[i_charm_meson_origin_mu2].Data() << " original pdg: " << original_pdg[i_charm_meson_origin_mu2] << endl;
            Double_t correction_mu1 = (BR_charm_hadrons2mu_MEAS[i_charm_meson_origin_mu1] * Frag_charm_hadrons_MEAS[i_charm_meson_origin_mu1]) / (BR_charm_hadrons2mu_PYTHIA[i_charm_meson_origin_mu1] * Frag_charm_hadrons_PYTHIA[i_charm_meson_origin_mu1]);

            Double_t correction_mu2 = (BR_charm_hadrons2mu_MEAS[i_charm_meson_origin_mu2] * Frag_charm_hadrons_MEAS[i_charm_meson_origin_mu2]) / (BR_charm_hadrons2mu_PYTHIA[i_charm_meson_origin_mu2] * Frag_charm_hadrons_PYTHIA[i_charm_meson_origin_mu2]);

            Double_t correction_dimu = correction_mu1 * correction_mu2;

            cout << "Correction Dimuon: " << correction_dimu << endl;
            cout << " ==================================== " << endl;

            //---Pt---//
            h_Pdg1Pdg2Pt_DiMuon_Gen_Original->GetYaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu2], charm_meson_origin[i_charm_meson_origin_mu2 + 1]);
            h_Pt_DiMuon_Gen_Original[index] = (TH1F *)h_Pdg1Pdg2Pt_DiMuon_Gen_Original->Project3D("NUFze");
            h_Pt_DiMuon_Gen_Original[index]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
            h_Pt_DiMuon_Gen_Original[index]->SetName(Form("charm_%s_%s_pt_%d_Original", charm_meson_origin_name[i_charm_meson_origin_mu1].Data(), charm_meson_origin_name[i_charm_meson_origin_mu2].Data(), index));
            h_Pt_DiMuon_Gen_Original[index]->SetTitle(Form("#mu #leftarrow %s, #mu #leftarrow %s", charm_meson_origin_title[i_charm_meson_origin_mu1].Data(), charm_meson_origin_title[i_charm_meson_origin_mu2].Data()));

            h_Pdg1Pdg2Pt_DiMuon_Gen_Corrected->GetYaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu2], charm_meson_origin[i_charm_meson_origin_mu2 + 1]);
            h_Pt_DiMuon_Gen_Corrected[index] = (TH1F *)h_Pdg1Pdg2Pt_DiMuon_Gen_Corrected->Project3D("z");
            h_Pt_DiMuon_Gen_Corrected[index]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
            h_Pt_DiMuon_Gen_Corrected[index]->SetName(Form("charm_%s_%s_pt_%d_Corrected", charm_meson_origin_name[i_charm_meson_origin_mu1].Data(), charm_meson_origin_name[i_charm_meson_origin_mu2].Data(), index));
            h_Pt_DiMuon_Gen_Corrected[index]->SetTitle(Form("#mu #leftarrow %s, #mu #leftarrow %s", charm_meson_origin_title[i_charm_meson_origin_mu1].Data(), charm_meson_origin_title[i_charm_meson_origin_mu2].Data()));

            if (index == 0)
            {
                h_Pt_DiMuon_Gen_Sum = (TH1F *)h_Pt_DiMuon_Gen_Original[index]->Clone("h_Pt_DiMuon_Gen_Sum");
                h_Pt_DiMuon_Gen_Sum_Corrected = (TH1F *)h_Pt_DiMuon_Gen_Corrected[index]->Clone("h_Pt_DiMuon_Gen_Sum");
            }
            else
            {
                h_Pt_DiMuon_Gen_Sum->Add(h_Pt_DiMuon_Gen_Original[index]);
                h_Pt_DiMuon_Gen_Sum_Corrected->Add(h_Pt_DiMuon_Gen_Corrected[index]);
            }

            //---M---//
            h_Pdg1Pdg2M_DiMuon_Gen_Original->GetYaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu2], charm_meson_origin[i_charm_meson_origin_mu2 + 1]);
            h_M_DiMuon_Gen_Original[index] = (TH1F *)h_Pdg1Pdg2M_DiMuon_Gen_Original->ProjectionZ();
            h_M_DiMuon_Gen_Original[index]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
            h_M_DiMuon_Gen_Original[index]->SetName(Form("charm_%s_%s_M_%d_Original", charm_meson_origin_name[i_charm_meson_origin_mu1].Data(), charm_meson_origin_name[i_charm_meson_origin_mu2].Data(), index));
            h_M_DiMuon_Gen_Original[index]->SetTitle(Form("#mu #leftarrow %s, #mu #leftarrow %s", charm_meson_origin_title[i_charm_meson_origin_mu1].Data(), charm_meson_origin_title[i_charm_meson_origin_mu2].Data()));

            h_Pdg1Pdg2M_DiMuon_Gen_Corrected->GetYaxis()->SetRangeUser(charm_meson_origin[i_charm_meson_origin_mu2], charm_meson_origin[i_charm_meson_origin_mu2 + 1]);
            h_M_DiMuon_Gen_Corrected[index] = (TH1F *)h_Pdg1Pdg2M_DiMuon_Gen_Corrected->Project3D("z");
            h_M_DiMuon_Gen_Corrected[index]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
            h_M_DiMuon_Gen_Corrected[index]->SetName(Form("charm_%s_%s_M_%d_Corrected", charm_meson_origin_name[i_charm_meson_origin_mu1].Data(), charm_meson_origin_name[i_charm_meson_origin_mu2].Data(), index));
            h_M_DiMuon_Gen_Corrected[index]->SetTitle(Form("#mu #leftarrow %s, #mu #leftarrow %s", charm_meson_origin_title[i_charm_meson_origin_mu1].Data(), charm_meson_origin_title[i_charm_meson_origin_mu2].Data()));
            cout << "index: " << index << endl;
            cout << "h_Pt_DiMuon_Gen_Original[index]->GetEntries() " << h_Pt_DiMuon_Gen_Original[index]->GetEntries() << endl;
            cout << "h_Pt_DiMuon_Gen_Corrected[index]->GetEntries() " << h_Pt_DiMuon_Gen_Corrected[index]->GetEntries() << endl;
            h_Pt_DiMuon_Gen_Original[index]->Sumw2();
            h_Pt_DiMuon_Gen_Corrected[index]->Sumw2();
            // hist1D_graphic_opt(h_Pt_DiMuon_Gen_Corrected[index], 1, 20, kBlack, 1. / h_Pt_DiMuon_Gen_Corrected[index]->GetEntries());
            // hist1D_graphic_opt(h_Pt_DiMuon_Gen_Original[index], 1, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Original[index]->GetEntries());

            TH1F *ratio = (TH1F *)h_Pt_DiMuon_Gen_Corrected[index]->Clone("ratio");
            ratio->Divide(h_Pt_DiMuon_Gen_Original[index]);

            // c[index] = new TCanvas(Form("c_%d", index), Form("c_%d", index), 1000, 1500);
            // c[index]->Divide(1, 2);
            // c[index]->cd(1);
            // h_Pt_DiMuon_Gen_Original[index]->Draw("PESAME");
            // h_Pt_DiMuon_Gen_Corrected[index]->Draw("PESAME");
            // c[index]->cd(2);
            // ratio->Draw("PE");

            index++;
        }
    }

    // TCanvas *test = new TCanvas("test", "test", 1000, 1000);
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Sum, 5, 20, kMagenta + 2, 1.0 / h_Pt_DiMuon_Gen_Sum->GetEntries());
    hist1D_graphic_opt(h_Pt_DiMuon_Gen_Sum_Corrected, 5, 20, kRed + 2, 1. / h_Pt_DiMuon_Gen_Sum_Corrected->GetEntries());

    TH1F *ratio = (TH1F *)h_Pt_DiMuon_Gen_Sum_Corrected->Clone("ratio");
    ratio->Divide(h_Pt_DiMuon_Gen_Sum);

    // TF1 *pdf_M_MONASH = new TF1("pdf_M_MONASH", FuncPtMass, 0, 40, 4);
    // pdf_M_MONASH->SetTitle("Fit");
    // pdf_M_MONASH->SetParameter(3, 1);
    // pdf_M_MONASH->SetParameter(0, 3.6);
    // pdf_M_MONASH->SetParameter(1, 2.81);
    // pdf_M_MONASH->SetParameter(2, 2.5);
    // pdf_M_MONASH->SetNpx(300);
    // pdf_M_MONASH->SetLineWidth(3);
    // h_Pt_DiMuon_Gen_DQcut->Fit(pdf_M_MONASH, "LR0I");
    h_Pt_DiMuon_Gen_Sum->SetTitle("Pre corr.");
    h_Pt_DiMuon_Gen_Sum_Corrected->SetTitle("After corr.");
    TCanvas *test=two_histo_ratio(h_Pt_DiMuon_Gen_Sum, h_Pt_DiMuon_Gen_Sum_Corrected, ratio, "test", "#mu#mu #leftarrow c, PYTHIA8 MNR",kTRUE);

    // test->Divide(1, 2);
    // test->cd(1);
    // h_Pt_DiMuon_Gen_Sum->Draw("PE");
    // h_Pt_DiMuon_Gen_Sum_Corrected->Draw("SAME");
    // test->cd(2);
    // ratio->Draw();
}

double FuncPtMass(double *x, double *par)
{
    return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
}
