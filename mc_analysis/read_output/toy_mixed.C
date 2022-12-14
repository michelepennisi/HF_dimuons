#if !defined(__CINT__) || defined(__CLING__)
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TGrid.h"
#include "TRandom3.h"
#include "TGraph.h"

#endif
void cutt(Int_t n_ev = 10000, Double_t cut = 0.5);

void result();

void toy_mixed()
{
    /*
    
    cutt(10000, 0.00);

    cutt(10000, 0.05);

    cutt(10000, 0.10);

    cutt(10000, 0.15);

    cutt(10000, 0.20);

    cutt(10000, 0.25);

    cutt(10000, 0.3);

    cutt(10000, 0.35);

    cutt(10000, 0.4);

    cutt(10000, 0.45);

    cutt(10000, 0.5);

    cutt(10000, 0.55);

    cutt(10000, 0.6);

    cutt(10000, 0.65);

    cutt(10000, 0.7);

    cutt(10000, 0.75);

    cutt(10000, 0.8);

    cutt(10000, 0.85);

    cutt(10000, 0.9);

    cutt(10000, 0.95);

    cutt(10000, 1.00);

    */

    result();
}

void result()
{
    printf("entro\n");
    const Int_t n_Fraction = 21;
    Double_t Fractions[n_Fraction] = {1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.00};

    TFile *fOut;
    Double_t def_frac_mub[n_Fraction];
    Double_t def_frac_muc[n_Fraction];

    Double_t def_frac_mub_ccut[n_Fraction];
    Double_t def_frac_muc_ccut[n_Fraction];
    Double_t def_frac_mub_bcut[n_Fraction];
    Double_t def_frac_muc_bcut[n_Fraction];

    Double_t def_frac_DiMu_b[n_Fraction];
    Double_t def_frac_DiMu_c[n_Fraction];
    Double_t def_frac_DiMu_m[n_Fraction];

    Double_t def_frac_DiMu_b_ccut[n_Fraction];
    Double_t def_frac_DiMu_c_ccut[n_Fraction];
    Double_t def_frac_DiMu_m_ccut[n_Fraction];

    Double_t def_frac_DiMu_b_bcut[n_Fraction];
    Double_t def_frac_DiMu_c_bcut[n_Fraction];
    Double_t def_frac_DiMu_m_bcut[n_Fraction];

    for (Int_t i = n_Fraction - 1; i >= 0; i--)
    {
        fOut = new TFile(Form("fOut_frac%0.2f.root", Fractions[i]), "READ");
        fOut->ls();

        TH1D *Muon_Fractions = (TH1D *)fOut->Get("Muon/Muon_Fractions");

        def_frac_mub[i] = (Muon_Fractions->GetBinContent(1) / (Muon_Fractions->GetBinContent(1) + Muon_Fractions->GetBinContent(2)));
        def_frac_muc[i] = (Muon_Fractions->GetBinContent(2) / (Muon_Fractions->GetBinContent(1) + Muon_Fractions->GetBinContent(2)));

        printf("mu b %0.3f || mu c %0.3f \n", def_frac_mub[i], def_frac_muc[i]);

        TH1D *Muon_Fractions_ccut = (TH1D *)fOut->Get("Muon/Muon_Fractions_ccut");

        def_frac_mub_ccut[i] = 100 * (Muon_Fractions_ccut->GetBinContent(1) / (Muon_Fractions_ccut->GetBinContent(1) + Muon_Fractions_ccut->GetBinContent(2)));
        def_frac_muc_ccut[i] = 100 * (Muon_Fractions_ccut->GetBinContent(2) / (Muon_Fractions_ccut->GetBinContent(1) + Muon_Fractions_ccut->GetBinContent(2)));

        printf("mu b ccut %0.3f || mu c ccut %0.3f \n", def_frac_mub_ccut[i], def_frac_muc_ccut[i]);

        TH1D *Muon_Fractions_bcut = (TH1D *)fOut->Get("Muon/Muon_Fractions_bcut");

        def_frac_mub_bcut[i] = 100 * (Muon_Fractions_bcut->GetBinContent(1) / (Muon_Fractions_bcut->GetBinContent(1) + Muon_Fractions_bcut->GetBinContent(2)));
        def_frac_muc_bcut[i] = 100 * (Muon_Fractions_bcut->GetBinContent(2) / (Muon_Fractions_bcut->GetBinContent(1) + Muon_Fractions_bcut->GetBinContent(2)));

        printf("mu b bcut %0.3f || mu c bcut %0.3f \n", def_frac_mub_bcut[i], def_frac_muc_bcut[i]);

        TH1D *DiMuon_Fractions = (TH1D *)fOut->Get("DiMuon/DiMuon_Fractions");

        def_frac_DiMu_b[i] = 100 * (DiMuon_Fractions->GetBinContent(1) / (DiMuon_Fractions->GetBinContent(1) + DiMuon_Fractions->GetBinContent(2) + DiMuon_Fractions->GetBinContent(3)));
        def_frac_DiMu_c[i] = 100 * (DiMuon_Fractions->GetBinContent(2) / (DiMuon_Fractions->GetBinContent(1) + DiMuon_Fractions->GetBinContent(2) + DiMuon_Fractions->GetBinContent(3)));
        def_frac_DiMu_m[i] = 100 * (DiMuon_Fractions->GetBinContent(3) / (DiMuon_Fractions->GetBinContent(1) + DiMuon_Fractions->GetBinContent(2) + DiMuon_Fractions->GetBinContent(3)));

        printf("DiMu b %0.3f || DiMu c %0.3f || DiMu m %0.3f \n", def_frac_DiMu_b[i], def_frac_DiMu_c[i], def_frac_DiMu_m[i]);

        TH1D *DiMuon_Fractions_ccut = (TH1D *)fOut->Get("DiMuon/DiMuon_Fractions_ccut");

        def_frac_DiMu_b_ccut[i] = 100 * (DiMuon_Fractions_ccut->GetBinContent(1) / (DiMuon_Fractions_ccut->GetBinContent(1) + DiMuon_Fractions_ccut->GetBinContent(2) + DiMuon_Fractions_ccut->GetBinContent(3)));
        def_frac_DiMu_c_ccut[i] = 100 * (DiMuon_Fractions_ccut->GetBinContent(2) / (DiMuon_Fractions_ccut->GetBinContent(1) + DiMuon_Fractions_ccut->GetBinContent(2) + DiMuon_Fractions_ccut->GetBinContent(3)));
        def_frac_DiMu_m_ccut[i] = 100 * (DiMuon_Fractions_ccut->GetBinContent(3) / (DiMuon_Fractions_ccut->GetBinContent(1) + DiMuon_Fractions_ccut->GetBinContent(2) + DiMuon_Fractions_ccut->GetBinContent(3)));

        printf("DiMu b %0.3f || DiMu c %0.3f || DiMu m %0.3f \n", def_frac_DiMu_b_ccut[i], def_frac_DiMu_c_ccut[i], def_frac_DiMu_m_ccut[i]);

        TH1D *DiMuon_Fractions_bcut = (TH1D *)fOut->Get("DiMuon/DiMuon_Fractions_bcut");

        def_frac_DiMu_b_bcut[i] = 100 * (DiMuon_Fractions_bcut->GetBinContent(1) / (DiMuon_Fractions_bcut->GetBinContent(1) + DiMuon_Fractions_bcut->GetBinContent(2) + DiMuon_Fractions_bcut->GetBinContent(3)));
        def_frac_DiMu_c_bcut[i] = 100 * (DiMuon_Fractions_bcut->GetBinContent(2) / (DiMuon_Fractions_bcut->GetBinContent(1) + DiMuon_Fractions_bcut->GetBinContent(2) + DiMuon_Fractions_bcut->GetBinContent(3)));
        def_frac_DiMu_m_bcut[i] = 100 * (DiMuon_Fractions_bcut->GetBinContent(3) / (DiMuon_Fractions_bcut->GetBinContent(1) + DiMuon_Fractions_bcut->GetBinContent(2) + DiMuon_Fractions_bcut->GetBinContent(3)));

        printf("DiMu b %0.3f || DiMu c %0.3f || DiMu m %0.3f \n", def_frac_DiMu_b_bcut[i], def_frac_DiMu_c_bcut[i], def_frac_DiMu_m_bcut[i]);
    }

    // DiMu_charm_ccut->Draw();
    TCanvas *canvas = new TCanvas("canvas", "canvas", 850, 650);
    canvas->SetTicks();
    canvas->cd();
    gPad->SetTicks();
    // gPad->SetLogy(1);
    gPad->SetTopMargin(0.10);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.15);

    // TLine* l = new TLine(55, 2.15, 115, 2.15);
    // l->SetLineWidth(4);
    // l->SetLineStyle(2);
    // l->SetLineColor(kMagenta);
    TH2D *h_grid_dimucharm = new TH2D("h_grid_dimucharm", "Dimu from Charm ", n_Fraction, 0.0, 1, 10, 0, 100);
    h_grid_dimucharm->SetStats(0);
    for (Int_t q = n_Fraction - 1; q >= 0; q--)
    {
        // printf("%s\n", Form("Frac %0.2f", Fractions[q]));
        h_grid_dimucharm->GetXaxis()->SetBinLabel(n_Fraction - q, Form("%0.0f", Fractions[q] * 100));
    }

    h_grid_dimucharm->LabelsDeflate();
    h_grid_dimucharm->GetYaxis()->SetTitle("%");
    h_grid_dimucharm->GetXaxis()->SetTitleOffset(1.3);
    h_grid_dimucharm->GetXaxis()->SetTitleSize(0.0475);
    h_grid_dimucharm->GetXaxis()->SetLabelSize(0.045);

    h_grid_dimucharm->GetYaxis()->SetNdivisions(505);
    h_grid_dimucharm->GetYaxis()->SetTitleOffset(1.1);
    h_grid_dimucharm->GetYaxis()->SetTitleSize(0.0475);
    h_grid_dimucharm->GetYaxis()->SetLabelSize(0.045);
    h_grid_dimucharm->Draw();

    TGraph *DiMu_charm_ccut = new TGraph(n_Fraction, Fractions, def_frac_DiMu_c_ccut);
    DiMu_charm_ccut->SetName("DiMu_charm_ccut");
    DiMu_charm_ccut->SetMarkerColor(kViolet);
    DiMu_charm_ccut->SetMarkerStyle(20);
    DiMu_charm_ccut->SetMarkerSize(1.5);

    DiMu_charm_ccut->Draw("PSAME");

    TGraph *DiMu_charm_bcut = new TGraph(n_Fraction, Fractions, def_frac_DiMu_c_bcut);
    DiMu_charm_bcut->SetName("DiMu_charm_bcut");
    DiMu_charm_bcut->SetMarkerColor(kGreen);
    DiMu_charm_bcut->SetMarkerStyle(20);
    DiMu_charm_bcut->SetMarkerSize(1.5);

    DiMu_charm_bcut->Draw("PSAME");

    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 850, 650);
    canvas1->SetTicks();
    canvas1->cd();
    gPad->SetTicks();
    // gPad->SetLogy(1);
    gPad->SetTopMargin(0.10);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.15);

    TH2D *h_grid_dimubeauty = new TH2D("h_grid_dimubeauty", "Dimu from Beauty ", n_Fraction, 0.0, 1, 10, 0, 100);
    h_grid_dimubeauty->SetStats(0);
    for (Int_t q = n_Fraction - 1; q >= 0; q--)
    {
        // printf("%s\n", Form("Frac %0.2f", Fractions[q]));
        h_grid_dimubeauty->GetXaxis()->SetBinLabel(n_Fraction - q, Form("%0.0f", Fractions[q] * 100));
    }

    h_grid_dimubeauty->LabelsDeflate();
    // h_grid_dimubeauty->GetXaxis()->SetTitle("Plane");
    h_grid_dimubeauty->GetYaxis()->SetTitle("%");
    h_grid_dimubeauty->GetXaxis()->SetTitleOffset(1.3);
    h_grid_dimubeauty->GetXaxis()->SetTitleSize(0.0475);
    h_grid_dimubeauty->GetXaxis()->SetLabelSize(0.045);

    h_grid_dimubeauty->GetYaxis()->SetNdivisions(505);
    h_grid_dimubeauty->GetYaxis()->SetTitleOffset(1.1);
    h_grid_dimubeauty->GetYaxis()->SetTitleSize(0.0475);
    h_grid_dimubeauty->GetYaxis()->SetLabelSize(0.045);
    h_grid_dimubeauty->Draw();

    TGraph *DiMu_beauty_ccut = new TGraph(n_Fraction, Fractions, def_frac_DiMu_b_ccut);
    DiMu_beauty_ccut->SetName("DiMu_beauty_ccut");
    DiMu_beauty_ccut->SetMarkerColor(kViolet);
    DiMu_beauty_ccut->SetMarkerStyle(20);
    DiMu_beauty_ccut->SetMarkerSize(1.5);

    DiMu_beauty_ccut->Draw("PSAME");

    TGraph *DiMu_beauty_bcut = new TGraph(n_Fraction, Fractions, def_frac_DiMu_b_bcut);
    DiMu_beauty_bcut->SetName("DiMu_beauty_bcut");
    DiMu_beauty_bcut->SetMarkerColor(kGreen);
    DiMu_beauty_bcut->SetMarkerStyle(20);
    DiMu_beauty_bcut->SetMarkerSize(1.5);

    DiMu_beauty_bcut->Draw("PSAME");

    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 850, 650);
    canvas2->SetTicks();
    canvas2->cd();
    gPad->SetTicks();
    // gPad->SetLogy(1);
    gPad->SetTopMargin(0.10);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.15);

    TH2D *h_grid_dimumixed = new TH2D("h_grid_dimumixed", "Mixed HF Dimuon from TOY ", n_Fraction, -0.025, 1, 10, 0, 7.75);
    h_grid_dimumixed->SetStats(0);
    for (Int_t q = n_Fraction - 1; q >= 0; q--)
    {
        // printf("%s\n", Form("Frac %0.2f", Fractions[q]));
        h_grid_dimumixed->GetXaxis()->SetBinLabel(n_Fraction - q, Form("%0.0f", Fractions[q] * 100));
    }

    h_grid_dimumixed->LabelsDeflate();
    h_grid_dimumixed->GetXaxis()->SetTitle("% Fraction cut");
    h_grid_dimumixed->GetYaxis()->SetTitle("% Mixed Dimuon over total");
    h_grid_dimumixed->GetXaxis()->SetTitleOffset(1.1);
    h_grid_dimumixed->GetXaxis()->SetTitleSize(0.0475);
    h_grid_dimumixed->GetXaxis()->SetLabelSize(0.045);

    h_grid_dimumixed->GetYaxis()->SetNdivisions(505);
    h_grid_dimumixed->GetYaxis()->SetTitleOffset(1.1);
    h_grid_dimumixed->GetYaxis()->SetTitleSize(0.0475);
    h_grid_dimumixed->GetYaxis()->SetLabelSize(0.045);
    h_grid_dimumixed->Draw();

    TGraph *DiMu_mixed_ccut = new TGraph(n_Fraction, Fractions, def_frac_DiMu_m_ccut);
    DiMu_beauty_ccut->SetName("DiMu_beauty_ccut");
    DiMu_mixed_ccut->SetMarkerColor(kViolet);
    DiMu_mixed_ccut->SetMarkerStyle(20);
    DiMu_mixed_ccut->SetMarkerSize(1.5);

    DiMu_mixed_ccut->Draw("PSAME");

    TGraph *DiMu_mixed_bcut = new TGraph(n_Fraction, Fractions, def_frac_DiMu_m_bcut);
    DiMu_beauty_bcut->SetName("DiMu_beauty_bcut");
    DiMu_mixed_bcut->SetMarkerColor(kGreen);
    DiMu_mixed_bcut->SetMarkerStyle(20);
    DiMu_mixed_bcut->SetMarkerSize(1.5);

    DiMu_mixed_bcut->Draw("PSAME");

    TLegend *legend = new TLegend(0.175, 0.575, 0.65, 0.775);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);
    legend->SetTextAlign(12);
    legend->AddEntry(DiMu_mixed_ccut, "Cut on charm", "LP");
    legend->AddEntry(DiMu_mixed_bcut, "Cut on beauty", "LP");
    legend->Draw();

    canvas2->SaveAs("TOY_mixed.pdf");
}

void cutt(Int_t n_ev = 10000, Double_t cut = 0.5)
{
    const Int_t n_Fraction = 12;
    TH1D *Muon_Fractions = new TH1D("Muon_Fractions", "Muon_Fractions", 2, -0.5, 1.5);
    // Muon_Fractions->SetStats(0);
    // Muon_Fractions->SetFillColor(38);
    Muon_Fractions->LabelsDeflate();

    TH1D *Muon_Fractions_ccut = new TH1D("Muon_Fractions_ccut", "Muon_Fractions_ccut", 2, -0.5, 1.5);
    // Muon_Fractions_ccut->SetStats(0);
    // Muon_Fractions_ccut->SetFillColor(38);
    Muon_Fractions_ccut->LabelsDeflate();

    TH1D *Muon_Fractions_bcut = new TH1D("Muon_Fractions_bcut", "Muon_Fractions_bcut", 2, -0.5, 1.5);
    // Muon_Fractions_bcut->SetStats(0);
    // Muon_Fractions_bcut->SetFillColor(38);
    Muon_Fractions_bcut->LabelsDeflate();

    TH1D *DiMuon_Fractions = new TH1D("DiMuon_Fractions", "DiMuon_Fractions", 3, -0.5, 2.5);
    // DiMuon_Fractions->SetStats(0);
    // DiMuon_Fractions->SetFillColor(38);
    DiMuon_Fractions->LabelsDeflate();

    TH1D *Control_DiMuon_Fractions = new TH1D("Control_DiMuon_Fractions", "Control_DiMuon_Fractions", 3, -0.5, 2.5);
    Control_DiMuon_Fractions->SetStats(0);
    Control_DiMuon_Fractions->SetFillColor(38);
    Control_DiMuon_Fractions->LabelsDeflate();

    TH1D *DiMuon_Fractions_ccut = new TH1D("DiMuon_Fractions_ccut", "DiMuon_Fractions_ccut", 3, -0.5, 2.5);
    // DiMuon_Fractions_ccut->SetStats(0);
    // DiMuon_Fractions_ccut->SetFillColor(38);
    DiMuon_Fractions_ccut->LabelsDeflate();

    TH1D *DiMuon_Fractions_bcut = new TH1D("DiMuon_Fractions_bcut", "DiMuon_Fractions_bcut", 3, -0.5, 2.5);
    // DiMuon_Fractions_bcut->SetStats(0);
    // DiMuon_Fractions_bcut->SetFillColor(38);
    DiMuon_Fractions_bcut->LabelsDeflate();

    Double_t tot_mu_charm = 0;
    Double_t tot_mu_beauty = 0;
    Double_t tot_mu = 0;

    Double_t tot_mu_charm_ccut = 0;
    Double_t tot_mu_beauty_ccut = 0;
    Double_t tot_mu_ccut = 0;

    Double_t tot_mu_charm_bcut = 0;
    Double_t tot_mu_beauty_bcut = 0;
    Double_t tot_mu_bcut = 0;

    Double_t tot_Dimu_charm = 0;
    Double_t tot_Dimu_beauty = 0;
    Double_t tot_Dimu_mixed = 0;
    Double_t tot_Dimu = 0;

    Double_t tot_Dimu_charm_ccut = 0;
    Double_t tot_Dimu_beauty_ccut = 0;
    Double_t tot_Dimu_mixed_ccut = 0;
    Double_t tot_Dimu_ccut = 0;

    Double_t tot_Dimu_charm_bcut = 0;
    Double_t tot_Dimu_beauty_bcut = 0;
    Double_t tot_Dimu_mixed_bcut = 0;
    Double_t tot_Dimu_bcut = 0;

    Double_t control_tot_Dimu_charm = 0;
    Double_t control_tot_Dimu_beauty = 0;
    Double_t control_tot_Dimu_mixed = 0;
    Double_t control_tot_Dimu = 0;

    // TH1D *Muon_Fractions_ccut[n_Fraction];
    // TH1D *Muon_Fractions_bcut[n_Fraction];

    // TH1D *DiMuon_Fractions_ccut[n_Fraction];
    // TH1D *DiMuon_Fractions_bcut[n_Fraction];

    Double_t Fractions[n_Fraction] = {1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5};

    Double_t Fraction_mu_beauty[n_Fraction][2] = {0};
    Double_t Fraction_mu_charm[n_Fraction][2] = {0};
    Double_t Total_mu[n_Fraction][2] = {0};

    Double_t Fraction_dimu_beauty[n_Fraction][2] = {0};
    Double_t Fraction_dimu_charm[n_Fraction][2] = {0};
    Double_t Fraction_dimu_mixed[n_Fraction][2] = {0};

    Double_t Total_dimu[n_Fraction][2] = {0};

    for (Int_t i = 0; i < n_ev; i++)
    {
        Double_t NMuons_event_double = gRandom->Rndm() * 20;

        Int_t NMuons_event = (Int_t)NMuons_event_double;

        // printf("Evento %d) NMuons_event %d \n", i, NMuons_event);

        if (i % 500 == 0)
        {
            printf("Evento %d\n", i);
        }
        vector<Int_t> muons;
        vector<Int_t> muons_bcut;
        vector<Int_t> muons_ccut;

        double Dimuon_id = gRandom->Rndm();

        if (Dimuon_id < 0.6051)
        {
            muons.push_back(5);
            muons_ccut.push_back(5);

            if (gRandom->Rndm() < cut)
            {
                muons_bcut.push_back(5);
            }
            muons.push_back(5);
            muons_ccut.push_back(5);

            if (gRandom->Rndm() < cut)
            {
                muons_bcut.push_back(5);
            }
            tot_Dimu_beauty++;
            tot_Dimu++;
        }
        else if (Dimuon_id > 0.6051 && Dimuon_id < 0.963837308)
        {
            muons.push_back(4);
            muons_bcut.push_back(4);
            if (gRandom->Rndm() < cut)
            {
                muons_ccut.push_back(4);
            }
            muons.push_back(4);
            muons_bcut.push_back(4);
            if (gRandom->Rndm() < cut)
            {
                muons_ccut.push_back(4);
            }
            tot_Dimu_charm++;
            tot_Dimu++;
        }
        else
        {
            muons.push_back(5);
            muons_ccut.push_back(5);
            if (gRandom->Rndm() < cut)
            {
                muons_ccut.push_back(4);
            }
            muons.push_back(4);
            muons_bcut.push_back(4);
            if (gRandom->Rndm() < cut)
            {
                muons_bcut.push_back(5);
            }
            tot_Dimu_mixed++;
            tot_Dimu++;
        }

        for (Int_t imu1 = 0; imu1 < muons.size(); imu1++)
        {
            if (muons[imu1] == 5)
            {
                tot_mu_beauty++;
                tot_mu++;
            }
            else if (muons[imu1] == 4)
            {
                tot_mu_charm++;
                tot_mu++;
            }

            for (Int_t imu2 = imu1 + 1; imu2 < muons.size(); imu2++)
            {
                if (muons[imu1] == 5 && muons[imu2] == 5)
                {
                    control_tot_Dimu_beauty++;
                    control_tot_Dimu++;
                }
                else if (muons[imu1] == 4 && muons[imu2] == 4)
                {
                    control_tot_Dimu_charm++;
                    control_tot_Dimu++;
                }
                else if ((muons[imu1] == 5 && muons[imu2] == 4) || (muons[imu1] == 4 && muons[imu2] == 5))
                {
                    control_tot_Dimu_mixed++;
                    control_tot_Dimu++;
                }
            }
        }

        for (Int_t imu1 = 0; imu1 < muons_ccut.size(); imu1++)
        {
            if (muons_ccut[imu1] == 5)
            {
                tot_mu_beauty_ccut++;
                tot_mu_ccut++;
            }
            else if (muons_ccut[imu1] == 4)
            {
                tot_mu_charm_ccut++;
                tot_mu_ccut++;
            }

            for (Int_t imu2 = imu1 + 1; imu2 < muons_ccut.size(); imu2++)
            {
                if (muons_ccut[imu1] == 5 && muons_ccut[imu2] == 5)
                {
                    tot_Dimu_beauty_ccut++;
                    tot_Dimu_ccut++;
                }
                else if (muons_ccut[imu1] == 4 && muons_ccut[imu2] == 4)
                {
                    tot_Dimu_charm_ccut++;
                    tot_Dimu_ccut++;
                }
                else if ((muons_ccut[imu1] == 5 && muons_ccut[imu2] == 4) || (muons_ccut[imu1] == 4 && muons_ccut[imu2] == 5))
                {
                    tot_Dimu_mixed_ccut++;
                    tot_Dimu_ccut++;
                }
            }
        }

        for (Int_t imu1 = 0; imu1 < muons_bcut.size(); imu1++)
        {
            if (muons_bcut[imu1] == 5)
            {
                tot_mu_beauty_bcut++;
                tot_mu_bcut++;
            }
            else if (muons_bcut[imu1] == 4)
            {
                tot_mu_charm_bcut++;
                tot_mu_bcut++;
            }

            for (Int_t imu2 = imu1 + 1; imu2 < muons_bcut.size(); imu2++)
            {
                if (muons_bcut[imu1] == 5 && muons_bcut[imu2] == 5)
                {
                    tot_Dimu_beauty_bcut++;
                    tot_Dimu_bcut++;
                }
                else if (muons_bcut[imu1] == 4 && muons_bcut[imu2] == 4)
                {
                    tot_Dimu_charm_bcut++;
                    tot_Dimu_bcut++;
                }
                else if ((muons_bcut[imu1] == 5 && muons_bcut[imu2] == 4) || (muons_bcut[imu1] == 4 && muons_bcut[imu2] == 5))
                {
                    tot_Dimu_mixed_bcut++;
                    tot_Dimu_bcut++;
                }
            }
        }

        // for (Int_t mu = 0; mu < muons.size(); mu++)
        // {
        //     if (muons[mu] == 4)
        //     {
        //         tot_mu_charm++;
        //         tot_mu++;
        //     }
        //     if (muons[mu] == 5)
        //     {
        //         tot_mu_beauty++;
        //         tot_mu++;
        //     }

        //     for (Int_t frac = 0; frac < n_Fraction; frac++)
        //     {
        //         muons_bcut[frac].push_back(4);
        //         muons_ccut[frac].push_back(5);
        //         if (muons[mu] == 4 && prob_cut < Fractions[frac])
        //         {
        //             muons_ccut[frac].push_back(4);
        //         }
        //         if (muons[mu] == 5 && prob_cut < Fractions[frac])
        //         {
        //             muons_bcut[frac].push_back(5);
        //         }
        //     }
        // }

        // for (Int_t frac = 0; frac < n_Fraction; frac++)
        // {
        //     for (Int_t mu_charm = 0; mu_charm < muons_ccut[frac].size(); mu_charm++)
        //     {
        //         if (muons_ccut[frac][mu_charm] == 5)
        //             Fraction_mu_beauty[frac][0]++;
        //         if (muons_ccut[frac][mu_charm] == 4)
        //             Fraction_mu_charm[frac][0]++;
        //     }
        //     for (Int_t mu_beauty = 0; mu_beauty < muons_bcut[frac].size(); mu_beauty++)
        //     {
        //         if (muons_bcut[frac][mu_beauty] == 5)
        //             Fraction_mu_beauty[frac][1]++;
        //         if (muons_bcut[frac][mu_beauty] == 4)
        //             Fraction_mu_charm[frac][1]++;
        //     }
        //     Total_mu[frac][0] = Total_mu[frac][0] + muons_ccut[frac].size();
        //     Total_mu[frac][1] = Total_mu[frac][1] + muons_bcut[frac].size();
        // }
    }
    TFile fOut(Form("fOut_frac%0.2f.root", cut), "recreate");
    fOut.cd();

    TString Mu_label[2];
    Mu_label[0].Form("#mu from b");
    Mu_label[1].Form("#mu from c");

    TString DiMu_label[2];
    DiMu_label[0].Form("#mu#mu from b");
    DiMu_label[1].Form("#mu#mu from c");
    DiMu_label[2].Form("#mu#mu from c,b");
    fOut.mkdir("DiMuon");
    fOut.cd("DiMuon");

    DiMuon_Fractions->Fill(DiMu_label[0], tot_Dimu_beauty);
    DiMuon_Fractions->Fill(DiMu_label[1], tot_Dimu_charm);
    DiMuon_Fractions->Fill(DiMu_label[2], tot_Dimu_mixed);
    DiMuon_Fractions->Write(0, 2, 0);

    Control_DiMuon_Fractions->Fill(DiMu_label[0], control_tot_Dimu_beauty / control_tot_Dimu);
    Control_DiMuon_Fractions->Fill(DiMu_label[1], control_tot_Dimu_charm / control_tot_Dimu);
    Control_DiMuon_Fractions->Fill(DiMu_label[2], control_tot_Dimu_mixed / control_tot_Dimu);
    Control_DiMuon_Fractions->Write(0, 2, 0);

    DiMuon_Fractions_ccut->Fill(DiMu_label[0], tot_Dimu_beauty_ccut);
    DiMuon_Fractions_ccut->Fill(DiMu_label[1], tot_Dimu_charm_ccut);
    DiMuon_Fractions_ccut->Fill(DiMu_label[2], tot_Dimu_mixed_ccut);
    DiMuon_Fractions_ccut->Write(0, 2, 0);

    DiMuon_Fractions_bcut->Fill(DiMu_label[0], tot_Dimu_beauty_bcut);
    DiMuon_Fractions_bcut->Fill(DiMu_label[1], tot_Dimu_charm_bcut);
    DiMuon_Fractions_bcut->Fill(DiMu_label[2], tot_Dimu_mixed_bcut);
    DiMuon_Fractions_bcut->Write(0, 2, 0);
    fOut.cd();
    fOut.mkdir("Muon");
    fOut.cd("Muon");
    Muon_Fractions->Fill(Mu_label[0], tot_mu_beauty);
    Muon_Fractions->Fill(Mu_label[1], tot_mu_charm);
    Muon_Fractions->Write(0, 2, 0);

    Muon_Fractions_ccut->Fill(Mu_label[0], tot_mu_beauty_ccut);
    Muon_Fractions_ccut->Fill(Mu_label[1], tot_mu_charm_ccut);
    Muon_Fractions_ccut->Write(0, 2, 0);

    Muon_Fractions_bcut->Fill(Mu_label[0], tot_mu_beauty_bcut);
    Muon_Fractions_bcut->Fill(Mu_label[1], tot_mu_charm_bcut);
    Muon_Fractions_bcut->Write(0, 2, 0);

    // for (Int_t frac = 0; frac < n_Fraction; frac++)
    // {
    //     Muon_Fractions_ccut[frac] = new TH1D(Form("Muon_Fractions_ccut_Frac%0.2f", Fractions[frac]), Form("Muon_Fractions_ccut_Frac%0.2f", Fractions[frac]), 2, -0.5, 1.5);
    //     Muon_Fractions_ccut[frac]->SetStats(0);
    //     Muon_Fractions_ccut[frac]->SetFillColor(38);
    //     Muon_Fractions_ccut[frac]->LabelsDeflate();

    //     Muon_Fractions_bcut[frac] = new TH1D(Form("Muon_Fractions_bcut_Frac%0.2f", Fractions[frac]), Form("Muon_Fractions_bcut_Frac%0.2f", Fractions[frac]), 2, -0.5, 1.5);
    //     Muon_Fractions_bcut[frac]->SetStats(0);
    //     Muon_Fractions_bcut[frac]->SetFillColor(38);
    //     Muon_Fractions_bcut[frac]->LabelsDeflate();

    //     Muon_Fractions_ccut[frac]->Fill(Mu_label[0], Fraction_mu_beauty[frac][0]);
    //     Muon_Fractions_ccut[frac]->Fill(Mu_label[1], Fraction_mu_charm[frac][0]);

    //     Muon_Fractions_bcut[frac]->Fill(Mu_label[0], Fraction_mu_beauty[frac][1]);
    //     Muon_Fractions_bcut[frac]->Fill(Mu_label[1], Fraction_mu_charm[frac][1]);

    //     // Muon_Fractions_ccut[frac]->SetBinContent(1, Fraction_mu_beauty[frac][0] /Total_mu[frac][0]);
    //     // Muon_Fractions_ccut[frac]->SetBinContent(2, Fraction_mu_charm[frac][0] /Total_mu[frac][0]);

    //     // Muon_Fractions_bcut[frac]->SetBinContent(1, Fraction_mu_beauty[frac][1] / Total_mu[frac][1]);
    //     // Muon_Fractions_bcut[frac]->SetBinContent(2, Fraction_mu_charm[frac][1] / Total_mu[frac][1]);

    //     // Control_DiMuon_ID_ccut[frac]->SetBinContent(1, control_dimu_beauty / (control_dimuons));
    //     // Control_DiMuon_ID_ccut[frac]->SetBinContent(2, total_dimu_charm_ccut[frac] / (control_dimuons));
    //     // Control_DiMuon_ID_ccut[frac]->SetBinContent(3, total_dimu_mixed_ccut[frac] / (control_dimuons));

    //     Muon_Fractions_ccut[frac]->Write(0, 2, 0);
    //     Muon_Fractions_bcut[frac]->Write(0, 2, 0);
    // }

    fOut.Close();
}