#include "/home/michele_pennisi/cernbox/common_include.h"
void muon_comparison(Int_t origin=4){

    Double_t Powheg_cs = 0.487e+4;       // Powheg cs in mub
    Double_t Pythia_Monash_cs = 7.85e+4; // Pythia cs in mub
    Double_t Pythia_MNR_cs = 5.42e+4;    // Pythia cs in mub
    Int_t Rebin_X=10;

    TString Origin_powheg;
    TString Origin;
    TString Header;
    if (origin==4){
        Header="Gen. #mu #leftarrow c";
        Origin_powheg="charm";
        Origin="Charm";
    }
    else if (origin==5){
        Header="Gen. #mu #leftarrow b";
        Origin_powheg="beauty";
        Origin="Beauty";
    }
    
    // Getting HF hadron hist from Powheg sim
    TFile *fIn_Powheg = new TFile(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Powheg_Sim/powheg_%s/Version1/save_mc_output/powheg_%s_MC_output_Hist_merged.root",Origin_powheg.Data(),Origin_powheg.Data()), "READ");
    TH2F *h_PtY_Muon_Gen_Charm_POWHEG=(TH2F *)fIn_Powheg->Get(Form("Muon_Gen/h_PtY_Muon_Gen_DQcut_%s",Origin.Data()));
    TH1F *h_Pt_Muon_Gen_Charm_POWHEG=(TH1F *)h_PtY_Muon_Gen_Charm_POWHEG->ProjectionX();
    h_Pt_Muon_Gen_Charm_POWHEG->SetName("h_Pt_Muon_Gen_Charm_POWHEG");
    TH1F *h_Y_Muon_Gen_Charm_POWHEG=(TH1F *)h_PtY_Muon_Gen_Charm_POWHEG->ProjectionY();
    h_Y_Muon_Gen_Charm_POWHEG->SetName("h_Y_Muon_Gen_Charm_POWHEG");

    TH1F *h_Nevents_Powheg = (TH1F *)fIn_Powheg->Get("h_Nevents");
    h_Nevents_Powheg->SetName("h_Nevents_Powheg");
    cout<<h_Nevents_Powheg->GetBinContent(2)<<endl;
    Double_t norm_factor_ALICE_POWHEG = 0.5 * Powheg_cs  / h_Nevents_Powheg->GetBinContent(2);
    h_Pt_Muon_Gen_Charm_POWHEG->Rebin(Rebin_X);
    h_Pt_Muon_Gen_Charm_POWHEG->Scale(1./h_Pt_Muon_Gen_Charm_POWHEG->GetEntries(), "width");
    h_Pt_Muon_Gen_Charm_POWHEG->SetMarkerColor(kCyan + 1);
    h_Pt_Muon_Gen_Charm_POWHEG->SetMarkerStyle(24);
    h_Pt_Muon_Gen_Charm_POWHEG->SetMarkerSize(1.5);
    h_Pt_Muon_Gen_Charm_POWHEG->SetLineColor(kCyan + 1);
    h_Pt_Muon_Gen_Charm_POWHEG->SetLineWidth(2);    

    TFile *fIn_Pythia_Monash = new TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Pythia_Sim/MB/Version3/save_mc_output/MB_MC_output_Hist_merged.root", "READ");

    TH2F *h_PtY_Muon_Gen_Charm_Pythia_Monash=(TH2F *)fIn_Pythia_Monash->Get(Form("Muon_Gen/h_PtY_Muon_Gen_DQcut_%s",Origin.Data()));
    TH1F *h_Pt_Muon_Gen_Charm_Pythia_Monash=(TH1F *)h_PtY_Muon_Gen_Charm_Pythia_Monash->ProjectionX();
    h_Pt_Muon_Gen_Charm_Pythia_Monash->SetName("h_Pt_Muon_Gen_Charm_Pythia_Monash");
    TH1F *h_Y_Muon_Gen_Charm_Pythia_Monash=(TH1F *)h_PtY_Muon_Gen_Charm_Pythia_Monash->ProjectionY();
    h_Y_Muon_Gen_Charm_Pythia_Monash->SetName("h_Y_Muon_Gen_Charm_Pythia_Monash");

    TH1F *h_Nevents_Pythia_Monash = (TH1F *)fIn_Pythia_Monash->Get("h_Nevents");
    h_Nevents_Pythia_Monash->SetName("h_Nevents_Pythia_Monash");
    cout<<h_Nevents_Pythia_Monash->GetBinContent(2)<<endl;
    Double_t norm_factor_ALICE_Pythia_Monash = (Pythia_Monash_cs) / (2. * h_Nevents_Pythia_Monash->GetBinContent(2));
    
    h_Pt_Muon_Gen_Charm_Pythia_Monash->Rebin(Rebin_X);
    h_Pt_Muon_Gen_Charm_Pythia_Monash->Scale(1./h_Pt_Muon_Gen_Charm_Pythia_Monash->GetEntries(), "width");
    h_Pt_Muon_Gen_Charm_Pythia_Monash->SetMarkerColor(kGreen + 2);
    h_Pt_Muon_Gen_Charm_Pythia_Monash->SetMarkerStyle(20);
    h_Pt_Muon_Gen_Charm_Pythia_Monash->SetMarkerSize(1.5);
    h_Pt_Muon_Gen_Charm_Pythia_Monash->SetLineColor(kGreen + 2);
    h_Pt_Muon_Gen_Charm_Pythia_Monash->SetLineWidth(2);

    TFile *fIn_Pythia_MNR = new TFile("test/HF_MC_output_Hist_294009.root", "READ");


    TH2F *h_PtY_Muon_Gen_Charm_Pythia_MNR=(TH2F *)fIn_Pythia_MNR->Get(Form("Muon_Gen/h_PtY_Muon_Gen_DQcut_%s",Origin.Data()));
    TH1F *h_Pt_Muon_Gen_Charm_Pythia_MNR=(TH1F *)h_PtY_Muon_Gen_Charm_Pythia_MNR->ProjectionX();
    h_Pt_Muon_Gen_Charm_Pythia_MNR->SetName("h_Pt_Muon_Gen_Charm_Pythia_MNR");
    TH1F *h_Y_Muon_Gen_Charm_Pythia_MNR=(TH1F *)h_PtY_Muon_Gen_Charm_Pythia_MNR->ProjectionY();
    h_Y_Muon_Gen_Charm_Pythia_MNR->SetName("h_Y_Muon_Gen_Charm_Pythia_MNR");
    
    TH1F *h_Nevents_Pythia_MNR = (TH1F *)fIn_Pythia_MNR->Get("h_Nevents");
    h_Nevents_Pythia_MNR->SetName("h_Nevents_Pythia_MNR");
    Double_t norm_factor_ALICE_Pythia_MNR = (h_Pt_Muon_Gen_Charm_Pythia_Monash->GetEntries()/h_Pt_Muon_Gen_Charm_Pythia_MNR->GetEntries())*(Pythia_Monash_cs) / (2. * h_Nevents_Pythia_Monash->GetBinContent(2));
    
    h_Pt_Muon_Gen_Charm_Pythia_MNR->Rebin(Rebin_X);
    h_Pt_Muon_Gen_Charm_Pythia_MNR->Scale(1./h_Pt_Muon_Gen_Charm_Pythia_MNR->GetEntries(), "width");
    h_Pt_Muon_Gen_Charm_Pythia_MNR->SetMarkerColor(kMagenta + 2);
    h_Pt_Muon_Gen_Charm_Pythia_MNR->SetMarkerStyle(47);
    h_Pt_Muon_Gen_Charm_Pythia_MNR->SetMarkerSize(1.5);
    h_Pt_Muon_Gen_Charm_Pythia_MNR->SetLineColor(kMagenta + 2);
    h_Pt_Muon_Gen_Charm_Pythia_MNR->SetLineWidth(2);

    h_Pt_Muon_Gen_Charm_Pythia_Monash->Draw("PE");
    h_Pt_Muon_Gen_Charm_Pythia_MNR->Draw("PESAME");
    h_Pt_Muon_Gen_Charm_POWHEG->Draw("PESAME");

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas(Form("canvas_mu_%s",Origin.Data()), Form("canvas_mu_%s",Origin.Data()), 1000, 1600);
    c->SetTicks();
    c->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.50, 1.0, 1.00);
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.50);
    
    pad1->Draw();
    pad2->Draw();
    
    pad1->cd();
    pad1->SetTicks();
    pad1->SetLogy(1);
    pad1->SetTopMargin(0.05);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.14);
    pad1->SetBottomMargin(0.0);
    h_Pt_Muon_Gen_Charm_Pythia_Monash->GetYaxis()->SetLabelSize(0.045);
    h_Pt_Muon_Gen_Charm_Pythia_Monash->GetYaxis()->SetTitleSize(0.055);
    h_Pt_Muon_Gen_Charm_Pythia_Monash->GetYaxis()->SetTitle("d^{2}#sigma/d#it{p}_{T}d#it{y} #mubGeV^{-1}#it{c}");
    h_Pt_Muon_Gen_Charm_POWHEG->SetTitle("PYTHIA8 MNR");
    h_Pt_Muon_Gen_Charm_Pythia_Monash->SetTitle("PYTHIA8 Monash");
    h_Pt_Muon_Gen_Charm_Pythia_MNR->SetTitle("POWHEG+PYTHIA6");

    h_Pt_Muon_Gen_Charm_Pythia_Monash->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h_Pt_Muon_Gen_Charm_Pythia_Monash->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    h_Pt_Muon_Gen_Charm_Pythia_Monash->GetXaxis()->SetRangeUser(0,20);
    h_Pt_Muon_Gen_Charm_Pythia_Monash->Draw();
    h_Pt_Muon_Gen_Charm_Pythia_MNR->Draw("SAME");
    h_Pt_Muon_Gen_Charm_POWHEG->Draw("SAME");
    TLegend *legend = new TLegend(0.575, 0.675, 0.825, 0.875);
    legend->SetHeader(Header);
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0425);

    legend->AddEntry(h_Pt_Muon_Gen_Charm_Pythia_Monash);
    legend->AddEntry(h_Pt_Muon_Gen_Charm_Pythia_MNR);
    legend->AddEntry(h_Pt_Muon_Gen_Charm_POWHEG);
    legend->Draw();

    pad2->cd();
    pad2->SetTicks();
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.03);
    pad2->SetLeftMargin(0.14);
    pad2->SetBottomMargin(0.175);
    TH1F *ratio_powheg_Monash = (TH1F *)h_Pt_Muon_Gen_Charm_Pythia_Monash->Clone("ratio_powheg_Monash");
    ratio_powheg_Monash->GetYaxis()->SetTitle(" ");
    ratio_powheg_Monash->GetXaxis()->SetTitleFont(42);
    ratio_powheg_Monash->GetXaxis()->SetLabelFont(42);
    ratio_powheg_Monash->GetXaxis()->SetLabelSize(0.05);
    ratio_powheg_Monash->GetXaxis()->SetTitleSize(0.06);
    ratio_powheg_Monash->GetXaxis()->SetTitleOffset(1.2);
    ratio_powheg_Monash->GetYaxis()->SetLabelFont(42);
    ratio_powheg_Monash->GetYaxis()->SetRangeUser(-0.5, 2.5);
    ratio_powheg_Monash->GetYaxis()->SetNdivisions(505);
    ratio_powheg_Monash->GetYaxis()->SetLabelOffset(0.0175);
    ratio_powheg_Monash->GetYaxis()->SetLabelSize(0.05);

    ratio_powheg_Monash->SetMarkerStyle(h_Pt_Muon_Gen_Charm_POWHEG->GetMarkerStyle());
    ratio_powheg_Monash->SetMarkerColor(h_Pt_Muon_Gen_Charm_POWHEG->GetMarkerColor());
    ratio_powheg_Monash->SetLineColor(h_Pt_Muon_Gen_Charm_POWHEG->GetLineColor());
    ratio_powheg_Monash->Divide(h_Pt_Muon_Gen_Charm_POWHEG);

    TH1F *ratio_MNR_Monash = (TH1F *)h_Pt_Muon_Gen_Charm_Pythia_Monash->Clone("ratio_MNR_Monash");
    ratio_MNR_Monash->GetYaxis()->SetTitle(" ");
    ratio_MNR_Monash->Divide(h_Pt_Muon_Gen_Charm_Pythia_MNR);
    ratio_MNR_Monash->SetMarkerStyle(h_Pt_Muon_Gen_Charm_Pythia_MNR->GetMarkerStyle());
    ratio_MNR_Monash->SetMarkerColor(h_Pt_Muon_Gen_Charm_Pythia_MNR->GetMarkerColor());
    ratio_MNR_Monash->SetLineColor(h_Pt_Muon_Gen_Charm_Pythia_MNR->GetLineColor());

    ratio_powheg_Monash->Draw("PESAME");
    ratio_MNR_Monash->Draw("PESAME");
    
    
    TLegend *legend1 = new TLegend(0.175, 0.2, 0.425, 0.4);
    legend1->SetFillStyle(0);
    legend1->SetLineColor(kWhite);
    legend1->SetBorderSize(0);
    legend1->SetTextSize(0.0425);

    legend1->AddEntry(ratio_powheg_Monash,"POWHEG+PYTHIA6/PYTHIA8 Monash");
    legend1->AddEntry(ratio_MNR_Monash,"PYTHIA8 MNR/PYTHIA8 Monash");
    legend1->Draw();

}