void dimuon_acc()
{
    const Int_t n_scan = 7;
    Double_t mass_cut[n_scan] = {4, 7, 10.0, 12.0, 15.0, 20.0, 30.0};

    Double_t pt_cut[n_scan] = {0.0, 5, 10.0, 12.0, 15.0, 20.0, 30.0};

    TFile *fIn = new TFile("/home/michele_pennisi/cernbox/data_HF_dimuons/mc_LHC18p/HistLite_HF_MCDimuHFTree_merged.root", "READ");

    TH2D *h_PtM_Gen = (TH2D *)fIn->Get("DiMuon/M4/Gen_DQ_cut/ULS/h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromHF");
    h_PtM_Gen->Draw("text");
    h_PtM_Gen->RebinX(30);
    h_PtM_Gen->RebinY(26);

    h_PtM_Gen->Integral(0, 10, 0, 10);

    printf("GEN)Integral pt<10 and m<10: %0.1f \n", h_PtM_Gen->Integral(0, 10, 0, 10));

    TH2D *h_PtM_Rec = (TH2D *)fIn->Get("DiMuon/M4/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromHF");
    // h_PtM_Rec->RebinX(30);
    // h_PtM_Rec->RebinY(26);
    // TH1D *h_Pt_Rec = (TH1D *)h_PtM_Rec->ProjectionX();
    // TH1D *h_M_Rec = (TH1D *)h_PtM_Rec->ProjectionY();
    TCanvas *c1[n_scan];
    TH2D *h_PtM_Rec_resize[n_scan];
    TH1D *h_Pt_Rec_resize[n_scan];
    TH1D *h_M_Rec_resize[n_scan];
    for (Int_t i = 0; i < n_scan - 1; i++)
    {
        c1[i] = new TCanvas(Form("canvas_opt%d", i), Form("canvas %0.0f < #it{m} < %0.0f & %0.0f < #it{p}_{T} < %0.0f", mass_cut[i], mass_cut[i + 1], pt_cut[i], pt_cut[i + 1]), 1600, 1600);
        c1[i]->Divide(3, 1);
        c1[i]->cd(1);
        h_PtM_Rec_resize[i] = (TH2D *)h_PtM_Rec->Clone(Form("h_PtM_Rec_resize_opt%d", i));
        h_PtM_Rec_resize[i]->SetTitle(Form("Rec %0.0f < #it{m} < %0.0f & %0.0f < #it{p}_{T} < %0.0f", mass_cut[i], mass_cut[i + 1], pt_cut[i], pt_cut[i + 1]));

        h_PtM_Rec_resize[i]->GetYaxis()->SetRangeUser(mass_cut[i], mass_cut[i + 1]);
        h_PtM_Rec_resize[i]->GetXaxis()->SetRangeUser(pt_cut[i], pt_cut[i + 1]);
        h_PtM_Rec_resize[i]->Draw("COLZ");

        printf("REC)Integral %0.0f < #it{m} < %0.0f & %0.0f < #it{p}_{T} < %0.0f: %0.0f (%0.1e %%) \n", mass_cut[i], mass_cut[i + 1], pt_cut[i], pt_cut[i + 1], h_PtM_Rec_resize[i]->Integral(), h_PtM_Rec_resize[i]->Integral()/ h_PtM_Rec->Integral());

        c1[i]->cd(2);
        h_Pt_Rec_resize[i] = (TH1D *)h_PtM_Rec_resize[i]->ProjectionX();
        h_Pt_Rec_resize[i]->Draw();
        c1[i]->cd(3);
        h_M_Rec_resize[i] = (TH1D *)h_PtM_Rec_resize[i]->ProjectionY();

        h_M_Rec_resize[i]->Draw();

        c1[i]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/plot/%s.png", c1[i]->GetName()));
    }
}