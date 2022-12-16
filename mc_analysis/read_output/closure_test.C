void closure_test() {
    const double maxPt = 30;
    const double maxM = 30;

    TFile *fIn = new TFile("~/cernbox/data_HF_dimuons/mc_LHC18p/HistLite_HF_MCDimuHFTree_merged.root", "READ");
    TH2D *h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm = (TH2D*) fIn -> Get("DiMuon/M4/Gen_DQ_cut/ULS/h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm");
    TH2D *h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm = (TH2D*) fIn -> Get("DiMuon/M4/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_fromCharm");

    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm -> RebinX(4);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm -> RebinY(4);

    h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm -> RebinX(4);
    h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm -> RebinY(4);

    const int nBinsPt = h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm -> GetNbinsX();
    const int nBinsM = h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm -> GetNbinsY();

    TH2D *h_PtMDiMu_M4_AxE_DQ_cut_ULS_fromCharm = new TH2D("h_PtMDiMu_M4_AxE_DQ_cut_ULS_fromCharm", "", nBinsPt, 0, 30, nBinsM, 4, 30);
    h_PtMDiMu_M4_AxE_DQ_cut_ULS_fromCharm -> Divide(h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm, h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm, 1, 1, "B");

    std::cout << nBinsPt << " " << nBinsM << std::endl;

    TCanvas *canvasComp_AxE = new TCanvas("canvasComp_AxE", "", 1800, 600);
    canvasComp_AxE -> Divide(3, 1);
    canvasComp_AxE -> cd(1);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm -> Draw("COLZ");

    canvasComp_AxE -> cd(2);
    h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm -> Draw("COLZ");

    canvasComp_AxE -> cd(3);
    h_PtMDiMu_M4_AxE_DQ_cut_ULS_fromCharm -> Draw("COLZ");
 

    TH1D *h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt = h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm -> ProjectionX("h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt");
    TH1D *h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjPt = h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm -> ProjectionX("h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjPt");

    TH1D *h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM = h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm -> ProjectionY("h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM");
    TH1D *h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjM = h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm -> ProjectionY("h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjM");

    TH1D *h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt_AxE = new TH1D("h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt_AxE", "", nBinsPt, 0, maxPt);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt_AxE -> Divide(h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjPt, h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt, 1, 1, "B");
    
    TH1D *h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM_AxE = new TH1D("h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM_AxE", "", nBinsM, 4, maxM);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM_AxE -> Divide(h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjM, h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM, 1, 1, "B");

    TLine *lineUnityPt = new TLine(0, 1., 30., 1.);
    lineUnityPt -> SetLineColor(kGray+1);
    lineUnityPt -> SetLineWidth(2);
    lineUnityPt -> SetLineStyle(kDashed);

    TLine *lineUnityM = new TLine(4, 1., 30., 1.);
    lineUnityM -> SetLineColor(kGray+1);
    lineUnityM -> SetLineWidth(2);
    lineUnityM -> SetLineStyle(kDashed);


    TCanvas *canvasCompProj_AxE = new TCanvas("canvasCompProj_AxE", "", 1200, 1200);
    canvasCompProj_AxE -> Divide(2, 2);
    
    canvasCompProj_AxE -> cd(1);
    gPad -> SetLogy(1);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt -> SetLineColor(kAzure+2);
    h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjPt -> SetLineColor(kRed+1);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt -> Draw();
    h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjPt -> Draw("same");

    canvasCompProj_AxE -> cd(2);
    gPad -> SetLogy(1);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM -> SetLineColor(kAzure+2);
    h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjM -> SetLineColor(kRed+1);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM -> Draw();
    h_PtMDiMu_M4_Rec_DQ_cut_ULS_fromCharm_ProjM -> Draw("same");


    canvasCompProj_AxE -> cd(3);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt_AxE -> GetYaxis() -> SetRangeUser(0, 1.2);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt_AxE -> SetLineColor(kBlack);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjPt_AxE -> Draw();
    lineUnityPt -> Draw("same");

    canvasCompProj_AxE -> cd(4);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM_AxE -> GetYaxis() -> SetRangeUser(0, 1.2);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM_AxE -> SetLineColor(kBlack);
    h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromCharm_ProjM_AxE -> Draw();
    lineUnityM -> Draw("same");

    TFile *fOut = new TFile("AxE_LUT.root", "RECREATE");
    h_PtMDiMu_M4_AxE_DQ_cut_ULS_fromCharm -> Write();
    fOut -> Close();
}