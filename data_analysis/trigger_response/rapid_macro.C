void rapid_macro()
{
    TFile *fIn = new TFile("HF_TriggeResponse_294009.root", "READ");

    TList *list = (TList *)fIn->Get("MCList");

    TH1D *h_allpt = (TH1D *)list->FindObject("fHistPt_allmuon_alleta");
    h_allpt->Rebin(5);
    

    TH1D *h_allpt_lowpt_thr = (TH1D *)list->FindObject("fHistPt_allmuon_alleta_lowpt_thr");
    h_allpt_lowpt_thr->Rebin(5);
    

    TH1D *trigger_response=(TH1D*)h_allpt_lowpt_thr->Clone("trigger_response");
    trigger_response->Divide(h_allpt);
    
    trigger_response->SetMarkerStyle(20);
    trigger_response->Draw("P");
}