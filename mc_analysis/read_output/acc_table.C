#include "/home/michele_pennisi/cernbox/common_include.h"

TFile fIn("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Version1/HF/Analysis_MCsim/HF_Analysis_MCsim_merged.root", "READ");

TFile fOut("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/Version1/HF/Analysis_MCsim/HF_Analysis_MCsim_Acc_Table_merged.root", "RECREATE");
void acc_table_muon()
{
    TString Origin[2] = {"Charm", "Beauty"};
    TString Type[2] = {"Gen", "Rec"};

    TH3F *h_PtYPdg_Muon_Charm_Gen = (TH3F *)fIn.Get(Form("Muon_%s/h_PtYPdg_Muon_%s%s", Type[0].Data(), Type[0].Data(), Origin[0].Data()));

    TH2F *h_PtY_Muon_Charm_Gen = (TH2F *)h_PtYPdg_Muon_Charm_Gen->Project3D("yxe");
    h_PtY_Muon_Charm_Gen->RebinX(15);
    h_PtY_Muon_Charm_Gen->RebinY(10);

    TH3F *h_PtYPdg_Muon_Charm_Rec = (TH3F *)fIn.Get(Form("Muon_%s/h_PtYPdg_Muon_%s%s", Type[1].Data(), Type[1].Data(), Origin[0].Data()));

    TH2F *h_PtY_Muon_Charm_Rec = (TH2F *)h_PtYPdg_Muon_Charm_Rec->Project3D("yxe");
    h_PtY_Muon_Charm_Rec->RebinX(15);
    h_PtY_Muon_Charm_Rec->RebinY(10);

    TH2F *h_PtY_Muon_Charm_Acc = (TH2F *)h_PtY_Muon_Charm_Rec->Clone("h_PtY_Muon_Charm_Acc");
    h_PtY_Muon_Charm_Acc->Divide(h_PtY_Muon_Charm_Gen);
    TCanvas *charm = new TCanvas("canvas_charm_muon", "canvas_charm_muon", 1200, 1200);
    charm->cd();
    Bool_t good_norm=kFALSE;
    Int_t Final=30;
    while (!good_norm)
    {
        if (h_PtY_Muon_Charm_Acc->GetMaximum()>1.0)
        {
            Final=Final-1;
            h_PtY_Muon_Charm_Acc->GetXaxis()->SetRangeUser(0,Final);
        }else 
            good_norm=kTRUE;
        
    }
    h_PtY_Muon_Charm_Acc->GetXaxis()->SetRangeUser(0,Final);
    h_PtY_Muon_Charm_Acc->Draw("COLZ");

    TH3F *h_PtYPdg_Muon_Beauty_Gen = (TH3F *)fIn.Get(Form("Muon_%s/h_PtYPdg_Muon_%s%s", Type[0].Data(), Type[0].Data(), Origin[1].Data()));

    TH2F *h_PtY_Muon_Beauty_Gen = (TH2F *)h_PtYPdg_Muon_Beauty_Gen->Project3D("yxe");
    h_PtY_Muon_Beauty_Gen->RebinX(15);
    h_PtY_Muon_Beauty_Gen->RebinY(10);

    TH3F *h_PtYPdg_Muon_Beauty_Rec = (TH3F *)fIn.Get(Form("Muon_%s/h_PtYPdg_Muon_%s%s", Type[1].Data(), Type[1].Data(), Origin[1].Data()));

    TH2F *h_PtY_Muon_Beauty_Rec = (TH2F *)h_PtYPdg_Muon_Beauty_Rec->Project3D("yxe");
    h_PtY_Muon_Beauty_Rec->RebinX(15);
    h_PtY_Muon_Beauty_Rec->RebinY(10);

    TH2F *h_PtY_Muon_Beauty_Acc = (TH2F *)h_PtY_Muon_Beauty_Rec->Clone("h_PtY_Muon_Beauty_Acc");
    h_PtY_Muon_Beauty_Acc->Divide(h_PtY_Muon_Beauty_Gen);
    TCanvas *beauty = new TCanvas("canvas_beauty_muon", "canvas_beauty_muon", 1200, 1200);
    beauty->cd();
    good_norm=kFALSE;
    Final=30;
    while (!good_norm)
    {
        if (h_PtY_Muon_Beauty_Acc->GetMaximum()>1.0)
        {
            Final=Final-1;
            h_PtY_Muon_Beauty_Acc->GetXaxis()->SetRangeUser(0,Final);
        }else 
            good_norm=kTRUE;
        
    }
    h_PtY_Muon_Beauty_Acc->GetXaxis()->SetRangeUser(0,Final);
    fOut.cd();
    if (!fOut.GetDirectory("Muon_Corr"))
        fOut.mkdir("Muon_Corr");
    fOut.cd("Muon_Corr");
    h_PtY_Muon_Charm_Acc->Write();
    h_PtY_Muon_Beauty_Acc->Write();
    
    TH1D *prova=(TH1D *)h_PtY_Muon_Beauty_Acc->ProjectionX();
    prova->Scale(1./15);
    prova->Draw("PE");
}

void acc_table_dimuon()
{
    TString Origin[2] = {"Charm", "Beauty"};
    TString Type[2] = {"Gen", "Rec"};

    fIn.cd("DiMu_Gen");
    fIn.ls();

    TH3F *h_PtYPdg_DiMu_Gen_ULS_M49_Pt010_Charm = (TH3F *)fIn.Get(Form("DiMu_%s/h_PtYPdg_DiMu_%s_ULS_M49_Pt010_%s", Type[0].Data(), Type[0].Data(), Origin[0].Data()));
    TH2F *h_PtY_DiMu_Gen_ULS_M49_Pt010_Charm = (TH2F *)h_PtYPdg_DiMu_Gen_ULS_M49_Pt010_Charm->Project3D("yxe");
    h_PtY_DiMu_Gen_ULS_M49_Pt010_Charm->RebinX(10);
    h_PtY_DiMu_Gen_ULS_M49_Pt010_Charm->RebinY(10);

    TH3F *h_PtYPdg_DiMu_Rec_ULS_M49_Pt010_Charm = (TH3F *)fIn.Get(Form("DiMu_%s/h_PtYPdg_DiMu_%s_ULS_M49_Pt010_%s", Type[1].Data(), Type[1].Data(), Origin[0].Data()));
    TH2F *h_PtY_DiMu_Rec_ULS_M49_Pt010_Charm = (TH2F *)h_PtYPdg_DiMu_Rec_ULS_M49_Pt010_Charm->Project3D("yxe");
    h_PtY_DiMu_Rec_ULS_M49_Pt010_Charm->RebinX(10);
    h_PtY_DiMu_Rec_ULS_M49_Pt010_Charm->RebinY(10);

    TH2F *h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm = (TH2F *)h_PtY_DiMu_Rec_ULS_M49_Pt010_Charm->Clone("h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm");

    h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm->Divide(h_PtY_DiMu_Gen_ULS_M49_Pt010_Charm);

    TCanvas *charm = new TCanvas("canvas_charm_DiMu", "canvas_charm_DiMu", 1200, 1200);
    charm->cd();
    h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm->Draw("COLZ");

    TH3F *h_PtYPdg_DiMu_Gen_ULS_M49_Pt010_Beauty = (TH3F *)fIn.Get(Form("DiMu_%s/h_PtYPdg_DiMu_%s_ULS_M49_Pt010_%s", Type[0].Data(), Type[0].Data(), Origin[1].Data()));
    TH2F *h_PtY_DiMu_Gen_ULS_M49_Pt010_Beauty = (TH2F *)h_PtYPdg_DiMu_Gen_ULS_M49_Pt010_Beauty->Project3D("yxe");
    h_PtY_DiMu_Gen_ULS_M49_Pt010_Beauty->RebinX(10);
    h_PtY_DiMu_Gen_ULS_M49_Pt010_Beauty->RebinY(10);

    TH3F *h_PtYPdg_DiMu_Rec_ULS_M49_Pt010_Beauty = (TH3F *)fIn.Get(Form("DiMu_%s/h_PtYPdg_DiMu_%s_ULS_M49_Pt010_%s", Type[1].Data(), Type[1].Data(), Origin[1].Data()));
    TH2F *h_PtY_DiMu_Rec_ULS_M49_Pt010_Beauty = (TH2F *)h_PtYPdg_DiMu_Rec_ULS_M49_Pt010_Beauty->Project3D("yxe");
    h_PtY_DiMu_Rec_ULS_M49_Pt010_Beauty->RebinX(10);
    h_PtY_DiMu_Rec_ULS_M49_Pt010_Beauty->RebinY(10);

    TH2F *h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty = (TH2F *)h_PtY_DiMu_Rec_ULS_M49_Pt010_Beauty->Clone("h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty");

    h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty->Divide(h_PtY_DiMu_Gen_ULS_M49_Pt010_Beauty);
    TCanvas *beauty = new TCanvas("canvas_beauty_DiMu", "canvas_beauty_DiMu", 1200, 1200);
    beauty->cd();
    h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty->Draw("COLZ");

    fOut.cd();
    if (!fOut.GetDirectory("DiMu_Corr"))
        fOut.mkdir("DiMu_Corr");
    fOut.cd("DiMu_Corr");
    h_PtY_DiMu_Acc_ULS_M49_Pt010_Charm->Write();
    h_PtY_DiMu_Acc_ULS_M49_Pt010_Beauty->Write();
}

void test()
{
    TString Origin[2] = {"Charm", "Beauty"};
    TString Type[2] = {"Meson", "Barion"};

    TH3F *h_PtYPdg_DiMu_Rec_ULS_M4cut = (TH3F *)fIn.Get(Form("DiMu_Rec/h_PtYPdg_DiMu_Rec_%s_ULS_M4cut", Type[0].Data()));
    TH3F *h_PtYPdg_DiMu_Rec_ULS = (TH3F *)fIn.Get(Form("DiMu_Rec/h_PtYPdg_DiMu_Rec_%s_ULS", Type[0].Data()));

    TH3F *h_PtYPdg_DiMu_Rec_ULS_Cloned_M4cut = (TH3F *)h_PtYPdg_DiMu_Rec_ULS_M4cut->Clone(Form("h_PtYPdg_DiMu_Rec_%s_%s_ULS_M4cut", Origin[0].Data(), Type[0].Data()));
    h_PtYPdg_DiMu_Rec_ULS_Cloned_M4cut->GetZaxis()->SetRange(1, 1);

    TH3F *h_PtYPdg_DiMu_Rec_ULS_Cloned = (TH3F *)h_PtYPdg_DiMu_Rec_ULS->Clone(Form("h_PtYPdg_DiMu_Rec_%s_%s_ULS", Origin[0].Data(), Type[0].Data()));
    h_PtYPdg_DiMu_Rec_ULS_Cloned->GetZaxis()->SetRange(1, 1);

    TH2F *h_PtY_DiMu_Rec_ULS_M4cut = (TH2F *)h_PtYPdg_DiMu_Rec_ULS_Cloned_M4cut->Project3D("xye");
    h_PtY_DiMu_Rec_ULS_M4cut->RebinX(15);
    h_PtY_DiMu_Rec_ULS_M4cut->RebinY(15);

    TH2F *h_PtY_DiMu_Rec_ULS = (TH2F *)h_PtYPdg_DiMu_Rec_ULS_Cloned->Project3D("xye");
    h_PtY_DiMu_Rec_ULS->RebinX(15);
    h_PtY_DiMu_Rec_ULS->RebinY(15);

    TH2F *ratio = (TH2F *)h_PtY_DiMu_Rec_ULS_M4cut->Clone("ratio");
    ratio->Divide(h_PtY_DiMu_Rec_ULS);
    ratio->Draw("COLZ");
}

void acc_table()
{
    acc_table_muon();
    acc_table_dimuon();
}