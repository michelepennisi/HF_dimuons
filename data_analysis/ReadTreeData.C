#include "/home/michele_pennisi/cernbox/common_include.h"
#include "TKey.h"
void ReadTreeData(TString Run_Period = "LHC18p")
{

    TFile *fIn = new TFile(Form("/home/michele_pennisi/Remote_FarmData/Run2_Data/%s/Tree_MassCut4_%s.root", Run_Period.Data(), Run_Period.Data()), "READ");
    TIter next(fIn->GetListOfKeys());
    TKey *key = new TKey();
    TObject *iobj = nullptr;
    while ((key = (TKey *)next()))
    {

        if (TString::Format("%s", key->GetClassName()).Contains("TTree"))
        {
            iobj = (TObject *)key->ReadObj();
            iobj->InheritsFrom(TString::Format("TTree::Class()"));
            ((TTree *)iobj)->SetDirectory(gROOT);
        }
    }

    bool IsPhysSelected;
    char Trigger[500];
    Int_t NDimu, NMuons;
    Int_t pDCA[20], DimuCharge[20], DimuMatch[20];
    Double_t DimuPt[20], DimuMass[20], DimuY[20], RAtAbsEnd[20], Eta[20];
    Int_t DimuMu[20][2];
    ((TTree *)iobj)->SetBranchAddress("FiredTriggerClasses", Trigger);
    ((TTree *)iobj)->SetBranchAddress("NDimu", &NDimu);
    ((TTree *)iobj)->SetBranchAddress("Eta", Eta);
    ((TTree *)iobj)->SetBranchAddress("RAtAbsEnd", RAtAbsEnd);
    ((TTree *)iobj)->SetBranchAddress("pDCA", pDCA);
    ((TTree *)iobj)->SetBranchAddress("DimuY", DimuY);
    ((TTree *)iobj)->SetBranchAddress("DimuMass", DimuMass);
    ((TTree *)iobj)->SetBranchAddress("DimuPt", DimuPt);
    ((TTree *)iobj)->SetBranchAddress("DimuCharge", DimuCharge);
    ((TTree *)iobj)->SetBranchAddress("DimuMatch", DimuMatch);
    ((TTree *)iobj)->SetBranchAddress("DimuMu", DimuMu);
    ((TTree *)iobj)->SetBranchAddress("IsPhysSelected", &IsPhysSelected);

    Int_t Nentries = ((TTree *)iobj)->GetEntries();
    Int_t DimuGood = 0;
    Double_t *Pt_Dimu_Rec = new Double_t;
    Double_t *M_Dimu_Rec = new Double_t;
    TFile *fout = new TFile(Form("/home/michele_pennisi/Remote_FarmData/Run2_Data/%s/Tree_MassPtCut4_%s.root", Run_Period.Data(), Run_Period.Data()), "RECREATE");

    TH1D *h_Nev = new TH1D("h_Nev", ";ev", 1, 0, 1);
    h_Nev->SetBinContent(1, Nentries);
    TTree *Tree_DiMuon_Rec = new TTree("rec_data_tree", "Data dimuons with 4 < m <30 Gev");
    Tree_DiMuon_Rec->Branch("m", M_Dimu_Rec, "m/D");
    Tree_DiMuon_Rec->Branch("pt", Pt_Dimu_Rec, "pt/D");
    for (Int_t entry = 0; entry < Nentries; entry++)
    {
        if (entry % (Int_t)(Nentries * 0.1) == 0)
            progress_status(entry, Nentries);

        ((TTree *)iobj)->GetEntry(entry);
        TString Ev_Trigger = Trigger;
        if (!(Ev_Trigger.Contains("CMUL7-B-NOPF-MUFAST") && IsPhysSelected))
            continue;

        for (Int_t i_DiMu = 0; i_DiMu < NDimu; i_DiMu++)
        {
            if (DimuMatch[i_DiMu] != 2)
                continue;
            Double_t Y_DiMu = DimuY[i_DiMu];
            if (Y_DiMu < -4 || Y_DiMu > -2.5)
                continue;
            Double_t Pt_DiMu = DimuPt[i_DiMu];
            Double_t Mass_DiMu = DimuMass[i_DiMu];
            if (Mass_DiMu > 30)
                continue;

            Int_t Charge_DiMu = DimuCharge[i_DiMu];
            if (Charge_DiMu != 0)
                continue;

            Double_t Mu0_pDCA = pDCA[DimuMu[i_DiMu][0]];
            Double_t Mu1_pDCA = pDCA[DimuMu[i_DiMu][1]];
            if (Mu0_pDCA != 1 || Mu1_pDCA != 1)
                continue;

            Double_t Mu0_RAbs = RAtAbsEnd[DimuMu[i_DiMu][0]];
            Double_t Mu1_RAbs = RAtAbsEnd[DimuMu[i_DiMu][1]];
            if ((Mu0_RAbs < 17.6 || Mu0_RAbs > 89.5) || (Mu1_RAbs < 17.6 || Mu1_RAbs > 89.5))
                continue;
            Double_t Mu0_Eta = Eta[DimuMu[i_DiMu][0]];
            Double_t Mu1_Eta = Eta[DimuMu[i_DiMu][1]];
            if ((Mu0_Eta < -4.0 || Mu0_Eta > -2.5) || (Mu1_Eta < -4.0 || Mu1_Eta > -2.5))
                continue;
            DimuGood++;
            *M_Dimu_Rec = Mass_DiMu;
            *Pt_Dimu_Rec = Pt_DiMu;
            Tree_DiMuon_Rec->Fill();
        }
    }

    h_Nev->Write();
    Tree_DiMuon_Rec->Write();

    delete iobj;
    delete Tree_DiMuon_Rec;
    delete h_Nev;
}
