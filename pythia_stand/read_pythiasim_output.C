#include "PYTHIA_Muon_Analisys.h"
void particle_output();
void Muons_Analisys();
void read_pythiasim_output()
{
    Muons_Analisys();
}

void particle_output()
{

    TChain *tree = nullptr;

    // if (!gSystem->OpenDirectory(dir_filename.Data())) {
    //     printf("Dir not found\n");
    //     return tree;
    // }
    tree = new TChain("T1");
    //    printf("%s \n",Form("%s/%s",dir_filename.Data(),filename.Data()));
    tree->AddFile("~/pythia_stand_provAll_HFtest_Sim_SoftQCD_Def_10000_2710_DefaultBR.root");

    TClonesArray *data_charm = new TClonesArray("TLorentzVector", 400);
    TClonesArray *data_beauty = new TClonesArray("TLorentzVector", 400);
    TClonesArray *data_charm_prompt = new TClonesArray("TLorentzVector", 400);
    TClonesArray *data_charm_from_beauty = new TClonesArray("TLorentzVector", 400);
    TClonesArray *datamuminus_charm = new TClonesArray("TLorentzVector", 400);
    TClonesArray *datamuplus_charm = new TClonesArray("TLorentzVector", 400);
    TClonesArray *datamuminus_beauty = new TClonesArray("TLorentzVector", 400);
    TClonesArray *datamuplus_beauty = new TClonesArray("TLorentzVector", 400);

    int mu_plus_mother_charm[15][2];
    int mu_minus_mother_charm[15][2];
    int mu_plus_mother_beauty[15][2];
    int mu_minus_mother_beauty[15][2];
    vector<int> *code_charm = 0;
    TBranch *bcode_charm = 0;
    vector<int> *code_charm_prompt = 0;
    TBranch *bcode_charm_prompt = 0;
    vector<int> *code_charm_beauty = 0;
    TBranch *bcode_charm_beauty = 0;
    vector<int> *code_beauty = 0;
    TBranch *bcode_beauty = 0;

    tree->SetBranchAddress("charm_particle", &data_charm);
    tree->SetBranchAddress("code_charm", &code_charm, &bcode_charm);
    tree->SetBranchAddress("beauty_particle", &data_beauty);
    tree->SetBranchAddress("code_beauty", &code_beauty, &bcode_beauty);
    tree->SetBranchAddress("charm_particle_prompt", &data_charm_prompt);
    tree->SetBranchAddress("code_charm_prompt", &code_charm_prompt, &bcode_charm_prompt);
    tree->SetBranchAddress("charm_particle_from_beauty", &data_charm_from_beauty);
    tree->SetBranchAddress("code_charm_beauty", &code_charm_beauty, &bcode_charm_beauty);

    tree->SetBranchAddress("muon_minus_charm", &datamuminus_charm);
    tree->SetBranchAddress("muminus_charm", mu_minus_mother_charm);
    tree->SetBranchAddress("muon_plus_charm", &datamuplus_charm);
    tree->SetBranchAddress("muplus_charm", mu_plus_mother_charm);

    tree->SetBranchAddress("muon_minus_beauty", &datamuminus_beauty);
    tree->SetBranchAddress("muminus_beauty", mu_minus_mother_beauty);
    tree->SetBranchAddress("muon_plus_beauty", &datamuplus_beauty);
    tree->SetBranchAddress("muplus_beauty", mu_plus_mother_beauty);

    Int_t n_data_charm = 0;
    Int_t n_data_beauty = 0;
    Int_t n_data_charm_prompt = 0;
    Int_t n_data_charm_from_beauty = 0;
    Int_t n_datamuminus_charm = 0;
    Int_t n_datamuplus_charm = 0;
    Int_t n_datamuminus_beauty = 0;
    Int_t n_datamuplus_beauty = 0;

    for (int i = 0; i < tree->GetEntries(); i++)
    {

        if (i % (50000) == 0)
        {
            Printf("Evento: %d ", i);
        }

        tree->GetEvent(i);
        n_data_charm = data_charm->GetEntries();
        n_data_beauty = data_beauty->GetEntries();
        n_data_charm_prompt = data_charm_prompt->GetEntries();
        n_data_charm_from_beauty = data_charm_from_beauty->GetEntries();
        n_datamuminus_charm = datamuminus_charm->GetEntries();
        n_datamuplus_charm = datamuplus_charm->GetEntries();
        n_datamuminus_beauty = datamuminus_beauty->GetEntries();
        n_datamuplus_beauty = datamuplus_beauty->GetEntries();

        printf("%d\n", n_datamuplus_charm);

        for (int q = 0; q < n_data_charm; q++)
        {
            TLorentzVector *rec_part = (TLorentzVector *)data_charm->At(q);

            if (TMath::Abs(code_charm->at(q)) == 421)
            {

                h_ptDzero->Fill(rec_part->Pt());
                h_yDzero->Fill(rec_part->Rapidity());
            }
            if (TMath::Abs(code_charm->at(q)) == 411)
            {

                h_ptDplus->Fill(rec_part->Pt());
                h_yDplus->Fill(rec_part->Rapidity());
            }

            if (TMath::Abs(code_charm->at(q)) == 431)
            {

                h_ptDstrangeplus->Fill(rec_part->Pt());
                h_yDstrangeplus->Fill(rec_part->Rapidity());
            }

            if (TMath::Abs(code_charm->at(q)) == 4122)
            {

                h_ptLambdacharmedplus->Fill(rec_part->Pt());
                h_yLambdacharmedplus->Fill(rec_part->Rapidity());
            }

            if (TMath::Abs(code_charm->at(q)) > 400 && TMath::Abs(code_charm->at(q)) < 500)
            {
                h_ptAllCharm->Fill(rec_part->Pt());
                h_yAllCharm->Fill(rec_part->Rapidity());
                h_ptCharmMeson->Fill(rec_part->Pt());
                h_yCharmMeson->Fill(rec_part->Rapidity());
            }
            if (TMath::Abs(code_charm->at(q)) > 4000 && TMath::Abs(code_charm->at(q)) < 5000)
            {
                h_ptAllCharm->Fill(rec_part->Pt());
                h_yAllCharm->Fill(rec_part->Rapidity());
                h_ptCharmBarion->Fill(rec_part->Pt());
                h_yCharmBarion->Fill(rec_part->Rapidity());
            }
        }
    }

    TFile *fOut = new TFile("prova.root", "UPDATE");
    fOut->cd();

    h_ptDzero->Write(0, 2, 0);

    h_yDzero->Write(0, 2, 0);

    fOut->Close();
}

void Muons_Analisys()
{

    double LowPt = 0.0;
    double HighPt = 100.0;
    double LowY = -4.0;
    double HighY = -2.5;

    TString fileout;

    TChain *tree = nullptr;

    // if (!gSystem->OpenDirectory(dir_filename.Data())) {
    //     printf("Dir not found\n");
    //     return tree;
    // }
    tree = new TChain("T1");
    //    printf("%s \n",Form("%s/%s",dir_filename.Data(),filename.Data()));
    tree->AddFile("~/pythia_stand_provAll_HFtest_Sim_SoftQCD_Def_10000_2710_DefaultBR.root");

    TClonesArray *data_charm = new TClonesArray("TLorentzVector", 400);
    
    TClonesArray *data_beauty = new TClonesArray("TLorentzVector", 400);
    TClonesArray *data_charm_prompt = new TClonesArray("TLorentzVector", 400);
    TClonesArray *data_charm_from_beauty = new TClonesArray("TLorentzVector", 400);
    TClonesArray *datamuminus_charm = new TClonesArray("TLorentzVector", 400);
    TClonesArray *datamuplus_charm = new TClonesArray("TLorentzVector", 400);
    TClonesArray *datamuminus_beauty = new TClonesArray("TLorentzVector", 400);
    TClonesArray *datamuplus_beauty = new TClonesArray("TLorentzVector", 400);

    int mu_plus_mother_charm[15][2];
    int mu_minus_mother_charm[15][2];
    int mu_plus_mother_beauty[15][2];
    int mu_minus_mother_beauty[15][2];
    vector<int> *code_charm = 0;
    TBranch *bcode_charm = 0;
    vector<int> *code_charm_prompt = 0;
    TBranch *bcode_charm_prompt = 0;
    vector<int> *code_charm_beauty = 0;
    TBranch *bcode_charm_beauty = 0;
    vector<int> *code_beauty = 0;
    TBranch *bcode_beauty = 0;

    //    tree->SetBranchAddress("charm_particle",&data_charm);
    //    tree->SetBranchAddress("code_charm",&code_charm,&bcode_charm);
    //    tree->SetBranchAddress("beauty_particle",&data_beauty);
    //    tree->SetBranchAddress("code_beauty",&code_beauty,&bcode_beauty);
    //    tree->SetBranchAddress("charm_particle_prompt",&data_charm_prompt);
    //    tree->SetBranchAddress("code_charm_prompt",&code_charm_prompt,&bcode_charm_prompt);
    //    tree->SetBranchAddress("charm_particle_from_beauty",&data_charm_from_beauty);
    //    tree->SetBranchAddress("code_charm_beauty",&code_charm_beauty,&bcode_charm_beauty);

    tree->SetBranchAddress("muon_minus_charm", &datamuminus_charm);
    tree->SetBranchAddress("muminus_charm", mu_minus_mother_charm);
    tree->SetBranchAddress("muon_plus_charm", &datamuplus_charm);
    tree->SetBranchAddress("muplus_charm", mu_plus_mother_charm);

    tree->SetBranchAddress("muon_minus_beauty", &datamuminus_beauty);
    tree->SetBranchAddress("muminus_beauty", mu_minus_mother_beauty);
    tree->SetBranchAddress("muon_plus_beauty", &datamuplus_beauty);
    tree->SetBranchAddress("muplus_beauty", mu_plus_mother_beauty);

    TFile *outFile = new TFile("prova.root", "UPDATE");

    int n_datamuminus_charm = 0;
    int n_datamuplus_charm = 0;

    int n_datamuminus_beauty = 0;
    int n_datamuplus_beauty = 0;

    int Dimu_charm = 0;
    int Dimu_beauty = 0;

    //    double Ptmuplus=999.0;
    //    double Ptmuminus=999.0;
    //
    //    double Ymuplus=999.0;
    //    double Ymuminus=999.0;
    //
    //    double Phimuplus=999.0;
    //    double Phimuminus=999.0;

    bool verbose = false;
    
    for (int i = 0; i < tree->GetEntries(); i++)
    {

        if (i % (500000) == 0)
        {
            Printf("Evento: %d ", i);
        }

        tree->GetEvent(i);
        n_datamuminus_charm = datamuminus_charm->GetEntries();
        n_datamuplus_charm = datamuplus_charm->GetEntries();

        n_datamuminus_beauty = datamuminus_beauty->GetEntries();
        n_datamuplus_beauty = datamuplus_beauty->GetEntries();

        int totmuminus_charm = 0;
        int totmuplus_charm = 0;
        int totmuminus_beauty = 0;
        int totmuplus_beauty = 0;

        // Archiviazione mu+ da ccbar
        for (int q = 0; q < n_datamuplus_charm; q++)
        {

            TLorentzVector *rec_muplus = (TLorentzVector *)datamuplus_charm->At(q);
            h_yptMuon_plus_Charm->Fill(rec_muplus->Rapidity(), rec_muplus->Pt());

            double Ptmuplus = 999.0;
            double Ymuplus = 999.0;
            double Phimuplus = 999.0;

            Ptmuplus = rec_muplus->Pt();
            Ymuplus = rec_muplus->Rapidity();
            Phimuplus = rec_muplus->Phi();

            if ((Ptmuplus > LowPt && Ptmuplus < HighPt) && (Ymuplus > LowY && Ymuplus < HighY))
            {

                totmuplus_charm++;
                h_ptMuon_plus->Fill(rec_muplus->Pt());
                h_yMuon_plus->Fill(rec_muplus->Rapidity());

                h_ptMuon_plus_Charm->Fill(rec_muplus->Pt());
                h_yMuon_plus_Charm->Fill(rec_muplus->Rapidity());

                if ((TMath::Abs(mu_plus_mother_charm[q][0]) > 400 && TMath::Abs(mu_plus_mother_charm[q][0]) < 500) || (TMath::Abs(mu_plus_mother_charm[q][1]) > 400 && TMath::Abs(mu_plus_mother_charm[q][1]) < 500))
                {
                    h_ptMuon_plus_Charm_Meson->Fill(rec_muplus->Pt());
                    h_yMuon_plus_Charm_Meson->Fill(rec_muplus->Rapidity());

                    if (TMath::Abs(mu_plus_mother_charm[q][0]) == 421 || TMath::Abs(mu_plus_mother_charm[q][1]) == 421)
                    {
                        h_ptMuon_plus_Dzero->Fill(rec_muplus->Pt());
                        h_yMuon_plus_Dzero->Fill(rec_muplus->Rapidity());
                    }
                    else if (TMath::Abs(mu_plus_mother_charm[q][0]) == 411 || TMath::Abs(mu_plus_mother_charm[q][1]) == 411)
                    {
                        h_ptMuon_plus_Dplus->Fill(rec_muplus->Pt());
                        h_yMuon_plus_Dplus->Fill(rec_muplus->Rapidity());
                    }
                    else if (TMath::Abs(mu_plus_mother_charm[q][0]) == 431 || TMath::Abs(mu_plus_mother_charm[q][1]) == 431)
                    {
                        h_ptMuon_plus_Dstrangeplus->Fill(rec_muplus->Pt());
                        h_yMuon_plus_Dstrangeplus->Fill(rec_muplus->Rapidity());
                    }
                    else
                    {
                        h_ptMuon_plus_Charmrest->Fill(rec_muplus->Pt());
                        h_yMuon_plus_Charmrest->Fill(rec_muplus->Rapidity());
                    }
                }
                if ((TMath::Abs(mu_plus_mother_charm[q][0]) > 4000 && TMath::Abs(mu_plus_mother_charm[q][0]) < 5000) || (TMath::Abs(mu_plus_mother_charm[q][1]) > 4000 && TMath::Abs(mu_plus_mother_charm[q][1]) < 5000))
                {
                    h_ptMuon_plus_Charm_Barion->Fill(rec_muplus->Pt());
                    h_yMuon_plus_Charm_Barion->Fill(rec_muplus->Rapidity());
                    if (TMath::Abs(mu_plus_mother_charm[q][0]) == 4122 || TMath::Abs(mu_plus_mother_charm[q][1]) == 4122)
                    {
                        h_ptMuon_plus_Lambda->Fill(rec_muplus->Pt());
                        h_yMuon_plus_Lambda->Fill(rec_muplus->Rapidity());
                    }
                    else
                    {
                        h_ptMuon_plus_Charmrest->Fill(rec_muplus->Pt());
                        h_yMuon_plus_Charmrest->Fill(rec_muplus->Rapidity());
                    }
                }
            }
        }

        h_mult_mu_plus_charm->Fill(totmuplus_charm);

        // Archiviazione mu- e Dimuoni da ccbar
        for (int q = 0; q < n_datamuminus_charm; q++)
        {

            TLorentzVector *rec_muminus = (TLorentzVector *)datamuminus_charm->At(q);

            double Ptmuminus = 999.0;
            double Ymuminus = 999.0;
            double Phimuminus = 999.0;

            Ptmuminus = rec_muminus->Pt();
            Ymuminus = rec_muminus->Rapidity();
            Phimuminus = rec_muminus->Phi();

            if ((Ptmuminus > LowPt && Ptmuminus < HighPt) && (Ymuminus > LowY && Ymuminus < HighY))
            {

                totmuminus_charm++;
                h_ptMuon_minus->Fill(rec_muminus->Pt());
                h_yMuon_minus->Fill(rec_muminus->Rapidity());

                h_ptMuon_minus_Charm->Fill(rec_muminus->Pt());
                h_yMuon_minus_Charm->Fill(rec_muminus->Rapidity());

                h_yptMuon_minus_Charm->Fill(rec_muminus->Rapidity(), rec_muminus->Pt());

                if ((TMath::Abs(mu_minus_mother_charm[q][0]) > 400 && TMath::Abs(mu_minus_mother_charm[q][0]) < 500) || (TMath::Abs(mu_minus_mother_charm[q][1]) > 400 && TMath::Abs(mu_minus_mother_charm[q][1]) < 500))
                {

                    h_ptMuon_minus_Charm_Meson->Fill(rec_muminus->Pt());
                    h_yMuon_minus_Charm_Meson->Fill(rec_muminus->Rapidity());
                    if (TMath::Abs(mu_minus_mother_charm[q][0]) == 421 || TMath::Abs(mu_minus_mother_charm[q][1]) == 421)
                    {

                        h_ptMuon_minus_Dzero->Fill(rec_muminus->Pt());
                        h_yMuon_minus_Dzero->Fill(rec_muminus->Rapidity());
                    }
                    else if (TMath::Abs(mu_minus_mother_charm[q][0]) == 411 || TMath::Abs(mu_minus_mother_charm[q][1]) == 411)
                    {

                        h_ptMuon_minus_Dplus->Fill(rec_muminus->Pt());
                        h_yMuon_minus_Dplus->Fill(rec_muminus->Rapidity());
                    }
                    else if (TMath::Abs(mu_minus_mother_charm[q][0]) == 431 || TMath::Abs(mu_minus_mother_charm[q][1]) == 431)
                    {
                        h_ptMuon_minus_Dstrangeplus->Fill(rec_muminus->Pt());
                        h_yMuon_minus_Dstrangeplus->Fill(rec_muminus->Rapidity());
                    }
                    else
                    {
                        h_ptMuon_minus_Charmrest->Fill(rec_muminus->Pt());
                        h_yMuon_minus_Charmrest->Fill(rec_muminus->Rapidity());
                    }
                }

                if ((TMath::Abs(mu_minus_mother_charm[q][0]) > 4000 && TMath::Abs(mu_minus_mother_charm[q][0]) < 5000) || (TMath::Abs(mu_minus_mother_charm[q][1]) > 4000 && TMath::Abs(mu_minus_mother_charm[q][1]) < 5000))
                {

                    h_ptMuon_minus_Charm_Barion->Fill(rec_muminus->Pt());
                    h_yMuon_minus_Charm_Barion->Fill(rec_muminus->Rapidity());

                    if (TMath::Abs(mu_minus_mother_charm[q][0]) == 4122 || TMath::Abs(mu_minus_mother_charm[q][1]) == 4122)
                    {
                        h_ptMuon_minus_Lambda->Fill(rec_muminus->Pt());
                        h_yMuon_minus_Lambda->Fill(rec_muminus->Rapidity());
                    }
                    else
                    {
                        h_ptMuon_minus_Charmrest->Fill(rec_muminus->Pt());
                        h_yMuon_minus_Charmrest->Fill(rec_muminus->Rapidity());
                    }
                }

                for (int w = 0; w < n_datamuplus_charm; w++)
                {

                    TLorentzVector *rec_muplus = (TLorentzVector *)datamuplus_charm->At(w);

                    double Ptmuplus = 999.0;
                    double Phimuplus = 999.0;
                    double Ymuplus = 999.0;

                    Ptmuplus = rec_muplus->Pt();
                    Ymuplus = rec_muplus->Rapidity();
                    Phimuplus = rec_muplus->Phi();

                    if (Ymuplus > LowY && Ymuplus < HighY)
                    {

                        if (Ptmuplus > LowPt && Ptmuplus < HighPt)
                        {

                            TLorentzVector DiMu(*rec_muplus + *rec_muminus);
                            if (false)
                            {
                                cout << "-----------------------------------------" << endl;
                                cout << "evento " << i << endl;
                                cout << "Mu+: " << q << " Mu-: " << w << endl;
                                printf("pt Mu- (%0.1f) + pt Mu+ (%0.1f) = pt DiMu (%0.1f)\n", rec_muminus->Pt(), rec_muplus->Pt(), DiMu.Pt());

                                printf("y Mu- (%0.1f) + y Mu+ (%0.1f) = y DiMu (%0.1f)\n", rec_muminus->Rapidity(), rec_muplus->Rapidity(), DiMu.Rapidity());

                                printf("phi Mu- (%0.1f) + phi Mu+ (%0.1f) = phi DiMu (%0.1f)\n", rec_muminus->Phi(), rec_muplus->Phi(), DiMu.Phi());
                            }

                            Dimu_charm++;
                            h_ptDiMu->Fill(DiMu.Pt());
                            h_yDiMu->Fill(DiMu.Rapidity());
                            h_phiDiMu->Fill(DiMu.Phi());
                            h_mDiMu->Fill(DiMu.M());

                            h_ptDiMu_Charm->Fill(DiMu.Pt());
                            h_yDiMu_Charm->Fill(DiMu.Rapidity());
                            h_phiDiMu_Charm->Fill(DiMu.Phi());
                            h_mDiMu_Charm->Fill(DiMu.M());

                            if (DiMu.M() > 4.0)
                            {

                                h_ptDiMu4MassCut->Fill(DiMu.Pt());
                                h_yDiMu4MassCut->Fill(DiMu.Rapidity());
                                h_phiDiMu4MassCut->Fill(DiMu.Phi());
                                h_mDiMu4MassCut->Fill(DiMu.M());

                                h_ptDiMu_Charm4MassCut->Fill(DiMu.Pt());
                                h_yDiMu_Charm4MassCut->Fill(DiMu.Rapidity());
                                h_phiDiMu_Charm4MassCut->Fill(DiMu.Phi());
                                h_mDiMu_Charm4MassCut->Fill(DiMu.M());
                            }
                            if (DiMu.M() > 4.0 && DiMu.M() < 8.0)
                            {

                                h_ptDiMu_light->Fill(DiMu.Pt());
                                h_yDiMu_light->Fill(DiMu.Rapidity());
                                h_phiDiMu_light->Fill(DiMu.Phi());
                                h_mDiMu_light->Fill(DiMu.M());

                                h_ptDiMu_Charm_light->Fill(DiMu.Pt());
                                h_yDiMu_Charm_light->Fill(DiMu.Rapidity());
                                h_phiDiMu_Charm_light->Fill(DiMu.Phi());
                                h_mDiMu_Charm_light->Fill(DiMu.M());
                            }
                        }
                    }
                }
            }
        }
        h_mult_mu_minus_charm->Fill(totmuminus_charm);

        if (totmuminus_charm > 0 && totmuplus_charm > 0)
        {
            h_mult_mu_OSpair_charm->Fill(totmuminus_charm + totmuplus_charm);
        }

        // Archiviazione mu+ da bbbar
        for (int q = 0; q < n_datamuplus_beauty; q++)
        {

            TLorentzVector *rec_muplus_beauty = (TLorentzVector *)datamuplus_beauty->At(q);
            double Ptmuplus_beauty = 999.0;
            double Ymuplus_beauty = 999.0;
            double Phimuplus_beauty = 999.0;

            Ptmuplus_beauty = rec_muplus_beauty->Pt();
            Ymuplus_beauty = rec_muplus_beauty->Rapidity();
            Phimuplus_beauty = rec_muplus_beauty->Phi();

            if ((Ptmuplus_beauty > LowPt && Ptmuplus_beauty < HighPt) && (Ymuplus_beauty > LowY && Ymuplus_beauty < HighY))
            {
                totmuplus_beauty++;
                h_ptMuon_plus->Fill(rec_muplus_beauty->Pt());
                h_yMuon_plus->Fill(rec_muplus_beauty->Rapidity());
                //            cout<<"Evento: "<<i<<" mu più: "<<q<<endl;
                //            cout<<"Madre 1: "<<mu_plus_mother_beauty[q][0]<<endl;
                //            cout<<"Madre 2: "<<mu_plus_mother_beauty[q][1]<<endl;
                h_ptMuon_plus_AllBeauty->Fill(rec_muplus_beauty->Pt());
                h_yMuon_plus_AllBeauty->Fill(rec_muplus_beauty->Rapidity());

                if ((TMath::Abs(mu_plus_mother_beauty[q][0]) > 400 && TMath::Abs(mu_plus_mother_beauty[q][0]) < 500) || (TMath::Abs(mu_plus_mother_beauty[q][1]) > 400 && TMath::Abs(mu_plus_mother_beauty[q][1]) < 500))
                {

                    h_ptMuon_plus_BeautyCharm->Fill(rec_muplus_beauty->Pt());
                    h_yMuon_plus_BeautyCharm->Fill(rec_muplus_beauty->Rapidity());
                }
                else if ((TMath::Abs(mu_plus_mother_beauty[q][0]) > 4000 && TMath::Abs(mu_plus_mother_beauty[q][0]) < 5000) || (TMath::Abs(mu_plus_mother_beauty[q][1]) > 4000 && TMath::Abs(mu_plus_mother_beauty[q][1]) < 5000))
                {
                    h_ptMuon_plus_BeautyCharm->Fill(rec_muplus_beauty->Pt());
                    h_yMuon_plus_BeautyCharm->Fill(rec_muplus_beauty->Rapidity());
                }
                else if ((TMath::Abs(mu_plus_mother_beauty[q][0]) > 500 && TMath::Abs(mu_plus_mother_beauty[q][0]) < 600) || (TMath::Abs(mu_plus_mother_beauty[q][1]) > 500 && TMath::Abs(mu_plus_mother_beauty[q][1]) < 600))
                {
                    h_ptMuon_plus_Beauty->Fill(rec_muplus_beauty->Pt());
                    h_yMuon_plus_Beauty->Fill(rec_muplus_beauty->Rapidity());
                }

                else if ((TMath::Abs(mu_plus_mother_beauty[q][0]) > 5000 && TMath::Abs(mu_plus_mother_beauty[q][0]) < 6000) || (TMath::Abs(mu_plus_mother_beauty[q][1]) > 5000 && TMath::Abs(mu_plus_mother_beauty[q][1]) < 6000))
                {
                    h_ptMuon_plus_Beauty->Fill(rec_muplus_beauty->Pt());
                    h_yMuon_plus_Beauty->Fill(rec_muplus_beauty->Rapidity());
                }
            }
        }

        h_mult_mu_plus_beauty->Fill(totmuplus_beauty);

        // Archiviazione mu- e Dimuoni da bbbar
        //        printf("eveto %d n_datamuminus_beauty %d\n",i,n_datamuminus_beauty);
        for (int q = 0; q < n_datamuminus_beauty; q++)
        {

            TLorentzVector *rec_muminus_beauty = (TLorentzVector *)datamuminus_beauty->At(q);

            double Ptmuminus_beauty = 999.0;
            double Ymuminus_beauty = 999.0;
            double Phimuminus_beauty = 999.0;

            Ptmuminus_beauty = rec_muminus_beauty->Pt();
            Ymuminus_beauty = rec_muminus_beauty->Rapidity();
            Phimuminus_beauty = rec_muminus_beauty->Phi();

            if ((Ptmuminus_beauty > LowPt && Ptmuminus_beauty < HighPt) && (Ymuminus_beauty > LowY && Ymuminus_beauty < HighY))
            {
                totmuminus_beauty++;
                h_ptMuon_minus->Fill(rec_muminus_beauty->Pt());
                h_yMuon_minus->Fill(rec_muminus_beauty->Rapidity());

                h_ptMuon_minus_AllBeauty->Fill(rec_muminus_beauty->Pt());
                h_yMuon_minus_AllBeauty->Fill(rec_muminus_beauty->Rapidity());

                if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 500))
                {
                    //                cout<<"Evento: "<<i<<" mu meno da charm: "<<q<<endl;
                    //                cout<<"Madre 1: "<<mu_minus_mother_beauty[q][0]<<endl;
                    //                cout<<"Madre 2: "<<mu_minus_mother_beauty[q][1]<<endl;
                    h_ptMuon_minus_BeautyCharm->Fill(rec_muminus_beauty->Pt());
                    h_yMuon_minus_BeautyCharm->Fill(rec_muminus_beauty->Rapidity());
                }
                else if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 5000))
                {
                    h_ptMuon_minus_BeautyCharm->Fill(rec_muminus_beauty->Pt());
                    h_yMuon_minus_BeautyCharm->Fill(rec_muminus_beauty->Rapidity());
                }
                else if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 600))
                {
                    //                cout<<"Evento: "<<i<<" mu meno da beauty: "<<q<<endl;
                    //                cout<<"Madre 1: "<<mu_minus_mother_beauty[q][0]<<endl;
                    //                cout<<"Madre 2: "<<mu_minus_mother_beauty[q][1]<<endl;
                    h_ptMuon_minus_Beauty->Fill(rec_muminus_beauty->Pt());
                    h_yMuon_minus_Beauty->Fill(rec_muminus_beauty->Rapidity());
                }
                else if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 6000))
                {
                    h_ptMuon_minus_Beauty->Fill(rec_muminus_beauty->Pt());
                    h_yMuon_minus_Beauty->Fill(rec_muminus_beauty->Rapidity());
                }

                for (int w = 0; w < n_datamuplus_beauty; w++)
                {
                    TLorentzVector *rec_muplus_beauty = (TLorentzVector *)datamuplus_beauty->At(w);

                    double Ptmuplus_beauty = 999.0;
                    double Ymuplus_beauty = 999.0;
                    double Phimuplus_beauty = 999.0;

                    Ptmuplus_beauty = rec_muplus_beauty->Pt();
                    Ymuplus_beauty = rec_muplus_beauty->Rapidity();
                    Phimuplus_beauty = rec_muplus_beauty->Phi();

                    //-----------Primo Taglio in Rapidity-----------//
                    if (Ymuplus_beauty > LowY && Ymuplus_beauty < HighY)
                    {

                        if (Ptmuplus_beauty > LowPt && Ptmuplus_beauty < HighPt)
                        {

                            Dimu_beauty++;
                            TLorentzVector DiMu(*rec_muplus_beauty + *rec_muminus_beauty);
                            if (false)
                            {
                                cout << "Mu+: " << q << " Mu-: " << w << endl;
                                printf("pt Mu- (%0.1f) + pt Mu+ (%0.1f) = pt DiMu (%0.1f)\n", rec_muminus_beauty->Pt(), rec_muplus_beauty->Pt(), DiMu.Pt());

                                printf("y Mu- (%0.1f) + y Mu+ (%0.1f) = y DiMu (%0.1f)\n", rec_muminus_beauty->Rapidity(), rec_muplus_beauty->Rapidity(), DiMu.Rapidity());

                                printf("phi Mu- (%0.1f) + phi Mu+ (%0.1f) = phi DiMu (%0.1f)\n", rec_muminus_beauty->Phi(), rec_muplus_beauty->Phi(), DiMu.Phi());
                            }

                            h_ptDiMu->Fill(DiMu.Pt());
                            h_yDiMu->Fill(DiMu.Rapidity());
                            h_phiDiMu->Fill(DiMu.Phi());
                            h_mDiMu->Fill(DiMu.M());

                            h_ptDiMu_AllBeauty->Fill(DiMu.Pt());
                            h_yDiMu_AllBeauty->Fill(DiMu.Rapidity());
                            h_phiDiMu_AllBeauty->Fill(DiMu.Phi());
                            h_mDiMu_AllBeauty->Fill(DiMu.M());

                            if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[q][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 600) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 6000))
                            {

                                if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 6000))
                                {

                                    h_ptDiMu_Beauty->Fill(DiMu.Pt());
                                    h_yDiMu_Beauty->Fill(DiMu.Rapidity());
                                    h_phiDiMu_Beauty->Fill(DiMu.Phi());
                                    h_mDiMu_Beauty->Fill(DiMu.M());
                                }

                                if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 5000))
                                {

                                    h_ptDiMu_BeautyCharm->Fill(DiMu.Pt());
                                    h_yDiMu_BeautyCharm->Fill(DiMu.Rapidity());
                                    h_phiDiMu_BeautyCharm->Fill(DiMu.Phi());
                                    h_mDiMu_BeautyCharm->Fill(DiMu.M());
                                }
                            }
                            else if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[q][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 500) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 5000))
                            {

                                if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 6000))
                                {

                                    h_ptDiMu_BeautyCharm->Fill(DiMu.Pt());
                                    h_yDiMu_BeautyCharm->Fill(DiMu.Rapidity());
                                    h_phiDiMu_BeautyCharm->Fill(DiMu.Phi());
                                    h_mDiMu_BeautyCharm->Fill(DiMu.M());
                                }

                                if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 5000))
                                {

                                    h_ptDiMu_BeautyCharm->Fill(DiMu.Pt());
                                    h_yDiMu_BeautyCharm->Fill(DiMu.Rapidity());
                                    h_phiDiMu_BeautyCharm->Fill(DiMu.Phi());
                                    h_mDiMu_BeautyCharm->Fill(DiMu.M());
                                }
                            }

                            if (DiMu.M() > 4.0)
                            {

                                h_ptDiMu4MassCut->Fill(DiMu.Pt());
                                h_yDiMu4MassCut->Fill(DiMu.Rapidity());
                                h_phiDiMu4MassCut->Fill(DiMu.Phi());
                                h_mDiMu4MassCut->Fill(DiMu.M());

                                h_ptDiMu_AllBeauty4MassCut->Fill(DiMu.Pt());
                                h_yDiMu_AllBeauty4MassCut->Fill(DiMu.Rapidity());
                                h_phiDiMu_AllBeauty4MassCut->Fill(DiMu.Phi());
                                h_mDiMu_AllBeauty4MassCut->Fill(DiMu.M());

                                if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[q][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 600) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 6000))
                                {

                                    if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 6000))
                                    {

                                        h_ptDiMu_Beauty4MassCut->Fill(DiMu.Pt());
                                        h_yDiMu_Beauty4MassCut->Fill(DiMu.Rapidity());
                                        h_phiDiMu_Beauty4MassCut->Fill(DiMu.Phi());
                                        h_mDiMu_Beauty4MassCut->Fill(DiMu.M());
                                    }

                                    if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 5000))
                                    {

                                        h_ptDiMu_BeautyCharm4MassCut->Fill(DiMu.Pt());
                                        h_yDiMu_BeautyCharm4MassCut->Fill(DiMu.Rapidity());
                                        h_phiDiMu_BeautyCharm4MassCut->Fill(DiMu.Phi());
                                        h_mDiMu_BeautyCharm4MassCut->Fill(DiMu.M());
                                    }
                                }
                                else if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[q][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 500) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 5000))
                                {

                                    if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 6000))
                                    {

                                        h_ptDiMu_BeautyCharm4MassCut->Fill(DiMu.Pt());
                                        h_yDiMu_BeautyCharm4MassCut->Fill(DiMu.Rapidity());
                                        h_phiDiMu_BeautyCharm4MassCut->Fill(DiMu.Phi());
                                        h_mDiMu_BeautyCharm4MassCut->Fill(DiMu.M());
                                    }

                                    if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 5000))
                                    {

                                        h_ptDiMu_BeautyCharm4MassCut->Fill(DiMu.Pt());
                                        h_yDiMu_BeautyCharm4MassCut->Fill(DiMu.Rapidity());
                                        h_phiDiMu_BeautyCharm4MassCut->Fill(DiMu.Phi());
                                        h_mDiMu_BeautyCharm4MassCut->Fill(DiMu.M());
                                    }
                                }
                            }

                            if (DiMu.M() > 4.0 && DiMu.M() < 8.0)
                            {

                                h_ptDiMu_light->Fill(DiMu.Pt());
                                h_yDiMu_light->Fill(DiMu.Rapidity());
                                h_phiDiMu_light->Fill(DiMu.Phi());
                                h_mDiMu_light->Fill(DiMu.M());

                                h_ptDiMu_AllBeauty_light->Fill(DiMu.Pt());
                                h_yDiMu_AllBeauty_light->Fill(DiMu.Rapidity());
                                h_phiDiMu_AllBeauty_light->Fill(DiMu.Phi());
                                h_mDiMu_AllBeauty_light->Fill(DiMu.M());

                                if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[q][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 600) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 6000))
                                {

                                    if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 6000))
                                    {

                                        h_ptDiMu_Beauty_light->Fill(DiMu.Pt());
                                        h_yDiMu_Beauty_light->Fill(DiMu.Rapidity());
                                        h_phiDiMu_Beauty_light->Fill(DiMu.Phi());
                                        h_mDiMu_Beauty_light->Fill(DiMu.M());
                                    }

                                    if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 5000))
                                    {

                                        h_ptDiMu_BeautyCharm_light->Fill(DiMu.Pt());
                                        h_yDiMu_BeautyCharm_light->Fill(DiMu.Rapidity());
                                        h_phiDiMu_BeautyCharm_light->Fill(DiMu.Phi());
                                        h_mDiMu_BeautyCharm_light->Fill(DiMu.M());
                                    }
                                }
                                else if ((TMath::Abs(mu_minus_mother_beauty[q][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[q][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[q][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 500) || (TMath::Abs(mu_minus_mother_beauty[q][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[q][1]) < 5000))
                                {

                                    if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 6000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 500 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 600) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 5000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 6000))
                                    {

                                        h_ptDiMu_BeautyCharm_light->Fill(DiMu.Pt());
                                        h_yDiMu_BeautyCharm_light->Fill(DiMu.Rapidity());
                                        h_phiDiMu_BeautyCharm_light->Fill(DiMu.Phi());
                                        h_mDiMu_BeautyCharm_light->Fill(DiMu.M());
                                    }

                                    if ((TMath::Abs(mu_minus_mother_beauty[w][0]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][0]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][0]) < 5000) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 400 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 500) || (TMath::Abs(mu_minus_mother_beauty[w][1]) > 4000 && TMath::Abs(mu_minus_mother_beauty[w][1]) < 5000))
                                    {

                                        h_ptDiMu_BeautyCharm_light->Fill(DiMu.Pt());
                                        h_yDiMu_BeautyCharm_light->Fill(DiMu.Rapidity());
                                        h_phiDiMu_BeautyCharm_light->Fill(DiMu.Phi());
                                        h_mDiMu_BeautyCharm_light->Fill(DiMu.M());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        h_mult_mu_minus_beauty->Fill(totmuminus_beauty);

        if (totmuminus_beauty > 0 && totmuplus_beauty > 0)
        {
            h_mult_mu_OSpair_beauty->Fill(totmuminus_beauty + totmuplus_beauty);
        }
    }

    printf("\n DiMuons from Charm %d\n", Dimu_charm);
    printf("\n DiMuons from Beauty %d\n", Dimu_beauty);

    outFile->cd();

 Writing_Muons_Analisys(outFile, "name_folder", "name_subfolder");

    outFile->Close();
}



void Writing_Muons_Analisys(TFile *outFile,TString name_folder,TString name_subfolder){
    
    TDirectory *subsub1 = outFile->GetDirectory(name_subfolder.Data());

    if (!subsub1) subsub1 = outFile->mkdir(name_subfolder.Data());
    else printf("\n %s già esistente\n",name_subfolder.Data());

    outFile->cd(name_subfolder.Data());
    
    h_mult_mu_plus_charm->GetXaxis()->SetTitle("N° #mu^{+}");
    h_mult_mu_plus_charm->Write(0,2,0);
    
    h_mult_mu_minus_charm->GetXaxis()->SetTitle("N° #mu^{-}");
    h_mult_mu_minus_charm->Write(0,2,0);
    
    h_ptMuon_plus_Charm->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_plus_Charm->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_plus_Charm->Write(0,2,0);
    
    h_ptMuon_plus_Charm_Meson->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_plus_Charm_Meson->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_plus_Charm_Meson->Write(0,2,0);
    
    h_ptMuon_plus_Charm_Barion->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_plus_Charm_Barion->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_plus_Charm_Barion->Write(0,2,0);
    
    h_ptMuon_plus_Charmrest->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_plus_Charmrest->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_plus_Charmrest->Write(0,2,0);
    
    h_ptMuon_plus_Dzero->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_plus_Dzero->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_plus_Dzero->Write(0,2,0);
    
    h_ptMuon_plus_Dplus->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_plus_Dplus->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_plus_Dplus->Write(0,2,0);
    
    h_ptMuon_plus_Dstrangeplus->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_plus_Dstrangeplus->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_plus_Dstrangeplus->Write(0,2,0);
    
    h_ptMuon_plus_Lambda->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_plus_Lambda->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_plus_Lambda->Write(0,2,0);
    
    h_yMuon_plus_Charm->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_Charm->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_Charm->Write(0,2,0);
    
    h_yMuon_plus_Charm_Meson->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_Charm_Meson->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_Charm_Meson->Write(0,2,0);
    
    h_yMuon_plus_Charm_Barion->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_Charm_Barion->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_Charm_Barion->Write(0,2,0);
    
    h_yMuon_plus_Charmrest->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_Charmrest->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_Charmrest->Write(0,2,0);
    
    h_yMuon_plus_Dzero->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_Dzero->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_Dzero->Write(0,2,0);
    
    h_yMuon_plus_Dplus->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_Dplus->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_Dplus->Write(0,2,0);
    
    h_yMuon_plus_Dstrangeplus->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_Dstrangeplus->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_Dstrangeplus->Write(0,2,0);
    
    h_yMuon_plus_Dstrangeplus->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_Dstrangeplus->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_Dstrangeplus->Write(0,2,0);
    
    h_ptMuon_minus_Charm->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_minus_Charm->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_minus_Charm->Write(0,2,0);
    
    h_ptMuon_minus_Charm_Meson->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_minus_Charm_Meson->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_minus_Charm_Meson->Write(0,2,0);
    
    h_ptMuon_minus_Charm_Barion->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_minus_Charm_Barion->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_minus_Charm_Barion->Write(0,2,0);
    
    h_ptMuon_minus_Charmrest->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_minus_Charmrest->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_minus_Charmrest->Write(0,2,0);
    
    h_ptMuon_minus_Dzero->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_minus_Dzero->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_minus_Dzero->Write(0,2,0);
    
    h_ptMuon_minus_Dplus->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_minus_Dplus->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_minus_Dplus->Write(0,2,0);
    
    h_ptMuon_minus_Dstrangeplus->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_minus_Dstrangeplus->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_minus_Dstrangeplus->Write(0,2,0);
    
    h_ptMuon_minus_Lambda->GetXaxis()->SetTitle("#it{p}_{t}");
    h_ptMuon_minus_Lambda->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_ptMuon_minus_Lambda->Write(0,2,0);
    
    h_yMuon_minus_Charm->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_Charm->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_Charm->Write(0,2,0);
    
    h_yMuon_minus_Charm_Meson->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_Charm_Meson->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_Charm_Meson->Write(0,2,0);
    
    h_yMuon_minus_Charm_Barion->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_Charm_Barion->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_Charm_Barion->Write(0,2,0);
    
    h_yMuon_minus_Charmrest->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_Charmrest->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_Charmrest->Write(0,2,0);
    
    h_yMuon_minus_Dzero->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_Dzero->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_Dzero->Write(0,2,0);
    
    h_yMuon_minus_Dplus->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_Dplus->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_Dplus->Write(0,2,0);
    
    h_yMuon_minus_Dstrangeplus->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_Dstrangeplus->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_Dstrangeplus->Write(0,2,0);
    
    h_yMuon_minus_Dstrangeplus->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_Dstrangeplus->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_Dstrangeplus->Write(0,2,0);
    
    name_subfolder.Form("%s/Muons_from_Beauty",name_folder.Data());

    TDirectory *subsub2 = outFile->GetDirectory(name_subfolder.Data());

    if (!subsub2) subsub2 = outFile->mkdir(name_subfolder.Data());
    else printf("\n %s già esistente\n",name_subfolder.Data());

    outFile->cd(name_subfolder.Data());
    h_mult_mu_plus_beauty->GetXaxis()->SetTitle("N° #mu^{+}");
    h_mult_mu_plus_beauty->Write(0,2,0);
    
    h_mult_mu_minus_beauty->GetXaxis()->SetTitle("N° #mu^{-}");
    h_mult_mu_minus_beauty->Write(0,2,0);
    
    h_ptMuon_plus_AllBeauty->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptMuon_plus_AllBeauty->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptMuon_plus_AllBeauty->Write(0,2,0);
    
    h_ptMuon_plus_Beauty->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptMuon_plus_Beauty->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptMuon_plus_Beauty->Write(0,2,0);
    
    h_ptMuon_plus_BeautyCharm->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptMuon_plus_BeautyCharm->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptMuon_plus_BeautyCharm->Write(0,2,0);
    
    h_yMuon_plus_AllBeauty->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_AllBeauty->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_AllBeauty->Write(0,2,0);
    
    h_yMuon_plus_Beauty->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_Beauty->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_Beauty->Write(0,2,0);
    
    h_yMuon_plus_BeautyCharm->GetXaxis()->SetTitle("Y");
    h_yMuon_plus_BeautyCharm->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_plus_BeautyCharm->Write(0,2,0);
    
    h_ptMuon_minus_AllBeauty->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptMuon_minus_AllBeauty->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptMuon_minus_AllBeauty->Write(0,2,0);
    
    h_ptMuon_minus_Beauty->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptMuon_minus_Beauty->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptMuon_minus_Beauty->Write(0,2,0);
    
    h_ptMuon_minus_BeautyCharm->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptMuon_minus_BeautyCharm->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptMuon_minus_BeautyCharm->Write(0,2,0);
    
    h_yMuon_minus_AllBeauty->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_AllBeauty->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_AllBeauty->Write(0,2,0);
    
    h_yMuon_minus_Beauty->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_Beauty->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_Beauty->Write(0,2,0);
    
    h_yMuon_minus_BeautyCharm->GetXaxis()->SetTitle("Y");
    h_yMuon_minus_BeautyCharm->GetYaxis()->SetTitle("dN/dY");
    h_yMuon_minus_BeautyCharm->Write(0,2,0);
    
    name_subfolder.Form("%s/DiMu",name_folder.Data());

    TDirectory *subsub3 = outFile->GetDirectory(name_subfolder.Data());

    if (!subsub3) subsub3 = outFile->mkdir(name_subfolder.Data());
    else printf("\n %s già esistente\n",name_subfolder.Data());

    outFile->cd(name_subfolder.Data());
    
    h_mult_mu_OSpair_charm->GetXaxis()->SetTitle("N° #mu^{+}+#mu^{-}");
    h_mult_mu_OSpair_charm->Write(0,2,0);
    
    h_mult_mu_OSpair_beauty->GetXaxis()->SetTitle("N° #mu^{+}+#mu^{-}");
    h_mult_mu_OSpair_beauty->Write(0,2,0);
    
    h_ptDiMu->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu->Write(0,2,0);
    
    h_yDiMu->Write(0,2,0);
    h_phiDiMu->Write(0,2,0);
    
    h_mDiMu->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu->Write(0,2,0);
    
    h_ptDiMu_Charm->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_Charm->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_Charm->Write(0,2,0);
    
    h_yDiMu_Charm->Write(0,2,0);
    h_phiDiMu_Charm->Write(0,2,0);
    
    h_mDiMu_Charm->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_Charm->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_Charm->Write(0,2,0);
    
    h_ptDiMu_AllBeauty->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_AllBeauty->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_AllBeauty->Write(0,2,0);
    
    h_yDiMu_AllBeauty->Write(0,2,0);
    h_phiDiMu_AllBeauty->Write(0,2,0);
    
    h_mDiMu_AllBeauty->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_AllBeauty->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_AllBeauty->Write(0,2,0);
    
    h_ptDiMu_Beauty->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_Beauty->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_Beauty->Write(0,2,0);
    
    h_yDiMu_Beauty->Write(0,2,0);
    h_phiDiMu_Beauty->Write(0,2,0);
    
    h_mDiMu_Beauty->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_Beauty->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_Beauty->Write(0,2,0);
    
    h_ptDiMu_BeautyCharm->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_BeautyCharm->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_BeautyCharm->Write(0,2,0);
    
    
    h_yDiMu_BeautyCharm->Write(0,2,0);
    h_phiDiMu_BeautyCharm->Write(0,2,0);
   
    h_mDiMu_BeautyCharm->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_BeautyCharm->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_BeautyCharm->Write(0,2,0);
   
    h_ptDiMu4MassCut->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu4MassCut->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu4MassCut->Write(0,2,0);
    
    h_yDiMu4MassCut->Write(0,2,0);
    h_phiDiMu4MassCut->Write(0,2,0);
    
    h_mDiMu4MassCut->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu4MassCut->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu4MassCut->Write(0,2,0);
    
    h_ptDiMu_Charm4MassCut->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_Charm4MassCut->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_Charm4MassCut->Write(0,2,0);
    
    h_yDiMu_Charm4MassCut->Write(0,2,0);
    h_phiDiMu_Charm4MassCut->Write(0,2,0);
    
    h_mDiMu_Charm4MassCut->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_Charm4MassCut->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_Charm4MassCut->Write(0,2,0);
    
    h_ptDiMu_AllBeauty4MassCut->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_AllBeauty4MassCut->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_AllBeauty4MassCut->Write(0,2,0);
    
    h_yDiMu_AllBeauty4MassCut->Write(0,2,0);
    h_phiDiMu_AllBeauty4MassCut->Write(0,2,0);
    
    h_mDiMu_AllBeauty4MassCut->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_AllBeauty4MassCut->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_AllBeauty4MassCut->Write(0,2,0);
    
    h_ptDiMu_Beauty4MassCut->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_Beauty4MassCut->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_Beauty4MassCut->Write(0,2,0);
    
    h_yDiMu_Beauty4MassCut->Write(0,2,0);
    h_phiDiMu_Beauty4MassCut->Write(0,2,0);
    
    h_mDiMu_Beauty4MassCut->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_Beauty4MassCut->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_Beauty4MassCut->Write(0,2,0);
    
    h_ptDiMu_BeautyCharm4MassCut->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_BeautyCharm4MassCut->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_BeautyCharm4MassCut->Write(0,2,0);
    
    h_yDiMu_BeautyCharm4MassCut->Write(0,2,0);
    h_phiDiMu_BeautyCharm4MassCut->Write(0,2,0);
    
    h_mDiMu_BeautyCharm4MassCut->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_BeautyCharm4MassCut->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_BeautyCharm4MassCut->Write(0,2,0);
    
    h_ptDiMu_light->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_light->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_light->Write(0,2,0);
    
    h_yDiMu_light->Write(0,2,0);
    h_phiDiMu_light->Write(0,2,0);
    
    h_mDiMu_light->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_light->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_light->Write(0,2,0);
    
    h_ptDiMu_Charm_light->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_Charm_light->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_Charm_light->Write(0,2,0);
    
    h_yDiMu_Charm_light->Write(0,2,0);
    h_phiDiMu_Charm_light->Write(0,2,0);
   
    h_mDiMu_Charm_light->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_Charm_light->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_Charm_light->Write(0,2,0);
   
    h_ptDiMu_AllBeauty_light->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_AllBeauty_light->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_AllBeauty_light->Write(0,2,0);
    
    h_yDiMu_AllBeauty_light->Write(0,2,0);
    h_phiDiMu_AllBeauty_light->Write(0,2,0);
    
    h_mDiMu_AllBeauty_light->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_AllBeauty_light->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_AllBeauty_light->Write(0,2,0);
    
    h_ptDiMu_Beauty_light->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_Beauty_light->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_Beauty_light->Write(0,2,0);
    
    h_yDiMu_Beauty_light->Write(0,2,0);
    h_phiDiMu_Beauty_light->Write(0,2,0);
    
    h_mDiMu_Beauty_light->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_Beauty_light->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_Beauty_light->Write(0,2,0);
    
    h_ptDiMu_BeautyCharm_light->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
    h_ptDiMu_BeautyCharm_light->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
    h_ptDiMu_BeautyCharm_light->Write(0,2,0);
    
    h_yDiMu_BeautyCharm_light->Write(0,2,0);
    h_phiDiMu_BeautyCharm_light->Write(0,2,0);
    
    h_mDiMu_BeautyCharm_light->GetXaxis()->SetTitle("M [GeV/c^{2}]");
    h_mDiMu_BeautyCharm_light->GetYaxis()->SetTitle("dN/dM [GeV/c^{2}]^{-1}");
    h_mDiMu_BeautyCharm_light->Write(0,2,0);
    
}

// void Muons(){
//     char text[400];
    
//     TH1D *Mu_Charm=onehist("/Volumes/REPOSITORY/Result/Hist_HFtest_Sim_SoftQCD_Def_DefaultBR_150M.root","PtCut_1.0_20.0_YCut_4.0_2.5/Muons_from_Charm/h_ptMuon_plus_Charm",true);
//     sprintf(text,"#mu^{+} from prompt charm");
//     Mu_Charm->SetTitle(text);
//     Mu_Charm->GetXaxis()->SetTitle("#it{p}_{t} [GeV/c]");
//     Mu_Charm->GetYaxis()->SetTitle("dN/d#it{p}_{t} [GeV/c]^{-1}");
//     TH1D *Mu_Beauty=onehist("/Volumes/REPOSITORY/Result/Hist_HFtest_Sim_SoftQCD_Def_DefaultBR_150M.root","PtCut_1.0_20.0_YCut_4.0_2.5/Muons_from_Beauty/h_ptMuon_plus_AllBeauty",true);
//     sprintf(text,"#mu^{+} from beauty");
//     Mu_Beauty->SetTitle(text);
    
//     sprintf(text,"#splitline{#splitline{p-p collision,#sqrt{s} 13=TeV}{150M PYTHIA8 SoftQCD ND Monash DefaultBR}}{#mu^{+} 1<#it{p}_{T}<20 GeV/c -4.0<y<-2.5}");
    
//     TCanvas *c=allprint2hist_marker(Mu_Charm,Mu_Beauty,text,false,true);
//     TFile* outFile=new TFile("pippo.root","UPDATE");
//     c->Write(0,2,0);
//     outFile->Close();
    
    
// }