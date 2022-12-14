#include "HistMCDimuonHF.h"

TChain *Getting_Tree(TString dir_filename, TString filename, bool from_data)
{
  TChain *tree = nullptr;

  // if (!gSystem->OpenDirectory(dir_filename.Data())) {
  //     printf("Dir not found\n");
  //     return tree;
  // }
  tree = new TChain("MCTree");
  //    printf("%s \n",Form("%s/%s",dir_filename.Data(),filename.Data()));
  tree->AddFile(Form("%s/%s", dir_filename.Data(), filename.Data()));

  tree->SetBranchAddress("N_HFquarks_gen", &N_HFquarks_gen);

  tree->SetBranchAddress("PDG_HFquark_gen", PDG_HFquark_gen);
  tree->SetBranchAddress("Pt_HFquark_gen", Pt_HFquark_gen);
  tree->SetBranchAddress("Y_HFquark_gen", Y_HFquark_gen);

  // tree->SetBranchAddress("N_HFquarks_rec", &N_HFquarks_rec);

  // tree->SetBranchAddress("PDG_HFquark_rec", PDG_HFquark_rec);
  // tree->SetBranchAddress("Pt_HFquark_rec", Pt_HFquark_rec);
  // tree->SetBranchAddress("Y_HFquark_rec", Y_HFquark_rec);

  tree->SetBranchAddress("NDimu_gen", &NDimu_gen);
  tree->SetBranchAddress("DimuMu_gen", DimuMu_gen);
  tree->SetBranchAddress("DimuPt_gen", DimuPt_gen);
  tree->SetBranchAddress("DimuPx_gen", DimuPx_gen);
  tree->SetBranchAddress("DimuPy_gen", DimuPy_gen);
  tree->SetBranchAddress("DimuPz_gen", DimuPz_gen);
  tree->SetBranchAddress("DimuY_gen", DimuY_gen);
  tree->SetBranchAddress("DimuMass_gen", DimuMass_gen);
  tree->SetBranchAddress("DimuCharge_gen", DimuCharge_gen);

  tree->SetBranchAddress("NMuons_gen", &NMuons_gen);
  tree->SetBranchAddress("PDGmum_gen", PDGmum_gen);
  tree->SetBranchAddress("Pt_gen", Pt_gen);
  tree->SetBranchAddress("E_gen", E_gen);
  tree->SetBranchAddress("Px_gen", Px_gen);
  tree->SetBranchAddress("Py_gen", Py_gen);
  tree->SetBranchAddress("Pz_gen", Pz_gen);
  tree->SetBranchAddress("Y_gen", Y_gen);
  tree->SetBranchAddress("Eta_gen", Eta_gen);
  tree->SetBranchAddress("Phi_gen", Phi_gen);
  tree->SetBranchAddress("Theta_gen", Theta_gen);

  tree->SetBranchAddress("NDimu_rec", &NDimu_rec);
  tree->SetBranchAddress("DimuMu_rec", DimuMu_rec);
  tree->SetBranchAddress("DimuPt_rec", DimuPt_rec);
  tree->SetBranchAddress("DimuPx_rec", DimuPx_rec);
  tree->SetBranchAddress("DimuPy_rec", DimuPy_rec);
  tree->SetBranchAddress("DimuPz_rec", DimuPz_rec);
  tree->SetBranchAddress("DimuY_rec", DimuY_rec);
  tree->SetBranchAddress("DimuMass_rec", DimuMass_rec);
  tree->SetBranchAddress("DimuCharge_rec", DimuCharge_rec);
  tree->SetBranchAddress("DimuMatch_rec", DimuMatch_rec);
  tree->SetBranchAddress("DimuPhi_rec", DimuPhi_rec);
  tree->SetBranchAddress("DimuTheta_rec", DimuTheta_rec);

  tree->SetBranchAddress("NMuons_rec", &NMuons_rec);
  tree->SetBranchAddress("PDGmum_rec", PDGmum_rec);
  tree->SetBranchAddress("E_rec", E_rec);
  tree->SetBranchAddress("Px_rec", Px_rec);
  tree->SetBranchAddress("Pt_rec", Pt_rec);
  tree->SetBranchAddress("Py_rec", Py_rec);
  tree->SetBranchAddress("Pz_rec", Pz_rec);
  tree->SetBranchAddress("Y_rec", Y_rec);
  tree->SetBranchAddress("Eta_rec", Eta_rec);
  tree->SetBranchAddress("MatchTrig_rec", MatchTrig_rec);
  tree->SetBranchAddress("TrackChi2_rec", TrackChi2_rec);
  tree->SetBranchAddress("MatchTrigChi2_rec", MatchTrigChi2_rec);
  tree->SetBranchAddress("Charge_rec", Charge_rec);
  tree->SetBranchAddress("RAtAbsEnd_rec", RAtAbsEnd_rec);
  tree->SetBranchAddress("pDCA_rec", pDCA_rec);
  tree->SetBranchAddress("Phi_rec", Phi_rec);
  tree->SetBranchAddress("Theta_rec", Theta_rec);

  if (from_data)
    printf("Data option not avaible\n");

  return tree;
}

TH1D *Get1Dhist(const char *filename, const char *hist_name, bool norm)
{

  TFile *myFile = TFile::Open(filename, "UPDATE");

  if (!myFile || myFile->IsZombie())
  {
    printf("%s non trovato \n", filename);
    //        cout<<"File non trovato"<<endl;
    return nullptr;
  }
  TH1D *h_hist = (TH1D *)myFile->Get(hist_name);
  if (!h_hist || h_hist->IsZombie())
  {
    printf("%s non trovato \n", hist_name);
    //        cout<<"Istogramma non presente o con nome errato"<<endl;
    return nullptr;
  }
  h_hist->SetDirectory(0);
  h_hist->SetMinimum(1);
  if (norm)
    h_hist->Scale(1, "width");
  myFile->Close();
  //    gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  h_hist->SetMarkerStyle(0);
  return h_hist;
}

void HistMCDimuonHF_lite(
    const char *RunMode = "HF",
    Int_t RunNumber = 294009,
    TString prefix_dir_filename = "30_09_2022",
    TString dir_fileout = "/home/michele_pennisi/dimuon_HF_pp/MC_Analysis/hist_mc_aod/test",
    TString prefix_filename = "MCDimuHFTree")
{
  Double_t Mass_cut = 0;
  SetHist();
  TString dir_filename;
  // dir_filename.Form("/home/michele_pennisi/dimuon_HF_pp/Grid_Sim/%s/RunMCDimuHF/%s", RunMode, prefix_dir_filename.Data());//Local test
  dir_filename.Form("/media/michele_pennisi/DataBoi/Grid_Sim/%s/RunMCDimuHF/%s", RunMode, prefix_dir_filename.Data()); // official
  printf("%s", dir_filename.Data());

  TString filename;
  filename.Form("%s_MCDimuHFTree_%d.root", RunMode, RunNumber);

  printf("Creating hist for %s/%s, with Mass Cut = %0.1f\n", dir_filename.Data(), filename.Data(), Mass_cut);
  TString fileout;
  fileout.Form("%s/%s/HistLite_%s", dir_fileout.Data(), RunMode, filename.Data());
  // fileout.Form("%sLite_Hist%s", RunMode, filename.Data()); // test
  // fileout.Form("Test_hist.root"); // test GRid

  TString fileout_tree;
  fileout_tree.Form("%s/%s/Tree_%s", dir_fileout.Data(), RunMode, filename.Data());
  // fileout_tree.Form("%sTree%s", RunMode, filename.Data()); // test
  // fileout_tree.Form("Test_tree.root"); // test Grid

  TChain *input_tree = Getting_Tree(dir_filename, filename, false);

  TFile fr(fileout_tree.Data(), "recreate");
  TTree *rec_tree_mucharm = new TTree("rec_tree_mucharm", "Tree pt dimuon charm with 4<M<30");
  TTree *rec_tree_mubeauty = new TTree("rec_tree_mubeauty", "Tree pt dimuon beauty with 4<M<30");
  TTree *rec_tree_mumixed = new TTree("rec_tree_mumixed", "Tree pt dimuon mixed with 4<M<30");

  TTree *rec_tree_mucharm_withcut = new TTree("rec_tree_mucharm_withcut", "Tree m dimuon charm with 4<M<30");
  TTree *rec_tree_mubeauty_withcut = new TTree("rec_tree_mubeauty_withcut", "Tree m dimuon beauty with 4<M<30");
  TTree *rec_tree_mumixed_withcut = new TTree("rec_tree_mumixed_withcut", "Tree m dimuon mixed with 4<M<30");

  TTree *rec_tree_mucharm_lowmass = new TTree("rec_tree_mucharm_LowMass", "Tree m dimuon charm with 4<M<9");
  TTree *rec_tree_mubeauty_lowmass = new TTree("rec_tree_mubeauty_LowMass", "Tree m dimuon beauty with 4<M<9");
  TTree *rec_tree_mumixed_lowmass = new TTree("rec_tree_mumixed_LowMass", "Tree m dimuon mixed with 4<M<9");

  TTree *rec_tree_mucharm_highmass = new TTree("rec_tree_mucharm_HighMass", "Tree m dimuon charm with 11<M<30");
  TTree *rec_tree_mubeauty_highmass = new TTree("rec_tree_mubeauty_HighMass", "Tree m dimuon beauty with 11<M<30");
  TTree *rec_tree_mumixed_highmass = new TTree("rec_tree_mumixed_HighMass", "Tree m dimuon mixed with 11<M<30");

  TTree *rec_tree_mucharm_lowmass_lowpt = new TTree("rec_tree_mucharm_LowMass_LowPt", "Tree m dimuon charm with 4<M<9 && pt<10");
  TTree *rec_tree_mubeauty_lowmass_lowpt = new TTree("rec_tree_mubeauty_LowMass_LowPt", "Tree m dimuon beauty with 4<M<9 && pt<10");
  TTree *rec_tree_mumixed_lowmass_lowpt = new TTree("rec_tree_mumixed_LowMass_LowPt", "Tree m dimuon mixed with 4<M<9 && pt<10");

  Double_t *pt_dimu_charm = new Double_t;
  Double_t *pt_dimu_beauty = new Double_t;
  Double_t *pt_dimu_mixed = new Double_t;

  Double_t *pt_dimu_charm_withcut = new Double_t;
  Double_t *pt_dimu_beauty_withcut = new Double_t;
  Double_t *pt_dimu_mixed_withcut = new Double_t;

  Double_t *pt_dimu_charm_lowmass = new Double_t;
  Double_t *pt_dimu_beauty_lowmass = new Double_t;
  Double_t *pt_dimu_mixed_lowmass = new Double_t;

  Double_t *pt_dimu_charm_highmass = new Double_t;
  Double_t *pt_dimu_beauty_highmass = new Double_t;
  Double_t *pt_dimu_mixed_highmass = new Double_t;

  Double_t *pt_dimu_charm_lowmass_lowpt = new Double_t;
  Double_t *pt_dimu_beauty_lowmass_lowpt = new Double_t;
  Double_t *pt_dimu_mixed_lowmass_lowpt = new Double_t;

  Double_t *m_dimu_charm = new Double_t;
  Double_t *m_dimu_beauty = new Double_t;
  Double_t *m_dimu_mixed = new Double_t;

  Double_t *m_dimu_charm_withcut = new Double_t;
  Double_t *m_dimu_beauty_withcut = new Double_t;
  Double_t *m_dimu_mixed_withcut = new Double_t;

  Double_t *m_dimu_charm_lowmass = new Double_t;
  Double_t *m_dimu_beauty_lowmass = new Double_t;
  Double_t *m_dimu_mixed_lowmass = new Double_t;

  Double_t *m_dimu_charm_highmass = new Double_t;
  Double_t *m_dimu_beauty_highmass = new Double_t;
  Double_t *m_dimu_mixed_highmass = new Double_t;

  Double_t *m_dimu_charm_lowmass_lowpt = new Double_t;
  Double_t *m_dimu_beauty_lowmass_lowpt = new Double_t;
  Double_t *m_dimu_mixed_lowmass_lowpt = new Double_t;

  rec_tree_mucharm->Branch("pt", pt_dimu_charm, "pt/D");
  rec_tree_mucharm->Branch("m", m_dimu_charm, "m/D");

  rec_tree_mubeauty->Branch("pt", pt_dimu_beauty, "pt/D");
  rec_tree_mubeauty->Branch("m", m_dimu_beauty, "m/D");

  rec_tree_mumixed->Branch("pt", pt_dimu_mixed, "pt/D");
  rec_tree_mumixed->Branch("m", m_dimu_mixed, "m/D");

  rec_tree_mucharm_withcut->Branch("pt", pt_dimu_charm_withcut, "pt/D");
  rec_tree_mucharm_withcut->Branch("m", m_dimu_charm_withcut, "m/D");

  rec_tree_mubeauty_withcut->Branch("pt", pt_dimu_beauty_withcut, "pt/D");
  rec_tree_mubeauty_withcut->Branch("m", m_dimu_beauty_withcut, "m/D");

  rec_tree_mumixed_withcut->Branch("pt", pt_dimu_mixed_withcut, "pt/D");
  rec_tree_mumixed_withcut->Branch("m", m_dimu_mixed_withcut, "m/D");

  rec_tree_mucharm_lowmass->Branch("pt", pt_dimu_charm_lowmass, "pt/D");
  rec_tree_mucharm_lowmass->Branch("m", m_dimu_charm_lowmass, "m/D");

  rec_tree_mubeauty_lowmass->Branch("pt", pt_dimu_beauty_lowmass, "pt/D");
  rec_tree_mubeauty_lowmass->Branch("m", m_dimu_beauty_lowmass, "m/D");

  rec_tree_mumixed_lowmass->Branch("pt", pt_dimu_mixed_lowmass, "pt/D");
  rec_tree_mumixed_lowmass->Branch("m", m_dimu_mixed_lowmass, "m/D");

  rec_tree_mucharm_highmass->Branch("pt", pt_dimu_charm_highmass, "pt/D");
  rec_tree_mucharm_highmass->Branch("m", m_dimu_charm_highmass, "m/D");

  rec_tree_mubeauty_highmass->Branch("pt", pt_dimu_beauty_highmass, "pt/D");
  rec_tree_mubeauty_highmass->Branch("m", m_dimu_beauty_highmass, "m/D");

  rec_tree_mumixed_highmass->Branch("pt", pt_dimu_mixed_highmass, "pt/D");
  rec_tree_mumixed_highmass->Branch("m", m_dimu_mixed_highmass, "m/D");

  rec_tree_mucharm_lowmass_lowpt->Branch("pt", pt_dimu_charm_lowmass_lowpt, "pt/D");
  rec_tree_mucharm_lowmass_lowpt->Branch("m", m_dimu_charm_lowmass_lowpt, "m/D");

  rec_tree_mubeauty_lowmass_lowpt->Branch("pt", pt_dimu_beauty_lowmass_lowpt, "pt/D");
  rec_tree_mubeauty_lowmass_lowpt->Branch("m", m_dimu_beauty_lowmass_lowpt, "m/D");

  rec_tree_mumixed_lowmass_lowpt->Branch("pt", pt_dimu_mixed_lowmass_lowpt, "pt/D");
  rec_tree_mumixed_lowmass_lowpt->Branch("m", m_dimu_mixed_lowmass_lowpt, "m/D");

  if (!input_tree)
    return;

  input_tree->ls();

  for (Int_t i = 0; i < input_tree->GetEntries(); i++)
  {
    Int_t Dimuon_rec_charm = 0;
    Int_t Dimuon_rec_beauty = 0;
    Int_t Dimuon_rec_mixed = 0;

    Int_t counter_charm = 0;
    Int_t counter_beauty = 0;
    Int_t counter_mixed = 0;

    Int_t N_charm_xevent = 0;
    Int_t N_beauty_xevent = 0;

    Int_t N_anticharm_xevent = 0;
    Int_t N_antibeauty_xevent = 0;

    Int_t N_charmpair_xevent = 0;
    Int_t N_beautypair_xevent = 0;

    Double_t Pt_dimu_Charm[100];
    Double_t M_dimu_Charm[100];

    Double_t Pt_dimu_Beauty[100];
    Double_t M_dimu_Beauty[100];

    Double_t Pt_dimu_Mixed[100];
    Double_t M_dimu_Mixed[100];

    input_tree->GetEntry(i);
    if (i % 500000 == 0)
      printf("Evento :i %d\n", i);

    for (Int_t i = 0; i < N_HFquarks_gen; i++)
    {
      Int_t PDG_quark_gen1 = PDG_HFquark_gen[i];  // single gen mu PDG mum
      Double_t Pt_quark_gen1 = Pt_HFquark_gen[i]; // single gen mu pT
      Double_t Y_quark_gen1 = Y_HFquark_gen[i];   // single gen mu E

      if (PDG_quark_gen1 == 4)
      {
        h_pt_charm_quark->Fill(Pt_quark_gen1);
        h_y_charm_quark->Fill(Y_quark_gen1);
        N_charm_xevent++;
      }
      else if (PDG_quark_gen1 == -4)
      {
        h_pt_charm_antiquark->Fill(Pt_quark_gen1);
        h_y_charm_antiquark->Fill(Y_quark_gen1);
        N_anticharm_xevent++;
      }
      else if (PDG_quark_gen1 == 5)
      {
        h_pt_beauty_quark->Fill(Pt_quark_gen1);
        h_y_beauty_quark->Fill(Y_quark_gen1);
        N_beauty_xevent++;
      }
      else if (PDG_quark_gen1 == 5)
      {
        h_pt_beauty_antiquark->Fill(Pt_quark_gen1);
        h_y_beauty_antiquark->Fill(Y_quark_gen1);
        N_antibeauty_xevent++;
      }

      for (Int_t j = i + 1; i < N_HFquarks_gen; i++)
      {
        Int_t PDG_quark_gen2 = PDG_HFquark_gen[j];  // sjngle gen mu PDG mum
        Double_t Pt_quark_gen2 = Pt_HFquark_gen[j]; // single gen mu pT
        Double_t Y_quark_gen2 = Y_HFquark_gen[j];   // single gen mu E

        if (TMath::Abs(PDG_quark_gen1) == 4 && TMath::Abs(PDG_quark_gen2) == 4)
        {
          if ((PDG_quark_gen1 + PDG_quark_gen2) == 0)
          {
            N_charmpair_xevent++;
          }
        }

        if (TMath::Abs(PDG_quark_gen1) == 5 && TMath::Abs(PDG_quark_gen2) == 5)
        {
          if ((PDG_quark_gen1 + PDG_quark_gen2) == 0)
          {
            N_beautypair_xevent++;
          }
        }
      }
    }
    h_Ncharm_pairs->Fill(N_charmpair_xevent);
    h_Nbeauty_pairs->Fill(N_beautypair_xevent);

    h_Ncharm_quark->Fill(N_charm_xevent);
    h_Nbeauty_quark->Fill(N_beauty_xevent);
    /*
      for (Int_t i = 0; i < NMuons_gen; i++)
      {
        Int_t PDG_Mu = PDGmum_gen[i];     // single gen mu PDG mum
        Double_t Pt_Mu = Pt_gen[i];       // single gen mu pT
        Double_t E_Mu = E_gen[i];         // single gen mu E
        Double_t Px_Mu = Px_gen[i];       // single gen mu px
        Double_t Py_Mu = Py_gen[i];       // single gen mu py
        Double_t Pz_Mu = Pz_gen[i];       // single gen mu pz
        Double_t Y_Mu = Y_gen[i];         // single gen mu y
        Double_t Eta_Mu = Eta_gen[i];     // single gen mu eta
        Double_t Phi_Mu = Phi_gen[i];     // single gen mu phi
        Double_t Theta_Mu = Theta_gen[i]; // single gen mu theta

        Bool_t Selection_Mu_gen[n_MuSelection] = {kFALSE};
        Bool_t Kin_Mu_gen[n_MuCut] = {kFALSE};

        Selection_Mu_gen[0] = kTRUE;
        Kin_Mu_gen[0] = kTRUE;

        if ((TMath::Abs(PDG_Mu) == 4) || (TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500) || (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000))
        {
          Selection_Mu_gen[1] = kTRUE;
          Selection_Mu_gen[2] = kTRUE;
        }
        else if ((TMath::Abs(PDG_Mu) == 5) || (TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600) || (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000))
        {
          Selection_Mu_gen[1] = kTRUE;
          Selection_Mu_gen[3] = kTRUE;
        }
        else if ((TMath::Abs(PDG_Mu) > 0 && (TMath::Abs(PDG_Mu) < 4)) || (TMath::Abs(PDG_Mu) > 100 && TMath::Abs(PDG_Mu) < 400) || (TMath::Abs(PDG_Mu) > 1000 && TMath::Abs(PDG_Mu) < 4000))
        {
          Selection_Mu_gen[4] = kTRUE;
        }
        else
        {
          Selection_Mu_gen[5] = kTRUE;
        }

        if ((Y_Mu > -4.0 && Y_Mu < -2.5))
        {
          Kin_Mu_gen[1] = kTRUE;

          if (Pt_Mu > 0.5)
          {
            Kin_Mu_gen[2] = kTRUE;
          }
          else
            Kin_Mu_gen[2] = kFALSE;
        }

        if ((Eta_Mu > -4.0 && Eta_Mu < -2.5))
        {
          Kin_Mu_gen[3] = kTRUE;
        }

        for (Int_t Mu_b = 0; Mu_b < n_MuCut - n_MuCut_rec; Mu_b++)
        {
          if (Kin_Mu_gen[Mu_b])
          {
            for (Int_t Mu_c = 0; Mu_c < n_MuSelection; Mu_c++)
            {
              if (Selection_Mu_gen[Mu_c])
              {
                h_PtYMu[0][Mu_b][Mu_c]->Fill(Pt_Mu, Y_Mu);
                h_pdgMu[0][Mu_b][Mu_c]->Fill(PDG_Mu);
                nMu_xevent[0][Mu_b][Mu_c]++;
              }
            } // End loop on muon selection
          }
        } // End loop on kinematic cuts

      } // End loop on generated muons

      for (Int_t i = 0; i < NMuons_rec; i++)
      {
        Int_t counter_rec = n_MuCut - n_MuCut_rec;

        Int_t PDG_Mu = PDGmum_rec[i];             // single rec mu PDG mum
        Double_t Pt_Mu = Pt_rec[i];               // single rec mu pT
        Double_t E_Mu = E_rec[i];                 // single rec mu E
        Double_t Px_Mu = Px_rec[i];               // single rec mu px
        Double_t Py_Mu = Py_rec[i];               // single rec mu py
        Double_t Pz_Mu = Pz_rec[i];               // single rec mu pz
        Double_t Y_Mu = Y_rec[i];                 // single rec mu y
        Double_t Eta_Mu = Eta_rec[i];             // single rec mu eta
        Double_t Phi_Mu = Phi_rec[i];             // single rec mu phi
        Double_t Theta_Mu = Theta_rec[i];         // single rec mu theta
        Double_t MatchTrig_Mu = MatchTrig_rec[i]; // single gen mu theta
        Double_t RAbs_Mu = RAtAbsEnd_rec[i];
        Double_t pDCA_Mu = pDCA_rec[i];

        Bool_t Selection_Mu_rec[n_MuSelection] = {kFALSE};
        Bool_t Kin_Mu_rec[n_DiMuCut] = {kFALSE};

        Selection_Mu_rec[0] = kTRUE;
        Kin_Mu_rec[counter_rec] = kTRUE;

        if ((TMath::Abs(PDG_Mu) == 4) || (TMath::Abs(PDG_Mu) > 400 && TMath::Abs(PDG_Mu) < 500) || (TMath::Abs(PDG_Mu) > 4000 && TMath::Abs(PDG_Mu) < 5000))
        {
          Selection_Mu_rec[1] = kTRUE;
          Selection_Mu_rec[2] = kTRUE;
        }
        else if ((TMath::Abs(PDG_Mu) == 5) || (TMath::Abs(PDG_Mu) > 500 && TMath::Abs(PDG_Mu) < 600) || (TMath::Abs(PDG_Mu) > 5000 && TMath::Abs(PDG_Mu) < 6000))
        {
          Selection_Mu_rec[1] = kTRUE;
          Selection_Mu_rec[3] = kTRUE;
        }
        else if ((TMath::Abs(PDG_Mu) > 0 && (TMath::Abs(PDG_Mu) < 4)) || (TMath::Abs(PDG_Mu) > 100 && TMath::Abs(PDG_Mu) < 400) || (TMath::Abs(PDG_Mu) > 1000 && TMath::Abs(PDG_Mu) < 4000))
        {
          Selection_Mu_rec[4] = kTRUE;
        }
        else
        {
          Selection_Mu_rec[5] = kTRUE;
        }
        // printf("%0.0f \n",MatchTrig_Mu);
        if ((RAbs_Mu > 17.6 && RAbs_Mu < 89.5) && pDCA_Mu == 1)
        {
          if ((Pt_Mu > 0.0 && Pt_Mu < 30.0))
          {

            if ((Y_Mu > -4.0 && Y_Mu < -2.5))

              if (MatchTrig_Mu >= 0)
              {

                Kin_Mu_rec[counter_rec + 1] = kTRUE;

                if (MatchTrig_Mu >= 2)
                {
                  Kin_Mu_rec[counter_rec + 2] = kTRUE;
                }
              }
            if ((Eta_Mu > -4.0 && Eta_Mu < -2.5))
            {
              if (MatchTrig_Mu >= 0)
              {

                Kin_Mu_rec[counter_rec + 3] = kTRUE;

                if (MatchTrig_Mu >= 2)
                {
                  Kin_Mu_rec[counter_rec + 4] = kTRUE;
                }
              }
            }
          }
        } // End if to apply the kinematic cut applied in the data

        for (Int_t Mu_b = counter_rec; Mu_b < n_MuCut; Mu_b++)
        {
          for (Int_t Mu_c = 0; Mu_c < n_MuSelection; Mu_c++)
          {
            if (Selection_Mu_rec[Mu_c] && Kin_Mu_rec[Mu_b])
            {
              h_PtYMu[0][Mu_b][Mu_c]->Fill(Pt_Mu, Y_Mu);
              h_pdgMu[0][Mu_b][Mu_c]->Fill(PDG_Mu);
              nMu_xevent[0][Mu_b][Mu_c]++;
            }
          } // End loop on muon selection
        }

      } // End loop on reconstructed muons
  */
    //----------------------Start Generated Dimuon---------------------------------//

    // for (Int_t i = 0; i < NDimu_gen; i++)
    // {

    //   Double_t Pt_DiMu = DimuPt_gen[i];      // gen dimuon pT
    //   Double_t Px_DiMu = DimuPx_gen[i];      // gen dimuon px
    //   Double_t Py_DiMu = DimuPy_gen[i];      // gen dimuon py
    //   Double_t Pz_DiMu = DimuPz_gen[i];      // gen dimuon pz
    //   Double_t Y_DiMu = DimuY_gen[i];        // gen dimuon y
    //   Double_t M_DiMu = DimuMass_gen[i];     // gen dimuon invariant mass
    //   Int_t Charge_DiMu = DimuCharge_gen[i]; // gen dimuon charge

    //   Double_t Pt_Mu0 = Pt_gen[DimuMu_gen[i][0]];
    //   Double_t Y_Mu0 = Y_gen[DimuMu_gen[i][0]];
    //   Double_t Eta_Mu0 = Eta_gen[DimuMu_gen[i][0]];

    //   Double_t Pt_Mu1 = Pt_gen[DimuMu_gen[i][1]];
    //   Double_t Y_Mu1 = Y_gen[DimuMu_gen[i][1]];
    //   Double_t Eta_Mu1 = Eta_gen[DimuMu_gen[i][1]];

    //   Int_t PDG_Mu0 = PDGmum_gen[DimuMu_gen[i][0]];
    //   Int_t PDG_Mu1 = PDGmum_gen[DimuMu_gen[i][1]];

    //   Bool_t Mass_DiMu_gen[n_Mass_Cut] = {kFALSE};
    //   /*
    //   0 For no cut on Dimu Mass
    //   1 For M>4 GeV/c^2
    //   */
    //   Bool_t Kin_DiMu_gen[n_DiMuCut] = {kFALSE};
    //   /*
    //   0 For Generated
    //   1 For Generated -4.0<Y<-2.5
    //   2 For Generated in ALICEacc
    //   4 For Reconstructed
    //   5 For Reconstructed with all cut applied
    //   */
    //   Bool_t Charge_DiMu_gen[n_DiMu_Charge] = {kFALSE};
    //   /*
    //   0 for ULS,
    //   1 for LS++,
    //   2 for LS--,
    //   3 LS total
    //   */

    //   Bool_t Selection_DiMu_gen[n_DiMuSelection] = {kFALSE};
    //   /*
    //   0 For All
    //   1 For HF
    //   2 For Charm
    //   3 For Beauty
    //   4 For HF Mixed (one muon from Charm, one muon from Beauty)
    //   5 For LF (two muons from LF)
    //   6 For LF Mixed (one muon from LH, one muon from HF)
    //   7 Others
    //   */

    //   Selection_DiMu_gen[0] = kTRUE;
    //   Kin_DiMu_gen[0] = kTRUE;
    //   if (M_DiMu > 0.0)
    //     Mass_DiMu_gen[0] = kTRUE;

    //   Bool_t Charm_mu0 = kFALSE;
    //   Bool_t Charm_mu1 = kFALSE;
    //   Bool_t Beauty_mu0 = kFALSE;
    //   Bool_t Beauty_mu1 = kFALSE;
    //   Bool_t HF_mu0 = kFALSE;
    //   Bool_t HF_mu1 = kFALSE;
    //   Bool_t LF_mu0 = kFALSE;
    //   Bool_t LF_mu1 = kFALSE;

    //   if (M_DiMu > 4.0 && M_DiMu <= 30.0)
    //     Mass_DiMu_gen[1] = kTRUE;

    //   if (M_DiMu > 4.0 && M_DiMu <= 9.0)
    //   {
    //     Mass_DiMu_gen[2] = kTRUE;
    //     Mass_DiMu_gen[6] = kTRUE;
    //     if (Pt_DiMu<=10)
    //     {
    //       Mass_DiMu_gen[8] = kTRUE;
    //     }

    //   }
    //   else if (M_DiMu > 9.0 && M_DiMu <= 11.0)
    //   {
    //     Mass_DiMu_gen[3] = kTRUE;
    //   }
    //   else if (M_DiMu > 11.0 && M_DiMu <= 15.0)
    //   {
    //     Mass_DiMu_gen[4] = kTRUE;
    //     Mass_DiMu_gen[6] = kTRUE;
    //     Mass_DiMu_gen[7] = kTRUE;
    //   }
    //   else if (M_DiMu > 15.0 && M_DiMu <= 30.0)
    //   {
    //     Mass_DiMu_gen[5] = kTRUE;
    //     Mass_DiMu_gen[6] = kTRUE;
    //     Mass_DiMu_gen[7] = kTRUE;
    //   }

    //   if ((TMath::Abs(PDG_Mu0) == 4) || (TMath::Abs(PDG_Mu0) > 400 && TMath::Abs(PDG_Mu0) < 500) || (TMath::Abs(PDG_Mu0) > 4000 && TMath::Abs(PDG_Mu0) < 5000))
    //   {
    //     HF_mu0 = kTRUE;
    //     Charm_mu0 = kTRUE;
    //   }
    //   else if ((TMath::Abs(PDG_Mu0) == 5) || (TMath::Abs(PDG_Mu0) > 500 && TMath::Abs(PDG_Mu0) < 600) || (TMath::Abs(PDG_Mu0) > 5000 && TMath::Abs(PDG_Mu0) < 6000))
    //   {
    //     HF_mu0 = kTRUE;
    //     Beauty_mu0 = kTRUE;
    //   }
    //   else if ((TMath::Abs(PDG_Mu0) > 0 && TMath::Abs(PDG_Mu0) < 4) || (TMath::Abs(PDG_Mu0) > 100 && TMath::Abs(PDG_Mu0) < 400) || (TMath::Abs(PDG_Mu0) > 1000 && TMath::Abs(PDG_Mu0) < 4000))
    //   {
    //     LF_mu0 = kTRUE;
    //     Charm_mu0 = kFALSE;
    //     Beauty_mu0 = kFALSE;
    //     HF_mu0 = kFALSE;
    //   }
    //   else
    //   {
    //     LF_mu0 = kFALSE;
    //     Charm_mu0 = kFALSE;
    //     Beauty_mu0 = kFALSE;
    //     HF_mu0 = kFALSE;
    //   }

    //   if ((TMath::Abs(PDG_Mu1) == 4) || (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500) || (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000))
    //   {
    //     HF_mu1 = kTRUE;
    //     Charm_mu1 = kTRUE;
    //   }
    //   else if ((TMath::Abs(PDG_Mu1) == 5) || (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600) || (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000))
    //   {
    //     HF_mu1 = kTRUE;
    //     Beauty_mu1 = kTRUE;
    //   }
    //   else if ((TMath::Abs(PDG_Mu1) > 0 && TMath::Abs(PDG_Mu1) < 4) || (TMath::Abs(PDG_Mu1) > 100 && TMath::Abs(PDG_Mu1) < 400) || (TMath::Abs(PDG_Mu1) > 1000 && TMath::Abs(PDG_Mu1) < 4000))
    //   {
    //     LF_mu1 = kTRUE;
    //     Charm_mu1 = kFALSE;
    //     Beauty_mu1 = kFALSE;
    //     HF_mu1 = kFALSE;
    //   }
    //   else
    //   {
    //     LF_mu1 = kFALSE;
    //     Charm_mu1 = kFALSE;
    //     Beauty_mu1 = kFALSE;
    //     HF_mu1 = kFALSE;
    //   }

    //   if (HF_mu0 && HF_mu1)
    //   {
    //     Selection_DiMu_gen[1] = kTRUE;
    //     if (Charm_mu0 && Charm_mu1)
    //       Selection_DiMu_gen[2] = kTRUE;
    //     else if (Beauty_mu0 && Beauty_mu1)
    //       Selection_DiMu_gen[3] = kTRUE;
    //     else if ((Charm_mu0 && Beauty_mu1) || (Beauty_mu0 && Charm_mu1))
    //       Selection_DiMu_gen[4] = kTRUE;
    //   }
    //   else if (LF_mu0 && LF_mu1)
    //     Selection_DiMu_gen[5] = kTRUE;
    //   else if ((LF_mu0 && HF_mu1) || (HF_mu0 && LF_mu1))
    //     Selection_DiMu_gen[6] = kTRUE;
    //   else
    //     Selection_DiMu_gen[7] = kTRUE;

    //   if ((Y_Mu0 > -4.0 && Y_Mu0 < -2.5) && (Y_Mu1 > -4.0 && Y_Mu1 < -2.5))
    //   {
    //     Kin_DiMu_gen[1] = kTRUE;

    //     if ((Pt_Mu0 > 0.5) && (Pt_Mu1 > 0.5))
    //       Kin_DiMu_gen[2] = kTRUE;
    //   }

    //   if ((Eta_Mu0 > -4.0 && Eta_Mu0 < -2.5) && (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5))
    //   {
    //     Kin_DiMu_gen[3] = kTRUE;
    //   }

    //   if (Charge_DiMu == 0)
    //   {
    //     Charge_DiMu_gen[0] = kTRUE;
    //     Charge_DiMu_gen[1] = kFALSE;
    //     Charge_DiMu_gen[2] = kFALSE;
    //     Charge_DiMu_gen[3] = kFALSE;
    //   }
    //   else if (Charge_DiMu == 2)
    //   {
    //     Charge_DiMu_gen[0] = kFALSE;
    //     Charge_DiMu_gen[1] = kTRUE;
    //     Charge_DiMu_gen[2] = kFALSE;
    //     Charge_DiMu_gen[3] = kTRUE;
    //   }
    //   else if (Charge_DiMu == -2)
    //   {
    //     Charge_DiMu_gen[0] = kFALSE;
    //     Charge_DiMu_gen[1] = kFALSE;
    //     Charge_DiMu_gen[2] = kTRUE;
    //     Charge_DiMu_gen[3] = kTRUE;
    //   }

    //   for (Int_t DiMu_a = 0; DiMu_a < n_Mass_Cut; DiMu_a++)
    //   {
    //     for (Int_t DiMu_b = 0; DiMu_b < n_DiMuCut - n_DiMuCut_rec; DiMu_b++)
    //     {
    //       for (Int_t DiMu_c = 0; DiMu_c < n_DiMu_Charge; DiMu_c++)
    //       {
    //         for (Int_t DiMu_d = 0; DiMu_d < n_DiMuSelection; DiMu_d++)
    //         {

    //           if (Mass_DiMu_gen[DiMu_a] && Kin_DiMu_gen[DiMu_b] && Charge_DiMu_gen[DiMu_c] && Selection_DiMu_gen[DiMu_d])
    //           {
    //             // printf("Filling %s with Charge %d\n",name_DiMu_Charge[DiMu_c].Data(),Charge_DiMu);
    //             h_PtYDiMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(Pt_DiMu, Y_DiMu);
    //             h_PtMDiMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(Pt_DiMu, M_DiMu);
    //             h_MDiMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(M_DiMu);
    //             h_pdgDimuMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(PDG_Mu0);
    //             h_pdgDimuMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(PDG_Mu1);
    //             nDiMu_xevent[DiMu_a][DiMu_b][DiMu_c][DiMu_d]++;
    //           }

    //         } // End definition over DiMu selection

    //       } // End definition over DiMu charge cut

    //     } // End definition over DiMu cut

    //   } // End loop over M_DiMu cut

    // } // End loop on generated Dimuons

    //----------------------Start Reconstructed Dimuon---------------------------------//

    for (Int_t i = 0; i < NDimu_rec; i++)
    {

      Int_t counter_rec = n_DiMuCut - n_DiMuCut_rec;

      Double_t Pt_DiMu = DimuPt_rec[i];      // rec dimuon pT
      Double_t Px_DiMu = DimuPx_rec[i];      // rec dimuon px
      Double_t Py_DiMu = DimuPy_rec[i];      // rec dimuon py
      Double_t Pz_DiMu = DimuPz_rec[i];      // rec dimuon pz
      Double_t Y_DiMu = DimuY_rec[i];        // rec dimuon y
      Double_t M_DiMu = DimuMass_rec[i];     // rec dimuon invariant mass
      Int_t Charge_DiMu = DimuCharge_rec[i]; // rec dimuon charge

      Double_t Pt_Mu0 = Pt_rec[DimuMu_rec[i][0]];
      Double_t Y_Mu0 = Y_rec[DimuMu_rec[i][0]];
      Int_t PDG_Mu0 = PDGmum_rec[DimuMu_rec[i][0]];
      Double_t Charge_Mu0 = Charge_rec[DimuMu_rec[i][0]];
      Double_t RAbs_Mu0 = RAtAbsEnd_rec[DimuMu_rec[i][0]];
      Double_t pDCA_Mu0 = pDCA_rec[DimuMu_rec[i][0]];
      Double_t Eta_Mu0 = Eta_rec[DimuMu_rec[i][0]];

      Double_t Y_Mu1 = Y_rec[DimuMu_rec[i][1]];
      Double_t Pt_Mu1 = Pt_rec[DimuMu_rec[i][1]];
      Int_t PDG_Mu1 = PDGmum_rec[DimuMu_rec[i][1]];
      Double_t Charge_Mu1 = Charge_rec[DimuMu_rec[i][1]];
      Double_t RAbs_Mu1 = RAtAbsEnd_rec[DimuMu_rec[i][1]];
      Double_t pDCA_Mu1 = pDCA_rec[DimuMu_rec[i][1]];
      Double_t Eta_Mu1 = Eta_rec[DimuMu_rec[i][1]];

      // Int_t controlcharge0 = TMath::Sign(1, PDG_Mu0);
      // Int_t controlcharge1 = TMath::Sign(1, PDG_Mu1);
      // printf("Charge DiMu %d Charg1 %0.0f Charg2 %0.0f || Charge tot %d = CHarge1 %d + CHarge2 %d || PDG1 %d PDG2 %d \n", Charge_DiMu, Charge_Mu0, Charge_Mu1, controlcharge0 + controlcharge1, controlcharge0, controlcharge1, PDG_Mu0, PDG_Mu1);

      Bool_t Mass_DiMu_rec[n_Mass_Cut] = {kFALSE};
      /*
      0 For no cut on Dimu Mass
      1 For M>4 GeV/c^2
      */
      Bool_t Kin_DiMu_rec[n_DiMuCut] = {kFALSE};
      /*
      0 For Generated,
      1 For Generated -4.0<Y<-2.5
      2 For Generated in ALICEacc
      3 For Generated with DQ Cut
      4 For Reconstructed
      5 For Reconstructed with OLD all cut applied
      6 For Reconstructed OLD with all cut,
      7 For Reconstructed with DQ CUT
      8 For Reconstructed with DQ CUT and match

      */
      Bool_t Charge_DiMu_rec[n_DiMu_Charge] = {kFALSE};
      /*
      0 for ULS,printf("Charge_Mu1 %d + Charge_Mu2 %d == %d || PDG1 %d PDG2 %d \n", Charge_Mu1, Charge_Mu2, Charge_Mu1 + Charge_Mu2, PDG_Mu1, PDG_Mu2);
      1 for LS++,
      2 for LS--,
      3 LS total
      */

      Bool_t Selection_DiMu_rec[n_DiMuSelection] = {kFALSE};
      /*
      0 For All
      1 For HF
      2 For Charm
      3 For Beauty
      4 For HF Mixed (one muon from Charm, one muon from Beauty)
      5 For LF (two muons from LF)
      6 For LF Mixed (one muon from LH, one muon from HF)
      7 Others
      */

      Selection_DiMu_rec[0] = kTRUE;
      Kin_DiMu_rec[counter_rec] = kTRUE;
      Mass_DiMu_rec[0] = kTRUE;

      Bool_t Charm_mu0 = kFALSE;
      Bool_t Charm_mu1 = kFALSE;
      Bool_t Beauty_mu0 = kFALSE;
      Bool_t Beauty_mu1 = kFALSE;
      Bool_t HF_mu0 = kFALSE;
      Bool_t HF_mu1 = kFALSE;
      Bool_t LF_mu0 = kFALSE;
      Bool_t LF_mu1 = kFALSE;

      if (M_DiMu > 4.0 && M_DiMu <= 30.0)
        Mass_DiMu_rec[1] = kTRUE;

      if (M_DiMu > 4.0 && M_DiMu <= 9.0)
      {
        Mass_DiMu_rec[2] = kTRUE;
        Mass_DiMu_rec[6] = kTRUE;
        if (Pt_DiMu <= 10.0)
        {
          Mass_DiMu_rec[8] = kTRUE;
        }
      }
      else if (M_DiMu > 9.0 && M_DiMu <= 11.0)
      {
        Mass_DiMu_rec[3] = kTRUE;
      }
      else if (M_DiMu > 11.0 && M_DiMu <= 15.0)
      {
        Mass_DiMu_rec[4] = kTRUE;
        Mass_DiMu_rec[6] = kTRUE;
        Mass_DiMu_rec[7] = kTRUE;
      }
      else if (M_DiMu > 15.0 && M_DiMu <= 30.0)
      {
        Mass_DiMu_rec[5] = kTRUE;
        Mass_DiMu_rec[6] = kTRUE;
        Mass_DiMu_rec[7] = kTRUE;
      }

      if ((TMath::Abs(PDG_Mu0) == 4) || (TMath::Abs(PDG_Mu0) > 400 && TMath::Abs(PDG_Mu0) < 500) || (TMath::Abs(PDG_Mu0) > 4000 && TMath::Abs(PDG_Mu0) < 5000))
      {
        HF_mu0 = kTRUE;
        Charm_mu0 = kTRUE;
      }
      else if ((TMath::Abs(PDG_Mu0) == 5) || (TMath::Abs(PDG_Mu0) > 500 && TMath::Abs(PDG_Mu0) < 600) || (TMath::Abs(PDG_Mu0) > 5000 && TMath::Abs(PDG_Mu0) < 6000))
      {
        HF_mu0 = kTRUE;
        Beauty_mu0 = kTRUE;
      }
      else if ((TMath::Abs(PDG_Mu0) > 0 && TMath::Abs(PDG_Mu0) < 4) || (TMath::Abs(PDG_Mu0) > 100 && TMath::Abs(PDG_Mu0) < 400) || (TMath::Abs(PDG_Mu0) > 1000 && TMath::Abs(PDG_Mu0) < 4000))
      {
        LF_mu0 = kTRUE;
        Charm_mu0 = kFALSE;
        Beauty_mu0 = kFALSE;
        HF_mu0 = kFALSE;
      }
      else
      {
        LF_mu0 = kFALSE;
        Charm_mu0 = kFALSE;
        Beauty_mu0 = kFALSE;
        HF_mu0 = kFALSE;
      }

      if ((TMath::Abs(PDG_Mu1) == 4) || (TMath::Abs(PDG_Mu1) > 400 && TMath::Abs(PDG_Mu1) < 500) || (TMath::Abs(PDG_Mu1) > 4000 && TMath::Abs(PDG_Mu1) < 5000))
      {
        HF_mu1 = kTRUE;
        Charm_mu1 = kTRUE;
      }
      else if ((TMath::Abs(PDG_Mu1) == 5) || (TMath::Abs(PDG_Mu1) > 500 && TMath::Abs(PDG_Mu1) < 600) || (TMath::Abs(PDG_Mu1) > 5000 && TMath::Abs(PDG_Mu1) < 6000))
      {
        HF_mu1 = kTRUE;
        Beauty_mu1 = kTRUE;
      }
      else if ((TMath::Abs(PDG_Mu1) > 0 && TMath::Abs(PDG_Mu1) < 4) || (TMath::Abs(PDG_Mu1) > 100 && TMath::Abs(PDG_Mu1) < 400) || (TMath::Abs(PDG_Mu1) > 1000 && TMath::Abs(PDG_Mu1) < 4000))
      {
        LF_mu1 = kTRUE;
        Charm_mu1 = kFALSE;
        Beauty_mu1 = kFALSE;
        HF_mu1 = kFALSE;
      }
      else
      {
        LF_mu1 = kFALSE;
        Charm_mu1 = kFALSE;
        Beauty_mu1 = kFALSE;
        HF_mu1 = kFALSE;
      }

      if (HF_mu0 && HF_mu1)
      {
        Selection_DiMu_rec[1] = kTRUE;
        if (Charm_mu0 && Charm_mu1)
          Selection_DiMu_rec[2] = kTRUE;
        else if (Beauty_mu0 && Beauty_mu1)
          Selection_DiMu_rec[3] = kTRUE;
        else if ((Charm_mu0 && Beauty_mu1) || (Beauty_mu0 && Charm_mu1))
          Selection_DiMu_rec[4] = kTRUE;
      }
      else if (LF_mu0 && LF_mu1)
        Selection_DiMu_rec[5] = kTRUE;
      else if ((LF_mu0 && HF_mu1) || (HF_mu0 && LF_mu1))
        Selection_DiMu_rec[6] = kTRUE;
      else
        Selection_DiMu_rec[7] = kTRUE;

      if (pDCA_Mu0 == 1 && pDCA_Mu1 == 1)
      {

        if ((RAbs_Mu0 > 17.6 && RAbs_Mu0 < 89.5) && (RAbs_Mu1 > 17.6 && RAbs_Mu1 < 89.5))
        {

          if ((Pt_Mu0 > 0.0 && Pt_Mu0 < 30.0) && (Pt_Mu1 > 0.0 && Pt_Mu1 < 30.0))
          {

            if ((Y_Mu0 > -4.0 && Y_Mu0 < -2.5) && (Y_Mu1 > -4.0 && Y_Mu1 < -2.5))
            {
              Kin_DiMu_rec[counter_rec + 1] = kTRUE;

              if (DimuMatch_rec[i] == 2)
              {
                Kin_DiMu_rec[counter_rec + 2] = kTRUE;
              }
            }
            if ((Y_DiMu > -4.0 && Y_DiMu < -2.5) && (Eta_Mu0 > -4.0 && Eta_Mu0 < -2.5) && (Eta_Mu1 > -4.0 && Eta_Mu1 < -2.5))
            {
              Kin_DiMu_rec[counter_rec + 3] = kTRUE;
              if (DimuMatch_rec[i] == 2)
              {
                Kin_DiMu_rec[counter_rec + 4] = kTRUE;
              }
            }
          }
        }
      } // End if to apply the kinematic cut applied in the data

      if (Charge_DiMu == 0)
      {
        Charge_DiMu_rec[0] = kTRUE;
      }
      else if (Charge_DiMu == 2)
      {
        Charge_DiMu_rec[1] = kTRUE;
        Charge_DiMu_rec[3] = kTRUE;
      }
      else if (Charge_DiMu == -2)
      {
        Charge_DiMu_rec[2] = kTRUE;
        Charge_DiMu_rec[3] = kTRUE;
      } // End selection over dimuons charge

      for (Int_t DiMu_a = 0; DiMu_a < n_Mass_Cut; DiMu_a++)
      {
        for (Int_t DiMu_b = n_DiMuCut - 1; DiMu_b < n_DiMuCut; DiMu_b++)
        {
          for (Int_t DiMu_c = 0; DiMu_c < n_DiMu_Charge; DiMu_c++)
          {
            for (Int_t DiMu_d = 0; DiMu_d < n_DiMuSelection; DiMu_d++)
            {

              if (Mass_DiMu_rec[DiMu_a] && Kin_DiMu_rec[DiMu_b] && Charge_DiMu_rec[DiMu_c] && Selection_DiMu_rec[DiMu_d])
              {
                h_PtYDiMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(Pt_DiMu, Y_DiMu);
                h_PtMDiMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(Pt_DiMu, M_DiMu);
                h_MDiMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(M_DiMu);
                h_pdgDimuMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(PDG_Mu0);
                h_pdgDimuMu[DiMu_a][DiMu_b][DiMu_c][DiMu_d]->Fill(PDG_Mu1);
                nDiMu_xevent[DiMu_a][DiMu_b][DiMu_c][DiMu_d]++;
              }

            } // End definition over DiMu selection

          } // End definition over DiMu charge cut

        } // End definition over DiMu cut

      } // End loop over M_DiMu cut

      if (Kin_DiMu_rec[8] && Charge_DiMu_rec[0])
      {
        /*
        2 For Charm
        3 For Beauty
        4 For HF Mixed (one muon from Charm, one muon from Beauty)
        */
        if (Selection_DiMu_rec[2])
        {
          // Pt_dimu_Charm[counter_charm] = Pt_DiMu;
          // M_dimu_Charm[counter_charm] = M_DiMu;
          // counter_charm++;
          if (Mass_DiMu_rec[1])
          {
            *pt_dimu_charm = Pt_DiMu;
            *m_dimu_charm = M_DiMu;
            rec_tree_mucharm->Fill();
            counter_charm += 2;
            Dimuon_rec_charm++;
          }

          if (Mass_DiMu_rec[2])
          {
            *pt_dimu_charm_lowmass = Pt_DiMu;
            *m_dimu_charm_lowmass = M_DiMu;
            rec_tree_mucharm_lowmass->Fill();
          }

          if (Mass_DiMu_rec[4] || Mass_DiMu_rec[5])
          {
            *pt_dimu_charm_highmass = Pt_DiMu;
            *m_dimu_charm_highmass = M_DiMu;
            rec_tree_mucharm_highmass->Fill();
          }

          if (Mass_DiMu_rec[6])
          {
            *pt_dimu_charm_withcut = Pt_DiMu;
            *m_dimu_charm_withcut = M_DiMu;
            rec_tree_mucharm_withcut->Fill();
          }

          if (Mass_DiMu_rec[8])
          {
            *pt_dimu_charm_lowmass_lowpt = Pt_DiMu;
            *m_dimu_charm_lowmass_lowpt = M_DiMu;
            rec_tree_mucharm_lowmass_lowpt->Fill();
          }
        }

        else if (Selection_DiMu_rec[3])
        {
          if (Mass_DiMu_rec[1])
          {

            *pt_dimu_beauty = Pt_DiMu;
            *m_dimu_beauty = M_DiMu;

            rec_tree_mubeauty->Fill();
            counter_beauty += 2;
            Dimuon_rec_beauty++;
          }
          if (Mass_DiMu_rec[2])
          {
            *pt_dimu_beauty_lowmass = Pt_DiMu;
            *m_dimu_beauty_lowmass = M_DiMu;
            rec_tree_mubeauty_lowmass->Fill();
          }

          if (Mass_DiMu_rec[4] || Mass_DiMu_rec[5])
          {
            *pt_dimu_beauty_highmass = Pt_DiMu;
            *m_dimu_beauty_highmass = M_DiMu;
            rec_tree_mubeauty_highmass->Fill();
          }

          if (Mass_DiMu_rec[6])
          {
            *pt_dimu_beauty_withcut = Pt_DiMu;
            *m_dimu_beauty_withcut = M_DiMu;

            rec_tree_mubeauty_withcut->Fill();
          }
          if (Mass_DiMu_rec[8])
          {
            *pt_dimu_beauty_lowmass_lowpt = Pt_DiMu;
            *m_dimu_beauty_lowmass_lowpt = M_DiMu;

            rec_tree_mubeauty_lowmass_lowpt->Fill();
          }
        }
        else if (Selection_DiMu_rec[4])
        {
          if (Mass_DiMu_rec[1])
          {
            *pt_dimu_mixed = Pt_DiMu;
            *m_dimu_mixed = M_DiMu;

            rec_tree_mumixed->Fill();
            counter_charm++;
            counter_beauty++;
            Dimuon_rec_mixed++;
          }
          if (Mass_DiMu_rec[2])
          {
            *pt_dimu_mixed_lowmass = Pt_DiMu;
            *m_dimu_mixed_lowmass = M_DiMu;
            rec_tree_mumixed_lowmass->Fill();
          }

          if (Mass_DiMu_rec[4] || Mass_DiMu_rec[5])
          {
            *pt_dimu_mixed_highmass = Pt_DiMu;
            *m_dimu_mixed_highmass = M_DiMu;
            rec_tree_mumixed_highmass->Fill();
          }

          if (Mass_DiMu_rec[6])
          {
            *pt_dimu_mixed_withcut = Pt_DiMu;
            *m_dimu_mixed_withcut = M_DiMu;

            rec_tree_mumixed_withcut->Fill();
          }

          if (Mass_DiMu_rec[8])
          {
            *pt_dimu_mixed_lowmass_lowpt = Pt_DiMu;
            *m_dimu_mixed_lowmass_lowpt = M_DiMu;

            rec_tree_mumixed_lowmass_lowpt->Fill();
          }
        }
      }

    } // End loop on reconstructed Dimuons
    h_nDiMu_xevent[1][8][0][2]->Fill(Dimuon_rec_charm);
    h_nDiMu_xevent[1][8][0][3]->Fill(Dimuon_rec_beauty);
    h_nDiMu_xevent[1][8][0][4]->Fill(Dimuon_rec_mixed);
    h_nMu_xevent_rec_charm->Fill(counter_charm);
    h_nMu_xevent_rec_beauty->Fill(counter_beauty);

  } // End Loop over Tree Entries

  fr.cd();
  h_Ncharm_quark->Write(0, 2, 0);
  h_pt_charm_quark->Write(0, 2, 0);
  h_y_charm_quark->Write(0, 2, 0);
  h_Nbeauty_quark->Write(0, 2, 0);
  h_pt_beauty_quark->Write(0, 2, 0);
  h_y_beauty_quark->Write(0, 2, 0);

  h_Ncharm_pairs->Write(0, 2, 0);
  h_Nbeauty_pairs->Write(0, 2, 0);

  rec_tree_mucharm->Write(0, 2, 0);
  rec_tree_mubeauty->Write(0, 2, 0);
  rec_tree_mumixed->Write(0, 2, 0);

  rec_tree_mucharm_withcut->Write(0, 2, 0);
  rec_tree_mubeauty_withcut->Write(0, 2, 0);
  rec_tree_mumixed_withcut->Write(0, 2, 0);

  rec_tree_mucharm_lowmass->Write(0, 2, 0);
  rec_tree_mubeauty_lowmass->Write(0, 2, 0);
  rec_tree_mumixed_lowmass->Write(0, 2, 0);

  rec_tree_mucharm_highmass->Write(0, 2, 0);
  rec_tree_mubeauty_highmass->Write(0, 2, 0);
  rec_tree_mumixed_highmass->Write(0, 2, 0);

  rec_tree_mucharm_lowmass_lowpt->Write(0, 2, 0);
  rec_tree_mubeauty_lowmass_lowpt->Write(0, 2, 0);
  rec_tree_mumixed_lowmass_lowpt->Write(0, 2, 0);

  fr.Close();

  Bool_t verbose = kFALSE;
  // TFile *f = TFile::Open(Form("%s/pyxsec_hists.root",dir_filename.Data()));
  // TList *list = (TList*)f->Get("cFilterList");
  // TH1F *h1Trials = (TH1F*)list->FindObject("h1Trials");
  // TH1F *h1Xsec = (TH1F*)list->FindObject("h1Xsec");

  TFile *outFile = new TFile(fileout.Data(), "UPDATE");
  outFile->cd();
  // h1Trials->Write(0,2,0);
  // h1Xsec->Write(0,2,0);
  outFile->ls();
  // h_PtYDiMu[0][0][0][0]->Write(0, 2, 0);
  // h_PtYDiMu[0][1][0][4]->Write(0, 2, 0);
  // h_PtYDiMu[0][1][1][4]->Write(0, 2, 0);
  // h_PtYDiMu[0][1][2][4]->Write(0, 2, 0);

  TString name_dir;

  // TDirectory *dir = outFile->GetDirectory("Muon");
  // if (!dir)
  //   dir = outFile->mkdir("Muon");
  // else
  //   printf("%s already exists \n", "Muon");
  // outFile->cd();

  TDirectory *dir = outFile->GetDirectory("DiMuon");
  if (!dir)
    dir = outFile->mkdir("DiMuon");
  else
    printf("%s already exists \n", "DiMuon");
  outFile->cd();

  for (Int_t a = 0; a < n_Mass_Cut; a++)
  {
    /*
    outFile->cd("Muon");

    for (Int_t Mu_b = 0; Mu_b < n_MuCut; Mu_b++)
    {
      dir = outFile->GetDirectory(Form("Muon/%s/%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data()));
      if (!dir)
        dir = outFile->mkdir(Form("Muon/%s/%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data()));
      else
        printf("%s already exists \n", Form("Muon/%s/%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data()));
      outFile->cd(Form("Muon/%s/%s", name_Mass_Cut[a].Data(), name_MuCut[Mu_b].Data()));
      for (Int_t Mu_c = 0; Mu_c < n_MuSelection; Mu_c++)
      {

        h_PtYMu[a][Mu_b][Mu_c]->Write(0, 2, 0);

        h_pdgMu[a][Mu_b][Mu_c]->Write(0, 2, 0);

        h_nMu_xevent[a][Mu_b][Mu_c]->Write(0, 2, 0);

      } // End saving over Mu charge

    } // End saving over Mu selection
    */
    outFile->cd();

    outFile->cd("DiMuon");

    for (Int_t DiMu_b = n_DiMuCut - 1; DiMu_b < n_DiMuCut; DiMu_b++)
    {
      for (Int_t DiMu_c = 0; DiMu_c < n_DiMu_Charge; DiMu_c++)
      {
        name_dir.Form("DiMuon/%s/%s/%s", name_Mass_Cut[a].Data(), name_DiMuCut[DiMu_b].Data(), name_DiMu_Charge[DiMu_c].Data());
        dir = outFile->GetDirectory(name_dir.Data());

        if (!dir)
          dir = outFile->mkdir(name_dir.Data());
        else
          printf("%s already exists \n", name_dir.Data());
        outFile->cd(name_dir.Data());
        for (Int_t DiMu_d = 0; DiMu_d < n_DiMuSelection; DiMu_d++)
        {

          h_PtYDiMu[a][DiMu_b][DiMu_c][DiMu_d]->Write(0, 2, 0);

          h_PtMDiMu[a][DiMu_b][DiMu_c][DiMu_d]->Write(0, 2, 0);

          h_MDiMu[a][DiMu_b][DiMu_c][DiMu_d]->Write(0, 2, 0);

          h_pdgDimuMu[a][DiMu_b][DiMu_c][DiMu_d]->Write(0, 2, 0);

          h_nDiMu_xevent[a][DiMu_b][DiMu_c][DiMu_d]->Write(0, 2, 0);

        } // End definition over DiMu selection

      } // End definition over DiMu cut
    }
  }
  printf("Finish Filling\n");
  h_nMu_xevent_rec_charm->Write(0, 2, 0);

  h_nMu_xevent_rec_beauty->Write(0, 2, 0);

  outFile->Close();
}