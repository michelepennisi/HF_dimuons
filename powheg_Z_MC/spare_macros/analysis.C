void analysis(string energy = "13TeV"){
    // CMS binning (https://arxiv.org/pdf/1812.10529.pdf) 13 TeV
    //const int n_bins = 8;
    //float mass_min[] = {20, 25, 30, 35, 40, 45, 50, 55};
    //float mass_max[] = {25, 30, 35, 40, 45, 50, 55, 60};
    //float mass_bins[] = {20, 25, 30, 35, 40, 45, 50, 55, 60};
    // ATLAS binning (https://arxiv.org/pdf/1404.1212.pdf) 7 TeV
    const int n_bins = 8;
    float mass_min[] = {26, 31, 36, 41, 46, 51, 56, 61};
    float mass_max[] = {31, 36, 41, 46, 51, 56, 61, 66};
    float mass_bins[] = {26, 31, 36, 41, 46, 51, 56, 61, 66};

    vector<int> *fPdgID = 0;
    vector<float> *fMass = 0;
    vector<float> *fPt = 0;
    vector<float> *fY = 0;

    vector<float> *fPx = 0;
    vector<float> *fPy = 0;
    vector<float> *fPz = 0;
    vector<float> *fE = 0;

    vector<int> *fMother1 = 0;
    vector<int> *fMother2 = 0;

    vector<int> index_mups;
    vector<int> index_mums;

    float px, py, pz, e = 0;
    float px_mup, py_mup, pz_mup, e_mup = 0;
    float px_mum, py_mum, pz_mum, e_mum = 0;


    TH1I *hist_mother_pdg = new TH1I("hist_mother_pdg", "; PDG", 1000, -0.5, 999.5);

    TH1F *hist_drell_yan_vs_mass[n_bins];
    TH1F *hist_drell_yan_vs_mass_ALICE[n_bins];
    TH1F *hist_drell_yan_vs_mass_CMS[n_bins];
    TH1F *hist_drell_yan_template_mass = new TH1F("hist_drell_yan_template_mass", "; #it{m} (GeV/c)", 30, 0, 30);
    TH1F *hist_drell_yan_template_pt = new TH1F("hist_drell_yan_template_pt", "; #it{p}_{T} (GeV/c^{2})", 30, 0, 30);
    for (int i = 0;i < n_bins;i++) {
        //hist_drell_yan_vs_mass[i] = new TH1F(Form("hist_drell_yan_vs_mass_%1.0f_%1.0f", mass_min[i], mass_max[i]), "; #it{m}_{T} (GeV/c^{2})", 8, 20, 60);
        //hist_drell_yan_vs_mass_ALICE[i] = new TH1F(Form("hist_drell_yan_vs_mass_ALICE_%1.0f_%1.0f", mass_min[i], mass_max[i]), "; #it{m}_{T} (GeV/c^{2})", 8, 20, 60);
        //hist_drell_yan_vs_mass_CMS[i] = new TH1F(Form("hist_drell_yan_vs_mass_CMS_%1.0f_%1.0f", mass_min[i], mass_max[i]), "; #it{m}_{T} (GeV/c^{2})", 8, 20, 60);
        hist_drell_yan_vs_mass[i] = new TH1F(Form("hist_drell_yan_vs_mass_%1.0f_%1.0f", mass_min[i], mass_max[i]), "; #it{m}_{T} (GeV/c^{2})", n_bins, mass_bins);
        hist_drell_yan_vs_mass_ALICE[i] = new TH1F(Form("hist_drell_yan_vs_mass_ALICE_%1.0f_%1.0f", mass_min[i], mass_max[i]), "; #it{m}_{T} (GeV/c^{2})", n_bins, mass_bins);
        hist_drell_yan_vs_mass_CMS[i] = new TH1F(Form("hist_drell_yan_vs_mass_CMS_%1.0f_%1.0f", mass_min[i], mass_max[i]), "; #it{m}_{T} (GeV/c^{2})", n_bins, mass_bins);
    }

    TH1F *hist_dimu_mass = new TH1F("hist_dimu_mass", "; #it{m}^{#mu#mu}", 1000, 0, 2000);
    TH1F *hist_dimu_pt = new TH1F("hist_dimu_pt", "; #it{p}_{T}^{#mu#mu}", 1000, 0, 100);
    TH1F *hist_dimu_y = new TH1F("hist_dimu_y", "; #it{y}^{#mu#mu}", 1000, -10, 10);
    TH1F *hist_dimu_y_CMS_cut = new TH1F("hist_dimu_y_CMS_cut", "; #it{y}^{#mu#mu}", 250, 0, 5);
    TH1F *hist_dimu_y_ALICE_cut = new TH1F("hist_dimu_y_ALICE_cut", "; #it{y}^{#mu#mu}", 250, 0, 5);
    TH2F *hist_dimu_y_mass = new TH2F("hist_dimu_y_mass", "; #it{y}^{#mu#mu}; #it{m}^{#mu#mu}", 1000, -10, 10, 1000, 0, 100);

    TH1F *hist_muon_pt = new TH1F("hist_muon_pt", "; #it{p}_{T}^{#mu}", 1000, 0, 100);
    TH1F *hist_muon_y = new TH1F("hist_muon_y", "; #it{y}^{#mu}", 1000, -10, 10);

    TH2F *hist_muon_pt_dimu_mass = new TH2F("hist_muon_pt_dimu_mass", "; #it{p}_{T}^{#mu}; #it{m}^{#mu#mu}", 1000, 0, 100, 1000, 0, 100);
    
    TFile *fIn = new TFile(Form("data/root_tree_%s.root", energy.c_str()), "READ");

    TTree *tree = (TTree*) fIn -> Get("mytree");
    tree -> SetBranchAddress("pdgID", &fPdgID);
    tree -> SetBranchAddress("Mass", &fMass);
    tree -> SetBranchAddress("Pt", &fPt);
    tree -> SetBranchAddress("Y", &fY);

    tree -> SetBranchAddress("px", &fPx);
    tree -> SetBranchAddress("py", &fPy);
    tree -> SetBranchAddress("pz", &fPz);
    tree -> SetBranchAddress("E", &fE);

    tree -> SetBranchAddress("mother1", &fMother1);
    tree -> SetBranchAddress("mother2", &fMother2);

    for (int iEntry = 0;iEntry < tree -> GetEntries();iEntry++) {
        tree -> GetEntry(iEntry);
        for (ULong_t iTrack = 0; iTrack < fPdgID -> size(); iTrack++) {
            if (fPdgID -> at(iTrack) == 22 || fPdgID -> at(iTrack) == 23) {
                hist_drell_yan_template_mass -> Fill(fMass -> at(iTrack));
                hist_drell_yan_template_pt -> Fill(fPt -> at(iTrack));
                for (int iMass = 0;iMass < 8;iMass++) {
                    px = fPx -> at(iTrack);
                    py = fPy -> at(iTrack);
                    pz = fPz -> at(iTrack);
                    e = fE -> at(iTrack);
                    TLorentzVector mother(px, py, pz, e);
                    if (fMass -> at(iTrack) > mass_min[iMass] && fMass -> at(iTrack) < mass_max[iMass]) {
                        hist_drell_yan_vs_mass[iMass] -> Fill(fMass -> at(iTrack));

                        //if (TMath::Abs(mother.Eta()) > 2.5 && TMath::Abs(mother.Eta()) < 4) {
                            //hist_drell_yan_vs_mass_ALICE[iMass] -> Fill(fMass -> at(iTrack));
                        //}
                        if (mother.Eta() > 2.5 && mother.Eta() < 4) {
                            hist_drell_yan_vs_mass_ALICE[iMass] -> Fill(fMass -> at(iTrack));
                        }
                        if (TMath::Abs(mother.Eta()) < 2.4) {
                            hist_drell_yan_vs_mass_CMS[iMass] -> Fill(fMass -> at(iTrack));
                        }
                    }
                }
            }
            if (fPdgID -> at(iTrack) == 13) {
                index_mups.push_back(iTrack);
                hist_muon_pt -> Fill(fPt -> at(iTrack));
                hist_muon_y -> Fill(fY -> at(iTrack));
                int mother_index1 = fMother1 -> at(iTrack);
                int mother_index2 = fMother2 -> at(iTrack);
                std::cout << "M1 = " << fPdgID -> at(mother_index1-3) << " M2 = " << fPdgID -> at(mother_index2-3) << std::endl;
                hist_mother_pdg -> Fill(fPdgID -> at(mother_index1-3));
            } 
            if (fPdgID -> at(iTrack) == -13) {
                index_mums.push_back(iTrack);
                hist_muon_pt -> Fill(fPt -> at(iTrack));
                hist_muon_y -> Fill(fY -> at(iTrack));
                int mother_index1 = fMother1 -> at(iTrack);
                int mother_index2 = fMother2 -> at(iTrack);
                std::cout << "M1 = " << fPdgID -> at(mother_index1-3) << " M2 = " << fPdgID -> at(mother_index2-3) << std::endl;
                hist_mother_pdg -> Fill(fPdgID -> at(mother_index1-3));
            }
            //std::cout << iEntry << " : " << iTrack << " : " << fMass -> at(iTrack) << std::endl;
        }
        std::cout << iEntry << " -> N part = " << fPdgID -> size() << "  : N mup = " << index_mups.size() << " : N mum = " << index_mups.size() << std::endl;

        for (ULong_t iTrack1 = 0; iTrack1 < index_mums.size(); iTrack1++) {
            int index_mum = index_mums.at(iTrack1);
            px_mum = fPx -> at(index_mum);
            py_mum = fPy -> at(index_mum);
            pz_mum = fPz -> at(index_mum);
            e_mum = fE -> at(index_mum);
            TLorentzVector mum(px_mum, py_mum, pz_mum, e_mum);

            for (ULong_t iTrack2 = 0; iTrack2 < index_mums.size(); iTrack2++) {
                int index_mup = index_mups.at(iTrack2);
                px_mup = fPx -> at(index_mup);
                py_mup = fPy -> at(index_mup);
                pz_mup = fPz -> at(index_mup);
                e_mup = fE -> at(index_mup);
                TLorentzVector mup(px_mup, py_mup, pz_mup, e_mup);

                TLorentzVector dimuon = mum + mup;
                hist_dimu_mass -> Fill(dimuon.M());
                hist_dimu_pt -> Fill(dimuon.Pt());
                hist_dimu_y -> Fill(dimuon.Rapidity());
                hist_dimu_y_mass -> Fill(dimuon.Rapidity(), dimuon.M());

                hist_muon_pt_dimu_mass -> Fill(mup.Pt(), dimuon.M());

                // CMS cuts (https://arxiv.org/pdf/1310.7291.pdf)
                if (dimuon.M() > 20 && dimuon.M() < 30) {
                    if (mum.Pt() > 10 && mup.Pt() > 10) {
                        if (TMath::Abs(mum.Eta()) < 2.4 && TMath::Abs(mup.Eta()) < 2.4) {
                            hist_dimu_y_CMS_cut -> Fill(TMath::Abs(dimuon.Rapidity()));
                        }
                    }
                }

                // ALICE cuts
                if (dimuon.M() > 20 && dimuon.M() < 30) {
                    if (mum.Pt() > 0.5 && mup.Pt() > 0.5) {
                        if (TMath::Abs(mum.Eta()) > 2.5 && TMath::Abs(mum.Eta()) < 4) {
                            if (TMath::Abs(mup.Eta()) > 2.5 && TMath::Abs(mup.Eta()) < 4) {
                                hist_dimu_y_ALICE_cut -> Fill(TMath::Abs(dimuon.Rapidity()));
                            }
                        }
                    }
                }
            }
        }
        index_mups.clear();
        index_mums.clear();
    }

    TFile *fOut = new TFile(Form("AnalysisResults_%s.root", energy.c_str()), "RECREATE");
    hist_mother_pdg -> Write();
    hist_dimu_mass -> Write();
    hist_dimu_pt -> Write();
    hist_dimu_y -> Write();
    hist_dimu_y_CMS_cut -> Write();
    hist_dimu_y_ALICE_cut -> Write();
    hist_dimu_y_mass -> Write();
    hist_muon_pt -> Write();
    hist_muon_y -> Write();
    hist_muon_pt_dimu_mass -> Write();

    hist_drell_yan_template_mass -> Write();
    hist_drell_yan_template_pt -> Write();

    for (int iMass = 0;iMass < n_bins;iMass++) {
        hist_drell_yan_vs_mass[iMass] -> Write();
        hist_drell_yan_vs_mass_ALICE[iMass] -> Write();
        hist_drell_yan_vs_mass_CMS[iMass] -> Write();
    }
    fOut -> Close();
}

void compare_cross_section(){
    const float drell_yan_cross_section_13TeV = 4122.4; // pb
    const float drell_yan_cross_section_7TeV =  2246.8; // pb

    float n_events = 1.e6; // number of events

    // https://arxiv.org/pdf/1812.10529.pdf
    const int n_bins_CMS = 8;
    float mass_min_CMS[] = {20, 25, 30, 35, 40, 45, 50, 55};
    float mass_max_CMS[] = {25, 30, 35, 40, 45, 50, 55, 60};
    float mass_bins_CMS[] = {20, 25, 30, 35, 40, 45, 50, 55, 60};
    float mass_width_CMS[] = {5, 5, 5, 5, 5, 5, 5, 5};
    float cross_section_CMS[] = {99, 53, 28, 17, 12, 8.5, 6.3, 5.3};

    TFile *fIn_13TeV = new TFile("AnalysisResults_13TeV.root", "READ");
    TH1F *hist_drell_yan_POWHEG_13TeV = new TH1F("hist_drell_yan_POWHEG_13TeV", "; #it{m}^{#mu#mu}", n_bins_CMS, mass_bins_CMS);
    hist_drell_yan_POWHEG_13TeV -> SetLineColor(kAzure+2);

    TH1F *hist_drell_yan_POWHEG_ALICE_cuts = new TH1F("hist_drell_yan_POWHEG_ALICE_cuts", "; #it{m}^{#mu#mu}", n_bins_CMS, mass_bins_CMS);
    hist_drell_yan_POWHEG_ALICE_cuts -> SetLineColor(kViolet+1);

    TH1F *hist_drell_yan_POWHEG_CMS_cuts = new TH1F("hist_drell_yan_POWHEG_CMS_cuts", "; #it{m}^{#mu#mu}", n_bins_CMS, mass_bins_CMS);
    hist_drell_yan_POWHEG_CMS_cuts -> SetLineColor(kOrange+1);

    for (int i = 0;i < n_bins_CMS;i++) {
        TH1F *hist = (TH1F*) fIn_13TeV -> Get(Form("hist_drell_yan_vs_mass_%1.0f_%1.0f", mass_min_CMS[i], mass_max_CMS[i]));
        //float entries = (hist -> GetEntries() / n_events) *  drell_yan_cross_section_13TeV;
        float entries = (hist -> GetEntries() / n_events) * (drell_yan_cross_section_13TeV / mass_width_CMS[i]);
        float err_entries = (1. / TMath::Sqrt(hist -> GetEntries())) * entries;
        //std::cout << mass_min_CMS[i] << " " << mass_max_CMS[i] << " " << entries << " " << err_entries << std::endl;
        hist_drell_yan_POWHEG_13TeV -> SetBinContent(i+1, entries);
        hist_drell_yan_POWHEG_13TeV -> SetBinError(i+1, err_entries);

        TH1F *hist_ALICE = (TH1F*) fIn_13TeV -> Get(Form("hist_drell_yan_vs_mass_ALICE_%1.0f_%1.0f", mass_min_CMS[i], mass_max_CMS[i]));
        entries = (hist_ALICE -> GetEntries() / n_events) * (drell_yan_cross_section_13TeV / mass_width_CMS[i]);
        err_entries = (1. / TMath::Sqrt(hist_ALICE -> GetEntries())) * entries;
        std::cout << mass_min_CMS[i] << " " << mass_max_CMS[i] << " " << entries << " " << err_entries << std::endl;
        hist_drell_yan_POWHEG_ALICE_cuts -> SetBinContent(i+1, entries);
        hist_drell_yan_POWHEG_ALICE_cuts -> SetBinError(i+1, err_entries);

        TH1F *hist_CMS = (TH1F*) fIn_13TeV -> Get(Form("hist_drell_yan_vs_mass_CMS_%1.0f_%1.0f", mass_min_CMS[i], mass_max_CMS[i]));
        entries = (hist_CMS -> GetEntries() / n_events) * (drell_yan_cross_section_13TeV / mass_width_CMS[i]);
        err_entries = (1. / TMath::Sqrt(hist_CMS -> GetEntries())) * entries;
        //std::cout << entries << " +/- " << err_entries << std::endl;
        hist_drell_yan_POWHEG_CMS_cuts -> SetBinContent(i+1, entries);
        hist_drell_yan_POWHEG_CMS_cuts -> SetBinError(i+1, err_entries);
    }

    return;

    TH1F *hist_drell_yan_CMS = new TH1F("hist_drell_yan_CMS", "", 8, 20, 60);
    hist_drell_yan_CMS -> SetLineColor(kRed+1);

    std::cout << "---------" << std::endl;
    for (int i = 0;i < n_bins_CMS;i++) {
        float entries = cross_section_CMS[i];
        float err_entries = 0;
        //hist_drell_yan_CMS -> SetBinContent(i+1, (cross_section_CMS[i] / 1000.) * mass_width_CMS[i]);
        hist_drell_yan_CMS -> SetBinContent(i+1, entries);
        hist_drell_yan_CMS -> SetBinError(i+1, err_entries);
        std::cout << mass_min_CMS[i] << " " << mass_max_CMS[i] << " " << entries << " " << err_entries << std::endl;
    }

    TCanvas *canvas_13TeV = new TCanvas("comparison_13TeV", "", 800, 600);
    gPad -> SetLogy(1);
    hist_drell_yan_POWHEG_13TeV -> Draw("EP");
    hist_drell_yan_CMS -> Draw("H SAME");
    hist_drell_yan_POWHEG_ALICE_cuts -> Draw("EP SAME");
    hist_drell_yan_POWHEG_CMS_cuts -> Draw("EP SAME");
    canvas_13TeV -> BuildLegend();
    

    TFile *fOut_13TeV = new TFile("comparison_13TeV.root", "RECREATE");
    hist_drell_yan_CMS -> Write();
    hist_drell_yan_POWHEG_13TeV -> Write();
    hist_drell_yan_POWHEG_ALICE_cuts -> Write();
    hist_drell_yan_POWHEG_CMS_cuts -> Write();
    canvas_13TeV -> Write();
    fOut_13TeV -> Close();

    //https://arxiv.org/pdf/1404.1212.pdf
    const int n_bins_ATLAS = 8;
    float mass_min_ATLAS[] = {26, 31, 36, 41, 46, 51, 56, 61};
    float mass_max_ATLAS[] = {31, 36, 41, 46, 51, 56, 61, 66};
    float mass_bins_ATLAS[] = {26, 31, 36, 41, 46, 51, 56, 61, 66};
    float mass_width_ATLAS[] = {5, 5, 5, 5, 5, 5, 5, 5, 5};
    //float cross_section_ATLAS[] = {1.80, 3.12, 2.64, 2.03, 1.54, 1.19, 1.00, 0.90}; // Table 8
    float cross_section_ATLAS[] = {27.9, 16.5, 9.5, 6.1, 4.2, 3.2, 2.4, 2.1}; // Fig. 5 

    TFile *fIn_7TeV = new TFile("AnalysisResults_7TeV.root", "READ");

    TH1F *hist_drell_yan_POWHEG_7TeV = new TH1F("hist_drell_yan_POWHEG_7TeV", "; #it{m}^{#mu#mu}", n_bins_ATLAS, mass_bins_ATLAS);
    hist_drell_yan_POWHEG_7TeV -> SetLineColor(kAzure+2);

    for (int i = 0;i < n_bins_ATLAS;i++) {
        TH1F *hist = (TH1F*) fIn_7TeV -> Get(Form("hist_drell_yan_vs_mass_%1.0f_%1.0f", mass_min_ATLAS[i], mass_max_ATLAS[i]));
        //float entries = (hist -> GetEntries() / n_events) *  drell_yan_cross_section_7TeV;
        float entries = (hist -> GetEntries() / n_events) *  (drell_yan_cross_section_7TeV / mass_width_ATLAS[i]);
        float err_entries = (1. / TMath::Sqrt(hist -> GetEntries())) * entries;
        std::cout << mass_min_ATLAS[i] << " " << mass_max_ATLAS[i] << " " << entries << " " << err_entries << std::endl;
        hist_drell_yan_POWHEG_7TeV -> SetBinContent(i+1, entries);
        hist_drell_yan_POWHEG_7TeV -> SetBinError(i+1, err_entries);
    }

    TH1F *hist_drell_yan_ATLAS = new TH1F("hist_drell_yan_ATLAS", "", n_bins_ATLAS, mass_bins_ATLAS);
    hist_drell_yan_ATLAS -> SetLineColor(kRed+1);

    std::cout << "---------" << std::endl;
    for (int i = 0;i < n_bins_ATLAS;i++) {
        float entries = cross_section_ATLAS[i];
        float err_entries = 0;
        //hist_drell_yan_ATLAS -> SetBinContent(i+1, (cross_section_ATLAS[i] / 1000.) * mass_width_ATLAS[i]);
        hist_drell_yan_ATLAS -> SetBinContent(i+1, entries);
        hist_drell_yan_ATLAS -> SetBinError(i+1, err_entries);
        std::cout << mass_min_ATLAS[i] << " " << mass_max_ATLAS[i] << " " << entries << " " << err_entries << std::endl;
    }

    TCanvas *canvas_7TeV = new TCanvas("comparison_7TeV", "", 800, 600);
    gPad -> SetLogy(1);
    hist_drell_yan_POWHEG_7TeV -> Draw("EP");
    hist_drell_yan_ATLAS -> Draw("H SAME");
    canvas_7TeV -> BuildLegend();

    TFile *fOut_7TeV = new TFile("comparison_7TeV.root", "RECREATE");
    hist_drell_yan_ATLAS -> Write();
    hist_drell_yan_POWHEG_7TeV -> Write();
    canvas_7TeV -> Write();
    fOut_7TeV -> Close();
}