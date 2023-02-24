#include <iostream>
#include "Riostream.h"
#include <fstream>
#include <vector>
#include "Pythia8/Pythia.h"
#include "Pythia8/Analysis.h"
#include "Pythia8Plugins/ColourReconnectionHooks.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "THnSparse.h"
#include "TMath.h"

#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TLorentzVector.h"
// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TLatex.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
// Versione aggiornata al 3/11/21
// scambio-> oggetti mu- contengono i mu+ e viceversa
// Salvataggio in due root file, il primo completo mentre il secondo contiene solo l'info relativi ai muoni
using namespace Pythia8;
Pythia pythia;

bool isCharm(int i, int iEvent, bool verbose);
bool IsChargedPhysicalPrimary(Particle);
bool IsPIDPhysicalPrimary(Particle);
void ConfigHeavyFlavor();
void AtlasTuning();

int main(int argc, char *argv[])
{
    char *selectedprocess[5] = {"SoftQCD", "ccbar", "bbbar", "SoftQCDall", "Drell-Yan"};
    char *change[4] = {"Def", "Config", "Config+Atlas", "Mode2"};
    char *selectedBR[2] = {"DefaultBR", "ForceBR"};
    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);
    // printf("argc = %d\n argv = %s\n",theApp.Argc(),theApp.Argv(0));
    int nevents = 1000;
    int seed = 9170;
    int mode = -1;  //-1==Monash  .  1 not allowed. Check JHEP08(2015)003
    int n_MPI = -1; //-1==we take everything --> >0--> register only larger nMPI events
    int chooseprocess = 1;
    int choosechange = 0; // 1 = SoftQCD, 2 = ccbar, 3 = bbbar
    int BR = 0;           // 0 = original BR, 1 = force semileptonic decay of charm

    // read arguments to change settings
    for (int ig = 0; ig < argc; ++ig)
    {
        TString str = argv[ig];
        printf("str %d %s\n", ig, argv[ig]);
        if (str.Contains("--events"))
        {
            sscanf(argv[ig + 1], "%d", &nevents);
            printf("nevents after setting: %d\n", nevents);
        }
        if (str.Contains("--seed"))
        {
            sscanf(argv[ig + 1], "%d", &seed);
            printf("seed after setting: %d\n", seed);
        }
        if (str.Contains("--mode"))
        {
            sscanf(argv[ig + 1], "%d", &mode);
            printf("mode after setting: %d\n", mode);
            if (mode != -1 && mode != 0 && mode != 2 && mode != 3)
            {
                printf("mode must be -1, 0, 2 or 3. Try again\n");
                return 0;
            }
        }
        if (str.Contains("--triggerMPI"))
        {
            sscanf(argv[ig + 1], "%d", &n_MPI);
            printf("we trigger on %d MPI %d\n", n_MPI);
        }
        if (str.Contains("--process"))
        {
            sscanf(argv[ig + 1], "%d", &chooseprocess);
            //            printf("generation type %d\n",chooseprocess);
            printf("process after setting %s\n", selectedprocess[chooseprocess - 1]);
        }
        if (str.Contains("--choosechange"))
        {
            sscanf(argv[ig + 1], "%d", &choosechange);
            //            printf("generation type %d\n",choosechange);
            printf("mode proceess after setting %s\n", change[choosechange]);
        }
        if (str.Contains("--BR"))
        {
            sscanf(argv[ig + 1], "%d", &BR);
            //            printf("force semileptonic D %d\n",BR);
            printf("BR after setting %s\n", selectedBR[BR]);
        }
    }

    Settings &settings = pythia.settings;
    cout << "PROBLEMA" << endl;
    // Basic settings (collision energ, processes....)

    if (chooseprocess == 1)
    {
        if (choosechange == 0)
        {
            cout << "choosechange 0" << endl;
            pythia.readString("SoftQCD:nonDiffractive = on ");
        }
        else if (choosechange == 1)
        {
            cout << "ConfigHeavyFlavor" << endl;
            pythia.readString("SoftQCD:nonDiffractive = on ");

            // ConfigHeavyFlavor
            ConfigHeavyFlavor();

            // Mix Set
            pythia.readString("SigmaProcess:factorMultFac = 1.");
            // Intrinsic <kT>
            pythia.readString("BeamRemnants:primordialKT = on");
            pythia.readString("BeamRemnants:primordialKTsoft = 0.");
            pythia.readString("BeamRemnants:primordialKThard = 1.0");
            pythia.readString("BeamRemnants:halfScaleForKT = 0.");
            pythia.readString("BeamRemnants:halfMassForKT = 0.");
            // Set c and b quark masses
            pythia.readString("ParticleData:mcRun = 1.20");
            pythia.readString("ParticleData:mbRun = 4.75");
        }
        else if (choosechange == 2)
        {
            cout << "ConfigHeavyFlavor+AtlasTuning" << endl;
            pythia.readString("SoftQCD:nonDiffractive = on ");
            ConfigHeavyFlavor();
            // MIX SET
            pythia.readString("SigmaProcess:factorMultFac = 1.");
            // Intrinsic <kT>
            pythia.readString("BeamRemnants:primordialKT = on");
            pythia.readString("BeamRemnants:primordialKTsoft = 0.");
            pythia.readString("BeamRemnants:primordialKThard = 1.0");
            pythia.readString("BeamRemnants:halfScaleForKT = 0.");
            pythia.readString("BeamRemnants:halfMassForKT = 0.");
            // Set c and b quark masses
            pythia.readString("ParticleData:mcRun = 1.20");
            pythia.readString("ParticleData:mbRun = 4.75");

            AtlasTuning();
        }
        else if (choosechange == 3)
        {
            cout << "Mode2" << endl;
            pythia.readString("SoftQCD:nonDiffractive = on");
            pythia.readString(Form("Tune:pp = %d", 14)); // 14: Monash 2013
            pythia.readString(Form("Beams:eCM = %f", (float)13000));
            // initialize random seed
            pythia.readString("Random:setSeed = on");
            pythia.readString(Form("Random:seed %d", seed));

            pythia.readString("ColourReconnection:mode = 1");
            cout << "ColourReconnection:mode = 1" << endl;
            pythia.readString("ColourReconnection:allowDoubleJunRem = off");
            pythia.readString("ColourReconnection:m0 = 0.3");
            pythia.readString("ColourReconnection:allowJunctions = on");
            pythia.readString("ColourReconnection:junctionCorrection = 1.20");
            pythia.readString("ColourReconnection:timeDilationMode = 2");
            pythia.readString("ColourReconnection:timeDilationPar = 0.18");
            pythia.readString("StringPT:sigma = 0.335");
            pythia.readString("StringZ:aLund = 0.36");
            pythia.readString("StringZ:bLund = 0.56");
            pythia.readString("StringFlav:probQQtoQ = 0.078");
            pythia.readString("StringFlav:ProbStoUD = 0.2");
            pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
            pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
            pythia.readString("BeamRemnants:remnantMode = 1");
            cout << "BeamRemnants:remnantMode = 1" << endl;
            pythia.readString("BeamRemnants:saturation =5");
        }
        //        pythia.readString("SoftQCD:nonDiffractive = on ");
        ////       // ConfigHeavyFlavor();
        ////        // QCD scales
        //        pythia.readString("SigmaProcess:factorMultFac = 1.");
        //        // Intrinsic <kT>
        //        pythia.readString("BeamRemnants:primordialKT = on");
        //        pythia.readString("BeamRemnants:primordialKTsoft = 0.");
        //        pythia.readString("BeamRemnants:primordialKThard = 1.0");
        //        pythia.readString("BeamRemnants:halfScaleForKT = 0.");
        //        pythia.readString("BeamRemnants:halfMassForKT = 0.");
        //        // Set c and b quark masses
        //        pythia.readString("ParticleData:mcRun = 1.20");
        //        pythia.readString("ParticleData:mbRun = 4.75");
    }
    else if (chooseprocess == 2)
        pythia.readString(" HardQCD:hardccbar = on ");
    else if (chooseprocess == 3)
        pythia.readString(" HardQCD:hardbbbar = on "); // forza ad avere almeno una coppia di b e bbar prodotti
    else if (chooseprocess == 4)
        pythia.readString(" SoftQCD:all = on ");
    else if (chooseprocess == 5)
    {
        printf("DRELL-YAN");
        pythia.readString("Beams:idA = 2212");   //       ! first beam, p = 2212, pbar = -2212
        pythia.readString("Beams:idB = 2212");   //       ! second beam, p = 2212, pbar = -2212
        pythia.readString("Beams:eCM = 13000."); //       ! CM energy of collision
        pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
        pythia.readString("WeakZ0:gmZmode = 1 ");
        pythia.readString("PhaseSpace:mHatMin = 4.");
    }
    else
    {
        printf("Invalid process!\n");
        exit(0);
    }
    // forza ad avere decadimento buono in coppie di ccbar e bbar
    char outfilename1[100];
    char outfilename2[100];

    // const char *path="/Volumes/REPOSITORY/Test/";

    const char *path = "/home/michele_pennisi/cernbox/HF_dimuons/pythia_stand/sim";
    sprintf(outfilename1, "%s/pythia_sim_%s_%s_%d_%d_%s.root", path, selectedprocess[chooseprocess - 1], change[choosechange], nevents, seed, selectedBR[BR]);

    if (choosechange != 3)
    {
        pythia.readString(Form("Tune:pp = %d", 14)); // 14: Monash 2013
        pythia.readString(Form("Beams:eCM = %f", (float)13000));
        // initialize random seed
        pythia.readString("Random:setSeed = on");
        pythia.readString(Form("Random:seed %d", seed));
    }

    //   if(mode==-1){ // Monash mode
    //     pythia.readString("ColourReconnection:mode = 0");
    //     pythia.readString("ColourReconnection:allowDoubleJunRem = on");
    //     pythia.readString("StringPT:sigma = 0.335");
    //     pythia.readString("StringZ:aLund = 0.68");
    //     pythia.readString("StringZ:bLund = 0.98");
    //     pythia.readString("StringFlav:probQQtoQ = 0.081");
    //     pythia.readString("StringFlav:ProbStoUD = 0.217");
    //     pythia.readString("StringFlav:probQQ1toQQ0join = 0.5,0.7,0.9,1.0");
    //     pythia.readString("MultiPartonInteractions:pT0Ref = 2.28");
    //     pythia.readString("BeamRemnants:remnantMode = 0");
    //   }

    pythia.init();
    // Switch off decays not involving muons

    if (chooseprocess == 5 && BR == 1)
    {
        pythia.readString("23:onMode = 0");
        pythia.readString("23:onIfAny = 13");
        pythia.particleData.list(23);
    }
    if (chooseprocess != 5 && BR == 1)
    {
        pythia.readString("411:onMode = 0");
        pythia.readString("411:onIfAny = 13");
        //        pythia.particleData.list(411);
        pythia.readString("421:onMode = 0");
        pythia.readString("421:onIfAny = 13");
        //        pythia.particleData.list(421);
        pythia.readString("431:onMode = 0");
        pythia.readString("431:onIfAny = 13");
        //        pythia.particleData.list(431);
        pythia.readString("4122:onMode = 0");
        pythia.readString("4122:onIfAny = 13");
        //        pythia.particleData.list(4122);
        pythia.readString("4132:onMode = 0");
        pythia.readString("4132:onIfAny = 13");
        //        pythia.particleData.list(4132);
        pythia.readString("4232:onMode = 0");
        pythia.readString("4232:onIfAny = 13");
        //        pythia.particleData.list(4232);

        //        pythia.particleData.list(5);
        //        pythia.readString("511:onMode = 0");
        //        pythia.readString("511:onIfAny = 13");
        //        pythia.particleData.list(511);
        //        pythia.readString("521:onMode = 0");
        //        pythia.readString("521:onIfAny = 13");
        //        pythia.particleData.list(521);
        //        pythia.readString("531:onMode = 0");
        //        pythia.readString("531:onIfAny = 13");
        //        pythia.particleData.list(531);
        //        pythia.readString("5122:onMode = 0");
        //        pythia.readString("5122:onIfAny = 13");
        //        pythia.particleData.list(5122);
        //        pythia.readString("5132:onMode = 0");
        //        pythia.readString("5132:onIfAny = 13");
        //        pythia.particleData.list(5132);
        //        pythia.readString("5232:onMode = 0");
        //        pythia.readString("5232:onIfAny = 13");
        //        pythia.particleData.list(5232);
    }

    Int_t fNMuons_gen = 0;     // gen muon in the event
    Int_t fNDimu_gen;          // gen dimuons in the event
    Int_t fN_HFquarks_gen = 0; // gen c/cbar or b/bar HFquarks in the event
    static const Int_t fMuons_dim = 500;
    static const Int_t fDimu_dim = 500;
    Int_t fPDG_HFquark_gen[fMuons_dim];           // single gen c/cbar PDG mum
    Int_t fPDG_HFquark_gen_daughter1[fMuons_dim]; // single gen c/cbar PDG mum
    Int_t fPDG_HFquark_gen_daughter2[fMuons_dim]; // single gen c/cbar PDG mum
    Double_t fPt_HFquark_gen[fMuons_dim];         // single gen c/cbar or b/bbar HFquark pT
    Double_t fY_HFquark_gen[fMuons_dim];          // single gen c/cbar or b/bbar HFquark y

    Int_t fPDGmum_gen[fMuons_dim];    // single gen mu PDG mum
    Int_t fPDG_gen[fMuons_dim];       // single gen mu PDG mum
    Int_t fPromptmum_gen[fMuons_dim]; // single gen mu PDG mum
    Double_t fPt_mum_gen[fMuons_dim]; // pt of single mu mum
    Double_t fY_mum_gen[fMuons_dim];  // Y of single mu mum
    Double_t fPt_gen[fMuons_dim];     // single gen mu pT
    Double_t fE_gen[fMuons_dim];      // single gen mu E
    Double_t fPx_gen[fMuons_dim];     // single gen mu px
    Double_t fPy_gen[fMuons_dim];     // single gen mu py
    Double_t fPz_gen[fMuons_dim];     // single gen mu pz
    Double_t fY_gen[fMuons_dim];      // single gen mu y
    Double_t fEta_gen[fMuons_dim];    // single gen mu eta
    Int_t fCharge_gen[fMuons_dim];    // single gen mu charge
    Double_t fPhi_gen[fMuons_dim];    // single gen mu phi
    Double_t fTheta_gen[fMuons_dim];  // single gen mu theta

    Int_t fNHadron_gen = 0;
    Int_t fPDGHadron_gen[fMuons_dim];
    Int_t fPromptHadron_gen[fMuons_dim];
    Double_t fHadron_Pt_gen[fMuons_dim];
    Double_t fHadron_E_gen[fMuons_dim];
    Double_t fHadron_Px_gen[fMuons_dim];
    Double_t fHadron_Py_gen[fMuons_dim];
    Double_t fHadron_Pz_gen[fMuons_dim];
    Double_t fHadron_Y_gen[fMuons_dim];
    Double_t fHadron_Eta_gen[fMuons_dim];
    Double_t fHadron_Phi_gen[fMuons_dim];
    Double_t fHadron_Theta_gen[fMuons_dim];

    Int_t fDimuMu_gen[fDimu_dim][2];   // reference to single gen mus
    Double_t fDimuPt_gen[fDimu_dim];   // gen dimuon pT
    Double_t fDimuPx_gen[fDimu_dim];   // gen dimuon px
    Double_t fDimuPy_gen[fDimu_dim];   // gen dimuon py
    Double_t fDimuPz_gen[fDimu_dim];   // gen dimuon pz
    Double_t fDimuY_gen[fDimu_dim];    // gen dimuon y
    Double_t fDimuMass_gen[fDimu_dim]; // gen dimuon invariant mass
    Int_t fDimuCharge_gen[fDimu_dim];  // gen dimuon charge

    for (Int_t i = 0; i < fMuons_dim; i++)
    {
        fPDG_HFquark_gen[i] = 9999.;
        fPt_HFquark_gen[i] = 9999.;
        fY_HFquark_gen[i] = 9999.;

        fPDGmum_gen[i] = 9999.;
        fPDG_gen[i] = 9999.;
        fPromptmum_gen[i] = 9999.;
        fPt_mum_gen[i] = 9999.;
        fY_mum_gen[i] = 9999.;
        fPt_gen[i] = 999.;
        fE_gen[i] = 999.;
        fPx_gen[i] = 999;
        fPy_gen[i] = 999;
        fPz_gen[i] = 999;
        fY_gen[i] = 999.;
        fEta_gen[i] = 999.;
        fPhi_gen[i] = 999.;
        fTheta_gen[i] = 999.;
        fCharge_gen[i] = 999.;

        fPDGHadron_gen[i] = 999.;
        fPromptHadron_gen[i] = 999.;

        fHadron_Pt_gen[i] = 999.;
        fHadron_E_gen[i] = 999.;
        fHadron_Px_gen[i] = 999.;
        fHadron_Py_gen[i] = 999.;
        fHadron_Pz_gen[i] = 999.;
        fHadron_Y_gen[i] = 999.;
        fHadron_Eta_gen[i] = 999.;
        fHadron_Phi_gen[i] = 999.;
        fHadron_Theta_gen[i] = 999.;
    }

    TFile *outFile1 = new TFile(outfilename1, "RECREATE");
    TTree *fOutputTree = new TTree("MCTree", "Data Tree");

    fOutputTree = new TTree("MCTree", "Data Tree");

    // fOutputTree->Branch("PercentV0M",&fPercentV0M,"PercentV0M/D");
    fOutputTree->Branch("N_HFquarks_gen", &fN_HFquarks_gen, "N_HFquarks_gen/I");
    fOutputTree->Branch("PDG_HFquark_gen", fPDG_HFquark_gen, "PDG_HFquark_gen[N_HFquarks_gen]/I");
    fOutputTree->Branch("PDG_HFquark_gen_daughter1", fPDG_HFquark_gen_daughter1, "PDG_HFquark_gen_daughter1[N_HFquarks_gen]/I");
    fOutputTree->Branch("PDG_HFquark_gen_daughter2", fPDG_HFquark_gen_daughter2, "PDG_HFquark_gen_daughter2[N_HFquarks_gen]/I");
    fOutputTree->Branch("Pt_HFquark_gen", fPt_HFquark_gen, "Pt_HFquark_gen[N_HFquarks_gen]/D");
    fOutputTree->Branch("Y_HFquark_gen", fY_HFquark_gen, "Y_HFquark_gen[N_HFquarks_gen]/D");

    // fOutputTree->Branch("NDimu_gen", &fNDimu_gen, "NDimu_gen/I");
    // fOutputTree->Branch("DimuMu_gen", fDimuMu_gen, "DimuMu_gen[NDimu_gen][2]/I");
    // fOutputTree->Branch("DimuPt_gen", fDimuPt_gen, "DimuPt_gen[NDimu_gen]/D");
    // fOutputTree->Branch("DimuPx_gen", fDimuPx_gen, "DimuPx_gen[NDimu_gen]/D");
    // fOutputTree->Branch("DimuPy_gen", fDimuPy_gen, "DimuPy_gen[NDimu_gen]/D");
    // fOutputTree->Branch("DimuPz_gen", fDimuPz_gen, "DimuPz_gen[NDimu_gen]/D");
    // fOutputTree->Branch("DimuY_gen", fDimuY_gen, "DimuY_gen[NDimu_gen]/D");
    // fOutputTree->Branch("DimuMass_gen", fDimuMass_gen, "DimuMass_gen[NDimu_gen]/D");
    // fOutputTree->Branch("DimuCharge_gen", fDimuCharge_gen, "DimuCharge_gen[NDimu_gen]/I");

    fOutputTree->Branch("NMuons_gen", &fNMuons_gen, "NMuons_gen/I");
    fOutputTree->Branch("PDGmum_gen", fPDGmum_gen, "PDGmum_gen[NMuons_gen]/I");
    fOutputTree->Branch("Promptmum_gen", fPromptmum_gen, "Promptmum_gen[NMuons_gen]/I");
    fOutputTree->Branch("PDG_gen", fPDG_gen, "PDG_gen[NMuons_gen]/I");
    fOutputTree->Branch("Pt_mum_gen", fPt_mum_gen, "PDGmum_gen[NMuons_gen]/D");
    fOutputTree->Branch("Y_mum_gen", fY_mum_gen, "PDGmum_gen[NMuons_gen]/D");
    fOutputTree->Branch("Pt_gen", fPt_gen, "Pt_gen[NMuons_gen]/D");
    fOutputTree->Branch("E_gen", fE_gen, "E_gen[NMuons_gen]/D");
    fOutputTree->Branch("Px_gen", fPx_gen, "Px_gen[NMuons_gen]/D");
    fOutputTree->Branch("Py_gen", fPy_gen, "Py_gen[NMuons_gen]/D");
    fOutputTree->Branch("Pz_gen", fPz_gen, "Pz_gen[NMuons_gen]/D");
    fOutputTree->Branch("Y_gen", fY_gen, "Y_gen[NMuons_gen]/D");
    fOutputTree->Branch("Eta_gen", fEta_gen, "Eta_gen[NMuons_gen]/D");
    fOutputTree->Branch("Phi_gen", fPhi_gen, "Phi_gen[NMuons_gen]/D");
    fOutputTree->Branch("Theta_gen", fTheta_gen, "Theta_gen[NMuons_gen]/D");
    fOutputTree->Branch("Charge_gen", fCharge_gen, "Charge_gen[NMuons_gen]/I");

    fOutputTree->Branch("NHadron_gen", &fNHadron_gen, "NHadron_gen/I");
    fOutputTree->Branch("PDGHadron_gen", fPDGHadron_gen, "PDGHadron_gen[NHadron_gen]/I");
    fOutputTree->Branch("PromptHadron_gen", fPromptHadron_gen, "PromptHadron_gen[NHadron_gen]/I");
    fOutputTree->Branch("Hadron_Pt_gen", fHadron_Pt_gen, "Hadron_Pt_gen[NHadron_gen]/D");
    fOutputTree->Branch("Hadron_E_gen", fHadron_E_gen, "Hadron_E_gen[NHadron_gen]/D");
    fOutputTree->Branch("Hadron_Px_gen", fHadron_Px_gen, "Hadron_Px_gen[NHadron_gen]/D");
    fOutputTree->Branch("Hadron_Py_gen", fHadron_Py_gen, "Hadron_Py_gen[NHadron_gen]/D");
    fOutputTree->Branch("Hadron_Pz_gen", fHadron_Pz_gen, "Hadron_Pz_gen[NHadron_gen]/D");
    fOutputTree->Branch("Hadron_Y_gen", fHadron_Y_gen, "Hadron_Y_gen[NHadron_gen]/D");
    fOutputTree->Branch("Hadron_Eta_gen", fHadron_Eta_gen, "Hadron_Eta_gen[NHadron_gen]/D");
    fOutputTree->Branch("Hadron_Phi_gen", fHadron_Phi_gen, "Hadron_Phi_gen[NHadron_gen]/D");
    fOutputTree->Branch("Hadron_Theta_gen", fHadron_Theta_gen, "Hadron_Theta_gen[NHadron_gen]/D");

    TH1D *h_charm = new TH1D("h_charm", "charm number", 10, -0.5, 9.5);
    TH1D *h_bbbar = new TH1D("h_bbbar", "bbbar number", 10, -0.5, 9.5);

    // Definzione contatori
    int charm_number = 0;
    int prompt_charm_number = 0;
    int beauty_number = 0;

    // Def contatori muoni
    int muplus_number = 0;
    int muminus_number = 0;

    int muplus_number_charm = 0;
    int muplus_number_beauty = 0;

    int muminus_number_charm = 0;
    int muminus_number_beauty = 0;

    // taking track of time spent
    TStopwatch t;
    t.Start();

    // Contatori per visualizzare il numero di tagli su particelle con charm
    int charm_from_bquark = 0;
    int charm_from_bmeson = 0;
    int charm_from_bbarion = 0;

    int muminus_bquark = 0;
    int muminus_bmeson = 0;
    int muminus_bbarion = 0;

    int muplus_bquark = 0;
    int muplus_bmeson = 0;
    int muplus_bbarion = 0;

    // Inizio del loop sul numero di eventi
    for (int iEvent = 0; iEvent < nevents; ++iEvent)
    {
        if (iEvent % (10000) == 0)
        {
            Printf("Evento: %d time: CPU %f (min) real %f", iEvent, t.CpuTime(), t.RealTime());
            t.Start(kFALSE);
        }
        // useful variables
        Int_t *return_mom;

        // generation
        if (!pythia.next())
            continue;
        Int_t nHFquark_gen = 0;
        Int_t nmu_gen = 0;
        Int_t ndimu_gen = 0;

        Int_t nHFHadron_gen = 0;

        for (int i = 0; i < pythia.event.size(); i++)
        {
            TLorentzVector Particle(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
            if (TMath::Abs(pythia.event[i].id()) > 400 && TMath::Abs(pythia.event[i].id()) < 500 || TMath::Abs(pythia.event[i].id()) > 4000 && TMath::Abs(pythia.event[i].id()) < 5000)
            {

                bool isPrompt = true;
                fPDGHadron_gen[nHFHadron_gen] = pythia.event[i].id();
                fHadron_Pt_gen[nHFHadron_gen] = Particle.Pt();
                fHadron_E_gen[nHFHadron_gen] = Particle.E();
                fHadron_Px_gen[nHFHadron_gen] = Particle.Px();
                fHadron_Py_gen[nHFHadron_gen] = Particle.Py();
                fHadron_Pz_gen[nHFHadron_gen] = Particle.Pz();
                fHadron_Y_gen[nHFHadron_gen] = Particle.Rapidity();
                fHadron_Eta_gen[nHFHadron_gen] = Particle.Eta();
                fHadron_Phi_gen[nHFHadron_gen] = Particle.Phi();
                fHadron_Theta_gen[nHFHadron_gen] = Particle.Theta();

                isPrompt = isCharm(i, iEvent, false);

                if (isPrompt == true)
                    fPromptHadron_gen[nHFHadron_gen] = 4;
                else
                    fPromptHadron_gen[nHFHadron_gen] = 5;

                nHFHadron_gen++;
            }
            if (TMath::Abs(pythia.event[i].id()) > 500 && TMath::Abs(pythia.event[i].id()) < 600 || TMath::Abs(pythia.event[i].id()) > 5000 && TMath::Abs(pythia.event[i].id()) < 6000)
            {

                bool isPrompt = true;
                fPDGHadron_gen[nHFHadron_gen] = pythia.event[i].id();
                fPromptHadron_gen[nHFHadron_gen] = 5;
                fHadron_Pt_gen[nHFHadron_gen] = Particle.Pt();
                fHadron_E_gen[nHFHadron_gen] = Particle.E();
                fHadron_Px_gen[nHFHadron_gen] = Particle.Px();
                fHadron_Py_gen[nHFHadron_gen] = Particle.Py();
                fHadron_Pz_gen[nHFHadron_gen] = Particle.Pz();
                fHadron_Y_gen[nHFHadron_gen] = Particle.Rapidity();
                fHadron_Eta_gen[nHFHadron_gen] = Particle.Eta();
                fHadron_Phi_gen[nHFHadron_gen] = Particle.Phi();
                fHadron_Theta_gen[nHFHadron_gen] = Particle.Theta();

                nHFHadron_gen++;
            }
            if (Particle.Rapidity() < -4.0 || Particle.Rapidity() > -2.5)
                continue;

            if (TMath::Abs(pythia.event[i].id()) == 4 || TMath::Abs(pythia.event[i].id()) == 5)
            {
                Int_t index_daughter1 = pythia.event[i].daughter1();
                Int_t index_daughter2 = pythia.event[i].daughter2();
                // pythia.event.list();

                Int_t pdg_daughter1 = pythia.event[index_daughter1].id();
                Int_t pdg_daughter2 = pythia.event[index_daughter2].id();

                if ((TMath::Abs(pdg_daughter1) > 400 && TMath::Abs(pdg_daughter1) < 600) || (TMath::Abs(pdg_daughter1) > 4000 && TMath::Abs(pdg_daughter1) < 6000) || (TMath::Abs(pdg_daughter2) > 400 && TMath::Abs(pdg_daughter2) < 600) || (TMath::Abs(pdg_daughter2) > 4000 && TMath::Abs(pdg_daughter2) < 6000))
                {
                    // printf("Array index %d || PDG %d || PT %f || Y %f \n", i, pythia.event[i].id(), Particle.Pt(), Particle.Rapidity());

                    // printf("Array index %d || PDG %d || Array index FIGLIA1 %d || PDG FIGLIA1 %d || Array index FIGLIA2 %d || PDG FIGLIA2 %d \n", i, pythia.event[i].id(), index_daughter1, pythia.event[index_daughter1].id(), index_daughter2, pythia.event[index_daughter2].id());
                    fPDG_HFquark_gen_daughter1[nHFquark_gen] = pdg_daughter1;
                    fPDG_HFquark_gen_daughter2[nHFquark_gen] = pdg_daughter2;

                    fPDG_HFquark_gen[nHFquark_gen] = pythia.event[i].id();
                    fPt_HFquark_gen[nHFquark_gen] = Particle.Pt();
                    fY_HFquark_gen[nHFquark_gen] = Particle.Rapidity();

                    nHFquark_gen++;
                }
            }

            // Ricerca singolo mu meno
            if (!(TMath::Abs(pythia.event[i].id()) == 13))
                continue;

            TLorentzVector Muon_Mum1(pythia.event[pythia.event[i].mother1()].px(), pythia.event[pythia.event[i].mother1()].py(), pythia.event[pythia.event[i].mother1()].pz(), pythia.event[pythia.event[i].mother1()].e());
            TLorentzVector Muon_Mum2(pythia.event[pythia.event[i].mother2()].px(), pythia.event[pythia.event[i].mother2()].py(), pythia.event[pythia.event[i].mother2()].pz(), pythia.event[pythia.event[i].mother2()].e());
            bool isMuCharm = false;
            bool isMuBeauty = false;
            Int_t PDG_Mu_Mum1=pythia.event[pythia.event[i].mother1()].id();
            Int_t PDG_Mu_Mum2=pythia.event[pythia.event[i].mother2()].id();

            if ((TMath::Abs(PDG_Mu_Mum1) > 400 && TMath::Abs(PDG_Mu_Mum1) < 500) || (TMath::Abs(PDG_Mu_Mum1) > 4000 && TMath::Abs(PDG_Mu_Mum1) < 5000))
            {

                isMuCharm = isCharm(pythia.event[i].mother1(), iEvent, kFALSE);
                if (isMuCharm)
                    fPromptmum_gen[nmu_gen] = 4;
                else
                    fPromptmum_gen[nmu_gen] = 5;

                fPDGmum_gen[nmu_gen] = PDG_Mu_Mum1;
                fPt_mum_gen[nmu_gen] = Muon_Mum1.Pt();
                fY_mum_gen[nmu_gen] = Muon_Mum1.Rapidity();
            }
            else if ( (TMath::Abs(PDG_Mu_Mum2) > 400 && TMath::Abs(PDG_Mu_Mum2) < 500) || (TMath::Abs(PDG_Mu_Mum2) > 4000 && TMath::Abs(PDG_Mu_Mum2) < 5000))
            {

                isMuCharm = isCharm(pythia.event[i].mother2(), iEvent, kFALSE);
                if (isMuCharm)
                    fPromptmum_gen[nmu_gen] = 4;
                else
                    fPromptmum_gen[nmu_gen] = 5;

                fPDGmum_gen[nmu_gen] = PDG_Mu_Mum2;
                fPt_mum_gen[nmu_gen] = Muon_Mum2.Pt();
                fY_mum_gen[nmu_gen] = Muon_Mum2.Rapidity();
            }
            else if ((TMath::Abs(PDG_Mu_Mum1) > 500 && TMath::Abs(PDG_Mu_Mum1) < 600) || (TMath::Abs(PDG_Mu_Mum1) > 5000 && TMath::Abs(PDG_Mu_Mum1) < 6000))
            {

                isMuBeauty = true;
                fPDGmum_gen[nmu_gen] = PDG_Mu_Mum1;
                fPt_mum_gen[nmu_gen] = Muon_Mum1.Pt();
                fY_mum_gen[nmu_gen] = Muon_Mum1.Rapidity();
            }
            else if ((TMath::Abs(PDG_Mu_Mum2) > 500 && TMath::Abs(PDG_Mu_Mum2) < 600) || (TMath::Abs(PDG_Mu_Mum2) > 5000 && TMath::Abs(PDG_Mu_Mum2) < 6000))
            {

                isMuBeauty = true;
                fPDGmum_gen[nmu_gen] = PDG_Mu_Mum2;
                fPt_mum_gen[nmu_gen] = Muon_Mum2.Pt();
                fY_mum_gen[nmu_gen] = Muon_Mum2.Rapidity();
            }
            else
                continue;

            if (isMuCharm || isMuBeauty)
            {
                fPDG_gen[nmu_gen] = pythia.event[i].id();
                fPt_gen[nmu_gen] = Particle.Pt();
                fE_gen[nmu_gen] = Particle.E();
                fPx_gen[nmu_gen] = Particle.Px();
                fPy_gen[nmu_gen] = Particle.Py();
                fPz_gen[nmu_gen] = Particle.Pz();
                fY_gen[nmu_gen] = Particle.Rapidity();
                fEta_gen[nmu_gen] = Particle.Eta();
                fPhi_gen[nmu_gen] = Particle.Phi();
                fTheta_gen[nmu_gen] = Particle.Theta();
                fCharge_gen[nmu_gen] = -TMath::Sign(1, pythia.event[i].id());
                nmu_gen++;
            }
        }
        fN_HFquarks_gen = nHFquark_gen;
        fNMuons_gen = nmu_gen;
        fNHadron_gen = nHFHadron_gen;
        fOutputTree->Fill();
    }

    // accessing stat info and writing output

    pythia.stat();
    // cout << "Number of total mu meno: " << muminus_number << endl;
    // cout << "Number of total mu meno prompt: " << muminus_number_charm << endl;
    // cout << "Number of mu minus from BQuark" << muminus_bquark << endl;
    // cout << "Number of mu minus from BMeson" << muminus_bmeson << endl;
    // cout << "Number of mu minus from BBarion" << muminus_bbarion << endl;
    // cout << "Number of mu meno from Beauty" << muminus_number_beauty << endl;
    // cout << "Number of total mu più: " << muplus_number << endl;
    // cout << "Number of mu plus from Charm Quark" << muplus_bquark << endl;
    // cout << "Number of mu plus from Charm Meson" << muplus_bmeson << endl;
    // cout << "Number of mu plus from Charm Barion" << muplus_bbarion << endl;
    // cout << "Number of mu plus from Beauty" << muplus_number_beauty << endl;
    // cout << "Tagli Particle quark: " << charm_from_bquark << endl;
    // cout << "Tagli Particle meson: " << charm_from_bmeson << endl;
    // cout << "Tagli Particle barion: " << charm_from_bbarion << endl;
    // cout << "Tagli TOTALI: " << charm_from_bquark + charm_from_bmeson + charm_from_bbarion << endl;
    // cout << "Number of total Charm Particle: " << charm_number << endl;
    // cout << "Number of total Beauty Particle: " << beauty_number << endl;
    // cout << "Number of prompt Prompt Charm Particle: " << prompt_charm_number << endl;

    // printf("Charm mode1 %d,Charm mode2 %d\n", charm_mode1, charm_mode2);
    t.Start(kFALSE);
    printf("Tempo di esecuzione: %0.01f\n", t.RealTime());
    printf("Salvataggio Completo in %s\n", outfilename1);

    outFile1->cd();
    h_charm->Write(0, 2, 0);
    h_bbbar->Write(0, 2, 0);
    fOutputTree->Write(0, 2, 0);

    outFile1->Close();
    // Done.
    return 0;
}

bool isCharm(int i, int iEvent, bool verbose)
{
    int *mothers = NULL;
    bool isPrompt = true;
    mothers = new int[200];
    for (int i = 0; i < 200; i++)
    {
        mothers[i] = 0;
    }

    vector<int> charm_ancestor;
    charm_ancestor.clear();
    charm_ancestor = pythia.event[i].motherList();
    bool testing = true;
    int m1 = i;
    int m2 = i;
    int mom1 = pythia.event[m1].mother1();
    int mom2 = pythia.event[m1].mother2();
    int momindex = 0;
    int giro = 1;
    mothers[momindex] = mom1;
    mothers[momindex + 1] = mom2;

    // Inizio il test sulla particella charmata

    while (testing)
    {

        int temp = momindex + giro + 1;

        if (mothers[momindex] <= 2 && mothers[momindex + 1] <= 2)
        {
            if (mothers[momindex + 2] <= 2)
            {
                if (mothers[momindex + 4] <= 2)
                {
                    testing = false;
                    break;
                }
                else
                    momindex = momindex + 4;
            }
            else
            {
                momindex = momindex + 2;
                testing = true;
            }
        }
        else if (mothers[momindex] == mothers[momindex + 1])
        {
            mothers[momindex] = mothers[momindex];
            mothers[momindex + 1] = 0;
        }
        else if (mothers[momindex + 1] == 0)
        {
            if (TMath::Abs(pythia.event[mothers[momindex]].id()) == 5)
            {
                if (verbose)
                {
                    pythia.event.list();
                    printf("Charm from bquark: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n", iEvent, i, pythia.event[i].id(), mothers[momindex], pythia.event[mothers[momindex]].id(), mothers[momindex + 1], pythia.event[mothers[momindex + 1]].id());
                }
                isPrompt = false;
                testing = false;
                break;
            }
            else if (TMath::Abs(pythia.event[mothers[momindex]].id()) > 500 && TMath::Abs(pythia.event[mothers[momindex]].id()) < 600)
            {
                if (verbose)
                {
                    pythia.event.list();
                    printf("Charm from bmeson: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n", iEvent, i, pythia.event[i].id(), mothers[momindex], pythia.event[mothers[momindex]].id(), mothers[momindex + 1], pythia.event[mothers[momindex + 1]].id());
                }
                testing = false;
                isPrompt = false;
                break;
            }
            else if (TMath::Abs(pythia.event[mothers[momindex]].id()) > 5000 && TMath::Abs(pythia.event[mothers[momindex]].id()) < 6000)
            {
                if (verbose)
                {
                    pythia.event.list();
                    printf("Charm from bbarion: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n", iEvent, i, pythia.event[i].id(), mothers[momindex], pythia.event[mothers[momindex]].id(), mothers[momindex + 1], pythia.event[mothers[momindex + 1]].id());
                }
                testing = false;
                isPrompt = false;
                break;
            }
            else
            {
                charm_ancestor.clear();
                charm_ancestor = pythia.event[mothers[momindex]].motherList();
                mothers[momindex] = pythia.event[mothers[momindex]].mother1();
                mothers[momindex + 1] = pythia.event[mothers[momindex]].mother2();
                if (charm_ancestor.size() > 1)
                {
                    mothers[momindex] = charm_ancestor[0];
                    mothers[momindex + 1] = charm_ancestor[charm_ancestor.size() - 1];
                }
                else
                {
                    mothers[momindex] = charm_ancestor[0];
                    mothers[momindex + 1] = 0;
                }
                testing = true;
            }
        }
        else
        {
            if (TMath::Abs(pythia.event[mothers[momindex]].id()) == 5 || TMath::Abs(pythia.event[mothers[momindex + 1]].id()) == 5)
            {
                if (verbose)
                {
                    pythia.event.list();
                    printf("Charm from bquark: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n", iEvent, i, pythia.event[i].id(), mothers[momindex], pythia.event[mothers[momindex]].id(), mothers[momindex + 1], pythia.event[mothers[momindex + 1]].id());
                }

                testing = false;
                isPrompt = false;
                break;
            }
            else if (TMath::Abs(pythia.event[mothers[momindex]].id()) > 500 && TMath::Abs(pythia.event[mothers[momindex]].id()) < 600 || TMath::Abs(pythia.event[mothers[momindex + 1]].id()) > 500 && TMath::Abs(pythia.event[mothers[momindex + 1]].id()) < 600)
            {
                if (verbose)
                {
                    pythia.event.list();
                    printf("Charm from bmeson: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n", iEvent, i, pythia.event[i].id(), mothers[momindex], pythia.event[mothers[momindex]].id(), mothers[momindex + 1], pythia.event[mothers[momindex + 1]].id());
                }
                testing = false;
                isPrompt = false;
                break;
            }
            else if (TMath::Abs(pythia.event[mothers[momindex]].id()) > 5000 && TMath::Abs(pythia.event[mothers[momindex]].id()) < 6000 || TMath::Abs(pythia.event[mothers[momindex + 1]].id()) > 5000 && TMath::Abs(pythia.event[mothers[momindex + 1]].id()) < 6000)
            {
                if (verbose)
                {
                    pythia.event.list();
                    printf("Charm from bbarion: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n", iEvent, i, pythia.event[i].id(), mothers[momindex], pythia.event[mothers[momindex]].id(), mothers[momindex + 1], pythia.event[mothers[momindex + 1]].id());
                }
                isPrompt = false;
                break;
            }
            else
            {
                mothers[temp] = pythia.event[mothers[momindex]].mother1();
                mothers[temp + 1] = pythia.event[mothers[momindex]].mother2();

                mothers[temp + 2] = pythia.event[mothers[momindex + 1]].mother1();
                mothers[temp + 3] = pythia.event[mothers[momindex + 1]].mother2();
                momindex = momindex + 2;
                giro = giro + 2;
                testing = true;
            }
        }
    }

    delete[] mothers;

    return isPrompt;
}

bool IsChargedPhysicalPrimary(Particle part)
{

    if (part.isFinal() && part.isCharged() && part.tau0() > 10)
    {
        if (part.status() / 10 == 8)
        {
            return kTRUE;
        }
        else if (TMath::Abs(part.status()) == 91)
        {
            int indmoth = part.mother1();
            if (pythia.event[indmoth].tau0() < 10)
            {
                return kTRUE;
            }
        }
    }
    return kFALSE;
}

bool IsPIDPhysicalPrimary(Particle part)
{

    if (part.status() / 10 == 8 && part.id() != 310)
        return kTRUE;
    else
    {
        int indmoth = part.mother1();
        if (pythia.event[indmoth].tau0() < 10)
            return kTRUE;
    }
    return kFALSE;
}
void ConfigHeavyFlavor()
{
    //
    // Default configuration for Heavy Flavor production
    //
    // All QCD processes
    //

    // No multiple interactions
    pythia.readString("PartonLevel:MPI = off");
    // multiple parton interaction
    pythia.readString("MultipartonInteractions:pTmin = 0.0");
    pythia.readString("MultipartonInteractions:pT0Ref = 0.0");

    // Initial/final parton shower on (Pythia default)
    pythia.readString("PartonLevel:ISR = on");
    pythia.readString("PartonLevel:FSR = on");

    // 2nd order alpha_s
    pythia.readString("SigmaProcess:alphaSorder = 2");

    // QCD scales processo di prudzione coppia ccbar + costante alphas al secondo ordine
    pythia.readString("SigmaProcess:renormScale2 = 2");
    pythia.readString("SigmaProcess:renormMultFac = 1.");
}
void AtlasTuning()
{
    //
    // Configuration for the ATLAS tuning
    //    ReadString(Form("PDF:LHAPDFset = %s", AliStructFuncType::PDFsetName(kCTEQ5L).Data()));
    pythia.readString("PartonLevel:MPI = on");
    pythia.readString("MultipartonInteractions:pTmin = 1.9");
    pythia.readString("MultipartonInteractions:pT0Ref = 1.8");
    pythia.readString("MultipartonInteractions:ecmRef = 1000.");
    pythia.readString("MultipartonInteractions:expPow = 0.16");
    pythia.readString("MultipartonInteractions:bProfile = 2");
    pythia.readString("MultipartonInteractions:coreFraction = 0.16");
    pythia.readString("MultipartonInteractions:coreRadius = 0.5");
    //    SetPARP(85,0.33);          // Regulates gluon prod. mechanism
    //    SetPARP(86,0.66);          // Regulates gluon prod. mechanism
    pythia.readString("SigmaProcess:factorMultFac = 1.");
}

// problemi nell'usare runDY.sh
