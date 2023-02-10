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
//Versione aggiornata al 3/11/21
//scambio-> oggetti mu- contengono i mu+ e viceversa
//Salvataggio in due root file, il primo completo mentre il secondo contiene solo l'info relativi ai muoni
using namespace Pythia8;
Pythia pythia;

bool isCharm(int,int ,int &,int &,int &,bool);
bool IsChargedPhysicalPrimary(Particle);
bool IsPIDPhysicalPrimary(Particle);
void ConfigHeavyFlavor();
void AtlasTuning();

int main(int argc, char* argv[]) {
    char *selectedprocess[5]={"SoftQCD","ccbar","bbbar","SoftQCDall","Drell-Yan"};
    char *change[4]={"Def","Config","Config+Atlas","Mode2"};
    char *selectedBR[2]={"DefaultBR","ForceBR"};
    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);
    //printf("argc = %d\n argv = %s\n",theApp.Argc(),theApp.Argv(0));
    int nevents=100000;
    int seed=2710;
    int mode=-1;    //-1==Monash  .  1 not allowed. Check JHEP08(2015)003
    int n_MPI=-1;    //-1==we take everything --> >0--> register only larger nMPI events
    int chooseprocess=1;
    int choosechange=0;    //1 = SoftQCD, 2 = ccbar, 3 = bbbar
    int BR=0;    //0 = original BR, 1 = force semileptonic decay of charm
    
    //read arguments to change settings
    for(int ig=0;ig<argc;++ig){
        TString str=argv[ig];
        printf("str %d %s\n",ig,argv[ig]);
        if(str.Contains("--events")){
            sscanf(argv[ig+1],"%d",&nevents);
            printf("nevents after setting: %d\n",nevents);
        }
        if(str.Contains("--seed")){
            sscanf(argv[ig+1],"%d",&seed);
            printf("seed after setting: %d\n",seed);
        }
        if(str.Contains("--mode")){
            sscanf(argv[ig+1],"%d",&mode);
            printf("mode after setting: %d\n",mode);
            if(mode !=-1 && mode!=0 && mode!=2 && mode!=3) {printf("mode must be -1, 0, 2 or 3. Try again\n"); return 0;}
        }
        if(str.Contains("--triggerMPI")){
            sscanf(argv[ig+1],"%d",&n_MPI);
            printf("we trigger on %d MPI %d\n",n_MPI);
        }
        if(str.Contains("--process")){
            sscanf(argv[ig+1],"%d",&chooseprocess);
//            printf("generation type %d\n",chooseprocess);
            printf("process after setting %s\n",selectedprocess[chooseprocess-1]);
        }
        if(str.Contains("--choosechange")){
            sscanf(argv[ig+1],"%d",&choosechange);
//            printf("generation type %d\n",choosechange);
            printf("mode proceess after setting %s\n",change[choosechange]);
        }
        if(str.Contains("--BR")){
            sscanf(argv[ig+1],"%d",&BR);
//            printf("force semileptonic D %d\n",BR);
            printf("BR after setting %s\n",selectedBR[BR]);
        }
        
    }
    
    Settings& settings=pythia.settings;
    cout<<"PROBLEMA"<<endl;
    // Basic settings (collision energ, processes....)
    
    if(chooseprocess==1)
    {
        if (choosechange==0) {
            cout<<"choosechange 0"<<endl;
            pythia.readString("SoftQCD:nonDiffractive = on ");
        }
        else if (choosechange==1){
            cout<<"ConfigHeavyFlavor"<<endl;
            pythia.readString("SoftQCD:nonDiffractive = on ");
            
            //ConfigHeavyFlavor
            ConfigHeavyFlavor();
            
            //Mix Set
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
        else if(choosechange==2){
            cout<<"ConfigHeavyFlavor+AtlasTuning"<<endl;
            pythia.readString("SoftQCD:nonDiffractive = on ");
            ConfigHeavyFlavor();
            //MIX SET
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
        else if(choosechange==3){
            cout<<"Mode2"<<endl;
            pythia.readString("SoftQCD:nonDiffractive = on");
            pythia.readString(Form("Tune:pp = %d",14));                       //14: Monash 2013
            pythia.readString(Form("Beams:eCM = %f",(float)13000));
            // initialize random seed
            pythia.readString("Random:setSeed = on");
            pythia.readString(Form("Random:seed %d",seed));
            
            pythia.readString("ColourReconnection:mode = 1");
            cout<<"ColourReconnection:mode = 1"<<endl;
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
            cout<<"BeamRemnants:remnantMode = 1"<<endl;
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
    else if (chooseprocess==2) pythia.readString(" HardQCD:hardccbar = on ");
    else if (chooseprocess==3) pythia.readString(" HardQCD:hardbbbar = on ");//forza ad avere almeno una coppia di b e bbar prodotti
    else if (chooseprocess==4) pythia.readString(" SoftQCD:all = on ");
    else if (chooseprocess==5){
        pythia.readString("Beams:idA = 2212");            //       ! first beam, p = 2212, pbar = -2212
        pythia.readString("Beams:idB = 2212");            //       ! second beam, p = 2212, pbar = -2212
        pythia.readString("Beams:eCM = 13000.");          //       ! CM energy of collision
        pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
        pythia.readString("WeakZ0:gmZmode = 1 ");
        pythia.readString("PhaseSpace:mHatMin = 4.");
    }
    else {
        printf("Invalid process!\n");
        exit(0);
    }
    //forza ad avere decadimento buono in coppie di ccbar e bbar
    char outfilename1[100];
    char outfilename2[100];

    // const char *path="/Volumes/REPOSITORY/Test/";

    const char *path="~/pythia_stand_prov";
    sprintf(outfilename1,"%sAll_HFtest_Sim_%s_%s_%d_%d_%s.root",path,selectedprocess[chooseprocess-1],change[choosechange],nevents,seed,selectedBR[BR]);
    sprintf(outfilename2,"%sMuon_HFtest_Sim_%s_%s_%d_%d_%s.root",path,selectedprocess[chooseprocess-1],change[choosechange],nevents,seed,selectedBR[BR]);
    
    if (choosechange !=3){
        pythia.readString(Form("Tune:pp = %d",14));                       //14: Monash 2013
        pythia.readString(Form("Beams:eCM = %f",(float)13000));
        // initialize random seed
        pythia.readString("Random:setSeed = on");
        pythia.readString(Form("Random:seed %d",seed));
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
    
    if(chooseprocess==5 && BR==1){
        pythia.readString("23:onMode = 0");
        pythia.readString("23:onIfAny = 13");
        pythia.particleData.list(23);
    }
    if(chooseprocess !=5 && BR==1){
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
    
    TH1D *h_charm = new TH1D("h_charm","charm number", 10, -0.5, 9.5);
    TH1D *h_bbbar = new TH1D("h_bbbar","bbbar number", 10, -0.5, 9.5);
    
    //Definzione contatori
    int charm_number=0;
    int prompt_charm_number=0;
    int beauty_number=0;


    //Def contatori muoni
    int muplus_number=0;
    int muminus_number=0;
    
    int muplus_number_charm=0;
    int muplus_number_beauty=0;

    int muminus_number_charm=0;
    int muminus_number_beauty=0;

    //taking track of time spent
    TStopwatch t;
    t.Start();
    
    //Definizione oggetti tree
    
    TFile* outFile1 = new TFile(outfilename1,"RECREATE");
    TFile* outFile2 = new TFile(outfilename2,"RECREATE");
    outFile1->cd();
    TTree *tree = new TTree("T1","Particle Data");
    outFile2->cd();
    TTree *tree2 = new TTree("T2","Muon Data");
    
    //Definizione TClonesArray per salvare la cinematica delle diverse particelle e vettori per salvare il PDG Code
    
    TClonesArray *charm_particle = new TClonesArray("TLorentzVector",100);
    TClonesArray& CHARM_PARTICLE = *charm_particle;
    tree->Branch("charm_particle",&charm_particle);
    
    vector <int> charm_code;
    tree->Branch("code_charm",&charm_code);
    
    TClonesArray *beauty_particle = new TClonesArray("TLorentzVector",100);
    TClonesArray& BEAUTY_PARTICLE = *beauty_particle;
    tree->Branch("beauty_particle",&beauty_particle);
    
    vector <int> beauty_code;
    tree->Branch("code_beauty",&beauty_code);
    
    TClonesArray *prompt_charm_particle = new TClonesArray("TLorentzVector",100);
    TClonesArray& PROMPT_CHARM_PARTICLE = *prompt_charm_particle;
    tree->Branch("charm_particle_prompt",&prompt_charm_particle);
    
    vector <int> charm_prompt_code;
    tree->Branch("code_charm_prompt",&charm_prompt_code);
    
    TClonesArray *beauty_charm_particle = new TClonesArray("TLorentzVector",100);
    TClonesArray& BEAUTY_CHARM_PARTICLE = *beauty_charm_particle;
    tree->Branch("charm_particle_from_beauty",&beauty_charm_particle);
    
    vector <int> charm_beauty_code;
    tree->Branch("code_charm_beauty",&charm_beauty_code);

    TClonesArray *mu_minus_charm = new TClonesArray("TLorentzVector",100);
    TClonesArray& MU_MINUS_CHARM = *mu_minus_charm;
    tree->Branch("muon_minus_charm",&mu_minus_charm);
    tree2->Branch("muon_minus_charm",&mu_minus_charm);

    TClonesArray *mu_plus_charm = new TClonesArray("TLorentzVector",100);
    TClonesArray& MU_PLUS_CHARM = *mu_plus_charm;
    tree->Branch("muon_plus_charm",&mu_plus_charm);
    tree2->Branch("muon_plus_charm",&mu_plus_charm);
        
    int mu_plus_mother_charm[15][2];
    
    for (int q=0; q<15; q++) {
        mu_plus_mother_charm[q][0]=99999;
        mu_plus_mother_charm[q][1]=99999;
    }
    
    tree->Branch("muplus_charm",mu_plus_mother_charm,"mu_plus_mother_charm[15][2]/I");
    tree2->Branch("muplus_charm",mu_plus_mother_charm,"mu_plus_mother_charm[15][2]/I");
    
    int mu_minus_mother_charm[15][2];
    
    for (int q=0; q<15; q++) {
        mu_minus_mother_charm[q][0]=99999;
        mu_minus_mother_charm[q][1]=99999;
    }
    
    tree->Branch("muminus_charm",mu_minus_mother_charm,"mu_minus_mother_charm[15][2]/I");
    tree2->Branch("muminus_charm",mu_minus_mother_charm,"mu_minus_mother_charm[15][2]/I");
    
    TClonesArray *mu_minus_beauty = new TClonesArray("TLorentzVector",100);
    TClonesArray& MU_MINUS_BEAUTY = *mu_minus_beauty;
    tree->Branch("muon_minus_beauty",&mu_minus_beauty);
    tree2->Branch("muon_minus_beauty",&mu_minus_beauty);

    TClonesArray *mu_plus_beauty = new TClonesArray("TLorentzVector",100);
    TClonesArray& MU_PLUS_BEAUTY = *mu_plus_beauty;
    tree->Branch("muon_plus_beauty",&mu_plus_beauty);
    tree2->Branch("muon_plus_beauty",&mu_plus_beauty);
    
    int mu_minus_mother_beauty[15][2];
    
    for (int q=0; q<15; q++) {
        mu_minus_mother_beauty[q][0]=99999;
        mu_minus_mother_beauty[q][1]=99999;
    }
    
    tree->Branch("muminus_beauty",mu_minus_mother_beauty,"mu_minus_mother_beauty[15][2]/I");
    tree2->Branch("muminus_beauty",mu_minus_mother_beauty,"mu_minus_mother_beauty[15][2]/I");
        
    int mu_plus_mother_beauty[15][2];
    
    for (int q=0; q<15; q++) {
        mu_plus_mother_beauty[q][0]=99999;
        mu_plus_mother_beauty[q][1]=99999;
    }
    
    tree->Branch("muplus_beauty",mu_plus_mother_beauty,"mu_plus_mother_beauty[15][2]/I");
    tree2->Branch("muplus_beauty",mu_plus_mother_beauty,"mu_plus_mother_beauty[15][2]/I");
    
    //Contatori per visualizzare il numero di tagli su particelle con charm
    int charm_from_bquark=0;
    int charm_from_bmeson=0;
    int charm_from_bbarion=0;
    
    int muminus_bquark=0;
    int muminus_bmeson=0;
    int muminus_bbarion=0;
    
    int muplus_bquark=0;
    int muplus_bmeson=0;
    int muplus_bbarion=0;
    
    int charm_mode1=0;
    int charm_mode2=0;
    
    //Vettori per salvare provvisoriamente le madri della singola particella nella selezione del prompt
    int *mothers1=NULL;
    //Vettori per salvare provvisoriamente le madri del singolo muone in modo da risalire alla provenienza
    int *mothers2=NULL;
    int *mothers3=NULL;
    
    //Vettori per salvare l'indice della singola particella e successivamente usato per riempire i TClonesArray al termine del singolo evento
    vector<int> charm_event;
    vector<int> beauty_event;
    vector<int> charm_prompt_event;
    vector<int> charm_beauty_event;
    vector<int> muplus_event;
    vector<int> muminus_event;
    vector<int> muplus_event_charm;
    vector<int> muplus_event_beauty;
    vector<int> muminus_event_charm;
    vector<int> muminus_event_beauty;
    
    // Inizio del loop sul numero di eventi
    for (int iEvent = 0; iEvent < nevents; ++iEvent) {
        if (iEvent%(10000)==0) {
            Printf("Evento: %d time: CPU %f (min) real %f",iEvent,t.CpuTime(),t.RealTime());
            t.Start(kFALSE);
        }
        //useful variables
        int n_ccbar=0, n_bbbar=0;
        
        int ncharm=0,ncharm_prompt=0,ncharm_frombeauty=0,nbeauty=0;
        int nmuplus=0, nmuminus=0;
        int nmuplus_prompt=0, nmuminus_prompt=0;
        int nmuplus_beauty=0, nmuminus_beauty=0;
        
        //generation
        if (!pythia.next()) continue;
        
        //if (n_MPI>0 && pythia.info.nMPI() < n_MPI ) {iEvent--; continue;}
        //          cout<<"Sta uscendo l'intero evento"<<endl;
                
        
        
//            if(iEvent%(1)==0){
//                Printf("Event %d time: CPU %f real %f",iEvent,t.CpuTime(),t.RealTime());
//
//
//              t.Start(kFALSE);
//            }
    
        charm_event.clear();
        beauty_event.clear();
        charm_code.clear();
        beauty_code.clear();
        charm_prompt_event.clear();
        charm_prompt_code.clear();
        charm_beauty_event.clear();
        charm_beauty_code.clear();
        muplus_event.clear();
        muminus_event.clear();
        
        muplus_event_charm.clear();
        muplus_event_beauty.clear();
        muminus_event_charm.clear();
        muminus_event_beauty.clear();
        
        bool testing;
        
        //Metodo di selezione 3.0
        //Loop sul numero di particelle generate per singolo evento
        
        for(int i=0;i<pythia.event.size();i++){
//            if (TMath::Abs(pythia.event[i].id())==4){
//                if (TMath::Abs(pythia.event[pythia.event[i].daughter1()].id())==4 || TMath::Abs(pythia.event[pythia.event[i].daughter2()].id())==4)
//                {
//                    //TMath::Abs(pythia.event[pythia.event[i].daughter1()]).id()!=4 && TMath::Abs(pythia.event[pythia.event[i].daughter2()]).id()!=4
////                    printf ("NO-> Evento: %d | Codice %d | ID %d | Figlio1 %d | Figlio2 %d\n",iEvent,i,pythia.event[i].id(),pythia.event[i].daughter1(),pythia.event[i].daughter2());
//                }else{
//                    if (TMath::Abs(pythia.event[pythia.event[i].daughter1()].id())>100) {
//                        n_ccbar++;
//                        charm_mode1=charm_mode1+n_ccbar;
//                    }
//                }
//            }
            
            //Ricerca particelle con quark charm
            if (TMath::Abs(pythia.event[i].id())==4){

                if (TMath::Abs(pythia.event[pythia.event[i].daughter1()].id())>100 || TMath::Abs(pythia.event[pythia.event[i].daughter2()].id())>100){
                    charm_mode2++;
                    bool isPrompt=true;
                    charm_number++;
                    ncharm++;
                    charm_event.push_back(i);
                    charm_code.push_back(pythia.event[i].id());
                    isPrompt=isCharm(i,iEvent,charm_from_bquark,charm_from_bmeson,charm_from_bbarion,false);

                    if (isPrompt==true) {
                        n_ccbar++;
                        prompt_charm_number++;
                        ncharm_prompt++;
                        charm_prompt_event.push_back(i);
                        charm_prompt_code.push_back(pythia.event[i].id());

                    }else {
                        ncharm_frombeauty++;
                        charm_beauty_event.push_back(i);
                        charm_beauty_code.push_back(pythia.event[i].id());
                    }
                }
            }
            //Ricerca particelle particelle charmate
            if (TMath::Abs(pythia.event[i].id())>400 && TMath::Abs(pythia.event[i].id())<500 || TMath::Abs(pythia.event[i].id())>4000 && TMath::Abs(pythia.event[i].id())<5000){

                    bool isPrompt=true;
                    charm_number++;
                    ncharm++;
                    charm_event.push_back(i);
                    charm_code.push_back(pythia.event[i].id());
                    isPrompt=isCharm(i,iEvent,charm_from_bquark,charm_from_bmeson,charm_from_bbarion,false);

                    if (isPrompt==true) {
                        prompt_charm_number++;
                        ncharm_prompt++;
                        charm_prompt_event.push_back(i);
                        charm_prompt_code.push_back(pythia.event[i].id());
                    }else {
                        ncharm_frombeauty++;
                        charm_beauty_event.push_back(i);
                        charm_beauty_code.push_back(pythia.event[i].id());
                    }

            }
            //Ricerca particelle con quark beauty
            if (TMath::Abs(pythia.event[i].id())==5){
                if (TMath::Abs(pythia.event[pythia.event[i].daughter1()].id())>100 || TMath::Abs(pythia.event[pythia.event[i].daughter2()].id())>100){

//                    printf ("bbbar Evento: %d | Codice %d | ID %d | Figlio1 %d | Figlio2 %d\n",iEvent,i,pythia.event[i].id(),pythia.event[i].daughter1(),pythia.event[i].daughter2());
                    n_bbbar++;
                    beauty_number++;
                    nbeauty++;
                    beauty_event.push_back(i);
                    beauty_code.push_back(pythia.event[i].id());

                }
            }
            //Ricerca particelle particelle con beauty
            if (TMath::Abs(pythia.event[i].id())>500 && TMath::Abs(pythia.event[i].id())<600 || TMath::Abs(pythia.event[i].id())>5000 && TMath::Abs(pythia.event[i].id())<6000) {
                beauty_number++;
                nbeauty++;
                beauty_event.push_back(i);
                beauty_code.push_back(pythia.event[i].id());
            }

            //Ricerca singolo mu meno
            if(pythia.event[i].id()==13){
                
                muminus_number++;
                nmuminus++;
                muminus_event.push_back(i);
                bool isMuCharm=false;
                bool isMuBeauty=false;
                
                if (TMath::Abs(pythia.event[pythia.event[i].mother1()].id())>400 && TMath::Abs(pythia.event[pythia.event[i].mother1()].id())<500 || TMath::Abs(pythia.event[pythia.event[i].mother1()].id())>4000 && TMath::Abs(pythia.event[pythia.event[i].mother1()].id())<5000) {
                    
                    isMuCharm=true;
                    isMuCharm=isCharm(pythia.event[i].mother1(),iEvent,muminus_bquark,muminus_bmeson,muminus_bbarion,false);
                    if(isMuCharm==false) isMuBeauty=true;
                    
                }else if (TMath::Abs(pythia.event[pythia.event[i].mother2()].id())>400 && TMath::Abs(pythia.event[pythia.event[i].mother2()].id())<500 || TMath::Abs(pythia.event[pythia.event[i].mother2()].id())>4000 && TMath::Abs(pythia.event[pythia.event[i].mother2()].id())<5000){
                    
                    isMuCharm=true;
                    isMuCharm=isCharm(pythia.event[i].mother2(),iEvent,muminus_bquark,muminus_bmeson,muminus_bbarion,false);
                    if(isMuCharm==false) isMuBeauty=true;
                    
                }else if (TMath::Abs(pythia.event[pythia.event[i].mother1()].id())>500 && TMath::Abs(pythia.event[pythia.event[i].mother1()].id())<600 || TMath::Abs(pythia.event[pythia.event[i].mother1()].id())>5000 && TMath::Abs(pythia.event[pythia.event[i].mother1()].id())<6000 || TMath::Abs(pythia.event[pythia.event[i].mother2()].id())>500 && TMath::Abs(pythia.event[pythia.event[i].mother2()].id())<600 || TMath::Abs(pythia.event[pythia.event[i].mother2()].id())>5000 && TMath::Abs(pythia.event[pythia.event[i].mother2()].id())<6000){
                    
                    isMuBeauty=true;
                }
                
                if (isMuCharm==true) {
                    
                    mu_minus_mother_charm[nmuminus_prompt][0]=pythia.event[pythia.event[i].mother1()].id();
                    mu_minus_mother_charm[nmuminus_prompt][1]=pythia.event[pythia.event[i].mother2()].id();
                    
                    nmuminus_prompt++;
                    muminus_number_charm++;
                    muminus_event_charm.push_back(i);
                    
                }else if(isMuBeauty==true){
                    if(false){
                        pythia.event.list();
                        printf ("Mu-: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n",iEvent,i,pythia.event[i].id());
                    }
                    
                    mu_minus_mother_beauty[nmuminus_beauty][0]=pythia.event[pythia.event[i].mother1()].id();
                    mu_minus_mother_beauty[nmuminus_beauty][1]=pythia.event[pythia.event[i].mother2()].id();
                    
                    nmuminus_beauty++;
                    muminus_number_beauty++;
                    muminus_event_beauty.push_back(i);
                }
            }
            //Ricerca singolo mu più
            if(pythia.event[i].id()==-13){
                
                muplus_number++;
                nmuplus++;
                muplus_event.push_back(i);
                bool isMuCharm=false;
                bool isMuBeauty=false;
                
                if (TMath::Abs(pythia.event[pythia.event[i].mother1()].id())>400 && TMath::Abs(pythia.event[pythia.event[i].mother1()].id())<500 || TMath::Abs(pythia.event[pythia.event[i].mother1()].id())>4000 && TMath::Abs(pythia.event[pythia.event[i].mother1()].id())<5000) {
                    
                    isMuCharm=true;
                    isMuCharm=isCharm(pythia.event[i].mother1(),iEvent,muplus_bquark,muplus_bmeson,muplus_bbarion,false);
                    if(isMuCharm==false) isMuBeauty=true;
                    
                }else if (TMath::Abs(pythia.event[pythia.event[i].mother2()].id())>400 && TMath::Abs(pythia.event[pythia.event[i].mother2()].id())<500 || TMath::Abs(pythia.event[pythia.event[i].mother2()].id())>4000 && TMath::Abs(pythia.event[pythia.event[i].mother2()].id())<5000){
                    
                    isMuCharm=true;
                    isMuCharm=isCharm(pythia.event[i].mother2(),iEvent,muplus_bquark,muplus_bmeson,muplus_bbarion,false);
                    if(isMuCharm==false) isMuBeauty=true;
                    
                }else if (TMath::Abs(pythia.event[pythia.event[i].mother1()].id())>500 && TMath::Abs(pythia.event[pythia.event[i].mother1()].id())<600 || TMath::Abs(pythia.event[pythia.event[i].mother1()].id())>5000 && TMath::Abs(pythia.event[pythia.event[i].mother1()].id())<6000 || TMath::Abs(pythia.event[pythia.event[i].mother2()].id())>500 && TMath::Abs(pythia.event[pythia.event[i].mother2()].id())<600 || TMath::Abs(pythia.event[pythia.event[i].mother2()].id())>5000 && TMath::Abs(pythia.event[pythia.event[i].mother2()].id())<6000){
                    
                    isMuBeauty=true;
                }
                
                if (isMuCharm==true) {
                    
                    mu_plus_mother_charm[nmuplus_prompt][0]=pythia.event[pythia.event[i].mother1()].id();
                    mu_plus_mother_charm[nmuplus_prompt][1]=pythia.event[pythia.event[i].mother2()].id();
                    
                    nmuplus_prompt++;
                    muplus_number_charm++;
                    muplus_event_charm.push_back(i);
                    
                }else if(isMuBeauty==true){
                    if(false){
                        pythia.event.list();
                        printf ("Mu+: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n",iEvent,i,pythia.event[i].id());
                    }
                    mu_plus_mother_beauty[nmuplus_beauty][0]=pythia.event[pythia.event[i].mother1()].id();
                    mu_plus_mother_beauty[nmuplus_beauty][1]=pythia.event[pythia.event[i].mother2()].id();
                    
                    nmuplus_beauty++;
                    muplus_number_beauty++;
                    muplus_event_beauty.push_back(i);
                }
            }
            
            
            
//            cout<<"Finito loop mu minus"<<endl;
        }
        
        if (n_ccbar%2!=0) {
            cout<<"n ccbar: "<<n_ccbar<<endl;
            pythia.event.list();
        }
        
        h_charm->Fill(n_ccbar);
        h_bbbar->Fill(n_bbbar);
        
        for (int i=0;i<ncharm;i++){
            new(CHARM_PARTICLE[i])TLorentzVector(pythia.event[charm_event[i]].px(),pythia.event[charm_event[i]].py(),pythia.event[charm_event[i]].pz(),pythia.event[charm_event[i]].e());
            TLorentzVector *tst=(TLorentzVector*)charm_particle->At(i);
        }
        for (int i=0;i<nbeauty;i++){
            new(BEAUTY_PARTICLE[i])TLorentzVector(pythia.event[beauty_event[i]].px(),pythia.event[beauty_event[i]].py(),pythia.event[beauty_event[i]].pz(),pythia.event[beauty_event[i]].e());
            TLorentzVector *tst_beauty=(TLorentzVector*)beauty_particle->At(i);
        }

        for (int i=0;i<ncharm_prompt;i++){
            new(PROMPT_CHARM_PARTICLE[i])TLorentzVector(pythia.event[charm_prompt_event[i]].px(),pythia.event[charm_prompt_event[i]].py(),pythia.event[charm_prompt_event[i]].pz(),pythia.event[charm_prompt_event[i]].e());
            TLorentzVector *tst_prompt=(TLorentzVector*)prompt_charm_particle->At(i);
        }
        for (int i=0;i<ncharm_frombeauty;i++){
            new(BEAUTY_CHARM_PARTICLE[i])TLorentzVector(pythia.event[charm_beauty_event[i]].px(),pythia.event[charm_beauty_event[i]].py(),pythia.event[charm_beauty_event[i]].pz(),pythia.event[charm_beauty_event[i]].e());
            TLorentzVector *tst_chb=(TLorentzVector*)beauty_charm_particle->At(i);
        }

        //Fill Branch with muon kinematic
        
        for (int i=0;i<nmuminus_prompt;i++){
            new(MU_MINUS_CHARM[i])TLorentzVector(pythia.event[muminus_event_charm[i]].px(),pythia.event[muminus_event_charm[i]].py(),pythia.event[muminus_event_charm[i]].pz(),pythia.event[muminus_event_charm[i]].e());
            TLorentzVector *tstmu_minus_prompt=(TLorentzVector*)mu_minus_charm->At(i);
        }

        for (int i=0;i<nmuplus_prompt;i++){
            new(MU_PLUS_CHARM[i])TLorentzVector(pythia.event[muplus_event_charm[i]].px(),pythia.event[muplus_event_charm[i]].py(),pythia.event[muplus_event_charm[i]].pz(),pythia.event[muplus_event_charm[i]].e());
            TLorentzVector *tstmu_plus_prompt=(TLorentzVector*)mu_plus_charm->At(i);
        }
        for (int i=0;i<nmuminus_beauty;i++){
            new(MU_MINUS_BEAUTY[i])TLorentzVector(pythia.event[muminus_event_beauty[i]].px(),pythia.event[muminus_event_beauty[i]].py(),pythia.event[muminus_event_beauty[i]].pz(),pythia.event[muminus_event_beauty[i]].e());
            TLorentzVector *tstmu_minus_beauty=(TLorentzVector*)mu_minus_beauty->At(i);
        }
        
        for (int i=0;i<nmuplus_beauty;i++){
            new(MU_PLUS_BEAUTY[i])TLorentzVector(pythia.event[muplus_event_beauty[i]].px(),pythia.event[muplus_event_beauty[i]].py(),pythia.event[muplus_event_beauty[i]].pz(),pythia.event[muplus_event_beauty[i]].e());
            TLorentzVector *tstmu_plus_beauty=(TLorentzVector*)mu_plus_beauty->At(i);
        }
      
        // Riempimento del tree e pulizia dei TClonesArray
        
        tree->Fill();
        
        tree2->Fill();
        
        charm_particle->Clear();
        beauty_particle->Clear();
        prompt_charm_particle->Clear();
        beauty_charm_particle->Clear();

        
        mu_minus_charm->Clear();
        mu_plus_charm->Clear();
        mu_minus_beauty->Clear();
        mu_plus_beauty->Clear();
        
    }
    
    
        //accessing stat info and writing output
    
    pythia.stat();
    cout<<"Number of total mu meno: "<<muminus_number<<endl;
    cout<<"Number of total mu meno prompt: "<<muminus_number_charm<<endl;
    cout<<"Number of mu minus from BQuark"<<muminus_bquark<<endl;
    cout<<"Number of mu minus from BMeson"<<muminus_bmeson<<endl;
    cout<<"Number of mu minus from BBarion"<<muminus_bbarion<<endl;
    cout<<"Number of mu meno from Beauty"<<muminus_number_beauty<<endl;
    cout<<"Number of total mu più: "<<muplus_number<<endl;
    cout<<"Number of mu plus from Charm Quark"<<muplus_bquark<<endl;
    cout<<"Number of mu plus from Charm Meson"<<muplus_bmeson<<endl;
    cout<<"Number of mu plus from Charm Barion"<<muplus_bbarion<<endl;
    cout<<"Number of mu plus from Beauty"<<muplus_number_beauty<<endl;
    cout<<"Tagli Particle quark: "<<charm_from_bquark<<endl;
    cout<<"Tagli Particle meson: "<<charm_from_bmeson<<endl;
    cout<<"Tagli Particle barion: "<<charm_from_bbarion<<endl;
    cout<<"Tagli TOTALI: "<<charm_from_bquark+charm_from_bmeson+charm_from_bbarion<<endl;
    cout<<"Number of total Charm Particle: "<<charm_number<<endl;
    cout<<"Number of total Beauty Particle: "<<beauty_number<<endl;
    cout<<"Number of prompt Prompt Charm Particle: "<<prompt_charm_number<<endl;

    printf("Charm mode1 %d,Charm mode2 %d\n",charm_mode1,charm_mode2);
    t.Start(kFALSE);
    printf("Tempo di esecuzione: %0.01f\n",t.RealTime());
    printf("Salvataggio Completo in %s\n",outfilename1);
    
    outFile1->cd();
    h_charm->Write(0,2,0);
    h_bbbar->Write(0,2,0);
    tree->Write(0,2,0);
    
    
    printf("Salvataggio Muoni in %s\n",outfilename2);
    outFile2->cd();
    tree2->Write(0,2,0);
    outFile1->Close();
    outFile2->Close();
    
        // Done.
        return 0;
}

bool isCharm(int i,int iEvent, int &charm_from_bquark,int &charm_from_bmeson, int &charm_from_bbarion, bool verbose){
    int *mothers=NULL;
    bool isPrompt=true;
    mothers=new int [200];
    for (int i=0; i<200; i++) {
            mothers[i]=0;
        }
        
        vector<int> charm_ancestor;
        charm_ancestor.clear();
        charm_ancestor=pythia.event[i].motherList();
        bool testing=true;
        int m1=i;
        int m2=i;
        int mom1=pythia.event[m1].mother1();
        int mom2=pythia.event[m1].mother2();
        int momindex=0;
        int giro=1;
        mothers[momindex]=mom1;
        mothers[momindex+1]=mom2;
        
        //Inizio il test sulla particella charmata
        
    while (testing) {

            int temp=momindex+giro+1;
            
            if (mothers[momindex]<=2 && mothers[momindex+1]<=2){
                if (mothers[momindex+2]<=2) {
                    if (mothers[momindex+4]<=2) {
                        testing=false;
                        break;
                    }
                    else momindex=momindex+4;

                } else {
                    momindex=momindex+2;
                    testing=true;
                }
            }
            else if (mothers[momindex]==mothers[momindex+1]){
                mothers[momindex]=mothers[momindex];
                mothers[momindex+1]=0;
                
            }
            else if (mothers[momindex+1]==0){
                if (TMath::Abs(pythia.event[mothers[momindex]].id())==5){
                    if(verbose){
                        pythia.event.list();
                        printf ("Charm from bquark: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n",iEvent,i,pythia.event[i].id(),mothers[momindex],pythia.event[mothers[momindex]].id(),mothers[momindex+1],pythia.event[mothers[momindex+1]].id());
                    }
                    charm_from_bquark++;
                    isPrompt=false;
                    testing=false;
                    break;
                }
                else if (TMath::Abs(pythia.event[mothers[momindex]].id())>500 && TMath::Abs(pythia.event[mothers[momindex]].id())<600){
                    if(verbose){
                        pythia.event.list();
                        printf ("Charm from bmeson: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n",iEvent,i,pythia.event[i].id(),mothers[momindex],pythia.event[mothers[momindex]].id(),mothers[momindex+1],pythia.event[mothers[momindex+1]].id());
                    }
                    charm_from_bmeson++;
                    testing=false;
                    isPrompt=false;
                    break;
                }
                else if (TMath::Abs(pythia.event[mothers[momindex]].id())>5000 && TMath::Abs(pythia.event[mothers[momindex]].id())<6000){
                    if(verbose){
                        pythia.event.list();
                        printf ("Charm from bbarion: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n",iEvent,i,pythia.event[i].id(),mothers[momindex],pythia.event[mothers[momindex]].id(),mothers[momindex+1],pythia.event[mothers[momindex+1]].id());
                    }
                    charm_from_bbarion++;
                    testing=false;
                    isPrompt=false;
                    break;
                }
                else{
                    charm_ancestor.clear();
                    charm_ancestor=pythia.event[mothers[momindex]].motherList();
                    mothers[momindex]=pythia.event[mothers[momindex]].mother1();
                    mothers[momindex+1]=pythia.event[mothers[momindex]].mother2();
                    if (charm_ancestor.size()>1) {
                        mothers[momindex]=charm_ancestor[0];
                        mothers[momindex+1]=charm_ancestor[charm_ancestor.size()-1];
                    }
                    else{
                        mothers[momindex]=charm_ancestor[0];
                        mothers[momindex+1]=0;
                    }
                    testing=true;
                }
            }
            else{
                if(TMath::Abs(pythia.event[mothers[momindex]].id())==5 || TMath::Abs(pythia.event[mothers[momindex+1]].id())==5){
                    if(verbose){
                        pythia.event.list();
                        printf ("Charm from bquark: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n",iEvent,i,pythia.event[i].id(),mothers[momindex],pythia.event[mothers[momindex]].id(),mothers[momindex+1],pythia.event[mothers[momindex+1]].id());
                    }

                    charm_from_bquark++;
                    testing=false;
                    isPrompt=false;
                    break;
                }
                else if (TMath::Abs(pythia.event[mothers[momindex]].id())>500 && TMath::Abs(pythia.event[mothers[momindex]].id())<600 || TMath::Abs(pythia.event[mothers[momindex+1]].id())>500 && TMath::Abs(pythia.event[mothers[momindex+1]].id())<600){
                    if(verbose){
                        pythia.event.list();
                        printf ("Charm from bmeson: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n",iEvent,i,pythia.event[i].id(),mothers[momindex],pythia.event[mothers[momindex]].id(),mothers[momindex+1],pythia.event[mothers[momindex+1]].id());
                    }
                    charm_from_bmeson++;
                    testing=false;
                    isPrompt=false;
                    break;
                }
                else if (TMath::Abs(pythia.event[mothers[momindex]].id())>5000 && TMath::Abs(pythia.event[mothers[momindex]].id())<6000 || TMath::Abs(pythia.event[mothers[momindex+1]].id())>5000 && TMath::Abs(pythia.event[mothers[momindex+1]].id())<6000){
                    if(verbose){
                        pythia.event.list();
                        printf ("Charm from bbarion: Evento: %d | Codice %d | ID %d | Codice Madre1 %d | Madre1 %d | Codice Madre2 %d | Madre2 %d\n",iEvent,i,pythia.event[i].id(),mothers[momindex],pythia.event[mothers[momindex]].id(),mothers[momindex+1],pythia.event[mothers[momindex+1]].id());
                    }
                    charm_from_bbarion++;
                    isPrompt=false;
                    break;
                }
                else {
                    mothers[temp]=pythia.event[mothers[momindex]].mother1();
                    mothers[temp+1]=pythia.event[mothers[momindex]].mother2();
                    
                    mothers[temp+2]=pythia.event[mothers[momindex+1]].mother1();
                    mothers[temp+3]=pythia.event[mothers[momindex+1]].mother2();
                    momindex=momindex+2;
                    giro=giro+2;
                    testing=true;
                    
                }
            }
        }
    
    delete [] mothers;
    
    return isPrompt;
}

bool IsChargedPhysicalPrimary(Particle part){
    
    if(part.isFinal() && part.isCharged()&&part.tau0()>10){
        if(part.status()/10==8){
            return kTRUE;
        }
        else if(TMath::Abs(part.status())==91){
            int indmoth=part.mother1();
            if(pythia.event[indmoth].tau0()<10){
                return kTRUE;
            }
        }
    }
    return kFALSE;
    
}


bool IsPIDPhysicalPrimary(Particle part){
    
    if( part.status()/10==8 && part.id()!=310 ) return kTRUE;
    else {
        int indmoth=part.mother1();
        if(pythia.event[indmoth].tau0()<10) return kTRUE;
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
    //multiple parton interaction
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

//problemi nell'usare runDY.sh
