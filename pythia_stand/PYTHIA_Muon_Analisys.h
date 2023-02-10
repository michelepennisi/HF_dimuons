
#ifndef PYTHIA_MUON_ANALISYS_H

#define PYTHIA_MUON_ANALISYS_H

//Macro riferita al file di simulazione DY_sigle_new.cc. Consente di analizzare le sim prodotte con questo + il confronto con il vecchio sistema di simulazione.

TH1D *h_charm = new TH1D("h_charm","c#bar{c} pair", 10, -0.5, 9.5);
TH1D *h_beauty = new TH1D("h_beauty","b#bar{b} pair", 10, -0.5, 9.5);

static const int NBINS = 43;
Double_t edges[NBINS + 1] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,10.0,11.5,13.0,14.5,16.0,18.5,20.0,22.5,25.0,28,31};

//Hist to store charm particle
TH1D *h_ptccbar =new TH1D("h_ptccbar","h_ptccbar",NBINS,edges);

TH1D *h_yccbar =new TH1D("h_yccbar","h_yccbar",100,-9.5,9.5);


TH1D *h_ptAllCharm =new TH1D("h_ptAllCharm","h_ptAllCharm",NBINS,edges);
TH1D *h_yAllCharm =new TH1D("h_yAllCharm","h_yAllCharm",100,-9.5,9.5);

TH1D *h_ptCharmMeson =new TH1D("h_ptCharmMeson","h_ptCharmMeson",NBINS,edges);
TH1D *h_ptCharmBarion =new TH1D("h_ptCharmBarion","h_ptCharmBarion",NBINS,edges);

TH1D *h_yCharmMeson =new TH1D("h_yCharmMeson","h_yCharmMeson",100,-9.5,9.5);
TH1D *h_yCharmBarion =new TH1D("h_yCharmBarion","h_yCharmBarion",100,-9.5,9.5);

//Hist to store beauty particle
TH1D *h_ptbbbar =new TH1D("h_ptbbbar","h_ptbbbar",NBINS,edges);

TH1D *h_ybbbar =new TH1D("h_ybbbar","h_ybbbar",100,-9.5,9.5);

TH1D *h_ptAllBeauty =new TH1D("h_ptAllBeauty","h_ptAllBeauty",NBINS,edges);
TH1D *h_yAllBeauty =new TH1D("h_yAllBeauty","h_yAllBeauty",100,-9.5,9.5);

TH1D *h_ptBeautyMeson =new TH1D("h_ptBeautyMeson","h_ptBeautyMeson",NBINS,edges);
TH1D *h_ptBzero =new TH1D("h_ptBzero","h_ptBzero",NBINS,edges);
TH1D *h_ptBplus =new TH1D("h_ptBplus","h_ptBplus",NBINS,edges);
TH1D *h_ptBstrangeplus =new TH1D("h_ptBstrangeplus","h_ptBstrangeplus",NBINS,edges);
TH1D *h_ptBeautyMeson_rest =new TH1D("h_ptBeautyMeson_rest","h_ptBeautyMeson_rest",NBINS,edges);

TH1D *h_ptBeautyBarion =new TH1D("h_ptBeautyBarion","h_ptBeautyBarion",NBINS,edges);
TH1D *h_ptLambdabeautyplus =new TH1D("h_ptLambdabeautyplus","h_ptLambdabeautyplus",NBINS,edges);
TH1D *h_ptSigmabeautyzero =new TH1D("h_ptSigmabeautyzero","h_ptSigmabeautyzero",NBINS,edges);
TH1D *h_ptSigmabeautyplus =new TH1D("h_ptSigmabeautyplus","h_ptSigmabeautyplus",NBINS,edges);
TH1D *h_ptBeautyBarion_rest =new TH1D("h_ptBeautyBarion_rest","h_ptBeautyBarion_rest",NBINS,edges);

TH1D *h_yBeautyMeson =new TH1D("h_yBeautyMeson","h_yBeautyMeson",100,-9.5,9.5);
TH1D *h_yBzero =new TH1D("h_yBzero","h_yBzero",100,-9.5,9.5);
TH1D *h_yBplus =new TH1D("h_yBplus","h_yBplus",100,-9.5,9.5);
TH1D *h_yBstrangeplus =new TH1D("h_yBstrangeplus","h_yBstrangeplus",100,-9.5,9.5);
TH1D *h_yBeautyMeson_rest =new TH1D("h_yBeautyMeson_rest","h_yBeautyMeson_rest",100,-9.5,9.5);

TH1D *h_yBeautyBarion =new TH1D("h_yBeautyBarion","h_yBeautyBarion",100,-9.5,9.5);
TH1D *h_yLambdabeautyplus =new TH1D("h_yLambdabeautyplus","h_yLambdabeautyplus",100,-9.5,9.5);
TH1D *h_ySigmabeautyzero =new TH1D("h_ySigmabeautyzero","h_ySigmabeautyzero",100,-9.5,9.5);
TH1D *h_ySigmabeautyplus =new TH1D("h_ySigmabeautyplus","h_ySigmabeautyplus",100,-9.5,9.5);
TH1D *h_yBeautyBarion_rest =new TH1D("h_yBeautyBarion_rest","h_yBeautyBarion_rest",100,-9.5,9.5);

TH1D *h_ptAntiBeauty =new TH1D("h_ptAntiBeauty","h_ptAntiBeauty",NBINS,edges);
TH1D *h_ptBeautyAntiMeson =new TH1D("h_ptBeautyAntiMeson","h_ptBeautyAntiMeson",NBINS,edges);
TH1D *h_ptBeautyAntiBarion =new TH1D("h_ptBeautyAntiBarion","h_ptBeautyAntiBarion",NBINS,edges);

TH1D *h_yAntiBeauty =new TH1D("h_yAntiBeauty","h_yAntiBeauty",100,-9.5,9.5);
TH1D *h_yBeautyAntiMeson =new TH1D("h_yBeautyAntiMeson","h_yBeautyAntiMeson",100,-9.5,9.5);
TH1D *h_yBeautyAntiBarion =new TH1D("h_yBeautyAntiBarion","h_yBeautyAntiBarion",100,-9.5,9.5);

//Hist to store prompt charm particle
TH1D *h_ptAllCharm_prompt =new TH1D("h_ptAllCharm_prompt","h_ptAllCharm_prompt",NBINS,edges);
TH1D *h_yAllCharm_prompt =new TH1D("h_yAllCharm_prompt","h_yAllCharm_prompt",100,-9.5,9.5);

TH1D *h_ptAllAntiCharm_prompt =new TH1D("h_ptAllAntiCharm_prompt","h_ptAllAntiCharm_prompt",NBINS,edges);
TH1D *h_yAllAntiCharm_prompt =new TH1D("h_yAllAntiCharm_prompt","h_yAllAntiCharm_prompt",100,-9.5,9.5);

TH1D *h_ptCharmMeson_prompt =new TH1D("h_ptCharmMeson_prompt","h_ptCharmMeson_prompt",NBINS,edges);
TH1D *h_ptCharmBarion_prompt =new TH1D("h_ptCharmBarion_prompt","h_ptCharmBarion_prompt",NBINS,edges);

TH1D *h_yCharmMeson_prompt =new TH1D("h_yCharmMeson_prompt","h_yCharmMeson_prompt",100,-9.5,9.5);
TH1D *h_yCharmBarion_prompt =new TH1D("h_yCharmBarion_prompt","h_yCharmBarion_prompt",100,-9.5,9.5);

TH1D *h_ptCharmAntiMeson_prompt =new TH1D("h_ptCharmAntiMeson_prompt","h_ptCharmAntiMeson_prompt",NBINS,edges);
TH1D *h_ptCharmAntiBarion_prompt =new TH1D("h_ptCharmAntiBarion_prompt","h_ptCharmAntiBarion_prompt",NBINS,edges);

TH1D *h_yCharmAntiMeson_prompt =new TH1D("h_yCharmAntiMeson_prompt","h_yCharmAntiMeson_prompt",100,-9.5,9.5);
TH1D *h_yCharmAntiBarion_prompt =new TH1D("h_yCharmAntiBarion_prompt","h_yCharmAntiBarion_prompt",100,-9.5,9.5);

//Hist to store charm particle from beauty
TH1D *h_ptAllCharm_beauty =new TH1D("h_ptAllCharm_beauty","h_ptAllCharm_beauty",NBINS,edges);
TH1D *h_yAllCharm_beauty =new TH1D("h_yAllCharm_beauty","h_yAllCharm_beauty",100,-9.5,9.5);

TH1D *h_ptAllAntiCharm_beauty =new TH1D("h_ptAllAntiCharm_beauty","h_ptAllAntiCharm_beauty",NBINS,edges);
TH1D *h_yAllAntiCharm_beauty =new TH1D("h_yAllAntiCharm_beauty","h_yAllAntiCharm_beauty",100,-9.5,9.5);

TH1D *h_ptCharmMeson_beauty =new TH1D("h_ptCharmMeson_beauty","h_ptCharmMeson_beauty",NBINS,edges);
TH1D *h_ptCharmBarion_beauty =new TH1D("h_ptCharmBarion_beauty","h_ptCharmBarion_beauty",NBINS,edges);

TH1D *h_yCharmMeson_beauty =new TH1D("h_yCharmMeson_beauty","h_yCharmMeson_beauty",100,-9.5,9.5);
TH1D *h_yCharmBarion_beauty =new TH1D("h_yCharmBarion_beauty","h_yCharmBarion_beauty",100,-9.5,9.5);

TH1D *h_ptCharmAntiMeson_beauty =new TH1D("h_ptCharmAntiMeson_beauty","h_ptCharmAntiMeson_beauty",NBINS,edges);
TH1D *h_ptCharmAntiBarion_beauty =new TH1D("h_ptCharmAntiBarion_beauty","h_ptCharmAntiBarion_beauty",NBINS,edges);

TH1D *h_yCharmAntiMeson_beauty =new TH1D("h_yCharmAntiMeson_beauty","h_yCharmAntiMeson_beauty",100,-9.5,9.5);
TH1D *h_yCharmAntiBarion_beauty =new TH1D("h_yCharmAntiBarion_beauty","h_yCharmAntiBarion_beauty",100,-9.5,9.5);

//Hist to store Charm Meson

TH1D *h_ptDzero =new TH1D("h_ptDzero","h_ptDzero",NBINS,edges);
TH1D *h_ptDzero_prompt =new TH1D("h_ptDzero_prompt","h_ptDzero_prompt",NBINS,edges);

TH1D *h_ptDplus =new TH1D("h_ptDplus","h_ptDplus",NBINS,edges);
TH1D *h_ptDplus_prompt =new TH1D("h_ptDplus_prompt","h_ptDplus_prompt",NBINS,edges);

TH1D *h_ptDstrangeplus =new TH1D("h_ptDstrangeplus","h_ptDstrangeplus",NBINS,edges);
TH1D *h_ptDstrangeplus_prompt =new TH1D("h_ptDstrangeplus_prompt","h_ptDstrangeplus_prompt",NBINS,edges);

TH1D *h_ptCharmMeson_rest_prompt =new TH1D("h_ptCharmMeson_rest_prompt","h_ptCharmMeson_rest_prompt",NBINS,edges);

TH1D *h_ptLambdacharmedplus=new TH1D("h_ptLambdacharmedplus","h_ptLambdacharmedplus",NBINS,edges);
TH1D *h_ptLambdacharmedplus_prompt=new TH1D("h_ptLambdacharmedplus_prompt","h_ptLambdacharmedplus_prompt",NBINS,edges);

TH1D *h_ptSigmacharmedzero=new TH1D("h_ptSigmacharmedzero","h_ptSigmacharmedzero",NBINS,edges);
TH1D *h_ptSigmacharmedzero_prompt=new TH1D("h_ptSigmacharmedzero_prompt","h_ptSigmacharmedzero_prompt",NBINS,edges);

TH1D *h_ptSigmacharmedplus=new TH1D("h_ptSigmacharmedplus","h_ptSigmacharmedplus",NBINS,edges);
TH1D *h_ptSigmacharmedplus_prompt=new TH1D("h_ptSigmacharmedplus_prompt","h_ptSigmacharmedplus_prompt",NBINS,edges);

TH1D *h_ptCharmBarion_rest_prompt=new TH1D("h_ptCharmBarion_rest_prompt","h_ptCharmBarion_rest_prompt",NBINS,edges);

TH1D *h_ptAntiDzero =new TH1D("h_ptAntiDzero","h_ptAntiDzero",NBINS,edges);
TH1D *h_ptAntiDzero_prompt =new TH1D("h_ptAntiDzero_prompt","h_ptAntiDzero_prompt",NBINS,edges);

TH1D *h_ptAntiDplus =new TH1D("h_ptAntiDplus","h_ptAntiDplus",NBINS,edges);
TH1D *h_ptAntiDplus_prompt =new TH1D("h_ptAntiDplus_prompt","h_ptAntiDplus_prompt",NBINS,edges);

TH1D *h_ptAntiDstrangeplus =new TH1D("h_ptAntiDstrangeplus","h_ptAntiDstrangeplus",NBINS,edges);
TH1D *h_ptAntiDstrangeplus_prompt =new TH1D("h_ptAntiDstrangeplus_prompt","h_ptAntiDstrangeplus_prompt",NBINS,edges);

TH1D *h_ptCharmAntiMeson_rest_prompt =new TH1D("h_ptCharmAntiMeson_rest_prompt","h_ptCharmAntiMeson_rest_prompt",NBINS,edges);

TH1D *h_ptAntiLambdacharmedplus=new TH1D("h_ptAntiLambdacharmedplus","h_ptAntiLambdacharmedplus",NBINS,edges);
TH1D *h_ptAntiLambdacharmedplus_prompt=new TH1D("h_ptAntiLambdacharmedplus_prompt","h_ptAntiLambdacharmedplus_prompt",NBINS,edges);

TH1D *h_ptCharmAntiBarion_rest_prompt=new TH1D("h_ptCharmAntiBarion_rest_prompt","h_ptCharmAntiBarion_rest_prompt",NBINS,edges);

TH1D *h_yDzero =new TH1D("h_yDzero","h_yDzero",100,-9.5,9.5);
TH1D *h_yDzero_prompt =new TH1D("h_yDzero_prompt","h_yDzero_prompt",100,-9.5,9.5);

TH1D *h_yDplus =new TH1D("h_yDplus","h_yDplus",NBINS,edges);
TH1D *h_yDplus_prompt =new TH1D("h_yDplus_prompt","h_yDplus_prompt",100,-9.5,9.5);

TH1D *h_yDstrangeplus =new TH1D("h_yDstrangeplus","h_yDstrangeplus",100,-9.5,9.5);
TH1D *h_yDstrangeplus_prompt =new TH1D("h_yDstrangeplus_prompt","h_yDstrangeplus_prompt",100,-9.5,9.5);

TH1D *h_yCharmMeson_rest_prompt =new TH1D("h_yCharmMeson_rest_prompt","h_yCharmMeson_rest_prompt",100,-9.5,9.5);

TH1D *h_yLambdacharmedplus=new TH1D("h_yLambdacharmedplus","h_yLambdacharmedplus",100,-9.5,9.5);
TH1D *h_yLambdacharmedplus_prompt=new TH1D("h_yLambdacharmedplus_prompt","h_yLambdacharmedplus_prompt",100,-9.5,9.5);

TH1D *h_ySigmacharmedzero=new TH1D("h_ySigmacharmedzero","h_ySigmacharmedzero",100,-9.5,9.5);
TH1D *h_ySigmacharmedzero_prompt=new TH1D("h_ySigmacharmedzero_prompt","h_ySigmacharmedzero_prompt",100,-9.5,9.5);

TH1D *h_ySigmacharmedplus=new TH1D("h_ySigmacharmedplus","h_ySigmacharmedplus",100,-9.5,9.5);
TH1D *h_ySigmacharmedplus_prompt=new TH1D("h_ySigmacharmedplus_prompt","h_ySigmacharmedplus_prompt",100,-9.5,9.5);

TH1D *h_yCharmBarion_rest_prompt=new TH1D("h_yCharmBarion_rest_prompt","h_yCharmBarion_rest_prompt",100,-9.5,9.5);

TH1D *h_yAntiDzero =new TH1D("h_yAntiDzero","h_yAntiDzero",100,-9.5,9.5);
TH1D *h_yAntiDzero_prompt =new TH1D("h_yAntiDzero_prompt","h_yAntiDzero_prompt",100,-9.5,9.5);

TH1D *h_yAntiDplus =new TH1D("h_yAntiDplus","h_yAntiDplus",100,-9.5,9.5);
TH1D *h_yAntiDplus_prompt =new TH1D("h_yAntiDplus_prompt","h_yAntiDplus_prompt",100,-9.5,9.5);

TH1D *h_yAntiDstrangeplus =new TH1D("h_yAntiDstrangeplus","h_yAntiDstrangeplus",100,-9.5,9.5);
TH1D *h_yAntiDstrangeplus_prompt =new TH1D("h_yAntiDstrangeplus_prompt","h_yAntiDstrangeplus_prompt",100,-9.5,9.5);

TH1D *h_yCharmAntiMeson_rest_prompt =new TH1D("h_yCharmAntiMeson_rest_prompt","h_yCharmAntiMeson_rest_prompt",100,-9.5,9.5);

TH1D *h_yAntiLambdacharmedplus=new TH1D("h_yAntiLambdacharmedplus","h_yAntiLambdacharmedplus",100,-9.5,9.5);
TH1D *h_yAntiLambdacharmedplus_prompt=new TH1D("h_yAntiLambdacharmedplus_prompt","h_yAntiLambdacharmedplus_prompt",100,-9.5,9.5);

TH1D *h_yCharmAntiBarion_rest_prompt=new TH1D("h_yCharmAntiBarion_rest_prompt","h_yCharmAntiBarion_rest_prompt",100,-9.5,9.5);

//Hist Muon PLUS from prompt particle
static const Int_t PT_NBINS_muon = 22;
Double_t PT_edges_muon[PT_NBINS_muon + 1] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,4.0,5.0,7.5,10.0,15};
static const Int_t Y_NBINS_muon = 26;
Double_t Y_edges_muon[Y_NBINS_muon + 1] = {-9.0,-7.0,-6.0,-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,9.0};

TH1D *h_ptMuon_plus= new TH1D("h_ptMuon_plus","h_ptMuon_plus",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_Charm= new TH1D("h_ptMuon_plus_Charm","h_ptMuon_plus_Charm",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_Charm_Meson= new TH1D("h_ptMuon_plus_Charm_Meson","h_ptMuon_plus_Charm_Meson",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_Charm_Barion= new TH1D("h_ptMuon_plus_Charm_Barion","h_ptMuon_plus_Charm_Barion",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_AllBeauty= new TH1D("h_ptMuon_plus_AllBeauty","h_ptMuon_plus_AllBeauty",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_Beauty= new TH1D("h_ptMuon_plus_Beauty","h_ptMuon_plus_Beauty",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_BeautyCharm= new TH1D("h_ptMuon_plus_BeautyCharm","h_ptMuon_plus_BeautyCharm",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_Dzero= new TH1D("h_ptMuon_plus_Dzero","h_ptMuon_plus_Dzero",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_Dplus= new TH1D("h_ptMuon_plus_Dplus","h_ptMuon_plus_Dplus",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_Dstrangeplus= new TH1D("h_ptMuon_plus_Dstrangeplus","h_ptMuon_plus_Dstrangeplus",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_Lambda= new TH1D("h_ptMuon_plus_Lambda","h_ptMuon_plus_Lambda",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_plus_Charmrest= new TH1D("h_ptMuon_plus_Charmrest","h_ptMuon_plus_Charmrest",PT_NBINS_muon,PT_edges_muon);

TH1D *h_yMuon_plus= new TH1D("h_yMuon_plus","h_yMuon_plus", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_Charm= new TH1D("h_yMuon_plus_Charm","h_yMuon_plus_Charm", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_Charm_Meson= new TH1D("h_yMuon_plus_Charm_Meson","h_yMuon_plus_Charm_Meson",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_Charm_Barion= new TH1D("h_yMuon_plus_Charm_Barion","h_yMuon_plus_Charm_Barion",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_AllBeauty= new TH1D("h_yMuon_plus_AllBeauty","h_yMuon_plus_AllBeauty", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_BeautyCharm= new TH1D("h_yMuon_plus_BeautyCharm","h_yMuon_plus_BeautyCharm", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_Beauty= new TH1D("h_yMuon_plus_Beauty","h_yMuon_plus_Beauty", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_Dzero= new TH1D("h_yMuon_plus_Dzero","h_yMuon_plus_Dzero",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_Dplus= new TH1D("h_yMuon_plus_Dplus","h_yMuon_plus_Dplus",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_Dstrangeplus= new TH1D("h_yMuon_plus_Dstrangeplus","h_yMuon_plus_Dstrangeplus",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_Lambda= new TH1D("h_yMuon_plus_Lambda","h_yMuon_plus_Lambda",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_plus_Charmrest= new TH1D("h_yMuon_plus_Charmrest","h_yMuon_plus_Charmrest",Y_NBINS_muon, Y_edges_muon);

TH1D *h_ptMuon_minus= new TH1D("h_ptMuon_minus","h_ptMuon_minus",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_Charm= new TH1D("h_ptMuon_minus_Charm","h_ptMuon_minus_Charm",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_Charm_Meson= new TH1D("h_ptMuon_minus_Charm_Meson","h_ptMuon_minus_Charm_Meson",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_Charm_Barion= new TH1D("h_ptMuon_minus_Charm_Barion","h_ptMuon_minus_Charm_Barion",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_AllBeauty= new TH1D("h_ptMuon_minus_AllBeauty","h_ptMuon_minus_AllBeauty",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_BeautyCharm= new TH1D("h_ptMuon_minus_BeautyCharm","h_ptMuon_minus_BeautyCharm",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_Beauty= new TH1D("h_ptMuon_minus_Beauty","h_ptMuon_minus_Beauty",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_Dzero= new TH1D("h_ptMuon_minus_Dzero","h_ptMuon_minus_Dzero",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_Dplus= new TH1D("h_ptMuon_minus_Dplus","h_ptMuon_minus_Dplus",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_Dstrangeplus= new TH1D("h_ptMuon_minus_Dstrangeplus","h_ptMuon_minus_Dstrangeplus",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_Lambda= new TH1D("h_ptMuon_minus_Lambda","h_ptMuon_minus_Lambda",PT_NBINS_muon,PT_edges_muon);

TH1D *h_ptMuon_minus_Charmrest= new TH1D("h_ptMuon_minus_Charmrest","h_ptMuon_minus_Charmrest",PT_NBINS_muon,PT_edges_muon);

TH1D *h_yMuon_minus= new TH1D("h_yMuon_minus","h_yMuon_minus", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_Charm= new TH1D("h_yMuon_minus_Charm","h_yMuon_minus_Charm", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_Charm_Meson= new TH1D("h_yMuon_minus_Charm_Meson","h_yMuon_minus_Charm_Meson",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_Charm_Barion= new TH1D("h_yMuon_minus_Charm_Barion","h_yMuon_minus_Charm_Barion",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_AllBeauty= new TH1D("h_yMuon_minus_AllBeauty","h_yMuon_minus_AllBeauty", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_Beauty= new TH1D("h_yMuon_minus_Beauty","h_yMuon_minus_Beauty", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_BeautyCharm= new TH1D("h_yMuon_minus_BeautyCharm","h_yMuon_minus_BeautyCharm", Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_Dzero= new TH1D("h_yMuon_minus_Dzero","h_yMuon_minus_Dzero",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_Dplus= new TH1D("h_yMuon_minus_Dplus","h_yMuon_minus_Dplus",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_Dstrangeplus= new TH1D("h_yMuon_minus_Dstrangeplus","h_yMuon_minus_Dstrangeplus",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_Lambda= new TH1D("h_yMuon_minus_Lambda","h_yMuon_minus_Lambda",Y_NBINS_muon, Y_edges_muon);

TH1D *h_yMuon_minus_Charmrest= new TH1D("h_yMuon_minus_rest","h_yMuon_minus_rest",Y_NBINS_muon, Y_edges_muon);

TH2D *h_yptMuon_plus_Charm=new TH2D("h_yptMuon_plus_Charm","2D Mu+ pt y",100,-9,9,20,0,20);

TH2D *h_yptMuon_minus_Charm=new TH2D("h_yptMuon_minus_Charm","2D Mu- pt y",100,-9,9,20,0,20);

//Hist to store muon multiplicity

TH1D *h_mult_mu_plus_charm = new TH1D("h_mult_mu_plus_charm","mu plus multiplicity from charm ", 10,-0.5, 9.5);

TH1D *h_mult_mu_minus_charm = new TH1D("h_mult_mu_minus_charm","mu minus multiplicity from charm", 10,-0.5, 9.5);

TH1D *h_mult_mu_OSpair_charm = new TH1D("h_mult_mu_OSpair_charm","mu OS pair from charm ", 10,-0.5, 9.5);

TH1D *h_mult_mu_plus_beauty = new TH1D("h_mult_mu_plus_beauty","mu plus multiplicity from beauty ", 10,-0.5, 9.5);

TH1D *h_mult_mu_minus_beauty = new TH1D("h_mult_mu_minus_beauty","mu minus multiplicity from beauty", 10,-0.5, 9.5);

TH1D *h_mult_mu_OSpair_beauty = new TH1D("h_mult_mu_OSpair_beauty","mu OS pair from beauty ", 10,-0.5, 9.5);

//hist to store DiMuon

static const Int_t PT_NBINS_Dimuon = 19;
Double_t PT_edges_Dimuon[PT_NBINS_Dimuon + 1]={0.0,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,9.0,12.5,17.5,22.5,27.5};

static const Int_t Mass_NBINS_Dimuon = 11;
Double_t Mass_edges_Dimuon[PT_NBINS_Dimuon + 1]={4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,9.0,12.5,17.5,22.5};

static const Int_t Mass_NBINS_Dimuon_Light = 8;
Double_t Mass_edges_Dimuon_Light[Mass_NBINS_Dimuon_Light + 1]={4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0};

TH1D *h_ptDiMu= new TH1D("h_ptDiMu","dN/dpt DiMuon single mu cut",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu= new TH1D("h_mDiMu","dN/dM DiMuon single mu cut",Mass_NBINS_Dimuon,Mass_edges_Dimuon);
TH1D *h_yDiMu= new TH1D("h_yDiMu","dN/dy DiMuon single mu cut", 100, -9.5, 9.5);
TH1D *h_phiDiMu= new TH1D("h_phiDiMu","dN/d#phi DiMuon single mu cut", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_Charm= new TH1D("h_ptDiMu_Charm","dN/dpt DiMuon from Charm single mu cut",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_Charm= new TH1D("h_mDiMu_Charm","dN/dM DiMuon from Charm single mu cut",44,-0.5,22.5);
TH1D *h_yDiMu_Charm= new TH1D("h_yDiMu_Charm","dN/dy DiMuon from Charm single mu cut", 100, -9.5, 9.5);
TH1D *h_phiDiMu_Charm= new TH1D("h_phiDiMu_Charm","dN/d#phi DiMuon from Charm single mu cut", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_AllBeauty= new TH1D("h_ptDiMu_AllBeauty","DiMuon from Beauty single mu cut",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_AllBeauty= new TH1D("h_mDiMu_AllBeauty","DiMuon from Beauty single mu cut",44,-0.5,22.5);
TH1D *h_yDiMu_AllBeauty= new TH1D("h_yDiMu_AllBeauty","DiMuon from Beauty single mu cut", 100, -9.5, 9.5);
TH1D *h_phiDiMu_AllBeauty= new TH1D("h_phiDiMu_AllBeauty","DiMuon from Beauty single mu cut", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_Beauty= new TH1D("h_ptDiMu_Beauty","DiMuon from direct Beauty single mu cut",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_Beauty= new TH1D("h_mDiMu_Beauty","DiMuon from direct Beauty single mu cut",44,-0.5,22.5);
TH1D *h_yDiMu_Beauty= new TH1D("h_yDiMu_Beauty","DiMuon from direct Beauty single mu cut", 100, -9.5, 9.5);
TH1D *h_phiDiMu_Beauty= new TH1D("h_phiDiMu_Beauty","DiMuon from direct Beauty single mu cut", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_BeautyCharm= new TH1D("h_ptDiMu_BeautyCharm","DiMuon from Beauty via Charm single mu cut",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_BeautyCharm= new TH1D("h_mDiMu_BeautyCharm","DiMuon from Beauty via Charm single mu cut",44,-0.5,22.5);
TH1D *h_yDiMu_BeautyCharm= new TH1D("h_yDiMu_BeautyCharm","DiMuon from Beauty via Charm single mu cut", 100, -9.5, 9.5);
TH1D *h_phiDiMu_BeautyCharm= new TH1D("h_phiDiMu_BeautyCharm","DiMuon from Beauty viaCharm single mu cut", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu4MassCut= new TH1D("h_ptDiMu4MassCut","DiMuon M>4 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu4MassCut= new TH1D("h_mDiMu4MassCut","DiMuon M>4 GeV/c^{2}",Mass_NBINS_Dimuon,Mass_edges_Dimuon);
TH1D *h_yDiMu4MassCut= new TH1D("h_yDiMu4MassCut","DiMuon M>4 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu4MassCut= new TH1D("h_phiDiMu4MassCut","DiMuon M>4 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_Charm4MassCut= new TH1D("h_ptDiMu_Charm4MassCut","DiMuon from Charm M>4 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_Charm4MassCut= new TH1D("h_mDiMu_Charm4MassCut","DiMuon from Charm M>4 GeV/c^{2}",Mass_NBINS_Dimuon,Mass_edges_Dimuon);
TH1D *h_yDiMu_Charm4MassCut= new TH1D("h_yDiMu_Charm4MassCut","DiMuon from Charm M>4 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu_Charm4MassCut= new TH1D("h_phiDiMu_Charm4MassCut","DiMuon from Charm M>4 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_AllBeauty4MassCut= new TH1D("h_ptDiMu_AllBeauty4MassCut","DiMuon from Beauty M>4 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_AllBeauty4MassCut= new TH1D("h_mDiMu_AllBeauty4MassCut","DiMuon from Beauty M>4 GeV/c^{2}",Mass_NBINS_Dimuon,Mass_edges_Dimuon);
TH1D *h_yDiMu_AllBeauty4MassCut= new TH1D("h_yDiMu_AllBeauty4MassCut","DiMuon from Beauty M>4 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu_AllBeauty4MassCut= new TH1D("h_phiDiMu_AllBeauty4MassCut","DiMuon from Beauty M>4 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_Beauty4MassCut= new TH1D("h_ptDiMu_Beauty4MassCut","DiMuon from direct Beauty M>4 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_Beauty4MassCut= new TH1D("h_mDiMu_Beauty4MassCut","DiMuon from direct Beauty M>4 GeV/c^{2}",Mass_NBINS_Dimuon,Mass_edges_Dimuon);
TH1D *h_yDiMu_Beauty4MassCut= new TH1D("h_yDiMu_Beauty4MassCut","DiMuon from direct Beauty M>4 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu_Beauty4MassCut= new TH1D("h_phiDiMu_Beauty4MassCut","DiMuon from direct Beauty M>4 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());


TH1D *h_ptDiMu_BeautyCharm4MassCut= new TH1D("h_ptDiMu_BeautyCharm4MassCut","DiMuon from Beauty via Charm M>4 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_BeautyCharm4MassCut= new TH1D("h_mDiMu_BeautyCharm4MassCut","DiMuon from Beauty via Charm M>4 GeV/c^{2}",Mass_NBINS_Dimuon,Mass_edges_Dimuon);
TH1D *h_yDiMu_BeautyCharm4MassCut= new TH1D("h_yDiMu_BeautyCharm4MassCut","DiMuon from Beauty via Charm M>4 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu_BeautyCharm4MassCut= new TH1D("h_phiDiMu_BeautyCharm4MassCut","DiMuon from Beauty via Charm M>4 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_light= new TH1D("h_ptDiMu_light","DiMuon 4<M<8 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_light= new TH1D("h_mDiMu_light","DiMuon 4<M<8 GeV/c^{2}",Mass_NBINS_Dimuon_Light,Mass_edges_Dimuon_Light);
TH1D *h_yDiMu_light= new TH1D("h_yDiMu_light","DiMuon 4<M<8 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu_light= new TH1D("h_phiDiMu_light","DiMuon 4<M<8 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_Charm_light= new TH1D("h_ptDiMu_Charm_light","DiMuon from Charm 4<M<8 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_Charm_light= new TH1D("h_mDiMu_Charm_light","DiMuon from Charm 4<M<8 GeV/c^{2}",Mass_NBINS_Dimuon_Light,Mass_edges_Dimuon_Light);
TH1D *h_yDiMu_Charm_light= new TH1D("h_yDiMu_Charm_light","DiMuon from Charm 4<M<8 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu_Charm_light= new TH1D("h_phiDiMu_Charm_light","DiMuon from Charm 4<M<8 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_AllBeauty_light= new TH1D("h_ptDiMu_AllBeauty_light","DiMuon from Beauty 4<M<8 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_AllBeauty_light= new TH1D("h_mDiMu_AllBeauty_light","DiMuon from Beauty 4<M<8 GeV/c^{2}",Mass_NBINS_Dimuon_Light,Mass_edges_Dimuon_Light);
TH1D *h_yDiMu_AllBeauty_light= new TH1D("h_yDiMu_AllBeauty_light","DiMuon from Beauty 4<M<8 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu_AllBeauty_light= new TH1D("h_phiDiMu_AllBeauty_light","DiMuon from Beauty 4<M<8 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_Beauty_light= new TH1D("h_ptDiMu_Beauty_light","DiMuon from Beauty 4<M<8 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_Beauty_light= new TH1D("h_mDiMu_Beauty_light","DiMuon from Beauty 4<M<8 GeV/c^{2}",Mass_NBINS_Dimuon_Light,Mass_edges_Dimuon_Light);
TH1D *h_yDiMu_Beauty_light= new TH1D("h_yDiMu_Beauty_light","DiMuon from Beauty 4<M<8 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu_Beauty_light= new TH1D("h_phiDiMu_Beauty_light","DiMuon from Beauty 4<M<8 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());

TH1D *h_ptDiMu_BeautyCharm_light= new TH1D("h_ptDiMu_BeautyCharm_light","DiMuon from Beauty via Charm 4<M<8 GeV/c^{2}",PT_NBINS_Dimuon,PT_edges_Dimuon);
TH1D *h_mDiMu_BeautyCharm_light= new TH1D("h_mDiMu_BeautyCharm_light","DiMuon from Beauty via Charm 4<M<8 GeV/c^{2}",Mass_NBINS_Dimuon_Light,Mass_edges_Dimuon_Light);
TH1D *h_yDiMu_BeautyCharm_light= new TH1D("h_yDiMu_BeautyCharm_light","DiMuon from Beauty via Charm 4<M<8 GeV/c^{2}", 100, -9.5, 9.5);
TH1D *h_phiDiMu_BeautyCharm_light= new TH1D("h_phiDiMu_BeautyCharm_light","DiMuon from Beauty via Charm 4<M<8 GeV/c^{2}", 100,-TMath::Pi(),TMath::Pi());

using namespace std;

//void nuovovsvecchio();
void Particle_Analisys(int argc, char* argv[]);
void Muons_Analisys(double LowPt,double HighPt,double LowY, double HighY,int argc, char* argv[],TString input_path, TString input_dir, TString saving_path,bool Pythia_standalone);
void Plus_Muon_origin(int argc, char* argv[]);
void Minus_Muon_origin(int argc, char* argv[]);
void Particle_Analisys(int argc, char* argv[],TString input_path,TString input_dir,TString saving_path,bool Pythia_standalone);
void AntiParticle_origin(int argc, char* argv[]);
void DiMuon(int argc, char* argv[],double LowPt,double HighPt,double LowY,double HighY);
void DefvsForce(int argc, char* argv[]);
void Muon_DefvsForce_Charm(int argc, char* argv[]);
void Muon_DefvsForce_BeautyCharm(int argc, char* argv[]);
void Muon_DefvsForce_Beauty(int argc, char* argv[]);
void Process(int argc, char* argv[]);
void ALICE(int codeParticle,int argc, char* argv[]);
void LHCb(int codeParticle,int argc, char* argv[]);
void Writing_Muons_Analisys(TFile *outFile,TString name_folder,TString name_subfolder);
void Muons();
TChain *Getting_File(int argc, char* argv[],TString &fileout, TString input_path,TString input_dir,TString saving_path,bool Pythia_standalone);

void input_for_analysis(int argc, char* argv[],int &func, int &nevents,int &choice_Seed, TString &Process, int &choice_Process, TString &Change, int &choice_Change, TString &BR, int &choice_BR, bool Pythia_standalone) {

    const char *selectedprocess[5]={"SoftQCD","ccbar","bbbar","SoftQCDall","Drell-Yan"};
    const char *selectedchange_PythiaSim[4]={"Def","Config","Config+Atlas","Mode2"};
    const char *selectedBR[2]={"DefaultBR","ForceBR"};
    
    const char *selectedchange_LxplusSim[4]={"Pythia_Only","Custom_ND","Custom_ND_All",""};
    
    for(int ig=0;ig<argc;++ig){
        TString str=argv[ig];
        printf("str %d %s\n",ig,argv[ig]);
        if(str.Contains("--Function")) {
            sscanf(argv[ig+1],"%d",&func);
            printf("func %d",func);
        }
        if(str.Contains("--Events")){
            sscanf(argv[ig+1],"%d",&nevents);
            printf("nevents after setting: %d\n",nevents);
        }
        if(str.Contains("--Seed")){
            sscanf(argv[ig+1],"%d",&choice_Seed);
            printf("Seed: %d\n",choice_Seed);
        }
        if(str.Contains("--Process")){
            sscanf(argv[ig+1],"%d",&choice_Process);
            Process.Form("%s",selectedprocess[choice_Process]);
            printf("Process: %d) => %s \n", choice_Process, Process.Data());
            if(choice_Process>4 && Pythia_standalone) printf("Uncorrect Process choice %d",choice_Change);
        }
        if(str.Contains("--BR")){
            sscanf(argv[ig+1],"%d",&choice_BR);
            BR.Form("%s",selectedBR[choice_BR]);
            printf("BR: %d) => %s \n", choice_BR, BR.Data());
            if(choice_BR>1 && Pythia_standalone) printf("Uncorrect BR choice %d",choice_Change);
        }
        if(str.Contains("--Change")){
            sscanf(argv[ig+1],"%d",&choice_Change);
            if (Pythia_standalone) Change.Form("%s",selectedchange_PythiaSim[choice_Change]);
            else Change.Form("%s",selectedchange_LxplusSim[choice_Change]);
            printf("Change: %d) => %s \n", choice_BR, BR.Data());
            if(choice_Change>4) printf("Uncorrect Change choice %d",choice_Change);
        }
        
    }
}

void SetLabelsXYandSave(TH1D *h_pt,TH1D *h_y){
    h_pt->GetXaxis()->SetTitle("#it{p}_{t}");
    h_pt->GetYaxis()->SetTitle("dN/d#it{p}_{t}");
    h_pt->Write(0,2,0);
    
    h_y->GetXaxis()->SetTitle("y");
    h_y->GetYaxis()->SetTitle("dN/dy");
    h_y->Write(0,2,0);
}
void SetLabelsCSandSave(TH1D *h_cs){
    h_cs->GetXaxis()->SetTitle("#it{p}_{t}");
    h_cs->GetYaxis()->SetTitle("d^{2}#sigma/dyd#it{p}_{T} (#mu b/GeV/c)");
    h_cs->Write(0,2,0);
}

#endif
