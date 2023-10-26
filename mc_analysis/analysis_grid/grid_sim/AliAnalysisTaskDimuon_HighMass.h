#ifndef AliAnalysisTaskDimuon_HighMass_H
#define AliAnalysisTaskDimuon_HighMass_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TObjArray;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;
class AliMuonTrackCuts;

class AliAnalysisTaskDimuon_HighMass : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskDimuon_HighMass();
    AliAnalysisTaskDimuon_HighMass(const char *name);
    virtual ~AliAnalysisTaskDimuon_HighMass();

    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *);
    void SetSaving(TString);
    virtual void NotifyRun();

    void SetBeamEnergy(Double_t en) { fBeamEnergy = en; }
    void SetAnalysisType(const char *type) { fkAnalysisType = type; }
    void SetPeriod(TString period) { fPeriod = period; }

    // typedef enum
    // {
    //     kSaveHF,
    //     kSaveDY,
    //     kSaveMuonOnly,
    // } Saving_opt_t;

    virtual void SetSaving_opt(TString opt)
    {
        fSaving_opt = opt;
    }

    virtual void SetAOD_origin(TString opt) { fAOD_origin = opt; }

    void Set_Verbosity(Bool_t verbose) { fVerbose = verbose; }

private:
    AliAnalysisTaskDimuon_HighMass(const AliAnalysisTaskDimuon_HighMass &);
    AliAnalysisTaskDimuon_HighMass &operator=(const AliAnalysisTaskDimuon_HighMass &);

    TString fSaving_opt;

    TString fAOD_origin;

    Bool_t fVerbose;

    // protected:

    Double_t fBeamEnergy;       // Energy of the beam (required for the CS angle)
    const char *fkAnalysisType; // ESD or AOD based analysis
    TString fPeriod;            // period
    AliAODEvent *fAODEvent;     //! AOD event  //tolgo !
    TTree *fOutputTree;         //! tree output
    AliMuonTrackCuts *fMuonTrackCuts;

    Int_t fNMuons_gen;     // gen muon in the event
    Int_t fNHadrons_gen;   // gen muon in the event
    Int_t fNDimu_gen;      // gen dimuons in the event
    Int_t fN_HFquarks_gen; // gen c/cbar or b/bar HFquarks in the event

    Int_t fN_gamma; // gen gamma* in the event

    Int_t fNMuons_rec; // rec muon tracks in the event
    Int_t fNDimu_rec;  // rec dimuons in the event

    Int_t fN_HFquarks_rec; // gen c/cbar or b/bar HFquarks in the event

    Double_t fPercentV0M;

    static const Int_t fMuons_dim = 5000;
    static const Int_t fDimu_dim = 5000;

    Double_t fPt_gamma[fMuons_dim]; // Gamma star pt
    Double_t fM_gamma[fMuons_dim];  // Gamma star M
    Double_t fY_gamma[fMuons_dim];  // Gamma star Y

    Int_t fPDGmum_rec[fMuons_dim];           // single rec mu PDG mum
    Double_t fPt_rec[fMuons_dim];            // single rec mu pT
    Double_t fE_rec[fMuons_dim];             // single rec mu E
    Double_t fPx_rec[fMuons_dim];            // single rec mu px
    Double_t fPy_rec[fMuons_dim];            // single rec mu py
    Double_t fPz_rec[fMuons_dim];            // single rec mu pz
    Double_t fY_rec[fMuons_dim];             // single rec mu y
    Double_t fEta_rec[fMuons_dim];           // single rec mu eta
    Int_t fMatchTrig_rec[fMuons_dim];        // single rec mu match trigger
    Double_t fTrackChi2_rec[fMuons_dim];     // single rec mu chi2 track
    Double_t fMatchTrigChi2_rec[fMuons_dim]; // single rec mu chi2 of match trigger
    Int_t fCharge_rec[fMuons_dim];           // single rec mu charge
    Double_t fRAtAbsEnd_rec[fMuons_dim];     // single rec mu distance from beam center at end abs
    Int_t fpDCA_rec[fMuons_dim];             // single rec mu charge
    Double_t fPhi_rec[fMuons_dim];           // single rec mu phi
    Double_t fTheta_rec[fMuons_dim];         // single rec mu theta
    Int_t fFrom_Powheg_rec[fMuons_dim];      // check muon rec origin
    Int_t fInitial_Parton_rec[fMuons_dim];   // check muon rec original parton

    Int_t fPDG_HFquark_gen[fMuons_dim];   // single gen c/cbar PDG mum
    Double_t fPx_HFquark_gen[fMuons_dim]; // single gen c/cbar or b/bbar HFquark pT
    Double_t fPy_HFquark_gen[fMuons_dim]; // single gen c/cbar or b/bbar HFquark y
    Double_t fPz_HFquark_gen[fMuons_dim]; // single gen c/cbar or b/bbar HFquark pT
    Double_t fPt_HFquark_gen[fMuons_dim]; // single gen c/cbar or b/bbar HFquark pT
    Double_t fY_HFquark_gen[fMuons_dim];  // single gen c/cbar or b/bbar HFquark y
    Int_t fMother_index[fMuons_dim];      // single gen c/cbar or b/bbar HFquark y

    Int_t fPDGmum_gen[fMuons_dim];   // single gen mu PDG mum
    Double_t fPt_gen[fMuons_dim];    // single gen mu pT
    Double_t fE_gen[fMuons_dim];     // single gen mu E
    Double_t fPx_gen[fMuons_dim];    // single gen mu px
    Double_t fPy_gen[fMuons_dim];    // single gen mu py
    Double_t fPz_gen[fMuons_dim];    // single gen mu pz
    Double_t fY_gen[fMuons_dim];     // single gen mu y
    Double_t fEta_gen[fMuons_dim];   // single gen mu eta
    Int_t fCharge_gen[fMuons_dim];   // single gen mu charge
    Double_t fPhi_gen[fMuons_dim];   // single gen mu phi
    Double_t fTheta_gen[fMuons_dim]; // single gen mu theta
    Int_t fFrom_Powheg_gen[fMuons_dim]; // check muon gen origin
    Int_t fInitial_Parton_gen[fMuons_dim]; // check muon gen original parton

    Int_t fPDGmum_Hadron_gen[fMuons_dim]; // gen Hadron PDG mum
    Int_t fPDG_Hadron_gen[fMuons_dim];    // gen Hadron PDG
    Double_t fPt_Hadron_gen[fMuons_dim];  // gen Hadron pT
    Double_t fE_Hadron_gen[fMuons_dim];   // gen Hadron E
    Double_t fPx_Hadron_gen[fMuons_dim];  // gen Hadron px
    Double_t fPy_Hadron_gen[fMuons_dim];  // gen Hadron py
    Double_t fPz_Hadron_gen[fMuons_dim];  // gen Hadron pz
    Double_t fY_Hadron_gen[fMuons_dim];   // gen Hadron y
    Double_t fEta_Hadron_gen[fMuons_dim]; // gen Hadron eta
    Int_t fHadronFrom_Powheg_gen[fMuons_dim]; // check muon gen origin

    Int_t fDimuMu_gen[fDimu_dim][2];   // reference to single gen mus
    Double_t fDimuPt_gen[fDimu_dim];   // gen dimuon pT
    Double_t fDimuPx_gen[fDimu_dim];   // gen dimuon px
    Double_t fDimuPy_gen[fDimu_dim];   // gen dimuon py
    Double_t fDimuPz_gen[fDimu_dim];   // gen dimuon pz
    Double_t fDimuY_gen[fDimu_dim];    // gen dimuon y
    Double_t fDimuMass_gen[fDimu_dim]; // gen dimuon invariant mass
    Int_t fDimuCharge_gen[fDimu_dim];  // gen dimuon charge

    Int_t fDimuMu_rec[fDimu_dim][2];    // reference to single rec mus
    Double_t fDimuPt_rec[fDimu_dim];    // rec dimuon pT
    Double_t fDimuPx_rec[fDimu_dim];    // rec dimuon px
    Double_t fDimuPy_rec[fDimu_dim];    // rec dimuon py
    Double_t fDimuPz_rec[fDimu_dim];    // rec dimuon pz
    Double_t fDimuY_rec[fDimu_dim];     // rec dimuon y
    Double_t fDimuMass_rec[fDimu_dim];  // rec dimuon invariant mass
    Int_t fDimuCharge_rec[fDimu_dim];   // rec dimuon charge
    Int_t fDimuMatch_rec[fDimu_dim];    // rec dimuon match
    Double_t fDimuPhi_rec[fDimu_dim];   // rec dimuon phi
    Double_t fDimuTheta_rec[fDimu_dim]; // rec dimuon theta

    ClassDef(AliAnalysisTaskDimuon_HighMass, 1);
};

#endif
