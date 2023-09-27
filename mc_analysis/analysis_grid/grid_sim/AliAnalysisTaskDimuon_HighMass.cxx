/**************************************************************************
 * Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Coiibutors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAnalysisTaskDimuon_HighMass.cxx $ */

// ROOT includes
#include "TROOT.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGrid.h"

#include "AliInputEventHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODMCParticle.h"
#include "AliMultSelection.h"
#include "AliMuonTrackCuts.h"

#include "AliAnalysisTaskDimuon_HighMass.h"

Int_t IsPrompt(AliAODMCParticle *mcp_mumum, TClonesArray *mcarray);
Int_t *Muon_ancestor_features(AliAODMCParticle *mcp_mumum, TClonesArray *mcarray, TString AOD_origin);

ClassImp(AliAnalysisTaskDimuon_HighMass)
    //__________________________________________________________________________
    AliAnalysisTaskDimuon_HighMass::AliAnalysisTaskDimuon_HighMass() : AliAnalysisTaskSE(),
                                                                       fBeamEnergy(0.),
                                                                       fkAnalysisType(0x0),
                                                                       fPeriod(0x0),
                                                                       fAODEvent(0x0),
                                                                       fOutputTree(0x0),
                                                                       fMuonTrackCuts(0x0),
                                                                       fNMuons_gen(0x0),
                                                                       fNHadrons_gen(0x0),
                                                                       fN_gamma(0x0),
                                                                       fNDimu_gen(0x0),
                                                                       fNMuons_rec(0x0),
                                                                       fPercentV0M(0x0),
                                                                       fVerbose(kFALSE)

{
    /// Default ctor.
    printf("fMuons_dim %d\n", fMuons_dim);
    for (Int_t i = 0; i < fMuons_dim; i++)
    {

        fPDG_HFquark_gen[i] = 9999.;
        fPx_HFquark_gen[i] = 9999.;
        fPy_HFquark_gen[i] = 9999.;
        fPz_HFquark_gen[i] = 9999.;
        fPt_HFquark_gen[i] = 9999.;
        fY_HFquark_gen[i] = 9999.;
        fMother_index[i] = 9999;

        fPt_gamma[i] = 9999.;
        fM_gamma[i] = 9999.;
        fY_gamma[i] = 9999.;

        fPDGmum_gen[i] = 9999.;
        fPt_gen[i] = 999.;
        fE_gen[i] = 999.;
        fPx_gen[i] = 999;
        fPy_gen[i] = 999;
        fPz_gen[i] = 999;
        fY_gen[i] = 999.;
        fEta_gen[i] = 999.;
        fPhi_gen[i] = 999.;
        fTheta_gen[i] = 999.;
        fFrom_Powheg_gen[i] = 999;

        fPDGmum_Hadron_gen[i] = 999; // gen Hadron PDG mum
        fPDG_Hadron_gen[i] = 999;    // gen Hadron PDG
        fPt_Hadron_gen[i] = 999;     // gen Hadron pT
        fE_Hadron_gen[i] = 999;      // gen Hadron E
        fPx_Hadron_gen[i] = 999;     // gen Hadron px
        fPy_Hadron_gen[i] = 999;     // gen Hadron py
        fPz_Hadron_gen[i] = 999;     // gen Hadron pz
        fY_Hadron_gen[i] = 999;      // gen Hadron y
        fEta_Hadron_gen[i] = 999;    // gen Hadron eta

        fPDGmum_rec[i] = 9999.;
        fPt_rec[i] = 999.;
        fE_rec[i] = 999.;
        fPx_rec[i] = 999;
        fPy_rec[i] = 999;
        fPz_rec[i] = 999;
        fY_rec[i] = 999.;
        fEta_rec[i] = 999.;
        fMatchTrig_rec[i] = 999.;
        fTrackChi2_rec[i] = 999.;
        fMatchTrigChi2_rec[i] = 999.;
        fCharge_rec[i] = 999;
        fRAtAbsEnd_rec[i] = 999;
        fpDCA_rec[i] = 999.;
        fPhi_rec[i] = 999.;
        fTheta_rec[i] = 999.;
        fFrom_Powheg_rec[i] = 999;
    }
    for (Int_t i = 0; i < fDimu_dim; i++)
    {
        fDimuPt_gen[i] = 999.;
        fDimuPx_gen[i] = 999.;
        fDimuPy_gen[i] = 999.;
        fDimuPz_gen[i] = 999.;
        fDimuY_gen[i] = 999.;
        fDimuMass_gen[i] = 999.;
        fDimuCharge_gen[i] = 999;

        fDimuPt_rec[i] = 999.;
        fDimuPx_rec[i] = 999.;
        fDimuPy_rec[i] = 999.;
        fDimuPz_rec[i] = 999.;
        fDimuY_rec[i] = 999.;
        fDimuMass_rec[i] = 999.;
        fDimuCharge_rec[i] = 999;
        fDimuMatch_rec[i] = 999;
        for (Int_t k = 0; k < 2; k++)
            fDimuMu_rec[i][k] = 999;
        fDimuPhi_rec[i] = 999.;
        fDimuTheta_rec[i] = 999.;
    }

    fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
    fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
    fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
}

//__________________________________________________________________________
AliAnalysisTaskDimuon_HighMass::AliAnalysisTaskDimuon_HighMass(const char *name) : AliAnalysisTaskSE(name),
                                                                                   fBeamEnergy(0.),
                                                                                   fkAnalysisType(0x0),
                                                                                   fPeriod(0x0),
                                                                                   fAODEvent(0x0),
                                                                                   fOutputTree(0x0),
                                                                                   fMuonTrackCuts(0x0),
                                                                                   fNMuons_gen(0x0),
                                                                                   fNHadrons_gen(0x0),
                                                                                   fN_gamma(0x0),
                                                                                   fNDimu_gen(0x0),
                                                                                   fNMuons_rec(0x0),
                                                                                   fPercentV0M(0x0),
                                                                                   fVerbose(kFALSE)
{
    //
    // Constructor. Initialization of Inputs and Outputs
    //

    Info("AliAnalysisTaskDimuon_HighMass", "Calling Constructor");

    printf("fMuons_dim %d\n", fMuons_dim);
    for (Int_t i = 0; i < fMuons_dim; i++)
    {
        fPDG_HFquark_gen[i] = 9999.;
        fPx_HFquark_gen[i] = 9999.;
        fPy_HFquark_gen[i] = 9999.;
        fPz_HFquark_gen[i] = 9999.;
        fPt_HFquark_gen[i] = 9999.;
        fY_HFquark_gen[i] = 9999.;
        fMother_index[i] = 9999;

        fPt_gamma[i] = 9999.;
        fM_gamma[i] = 9999.;
        fY_gamma[i] = 9999.;

        fPDGmum_gen[i] = 9999.;
        fPt_gen[i] = 999.;
        fE_gen[i] = 999.;
        fPx_gen[i] = 999;
        fPy_gen[i] = 999;
        fPz_gen[i] = 999;
        fY_gen[i] = 999.;
        fEta_gen[i] = 999.;
        fPhi_gen[i] = 999.;
        fTheta_gen[i] = 999.;
        fFrom_Powheg_gen[i] = 999;

        fPDGmum_Hadron_gen[i] = 999; // gen Hadron PDG mum
        fPDG_Hadron_gen[i] = 999;    // gen Hadron PDG
        fPt_Hadron_gen[i] = 999;     // gen Hadron pT
        fE_Hadron_gen[i] = 999;      // gen Hadron E
        fPx_Hadron_gen[i] = 999;     // gen Hadron px
        fPy_Hadron_gen[i] = 999;     // gen Hadron py
        fPz_Hadron_gen[i] = 999;     // gen Hadron pz
        fY_Hadron_gen[i] = 999;      // gen Hadron y
        fEta_Hadron_gen[i] = 999;    // gen Hadron eta

        fPDGmum_rec[i] = 9999.;
        fPt_rec[i] = 999.;
        fE_rec[i] = 999.;
        fPx_rec[i] = 999;
        fPy_rec[i] = 999;
        fPz_rec[i] = 999;
        fY_rec[i] = 999.;
        fEta_rec[i] = 999.;
        fMatchTrig_rec[i] = 999.;
        fTrackChi2_rec[i] = 999.;
        fMatchTrigChi2_rec[i] = 999.;
        fCharge_rec[i] = 999;
        fRAtAbsEnd_rec[i] = 999;
        fpDCA_rec[i] = 999.;
        fPhi_rec[i] = 999.;
        fTheta_rec[i] = 999.;
        fFrom_Powheg_rec[i] = 999;
    }
    for (Int_t i = 0; i < fDimu_dim; i++)
    {
        fDimuPt_gen[i] = 999.;
        fDimuPx_gen[i] = 999.;
        fDimuPy_gen[i] = 999.;
        fDimuPz_gen[i] = 999.;
        fDimuY_gen[i] = 999.;
        fDimuMass_gen[i] = 999.;
        fDimuCharge_gen[i] = 999;

        fDimuPt_rec[i] = 999.;
        fDimuPx_rec[i] = 999.;
        fDimuPy_rec[i] = 999.;
        fDimuPz_rec[i] = 999.;
        fDimuY_rec[i] = 999.;
        fDimuMass_rec[i] = 999.;
        fDimuCharge_rec[i] = 999;
        fDimuMatch_rec[i] = 999;
        for (Int_t k = 0; k < 2; k++)
            fDimuMu_rec[i][k] = 999;
        fDimuPhi_rec[i] = 999.;
        fDimuTheta_rec[i] = 999.;
    }

    fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
    fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
    fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
    DefineOutput(1, TTree::Class());
}

//___________________________________________________________________________
AliAnalysisTaskDimuon_HighMass &AliAnalysisTaskDimuon_HighMass::operator=(const AliAnalysisTaskDimuon_HighMass &c)
{
    //
    // Assignment operator
    //
    if (this != &c)
    {
        AliAnalysisTaskSE::operator=(c);
    }
    return *this;
}

//___________________________________________________________________________
AliAnalysisTaskDimuon_HighMass::AliAnalysisTaskDimuon_HighMass(const AliAnalysisTaskDimuon_HighMass &c) : AliAnalysisTaskSE(c),
                                                                                                          fBeamEnergy(c.fBeamEnergy),
                                                                                                          fkAnalysisType(c.fkAnalysisType),
                                                                                                          fPeriod(c.fPeriod),
                                                                                                          fAODEvent(c.fAODEvent),
                                                                                                          fOutputTree(c.fOutputTree),
                                                                                                          fMuonTrackCuts(c.fMuonTrackCuts),
                                                                                                          fNMuons_gen(c.fNMuons_gen),
                                                                                                          fNHadrons_gen(c.fNHadrons_gen),
                                                                                                          fNMuons_rec(c.fNMuons_rec),
                                                                                                          fNDimu_rec(c.fNDimu_rec),
                                                                                                          fPercentV0M(c.fPercentV0M),
                                                                                                          fVerbose(c.fVerbose)
{
    //
    // Copy Constructor
    //
}

//___________________________________________________________________________
AliAnalysisTaskDimuon_HighMass::~AliAnalysisTaskDimuon_HighMass()
{
    //
    // destructor
    //
    Info("~AliAnalysisTaskDimuon_HighMass", "Calling Destructor");
    if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis)
        delete fOutputTree;
}

//___________________________________________________________________________
void AliAnalysisTaskDimuon_HighMass::NotifyRun()
{
    fMuonTrackCuts->SetRun(fInputHandler);
}

//___________________________________________________________________________
void AliAnalysisTaskDimuon_HighMass::UserCreateOutputObjects()
{

    if (fOutputTree)
        return;

    OpenFile(1, "RECREATE");

    fOutputTree = new TTree("MCTree", "Data Tree");
    switch (fSaving_opt)
    {
    case kSaveHF:
        fOutputTree->Branch("N_HFquarks_gen", &fN_HFquarks_gen, "N_HFquarks_gen/I");
        fOutputTree->Branch("PDG_HFquark_gen", fPDG_HFquark_gen, "PDG_HFquark_gen[N_HFquarks_gen]/I");
        fOutputTree->Branch("Px_HFquark_gen", fPx_HFquark_gen, "Px_HFquark_gen[N_HFquarks_gen]/D");
        fOutputTree->Branch("Py_HFquark_gen", fPy_HFquark_gen, "Py_HFquark_gen[N_HFquarks_gen]/D");
        fOutputTree->Branch("Pz_HFquark_gen", fPz_HFquark_gen, "Pz_HFquark_gen[N_HFquarks_gen]/D");
        fOutputTree->Branch("Pt_HFquark_gen", fPt_HFquark_gen, "Pt_HFquark_gen[N_HFquarks_gen]/D");
        fOutputTree->Branch("Y_HFquark_gen", fY_HFquark_gen, "Y_HFquark_gen[N_HFquarks_gen]/D");
        fOutputTree->Branch("Mother_index", fMother_index, "Mother_index[N_HFquarks_gen]/I");

        fOutputTree->Branch("NHadrons_gen", &fNHadrons_gen, "NHadrons_gen/I");
        fOutputTree->Branch("PDGmum_Hadron_gen", fPDGmum_Hadron_gen, "PDGmum_Hadron_gen[NHadrons_gen]/I");
        fOutputTree->Branch("PDG_Hadron_gen", fPDG_Hadron_gen, "PDG_Hadron_gen[NHadrons_gen]/I");
        fOutputTree->Branch("Pt_Hadron_gen", fPt_Hadron_gen, "Pt_Hadron_gen[NHadrons_gen]/D");
        fOutputTree->Branch("E_Hadron_gen", fE_Hadron_gen, "E_Hadron_gen[NHadrons_gen]/D");
        fOutputTree->Branch("Px_Hadron_gen", fPx_Hadron_gen, "Px_Hadron_gen[NHadrons_gen]/D");
        fOutputTree->Branch("Py_Hadron_gen", fPy_Hadron_gen, "Py_Hadron_gen[NHadrons_gen]/D");
        fOutputTree->Branch("Pz_Hadron_gen", fPz_Hadron_gen, "Pz_Hadron_gen[NHadrons_gen]/D");
        fOutputTree->Branch("Y_Hadron_gen", fY_Hadron_gen, "Y_Hadron_gen[NHadrons_gen]/D");
        fOutputTree->Branch("Eta_Hadron_gen", fEta_Hadron_gen, "Eta_Hadron_gen[NHadrons_gen]/D");
        break;

    case kSaveDY:
        fOutputTree->Branch("N_gamma_gen", &fN_gamma, "N_gamma_gen/I");
        fOutputTree->Branch("Pt_gamma_gen", fPt_gamma, "Pt_gamma_gen[N_gamma_gen]/D");
        fOutputTree->Branch("M_gamma_gen", fM_gamma, "M_gamma_gen[N_gamma_gen]/D");
        fOutputTree->Branch("Y_gamma_gen", fY_gamma, "Y_gamma_gen[N_gamma_gen]/D");
        break;
    }

    fOutputTree->Branch("NMuons_gen", &fNMuons_gen, "NMuons_gen/I");
    fOutputTree->Branch("PDGmum_gen", fPDGmum_gen, "PDGmum_gen[NMuons_gen]/I");
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
    fOutputTree->Branch("From_Powheg_gen", fFrom_Powheg_gen, "From_Powheg_gen[NMuons_gen]/I");

    fOutputTree->Branch("NMuons_rec", &fNMuons_rec, "NMuons_rec/I");
    fOutputTree->Branch("PDGmum_rec", fPDGmum_rec, "PDGmum_rec[NMuons_rec]/I");
    fOutputTree->Branch("Pt_rec", fPt_rec, "Pt_rec[NMuons_rec]/D");
    fOutputTree->Branch("E_rec", fE_rec, "E_rec[NMuons_rec]/D");
    fOutputTree->Branch("Px_rec", fPx_rec, "Px_rec[NMuons_rec]/D");
    fOutputTree->Branch("Py_rec", fPy_rec, "Py_rec[NMuons_rec]/D");
    fOutputTree->Branch("Pz_rec", fPz_rec, "Pz_rec[NMuons_rec]/D");
    fOutputTree->Branch("Y_rec", fY_rec, "Y_rec[NMuons_rec]/D");
    fOutputTree->Branch("Eta_rec", fEta_rec, "Eta_rec[NMuons_rec]/D");
    fOutputTree->Branch("MatchTrig_rec", fMatchTrig_rec, "MatchTrig_rec[NMuons_rec]/I");
    fOutputTree->Branch("TrackChi2_rec", fTrackChi2_rec, "TrackChi2_rec[NMuons_rec]/D");
    fOutputTree->Branch("MatchTrigChi2_rec", fMatchTrigChi2_rec, "MatchTrigChi2_rec[NMuons_rec]/D");
    fOutputTree->Branch("Charge_rec", fCharge_rec, "Charge_rec[NMuons_rec]/I");
    fOutputTree->Branch("RAtAbsEnd_rec", fRAtAbsEnd_rec, "RAtAbsEnd_rec[NMuons_rec]/D");
    fOutputTree->Branch("pDCA_rec", fpDCA_rec, "pDCA[NMuons_rec]/I");
    fOutputTree->Branch("Phi_rec", fPhi_rec, "Phi_rec[NMuons_rec]/D");
    fOutputTree->Branch("Theta_rec", fTheta_rec, "Theta_rec[NMuons_rec]/D");
    fOutputTree->Branch("From_Powheg_rec", fFrom_Powheg_rec, "From_Powheg_rec[NMuons_rec]/I");

    fOutputTree->Branch("NDimu_gen", &fNDimu_gen, "NDimu_gen/I");
    fOutputTree->Branch("DimuMu_gen", fDimuMu_gen, "DimuMu_gen[NDimu_gen][2]/I");
    fOutputTree->Branch("DimuPt_gen", fDimuPt_gen, "DimuPt_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuPx_gen", fDimuPx_gen, "DimuPx_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuPy_gen", fDimuPy_gen, "DimuPy_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuPz_gen", fDimuPz_gen, "DimuPz_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuY_gen", fDimuY_gen, "DimuY_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuMass_gen", fDimuMass_gen, "DimuMass_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuCharge_gen", fDimuCharge_gen, "DimuCharge_gen[NDimu_gen]/I");

    fOutputTree->Branch("NDimu_rec", &fNDimu_rec, "NDimu_rec/I");
    fOutputTree->Branch("DimuMu_rec", fDimuMu_rec, "DimuMu_rec[NDimu_rec][2]/I");
    fOutputTree->Branch("DimuPt_rec", fDimuPt_rec, "DimuPt_rec[NDimu_rec]/D");
    fOutputTree->Branch("DimuPx_rec", fDimuPx_rec, "DimuPx_rec[NDimu_rec]/D");
    fOutputTree->Branch("DimuPy_rec", fDimuPy_rec, "DimuPy_rec[NDimu_rec]/D");
    fOutputTree->Branch("DimuPz_rec", fDimuPz_rec, "DimuPz_rec[NDimu_rec]/D");
    fOutputTree->Branch("DimuY_rec", fDimuY_rec, "DimuY_rec[NDimu_rec]/D");
    fOutputTree->Branch("DimuMass_rec", fDimuMass_rec, "DimuMass_rec[NDimu_rec]/D");
    fOutputTree->Branch("DimuCharge_rec", fDimuCharge_rec, "DimuCharge_rec[NDimu_rec]/I");
    fOutputTree->Branch("DimuMatch_rec", fDimuMatch_rec, "DimuMatch_rec[NDimu_rec]/I");
    fOutputTree->Branch("DimuPhi_rec", fDimuPhi_rec, "DimuPhi_rec[NDimu_rec]/D");
    fOutputTree->Branch("DimuTheta_rec", fDimuTheta_rec, "DimuTheta_rec[NDimu_rec]/D");

    fOutputTree->ls();

    PostData(1, fOutputTree);
}

//_________________________________________________
void AliAnalysisTaskDimuon_HighMass::UserExec(Option_t *)
{

    // printf("Entro in UserExec\n" );
    fNMuons_gen = 0;
    fNHadrons_gen = 0;
    fN_gamma = 0;
    fNDimu_gen = 0;
    fNMuons_rec = 0;
    fNDimu_rec = 0;
    fN_HFquarks_gen = 0;
    fN_HFquarks_rec = 0;
    for (Int_t i = 0; i < fMuons_dim; i++)
    {
        fPDG_HFquark_gen[i] = 9999.;
        fPx_HFquark_gen[i] = 9999.;
        fPy_HFquark_gen[i] = 9999.;
        fPz_HFquark_gen[i] = 9999.;
        fPt_HFquark_gen[i] = 9999.;
        fY_HFquark_gen[i] = 9999.;
        fMother_index[i] = 9999;

        fPt_gamma[i] = 9999.;
        fM_gamma[i] = 9999.;
        fY_gamma[i] = 9999.;

        fPDGmum_gen[i] = 9999.;
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
        fFrom_Powheg_gen[i] = 999;

        fPDGmum_Hadron_gen[i] = 999; // gen Hadron PDG mum
        fPDG_Hadron_gen[i] = 999;    // gen Hadron PDG
        fPt_Hadron_gen[i] = 999;     // gen Hadron pT
        fE_Hadron_gen[i] = 999;      // gen Hadron E
        fPx_Hadron_gen[i] = 999;     // gen Hadron px
        fPy_Hadron_gen[i] = 999;     // gen Hadron py
        fPz_Hadron_gen[i] = 999;     // gen Hadron pz
        fY_Hadron_gen[i] = 999;      // gen Hadron y
        fEta_Hadron_gen[i] = 999;    // gen Hadron eta

        fPDGmum_rec[i] = 9999.;
        fPt_rec[i] = 999.;
        fE_rec[i] = 999.;
        fPx_rec[i] = 999;
        fPy_rec[i] = 999;
        fPz_rec[i] = 999;
        fY_rec[i] = 999.;
        fEta_rec[i] = 999.;
        fMatchTrig_rec[i] = 999.;
        fTrackChi2_rec[i] = 999.;
        fMatchTrigChi2_rec[i] = 999.;
        fCharge_rec[i] = 999;
        fRAtAbsEnd_rec[i] = 999;
        fpDCA_rec[i] = 999.;
        fPhi_rec[i] = 999.;
        fTheta_rec[i] = 999.;
        fFrom_Powheg_rec[i] = 999;
    }
    for (Int_t i = 0; i < fDimu_dim; i++)
    {
        fDimuPt_gen[i] = 999.;
        fDimuPx_gen[i] = 999.;
        fDimuPy_gen[i] = 999.;
        fDimuPz_gen[i] = 999.;
        fDimuY_gen[i] = 999.;
        fDimuMass_gen[i] = 999.;
        fDimuCharge_gen[i] = 999.;

        fDimuPt_rec[i] = 999.;
        fDimuPx_rec[i] = 999.;
        fDimuPy_rec[i] = 999.;
        fDimuPz_rec[i] = 999.;
        fDimuY_rec[i] = 999.;
        fDimuMass_rec[i] = 999.;
        fDimuCharge_rec[i] = 999.;
        for (Int_t k = 0; k < 2; k++)
            fDimuMu_rec[i][k] = 999;
        fDimuPhi_rec[i] = 999.;
        fDimuTheta_rec[i] = 999.;
    }

    //---------------------------------------------------
    // Execute analysis for current event
    //---------------------------------------------------

    fAODEvent = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!fAODEvent)
    {
        AliError("AOD event not found. Nothing done!");
        return;
    }

    //-----------------------------------------------
    // loop on MC generated event
    //-----------------------------------------------
    TClonesArray *mcarray = dynamic_cast<TClonesArray *>(fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    AliAODHeader *aodheader = dynamic_cast<AliAODHeader *>(fAODEvent->GetHeader());
    TString firedtrigger = aodheader->GetFiredTriggerClasses();
    printf("%s\n", firedtrigger.Data());

    Int_t nHFquark_gen = 0;
    Int_t nmu_gen = 0;
    Int_t nHadron_gen = 0;
    Int_t ndimu_gen = 0;
    Int_t n_gamma = 0;

    Int_t n_MC_part = mcarray->GetEntries();
    Int_t *LabelOld1 = new Int_t[n_MC_part];
    Int_t *LabelOld2 = new Int_t[n_MC_part];
    Bool_t *GoodMuon = new Bool_t[n_MC_part];
    for (Int_t i_nMCpart = 0; i_nMCpart < n_MC_part; i_nMCpart++)
    {
        LabelOld1[i_nMCpart] = 999;
        LabelOld2[i_nMCpart] = 999;
        GoodMuon[i_nMCpart] = kFALSE;
    }

    for (Int_t i_nMCpart = 0; i_nMCpart < n_MC_part; i_nMCpart++)
    {

        AliAODMCParticle *MC_part0 = (AliAODMCParticle *)mcarray->At(i_nMCpart);
        Int_t PDG_MC_part0 = MC_part0->GetPdgCode();
        if (MC_part0->Y() < -15.0 || MC_part0->Y() > 15.0)
            continue;

        Bool_t Muon = kFALSE;
        Bool_t Gamma_star = kFALSE;
        Bool_t HF_quark = kFALSE;
        Bool_t Charm_Hadron = kFALSE;
        Bool_t Beauty_Hadron = kFALSE;

        if (TMath::Abs(PDG_MC_part0) == 13)
            Muon = kTRUE;
        else if (TMath::Abs(PDG_MC_part0) == 4 || TMath::Abs(PDG_MC_part0) == 5)
            HF_quark = kTRUE;
        else if (TMath::Abs(PDG_MC_part0) == 23)
            Gamma_star = kTRUE;
        else if ((TMath::Abs(PDG_MC_part0) > 400 && TMath::Abs(PDG_MC_part0) < 500) || (TMath::Abs(PDG_MC_part0) > 4000 && TMath::Abs(PDG_MC_part0) < 5000))
            Charm_Hadron = kTRUE;
        else if ((TMath::Abs(PDG_MC_part0) > 500 && TMath::Abs(PDG_MC_part0) < 600) || (TMath::Abs(PDG_MC_part0) > 5000 && TMath::Abs(PDG_MC_part0) < 6000))
            Beauty_Hadron = kTRUE;

        if (!Muon && !Charm_Hadron && !Beauty_Hadron && !HF_quark & !Gamma_star)
            continue;
        Int_t index_Mum_MC_part0 = MC_part0->GetMother();
        AliAODMCParticle *Mum_MC_part0;
        Int_t PDG_Mum_MC_part0 = 999;
        Int_t final_mom_mu0 = 999;

        if (Charm_Hadron && index_Mum_MC_part0 == -1)
            PDG_Mum_MC_part0 = -TMath::Sign(1, PDG_MC_part0) * 4;
        else if (Beauty_Hadron && index_Mum_MC_part0 == -1)
            PDG_Mum_MC_part0 = -TMath::Sign(1, PDG_MC_part0) * 5;
        else if (HF_quark && index_Mum_MC_part0 == -1)
            PDG_Mum_MC_part0 = 999;
        else if (Gamma_star && index_Mum_MC_part0 == -1)
            PDG_Mum_MC_part0 = 999;
        else
        {
            Mum_MC_part0 = (AliAODMCParticle *)mcarray->At(index_Mum_MC_part0);
            PDG_Mum_MC_part0 = Mum_MC_part0->GetPdgCode();
        }

        if (TMath::Abs(PDG_Mum_MC_part0) == TMath::Abs(PDG_MC_part0)) // Remove particles formally produced by the same particle (Powheg transport)
            continue;
        TLorentzVector vector_mcp0;

        if (Muon)
        {
            if ((PDG_Mum_MC_part0 != 23) && !(TMath::Abs(PDG_Mum_MC_part0) > 400 && TMath::Abs(PDG_Mum_MC_part0) < 500) && !(TMath::Abs(PDG_Mum_MC_part0) > 4000 && TMath::Abs(PDG_Mum_MC_part0) < 5000) && !(TMath::Abs(PDG_Mum_MC_part0) > 500 && TMath::Abs(PDG_Mum_MC_part0) < 600) && !(TMath::Abs(PDG_Mum_MC_part0) > 5000 && TMath::Abs(PDG_Mum_MC_part0) < 6000))
                continue;

            if (fVerbose)
            {
                printf("Muon mother: %d\n", PDG_Mum_MC_part0);
                printf("i_nMCpart %d) PDG: %d | MCParticle Pt %0.1f | Px %0.1f | Py %0.1f | Pz %0.1f | Y %0.1f || Mother %d , PDG Mother %d \n", i_nMCpart, PDG_MC_part0, MC_part0->Pt(), MC_part0->Px(), MC_part0->Py(), MC_part0->Pz(), MC_part0->Y(), index_Mum_MC_part0, PDG_Mum_MC_part0);
            }
            fY_gen[nmu_gen] = MC_part0->Y();
            fPt_gen[nmu_gen] = MC_part0->Pt();
            fE_gen[nmu_gen] = MC_part0->E();
            fPx_gen[nmu_gen] = MC_part0->Px();
            fPy_gen[nmu_gen] = MC_part0->Py();
            fPz_gen[nmu_gen] = MC_part0->Pz();
            fEta_gen[nmu_gen] = MC_part0->Eta();
            fPhi_gen[nmu_gen] = MC_part0->Phi();
            fCharge_gen[nmu_gen] = -TMath::Sign(1, PDG_MC_part0); // Sign minus necessary
            fTheta_gen[nmu_gen] = MC_part0->Theta();
            MC_part0->Momentum(vector_mcp0);
            Int_t Final_PDG = 999;
            Bool_t isFromPowheg = kFALSE;
            // Int_t final_mom_mu0 = IsPrompt(Mum_MC_part0, mcarray);
            // if ((TMath::Abs(PDG_Mum_MC_part0) > 400 && TMath::Abs(PDG_Mum_MC_part0) < 500) || (TMath::Abs(PDG_Mum_MC_part0) > 4000 && TMath::Abs(PDG_Mum_MC_part0) < 5000))
            //     fPDGmum_gen[nmu_gen] = IsPrompt(Mum_MC_part0, mcarray);
            // else
            //     fPDGmum_gen[nmu_gen] = PDG_Mum_MC_part0;
            Int_t *muon_properties = Muon_ancestor_features(Mum_MC_part0, mcarray, fAOD_origin);

            printf("Gen Muon properties: PDG mum %d , is from powheg %d\n", muon_properties[0], muon_properties[1]);
            fPDGmum_gen[nmu_gen] = muon_properties[0];
            fFrom_Powheg_gen[nmu_gen] = muon_properties[1];

            // printf("Final Muon mother: %d\n", fPDGmum_gen[nmu_gen]);
            nmu_gen++;
        }
        else if (Gamma_star)
        {
            if (fVerbose)
            {
                printf("i_nMCpart %d) PDG: %d | MCParticle Pt %0.1f | Px %0.1f | Py %0.1f | Pz %0.1f | Y %0.1f || Mother %d , PDG Mother %d \n", i_nMCpart, PDG_MC_part0, MC_part0->Pt(), MC_part0->Px(), MC_part0->Py(), MC_part0->Pz(), MC_part0->Y(), index_Mum_MC_part0, PDG_Mum_MC_part0);
            }

            fPt_gamma[n_gamma] = MC_part0->Pt();
            fM_gamma[n_gamma] = MC_part0->GetCalcMass();
            fY_gamma[n_gamma] = MC_part0->Y();
            if (index_Mum_MC_part0 == -1 && fSaving_opt == kSaveDY)
            {
                printf("\n");
                printf("index %d\n", i_nMCpart);
                MC_part0->Print();
                printf("\n");
            }
            n_gamma++;
        }
        else if (HF_quark)
        {
            if (fVerbose)
            {
                printf("Flag %u Gen status %d\n", MC_part0->GetFlag(), MC_part0->GetGeneratorIndex());
                printf("i_nMCpart %d) PDG: %d | MCParticle Pt %0.1f | Px %0.1f | Py %0.1f | Pz %0.1f | Y %0.1f || Mother %d , PDG Mother %d \n", i_nMCpart, PDG_MC_part0, MC_part0->Pt(), MC_part0->Px(), MC_part0->Py(), MC_part0->Pz(), MC_part0->Y(), index_Mum_MC_part0, PDG_Mum_MC_part0);
            }

            fPDG_HFquark_gen[nHFquark_gen] = PDG_MC_part0;
            fPx_HFquark_gen[nHFquark_gen] = MC_part0->Px();
            fPy_HFquark_gen[nHFquark_gen] = MC_part0->Py();
            fPz_HFquark_gen[nHFquark_gen] = MC_part0->Pz();
            fPt_HFquark_gen[nHFquark_gen] = MC_part0->Pt();
            fY_HFquark_gen[nHFquark_gen] = MC_part0->Y();
            fMother_index[nHFquark_gen] = index_Mum_MC_part0;

            if (index_Mum_MC_part0 == -1 && fSaving_opt == kSaveHF)
            {
                printf("\n");
                printf("index %d\n", i_nMCpart);
                MC_part0->Print();
                printf("\n");
            }

            nHFquark_gen++;
        }

        else if (Charm_Hadron || Beauty_Hadron)
        {
            // MC_part0->Print();
            // printf("i_nMCpart %d) PDG: %d | MCParticle Pt %0.1f | Px %0.1f | Py %0.1f | Pz %0.1f | Y %0.1f || Mother %d , PDG Mother %d \n", i_nMCpart, PDG_MC_part0, MC_part0->Pt(), MC_part0->Px(), MC_part0->Py(), MC_part0->Pz(), MC_part0->Y(), index_Mum_MC_part0, PDG_Mum_MC_part0);

            if (Charm_Hadron && TMath::Abs(PDG_Mum_MC_part0) != 4)
                final_mom_mu0 = IsPrompt(Mum_MC_part0, mcarray);
            else
                final_mom_mu0 = PDG_Mum_MC_part0;
            // printf("final mum %i \n", final_mom_mu0);
            fPDGmum_Hadron_gen[nHadron_gen] = final_mom_mu0; // gen Hadron PDG mum
            fPDG_Hadron_gen[nHadron_gen] = PDG_MC_part0;     // gen Hadron PDG
            fPt_Hadron_gen[nHadron_gen] = MC_part0->Pt();
            fE_Hadron_gen[nHadron_gen] = MC_part0->E();
            fPx_Hadron_gen[nHadron_gen] = MC_part0->Px();
            fPy_Hadron_gen[nHadron_gen] = MC_part0->Py();
            fPz_Hadron_gen[nHadron_gen] = MC_part0->Pz();
            fY_Hadron_gen[nHadron_gen] = MC_part0->Y();
            fEta_Hadron_gen[nHadron_gen] = MC_part0->Eta();
            nHadron_gen++;
        }

        if (TMath::Abs(PDG_MC_part0) != 13)
            continue;
        for (Int_t j_nMCpart1 = i_nMCpart + 1; j_nMCpart1 < n_MC_part; j_nMCpart1++)
        {
            AliAODMCParticle *MC_part1 = (AliAODMCParticle *)mcarray->At(j_nMCpart1);
            Int_t PDG_MC_part1 = MC_part1->GetPdgCode();
            if (TMath::Abs(PDG_MC_part1) != 13)
                continue;
            if (MC_part1->Y() < -15.0 || MC_part1->Y() > 15.0)
                continue;
            Int_t index_Mum_MC_part1 = MC_part1->GetMother();
            AliAODMCParticle *Mum_MC_part1 = (AliAODMCParticle *)mcarray->At(index_Mum_MC_part1);
            Int_t PDG_Mum_MC_part1 = Mum_MC_part1->GetPdgCode();
            Int_t final_Mum_MC_part1 = PDG_Mum_MC_part1;

            if ((PDG_Mum_MC_part1 != 23) && !(TMath::Abs(PDG_Mum_MC_part1) > 400 && TMath::Abs(PDG_Mum_MC_part1) < 500) && !(TMath::Abs(PDG_Mum_MC_part1) > 4000 && TMath::Abs(PDG_Mum_MC_part1) < 5000) && !(TMath::Abs(PDG_Mum_MC_part1) > 500 && TMath::Abs(PDG_Mum_MC_part1) < 600) && !(TMath::Abs(PDG_Mum_MC_part1) > 5000 && TMath::Abs(PDG_Mum_MC_part1) < 6000))
                continue;
            TLorentzVector vector_mcp1;
            MC_part1->Momentum(vector_mcp1);

            TLorentzVector dimu = vector_mcp0 + vector_mcp1;

            fDimuPt_gen[ndimu_gen] = dimu.Pt();
            fDimuPx_gen[ndimu_gen] = dimu.Px();
            fDimuPy_gen[ndimu_gen] = dimu.Py();
            fDimuPz_gen[ndimu_gen] = dimu.Pz();
            fDimuY_gen[ndimu_gen] = dimu.Rapidity();
            fDimuMass_gen[ndimu_gen] = dimu.M();
            fDimuCharge_gen[ndimu_gen] = -(TMath::Sign(1, PDG_MC_part0) + TMath::Sign(1, PDG_MC_part1));
            LabelOld1[ndimu_gen] = i_nMCpart;
            LabelOld2[ndimu_gen] = j_nMCpart1;
            GoodMuon[i_nMCpart] = kTRUE;
            GoodMuon[j_nMCpart1] = kTRUE;
            ndimu_gen++;
        }
    }
    fN_gamma = n_gamma;
    fN_HFquarks_gen = nHFquark_gen;
    fNMuons_gen = nmu_gen;
    fNHadrons_gen = nHadron_gen;
    fNDimu_gen = ndimu_gen;
    for (Int_t i_nDimu_Gen = 0; i_nDimu_Gen < ndimu_gen; i_nDimu_Gen++)
    {
        Int_t LabelNew1 = 0;
        Int_t LabelNew2 = 0;
        for (int j = 0; j < LabelOld1[i_nDimu_Gen]; j++)
        {
            if (GoodMuon[j])
                LabelNew1++;
        }
        for (int j = 0; j < LabelOld2[i_nDimu_Gen]; j++)
        {
            if (GoodMuon[j])
                LabelNew2++;
        }
        fDimuMu_gen[i_nDimu_Gen][0] = LabelNew1;
        fDimuMu_gen[i_nDimu_Gen][1] = LabelNew2;
    }

    Int_t numtracks = fAODEvent->GetNumberOfTracks();
    // printf("Number of tracks per event %i\n", numtracks);
    Int_t nHFquark_rec = 0;

    Int_t nmu_rec = 0;
    Int_t ndimu_rec = 0;

    Int_t *LabelOld1_rec = new Int_t[numtracks];
    Int_t *LabelOld2_rec = new Int_t[numtracks];
    Bool_t *GoodMuon_rec = new Bool_t[numtracks];

    for (Int_t i_Track0 = 0; i_Track0 < numtracks; i_Track0++)
    {
        LabelOld1_rec[i_Track0] = 999;
        LabelOld2_rec[i_Track0] = 999;
        GoodMuon_rec[i_Track0] = kFALSE;
    }

    for (Int_t i_Track0 = 0; i_Track0 < numtracks; i_Track0++)
    {
        GoodMuon_rec[i_Track0] = kFALSE;
        AliAODTrack *Track0 = (AliAODTrack *)fAODEvent->GetTrack(i_Track0);

        if (Track0->GetLabel() < 0)
            continue;

        if (Track0->Eta() < -4.5 || Track0->Eta() > -2.0) // Condition on muon acceptance
            continue;

        AliAODMCParticle *MCpart_Track0 = (AliAODMCParticle *)mcarray->At(Track0->GetLabel());
        Int_t PDG_MCpart_Track0 = MCpart_Track0->GetPdgCode();

        if (TMath::Abs(PDG_MCpart_Track0) != 13) // request of a muon
            continue;

        Int_t index_Mum_MCpart_Track0 = MCpart_Track0->GetMother();
        AliAODMCParticle *Mum_MCpart_Track0 = (AliAODMCParticle *)mcarray->At(index_Mum_MCpart_Track0);
        Int_t PDG_Mum_MCpart_Track0 = Mum_MCpart_Track0->GetPdgCode();

        // MCpart_Track0->Print();

        Int_t round = 0;
        TString str = Form("pos %d) PDG particle rec %d [Pt %0.03f,Y %0.03f]==>> (round %d)PDG MUM particle rec %d [pos %d, Pt %0.3f, Y %0.03f]", i_Track0, PDG_MCpart_Track0, MCpart_Track0->Pt(), MCpart_Track0->Y(), round, PDG_Mum_MCpart_Track0, index_Mum_MCpart_Track0, Mum_MCpart_Track0->Pt(), Mum_MCpart_Track0->Y());

        if (TMath::Abs(PDG_Mum_MCpart_Track0) == TMath::Abs(PDG_MCpart_Track0))
        {
            Bool_t testing_Track0 = kTRUE;
            Int_t final_mum_Track0 = 999;
            while (testing_Track0)
            {
                //        printf("Round %d | Index Muon Ancestor %d | PDG Muon Ancestor %d | \n",round,ancestor_mum0, PDGcode_ancestor0);
                round++;
                if (index_Mum_MCpart_Track0 < 2)
                {
                    testing_Track0 = kFALSE;
                    break;
                }

                else
                {
                    index_Mum_MCpart_Track0 = Mum_MCpart_Track0->GetMother();
                    if (index_Mum_MCpart_Track0 > 0)
                        Mum_MCpart_Track0 = (AliAODMCParticle *)mcarray->At(index_Mum_MCpart_Track0); // Salto gen Muone

                    final_mum_Track0 = Mum_MCpart_Track0->GetPdgCode();
                    str += Form(" ==>> (round %d)PDG MUM particle rec %d [pos %d, Pt %0.3f, Y %0.03f]", round, PDG_Mum_MCpart_Track0, index_Mum_MCpart_Track0, Mum_MCpart_Track0->Pt(), Mum_MCpart_Track0->Y());
                    if (TMath::Abs(final_mum_Track0) != 13)
                    {
                        PDG_Mum_MCpart_Track0 = final_mum_Track0;
                        testing_Track0 = kFALSE;
                    }
                }
            }
        } // Remove particles formally produced by the same particle (Powheg transport)
        if (fVerbose)
            printf("%s\n", str.Data());
        if ((PDG_Mum_MCpart_Track0 == 23) || (TMath::Abs(PDG_Mum_MCpart_Track0) > 400 && TMath::Abs(PDG_Mum_MCpart_Track0) < 500) || (TMath::Abs(PDG_Mum_MCpart_Track0) > 4000 && TMath::Abs(PDG_Mum_MCpart_Track0) < 5000) || (TMath::Abs(PDG_Mum_MCpart_Track0) > 500 && TMath::Abs(PDG_Mum_MCpart_Track0) < 600) || (TMath::Abs(PDG_Mum_MCpart_Track0) > 5000 && TMath::Abs(PDG_Mum_MCpart_Track0) < 6000))
            GoodMuon_rec[i_Track0] = kTRUE;

        if (!GoodMuon_rec[i_Track0])
            continue;

        fPt_rec[nmu_rec] = Track0->Pt();
        fPx_rec[nmu_rec] = Track0->Px();
        fPy_rec[nmu_rec] = Track0->Py();
        fPz_rec[nmu_rec] = Track0->Pz();
        fY_rec[nmu_rec] = Track0->Y();
        fEta_rec[nmu_rec] = Track0->Eta();
        fE_rec[nmu_rec] = Track0->E();
        fPhi_rec[nmu_rec] = Track0->Phi();
        fTheta_rec[nmu_rec] = Track0->Theta();
        fCharge_rec[nmu_rec] = Track0->Charge();
        fMatchTrig_rec[nmu_rec] = Track0->GetMatchTrigger();
        fMatchTrigChi2_rec[nmu_rec] = Track0->GetChi2MatchTrigger();
        fRAtAbsEnd_rec[nmu_rec] = Track0->GetRAtAbsorberEnd();

        // if ((TMath::Abs(PDG_Mum_MCpart_Track0) > 400 && TMath::Abs(PDG_Mum_MCpart_Track0) < 500) || (TMath::Abs(PDG_Mum_MCpart_Track0) > 4000 && TMath::Abs(PDG_Mum_MCpart_Track0) < 5000))
        // {
        //     fPDGmum_rec[nmu_rec] = IsPrompt(Mum_MCpart_Track0, mcarray);
        // }
        // else
        //     fPDGmum_rec[nmu_rec] = PDG_Mum_MCpart_Track0;
        Int_t *muon_properties = Muon_ancestor_features(Mum_MCpart_Track0, mcarray, fAOD_origin);

        printf("Rec Muon properties: PDG mum %d , is from powheg %d\n", muon_properties[0], muon_properties[1]);
        fPDGmum_rec[nmu_rec] = muon_properties[0];
        fFrom_Powheg_rec[nmu_rec] = muon_properties[1];

        if (fMuonTrackCuts->IsSelected(Track0))
            fpDCA_rec[nmu_rec] = 1;

        nmu_rec++;
        for (Int_t j_Track1 = i_Track0 + 1; j_Track1 < numtracks; j_Track1++)
        {
            AliAODTrack *Track1 = (AliAODTrack *)fAODEvent->GetTrack(j_Track1);
            if (Track1->GetLabel() < 0)
                continue;
            AliAODMCParticle *MCpart_Track1 = (AliAODMCParticle *)mcarray->At(Track1->GetLabel());
            Int_t PDG_mctrack1 = MCpart_Track1->GetPdgCode();
            if (TMath::Abs(PDG_mctrack1) != 13)
                continue;
            if (Track1->Eta() < -4.5 || Track1->Eta() > -2.0)
                continue;

            Int_t index_Mum_MCpart_Track1 = MCpart_Track1->GetMother();
            AliAODMCParticle *Mum_MCpart_Track1 = (AliAODMCParticle *)mcarray->At(index_Mum_MCpart_Track1);
            Int_t PDG_Mum_MCpart_Track1 = Mum_MCpart_Track1->GetPdgCode();

            if (TMath::Abs(PDG_Mum_MCpart_Track1) == TMath::Abs(PDG_mctrack1))
            {
                Int_t final_mum_Track1 = 999;
                Bool_t testing_Track1 = kTRUE;
                while (testing_Track1)
                {
                    //        printf("Round %d | Index Muon Ancestor %d | PDG Muon Ancestor %d | \n",round,ancestor_mum0, PDGcode_ancestor0);
                    round++;
                    if (index_Mum_MCpart_Track1 < 2)
                    {
                        testing_Track1 = kFALSE;
                        break;
                    }

                    else
                    {
                        index_Mum_MCpart_Track1 = Mum_MCpart_Track1->GetMother();
                        if (index_Mum_MCpart_Track1 > 0)
                            Mum_MCpart_Track1 = (AliAODMCParticle *)mcarray->At(index_Mum_MCpart_Track1); // Salto gen Muone

                        final_mum_Track1 = Mum_MCpart_Track1->GetPdgCode();
                        if (TMath::Abs(final_mum_Track1) != 13)
                        {
                            PDG_Mum_MCpart_Track1 = final_mum_Track1;
                            testing_Track1 = kFALSE;
                        }
                    }
                }
            } // Remove particles formally produced by the same particle (Powheg transport)

            if ((PDG_Mum_MCpart_Track1 == 23) || (TMath::Abs(PDG_Mum_MCpart_Track1) > 400 && TMath::Abs(PDG_Mum_MCpart_Track1) < 500) || (TMath::Abs(PDG_Mum_MCpart_Track1) > 4000 && TMath::Abs(PDG_Mum_MCpart_Track1) < 5000) || (TMath::Abs(PDG_Mum_MCpart_Track1) > 500 && TMath::Abs(PDG_Mum_MCpart_Track1) < 600) || (TMath::Abs(PDG_Mum_MCpart_Track1) > 5000 && TMath::Abs(PDG_Mum_MCpart_Track1) < 6000))
                GoodMuon_rec[j_Track1] = kTRUE;

            if (!GoodMuon_rec[j_Track1])
                continue;
            AliAODDimuon *dimu = new AliAODDimuon(Track0, Track1);
            fDimuMass_rec[ndimu_rec] = dimu->Mass();
            fDimuPt_rec[ndimu_rec] = dimu->Pt();
            fDimuPx_rec[ndimu_rec] = dimu->Px();
            fDimuPy_rec[ndimu_rec] = dimu->Py();
            fDimuPz_rec[ndimu_rec] = dimu->Pz();
            fDimuY_rec[ndimu_rec] = dimu->Y();
            fDimuCharge_rec[ndimu_rec] = dimu->Charge();
            fDimuPhi_rec[ndimu_rec] = dimu->Phi();
            fDimuTheta_rec[ndimu_rec] = dimu->Theta();

            // printf("Charge Dimu %d || ChargeMu0 %d ChargeMu1 %d \n",dimu->Charge(),mu0->Charge(),mu1->Charge());
            if (Track0->GetMatchTrigger() > 1 || Track1->GetMatchTrigger() > 1)
                fDimuMatch_rec[ndimu_rec] = 1;
            if (Track0->GetMatchTrigger() > 1 && Track1->GetMatchTrigger() > 1)
                fDimuMatch_rec[ndimu_rec] = 2;

            LabelOld1_rec[ndimu_rec] = i_Track0;
            LabelOld2_rec[ndimu_rec] = j_Track1;
            delete dimu;
            ndimu_rec++;
        }
    }
    fN_HFquarks_rec = nHFquark_rec;
    fNMuons_rec = nmu_rec;
    fNDimu_rec = ndimu_rec;
    for (Int_t i = 0; i < ndimu_rec; i++)
    {
        Int_t LabelNew1 = 0;
        Int_t LabelNew2 = 0;

        for (Int_t j = 0; j < LabelOld1_rec[i]; j++)
        {
            if (GoodMuon_rec[j])
                LabelNew1++;
        }
        for (Int_t j_LabelOld2_rec = 0; j_LabelOld2_rec < LabelOld2_rec[i]; j_LabelOld2_rec++)
        {
            if (GoodMuon_rec[j_LabelOld2_rec])
                LabelNew2++;
        }
        fDimuMu_rec[i][0] = LabelNew1;
        fDimuMu_rec[i][1] = LabelNew2;
    }

    fOutputTree->Fill();
    PostData(1, fOutputTree);
}

//________________________________________________________________________

Int_t *Muon_ancestor_features(AliAODMCParticle *mcp_mumum, TClonesArray *mcarray, TString AOD_origin)
{
    static Int_t Muon_properties[2];
    Muon_properties[0] = mcp_mumum->GetPdgCode();
    Muon_properties[1] = 999;
    Bool_t Verbosity = kTRUE;
    Int_t index_ancestor_mum = mcp_mumum->GetMother();
    if (index_ancestor_mum == -1)
    {
        if (AOD_origin.Contains("Powheg"))
            Muon_properties[1] = index_ancestor_mum;
        else
            Muon_properties[1] = 999;
        return Muon_properties;
    }
    Int_t PDG_muon_mum = mcp_mumum->GetPdgCode();

    if ((TMath::Abs(PDG_muon_mum) > 400 && TMath::Abs(PDG_muon_mum) < 500) || (TMath::Abs(PDG_muon_mum) > 4000 && TMath::Abs(PDG_muon_mum) < 5000))
    {
        if (Verbosity)
            printf("Starting: Index Muon Ancestor %d | PDG Muon Ancestor %d | \n", index_ancestor_mum, PDG_muon_mum);
        Bool_t testing = kTRUE;
        Int_t round = 0;

        AliAODMCParticle *mcp_ancestor0 = (AliAODMCParticle *)mcarray->At(index_ancestor_mum); // Muon ancestor definition
        Int_t PDGcode_ancestor0 = mcp_ancestor0->GetPdgCode();
        if (Verbosity)
            printf("Round: %d Index Muon Ancestor %d | PDG Muon Ancestor %d | \n", round, index_ancestor_mum, PDGcode_ancestor0);

        while (testing)
        {
            if (Verbosity)
                printf("Round %d | Index Muon Ancestor %d | PDG Muon Ancestor %d | \n", round, index_ancestor_mum, PDGcode_ancestor0);
            round++;
            if (index_ancestor_mum < 2)
            {
                testing = kFALSE;
                break;
            }
            else if (TMath::Abs(PDGcode_ancestor0) == 5)
            {
                if (Verbosity)
                    printf("Muon from Beauty Quark \n");
                testing = kFALSE;
                break;
            }
            else if (TMath::Abs(PDGcode_ancestor0) > 500 && TMath::Abs(PDGcode_ancestor0) < 600)
            {
                if (Verbosity)
                    printf("Muon from Beauty Meson \n");
                PDG_muon_mum = PDGcode_ancestor0; // if the charm hadron ancestor had a beauty hadron as ancestor, save this last one as the muon ancestor
                testing = kFALSE;
                break;
            }
            else if (TMath::Abs(PDGcode_ancestor0) > 5000 && TMath::Abs(PDGcode_ancestor0) < 6000)
            {
                if (Verbosity)
                    printf("Muon from Beauty Barion \n");
                PDG_muon_mum = PDGcode_ancestor0; // if the charm hadron ancestor had a beauty hadron as ancestor, save this last one as the muon ancestor
                testing = kFALSE;
                break;
            }
            else
            {
                index_ancestor_mum = mcp_ancestor0->GetMother();
                if (index_ancestor_mum > 0)
                {
                    mcp_ancestor0 = (AliAODMCParticle *)mcarray->At(index_ancestor_mum); // Salto gen Muone
                    PDGcode_ancestor0 = mcp_ancestor0->GetPdgCode();
                }
                else
                {
                    mcp_ancestor0 = (AliAODMCParticle *)mcarray->At(mcp_mumum->GetMother());
                    testing = kFALSE;
                }
            }
        }
    }
    Muon_properties[0] = PDG_muon_mum;

    if (AOD_origin != "Powheg")
        return Muon_properties;

    if (Verbosity)
        printf("Looking for powheg HF quark\n");

    if (Verbosity)
        printf("Starting: Index Muon Ancestor %d | PDG Muon Ancestor %d | \n", index_ancestor_mum, PDG_muon_mum);
    Bool_t testing = kTRUE;
    Int_t round = 0;

    if (index_ancestor_mum < 2)
        Muon_properties[1] = index_ancestor_mum;
    else
    {

        AliAODMCParticle *mcp_ancestor0 = (AliAODMCParticle *)mcarray->At(index_ancestor_mum); // Muon ancestor definition
        Int_t PDGcode_ancestor0 = mcp_ancestor0->GetPdgCode();
        if (Verbosity)
            printf("Round: %d Index Muon Ancestor %d | PDG Muon Ancestor %d | \n", round, index_ancestor_mum, PDGcode_ancestor0);

        while (testing)
        {
            index_ancestor_mum = mcp_ancestor0->GetMother();
            if (index_ancestor_mum < 2)
            {
                testing = kFALSE;
                printf("HF from Powheg\n");
                break;
            }
            mcp_ancestor0 = (AliAODMCParticle *)mcarray->At(index_ancestor_mum); // Salto gen Muone
            PDGcode_ancestor0 = mcp_ancestor0->GetPdgCode();

            if (Verbosity)
                printf("Round %d | Index Mum %d | PDG Muon Ancestor %d | \n", round, index_ancestor_mum, PDGcode_ancestor0);
            round++;
        }

        Muon_properties[1] = index_ancestor_mum;
    }

    return Muon_properties;
}

Int_t IsPrompt(AliAODMCParticle *mcp_mumum, TClonesArray *mcarray)
{
    Bool_t Verbosity = kFALSE;
    Int_t Final_PDG = mcp_mumum->GetPdgCode();
    Int_t ancestor_mum0 = mcp_mumum->GetMother();
    if (ancestor_mum0 == -1)
        return Final_PDG;
    AliAODMCParticle *mcp_ancestor0 = (AliAODMCParticle *)mcarray->At(ancestor_mum0); // Definizione Antenato Muone
    Int_t PDGcode_ancestor0 = mcp_ancestor0->GetPdgCode();
    if (Verbosity)
        printf("Index Muon Ancestor %d | PDG Muon Ancestor %d | \n", ancestor_mum0, PDGcode_ancestor0);

    Int_t round = 0;
    Bool_t Prompt = kTRUE;
    Bool_t testing = kTRUE;

    while (testing)
    {
        if (Verbosity)
            printf("Round %d | Index Muon Ancestor %d | PDG Muon Ancestor %d | \n", round, ancestor_mum0, PDGcode_ancestor0);
        round++;
        if (ancestor_mum0 < 2)
        {
            testing = kFALSE;
            break;
        }
        else if (TMath::Abs(PDGcode_ancestor0) == 5)
        {
            if (Verbosity)
                printf("Muon from Beauty Quark \n");
            Final_PDG = PDGcode_ancestor0;
            testing = kFALSE;
            Prompt = kFALSE;
            break;
        }
        else if (TMath::Abs(PDGcode_ancestor0) > 500 && TMath::Abs(PDGcode_ancestor0) < 600)
        {
            if (Verbosity)
                printf("Muon from Beauty Meson \n");
            Final_PDG = PDGcode_ancestor0;
            testing = kFALSE;
            Prompt = kFALSE;
            break;
        }
        else if (TMath::Abs(PDGcode_ancestor0) > 5000 && TMath::Abs(PDGcode_ancestor0) < 6000)
        {
            if (Verbosity)
                printf("Muon from Beauty Barion \n");
            Final_PDG = PDGcode_ancestor0;
            testing = kFALSE;
            Prompt = kFALSE;
            break;
        }
        else
        {
            ancestor_mum0 = mcp_ancestor0->GetMother();
            if (ancestor_mum0 > 0)
                mcp_ancestor0 = (AliAODMCParticle *)mcarray->At(ancestor_mum0); // Salto gen Muone
            else
            {
                mcp_ancestor0 = (AliAODMCParticle *)mcarray->At(mcp_mumum->GetMother());
                PDGcode_ancestor0 = mcp_ancestor0->GetPdgCode();
                testing = kFALSE;
            }
            PDGcode_ancestor0 = mcp_ancestor0->GetPdgCode();
        }
    }

    return Final_PDG;
}

void AliAnalysisTaskDimuon_HighMass::Terminate(Option_t *)
{
    printf("Terminate Task\n");
}
