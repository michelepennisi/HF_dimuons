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

/* $Id: AliAnalysisTaskDimuonPythia.cxx $ */

//-----------------------------------------------------------------------------
// Analysis task to create a tree containing MC embedding information
// for both generated and reconstructed particles
// R. Arnaldi
//-----------------------------------------------------------------------------

// ROOT includes
#include "TROOT.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "AliInputEventHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODMCParticle.h"
#include "AliMultSelection.h"
#include "AliMuonTrackCuts.h"

#include "AliAnalysisTaskDimuonPythia.h"

Int_t IsPrompt(AliAODMCParticle *mcp_mumum, TClonesArray *mcarray);

ClassImp(AliAnalysisTaskDimuonPythia)
    //__________________________________________________________________________
    AliAnalysisTaskDimuonPythia::AliAnalysisTaskDimuonPythia() : AliAnalysisTaskSE(),
                                                         fBeamEnergy(0.),
                                                         fkAnalysisType(0x0),
                                                         fPeriod(0x0),
                                                         fAODEvent(0x0),
                                                         fOutputTree(0x0),
                                                         fMuonTrackCuts(0x0),
                                                         fNMuons_gen(0x0),
                                                         fNDimu_gen(0x0),
                                                         fNMuons_rec(0x0),
                                                         fPercentV0M(0x0)

{
    /// Default ctor.
    printf("fMuons_dim %d\n", fMuons_dim);
    for (Int_t i = 0; i < fMuons_dim; i++)
    {

        fPDG_HFquark_rec[i] = 9999.;
        fPt_HFquark_rec[i] = 9999.;
        fY_HFquark_rec[i] = 9999.;

        fPDG_HFquark_gen[i] = 9999.;
        fPt_HFquark_gen[i] = 9999.;
        fY_HFquark_gen[i] = 9999.;

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
AliAnalysisTaskDimuonPythia::AliAnalysisTaskDimuonPythia(const char *name) : AliAnalysisTaskSE(name),
                                                                     fBeamEnergy(0.),
                                                                     fkAnalysisType(0x0),
                                                                     fPeriod(0x0),
                                                                     fAODEvent(0x0),
                                                                     fOutputTree(0x0),
                                                                     fMuonTrackCuts(0x0),
                                                                     fNMuons_gen(0x0),
                                                                     fNDimu_gen(0x0),
                                                                     fNMuons_rec(0x0),
                                                                     fPercentV0M(0x0)
{
    //
    // Constructor. Initialization of Inputs and Outputs
    //

    Info("AliAnalysisTaskDimuonPythia", "Calling Constructor");

    printf("fMuons_dim %d\n", fMuons_dim);
    for (Int_t i = 0; i < fMuons_dim; i++)
    {
        fPDG_HFquark_rec[i] = 9999.;
        fPt_HFquark_rec[i] = 9999.;
        fY_HFquark_rec[i] = 9999.;

        fPDG_HFquark_gen[i] = 9999.;
        fPt_HFquark_gen[i] = 9999.;
        fY_HFquark_gen[i] = 9999.;

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
AliAnalysisTaskDimuonPythia &AliAnalysisTaskDimuonPythia::operator=(const AliAnalysisTaskDimuonPythia &c)
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
AliAnalysisTaskDimuonPythia::AliAnalysisTaskDimuonPythia(const AliAnalysisTaskDimuonPythia &c) : AliAnalysisTaskSE(c),
                                                                                     fBeamEnergy(c.fBeamEnergy),
                                                                                     fkAnalysisType(c.fkAnalysisType),
                                                                                     fPeriod(c.fPeriod),
                                                                                     fAODEvent(c.fAODEvent),
                                                                                     fOutputTree(c.fOutputTree),
                                                                                     fMuonTrackCuts(c.fMuonTrackCuts),
                                                                                     fNMuons_gen(c.fNMuons_gen),
                                                                                     fNMuons_rec(c.fNMuons_rec),
                                                                                     fNDimu_rec(c.fNDimu_rec),
                                                                                     fPercentV0M(c.fPercentV0M)
{
    //
    // Copy Constructor
    //
}

//___________________________________________________________________________
AliAnalysisTaskDimuonPythia::~AliAnalysisTaskDimuonPythia()
{
    //
    // destructor
    //
    Info("~AliAnalysisTaskDimuonPythia", "Calling Destructor");
    if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis)
        delete fOutputTree;
}

//___________________________________________________________________________
void AliAnalysisTaskDimuonPythia::NotifyRun()
{
    fMuonTrackCuts->SetRun(fInputHandler);
}

//___________________________________________________________________________
void AliAnalysisTaskDimuonPythia::UserCreateOutputObjects()
{

    if (fOutputTree)
        return;

    OpenFile(1, "RECREATE");

    fOutputTree = new TTree("MCTree", "Data Tree");

    // fOutputTree->Branch("PercentV0M",&fPercentV0M,"PercentV0M/D");
    fOutputTree->Branch("N_HFquarks_gen", &fN_HFquarks_gen, "N_HFquarks_gen/I");

    fOutputTree->Branch("PDG_HFquark_gen", fPDG_HFquark_gen, "PDG_HFquark_gen[N_HFquarks_gen]/I");
    fOutputTree->Branch("Pt_HFquark_gen", fPt_HFquark_gen, "Pt_HFquark_gen[N_HFquarks_gen]/D");
    fOutputTree->Branch("Y_HFquark_gen", fY_HFquark_gen, "Y_HFquark_gen[N_HFquarks_gen]/D");

    // fOutputTree->Branch("N_HFquarks_rec", &fN_HFquarks_rec, "N_HFquarks_rec/I");

    // fOutputTree->Branch("PDG_HFquark_rec", fPDG_HFquark_rec, "PDG_HFquark_rec[N_HFquarks_rec]/I");
    // fOutputTree->Branch("Pt_HFquark_rec", fPt_HFquark_rec, "Pt_HFquark_rec[N_HFquarks_rec]/D");
    // fOutputTree->Branch("Y_HFquark_rec", fY_HFquark_rec, "Y_HFquark_rec[N_HFquarks_rec]/D");

    fOutputTree->Branch("NDimu_gen", &fNDimu_gen, "NDimu_gen/I");
    fOutputTree->Branch("DimuMu_gen", fDimuMu_gen, "DimuMu_gen[NDimu_gen][2]/I");
    fOutputTree->Branch("DimuPt_gen", fDimuPt_gen, "DimuPt_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuPx_gen", fDimuPx_gen, "DimuPx_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuPy_gen", fDimuPy_gen, "DimuPy_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuPz_gen", fDimuPz_gen, "DimuPz_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuY_gen", fDimuY_gen, "DimuY_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuMass_gen", fDimuMass_gen, "DimuMass_gen[NDimu_gen]/D");
    fOutputTree->Branch("DimuCharge_gen", fDimuCharge_gen, "DimuCharge_gen[NDimu_gen]/I");

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

    fOutputTree->ls();

    PostData(1, fOutputTree);
}

//_________________________________________________
void AliAnalysisTaskDimuonPythia::UserExec(Option_t *)
{

    // printf("Entro in UserExec\n" );
    fNMuons_gen = 0;
    fNDimu_gen = 0;
    fNMuons_rec = 0;
    fNDimu_rec = 0;
    fN_HFquarks_gen = 0;
    fN_HFquarks_rec = 0;
    for (Int_t i = 0; i < fMuons_dim; i++)
    {
        fPDG_HFquark_rec[i] = 9999.;
        fPt_HFquark_rec[i] = 9999.;
        fY_HFquark_rec[i] = 9999.;

        fPDG_HFquark_gen[i] = 9999.;
        fPt_HFquark_gen[i] = 9999.;
        fY_HFquark_gen[i] = 9999.;

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
    printf("%s\n",firedtrigger.Data());
    // if (firedtrigger.Contains("CINT7-B-NOPF-MUFAST"))
    // {
    //     // TriggerSelected_CINT7_CENT = kTRUE;
    //     printf("firedtrigger class %s\n", firedtrigger.Data());
    // }
    // AliMultSelection *MultSelection = (AliMultSelection * ) fAODEvent->FindListObject("MultSelection");
    // fPercentV0M = (Double_t) MultSelection->GetMultiplicityPercentile("V0M");

    Int_t nHFquark_gen = 0;

    Int_t nmu_gen = 0;
    Int_t ndimu_gen = 0;

    Int_t LabelOld1[10000] = {999};
    Int_t LabelOld2[10000] = {999};
    Bool_t GoodMuon[10000] = {kFALSE};

    for (Int_t i = 0; i < mcarray->GetEntries(); i++)
    {
        AliAODMCParticle *mcp0 = (AliAODMCParticle *)mcarray->At(i);
        Int_t PDG_mcp0 = mcp0->GetPdgCode();
        if (mcp0->Y() < -4.0 || mcp0->Y() > -2.5)
            continue;
        if (TMath::Abs(mcp0->Y()) > 15.0)
            continue;

        if (TMath::Abs(PDG_mcp0) == 4 || TMath::Abs(PDG_mcp0) == 5)
        {
            // printf("Array index %d || PDG %d || Array index FIGLIA1 %d || Array index FIGLIA2 %d\n", i, mcp0->GetPdgCode(), mcp0->GetDaughterFirst(), mcp0->GetDaughterLast());
            AliAODMCParticle *mcp0_daughter1;
            AliAODMCParticle *mcp0_daughter2;
            if (mcp0->GetDaughterFirst() != -1)
            {
                mcp0_daughter1 = (AliAODMCParticle *)mcarray->At(mcp0->GetDaughterFirst());
            }
            else
                mcp0_daughter1 = (AliAODMCParticle *)mcarray->At(i);

            if (mcp0->GetDaughterLast() != -1)
            {
                mcp0_daughter2 = (AliAODMCParticle *)mcarray->At(mcp0->GetDaughterLast());
            }
            else
                mcp0_daughter2 = (AliAODMCParticle *)mcarray->At(i);

            // printf("Array index %d || PDG %d || PT %f || Y %f || PDG FIGLIA1 %d || PDG FIGLIA2 %d\n", i, mcp0->GetPdgCode(), mcp0->Pt(), mcp0->Y(), mcp0_daughter1->GetPdgCode(), mcp0_daughter2->GetPdgCode());

            Int_t PDG_daughter1 = mcp0_daughter1->GetPdgCode();
            Int_t PDG_daughter2 = mcp0_daughter2->GetPdgCode();

            if ((TMath::Abs(PDG_daughter1) > 400 && TMath::Abs(PDG_daughter1) < 600) || (TMath::Abs(PDG_daughter1) > 4000 && TMath::Abs(PDG_daughter1) < 6000) || (TMath::Abs(PDG_daughter2) > 400 && TMath::Abs(PDG_daughter2) < 600) || (TMath::Abs(PDG_daughter2) > 4000 && TMath::Abs(PDG_daughter2) < 6000))
            {
                // printf("Array index %d || PDG %d || Array index FIGLIA1 %d || Array index FIGLIA2 %d\n", i, mcp0->GetPdgCode(), mcp0->GetDaughterFirst(), mcp0->GetDaughterLast());

                // printf("Array index %d || PDG %d || PT %f || Y %f || PDG FIGLIA1 %d || PDG FIGLIA2 %d\n", i, mcp0->GetPdgCode(), mcp0->Pt(), mcp0->Y(), mcp0_daughter1->GetPdgCode(), mcp0_daughter2->GetPdgCode());

                fPDG_HFquark_gen[nHFquark_gen] = mcp0->GetPdgCode();
                fPt_HFquark_gen[nHFquark_gen] = mcp0->Pt();
                fY_HFquark_gen[nHFquark_gen] = mcp0->Y();
                nHFquark_gen++;
            }
        }
        if (TMath::Abs(PDG_mcp0) != 13)
            continue;

        TLorentzVector vector_mcp0;
        //        printf("Find a Muon %d \n",PDG_mcp0);
        fY_gen[nmu_gen] = mcp0->Y();
        fPt_gen[nmu_gen] = mcp0->Pt();
        fE_gen[nmu_gen] = mcp0->E();
        fPx_gen[nmu_gen] = mcp0->Px();
        fPy_gen[nmu_gen] = mcp0->Py();
        fPz_gen[nmu_gen] = mcp0->Pz();
        fEta_gen[nmu_gen] = mcp0->Eta();
        fPhi_gen[nmu_gen] = mcp0->Phi();
        // printf("Charge a mano: %d || Charge metodo: %d",-TMath::Sign(1, PDG_mcp0),mcp0->Charge());
        fCharge_gen[nmu_gen] = -TMath::Sign(1, PDG_mcp0); // Sign minus necessary
                                                          //        printf("PDG_mcp0 %d, TMath::Sign(-1,PDG_mcp0) %d \n",PDG_mcp0,fCharge_gen[nmu_gen]);
        fTheta_gen[nmu_gen] = mcp0->Theta();

        mcp0->Momentum(vector_mcp0);

        //        printf("nmu_gen %d) MCParticle Pt %0.1f, Px %0.1f, Py %0.1f, Pz %0.1f, Y %0.1f\n",nmu_gen,fPt_gen[nmu_gen],fPx_gen[nmu_gen],fPy_gen[nmu_gen],fPz_gen[nmu_gen],fY_gen[nmu_gen]);
        //        printf("Pippo Pt %0.1f, Px %0.1f, Py %0.1f, Pz %0.1f, Y %0.1f\n",vector_mcp0.Pt(),vector_mcp0.Px(),vector_mcp0.Py(),vector_mcp0.Pz(),vector_mcp0.Y());

        Int_t mum_mcp0 = mcp0->GetMother();
        AliAODMCParticle *mcp0_mother = (AliAODMCParticle *)mcarray->At(mum_mcp0);
        Int_t PDG_mum = mcp0_mother->GetPdgCode();

        //        if (mcp0->Y()>-2.0 && mcp0->Y()<2.0) {
        //            printf("--------Muon-------\n");
        //            mcp0->Print();
        //            printf("PDGcode mother: %d \n",PDG_mum);
        //            printf("-------Mother------\n");
        //            mcp0_mother->Print();
        //            AliAODMCParticle *mcp0_grandmother = (AliAODMCParticle *)mcarray->At(mcp0_mother->GetMother());
        //            printf("----GrandMother----\n");
        //            mcp0_grandmother->Print();
        //            Int_t PDG_grandmum=mcp0_grandmother->GetPdgCode();
        //            printf("PDGcode grandmother: %d \n",PDG_grandmum);
        //        }

        if (TMath::Abs(PDG_mum) == 4)
        {
            //            printf ("Muon from Charm Quark \n");
            fPDGmum_gen[nmu_gen] = IsPrompt(mcp0_mother, mcarray);
        }
        else if (TMath::Abs(PDG_mum) > 400 && TMath::Abs(PDG_mum) < 500)
        {
            //            printf ("Muon from Charm Meson \n");
            fPDGmum_gen[nmu_gen] = IsPrompt(mcp0_mother, mcarray);
        }
        else if (TMath::Abs(PDG_mum) > 4000 && TMath::Abs(PDG_mum) < 5000)
        {
            //            printf ("Muon from Charm Barion \n");
            fPDGmum_gen[nmu_gen] = IsPrompt(mcp0_mother, mcarray);
        }
        else
            fPDGmum_gen[nmu_gen] = PDG_mum;

        nmu_gen++;
        for (Int_t j = i + 1; j < mcarray->GetEntries(); j++)
        {
            AliAODMCParticle *mcp1 = (AliAODMCParticle *)mcarray->At(j);
            Int_t PDG_mcp1 = mcp1->GetPdgCode();
            if (TMath::Abs(PDG_mcp1) != 13)
                continue;
            if (TMath::Abs(mcp1->Y()) > 15.0)
                continue;
            if (mcp1->Y() < -4.0 || mcp1->Y() > -2.5)
                continue;
            //            if(mcp1->IsSecondaryFromMaterial()) continue;
            TLorentzVector vector_mcp1;
            mcp1->Momentum(vector_mcp1);

            //            printf("j %d) MCParticle Pt %0.1f, Px %0.1f, Py %0.1f, Pz %0.1f, Y %0.1f\n",j,mcp1->Pt(),mcp1->Px(),mcp1->Py(),mcp1->Pz(),mcp0->Y());
            //            printf("vector_mcp1, %0.1f, Px %0.1f, Py %0.1f, Pz %0.1f, Y %0.1f\n",vector_mcp1.Pt(),vector_mcp1.Px(),vector_mcp1.Py(),vector_mcp1.Pz(),vector_mcp1.Y());

            TLorentzVector dimu = vector_mcp0 + vector_mcp1;

            fDimuPt_gen[ndimu_gen] = dimu.Pt();
            fDimuPx_gen[ndimu_gen] = dimu.Px();
            fDimuPy_gen[ndimu_gen] = dimu.Py();
            fDimuPz_gen[ndimu_gen] = dimu.Pz();
            fDimuY_gen[ndimu_gen] = dimu.Rapidity();
            fDimuMass_gen[ndimu_gen] = dimu.M();
            fDimuCharge_gen[ndimu_gen] = -(TMath::Sign(1, PDG_mcp0) + TMath::Sign(1, PDG_mcp1));
            //            printf("PDG_mcp0 %d | TMath::Sign(-1,PDG_mcp0) %d | PDG_mcp1 %d | TMath::Sign(-1,PDG_mcp0) %d | fDimuCharge_gen[ndimu_gen] %d\n",PDG_mcp0,-TMath::Sign(1,PDG_mcp0),PDG_mcp0,-TMath::Sign(1,PDG_mcp1),fDimuCharge_gen[ndimu_gen]);
            LabelOld1[ndimu_gen] = i;
            LabelOld2[ndimu_gen] = j;
            GoodMuon[i] = kTRUE;
            GoodMuon[j] = kTRUE;
            ndimu_gen++;
        }
    }

    printf("NUMERO DI HF quark per evento %d\n", nHFquark_gen);
    fN_HFquarks_gen = nHFquark_gen;
    fNMuons_gen = nmu_gen;
    fNDimu_gen = ndimu_gen;

    for (Int_t i = 0; i < ndimu_gen; i++)
    {
        Int_t LabelNew1 = 0;
        Int_t LabelNew2 = 0;
        for (int j = 0; j < LabelOld1[i]; j++)
        {
            if (GoodMuon[j])
                LabelNew1++;
        }
        for (int j = 0; j < LabelOld2[i]; j++)
        {
            if (GoodMuon[j])
                LabelNew2++;
        }
        fDimuMu_gen[i][0] = LabelNew1;
        fDimuMu_gen[i][1] = LabelNew2;
    }

    Int_t numtracks = fAODEvent->GetNumberOfTracks();

    Int_t nHFquark_rec = 0;

    Int_t nmu_rec = 0;
    Int_t ndimu_rec = 0;

    Int_t LabelOld1_rec[10000] = {999};
    Int_t LabelOld2_rec[10000] = {999};
    Int_t GoodMuon_rec[10000] = {kFALSE};

    for (Int_t j = 0; j < numtracks; j++)
    {
        AliAODTrack *mu0 = (AliAODTrack *)fAODEvent->GetTrack(j);
        if (mu0->GetLabel() == -1)
            continue;
        AliAODMCParticle *mctrack0 = (AliAODMCParticle *)mcarray->At(mu0->GetLabel());
        Int_t PDG_mctrack0 = mctrack0->GetPdgCode();
        // if (TMath::Abs(PDG_mctrack0) == 4 || TMath::Abs(PDG_mctrack0) == 5)
        // {
        //     if (true)
        //     {
        //         fPDG_HFquark_rec[nHFquark_rec] = PDG_mctrack0;
        //         fPt_HFquark_rec[nHFquark_rec] = mu0->Pt();
        //         fY_HFquark_rec[nHFquark_rec] = mu0->Y();
        //         nHFquark_rec++;
        //     }
        // }

        if (TMath::Abs(PDG_mctrack0) != 13)
            continue;
        //        if(mctrack0->IsSecondaryFromMaterial()) continue;

        if (mu0->Y() > -4.0 && mu0->Y() < -2.5)
            GoodMuon_rec[j] = kTRUE;

        if (!GoodMuon_rec[j])
            continue;

        fPt_rec[nmu_rec] = mu0->Pt();
        fPx_rec[nmu_rec] = mu0->Px();
        fPy_rec[nmu_rec] = mu0->Py();
        fPz_rec[nmu_rec] = mu0->Pz();
        fY_rec[nmu_rec] = mu0->Y();
        fEta_rec[nmu_rec] = mu0->Eta();
        fE_rec[nmu_rec] = mu0->E();
        fMatchTrig_rec[nmu_rec] = mu0->GetMatchTrigger();
        fMatchTrigChi2_rec[nmu_rec] = mu0->GetChi2MatchTrigger();
        fRAtAbsEnd_rec[nmu_rec] = mu0->GetRAtAbsorberEnd();
        if (fMuonTrackCuts->IsSelected(mu0))
            fpDCA_rec[nmu_rec] = 1;
        fPhi_rec[nmu_rec] = mu0->Phi();
        fTheta_rec[nmu_rec] = mu0->Theta();
        fCharge_rec[nmu_rec] = mu0->Charge();

        Int_t mum_mctrack0 = mctrack0->GetMother();
        AliAODMCParticle *mctrack0_mother = (AliAODMCParticle *)mcarray->At(mum_mctrack0);
        Int_t PDG_mum_mctrack0 = mctrack0_mother->GetPdgCode();

        if (TMath::Abs(PDG_mum_mctrack0) == 4)
        {
            //            printf ("Muon from Charm Quark \n");
            fPDGmum_rec[nmu_rec] = IsPrompt(mctrack0_mother, mcarray);
        }
        else if (TMath::Abs(PDG_mum_mctrack0) > 400 && TMath::Abs(PDG_mum_mctrack0) < 500)
        {
            //            printf ("Muon from Charm Meson \n");
            fPDGmum_rec[nmu_rec] = IsPrompt(mctrack0_mother, mcarray);
        }
        else if (TMath::Abs(PDG_mum_mctrack0) > 4000 && TMath::Abs(PDG_mum_mctrack0) < 5000)
        {
            //            printf ("Muon from Charm Barion \n");
            fPDGmum_rec[nmu_rec] = IsPrompt(mctrack0_mother, mcarray);
        }
        else
        {
            fPDGmum_rec[nmu_rec] = PDG_mum_mctrack0;
        }

        nmu_rec++;

        for (Int_t k = j + 1; k < numtracks; k++)
        {
            AliAODTrack *mu1 = (AliAODTrack *)fAODEvent->GetTrack(k);
            if (mu1->GetLabel() == -1)
                continue;
            AliAODMCParticle *mctrack1 = (AliAODMCParticle *)mcarray->At(mu1->GetLabel());
            Int_t PDG_mctrack1 = mctrack1->GetPdgCode();
            if (TMath::Abs(PDG_mctrack1) != 13)
                continue;
            //            if(mctrack1->IsSecondaryFromMaterial()) continue;

            if (mu1->Y() > -4.0 && mu1->Y() < -2.5)
                GoodMuon_rec[k] = kTRUE;

            if (!GoodMuon_rec[k])
                continue;

            AliAODDimuon *dimu = new AliAODDimuon(mu0, mu1);
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
            if (mu0->GetMatchTrigger() > 1 || mu1->GetMatchTrigger() > 1)
                fDimuMatch_rec[ndimu_rec] = 1;
            if (mu0->GetMatchTrigger() > 1 && mu1->GetMatchTrigger() > 1)
                fDimuMatch_rec[ndimu_rec] = 2;

            LabelOld1_rec[ndimu_rec] = j;
            LabelOld2_rec[ndimu_rec] = k;

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
        for (int j = 0; j < LabelOld1_rec[i]; j++)
        {
            if (GoodMuon_rec[j])
                LabelNew1++;
        }
        for (int j = 0; j < LabelOld2_rec[i]; j++)
        {
            if (GoodMuon_rec[j])
                LabelNew2++;
        }
        fDimuMu_rec[i][0] = LabelNew1;
        fDimuMu_rec[i][1] = LabelNew2;
    }

    fOutputTree->Fill();
    PostData(1, fOutputTree);
}

//________________________________________________________________________

Int_t IsPrompt(AliAODMCParticle *mcp_mumum, TClonesArray *mcarray)
{

    Int_t Final_PDG = mcp_mumum->GetPdgCode();
    //    printf("Canditate Muon from Charm Particle \n");
    //    printf("Index Muon Mum %d | PDG Muon Mum %d | \n",mcp_mu_charm_canditate->GetPdgCode());
    //    mcp_mumum->Print();

    Int_t ancestor_mum0 = mcp_mumum->GetMother();
    if (ancestor_mum0 == -1)
        return Final_PDG;
    AliAODMCParticle *mcp_ancestor0 = (AliAODMCParticle *)mcarray->At(ancestor_mum0); // Definizione Antenato Muone
    Int_t PDGcode_ancestor0 = mcp_ancestor0->GetPdgCode();
    //    printf("Index Muon Ancestor %d | PDG Muon Ancestor %d | \n",ancestor_mum0, PDGcode_ancestor0);

    Int_t round = 0;
    Bool_t Prompt = kTRUE;
    Bool_t testing = kTRUE;

    while (testing)
    {
        //        printf("Round %d | Index Muon Ancestor %d | PDG Muon Ancestor %d | \n",round,ancestor_mum0, PDGcode_ancestor0);
        round++;
        if (ancestor_mum0 < 2)
        {
            testing = kFALSE;
            break;
        }
        else if (TMath::Abs(PDGcode_ancestor0) == 5)
        {
            //            printf("Muon from Beauty Quark \n");
            Final_PDG = PDGcode_ancestor0;
            testing = kFALSE;
            Prompt = kFALSE;
            break;
        }
        else if (TMath::Abs(PDGcode_ancestor0) > 500 && TMath::Abs(PDGcode_ancestor0) < 600)
        {
            //            printf("Muon from Beauty Meson \n");
            Final_PDG = PDGcode_ancestor0;
            testing = kFALSE;
            Prompt = kFALSE;
            break;
        }
        else if (TMath::Abs(PDGcode_ancestor0) > 5000 && TMath::Abs(PDGcode_ancestor0) < 6000)
        {
            //            printf("Muon from Beauty Barion \n");
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

void AliAnalysisTaskDimuonPythia::Terminate(Option_t *)
{
    printf("Terminate Task\n");
}
