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

/* $Id: AliAnalysisTaskDimuonPowheg.cxx $ */

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

#include "AliAnalysisTaskDimuonPowheg.h"

Int_t IsPrompt(AliAODMCParticle *mcp_mumum, TClonesArray *mcarray);

ClassImp(AliAnalysisTaskDimuonPowheg)
    //__________________________________________________________________________
    AliAnalysisTaskDimuonPowheg::AliAnalysisTaskDimuonPowheg() : AliAnalysisTaskSE(),
                                                                 fBeamEnergy(0.),
                                                                 fkAnalysisType(0x0),
                                                                 fPeriod(0x0),
                                                                 fAODEvent(0x0),
                                                                 fOutputTree(0x0),
                                                                 fMuonTrackCuts(0x0),
                                                                 fNMuons_gen(0x0),
                                                                 fNHadrons_gen(0x0),
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
AliAnalysisTaskDimuonPowheg::AliAnalysisTaskDimuonPowheg(const char *name) : AliAnalysisTaskSE(name),
                                                                             fBeamEnergy(0.),
                                                                             fkAnalysisType(0x0),
                                                                             fPeriod(0x0),
                                                                             fAODEvent(0x0),
                                                                             fOutputTree(0x0),
                                                                             fMuonTrackCuts(0x0),
                                                                             fNMuons_gen(0x0),
                                                                             fNHadrons_gen(0x0),
                                                                             fNDimu_gen(0x0),
                                                                             fNMuons_rec(0x0),
                                                                             fPercentV0M(0x0)
{
    //
    // Constructor. Initialization of Inputs and Outputs
    //

    Info("AliAnalysisTaskDimuonPowheg", "Calling Constructor");

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
AliAnalysisTaskDimuonPowheg &AliAnalysisTaskDimuonPowheg::operator=(const AliAnalysisTaskDimuonPowheg &c)
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
AliAnalysisTaskDimuonPowheg::AliAnalysisTaskDimuonPowheg(const AliAnalysisTaskDimuonPowheg &c) : AliAnalysisTaskSE(c),
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
                                                                                                 fPercentV0M(c.fPercentV0M)
{
    //
    // Copy Constructor
    //
}

//___________________________________________________________________________
AliAnalysisTaskDimuonPowheg::~AliAnalysisTaskDimuonPowheg()
{
    //
    // destructor
    //
    Info("~AliAnalysisTaskDimuonPowheg", "Calling Destructor");
    if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis)
        delete fOutputTree;
}

//___________________________________________________________________________
void AliAnalysisTaskDimuonPowheg::NotifyRun()
{
    fMuonTrackCuts->SetRun(fInputHandler);
}

//___________________________________________________________________________
void AliAnalysisTaskDimuonPowheg::UserCreateOutputObjects()
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
void AliAnalysisTaskDimuonPowheg::UserExec(Option_t *)
{

    // printf("Entro in UserExec\n" );
    fNMuons_gen = 0;
    fNHadrons_gen = 0;
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
    printf("BEAUTY SIM \n");
    AliAODHeader *aodheader = dynamic_cast<AliAODHeader *>(fAODEvent->GetHeader());
    TString firedtrigger = aodheader->GetFiredTriggerClasses();
    printf("%s\n", firedtrigger.Data());

    Int_t nHFquark_gen = 0;
    Int_t nmu_gen = 0;
    Int_t nHadron_gen = 0;
    Int_t ndimu_gen = 0;

    Int_t LabelOld1[10000] = {999};
    Int_t LabelOld2[10000] = {999};
    Bool_t GoodMuon[10000] = {kFALSE};

    printf("Numero part %d \n", mcarray->GetEntries());
    for (Int_t i_nMCpart = 0; i_nMCpart < mcarray->GetEntries(); i_nMCpart++)
    {
        AliAODMCParticle *MC_part0 = (AliAODMCParticle *)mcarray->At(i_nMCpart);
        Int_t PDG_MC_part0 = MC_part0->GetPdgCode();
        if (MC_part0->Y() < -15.0 || MC_part0->Y() > 15.0)
            continue;
        Bool_t Muon = kFALSE;
        Bool_t Charm_Hadron = kFALSE;
        Bool_t Beauty_Hadron = kFALSE;

        if (TMath::Abs(PDG_MC_part0) == 13)
            Muon = kTRUE;
        else if ((TMath::Abs(PDG_MC_part0) > 400 && TMath::Abs(PDG_MC_part0) < 500) || (TMath::Abs(PDG_MC_part0) > 4000 && TMath::Abs(PDG_MC_part0) < 5000))
            Charm_Hadron = kTRUE;
        else if ((TMath::Abs(PDG_MC_part0) > 500 && TMath::Abs(PDG_MC_part0) < 600) || (TMath::Abs(PDG_MC_part0) > 5000 && TMath::Abs(PDG_MC_part0) < 6000))
            Beauty_Hadron = kTRUE;

        if (!Muon && !Charm_Hadron && !Beauty_Hadron)
            continue;

        Int_t index_Mum_MC_part0 = MC_part0->GetMother();
        AliAODMCParticle *Mum_MC_part0 = (AliAODMCParticle *)mcarray->At(index_Mum_MC_part0);

        Int_t PDG_Mum_MC_part0 = Mum_MC_part0->GetPdgCode();
        if (TMath::Abs(PDG_Mum_MC_part0) == TMath::Abs(PDG_MC_part0)) // Remove particles formally produced by the same particle (PYTHIA transport)
            continue;
        TLorentzVector vector_mcp0;

        if (Muon)
        {
            if ((PDG_Mum_MC_part0 != 23) && !(TMath::Abs(PDG_Mum_MC_part0) > 400 && TMath::Abs(PDG_Mum_MC_part0) < 500) && !(TMath::Abs(PDG_Mum_MC_part0) > 4000 && TMath::Abs(PDG_Mum_MC_part0) < 5000) && !(TMath::Abs(PDG_Mum_MC_part0) > 500 && TMath::Abs(PDG_Mum_MC_part0) < 600) && !(TMath::Abs(PDG_Mum_MC_part0) > 5000 && TMath::Abs(PDG_Mum_MC_part0) < 6000))
                continue;
            MC_part0->Print();
            printf("i_nMCpart %d) PDG: %d | MCParticle Pt %0.1f | Px %0.1f | Py %0.1f | Pz %0.1f | Y %0.1f || Mother %d , PDG Mother %d \n", i_nMCpart, PDG_MC_part0, MC_part0->Pt(), MC_part0->Px(), MC_part0->Py(), MC_part0->Pz(), MC_part0->Y(), index_Mum_MC_part0, PDG_Mum_MC_part0);
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

            Int_t final_mom_mu0 = IsPrompt(Mum_MC_part0, mcarray);
            printf("Final Muon mother: %d\n", final_mom_mu0);
            fPDGmum_gen[nmu_gen] = final_mom_mu0;
            nmu_gen++;
        }
        else if (Charm_Hadron || Beauty_Hadron)
        {
            MC_part0->Print();
            printf("i_nMCpart %d) PDG: %d | MCParticle Pt %0.1f | Px %0.1f | Py %0.1f | Pz %0.1f | Y %0.1f || Mother %d , PDG Mother %d \n", i_nMCpart, PDG_MC_part0, MC_part0->Pt(), MC_part0->Px(), MC_part0->Py(), MC_part0->Pz(), MC_part0->Y(), index_Mum_MC_part0, PDG_Mum_MC_part0);
            Int_t final_mom_mu0 = IsPrompt(Mum_MC_part0, mcarray);
            printf("Final Hadron mother: %d\n", final_mom_mu0);
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

        for (Int_t j_nMCpart1 = i_nMCpart + 1; j_nMCpart1 < mcarray->GetEntries(); j_nMCpart1++)
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
            if (TMath::Abs(PDG_MC_part1) == TMath::Abs(PDG_Mum_MC_part1)) // Remove particles formally produced by the same particle (PYTHIA transport)
                continue;

            Int_t final_Mum_MC_part1 = PDG_Mum_MC_part1;

            if ((PDG_Mum_MC_part1 != 23) && !(TMath::Abs(PDG_Mum_MC_part1) > 400 && TMath::Abs(PDG_Mum_MC_part1) < 500) && !(TMath::Abs(PDG_Mum_MC_part1) > 4000 && TMath::Abs(PDG_Mum_MC_part1) < 5000))
                continue;
            printf("Find A DIMUON");
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

    printf("NUMERO DI HF quark per evento %d\n", nHFquark_gen);
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

    Int_t nHFquark_rec = 0;

    Int_t nmu_rec = 0;
    Int_t ndimu_rec = 0;

    Int_t LabelOld1_rec[10000] = {999};
    Int_t LabelOld2_rec[10000] = {999};
    Int_t GoodMuon_rec[10000] = {kFALSE};

    for (Int_t i_Track0 = 0; i_Track0 < numtracks; i_Track0++)
    {
        AliAODTrack *Track0 = (AliAODTrack *)fAODEvent->GetTrack(i_Track0);
        if (Track0->GetLabel() == -1)
            continue;
        //        if(mctrack0->IsSecondaryFromMaterial()) continue;

        if (Track0->Y() > -15.0 && Track0->Y() < 15.0)
            GoodMuon_rec[i_Track0] = kTRUE;

        if (!GoodMuon_rec[i_Track0])
            continue;

        AliAODMCParticle *MCpart_Track0 = (AliAODMCParticle *)mcarray->At(Track0->GetLabel());
        Int_t PDG_MCpart_Track0 = MCpart_Track0->GetPdgCode();

        if (TMath::Abs(PDG_MCpart_Track0) != 13)
            continue;

        Int_t index_Mum_MCpart_Track0 = MCpart_Track0->GetMother();
        AliAODMCParticle *Mum_MCpart_Track0 = (AliAODMCParticle *)mcarray->At(index_Mum_MCpart_Track0);
        Int_t PDG_Mum_MCpart_Track0 = Mum_MCpart_Track0->GetPdgCode();

        if ((PDG_Mum_MCpart_Track0 != 23) && !(TMath::Abs(PDG_Mum_MCpart_Track0) > 400 && TMath::Abs(PDG_Mum_MCpart_Track0) < 500) && !(TMath::Abs(PDG_Mum_MCpart_Track0) > 4000 && TMath::Abs(PDG_Mum_MCpart_Track0) < 5000) && !(TMath::Abs(PDG_Mum_MCpart_Track0) > 500 && TMath::Abs(PDG_Mum_MCpart_Track0) < 600) && !(TMath::Abs(PDG_Mum_MCpart_Track0) > 5000 && TMath::Abs(PDG_Mum_MCpart_Track0) < 6000))
            continue;

        MCpart_Track0->Print();
        
        fPt_rec[nmu_rec] = Track0->Pt();
        fPx_rec[nmu_rec] = Track0->Px();
        fPy_rec[nmu_rec] = Track0->Py();
        fPz_rec[nmu_rec] = Track0->Pz();
        fY_rec[nmu_rec] = Track0->Y();
        fEta_rec[nmu_rec] = Track0->Eta();
        fE_rec[nmu_rec] = Track0->E();
        fMatchTrig_rec[nmu_rec] = Track0->GetMatchTrigger();
        fMatchTrigChi2_rec[nmu_rec] = Track0->GetChi2MatchTrigger();
        fRAtAbsEnd_rec[nmu_rec] = Track0->GetRAtAbsorberEnd();
        if (fMuonTrackCuts->IsSelected(Track0))
            fpDCA_rec[nmu_rec] = 1;
        fPhi_rec[nmu_rec] = Track0->Phi();
        fTheta_rec[nmu_rec] = Track0->Theta();
        fCharge_rec[nmu_rec] = Track0->Charge();
        fPDGmum_rec[nmu_rec] = IsPrompt(Mum_MCpart_Track0, mcarray);;
        nmu_rec++;

        for (Int_t j_Track1 = i_Track0 + 1; j_Track1 < numtracks; j_Track1++)
        {
            AliAODTrack *Track1 = (AliAODTrack *)fAODEvent->GetTrack(j_Track1);
            if (Track1->GetLabel() == -1)
                continue;
            AliAODMCParticle *MCpart_Track1 = (AliAODMCParticle *)mcarray->At(Track1->GetLabel());
            Int_t PDG_mctrack1 = MCpart_Track1->GetPdgCode();
            if (TMath::Abs(PDG_mctrack1) != 13)
                continue;
            if (Track1->Y() > -15.0 && Track1->Y() < 15.0)
                GoodMuon_rec[j_Track1] = kTRUE;

            if (!GoodMuon_rec[j_Track1])
                continue;
            Int_t index_Mum_MCpart_Track1 = MCpart_Track1->GetMother();
            AliAODMCParticle *Mum_MCpart_Track1 = (AliAODMCParticle *)mcarray->At(index_Mum_MCpart_Track1);
            Int_t PDG_Mum_MCpart_Track1 = Mum_MCpart_Track1->GetPdgCode();

            if ((PDG_Mum_MCpart_Track1 != 23) && !(TMath::Abs(PDG_Mum_MCpart_Track1) > 400 && TMath::Abs(PDG_Mum_MCpart_Track1) < 500) && !(TMath::Abs(PDG_Mum_MCpart_Track1) > 4000 && TMath::Abs(PDG_Mum_MCpart_Track1) < 5000))
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
    Bool_t Verbosity = kTRUE;
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

void AliAnalysisTaskDimuonPowheg::Terminate(Option_t *)
{
    printf("Terminate Task\n");
}
