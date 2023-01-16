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

/* $Id: AliAnalysisTaskTriggerEff.cxx $ */

//-----------------------------------------------------------------------------
// Analysis task to study trigger efficiency response function of single muons
// M. Pennisi
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

#include "AliAnalysisTaskTriggerEff.h"

ClassImp(AliAnalysisTaskTriggerEff)
    //__________________________________________________________________________
    AliAnalysisTaskTriggerEff::AliAnalysisTaskTriggerEff() : AliAnalysisTaskSE(),
                                                             fBeamEnergy(0.),
                                                             fkAnalysisType(0x0),
                                                             fPeriod(0x0),
                                                             fAODEvent(0x0),
                                                             fMuonTrackCuts(0x0),
                                                             fPercentV0M(0x0),
                                                             fOutputList{0},
                                                             fHistPt_allmuon_alleta{0},
                                                             fHistPt_allmuon_eta1_allpt{0},
                                                             fHistPt_allmuon_eta2_allpt{0},
                                                             fHistPt_allmuon_eta3_allpt{0},
                                                             fHistPt_allmuon_alleta_lowpt_thr{0},
                                                             fHistPt_allmuon_eta1_lowpt_thr{0},
                                                             fHistPt_allmuon_eta2_lowpt_thr{0},
                                                             fHistPt_allmuon_eta3_lowpt_thr{0}

{
    fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
    fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
    fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
}

//__________________________________________________________________________
AliAnalysisTaskTriggerEff::AliAnalysisTaskTriggerEff(const char *name) : AliAnalysisTaskSE(name),
                                                                         fBeamEnergy(0.),
                                                                         fkAnalysisType(0x0),
                                                                         fPeriod(0x0),
                                                                         fAODEvent(0x0),
                                                                         fMuonTrackCuts(0x0),
                                                                         fPercentV0M(0x0),
                                                                         fOutputList{0},
                                                                         fHistPt_allmuon_alleta{0},
                                                                         fHistPt_allmuon_eta1_allpt{0},
                                                                         fHistPt_allmuon_eta2_allpt{0},
                                                                         fHistPt_allmuon_eta3_allpt{0},
                                                                         fHistPt_allmuon_alleta_lowpt_thr{0},
                                                                         fHistPt_allmuon_eta1_lowpt_thr{0},
                                                                         fHistPt_allmuon_eta2_lowpt_thr{0},
                                                                         fHistPt_allmuon_eta3_lowpt_thr{0}
{
    //
    // Constructor. Initialization of Inputs and Outputs
    //
    fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
    fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
    fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
    DefineOutput(1, TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskTriggerEff &AliAnalysisTaskTriggerEff::operator=(const AliAnalysisTaskTriggerEff &c)
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
AliAnalysisTaskTriggerEff::AliAnalysisTaskTriggerEff(const AliAnalysisTaskTriggerEff &c) : AliAnalysisTaskSE(c),
                                                                                           fBeamEnergy(c.fBeamEnergy),
                                                                                           fkAnalysisType(c.fkAnalysisType),
                                                                                           fPeriod(c.fPeriod),
                                                                                           fAODEvent(c.fAODEvent),
                                                                                           fMuonTrackCuts(c.fMuonTrackCuts),
                                                                                           fPercentV0M(c.fPercentV0M),
                                                                                           fOutputList(c.fOutputList),
                                                                                           fHistPt_allmuon_alleta(c.fHistPt_allmuon_alleta),
                                                                                           fHistPt_allmuon_eta1_allpt{c.fHistPt_allmuon_eta1_allpt},
                                                                                           fHistPt_allmuon_eta2_allpt{c.fHistPt_allmuon_eta2_allpt},
                                                                                           fHistPt_allmuon_eta3_allpt{c.fHistPt_allmuon_eta3_allpt},
                                                                                           fHistPt_allmuon_alleta_lowpt_thr{c.fHistPt_allmuon_alleta_lowpt_thr},
                                                                                           fHistPt_allmuon_eta1_lowpt_thr{c.fHistPt_allmuon_eta1_lowpt_thr},
                                                                                           fHistPt_allmuon_eta2_lowpt_thr{c.fHistPt_allmuon_eta2_lowpt_thr},
                                                                                           fHistPt_allmuon_eta3_lowpt_thr{c.fHistPt_allmuon_eta3_lowpt_thr}
{
    //
    // Copy Constructor
    //
}

//___________________________________________________________________________
AliAnalysisTaskTriggerEff::~AliAnalysisTaskTriggerEff()
{
    //
    // destructor
    //
    Info("~AliAnalysisTaskTriggerEff", "Calling Destructor");
    if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis)
        delete fOutputList;
}

//___________________________________________________________________________
void AliAnalysisTaskTriggerEff::NotifyRun()
{
    fMuonTrackCuts->SetRun(fInputHandler);
}

//___________________________________________________________________________
void AliAnalysisTaskTriggerEff::UserCreateOutputObjects()
{

    if (fOutputList)
        return;
    OpenFile(1, "RECREATE");

    fOutputList = new TList();
    fOutputList->SetOwner(true);

    fHistPt_allmuon_alleta = new TH1D("fHistPt_allmuon_alleta", "fHistPt_allmuon_alleta", 3000, 0, 30);
    fHistPt_allmuon_eta1_allpt = new TH1D("fHistPt_allmuon_eta1_allpt", "fHistPt_allmuon_eta1_allpt", 3000, 0, 30);
    fHistPt_allmuon_eta2_allpt = new TH1D("fHistPt_allmuon_eta2_allpt", "fHistPt_allmuon_eta2_allpt", 3000, 0, 30);
    fHistPt_allmuon_eta3_allpt = new TH1D("fHistPt_allmuon_eta3_allpt", "fHistPt_allmuon_eta3_allpt", 3000, 0, 30);

    fHistPt_allmuon_alleta_lowpt_thr = new TH1D("fHistPt_allmuon_alleta_lowpt_thr", "fHistPt_allmuon_alleta_lowpt_thr", 3000, 0, 30);
    fHistPt_allmuon_eta1_lowpt_thr = new TH1D("fHistPt_allmuon_eta1_lowpt_thr", "fHistPt_allmuon_eta1_lowpt_thr", 3000, 0, 30);
    fHistPt_allmuon_eta2_lowpt_thr = new TH1D("fHistPt_allmuon_eta2_lowpt_thr", "fHistPt_allmuon_eta2_lowpt_thr", 3000, 0, 30);
    fHistPt_allmuon_eta3_lowpt_thr = new TH1D("fHistPt_allmuon_eta3_lowpt_thr", "fHistPt_allmuon_eta3_lowpt_thr", 3000, 0, 30);

    fOutputList->Add(fHistPt_allmuon_alleta);
    fOutputList->Add(fHistPt_allmuon_eta1_allpt);
    fOutputList->Add(fHistPt_allmuon_eta2_allpt);
    fOutputList->Add(fHistPt_allmuon_eta3_allpt);

    fOutputList->Add(fHistPt_allmuon_alleta_lowpt_thr);
    fOutputList->Add(fHistPt_allmuon_eta1_lowpt_thr);
    fOutputList->Add(fHistPt_allmuon_eta2_lowpt_thr);
    fOutputList->Add(fHistPt_allmuon_eta3_lowpt_thr);

    fOutputList->ls();

    PostData(1, fOutputList);
}

//_________________________________________________
void AliAnalysisTaskTriggerEff::UserExec(Option_t *)
{
    fAODEvent = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!fAODEvent)
    {
        AliError("AOD event not found. Nothing done!");
        return;
    }

    TClonesArray *mcarray = dynamic_cast<TClonesArray *>(fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    //-----------------------------------------------
    // loop on MC generated event
    //-----------------------------------------------
    Int_t numtracks = fAODEvent->GetNumberOfTracks();
    Int_t nmu_rec = 0;
    Int_t GoodMuon_rec[10000];

    for (Int_t j = 0; j < numtracks; j++)
    {
        AliAODTrack *mu0 = (AliAODTrack *)fAODEvent->GetTrack(j);
        if (mu0->GetLabel() == -1)
            continue;
        AliAODMCParticle *mctrack0 = (AliAODMCParticle *)mcarray->At(mu0->GetLabel());
        if (mu0->GetLabel() != mctrack0->GetLabel())
            continue;
        Int_t PDG_mctrack0 = mctrack0->GetPdgCode();

        if (TMath::Abs(PDG_mctrack0) != 13)
            continue;

        if (mu0->Y() > -4.0 && mu0->Y() < -2.5)
            GoodMuon_rec[j] = kTRUE;

        if (!GoodMuon_rec[j])
            continue;

        // fPt_rec[nmu_rec] = mu0->Pt();
        // fY_rec[nmu_rec] = mu0->Y();
        // fEta_rec[nmu_rec] = mu0->Eta();
        // fMatchTrig_rec[nmu_rec] = mu0->GetMatchTrigger();
        // fMatchTrigChi2_rec[nmu_rec] = mu0->GetChi2MatchTrigger();
        // if (fMuonTrackCuts->IsSelected(mu0))
        //     fpDCA_rec[nmu_rec] = 1;q
        // fPhi_rec[nmu_rec] = mu0->Phi();

        if (GoodMuon_rec[j])
        {
            fHistPt_allmuon_alleta->Fill(mu0->Pt());

            //-------------------------All THR---------------------------------//

            if (mu0->Y() > -4.0 && mu0->Y() < -3.5)
            {
                fHistPt_allmuon_eta1_allpt->Fill(mu0->Pt());
            }
            else if (mu0->Y() > -3.5 && mu0->Y() < -3.0)
            {
                fHistPt_allmuon_eta2_allpt->Fill(mu0->Pt());
            }
            else if (mu0->Y() > -3.0 && mu0->Y() < -2.5)
            {
                fHistPt_allmuon_eta3_allpt->Fill(mu0->Pt());
            }

            //-------------------------Low THR---------------------------------//
            if (mu0->GetMatchTrigger() > 1)
            {
                fHistPt_allmuon_alleta_lowpt_thr->Fill(mu0->Pt());
                if ((mu0->Y() > -4.0 && mu0->Y() < -3.5))
                {
                    fHistPt_allmuon_eta1_lowpt_thr->Fill(mu0->Pt());
                }
                else if ((mu0->Y() > -3.5 && mu0->Y() < -3.0))
                {
                    fHistPt_allmuon_eta2_lowpt_thr->Fill(mu0->Pt());
                }
                else if ((mu0->Y() > -3.0 && mu0->Y() < -2.5))
                {
                    fHistPt_allmuon_eta3_lowpt_thr->Fill(mu0->Pt());
                }
            }
        }

        // fTheta_rec[nmu_rec] = mu0->Theta();
        // fCharge_rec[nmu_rec] = mu0->Charge();

        // Int_t mum_mctrack0 = mctrack0->GetMother();
        // AliAODMCParticle *mctrack0_mother = (AliAODMCParticle *)mcarray->At(mum_mctrack0);
        // Int_t PDG_mum_mctrack0 = mctrack0_mother->GetPdgCode();

        // if (TMath::Abs(PDG_mum_mctrack0) == 4)
        // {
        //     //            printf ("Muon from Charm Quark \n");
        //     fPDGmum_rec[nmu_rec] = IsPrompt(mctrack0_mother, mcarray);
        // }
        // else if (TMath::Abs(PDG_mum_mctrack0) > 400 && TMath::Abs(PDG_mum_mctrack0) < 500)
        // {
        //     //            printf ("Muon from Charm Meson \n");
        //     fPDGmum_rec[nmu_rec] = IsPrompt(mctrack0_mother, mcarray);
        // }
        // else if (TMath::Abs(PDG_mum_mctrack0) > 4000 && TMath::Abs(PDG_mum_mctrack0) < 5000)
        // {
        //     //            printf ("Muon from Charm Barion \n");
        //     fPDGmum_rec[nmu_rec] = IsPrompt(mctrack0_mother, mcarray);
        // }
        // else
        // {
        //     fPDGmum_rec[nmu_rec] = PDG_mum_mctrack0;
        // }
    }

    PostData(1, fOutputList);
}

//________________________________________________________________________

void AliAnalysisTaskTriggerEff::Terminate(Option_t *)
{
    printf("Terminate Task\n");
}
