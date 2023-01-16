#ifndef AliAnalysisTaskTriggerEff_H
#define AliAnalysisTaskTriggerEff_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TObjArray;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;
class AliMuonTrackCuts;

class AliAnalysisTaskTriggerEff : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskTriggerEff();
    AliAnalysisTaskTriggerEff(const char *name);
    virtual ~AliAnalysisTaskTriggerEff();

    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *);
    virtual void NotifyRun();

    void SetBeamEnergy(Double_t en) { fBeamEnergy = en; }
    void SetAnalysisType(const char *type) { fkAnalysisType = type; }
    void SetPeriod(TString period) { fPeriod = period; }

private:
    AliAnalysisTaskTriggerEff(const AliAnalysisTaskTriggerEff &);
    AliAnalysisTaskTriggerEff &operator=(const AliAnalysisTaskTriggerEff &);

    // protected:

    Double_t fBeamEnergy;       // Energy of the beam (required for the CS angle)
    const char *fkAnalysisType; // ESD or AOD based analysis
    TString fPeriod;            // period
    AliAODEvent *fAODEvent;     //! AOD event  //tolgo !
    AliMuonTrackCuts *fMuonTrackCuts;
    Double_t fPercentV0M;

    TList *fOutputList;
    TH1D *fHistPt_allmuon_alleta;
    TH1D *fHistPt_allmuon_eta1_allpt;
    TH1D *fHistPt_allmuon_eta2_allpt;
    TH1D *fHistPt_allmuon_eta3_allpt;

    TH1D *fHistPt_allmuon_alleta_lowpt_thr;
    TH1D *fHistPt_allmuon_eta1_lowpt_thr;
    TH1D *fHistPt_allmuon_eta2_lowpt_thr;
    TH1D *fHistPt_allmuon_eta3_lowpt_thr;

    ClassDef(AliAnalysisTaskTriggerEff, 1);
};

#endif
