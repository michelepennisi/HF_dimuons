#if !defined(__CINT__) || defined(__CLING__)
#include "TROOT.h"
#include "AliAnalysisAlien.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TChain.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskDimuon_HighMass.h"
#endif
void ReadMCDimuon_HighMass_Local(const char *RunMode = "terminate",
                                 Int_t RunNumber = 294154,
                                 TString Version = "Version_5_AliAOD_skimmed_fwd",
                                 TString MC_type = "powheg_beauty_nocut_Version_5_AliAOD_withHF_Q",
                                 TString AOD_origin = "Powheg")
{
    // header location
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/include ");
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY");
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("analysize high mass dimuon MC on Grid");
    mgr->SetDebugLevel(AliLog::kError);
    AliLog::SetGlobalLogLevel(AliLog::kDebug);

    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    // compile the class (locally) with debug symbols
    gInterpreter->LoadMacro("AliAnalysisTaskDimuon_HighMass.cxx++g");
    AliAnalysisTaskDimuon_HighMass *MCTask = reinterpret_cast<AliAnalysisTaskDimuon_HighMass *>(gInterpreter->ProcessLine(Form(".x %s(\"%s\",%d)", "AddTaskDimuon_HighMass_Local.C", MC_type.Data(), RunNumber)));
    MCTask->SetSaving_opt(Version);
    MCTask->SetAOD_origin(AOD_origin);
    // if you want to run locally, we need to define some input
    TChain *chain = new TChain("aodTree");
    // chain->Add("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/LHC23i1_data/AliAOD.root");
    for (Int_t i = 0; i < 100; i++)
    {
        chain->Add(Form("/home/michele_pennisi/cernbox/HF_dimuons/mc_analysis/analysis_grid/grid_sim/%s/AOD/%d/AliAOD.root", MC_type.Data(), i));
    }

    // start the analysis locally
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
}