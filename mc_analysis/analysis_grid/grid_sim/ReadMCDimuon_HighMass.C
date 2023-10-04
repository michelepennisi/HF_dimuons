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
#include "AliAnalysisTaskDimuon_HighMass.h"
#endif
/// alice/sim/2022/LHC22b3/294925/AOD/
/// alice/cern.ch/user/m/mpennisi/powheg_jdl_sub_test/LHC18p
void ReadMCDimuon_HighMass(
    const char *RunMode = "full",
    Int_t RunNumber = 294013,
    TString Version = "Version_3_AliAODMuons",
    TString MC_type = "LHC23i1",
    TString GridDir = "/alice/sim/2023",
    // TString GridDir = "/alice/cern.ch/user/m/mpennisi/jira_test_charm",
    // TString GridDir = "/alice/cern.ch/user/m/mpennisi",
    TString AOD_origin = "Powheg",
    Bool_t usePhysicsSelection = kFALSE,
    TString DataPattern = "/AOD/*/AliAOD.Muons.root",
    TString AliPhysicsVersion = "vAN-20220204_ROOT6-1",
    Bool_t gridMerge = kTRUE)
{

  // LHC22b3 HF
  // LHC22c1 MB
  //  RunMode can be set to "test" to test running of the code
  //  RunMode can be set to "full" to produce Trees
  //  RunMode can be set to "terminate" to do the merging of the Trees

  //------------------------------------------------------------------------------------
  // Load paths
  //------------------------------------------------------------------------------------
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include ");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY");
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");

  //------------------------------------------------------------------------------------
  // Handler and analysis manager
  //------------------------------------------------------------------------------------

  // AliAnalysisManager *mgr = new AliAnalysisManager(Form("analysize MC from %s Simulation", MC_type.Data()));
  AliAnalysisManager *mgr = new AliAnalysisManager("analysize high mass dimuon MC on Grid");
  mgr->SetDebugLevel(AliLog::kError);
  AliLog::SetGlobalLogLevel(AliLog::kDebug);

  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  gInterpreter->LoadMacro("../CreateAlienHandler_HighMass.C");

  AliAnalysisAlien *alienHandler = reinterpret_cast<AliAnalysisAlien *>(gInterpreter->ProcessLine(Form("CreateAlienHandler_HighMass(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d,\"%s\",%d)", RunMode, Version.Data(), GridDir.Data(), MC_type.Data(), DataPattern.Data(), RunNumber, AliPhysicsVersion.Data(), gridMerge)));
  if (!alienHandler)
    return;
  mgr->SetGridHandler(alienHandler);

  AliAnalysisDataContainer *cinput1 = NULL;
  cinput1 = mgr->GetCommonInputContainer();

  //------------------------------------------------------------------------------------
  // Add task
  //------------------------------------------------------------------------------------
  gROOT->LoadMacro("AliAnalysisTaskDimuon_HighMass.cxx++g");
  AliAnalysisTaskDimuon_HighMass *MCTask = reinterpret_cast<AliAnalysisTaskDimuon_HighMass *>(gInterpreter->ProcessLine(Form(".x %s(\"%s\",%d)", "../AddTaskDimuon_HighMass.C", MC_type.Data(), RunNumber))); // I set by hand usePhysicsSelection=kTRUE
  MCTask->SetSaving_opt(Version);
  MCTask->SetAOD_origin(AOD_origin);
  mgr->AddTask(MCTask);
  if (usePhysicsSelection)
  {
    MCTask->SelectCollisionCandidates(AliVEvent::kAny);
  }

  //------------------------------------------------------------------------------------
  // Init analysis on GRID
  //------------------------------------------------------------------------------------
  if (!mgr->InitAnalysis())
    return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");
}
