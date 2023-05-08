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
#include "AliAnalysisTaskDimuonHF.h"
#endif
/// alice/sim/2022/LHC22b3/294925/AOD/
void ReadMCDimuonHF_Grid(
    const char *RunMode = "test",
    Int_t RunNumber = 294743,
    TString MC_type = "Z",
    Bool_t usePhysicsSelection = kFALSE,
    TString GridDir = "/alice/cern.ch/user/m/mpennisi/powheg_test/LHC18p",
    TString DataPattern = "/*/AliAOD.root",
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
  AliAnalysisManager *mgr = new AliAnalysisManager("analysize MC from Powheg Sim Simulation for DY process");
  mgr->SetDebugLevel(AliLog::kError);
  AliLog::SetGlobalLogLevel(AliLog::kDebug);

  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  gInterpreter->LoadMacro("CreateAlienHandler_Grid.C");

  AliAnalysisGrid *alienHandler = reinterpret_cast<AliAnalysisGrid *>(gInterpreter->ProcessLine(Form("CreateAlienHandler_Grid(\"%s\",\"%s\",\"%s\",\"%s\",%d,\"%s\",%d)", RunMode, GridDir.Data(), MC_type.Data(), DataPattern.Data(), RunNumber, AliPhysicsVersion.Data(), gridMerge)));
  if (!alienHandler)
    return;
  mgr->SetGridHandler(alienHandler);

  AliAnalysisDataContainer *cinput1 = NULL;
  cinput1 = mgr->GetCommonInputContainer();

  //------------------------------------------------------------------------------------
  // Add task
  //------------------------------------------------------------------------------------
  gROOT->LoadMacro("AliAnalysisTaskDimuonHF.cxx++g");
  AliAnalysisTaskDimuonHF *MCEmbeddingTask = reinterpret_cast<AliAnalysisTaskDimuonHF *>(gInterpreter->ProcessLine(Form(".x %s(\"%s\",%d)", "./AddTaskDimuonHF.C", MC_type.Data(),RunNumber))); // I set by hand usePhysicsSelection=kTRUE
  mgr->AddTask(MCEmbeddingTask);
  if (usePhysicsSelection)
  {
    MCEmbeddingTask->SelectCollisionCandidates(AliVEvent::kAny);
  }

  //------------------------------------------------------------------------------------
  // Init analysis on GRID
  //------------------------------------------------------------------------------------
  if (!mgr->InitAnalysis())
    return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");
}
