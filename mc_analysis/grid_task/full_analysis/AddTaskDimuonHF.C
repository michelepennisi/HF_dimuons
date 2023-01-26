AliAnalysisTaskDimuonHF *AddTaskDimuonHF(TString MC_type, Int_t RunNumber = 294009)
{

   //****************************************************************************************
   // Add task class.
   // The attached class prepares a MC tree (embedding) with generated and reconstructed muons/dimuons
   // Michele Pennisi
   //****************************************************************************************

   printf("Creating Task to creat a Muon/Dimuon tree from Embedding\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr)
   {
      ::Error("AddTaskDimuonHF", "No analysis manager to connect to.");
      return NULL;
   }
   TString fnameout;
   fnameout.Form("%s_MCDimuHFTree_%d.root", MC_type.Data(), RunNumber);
   printf("Fnameout = %s\n", fnameout.Data());

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("MCTree", TTree::Class(), AliAnalysisManager::kOutputContainer, fnameout);

   AliAnalysisTaskDimuonHF *MCQuarkoniumTask = new AliAnalysisTaskDimuonHF("AliAnalysisTaskDimuHFTreeMC");

   // MCQuarkoniumTask->SetBeamEnergy(5.02);
   // MCQuarkoniumTask->SetResonance(resonance);
   mgr->AddTask(MCQuarkoniumTask);

   mgr->ConnectInput(MCQuarkoniumTask, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(MCQuarkoniumTask, 1, coutput1);

   return MCQuarkoniumTask;
}
