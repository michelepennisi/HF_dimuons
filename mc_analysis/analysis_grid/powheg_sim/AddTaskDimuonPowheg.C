AliAnalysisTaskDimuonPowheg *AddTaskDimuonPowheg(TString MC_type, Int_t RunNumber = 294009)
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
      ::Error("AddTaskDimuonPowheg", "No analysis manager to connect to.");
      return NULL;
   }
   TString fnameout;
   fnameout.Form("%s_MCDimuHFTree_%d.root", MC_type.Data(), RunNumber);   
   // fnameout.Form("pow_Z_Sim_MCDimuHFTree_%d.root", RunNumber);   
   printf("Fnameout = %s\n", fnameout.Data());

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("MCTree", TTree::Class(), AliAnalysisManager::kOutputContainer, fnameout);

   AliAnalysisTaskDimuonPowheg *MC_DiMu = new AliAnalysisTaskDimuonPowheg("AliAnalysisTaskDimuHFTreeMC");

   mgr->AddTask(MC_DiMu);

   mgr->ConnectInput(MC_DiMu, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(MC_DiMu, 1, coutput1);

   return MC_DiMu;
}
