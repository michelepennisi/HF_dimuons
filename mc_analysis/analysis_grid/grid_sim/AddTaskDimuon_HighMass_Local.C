AliAnalysisTaskDimuon_HighMass *AddTaskDimuon_HighMass_Local(TString MC_type, Int_t RunNumber = 294009)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    // resolve the name of the output file
    TString fileName;
    fileName.Form("%s_MCDimuHFTree_%d.root", MC_type.Data(), RunNumber);

    // now we create an instance of your task
    AliAnalysisTaskDimuon_HighMass *task = new AliAnalysisTaskDimuon_HighMass("AliAnalysisTaskDimuHFTreeMC");

    // add your task to the manager
    mgr->AddTask(task);

    // connect the manager to your task
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task, 1, mgr->CreateContainer("MyOutputContainer", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    // important: return a pointer to your task
    return task;
}