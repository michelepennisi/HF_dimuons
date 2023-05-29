#include "save_mc_output.h"

void save_mc_output(
    TString RunMode = "powheg_DY_mass_3_35",
    Int_t RunNumber = 294009,
    TString Task_Version = "Version1",
    Bool_t test = kFALSE,
    TString prefix_filename = "MCDimuHFTree")
{

    Set_Histograms();
    h_MPdg1Pdg2_DiMuon_Rec[0]->Draw();
}
