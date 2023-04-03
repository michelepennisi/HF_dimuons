import ROOT as rt

def acc_table():
    fIn=rt.TFile("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/root_files/test/HF_Analysis_MCsim_294241.root","READ")
    fIn.cd()

    h_PtYPdg_DiMu_Rec_Meson_ULS=fIn.Get("DiMu_Rec/h_PtYPdg_DiMu_Rec_Meson_ULS_M4cut")
    h_PtYPdg_DiMu_Rec_Meson_ULS.GetZaxis().SetRange(1,1)
    h_PtY_DiMu_Rec_Charm_Meson_ULS=h_PtYPdg_DiMu_Rec_Meson_ULS.Project3D("xye")
    h_PtY_DiMu_Rec_Charm_Meson_ULS.RebinX(20)
    h_PtY_DiMu_Rec_Charm_Meson_ULS.RebinY(20)

    h_PtYPdg_DiMu_Gen_Meson_ULS=fIn.Get("DiMu_Gen/h_PtYPdg_DiMu_Gen_Meson_ULS_M4cut")
    h_PtYPdg_DiMu_Gen_Meson_ULS.GetZaxis().SetRange(1,1)
    h_PtY_DiMu_Gen_Charm_Meson_ULS=h_PtYPdg_DiMu_Gen_Meson_ULS.Project3D("xye")
    h_PtY_DiMu_Gen_Charm_Meson_ULS.RebinX(20)
    h_PtY_DiMu_Gen_Charm_Meson_ULS.RebinY(20)

    h_PtYPdg_DiMu_Rec_Meson_ULS=fIn.Get("DiMu_Rec/h_PtYPdg_DiMu_Rec_Meson_ULS_M4cut")
    h_PtYPdg_DiMu_Rec_Meson_ULS.GetZaxis().SetRange(2,2)
    h_PtY_DiMu_Rec_Beauty_Meson_ULS=h_PtYPdg_DiMu_Rec_Meson_ULS.Project3D("xye")
    h_PtY_DiMu_Rec_Beauty_Meson_ULS.RebinX(20)
    h_PtY_DiMu_Rec_Beauty_Meson_ULS.RebinY(20)

    h_PtYPdg_DiMu_Gen_Meson_ULS=fIn.Get("DiMu_Gen/h_PtYPdg_DiMu_Gen_Meson_ULS_M4cut")
    h_PtYPdg_DiMu_Gen_Meson_ULS.GetZaxis().SetRange(2,2)
    h_PtY_DiMu_Gen_Beauty_Meson_ULS=h_PtYPdg_DiMu_Gen_Meson_ULS.Project3D("xye")
    h_PtY_DiMu_Gen_Beauty_Meson_ULS.RebinX(20)
    h_PtY_DiMu_Gen_Beauty_Meson_ULS.RebinY(20)

    canvas_charm=rt.TCanvas("canvas_charm","canvas_charm",1000,1000)
    canvas_charm.Divide(3,1)
    canvas_charm.cd(1)

    h_PtY_DiMu_Rec_Charm_Meson_ULS.Draw("COLZ")
    canvas_charm.cd(2)
    h_PtY_DiMu_Gen_Charm_Meson_ULS.Draw("COLZ")

    canvas_charm.cd(3)
    h_PtY_DiMu_ACC_Charm_Meson_ULS=h_PtY_DiMu_Rec_Charm_Meson_ULS.Clone("h_PtY_DiMu_ACC_Charm_Meson_ULS")
    h_PtY_DiMu_ACC_Charm_Meson_ULS.Divide(h_PtY_DiMu_Gen_Charm_Meson_ULS)
    h_PtY_DiMu_ACC_Charm_Meson_ULS.Draw("COLZ")

    canvas_Beauty=rt.TCanvas("canvas_Beauty","canvas_Beauty",1000,1000)
    canvas_Beauty.Divide(3,1)
    canvas_Beauty.cd(1)

    h_PtY_DiMu_Rec_Beauty_Meson_ULS.Draw("COLZ")
    canvas_Beauty.cd(2)
    h_PtY_DiMu_Gen_Beauty_Meson_ULS.Draw("COLZ")

    canvas_Beauty.cd(3)
    h_PtY_DiMu_ACC_Beauty_Meson_ULS=h_PtY_DiMu_Rec_Beauty_Meson_ULS.Clone("h_PtY_DiMu_ACC_Beauty_Meson_ULS")
    h_PtY_DiMu_ACC_Beauty_Meson_ULS.Divide(h_PtY_DiMu_Gen_Beauty_Meson_ULS)
    h_PtY_DiMu_ACC_Beauty_Meson_ULS.Draw("COLZ")

    input()

acc_table()