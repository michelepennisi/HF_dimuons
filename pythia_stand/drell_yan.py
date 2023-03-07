import ROOT as rt
def LoadStyle():
    rt.gStyle.SetImageScaling(3.)
    font = 42
    rt.gStyle.SetFrameBorderMode(0)
    rt.gStyle.SetFrameFillColor(0)
    rt.gStyle.SetCanvasBorderMode(0)
    rt.gStyle.SetPadBorderMode(0)
    rt.gStyle.SetPadColor(10)
    rt.gStyle.SetCanvasColor(10)
    # rt.gStyle.SetTitleFillColor(10)
    # rt.gStyle.SetTitleBorderSize(1)
    rt.gStyle.SetStatColor(10)
    rt.gStyle.SetStatBorderSize(1)
    rt.gStyle.SetLegendBorderSize(1)
    rt.gStyle.SetDrawBorder(0)
    rt.gStyle.SetTextFont(font)
    rt.gStyle.SetStatFontSize(0.05)
    rt.gStyle.SetStatX(0.97)
    rt.gStyle.SetStatY(0.98)
    rt.gStyle.SetStatH(0.03)
    rt.gStyle.SetStatW(0.3)
    rt.gStyle.SetTickLength(0.02, "y")
    rt.gStyle.SetEndErrorSize(3)
    rt.gStyle.SetLabelSize(0.04, "xyz")
    rt.gStyle.SetLabelFont(font, "xyz")
    rt.gStyle.SetLabelOffset(0.01, "xyz")
    rt.gStyle.SetTitleFont(font, "xyz")
    rt.gStyle.SetTitleOffset(0.9, "x")
    rt.gStyle.SetTitleOffset(1.02, "y")
    rt.gStyle.SetTitleSize(0.04, "xyz")
    rt.gStyle.SetMarkerSize(1.3)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetEndErrorSize(0)
    rt.gStyle.SetCanvasPreferGL(rt.kTRUE)
    rt.gStyle.SetHatchesSpacing(0.5)



def canvas_style():
    canvas = rt.TCanvas("c", "c", 1000, 1100)
    canvas.cd()
    canvas.SetRightMargin(0.05)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetBorderSize(0)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.SetLeftMargin(0.126)
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.12)
    canvas.SetFrameBorderMode(0)
    canvas.SetFrameBorderMode(0)

    return canvas

def main(var):
    LoadStyle()
    fIn_drell=rt.TFile("sim/Zperevent_mgt4.root","READ")   
    fIn_drell.ls()

    h_MDimu_DY=fIn_drell.Get("h_dimu%s_OS_Z" % var)
    # h_PtDimu_DY.Rebin(10) 
    if(var.find("mass")!=-1):
        h_MDimu_DY.GetXaxis().SetRangeUser(4,30)
    
    # xbins = np.array((0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7,8,9,10,11.5,13,15,22.5,30.0))
    h_MDimu_DY.Rebin(2)
    h_MDimu_DY.SetMarkerStyle(20)
    h_MDimu_DY.SetMarkerColor(rt.kOrange+2)
    h_MDimu_DY.SetLineColor(rt.kOrange+2)
    h_MDimu_DY.SetLineWidth(3)
    h_MDimu_DY.Scale(1.,"width")

    fIn_HFsim = rt.TFile("~/cernbox/output_HF_dimuons/mc_analysis_output/Hist_fromSim/Version1/HF/HistLite_HF_MCDimuHFTree_merged.root", "READ")
    h_Nevents_HFsim = fIn_HFsim.Get("h_Nevents")
    nEvents_HFsim = h_Nevents_HFsim.GetBinContent(2)

    h_PtM_DimuHF_HFsim = fIn_HFsim.Get("DiMuon/M4/Gen_DQ_cut/ULS/h_PtMDiMu_M4_Gen_DQ_cut_ULS_fromHF")
    if(var.find("mass")!=-1):
        h_M_DimuHF_HFsim = h_PtM_DimuHF_HFsim.ProjectionY()
        h_M_DimuHF_HFsim.SetName("h_M_DimuHF_HFsim")
        h_M_DimuHF_HFsim.Rebin(20)
    elif(var.find("pt")!=-1):
        h_M_DimuHF_HFsim = h_PtM_DimuHF_HFsim.ProjectionX()
        h_M_DimuHF_HFsim.SetName("h_Pt_DimuHF_HFsim")
        h_M_DimuHF_HFsim.Rebin(20)

    print("N dimu from HF %0.2f\n",h_M_DimuHF_HFsim.GetEntries())
    
    h_M_DimuHF_HFsim.Scale(1. / (216*nEvents_HFsim), "width")
    h_M_DimuHF_HFsim.SetMarkerStyle(24)
    h_M_DimuHF_HFsim.SetMarkerColor(rt.kBlue+2)
    h_M_DimuHF_HFsim.SetLineColor(rt.kBlue+2)
    h_M_DimuHF_HFsim.SetLineWidth(3)
    # h_Pt_DimuHF_HFsim.SetMinimum(h_PtDimu_DY.GetMinimum()*0.02)
    
    canvas_comparison=canvas_style()
    canvas_comparison.cd()
    rt.gStyle.SetOptTitle(0)
    rt.gPad.SetLogy()
    if(var.find("mass")!=-1):
        h_M_DimuHF_HFsim.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})")
        h_M_DimuHF_HFsim.GetYaxis().SetTitle("d#it{N}_{#mu^{#plus}#mu^{#minus}}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}")
    elif(var.find("pt")!=-1):
        h_M_DimuHF_HFsim.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
        h_M_DimuHF_HFsim.GetYaxis().SetTitle("d#it{N}_{#mu^{#plus}#mu^{#minus}}/d#it{p}_{T} (GeV/#it{c})^{-1}")
    
    h_M_DimuHF_HFsim.GetXaxis().SetTitleOffset(1.6)
    h_M_DimuHF_HFsim.GetXaxis().SetTitleSize(0.035)
    h_M_DimuHF_HFsim.GetXaxis().SetLabelSize(0.03)
    
    h_M_DimuHF_HFsim.GetYaxis().SetTitleOffset(1.6)
    h_M_DimuHF_HFsim.GetYaxis().SetTitleSize(0.035)
    h_M_DimuHF_HFsim.GetYaxis().SetLabelSize(0.035)
    h_M_DimuHF_HFsim.Draw("PE")
    h_MDimu_DY.Draw("PESAME")

    Leged = rt.TLegend(0.18, 0.18, 0.325, 0.425, " ", "brNDC")
    

    Leged.AddEntry(h_M_DimuHF_HFsim, "#mu^{#plus}#mu^{#minus} #leftarrow HF", "EP")
    Leged.AddEntry(h_MDimu_DY, "#mu^{#plus}#mu^{#minus} #leftarrow DY", "EP")
    
    Leged.SetBorderSize(0)
    Leged.SetFillColor(10)
    Leged.SetFillStyle(1)
    Leged.SetLineStyle(0)
    Leged.SetLineColor(0)
    Leged.SetTextFont(42)
    Leged.SetTextSize(0.035)
    Leged.Draw("SAME")
    rt.gPad.Update()
    canvas_comparison.SaveAs("DY_HF_comp_%s.png" % var)
    input()

main("mass")
main("pt")
