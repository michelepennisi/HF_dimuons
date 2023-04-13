import numpy as np
import ROOT as rt
import sys
sys.path.insert(0, '/home/michele_pennisi/cernbox/HF_dimuons')

import graphic_functions as gf

def m_distrib_data():
    fIn = rt.TFile("~/dimuon_HF_pp/data/LHC18p/Hist_AOD/3_11_2022/HistResults_merged.root","READ")
    h_PtMdimu=fIn.Get("Dimuon/CMUL7/DQ_cut_match_LT_ULS/h_PtMdiMu_CMUL7_DQ_cut_match_LT_ULS_norm")
    h_Mdimu=h_PtMdimu.ProjectionY()
    h_Mdimu.Rebin(1)
    h_Mdimu.SetMarkerStyle(20)
    h_Mdimu.SetMarkerSize(1.3)
    h_Mdimu.SetLineWidth(2)
    h_Mdimu.SetLineColor(rt.kBlack)
    h_Mdimu.GetXaxis().SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})")
    h_Mdimu.GetYaxis().SetTitle("d#it{N}/d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}")
    # h_Mdimu.Scale(1.,"width")
    # for x in range(h_Mdimu.GetNbinsX()):
    #     print("%0.1f," % h_Mdimu.GetBinLowEdge(x+1),end="")

    array_low_bin=np.array([4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,11.5,12.0,12.5,13.0,13.5,14.0,15.0,16.0,18.0,20.5,22.0,24.0,27.0,30.0])

    h_Mdimu_REBINNED=h_Mdimu.Rebin(array_low_bin.size-1,("h_Mdimu_REBINNED"),array_low_bin)
    h_Mdimu_REBINNED.Scale(1.,"width")
    h_Mdimu_REBINNED.Smooth()
    # h_Mdimu_REBINNED.Draw()

    
    h_Mdimu_Upsilon=h_Mdimu.Clone("h_Mdimu_Upsilon")
    h_Mdimu_Upsilon.GetXaxis().SetRangeUser(9,11)
    h_Mdimu_Upsilon.Scale(1.,"width")
    h_Mdimu_Upsilon.SetMarkerColor(rt.kGray)
    h_Mdimu_Upsilon.SetLineColor(rt.kGray)
    h_Mdimu_Upsilon.Smooth()
    
    fill_histo=rt.TH1D("fill_histo","fill_histo",1,9,11)
    fill_histo.SetBinContent(1,1)
    max=h_Mdimu_REBINNED.GetMaximum()*20
    min=h_Mdimu_REBINNED.GetMaximum()*0.000002
    # fill_histo.SetBinError(1,5)
    fill_histo.SetFillColorAlpha(rt.kRed, 0.15);
    
    h_grid=h_Mdimu.Clone("h_grid")
    h_grid.Reset()

    gf.LoadStyle()
    h_grid.SetMaximum(max)
    h_grid.SetMinimum(min)
    canvas=gf.draw_single_hist(h_grid)
    canvas.cd()
    h_grid.Draw()
    h_Mdimu_REBINNED.Draw("PESAME")
    
    fill_histo.Draw("HistSAME")
    h_Mdimu_Upsilon.Draw("PESAME")

    letexTitle = rt.TLatex()
    letexTitle.SetNDC()
    letexTitle.SetTextFont(42)
    letexTitle.SetTextSize(0.035)
    letexTitle.DrawLatex(0.525,0.86,"ALICE Preliminary, pp #sqrt{#it{s}} = 13 TeV")
    letexTitle.DrawLatex(0.525,0.80,"2.5 < #it{#eta}_{#mu} < 4.0, 2.5 < #it{y}_{#mu#mu} < 4.0")

    l_low=rt.TLine(9,min,9,max)
    l_low.SetLineColor(rt.kRed)
    l_low.SetLineStyle(rt.kDashed)
    l_low.SetLineWidth(2)
    l_low.Draw()

    l_highw=rt.TLine(11,min,11,max)
    l_highw.SetLineColor(rt.kRed)
    l_highw.SetLineStyle(rt.kDashed)
    l_highw.SetLineWidth(2)
    l_highw.Draw()

    canvas.SaveAs("m_distrib_data.pdf")
    canvas.SaveAs("m_distrib_data.png")
    input()

m_distrib_data()