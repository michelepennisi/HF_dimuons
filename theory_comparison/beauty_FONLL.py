import pandas as pd
import numpy as np

import ROOT


def LoadStyle_canvas():
    ROOT.gStyle.SetImageScaling(3.)
    font = 42
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetFrameFillColor(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(10)
    ROOT.gStyle.SetCanvasColor(10)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleBorderSize(1)
    ROOT.gStyle.SetStatColor(10)
    ROOT.gStyle.SetStatBorderSize(1)
    ROOT.gStyle.SetLegendBorderSize(1)
    ROOT.gStyle.SetDrawBorder(0)
    ROOT.gStyle.SetTextFont(font)
    ROOT.gStyle.SetStatFontSize(0.05)
    ROOT.gStyle.SetStatX(0.97)
    ROOT.gStyle.SetStatY(0.98)
    ROOT.gStyle.SetStatH(0.03)
    ROOT.gStyle.SetStatW(0.3)
    ROOT.gStyle.SetTickLength(0.02, "y")
    ROOT.gStyle.SetEndErrorSize(3)
    ROOT.gStyle.SetLabelSize(0.04, "xyz")
    ROOT.gStyle.SetLabelFont(font, "xyz")
    ROOT.gStyle.SetLabelOffset(0.01, "xyz")
    ROOT.gStyle.SetTitleFont(font, "xyz")
    ROOT.gStyle.SetTitleOffset(0.9, "x")
    ROOT.gStyle.SetTitleOffset(1.02, "y")
    ROOT.gStyle.SetTitleSize(0.04, "xyz")
    ROOT.gStyle.SetMarkerSize(1.3)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetEndErrorSize(0)
    ROOT.gStyle.SetCanvasPreferGL(ROOT.kTRUE)
    ROOT.gStyle.SetHatchesSpacing(0.5)

    canvas = ROOT.TCanvas("c", "c", 1600, 1200)
    canvas.cd()
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.10)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetBorderSize(0)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.1518219)
    canvas.SetFrameBorderMode(0)
    canvas.SetFrameBorderMode(0)

    return canvas


def main():

    y_fwd = np.array([-3.25, 3.25])
    dy_fwd = np.array([0.75, 0.75])

    ds_bb_dy_PYTHIA_fwd = np.array([24.6, 24.6])
    stat_ds_bb_dy_PYTHIA_fwd = np.array([0.9, 0.9])
    syst_ds_bb_dy_PYTHIA_fwd = np.array([7.5, 7.5])

    bb_cs_PYTHIA_fwd = ROOT.TGraphMultiErrors("bb_cs_PYTHIA_fwd", "TGraphMultiErrors Example",
                                              2, y_fwd, ds_bb_dy_PYTHIA_fwd, dy_fwd, dy_fwd, stat_ds_bb_dy_PYTHIA_fwd, stat_ds_bb_dy_PYTHIA_fwd)
    bb_cs_PYTHIA_fwd.AddYError(
        2, syst_ds_bb_dy_PYTHIA_fwd, syst_ds_bb_dy_PYTHIA_fwd)
    bb_cs_PYTHIA_fwd.SetMarkerStyle(20)
    bb_cs_PYTHIA_fwd.SetMarkerColor(ROOT.kRed)
    bb_cs_PYTHIA_fwd.SetLineColor(ROOT.kRed)
    bb_cs_PYTHIA_fwd.SetLineWidth(3)
    bb_cs_PYTHIA_fwd.GetAttLine(0).SetLineColor(ROOT.kRed)
    bb_cs_PYTHIA_fwd.GetAttLine(0).SetLineWidth(3)
    bb_cs_PYTHIA_fwd.GetAttLine(1).SetLineColor(ROOT.kRed)
    bb_cs_PYTHIA_fwd.GetAttLine(1).SetLineWidth(3)
    bb_cs_PYTHIA_fwd.GetAttFill(1).SetFillStyle(0)

    y_mid = np.array([0, 0.000001])
    dy_mid = np.array([0.5, 0])

    ds_bb_dy_PYTHIA_mid = np.array([79, 0.0])
    stat_ds_bb_dy_PYTHIA_mid = np.array([14, 0.0])
    syst_ds_bb_dy_PYTHIA_mid = np.array([11, 0.0])

    bb_cs_PYTHIA_mid = ROOT.TGraphMultiErrors("bb_cs_PYTHIA_mid", "TGraphMultiErrors Example",2, y_mid, ds_bb_dy_PYTHIA_mid, dy_mid, dy_mid, stat_ds_bb_dy_PYTHIA_mid, stat_ds_bb_dy_PYTHIA_mid)
    bb_cs_PYTHIA_mid.AddYError(2, syst_ds_bb_dy_PYTHIA_mid, syst_ds_bb_dy_PYTHIA_mid)
    bb_cs_PYTHIA_mid.SetMarkerStyle(24)
    bb_cs_PYTHIA_mid.SetMarkerColor(ROOT.kMagenta+2)
    bb_cs_PYTHIA_mid.SetLineColor(ROOT.kMagenta+2)
    bb_cs_PYTHIA_mid.SetLineWidth(3)
    bb_cs_PYTHIA_mid.GetAttLine(0).SetLineColor(ROOT.kMagenta+2)
    bb_cs_PYTHIA_mid.GetAttLine(0).SetLineWidth(3)
    bb_cs_PYTHIA_mid.GetAttLine(1).SetLineColor(ROOT.kMagenta+2)
    bb_cs_PYTHIA_mid.GetAttLine(1).SetLineWidth(3)
    bb_cs_PYTHIA_mid.GetAttFill(1).SetFillStyle(0)

    df_FONLL_bb_NNPDF = pd.read_csv('FONLL_bb_Pt_0_30_NNPDF.txt', sep=' ')
    y_FONLL_NNPDF = df_FONLL_bb_NNPDF["y"]
    low_ey_FONLL_NNPDF = df_FONLL_bb_NNPDF["low_y"]
    up_ey_FONLL_NNPDF = df_FONLL_bb_NNPDF["up_y"]
    ds_dy_bb_FONLL_NNPDF = df_FONLL_bb_NNPDF["central"] * 1e-06
    min_ds_dy_bb_FONLL_NNPDF = df_FONLL_bb_NNPDF["min"] * 1e-06
    max_ds_dy_bb_FONLL_NNPDF = df_FONLL_bb_NNPDF["max"] * 1e-06
    min_sc_ds_dy_bb_FONLL_NNPDF = df_FONLL_bb_NNPDF["min_sc"] * 1e-06
    max_sc_ds_dy_bb_FONLL_NNPDF = df_FONLL_bb_NNPDF["max_sc"] * 1e-06
    min_mass_ds_dy_bb_FONLL_NNPDF = df_FONLL_bb_NNPDF["min_mass"] * 1e-06
    max_mass_ds_dy_bb_FONLL_NNPDF = df_FONLL_bb_NNPDF["max_mass"] * 1e-06
    min_pdf_ds_dy_bb_FONLL_NNPDF = df_FONLL_bb_NNPDF["min_pdf"] * 1e-06
    max_pdf_ds_dy_bb_FONLL_NNPDF = df_FONLL_bb_NNPDF["max_pdf"] * 1e-06

    Arr_y_FONLL_NNPDF = y_FONLL_NNPDF.to_numpy()
    Arr_low_ey_FONLL_NNPDF = low_ey_FONLL_NNPDF.to_numpy()
    Arr_up_ey_FONLL_NNPDF = up_ey_FONLL_NNPDF.to_numpy()
    Arr_ds_dy_bb_FONLL_NNPDF = ds_dy_bb_FONLL_NNPDF.to_numpy()
    Arr_min_ds_dy_bb_FONLL_NNPDF = min_ds_dy_bb_FONLL_NNPDF.to_numpy()
    Arr_max_ds_dy_bb_FONLL_NNPDF = max_ds_dy_bb_FONLL_NNPDF.to_numpy()
    Arr_min_sc_ds_dy_bb_FONLL_NNPDF = min_sc_ds_dy_bb_FONLL_NNPDF.to_numpy()
    Arr_max_sc_ds_dy_bb_FONLL_NNPDF = max_sc_ds_dy_bb_FONLL_NNPDF.to_numpy()
    Arr_min_mass_ds_dy_bb_FONLL_NNPDF = min_mass_ds_dy_bb_FONLL_NNPDF.to_numpy()
    Arr_max_mass_ds_dy_bb_FONLL_NNPDF = max_mass_ds_dy_bb_FONLL_NNPDF.to_numpy()
    Arr_min_pdf_ds_dy_bb_FONLL_NNPDF = min_pdf_ds_dy_bb_FONLL_NNPDF.to_numpy()
    Arr_max_pdf_ds_dy_bb_FONLL_NNPDF = max_pdf_ds_dy_bb_FONLL_NNPDF.to_numpy()

    df_FONLL_bb_CTEQ6 = pd.read_csv('FONLL_bb_Pt_0_30_CTEQ6.txt', sep=' ')
    y_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["y"]
    low_ey_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["low_y"]
    up_ey_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["up_y"]
    ds_dy_bb_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["central"] * 1e-06
    min_ds_dy_bb_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["min"] * 1e-06
    max_ds_dy_bb_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["max"] * 1e-06
    min_sc_ds_dy_bb_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["min_sc"] * 1e-06
    max_sc_ds_dy_bb_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["max_sc"] * 1e-06
    min_mass_ds_dy_bb_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["min_mass"] * 1e-06
    max_mass_ds_dy_bb_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["max_mass"] * 1e-06
    min_pdf_ds_dy_bb_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["min_pdf"] * 1e-06
    max_pdf_ds_dy_bb_FONLL_CTEQ6 = df_FONLL_bb_CTEQ6["max_pdf"] * 1e-06

    Arr_y_FONLL_CTEQ6 = y_FONLL_CTEQ6.to_numpy()
    Arr_ds_dy_bb_FONLL_CTEQ6 = ds_dy_bb_FONLL_CTEQ6.to_numpy()
    Arr_low_ey_FONLL_CTEQ6 = low_ey_FONLL_CTEQ6.to_numpy()
    Arr_up_ey_FONLL_CTEQ6 = up_ey_FONLL_CTEQ6.to_numpy()
    Arr_min_ds_dy_bb_FONLL_CTEQ6 = min_ds_dy_bb_FONLL_CTEQ6.to_numpy()
    Arr_max_ds_dy_bb_FONLL_CTEQ6 = max_ds_dy_bb_FONLL_CTEQ6.to_numpy()
    Arr_min_sc_ds_dy_bb_FONLL_CTEQ6 = min_sc_ds_dy_bb_FONLL_CTEQ6.to_numpy()
    Arr_max_sc_ds_dy_bb_FONLL_CTEQ6 = max_sc_ds_dy_bb_FONLL_CTEQ6.to_numpy()
    Arr_min_mass_ds_dy_bb_FONLL_CTEQ6 = min_mass_ds_dy_bb_FONLL_CTEQ6.to_numpy()
    Arr_max_mass_ds_dy_bb_FONLL_CTEQ6 = max_mass_ds_dy_bb_FONLL_CTEQ6.to_numpy()
    Arr_min_pdf_ds_dy_bb_FONLL_CTEQ6 = min_pdf_ds_dy_bb_FONLL_CTEQ6.to_numpy()
    Arr_max_pdf_ds_dy_bb_FONLL_CTEQ6 = max_pdf_ds_dy_bb_FONLL_CTEQ6.to_numpy()

    Error_y = np.empty(16, dtype=object)

    canvas_bb_cs_NNPDF = LoadStyle_canvas()
    canvas_bb_cs_NNPDF.SetName("canvas_bb_cs_NNPDF")
    canvas_bb_cs_NNPDF.cd()
    FONLL_bb_cs_NNPDF_minmaxerror = ROOT.TGraphAsymmErrors(17, Arr_y_FONLL_NNPDF, Arr_ds_dy_bb_FONLL_NNPDF, Arr_low_ey_FONLL_NNPDF,
                                                           Arr_up_ey_FONLL_NNPDF, Arr_min_ds_dy_bb_FONLL_NNPDF, Arr_max_ds_dy_bb_FONLL_NNPDF)
    FONLL_bb_cs_NNPDF_minmaxerror.SetName("FONLL_bb_cs_NNPDF_minmaxerror")
    FONLL_bb_cs_NNPDF_minmaxerror.SetLineWidth(2)
    FONLL_bb_cs_NNPDF_minmaxerror.SetFillColorAlpha(ROOT.kGray, 0.7)
    FONLL_bb_cs_NNPDF_minmaxerror.SetFillColorAlpha(ROOT.kGray, 0.7)
    
    FONLL_bb_cs_NNPDF_scalerror = ROOT.TGraphAsymmErrors(17, Arr_y_FONLL_NNPDF, Arr_ds_dy_bb_FONLL_NNPDF, Arr_low_ey_FONLL_NNPDF,
                                                         Arr_up_ey_FONLL_NNPDF, Arr_min_sc_ds_dy_bb_FONLL_NNPDF, Arr_max_sc_ds_dy_bb_FONLL_NNPDF)
    FONLL_bb_cs_NNPDF_scalerror.SetName("FONLL_bb_cs_NNPDF_scalerror")
    FONLL_bb_cs_NNPDF_scalerror.SetFillColorAlpha(ROOT.kCyan-10, 0.2)
    FONLL_bb_cs_NNPDF_scalerror.SetFillColorAlpha(ROOT.kCyan-10, 0.2)
    FONLL_bb_cs_NNPDF_scalerror.SetLineWidth(2)
    

    FONLL_bb_cs_NNPDF_masserror = ROOT.TGraphAsymmErrors(17, Arr_y_FONLL_NNPDF, Arr_ds_dy_bb_FONLL_NNPDF, Arr_low_ey_FONLL_NNPDF,
                                                         Arr_up_ey_FONLL_NNPDF, Arr_min_mass_ds_dy_bb_FONLL_NNPDF, Arr_max_mass_ds_dy_bb_FONLL_NNPDF)
    FONLL_bb_cs_NNPDF_masserror.SetName("FONLL_bb_cs_NNPDF_masserror")
    FONLL_bb_cs_NNPDF_masserror.SetMarkerStyle(20)
    FONLL_bb_cs_NNPDF_masserror.SetFillColorAlpha(ROOT.kGreen, 0.2)
    FONLL_bb_cs_NNPDF_masserror.SetFillColorAlpha(ROOT.kGreen, 0.2)
    FONLL_bb_cs_NNPDF_masserror.SetLineWidth(2)
    

    FONLL_bb_cs_NNPDF_pdferror = ROOT.TGraphAsymmErrors(17, Arr_y_FONLL_NNPDF, Arr_ds_dy_bb_FONLL_NNPDF, Arr_low_ey_FONLL_NNPDF,
                                                         Arr_up_ey_FONLL_NNPDF, Arr_min_pdf_ds_dy_bb_FONLL_NNPDF, Arr_max_pdf_ds_dy_bb_FONLL_NNPDF)
    FONLL_bb_cs_NNPDF_pdferror.SetName("FONLL_bb_cs_NNPDF_pdferror")
    FONLL_bb_cs_NNPDF_pdferror.SetMarkerStyle(20)
    FONLL_bb_cs_NNPDF_pdferror.SetFillColorAlpha(ROOT.kMagenta-9, 0.2)
    FONLL_bb_cs_NNPDF_pdferror.SetFillColorAlpha(ROOT.kMagenta-9, 0.2)
    FONLL_bb_cs_NNPDF_pdferror.SetLineWidth(2)
    

    bb_bar_cs_NNPDF = ROOT.TMultiGraph()
    bb_bar_cs_NNPDF. Add(FONLL_bb_cs_NNPDF_scalerror)
    bb_bar_cs_NNPDF. Add(FONLL_bb_cs_NNPDF_masserror)
    bb_bar_cs_NNPDF. Add(FONLL_bb_cs_NNPDF_pdferror)
    bb_bar_cs_NNPDF. Add(FONLL_bb_cs_NNPDF_minmaxerror)
    ROOT.gPad.SetLogy()
    
    bb_bar_cs_NNPDF.Add(bb_cs_PYTHIA_fwd, "APS; Z ; 5 s=0.5")
    bb_bar_cs_NNPDF.Add(bb_cs_PYTHIA_mid, "APS; Z ; 5 s=0.5")
    bb_bar_cs_NNPDF.Draw("A3")
    # FONLL_bb_cs_NNPDF_minmaxerror.SetLineColorAlpha(ROOT.kGray, 1)
    # FONLL_bb_cs_NNPDF_minmaxerror_bis.Draw("C PLC PFC SAME")
    bb_bar_cs_NNPDF.GetXaxis().SetTitle("#it{y}")
    bb_bar_cs_NNPDF.GetXaxis().SetTitleOffset(1.2)
    bb_bar_cs_NNPDF.GetXaxis().SetTitleSize(0.05)
    bb_bar_cs_NNPDF.GetXaxis().SetLabelSize(0.05)
    bb_bar_cs_NNPDF.GetYaxis().SetTitle("d#sigma_{b#bar{b}} / d#it{y} (#mub)")
    bb_bar_cs_NNPDF.GetYaxis().SetTitleOffset(1.2)
    bb_bar_cs_NNPDF.GetYaxis().SetTitleSize(0.055)
    bb_bar_cs_NNPDF.GetYaxis().SetLabelSize(0.05)

    bb_bar_cs_NNPDF.SetMaximum(1400)
    bb_bar_cs_NNPDF.SetMinimum(0.075)
    ROOT.gPad.Modified()
    ROOT.gPad.Update()

    Leged_bb_cs_NNPDF = ROOT.TLegend(0.18, 0.18, 0.85, 0.425, " ", "brNDC")
    Leged_bb_cs_NNPDF.SetNColumns(2)

    Leged_bb_cs_NNPDF.AddEntry(FONLL_bb_cs_NNPDF_minmaxerror, "FONLL unc. tot.", "F")
    Leged_bb_cs_NNPDF.AddEntry(FONLL_bb_cs_NNPDF_scalerror, "FONLL unc. scale", "F")
    Leged_bb_cs_NNPDF.AddEntry(FONLL_bb_cs_NNPDF_masserror, "FONLL unc. mass", "F")
    Leged_bb_cs_NNPDF.AddEntry(FONLL_bb_cs_NNPDF_pdferror, "FONLL unc. pdf", "F")
    # Leged_bb_cs_NNPDF.AddEntry(FONLL_bb_cs_NNPDF_minmaxerror_bis, "FONLL", "CF")
    Leged_bb_cs_NNPDF.AddEntry(bb_cs_PYTHIA_fwd, "#mu^{#plus}#mu^{#minus} #leftarrow b,b", "EP")
    Leged_bb_cs_NNPDF.AddEntry(bb_cs_PYTHIA_mid, "e^{#plus} e^{#minus} #leftarrow b,b", "EP")
    

    Leged_bb_cs_NNPDF.SetBorderSize(0)
    Leged_bb_cs_NNPDF.SetFillColor(10)
    Leged_bb_cs_NNPDF.SetFillStyle(1)
    Leged_bb_cs_NNPDF.SetLineStyle(0)
    Leged_bb_cs_NNPDF.SetLineColor(0)
    Leged_bb_cs_NNPDF.SetTextFont(42)
    Leged_bb_cs_NNPDF.SetTextSize(0.05)
    Leged_bb_cs_NNPDF.Draw()
    
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    letexTitle = ROOT.TLatex()
    letexTitle.SetTextSize(0.055)
    letexTitle.SetNDC()
    letexTitle.SetTextFont(42)
    letexTitle.DrawLatex(0.2, 0.82, "ALICE Preliminary, pp#sqrt{#it{s}} = 13 TeV")
    letexTitle.SetTextSize(0.0375)
    letexTitle.DrawLatex(0.2, 0.74, "FONLL NNPDF")
    ROOT.gPad.Modified()
    ROOT.gPad.Update()

    canvas_bb_cs_NNPDF.SaveAs("canvas_bb_cs_NNPDF.png")
    canvas_bb_cs_NNPDF.SaveAs("canvas_bb_cs_NNPDF.pdf")

    ##CPTQE PDF CHARM

    canvas_bb_cs_CTEQ6 = LoadStyle_canvas()
    canvas_bb_cs_CTEQ6.SetName("canvas_bb_cs_CTEQ6")
    canvas_bb_cs_CTEQ6.cd()
    FONLL_bb_cs_CTEQ6_minmaxerror = ROOT.TGraphAsymmErrors(17, Arr_y_FONLL_CTEQ6, Arr_ds_dy_bb_FONLL_CTEQ6, Arr_low_ey_FONLL_CTEQ6,
                                                           Arr_up_ey_FONLL_CTEQ6, Arr_min_ds_dy_bb_FONLL_CTEQ6, Arr_max_ds_dy_bb_FONLL_CTEQ6)
    FONLL_bb_cs_CTEQ6_minmaxerror.SetName("FONLL_bb_cs_CTEQ6_minmaxerror")
    FONLL_bb_cs_CTEQ6_minmaxerror.SetLineWidth(2)
    FONLL_bb_cs_CTEQ6_minmaxerror.SetFillColorAlpha(ROOT.kGray, 0.7)
    FONLL_bb_cs_CTEQ6_minmaxerror.SetFillColorAlpha(ROOT.kGray, 0.7)
    
    FONLL_bb_cs_CTEQ6_scalerror = ROOT.TGraphAsymmErrors(17, Arr_y_FONLL_CTEQ6, Arr_ds_dy_bb_FONLL_CTEQ6, Arr_low_ey_FONLL_CTEQ6,
                                                         Arr_up_ey_FONLL_CTEQ6, Arr_min_sc_ds_dy_bb_FONLL_CTEQ6, Arr_max_sc_ds_dy_bb_FONLL_CTEQ6)
    FONLL_bb_cs_CTEQ6_scalerror.SetName("FONLL_bb_cs_CTEQ6_scalerror")
    FONLL_bb_cs_CTEQ6_scalerror.SetFillColorAlpha(ROOT.kCyan-10, 0.2)
    FONLL_bb_cs_CTEQ6_scalerror.SetFillColorAlpha(ROOT.kCyan-10, 0.2)
    FONLL_bb_cs_CTEQ6_scalerror.SetLineWidth(2)
    

    FONLL_bb_cs_CTEQ6_masserror = ROOT.TGraphAsymmErrors(17, Arr_y_FONLL_CTEQ6, Arr_ds_dy_bb_FONLL_CTEQ6, Arr_low_ey_FONLL_CTEQ6,
                                                         Arr_up_ey_FONLL_CTEQ6, Arr_min_mass_ds_dy_bb_FONLL_CTEQ6, Arr_max_mass_ds_dy_bb_FONLL_CTEQ6)
    FONLL_bb_cs_CTEQ6_masserror.SetName("FONLL_bb_cs_CTEQ6_masserror")
    FONLL_bb_cs_CTEQ6_masserror.SetMarkerStyle(20)
    FONLL_bb_cs_CTEQ6_masserror.SetFillColorAlpha(ROOT.kGreen, 0.2)
    FONLL_bb_cs_CTEQ6_masserror.SetFillColorAlpha(ROOT.kGreen, 0.2)
    FONLL_bb_cs_CTEQ6_masserror.SetLineWidth(2)
    

    FONLL_bb_cs_CTEQ6_pdferror = ROOT.TGraphAsymmErrors(17, Arr_y_FONLL_CTEQ6, Arr_ds_dy_bb_FONLL_CTEQ6, Arr_low_ey_FONLL_CTEQ6,
                                                         Arr_up_ey_FONLL_CTEQ6, Arr_min_pdf_ds_dy_bb_FONLL_CTEQ6, Arr_max_pdf_ds_dy_bb_FONLL_CTEQ6)
    FONLL_bb_cs_CTEQ6_pdferror.SetName("FONLL_bb_cs_CTEQ6_pdferror")
    FONLL_bb_cs_CTEQ6_pdferror.SetMarkerStyle(20)
    FONLL_bb_cs_CTEQ6_pdferror.SetFillColorAlpha(ROOT.kMagenta-9, 0.2)
    FONLL_bb_cs_CTEQ6_pdferror.SetFillColorAlpha(ROOT.kMagenta-9, 0.2)
    FONLL_bb_cs_CTEQ6_pdferror.SetLineWidth(2)
    

    bb_bar_cs_CTEQ6 = ROOT.TMultiGraph()
    bb_bar_cs_CTEQ6. Add(FONLL_bb_cs_CTEQ6_scalerror)
    bb_bar_cs_CTEQ6. Add(FONLL_bb_cs_CTEQ6_masserror)
    bb_bar_cs_CTEQ6. Add(FONLL_bb_cs_CTEQ6_pdferror)
    bb_bar_cs_CTEQ6. Add(FONLL_bb_cs_CTEQ6_minmaxerror)
    ROOT.gPad.SetLogy()
    
    bb_bar_cs_CTEQ6.Add(bb_cs_PYTHIA_fwd, "APS; Z ; 5 s=0.5")
    bb_bar_cs_CTEQ6.Add(bb_cs_PYTHIA_mid, "APS; Z ; 5 s=0.5")
    bb_bar_cs_CTEQ6.Draw("A3")
    # FONLL_bb_cs_CTEQ6_minmaxerror.SetLineColorAlpha(ROOT.kGray, 1)
    # FONLL_bb_cs_CTEQ6_minmaxerror_bis.Draw("C PLC PFC SAME")
    bb_bar_cs_CTEQ6.GetXaxis().SetTitle("#it{y}")
    bb_bar_cs_CTEQ6.GetXaxis().SetTitleOffset(1.2)
    bb_bar_cs_CTEQ6.GetXaxis().SetTitleSize(0.05)
    bb_bar_cs_CTEQ6.GetXaxis().SetLabelSize(0.05)
    bb_bar_cs_CTEQ6.GetYaxis().SetTitle("d#sigma_{b#bar{b}} / d#it{y} (#mub)")
    bb_bar_cs_CTEQ6.GetYaxis().SetTitleOffset(1.2)
    bb_bar_cs_CTEQ6.GetYaxis().SetTitleSize(0.055)
    bb_bar_cs_CTEQ6.GetYaxis().SetLabelSize(0.05)

    bb_bar_cs_CTEQ6.SetMaximum(1400)
    bb_bar_cs_CTEQ6.SetMinimum(0.075)
    ROOT.gPad.Modified()
    ROOT.gPad.Update()

    Leged_bb_cs_CTEQ6 = ROOT.TLegend(0.18, 0.18, 0.85, 0.425, " ", "brNDC")
    Leged_bb_cs_CTEQ6.SetNColumns(2)

    Leged_bb_cs_CTEQ6.AddEntry(FONLL_bb_cs_CTEQ6_minmaxerror, "FONLL unc. tot.", "F")
    Leged_bb_cs_CTEQ6.AddEntry(FONLL_bb_cs_CTEQ6_scalerror, "FONLL unc. scale", "F")
    Leged_bb_cs_CTEQ6.AddEntry(FONLL_bb_cs_CTEQ6_masserror, "FONLL unc. mass", "F")
    Leged_bb_cs_CTEQ6.AddEntry(FONLL_bb_cs_CTEQ6_pdferror, "FONLL unc. pdf", "F")
    # Leged_bb_cs_CTEQ6.AddEntry(FONLL_bb_cs_CTEQ6_minmaxerror_bis, "FONLL", "CF")
    Leged_bb_cs_CTEQ6.AddEntry(bb_cs_PYTHIA_fwd, "#mu^{#plus}#mu^{#minus} #leftarrow b,b", "EP")
    Leged_bb_cs_CTEQ6.AddEntry(bb_cs_PYTHIA_mid, "e^{#plus} e^{#minus} #leftarrow b,b", "EP")
    

    Leged_bb_cs_CTEQ6.SetBorderSize(0)
    Leged_bb_cs_CTEQ6.SetFillColor(10)
    Leged_bb_cs_CTEQ6.SetFillStyle(1)
    Leged_bb_cs_CTEQ6.SetLineStyle(0)
    Leged_bb_cs_CTEQ6.SetLineColor(0)
    Leged_bb_cs_CTEQ6.SetTextFont(42)
    Leged_bb_cs_CTEQ6.SetTextSize(0.05)
    Leged_bb_cs_CTEQ6.Draw()
    
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    letexTitle = ROOT.TLatex()
    letexTitle.SetTextSize(0.055)
    letexTitle.SetNDC()
    letexTitle.SetTextFont(42)
    letexTitle.DrawLatex(0.2, 0.82, "ALICE Preliminary, pp#sqrt{#it{s}} = 13 TeV")
    letexTitle.SetTextSize(0.0375)
    letexTitle.DrawLatex(0.2, 0.74, "FONLL CTEQ6")
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    canvas_bb_cs_CTEQ6.SaveAs("canvas_bb_cs_CTEQ6.png")
    canvas_bb_cs_CTEQ6.SaveAs("canvas_bb_cs_CTEQ6.pdf")

    input()


main()
