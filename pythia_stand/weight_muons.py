import numpy as np

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

def ratio(mode2,monash,ratio,string,legend):
    canvas = rt.TCanvas("canvas", "canvas", 1000, 1100)
    canvas.SetTicks()
    canvas.cd()
    pad1 =rt.TPad("pad1", "pad1", 0, 0.375, 1.0, 1.0)
    pad1.SetTicks()
    pad1.SetLogy(1)
    pad1.SetTopMargin(0.05)
    pad1.SetRightMargin(0.03)
    pad1.SetLeftMargin(0.14)
    pad1.SetBottomMargin(0.0)
    pad1.Draw()

    pad2 = rt.TPad("pad2", "pad2", 0, 0.0, 1, 0.375)
    pad2.SetTicks()
    pad2.SetTopMargin(0.0)
    pad2.SetRightMargin(0.03)
    pad2.SetLeftMargin(0.14)
    pad2.SetBottomMargin(0.25)
    pad2.Draw()

    pad1.cd()
    # TH2 *h_grid = new TH2D(Form("grid_%s", data->GetName()), " ", data->GetXaxis()->GetNbins(), data->GetXaxis()->GetBinLowEdge(1), data->GetXaxis()->GetBinLowEdge(data->GetXaxis()->GetNbins() + 1), 1000, data->GetMinimum() * 0.025, data->GetMaximum() * 2050);
    # h_grid->GetYaxis()->SetTitle("d^{2}#sigma/d#it{p}_{T}d#it{y} #mub/(GeV/#it{c})");
    # h_grid->GetXaxis()->SetTitleOffset(1.3);
    # h_grid->GetXaxis()->SetTitleSize(0.0475);
    # h_grid->GetXaxis()->SetLabelSize(0.045);

    # h_grid->GetYaxis()->SetNdivisions(505);
    # h_grid->GetYaxis()->SetTitleOffset(1.1);
    # h_grid->GetYaxis()->SetTitleSize(0.05);
    # h_grid->GetYaxis()->SetLabelSize(0.05);

    # h_grid->Draw();

    monash.Draw("PE")
    mode2.Draw("PEsame")

    
    # // legend->SetNColumns(2);
    # // legend->AddEntry((TObject*)0, "", "");
    legend.SetFillStyle(0)
    legend.SetLineColor(rt.kWhite)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.0425)
    # // legend->SetHeader("Data");
    legend.SetTextAlign(12)
    # if (MC_name.Contains("Def"))
    # {
    #     legend->AddEntry(MC, "Monash", "LP");
    # }
    # else if (MC_name.Contains("Config+Atlas"))
    # {
    #     legend->AddEntry(MC, "MNR", "LP");
    # }
    # else if (MC_name.Contains("Mode2"))
    # {
    #     legend->AddEntry(MC, "CR mode 2", "LP");
    # }

    legend.Draw("SAME")
    
    letexTitle = rt.TLatex()
    letexTitle.SetNDC()
    letexTitle.SetTextFont(42)
    letexTitle.SetTextSize(0.0475)
    # // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle.DrawLatex(0.275, 0.875, "PYTHIA8, pp #sqrt{#it{s}} = 13 TeV")
    letexTitle.DrawLatex(0.275, 0.805, string)
    pad2.cd()
    pad2.SetTicks()
    l = rt.TLine(monash.GetXaxis().GetBinLowEdge(1), 1.0, monash.GetXaxis().GetBinLowEdge(monash.GetXaxis().GetNbins() + 1), 1.0)
    l.SetLineWidth(3)
    l.SetLineStyle(2)
    l.SetLineColor(rt.kRed)

    ratio.SetLineColor(rt.kBlack)
    ratio.SetMarkerColor(rt.kBlack)
    ratio.SetMarkerStyle(20)

    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetNdivisions(504)
    ratio.GetYaxis().SetTitleSize(0.075)
    ratio.GetYaxis().SetTitleOffset(0.9)
    ratio.GetYaxis().SetLabelOffset(0.02)
    ratio.GetYaxis().SetLabelSize(0.07)
    
    ratio.GetXaxis().SetTitleSize(0.075)
    ratio.GetXaxis().SetTitleOffset(1.3)
    ratio.GetXaxis().SetLabelSize(0.07)
    # h_grid_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    # h_grid_ratio->Draw();
    ratio.Draw()
    l.Draw("SAME")
    return canvas

def weight_muons(string):
    LoadStyle()
    fIn_monash=rt.TFile("test/Hist_monash_sim_fixed_10Mev.root")
    # fIn_monash=rt.TFile("test/Hist_monash_sim_new_fixed.root")
    # fIn_monash=rt.TFile("test/Hist_monash_100Mev.root")
    h_PtY_monash=fIn_monash.Get("muon/All/h_PtYMu_from%s_All" % string)
    
    h_PtY_monash.SetName(("h_pt_y_%s_monash" % string))
    h_Pt_monash=h_PtY_monash.ProjectionX()
    h_Pt_monash.SetName(("h_pt_%s_monash" % string))
    h_Y_monash=h_PtY_monash.ProjectionY()
    h_Y_monash.SetName(("h_y_%s_monash" % string))


    fIn_mode2=rt.TFile("test/Hist_mode2_sim_fixed_10Mev.root")
    # fIn_mode2=rt.TFile("test/Hist_mode2_sim_new_fixed.root")
    # fIn_mode2=rt.TFile("test/Hist_mode2_100Mev.root")
    h_PtY_mode2=fIn_mode2.Get("muon/All/h_PtYMu_from%s_All" % string)
    
    h_PtY_mode2.SetName(("h_pt_y_%s_mode2" % string))
    h_Pt_mode2=h_PtY_mode2.ProjectionX()
    h_Pt_mode2.SetName(("h_pt_%s_mode2" % string))
    h_Y_mode2=h_PtY_mode2.ProjectionY()
    h_Y_mode2.SetName(("h_y_%s_mode2" % string))

    xbins = np.array((0.0,0.5,1.0,2.0,3.0,5,7,9,12,15,20,30.0))
    # print(xbins.size)

    h_Pt_mode2_rebinned=h_Pt_mode2.Rebin(xbins.size - 1 ,("h_Pt_%s_mode2_rebinned" % string),xbins)
    # h_Pt_mode2_rebinned=h_Pt_mode2.Rebin(10)
    # h_Pt_mode2_rebinned.Scale(1.,"width")

    h_Pt_mode2_rebinned.Scale(1./h_Pt_mode2.Integral(),"width")
    h_Pt_mode2_rebinned.Sumw2()
    h_Pt_mode2_rebinned.SetMarkerColor(rt.kRed)
    h_Pt_mode2_rebinned.SetMarkerSize(1.2)
    h_Pt_mode2_rebinned.SetMarkerStyle(20)
    h_Pt_mode2_rebinned.SetLineWidth(2)
    h_Pt_mode2_rebinned.SetLineColor(rt.kRed)
    
    h_Pt_monash_rebinned=h_Pt_monash.Rebin(xbins.size - 1 ,("h_Pt_%s_monash_rebinned" % string),xbins)
    # h_Pt_monash_rebinned=h_Pt_monash.Rebin(10)
    # h_Pt_monash_rebinned.Scale(1.,"width")

    h_Pt_monash_rebinned.Scale(1./h_Pt_monash.Integral(),"width")
    h_Pt_monash_rebinned.Sumw2()
    h_Pt_monash_rebinned.SetMarkerColor(rt.kBlack)
    h_Pt_monash_rebinned.SetMarkerSize(1.2)
    h_Pt_monash_rebinned.SetMarkerStyle(20)
    h_Pt_monash_rebinned.SetLineWidth(2)
    h_Pt_monash_rebinned.SetLineColor(rt.kBlack)

    h_Pt_weight=h_Pt_mode2_rebinned.Clone(("h_Pt_%s_weight" % string))
    h_Pt_weight.Divide(h_Pt_monash_rebinned)
    h_Pt_weight.GetYaxis().SetTitle("#frac{mode2}{monash}")
    h_Pt_weight.SetMaximum(1.75)
    h_Pt_weight.SetMinimum(-0.25)
    
    if(string.find("Charm")!=-1):
        info="#mu #leftarrow c, -4.0 < #it{y} < -2.5"
    elif(string.find("Beauty")!=-1):
        info="#mu #leftarrow b, -4.0 < #it{y} < -2.5"
    
    legend=rt.TLegend(0.6, 0.275, 0.9, 0.495)
    legend.AddEntry(h_Pt_mode2_rebinned,"Mode2 Tune","LP")
    legend.AddEntry(h_Pt_monash_rebinned,"Monash Tune","LP")

    canvas_pt=ratio(h_Pt_mode2_rebinned,h_Pt_monash_rebinned,h_Pt_weight,info,legend)
    if(string.find("Charm")!=-1):
        canvas_pt.SetName("canvas_pt_charm")
    elif(string.find("Beauty")!=-1):
        canvas_pt.SetName("canvas_pt_beauty")

    canvas_pt.SaveAs("plot/%s.pdf" % canvas_pt.GetName())
    canvas_pt.SaveAs("plot/%s.png" % canvas_pt.GetName())

    

    # canvas_y=rt.TCanvas("canvas_y","canvas_y",1600,1200)
    # canvas_y.cd()
    h_Y_mode2.Rebin(10)
    h_Y_mode2.Scale(1./h_Y_mode2.Integral(),"width")
    h_Y_mode2.SetMarkerColor(rt.kRed)
    h_Y_mode2.SetMarkerSize(1.2)
    h_Y_mode2.SetMarkerStyle(20)
    h_Y_mode2.SetLineWidth(2)
    h_Y_mode2.SetLineColor(rt.kRed)
    
    h_Y_monash.Rebin(10)
    h_Y_monash.Scale(1./h_Y_monash.Integral(),"width")
    h_Y_monash.SetMarkerColor(rt.kBlack)
    h_Y_monash.SetMarkerSize(1.2)
    h_Y_monash.SetMarkerStyle(20)
    h_Y_monash.SetLineWidth(2)
    h_Y_monash.SetLineColor(rt.kBlack)

    h_Y_weight=h_Y_mode2.Clone("h_Y_%s_weight" % string)    
    h_Y_weight.Divide(h_Y_monash)
    h_Y_weight.GetYaxis().SetTitle("#frac{mode2}{monash}")
    h_Y_weight.SetMaximum(1.75)
    h_Y_weight.SetMinimum(0.25)

    if(string.find("Charm")!=-1):
        info="#mu #leftarrow c"
    elif(string.find("Beauty")!=-1):
        info="#mu #leftarrow b"
    
    canvas_y=ratio(h_Y_monash,h_Y_mode2,h_Y_weight,info,legend)
    if(string.find("Charm")!=-1):
        canvas_y.SetName("canvas_y_charm")
    elif(string.find("Beauty")!=-1):
        canvas_y.SetName("canvas_y_beauty")
    canvas_y.SaveAs("plot/%s.pdf" % canvas_y.GetName())
    canvas_y.SaveAs("plot/%s.png" % canvas_y.GetName())
    input()

    # fOut=rt.TFile("test/muon_mode2_correction.root","UPDATE")
    # fOut.cd()
    # h_Pt_monash_rebinned.Write()
    # h_Pt_mode2_rebinned.Write()
    # h_Pt_weight.Write()
    
    # h_Y_monash.Write()
    # h_Y_mode2.Write()
    # h_Y_weight.Write()

    # fOut.Close()
    
    

weight_muons("Charm")
weight_muons("Beauty")


