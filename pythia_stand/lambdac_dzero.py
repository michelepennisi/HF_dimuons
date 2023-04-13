import ROOT as rt

fIn_monash=rt.TFile("test/Hist_monash_sim_fixed_10Mev.root","READ")
fIn_mode2=rt.TFile("test/Hist_mode2_sim_fixed_10Mev.root","READ")

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
    # pad1.SetLogy(1)
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


    # h_grid->Draw();

    mode2.GetYaxis().SetNdivisions(505)
    mode2.GetYaxis().SetTitleOffset(1.1)
    mode2.GetYaxis().SetTitleSize(0.045)
    mode2.GetYaxis().SetLabelSize(0.0425)
    mode2.Draw("PE")
    monash.Draw("PESAME")

    
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
    ratio.GetYaxis().SetTitleOffset(0.75)
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

def lambdac_dzero():
    LoadStyle()
    rebin=20
    
    lambdac_monash_midy=fIn_monash.Get("hf_hadron/prompt/h_ptLambda_prompt_MidY")
    lambdac_monash_midy.SetName("lambdac_monash_midy")
    lambdac_monash_midy.Rebin(rebin)
    lambdac_monash_midy.Scale(1.,"width")
    lambdac_monash_midy.Sumw2()
    
    dzero_monash_midy=fIn_monash.Get("hf_hadron/prompt/h_ptDzero_prompt_MidY")
    dzero_monash_midy.SetName("dzero_monash_midy")
    dzero_monash_midy.Rebin(rebin)
    dzero_monash_midy.Scale(1.,"width")
    dzero_monash_midy.Sumw2()

    lambdac_over_dzero_monash_midy=lambdac_monash_midy.Clone("lambdac_over_dzero")
    lambdac_over_dzero_monash_midy.Divide(dzero_monash_midy)
    # lambdac_over_dzero_monash.Draw()
    lambdac_over_dzero_monash_midy.Sumw2()

    lambdac_mode2_midy=fIn_mode2.Get("hf_hadron/prompt/h_ptLambda_prompt_MidY")
    lambdac_mode2_midy.SetName("lambdac_mode2_midy")
    lambdac_mode2_midy.Rebin(rebin)
    lambdac_mode2_midy.Scale(1.,"width")
    lambdac_mode2_midy.Sumw2()
    
    dzero_mode2_midy=fIn_mode2.Get("hf_hadron/prompt/h_ptDzero_prompt_MidY")
    dzero_mode2_midy.SetName("dzero_mode2_midy")
    dzero_mode2_midy.Rebin(rebin)
    dzero_mode2_midy.Scale(1.,"width")
    dzero_mode2_midy.Sumw2()

    lambdac_over_dzero_mode2_midy=lambdac_mode2_midy.Clone("lambdac_over_dzero")
    lambdac_over_dzero_mode2_midy.Divide(dzero_mode2_midy)
    # lambdac_over_dzero_mode2.Draw()
    lambdac_over_dzero_mode2_midy.Sumw2()

    lambdac_over_dzero_mode2_over_monash_midy=lambdac_over_dzero_mode2_midy.Clone("lambdac_over_dzero_mode2_over_monash")
    lambdac_over_dzero_mode2_over_monash_midy.Divide(lambdac_over_dzero_monash_midy)

    lambdac_over_dzero_mode2_midy.SetMarkerColor(rt.kRed)
    lambdac_over_dzero_mode2_midy.SetMarkerStyle(24)
    lambdac_over_dzero_mode2_midy.SetMarkerSize(1.3)
    lambdac_over_dzero_mode2_midy.SetLineWidth(2)
    lambdac_over_dzero_mode2_midy.SetLineColor(rt.kRed)
    
    lambdac_over_dzero_monash_midy.SetMarkerColor(rt.kBlack)
    lambdac_over_dzero_monash_midy.SetMarkerStyle(20)
    lambdac_over_dzero_monash_midy.SetMarkerSize(1.3)
    lambdac_over_dzero_monash_midy.SetLineWidth(2)
    lambdac_over_dzero_monash_midy.SetLineColor(rt.kBlack)

    x_axis_name="#it{p}_{T} (GeV/#it{c})"
    y_axis_name="d#it{N}/d#it{p}_{T} (GeV/#it{c})^-1"
    
    lambdac_over_dzero_mode2_midy.GetXaxis().SetTitle(x_axis_name)
    lambdac_over_dzero_mode2_midy.GetYaxis().SetTitle(y_axis_name)

    lambdac_over_dzero_mode2_over_monash_midy.GetXaxis().SetTitle(x_axis_name)
    lambdac_over_dzero_mode2_over_monash_midy.GetYaxis().SetTitle("#frac{CR mode2}{Monash}")
    
    

    legend=rt.TLegend(0.6, 0.675, 0.9, 0.895)
    legend.AddEntry(lambdac_over_dzero_mode2_midy,"#Lambda_{c}/D^{0} with CR mode2","LP")
    legend.AddEntry(lambdac_over_dzero_monash_midy,"#Lambda_{c}/D^{0} with Monash","LP")

    canvas=ratio(lambdac_over_dzero_mode2_midy,lambdac_over_dzero_monash_midy,lambdac_over_dzero_mode2_over_monash_midy,"|y|<0.5",legend)
    canvas.SaveAs("lambdac_over_dzero_mode2_over_monash.pdf")
    canvas.SaveAs("lambdac_over_dzero_mode2_over_monash.png")


    lambdac_monash_fwdy2=fIn_monash.Get("hf_hadron/prompt/h_ptLambda_prompt_FwdY2")
    lambdac_monash_fwdy2.SetName("lambdac_monash_FwdY2")

    lambdac_monash_fwdy3=fIn_monash.Get("hf_hadron/prompt/h_ptLambda_prompt_FwdY3")
    lambdac_monash_fwdy3.SetName("lambdac_monash_FwdY3")

    lambdac_monash_fwdy4=fIn_monash.Get("hf_hadron/prompt/h_ptLambda_prompt_FwdY4")
    lambdac_monash_fwdy4.SetName("lambdac_monash_FwdY4")

    dzero_monash_fwdy2=fIn_monash.Get("hf_hadron/prompt/h_ptDzero_prompt_FwdY2")
    dzero_monash_fwdy2.SetName("dzero_monash_FwdY2")

    dzero_monash_fwdy3=fIn_monash.Get("hf_hadron/prompt/h_ptDzero_prompt_FwdY3")
    dzero_monash_fwdy3.SetName("dzero_monash_FwdY3")

    dzero_monash_fwdy4=fIn_monash.Get("hf_hadron/prompt/h_ptDzero_prompt_FwdY4")
    dzero_monash_fwdy4.SetName("dzero_monash_FwdY4")

    lambdac_monash_fwdy=lambdac_monash_fwdy2.Clone("lambdac_monash_fwdy")
    lambdac_monash_fwdy.Add(lambdac_monash_fwdy3)
    lambdac_monash_fwdy.Add(lambdac_monash_fwdy4)
    lambdac_monash_fwdy.Rebin(rebin)
    lambdac_monash_fwdy.Scale(1.,"width")
    lambdac_monash_fwdy.Sumw2()

    dzero_monash_fwdy=dzero_monash_fwdy2.Clone("dzero_monash_fwdy")
    dzero_monash_fwdy.Add(dzero_monash_fwdy3)
    dzero_monash_fwdy.Add(dzero_monash_fwdy4)
    dzero_monash_fwdy=fIn_monash.Get("hf_hadron/prompt/h_ptDzero_prompt_MidY")
    dzero_monash_fwdy.SetName("dzero_monash_midy")
    dzero_monash_fwdy.Rebin(rebin)
    dzero_monash_fwdy.Scale(1.,"width")
    dzero_monash_fwdy.Sumw2()

    input()

    
lambdac_dzero()
