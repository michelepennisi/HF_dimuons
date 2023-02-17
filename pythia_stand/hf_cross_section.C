
TCanvas *compare_cs(TH1D *MC, TH1D *data, TString part, TString rapidity);
TCanvas *compare_cs_LHCb_ratio(TH1D *MC[4][6], TH1D *data[4][6], TString part[4], TString rapidity[6], Int_t i_hadron = 0);
TCanvas *compare_cs_LHCb(TH1D *MC[4][6], TH1D *data[4][6], TString part[4], TString rapidity[6], Int_t i_hadron = 0);
void hf_cross_section()
{
    TString output_dir("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/hf_hadron_cs");

    TFile *fIn = new TFile("Hist_HFHadronpythia_sim_SoftQCD_Def_1000000_2710_DefaultBR_HFcount_HFhadron.root", "READ");
    fIn->ls();
    double pythia_cs = 5.642e+4; // pythia cs in mb

    const Int_t n_Hadron_studied = 4;
    const Int_t n_Rapidity_studied = 6;
    Double_t norm_factor[n_Hadron_studied];
    norm_factor[0] = pythia_cs / (2. * 1000000);
    norm_factor[1] = pythia_cs * 0.0913 / (2.);
    norm_factor[2] = pythia_cs * 0.0228 / (2.);
    norm_factor[3] = pythia_cs / (2. * 1000000);

    TString name_root_files[n_Hadron_studied + 1];
    name_root_files[0].Form("%s/root_files/D0-Lc-Xsec-pp13.root", output_dir.Data());
    name_root_files[1].Form("%s/root_files/HFPtSpectrum_13TeV_Dplus_withOldSys_RenuBala.root", output_dir.Data());
    name_root_files[2].Form("%s/root_files/HFPtSpectrum_13TeV_Ds_JulienHamon.root", output_dir.Data());
    name_root_files[3].Form("%s/root_files/D0-Lc-Xsec-pp13.root", output_dir.Data());
    name_root_files[4].Form("%s/root_files/HEPData-ins1396331-v2-root.root", output_dir.Data());

    TString name_histo[n_Hadron_studied + 1];
    name_histo[0].Form("hxsecD0");
    name_histo[1].Form("histoSigmaCorr");
    name_histo[2].Form("histoSigmaCorr");
    name_histo[3].Form("hxsecLc");
    // name_root_files[4].Form("HEPData-ins1396331-v2-root.root");

    TString name_Hadron_studied[n_Hadron_studied];
    TString name_rapidity_studied[n_Rapidity_studied];

    name_Hadron_studied[0].Form("Dzero");
    name_Hadron_studied[1].Form("Dplus");
    name_Hadron_studied[2].Form("Dstrange");
    name_Hadron_studied[3].Form("Lambda");

    name_rapidity_studied[0].Form("MidY");
    name_rapidity_studied[1].Form("FwdY1");
    name_rapidity_studied[2].Form("FwdY2");
    name_rapidity_studied[3].Form("FwdY3");
    name_rapidity_studied[4].Form("FwdY4");
    name_rapidity_studied[5].Form("FwdY5");

    TH1D *h_ptHadron_prompt[n_Hadron_studied][n_Rapidity_studied];
    TH1D *h_yHadron_prompt[n_Hadron_studied][n_Rapidity_studied];
    TH1D *h_pdgHadron_prompt[n_Hadron_studied][n_Rapidity_studied];

    TH1D *h_cross_section_Hadron_prompt[n_Hadron_studied][n_Rapidity_studied];

    TH1D *h_ptHadron_Notprompt[n_Hadron_studied][n_Rapidity_studied];
    TH1D *h_yHadron_Notprompt[n_Hadron_studied][n_Rapidity_studied];
    TH1D *h_pdgHadron_Notprompt[n_Hadron_studied][n_Rapidity_studied];

    TCanvas *cs_compare[n_Hadron_studied][n_Rapidity_studied];
    TH1D *hist_data_LHCb[n_Hadron_studied][n_Rapidity_studied];

    for (size_t i_hadron = 0; i_hadron < n_Hadron_studied; i_hadron++)
    {
        for (size_t i_rapidity = 0; i_rapidity < n_Rapidity_studied; i_rapidity++)
        {
            h_ptHadron_prompt[i_hadron][i_rapidity] = (TH1D *)fIn->Get(Form("%s/h_pt%s_prompt_%s", name_rapidity_studied[i_rapidity].Data(), name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
            h_yHadron_prompt[i_hadron][i_rapidity] = (TH1D *)fIn->Get(Form("%s/h_y%s_prompt_%s", name_rapidity_studied[i_rapidity].Data(), name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
            h_pdgHadron_prompt[i_hadron][i_rapidity] = (TH1D *)fIn->Get(Form("%s/h_pdg%s_prompt_%s", name_rapidity_studied[i_rapidity].Data(), name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));

            h_ptHadron_Notprompt[i_hadron][i_rapidity] = (TH1D *)fIn->Get(Form("%s/h_pt%s_Notprompt_%s", name_rapidity_studied[i_rapidity].Data(), name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
            h_yHadron_Notprompt[i_hadron][i_rapidity] = (TH1D *)fIn->Get(Form("%s/h_y%s_Notprompt_%s", name_rapidity_studied[i_rapidity].Data(), name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
            h_pdgHadron_Notprompt[i_hadron][i_rapidity] = (TH1D *)fIn->Get(Form("%s/h_pdg%s_Notprompt_%s", name_rapidity_studied[i_rapidity].Data(), name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
            // if (i_rapidity == 0)
            if (i_rapidity == 0)
            {
                TFile *fIn_data = new TFile(Form("%s", name_root_files[i_hadron].Data()), "READ");
                TH1D *hist_data = (TH1D *)fIn_data->Get(Form("%s", name_histo[i_hadron].Data()));
                hist_data->SetDirectory(0);
                h_ptHadron_prompt[i_hadron][i_rapidity]->Scale(norm_factor[i_hadron], "width");
                if (i_hadron == 1 || i_hadron == 2)
                {
                    h_ptHadron_prompt[i_hadron][i_rapidity]->Scale(1. / 1000000);
                    hist_data->Scale(1. / 1000000);
                }

                cs_compare[i_hadron][i_rapidity] = compare_cs(h_ptHadron_prompt[i_hadron][i_rapidity], hist_data, name_Hadron_studied[i_hadron], name_rapidity_studied[i_rapidity]);
                cs_compare[i_hadron][i_rapidity]->SetName(Form("cs_compariso_%s_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
            }
            else
            {
                if (i_hadron < n_Hadron_studied - 1)
                {
                    TFile *fIn_data = new TFile(Form("%s", name_root_files[n_Hadron_studied].Data()), "READ");
                    hist_data_LHCb[i_hadron][i_rapidity] = (TH1D *)fIn_data->Get(Form("Table %zu/Hist1D_y%zu", i_hadron + 1, i_rapidity));
                    TGraphAsymmErrors *graph_data = (TGraphAsymmErrors *)fIn_data->Get(Form("Table %zu/Graph1D_y%zu", i_hadron + 1, i_rapidity));
                    hist_data_LHCb[i_hadron][i_rapidity]->SetDirectory(0);
                    hist_data_LHCb[i_hadron][i_rapidity]->Reset();
                    for (size_t i_bin = 0; i_bin < hist_data_LHCb[i_hadron][i_rapidity]->GetXaxis()->GetNbins(); i_bin++)
                    {
                        hist_data_LHCb[i_hadron][i_rapidity]->SetBinContent(i_bin + 1, graph_data->Eval(hist_data_LHCb[i_hadron][i_rapidity]->GetBinCenter(i_bin + 1)));
                        hist_data_LHCb[i_hadron][i_rapidity]->SetBinError(i_bin + 1, graph_data->GetErrorY(i_bin + 1));
                    }
                    // graph_data->Draw();
                    // hist_data_LHCb[i_hadron][i_rapidity]->SetLineColor(kRed);
                    // hist_data_LHCb[i_hadron][i_rapidity]->Draw("PESAME");
                    // return;

                    h_ptHadron_prompt[i_hadron][i_rapidity]->Scale(2 * pythia_cs / (1000000), "width");

                    // cs_compare[i_hadron][i_rapidity] = compare_cs(h_ptHadron_prompt[i_hadron][i_rapidity], hist_data_LHCb[i_hadron][i_rapidity], name_Hadron_studied[i_hadron], name_rapidity_studied[i_rapidity]);
                    // cs_compare[i_hadron][i_rapidity]->SetName(Form("cs_compariso_%s_%s", name_Hadron_studied[i_hadron].Data(), name_rapidity_studied[i_rapidity].Data()));
                }
                else
                    continue;

                // return;s
            }
        }
    }
    TCanvas *cs_compare_LHCb_Dzero = compare_cs_LHCb_ratio(h_ptHadron_prompt, hist_data_LHCb, name_Hadron_studied, name_rapidity_studied, 0);
    cs_compare_LHCb_Dzero->SetName("cs_compare_LHCb_Dzero");
    cs_compare_LHCb_Dzero->SaveAs(Form("%s/plot/cs_compare_LHCb_Dzero_Monashstd.pdf", output_dir.Data()));
    cs_compare_LHCb_Dzero->SaveAs(Form("%s/plot/cs_compare_LHCb_Dzero_Monashstd.png", output_dir.Data()));

    TCanvas *cs_compare_LHCb_Dplus = compare_cs_LHCb_ratio(h_ptHadron_prompt, hist_data_LHCb, name_Hadron_studied, name_rapidity_studied, 1);
    cs_compare_LHCb_Dplus->SetName("cs_compare_LHCb_Dplus");
    cs_compare_LHCb_Dplus->SaveAs(Form("%s/plot/cs_compare_LHCb_Dplus_Monashstd.pdf", output_dir.Data()));
    cs_compare_LHCb_Dplus->SaveAs(Form("%s/plot/cs_compare_LHCb_Dplus_Monashstd.png", output_dir.Data()));

    TCanvas *cs_compare_LHCb_Dstrage = compare_cs_LHCb_ratio(h_ptHadron_prompt, hist_data_LHCb, name_Hadron_studied, name_rapidity_studied, 2);
    cs_compare_LHCb_Dstrage->SetName("cs_compare_LHCb_Dstrage");
    cs_compare_LHCb_Dstrage->SaveAs(Form("%s/plot/cs_compare_LHCb_Dstrage_Monashstd.pdf", output_dir.Data()));
    cs_compare_LHCb_Dstrage->SaveAs(Form("%s/plot/cs_compare_LHCb_Dstrage_Monashstd.png", output_dir.Data()));

    // h_ptHadron_prompt[0][0]->SetLineColor(kRed);
    // h_ptHadron_prompt[0][0]->Draw();

    // h_cross_section_Hadron_prompt[0][0]->Scale(1., "width");
    // h_cross_section_Hadron_prompt[0][0]->Draw("same");
    // TH1D *prova = (TH1D *)h_ptHadron_prompt[0][5]->Rebin(fPtBinning->GetNbins(), "prova", array);

    // prova->Draw("same");
}

TCanvas *compare_cs(TH1D *MC, TH1D *data, TString part, TString rapidity)
{
    gStyle->SetImageScaling(3.);
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1100);
    canvas->SetTicks();
    canvas->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.375, 1.0, 1.0);
    pad1->SetTicks();
    pad1->SetLogy(1);
    pad1->SetTopMargin(0.05);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.14);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.375);
    pad2->SetTicks();
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.03);
    pad2->SetLeftMargin(0.14);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();

    pad1->cd();
    // if (norm)
    // {
    //   frame->SetMaximum(1.5e-0);
    //   frame->SetMinimum(1.5e-12);
    // }
    // else
    // {
    //   frame->SetMaximum(1.5e+7);
    //   frame->SetMinimum(1.5e-2);
    // }
    TH2 *h_grid = new TH2D(Form("grid_%s", data->GetName()), " ", data->GetXaxis()->GetNbins(), data->GetXaxis()->GetBinLowEdge(1), data->GetXaxis()->GetBinLowEdge(data->GetXaxis()->GetNbins() + 1), 1000, data->GetMinimum() * 0.025, data->GetMaximum() * 2050);
    h_grid->GetYaxis()->SetTitle("d^{2}#sigma/d#it{p}_{T}d#it{y} #mub/(GeV/#it{c})");
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.0475);
    h_grid->GetXaxis()->SetLabelSize(0.045);

    h_grid->GetYaxis()->SetNdivisions(505);
    h_grid->GetYaxis()->SetTitleOffset(1.1);
    h_grid->GetYaxis()->SetTitleSize(0.05);
    h_grid->GetYaxis()->SetLabelSize(0.05);

    h_grid->Draw();

    data->SetLineColor(kBlack);
    data->SetLineWidth(2.0);
    data->SetMarkerColor(kBlack);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(1.2);

    data->Draw("PEsame");
    MC->SetLineColor(kRed);
    MC->SetMarkerColor(kRed);
    MC->SetMarkerStyle(24);
    MC->SetLineWidth(2.0);
    MC->SetMarkerSize(1.2);

    MC->Draw("PEsame");

    TLegend *legend = new TLegend(0.675, 0.375, 1.0, 0.595);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0425);
    // legend->SetHeader("Data");
    legend->SetTextAlign(12);
    legend->AddEntry(data, "Data", "LP");
    legend->AddEntry(MC, "Monash", "LP");

    legend->Draw();
    TLatex *letexTitle = new TLatex();
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);
    letexTitle->SetTextSize(0.0475);
    // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.175, 0.875, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    if (part.Contains("Dzero"))
    {
        letexTitle->DrawLatex(0.175, 0.775, "D^{0}#rightarrow K^{#minus}#pi^{#plus} #plus c.c");
    }
    else if (part.Contains("Dplus"))
    {
        letexTitle->DrawLatex(0.175, 0.775, "D^{#plus}#rightarrow K^{#plus}#pi^{#minus}#pi^{#plus} #plus c.c");
    }
    else if (part.Contains("Dstrange"))
    {
        letexTitle->DrawLatex(0.175, 0.775, "D^{#plus}_{s}#rightarrow K^{#plus}K^{#minus}#pi^{#plus} #plus c.c");
    }
    else if (part.Contains("Lambda"))
    {
        letexTitle->DrawLatex(0.175, 0.775, "#Lambda^{#plus}_{c}#rightarrow p K^{#minus}#pi^{#plus} #plus c.c");
    }

    if (rapidity.Contains("MidY"))
        letexTitle->DrawLatex(0.825, 0.875, "|#it{y}|<0.5");
    else if (rapidity.Contains("FwdY1"))
        letexTitle->DrawLatex(0.805, 0.875, "2.0< #it{y}<2.5");
    else if (rapidity.Contains("FwdY2"))
        letexTitle->DrawLatex(0.805, 0.875, "2.5< #it{y}<3.0");
    else if (rapidity.Contains("FwdY3"))
        letexTitle->DrawLatex(0.805, 0.875, "3.0< #it{y}<3.5");
    else if (rapidity.Contains("FwdY4"))
        letexTitle->DrawLatex(0.805, 0.875, "3.5< #it{y}<4.0");
    else if (rapidity.Contains("FwdY5"))
        letexTitle->DrawLatex(0.805, 0.875, "4.0< #it{y}<4.5");

    pad2->cd();
    pad2->SetTicks();
    TLine *l = new TLine(data->GetXaxis()->GetBinLowEdge(1), 1.0, data->GetXaxis()->GetBinLowEdge(data->GetXaxis()->GetNbins() + 1), 1.0);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->SetLineColor(kRed);

    TH1D *ratio = (TH1D *)data->Clone("ratio");
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetMarkerStyle(20);
    // ratioata->Rebin(15);
    ratio->Divide(MC);

    TH2D *h_grid_ratio = new TH2D(Form("grid_ratio_%s", data->GetName()), " ", data->GetXaxis()->GetNbins(), data->GetXaxis()->GetBinLowEdge(1), data->GetXaxis()->GetBinLowEdge(data->GetXaxis()->GetNbins() + 1), 1000, 0.002, ratio->GetMaximum() * 1.45);
    h_grid_ratio->SetTitle("");

    // TH1D *c_data = (TH1D *)data->Clone("c_data");

    h_grid_ratio->GetYaxis()->SetTitle(Form("#frac{Data}{Monash}"));
    h_grid_ratio->GetYaxis()->CenterTitle();
    h_grid_ratio->GetYaxis()->SetNdivisions(504);
    h_grid_ratio->GetYaxis()->SetTitleSize(0.07);
    // h_grid_ratio->GetYaxis()->SetTitleOffset(0.8);
    h_grid_ratio->GetYaxis()->SetLabelOffset(0.02);
    h_grid_ratio->GetYaxis()->SetLabelSize(0.09);

    h_grid_ratio->GetXaxis()->SetTitleSize(0.09);
    h_grid_ratio->GetXaxis()->SetTitleOffset(1.1);
    h_grid_ratio->GetXaxis()->SetLabelSize(0.09);
    h_grid_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h_grid_ratio->Draw();
    ratio->Draw("PESAME");
    l->Draw();
    // l1->Draw();
    // l2->Draw();

    return canvas;
}
TCanvas *compare_cs_LHCb_ratio(TH1D *MC[4][6], TH1D *data[4][6], TString part[4], TString rapidity[6], Int_t i_hadron = 0)
{
    gStyle->SetImageScaling(3.);
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 2200);
    canvas->SetTicks();

    canvas->Divide(3, 2, 0.0001, 0.00001, 19);
    TLatex *letexTitle = new TLatex();
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);
    letexTitle->SetTextSize(0.0475);
    for (Int_t i_rapidity = 1; i_rapidity < 6; i_rapidity++)
    {
        data[i_hadron][i_rapidity]->SetLineColor(kBlack);
        data[i_hadron][i_rapidity]->SetLineWidth(2.0);
        data[i_hadron][i_rapidity]->SetMarkerColor(kBlack);
        data[i_hadron][i_rapidity]->SetMarkerStyle(20);
        data[i_hadron][i_rapidity]->SetMarkerSize(1.2);
        if (rapidity[i_rapidity].Contains("FwdY1"))
            data[i_hadron][i_rapidity]->SetTitle(" 2.0< #it{y}<2.5");
        else if (rapidity[i_rapidity].Contains("FwdY2"))
            data[i_hadron][i_rapidity]->SetTitle(" 2.5< #it{y}<3.0");
        else if (rapidity[i_rapidity].Contains("FwdY3"))
            data[i_hadron][i_rapidity]->SetTitle(" 3.0< #it{y}<3.5");
        else if (rapidity[i_rapidity].Contains("FwdY4"))
            data[i_hadron][i_rapidity]->SetTitle(" 3.5< #it{y}<4.0");
        else if (rapidity[i_rapidity].Contains("FwdY5"))
            data[i_hadron][i_rapidity]->SetTitle(" 4.0< #it{y}<4.5");

        data[i_hadron][i_rapidity]->GetYaxis()->SetTitle("d^{2}#sigma/d#it{p}_{T}d#it{y} #mub/(GeV/#it{c})");
        data[i_hadron][i_rapidity]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        data[i_hadron][i_rapidity]->GetXaxis()->SetTitleOffset(1.5);
        data[i_hadron][i_rapidity]->GetXaxis()->SetTitleSize(0.035);
        data[i_hadron][i_rapidity]->GetXaxis()->SetLabelSize(0.035);

        data[i_hadron][i_rapidity]->GetYaxis()->SetNdivisions(505);
        data[i_hadron][i_rapidity]->GetYaxis()->SetTitleOffset(1.75);
        data[i_hadron][i_rapidity]->GetYaxis()->SetTitleSize(0.035);
        data[i_hadron][i_rapidity]->GetYaxis()->SetLabelSize(0.035);

        MC[i_hadron][i_rapidity]->SetLineColor(kRed);
        MC[i_hadron][i_rapidity]->SetMarkerColor(kRed);
        MC[i_hadron][i_rapidity]->SetMarkerStyle(24);
        MC[i_hadron][i_rapidity]->SetLineWidth(2.0);
        MC[i_hadron][i_rapidity]->SetMarkerSize(1.2);

        canvas->cd(i_rapidity);
        gPad->SetLogy();

        auto rp = new TRatioPlot(data[i_hadron][i_rapidity], MC[i_hadron][i_rapidity]);

        rp->SetLeftMargin(0.14);
        rp->SetUpTopMargin(0.05);
        rp->SetUpBottomMargin(0.00);
        rp->SetLowTopMargin(0.00);
        rp->SetRightMargin(0.03);
        rp->SetH1DrawOpt("PE");
        rp->SetH2DrawOpt("PE");
        rp->SetGraphDrawOpt("AB");
        rp->Draw();
        rp->GetLowerRefGraph()->SetMinimum(-0.2);
        rp->GetLowerRefGraph()->SetMaximum(2.5);
        rp->GetLowYaxis()->SetNdivisions(508);
        rp->GetLowerRefYaxis()->SetTitle("#frac{Data}{Monash}");
        rp->GetLowerRefYaxis()->CenterTitle();
        rp->GetLowerRefYaxis()->SetTitleOffset(1.75);
        canvas->Update();
    }
    canvas->cd(6);
    TLegend *legend = new TLegend(0.175, 0.375, 0.35, 0.595);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0425);
    // legend->SetHeader("Data");
    legend->SetTextAlign(12);
    legend->AddEntry(data[i_hadron][1], "Data", "LP");
    legend->AddEntry(MC[i_hadron][1], "PYTHIA8 Monash", "LP");
    legend->Draw();

    // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.175, 0.875, "LHCb, pp #sqrt{#it{s}} = 13 TeV");
    if (part[i_hadron].Contains("Dzero"))
    {
        letexTitle->DrawLatex(0.175, 0.775, "D^{0}#rightarrow K^{#minus}#pi^{#plus} #plus c.c");
    }
    else if (part[i_hadron].Contains("Dplus"))
    {
        letexTitle->DrawLatex(0.175, 0.775, "D^{#plus}#rightarrow K^{#plus}#pi^{#minus}#pi^{#plus} #plus c.c");
    }
    else if (part[i_hadron].Contains("Dstrange"))
    {
        letexTitle->DrawLatex(0.175, 0.775, "D^{#plus}_{s}#rightarrow K^{#plus}K^{#minus}#pi^{#plus} #plus c.c");
    }
    else if (part[i_hadron].Contains("Lambda"))
    {
        letexTitle->DrawLatex(0.175, 0.775, "#Lambda^{#plus}_{c}#rightarrow p K^{#minus}#pi^{#plus} #plus c.c");
    }
    return canvas;
}

TCanvas *compare_cs_LHCb(TH1D *MC[4][6], TH1D *data[4][6], TString part[4], TString rapidity[6], Int_t i_hadron = 0)
{
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1100);
    canvas->SetTicks();
    canvas->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.375, 1.0, 1.0);
    pad1->SetTicks();
    pad1->SetLogy(1);
    pad1->SetTopMargin(0.05);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.14);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();
    TH1D *ratio[4][6];
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.375);
    pad2->SetTicks();
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.03);
    pad2->SetLeftMargin(0.14);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();

    TH2 *h_grid = new TH2D(Form("grid_%s", data[i_hadron][1]->GetName()), " ", data[i_hadron][1]->GetXaxis()->GetNbins(), data[i_hadron][1]->GetXaxis()->GetBinLowEdge(1), data[i_hadron][1]->GetXaxis()->GetBinLowEdge(data[i_hadron][1]->GetXaxis()->GetNbins() + 1), 1000, data[i_hadron][5]->GetMinimum() * 0.0000025, data[i_hadron][1]->GetMaximum() * 2050);
    h_grid->GetYaxis()->SetTitle("d^{2}#sigma/d#it{p}_{T}d#it{y} #mub/(GeV/#it{c})");
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.0475);
    h_grid->GetXaxis()->SetLabelSize(0.045);

    h_grid->GetYaxis()->SetNdivisions(505);
    h_grid->GetYaxis()->SetTitleOffset(1.1);
    h_grid->GetYaxis()->SetTitleSize(0.05);
    h_grid->GetYaxis()->SetLabelSize(0.05);
    pad1->cd();
    h_grid->Draw();
    TLegend *legend = new TLegend(0.675, 0.375, 1.0, 0.595);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0425);
    // legend->SetHeader("Data");
    legend->SetTextAlign(12);
    legend->AddEntry(data[i_hadron][1], "Data", "LP");
    legend->AddEntry(MC[i_hadron][1], "Monash", "LP");

    for (Int_t i_rapidity = 1; i_rapidity < 6; i_rapidity++)
    {
        data[i_hadron][i_rapidity]->SetLineColor(kBlack);
        data[i_hadron][i_rapidity]->SetLineWidth(2.0);
        data[i_hadron][i_rapidity]->SetMarkerColor(kBlack);
        data[i_hadron][i_rapidity]->SetMarkerStyle(20);
        data[i_hadron][i_rapidity]->SetMarkerSize(1.2);

        data[i_hadron][i_rapidity]->Scale(TMath::Power(10., -(i_rapidity - 1)));
        data[i_hadron][i_rapidity]->Draw("PEsame");
        MC[i_hadron][i_rapidity]->SetLineColor(kRed);
        MC[i_hadron][i_rapidity]->SetMarkerColor(kRed);
        MC[i_hadron][i_rapidity]->SetMarkerStyle(24);
        MC[i_hadron][i_rapidity]->SetLineWidth(2.0);
        MC[i_hadron][i_rapidity]->SetMarkerSize(1.2);

        MC[i_hadron][i_rapidity]->Scale(TMath::Power(10., -(i_rapidity - 1)));
        MC[i_hadron][i_rapidity]->Draw("PEsame");

        TLatex *letexTitle = new TLatex();
        letexTitle->SetNDC();
        letexTitle->SetTextFont(42);
        letexTitle->SetTextSize(0.0475);
        // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
        letexTitle->DrawLatex(0.175, 0.875, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
        if (part[i_hadron].Contains("Dzero"))
        {
            letexTitle->DrawLatex(0.175, 0.775, "D^{0}#rightarrow K^{#minus}#pi^{#plus} #plus c.c");
        }
        else if (part[i_hadron].Contains("Dplus"))
        {
            letexTitle->DrawLatex(0.175, 0.775, "D^{#plus}#rightarrow K^{#plus}#pi^{#minus}#pi^{#plus} #plus c.c");
        }
        else if (part[i_hadron].Contains("Dstrange"))
        {
            letexTitle->DrawLatex(0.175, 0.775, "D^{#plus}_{s}#rightarrow K^{#plus}K^{#minus}#pi^{#plus} #plus c.c");
        }
        else if (part[i_hadron].Contains("Lambda"))
        {
            letexTitle->DrawLatex(0.175, 0.775, "#Lambda^{#plus}_{c}#rightarrow p K^{#minus}#pi^{#plus} #plus c.c");
        }

        if (rapidity[i_rapidity].Contains("MidY"))
            letexTitle->DrawLatex(0.825, 0.875, "|#it{y}|<0.5");
        else if (rapidity[i_rapidity].Contains("FwdY1"))
            letexTitle->DrawLatex(0.805, 0.875, "2.0< #it{y}<2.5");
        else if (rapidity[i_rapidity].Contains("FwdY2"))
            letexTitle->DrawLatex(0.805, 0.875, "2.5< #it{y}<3.0");
        else if (rapidity[i_rapidity].Contains("FwdY3"))
            letexTitle->DrawLatex(0.805, 0.875, "3.0< #it{y}<3.5");
        else if (rapidity[i_rapidity].Contains("FwdY4"))
            letexTitle->DrawLatex(0.805, 0.875, "3.5< #it{y}<4.0");
        else if (rapidity[i_rapidity].Contains("FwdY5"))
            letexTitle->DrawLatex(0.805, 0.875, "4.0< #it{y}<4.5");

        // TLine *l = new TLine(data->GetXaxis()->GetBinLowEdge(1), 1.0, data->GetXaxis()->GetBinLowEdge(data->GetXaxis()->GetNbins() + 1), 1.0);
        // l->SetLineWidth(3);
        // l->SetLineStyle(2);
        // l->SetLineColor(kRed);

        // l->Draw();
        // l1->Draw();
        // l2->Draw();
        ratio[i_hadron][i_rapidity] = (TH1D *)data[i_hadron][i_rapidity]->Clone(Form("ratio_y%d", i_rapidity));
        ratio[i_hadron][i_rapidity]->SetLineColor(kBlack);
        ratio[i_hadron][i_rapidity]->SetMarkerColor(kBlack);
        ratio[i_hadron][i_rapidity]->SetMarkerStyle(20);
        // ratioata->Rebin(15);
        ratio[i_hadron][i_rapidity]->Divide(MC[i_hadron][i_rapidity]);
        ratio[i_hadron][i_rapidity]->Scale(TMath::Power(10., -(i_rapidity - 1)));
    }
    legend->Draw();

    pad2->cd();
    TH2D *h_grid_ratio = new TH2D(Form("grid_ratio_%s", data[i_hadron][1]->GetName()), " ", data[i_hadron][1]->GetXaxis()->GetNbins(), data[i_hadron][1]->GetXaxis()->GetBinLowEdge(1), data[i_hadron][1]->GetXaxis()->GetBinLowEdge(data[i_hadron][1]->GetXaxis()->GetNbins() + 1), 10000, ratio[i_hadron][5]->GetMinimum() * 0.0000025, ratio[i_hadron][5]->GetMaximum() * 15);
    // TH2D *h_grid_ratio = new TH2D(Form("grid_ratio_%s", "ciao"), " ", 100, 0, 30, 1000, 0.00000000002, 1.45);
    h_grid_ratio->SetTitle("");

    // TH1D *c_data = (TH1D *)data->Clone("c_data");

    h_grid_ratio->GetYaxis()->SetTitle(Form("#frac{Data}{Monash}"));
    h_grid_ratio->GetYaxis()->CenterTitle();
    h_grid_ratio->GetYaxis()->SetNdivisions(504);
    h_grid_ratio->GetYaxis()->SetTitleSize(0.07);
    // h_grid_ratio->GetYaxis()->SetTitleOffset(0.8);
    h_grid_ratio->GetYaxis()->SetLabelOffset(0.02);
    h_grid_ratio->GetYaxis()->SetLabelSize(0.09);

    h_grid_ratio->GetXaxis()->SetTitleSize(0.09);
    h_grid_ratio->GetXaxis()->SetTitleOffset(1.1);
    h_grid_ratio->GetXaxis()->SetLabelSize(0.09);
    h_grid_ratio->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h_grid_ratio->Draw();

    for (Int_t i_rapidity = 1; i_rapidity < 6; i_rapidity++)
    {

        ratio[i_hadron][i_rapidity]->Draw("PESAME");
    }

    return canvas;
}