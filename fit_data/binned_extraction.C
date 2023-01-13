double FuncPtMass(double *x, double *par);
TCanvas *printMC_ratio(TH1D *data, TF1 *pdf, TH1D *hint, Color_t color, Int_t minx = 0, Int_t max_x = 30);

TCanvas *linear_error(Int_t i, TH1D *up_hint_error_cl95, TH1D *lo_hint_error_cl95, TH1D *up_hint_error_cl99, TH1D *lo_hint_error_cl99, TH1D *hint_cl95, TH1D *hint_cl99, Color_t color, Int_t minx = 0, Int_t max_x = 30);
void cl_95();

const Int_t n_DiMuSelection = 3;
TString name_DiMuSelection[n_DiMuSelection];

void binned_extraction()
{
    cl_95();
}
void cl_95()
{
    TFile *fOut = new TFile("cl_95.root", "UPDATE");
    TFile *fIn = new TFile("/home/michele_pennisi/dimuon_HF_pp/data/LHC18p/Hist_MC/3_11_22/HF/HistLite_HF_MCDimuHFTree_merged.root", "READ");

    name_DiMuSelection[0].Form("Charm");
    name_DiMuSelection[1].Form("Beauty");
    name_DiMuSelection[2].Form("Mixed");

    TH2D *h_PtM_MC[n_DiMuSelection];

    TH1D *h_Pt_MC[n_DiMuSelection];
    TH1D *h_M_MC[n_DiMuSelection];

    TF1 *pdf_M[n_DiMuSelection];
    TF1 *pdf_Pt[n_DiMuSelection];

    TH1D *mass_hint_cl95[n_DiMuSelection];

    TH1D *pt_hint_cl95[n_DiMuSelection];

    TH1D *mass_hint_cl99[n_DiMuSelection];
    
    TH1D *pt_hint_cl99[n_DiMuSelection];

    TH1D *up_mass_hint_error_cl95[n_DiMuSelection];

    TH1D *lo_mass_hint_error_cl95[n_DiMuSelection];

    TH1D *up_mass_hint_error_cl99[n_DiMuSelection];

    TH1D *lo_mass_hint_error_cl99[n_DiMuSelection];

    TH1D *up_pt_hint_error_cl95[n_DiMuSelection];

    TH1D *lo_pt_hint_error_cl95[n_DiMuSelection];

    TH1D *up_pt_hint_error_cl99[n_DiMuSelection];

    TH1D *lo_pt_hint_error_cl99[n_DiMuSelection];

    Color_t color[n_DiMuSelection+1] = {kMagenta + 2, kSpring - 6, kAzure + 9,kMagenta + 2};
    Color_t fillcolor[n_DiMuSelection] = {kMagenta - 10, kGreen - 10, kCyan - 10};

    TCanvas *mass_canvas[n_DiMuSelection];
    TCanvas *pt_canvas[n_DiMuSelection];

    TCanvas *mass_canvas_error[n_DiMuSelection];

    TCanvas *pt_canvas_error[n_DiMuSelection];

    Double_t B[n_DiMuSelection] = {2.5, 2.65, 3.43};
    Double_t n1[n_DiMuSelection] = {2.81, 2.81, 1.00};
    Double_t n2[n_DiMuSelection] = {2.5, 2.5, 7.5};

    for (Int_t i = 0; i < n_DiMuSelection; i++)
    {
        h_PtM_MC[i] = (TH2D *)fIn->Get(Form("DiMuon/M4/Rec_DQ_cut_match_LT/ULS/h_PtMDiMu_M4_Rec_DQ_cut_match_LT_ULS_from%s", name_DiMuSelection[i].Data()));

        h_Pt_MC[i] = (TH1D *)h_PtM_MC[i]->ProjectionX();
        h_Pt_MC[i]->Rebin(10);
        h_Pt_MC[i]->Scale(1., "width");
        h_Pt_MC[i]->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");

        pdf_Pt[i] = new TF1(Form("pt_pdf_%s", name_DiMuSelection[i].Data()), FuncPtMass, 0, 30, 4);
        pdf_Pt[i]->SetParameter(3, 20000);
        pdf_Pt[i]->SetParameter(0, B[i]);
        pdf_Pt[i]->SetParameter(1, n1[i]);
        pdf_Pt[i]->SetParameter(2, n2[i]);

        h_Pt_MC[i]->Fit(pdf_Pt[i], "LR0I");

        // mass_hint_cpti] = new TH1D(Form("pt_hint_%s", name_DiMuSelection[i].Data()), "Fitted Gaussian with .95 conf.band", 30, 0, 30);

        pt_hint_cl95[i] = (TH1D *)h_Pt_MC[i]->Clone(Form("pt_hint_%s_cl95", name_DiMuSelection[i].Data()));
        pt_hint_cl95[i]->SetTitle(Form("pt_hint_%s_cl95", name_DiMuSelection[i].Data()));

        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(pt_hint_cl95[i], 0.95);

        pt_hint_cl99[i] = (TH1D *)h_Pt_MC[i]->Clone(Form("pt_hint_%s_cl99", name_DiMuSelection[i].Data()));
        pt_hint_cl99[i]->SetTitle(Form("pt_hint_%s_cl99", name_DiMuSelection[i].Data()));

        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(pt_hint_cl99[i], 0.99);

        pt_hint_cl95[i]->SetStats(false);
        pt_hint_cl95[i]->SetFillColor(fillcolor[i]);

        // pt_canvas[i] = new TCanvas(Form("pt_canvas_%s", name_DiMuSelection[i].Data()), Form("pt_canvas_%s", name_DiMuSelection[i].Data()), 1000, 1200);
        pt_canvas[i] = printMC_ratio(h_Pt_MC[i], pdf_Pt[i], pt_hint_cl95[i], color[i], 0, 30);
        pt_canvas[i]->SetName(Form("pt_canvas_%s", name_DiMuSelection[i].Data()));
        pt_canvas[i]->SetTitle(Form("pt_canvas_%s", name_DiMuSelection[i].Data()));
        pt_canvas[i]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/%s.png", pt_canvas[i]->GetName()));
        up_pt_hint_error_cl95[i] = (TH1D *)pt_hint_cl95[i]->Clone(Form("up_pt_hint_error_%s_cl95", name_DiMuSelection[i].Data()));
        up_pt_hint_error_cl95[i]->Reset();
        lo_pt_hint_error_cl95[i] = (TH1D *)pt_hint_cl95[i]->Clone(Form("lo_pt_hint_error_%s_cl95", name_DiMuSelection[i].Data()));
        lo_pt_hint_error_cl95[i]->Reset();

        up_pt_hint_error_cl99[i] = (TH1D *)pt_hint_cl99[i]->Clone(Form("up_pt_hint_error_%s_cl99", name_DiMuSelection[i].Data()));
        up_pt_hint_error_cl99[i]->Reset();
        lo_pt_hint_error_cl99[i] = (TH1D *)pt_hint_cl99[i]->Clone(Form("lo_pt_hint_error_%s_cl99", name_DiMuSelection[i].Data()));
        lo_pt_hint_error_cl99[i]->Reset();

        pt_canvas_error[i] = linear_error(i, up_pt_hint_error_cl95[i], lo_pt_hint_error_cl95[i], up_pt_hint_error_cl99[i], lo_pt_hint_error_cl99[i], pt_hint_cl95[i], pt_hint_cl99[i], color[i], 0, 30);
        pt_canvas_error[i]->SetName(Form("pt_canvas_error_%s", name_DiMuSelection[i].Data()));
        pt_canvas_error[i]->SetTitle(Form("pt_canvas_error_%s", name_DiMuSelection[i].Data()));
        pt_canvas_error[i]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/%s.png", pt_canvas_error[i]->GetName()));
        // pt_canvas[i]->cd();
        // gPad->SetLogy();
        // h_Pt_MC[i]->Draw("PE");
        // pdf_Pt[i]->Draw("LSAME");
        // mass_hint_cl95[i]->Draw("e3 same");
        // h_Pt_MC[i]->Draw("PESAME");
        // pdf_Pt[i]->Draw("LSAME");

        h_M_MC[i] = (TH1D *)h_PtM_MC[i]->ProjectionY();
        h_M_MC[i]->Rebin(10);
        h_M_MC[i]->Scale(1., "width");
        h_M_MC[i]->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu #mu} (GeV/#it{c}^2)^{-1}");

        pdf_M[i] = new TF1(Form("mass_pdf_%s", name_DiMuSelection[i].Data()), FuncPtMass, 4, 30, 4);
        pdf_M[i]->SetParameter(3, 20000);
        pdf_M[i]->SetParameter(0, B[i]);
        pdf_M[i]->SetParameter(1, n1[i]);
        pdf_M[i]->SetParameter(2, n2[i]);

        h_M_MC[i]->Fit(pdf_M[i], "LR0I");

        mass_hint_cl95[i] = (TH1D *)h_M_MC[i]->Clone(Form("mass_hint_%s_cl95", name_DiMuSelection[i].Data()));
        mass_hint_cl95[i]->SetTitle(Form("mass_hint_%s_cl95", name_DiMuSelection[i].Data()));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mass_hint_cl95[i], 0.95);

        mass_hint_cl99[i] = (TH1D *)h_M_MC[i]->Clone(Form("mass_hint_%s_cl99", name_DiMuSelection[i].Data()));
        mass_hint_cl99[i]->SetTitle(Form("mass_hint_%s_cl99", name_DiMuSelection[i].Data()));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mass_hint_cl99[i], 0.99);

        mass_hint_cl95[i]->SetStats(false);
        mass_hint_cl95[i]->SetFillColor(fillcolor[i]);

        // mass_canvas[i] = new TCanvas(Form("mass_canvas_%s", name_DiMuSelection[i].Data()), Form("mass_canvas_%s", name_DiMuSelection[i].Data()), 1000, 1200);

        // mass_canvas[i]->cd();
        mass_canvas[i] = printMC_ratio(h_M_MC[i], pdf_M[i], mass_hint_cl95[i], color[i], 4, 30);
        mass_canvas[i]->SetName(Form("mass_canvas_%s", name_DiMuSelection[i].Data()));
        mass_canvas[i]->SetTitle(Form("mass_canvas_%s", name_DiMuSelection[i].Data()));
        mass_canvas[i]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/%s.png", mass_canvas[i]->GetName()));
        // gPad->SetLogy();
        // h_M_MC[i]->Draw("PE");
        // pdf_M[i]->Draw("LSAME");
        // mass_hint_cl95[i]->Draw("e3 same");
        // h_M_MC[i]->Draw("PESAME");
        // pdf_M[i]->Draw("LSAME");

        // mass_canvas_error[i] = new TCanvas(Form("mass_canvas_error_%s", name_DiMuSelection[i].Data()), Form("mass_canvas_error_%s", name_DiMuSelection[i].Data()), 1000, 1200);
        // mass_canvas_error[i]->cd();
        up_mass_hint_error_cl95[i] = (TH1D *)mass_hint_cl95[i]->Clone(Form("up_mass_hint_error_%s_cl95", name_DiMuSelection[i].Data()));
        up_mass_hint_error_cl95[i]->Reset();

        lo_mass_hint_error_cl95[i] = (TH1D *)mass_hint_cl95[i]->Clone(Form("lo_mass_hint_error_%s_cl95", name_DiMuSelection[i].Data()));
        lo_mass_hint_error_cl95[i]->Reset();

        up_mass_hint_error_cl99[i] = (TH1D *)mass_hint_cl99[i]->Clone(Form("up_mass_hint_error_%s_cl99", name_DiMuSelection[i].Data()));
        up_mass_hint_error_cl99[i]->Reset();

        lo_mass_hint_error_cl99[i] = (TH1D *)mass_hint_cl99[i]->Clone(Form("lo_mass_hint_error_%s_cl99", name_DiMuSelection[i].Data()));
        lo_mass_hint_error_cl99[i]->Reset();

        mass_canvas_error[i] = linear_error(i, up_mass_hint_error_cl95[i], lo_mass_hint_error_cl95[i], up_mass_hint_error_cl99[i], lo_mass_hint_error_cl99[i], mass_hint_cl95[i], mass_hint_cl99[i], color[i], 4, 30);
        mass_canvas_error[i]->SetName(Form("mass_canvas_error_%s", name_DiMuSelection[i].Data()));
        mass_canvas_error[i]->SetTitle(Form("mass_canvas_error_%s", name_DiMuSelection[i].Data()));

        mass_canvas_error[i]->SaveAs(Form("/home/michele_pennisi/cernbox/output_HF_dimuons/fit_data_output/plot/%s.png", mass_canvas_error[i]->GetName()));

        
        fOut->cd();
        printf("-------------%s-----------\n",name_DiMuSelection[i].Data());
        printf("name %s\n",up_pt_hint_error_cl95[i]->GetName());
        printf("name %s\n",up_mass_hint_error_cl95[i]->GetName());
        up_pt_hint_error_cl95[i]->Write(0, 2, 0);
        lo_pt_hint_error_cl95[i]->Write(0, 2, 0);

        up_mass_hint_error_cl95[i]->Write(0, 2, 0);
        lo_mass_hint_error_cl95[i]->Write(0, 2, 0);

        // up_mass_hint_error_cl95[i] = (TH1D *)h_M_MC[i]->Clone(Form("up_mass_hint_error_%s", name_DiMuSelection[i].Data()));
        // up_mass_hint_error_cl95[i]->Reset();
        // up_mass_hint_error_cl95[i]->GetYaxis()->SetRangeUser(-3.5, 3.5);

        // lo_mass_hint_error_cl95[i] = (TH1D *)h_M_MC[i]->Clone(Form("lo_mass_hint_error_%s", name_DiMuSelection[i].Data()));
        // lo_mass_hint_error_cl95[i]->Reset();

        // for (Int_t q = 0; q < up_mass_hint_error_cl95[i]->GetNbinsX(); q++)
        // {
        //     // cout << "hint[i]->GetBinError(q+1) " << mass_hint_cl95[i]->GetBinError(q + 1) << " hint[i]->GetBinContent(q+1)" << mass_hint_cl95[i]->GetBinContent(q + 1) << endl;
        //     if (mass_hint_cl95[i]->GetBinContent(q + 1) != 0)
        //     {

        //         up_mass_hint_error_cl95[i]->SetBinContent(q + 1, (Double_t)mass_hint_cl95[i]->GetBinError(q + 1) / mass_hint_cl95[i]->GetBinContent(q + 1));
        //         up_mass_hint_error_cl95[i]->SetBinError(q + 1, 0.1 * (Double_t)mass_hint_cl95[i]->GetBinError(q + 1) / mass_hint_cl95[i]->GetBinContent(q + 1));

        //         lo_mass_hint_error_cl95[i]->SetBinContent(q + 1, -(Double_t)mass_hint_cl95[i]->GetBinError(q + 1) / mass_hint_cl95[i]->GetBinContent(q + 1));
        //         lo_mass_hint_error_cl95[i]->SetBinError(q + 1, 0.1 * (Double_t)mass_hint_cl95[i]->GetBinError(q + 1) / mass_hint_cl95[i]->GetBinContent(q + 1));
        //     }
        //     else
        //     {

        //         up_mass_hint_error_cl95[i]->SetBinContent(q + 1, 0);
        //         up_mass_hint_error_cl95[i]->SetBinError(q + 1, 0);

        //         lo_mass_hint_error_cl95[i]->SetBinContent(q + 1, 0);
        //         lo_mass_hint_error_cl95[i]->SetBinError(q + 1, 0);
        //     }
        // }
        // cout << "hint[i]->GetBinError(q+1) " << up_mass_hint_error_cl95[i]->GetBinError(2) << " hint[i]->GetBinContent(q+1)" << up_mass_hint_error_cl95[i]->GetBinContent(2) << endl;
        // TLine *lo2up_line = new TLine(4, lo_mass_hint_error_cl95[i]->GetBinContent(1), 30, up_mass_hint_error_cl95[i]->GetBinContent(up_mass_hint_error_cl95[i]->GetNbinsX()));

        // TLine *up2lo_line = new TLine(4, up_mass_hint_error_cl95[i]->GetBinContent(1), 30, lo_mass_hint_error_cl95[i]->GetBinContent(up_mass_hint_error_cl95[i]->GetNbinsX()));

        // lo2up_line->SetLineColor(color[i]);
        // lo2up_line->SetLineWidth(3);
        // lo2up_line->SetLineStyle(2);
        // up2lo_line->SetLineColor(color[i]);
        // up2lo_line->SetLineWidth(3);
        // up2lo_line->SetLineStyle(2);

        // up_mass_hint_error_cl95[i]->SetMarkerStyle(20);
        // up_mass_hint_error_cl95[i]->SetMarkerSize(1.25);
        // up_mass_hint_error_cl95[i]->SetMarkerColor(color[i]);

        // lo_mass_hint_error_cl95[i]->SetMarkerStyle(20);
        // lo_mass_hint_error_cl95[i]->SetMarkerSize(1.25);
        // lo_mass_hint_error_cl95[i]->SetMarkerColor(color[i]);

        // up_mass_hint_error_cl95[i]->Draw("PE");
        // lo_mass_hint_error_cl95[i]->Draw("PESAME");

        // lo2up_line->Draw("same");
        // up2lo_line->Draw("same");
    }

    // break;
}

double FuncPtMass(double *x, double *par)
{
    return par[3] * (x[0] / TMath::Power(1 + TMath::Power(x[0] / par[0], par[1]), par[2]));
}

TCanvas *linear_error(Int_t i, TH1D *up_hint_error_cl95, TH1D *lo_hint_error_cl95, TH1D *up_hint_error_cl99, TH1D *lo_hint_error_cl99, TH1D *hint_cl95, TH1D *hint_cl99, Color_t color, Int_t minx = 0, Int_t max_x = 30)
{
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 1000);
    canvas->SetTicks();
    canvas->cd();

    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.11);

    for (Int_t q = 0; q < up_hint_error_cl99->GetNbinsX(); q++)
    {
        // cout << "hint_cl95->GetBinError(q+1) " << hint_cl95->GetBinError(q + 1) << " hint_cl95->GetBinContent(q+1)" << hint_cl95->GetBinContent(q + 1) << endl;
        if (hint_cl95->GetBinContent(q + 1) != 0)
        {

            up_hint_error_cl95->SetBinContent(q + 1, (Double_t)hint_cl95->GetBinError(q + 1) / hint_cl95->GetBinContent(q + 1));
            up_hint_error_cl95->SetBinError(q + 1, 0.1 * (Double_t)hint_cl95->GetBinError(q + 1) / hint_cl95->GetBinContent(q + 1));

            lo_hint_error_cl95->SetBinContent(q + 1, -(Double_t)hint_cl95->GetBinError(q + 1) / hint_cl95->GetBinContent(q + 1));
            lo_hint_error_cl95->SetBinError(q + 1, 0.1 * (Double_t)hint_cl95->GetBinError(q + 1) / hint_cl95->GetBinContent(q + 1));
        }
        else
        {

            up_hint_error_cl95->SetBinContent(q + 1, 0);
            up_hint_error_cl95->SetBinError(q + 1, 0);

            lo_hint_error_cl95->SetBinContent(q + 1, 0);
            lo_hint_error_cl95->SetBinError(q + 1, 0);
        }
    }

    for (Int_t q = 0; q < up_hint_error_cl99->GetNbinsX(); q++)
    {
        // cout << "hint_cl95->GetBinError(q+1) " << hint_cl95->GetBinError(q + 1) << " hint_cl95->GetBinContent(q+1)" << hint_cl95->GetBinContent(q + 1) << endl;
        if (hint_cl99->GetBinContent(q + 1) != 0)
        {

            up_hint_error_cl99->SetBinContent(q + 1, (Double_t)hint_cl99->GetBinError(q + 1) / hint_cl99->GetBinContent(q + 1));
            up_hint_error_cl99->SetBinError(q + 1, 0.1 * (Double_t)hint_cl99->GetBinError(q + 1) / hint_cl99->GetBinContent(q + 1));

            lo_hint_error_cl99->SetBinContent(q + 1, -(Double_t)hint_cl99->GetBinError(q + 1) / hint_cl99->GetBinContent(q + 1));
            lo_hint_error_cl99->SetBinError(q + 1, 0.1 * (Double_t)hint_cl99->GetBinError(q + 1) / hint_cl99->GetBinContent(q + 1));
        }
        else
        {

            up_hint_error_cl99->SetBinContent(q + 1, 0);
            up_hint_error_cl99->SetBinError(q + 1, 0);

            lo_hint_error_cl99->SetBinContent(q + 1, 0);
            lo_hint_error_cl99->SetBinError(q + 1, 0);
        }
    }

    up_hint_error_cl95->SetTitle(" ");
    up_hint_error_cl95->GetYaxis()->SetTitle("error/value pdf");
    up_hint_error_cl95->GetXaxis()->SetTitleOffset(1.2);
    up_hint_error_cl95->GetYaxis()->SetRangeUser(-up_hint_error_cl95->GetMaximum()*2, up_hint_error_cl95->GetMaximum()*2.);

    cout << "hint_cl95->GetBinError(q+1) " << up_hint_error_cl95->GetBinError(2) << " hint_cl95->GetBinContent(q+1)" << up_hint_error_cl95->GetBinContent(2) << endl;
    TLine *lo2up_line = new TLine(minx, lo_hint_error_cl95->GetBinContent(1), 30, up_hint_error_cl95->GetBinContent(up_hint_error_cl95->GetNbinsX()));

    TLine *up2lo_line = new TLine(minx, up_hint_error_cl95->GetBinContent(1), 30, lo_hint_error_cl95->GetBinContent(up_hint_error_cl95->GetNbinsX()));

    lo2up_line->SetLineColor(kRed);
    lo2up_line->SetLineWidth(3);
    lo2up_line->SetLineStyle(2);
    up2lo_line->SetLineColor(kRed);
    up2lo_line->SetLineWidth(3);
    up2lo_line->SetLineStyle(2);

    up_hint_error_cl95->SetMarkerStyle(20);
    up_hint_error_cl95->SetMarkerSize(1.25);
    up_hint_error_cl95->SetMarkerColor(color);

    lo_hint_error_cl95->SetMarkerStyle(20);
    lo_hint_error_cl95->SetMarkerSize(1.25);
    lo_hint_error_cl95->SetMarkerColor(color);

    up_hint_error_cl99->SetMarkerStyle(24);
    up_hint_error_cl99->SetMarkerSize(1.25);
    up_hint_error_cl99->SetMarkerColor(color);

    lo_hint_error_cl99->SetMarkerStyle(24);
    lo_hint_error_cl99->SetMarkerSize(1.25);
    lo_hint_error_cl99->SetMarkerColor(color);

    canvas->cd();

    up_hint_error_cl95->Draw("PE");
    lo_hint_error_cl95->Draw("PESAME");

    up_hint_error_cl99->Draw("PESAME");
    lo_hint_error_cl99->Draw("PESAME");

    lo2up_line->Draw("same");
    up2lo_line->Draw("same");

    TLegend *legend = new TLegend(0.215, 0.15, 0.45, 0.4);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0325);

    legend->SetTextAlign(12);
    legend->AddEntry(up_hint_error_cl95, "CL 95% ", "PE");
    legend->AddEntry(lo_hint_error_cl99, "CL 99% ", "PE");

    legend->Draw();

    return canvas;
}

TCanvas *printMC_ratio(TH1D *data, TF1 *pdf, TH1D *hint, Color_t color, Int_t minx = 0, Int_t max_x = 30)
{
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 900, 1000);
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
    TString str;
    str.Form("%s", data->GetTitle());

    //   if (str.Contains("h_M"))
    //   {
    //     data->GetXaxis()->SetTitle("#it{m}_{#mu^{#plus}#mu^{#plus}} (GeV/#it{c}^{2})");
    //     data->GetYaxis()->SetTitle("d#it{N}/d#it{m}_{#mu^{#plus}#mu^{#plus}} (GeV/#it{c}^{2})^{-1}");
    //   }
    //   else
    //   {
    //     data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    //     data->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}");
    //   }

    TH2D *h_grid = new TH2D("h_grid", "", 100, minx, max_x, 10000, 0.05, 120000);
    h_grid->SetTitle("");

    h_grid->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());
    h_grid->GetYaxis()->SetTitle(data->GetYaxis()->GetTitle());
    h_grid->GetXaxis()->SetTitleOffset(1.3);
    h_grid->GetXaxis()->SetTitleSize(0.0475);
    h_grid->GetXaxis()->SetLabelSize(0.045);

    h_grid->GetYaxis()->SetNdivisions(505);
    h_grid->GetYaxis()->SetTitleOffset(0.9);
    h_grid->GetYaxis()->SetTitleSize(0.065);
    h_grid->GetYaxis()->SetLabelSize(0.055);

    data->SetMarkerStyle(24);
    data->SetMarkerSize(1.25);
    data->SetMarkerColor(color);

    pdf->SetLineColor(color);
    pdf->SetLineWidth(2);
    pdf->SetLineStyle(9);
    pdf->SetNpx(10000);

    h_grid->Draw();
    data->Draw("PESAME");
    pdf->Draw("LSAME");
    hint->Draw("e3 same");
    data->Draw("PESAME");
    pdf->Draw("LSAME");

    // data->SetMarkerSize(1.75);
    // data->SetMarkerStyle(24);
    // data->SetMarkerColor(color);
    // data->SetLineColor(color);
    // data->SetLineWidth(2);
    // data->Scale(1. / data->Integral(), "width");
    // data->Draw("PESAME");

    // pdf->SetLineColor(color);
    // pdf->SetLineWidth(2);

    // pdf->Draw("LSAME");

    TLegend *legend = new TLegend(0.315, 0.15, 0.55, 0.4);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject*)0, "", "");
    legend->SetFillStyle(0);
    legend->SetLineColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0525);
    legend->SetHeader("MC");
    legend->SetTextAlign(11);
    legend->AddEntry(data, " ", "PE");

    TLegend *pdf_legend = new TLegend(0.375, 0.15, 0.635, 0.4);
    pdf_legend->SetTextAlign(11);
    pdf_legend->SetFillStyle(0);
    pdf_legend->SetBorderSize(0);
    pdf_legend->SetTextSize(0.0525);
    pdf_legend->SetHeader("PDF");
    pdf_legend->AddEntry(pdf, " ", "L");

    legend->Draw();
    pdf_legend->Draw();

    TLatex *letexTitle = new TLatex();
    letexTitle->SetNDC();
    letexTitle->SetTextFont(42);
    letexTitle->SetTextSize(0.05);
    if (str.Contains("Charm"))
        letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow c,c");
    if (str.Contains("Beauty"))
        letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow b,b");
    if (str.Contains("Mixed"))
        letexTitle->DrawLatex(0.16, 0.202, "#mu^{#plus}#mu^{#minus} #leftarrow b,c");
    // if (roohist_name.Contains("Charm")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow c,c");
    // if (roohist_name[option].Contains("Beauty")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,b");
    // if (roohist_name[option].Contains("Mixed")) letexTitle -> DrawLatex(0.20,0.202,"#mu^{#plus}#mu^{#minus} #leftarrow b,c");

    letexTitle->SetTextSize(0.06);
    // letexTitle -> DrawLatex(0.405,0.86,"ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    // printf("WOW %s\n",roohist_name[option].Data() );
    // letexTitle->SetTextSize(0.055);
    // letexTitle->DrawLatex(0.675, 0.825, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c} = %0.3f #pm %0.3f", fit_output[0]->getVal(), fit_output[0]->getError()));
    // letexTitle->DrawLatex(0.675, 0.725, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow b} = %0.3f #pm %0.3f", fit_output[1]->getVal(), fit_output[1]->getError()));
    // letexTitle->DrawLatex(0.675, 0.625, Form("#it{f}^{Output}_{#mu^{#plus}#mu^{#minus} #leftarrow c,b} = %0.3f", 1 - fit_output[0]->getVal() - fit_output[1]->getVal()));

    letexTitle->SetTextSize(0.055);
    letexTitle->DrawLatex(0.355, 0.88, "ALICE Simulation, pp #sqrt{#it{s}} = 13 TeV");
    letexTitle->DrawLatex(0.355, 0.81, "PYTHIA8 Monash Tune, N_{ev} = 2 #upoint 10^{8}");
    if (str.Contains("h_Pt"))
    {
        letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{#eta}_{#mu} < 4.0");
        letexTitle->DrawLatex(0.355, 0.66, Form("2.5 < #it{y}_{#mu^{#plus}#mu^{#minus}} < 4.0, 4 < #it{m}_{#mu^{#plus}#mu^{#minus}} < 30 GeV/#it{c}^{2}"));
    }
    else
    {
        letexTitle->DrawLatex(0.355, 0.74, "Reconstructed #mu^{#plus}#mu^{#minus}, 2.5 < #it{#eta}_{#mu} < 4.0");
        letexTitle->DrawLatex(0.355, 0.66, Form("2.5 < #it{y}_{#mu^{#plus}#mu^{#minus}} < 4.0"));
    }

    pad2->cd();
    pad2->SetTicks();
    TLine *l = new TLine(minx, 1.0, 30.0, 1.0);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->SetLineColor(kRed);

    TLine *l1 = new TLine(minx, 1.5, 30.0, 1.5);
    l1->SetLineWidth(2);
    l1->SetLineStyle(9);
    l1->SetLineColor(kGray + 2);

    TLine *l2 = new TLine(minx, 0.5, 30.0, 0.5);
    l2->SetLineWidth(2);
    l2->SetLineStyle(9);
    l2->SetLineColor(kGray + 2);

    TH2D *h_grid_ratio = new TH2D("h_grid_ratio", "", 100, minx, max_x, 100, 0.25, 2.55);
    h_grid_ratio->SetTitle("");

    TH1D *c_data = (TH1D *)data->Clone("c_data");
    TH1D *c_data_band = (TH1D *)data->Clone("c_data_band");

    h_grid_ratio->GetYaxis()->SetTitle(Form("#frac{MC}{PDF}"));
    h_grid_ratio->GetYaxis()->CenterTitle();
    h_grid_ratio->GetYaxis()->SetNdivisions(504);
    h_grid_ratio->GetYaxis()->SetTitleSize(0.08);
    // h_grid_ratio->GetYaxis()->SetTitleOffset(0.8);
    h_grid_ratio->GetYaxis()->SetLabelOffset(0.02);
    h_grid_ratio->GetYaxis()->SetLabelSize(0.1);

    h_grid_ratio->GetXaxis()->SetTitleSize(0.1);
    h_grid_ratio->GetXaxis()->SetTitleOffset(1.1);
    h_grid_ratio->GetXaxis()->SetLabelSize(0.1);
    h_grid_ratio->GetXaxis()->SetTitle(c_data->GetXaxis()->GetTitle());

    // c_data->SetLineColor(kBlack);
    // c_data->SetMarkerColor(kBlack);
    // c_data->SetMarkerStyle(20);

    // c_data->Rebin(10);
    // c_data->Scale(1. / c_data->Integral(), "width");
    c_data->Divide(pdf);
    c_data_band->Divide(hint);

    h_grid_ratio->Draw();
    c_data_band->SetFillColor(hint->GetFillColor());
    c_data_band->SetMinimum(0.2);
    c_data_band->Draw("e3gsame");
    l->Draw();
    l1->Draw();
    l2->Draw();
    c_data->Draw("PESAME");
    return canvas;
}
