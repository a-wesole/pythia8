// old version was called SignalSwapOnly.C
// updated to include all bins
#include <iostream>
#include "TRandom.h"
#include "TFile.h"
#include <iostream>
#include <random>
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include <TBranch.h>
#include <vector>
#include <TRandom3.h>
#include "TROOT.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TMath.h"
#include <TPDF.h>
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"

double histo_entries = 3.0e+08;
const double avg_mass = 1.865, std_mass = 1.902e-02, fg_lower = avg_mass - 3 * std_mass, fg_upper = avg_mass + 3 * std_mass; // July 30 updates
// const double fit_range_low = 1.82, fit_range_high = 1.9, D0_mass = 1.8648;
const double D0_mass = 1.8648;
TH2D *S1S2Template, *B1B2Template, *SW1SW2Template, *S1B2Template, *B1S2Template, *SW1B2Template, *B1SW2Template, *S1SW2Template, *SW1S2Template;
const int rebin_factor = 10;
const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / 5);
const float phi_array[5] = {in_val + sc_val, in_val + 2 * sc_val, in_val + 3 * sc_val, in_val + 4 * sc_val};
const bool check_templates = false;
// const float fit_range_low = 1.55,  fit_range_high = 2.2;
using namespace std;

Double_t MyCustomFunction(Double_t *x, Double_t *par)
{
    Double_t xx = x[0]; // x-coordinate
    Double_t yy = x[1]; // y-coordinate

    // Find bin indices once using any of the histograms, as they have the same binning
    Int_t binX = S1S2Template->GetXaxis()->FindBin(xx);
    Int_t binY = S1S2Template->GetYaxis()->FindBin(yy);

    // Retrieve bin contents using the same bin indices for all templates
    Double_t content_S1S2 = S1S2Template->GetBinContent(binX, binY);
    Double_t content_B1B2 = B1B2Template->GetBinContent(binX, binY);
    Double_t content_SW1SW2 = SW1SW2Template->GetBinContent(binX, binY);
    Double_t content_S1B2 = S1B2Template->GetBinContent(binX, binY);
    Double_t content_B1S2 = B1S2Template->GetBinContent(binX, binY);
    Double_t content_SW1B2 = SW1B2Template->GetBinContent(binX, binY);
    Double_t content_B1SW2 = B1SW2Template->GetBinContent(binX, binY);
    Double_t content_S1SW2 = S1SW2Template->GetBinContent(binX, binY);
    Double_t content_SW1S2 = SW1S2Template->GetBinContent(binX, binY);

    // Linear combination of templates with parameters
    return par[0] * content_S1S2 + par[1] * content_B1B2 + par[2] * content_SW1SW2 + par[3] * content_S1B2 + par[4] * content_B1S2 + par[5] * content_SW1S2 + par[6] * content_S1SW2 + par[7] * content_B1SW2 + par[8] * content_SW1B2;
}

TH2D *CreateTemplate(TF1 *fit1, TF1 *fit2, TH2F *reference_histo)
{
    // your function here
    TString newTitle = TString(reference_histo->GetTitle()) + "Template";

    TH2D *Template = new TH2D(newTitle, newTitle, 1000 / rebin_factor, 1.55, 2.2, 1000 / rebin_factor, 1.55, 2.2);

    for (int i = 0; i < histo_entries; i++)
    {
        float x = fit1->GetRandom();
        float y = fit2->GetRandom();
        Template->Fill(x, y);
    }
    // Template->Scale(reference_histo->GetEntries() / histo_entries);
    Template->Scale(1.0 / Template->GetEntries());
    return Template;
}

TH2D *CreateClone(TH2D *inputHisto, TH2D *templato)
{
    TString newTitle = TString(inputHisto->GetTitle()) + " - Template";
    TH2D *histoClone = new TH2D("histoClone", newTitle,
                                inputHisto->GetNbinsX(), inputHisto->GetXaxis()->GetXmin(), inputHisto->GetXaxis()->GetXmax(),
                                inputHisto->GetNbinsY(), inputHisto->GetYaxis()->GetXmin(), inputHisto->GetYaxis()->GetXmax());
    inputHisto->Copy(*histoClone);
    histoClone->Scale(1.0 / histoClone->GetEntries());
    histoClone->SetTitle(newTitle);
    histoClone->Add(templato, -1);
    histoClone->SetLineColor(1);
    histoClone->SetMaximum(templato->GetMaximum());
    return histoClone;
}

void TemplatesPhiCorr_differentialBinning(int bin)
{

   // TString outfile = TString(Form("templates_01_bin%d.root", bin));
    TFile *results = new TFile(outfile, "recreate");

    //~~First open the file nad retrieve all histograms & their x and y projections

    cout << "working to open file..." << endl;
    // TString inputfile1 = "TH2F_output_1k.root";
    // TString inputfile1 = "TH2F_08Nov_yesskip.root";
    // TString inputfile1 = "TH2F_12Nov_02.root";
    TString inputfile1 = "TH2F_08Nov_yesskip.root";
    TFile *inf1 = TFile::Open(inputfile1);
    if (!inf1)
        cout << "hey no inf1 !!!!!!!!!!!!!!!!" << endl;
    TString fits_file_path = "TF1_outputs_12Nov_01_allbins.root";
    // TString fits_file_path = "TF1_outputs_28Oct_swap_mean.root";
    TFile *inf2 = TFile::Open(fits_file_path);
    TString templates_path = Form("templates_01_bin%d.root", bin); // actual one we need
    TFile *inf3 = TFile::Open(templates_path);

        TH2F *M1M2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("M1M2Mass_%d", bin)));
        if (!M1M2Mass)
            cout << "hey no m1m2 mass .. !!!!!" << endl;
        M1M2Mass->RebinX(rebin_factor);
        M1M2Mass->RebinY(rebin_factor);
        TH1D *M1M2hx = M1M2Mass->ProjectionX("M1M2hx");
        TH1D *M1M2hy = M1M2Mass->ProjectionY("M1M2hy");
        // TH2F *M1M2Mass = new TH2D("M1M2Mass", "M1M2Mass", 100, 1.6, 2.1, 100, 1.6, 2.1); // Full Range of F1F2 Mass

        TH2F *S1S2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("S1S2Mass_%d", bin)));
        S1S2Mass->RebinX(rebin_factor);
        S1S2Mass->RebinY(rebin_factor);
        S1S2Mass->SetMinimum(0);
        TH1D *S1S2hx = S1S2Mass->ProjectionX("S1S2hx");
        TH1D *S1S2hy = S1S2Mass->ProjectionY("S1S2hy");

        TH2F *SignalSwap12Mass = dynamic_cast<TH2F *>(inf1->Get(Form("SignalSwapMass_%d", bin)));
        SignalSwap12Mass->SetMinimum(0);
        SignalSwap12Mass->RebinX(rebin_factor);
        SignalSwap12Mass->RebinY(rebin_factor);

        TH1D *SSWhx = SignalSwap12Mass->ProjectionX("SSWhx");
        TH1D *SSWhy = SignalSwap12Mass->ProjectionY("SSWhy");

        TH2F *SW1SW2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("SW1SW2Mass_%d", bin)));
        SW1SW2Mass->SetMinimum(0);
        SW1SW2Mass->RebinX(rebin_factor);
        SW1SW2Mass->RebinY(rebin_factor);
        TH1D *SWhx = SW1SW2Mass->ProjectionX("SWhx");
        TH1D *SWhy = SW1SW2Mass->ProjectionY("SWhy");

        TH2F *B1B2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("B1B2Mass_%d", bin)));
        B1B2Mass->SetMinimum(0);
        B1B2Mass->RebinX(rebin_factor);
        B1B2Mass->RebinY(rebin_factor);
        TH1D *Bkghx = B1B2Mass->ProjectionX("Bkghx");
        TH1D *Bkghy = B1B2Mass->ProjectionY("Bkghy");

        TH2F *S1SW2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("S1SW2Mass_%d", bin)));
        S1SW2Mass->SetMinimum(0);
        S1SW2Mass->RebinX(rebin_factor);
        S1SW2Mass->RebinY(rebin_factor);

        TH2F *SW1S2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("SW1S2Mass_%d", bin)));
        SW1S2Mass->SetMinimum(0);
        SW1S2Mass->RebinX(rebin_factor);
        SW1S2Mass->RebinY(rebin_factor);

        TH2F *S1B2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("S1B2Mass_%d", bin)));
        S1B2Mass->SetMinimum(0);
        S1B2Mass->RebinX(rebin_factor);
        S1B2Mass->RebinY(rebin_factor);

        TH2F *B1S2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("B1S2Mass_%d", bin)));
        B1S2Mass->SetMinimum(0);
        B1S2Mass->RebinX(rebin_factor);
        B1S2Mass->RebinY(rebin_factor);

        TH2F *SW1B2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("SW1B2Mass_%d", bin)));
        SW1B2Mass->SetMinimum(0);
        SW1B2Mass->RebinX(rebin_factor);
        SW1B2Mass->RebinY(rebin_factor);

        TH2F *B1SW2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("B1SW2Mass_%d", bin)));
        B1SW2Mass->SetMinimum(0);
        B1SW2Mass->RebinX(rebin_factor);
        B1SW2Mass->RebinY(rebin_factor);

        cout << "File Successfully Opened!" << endl;
        /*
        M1M2Mass->Add(S1S2Mass);
        //M1M2Mass->Add(S1B2Mass);
        //M1M2Mass->Add(B1S2Mass);
        M1M2Mass->Add(BkgOnlyMass);
        M1M2Mass->Add(SwapOnlyMass);
        */
        // histo_entries = M1M2Mass->GetEntries();

        //~~Next define a signal swap and bkg function for each d0 and dbar candidiate
        //~~Signal is double gaus, swap is double gaus (see notes) bkg is linear pol

        double fit_range_low = 1.55, fit_range_high = 2.2;

        
        S1S2Template = (TH2D *)inf3->Get("S1S2Mass_6Template");
        if (!S1S2Template)
            cout << "NO S1S2temp!!!!!!!!!!!!!!" << endl;
        S1SW2Template = (TH2D *)inf3->Get("S1SW2Mass_6Template");
        S1B2Template = (TH2D *)inf3->Get("S1B2Mass_6Template");
        B1S2Template = (TH2D *)inf3->Get("B1S2Mass_6Template");
        B1SW2Template = (TH2D *)inf3->Get("B1SW2Mass_6Template");
        B1B2Template = (TH2D *)inf3->Get("B1B2Mass_6Template");
        SW1S2Template = (TH2D *)inf3->Get("SW1S2Mass_6Template");
        SW1SW2Template = (TH2D *)inf3->Get("SW1SW2Mass_6Template");
        SW1B2Template = (TH2D *)inf3->Get("SW1B2Mass_6Template");

        // following is code to check each created template agaist data
        /*
        TH2D *cloneS1S2 = CreateClone(S1S2Mass, S1S2Template);
        TH2D *cloneS1SW2 = CreateClone(S1SW2Mass, S1SW2Template);
        TH2D *cloneS1B2 = CreateClone(S1B2Mass, S1B2Template);

        TH2D *cloneSW1S2 = CreateClone(SW1S2Mass, SW1S2Template);
        TH2D *cloneSW1SW2 = CreateClone(SW1SW2Mass, SW1SW2Template);
        TH2D *cloneSW1B2 = CreateClone(SW1B2Mass, SW1B2Template);

        TH2D *cloneB1S2 = CreateClone(B1S2Mass, B1S2Template);
        TH2D *cloneB1SW2 = CreateClone(B1SW2Mass, B1SW2Template);
        TH2D *cloneB1B2 = CreateClone(B1B2Mass, B1B2Template);

        TCanvas *temp_check = new TCanvas("temp_check", "temp_check", 1000, 1000);
        temp_check->Divide(3,3);
        temp_check->cd(1);
        cloneS1S2->Draw("surf1");
        temp_check->cd(2);
        cloneS1SW2->Draw("surf1");
        temp_check->cd(3);
        cloneS1B2->Draw("surf1");
        temp_check->cd(4);
        cloneSW1S2->Draw("surf1");
        temp_check->cd(5);
        cloneSW1SW2->Draw("surf1");
        temp_check->cd(6);
        cloneSW1B2->Draw("surf1");
        temp_check->cd(7);
        cloneB1S2->Draw("surf1");
        temp_check->cd(8);
        cloneB1SW2->Draw("surf1");
        temp_check->cd(9);
        cloneB1B2->Draw("surf1");
        //temp_check->SaveAs("pdfs/temp_check_08Nov.pdf");

        TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
        c1->cd();
        Bkghx->Scale(1.0 / Bkghx->GetEntries());
        Bkghx->Draw();
        B1B2Template->ProjectionX()->SetLineColor(kRed);
        B1B2Template->ProjectionX()->Draw("same");
        */

        // the following code only needs to be run if creating templates from scratch. otherwise read them in from the input file
        /*

        // TF1 *F1 = (TF1 *)inf2->Get("F1_6");
        cout << " <<<<<<<<<<<<" << endl;
        cout << " << BIN = " << j << " << " << endl;
        cout << " <<<<<<<<<<<<" << endl;
        TF1 *F1 = (TF1 *)inf2->Get(Form("F1_%d", j));

        TF1 *signal1 = new TF1("signal1", "([6]*[7]*([0]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [5]*((1-[0]))*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]) + (1-[5])* (1-[0])*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4])))", fit_range_low, fit_range_high);
        signal1->FixParameter(0, F1->GetParameter(4));  // ratio btw 1 and 2
        signal1->FixParameter(1, F1->GetParameter(1));  // mean
        signal1->FixParameter(2, F1->GetParameter(2));  // sigma 1
        signal1->FixParameter(3, F1->GetParameter(3));  // sihgma 2
        signal1->FixParameter(4, F1->GetParameter(14)); // sigma 3
        signal1->FixParameter(5, F1->GetParameter(12)); // ratio btw gaus 2 & 3
        signal1->FixParameter(6, F1->GetParameter(0));  // scaling
        signal1->FixParameter(7, F1->GetParameter(5));  // signal fraction
        // signal1->Draw();

        TF1 *swap1 = new TF1("swap1", "[0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
        swap1->FixParameter(0, F1->GetParameter(0));  // norm
        swap1->FixParameter(1, F1->GetParameter(16)); // mean mass
        swap1->FixParameter(2, F1->GetParameter(7));  // sigma1
        swap1->FixParameter(3, F1->GetParameter(15)); // sigma2
        swap1->FixParameter(4, F1->GetParameter(13)); // ratio
        swap1->FixParameter(5, F1->GetParameter(5));  // signal fraction
        swap1->FixParameter(6, F1->GetParameter(6));  // smearing

        TF1 *background1 = new TF1("background1", "[8] + [9]*x + [10]*x*x + [11]*x*x*x", fit_range_low, fit_range_high);
        background1->FixParameter(8, F1->GetParameter(8));
        background1->FixParameter(9, F1->GetParameter(9));
        background1->FixParameter(10, F1->GetParameter(10));
        background1->FixParameter(11, F1->GetParameter(11));

        TF1 *F2 = (TF1 *)inf2->Get(Form("F2_%d", j));

        TF1 *signal2 = new TF1("signal2", "([6]*[7]*([0]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [5]*((1-[0]))*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]) + (1-[5])* (1-[0])*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4])))", fit_range_low, fit_range_high);
        signal2->FixParameter(0, F2->GetParameter(4));  // ratio btw 1 and 2
        signal2->FixParameter(1, F2->GetParameter(1));  // mean
        signal2->FixParameter(2, F2->GetParameter(2));  // sigma 1
        signal2->FixParameter(3, F2->GetParameter(3));  // sihgma 2
        signal2->FixParameter(4, F2->GetParameter(14)); // sigma 3
        signal2->FixParameter(5, F2->GetParameter(12)); // ratio btw gaus 2 & 3
        signal2->FixParameter(6, F2->GetParameter(0));  // scaling
        signal2->FixParameter(7, F2->GetParameter(5));  // signal fraction

        TF1 *swap2 = new TF1("swap2", "[0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
        swap2->FixParameter(0, F2->GetParameter(0));  // norm
        swap2->FixParameter(1, F2->GetParameter(16)); // mean mass
        swap2->FixParameter(2, F2->GetParameter(7));  // sigma1
        swap2->FixParameter(3, F2->GetParameter(15)); // sigma2
        swap2->FixParameter(4, F2->GetParameter(13)); // ratio
        swap2->FixParameter(5, F2->GetParameter(5));  // signal fraction
        swap2->FixParameter(6, F2->GetParameter(6));  // smearing

        TF1 *background2 = new TF1("background2", "[8] + [9]*x + [10]*x*x + [11]*x*x*x", fit_range_low, fit_range_high);
        background2->FixParameter(8, F2->GetParameter(8));
        background2->FixParameter(9, F2->GetParameter(9));
        background2->FixParameter(10, F2->GetParameter(10));
        background2->FixParameter(11, F2->GetParameter(11));

        cout << "Generating template (1) - S1S2" << endl;
        TH2D *S1S2Template = CreateTemplate(signal1, signal2, S1S2Mass);
        cout << "Generating template (2) - S1SW2" << endl;
        TH2D *S1SW2Template = CreateTemplate(signal1, swap2, S1SW2Mass);
        cout << "Generating template (3) - S1B2" << endl;
        TH2D *S1B2Template = CreateTemplate(signal1, background2, S1B2Mass);

        cout << "Generating template (4) - SW1S2" << endl;
        TH2D *SW1S2Template = CreateTemplate(swap1, signal2, SW1S2Mass);
        cout << "Generating template (5) - SW1SW2" << endl;
        TH2D *SW1SW2Template = CreateTemplate(swap1, swap2, SW1SW2Mass);
        cout << "Generating template (6) - SW1B2" << endl;
        TH2D *SW1B2Template = CreateTemplate(swap1, background2, SW1B2Mass);

        cout << "Generating template (7) - B1S2" << endl;
        TH2D *B1S2Template = CreateTemplate(background1, signal2, B1S2Mass);
        cout << "Generating template (8) - B1SW2" << endl;
        TH2D *B1SW2Template = CreateTemplate(background1, swap2, B1SW2Mass);
        cout << "Generating template (9) - B1B2" << endl;
        TH2D *B1B2Template = CreateTemplate(background1, background2, B1B2Mass);

        TCanvas *temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1000, 1000);
        temp_canvas->Divide(3, 3);

        temp_canvas->cd(1);
        S1S2Template->Draw("surf1");
        temp_canvas->cd(2);
        S1SW2Template->Draw("surf1");
        temp_canvas->cd(3);
        S1B2Template->Draw("surf1");
        temp_canvas->cd(4);
        SW1S2Template->Draw("surf1");
        temp_canvas->cd(5);
        SW1SW2Template->Draw("surf1");
        temp_canvas->cd(6);
        SW1B2Template->Draw("surf1");
        temp_canvas->cd(7);
        B1S2Template->Draw("surf1");
        temp_canvas->cd(8);
        B1SW2Template->Draw("surf1");
        temp_canvas->cd(9);
        B1B2Template->Draw("surf1");

        temp_canvas->SaveAs(Form("pdfs/templates_21Nov_bin_%d.pdf", j));
        */

        TCanvas *cg = new TCanvas("cg", "cg", 800, 1200);
        cg->Divide(2, 3);
        cg->cd(1);
        TCanvas *can = new TCanvas("can", "", 1000, 1000);
        can->Divide(2, 2);
        can->cd(1);

        cout << "defining fit function" << endl;

        TF2 *fitF1F2 = new TF2("fitF1F2", MyCustomFunction, 1.55, 2.2, 1.55, 2.2, 9);
        for (int i = 0; i < 9; i++)
        {
            fitF1F2->SetParameter(i, 100.);
        }
        fitF1F2->SetParNames("s1s2", "b1b2", "sw1sw2", "s1b2", "b1s2", "sw1s2", "s1sw2", "b1sw2", "sw1b2");
        fitF1F2->SetNpx(100);
        fitF1F2->SetNpy(100);

        for (int i = 0; i < 9; i++)
        {
            fitF1F2->SetParLimits(i, 000.0, M1M2Mass->GetEntries());
        }
        TH2D *M1M2MassClone = (TH2D *)M1M2Mass->Clone("M1M2MassClone");

        for (int i = 0; i < 9; i++)
        {
            fitF1F2->SetParameter(i, 100.0);
        }

        for (int i = 0; i < 1; i++)
        {
            // fitF1F2->SetParameter(0, S1S2Mass->GetEntries());     // s1s2
            fitF1F2->SetParameter(0, 100.);                     // s1s2
            /*
            fitF1F2->FixParameter(1, B1B2Mass->GetEntries());   // b1b2
            fitF1F2->FixParameter(2, SW1SW2Mass->GetEntries()); // sw1sw2
            fitF1F2->FixParameter(3, S1B2Mass->GetEntries());   // s1b2
            fitF1F2->FixParameter(4, B1S2Mass->GetEntries());   // b1s2
            fitF1F2->FixParameter(5, SW1S2Mass->GetEntries());  // sw1s2
            fitF1F2->FixParameter(6, S1SW2Mass->GetEntries());  // s1sw2
            fitF1F2->FixParameter(7, B1SW2Mass->GetEntries());  // b1sw2
            fitF1F2->FixParameter(8, SW1B2Mass->GetEntries());  // sw1b2
            */
            can->cd(1);
            M1M2Mass->Draw("surf1");

            can->cd(2);
            cout << "---------" << endl;
            cout << "iteration = " << endl;
            for (int i = 0; i < 10; i++)
            {
                cout << i << endl;
                M1M2MassClone->Fit(fitF1F2, "LQ");
            }
             M1M2MassClone->Fit(fitF1F2, "LM");
            can->cd(2);
            M1M2MassClone->Draw("surf1");

            can->cd(3);
            TH2D *FitCheck = new TH2D("FitCheck", "Data-Fit", 1000 / rebin_factor, 1.55, 2.2, 1000 / rebin_factor, 1.55, 2.2);
            // FitCheck->Add(M1M2MassClone);
            FitCheck->Add(M1M2Mass);
            FitCheck->Add(fitF1F2, -1);
            FitCheck->RebinX(2);
            FitCheck->RebinY(2);
            FitCheck->SetMaximum(M1M2Mass->GetMaximum());
            FitCheck->Draw("surf1");

            can->cd(4);
            fitF1F2->Draw("surf1");

            can->Update();
            can->cd();

            double fract = fitF1F2->GetParameter(0) / S1S2Mass->GetEntries();
            cout << " <<<<<<<<<<<<" << endl;
            cout << " << BIN = " << bin << " << " << endl;
            cout << " <<<<<<<<<<<<" << endl;
            cout << "fit S1S2 yield = " << fitF1F2->GetParameter(0) << endl;
            cout << "fit yield / data yield for S1S2 = " << fitF1F2->GetParameter(0) << " / " << S1S2Mass->GetEntries() << endl;
            cout << "= " << fract << endl;
            cout << " <<<<<<<<<<<<" << endl;
        }

        // fit_range_low -= 0.02;
        // fit_range_high +=0.02;

        // getchar();

        /*
        results->cd();

        S1S2Template->Write();
        S1SW2Template->Write();
        S1B2Template->Write();
        SW1S2Template->Write();
        SW1SW2Template->Write();
        SW1B2Template->Write();
        B1S2Template->Write();
        B1SW2Template->Write();
        B1B2Template->Write();
        can->Write();
        */
    /*
    M1M2MassClone->Write();
    M1M2Mass->Write();
    fitF1F2->Write();
    */

    results->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) { // Ensure that exactly one argument is provided
        std::cerr << "Usage: " << argv[0] << " <bin>" << std::endl;
        return 1;
    }

    int bin = std::atoi(argv[1]); // Convert the command-line argument to an integer
    TemplatesPhiCorr_differentialBinning(bin); // Call the function with the provided bin number

    return 0;
}

