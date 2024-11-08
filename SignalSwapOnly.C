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

double histo_entries = 1.0e+08;
const double avg_mass = 1.865, std_mass = 1.902e-02, fg_lower = avg_mass - 3 * std_mass, fg_upper = avg_mass + 3 * std_mass; // July 30 updates
// const double fit_range_low = 1.82, fit_range_high = 1.9, D0_mass = 1.8648;
const double D0_mass = 1.8648;
TH2D *S1S2Template, *B1B2Template, *SW1SW2Template, *S1B2Template, *B1S2Template, *SW1B2Template, *B1SW2Template, *S1SW2Template, *SW1S2Template;
const int rebin_factor = 10;
const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / 5);
const float phi_array[5] = {in_val + sc_val, in_val + 2 * sc_val, in_val + 3 * sc_val, in_val + 4 * sc_val};
const bool check_templates = false;
const float fit_range_low = 1.55,  fit_range_high = 2.2;
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

TH2D *CreateTemplate(TF1 *fit1, TF1 *fit2, TH2D *reference_histo)
{
    // your function here
    TString newTitle = TString(reference_histo->GetTitle()) + "Template";

    TH2D *Template = new TH2D(newTitle, newTitle, 1000 / rebin_factor, 1.6, 2.1, 1000 / rebin_factor, 1.6, 2.1);

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
    histoClone->SetTitle(newTitle);
    histoClone->Add(templato, -1);
    histoClone->SetLineColor(1);
    histoClone->SetMaximum(inputHisto->GetMaximum());
    return histoClone;
}


void SignalSwapOnly()
{
    TFile *outfile = new TFile("templates_out.root", "recreate");

    //~~First open the file nad retrieve all histograms & their x and y projections

    cout << "working to open file..." << endl;
    //TString inputfile1 = "TH2F_output_1k.root";
    TString inputfile1 = "data_files/TH2F_06Nov.root";
    TFile *inf1 = TFile::Open(inputfile1);
    TString fits_file_path = "TF1_outputs_06Nov_swap_mean.root";
    //TString fits_file_path = "TF1_outputs_28Oct_swap_mean.root";
    TFile *inf2 = TFile::Open(fits_file_path);
    /*
    TString templates_path = "templates_out.root"; // actual one we need
    TFile *inf3 = TFile::Open(templates_path);
    */

    //TH2D *M1M2Mass = (TH2D *)inf1->Get("M1M2Mass_6"); // Full Range of F1F2 Mass
    TH2D *M1M2Mass = new TH2D("M1M2Mass", "M1M2Mass", 100, 1.6, 2.1, 100, 1.6, 2.1); // Full Range of F1F2 Mass
    //M1M2Mass->RebinX(rebin_factor);
    //M1M2Mass->RebinY(rebin_factor);
    TH1D *M1M2hx = M1M2Mass->ProjectionX("M1M2hx");
    TH1D *M1M2hy = M1M2Mass->ProjectionY("M1M2hy");

    TH2D *S1S2Mass = (TH2D *)inf1->Get("S1S2Mass_6"); // Full Range of F1F2 Mass
    S1S2Mass->RebinX(rebin_factor);
    S1S2Mass->RebinY(rebin_factor);
    S1S2Mass->SetMinimum(0);
    TH1D *S1S2hx = S1S2Mass->ProjectionX("S1S2hx");
    TH1D *S1S2hy = S1S2Mass->ProjectionY("S1S2hy");

    TH2D *SignalSwap12Mass = (TH2D *)inf1->Get("SignalSwapMass_6"); // Full Range of F1F2 Mass
    SignalSwap12Mass->SetMinimum(0);
    SignalSwap12Mass->RebinX(rebin_factor);
    SignalSwap12Mass->RebinY(rebin_factor);

    TH1D *SSWhx = SignalSwap12Mass->ProjectionX("SSWhx");
    TH1D *SSWhy = SignalSwap12Mass->ProjectionY("SSWhy");

    TH2D *SwapOnlyMass = (TH2D *)inf1->Get("SW1SW2Mass_6"); // Full Range of F1F2 Mass
    SwapOnlyMass->SetMinimum(0);
    SwapOnlyMass->RebinX(rebin_factor);
    SwapOnlyMass->RebinY(rebin_factor);
    TH1D *SWhx = SwapOnlyMass->ProjectionX("SWhx");
    TH1D *SWhy = SwapOnlyMass->ProjectionY("SWhy");

    TH2D *BkgOnlyMass = (TH2D *)inf1->Get("B1B2Mass_6"); // Full Range of F1F2 Mass
    BkgOnlyMass->SetMinimum(0);
    BkgOnlyMass->RebinX(rebin_factor);
    BkgOnlyMass->RebinY(rebin_factor);
    TH1D *Bkghx = BkgOnlyMass->ProjectionX("Bkghx");
    TH1D *Bkghy = BkgOnlyMass->ProjectionY("Bkghy");

    TH2D *S1SW2Mass = (TH2D *)inf1->Get("S1SW2Mass_6"); // Full Range of F1F2 Mass
    S1SW2Mass->SetMinimum(0);
    S1SW2Mass->RebinX(rebin_factor);
    S1SW2Mass->RebinY(rebin_factor);

    TH2D *SW1S2Mass = (TH2D *)inf1->Get("SW1S2Mass_6"); // Full Range of F1F2 Mass
    SW1S2Mass->SetMinimum(0);
    SW1S2Mass->RebinX(rebin_factor);
    SW1S2Mass->RebinY(rebin_factor);
    TH2D *S1B2Mass = (TH2D *)inf1->Get("S1B2Mass_6"); // Full Range of F1F2 Mass
    S1B2Mass->SetMinimum(0);
    S1B2Mass->RebinX(rebin_factor);
    S1B2Mass->RebinY(rebin_factor);
    TH2D *B1S2Mass = (TH2D *)inf1->Get("B1S2Mass_6"); // Full Range of F1F2 Mass
    B1S2Mass->SetMinimum(0);
    B1S2Mass->RebinX(rebin_factor);
    B1S2Mass->RebinY(rebin_factor);
    TH2D *SW1B2Mass = (TH2D *)inf1->Get("SW1B2Mass_6"); // Full Range of F1F2 Mass
    SW1B2Mass->SetMinimum(0);
    SW1B2Mass->RebinX(rebin_factor);
    SW1B2Mass->RebinY(rebin_factor);
    TH2D *B1SW2Mass = (TH2D *)inf1->Get("B1SW2Mass_6"); // Full Range of F1F2 Mass
    B1SW2Mass->SetMinimum(0);
    B1SW2Mass->RebinX(rebin_factor);
    B1SW2Mass->RebinY(rebin_factor);

    cout << "File Successfully Opened!" << endl;
    M1M2Mass->Add(S1S2Mass);
    //M1M2Mass->Add(S1B2Mass);
    //M1M2Mass->Add(B1S2Mass);
    M1M2Mass->Add(BkgOnlyMass);
    M1M2Mass->Add(SwapOnlyMass);
    // histo_entries = M1M2Mass->GetEntries();

    //~~Next define a signal swap and bkg function for each d0 and dbar candidiate
    //~~Signal is double gaus, swap is double gaus (see notes) bkg is linear pol

    TH2D *S1S2Template = new TH2D("S1S2Template", "S1S2Template", 100, fit_range_low, fit_range_high, 100, fit_range_low, fit_range_high);
    TH2D *S1SW2Template = new TH2D("S1SW2Template", "S1SW2Template", 100, fit_range_low, fit_range_high, 100, fit_range_low, fit_range_high);
    TH2D *S1B2Template = new TH2D("S1B2Template", "S1B2Template", 100, fit_range_low, fit_range_high, 100, fit_range_low, fit_range_high);
    TH2D *B1S2Template = new TH2D("B1S2Template", "B1S2Template", 100, fit_range_low, fit_range_high, 100, fit_range_low, fit_range_high);
    TH2D *B1SW2Template = new TH2D("B1SW2Template", "B1SW2Template", 100, fit_range_low, fit_range_high, 100, fit_range_low, fit_range_high);
    TH2D *B1B2Template = new TH2D("B1B2Template", "B1B2Template", 100, fit_range_low, fit_range_high, 100, fit_range_low, fit_range_high);
    TH2D *SW1S2Template = new TH2D("SW1S2Template", "SW1S2Template", 100, fit_range_low, fit_range_high, 100, fit_range_low, fit_range_high);
    TH2D *SW1SW2Template = new TH2D("SW1SW2Template", "SW1SW2Template", 100, fit_range_low, fit_range_high, 100, fit_range_low, fit_range_high);
    TH2D *SW1B2Template = new TH2D("SW1B2Template", "SW1B2Template", 100, fit_range_low, fit_range_high, 100, fit_range_low, fit_range_high);


    //double fit_range_low = 1.8, fit_range_high = 1.92;
    /*

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
    */


    TF1 *F1 = (TF1 *)inf2->Get("F1_6");

    TF1 *signal1 = new TF1("signal1", "([6]*[7]*([0]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [5]*((1-[0]))*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]) + (1-[5])* (1-[0])*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4])))", fit_range_low, fit_range_high);
    signal1->FixParameter(0, F1->GetParameter(4));  // ratio btw 1 and 2
    signal1->FixParameter(1, F1->GetParameter(1));  // mean
    signal1->FixParameter(2, F1->GetParameter(2));  // sigma 1
    signal1->FixParameter(3, F1->GetParameter(3));  // sihgma 2
    signal1->FixParameter(4, F1->GetParameter(14)); // sigma 3
    signal1->FixParameter(5, F1->GetParameter(12)); // ratio btw gaus 2 & 3
    signal1->FixParameter(6, F1->GetParameter(0));  // scaling
    signal1->FixParameter(7, F1->GetParameter(5));  // signal fraction

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

    TF1 *F2 = (TF1 *)inf2->Get("F2_6");

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

   *CreateTemplate(TF1 *fit1, TF1 *fit2, TH2D *reference_histo)




    TCanvas *cg = new TCanvas("cg", "cg", 800, 1200);
    cg->Divide(2, 3);
    cg->cd(1);
    TCanvas *can = new TCanvas("can", "", 1000, 1000);
    can->Divide(2, 2);
    can->cd(1);

    cout << "defining fit function" << endl;

    TF2 *fitF1F2 = new TF2("fitF1F2", MyCustomFunction, 1.6, 2.1, 1.6, 2.1, 9);
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
    for (int i = 3; i < 4; i++)
    {
        fitF1F2->FixParameter(0, S1S2Mass->GetEntries());     // s1s2
        /*
        fitF1F2->FixParameter(1, BkgOnlyMass->GetEntries());  // b1b2
        fitF1F2->FixParameter(2, SwapOnlyMass->GetEntries()); // sw1sw2
        fitF1F2->FixParameter(3, S1B2Mass->GetEntries());     // s1b2
        fitF1F2->FixParameter(4, B1S2Mass->GetEntries());     // b1s2
        fitF1F2->FixParameter(5, SW1S2Mass->GetEntries());    // sw1s2
        fitF1F2->FixParameter(6, S1SW2Mass->GetEntries());    // s1sw2
        fitF1F2->FixParameter(7, B1SW2Mass->GetEntries());    // b1sw2
        fitF1F2->FixParameter(8, SW1B2Mass->GetEntries());    // sw1b2
        */
        for (int i=1; i<9; i++){
            fitF1F2->FixParameter(i, 0.0);
        }

        can->cd(1);
        M1M2Mass->Draw("surf1");

        can->cd(2);
        cout << "---------" << endl;
        cout << "Released parameters are:" << endl;
        if (i >= 0)
        {
            fitF1F2->ReleaseParameter(0); // s1s2
            fitF1F2->SetParameter(0, S1S2Mass->GetEntries());
            fitF1F2->SetParLimits(0, 0.00, M1M2Mass->GetEntries());
            cout << "par(0) - S1S2" << endl;
        }
        if (i >= 6)
        {
            fitF1F2->ReleaseParameter(3); // s1b2
            fitF1F2->ReleaseParameter(4); // b1s2
            fitF1F2->SetParLimits(3, 0.00, M1M2Mass->GetEntries());
            fitF1F2->SetParLimits(4, 0.00, M1M2Mass->GetEntries());
            cout << "par(3) - S1B2" << endl;
            cout << "par(4) - B1S2" << endl;
        }
        if (i == 3)
        {
            fitF1F2->ReleaseParameter(1); // b1b2
            fitF1F2->SetParLimits(1, 0.00, M1M2Mass->GetEntries());
            cout << "par(1) - B1B2" << endl;
        }
        if (i == 3)
        {
            fitF1F2->ReleaseParameter(2); // sw1sw2
            fitF1F2->SetParLimits(2, 0.00, M1M2Mass->GetEntries());
            cout << "par(2) - SW1SW2" << endl;
        }
        if (i >= 4)
        {
            fitF1F2->ReleaseParameter(5); // sw1s2
            fitF1F2->ReleaseParameter(6); // s1sw2
            fitF1F2->SetParLimits(5, 0.00, M1M2Mass->GetEntries());
            fitF1F2->SetParLimits(6, 0.00, M1M2Mass->GetEntries());
            cout << "par(5) - SW1S2" << endl;
            cout << "par(6) - S1SW2" << endl;
        }
        if (i >= 5)
        {
            fitF1F2->ReleaseParameter(7); // b1sw2
            fitF1F2->ReleaseParameter(8); // sw1b2
            fitF1F2->SetParLimits(7, 0.00, M1M2Mass->GetEntries());
            fitF1F2->SetParLimits(8, 0.00, M1M2Mass->GetEntries());
            cout << "par(7) - B1SW2" << endl;
            cout << "par(8) - SW1B2" << endl;
        }
        cout << "iteration = " << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << i << endl;
            M1M2MassClone->Fit(fitF1F2, "Q");
        }
        M1M2MassClone->Fit(fitF1F2, "M");
        can->cd(2);
        M1M2MassClone->Draw("surf1");

        can->cd(3);
        TH2D *FitCheck = new TH2D("FitCheck", "Data-Fit", 1000 / rebin_factor, 1.6, 2.1, 1000 / rebin_factor, 1.6, 2.1);
        //FitCheck->Add(M1M2MassClone);
        FitCheck->Add(M1M2Mass);
        FitCheck->Add(fitF1F2, -1);
        FitCheck->RebinX(2);
        FitCheck->RebinY(2);
        FitCheck->Draw("surf1");

        can->cd(4);
        fitF1F2->Draw("surf1");

        can->Update();
        can->cd();

        double fract = fitF1F2->GetParameter(0) / S1S2Mass->GetEntries();
        cout << " <<<<<<<<<<<<" << endl;
        cout << "fit S1S2 yield = " << fitF1F2->GetParameter(0) << endl;
        cout << "fit yield / data yield for S1S2 = " << fitF1F2->GetParameter(0) << " / " << S1S2Mass->GetEntries() << endl;
        cout << "= " << fract << endl;
        cout << " <<<<<<<<<<<<" << endl;

        getchar();
    }


    outfile->cd();

    M1M2MassClone->Write();
    M1M2Mass->Write();
    fitF1F2->Write();

    outfile->Close();
}
