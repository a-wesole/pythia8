#include <iostream>
#include <random>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TNtuple.h"
#include <TTree.h>
#include "TTree.h"
#include <TBranch.h>
#include <vector>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom3.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include <TF1.h>

//const char *file_path = "/scratch/bell/awesole/mergedfile29Sept2024_cuts.root";
const char *file_path = "data_files/mergedfile.root";
//const char *limits_path = "sidebands_001.root";
//const char *limits_path = "test.root";
// const char *file_path = "/scratch/bell/awesole/test_57.root";
const float upper_limit = 3.14159, lower_limit = -3.14159;                                             // bounds to transition delta phi from (-2pi, 2pi) to (-0.5pi, 1.5*pi)
const float nbins_mass = 1.0e+04, xmin_mass = 1.55, xmax_mass = 2.2;                                    // for mass plots
const float nbins = 25.0, xmin = -0.5 * TMath::Pi(), xmax = 1.5 * TMath::Pi(), ymin = 0.0, ymax = 1.5; // for dphi plots
// const double avg_mass = 1.865, std_mass = 1.902e-02;                       // July 30 updates
// const double avg_mass = 1.86502, std_mass = 0.0147844;                       // July 30 updates
//  below are the limits for sidebands...
const float pT_minimum = 0.5; //   GeV/c
/*
const float avg_massD[6] = {1.806146, 1.804443, 1.805291, 1.802582, 1.805001, 1.805216};
const float std_massD[6] = {1.923592, 1.925303, 1.924445, 1.927122, 1.924667, 1.924508};
const float sb_l1D[6] = {1.746271, 1.742568, 1.744541, 1.738957, 1.743876, 1.744341}; // 99.7 closure   fraction=99.7 data style
const float sb_u2D[6] = {1.983467, 1.987178, 1.985195, 1.990747, 1.985792, 1.985383}; // 99.7% closure  fraction=99.7 data style

const float avg_massDbar[6] = {1.803993, 1.803237, 1.806124, 1.800378, 1.801654, 1.806551};
const float std_massDbar[6] = {1.925499, 1.926383, 1.923822, 1.929234, 1.927730, 1.923143};
const float sb_l1Dbar[6] = {1.741868, 1.740362, 1.746124, 1.734503, 1.737279, 1.747176}; // 99.7 closure   fraction=99.7 data style
const float sb_u2Dbar[6] = {1.987624, 1.989258, 1.983822, 1.995109, 1.992105, 1.982518}; // 99.7% closure  fraction=99.7 data style
*/

const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / 5);
const float phi_array[6] = {in_val, in_val + sc_val, in_val + 2 * sc_val, in_val + 3 * sc_val, in_val + 4 * sc_val, in_val + 5 * sc_val};

int bin = 0;
//int goodbin = 6;

double step_size = 0.00001;
bool fisB1, fisB2, fisS1, fisS2, fisF1, fisF2, fisSB1, fisSB2;
bool gisB1, gisB2, gisS1, gisS2, gisF1, gisF2, gisSB1, gisSB2;
bool fisD0candidate, fisDbarcandidate, gisD0candidate, gisDbarcandidate;
bool f_foreground, f_sidebands, g_foreground, g_sidebands;

// const float prob_cutoff = 0.1;
//const float prob_cutoff = 0.01;

using namespace std;

float transition_phi(float &D0del_phi)
{
    // this function reads in an angle d0del_phi and changes it so that the range is confined to (-pi/2 to 3pi/2)
    if (D0del_phi > upper_limit)
    {
        D0del_phi = D0del_phi - 2 * TMath::Pi();
    }
    if (D0del_phi < lower_limit)
    {
        D0del_phi = 2 * TMath::Pi() + D0del_phi;
    }
    if (D0del_phi < -0.5 * TMath::Pi() && D0del_phi > -1 * TMath::Pi())
    {
        D0del_phi = 2 * TMath::Pi() + D0del_phi;
    }
    return D0del_phi;
}

int temp_5bins(int goodbin, float prob_cutoff, TString outfile, TString limits_path)
{
    // TString outfile = TString(Form("output/bi6_noskip.root", goodbin));
    //TString outfile = TString("newCode_oldData/SB_001.root");
    // TString outfile = TString("output/albins_triplefit_09Oct.root");

    std::vector<float> *mass = nullptr;
    std::vector<float> *phi = nullptr;
    std::vector<float> *prob = nullptr;
    std::vector<float> *vertex_no = nullptr;

    std::vector<float> *iparticle_px = nullptr;
    std::vector<float> *iparticle_py = nullptr;
    std::vector<float> *iparticle_pz = nullptr;

    std::vector<float> *jparticle_px = nullptr;
    std::vector<float> *jparticle_py = nullptr;
    std::vector<float> *jparticle_pz = nullptr;
    std::vector<float> *phy_process = nullptr;
    std::vector<float> *TruePdg_i = nullptr;
    std::vector<float> *TruePdg_j = nullptr;
    std::vector<float> *AssignedPdg_i = nullptr;
    std::vector<float> *AssignedPdg_j = nullptr;
    std::vector<float> *fromD0 = nullptr;
    // std::vector<float> F1Mass, F2Mass, SB1Mass, SB2Mass, S1Mass, S2Mass, FullRangeF1, FullRangeF2, FullRangeF1Mass, FullRangeF2Mass;

    /*
    TH1F *inclusiveSignal = new TH1F("inclusiveSignal", "inclusive signal dphi", nbins, xmin, xmax);
    TH1F *correlated_dphi = new TH1F("correlated_dphi", "correlated dphi", nbins, xmin, xmax);
    TH1F *uncorrelated_dphi = new TH1F("uncorrelated_dphi", "uncorrelated dphi", nbins, xmin, xmax);
    TH1F *background_dphi = new TH1F("background_dphi", "background dphi", nbins, xmin, xmax);
    TH1F *Gluon_splitting = new TH1F("Gluon_splitting", "gluon splitting", nbins, xmin, xmax);
    TH1F *Gluon_fusion = new TH1F("Gluon_fusion", "gluon fusion", nbins, xmin, xmax);
    TH1F *flavor_excitation_quark = new TH1F("flavor_excitation_quark", "flavor_excitation_quark ", nbins, xmin, xmax);
    TH1F *flavor_excitation_gluon = new TH1F("flavor_excitation_gluon", "flavor_excitation_gluon ", nbins, xmin, xmax);
    TH1F *Quark_annihlation = new TH1F("Quark_annihlation", "Quark_annihlation ", nbins, xmin, xmax);
    TH1F *nonprompt_phi = new TH1F("nonprompt_dphi", "nonprompt dphi", nbins, xmin, xmax);

    */

    TH1F *F1Mass = new TH1F("F1", "F1", nbins_mass, xmin_mass, xmax_mass); // Foreground of D0 candidates
    TH1F *F2Mass = new TH1F("F2", "F2", nbins_mass, xmin_mass, xmax_mass); // foreground D0bar candidates
    TH1F *inclusiveSignal = new TH1F("InclusiveSignal", "", nbins, xmin, xmax);
    // TH1F *B1Mass = new TH1F("B1Mass", "B1 Mass", nbins_mass, xmin_mass, xmax_mass);    // Background of D0 candidates
    // TH1F *B2Mass = new TH1F("B2Mass", "B2 Mass", nbins_mass, xmin_mass, xmax_mass);    // Background D0bar candidates
    // TH1F *S1Mass = new TH1F("S1Mass", "S1 Mass", nbins_mass, xmin_mass, xmax_mass);    // D0 signal
    // TH1F *S2Mass = new TH1F("S2Mass", "S2 Mass", nbins_mass, xmin_mass, xmax_mass);    // D0bar signal
    // TH1F *SW1Mass = new TH1F("SW1Mass", "SW1 Mass", nbins_mass, xmin_mass, xmax_mass); // SW1 mass in 3 sigma
    // TH1F *SW2Mass = new TH1F("SW2Mass", "SW2 Mass", nbins_mass, xmin_mass, xmax_mass); // SW2 mass in 3 sigma
    // TH1F *SBSW1Mass = new TH1F("SBSW1Mass", "SBSW1Mass", nbins_mass, xmin_mass, xmax_mass);
    // TH1F *SBSW2Mass = new TH1F("SBSW2Mass", "SBSW2Mass", nbins_mass, xmin_mass, xmax_mass);
    // TH1F *SB1Mass = new TH1F("SB1Mass", "SB1Mass", nbins_mass, xmin_mass, xmax_mass);
    // TH1F *SB2Mass = new TH1F("SB2Mass", "SB2Mass", nbins_mass, xmin_mass, xmax_mass);

    TH1F *F1F2 = new TH1F("F1F2", "", nbins, xmin, xmax);
    // F1F2->Sumw2();
    TH1F *F1SB2 = new TH1F("F1SB2", "", nbins, xmin, xmax);
    // F1SB2->Sumw2();
    TH1F *SB1F2 = new TH1F("SB1F2", "", nbins, xmin, xmax);
    // SB1F2->Sumw2();
    TH1F *SB1SB2 = new TH1F("SB1SB2", "", nbins, xmin, xmax);
    // SB1SB2->Sumw2();
    // TH1F *S1S2 = new TH1F("S1S2", "S1S2", nbins, xmin, xmax);
    TH1F *B1B2 = new TH1F("B1B2", "", nbins, xmin, xmax);
    // B1B2->Sumw2();
    TH1F *F1B2 = new TH1F("F1B2", "", nbins, xmin, xmax);
    // F1B2->Sumw2();
    TH1F *B1F2 = new TH1F("B1F2", "", nbins, xmin, xmax);
    // B1F2->Sumw2();

    TH1F *S1S2 = new TH1F("S1S2", "", nbins, xmin, xmax);
    // S1S2->Sumw2();
    TH1F *S1B2 = new TH1F("S1B2", "", nbins, xmin, xmax);
    // S1B2->Sumw2();
    TH1F *B1S2 = new TH1F("B1S2", "", nbins, xmin, xmax);
    // B1S2->Sumw2();

    TH2F *FullRangeM1M2Mass = new TH2F("FullRangeM1M2", "", 100.0, 1.55, 2.2, 100.0, 1.55, 2.2);
    FullRangeM1M2Mass->SetXTitle("M1 Mass");
    FullRangeM1M2Mass->SetYTitle("M2 Mass");
    FullRangeM1M2Mass->SetOption("SURF1");

    TH2F *S1S2inclusiveSignalMass = new TH2F("InclusiveSIgnalMass", "", 16.0, 1.55, 2.2, 16.0, 1.55, 2.2);
    S1S2inclusiveSignalMass->SetXTitle("M1 Mass");
    S1S2inclusiveSignalMass->SetYTitle("M2 Mass");
    S1S2inclusiveSignalMass->SetOption("SURF1");

    float D_l1, D_l2, D_l3, D_l4;
    float Dbar_l1, Dbar_l2, Dbar_l3, Dbar_l4;

    // bool debug = false;
    float Ifile, event_no, error_files = 0.;
    float pT_i, pT_j;

    TFile *inf1 = TFile::Open(limits_path);
    TTree *D_limits = (TTree *)inf1->Get("D_mass_values");
    D_limits->SetBranchAddress("l1", &D_l1);
    D_limits->SetBranchAddress("l2", &D_l2);
    D_limits->SetBranchAddress("l3", &D_l3);
    D_limits->SetBranchAddress("l4", &D_l4);
    TTree *Dbar_limits = (TTree *)inf1->Get("Dbar_mass_values");
    Dbar_limits->SetBranchAddress("l1", &Dbar_l1);
    Dbar_limits->SetBranchAddress("l2", &Dbar_l2);
    Dbar_limits->SetBranchAddress("l3", &Dbar_l3);
    Dbar_limits->SetBranchAddress("l4", &Dbar_l4);

    TFile *infile = new TFile(file_path);
    TTree *t = (TTree *)infile->Get("event_tree");

    t->SetBranchAddress("mass", &mass);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("prob", &prob);
    t->SetBranchAddress("vertex_no", &vertex_no);
    t->SetBranchAddress("iparticle_px", &iparticle_px);
    t->SetBranchAddress("iparticle_py", &iparticle_py);
    t->SetBranchAddress("iparticle_pz", &iparticle_pz);
    t->SetBranchAddress("jparticle_px", &jparticle_px);
    t->SetBranchAddress("jparticle_py", &jparticle_py);
    t->SetBranchAddress("jparticle_pz", &jparticle_pz);
    t->SetBranchAddress("phy_process", &phy_process);
    t->SetBranchAddress("TruePdg_i", &TruePdg_i);
    t->SetBranchAddress("TruePdg_j", &TruePdg_j);
    t->SetBranchAddress("AssignedPdg_i", &AssignedPdg_i);
    t->SetBranchAddress("AssignedPdg_j", &AssignedPdg_j);
    t->SetBranchAddress("event_no", &event_no);
    t->SetBranchAddress("Ifile", &Ifile);
    t->SetBranchAddress("fromD0", &fromD0);

    for (int i = 0; i < t->GetEntries(); i++)
    // for (int i = 5000; i < 10000; i++)
    {
        t->GetEntry(i);
        /*
        if (Ifile != 2)
            continueVVVVVVVV;
        if (event_no < 13)
            continue;
        if (event_no > 14)
            break;
            */
        // cout << "Ifile =" << Ifile << " and event_no =" << event_no << endl;
        if (i % 6000 == 0)
            // if (i % 1 == 0)
            cout << i << " / " << t->GetEntries() << "  " << 100 * i / t->GetEntries() << "%" << endl;

        float N_SB1F2 = 0.0, N_B1F2 = 0.0, N_SB1SB2 = 0.0, N_B1B2 = 0.0, N_F1B2 = 0.0, N_F1SB2 = 0.0;

        for (int f = 0; f < phy_process->size() - 1; f++)
        { // for first kpi pair

            /*
            cout << "++++++++++++" << endl;
            cout << "f = " << f << endl;
            cout << "SB1F2->GetEntries()=" << SB1F2->GetEntries() << endl;
            */

            fisB1 = fisB2 = fisS1 = fisS2 = fisF1 = fisF2 = fisSB1 = fisSB2 = fisD0candidate = fisDbarcandidate = f_foreground = f_sidebands = false;

            pT_i = 0.0;
            pT_j = 0.0;
            if (prob->at(f) > prob_cutoff && fromD0->at(f) == 0)
                continue; // background cuts
            if (phy_process->at(f) > 5)
                continue; // skip all phy_process that are not 1-5

            pT_i = std::sqrt(iparticle_px->at(f) * iparticle_px->at(f) + iparticle_py->at(f) * iparticle_py->at(f));
            pT_j = std::sqrt(jparticle_px->at(f) * jparticle_px->at(f) + jparticle_py->at(f) * jparticle_py->at(f));
            if (pT_i < pT_minimum || pT_j < pT_minimum)
                continue; // apply pT cuts tp each daughter (iparticle and j particle)

            // if(fromD0->at(f)!=0 && AssignedPdg_i->at(f) == TruePdg_i->at(f) && AssignedPdg_j->at(f) == TruePdg_j->at(f)) continue;//skip all signal
            //if (AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) continue; // skip all swap
            // if(fromD0->at(f)!=0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) continue;//skip all swap
            // if(fromD0->at(f) ==0) continue;

            // determine D0 or Dbar candidate
            if ((AssignedPdg_i->at(f) == 211 && AssignedPdg_j->at(f) == -321) || (AssignedPdg_i->at(f) == -321 && AssignedPdg_j->at(f) == 211))
                fisD0candidate = true;
            if ((AssignedPdg_i->at(f) == -211 && AssignedPdg_j->at(f) == 321) || (AssignedPdg_i->at(f) == 321 && AssignedPdg_j->at(f) == -211))
                fisDbarcandidate = true;

            ///////////////////// proceed to second kpi pair //////////////////////////////////////////////////////////////////////////////

            // for (int g = f + 1; g < phy_process->size(); g++)
            for (int g = f + 1; g < phy_process->size(); g++)
            { // for second kpi pair

                bin = 0;

                f_foreground = f_sidebands = false;
                gisB1 = gisB2 = gisS1 = gisS2 = gisF1 = gisF2 = gisSB1 = gisSB2 = gisD0candidate = gisDbarcandidate = g_foreground = g_sidebands = false;

                pT_i = 0.0;
                pT_j = 0.0;
                if (prob->at(g) > prob_cutoff && fromD0->at(g) == 0)
                    continue; // background cuts
                if (phy_process->at(g) > 5)
                    continue; // skip all phy_process that are not 1-5

                // if(fromD0->at(g)!=0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_j->at(g) == TruePdg_j->at(g)) continue;//skip all signal
                //if (AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) continue;
                // if(fromD0->at(g)!=0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g))
                // if(fromD0->at(g) ==0) continue;

                pT_i = std::sqrt(iparticle_px->at(g) * iparticle_px->at(g) + iparticle_py->at(g) * iparticle_py->at(g));
                pT_j = std::sqrt(jparticle_px->at(g) * jparticle_px->at(g) + jparticle_py->at(g) * jparticle_py->at(g));
                if (pT_i < pT_minimum || pT_j < pT_minimum)
                    continue; // apply pT cuts to each daughter (iparticle & jparticle)
                if (iparticle_px->at(f) == iparticle_px->at(g) && iparticle_py->at(f) == iparticle_py->at(g) && iparticle_pz->at(f) == iparticle_pz->at(g) &&
                    jparticle_px->at(f) == jparticle_px->at(g) && jparticle_py->at(f) == jparticle_py->at(g) && jparticle_pz->at(f) == jparticle_pz->at(g)){
                    continue;
                }

                // determine D0 and Dbar candidate:
                if ((AssignedPdg_i->at(g) == 211 && AssignedPdg_j->at(g) == -321) || (AssignedPdg_i->at(g) == -321 && AssignedPdg_j->at(g) == 211))
                    gisD0candidate = true;
                if ((AssignedPdg_i->at(g) == -211 && AssignedPdg_j->at(g) == 321) || (AssignedPdg_i->at(g) == 321 && AssignedPdg_j->at(g) == -211))
                    gisDbarcandidate = true;

                /////////////////////////////////////////////////////////////////////////////
                // begin to classify the DDbar phi
                //  fromD0==0 : kpi does not come from D0
                //  fromD0 !=0 : kpi comes from D0
                // vertex_no ==0 : not signal
                // vertex_no !=0 : true signal
                /////////////////////////////////////////////////////////////////////////////

                if (fisD0candidate && gisDbarcandidate)
                { // d0 then d0bar
                    // D0candidate then D0barcandidate
                    float delta_phi = phi->at(f) - phi->at(g);
                    transition_phi(delta_phi);
                    // if (abs(delta_phi)<1.0e-6) continue;


                    if (delta_phi >= phi_array[0] && delta_phi < phi_array[1])
                        bin = 1;
                    if (delta_phi >= phi_array[1] && delta_phi < phi_array[2])
                        bin = 2;
                    if (delta_phi >= phi_array[2] && delta_phi < phi_array[3])
                        bin = 3;
                    if (delta_phi >= phi_array[3] && delta_phi < phi_array[4])
                        bin = 4;
                    if (delta_phi >= phi_array[4] && delta_phi <= phi_array[5])
                        bin = 5;

                    if (bin != goodbin){
                        //cout << "Error BREAK  fD0 && bin = " << bin << endl;
                        continue;
                    }
                    // cout << "f d  bin=" << bin << endl;
                    /*
                    if (delta_phi < phi_array[3] && bin !=0) {
                          cout << "uh oh!  bin should=0 but instead bin = " << bin << endl;
                    }
                    if(bin ==0) continue;
                    */

                     D_limits->GetEntry(bin-1);
                     Dbar_limits->GetEntry(bin-1);
                    // determine if f is in foreground or sidebands ( mass regions)
                    // cout << "look here  avg_massD[bin-1] = " << avg_massD[bin-1] << endl;
                    f_foreground = f_sidebands = false;
                    if (f_foreground && f_sidebands)
                        cout << "BAD BAD BAD BAD BAD fFG 1 & SB" << endl;
                    if (mass->at(f) > D_l2 && mass->at(f) < D_l3)
                        f_foreground = true;

                    else if ((mass->at(f) <= D_l2 && mass->at(f) >= D_l1) || (mass->at(f) >= D_l3 && mass->at(f) <= D_l4))
                        f_sidebands = true;

                    /////////////////////determine background, swap, signal and sidebands components for FIRST kpi pair////////////////////////////

                    if (f_foreground)
                    { // mass is within foreground region (3 sigma)
                        fisF1 = fisF2 = fisB1 = fisB2 = fisS1 = fisS2 = false;
                        if (fisD0candidate)
                            fisF1 = true;
                        if (fisDbarcandidate)
                            fisF2 = true;
                        if (fromD0->at(f) == 0) // fromD0==0 means kpi do not go to D0
                        {
                            if (fisD0candidate)
                                fisB1 = true;
                            if (fisDbarcandidate)
                                fisB2 = true;
                        }
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == TruePdg_i->at(f) && AssignedPdg_j->at(f) == TruePdg_j->at(f)) // True Signal
                        {
                            if (fisD0candidate)
                                fisS1 = true;
                            if (fisDbarcandidate)
                                fisS2 = true;
                        }
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) // Swap
                        {
                            if (fisD0candidate)
                                fisB1 = true;
                            if (fisDbarcandidate)
                                fisB2 = true;
                        }
                    }

                    else if (f_sidebands)
                    { // mass in within sidebands region
                        fisSB1 = fisSB2 = false;
                        if (fisD0candidate)
                            fisSB1 = true;
                        if (fisDbarcandidate)
                            fisSB2 = true;
                        /*
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) // Sidebands - swap
                        {
                            if (fisD0candidate)
                                fisSB1 = true;
                            if (fisDbarcandidate)
                                fisSB2 = true;
                        }
                        if (fromD0->at(f) == 0) // sidebands - background
                        {
                            if (fisD0candidate)
                                fisSB1 = true;
                            if (fisDbarcandidate)
                                fisSB2 = true;
                        }
                        */
                    }
                    else
                        continue;
                    // determine if mass is foreground or sidebands region, else continue;
                    g_foreground = g_sidebands = false;

                    if (mass->at(g) > Dbar_l2 && mass->at(g) < Dbar_l3)
                        g_foreground = true;
                    else if ((mass->at(g) <= Dbar_l2 && mass->at(g) >= Dbar_l1) || (mass->at(g) >= Dbar_l3 && mass->at(g) <= Dbar_l4))
                        g_sidebands = true;

                    /////////////////////determine background, swap, signal and sidebands components for SECOND kpi pair////////////////////////////

                    if (g_foreground && g_sidebands)
                        cout << "BAD BAD BAD BAD BAD GFG & SB" << endl;
                    if (g_foreground)
                    { // mass is within foreground region (3 sigma)
                        gisF1 = gisF2 = gisB1 = gisB2 = gisS1 = gisS2 = false;
                        if (gisD0candidate)
                            gisF1 = true;
                        if (gisDbarcandidate)
                            gisF2 = true;
                        if (fromD0->at(g) == 0) // from D0 ==0 means background
                        {
                            if (gisD0candidate)
                                gisB1 = true;
                            if (gisDbarcandidate)
                                gisB2 = true;
                        }
                        if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_j->at(g) == TruePdg_j->at(g)) // Signal
                        {
                            if (gisD0candidate)
                                gisS1 = true;
                            if (gisDbarcandidate)
                                gisS2 = true;
                        }
                        if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) // Swap
                        {
                            if (gisD0candidate)
                                gisB1 = true;
                            if (gisDbarcandidate)
                                gisB2 = true;
                        }
                    }

                    else if (g_sidebands)
                    { // mass in within sidebands region
                        gisSB1 = gisSB2 = false;
                        if (gisD0candidate)
                            gisSB1 = true;
                        if (gisDbarcandidate)
                            gisSB2 = true;
                        /*
                    if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) // SB - swap
                    {
                        if (gisD0candidate)
                            gisSB1 = true;
                        if (gisDbarcandidate)
                            gisSB2 = true;
                    }
                    if (fromD0->at(g) == 0) // SB - background
                    {
                        if (gisD0candidate)
                            gisSB1 = true;
                        if (gisDbarcandidate)
                            gisSB2 = true;
                    }
                    */
                    }
                    else
                        continue;

                    if (fromD0->at(f) != 0 && fromD0->at(g) != 0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_i->at(f) == TruePdg_i->at(f) &&
                        AssignedPdg_j->at(g) == TruePdg_j->at(g) && AssignedPdg_j->at(f) == TruePdg_j->at(f) && f_foreground && g_foreground)
                    { // any kpi pair from DDbar pair

                        // inclusiveSignal.at(0)->Fill(delta_phi); // any ddbar pair from D0s -- inclusive signal
                        // S1S2inclusiveSignalMass.at(0)->Fill(mass->at(f), mass->at(g));
                        inclusiveSignal->Fill(delta_phi); // any ddbar pair from D0s -- inclusive signal
                        S1S2inclusiveSignalMass->Fill(mass->at(f), mass->at(g));

                        /*
                        if (fromD0->at(f) == -101 && fromD0->at(g) == -101)
                            nonprompt_phi->Fill(delta_phi);
                        if (((vertex_no->at(f) != vertex_no->at(g)) || (vertex_no->at(f) == 0 && vertex_no->at(g) == 0)) && (fromD0->at(f) != -101 || fromD0->at(g) != -101))
                        {
                            uncorrelated_dphi->Fill(delta_phi); // dphi from different origins, unocrrelated
                        }
                        if (vertex_no->at(f) == vertex_no->at(g) && vertex_no->at(f) < 0 && vertex_no->at(g) < 0)
                        {
                            correlated_dphi->Fill(delta_phi); // correlated dphi pair, true signal
                            if (phy_process->at(f) == 1)
                                Gluon_splitting->Fill(delta_phi);
                            if (phy_process->at(f) == 2)
                                Gluon_fusion->Fill(delta_phi);
                            if (phy_process->at(f) == 3)
                                flavor_excitation_gluon->Fill(delta_phi);
                            if (phy_process->at(f) == 4)
                                flavor_excitation_quark->Fill(delta_phi);
                            if (phy_process->at(f) == 5)
                                Quark_annihlation->Fill(delta_phi);
                        }*/
                    }
                    if (mass->at(f) > 1.55 && mass->at(f) < 2.2 && mass->at(g) > 1.55 && mass->at(g) < 2.2)
                    {
                        F1Mass->Fill(mass->at(f));
                        F2Mass->Fill(mass->at(g));
                        FullRangeM1M2Mass->Fill(mass->at(f), mass->at(g));
                    }

                    /*
                    if ((fromD0->at(f) == 0 || fromD0->at(g) == 0))
                        background_dphi.at(bin-1)->Fill(delta_phi); // background dphi
                        */

                    // for dphi correlation below!!!!!!!!!!!!!!!!!!!!//

                    /*
                    cout << "" << endl;
                    cout << "-------------- new particle pair--------------" << endl;
                    cout << "-----------------f is D0 cand-----------------" << endl;
                    cout << "particle f mass = " << mass->at(f) << " & g mass = " << mass->at(g) << endl;
                    cout << "delta phi = " << delta_phi << endl;
                    cout << "code determined the dphi correlation component to be " << endl;
                    */
                    if (fisF1 && gisF2)
                    {
                        F1F2->Fill(delta_phi);
                        // cout << "mass f1 = " << mass->at(f) << " and f2 =" << mass->at(g) << endl;
                        //  cout << "~F1F2" << endl;
                    }
                    if (fisF1 && gisSB2)
                    {
                        F1SB2->Fill(delta_phi);
                        N_F1SB2 += 1;
                        // cout << "mass f1 = " << mass->at(f) << " and sb2 =" << mass->at(g) << endl;
                        //  cout << "~F1SB2" << endl;
                    }
                    if (fisSB1 && gisF2)
                    {
                        SB1F2->Fill(delta_phi);
                        N_SB1F2 += 1;
                        // cout << "mass sb1 = " << mass->at(f) << " and f2 =" << mass->at(g) << endl;
                        //  cout << "~SB1F2" << endl;
                    }
                    if (fisSB1 && gisSB2)
                    {
                        // cout << "~SB1SB2" << endl;
                        SB1SB2->Fill(delta_phi);
                        N_SB1SB2 += 1;
                        // cout << "mass sb1 = " << mass->at(f) << " and sb2 =" << mass->at(g) << endl;
                        ////SB1Mass->Fill(mass->at(f));
                        // SB2Mass->Fill(mass->at(g));
                    }
                    if (fisB1 && gisB2)
                    {
                        ////cout << "~B1B2" << endl;
                        B1B2->Fill(delta_phi);
                        N_B1B2 += 1;
                        // cout << "mass b1 = " << mass->at(f) << " and b2 =" << mass->at(g) << endl;
                        //  B1Mass->Fill(mass->at(f));
                        //  B2Mass->Fill(mass->at(g));
                    }
                    if (fisB1 && gisF2)
                    {
                        // cout << "~B1F2" << endl;
                        N_B1F2 += 1;
                        B1F2->Fill(delta_phi);
                        // cout << "mass b1 = " << mass->at(f) << " and f2 =" << mass->at(g) << endl;
                    }
                    if (fisF1 && gisB2)
                    {
                        // cout << "~F1B2" << endl;
                        F1B2->Fill(delta_phi);
                        N_F1B2 += 1;
                        // cout << "mass f1 = " << mass->at(f) << " and b2 =" << mass->at(g) << endl;
                    }

                    // if(fisS1 && gisS2) S1S2->Fill(delta_phi);
                    if (fisS1 && gisB2)
                    {
                        ////cout << "~S1B2" << endl;
                        S1B2->Fill(delta_phi);
                    }
                    if (fisB1 && gisS2)
                    {
                        ////cout << "~B1S2" << endl;
                        B1S2->Fill(delta_phi);
                    }

                } // D0can then D0barcan

                else if (fisDbarcandidate && gisD0candidate)
                { // d0bar then d0
                  // if D0barcand then D0cand

                    float delta_phi = phi->at(g) - phi->at(f);
                    transition_phi(delta_phi);
                    // if(abs(delta_phi)<1.0e-6) continue;

                    if (delta_phi >= phi_array[0] && delta_phi < phi_array[1])
                        bin = 1;
                    if (delta_phi >= phi_array[1] && delta_phi < phi_array[2])
                        bin = 2;
                    if (delta_phi >= phi_array[2] && delta_phi < phi_array[3])
                        bin = 3;
                    if (delta_phi >= phi_array[3] && delta_phi < phi_array[4])
                        bin = 4;
                    if (delta_phi >= phi_array[4] && delta_phi <= phi_array[5])
                        bin = 5;

                    if (bin != goodbin) {
                        //cout << "Error BREAK  fDbar && bin = " << bin << endl;
                        continue;
                    }
                    D_limits->GetEntry(bin-1);
                    Dbar_limits->GetEntry(bin-1);
                    //D_limits->GetEntry(0);
                    //Dbar_limits->GetEntry(0);
                    // if (abs(delta_phi) < 1.0e-6) continue;
                    // cout << "f dbar  bin=" << bin << endl;
                    /*
                       if (delta_phi < phi_array[3] && bin !=0) {
                       cout << "uh oh!  bin should=0 but instead bin = " << bin << endl;
                       }
                       if(bin ==0) continue;
                       */

                    // determine if f is in foreground or sidebands ( mass regions)
                    f_foreground = f_sidebands = false;
                    if (mass->at(f) > Dbar_l2 && mass->at(f) < Dbar_l3)
                        f_foreground = true;
                    else if ((mass->at(f) <= Dbar_l2 && mass->at(f) >= Dbar_l1) || (mass->at(f) >= Dbar_l3 && mass->at(f) <= Dbar_l4))
                        f_sidebands = true;

                    if (f_foreground && f_sidebands)
                        cout << "BAD BAD BAD BAD BAD fFG 2 & SB" << endl;
                    if (f_foreground)
                    { // mass is within foreground region (3 sigma)
                        fisF1 = fisF2 = fisB1 = fisB2 = fisS1 = fisS2 = false;
                        if (fisD0candidate)
                            fisF1 = true;
                        if (fisDbarcandidate)
                            fisF2 = true;
                        if (fromD0->at(f) == 0) // fromD0==0 means kpi do not go to D0
                        {
                            if (fisD0candidate)
                                fisB1 = true;
                            if (fisDbarcandidate)
                                fisB2 = true;
                        }
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == TruePdg_i->at(f) && AssignedPdg_j->at(f) == TruePdg_j->at(f)) // True Signal
                        {
                            if (fisD0candidate)
                                fisS1 = true;
                            if (fisDbarcandidate)
                                fisS2 = true;
                        }
                        if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) // Swap
                        {
                            if (fisD0candidate)
                                fisB1 = true;
                            if (fisDbarcandidate)
                                fisB2 = true;
                        }
                    }

                    else if (f_sidebands)
                    { // mass in within sidebands region
                        fisSB1 = fisSB2 = false;
                        if (fisD0candidate)
                            fisSB1 = true;
                        if (fisDbarcandidate)
                            fisSB2 = true;

                        /*
                           if (fromD0->at(f) != 0 && AssignedPdg_i->at(f) == -TruePdg_j->at(f) && AssignedPdg_j->at(f) == -TruePdg_i->at(f)) // Sidebands - swap
                           {
                           if (fisD0candidate)
                           fisSB1 = true;
                           if (fisDbarcandidate)
                           fisSB2 = true;
                           }
                           if (fromD0->at(f) == 0) // sidebands - background
                           {
                           if (fisD0candidate)
                           fisSB1 = true;
                           if (fisDbarcandidate)
                           fisSB2 = true;
                           }
                           */
                    }
                    else
                        continue;
                    // determine if mass is foreground or sidebands region, else continue;
                    g_foreground = g_sidebands = false;
                    if (mass->at(g) > D_l2 && mass->at(g) < D_l3)
                        g_foreground = true;
                    else if ((mass->at(g) <= D_l2 && mass->at(g) >= D_l1) || (mass->at(g) >= D_l3 && mass->at(g) <= D_l4))
                        g_sidebands = true;

                    /////////////////////determine background, swap, signal and sidebands components for SECOND kpi pair////////////////////////////

                    if (g_foreground && g_sidebands)
                        cout << "BAD BAD BAD BAD BAD GFG & SB" << endl;
                    if (g_foreground)
                    { // mass is within foreground region (3 sigma)
                        gisF1 = gisF2 = gisB1 = gisB2 = gisS1 = gisS2 = false;
                        if (gisD0candidate)
                            gisF1 = true;
                        if (gisDbarcandidate)
                            gisF2 = true;
                        if (fromD0->at(g) == 0) // from D0 ==0 means background
                        {
                            if (gisD0candidate)
                                gisB1 = true;
                            if (gisDbarcandidate)
                                gisB2 = true;
                        }
                        if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_j->at(g) == TruePdg_j->at(g)) // Signal
                        {
                            if (gisD0candidate)
                                gisS1 = true;
                            if (gisDbarcandidate)
                                gisS2 = true;
                        }
                        if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) // Swap
                        {
                            if (gisD0candidate)
                                gisB1 = true;
                            if (gisDbarcandidate)
                                gisB2 = true;
                        }
                    }

                    else if (g_sidebands)
                    { // mass in within sidebands region
                        gisSB1 = gisSB2 = false;
                        if (gisD0candidate)
                            gisSB1 = true;
                        if (gisDbarcandidate)
                            gisSB2 = true;
                        /*
                           if (fromD0->at(g) != 0 && AssignedPdg_i->at(g) == -TruePdg_j->at(g) && AssignedPdg_j->at(g) == -TruePdg_i->at(g)) // SB - swap
                           {
                           if (gisD0candidate)
                           gisSB1 = true;
                           if (gisDbarcandidate)
                           gisSB2 = true;
                           }
                           if (fromD0->at(g) == 0) // SB - background
                           {
                           if (gisD0candidate)
                           gisSB1 = true;
                           if (gisDbarcandidate)
                           gisSB2 = true;
                           }
                           */
                    }
                    else
                        continue;

                    if (fromD0->at(f) != 0 && fromD0->at(g) != 0 && AssignedPdg_i->at(g) == TruePdg_i->at(g) && AssignedPdg_i->at(f) == TruePdg_i->at(f) &&
                        AssignedPdg_j->at(g) == TruePdg_j->at(g) && AssignedPdg_j->at(f) == TruePdg_j->at(f) && f_foreground && g_foreground)
                    {
                        inclusiveSignal->Fill(delta_phi);
                        S1S2inclusiveSignalMass->Fill(mass->at(g), mass->at(f));

                        /*
                           if (fromD0->at(f) == -101 && fromD0->at(g) == -101)
                           nonprompt_phi->Fill(delta_phi);
                           if (((vertex_no->at(f) != vertex_no->at(g)) || (vertex_no->at(f) == 0 && vertex_no->at(g) == 0)) && (fromD0->at(f) != -101 || fromD0->at(g) != -101))
                           uncorrelated_dphi->Fill(delta_phi);
                           if (vertex_no->at(f) == vertex_no->at(g) && vertex_no->at(f) < 0 && vertex_no->at(g) < 0)
                           {
                           correlated_dphi->Fill(delta_phi);
                           if (phy_process->at(f) == 1)
                           Gluon_splitting->Fill(delta_phi);
                           if (phy_process->at(f) == 2)
                           Gluon_fusion->Fill(delta_phi);
                           if (phy_process->at(f) == 3)
                           flavor_excitation_gluon->Fill(delta_phi);
                           if (phy_process->at(f) == 4)
                           flavor_excitation_quark->Fill(delta_phi);
                           if (phy_process->at(f) == 5)
                           Quark_annihlation->Fill(delta_phi);
                           }
                           */
                    }
                    if (mass->at(f) > 1.55 && mass->at(f) < 2.2 && mass->at(g) > 1.55 && mass->at(g) < 2.2)
                    {
                        F1Mass->Fill(mass->at(g));
                        F2Mass->Fill(mass->at(f));
                        FullRangeM1M2Mass->Fill(mass->at(g), mass->at(f));
                    }

                    /*
                       if ((fromD0->at(f) == 0 || fromD0->at(g) == 0))
                       {
                       background_dphi.at(bin-1)->Fill(delta_phi); // background dphi
                       }
                       */

                    //////for dphi correlation below!!!!!!!!!!!!!!!!!!!!//

                    /*
                       cout << "" << endl;
                       cout << "-------------- new particle pair--------------" << endl;
                       cout << "-----------------g is D0 cand-----------------" << endl;
                       cout << "particle f mass = " << mass->at(f) << " and particle g mass = " << mass->at(g) << endl;
                       cout << "delta phi = " << delta_phi << endl;
                       cout << "code determined the dphi correlation component to be " << endl;
                       */
                    if (fisF2 && gisF1)
                    {
                        F1F2->Fill(delta_phi);
                        // cout << "mass f2 = " << mass->at(f) << " and f1 =" << mass->at(g) << endl;
                        //  cout << "~F1F2" << endl;
                    }
                    if (fisSB2 && gisF1)
                    {
                        F1SB2->Fill(delta_phi);
                        N_F1SB2 += 1;
                        // cout << "mass f1 = " << mass->at(g) << " and sb2 =" << mass->at(f) << endl;
                        // cout << "delta_phi = " << endl;
                        //  cout << "~F1SB2" << endl;
                    }
                    if (fisF2 && gisSB1)
                    {
                        SB1F2->Fill(delta_phi);
                        N_SB1F2 += 1;
                        // cout << "mass f2 = " << mass->at(f) << " and sb1 =" << mass->at(g) << endl;
                        //  cout << "~SB1F2" << endl;
                    }
                    if (fisSB2 && gisSB1)
                    {
                        N_SB1SB2 += 1;
                        SB1SB2->Fill(delta_phi);
                        // cout << "mass sb2 = " << mass->at(f) << " and sb1 =" << mass->at(g) << endl;
                        //  cout << "~SB1SB2" << endl;
                        //   SB1Mass->Fill(mass->at(g));
                        //   SB2Mass->Fill(mass->at(f));
                    }
                    if (gisB1 && fisB2)
                    {
                        N_B1B2 += 1;
                        B1B2->Fill(delta_phi);
                        // cout << "mass b2 = " << mass->at(f) << " and b1 =" << mass->at(g) << endl;
                        ////cout << "~B1B2" << endl;
                        // B1Mass->Fill(mass->at(g));
                        // B2Mass->Fill(mass->at(f));
                    }
                    if (gisB1 && fisF2)
                    {
                        B1F2->Fill(delta_phi);
                        N_B1F2 += 1;
                        // cout << "mass f2 = " << mass->at(f) << " and b1 =" << mass->at(g) << endl;
                        //  cout << "~B1F2" << endl;
                    }
                    if (gisF1 && fisB2)
                    {
                        F1B2->Fill(delta_phi);
                        N_F1B2 += 1;
                        // cout << "mass b2 = " << mass->at(f) << " and f1 =" << mass->at(g) << endl;
                        //  cout << "~F1B2" << endl;
                    }

                    // if(gisS1 && fisS2) S1S2->Fill(delta_phi);
                    if (gisS1 && fisB2)
                    {
                        S1B2->Fill(delta_phi);
                        // cout << "~S1B2" << endl;
                    }
                    if (gisB1 && fisS2)
                    {
                        B1S2->Fill(delta_phi);
                        // cout << "~B1S2" << endl;
                    }

                } // D0barcand then D0can

                bin = 0;
            } // g loop

        } // f loop

        /* 
        /////
        if ((N_SB1F2 * N_F1SB2 * N_B1B2) != (N_SB1SB2 * N_B1F2 * N_F1B2))
        {
            {
                cout << "ERROR WITH RATIOS ___ BREAK" << endl;
                cout << "Ifile =" << Ifile << " and event_no =" << event_no << endl;
                cout << "SB1SB2 / B1B2 = " << N_SB1SB2 / N_B1B2 << endl;
                cout << "SB1F2 / B1F2 = " << N_SB1F2 / N_B1F2 << endl;
                cout << "F1SB2 / F1B2 = " << N_F1SB2 / N_F1B2 << endl;
                error_files += 1;
            }
        }
        */
        /*
        else
        {
            cout << "SB1F2->GetEntries()=" << SB1F2->GetEntries() << endl;
            cout << "a =  " << a << " b = " << b << " c = " << c << endl;
        }
        */
    } // for loop

    /*
    TCanvas *c = new TCanvas("c", "", 1200, 800);
    c->Divide(3, 2);
    c->cd(1);
    inclusiveSignal->Draw();
    c->cd(2);
    correlated_dphi->Draw();
    c->cd(3);
    uncorrelated_dphi->Draw();
    c->cd(4);
    background_dphi->Draw();
    c->cd(5);
    nonprompt_phi->Draw();
    c->cd(6);
    auto hs = new THStack("hs", "");
    hs->SetTitle("Dphi by Physics Process");
    Gluon_splitting->SetFillColor(kBlue);
    Gluon_fusion->SetFillColor(kRed);
    flavor_excitation_gluon->SetFillColor(kGreen);
    flavor_excitation_quark->SetFillColor(kCyan);
    Quark_annihlation->SetFillColor(kYellow);
    hs->Add(Quark_annihlation);
    hs->Add(flavor_excitation_quark);
    hs->Add(flavor_excitation_gluon);
    hs->Add(Gluon_fusion);
    hs->Add(Gluon_splitting);
    hs->Draw();
    c->SaveAs("test_histograms_massplots.pdf");

    TCanvas *c1 = new TCanvas("c1", "", 1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    B1B2->Draw();
    c1->cd(2);
    B1Mass->Draw();
    c1->cd(3);
    B2Mass->Draw();
    c1->cd(4);
    SB1SB2->Draw();
    c1->cd(5);
    SB1Mass->Draw();
    c1->cd(6);
    SB2Mass->Draw();
    c1->SaveAs("sidebands_background_check_after_cuts.pdf");
    */

    TCanvas *c2 = new TCanvas("c2", "", 1200, 1200);
    c2->Divide(2, 2);
    c2->cd(1);
    F1F2->Draw();
    c2->cd(2);
    F1SB2->Draw();
    c2->cd(3);
    SB1F2->Draw();
    c2->cd(4);
    SB1SB2->Draw();
    // c2->SaveAs("dphi_closure_test_components.pdf");

    TFile *results = new TFile(outfile, "recreate");
    results->cd();
    // SB1Mass->Write();
    // SB2Mass->Write();
    // SW1Mass->Write();
    // SW2Mass->Write();
    // SBSW1Mass->Write();
    // SBSW2Mass->Write();
    // S1Mass->Write();
    // S2Mass->Write();
    F1Mass->Write();
    F2Mass->Write();
    FullRangeM1M2Mass->Write();
    // S1S2inclusiveSignalMass->Write();
    // B1Mass->Write();
    // B2Mass->Write();
    S1B2->Write();
    B1S2->Write();
    F1F2->Write();
    F1SB2->Write();
    SB1F2->Write();
    SB1SB2->Write();
    B1B2->Write();
    B1F2->Write();
    F1B2->Write();
    // S1S2.at(0)->Add(F1F2.at(0));
    // S1S2.at(0)->Add(SB1SB2.at(0));
    // S1S2.at(0)->Add(F1SB2.at(0), -1);
    // S1S2.at(0)->Add(SB1F2.at(0), -1);
    S1S2->Add(F1F2);
    S1S2->Add(SB1SB2);
    S1S2->Add(F1SB2, -1);
    S1S2->Add(SB1F2, -1);
    S1S2->Write();
    inclusiveSignal->Write();

    TH1F *S1S2e = new TH1F("S1S2e", "S1S2e", nbins, xmin, xmax);
    S1S2e->Sumw2();
    S1S2e->Add(F1F2);
    S1S2e->Add(F1SB2, -1);
    S1S2e->Add(SB1F2, -1);
    S1S2e->Add(SB1SB2);
    S1S2e->Write();

    TCanvas *c3 = new TCanvas("c3", "", 1200, 1200);
    c3->Divide(2, 2);
    c3->cd(1);
    inclusiveSignal->SetMinimum(0);
    inclusiveSignal->SetMaximum(30000);
    inclusiveSignal->SetLineColor(kBlack);
    inclusiveSignal->Draw();
    c3->cd(2);
    S1S2e->SetMinimum(0);
    S1S2e->SetTitle("S1S2 via SB Method");
    S1S2e->SetLineColor(kRed);
    S1S2e->SetMaximum(30000);
    inclusiveSignal->SetMaximum(30000);
    S1S2e->Draw();
    c3->cd(3);

    TH1F *hist = new TH1F("hist", "S1S2 Template Fit", nbins, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    double values[5] = {21780.5, 23632.3, 22999.0, 26579.6, 23930.4};
    double errors[5] = {465.707, 440.069, 431.497, 536.060, 560.55};

    for (int i = 1; i < 6; i++)
    {
        hist->SetBinContent(i, values[i - 1]);
        hist->SetBinError(i, errors[i - 1]);
    }
    hist->SetMinimum(0);
    hist->SetMaximum(30000);
    hist->SetLineColor(kBlue);
    hist->Draw("E");
    c3->cd(4);
    inclusiveSignal->Draw();
    S1S2e->Draw("same");
    hist->Draw("same");

    cout << "---" << endl;
    for (int i = 1; i < 6; i++)
    {
        cout << "for bin " << i << " fraction of S1S2 / inclusive signal = " << S1S2->GetBinContent(i) / inclusiveSignal->GetBinContent(i) << endl;
        cout << "---" << endl;
    }
    cout << "overall closure = " << S1S2->GetEntries() << "/" << inclusiveSignal->GetEntries() << endl;
    cout << "overall closure = " << S1S2->GetEntries() / inclusiveSignal->GetEntries() << endl;
    cout << "---" << endl;
    /*
    correlated_dphi->Write();
    uncorrelated_dphi->Write();
    background_dphi->Write();
    nonprompt_phi->Write();
    Gluon_splitting->Write();
    Gluon_fusion->Write();
    flavor_excitation_gluon->Write();
    flavor_excitation_quark->Write();
    Quark_annihlation->Write();
    */
    cout << "VVVVVVVV error files = " << error_files << " VVVVVVVV" << endl;
    inf1->Close();
    infile->Close();
    // results->Close();

    return 0;
}

int main(int argc, char* argv[]) {
        if (argc != 5) {
                std::cerr << "Usage: " << argv[0] << " <goodbin> <prob_cutoff> <outfile> <limits_path>" << std::endl;
                return 1;
        }
        float prob_cutoff = std::atof(argv[2]);
        TString outfile = argv[3];
        TString limits_path = argv[4];
        int goodbin = std::atoi(argv[1]);

        temp_5bins(goodbin, prob_cutoff, outfile, limits_path);
        return 0;
}
