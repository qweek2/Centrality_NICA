#include <fstream>
#include <iostream>
#include <TNtuple.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <Rtypes.h>
#include <TBrowser.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMath.h>
#include <TKey.h>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TCutG.h"
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TF2.h>
#include <vector>

//git version

Double_t g2(Double_t *x, Double_t *par)
{
    Double_t r1 = Double_t((x[0] - par[1]) / par[2]);
    Double_t r2 = Double_t((x[1] - par[3]) / par[4]);
    return par[0] * TMath::Exp(-0.5 * (r1 * r1 + r2 * r2));
}

Double_t fun2(Double_t *x, Double_t *par) //for signal
{
    Double_t z2 = x[0] * x[0] / par[0] / par[0] + x[1] * x[1] / par[0] / par[0];
    z2 *= par[1] * par[1];
    return -sqrt(z2) + par[2];
}

/*Double_t fun2(Double_t *x, Double_t *par) //for signal
{
	Double_t z2 = x[0] * x[0] / par[0] / par[0] + x[1] * x[1] / par[1] / par[1];
	z2 *= par[2] * par[2];
	return -sqrt(z2) + par[3];
}*/

// c/a=ctg(φ)
Double_t fun3(Double_t *x, Double_t *par) //for pions
{
    //0.000142143 (result in radians) // 0.008144°
    double phi = 0.000142143;
    Double_t z2 = x[0] * x[0] / ((par[0] * tan(phi)) * (par[0] * tan(phi))) + x[1] * x[1] / ((par[0] * tan(phi)) * (par[0] * tan(phi)));
    z2 *= par[0] * par[0];
    return -sqrt(z2) + par[1];
}

//alternative
Double_t fun4(Double_t *x, Double_t *par) //for pions 2
{
    //0.000142143 (result in radians) // 0.008144°
    //double phi = 0.000142143;
    //t = cos(phi)/sin(phi);
    //f->FixParameter(1, 2.); // to try
    Double_t z2 = (x[0] * x[0] * par[0] * par[0]) / (par[1] * par[1]) + (x[1] * x[1] * par[0] * par[0]) / (par[1] * par[1]);
    z2 *= par[1] * par[1];
    return -sqrt(z2) + par[2];
}

void twod_gauss_energy()
{
    ofstream myfile;
    myfile.open("SMM_11_Edep_Emax.txt");
    TGraph *g = new TGraph();  // using the blank constructor
    TGraph *g1 = new TGraph(); // using the blank constructor
    TGraph *g2 = new TGraph(); // using the blank constructor
    TGraph *g3 = new TGraph(); // using the blank constructor
    TGraph *g4 = new TGraph(); // using the blank constructor
    TGraph *g5 = new TGraph(); // using the blank constructor
    TGraph *g6 = new TGraph(); // using the blank constructor

    Double_t edep_7sect_1, edep_7sect_2;
    Double_t edepMod_1[91], edepMod_2[91], edepMod_3[91], edepMod_4[91], edepMod_5[91], edepMod_6[91], edepMod_7[91];
    Double_t edepMod_1_p[91], edepMod_2_p[91], edepMod_3_p[91], edepMod_4_p[91], edepMod_5_p[91], edepMod_6_p[91], edepMod_7_p[91];
    Double_t impPar, impPar_p;

    //TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_full_geom_QGSM_AuAu_11_mb_500_99500ev.root");
    //TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_full_geom_SMM_AuAu_MB_s5GeV_1_2_98977ev.root");
    //TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_full_geom_QGSM_AuAuss5mb_99500ev.root");
    TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_full_geom_SMM_AuAu_MB_s11GeV_1_2_3_4_5_99483ev.root");
    //analzdc_full_geom_SMM_AuAu_MB_s11GeV_1_2_3_4_5_99483ev.root
    TTree *nt1 = (TTree *)_file0->Get("nt1");
    nt1->SetBranchAddress("edepMod_1", &edepMod_1);
    nt1->SetBranchAddress("edepMod_2", &edepMod_2);
    nt1->SetBranchAddress("edepMod_3", &edepMod_3);
    nt1->SetBranchAddress("edepMod_4", &edepMod_4);
    nt1->SetBranchAddress("edepMod_5", &edepMod_5);
    nt1->SetBranchAddress("edepMod_6", &edepMod_6);
    nt1->SetBranchAddress("edepMod_7", &edepMod_7);
    nt1->SetBranchAddress("impPar", &impPar);
    nt1->SetBranchAddress("edep_7sect_1", &edep_7sect_1);
    nt1->SetBranchAddress("edep_7sect_2", &edep_7sect_2);

    TFile *_file1 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_QGSM_AuAu_11_mb_only_pions_98500ev.root"); //pion
    TTree *nt1_pion = (TTree *)_file1->Get("nt1");

    nt1_pion->SetBranchAddress("edepMod_1", &edepMod_1_p);
    nt1_pion->SetBranchAddress("edepMod_2", &edepMod_2_p);
    nt1_pion->SetBranchAddress("edepMod_3", &edepMod_3_p);
    nt1_pion->SetBranchAddress("edepMod_4", &edepMod_4_p);
    nt1_pion->SetBranchAddress("edepMod_5", &edepMod_5_p);
    nt1_pion->SetBranchAddress("edepMod_6", &edepMod_6_p);
    nt1_pion->SetBranchAddress("edepMod_7", &edepMod_7_p);
    nt1_pion->SetBranchAddress("impPar", &impPar_p);

    TString result_filename = "Edep_Emax_SMM_min.root";
    TFile *result_file;
    result_file = TFile::Open(result_filename, "RECREATE");

    const int nbCBMmods = 90; //NICA 45+45 mods (VETO1 - 1-45, VETO 2 -46-89)
    double xxxxNICA[nbCBMmods], yyyyNICA[nbCBMmods];

    Int_t modNb = -1;
    double xCur = 60;
    double yCur = 45;
    for (Int_t iii = 0; iii < 5; iii++)
    {
        for (Int_t jjj = 0; jjj < 7; jjj++)
        {
            modNb = 7 * iii + jjj + 5;
            xxxxNICA[modNb] = xCur - (jjj + 1) * 15.;
            yyyyNICA[modNb] = yCur - (iii + 1) * 15.;
        }
    }

    xCur = 45;
    for (Int_t jjj = 0; jjj < 2; jjj++)
    {
        for (Int_t iii = 0; iii < 5; iii++)
        {
            if (jjj == 0)
            {
                modNb = iii;
                yyyyNICA[modNb] = yCur;
            }
            if (jjj == 1)
            {
                modNb = 40 + iii;
                yyyyNICA[modNb] = -yCur;
            }
            xxxxNICA[modNb] = xCur - (iii + 1) * 15.;
            //cout << "xx " << iii << " " << jjj << " " << modNb + 1 << " " << xxxxNICA[modNb] << " " << yyyyNICA[modNb] << endl;
        }
    }

    //NICA 45 mods, 15x15 hole 10cm in mod 68 VETO 2
    for (Int_t jjj = 0; jjj < nbCBMmods / 2; jjj++)
    {
        xxxxNICA[jjj + 45] = -xxxxNICA[jjj];
        //xxxxNICA[jjj+45] = xxxxNICA[jjj];
        yyyyNICA[jjj + 45] = yyyyNICA[jjj];
    }

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    TH1 *th1_hist_ptr = NULL;
    TH2 *th2_hist_ptr = NULL;
    Int_t module_binsX = 7;
    Int_t module_binsY = 7;
    Float_t max_energy = 50;
    Int_t energy_bins = 200;
    Float_t max_imp_par = 15.;
    Int_t imp_par_bins = 200;

    vector<float> ring_radius;
    for (Int_t i = 1; i < 4; i++)
    {
        for (Int_t j = 0; j < 4; j++)
        {
            float radius = sqrt(pow(15 * i, 2) + pow(15 * j, 2));
            bool IsNew = true;
            for (uint8_t k = 0; k < ring_radius.size(); k++)
                if (abs(radius - ring_radius.at(k)) < 0.1)
                {
                    IsNew = false;
                    break;
                }
            if (IsNew)
                ring_radius.push_back(radius);
        }
    }
    std::sort(ring_radius.begin(), ring_radius.end());

    int en_rings_nmbr = ring_radius.size();
    for (int i = 0; i < en_rings_nmbr; i++)
        cout << " " << ring_radius.at(i) << endl;
    float energy_in_ring[en_rings_nmbr];
    int entries_in_ring[en_rings_nmbr];

    TH2F *EdepEcalc = new TH2F("EdepEcalc", "EdepEcalc", 1000, 0., 20., 1000, 0, 20.);

    const Int_t npar = 4;
    const Int_t npar2 = 2;
    const Int_t npar3 = 3;

    TF2 *f2 = new TF2("f2", fun2, -52.5, 52.5, -52.5, 52.5, npar);
    TF2 *f5 = new TF2("f5", fun2, -52.5, 52.5, -52.5, 52.5, npar);
    TF2 *f_diff = new TF2("f_diff", fun2, -52.5, 52.5, -52.5, 52.5, npar);
    TF2 *f6 = new TF2("f6", fun2, -52.5, 52.5, -52.5, 52.5, npar3);
    TF2 *f7 = new TF2("f7", fun2, -52.5, 52.5, -52.5, 52.5, npar3);
    TF2 *f8 = new TF2("f8", fun2, -52.5, 52.5, -52.5, 52.5, npar3);
    TF2 *f_p = new TF2("f_p", fun2, -52.5, 52.5, -52.5, 52.5, npar3);
    TF2 *f_p2 = new TF2("f_p2", fun4, -52.5, 52.5, -52.5, 52.5, npar3);
    f_p2->SetParLimits(0, 0.4, 100000); //set range for par[0] = c/a=cos/sin 1000 to 100000 default

    //TF2 *f2 = new TF2("f2", fun2, -90, 90, -70, 70, npar);

    th1_hist_ptr = new TH1F("Impact_parameter", "Impact_parameter", 500, 0, 20);

    th2_hist_ptr = new TH2F("PSD_energy_surface", "PSD_energy_surface", module_binsX, -52.5, 52.5, module_binsY, -52.5, 52.5);           //+
    th2_hist_ptr = new TH2F("PSD_energy_surface_pion", "PSD_energy_surface_pion", module_binsX, -52.5, 52.5, module_binsY, -52.5, 52.5); //pion
    TH2F *hDiff = new TH2F("hDiff", "hDiff", module_binsX, -52.5, 52.5, module_binsY, -52.5, 52.5);                                      //+
                                                                                                                                         //+
    TH2F *hPionsFit = new TH2F("hPionsFit", "hPionsFit", module_binsX, -52.5, 52.5, module_binsY, -52.5, 52.5);                          //+
    TH2F *hFin = new TH2F("hFin", "hFin", module_binsX, -52.5, 52.5, module_binsY, -52.5, 52.5);                                         //+
    TH2F *hFin2 = new TH2F("hFin2", "hFin2", module_binsX, -52.5, 52.5, module_binsY, -52.5, 52.5);                                      //+
    TH2F *hEnergy = new TH2F("hEnergy", "hEnergy", module_binsX, -52.5, 52.5, module_binsY, -52.5, 52.5);                                //+
    TH2F *hImpAng = new TH2F("hImpAng", "hImpAng", 500, 0, 20, 500, 0, 20);                                                              //+
    TH2F *hImp_E_p = new TH2F("hImpEp", "hImpEp", 500, 0, 25, 500, 0, 25);                                                               //+
    TH2F *hRadEmax = new TH2F("hRadEmax", "hRadEmax", 2000, 0, 100, 2000, 0, 3.5);                                                       //+
    TH2F *hEdepEmax = new TH2F("EdepEmax", "EdepEmax", 1000, 0, 1, 1000, 0, 1);                                                          //+
    TH2F *hPions = new TH2F("hPions", "hPions", 500, 0, 16, 500, 0, 10);                                                                 //+

    th1_hist_ptr = new TH1F("PSD_energy_in_ring", "PSD_energy_in_ring", 100, 0, ring_radius.back() + 1);
    th2_hist_ptr = new TH2F("PSD_total_energy_vs_imp_par", "PSD_total_energy_vs_imp_par", imp_par_bins, 0, max_imp_par, energy_bins, 0, max_energy);
    th2_hist_ptr = new TH2F("PSD_surface_integral_vs_imp_par", "PSD_surface_integral_vs_imp_par", imp_par_bins, 0, max_imp_par, energy_bins, 0, max_energy);
    th2_hist_ptr = new TH2F("PSD_surface_sigma_vs_imp_par", "PSD_surface_sigma_vs_imp_par", imp_par_bins, 0, max_imp_par, 200, 0, 100);
    th2_hist_ptr = new TH2F("PSD_surface_sigma_vs_assymetry", "PSD_surface_sigma_vs_assymetry", 200, -1, 1, 200, 0, 100);

    time_t start_time = time(NULL);
    printf("Filling initial histos\n");

    ///////////////////////////////////
    ///////////////////////////////////
    /////////////  data  //////////////
    ///////////////////////////////////
    ///////////////////////////////////
    for (Int_t entry = 0; entry < 50000; entry++) //nt1->GetEntries() // 5 vs 8 1700 bad
    {
        th2_hist_ptr->Reset();
        th2_hist_ptr = ((TH2 *)(gDirectory->FindObjectAny("PSD_energy_surface")));

        hDiff->Reset();
        hEnergy->Reset();

        //hPionsFit->Reset();
        for (int iBin = 1; iBin < th2_hist_ptr->GetNbinsX() * th2_hist_ptr->GetNbinsY(); iBin++)
        {
            if (th2_hist_ptr->GetBinContent(iBin) > 0)
                th2_hist_ptr->SetBinError(iBin, 1. / sqrt(th2_hist_ptr->GetBinContent(iBin)));
            else
                th2_hist_ptr->SetBinError(iBin, numeric_limits<double>::max());
        }
        for (Int_t ring_iter = 0; ring_iter < en_rings_nmbr; ring_iter++)
        {
            energy_in_ring[ring_iter] = 0.;
            entries_in_ring[ring_iter] = 0;
        }
        nt1->GetEntry(entry);

        if ((entry % 10) == 0)
        {
            time_t current_time = time(NULL);
            Int_t proc_sec = difftime(current_time, start_time);
            Float_t percents = (float)entry / (nt1->GetEntries());
            Int_t time_est = (percents == 0) ? 0 : proc_sec / percents * (1. - percents);
            Float_t proc_rate = (float)proc_sec / entry * 1000000. / 60.;

            printf("Processed events: %i (%5.1f%%); [pas %3.0dm %2.0is] [est %3.0dm %2.0is] [%.1f min [%.1fh]/1M ev]\r",
                   entry, (percents * 100.),
                   (proc_sec / 60), proc_sec % 60,
                   (time_est / 60), time_est % 60,
                   proc_rate, (proc_rate / 60.));
            cout << flush;
        }

        for (Int_t m_i = 0; m_i < 90; m_i++)
        //NICA 45+45 mods (VETO1 - 1-45 (from 0 to 44), VETO 2 -46-89 (from 45 to 89))
        {
            if (m_i == 22 || m_i == 67)
                continue;
            Float_t Xposition = xxxxNICA[m_i];
            Float_t Yposition = yyyyNICA[m_i];
            Double_t Energy = edepMod_7[m_i + 1];

            th2_hist_ptr = ((TH2 *)(gDirectory->FindObjectAny("PSD_energy_surface")));
            th2_hist_ptr->Fill(Xposition, Yposition, Energy * 1000);
            // hist with data
            hDiff->Fill(Xposition, Yposition, 1000 * Energy);
            ////////////////////////////
            //////uniformation begin////
            ////////////////////////////
            if (m_i == 15 || m_i == 23 || m_i == 21 || m_i == 29 || m_i == 15 + 45 || m_i == 23 + 45 || m_i == 21 + 45 || m_i == 29 + 45) //1
            {
                hEnergy->Fill(xxxxNICA[15], yyyyNICA[15], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[21], yyyyNICA[21], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[23], yyyyNICA[23], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[29], yyyyNICA[29], 1000 * Energy / 4);
            }
            else if (m_i == 16 || m_i == 14 || m_i == 30 || m_i == 28 || m_i == 16 + 45 || m_i == 14 + 45 || m_i == 30 + 45 || m_i == 28 + 45) //2
            {
                hEnergy->Fill(xxxxNICA[16], yyyyNICA[16], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[14], yyyyNICA[14], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[30], yyyyNICA[30], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[28], yyyyNICA[28], 1000 * Energy / 4);
            }
            else if (m_i == 8 || m_i == 20 || m_i == 24 || m_i == 36 || m_i == 8 + 45 || m_i == 20 + 45 || m_i == 24 + 45 || m_i == 36 + 45) //3
            {
                hEnergy->Fill(xxxxNICA[8], yyyyNICA[8], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[20], yyyyNICA[20], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[36], yyyyNICA[36], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[24], yyyyNICA[24], 1000 * Energy / 4);
            }
            else if (m_i == 9 || m_i == 7 || m_i == 17 || m_i == 13 || m_i == 31 || m_i == 27 || m_i == 37 || m_i == 35 || m_i == 9 + 45 || m_i == 7 + 45 || m_i == 17 + 45 || m_i == 13 + 45 || m_i == 31 + 45 || m_i == 27 + 45 || m_i == 37 + 45 || m_i == 35 + 45)
            { //4
                hEnergy->Fill(xxxxNICA[9], yyyyNICA[9], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[7], yyyyNICA[7], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[17], yyyyNICA[17], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[13], yyyyNICA[13], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[31], yyyyNICA[31], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[27], yyyyNICA[27], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[37], yyyyNICA[37], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[35], yyyyNICA[35], 1000 * Energy / 8);
            }
            else if (m_i == 10 || m_i == 6 || m_i == 34 || m_i == 38 || m_i == 10 + 45 || m_i == 6 + 45 || m_i == 34 + 45 || m_i == 38 + 45) //5
            {
                hEnergy->Fill(xxxxNICA[10], yyyyNICA[10], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[6], yyyyNICA[6], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[34], yyyyNICA[34], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[38], yyyyNICA[38], 1000 * Energy / 4);
            }
            else if (m_i == 2 || m_i == 25 || m_i == 19 || m_i == 42 || m_i == 2 + 45 || m_i == 25 + 45 || m_i == 19 + 45 || m_i == 42 + 45) //6
            {
                hEnergy->Fill(xxxxNICA[25], yyyyNICA[25], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[2], yyyyNICA[2], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[19], yyyyNICA[19], 1000 * Energy / 4);
                hEnergy->Fill(xxxxNICA[42], yyyyNICA[42], 1000 * Energy / 4);
            }
            else if (m_i == 3 || m_i == 1 || m_i == 18 || m_i == 12 || m_i == 32 || m_i == 26 || m_i == 43 || m_i == 41 || m_i == 3 + 45 || m_i == 1 + 45 || m_i == 18 + 45 || m_i == 12 + 45 || m_i == 32 + 45 || m_i == 26 + 45 || m_i == 43 + 45 || m_i == 41 + 45)
            { //7
                hEnergy->Fill(xxxxNICA[3], yyyyNICA[3], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[1], yyyyNICA[1], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[18], yyyyNICA[18], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[12], yyyyNICA[12], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[32], yyyyNICA[32], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[26], yyyyNICA[26], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[43], yyyyNICA[43], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[41], yyyyNICA[41], 1000 * Energy / 8);
            }
            else if (m_i == 44 || m_i == 40 || m_i == 11 || m_i == 5 || m_i == 4 || m_i == 0 || m_i == 39 || m_i == 33 || m_i == 44 + 45 || m_i == 40 + 45 || m_i == 11 + 45 || m_i == 5 + 45 || m_i == 4 + 45 || m_i == 0 + 45 || m_i == 39 + 45 || m_i == 33 + 45)
            { //8 last point
                hEnergy->Fill(xxxxNICA[40], yyyyNICA[40], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[44], yyyyNICA[44], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[11], yyyyNICA[11], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[5], yyyyNICA[5], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[4], yyyyNICA[4], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[0], yyyyNICA[0], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[39], yyyyNICA[39], 1000 * Energy / 8);
                hEnergy->Fill(xxxxNICA[33], yyyyNICA[33], 1000 * Energy / 8);
            }

            ////////////////////////////
            //////uniformation end//////
            ////////////////////////////
        }

        /*double E_sum = hDiff->Integral();
		for (Int_t module_iter = 0; module_iter < 90; module_iter++) //NICA 45+45 mods (VETO1 - 1-45 (from 0 to 44), VETO 2 -46-89 (from 45 to 89))
		{
			if (module_iter == 22 || module_iter == 67)
				continue;
			Float_t Xposition = xxxxNICA[module_iter];
			Float_t Yposition = yyyyNICA[module_iter];
			Double_t Energy = edepMod_7[module_iter + 1];
			// hist with normalized data
			hDiff->Fill(Xposition, Yposition, 100 * Energy/E_sum);
		}*/

        th1_hist_ptr = ((TH1 *)(gDirectory->FindObjectAny("PSD_energy_in_ring")));
        th1_hist_ptr->Reset();
        for (int i = 0; i < en_rings_nmbr; i++)
            if (entries_in_ring[i] > 0)
                th1_hist_ptr->Fill(ring_radius.at(i), energy_in_ring[i] / entries_in_ring[i]);

        //weights for signal
        for (int iBin = 1; iBin <= hEnergy->GetNbinsX() * hEnergy->GetNbinsY(); iBin++)
        {
            if (hEnergy->GetBinContent(iBin) > 0)
                hEnergy->SetBinError(iBin, 1. / sqrt(hEnergy->GetBinContent(iBin))); //true from physics
                                                                                     //hEnergy->SetBinError(iBin, 1/hEnergy->GetBinContent(iBin));
            else
                hEnergy->SetBinError(iBin, numeric_limits<double>::max());
        }

        //weights for pions
        for (int iBin = 1; iBin <= hPionsFit->GetNbinsX() * hPionsFit->GetNbinsY(); iBin++)
        {
            if (hPionsFit->GetBinContent(iBin) > 0)
                hPionsFit->SetBinError(iBin, 1. / sqrt(hPionsFit->GetBinContent(iBin)));
            else
                hPionsFit->SetBinError(iBin, numeric_limits<double>::max());
        }
        /*cout << "" << endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!! next fit  !!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "" << endl;*/
        Double_t f6params[npar3] = {35, 1000, 1000};
        f6->SetParameters(f6params);
        hEnergy->Fit("f6", "QR M"); //VR - more words
        //cout << f6->Eval(0, 0) << endl;

        hFin->Reset();
        //double Ecalc = 0; // = hFin->Integral();
        for (int i = 1; i < hFin->GetXaxis()->GetNbins() + 1; i++)
        {
            for (int j = 1; j < hFin->GetYaxis()->GetNbins() + 1; j++)
            {
                if (i == 4 && j == 4)
                    continue;
                double x = hFin->GetXaxis()->GetBinCenter(i);
                double y = hFin->GetYaxis()->GetBinCenter(j);
                double z = f6->Eval(x, y);
                if (z < 0) // ??? due to conic equation or fit
                    z = 0;
                //Ecalc += z;
                hFin->Fill(x, y, z);
            }
        }
        //weights for signal
        for (int iBin = 1; iBin <= hFin->GetNbinsX() * hFin->GetNbinsY(); iBin++)
        {
            if (hFin->GetBinContent(iBin) > 0)
                hFin->SetBinError(iBin, 1. / sqrt(hFin->GetBinContent(iBin))); //true from physics
                                                                               //hFin->SetBinError(iBin, 1/hFin->GetBinContent(iBin));
            else
                hFin->SetBinError(iBin, numeric_limits<double>::max());
        }
        Double_t f7params[npar3] = {35, 1000, 1000};
        f7->SetParameters(f7params);
        hFin->Fit("f7", "QR M");
        hPionsFit->Reset();
        //pions from fit
        /*double pion_energy = f6->Eval(xxxxNICA[0], yyyyNICA[0]) / 4;
		pion_energy += f6->Eval(xxxxNICA[4], yyyyNICA[4]) / 4;
		pion_energy += f6->Eval(xxxxNICA[40], yyyyNICA[40]) / 4;
		pion_energy += f6->Eval(xxxxNICA[44], yyyyNICA[44]) / 4;*/

        //pions from hist
        double pion_energy = 0;
        pion_energy = hEnergy->GetBinContent(2, 1) / 8;
        pion_energy += hEnergy->GetBinContent(6, 1) / 8;
        pion_energy += hEnergy->GetBinContent(1, 2) / 8;
        pion_energy += hEnergy->GetBinContent(7, 2) / 8;
        pion_energy += hEnergy->GetBinContent(1, 6) / 8;
        pion_energy += hEnergy->GetBinContent(7, 6) / 8;
        pion_energy += hEnergy->GetBinContent(2, 7) / 8;
        pion_energy += hEnergy->GetBinContent(6, 7) / 8;
        //myfile << -1.5 * pion_energy / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) << endl;
        for (Int_t module_iter = 0; module_iter < 45; module_iter++)
        {
            if (module_iter == 22 || module_iter == 67)
                continue;
            //myfile << entry << " " <<  f7->Eval(xxxxNICA[0], yyyyNICA[0]) << endl;
            if (pion_energy >= 0)
            {
                if ((module_iter == 15 || module_iter == 23 || module_iter == 21 || module_iter == 29))                                                                                                                                                                                                                                                                                                                     //1
                    hPionsFit->Fill(xxxxNICA[module_iter], yyyyNICA[module_iter], 1.5 * pion_energy * (sqrt(xxxxNICA[module_iter] * xxxxNICA[module_iter] + yyyyNICA[module_iter] * yyyyNICA[module_iter]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) + pion_energy); //1.5 * pion_energy) * (sqrt(xxxxNICA[module_iter] * xxxxNICA[module_iter] + yyyyNICA[module_iter] * yyyyNICA[module_iter]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) + pion_energy));

                else if ((module_iter == 16 || module_iter == 14 || module_iter == 30 || module_iter == 28)) //2
                    hPionsFit->Fill(xxxxNICA[module_iter], yyyyNICA[module_iter], 1.5 * pion_energy * (sqrt(xxxxNICA[module_iter] * xxxxNICA[module_iter] + yyyyNICA[module_iter] * yyyyNICA[module_iter]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) + pion_energy);

                else if ((module_iter == 8 || module_iter == 20 || module_iter == 24 || module_iter == 36)) //3
                    hPionsFit->Fill(xxxxNICA[module_iter], yyyyNICA[module_iter], 1.5 * pion_energy * (sqrt(xxxxNICA[module_iter] * xxxxNICA[module_iter] + yyyyNICA[module_iter] * yyyyNICA[module_iter]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) + pion_energy);

                else if ((module_iter == 9 || module_iter == 7 || module_iter == 17 || module_iter == 13 || module_iter == 31 || module_iter == 27 || module_iter == 37 || module_iter == 35)) //4
                    hPionsFit->Fill(xxxxNICA[module_iter], yyyyNICA[module_iter], 1.5 * pion_energy * (sqrt(xxxxNICA[module_iter] * xxxxNICA[module_iter] + yyyyNICA[module_iter] * yyyyNICA[module_iter]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) + pion_energy);

                else if ((module_iter == 10 || module_iter == 6 || module_iter == 34 || module_iter == 38)) //5
                    hPionsFit->Fill(xxxxNICA[module_iter], yyyyNICA[module_iter], 1.5 * pion_energy * (sqrt(xxxxNICA[module_iter] * xxxxNICA[module_iter] + yyyyNICA[module_iter] * yyyyNICA[module_iter]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) + pion_energy);

                else if ((module_iter == 2 || module_iter == 25 || module_iter == 19 || module_iter == 42)) //6
                    hPionsFit->Fill(xxxxNICA[module_iter], yyyyNICA[module_iter], 1.5 * pion_energy * (sqrt(xxxxNICA[module_iter] * xxxxNICA[module_iter] + yyyyNICA[module_iter] * yyyyNICA[module_iter]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) + pion_energy);

                else if ((module_iter == 3 || module_iter == 1 || module_iter == 18 || module_iter == 12 || module_iter == 32 || module_iter == 26 || module_iter == 43 || module_iter == 41))
                    hPionsFit->Fill(xxxxNICA[module_iter], yyyyNICA[module_iter], 1.5 * pion_energy * (sqrt(xxxxNICA[module_iter] * xxxxNICA[module_iter] + yyyyNICA[module_iter] * yyyyNICA[module_iter]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) + pion_energy);

                else if ((module_iter == 44 || module_iter == 40 || module_iter == 11 || module_iter == 5 || module_iter == 4 || module_iter == 0 || module_iter == 39 || module_iter == 33))
                    //8 last point
                    hPionsFit->Fill(xxxxNICA[module_iter], yyyyNICA[module_iter], 1.5 * pion_energy * (sqrt(xxxxNICA[module_iter] * xxxxNICA[module_iter] + yyyyNICA[module_iter] * yyyyNICA[module_iter]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) / (sqrt(xxxxNICA[15] * xxxxNICA[15] + yyyyNICA[15] * yyyyNICA[15]) - sqrt(xxxxNICA[0] * xxxxNICA[0] + yyyyNICA[0] * yyyyNICA[0])) + pion_energy);
            }
        }

        /*if (atan(f7->GetParameter(2) / f7->GetParameter(0)) * 180 / 3.1415 > 0)
		{
			//myfile << atan(f7->GetParameter(2) / f7->GetParameter(0)) * 180 / 3.1415 << " " << f7->GetParameter(2) << endl;
			g->SetMarkerColor(1);
			g->SetPoint(g->GetN(), , edep_7sect_1 + edep_7sect_2);
		}*/

        Double_t fpparams[npar3] = {50, 500, 500};
        f_p->SetParameters(fpparams);
        hPionsFit->Fit("f_p", "QR");

        ///////////////////////
        ////        ////    ///
        /////     ///////  ////
        /////     //////   ////
        /////     //////   ////
        /////     /////    ////
        /////      ///  ///////
        ///////////////////////
        ///////////////////////
        //hFin->Add(hPionsFit, -1.);
        //fit w/o pions
        for (int iBin = 1; iBin <= hFin->GetNbinsX() * hFin->GetNbinsY(); iBin++)
        {
            if (hFin->GetBinContent(iBin) > 0)
                hFin->SetBinError(iBin, 1. / sqrt(hFin->GetBinContent(iBin))); //true from physics
                                                                               //hFin->SetBinError(iBin, 1/hFin->GetBinContent(iBin));
            else
                hFin->SetBinError(iBin, numeric_limits<double>::max());
        }
        f8->SetParameters(f7params);
        /*f8->SetParLimits(0, 0, 100);
        f8->SetParLimits(1, 0, 100);
        f8->SetParLimits(2, 0, 5000);*/

        hFin->Fit("f8", "QR M");
        double par2_max, edep_max;
        if (f8->GetParameter(2) > par2_max)
            par2_max = f8->GetParameter(2);
        if (edep_7sect_1 + edep_7sect_2 > edep_max)
            edep_max = edep_7sect_1 + edep_7sect_2;
        if (entry == 49999)
            cout << par2_max / 1000 << "  " << edep_max << endl;

        /*hFin2->Reset();
        for (int i = 1; i < hFin2->GetXaxis()->GetNbins() + 1; i++)
        {
            for (int j = 1; j < hFin2->GetYaxis()->GetNbins() + 1; j++)
            {
                if (i == 4 && j == 4)
                    continue;
                double x = hFin2->GetXaxis()->GetBinCenter(i);
                double y = hFin2->GetYaxis()->GetBinCenter(j);
                double z = f8->Eval(x, y);
                if (z < 0) // ??? due to conic equation
                    z = 0;
                //Ecalc += z;
                hFin2->Fill(x, y, z);
            }
        }*/

        double radius_con = f8->GetParameter(2) * f8->GetParameter(0) / f8->GetParameter(1);
        double Ecalc = (3.1415 * radius_con * radius_con * f8->GetParameter(2)) / 675; //из геометрии V=1/3pi*h*r^2
        double radius_con_p = f_p->GetParameter(2) * f_p->GetParameter(0) / f_p->GetParameter(1);
        double E_p = (3.1415 * radius_con_p * radius_con_p * f_p->GetParameter(2)) / 675;
        //double radius_con = f6->GetParameter(2) * f6->GetParameter(0) / f6->GetParameter(1);
        //myfile << "rad " << radius_con << " par 3 " << f8->GetParameter(2) << " par 0 " << f8->GetParameter(0) << " par 2 " << f8->GetParameter(2) << endl;
        //myfile << atan(f8->GetParameter(2) / f8->GetParameter(0)) * 180 / 3.1415 << " " << f8->GetParameter(2) << endl;
        g->SetMarkerColor(1);
        g1->SetMarkerColor(1);
        g2->SetMarkerColor(1);
        g3->SetMarkerColor(1);
        g4->SetMarkerColor(1);

        if (radius_con > 0 && radius_con < 100 && f8->GetChisquare() < 1000)
            g->SetPoint(g->GetN(), atan(radius_con / 320) * 180 / 3.1415, f8->GetParameter(2) / 1000); // angle vs Emax

        if (radius_con > 0 && radius_con < 100 && f8->GetChisquare() < 1000)
        {
            g1->SetPoint(g1->GetN(), radius_con, f8->GetParameter(2) / 1000); //rad vs Emax
            hRadEmax->Fill(radius_con, f8->GetParameter(2) / 1000);
            //myfile << radius_con << " " << f8->GetParameter(2) / 1000 << " " << impPar << endl;
        }
        //if (radius_con > 0 && radius_con < 100 && f8->GetChisquare() < 1000)
        g2->SetPoint(g2->GetN(), impPar, f8->GetParameter(2) / 1000); // impact vs emax

        if (radius_con > 0 /*&& radius_con < 100 && f8->GetChisquare() < 6000*/)
        {
            //g3->SetPoint(g3->GetN(), f8->GetParameter(2) / 1000, edep_7sect_1 + edep_7sect_2); // emax vs edep ORIGINAL
            g3->SetPoint(g3->GetN(), f8->GetParameter(2) / 1000, edep_7sect_1 + edep_7sect_2); // emax vs edep
            //hEdepEmax->Fill(f8->GetParameter(2) / 1000, edep_7sect_1 + edep_7sect_2);          //free scale
            //myfile << f8->GetParameter(2) / 1000 << " " << edep_7sect_1 + edep_7sect_2 << " " << impPar << endl;
            hEdepEmax->Fill(f8->GetParameter(2) / 1000 / 2.2832, (edep_7sect_1 + edep_7sect_2) / 18.5943); //scale 1:1
            myfile << f8->GetParameter(2) / 1000 / 2.2832 << " " << (edep_7sect_1 + edep_7sect_2) / 18.5943 << " " << impPar << endl;
        }
        if (radius_con > 0 && Ecalc > 0 && Ecalc < 30000 && f8->GetParameter(2) > 0 && f8->GetParameter(2) < 30000)
        {
            g4->SetPoint(g4->GetN(), f8->GetParameter(2) / 1000, Ecalc / 1000); // emax vs ecalc
            hImp_E_p->Fill(Ecalc / 1000, E_p / 1000);
        }

        if (radius_con > 0 && radius_con < 100 && f8->GetChisquare() < 1000)
            g5->SetPoint(g5->GetN(), atan(radius_con / 320) * 180 / 3.1415, edep_7sect_1 + edep_7sect_2); // angle vs edep

        if (true)
        {
            g6->SetPoint(g6->GetN(), impPar, atan(radius_con / 320) * 180 / 3.1415); // impact vs angle
            hImpAng->Fill(impPar, atan(radius_con / 320) * 180 / 3.1415);
        }

        hPions->Fill(impPar, hPionsFit->Integral() / 1000); //imp vs E_pions
        //cout << impPar << " " << hPionsFit->Integral() << endl;

        //double Ecalc1 = (3.1415 * radius_con * radius_con * f6->GetParameter(2)) / 675; //из геометрии V=1/3pi*h*r^2
        // double Ecalc2 = (3.1415 * radius_con_p * radius_con_p * f_p->GetParameter(2)) / 675;                                         //из геометрии V=1/3pi*h*r^2
        //double Ecalc2 = (3.1415 * radius_con_p * radius_con_p * fabs(f_p->GetParameter(2))) / 675;											  //из геометрии V=1/3pi*h*r^2
        //double Ecalc = Ecalc1 - Ecalc2;
        //myfile << entry << " geo " << Ecalc << endl;

        //TCanvas *c1 = new TCanvas("Surface fitted pions", "Surface fitted pions");
        //hPionsFit->Draw();
        //c1->Print("h1.png(", "png");
        //c1->SaveAs(Form("%e.png", entry));

        double Edep = edep_7sect_1 + edep_7sect_2;
        if (radius_con > 0)
        {
            //myfile << Ecalc / 1000 << " " << Edep << " " << impPar << endl;
            EdepEcalc->Fill(Ecalc / 1000, Edep);
        }
        //myfile << entry << " " << Ecalc << " " << radius_con << " " << pion_energy << " " << f8->GetParameter(2) << "  " << f8->GetParameter(0) << " " << f8->GetParameter(2) << " " << hPionsFit->Integral() << endl;
        //EdepEcalc->Fill(Ecalc / 1000, Edep_2);

        ///////////////////// paste out
    }
    cout << endl;

    ///////////////////////////////////////
    ///////////////////////////////////////
    //////////////  draw  /////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////

    TCanvas *canv_surface = new TCanvas("Surface", "Surface");
    th2_hist_ptr = ((TH2 *)(gDirectory->FindObjectAny("PSD_energy_surface")));
    th2_hist_ptr->SetTitle("PSD energy surface; X [cm]; Y [cm]");
    th2_hist_ptr->GetZaxis()->SetTitle("E [MeV]");
    th2_hist_ptr->GetXaxis()->SetTitleSize(0.05);
    th2_hist_ptr->GetYaxis()->SetTitleSize(0.05);
    th2_hist_ptr->GetZaxis()->SetTitleSize(0.05);
    th2_hist_ptr->Draw("lego2z");

    TCanvas *canv_energy_sub = new TCanvas(" uniform distr ", " uniform distr ");
    gStyle->SetOptFit(1111);
    hEnergy->SetTitle("Energy per rings uniform distr");
    hEnergy->GetXaxis()->SetTitle("X [cm]");
    hEnergy->GetYaxis()->SetTitle("Y [cm]");
    hEnergy->GetZaxis()->SetTitle("E [GeV]");
    hEnergy->GetXaxis()->SetTitleSize(0.05);
    hEnergy->GetYaxis()->SetTitleSize(0.05);
    hEnergy->Draw("lego2z");

    //f2->Draw("surf same bb");
    auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
    /*TCanvas *canv_surface_diff = new TCanvas("Surface diff", "Surface diff");
    gStyle->SetOptFit(1111);

    hDiff->Draw("lego2z");
    legend->AddEntry((TObject *)0, TString::Format("Chi = %g", f6->GetChisquare()), "");
    cout << " Chi = " << f6->GetChisquare() << endl;*/

    /* TCanvas *canv_surface_fitted_fill = new TCanvas("Surface fit fill", "Surface fitted fill");
    hFin2->SetTitle("hFin2");
    hFin2->Draw("lego2z");*/

    TCanvas *canv_surface_fin = new TCanvas("fin", "fin");
    hFin->SetTitle("final");
    hFin->GetXaxis()->SetTitle("X [cm]");
    hFin->GetYaxis()->SetTitle("Y [cm]");
    hFin->GetXaxis()->SetTitleSize(0.05);
    hFin->GetYaxis()->SetTitleSize(0.05);
    hFin->Draw("lego2z");

    TCanvas *canv_surface_fitted = new TCanvas("Surface fitted", "Surface fitted");
    hPionsFit->SetTitle("hPionsFit");
    hPionsFit->GetXaxis()->SetTitle("X [cm]");
    hPionsFit->GetYaxis()->SetTitle("Y [cm]");
    hPionsFit->GetXaxis()->SetTitleSize(0.05);
    hPionsFit->GetYaxis()->SetTitleSize(0.05);
    hPionsFit->Draw("lego2z");

    //const int n_lev = 20;
    TCanvas *canv_surface_edepecalc = new TCanvas(" c_EdepEcalc", " c_EdepEcalc");
    gStyle->SetPalette(55);
    EdepEcalc->SetTitle("EdepEcalc");
    EdepEcalc->GetXaxis()->SetTitle("E_fit [GeV]");
    EdepEcalc->GetYaxis()->SetTitle("E_dep [GeV]");
    EdepEcalc->GetXaxis()->SetTitleSize(0.05);
    EdepEcalc->GetYaxis()->SetTitleSize(0.05);
    EdepEcalc->SetContour(200);
    EdepEcalc->Draw("colz");
    //TFile histoFileFull2("EdepEcalc_SMM_min.root", "RECREATE");
    //EdepEcalc->Write();

    TCanvas *canv_edep1_vs_angle = new TCanvas(" edep_v_angle", " edep_v_angle");
    g->GetYaxis()->SetTitle("E_max [GeV]");
    g->GetXaxis()->SetTitle("angle [deg]");
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetTitleSize(0.05);
    g->Draw("AP");

    TCanvas *canv_r_vs_emax = new TCanvas(" rad vs emax", " rad vs emax");
    g1->GetYaxis()->SetTitle("E_max [GeV]");
    g1->GetXaxis()->SetTitle("radius [cm]");
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetYaxis()->SetTitleSize(0.05);
    g1->Draw("AP");

    TCanvas *canv_imp_p3 = new TCanvas(" imp_p3", " imp_p3");
    g2->GetYaxis()->SetTitle("E_max [GeV]");
    g2->GetXaxis()->SetTitle("Impact parameter [fm]");
    g2->GetXaxis()->SetTitleSize(0.05);
    g2->GetYaxis()->SetTitleSize(0.05);
    g2->Draw("AP");

    TCanvas *canv_p3_edep = new TCanvas(" edep_v_p3", " edep_v_p3");
    g3->GetYaxis()->SetTitle("Edep [GeV]");
    g3->GetXaxis()->SetTitle("E_max [GeV]");
    g3->GetXaxis()->SetTitleSize(0.05);
    g3->GetYaxis()->SetTitleSize(0.05);
    g3->Draw("AP");

    TCanvas *canv_p3_ecalc = new TCanvas(" emax vs ecalc", " emax vs ecalc");
    g4->GetXaxis()->SetTitle("E_max [GeV]");
    g4->GetYaxis()->SetTitle("E_rec [GeV]");
    g4->GetXaxis()->SetTitleSize(0.05);
    g4->GetYaxis()->SetTitleSize(0.05);
    g4->Draw("AP");

    TCanvas *canv_an_vs_edep = new TCanvas(" angle vs edep ", " angle vs edep ");
    g5->GetYaxis()->SetTitle("Edep [GeV]");
    g5->GetXaxis()->SetTitle("angle [deg] ");
    g5->GetXaxis()->SetTitleSize(0.05);
    g5->GetYaxis()->SetTitleSize(0.05);
    g5->Draw("AP");

    TCanvas *canv_imp_vs_angle = new TCanvas(" impact vs angle ", " impact vs angle ");
    g6->GetYaxis()->SetTitle("angle [deg]");
    g6->GetXaxis()->SetTitle("imp [fm]");
    g6->GetXaxis()->SetTitleSize(0.05);
    g6->GetYaxis()->SetTitleSize(0.05);
    g6->Draw("AP");

    TCanvas *canv_imp_vs_angle_hist = new TCanvas(" imp_vs_angle", " imp_vs_angle");
    hImpAng->SetTitle("hImpAng");
    hImpAng->GetXaxis()->SetTitle("Imp [fm]");
    hImpAng->GetYaxis()->SetTitle("angle [deg]");
    hImpAng->GetXaxis()->SetTitleSize(0.05);
    hImpAng->GetYaxis()->SetTitleSize(0.05);
    hImpAng->SetContour(5);
    hImpAng->Draw("colz");
    //TFile histoFileFull("imp_angle_QGSM_2000.root", "RECREATE");
    //hImpAng->Write();

    TCanvas *canv_imp_vs_e_p = new TCanvas(" imp_vs_ep", " imp_vs_ep");
    hImp_E_p->SetTitle("hImp_E_p");
    hImp_E_p->GetYaxis()->SetTitle("E_rec [GeV]");
    hImp_E_p->GetXaxis()->SetTitle("E_pions [GeV]");
    hImp_E_p->GetXaxis()->SetTitleSize(0.05);
    hImp_E_p->GetYaxis()->SetTitleSize(0.05);
    //hImp_E_p->SetContour(25);
    hImp_E_p->Draw("colz");

    TCanvas *canv_edep_emax = new TCanvas(" edep_emax", "edep_emax");
    hEdepEmax->SetTitle("hEdepEmax");
    hEdepEmax->GetXaxis()->SetTitle("Emax [GeV]");
    hEdepEmax->GetYaxis()->SetTitle("Edep [GeV]");
    hEdepEmax->GetXaxis()->SetTitleSize(0.05);
    hEdepEmax->GetYaxis()->SetTitleSize(0.05);
    hEdepEmax->SetContour(500);
    hEdepEmax->Draw("colz");
    TFile histoFileFull2("EdepEmax_SMM_full_1_to_1.root", "RECREATE");
    hEdepEmax->Write();

    TCanvas *canv_rad_vs_emax = new TCanvas(" r_vs_emax", " r_vs_emax");
    hRadEmax->SetTitle("hRadEmax");
    hRadEmax->GetXaxis()->SetTitle("radius [cm]");
    hRadEmax->GetYaxis()->SetTitle("E_max [GeV]");
    hRadEmax->GetXaxis()->SetTitleSize(0.05);
    hRadEmax->GetYaxis()->SetTitleSize(0.05);
    hRadEmax->SetContour(5);
    hRadEmax->Draw("colz");
    //TFile histoFileFull("rad_emax_SMM.root", "RECREATE");
    //hRadEmax->Write();

    TCanvas *canv_imp_v_E_pions = new TCanvas(" imp_v_E_pions", " imp_v_E_pions");
    hPions->SetTitle("Pions Integral");
    hPions->GetXaxis()->SetTitle("Impact parameter [fm]");
    hPions->GetYaxis()->SetTitle("E_pions [GeV]");
    hPions->GetXaxis()->SetTitleSize(0.05);
    hPions->GetYaxis()->SetTitleSize(0.05);
    //Pions->SetContour(5);
    hPions->Draw("colz");

    myfile.close();
}