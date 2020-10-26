#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <cmath>
#include <iostream>
#include <math.h>
using namespace std;
using namespace TMath;

void ImpactIt()
{
    ofstream myfile;
    myfile.open("data_check.txt");
    TGraph *g1 = new TGraph(); // using the blank constructor
    TGraph *g2 = new TGraph(); // using the blank constructor
                               // 1, 2 st
    //TFile *f_input = new TFile("/mnt/d/Work/root/builddir/macros/EdepEmax_QGSM_minus.root");
    TFile *f_input = new TFile("/mnt/d/Work/root/builddir/macros/EdepEmax_QGSM_full_1_to_1.root");
    //TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/EdepEmax_QGSM_minus.root");
    //TH2F *hist = (TH2F *)f_input->Get("pEcalcEdep_nt");
    TH2F *hist = (TH2F *)f_input->Get("EdepEmax");
    //hist->RebinX(2);
    //hist->RebinY(2);
    int Bin_Y = hist->GetYaxis()->GetNbins();
    int Bin_X = hist->GetXaxis()->GetNbins();
    cout << Bin_X << "  " << Bin_Y << endl;
    //TCutG* graph_cut = new TCutG("mycut", 4);
    //TCanvas *canvas = new TCanvas("c_Et_El_a", "c_Et_El_a", 600, 450);

    TH2F *pElEt[24];
    for (int i_pelet = 0; i_pelet < 25; i_pelet++)
    {
        pElEt[i_pelet] = new TH2F("E", "E", 1000, 0, 1, 1000, 0, 1);
    }

    TH1F *ImpPar[24];
    for (int i_impact = 0; i_impact < 25; i_impact++)
    {
        ImpPar[i_impact] = new TH1F("hImpPar", "hImpPar", 500, 0., 16.);
    }
    TH1F *hImpPar = new TH1F("ImpPar", "ImpPar", 500, 0., 16.);

    int count = 0.;
    int k = 0;
    int sum = 0;
    int i_cut;

    cout << hist->GetYaxis()->GetNbins() << " " << hist->GetXaxis()->GetNbins() << endl;
    //подсчет кол-ва событий

    /*for (int i = 0; i < hist->GetXaxis()->GetNbins(); ++i) // integral or i = 1 СДЕЛАЙ С 0 !!!!!!!!!!!!!!1
    {
        for (int j = 0; j < hist->GetYaxis()->GetNbins(); ++j)
        {
            count = hist->GetBinContent(i, j);
            if (count > 0)
            {
                sum = sum + count;
            }
        }
    }*/
    sum = hist->Integral();
    int arr_cuts[24];
    int n_cuts = 20;
    for (int init_cuts_array = 0; init_cuts_array <= n_cuts; init_cuts_array++)
    {
        arr_cuts[init_cuts_array] = std::round(sum / n_cuts * init_cuts_array);
        cout << init_cuts_array << " " << arr_cuts[init_cuts_array] << endl;
    }
    cout << sum << endl;

    double xs, ys, xf, yf, koef, xs_d, ys_d, xf_d, yf_d, koef_d, yorth, xorth, ycheck, xs_rot, ys_rot, xf_rot, yf_rot, xs_rot_d, ys_rot_d, xf_rot_d, yf_rot_d;
    double cut_point_d[24][5];
    double xorth_prev;
    const double y0 = 0.17992;
    const double a = 0.891876;
    const double b = 0.127786;
    const double th = 0.754254; //radians
    const double x0 = 0.0005;

    bool BinDup[1000][1000];
    for (int initBD = 0; initBD < 1000; initBD++)
    {
        for (int initBD2 = 0; initBD2 < 1000; initBD2++)
        {
            BinDup[initBD][initBD2] = false;
        }
    }

    int points_k = 0;
    int binnumX = 0;
    int binnumY = 0;
    bool DupFlag;
    int init_cuts_array = 0;
    i_cut = 0;

    /////////////////////////////////////
    ////         //                  ////
    ////      ////                   ////
    ////       //                    ////
    ////     ///// st stage          ////
    /////////////////////////////////////

    cout << "***************1**************" << endl;
    for (int i = -200000; i < 200000; i++) // along x
    {
        // ===========================================
        // под

        //////////////////////////////////////////////
        xs_d = i * 0.00001; //считaем точки на эллипсе 1
        ys_d = y0 - (b * sqrt(a * a - xs_d * xs_d + 2 * xs_d * x0 - x0 * x0)) / a;

        xf_d = xs_d + 0.000001; //считaем точки на эллипсе 2 (смещение)
        yf_d = y0 - (b * sqrt(a * a - xf_d * xf_d + 2 * xf_d * x0 - x0 * x0)) / a;

        // под
        xs_rot_d = x0 + (xs_d - x0) * cos(th) - (ys_d - y0) * sin(th); // получили новые коорд в повернутой СК
        ys_rot_d = y0 + (xs_d - x0) * sin(th) + (ys_d - y0) * cos(th);
        xf_rot_d = x0 + (xf_d - x0) * cos(th) - (yf_d - y0) * sin(th); // получили новые коорд в повернутой СК (смещенные)
        yf_rot_d = y0 + (xf_d - x0) * sin(th) + (yf_d - y0) * cos(th);
        koef_d = (yf_rot_d - ys_rot_d) / (xf_rot_d - xs_rot_d); // коэф наклона прямой через 2 точки на эллипсе под
        g1->SetMarkerColor(1);
        //g1->SetPoint(g1->GetN(), xs_rot_d, ys_rot_d);
        //g1->SetPoint(g1->GetN(), xs_rot_d, 0.939572 * xs_rot_d + 0.17945);

        //под осью
        for (int l = 0; l < 1000; l++) //идем по перпендикуляру
        {
            xorth = l * 0.001;
            yorth = ys_rot_d - (1 / koef_d) * (xorth - xs_rot_d); //получили ур-е перпендикуляра
            //cout << i << " | points | " << points_k << endl;
            if (points_k > arr_cuts[init_cuts_array] && ((yorth - (0.939572 * xorth + 0.17945)) < 0))
            {

                init_cuts_array = init_cuts_array + 1;
                cut_point_d[i_cut][1] = (ys_rot_d + xs_rot_d / koef_d - 0.17945) / (0.939572 + 1 / koef_d); //? +-
                cut_point_d[i_cut][2] = 0.939572 * cut_point_d[i_cut][1] + 0.17945;
                cut_point_d[i_cut][4] = 0;
                cut_point_d[i_cut][3] = koef_d * (ys_rot_d - cut_point_d[i_cut][4]) + xs_rot_d;
                if (cut_point_d[i_cut][3] < 0)
                {
                    cut_point_d[i_cut][3] = 10;
                    cut_point_d[i_cut][4] = ys_rot_d - (1 / koef_d) * (cut_point_d[i_cut][3] - xs_rot_d);
                }
                myfile << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2] << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;
                cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
                     << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << points_k << endl;
            }

            if (((yorth - (0.939572 * xorth + 0.17945)) < 0))
            {
                binnumX = hist->GetXaxis()->FindBin(xorth);
                binnumY = hist->GetYaxis()->FindBin(yorth);
                count = hist->GetBinContent(binnumX, binnumY);
                if (BinDup[binnumX][binnumY])
                {
                    DupFlag = false;
                }
                else
                {
                    BinDup[binnumX][binnumY] = true;
                    DupFlag = true;
                }

                if (DupFlag && count > 0)
                {
                    points_k = points_k + count;
                }
            }
        }
    }

    cut_point_d[i_cut][1] = 1; // x
    cut_point_d[i_cut][2] = 0.939572 * 1 + 0.17945;
    cut_point_d[i_cut][3] = 1;
    cut_point_d[i_cut][4] = 0;
    myfile << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2] << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;

    cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
         << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << points_k << endl;

    cut_point_d[i_cut][1] = 0; // x
    cut_point_d[i_cut][2] = 1;
    cut_point_d[i_cut][3] = 0;
    cut_point_d[i_cut][4] = 0.2;
    myfile << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2] << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;

    cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
         << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << points_k << endl;

    /////////////////////////////////////
    ////    ////////                 ////
    ////        ///                  ////
    ////    ///                      ////
    ////   //////// nd stage         ////
    /////////////////////////////////////

    points_k = arr_cuts[init_cuts_array - 1];
    cout << "2nd stage " << endl;
    init_cuts_array--;
    cout << "Points k " << points_k << " current array " << arr_cuts[init_cuts_array] << endl;
    bool flag_v2 = true;
    xorth_prev = 0;
    for (int i = -200000; i < 200000; i++)
    {
        // ===========================================
        //над осью
        //cout << i << "\r";
        xs = i * 0.00001; //считaем точки на эллипсе 1 над осью
        ys = (b * sqrt(a * a - xs * xs + 2 * xs * x0 - x0 * x0)) / a + y0;

        xf = xs + 0.000001; //считaем точки на эллипсе 2 (смещение) над осью
        yf = (b * sqrt(a * a - xf * xf + 2 * xf * x0 - x0 * x0)) / a + y0;

        //над
        xs_rot = x0 + (xs - x0) * cos(th) - (ys - y0) * sin(th); // получили новые коорд в повернутой СК
        ys_rot = y0 + (xs - x0) * sin(th) + (ys - y0) * cos(th);
        xf_rot = x0 + (xf - x0) * cos(th) - (yf - y0) * sin(th); // получили новые коорд в повернутой СК (смещенные)
        yf_rot = y0 + (xf - x0) * sin(th) + (yf - y0) * cos(th);
        g1->SetMarkerColor(1);
        //g1->SetPoint(g1->GetN(), xs_rot, ys_rot);

        koef = (yf_rot - ys_rot) / (xf_rot - xs_rot); // коэф наклона прямой через 2 точки на эллипсе над

        //полуось y=0,065x+0.253
        //над большой полуосью
        for (int l = 0; l < 1000; l++) //идем по перпендикуляру l = 8000
        {
            xorth = l * 0.001;
            yorth = ys_rot - (1 / koef) * (xorth - xs_rot);
            //if (l % 5 == 0 && (yorth - (0.939572 * xorth + 0.17945)) > 0)
            //g1->SetPoint(g1->GetN(), xorth, yorth);

            if (points_k > arr_cuts[init_cuts_array] && ((yorth - (0.939572 * xorth + 0.17945)) > 0))
            {
                init_cuts_array = init_cuts_array + 1;
                cut_point_d[i_cut][3] = (ys_rot + xs_rot / koef - 0.17945) / (0.939572 + 1 / koef);
                cut_point_d[i_cut][4] = 0.939572 * cut_point_d[i_cut][3] + 0.17945;
                cut_point_d[i_cut][2] = 1;
                cut_point_d[i_cut][1] = koef * (ys_rot - cut_point_d[i_cut][2]) + xs_rot;
                myfile << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2] << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;

                cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
                     << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << points_k << endl;
                xorth_prev = xorth;
            }

            binnumX = hist->GetXaxis()->FindBin(xorth);
            binnumY = hist->GetYaxis()->FindBin(yorth);
            if (((yorth - (0.939572 * xs_rot + 0.17945)) > 0) && xorth > 0 && yorth > 0)
            {
                binnumX = hist->GetXaxis()->FindBin(xorth);
                binnumY = hist->GetYaxis()->FindBin(yorth);
                count = hist->GetBinContent(binnumX, binnumY);
                if (BinDup[binnumX][binnumY])
                {
                    DupFlag = false;
                }
                else
                {
                    BinDup[binnumX][binnumY] = true;
                    DupFlag = true;
                }

                if (DupFlag && count > 0 /* && (koef * (ys_rot - 30) + xs_rot) > 0*/)
                {
                    points_k = points_k + count;
                }
            }
        }
    }

    cut_point_d[i_cut][1] = 1; // x
    cut_point_d[i_cut][2] = 1; // y
    cut_point_d[i_cut][3] = 1;
    cut_point_d[i_cut][4] = 0.939572 * 1 + 0.17945; // y=0,065x+0.253;
    myfile << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2] << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;

    cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
         << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << points_k << " | " << init_cuts_array + 1
         << " |arr_cuts| " << arr_cuts[init_cuts_array + 1] << " | Delta " << points_k - arr_cuts[init_cuts_array + 1] << endl;

    ifstream fp("data_check.txt");
    for (int row = 0; row < 24; row++)
    {
        for (int col = 1; col < 5; col++)
        {
            fp >> cut_point_d[row][col];
            if (!fp)
            {
                cout << "error" << endl;
            }
        }
    }

    for (i_cut = 0; i_cut < 24; i_cut++)
    {
        cout << " CUTS " << i_cut << " | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
             << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;
    }
    //sleep(5000);

    /////////////////////////////////////
    ////    ////////                 ////
    ////        ///                  ////
    ////      ///                    ////
    ////        ///                  ////
    ////   //////// rd stage (events)////
    /////////////////////////////////////

    TCutG *graph_cut = new TCutG("mycut", 5);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gStyle->SetStatW(0.3);
    gStyle->SetStatH(0.2);
    gStyle->SetOptStat(0000);
    gStyle->SetOptFit(00110);
    //root6 kBird palette

    Char_t buf[1000], buf1[1000], buff[1000], buf2[1000], buf3[1000];
    Double_t edep_7sect_1, edep_7sect_2, edep_7sect_mod15_31_1, edep_7sect_mod15_31_2, edep_7sect_mod7_39_1, edep_7sect_mod7_39_2;
    Double_t edep;
    Double_t edepMod_1[91], edepMod_2[91], edepMod_3[91], edepMod_4[91], edepMod_5[91], edepMod_6[91], edepMod_7[91];

    ifstream fp2("QGSM_11_Edep_Emax.txt");
    double edep_read, emax, impPar;
    if (fp2.fail()) // checks to see if file opended
    {
        cout << "error fp2" << endl;
        return 1; // no point continuing if the file didn't open...
    }
    while (!fp2.eof()) // reads file to end of *file*, not line
    {
        fp2 >> emax;      // read first column number
        fp2 >> edep_read; // read second column number
        fp2 >> impPar;    // read third column number
        bool IsPointInside = false;
        hImpPar->Fill(impPar);
        for (int iii = 0; iii < 23; iii++)
        { //проверяем, к какому сектору относится точка
            if (iii == 14)
                continue;
            graph_cut->SetPoint(0, cut_point_d[iii][1], cut_point_d[iii][2]);
            graph_cut->SetPoint(1, cut_point_d[iii][3], cut_point_d[iii][4]);
            graph_cut->SetPoint(2, cut_point_d[iii + 1][3], cut_point_d[iii + 1][4]);
            graph_cut->SetPoint(3, cut_point_d[iii + 1][1], cut_point_d[iii + 1][2]);
            graph_cut->SetPoint(4, cut_point_d[iii][1], cut_point_d[iii][2]);

            if (graph_cut->IsInside(emax, edep_read))
            {
                if (iii == 22)
                {
                    ImpPar[13]->Fill(impPar);
                    pElEt[13]->Fill(emax, edep_read);
                }
                else
                {
                    ImpPar[iii]->Fill(impPar);
                    pElEt[iii]->Fill(emax, edep_read);
                }
            }
        } //for (int iii = 0; iii < 13; iii++)
    }     //while

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    TF1 *fit1 = new TF1("fit1", "gaus");
    myfile << "Mean Sigma " << endl;

    for (int i_impact = 0; i_impact < 23; i_impact++)
    {
        if (i_impact != 14)
        {
            ImpPar[i_impact]->Fit(fit1, "Q");
            double sigma = fit1->GetParameter(2);
            double mean = fit1->GetParameter(1);
            cout << mean << " " << sigma << endl;
            myfile << mean << " " << sigma << " " << i_impact << endl;
        }
    }

    TCanvas *canv_imp = new TCanvas(" impact", "impact");
    hImpPar->Draw();
    hImpPar->GetXaxis()->SetTitle("Impact parameter [fm]");
    hImpPar->GetYaxis()->SetTitle("counts");
    hImpPar->GetXaxis()->SetTitleSize(0.05);
    hImpPar->GetYaxis()->SetTitleSize(0.05);
    for (int i_imp = 0; i_imp < 23; i_imp++)
    {
        if (i_imp == 14)
            continue;
        if (i_imp % 10 == 0)
        {
            ImpPar[i_imp]->SetLineColor(i_imp + 15);
            ImpPar[i_imp]->Draw("same");
        }
        else
        {
            ImpPar[i_imp]->SetLineColor(i_imp % 10);
            ImpPar[i_imp]->Draw("same");
        }
    }

    TCanvas *canv_edep_emax = new TCanvas(" edep_emax", "edep_emax");
    for (int i_pelet = 0; i_pelet < 23; i_pelet++)
    {
        if (i_pelet != 14)
        {
            if (i_pelet % 10 == 0)
            {
                pElEt[i_pelet]->SetMarkerColor(i_pelet + 15);
                pElEt[i_pelet]->Draw("same");
                pElEt[i_pelet]->GetXaxis()->SetTitle("E_max [a.u.]");
                pElEt[i_pelet]->GetYaxis()->SetTitle("E_dep [a.u.]");
                pElEt[i_pelet]->GetXaxis()->SetTitleSize(0.05);
                pElEt[i_pelet]->GetYaxis()->SetTitleSize(0.05);
            }
            else
            {
                pElEt[i_pelet]->SetMarkerColor(i_pelet % 10);
                pElEt[i_pelet]->Draw("same");
            }
        }
    }

    myfile.close();

    //g1->Draw("P");
    //g2->Draw("P");
    /*TMultiGraph *mg = new TMultiGraph();
	mg->Add(g1,"P");
	mg->Add(g2,"P");
	mg->Draw("P");
*/
    //canvas->SaveAs("test.png");
}