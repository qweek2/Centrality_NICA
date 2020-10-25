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
#include <chrono>
#include <iostream>
using namespace std;
using namespace TMath;

void ImpactIt()
{

	ofstream myfile;
	myfile.open("data_11_GeV_QGSM.txt");
	gSystem->Beep(650, 300);
	using clock = std::chrono::steady_clock;
	clock::time_point start = clock::now();
	gSystem->Beep(500, 300);

	// 1, 2 st
	//TFile *f_input = new TFile("/mnt/d/Work/root/builddir/macros/pelet_11gev_SMM_full.root");
	TFile *f_input = new TFile("/mnt/d/Work/root/builddir/macros/pelet_11gev_QGSM_full.root");
	//TH2F *hist = (TH2F *)f_input->Get("pEcalcEdep_nt");
	TH2F *hist = (TH2F *)f_input->Get("pElEt_11GeV");

	int Bin_Y = hist->GetYaxis()->GetNbins();
	int Bin_X = hist->GetXaxis()->GetNbins();
	cout << Bin_X << "  " << Bin_Y << endl;
	//TCutG* graph_cut = new TCutG("mycut", 4);
	TCanvas *canvas = new TCanvas("c_Et_El_a", "c_Et_El_a", 600, 450);
	TCanvas *canvas2 = new TCanvas("2", "2", 600, 450);

	// 3 st
	//double fNevents = 0, gdTEvents, gdTEvents_1 = 0;
	//double edep = 0, impPar_scint = 0;
	TH2F *pElEt2 = new TH2F("Pelet_5_GeV", "Pelet_5_GeV", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_1 = new TH2F("Pelet_5_GeV", "Pelet_5_GeV", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_2 = new TH2F("pElEt_2", "pElEt_2", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_3 = new TH2F("pElEt_3", "pElEt_3", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_4 = new TH2F("pElEt_4", "pElEt_4", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_5 = new TH2F("pElEt_5", "pElEt_5", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_6 = new TH2F("pElEt_6", "pElEt_6", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_7 = new TH2F("pElEt_7", "pElEt_7", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_8 = new TH2F("pElEt_8", "pElEt_8", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_9 = new TH2F("pElEt_9", "pElEt_9", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_10 = new TH2F("pElEt_10", "pElEt_10", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_11 = new TH2F("pElEt_11", "pElEt_11", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_12 = new TH2F("pElEt_12", "pElEt_12", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_13 = new TH2F("pElEt_13", "pElEt_13", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_14 = new TH2F("pElEt_14", "pElEt_14", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_15 = new TH2F("pElEt_15", "pElEt_15", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_16 = new TH2F("pElEt_16", "pElEt_16", 2000, 0, 25, 2000, 0, 3);
	TH2F *pElEt_17 = new TH2F("pElEt_17", "pElEt_17", 2000, 0, 25, 2000, 0, 1);
	TH2F *pElEt_18 = new TH2F("pElEt_18", "pElEt_18", 2000, 0, 12, 2000, 0, 1);
	TH2F *pElEt_19 = new TH2F("pElEt_19", "pElEt_19", 2000, 0, 12, 2000, 0, 1);
	TH2F *pElEt_20 = new TH2F("pElEt_20", "pElEt_20", 2000, 0, 12, 2000, 0, 1);
	TH2F *pElEt_21 = new TH2F("pElEt_21", "pElEt_21", 2000, 0, 12, 2000, 0, 1);
	TH2F *pElEt_22 = new TH2F("pElEt_22", "pElEt_22", 2000, 0, 12, 2000, 0, 1);
	TH2F *pElEt_23 = new TH2F("pElEt_23", "pElEt_23", 2000, 0, 12, 2000, 0, 1);

	//TFile *analzdc = new TFile();
	//analzdc = TFile::Open("analzdc_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root");
	//TTree *theTree = (TTree*)analzdc->Get("nt1");
	//gdTEvents=theTree->GetEntries();
	//theTree = SetBranchAddress("impPar",&impPar);

	int count = 0.;
	int k = 0;
	int sum = 0;
	int i_cut = 0;
	double fNevents = 0, gdTEvents, gdTEvents_1 = 0;

	TMarker *m;

	cout << hist->GetYaxis()->GetNbins() << " " << hist->GetXaxis()->GetNbins() << endl;
	int s = 0;

	//подсчет кол-ва событий

	for (int i = 0; i < hist->GetXaxis()->GetNbins(); i++) // integral or i = 1 СДЕЛАЙ С 0 !!!!!!!!!!!!!!1
	{
		for (int j = 0; j < hist->GetYaxis()->GetNbins(); j++)
		{
			count = hist->GetBinContent(i, j);
			//cout << i<< " " << j << " "<< count << endl;
			if (count > 0)
			{
				sum = sum + count;
				//cout << sum << "\r";
			}
		}
	}

	int arr_cuts[10];
	int n_cuts = 10;
	for (int init_cuts_array = 0; init_cuts_array < n_cuts + 1; init_cuts_array++)
	{
		arr_cuts[init_cuts_array] = round(sum / n_cuts * init_cuts_array);
		cout << init_cuts_array << " " << arr_cuts[init_cuts_array] << endl;
	}
	cout << sum << endl;
	//sleep(1000);
	//уравнение полуоси y=0.606x+0.222

	double xs, ys, xf, yf, koef, xs_d, ys_d, xf_d, yf_d, koef_d, yorth, xorth, ycheck, xs_rot, ys_rot, xf_rot, yf_rot, xs_rot_d, ys_rot_d, xf_rot_d, yf_rot_d;
	double cut_point_d[12][5];
	double delta_x = 0;
	double xorth_prev;
	double a = 3.85052;
	double b = 0.0661745;
	double th = 0.06981; //наклон эллипса 0,769760013422076672 // 44.1* градуса
	double x0 = 4.77954;
	double y0 = 0.402009;
	/*
	double a = 23.9928;
	double b = 2.64050;
	double th = 0.7698; //наклон эллипса 0,769760013422076672 // 44.1* градуса
	double x0 = 1.195;
	double y0 = 5.0276;
	*/
	//b≈149.46 and k≈0.53610 qgsm 11 gev
	bool BinDup[2000][2000];
	//double means[11];
	for (int initBD = 0; initBD < 2000; initBD++)
	{
		for (int initBD2 = 0; initBD2 < 2000; initBD2++)
		{
			BinDup[initBD][initBD2] = false;
		}
	}
	int points_k = 0;
	int points_k2 = 0;

	int binnumX = 0;
	int binnumY = 0;
	bool DupFlag;
	bool flag_v = true;
	int init_cuts_array = 0;

	binnumX = hist->GetXaxis()->FindBin(x0);
	binnumY = hist->GetYaxis()->FindBin(y0);
	cout << binnumX << " " << binnumY << endl;
	binnumX = hist->GetXaxis()->FindBin(8.61794);
	binnumY = hist->GetYaxis()->FindBin(0.707497);
	cout << binnumX << " " << binnumY << endl;

	/////////////////////////////////////
	////                             ////
	////         //                  ////
	////      ////                   ////
	////       //                    ////
	////     ///// st stage          ////
	/////////////////////////////////////
	//b = 0.222 and k≈0.064555

	for (int i = 0; i < hist->GetXaxis()->GetNbins(); i++)
	{
		for (int j = 0; j < hist->GetYaxis()->GetNbins(); j++)
		{
			if (j < 0.53610 * i + 149)
			{
				count = hist->GetBinContent(i, j);
				k = k + count;
				s = k;
				//cout << i<< " " << j << " "<< count << endl;
				if (k > arr_cuts[init_cuts_array])
				{
					init_cuts_array++;
					double x, y;
					y = ((TAxis *)hist->GetYaxis())->GetBinCenter(j);
					x = ((TAxis *)hist->GetXaxis())->GetBinCenter(i);

					cut_point_d[i_cut][1] = x;
					cut_point_d[i_cut][2] = 0.064555 * cut_point_d[i_cut][1] + 0.222;
					cut_point_d[i_cut][4] = 0;
					cut_point_d[i_cut][3] = x;
					myfile << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2] << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;

					cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
						 << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << k << endl;
					i_cut++;
				}
			}
		}
	}

	cut_point_d[i_cut][1] = 25;
	cut_point_d[i_cut][2] = 0.064555 * cut_point_d[i_cut][1] + 0.222;
	cut_point_d[i_cut][4] = 0;
	cut_point_d[i_cut][3] = 25;
	myfile << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2] << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;
	cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
		 << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << k << endl;
	i_cut++;

	/////////////////////////////////////
	////    ////////                 ////
	////        ///                  ////
	////    ///                      ////
	////   //////// nd stage         ////
	/////////////////////////////////////
	cout << " %%%%%%%%% 2 %%%%%%%%%" << endl;
	k = arr_cuts[init_cuts_array];

	TH1F *ImpPar[30];
	for (int i_impact = 0; i_impact < 30; i_impact++)
	{
		ImpPar[i_impact] = new TH1F("hImpPar", "hImpPar", 400, 0., 17.);
	}
	TH1F *hImpPar = new TH1F("hImpPar", "hImpPar", 400, 0., 17.);

	for (int i = 0; i < hist->GetXaxis()->GetNbins(); i++) //
	{
		for (int j = 0; j < hist->GetYaxis()->GetNbins(); j++)
		{
			if (j > 0.53610 * i + 149)
			{
				count = hist->GetBinContent(i, j);
				k = k + count;
				s = k;
				//cout << i<< " " << j << " "<< count << endl;
				if (k > arr_cuts[init_cuts_array] && init_cuts_array < 11)
				{
					init_cuts_array++;
					double x, y;
					y = ((TAxis *)hist->GetYaxis())->GetBinCenter(j);
					x = ((TAxis *)hist->GetXaxis())->GetBinCenter(i);

					cut_point_d[i_cut][3] = x;
					cut_point_d[i_cut][4] = 0.064555 * cut_point_d[i_cut][3] + 0.222;
					cut_point_d[i_cut][2] = 10;
					cut_point_d[i_cut][1] = x;
					myfile << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2] << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;

					cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
						 << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << k << endl;
					i_cut++;
				}
			}
		}
	}

	cut_point_d[i_cut][1] = 25; // x
	cut_point_d[i_cut][2] = 3;	// y
	cut_point_d[i_cut][3] = 25;
	cut_point_d[i_cut][4] = 0.064555 * 25 + 0.222;
	myfile << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2] << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;

	cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
		 << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << k << endl;

	cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	for (i_cut = 0; i_cut < 13; i_cut++)
	{
		cout << " CUTS " << i_cut << " | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
			 << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;
	}

	/////////////////////////////////////
	////    ////////                 ////
	////        ///                  ////
	////      ///                    ////
	////        ///                  ////
	////   //////// rd stage (events)////
	/////////////////////////////////////

	//TFile *_file0 = TFile::Open("analzdc_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root");
	TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_full_geom_QGSM_AuAu_11_mb_500_99500ev.root");
	//TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root")
	//TFile *_file0 = TFile::Open("data.root");
	TTree *nt1 = (TTree *)_file0->Get("nt1");
	gdTEvents = nt1->GetEntries();
	cout << "loop 1 -> gdTEvents " << gdTEvents << endl;
	TFile histoFile("hists_new.root", "RECREATE");
	TCutG *graph_cut = new TCutG("mycut", 5);
	TH2F *pElEt_5Gev = new TH2F("pElEt_5GeV", "pElEt_5GeV", 2000, 0, 25, 2000, 0, 3);

	const int nbCBMmods = 90; //NICA 45+45 mods (VETO1 - 1-45, VETO 2 -46-89)
	const int nbCBMsect = 7;
	Float_t zzzzNICA[nbCBMsect];
	Int_t modNb = -1, rNICA_int, rsumNICA_int;
	double xxxxNICA[nbCBMmods], yyyyNICA[nbCBMmods];
	Double_t fphi_mod[nbCBMmods] = {0.};
	Double_t fphi_sect[nbCBMmods][nbCBMsect] = {0.};
	Double_t radTodeg = 180. / TMath::Pi();
	Double_t degToRad = TMath::Pi() / 180.;
	Double_t edepMod_1[91], edepMod_2[91], edepMod_3[91], edepMod_4[91], edepMod_5[91], edepMod_6[91], edepMod_7[91];
	Double_t impPar, edep_7sect_1, edep_7sect_2, edep_7sect_mod15_31_1, edep_7sect_mod15_31_2, edep_7sect_mod7_39_1, edep_7sect_mod7_39_2;
	Double_t Et, El, Et_test, E_test, edepSect[nbCBMmods][nbCBMsect];

	//NICA 45 mods, 15x15 hole 10cm in mod 23 VETO 1
	double xCur = 60;
	double yCur = 45;
	for (Int_t iii = 0; iii < 5; iii++)
	{
		for (Int_t jjj = 0; jjj < 7; jjj++)
		{
			modNb = 7 * iii + jjj + 5;
			xxxxNICA[modNb] = xCur - (jjj + 1) * 15.;
			yyyyNICA[modNb] = yCur - (iii + 1) * 15.;
			//cout << "xx " << iii << " " << jjj << " " << modNb + 1 << " " << xxxxNICA[modNb] << " " << yyyyNICA[modNb] << endl;
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

	//NICA sections
	zzzzNICA[0] = 319.5 + 6.32;	 //sect 1
	zzzzNICA[1] = 319.5 + 19.76; //sect 2
	zzzzNICA[2] = 319.5 + 31.00; //sect 3
	zzzzNICA[3] = 319.5 + 43.24; //sect 4
	zzzzNICA[4] = 319.5 + 55.48; //sect 5
	zzzzNICA[5] = 319.5 + 67.72; //sect 6
	zzzzNICA[6] = 319.5 + 79.96; //sect 7

	//NICA
	Double_t cs, sn, r = 0;

	for (Int_t im = 0; im < nbCBMmods; im++)
	{
		if ((im + 1) != 23 && (im + 1) != 68)
		{
			r = TMath::Sqrt(yyyyNICA[im] * yyyyNICA[im] + xxxxNICA[im] * xxxxNICA[im]);
			for (Int_t is = 0; is < nbCBMsect; is++)
			{
				fphi_sect[im][is] = TMath::ATan2(r, zzzzNICA[is]);
				//if ((im + 1) == 20)
				//cout << "fphi_sect " << im + 1 << " " << is + 1 << " " << fphi_sect[im][is] << " " << fphi_sect[im][is]
				// << " " << TMath::Cos(fphi_sect[im][is]) << " " << TMath::Sin(fphi_sect[im][is]) << endl;
			}
		}
	}

	//TTree *nt1=(TTree*)_file0->Get("nt1");

	nt1->SetBranchAddress("edepMod_1", &edepMod_1);
	nt1->SetBranchAddress("edepMod_2", &edepMod_2);
	nt1->SetBranchAddress("edepMod_3", &edepMod_3);
	nt1->SetBranchAddress("edepMod_4", &edepMod_4);
	nt1->SetBranchAddress("edepMod_5", &edepMod_5);
	nt1->SetBranchAddress("edepMod_6", &edepMod_6);
	nt1->SetBranchAddress("edepMod_7", &edepMod_7);
	nt1->SetBranchAddress("impPar", &impPar);

	for (Int_t iEventN = 0; iEventN < 54440; iEventN++)
	{ // loop on events
		//cout << "iEventN = " << iEventN << endl;
		nt1->GetEntry(iEventN);
		hImpPar->Fill(impPar);

		Et = 0;
		El = 0;
		bool IsPointInside = false;

		for (Int_t im = 0; im < nbCBMmods; im++)
		{
			//cout << "edepMod_7 " << im + 1 << " " << edepMod_7[im + 1] << endl;
			for (Int_t is = 0; is < nbCBMsect; is++)
			{
				if (is == 0)
					edepSect[im][is] = edepMod_1[im + 1];
				if (is == 1)
					edepSect[im][is] = edepMod_2[im + 1] - edepMod_1[im + 1];
				if (is == 2)
					edepSect[im][is] = edepMod_3[im + 1] - edepMod_2[im + 1];
				if (is == 3)
					edepSect[im][is] = edepMod_4[im + 1] - edepMod_3[im + 1];
				if (is == 4)
					edepSect[im][is] = edepMod_5[im + 1] - edepMod_4[im + 1];
				if (is == 5)
					edepSect[im][is] = edepMod_6[im + 1] - edepMod_5[im + 1];
				if (is == 6)
					edepSect[im][is] = edepMod_7[im + 1] - edepMod_6[im + 1];
				Et += edepSect[im][is] * TMath::Sin(fphi_sect[im][is]);
				El += edepSect[im][is] * TMath::Cos(fphi_sect[im][is]);

			} //sect
		}	  //mod
		for (int iii = 0; iii < 13; iii++)
		{ //проверяем, к какому сектору относится точка
			graph_cut->SetPoint(0, cut_point_d[iii][1], cut_point_d[iii][2]);
			graph_cut->SetPoint(1, cut_point_d[iii][3], cut_point_d[iii][4]);
			graph_cut->SetPoint(2, cut_point_d[iii + 1][3], cut_point_d[iii + 1][4]);
			graph_cut->SetPoint(3, cut_point_d[iii + 1][1], cut_point_d[iii + 1][2]);
			graph_cut->SetPoint(4, cut_point_d[iii][1], cut_point_d[iii][2]);

			if (graph_cut->IsInside(El, Et) && iii == 0)
			{
				ImpPar[1]->Fill(impPar);
				pElEt_1->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 1)
			{
				ImpPar[2]->Fill(impPar);
				pElEt_2->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 2)
			{
				ImpPar[3]->Fill(impPar);
				pElEt_3->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 3)
			{
				ImpPar[4]->Fill(impPar);
				pElEt_4->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 4)
			{
				ImpPar[5]->Fill(impPar);
				pElEt_5->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 5)
			{
				ImpPar[6]->Fill(impPar);
				pElEt_6->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 6)
			{
				ImpPar[7]->Fill(impPar);
				pElEt_7->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 7)
			{
				ImpPar[8]->Fill(impPar);
				pElEt_8->Fill(El, Et);
			}
			/*else if (graph_cut->IsInside(El, Et) && iii == 8)
			{
				ImpPar[9]->Fill(impPar);
				pElEt_9->Fill(El, Et);
			}*/
			else if (graph_cut->IsInside(El, Et) && iii == 9)
			{
				ImpPar[10]->Fill(impPar);
				pElEt_10->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 10)
			{
				ImpPar[11]->Fill(impPar);
				pElEt_11->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 11)
			{
				ImpPar[8]->Fill(impPar);
				pElEt_8->Fill(El, Et);
			}
			/*else if (graph_cut->IsInside(El, Et) && iii == 12)
			{
				hImpPar13->Fill(impPar);
				pElEt_13->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 13)
			{
				hImpPar14->Fill(impPar);
				pElEt_14->Fill(El, Et);
			}
			else if (graph_cut->IsInside(El, Et) && iii == 14)
			{
				hImpPar15->Fill(impPar);
				pElEt_15->Fill(El, Et);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(El, Et) && iii == 15)
			{
				hImpPar16->Fill(impPar);
				pElEt_16->Fill(El, Et);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(El, Et) && iii == 16)
			{
				hImpPar17->Fill(impPar);
				pElEt_17->Fill(El, Et);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(El, Et) && iii == 17)
			{
				hImpPar18->Fill(impPar);
				pElEt_18->Fill(El, Et);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(El, Et) && iii == 18)
			{
				hImpPar19->Fill(impPar);
				pElEt_19->Fill(El, Et);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(El, Et) && iii == 19)
			{
				hImpPar20->Fill(impPar);
				pElEt_20->Fill(El, Et);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(El, Et) && iii == 20)
			{
				hImpPar21->Fill(impPar);
				pElEt_21->Fill(El, Et);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(El, Et) && iii == 21)
			{
				hImpPar22->Fill(impPar);
				pElEt_22->Fill(El, Et);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(El, Et) && iii == 22)
			{
				hImpPar23->Fill(impPar);
				pElEt_23->Fill(El, Et);
				cout << "++11+++++++++++++++" << endl;
			}*/
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	//cout << Et << " " << El << endl;
	//pElEt2->Fill(El, Et);
	//hImpPar->Write();

	//hImpPar1->Write();

	TF1 *fit1 = new TF1("fit1", "gaus");
	//myfile << "|  Mean  |  Sigma  |" << endl;
	for (int i_impact = 1; i_impact < 12; i_impact++)
	{
		if (i_impact != 9)
		{
			ImpPar[i_impact]->Fit(fit1);
			double sigma = fit1->GetParameter(2);
			double mean = fit1->GetParameter(1);
			//myfile << mean << " " << sigma << endl;
		}
	}
	/*
	hImpPar1->Fit("gaus");
	hImpPar2->Fit("gaus");
	hImpPar3->Fit("gaus");
	hImpPar4->Fit("gaus");
	hImpPar5->Fit("gaus");
	hImpPar6->Fit("gaus");
	hImpPar7->Fit("gaus"); //,"","",5,11);
	//hImpPar8->Fit("gaus"); //,"","",4,10);
	hImpPar9->Fit("gaus");
	hImpPar10->Fit("gaus"); //,"","",0,5);
	hImpPar11->Fit("gaus");*/
	//hImpPar12->Fit("gaus");
	//hImpPar13->Fit("gaus");
	/*hImpPar14->Fit("gaus");
	hImpPar15->Fit("gaus");
	hImpPar16->Fit("gaus");
	hImpPar17->Fit("gaus");
	hImpPar18->Fit("gaus"); 
	hImpPar19->Fit("gaus");
	hImpPar20->Fit("gaus");
	hImpPar21->Fit("gaus");
	hImpPar22->Fit("gaus");*/
	//hImpPar23->Fit("gaus");
	TCanvas *canv_imp = new TCanvas(" impPar", "impPar");
	hImpPar->GetXaxis()->SetTitle("Impact parameter [fm]");
	hImpPar->GetYaxis()->SetTitle("counts");
	hImpPar->GetXaxis()->SetTitleSize(0.05);
	hImpPar->GetYaxis()->SetTitleSize(0.05);
	hImpPar->SetLineColor(1);
	hImpPar->Draw();
	ImpPar[1]->SetLineColor(1);
	ImpPar[1]->Draw("same");
	ImpPar[2]->SetLineColor(2);
	ImpPar[2]->Draw("same");
	ImpPar[3]->SetLineColor(3);
	ImpPar[3]->Draw("same");
	ImpPar[4]->SetLineColor(4);
	ImpPar[4]->Draw("same");
	ImpPar[5]->SetLineColor(5);
	ImpPar[5]->Draw("same");
	ImpPar[6]->SetLineColor(6);
	ImpPar[6]->Draw("same");
	ImpPar[7]->SetLineColor(14);
	ImpPar[7]->Draw("same");
	ImpPar[8]->SetLineColor(8);
	ImpPar[8]->Draw("same");
	ImpPar[9]->SetLineColor(9);
	ImpPar[9]->Draw("same");
	ImpPar[10]->SetLineColor(2);
	ImpPar[10]->Draw("same");
	ImpPar[11]->SetLineColor(2);
	ImpPar[11]->Draw("same");
	/*hImpPar12->SetLineColor(8);
	hImpPar12->Draw("same");
	hImpPar13->SetLineColor(5);
	hImpPar13->Draw("same");
	hImpPar14->SetLineColor(6);
	hImpPar14->Draw("same");
	hImpPar15->SetLineColor(7);
	hImpPar15->Draw("same");
	hImpPar16->SetLineColor(8);
	hImpPar16->Draw("same");
*/
	//hImpPar4->Draw("same");
	//hImpPar5->Draw("same");
	//pElEt2->Draw();

	TCanvas *canv_pelet = new TCanvas(" Et El", " Et El");
	pElEt_1->GetXaxis()->SetTitle("El [GeV]");
	pElEt_1->GetYaxis()->SetTitle("Et [GeV]");
	pElEt_1->GetXaxis()->SetTitleSize(0.05);
	pElEt_1->GetYaxis()->SetTitleSize(0.05);
	pElEt_1->SetMarkerColor(1);
	pElEt_1->Draw();
	//TFile histoFileFull("hists_NICA.root","RECREATE");
	//pElEt_1->Write();
	pElEt_2->SetMarkerColor(2);
	pElEt_2->Draw("same");
	pElEt_3->SetMarkerColor(3);
	pElEt_3->Draw("same");
	pElEt_4->SetMarkerColor(4);
	pElEt_4->Draw("same");
	pElEt_5->SetMarkerColor(5);
	pElEt_5->Draw("same");
	pElEt_6->SetMarkerColor(6);
	pElEt_6->Draw("same");
	pElEt_7->SetMarkerColor(14);
	pElEt_7->Draw("same");
	pElEt_8->SetMarkerColor(8);
	pElEt_8->Draw("same");
	pElEt_9->SetMarkerColor(9);
	pElEt_9->Draw("same");
	pElEt_10->SetMarkerColor(1);
	pElEt_10->Draw("same");
	pElEt_11->SetMarkerColor(2);
	pElEt_11->Draw("same");
	pElEt_12->SetMarkerColor(3);
	pElEt_12->Draw("same");
	/*pElEt_13->SetMarkerColor(4);
	pElEt_13->Draw("same");
	//pElEt_14->SetMarkerColor(5);
	//pElEt_14->Draw("same");
	pElEt_15->SetMarkerColor(6);
	pElEt_15->Draw("same");
	pElEt_16->SetMarkerColor(8);
	pElEt_16->Draw("same");
	pElEt_17->SetMarkerColor(9);
	pElEt_17->Draw("same");
	pElEt_18->SetMarkerColor(1);
	pElEt_18->Draw("same");
	pElEt_19->SetMarkerColor(2);
	pElEt_19->Draw("same");
	pElEt_20->SetMarkerColor(3);
	pElEt_20->Draw("same");
	pElEt_21->SetMarkerColor(4);
	pElEt_21->Draw("same");
	pElEt_22->SetMarkerColor(5);
	pElEt_22->Draw("same");
	*/

	//canvas->SetCansvasSize(600,600);
	//canvas->SetWindowSize(500, 500);

	//g1->Draw("AP");
	//g2->Draw("AP");
	myfile.close();

	//canvas->SaveAs("test.png");
	clock::time_point end = clock::now();
	clock::duration execution_time = end - start;
	gSystem->Beep(432, 1000);
	cout << "Complete in " << (execution_time.count() / (1000000000)) << " s" << endl;
}