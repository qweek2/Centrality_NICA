#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <cmath>
#include <iostream>
#include "EllipseTGraphRMMSMM.cpp"
using namespace std;

//run RunsFit.sh

void FitIt(int run)
{
	
		
	//TFile* f_input = new TFile("hists_BMN_ZDC_Wall_AuAu4_5mb_pEdepSumZ_pEdepSquareZ.root");
	//TFile* f_input = new TFile("hists_NICA_5000bins.root");
	//TFile* f_input = new TFile("hists_NICA_EtEl_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root");
	//TFile* f_input = new TFile("hists_BMN_AuAu4_5mb.root");
	//TFile* f_input = new TFile("hists_NICA_11mb_edepecalc_1000bin_SMM_new.root");
	TFile* f_input = new TFile("EdepEmax_SMM_full_1_to_1.root");


	//TH2F* hist = (TH2F*)f_input->Get("pEdepLnSumZ_nt_Z_GT1");
	//TH2F* hist = (TH2F*)f_input->Get("pElEt_1");
	//hist->RebinX(2);
	//hist->RebinY(2);
	TH2F* hist = (TH2F*)f_input->Get("EdepEmax");
	//TH2F* hist = (TH2F*)f_input->Get("EdepEcalc");
	//TH2F* hist = (TH2F*)f_input->Get("pEdepLnSquareZ_nt_Z_GT1");

	TCanvas *canvas = new TCanvas();
	TGraph* g = new TGraph(); // using the blank constructor
	
	//auto nPoints = graph.GetN(); // number of points in your TGraph

	double count =0.;
	int Bin_Y = hist->GetYaxis()->GetNbins();
	int Bin_X = hist->GetXaxis()->GetNbins();	
	   TMarker *m;


/*
	for(int k=0; k < hist->GetYaxis()->GetNbins(); ++k) 
	{
		for (int p=0; p < hist->GetXaxis()->GetNbins(); ++p)
		{
				count = hist->GetBinContent(k,p);
				if (count > 0) 
				{
					new_hist->Fill(k,p);
				}
		}
	}

	new_hist->Draw();
	c->Print("hist_new.png");
*/
	cout << hist->GetYaxis()->GetNbins() << " "<< hist->GetXaxis()->GetNbins() << endl;
	//sleep(2000);
	//ij 120 for edepecalc 5000 bins
	//ij = 40 for ecalcedep 1200 bins // 40 bins = 1 step on x
	for(int i=run; i < hist->GetXaxis()->GetNbins(); ++i) 
	{
		for (int j=run; j < hist->GetYaxis()->GetNbins(); ++j)
		{
			count = hist->GetBinContent(i,j);
			//cout << i<< " " << j << " "<< count << endl;
			if (count > 0)
			{
				double x,y;
				y = ((TAxis*)hist->GetYaxis())->GetBinCenter(j);
				x = ((TAxis*)hist->GetXaxis())->GetBinCenter(i);
				g->SetMarkerColor(1);
				//g->SetMarkerStyle(1);
				g->SetPoint(g->GetN(),x,y);
				cout << count << " " << i << "/" << Bin_X << " " << j << "/" << Bin_Y << " " << "\r";
			}
		}
	}

	g->Draw("AP");
	EllipseTGraphRMMSMM(run,g);
	
}