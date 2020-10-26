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
#include "Ellipse.cpp"
using namespace std;

//run RunsFit.sh

void FitIt(int run)
{
	TFile* f_input = new TFile("EdepEmax_SMM_full_1_to_1.root");
	TH2F* hist = (TH2F*)f_input->Get("EdepEmax");
	TCanvas *canvas = new TCanvas();
	TGraph* g = new TGraph(); // using the blank constructor
	
	double count =0.;
	int Bin_Y = hist->GetYaxis()->GetNbins();
	int Bin_X = hist->GetXaxis()->GetNbins();	
	   TMarker *m;

	cout << hist->GetYaxis()->GetNbins() << " "<< hist->GetXaxis()->GetNbins() << endl;

	for(int i=run; i < hist->GetXaxis()->GetNbins(); ++i) 
	{
		for (int j=run; j < hist->GetYaxis()->GetNbins(); ++j)
		{
			count = hist->GetBinContent(i,j);
			if (count > 0)
			{
				double x,y;
				y = ((TAxis*)hist->GetYaxis())->GetBinCenter(j);
				x = ((TAxis*)hist->GetXaxis())->GetBinCenter(i);
				g->SetMarkerColor(1);
				g->SetPoint(g->GetN(),x,y);
				cout << count << " " << i << "/" << Bin_X << " " << j << "/" << Bin_Y << " " << "\r";
			}
		}
	}

	g->Draw("AP");
	Ellipse(run,g);
}