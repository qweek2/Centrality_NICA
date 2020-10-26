# Centrality_NICA

## Files description.

twod_gauss_energy_uni_distr.cpp -- file for creating different histograms with dependencies (e.g. impact vs. angle, edep vs emax etc.).

FitIt.cpp is for fitting histograms in order to divide it into sectors (10%, 5% etc.)

RunFits.sh is needed if the defaults do not result in a fitting with the 1st time.

EllipseTGraphRMM.cpp. FitIt uses this file to fit with the ellipse.

EdepEmax_ell_QGSM_high_bin_2_percent.cpp is a macro for splitting the histograms into central classes.




1. FitIt is for Fitting. Set path to the hist you want to fit with ellipse.

Example:
        
        TFile* f_input = new TFile("path/file_name.root");
        TH2F* hist = (TH2F*)f_input->Get("hist_name");
To use: 

        .L FitData.cpp
        FitIt()

in case the program does not work correctly cut off the noise using the i and j parameters in the loop for (line 61, 63).



	
2. Impact.cpp is for colorizing hists, dividing it by events %. Impact factors visualising.

Set path to file and histo:

	TFile *f_input = new TFile("path/file_name.root");
	TH2F *hist = (TH2F *)f_input->Get("hist_name");
To use:

        .L impact.cpp
        ImpactIt()
 
