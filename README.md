# Centrality_NICA
NICA

1. FitData is for Fitting. Set path to the hist you want to fit with ellipse.

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
 
