# Centrality_NICA

This files implemented to determine centrality at MPD/NICA experiment with FHCal.

## Files description

**twod_gauss_energy_uni_distr.cpp** --- file for creating different histograms with dependencies (e.g. impact vs. angle, edep vs emax etc.).

**FitIt.cpp** is for fitting histograms in order to divide it into sectors (10%, 5% etc.)

**RunFits.sh** is needed if the defaults do not result in a fitting with the 1st time.

**EllipseTGraphRMM.cpp**. FitIt uses this file to fit with the ellipse.

**EdepEmax_ell_QGSM_high_bin_2_percent.cpp** is a macro for splitting the histograms into central classes.

## How to

1. **twod_gauss_energy_uni_distr.cpp**

Set path to your data file here:

        TFile *_file0 = TFile::Open("/path/datafile.root");
        
Uncomment this line if you want to try pion subtraction:

        //hFin->Add(hPionsFit, -1.);
        
You can draw a lot of histograms, some of them will be already available after the first run, you can see the list in the "draw" part of the code. For futher steps we will use Edep_Emax histo as an example.

        if (radius_con > 0) //here you can set any cut you need
        {
            hEdepEmax->Fill(f8->GetParameter(2) / 1000, edep_7sect_1 + edep_7sect_2); //f8->GetParameter(2) is a cone height
            myfile << f8->GetParameter(2) / 1000 << " " << edep_7sect_1 + edep_7sect_2 << " " << impPar << endl; //write data to the .txt
        }
        
You'll need the .txt file for the next step, its name can be changed here:

            myfile.open("file_name.txt");

2. 

 
