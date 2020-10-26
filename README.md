# Centrality_NICA

This files implemented to determine centrality at MPD/NICA experiment with FHCal.

## Files description

**twod_gauss_energy_uni_distr.cpp** --- file for creating different histograms with dependencies (e.g. impact vs. angle, edep vs emax etc.).

**FitIt.cpp** is for fitting histograms in order to divide it into sectors (10%, 5% etc.)

**RunFits.sh** is needed if the defaults do not result in a fitting with the 1st time.

**Ellipse.cpp**. FitIt uses this file to fit with the ellipse.

**EdepEmax_ell_QGSM_high_bin_2_percent.cpp** is a macro for splitting the histograms into central classes.

## How to

**1. twod_gauss_energy_uni_distr.cpp**

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
            
To run:
        
        root [0] .L twod_gauss_energy_uni_distr.cpp  
        root [1] twod_gauss_energy()

**2. Ellipse.cpp**

Now you have a histogram, you need to roughly determine the center of the ellipse and the size of the axes. You need to paste these estimates into the FitIt.cpp file to make the fit more accurate. Put this values (x0, y0, a, b) into the lines 131-134. 

Specify the path to the folder (or create it) where the pictures with the fit results will be saved:

    gSystem->cd("/mnt/d/Work/INR/centrality/pics");
    gSystem->Exec("mkdir fits_ell");
    gSystem->cd("/mnt/d/Work/INR/centrality/pics/fits_ell");

**3. RunFits.sh**

Here just set number of FitIt.cpp iterations. The difference between iterations lies in the fact that probably at a single iteration fit may not work, so each run 1 bin in the x and y axis (from the left bottom) will be cutted. You can monitor the quality of the fit "online" by watching the output of Minuit in the terminal, or just look at the pictures and select the fit after finishing.

To run:

        ./RunFits.sh

**4.  **

