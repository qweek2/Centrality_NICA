# Centrality_NICA

This files implemented to determine centrality at MPD/NICA experiment with FHCal.

**The preliminary version of the code, requires more automatization, however, is currently working.**

## Files description

**twod_gauss_energy_uni_distr.cpp**. Creating different histograms with dependencies (e.g. impact vs. angle, edep vs emax etc.).

**FitIt.cpp** is for fitting histograms in order to divide it into sectors (10%, 5% etc.)

**RunFits.sh** is needed if the defaults do not result in a fitting with the 1st time.

**Ellipse.cpp**. FitIt uses this file to fit with an ellipse.

**Impact.cpp** is a macro for splitting the histograms into centrality classes.

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

**4. Impact.cpp**

Final stage. Not very user-friendly and requires manual actions. Once you have finished fitting and selected the right fit, you have the parameters of a curve that envelope the histogram data. The output looks like this:

        EXT PARAMETER                                   STEP         FIRST   
        NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
        1  x0           5.00000e-04   2.51883e-04   3.74524e-05** at limit **
        2  y0           2.37348e-01   3.36057e-03  -4.11668e-06  -6.81629e-02
        3  a            7.54648e-01   6.39113e-03  -1.85627e-06   1.10214e-02
        4  b            1.75969e-01   2.27828e-03  -4.68362e-06   5.60232e-02
        5  theta        3.80269e+01   2.85929e-01   5.07917e-06  -1.90082e-02
        top x 0.594961
        top y 0.702224
        theta deg 0.663675

Example of the fit picture:

<img src="https://raw.githubusercontent.com/qweek2/Centrality_NICA/master/pics_for_readme/7_run%20Ellipse%20Fit%20.png" width="400">

This way you have all the data to go on dividing the histogram into sectors (in this code 5% sectors). However, there are still a couple of steps left, it is necessary to make changes (enter these parameters) in the program code (at this stage it is organized manually, in the future everything will be automated). In case your fit is obtained at zero iteration of FitIt.cpp you should do the following:

Set paths and name of your .txt with data from step 1:

    TFile *f_input = new TFile("/mnt/d/Work/root/builddir/macros/EdepEmax_QGSM_full_1_to_1.root");
    TH2F *hist = (TH2F *)f_input->Get("EdepEmax");
    
    ifstream fp2("QGSM_11_Edep_Emax.txt");
    
Set parameters of the curve in lines 64 - 68. Example:
    
    const double y0 = 0.17992;
    const double x0 = 0.0005;
    const double a = 0.891876;
    const double b = 0.127786;
    const double th = 0.754254;
    
Then you have to determine k and b parameters from the y=kx+b equation. You should solve the system of 2 equations putting x0, y0 and top_x, top_y. Then replace all the k and b in program to your own (ctrl+H).

Set the range in for cycles:
   
    for (int i = -200000; i < 200000; i++) // along x
    
and further.

The last thing to do before starting is to determine at what point the transition from the edge of the lower branch to the edge of the upper branch occurs (see the picture illustrating this transition).
<img src="https://raw.githubusercontent.com/qweek2/Centrality_NICA/master/pics_for_readme/cuts.png" width="400">



