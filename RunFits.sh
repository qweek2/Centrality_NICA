
for run in {0..10}; 
do  
    echo -------------------------
    echo --------- $run -------------
    echo -------------------------


    root -l -b -q 'FitIt.cpp('$run')'
done