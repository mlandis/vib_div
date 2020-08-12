RevBayes scripts to run the main analysis in

> MJ Landis, DAR Eaton, WL Clement, B Park, EL Spriggs, PW Sweeney, EJ Edwards, & MJ Donoghue. Joint phylogenetic estimation of geographic movements and biome shifts during the global diversification of Viburnum. Systematic Biology (advance access), https://doi.org/10.1093/sysbio/syaa027.

To run the analysis, call `rb run_quick.Rev` from within this directory. That script will in turn call all the necessary model scripts. There are many settings you can tweak within `run_quick.Rev`. I've done my best to comment the scripts so that others may modify them freely. Feel welcome to contact me with any quesitons.

The `plot/` directory contains R scripts to generate the plots in the manuscript using the files in `../output/`.
