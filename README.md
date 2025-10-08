# SCLC-Commentary

*Accompanying paper:*
Evaluating the prognostic value of mutational signatures in small-cell lung cancer through data-driven threshold optimization and signature assignment

*Order of running code:*
[Data-Import-and-Cleaning.R](https://github.com/rishug123/SCLC-Commentary/blob/main/Data-Import-and-Cleaning.R) -> [Functions.R](https://github.com/rishug123/SCLC-Commentary/blob/main/Functions.R) -> [Main-Analysis.R](https://github.com/rishug123/SCLC-Commentary/blob/main/Main-Analysis.R) -> [Plots.R](https://github.com/rishug123/SCLC-Commentary/blob/main/plots/Plots.R) -> [Figures.R](https://github.com/rishug123/SCLC-Commentary/blob/main/plots/Figures.R)

*Needed files:*
Download deconstructSigs_1.8.0.tar.gz from https://github.com/raerose01/deconstructSigs/releases in order to assign signatures with deconstructSigs.

Run Bootstrap_calculations.R with num_bs = 1000 in order to replicate the bootstrap optimal cutoffs. 

