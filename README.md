# SCLC-Commentary

*Accompanying paper:*
Evaluating the prognostic value of mutational signatures in small-cell lung cancer through data-driven threshold optimization and signature assignment

*Order of running code:*
Data-Import-and-Cleaning.R -> Functions.R -> Main-Analysis.R -> Plots.R -> Figures.R

*Needed files:*
Download deconstructSigs_1.8.0.tar.gz from https://github.com/raerose01/deconstructSigs/releases in order to assign signatures with deconstructSigs.

Run bootstrap_analysis_SBS4_survival_optimal_cutoff_paper, bootstrap_analysis_SBS13_survival_optimal_cutoff_paper, bootstrap_analysis_SBS4_TMB_optimal_cutoff_paper, bootstrap_analysis_SBS13_TMB_optimal_cutoff_paper with num_bs = 1000 in order to replicate the bootstrap optimal cutoffs. 