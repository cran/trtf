################################################################################

# Model-based Random Forests for Ordinal Regression

# Muriel Buri and Torsten Hothorn

# SUBFOLDER: "ALS"

################################################################################

This folder contains simplified and commented R-scripts used for the re-analysis
of the Respiratory sub-item of the ALS Functional Rating Scale Revised (ALSFRS-R).

Please refer to the manuscript for details on the functions, models and
simulation procedures.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Reproducing the re-analysis:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
- As written in the manuscript, the patient-level data is available to registered
  users from https://nctu.partners.org/ProACT.
  The code for data preprocessing is available in the \pkg{TH.data} add-on
  R package.

- The R-script '01_ALS_sample.R' generates the 100 random splits of the data set
  and saves them as 'ALS_samples.rda' to the same directory.

- To compare the seven methods that predict the “Respiratory” ALSFRS-R score at
  half a year after diagnosis based on the 100 random splits are estimated while
  running the file '02_run_1_100.R'. For this, the file 'ALS_analysis.R' is ran
  within a for-loop using the functions of the 'ALS_competitors.R' file.
  The out-of-sample log-likelihoods and the quadratic weighted Cohen's Kappa are
  saved to the file 'ALS_results_combined.rda' in the same directory.

- The R-script '03_summary_plots.R' allows you to plot the results according to
  Figure 2 and 3 of the manuscript.


