################################################################################

# Model-based Random Forests for Ordinal Regression

# Muriel Buri and Torsten Hothorn

# SUBFOLDER: "empeval"

################################################################################

This folder contains simplified and commented R-scripts used for the
empirical evaluation of the manuscript.

Please refer to the manuscript for details on the functions, models and
simulation procedures.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Reproducing the empirical evaluation results:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- To generate the data, run '01_data_simulation.R'. The resulting simulated data
  set called 'simulated_data_nsim100.rda' will be saved in the (automatically
  generated) sub-directory called 'rda'.
  (Note that the whole simulation study with the original experimental conditions
  is quite time-consuming. You may want to run it on a multi-core server.)

- All forest-based prognostic methods with the corresponding predict functions
  as well as the log-likelihood and Kullback-Leibler divergence calculations are
  realized in 'competitors.R'.

- To run the empirical evaluation on the simulated data sets, run '02_empeval.R'.
  At the end of this, the results will again be written and saved to the
  sub-directory 'rda'. The resulting files are called:
    o 'results_empeval_nsim100_ntree250.rda', and
    o 'results_empeval_nsim100_ntree2000.rda'.

- The R-script '03_summary_plots.R' allows you to plot the results according to
  Figure 1, 4, 5, 6, 7, 8, 9, 10, 11, and 12 of the manuscript.

