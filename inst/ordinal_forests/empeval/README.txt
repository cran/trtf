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

- To generate the data, run '01_data_simulation.R'. As a result, in the same
  directory, where '01_data_simulation.R' is located, a new file called
  'simulated_data.rda' will be created.
  (Note that the whole simulation study with the original experimental conditions
  is quite time-consuming. You may want to run it on a multi-core server
  and/or parallelize it.)

- All forest-based prognostic methods with the corresponding predict functions
  as well as log-likelihood calculations are realized in 'competitors.R'.

- To run the empirical evaluation on the simulated data sets, run '02_empeval.R'.
  At the end of this, a new RDA-file called 'results_empeval_250.rda' will be
  written and hence saved in the same directory.

- The R-script '03_summary_plot.R' allows you to plot the results according to
  Figure 1 of the manuscript.

