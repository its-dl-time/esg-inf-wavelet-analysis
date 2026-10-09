#!/usr/bin/env Rscript
# Run from the repository root. This reproduces the archived smoothing=5/0.25 analysis.
source('robustness/legacy/00_install_packages.R')
source('R/prepare_data.R')
source('robustness/legacy/02_run_wavelet.R')
source('robustness/legacy/03_diagnostics.R')
source('robustness/legacy/06_stationarity_robustness.R')
source('robustness/legacy/07_compare_robustness.R')
source('robustness/legacy/04_figures_and_tables.R')
cat('Legacy analysis completed.\n')
