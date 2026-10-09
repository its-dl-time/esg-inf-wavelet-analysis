# Earlier full analysis and optional sensitivities

Run `Rscript robustness/legacy/run_all.R` from the repository root to reproduce the original 149-job workflow. Its main smoother is 5/0.25; its `smooth_time_9` sensitivity does **not** use scale=1 and does **not** recalibrate p-values. It also checks seasonal alternatives, AR/Fourier nulls, COI, drift, thresholds and biwavelet. It writes `results/` and `figures/`; existing `results/legacy/` tables are the archived subset for quick review.

For comparable cluster/episode inference under alternative nulls and smoothing, rerun `R/core_wtc.R` with `WTC_NULL=ar1` or `fourier`, or `WTC_TIME=5 WTC_SCALE=0.25`, each into its own `WTC_OUT`. Keep the calibration and test draw counts unchanged for a final comparison. The original pipeline's pointwise p-values and the newer cluster null are different estimands.
