# EPU–inflation wavelet coherence: 18 countries

This repository freezes the 18-country monthly EPU/CPI panel and separates the chosen WTC inference from earlier sensitivity analyses. It is a research workflow, not a causal shock identification model. The public repository name contains `esg`, while the variables analyzed here are **EPU** and inflation.

## Quick start

```bash
git clone https://github.com/its-dl-time/esg-inf-wavelet-analysis.git
cd esg-inf-wavelet-analysis
python scripts/unpack_data.py
Rscript -e 'install.packages(c("WaveletComp","biwavelet","urca","jsonlite","viridisLite"), repos="https://cloud.r-project.org")'
Rscript R/core_wtc.R
```

Run the script from the repository root. Default: 18 countries; Morlet `dj=1/20`; WaveletComp Bartlett smoothing `window.size.t=9`, `window.size.s=1`; paired independent, locally heteroskedastic AR(1) surrogate null; 250 independent calibration maps to define per-cell cluster thresholds; 999 fresh maps to test band means, maximum 8-connected cluster excess, and largest qualifying six-month episode mean. Three disjoint bands are `[2,8)`, `[8,32)` and `[32, min(128, floor(N/3))]` months. Calculations use the longest observed uninterrupted EPU–INF overlap per country, inside the native cone of influence. Raw p-values use `(1+exceedances)/(B+1)`.

The full run can take many hours; each country is saved in `results/core/checkpoints/`. A resumed run reuses checkpoints only when the data and configuration bytes match. To test a small set first:

```bash
WTC_COUNTRIES=US,China,Canada WTC_CAL=25 WTC_TEST=49 WTC_OUT=results/smoke Rscript R/core_wtc.R
```

That small run only checks execution. Use the defaults for reported inference. `WTC_NULL=ar1` and `WTC_NULL=fourier` are alternative nulls; give each run a different `WTC_OUT`. `WTC_TIME=5 WTC_SCALE=0.25` restores the original smoother. `results/core/summary.csv` includes p, BH/BY q and Holm p for each metric separately. The adjustments only span countries actually processed; a subset run does **not** represent the 18-country testing family. Do not combine p-values from separate partial runs without recalculating multiplicity across the entire intended family.

## Contents

| Path | Purpose |
|---|---|
| `data/raw/*.csv.gz` | Frozen original EPU/CPI inputs; include excluded countries and Chinese proxies for traceability. |
| `data/processed/core_18.csv.gz` | All 7,290 rows for the 18 retained countries; source columns and calendar dates retained. |
| `data/processed/processed_panel.csv.gz` | Full processed input including alternative Chinese indexes and two excluded countries. |
| `R/core_wtc.R`, `R/wtc_functions.R` | Chosen WTC inference with independent calibration/test surrogates and checkpoints. |
| `R/prepare_data.R` | Rebuilds processed panel from the frozen raw inputs and records the seasonal decision. |
| `robustness/legacy/` | Original 18-country pipeline (smoothing 5/0.25, AR(1), 149 jobs), kept separately. |
| `results/legacy/` | Previously computed results for the original smoother; **not** the new 9/1 core. |
| `extensions/uct/` | US and China conditional WTC with UCT and Brent controls; descriptive extension. |
| `docs/` | Plan, work log, results, limitations and next steps. |

The complete reproduction of the original 18-country baseline is `Rscript robustness/legacy/run_all.R`. It regenerates `results/prepared.rds` and many files under `results/` (including figures); it is computationally expensive. Run `python scripts/validate_data.py` after unpacking to verify frozen snapshots. To rebuild `data/processed/core_18.csv` after updating raw data, rerun `R/prepare_data.R` and filter its `results/processed_panel.csv` to the 18 `series` listed in `docs/data_manifest.csv`; then rerun core. Do not reuse the frozen core CSV after altering the source.

The locally varying volatility surrogate is an AR(1) recursion with a 25-month rolling residual variance profile. It preserves an approximate local volatility envelope, but is **not** an exact nonstationary wavelet phase randomization. Cluster thresholds are calibrated under the same null and smoothing as the test maps. An eight-connected cluster is within one band only; no cluster crosses band boundaries. The selected episode statistic is the maximum mean WTC among contiguous runs at one scale with at least six calendar months above the calibrated cell threshold (also at least 0.70). The null replicates that selection.

For this retrospective dataset, BH across 18×3 country-band hypotheses per statistic is a declared reporting convention for the forthcoming complete run. BY and Holm remain sensitivity checks; correlated tests do not automatically meet BH's positive-dependence conditions. A nominal within-map cluster p is not evidence across the 54 tests by itself. Maps, arrows and episode dates require separate interpretation. Phase lead is modulo the period and does not identify structural shocks.

See [results](docs/RESULTS.md) for what was actually computed and [next steps](docs/NEXT_STEPS.md) for what is pending.
