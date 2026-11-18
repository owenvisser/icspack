# icspack analysis — rebuilt around `clustPear`

This is a repository for the analysis done in the manuscript titled *Estimating the marginal association between paired outcomes in clustered data with application to caries and periodontal disease from the National Health and Nutrition Examination Survey* by Owen Visser, Bruce Dye, Shaoping Zhang, Victor Agboli, Steven Levy, and Somnath Datta. The correspondent for this repository is Owen Visser (email: owviss@gmail.com). `clustPear` is the sole computational engine for NW, CW, EQ and ER Pearson correlations and their cluster-sandwich delta-method standard errors.

## Statistical specification

The source of truth of methods is supplied by the author via email or can be seen within the code. Briefly:

- NW: unit weights.
- CW: inverse observed cluster-size weights.
- OPW: within each cluster, observed joint dichotomy `p` receives weight `1 / (N_iP * n_i(p))` per observation.
- ER: mean of OPW-weighted Pearson correlations over all distinct paired dichotomies, or the stated deterministic capped reference if exhaustive evaluation exceeds 100,000 pairs.
- EQ: mean over the selected distinct-rank quantile grid, with `K = 2, 3, 5, 10, 20, 30, 50, 100` in simulations and `K = 20` for NHANES.
- Uncertainty: `clustPear` supplies independent-cluster delta-method sandwich SEs. For EQ and ER, influences are averaged before squaring, treating observed cutpoints as fixed. Normal bounds are the package's `estimate +/- 1.96 * SE`.

### NHANES note on survey weights and standard errors

Four-year MEC weights are sampling weights in the estimating equations through `survey_weights = "WTMEC4YR"`. The reported SEs treat NHANES participants as independent clusters and do not adjust for the NHANES masked survey strata and PSU design. They therefore must not be described as NHANES-complex-design-valid SEs. Document this limitation explicitly if using the intervals or standard errors for inference about the NHANES target population.

## Installation

Install the author's tested `clustPear` 0.1.2 package separately, plus the necessary R packages:

```r
install.packages(c("Rcpp", "tibble", "here", "foreign", "dplyr", "tidyr", "stringr", "ggplot2", "remotes"))
remotes::install_github("owenvisser/clustPear")
```

The C++ statistical implementation comes only from installed `clustPear`, so it is kept fast and generic.

Run all scripts from the root of this repository. The empty `.here` marker pins relative paths.

## NHANES Analysis

No raw data are distributed in this repository. All NHANES data is publicly available on their website linked here: `https://wwwn.cdc.gov/nchs/nhanes/default.aspx`

```r
source("scripts/01_import_nhanes_raw.R")
source("scripts/02_score_nhanes_teeth.R")
source("scripts/03_prepare_nhanes_analysis.R")
source("scripts/04_nhanes_clustPear_delta.R")
source("scripts/08_make_nhanes_tables.R")
source("scripts/10_plot_NHANES_F2.R")
```

The tooth-level processing is the dental scoring (age `RIDAGEYR >= 30`, full oral examination, at least 10 present teeth, third molars excluded). It scores DS/DSI/FS/FSI as an average of available surfaces (with I being interproximal); after conversion, `DFS=(DS+FS)/2` and `DFSI=(DSI+FSI)/2` as used in the manuscript application.

The NHANES analysis is restricted to adults aged 30 years or older who completed the oral health examination and had at least 10 present teeth, excluding third molars. Tooth-level periodontal and caries measurements are obtained using the original icspack dental scoring procedures. For each periodontal–caries outcome pairing, observations with incomplete paired measurements are excluded, and participants with fewer than two complete tooth pairs for that pairing are omitted. NW, CW, EQ, and ER are then estimated using `clustPear`, with four-year MEC sampling weights and cluster-sandwich delta-method standard errors.

The processed data are in `data/nhanes_teeth.rds`; the estimates (including SEs) are in `results/nhanes/correlations_delta.csv`.

## Simulation Analysis

The package's `simulate_clustered()` provides the manuscript bivariate random-effects outcomes and the retention probability `exp[-eta * expit(-X*) * expit(-Y*)]`; only observed pairs are returned. Clusters with fewer than two observed units are excluded before estimation. All four Pearson methods use `clustPear` directly to simulate and compute the data.

To run a very small pilot in a separate folder:

```r
Sys.setenv(SIM_MODE = "pilot")
source("scripts/05_run_simulation.R")
source("scripts/06_summarize_simulation.R")
```

For the manuscript's 9 settings × 1,000 replicates × 8 K-values with 100 clusters and 100 potential units per cluster:

```r
Sys.setenv(SIM_MODE = "full")
source("scripts/05_run_simulation.R")
source("scripts/06_summarize_simulation.R")
source("scripts/07_make_simulation_tables.R")
source("scripts/09_plot_simulation_F1_S1.R")
```

The simulation writes checkpoints under `results/simulation/full/chunks/` and provides an update every 5 chunks
```r
Sys.setenv(SIM_MODE = "full")
source("scripts/05_monitor_simulation.R")
```

The simulation follows the seed `20261103L` as in R/simulation.R

## Outputs used in manuscript

| Item | Exact output |
|---|---|
| **T1** simulation mean Pearson (mean delta SE), K=20 | `manuscript/generated/full/T1_simulation_estimate_meanSE.tex` |
| **T2** CAL vs DFS, survey-weighted estimate (delta SE) | `manuscript/generated/nhanes/T2_NHANES_CAL_DFS.tex` |
| **T3** PD vs DFS, survey-weighted estimate (delta SE) | `manuscript/generated/nhanes/T3_NHANES_PD_DFS.tex` |
| **F1** pooled MAE vs K | `figures/full/F1_pooled_MAE.pdf` |
| **S1** nine-panel setting-specific MAE | `figures/full/S1_MAE_3x3_grid.pdf` |
| **F2** black-and-white CAL/PD distributions by DFS presence | `figures/nhanes/F2_NHANES_DFS_only_CAL_PD_pattern.pdf` |
| Supplement: NW/CW/ER MAE, nine scenarios | `manuscript/generated/full/S_MAE_reference_values.tex` |
| Supplement: EQ MAE per K, nine scenarios | `manuscript/generated/full/S_MAE_EQ_values.tex` |
| Supplement: numerical F1 pooled MAE | `manuscript/generated/full/S_MAE_pooled_values.tex` |
| Supplement: CAL associations for DS, DSI, FS, FSI and DFSI | `manuscript/generated/nhanes/S_NHANES_CAL_other_caries.tex` |
| Supplement: PD associations for DS, DSI, FS, FSI and DFSI | `manuscript/generated/nhanes/S_NHANES_PD_other_caries.tex` |

All figures have PNG previews as well as vector PDF for LaTeX. Tables are ordinary `table`/`table*` and `tabular` environments (no `longtable`). LaTeX tables need `booktabs` and `multirow` and the figure snippets need `graphicx` in the manuscript preamble.

