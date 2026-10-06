# icspack

`icspack` contains the code used for simulation studies and the NHANES oral-health application for evaluating marginal association methods under informative cluster size and informative subgroup size.

## Contents

- `R/` — functions for correlation estimation and weighting methods
- `data/` — processed NHANES and simulated data
- `results/` — simulation and NHANES analysis results
- `figures/` — figures and tables generated for the manuscript
- `scripts/` — data processing, simulation, association testing, correlation analysis, and visualization scripts

The analyses compare standard cluster weighting with subgroup-based weighting approaches, including PPW and OPW, for estimating marginal association between paired clustered outcomes.

## Reproducibility

Run the scripts in numerical order from the root `icspack` directory. Package dependencies are managed using `renv`.