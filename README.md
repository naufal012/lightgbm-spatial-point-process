# LightGBM for Spatial Point Processes (LightGBMPP)

A research codebase implementing **LightGBM-based intensity estimation** for spatial point processes (SPP) and linear point processes (LPP). This project covers both simulation studies and real-world applications using custom Poisson and logistic log-likelihood objectives within the LightGBM framework, with XGBoost available as a comparison model.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Key Differences from XGBoostPP](#key-differences-from-xgboostpp)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Simulation Study](#simulation-study)
  - [SPP Simulation (Planar)](#spp-simulation-planar)
  - [LPP Simulation (Linear Network)](#lpp-simulation-linear-network)
- [Real-World Application — Traffic Accidents (Waejuc/Nganjuk)](#real-world-application--traffic-accidents-waejucnganjuk)
  - [Dataset Description](#dataset-description)
  - [Running the Poisson Application](#running-the-poisson-application)
  - [Running the Logistic Application](#running-the-logistic-application)
- [Loss Functions](#loss-functions)
- [Hyperparameter Tuning](#hyperparameter-tuning)
- [Outputs](#outputs)
- [Data Description](#data-description)
- [Citation](#citation)

---

## Overview

This project estimates the **intensity function** λ(u) of a spatial point process using gradient boosting (LightGBM / XGBoost) with custom objectives derived from:

- **Poisson log-likelihood** — standard Poisson process likelihood
- **Weighted Poisson log-likelihood** — accounts for clustering via an F_prime correction
- **Logistic log-likelihood** — Baddeley-style logistic regression for point processes
- **Weighted logistic log-likelihood** — weighted variant with cluster correction

The model is trained on a **case-control** dataset: observed event points (label = +1) mixed with randomly generated dummy/quadrature points (label = −1). Intensity predictions at dummy locations reconstruct the continuous intensity surface.

Key features:
- **Primary model:** LightGBM with GOSS boosting (`boosting_type = "goss"`)
- **Comparison model:** XGBoost (shared Python module)
- Works on **planar windows** (SPP) and **linear networks** (LPP)
- Supports **fixed parameters** and **Bayesian hyperparameter tuning** via Optuna
- Seamlessly integrates R (`spatstat`) with Python (`lightgbm`, `xgboost`) via `reticulate`
- Real-world application on **traffic accident data** from Nganjuk, East Java (2022–2025)
- Includes MISE and log-likelihood evaluation metrics

---

## Repository Structure

```
lightgbm-spatial-point-process/
│
├── python/                           # Core Python modules
│   ├── lgbpp.py                      # LightGBMPP: custom objectives + Optuna tuner
│   └── xgbpp.py                      # XGBoostPP: custom objectives + Optuna tuner (comparison)
│
├── R/
│   ├── simulation/                   # Simulation study (planar SPP)
│   │   ├── simulate_processes.R      # Simulate Poisson, Thomas, LGCP, Strauss processes
│   │   ├── wrapper_spp.R             # run_analysis(): train + evaluate + visualize per run
│   │   └── run_simulation_spp.R      # Main entry point for SPP simulation loop
│   │
│   ├── lpp/                          # Linear Point Process (LPP) simulation
│   │   ├── simulate_lpp.R            # Simulate LPP on a linear network
│   │   ├── wrapper_lpp.R             # run_analysis_LPP(): LPP-specific train + evaluate
│   │   ├── run_simulation_lpp.R      # Main entry point for LPP simulation loop
│   │   └── linquad.R                 # Linear quadrature weight helpers
│   │
│   └── application/                  # Real-world application pipeline
│       ├── application_poisson.R     # Poisson loss pipeline (Waejuc/Nganjuk accident data)
│       └── application_logistic.R    # Logistic loss pipeline (same dataset)
│
├── data/
│   ├── bci/                          # BCI Forest Plot data
│   │   ├── bci.covars.Rda            # BCI covariate images (elevation, slope, etc.)
│   │   └── tree/                     # BCI tree species point patterns (tree1–tree8)
│   │       ├── bci.tree1.rdata
│   │       ├── bci.tree2.rdata
│   │       └── ...
│   └── accident/                     # Nganjuk traffic accident data (LPP application)
│       ├── nganjuk_linnet_rescaled.rds   # Nganjuk road network as linnet object
│       ├── dataall.csv                   # Cleaned accident records (all years)
│       ├── tanggal1.csv                  # Accident date/time table
│       ├── trafficlights.csv             # Traffic light coordinates
│       ├── bmkg_2018_2020.csv            # BMKG weather data (cleaned)
│       ├── pcf200m_h.rds                 # Pre-computed PCF at 200m bandwidth
│       ├── accident_2022.xlsx            # Raw accident records 2022
│       ├── accident_2023.xlsx            # Raw accident records 2023
│       ├── accident_2024.xlsx            # Raw accident records 2024
│       └── accident_2025.xlsx            # Raw accident records 2025
│
├── output/                           # Auto-generated outputs (gitignored)
│   ├── figures/
│   ├── tables/
│   └── models/
│
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Key Differences from XGBoostPP

| Aspect | LightGBMPP (this repo) | XGBoostPP |
|---|---|---|
| **Primary model** | LightGBM (GOSS boosting) | XGBoost (hist method) |
| **Tree complexity param** | `num_leaves` (leaf-wise) | `max_depth` (level-wise) |
| **Regularisation** | `lambda_l1`, `lambda_l2` | `alpha` (L1), `lambda` (L2) |
| **Speed** | Faster on large datasets | Slightly slower |
| **Application dataset** | Nganjuk traffic accidents (LPP, 2022–2025) | BEI/Ryno trees, Crime, Accidents |
| **Application scripts** | Two separate scripts (Poisson / Logistic) | One universal pipeline |
| **BCI tree data** | 8 species (`bci.tree1`–`bci.tree8`) | Single species (`bci.tree1`) |
| **Metric return** | `(name, value, is_higher_better)` | `(name, value)` |

---

## Requirements

### Python (≥ 3.9)

```
lightgbm>=3.3
xgboost>=1.7
optuna>=3.0
numpy>=1.23
pandas>=1.5
tqdm
```

### R (≥ 4.1)

```r
spatstat
spatstat.geom
spatstat.explore
spatstat.linnet
spatstat.random
reticulate
lightgbm
xgboost
ggplot2
viridis
RColorBrewer
openxlsx
tictoc
dplyr
readr
lubridate
sf
sp
```

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/yourusername/lightgbm-spatial-point-process.git
cd lightgbm-spatial-point-process
```

### 2. Set up a Python environment

```bash
conda create -n lgbm-env python=3.10
conda activate lgbm-env
pip install -r requirements.txt
```

### 3. Install R packages

```r
install.packages(c(
  "reticulate", "spatstat", "spatstat.geom", "spatstat.explore",
  "spatstat.linnet", "spatstat.random", "lightgbm", "xgboost",
  "ggplot2", "viridis", "RColorBrewer", "openxlsx", "tictoc",
  "dplyr", "readr", "lubridate", "sf", "sp"
))
```

### 4. Link R to your conda environment

Add this at the top of any R script before loading `reticulate`:

```r
library(reticulate)
use_condaenv("lgbm-env", required = TRUE)
```

---

## Quick Start

### Minimal example — train LightGBMPP on simulated data

```r
library(reticulate)
use_condaenv("lgbm-env", required = TRUE)

source_python("python/lgbpp.py")
source("R/simulation/simulate_processes.R")
source("R/simulation/wrapper_spp.R")

lgb <- import("lightgbm")
pd  <- import("pandas")

# Load BCI covariates
load("data/bci/bci.covars.Rda")
data(bei)

# Simulate a Poisson process
sim <- simulate_poisson_process(
  covariate_names = c("elev", "grad"),
  intercept    = 0,
  coefficients = c(1, -1),
  bci_covars   = bci.covars,
  bei_window   = Window(bci.covars[[1]]),
  scale_factor = 500,
  n_points     = 2000
)

# Train with fixed parameters
train_lgbpp_fixed <- function(X, y, vol, loss, F_prime, base_params) {
  final_params <- c(base_params, list(
    learning_rate = 0.001, lambda_l1 = 0, lambda_l2 = 0, num_leaves = 63L
  ))
  train_set <- lgb$Dataset(data = as.matrix(X), label = pd$Series(y))
  lgbpp_py(
    data = train_set, vol = pd$Series(vol), params = final_params,
    loss = loss, F_prime = F_prime, num_boost_round = 5000L,
    valid_sets = list(train_set),
    callbacks = list(lgb$early_stopping(stopping_rounds = 50L, verbose = FALSE))
  )
}
```

---

## Simulation Study

### SPP Simulation (Planar)

Evaluates LightGBMPP on three spatial point process types over the BCI forest plot window (1000×500m rescaled to km).

**Entry point:** `R/simulation/run_simulation_spp.R`

#### Supported process types

| Process | Description |
|---|---|
| `Poisson` | Inhomogeneous Poisson with covariate-driven log-intensity |
| `Thomas` | Neyman–Scott cluster process |
| `LGCP` | Log-Gaussian Cox Process |
| `Strauss` | Inhibition process via Metropolis-Hastings |

#### How to run

1. Open `R/simulation/run_simulation_spp.R`
2. Update the conda path to match your system:
   ```r
   Sys.setenv(RETICULATE_PYTHON = "/path/to/envs/lgbm-env/python.exe")
   use_condaenv("lgbm-env", conda = "/path/to/conda.exe", required = TRUE)
   ```
3. Set simulation parameters:
   ```r
   N_SIMULATIONS  <- 50     # Monte Carlo replicates
   num_sim_points <- 2000   # Points per simulation
   scale_factor   <- 500
   ```
4. Run the script. Output is saved per replicate under `output/`.

#### Key functions (from `R/simulation/`)

```r
# Simulate processes
simulate_poisson_process(covariate_names, intercept, coefficients,
                         bci_covars, bei_window, scale_factor, n_points)
simulate_thomas_process(...)
simulate_lgcp_process(...)
simulate_strauss_process(...)

# Run analysis for one replicate
run_analysis(model_type,       # "lgb" or "xgb"
             loss_type,        # "poisson", "logistic", "weighted_poisson", ...
             sim_data, sim_intensity, sim_points,
             base_output_dir, run_number, base_params,
             analysis_type,    # "fixed" or "tuned"
             scale_factor)
```

---

### LPP Simulation (Linear Network)

Extends the framework to **linear point processes** on the Nganjuk road network.

**Entry point:** `R/lpp/run_simulation_lpp.R`

#### How to run

1. Open `R/lpp/run_simulation_lpp.R` and update the conda path
2. The network is loaded from `data/accident/nganjuk_linnet_rescaled.rds`
3. Run the script; outputs are saved per replicate

#### Key differences from SPP

- Observation window is a `linnet` object (road network), not a polygon
- Quadrature weights (`vol`) computed using `linquad.R` for network geometry
- MISE is computed over **total network length**, not area
- Intensity maps are `linim` objects from `spatstat.linnet`
- `F_prime` estimated using `Knet()` instead of `Kinhom()`

---

## Real-World Application — Traffic Accidents (Waejuc/Nganjuk)

The application analyses traffic accident events on the **Nganjuk regency road network** (East Java, Indonesia), spanning 2022–2025. Two separate pipelines are provided for the two loss functions.

### Dataset Description

| File | Description |
|---|---|
| `nganjuk_linnet_rescaled.rds` | Road network as a rescaled `linnet` object |
| `dataall.csv` | Cleaned and merged accident records across all years |
| `tanggal1.csv` | Accident date/time index table |
| `trafficlights.csv` | Traffic signal locations (lon/lat, WGS84) |
| `bmkg_2018_2020.csv` | BMKG weather station data (cleaned) |
| `pcf200m_h.rds` | Pre-computed pair correlation function at 200m bandwidth |
| `accident_20XX.xlsx` | Raw yearly accident location files (2022–2025) |

The data pipeline:
1. Reads and cleans raw accident records from `.xlsx` files
2. Parses Indonesian datetime strings (e.g. "Senin, 3 Januari 2022 Jam 08.00 Wib")
3. Projects coordinates from WGS84 to UTM Zone 49
4. Constructs the case-control dataset on the `linnet` object
5. Computes quadrature weights via linear quadrature

### Running the Poisson Application

```r
# Open R/application/application_poisson.R
# Update the conda path at the top, then source:
source("R/application/application_poisson.R")
```

This pipeline:
- Fits LightGBMPP with **Poisson log-likelihood**
- Optionally also fits XGBoostPP for direct comparison
- Evaluates with Poisson log-likelihood and predicted event count
- Produces intensity maps on the linear network and feature importance plots

### Running the Logistic Application

```r
source("R/application/application_logistic.R")
```

This pipeline:
- Fits LightGBMPP with **logistic log-likelihood** (Baddeley formulation)
- Uses `linquad.R` for precise quadrature weight computation
- Reports both logistic and Poisson log-likelihoods for comparison

---

## Loss Functions

All objectives are implemented in `python/lgbpp.py`. The model predicts `f(u)` on the log-intensity scale, so the estimated intensity is `λ̂(u) = exp(f(u))`.

> **Note:** LightGBM's `feval` callback requires a 3-tuple return `(name, value, is_higher_better)`. This differs from XGBoost's 2-tuple — the two metric implementations are therefore not interchangeable.

| Loss | Function | Best for |
|---|---|---|
| `"poisson"` | Poisson log-likelihood with case-control normalisation | Standard IPP |
| `"weighted_poisson"` | Poisson + F_prime cluster correction | Clustered patterns |
| `"logistic"` | Baddeley logistic log-likelihood | Alternative to Poisson |
| `"weighted_logistic"` | Weighted logistic with cluster correction | Clustered + logistic |

The `F_prime` correction is estimated from the inhomogeneous K-function:

```r
# Planar SPP
k_func  <- Kinhom(points, lambda = intensity, correction = "translation")
r_med   <- median(nndist(points))
F_prime <- with(k_func[which.min(abs(k_func$r - r_med)), ], trans - theo)
F_prime <- max(F_prime, 0)   # clamp to zero if negative

# Linear network (LPP)
k_func  <- Knet(lpp_points, lambda = intensity_linim, correction = "translation")
```

---

## Hyperparameter Tuning

### Fixed parameters (`analysis_type = "fixed"`)

Fast defaults for exploration:

```python
params = {
    "boosting_type"        : "goss",
    "top_rate"             : 0.2,
    "other_rate"           : 0.2,
    "feature_fraction"     : 1/3,
    "max_bin"              : 255,
    "deterministic"        : True,
    "num_threads"          : 1,
    "learning_rate"        : 0.001,
    "lambda_l1"            : 0,
    "lambda_l2"            : 0,
    "num_leaves"           : 63
}
```

### Bayesian tuning (`analysis_type = "tuned"`)

Uses **Optuna** to search `learning_rate`, `lambda_l1`, `lambda_l2`, and `num_leaves` with up to 5000 boosting rounds and early stopping (50 rounds):

```r
# In R (via reticulate)
tuning_results <- tune_lgbpp(
  X_df              = X,
  y_series          = y,
  vol_series        = vol,
  loss              = "poisson",   # or "logistic"
  F_prime           = 0,
  n_trials          = 100L,
  constrain_events  = TRUE,        # penalise if predicted N ≠ observed N
  constraint_strength = 1.0
)

final_model    <- tuning_results$final_model
log_likelihood <- tuning_results$best_log_likelihood
num_events     <- tuning_results$num_events
```

All trial results are saved to `lightgbm/<loss>/optuna_search_results.csv`. Final model training time is saved to `lightgbm/<loss>/final_model_time.txt`.

---

## Outputs

```
output/
├── figures/
│   ├── intensity_<run>.png             # Predicted intensity heatmap (linim / im)
│   ├── simulated_intensity_run_*.png   # Per-run simulation intensity
│   └── importance_*.png                # Feature importance bar chart
├── tables/
│   ├── summary_results.xlsx            # Per-run metrics (log-lik, MISE, time)
│   └── results_dataframe.csv           # Point-level predictions + covariates
└── models/
    └── (saved model objects if applicable)
```

Per-run Excel summaries (`summary_results.xlsx`) contain:

| Metric | Description |
|---|---|
| Log-Likelihood | Model Poisson log-likelihood |
| Logistic Log-Likelihood | Model logistic log-likelihood |
| True Log-Likelihood | Oracle log-likelihood from true intensity |
| MISE | Mean Integrated Squared Error (scaled & unscaled) |
| Predicted Events | ∫ λ̂(u) du over observation window |
| True Events | Observed event count |
| Computation Time | Wall-clock time in minutes |

---

## Data Description

### BCI Forest Plot (`data/bci/`)

- `bci.covars.Rda` — covariate rasters: `elev`, `grad`, `aspect`, `convex`, `beers`
- `bci.tree1.rdata` – `bci.tree8.rdata` — point patterns for 8 different BCI tree species

### Nganjuk Traffic Accidents (`data/accident/`)

Traffic accidents recorded by the Nganjuk Police (Polres Nganjuk), East Java. Records span 2022–2025 and include:
- Accident location (latitude/longitude, then projected to UTM Zone 49S)
- Date and time (parsed from Indonesian-language strings)
- Linked to the Nganjuk road network (`nganjuk_linnet_rescaled.rds`)

Weather covariates from `bmkg_2018_2020.csv` (BMKG meteorological stations) and road infrastructure features from `trafficlights.csv` are incorporated as spatial covariates.

---

## Citation

If you use this code in your research, please cite:

```
[Your paper citation here]
```

---

## Notes

- All Indonesian-language comments (`#` lines in Bahasa Indonesia) in the source code are original annotations from the development phase.
- The `output/` directory is gitignored; all result files are generated locally at runtime.
- Both `lgbpp.py` and `xgbpp.py` must be present for the R scripts to run, since simulation wrappers train and compare both models.
