# Target-Calibrated Selective Transfer Learning

This repository contains the simulation and real-data analysis code for
**Target-Calibrated Selective Transfer Learning for High-Dimensional
Prediction Intervals**.

The proposed procedure constructs a feature-wise adaptive transfer candidate
and a Target-only candidate, compares them using cross-fitted target residual
quantiles, and calibrates the selected predictor with an independent target
calibration sample.

## Repository structure

- `simulation/`: data-generating mechanisms, estimators, formal simulation
  runners, validation checks, final replication-level results, and figures.
- `real_data/`: preparation, analysis, and plotting scripts for the nine ADME
  endpoints from Therapeutics Data Commons (TDC), together with the final
  numerical results and figures.
- `external/`: instructions for obtaining the official TransLasso comparison
  code. That code is not redistributed here because its source repository
  does not provide an explicit license.

Large downloaded or generated data objects are intentionally excluded. The
ADME preparation scripts download the public TDC datasets and reconstruct all
processed inputs.

## Software requirements

The analyses were developed with R and Python. Install the required R
packages with:

```r
install.packages(c(
  "glmnet", "Matrix", "data.table", "ggplot2", "gridExtra", "scales"
))
```

Install the Python dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

## Validation and smoke tests

From the repository root, run:

```bash
Rscript simulation/run_validation_tests.R
Rscript simulation/run_smoke_test.R
```

These checks exercise the main estimators and a small end-to-end simulation
before a formal Monte Carlo run is started.

## Formal simulation studies

Simulations 1--3 use fixed penalty multipliers selected before the formal
replications:

```bash
Rscript simulation/run_simulation.R --simulation=1 --mode=full \
  --replications=200 --workers=8 --grid=fixed \
  --seed_offset=10000000 --run_id=fixed_final200

Rscript simulation/run_simulation.R --simulation=2 --mode=full \
  --replications=200 --workers=8 --grid=fixed \
  --seed_offset=20000000 --run_id=fixed_final200

Rscript simulation/run_simulation.R --simulation=3 --mode=full \
  --replications=200 --workers=8 --grid=fixed \
  --seed_offset=30000000 --run_id=fixed_final200
```

Simulations 4--5 use method-specific parameters fixed from the first ten
pilot replications. Obtain `TransLasso-functions.R` as described in
`external/README.md`, then run:

```bash
Rscript simulation/run_simulation.R --simulation=4 --mode=full \
  --replications=200 --workers=8 --seed_offset=40000000 \
  --fixed_params=simulation/results/simulation_4_pilot10_fixed_parameters.csv \
  --translasso=/absolute/path/to/TransLasso-functions.R \
  --run_id=pilotfixed_final200

Rscript simulation/run_simulation.R --simulation=5 --mode=full \
  --replications=200 --workers=8 --seed_offset=60000000 \
  --fixed_params=simulation/results/simulation_5_pilot10_fixed_parameters.csv \
  --translasso=/absolute/path/to/TransLasso-functions.R \
  --run_id=pilotfixed_final200
```

The `simulation/results/` directory contains the final replication-level CSV
files and the figures reported in the manuscript.

## ADME real-data analysis

The primary application uses nine continuous ADME endpoints from TDC. From
the repository root, run:

```bash
cd real_data
python3 22_prepare_tdc_adme.py
Rscript 23_finalize_tdc_adme.R
Rscript 02_run_pilot.R 10 300 adme 40
Rscript 14_run_gas_sensor_translasso.R 10 300 40 \
  /absolute/path/to/TransLasso-functions.R adme
Rscript 26_summarize_adme_analysis.R
Rscript 27_plot_adme_analysis.R
Rscript 28_plot_adme_prediction_intervals.R
```

The preparation stage canonicalizes SMILES, aggregates repeated records for
the same compound, removes target compounds from every source, constructs
1024-bit Morgan fingerprints, and locks the 18 target-source tasks before
model fitting. The reported analysis uses 40 target fitting compounds, 300
target-variance-selected fingerprint features, and ten repetitions.

## Data availability

The ADME datasets are publicly available through
[Therapeutics Data Commons](https://tdcommons.ai/). The script
`real_data/22_prepare_tdc_adme.py` downloads the named datasets through the
PyTDC interface and reconstructs the analysis inputs.

## Code availability

The code used to reproduce the simulation studies and real-data analysis is
publicly available in this repository.

## Contact

Xun Zhao, School of Science, Jimei University, Xiamen 361021, China  
Email: zhaoxun@jmu.edu.cn

