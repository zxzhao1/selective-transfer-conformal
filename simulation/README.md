# Simulation code

This directory implements the simulation design in
`model/selective_transfer_conformal_prediction.tex`.

## Implemented components

- target/source data-generating mechanism;
- target-only Lasso;
- complete pooling;
- F-AdaTrans with SCAD-derived adaptive weights;
- Ora-F-Ada with known target and contrast supports;
- the proposed outer target-based selection, refitting, and target-only
  conformal calibration;
- empirical coverage, interval length, prediction MSE, slope estimation
  error, and transfer-selection rate;
- scenario grids for Simulations 1--5.

Ora-F-Ada is evaluated only when every true contrast support satisfies
`s_delta < n_source`, the condition required for the oracle estimator in
the reference method. It is omitted from the mixed and dense full-size
Simulation 4 settings, where this condition fails.

## Required packages

```r
install.packages(c("glmnet", "Matrix", "data.table", "ggplot2", "scales"))
```

## Checks

Run the small end-to-end test before any full simulation:

```bash
Rscript simulation/run_smoke_test.R
```

Run one small Simulation 1 scenario:

```bash
Rscript simulation/run_simulation.R \
  --simulation=1 --mode=smoke --replications=1 --scenario=1
```

Every run writes the raw replication-level CSV together with:

- long- and wide-format Monte Carlo summary tables;
- a LaTeX table containing means and Monte Carlo standard errors;
- PDF and PNG figures for coverage, interval length, prediction MSE,
  slope estimation error, and the transfer-selection rate when available.

Existing raw results can be summarized without rerunning a simulation:

```bash
Rscript simulation/summarize_results.R \
  simulation/results/simulation_1_smoke.csv 0.1
```

## Full runs

For example:

```bash
Rscript simulation/run_simulation.R \
  --simulation=1 --mode=full --replications=200 \
  --workers=8 --grid=fixed --seed_offset=10000000 \
  --run_id=fixed_final200
```

For the reference-style pilot that selects the target initial Lasso
penalty first and then tunes the linked transfer penalty on a wider grid:

```bash
Rscript simulation/run_simulation.R \
  --simulation=1 --mode=full --replications=5 --workers=6 \
  --initial_tuning=separate --grid=wide --run_id=separate_wide_pilot
```

For a boundary diagnostic, `--grid=extra_wide` uses multipliers
`0.125, 0.25, 0.5, 1, 2, 4, 8`. This option is intended for pilot runs;
the final grid should be frozen before a formal Monte Carlo study.

Replication-level checkpoints are written under `results/checkpoints`.
Rerunning the same command reuses completed checkpoints by default.

Simulations 4 and 5 require the original `TransLasso-functions.R` from
<https://github.com/saili0103/TransLasso>. The repository does not display
an explicit license, so the author file is not copied into this project.
The implementation was validated against repository commit
`285009edf3e638567474108f823adcc1dc6c7c4a`; the MD5 checksum of the
supplied author file is stored with the fitted TransLasso object.
Pass its local path explicitly:

```bash
Rscript simulation/run_simulation.R \
  --simulation=4 --mode=full --replications=200 \
  --translasso=/absolute/path/to/TransLasso-functions.R
```

## Tuning-grid status

The formal Simulation 1 run uses `--grid=fixed`. It sets the target and
source initial Lasso penalties to one times their theoretical scales and
sets `lambda0` and `lambda1` to one half of their respective theoretical
scales. These constants were selected before the formal run using separate
pilot replications. No inner tuning-grid search is performed when
`--grid=fixed` is supplied. The `standard`, `wide`, and `extra_wide`
options remain available for diagnostics and for later simulations.

The `narrow` option uses multipliers `0.5, 1`. It retains target-based
validation while reducing the unlinked joint transfer grid from 27 to
8 combinations. The Simulation 4 boundary diagnostic showed that this
range is too narrow for heterogeneous and dense contrasts.

The `compact` option uses initial-penalty multipliers `0.5, 1` and a
common transfer-penalty multiplier from `0.5, 1, 2, 4, 8`. The resulting
prespecified joint grid has 10 candidates and retains the broad transfer
range needed by Simulation 4. Target-only uses the two initial-scale
multipliers, whereas Complete pooling uses the five transfer-scale
multipliers.

The transfer-penalty scales are
`sqrt(log(p) / N)` for `lambda0` and
`(n_source / N) * sqrt(log(p) / n_target)` for `lambda1`, matching the
reference implementation's target-sample scaling. Checkpoints include an
implementation version and are recomputed when that version changes.
