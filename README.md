# Pressure-driven multispecies tumor growth

**A decoupled structure-preserving scheme for pressure-driven multispecies tumor growth**  
Weijie Huang and Xinran Ruan.

MATLAB code and numerical data accompanying the article.

## Setup

Use MATLAB R2023a (9.14.0.2206163). Keep `code/`, `data/` and this README
together, and set the MATLAB Current Folder to the repository root.

```text
README.md
code/
    redraw_paper.m
    run_paper_case.m
    startup_paper.m
    validate_release.m
    +model/                 Model parameters and initial conditions
    +src/                   Numerical solvers and experiment drivers
    +utils/                 Plotting, data I/O and postprocessing
data/
    Figure1/ ... Figure6/
    Table1/  ... Table3/
```

The code uses MATLAB packages: `model` defines the model and cases, `src`
advances the numerical solution, and `utils` processes and plots the results.
`data/` contains numerical MAT and CSV files with parameter JSON files.

## Reproduce the figures and tables

From the repository root, run:

```matlab
addpath(fullfile(pwd, 'code'));
out = fullfile(tempdir, 'two-species-results');
redraw_paper('all', out);
```

The command generates all six figures and three tables from the supplied
data. Results are saved under `out`. Individual results can be selected by
number or name:

```matlab
redraw_paper(1, out);
redraw_paper('Figure6', out);
redraw_paper('Table1', out);
redraw_paper('tables', out);
```

| Result | Data files, relative to `data/` | Content | Command |
| --- | --- | --- | --- |
| Figure 1 | `Figure1/plot_data.mat`, `temporal_errors.csv`, `spatial_errors.csv` | Temporal and spatial convergence for smooth and near-vacuum data | `redraw_paper(1,out)` |
| Figure 2 | `Figure2/plot_data.mat`, `structure_history.mat` | Total density, common pressure, modified energy and total mass | `redraw_paper(2,out)` |
| Figure 3 | `Figure3/plot_data.mat`, `production_history_m2.mat`, `production_history_m3.mat` | Two- and three-species density profiles and mass fractions | `redraw_paper(3,out)` |
| Figure 4 | `Figure4/switching_data.mat` | Space-time evolution of the type-2 fraction | `redraw_paper(4,out)` |
| Figure 5 | `Figure5/stiff_pressure_data.mat` | Density and pressure profiles with Hele--Shaw references | `redraw_paper(5,out)` |
| Figure 6 | `Figure6/gamma20_data.mat`, `gamma3_data.mat` | Three-species composition for `gamma=20` and `3` | `redraw_paper(6,out)` |
| Table 1 | `Table1/cost.csv`, `cost_g3.csv`, `cost_g20.csv` | Temporal errors, median computation times and speedup | `redraw_paper('Table1',out)` |
| Table 2 | `Table2/stiff_pressure_data.mat`, `table2_errors.csv` | Field errors and mass-based boundary/interface errors | `redraw_paper('Table2',out)` |
| Table 3 | `Table3/stiff_pressure_data.mat`, `table3_structure.csv` | Complementarity and interior pressure-equation defects | `redraw_paper('Table3',out)` |

Figures are exported as PNG and PDF. Figure 6 is exported as separate
`gamma=20` and `gamma=3` rows. Tables are exported as CSV and LaTeX.
Tables 2 and 3 are calculated from the saved numerical fields. Table 1
uses median times from three trials per method and exponent.

An alternative data directory can be supplied as the third argument:

```matlab
redraw_paper('all', out, another_data_root);
```

## Run the numerical experiments

Use `run_paper_case` to select an experiment and an output directory:

```matlab
addpath(fullfile(pwd, 'code'));
runs = fullfile(tempdir, 'two-species-runs');
out = fullfile(tempdir, 'two-species-results');
run_paper_case('Figure3', runs);
redraw_paper('Figure3', out, runs);
```

| Case ID | Experiment | Outputs |
| --- | --- | --- |
| `Figure1` | Temporal and spatial accuracy tests | Reference/test states, error CSV files and plotting data |
| `Figure2` | Two-species structure test | Density, pressure, diagnostic histories and plotting data |
| `Figure3` | Two- and three-species comparison | Field snapshots, full-time diagnostics and plotting data |
| `Figure4` | Pressure-induced switching | Species, pressure and composition histories |
| `Figure5` | Fixed-grid stiff-pressure test | Fields and diagnostics for Figure 5 and Tables 2--3 |
| `Figure6Gamma3` | Two-dimensional growth with `gamma=3` | Three-species fields and snapshots |
| `Table1` | Efficiency benchmark | Three timing trials per method and exponent |

Each case writes to its corresponding folder under `runs`. Figures 1--3
include `plot_data.mat` for direct use with `redraw_paper`. Figure 5 also
writes the shared fields to `Table2/` and `Table3/`. `Figure6Gamma3` writes
the `gamma=3` data for the corresponding row of Figure 6. New timing
measurements can be formatted with:

```matlab
run_paper_case('Table1', runs);
redraw_paper('Table1', out, runs);
```

The available cases and their settings are defined in `model.paper_cases`.

## Parameters and data format

All variables are dimensionless. The common homeostatic pressure is
`p_H=1`, and the common mechanical pressure is `p=(sum_i n_i)^gamma`.

- **Figure 1:** `gamma=20`, `T=0.2`. Temporal tests use `N=1000`,
  `dt=0.008/2^k` for `k=0,...,6`, and a PC reference with `dt=5e-6`.
  Spatial tests use the first-order scheme on
  `N=100,200,400,800,1600` with `dt=1e-5`. Fine-grid states are
  conservatively averaged onto the coarse grid for spatial comparisons.
- **Figures 2 and 3:** domain `(-10,10)`, `gamma=20`, `N=2000`,
  `dt=1e-4`, `T=2`. Growth coefficients are `10,2` for two species and
  `10,10,2` for the three-species split. Figure 2 uses
  `t=0,0.8,1.5,2`; Figure 3 shows density profiles at `t=1` and
  species mass fractions over `0<=t<=2`.
- **Figure 4:** domain `(-10,10)`, `gamma=20`, `N=2000`, `dt=1e-4`,
  `T=2`. Both species have net growth `6(1-p)`; switching rates are
  `p` and `0.25(1-p)`. Initial fractions are `1/2`. Composition is
  displayed where `n>1e-4`, with the pressure contour `p=0.8`.
- **Figure 5 and Tables 2--3:** domain `(-5,5)`, `N=6400`, `dt=1e-4`,
  `T=1`, growth `1-p` for both species, and no switching. Saved exponents
  are `10,20,40,80,160,320,640`; Figure 5 selects `10,20,160`.
  The reference uses `R0=1`, `X0=0.25`. The mass-based estimates are
  `R_gamma=(h/2)*sum(n)` and `X_gamma=(h/2)*sum(n1-n2)`.
  The interior pressure defect uses the region where the limiting
  reference pressure is at least `0.2`.
- **Figure 6:** domain `(-10,10)^2`, a `500 x 500` grid, `dt=0.01`,
  and growth coefficients `8,4,1.5`. The two rows use `gamma=20` and `3`
  at `t=0,1,2.5`. Color indicates species fractions, opacity is `n^0.75`,
  white indicates `n<0.02`, and the black contour is `n=0.5`.

`plot_data.mat` contains panel labels, curve coordinates and plotting
styles. The `panelNN_curveNN.csv` files provide individual numerical
curves as two columns, `x,y`, without a header. Other MAT files contain
grids, numerical fields, snapshots, time histories and parameter structs.
Figure 3 includes 21 field snapshots on `0:0.1:2` and diagnostics at every
time step. Parameter JSON files provide the case settings in text form.
The identifier `dvd` in internal solver and data fields denotes the
predictor--corrector (PC) scheme. Computation times depend on hardware
and software.
