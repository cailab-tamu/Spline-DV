# MATLAB

To use the functions, please first install [scGEAToolbox](https://www.mathworks.com/matlabcentral/fileexchange/72917-scgeatoolbox-single-cell-gene-expression-analysis-toolbox).

Then, from this folder:

```matlab
test          % runs both methods on the bundled data and compares them
```

## Scripts

| | |
|---|---|
| [`test.m`](test.m) | The main demo, ~5 s. Splits the bundled data into its two samples, builds a reference curve for each with **both** methods, scores every gene, and compares the two rankings. Ends by printing the `test_gui` call. |
| [`test_gui.m`](test_gui.m) | Interactive view of a `test.m` result. Left panel is the two gene clouds with their reference curves; click a gene to see its expression profile across the cells of each sample. Call it with the arguments `test.m` prints. |
| [`make_analytic_figures.m`](make_analytic_figures.m) | Rebuilds the two figures in `../docs/`, ~30 s. Prints $\alpha$ and both estimates of $\phi$ as it goes. `PhiSource="genes"` switches to the estimate `sc_analyticfit` itself uses. |
| [`depth_matching_experiment.m`](depth_matching_experiment.m) | Answers whether the deeper sample should be downsampled so the two reference curves coincide, ~11 s. Answer: no — the thinning injects more noise than the bias it removes. Two seeds, so the noise floor is measurable. `help` gives the finding without running it. |
| `GSM3308547_GSM3308548.mat` | The bundled two-sample dataset every script above uses. GEO samples GSM3308547 and GSM3308548; 11,235 genes, 1,869 cells; mouse pancreas, with both islet (`Ins1`, `Ins2`, `Gcg`, `Sst`) and acinar (`Try4`, `Ctrb1`, `Clps`) populations. Loads as a `SingleCellExperiment`. |

## Where the fitting functions live

Not here — they are part of scGEAToolbox:

| | |
|---|---|
| `sc_splinefit`, `sc_splinefit2` | the fitted smoothing spline, as published |
| `sc_analyticfit`, `sc_analyticfit2` | the closed-form curve, a drop-in replacement returning identically shaped and named tables |

The derivation is in [`../docs/analytic-curve.md`](../docs/analytic-curve.md).
