# Spline-DV

A spline-based scRNA-seq method for identifying differentially variable (DV) genes across two experimental conditions.

Read the paper here: [https://www.nature.com/articles/s41540-025-00507-z](https://www.nature.com/articles/s41540-025-00507-z)

## The reference curve: fitted, or analytic

Spline-DV places every gene at $(\log(1+\mu),\ \log(1+C_\sigma),\ D_r)$ — log mean, log CV, dropout rate — and measures how far it sits from a reference curve drawn through that cloud. The published method fits that curve as a smoothing spline. It also has a closed form.

Under a gamma-Poisson count model, with $c$ the normalization scale and $L_j$ the raw library sizes, the curve is, parameterized by the normalized mean $\mu$:

$$x(\mu)=\log(1+\mu),\qquad
y(\mu)=\log\big(1+\sqrt{\alpha/\mu+\phi}\big),\qquad
z(\mu)=\frac1n\sum_j e^{-L_j\mu/c}$$

$\alpha = \frac{c}{n}\sum_j 1/L_j$ is fixed by the library sizes, not fitted. The dropout coordinate is exactly the **empirical Laplace transform of the library sizes** and carries no free parameter at all, so the entire reference curve costs a single fitted scalar: the overdispersion $\phi$.

On the bundled data it reproduces the fitted spline at $R^2 = 0.999$ in $\log(1+C_\sigma)$ and 0.997–0.999 in dropout rate, and the DV gene ranking at Spearman $\rho = 0.93$. Where the two disagree the spline is the one at fault — its dropout coordinate descends to $-0.086$ and $-0.165$, outside the range a rate can take — and the analytic curve needs none of the out-of-range patching the fitted one requires.

Full derivation, validation and figures: **[docs/analytic-curve.md](docs/analytic-curve.md)**, also as Word with typeset equations in [docs/Spline-DV_analytic_curve.docx](docs/Spline-DV_analytic_curve.docx).

## Getting started

Install [scGEAToolbox](https://www.mathworks.com/matlabcentral/fileexchange/72917-scgeatoolbox-single-cell-gene-expression-analysis-toolbox), then from `MATLAB/`:

```matlab
test          % runs both methods on the bundled data and compares them
```

## Contents

| | |
|---|---|
| [`MATLAB/test.m`](MATLAB/test.m) | runs both methods side by side and compares the gene rankings |
| [`MATLAB/test_gui.m`](MATLAB/test_gui.m) | interactive view: click a gene to see its expression profile in both conditions |
| [`MATLAB/make_analytic_figures.m`](MATLAB/make_analytic_figures.m) | regenerates the figures in `docs/` |
| [`MATLAB/depth_matching_experiment.m`](MATLAB/depth_matching_experiment.m) | should the deeper sample be downsampled first? (no — and why) |
| `example_data/` | DE and DV result tables for the paper's three case studies (obesity, fibrosis, colorectal cancer) |
| `R/` | pointer to the [R implementation](https://github.com/Xenon8778/SplineDV) |

The fitting functions themselves live in scGEAToolbox: `sc_splinefit` and `sc_splinefit2` for the fitted spline, `sc_analyticfit` and `sc_analyticfit2` for the closed form. The analytic pair is a drop-in replacement, returning identically shaped and named tables.
