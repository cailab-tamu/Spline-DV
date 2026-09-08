# An analytic form for the Spline-DV 3-D curve

The spline-DV method places every gene at
$(x,y,z) = (\log(1+\mu),\ \log(1+C_\sigma),\ D_r)$ — log mean, log CV, dropout
rate — and fits a smoothing spline through the resulting point cloud
(`sc_splinefit`: 15 pieces, robustness 0.75, arclength parameterization).
Deviation from that curve is the per-condition statistic; the difference of the
two deviation vectors is the DV score.

The curve does not have to be estimated nonparametrically. Under the standard
count model for droplet scRNA-seq it has a closed form with **one** fitted
parameter, and its dropout coordinate has a closed form with **none**.

## How it works, step by step

*A walkthrough from raw counts, using sample 1 of the bundled data. The formal
derivation begins at Section 1.*

**Step 0 — what you have.** A UMI matrix $K$, integers, genes by cells. After
QC, $9{,}654 \times 994$. Nothing else.

**Step 1 — the column sums are the important part.** $L_j=\sum_i K_{ij}$, the
library size of each cell; median 5,011 here. This one vector determines almost
the whole reference curve, which is the part that is not obvious until the
derivation is done.

**Step 2 — normalize, and notice what it does not do.** $X_{ij}=c\,K_{ij}/L_j$
puts deep and shallow cells on a common scale. But **zeros are untouched** —
dividing by a positive number can neither create nor destroy one. That is what
makes the dropout coordinate free.

**Step 3 — three numbers per gene.** Mean, $C_\sigma=\text{std}/\text{mean}$,
and the fraction of cells with a zero count. That is the cloud the method works
in.

**Step 4 — ask what an *ordinary* gene would look like.** Gene $i$ has a
transcriptome share $\lambda_i$, cell $j$ samples $L_j$ molecules:

$$\mathbb E[K_{ij}] = L_j\lambda_i,\qquad
\operatorname{Var}(K_{ij}) = L_j\lambda_i + \phi\,(L_j\lambda_i)^2$$

One number per gene, $\lambda_i$; one number for the whole dataset, $\phi$ —
how much more variable than pure sampling noise.

**Step 5 — the mean coordinate collapses.** $\mathbb E[X_{ij}]=c\lambda_i$, the
*same in every cell*, deep or shallow. Normalization has done its job, and
$\mu_i=c\lambda_i$ simply re-reads $\lambda$.

**Step 6 — the CV coordinate, the key step.** Every cell now has the same mean
but a *different variance*, because a deeper cell measures the gene more
precisely. So the pooled variance is the average of the per-cell variances,
and dividing by $\mu^2$ makes the $\lambda$'s cancel:

$$C_\sigma^2 = \frac{\alpha}{\mu} + \phi,
\qquad \alpha = c\,\overline{(1/L)} = \frac{c}{\text{harmonic mean of }L_j}$$

Read it as two terms. $\alpha/\mu$ is the **counting-statistics floor** — a
rare gene cannot be measured precisely, whatever the biology — and $\phi$ is
the biology on top. Here $\alpha = 1.9973$, straight from the library sizes,
nothing fitted.

**Step 7 — the dropout coordinate is free.** Zeros live on the raw counts, so
$\Pr(K_{ij}=0)=e^{-L_j\lambda_i}$ under Poisson sampling. Averaging over the
cells you actually have,

$$D_r(\mu) = \frac1n\sum_j e^{-L_j\mu/c},$$

the **empirical Laplace transform of your library sizes**. No free parameter at
all: every $L_j$ is already known, so the mean–dropout relation of a dataset is
fixed the moment it is loaded.

**Step 8 — sweep $\mu$ and the curve appears.** Total cost: one fitted number,
$\phi = 0.152$ here, from a robust regression on the gene cloud.

**Step 9 — check it on real genes.**

| gene | mean | $C_\sigma$ | predicted $C_\sigma$ | dropout | predicted dropout |
|---|---|---|---|---|---|
| Try4 | 0.206 | 3.295 | 3.138 | 0.886 | **0.882** |
| Sparc | 0.945 | 1.663 | 1.506 | 0.612 | 0.584 |
| Malat1 | 472.6 | 0.532 | 0.396 | 0.000 | 0.000 |
| Actb | 18.24 | 0.895 | 0.512 | 0.009 | 0.001 |

Dropout is predicted to three decimals with nothing fitted. Try4 sits
essentially *on* the curve in this sample — it is a DV gene because it moves
*off* the curve in sample 2. Actb's observed CV, 0.895, sits well above its
prediction of 0.512: that gap is real biological variability, and finding such
gaps is what the method is for.

**What this buys you.** The spline asks the 9,654 genes where the curve should
go; the analytic form asks the count model, and the genes contribute only one
scalar. So the reference cannot be dragged around by the very genes being
measured, it is defined at every mean rather than only where genes happen to
sit, and it cannot leave the feasible region the way the fitted spline does
when its dropout coordinate goes negative.

## 1. Setup

Let $K_{ij}$ be the raw count of gene $i$ in cell $j$, $L_j=\sum_i K_{ij}$ the
library size, $n$ the number of cells, and $c$ the normalization scale
(`pkg.norm_libsize` uses $c=10^4$, i.e. CP10K). The pipeline analyses

$$X_{ij} = c\,K_{ij}/L_j .$$

Model the counts as gamma-Poisson (negative binomial) with cell-independent
gene rate $\lambda_i$ (an expression *fraction*) and dispersion $\phi$:

$$\mathbb{E}[K_{ij}] = L_j\lambda_i,\qquad
\operatorname{Var}(K_{ij}) = L_j\lambda_i + \phi (L_j\lambda_i)^2 .$$

$\phi\to 0$ is the pure-Poisson (technical-noise-only) limit.

## 2. The three coordinates

**Mean.** $\mathbb{E}[X_{ij}] = c\lambda_i$ for every cell, so

$$\mu_i = c\lambda_i .$$

Library-size normalization has removed the cell-to-cell mean variation — this
is what makes the rest of the derivation clean.

**CV.** Because every cell has the same mean, the pooled variance is the
average of the per-cell variances:

$$\operatorname{Var}(X_{i\cdot})
= \frac1n\sum_j \frac{c^2}{L_j^2}\big(L_j\lambda_i + \phi L_j^2\lambda_i^2\big)
= c^2\lambda_i \overline{\left(\tfrac1L\right)} + \phi\,c^2\lambda_i^2 ,$$

hence

$$\boxed{\;C_{\sigma}^2(\mu) = \frac{\alpha}{\mu} + \phi,
\qquad \alpha \equiv c\,\overline{\left(\tfrac1L\right)} = \frac{c}{n}\sum_j \frac{1}{L_j}\;}$$

$\alpha$ is **not** fitted: it is $c$ divided by the harmonic mean of the
library sizes. On the bundled GSM3308547/8 data, $\alpha = 1.997$ and $1.597$.

*On the divisor.* The step above uses the population variance $\frac1n\sum_j$,
while the pipeline calls `std(X,0,·)`, which divides by $N-1$. That is exactly
what makes $\alpha$ correct rather than approximately correct. For independent
$X_j$ with a common mean but per-cell variances $\sigma_j^2$, and $m$ the
sample mean,

$$\mathbb{E}\Big[\textstyle\sum_j (X_j-m)^2\Big] = \Big(1-\tfrac1n\Big)\sum_j \sigma_j^2
\quad\Longrightarrow\quad
\mathbb{E}\Big[\tfrac{1}{n-1}\textstyle\sum_j (X_j-m)^2\Big] = \frac1n\sum_j \sigma_j^2 ,$$

with no $O(1/n)$ remainder. Simulated on a synthetic Poisson gene at
$\mu=1$ CP10K over 2000 replicates, the mean sample variance is 1.99850
against $c^2\lambda\,\overline{(1/L)} = 1.99729$ — a 0.06% gap that is Monte
Carlo noise. So $\alpha$ needs no finite-sample correction.

**Dropout.** Zeros are invariant to a positive per-cell rescaling, so
$D_r$ can be computed on raw counts. For the gamma-Poisson,
$\Pr(K_{ij}=0) = (1+\phi L_j\lambda_i)^{-1/\phi}$, so

$$\boxed{\;D_r(\mu) = \frac1n\sum_j\Big(1 + \phi_z\,\frac{L_j\mu}{c}\Big)^{-1/\phi_z}
\;\xrightarrow[\ \phi_z\to0\ ]{}\;
\frac1n\sum_j e^{-L_j\mu/c}\;}$$

The Poisson limit is exactly the **empirical Laplace transform of the library
sizes**, evaluated at the gene's rate $\lambda=\mu/c$. It contains no free
parameter at all — it is a deterministic function of the cell library sizes,
which are already known.

## 3. The curve

Parameterized by $\mu>0$:

$$\mathbf{r}(\mu) = \Big(\ \log(1+\mu),\ \
\log\!\big(1+\sqrt{\alpha/\mu+\phi}\big),\ \
\tfrac1n\textstyle\sum_j (1+\phi_z L_j\mu/c)^{-1/\phi_z}\ \Big)$$

with $\alpha$ given in closed form, $\phi$ the single fitted overdispersion,
and $\phi_z=0$ (Poisson dropout) working best in practice — see §5.

**Implicit form.** The mean–CV relation is the cylinder

$$(e^{y}-1)^2\,(e^{x}-1) = \alpha + \phi\,(e^x-1),$$

and in the idealized case of equal library sizes ($L_j\equiv L$, $\phi=0$,
$\alpha=c/L$) the dropout coordinate collapses to a relation in $y$ alone,
free of every parameter:

$$z = \exp\!\big(-(e^{x}-1)/\alpha\big) = \exp\!\Big(-\big(e^{y}-1\big)^{-2}\Big).$$

So the "3-D curve" is, structurally, the intersection of two surfaces; only
library-size heterogeneity and biological overdispersion perturb it.

**Inputs and outputs.** This is a parametric space curve, not a surface: one
independent variable, three dependent ones.

| | |
|---|---|
| independent variable | $\mu>0$, the gene's normalized mean (equivalently $\lambda=\mu/c$) |
| dependent variables | $x(\mu),y(\mu),z(\mu)$ — log mean, log CV, dropout rate |

To instantiate the curve you need the raw library sizes $\{L_j\}$, the scale
$c$, and the two dispersions. Of those, $\alpha$ is a closed form of $\{L_j\}$,
$\phi_z$ is fixed at 0, and **$\phi$ is the only fitted quantity** — the whole
reference curve costs one scalar, and its dropout coordinate costs none.

Per gene, projecting $\mathbf p_i=(x_i,y_i,z_i)$ onto the curve returns
$\hat\mu_i$ (the fitted rate, an interpretable by-product), the foot point
$\mathbf r(\hat\mu_i)$, the deviation $\vec v_i=\mathbf p_i-\mathbf r(\hat\mu_i)$
with $d_i=\|\vec v_i\|$, and hence the DV score
$\|\vec v_i^{(2)}-\vec v_i^{(1)}\|$.

Note the contrast: with a spline, all 9,654 genes are collectively inputs to
the reference curve, so it depends on the sample in an uncontrolled way. Here
the genes enter the reference only through $\phi$.

**Derivatives** (for exact Newton projection, no discretization):

$$x'=\frac{1}{1+\mu},\qquad
y'=\frac{w'}{1+w}\ \ \text{with}\ \ w=\sqrt{\alpha/\mu+\phi},\ w'=-\frac{\alpha}{2\mu^2 w},$$
$$z'=-\frac1n\sum_j \frac{L_j}{c}\Big(1+\phi_z\frac{L_j\mu}{c}\Big)^{-1/\phi_z-1}
\ \xrightarrow[\phi_z\to0]{}\ -\frac1n\sum_j \frac{L_j}{c}e^{-L_j\mu/c}.$$

The foot point of a gene $\mathbf{p}$ solves
$\big(\mathbf{r}(\mu)-\mathbf{p}\big)\cdot\mathbf{r}'(\mu)=0$.

## 4. Estimating $\phi$

$\phi$ can be fitted with a one-dimensional robust regression directly on the
gene cloud — no spline needed anywhere:

$$\hat\phi = \arg\min_\phi \sum_i \big|\log(1+\sqrt{\alpha/\mu_i+\phi}) - \log(1+C_{\sigma,i})\big| .$$

L1 keeps genuinely variable genes from dragging the reference upward, playing
the role the spline's robustness weight (0.75) plays now.

## 5. How well it reproduces the spline

Bundled data (`MATLAB/GSM3308547_GSM3308548.mat`), QC-filtered, gene sets
intersected: 9,654 genes; 994 and 875 cells.

| | sample 1 | sample 2 |
|---|---|---|
| $\alpha$ (closed form) | 1.9973 | 1.5970 |
| $\hat\phi$ (L1 on genes / fitted to spline) | 0.152 / 0.129 | 0.200 / 0.159 |
| $\log(1+C_\sigma)$: RMSE vs. spline (range ≈ 0.3–2.3) | 0.0110 | 0.0132 |
| $\log(1+C_\sigma)$: $R^2$ vs. spline | 0.99934 | 0.99896 |
| $D_r$, Poisson, **0 fitted parameters**: RMSE vs. spline | 0.0071 | 0.0116 |
| $D_r$, Poisson: $R^2$ vs. spline | 0.99868 | 0.99688 |
| $D_r$, Poisson: $R^2$ vs. the 9,654 **genes** | 0.9944 | 0.9915 |

Downstream, replacing the spline with the analytic curve throughout the
`test.m` DV pipeline gives Spearman $\rho = 0.93$ with the spline DV score,
73/100 top-100 overlap and 388/500 top-500 overlap. Runtime is 2–3 s per
sample.

![analytic vs spline, 3-D](analytic_vs_spline_3d.png)

*The gene cloud in the 3-D space the method works in, with the fitted spline
(solid) and the analytic curve (dashed). The shaded surface is $z=0$. The two
curves are indistinguishable over the body of the data; the spline then stops
where the genes thin out, having already left the feasible region, while the
analytic curve stays on the plane and continues across the full range of gene
means.*

![analytic vs spline, coordinate by coordinate](analytic_vs_spline.png)

*The same comparison coordinate by coordinate, which makes the size of the
discrepancy readable.*

## 6. Why the analytic curve is preferable where the two disagree

* **The spline leaves the feasible region.** Its dropout coordinate reaches
  $-0.086$ (sample 1) and $-0.165$ (sample 2) — a negative rate. The largest
  analytic-vs-spline discrepancies are exactly there. The analytic curve is
  a mixture of survival functions, so $0<z\le1$ by construction, and $y$ is
  monotone decreasing as it must be.
* **The curve is defined for all $\mu>0$.** The spline only spans the
  arclength of the observed genes, which is why `sc_splinefit` rescales the
  distance by 1/100 and 1/10 for genes past either end, and why `test.m`
  zeroes out genes whose `nearidx` is the first or last knot. Those patches
  are unnecessary with a closed form. (On these data the spline stops near
  $\log(1+\mu)=3.3$ while genes extend to 7.7.)
* **No hyperparameters and no rank-deficiency.** The 15 pieces, the 0.75
  robustness weight, the arclength parameterization and the
  `MATLAB:rankDeficientMatrix` warning suppression all go away.
* **The reference is comparable across conditions.** Today $\vec v_1$ and
  $\vec v_2$ are measured against two independently estimated splines, so the
  DV score mixes per-gene deviation with sampling noise in the two reference
  curves. Analytically the two references differ only through
  $(\alpha,\phi,\{L_j\})$, all of which are interpretable and estimable.
  This shows up directly in the Procrustes dissimilarity between the two
  conditions' curves, **0.0945 analytically against 0.9776 for the splines**:
  two analytic references stay inside a two-parameter family whatever the data
  does, whereas two independently fitted splines are free to take different
  shapes, and do. Both figures are computed with the two curves paired by gene
  name — Procrustes on the raw fits compares mismatched genes, because each
  fit sorts genes by its own condition's statistics.
* **The residual becomes interpretable.** $\hat\mu$ at the foot point
  estimates the gene's rate, and the CV residual is excess variability over
  the gamma-Poisson expectation rather than over an empirical trend.

## 7. Usage

`sc_analyticfit` is a drop-in replacement for `sc_splinefit`. Both it and
the two-condition wrapper `sc_analyticfit2` live in scGEAToolbox, beside
`sc_splinefit` and `sc_splinefit2`:

```matlab
[T, xyzFit, curve, params] = sc_analyticfit(Xraw, genelist);
v = [T.lgu T.lgcv T.dropr] - xyzFit(T.nearidx, :);   % same pattern as test.m
```

For a two-condition comparison in one call, matching `sc_splinefit2`:

```matlab
[T, sx, sy, sz, d, genesfit, params] = sc_analyticfit2(Xraw, Yraw, genelist);
```

Pass `IsNormalized=true` together with `LibSize=` if the matrix has already
been through `sc_norm`. `Dispersion=0` gives the pure-Poisson technical null;
`Dispersion=` a scalar forces one $\phi$ on both conditions; the default fits
$\phi$ robustly per condition. Projection is a grid search refined
parabolically in $\log\mu$, accurate to $\sim10^{-6}$ against a 150,000-point
brute-force search.

`sc_analyticfit2` returns `sx`, `sy` and `sz` in the shared gene order
`genesfit`, with `d` the Procrustes dissimilarity, matching `sc_splinefit2`.
Prefer `params` to `d` for judging how far apart the two references are:
Procrustes may rotate and rescale, so it measures shape similarity after
alignment rather than separation in the space the DV score lives in. The
direct quantity is

```matlab
mu  = logspace(-3, 2, 500)';
gap = max(abs(params.curve1.y(mu) - params.curve2.y(mu)));   % 0.1092
```

## 8. Should the deeper sample be downsampled first?

The two conditions' reference curves differ mainly through $\alpha$, which is
fixed by sequencing depth (harmonic mean library 5,007 vs 6,262 here). It is
natural to ask whether thinning the deeper sample to match would help.
Running `depth_matching_experiment` (~11 s) answers it — **no**:

* Thinning by $p=\text{HM}_1/\text{HM}_2 = 0.7996$ does converge $\alpha$
  (1.5970 → 1.998 against 1.9973) and shrinks the gap between the two
  references by about 60%, $\max|\Delta y|$ 0.109 → 0.042.
* But two runs of the same procedure, differing only in the random thinning,
  agree at only $\rho = 0.774$, while a thinned run agrees with the
  full-depth run at $\rho = 0.848$. The thinning injects more noise than the
  bias it removes.
* The top 100 genes are stable either way (92–93 shared), so the damage is in
  the bulk of the ranking, where scores are near zero.

The deviation-vector construction already corrects for depth analytically;
downsampling redoes that correction destructively.

The by-product is more interesting than the answer: once $\alpha$ matches,
$\phi$ still does not — 0.152 for sample 1 against 0.204 for sample 2. That
residual can no longer be attributed to sequencing depth, so it is a real
global difference in biological overdispersion between the two conditions.
A fitted spline smears the same information into the shape of a curve, where
nothing can read it back out.

## 9. Reproducing the figures

Both figures above are regenerated from the bundled data by

```matlab
make_analytic_figures            % ~30 s; needs scGEAToolbox on the path
```

which prints $\alpha$ and both estimates of $\phi$, and writes the two PNGs
into this folder. `PhiSource="genes"` switches the figure to the robust L1
estimate `sc_analyticfit` uses (0.152 / 0.200) instead of the spline-matched
value (0.129 / 0.159) plotted here.
