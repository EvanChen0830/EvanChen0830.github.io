---
title: "Cross-Sectional Predictive Geometry, Part 2: Feature/Label Preprocessing and Horizon Bands"
description: "Building the label mechanism space starts here: winsorizing and standardizing features so quality scores are comparable across the library, the transform-family trick for seeing nonlinear relationships through a linear tool, and why horizon bands should default to non-overlapping incremental returns."
date: 2026-09-14
series: ["Cross-Sectional Predictive Geometry"]
part: 2
tags: ["cross-sectional", "methodology", "quant-research", "predictive-geometry", "feature-engineering"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Predictive Geometry" >}}

# Part Two: Building the Label Mechanism Space

This part builds $\Phi$ — the coordinate system every later chapter works
in. Build order: feature preprocessing (§4) → label preprocessing (§5) →
horizon banding (§6) → whitening and rank control (§7) → the three liftings
(§8) → hierarchical orthogonalization (§9) → the fit/transform boundary
(§10).

## §4. Feature Preprocessing

### Motivation

Feature-side preprocessing serves two purposes that need to be understood
separately:

1. **Robustness**: removing obvious bad prints and numerical errors so
   downstream moment estimates aren't dominated by a single data error.
2. **Comparability**: making $Q_p=\|c_p\|^2$ directly comparable across
   features. This requires every feature to have the same scale under the
   weighted inner product. Without this step, a feature priced in basis
   points and one priced in percent would differ in $Q$ by a factor of
   $10^4$ — which has nothing to do with predictive power.

### §4.1 Definition

For each date $d$ and feature $p$, take the valid observations
$\{X_{idp}:i\in\mathcal U_d\}$.

**Step one: winsorize.** Default is cross-sectional 1% / 99%:

$$X'_{idp}=\operatorname{clip}\big(X_{idp},\ q_{0.01,d,p},\ q_{0.99,d,p}\big).$$

**Step two: within-date standardization.**

$$\mu_{dp}=\frac1{N_d}\sum_iX'_{idp},\qquad \sigma_{dp}=\sqrt{\frac1{N_d}\sum_i(X'_{idp}-\mu_{dp})^2},$$

$$\boxed{X^{XS}_{idp}=\frac{X'_{idp}-\mu_{dp}}{\sigma_{dp}+\epsilon}.}$$

**Step three: global normalization.** Missing-value imputation and coverage
differences will leave different features with different actual weighted
norms (see §4.2), so a global rescale is applied before entering §11:

$$\|X_p\|_W^2=\frac1T\sum_d\frac1{N_d}\sum_i\big(X^{XS}_{idp}\big)^2,\qquad\boxed{X_p\leftarrow\frac{X_p}{\|X_p\|_W}.}$$

After step three,

$$\boxed{X_p^\top WX_p=1\quad\text{(exactly)}}$$

holds for every $p$ — this is the precondition for $Q_p\le1$ becoming a
strict inequality in §16.

### §4.2 Missing Values

Record daily coverage

$$m_{dp}=\frac{\#\{\text{valid names}\}}{N_d}.$$

Rules:

- $m_{dp}<0.8$: the feature does not contribute on that date; the date is
  recorded in $\mathcal D_p^c$ (the absence set).
- Other missing values are filled with $0$ in standardized space
  (cross-sectionally neutral).
- The report records $\text{coverage}_p=\operatorname{median}_dm_{dp}$ and
  $\text{active\_ratio}_p=T_p/T$, where $T_p=|\mathcal D_p|$ is the number
  of valid days.

**Average over valid days, not the full period.** §11's definition is

$$C_{p,:}=\frac1{T_p}\sum_{d\in\mathcal D_p}C^{(d)}_{p,:},$$

i.e. absent days are not counted as zero in the numerator. The reason:
$Q_p$ should measure "how strong this feature is when it's available," and
"how often it's available" is a separate operational question — the two
should not be multiplied together and conflated. The cost is that different
rows of $C$ are averaged over different date sets, so technically the SVD
in §12 is not the decomposition of a single linear operator. When the
library's coverage distribution is very dispersed (e.g. half the features
have 100% coverage and half have 50%), low-coverage features should be
pooled separately and run through the full pipeline on their own.

**Coverage's hidden screening effect deserves caution.** A feature like
"20-day analyst estimate revision" might have analyst coverage for only
about 60% of the universe, and under the $m_{dp}<0.8$ rule it gets excluded
almost every day, giving $Q_p\to0$ — not because it lacks predictive power,
but because of the threshold. Features like this should go into a
"restricted universe" sub-study before entering the pipeline, rather than
compete in the same pool as full-universe features.

### §4.3 Transform Family

**Motivation.** $c_p=X_p^\top W\Phi$ is linear in $X_p$. So if a feature's
relationship to future returns is U-shaped (e.g. extreme valuation
deviations in either direction predict reversion), then
$\mathbb E[X\cdot z]\approx0$, and this method reports $Q_p\approx0$ despite
the feature carrying complete predictive information.

The fix is not to abandon linearity (that would destroy the entire
geometric architecture) — it's to **move the nonlinearity to the feature
side**: every semantic feature is expanded into a set of transforms, each
of which runs through the full pipeline as its own library entry.

**Definition.** For each source feature $p$, define the transform family

$$\boxed{\mathcal T_p=\{t_1,t_2,\dots\},\qquad X_{p,k}=\text{standardize}\big(t_k(X^{XS}_p)\big).}$$

The suggested minimal starting family has three:

$$t_1(x)=x,\qquad t_2(x)=|x|,\qquad t_3(x)=x|x|.$$

Their roles divide as: $t_1$ captures a monotonic linear relationship; $t_2$
captures "degree of deviation" (symmetric, directionless) — exactly the
linearization of a U-shaped relationship; $t_3$ is an odd function with
amplified convexity, capturing "the directional effect of extreme values is
stronger than mid-range values." Note that $t_2$'s output must be
**re-demeaned and re-standardized within-date**, or it will carry a large
day-level common component (that day's cross-sectional dispersion), which
by Lemma 4.1 contributes nothing to $C$ but wastes $X$'s variance.

When a richer function family is needed, a rank-space B-spline basis (3–5
knots) can be used instead, at the cost of each source feature producing
3–5 entries, inflating library size and multiple-testing burden
accordingly.

**Indexing and aggregation.** Library entries are indexed by the pair
$(p,k)$; all statistics ($Q,A^{\rm tgt},\pi,M$, etc.) are computed at the
entry level. Source-level aggregation is only used for reporting and final
narrative:

$$Q_p^{\rm src}=\max_kQ_{p,k},\qquad k_p^*=\arg\max_kQ_{p,k}.$$

**Selection rule.** Different transforms of the same source are naturally
highly correlated, so:

> **Rule 4.1**: all transform entries of the same source feature are
> **forced into the same family** during the raw-redundancy clustering of
> §36, unless they fall into different mechanism groups (§35) *and* each
> independently passes eligibility (§34). In other words, a source feature
> contributes at most one selected entry in principle; only when $t_1$ and
> $t_2$ clearly predict **different mechanisms** (e.g. $t_1$ predicts
> direction, $t_2$ predicts activity) can both be selected simultaneously.

This rule keeps transform-family cost under control: the library's nominal
size is multiplied by $|\mathcal T|$, but the final selected count doesn't
inflate.

### §4.4 Lemma 4.1 (the absorption property of daily demeaning)

Since $X^{XS}_d$ has zero mean each period, for any matrix $\Phi_d$,

$$\frac1{N_d}X_d^{XS\top}\Phi_d=\frac1{N_d}X_d^{XS\top}\big(\Phi_d-\mathbf 1_{N_d}\bar\Phi_d\big),$$

where $\bar\Phi_d$ is $\Phi_d$'s that-period row average.

**Proof**: $X_d^{XS\top}\mathbf 1_{N_d}=0$, so $X_d^{XS\top}\mathbf
1\bar\Phi_d=0$. $\blacksquare$

**Corollary 4.2**: $C^{(d)}$ is automatically a **that-period
cross-sectional covariance**, not an uncentered cross-moment; and any
"that-period-constant" component of $\Phi$ contributes exactly zero to $C$.

This lemma gets used in §8, §9, §13, and §29. It explains two things at
once: why the $-I_K$ term in §8.2 doesn't affect any result, and why taking
the ticker intersection under permutation **requires re-centering**
(because on the intersection $X^\top\mathbf 1\ne0$, breaking the lemma's
precondition).

---

## §5. Label Preprocessing

### Definition

For each horizon:

$$\boxed{Y^{XS}_{id,h}=Y^{raw}_{id,h}-\frac1{N_d}\sum_{j\in\mathcal U_d}Y^{raw}_{jd,h}.}$$

So $\sum_iY^{XS}_{id,h}=0$.

### Two design choices and their consequences

**Choice one: demean, don't beta-neutralize.** Demeaning removes that
period's equal-weighted market return. This fixes the entire analysis's
meaning to "**cross-sectional contrast relative to the equal-weighted
market**." This means:

- A pure market-timing signal (predicting overall market direction) gets
  $Q=0$ under this method. That's correct — it isn't a cross-sectional
  signal.
- If the universe's sector composition is systematically skewed, demeaning
  does not remove sector-level common movement. To analyze purely
  idiosyncratic predictive power, run a cross-sectional regression against
  known factor exposures after this step and take the residual. This is an
  optional additional preprocessing layer that doesn't change any formula
  downstream.

**Choice two: no daily volatility normalization.** I.e. we do not set
$Y_d\leftarrow Y_d/\hat\sigma_d$.

The consequence is that $\Sigma_Y$ (§7) is the full-sample covariance, and
$\|z_{id}\|$ is systematically larger on high-volatility days. As discussed
in §3, this makes $Q$ dispersion-weighted. This choice is deliberate:
**the cross-sectional dispersion of future returns is itself the
information the $\mathsf V$ and $\mathsf L$ blocks are meant to analyze.**
Normalizing it away first would leave the $\mathsf L$ block (directionality
during large moves) with nothing to be defined over. The cost is that $Q$
is not regime-immune, and must be read alongside $T_{\rm eff}$
(Definition 3.1).

---

## §6. Horizon Bands

### Motivation

Prediction problems at different horizons are economically different
problems. 5-minute returns are mainly governed by microstructure and
liquidity supply; 24-hour returns start to reflect information diffusion
and overnight risk premia. Mixing them into one whitening and one SVD
produces a coordinate system dominated by the highest-variance horizon,
meaningless for the others. Hence: band first, then analyze independently
within each band.

### Definition

Partition the $H$ horizons into bands $H_1,\dots,H_{B_{\rm band}}$. For
band $\beta$:

$$K=|H_\beta|,\qquad Y_d^{(\beta)}=Y_d^{XS}[:,H_\beta].$$

Everything downstream runs independently within each band.

### Overlapping horizons vs. an incremental basis

This is a specification point that must be chosen explicitly. If
$Y_{\cdot,h}$ are **cumulative** returns measured from the same starting
point, different horizons are naturally nested: under a random-walk
assumption,

$$\operatorname{corr}(r_{0\to1h},r_{0\to2h})=\sqrt{1/2}\approx0.71,\qquad\operatorname{corr}(r_{0\to1h},r_{0\to4h})=\sqrt{1/4}=0.5.$$

A $K=4$ cumulative band has a correlation matrix with condition number
roughly 10–20, and whitening amplifies the smallest-eigenvalue direction by
$\lambda_{\min}^{-1/2}$ — amplifying estimation noise right along with it.

> **Specification 6.1**: within a band, prefer **non-overlapping incremental
> returns**
> $$r_{0\to1h},\quad r_{1h\to2h},\quad r_{2h\to4h},\quad r_{4h\to8h},$$
> whose correlation matrix under a random walk is the identity, with
> condition number near 1. If the research question genuinely requires
> cumulative returns (e.g. execution holds to a $4h$ close), keep the
> cumulative form, but §7's rank control becomes mandatory rather than
> optional.
>
> Either way, the band definition (cumulative or incremental) must be
> recorded in configuration and appear in the report, because $c_p$'s
> coordinate meaning depends entirely on it.

### A sense of dimension

This book uses the following band split as its running example:

$$H_1=\{5m,10m,15m,30m\},\quad H_2=\{1h,2h,4h\},\quad H_3=\{6h,12h,24h\},$$

| band | $K$ | $d_{\mathsf D}$ | $d_{\mathsf V}=K(K{+}1)/2$ | $d_{\mathsf L}$ | $d_\Phi$ |
|---|---:|---:|---:|---:|---:|
| $H_1$ | 4 | 4 | 10 | 4 | 18 |
| $H_2$ | 3 | 3 | 6 | 3 | 12 |
| $H_3$ | 3 | 3 | 6 | 3 | 12 |

Remember this table. $d_{\mathsf V}$ grows quadratically in $K$ while
$d_{\mathsf D}$ only grows linearly: at $K=4$, the Activity block already
takes up $10/18=56\%$ of all coordinates. This fact has decisive
consequences in two later places: §19 (the fingerprint's dimensional
baseline) and §38 (coverage's rank ceiling).
