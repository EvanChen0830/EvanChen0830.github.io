---
title: "Factor Definition and Evaluation"
description: "How the factor set used across the pipeline is defined and judged: centered point-in-time exposures, characteristic factors and their returns, PCA factors and why they rotate, what stays stable under rotation, how to align and evaluate a factor set, and what each later chapter takes from it."
date: 2026-11-08
tags: ["factor", "risk-model", "pca", "methodology", "quant-research"]
math: true
aliases:
  - "/quant/factor/factor-research-framework/"
---

## The problem this chapter solves

Universe, Feature and Alpha research each ask a factor question:

- Universe: is this candidate sub-universe just a factor tilt?
  ([Universe §21](/quant/universe/multi-sub-universe-directional-coupling-geometry/#21-beta-demeaning-and-the-factor-neutral-variant))
- Feature: how much of this feature's predictive power is a known factor, and
  **which** factor does it predict?
  ([Predictive Geometry §35](/quant/feature-evaluation/cross-sectional-predictive-geometry/#35-factor-decomposition-of-a-feature))
- Alpha: on **which** factor does the model's output have predictive power, and
  does that survive the portfolio's neutrality rules?
  ([Label–Alpha §28](/quant/alpha-evaluation/label-alpha-geometry/#28-factor-aligned-and-specific-alpha))

All three answers are only as good as the factor set they are measured
against. This chapter fixes that object: how exposures are defined, how factor
returns are obtained, how statistical (PCA) factors are made usable despite
their rotation, and how a factor set is evaluated before any other chapter
relies on it.

> **Which factor set, built how, is stable and accurate enough that "this
> alpha predicts factor $k$" means the same thing next month?**

---

## 1. Notation

Symbols from
[Cross-Sectional Predictive Geometry §1](/quant/feature-evaluation/cross-sectional-predictive-geometry/#1-notation)
keep their meaning: $d$ is the decision period, $i$ the ticker, $h$ the
horizon, $W=\operatorname{diag}(1/(TN_d))$, $\Phi=[\Phi_{\mathsf D}\,|\,\Phi_{\mathsf V}\,|\,\Phi_{\mathsf T}]$,
$c_s=s^\top W\Phi$. Factor objects use letters that are not already taken:

| Symbol | Shape | Meaning |
|---|---|---|
| $m$ | scalar | number of factors |
| $F_d$ | $N_d\times m$ | exposures known at $d$, **within-date centered** (Section 3) |
| $\bar F_d$ | $N_d\times(m+1)$ | $[F_d\ \mathbf1]$, exposures plus the dollar-neutral direction |
| $y_{d,h}$, $Y_d$ | $N_d$, $N_d\times K$ | within-date demeaned future returns (one horizon / one band) |
| $\hat\psi_{d,h}$, $\hat\Psi_d$ | $m$, $m\times K$ | factor returns estimated from the cross-section at $d$ |
| $\hat u_{d,h}$ | $N_d$ | specific (residual) returns |
| $\Sigma^{\psi}_{d}$, $\Sigma^{u}_{d}$ | $m\times m$, $N_d\times N_d$ | factor-return covariance, specific covariance |
| $P_{F_d}$, $M_{F_d}$ | $N_d\times N_d$ | projector onto $\operatorname{span}(F_d)$ and its complement |
| $\theta_j$ | angle | principal angles between two factor subspaces (Section 5.3) |

$B$, $f$, $D$, $K$ and $w$ are deliberately not used for factor objects: in
this framework they already mean the whitening matrix or bootstrap count, the
mechanism fingerprint, the Direction block, the label dimension and the row
weights.

---

## 2. Two kinds of factor

| | Characteristic factors | Statistical (PCA) factors |
|---|---|---|
| Exposure is | a measured attribute: beta, size, sector, momentum, funding rate, volatility | a loading estimated from the joint behaviour of returns |
| Defined before seeing returns? | yes | no: the exposures **are** an estimate |
| Interpretable column by column? | yes | only after alignment (Section 5.4) |
| Stable column by column? | usually, if the attribute is persistent | no: signs flip, columns swap and rotate (Section 5.2) |
| What it is good for | naming what an alpha predicts; risk budgets by economic driver | capturing common movement that no named attribute describes; detecting omitted factors |

The two are complements, not alternatives. The default in this framework is a
**characteristic set** for attribution and a **PCA set on the characteristic
residuals** (Section 5.6) to catch what the named factors miss.

---

## 3. The exposure contract

Every factor set, of either kind, is delivered to the other chapters as
$F_d$ satisfying four rules.

**Point-in-time.** Every input to $F_d$ must have been observable at $d$. A
sector classification revised later, or a PCA fit on a window that reaches
past $d$, is look-ahead.

**Within-date centered.**

$$
\boxed{
F_d\leftarrow\Big(I-\tfrac1{N_d}\mathbf1\mathbf1^\top\Big)F_d^{raw},
\qquad
\mathbf1^\top F_d=0.
}
$$

Every label in the framework is already within-date demeaned, so the market's
common move has been removed from $Y$. Centering the exposures keeps the two
consistent:

- The market factor survives only as **beta dispersion**,
  $\beta_i-\bar\beta$. That is exactly the part demeaning cannot remove
  ([Universe §21](/quant/universe/multi-sub-universe-directional-coupling-geometry/#21-beta-demeaning-and-the-factor-neutral-variant)).
  A pure intercept column becomes zero after centering and is dropped.
- With centered exposures, residualizing a demeaned vector against $F_d$
  leaves it demeaned, under any weighting. Each chapter relies on this.
- Nothing is lost for **dollar-neutral** portfolios: their exposures, returns
  and risk are identical whether computed with raw or centered exposures. A
  portfolio that is not dollar-neutral needs the market intercept and a
  separate market model.

**Analysis cross-section.** Wherever a signal is defined only on some names
(missing values, matched thinning in a null), restrict $F_d$ to those names
and re-center it there. Otherwise "centered" no longer holds on the set that
is actually used.

**Rank is logged.** Non-overlapping sector dummies sum to the constant vector,
so after centering they lose one rank. This is harmless for every projection
($P_{F_d}$ depends only on the span of $F_d$) and matters only when sector
returns are reported (Section 4.2).

---

## 4. Characteristic factors

### 4.1 Construction

For each characteristic, each date:

1. winsorize cross-sectionally (default 1% / 99%, as for features);
2. standardize within date;
3. center (already true after step 2 for continuous attributes; required
   explicitly for dummies).

Group memberships enter as dummies $H_d$ (sector, industry, asset class, or a
crypto sector taxonomy), centered the same way.

### 4.2 Factor returns

For each date and horizon, regress the demeaned returns on the exposures
**plus an intercept column**, and keep the factor coefficients:

$$
\begin{bmatrix}\hat\psi_{d,h}\\ \hat a_{d,h}\end{bmatrix}
=\big(\bar F_d^\top\Pi_d\bar F_d\big)^{+}\bar F_d^\top\Pi_d\,y_{d,h},
\qquad
\hat u_{d,h}=y_{d,h}-F_d\hat\psi_{d,h}-\hat a_{d,h}\mathbf1 .
$$

$\Pi_d$ is the regression weight: $I$ by default, or the inverse specific
variance when a risk model exists.

**Why the intercept column, when exposures are centered.** With $\Pi=I$ it
changes nothing: centered exposures with no intercept give the same
$\hat\psi$. With any other weight it matters, because a weighted regression
removes the *weighted* mean, not the equal mean. Keeping $\mathbf1$ in the
basis makes $\hat\psi$ correct under every weight, and the same whether the
labels were demeaned or not.

**Reading $\hat\psi$.** Row $k$ of $(\bar F^\top\Pi\bar F)^{-1}\bar F^\top\Pi$
is a portfolio with unit exposure to factor $k$, zero exposure to every other
factor, and zero net weight. $\hat\psi_{d,h,k}$ is that **pure factor
portfolio's** return.

**Sector returns need a stated restriction.** Because centered dummies lose
one rank, individual sector returns are defined only relative to a choice:
drop one sector, or require the sector returns to sum to zero. A
pseudo-inverse silently picks the equal-weight sum-to-zero version. The
choice changes how each sector's return reads, never the fitted values, never
any projection.

---

## 5. Statistical (PCA) factors

### 5.1 Definition

On each date $d$, take a trailing window of demeaned returns that **ends
before $d$**, compute the top $m$ principal components of the
cross-sectional covariance, and use the loadings as exposures:

$$
\widehat{\operatorname{Cov}}_{[d-L,\,d-1]}=V\Lambda V^\top,\qquad
F^{\rm pca}_d=\text{center}\big(V_{:,1:m}\big).
$$

Because the returns are within-date demeaned, $\mathbf1$ lies in the null
space of their cross-sectional covariance, so the retained loadings are
already centered and orthonormal; the centering step only matters if a
different return definition is used. Names that enter the universe without
enough history get no loading; treat them as missing (Section 3, analysis
cross-section).

### 5.2 Why PCA factors rotate

An eigenvector is defined only up to the eigenvalue it belongs to. Three
things follow, and all three happen in real data:

| Instability | Cause | Symptom |
|---|---|---|
| **Sign flip** | $v$ and $-v$ are equally valid | factor $k$'s exposure and return change sign between windows |
| **Swap** | two eigenvalues cross as the window rolls | "factor 2" this month is "factor 3" last month |
| **Rotation** | two or more eigenvalues are close (small eigengap) | the columns mix: each window returns a different orthogonal basis of the **same** plane |

Rotation is the dangerous one, because nothing flags it: each window's PCA is
correct, and yet "factor 2" means a different combination of names every
month. Any number reported **per factor** — an exposure, a factor return, an
alpha's predictive power on factor 2 — inherits that instability.

### 5.3 What does not rotate: the subspace

Rotation, sign flips and swaps among the retained components all leave
$\operatorname{span}(F_d)$ unchanged, as long as the eigengap at the cut
($\lambda_m-\lambda_{m+1}$) is clear. Stability is therefore measured between
**subspaces**, with principal angles: for two orthonormal loading bases
$V_a$, $V_b$ ($N\times m$),

$$
\cos\theta_j=\sigma_j\big(V_a^\top V_b\big),\qquad
\text{subspace similarity}=\frac1m\sum_j\cos^2\theta_j\in[0,1].
$$

This is the same principal-angle similarity the Universe chapter uses for
canonical-subspace stability.

Every quantity in the other chapters falls into one of two groups:

| Depends only on the subspace (safe for PCA factors) | Depends on the basis (needs Section 5.4) |
|---|---|
| residualization $M_{F_d}$; factor-neutral labels (Universe §21.2) | individual exposures $F_{d,k}$ |
| a signal's factor share $\rho^2_F$ and its factor / specific split (Predictive §35, Label–Alpha §28) | individual factor returns $\hat\psi_{d,k}$ |
| membership factor share $\rho^2_F(g)$ (Universe §21.3) | "which factor" attribution, factor by factor (Predictive §35.7, Label–Alpha §28.6) |
| the **total** factor part of a signal's predictive power | exposure strata for the exposure-matched null (Universe §21.4) |
| portfolio neutrality to the whole factor set | budgets or limits on one factor |

The left column can use raw PCA factors directly. The right column cannot.

### 5.4 Making PCA factors usable column by column

Apply these in order; each one removes one source of instability.

1. **Eigengap blocks.** Group consecutive components whose eigenvalues are
   within a pre-set ratio (default: $\lambda_{j+1}/\lambda_j>0.9$) into one
   block. A block is reported as one unit; only blocks of size 1 are reported
   as single factors.
2. **Procrustes alignment.** Rotate each window's loadings to best match a
   reference basis $V_{\rm ref}$ (the previous window, or a fixed reference
   window):
   $$
   Q=\arg\min_{Q^\top Q=I}\|V_dQ-V_{\rm ref}\|_F=U\,\tilde V^\top,
   \qquad U\Sigma\tilde V^\top=\operatorname{svd}\big(V_d^\top V_{\rm ref}\big).
   $$
   Restrict $Q$ to be block-diagonal over the eigengap blocks if the blocks
   should stay separate. Alignment changes the columns, never the subspace.
3. **Sign convention.** After alignment, fix each column's sign (e.g. positive
   loading-weighted average of a reference attribute such as beta).
4. **Target rotation for names.** To read PCA factors in economic terms,
   rotate the subspace toward the characteristic exposures (Procrustes with
   $V_{\rm ref}$ = normalized characteristic exposures). The result is
   "the part of the PCA subspace that looks like beta", etc. — interpretable,
   but only as good as the match, which must be reported.

Chained alignment (each window to the previous) can drift slowly; alignment to
a fixed reference cannot drift but degrades if the reference period is
unrepresentative. Report which was used.

### 5.5 How many factors

Choose $m$ by comparing eigenvalues to a null, not by a fixed count:

- **Permutation parallel analysis** — the same "$j$-th against $j$-th"
  comparison the Predictive Geometry chapter uses for $r_{\rm sig}$: shuffle
  each name's return series in time independently (block-wise, to keep
  autocorrelation), recompute the spectrum, keep components whose eigenvalue
  exceeds the null's 95th percentile at the same rank.
- **Marchenko–Pastur edge** as a quick reference when returns are close to
  independent across names.

Then check the eigengap at the cut. A cut inside a near-degenerate block makes
$m$ itself unstable: the $m$-th component rotates in and out of the retained
set. Move the cut to a block boundary.

### 5.6 PCA on characteristic residuals

Running PCA on $\hat u_{d,h}$ (returns after the characteristic factors)
instead of raw returns gives **residual statistical factors**: common
movement that the named factors do not explain. Two uses:

- appended to the characteristic set, they complete the risk model;
- each retained residual component is a candidate **omitted factor**, to be
  named by target rotation (Section 5.4, step 4) or by correlating its
  loadings with candidate attributes.

---

## 6. Evaluating a factor set

The same five checks apply to characteristic, PCA and mixed sets. All are
computed out of sample: exposures from $d$, returns realized after $d$.

### 6.1 Explanatory power

Out-of-sample cross-sectional $R^2$ of returns on exposures, averaged with
equal date weights (as everywhere in this framework):

$$
R^2_{F}=\frac1T\sum_d\frac{\|P_{F_d}y_{d,h}\|^2}{\|y_{d,h}\|^2}.
$$

More factors always raise in-sample $R^2$; only the out-of-sample version is
informative.

### 6.2 Risk forecast accuracy

For a set of test portfolios (pure factor portfolios, random dollar-neutral
portfolios, and the strategies actually traded), standardize each realized
return by its predicted volatility and look at the dispersion:

$$
z_{d}=\frac{\omega_d^\top y_{d,h}}{\sqrt{\omega_d^\top\big(F_d\Sigma^{\psi}_dF_d^\top+\Sigma^{u}_d\big)\omega_d}},
\qquad
b=\operatorname{sd}_d(z_d).
$$

$b\approx1$ means risk is forecast correctly; $b>1$ means it is understated.
A good $R^2$ with $b$ far from 1 means the factors explain returns but the
covariance model on top of them is wrong.

### 6.3 Residual dependence

Compute the eigenvalues of the residual correlation matrix of $\hat u$ and
compare them with the same null as Section 5.5. Eigenvalues above the null
mean common structure remains in the residuals: an omitted factor (Section
5.6). Test cross-sectional residual dependence and time-series residual
autocorrelation separately; they are different failures.

### 6.4 Stability

| Quantity | Definition | Reading |
|---|---|---|
| exposure persistence | cross-sectional correlation of $F_{d,k}$ with $F_{d-1,k}$, averaged | low means turnover in what "factor $k$" means; check the attribute or the alignment |
| subspace stability | similarity (Section 5.3) between consecutive windows, and across block-bootstrap resamples | the primary stability measure for PCA sets |
| column stability after alignment | per-column $\cos$ between aligned consecutive windows | must be high for every factor that will be attributed individually |
| eigengap at the cut | $\lambda_m/\lambda_{m+1}$ | close to 1 means $m$ is not well defined |

### 6.5 Identifiability

Report $\operatorname{rank}F_d$, the condition number of $F_d^\top F_d$, and
the restriction used for group dummies. A poorly conditioned set gives
unstable individual factor returns even when the subspace is fine.

### 6.6 Factor-set report

| Field | Content |
|---|---|
| `factor_basis_id` | version, list of factors, kind (characteristic / PCA / residual PCA) |
| `m`, `blocks` | number of factors; eigengap blocks for PCA |
| `align` | none / chained Procrustes / fixed reference / target rotation |
| `R2_oos` | Section 6.1 |
| `bias_stat_median`, `bias_stat_p90` | Section 6.2 across test portfolios |
| `resid_eig_above_null` | Section 6.3 |
| `subspace_sim_med`, `subspace_sim_p10` | Section 6.4 |
| `col_stability_min` | worst aligned column |
| `cond_FtF_med`, `rank_min`, `dummy_restriction` | Section 6.5 |
| `usable_per_factor` | the list of factors that pass the column-stability gate and may be used for per-factor attribution |

**Gates.** A set is usable for subspace quantities when its subspace
similarity and explanatory power pass. It is usable for **per-factor**
quantities only for the factors listed in `usable_per_factor`; every other
factor is reported at block level or not at all.

---

## 7. What each chapter takes from the factor set

| Chapter | Uses | Needs per-factor stability? |
|---|---|---|
| [Universe §21](/quant/universe/multi-sub-universe-directional-coupling-geometry/#21-beta-demeaning-and-the-factor-neutral-variant) | factor-neutral label variant, membership factor share | no |
| | exposure-matched membership null | yes (strata are built per factor) |
| [Predictive Geometry §35](/quant/feature-evaluation/cross-sectional-predictive-geometry/#35-factor-decomposition-of-a-feature) | factor share and factor / specific split of each feature | no |
| | which factor a feature predicts (§35.7) | yes |
| [Selection Bridge](/quant/feature-evaluation/selection-bridge-eligibility-grouping-coverage/#factor-aware-eligibility) | mandate-dependent eligibility, redundancy, Score | no |
| [Label–Alpha §28](/quant/alpha-evaluation/label-alpha-geometry/#28-factor-aligned-and-specific-alpha) | factor / specific split, four-cell decomposition, decision alpha | no |
| | which factor an alpha predicts (§28.6) | yes |

---

## 8. Execution order

1. Choose the characteristic attributes and group taxonomy; lock definitions.
2. Build $F_d$ point-in-time; winsorize, standardize, center; log rank.
3. Estimate $\hat\psi$, $\hat u$ by the augmented regression (Section 4.2).
4. On $\hat u$, run rolling PCA; choose $m$ by parallel analysis; form eigengap
   blocks; align (Section 5.4).
5. Evaluate the combined set (Section 6) on a period not used to design it.
6. Publish `factor_basis_id`, the report, and `usable_per_factor`. Every other
   chapter cites this ID with every factor-related number.

## 9. Fit / transform boundary

| Object | Fit on | Applied to |
|---|---|---|
| characteristic definitions, winsorization rule | design period | all dates |
| $F_d$ (characteristic) | data observable at $d$ | date $d$ |
| PCA loadings | trailing window ending at $d-1$ | date $d$ |
| $m$, eigengap ratio, alignment reference | design period | all dates |
| $\hat\psi_{d,h}$, $\hat u_{d,h}$ | the cross-section at $d$ and returns over $(d,d+h]$ | attribution only, never as an input at $d$ |
| $\Sigma^{\psi}_d$, $\Sigma^{u}_d$ | data strictly before $d$ | date $d$ |

## 10. Numerical checks

`scripts/verify_factor_math.py` checks every identity used by this chapter
and by the factor sections of the other chapters: centering is lossless for
dollar-neutral portfolios, residuals stay centered, the augmented regression
is invariant to demeaning while the centered one is not under weighting,
projectors depend only on the span, Procrustes alignment leaves the projector
unchanged, and the attribution identities of Predictive Geometry §35 and
Label–Alpha §28 hold exactly.

## 11. Open research questions

| ID | Question |
|---|---|
| IR-23 | Which characteristic set and taxonomy work across equities, futures, FX and crypto, and which role does each factor serve (explain, earn premium, constrain)? |
| IR-24 | How many factors are estimable, and does a mixed characteristic + residual-PCA set beat either alone on Section 6? |
| IR-25 | Characteristic exposures versus estimated betas versus state-space betas: which are most stable and best for attribution? |
| IR-26 | Are factor premia and factor timing predictable, and by which features? (the attribution of Predictive Geometry §35.7 supplies candidates) |
| IR-27 | Is specific risk diagonal, block, sparse or heavy-tailed; which residual PCA components are real omitted factors? |
| IR-28 | Which factors should be free, budgeted or hard-neutral, given what the alpha attribution (Label–Alpha §28.6) shows? |
| IR-29 | How sensitive are attributions and portfolios to the choice among legitimate factor sets? |

## 12. What this chapter can and cannot answer

**Can answer**

- Whether a factor set is accurate, stable and identifiable enough to be used.
- Which factor-related numbers are trustworthy for PCA factors (subspace
  quantities) and which need alignment first (per-factor quantities).
- How factor returns are obtained consistently with demeaned labels.

**Cannot answer**

- Which factor set is "true". Every attribution in the other chapters is
  relative to the published `factor_basis_id`.
- Whether a factor is priced. That is a question about $\hat\psi$ over time
  (IR-26), not about exposures.
