---
title: "Cross-Sectional Predictive Geometry, Part 4: The Predictive Operator, SVD, and Permutation Nulls"
description: "C = X^T W Phi turns out to be literally a matrix of ICs, bounded by at most d_Phi detectable modes regardless of library size. Calibrating its null distribution correctly means a circular date-shift permutation with matched thinning for universe mismatch, plus a stratified alternative for slow-moving features that would otherwise defeat the test entirely."
date: 2026-09-14
series: ["Cross-Sectional Predictive Geometry"]
part: 4
tags: ["cross-sectional", "methodology", "quant-research", "predictive-geometry", "permutation-test"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Predictive Geometry" >}}

# Part Three: Feature–Label Geometry

## §11. The Predictive Operator

### Definition

For each period,

$$\boxed{C^{(d)}=\frac1{N_d}X_d^{XS\top}\Phi_d\ \in\mathbb R^{P\times d_\Phi},}$$

over the full period (each feature averaged over its own valid-day set
$\mathcal D_p$, per §4.2),

$$\boxed{C_{p,:}=\frac1{T_p}\sum_{d\in\mathcal D_p}C^{(d)}_{p,:}.}$$

When every feature has 100% coverage, this reduces to
$C=\frac1T\sum_dC^{(d)}=X^\top W\Phi$.

$C$ is the **average daily cross-sectional predictive operator**. Row $p$,
$C_{p,:}$, is feature $p$'s average cross-sectional association with every
future-return mechanism coordinate.

### Scale: Every Element of $C$ Is Literally an IC

By §4, $X^{XS}_{\cdot p}$ has zero mean and unit variance each period; by
§9, each column of $\Phi$ has unit norm under the weighted inner product.
So

$$C_{p,k}\approx\overline{\operatorname{corr}_d\big(X_{\cdot dp},\Phi_{\cdot dk}\big)}$$

— feature $p$'s **average per-period cross-sectional correlation** with the
$k$-th mechanism coordinate, i.e. literally an IC. Hence

$$Q_p=\|c_p\|^2\approx\sum_k\overline{IC}_{p,k}^2.$$

This gives the scale for reading $Q$:

| Average per-period IC | $Q$ contribution of one coordinate | Meaning in equity markets |
|---:|---:|---|
| 0.005 | $2.5\times10^{-5}$ | almost undetectable |
| 0.01 | $1\times10^{-4}$ | weak, but tradeable at high frequency |
| 0.03 | $9\times10^{-4}$ | a good single feature |
| 0.05 | $2.5\times10^{-3}$ | very strong |
| 0.10 | $1\times10^{-2}$ | extremely rare; suspect data leakage first |
| 0.20 | $4\times10^{-2}$ | almost certainly look-ahead bias |

**Every numeric example in this book falls in the band $Q^{\mathsf
D}\in[10^{-5},10^{-3}]$.** Getting this scale wrong by one or two orders of
magnitude will mis-set every threshold downstream, deleting all real
features or keeping all noise.

### Scale of Computation

With $P=400$, $d_\Phi=18$, $C$ is a $400\times18$ matrix — about 7200
numbers. The entire raw material of a feature study is these 7200 numbers
plus their null distribution.

---

## §12. Singular Value Decomposition

### Definition

$$\boxed{C=USV^\top,}\qquad S=\operatorname{diag}(s_1,\dots),\ s_1\ge s_2\ge\cdots$$

$V$'s columns are label-mechanism directions; $U$'s columns are feature-side
combinations.

### Interpretation

$C$'s SVD answers: **across the $d_\Phi$-dimensional label-mechanism space,
which directions does the entire feature library collectively
"illuminate"?**

- $v_1$: the whole library's most-agreed-upon predictive direction. On most
  equity data this will be a direction dominated by the $\mathsf D$ block,
  close to "all horizons point the same way" — the classic
  momentum/reversal axis.
- $u_1$: the feature-side linear combination achieving this prediction — can
  be thought of as an automatically-extracted super-feature.
- $s_1$: that direction's strength.

**Note the asymmetry between $U$ and $V$.** $V$ lives in label space
(dimension $d_\Phi$, small); $U$ lives in feature space (dimension $P$,
large). Therefore

$$r\le\min(P,d_\Phi)=d_\Phi.$$

**The ceiling on the number of detectable modes is set by the label space's
dimension, not the number of features.** At $K=4$, $r\le18$; in practice
$r_{\rm sig}$ almost always lands between 2 and 6. This fact has decisive
consequences in §38.

---

## §13. The Permutation Null

### Motivation

$Q_p$ needs a reference distribution. It cannot be given in closed form,
because $C$'s null distribution depends simultaneously on the feature's
cross-sectional distribution, $\Phi$'s cross-sectional distribution, and
both of their own temporal autocorrelation structure. The only reliable
approach is a resampling scheme that preserves all of that structure and
destroys only "which day pairs with which."

### §13.1 The Null Hypothesis

The original pairing is $X_d\leftrightarrow\Phi_d$. The null hypothesis is:

> the feature's cross-sectional ranking has no predictive relationship with
> the future-return mechanism **at the correct time**.

So we break date alignment: $X_{\pi(d)}\leftrightarrow\Phi_d$, while
preserving

1. each period's full cross-sectional structure in $X$ (including sector,
   market-cap clustering, etc.);
2. each period's full cross-sectional structure in $\Phi$;
3. ticker identity;
4. $X$'s own temporal autocorrelation;
5. the label's own temporal autocorrelation.

Preserving these five is what keeps the null from being artificially small
just because cross-sectional structure got scrambled — that's the most
common error, and it makes everything look significant.

### §13.2 Primary Scheme: Circular Date Shift

Given ordered dates, pick a shift $s$, and set

$$\pi_s(k)=1+\big((k+s-1)\bmod T\big).$$

The set of valid shifts is

$$\mathcal S=\{s:\Delta\le\min(s,T-s)\},\qquad\boxed{\Delta=H_{\max}+\ell_X,}$$

where $H_{\max}$ is the longest horizon in the band (in periods), and
$\ell_X$ is the feature's memory length (e.g. $\ell_X=20$ for a 20-day
moving average). This minimum offset ensures the paired $X$ and $Y$ never
share an overlapping time window. **$\Delta$ and the purge gap in §41 must
be the same number** — otherwise the same data is held to different leakage
standards in two places.

In practice, use $B_{\rm null}$ of the $|\mathcal S|$ possible shifts
(default: take all, capped at 2000).

### §13.3 Universe Mismatch: the Matched-Thinning Rule

$N_d$ changes over time, so $\mathcal U_{\pi_s(d)}\ne\mathcal U_d$. Three
things must be handled together.

**(a) Take the intersection.** Let

$$\mathcal I_{d,s}=\mathcal U_{\pi_s(d)}\cap\mathcal U_d.$$

If $|\mathcal I_{d,s}|<\rho_{\min}\cdot\min(N_{\pi_s(d)},N_d)$ (default
$\rho_{\min}=0.9$), discard the $(d,s)$ pair.

**(b) Re-standardize on the intersection.** By Lemma 4.1, $C^{(d)}$ is a
genuine covariance *provided* $X$ has zero mean on the universe used.
Restricting to $\mathcal I_{d,s}$ breaks this precondition: if a feature
correlates with "being present on both dates" (market cap, liquidity,
listing age naturally do), the intersection's mean is nonzero, and the
null's $C$ gets systematically contaminated by a nonzero mean term,
inflating or deflating the null for that class of feature. So **$X$ must be
re-demeaned and re-standardized on $\mathcal I_{d,s}$**:

$$\check X_{idp}=\frac{X'_{idp}-\mu_{dp}^{\mathcal I}}{\sigma_{dp}^{\mathcal I}+\epsilon},\qquad i\in\mathcal I_{d,s}.$$

$\Phi$ needs no reprocessing: by Lemma 4.1, as long as $\check X$ has zero
mean on $\mathcal I$, any that-period-constant component of $\Phi$ is
automatically absorbed.

**(c) Apply the same sampling rule to the observed statistic.** If the
observed statistic is computed on the full $\mathcal U_d$ while the null is
computed on a smaller intersection, the two have different effective sample
sizes, and the null's variance will be over- or under-estimated depending
on whether the intersection is selective. The symmetric rule is:

> **Specification 13.1 (matched thinning)**: for each date $d$, draw $s_d
> \sim\text{Uniform}(\mathcal S)$ independently, define the observation
> analysis universe
> $$\mathcal A_d=\mathcal U_d\cap\mathcal U_{\pi_{s_d}(d)},$$
> and re-standardize on $\mathcal A_d$ per (b) before computing the observed
> statistic. Repeat $R_{\rm thin}=5$ independent draws and average, as the
> final observed value.

This way, the observed and null statistics have identically distributed
universe sizes.

**Simplifying condition**: if $\operatorname{median}_{d,s}|\mathcal
I_{d,s}|/N_d\ge0.98$, thinning's effect is smaller than numerical error,
and observed values can be computed directly on $\mathcal U_d$. This ratio
must be computed and recorded, never assumed.

### §13.4 Slow-Moving Features and an Alternative Null

The circular-shift test's power depends on whether, after shifting by
$\Delta$ periods, a feature's cross-sectional ranking has actually changed.
For a feature with memory much longer than $\Delta$ (e.g. 250-day momentum,
market cap, sector dummies),

$$X_{d+\Delta,p}\approx X_{d,p},$$

so $C^{\rm null}\approx C^{\rm obs}$, and the test has almost no power.
**This systematically biases the method toward favoring fast-rotating
signals** — if left unhandled, slow but genuine signals get judged
not-significant across the board.

**Persistence diagnostic.** For each feature, define

$$\boxed{\kappa_p=\operatorname{median}_d\ \operatorname{corr}_{i\in\mathcal U_d\cap\mathcal U_{d+\Delta}}\big(X_{idp},X_{i,d+\Delta,p}\big).}$$

$\kappa_p$ must appear in the report. Rule of thumb: the circular-shift
null's power has dropped noticeably once $\kappa_p>0.5$, and is near zero
by $\kappa_p>0.8$.

**Alternative null (N2): stratified within-date permutation.** For features
with $\kappa_p>0.5$, replace the date shift with

$$X_{idp}\longrightarrow X_{\pi_d(i),dp},$$

where $\pi_d$ is a random permutation **within strata**:

$$\mathcal U_d=\bigsqcup_g\mathcal U_{d,g},\qquad g=\text{sector}\times\text{size bucket},$$

permuted independently within each $g$.

**Stratification is necessary, not a refinement.** Permuting freely across
the whole universe would destroy both sector and market-cap clustering at
once, making the null for a sector-level feature far too small, and the
test severely anti-conservative. Stratification fixes these known
exposures' cross-sectional composition in place, so the null only
destroys "individual-stock-level pairing."

A cleaner variant residualizes known exposures first: let $F$ be the known
factor-exposure matrix (sector, size, beta),

$$X^\perp_{d,p}=(I-P_{F_d})X_{d,p},\qquad P_{F_d}=F_d(F_d^\top F_d)^{-1}F_d^\top,$$

then run free within-date permutation on $X^\perp$. The null hypothesis
becomes "after controlling for known exposures, individual-stock-level
feature ranking has no predictive power" — more precise semantics, at the
cost that the conclusion covers only the residual portion.

**Reporting rule**: for features with $\kappa_p>0.5$, the primary conclusion
uses N2, with N1's result and $\kappa_p$ also listed; for $\kappa_p\le0.5$,
the primary conclusion uses N1. Which null was used must be marked per
feature in the report.

### §13.5 Implementation Flow

```text
Input: {X_d}, {Phi_d}, valid shift set S, number of resamples B_null
Precompute: kappa_p (persistence diagnostic) for each feature

# ---- observed statistic (matched thinning) ----
for r in 1..R_thin:
    for d in 1..T:
        s_d ~ Uniform(S)
        A_d = U_d ∩ U_{pi_{s_d}(d)}
        if |A_d| < rho_min * min(N_d, N_{pi_{s_d}(d)}): skip this d
        X_check = X_d re-standardized on A_d
        C_r += (1/|A_d|) * X_check^T Phi_d[A_d]
    C_obs_r = C_r / (number of valid days)
C_obs = mean_r C_obs_r

# ---- N1: circular-shift null ----
for b in 1..B_null:
    s = S[b]
    for d in 1..T:
        I = U_{pi_s(d)} ∩ U_d
        if |I| < rho_min * min(...): skip
        X_check = X_{pi_s(d)} re-standardized on I
        C_b += (1/|I|) * X_check^T Phi_d[I]
    C_null[b] = C_b / (number of valid days)

# ---- N2: stratified within-date permutation (slow-moving features only) ----
for b in 1..B_null:
    for d in 1..T:
        for g in strata(d):
            randomly permute X's rows within U_{d,g}
        C_b += (1/N_d) * X_perm_d^T Phi_d
    C_null2[b] = C_b / T
```

Every downstream null-based quantity ($r_{\rm sig}$, feature p-values,
similarity thresholds, fingerprint baselines, breadth and cancellation
baselines) is derived from this same $\{C^{\rm null}_b\}$, so this loop
only needs to run once.
