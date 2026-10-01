---
title: "Cross-Sectional Predictive Geometry, Part 8: Temporal and Cross-Sectional Stability, Breadth, and Cancellation"
description: "A feature can look identical on Q while hiding completely different failure modes: strength decay versus mechanism drift, both invisible to a p-value. And three diagnostics whose null baseline is never where intuition says it is -- breadth's noise floor is 0.405, not 1; cancellation's noise floor shrinks with the universe size, not fixed at all."
date: 2026-09-14
series: ["Cross-Sectional Predictive Geometry"]
part: 8
tags: ["cross-sectional", "methodology", "quant-research", "predictive-geometry", "stability"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Predictive Geometry" >}}

## §26. Temporal Stability: Bootstrap Design

### Motivation

$Q_p$ is a full-period average. A feature that's very strong in the first
half of the sample and vanishes in the second half can have the same $Q_p$
as a stably-moderate feature, but the former is unusable. We need to
measure whether a conclusion reproduces across different time sub-samples.

### Specification

**Moving block bootstrap (over dates).** Set block length
$\ell_{\rm blk}\ge H_{\max}$; each resample draws
$\lceil T/\ell_{\rm blk}\rceil$ starting points and concatenates them into a
length-$T$ date sequence (with replacement), giving

$$\boxed{C^{(b)}=\frac1{T}\sum_{d\in\text{bootstrap}_b}C^{(d)},\qquad b=1,\dots,B_{\rm time}.}$$

Default $B_{\rm time}=500$–$1000$.

**Must use moving block, not iid bootstrap**: $C^{(d)}$ is autocorrelated in
time (features persist, return volatility clusters). An iid resample would
destroy this structure and understate variance, making every $\pi$ look
over-optimistic. The block length must be at least $H_{\max}$ to ensure
overlapping label windows fall within the same block.

**$\Phi$ and $V_r$ are held fixed across all resamples.** This is
deliberate: this bootstrap measures "**given the coordinate system**, how
stable is the feature-level association." The coordinate system's own
uncertainty is measured by two other mechanisms: §7.2's spectral bootstrap
(stability of $\Sigma_Y$) and §10's $\delta_{\rm drift}$ (cross-period
drift).

> **Rule 26.1**: $\pi$ and $M$ are **lower bounds** on total uncertainty. The
> report must present $\delta_{\rm drift}$ and $K_{\rm eff}$'s bootstrap
> stability alongside them, or a reader will mistake "stable given the
> coordinate system" for "stable overall."

---

## §27. Temporal Stability Metrics

### Definitions

All metrics are **block-specific**: $\mathsf B$ must match the research
target (directional research always uses $\mathsf D$).

**Survival rate (strength)**:

$$\boxed{\pi^{{\rm time},\mathsf B}_p=\frac1{B_{\rm time}}\sum_b\mathbf 1\left[Q^{\mathsf B,(b)}_p>Q^{\mathsf B,{\rm null}}_{p,0.95}\right],}$$

where $Q^{\mathsf B,\rm null}_{p,0.95}$ is that feature-block's null
distribution's 95th percentile.

**Sign stability rate**:

$$\boxed{\pi^{\rm time,sign}_p=\frac1{B_{\rm time}}\sum_b\mathbf 1\left[\operatorname{sign}\big(A^{{\rm tgt},(b)}_p\big)=\operatorname{sign}\big(A^{\rm tgt}_p\big)\right].}$$

**Mechanism stability**:

$$M^{{\rm time},\mathsf B}_{p,b}=\frac{c_p^{\mathsf B,(b)\top}c_p^{\mathsf B}}{\|c_p^{\mathsf B,(b)}\|\,\|c_p^{\mathsf B}\|+\epsilon},\qquad\boxed{M^{{\rm time},\mathsf B}_p=\operatorname{median}_bM^{{\rm time},\mathsf B}_{p,b}.}$$

**Coefficient of variation of strength**:

$$\boxed{CV^{Q,\mathsf B}_p=\frac{\operatorname{sd}_b\big(Q^{\mathsf B,(b)}_p\big)}{\operatorname{mean}_b\big(Q^{\mathsf B,(b)}_p\big)+\epsilon}.}$$

### What Each of the Four Quantities Answers

| Quantity | Question | What failure means |
|---|---|---|
| $\pi^{\rm time,\mathsf B}$ | does strength persistently clear noise across sub-samples | effective only in a minority of periods |
| $\pi^{\rm time,sign}$ | is the direction stable | direction can flip — the most dangerous failure |
| $M^{\rm time,\mathsf B}$ | is the **mechanism** (direction of $c$) stable | what's being predicted is itself changing |
| $CV^{Q,\mathsf B}$ | how much does strength fluctuate | high $CV$ with high $\pi$: strength swings but never disappears |

The distinction between $M$ and $\pi$ matters especially. $\pi^{\rm time}$
only checks "is strength big enough"; $M^{\rm time}$ checks "has the
direction changed." A feature can have $\pi=0.95$ (always strong) but
$M=0.2$ (predicting a different mechanism in every sub-sample) — meaning
it's an unstable composite signal whose large $Q$ is only because it
associates with *something different* each period. Such a feature is
entirely unusable, and looking only at $\pi$ would miss it.

### Relationship to the P-Value

$\pi^{\rm time,\mathsf B}$ and $\tilde p_p$ look redundant but answer
different questions: $\tilde p_p$ asks "could the full-sample point
estimate plausibly be chance," $\pi$ asks "how often does this association
show up across time sub-samples." A feature that's genuine but only
effective 30% of the time can have a tiny $\tilde p$ and a $\pi$ of just
0.3. Both must be part of eligibility.

### Note: the $\pi^{\rm time,\mathsf B}$ Threshold Comparison

In $\pi$'s definition, the bootstrap sample's $Q$ is compared against the
**full-sample** null's 95th percentile. Strictly, a bootstrap sub-sample has
fewer effective periods, so its $Q$'s variance should be somewhat larger
than the full sample's — making this comparison slightly conservative.
This is acceptable: $\pi\ge0.8$ is already a screening threshold, not an
exact test. A more precise version would resample the null alongside every
bootstrap draw (multiplying cost by $B_{\rm time}$), but in practice the
gain doesn't justify the cost.

---
## §28. Two Failure Modes of Instability

The combination of $\pi$ and $M$ splits instability into two categories,
requiring completely different treatment.

### Mode One: Strength Decay (low $\pi^{\rm time,\mathsf D}$, high $M^{\rm time,\mathsf D}$)

The feature's predicted **mechanism hasn't changed** — its strength just
drops below noise in some periods.

Typical causes: a market regime change (volatility, liquidity) altered the
signal's detectability; or the signal is being gradually arbitraged away
(alpha decay).

Treatment: such a feature can still be used, but needs more information to
decide how. The key diagnostic is **whether the decay is monotonic**:
regress $Q^{\mathsf D,(b)}_p$ against the bootstrap samples' average date;
if the slope is significantly negative, this is alpha decay and should be
downweighted or retired; if strength correlates with an observable regime
variable (VIX, volume, spread), this is a conditional signal that should be
paired with a regime filter.

### Mode Two: Mechanism Drift (low $M^{\rm time,\mathsf D}$)

The feature predicts **different things** in different periods. $c_p$'s
direction is itself changing.

Typical causes: the feature's definition changed with the data supply
(e.g. a data field's coverage expanded partway through); or the feature is
a mixture of multiple underlying effects whose relative strength shifts
over time.

Treatment: such a feature **should not be used**, even if $Q_p$ is large and
$\tilde p_p$ is tiny. The reason: what it will predict during OOS is
completely unpredictable. The correct fix is to go back to the feature
construction layer and decompose it into more basic components, testing
each separately.

$M^{\rm time,\mathsf D}<0.3$ should be treated as a hard exclusion
condition, not a scoring input.

### Sign Flips ($\pi^{\rm time,sign}$ low)

This is a special case of mode two, but flagged separately because it's the
most dangerous. A feature whose sign can flip turns profit directly into
loss in a position. Features with $\pi^{\rm time,sign}<0.8$ are always
retired, and **cannot** be rescued by "flipping the sign back in OOS" (§21).

---

## §29. Cross-Sectional Stability

### Motivation

Temporal stability asks "does this hold across different periods."
Cross-sectional stability asks "is it driven by only a handful of stocks."
The two are independent failure dimensions: a feature can be effective
every single day ($\pi^{\rm time}=1$) while every day's effect comes from
the same 20 small-cap names.

### Specification

Keep all dates. For each date $d$, re-bootstrap tickers within
$\mathcal U_d$. If sector/cap metadata is available, use a **stratified**
bootstrap:

$$\mathcal U_d=\bigsqcup_g\mathcal U_{d,g},$$

drawing with replacement $N_{d,g}$ times within each stratum, ending with
$N_d$ rows total.

**Stratification is necessary**: without it, some bootstrap samples would
randomly over-represent a sector, so $C^{(d,b)}$'s variance would reflect
both "sector-composition sampling noise" and "the feature's cross-sectional
concentration" at once. Stratifying fixes sector composition in place, so
the result reflects only the latter.

### Every Resample Must Be Re-Centered

After sampling with replacement, $X$'s mean on the resampled set is
**no longer zero** (some names get drawn multiple times). By Lemma 4.1,
$C^{(d)}$ being a clean covariance depended exactly on $X^\top\mathbf 1=0$.
Once that precondition breaks,

$$\frac1{N_d}X_{d,b}^\top\Phi_{d,b}=\underbrace{\operatorname{Cov}_b(X,\Phi)}_{\text{wanted}}+\underbrace{\bar X_{d,b}\,\bar\Phi_{d,b}^\top}_{\text{spurious signal}}.$$

For a feature unrelated to sampling probability, the spurious term is just
extra noise; but for a feature correlated with stratum size (market cap,
liquidity, sector dummies), it produces a systematic bias.

> **Specification 29.1**: after every bootstrap resample, $X$ must be
> re-demeaned and re-standardized on the resampled sample:
> $$\check X_{idp}=\frac{X_{idp}-\bar X^{(b)}_{dp}}{\sigma^{(b)}_{dp}+\epsilon}.$$

The same requirement applies to §13.3(b)'s permutation intersection. Both
are applications of the same principle.

### Output

$$C^{(b)}_{\rm XS}=\frac1T\sum_d\frac1{N_d}\check X_{d,b}^\top\Phi_{d,b},$$

holding $V_r$ fixed, giving — exactly as in §27 —

$$\pi^{{\rm XS},\mathsf B}_p,\qquad \pi^{\rm XS,sign}_p,\qquad M^{{\rm XS},\mathsf B}_p.$$

Default $B_{\rm XS}=500$.

---

## §30. Predictive Breadth

### Motivation

$Q_p$ is an average, and doesn't tell us how many stocks contributed to it.
An IC contributed to evenly by 3000 stocks and an IC contributed to by 20
stocks (with the other 2980 contributing nothing) are completely different
objects in trading terms: the latter has tiny capacity, is extremely
sensitive to a single event, and is very likely to be eaten by trading
costs in live trading.

### Definition

For feature $p$, let $\hat c_p=c_p/\|c_p\|$. Stock $i$'s **signed
contribution** to the predictive relationship on date $d$ is

$$\boxed{g_{id,p}=X^{XS}_{id,p}\big\langle\Phi_{id,:},\hat c_p\big\rangle.}$$

This has an exact meaning: $\frac1{N_d}\sum_ig_{id,p}$ is precisely $c_p$'s
that-period contribution along the $\hat c_p$ direction, so

$$\frac{1}{T}\sum_d\frac1{N_d}\sum_ig_{id,p}=\hat c_p^\top\,c_p^{\rm full}\approx\|c_p\|.$$

Absolute contribution weight and effective number of contributors:

$$u_{id,p}=\frac{|g_{id,p}|}{\sum_j|g_{jd,p}|},\qquad\boxed{N^{\rm eff}_{d,p}=\frac1{\sum_iu_{id,p}^2},\qquad 1\le N^{\rm eff}_{d,p}\le N_d,}$$

$$\boxed{b_{d,p}=\frac{N^{\rm eff}_{d,p}}{N_d}\in(0,1].}$$

The report uses

$$b^{\rm med}_p=\operatorname{median}_db_{d,p},\qquad b^{\rm p10}_p=\text{the 10th percentile}_d(b_{d,p}).$$

$N^{\rm eff}$ is the inverse-Herfindahl of the $\{|g|\}$ distribution — the
"effective number of contributing stocks."

---

## §31. Reading Breadth and Its Null Baseline

### §31.1 Why $b=1$ Is Not the Reference Point

Intuitively, $b=1$ (every stock contributes equally) feels like the
reference point for "fully diversified." That intuition is wrong.

Consider a feature with **zero** predictive power. Then
$g_{id}=X_{id}\cdot\langle\Phi_{id},\hat c\rangle$ is the product of two
approximately independent, approximately standard-normal quantities. As
$N\to\infty$,

$$b=\frac{\big(\frac1N\sum|g|\big)^2}{\frac1N\sum g^2}\longrightarrow\frac{(\mathbb E|g|)^2}{\mathbb E[g^2]}=\frac{(\mathbb E|u|\cdot\mathbb E|v|)^2}{\mathbb E[u^2]\mathbb E[v^2]}=\frac{(2/\pi)^2}{1}\approx0.405.$$

(Using $\mathbb E|u|=\sqrt{2/\pi}$ for a standard normal.)

So **pure noise has breadth around 0.405, not 1**. (If the contribution
were a single normal quantity rather than a product of two, the baseline
would be $2/\pi\approx0.637$.) Seeing $b=0.35$ and calling it "highly
concentrated" is wrong — that's essentially the noise level.

### §31.2 Relative Breadth

Real $X$ and $\Phi$ are not independent normals ($X$ usually has fat tails,
$\Phi_{\mathsf V}$'s components are quadratic forms and thus heavily
skewed), so 0.405 is only an order-of-magnitude anchor, not a hard
threshold. The report uses a **relative** quantity instead:

$$\boxed{\tilde b_p=\frac{b^{\rm med}_p}{\operatorname{median}_b\big(b^{{\rm med},{\rm null}}_{p,b}\big)},}$$

with the null coming from §13's same resamples. $\tilde b_p\approx1$ means
concentration is no different from noise; $\tilde b_p\ll1$ means genuine
concentration.

**Why only the relative quantity is meaningful**: $b$'s null baseline
differs by block. The $\mathsf V$ block's $\Phi$ components are quadratic
forms of $z$, with much higher kurtosis than normal, so its null breadth is
notably below 0.405. Using the same absolute threshold for $\mathsf D$ and
$\mathsf V$ block features would systematically flag every $\mathsf V$
block feature as "over-concentrated."

### §31.3 Concentration Is Not Necessarily a Defect

> **Rule 31.1**: breadth is a soft diagnostic, not a hard gate.

Concentration itself isn't necessarily a flaw — some genuine mechanisms
(event-driven, earnings surprise, index reconstitution) naturally affect
only a handful of names. Hard-excluding them systematically biases toward
mediocre, broad-based signals.

Flagging rules:

- $\tilde b_p<0.5$ → `CONCENTRATED`
- $\tilde b^{\rm p10}_p<0.2$ → `EXTREME_CONCENTRATION`

These flags indicate **what to check next**, not an automatic exclusion:

1. Check whether the top-contributing names share a common trait (same
   sector? same event type? same data-vendor anomaly?).
2. Check whether $Q_p$ remains significant after excluding the top 1% of
   contributing names.
3. Estimate the actually tradeable capacity.

### Worked Example

| feature | $b^{\rm med}$ | null $b$ | $\tilde b$ | reading |
|---|---:|---:|---:|---|
| momentum | 0.42 | 0.41 | 1.02 | broad-based, indistinguishable from noise baseline |
| earnings surprise | 0.06 | 0.40 | 0.15 | highly concentrated, flag and check capacity |
| volatility prediction | 0.20 | 0.22 | 0.91 | looks concentrated, actually just the $\mathsf V$ block's natural baseline |

The third row demonstrates why relativization is necessary: at a 0.405
threshold this feature would be mislabeled as concentrated.

---

## §32. Cancellation Ratio

### Definition

$$\boxed{R^{\rm canc}_{d,p}=\frac{\big|\sum_ig_{id,p}\big|}{\sum_i|g_{id,p}|+\epsilon}\in[0,1].}$$

### Division of Labor With Breadth

Both are derived from $\{g_{id}\}$, but ask different questions:

- **breadth**: how spread out are contributions across stocks? (the
  concentration of $|g|$)
- **cancellation**: do the contributions cancel each other out? (the sign
  consistency of $g$)

A feature can have high breadth but low cancellation — all 3000 stocks
contribute, but half positive and half negative, with the net effect coming
from a small imbalance. Such a signal's defining feature: trading costs
will eat almost all of it, because you have to trade all 3000 names to
capture that tiny net value.

### The Null Scale

As with breadth, $R^{\rm canc}=1$ is not the reference point. Under the
null, $g_{id}$ is symmetric with zero mean, so

$$\Big|\sum_ig_{id}\Big|\sim\big|\mathcal N(0,N\sigma_g^2)\big|\ \Longrightarrow\ \mathbb E\Big|\sum_ig\Big|=\sigma_g\sqrt{\frac{2N}\pi},$$

and $\sum_i|g|\approx N\,\mathbb E|g|$. So

$$\mathbb E\big[R^{\rm canc}\big]\approx\frac{\sigma_g}{\mathbb E|g|}\sqrt{\frac2{\pi N}}=\frac\pi2\sqrt{\frac2{\pi N}}=\sqrt{\frac{\pi}{2N}},$$

(the last step uses the product case $g=uv$: $\sigma_g=1$,
$\mathbb E|g|=2/\pi$).

$$N=500\Rightarrow0.056,\qquad N=1500\Rightarrow0.032,\qquad N=3000\Rightarrow0.023.$$

So the raw $R^{\rm canc}$ isn't just far below 1 — it **shrinks with
$N_d$**, and isn't comparable across dates or universes.

### Relative Cancellation

$$\boxed{\tilde R_p=\frac{\operatorname{median}_dR^{\rm canc}_{d,p}}{\operatorname{median}_b\big(\operatorname{median}_dR^{{\rm canc},{\rm null}}_{d,p,b}\big)}.}$$

$\tilde R_p\approx1$ means the net effect is indistinguishable from pure
noise — **a serious warning sign**: it means $Q_p$ may come from a handful
of extreme names rather than a systematic cross-sectional relationship.
$\tilde R_p\ge3$ means there's a genuine, sign-consistent cross-sectional
effect.

The report lists both the raw $R^{\rm canc}$ (for capacity and
cost estimation, since absolute value is what matters for trading
economics) and $\tilde R_p$ (for statistical interpretation).

---

## §33. Effective Spectral Rank

### Definition

$$e_j=\frac{s_j^2}{\sum_ks_k^2},\qquad\boxed{r_{\rm eff}=\frac1{\sum_je_j^2}.}$$

This is the participation ratio of the singular-value energy distribution:
if energy concentrates in a single mode, $r_{\rm eff}\to1$; if $r$ modes
share equal energy, $r_{\rm eff}=r$.

### Use

$r_{\rm eff}$ and $r_{\rm sig}$ answer different questions:

- $r_{\rm sig}$: **how many directions are statistically above noise** (a
  test result).
- $r_{\rm eff}$: **how many directions the energy is actually spread
  across** (a descriptive quantity).

Comparing the two is itself informative:

| Case | Meaning |
|---|---|
| $r_{\rm eff}\approx1\ll r_{\rm sig}$ | the whole library essentially illuminates one mechanism; the rest of the significant modes are weak. The library is too homogeneous — add new data sources. |
| $r_{\rm eff}\approx r_{\rm sig}$ | several mechanisms of comparable strength; good library diversity. |
| $r_{\rm eff}\gg r_{\rm sig}$ | energy spread across many weak directions. Usually means the overall signal is weak, $V_r$'s estimate is unstable, and downstream conclusions should be downgraded. |

$r_{\rm eff}$ also gives a lower-bound reference for feature budget: if
$r_{\rm eff}=3$, selecting 40 features is necessarily highly redundant at
the mechanism level — redundancy itself isn't wrong (spreading noise has
value), but it should be a conscious decision, not a mistaken belief that
40 genuinely different signals were found.
