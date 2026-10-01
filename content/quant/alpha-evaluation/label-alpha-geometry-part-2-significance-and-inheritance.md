---
title: "Label–Alpha Geometry, Part 2: Significance, Breadth, and the Two-Layer Feature–Alpha Comparison"
description: "Significance testing for alpha and the hidden multiplicity cost of trying many models, two different breadth diagnostics and why their divergence is the most actionable finding in this section, why traditional IC must still be reported, and the two-layer framework for asking what a model actually inherited from the feature library."
date: 2026-09-14
series: ["Label-Alpha Geometry"]
part: 2
tags: ["alpha-evaluation", "cross-sectional", "methodology", "quant-research"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Label-Alpha Geometry" >}}

Continues directly from
[Part 1](/quant/alpha-evaluation/label-alpha-geometry-part-1-representation-and-quality/),
which defines $c^{\rm full}_\alpha$, $Q_\alpha$, $f_\alpha$, and the
inherited/novel decomposition used throughout this post.

## §46. Alpha significance

### Specification

Reuses the §13 permutation exactly: replace $X_d$ with $\tilde\alpha_d$,
destroy pairing with a circular date shift, apply the §13.3 matched
thinning and intersection re-standardization, and obtain
$\{c^{\rm full,null}_{\alpha,b}\}$.

From this, compute all six p-values:

$$p^{\rm full}_\alpha,\quad p^{\mathsf D}_\alpha,\quad p^{\mathsf V}_\alpha,\quad p^{\mathsf L}_\alpha,\quad p^{\rm tgt}_\alpha\ \big((A^{\rm tgt}_\alpha)^2\big),\quad p^{\rm nov}_\alpha\ \big(Q^{\rm nov}_\alpha\big).$$

Testing $Q^{\rm nov}_\alpha$ **requires $V_r$ to be independent of the
evaluation data** (the §15 split specification) — otherwise the direction
of $I-P_V$ was itself chosen using the leftover of the same data, and the
null will be skewed. This is where the split specification directly pays
off in alpha analysis.

### Model multiplicity

A **single** model needs no multiple-testing adjustment. But in practice a
researcher almost always compares several models or hyperparameter
settings, in which case:

> **Rule 46.1.** Log `n_model_variants_evaluated` — the total number of
> model variants (including abandoned ones) for which $Q_\alpha$ was
> computed on the same evaluation data. If $>1$, either apply a
> Westfall–Young adjustment across the model family, or hold out a
> completely independent final test set and use it exactly once.

This field is routinely omitted, and it is the single largest hidden source
of multiple testing in alpha research. A researcher who tries 50
hyperparameter combinations and reports the best one has an effective
significance level one to two orders of magnitude worse than the nominal
one.

## §47. Alpha signal breadth

### Definition

This section measures the concentration of **the model's output itself**,
unrelated to predictive power:

$$u^{\rm sig}_{id}=\frac{|\tilde\alpha_{id}|}{\sum_j|\tilde\alpha_{jd}|},\qquad\boxed{N^{\rm eff,sig}_d=\frac1{\sum_i(u^{\rm sig}_{id})^2},\qquad b^{\rm sig}_d=\frac{N^{\rm eff,sig}_d}{N_d}.}$$

### Interpretation

A low $b^{\rm sig}$ means the model's output is highly concentrated on a
handful of names — most stocks get a near-zero prediction. This isn't
inherently good or bad, but it determines position structure:

- $b^{\rm sig}\approx0.4$ (typical level): the model has an opinion on the
  whole universe; positions will be broadly diversified.
- $b^{\rm sig}<0.1$: the model is effectively a "stock picker" rather than
  a "ranker," with strong opinions on only a few names. Consider
  constructing positions directly from the top/bottom-$k$ rather than
  linear weighting.

As in §31, the null baseline is not 1. If $\tilde\alpha$ is approximately
cross-sectionally normal, $b^{\rm sig}\to2/\pi\approx0.637$ (a single normal
quantity, not a product of two).

## §48. Alpha predictive breadth

### Definition

Exactly isomorphic to §30, applied to alpha:

$$g_{id,\alpha}=\tilde\alpha_{id}\big\langle\Phi_{id,:},\hat c^{\rm full}_\alpha\big\rangle,\qquad\boxed{N^{\rm eff,pred}_{d,\alpha}=\frac1{\sum_i(u^{\rm pred}_{id,\alpha})^2},\quad b^{\rm pred}_{d,\alpha}=\frac{N^{\rm eff,pred}_{d,\alpha}}{N_d}.}$$

Can also be restricted to the $\mathsf D$ block, giving $b^{{\rm pred},\mathsf D}_{d,\alpha}$.

### The key contrast

The gap between $b^{\rm sig}$ and $b^{\rm pred}$ is an important diagnostic:

| $b^{\rm sig}$ | $b^{\rm pred}$ | Meaning |
|---|---|---|
| High | High | Model bets broadly, and predictive power is broadly distributed too. Good capacity, robust. |
| High | Low | **Model bets broadly, but predictive power comes from only a few names.** Most of the position is dead weight — it contributes cost and risk without contributing return. |
| Low | Low | Model bets narrowly and predictive power is narrow. Consistent, but low capacity and sensitive to a single event. |
| Low | High | Rare. Model bets on only a few names, but predictive power is spread out — usually means $\hat c_\alpha$ is unstably estimated. |

The second row is the most valuable finding, because it points to a
directly fixable problem: applying a position cap or a score threshold to
names contributing near-zero predictive power can cut cost and turnover
substantially with almost no return given up.

## §49. Traditional IC must still be reported

### Specification

Compute the target signal $y^*_d=Y_d^{(\beta)}a$, daily

$$IC_d=\operatorname{Spearman}(\tilde\alpha_d,y^*_d),$$

and report

$$\overline{IC},\quad \operatorname{median}(IC_d),\quad P(IC_d>0),\quad ICIR=\frac{\operatorname{mean}(IC_d)}{\operatorname{std}(IC_d)}.$$

### Why this can't be skipped

This method's geometric quantities are complementary to IC, not a
replacement for it:

1. **IC is rank-based, immune to monotone transforms.** The geometric
   quantities are linear. Disagreement between the two is itself
   informative: if $A^{\rm tgt}_\alpha$ is strong but IC is weak, predictive
   power is concentrated in the tails (linear correlation propped up by a
   few extreme values, without an improvement in ranking); the reverse
   means the relationship is monotone but nonlinear.
2. **IC maps directly to tradeable ranking performance**, and is the common
   language for communicating with outside teams or other methodologies.
3. **IC has a large base of existing empirical benchmarks**, while $Q$'s
   benchmark has to be built within this method itself.

$ICIR$ also supplies something $Q$ cannot: the **temporal stability** of
predictive power on a daily basis, complementary to the bootstrap stability
of §27.

### Technical note

The $IC_d$ series is highly autocorrelated (overlapping label windows), so
$ICIR\cdot\sqrt{T}$ is **not** a valid t-statistic. If IC significance is
needed, use the same §13 permutation set, not a closed-form formula.

## §50. Feature and alpha must be compared in the same OOS window

### Specification

For the OOS period, using the **same** $\Phi_{\rm OOS}$ (obtained via the
§10 transform) and the **same** $V_r$, recompute

$$\boxed{c_p^{\rm OOS}=\left[\frac1{T_{\rm OOS}}\sum_{d}\frac1{N_d}X^{\rm OOS\top}_{d,p}\Phi^{\rm OOS}_d\right],\qquad c_\alpha^{\rm OOS}=\tilde\alpha^{\rm OOS\top}W_{\rm OOS}\Phi^{\rm OOS}.}$$

Every feature-alpha comparison (§51–§53) must use these two OOS
representations.

### Why this rule is non-negotiable

Comparing a research-period $c_p$ against an OOS $c_\alpha$ mixes two
different market periods, different realized returns, different regimes.
The answer to "which feature's mechanism did the model inherit" would then
be contaminated by two factors unrelated to what the model actually
learned:

1. drift in the feature's own mechanism between the two periods;
2. drift in the label covariance structure ($\delta_{\rm drift}$).

Neither has anything to do with what the model learned, yet both change the
angle between the vectors. Computing both on the same window eliminates
both contamination sources entirely — whatever difference remains can only
come from the model itself.

### A free feature-decay diagnostic

The same computation gives an important result for free: comparing
$c_p^{\mathcal E}$ against $c_p^{\rm OOS}$ yields each selected feature's
OOS mechanism stability,

$$M^{\rm OOS}_p=\frac{c_p^{{\rm OOS},\mathsf D\top}c_p^{\mathcal E,\mathsf D}}{\|c_p^{{\rm OOS},\mathsf D}\|\|c_p^{\mathcal E,\mathsf D}\|}.$$

If model performance has degraded and $M^{\rm OOS}_p$ is broadly low, the
problem is in the features, not the model. This distinction is decisive for
where to direct further research.

## §51. Feature–alpha layer one: raw output exposure

### Definition

$$\rho_{p,\alpha,d}=\operatorname{Spearman}\big(X_{d,p},\tilde\alpha_d\big),\qquad\boxed{\rho^{\rm rank}_{p,\alpha}=\operatorname{median}_d\rho_{p,\alpha,d}.}$$

A Pearson version $\rho^{\rm lin}_{p,\alpha}$ can also be computed.

### What this layer asks

> How close is the model's output ranking to some raw feature's ranking?

This is an **output-level** similarity, unrelated to predictive mechanism.
Its uses:

- Detecting whether the model has degenerated into a wrapper around a
  single feature ($\rho^{\rm rank}>0.9$).
- Understanding the position's actual composition: if the model's output is
  highly correlated with market cap, positions carry a systematic size
  exposure regardless of the model's predictive mechanism.

$\rho^{\rm rank}$ and the next layer ask **different questions** and must
not be conflated (§53).

## §52. Feature–alpha layer two: mechanism inheritance

### Definition

Per Convention 43.2, compare inside $\operatorname{span}(V_r)$:

$$\hat c_p=\frac{c_p}{\|c_p\|},\qquad \hat c^{\rm inh}_\alpha=\frac{c^{\rm inh}_\alpha}{\|c^{\rm inh}_\alpha\|},$$

$$\boxed{A^{\rm mech}_{p,\alpha}=\hat c_p^\top\hat c^{\rm inh}_\alpha,\qquad S^{\rm mech}_{p,\alpha}=\big|A^{\rm mech}_{p,\alpha}\big|,}$$

and the block-specific version

$$\boxed{A^{\mathsf D}_{p,\alpha}=\frac{c_p^{\mathsf D\top}c_\alpha^{{\rm inh},\mathsf D}}{\|c_p^{\mathsf D}\|\,\|c_\alpha^{{\rm inh},\mathsf D}\|}.}$$

Directional research always leads with $A^{\mathsf D}_{p,\alpha}$ (Rule 24.1).

### What this layer asks

> Which feature's mechanism is closest to the mechanism the model inherited?

Of the three layers, this is the one that most directly answers "which
feature did the model learn."

### Three things that must be presented alongside it

**(a) Coverage fraction.** $A^{\rm mech}_{p,\alpha}$ only describes alpha's
inherited component. $\tilde\nu_\alpha$ must be reported alongside it,
otherwise a reader will misread "the model's mechanism resembles feature X"
as "the model simply is feature X," when in fact most of the predictive
geometry may sit in the novel component.

**(b) Dimensional calibration.** The low-dimension caveat of §23.3 applies
fully here. At $m_{\mathsf D}=3$, a random direction has a 10% chance of
$|\cos|>0.9$ with any given feature; across $P$ features, the chance that
at least one exceeds 0.9 is close to 1. So:

> **The number $A^{\mathsf D}_{p,\alpha}=0.92$ on its own is not sufficient
> to claim "the model inherited feature $p$."**

What must be compared is the **relative ranking and gap**: rank $|A^{\mathsf D}|$
across all selected features, and look at the gap between first and second
place. If the top five all sit between 0.85 and 0.93, that's the inevitable
consequence of a low-dimensional space, not an inheritance relationship.
Also report the quantile of $\max_p|A^{\mathsf D}_{p,\alpha}|$ under the null
(generated from $\alpha$ via the independent shifts of §23.2).

**(c) Combinations, not a single feature.** A model usually inherits a
**combination** of features, not one. A more informative quantity is how
much of alpha's inherited component can be spanned by the selected
features' mechanisms:

$$\boxed{\chi_\alpha=\frac{\big\|c^{{\rm inh},\mathsf D}_\alpha P^{\mathsf D}_{\mathcal F}\big\|^2}{\big\|c^{{\rm inh},\mathsf D}_\alpha\big\|^2}\in[0,1],}$$

where $P^{\mathsf D}_{\mathcal F}$ is the projection onto the subspace
spanned by the selected features' $\mathsf D$ mechanism vectors (§38.1).
Note that by Proposition 38.2, once $|\mathcal F|\ge K$ with linearly
independent mechanism vectors, $\chi_\alpha\equiv1$ — so $\chi_\alpha$ is
only informative when $|\mathcal F|<K$ or the mechanism vectors are
degenerate. This is exactly why "what did the model inherit" ultimately has
to be answered by the linear surrogate of §54–§57, not by geometric
projection alone.

## §53. Joint interpretation of layers one and two

$\rho^{\rm rank}_{p,\alpha}$ (output similarity) and $A^{\mathsf D}_{p,\alpha}$
(mechanism similarity) form a $2\times2$:

| | $A^{\mathsf D}$ high | $A^{\mathsf D}$ low |
|---|---|---|
| $\rho^{\rm rank}$ high | **Model ≈ a monotone transform of this feature.** Check whether the model is just wrapping a single feature; if so, its complexity isn't earning its keep. | **Output looks similar, mechanism doesn't.** Usually means $\rho$ is propped up by shared cross-sectional structure (industry, market cap) while the actual predictive content differs. Recompute $\rho$ after controlling for known exposures. |
| $\rho^{\rm rank}$ low | **Same mechanism, different output.** The most valuable case: the model reaches the same predictive mechanism by a different route (different name weighting, different nonlinearity), possibly with better noise or cost characteristics. | **The two are unrelated.** If this holds for every feature and $\tilde\nu_\alpha$ is also high, the model's predictive power comes from outside the library — proceed to the residual analysis of §58. |

The bottom-left cell deserves special note. Its typical cause: the model
learned to "only use feature $p$ under certain conditions." Output ranking
then differs from $X_p$ (the model doesn't bet when the condition isn't
met), but the predictive mechanism is the same (it predicts the same thing
when the condition is met). Models like this often have a better
$b^{\rm pred}/b^{\rm sig}$ ratio (§48) — a case where complexity genuinely
earns value.
