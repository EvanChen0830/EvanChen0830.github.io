---
title: "Label–Alpha Geometry"
description: "Judging alpha quality after it leaves the model: leakage control, the full/inherited/novel decomposition, significance and breadth, whether a model's complexity earned anything over a linear surrogate, the five diagnostic cases, and complete execution order."
date: 2026-09-14
tags: ["alpha-evaluation", "methodology", "quant-research", "predictive-geometry", "linear-surrogate"]
math: true
aliases:
  - "/quant/alpha-evaluation/label-alpha-geometry-part-1-representation-and-quality/"
  - "/quant/alpha-evaluation/label-alpha-geometry-part-2-significance-and-inheritance/"
  - "/quant/alpha-evaluation/label-alpha-geometry-part-3-linear-surrogate-and-residual-alpha/"
  - "/quant/alpha-evaluation/label-alpha-geometry-part-4-five-diagnostic-cases/"
  - "/quant/alpha-evaluation/label-alpha-geometry-part-5-execution-order-and-validation/"
---

Feature evaluation asks whether a feature carries signal. This series asks
the same question one layer downstream: once a model has turned selected
features into an output, how much of that output is real, what kind of
prediction is it actually making, and how much of it is genuinely new
versus a repackaging of what the feature layer already found.

## 1. Model training and leakage control

Once feature selection is finished, the model is trained using only
information available in the training period, producing

$$\boxed{\hat\alpha_{id}\quad\text{(must be a strict OOS prediction).}}$$

Under walk-forward evaluation, each test period's predictions can only come
from a model fit before that period.

### Six objects that must be inside the OOS boundary

"Strict OOS" isn't only about the model's own parameters. At least six
objects in this method are fit on the research period, and every one of
them can leak:

1. $B_Y$ (label whitening) and $K_{\rm eff}$;
2. $\Theta_{\mathsf V},\Theta_{\mathsf T},\Gamma_{\mathsf{DV}},\Gamma_{\mathsf{DT}},\Gamma_{\mathsf{VT}}$ (orthogonalization coefficients);
3. $V_r$ and $r_{\rm sig}$;
4. The feature eligibility list;
5. Mechanism-group and raw-family structure;
6. The selected feature set.

> **Rule 41.1.** Even if the model itself is completely OOS, if the selected
> feature set was chosen using data that includes the test period, the
> entire evaluation is invalid.

This is the most common leakage path in iterative research: a researcher
runs feature selection once on the full sample, then starts walk-forward
modeling afterward. The §15 discovery/evaluation split handles the
statistical-inference layer of selection effects, but it does **not**
automatically handle this layer — if the final model's test period overlaps
$\mathcal E$, the selection decisions made on $\mathcal E$ are themselves a
leak.

Two correct approaches, chosen by cost:

- **Full version:** fold feature selection itself into walk-forward,
  re-running §13–§39 independently for every fold. Expensive, but the only
  strictly correct option.
- **Practical version:** hold out a final segment that never participates
  in any selection or modeling decision, and use it **exactly once**. The
  number of times it's used must be logged in the report.

### Purging and embargo

Because labels carry a lookahead window of $H_{\max}$ and features carry
memory $\ell_X$, the gap between train and test must be at least

$$\boxed{\Delta=H_{\max}+\ell_X}$$

remove training samples whose label window overlaps the test period
(purge), and leave an additional $\ell_X$ embargo before the test period
begins.

**This is the same number as the minimum shift in §13.2.** It must be read
from the same configuration constant in implementation — otherwise the same
data gets treated by two different leakage standards in two places, and
that inconsistency is invisible in any report.

## 2. Alpha preprocessing

### Definition

$$\boxed{\tilde\alpha_{id}=\frac{\hat\alpha_{id}-\bar\alpha_d}{s_{\alpha,d}+\epsilon},}\qquad \bar\alpha_d=\frac1{N_d}\sum_i\hat\alpha_{id},\quad s_{\alpha,d}^2=\frac1{N_d}\sum_i(\hat\alpha_{id}-\bar\alpha_d)^2.$$

Then, same as a feature, a single global normalization is applied so that
$\|\tilde\alpha\|_W=1$.

### What this step discards

Within-date standardization removes two dimensions of information:

- **The absolute magnitude of the model's prediction** (e.g. "the model
  thinks today's overall opportunity is large"). This is a market-timing
  component, outside this method's scope (§5).
- **Conviction changes across periods** (e.g. "the model's cross-sectional
  dispersion today is three times its usual level").

The second point is practically valuable and should be logged separately as
a time series $s_{\alpha,d}$ and compared against realized return
dispersion: if the two correlate, the model has a conditional
conviction-calibration ability — an additional benefit this method doesn't
otherwise measure.

Standardization is necessary because it's what makes $Q_\alpha$ and $Q_p$
directly comparable — both become "predictive energy per unit of signal
variance."

## 3. Representation convention: full / inherited / novel

**This section is the foundation for all of Parts 6 and 7.** It defines how
every signal — feature, alpha, residual — is represented in $\Phi$.

### 3.1 Three vectors

For any within-date-standardized signal $s$ with $\|s\|_W=1$:

$$\boxed{c^{\rm full}_s=s^\top W\Phi\ \in\mathbb R^{d_\Phi}}$$

This is $s$'s **full representation**, with no projection applied.

Let $P_V=V_rV_r^\top$ be the projection onto the feature library's
significant mechanism subspace (§15), and define

$$\boxed{c^{\rm inh}_s=c^{\rm full}_sP_V,\qquad c^{\rm nov}_s=c^{\rm full}_s(I-P_V).}$$

with corresponding energies

$$Q^{\rm full}_s=\|c^{\rm full}_s\|^2,\qquad Q^{\rm inh}_s=\|c^{\rm inh}_s\|^2,\qquad Q^{\rm nov}_s=\|c^{\rm nov}_s\|^2.$$

**Proposition 43.1 (orthogonal decomposition).**

$$Q^{\rm full}_s=Q^{\rm inh}_s+Q^{\rm nov}_s.$$

*Proof.* $P_V$ and $I-P_V$ are complementary orthogonal projections, so
$c^{\rm inh}\perp c^{\rm nov}$, and Pythagoras gives the result. $\blacksquare$

### 3.2 Why alpha can't be judged on its projected version alone

$V_r$'s content is "**the label mechanisms the original feature library can
see linearly.**" If a nonlinear model learns

$$\hat\alpha=f(X)$$

and this nonlinear combination produces a new predictive mechanism

$$v_{\rm new}\perp\operatorname{span}(V_r),$$

then computing $c_\alpha^{\rm full}P_V$ **kills it outright.**

Asking afterward "did the model learn any new predictive structure" is then
unanswerable on a representation that has already projected the new
structure away — a logical contradiction.

It's worse for residual alpha (Section 18): the entire **purpose** of
$\varepsilon_\alpha=\tilde\alpha-\text{Linear}(X)$ is to find whatever lies
outside the linear feature span. Projecting it back onto the feature-derived
subspace negates its own definition.

### 3.3 Novelty fraction

$$\boxed{\nu_s=\frac{Q^{\rm nov}_s}{Q^{\rm full}_s}\in[0,1].}$$

$\nu_\alpha\approx0$: the model is almost entirely recombining predictive
mechanisms the feature library already knows about.
$\nu_\alpha=0.4$: 40% of the model's predictive geometry sits in directions
the linear feature analysis couldn't originally see.

### 3.4 Novelty needs null calibration

The exact dimensional effect from fingerprints (§19.2) reappears here, and
worse.

Under the null ($\alpha$ unrelated to $\Phi$), $c^{\rm full}_\alpha$ is
approximately isotropic in the $d_\Phi$-dimensional space, so

$$\mathbb E_0[\nu_\alpha]\approx\frac{d_\Phi-r_{\rm sig}}{d_\Phi}.$$

At $d_\Phi=18$, $r_{\rm sig}=3$, $\mathbb E_0[\nu_\alpha]\approx0.83$. That
is: **a completely non-predictive random signal has a novelty fraction of
0.83.** Reading $\nu_\alpha=0.4$ as "the model has 40% new content" is a
serious misread — 0.4 is in fact well *below* the null baseline, meaning
the model's geometry is **more concentrated in the known subspace than
random.**

So define the enrichment version:

$$e^{\rm inh}_s=\frac{Q^{\rm inh}_s}{\operatorname{median}_bQ^{{\rm inh},{\rm null}}_{s,b}},\qquad e^{\rm nov}_s=\frac{Q^{\rm nov}_s}{\operatorname{median}_bQ^{{\rm nov},{\rm null}}_{s,b}},$$

$$\boxed{\tilde\nu_s=\frac{e^{\rm nov}_s}{e^{\rm inh}_s+e^{\rm nov}_s}.}$$

Under the null, $\mathbb E[\tilde\nu_s]\approx0.5$.

**Reports always use $\tilde\nu$ as the basis for interpretation**; $\nu$ is
only a descriptive quantity for the energy decomposition.

| $\tilde\nu_\alpha$ | Interpretation |
|---:|---|
| $<0.3$ | The model's predictive geometry is concentrated in mechanisms the feature library already knows |
| $\approx0.5$ | New and known components have comparable enrichment |
| $>0.7$ | The model's predictive power mostly comes from directions **outside** known mechanisms |

A high $\tilde\nu$ must always be read alongside the significance test on
$Q^{\rm nov}_\alpha$ (Section 6): high $\tilde\nu$ with a non-significant
$Q^{\rm nov}$ means the new component is just noise.

### 3.5 Summary table of representation conventions

> **Convention 43.2:**
> - **Feature-level selection statistics** (§16–§39) use the **projected**
>   $c_p=c_p^{\rm full}P_V$. Reason: the point of selection is to compare
>   features along statistically validated mechanism directions, and $P_V$
>   is a necessary denoising step. The report also lists $Q_p^{\rm full}$
>   for detecting isolated mechanisms.
> - **Alpha and residuals** (Section 3–Section 20) always use $c^{\rm full}$ as the
>   primary result, reporting the inherited/novel decomposition and
>   $\tilde\nu$ alongside it.
> - **Feature-alpha cross-comparison** (Section 12) uses $c_p$ and
>   $c^{\rm inh}_\alpha$ — both live inside $\operatorname{span}(V_r)$, so
>   the angle between them is meaningful — while also reporting
>   $\tilde\nu_\alpha$ to show how much of alpha this comparison actually
>   covers.
> - **Target alignment** (Section 5) uses the $\mathsf D$ component of
>   $c^{\rm full}_\alpha$, because the target direction $r_{\mathsf D}$ is
>   defined in the full $\mathsf D$ space, independent of $V_r$.

### 3.6 The boundary of this decomposition

$Q^{\rm nov}$ measures what's "inside $\operatorname{span}(\Phi)$ but
outside $\operatorname{span}(V_r)$." **Anything entirely outside
$\operatorname{span}(\Phi)$ is still invisible** — that would require
extending the §8 lifting family, which this decomposition cannot address.

So the correct statement is:

$$\nu_\alpha\ \text{measures "novelty relative to the feature library," not "absolute novelty."}$$

## 4. Alpha Quality and fingerprint

### Definition

$$\boxed{Q_\alpha=Q^{\rm full}_\alpha=\|c^{\rm full}_\alpha\|^2,}$$

with block components $Q^{\mathsf D}_\alpha,Q^{\mathsf V}_\alpha,Q^{\mathsf T}_\alpha$
also taken from $c^{\rm full}_\alpha$, satisfying the additive identity.

$$\boxed{f_\alpha=\left(\frac{Q^{\mathsf D}_\alpha}{Q_\alpha},\frac{Q^{\mathsf V}_\alpha}{Q_\alpha},\frac{Q^{\mathsf T}_\alpha}{Q_\alpha}\right),}$$

along with the enrichment version $\tilde f_\alpha$ (§19.3, null baseline
from the alpha permutation in Section 6).

### $\tilde f_\alpha$ is the single most important diagnostic

$\tilde f_\alpha$ answers a question traditional IC analysis cannot see at
all: **what kind of predictive power does this model have?**

| Dominant component of $\tilde f_\alpha$ | Meaning | Consequence for position design |
|---|---|---|
| $\mathsf D$ | Model predicts direction | Can directly build a long-short ranking position |
| $\mathsf V$ | Model predicts activity | A long-short ranking becomes a long-gamma position — entirely different return character |
| $\mathsf T$ | Model predicts direction during large moves (tail regime) | Needs conditional position sizing; exposure should be reduced in calm periods |

A model with $\tilde f_\alpha=(0.2,0.7,0.1)$ can look fine on traditional IC
(because volatility correlates with return skew, leaking a positive rank
IC), yet it isn't actually a directional model. Using it as one will
realize a return distribution that seriously diverges from the backtest.

**This is this method's most direct practical value:** it diagnoses, before
deployment, the common failure where "the model looks like it has IC but
realized positions underperform expectations."

### Reference scale

$Q_\alpha$ and $Q_p$ share the same scale (§11). For a good ensemble model,
$Q^{\mathsf D}_\alpha$ should be clearly higher than the best single
feature's $Q^{\mathsf D}_p$ — if not, the model gained nothing from
combining features, which is itself a conclusion. Typical values: best
single feature $Q^{\mathsf D}_p\approx5\times10^{-4}$ (IC 0.022), a good
model $Q^{\mathsf D}_\alpha\approx2\times10^{-3}$ (IC 0.045).

## 5. Alpha target alignment

### Definition

$$\boxed{A^{\rm tgt}_\alpha=\big(c^{{\rm full},\mathsf D}_\alpha\big)^\top r_{\mathsf D},\qquad \eta^{\mathsf D}_\alpha=\frac{(A^{\rm tgt}_\alpha)^2}{Q^{\mathsf D}_\alpha+\epsilon}\in[0,1].}$$

### A complete directional diagnostic, made of three numbers

$$Q^{\mathsf D}_\alpha\ \text{(how much directional information)},\quad A^{\rm tgt}_\alpha\ \text{(is the direction right)},\quad \eta^{\mathsf D}_\alpha\ \text{(how much of it is aimed at the target)}.$$

This triple distinguishes three failure modes that look identical under a
single IC number:

| Case | $Q^{\mathsf D}_\alpha$ | $\eta^{\mathsf D}_\alpha$ | $A^{\rm tgt}_\alpha$ | Diagnosis |
|---|---:|---:|---:|---|
| Healthy | High | High | Large, positive | Model is strongly directional and aimed at the target |
| Misaligned | High | Low | Small | Model is directional, but predicting a different horizon combination |
| Inverted | High | High | Large, negative | Model is systematically inverted — check label sign or training-target setup |

The second case is common in multi-horizon training: if the loss weights
every horizon equally while $\Sigma_Y$'s spectrum varies a lot across
horizons, the model automatically skews toward the highest-variance
horizon, which may not be the one the user actually wants to trade. The fix
is to re-weight the training target's horizons to match $a$.

---

Continues directly from
[Section 1–45](#41-model-training-and-leakage-control),
which defines $c^{\rm full}_\alpha$, $Q_\alpha$, $f_\alpha$, and the
inherited/novel decomposition used throughout this post.

## 6. Alpha significance

### Specification

Reuses the §13 permutation exactly: replace $X_d$ with $\tilde\alpha_d$,
destroy pairing with a circular date shift, apply the §13.3 matched
thinning and intersection re-standardization, and obtain
$\{c^{\rm full,null}_{\alpha,b}\}$.

From this, compute all six p-values:

$$p^{\rm full}_\alpha,\quad p^{\mathsf D}_\alpha,\quad p^{\mathsf V}_\alpha,\quad p^{\mathsf T}_\alpha,\quad p^{\rm tgt}_\alpha\ \big((A^{\rm tgt}_\alpha)^2\big),\quad p^{\rm nov}_\alpha\ \big(Q^{\rm nov}_\alpha\big).$$

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

## 7. Alpha signal breadth

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

## 8. Alpha predictive breadth

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

## 9. Traditional IC must still be reported

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

## 10. Feature and alpha must be compared in the same OOS window

### Specification

For the OOS period, using the **same** $\Phi_{\rm OOS}$ (obtained via the
§10 transform) and the **same** $V_r$, recompute

$$\boxed{c_p^{\rm OOS}=\left[\frac1{T_{\rm OOS}}\sum_{d}\frac1{N_d}X^{\rm OOS\top}_{d,p}\Phi^{\rm OOS}_d\right],\qquad c_\alpha^{\rm OOS}=\tilde\alpha^{\rm OOS\top}W_{\rm OOS}\Phi^{\rm OOS}.}$$

Every feature-alpha comparison (Section 11–Section 13) must use these two OOS
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

## 11. Feature–alpha layer one: raw output exposure

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
not be conflated (Section 13).

## 12. Feature–alpha layer two: mechanism inheritance

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
only informative when $|\mathcal F|\lt K$ or the mechanism vectors are
degenerate. This is exactly why "what did the model inherit" ultimately has
to be answered by the linear surrogate of Section 14–Section 17, not by geometric
projection alone.

## 13. Joint interpretation of layers one and two

$\rho^{\rm rank}_{p,\alpha}$ (output similarity) and $A^{\mathsf D}_{p,\alpha}$
(mechanism similarity) form a $2\times2$:

| | $A^{\mathsf D}$ high | $A^{\mathsf D}$ low |
|---|---|---|
| $\rho^{\rm rank}$ high | **Model ≈ a monotone transform of this feature.** Check whether the model is just wrapping a single feature; if so, its complexity isn't earning its keep. | **Output looks similar, mechanism doesn't.** Usually means $\rho$ is propped up by shared cross-sectional structure (industry, market cap) while the actual predictive content differs. Recompute $\rho$ after controlling for known exposures. |
| $\rho^{\rm rank}$ low | **Same mechanism, different output.** The most valuable case: the model reaches the same predictive mechanism by a different route (different name weighting, different nonlinearity), possibly with better noise or cost characteristics. | **The two are unrelated.** If this holds for every feature and $\tilde\nu_\alpha$ is also high, the model's predictive power comes from outside the library — proceed to the residual analysis of Section 18. |

The bottom-left cell deserves special note. Its typical cause: the model
learned to "only use feature $p$ under certain conditions." Output ranking
then differs from $X_p$ (the model doesn't bet when the condition isn't
met), but the predictive mechanism is the same (it predicts the same thing
when the condition is met). Models like this often have a better
$b^{\rm pred}/b^{\rm sig}$ ratio (Section 8) — a case where complexity genuinely
earns value.

---

## 14. The third-layer problem

{{< svgfig src="alpha-eval/fig-alpha-three-layers.svg" >}}
**Figure 1.** The three-layer feature–alpha attribution framework. Layer 1 measures surface-level output correlation; Layer 2 measures geometric mechanism inheritance; Layer 3 fits an out-of-sample linear surrogate and tests whether the residual $\varepsilon_\alpha = \alpha - \hat\alpha$ carries genuine novel signal or merely repackages existing features.
{{< /svgfig >}}

The first two layers (Section 11 output similarity, Section 12 mechanism similarity) are
both geometric, correlational. Neither can answer the question that
actually matters:

> **Did the model's complexity buy anything?**

The right control for this question is: **the same inputs, a linear
model.** If a ridge regression can reconstruct 95% of the model's output
OOS, the model's nonlinearity is just reproducing a linear function — the
complexity bought nothing, and deploying ridge instead would be far better
on operating cost, interpretability, and stability.

So the third layer's design is: build an OOS linear surrogate, measure how
much of the model it reconstructs, and analyze whether **the part it
couldn't reconstruct** has any predictive power.

## 15. Walk-forward linear surrogate

### Specification

Split the OOS alpha period into consecutive folds over time. For fold $k$,
fit on $Train_k^{\rm sur}$ and predict on $Test_k^{\rm sur}$:

$$\hat\beta_k=\arg\min_\beta\left[\sum_{(d,i)\in Train_k^{\rm sur}}w_{id}\big(\tilde\alpha_{id}-X_{id}\beta\big)^2+\lambda\|\beta\|_2^2\right],\qquad \hat\alpha^{\rm lin}_{id}=X_{id}\hat\beta_k\ \ \text{for }(d,i)\in Test_k^{\rm sur}.$$

Concatenating every fold's predictions gives the complete OOS linear
surrogate $\hat\alpha^{\rm lin}$.

### Four specification points that must be made explicit

**(a) Which $X$ to use.** The primary report uses **the model's actual
input set** (the `selected` set from §39). Reason: the question is "did the
model's nonlinearity add value," and its correct control is "same inputs, a
linear model."

Also report, as a secondary figure, the version using **every eligible
feature**, which answers a different question: "can this alpha be linearly
reproduced by the entire library?" Name the two

$$R^2_{\rm lin,selected},\qquad R^2_{\rm lin,library}.$$

$R^2_{\rm lin,library}\gg R^2_{\rm lin,selected}$ means the model is using
information carried by features that weren't selected — usually meaning §39
was too strict, or the model's actual inputs don't match the selected set
(an implementation bug that should be excluded).

**(b) Standardization.** $X$ uses the §4.1-processed $X^{XS}$ (already
within-date standardized and globally normalized), consistent with the rest
of the pipeline. $\tilde\alpha$ uses the Section 2 version. Standardizing again
inside the surrogate fit is **not allowed** — it would break the
correspondence between $R^2$'s scale and $Q$'s scale.

**(c) Intercept.** Both $X^{XS}$ and $\tilde\alpha$ are already within-date
demeaned, so the weighted-least-squares intercept is identically zero. **Do
not fit an intercept.** This must be explicitly turned off in code (most
libraries default it on), or it introduces a degree of freedom that absorbs
day-level common components.

**(d) Choosing $\lambda$.** Determine it **inside** $Train_k^{\rm sur}$ via
time-ordered CV (split by time, never randomly), never using any
information from $Test_k^{\rm sur}$. Each fold determines its own $\lambda_k$,
and $\{\lambda_k\}$ is logged in the report — if $\lambda_k$ varies by more
than an order of magnitude across folds, the surrogate itself is unstable,
and subsequent $R^2$ readings should be downgraded accordingly.

### Why ridge, not OLS

Selected features are necessarily correlated (§36 only removes duplicates
above $S^{\rm raw}>0.9$; correlations in the $0.6$–$0.9$ range are common).
OLS coefficients are highly unstable under these conditions — $\hat\beta_k$
would swing wildly between folds, so the surrogate's OOS performance would
be driven by coefficient noise rather than model structure. Ridge trades a
small amount of bias for a large reduction in variance, giving a more
stable answer to "can this be reconstructed."

Note the surrogate's purpose is **not** to estimate each feature's
contribution (that's Section 17's job), but to estimate **the reconstructible
fraction**. The coefficients themselves are not interpreted.

## 16. Reconstruction metrics

### Definition

$$\boxed{R^2_{\rm lin,OOS}=1-\frac{\sum_{d,i}w_{id}\big(\tilde\alpha_{id}-\hat\alpha^{\rm lin}_{id}\big)^2}{\sum_{d,i}w_{id}\tilde\alpha_{id}^2}=1-\big\|\varepsilon_\alpha\big\|_W^2,}$$

(the last equality uses $\|\tilde\alpha\|_W=1$), and at the ranking level,

$$\rho^{\rm rank}_{\rm lin}=\operatorname{median}_d\operatorname{Spearman}\big(\tilde\alpha_d,\hat\alpha^{\rm lin}_d\big).$$

### The division of labor between the two metrics

$R^2$ measures **numerical** reconstruction; $\rho^{\rm rank}$ measures
**ranking** reconstruction. The two can differ substantially:

- $R^2=0.5$ but $\rho^{\rm rank}=0.9$: the model and the surrogate produce
  almost the same ranking, differing mainly in tail magnitude. Since
  positions are usually determined by ranking, the model's real-world value
  is close to the surrogate's in this case.
- $R^2=0.8$ but $\rho^{\rm rank}=0.6$: numerically close but ranked
  differently — usually means the model and surrogate diverge in the
  middle of the distribution (most names), which happens to be where
  position weight is densest.

**The report must show both side by side**, with $\rho^{\rm rank}_{\rm lin}$
as the primary criterion for "can the model be replaced by the surrogate,"
since it's closer to trading semantics.

### OOS $R^2$ can be negative

If the surrogate fails badly on some folds, $R^2_{\rm lin,OOS}<0$ is
possible, meaning the surrogate is worse than predicting zero everywhere.
This should not be read as "the model is highly nonlinear" — first check
the surrogate's implementation ($\lambda$ too small, leakage or
discontinuity in the fold split).

## 17. Which feature families reconstruct alpha

### Specification

A single ridge coefficient is not a suitable importance measure (correlated
features make coefficients unstable and arbitrarily interchangeable).
Instead use two complementary quantities, **both of which must be
reported.**

**(a) Conditional contribution (drop-family):** for redundancy family
$R_g$, remove every feature in that family, rerun the surrogate through the
identical pipeline, and get $R^2_{-g}$:

$$\boxed{\Delta R^2_g=R^2_{\rm full}-R^2_{-g}.}$$

Meaning: this family's marginal contribution **given that every other
family is present.**

**(b) Standalone contribution:** fit the surrogate using only family
$R_g$'s features:

$$\boxed{R^2_{g,{\rm own}}.}$$

Meaning: how much this family can reconstruct on its own.

### Why both are needed

Two families that are highly correlated with each other (even if
$S^{\rm raw}<0.9$ and so weren't merged) will mask each other:

$$\Delta R^2_g\approx0\quad\text{and}\quad\Delta R^2_{g'}\approx0,\qquad\text{yet}\quad R^2_{g,{\rm own}}=R^2_{g',{\rm own}}=0.6.$$

Looking only at $\Delta R^2$ would conclude "neither family matters," when
the truth is "either one alone suffices, but you can't drop both."

Interpretation rule:

| $\Delta R^2_g$ | $R^2_{g,\rm own}$ | Interpretation |
|---|---|---|
| High | High | Core family, irreplaceable |
| Low | High | Redundant with another family; usable as a backup/substitute |
| High | Low | Only useful in combination with other families (interaction or residual correction) |
| Low | Low | No contribution to linear reconstruction |

### A more complete version

If there aren't too many families ($\le12$), compute family-level Shapley
values: compute $R^2$ for every subset, average the marginal contributions.
Cost is $2^G$ surrogate fits — about 1000 at $G=10$, practically
affordable. Shapley values resolve the masking problem above and give the
cleanest attribution, but don't scale to large $G$.

## 18. Residual alpha

### Definition

$$\boxed{\varepsilon_\alpha=\tilde\alpha-\hat\alpha^{\rm lin}.}$$

Since both $\tilde\alpha$ and $X^{XS}$ are already within-date demeaned,
$\varepsilon_\alpha$ automatically satisfies $\sum_i\varepsilon_{\alpha,id}=0$.

**Two scales must be recorded separately.**

**Absolute residual energy:**

$$s^2_\varepsilon=\|\varepsilon_\alpha\|_W^2=1-R^2_{\rm lin,OOS}.$$

**Standardized residual signal:** within-date standardize $\varepsilon_\alpha$
and globally normalize, giving $\tilde\varepsilon_\alpha$ with
$\|\tilde\varepsilon_\alpha\|_W=1$, then compute per Section 3

$$c^{\rm full}_\varepsilon=\tilde\varepsilon_\alpha^\top W\Phi,\qquad c^{\rm inh}_\varepsilon=c^{\rm full}_\varepsilon P_V,\qquad c^{\rm nov}_\varepsilon=c^{\rm full}_\varepsilon(I-P_V),$$

$$Q_\varepsilon=\|c^{\rm full}_\varepsilon\|^2,\quad Q^{\mathsf D}_\varepsilon,\ Q^{\mathsf V}_\varepsilon,\ Q^{\mathsf T}_\varepsilon,\quad A^{\rm tgt}_\varepsilon,\ \eta^{\mathsf D}_\varepsilon,\quad\tilde\nu_\varepsilon.$$

### Why $Q_\varepsilon$ must be computed after standardizing

Using $\varepsilon_\alpha$ at its raw scale gives $\|\varepsilon\|_W^2=1-R^2$, so

$$Q_\varepsilon^{\rm raw}\le1-R^2.$$

At $R^2=0.9$, $Q_\varepsilon^{\rm raw}$'s upper bound is only 0.1 times
$Q_\alpha$'s upper bound. Comparing the two numbers directly would conclude
"the residual has no predictive power," when that's purely scale
arithmetic, unrelated to the residual's actual quality.

After standardizing, $Q_\varepsilon$ and $Q_\alpha$ sit on **the same
scale** (the IC scale of §11) and can be compared directly. The relationship
between them is

$$\boxed{\underbrace{s^2_\varepsilon\cdot Q_\varepsilon}_{\text{residual's absolute contribution to total predictive energy}}\quad\text{vs}\quad\underbrace{Q_\alpha}_{\text{total predictive energy}}.}$$

**The report must list both $Q_\varepsilon$ (quality) and
$s^2_\varepsilon Q_\varepsilon$ (absolute contribution).** The former
answers "is the residual signal good"; the latter answers "how much of the
total it accounts for."

### Significance

Test $\tilde\varepsilon_\alpha$ using the same §13 permutation set, giving
$p^{\mathsf D}_\varepsilon$, $p^{\rm nov}_\varepsilon$, etc.

Worth noting: even though standardization is what makes $Q_\varepsilon$'s
**magnitude** comparable, the permutation p-value itself was always
scale-invariant — it compares the observed value of the same standardized
signal against its own null distribution. So the significance conclusion is
identical before and after standardizing; standardizing only changes which
magnitudes are comparable.

### This section is this method's final product

The combination of $Q^{\mathsf D}_\varepsilon$ and $\tilde\nu_\varepsilon$
answers a question that's extremely hard to answer in quantitative research
generally, and that this method answers directly:

> **Does the part of the model's output that lies outside the linear
> feature span have any predictive value at all? And if it does, is it a
> re-weighting of known mechanisms, or an entirely new one?**

$Q^{\mathsf D}_\varepsilon$ significance answers the first half;
$\tilde\nu_\varepsilon$ answers the second.

---

Continues from
[Section 14–58](#54-the-third-layer-problem),
which defines $R^2_{\rm lin,OOS}$, $\rho^{\rm rank}_{\rm lin}$, and
$Q_\varepsilon$/$\tilde\nu_\varepsilon$ used throughout this post.

## 19. How the three layers relate logically

$$\underbrace{\rho^{\rm rank}_{p,\alpha}}_{\text{layer 1: output}}\ \longrightarrow\ \underbrace{A^{\mathsf D}_{p,\alpha},\ \tilde\nu_\alpha}_{\text{layer 2: mechanism}}\ \longrightarrow\ \underbrace{R^2_{\rm lin},\ Q_\varepsilon,\ \tilde\nu_\varepsilon}_{\text{layer 3: beyond-linear structure}}$$

The three layers **tighten progressively**:

- Layer one can be computed without $\Phi$ at all, but only sees surface
  similarity.
- Layer two needs $\Phi$ and $V_r$, and can see predictive content, but
  cannot see whether the complexity was necessary.
- Layer three needs an additional surrogate fit, and can see the marginal
  value of complexity, but is only meaningful relative to a linear control.

**No layer substitutes for another.** A common mistake is looking only at
layer three: $R^2_{\rm lin}=0.3$ gets read as "the model is highly
nonlinear and valuable." But if $Q_\varepsilon$ isn't significant, that 70%
unreconstructed part is just noise, and the model is in fact *worse* than
the surrogate (Section 20, Case C).

Another common mistake is looking only at layer two: $A^{\mathsf D}_{p,\alpha}=0.95$
gets read as "the model simply learned feature $p$," ignoring that
$\tilde\nu_\alpha$ might be high (most of the model's geometry isn't even
inside $V_r$) and the coincidence effect of low-dimensional spaces (Section 12(b)).

### Two entirely different meanings of "nonlinear"

This method splits a chronically muddled concept into two measurable
quantities:

$$\boxed{\text{output nonlinearity}:\ 1-R^2_{\rm lin,OOS}\ \text{large}}\qquad\text{vs}\qquad\boxed{\text{predictively novel nonlinearity}:\ Q^{\mathsf D}_\varepsilon\ \text{significant}}$$

The first only says "the model's output isn't a linear function of its
inputs" — automatically true for any tree model or neural network, and
carries no information. The second is what we actually want to know:
**did that nonlinearity buy any new predictive power.**

## 20. Five diagnostic cases

The following five cases cover the large majority of model diagnostic
conclusions encountered in practice.

### Case A: the model is equivalent to its linear surrogate

$$R^2_{\rm lin,OOS}=0.90,\quad \rho^{\rm rank}_{\rm lin}=0.94,\quad Q^{\mathsf D}_\varepsilon\ \text{not significant}.$$

**Conclusion:** alpha can almost entirely be reconstructed by a stable
linear combination of features; whatever lies outside the linear span has
no significant directional predictive value.

**Action:** deploy the ridge surrogate, not the original model.

This is **not** a failure of the model — it's an extremely valuable, clear
answer. It saves the team substantial operating cost, removes a source of
overfitting risk, and greatly improves interpretability and
monitorability. Being able to make this call with confidence is itself one
of this method's main payoffs.

### Case B: the nonlinearity is valuable

$$R^2_{\rm lin,OOS}=0.40,\quad Q^{\mathsf D}_\varepsilon\ \text{significant},\quad A^{\rm tgt}_\varepsilon>0,\quad \tilde\nu_\varepsilon\approx0.4.$$

**Conclusion:** the model produces substantial structure that a stable
linear feature span cannot reconstruct, and that structure still has
positive directional predictive value toward the intended target; it
mostly sits in a re-weighting of known mechanisms (since $\tilde\nu_\varepsilon$
isn't high).

**Action:** deploy the original model. Then use the Section 17 family attribution
and conditional analysis to try to understand the *form* of this nonlinear
structure — if it can be understood (e.g. "the model only uses momentum
features on high-volatility days"), it can often be hard-coded as a new
feature, making it linearly available in the next round. This is the
standard path for converting model capability into a permanent library
asset.

### Case C: the nonlinearity is pure cost

$$R^2_{\rm lin,OOS}=0.30,\quad Q^{\mathsf D}_\varepsilon\ \text{not significant}.$$

**Conclusion:** the model's output has substantial nonlinear/non-stationary
structure, but that structure provides no statistically supportable
directional predictive information.

**Action:** this is the worst case — the model has introduced a large
amount of useless variance, raising turnover, transaction cost, and
unpredictability with no compensating return. Strengthen regularization, or
fall back to the surrogate directly.

**Further diagnosis:** check $Q^{\mathsf V}_\varepsilon$. If
$Q^{\mathsf V}_\varepsilon$ is significant while $Q^{\mathsf D}_\varepsilon$
isn't, the model's nonlinear component is predicting volatility. This is
usually caused by the loss function: a model trained with MSE on raw
returns automatically devotes more capacity to high-volatility names,
because they contribute the most loss. The fix is to switch to a
volatility-normalized training target, or train directly on $y^*$ (§20) as
the label.

### Case D: strong mechanism inheritance

$$A^{\mathsf D}_{p,\alpha}=0.92,\quad\text{clearly higher than the runner-up},\quad \Delta R^2_g\ \text{large},\quad R^2_{g,\rm own}\ \text{large}.$$

**Conclusion:** the final alpha strongly inherits this feature family's
predictive mechanism, and the family also makes a substantial independent
contribution to alpha's linear reconstruction. This is strong evidence "the
model is using this mechanism."

**Necessary verification (Section 12(b)):** $A^{\mathsf D}=0.92$ alone isn't rare
at $m_{\mathsf D}=3$. This conclusion only holds if three conditions are
met: (i) the gap to the runner-up is significant; (ii) it exceeds the 95th
percentile of the null distribution of $\max_p|A^{\mathsf D}_{p,\alpha}|$;
(iii) both $\Delta R^2_g$ and $R^2_{g,\rm own}$ support it. The geometric
angle alone is not sufficient to draw this conclusion.

### Case E: a new mechanism

$$R^2_{\rm lin,OOS}=0.50,\quad Q^{\mathsf D}_\varepsilon\ \text{significant},\quad \tilde\nu_\varepsilon>0.7,\quad p^{\rm nov}_\varepsilon<0.05.$$

**Conclusion:** the model has not only created directional predictive
structure beyond the linear span, but most of that structure sits in a
label direction the existing feature library **cannot see at all.**

**Action:** this is the strongest positive result, and the one that demands
the most caution, because it claims to have found a new mechanism.
Required follow-up verification:

1. Examine $c^{\rm nov}_\varepsilon$'s actual loadings in $\Phi$: which
   block, which horizon combination? Translate it into economic language.
2. Re-verify on a completely independent holdout that never participated in
   any decision (and use it only once).
3. Check $\delta_{\rm drift}$: if the coordinate system itself is drifting,
   the direction of $I-P_V$ might just reflect instability in $\Phi$.
4. Try to construct a new feature that captures this direction. If
   successful, the mechanism converts from "a black-box model's exclusive
   capability" into "a permanent library asset" — it will appear inside
   $V_r$ in the next round.

Step 4 is the ultimate purpose of this method's design: **make whatever a
model discovers systematically recoverable.**

## 21. Model report

One row per model × horizon band.

**Strength and mechanism**

| Field | Description |
|---|---|
| `Q_alpha_full`, `Q_alpha_D/V/L` | Section 4 |
| `Q_alpha_inh`, `Q_alpha_nov` | Section 3.1 |
| `f_alpha_D/V/L` | Raw fingerprint |
| `ft_alpha_D/V/L` | Enrichment fingerprint |
| `nu_alpha`, `nu_tilde_alpha` | Section 3.3, Section 3.4 |
| `A_target_alpha`, `eta_D_alpha` | Section 5 |

**Significance**

| Field | Description |
|---|---|
| `p_alpha_full`, `p_alpha_D/V/L` | Section 6 |
| `p_alpha_target`, `p_alpha_nov` | Section 6 |
| `n_model_variants_evaluated` | Rule 46.1 |
| `holdout_used_count` | Rule 41.1 |

**Cross-sectional support**

| Field | Description |
|---|---|
| `b_signal_med` | Section 7 |
| `b_pred_med`, `b_pred_D_med` | Section 8 |
| `b_pred/b_signal` | Section 8's key ratio |
| `R_cancel_alpha`, `R_tilde_alpha` | §32 applied to alpha |

**Traditional metrics**

| Field | Description |
|---|---|
| `IC_mean`, `IC_median`, `IC_hit_rate`, `ICIR` | Section 9 |

**Linear surrogate and residual**

| Field | Description |
|---|---|
| `R2_lin_selected`, `R2_lin_library` | Section 15(a) |
| `rho_rank_lin` | Section 16 |
| `lambda_by_fold` | Section 15(d) |
| `s2_resid` ($=1-R^2$) | Section 18 |
| `Q_resid`, `Q_resid_D/V/L` | Section 18 (post-standardization) |
| `Q_resid_abs` ($=s^2_\varepsilon Q_\varepsilon$) | Section 18 |
| `A_target_resid`, `eta_D_resid` | Section 18 |
| `nu_tilde_resid` | Section 18 |
| `p_resid_D`, `p_resid_nov` | Section 18 |

**Coordinate-system health**

| Field | Description |
|---|---|
| `delta_drift`, `delta_D/V/L` | §10 |
| `T_eff_D`, `T_eff_V` (OOS period) | §3 |
| `diagnosis_case` | A–E from Section 20 |

## 22. Feature–alpha report

One row per selected feature (all computed in the same OOS window, Section 10):

| Field | Description |
|---|---|
| `rho_rank_p_alpha`, `rho_lin_p_alpha` | Section 11 |
| `A_mech_p_alpha`, `S_mech_p_alpha` | Section 12 |
| `A_D_p_alpha` | Section 12 (primary field for directional research) |
| `A_D_rank` | Rank among all selected features |
| `A_D_gap_to_next` | Gap to the runner-up |
| `A_D_null_q95` | 95th percentile of $\max_p|A^{\mathsf D}_{p,\alpha}|$ under the null |
| `delta_R2_family` | Section 17(a) |
| `R2_own_family` | Section 17(b) |
| `M_OOS_p` | Feature mechanism OOS stability from Section 10 |
| `Q_p_OOS`, `A_target_p_OOS` | Feature representation in the OOS period |
| `Q_decay` ($=Q_p^{\rm OOS}/Q_p^{\mathcal E}$) | Feature decay |
| `cell_2x2` | The four-cell classification from Section 13 |

**Global fields:** $\tilde\nu_\alpha$, $\chi_\alpha$, $m_{\mathsf D}$ (a
reminder of the dimensional limit), and the OOS window's start/end dates.

---

This closing post is the capstone for the entire geometric framework — both
the feature-evaluation layer
([Cross-Sectional Predictive Geometry](/quant/feature-evaluation/) /
[Cross-Sectional Feature Selection](/quant/feature-evaluation/cross-sectional-feature-selection/))
and this alpha-evaluation series. It's reference material: an execution
checklist, a boundary table, a validation-test suite, and a default
parameter table meant to be checked against directly, not read start to
end.

## 23. Complete execution order

The following order is complete and can be used directly as an
implementation checklist. Bracketed references point to the relevant
sections.

**Stage 0: Setup and data**

```text
 0. Lock the configuration: band definition (cumulative/incremental),
    target weight a, Delta = H_max + l_X, split ratios, B_null / B_time /
    B_XS, every threshold                                    [§6,§13,§20,Section 26]
 1. Align by (time, ticker) key; verify the point-in-time contract  [§2]
 2. Split discovery D and evaluation E, with a purge gap of Delta   [§15.2]
 3. Feature preprocessing: winsorize -> within-date standardize
    -> global normalization so ||X_p||_W = 1                       [§4.1]
 4. Expand the transform family; record (source, transform) index  [§4.3]
 5. Handle missing values; record coverage_p and active_ratio_p    [§4.2]
 6. Within-date demean labels; cut horizon bands                   [§5,§6]
```

**Stage 1: Build the coordinate system on discovery**

```text
 7. Sigma_Y; spectral bootstrap; decide K_eff; record lambda_j/lambda_1
    and cond(Sigma_Y)                                              [§7.1,§7.2]
 8. B_Y -> Z; verify Z^T W Z = I holds exactly                      [§7.3,§7.4]
 9. Phi_D := Z (no orthogonalization)                               [§9.2]
10. G_V, G_T -> within-date centering -> Gram-Schmidt -> Phi         [Feature Evaluation §8.1, §9.3]
11. Store B_Y, Theta_V, Theta_T, Gamma_DV, Gamma_DT, Gamma_VT        [Feature Evaluation §10]
12. Compute T_eff^D, T_eff^V                                        [§3.1]
13. C_D; permutation null (common shift, matched thinning,
    intersection re-standardization)                               [§11,§13]
14. Parallel analysis -> r_sig, V_r (fixed from here on)            [§14]
15. r_eff; m_D (effective mechanism dimension)                      [§33,§23.4]
```

**Stage 2: Do all inference on evaluation**

```text
16. Transform the label basis onto E; compute delta_drift(D->E)     [§10]
17. C_E -> c_p = C_E V_r V_r^T -> Q, Q_D, Q_V, Q_T, Q_full           [§15.2,§16]
18. kappa_p; decide whether each feature uses the N1 or N2 null     [§13.4]
19. Permutation on E -> raw p-values -> Westfall-Young adjustment   [§18]
20. Fingerprint f and enrichment f_tilde                            [§19]
21. r_D = normalize(Phi_D^T W y*); A_target; eta_D;
    significance of (A_target)^2                                   [§20,§21,§22]
22. Independent-shift null -> similarity null distribution ->
    zeta(tau); S^D                                                  [§23,§24]
23. Temporal moving-block bootstrap -> pi, M, CV (block-specific)   [§26,§27]
24. Cross-sectional stratified bootstrap (re-center after each
    resample) -> pi^XS, M^XS                                        [§29]
25. Breadth and cancellation; relativize to b_tilde, R_tilde         [§30,§31,§32]
```

**Stage 3: Selection**

```text
26. Eligibility E1-E5                                                [§34]
27. Mechanism grouping (S^D, average linkage, tau=0.90);
    report zeta and m_D                                              [§35]
28. Raw redundancy families (complete linkage);
    apply the transform-family merge rule                            [§36,§4.3]
29. Score_p = median_b (A_target^(b))^2 -> family representatives    [§37]
30. Greedily maximize Cov^feat_D to 0.90 or budget; then fill V, L   [§39]
31. Record Cov^span as a diagnostic                                  [§38.1]
32. Output the feature report                                        [§40]
```

**Stage 4: Model and attribution**

```text
33. Train the model (purge + embargo = Delta; the selected set must
    not contain test information)                                   [Section 1]
34. Within-date standardize alpha + global normalization             [Section 2]
35. OOS: transform Phi; record delta_drift                           [§10]
36. c_alpha^full -> inherited/novel decomposition; Q, f, f_tilde,
    nu_tilde                                                         [Section 3,Section 4]
37. A_target_alpha, eta_D_alpha                                      [Section 5]
38. Alpha permutation -> p_full, p_D/V/L, p_target, p_nov             [Section 6]
39. b_signal, b_pred; IC / ICIR                                       [Section 7,Section 8,Section 9]
40. Recompute c_p^OOS in the same OOS window; M_OOS_p; Q_decay        [Section 10]
41. rho_rank; A^D_{p,alpha} + null max quantile; 2x2 interpretation   [Section 11,Section 12,Section 13]
42. Walk-forward ridge surrogate -> R2_selected, R2_library,
    rho_rank_lin                                                      [Section 15,Section 16]
43. Family attribution: delta_R2 and R2_own                           [Section 17]
44. Residual: s2_resid, Q_resid (post-standardization), Q_resid_abs,
    nu_tilde_resid, p_resid                                           [Section 18]
45. Determine Case A-E                                                [Section 20]
46. Output the model report and feature-alpha report                  [Section 21,Section 22]
```

**Stage 5 (optional): production refit**

```text
47. Refit the label basis and V_r on D union E, for production use
    only. No new decisions may be made after refitting.               [§15.4]
```

## 24. Fit/transform boundary table

The table below lists every fitted object, what data it's fit on, and how
it's applied. Violating any row constitutes a leak.

| Object | Fit on | Applied to | Section |
|---|---|---|---|
| Winsorization quantiles | Computed per period (no cross-period fit) | The current period | §4.1 |
| Within-date $\mu,\sigma$ | Computed per period | The current period | §4.1 |
| Global normalization $\|X_p\|_W$ | $\mathcal D$ | The entire period | §4.1 |
| $\Sigma_Y$, $K_{\rm eff}$, $B_Y$ | $\mathcal D$ | The entire period | §7 |
| Within-date centering of $G_\bullet$ | Computed per period | The current period | §8.1 |
| $\Theta_{\mathsf V},\Theta_{\mathsf T},\Gamma_\bullet$ | $\mathcal D$ | The entire period | §9 |
| $r_{\rm sig}$, $V_r$ | $\mathcal D$ | The entire period | §14, §15 |
| $r_{\mathsf D}$ | $\mathcal D$ (or $\mathcal E$, must be consistent and logged) | The entire period | §20 |
| Eligibility list | $\mathcal E$ | Model input | §34 |
| Grouping / families | $\mathcal E$ | Model input | §35, §36 |
| Selected set | $\mathcal E$ | Model input | §39 |
| Model parameters | $\subseteq\mathcal D\cup\mathcal E$, strictly before OOS | OOS | Section 1 |
| Surrogate $\hat\beta_k,\lambda_k$ | Each fold's train | That fold's test | Section 15 |

### Vectors across folds cannot be averaged directly

Under the full walk-forward feature-selection variant (Rule 41.1's complete
version), each fold has its own $\Phi^{(k)}$ and $V_r^{(k)}$. In that case

$$\frac1K\sum_kc_p^{(k)}\quad\text{is meaningless}$$

— different folds have different coordinate systems, and adding
component-wise mixes unrelated coordinates together.

There are three things you can actually do:

1. **Average scalars:** $\overline{Q_p}=\frac1K\sum_kQ_p^{(k)}$,
   $\overline{|A^{\rm tgt}_p|}$, etc. Scalars are coordinate-independent and
   can be averaged.
2. **Average angles:** $\overline{M_p}=\frac1K\sum_k\cos\big(c_p^{(k)},c_p^{(\rm ref)}\big)$,
   where $c_p^{(\rm ref)}$ comes from a designated reference fold.
3. **Project onto a common coordinate system:** pick one fold's
   $\Phi^{(\rm ref)}$ as reference, use the §10 transform to project every
   fold's data onto it, then average. This is the only correct way to get
   an "average mechanism vector," at the cost of introducing an arbitrary
   choice of reference fold.

## 25. Numerical validation

Split into two groups: **identity tests** (verify the implementation is
correct) and **calibration tests** (verify the statistics are valid). Both
groups must pass.

### 25.1 Identity tests

| Test | Criterion |
|---|---|
| Feature demeaning | $\big|\frac1{N_d}\sum_iX^{XS}_{idp}\big|<10^{-8}$ |
| Feature normalization | $\big|\|X_p\|_W^2-1\big|<10^{-8}$ |
| Label demeaning | $\big|\sum_iY^{XS}_{id,h}\big|<10^{-8}\cdot$ scale |
| Exact whitening | $\|Z^\top WZ-I\|_F<10^{-10}$ (**must hold exactly**, see Proposition 7.1) |
| Basis orthogonality | $\|\Phi^\top W\Phi-I\|_F<10^{-6}$ (in-sample) |
| $C$ computed two ways agrees | Direct long-form vs. $\frac1T\sum_dC^{(d)}$ |
| Quality additivity | $Q_p=Q_p^{\mathsf D}+Q_p^{\mathsf V}+Q_p^{\mathsf T}$ |
| Novelty additivity | $Q^{\rm full}_s=Q^{\rm inh}_s+Q^{\rm nov}_s$ |
| Efficiency bound | $0\le\eta^{\mathsf D}\le1+\epsilon$ |
| Sign flip | $X_q=-X_p\Rightarrow S^{\rm mech}_{pq}=1,\ A^{\rm mech}_{pq}=-1$ |
| Coverage monotonicity | ${\rm Cov}^{\rm feat}(\mathcal F\cup\{p\})\ge{\rm Cov}^{\rm feat}(\mathcal F)$ |
| Coverage submodularity | Draw 100 random $\mathcal A\subset\mathcal B$ pairs, verify diminishing marginal gain (Theorem 38.3) |

Note the quality-additivity test is **tautologically true** (Proposition
16.1, a consequence of coordinate decomposition) — it can only catch
implementation slicing bugs. What actually verifies orthogonality is the
exact-whitening and basis-orthogonality tests.

### 25.2 Calibration tests

Run on synthetic data. This group is the direct evidence for whether
statistical conclusions can be trusted.

**T1 Null calibration ($Q$).** Generate independent $X$ and $Y$ (preserving
real cross-sectional correlation structure and temporal autocorrelation,
but independent of each other), run the full pipeline.

Criterion: p-values roughly uniform; false-positive rate $\le5\%$ at
$\gamma_{\rm FWER}=0.05$.

*This is the single most important test.* If the §15 split or refit
specification is implemented incorrectly (e.g. the null still uses the
observed $V_r$), this test's false-positive rate will be significantly
above 5%.

**T2 Null calibration (similarity).** On the same synthetic data, compute
$\zeta(\tau_{\mathsf D})$.

Criterion: $\zeta\approx0$ (within Monte Carlo error).

**T3 Fingerprint calibration.** $\tilde f_p$ of a pure-noise feature.

Criterion: close to $(1/3,1/3,1/3)$; average deviation of the three
components $<0.05$.

(Cross-check: the same features' raw $f_p$ should be close to the dimension
ratio $(4,10,4)/18$, which is exactly what the §19.2 analysis predicts.)

**T4 Rotation invariance.** Rerun the entire pipeline with 10 different BLAS
thread counts and random seeds.

Criterion: every quantity listed in Proposition 9.2 agrees exactly within
numerical error ($10^{-10}$ relative).

*Must be repeated multiple times.* This test's failure mode is
**intermittent** — passing once does not mean it's correct.

**T5 Novelty calibration.** Generate random $\alpha$ unrelated to $\Phi$.

Criterion: $\tilde\nu_\alpha\approx0.5$; $p^{\rm nov}_\alpha$ uniform.

(Cross-check: raw $\nu_\alpha$ should be close to $(d_\Phi-r_{\rm sig})/d_\Phi$,
verifying Section 3.4.)

**T6 Breadth/cancellation calibration.** Pure-noise feature.

Criterion: $\tilde b_p\approx1$, $\tilde R_p\approx1$.

(Cross-check: raw $b$ should be close to 0.405, raw $R^{\rm canc}$ should be
close to $\sqrt{\pi/2N}$.)

**T7 Target recovery.** Set $\alpha=y^*$.

Criterion: $A^{\rm tgt}_\alpha>0$ and $\eta^{\mathsf D}_\alpha\approx1$;
setting $\alpha=-y^*$ gives $A^{\rm tgt}_\alpha<0$.

**T8 Injected-signal recovery.** Construct
$X_p=\text{IC}\cdot y^*/\|y^*\|+\sqrt{1-\text{IC}^2}\cdot\text{noise}$, with
$\text{IC}\in\{0.01,0.03,0.05\}$.

Criterion: the estimated $|A^{\rm tgt}_p|$ matches the injected IC within
$\pm20\%$. This also validates §11's magnitude scale.

**T9 Dynamic universe.** Randomly vary $N_d$ within $\pm40\%$, randomly add
and remove tickers.

Criterion: the whole pipeline doesn't depend on equal row counts; T1's
calibration still holds under a dynamic universe (this validates §13.3's
matched thinning).

**T10 Leakage detection.** Shift the label forward in time by $2\Delta$
periods as a whole (destroying the real association while keeping all other
structure intact).

Criterion: every $Q$ drops to the null level. If not, there's a leak
somewhere in the pipeline.

**T11 Slow-moving feature.** Construct a persistent feature with
$\kappa_p=0.95$ that is **unrelated to the label**.

Criterion: N1's (circular shift) p-value visibly deviates from uniform
(exactly the power loss §13.4 predicts); N2's (stratified within-date
permutation) p-value is uniform. This contrast validates the necessity of
the dual-null design.

### 25.3 Acceptance criteria

> An implementation's output should not be trusted for any significance
> conclusion until **T1, T2, T3, T4, and T5** all pass. T1 validates the
> inference specification, T2 validates similarity calibration, T3 and T5
> validate dimensional calibration, T4 validates that the coordinate system
> is well-defined. Together these five cover every part of the method that
> can systematically change a conclusion.

## 26. Default parameters

| Parameter | Default | Section |
|---|---:|---|
| $X$ winsorization | 1% / 99% cross-sectional | §4.1 |
| Minimum daily feature coverage | 0.80 | §4.2 |
| Transform family | $\{x,|x|,x|x|\}$ | §4.3 |
| Horizon band convention | Non-overlapping incremental (preferred) | §6 |
| $\kappa_{\max}$ (whitening condition number ceiling) | 100 | §7.2 |
| $\tau_\lambda$ (eigenvalue bootstrap relative SE ceiling) | 0.25 | §7.2 |
| $\epsilon_{\rm rank}$ ($\operatorname{orth}_W$ numerical cutoff) | $10^{-10}$ | §9.1 |
| Discovery/evaluation ratio | 60/40 | §15.2 |
| $\Delta$ (minimum permutation shift = purge gap) | $H_{\max}+\ell_X$ | §13.2, Section 1 |
| $\rho_{\min}$ (minimum intersection ratio) | 0.90 | §13.3 |
| $R_{\rm thin}$ (matched thinning repeats) | 5 | §13.3 |
| $\kappa_p$ threshold (switch to N2 null) | 0.50 | §13.4 |
| $B_{\rm null}$ | 2000 (capped at all legal shifts) | §13.2, §18.3 |
| $\gamma_{\rm mode}$ (mode significance level) | 0.05 | §14 |
| Primary multiple-testing scheme | Westfall–Young, $\gamma_{\rm FWER}=0.05$ | §18.3 |
| Secondary multiple-testing scheme | Permutation FDR, $\gamma_{\rm FDR}=0.05$ | §18.4 |
| $B_{\rm time}$ (temporal bootstrap) | 1000 | §26 |
| Block length $\ell_{\rm blk}$ | $\ge H_{\max}$ | §26 |
| $B_{\rm XS}$ (cross-sectional bootstrap) | 500 | §29 |
| $\tau_{\mathsf D}$ (mechanism-grouping cut point) | 0.90 (= accepting 19% energy loss) | §23.1, §35 |
| Raw redundancy cut point | 0.90 | §36 |
| $\pi$ threshold (three terms) | 0.80 | §34 |
| $M$ threshold (two terms) | 0.60 | §34 |
| Breadth warning | $\tilde b<0.5$ / $\tilde b^{\rm p10}<0.2$ | §31.3 |
| ${\rm Cov}^{\rm feat}_{\mathsf D}$ target | 0.90 | §39 |
| ${\rm Cov}^{\rm feat}_{\mathsf V},{\rm Cov}^{\rm feat}_{\mathsf T}$ target | 0.80 | §39 |
| $m_{\mathsf D}$ warning threshold | $<4$ → grouping has no discriminative power | §23.4 |
| $\delta_{\rm drift}$ warning threshold | $>0.3$ | §10 |
| $T_{\rm eff}/T$ warning threshold | $<0.1$ | §3.1 |
| $\epsilon$ (numerical floor) | $10^{-12}$ | §0.3 |

**These thresholds are operational filtering parameters, not statistical
criteria.** Statistical significance is always determined by permutation
(§14, §18, Section 6). The thresholds exist to filter out features that are
"significant but unusable"; their values can be adjusted for the research
target, but **must be locked and logged before research begins** — adjusting
thresholds after the fact until results look good is a form of overfitting
this method cannot detect.

## 27. The unified representation

### One object

The entire method ultimately reduces to this: for any within-date
standardized cross-sectional signal $s$ with $\|s\|_W=1$,

$$\boxed{c^{\rm full}_s=s^\top W\Phi\ \in\mathbb R^{d_\Phi},\qquad c^{\rm inh}_s=c^{\rm full}_sP_V,\qquad c^{\rm nov}_s=c^{\rm full}_s(I-P_V).}$$

From this single object, every geometric property of the signal can be read
off:

$$\begin{aligned}
&Q_s=\|c^{\rm full}_s\|^2 && \text{strength}\\
&\big(Q^{\mathsf D}_s,Q^{\mathsf V}_s,Q^{\mathsf T}_s\big) && \text{mechanism decomposition}\\
&\tilde f_s && \text{mechanism fingerprint (dimension-calibrated)}\\
&A^{\rm tgt}_s=\big(c^{{\rm full},\mathsf D}_s\big)^\top r_{\mathsf D} && \text{direction}\\
&\eta^{\mathsf D}_s && \text{target precision}\\
&\tilde\nu_s && \text{novelty relative to the feature library}\\
&\pi_s,\ M_s && \text{temporal and cross-sectional stability}\\
&\tilde b_s,\ \tilde R_s && \text{cross-sectional support}
\end{aligned}$$

and the relationship between two signals is

$$A^{\rm mech}_{s_1,s_2}=\hat c_{s_1}^\top\hat c_{s_2},\qquad S^{\rm mech}_{s_1,s_2}=\big|A^{\rm mech}_{s_1,s_2}\big|.$$

### Three applications

**Feature layer:**

$$X_p\ \xrightarrow{\ \Phi(Y)\ }\ c_p\quad\Longrightarrow\quad\text{select features that are real, stable, and mechanism-complementary}.$$

**Alpha layer:**

$$\hat\alpha\ \xrightarrow{\ \Phi(Y)\ }\ c_\alpha\quad\Longrightarrow\quad\text{is the model predicting direction, activity, or large-move direction}.$$

**Cross layer:**

$$c_p\ \longleftrightarrow\ c^{\rm inh}_\alpha\quad\Longrightarrow\quad\text{which known mechanisms did the model inherit},$$

$$\varepsilon_\alpha=\tilde\alpha-\text{OOSLinear}(X)\ \xrightarrow{\ \Phi(Y)\ }\ c_\varepsilon\quad\Longrightarrow\quad\text{what did the model create that's genuinely new}.$$

### The cost of this unification

Putting three questions into the same coordinate system buys comparability,
but it also means **all three questions share the same assumptions:**

1. The space spanned by the three liftings is sufficient to describe the
   predictive patterns we care about (§8);
2. The mapping from signal to label is linear on the feature side, or can
   be covered by the transform family (§4.3);
3. $\Phi$ is sufficiently stable between the research period and the OOS
   period (monitored via $\delta_{\rm drift}$, §10).

If any one of the three fails, the answers to all three questions fail
**simultaneously**, and in the same way. This is the structural risk of a
unified coordinate system, and it's why the three coordinate-system health
indicators — $\delta_{\rm drift}$, the stability of $K_{\rm eff}$, and
$T_{\rm eff}$ — must travel with every single report.

What it buys in return is an ability rare in quantitative research: **the
same ruler, measuring features, measuring the model, and measuring what's
left over once the model's output has the features subtracted out.**

## Appendix A: Method boundaries

This appendix states plainly the questions this method **cannot** answer.
These aren't flaws — they're the scope of applicability, stated explicitly
so that "this method can't see it" isn't misread as "it doesn't exist."

### A.1 The label-side function family is finite

$\Phi$ is spanned by three liftings: $z$, $\operatorname{svec}(zz^\top)$,
$\|z\|z$. Other functions of future returns are not covered, for example:

- The sign of extreme tails (indicator functions like $\mathbf 1[z_h>3]$);
- Third-order-and-above structure (cross-sectional differences in skew or
  kurtosis);
- Cross-day path shape (a band only contains several horizons computed
  from the same starting point).

A signal that predicts only these will get $Q\approx0$ from this method.

The way to extend coverage is to add more liftings (e.g. a third-order
version of $\operatorname{svec}$, or quantile indicator functions), at the
cost of increasing $d_\Phi$, which dilutes the fingerprint's denominator and
reduces statistical power in every direction. This is an explicit
trade-off: **more dimensions means more patterns visible, but each one seen
less clearly.**

### A.2 The signal-side mapping is linear

$c_s=s^\top W\Phi$ is linear in $s$. If $Y=X^2+\text{noise}$ and $X$ is
symmetrically distributed, then $\mathbb E[X\cdot z]=\mathbb E[X^3]=0$, and
the raw $X$'s $Q\approx0$.

The §4.3 transform family pushes the boundary from "is $X$ itself linearly
related" to "is any member of $X$'s specified transform family linearly
related." But it's still not "any function of $X$."

So the correct statement of this method is:

> **It answers "do these representations of this feature carry predictive
> information," not "the raw variable $X$ has no information under any
> possible nonlinear transform."**

This distinction matters especially when the final model is itself
nonlinear: the model may use a transform the transform family doesn't
cover, in which case Section 12's mechanism-inheritance analysis will
underestimate the degree of inheritance, and Section 18's residual analysis will
classify it as a "new mechanism." Step 4 of Case E's verification in Section 20
(trying to construct a new feature that captures that direction) exists
precisely to distinguish these two situations.

### A.3 Circular-shift null power against persistent signals

Fully addressed in §13.4: diagnosed via $\kappa_p$, with stratified
within-date permutation as an alternative null. But it's important to keep
in mind the two nulls test **different hypotheses**:

- N1: "the feature's cross-sectional ranking is unrelated to the return
  mechanism at the correct time";
- N2: "after controlling for strata (industry × market cap), the
  stock-level feature ranking is unrelated to the return mechanism."

N2 is weaker: it cannot detect a purely industry-level or market-cap-level
effect (those effects are preserved within strata). So for a slow-moving
feature, the strongest conclusion this method can offer is "there is
stock-level predictive power," not "there is predictive power."

### A.4 Equal weighting is equal-row weighting, not equal-influence weighting

Fully addressed in §3: diagnosed via $T_{\rm eff}$. The practical
consequence is that $Q$ is not immune to market regime, especially for the
$\mathsf V$ and $\mathsf T$ blocks (whose contribution scales like
$\sim\|z\|^2$).

### A.5 Discriminative power of a low-dimensional mechanism space

Fully addressed in §23.3. When $m_{\mathsf D}\le3$, the question "are these
two features the same mechanism" cannot be answered by the angle between
them alone. This method's design therefore bases selection on a coverage
criterion (§38.2) rather than pairwise thresholds, but **descriptive**
conclusions like "the model inherited feature $p$" are still bound by this
limit, and only hold up when paired with Section 12's null quantile and Section 17's
linear attribution.

### A.6 Association is not causation

A large $Q_p$ only means the feature has a stable cross-sectional
association with the future return mechanism. It may come from a genuine
pricing mechanism, or from a shared exposure to an unmodeled risk factor.
This method does not distinguish the two.

To distinguish them requires adding, after §5, a preprocessing layer that
cross-sectionally regresses out exposure to known factors and keeps the
residual. This doesn't change any downstream formula, but it changes the
meaning of every conclusion: from "this feature predicts returns" to "this
feature predicts returns beyond known factors."

### A.7 $Q$ is not return

$Q$ does not account for transaction cost, capacity, turnover, borrow
availability, or market impact. A feature with high $Q$ but extremely low
$\tilde b$ and $\tilde R$ close to 1 may be entirely untradeable in
practice. The §30–§32 diagnostics exist for exactly this reason, but they
are only warnings — they cannot substitute for an actual cost model and
capacity analysis.

## Appendix B: Interpretation quick-reference

| Observation | Meaning | What to do |
|---|---|---|
| $K_{\rm eff}\lt K$ | Whitening hit an ill-conditioned case; some horizon directions are unstable | Go back to §6, use the incremental basis or merge horizons; truncation is a last resort |
| $T_{\rm eff}/T<0.1$ | Conclusions are driven by a handful of event days | Downgrade all of that block's conclusions to exploratory |
| $\delta_{\rm drift}>0.3$ | The coordinate system itself is drifting | Cross-period comparisons are invalid; rule out coordinate issues before discussing model decay |
| $r_{\rm sig}=1$ | The library illuminates only one mechanism | Add a new data source, not more features of the same type |
| $r_{\rm eff}\gg r_{\rm sig}$ | Energy is spread across weak directions | $V_r$ is unstable; downgrade downstream conclusions |
| $m_{\mathsf D}<4$ | The mechanism space is too small for angles to discriminate | Use grouping only for report organization; let the coverage criterion drive selection |
| $\zeta(\tau)\le0$ | The observed similarity structure is indistinguishable from the null | Same as above; conclusions must not be drawn at the cluster level |
| No feature is significant | Usually $B_{\rm null}$ is too small | Check whether $1/(B_{\rm null}+1)$ is small enough (§18.2) |
| $\tilde f_p\approx(1/3,1/3,1/3)$ | This feature is noise | Regardless of how structured the raw $f$ looks |
| $\tilde b_p\approx1,\tilde R_p\approx1$ | Concentration and net effect are both indistinguishable from noise | Strongly suspect $Q_p$ comes from a handful of extreme names |
| $\pi^{\rm time}$ high but $M^{\rm time}$ low | Always strong, but predicting something different each period | Hard-reject; go back and decompose at the feature-construction layer |
| $\pi^{\rm time,sign}<0.8$ | Sign flips | Reject; never flip sign on OOS |
| $\eta^{\mathsf D}_p\approx0$ but $Q^{\mathsf D}_p$ large | Directional information is pointed at a different horizon combination | Flag `TARGET_MISMATCH`, route to another band |
| ${\rm Cov}^{\rm span}_{\mathsf D}=1$ with only 3 features selected | Span coverage has already saturated | Inevitable consequence of $d_{\mathsf D}=K$; use ${\rm Cov}^{\rm feat}$ instead |
| $\tilde f_\alpha$ dominated by $\mathsf V$ | The model is predicting volatility, not direction | A long-short rank becomes a long-gamma position; change the training target |
| $\nu_\alpha=0.4$ | Below the null baseline $(d_\Phi-r_{\rm sig})/d_\Phi$ | Use $\tilde\nu$ for interpretation; 0.4 actually means "not novel" |
| $b^{\rm sig}$ high but $b^{\rm pred}$ low | Bets broadly, predictive power is concentrated | Threshold low-contribution names to cut cost with almost no return given up |
| $R^2_{\rm lin}$ low but $Q_\varepsilon$ not significant | The nonlinearity is pure cost (Case C) | Strengthen regularization or fall back to the surrogate |
| $R^2_{\rm lin,library}\gg R^2_{\rm lin,selected}$ | The model is using information outside the selected set | Check whether the model's actual inputs match `selected` |
| $\Delta R^2_g\approx0$ but $R^2_{g,\rm own}$ high | Masked by another family | Either one alone suffices, but you can't drop both |
| $A^{\mathsf D}_{p,\alpha}=0.92$ with the top five close together | Inevitable in a low-dimensional space, not inheritance | Needs the null max quantile and Section 17 support together |
| $\tilde\nu_\varepsilon>0.7$ and $p^{\rm nov}_\varepsilon<0.05$ | Possibly a new mechanism (Case E) | Four-step verification; the ultimate goal is building a new feature that captures it |

---

*Every definition, formula, default, and report field in this book forms
one self-consistent specification. Section 25's identity tests and calibration
tests are the acceptance criteria: no significance conclusion should be
trusted until all five core calibration tests (T1–T5) pass.*
