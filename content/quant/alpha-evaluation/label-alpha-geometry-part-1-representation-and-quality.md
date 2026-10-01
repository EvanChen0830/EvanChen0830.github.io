---
title: "Label–Alpha Geometry, Part 1: Leakage Control, Representation Conventions, and Alpha Quality"
description: "Six objects that can leak beyond the model's own parameters, the full/inherited/novel decomposition that every alpha diagnostic is built on, and why a model's fingerprint matters more than its headline IC."
date: 2026-09-14
series: ["Label-Alpha Geometry"]
part: 1
tags: ["alpha-evaluation", "cross-sectional", "methodology", "quant-research", "leakage-control"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Label-Alpha Geometry" >}}

Feature evaluation asks whether a feature carries signal. This series asks
the same question one layer downstream: once a model has turned selected
features into an output, how much of that output is real, what kind of
prediction is it actually making, and how much of it is genuinely new
versus a repackaging of what the feature layer already found.

## §41. Model training and leakage control

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
2. $\Theta_{\mathsf V},\Theta_{\mathsf L},\Gamma_{\mathsf{DV}},\Gamma_{\mathsf{DL}},\Gamma_{\mathsf{VL}}$ (orthogonalization coefficients);
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

## §42. Alpha preprocessing

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

## §43. Representation convention: full / inherited / novel

**This section is the foundation for all of Parts 6 and 7.** It defines how
every signal — feature, alpha, residual — is represented in $\Phi$.

### §43.1 Three vectors

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

### §43.2 Why alpha can't be judged on its projected version alone

$V_r$'s content is "**the label mechanisms the original feature library can
see linearly.**" If a nonlinear model learns

$$\hat\alpha=f(X)$$

and this nonlinear combination produces a new predictive mechanism

$$v_{\rm new}\perp\operatorname{span}(V_r),$$

then computing $c_\alpha^{\rm full}P_V$ **kills it outright.**

Asking afterward "did the model learn any new predictive structure" is then
unanswerable on a representation that has already projected the new
structure away — a logical contradiction.

It's worse for residual alpha (§58): the entire **purpose** of
$\varepsilon_\alpha=\tilde\alpha-\text{Linear}(X)$ is to find whatever lies
outside the linear feature span. Projecting it back onto the feature-derived
subspace negates its own definition.

### §43.3 Novelty fraction

$$\boxed{\nu_s=\frac{Q^{\rm nov}_s}{Q^{\rm full}_s}\in[0,1].}$$

$\nu_\alpha\approx0$: the model is almost entirely recombining predictive
mechanisms the feature library already knows about.
$\nu_\alpha=0.4$: 40% of the model's predictive geometry sits in directions
the linear feature analysis couldn't originally see.

### §43.4 Novelty needs null calibration

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
$Q^{\rm nov}_\alpha$ (§46): high $\tilde\nu$ with a non-significant
$Q^{\rm nov}$ means the new component is just noise.

### §43.5 Summary table of representation conventions

> **Convention 43.2:**
> - **Feature-level selection statistics** (§16–§39) use the **projected**
>   $c_p=c_p^{\rm full}P_V$. Reason: the point of selection is to compare
>   features along statistically validated mechanism directions, and $P_V$
>   is a necessary denoising step. The report also lists $Q_p^{\rm full}$
>   for detecting isolated mechanisms.
> - **Alpha and residuals** (§43–§60) always use $c^{\rm full}$ as the
>   primary result, reporting the inherited/novel decomposition and
>   $\tilde\nu$ alongside it.
> - **Feature-alpha cross-comparison** (§52) uses $c_p$ and
>   $c^{\rm inh}_\alpha$ — both live inside $\operatorname{span}(V_r)$, so
>   the angle between them is meaningful — while also reporting
>   $\tilde\nu_\alpha$ to show how much of alpha this comparison actually
>   covers.
> - **Target alignment** (§45) uses the $\mathsf D$ component of
>   $c^{\rm full}_\alpha$, because the target direction $r_{\mathsf D}$ is
>   defined in the full $\mathsf D$ space, independent of $V_r$.

### §43.6 The boundary of this decomposition

$Q^{\rm nov}$ measures what's "inside $\operatorname{span}(\Phi)$ but
outside $\operatorname{span}(V_r)$." **Anything entirely outside
$\operatorname{span}(\Phi)$ is still invisible** — that would require
extending the §8 lifting family, which this decomposition cannot address.

So the correct statement is:

$$\nu_\alpha\ \text{measures "novelty relative to the feature library," not "absolute novelty."}$$

## §44. Alpha Quality and fingerprint

### Definition

$$\boxed{Q_\alpha=Q^{\rm full}_\alpha=\|c^{\rm full}_\alpha\|^2,}$$

with block components $Q^{\mathsf D}_\alpha,Q^{\mathsf V}_\alpha,Q^{\mathsf L}_\alpha$
also taken from $c^{\rm full}_\alpha$, satisfying the additive identity.

$$\boxed{f_\alpha=\left(\frac{Q^{\mathsf D}_\alpha}{Q_\alpha},\frac{Q^{\mathsf V}_\alpha}{Q_\alpha},\frac{Q^{\mathsf L}_\alpha}{Q_\alpha}\right),}$$

along with the enrichment version $\tilde f_\alpha$ (§19.3, null baseline
from the alpha permutation in §46).

### $\tilde f_\alpha$ is the single most important diagnostic

$\tilde f_\alpha$ answers a question traditional IC analysis cannot see at
all: **what kind of predictive power does this model have?**

| Dominant component of $\tilde f_\alpha$ | Meaning | Consequence for position design |
|---|---|---|
| $\mathsf D$ | Model predicts direction | Can directly build a long-short ranking position |
| $\mathsf V$ | Model predicts activity | A long-short ranking becomes a long-gamma position — entirely different return character |
| $\mathsf L$ | Model predicts direction during large moves | Needs conditional position sizing; exposure should be reduced in calm periods |

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

## §45. Alpha target alignment

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
