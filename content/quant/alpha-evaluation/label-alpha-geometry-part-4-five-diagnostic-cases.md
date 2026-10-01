---
title: "Label–Alpha Geometry, Part 4: Five Diagnostic Cases a Model Can Fall Into"
description: "How the three comparison layers chain together, why nonlinearity has two completely different meanings that this method finally separates, and the five case patterns — from 'deploy the ridge surrogate instead' to 'you may have found a new mechanism' — that cover almost every real model diagnosis."
date: 2026-09-14
series: ["Label-Alpha Geometry"]
part: 4
tags: ["alpha-evaluation", "cross-sectional", "methodology", "quant-research", "model-diagnostics"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Label-Alpha Geometry" >}}

Continues from
[Part 3](/quant/alpha-evaluation/label-alpha-geometry-part-3-linear-surrogate-and-residual-alpha/),
which defines $R^2_{\rm lin,OOS}$, $\rho^{\rm rank}_{\rm lin}$, and
$Q_\varepsilon$/$\tilde\nu_\varepsilon$ used throughout this post.

## §59. How the three layers relate logically

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
the surrogate (§60, Case C).

Another common mistake is looking only at layer two: $A^{\mathsf D}_{p,\alpha}=0.95$
gets read as "the model simply learned feature $p$," ignoring that
$\tilde\nu_\alpha$ might be high (most of the model's geometry isn't even
inside $V_r$) and the coincidence effect of low-dimensional spaces (§52(b)).

### Two entirely different meanings of "nonlinear"

This method splits a chronically muddled concept into two measurable
quantities:

$$\boxed{\text{output nonlinearity}:\ 1-R^2_{\rm lin,OOS}\ \text{large}}\qquad\text{vs}\qquad\boxed{\text{predictively novel nonlinearity}:\ Q^{\mathsf D}_\varepsilon\ \text{significant}}$$

The first only says "the model's output isn't a linear function of its
inputs" — automatically true for any tree model or neural network, and
carries no information. The second is what we actually want to know:
**did that nonlinearity buy any new predictive power.**

## §60. Five diagnostic cases

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

**Action:** deploy the original model. Then use the §57 family attribution
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

**Necessary verification (§52(b)):** $A^{\mathsf D}=0.92$ alone isn't rare
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

## §61. Model report

One row per model × horizon band.

**Strength and mechanism**

| Field | Description |
|---|---|
| `Q_alpha_full`, `Q_alpha_D/V/L` | §44 |
| `Q_alpha_inh`, `Q_alpha_nov` | §43.1 |
| `f_alpha_D/V/L` | Raw fingerprint |
| `ft_alpha_D/V/L` | Enrichment fingerprint |
| `nu_alpha`, `nu_tilde_alpha` | §43.3, §43.4 |
| `A_target_alpha`, `eta_D_alpha` | §45 |

**Significance**

| Field | Description |
|---|---|
| `p_alpha_full`, `p_alpha_D/V/L` | §46 |
| `p_alpha_target`, `p_alpha_nov` | §46 |
| `n_model_variants_evaluated` | Rule 46.1 |
| `holdout_used_count` | Rule 41.1 |

**Cross-sectional support**

| Field | Description |
|---|---|
| `b_signal_med` | §47 |
| `b_pred_med`, `b_pred_D_med` | §48 |
| `b_pred/b_signal` | §48's key ratio |
| `R_cancel_alpha`, `R_tilde_alpha` | §32 applied to alpha |

**Traditional metrics**

| Field | Description |
|---|---|
| `IC_mean`, `IC_median`, `IC_hit_rate`, `ICIR` | §49 |

**Linear surrogate and residual**

| Field | Description |
|---|---|
| `R2_lin_selected`, `R2_lin_library` | §55(a) |
| `rho_rank_lin` | §56 |
| `lambda_by_fold` | §55(d) |
| `s2_resid` ($=1-R^2$) | §58 |
| `Q_resid`, `Q_resid_D/V/L` | §58 (post-standardization) |
| `Q_resid_abs` ($=s^2_\varepsilon Q_\varepsilon$) | §58 |
| `A_target_resid`, `eta_D_resid` | §58 |
| `nu_tilde_resid` | §58 |
| `p_resid_D`, `p_resid_nov` | §58 |

**Coordinate-system health**

| Field | Description |
|---|---|
| `delta_drift`, `delta_D/V/L` | §10 |
| `T_eff_D`, `T_eff_V` (OOS period) | §3 |
| `diagnosis_case` | A–E from §60 |

## §62. Feature–alpha report

One row per selected feature (all computed in the same OOS window, §50):

| Field | Description |
|---|---|
| `rho_rank_p_alpha`, `rho_lin_p_alpha` | §51 |
| `A_mech_p_alpha`, `S_mech_p_alpha` | §52 |
| `A_D_p_alpha` | §52 (primary field for directional research) |
| `A_D_rank` | Rank among all selected features |
| `A_D_gap_to_next` | Gap to the runner-up |
| `A_D_null_q95` | 95th percentile of $\max_p|A^{\mathsf D}_{p,\alpha}|$ under the null |
| `delta_R2_family` | §57(a) |
| `R2_own_family` | §57(b) |
| `M_OOS_p` | Feature mechanism OOS stability from §50 |
| `Q_p_OOS`, `A_target_p_OOS` | Feature representation in the OOS period |
| `Q_decay` ($=Q_p^{\rm OOS}/Q_p^{\mathcal E}$) | Feature decay |
| `cell_2x2` | The four-cell classification from §53 |

**Global fields:** $\tilde\nu_\alpha$, $\chi_\alpha$, $m_{\mathsf D}$ (a
reminder of the dimensional limit), and the OOS window's start/end dates.
