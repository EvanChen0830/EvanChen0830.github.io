---
title: "Multi-Sub-Universe Directional Coupling Geometry, Part 4: Stability, OOS Testing, and Redundancy"
description: "Does a candidate universe's geometry survive contact with time it didn't see during discovery? Bootstrap persistence, the frozen train/OOS canonical test, redundancy between candidate universes, and the local refit that lets an accepted universe become its own market domain."
date: 2026-09-14
series: ["Multi-Sub-Universe Directional Coupling Geometry"]
part: 4
tags: ["universe", "cross-sectional", "methodology", "quant-research", "out-of-sample"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Multi-Sub-Universe Directional Coupling Geometry" >}}

## §83. Temporal stability

### §83.1 Moving-block bootstrap

Run a moving-block bootstrap over dates within the evaluation period, with block length

$$
\ell_{blk}\ge H_{\max}.
$$

Recompute on every bootstrap draw:

$$
q_D^{(b)},
\quad
\rho_k^{(b)},
\quad
Q^{(b)},
\quad
r_{eff}^{(b)},
\quad
V_{D,r}^{(b)}.
$$

### §83.2 Enrichment persistence

Define:

$$
\boxed{
\pi_{q_D}
=\Pr_b(q_D^{(b)}>1).
}
$$

### §83.3 Coupling persistence

For the top canonical correlation:

$$
\boxed{
\pi_{\rho_1}
=\Pr_b(\rho_1^{(b)}>\rho_{\min}^{can}).
}
$$

Also report the median and the 10th/90th percentiles.

### §83.4 Directional-mode stability

If the top $r$ canonical D modes are retained,

$$
V_{g,r}=[v_{g,1},\dots,v_{g,r}],
$$

then for bootstrap draw $b$:

$$
S_{g,r}^{(b)} =
\frac1r
\|V_{g,r}^{\top}V_{g,r}^{(b)}\|_F^2.
$$

Range $[0,1]$; $1$ means the canonical directional subspace is completely stable.

Looking at only the first mode:

$$
M_{g,1}^{(b)}=|v_{g,1}^\top v_{g,1}^{(b)}|.
$$

Absolute value is used because a singular vector's sign is not identifiable.

---

## §84. The train / freeze / OOS canonical test

### Motivation

A bootstrap still resamples within the same evaluation-period distribution. The test that actually matters is whether the canonical channel survives on a future period that genuinely did not participate in discovery.

### §84.1 Freeze

On the train window, obtain:

$$
\mathcal S_g,
\quad
G_{VV}^{train},
\quad
G_{DD}^{train},
\quad
u_{g,k}^{train},
\quad
v_{g,k}^{train},
\quad
\rho_{g,k}^{train}.
$$

For the $k$-th canonical pair, freeze the normalized coefficients:

$$
a_{g,k}^{train}
=(G_{VV}^{train})^{-1/2}u_{g,k}^{train},
$$

$$
b_{g,k}^{train}
=(G_{DD}^{train})^{-1/2}v_{g,k}^{train}.
$$

### §84.2 OOS evaluation

On OOS data, only transform the reference $\Phi$ — compute the frozen signals:

$$
v_{g,k,t}^{OOS}
=\Phi_{\mathsf V,t}a_{g,k}^{train},
$$

$$
d_{g,k,t}^{OOS}
=\Phi_{\mathsf D,t}b_{g,k}^{train}.
$$

then compute directly:

$$
\boxed{
\rho_{g,k}^{OOS} =
\operatorname{Corr}_{W_g^{OOS}}
(v_{g,k}^{OOS},d_{g,k}^{OOS}).
}
$$

You may not rerun CCA on OOS data and then report its own top singular value — that would reintroduce selection.

### §84.3 Decay

$$
\boxed{
D_{\rho,g,k} =
\frac{|\rho_{g,k}^{OOS}|}
{|\rho_{g,k}^{train}|+\epsilon}.
}
$$

and sign consistency:

$$
\boxed{
S_{\rho,g,k}
=\operatorname{sign}(\rho_{train})
\operatorname{sign}(\rho_{OOS}).
}
$$

Because a CCA singular value carries no sign of its own, the implementation must fix the canonical pair's orientation on train, then apply that same orientation to the OOS bilinear covariance.

---

## §85. Redundancy and retention across multiple universes

### Motivation

Several candidate universes can all pass eligibility simultaneously while really just being different cuts of the same tickers, or the same D mode. Membership redundancy and mechanism redundancy must be kept separate.

### §85.1 Membership overlap

Jaccard index:

$$
\boxed{
J_{gh} =
\frac{|\mathcal S_g\cap\mathcal S_h|}
{|\mathcal S_g\cup\mathcal S_h|}.
}
$$

For dynamic universes, report the median-date Jaccard along with p10/p90.

### §85.2 Canonical directional similarity

If two universes each retain their top $r_g, r_h$ stable D canonical modes,

$$
V_g\in\mathbb R^{d_D\times r_g},
\qquad
V_h\in\mathbb R^{d_D\times r_h},
$$

the principal-angle similarity is

$$
\boxed{
S_{gh}^{D} =
\frac1{r_*}
\|V_g^\top V_h\|_F^2,
\qquad
r_*=\min(r_g,r_h).
}
$$

If only the top mode is kept in each:

$$
S_{gh}^{D}=|v_{g,1}^\top v_{h,1}|^2.
$$

### §85.3 Retention rule

A universe is not automatically dropped just for high membership overlap. It is treated as a redundant version within the same universe family only when **both**:

$$
J_{gh}>\tau_J
$$

**and**

$$
S_{gh}^{D}>\tau_D
$$

hold. Within the same family, priority for retention goes, in order:

1. stronger evaluation-period adjusted significance;
2. smaller $\rho^{OOS}$ decay;
3. lower membership turnover;
4. if still tied, the larger universe.

This is a lexicographic rule — no arbitrary weighted composite score is built.

---

## §86. Local refit of an accepted universe

### Motivation

The shared $\Phi_{\mathcal U}$'s role is discovery and cross-universe comparison. Once $\mathcal S_g$ is formally accepted, the downstream feature-label-alpha pipeline should treat it as its own market domain, and rebuild a local basis for it.

### §86.1 The local label basis

On universe $g$'s own model-training period:

$$
Y_g
\rightarrow
\Sigma_{Y,g}
\rightarrow
B_{Y,g}
\rightarrow
Z_g
\rightarrow
\Phi_g.
$$

satisfying

$$
\boxed{
\Phi_g^\top W_g\Phi_g=I.
}
$$

From here on, universe $g$'s feature geometry uses

$$
C_g=X_g^\top W_g\Phi_g,
$$

not a permanent dependence on the reference $\Phi_{\mathcal U}$.

### §86.2 Reference-to-local diagnostic

Restrict the reference basis to the subset and orthonormalize it under $W_g$:

$$
Q_g
=\operatorname{orth}_{W_g}(\Phi_{\mathcal U}|_{\mathcal S_g}).
$$

Compare against the local basis:

$$
A_g=Q_g^\top W_g\Phi_g.
$$

Take its SVD:

$$
A_g=U\Sigma V^\top,
$$

$$
\sigma_k=\cos\theta_k.
$$

$\theta_k$ is the reference-to-local principal angle. It's a diagnostic only, not a universe-selection objective.

### §86.3 Interpretation

- small angle: the subset's local mechanism space is roughly a reweighting of the full-universe mechanism;
- large angle: the subset has meaningfully different local label geometry;
- either outcome is acceptable to retain — it just means downstream feature interpretation must use the local $\Phi_g$, not the reference basis.

---

## §87. From canonical D modes to feature research

### Motivation

Up to this point, universe geometry has not used any engineered feature at all. An accepted universe's canonical D modes can serve as an explicit target for subsequent feature research, but they cannot replace full mechanism analysis.

### §87.1 The reference canonical target

For universe $g$'s stable top D canonical coefficient $b_{g,1}$:

$$
\boxed{
d_g^*
=\Phi_{\mathsf D}b_{g,1}.
}
$$

This represents the strongest V-coupled directional combination in this universe.

### §87.2 Feature alignment

For a point-in-time feature $x_p$:

$$
\boxed{
A_{p,g}^{can}
=x_p^\top W_gd_g^*.
}
$$

This must go through the same permutation significance, temporal stability, and OOS decay checks as any other target alignment elsewhere on this site.

### §87.3 The complete downstream pipeline

The canonical target only answers one additional question:

> can a feature predict the low-dimensional D channel that universe discovery found?

Full feature analysis still has to compute

$$
c_{p,g}=x_p^\top W_g\Phi_g,
$$

$$
Q_{p,g}^{D},Q_{p,g}^{V},Q_{p,g}^{L},
$$

fingerprint, stability, breadth, redundancy, coverage, and alpha attribution — the canonical target is a useful extra lens, not a shortcut around the rest of [feature evaluation](/quant/feature-evaluation/).
