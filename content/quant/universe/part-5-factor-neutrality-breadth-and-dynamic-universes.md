---
title: "Multi-Sub-Universe Directional Coupling Geometry, Part 5: Factor-Neutrality, Breadth, and Dynamic Universes"
description: "Why demeaning doesn't remove heterogeneous beta dispersion and when a factor-neutral variant is needed, how to check whether a universe's coupling is driven by a handful of names, and how to update membership over time without leaking the future into the selection rule."
date: 2026-09-14
series: ["Multi-Sub-Universe Directional Coupling Geometry"]
part: 5
tags: ["universe", "cross-sectional", "methodology", "quant-research"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Multi-Sub-Universe Directional Coupling Geometry" >}}

## §88. Beta, demeaning, and the factor-neutral variant

### Baseline

The existing label preprocessing is within-date demeaning:

$$
Y^{XS}_{id,h}
=Y^{raw}_{id,h}-\bar Y_{d,h}.
$$

This is equivalent to projecting out

$$
\operatorname{span}\{\mathbf1\}.
$$

So it's a rank-one special case of unit-exposure common-factor neutralization.

### It does not remove heterogeneous beta dispersion

If

$$
r_i=\beta_i r_m+\epsilon_i,
$$

then after demeaning,

$$
r_i-\bar r
=(\beta_i-\bar\beta)r_m
+(\epsilon_i-\bar\epsilon).
$$

So a subset with high beta dispersion can inflate $q_D$ or produce V-D coupling purely through that mechanism.

### The factor-neutral variant

If the goal is to study idiosyncratic geometry specifically, first apply

$$
Y_d^{idio}
=M_{F_d}Y_d,
$$

$$
M_{F_d}
=I-F_d(F_d^\top F_d)^{-1}F_d^\top,
$$

then rebuild the entire chain:

$$
\Phi_{\mathcal U}^{idio},
\quad
G_g^{idio},
\quad
R_{VD,g}^{idio}.
$$

Mixing blocks from the raw basis and the idiosyncratic basis in the same calculation is not permitted.

### Reporting

When a factor-neutral variant exists, report both sets side by side:

$$
q_D^{raw},Q_{VD}^{raw},\rho_1^{raw},
$$

$$
q_D^{idio},Q_{VD}^{idio},\rho_1^{idio}.
$$

Strong raw / weak idio means the structure is driven mainly by known common exposures. Both strong means the structure survives factor neutralization.

---

## §89. Cross-sectional breadth and single-ticker dominance

### Motivation

A universe's high $q_D$ or high $\rho_1$ might be determined by a tiny handful of names. Cross-sectional breadth has to be checked explicitly.

### §89.1 D-energy contribution

For date $d$:

$$
e_{id}^{D}
=\|\phi_{id}^{D}\|^2.
$$

Normalized share:

$$
a_{id}^{D}
=\frac{e_{id}^{D}}
{\sum_{j\in\mathcal S_{g,d}}e_{jd}^{D}}.
$$

Effective contributor count:

$$
\boxed{
N_{g,d}^{eff,D}
=\frac1{\sum_i(a_{id}^{D})^2}.
}
$$

Relative breadth:

$$
\boxed{
b_{g,d}^{D}
=\frac{N_{g,d}^{eff,D}}{N_{g,d}}.
}
$$

### §89.2 Canonical-channel breadth

For a frozen canonical pair $(a_g,b_g)$, each observation's coupling contribution is

$$
z_{id}^{can}
=(\phi_{id}^{V\top}a_g)
(\phi_{id}^{D\top}b_g).
$$

Define share by absolute contribution:

$$
c_{id}^{can}
=\frac{|z_{id}^{can}|}
{\sum_j|z_{jd}^{can}|}.
$$

$$
\boxed{
N_{g,d}^{eff,can}
=\frac1{\sum_i(c_{id}^{can})^2}.
}
$$

A high $\rho_1$ paired with very low canonical breadth should be flagged as concentrated — not automatically excluded, but flagged, since it changes what the universe's tradability and capacity actually look like.

---

## §90. The dynamic universe

### Motivation

If universe membership updates over time, the selection rule must depend only on label history already realized strictly before that point in time.

### §90.1 Rolling ticker operators

For ticker $i$:

$$
\bar K_{i,t}^{D}
=
\sum_{\tau<t}\omega_{t-\tau}
K_{i\tau}^{D},
$$

$$
\bar K_{i,t}^{VD}
=
\sum_{\tau<t}\omega_{t-\tau}
K_{i\tau}^{VD}.
$$

$\omega$ can be a fixed rolling window or exponential decay, and must be locked in advance.

### §90.2 Membership update

$$
\mathcal S_{g,t}
=\mathcal A_g
(\{\bar K_{i,t}^{D},\bar K_{i,t}^{VD}\}_i).
$$

Update frequency should default to weekly or monthly — not updating membership on every 5-minute bar.

### §90.3 Turnover

Universe turnover:

$$
\boxed{
\operatorname{Turn}_g(t)
=1-
\frac{|\mathcal S_{g,t}\cap\mathcal S_{g,t-1}|}
{|\mathcal S_{g,t}\cup\mathcal S_{g,t-1}|}.
}
$$

Report median and p90 turnover. High turnover means the universe's identity is itself unstable, and downstream feature training will be affected by distribution drift.

### §90.4 Freeze boundary

The dynamic rule's functional form, window, cluster count, size grid, and thresholds are all frozen once discovery is complete. The evaluation period can update the *statistics*, but it cannot modify the *rule*.
