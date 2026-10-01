---
title: "Cross-Sectional Predictive Geometry, Part 6: Target Alignment and Directional Efficiency"
description: "Having directional information isn't the same as having it aimed at the horizon combination you actually trade. Defining the intended target direction, the signed alignment that tells you which way and how hard, and the efficiency ratio that separates 'how much direction' from 'how well aimed'."
date: 2026-09-14
series: ["Cross-Sectional Predictive Geometry"]
part: 6
tags: ["cross-sectional", "methodology", "quant-research", "predictive-geometry"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Predictive Geometry" >}}

# Part Four: Target Direction, Similarity, and Stability

## §20. Intended Target Direction

### Motivation

$Q^{\mathsf D}_p$ tells us how much directional information a feature has,
but not whether **that direction is the one we want**. In a $K=3$ band, the
$\mathsf D$ block is three-dimensional; a feature could be very strong along
some direction orthogonal to the horizon combination we actually trade.
That's perfectly valid mathematically and entirely useless for trading.

So "the direction we want" has to be defined explicitly.

### Definition

The user specifies horizon weights $a\in\mathbb R^K$ (in **raw return
units**, before whitening), defining the target signal

$$\boxed{y^*_d=Y_d^{(\beta)}a\ \in\mathbb R^{N_d}.}$$

For example, in the $H_2=\{1h,2h,4h\}$ incremental band, $a=(1,1,1)^\top$
means "I care about the total return held to 4 hours"; $a=(1,0,0)^\top$
means "I only care about the first hour."

The target direction in $\mathsf D$-block coordinates is defined by
projecting $y^*$ onto $\Phi_{\mathsf D}$ and normalizing:

$$\boxed{r_{\mathsf D}=\frac{\Phi_{\mathsf D}^\top Wy^*}{\big\|\Phi_{\mathsf D}^\top Wy^*\big\|},\qquad\|r_{\mathsf D}\|=1.}$$

**This is the implementation spec.** It's a one-line computation, and it's
immune to any rank truncation in §7 or any coordinate convention — because
it measures directly "where the target signal itself lands in
$\Phi_{\mathsf D}$," rather than relying on any closed-form inversion.

### Proposition 20.1 (closed-form equivalence)

**No-truncation case ($K_{\rm eff}=K$)**: $\Phi_{\mathsf
D}=Z=Y^{(\beta)}\Sigma_Y^{-1/2}$, so

$$\Phi_{\mathsf D}^\top Wy^*=\Sigma_Y^{-1/2}\big(Y^{\top}WY\big)a=\Sigma_Y^{-1/2}\Sigma_Ya=\Sigma_Y^{1/2}a,$$

$$\boxed{r_{\mathsf D}=\frac{\Sigma_Y^{1/2}a}{\sqrt{a^\top\Sigma_Ya}}.}$$

**Truncation case**: $B_Y=E_{1:K_{\rm eff}}\Lambda_{1:K_{\rm eff}}^{-1/2}$,
so

$$\Phi_{\mathsf D}^\top Wy^*=B_Y^\top\Sigma_Ya=\Lambda^{-1/2}E^\top\Sigma_Ya=\Lambda^{1/2}E^\top a\quad(\text{first }K_{\rm eff}\text{ components}),$$

$$r_{\mathsf D}=\frac{\Lambda_{1:K_{\rm eff}}^{1/2}E_{1:K_{\rm eff}}^\top a}{\big\|\Lambda_{1:K_{\rm eff}}^{1/2}E_{1:K_{\rm eff}}^\top a\big\|}.$$

$\blacksquare$

The closed form is for understanding only: it shows $r_{\mathsf D}$'s
direction is "reweighted" by $\Sigma_Y$'s spectrum — higher-variance return
directions get proportionally more weight in the target. This is correct:
in whitened coordinates, achieving one unit of $y^*$ requires more
contribution from the high-variance directions.

**In implementation, always use the projection formula.** This eliminates
an entire class of the most hard-to-spot bugs: any rotation, truncation, or
reordering of the coordinate system is automatically followed by the
projection formula, and is not by the closed form.

### Worked Example

$H_2=\{1h,2h,4h\}$ incremental band, $\Sigma_Y=\operatorname{diag}(\sigma_1^2,
\sigma_2^2,\sigma_3^2)$, $\sigma=(0.004,0.005,0.007)$ (per-period standard
deviations). Take $a=(1,1,1)$:

$$\Sigma_Y^{1/2}a=(0.004,0.005,0.007),\qquad \|\cdot\|=0.00949,\qquad r_{\mathsf D}=(0.421,0.527,0.738).$$

The target direction loads most heavily on the longest horizon, because it
contributes the most return variance. Taking $a=(1,0,0)$ instead gives
$r_{\mathsf D}=(1,0,0)$.

---

## §21. Signed Target Alignment

### Definition

$$\boxed{A_p^{\rm tgt}=(c_p^{\mathsf D})^\top r_{\mathsf D}.}$$

This is a **signed** quantity. Its meaning: feature $p$'s average per-period
cross-sectional IC against the target signal $y^*$ (under $\Phi_{\mathsf
D}$'s metric).

### Interpretation

- $A^{\rm tgt}_p>0$: higher feature value implies higher target return.
- $A^{\rm tgt}_p<0$: an inverted signal. This is **not a defect** — as long
  as the sign is stable over time (§27's $\pi^{\rm time,sign}$), an inverted
  signal is just as usable as a direct one; simply negate it when
  constructing positions.
- $A^{\rm tgt}_p\approx0$ but $Q^{\mathsf D}_p$ large: the feature has
  directional information, but it's aimed at a direction orthogonal to the
  target. Flag it `TARGET_MISMATCH` and re-evaluate it under a different
  band or a different $a$.

### Pitfall: the Sign Must Be Frozen

Once $A^{\rm tgt}$'s sign is determined in the research period, it must be
treated as part of the feature's definition and locked. If the sign flips
during OOS and gets "conveniently" flipped back, that's a decision made
using OOS data, invalidating the whole OOS evaluation. A sign flip is a
**result**, not a bug to be patched (§28).

### Significance

$(A^{\rm tgt}_p)^2$ is permutation-tested and Westfall–Young-adjusted
exactly like $Q$, giving $\tilde p^{\rm tgt}_p$. This quantity plays a key
role in §34's eligibility, because it merges "how much directional
information" and "how well aimed" into a single testable quantity — it's
precisely "the average squared IC against the intended target."

---

## §22. Directional Target Efficiency

### Definition

$$\boxed{\eta_p^{\mathsf D}=\frac{(A_p^{\rm tgt})^2}{Q_p^{\mathsf D}+\epsilon}.}$$

### Proposition 22.1 (bound)

$$0\le\eta^{\mathsf D}_p\le1.$$

**Proof**: by Cauchy–Schwarz, $(c_p^{\mathsf D\top}r_{\mathsf D})^2\le
\|c_p^{\mathsf D}\|^2\|r_{\mathsf D}\|^2=Q_p^{\mathsf D}$. $\blacksquare$

### Geometric Interpretation

$\eta^{\mathsf D}_p=\cos^2\theta$, where $\theta$ is the angle between
$c_p^{\mathsf D}$ and $r_{\mathsf D}$.

- $\eta\approx1$: the feature's directional information is **almost
  entirely** aimed at the target.
- $\eta\approx0.3$: only 30% of the directional energy is on the target
  direction; the rest points at other horizon combinations within the band.
- $\eta\approx0$: completely misaligned.

### How Three Quantities Divide Labor

$$Q^{\mathsf D}_p\ \text{(how much directional information)},\quad A^{\rm tgt}_p\ \text{(which way, how precisely)},\quad \eta^{\mathsf D}_p\ \text{(precision ratio)}.$$

The three satisfy $(A^{\rm tgt}_p)^2=Q^{\mathsf D}_p\cdot\eta^{\mathsf
D}_p$, so there are only two degrees of freedom. All three are reported
because they answer different practical questions: $Q^{\mathsf D}$ decides
"is this feature worth keeping at all," $\eta^{\mathsf D}$ decides "is it
worth keeping for the *current* target," and $A^{\rm tgt}$ decides "how to
use it."

### Worked Example

$$Q^{\mathsf D}_p=4.0\times10^{-4},\quad A^{\rm tgt}_p=-1.6\times10^{-2},\quad \eta^{\mathsf D}_p=\frac{2.56\times10^{-4}}{4.0\times10^{-4}}=0.64.$$

Reading: this is an inverted directional feature, with average IC of
$-0.016$ against the target, and 64% of its directional energy aimed at the
target. This is a fully usable signal: feed $-X_p$ into the long/short
ranking.

Another feature:

$$Q^{\mathsf D}_q=9.0\times10^{-4},\quad A^{\rm tgt}_q=3.0\times10^{-3},\quad \eta^{\mathsf D}_q=0.01.$$

$Q^{\mathsf D}$ is larger, but $\eta=0.01$: almost all of its directional
information points at a direction orthogonal to the target. If $a=(1,1,1)$,
this feature might actually be predicting "1h up but 4h down" — a
horizon-hedge structure. Useless for the current target, but potentially
very strong in a study with $a=(1,0,-1)$. Flag `TARGET_MISMATCH`, keep in
the library.

### The Role of $\epsilon$

The denominator's $\epsilon$ is purely numerical protection. When
$Q^{\mathsf D}_p$ is itself at noise level (failed §18's significance
test), $\eta$ is a ratio of two noise quantities, numerically anywhere
between 0 and 1, **meaningless**. $\eta$ should therefore only be read for
features that already passed the significance test; the report should show
a blank, not a number, for the rest.
