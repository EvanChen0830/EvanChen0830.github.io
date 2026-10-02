---
title: "Multi-Sub-Universe Directional Coupling Geometry, Part 1: Research Goals and the Reference Basis"
description: "Universe discovery as its own geometric research layer, before feature engineering ever starts: what a sub-universe family is, the three separate claims a good sub-universe must satisfy, the discovery/evaluation boundary, and the shared reference basis that makes different candidate universes comparable."
date: 2026-09-14
series: ["Multi-Sub-Universe Directional Coupling Geometry"]
part: 1
tags: ["universe", "cross-sectional", "methodology", "quant-research"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Multi-Sub-Universe Directional Coupling Geometry" >}}

## The problem this chapter solves

The [Cross-Sectional Predictive Geometry](/quant/feature-evaluation/) framework starts from a *fixed* universe: build a label mechanism basis from future returns,

$$
\Phi=[\Phi_{\mathsf D}\mid\Phi_{\mathsf V}\mid\Phi_{\mathsf L}],
\qquad
\Phi^\top W\Phi=I,
$$

then project features, alphas, and residuals into that one coordinate system. That setup answers "given this universe, which features have predictive power, and what did the model learn" — but it does not answer a more upstream question:

> **Which group of tickers should features and models even be built on? Do multiple sub-universes exist in which directional label structure is more concentrated, lower-dimensional, and forms a stable low-rank coupling with activity structure than it does in the parent universe?**

This chapter promotes the universe itself to an analyzable geometric object. For any candidate sub-universe $\mathcal S_g$, membership is not treated as an ordinary feature — it's treated as a restriction of the observation measure. That gives

$$
G_g=\Phi_{\mathcal U}^\top W_g\Phi_{\mathcal U},
$$

where $\Phi_{\mathcal U}$ is the shared reference basis fit on the parent universe $\mathcal U$, and $W_g$ is the row weight renormalized to live only on $\mathcal S_g$.

$G_g$ is the Gram operator of the reference label basis, restricted to the sub-universe. It answers: among the label mechanisms that are mutually orthogonal and equally scaled in the parent universe, which get amplified, which get suppressed, and which start coupling to each other once you restrict to $\mathcal S_g$?

The core object of this chapter isn't just directional energy $G_{g,\mathsf{DD}}$ on its own — it's the canonical coupling between activity and direction:

$$
R_{g,\mathsf{VD}} =
G_{g,\mathsf{VV}}^{-1/2}
G_{g,\mathsf{VD}}
G_{g,\mathsf{DD}}^{-1/2},
$$

and the directional explained-variance operator built from it,

$$
\boxed{
H^{\mathsf D\leftarrow\mathsf V}_g =
R_{g,\mathsf{VD}}^\top R_{g,\mathsf{VD}} =
G_{g,\mathsf{DD}}^{-1/2}
G_{g,\mathsf{DV}}
G_{g,\mathsf{VV}}^{-1}
G_{g,\mathsf{VD}}
G_{g,\mathsf{DD}}^{-1/2}.
}
$$

Its eigenvalues are squared canonical correlations; its eigenvectors are "the directional modes best represented by activity geometry in this sub-universe."

> **A note on notation.** This book, like the [Cross-Sectional Predictive Geometry](/quant/feature-evaluation/) series, decomposes predictive structure into three blocks: $\mathsf D$ (Direction), $\mathsf V$ (Activity), and $\mathsf L$ — the same block the companion series calls $T$ (Tail), $G_{\mathsf L}(z)=\|z\|z$. Same definition, different letter. Here that machinery is applied one level up: not to whether a *feature* predicts direction, but to whether a *sub-universe's membership* concentrates and structures the direction block itself.

So this chapter does not output a single best universe. It outputs a family of possibly-overlapping universes,

$$
\boxed{
\mathscr S =
\{\mathcal S_1,\dots,\mathcal S_G\},
}
$$

and reports four separate facets for each one:

1. **Directional enrichment** — is $\mathsf D$ energy higher than in the parent universe?
2. **Directional compressibility** — is $\mathsf D$ energy concentrated in a few directional modes?
3. **$\mathsf V \leftrightarrow \mathsf D$ coupling** — do activity and direction form strong, low-rank, repeatable canonical channels?
4. **Out-of-sample stability** — does the above geometry still hold on future data that didn't participate in discovery?

All four are purely **label-side structure**. Actual learnability still has to be confirmed afterward, on point-in-time features and alpha OOS results, once the universe has been frozen.

---

## §68. Research goal and output objects

### Motivation

Let the parent universe be $\mathcal U_d$, allowed to change with date $d$. The conventional approach usually cuts the universe first — by industry, market cap, liquidity, or hand-written rules — and only then builds features inside each slice. That treats the universe as an exogenous setting, and can't answer:

- which tickers have directional label geometry that is inherently more concentrated;
- which sub-universes have activity and direction forming a low-dimensional coupling;
- whether several distinct market domains exist, each with its own clean directional representation;
- whether a universe's label-side structure is merely a byproduct of beta dispersion, volatility scale, or sample mining.

So universe discovery is defined explicitly as an independent research layer, prior to feature engineering.

### Definition 68.1 (Parent universe)

For every observation time $d = 1, \dots, T$, the parent universe is

$$
\mathcal U_d,
\qquad
N_d=|\mathcal U_d|.
$$

Every candidate universe must satisfy

$$
\mathcal S_{g,d}\subseteq\mathcal U_d.
$$

### Definition 68.2 (Universe family)

Up to $G$ sub-universes may be retained simultaneously:

$$
\boxed{
\mathscr S =
\{\mathcal S_1,\dots,\mathcal S_G\}.
}
$$

Mutual exclusivity is not required:

$$
\mathcal S_g\cap\mathcal S_h\neq\varnothing
$$

is allowed; nested families are also allowed:

$$
\mathcal S_1\subset\mathcal S_2\subset\cdots.
$$

The goal of this chapter is not to solve

$$
\arg\max_{\mathcal S}J(\mathcal S)
$$

for a unique winner — it's to retain, from a pre-specified candidate generator, the set of universes that clear statistical and stability tests and carry non-redundant structure in their directional mechanism.

### Three separate levels

A "good" sub-universe splits into three distinct claims:

**U1: D-enriched.**

$$
q_g^{\mathsf D}>1.
$$

Average label energy per directional dimension in this subset is higher than in the parent universe.

**U2: D-compressible.**

$$
r_{g,\rm eff}^{\mathsf D}\ll d_{\mathsf D}.
$$

Directional energy is concentrated in a small number of horizon combinations.

**U3: D-structurally represented.**

$$
Q_g^{\mathsf D\leftarrow\mathsf V}
\text{ large},
\qquad
r_{g,\rm eff}^{\mathsf D\leftarrow\mathsf V}
\text{ small}.
$$

That is: a small number of activity modes form a strong canonical coupling with a small number of directional modes.

These three are not equivalent to each other, and a report must keep them separate.

---

## §69. The point-in-time and discovery/evaluation boundary

### Principle

Universe discovery uses future labels, so it is itself a supervised selection procedure. If, over the same period, you simultaneously:

1. use the label to choose membership;
2. use that same label to evaluate that membership's geometry;
3. treat the result as statistical evidence;

the result is necessarily selection bias.

### Definition 69.1 (Time split)

Split the universe-research data into

$$
\mathcal D_U
\quad\text{(discovery)},
$$

$$
\mathcal E_U
\quad\text{(evaluation)},
$$

with a purge gap between them,

$$
\Delta_U\ge H_{\max}.
$$

If the candidate generator uses a rolling history of length $\ell_U$, then

$$
\boxed{
\Delta_U\ge H_{\max}+\ell_U^{\rm overlap},
}
$$

where $\ell_U^{\rm overlap}$ is any additional lookback overlap that would cross the evaluation boundary.

### Rule 69.1 (Freeze)

Decide, on $\mathcal D_U$ only:

- the candidate-generation rule;
- membership $\mathcal S_g$ or the membership rule;
- the reference basis $\Phi_{\mathcal U}$;
- the retained canonical rank;
- the universe-family pruning rule;
- every threshold.

Once $\mathcal E_U$ begins, all of the above are frozen. Membership or thresholds may not be re-adjusted because of evaluation-period results.

### Static vs. dynamic universes

**Static research universe**: membership is decided only by the discovery period, and held fixed across the entire evaluation period.

**Dynamic research universe**: at every evaluation time $t$, membership may only be decided by rolling geometry already realized before $t$:

$$
\mathcal S_{g,t} =
\mathcal A_g(\mathcal F_{t^-}),
$$

where $\mathcal F_{t^-}$ denotes information available strictly before time $t$.

Any rule that uses $Y_{t\to t+h}$ directly to decide $\mathcal S_{g,t}$ constitutes leakage.

---

## §70. The full-universe reference basis

### Motivation

If every candidate universe independently refit its own whitening and $\Phi_g$ from scratch, each universe's coordinate system would differ, and $q_D$, coupling directions, and eigenvectors could not be compared directly across universes. So the discovery stage first establishes one shared reference coordinate system.

### Definition

On the parent universe, during the discovery period, fit exactly the existing framework end to end:

$$
Y^{XS}_{\mathcal U}
\longrightarrow
Z_{\mathcal U}
\longrightarrow
\Phi_{\mathcal U} =
[\Phi_{\mathsf D}\mid\Phi_{\mathsf V}\mid\Phi_{\mathsf L}].
$$

satisfying

$$
\boxed{
\Phi_{\mathcal U}^\top
W_{\mathcal U}
\Phi_{\mathcal U}
=I_{d_\Phi}.
}
$$

where

$$
W_{\mathcal U} =
\operatorname{diag}
\left(
\frac1{TN_d}
\right).
$$

### Reading

This identity fixes the parent universe as the geometric benchmark:

- every reference mechanism direction has energy exactly 1;
- different reference directions are mutually orthogonal;
- the three blocks are mutually orthogonal in the parent universe.

So any deviation observed on a candidate subset can be read directly as concentration, suppression, or coupling caused by the restriction itself — not as an artifact of a different coordinate system.

### Fit / transform

$B_Y$, $\Theta_{\mathsf V}$, $\Theta_{\mathsf L}$, and $\Gamma_\bullet$ are all fit on $\mathcal D_U$ only. Evaluation-period labels are only ever transformed into the already-frozen reference basis — never refit.
