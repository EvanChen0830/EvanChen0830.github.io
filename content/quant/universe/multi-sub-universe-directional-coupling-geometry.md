---
title: "Multi-Sub-Universe Directional Coupling Geometry"
description: "Treating the universe itself as an analyzable geometric object: restriction Gram operators, V↔D coupling, candidate generation, stratified permutation nulls, OOS freeze tests, breadth, and local refit."
date: 2026-09-14
tags: ["universe", "geometry", "methodology", "quant-research", "coupling-operator"]
math: true
aliases:
  - "/quant/universe/part-1-research-goals-and-reference-basis/"
  - "/quant/universe/part-2-restriction-measures-and-the-coupling-operator/"
  - "/quant/universe/part-3-candidate-generation-and-the-universe-null/"
  - "/quant/universe/part-4-stability-oos-testing-and-redundancy/"
  - "/quant/universe/part-5-factor-neutrality-breadth-and-dynamic-universes/"
  - "/quant/universe/part-6-reports-implementation-order-and-validation/"
  - "/quant/universe/part-7-unified-viewpoint-and-quick-reference/"
---

## The problem this chapter solves

The [Cross-Sectional Predictive Geometry](/quant/feature-evaluation/) framework starts from a *fixed* universe: build a label mechanism basis from future returns,

$$
\Phi=[\Phi_{\mathsf D}\mid\Phi_{\mathsf V}\mid\Phi_{\mathsf T}],
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

> **A note on notation.** This book, like the [Cross-Sectional Predictive Geometry](/quant/feature-evaluation/) series, decomposes predictive structure into three blocks: $\mathsf D$ (Direction), $\mathsf V$ (Activity), and $\mathsf T$ (Tail), where $G_{\mathsf T}(z)=\|z\|_2\,z$. Here that machinery is applied one level up: not to whether a *feature* predicts direction, but to whether a *sub-universe's membership* concentrates and structures the direction block itself.

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

## 1. Research goal and output objects

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

## 2. The point-in-time and discovery/evaluation boundary

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

## 3. The full-universe reference basis

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
[\Phi_{\mathsf D}\mid\Phi_{\mathsf V}\mid\Phi_{\mathsf T}].
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

$B_Y$, $\Theta_{\mathsf V}$, $\Theta_{\mathsf T}$, and $\Gamma_\bullet$ are all fit on $\mathcal D_U$ only. Evaluation-period labels are only ever transformed into the already-frozen reference basis — never refit.

---

## 4. Sub-universe weight and the restriction measure

### 4.1 Definition

For a candidate universe $\mathcal S_g$, let

$$
N_{g,d}=|\mathcal S_{g,d}|.
$$

Define the renormalized row weight:

$$
\boxed{
w^{(g)}_{id} =
\frac{\mathbf1\{i\in\mathcal S_{g,d}\}}
{T_gN_{g,d}},
}
$$

where $T_g$ is the number of observation times for which this universe is effective.

In implementation no diagonal matrix is actually built; every quantity is computed as

$$
\frac1{T_g}
\sum_{d\in\mathcal D_g}
\frac1{N_{g,d}}
(\cdot).
$$

### Why renormalize

If you only mask,

$$
\widetilde W_g=M_gW_{\mathcal U},
$$

then every energy shrinks linearly with subset size — conflating two different things: "how many tickers are in the subset" and "what the subset's own label geometry looks like."

The purpose of $W_g$ here is to ask:

> **If $\mathcal S_g$ is treated as a complete cross-section in its own right, how does its average label geometry change relative to the parent universe?**

So every date is still equally weighted, and every ticker within the subset is equally weighted.

### Minimum size

If $N_{g,d}\lt d_\Phi+5$, the Gram matrices may be numerically computable but statistically very unstable. This chapter requires at minimum:

$$
\operatorname{median}_dN_{g,d}\ge N_{\min},
$$

where $N_{\min}$ should be much larger than $d_\Phi$; see Section 29 for the default.

---

## 5. The restriction Gram operator

### Definition 72.1

For any candidate universe $g$:

$$
\boxed{
G_g =
\Phi_{\mathcal U}^\top
W_g
\Phi_{\mathcal U}
\in\mathbb R^{d_\Phi\times d_\Phi}.
}
$$

### Proposition 72.1 (Positive semi-definite)

$$
G_g\succeq0.
$$

**Proof.** For any $v\in\mathbb R^{d_\Phi}$,

$$
v^\top G_gv =
(\Phi_{\mathcal U}v)^\top
W_g
(\Phi_{\mathcal U}v)
\ge0.
$$

$\blacksquare$

### 5.1 Physical meaning of the Rayleigh quotient

Take any unit mechanism direction $\|v\|=1$. In the parent universe:

$$
\|\Phi_{\mathcal U}v\|^2_{W_{\mathcal U}}
=v^\top Iv=1.
$$

On the subset:

$$
\boxed{
\|\Phi_{\mathcal U}v\|^2_{W_g}
=v^\top G_gv.
}
$$

So $v^\top G_gv$ is exactly the **relative energy** of that full-universe mechanism combination, measured inside $\mathcal S_g$.

If

$$
G_gv_k=\lambda_{g,k}v_k,
$$

then:

- $\lambda_{g,k}>1$: this direction is concentrated in the subset relative to the parent universe;
- $\lambda_{g,k}<1$: it's suppressed;
- $\lambda_{g,k}=1$: average energy is unchanged.

### 5.2 Not a learnability operator

$G_g$ is built entirely from future labels. A large eigenvalue only means label variation / mechanism energy is concentrated — it says nothing about whether any point-in-time $X$ exists that could actually predict that direction.

So $G_g$ is a **label-geometry operator**, not a predictive operator.

---

## 6. Block decomposition and directional enrichment

By the block structure of the reference basis:

$$
G_g=
\begin{pmatrix}
G_{g,\mathsf{DD}} & G_{g,\mathsf{DV}} & G_{g,\mathsf{DL}}\\
G_{g,\mathsf{VD}} & G_{g,\mathsf{VV}} & G_{g,\mathsf{VL}}\\
G_{g,\mathsf{LD}} & G_{g,\mathsf{LV}} & G_{g,\mathsf{LL}}
\end{pmatrix}.
$$

### 6.1 Average mechanism energy per dimension

Because block dimensions differ, traces can't be compared directly. Define

$$
\boxed{
q_g^{\mathsf B} =
\frac1{d_{\mathsf B}}
\operatorname{tr}(G_{g,\mathsf{BB}}),
\qquad
\mathsf B\in\{\mathsf D,\mathsf V,\mathsf T\}.
}
$$

The parent universe satisfies, exactly:

$$
q_{\mathcal U}^{\mathsf D}
=q_{\mathcal U}^{\mathsf V}
=q_{\mathcal U}^{\mathsf T}
=1.
$$

So $q_g^{\mathsf D}=1.7$ means:

> this sub-universe's average squared amplitude per reference directional dimension is $1.7\times$ the parent universe's.

### 6.2 Universe mechanism fingerprint

Define

$$
\boxed{
f_g^{\mathsf B} =
\frac{q_g^{\mathsf B}}
{q_g^{\mathsf D}+q_g^{\mathsf V}+q_g^{\mathsf T}}.
}
$$

The parent-universe baseline is

$$
f_{\mathcal U}=(1/3,1/3,1/3).
$$

$q_g^{\mathsf D}$ is absolute enrichment; $f_g^{\mathsf D}$ is relative purity. Both must be reported together.

### Example

Universe A:

$$
(q_D,q_V,q_L)=(0.8,0.15,0.15).
$$

Then $f_D\approx0.73$, but $q_D<1$: it's direction-*pure* but *weak*.

Universe B:

$$
(q_D,q_V,q_L)=(1.9,0.9,0.8).
$$

This one has both absolute D enrichment and relative D enrichment at once.

---

## 7. The directional spectrum and compressibility

### Motivation

A large $q_g^{\mathsf D}$ doesn't by itself mean the directional representation is simple. If all $d_{\mathsf D}$ directions are scaled up together proportionally, there is no low-dimensional structure.

### 7.1 The directional eigenspectrum

Eigendecompose

$$
G_{g,\mathsf{DD}}
=V_g^{\mathsf D}
\Lambda_g^{\mathsf D}
V_g^{\mathsf D\top}
$$

to get

$$
\lambda_{g,1}^{\mathsf D}
\ge\cdots\ge
\lambda_{g,d_{\mathsf D}}^{\mathsf D}\ge0.
$$

$v_{g,k}^{\mathsf D}$ is the $k$-th most concentrated reference directional combination in this universe.

### 7.2 Effective directional rank

Let

$$
p_{g,k}^{\mathsf D} =
\frac{\lambda_{g,k}^{\mathsf D}}
{\sum_j\lambda_{g,j}^{\mathsf D}}.
$$

Define the entropy effective rank:

$$
\boxed{
r_{g,\rm eff}^{\mathsf D} =
\exp\left(
-\sum_kp_{g,k}^{\mathsf D}
\log p_{g,k}^{\mathsf D}
\right).
}
$$

$$
1\le r_{g,\rm eff}^{\mathsf D}\le d_{\mathsf D}.
$$

### 7.3 Top-$m$ directional concentration

$$
\boxed{
C_{g,m}^{\mathsf D} =
\frac{\sum_{k=1}^m\lambda_{g,k}^{\mathsf D}}
{\sum_j\lambda_{g,j}^{\mathsf D}}.
}
$$

and

$$
\boxed{
r_{g,\eta}^{\mathsf D} =
\min\left\{m:C_{g,m}^{\mathsf D}\ge\eta\right\}.
}
$$

### Reading

- high $q_D$ + high $r_{\rm eff}^{D}$: directional energy is strong but spread out;
- high $q_D$ + low $r_{\rm eff}^{D}$: directional energy is strong *and* low-dimensional;
- low $q_D$ + low $r_{\rm eff}^{D}$: low-dimensional but weak;
- low $q_D$ + high $r_{\rm eff}^{D}$: weak and spread out.

This section only answers **whether D representation is concentrated** — it does not yet say whether magnitude/activity structure actually represents that D structure. That's next.

---

## 8. From $|Z|^\top WZ$ to $\mathsf V \leftrightarrow \mathsf D$ coupling

### Motivation

The most intuitive magnitude-direction cross-moment is

$$
A_{MD}=|Z|^\top WZ,
$$

whose entries

$$
(A_{MD})_{hk}
=\mathbb E_W[|Z_h|Z_k]
$$

answer: is horizon $h$'s move magnitude related to horizon $k$'s signed direction?

The intuition is right, but this object has three limitations:

1. $|Z|$ only keeps per-horizon absolute magnitude — it has no complete second-order activity structure;
2. its singular values are simultaneously driven by magnitude variance and direction variance;
3. the existing framework already has a more complete activity block $\Phi_{\mathsf V}$, already orthogonalized against $\mathsf D$.

So the formal method replaces $|Z|$ with $\Phi_{\mathsf V}$, and poses the question as $\mathsf V \leftrightarrow \mathsf D$ coupling within the subset.

### 8.1 The cross-block Gram

Define

$$
\boxed{
G_{g,\mathsf{VD}} =
\Phi_{\mathsf V}^\top
W_g
\Phi_{\mathsf D}.
}
$$

In the parent universe,

$$
\boxed{
G_{\mathcal U,\mathsf{VD}}=0
}
$$

holds exactly, by the block orthogonality of the reference basis.

On a subset it generally is no longer zero:

$$
G_{g,\mathsf{VD}}\neq0.
$$

The physical meaning:

> Activity and direction, already orthogonally separated in the parent universe, exhibit systematic covariation under the restriction measure $\mathcal S_g$.

This is not a basis-construction error — it's a signal about the subset's geometry.

---

## 9. The canonical $\mathsf V \leftrightarrow \mathsf D$ coupling operator

### Motivation

Looking at $G_{g,\mathsf{VD}}$ directly still conflates the coupling with $\mathsf V$'s and $\mathsf D$'s own scale within the subset. Both sides need to be whitened first.

### 9.1 Regularized inverse square root

For

$$
G_{g,\mathsf{VV}}
=E_V\Lambda_VE_V^\top,
$$

$$
G_{g,\mathsf{DD}}
=E_D\Lambda_DE_D^\top,
$$

retain only

$$
\lambda_j\ge\epsilon_{\rm can}\lambda_{\max}.
$$

Define the pseudo-inverse square root:

$$
G_{g,\mathsf{VV}}^{-1/2}
=E_{V,r}\Lambda_{V,r}^{-1/2}E_{V,r}^\top,
$$

$$
G_{g,\mathsf{DD}}^{-1/2}
=E_{D,r}\Lambda_{D,r}^{-1/2}E_{D,r}^\top.
$$

### 9.2 Definition

$$
\boxed{
R_{g,\mathsf{VD}} =
G_{g,\mathsf{VV}}^{-1/2}
G_{g,\mathsf{VD}}
G_{g,\mathsf{DD}}^{-1/2}.
}
$$

Take its SVD:

$$
\boxed{
R_{g,\mathsf{VD}}
=U_g
\operatorname{diag}(\rho_{g,1},\dots,\rho_{g,r_g})
V_g^\top,
}
$$

where

$$
1\ge\rho_{g,1}\ge\cdots\ge0.
$$

### Proposition 76.1 (Equivalence to CCA)

$\rho_{g,k}$ is the $k$-th canonical correlation between activity coordinates and directional coordinates on subset $g$.

**Derivation.** Consider

$$
\max_{a,b}
\frac{a^\top G_{VD}b}
{\sqrt{a^\top G_{VV}a}\sqrt{b^\top G_{DD}b}}.
$$

Let

$$
\tilde a=G_{VV}^{1/2}a,
\qquad
\tilde b=G_{DD}^{1/2}b.
$$

The objective becomes

$$
\max_{\|\tilde a\|=\|\tilde b\|=1}
\tilde a^\top
R_{VD}
\tilde b,
$$

whose maximum is $\sigma_1(R_{VD})=\rho_1$; the remaining singular values follow in order under orthogonality constraints. $\blacksquare$

### 9.3 Canonical modes

The $k$-th activity mode's coefficient:

$$
a_{g,k}
=G_{VV}^{-1/2}u_{g,k}.
$$

The $k$-th directional mode:

$$
b_{g,k}
=G_{DD}^{-1/2}v_{g,k}.
$$

The corresponding observation-space signals:

$$
v^{*}_{g,k}
=\Phi_{\mathsf V}a_{g,k},
$$

$$
d^{*}_{g,k}
=\Phi_{\mathsf D}b_{g,k}.
$$

Under $W_g$:

$$
\|v^*_{g,k}\|_{W_g}
=\|d^*_{g,k}\|_{W_g}=1,
$$

and

$$
\boxed{
\langle v^*_{g,k},d^*_{g,k}\rangle_{W_g}
=\rho_{g,k}.
}
$$

So every $\rho_{g,k}$ corresponds to one concrete

$$
\mathsf V_k\leftrightarrow\mathsf D_k
$$

channel.

---

## 10. The directional explained-variance operator

### Definition

From Direction's point of view, collapse all canonical coupling into a single PSD operator:

$$
\boxed{
H_g^{\mathsf D\leftarrow\mathsf V} =
R_{g,\mathsf{VD}}^\top
R_{g,\mathsf{VD}}.
}
$$

Expanded:

$$
\boxed{
H_g^{\mathsf D\leftarrow\mathsf V} =
G_{g,\mathsf{DD}}^{-1/2}
G_{g,\mathsf{DV}}
G_{g,\mathsf{VV}}^{-1}
G_{g,\mathsf{VD}}
G_{g,\mathsf{DD}}^{-1/2}.
}
$$

### 10.1 Regression interpretation

Let the subset's activity coordinate matrix be $V_g$ and its direction coordinate matrix be $D_g$. Consider the population linear regression:

$$
D_g=V_gB+\varepsilon.
$$

The least-squares solution is

$$
B^*=G_{VV}^{-1}G_{VD}.
$$

The directional covariance explained linearly by activity space is

$$
\boxed{
\Sigma_{D,\rm expl}
=G_{DV}G_{VV}^{-1}G_{VD}.
}
$$

Whitening this under D-space's own covariance:

$$
G_{DD}^{-1/2}
\Sigma_{D,\rm expl}
G_{DD}^{-1/2}
=H^{D\leftarrow V}.
$$

So the eigenvalues of $H^{D\leftarrow V}$,

$$
\boxed{
\lambda_k(H^{D\leftarrow V})=\rho_k^2,
}
$$

can be read as the normalized fraction of the $k$-th canonical directional mode's covariance explained by activity space.

### 10.2 Boundary

"Explained" here means **linear geometric explanation within the same future-label system** — it is not causal explanation, and not point-in-time predictability.

Both $\mathsf V$ and $\mathsf D$ are constructed from future labels, so $H^{D\leftarrow V}$ can only describe label-side representation.

---

## 11. Coupling quality, rank, and concentration

### 11.1 Total normalized coupling

Define

$$
\boxed{
Q_g^{\mathsf D\leftarrow\mathsf V} =
\frac1{d_{D,g}^{\rm can}}
\operatorname{tr}
H_g^{\mathsf D\leftarrow\mathsf V} =
\frac1{d_{D,g}^{\rm can}}
\sum_k\rho_{g,k}^2,
}
$$

where $d_{D,g}^{\rm can}$ is the effective D dimension that survives the canonical-rank cutoff.

Range:

$$
0\le Q_g^{\mathsf D\leftarrow\mathsf V}\le1.
$$

A large $Q$ means a large fraction of normalized D-space energy has a systematic linear coupling with V-space.

### 11.2 Effective coupling rank

If $\sum_k\rho_k^2>0$, let

$$
p_{g,k}^{\rm can} =
\frac{\rho_{g,k}^2}
{\sum_j\rho_{g,j}^2}.
$$

Define

$$
\boxed{
r_{g,\rm eff}^{\mathsf D\leftarrow\mathsf V} =
\exp\left(
-\sum_kp_{g,k}^{\rm can}
\log p_{g,k}^{\rm can}
\right).
}
$$

### 11.3 Top-$m$ coupling concentration

$$
\boxed{
C_{g,m}^{\rm can} =
\frac{\sum_{k=1}^m\rho_{g,k}^2}
{\sum_j\rho_{g,j}^2}.
}
$$

$C_{g,1}^{\rm can}$ is one of the most important operational statistics here: if it's close to 1, most of the V-D relationship can be compressed into a single canonical channel.

### Four regimes

| $Q^{D\leftarrow V}$ | $r_{eff}^{D\leftarrow V}$ | Reading |
|---|---|---|
| high | low | strong, low-dimensional coupling — prioritize this universe |
| high | high | strong but multi-dimensional — needs a more complex representation |
| low | low | only one weak channel — probably not practically valuable |
| low | high | close to diffuse noise |

### No single composite score

This chapter deliberately does not define

$$
J=q_D+\alpha Q-\beta r_{eff}
$$

or any similar arbitrary single score. $q_D$, $Q$, rank, and OOS stability each represent a different property, and must be kept separate for eligibility and ranking diagnostics rather than collapsed into one number.

---

## 12. Ticker-level directional and coupling contribution

### Motivation

An arbitrary subset search has $2^N$ possibilities; directly optimizing $q_D$ or $Q^{D\leftarrow V}$ over that space is guaranteed to overfit. Candidate universes have to come from a pre-specified, low-complexity generator.

### 12.1 The directional contribution operator

For observation $(i,d)$, let the reference D row be

$$
\phi^{D}_{id}\in\mathbb R^{d_D}.
$$

Define the rank-one operator:

$$
\boxed{
K^{D}_{id}
=\phi^D_{id}\phi_{id}^{D\top}.
}
$$

Average over ticker $i$'s history in the discovery period:

$$
\boxed{
\bar K_i^D =
\frac1{T_i}
\sum_{d:i\in\mathcal U_d}
K_{id}^D.
}
$$

Its trace

$$
s_i^D
=\frac1{d_D}\operatorname{tr}(\bar K_i^D)
$$

is the ticker-level directional energy.

Define $s_i^V,s_i^L$ analogously, and directional purity:

$$
\boxed{
p_i^D =
\frac{s_i^D}{s_i^D+s_i^V+s_i^L+\epsilon}.
}
$$

### 12.2 Directional shape

To separate magnitude from shape:

$$
\boxed{
\widetilde K_i^D =
\frac{\bar K_i^D}
{\operatorname{tr}(\bar K_i^D)+\epsilon}.
}
$$

$\widetilde K_i^D$ describes *which* horizon combinations a ticker's directional energy lives in, independent of the total amount.

### 12.3 Cross contribution

An observation-level V-D outer product can also be defined:

$$
K_{id}^{VD} =
\phi^V_{id}\phi_{id}^{D\top}.
$$

with ticker-level average

$$
\bar K_i^{VD} =
\frac1{T_i}\sum_dK_{id}^{VD}.
$$

Its singular values cannot be treated directly as canonical coupling, since it hasn't yet been normalized by the subset-level $G_{VV}, G_{DD}$ — it only serves as a low-order descriptor for candidate generation.

---

## 13. Multi-universe candidate generation

### Principle

Candidate generation and candidate evaluation must stay separate. The generator only produces a finite, pre-specified family; statistical and OOS validation happen afterward, in Section 15-84.

### 13.1 D-purity / D-energy ladder

Rank discovery-period ticker statistics by

$$
\text{rank by }p_i^D
$$

or by

$$
\text{rank by }u_i^D =
\log(s_i^D+\epsilon)
-\frac12\log(s_i^V+\epsilon)
-\frac12\log(s_i^L+\epsilon).
$$

Pre-specify a size grid,

$$
\mathcal M=\{M_1,\dots,M_J\},
$$

and generate nested universes:

$$
\mathcal S_j
=\operatorname{TopM_j}(u_i^D).
$$

This line finds D-enriched candidate families.

### 13.2 Directional-shape clustering

Using the Frobenius distance between $\widetilde K_i^D$ matrices,

$$
d_D(i,j)
=\|\widetilde K_i^D-\widetilde K_j^D\|_F,
$$

cluster with a fixed, pre-specified method (e.g. average linkage or spectral clustering). Each cluster becomes a candidate universe.

This line finds non-nested universes that differ in their D *shape*.

### 13.3 VD-shape clustering

Concatenate

$$
\operatorname{vec}(\widetilde K_i^D),
\quad
\operatorname{vec}(\bar K_i^{VD})
$$

run discovery-only PCA on the concatenation, then cluster. This produces candidate groups that share similar raw V-D cross-structure.

### 13.4 Hybrid family

Within a directional-shape cluster $c$, further split by top fraction of $u_i^D$:

$$
\mathcal S_{c,q} =
\{i\in c:u_i^D\ge Q_q(u^D\mid c)\}.
$$

This lets you retain, simultaneously:

- different D directions;
- different enrichment levels within the same direction.

### Rule 80.1 (No arbitrary subset search)

Genetic search, beam search, or greedily maximizing evaluation-period $Q$ ticker by ticker are not permitted, unless the entire optimizer is fully embedded inside a permutation null with an independent OOS test. The baseline specification does not allow it.

---

## 14. Universe eligibility: the geometric layer

Geometric eligibility only defines "worth submitting to statistical validation" — it is not a significance criterion.

**U1 — size**

$$
\operatorname{median}_dN_{g,d}\ge N_{\min}.
$$

**U2 — directional enrichment**

$$
q_g^D>1.
$$

**U3 — directional purity**

$$
f_g^D>1/3.
$$

**U4 — coupling presence**

$$
\rho_{g,1}>\rho_{\min}^{can}.
$$

**U5 — coupling concentration**

$$
C_{g,1}^{can}>C_{\min}^{can}
$$

or

$$
r_{g,eff}^{D\leftarrow V}\le r_{\max}^{can}.
$$

All thresholds are locked before research begins. Clearing U1-U5 still does not mean statistically significant.

---

## 15. The universe null: is membership actually special?

### 15.1 Null hypothesis

For a fixed candidate $\mathcal S_g$, the null is:

> Under the same dates, the same subset size, and the same pre-locked liquidity/availability strata, random membership can produce $q_D$, $Q^{D\leftarrow V}$, or concentration that is equally or more extreme.

### 15.2 Size-matched stratified membership permutation

For every date $d$:

1. fix $N_{g,d}$;
2. stratify the parent universe by **pre-locked** strata such as liquidity, market cap, or listing age;
3. within each stratum, draw the same count as the real candidate;
4. obtain $\mathcal S_g^{(b)}$;
5. recompute the full chain:
   $$
   G_g^{(b)},
   R_{VD}^{(b)},
   \rho_k^{(b)},
   Q^{(b)},
   r_{eff}^{(b)}.
   $$

You cannot just permute $G_{VD}$ alone, because the canonical normalizers $G_{VV}, G_{DD}$ also change with membership.

### 15.3 Primary test statistics

Report a p-value for at least these three:

$$
T_1=q_g^D,
$$

$$
T_2=Q_g^{D\leftarrow V},
$$

$$
T_3=\rho_{g,1}^2.
$$

If the research question is specifically about low-rank concentration, add

$$
T_4=C_{g,1}^{can}.
$$

For a "larger is more special" statistic:

$$
\boxed{
p_g
=\frac{1+\sum_b\mathbf1\{T_g^{(b)}\ge T_g^{obs}\}}
{B_U+1}.
}
$$

### 15.4 Multiple testing

Because multiple universes are evaluated at once, the primary approach uses Westfall-Young max-T: each membership permutation produces the statistic for every candidate at once, then take

$$
M^{(b)}=\max_g T_g^{(b)}.
$$

Family-wise error is then calibrated from the distribution of $M^{(b)}$.

This is especially necessary for a nested ladder, since adjacent universes are highly correlated.

### 15.5 When the candidate generator is itself derived from labels

If membership itself was generated from label-derived statistics such as $u_i^D$ or $\widetilde K_i^D$, a fixed-membership null only answers "is this already-chosen set special" — it does not account for search cost.

There are two legitimate ways to do formal inference here:

1. **Discovery/evaluation split (primary approach)**: the generator only runs on $\mathcal D_U$; membership is fixed and evaluated on $\mathcal E_U$.
2. **Full rerun null (expensive alternative)**: every permutation reruns the entire chain — ticker scores, clustering, candidate generation, and pruning — from scratch.

The baseline specification uses approach 1.

---

## 16. Temporal stability

### 16.1 Moving-block bootstrap

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

### 16.2 Enrichment persistence

Define:

$$
\boxed{
\pi_{q_D}
=\Pr_b(q_D^{(b)}>1).
}
$$

### 16.3 Coupling persistence

For the top canonical correlation:

$$
\boxed{
\pi_{\rho_1}
=\Pr_b(\rho_1^{(b)}>\rho_{\min}^{can}).
}
$$

Also report the median and the 10th/90th percentiles.

### 16.4 Directional-mode stability

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

## 17. The train / freeze / OOS canonical test

### Motivation

A bootstrap still resamples within the same evaluation-period distribution. The test that actually matters is whether the canonical channel survives on a future period that genuinely did not participate in discovery.

### 17.1 Freeze

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

### 17.2 OOS evaluation

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

### 17.3 Decay

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

## 18. Redundancy and retention across multiple universes

### Motivation

Several candidate universes can all pass eligibility simultaneously while really just being different cuts of the same tickers, or the same D mode. Membership redundancy and mechanism redundancy must be kept separate.

### 18.1 Membership overlap

Jaccard index:

$$
\boxed{
J_{gh} =
\frac{|\mathcal S_g\cap\mathcal S_h|}
{|\mathcal S_g\cup\mathcal S_h|}.
}
$$

For dynamic universes, report the median-date Jaccard along with p10/p90.

### 18.2 Canonical directional similarity

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

### 18.3 Retention rule

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

## 19. Local refit of an accepted universe

### Motivation

The shared $\Phi_{\mathcal U}$'s role is discovery and cross-universe comparison. Once $\mathcal S_g$ is formally accepted, the downstream feature-label-alpha pipeline should treat it as its own market domain, and rebuild a local basis for it.

### 19.1 The local label basis

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

### 19.2 Reference-to-local diagnostic

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

### 19.3 Interpretation

- small angle: the subset's local mechanism space is roughly a reweighting of the full-universe mechanism;
- large angle: the subset has meaningfully different local label geometry;
- either outcome is acceptable to retain — it just means downstream feature interpretation must use the local $\Phi_g$, not the reference basis.

---

## 20. From canonical D modes to feature research

### Motivation

Up to this point, universe geometry has not used any engineered feature at all. An accepted universe's canonical D modes can serve as an explicit target for subsequent feature research, but they cannot replace full mechanism analysis.

### 20.1 The reference canonical target

For universe $g$'s stable top D canonical coefficient $b_{g,1}$:

$$
\boxed{
d_g^*
=\Phi_{\mathsf D}b_{g,1}.
}
$$

This represents the strongest V-coupled directional combination in this universe.

### 20.2 Feature alignment

For a point-in-time feature $x_p$:

$$
\boxed{
A_{p,g}^{can}
=x_p^\top W_gd_g^*.
}
$$

This must go through the same permutation significance, temporal stability, and OOS decay checks as any other target alignment elsewhere on this site.

### 20.3 The complete downstream pipeline

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

---

## 21. Beta, demeaning, and the factor-neutral variant

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

## 22. Cross-sectional breadth and single-ticker dominance

### Motivation

A universe's high $q_D$ or high $\rho_1$ might be determined by a tiny handful of names. Cross-sectional breadth has to be checked explicitly.

### 22.1 D-energy contribution

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

### 22.2 Canonical-channel breadth

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

## 23. The dynamic universe

### Motivation

If universe membership updates over time, the selection rule must depend only on label history already realized strictly before that point in time.

### 23.1 Rolling ticker operators

For ticker $i$:

$$
\bar K_{i,t}^{D} =
\sum_{\tau\lt t}\omega_{t-\tau}
K_{i\tau}^{D},
$$

$$
\bar K_{i,t}^{VD} =
\sum_{\tau\lt t}\omega_{t-\tau}
K_{i\tau}^{VD}.
$$

$\omega$ can be a fixed rolling window or exponential decay, and must be locked in advance.

### 23.2 Membership update

$$
\mathcal S_{g,t}
=\mathcal A_g
(\{\bar K_{i,t}^{D},\bar K_{i,t}^{VD}\}_i).
$$

Update frequency should default to weekly or monthly — not updating membership on every 5-minute bar.

### 23.3 Turnover

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

### 23.4 Freeze boundary

The dynamic rule's functional form, window, cluster count, size grid, and thresholds are all frozen once discovery is complete. The evaluation period can update the *statistics*, but it cannot modify the *rule*.

---

## 24. Universe report format

Every candidate — accepted or rejected — gets a row. Only saving the winners is not permitted.

### A. Identity

| Field | Content |
|---|---|
| `universe_id` | unique ID |
| `generator` | ladder / D-shape / VD-shape / hybrid |
| `parent_cluster` | if applicable |
| `size_median` | median $N_{g,d}$ |
| `size_p10`, `size_p90` | dynamic universe size range |
| `turnover_median` | dynamic membership turnover |

### B. Mechanism energy

| Field | Definition |
|---|---|
| `q_D` | $\operatorname{tr}(G_{DD})/d_D$ |
| `q_V` | $\operatorname{tr}(G_{VV})/d_V$ |
| `q_L` | $\operatorname{tr}(G_{LL})/d_L$ |
| `f_D,f_V,f_L` | dimension-adjusted fingerprint |
| `r_eff_D` | D spectrum effective rank |
| `C_D_1`, `C_D_2` | top-1 / top-2 D concentration |

### C. Canonical coupling

| Field | Definition |
|---|---|
| `rho_1...rho_r` | canonical correlations |
| `Q_D_from_V` | mean $\rho_k^2$ |
| `r_eff_can` | coupling effective rank |
| `C_can_1`, `C_can_2` | top coupling concentration |
| `can_D_loadings_1` | top D canonical coefficients |
| `can_V_loadings_1` | top V canonical coefficients |

### D. Inference

| Field | Content |
|---|---|
| `p_qD_raw`, `p_qD_adj` | D-enrichment null |
| `p_Q_raw`, `p_Q_adj` | coupling null |
| `p_rho1_raw`, `p_rho1_adj` | top canonical correlation null |
| `pi_qD` | temporal persistence |
| `pi_rho1` | coupling persistence |
| `subspace_stability` | principal-angle bootstrap stability |

### E. OOS

| Field | Content |
|---|---|
| `rho1_train` | frozen train canonical correlation |
| `rho1_oos` | frozen OOS canonical covariance/correlation |
| `rho_decay` | $|\rho_{OOS}|/(|\rho_{train}|+\epsilon)$ |
| `rho_sign_consistency` | orientation-consistent sign |
| `qD_oos` | OOS D enrichment |

### F. Breadth / factor diagnostics

| Field | Content |
|---|---|
| `breadth_D_median` | D contributor breadth |
| `breadth_can_median` | canonical contribution breadth |
| `qD_idio` | if a factor-neutral variant exists |
| `rho1_idio` | if a factor-neutral variant exists |

### G. Decision

| Field | Content |
|---|---|
| `eligible_geometry` | U1-U5 |
| `significant` | permutation / max-T |
| `stable` | bootstrap + OOS |
| `family_id` | universe redundancy family |
| `accepted` | whether it ultimately enters the local feature pipeline |
| `reject_reason` | size / null / instability / redundancy / breadth, etc. |

---

## 25. Complete implementation order

The following can be used directly as a pipeline checklist.

### Stage 0: Setup

```text
 0. Lock the horizon band, parent universe, and discovery/evaluation/OOS windows.
 1. Lock candidate generators, size grid, clustering method, and cluster-count range.
 2. Lock the canonical rank cutoff, B_U, B_time, and all eligibility thresholds.
 3. Lock the dynamic universe's rolling window / decay and update frequency (if used).
```

### Stage 1: Reference basis

```text
 4. On the parent universe's discovery period, within-date demean the labels.
 5. Fit Sigma_Y, K_eff, B_Y; transform -> Z.
 6. Phi_D := Z; build Phi_V, Phi_T; verify Phi^T W_U Phi = I.
 7. Freeze the reference fit objects.
```

### Stage 2: Candidate generation

```text
 8. Compute ticker-level K_i^D, s_i^D/V/L, p_i^D, u_i^D.
 9. Compute the normalized D shape Ktilde_i^D; the K_i^VD descriptor if needed.
10. Generate ladder / D-shape / VD-shape / hybrid candidate families per the pre-specification.
11. Drop candidates below N_min; do not adjust the family based on evaluation performance.
```

### Stage 3: Discovery geometry

```text
12. Build W_g for each candidate.
13. G_g = Phi^T W_g Phi; split into DD/VV/VD etc. blocks.
14. Compute q_D/q_V/q_L, fingerprint, D eigenspectrum, r_eff_D.
15. Build G_VV^{-1/2}, G_DD^{-1/2} (rank-controlled).
16. R_VD; SVD -> rho_k, U_k, V_k.
17. H_D<-V = R^T R; Q_D<-V, r_eff_can, C_can_m.
18. Apply geometric eligibility U1-U5.
```

### Stage 4: Evaluation inference

```text
19. Freeze membership / the membership rule and canonical rank choices.
20. Transform evaluation labels into the frozen reference Phi.
21. Recompute G_g, q_D, R_VD, rho, Q, rank.
22. Size-matched stratified membership null; recompute every canonical object on each draw.
23. Westfall-Young max-T across the candidate family.
24. Moving-block bootstrap -> pi_qD, pi_rho1, canonical subspace stability.
25. D / canonical breadth diagnostics.
```

### Stage 5: The real OOS freeze test

```text
26. At the end of train/evaluation, freeze canonical coefficients a_gk, b_gk.
27. On OOS, only transform Phi; CCA may not be rerun.
28. Compute frozen rho_OOS, rho_decay, sign consistency, qD_OOS.
29. Decide stable / unstable per the pre-registered stability rule.
```

### Stage 6: Multi-universe pruning

```text
30. Compute pairwise Jaccard membership overlap.
31. Compute canonical D-subspace principal-angle similarity.
32. Merge into the same family only when both membership and D mechanism are highly duplicated.
33. Use the lexicographic rule to keep a representative within each family; multiple families may coexist.
```

### Stage 7: Local refit and the feature pipeline

```text
34. For each accepted universe, rebuild a local Phi_g on its own model-training window.
35. Report reference-to-local principal angles.
36. Start the feature library, C_g, Q, fingerprint, stability, and selection on the local Phi_g.
37. Additionally check feature alignment against the frozen canonical D target.
38. Model and alpha work from here on follows the existing Feature-Label-Alpha framework unchanged.
```

---

## 26. The fit/transform boundary table

| Object | Fit on | Applied to | Refittable OOS? |
|---|---|---|---|
| $\Sigma_Y,B_Y,\Phi_{\mathcal U}$ | universe-discovery train | evaluation / OOS transform | no |
| ticker $\bar K_i^D,\bar K_i^{VD}$ | discovery history | candidate generation | no; dynamic version only uses past rolling updates |
| candidate generator | discovery | evaluation / OOS | no |
| static membership $\mathcal S_g$ | discovery | evaluation / OOS | no |
| dynamic membership rule | discovery | future times | rule is immutable, state may update |
| canonical rank cutoff | discovery | evaluation / OOS | no |
| $U_g,V_g,\rho_g$ for selection | train/evaluation | diagnosis | may not be redone on OOS to re-pick the top |
| frozen canonical coefficients $a_g,b_g$ | pre-OOS | OOS | no |
| local $\Phi_g$ | accepted-universe model train | downstream model OOS | per each model's own fold specification |
| factor-neutral basis | its corresponding train window | its corresponding evaluation | kept separate from the raw basis |

---

## 27. Numerical identity tests

All of the following must be automated in the production implementation.

| Test | Criterion |
|---|---|
| Reference orthogonality | $\|\Phi_U^TW_U\Phi_U-I\|_F<10^{-6}$ |
| Mother-universe energy | $|q_U^D-1|,|q_U^V-1|,|q_U^L-1|<10^{-8}$ |
| Mother cross-block | $\|G_{U,VD}\|_F<10^{-8}$ |
| PSD of $G_g$ | $\lambda_{min}(G_g)>-10^{-10}$ |
| Canonical bound | $\rho_k\le1+10^{-8}$ |
| $H$ identity | $\|H-R^TR\|_F<10^{-10}$ |
| Eigen/SVD match | eigenvalues$(H)$ = $\rho_k^2$ within $10^{-10}$ |
| CCA normalization | canonical V/D signals have unit $W_g$ norm |
| CCA covariance | $\langle v_k^*,d_k^*\rangle_{W_g}=\rho_k$ |
| Effective rank bound | $1\le r_{eff}\le$ retained rank |
| Top concentration monotonic | $C_1\le C_2\le\cdots\le1$ |
| Sign invariance | flipping singular-vector signs leaves every scalar metric unchanged |
| Membership identity | $\mathcal S_g=\mathcal U\Rightarrow G_g=I$ |

---

## 28. Synthetic calibration tests

### T1: Random membership null

Generate membership of the same size as a candidate, randomly, on real label geometry.

Criterion: with no planted special subset, $p_{qD},p_Q,p_{\rho_1}$ are roughly uniform; at FWER 0.05, the family false-positive rate is $\le5\%$.

### T2: Planted D-enrichment

Generate a subset whose D-block covariance is multiplied by $c>1$, with V-D coupling left unchanged.

Expected:

$$
q_D\uparrow,
\quad
Q^{D\leftarrow V}\approx\text{baseline}.
$$

This verifies that D enrichment and coupling are correctly separated from each other.

### T3: Planted rank-one V-D channel

Generate:

$$
d_1=\rho v_1+\sqrt{1-\rho^2}\epsilon,
$$

with the remaining D directions independent of V.

Expected:

$$
\hat\rho_1\approx\rho,
$$

$$
\rho_{k>1}\approx0,
$$

$$
r_{eff}^{can}\approx1.
$$

### T4: Diffuse coupling

Plant $r$ equal canonical correlations.

Expected:

$$
r_{eff}^{can}\approx r,
$$

and $C_1^{can}\approx1/r$.

### T5: Subset-scale only

Raise only the subset's D and V variance, without changing the normalized correlation.

Expected: raw $\|G_{VD}\|$ can rise, but the normalized $\rho_k$ should not rise systematically. This verifies that canonical whitening is actually necessary.

### T6: Search-cost calibration

Use discovery labels to generate multiple ladders/clusters, with evaluation labels independent of discovery.

Expected: using the discovery/evaluation split, the adjusted p-value stays calibrated; if selection and evaluation happen on the same period, the false-positive rate inflates noticeably.

### T7: Frozen OOS channel

Plant the same canonical channel in both train and test.

Criterion: the frozen $\rho_{OOS}$ recovers the true value; if the test-period channel is rotated 90°, the frozen $\rho_{OOS}$ should be close to 0 — even though rerunning CCA on the test period would still find a high top singular value. This test specifically verifies "OOS may not re-pick a canonical mode."

### T8: Beta-driven subset

Generate

$$
r_i=\beta_ir_m+\epsilon_i
$$

and give the candidate a larger $\operatorname{Var}(\beta)$.

Expected: raw $q_D$ rises; the factor-neutral variant's enrichment disappears or drops sharply.

### T9: Dynamic universe

Let $N_d$ and membership change over time.

Criterion: $W_g$'s total weight per date is exactly $1/T_g$; null calibration is not broken by changes in $N_d$.

### T10: Single-name domination

Plant one extreme ticker that produces most of the coupling.

Expected: $\rho_1$ can be very high, but canonical breadth drops sharply. This verifies the breadth diagnostic can actually detect concentration.

---

## 29. Default operating parameters

| Parameter | Default | Note |
|---|---:|---|
| discovery / evaluation | 60 / 40 | universe-discovery inference |
| purge gap | $\ge H_{max}$ | label overlap |
| $N_{min}$ | $\max(100,10d_\Phi)$ | minimum candidate median size |
| ladder size grid | 100, 200, 300, 500, 800, 1200 | adjust to parent universe, lock in advance |
| canonical eigen cutoff $\epsilon_{can}$ | $10^{-6}$ relative | more conservative than a pure numerical orth cutoff |
| $B_U$ membership null | 2000 | at least 1000 if many candidates |
| multiple testing | Westfall-Young max-T | family-wise control |
| $B_{time}$ | 1000 | temporal bootstrap |
| block length | $\ge H_{max}$ | overlapping label dependence |
| $q_D$ geometry eligibility | $>1$ | not a significance criterion |
| $f_D$ geometry eligibility | $>1/3$ | not a significance criterion |
| $\rho_{min}^{can}$ | 0.20 | initial operating threshold, needs calibration |
| $C_{min}^{can}$ | 0.50 | top canonical pair explains at least half the coupling energy |
| $\pi_{qD}$ | $\ge0.80$ | bootstrap persistence |
| $\pi_{\rho_1}$ | $\ge0.80$ | bootstrap persistence |
| D-subspace stability | $\ge0.70$ | bootstrap principal-angle similarity |
| Jaccard redundancy | $>0.80$ | used jointly with D similarity |
| D-subspace redundancy | $>0.90$ | used jointly with Jaccard |
| dynamic update | weekly / monthly | reduces membership churn |
| numerical $\epsilon$ | $10^{-12}$ | denominator protection |

All of these thresholds are operating/screening parameters, not statistical criteria. Formal significance comes from the Section 15 null; stability comes from Section 16-84; actual feature learnability comes from the downstream pipeline starting at Section 20.

---

## 30. The four core objects of a universe

For every candidate universe $g$, ultimately only four objects need to be remembered.

### 1. Restriction geometry

$$
\boxed{
G_g=\Phi_U^TW_g\Phi_U.
}
$$

Answers: how are the reference label mechanisms reweighted on this subset?

### 2. The directional block

$$
\boxed{
G_{g,DD}.
}
$$

Answers: how strong is D energy, and which horizon combinations does it concentrate in?

### 3. Canonical coupling

$$
\boxed{
R_{g,VD}
=G_{VV}^{-1/2}G_{VD}G_{DD}^{-1/2}.
}
$$

Answers: what is the normalized cross-geometry between V and D?

### 4. The directional explained-variance operator

$$
\boxed{
H_g^{D\leftarrow V}=R_{g,VD}^TR_{g,VD}.
}
$$

Answers: which D directions have a low-dimensional activity-side representation in this subset?

These four objects are not interchangeable — each answers a question the others cannot.

---

## 31. The parallel structure across universe, feature, and alpha

The existing framework, for a signal $s$:

$$
c_s=s^TW\Phi.
$$

A universe changes the *measure* rather than providing a one-dimensional signal, so its natural object is:

$$
G_g=\Phi^TW_g\Phi.
$$

The two form the following correspondence:

| Feature / Alpha Geometry | Universe Geometry |
|---|---|
| $c_s=s^TW\Phi$ | $G_g=\Phi^TW_g\Phi$ |
| $Q_s=\|c_s\|^2$ | $q_g^B=\operatorname{tr}(G_{BB})/d_B$ |
| mechanism fingerprint | universe fingerprint |
| $c_s/\|c_s\|$ | enriched D eigenvectors |
| feature mechanism rank / coverage | D effective rank |
| feature-label alignment | V-D canonical correlation |
| feature stability | universe canonical-subspace stability |
| feature breadth | D / canonical contributor breadth |

So the complete research order is:

$$
\boxed{
\text{Universe Geometry}
\longrightarrow
\text{Local Feature Geometry}
\longrightarrow
\text{Alpha Geometry.}
}
$$

---

## 32. What this chapter can and cannot answer

### Can answer

1. Does a given subset have higher directional label energy?
2. Is directional energy concentrated in a small number of modes?
3. Do activity and direction exhibit normalized canonical coupling within the subset?
4. Is that coupling low-rank?
5. Is the canonical D subspace stable across time?
6. Do different universes concentrate different D mechanisms?
7. Is the structure driven mainly by beta / common-factor dispersion?

### Cannot answer

1. Whether any point-in-time feature exists that can predict these D modes;
2. causal mechanism;
3. whether there is positive net Sharpe after trading costs;
4. whether a given universe is "optimal";
5. the optimal architecture for the local model;
6. whether the canonical V-D relation represents V preceding D in economic time ordering.

That last point matters especially: $\Phi_V$ and $\Phi_D$ are both constructed from the same future-label system. $V \leftrightarrow D$ coupling is a representation relation, not a temporal causal arrow.

---

## 33. Minimum viable version (MVP)

A first implementation does only the following ten steps:

1. Build a fixed $\Phi_U$ on the parent universe's discovery period.
2. Compute $s_i^D,s_i^V,s_i^L,p_i^D$ per ticker.
3. Build 5-8 ladder universes using a pre-specified top-$M$ grid.
4. For each universe, compute
   $$G_g,\quad q_D,q_V,q_L,\quad f_D,\quad r_{eff}^D.$$
5. Compute
   $$R_{VD}=G_{VV}^{-1/2}G_{VD}G_{DD}^{-1/2}.$$
6. Report
   $$\rho_1,\rho_2,\quad Q^{D\leftarrow V},\quad r_{eff}^{can},\quad C_1^{can}.$$
7. Run the size-matched membership null + Westfall-Young on the evaluation period.
8. Use a moving-block bootstrap to check $q_D$, $\rho_1$, and top-D-canonical-mode stability.
9. Freeze the top canonical coefficients, and test $\rho_1^{OOS}$ on genuine OOS data.
10. For each universe that passes, rebuild its own local $\Phi_g$, and only then begin feature engineering.

The MVP's accept condition is deliberately not a single score — it's four gates:

$$
\boxed{
\text{D enrichment}
\land
\text{canonical coupling}
\land
\text{low-rank / stable representation}
\land
\text{OOS persistence}.
}
$$

Only a universe that clears all four gates proceeds to downstream feature research — and "actually learnable" still has to be confirmed by the feature-label-alpha OOS results that follow.

---

## Appendix U-A: Core formula reference

**Reference**

$$
\Phi_U^TW_U\Phi_U=I.
$$

**Subset measure**

$$
w_{id}^{(g)}
=\frac{\mathbf1\{i\in\mathcal S_{g,d}\}}{T_gN_{g,d}}.
$$

**Restriction Gram**

$$
G_g=\Phi_U^TW_g\Phi_U.
$$

**D enrichment**

$$
q_g^D=\frac1{d_D}\operatorname{tr}(G_{g,DD}).
$$

**D effective rank**

$$
r_{eff}^D
=\exp\left(-\sum_kp_k^D\log p_k^D\right),
\qquad
p_k^D=\frac{\lambda_k^D}{\sum_j\lambda_j^D}.
$$

**Canonical V-D operator**

$$
R_{g,VD}
=G_{g,VV}^{-1/2}
G_{g,VD}
G_{g,DD}^{-1/2}.
$$

**Canonical correlations**

$$
R_{g,VD}=U\operatorname{diag}(\rho_k)V^T.
$$

**D explained-variance operator**

$$
H_g^{D\leftarrow V}
=R_{g,VD}^TR_{g,VD}.
$$

**Coupling quality**

$$
Q_g^{D\leftarrow V}
=\frac1{d_{D,g}^{can}}\sum_k\rho_k^2.
$$

**Coupling effective rank**

$$
r_{eff}^{D\leftarrow V}
=\exp\left(-\sum_kp_k^{can}\log p_k^{can}\right),
$$

$$
p_k^{can}=\frac{\rho_k^2}{\sum_j\rho_j^2}.
$$

**Top-$m$ concentration**

$$
C_{g,m}^{can}
=\frac{\sum_{k=1}^m\rho_k^2}{\sum_j\rho_j^2}.
$$

**Universe mechanism similarity**

$$
S_{gh}^{D}
=\frac1{r_*}\|V_g^TV_h\|_F^2.
$$

**OOS frozen channel**

$$
\rho_{g,k}^{OOS}
=\operatorname{Corr}_{W_g^{OOS}}
(\Phi_Va_{g,k}^{train},\Phi_Db_{g,k}^{train}).
$$

---

## Appendix U-B: Interpretation quick-reference

| Pattern | Geometric meaning | Next step |
|---|---|---|
| $q_D\gg1$, $Q_{D\leftarrow V}\approx0$ | D energy is large but has no low-dimensional coupling with activity | can retain as a D-enriched domain; don't call it structural coupling |
| $q_D>1$, $Q$ high, $r_{eff}^{can}\approx1$ | D is strong and almost entirely represented by a single V-D channel | prioritize the frozen OOS canonical test |
| $q_D\approx1$, $Q$ high | D isn't especially strong, but the V-D relation is especially coherent | can be kept as a specialized universe rather than a D-enrichment universe |
| high $q_D$, high $r_{eff}^D$ | D is strong but diffuse | needs multi-dimensional features/model; not a low-dimensional representation |
| high $\rho_1^{train}$, $\rho_1^{OOS}\approx0$ | canonical-mode selection overfit / regime-specific | reject on stability |
| strong raw, weak idio | driven by common-factor dispersion | flag as factor-driven |
| high membership overlap, low D similarity | same names, different mechanism | can keep both |
| low membership overlap, high D similarity | different names concentrating the same D mode | depending on research goal, keep both domains or merge the family |
| high $\rho_1$ but very low canonical breadth | a handful of tickers dominate | concentration warning |
| large local-reference angle | the subset has its own label basis | downstream must use the local $\Phi_g$ |

---

## Appendix U-C: One-line summary

The entire logic of this chapter compresses to:

$$
\boxed{
\mathcal U
\xrightarrow{\text{reference }\Phi}
\mathcal S_g
\xrightarrow{G_g}
G_{DD},G_{VV},G_{VD}
\xrightarrow{\text{CCA}}
H_g^{D\leftarrow V}
\xrightarrow{\text{null + stability + OOS}}
\mathscr S_{accepted}
\xrightarrow{\text{local refit}}
\Phi_g
\xrightarrow{X,\alpha}
\text{true predictive learnability.}
}
$$

The central question of universe discovery is not "which group of stocks will have the biggest future return" — it's:

> **Which market sub-domains have directional label structure that concentrates into stable, low-dimensional canonical modes representable by activity geometry — and does that structure still hold on genuine out-of-sample data?**
