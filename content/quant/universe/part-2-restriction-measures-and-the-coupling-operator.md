---
title: "Multi-Sub-Universe Directional Coupling Geometry, Part 2: Restriction Measures and the Coupling Operator"
description: "How restricting to a sub-universe reshapes the reference label geometry: the restriction Gram operator, directional enrichment and compressibility, and the canonical Activity-Direction coupling operator that is the core object of this chapter."
date: 2026-09-14
series: ["Multi-Sub-Universe Directional Coupling Geometry"]
part: 2
tags: ["universe", "cross-sectional", "methodology", "quant-research"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Multi-Sub-Universe Directional Coupling Geometry" >}}

## §71. Sub-universe weight and the restriction measure

### §71.1 Definition

For a candidate universe $\mathcal S_g$, let

$$
N_{g,d}=|\mathcal S_{g,d}|.
$$

Define the renormalized row weight:

$$
\boxed{
w^{(g)}_{id}
=
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

If $N_{g,d}<d_\Phi+5$, the Gram matrices may be numerically computable but statistically very unstable. This chapter requires at minimum:

$$
\operatorname{median}_dN_{g,d}\ge N_{\min},
$$

where $N_{\min}$ should be much larger than $d_\Phi$; see §96 for the default.

---

## §72. The restriction Gram operator

### Definition 72.1

For any candidate universe $g$:

$$
\boxed{
G_g
=
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
v^\top G_gv
=
(\Phi_{\mathcal U}v)^\top
W_g
(\Phi_{\mathcal U}v)
\ge0.
$$

$\blacksquare$

### §72.1 Physical meaning of the Rayleigh quotient

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

### §72.2 Not a learnability operator

$G_g$ is built entirely from future labels. A large eigenvalue only means label variation / mechanism energy is concentrated — it says nothing about whether any point-in-time $X$ exists that could actually predict that direction.

So $G_g$ is a **label-geometry operator**, not a predictive operator.

---

## §73. Block decomposition and directional enrichment

By the block structure of the reference basis:

$$
G_g=
\begin{pmatrix}
G_{g,\mathsf{DD}} & G_{g,\mathsf{DV}} & G_{g,\mathsf{DL}}\\
G_{g,\mathsf{VD}} & G_{g,\mathsf{VV}} & G_{g,\mathsf{VL}}\\
G_{g,\mathsf{LD}} & G_{g,\mathsf{LV}} & G_{g,\mathsf{LL}}
\end{pmatrix}.
$$

### §73.1 Average mechanism energy per dimension

Because block dimensions differ, traces can't be compared directly. Define

$$
\boxed{
q_g^{\mathsf B}
=
\frac1{d_{\mathsf B}}
\operatorname{tr}(G_{g,\mathsf{BB}}),
\qquad
\mathsf B\in\{\mathsf D,\mathsf V,\mathsf L\}.
}
$$

The parent universe satisfies, exactly:

$$
q_{\mathcal U}^{\mathsf D}
=q_{\mathcal U}^{\mathsf V}
=q_{\mathcal U}^{\mathsf L}
=1.
$$

So $q_g^{\mathsf D}=1.7$ means:

> this sub-universe's average squared amplitude per reference directional dimension is $1.7\times$ the parent universe's.

### §73.2 Universe mechanism fingerprint

Define

$$
\boxed{
f_g^{\mathsf B}
=
\frac{q_g^{\mathsf B}}
{q_g^{\mathsf D}+q_g^{\mathsf V}+q_g^{\mathsf L}}.
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

## §74. The directional spectrum and compressibility

### Motivation

A large $q_g^{\mathsf D}$ doesn't by itself mean the directional representation is simple. If all $d_{\mathsf D}$ directions are scaled up together proportionally, there is no low-dimensional structure.

### §74.1 The directional eigenspectrum

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

### §74.2 Effective directional rank

Let

$$
p_{g,k}^{\mathsf D}
=
\frac{\lambda_{g,k}^{\mathsf D}}
{\sum_j\lambda_{g,j}^{\mathsf D}}.
$$

Define the entropy effective rank:

$$
\boxed{
r_{g,\rm eff}^{\mathsf D}
=
\exp\left(
-\sum_kp_{g,k}^{\mathsf D}
\log p_{g,k}^{\mathsf D}
\right).
}
$$

$$
1\le r_{g,\rm eff}^{\mathsf D}\le d_{\mathsf D}.
$$

### §74.3 Top-$m$ directional concentration

$$
\boxed{
C_{g,m}^{\mathsf D}
=
\frac{\sum_{k=1}^m\lambda_{g,k}^{\mathsf D}}
{\sum_j\lambda_{g,j}^{\mathsf D}}.
}
$$

and

$$
\boxed{
r_{g,\eta}^{\mathsf D}
=
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

## §75. From $|Z|^\top WZ$ to $\mathsf V \leftrightarrow \mathsf D$ coupling

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

### §75.1 The cross-block Gram

Define

$$
\boxed{
G_{g,\mathsf{VD}}
=
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

## §76. The canonical $\mathsf V \leftrightarrow \mathsf D$ coupling operator

### Motivation

Looking at $G_{g,\mathsf{VD}}$ directly still conflates the coupling with $\mathsf V$'s and $\mathsf D$'s own scale within the subset. Both sides need to be whitened first.

### §76.1 Regularized inverse square root

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

### §76.2 Definition

$$
\boxed{
R_{g,\mathsf{VD}}
=
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

### §76.3 Canonical modes

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

## §77. The directional explained-variance operator

### Definition

From Direction's point of view, collapse all canonical coupling into a single PSD operator:

$$
\boxed{
H_g^{\mathsf D\leftarrow\mathsf V}
=
R_{g,\mathsf{VD}}^\top
R_{g,\mathsf{VD}}.
}
$$

Expanded:

$$
\boxed{
H_g^{\mathsf D\leftarrow\mathsf V}
=
G_{g,\mathsf{DD}}^{-1/2}
G_{g,\mathsf{DV}}
G_{g,\mathsf{VV}}^{-1}
G_{g,\mathsf{VD}}
G_{g,\mathsf{DD}}^{-1/2}.
}
$$

### §77.1 Regression interpretation

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

### §77.2 Boundary

"Explained" here means **linear geometric explanation within the same future-label system** — it is not causal explanation, and not point-in-time predictability.

Both $\mathsf V$ and $\mathsf D$ are constructed from future labels, so $H^{D\leftarrow V}$ can only describe label-side representation.

---

## §78. Coupling quality, rank, and concentration

### §78.1 Total normalized coupling

Define

$$
\boxed{
Q_g^{\mathsf D\leftarrow\mathsf V}
=
\frac1{d_{D,g}^{\rm can}}
\operatorname{tr}
H_g^{\mathsf D\leftarrow\mathsf V}
=
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

### §78.2 Effective coupling rank

If $\sum_k\rho_k^2>0$, let

$$
p_{g,k}^{\rm can}
=
\frac{\rho_{g,k}^2}
{\sum_j\rho_{g,j}^2}.
$$

Define

$$
\boxed{
r_{g,\rm eff}^{\mathsf D\leftarrow\mathsf V}
=
\exp\left(
-\sum_kp_{g,k}^{\rm can}
\log p_{g,k}^{\rm can}
\right).
}
$$

### §78.3 Top-$m$ coupling concentration

$$
\boxed{
C_{g,m}^{\rm can}
=
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
