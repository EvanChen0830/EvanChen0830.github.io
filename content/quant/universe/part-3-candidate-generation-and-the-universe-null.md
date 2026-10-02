---
title: "Multi-Sub-Universe Directional Coupling Geometry, Part 3: Candidate Generation and the Universe Null"
description: "Ticker-level directional contribution statistics, four disciplined ways to generate candidate sub-universes without overfitting, a geometric eligibility floor, and the size-matched stratified permutation null that decides whether a candidate's geometry is actually special."
date: 2026-09-14
series: ["Multi-Sub-Universe Directional Coupling Geometry"]
part: 3
tags: ["universe", "cross-sectional", "methodology", "quant-research", "permutation-test"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Multi-Sub-Universe Directional Coupling Geometry" >}}

## §79. Ticker-level directional and coupling contribution

### Motivation

An arbitrary subset search has $2^N$ possibilities; directly optimizing $q_D$ or $Q^{D\leftarrow V}$ over that space is guaranteed to overfit. Candidate universes have to come from a pre-specified, low-complexity generator.

### §79.1 The directional contribution operator

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

### §79.2 Directional shape

To separate magnitude from shape:

$$
\boxed{
\widetilde K_i^D =
\frac{\bar K_i^D}
{\operatorname{tr}(\bar K_i^D)+\epsilon}.
}
$$

$\widetilde K_i^D$ describes *which* horizon combinations a ticker's directional energy lives in, independent of the total amount.

### §79.3 Cross contribution

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

## §80. Multi-universe candidate generation

### Principle

Candidate generation and candidate evaluation must stay separate. The generator only produces a finite, pre-specified family; statistical and OOS validation happen afterward, in §82-84.

### §80.1 D-purity / D-energy ladder

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

### §80.2 Directional-shape clustering

Using the Frobenius distance between $\widetilde K_i^D$ matrices,

$$
d_D(i,j)
=\|\widetilde K_i^D-\widetilde K_j^D\|_F,
$$

cluster with a fixed, pre-specified method (e.g. average linkage or spectral clustering). Each cluster becomes a candidate universe.

This line finds non-nested universes that differ in their D *shape*.

### §80.3 VD-shape clustering

Concatenate

$$
\operatorname{vec}(\widetilde K_i^D),
\quad
\operatorname{vec}(\bar K_i^{VD})
$$

run discovery-only PCA on the concatenation, then cluster. This produces candidate groups that share similar raw V-D cross-structure.

### §80.4 Hybrid family

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

## §81. Universe eligibility: the geometric layer

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

## §82. The universe null: is membership actually special?

### §82.1 Null hypothesis

For a fixed candidate $\mathcal S_g$, the null is:

> Under the same dates, the same subset size, and the same pre-locked liquidity/availability strata, random membership can produce $q_D$, $Q^{D\leftarrow V}$, or concentration that is equally or more extreme.

### §82.2 Size-matched stratified membership permutation

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

### §82.3 Primary test statistics

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

### §82.4 Multiple testing

Because multiple universes are evaluated at once, the primary approach uses Westfall-Young max-T: each membership permutation produces the statistic for every candidate at once, then take

$$
M^{(b)}=\max_g T_g^{(b)}.
$$

Family-wise error is then calibrated from the distribution of $M^{(b)}$.

This is especially necessary for a nested ladder, since adjacent universes are highly correlated.

### §82.5 When the candidate generator is itself derived from labels

If membership itself was generated from label-derived statistics such as $u_i^D$ or $\widetilde K_i^D$, a fixed-membership null only answers "is this already-chosen set special" — it does not account for search cost.

There are two legitimate ways to do formal inference here:

1. **Discovery/evaluation split (primary approach)**: the generator only runs on $\mathcal D_U$; membership is fixed and evaluated on $\mathcal E_U$.
2. **Full rerun null (expensive alternative)**: every permutation reruns the entire chain — ticker scores, clustering, candidate generation, and pruning — from scratch.

The baseline specification uses approach 1.
