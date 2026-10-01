---
title: "Cross-Sectional Feature Selection, Part 1: Setup & Matrix Mechanics"
description: "The full matrix-level pipeline: within-date demeaning, horizon-band slicing, whitening, orthogonal lifting into Direction/Activity/Tail blocks, the predictive cross-covariance operator, null calibration, SVD truncation, and the two-stage selection rule."
date: 2026-09-09
series: ["Cross-Sectional Feature Selection"]
part: 1
tags: ["feature-selection", "cross-sectional", "methodology", "quant-research", "matrix-algebra"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Feature Selection" >}}

## Part A — Setup Recap

You have:

- $X\in\mathbb R^{L\times P}$ — feature matrix, long format. Each row is an $(i,d)$ pair (asset $i$ on date $d$); $L=\sum_{d=1}^TN_d$.
- $Y_{\rm raw}\in\mathbb R^{L\times H}$ — raw multi-horizon target (not de-meaned).
- $Y^{XS}_{\rm raw}\in\mathbb R^{L\times H}$ — the same target, already stripped of the market/common component.
- A date index $d(\cdot):\{1,\dots,L\}\to\{1,\dots,T\}$ recording which row belongs to which date.

Goal: produce $n=B$ subsets $\mathcal F^{(1)},\dots,\mathcal F^{(B)}$, one per horizon band, each of which is good at predicting **cross-sectional relative strength** — and to make sure the selection procedure doesn't accidentally dilute the directional (D-block) content that portfolio construction actually depends on.

Everything below is **two-dimensional matrix algebra**. No tensor is ever required; the asset and date axes are pre-flattened into $L$, and "which day" is carried by auxiliary matrices, never by a third array axis.

---

## Part B — End-to-End Matrix Implementation

### B0. Data layout

Fix a row order once (e.g., sort by date, then by asset id within date) and keep $d(\cdot)$ as a plain integer array of length $L$. Every operation below either (a) treats $X,Y,\Phi,\dots$ as ordinary $L\times(\cdot)$ matrices, or (b) uses an auxiliary $L\times L$ or $L\times T$ matrix built from $d(\cdot)$ to encode grouping. There is no step that requires $N_d$ to be constant across days.

### B1. Within-date demeaning

Build the date indicator (one-hot) matrix

$$
\mathbf D\in\{0,1\}^{L\times T},\qquad \mathbf D_{i,d}=\mathbf 1[d(i)=d].
$$

$\mathbf D$ is sparse (exactly one 1 per row) — in practice this is implemented as a `groupby(date).transform('mean')`, not a dense matrix multiply, but the two are algebraically identical:

$$
X^{XS}=X-\mathbf D(\mathbf D^\top\mathbf D)^{-1}\mathbf D^\top X,
\qquad
Y^{XS}=Y_{\rm raw}-\mathbf D(\mathbf D^\top\mathbf D)^{-1}\mathbf D^\top Y_{\rm raw}.
\tag{B1}
$$

$(\mathbf D^\top\mathbf D)^{-1}\mathbf D^\top X\in\mathbb R^{T\times P}$ is just the matrix of per-date column means; $\mathbf D$ broadcasts each date's mean back down to its rows. If you already have $Y^{XS}_{\rm raw}$ from upstream, skip that half — but **always** apply this to $X$. This is not optional: without it, the Activity ($V$) and Tail ($T$) blocks of $\Phi$ are not protected from market-level contamination (see B4 and the note at the end of this section).

**Why demeaning $X$ alone is enough.** For any column of $\Phi$ that has a day-level constant component $\gamma_d$ (e.g., "today is a high-vol day," embedded identically for every stock that day):

$$
\sum_{i\in d}X^{XS}_{i,p}\cdot\gamma_d=\gamma_d\underbrace{\sum_{i\in d}X^{XS}_{i,p}}_{=0}=0.
$$

The day-level component cancels automatically once $X$ is demeaned, regardless of whether $\Phi$ itself is demeaned. You get this protection for free without touching $\Phi_V,\Phi_T$.

### B2. Horizon-band slicing

Choose bands $H_1,\dots,H_B$ partitioning $\{1,\dots,H\}$ (need not be equal width — pick boundaries where the label's autocorrelation structure changes, e.g., 1–5 / 6–20 / 21–60 days). For each band:

$$
Y^{(b)}=Y^{XS}[:,H_b]\in\mathbb R^{L\times K},\qquad K=|H_b|.
\tag{B2}
$$

Everything from here on is run **once per band**, independently, with its own $K$.

### B3. Whitening

$$
\hat\Sigma_Z=\frac1L Y^{(b)\top}Y^{(b)}\ (\text{shrunk, e.g. Ledoit–Wolf}),
\qquad
Z=Y^{(b)}\big(\hat\Sigma_Z^{-1/2}\big)^\top\in\mathbb R^{L\times K}.
\tag{B3}
$$

Ordinary $(L\times K)(K\times K)$ matrix product. Because $Y^{(b)}$ inherits within-date zero-sum from B1, and $\hat\Sigma_Z^{-1/2}$ is applied identically to every row, $Z$ also has within-date zero-sum:

$$
\sum_{i\in d}Z_{i,:}=\Big(\hat\Sigma_Z^{-1/2}\Big)^{\!\top}\!\!\sum_{i\in d}Y^{(b)}_{i,:}=0\quad\forall d.
\tag{B3'}
$$

This is what later makes $\Phi_D$ immune to market contamination by construction.

### B4. Orthogonal lifting — row-wise, not a matrix product

For row $t$ (i.e., one $(i,d)$ pair), with $z_t\in\mathbb R^K$:

$$
G_D(z_t)=z_t\in\mathbb R^K,
\qquad
G_V(z_t)=\operatorname{svec}(z_tz_t^\top-I_K)\in\mathbb R^{K(K+1)/2},
\qquad
G_T(z_t)=\|z_t\|_2\,z_t\in\mathbb R^K.
\tag{B4a}
$$

This is a fixed nonlinear map applied identically and independently to each row — vectorize with broadcasting (e.g. `np.einsum('ti,tj->tij', z, z)` then take the upper triangle), never a for-loop, but conceptually it touches one row at a time. Stack: $G=[G_D\mid G_V\mid G_T]\in\mathbb R^{L\times d_\Phi^{\rm raw}}$, $d_\Phi^{\rm raw}=K+\tfrac{K(K+1)}2+K$.

Then Gram–Schmidt, in order $D\to V\to T$, is ordinary 2D algebra:

$$
B_{aq}=\frac1L\Phi_a^\top G_q\in\mathbb R^{d_a\times d_q},
\qquad
R_q=G_q-\Phi_{<q}B_{<q,q},
\qquad
\Phi_q=R_q\Big(\tfrac1LR_q^\top R_q\Big)_{+}^{-1/2}.
\tag{B4b}
$$

$\Phi_a^\top G_q$ is a $(d_a\times L)(L\times d_q)$ product; $R_q^\top R_q$ is a small $d_q\times d_q$ matrix whose eigendecomposition gives the inverse-square-root. Result: $\Phi=[\Phi_D\mid\Phi_V\mid\Phi_T]\in\mathbb R^{L\times d_\Phi}$, $\tfrac1L\Phi^\top\Phi=I_{d_\Phi}$.

### B5. The predictive cross-covariance operator, and its day-by-day decomposition

$$
C=\frac1LX^{XS\top}\Phi\in\mathbb R^{P\times d_\Phi}.
\tag{B5}
$$

This single $(P\times L)(L\times d_\Phi)$ product **is simultaneously a time average and a cross-sectional statistic**. Splitting the sum over rows into a sum over dates, then over stocks within each date:

$$
\boxed{\,C=\sum_{d=1}^Tw_d\,C^{(d)},\qquad w_d=\frac{N_d}{L},\qquad C^{(d)}=\frac1{N_d}\sum_{i\in d}X^{XS}_{i,:}{}^\top\Phi_{i,:}\,}
\tag{B5'}
$$

$C^{(d)}\in\mathbb R^{P\times d_\Phi}$ is a **pure, single-day, cross-sectional** covariance (who's ranked where today vs. who wins later); $\sum_dw_d(\cdot)$ is the **time-pooling** that gives you enough samples to estimate and test it. Cross-section defines what the number *means*; time is only there to make the estimate stable enough to trust. This decomposition is the basis for both the null calibration below and the diagnostic in B12.

### B6. Null calibration — block permutation matrix

Represent the circular date shift as a permutation matrix $\Pi^{(b)}\in\{0,1\}^{L\times L}$ built to move entire date-blocks together (never splits a day):

$$
X^{(b)}_{\rm null}=\Pi^{(b)}X^{XS},
\qquad
C^{(b)}_{\rm null}=\frac1L X^{(b)\top}_{\rm null}\Phi.
\tag{B6}
$$

Draw $B_{\rm null}\ (\sim\!1000)$ such permutations with offset $\ge\Delta$ (bigger than horizon + feature memory). Get $s^{(b)}_{j,\rm null}$ from each. Use the **sequential deflation test**: at step $j$, project both the candidate mode and the null draws onto the space orthogonal to the first $j-1$ accepted modes, compare $s_j$ against the resulting null distribution, stop at the first non-significant $j$. Set $r=j-1$.

### B7. SVD & truncation

$$
C=USV^\top,\qquad C_{\rm sig}=U_rS_rV_r^\top=CV_rV_r^\top\in\mathbb R^{P\times d_\Phi}.
\tag{B7}
$$

Standard SVD of a small $P\times d_\Phi$ matrix. $c_p=C_{\rm sig}[p,:]$ is feature $p$'s signal row, in fixed $\Phi$-coordinates — stable across refits even if $r$, $U$, $V$ change.

### B8. Quality, fingerprint, mechanism

$$
Q_p=\|c_p\|_2^2,\qquad
Q_p^q=\|c_p^q\|_2^2\ (q\in\{D,V,T\}),\qquad
f_p=\frac{(Q_p^D,Q_p^V,Q_p^T)}{Q_p},\qquad
\hat c_p=\frac{c_p}{\|c_p\|_2}.
\tag{B8}
$$

All scalar/vector arithmetic on rows already extracted — no new matrix shapes.

### B9. Eligibility & mechanism clustering

Null-based quality floor $Q_{\min}$ (quantile of $\{Q^{(b)}_{p,\rm null}\}$ from B6, projected through $V_r$); $\mathcal E=\{p:Q_p>Q_{\min}\}$.

Mechanism similarity is a Gram matrix on unit vectors:

$$
\hat C=[\hat c_p]_{p\in\mathcal E}\in\mathbb R^{|\mathcal E|\times d_\Phi},
\qquad
S^{\rm pred}=|\hat C\hat C^\top|\ (\text{elementwise abs}),
\qquad
D^{\rm pred}=1-S^{\rm pred}.
\tag{B9}
$$

Agglomerative clustering (average linkage) on $D^{\rm pred}$, cut at $\tau_{\rm pred}$, gives groups $G_1,\dots,G_M$.

### B10. Raw redundancy clustering

Within each $G_m$, build the Spearman correlation Gram matrix on the *raw* (standardized) features, $S^{\rm raw}_{pq}=|\bar\rho_{pq}|$ (B1's per-date rank correlation averaged over $d$), $D^{\rm raw}=1-S^{\rm raw}$. Complete-linkage clustering, cut at $\tau_{\rm raw}$, gives families $R_{m1},R_{m2},\dots$

### B11. Temporal stability — block bootstrap matrix

A resample is a block-selection matrix $\mathbf S^{(b)}\in\{0,1,2,\dots\}^{L\times L}$ built by drawing dates with replacement and stacking their full day-blocks (row sums need not be 1). Apply the **same** $\mathbf S^{(b)}$ to $X^{XS}$ and $Y^{(b)}$ to preserve pairing:

$$
X^{(b)}=\mathbf S^{(b)}X^{XS},\qquad Y^{(b)}_{\rm rs}=\mathbf S^{(b)}Y^{(b)}.
\tag{B11}
$$

Refit B3–B8 on each ($\Phi$'s lifting coefficients frozen from the reference fit) to get $Q_p^{(b)},\hat c_p^{(b)}$, hence $\bar Q_p,\mathrm{CV}_p,\pi_p,M_p$.

### B12. Cross-sectional stability and concentration ($L_{\rm eff}$)

**Stratified asset bootstrap**: build $\mathbf S^{(b)}_{XS}$ by resampling *within* each date (by sector/cap bucket) rather than resampling dates — gives $\pi_p^{XS},M_p^{XS}$, same formulas as B11, different resampling axis.

**Concentration diagnostic** (no refit needed — computed directly from B5/B8): each row's signed contribution to feature $p$'s own mechanism direction is

$$
\kappa_{ip}=\frac{X^{XS}_{i,p}\,\big\langle\Phi_{i,:},\hat c_p\big\rangle}{\sum_{i'=1}^LX^{XS}_{i',p}\,\big\langle\Phi_{i',:},\hat c_p\big\rangle},
\qquad
L_{\rm eff}(p)=\frac1{\sum_i\kappa_{ip}^2}.
\tag{B12}
$$

This is an inverse-Herfindahl "effective number of observations" behind $c_p$ — the same construction as $r_{\rm eff}$ (ERANK) in the spectral audit, but at the observation level rather than the mode level. $L_{\rm eff}(p)\ll L/T$ (i.e., far fewer than a typical day's stock count) means the signal is carried by a handful of names, not a broad cross-section.

### B13. Selection — Stage 1 + Stage 2

**Score** (per candidate, within a family):
$$
\mathrm{Score}_p=\bar Q_p\cdot\pi_p^{w_\pi}\cdot M_p^{w_M}\cdot g(E_p),\qquad g(E)=\tfrac12[1+\tanh(E/E_0)].
\tag{B13a}
$$

**Stage 1** — one representative per family, closed-form, no search:
$$
p^*_{mk}=\arg\max_{p\in R_{mk}}\mathrm{Score}_p,\qquad \mathcal F_{\rm rep}=\{p^*_{mk}\}_{m,k}.
\tag{B13b}
$$

**Stage 2** — group-quota floor, then greedy fill on block-normalized marginal coverage gain:
$$
\mathcal F\supseteq\Big\{\arg\max_{p\in G_m}\mathrm{Score}_p:m=1,\dots,M\Big\},
\qquad
\Delta\tilde{\mathrm{Cov}}(p)=\sum_{q\in\{D,V,T\}}w_q\cdot\frac{\Delta\mathrm{Cov}_q(p)}{1-\mathrm{Cov}_q(\mathcal F)+\epsilon}.
\tag{B13c}
$$

Add candidates by descending $\mathrm{Score}_p$, accept if $\Delta\tilde{\mathrm{Cov}}(p)>\delta$, stop at $P_{\max}$ or exhaustion. Output: $\mathcal F^{(b)}$ for band $b$.

