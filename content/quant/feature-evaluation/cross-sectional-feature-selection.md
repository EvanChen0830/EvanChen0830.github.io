---
title: "Cross-Sectional Feature Selection: Setup, Matrix Mechanics & Worked Example"
description: "The full matrix-level pipeline: within-date demeaning, horizon-band slicing, whitening, orthogonal lifting into Direction/Activity/Tail blocks, the predictive cross-covariance operator, null calibration, SVD truncation, worked toy numbers, and protecting directional exposure."
date: 2026-09-09
tags: ["feature-selection", "cross-sectional", "methodology", "quant-research", "matrix-algebra"]
math: true
aliases:
  - "/quant/feature-evaluation/part-1-setup-and-matrix-mechanics/"
  - "/quant/feature-evaluation/part-2-worked-toy-example/"
  - "/quant/feature-evaluation/part-3-reading-the-subsets/"
  - "/quant/feature-evaluation/part-4-protecting-directional-exposure/"
---

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
R_q=G_q-\Phi_{\lt q}B_{\lt q,q},
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

---

## Part C — Toy Walkthrough (end-to-end, with shapes and numbers)

Setup: $T=2$ days, $N=3$ stocks/day (balanced, for hand-calculation), $P=2$ features, one horizon band with $K=2$. $L=6$. Row order: $(A,1),(B,1),(C,1),(A,2),(B,2),(C,2)$.

**Raw inputs**

$$
X=\begin{pmatrix}1.0&0.2\\2.0&0.4\\3.0&0.3\\1.5&0.1\\2.5&0.5\\2.0&0.6\end{pmatrix},\qquad
Y_{\rm raw}=\begin{pmatrix}0.05&0.10\\0.02&0.01\\-0.01&-0.02\\0.04&0.06\\0.00&-0.01\\0.03&0.02\end{pmatrix}
$$

**B1 — demean.** Day averages: $\bar X_{d1}=(2.0,0.3)$, $\bar X_{d2}=(2.0,0.4)$.

$$
X^{XS}=\begin{pmatrix}-1.0&-0.1\\0.0&0.1\\1.0&0.0\\-0.5&-0.3\\0.5&0.1\\0.0&0.2\end{pmatrix}
$$

(each day's rows sum to zero — check: $-1.0+0.0+1.0=0$, $-0.1+0.1+0.0=0$ ✓). $Y^{XS}$ is computed the same way (omitted, same mechanics).

**B3 — whiten.** Skipping actual $\hat\Sigma_Z$ estimation for brevity, take row 1 of $Z\approx(0.03,0.07)$ as illustration.

**B4 — lift**, row 1:
$$
G_D=(0.03,\,0.07),\quad
G_V=(0.03^2{-}1,\ \sqrt2(0.03)(0.07),\ 0.07^2{-}1)=(-0.999,\,0.00297,\,-0.995),\quad
G_T=0.0762\times(0.03,0.07)
$$
Stack over all 6 rows $\to G\in\mathbb R^{6\times7}$; Gram–Schmidt $\to\Phi\in\mathbb R^{6\times7}$ (in a real run $L\gg d_\Phi$; here $L=6\lt d_\Phi=7$ is only illustrative of the mechanics, not a usable sample size).

**B5 — cross-covariance.** $C=X^{XS\top}\Phi/6\in\mathbb R^{2\times7}$. Day-decomposition check on $(f_1,\phi_{D_1})$: $C_D^{(1)}[f_1,1]=\tfrac13\big[(-1.0)(0.03)+0+(1.0)(-0.03)\big]=-0.02$, $C_D^{(2)}[f_1,1]=-0.00667$; pooled with $w_1=w_2=0.5$ gives $-0.01333$, matching the direct 6-row computation exactly — confirming (B5').

**B7–B8 — spectral & quality (illustrative numbers for a 6-feature long-band fit).**

| feature | $Q_p$ | $f_p=(D,V,T)$ |
|---|---|---|
| $\text{mom}_{20}$ | 0.0062 | (0.60, 0.02, 0.38) |
| $\text{mom}_{60}$ | 0.0135 | (0.62, 0.04, 0.34) |
| $\text{vol}_5$ | 0.0069 | (0.02, 0.94, 0.04) |
| $\text{turn}_5$ | 0.0050 | (0.00, 0.97, 0.03) |
| $\text{ma}_{10}$ | 0.0086 | (0.61, 0.02, 0.37) |
| $\text{ma}_{12}$ | 0.0074 | (0.62, 0.02, 0.36) |

**B9 — mechanism clustering.** $S^{\rm pred}(\text{mom}_{20},\text{ma}_{10})=0.98$, $S^{\rm pred}(\text{mom}_{20},\text{vol}_5)=0.08$. Cut at 0.8:
$$G_1=\{\text{mom}_{20},\text{mom}_{60},\text{ma}_{10},\text{ma}_{12}\}\ (D\text{-mechanism}),\qquad G_2=\{\text{vol}_5,\text{turn}_5\}\ (V\text{-mechanism}).$$

**B10 — redundancy families.** Within $G_1$: $S^{\rm raw}(\text{ma}_{10},\text{ma}_{12})=0.97$, $S^{\rm raw}(\text{mom}_{20},\text{ma}_{10})=0.93$, $S^{\rm raw}(\text{mom}_{20},\text{mom}_{60})=0.55$. Cut at 0.9: $R_{11}=\{\text{ma}_{10},\text{ma}_{12},\text{mom}_{20}\}$, $R_{12}=\{\text{mom}_{60}\}$. Within $G_2$: $S^{\rm raw}(\text{vol}_5,\text{turn}_5)=0.40$ → each its own family.

**B11–B12 — stability.**

| feature | $\pi_p$ | $M_p$ | $\pi_p^{XS}$ | $L_{\rm eff}$ (of $N{=}500$) |
|---|---|---|---|---|
| mom20 | 0.88 | 0.91 | 0.85 | 210 |
| mom60 | 0.90 | 0.93 | 0.82 | 190 |
| ma10 | 0.85 | 0.89 | 0.60 | **35** |
| vol5 | 0.79 | 0.86 | 0.83 | 240 |
| turn5 | 0.74 | 0.80 | 0.77 | 205 |

**B13 — selection.** In $R_{11}$, ma10's very low $L_{\rm eff}$ marks it as a few-large-caps bet; if hard-filtered, mom20 wins the family by Score. Group quota guarantees $G_2$ keeps at least one representative regardless of $G_1$'s higher raw scores.

$$
\mathcal F^{(\rm long)}_{\rm rep}=\{\text{mom}_{20},\ \text{mom}_{60},\ \text{vol}_5,\ \text{turn}_5\}
$$

Repeating B2–B13 for the short band (independent whitening, independent null block length) might give, e.g., $\mathcal F^{(\rm short)}=\{\text{turn}_5,\text{vol}_5,\text{mom}_{20}\}$ with a more $V$-tilted aggregate fingerprint.

---

## Part D — Reading the $n$ Subsets

Each $\mathcal F^{(b)}$ comes with its own audit block (Coverage, block coverage, $r_{\rm eff}$, stability). Use these to characterize what each subset is actually good for:

| Diagnostic | What it tells you |
|---|---|
| $\operatorname{Cov}_D(\mathcal F^{(b)})$ high, $\operatorname{Cov}_V,\operatorname{Cov}_T$ low | Subset is **directional**: predicts the sign and magnitude of relative outperformance. Suited to long-short position sizing — the classic "relative strength" portfolio. |
| $\operatorname{Cov}_V(\mathcal F^{(b)})$ high | Subset predicts future dispersion/co-movement structure, not who wins. Useful for risk sizing, pairs/dispersion trades, not a directional bet. |
| $\operatorname{Cov}_T(\mathcal F^{(b)})$ high | Subset predicts direction conditional on large moves only. Useful for tail-hedge timing or crash-alpha, weak signal in calm regimes. |
| $r^{(b)}$, $r_{\rm eff}^{(b)}$ across bands | How many independent mechanisms the library supports at that horizon. A band with $r_{\rm eff}\approx1$ has one alpha direction no matter how many features are in $\mathcal F^{(b)}$; compare across bands to see whether independent mechanisms emerge only at certain horizons. |
| Which $G_m$ are represented in $\mathcal F^{(b)}$ vs. which existed in $\mathcal E^{(b)}$ but got squeezed out | Tells you whether the band's budget was directional-sufficient or diluted — see Part E. |
| $L_{\rm eff}(p)$ of the features actually selected | Whether the band's directional power rests on a broad cross-section or a handful of large caps — matters directly for tradability/capacity. |

**Reading example** (continuing the toy numbers): a long-band subset with fingerprint aggregate near $(0.55,0.10,0.35)$ is directional-with-a-tail-tilt — good for a long-short book that also wants to lean into big moves; a short-band subset near $(0.05,0.85,0.10)$ is not a relative-strength book at all — it's an activity/volatility-timing signal, and using it for position sizing on direction would be a category error even though its $Q_p$ values look respectable.

**Comparing bands.** If $\operatorname{Cov}_D$ decays sharply from short to long bands while $\operatorname{Cov}_V$ or $\operatorname{Cov}_T$ holds up, it means your directional relative-strength edge decays faster than your activity/tail edge — a genuinely useful finding for holding-period design, not an artifact to explain away.

---

## Part E — Protecting Directional (Relative-Strength) Exposure in Selection

**First, confirm what "relative strength magnitude" maps to in this framework.** $\Phi_D=z$ is the whitened, within-date-demeaned, *continuous* future return — not a rank transform. So $Q_p^D$ already measures a feature's ability to predict the **size** of relative over/underperformance, not merely its sign. $\operatorname{Cov}_D(\mathcal F)$ at the set level is therefore the single number that most directly answers "how much of the cross-sectional relative-strength signal did I keep." Treat it as the primary KPI in every report (Section 26/31 of the main document), not merely one of three equally-weighted blocks.

The risk you're worried about is specifically in **Stage 2** (B13c): the block-normalized greedy fill was designed to *stop* Direction from crowding out Activity/Tail. That is exactly the mechanism to relax when Direction is what you actually want to protect. Group quota (the floor in B13c) is **not** the culprit — it only guarantees a minimum for $G_2,G_3,\dots$; it never caps how many representatives $G_1$ (the D-mechanism group) can keep.

Concrete levers, roughly in order of how directly they target the problem:

**1. Lexicographic budget allocation.** Replace the single water-filling pass with two passes: first spend budget purely on maximizing $\Delta\mathrm{Cov}_D(p)$ among D-eligible candidates until $\operatorname{Cov}_D(\mathcal F)\ge\gamma_D$ (e.g. 90–95%) or the sub-budget is exhausted; only then spend remaining slots on $V,T$ with the original water-filling logic. This guarantees Direction coverage is never sacrificed to satisfy fairness toward the other two channels — fairness only governs the leftover budget.

**2. Direction-specific scoring track.** Within a D-mechanism group, rank candidates by $\bar Q_p^D$ (block-specific quality) rather than total $\bar Q_p$:
$$
\mathrm{Score}_p^D=\bar Q_p^D\cdot\pi_p^{w_\pi}\cdot M_p^{w_M}\cdot g(E_p).
$$
This prevents a feature with strong but incidental $V$/$T$ content from outscoring a purer directional feature purely because its *total* $Q_p$ happens to be larger.

**3. Loosen the redundancy threshold inside D-mechanism groups.** Raising $\tau_{\rm raw}$ specifically for families formed within $G_1$ (e.g., require $S^{\rm raw}\ge0.95$ instead of $0.90$ before two features merge into one family) produces *more*, smaller families — hence more surviving representatives after Stage 1. This is the direct lever for "don't prune my directional features so aggressively," distinct from and complementary to lever 1.

**4. Turn $L_{\rm eff}$ into a soft discount, not a hard filter, when direction is the priority.** A hard cutoff $L_{\rm eff}(p)>L_{\min}$ can throw away a concentrated-but-genuinely-strong directional feature even when nothing else in its family comes close. Prefer
$$
g_L(p)=\min\!\left(1,\ L_{\rm eff}(p)/L_0\right),
$$
multiplied into $\mathrm{Score}_p$ (or $\mathrm{Score}_p^D$), so concentration discounts a candidate proportionally rather than disqualifying it outright — unless your capacity constraints genuinely require a hard floor, in which case keep it, but set $L_{\min}$ relative to your own tradable universe size rather than an arbitrary global number.

**5. Report $\operatorname{Cov}_D$ per feature-slot spent.** $\operatorname{Cov}_D(\mathcal F)/|\mathcal F|$ (or the marginal $\Delta\mathrm{Cov}_D$ of the last few features admitted) tells you whether your budget is still buying meaningful directional coverage or whether you're spending slots on diminishing returns — useful for deciding how big $P_{\max}$ needs to be before Direction saturates.

**Recommended default combination:** run lever 1 (lexicographic budget split) as the structural fix, lever 3 (looser $\tau_{\rm raw}$ inside D groups) so Stage 1 doesn't over-consolidate before Stage 2 even runs, and lever 4 (soft $L_{\rm eff}$ discount) so a single concentration diagnostic can't unilaterally veto your strongest directional candidate. Use lever 2 as a secondary tie-break inside $G_1$ rather than replacing $\mathrm{Score}_p$ everywhere — total-quality scoring is still the right criterion for $G_2,G_3,\dots$, since there Direction *isn't* what you're optimizing for.
