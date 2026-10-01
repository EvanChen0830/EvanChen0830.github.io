---
title: "Multi-Sub-Universe Directional Coupling Geometry, Part 6: Reports, Implementation Order, and Validation"
description: "The complete universe report schema (every candidate, not just winners), a 38-step implementation checklist from reference basis to local refit, the fit/transform boundary table, numerical identity tests, ten synthetic calibration tests, and default operating parameters."
date: 2026-09-14
series: ["Multi-Sub-Universe Directional Coupling Geometry"]
part: 6
tags: ["universe", "cross-sectional", "methodology", "quant-research", "reproducibility"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Multi-Sub-Universe Directional Coupling Geometry" >}}

## §91. Universe report format

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

## §92. Complete implementation order

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
 6. Phi_D := Z; build Phi_V, Phi_L; verify Phi^T W_U Phi = I.
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

## §93. The fit/transform boundary table

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

## §94. Numerical identity tests

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

## §95. Synthetic calibration tests

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

## §96. Default operating parameters

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

All of these thresholds are operating/screening parameters, not statistical criteria. Formal significance comes from the §82 null; stability comes from §83-84; actual feature learnability comes from the downstream pipeline starting at §87.
