---
title: "Label–Alpha Geometry, Part 5: Execution Order, Validation, and the Unified Representation"
description: "The complete execution checklist from raw data to model report, a fit/transform boundary table where any violated row means a leak, eleven numerical validation tests that must pass before any significance conclusion is trusted, the default parameter table, the method's single unifying object, and an honest appendix on what it cannot answer."
date: 2026-09-14
series: ["Label-Alpha Geometry"]
part: 5
tags: ["alpha-evaluation", "feature-selection", "cross-sectional", "methodology", "quant-research", "validation"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Label-Alpha Geometry" >}}

This closing post is the capstone for the entire geometric framework — both
the feature-evaluation layer
([Cross-Sectional Predictive Geometry](/quant/feature-evaluation/) /
[Cross-Sectional Feature Selection](/quant/feature-evaluation/part-1-setup-and-matrix-mechanics/))
and this alpha-evaluation series. It's reference material: an execution
checklist, a boundary table, a validation-test suite, and a default
parameter table meant to be checked against directly, not read start to
end.

## §63. Complete execution order

The following order is complete and can be used directly as an
implementation checklist. Bracketed references point to the relevant
sections.

**Stage 0: Setup and data**

```text
 0. Lock the configuration: band definition (cumulative/incremental),
    target weight a, Delta = H_max + l_X, split ratios, B_null / B_time /
    B_XS, every threshold                                    [§6,§13,§20,§66]
 1. Align by (time, ticker) key; verify the point-in-time contract  [§2]
 2. Split discovery D and evaluation E, with a purge gap of Delta   [§15.2]
 3. Feature preprocessing: winsorize -> within-date standardize
    -> global normalization so ||X_p||_W = 1                       [§4.1]
 4. Expand the transform family; record (source, transform) index  [§4.3]
 5. Handle missing values; record coverage_p and active_ratio_p    [§4.2]
 6. Within-date demean labels; cut horizon bands                   [§5,§6]
```

**Stage 1: Build the coordinate system on discovery**

```text
 7. Sigma_Y; spectral bootstrap; decide K_eff; record lambda_j/lambda_1
    and cond(Sigma_Y)                                              [§7.1,§7.2]
 8. B_Y -> Z; verify Z^T W Z = I holds exactly                      [§7.3,§7.4]
 9. Phi_D := Z (no orthogonalization)                               [§9.2]
10. G_V, G_L -> within-date centering -> Gram-Schmidt -> Phi         [§8.1,§9.3]
11. Store B_Y, Theta_V, Theta_L, Gamma_DV, Gamma_DL, Gamma_VL        [§10]
12. Compute T_eff^D, T_eff^V                                        [§3.1]
13. C_D; permutation null (common shift, matched thinning,
    intersection re-standardization)                               [§11,§13]
14. Parallel analysis -> r_sig, V_r (fixed from here on)            [§14]
15. r_eff; m_D (effective mechanism dimension)                      [§33,§23.4]
```

**Stage 2: Do all inference on evaluation**

```text
16. Transform the label basis onto E; compute delta_drift(D->E)     [§10]
17. C_E -> c_p = C_E V_r V_r^T -> Q, Q_D, Q_V, Q_L, Q_full           [§15.2,§16]
18. kappa_p; decide whether each feature uses the N1 or N2 null     [§13.4]
19. Permutation on E -> raw p-values -> Westfall-Young adjustment   [§18]
20. Fingerprint f and enrichment f_tilde                            [§19]
21. r_D = normalize(Phi_D^T W y*); A_target; eta_D;
    significance of (A_target)^2                                   [§20,§21,§22]
22. Independent-shift null -> similarity null distribution ->
    zeta(tau); S^D                                                  [§23,§24]
23. Temporal moving-block bootstrap -> pi, M, CV (block-specific)   [§26,§27]
24. Cross-sectional stratified bootstrap (re-center after each
    resample) -> pi^XS, M^XS                                        [§29]
25. Breadth and cancellation; relativize to b_tilde, R_tilde         [§30,§31,§32]
```

**Stage 3: Selection**

```text
26. Eligibility E1-E5                                                [§34]
27. Mechanism grouping (S^D, average linkage, tau=0.90);
    report zeta and m_D                                              [§35]
28. Raw redundancy families (complete linkage);
    apply the transform-family merge rule                            [§36,§4.3]
29. Score_p = median_b (A_target^(b))^2 -> family representatives    [§37]
30. Greedily maximize Cov^feat_D to 0.90 or budget; then fill V, L   [§39]
31. Record Cov^span as a diagnostic                                  [§38.1]
32. Output the feature report                                        [§40]
```

**Stage 4: Model and attribution**

```text
33. Train the model (purge + embargo = Delta; the selected set must
    not contain test information)                                   [§41]
34. Within-date standardize alpha + global normalization             [§42]
35. OOS: transform Phi; record delta_drift                           [§10]
36. c_alpha^full -> inherited/novel decomposition; Q, f, f_tilde,
    nu_tilde                                                         [§43,§44]
37. A_target_alpha, eta_D_alpha                                      [§45]
38. Alpha permutation -> p_full, p_D/V/L, p_target, p_nov             [§46]
39. b_signal, b_pred; IC / ICIR                                       [§47,§48,§49]
40. Recompute c_p^OOS in the same OOS window; M_OOS_p; Q_decay        [§50]
41. rho_rank; A^D_{p,alpha} + null max quantile; 2x2 interpretation   [§51,§52,§53]
42. Walk-forward ridge surrogate -> R2_selected, R2_library,
    rho_rank_lin                                                      [§55,§56]
43. Family attribution: delta_R2 and R2_own                           [§57]
44. Residual: s2_resid, Q_resid (post-standardization), Q_resid_abs,
    nu_tilde_resid, p_resid                                           [§58]
45. Determine Case A-E                                                [§60]
46. Output the model report and feature-alpha report                  [§61,§62]
```

**Stage 5 (optional): production refit**

```text
47. Refit the label basis and V_r on D union E, for production use
    only. No new decisions may be made after refitting.               [§15.4]
```

## §64. Fit/transform boundary table

The table below lists every fitted object, what data it's fit on, and how
it's applied. Violating any row constitutes a leak.

| Object | Fit on | Applied to | Section |
|---|---|---|---|
| Winsorization quantiles | Computed per period (no cross-period fit) | The current period | §4.1 |
| Within-date $\mu,\sigma$ | Computed per period | The current period | §4.1 |
| Global normalization $\|X_p\|_W$ | $\mathcal D$ | The entire period | §4.1 |
| $\Sigma_Y$, $K_{\rm eff}$, $B_Y$ | $\mathcal D$ | The entire period | §7 |
| Within-date centering of $G_\bullet$ | Computed per period | The current period | §8.1 |
| $\Theta_{\mathsf V},\Theta_{\mathsf L},\Gamma_\bullet$ | $\mathcal D$ | The entire period | §9 |
| $r_{\rm sig}$, $V_r$ | $\mathcal D$ | The entire period | §14, §15 |
| $r_{\mathsf D}$ | $\mathcal D$ (or $\mathcal E$, must be consistent and logged) | The entire period | §20 |
| Eligibility list | $\mathcal E$ | Model input | §34 |
| Grouping / families | $\mathcal E$ | Model input | §35, §36 |
| Selected set | $\mathcal E$ | Model input | §39 |
| Model parameters | $\subseteq\mathcal D\cup\mathcal E$, strictly before OOS | OOS | §41 |
| Surrogate $\hat\beta_k,\lambda_k$ | Each fold's train | That fold's test | §55 |

### Vectors across folds cannot be averaged directly

Under the full walk-forward feature-selection variant (Rule 41.1's complete
version), each fold has its own $\Phi^{(k)}$ and $V_r^{(k)}$. In that case

$$\frac1K\sum_kc_p^{(k)}\quad\text{is meaningless}$$

— different folds have different coordinate systems, and adding
component-wise mixes unrelated coordinates together.

There are three things you can actually do:

1. **Average scalars:** $\overline{Q_p}=\frac1K\sum_kQ_p^{(k)}$,
   $\overline{|A^{\rm tgt}_p|}$, etc. Scalars are coordinate-independent and
   can be averaged.
2. **Average angles:** $\overline{M_p}=\frac1K\sum_k\cos\big(c_p^{(k)},c_p^{(\rm ref)}\big)$,
   where $c_p^{(\rm ref)}$ comes from a designated reference fold.
3. **Project onto a common coordinate system:** pick one fold's
   $\Phi^{(\rm ref)}$ as reference, use the §10 transform to project every
   fold's data onto it, then average. This is the only correct way to get
   an "average mechanism vector," at the cost of introducing an arbitrary
   choice of reference fold.

## §65. Numerical validation

Split into two groups: **identity tests** (verify the implementation is
correct) and **calibration tests** (verify the statistics are valid). Both
groups must pass.

### §65.1 Identity tests

| Test | Criterion |
|---|---|
| Feature demeaning | $\big|\frac1{N_d}\sum_iX^{XS}_{idp}\big|<10^{-8}$ |
| Feature normalization | $\big|\|X_p\|_W^2-1\big|<10^{-8}$ |
| Label demeaning | $\big|\sum_iY^{XS}_{id,h}\big|<10^{-8}\cdot$ scale |
| Exact whitening | $\|Z^\top WZ-I\|_F<10^{-10}$ (**must hold exactly**, see Proposition 7.1) |
| Basis orthogonality | $\|\Phi^\top W\Phi-I\|_F<10^{-6}$ (in-sample) |
| $C$ computed two ways agrees | Direct long-form vs. $\frac1T\sum_dC^{(d)}$ |
| Quality additivity | $Q_p=Q_p^{\mathsf D}+Q_p^{\mathsf V}+Q_p^{\mathsf L}$ |
| Novelty additivity | $Q^{\rm full}_s=Q^{\rm inh}_s+Q^{\rm nov}_s$ |
| Efficiency bound | $0\le\eta^{\mathsf D}\le1+\epsilon$ |
| Sign flip | $X_q=-X_p\Rightarrow S^{\rm mech}_{pq}=1,\ A^{\rm mech}_{pq}=-1$ |
| Coverage monotonicity | ${\rm Cov}^{\rm feat}(\mathcal F\cup\{p\})\ge{\rm Cov}^{\rm feat}(\mathcal F)$ |
| Coverage submodularity | Draw 100 random $\mathcal A\subset\mathcal B$ pairs, verify diminishing marginal gain (Theorem 38.3) |

Note the quality-additivity test is **tautologically true** (Proposition
16.1, a consequence of coordinate decomposition) — it can only catch
implementation slicing bugs. What actually verifies orthogonality is the
exact-whitening and basis-orthogonality tests.

### §65.2 Calibration tests

Run on synthetic data. This group is the direct evidence for whether
statistical conclusions can be trusted.

**T1 Null calibration ($Q$).** Generate independent $X$ and $Y$ (preserving
real cross-sectional correlation structure and temporal autocorrelation,
but independent of each other), run the full pipeline.

Criterion: p-values roughly uniform; false-positive rate $\le5\%$ at
$\gamma_{\rm FWER}=0.05$.

*This is the single most important test.* If the §15 split or refit
specification is implemented incorrectly (e.g. the null still uses the
observed $V_r$), this test's false-positive rate will be significantly
above 5%.

**T2 Null calibration (similarity).** On the same synthetic data, compute
$\zeta(\tau_{\mathsf D})$.

Criterion: $\zeta\approx0$ (within Monte Carlo error).

**T3 Fingerprint calibration.** $\tilde f_p$ of a pure-noise feature.

Criterion: close to $(1/3,1/3,1/3)$; average deviation of the three
components $<0.05$.

(Cross-check: the same features' raw $f_p$ should be close to the dimension
ratio $(4,10,4)/18$, which is exactly what the §19.2 analysis predicts.)

**T4 Rotation invariance.** Rerun the entire pipeline with 10 different BLAS
thread counts and random seeds.

Criterion: every quantity listed in Proposition 9.2 agrees exactly within
numerical error ($10^{-10}$ relative).

*Must be repeated multiple times.* This test's failure mode is
**intermittent** — passing once does not mean it's correct.

**T5 Novelty calibration.** Generate random $\alpha$ unrelated to $\Phi$.

Criterion: $\tilde\nu_\alpha\approx0.5$; $p^{\rm nov}_\alpha$ uniform.

(Cross-check: raw $\nu_\alpha$ should be close to $(d_\Phi-r_{\rm sig})/d_\Phi$,
verifying §43.4.)

**T6 Breadth/cancellation calibration.** Pure-noise feature.

Criterion: $\tilde b_p\approx1$, $\tilde R_p\approx1$.

(Cross-check: raw $b$ should be close to 0.405, raw $R^{\rm canc}$ should be
close to $\sqrt{\pi/2N}$.)

**T7 Target recovery.** Set $\alpha=y^*$.

Criterion: $A^{\rm tgt}_\alpha>0$ and $\eta^{\mathsf D}_\alpha\approx1$;
setting $\alpha=-y^*$ gives $A^{\rm tgt}_\alpha<0$.

**T8 Injected-signal recovery.** Construct
$X_p=\text{IC}\cdot y^*/\|y^*\|+\sqrt{1-\text{IC}^2}\cdot\text{noise}$, with
$\text{IC}\in\{0.01,0.03,0.05\}$.

Criterion: the estimated $|A^{\rm tgt}_p|$ matches the injected IC within
$\pm20\%$. This also validates §11's magnitude scale.

**T9 Dynamic universe.** Randomly vary $N_d$ within $\pm40\%$, randomly add
and remove tickers.

Criterion: the whole pipeline doesn't depend on equal row counts; T1's
calibration still holds under a dynamic universe (this validates §13.3's
matched thinning).

**T10 Leakage detection.** Shift the label forward in time by $2\Delta$
periods as a whole (destroying the real association while keeping all other
structure intact).

Criterion: every $Q$ drops to the null level. If not, there's a leak
somewhere in the pipeline.

**T11 Slow-moving feature.** Construct a persistent feature with
$\kappa_p=0.95$ that is **unrelated to the label**.

Criterion: N1's (circular shift) p-value visibly deviates from uniform
(exactly the power loss §13.4 predicts); N2's (stratified within-date
permutation) p-value is uniform. This contrast validates the necessity of
the dual-null design.

### §65.3 Acceptance criteria

> An implementation's output should not be trusted for any significance
> conclusion until **T1, T2, T3, T4, and T5** all pass. T1 validates the
> inference specification, T2 validates similarity calibration, T3 and T5
> validate dimensional calibration, T4 validates that the coordinate system
> is well-defined. Together these five cover every part of the method that
> can systematically change a conclusion.

## §66. Default parameters

| Parameter | Default | Section |
|---|---:|---|
| $X$ winsorization | 1% / 99% cross-sectional | §4.1 |
| Minimum daily feature coverage | 0.80 | §4.2 |
| Transform family | $\{x,|x|,x|x|\}$ | §4.3 |
| Horizon band convention | Non-overlapping incremental (preferred) | §6 |
| $\kappa_{\max}$ (whitening condition number ceiling) | 100 | §7.2 |
| $\tau_\lambda$ (eigenvalue bootstrap relative SE ceiling) | 0.25 | §7.2 |
| $\epsilon_{\rm rank}$ ($\operatorname{orth}_W$ numerical cutoff) | $10^{-10}$ | §9.1 |
| Discovery/evaluation ratio | 60/40 | §15.2 |
| $\Delta$ (minimum permutation shift = purge gap) | $H_{\max}+\ell_X$ | §13.2, §41 |
| $\rho_{\min}$ (minimum intersection ratio) | 0.90 | §13.3 |
| $R_{\rm thin}$ (matched thinning repeats) | 5 | §13.3 |
| $\kappa_p$ threshold (switch to N2 null) | 0.50 | §13.4 |
| $B_{\rm null}$ | 2000 (capped at all legal shifts) | §13.2, §18.3 |
| $\gamma_{\rm mode}$ (mode significance level) | 0.05 | §14 |
| Primary multiple-testing scheme | Westfall–Young, $\gamma_{\rm FWER}=0.05$ | §18.3 |
| Secondary multiple-testing scheme | Permutation FDR, $\gamma_{\rm FDR}=0.05$ | §18.4 |
| $B_{\rm time}$ (temporal bootstrap) | 1000 | §26 |
| Block length $\ell_{\rm blk}$ | $\ge H_{\max}$ | §26 |
| $B_{\rm XS}$ (cross-sectional bootstrap) | 500 | §29 |
| $\tau_{\mathsf D}$ (mechanism-grouping cut point) | 0.90 (= accepting 19% energy loss) | §23.1, §35 |
| Raw redundancy cut point | 0.90 | §36 |
| $\pi$ threshold (three terms) | 0.80 | §34 |
| $M$ threshold (two terms) | 0.60 | §34 |
| Breadth warning | $\tilde b<0.5$ / $\tilde b^{\rm p10}<0.2$ | §31.3 |
| ${\rm Cov}^{\rm feat}_{\mathsf D}$ target | 0.90 | §39 |
| ${\rm Cov}^{\rm feat}_{\mathsf V},{\rm Cov}^{\rm feat}_{\mathsf L}$ target | 0.80 | §39 |
| $m_{\mathsf D}$ warning threshold | $<4$ → grouping has no discriminative power | §23.4 |
| $\delta_{\rm drift}$ warning threshold | $>0.3$ | §10 |
| $T_{\rm eff}/T$ warning threshold | $<0.1$ | §3.1 |
| $\epsilon$ (numerical floor) | $10^{-12}$ | §0.3 |

**These thresholds are operational filtering parameters, not statistical
criteria.** Statistical significance is always determined by permutation
(§14, §18, §46). The thresholds exist to filter out features that are
"significant but unusable"; their values can be adjusted for the research
target, but **must be locked and logged before research begins** — adjusting
thresholds after the fact until results look good is a form of overfitting
this method cannot detect.

## §67. The unified representation

### One object

The entire method ultimately reduces to this: for any within-date
standardized cross-sectional signal $s$ with $\|s\|_W=1$,

$$\boxed{c^{\rm full}_s=s^\top W\Phi\ \in\mathbb R^{d_\Phi},\qquad c^{\rm inh}_s=c^{\rm full}_sP_V,\qquad c^{\rm nov}_s=c^{\rm full}_s(I-P_V).}$$

From this single object, every geometric property of the signal can be read
off:

$$\begin{aligned}
&Q_s=\|c^{\rm full}_s\|^2 && \text{strength}\\
&\big(Q^{\mathsf D}_s,Q^{\mathsf V}_s,Q^{\mathsf L}_s\big) && \text{mechanism decomposition}\\
&\tilde f_s && \text{mechanism fingerprint (dimension-calibrated)}\\
&A^{\rm tgt}_s=\big(c^{{\rm full},\mathsf D}_s\big)^\top r_{\mathsf D} && \text{direction}\\
&\eta^{\mathsf D}_s && \text{target precision}\\
&\tilde\nu_s && \text{novelty relative to the feature library}\\
&\pi_s,\ M_s && \text{temporal and cross-sectional stability}\\
&\tilde b_s,\ \tilde R_s && \text{cross-sectional support}
\end{aligned}$$

and the relationship between two signals is

$$A^{\rm mech}_{s_1,s_2}=\hat c_{s_1}^\top\hat c_{s_2},\qquad S^{\rm mech}_{s_1,s_2}=\big|A^{\rm mech}_{s_1,s_2}\big|.$$

### Three applications

**Feature layer:**

$$X_p\ \xrightarrow{\ \Phi(Y)\ }\ c_p\quad\Longrightarrow\quad\text{select features that are real, stable, and mechanism-complementary}.$$

**Alpha layer:**

$$\hat\alpha\ \xrightarrow{\ \Phi(Y)\ }\ c_\alpha\quad\Longrightarrow\quad\text{is the model predicting direction, activity, or large-move direction}.$$

**Cross layer:**

$$c_p\ \longleftrightarrow\ c^{\rm inh}_\alpha\quad\Longrightarrow\quad\text{which known mechanisms did the model inherit},$$

$$\varepsilon_\alpha=\tilde\alpha-\text{OOSLinear}(X)\ \xrightarrow{\ \Phi(Y)\ }\ c_\varepsilon\quad\Longrightarrow\quad\text{what did the model create that's genuinely new}.$$

### The cost of this unification

Putting three questions into the same coordinate system buys comparability,
but it also means **all three questions share the same assumptions:**

1. The space spanned by the three liftings is sufficient to describe the
   predictive patterns we care about (§8);
2. The mapping from signal to label is linear on the feature side, or can
   be covered by the transform family (§4.3);
3. $\Phi$ is sufficiently stable between the research period and the OOS
   period (monitored via $\delta_{\rm drift}$, §10).

If any one of the three fails, the answers to all three questions fail
**simultaneously**, and in the same way. This is the structural risk of a
unified coordinate system, and it's why the three coordinate-system health
indicators — $\delta_{\rm drift}$, the stability of $K_{\rm eff}$, and
$T_{\rm eff}$ — must travel with every single report.

What it buys in return is an ability rare in quantitative research: **the
same ruler, measuring features, measuring the model, and measuring what's
left over once the model's output has the features subtracted out.**

## Appendix A: Method boundaries

This appendix states plainly the questions this method **cannot** answer.
These aren't flaws — they're the scope of applicability, stated explicitly
so that "this method can't see it" isn't misread as "it doesn't exist."

### A.1 The label-side function family is finite

$\Phi$ is spanned by three liftings: $z$, $\operatorname{svec}(zz^\top)$,
$\|z\|z$. Other functions of future returns are not covered, for example:

- The sign of extreme tails (indicator functions like $\mathbf 1[z_h>3]$);
- Third-order-and-above structure (cross-sectional differences in skew or
  kurtosis);
- Cross-day path shape (a band only contains several horizons computed
  from the same starting point).

A signal that predicts only these will get $Q\approx0$ from this method.

The way to extend coverage is to add more liftings (e.g. a third-order
version of $\operatorname{svec}$, or quantile indicator functions), at the
cost of increasing $d_\Phi$, which dilutes the fingerprint's denominator and
reduces statistical power in every direction. This is an explicit
trade-off: **more dimensions means more patterns visible, but each one seen
less clearly.**

### A.2 The signal-side mapping is linear

$c_s=s^\top W\Phi$ is linear in $s$. If $Y=X^2+\text{noise}$ and $X$ is
symmetrically distributed, then $\mathbb E[X\cdot z]=\mathbb E[X^3]=0$, and
the raw $X$'s $Q\approx0$.

The §4.3 transform family pushes the boundary from "is $X$ itself linearly
related" to "is any member of $X$'s specified transform family linearly
related." But it's still not "any function of $X$."

So the correct statement of this method is:

> **It answers "do these representations of this feature carry predictive
> information," not "the raw variable $X$ has no information under any
> possible nonlinear transform."**

This distinction matters especially when the final model is itself
nonlinear: the model may use a transform the transform family doesn't
cover, in which case §52's mechanism-inheritance analysis will
underestimate the degree of inheritance, and §58's residual analysis will
classify it as a "new mechanism." Step 4 of Case E's verification in §60
(trying to construct a new feature that captures that direction) exists
precisely to distinguish these two situations.

### A.3 Circular-shift null power against persistent signals

Fully addressed in §13.4: diagnosed via $\kappa_p$, with stratified
within-date permutation as an alternative null. But it's important to keep
in mind the two nulls test **different hypotheses**:

- N1: "the feature's cross-sectional ranking is unrelated to the return
  mechanism at the correct time";
- N2: "after controlling for strata (industry × market cap), the
  stock-level feature ranking is unrelated to the return mechanism."

N2 is weaker: it cannot detect a purely industry-level or market-cap-level
effect (those effects are preserved within strata). So for a slow-moving
feature, the strongest conclusion this method can offer is "there is
stock-level predictive power," not "there is predictive power."

### A.4 Equal weighting is equal-row weighting, not equal-influence weighting

Fully addressed in §3: diagnosed via $T_{\rm eff}$. The practical
consequence is that $Q$ is not immune to market regime, especially for the
$\mathsf V$ and $\mathsf L$ blocks (whose contribution scales like
$\sim\|z\|^2$).

### A.5 Discriminative power of a low-dimensional mechanism space

Fully addressed in §23.3. When $m_{\mathsf D}\le3$, the question "are these
two features the same mechanism" cannot be answered by the angle between
them alone. This method's design therefore bases selection on a coverage
criterion (§38.2) rather than pairwise thresholds, but **descriptive**
conclusions like "the model inherited feature $p$" are still bound by this
limit, and only hold up when paired with §52's null quantile and §57's
linear attribution.

### A.6 Association is not causation

A large $Q_p$ only means the feature has a stable cross-sectional
association with the future return mechanism. It may come from a genuine
pricing mechanism, or from a shared exposure to an unmodeled risk factor.
This method does not distinguish the two.

To distinguish them requires adding, after §5, a preprocessing layer that
cross-sectionally regresses out exposure to known factors and keeps the
residual. This doesn't change any downstream formula, but it changes the
meaning of every conclusion: from "this feature predicts returns" to "this
feature predicts returns beyond known factors."

### A.7 $Q$ is not return

$Q$ does not account for transaction cost, capacity, turnover, borrow
availability, or market impact. A feature with high $Q$ but extremely low
$\tilde b$ and $\tilde R$ close to 1 may be entirely untradeable in
practice. The §30–§32 diagnostics exist for exactly this reason, but they
are only warnings — they cannot substitute for an actual cost model and
capacity analysis.

## Appendix B: Interpretation quick-reference

| Observation | Meaning | What to do |
|---|---|---|
| $K_{\rm eff}\lt K$ | Whitening hit an ill-conditioned case; some horizon directions are unstable | Go back to §6, use the incremental basis or merge horizons; truncation is a last resort |
| $T_{\rm eff}/T<0.1$ | Conclusions are driven by a handful of event days | Downgrade all of that block's conclusions to exploratory |
| $\delta_{\rm drift}>0.3$ | The coordinate system itself is drifting | Cross-period comparisons are invalid; rule out coordinate issues before discussing model decay |
| $r_{\rm sig}=1$ | The library illuminates only one mechanism | Add a new data source, not more features of the same type |
| $r_{\rm eff}\gg r_{\rm sig}$ | Energy is spread across weak directions | $V_r$ is unstable; downgrade downstream conclusions |
| $m_{\mathsf D}<4$ | The mechanism space is too small for angles to discriminate | Use grouping only for report organization; let the coverage criterion drive selection |
| $\zeta(\tau)\le0$ | The observed similarity structure is indistinguishable from the null | Same as above; conclusions must not be drawn at the cluster level |
| No feature is significant | Usually $B_{\rm null}$ is too small | Check whether $1/(B_{\rm null}+1)$ is small enough (§18.2) |
| $\tilde f_p\approx(1/3,1/3,1/3)$ | This feature is noise | Regardless of how structured the raw $f$ looks |
| $\tilde b_p\approx1,\tilde R_p\approx1$ | Concentration and net effect are both indistinguishable from noise | Strongly suspect $Q_p$ comes from a handful of extreme names |
| $\pi^{\rm time}$ high but $M^{\rm time}$ low | Always strong, but predicting something different each period | Hard-reject; go back and decompose at the feature-construction layer |
| $\pi^{\rm time,sign}<0.8$ | Sign flips | Reject; never flip sign on OOS |
| $\eta^{\mathsf D}_p\approx0$ but $Q^{\mathsf D}_p$ large | Directional information is pointed at a different horizon combination | Flag `TARGET_MISMATCH`, route to another band |
| ${\rm Cov}^{\rm span}_{\mathsf D}=1$ with only 3 features selected | Span coverage has already saturated | Inevitable consequence of $d_{\mathsf D}=K$; use ${\rm Cov}^{\rm feat}$ instead |
| $\tilde f_\alpha$ dominated by $\mathsf V$ | The model is predicting volatility, not direction | A long-short rank becomes a long-gamma position; change the training target |
| $\nu_\alpha=0.4$ | Below the null baseline $(d_\Phi-r_{\rm sig})/d_\Phi$ | Use $\tilde\nu$ for interpretation; 0.4 actually means "not novel" |
| $b^{\rm sig}$ high but $b^{\rm pred}$ low | Bets broadly, predictive power is concentrated | Threshold low-contribution names to cut cost with almost no return given up |
| $R^2_{\rm lin}$ low but $Q_\varepsilon$ not significant | The nonlinearity is pure cost (Case C) | Strengthen regularization or fall back to the surrogate |
| $R^2_{\rm lin,library}\gg R^2_{\rm lin,selected}$ | The model is using information outside the selected set | Check whether the model's actual inputs match `selected` |
| $\Delta R^2_g\approx0$ but $R^2_{g,\rm own}$ high | Masked by another family | Either one alone suffices, but you can't drop both |
| $A^{\mathsf D}_{p,\alpha}=0.92$ with the top five close together | Inevitable in a low-dimensional space, not inheritance | Needs the null max quantile and §57 support together |
| $\tilde\nu_\varepsilon>0.7$ and $p^{\rm nov}_\varepsilon<0.05$ | Possibly a new mechanism (Case E) | Four-step verification; the ultimate goal is building a new feature that captures it |

---

*Every definition, formula, default, and report field in this book forms
one self-consistent specification. §65's identity tests and calibration
tests are the acceptance criteria: no significance conclusion should be
trusted until all five core calibration tests (T1–T5) pass.*
