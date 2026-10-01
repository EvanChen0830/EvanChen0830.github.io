---
title: "Cross-Sectional Predictive Geometry, Part 5: Significant Rank, the Discovery/Evaluation Split, and Feature Quality"
description: "Permutation parallel analysis to find how many modes are real, why comparing the j-th observed singular value against the null's j-th (not its first) matters, the selective-inference problem that forces a discovery/evaluation split, and the precise, three-qualifier meaning of the quality score Q."
date: 2026-09-14
series: ["Cross-Sectional Predictive Geometry"]
part: 5
tags: ["cross-sectional", "methodology", "quant-research", "predictive-geometry", "statistical-significance"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Predictive Geometry" >}}

## §14. Significant Rank $r_{\rm sig}$

### Method: Permutation Parallel Analysis

For each $b$, compute the null matrix's **full singular-value spectrum**
$s_1(C_b^{\rm null})\ge s_2(C_b^{\rm null})\ge\cdots$. For each $j$:

$$\boxed{p_j=\frac{1+\#\{b:\ s_j(C_b^{\rm null})\ge s_j(C)\}}{B_{\rm null}+1}.}$$

Accept step-down: starting at $j=1$, keep accepting as long as
$p_j<\gamma_{\rm mode}$ (default 0.05); stop at the first failure.

$$\boxed{r_{\rm sig}=\#\text{accepted modes}.}$$

(Adding 1 to both numerator and denominator is the standard Monte Carlo
p-value correction, guaranteeing $p\ge1/(B_{\rm null}+1)>0$ so the test is
valid at a finite $B$.)

### Why Compare "$j$-th Against $j$-th"

An intuitively natural but actually biased alternative is **sequential
deflation**: after accepting the first $j-1$ modes, let

$$R_j=P_U^\perp CP_V^\perp,\qquad P_U^\perp=I-U_{<j}U_{<j}^\top,\quad P_V^\perp=I-V_{<j}V_{<j}^\top,$$

and use $\|R_j\|_2$ as the observed statistic, with $\|P_U^\perp
C_b^{\rm null}P_V^\perp\|_2$ as the null.

The problem is an asymmetry between the two sides. For the observed matrix,
$P_U^\perp CP_V^\perp=\sum_{k\ge j}s_ku_kv_k^\top$, so $\|R_j\|_2$ is
**exactly equal to** $C$'s $j$-th singular value — the deflation precisely
removes the first $j-1$ principal directions. But for the null matrix, the
projector is determined by the *observed* data and has nothing to do with
$C_b^{\rm null}$'s own principal directions; removing $j-1$ randomly
oriented directions from a $P$-dimensional space barely affects an
isotropic random matrix's largest singular value. So

$$\|P_U^\perp C_b^{\rm null}P_V^\perp\|_2\approx s_1(C_b^{\rm null})\quad\text{for all }j.$$

The test is then effectively comparing the observed **$j$-th** singular
value against the null's **first**. A random matrix's singular-value
spectrum is fairly flat (typical ratio $s_4/s_1\approx0.8$), so this
comparison gets more and more unfair as $j$ grows, and $r_{\rm sig}$ is
systematically underestimated.

Concretely: suppose the null spectrum is $(1.00,0.92,0.86,0.81,\dots)$ and
the observed spectrum is $(3.0,1.5,1.05,0.95,\dots)$. Sequential deflation
gives $r_{\rm sig}=2$ (the third, $1.05$ vs. $1.00$, is marginal; the
fourth, $0.95$ vs. $1.00$, fails). Parallel analysis correctly compares the
fourth observed value, $0.95$, against the null's fourth, $0.81$, giving
$r_{\rm sig}=4$.

### A Note

Step-down "stop at the first failure" does not strictly control overall
FWER. This is acceptable in practice, since $r_{\rm sig}$ is a hyperparameter
for downstream analysis, not a final conclusion — but it should not be
cited as rigorous multiple-testing control. If $r_{\rm sig}$ lands at 1 or
exactly at $d_\Phi$, treat both as warning signs: the former means the
entire library illuminates only one direction (possibly built too
homogeneously), the latter means null calibration is likely wrong (all
$d_\Phi$ directions being significant on real data is extremely rare).

---

## §15. The Significant Subspace and the Discovery/Evaluation Split

### §15.1 Why a Split Is Needed

$V_r$ is chosen from **observed data** as the "strongest directions." If the
same observed data is then used to compute

$$Q_p=\|C_{p,:}V_rV_r^\top\|^2,$$

and tested against a null distribution obtained by "projecting the null
onto that same $V_r$," the two sides are treated unequally:

$$\text{observed}:\ \text{data}\to\text{pick best subspace}\to Q,\qquad \text{null}:\ \text{null data}\to\text{project onto someone else's best subspace}\to Q.$$

The null never gets to go through "picking," so it's systematically too
small, and p-values are optimistically biased. This is a textbook selective
inference problem, and the consequence isn't just theoretical tidiness: at
$d_\Phi=18$, $r_{\rm sig}=3$, the selection effect can underestimate a
marginal feature's p-value by several times.

There are two correct ways to handle this. This book uses the first as its
primary specification.

### §15.2 Primary Specification: the Discovery/Evaluation Split

Split the research period in time, with a purge gap of $\Delta=H_{\max}+
\ell_X$ in between:

$$\underbrace{\mathcal D}_{\text{discovery,}\ \approx60\%}\ \ \|\ \Delta\ \|\ \ \underbrace{\mathcal E}_{\text{evaluation,}\ \approx40\%}$$

**Decide the coordinate system and subspace on $\mathcal D$:**

1. Fit the label basis ($B_Y,\Theta,\Gamma$) per §7–§9.
2. Compute $C_{\mathcal D}$, and use §13's null and §14's parallel analysis
   to determine $r_{\rm sig}$ and

$$\boxed{V_r=[v_1,\dots,v_{r_{\rm sig}}]\quad\text{(determined by }\mathcal D\text{, fixed from here on)}.}$$

**Do all inference on $\mathcal E$:**

3. Apply the label basis to $\mathcal E$ per §10's transform rule and
   compute $C_{\mathcal E}$.
4. Define

$$\boxed{C_{\rm sig}=C_{\mathcal E}V_rV_r^\top,\qquad c_p=C_{\rm sig}[p,:].}$$

5. §13's permutation runs **only on $\mathcal E$**, with $V_r$ held fixed.
6. The bootstraps of §26–§29, the similarity-threshold calibration of §23,
   and every screening decision in §34–§39 all run on $\mathcal E$.

**Why this is valid.** $V_r$ is simply a matrix fixed outside of
$\mathcal E$'s data, so under the null hypothesis "no association between
$X$ and $\Phi$ on $\mathcal E$," the permutation distribution of
$C_{\mathcal E}V_rV_r^\top$ is exact (limited only by residual temporal
dependence between $\mathcal D$ and $\mathcal E$, which the purge gap
already handles). The selection effect is entirely isolated inside
$\mathcal D$.

**A side benefit.** The split also produces a free diagnostic: the
$\delta_{\rm drift}$ (§10) between $\mathcal D$ and $\mathcal E$. If it's
large, the coordinate system itself is already unstable between the two
segments, and any OOS conclusion should be downgraded first.

**Split ratio.** Too short a $\mathcal D$ gives an unstable $V_r$ estimate
and underestimates $r_{\rm sig}$; too short an $\mathcal E$ reduces every
p-value's resolution. $60/40$ is a reasonable starting point. The criterion:
$\mathcal D$ must be long enough for $\Sigma_Y$'s spectrum to be stable
(§7.2's bootstrap already provides this check); $\mathcal E$ must be long
enough that $T^{\mathsf D}_{\rm eff}(\mathcal E)\ge200$.

### §15.3 Alternative Specification: Rerun Selection Inside Each Permutation

When sample length doesn't support a split, instead:

> For every null resample $b$, compute the SVD of $C_b^{\rm null}$
> **itself**, take its first $r_{\rm sig}$ right-singular vectors
> $V_r^{(b)}$, and define
> $$Q_{p,b}^{\rm null}=\big\|C^{\rm null}_{b,p,:}V_r^{(b)}V_r^{(b)\top}\big\|^2.$$

I.e. every null resample fully re-enacts the "pick the best subspace" step.
$r_{\rm sig}$ is held fixed at the observed value (re-estimating $r_{\rm
sig}$ too would make the test excessively conservative and hard to
interpret).

The compute cost is negligible: $C_b^{\rm null}$ has to be computed anyway,
and an SVD of a $400\times18$ matrix is essentially free.

The two specifications differ in what they calibrate against: the split
moves the selection effect entirely out of inference; the rerun approach
folds the selection effect into the null, making both sides equal. Both are
valid, but **cannot be mixed** — every p-value in a single report must come
from the same specification, noted in the report header.

### §15.4 Final Refit

Once every screening decision is locked in, the label basis and $V_r$ can be
refit on $\mathcal D\cup\mathcal E$ for production use.

> **Rule 15.1**: no new decisions may be made after refitting. No
> re-screening features, no re-adjusting thresholds, no adding a feature
> back because it became significant after the refit.

Violating this rule silently invalidates the entire split design, in a way
that doesn't show up in the report.

---

## §16. Feature Quality

### Definition

$$\boxed{Q_p=\|c_p\|^2,}\qquad Q_p^{\mathsf B}=\big\|c_p[\mathsf B]\big\|^2,\ \ \mathsf B\in\{\mathsf D,\mathsf V,\mathsf L\},$$

where $c_p[\mathsf B]$ denotes $c_p$'s sub-vector on block $\mathsf B$'s
coordinates.

### Proposition 16.1 (additivity)

$$Q_p=Q_p^{\mathsf D}+Q_p^{\mathsf V}+Q_p^{\mathsf L}.$$

**Proof**: $\{1,\dots,d_\Phi\}$ is partitioned by the three blocks into
disjoint coordinate subsets; a vector's squared norm equals the sum of
squares over disjoint subsets. $\blacksquare$

This identity is **always true**, regardless of whether $\Phi$ is
orthogonal. So the corresponding unit test in §65 can only catch slicing
bugs, not verify orthogonality; what actually verifies orthogonality is the
$\|\Phi^\top W\Phi-I\|_F$ test.

### Proposition 16.2 ($Q$ Is a Weighted $R^2$)

Let $c_p^{\rm full}=X_p^\top W\Phi$ (unprojected). Then

$$\|c_p^{\rm full}\|^2=\big\|P_\Phi X_p\big\|_W^2,\qquad P_\Phi=\Phi\Phi^\top W,$$

i.e. the $R^2$ of a weighted least-squares regression of $X_p$ on $\Phi$'s
columns (since $\|X_p\|_W^2=1$).

**Proof**: $P_\Phi$ is the $W$-inner-product orthogonal projector onto
$\operatorname{span}(\Phi)$ ($P_\Phi^2=P_\Phi$, and $\langle P_\Phi
u,v\rangle_W=\langle u,P_\Phi v\rangle_W$, verified directly from
$\Phi^\top W\Phi=I$). So $\|P_\Phi X_p\|_W^2=X_p^\top W\Phi\Phi^\top
WX_p=\|c_p^{\rm full}\|^2$. $\blacksquare$

### Corollary 16.3 (bound)

$$\boxed{0\le Q_p\le\|c_p^{\rm full}\|^2\le\|X_p\|_W^2=1.}$$

The first inequality holds because $P_V$ is an orthogonal projection, the
second because $P_\Phi$ is an orthogonal projection and §4.1's step three
guarantees $\|X_p\|_W=1$.

### Interpretation

$Q_p$'s precise meaning is:

> **the fraction of feature $p$'s cross-sectional variance that can be
> linearly explained by a future-return function lying in the
> statistically significant mechanism subspace.**

None of the three qualifiers can be dropped:

- "cross-sectional variance" — contains no time-series predictive power;
- "linearly explained" — see §1's boundary two and §4.3's transform family;
- "lying in the significant subspace" — $Q_p$ has already been projected
  through $P_V$, so a feature predicting only directions outside $V_r$ gets
  $Q_p\approx0$.

The last point has a consequence that must be kept in mind: $V_r$ is
determined by **consensus across the whole library**. A genuine but
isolated mechanism — predicted by only one feature — won't appear in $V_r$,
and that feature will be judged invalid. This is the price this method pays
to control multiple testing. The mitigation is to also report the
unprojected $Q_p^{\rm full}=\|c_p^{\rm full}\|^2$ and its permutation
p-value; if a feature's $Q^{\rm full}$ is significant while $Q$ is not, it's
a candidate for this "isolated mechanism" case and should be examined
individually rather than discarded outright.

---

## §17. What $Q$ Can and Cannot Answer

$Q_p$ can answer:

- How much predictive information does this feature have? (magnitude)
- What kind of mechanism is it? ($Q^{\mathsf D}$ vs. $Q^{\mathsf V}$ vs.
  $Q^{\mathsf L}$)
- How does it compare to other features? (directly comparable, since
  $\|X_p\|_W=1$)

$Q_p$ **cannot** answer:

- How much money this feature can make. $Q$ doesn't account for trading
  costs, capacity, turnover, or borrow availability.
- The feature's marginal value in a portfolio. $Q$ is univariate; a
  medium-$Q$ signal orthogonal to everything else can have higher portfolio
  value than a high-$Q$ signal redundant with existing positions. This is
  exactly why the coverage criterion in §38 exists.
- Whether the direction is right. $Q$ is a squared quantity with no sign.
  Direction is answered by §21's $A^{\rm tgt}$.

### Worked Example

A feature's report shows

$$Q_p=4.0\times10^{-4},\quad Q_p^{\mathsf D}=3.6\times10^{-4},\quad Q_p^{\mathsf V}=3\times10^{-5},\quad Q_p^{\mathsf L}=1\times10^{-5}.$$

Reading: total IC magnitude is about $\sqrt{4\times10^{-4}}=0.02$, with 90%
of the energy on direction. This is a clean, moderately-strong directional
signal, usable directly for long/short ranking.

Another feature:

$$Q_q=9.0\times10^{-4},\quad Q_q^{\mathsf D}=5\times10^{-5},\quad Q_q^{\mathsf V}=8.0\times10^{-4},\quad Q_q^{\mathsf L}=5\times10^{-5}.$$

Larger overall magnitude (IC around 0.03), but 89% of the energy is in the
$\mathsf V$ block. This is **not** a directional signal — it's a
volatility/activity predictor. Put into a long/short ranking, it would
produce a position that's "long stocks about to move violently" —
economically like a long-gamma position, nothing like directional alpha.
In a directional-only study this feature should be excluded by §34's
eligibility (its $Q^{\mathsf D}$ is too small), but kept in the library for
volatility research.

---

## §18. Feature-Level Significance and Multiple Testing

### §18.1 Single-Feature P-Value

For each feature and block, the permutation p-value is

$$\boxed{p_p^{\mathsf B}=\frac{1+\#\{b:\ Q^{{\mathsf B},{\rm null}}_{p,b}\ge Q_p^{\mathsf B}\}}{B_{\rm null}+1},}$$

where $Q^{\rm null}_{p,b}$ is computed per §15's specification (split:
$V_r$ fixed; rerun: $V_r^{(b)}$ self-re-estimated). The same formula
applies to $Q_p$ (total) and to §21's $(A^{\rm tgt}_p)^2$.

### §18.2 Multiple Testing: a Resolution Problem

With $P=400$ features tested simultaneously, multiple-testing control is
mandatory. There's a mechanical limit that's easy to overlook:

$$\min_b p=\frac1{B_{\rm null}+1}.$$

Using Benjamini–Hochberg at $\gamma_{\rm FDR}=0.05$, the strongest feature
needs $p\le\gamma_{\rm FDR}/P=1.25\times10^{-4}$, requiring

$$B_{\rm null}\ge\frac{P}{\gamma_{\rm FDR}}-1=7999.$$

With only $B_{\rm null}=1000$, **not even the single strongest feature in
the world could pass BH** — not because it isn't significant, but because
p-value resolution is insufficient. This failure mode shows up on the
report as "no feature is significant," and is very easily misdiagnosed as
"no signal exists."

### §18.3 Primary Approach: Westfall–Young Max-T Step-Down

The Westfall–Young procedure controls FWER, and **does not need** p-value
resolution down to $1/P$, because it compares "the observed statistic"
against "the null's cross-feature maximum" — which is itself an extreme
value distribution.

**Step one: studentize.** Different features' null scales aren't identical
(coverage, persistence $\kappa_p$, cross-sectional distribution shape all
affect it). So standardize first:

$$t_p=\frac{Q_p^{\mathsf B}-\operatorname{mean}_b Q^{\rm null}_{p,b}}{\operatorname{sd}_bQ^{\rm null}_{p,b}+\epsilon},\qquad t_{p,b}=\frac{Q^{\rm null}_{p,b}-\operatorname{mean}_{b'}Q^{\rm null}_{p,b'}}{\operatorname{sd}_{b'}Q^{\rm null}_{p,b'}+\epsilon}.$$

**Step two: step-down maximum.** Sort observed values $t_{(1)}\ge
t_{(2)}\ge\cdots\ge t_{(P)}$; for each resample $b$ compute the suffix
maximum

$$m_{b,k}=\max_{j\ge k}t_{(j),b}.$$

**Step three: adjusted p-values.**

$$\tilde p_{(k)}=\frac{1+\#\{b:\ m_{b,k}\ge t_{(k)}\}}{B_{\rm null}+1},\qquad \tilde p_{(k)}\leftarrow\max\big(\tilde p_{(k-1)},\tilde p_{(k)}\big).$$

That last step enforces monotonicity, a required part of step-down
procedures.

**Why this solves the resolution problem**: the strongest feature's
$\tilde p$ can go as low as $1/(B_{\rm null}+1)$, and the criterion is
$\tilde p<\gamma_{\rm FWER}=0.05$. So $B_{\rm null}=2000$ is already more
than enough.

**Why Westfall–Young fits this problem better than BH**: BH assumes the
test statistics are positively dependent or independent. Here $Q_p$ values
are highly interdependent — every feature shares the same $\Phi$, the same
dates, and features correlate heavily with each other — and the dependence
structure is unknown and complex. Westfall–Young **reproduces this
structure directly** via permutation, requiring no assumptions. The cost is
that it controls FWER rather than FDR, which is stricter.

### §18.4 Secondary Approach: Permutation FDR

For a large library wanting a looser screen, a permutation-estimated FDR
can be used: at threshold $\tau$,

$$\widehat{\rm FDR}(\tau)=\frac{\frac1{B_{\rm null}}\sum_b\#\{p:t_{p,b}\ge\tau\}}{\#\{p:t_p\ge\tau\}},$$

take the smallest $\tau$ satisfying $\widehat{\rm FDR}(\tau)\le\gamma_{\rm
FDR}$. This is also entirely permutation-based and doesn't rely on
dependency assumptions.

**Reporting rule**: list both $\tilde p^{\rm WY}_p$ and the selection
produced by the $\widehat{\rm FDR}$ threshold. The primary screen (§34)
uses Westfall–Young; the FDR version serves as a sensitivity analysis. If
the two give feature counts differing by more than 3×, that means the
signal is concentrated in a handful of very strong features with the rest
marginal — a structural fact worth noting in itself.

---

## §19. Mechanism Fingerprint

### §19.1 Raw Fingerprint

$$\boxed{f_p=\left(\frac{Q_p^{\mathsf D}}{Q_p},\ \frac{Q_p^{\mathsf V}}{Q_p},\ \frac{Q_p^{\mathsf L}}{Q_p}\right),\qquad f_{\mathsf D}+f_{\mathsf V}+f_{\mathsf L}=1.}$$

This is a **compositional decomposition** of the feature's predictive
energy — a literal description of which coordinates its association falls
on.

### §19.2 Why the Raw Fingerprint Cannot Be Read as the Mechanism Directly

The blocks have different dimensions: at $K=4$, $(d_{\mathsf
D},d_{\mathsf V},d_{\mathsf L})=(4,10,4)$.

Consider a feature totally unrelated to future returns. Its $c_p^{\rm
full}$ is roughly isotropic in the $d_\Phi=18$-dimensional space, so the
expected energy allocation is

$$\mathbb E[f_{\mathsf D}]\approx\frac4{18}=0.222,\qquad \mathbb E[f_{\mathsf V}]\approx\frac{10}{18}=0.556,\qquad \mathbb E[f_{\mathsf L}]\approx\frac4{18}=0.222.$$

In other words, **a pure-noise feature's raw fingerprint is naturally
$(0.22,0.56,0.22)$**. Seeing $f_{\mathsf V}=0.55$ and declaring "this
feature mainly predicts activity" is a mistaken inference — that's just
dimensional arithmetic.

More precisely, after projecting onto $V_r$, the null expectation is

$$\mathbb E_0\big[Q_p^{\mathsf B}\big]\approx\sigma^2\sum_{j\le r_{\rm sig}}\big\|v_j[\mathsf B]\big\|^2,$$

i.e. the energy share each block occupies within $V_r$ — which depends on
$V_r$'s actual composition and is generally **not** equal to the dimension
ratio. So the analytic approximation is unreliable; it must be estimated
directly by permutation.

### §19.3 Null-Adjusted Enrichment Fingerprint

$$\boxed{e_p^{\mathsf B}=\frac{Q_p^{\mathsf B}}{\operatorname{median}_b\,Q^{{\mathsf B},{\rm null}}_{p,b}+\epsilon},\qquad \tilde f_p^{\mathsf B}=\frac{e_p^{\mathsf B}}{e_p^{\mathsf D}+e_p^{\mathsf V}+e_p^{\mathsf L}}.}$$

Under the null, $\mathbb E[\tilde f_p]\approx(1/3,1/3,1/3)$, so $\tilde f$
has a dimension-independent reference point.

Note the denominator uses **that feature's own** null median, not a common
baseline across the library. This is necessary: coverage, persistence
$\kappa_p$, and cross-sectional distribution shape all affect a single
feature's null scale.

**Reporting rule**: $f_p$ and $\tilde f_p$ are shown side by side. $f_p$
answers "where is the energy" (e.g. for position attribution); $\tilde f_p$
answers "what is the mechanism" (e.g. for clustering and interpretation).
§35's clustering and §53's reading always use $\tilde f$.

### §19.4 Worked Example: How Fingerprint Changes a Decision

Two features with identical $Q=5\times10^{-4}$:

| | $Q^{\mathsf D}$ | $Q^{\mathsf V}$ | $Q^{\mathsf L}$ | $f$ | $\tilde f$ |
|---|---:|---:|---:|---|---|
| A | $1.1\times10^{-4}$ | $2.8\times10^{-4}$ | $1.1\times10^{-4}$ | $(0.22,0.56,0.22)$ | $(0.33,0.34,0.33)$ |
| B | $3.5\times10^{-4}$ | $1.0\times10^{-4}$ | $5\times10^{-5}$ | $(0.70,0.20,0.10)$ | $(0.79,0.13,0.08)$ |

A's raw fingerprint looks "activity-dominant," but $\tilde f$ shows it's at
noise level on all three blocks — its entire $Q$ is a dimensional artifact,
and it shouldn't even clear §34's significance gate. B is a genuinely
directional feature.

Looking only at $f$, the two would be sorted into different mechanism
groups and both retained. Looking at $\tilde f$ reveals A shouldn't be in
the pipeline at all.
