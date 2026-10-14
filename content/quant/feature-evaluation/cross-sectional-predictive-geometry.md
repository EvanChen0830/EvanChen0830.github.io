---
title: "Cross-Sectional Predictive Geometry"
description: "Lifting cross-sectional returns into Direction, Activity, and Tail blocks: rank-controlled whitening, predictive cross-covariance operators, SVD, permutation nulls, feature quality, signed target alignment, mechanism redundancy, and temporal/cross-sectional stability."
date: 2026-09-14
tags: ["cross-sectional", "methodology", "quant-research", "predictive-geometry", "feature-selection"]
math: true
aliases:
  - "/quant/feature-evaluation/geometry-part-1-notation-and-data-contract/"
  - "/quant/feature-evaluation/geometry-part-2-preprocessing-and-horizon-bands/"
  - "/quant/feature-evaluation/geometry-part-3-whitening-and-the-three-blocks/"
  - "/quant/feature-evaluation/geometry-part-4-predictive-operator-svd-and-nulls/"
  - "/quant/feature-evaluation/geometry-part-5-significant-rank-and-feature-quality/"
  - "/quant/feature-evaluation/geometry-part-6-target-alignment-and-efficiency/"
  - "/quant/feature-evaluation/geometry-part-7-mechanism-similarity-and-redundancy/"
  - "/quant/feature-evaluation/geometry-part-8-stability-and-breadth/"
---

> This series is a from-scratch formal treatment of the same underlying
> framework as the [Cross-Sectional Feature Selection](/quant/feature-evaluation/cross-sectional-feature-selection/)
> series — same Direction/Activity/Tail decomposition, same discovery/evaluation
> discipline — but derived as formal propositions with proofs rather than
> presented as an implementation guide. Where the two overlap, this series is
> the more rigorous statement; the implementation-guide series is the faster
> path to code. One terminology note up front: this book labels the third
> predictive-mechanism block $\mathsf L$ where the implementation-guide series
> calls it $T$ (Tail) — identical definition, $G_{\mathsf L}(z) = \|z\|_2\,z$
> either way. The letter changed because $T$ is already used here for *total
> number of periods*.

## The problem this book solves

In cross-sectional quantitative research, a signal's relationship to future
returns is almost always compressed into a single number: rank IC. That
number has a fundamental problem — **it blends together three kinds of
predictive power that are completely different from a trading standpoint**:

1. "Who will outperform whom" (direction);
2. "Who will move more violently than whom" (magnitude, volatility, activity);
3. "Whose direction is right when a large move actually happens" (conditional
   direction).

A feature with an IC of 0.02 could be a clean directional signal, or it
could be a volatility predictor leaking through return skew. The two require
completely different position sizing and risk budgeting, but they're
indistinguishable at the level of IC.

This method's approach: **don't treat future returns as a scalar — lift them
into a finite-dimensional geometric space**, and assign every signal a
coordinate vector in that space. With a coordinate vector, questions like
"what does this feature actually predict," "are these two features the same
thing," and "did the model produce anything a linear combination of known
features couldn't" all become geometric operations on vectors: length,
angle, projection, span.

The entire method reduces to one object. For any cross-sectional signal $s$
(a raw feature, a model's alpha, or an alpha's residual after removing a
linear component), define

$$c_s=s^\top W\Phi,$$

where $\Phi$ is an orthonormal label-mechanism basis built from future
returns, and $W$ gives every cross-section equal weight. **Everything else in
this book is about how to correctly build $\Phi$, how to tell which
components of $c_s$ are not noise, and how to read it.**

## Who this is for, and how to use it

This book assumes familiarity with linear algebra (particularly SVD and
projection operators), basic permutation tests and the bootstrap, and the
general workflow of cross-sectional quant research. It does not assume
familiarity with any particular factor model.

Chapters are ordered by logical dependency; read in order. If time is short:

| Goal | Minimum chapters |
|---|---|
| Understand the overall idea | §1, §8, §11, §16, §67 |
| Implement the label basis | §4–§10 |
| Implement significance testing | §13–§18 |
| Implement feature selection | §26–§27, §34–§39 |
| Model diagnostics only (feature side already done) | §41–§46, §54–§60 |

Every section has a fixed structure: **motivation** (why this object is
needed) → **definition** (exact specification) → **derivation or properties**
(why it's well-defined, what its bounds are) → **worked example** (concrete
numbers) → **interpretation** (how to read the number on a report).

Every numeric example in this book is calibrated to real equity-market
magnitudes. This is worth stressing on its own: the method's core statistic
$Q$ has the magnitude of "an average per-period IC squared," so a genuinely
good feature corresponds to $Q$ on the order of $10^{-3}$, not $10^{-1}$.
Getting this scale wrong by one or two orders of magnitude will mis-set
every threshold downstream, deleting all real signal or keeping all noise.

## §0. Notation

This book uses a disambiguated symbol set. The correspondence to earlier
specification versions is listed at the end of this section, for comparison
against existing implementation code.

### §0.1 Indices and dimensions

| Symbol | Meaning | Range |
|---|---|---|
| $d$ | cross-sectional observation time ("day" by convention, can be any frequency) | $1,\dots,T$ |
| $i,j$ | ticker / instrument | $i\in\mathcal U_d$ |
| $p,q$ | feature index | $1,\dots,P$ |
| $h$ | horizon index | $1,\dots,H$ |
| $\beta$ | horizon band index | $1,\dots,B_{\rm band}$ |
| $b$ | permutation / bootstrap resample index | $1,\dots,B_\bullet$ |
| $k$ | feature transform index (§4.3) | $1,\dots,|\mathcal T_p|$ |
| $\mathcal U_d$ | universe at period $d$, $N_d=|\mathcal U_d|$ | $N_d$ may vary with $d$ |
| $T$ | total number of periods (**this** meaning only) | |
| $K$ | label dimension retained by the current band (becomes $K_{\rm eff}$ after §7) | |
| $L$ | total rows in long format, $\sum_dN_d$ | |
| $d_{\mathsf D},d_{\mathsf V},d_{\mathsf L}$ | dimensions of the three blocks | $K,\ K(K{+}1)/2,\ K$ |
| $d_\Phi$ | $d_{\mathsf D}+d_{\mathsf V}+d_{\mathsf L}$ (actual value after rank truncation) | |

The three predictive-mechanism blocks are always written in **sans-serif**
to distinguish them from other uses of the same letters:

$$
\mathsf D\ (\text{Direction}),\qquad
\mathsf V\ (\text{Volatility / Activity}),\qquad
\mathsf L\ (\text{Large-move directional}).
$$

### §0.2 Data, bases, and transforms

| Symbol | Meaning |
|---|---|
| $X_d\in\mathbb R^{N_d\times P}$ | raw feature panel |
| $X^{XS}$ | feature after within-date standardization and global normalization (§4) |
| $Y_d^{raw}\in\mathbb R^{N_d\times H}$ | raw future return |
| $Y^{XS}$ | label after within-date demeaning (§5) |
| $W=\operatorname{diag}(w_{id})$ | row weights, $w_{id}=1/(TN_d)$ (§3); never materialized |
| $\Sigma_Y$ | weighted label covariance (§7) |
| $E,\Lambda$ | eigendecomposition $\Sigma_Y=E\Lambda E^\top$ |
| $B_Y$ | whitening matrix (§7), $Z_d=Y_d^{(\beta)}B_Y$ |
| $Z$ | whitened label, satisfying $Z^\top WZ=I_K$ (**exactly**) |
| $G_{\mathsf D},G_{\mathsf V},G_{\mathsf L}$ | outputs of the three lifting maps (§8), all within-date demeaned |
| $\Theta_{\mathsf V},\Theta_{\mathsf L}$ | orthogonalization matrices (§9) |
| $\Gamma_{\mathsf{DV}},\Gamma_{\mathsf{DL}},\Gamma_{\mathsf{VL}}$ | Gram–Schmidt projection coefficients (§9) |
| $\Phi=[\Phi_{\mathsf D}\,|\,\Phi_{\mathsf V}\,|\,\Phi_{\mathsf L}]$ | label-mechanism basis, $\Phi^\top W\Phi=I_{d_\Phi}$ |

### §0.3 Core statistics

| Symbol | Meaning |
|---|---|
| $C=X^\top W\Phi$ | predictive operator (§11) |
| $C=USV^\top$ | SVD; $V_r$ is the significant right-singular subspace (§12, §14) |
| $P_V=V_rV_r^\top$ | projection onto the significant mechanism subspace |
| $r_{\rm sig},r_{\rm eff}$ | significant rank, effective spectral rank |
| $c_s^{\rm full}=s^\top W\Phi$ | full representation of signal $s$ (§43) |
| $c_s=c_s^{\rm full}P_V$ | representation projected onto the significant subspace |
| $Q_s=\|c_s\|^2$ | quality (energy); $Q^{\mathsf D},Q^{\mathsf V},Q^{\mathsf L}$ are the per-block components |
| $f_s,\tilde f_s$ | raw fingerprint, null-adjusted enrichment fingerprint (§19) |
| $a,\ y^*$ | intended target's horizon weights, target signal $y_d^*=Y_d^{(\beta)}a$ |
| $r_{\mathsf D}$ | the target's unit direction in $\Phi_{\mathsf D}$ coordinates (§20) |
| $A_s^{\rm tgt},\ \eta_s^{\mathsf D}$ | signed target alignment, directional target efficiency |
| $A^{\rm mech}_{pq},S^{\rm mech}_{pq}$ | signed / unsigned mechanism similarity (§23) |
| $S^{raw}_{pq}$ | raw Spearman redundancy (§25) |
| $\pi^{\rm time},\pi^{\rm XS},M^{\rm time},M^{\rm XS}$ | stability measures (§27, §29) |
| $N^{\rm eff}_{d,p},\ b_{d,p},\ \tilde b_p$ | effective number of contributors, normalized breadth, relative breadth (§30–§31) |
| $R^{\rm canc}_{d,p},\ \tilde R_p$ | cancellation ratio, relative cancellation (§32) |
| $\tilde p_\bullet$ | multiple-testing-adjusted p-value (§18) |
| $\gamma_{\rm mode},\gamma_{\rm FWER},\gamma_{\rm FDR}$ | significance levels for each testing layer |
| $\hat\alpha,\tilde\alpha$ | model's raw / standardized prediction |
| $\varepsilon_\alpha$ | OOS linear-surrogate residual (§58) |
| $\nu_s$ | novelty fraction (§43) |
| $\epsilon$ | pure numerical floor (denominator protection), default $10^{-12}$ |

### §0.4 Correspondence to earlier notation

| This book | Earlier spec | Reason for change |
|---|---|---|
| $\mathsf L$ block | $T$ block | $T$ is now used for total period count and transpose |
| $E$ (eigenvectors of $\Sigma_Y$) | $V$ | $V$ is now used for SVD right-singular vectors |
| $\mathsf V$ block | $V$ block | same reason |
| $\Theta_\bullet$ | $A_D,A_V,A_T$ | $A$ is now used for alignment |
| $\Gamma_\bullet$ | $B_{DV},B_{DT},B_{VT}$ | $B$ is now used for bootstrap resample count |
| $N^{\rm eff}_{d,p}$ | $B_{d,p}$ | same, and clearer semantics |
| $\tilde p_p$ | $q_p$ | $q$ stays reserved for a feature index |
| $\gamma_\bullet$ | $\alpha_{\rm mode}$ | $\alpha$ stays reserved for the model's prediction |
| $S^{\rm mech}$ | $S^{pred}$ | distinguished from "prediction" |

---

# Part One: Problem Setup and Data Contract

## §1. Research Objective

### Motivation

Any cross-sectional alpha research program is actually answering three
independent questions simultaneously, and in practice they're routinely
conflated:

1. **Which features actually have predictive power, and what do they
   predict?**
2. **What does the model's output (the alpha) actually predict? How does it
   relate to what the features predict?**
3. **Did the model produce anything a feature library can't see linearly?**

The third question is the hardest and the most important. A nonlinear model
will obviously produce nonlinear output, but "the output is nonlinear" and
"the content is new" are two completely different things. A model that
builds a complicated interaction out of three momentum features and still
only ends up predicting the same momentum direction has highly nonlinear
output but has discovered nothing new. This method provides exact tools for
separating the two.

### Three objects

$$
X_{id,p}\ \text{(feature)},\qquad
Y_{id,h}\ \text{(future return)},\qquad
\hat\alpha_{id}\ \text{(model prediction)}.
$$

### Three analytical questions

**Question one: Feature–Label.**

$$X\longrightarrow\Phi(Y).$$

Does each feature carry information about future returns? Which kind
(direction / activity / large-move direction)? Is it stable? Is it
redundant with other features?

**Question two: Label–Alpha.**

$$\hat\alpha\longrightarrow\Phi(Y).$$

What does the model actually predict? Is its predictive power directional?
Is its direction aimed at the horizon combination we actually want to
trade?

**Question three: Feature–Alpha.**

$$X\longleftrightarrow\hat\alpha.$$

Is the model recombining known features, or has it built a new mechanism?
This question has three layers of answer (§51–§59), corresponding to "does
the output look similar," "does the mechanism look similar," and "is there
anything beyond what's linearly reconstructible."

### A unified coordinate system

The three questions must be answered **in the same coordinate system**, or
the answers can't be compared against each other. This method's core design
is exactly that: build a fixed orthonormal basis $\Phi$ from $Y$, then
project features, alphas, and residuals onto $\Phi$:

$$X_p\longmapsto c_p,\qquad \tilde\alpha\longmapsto c_\alpha,\qquad \varepsilon_\alpha\longmapsto c_\varepsilon.$$

The three are vectors in the same space, so we can compare their lengths,
their angles, and whether one lies in the span of another.

### Boundaries of the method

This method has three boundaries that need to be stated up front. They are
not defects — they're the scope within which the tool is valid. Full
discussion is in Appendix A.

**Boundary one: the label-side function family is finite and fixed.** $\Phi$
is spanned by three lifting maps (§8). Some functions of future returns
(e.g. the sign of an extreme tail event, cross-day path shape) are not in
this space; any signal that predicts only those gets $Q\approx0$ under this
method.

**Boundary two: the signal-side map is linear.** $c_s=s^\top W\Phi$ is
linear in $s$. If $Y=X^2+\text{noise}$ and $X$ is symmetrically
distributed, then $\mathbb E[XY]=\mathbb E[X^3]=0$, so raw $X$ gets
$Q\approx0$ despite carrying complete information. This method's formal
treatment is the **transform family** of §4.3: every semantic feature is
expanded into a set of transforms, each analyzed separately and then
aggregated at the source level. This pushes the boundary from "is $X$ itself
linearly related" to "is any transform in $X$'s specified transform family
linearly related" — still not "any function of $X$," but strictly weaker.

**Boundary three: every conclusion is association, not causation.** A large
$Q_p$ only means the feature has a stable cross-sectional association with
a future-return mechanism. It could come from a genuine pricing mechanism,
or from shared exposure to some unmodeled risk factor. This method does not
distinguish between the two.

---

## §2. Data Structure

### Definition

Let $d=1,\dots,T$ index cross-sectional observation time. It can be a day,
5 minutes, an hour — any unit, as long as each $d$ represents one complete
cross-section.

The universe at period $d$ is $\mathcal U_d$, $N_d=|\mathcal U_d|$, allowed
to change over time. Every observation must be uniquely identified by
$(d,i)$. The data is

$$X_d\in\mathbb R^{N_d\times P},\qquad Y_d^{raw}\in\mathbb R^{N_d\times H},\qquad \hat\alpha_d\in\mathbb R^{N_d}.$$

### Data contract

> **Contract 2.1 (keyed alignment)**: every matrix operation must align by
> the `(time, ticker)` key. Never assume that row position on different
> dates corresponds to the same ticker.

This contract looks trivial but is the single easiest thing to break in
implementation, and the hardest break to notice. The classic failure mode:
the universe file is sorted by market cap, and the market-cap ranking
changes daily, so `X[d][0]` is AAPL today and MSFT tomorrow. Every
downstream statistic still computes cleanly and looks plausible — it's just
entirely meaningless. The Dynamic Universe test in §65 exists specifically
to catch this.

> **Contract 2.2 (point-in-time)**: $X_{id}$ may only use information
> publicly available at or before time $d$; $Y_{id,h}$ may only use returns
> after $d$. Any restatement, any field backfilled with final-version data,
> constitutes leakage.

### A sense of scale

Take US equity intraday research as an example: $d$ is a 5-minute bar, 78
bars per day, roughly $T\approx39{,}000$ over two years; $N_d$ floats
between 1200 and 3000; $P=400$; $H=10$. Then $L\approx7\times10^7$.

But this method **never needs $O(L^2)$ compute or memory**. Every
"full-sample" quantity is expressed as a weighted average of per-period
quantities — a direct consequence of the weighting design in §3. In
practice, the entire pipeline's compute bottleneck is a single pass over
the panel (computing $\{C^{(d)}\}$); everything after that — SVD,
permutation, bootstrap, clustering, selection — operates on a small
$P\times d_\Phi$ matrix (typically $400\times18$). This is why the method
is computationally so light, and why the permutation in §13 must be
implemented as a date-shift rather than a row-level permutation matrix.

---

## §3. Date Weighting: Every Cross-Section Gets Equal Weight

### Motivation

Before computing "the full-sample feature–label association," we first have
to answer: **what is the basic unit of observation?**

Two natural answers exist. One treats every $(d,i)$ as an observation, so
dates with more stocks naturally carry more weight. The other treats "one
day's cross-section" as the observation, so every day is weighted equally.

This method uses the latter. Reason: the true source of independent
repetition in cross-sectional research is **time**, not individual stocks.
3000 stocks' returns on the same day share the same market shock and the
same sector rotation — their effective independent sample size is far below
3000. Under row-level equal weighting, periods with a larger universe
(typically late-bull-market, IPO booms) would dominate the full-sample
statistic, which has nothing to do with signal quality and is purely an
artifact of data supply.

### Definition

$$\omega_d=\frac1T,\qquad\boxed{w_{id}=\frac1{TN_d}.}$$

Therefore

$$\sum_{i\in\mathcal U_d}w_{id}=\frac1T,\qquad\sum_{d,i}w_{id}=1.$$

Define the weighted inner product

$$\boxed{\langle A,B\rangle_W=A^\top WB,}\qquad W=\operatorname{diag}(w_{id}).$$

In implementation, $W$ is never materialized. Everything is computed as

$$\frac1T\sum_d\frac1{N_d}(\cdot).$$

### Worked example

Let $T=2$: day one has 3000 names, day two has 1200. Under row-level equal
weighting, day one gets $3000/4200=71.4\%$. Under this method, both days get
$50\%$. If day one happens to be a crash and the feature behaves anomalously
that day, row-level weighting amplifies that single event's influence by
1.43×.

### The exact scope of equal weighting, and the $T_{\rm eff}$ diagnostic

The weighting above guarantees **equal weight by row count**: every period
contributes equally to the sum. It does **not** guarantee **equal
influence**, because the label side is not normalized by daily volatility
(a design choice made in §5).

Concretely: $\Phi_{\mathsf D}\propto Z$, and $Z$'s scale is set by the
full-sample covariance, so $\|z_{id}\|$ is systematically larger on
high-volatility days. Since

$$C^{(d)}=\frac1{N_d}X_d^{\top}\Phi_d,$$

a day with twice the usual cross-sectional dispersion contributes roughly
twice the magnitude to $C$. It's stronger for the $\mathsf V$ and $\mathsf
L$ blocks: $G_{\mathsf V}\sim\|z\|^2$, $G_{\mathsf L}\sim\|z\|^2$, so the
contribution is **four times** the magnitude.

In other words, $Q^{\mathsf D}$ is an **average squared IC weighted by that
day's dispersion**, and $Q^{\mathsf V},Q^{\mathsf L}$ are weighted by
dispersion squared. This is deliberate: the cross-sectional dispersion of
future returns is itself information we want to analyze — normalizing it
away would make "directionality during large moves" a meaningless concept.
But precisely because of that, we need a diagnostic telling us *how many
days actually determined the result*.

> **Definition 3.1 (effective number of days)**: let
> $$\theta_d^{\mathsf B}=\frac{\frac1{N_d}\|\Phi_{\mathsf B,d}\|_F^2}{\sum_{d'}\frac1{N_{d'}}\|\Phi_{\mathsf B,d'}\|_F^2},\qquad\boxed{T^{\mathsf B}_{\rm eff}=\frac1{\sum_d(\theta_d^{\mathsf B})^2}}$$
> where $\mathsf B\in\{\mathsf D,\mathsf V,\mathsf L\}$.

This is the inverse-Herfindahl of the $\{\theta_d\}$ distribution — the
"effective number of contributing days." If $T=1000$ and $T^{\mathsf
D}_{\rm eff}=180$, the entire $Q^{\mathsf D}$ is actually determined by
roughly 180 "effective days." All permutation and bootstrap effective
sample sizes must be reinterpreted accordingly: a nominally 1000-day moving
block bootstrap really only carries about 180 days of independent
information.

$T_{\rm eff}$ must appear in the §40 report. A typical healthy range is
$T_{\rm eff}/T\in[0.3,0.7]$; $T^{\mathsf V}_{\rm eff}/T<0.1$ means every
conclusion about the $\mathsf V$ block is actually driven by a handful of
volatility events, and should be downgraded to exploratory.

### Implementation note

If the research goal is *not* to analyze large moves, and strict daily
equal-influence is required instead, an additional per-date scaling
$Z_d\leftarrow Z_d/\hat\sigma_d$ can be applied after §7, where
$\hat\sigma_d^2=\frac1{KN_d}\|Z_d\|_F^2$. The cost is that the $\mathsf L$
block's meaning degrades from "directionality during large moves" to
"directionality during *relatively* large moves," and the $\mathsf V$ block
is left with only the relative structure across horizons. Both are
reasonable designs, but one must be chosen and recorded in the §66
configuration table; the rest of this book assumes this scaling is **not**
applied.

---

# Part Two: Building the Label Mechanism Space

This part builds $\Phi$ — the coordinate system every later chapter works
in. Build order: feature preprocessing (§4) → label preprocessing (§5) →
horizon banding (§6) → whitening and rank control (§7) → the three liftings
(§8) → hierarchical orthogonalization (§9) → the fit/transform boundary
(§10).

## §4. Feature Preprocessing

### Motivation

Feature-side preprocessing serves two purposes that need to be understood
separately:

1. **Robustness**: removing obvious bad prints and numerical errors so
   downstream moment estimates aren't dominated by a single data error.
2. **Comparability**: making $Q_p=\|c_p\|^2$ directly comparable across
   features. This requires every feature to have the same scale under the
   weighted inner product. Without this step, a feature priced in basis
   points and one priced in percent would differ in $Q$ by a factor of
   $10^4$ — which has nothing to do with predictive power.

### §4.1 Definition

For each date $d$ and feature $p$, take the valid observations
$\{X_{idp}:i\in\mathcal U_d\}$.

**Step one: winsorize.** Default is cross-sectional 1% / 99%:

$$X'_{idp}=\operatorname{clip}\big(X_{idp},\ q_{0.01,d,p},\ q_{0.99,d,p}\big).$$

**Step two: within-date standardization.**

$$\mu_{dp}=\frac1{N_d}\sum_iX'_{idp},\qquad \sigma_{dp}=\sqrt{\frac1{N_d}\sum_i(X'_{idp}-\mu_{dp})^2},$$

$$\boxed{X^{XS}_{idp}=\frac{X'_{idp}-\mu_{dp}}{\sigma_{dp}+\epsilon}.}$$

**Step three: global normalization.** Missing-value imputation and coverage
differences will leave different features with different actual weighted
norms (see §4.2), so a global rescale is applied before entering §11:

$$\|X_p\|_W^2=\frac1T\sum_d\frac1{N_d}\sum_i\big(X^{XS}_{idp}\big)^2,\qquad\boxed{X_p\leftarrow\frac{X_p}{\|X_p\|_W}.}$$

After step three,

$$\boxed{X_p^\top WX_p=1\quad\text{(exactly)}}$$

holds for every $p$ — this is the precondition for $Q_p\le1$ becoming a
strict inequality in §16.

### §4.2 Missing Values

Record daily coverage

$$m_{dp}=\frac{\#\{\text{valid names}\}}{N_d}.$$

Rules:

- $m_{dp}<0.8$: the feature does not contribute on that date; the date is
  recorded in $\mathcal D_p^c$ (the absence set).
- Other missing values are filled with $0$ in standardized space
  (cross-sectionally neutral).
- The report records $\text{coverage}_p=\operatorname{median}_dm_{dp}$ and
  $\text{active\_ratio}_p=T_p/T$, where $T_p=|\mathcal D_p|$ is the number
  of valid days.

**Average over valid days, not the full period.** §11's definition is

$$C_{p,:}=\frac1{T_p}\sum_{d\in\mathcal D_p}C^{(d)}_{p,:},$$

i.e. absent days are not counted as zero in the numerator. The reason:
$Q_p$ should measure "how strong this feature is when it's available," and
"how often it's available" is a separate operational question — the two
should not be multiplied together and conflated. The cost is that different
rows of $C$ are averaged over different date sets, so technically the SVD
in §12 is not the decomposition of a single linear operator. When the
library's coverage distribution is very dispersed (e.g. half the features
have 100% coverage and half have 50%), low-coverage features should be
pooled separately and run through the full pipeline on their own.

**Coverage's hidden screening effect deserves caution.** A feature like
"20-day analyst estimate revision" might have analyst coverage for only
about 60% of the universe, and under the $m_{dp}<0.8$ rule it gets excluded
almost every day, giving $Q_p\to0$ — not because it lacks predictive power,
but because of the threshold. Features like this should go into a
"restricted universe" sub-study before entering the pipeline, rather than
compete in the same pool as full-universe features.

### §4.3 Transform Family

**Motivation.** $c_p=X_p^\top W\Phi$ is linear in $X_p$. So if a feature's
relationship to future returns is U-shaped (e.g. extreme valuation
deviations in either direction predict reversion), then
$\mathbb E[X\cdot z]\approx0$, and this method reports $Q_p\approx0$ despite
the feature carrying complete predictive information.

The fix is not to abandon linearity (that would destroy the entire
geometric architecture) — it's to **move the nonlinearity to the feature
side**: every semantic feature is expanded into a set of transforms, each
of which runs through the full pipeline as its own library entry.

**Definition.** For each source feature $p$, define the transform family

$$\boxed{\mathcal T_p=\{t_1,t_2,\dots\},\qquad X_{p,k}=\text{standardize}\big(t_k(X^{XS}_p)\big).}$$

The suggested minimal starting family has three:

$$t_1(x)=x,\qquad t_2(x)=|x|,\qquad t_3(x)=x|x|.$$

Their roles divide as: $t_1$ captures a monotonic linear relationship; $t_2$
captures "degree of deviation" (symmetric, directionless) — exactly the
linearization of a U-shaped relationship; $t_3$ is an odd function with
amplified convexity, capturing "the directional effect of extreme values is
stronger than mid-range values." Note that $t_2$'s output must be
**re-demeaned and re-standardized within-date**, or it will carry a large
day-level common component (that day's cross-sectional dispersion), which
by Lemma 4.1 contributes nothing to $C$ but wastes $X$'s variance.

When a richer function family is needed, a rank-space B-spline basis (3–5
knots) can be used instead, at the cost of each source feature producing
3–5 entries, inflating library size and multiple-testing burden
accordingly.

**Indexing and aggregation.** Library entries are indexed by the pair
$(p,k)$; all statistics ($Q,A^{\rm tgt},\pi,M$, etc.) are computed at the
entry level. Source-level aggregation is only used for reporting and final
narrative:

$$Q_p^{\rm src}=\max_kQ_{p,k},\qquad k_p^*=\arg\max_kQ_{p,k}.$$

**Selection rule.** Different transforms of the same source are naturally
highly correlated, so:

> **Rule 4.1**: all transform entries of the same source feature are
> **forced into the same family** during the raw-redundancy clustering of
> §36, unless they fall into different mechanism groups (§35) *and* each
> independently passes eligibility (§34). In other words, a source feature
> contributes at most one selected entry in principle; only when $t_1$ and
> $t_2$ clearly predict **different mechanisms** (e.g. $t_1$ predicts
> direction, $t_2$ predicts activity) can both be selected simultaneously.

This rule keeps transform-family cost under control: the library's nominal
size is multiplied by $|\mathcal T|$, but the final selected count doesn't
inflate.

### §4.4 Lemma 4.1 (the absorption property of daily demeaning)

Since $X^{XS}_d$ has zero mean each period, for any matrix $\Phi_d$,

$$\frac1{N_d}X_d^{XS\top}\Phi_d=\frac1{N_d}X_d^{XS\top}\big(\Phi_d-\mathbf 1_{N_d}\bar\Phi_d\big),$$

where $\bar\Phi_d$ is $\Phi_d$'s that-period row average.

**Proof**: $X_d^{XS\top}\mathbf 1_{N_d}=0$, so $X_d^{XS\top}\mathbf
1\bar\Phi_d=0$. $\blacksquare$

**Corollary 4.2**: $C^{(d)}$ is automatically a **that-period
cross-sectional covariance**, not an uncentered cross-moment; and any
"that-period-constant" component of $\Phi$ contributes exactly zero to $C$.

This lemma gets used in §8, §9, §13, and §29. It explains two things at
once: why the $-I_K$ term in §8.2 doesn't affect any result, and why taking
the ticker intersection under permutation **requires re-centering**
(because on the intersection $X^\top\mathbf 1\ne0$, breaking the lemma's
precondition).

---

## §5. Label Preprocessing

### Definition

For each horizon:

$$\boxed{Y^{XS}_{id,h}=Y^{raw}_{id,h}-\frac1{N_d}\sum_{j\in\mathcal U_d}Y^{raw}_{jd,h}.}$$

So $\sum_iY^{XS}_{id,h}=0$.

### Two design choices and their consequences

**Choice one: demean, don't beta-neutralize.** Demeaning removes that
period's equal-weighted market return. This fixes the entire analysis's
meaning to "**cross-sectional contrast relative to the equal-weighted
market**." This means:

- A pure market-timing signal (predicting overall market direction) gets
  $Q=0$ under this method. That's correct — it isn't a cross-sectional
  signal.
- If the universe's sector composition is systematically skewed, demeaning
  does not remove sector-level common movement. To analyze purely
  idiosyncratic predictive power, run a cross-sectional regression against
  known factor exposures after this step and take the residual. This is an
  optional additional preprocessing layer that doesn't change any formula
  downstream.

**Choice two: no daily volatility normalization.** I.e. we do not set
$Y_d\leftarrow Y_d/\hat\sigma_d$.

The consequence is that $\Sigma_Y$ (§7) is the full-sample covariance, and
$\|z_{id}\|$ is systematically larger on high-volatility days. As discussed
in §3, this makes $Q$ dispersion-weighted. This choice is deliberate:
**the cross-sectional dispersion of future returns is itself the
information the $\mathsf V$ and $\mathsf L$ blocks are meant to analyze.**
Normalizing it away first would leave the $\mathsf L$ block (directionality
during large moves) with nothing to be defined over. The cost is that $Q$
is not regime-immune, and must be read alongside $T_{\rm eff}$
(Definition 3.1).

---

## §6. Horizon Bands

### Motivation

Prediction problems at different horizons are economically different
problems. 5-minute returns are mainly governed by microstructure and
liquidity supply; 24-hour returns start to reflect information diffusion
and overnight risk premia. Mixing them into one whitening and one SVD
produces a coordinate system dominated by the highest-variance horizon,
meaningless for the others. Hence: band first, then analyze independently
within each band.

### Definition

Partition the $H$ horizons into bands $H_1,\dots,H_{B_{\rm band}}$. For
band $\beta$:

$$K=|H_\beta|,\qquad Y_d^{(\beta)}=Y_d^{XS}[:,H_\beta].$$

Everything downstream runs independently within each band.

### Overlapping horizons vs. an incremental basis

This is a specification point that must be chosen explicitly. If
$Y_{\cdot,h}$ are **cumulative** returns measured from the same starting
point, different horizons are naturally nested: under a random-walk
assumption,

$$\operatorname{corr}(r_{0\to1h},r_{0\to2h})=\sqrt{1/2}\approx0.71,\qquad\operatorname{corr}(r_{0\to1h},r_{0\to4h})=\sqrt{1/4}=0.5.$$

A $K=4$ cumulative band has a correlation matrix with condition number
roughly 10–20, and whitening amplifies the smallest-eigenvalue direction by
$\lambda_{\min}^{-1/2}$ — amplifying estimation noise right along with it.

> **Specification 6.1**: within a band, prefer **non-overlapping incremental
> returns**
> $$r_{0\to1h},\quad r_{1h\to2h},\quad r_{2h\to4h},\quad r_{4h\to8h},$$
> whose correlation matrix under a random walk is the identity, with
> condition number near 1. If the research question genuinely requires
> cumulative returns (e.g. execution holds to a $4h$ close), keep the
> cumulative form, but §7's rank control becomes mandatory rather than
> optional.
>
> Either way, the band definition (cumulative or incremental) must be
> recorded in configuration and appear in the report, because $c_p$'s
> coordinate meaning depends entirely on it.

### A sense of dimension

This book uses the following band split as its running example:

$$H_1=\{5m,10m,15m,30m\},\quad H_2=\{1h,2h,4h\},\quad H_3=\{6h,12h,24h\},$$

| band | $K$ | $d_{\mathsf D}$ | $d_{\mathsf V}=K(K{+}1)/2$ | $d_{\mathsf L}$ | $d_\Phi$ |
|---|---:|---:|---:|---:|---:|
| $H_1$ | 4 | 4 | 10 | 4 | 18 |
| $H_2$ | 3 | 3 | 6 | 3 | 12 |
| $H_3$ | 3 | 3 | 6 | 3 | 12 |

Remember this table. $d_{\mathsf V}$ grows quadratically in $K$ while
$d_{\mathsf D}$ only grows linearly: at $K=4$, the Activity block already
takes up $10/18=56\%$ of all coordinates. This fact has decisive
consequences in two later places: §19 (the fingerprint's dimensional
baseline) and §38 (coverage's rank ceiling).

---

## §7. Label Whitening and Rank Control

### Motivation

Whitening's purpose is to put different horizons' coordinates on the same
scale and make them mutually uncorrelated, so that "which horizon
combination did this feature predict" has a unique, scale-free answer.
Without whitening, the horizon with the largest variance in $\Sigma_Y$
(usually the longest) would dominate every inner product.

But whitening is also the **only step in the entire method that amplifies
noise**. Along $\Sigma_Y$'s $j$-th eigendirection, the data is multiplied by
$\lambda_j^{-1/2}$. If $\lambda_K/\lambda_1=10^{-2}$, the weakest direction
is amplified 10×, and so is the estimation error in that direction.
Whitening must therefore be paired with rank control, and the rank-control
threshold must be set by statistical stability, not an arbitrary numeric
floor.

### §7.1 Covariance and Spectrum

$$\boxed{\Sigma_Y=\frac1T\sum_{d=1}^T\frac1{N_d}Y_d^{(\beta)\top}Y_d^{(\beta)}\in\mathbb R^{K\times K}.}$$

(By the daily demeaning in §5, this genuinely is a covariance, not an
uncentered moment.) Eigendecomposition:

$$\Sigma_Y=E\Lambda E^\top,\qquad \lambda_1\ge\lambda_2\ge\cdots\ge\lambda_K>0.$$

### §7.2 Retained Dimension $K_{\rm eff}$

Moving block bootstrap over dates (block length $\ge H_{\max}$,
$B_\lambda=200$ draws), giving $\{\lambda_j^{(b)}\}$. Define

$$\boxed{K_{\rm eff}=\max\left\{k:\ \frac{\lambda_k}{\lambda_1}\ge\frac1{\kappa_{\max}}\ \ \text{and}\ \ \frac{\operatorname{sd}_b(\lambda_k^{(b)})}{\lambda_k}\le\tau_\lambda\right\},}$$

default $\kappa_{\max}=100$ (max 10× amplification), $\tau_\lambda=0.25$.

Both conditions are necessary: the first controls numerical amplification,
the second controls statistical instability. A direction with
$\lambda_k/\lambda_1=0.05$ but 60% bootstrap relative standard deviation has
a whitened coordinate that's almost entirely estimation noise, even though
it clears the condition-number threshold.

**The report must include**: the full $\{\lambda_j/\lambda_1\}$ spectrum,
$\kappa=\lambda_1/\lambda_{K_{\rm eff}}$, $K_{\rm eff}$, and the horizon
loadings $E[:,K_{\rm eff}{+}1:K]$ of the discarded directions (to see what
was actually discarded).

If $K_{\rm eff}\lt K$ and the discarded directions' loadings concentrate on
one or two horizons, the right fix is usually not to accept the
truncation — it's to go back to §6 and switch to an incremental basis or
merge horizons. Truncation is the last resort.

### §7.3 The Whitening Matrix

**Case one (no truncation, $K_{\rm eff}=K$)**: use the symmetric inverse
square root

$$B_Y=E\Lambda^{-1/2}E^\top=\Sigma_Y^{-1/2}\in\mathbb R^{K\times K}.$$

The advantage of the symmetric square root (over Cholesky) is that it
**doesn't depend on horizon ordering**, so $Z$'s $k$-th coordinate still
roughly corresponds to the $k$-th horizon when there's no truncation —
good interpretability. Cholesky whitening would turn the coordinates into a
triangular structure ("horizon 1," "horizon 2 net of horizon 1," ...) that
isn't comparable across bands.

**Case two (with truncation, $K_{\rm eff}\lt K$)**: use the reduced-dimension
form

$$B_Y=E_{1:K_{\rm eff}}\Lambda_{1:K_{\rm eff}}^{-1/2}\in\mathbb R^{K\times K_{\rm eff}}.$$

Here the coordinates' meaning shifts to "the leading $K_{\rm eff}$ principal
directions of future returns" and no longer correspond to individual
horizons. This semantic change must be stated explicitly in the report,
since each component of $c_p$ will be read against $E$'s loadings.

Both cases define

$$\boxed{Z_d=Y_d^{(\beta)}B_Y.}$$

### §7.4 Proposition 7.1 (exact whitening)

$$Z^\top WZ=I_{K_{\rm eff}}\quad\text{(exact, not approximate)}.$$

**Proof**: $Z^\top WZ=B_Y^\top\big(\frac1T\sum_d\frac1{N_d}Y_d^\top
Y_d\big)B_Y=B_Y^\top\Sigma_YB_Y$. In case two, $B_Y^\top\Sigma_YB_Y=
\Lambda^{-1/2}E^\top E\Lambda E^\top E\Lambda^{-1/2}$ (taking the first
$K_{\rm eff}$ columns) $=\Lambda^{-1/2}\Lambda\Lambda^{-1/2}=I$. Case one is
the special case. $\blacksquare$

This proposition is the foundation of §9: precisely because the whitening
is exact, $\Phi_{\mathsf D}$ can be taken directly as $Z$ with no further
orthogonalization — and therefore introduces no coordinate ambiguity.

Also, from §5's daily demeaning, $\sum_iZ_{id,:}=0$.

> **Convention**: after §7, the symbol $K$ always means $K_{\rm eff}$, the
> actually-retained label dimension.

## §8. The Three Predictive-Mechanism Blocks

### Overall Motivation

Given a whitened future-return vector $z\in\mathbb R^K$, we ask: which
**functions** of $z$ might a cross-sectional signal plausibly relate to?

In principle the answer is "any function of $z$," but that's an
infinite-dimensional space. This method picks three concrete lifting maps,
because they correspond to three prediction patterns that are completely
different in trading terms — requiring completely different handling:

| block | lifting | trading meaning | if a model only learned this |
|---|---|---|---|
| $\mathsf D$ | $z$ | who outperforms whom | long/short directly |
| $\mathsf V$ | $\operatorname{svec}(zz^\top)$ | who moves more violently, how horizons co-move | volatility/gamma positions only, no direction |
| $\mathsf L$ | $\|z\|z$ | whose direction is right during large moves | directional but heavily conditional; sizing and risk management are entirely different |

This table is the reason the whole method exists: under ordinary rank IC
these three contaminate each other; in this coordinate system they are
strictly separated.

### §8.1 The Three Liftings and Within-Date Centering

For each observation's $z_{id}\in\mathbb R^K$:

$$G_{\mathsf D}(z)=z,\qquad G_{\mathsf V}(z)=\operatorname{svec}(zz^\top),\qquad G_{\mathsf L}(z)=\|z\|_2\,z.$$

Dimensions are $K$, $K(K{+}1)/2$, $K$ respectively.

Then **all three are always within-date centered**:

$$\boxed{G_{\mathsf B,d}\leftarrow G_{\mathsf B,d}-\mathbf 1_{N_d}\,\overline{G_{\mathsf B,d}},\qquad \overline{G_{\mathsf B,d}}=\frac1{N_d}\sum_iG_{\mathsf B}(z_{id}).}$$

$G_{\mathsf D}$ already satisfies this (since $\sum_iz_{id}=0$); $G_{\mathsf
V}$ and $G_{\mathsf L}$ don't, and it must be done explicitly.

**Why this is mandatory.** By Lemma 4.1, these "that-period-constant"
components always contribute zero to $C$, so skipping centering **does not
produce a wrong $C$**. But it causes two real losses:

1. **Wasted basis dimension.** $\operatorname{orth}_W$ acting on a matrix
   with a day-common component produces some $\Phi_{\mathsf V}$ directions
   that are mostly "the market volatility time series." These directions
   are orthogonal to any demeaned $X$, forever contributing $Q=0$, yet
   occupy a $d_{\mathsf V}$ slot — diluting §19's fingerprint denominator
   and creating meaningless small singular values in §12's SVD.
2. **Distorted orthogonalization scale.** $\Xi_R=R^\top WR$ gets dominated
   by day-common variance, so $\Lambda_R^{-1/2}$ compresses genuine
   cross-sectional variance to the wrong scale, making the
   $\mathsf D\to\mathsf V\to\mathsf L$ Gram–Schmidt numerically less stable.

**Centering simplifies the definition.** Note the raw Activity lifting is
often written $\operatorname{svec}(zz^\top-I_K)$. After centering, $-I_K$ is
a global constant absorbed entirely by the that-period mean, so

$$\text{center}\big(\operatorname{svec}(zz^\top-I_K)\big)=\text{center}\big(\operatorname{svec}(zz^\top)\big).$$

This book therefore writes $G_{\mathsf V}(z)=\operatorname{svec}(zz^\top)$
directly.

**Semantic statement.** After this step, all of $\Phi$ lives in the
**within-date cross-sectional function space**. That's correct — this
method, start to finish, only cares about cross-sectional contrasts;
anything "common to all stocks that day" is outside the scope of study.

### §8.2 Definition and Verification of svec

`svec` flattens a symmetric matrix into a vector: diagonal entries unchanged,
off-diagonal entries multiplied by $\sqrt2$, upper triangle only.

**Property**: for any symmetric $A,B$,

$$\langle\operatorname{svec}A,\operatorname{svec}B\rangle=\operatorname{tr}(A^\top B).$$

**Verification**: $\operatorname{tr}(A^\top
B)=\sum_hA_{hh}B_{hh}+2\sum_{h<k}A_{hk}B_{hk}$, while
$\langle\operatorname{svec}A,\operatorname{svec}B\rangle=\sum_hA_{hh}B_{hh}+
\sum_{h<k}(\sqrt2A_{hk})(\sqrt2B_{hk})$ — the two are equal. $\blacksquare$

This property guarantees that the $\mathsf V$ block's Euclidean geometry
matches the matrix's Frobenius geometry, so $\|c_p^{\mathsf V}\|^2$ can be
read directly as "the strength of the feature's association with
second-order structure," with no coordinate artifact.

### §8.3 Content of the $\mathsf V$ Block

**Worked example ($K=2$, horizons being the 1h and 4h increments)**:
$G_{\mathsf V}$'s three coordinates are $(z_1^2,\ \sqrt2z_1z_2,\ z_2^2)$,
each minus its that-period mean after centering.

- Coordinates one and three: is this stock's 1h / 4h **magnitude** above
  the typical cross-sectional value that period. This is "activity."
- Coordinate two: the **co-directionality** of this stock's 1h and 4h
  returns. If 1h is up and 4h is also up, $z_1z_2>0$; if 1h is up and 4h
  gives it back, $z_1z_2<0$.

So the second coordinate is actually measuring **trend continuation vs.
mean reversion** — a feature that predicts it is, in substance, predicting
whether this stock's near-term return structure will be continuation-type
or reversal-type. This is valuable information, and it's completely
invisible in any single-horizon IC. This dimension is also the most
commonly misread: it is not volatility, it is **sign consistency across
horizons**.

### §8.4 Content of the $\mathsf L$ Block

$G_{\mathsf L}(z)=\|z\|z$ and $G_{\mathsf D}(z)=z$ are obviously highly
correlated (same direction, just reweighted by $\|z\|$). The $\mathsf L$
block's actual content is $G_{\mathsf L}$'s **residual** after projecting
out $\mathsf D$ and $\mathsf V$, intuitively approximately

$$R_{\mathsf L}\approx\big(\|z\|-\mathbb E\|z\|\big)\,z\ \text{, orthogonalized},$$

i.e. "**magnitude surprise × direction**."

A feature loading on the $\mathsf L$ block means: it has no directional
predictive power during calm periods, but it correctly calls direction on
the names that experience large moves. Such a feature's role in a portfolio
is a **conditional sizing signal**, not a baseline ranking signal. Its
position design in practice is completely different: it shouldn't get
equal-weighted long/short treatment — it should be used as a size modifier
on an existing directional position.

---

## §9. Hierarchical Orthogonalization

### Motivation

The three liftings' image spaces overlap. Used directly, $Q^{\mathsf D}+Q^{\mathsf
V}+Q^{\mathsf L}$ would double-count the same information, and the
fingerprint would become meaningless. Hierarchical orthogonalization
resolves this, at the cost of introducing an **order dependency**.

> **Principle 9.1 (attribution order)**: fixing the order
> $\mathsf D\to\mathsf V\to\mathsf L$ means "anything explainable by
> first-order direction is always credited to $\mathsf D$." So the correct
> reading of $Q^{\mathsf V}$ is **not** "how much second-order information
> does this feature have," but "how much second-order information does this
> feature have **that cannot be explained by first-order direction**."

This distinction is critical when reading a report. A feature that predicts
both direction and volatility will have a lower $Q^{\mathsf V}$ — not
because it doesn't predict volatility, but because its volatility
prediction is collinear with its directional prediction. The order is fixed
this way because directional information is the most directly tradeable; by
crediting the shared component to it, we ensure that whatever loads on
$\mathsf V$ and $\mathsf L$ is guaranteed to be "additional" information
that genuinely needs different handling.

### §9.1 The Weighted Orthogonalization Operator

Given a matrix $R$ (rows are observations), compute
$\Xi_R=R^\top WR$, eigendecompose $\Xi_R=U_R\Lambda_RU_R^\top$, retain
$\lambda_j>\epsilon_{\rm rank}\lambda_{\max}$ ($\epsilon_{\rm rank}=10^{-10}$,
purely numerical), and define

$$\Theta_R=U_R\Lambda_R^{-1/2},\qquad \operatorname{orth}_W(R)=R\,\Theta_R.$$

**Verification**: $(R\Theta_R)^\top W(R\Theta_R)=\Lambda_R^{-1/2}U_R^\top
\Xi_RU_R\Lambda_R^{-1/2}=I$. $\blacksquare$

If any dimension gets truncated by $\epsilon_{\rm rank}$, it must be
recorded and reported, since that means the lifting's image space actually
has lower rank than its nominal dimension (which can happen when $K$ is
large and $T$ is small).

### §9.2 The $\mathsf D$ Block: Take $Z$ Directly

$$\boxed{\Phi_{\mathsf D}:=Z.}$$

**No orthogonalization is applied.** By Proposition 7.1, $Z^\top WZ=I_K$
already holds exactly, so $\operatorname{orth}_W$ would be redundant — and
not merely redundant, actively harmful. Here's why.

If $\operatorname{orth}_W$ were applied to $Z$, then $\Xi_R=I_K$, whose
eigendecomposition is $I=U_RIU_R^\top$ — and **any orthogonal matrix $U_R$
is a valid eigenvector matrix**. So $\Theta_R=U_R$, giving

$$\Phi_{\mathsf D}=ZU_R,$$

where $U_R$ is determined by LAPACK implementation details, thread count,
even memory alignment. This is harmless for $Q^{\mathsf D}=\|c_p^{\mathsf
D}\|^2$ and any angle (orthogonal transforms preserve norms and inner
products), but it's fatal for §21's $A^{\rm tgt}_p=(c_p^{\mathsf
D})^\top r_{\mathsf D}$: if $r_{\mathsf D}$ is defined in one coordinate
system and $c_p^{\mathsf D}$ in another, the two differ by an unknown
rotation, and the alignment value would vary across execution environments.

Setting $\Phi_{\mathsf D}:=Z$ eliminates this ambiguity once and for all.
The $\mathsf D$ block's coordinates therefore have a clear economic
interpretation: the $k$-th coordinate is simply "the future-return vector's
component along the $k$-th whitened direction," roughly corresponding to
the $k$-th horizon when there's no truncation.

### §9.3 The $\mathsf V$ and $\mathsf L$ Blocks

$$\Gamma_{\mathsf{DV}}=\Phi_{\mathsf D}^\top WG_{\mathsf V},\qquad R_{\mathsf V}=G_{\mathsf V}-\Phi_{\mathsf D}\Gamma_{\mathsf{DV}},\qquad\boxed{\Phi_{\mathsf V}=\operatorname{orth}_W(R_{\mathsf V}).}$$

$$\Gamma_{\mathsf{DL}}=\Phi_{\mathsf D}^\top WG_{\mathsf L},\quad\Gamma_{\mathsf{VL}}=\Phi_{\mathsf V}^\top WG_{\mathsf L},\quad R_{\mathsf L}=G_{\mathsf L}-\Phi_{\mathsf D}\Gamma_{\mathsf{DL}}-\Phi_{\mathsf V}\Gamma_{\mathsf{VL}},$$

$$\boxed{\Phi_{\mathsf L}=\operatorname{orth}_W(R_{\mathsf L}).}$$

Finally

$$\boxed{\Phi=[\Phi_{\mathsf D}\,|\,\Phi_{\mathsf V}\,|\,\Phi_{\mathsf L}],\qquad \Phi^\top W\Phi=I_{d_\Phi}.}$$

### §9.4 Proposition 9.2 (within-block rotation invariance)

$\Phi_{\mathsf V}$ and $\Phi_{\mathsf L}$ still have internal rotation
freedom (eigenvector ordering and sign in $\operatorname{orth}_W$). This
causes no problem, because:

> **Proposition 9.2**: let $O=\operatorname{diag}(I_{d_{\mathsf
> D}},O_{\mathsf V},O_{\mathsf L})$ be any block-diagonal orthogonal matrix.
> Replacing $\Phi$ with $\Phi O$,
> $$C\mapsto CO,\quad V_r\mapsto O^\top V_r,\quad c_s\mapsto O^\top c_s,$$
> while the following are **completely invariant**: $Q_s$, $Q_s^{\mathsf
> D}$, $Q_s^{\mathsf V}$, $Q_s^{\mathsf L}$, $f_s$, $\tilde f_s$, $S^{\rm
> mech}_{pq}$, $A^{\rm mech}_{pq}$, $A^{\rm tgt}_s$, $\eta^{\mathsf D}_s$,
> $N^{\rm eff}$, $R^{\rm canc}$, $r_{\rm sig}$, $r_{\rm eff}$, and every
> $\pi$ and $M$.

**Proof**: $C\mapsto CO$ follows directly from the definition. The SVD
becomes $CO=US(O^\top V)^\top$, so singular values are unchanged,
$V\mapsto O^\top V$, hence $r_{\rm sig}$ and $r_{\rm eff}$ are unchanged;
$C_{\rm sig}=CV_rV_r^\top\mapsto CO\cdot O^\top V_rV_r^\top O=C_{\rm sig}O$,
so $c_s\mapsto O^\top c_s$, and norms and any pairwise inner product are
unchanged. Since $O$ is the identity on the $\mathsf D$ block, $c_s^{\mathsf
D}$ is **unchanged componentwise**, so $A^{\rm tgt}$ and $\eta^{\mathsf D}$
are unchanged. Each block's energy $\|c_s^{\mathsf B}\|^2$ is preserved
individually because $O$ is block-diagonal. In $g_{id}=X_{id}\langle
\Phi_{id},\hat c\rangle$, $\langle\Phi_{id}O,O^\top\hat c\rangle=\langle
\Phi_{id},\hat c\rangle$, so breadth and cancellation are unchanged.
$\blacksquare$

This proposition is the formal statement of the whole method being
"well-defined": **every reported quantity depends only on the subspace
spanned by each of $\Phi$'s three blocks, not on the choice of basis within
a block; and the $\mathsf D$ block's basis is pinned down by
$\Phi_{\mathsf D}=Z$, which in turn pins down the one quantity that has
directional meaning (target alignment).**

In implementation this gives a free unit test: rerun the whole pipeline
with different BLAS thread counts or random seeds — every quantity above
must come out identical up to numerical error (§65).

---

## §10. The Fit/Transform Boundary

### Motivation

This is the key mechanism separating "research" from "backtest." If the
label basis is refit at every evaluation period, then $c_p$ for different
periods lives in different coordinate systems, and any cross-period
comparison, any feature-vs-alpha comparison, becomes meaningless.

### Specification

During the training/research period, run `fit_label_basis(Y_train)` and
store

$$B_Y;\qquad \Theta_{\mathsf V},\ \Theta_{\mathsf L};\qquad \Gamma_{\mathsf{DV}},\ \Gamma_{\mathsf{DL}},\ \Gamma_{\mathsf{VL}}.$$

When new data $Y_{\rm new}$ arrives:

$$Z_{\rm new}=Y_{\rm new}B_Y,\qquad \Phi_{\mathsf D}^{\rm new}=Z_{\rm new},$$

$$G_{\mathsf V}^{\rm new},G_{\mathsf L}^{\rm new}\ \text{built per §8 and within-date centered},$$

$$R_{\mathsf V}^{\rm new}=G_{\mathsf V}^{\rm new}-\Phi_{\mathsf D}^{\rm new}\Gamma_{\mathsf{DV}},\qquad \Phi_{\mathsf V}^{\rm new}=R_{\mathsf V}^{\rm new}\Theta_{\mathsf V},$$

$$R_{\mathsf L}^{\rm new}=G_{\mathsf L}^{\rm new}-\Phi_{\mathsf D}^{\rm new}\Gamma_{\mathsf{DL}}-\Phi_{\mathsf V}^{\rm new}\Gamma_{\mathsf{VL}},\qquad \Phi_{\mathsf L}^{\rm new}=R_{\mathsf L}^{\rm new}\Theta_{\mathsf L}.$$

Note that within-date centering is **not a fitted parameter** — it's
applied in place on each date, so it does not constitute leakage.

### The Key Property: $\delta_{\rm drift}$

During OOS, $\Phi^{\rm new}$ no longer satisfies $\Phi^\top W\Phi=I$. This
isn't an error — it's a free regime-drift detector:

$$\boxed{\delta_{\rm drift}=\big\|\Phi_{\rm OOS}^\top W_{\rm OOS}\Phi_{\rm OOS}-I\big\|_F,}$$

further split into the three diagonal blocks' $\delta^{\mathsf
D},\delta^{\mathsf V},\delta^{\mathsf L}$.

If $\delta_{\rm drift}$ is small, future returns' covariance structure
matches the research period, and every quantity expressed via $c$ can be
compared directly across periods. If $\delta_{\rm drift}$ is large (say
$>0.3$), market structure has shifted, and a drop in $Q_\alpha$ **could come
from coordinate-system invalidity rather than model failure** — an
alternative explanation that must be ruled out first whenever reading any
OOS decay.

### Worked Example

Research period 2022–2023, OOS is the first half of 2024. If 2024H1's
market volatility is half the research period's, $\|z\|$ is systematically
smaller; since $G_{\mathsf V}\sim\|z\|^2$,

$$\Phi_{\mathsf V}^{\rm OOS\top}W\Phi_{\mathsf V}^{\rm OOS}\approx0.25\,I\ \Longrightarrow\ \delta^{\mathsf V}\approx0.75\sqrt{d_{\mathsf V}}\approx2.4\ (d_{\mathsf V}=10).$$

$Q_\alpha^{\mathsf V}$ will mechanically drop by roughly 4×, which
**cannot be read as decay in the model's activity-predicting power.** The
correct handling is to report $\delta^{\mathsf V}$ alongside it and, when
necessary, make cross-period comparisons on a $\delta$-corrected basis.

---

# Part Three: Feature–Label Geometry

## §11. The Predictive Operator

### Definition

For each period,

$$\boxed{C^{(d)}=\frac1{N_d}X_d^{XS\top}\Phi_d\ \in\mathbb R^{P\times d_\Phi},}$$

over the full period (each feature averaged over its own valid-day set
$\mathcal D_p$, per §4.2),

$$\boxed{C_{p,:}=\frac1{T_p}\sum_{d\in\mathcal D_p}C^{(d)}_{p,:}.}$$

When every feature has 100% coverage, this reduces to
$C=\frac1T\sum_dC^{(d)}=X^\top W\Phi$.

$C$ is the **average daily cross-sectional predictive operator**. Row $p$,
$C_{p,:}$, is feature $p$'s average cross-sectional association with every
future-return mechanism coordinate.

### Scale: Every Element of $C$ Is Literally an IC

By §4, $X^{XS}_{\cdot p}$ has zero mean and unit variance each period; by
§9, each column of $\Phi$ has unit norm under the weighted inner product.
So

$$C_{p,k}\approx\overline{\operatorname{corr}_d\big(X_{\cdot dp},\Phi_{\cdot dk}\big)}$$

— feature $p$'s **average per-period cross-sectional correlation** with the
$k$-th mechanism coordinate, i.e. literally an IC. Hence

$$Q_p=\|c_p\|^2\approx\sum_k\overline{IC}_{p,k}^2.$$

This gives the scale for reading $Q$:

| Average per-period IC | $Q$ contribution of one coordinate | Meaning in equity markets |
|---:|---:|---|
| 0.005 | $2.5\times10^{-5}$ | almost undetectable |
| 0.01 | $1\times10^{-4}$ | weak, but tradeable at high frequency |
| 0.03 | $9\times10^{-4}$ | a good single feature |
| 0.05 | $2.5\times10^{-3}$ | very strong |
| 0.10 | $1\times10^{-2}$ | extremely rare; suspect data leakage first |
| 0.20 | $4\times10^{-2}$ | almost certainly look-ahead bias |

**Every numeric example in this book falls in the band $Q^{\mathsf
D}\in[10^{-5},10^{-3}]$.** Getting this scale wrong by one or two orders of
magnitude will mis-set every threshold downstream, deleting all real
features or keeping all noise.

### Scale of Computation

With $P=400$, $d_\Phi=18$, $C$ is a $400\times18$ matrix — about 7200
numbers. The entire raw material of a feature study is these 7200 numbers
plus their null distribution.

---

## §12. Singular Value Decomposition

### Definition

$$\boxed{C=USV^\top,}\qquad S=\operatorname{diag}(s_1,\dots),\ s_1\ge s_2\ge\cdots$$

$V$'s columns are label-mechanism directions; $U$'s columns are feature-side
combinations.

### Interpretation

$C$'s SVD answers: **across the $d_\Phi$-dimensional label-mechanism space,
which directions does the entire feature library collectively
"illuminate"?**

- $v_1$: the whole library's most-agreed-upon predictive direction. On most
  equity data this will be a direction dominated by the $\mathsf D$ block,
  close to "all horizons point the same way" — the classic
  momentum/reversal axis.
- $u_1$: the feature-side linear combination achieving this prediction — can
  be thought of as an automatically-extracted super-feature.
- $s_1$: that direction's strength.

**Note the asymmetry between $U$ and $V$.** $V$ lives in label space
(dimension $d_\Phi$, small); $U$ lives in feature space (dimension $P$,
large). Therefore

$$r\le\min(P,d_\Phi)=d_\Phi.$$

**The ceiling on the number of detectable modes is set by the label space's
dimension, not the number of features.** At $K=4$, $r\le18$; in practice
$r_{\rm sig}$ almost always lands between 2 and 6. This fact has decisive
consequences in §38.

---

## §13. The Permutation Null

### Motivation

$Q_p$ needs a reference distribution. It cannot be given in closed form,
because $C$'s null distribution depends simultaneously on the feature's
cross-sectional distribution, $\Phi$'s cross-sectional distribution, and
both of their own temporal autocorrelation structure. The only reliable
approach is a resampling scheme that preserves all of that structure and
destroys only "which day pairs with which."

### §13.1 The Null Hypothesis

The original pairing is $X_d\leftrightarrow\Phi_d$. The null hypothesis is:

> the feature's cross-sectional ranking has no predictive relationship with
> the future-return mechanism **at the correct time**.

So we break date alignment: $X_{\pi(d)}\leftrightarrow\Phi_d$, while
preserving

1. each period's full cross-sectional structure in $X$ (including sector,
   market-cap clustering, etc.);
2. each period's full cross-sectional structure in $\Phi$;
3. ticker identity;
4. $X$'s own temporal autocorrelation;
5. the label's own temporal autocorrelation.

Preserving these five is what keeps the null from being artificially small
just because cross-sectional structure got scrambled — that's the most
common error, and it makes everything look significant.

### §13.2 Primary Scheme: Circular Date Shift

Given ordered dates, pick a shift $s$, and set

$$\pi_s(k)=1+\big((k+s-1)\bmod T\big).$$

The set of valid shifts is

$$\mathcal S=\{s:\Delta\le\min(s,T-s)\},\qquad\boxed{\Delta=H_{\max}+\ell_X,}$$

where $H_{\max}$ is the longest horizon in the band (in periods), and
$\ell_X$ is the feature's memory length (e.g. $\ell_X=20$ for a 20-day
moving average). This minimum offset ensures the paired $X$ and $Y$ never
share an overlapping time window. **$\Delta$ and the purge gap in §41 must
be the same number** — otherwise the same data is held to different leakage
standards in two places.

In practice, use $B_{\rm null}$ of the $|\mathcal S|$ possible shifts
(default: take all, capped at 2000).

### §13.3 Universe Mismatch: the Matched-Thinning Rule

$N_d$ changes over time, so $\mathcal U_{\pi_s(d)}\ne\mathcal U_d$. Three
things must be handled together.

**(a) Take the intersection.** Let

$$\mathcal I_{d,s}=\mathcal U_{\pi_s(d)}\cap\mathcal U_d.$$

If $|\mathcal I_{d,s}|<\rho_{\min}\cdot\min(N_{\pi_s(d)},N_d)$ (default
$\rho_{\min}=0.9$), discard the $(d,s)$ pair.

**(b) Re-standardize on the intersection.** By Lemma 4.1, $C^{(d)}$ is a
genuine covariance *provided* $X$ has zero mean on the universe used.
Restricting to $\mathcal I_{d,s}$ breaks this precondition: if a feature
correlates with "being present on both dates" (market cap, liquidity,
listing age naturally do), the intersection's mean is nonzero, and the
null's $C$ gets systematically contaminated by a nonzero mean term,
inflating or deflating the null for that class of feature. So **$X$ must be
re-demeaned and re-standardized on $\mathcal I_{d,s}$**:

$$\check X_{idp}=\frac{X'_{idp}-\mu_{dp}^{\mathcal I}}{\sigma_{dp}^{\mathcal I}+\epsilon},\qquad i\in\mathcal I_{d,s}.$$

$\Phi$ needs no reprocessing: by Lemma 4.1, as long as $\check X$ has zero
mean on $\mathcal I$, any that-period-constant component of $\Phi$ is
automatically absorbed.

**(c) Apply the same sampling rule to the observed statistic.** If the
observed statistic is computed on the full $\mathcal U_d$ while the null is
computed on a smaller intersection, the two have different effective sample
sizes, and the null's variance will be over- or under-estimated depending
on whether the intersection is selective. The symmetric rule is:

> **Specification 13.1 (matched thinning)**: for each date $d$, draw $s_d
> \sim\text{Uniform}(\mathcal S)$ independently, define the observation
> analysis universe
> $$\mathcal A_d=\mathcal U_d\cap\mathcal U_{\pi_{s_d}(d)},$$
> and re-standardize on $\mathcal A_d$ per (b) before computing the observed
> statistic. Repeat $R_{\rm thin}=5$ independent draws and average, as the
> final observed value.

This way, the observed and null statistics have identically distributed
universe sizes.

**Simplifying condition**: if $\operatorname{median}_{d,s}|\mathcal
I_{d,s}|/N_d\ge0.98$, thinning's effect is smaller than numerical error,
and observed values can be computed directly on $\mathcal U_d$. This ratio
must be computed and recorded, never assumed.

### §13.4 Slow-Moving Features and an Alternative Null

The circular-shift test's power depends on whether, after shifting by
$\Delta$ periods, a feature's cross-sectional ranking has actually changed.
For a feature with memory much longer than $\Delta$ (e.g. 250-day momentum,
market cap, sector dummies),

$$X_{d+\Delta,p}\approx X_{d,p},$$

so $C^{\rm null}\approx C^{\rm obs}$, and the test has almost no power.
**This systematically biases the method toward favoring fast-rotating
signals** — if left unhandled, slow but genuine signals get judged
not-significant across the board.

**Persistence diagnostic.** For each feature, define

$$\boxed{\kappa_p=\operatorname{median}_d\ \operatorname{corr}_{i\in\mathcal U_d\cap\mathcal U_{d+\Delta}}\big(X_{idp},X_{i,d+\Delta,p}\big).}$$

$\kappa_p$ must appear in the report. Rule of thumb: the circular-shift
null's power has dropped noticeably once $\kappa_p>0.5$, and is near zero
by $\kappa_p>0.8$.

**Alternative null (N2): stratified within-date permutation.** For features
with $\kappa_p>0.5$, replace the date shift with

$$X_{idp}\longrightarrow X_{\pi_d(i),dp},$$

where $\pi_d$ is a random permutation **within strata**:

$$\mathcal U_d=\bigsqcup_g\mathcal U_{d,g},\qquad g=\text{sector}\times\text{size bucket},$$

permuted independently within each $g$.

**Stratification is necessary, not a refinement.** Permuting freely across
the whole universe would destroy both sector and market-cap clustering at
once, making the null for a sector-level feature far too small, and the
test severely anti-conservative. Stratification fixes these known
exposures' cross-sectional composition in place, so the null only
destroys "individual-stock-level pairing."

A cleaner variant residualizes known exposures first: let $F$ be the known
factor-exposure matrix (sector, size, beta),

$$X^\perp_{d,p}=(I-P_{F_d})X_{d,p},\qquad P_{F_d}=F_d(F_d^\top F_d)^{-1}F_d^\top,$$

then run free within-date permutation on $X^\perp$. The null hypothesis
becomes "after controlling for known exposures, individual-stock-level
feature ranking has no predictive power" — more precise semantics, at the
cost that the conclusion covers only the residual portion.

**Reporting rule**: for features with $\kappa_p>0.5$, the primary conclusion
uses N2, with N1's result and $\kappa_p$ also listed; for $\kappa_p\le0.5$,
the primary conclusion uses N1. Which null was used must be marked per
feature in the report.

### §13.5 Implementation Flow

```text
Input: {X_d}, {Phi_d}, valid shift set S, number of resamples B_null
Precompute: kappa_p (persistence diagnostic) for each feature

# ---- observed statistic (matched thinning) ----
for r in 1..R_thin:
    for d in 1..T:
        s_d ~ Uniform(S)
        A_d = U_d ∩ U_{pi_{s_d}(d)}
        if |A_d| < rho_min * min(N_d, N_{pi_{s_d}(d)}): skip this d
        X_check = X_d re-standardized on A_d
        C_r += (1/|A_d|) * X_check^T Phi_d[A_d]
    C_obs_r = C_r / (number of valid days)
C_obs = mean_r C_obs_r

# ---- N1: circular-shift null ----
for b in 1..B_null:
    s = S[b]
    for d in 1..T:
        I = U_{pi_s(d)} ∩ U_d
        if |I| < rho_min * min(...): skip
        X_check = X_{pi_s(d)} re-standardized on I
        C_b += (1/|I|) * X_check^T Phi_d[I]
    C_null[b] = C_b / (number of valid days)

# ---- N2: stratified within-date permutation (slow-moving features only) ----
for b in 1..B_null:
    for d in 1..T:
        for g in strata(d):
            randomly permute X's rows within U_{d,g}
        C_b += (1/N_d) * X_perm_d^T Phi_d
    C_null2[b] = C_b / T
```

Every downstream null-based quantity ($r_{\rm sig}$, feature p-values,
similarity thresholds, fingerprint baselines, breadth and cancellation
baselines) is derived from this same $\{C^{\rm null}_b\}$, so this loop
only needs to run once.

---

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

$$R_j=P_U^\perp CP_V^\perp,\qquad P_U^\perp=I-U_{\lt j}U_{\lt j}^\top,\quad P_V^\perp=I-V_{\lt j}V_{\lt j}^\top,$$

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

---

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

---

## §23. Mechanism Similarity

### Definition

$$\hat c_p=\frac{c_p}{\|c_p\|+\epsilon},\qquad\boxed{A^{\rm mech}_{pq}=\hat c_p^\top\hat c_q\in[-1,1],\qquad S^{\rm mech}_{pq}=|A^{\rm mech}_{pq}|\in[0,1].}$$

$A^{\rm mech}$ keeps the sign, used to tell whether two features predict the
same mechanism in the same or opposite direction. $S^{\rm mech}$ drops the
sign, used for clustering — since a feature and its negated version carry
identical information and must be treated as the same thing when
clustering.

### §23.1 The Energy Interpretation of Similarity

This quantity has a precise decision-theoretic meaning. If feature $p$ is
used as feature $q$'s representative (i.e. $p$ is kept, $q$ is dropped in
selection), $q$'s mechanism vector is approximated by $p$'s direction, and
the fraction of energy retained is

$$\frac{\|\text{proj}_{\hat c_p}c_q\|^2}{\|c_q\|^2}=\big(S^{\rm mech}_{pq}\big)^2,$$

with fraction lost $1-(S^{\rm mech}_{pq})^2$.

| $S^{\rm mech}$ | Energy retained $S^2$ | Loss |
|---:|---:|---:|
| 0.95 | 0.90 | 10% |
| 0.90 | 0.81 | 19% |
| 0.80 | 0.64 | 36% |
| 0.70 | 0.49 | 51% |

**This is the correct basis for the selection threshold $\tau$**: $\tau$
directly expresses "how much mechanism energy am I willing to lose by
substituting the representative for the member." $\tau=0.90$ means
accepting up to 19% loss. This criterion is independent of dimension and of
library size.

### §23.2 Null Calibration: How Much Similarity Comes From Chance

The decision criterion above answers "should these be merged," but not
"is the observed similarity structure even real." The latter needs null
calibration.

**The null used for calibration differs from §13's.** §13 uses a
**shared** shift (all features shifted together), to preserve inter-feature
dependence — needed for the $Q_p$ test. But similarity calibration asks "how
similar can two features with **no shared mechanism** be," so it must use
**independent per-feature shifts**:

$$X_{d,p}\longrightarrow X_{\pi_{s_p}(d),p},\qquad s_p\ \text{drawn independently for each feature}.$$

This makes $c_p^{\rm null}$ and $c_q^{\rm null}$ independent, while each
still retains its true marginal geometry (coverage, persistence,
cross-sectional distribution).

This gives the null similarity distribution $\{S^{\rm null}_{pq,b}\}$, from
which we define the **excess pairing rate**

$$\boxed{\zeta(\tau)=\underbrace{\frac{\#\{(p,q):S^{\rm mech}_{pq}>\tau\}}{P(P-1)/2}}_{\text{observed}}-\underbrace{\mathbb E_b\left[\frac{\#\{(p,q):S^{\rm null}_{pq,b}>\tau\}}{P(P-1)/2}\right]}_{\text{null}}.}$$

### §23.3 A Fundamental Limit of Low-Dimensional Spaces

The following calculation must be understood, or clustering results will be
badly misread.

If $c_p$ and $c_q$ are two independent, uniformly distributed unit vectors
in an $m$-dimensional space, the tail probability of $|\cos|$ is

| $m$ | $P(|\cos|>0.9)$ | $P(|\cos|>0.95)$ |
|---:|---:|---:|
| 2 | 28.7% | 20.3% |
| 3 | 10.0% | 5.0% |
| 4 | 3.74% | 1.33% |
| 5 | 1.45% | 0.36% |
| 6 | 0.58% | 0.10% |

($m=3$ is the easiest to remember: by Archimedes' hat-box theorem, a
uniformly distributed vector on the 3D unit sphere has a cosine with any
fixed direction that is **uniformly distributed** on $[-1,1]$, so
$P(|\cos|>\tau)=1-\tau$.)

**What this means.** The $\mathsf D$ block's similarity $S^{\mathsf D}$
(§24) lives in a space of dimension

$$m_{\mathsf D}=\dim\big(\text{span}(V_r)\cap\mathsf D\text{ coordinates}\big)\le\min(K,r_{\rm sig}),$$

typically 2 or 3. At $m_{\mathsf D}=3$, $P=400$ (about $8\times10^4$ pairs),
the null already expects **roughly 8000 pairs with $S>0.9$**. Any FWER- or
Bonferroni-calibrated threshold would converge to 1 and become useless.

This isn't a calibration failure — it's a genuine structural fact: **in a
three-dimensional space, four hundred directions cannot help being "close"
to each other.**

### §23.4 Design Conclusions That Follow

> **Rule 23.1**: mechanism clustering (§35) is a **descriptive and
> organizational** tool, for grouping in reports and for human
> interpretation; **it is not the feature-selection decision mechanism.**
> Selection is decided by the coverage criterion of §38–§39, which directly
> measures "does the selected set represent every eligible feature" — and
> therefore doesn't depend on the significance of any pairwise threshold.
>
> **Rule 23.2**: the clustering threshold $\tau$ is set by §23.1's
> energy-loss criterion (default $\tau=0.90$, corresponding to 19% loss),
> and **must** be reported alongside $\zeta(\tau)$ and the effective
> mechanism dimension
> $$m_{\mathsf B}=\frac{\big(\sum_j\lambda_j^{\mathsf B}\big)^2}{\sum_j(\lambda_j^{\mathsf B})^2},\qquad \lambda^{\mathsf B}_j=\big\|v_j[\mathsf B]\big\|^2,$$
> i.e. $V_r$'s participation ratio of energy on block $\mathsf B$.
>
> If $\zeta(\tau)\le0$ or $m_{\mathsf D}<4$, the report must explicitly flag
> "mechanism clustering has no discriminative power under this setting," and
> **no cluster-level statement may be used as a conclusion.**

This is an honest and important conclusion: in a low-dimensional mechanism
space, "are these two features the same thing" cannot be answered by angle.
What can be answered is "which features need to be kept to represent the
whole" — and that's §38's question.

---

## §24. Block-Specific Mechanism Similarity

### Definition

$$\boxed{A^{\mathsf B}_{pq}=\frac{c_p^{\mathsf B\top}c_q^{\mathsf B}}{\|c_p^{\mathsf B}\|\|c_q^{\mathsf B}\|+\epsilon},\qquad S^{\mathsf B}_{pq}=\big|A^{\mathsf B}_{pq}\big|,\qquad \mathsf B\in\{\mathsf D,\mathsf V,\mathsf L\}.}$$

### Why This Section Can't Be Skipped

Full-space $S^{\rm mech}_{pq}$ can be high without the two features sharing
the same directional mechanism at all. Consider two features whose $c$ is
both dominated by the $\mathsf V$ block ($f_{\mathsf V}\approx0.8$) with
nearly identical $\mathsf V$ components, but orthogonal $\mathsf D$
components. Then

$$S^{\rm mech}_{pq}\approx0.8,\qquad S^{\mathsf D}_{pq}\approx0.$$

For **directional research**, these two features are not interchangeable
at all.

> **Rule 24.1**: every similarity computation involved in directional
> selection (§35's clustering, §38's coverage, §52's inheritance reading)
> uses $S^{\mathsf D}$, never $S^{\rm mech}$. The full-space $S^{\rm mech}$
> is only for the descriptive question "do these two features look alike
> overall."

The corollary: if the research target is volatility signals, switch
entirely to $S^{\mathsf V}$. The choice of block must match the research
target — the same applies to §27's stability and §34's thresholds.

---

## §25. Raw Feature Redundancy

### Motivation

$S^{\mathsf D}$ measures "is the predicted mechanism the same." But two
features might predict the same mechanism while carrying different raw
information (e.g. one from price, one from volume, that happen to both
predict short-term reversal). Keeping both may have value — their noise
differs, and combining them improves the signal-to-noise ratio.

Conversely, two raw-nearly-identical features (e.g. 20-day and 21-day
momentum) are just duplicates; keeping both is pure redundancy.

So a second layer of redundancy, measured directly on the raw data, is
needed.

### Definition

Compute daily Spearman correlation

$$\rho_{pq,d}=\operatorname{Spearman}\big(X_{d,p},X_{d,q}\big),$$

$$\boxed{S^{\rm raw}_{pq}=\operatorname{median}_d\big|\rho_{pq,d}\big|.}$$

Spearman rather than Pearson because raw features often have fat tails;
median rather than mean because single-date extremes (e.g. an earnings-day
data anomaly) shouldn't dominate the redundancy judgment.

### Relationship to $S^{\mathsf D}$

The two measure different things, and **there's no necessary ordering
between them**:

| | $S^{\rm raw}$ high | $S^{\rm raw}$ low |
|---|---|---|
| $S^{\mathsf D}$ high | genuine duplicate, keep only one | different source, same mechanism: possibly keep both (different noise) |
| $S^{\mathsf D}$ low | raw-correlated but predicts different mechanisms: keep both | completely unrelated: two different things |

The upper-right case is the most valuable in building a feature library:
two independent data sources pointing at the same predictive mechanism is
strong evidence the mechanism is real, and combining them can reduce
measurement noise.

---

## §26. Temporal Stability: Bootstrap Design

### Motivation

$Q_p$ is a full-period average. A feature that's very strong in the first
half of the sample and vanishes in the second half can have the same $Q_p$
as a stably-moderate feature, but the former is unusable. We need to
measure whether a conclusion reproduces across different time sub-samples.

### Specification

**Moving block bootstrap (over dates).** Set block length
$\ell_{\rm blk}\ge H_{\max}$; each resample draws
$\lceil T/\ell_{\rm blk}\rceil$ starting points and concatenates them into a
length-$T$ date sequence (with replacement), giving

$$\boxed{C^{(b)}=\frac1{T}\sum_{d\in\text{bootstrap}_b}C^{(d)},\qquad b=1,\dots,B_{\rm time}.}$$

Default $B_{\rm time}=500$–$1000$.

**Must use moving block, not iid bootstrap**: $C^{(d)}$ is autocorrelated in
time (features persist, return volatility clusters). An iid resample would
destroy this structure and understate variance, making every $\pi$ look
over-optimistic. The block length must be at least $H_{\max}$ to ensure
overlapping label windows fall within the same block.

**$\Phi$ and $V_r$ are held fixed across all resamples.** This is
deliberate: this bootstrap measures "**given the coordinate system**, how
stable is the feature-level association." The coordinate system's own
uncertainty is measured by two other mechanisms: §7.2's spectral bootstrap
(stability of $\Sigma_Y$) and §10's $\delta_{\rm drift}$ (cross-period
drift).

> **Rule 26.1**: $\pi$ and $M$ are **lower bounds** on total uncertainty. The
> report must present $\delta_{\rm drift}$ and $K_{\rm eff}$'s bootstrap
> stability alongside them, or a reader will mistake "stable given the
> coordinate system" for "stable overall."

---

## §27. Temporal Stability Metrics

### Definitions

All metrics are **block-specific**: $\mathsf B$ must match the research
target (directional research always uses $\mathsf D$).

**Survival rate (strength)**:

$$\boxed{\pi^{{\rm time},\mathsf B}_p=\frac1{B_{\rm time}}\sum_b\mathbf 1\left[Q^{\mathsf B,(b)}_p>Q^{\mathsf B,{\rm null}}_{p,0.95}\right],}$$

where $Q^{\mathsf B,\rm null}_{p,0.95}$ is that feature-block's null
distribution's 95th percentile.

**Sign stability rate**:

$$\boxed{\pi^{\rm time,sign}_p=\frac1{B_{\rm time}}\sum_b\mathbf 1\left[\operatorname{sign}\big(A^{{\rm tgt},(b)}_p\big)=\operatorname{sign}\big(A^{\rm tgt}_p\big)\right].}$$

**Mechanism stability**:

$$M^{{\rm time},\mathsf B}_{p,b}=\frac{c_p^{\mathsf B,(b)\top}c_p^{\mathsf B}}{\|c_p^{\mathsf B,(b)}\|\,\|c_p^{\mathsf B}\|+\epsilon},\qquad\boxed{M^{{\rm time},\mathsf B}_p=\operatorname{median}_bM^{{\rm time},\mathsf B}_{p,b}.}$$

**Coefficient of variation of strength**:

$$\boxed{CV^{Q,\mathsf B}_p=\frac{\operatorname{sd}_b\big(Q^{\mathsf B,(b)}_p\big)}{\operatorname{mean}_b\big(Q^{\mathsf B,(b)}_p\big)+\epsilon}.}$$

### What Each of the Four Quantities Answers

| Quantity | Question | What failure means |
|---|---|---|
| $\pi^{\rm time,\mathsf B}$ | does strength persistently clear noise across sub-samples | effective only in a minority of periods |
| $\pi^{\rm time,sign}$ | is the direction stable | direction can flip — the most dangerous failure |
| $M^{\rm time,\mathsf B}$ | is the **mechanism** (direction of $c$) stable | what's being predicted is itself changing |
| $CV^{Q,\mathsf B}$ | how much does strength fluctuate | high $CV$ with high $\pi$: strength swings but never disappears |

The distinction between $M$ and $\pi$ matters especially. $\pi^{\rm time}$
only checks "is strength big enough"; $M^{\rm time}$ checks "has the
direction changed." A feature can have $\pi=0.95$ (always strong) but
$M=0.2$ (predicting a different mechanism in every sub-sample) — meaning
it's an unstable composite signal whose large $Q$ is only because it
associates with *something different* each period. Such a feature is
entirely unusable, and looking only at $\pi$ would miss it.

### Relationship to the P-Value

$\pi^{\rm time,\mathsf B}$ and $\tilde p_p$ look redundant but answer
different questions: $\tilde p_p$ asks "could the full-sample point
estimate plausibly be chance," $\pi$ asks "how often does this association
show up across time sub-samples." A feature that's genuine but only
effective 30% of the time can have a tiny $\tilde p$ and a $\pi$ of just
0.3. Both must be part of eligibility.

### Note: the $\pi^{\rm time,\mathsf B}$ Threshold Comparison

In $\pi$'s definition, the bootstrap sample's $Q$ is compared against the
**full-sample** null's 95th percentile. Strictly, a bootstrap sub-sample has
fewer effective periods, so its $Q$'s variance should be somewhat larger
than the full sample's — making this comparison slightly conservative.
This is acceptable: $\pi\ge0.8$ is already a screening threshold, not an
exact test. A more precise version would resample the null alongside every
bootstrap draw (multiplying cost by $B_{\rm time}$), but in practice the
gain doesn't justify the cost.

---
## §28. Two Failure Modes of Instability

The combination of $\pi$ and $M$ splits instability into two categories,
requiring completely different treatment.

### Mode One: Strength Decay (low $\pi^{\rm time,\mathsf D}$, high $M^{\rm time,\mathsf D}$)

The feature's predicted **mechanism hasn't changed** — its strength just
drops below noise in some periods.

Typical causes: a market regime change (volatility, liquidity) altered the
signal's detectability; or the signal is being gradually arbitraged away
(alpha decay).

Treatment: such a feature can still be used, but needs more information to
decide how. The key diagnostic is **whether the decay is monotonic**:
regress $Q^{\mathsf D,(b)}_p$ against the bootstrap samples' average date;
if the slope is significantly negative, this is alpha decay and should be
downweighted or retired; if strength correlates with an observable regime
variable (VIX, volume, spread), this is a conditional signal that should be
paired with a regime filter.

### Mode Two: Mechanism Drift (low $M^{\rm time,\mathsf D}$)

The feature predicts **different things** in different periods. $c_p$'s
direction is itself changing.

Typical causes: the feature's definition changed with the data supply
(e.g. a data field's coverage expanded partway through); or the feature is
a mixture of multiple underlying effects whose relative strength shifts
over time.

Treatment: such a feature **should not be used**, even if $Q_p$ is large and
$\tilde p_p$ is tiny. The reason: what it will predict during OOS is
completely unpredictable. The correct fix is to go back to the feature
construction layer and decompose it into more basic components, testing
each separately.

$M^{\rm time,\mathsf D}<0.3$ should be treated as a hard exclusion
condition, not a scoring input.

### Sign Flips ($\pi^{\rm time,sign}$ low)

This is a special case of mode two, but flagged separately because it's the
most dangerous. A feature whose sign can flip turns profit directly into
loss in a position. Features with $\pi^{\rm time,sign}<0.8$ are always
retired, and **cannot** be rescued by "flipping the sign back in OOS" (§21).

---

## §29. Cross-Sectional Stability

### Motivation

Temporal stability asks "does this hold across different periods."
Cross-sectional stability asks "is it driven by only a handful of stocks."
The two are independent failure dimensions: a feature can be effective
every single day ($\pi^{\rm time}=1$) while every day's effect comes from
the same 20 small-cap names.

### Specification

Keep all dates. For each date $d$, re-bootstrap tickers within
$\mathcal U_d$. If sector/cap metadata is available, use a **stratified**
bootstrap:

$$\mathcal U_d=\bigsqcup_g\mathcal U_{d,g},$$

drawing with replacement $N_{d,g}$ times within each stratum, ending with
$N_d$ rows total.

**Stratification is necessary**: without it, some bootstrap samples would
randomly over-represent a sector, so $C^{(d,b)}$'s variance would reflect
both "sector-composition sampling noise" and "the feature's cross-sectional
concentration" at once. Stratifying fixes sector composition in place, so
the result reflects only the latter.

### Every Resample Must Be Re-Centered

After sampling with replacement, $X$'s mean on the resampled set is
**no longer zero** (some names get drawn multiple times). By Lemma 4.1,
$C^{(d)}$ being a clean covariance depended exactly on $X^\top\mathbf 1=0$.
Once that precondition breaks,

$$\frac1{N_d}X_{d,b}^\top\Phi_{d,b}=\underbrace{\operatorname{Cov}_b(X,\Phi)}_{\text{wanted}}+\underbrace{\bar X_{d,b}\,\bar\Phi_{d,b}^\top}_{\text{spurious signal}}.$$

For a feature unrelated to sampling probability, the spurious term is just
extra noise; but for a feature correlated with stratum size (market cap,
liquidity, sector dummies), it produces a systematic bias.

> **Specification 29.1**: after every bootstrap resample, $X$ must be
> re-demeaned and re-standardized on the resampled sample:
> $$\check X_{idp}=\frac{X_{idp}-\bar X^{(b)}_{dp}}{\sigma^{(b)}_{dp}+\epsilon}.$$

The same requirement applies to §13.3(b)'s permutation intersection. Both
are applications of the same principle.

### Output

$$C^{(b)}_{\rm XS}=\frac1T\sum_d\frac1{N_d}\check X_{d,b}^\top\Phi_{d,b},$$

holding $V_r$ fixed, giving — exactly as in §27 —

$$\pi^{{\rm XS},\mathsf B}_p,\qquad \pi^{\rm XS,sign}_p,\qquad M^{{\rm XS},\mathsf B}_p.$$

Default $B_{\rm XS}=500$.

---

## §30. Predictive Breadth

### Motivation

$Q_p$ is an average, and doesn't tell us how many stocks contributed to it.
An IC contributed to evenly by 3000 stocks and an IC contributed to by 20
stocks (with the other 2980 contributing nothing) are completely different
objects in trading terms: the latter has tiny capacity, is extremely
sensitive to a single event, and is very likely to be eaten by trading
costs in live trading.

### Definition

For feature $p$, let $\hat c_p=c_p/\|c_p\|$. Stock $i$'s **signed
contribution** to the predictive relationship on date $d$ is

$$\boxed{g_{id,p}=X^{XS}_{id,p}\big\langle\Phi_{id,:},\hat c_p\big\rangle.}$$

This has an exact meaning: $\frac1{N_d}\sum_ig_{id,p}$ is precisely $c_p$'s
that-period contribution along the $\hat c_p$ direction, so

$$\frac{1}{T}\sum_d\frac1{N_d}\sum_ig_{id,p}=\hat c_p^\top\,c_p^{\rm full}\approx\|c_p\|.$$

Absolute contribution weight and effective number of contributors:

$$u_{id,p}=\frac{|g_{id,p}|}{\sum_j|g_{jd,p}|},\qquad\boxed{N^{\rm eff}_{d,p}=\frac1{\sum_iu_{id,p}^2},\qquad 1\le N^{\rm eff}_{d,p}\le N_d,}$$

$$\boxed{b_{d,p}=\frac{N^{\rm eff}_{d,p}}{N_d}\in(0,1].}$$

The report uses

$$b^{\rm med}_p=\operatorname{median}_db_{d,p},\qquad b^{\rm p10}_p=\text{the 10th percentile}_d(b_{d,p}).$$

$N^{\rm eff}$ is the inverse-Herfindahl of the $\{|g|\}$ distribution — the
"effective number of contributing stocks."

---

## §31. Reading Breadth and Its Null Baseline

### §31.1 Why $b=1$ Is Not the Reference Point

Intuitively, $b=1$ (every stock contributes equally) feels like the
reference point for "fully diversified." That intuition is wrong.

Consider a feature with **zero** predictive power. Then
$g_{id}=X_{id}\cdot\langle\Phi_{id},\hat c\rangle$ is the product of two
approximately independent, approximately standard-normal quantities. As
$N\to\infty$,

$$b=\frac{\big(\frac1N\sum|g|\big)^2}{\frac1N\sum g^2}\longrightarrow\frac{(\mathbb E|g|)^2}{\mathbb E[g^2]}=\frac{(\mathbb E|u|\cdot\mathbb E|v|)^2}{\mathbb E[u^2]\mathbb E[v^2]}=\frac{(2/\pi)^2}{1}\approx0.405.$$

(Using $\mathbb E|u|=\sqrt{2/\pi}$ for a standard normal.)

So **pure noise has breadth around 0.405, not 1**. (If the contribution
were a single normal quantity rather than a product of two, the baseline
would be $2/\pi\approx0.637$.) Seeing $b=0.35$ and calling it "highly
concentrated" is wrong — that's essentially the noise level.

### §31.2 Relative Breadth

Real $X$ and $\Phi$ are not independent normals ($X$ usually has fat tails,
$\Phi_{\mathsf V}$'s components are quadratic forms and thus heavily
skewed), so 0.405 is only an order-of-magnitude anchor, not a hard
threshold. The report uses a **relative** quantity instead:

$$\boxed{\tilde b_p=\frac{b^{\rm med}_p}{\operatorname{median}_b\big(b^{{\rm med},{\rm null}}_{p,b}\big)},}$$

with the null coming from §13's same resamples. $\tilde b_p\approx1$ means
concentration is no different from noise; $\tilde b_p\ll1$ means genuine
concentration.

**Why only the relative quantity is meaningful**: $b$'s null baseline
differs by block. The $\mathsf V$ block's $\Phi$ components are quadratic
forms of $z$, with much higher kurtosis than normal, so its null breadth is
notably below 0.405. Using the same absolute threshold for $\mathsf D$ and
$\mathsf V$ block features would systematically flag every $\mathsf V$
block feature as "over-concentrated."

### §31.3 Concentration Is Not Necessarily a Defect

> **Rule 31.1**: breadth is a soft diagnostic, not a hard gate.

Concentration itself isn't necessarily a flaw — some genuine mechanisms
(event-driven, earnings surprise, index reconstitution) naturally affect
only a handful of names. Hard-excluding them systematically biases toward
mediocre, broad-based signals.

Flagging rules:

- $\tilde b_p<0.5$ → `CONCENTRATED`
- $\tilde b^{\rm p10}_p<0.2$ → `EXTREME_CONCENTRATION`

These flags indicate **what to check next**, not an automatic exclusion:

1. Check whether the top-contributing names share a common trait (same
   sector? same event type? same data-vendor anomaly?).
2. Check whether $Q_p$ remains significant after excluding the top 1% of
   contributing names.
3. Estimate the actually tradeable capacity.

### Worked Example

| feature | $b^{\rm med}$ | null $b$ | $\tilde b$ | reading |
|---|---:|---:|---:|---|
| momentum | 0.42 | 0.41 | 1.02 | broad-based, indistinguishable from noise baseline |
| earnings surprise | 0.06 | 0.40 | 0.15 | highly concentrated, flag and check capacity |
| volatility prediction | 0.20 | 0.22 | 0.91 | looks concentrated, actually just the $\mathsf V$ block's natural baseline |

The third row demonstrates why relativization is necessary: at a 0.405
threshold this feature would be mislabeled as concentrated.

---

## §32. Cancellation Ratio

### Definition

$$\boxed{R^{\rm canc}_{d,p}=\frac{\big|\sum_ig_{id,p}\big|}{\sum_i|g_{id,p}|+\epsilon}\in[0,1].}$$

### Division of Labor With Breadth

Both are derived from $\{g_{id}\}$, but ask different questions:

- **breadth**: how spread out are contributions across stocks? (the
  concentration of $|g|$)
- **cancellation**: do the contributions cancel each other out? (the sign
  consistency of $g$)

A feature can have high breadth but low cancellation — all 3000 stocks
contribute, but half positive and half negative, with the net effect coming
from a small imbalance. Such a signal's defining feature: trading costs
will eat almost all of it, because you have to trade all 3000 names to
capture that tiny net value.

### The Null Scale

As with breadth, $R^{\rm canc}=1$ is not the reference point. Under the
null, $g_{id}$ is symmetric with zero mean, so

$$\Big|\sum_ig_{id}\Big|\sim\big|\mathcal N(0,N\sigma_g^2)\big|\ \Longrightarrow\ \mathbb E\Big|\sum_ig\Big|=\sigma_g\sqrt{\frac{2N}\pi},$$

and $\sum_i|g|\approx N\,\mathbb E|g|$. So

$$\mathbb E\big[R^{\rm canc}\big]\approx\frac{\sigma_g}{\mathbb E|g|}\sqrt{\frac2{\pi N}}=\frac\pi2\sqrt{\frac2{\pi N}}=\sqrt{\frac{\pi}{2N}},$$

(the last step uses the product case $g=uv$: $\sigma_g=1$,
$\mathbb E|g|=2/\pi$).

$$N=500\Rightarrow0.056,\qquad N=1500\Rightarrow0.032,\qquad N=3000\Rightarrow0.023.$$

So the raw $R^{\rm canc}$ isn't just far below 1 — it **shrinks with
$N_d$**, and isn't comparable across dates or universes.

### Relative Cancellation

$$\boxed{\tilde R_p=\frac{\operatorname{median}_dR^{\rm canc}_{d,p}}{\operatorname{median}_b\big(\operatorname{median}_dR^{{\rm canc},{\rm null}}_{d,p,b}\big)}.}$$

$\tilde R_p\approx1$ means the net effect is indistinguishable from pure
noise — **a serious warning sign**: it means $Q_p$ may come from a handful
of extreme names rather than a systematic cross-sectional relationship.
$\tilde R_p\ge3$ means there's a genuine, sign-consistent cross-sectional
effect.

The report lists both the raw $R^{\rm canc}$ (for capacity and
cost estimation, since absolute value is what matters for trading
economics) and $\tilde R_p$ (for statistical interpretation).

---

## §33. Effective Spectral Rank

### Definition

$$e_j=\frac{s_j^2}{\sum_ks_k^2},\qquad\boxed{r_{\rm eff}=\frac1{\sum_je_j^2}.}$$

This is the participation ratio of the singular-value energy distribution:
if energy concentrates in a single mode, $r_{\rm eff}\to1$; if $r$ modes
share equal energy, $r_{\rm eff}=r$.

### Use

$r_{\rm eff}$ and $r_{\rm sig}$ answer different questions:

- $r_{\rm sig}$: **how many directions are statistically above noise** (a
  test result).
- $r_{\rm eff}$: **how many directions the energy is actually spread
  across** (a descriptive quantity).

Comparing the two is itself informative:

| Case | Meaning |
|---|---|
| $r_{\rm eff}\approx1\ll r_{\rm sig}$ | the whole library essentially illuminates one mechanism; the rest of the significant modes are weak. The library is too homogeneous — add new data sources. |
| $r_{\rm eff}\approx r_{\rm sig}$ | several mechanisms of comparable strength; good library diversity. |
| $r_{\rm eff}\gg r_{\rm sig}$ | energy spread across many weak directions. Usually means the overall signal is weak, $V_r$'s estimate is unstable, and downstream conclusions should be downgraded. |

$r_{\rm eff}$ also gives a lower-bound reference for feature budget: if
$r_{\rm eff}=3$, selecting 40 features is necessarily highly redundant at
the mechanism level — redundancy itself isn't wrong (spreading noise has
value), but it should be a conscious decision, not a mistaken belief that
40 genuinely different signals were found.
