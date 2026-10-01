---
title: "Cross-Sectional Predictive Geometry, Part 1: Problem Statement, Notation, and Data Contract"
description: "Why collapsing a signal's relationship to future returns into a single rank IC hides three economically different kinds of predictive power — and the geometric fix: lift returns into a finite-dimensional space and give every signal a coordinate vector in it."
date: 2026-09-14
series: ["Cross-Sectional Predictive Geometry"]
part: 1
tags: ["cross-sectional", "methodology", "quant-research", "predictive-geometry", "feature-selection"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Predictive Geometry" >}}

> This series is a from-scratch formal treatment of the same underlying
> framework as the [Cross-Sectional Feature Selection](/quant/feature-evaluation/part-1-setup-and-matrix-mechanics/)
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
