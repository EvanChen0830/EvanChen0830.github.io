---
title: "Eligibility, Grouping, and Submodular Coverage: the Formal Selection Layer"
description: "The formal specification behind feature selection: five independence-aware eligibility criteria, a two-tier mechanism/redundancy grouping scheme, a submodularity theorem that justifies greedy coverage, and the direction-first selection algorithm."
date: 2026-09-14
tags: ["feature-selection", "cross-sectional", "methodology", "quant-research", "submodularity"]
math: true
ShowToc: true
TocOpen: true
---

This post is the formal counterpart to the
[Cross-Sectional Feature Selection](/quant/feature-evaluation/cross-sectional-feature-selection/)
guide. That guide is an implementation reference: it states the selection
procedure (eligibility floor, mechanism clustering, redundancy families,
coverage-based greedy fill, direction-first budget protection) at the level
a practitioner needs to build it. What follows is the formal specification
those steps are built on — five independence-aware eligibility criteria
instead of a quality floor, a two-tier grouping scheme with an explicit
energy-loss justification for its threshold, a submodularity theorem that
turns "greedily add features" from a heuristic into a certified
approximation algorithm, and the exact direction-first selection algorithm
with its termination conditions spelled out. Where the correspondence is
exact, this post says so explicitly rather than restating the informal
version; where it goes further, that's the new material.

**Correspondence at a glance:**

| Here | Cross-Sectional Feature Selection series |
|---|---|
| §34 Eligibility | Part 1, step B9 (quality floor) — this is the rigorous version: five named criteria with an explicit "why it can't be dropped" for each, not one floor |
| §35 Mechanism grouping | Part 1, step B9 (mechanism clustering) |
| §36 Raw redundancy families | Part 1, step B10 |
| §37 Family representative | Part 1, step B13a-b (Score, Stage 1) |
| §38 Coverage | Part 1, step B13c — this is where the real gap is: the series states the greedy rule, this section proves it's a certified approximation |
| §39 Direction-first selection | Part 4 (Protecting Directional Exposure) — the series gives five practitioner levers; this is the single canonical algorithm those levers are informal variations of |
| §40 Feature report format | Not previously published — the complete field schema for a production feature report |

Three stages, three independent questions, never mixed:

1. **Eligibility (§34)** — is this feature real? A binary determination from
   significance and stability.
2. **Grouping (§35–§36)** — which features are saying the same thing? A
   descriptive organization for human interpretation.
3. **Coverage selection (§37–§39)** — which subset should represent the
   whole? An optimization problem, decided by submodular coverage.

The three must not be conflated. In particular: **grouping does not decide
selection** — selection is decided entirely by the coverage criterion in
stage three.

## §34. Eligibility

A feature enters selection only if all five of the following hold (for a
directional research target):

$$\boxed{
\begin{aligned}
&\text{(E1) } \tilde p^{\,\mathsf D}_p<0.05 &&\text{(Westfall–Young adjusted, §18)}\\
&\text{(E2) } \tilde p^{\,\rm tgt}_p<0.05 &&\text{(test on }(A^{\rm tgt}_p)^2\text{, §21)}\\
&\text{(E3) } \pi^{\rm time,\mathsf D}_p\ge0.8\ \text{and}\ \pi^{\rm time,sign}_p\ge0.8 &&\text{(§27)}\\
&\text{(E4) } \pi^{\rm XS,\mathsf D}_p\ge0.8 &&\text{(§29)}\\
&\text{(E5) } M^{\rm time,\mathsf D}_p\ge0.6\ \text{and}\ M^{\rm XS,\mathsf D}_p\ge0.6 &&\text{(§27, §29)}
\end{aligned}}$$

Breadth and cancellation are soft diagnostics: flagged, not excluded (§31.3).

**What each threshold excludes:**

| Threshold | Failure mode excluded |
|---|---|
| E1 | Purely coincidental association |
| E2 | Has directional information, but pointed at a target orthogonal to the one in question |
| E3 | Only works in a few periods; sign flips |
| E4 | Driven by a handful of names |
| E5 | The mechanism itself is drifting (§28, failure mode two) |

### Why E2 cannot be dropped

E1 only requires "has directional predictive power." A feature that is
significant, stable, and broad on the $\mathsf D$ block but has
$\eta^{\mathsf D}_p\approx0$ passes E1, E3, E4, and E5 cleanly — it could
even become a family representative. What it predicts is a different
horizon combination within the band: mathematically legitimate, useless for
the actual target at hand.

E2 tests $(A^{\rm tgt}_p)^2=Q^{\mathsf D}_p\cdot\eta^{\mathsf D}_p$ directly,
which simultaneously covers "how much directional information" and "how
well-aimed is it." This is better than a fixed threshold like
$\eta^{\mathsf D}_p\ge\eta_{\min}$, because $\eta$ itself is unstable for
weak features (§22: $\eta$ is meaningless once $Q^{\mathsf D}$ approaches
noise level), while E2's test statistic has its own null distribution by
construction.

**Features excluded by E2 should not be discarded.** Flag them
`TARGET_MISMATCH` and route them to evaluation against other bands or other
alignment targets $a$. They're still library assets — just not matched to
the current target.

### The thresholds are not independent

E1 through E5 are highly correlated: a feature with a tiny $\tilde p$ will
almost certainly have high $\pi$. The **joint** filtering power of these
five thresholds is therefore much weaker than the product of five
independent thresholds would suggest. This isn't a flaw, but it has to be
kept in mind when reading a pass rate: if 60 of 400 features pass all five,
that does not mean "60 features each independently survived five separate
tests." Genuine multiple-testing control only happens at E1 and E2 (both
Westfall–Young adjusted).

### Other research targets

For a volatility-signal research target, replace every $\mathsf D$ with
$\mathsf V$ and drop E2 (the $\mathsf V$ block has no notion of target
direction — it isn't a signed directional space). For large-move
directionality, replace with $\mathsf L$ and keep E2 ($\mathsf L$ *is* a
signed directional space; the analogue of $r_{\mathsf D}$ is
$r_{\mathsf L}=\text{normalize}(\Phi_{\mathsf L}^\top Wy^*)$).

## §35. Mechanism grouping

For the eligible set $\mathcal E$, run **average-linkage** hierarchical
clustering with distance $1-S^{\mathsf D}_{pq}$, cut at

$$\tau_{\mathsf D}=0.90\qquad\text{(i.e. }S^{\mathsf D}\ge0.90\text{ counts as the same group)},$$

to get mechanism groups $G_1,\dots,G_M$.

The $\tau_{\mathsf D}=0.90$ choice follows the energy-loss criterion of
§23.1: replacing a member with its representative loses at most
$1-0.9^2=19\%$ of mechanism energy.

**Why average linkage, not complete or single:** single linkage produces
chaining ($p\sim q\sim r$ but $p\not\sim r$), especially severe in a
low-dimensional space; complete linkage is overly sensitive to a single
outlier pair. Average linkage sits between the two and is the most robust
to this problem's distance structure.

### Calibration quantities that must accompany the grouping

By Rule 23.2, a grouping result **must** be reported alongside two
quantities:

$$\zeta(\tau_{\mathsf D})\ \text{(excess-pairing rate, §23.2)},\qquad m_{\mathsf D}\ \text{(effective mechanism dimension, §23.4)}.$$

If $\zeta\le0$ or $m_{\mathsf D}<4$, the report must carry this note:

> **Mechanism grouping has no discriminative power under this setting; the
> group boundaries are arbitrary and must not be cited as a conclusion.**

This is the norm, not the exception, whenever $r_{\rm sig}\le3$. Grouping
can still be used for visual organization of the report in that case, but
every substantive decision is deferred to §38–§39.

### Group structure is itself a library diagnostic

- Group count $M$ much smaller than $|\mathcal E|$: the library is severely
  over-sampling a small number of mechanisms. This is a warning sign for
  feature construction — add a new data source rather than more variants of
  the same feature type.
- A singleton group: that feature carries a unique mechanism, and will
  automatically get high priority in the §39 coverage selection.
- The largest group exceeds half of $|\mathcal E|$: usually means the
  library is dominated by a single data vendor or a single feature type
  (e.g. everything is price-derived).

## §36. Raw redundancy families

**Within** each mechanism group $G_m$, run **complete-linkage** clustering
on $S^{\rm raw}_{pq}$ (§25), cut at $S^{\rm raw}=0.90$, to get families
$R_{m1}, R_{m2}, \dots$

A family therefore means: **close in predictive mechanism, and also highly
redundant in raw cross-sectional information.**

### The logic of the two-tier structure

- **Raw redundancy is only checked within a mechanism group.** Two features
  that are highly correlated in raw terms but predict different mechanisms
  (e.g. the raw value predicts direction, its absolute value predicts
  activity) should both be kept.
- **Never merge across groups.** Even at $S^{\rm raw}=0.95$, different
  mechanisms means different families.

This order is correct because the actual goal is preserving **diversity of
predictive mechanism**; raw redundancy is only a secondary criterion for
picking a representative once the mechanism is already shared.

Using complete linkage here (versus average linkage in §35) is deliberate:
raw-redundancy judgments should be conservative — complete linkage requires
**every** pair within a group to clear the threshold, which avoids merging
features that are only indirectly related.

### The transform-family rule

By Rule 4.1, every transform of the same source feature,
$(p,1),(p,2),\dots$, is **forced into the same raw family**, unless they
fall into different mechanism groups.

In implementation this is a post-processing step: after §36's clustering,
scan each source feature, and if its transforms land in different families
within the same group, merge those families.

### Reporting

Record each feature's `mechanism_group`, `raw_family`, and the sizes of
both tiers. If a mechanism group contains 30 features but splits into only
2 raw families, that mechanism is being severely over-sampled.

## §37. Family representative

Each raw family $R_{mk}$ selects one representative:

$$\boxed{\text{Score}_p=\operatorname{median}_b\left[\big(A^{{\rm tgt},(b)}_p\big)^2\right],\qquad p^*=\arg\max_{p\in R_{mk}}\text{Score}_p,}$$

where $b$ ranges over the temporal bootstrap resamples of §26.

If the research target doesn't specify a direction $a$ (exploratory
research), use $\text{Score}_p=\operatorname{median}_b\big[Q^{\mathsf D,(b)}_p\big]$
instead.

### Why a bootstrap median, not a point estimate

Using the point estimate of $(A^{\rm tgt}_p)^2$ directly runs into
**winner's curse**: within a family, the feature with the largest point
estimate tends to be the one that's "genuinely above-average, and
especially lucky." The bootstrap median has a built-in shrinkage against
this effect — a feature propped up by a handful of periods has a heavily
right-skewed bootstrap distribution, so its median sits well below its
point estimate.

### Why not a multi-factor product score

An intuitively natural approach is to multiply every quality metric
together:

$$Q^{\mathsf D}_p\cdot\pi^{\rm time}_p\cdot\pi^{\rm sign}_p\cdot M^{\rm time}_p\cdot\pi^{\rm XS}_p\cdot M^{\rm XS}_p\cdot h(b_p).$$

This has two problems.

**First, it double-counts eligibility.** Anything reaching §37 already
satisfies $\pi\ge0.8$ (three terms) and $M\ge0.6$ (two terms). So these
factors only actually range over $[0.8,1]$ and $[0.6,1]$ respectively — the
six factors combined can move the product by at most roughly
$(1/0.8)^3(1/0.6)^2\approx5.4\times$.

**Second, $Q^{\mathsf D}$ ranges over more than 100×** (from $10^{-5}$ to
$10^{-3}$). So the product score is numerically dominated almost entirely
by $Q^{\mathsf D}$, and the other six factors are little more than
noise-level adjustments. Writing them into the score only manufactures an
appearance of thoroughness, while actually adding untraceable degrees of
freedom.

The correct design: **put quality thresholds in eligibility (a hard,
traceable gate), and put ranking in a single quantity with unambiguous
meaning (Score).**

### The role of breadth

$\tilde b_p$ does not enter the Score. It's used as a **tiebreaker** and a
**post-hoc check**: if the top two Scores within a family are within 20% of
each other, pick the one with higher $\tilde b$. The reasoning (§31.3):
concentration isn't itself a defect, but all else being roughly equal, a
broad-based signal has better capacity and robustness.

## §38. Coverage

This section defines two distinct notions of coverage. **They answer
different questions and are not interchangeable.**

### §38.1 Span coverage (diagnostic only)

$$P^{\mathsf B}_{\mathcal F}=C^{\mathsf B\top}_{\mathcal F}\big(C^{\mathsf B}_{\mathcal F}C^{\mathsf B\top}_{\mathcal F}\big)^+C^{\mathsf B}_{\mathcal F},\qquad\boxed{{\rm Cov}^{\rm span}_{\mathsf B}(\mathcal F)=\frac{\big\|C^{\mathsf B}_{\mathcal E}P^{\mathsf B}_{\mathcal F}\big\|_F^2}{\big\|C^{\mathsf B}_{\mathcal E}\big\|_F^2}\in[0,1],}$$

where $C^{\mathsf B}_{\mathcal F}$ stacks the mechanism vectors (block
$\mathsf B$ part) of the selected set $\mathcal F$, and $C^{\mathsf B}_{\mathcal E}$
does the same for every eligible feature.

**Proposition 38.1 (monotonicity).** $\mathcal F\subseteq\mathcal F'\Rightarrow{\rm Cov}^{\rm span}_{\mathsf B}(\mathcal F)\le{\rm Cov}^{\rm span}_{\mathsf B}(\mathcal F')$.

*Proof.* $\operatorname{span}(C_{\mathcal F})\subseteq\operatorname{span}(C_{\mathcal F'})$,
and the norm of a projection onto a larger subspace does not decrease. $\blacksquare$

**Proposition 38.2 (saturation).** $c_p^{\mathsf D}\in\mathbb R^{d_{\mathsf D}}=\mathbb R^K$.
So as soon as $\mathcal F$ contains $K$ linearly independent $c^{\mathsf D}$
vectors, $P^{\mathsf D}_{\mathcal F}=I$, and

$$\boxed{{\rm Cov}^{\rm span}_{\mathsf D}(\mathcal F)=1.}$$

At $K=3$, **any three linearly independent features already reach 100%
coverage.**

Proposition 38.2 makes span coverage entirely unsuitable as a selection
criterion: a 90% target is already met at the third or fourth feature
chosen, the greedy step stops immediately, and the three features selected
could be nearly arbitrary. Worse, this criterion structurally **cannot**
select more than $d_{\mathsf D}$ features, even if the library has 50
genuinely distinct, independently worth-keeping mechanisms.

> **Rule 38.1.** ${\rm Cov}^{\rm span}$ is diagnostic-only: it answers
> "does the selected set span every mechanism direction." It **does not
> enter** the §39 selection loop.

It remains useful diagnostically: if ${\rm Cov}^{\rm span}_{\mathsf D}(\mathcal F)<1$,
the selected set doesn't even span the $\mathsf D$ space — usually meaning
every selected feature points at the same horizon combination, which is a
problem worth investigating.

### §38.2 Feature coverage (selection criterion)

The actual question selection needs to answer is: **does every feature
worth keeping have a selected representative that's "close enough" to it?**

$$\boxed{{\rm Cov}^{\rm feat}_{\mathsf B}(\mathcal F)=\frac{\sum_{p\in\mathcal E}w_p\max_{j\in\mathcal F}S^{\mathsf B}_{pj}}{\sum_{p\in\mathcal E}w_p},\qquad w_p=\text{Score}_p\ \text{(§37)}.}$$

(By convention $\max_{j\in\varnothing}=0$, and $S^{\mathsf B}_{pp}=1$ when $p\in\mathcal F$.)

Interpretation: the Score-weighted average "representedness." ${\rm Cov}^{\rm feat}=0.9$
means that, weighted by importance, the eligible features' mechanisms have
on average 0.9 cosine similarity captured by some selected representative
— corresponding to $0.9^2=81\%$ of mechanism energy.

### §38.3 Theorem 38.3 (submodularity)

${\rm Cov}^{\rm feat}_{\mathsf B}$ is a monotone non-decreasing,
**submodular** set function.

*Proof.* For fixed $p$, let $h_p(\mathcal F)=\max_{j\in\mathcal F}S^{\mathsf B}_{pj}$.

*Monotonicity:* $\mathcal F\subseteq\mathcal F'$ widens the range the maximum
is taken over, so $h_p(\mathcal F)\le h_p(\mathcal F')$.

*Submodularity:* let $\mathcal A\subseteq\mathcal B$, $j\notin\mathcal B$. Then

$$h_p(\mathcal A\cup\{j\})-h_p(\mathcal A)=\max\Big(0,\ S^{\mathsf B}_{pj}-\max_{k\in\mathcal A}S^{\mathsf B}_{pk}\Big)\ \ge\ \max\Big(0,\ S^{\mathsf B}_{pj}-\max_{k\in\mathcal B}S^{\mathsf B}_{pk}\Big)=h_p(\mathcal B\cup\{j\})-h_p(\mathcal B),$$

which holds because $\max_{k\in\mathcal A}\le\max_{k\in\mathcal B}$.

Finally, ${\rm Cov}^{\rm feat}=\frac1{\sum w}\sum_pw_ph_p$ is a non-negative
weighted sum, and monotone submodular functions are closed under
non-negative linear combination. $\blacksquare$

**Corollary 38.4 (greedy guarantee).** Under a cardinality constraint
$|\mathcal F|\le n$, the greedy algorithm attains at least
$1-1/e\approx63\%$ of the optimum (Nemhauser–Wolsey–Fisher). In the "reach
coverage level $\theta$" formulation, greedy needs no more than
$O(\log\frac1{1-\theta})$ times the optimal number of elements (Wolsey's
submodular cover theorem).

This guarantee is ${\rm Cov}^{\rm feat}$'s substantive advantage over
${\rm Cov}^{\rm span}$: it turns "greedily add one at a time" from a
plausible heuristic into something with a formal quality guarantee.

### §38.4 Comparing the two coverage notions

| | ${\rm Cov}^{\rm span}$ | ${\rm Cov}^{\rm feat}$ |
|---|---|---|
| Question | Does it span the mechanism space? | Does every feature have a representative? |
| Upper bound set by | $d_{\mathsf B}$ ($K$ for $\mathsf D$) | No structural bound |
| Features needed to reach 1 | $K$ | All of them |
| Submodular | No (matroid-rank type, only monotone) | Yes |
| Use | Diagnostic | Selection |

## §39. Direction-first selection

**Step 1: seed.** $\mathcal F\leftarrow\{$the representative $p^*$ of every
raw family$\}$ (§37).

**Step 2: directional coverage.** Greedily add from the remaining eligible
features

$$p^\dagger=\arg\max_{p\notin\mathcal F}\Big[{\rm Cov}^{\rm feat}_{\mathsf D}(\mathcal F\cup\{p\})-{\rm Cov}^{\rm feat}_{\mathsf D}(\mathcal F)\Big],$$

until

$$\boxed{{\rm Cov}^{\rm feat}_{\mathsf D}(\mathcal F)\ge0.90}$$

or the feature budget $n_{\max}$ is reached.

**Step 3: fill $\mathsf V$ and $\mathsf L$.** If budget remains, repeat step
2 for the $\mathsf V$ and $\mathsf L$ blocks (using $S^{\mathsf V}$,
$S^{\mathsf L}$ and their respective Scores), target $0.80$.

**Step 4: diagnose.** Compute and record ${\rm Cov}^{\rm span}_{\mathsf D}(\mathcal F)$,
${\rm Cov}^{\rm span}_{\mathsf V}(\mathcal F)$, ${\rm Cov}^{\rm span}_{\mathsf L}(\mathcal F)$.
If ${\rm Cov}^{\rm span}_{\mathsf D}<1$, examine which direction isn't spanned.

### Why direction-first

The ordering principle — satisfy direction first, fill the other blocks
after — rests on three points:

1. The downstream objective is cross-sectional relative-strength ranking,
   which depends directly on $\mathsf D$.
2. $\mathsf V$ and $\mathsf L$ information is valuable at the portfolio
   level (position sizing, risk budgeting, conditional overlays), but
   **cannot substitute for direction.** A portfolio with only $\mathsf V$
   information cannot construct a long-short position.
3. If all three blocks compete within a single coverage objective,
   $\mathsf V$ mechanically claims more slots, since $d_{\mathsf V}>d_{\mathsf D}$
   (the dimension table in §6).

### Setting the budget

$n_{\max}$ should be set by the downstream model's sample size, not by this
method. A practical reference: for a tree model, $n_{\max}\approx\sqrt{T^{\mathsf D}_{\rm eff}}$
typically avoids severe overfitting; for a linear or shallow model, more is
fine.

Comparing $n_{\max}$ against $r_{\rm eff}$ (§33) is an important check: if
$r_{\rm eff}=3$ but $n_{\max}=40$, those 40 features are necessarily highly
redundant at the mechanism level. That can be intentional (multiple
measurements to reduce noise), but it should be a conscious decision.

### When the termination condition can't be met

If greedy runs through every eligible feature and still hasn't reached
0.90, the eligible set has genuinely high mechanism diversity — every
feature carries some part that nothing else represents. In that case,
record the coverage actually reached and accept it, while checking whether
$\tau_{\mathsf D}$ (§35) was set too high.

## §40. Feature report format

One row per feature (one row per transform entry $(p,k)$ within a transform
family):

**Identity and data quality**

| Field | Description |
|---|---|
| `feature_id`, `source_id`, `transform` | Entry, source feature, transform |
| `coverage`, `active_ratio` | $\operatorname{median}_dm_{dp}$, $T_p/T$ |
| `norm_before_renorm` | $\|X_p\|_W^2$ before step three of §4.1 |
| `kappa` | Persistence $\kappa_p$ (§13.4) |
| `null_used` | `N1_shift` or `N2_stratified` |

**Strength and mechanism**

| Field | Description |
|---|---|
| `Q`, `Q_D`, `Q_V`, `Q_L` | §16 |
| `Q_full` | $\|c_p^{\rm full}\|^2$ not projected onto $V_r$ (§16, closing note) |
| `f_D`, `f_V`, `f_L` | Raw fingerprint (§19.1) |
| `ft_D`, `ft_V`, `ft_L` | Enrichment fingerprint $\tilde f$ (§19.3) |
| `A_target`, `eta_D` | §21, §22 |

**Significance**

| Field | Description |
|---|---|
| `p_D_raw`, `p_target_raw` | Raw permutation p-value |
| `p_D_wy`, `p_target_wy` | Westfall–Young adjusted (§18.3) |
| `fdr_selected` | Whether selected under the permutation-FDR variant (§18.4) |
| `p_full` | p-value of $Q^{\rm full}$ |

**Stability**

| Field | Description |
|---|---|
| `pi_time_D`, `pi_time_sign`, `M_time_D`, `CV_Q_D` | §27 |
| `pi_XS_D`, `M_XS_D` | §29 |

**Cross-sectional support**

| Field | Description |
|---|---|
| `b_med`, `b_p10` | §30 |
| `b_tilde` | Relative breadth (§31.2) |
| `R_cancel_med`, `R_tilde` | §32 |
| `flags` | `CONCENTRATED` / `EXTREME_CONCENTRATION` / `TARGET_MISMATCH` |

**Structure and selection**

| Field | Description |
|---|---|
| `mechanism_group`, `raw_family` | §35, §36 |
| `is_representative` | §37 |
| `score` | Score from §37 |
| `selected` | Final outcome of §39 |
| `marginal_cov_gain` | $\Delta{\rm Cov}^{\rm feat}_{\mathsf D}$ at the moment of selection |

**Global fields (one row per run)**

$T$, $T^{\mathsf D}_{\rm eff}$, $T^{\mathsf V}_{\rm eff}$, $K$, $K_{\rm eff}$,
$\kappa(\Sigma_Y)$, $\{\lambda_j/\lambda_1\}$, band definition (cumulative or
incremental), $r_{\rm sig}$, $r_{\rm eff}$, $m_{\mathsf D}$,
$\zeta(\tau_{\mathsf D})$, $\tau_{\mathsf D}$, $B_{\rm null}$, $B_{\rm time}$,
$B_{\rm XS}$, inference specification (split or refit method), date ranges
of $\mathcal D$ and $\mathcal E$, $\delta_{\rm drift}(\mathcal D\to\mathcal E)$,
thinning ratio, ${\rm Cov}^{\rm feat}_{\mathsf D}$, ${\rm Cov}^{\rm span}_{\mathsf D}$.
