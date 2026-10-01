---
title: "Cross-Sectional Predictive Geometry, Part 7: Mechanism Similarity and Raw Redundancy"
description: "Cosine similarity between feature coordinate vectors has a precise energy-loss interpretation that sets the right clustering threshold -- but in a 2-4 dimensional block, hundreds of features can't help being close to each other, which is a hard structural limit on what clustering can tell you, not a calibration failure."
date: 2026-09-14
series: ["Cross-Sectional Predictive Geometry"]
part: 7
tags: ["cross-sectional", "methodology", "quant-research", "predictive-geometry", "redundancy"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Predictive Geometry" >}}

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
