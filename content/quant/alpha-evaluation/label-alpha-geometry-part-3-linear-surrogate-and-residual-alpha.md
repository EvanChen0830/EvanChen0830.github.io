---
title: "Label–Alpha Geometry, Part 3: The Linear Surrogate and Residual Alpha"
description: "The real test of whether a model's complexity earned anything: a walk-forward ridge surrogate, two reconstruction metrics that can disagree in informative ways, family-level attribution that avoids the pitfalls of raw regression coefficients, and the residual-alpha diagnostic that this whole method ultimately builds toward."
date: 2026-09-14
series: ["Label-Alpha Geometry"]
part: 3
tags: ["alpha-evaluation", "cross-sectional", "methodology", "quant-research", "model-diagnostics"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Label-Alpha Geometry" >}}

Continues from
[Part 2](/quant/alpha-evaluation/label-alpha-geometry-part-2-significance-and-inheritance/).

## §54. The third-layer problem

The first two layers (§51 output similarity, §52 mechanism similarity) are
both geometric, correlational. Neither can answer the question that
actually matters:

> **Did the model's complexity buy anything?**

The right control for this question is: **the same inputs, a linear
model.** If a ridge regression can reconstruct 95% of the model's output
OOS, the model's nonlinearity is just reproducing a linear function — the
complexity bought nothing, and deploying ridge instead would be far better
on operating cost, interpretability, and stability.

So the third layer's design is: build an OOS linear surrogate, measure how
much of the model it reconstructs, and analyze whether **the part it
couldn't reconstruct** has any predictive power.

## §55. Walk-forward linear surrogate

### Specification

Split the OOS alpha period into consecutive folds over time. For fold $k$,
fit on $Train_k^{\rm sur}$ and predict on $Test_k^{\rm sur}$:

$$\hat\beta_k=\arg\min_\beta\left[\sum_{(d,i)\in Train_k^{\rm sur}}w_{id}\big(\tilde\alpha_{id}-X_{id}\beta\big)^2+\lambda\|\beta\|_2^2\right],\qquad \hat\alpha^{\rm lin}_{id}=X_{id}\hat\beta_k\ \ \text{for }(d,i)\in Test_k^{\rm sur}.$$

Concatenating every fold's predictions gives the complete OOS linear
surrogate $\hat\alpha^{\rm lin}$.

### Four specification points that must be made explicit

**(a) Which $X$ to use.** The primary report uses **the model's actual
input set** (the `selected` set from §39). Reason: the question is "did the
model's nonlinearity add value," and its correct control is "same inputs, a
linear model."

Also report, as a secondary figure, the version using **every eligible
feature**, which answers a different question: "can this alpha be linearly
reproduced by the entire library?" Name the two

$$R^2_{\rm lin,selected},\qquad R^2_{\rm lin,library}.$$

$R^2_{\rm lin,library}\gg R^2_{\rm lin,selected}$ means the model is using
information carried by features that weren't selected — usually meaning §39
was too strict, or the model's actual inputs don't match the selected set
(an implementation bug that should be excluded).

**(b) Standardization.** $X$ uses the §4.1-processed $X^{XS}$ (already
within-date standardized and globally normalized), consistent with the rest
of the pipeline. $\tilde\alpha$ uses the §42 version. Standardizing again
inside the surrogate fit is **not allowed** — it would break the
correspondence between $R^2$'s scale and $Q$'s scale.

**(c) Intercept.** Both $X^{XS}$ and $\tilde\alpha$ are already within-date
demeaned, so the weighted-least-squares intercept is identically zero. **Do
not fit an intercept.** This must be explicitly turned off in code (most
libraries default it on), or it introduces a degree of freedom that absorbs
day-level common components.

**(d) Choosing $\lambda$.** Determine it **inside** $Train_k^{\rm sur}$ via
time-ordered CV (split by time, never randomly), never using any
information from $Test_k^{\rm sur}$. Each fold determines its own $\lambda_k$,
and $\{\lambda_k\}$ is logged in the report — if $\lambda_k$ varies by more
than an order of magnitude across folds, the surrogate itself is unstable,
and subsequent $R^2$ readings should be downgraded accordingly.

### Why ridge, not OLS

Selected features are necessarily correlated (§36 only removes duplicates
above $S^{\rm raw}>0.9$; correlations in the $0.6$–$0.9$ range are common).
OLS coefficients are highly unstable under these conditions — $\hat\beta_k$
would swing wildly between folds, so the surrogate's OOS performance would
be driven by coefficient noise rather than model structure. Ridge trades a
small amount of bias for a large reduction in variance, giving a more
stable answer to "can this be reconstructed."

Note the surrogate's purpose is **not** to estimate each feature's
contribution (that's §57's job), but to estimate **the reconstructible
fraction**. The coefficients themselves are not interpreted.

## §56. Reconstruction metrics

### Definition

$$\boxed{R^2_{\rm lin,OOS}=1-\frac{\sum_{d,i}w_{id}\big(\tilde\alpha_{id}-\hat\alpha^{\rm lin}_{id}\big)^2}{\sum_{d,i}w_{id}\tilde\alpha_{id}^2}=1-\big\|\varepsilon_\alpha\big\|_W^2,}$$

(the last equality uses $\|\tilde\alpha\|_W=1$), and at the ranking level,

$$\rho^{\rm rank}_{\rm lin}=\operatorname{median}_d\operatorname{Spearman}\big(\tilde\alpha_d,\hat\alpha^{\rm lin}_d\big).$$

### The division of labor between the two metrics

$R^2$ measures **numerical** reconstruction; $\rho^{\rm rank}$ measures
**ranking** reconstruction. The two can differ substantially:

- $R^2=0.5$ but $\rho^{\rm rank}=0.9$: the model and the surrogate produce
  almost the same ranking, differing mainly in tail magnitude. Since
  positions are usually determined by ranking, the model's real-world value
  is close to the surrogate's in this case.
- $R^2=0.8$ but $\rho^{\rm rank}=0.6$: numerically close but ranked
  differently — usually means the model and surrogate diverge in the
  middle of the distribution (most names), which happens to be where
  position weight is densest.

**The report must show both side by side**, with $\rho^{\rm rank}_{\rm lin}$
as the primary criterion for "can the model be replaced by the surrogate,"
since it's closer to trading semantics.

### OOS $R^2$ can be negative

If the surrogate fails badly on some folds, $R^2_{\rm lin,OOS}<0$ is
possible, meaning the surrogate is worse than predicting zero everywhere.
This should not be read as "the model is highly nonlinear" — first check
the surrogate's implementation ($\lambda$ too small, leakage or
discontinuity in the fold split).

## §57. Which feature families reconstruct alpha

### Specification

A single ridge coefficient is not a suitable importance measure (correlated
features make coefficients unstable and arbitrarily interchangeable).
Instead use two complementary quantities, **both of which must be
reported.**

**(a) Conditional contribution (drop-family):** for redundancy family
$R_g$, remove every feature in that family, rerun the surrogate through the
identical pipeline, and get $R^2_{-g}$:

$$\boxed{\Delta R^2_g=R^2_{\rm full}-R^2_{-g}.}$$

Meaning: this family's marginal contribution **given that every other
family is present.**

**(b) Standalone contribution:** fit the surrogate using only family
$R_g$'s features:

$$\boxed{R^2_{g,{\rm own}}.}$$

Meaning: how much this family can reconstruct on its own.

### Why both are needed

Two families that are highly correlated with each other (even if
$S^{\rm raw}<0.9$ and so weren't merged) will mask each other:

$$\Delta R^2_g\approx0\quad\text{and}\quad\Delta R^2_{g'}\approx0,\qquad\text{yet}\quad R^2_{g,{\rm own}}=R^2_{g',{\rm own}}=0.6.$$

Looking only at $\Delta R^2$ would conclude "neither family matters," when
the truth is "either one alone suffices, but you can't drop both."

Interpretation rule:

| $\Delta R^2_g$ | $R^2_{g,\rm own}$ | Interpretation |
|---|---|---|
| High | High | Core family, irreplaceable |
| Low | High | Redundant with another family; usable as a backup/substitute |
| High | Low | Only useful in combination with other families (interaction or residual correction) |
| Low | Low | No contribution to linear reconstruction |

### A more complete version

If there aren't too many families ($\le12$), compute family-level Shapley
values: compute $R^2$ for every subset, average the marginal contributions.
Cost is $2^G$ surrogate fits — about 1000 at $G=10$, practically
affordable. Shapley values resolve the masking problem above and give the
cleanest attribution, but don't scale to large $G$.

## §58. Residual alpha

### Definition

$$\boxed{\varepsilon_\alpha=\tilde\alpha-\hat\alpha^{\rm lin}.}$$

Since both $\tilde\alpha$ and $X^{XS}$ are already within-date demeaned,
$\varepsilon_\alpha$ automatically satisfies $\sum_i\varepsilon_{\alpha,id}=0$.

**Two scales must be recorded separately.**

**Absolute residual energy:**

$$s^2_\varepsilon=\|\varepsilon_\alpha\|_W^2=1-R^2_{\rm lin,OOS}.$$

**Standardized residual signal:** within-date standardize $\varepsilon_\alpha$
and globally normalize, giving $\tilde\varepsilon_\alpha$ with
$\|\tilde\varepsilon_\alpha\|_W=1$, then compute per §43

$$c^{\rm full}_\varepsilon=\tilde\varepsilon_\alpha^\top W\Phi,\qquad c^{\rm inh}_\varepsilon=c^{\rm full}_\varepsilon P_V,\qquad c^{\rm nov}_\varepsilon=c^{\rm full}_\varepsilon(I-P_V),$$

$$Q_\varepsilon=\|c^{\rm full}_\varepsilon\|^2,\quad Q^{\mathsf D}_\varepsilon,\ Q^{\mathsf V}_\varepsilon,\ Q^{\mathsf L}_\varepsilon,\quad A^{\rm tgt}_\varepsilon,\ \eta^{\mathsf D}_\varepsilon,\quad\tilde\nu_\varepsilon.$$

### Why $Q_\varepsilon$ must be computed after standardizing

Using $\varepsilon_\alpha$ at its raw scale gives $\|\varepsilon\|_W^2=1-R^2$, so

$$Q_\varepsilon^{\rm raw}\le1-R^2.$$

At $R^2=0.9$, $Q_\varepsilon^{\rm raw}$'s upper bound is only 0.1 times
$Q_\alpha$'s upper bound. Comparing the two numbers directly would conclude
"the residual has no predictive power," when that's purely scale
arithmetic, unrelated to the residual's actual quality.

After standardizing, $Q_\varepsilon$ and $Q_\alpha$ sit on **the same
scale** (the IC scale of §11) and can be compared directly. The relationship
between them is

$$\boxed{\underbrace{s^2_\varepsilon\cdot Q_\varepsilon}_{\text{residual's absolute contribution to total predictive energy}}\quad\text{vs}\quad\underbrace{Q_\alpha}_{\text{total predictive energy}}.}$$

**The report must list both $Q_\varepsilon$ (quality) and
$s^2_\varepsilon Q_\varepsilon$ (absolute contribution).** The former
answers "is the residual signal good"; the latter answers "how much of the
total it accounts for."

### Significance

Test $\tilde\varepsilon_\alpha$ using the same §13 permutation set, giving
$p^{\mathsf D}_\varepsilon$, $p^{\rm nov}_\varepsilon$, etc.

Worth noting: even though standardization is what makes $Q_\varepsilon$'s
**magnitude** comparable, the permutation p-value itself was always
scale-invariant — it compares the observed value of the same standardized
signal against its own null distribution. So the significance conclusion is
identical before and after standardizing; standardizing only changes which
magnitudes are comparable.

### This section is this method's final product

The combination of $Q^{\mathsf D}_\varepsilon$ and $\tilde\nu_\varepsilon$
answers a question that's extremely hard to answer in quantitative research
generally, and that this method answers directly:

> **Does the part of the model's output that lies outside the linear
> feature span have any predictive value at all? And if it does, is it a
> re-weighting of known mechanisms, or an entirely new one?**

$Q^{\mathsf D}_\varepsilon$ significance answers the first half;
$\tilde\nu_\varepsilon$ answers the second.
