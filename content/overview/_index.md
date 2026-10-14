---
title: "Overview"
description: "The mindset behind this research, and a high-level map of the methodology that the rest of the site goes into section by section."
math: true
ShowToc: true
TocOpen: true
---

## Mindset

Most of what's wrong with a trading signal isn't visible in its headline
number. A feature with a respectable IC can be carried by three large-cap
names. A model can look better on a 42-label objective purely because the
auxiliary labels are easier to rank, not because it learned anything more
useful. A backtest can look clean because the same data that chose the
features also graded them. None of these failures show up unless you go
looking for them specifically — the summary statistic is, by construction,
the thing that hides them.

So the methodology here is built around a small number of standing habits,
applied at every stage from raw features to the final portfolio:

- **Discovery and evaluation are different datasets.** Anything learned,
  chosen, or tuned on one window is not allowed to grade itself — a
  validation split that participated in a selection decision is no longer
  an independent estimate of that decision's quality.
- **Every signal gets a null.** A statistic is not "significant" because
  it's positive; it's significant because it clears a distribution of what
  the same procedure produces on data with the effect deliberately
  destroyed (block-permuted dates, shuffled universes). The null is
  computed with the same pipeline, not a generic significance table.
- **Effective sample size, not raw sample size.** A cross-section of 500
  names on a day when one sector moves together doesn't contain 500
  independent observations. Quantities like $L_{\rm eff}$ and $T_{\rm eff}$
  (inverse-Herfindahl concentration of a statistic's own contributions) are
  tracked specifically so a large-$N$ result can't disguise a small-$N$
  effect.
- **Ask what a number is actually measuring, not just whether it's big.**
  "Explains variance" and "predicts direction" are different claims, and a
  feature can be strong on one and worthless on the other. The decomposition
  below exists to keep those claims from being conflated.
- **The same discipline applies regardless of model family.** A linear
  factor, a gradient-boosted tree, a RankIC transformer, and a diffusion
  model are different ways of answering the same question, and get held to
  the same discovery/evaluation boundary, the same null calibration, and the
  same stability checks.

## The pipeline

Everything on this site sits at one stage of the same pipeline. A signal is
only as good as its weakest stage, which is why each stage is evaluated on
its own terms rather than judged only by the final portfolio result.

{{< pipeline >}}

- **Universe** — which assets are even eligible, and whether that
  eligibility is defined without hindsight.
- **Feature evaluation** — which features actually carry
  cross-sectional predictive signal, survive a null, and aren't redundant
  with something already kept.
- **Model** — generative modeling of multi-horizon cross-sectional return
  distributions (e.g. conditional diffusion) rather than point predictions.
- **Alpha evaluation** — evaluating model predictions with the same
  scrutiny as raw features: significance, stability, and residual novelty
  beyond linear factor surrogates.

## One picture: direction, activity, tail

A large fraction of the methodology reduces to one geometric idea, reused at
the feature layer and the alpha layer alike: **any cross-sectional
predictive relationship decomposes into three orthogonal mechanisms.**

Take a future return $Y$ at a given horizon, stripped of its market-wide
component by within-date demeaning, and whitened against its own
covariance to get $Z$. For a single cross-section (one asset $i$ on one
date), lift $Z$ into three blocks:

$$
\Phi_{\mathsf D} = Z, \qquad
\Phi_{\mathsf V} = \operatorname{svec}\!\left(ZZ^\top - I\right), \qquad
\Phi_{\mathsf T} = \|Z\|_2\, Z.
$$

- $\Phi_{\mathsf D}$ (**Direction**) is just the whitened return itself — it asks
  whether a feature predicts *who outperforms*.
- $\Phi_{\mathsf V}$ (**Activity**) is the outer-product residual — it asks whether a
  feature predicts *dispersion and co-movement structure*, independent of
  which side wins.
- $\Phi_{\mathsf T}$ (**Tail**) scales direction by its own magnitude — it asks
  whether a feature's predictive power is conditional on *large moves*.

{{< svgfig src="feature-eval/fig-mechanism-blocks.svg" >}}
**Figure 1.** The cross-sectional mechanism decomposition. Raw predictive signal separates additively into Direction ($\mathsf D$), Activity ($\mathsf V$), and Tail ($\mathsf T$) channels via sequential Gram-Schmidt orthogonalization.
{{< /svgfig >}}

After orthogonalizing these three blocks against each other, a feature's
(or an alpha's) entire predictive relationship with the target is a single
cross-covariance,

$$
C = X^{\top} W \Phi \in \mathbb{R}^{P \times d_\Phi},
$$

and its quality splits additively across the three mechanisms:

$$
Q = Q^{\mathsf D} + Q^{\mathsf V} + Q^{\mathsf T}, \qquad Q^{\mathsf B} = \|c^{\mathsf B}\|_2^2 \ \ (\mathsf B \in \{\mathsf D, \mathsf V, \mathsf T\}).
$$

That single decomposition answers the question that headline metrics hide:
*not* "is this feature good," but "what is it actually good at, and is that
what the portfolio downstream of it needs." A feature with high $Q_V$ and
negligible $Q_D$ is a volatility-timing signal dressed up as alpha; using it
for directional position sizing is a category error even though its raw
quality score looks fine. The [feature evaluation](/quant/feature-evaluation/)
and model-diagnostic write-ups go through exactly this decomposition in full
mechanical detail.

## Where the generative/diffusion approach fits

Every method above ultimately estimates a point summary: a rank, an
expected return, a quality score. The [diffusion](/quant/model/) line of
work asks a different question: instead of predicting one number, learn the
*entire conditional return distribution* $p(Y \mid X)$, and sample from it.

That turns out to be the same three questions in a different form. Once you
have a learned distribution instead of a point estimate, its first moment
*is* the direction signal, its spread and co-movement structure *is* the
activity signal, and its tail mass *is* the tail signal — they're no longer
separate diagnostics bolted onto a point prediction, they're different
slices of one object you can actually look at. That's the sense in which
"seeing the return distribution" is a continuation of the same mindset
rather than a separate one.
