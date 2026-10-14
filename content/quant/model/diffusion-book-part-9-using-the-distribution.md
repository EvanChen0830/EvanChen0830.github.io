---
title: "Diffusion Implementation Book, Part 9: Using the Distribution"
description: "Turning generated samples into trading decisions: conditional mean–variance optimization, exponential utility over the full empirical distribution, and clustering tickers by their return distributions."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 9
weight: 19
tags: ["diffusion", "generative-model", "portfolio-construction", "clustering", "implementation", "quant-research"]
math: true
---

Two portfolio methods follow directly. The first estimates mean and covariance from the generated distribution and uses them in mean–variance. The second maximizes exponential utility over the whole generated distribution; under a Gaussian assumption it reduces to mean–variance.

{{< svgfig src="diffusion-book/fig11.svg" >}}
**Figure 11.** Three uses of one date's samples at one horizon.
{{< /svgfig >}}

Choose a horizon and take the $K$ generated joint return vectors for a date, each of length $N$. Portfolio return is linear in simple returns, so if the labels are log returns, convert each sample first.

### Mean–variance

Estimate the mean vector as the average of the samples and the covariance as their sample covariance, then solve:

$$
\max_{w \in \mathcal{W}} \; w^\top\mu - \frac{\gamma}{2} \, w^\top\Sigma w
$$

$\gamma > 0$ is risk aversion. The allowed set $\mathcal{W}$ holds your portfolio constraints.

The key point is the covariance: **$\Sigma$ comes from the generated conditional distribution.** It doesn't have to be estimated from historical volatility. The features condition the distribution, and the distribution implies today's covariance. This only works if samples stay paired across tickers; sampling tickers independently throws the covariance away.

### Exponential utility over the full distribution

$$
\max_{w \in \mathcal{W}} \; \frac{1}{K}\sum_{k=1}^K \big[-\exp(-\gamma \cdot w^\top \hat y^{(k)})\big]
$$

$$
\text{equivalently, minimize:} \quad \frac{1}{K}\sum_{k=1}^K \exp(-\gamma \cdot w^\top \hat y^{(k)})
$$

This is what "take the exponential of minus gamma times the distribution" means precisely: evaluate the exponential on every sampled portfolio return, then average. Skewness, multiple peaks and tail outcomes can change the answer even when two portfolios have similar mean and variance. For numerical stability, optimize the equivalent log objective with log-sum-exp. The underlying expectation must exist: a heavy-tailed population may lack a finite exponential moment even though any finite sample gives a finite number.

**Why it becomes mean–variance.** If returns are Gaussian with mean $\mu$ and covariance $\Sigma$, then the portfolio return is Gaussian with mean $w^\top\mu$ and variance $w^\top\Sigma w$, and

$$
\mathbb{E}\big[\exp(-\gamma \cdot w^\top Y)\big] = \exp\!\left(-\gamma\, w^\top\mu + \frac{\gamma^2}{2}\, w^\top\Sigma w\right)
$$

Minimizing that is the same as maximizing $w^\top\mu - \frac{\gamma}{2}\, w^\top\Sigma w$. So mean–variance is exactly the Gaussian special case; for a non-Gaussian generated distribution, exponential utility keeps the extra information. Constraints, risk aversion and cost treatment are your choices.

### Cluster tickers by their return distributions

For each ticker, collect its generated returns across samples and horizons. That is the ticker's multi-horizon return distribution given today's features, and it can serve as the ticker's representation: raw samples, quantiles, moments, density values, horizon profiles, or a learned embedding. Define a distance between two tickers' distributions (Wasserstein, energy distance, MMD, or a quantile-function distance), build the ticker-by-ticker distance matrix, and cluster it, for example hierarchically or spectrally.

{{< callout type="feat" >}}
**The concept.** Tickers are grouped by the information in their output return distributions. The representation and the clustering method are yours to experiment with.
{{< /callout >}}

There are two different readings, and they must stay separate:

| Compare | Groups tickers that… |
|---|---|
| marginal distribution shapes | have similar location, scale, tails or multi-horizon profile |
| movement within the same joint draws | co-move in the same scenarios (paired-sample correlation, or the root-mean-square gap between two tickers across paired draws) |

Two tickers can have identical marginal distributions while moving in opposite directions in every paired draw: same shape, negative co-movement. Because the distributions are conditional, the groups change as the features change, which gives a way to study shifting ticker relationships. Note that it is the tickers being clustered, not the generated scenarios.

### Everything else is yours

Once you have the joint distribution there is a lot to extract. Mean and median describe direction. Standard deviation gives output-return volatility at a chosen horizon, and mean over volatility gives a horizon-specific risk-adjusted statistic. Skewness, quantiles, tails and relationships across horizons add more. These are examples, not a prescribed strategy. The book gives you the features, the model, the generated distribution and ways to evaluate it; the rest of the research space is deliberately open. Questions to guide it:

1. What property of the return distribution do I want to study?
2. Which tickers and horizons should be considered together?
3. Which statistic, transformation, grouping or objective represents that property?
4. How will I compare its usefulness with realized outcomes?
