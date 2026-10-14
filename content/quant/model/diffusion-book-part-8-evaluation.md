---
title: "Diffusion Implementation Book, Part 8: Evaluation"
description: "Evaluating distributional forecasts when only one label is realized: fixed unconditional bins, binned log-likelihood scoring, comparing against conditional and unconditional Gaussians (ΔLS), and evaluating joint structure."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 8
weight: 18
tags: ["diffusion", "generative-model", "evaluation", "log-likelihood", "calibration", "implementation", "quant-research"]
math: true
---

Evaluation is harder because each date has exactly one realized label vector. We never see the true conditional distribution for that exact set of features. So we score each realized label under its predicted distribution and aggregate across the test set, comparing likelihoods against a Gaussian.

{{< svgfig src="diffusion-book/fig10.svg" >}}
**Figure 10.** Scoring one realized label. Both models are scored on the same fixed bins. A two-peaked generated distribution can beat a Gaussian with the same mean and variance.
{{< /svgfig >}}

### Cut bins from the unconditional distribution

Define one fixed bin grid from the unconditional label distribution of the training data, then use the same edges for every conditional output and every baseline. Never re-cut bins per date. This matters for a two-peaked distribution: on a common grid of about 100 bins, the model can put probability in two separate regions, and those regions mean the same thing across dates. Re-cutting changes the coordinate system and makes scores incomparable.

Equal-width or unconditional-quantile bins both work. Include lower and upper overflow bins so future labels outside the training range still land somewhere. "Full dataset distribution" means the full training-label distribution for that research split; using test or future labels to set the grid leaks information.

### Likelihood of each realized label

$$
\text{generated mass in bin } b = \frac{\text{samples in bin } b + \kappa}{K + M \cdot \kappa}
$$

$$
\text{score for a date} = \log\!\big(\text{generated mass of the bin containing the realized label}\big)
$$

$M$ is the number of bins and $\kappa$ a small smoothing pseudo-count, so empty bins stay usable.

This is a binned likelihood score: the log probability the model gave the bin where the outcome landed. For many tickers and horizons, average the scalar scores within each date (with your horizon weights), then over dates. That keeps a date with a bigger universe from counting more just because it has more tickers.

### Compare with a Gaussian

Fit a Gaussian baseline on training data, either unconditional or conditional (mean and variance depending on features). Its probability for a bin is the normal CDF at the upper edge minus the CDF at the lower edge. Score the realized label in the same bin, and average the difference:

$$
\Delta\text{LS} = \operatorname{average}_{\text{test dates}}\big[\log(\text{generated mass}) - \log(\text{Gaussian mass})\big]
$$

*Example:* masses of $0.03$ vs $0.02$ in the realized bin contribute $\log(1.5) \approx 0.405\text{ nats}$ of improvement.

A positive value means the generated distribution, on average, put more probability where the outcomes actually landed. Both models must be scored the same way: don't subtract a Gaussian point density from a histogram mass; for this method both are bin masses. The conditional Gaussian is the more informative comparison, because it already tracks changing means and volatilities, so what's left over reflects the flexible shape the generative model learned.

### Keep the meaning of the score clear

- **It scores marginals.** Summing ticker–horizon scores evaluates marginal distributions, not the joint. A full histogram with 100 bins per coordinate would need $100^{N \cdot H}$ cells, which is impossible. Add likelihoods of selected spreads or portfolios, and sample-based joint scores.
- **Complementary tools.** CRPS for individual distributions; energy and variogram scores for aspects of the joint; quantile coverage, PIT diagnostics and regime-by-regime comparisons for calibration. The fixed-bin Gaussian comparison stays the main method.
- **Unequal-width bins.** If you want a continuous density estimate, divide a finite bin's mass by its width. Infinite overflow bins are scored as masses, not as a density over infinite width.
- **Uncertainty.** Keep the per-date model-minus-baseline differences for error bars. Overlapping horizons and serial market dependence make test dates dependent. The $K$ generated draws are samples from the model, not $K$ independent market outcomes.
