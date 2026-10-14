---
title: "Diffusion Implementation Book, Part 4: Label Preprocessing"
description: "Bringing heavy-tailed return labels to a stable scale: compress tails with an invertible transform, standardize, generate in that space and transform back — fitted on training dates only, with every sample inverse-transformed before any statistic."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 4
weight: 14
tags: ["diffusion", "generative-model", "label-preprocessing", "implementation", "quant-research"]
math: true
---

Before training, bring labels to a stable scale: roughly zero mean and unit
variance. Trading labels need care, because long tails make ordinary
mean–standard-deviation normalization unstable. Apply a nonlinear transform
first to compress the tails, then normalize.

{{< svgfig src="diffusion-book/fig04.svg" >}}
**Figure 4.** Synthetic heavy-tailed returns before and after an asinh
transform with standardization.
{{< /svgfig >}}

$$
\text{signed log:}\quad g(y) = \operatorname{sign}(y)\,\log\!\left(1 + \frac{|y|}{c}\right)
$$

$$
\text{inverse hyperbolic sine:}\quad g(y) = \operatorname{asinh}\!\left(\frac{y}{c}\right), \qquad c > 0
$$

$$
\text{then:}\quad Z_0 = \frac{g(Y) - m}{s}
$$

$c$, $m$ and $s$ are fitted on the training labels, per ticker and horizon,
or shared within an explicitly defined group.

The goal is a better-behaved label space. Standardizing doesn't make a
distribution Gaussian, and diffusion doesn't need Gaussian clean returns. The
near-Gaussian end point belongs to the noise the diffusion process adds, not
to the data.

### Transform, generate in that space, transform back

Save the fitted transform and use the same one everywhere:

$$
Y \;\longrightarrow\; Z \;\longrightarrow\; \text{diffusion training and sampling} \;\longrightarrow\; \hat Z \;\longrightarrow\; \hat Y
$$

1. **Inverse-transform every sample** before computing the mean, volatility,
   covariance or a portfolio objective. With a nonlinear inverse, the average
   of the inverted samples is not the inverse of the averaged sample, so the
   mean in return units must come from the inverted samples.
2. **Fit on training dates only.** Long tails are not leakage. Leakage
   happens when future labels influence the transform parameters, quantiles,
   clipping bounds or normalization statistics. Fit on the chronological
   training set, freeze, and reuse for validation, test and inference.
3. **Prefer invertible transforms.** When the full distribution matters,
   clipping is a bad idea: it changes the distribution and cannot be undone.
4. **Densities need the Jacobian.** If you compute continuous densities in
   transformed space and report them in return units, add the log-Jacobian of
   the transform. Probability masses of matching bins need no correction.
