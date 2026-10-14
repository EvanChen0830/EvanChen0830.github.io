---
title: "Diffusion Implementation Book, Part 7: Sampling"
description: "Generating the empirical conditional distribution: sampling K joint scenarios per date, choosing sample size K and monitoring stability, faster ODE solvers (DPM-Solver, UniPC), consistency distillation, and encoder reuse."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 7
weight: 17
tags: ["diffusion", "generative-model", "sampling", "dpm-solver", "distillation", "implementation", "quant-research"]
math: true
---

For each date we want a distribution. Hold one set of features fixed and generate $K$ samples. Each sample contains every ticker and every horizon together, and the collection of $K$ samples forms the empirical conditional distribution.

{{< svgfig src="diffusion-book/fig09.svg" >}}
**Figure 9.** Sampling pipeline and output. Each of the $K$ slices is a full $N \times H$ scenario for the same date.
{{< /svgfig >}}

In practice: repeat the condition across $K$ independent noise starts and run the sampler, batch-wise, giving output of shape `[B, K, N, H]`. The features are the same across the $K$ draws for one date; only the starting noise differs. Inside each draw, tickers are generated together. Generating each ticker on its own would destroy the joint dependence we set out to model.

### You choose K

$K$ might be 1,000, 10,000 or 100,000. Compare these on the same date: if the estimated mean or standard deviation moves materially, you need more samples. The same applies to covariance, quantiles, or whatever statistic you plan to use.

Use nested samples, so the first 1,000 draws are part of the 10,000 and 100,000 sets; the comparison then isolates the effect of adding samples. Repeating with independent seeds shows the remaining Monte Carlo noise. For the mean, the standard error is about $\sigma / \sqrt{K}$:

| $K$ | Standard error of the mean ($\sigma = 20\text{ bp}$) |
|---|---|
| 1,000 | $\approx 0.63\text{ bp}$ |
| 10,000 | $\approx 0.20\text{ bp}$ |
| 100,000 | $\approx 0.063\text{ bp}$ |

The real question is whether the quantities you need have stabilized. A $K$ that's enough for the mean can still be too small for a rare tail or an exponential-utility objective. More samples estimate the *model's* distribution better; how accurate that distribution is gets checked separately against realized labels.

### Faster inference: solvers and distillation

**A different solver.** Keep the trained model and cross the reverse process in fewer evaluations: DDIM, Euler, Heun, DPM-Solver, DPM-Solver++, UniPC, or other compatible higher-order or adaptive methods. A model trained with 100 noise levels can often be sampled in 10 or 20 steps if the solver is accurate enough.

**Distillation.** A slow many-step teacher trains a fast student: a consistency model, a progressively distilled diffusion model, or another teacher–student method that matches the teacher's predictions, transitions, trajectories or outputs. Consistency distillation uses a teacher; direct consistency training is another way to get a consistency model and is not teacher–student distillation.

**Reuse.** Encode the features once per date, reuse that encoding across all draws, and batch the $K$ samples. Where the architecture allows, cache the cross-attention keys and values from the encoder; the noisy-label states change at every step and can't be cached.

| Parameter | Meaning |
|---|---|
| $T_{\text{diff}}$ | number of forward noise levels |
| $S$ | number of reverse solver steps |
| $K$ | number of generated samples per date |

Cost depends on network evaluations, not just steps: guidance evaluates two branches, and some solvers evaluate several times per step. Raising $K$ estimates the chosen sampler's distribution more precisely but cannot fix an inaccurate solver. Whatever you speed up, the fast output must still represent the intended multi-ticker, multi-horizon distribution: compare mean, variance, dependence, tails and the downstream numbers you actually use.
