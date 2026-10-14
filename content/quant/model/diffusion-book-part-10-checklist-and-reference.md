---
title: "Diffusion Implementation Book, Part 10: Checklist and Reference"
description: "The 11-step implementation checklist, alternative output parameterizations (noise, velocity, clean data, score, flow matching, consistency models), and literature references."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 10
weight: 20
tags: ["diffusion", "generative-model", "checklist", "velocity-prediction", "flow-matching", "implementation", "quant-research"]
math: true
---

### Checklist

1. **Tickers are tokens; horizons are coordinates.** IDs and masks travel with every batch.
2. **Label transform fitted on training dates only**, invertible, frozen across splits.
3. **Final noise step is actually pure noise** ($\bar\alpha \approx 0$).
4. **Market features reach the decoder only through gated paths** (cross-attention).
5. **Both branches trained** (conditional and unconditional); loss gap logged per noise band.
6. **IC loss in return units**, weighted by SNR, skipped where undefined.
7. **Every sample inverse-transformed** before calculating any summary statistic.
8. **Sample size $K$ chosen by nested-sample stability** of the specific statistic used.
9. **Fast solvers checked against the slow sampler** on the downstream numbers you trade.
10. **Evaluation on one fixed bin grid**, against unconditional and conditional Gaussians, with serial-dependence-aware error bars.
11. **Samples kept paired across tickers** for covariance, portfolios, and co-movement analysis.

### Other output types

Noise prediction is the primary example in this guide. The same token and conditioning structure works with other parameterizations, as long as the training target and sampler agree.

| Output | What the model predicts | Clean labels from it |
|---|---|---|
| noise ($\varepsilon$) | the added noise | $\hat Z_0 = (Z_\tau - \sqrt{1 - \bar\alpha}\,\hat\varepsilon) / \sqrt{\bar\alpha}$ |
| diffusion velocity ($v$) | $v = \sqrt{\bar\alpha}\,\varepsilon - \sqrt{1 - \bar\alpha}\,Z_0$ | $\hat Z_0 = \sqrt{\bar\alpha}\,Z_\tau - \sqrt{1 - \bar\alpha}\,\hat v$ (defined even at zero final signal) |
| clean data | $Z_0$ directly | the prediction itself |
| score | gradient of log density of noisy state $\nabla_{Z_\tau} \log p(Z_\tau)$ | via $\hat\varepsilon = -\sqrt{1 - \bar\alpha} \cdot \text{score}$ |
| flow matching | path velocity; on straight path from noise $\xi$ to $Z_0$, target is $Z_0 - \xi$ | integrate velocity ODE; not the same target as diffusion velocity |
| consistency model | consistent clean output from any point on generation trajectory | one or a few steps |

Across these families the common object is the conditional joint return distribution; only the forward trajectory construction and the training objective differ.

### References

These papers represent the foundational generative methods adapted in this work. The volatility-only guider, the cross-sectional adaptations, and the financial feature architecture are the trading-specific ideas developed throughout this guide.

- [Attention Is All You Need](https://arxiv.org/abs/1706.03762) (Vaswani et al., 2017)
- [Denoising Diffusion Probabilistic Models](https://arxiv.org/abs/2006.11239) (Ho et al., 2020)
- [Improved Denoising Diffusion Probabilistic Models](https://arxiv.org/abs/2102.09672) (Nichol & Dhariwal, 2021)
- [Classifier-Free Diffusion Guidance](https://arxiv.org/abs/2207.12598) (Ho & Salimans, 2022)
- [Guiding a Diffusion Model with a Bad Version of Itself](https://arxiv.org/abs/2406.02507) (Karras et al., NeurIPS 2024), with the [official AutoGuidance implementation](https://github.com/NVlabs/edm2)
- [Common Diffusion Noise Schedules and Sample Steps are Flawed](https://arxiv.org/abs/2305.08891) (Lin et al., 2023)
- [Denoising Diffusion Implicit Models](https://arxiv.org/abs/2010.02502) (Song et al., 2020)
- [DPM-Solver: A Fast ODE Solver for Diffusion Probabilistic Model Sampling](https://arxiv.org/abs/2206.00927) (Lu et al., 2022)
- [DPM-Solver++: Fast Solver for Guided Sampling of Diffusion Probabilistic Models](https://arxiv.org/abs/2211.01095) (Lu et al., 2022)
- [UniPC: A Unified Predictor-Corrector Framework for Fast Sampling of Diffusion Models](https://arxiv.org/abs/2302.04867) (Zhao et al., 2023)
- [Progressive Distillation for Fast Sampling of Diffusion Models](https://arxiv.org/abs/2202.00512) (Salimans & Ho, 2022)
- [Consistency Models](https://arxiv.org/abs/2303.01469) (Song et al., 2023)
- [Flow Matching for Generative Modeling](https://arxiv.org/abs/2210.02747) (Lipman et al., 2022)
