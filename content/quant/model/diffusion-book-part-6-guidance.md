---
title: "Diffusion Implementation Book, Part 6: Guidance"
description: "Amplifying signal in noisy markets: classifier-free guidance, AutoGuidance using a weaker model as reference, volatility-only guiders to amplify direction, and reward-guided inference."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 6
weight: 16
tags: ["diffusion", "generative-model", "guidance", "cfg", "autoguidance", "implementation", "quant-research"]
math: true
---

Trading labels are noisy. A feature can carry real information while its effect stays small compared with the spread of future returns. Guidance runs two versions of the model at each sampling step and pushes the result further in the direction where they disagree. Which reference branch you choose decides what gets amplified.

{{< svgfig src="diffusion-book/fig08.svg" >}}
**Figure 8.** Three ways to build the reference branch. The guided prediction extrapolates away from the reference.
{{< /svgfig >}}

### Classifier-free guidance: amplify the condition

$$
\varepsilon_{\text{CFG}} = \varepsilon_{\text{uncond}} + \eta \cdot (\varepsilon_{\text{cond}} - \varepsilon_{\text{uncond}})
$$

$\eta = 1$ gives the plain conditional prediction. $\eta > 1$ amplifies the difference the condition makes.

**Classifier-free guidance amplifies the condition in the output distribution.** In trading experiments with noisy labels, this amplification has helped. The guidance scale changes the sampling distribution, so the guided output is evaluated with the same distribution and trading metrics as any other model.

The branches come from the gate in Figure 3 ([Part 3](/quant/model/diffusion-book-part-3-model-structure/)): cross-attention on for the conditional branch, off for the unconditional one, with every other feature path switched the same way. The features stay intact; only their use is switched, so no feature vectors are built by arbitrarily deleting components.

### AutoGuidance: use a bad guider

AutoGuidance comes from *Guiding a Diffusion Model with a Bad Version of Itself* (NeurIPS 2024, oral). It extends the CFG idea: instead of an unconditional model, the reference is a deliberately weaker version of the model under the same condition.

$$
\varepsilon_{\text{AG}} = \varepsilon_{\text{good}} + \eta \cdot (\varepsilon_{\text{good}} - \varepsilon_{\text{bad}})
$$

$\eta = 0$ gives the good model with no extra guidance.

The bad model can be smaller, trained for fewer iterations, or taken from an earlier checkpoint. Both models receive the same features and model the same label distribution, with compatible noise levels and the same type of output. Guidance then amplifies what the good model learned beyond the weaker one. The bad guider must still be a meaningful model of the same problem; an unrelated or collapsed model gives a meaningless contrast.

| Method | What the difference between branches represents |
|---|---|
| Classifier-free guidance | conditioned prediction vs. prediction without market features |
| AutoGuidance | better approximation vs. weaker approximation, same condition |
| Selective-conditioning guider | full information vs. information with one chosen source removed |

The AutoGuidance paper establishes the method for generative modelling in general. Using it on joint financial returns is a trading research choice.

### Choose the guider by what you want to amplify

The reference doesn't have to be a generic weaker model. You can choose what information it contains, and that choice defines what the difference between branches means.

**A volatility-only guider amplifies direction.** Split the features into volatility features and directional features. The full model sees both; the reference sees only volatility features.

$$
\varepsilon_{\text{guided}} = \varepsilon_{\text{full}} + \eta \cdot (\varepsilon_{\text{full}} - \varepsilon_{\text{vol-only}})
$$

The volatility-only model describes the distribution given volatility; the full model adds the directional features. Their difference, in score terms, is the extra evidence the directional features add beyond what volatility explains, so amplifying it amplifies directional information. This can be two models, or one architecture with separately gated volatility and directional pathways; either way the volatility-only branch must be trained for its task.

This is a selective-conditioning extension inspired by guidance, and it differs from the original AutoGuidance: the reference has less information, rather than being a worse approximation of the same thing. Its meaning depends on the features. If the directional features also predict scale or tails, those effects are amplified too, and keeping volatility in both branches does not fix the variance of the guided distribution.

**The same idea can target cross-sectional information.** If the downstream metric is cross-sectional IC, you may want to amplify relative ordering. Let the full branch use the cross-sectional alpha features and the reference keep volatility and market-wide information but drop those alpha features. This strengthens the ranking information without requiring the whole market to shift. A long–short book can profit from A beating B even when both are expected to fall.

**Or guide toward an objective directly.** Take an alpha score available at inference and define a reward: the horizon-weighted cross-sectional correlation between that score and the reconstructed returns. At each sampling step, add the gradient of that reward to the model's score (converted to a noise adjustment if the sampler works with predicted noise). Generated outcomes then lean toward the alpha's ordering. That agreement is built in by design, so only held-out labels can tell you whether realized IC improved. The training-time IC loss from [Part 5](/quant/model/diffusion-book-part-5-training/) targets the observed labels directly.
