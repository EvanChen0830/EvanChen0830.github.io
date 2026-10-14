---
title: "Diffusion Implementation Book, Part 5: Training"
description: "The forward noising process, why 100 steps works if the beta schedule reaches pure noise, training joint conditional and unconditional branches, the reconstructed-clean-label IC loss, and SNR-based step weighting."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 5
weight: 15
tags: ["diffusion", "generative-model", "training", "ic-loss", "implementation", "quant-research"]
math: true
---

The forward process adds noise to the normalized labels. For a diffusion step $\tau$:

$$
Z_\tau = \sqrt{\bar\alpha_\tau}\, Z_0 + \sqrt{1 - \bar\alpha_\tau}\,\varepsilon, \qquad \varepsilon \sim \mathcal{N}(0, I)
$$

$$
\alpha_\tau = 1 - \beta_\tau, \qquad \bar\alpha_\tau = \alpha_1 \cdot \alpha_2 \cdots \alpha_\tau
$$

$$
\text{loss} = \text{average of } \|\varepsilon - \hat\varepsilon\|^2 \text{ over dates, diffusion steps and noise draws}
$$

The model usually predicts the noise. In each training pass, sample one diffusion step per example and build that noisy state directly from the formula; the forward noising steps are never run one by one. This is the standard DDPM training setup.

{{< svgfig src="diffusion-book/fig05.svg" >}}
**Figure 5.** One training step. Both branches see identical noisy labels, diffusion step and noise, so their losses can be compared directly. The IC loss uses the conditional branch.
{{< /svgfig >}}

### The forward process can have 100 steps

A common setup uses 1,000 noise levels, but about 100 usually works. The important part is the beta schedule: when you change the number of steps, the total noise must still be right. At the last step the labels must be effectively pure noise, meaning $\bar\alpha$ at the final step is close to 0 (equivalently, the final signal-to-noise ratio is close to 0), so that training ends where sampling begins.

Two shortcuts break this. Truncating a 1,000-step schedule to its first 100 values leaves most of the signal in. Keeping the usual linear endpoints ($\beta$ from $10^{-4}$ to $0.02$) over only 100 steps still leaves $\bar\alpha \approx 0.36$. So going to 100 steps means redesigning or rescaling the schedule until its end point matches the Gaussian noise used at inference. The paper *Common Diffusion Noise Schedules and Sample Steps are Flawed* discusses this end-point problem.

{{< svgfig src="diffusion-book/fig06.svg" >}}
**Figure 6.** Signal retained over the forward process. The dashed schedule never reaches pure noise. The cosine schedule does, and spends more of the trajectory at informative noise levels.
{{< /svgfig >}}

### The scheduler changes how noise is added

Linear and cosine are the familiar choices; other schedules from diffusion research can be used or customized. A linear schedule raises $\beta$ step by step in equal increments, but because $\bar\alpha$ is a running product, the signal does not fall in a straight line. A cosine schedule instead defines the signal curve $\bar\alpha$ directly and derives $\beta$ from it, spreading noise levels so more of the trajectory keeps informative signal (Improved DDPM).

The reason this matters: on nearly pure-noise inputs, predicting the noise is easy and teaches the model little about returns. A low average loss can therefore hide a weak model. Look at the schedule and at the loss per noise band. Still, the scheduler is just one parameter to tune; the guidance methods and label-space objectives that follow are the more substantive trading tools.

### Train a conditional and an unconditional branch

The same network is trained with the gate open (conditional) and closed (unconditional). Either choose a branch at random for each training example, or run both branches on every example and optimize both losses. The financial features are never altered: we switch whether the decoder uses them, rather than deleting individual features to fake a missing condition. "Unconditional" means without market features only; ticker identity and diffusion step stay in.

Running both branches on the same clean labels, diffusion step and noise gives two directly comparable losses:

$$
L_{\text{cond}} = \|\varepsilon - \hat\varepsilon_{\text{cond}}\|^2, \qquad L_{\text{uncond}} = \|\varepsilon - \hat\varepsilon_{\text{uncond}}\|^2
$$

$$
L_{\text{base}} = \lambda_{\text{cond}} \cdot L_{\text{cond}} + \lambda_{\text{uncond}} \cdot L_{\text{uncond}}
$$

Log both separately, including by diffusion-step band. The gap $L_{\text{uncond}} - L_{\text{cond}}$ shows how much the features help denoising, and pairing the noise and step makes that comparison much cleaner than comparing unrelated passes. It is a model diagnostic only. Whether the features are useful is decided by evaluating generated forecasts against future labels.

### IC loss on reconstructed clean labels

The model predicts noise, but from the noisy labels and the predicted noise you can rebuild an estimate of the clean labels: noisy labels minus the scaled noise, rescaled by the signal factor. Map that back through the inverse transform and you have return-space predictions a trading loss can use.

$$
\hat Z_0 = \frac{Z_\tau - \sqrt{1 - \bar\alpha_\tau}\,\hat\varepsilon}{\sqrt{\bar\alpha_\tau}}, \qquad \hat Y_0 = \text{inverse transform of } \hat Z_0
$$

This reconstruction starts from a noisy state built from a real label. It is not a one-step sampler from pure noise.

For stat-arb, score the reconstruction with a cross-sectional IC against the actual labels, computed within the same date and horizon across valid tickers:

$$
L_{\text{IC},h} = 1 - \operatorname{Corr}_{\text{tickers}}(\hat Y_{0,h}, Y_h)
$$

$$
L_{\text{IC}} = \sum_h \omega_h \cdot L_{\text{IC},h} \quad (\text{horizon weights})
$$

$$
\text{total loss} = L_{\text{base}} + \lambda_{\text{IC}} \cdot q(\tau) \cdot L_{\text{IC}}
$$

- Pearson IC is differentiable. Rank IC needs a differentiable rank approximation.
- The diffusion loss stays the base objective, because IC alone does not pin down a distribution. The IC term adds a preference for cross-sectional usefulness.
- If labels are normalized per ticker, IC in normalized space can differ from the economic IC. Choose the loss space deliberately; computing it after the inverse transform keeps it in return units.

### Different diffusion steps get different weight

This is the key training detail. At a late, noisy step the reconstruction is unreliable: a noise-prediction error of size $e$ becomes a clean-label error of size $e / \sqrt{\text{SNR}_\tau}$. The same small noise error turns into a very large label error when the signal-to-noise ratio is low. So weight the IC loss by how informative the step is, for example:

$$
q(\tau) = \frac{\text{SNR}_\tau}{\text{SNR}_\tau + \rho}, \qquad \rho > 0
$$

Any other schedule works, as does skipping the IC loss below a chosen SNR threshold.

{{< svgfig src="diffusion-book/fig07.svg" >}}
**Figure 7.** The weight $q(\tau)$ on a 100-step cosine schedule. Low- and mid-noise steps drive the IC loss; near-pure-noise steps contribute little.
{{< /svgfig >}}

If the last step has exactly zero signal, the reconstruction formula divides by zero. Skip the IC loss there, or use a clean-data or velocity output (see the reference in [Part 10](/quant/model/diffusion-book-part-10-checklist-and-reference/)). Multiplying an undefined value by a zero weight does not fix it.

{{< callout >}}
**Two different ICs.** *Training reconstruction IC* uses a noisy state built from a real label. *Inference IC* uses a signal from samples generated from fresh noise. Only the second tells you whether training improved the actual forecast.
{{< /callout >}}
