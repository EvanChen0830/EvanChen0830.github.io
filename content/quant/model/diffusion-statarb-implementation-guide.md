---
title: "Diffusion Models for Statistical Arbitrage: An Implementation Guide"
description: "An end-to-end implementation guide to learning the joint distribution of multi-ticker, multi-horizon returns: model structure, tensor shapes, training, guidance, sampling, evaluation, and portfolio construction."
date: 2026-10-02
tags: ["diffusion", "generative-model", "implementation", "quant-research", "stat-arb"]
math: true
aliases:
  - "/diffusion-statarb-implementation-book.html"
  - "/quant/model/diffusion-book-part-1-the-idea/"
  - "/quant/model/diffusion-book-part-2-data-and-shapes/"
  - "/quant/model/diffusion-book-part-3-model-structure/"
  - "/quant/model/diffusion-book-part-4-label-preprocessing/"
  - "/quant/model/diffusion-book-part-5-training/"
  - "/quant/model/diffusion-book-part-6-guidance/"
  - "/quant/model/diffusion-book-part-7-sampling/"
  - "/quant/model/diffusion-book-part-8-evaluation/"
  - "/quant/model/diffusion-book-part-9-using-the-distribution/"
  - "/quant/model/diffusion-book-part-10-checklist-and-reference/"
---

This book is an implementation guide to learning the joint distribution of
multi-ticker, multi-horizon returns: model structure, tensor shapes,
training, guidance, sampling, evaluation, and how to use the samples. For a
shorter conceptual introduction, start with
[Diffusion Models for Statistical Arbitrage](/quant/model/diffusion-models-for-statistical-arbitrage/).

{{< figkey >}}

## 1. The Idea

Start from the return distribution itself. A deep generative model can learn
how future returns are distributed given the features available now.
Diffusion, flow matching and consistency models are different ways to build
such a model. The question this book studies is what becomes available to a
trader once the model produces a whole conditional distribution rather than a
single forecast.

For statistical arbitrage there are many tickers, and each ticker has labels
at several horizons. We model all of them together: the model learns
$p(Y \mid X)$, where $X$ holds the features available at market time $t$ and
$Y$ is an $N \times H$ block of returns for $N$ tickers and $H$ horizons.

Two structures come out of that joint view.

**Across tickers.** The joint distribution shows how tickers move together,
how their returns differ, and how they form groups. It can contain common
movement, relative movement, and more than one peak. For stat-arb, these
relationships are part of what we want to learn, not a nuisance.

**Across horizons.** A ticker's label is a vector of returns at several
horizons. Its distribution says how the future unfolds: whether a direction
persists or reverses, and whether the features point to momentum-like or
mean-reverting behaviour.

{{< callout type="feat" >}}
**The concept in one line.** Features condition a joint distribution of
multi-ticker, multi-horizon returns. The distribution is the output we study.
{{< /callout >}}

### Cumulative or incremental labels?

If horizon labels are cumulative returns from the same start time, their
intervals overlap, so they are positively correlated even when successive
increments are independent. To read persistence or reversal, know which kind
of label you have, and compare the conditional horizon structure with the
unconditional one.

$$
\text{cumulative log return:}\quad R_{t,h} = \log P(t + \Delta_h) - \log P(t)
$$

$$
\text{non-overlapping increment:}\quad I_{t,h} = R_{t,h} - R_{t,h-1}, \qquad R_{t,0} = 0
$$

Both are useful; pick the one that matches the structure you want to
investigate.

The worked scenario is statistical arbitrage. The same idea works at high
frequency, where the feature set, label timing and compute budget change with
the trading setup.

---

## 2. Data and Shapes

Before any training technique, be exact about what goes in and what comes
out. One training example is one date.

{{< svgfig src="diffusion-book/fig01.svg" >}}
**Figure 1.** What goes in and what comes out. Each ticker row becomes one
token; its $H$ label values are coordinates of that token, not separate
tokens.
{{< /svgfig >}}

For each ticker the model uses four things: that ticker's features, a ticker
embedding, its noisy multi-horizon label vector, and an embedding of the
diffusion step. It outputs one multi-horizon vector per ticker. With noise
(epsilon) prediction, that vector is the predicted noise on the ticker's
labels. After the full reverse process, what comes out is the ticker's
generated multi-horizon return vector.

| Tensor | Shape | Notes |
|---|---|---|
| features $X$ | `[B, N, F]` | or `[B, N, L, F]` with $L$ past timestamps |
| labels $Y$ | `[B, N, H]` | multi-horizon returns, same ticker order as $X$ |
| normalized clean labels $Z_0$ | `[B, N, H]` | labels after the transform in [Chapter 4: Label Preprocessing](#4-label-preprocessing) |
| noisy labels $Z_\tau$ | `[B, N, H]` | $Z_0$ with noise added at step $\tau$ |
| ticker ids | `[B, N]` | index into the ticker embedding table |
| valid mask | `[B, N]` | marks tickers not in the universe on that date |
| diffusion step $\tau$ | `[B]` | one noise level per example |
| model output $\hat\varepsilon$ | `[B, N, H]` | predicted noise, same shape as the labels |
| generated samples | `[B, K, N, H]` | $K$ joint scenarios per date ([Chapter 7: Sampling](#7-sampling)) |

**The ticker is the Transformer token. The horizons are coordinates inside
each ticker's label vector.** Making every ticker–horizon pair its own token
isn't needed for the structures here. With 200 tickers and four horizons, the
sequence has 200 tokens, each with a four-dimensional noisy-label input and a
four-dimensional noise-prediction output.

### Four kinds of time

These get mixed up easily. The time embedding inside the denoiser is
diffusion time. It is never a forecast horizon.

| Quantity | What it means |
|---|---|
| market time $t$ | when the conditioning observation is available |
| feature-history index | which past observation a feature comes from (the $L$ axis) |
| label horizon $\Delta_h$ | how far into the future the return is measured (the $H$ axis) |
| diffusion time $\tau$ | how much noise has been added to the labels |

### One timestamp or a history of features

Features can be a single snapshot, shape `[N, F]`, or a history, shape
`[N, L, F]`. With a history, either encode each ticker's history into one
vector first and pass that into the joint model, or keep both the ticker axis
and the history axis inside the model. This is an architecture choice only.
The target doesn't change: each ticker still has its own multi-horizon label
vector, and all tickers are generated jointly.

{{< callout >}}
**Universe changes.** Features, ticker embeddings, labels and outputs must
all refer to the same ticker identities. When the universe changes, carry ids
and masks instead of assuming a row position always means the same ticker.
{{< /callout >}}

---

## 3. Model Structure

In its simplest form, each ticker's token is just a sum of four pieces, all
projected to the same hidden width $d$:

$$
\text{token}_i = P_X(\text{features}_i) + \text{ticker embedding}_i + P_Z(\text{noisy labels}_i) + \text{diffusion-step embedding}
$$

$P_X$ and $P_Z$ are linear projections. A Transformer encoder processes the
$N$ tokens jointly and a head outputs $\hat\varepsilon$ of shape `[N, H]`.

For the guidance methods in [Chapter 6: Guidance](#6-guidance)
we separate the two sides instead: an encoder reads the features, and a
decoder denoises the labels. This is the encoder–decoder Transformer from
*Attention Is All You Need*, with three changes: the inputs are tickers
instead of words, the ticker embedding plays the role of positional encoding,
and the decoder's self-attention has no causal mask.

{{< svgfig src="diffusion-book/fig02.svg" >}}
**Figure 2.** Token construction for the encoder–decoder version. Features
and noisy labels are projected to the same width $d$ and summed with a shared
ticker embedding. The decoder token also gets a diffusion-step embedding.
{{< /svgfig >}}

{{< svgfig src="diffusion-book/fig03.svg" >}}
**Figure 3.** The denoiser, drawn in the layout of the original Transformer.
Encoder on the left, decoder on the right, residual connections around every
sublayer. The encoder output $C$ enters the decoder only through
cross-attention, and a gate $G$ switches that path on or off. The output is a
linear regression head, with no softmax.
{{< /svgfig >}}

### How it differs from the original Transformer

| Original Transformer | This model |
|---|---|
| tokens are words in a sentence | tokens are tickers on one date |
| positional encoding | ticker embedding, shared by encoder and decoder |
| decoder input: previous output tokens | noisy labels $Z_\tau$ plus a diffusion-step embedding |
| masked (causal) decoder self-attention | full self-attention across tickers, no mask |
| cross-attention always on | cross-attention gated by $G \in \{0, 1\}$ |
| linear + softmax over a vocabulary | linear $d \to H$, predicted noise per horizon |
| generates one token at a time | denoises all tickers at once, over several diffusion steps |

"Decoder" here only names the denoising half. Nothing is generated left to
right.

### Shape walk-through

| Stage | Output shape |
|---|---|
| feature embedding + ticker embedding | `[B, N, d]` |
| encoder stack | `C: [B, N, d]` |
| label embedding + ticker + diffusion-step embedding | `[B, N, d]` |
| decoder block (self-attn, gated cross-attn, feed forward) | `[B, N, d]` |
| output head | `ε̂: [B, N, H]` |

Setting $G = 0$ gives a decoder that works on the noisy labels, ticker
embedding and diffusion step alone. The market features stay in the encoder;
the decoder simply stops reading them. This is the clean on/off switch that
classifier-free guidance needs.

{{< callout >}}
**Gate every path.** If market features reach the decoder any other way,
such as a shortcut or a conditional normalization layer, that path must
switch off with $G$ too. Otherwise the "unconditional" branch is not
unconditional.
{{< /callout >}}

---

## 4. Label Preprocessing

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

---

## 5. Training

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

If the last step has exactly zero signal, the reconstruction formula divides by zero. Skip the IC loss there, or use a clean-data or velocity output (see the reference in [Chapter 10: Checklist and Reference](#10-checklist-and-reference)). Multiplying an undefined value by a zero weight does not fix it.

{{< callout >}}
**Two different ICs.** *Training reconstruction IC* uses a noisy state built from a real label. *Inference IC* uses a signal from samples generated from fresh noise. Only the second tells you whether training improved the actual forecast.
{{< /callout >}}

---

## 6. Guidance

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

The branches come from the gate in Figure 3 ([Chapter 3: Model Structure](#3-model-structure)): cross-attention on for the conditional branch, off for the unconditional one, with every other feature path switched the same way. The features stay intact; only their use is switched, so no feature vectors are built by arbitrarily deleting components.

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

**Or guide toward an objective directly.** Take an alpha score available at inference and define a reward: the horizon-weighted cross-sectional correlation between that score and the reconstructed returns. At each sampling step, add the gradient of that reward to the model's score (converted to a noise adjustment if the sampler works with predicted noise). Generated outcomes then lean toward the alpha's ordering. That agreement is built in by design, so only held-out labels can tell you whether realized IC improved. The training-time IC loss from [Chapter 5: Training](#5-training) targets the observed labels directly.

---

## 7. Sampling

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

---

## 8. Evaluation

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

---

## 9. Using the Distribution

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

## 10. Checklist and Reference

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
