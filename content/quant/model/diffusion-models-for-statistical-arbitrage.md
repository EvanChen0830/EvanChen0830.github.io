---
title: "Diffusion Models for Statistical Arbitrage: Seeing the Full Return Distribution"
description: "An introduction to the generative-model line of research: instead of predicting one number, learn the entire conditional return distribution and sample from it. Covers the diffusion mechanics, how inference and evaluation change when the output is a distribution, and two live applications."
date: 2026-09-14
tags: ["diffusion", "generative-model", "model", "alpha", "quant-research"]
math: true
ShowToc: true
TocOpen: true
---

Every method described elsewhere on this site — RankIC transformers, linear
factor models, the cross-sectional geometry framework — ultimately produces
a point summary: a rank, an expected return, a quality score. This line of
work asks a different question: given a signal, can we actually **see** the
return probability, rather than collapse it to one number?

{{< pipeline highlight="model,alpha,alpha-evaluation" >}}

## The idea

A generative model learns how future returns are distributed *given* the
features available right now, rather than learning a single expected
value. Feed it a condition — a trading signal, an alpha, any decision-time
information — and instead of a point forecast it hands back a full
distribution: possibly multimodal, possibly skewed, possibly just a vanilla
Gaussian when the condition carries no information about shape.

For statistical arbitrage specifically, there are many tickers, each with
labels at several horizons. The model learns $p(Y \mid X)$ jointly across
all of them, where $X$ holds the features available at decision time and
$Y$ is an $N \times H$ block of returns for $N$ tickers and $H$ horizons —
not $N \times H$ independent point predictions, one joint object.

## The mechanics: diffusion in one page

Diffusion models learn this distribution by a forward process that
destroys structure and a backward process that learns to undo it.

**Forward process.** Starting from real data $x_0$, noise is added over a
schedule of steps until nothing recognizable remains:

$$
\sqrt{\alpha_t}\, x_0 + \sqrt{1-\alpha_t}\,\epsilon \;=\; x_t,
\qquad \epsilon \sim \mathcal N(0, I).
$$

In continuous time this is a stochastic differential equation,

$$
dx = f(x,t)\,dt + g(t)\,dw,
$$

and the corresponding **backward process** — the one that actually
generates samples — runs it in reverse using the *score function*
$\nabla_{x_t} \log p(x_t)$, the gradient of the log-density of the noised
data:

$$
dx = \left[f(x,t) - g(t)^2 \nabla_{x_t}\log p(x_t)\right]d\tilde t + g(t)\,d\tilde w.
$$

**What the network actually learns.** The model never sees the label
directly as a regression target. It sees a noised version $x_t$ and the
step $t$, and predicts the noise $\epsilon_\theta$ that was added:

$$
\min_\theta\; \mathbb E\left[\,\lVert \epsilon - \epsilon_\theta(x_t, t)\rVert^2\,\right].
$$

This predicted noise *is* an estimate of the score function — hence
"score-based generative model." The network is fitting the data
distribution indirectly, through its gradient, one noise level at a time.

## Inference: sampling to build a distribution

A single forward pass doesn't give you "the prediction." It gives you one
sample. The distribution only appears once you sample repeatedly from the
same condition:

$$
x_0^{(1)}, x_0^{(2)}, \dots, x_0^{(n)} \;\sim\; p_\theta(\,\cdot \mid \text{condition}),
$$

and the empirical collection $\{x_0^{(i)}\}$ *is* the predicted return
distribution for that signal. This is the same operation as prompting an
image model for "a picture of a 5" a thousand times: no single output is
"the" answer, but the collection of outputs traces out the model's
learned concept of what a 5 looks like. Here the concept is "what returns
look like given this condition."

## Evaluation: you never observe the true distribution

This is where the pipeline's evaluation discipline has to be rebuilt from
scratch, because the usual move — compare prediction to realized label —
doesn't directly apply to a distributional output. For any one signal you
observe exactly one realized label $y_n$, never the distribution it was
drawn from.

The fix is to build a **proxy distribution**: pool the realized labels from
other signals whose *predicted* distributions are close to signal $n$'s
predicted distribution, and compare that pooled empirical distribution
against $P_n$ directly,

$$
\min_{} \; d\!\left(\text{Proxy Distribution}_n,\; P_n\right).
$$

If the model is well calibrated, signals that it believes have similar
return distributions should, in fact, realize similarly-distributed
outcomes — that agreement is what the distance $d(\cdot,\cdot)$ checks.
This is the alpha-evaluation question in distributional form: not "was the
point prediction close," but "was the *shape* of uncertainty honest."

## Two live applications

**Market making.** A single ticker has labels at multiple horizons, and
those horizon-specific predicted distributions don't have to agree. When
the short-horizon prediction moves against the long-horizon trend, the
honest response is to widen spreads and lean on inventory reversion rather
than take a view. When short- and long-horizon predictions agree, that's a
reason to shift weight from pure market-making toward directional skew —
and when only the short horizon carries a sharp, confident shape, that can
justify a direct directional take sized to that horizon alone. The point
is that the *joint* shape across horizons — not any single-horizon point
forecast — determines which of these regimes you're in.

**Cross-equity co-movement.** The same joint-distribution idea extends
across tickers: modeling $[\text{label}_A, \text{label}_B]$ together rather
than independently exposes co-movement structure directly in the sample
distribution, without a separate correlation-estimation step bolted on
afterward.

## Where this fits the pipeline

The distributional output doesn't replace the [direction / activity /
tail](/overview/) decomposition used everywhere else on this site — it
subsumes it. A learned distribution's first moment is the direction signal,
its dispersion and co-movement structure is the activity signal, and its
tail mass is the tail signal: three separate diagnostics bolted onto a
point prediction elsewhere become three different slices of one object you
can actually look at here.

The full derivation — model structure, tensor shapes, label preprocessing,
training objective weighting across diffusion steps, classifier-free and
auto-guidance, fast samplers, and the complete evaluation and checklist — is
written up in the companion implementation guide:

**[Diffusion models for stat-arb: an implementation guide →](/quant/model/diffusion-statarb-implementation-guide/)**

