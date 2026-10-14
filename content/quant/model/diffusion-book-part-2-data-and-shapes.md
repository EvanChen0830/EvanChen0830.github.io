---
title: "Diffusion Implementation Book, Part 2: Data and Shapes"
description: "Exactly what goes in and what comes out: one training example is one date, the ticker is the Transformer token, horizons are coordinates inside each token, the four kinds of time that get confused, and how to handle feature history and universe changes."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 2
weight: 12
tags: ["diffusion", "generative-model", "implementation", "quant-research"]
math: true
---

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
| normalized clean labels $Z_0$ | `[B, N, H]` | labels after the transform in [Part 4](/quant/model/diffusion-book-part-4-label-preprocessing/) |
| noisy labels $Z_\tau$ | `[B, N, H]` | $Z_0$ with noise added at step $\tau$ |
| ticker ids | `[B, N]` | index into the ticker embedding table |
| valid mask | `[B, N]` | marks tickers not in the universe on that date |
| diffusion step $\tau$ | `[B]` | one noise level per example |
| model output $\hat\varepsilon$ | `[B, N, H]` | predicted noise, same shape as the labels |
| generated samples | `[B, K, N, H]` | $K$ joint scenarios per date ([Part 7](/quant/model/diffusion-book-part-7-sampling/)) |

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
