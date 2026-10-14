---
title: "Diffusion Implementation Book, Part 3: Model Structure"
description: "The denoiser as an encoder–decoder Transformer over tickers: how encoder and decoder tokens are built, the gated cross-attention that makes classifier-free guidance possible, how it differs from the original Transformer, and a shape walk-through."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 3
weight: 13
tags: ["diffusion", "generative-model", "transformer", "implementation", "quant-research"]
math: true
---

In its simplest form, each ticker's token is just a sum of four pieces, all
projected to the same hidden width $d$:

$$
\text{token}_i = P_X(\text{features}_i) + \text{ticker embedding}_i + P_Z(\text{noisy labels}_i) + \text{diffusion-step embedding}
$$

$P_X$ and $P_Z$ are linear projections. A Transformer encoder processes the
$N$ tokens jointly and a head outputs $\hat\varepsilon$ of shape `[N, H]`.

For the guidance methods in [Part 6](/quant/model/diffusion-book-part-6-guidance/)
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
