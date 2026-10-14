---
title: "Diffusion Implementation Book, Part 1: The Idea"
description: "An implementation guide to learning the joint distribution of multi-ticker, multi-horizon returns. Part 1: why model the whole conditional distribution, the two structures it exposes (across tickers and across horizons), and cumulative vs. incremental labels."
date: 2026-10-02
series: ["Diffusion Implementation Book"]
part: 1
weight: 11
tags: ["diffusion", "generative-model", "implementation", "quant-research"]
math: true
aliases: ["/diffusion-statarb-implementation-book.html"]
---

This book is an implementation guide to learning the joint distribution of
multi-ticker, multi-horizon returns: model structure, tensor shapes,
training, guidance, sampling, evaluation, and how to use the samples. For a
shorter conceptual introduction, start with
[Diffusion Models for Statistical Arbitrage](/quant/model/diffusion-models-for-statistical-arbitrage/).

{{< figkey >}}

## The idea

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
