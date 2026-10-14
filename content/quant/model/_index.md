---
title: "Model"
description: "Model architecture, fitting, and validation: from point predictions to generative diffusion models learning the joint conditional return distribution."
---

{{< pipeline highlight="model" >}}

The Model stage converts engineered and evaluated features into forecasts. Rather than collapsing predictions into single scalar points, the primary modeling line here uses **deep generative diffusion models** to learn the entire conditional joint distribution $p(Y \mid X)$ across multiple tickers and multiple horizons simultaneously.

## Highlights in this section

- **[Diffusion Models for Statistical Arbitrage: An Implementation Guide](/quant/model/diffusion-statarb-implementation-guide/)**  
  A complete, end-to-end reference covering the conceptual foundations of conditional return distributions, tensor contracts, Transformer denoiser architecture, label preprocessing, training losses, guidance, fast sampling, binned likelihood evaluation, and portfolio construction with downstream market making and cross-equity applications.
