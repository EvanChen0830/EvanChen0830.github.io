---
title: "Model"
description: "Model architecture, fitting, and validation: from point predictions to generative diffusion models learning the joint conditional return distribution."
---

{{< pipeline highlight="model" >}}

The Model stage converts engineered and evaluated features into forecasts. Rather than collapsing predictions into single scalar points, the primary modeling line here uses **deep generative diffusion models** to learn the entire conditional joint distribution $p(Y \mid X)$ across multiple tickers and multiple horizons simultaneously.

## Highlights in this section

- **[Diffusion Models for Statistical Arbitrage: Seeing the Full Return Distribution](/quant/model/diffusion-models-for-statistical-arbitrage/)**  
  Conceptual overview of generative modeling for stat-arb: the diffusion mechanics, why point forecasts discard shape, proxy-distribution evaluation, and applications in market-making and cross-equity co-movement.

- **[Diffusion Models for Statistical Arbitrage: An Implementation Guide](/quant/model/diffusion-statarb-implementation-guide/)**  
  A complete, end-to-end implementation reference covering tensor contracts, Transformer denoiser structure, label preprocessing, training losses, guidance, fast sampling, binned likelihood evaluation, and portfolio construction.
