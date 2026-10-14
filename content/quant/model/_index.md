---
title: "Model"
description: "Model architecture, fitting, and validation: from point predictions to generative diffusion models learning the joint conditional return distribution."
---

{{< pipeline highlight="model" >}}

The Model stage converts engineered and evaluated features into forecasts. Rather than collapsing predictions into single scalar points, the primary modeling line here uses **deep generative diffusion models** to learn the entire conditional joint distribution $p(Y \mid X)$ across multiple tickers and multiple horizons simultaneously.

## Highlights in this section

- **[Diffusion Models for Statistical Arbitrage: Seeing the Full Return Distribution](/quant/model/diffusion-models-for-statistical-arbitrage/)**  
  Conceptual overview of generative modeling for stat-arb: the diffusion mechanics, why point forecasts discard shape, proxy-distribution evaluation, and applications in market-making and cross-equity co-movement.

- **[Diffusion Models for Stat-Arb: An Implementation Book (10 Parts)](/quant/model/diffusion-book-part-1-the-idea/)**  
  A complete implementation reference covering tensor contracts, Transformer denoiser structure, label preprocessing, training losses, guidance, fast sampling, binned likelihood evaluation, and portfolio construction:
  1. [The Idea](/quant/model/diffusion-book-part-1-the-idea/)
  2. [Data and Shapes](/quant/model/diffusion-book-part-2-data-and-shapes/)
  3. [Model Structure](/quant/model/diffusion-book-part-3-model-structure/)
  4. [Label Preprocessing](/quant/model/diffusion-book-part-4-label-preprocessing/)
  5. [Training](/quant/model/diffusion-book-part-5-training/)
  6. [Guidance](/quant/model/diffusion-book-part-6-guidance/)
  7. [Sampling](/quant/model/diffusion-book-part-7-sampling/)
  8. [Evaluation](/quant/model/diffusion-book-part-8-evaluation/)
  9. [Using the Distribution](/quant/model/diffusion-book-part-9-using-the-distribution/)
  10. [Checklist and Reference](/quant/model/diffusion-book-part-10-checklist-and-reference/)
