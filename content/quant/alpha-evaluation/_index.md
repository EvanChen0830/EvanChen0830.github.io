---
title: "Alpha Evaluation"
description: "Judging alpha quality after it leaves the model: leakage control, the full/inherited/novel decomposition, significance and breadth, whether a model's complexity earned anything over a linear surrogate, and the five diagnostic cases that decide what to deploy."
---

{{< pipeline highlight="alpha-evaluation" >}}

Alpha signals output by machine learning or generative models must face the same rigorous scrutiny as raw features. This stage evaluates whether model forecasts represent statistically significant trading power, diagnoses directional breadth across assets, and tests whether the model's non-linear complexity bought genuinely novel alpha or merely repackaged existing feature exposures.

## Guides in this section

- **[Label-Alpha Geometry: Significance, Mechanism Inheritance, and Residual Novelty](/quant/alpha-evaluation/label-alpha-geometry/)**  
  A complete end-to-end evaluation reference: point-in-time leakage control, the full/inherited/novel signal decomposition, walk-forward linear surrogate benchmarks, residual alpha testing, and the five diagnostic cases that decide production deployment.
