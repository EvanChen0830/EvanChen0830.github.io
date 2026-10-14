---
title: "Quant Research"
description: "Methodology notes, in the order a signal actually moves through the pipeline: universe, feature engineering, feature evaluation, model, alpha, alpha evaluation, final ensemble."
---

Every write-up below sits at one stage of the same pipeline. A signal is only as good as the weakest stage it passed through, so each stage gets its own evaluation step rather than being judged only by the final result.

{{< pipeline >}}

## Pipeline Sections & Featured Guides

### 1. [Universe](/quant/universe/)
Tests whether an asset universe carries authentic geometric structure rather than hindsight or random portfolio variance.
- **[Multi-Sub-Universe Directional Coupling Geometry](/quant/universe/multi-sub-universe-directional-coupling-geometry/)** — The restriction Gram operator, canonical $\mathsf V \leftrightarrow \mathsf D$ coupling, and size-matched membership null tests.

### 2. [Feature Engineering](/quant/feature-engineering/)
Constructing raw signals from market tick data, order book dynamics, fundamental ratios, and cross-asset interactions.

### 3. [Feature Evaluation](/quant/feature-evaluation/)
Disentangling raw predictive signal into orthogonal Direction, Activity, and Tail channels, calibrated against permutation nulls.
- **[Cross-Sectional Feature Selection: An Implementation Guide](/quant/feature-evaluation/cross-sectional-feature-selection/)** — Matrix pipeline from whitening to greedy coverage fill and protecting directional exposure.
- **[Cross-Sectional Predictive Geometry: Formal Theory & Operator SVD](/quant/feature-evaluation/cross-sectional-predictive-geometry/)** — Rigorous mathematical propositions, within-date centering, operator SVD, and multiple testing control.
- **[Eligibility, Grouping, and Submodular Coverage: The Selection Bridge](/quant/feature-evaluation/selection-bridge-eligibility-grouping-coverage/)** — Five eligibility gates, mechanism grouping, and submodular coverage optimization.

### 4. [Model](/quant/model/)
Generative modeling of multi-horizon, cross-sectional return distributions rather than point predictions.
- **[Diffusion Models for Statistical Arbitrage: An Implementation Guide](/quant/model/diffusion-statarb-implementation-guide/)** — End-to-end reference covering Transformer denoisers, tensor shapes, guidance, fast sampling, binned likelihoods, and portfolio construction.

### 5. [Alpha](/quant/alpha/)
Translating raw conditional distributions and model forecasts into tradeable, risk-adjusted alpha signals.

### 6. [Alpha Evaluation](/quant/alpha-evaluation/)
Evaluating model predictions with the same scrutiny as raw features.
- **[Label-Alpha Geometry: Significance, Mechanism Inheritance, and Residual Novelty](/quant/alpha-evaluation/label-alpha-geometry/)** — Walk-forward linear surrogates, residual alpha decomposition, and the five diagnostic cases deciding deployment.

### 7. [Final Ensemble](/quant/final-ensemble/)
Combining alphas across models, horizons, and sub-universes into executed positions with transaction-cost discipline.
