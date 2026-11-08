---
title: "Quant Research"
description: "Methodology notes across the quantitative pipeline: universe geometry, factor research, feature evaluation, generative models, and alpha evaluation."
---

Every write-up below sits at a key stage of the quantitative pipeline. Rather than judging methodology solely by final portfolio metrics, each stage undergoes rigorous geometric testing, permutation null calibration, and stability verification.

{{< pipeline >}}

## Pipeline Sections & Featured Guides

### 1. [Universe](/quant/universe/)
Tests whether an asset universe carries authentic geometric structure rather than hindsight or random portfolio variance.
- **[Multi-Sub-Universe Directional Coupling Geometry](/quant/universe/multi-sub-universe-directional-coupling-geometry/)** — The restriction Gram operator, canonical $\mathsf V \leftrightarrow \mathsf D$ coupling, and size-matched membership null tests.

### 2. [Factor](/quant/factor/)
Separates exposures, factor returns, factor and specific risk, and each factor's decision role into distinct research objects.
- **[Factor Research: Representation, Returns, Risk, and Neutralization](/quant/factor/factor-research-framework/)** — Exact feature-side factor decomposition in the D/V/T basis, signal- vs portfolio-level neutralization, and factor-return identification under demeaned labels.

### 3. [Feature Evaluation](/quant/feature-evaluation/)
Disentangling raw predictive signal into orthogonal Direction, Activity, and Tail channels, calibrated against permutation nulls.
- **[Cross-Sectional Feature Selection: An Implementation Guide](/quant/feature-evaluation/cross-sectional-feature-selection/)** — Matrix pipeline from whitening to greedy coverage fill and protecting directional exposure.
- **[Cross-Sectional Predictive Geometry: Formal Theory & Operator SVD](/quant/feature-evaluation/cross-sectional-predictive-geometry/)** — Rigorous mathematical propositions, within-date centering, operator SVD, and multiple testing control.
- **[Eligibility, Grouping, and Submodular Coverage: The Selection Bridge](/quant/feature-evaluation/selection-bridge-eligibility-grouping-coverage/)** — Five eligibility gates, mechanism grouping, and submodular coverage optimization.

### 4. [Model](/quant/model/)
Generative modeling of multi-horizon, cross-sectional return distributions rather than point predictions.
- **[Diffusion Models for Statistical Arbitrage: An Implementation Guide](/quant/model/diffusion-statarb-implementation-guide/)** — End-to-end reference covering Transformer denoisers, tensor shapes, guidance, fast sampling, binned likelihoods, and portfolio construction.

### 5. [Alpha Evaluation](/quant/alpha-evaluation/)
Evaluating model predictions with the same scrutiny as raw features.
- **[Label-Alpha Geometry: Significance, Mechanism Inheritance, and Residual Novelty](/quant/alpha-evaluation/label-alpha-geometry/)** — Walk-forward linear surrogates, residual alpha decomposition, and the five diagnostic cases deciding deployment.
