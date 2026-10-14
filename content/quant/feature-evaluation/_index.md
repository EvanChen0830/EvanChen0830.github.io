---
title: "Feature Evaluation"
description: "Judging which engineered features actually carry cross-sectional predictive signal, and selecting the subset that earns a place in the model — significance, fingerprinting, stability, and horizon-banded selection."
---

{{< pipeline highlight="feature-evaluation" >}}

The Feature Evaluation stage filters and selects features before they enter the model. Rather than relying on univariate rank IC or simple correlation floors, this section evaluates features against rigorous permutation null baselines, decomposes predictive energy into orthogonal mechanisms, and solves a submodular coverage optimization to select a maximally diverse, non-redundant feature subset.

## Guides in this section

- **[Cross-Sectional Feature Selection: An Implementation Guide](/quant/feature-evaluation/cross-sectional-feature-selection/)**  
  An end-to-end matrix implementation guide: within-date demeaning, horizon-band slicing, whitening, lifting into Direction/Activity/Tail blocks, the predictive cross-covariance operator, null calibration, SVD truncation, worked toy numbers, and protecting directional exposure.

- **[Cross-Sectional Predictive Geometry: Formal Theory & Operator SVD](/quant/feature-evaluation/cross-sectional-predictive-geometry/)**  
  A complete formal mathematical treatment: data stacking, within-date demeaning absorption, exact whitening, hierarchical Gram-Schmidt orthogonalization, circular date-permutation null distributions, multiple testing control (Westfall–Young / Permutation FDR), and effective spectral rank.

- **[Eligibility, Grouping, and Submodular Coverage: The Selection Bridge](/quant/feature-evaluation/selection-bridge-eligibility-grouping-coverage/)**  
  The formal bridge between evaluation and model training: five independence-aware eligibility gates, hierarchical clustering for mechanism grouping, a submodularity theorem proving the optimality of greedy coverage, and the direction-first selection algorithm.
