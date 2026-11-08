---
title: "Factor"
description: "How the factor set used across the pipeline is defined and evaluated: centered point-in-time exposures, characteristic factors, PCA factors and their rotation, alignment, and the checks a factor set must pass before Universe, Feature and Alpha research rely on it."
---

{{< pipeline highlight="factor" >}}

Universe, Feature and Alpha research each ask a factor question in their own
chapter: is this sub-universe a factor tilt, which factor does this feature
predict, on which factor does this alpha have predictive power. This section
defines the factor set those answers are measured against, and how to tell
whether it is stable enough to use. Statistical (PCA) factors get particular
attention: their columns rotate between windows, so only some quantities can
use them directly.

## Guides in this section

- **[Factor Definition and Evaluation](/quant/factor/factor-definition-and-evaluation/)**  
  Centered point-in-time exposures, characteristic factors and their returns, why PCA factors rotate and what stays stable (the subspace), eigengap blocks and Procrustes alignment, choosing the number of factors, the five evaluation checks, and which numbers each later chapter may take from a factor set.
