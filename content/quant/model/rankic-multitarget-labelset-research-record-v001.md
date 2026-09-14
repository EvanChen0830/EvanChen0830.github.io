---
title: "Research Record: Multi-Target Label-Set RankIC Benchmark (QD08 v001)"
description: "Does training one cross-sectional RankIC model on a 42-label target set improve its prediction of the tradeable raw return, compared with training on that raw-return label alone? A completed, reproducible research record: frozen data contract, implementation, results, and limitations."
date: 2026-09-14
tags: ["model", "rankic", "multi-task-learning", "research-record", "quant-research"]
math: true
ShowToc: true
TocOpen: true
---

> This is the reference format for a completed model-research report. It states
> the question, frozen data contract, implementation, selection boundary,
> metrics, results, limitations, and exact artifacts required to reproduce the
> conclusion. It is a benchmark record for the RankIC transformer; it does not
> replace the project's production diffusion methodology.

## 1. Decision question

Does training one cross-sectional RankIC model on the complete label set improve
its prediction of the tradeable raw return at a given horizon, compared with a
model trained only on that raw-return label?

The complete target set has seven label families at six horizons, giving 42
outputs:

| Family | Meaning | Horizons |
|---|---|---|
| `raw` | forward return | 5m, 15m, 30m, 1h, 2h, 4h |
| `cs` | cross-sectionally demeaned forward return | same |
| `raw_riskadj` | raw return divided by causal volatility forecast | same |
| `cs_riskadj` | cross-sectionally demeaned return divided by causal volatility forecast | same |
| `pca` | residual after a causal 30-day, two-component PCA factor fit | same |
| `pca_riskadj` | PCA-residual target divided by causal volatility forecast | same |
| `ts_demean` | return after causal one-day per-ticker time-series detrending | same |

The target names are `fwd_ret_{family}_{horizon}`. Their transform definitions
are governed by an internal label specification (`labels.md`). Cross-sectional targets are valid for
cross-sectional prediction: they use only contemporaneous future outcomes when
computing the realized label and are never supplied as a decision-time input.
Risk scaling, PCA state, and time-series detrending use state available at the
decision time.

## 2. Hypothesis and comparison

Auxiliary targets can help when they expose a common cross-sectional predictive
component while reducing target noise. They can hurt when their gradients point
the shared representation toward transformations that are easier to predict but
less aligned with the raw return.

For every primary horizon \(h\), compare:

\[
\begin{aligned}
M_h^{\mathrm{single}} &: X_{t,h}^{(\mathrm{selected})}
  \longrightarrow Y_{t}^{\mathrm{raw},h},\\
M_h^{42} &: X_{t,h}^{(\mathrm{selected})}
  \longrightarrow \left[Y_t^{\ell}\right]_{\ell=1}^{42}.
\end{aligned}
\]

The reported primary metric for both models is the Pearson IC of the raw output
\(\hat Y_t^{\mathrm{raw},h}\) against \(Y_t^{\mathrm{raw},h}\). Therefore the
comparison changes the training target set only; it does not compare different
features, universe, split, architecture, or raw-return metric.

## 3. Data and split contract

| Item | Value |
|---|---|
| Feature set | `feature_book_1m_streaming_full_v003_ensemble_current_v001` |
| Feature inventory | 988 features |
| Universe | 20 assets, identical across compared runs |
| Bar frequency | one minute inputs; decision timestamps in cached panel |
| Train | 2024-07-01 through 2025-07-16 |
| Validation | 2025-07-17 through 2025-10-12 |
| Test | 2025-10-13 through 2026-09-07 22:30 UTC |
| Purge between splits | one day |
| Golden OOS | not accessed |

All tensors were loaded once into
`artifacts/model_inputs/qd08_rankic_full_pca_v001/`. The complete cache stores
`all_labels.npy` and `all_labels_meta.json`, containing the 42 labels. Training
uses memory-mapped local arrays and does not read the database. The cache has
the same timestamp and asset order as `features.npy`.

## 4. Feature subset contract

The previous spectral/mechanism experiment formed feature groups using training
data, then used validation IC to choose one subgroup for each raw horizon. This
benchmark freezes those choices before fitting either side of the comparison.

| Primary horizon | Frozen subgroup | Features | Source mechanism groups |
|---|---:|---:|---|
| 5m | 3 | 226 | 6 |
| 15m | 2 | 177 | 4 |
| 30m | 1 | 71 | 1 |
| 1h | 2 | 41 | 2 |
| 2h | 1 | 15 | 1, 2 |
| 4h | 1 | 38 | 3, 1, 2 |

The manifests are under
`artifacts/feature_selections/qd08_standard_cs_raw_{horizon}_v001/` and the
summary is `artifacts/feature_selections/qd08_standard_selection_summary.json`.

This keeps the test comparison clean. Validation IC remains useful for model
selection, but it is not an independent estimate because it participated in
the earlier subgroup choice. Test IC is the selection-independent comparison
for this record.

## 5. Model and objective

This is a point-prediction diagnostic model, not the production diffusion
model. Each timestamp is a cross-section of 20 assets. A ticker embedding and
per-ticker feature projection feed a three-layer transformer that attends across
tickers, followed by one output head per target.

| Parameter | Value |
|---|---|
| Hidden width | 128 |
| Attention heads | 4 |
| Transformer layers | 3 |
| Dropout | 0.05 |
| Optimizer | AdamW |
| Learning rate | 0.0003 |
| Weight decay | 0.0001 |
| Batch size | 256 timestamps |
| Maximum epochs | 30 |
| Early-stopping patience | 5 epochs |
| Random seed | 7 |

Features are cross-sectionally demeaned at each timestamp and standardized
using train-split moments only. For a score matrix \(S\in\mathbb R^{N\times L}\)
and target matrix \(Y\in\mathbb R^{N\times L}\), the loss for a timestamp is:

\[
\mathcal L_t = 1 - \frac{1}{L}\sum_{\ell=1}^{L}
\operatorname{corr}_{i=1}^{N}\left(S_{t,i,\ell},Y_{t,i,\ell}\right).
\]

The training loss averages \(\mathcal L_t\) over timestamps and all 42 output
labels. Centering is performed along the asset axis, so the loss is exactly a
cross-sectional Pearson objective. It does not optimize hit rate or a
portfolio objective.

After training, each output is calibrated using a zero-intercept OLS slope on
**train data only**, after the same cross-sectional centering. Calibration
cannot affect Pearson IC, but it retains a consistent score-to-label contract.

## 6. Implementation changes

The RankIC runner was extended in
`packages/model/src/quantlab/model/rankic.py` with two generic capabilities:

1. A label-cache resolver now supports the complete 42-label cache in addition
   to the existing model-label and raw-label caches. It checks requested label
   names against cache metadata before opening the array.
2. `inference.write_scenarios: false` suppresses residual-bootstrap scenario
   materialization for IC-only experiments. A 42-label scenario panel would be
   large and does not contribute to this comparison. Checkpoints, configuration
   snapshots, metadata, and `signal_ic.json` are always written.

These are execution changes, not changes to data or model mathematics. The
Pearson objective and cross-sectional calibration were already the active
RankIC methodology before this experiment.

## 7. Results

### 7.1 Primary raw-return result

| Horizon | Single-label valid IC | 42-label valid IC | Δ valid | Single-label test IC | 42-label test IC | Δ test |
|---|---:|---:|---:|---:|---:|---:|
| 5m | 0.031971 | 0.031254 | -0.000718 | 0.034400 | 0.031974 | -0.002426 |
| 15m | 0.039335 | 0.034026 | -0.005309 | 0.032026 | 0.027963 | -0.004062 |
| 30m | 0.036359 | 0.036477 | +0.000119 | 0.028030 | 0.026745 | -0.001286 |
| 1h | 0.034511 | 0.037934 | +0.003422 | 0.021114 | 0.022782 | +0.001668 |
| 2h | 0.036176 | 0.033934 | -0.002242 | 0.021824 | 0.022244 | +0.000420 |
| 4h | 0.021688 | 0.019393 | -0.002294 | 0.005815 | 0.007834 | +0.002019 |

The 42-label objective improves test raw IC at 1h, 2h, and 4h. It degrades
5m, 15m, and 30m. The result is horizon-dependent; a single global choice of
single-target or all-target training is not supported by these data.

### 7.2 Auxiliary-label diagnostic

The following are the direct test ICs of outputs at the same primary horizon
from the 42-label model. Values across label families answer different target
questions and should not be treated as a common economic scale.

| Model horizon | Raw | CS | Raw risk-adj. | CS risk-adj. | PCA | PCA risk-adj. | TS demeaned |
|---|---:|---:|---:|---:|---:|---:|---:|
| 5m | 0.031974 | 0.032412 | 0.029009 | 0.031905 | 0.029720 | 0.028473 | 0.047074 |
| 15m | 0.027963 | 0.027630 | 0.022482 | 0.026721 | 0.023238 | 0.022652 | 0.059188 |
| 30m | 0.026745 | 0.027001 | 0.022512 | 0.026693 | 0.021804 | 0.023112 | 0.054385 |
| 1h | 0.022782 | 0.023140 | 0.012916 | 0.021141 | 0.014455 | 0.013425 | 0.048497 |
| 2h | 0.022244 | 0.021828 | 0.011893 | 0.020392 | 0.011779 | 0.011331 | 0.028585 |
| 4h | 0.007834 | 0.010032 | 0.000018 | 0.007709 | 0.001139 | 0.001063 | 0.173814 |

The high time-series-demeaned IC, especially at 4h, says that this transformed
target is easier for this architecture to rank. It does **not** establish that
it is a stronger tradeable alpha than raw return. That requires a separate,
properly selected portfolio experiment, which this record deliberately does
not run.

## 8. Decision and next experiment

Use the all-label RankIC variant as a candidate for the 1h, 2h, and 4h raw
signals. Retain the single-label variant as the current candidate for 5m, 15m,
and 30m. This is a research decision based on the untouched standard test
split; it does not promote a model to production or claim PnL improvement.

The next controlled experiment should replace equal weighting across the 42
Pearson losses with weights selected on training/validation only. Candidate
families are:

- horizon-local auxiliary targets, to test whether near-horizon labels explain
  the 1h+ improvement without short-horizon interference;
- family weights learned on validation and frozen before the test evaluation;
- a Pareto comparison of raw IC and portfolio utility after a separately
  authorized portfolio run.

Do not use this test result to choose these weights. Create a new versioned
configuration and reserve an untouched evaluation window for that decision.

## 9. Reproducibility manifest

The six complete-model configurations are:

```text
configs/models/qd08_standard_rankic_5m_labelset42_v001.yaml
configs/models/qd08_standard_rankic_15m_labelset42_v001.yaml
configs/models/qd08_standard_rankic_30m_labelset42_v001.yaml
configs/models/qd08_standard_rankic_1h_labelset42_v001.yaml
configs/models/qd08_standard_rankic_2h_labelset42_v001.yaml
configs/models/qd08_standard_rankic_4h_labelset42_v001.yaml
```

Each run persists `best.pt`, `meta.json`, its copied YAML config, and
`signal_ic.json` under `artifacts/models/qd08_standard_rankic_{horizon}_labelset42_v001/`.
The compact machine-generated result table is
`artifacts/reports/qd08_standard_labelset42_comparison_v001.md`.

To reproduce one run from the project root:

```bash
uv run --frozen python -c \
  "from quantlab.model import RankICModelPipeline; \
   RankICModelPipeline().run('configs/models/qd08_standard_rankic_1h_labelset42_v001.yaml')"
```

## 10. Report-writing checklist

A future research report must include all of the following before it is used
for a modeling decision:

1. The concrete decision question and a falsifiable comparison.
2. The exact data, universe, labels, split boundaries, purge, and whether any
   golden OOS data was accessed.
3. The selection boundary: what was learned on train, chosen on validation,
   and evaluated once on test.
4. All changes to feature inputs, target definitions, architecture, objective,
   calibration, and inference output.
5. Primary metrics with both validation and test values, plus uncertainty when
   the metric supplies it.
6. A conclusion restricted to what the comparison identifies, with limitations
   and the next experiment stated explicitly.
7. Versioned configs, artifacts, code locations, seed, and a command that can
   recreate the run without an unrecorded database read.
