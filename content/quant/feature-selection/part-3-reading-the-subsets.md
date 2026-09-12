---
title: "Cross-Sectional Feature Selection, Part 3: Reading the Resulting Subsets"
description: "How to read the audit block attached to each horizon-band subset — Direction/Activity/Tail coverage, effective rank, and per-feature concentration — and what each pattern implies for how the subset should actually be used."
date: 2026-09-11
series: ["Cross-Sectional Feature Selection"]
part: 3
tags: ["feature-selection", "cross-sectional", "methodology", "quant-research", "diagnostics"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Feature Selection" >}}

## Part D — Reading the $n$ Subsets

Each $\mathcal F^{(b)}$ comes with its own audit block (Coverage, block coverage, $r_{\rm eff}$, stability). Use these to characterize what each subset is actually good for:

| Diagnostic | What it tells you |
|---|---|
| $\operatorname{Cov}_D(\mathcal F^{(b)})$ high, $\operatorname{Cov}_V,\operatorname{Cov}_T$ low | Subset is **directional**: predicts the sign and magnitude of relative outperformance. Suited to long-short position sizing — the classic "relative strength" portfolio. |
| $\operatorname{Cov}_V(\mathcal F^{(b)})$ high | Subset predicts future dispersion/co-movement structure, not who wins. Useful for risk sizing, pairs/dispersion trades, not a directional bet. |
| $\operatorname{Cov}_T(\mathcal F^{(b)})$ high | Subset predicts direction conditional on large moves only. Useful for tail-hedge timing or crash-alpha, weak signal in calm regimes. |
| $r^{(b)}$, $r_{\rm eff}^{(b)}$ across bands | How many independent mechanisms the library supports at that horizon. A band with $r_{\rm eff}\approx1$ has one alpha direction no matter how many features are in $\mathcal F^{(b)}$; compare across bands to see whether independent mechanisms emerge only at certain horizons. |
| Which $G_m$ are represented in $\mathcal F^{(b)}$ vs. which existed in $\mathcal E^{(b)}$ but got squeezed out | Tells you whether the band's budget was directional-sufficient or diluted — see Part E. |
| $L_{\rm eff}(p)$ of the features actually selected | Whether the band's directional power rests on a broad cross-section or a handful of large caps — matters directly for tradability/capacity. |

**Reading example** (continuing the toy numbers): a long-band subset with fingerprint aggregate near $(0.55,0.10,0.35)$ is directional-with-a-tail-tilt — good for a long-short book that also wants to lean into big moves; a short-band subset near $(0.05,0.85,0.10)$ is not a relative-strength book at all — it's an activity/volatility-timing signal, and using it for position sizing on direction would be a category error even though its $Q_p$ values look respectable.

**Comparing bands.** If $\operatorname{Cov}_D$ decays sharply from short to long bands while $\operatorname{Cov}_V$ or $\operatorname{Cov}_T$ holds up, it means your directional relative-strength edge decays faster than your activity/tail edge — a genuinely useful finding for holding-period design, not an artifact to explain away.

