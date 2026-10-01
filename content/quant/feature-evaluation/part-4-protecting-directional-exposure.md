---
title: "Cross-Sectional Feature Selection, Part 4: Protecting Directional Exposure"
description: "Why the fairness mechanism in Stage 2 selection can dilute relative-strength predictive power, and five concrete levers to protect it when Direction is the channel that actually matters for the portfolio."
date: 2026-09-12
series: ["Cross-Sectional Feature Selection"]
part: 4
tags: ["feature-selection", "cross-sectional", "methodology", "quant-research", "portfolio-construction"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Feature Selection" >}}

## Part E — Protecting Directional (Relative-Strength) Exposure in Selection

**First, confirm what "relative strength magnitude" maps to in this framework.** $\Phi_D=z$ is the whitened, within-date-demeaned, *continuous* future return — not a rank transform. So $Q_p^D$ already measures a feature's ability to predict the **size** of relative over/underperformance, not merely its sign. $\operatorname{Cov}_D(\mathcal F)$ at the set level is therefore the single number that most directly answers "how much of the cross-sectional relative-strength signal did I keep." Treat it as the primary KPI in every report (Section 26/31 of the main document), not merely one of three equally-weighted blocks.

The risk you're worried about is specifically in **Stage 2** (B13c): the block-normalized greedy fill was designed to *stop* Direction from crowding out Activity/Tail. That is exactly the mechanism to relax when Direction is what you actually want to protect. Group quota (the floor in B13c) is **not** the culprit — it only guarantees a minimum for $G_2,G_3,\dots$; it never caps how many representatives $G_1$ (the D-mechanism group) can keep.

Concrete levers, roughly in order of how directly they target the problem:

**1. Lexicographic budget allocation.** Replace the single water-filling pass with two passes: first spend budget purely on maximizing $\Delta\mathrm{Cov}_D(p)$ among D-eligible candidates until $\operatorname{Cov}_D(\mathcal F)\ge\gamma_D$ (e.g. 90–95%) or the sub-budget is exhausted; only then spend remaining slots on $V,T$ with the original water-filling logic. This guarantees Direction coverage is never sacrificed to satisfy fairness toward the other two channels — fairness only governs the leftover budget.

**2. Direction-specific scoring track.** Within a D-mechanism group, rank candidates by $\bar Q_p^D$ (block-specific quality) rather than total $\bar Q_p$:
$$
\mathrm{Score}_p^D=\bar Q_p^D\cdot\pi_p^{w_\pi}\cdot M_p^{w_M}\cdot g(E_p).
$$
This prevents a feature with strong but incidental $V$/$T$ content from outscoring a purer directional feature purely because its *total* $Q_p$ happens to be larger.

**3. Loosen the redundancy threshold inside D-mechanism groups.** Raising $\tau_{\rm raw}$ specifically for families formed within $G_1$ (e.g., require $S^{\rm raw}\ge0.95$ instead of $0.90$ before two features merge into one family) produces *more*, smaller families — hence more surviving representatives after Stage 1. This is the direct lever for "don't prune my directional features so aggressively," distinct from and complementary to lever 1.

**4. Turn $L_{\rm eff}$ into a soft discount, not a hard filter, when direction is the priority.** A hard cutoff $L_{\rm eff}(p)>L_{\min}$ can throw away a concentrated-but-genuinely-strong directional feature even when nothing else in its family comes close. Prefer
$$
g_L(p)=\min\!\left(1,\ L_{\rm eff}(p)/L_0\right),
$$
multiplied into $\mathrm{Score}_p$ (or $\mathrm{Score}_p^D$), so concentration discounts a candidate proportionally rather than disqualifying it outright — unless your capacity constraints genuinely require a hard floor, in which case keep it, but set $L_{\min}$ relative to your own tradable universe size rather than an arbitrary global number.

**5. Report $\operatorname{Cov}_D$ per feature-slot spent.** $\operatorname{Cov}_D(\mathcal F)/|\mathcal F|$ (or the marginal $\Delta\mathrm{Cov}_D$ of the last few features admitted) tells you whether your budget is still buying meaningful directional coverage or whether you're spending slots on diminishing returns — useful for deciding how big $P_{\max}$ needs to be before Direction saturates.

**Recommended default combination:** run lever 1 (lexicographic budget split) as the structural fix, lever 3 (looser $\tau_{\rm raw}$ inside D groups) so Stage 1 doesn't over-consolidate before Stage 2 even runs, and lever 4 (soft $L_{\rm eff}$ discount) so a single concentration diagnostic can't unilaterally veto your strongest directional candidate. Use lever 2 as a secondary tie-break inside $G_1$ rather than replacing $\mathrm{Score}_p$ everywhere — total-quality scoring is still the right criterion for $G_2,G_3,\dots$, since there Direction *isn't* what you're optimizing for.
