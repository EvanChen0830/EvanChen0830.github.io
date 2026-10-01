---
title: "Cross-Sectional Feature Selection, Part 2: A Worked Toy Example"
description: "End-to-end walkthrough of the Part 1 pipeline on a tiny hand-computable dataset — demeaning, lifting, the cross-covariance operator, mechanism clustering, redundancy families, and stability diagnostics, with actual numbers at every step."
date: 2026-09-10
series: ["Cross-Sectional Feature Selection"]
part: 2
tags: ["feature-selection", "cross-sectional", "methodology", "quant-research", "worked-example"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Feature Selection" >}}

## Part C — Toy Walkthrough (end-to-end, with shapes and numbers)

Setup: $T=2$ days, $N=3$ stocks/day (balanced, for hand-calculation), $P=2$ features, one horizon band with $K=2$. $L=6$. Row order: $(A,1),(B,1),(C,1),(A,2),(B,2),(C,2)$.

**Raw inputs**

$$
X=\begin{pmatrix}1.0&0.2\\2.0&0.4\\3.0&0.3\\1.5&0.1\\2.5&0.5\\2.0&0.6\end{pmatrix},\qquad
Y_{\rm raw}=\begin{pmatrix}0.05&0.10\\0.02&0.01\\-0.01&-0.02\\0.04&0.06\\0.00&-0.01\\0.03&0.02\end{pmatrix}
$$

**B1 — demean.** Day averages: $\bar X_{d1}=(2.0,0.3)$, $\bar X_{d2}=(2.0,0.4)$.

$$
X^{XS}=\begin{pmatrix}-1.0&-0.1\\0.0&0.1\\1.0&0.0\\-0.5&-0.3\\0.5&0.1\\0.0&0.2\end{pmatrix}
$$

(each day's rows sum to zero — check: $-1.0+0.0+1.0=0$, $-0.1+0.1+0.0=0$ ✓). $Y^{XS}$ is computed the same way (omitted, same mechanics).

**B3 — whiten.** Skipping actual $\hat\Sigma_Z$ estimation for brevity, take row 1 of $Z\approx(0.03,0.07)$ as illustration.

**B4 — lift**, row 1:
$$
G_D=(0.03,\,0.07),\quad
G_V=(0.03^2{-}1,\ \sqrt2(0.03)(0.07),\ 0.07^2{-}1)=(-0.999,\,0.00297,\,-0.995),\quad
G_T=0.0762\times(0.03,0.07)
$$
Stack over all 6 rows $\to G\in\mathbb R^{6\times7}$; Gram–Schmidt $\to\Phi\in\mathbb R^{6\times7}$ (in a real run $L\gg d_\Phi$; here $L=6<d_\Phi=7$ is only illustrative of the mechanics, not a usable sample size).

**B5 — cross-covariance.** $C=X^{XS\top}\Phi/6\in\mathbb R^{2\times7}$. Day-decomposition check on $(f_1,\phi_{D_1})$: $C_D^{(1)}[f_1,1]=\tfrac13\big[(-1.0)(0.03)+0+(1.0)(-0.03)\big]=-0.02$, $C_D^{(2)}[f_1,1]=-0.00667$; pooled with $w_1=w_2=0.5$ gives $-0.01333$, matching the direct 6-row computation exactly — confirming (B5').

**B7–B8 — spectral & quality (illustrative numbers for a 6-feature long-band fit).**

| feature | $Q_p$ | $f_p=(D,V,T)$ |
|---|---|---|
| $\text{mom}_{20}$ | 0.0062 | (0.60, 0.02, 0.38) |
| $\text{mom}_{60}$ | 0.0135 | (0.62, 0.04, 0.34) |
| $\text{vol}_5$ | 0.0069 | (0.02, 0.94, 0.04) |
| $\text{turn}_5$ | 0.0050 | (0.00, 0.97, 0.03) |
| $\text{ma}_{10}$ | 0.0086 | (0.61, 0.02, 0.37) |
| $\text{ma}_{12}$ | 0.0074 | (0.62, 0.02, 0.36) |

**B9 — mechanism clustering.** $S^{\rm pred}(\text{mom}_{20},\text{ma}_{10})=0.98$, $S^{\rm pred}(\text{mom}_{20},\text{vol}_5)=0.08$. Cut at 0.8:
$$G_1=\{\text{mom}_{20},\text{mom}_{60},\text{ma}_{10},\text{ma}_{12}\}\ (D\text{-mechanism}),\qquad G_2=\{\text{vol}_5,\text{turn}_5\}\ (V\text{-mechanism}).$$

**B10 — redundancy families.** Within $G_1$: $S^{\rm raw}(\text{ma}_{10},\text{ma}_{12})=0.97$, $S^{\rm raw}(\text{mom}_{20},\text{ma}_{10})=0.93$, $S^{\rm raw}(\text{mom}_{20},\text{mom}_{60})=0.55$. Cut at 0.9: $R_{11}=\{\text{ma}_{10},\text{ma}_{12},\text{mom}_{20}\}$, $R_{12}=\{\text{mom}_{60}\}$. Within $G_2$: $S^{\rm raw}(\text{vol}_5,\text{turn}_5)=0.40$ → each its own family.

**B11–B12 — stability.**

| feature | $\pi_p$ | $M_p$ | $\pi_p^{XS}$ | $L_{\rm eff}$ (of $N{=}500$) |
|---|---|---|---|---|
| mom20 | 0.88 | 0.91 | 0.85 | 210 |
| mom60 | 0.90 | 0.93 | 0.82 | 190 |
| ma10 | 0.85 | 0.89 | 0.60 | **35** |
| vol5 | 0.79 | 0.86 | 0.83 | 240 |
| turn5 | 0.74 | 0.80 | 0.77 | 205 |

**B13 — selection.** In $R_{11}$, ma10's very low $L_{\rm eff}$ marks it as a few-large-caps bet; if hard-filtered, mom20 wins the family by Score. Group quota guarantees $G_2$ keeps at least one representative regardless of $G_1$'s higher raw scores.

$$
\mathcal F^{(\rm long)}_{\rm rep}=\{\text{mom}_{20},\ \text{mom}_{60},\ \text{vol}_5,\ \text{turn}_5\}
$$

Repeating B2–B13 for the short band (independent whitening, independent null block length) might give, e.g., $\mathcal F^{(\rm short)}=\{\text{turn}_5,\text{vol}_5,\text{mom}_{20}\}$ with a more $V$-tilted aggregate fingerprint.

