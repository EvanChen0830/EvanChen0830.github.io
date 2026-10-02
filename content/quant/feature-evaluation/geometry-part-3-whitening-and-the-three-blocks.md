---
title: "Cross-Sectional Predictive Geometry, Part 3: Whitening, the Three Predictive-Mechanism Blocks, and Orthogonalization"
description: "How the label mechanism basis gets built: rank-controlled whitening, the Direction/Activity/Tail lifting maps and why each coordinate means what it means, hierarchical Gram-Schmidt orthogonalization in a fixed D-then-V-then-L order, and the fit/transform boundary that separates research from backtest."
date: 2026-09-14
series: ["Cross-Sectional Predictive Geometry"]
part: 3
tags: ["cross-sectional", "methodology", "quant-research", "predictive-geometry", "whitening"]
math: true
ShowToc: true
TocOpen: true
---

{{< series name="Cross-Sectional Predictive Geometry" >}}

## §7. Label Whitening and Rank Control

### Motivation

Whitening's purpose is to put different horizons' coordinates on the same
scale and make them mutually uncorrelated, so that "which horizon
combination did this feature predict" has a unique, scale-free answer.
Without whitening, the horizon with the largest variance in $\Sigma_Y$
(usually the longest) would dominate every inner product.

But whitening is also the **only step in the entire method that amplifies
noise**. Along $\Sigma_Y$'s $j$-th eigendirection, the data is multiplied by
$\lambda_j^{-1/2}$. If $\lambda_K/\lambda_1=10^{-2}$, the weakest direction
is amplified 10×, and so is the estimation error in that direction.
Whitening must therefore be paired with rank control, and the rank-control
threshold must be set by statistical stability, not an arbitrary numeric
floor.

### §7.1 Covariance and Spectrum

$$\boxed{\Sigma_Y=\frac1T\sum_{d=1}^T\frac1{N_d}Y_d^{(\beta)\top}Y_d^{(\beta)}\in\mathbb R^{K\times K}.}$$

(By the daily demeaning in §5, this genuinely is a covariance, not an
uncentered moment.) Eigendecomposition:

$$\Sigma_Y=E\Lambda E^\top,\qquad \lambda_1\ge\lambda_2\ge\cdots\ge\lambda_K>0.$$

### §7.2 Retained Dimension $K_{\rm eff}$

Moving block bootstrap over dates (block length $\ge H_{\max}$,
$B_\lambda=200$ draws), giving $\{\lambda_j^{(b)}\}$. Define

$$\boxed{K_{\rm eff}=\max\left\{k:\ \frac{\lambda_k}{\lambda_1}\ge\frac1{\kappa_{\max}}\ \ \text{and}\ \ \frac{\operatorname{sd}_b(\lambda_k^{(b)})}{\lambda_k}\le\tau_\lambda\right\},}$$

default $\kappa_{\max}=100$ (max 10× amplification), $\tau_\lambda=0.25$.

Both conditions are necessary: the first controls numerical amplification,
the second controls statistical instability. A direction with
$\lambda_k/\lambda_1=0.05$ but 60% bootstrap relative standard deviation has
a whitened coordinate that's almost entirely estimation noise, even though
it clears the condition-number threshold.

**The report must include**: the full $\{\lambda_j/\lambda_1\}$ spectrum,
$\kappa=\lambda_1/\lambda_{K_{\rm eff}}$, $K_{\rm eff}$, and the horizon
loadings $E[:,K_{\rm eff}{+}1:K]$ of the discarded directions (to see what
was actually discarded).

If $K_{\rm eff}\lt K$ and the discarded directions' loadings concentrate on
one or two horizons, the right fix is usually not to accept the
truncation — it's to go back to §6 and switch to an incremental basis or
merge horizons. Truncation is the last resort.

### §7.3 The Whitening Matrix

**Case one (no truncation, $K_{\rm eff}=K$)**: use the symmetric inverse
square root

$$B_Y=E\Lambda^{-1/2}E^\top=\Sigma_Y^{-1/2}\in\mathbb R^{K\times K}.$$

The advantage of the symmetric square root (over Cholesky) is that it
**doesn't depend on horizon ordering**, so $Z$'s $k$-th coordinate still
roughly corresponds to the $k$-th horizon when there's no truncation —
good interpretability. Cholesky whitening would turn the coordinates into a
triangular structure ("horizon 1," "horizon 2 net of horizon 1," ...) that
isn't comparable across bands.

**Case two (with truncation, $K_{\rm eff}\lt K$)**: use the reduced-dimension
form

$$B_Y=E_{1:K_{\rm eff}}\Lambda_{1:K_{\rm eff}}^{-1/2}\in\mathbb R^{K\times K_{\rm eff}}.$$

Here the coordinates' meaning shifts to "the leading $K_{\rm eff}$ principal
directions of future returns" and no longer correspond to individual
horizons. This semantic change must be stated explicitly in the report,
since each component of $c_p$ will be read against $E$'s loadings.

Both cases define

$$\boxed{Z_d=Y_d^{(\beta)}B_Y.}$$

### §7.4 Proposition 7.1 (exact whitening)

$$Z^\top WZ=I_{K_{\rm eff}}\quad\text{(exact, not approximate)}.$$

**Proof**: $Z^\top WZ=B_Y^\top\big(\frac1T\sum_d\frac1{N_d}Y_d^\top
Y_d\big)B_Y=B_Y^\top\Sigma_YB_Y$. In case two, $B_Y^\top\Sigma_YB_Y=
\Lambda^{-1/2}E^\top E\Lambda E^\top E\Lambda^{-1/2}$ (taking the first
$K_{\rm eff}$ columns) $=\Lambda^{-1/2}\Lambda\Lambda^{-1/2}=I$. Case one is
the special case. $\blacksquare$

This proposition is the foundation of §9: precisely because the whitening
is exact, $\Phi_{\mathsf D}$ can be taken directly as $Z$ with no further
orthogonalization — and therefore introduces no coordinate ambiguity.

Also, from §5's daily demeaning, $\sum_iZ_{id,:}=0$.

> **Convention**: after §7, the symbol $K$ always means $K_{\rm eff}$, the
> actually-retained label dimension.

## §8. The Three Predictive-Mechanism Blocks

### Overall Motivation

Given a whitened future-return vector $z\in\mathbb R^K$, we ask: which
**functions** of $z$ might a cross-sectional signal plausibly relate to?

In principle the answer is "any function of $z$," but that's an
infinite-dimensional space. This method picks three concrete lifting maps,
because they correspond to three prediction patterns that are completely
different in trading terms — requiring completely different handling:

| block | lifting | trading meaning | if a model only learned this |
|---|---|---|---|
| $\mathsf D$ | $z$ | who outperforms whom | long/short directly |
| $\mathsf V$ | $\operatorname{svec}(zz^\top)$ | who moves more violently, how horizons co-move | volatility/gamma positions only, no direction |
| $\mathsf L$ | $\|z\|z$ | whose direction is right during large moves | directional but heavily conditional; sizing and risk management are entirely different |

This table is the reason the whole method exists: under ordinary rank IC
these three contaminate each other; in this coordinate system they are
strictly separated.

### §8.1 The Three Liftings and Within-Date Centering

For each observation's $z_{id}\in\mathbb R^K$:

$$G_{\mathsf D}(z)=z,\qquad G_{\mathsf V}(z)=\operatorname{svec}(zz^\top),\qquad G_{\mathsf L}(z)=\|z\|_2\,z.$$

Dimensions are $K$, $K(K{+}1)/2$, $K$ respectively.

Then **all three are always within-date centered**:

$$\boxed{G_{\mathsf B,d}\leftarrow G_{\mathsf B,d}-\mathbf 1_{N_d}\,\overline{G_{\mathsf B,d}},\qquad \overline{G_{\mathsf B,d}}=\frac1{N_d}\sum_iG_{\mathsf B}(z_{id}).}$$

$G_{\mathsf D}$ already satisfies this (since $\sum_iz_{id}=0$); $G_{\mathsf
V}$ and $G_{\mathsf L}$ don't, and it must be done explicitly.

**Why this is mandatory.** By Lemma 4.1, these "that-period-constant"
components always contribute zero to $C$, so skipping centering **does not
produce a wrong $C$**. But it causes two real losses:

1. **Wasted basis dimension.** $\operatorname{orth}_W$ acting on a matrix
   with a day-common component produces some $\Phi_{\mathsf V}$ directions
   that are mostly "the market volatility time series." These directions
   are orthogonal to any demeaned $X$, forever contributing $Q=0$, yet
   occupy a $d_{\mathsf V}$ slot — diluting §19's fingerprint denominator
   and creating meaningless small singular values in §12's SVD.
2. **Distorted orthogonalization scale.** $\Xi_R=R^\top WR$ gets dominated
   by day-common variance, so $\Lambda_R^{-1/2}$ compresses genuine
   cross-sectional variance to the wrong scale, making the
   $\mathsf D\to\mathsf V\to\mathsf L$ Gram–Schmidt numerically less stable.

**Centering simplifies the definition.** Note the raw Activity lifting is
often written $\operatorname{svec}(zz^\top-I_K)$. After centering, $-I_K$ is
a global constant absorbed entirely by the that-period mean, so

$$\text{center}\big(\operatorname{svec}(zz^\top-I_K)\big)=\text{center}\big(\operatorname{svec}(zz^\top)\big).$$

This book therefore writes $G_{\mathsf V}(z)=\operatorname{svec}(zz^\top)$
directly.

**Semantic statement.** After this step, all of $\Phi$ lives in the
**within-date cross-sectional function space**. That's correct — this
method, start to finish, only cares about cross-sectional contrasts;
anything "common to all stocks that day" is outside the scope of study.

### §8.2 Definition and Verification of svec

`svec` flattens a symmetric matrix into a vector: diagonal entries unchanged,
off-diagonal entries multiplied by $\sqrt2$, upper triangle only.

**Property**: for any symmetric $A,B$,

$$\langle\operatorname{svec}A,\operatorname{svec}B\rangle=\operatorname{tr}(A^\top B).$$

**Verification**: $\operatorname{tr}(A^\top
B)=\sum_hA_{hh}B_{hh}+2\sum_{h<k}A_{hk}B_{hk}$, while
$\langle\operatorname{svec}A,\operatorname{svec}B\rangle=\sum_hA_{hh}B_{hh}+
\sum_{h<k}(\sqrt2A_{hk})(\sqrt2B_{hk})$ — the two are equal. $\blacksquare$

This property guarantees that the $\mathsf V$ block's Euclidean geometry
matches the matrix's Frobenius geometry, so $\|c_p^{\mathsf V}\|^2$ can be
read directly as "the strength of the feature's association with
second-order structure," with no coordinate artifact.

### §8.3 Content of the $\mathsf V$ Block

**Worked example ($K=2$, horizons being the 1h and 4h increments)**:
$G_{\mathsf V}$'s three coordinates are $(z_1^2,\ \sqrt2z_1z_2,\ z_2^2)$,
each minus its that-period mean after centering.

- Coordinates one and three: is this stock's 1h / 4h **magnitude** above
  the typical cross-sectional value that period. This is "activity."
- Coordinate two: the **co-directionality** of this stock's 1h and 4h
  returns. If 1h is up and 4h is also up, $z_1z_2>0$; if 1h is up and 4h
  gives it back, $z_1z_2<0$.

So the second coordinate is actually measuring **trend continuation vs.
mean reversion** — a feature that predicts it is, in substance, predicting
whether this stock's near-term return structure will be continuation-type
or reversal-type. This is valuable information, and it's completely
invisible in any single-horizon IC. This dimension is also the most
commonly misread: it is not volatility, it is **sign consistency across
horizons**.

### §8.4 Content of the $\mathsf L$ Block

$G_{\mathsf L}(z)=\|z\|z$ and $G_{\mathsf D}(z)=z$ are obviously highly
correlated (same direction, just reweighted by $\|z\|$). The $\mathsf L$
block's actual content is $G_{\mathsf L}$'s **residual** after projecting
out $\mathsf D$ and $\mathsf V$, intuitively approximately

$$R_{\mathsf L}\approx\big(\|z\|-\mathbb E\|z\|\big)\,z\ \text{, orthogonalized},$$

i.e. "**magnitude surprise × direction**."

A feature loading on the $\mathsf L$ block means: it has no directional
predictive power during calm periods, but it correctly calls direction on
the names that experience large moves. Such a feature's role in a portfolio
is a **conditional sizing signal**, not a baseline ranking signal. Its
position design in practice is completely different: it shouldn't get
equal-weighted long/short treatment — it should be used as a size modifier
on an existing directional position.

---

## §9. Hierarchical Orthogonalization

### Motivation

The three liftings' image spaces overlap. Used directly, $Q^{\mathsf D}+Q^{\mathsf
V}+Q^{\mathsf L}$ would double-count the same information, and the
fingerprint would become meaningless. Hierarchical orthogonalization
resolves this, at the cost of introducing an **order dependency**.

> **Principle 9.1 (attribution order)**: fixing the order
> $\mathsf D\to\mathsf V\to\mathsf L$ means "anything explainable by
> first-order direction is always credited to $\mathsf D$." So the correct
> reading of $Q^{\mathsf V}$ is **not** "how much second-order information
> does this feature have," but "how much second-order information does this
> feature have **that cannot be explained by first-order direction**."

This distinction is critical when reading a report. A feature that predicts
both direction and volatility will have a lower $Q^{\mathsf V}$ — not
because it doesn't predict volatility, but because its volatility
prediction is collinear with its directional prediction. The order is fixed
this way because directional information is the most directly tradeable; by
crediting the shared component to it, we ensure that whatever loads on
$\mathsf V$ and $\mathsf L$ is guaranteed to be "additional" information
that genuinely needs different handling.

### §9.1 The Weighted Orthogonalization Operator

Given a matrix $R$ (rows are observations), compute
$\Xi_R=R^\top WR$, eigendecompose $\Xi_R=U_R\Lambda_RU_R^\top$, retain
$\lambda_j>\epsilon_{\rm rank}\lambda_{\max}$ ($\epsilon_{\rm rank}=10^{-10}$,
purely numerical), and define

$$\Theta_R=U_R\Lambda_R^{-1/2},\qquad \operatorname{orth}_W(R)=R\,\Theta_R.$$

**Verification**: $(R\Theta_R)^\top W(R\Theta_R)=\Lambda_R^{-1/2}U_R^\top
\Xi_RU_R\Lambda_R^{-1/2}=I$. $\blacksquare$

If any dimension gets truncated by $\epsilon_{\rm rank}$, it must be
recorded and reported, since that means the lifting's image space actually
has lower rank than its nominal dimension (which can happen when $K$ is
large and $T$ is small).

### §9.2 The $\mathsf D$ Block: Take $Z$ Directly

$$\boxed{\Phi_{\mathsf D}:=Z.}$$

**No orthogonalization is applied.** By Proposition 7.1, $Z^\top WZ=I_K$
already holds exactly, so $\operatorname{orth}_W$ would be redundant — and
not merely redundant, actively harmful. Here's why.

If $\operatorname{orth}_W$ were applied to $Z$, then $\Xi_R=I_K$, whose
eigendecomposition is $I=U_RIU_R^\top$ — and **any orthogonal matrix $U_R$
is a valid eigenvector matrix**. So $\Theta_R=U_R$, giving

$$\Phi_{\mathsf D}=ZU_R,$$

where $U_R$ is determined by LAPACK implementation details, thread count,
even memory alignment. This is harmless for $Q^{\mathsf D}=\|c_p^{\mathsf
D}\|^2$ and any angle (orthogonal transforms preserve norms and inner
products), but it's fatal for §21's $A^{\rm tgt}_p=(c_p^{\mathsf
D})^\top r_{\mathsf D}$: if $r_{\mathsf D}$ is defined in one coordinate
system and $c_p^{\mathsf D}$ in another, the two differ by an unknown
rotation, and the alignment value would vary across execution environments.

Setting $\Phi_{\mathsf D}:=Z$ eliminates this ambiguity once and for all.
The $\mathsf D$ block's coordinates therefore have a clear economic
interpretation: the $k$-th coordinate is simply "the future-return vector's
component along the $k$-th whitened direction," roughly corresponding to
the $k$-th horizon when there's no truncation.

### §9.3 The $\mathsf V$ and $\mathsf L$ Blocks

$$\Gamma_{\mathsf{DV}}=\Phi_{\mathsf D}^\top WG_{\mathsf V},\qquad R_{\mathsf V}=G_{\mathsf V}-\Phi_{\mathsf D}\Gamma_{\mathsf{DV}},\qquad\boxed{\Phi_{\mathsf V}=\operatorname{orth}_W(R_{\mathsf V}).}$$

$$\Gamma_{\mathsf{DL}}=\Phi_{\mathsf D}^\top WG_{\mathsf L},\quad\Gamma_{\mathsf{VL}}=\Phi_{\mathsf V}^\top WG_{\mathsf L},\quad R_{\mathsf L}=G_{\mathsf L}-\Phi_{\mathsf D}\Gamma_{\mathsf{DL}}-\Phi_{\mathsf V}\Gamma_{\mathsf{VL}},$$

$$\boxed{\Phi_{\mathsf L}=\operatorname{orth}_W(R_{\mathsf L}).}$$

Finally

$$\boxed{\Phi=[\Phi_{\mathsf D}\,|\,\Phi_{\mathsf V}\,|\,\Phi_{\mathsf L}],\qquad \Phi^\top W\Phi=I_{d_\Phi}.}$$

### §9.4 Proposition 9.2 (within-block rotation invariance)

$\Phi_{\mathsf V}$ and $\Phi_{\mathsf L}$ still have internal rotation
freedom (eigenvector ordering and sign in $\operatorname{orth}_W$). This
causes no problem, because:

> **Proposition 9.2**: let $O=\operatorname{diag}(I_{d_{\mathsf
> D}},O_{\mathsf V},O_{\mathsf L})$ be any block-diagonal orthogonal matrix.
> Replacing $\Phi$ with $\Phi O$,
> $$C\mapsto CO,\quad V_r\mapsto O^\top V_r,\quad c_s\mapsto O^\top c_s,$$
> while the following are **completely invariant**: $Q_s$, $Q_s^{\mathsf
> D}$, $Q_s^{\mathsf V}$, $Q_s^{\mathsf L}$, $f_s$, $\tilde f_s$, $S^{\rm
> mech}_{pq}$, $A^{\rm mech}_{pq}$, $A^{\rm tgt}_s$, $\eta^{\mathsf D}_s$,
> $N^{\rm eff}$, $R^{\rm canc}$, $r_{\rm sig}$, $r_{\rm eff}$, and every
> $\pi$ and $M$.

**Proof**: $C\mapsto CO$ follows directly from the definition. The SVD
becomes $CO=US(O^\top V)^\top$, so singular values are unchanged,
$V\mapsto O^\top V$, hence $r_{\rm sig}$ and $r_{\rm eff}$ are unchanged;
$C_{\rm sig}=CV_rV_r^\top\mapsto CO\cdot O^\top V_rV_r^\top O=C_{\rm sig}O$,
so $c_s\mapsto O^\top c_s$, and norms and any pairwise inner product are
unchanged. Since $O$ is the identity on the $\mathsf D$ block, $c_s^{\mathsf
D}$ is **unchanged componentwise**, so $A^{\rm tgt}$ and $\eta^{\mathsf D}$
are unchanged. Each block's energy $\|c_s^{\mathsf B}\|^2$ is preserved
individually because $O$ is block-diagonal. In $g_{id}=X_{id}\langle
\Phi_{id},\hat c\rangle$, $\langle\Phi_{id}O,O^\top\hat c\rangle=\langle
\Phi_{id},\hat c\rangle$, so breadth and cancellation are unchanged.
$\blacksquare$

This proposition is the formal statement of the whole method being
"well-defined": **every reported quantity depends only on the subspace
spanned by each of $\Phi$'s three blocks, not on the choice of basis within
a block; and the $\mathsf D$ block's basis is pinned down by
$\Phi_{\mathsf D}=Z$, which in turn pins down the one quantity that has
directional meaning (target alignment).**

In implementation this gives a free unit test: rerun the whole pipeline
with different BLAS thread counts or random seeds — every quantity above
must come out identical up to numerical error (§65).

---

## §10. The Fit/Transform Boundary

### Motivation

This is the key mechanism separating "research" from "backtest." If the
label basis is refit at every evaluation period, then $c_p$ for different
periods lives in different coordinate systems, and any cross-period
comparison, any feature-vs-alpha comparison, becomes meaningless.

### Specification

During the training/research period, run `fit_label_basis(Y_train)` and
store

$$B_Y;\qquad \Theta_{\mathsf V},\ \Theta_{\mathsf L};\qquad \Gamma_{\mathsf{DV}},\ \Gamma_{\mathsf{DL}},\ \Gamma_{\mathsf{VL}}.$$

When new data $Y_{\rm new}$ arrives:

$$Z_{\rm new}=Y_{\rm new}B_Y,\qquad \Phi_{\mathsf D}^{\rm new}=Z_{\rm new},$$

$$G_{\mathsf V}^{\rm new},G_{\mathsf L}^{\rm new}\ \text{built per §8 and within-date centered},$$

$$R_{\mathsf V}^{\rm new}=G_{\mathsf V}^{\rm new}-\Phi_{\mathsf D}^{\rm new}\Gamma_{\mathsf{DV}},\qquad \Phi_{\mathsf V}^{\rm new}=R_{\mathsf V}^{\rm new}\Theta_{\mathsf V},$$

$$R_{\mathsf L}^{\rm new}=G_{\mathsf L}^{\rm new}-\Phi_{\mathsf D}^{\rm new}\Gamma_{\mathsf{DL}}-\Phi_{\mathsf V}^{\rm new}\Gamma_{\mathsf{VL}},\qquad \Phi_{\mathsf L}^{\rm new}=R_{\mathsf L}^{\rm new}\Theta_{\mathsf L}.$$

Note that within-date centering is **not a fitted parameter** — it's
applied in place on each date, so it does not constitute leakage.

### The Key Property: $\delta_{\rm drift}$

During OOS, $\Phi^{\rm new}$ no longer satisfies $\Phi^\top W\Phi=I$. This
isn't an error — it's a free regime-drift detector:

$$\boxed{\delta_{\rm drift}=\big\|\Phi_{\rm OOS}^\top W_{\rm OOS}\Phi_{\rm OOS}-I\big\|_F,}$$

further split into the three diagonal blocks' $\delta^{\mathsf
D},\delta^{\mathsf V},\delta^{\mathsf L}$.

If $\delta_{\rm drift}$ is small, future returns' covariance structure
matches the research period, and every quantity expressed via $c$ can be
compared directly across periods. If $\delta_{\rm drift}$ is large (say
$>0.3$), market structure has shifted, and a drop in $Q_\alpha$ **could come
from coordinate-system invalidity rather than model failure** — an
alternative explanation that must be ruled out first whenever reading any
OOS decay.

### Worked Example

Research period 2022–2023, OOS is the first half of 2024. If 2024H1's
market volatility is half the research period's, $\|z\|$ is systematically
smaller; since $G_{\mathsf V}\sim\|z\|^2$,

$$\Phi_{\mathsf V}^{\rm OOS\top}W\Phi_{\mathsf V}^{\rm OOS}\approx0.25\,I\ \Longrightarrow\ \delta^{\mathsf V}\approx0.75\sqrt{d_{\mathsf V}}\approx2.4\ (d_{\mathsf V}=10).$$

$Q_\alpha^{\mathsf V}$ will mechanically drop by roughly 4×, which
**cannot be read as decay in the model's activity-predicting power.** The
correct handling is to report $\delta^{\mathsf V}$ alongside it and, when
necessary, make cross-period comparisons on a $\delta$-corrected basis.
