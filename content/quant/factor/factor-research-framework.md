---
title: "Factor Research: Representation, Returns, Risk, and Neutralization"
description: "Factor research as an independent pillar on top of the Universe layer: a notation that fits the existing D/V/T geometry, exact feature-side factor decomposition, signal- vs portfolio-level neutralization, factor-return identification under within-date demeaning, and seven research tracks (IR-23 to IR-29)."
date: 2026-11-08
tags: ["factor", "risk-model", "neutralization", "methodology", "quant-research"]
math: true
---

## The problem this chapter solves

The word *factor* is overloaded. A risk exposure, a rewarded characteristic, a
model feature, a tradable factor portfolio and a latent statistical component
can coincide, but they are not the same object by definition. Research that
lets them blur produces three recurring mistakes:

1. A feature looks predictive, but it is only a proxy for a known common
   exposure.
2. A good ex-post factor $R^2$ is taken as evidence of a good *covariance
   forecast*.
3. "Neutralize every risk factor" is treated as a rule, when it is a mandate
   choice that can delete real premium.

This chapter separates four objects that each need their own estimator,
null, and acceptance criterion:

| Object | Question | Typical failure if conflated |
|---|---|---|
| **Exposure** $F_d$ | which common directions does each asset load on, as known at $d$? | look-ahead through reconstituted taxonomies |
| **Factor return** $\psi_{d,h}$ | what did each common direction earn over horizon $h$? | rank-deficient regressions read as stable returns |
| **Factor risk** $\Sigma^{\psi}_{d,h},\ \Sigma^{u}_{d,h}$ | how do common and specific parts co-vary? | in-sample fit mistaken for forecast quality |
| **Decision role** | explain, earn premium, or constrain? | neutralizing away priced exposure |

Factor structure belongs to the Universe layer as a shared, point-in-time
service: the same $F_d$ is consumed by feature, alpha, covariance and
portfolio research. Its estimation and validation are a separate research
pillar, which this chapter sets out.

{{< callout type="feat" >}}
**Scope.** This is a research agenda with the mathematical contracts it needs.
Every identity stated below is verified numerically (Section 9). The research
tracks (Section 8) are open questions, not results.
{{< /callout >}}

---

## 1. Notation

This chapter extends the notation of
[Cross-Sectional Predictive Geometry §1](/quant/feature-evaluation/cross-sectional-predictive-geometry/#1-notation).
Every symbol there keeps its meaning: $d$ is the decision period, $i$ the
ticker, $h$ the horizon, $W=\operatorname{diag}(w_{id})$ with
$w_{id}=1/(TN_d)$, $\Phi=[\Phi_{\mathsf D}\,|\,\Phi_{\mathsf V}\,|\,\Phi_{\mathsf T}]$
the label-mechanism basis, $c_s=s^\top W\Phi$ a signal's representation,
$Q$ its energy.

Factor objects therefore use letters that are **not** already taken:

| Symbol | Shape | Meaning |
|---|---|---|
| $m$ | scalar | number of factors (not $K$, which is the label dimension) |
| $F_d$ | $N_d\times m$ | factor exposures known at $d$, **within-date centered** (Section 2) |
| $y_{d,h}$ | $N_d$ | within-date demeaned future return at horizon $h$, i.e. column $h$ of $Y^{XS}_d$ |
| $\psi_{d,h}$ | $m$ | factor returns realized over horizon $h$ |
| $u_{d,h}$ | $N_d$ | specific (residual) returns |
| $\Sigma^{\psi}_{d,h}$ | $m\times m$ | conditional factor-return covariance |
| $\Sigma^{u}_{d,h}$ | $N_d\times N_d$ | conditional specific covariance (diagonal only if validated) |
| $\Sigma^{r}_{d,h}$ | $N_d\times N_d$ | conditional asset-return covariance |
| $\Pi_d$ | $N_d\times N_d$ | positive-definite weight for factor projections (Section 3) |
| $\bar F_d$ | $N_d\times(m+1)$ | $[F_d\ \mathbf 1]$: exposures plus the dollar-neutrality direction (Sections 5–6) |
| $P^{\Pi}_{F},\ M^{\Pi}_{F}$ | $N_d\times N_d$ | $\Pi$-orthogonal projector onto $\operatorname{span}(F_d)$ and its complement |
| $\mu_{d,h}$ | $N_d$ | expected-return forecast |
| $\omega_d$ | $N_d$ | portfolio weights |
| $\gamma$ | scalar | risk aversion (as in the [diffusion guide §9](/quant/model/diffusion-statarb-implementation-guide/#9-using-the-distribution)) |
| $\operatorname{erank}(\cdot)$ | scalar | entropy effective rank of a PSD matrix |

Letters deliberately avoided, and why:

| Avoided | Already means | Used here instead |
|---|---|---|
| $B$ for exposures | whitening $B_Y$, bootstrap count $B_\bullet$, band count $B_{\rm band}$ | $F$ (as in Universe §21) |
| $f$ for factor returns | fingerprint $f_s,\tilde f_s$ | $\psi$ |
| $D$ for specific covariance | the $\mathsf D$ (Direction) block | $\Sigma^{u}$ |
| $\Omega$, $\Delta$ for covariances | (kept free; $\Delta$ appears in Universe §2) | $\Sigma^{\psi}$, $\Sigma^{u}$ |
| $K$ for factor count | label dimension | $m$ |
| $C$ for crowding / constraints | predictive operator $C=X^\top W\Phi$ | not used in this chapter |
| $r_{\rm eff}$ for covariance rank | effective spectral rank of $C$ | $\operatorname{erank}$ |
| $w$ for portfolio weights | row weights $w_{id}$ | $\omega$ |
| $\epsilon$ for specific returns | numerical floor; $\varepsilon_\alpha$ residual alpha | $u$ |

---

## 2. The demeaning convention and centered exposures

Every label in this framework is within-date demeaned:
$y_{d,h}=Y^{raw}_{d,\cdot h}-\mathbf 1\bar Y_{d,h}$. That has two consequences
for factor models that must be fixed before any estimation.

### Lemma F.0 (centering is lossless for dollar-neutral portfolios)

Let the raw model be $y^{raw}=F^{raw}\psi+u$ with covariance
$F^{raw}\Sigma^{\psi}F^{raw\top}+\Sigma^{u}$, and let
$C_N=I-\tfrac1N\mathbf 1\mathbf 1^\top$. Define $F=C_NF^{raw}$. For any
portfolio with $\mathbf 1^\top\omega=0$:

$$
\omega^\top y^{raw}=\omega^\top y,\qquad
\omega^\top F^{raw}=\omega^\top F,\qquad
\omega^\top\big(F^{raw}\Sigma^{\psi}F^{raw\top}+\Sigma^{u}\big)\omega
=\omega^\top\big(F\Sigma^{\psi}F^\top+\Sigma^{u}\big)\omega .
$$

*Proof.* $\mathbf 1^\top\omega=0$ gives $C_N\omega=\omega$, and $C_N$ is
symmetric, so $\omega^\top C_N z=\omega^\top z$ for any $z$. $\blacksquare$

So for dollar-neutral research nothing is lost by working with centered
exposures and demeaned labels. Two practical rules follow:

- **The market factor appears as $\beta-\bar\beta$.** A pure intercept column
  becomes identically zero after centering and must be dropped; keeping it
  makes $F^\top\Pi F$ singular.
- **This is exactly why demeaning does not remove beta.** Universe §21 shows
  $r_i-\bar r=(\beta_i-\bar\beta)r_m+(\epsilon_i-\bar\epsilon)$: the
  centered beta column is the part of the market that survives demeaning.

### Lemma F.1 (centered exposures keep residuals centered)

If $\mathbf 1^\top F_d=0$, then for **any** weight $\Pi_d$ and any centered
$y$, $\mathbf 1^\top M^{\Pi}_{F}y=0$.

*Proof.* $\mathbf 1^\top M^{\Pi}_{F}y=\mathbf 1^\top y-\mathbf 1^\top F(\cdots)=0-0$. $\blacksquare$

The alternative condition "$\mathbf 1\in\operatorname{span}(F)$" does
**not** give this once $\Pi\neq I$: it only gives
$\mathbf 1^\top\Pi M^{\Pi}_{F}y=0$, because $M^{\Pi}_F y$ is
$\Pi$-orthogonal to $\operatorname{span}(F)$, not Euclidean-orthogonal. That
is why the geometry uses centered exposures rather than an intercept column.
(The decision side does add $\mathbf 1$ to the basis, deliberately, for a
different reason: dollar neutrality, Section 5.)

---

## 3. Factor projectors and the choice of weight

### Definition

$$
\boxed{
P^{\Pi}_{F}=F_d\big(F_d^\top\Pi_dF_d\big)^{+}F_d^\top\Pi_d,\qquad
M^{\Pi}_{F}=I-P^{\Pi}_{F}.
}
$$

$P^{\Pi}_F$ is idempotent and self-adjoint under
$\langle x,z\rangle_{\Pi}=x^\top\Pi z$, so for any $x,z$:

$$
\langle x,z\rangle_{\Pi}
=\langle P^{\Pi}_Fx,P^{\Pi}_Fz\rangle_{\Pi}
+\langle M^{\Pi}_Fx,M^{\Pi}_Fz\rangle_{\Pi}.
$$

*Proof.* With $G=F^\top\Pi F$: $P^2=FG^{+}GG^{+}F^\top\Pi=P$ because
$G^{+}GG^{+}=G^{+}$; $\Pi P=\Pi FG^{+}F^\top\Pi$ is symmetric; and
$PF=F$ because $G^{+}G$ projects onto the row space of $F$ ($\Pi\succ0$). The
cross term $\langle Px,Mz\rangle_\Pi=x^\top\Pi P(I-P)z=0$. $\blacksquare$

The projector depends only on $\operatorname{span}(F_d)$, so it is unchanged
by rank deficiency or by any invertible re-parameterization $F\mapsto FA$.

This is an **inner-product** identity. With $\Pi=I$ and centered vectors
(Lemma F.1) it is a cross-sectional covariance identity; with $\Pi\neq I$ it
is a weighted one. It does **not** make IC, rank IC, squared geometry scores
or nonlinear importance additive.

### Two weights, two questions

| Weight | Definition | Use |
|---|---|---|
| $\Pi^{\rm geom}_d$ | $I_{N_d}$ | **geometry**: within a date, $W$ is a scalar multiple of $I$, so this projector is $W$-self-adjoint and compatible with $\Phi$. Universe §21's $M_{F_d}$ is this case. |
| $\Pi^{\rm GLS}_d$ | $\big(\Sigma^{u}_{d,h}\big)^{-1}$ | **decision**: the weight implied by a factor risk model in portfolio construction (Proposition F.3). |

The two coincide only when $\Sigma^{u}\propto I$. Reports must say which
weight a "factor-neutral" quantity used.

---

## 4. Feature-side factor decomposition in the mechanism space

Universe §21 handles factors on the **label** side:
$Y^{idio}=M_FY$, after which the whole chain $\Phi^{idio}$ is rebuilt. That
is a valid robustness variant, but $\Phi^{idio}$ is a different coordinate
system: raw and idiosyncratic $Q$ values cannot be subtracted, and blocks from
the two bases cannot be mixed.

The primary definition keeps the raw $\Phi$ fixed and decomposes the
**feature**.

### Proposition F.2 (exact feature-side decomposition)

Let $x$ be a within-date standardized feature with $\|x\|_W=1$, and apply
$P_F=P^{\Pi^{\rm geom}}_{F}$ date by date. Then

$$
\boxed{
c^{\rm full}_x=c^{\rm full}_{P_Fx}+c^{\rm full}_{M_Fx},
\qquad
\|P_Fx\|_W^2+\|M_Fx\|_W^2=1,
}
$$

and the energies satisfy

$$
Q^{\rm full}_x=Q^{\rm full}_{P_Fx}+Q^{\rm full}_{M_Fx}
+2\big\langle c^{\rm full}_{P_Fx},\,c^{\rm full}_{M_Fx}\big\rangle .
$$

The same holds block by block ($\mathsf D,\mathsf V,\mathsf T$) and after
projection by $P_V$.

*Proof.* $c_s=s^\top W\Phi$ is linear in $s$, which gives the first identity.
Within each date $W$ is a scalar multiple of $I$, so $P_F$ is $W$-self-adjoint,
and Pythagoras applied date by date gives the norm identity. Expanding the
squared norm of a sum gives the energy identity. $\blacksquare$

**Reading it.**

- $\rho^2_F(x)=\|P_Fx\|_W^2$ is the share of the feature's cross-sectional
  variance explained by exposures. It is a property of the **feature**, not of
  its predictive power.
- The cross term is generally nonzero. Its sign says whether the factor and
  idiosyncratic parts reinforce or cancel each other along the predictive
  directions. Reporting only $Q_{P_Fx}$ and $Q_{M_Fx}$ hides that cancellation.
- **Missing values.** Where a feature is zero-filled for missing names
  (Predictive Geometry §5.2), project on the valid names only, with $F_d$
  restricted to them and re-centered there. Projecting the full cross-section
  would assign factor-implied values to names that have no observation. The
  restricted projector is still a per-date $W$-self-adjoint projector, so the
  proposition is unchanged.
- To compare idiosyncratic quality on the IC scale, standardize the residual
  exactly as Label–Alpha Geometry §18 does for residual alpha:
  $\tilde x^{M}=M_Fx/\|M_Fx\|_W$ and
  $Q^{\rm full}_{\tilde x^M}=Q^{\rm full}_{M_Fx}/(1-\rho^2_F)$. Report both
  $Q^{\rm full}_{\tilde x^M}$ (quality) and $(1-\rho^2_F)\,Q^{\rm full}_{\tilde x^M}$
  (absolute contribution).

**What it cannot say.** The split is linear and depends on the chosen basis
$F$. A different legitimate factor basis gives a different split
(IR-29). And "explained by exposures" does not mean "priced": that is a
question about $\psi$, not about $F$.

---

## 5. Signal-level versus portfolio-level neutralization

### Setting

Labels are demeaned and exposures centered, so by Lemma F.0 the risk model

$$
\Sigma^{r}=F\Sigma^{\psi}F^\top+\Sigma^{u}
$$

describes the risk of **dollar-neutral** portfolios only. Every optimization
below therefore includes $\mathbf 1^\top\omega=0$. Write

$$
\bar F_d=[\,F_d\ \ \mathbf 1\,],\qquad A=\Pi^{\rm GLS}=(\Sigma^{u})^{-1},\qquad S=(\Sigma^{\psi})^{-1},
$$

with $\Sigma^{u}\succ0$ and $\Sigma^{\psi}\succ0$.

### Proposition F.3

**(a) Hard neutrality is a GLS residualization on $[F,\mathbf 1]$.** The
solution of

$$
\max_{\omega}\ \mu^\top\omega-\frac{\gamma}{2}\,\omega^\top\Sigma^{r}\omega
\quad\text{s.t.}\quad \bar F^\top\omega=0
$$

is

$$
\boxed{
\omega^*=\gamma^{-1}A\,M^{A}_{\bar F}\,\mu,\qquad
M^{A}_{\bar F}=I-\bar F\big(\bar F^\top A\bar F\big)^{+}\bar F^\top A .
}
$$

**(b) The portfolio is invariant to anything in $\operatorname{span}(\bar F)$.**
For any $b$ and any weight $\Pi$, replacing $\mu$ with $\mu+\bar Fb$, with
$M^{\Pi}_{F}\mu$, or with $M^{\Pi}_{\bar F}\mu$ leaves $\omega^*$ unchanged. In
particular the common shift $c\,\mathbf 1$ that demeaned labels cannot
identify (Section 6) has no effect, and neither does any signal-level
pre-neutralization.

**(c) With dollar neutrality only, the factor exposure is**

$$
F^\top\omega^{\rm dn}
=\gamma^{-1}S\big(S+F^\top AF\big)^{-1}F^\top A\,(\mu-\nu\mathbf 1),
\qquad
\nu=\frac{\mathbf 1^\top(\Sigma^{r})^{-1}\mu}{\mathbf 1^\top(\Sigma^{r})^{-1}\mathbf 1},
$$

which is zero if and only if $F^\top A(\mu-\nu\mathbf 1)=0$. It is zero when
$\mu$ is $A$-orthogonal to $\operatorname{span}(\bar F)$ (i.e.
$M^{A}_{\bar F}\mu=\mu$). It is **not** zero in general for an
OLS-neutralized signal ($F^\top\mu=0$, $\mathbf 1^\top\mu=0$), unless
$\Sigma^{u}\propto I$.

*Proof.* (a) On $\bar F^\top\omega=0$ the factor term vanishes, so the
objective is $\mu^\top\omega-\frac\gamma2\omega^\top\Sigma^{u}\omega$. KKT gives
$\gamma\Sigma^{u}\omega=\mu-\bar F\lambda$ with $\bar F^\top\omega=0$, so
$\bar F^\top A\bar F\lambda=\bar F^\top A\mu$ and
$\mu-\bar F\lambda=M^{A}_{\bar F}\mu$ (the product $\bar F\lambda$ is the same
for every solution $\lambda$).
(b) $M^{A}_{\bar F}\bar F=0$, and $M^{\Pi}_F\mu$, $M^{\Pi}_{\bar F}\mu$ differ
from $\mu$ by an element of $\operatorname{span}(\bar F)$.
(c) The dollar-neutral solution is
$\omega^{\rm dn}=\gamma^{-1}(\Sigma^{r})^{-1}(\mu-\nu\mathbf 1)$. Woodbury gives
$F^\top(\Sigma^{r})^{-1}=\big[I-F^\top AF(S+F^\top AF)^{-1}\big]F^\top A=S(S+F^\top AF)^{-1}F^\top A$,
and $S(S+F^\top AF)^{-1}$ is invertible. If $\bar F^\top A\mu=0$, then
$\mathbf 1^\top(\Sigma^{r})^{-1}\mu=\mathbf 1^\top A\mu-\mathbf 1^\top AF(S+F^\top AF)^{-1}F^\top A\mu=0$,
so $\nu=0$ and $F^\top A\mu=0$. If $A=aI$, then $F^\top A\mathbf 1=aF^\top\mathbf 1=0$
and $F^\top A(\mu-\nu\mathbf 1)=aF^\top\mu$, which vanishes for an OLS-neutral
$\mu$. $\blacksquare$

**Consequences.**

- A signal that looks factor-neutral in the geometry ($\Pi^{\rm geom}=I$)
  is not automatically neutral after optimization. The two agree when
  specific variances are equal across names.
- Under a hard constraint the choice of signal-level neutralization does not
  change the portfolio; it only changes the diagnostics. Effort belongs in the
  constraint and in $\Sigma^{u}$.
- Under a soft treatment (dollar-neutral optimization with a finite factor
  covariance, a quadratic exposure penalty, which acts like a larger
  $\Sigma^{\psi}$), the portfolio is factor-neutral whenever the signal is
  $A$-orthogonal to $[F,\mathbf 1]$; that is sufficient, and the exact
  condition is (c). An OLS-neutral signal generally does not satisfy it.
  Signal-side neutralization must use that weight and that basis for the
  signal-side and portfolio-side definitions to agree. (A signal that
  satisfies it also leaves any exposure budget non-binding.)

---

## 6. Factor returns and their identification

### Cross-sectional estimator

For each $(d,h)$, regress on the augmented basis $\bar F_d=[F_d\ \mathbf 1]$
and keep the first $m$ coefficients:

$$
\begin{bmatrix}\hat\psi_{d,h}\\ \hat a_{d,h}\end{bmatrix}
=\big(\bar F_d^\top\Pi_d\bar F_d\big)^{+}\bar F_d^\top\Pi_d\,y_{d,h}.
$$

Why the intercept column, when the exposures are already centered:

- With $\Pi=I$ it makes no difference. By Frisch–Waugh–Lovell, regressing
  demeaned $y$ on centered $F$ gives the same $\hat\psi$ as regressing raw $y$
  on $[F^{raw},\mathbf 1]$.
- With $\Pi\neq I$ (for example $\Pi^{\rm GLS}$) it does. Partialling out the
  intercept under weight $\Pi$ removes the $\Pi$-weighted mean, not the equal
  mean, so the centered regression is no longer the GLS estimator of the raw
  model. With $\mathbf 1$ in the basis, demeaning $y$ only changes $\hat a$, so
  $\hat\psi$ is the same for raw and demeaned labels.

Row $k\le m$ of $R_d=(\bar F_d^\top\Pi_d\bar F_d)^{-1}\bar F_d^\top\Pi_d$ is a
**pure factor portfolio**: since $R_d\bar F_d=I_{m+1}$, it has unit exposure to
factor $k$, zero exposure to every other factor, and zero net weight
($R_{d,k}\mathbf 1=0$). $\hat\psi_{d,h,k}$ is that portfolio's return.
Time-series summaries of $\hat\psi$ under the chosen $\Pi$ are the
Fama–MacBeth estimates.

### Lemma F.4 (centered group dummies lose one rank)

Let $H_d\in\{0,1\}^{N_d\times G}$ be **non-overlapping** group dummies, so
$H_d\mathbf 1_G=\mathbf 1_{N_d}$. Then $C_NH_d\mathbf 1_G=0$, so the
centered dummies have rank at most $G-1$.

Identification needs one linear restriction, chosen and reported in advance,
for example:

- drop one group (returns are then relative to the dropped group), or
- impose $\sum_g s_{g,d}\,\psi_{g,d,h}=0$ with group weights $s_{g,d}$
  (returns are then relative to the weighted average group).

The restriction changes how individual group returns are read, but not the
fitted values $F_d\hat\psi_{d,h}$, and not any projector $P^{\Pi}_F$ (which
depends only on the span). Note that a pseudo-inverse silently picks the
minimum-norm solution, which for centered dummies is the equal-weight
restriction $\sum_g\psi_g=0$; if a different restriction is intended, impose
it explicitly. Overlapping or soft memberships (IR-06) do
not satisfy $H\mathbf 1=\mathbf 1$; their rank must be checked numerically
rather than assumed.

### Expected return: factor part and specific part (IR-21)

$$
\mu_{d,h}=F_d\,\bar\psi_{d,h}+\bar u_{d,h},
$$

where $\bar\psi_{d,h}$ is a **forecast** of factor returns and $\bar u_{d,h}$
a **forecast** of specific returns. Two cautions:

- **Identification.** Every model trained on demeaned labels identifies
  $\mu_{d,h}$ only up to a common shift $c\,\mathbf 1$. Under
  $\mathbf 1^\top\omega=0$, which the centered risk model requires anyway
  (Section 5), the shift has no effect (Proposition F.3b). A portfolio that is
  not dollar-neutral needs the common component from a separate model, and a
  risk model with the market intercept.
- **Ex-post is not ex-ante.** Regressing a model's score on $F_d$ after the
  fact gives a decomposition of the score, not two independent forecasts.
  Calling the projected part "factor timing" double counts whenever the factor
  premium is also forecast elsewhere.

---

## 7. Covariance and attribution identities

### Risk model (IR-32)

$$
\Sigma^{r}_{d,h}=F_d\,\Sigma^{\psi}_{d,h}\,F_d^\top+\Sigma^{u}_{d,h},
$$

valid under three assumptions that each need a test:

| Assumption | Test |
|---|---|
| A1: $F_d$ is known at $d$ | point-in-time audit of every exposure input |
| A2: $\operatorname{Cov}(\psi_{d,h},u_{d,h})=0$ | correlation of $\hat\psi$ with cross-sectional residual moments |
| A3: the chosen structure of $\Sigma^{u}$ (diagonal, block, sparse) | residual correlation scans against a matched null (IR-27) |

$\Sigma^{\psi}$, $\Sigma^{u}$ and the label horizon $h$ must refer to the same
horizon. A 4h mean paired with a daily covariance is a contract violation, not
a modeling choice.

### Portfolio attribution (IR-44)

With exposures fixed at $d$, returns realized over $h$, and a dollar-neutral
portfolio ($\mathbf 1^\top\omega_d=0$):

$$
\omega_d^\top y_{d,h}
=\underbrace{\big(F_d^\top\omega_d\big)^\top\psi_{d,h}}_{\text{factor}}
+\underbrace{\omega_d^\top u_{d,h}}_{\text{specific}},
$$

$$
\operatorname{Var}_d\big(\omega_d^\top y_{d,h}\big)
=\omega_d^\top F_d\Sigma^{\psi}_{d,h}F_d^\top\omega_d
+\omega_d^\top\Sigma^{u}_{d,h}\omega_d
\quad\text{under A2.}
$$

The return identity holds by construction for any $\hat\psi$ once
$\hat u=y-F\hat\psi$ is defined as the residual. The variance split holds
only under A2, and only for dollar-neutral portfolios: the demeaned specific
return is $C_Nu$, and $\omega^\top C_N\Sigma^{u}C_N\omega=\omega^\top\Sigma^{u}\omega$
needs $C_N\omega=\omega$ (Lemma F.0). A large realized cross term is a
diagnostic signal, not noise to be allocated.

---

## 8. Research tracks

Each track follows the research-note template: a falsifiable question, the
objects it estimates, baselines and nulls, evaluation, and a versioned
deliverable. Priorities are dependency order, not expected PnL. Each track is
tagged with the **role** it serves: *explain* (common variation), *premium*
(forecastable return), or *control* (constraints and risk budgets).

### IR-23 · Factor ontology and taxonomy — P1 · explain / premium / control

- **Question.** Which factors are market, asset-class, sector, country,
  currency, style or latent, and which role is each meant to serve?
- **Methods.** Curated economic exposures, point-in-time classifications,
  nested taxonomy, an explicit role label per factor.
- **Evaluation.** Interpretability, identifiability (Lemma F.4), exposure
  coverage, applicability across equities, futures, FX and crypto.
- **Deliverable.** `factor_catalog` with type, horizon, definition, estimation
  method and intended role.

### IR-24 · Factor discovery, number and selection — P1 · explain

- **Question.** What is a useful factor basis, and how many dimensions are
  actually estimable?
- **Methods.** Characteristic models, rolling and weighted PCA, dynamic
  factors, low-rank-plus-sparse, latent factors; random-matrix (Marchenko–Pastur)
  baselines for the number of factors.
- **Evaluation.** Out-of-sample explained variance **and** covariance
  forecast quality; $\operatorname{erank}$, subspace stability, residual
  dependence, parameter uncertainty.
- **Deliverable.** Candidate-basis comparison and a factor-selection protocol.
- **Caution.** Statistical factors have arbitrary sign and rotation. Compare
  subspaces with principal angles under $\langle\cdot,\cdot\rangle_{\Pi}$, never
  column by column. More factors can improve in-sample fit without improving
  risk.

### IR-25 · Exposure estimation — P1 · explain / control

- **Question.** Should exposures be characteristics, trailing return betas,
  conditional states, or outputs of a hierarchical model?
- **Methods.** Rolling and weighted regressions, shrinkage, state-space betas
  and Kalman filters, characteristic exposures. All columns centered within
  date (Section 2).
- **Evaluation.** Exposure stability, out-of-sample return decomposition,
  exposure uncertainty, behaviour under missing assets.
- **Deliverable.** Timestamped $F_d$ with exposure confidence estimates.

### IR-26 · Factor returns and premia — P1 · premium

- **Question.** Can cross-sectional factor returns be estimated reliably, and
  are expected premia predictable?
- **Methods.** Constrained cross-sectional regressions (Section 6),
  Fama–MacBeth summaries, rolling conditional forecasts $\bar\psi_{d,h}$.
- **Evaluation.** Identification (rank and restriction reported),
  stability of $\hat\psi$, conditional calibration of $\bar\psi$.
- **Deliverable.** Historical $\hat\psi$ and a **separate** forecast
  $\bar\psi$; the two are never the same artifact.

### IR-27 · Specific risk and omitted factors — P1 · explain / control

- **Question.** Is $\Sigma^{u}$ diagonal, sparse, clustered or heavy-tailed?
  What common structure was omitted?
- **Methods.** Residual correlation scans; residual eigenvalues against the
  Marchenko–Pastur edge as candidate omitted factors; conditional residual
  covariance; sparse corrections. Test cross-sectional residual dependence and
  time-series residual autocorrelation **separately**: they are different
  failures.
- **Evaluation.** Out-of-sample specific-risk calibration and residual
  dependence, not explanatory $R^2$.
- **Deliverable.** Validated $\Sigma^{u}$ and a failure-mode report.

### IR-28 · Neutralization and constraint choice — P1 · control

- **Question.** Which exposures should be free, penalized, budgeted or
  hard-neutralized?
- **Methods.** Compare raw, sector-neutral, style-neutral and hierarchical
  partial-neutral versions, holding $\mu$ and $\Sigma^{r}$ fixed; use
  Proposition F.3 to separate signal-side from portfolio-side effects.
- **Evaluation.** Retained forecast information, feasible-set reduction, risk
  concentration and decision value.
- **Deliverable.** A factor policy report stating, per factor, the role
  (IR-23), the treatment, and the weight $\Pi$ used. No blanket rule that all
  risk factors must be neutralized: a factor with a forecastable premium is
  also alpha.

### IR-29 · Factor model uncertainty — P2 · control

- **Question.** How sensitive are portfolio decisions to alternative
  legitimate factor bases and noisy exposures?
- **Methods.** Competing risk-model ensembles, exposure perturbation,
  Bayesian or robust alternatives, alignment across universes.
- **Evaluation.** Change in weights and risk under basis changes; the spread
  of Proposition F.2's split across bases; calibration and stress sensitivity.
- **Deliverable.** A factor-model uncertainty envelope.

### Interfaces with other layers

| Layer | Track | What it consumes from this chapter | Contract | Where it is applied |
|---|---|---|---|---|
| Universe | IR-06, IR-08 | $F_d$, group dummies $H_d$ | point-in-time membership; rank checked (Lemma F.4) | [Universe §21](/quant/universe/multi-sub-universe-directional-coupling-geometry/#21-beta-demeaning-and-the-factor-neutral-variant): centered exposures, membership factor share, exposure-matched null; [§13.5](/quant/universe/multi-sub-universe-directional-coupling-geometry/#135-within-group-family) within-group candidates |
| Feature | IR-11, IR-12 | $P^{\Pi^{\rm geom}}_F$ | Proposition F.2 on fixed raw $\Phi$; label-side variant only as robustness | [Predictive Geometry §35](/quant/feature-evaluation/cross-sectional-predictive-geometry/#35-factor-decomposition-of-a-feature); [Selection Bridge](/quant/feature-evaluation/selection-bridge-eligibility-grouping-coverage/#factor-aware-eligibility) mandate-dependent eligibility, redundancy and Score; [Selection guide B5b](/quant/feature-evaluation/cross-sectional-feature-selection/#b5b-factor-split-of-the-operator-when-a-factor-model-is-available) |
| Alpha | IR-21 | $F_d$, $\bar\psi_{d,h}$, $\Sigma^u$ | demeaned $\mu$ identified up to a shift; ex-post split is not a forecast | [Label–Alpha §28](/quant/alpha-evaluation/label-alpha-geometry/#28-factor-aligned-and-specific-alpha): factor/specific split, four-cell decomposition, decision alpha |
| Covariance | IR-32 | $\Sigma^{\psi}$, $\Sigma^{u}$ | same horizon as $\mu$; PSD | planned | planned |
| Portfolio | IR-28, IR-39 | $F_d$, $\Pi^{\rm GLS}$ | Proposition F.3; $\mathbf 1^\top\omega=0$ when labels are demeaned | planned | planned |
| Attribution | IR-44 | $F_d$, $\hat\psi_{d,h}$ | exposures lagged to $d$; variance split only under A2 | planned | planned |

---

## 9. Numerical identity tests

In the same spirit as
[Label–Alpha Geometry §25](/quant/alpha-evaluation/label-alpha-geometry/#25-numerical-validation),
each identity in this chapter has a unit test. All must pass to floating-point
tolerance on random inputs before any empirical result is reported.

| Test | Check |
|---|---|
| Centering (Lemma F.0) | $\omega^\top F^{raw}=\omega^\top F$ and equal quadratic risk for random $\omega$ with $\mathbf 1^\top\omega=0$; intercept column becomes zero |
| Residual centering (Lemma F.1) | $\mathbf 1^\top M^{\Pi}_Fy=0$ for $\Pi=I$ and $\Pi=(\Sigma^u)^{-1}$ |
| Projector | $P^2=P$ and $\Pi P=(\Pi P)^\top$ |
| Feature decomposition (Prop. F.2) | $c_x=c_{P_Fx}+c_{M_Fx}$; norms sum to 1; energy identity including the cross term |
| Hard neutrality (Prop. F.3a) | KKT solution with $[F,\mathbf 1]^\top\omega=0$ equals $\gamma^{-1}AM^{A}_{\bar F}\mu$, and is dollar-neutral |
| Invariance (Prop. F.3b) | KKT solution unchanged under $\mu\to\mu+c\mathbf 1$, $\mu\to M^{I}_F\mu$, $\mu\to M^{I}_{\bar F}\mu$ |
| Dollar-neutral exposure (Prop. F.3c) | closed form matches $F^\top\omega^{\rm dn}$; zero for $\mu=M^{A}_{\bar F}\mu$; nonzero for OLS-neutral $\mu$ when $\Sigma^u\not\propto I$; zero for OLS-neutral $\mu$ when $\Sigma^u\propto I$ |
| Intercept handling (Section 6) | $\Pi=I$: centered regression equals regression with intercept; $\Pi\neq I$: they differ, and the augmented regression is invariant to demeaning $y$ |
| Missing values (Prop. F.2) | valid-only projection leaves zero-filled names at zero and keeps both pieces centered and orthogonal |
| Group rank (Lemma F.4) | centered non-overlapping dummies have rank $G-1$ |
| Pure factor portfolios | $R_d\bar F_d=I_{m+1}$ and $R_{d,k}\mathbf 1=0$ for $k\le m$ |
| Pseudo-inverse restriction | min-norm sector returns from centered dummies sum to zero |
| Rotation invariance | $P^{\Pi}_F$, $\rho^2_F$ and Prop. F.2's split unchanged under $F\mapsto FA$ for invertible $A$ |

---

## 10. What this chapter can and cannot answer

**Can answer**

- How much of a feature's cross-sectional variation, and of its predictive
  representation, lies along known exposures (Proposition F.2).
- Whether a "neutral" signal stays neutral after optimization, and when
  signal-side neutralization is irrelevant (Proposition F.3).
- How factor returns are identified under demeaned labels and group dummies.
- Which factor-model assumptions each downstream number depends on.

**Cannot answer**

- Whether a factor is priced in general. Proposition F.2 measures alignment
  with exposures, not premia; premia are IR-26.
- Which factor basis is "true". Every split here is relative to a chosen
  $F_d$ (IR-29).
- Anything about trading costs, capacity or crowding; those are separate
  layers.
