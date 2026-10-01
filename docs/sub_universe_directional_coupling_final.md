# Multi-Sub-Universe Directional Coupling Geometry

### 一套用於 Sub-Universe 發現、Directional Representation 壓縮與後續 Feature–Alpha 研究的統一幾何方法

---

## 本章要解決的問題

既有的 Cross-Sectional Predictive Geometry 從固定 universe 出發：先由 future return 建構 label mechanism basis

$$
\Phi=[\Phi_{\mathsf D}\mid\Phi_{\mathsf V}\mid\Phi_{\mathsf L}],
\qquad
\Phi^\top W\Phi=I,
$$

再把 feature、alpha 與 residual 投影到同一個座標系中。這個設定回答「在給定 universe 裡，哪些 feature 有預測力、模型學到了什麼」，但沒有回答更前置的問題：

> **應該在哪一群 ticker 上建立 feature 與模型？是否存在多個 sub-universe，使 directional label structure 比母 universe 更集中、更低維、且和 activity structure 形成穩定的低秩耦合？**

本章把 universe 本身提升成一個可分析的幾何物件。對任意 candidate sub-universe $\mathcal S_g$，不把 membership indicator 當成普通 feature，而是把它視為對 observation measure 的 restriction。由此定義

$$
G_g=\Phi_{\mathcal U}^\top W_g\Phi_{\mathcal U},
$$

其中 $\Phi_{\mathcal U}$ 是在母 universe $\mathcal U$ 上擬合的共同 reference basis，$W_g$ 是只在 $\mathcal S_g$ 上重新正規化的 row weight。

$G_g$ 是 reference label basis 在 sub-universe 上的 Gram operator。它回答：母 universe 中彼此正交、等尺度的 label mechanisms，進入 $\mathcal S_g$ 後哪些被放大、哪些被抑制、哪些開始互相耦合。

本章的核心進一步不是只看 directional energy $G_{g,\mathsf{DD}}$，而是研究 activity 與 direction 的 canonical coupling：

$$
R_{g,\mathsf{VD}}
=
G_{g,\mathsf{VV}}^{-1/2}
G_{g,\mathsf{VD}}
G_{g,\mathsf{DD}}^{-1/2},
$$

以及 directional explained-variance operator

$$
\boxed{
H^{\mathsf D\leftarrow\mathsf V}_g
=
R_{g,\mathsf{VD}}^\top R_{g,\mathsf{VD}}
=
G_{g,\mathsf{DD}}^{-1/2}
G_{g,\mathsf{DV}}
G_{g,\mathsf{VV}}^{-1}
G_{g,\mathsf{VD}}
G_{g,\mathsf{DD}}^{-1/2}.
}
$$

它的 eigenvalues 為 canonical correlations 的平方；eigenvectors 則是「在該 sub-universe 中最能由 activity geometry 表徵的 directional modes」。

因此本章最終不輸出唯一最佳 universe，而輸出一族可重疊的 universes

$$
\boxed{
\mathscr S
=
\{\mathcal S_1,\dots,\mathcal S_G\},
}
$$

並為每個 universe 報告四個彼此分離的面向：

1. **Directional enrichment**：$\mathsf D$ energy 是否高於母 universe；
2. **Directional compressibility**：$\mathsf D$ energy 是否集中在少數 directional modes；
3. **$\mathsf V\leftrightarrow\mathsf D$ coupling**：activity 與 direction 是否存在強、低秩、可重複的 canonical channels；
4. **Out-of-sample stability**：上述 geometry 是否在未參與 discovery 的未來資料仍成立。

這四項都只是 **label-side structure**。真正的可學習性仍需在 universe 被 freeze 後，由 point-in-time feature 與 alpha 的 OOS 結果確認。

---

# 第九部分：Multi-Sub-Universe Directional Geometry

## §68. 研究目標與輸出物件

### 動機

令母 universe 為 $\mathcal U_d$，允許隨日期 $d$ 改變。傳統做法通常先以產業、市值、流動性或人工規則切 universe，再在每個 universe 裡建立 feature。這種流程把 universe 當作外生設定，無法回答：

- 哪些 ticker 的 directional label geometry 本身比較集中；
- 哪些 sub-universe 中 activity 與 direction 形成低維耦合；
- 是否存在多個彼此不同、但各自有清楚 directional representation 的 market domains；
- 一個 universe 的 label-side 結構是否只是 beta dispersion、volatility scale 或 sample mining 的副產品。

因此 universe discovery 被明確定義成 feature engineering 之前的一個獨立研究層。

### 定義 68.1（母 universe）

對每個 observation time $d=1,\dots,T$，母 universe 為

$$
\mathcal U_d,
\qquad
N_d=|\mathcal U_d|.
$$

所有 candidate universes 都必須滿足

$$
\mathcal S_{g,d}\subseteq\mathcal U_d.
$$

### 定義 68.2（Universe family）

允許同時保留 $G$ 個 sub-universes：

$$
\boxed{
\mathscr S
=
\{\mathcal S_1,\dots,\mathcal S_G\}.
}
$$

不要求彼此互斥：

$$
\mathcal S_g\cap\mathcal S_h\neq\varnothing
$$

是合法的；也允許 nested family：

$$
\mathcal S_1\subset\mathcal S_2\subset\cdots.
$$

本章的目標不是求解

$$
\arg\max_{\mathcal S}J(\mathcal S)
$$

得到唯一 winner，而是從預先限定的 candidate generator 中保留一組通過統計與穩定度檢驗、且在 directional mechanism 上具有非冗餘結構的 universes。

### 三個層次

本章把「好的 sub-universe」拆成三個不同命題：

**U1：D-enriched。**

$$
q_g^{\mathsf D}>1.
$$

該 subset 每個 directional dimension 的平均 label energy 高於母 universe。

**U2：D-compressible。**

$$
r_{g,\rm eff}^{\mathsf D}\ll d_{\mathsf D}.
$$

Directional energy 集中在少數 horizon combinations。

**U3：D-structurally represented。**

$$
Q_g^{\mathsf D\leftarrow\mathsf V}
\text{ large},
\qquad
r_{g,\rm eff}^{\mathsf D\leftarrow\mathsf V}
\text{ small}.
$$

即少數 activity modes 與少數 directional modes 形成強 canonical coupling。

三者互不等價，報表必須分開。

---

## §69. Point-in-Time 與 Discovery / Evaluation 邊界

### 原則

Universe discovery 使用 future label，因此它本身是 supervised selection。若在同一期間同時：

1. 用 label 選 membership；
2. 用同一 label 評估該 membership 的 geometry；
3. 把結果當作統計證據；

則必然產生 selection bias。

### 定義 69.1（時間切分）

將 universe research 的資料切成：

$$
\mathcal D_U
\quad\text{(discovery)},
$$

$$
\mathcal E_U
\quad\text{(evaluation)},
$$

兩者中間留 purge gap

$$
\Delta_U\ge H_{\max}.
$$

若 candidate generator 使用長度為 $\ell_U$ 的 rolling history，則

$$
\boxed{
\Delta_U\ge H_{\max}+\ell_U^{\rm overlap},
}
$$

其中 $\ell_U^{\rm overlap}$ 是任何會跨越 evaluation 邊界的額外 lookback overlap。

### 規則 69.1（Freeze）

在 $\mathcal D_U$ 上決定：

- candidate generation rule；
- membership $\mathcal S_g$ 或 membership rule；
- reference basis $\Phi_{\mathcal U}$；
- retained canonical rank；
- universe family pruning rule；
- 所有 threshold。

進入 $\mathcal E_U$ 後全部 freeze。不得因 evaluation 結果重新調整 membership 或 thresholds。

### Static 與 dynamic universe

**Static research universe**：membership 只由 discovery period 決定，evaluation 全期固定。

**Dynamic research universe**：對每個 evaluation time $t$，membership 只能由 $t$ 之前已實現的 rolling geometry 決定：

$$
\mathcal S_{g,t}
=
\mathcal A_g(\mathcal F_{t^-}),
$$

其中 $\mathcal F_{t^-}$ 表示 $t$ 時點以前可得資訊。

任何使用 $Y_{t\to t+h}$ 直接決定 $\mathcal S_{g,t}$ 的規則都構成 leakage。

---

## §70. Full-Universe Reference Basis

### 動機

若每個 candidate universe 一開始就各自重新 fit whitening 與 $\Phi_g$，不同 universe 的座標系會不同，$q_D$、coupling direction 與 eigenvector 無法直接比較。因此 discovery 階段先建立一個共同 reference coordinate system。

### 定義

在 discovery period 的母 universe 上，完全依既有 framework 擬合：

$$
Y^{XS}_{\mathcal U}
\longrightarrow
Z_{\mathcal U}
\longrightarrow
\Phi_{\mathcal U}
=
[\Phi_{\mathsf D}\mid\Phi_{\mathsf V}\mid\Phi_{\mathsf L}].
$$

並滿足

$$
\boxed{
\Phi_{\mathcal U}^\top
W_{\mathcal U}
\Phi_{\mathcal U}
=I_{d_\Phi}.
}
$$

其中

$$
W_{\mathcal U}
=
\operatorname{diag}
\left(
\frac1{TN_d}
\right).
$$

### 解讀

這條恆等式把母 universe 定成幾何基準：

- 每個 reference mechanism direction 的 energy 為 1；
- 不同 reference directions 彼此 orthogonal；
- 三個 block 在母 universe 上互相 orthogonal。

因此任何 candidate subset 上的偏離都可以直接解讀為 restriction 造成的 concentration、suppression 或 coupling。

### Fit / transform

$B_Y$、$\Theta_{\mathsf V}$、$\Theta_{\mathsf L}$、$\Gamma_\bullet$ 全部只在 $\mathcal D_U$ fit。Evaluation period 的 label 只能 transform 到已凍結的 reference basis，不重新 fit。

---

## §71. Sub-Universe Weight 與 Restriction Measure

### §71.1 定義

對 candidate universe $\mathcal S_g$，令

$$
N_{g,d}=|\mathcal S_{g,d}|.
$$

定義重新正規化的 row weight：

$$
\boxed{
w^{(g)}_{id}
=
\frac{\mathbf1\{i\in\mathcal S_{g,d}\}}
{T_gN_{g,d}},
}
$$

其中 $T_g$ 是該 universe 有效的 observation times 數量。

實作上不建立 diagonal matrix，所有量以

$$
\frac1{T_g}
\sum_{d\in\mathcal D_g}
\frac1{N_{g,d}}
(\cdot)
$$

計算。

### 為什麼要重新 normalize

若只做 mask：

$$
\widetilde W_g=M_gW_{\mathcal U},
$$

則所有 energy 會隨 subset size 線性縮小；這混入了「有多少 ticker」與「subset 本身的 label geometry」兩件不同的事。

本章使用 $W_g$ 的目的，是問：

> **若把 $\mathcal S_g$ 本身視為一個完整 cross-section，它的平均 label geometry 相對於母 universe 如何改變？**

因此每個 date 仍等權，每個 subset 內 ticker 等權。

### 最低規模

若 $N_{g,d}<d_\Phi+5$，Gram matrices 可能數值上可算但統計上極不穩定。本章要求至少：

$$
\operatorname{median}_dN_{g,d}\ge N_{\min},
$$

其中 $N_{\min}$ 應遠大於 $d_\Phi$；預設見 §88。

---

## §72. Restriction Gram Operator

### 定義 72.1

對任意 candidate universe $g$：

$$
\boxed{
G_g
=
\Phi_{\mathcal U}^\top
W_g
\Phi_{\mathcal U}
\in\mathbb R^{d_\Phi\times d_\Phi}.
}
$$

### 命題 72.1（正半定）

$$
G_g\succeq0.
$$

**證明。** 對任意 $v\in\mathbb R^{d_\Phi}$，

$$
v^\top G_gv
=
(\Phi_{\mathcal U}v)^\top
W_g
(\Phi_{\mathcal U}v)
\ge0.
$$

∎

### §72.1 Rayleigh quotient 的物理意義

取任意 unit mechanism direction $\|v\|=1$。在母 universe：

$$
\|\Phi_{\mathcal U}v\|^2_{W_{\mathcal U}}
=v^\top Iv=1.
$$

在 subset：

$$
\boxed{
\|\Phi_{\mathcal U}v\|^2_{W_g}
=v^\top G_gv.
}
$$

因此 $v^\top G_gv$ 就是該 full-universe mechanism combination 在 $\mathcal S_g$ 中的 **relative energy**。

若

$$
G_gv_k=\lambda_{g,k}v_k,
$$

則：

- $\lambda_{g,k}>1$：該 direction 在 subset 中相對母 universe 被濃縮；
- $\lambda_{g,k}<1$：被抑制；
- $\lambda_{g,k}=1$：平均 energy 不變。

### §72.2 不是 learnability operator

$G_g$ 完全由 future labels 構成。大 eigenvalue 只代表 label variation / mechanism energy 集中，不代表存在任何 point-in-time $X$ 可以預測該 direction。

因此 $G_g$ 是 **label-geometry operator**，不是 predictive operator。

---

## §73. Block Decomposition 與 Directional Enrichment

依 reference basis 的 block 結構：

$$
G_g=
\begin{pmatrix}
G_{g,\mathsf{DD}} & G_{g,\mathsf{DV}} & G_{g,\mathsf{DL}}\\
G_{g,\mathsf{VD}} & G_{g,\mathsf{VV}} & G_{g,\mathsf{VL}}\\
G_{g,\mathsf{LD}} & G_{g,\mathsf{LV}} & G_{g,\mathsf{LL}}
\end{pmatrix}.
$$

### §73.1 每維平均 mechanism energy

由於 block dimensions 不同，不能直接比較 trace。定義

$$
\boxed{
q_g^{\mathsf B}
=
\frac1{d_{\mathsf B}}
\operatorname{tr}(G_{g,\mathsf{BB}}),
\qquad
\mathsf B\in\{\mathsf D,\mathsf V,\mathsf L\}.
}
$$

母 universe 精確滿足：

$$
q_{\mathcal U}^{\mathsf D}
=q_{\mathcal U}^{\mathsf V}
=q_{\mathcal U}^{\mathsf L}
=1.
$$

所以 $q_g^{\mathsf D}=1.7$ 的意思就是：

> 該 sub-universe 每個 reference directional dimension 的平均 squared amplitude 是母 universe 的 $1.7\times$。

### §73.2 Universe mechanism fingerprint

定義

$$
\boxed{
f_g^{\mathsf B}
=
\frac{q_g^{\mathsf B}}
{q_g^{\mathsf D}+q_g^{\mathsf V}+q_g^{\mathsf L}}.
}
$$

母 universe baseline 為

$$
f_{\mathcal U}=(1/3,1/3,1/3).
$$

$q_g^{\mathsf D}$ 是 absolute enrichment；$f_g^{\mathsf D}$ 是 relative purity。兩者必須同時報告。

### 實例

Universe A：

$$
(q_D,q_V,q_L)=(0.8,0.15,0.15).
$$

則 $f_D\approx0.73$，但 $q_D<1$。它是 direction-pure but weak。

Universe B：

$$
(q_D,q_V,q_L)=(1.9,0.9,0.8).
$$

同時具有 absolute D enrichment 與 relative D enrichment。

---

## §74. Directional Spectrum 與 Compressibility

### 動機

$q_g^{\mathsf D}$ 大並不代表 directional representation 簡單。若所有 $d_{\mathsf D}$ 個 directions 一起等比例放大，則沒有低維 structure。

### §74.1 Directional eigenspectrum

對

$$
G_{g,\mathsf{DD}}
=V_g^{\mathsf D}
\Lambda_g^{\mathsf D}
V_g^{\mathsf D\top}
$$

做 eigendecomposition：

$$
\lambda_{g,1}^{\mathsf D}
\ge\cdots\ge
\lambda_{g,d_{\mathsf D}}^{\mathsf D}\ge0.
$$

$v_{g,k}^{\mathsf D}$ 是在該 universe 中第 $k$ 個最被濃縮的 reference directional combination。

### §74.2 Effective directional rank

令

$$
p_{g,k}^{\mathsf D}
=
\frac{\lambda_{g,k}^{\mathsf D}}
{\sum_j\lambda_{g,j}^{\mathsf D}}.
$$

定義 entropy effective rank：

$$
\boxed{
r_{g,\rm eff}^{\mathsf D}
=
\exp\left(
-\sum_kp_{g,k}^{\mathsf D}
\log p_{g,k}^{\mathsf D}
\right).
}
$$

$$
1\le r_{g,\rm eff}^{\mathsf D}\le d_{\mathsf D}.
$$

### §74.3 Top-$m$ directional concentration

$$
\boxed{
C_{g,m}^{\mathsf D}
=
\frac{\sum_{k=1}^m\lambda_{g,k}^{\mathsf D}}
{\sum_j\lambda_{g,j}^{\mathsf D}}.
}
$$

並定義

$$
\boxed{
r_{g,\eta}^{\mathsf D}
=
\min\left\{m:C_{g,m}^{\mathsf D}\ge\eta\right\}.
}
$$

### 解讀

- 高 $q_D$ + 高 $r_{\rm eff}^{D}$：方向能量強但分散；
- 高 $q_D$ + 低 $r_{\rm eff}^{D}$：方向能量強且低維；
- 低 $q_D$ + 低 $r_{\rm eff}^{D}$：低維但弱；
- 低 $q_D$ + 高 $r_{\rm eff}^{D}$：弱且分散。

這一節只回答 **D representation 是否集中**，還沒有回答 magnitude / activity 是否提供對該 D structure 的表徵。

---

## §75. 從 $|Z|^\top WZ$ 到 $\mathsf V\leftrightarrow\mathsf D$ Coupling

### 動機

一個最直觀的 magnitude-direction cross-moment 是

$$
A_{MD}=|Z|^\top WZ.
$$

其元素

$$
(A_{MD})_{hk}
=\mathbb E_W[|Z_h|Z_k]
$$

回答：horizon $h$ 的 move magnitude 與 horizon $k$ 的 signed direction 是否相關。

這個物件的直覺正確，但它有三個限制：

1. $|Z|$ 只保留逐 horizon absolute magnitude，沒有完整二階 activity structure；
2. singular values 同時受 magnitude variance 與 direction variance 影響；
3. 既有 framework 已經有一個更完整且與 $\mathsf D$ 正交化過的 activity block $\Phi_{\mathsf V}$。

因此正式方法以 $\Phi_{\mathsf V}$ 取代 $|Z|$，把問題寫成 subset 中的 $\mathsf V\leftrightarrow\mathsf D$ coupling。

### §75.1 Cross-block Gram

定義

$$
\boxed{
G_{g,\mathsf{VD}}
=
\Phi_{\mathsf V}^\top
W_g
\Phi_{\mathsf D}.
}
$$

在母 universe：

$$
\boxed{
G_{\mathcal U,\mathsf{VD}}=0
}
$$

由 reference basis 的 block orthogonality 精確成立。

在 subset 上一般不再為零：

$$
G_{g,\mathsf{VD}}\neq0.
$$

其物理意義是：

> 母 universe 中已被正交分離的 activity 與 direction，在 $\mathcal S_g$ 這個 restriction measure 下出現系統性共變。

這不是 basis construction error；它是 subset geometry 的訊號。

---

## §76. Canonical $\mathsf V\leftrightarrow\mathsf D$ Coupling Operator

### 動機

直接看 $G_{g,\mathsf{VD}}$ 仍會混入 subset 中 $\mathsf V$ 和 $\mathsf D$ 自身的 scale。需要把兩邊各自 whiten。

### §76.1 Regularized inverse square root

對

$$
G_{g,\mathsf{VV}}
=E_V\Lambda_VE_V^\top,
$$

$$
G_{g,\mathsf{DD}}
=E_D\Lambda_DE_D^\top,
$$

只保留

$$
\lambda_j\ge\epsilon_{\rm can}\lambda_{\max}.
$$

定義 pseudo-inverse square root：

$$
G_{g,\mathsf{VV}}^{-1/2}
=E_{V,r}\Lambda_{V,r}^{-1/2}E_{V,r}^\top,
$$

$$
G_{g,\mathsf{DD}}^{-1/2}
=E_{D,r}\Lambda_{D,r}^{-1/2}E_{D,r}^\top.
$$

### §76.2 定義

$$
\boxed{
R_{g,\mathsf{VD}}
=
G_{g,\mathsf{VV}}^{-1/2}
G_{g,\mathsf{VD}}
G_{g,\mathsf{DD}}^{-1/2}.
}
$$

做 SVD：

$$
\boxed{
R_{g,\mathsf{VD}}
=U_g
\operatorname{diag}(\rho_{g,1},\dots,\rho_{g,r_g})
V_g^\top,
}
$$

其中

$$
1\ge\rho_{g,1}\ge\cdots\ge0.
$$

### 命題 76.1（CCA 等價）

$\rho_{g,k}$ 是 subset $g$ 上 activity coordinates 與 directional coordinates 的第 $k$ 個 canonical correlation。

**推導。** 考慮

$$
\max_{a,b}
\frac{a^\top G_{VD}b}
{\sqrt{a^\top G_{VV}a}\sqrt{b^\top G_{DD}b}}.
$$

令

$$
\tilde a=G_{VV}^{1/2}a,
\qquad
\tilde b=G_{DD}^{1/2}b.
$$

則 objective 化為

$$
\max_{\|\tilde a\|=\|\tilde b\|=1}
\tilde a^\top
R_{VD}
\tilde b,
$$

其最大值為 $\sigma_1(R_{VD})=\rho_1$，其後在正交限制下依序得到其餘 singular values。∎

### §76.3 Canonical modes

第 $k$ 個 activity mode 的 coefficient：

$$
a_{g,k}
=G_{VV}^{-1/2}u_{g,k}.
$$

第 $k$ 個 directional mode：

$$
b_{g,k}
=G_{DD}^{-1/2}v_{g,k}.
$$

對應 observation-space signals：

$$
v^{*}_{g,k}
=\Phi_{\mathsf V}a_{g,k},
$$

$$
d^{*}_{g,k}
=\Phi_{\mathsf D}b_{g,k}.
$$

在 $W_g$ 下：

$$
\|v^*_{g,k}\|_{W_g}
=\|d^*_{g,k}\|_{W_g}=1,
$$

且

$$
\boxed{
\langle v^*_{g,k},d^*_{g,k}\rangle_{W_g}
=\rho_{g,k}.
}
$$

因此每個 $\rho_{g,k}$ 都對應一條具體的

$$
\mathsf V_k\leftrightarrow\mathsf D_k
$$

channel。

---

## §77. Directional Explained-Variance Operator

### 定義

從 D 的角度，把所有 canonical coupling 收斂成一個 PSD operator：

$$
\boxed{
H_g^{\mathsf D\leftarrow\mathsf V}
=
R_{g,\mathsf{VD}}^\top
R_{g,\mathsf{VD}}.
}
$$

展開為

$$
\boxed{
H_g^{\mathsf D\leftarrow\mathsf V}
=
G_{g,\mathsf{DD}}^{-1/2}
G_{g,\mathsf{DV}}
G_{g,\mathsf{VV}}^{-1}
G_{g,\mathsf{VD}}
G_{g,\mathsf{DD}}^{-1/2}.
}
$$

### §77.1 Regression interpretation

令 subset 上 activity coordinate matrix 為 $V_g$，direction coordinate matrix 為 $D_g$。考慮 population linear regression：

$$
D_g=V_gB+\varepsilon.
$$

最小平方法：

$$
B^*=G_{VV}^{-1}G_{VD}.
$$

被 activity space 線性解釋的 directional covariance 為

$$
\boxed{
\Sigma_{D,\rm expl}
=G_{DV}G_{VV}^{-1}G_{VD}.
}
$$

在 D-space 自身 covariance 下 whitening：

$$
G_{DD}^{-1/2}
\Sigma_{D,\rm expl}
G_{DD}^{-1/2}
=H^{D\leftarrow V}.
$$

因此 $H^{D\leftarrow V}$ 的 eigenvalues

$$
\boxed{
\lambda_k(H^{D\leftarrow V})=\rho_k^2
}
$$

可解讀成第 $k$ 個 canonical directional mode 與 activity space 共變的 normalized fraction。

### §77.2 邊界

這裡的「explained」是 **同一 future-label system 內的線性幾何解釋**，不是 causal explanation，也不是 point-in-time predictability。

$\mathsf V$ 與 $\mathsf D$ 都是由 future label 建構，因此 $H^{D\leftarrow V}$ 只能說明 label-side representation。

---

## §78. Coupling Quality、Rank 與 Concentration

### §78.1 Total normalized coupling

定義

$$
\boxed{
Q_g^{\mathsf D\leftarrow\mathsf V}
=
\frac1{d_{D,g}^{\rm can}}
\operatorname{tr}
H_g^{\mathsf D\leftarrow\mathsf V}
=
\frac1{d_{D,g}^{\rm can}}
\sum_k\rho_{g,k}^2,
}
$$

其中 $d_{D,g}^{\rm can}$ 是通過 canonical rank cutoff 的有效 D dimension。

範圍：

$$
0\le Q_g^{\mathsf D\leftarrow\mathsf V}\le1.
$$

$Q$ 大表示 normalized D space 中有較多 energy 和 V-space 有系統性線性 coupling。

### §78.2 Effective coupling rank

若 $\sum_k\rho_k^2>0$，令

$$
p_{g,k}^{\rm can}
=
\frac{\rho_{g,k}^2}
{\sum_j\rho_{g,j}^2}.
$$

定義

$$
\boxed{
r_{g,\rm eff}^{\mathsf D\leftarrow\mathsf V}
=
\exp\left(
-\sum_kp_{g,k}^{\rm can}
\log p_{g,k}^{\rm can}
\right).
}
$$

### §78.3 Top-$m$ coupling concentration

$$
\boxed{
C_{g,m}^{\rm can}
=
\frac{\sum_{k=1}^m\rho_{g,k}^2}
{\sum_j\rho_{g,j}^2}.
}
$$

$C_{g,1}^{\rm can}$ 是最重要的 operational statistic 之一：若接近 1，表示大部分 V–D relation 可以壓縮成一條 canonical channel。

### 四種型態

| $Q^{D\leftarrow V}$ | $r_{eff}^{D\leftarrow V}$ | 解讀 |
|---|---|---|
| 高 | 低 | 強且低維 coupling；優先研究 |
| 高 | 高 | 強但多維；需要較複雜 representation |
| 低 | 低 | 只有一條弱 channel；可能無實際價值 |
| 低 | 高 | 接近 diffuse noise |

### 不建立單一 composite score

本章不定義

$$
J=q_D+\alpha Q-\beta r_{eff}
$$

之類的任意 single score。$q_D$、$Q$、rank、OOS stability 分別代表不同性質，必須分開做 eligibility 與排序診斷。

---

## §79. Ticker-Level Directional 與 Coupling Contribution

### 動機

任意 subset 搜尋有 $2^N$ 種可能，若直接 optimize $q_D$ 或 $Q^{D\leftarrow V}$ 必然 overfit。Candidate universe 必須來自預先限定的低複雜度 generator。

### §79.1 Directional contribution operator

對 observation $(i,d)$，令 reference D row 為

$$
\phi^{D}_{id}\in\mathbb R^{d_D}.
$$

定義 rank-one operator：

$$
\boxed{
K^{D}_{id}
=\phi^D_{id}\phi_{id}^{D\top}.
}
$$

對 ticker $i$ 在 discovery period 的歷史平均：

$$
\boxed{
\bar K_i^D
=
\frac1{T_i}
\sum_{d:i\in\mathcal U_d}
K_{id}^D.
}
$$

其 trace

$$
s_i^D
=\frac1{d_D}\operatorname{tr}(\bar K_i^D)
$$

是 ticker-level directional energy。

同理定義 $s_i^V,s_i^L$，以及 directional purity：

$$
\boxed{
p_i^D
=
\frac{s_i^D}{s_i^D+s_i^V+s_i^L+\epsilon}.
}
$$

### §79.2 Directional shape

為了把 magnitude 與 shape 分開：

$$
\boxed{
\widetilde K_i^D
=
\frac{\bar K_i^D}
{\operatorname{tr}(\bar K_i^D)+\epsilon}.
}
$$

$\widetilde K_i^D$ 描述 ticker 的 directional energy 分布在哪些 horizon combinations，而不看總量。

### §79.3 Cross contribution

可另外定義 observation-level V–D outer product：

$$
K_{id}^{VD}
=
\phi^V_{id}\phi_{id}^{D\top}.
$$

ticker-level average：

$$
\bar K_i^{VD}
=
\frac1{T_i}\sum_dK_{id}^{VD}.
$$

它不能直接取 singular values當作 canonical coupling，因為尚未用 subset-level $G_{VV},G_{DD}$ normalize；它只作 candidate generation 的低階 descriptor。

---

## §80. 多 Universe Candidate Generation

### 原則

Candidate generation 與 candidate evaluation 必須分開。Generator 只產生有限的、預先規格化的 family；後續 §82–§84 再做統計與 OOS 驗證。

### §80.1 D-purity / D-energy ladder

以 discovery-period ticker statistics 排序：

$$
\text{rank by }p_i^D
$$

或

$$
\text{rank by }u_i^D
=
\log(s_i^D+\epsilon)
-\frac12\log(s_i^V+\epsilon)
-\frac12\log(s_i^L+\epsilon).
$$

預先指定 size grid：

$$
\mathcal M=\{M_1,\dots,M_J\}.
$$

生成 nested universes：

$$
\mathcal S_j
=\operatorname{TopM_j}(u_i^D).
$$

這條線找的是 D-enriched candidate family。

### §80.2 Directional-shape clustering

用 $\widetilde K_i^D$ 的 Frobenius distance：

$$
d_D(i,j)
=\|\widetilde K_i^D-\widetilde K_j^D\|_F.
$$

做固定規格的 clustering（例如 average linkage 或 spectral clustering）。每個 cluster 形成 candidate universe。

這條線找的是具有不同 D shape 的非巢狀 universes。

### §80.3 VD-shape clustering

把

$$
\operatorname{vec}(\widetilde K_i^D),
\quad
\operatorname{vec}(\bar K_i^{VD})
$$

拼接後先做 discovery-only PCA，再 clustering。這可產生具有相似 raw V–D cross-structure 的 candidate groups。

### §80.4 Hybrid family

對 directional-shape cluster $c$，再用 $u_i^D$ 分 top fraction：

$$
\mathcal S_{c,q}
=
\{i\in c:u_i^D\ge Q_q(u^D\mid c)\}.
$$

因此可以同時保留：

- 不同 D directions；
- 同一 direction 的不同 enrichment levels。

### 規則 80.1（禁止 arbitrary subset search）

不得使用 genetic search、beam search 或逐 ticker greedily maximize evaluation-period $Q$，除非整個 optimizer 被完整嵌入 permutation null 且有獨立 OOS。Baseline 規格不允許。

---

## §81. Universe Eligibility：幾何層

幾何 eligibility 只定義「值得進入統計驗證的 candidate」，不是顯著性判準。

建議條件：

**U1 規模**

$$
\operatorname{median}_dN_{g,d}\ge N_{\min}.
$$

**U2 Directional enrichment**

$$
q_g^D>1.
$$

**U3 Directional purity**

$$
f_g^D>1/3.
$$

**U4 Coupling presence**

$$
\rho_{g,1}>\rho_{\min}^{can}.
$$

**U5 Coupling concentration**

$$
C_{g,1}^{can}>C_{\min}^{can}
$$

或

$$
r_{g,eff}^{D\leftarrow V}\le r_{\max}^{can}.
$$

所有 threshold 在研究開始前鎖定。通過 U1–U5 仍不表示 statistically significant。

---

## §82. Universe Null：Membership 是否真的特殊

### §82.1 虛無假設

對固定 candidate $\mathcal S_g$，虛無為：

> 在相同日期、相同 subset size、相同 liquidity / availability strata 下，隨機 membership 可以產生同等或更極端的 $q_D$、$Q^{D\leftarrow V}$ 或 concentration。

### §82.2 Size-matched stratified membership permutation

對每個 date $d$：

1. 固定 $N_{g,d}$；
2. 在母 universe 中按 liquidity / market-cap / listing-age 等**事前鎖定** strata 分層；
3. 在每層抽取與真實 candidate 相同數量；
4. 得到 $\mathcal S_g^{(b)}$；
5. 重算完整
   $$
   G_g^{(b)},
   R_{VD}^{(b)},
   \rho_k^{(b)},
   Q^{(b)},
   r_{eff}^{(b)}.
   $$

不能只 permutation $G_{VD}$，因為 canonical normalization 的 $G_{VV},G_{DD}$ 也會隨 membership 改變。

### §82.3 主要檢定量

至少對以下三個量給 p-value：

$$
T_1=q_g^D,
$$

$$
T_2=Q_g^{D\leftarrow V},
$$

$$
T_3=\rho_{g,1}^2.
$$

若研究重點是 low-rank concentration，再加入

$$
T_4=C_{g,1}^{can}.
$$

對「越大越特殊」的量：

$$
\boxed{
p_g
=\frac{1+\sum_b\mathbf1\{T_g^{(b)}\ge T_g^{obs}\}}
{B_U+1}.
}
$$

### §82.4 Multiple testing

因同時評估多個 universes，主方案使用 Westfall–Young max-T：每次 membership permutation 產生全部 candidates 的 statistic，取

$$
M^{(b)}=\max_g T_g^{(b)}.
$$

再由 max distribution 校準 family-wise error。

對 nested ladder 尤其必要，因相鄰 universes 高度相關。

### §82.5 Candidate generator 由 label 產生時

若 membership itself 是由 $u_i^D$、$\widetilde K_i^D$ 等 label-derived statistics 生成，則固定-membership null 只回答「這個已選好的集合是否特殊」，沒有包含 search cost。

正式 inference 有兩種合法方案：

1. **Discovery / evaluation split（主方案）**：generator 只在 $\mathcal D_U$ 跑，$\mathcal E_U$ 固定 membership 評估；
2. **Full rerun null（昂貴替代方案）**：每次 permutation 從 ticker scores、clustering、candidate generation 到 pruning 全部重跑。

Baseline 採方案 1。

---

## §83. Temporal Stability

### §83.1 Moving-block bootstrap

在 evaluation period 對日期做 moving-block bootstrap，block length

$$
\ell_{blk}\ge H_{\max}.
$$

每次 bootstrap 重算：

$$
q_D^{(b)},
\quad
\rho_k^{(b)},
\quad
Q^{(b)},
\quad
r_{eff}^{(b)},
\quad
V_{D,r}^{(b)}.
$$

### §83.2 Enrichment persistence

定義：

$$
\boxed{
\pi_{q_D}
=\Pr_b(q_D^{(b)}>1).
}
$$

### §83.3 Coupling persistence

對 top canonical correlation：

$$
\boxed{
\pi_{\rho_1}
=\Pr_b(\rho_1^{(b)}>\rho_{\min}^{can}).
}
$$

也報告 median 與 10/90 percentile。

### §83.4 Directional-mode stability

若保留前 $r$ 個 canonical D modes：

$$
V_{g,r}=[v_{g,1},\dots,v_{g,r}].
$$

對 bootstrap $b$：

$$
S_{g,r}^{(b)}
=
\frac1r
\|V_{g,r}^{\top}V_{g,r}^{(b)}\|_F^2.
$$

範圍 $[0,1]$。$1$ 表示 canonical directional subspace 完全穩定。

若只看第一 mode：

$$
M_{g,1}^{(b)}=|v_{g,1}^\top v_{g,1}^{(b)}|.
$$

使用絕對值是因 singular vector 的 sign 不具識別性。

---

## §84. Train / Freeze / OOS Canonical Test

### 動機

Bootstrap 仍在同一 evaluation distribution 內重抽；最重要的測試是 canonical channel 在真正未參與 discovery 的未來期間是否保持。

### §84.1 Freeze

在 train window 得到：

$$
\mathcal S_g,
\quad
G_{VV}^{train},
\quad
G_{DD}^{train},
\quad
u_{g,k}^{train},
\quad
v_{g,k}^{train},
\quad
\rho_{g,k}^{train}.
$$

對第 $k$ 個 canonical pair，freeze normalized coefficient：

$$
a_{g,k}^{train}
=(G_{VV}^{train})^{-1/2}u_{g,k}^{train},
$$

$$
b_{g,k}^{train}
=(G_{DD}^{train})^{-1/2}v_{g,k}^{train}.
$$

### §84.2 OOS evaluation

在 OOS 只 transform reference $\Phi$，計算 frozen signals：

$$
v_{g,k,t}^{OOS}
=\Phi_{\mathsf V,t}a_{g,k}^{train},
$$

$$
d_{g,k,t}^{OOS}
=\Phi_{\mathsf D,t}b_{g,k}^{train}.
$$

然後直接計算

$$
\boxed{
\rho_{g,k}^{OOS}
=
\operatorname{Corr}_{W_g^{OOS}}
(v_{g,k}^{OOS},d_{g,k}^{OOS}).
}
$$

不得在 OOS 重新做 CCA 後再報 top singular value；那會重新引入 selection。

### §84.3 Decay

$$
\boxed{
D_{\rho,g,k}
=
\frac{|\rho_{g,k}^{OOS}|}
{|\rho_{g,k}^{train}|+\epsilon}.
}
$$

以及 sign consistency：

$$
\boxed{
S_{\rho,g,k}
=\operatorname{sign}(\rho_{train})
\operatorname{sign}(\rho_{OOS}).
}
$$

由於 CCA singular value 本身無 sign，實作上必須在 train 固定 canonical pair orientation，再對 OOS bilinear covariance 使用同一 orientation。

---

## §85. 多 Universe 的 Redundancy 與保留

### 動機

多個 candidate universe 可以同時通過 eligibility，但可能只是同一組 ticker 或同一 D mode 的不同切法。Membership redundancy 與 mechanism redundancy 必須分開。

### §85.1 Membership overlap

Jaccard：

$$
\boxed{
J_{gh}
=
\frac{|\mathcal S_g\cap\mathcal S_h|}
{|\mathcal S_g\cup\mathcal S_h|}.
}
$$

Dynamic universe 報告 median-date Jaccard 與 p10/p90。

### §85.2 Canonical directional similarity

若兩個 universes 各保留前 $r_g,r_h$ 個 stable D canonical modes：

$$
V_g\in\mathbb R^{d_D\times r_g},
\qquad
V_h\in\mathbb R^{d_D\times r_h}.
$$

principal-angle similarity：

$$
\boxed{
S_{gh}^{D}
=
\frac1{r_*}
\|V_g^\top V_h\|_F^2,
\qquad
r_*=\min(r_g,r_h).
}
$$

若都只保留 top mode：

$$
S_{gh}^{D}=|v_{g,1}^\top v_{h,1}|^2.
$$

### §85.3 保留規則

Universe 不因高 membership overlap 自動刪除。只有同時：

$$
J_{gh}>\tau_J
$$

且

$$
S_{gh}^{D}>\tau_D
$$

時，才視為同一 universe family 的冗餘版本。

在同 family 中，優先保留：

1. evaluation-period adjusted significance 較強者；
2. $\rho^{OOS}$ 衰減較小者；
3. membership turnover 較低者；
4. 若以上接近，保留較大 universe。

這是一個 lexicographic rule，不建立任意加權總分。

---

## §86. Accepted Universe 的 Local Refit

### 動機

共同 $\Phi_{\mathcal U}$ 的角色是 discovery 與跨-universe 比較。當 $\mathcal S_g$ 被正式接受後，後續 feature–label–alpha pipeline 應把它視為自己的 market domain，重新建立 local basis。

### §86.1 Local label basis

在 universe $g$ 自己的 model-training period：

$$
Y_g
\rightarrow
\Sigma_{Y,g}
\rightarrow
B_{Y,g}
\rightarrow
Z_g
\rightarrow
\Phi_g.
$$

並滿足

$$
\boxed{
\Phi_g^\top W_g\Phi_g=I.
}
$$

此後 universe $g$ 的 feature geometry 使用

$$
C_g=X_g^\top W_g\Phi_g,
$$

而不是永久使用 reference $\Phi_{\mathcal U}$。

### §86.2 Reference-to-local diagnostic

將 reference basis restriction 到 subset，先在 $W_g$ 下 orthonormalize：

$$
Q_g
=\operatorname{orth}_{W_g}(\Phi_{\mathcal U}|_{\mathcal S_g}).
$$

比較 local basis：

$$
A_g=Q_g^\top W_g\Phi_g.
$$

SVD：

$$
A_g=U\Sigma V^\top,
$$

$$
\sigma_k=\cos\theta_k.
$$

$\theta_k$ 是 reference-to-local principal angle，只作診斷，不是 universe selection objective。

### §86.3 Interpretation

- 小 angle：subset 的 local mechanism space 大致是 full-universe mechanism 的重新加權；
- 大 angle：subset 有顯著不同的 local label geometry；
- 兩者都可以保留，只是下游 feature interpretation 必須使用 local $\Phi_g$。

---

## §87. 從 Canonical D Mode 到 Feature Research

### 動機

Universe geometry 到這一步仍然沒有使用任何 engineered feature。Accepted universe 的 canonical D modes 可以作為後續 feature research 的明確 target，但不能取代 full mechanism analysis。

### §87.1 Reference canonical target

對 universe $g$ 的 stable top D canonical coefficient $b_{g,1}$：

$$
\boxed{
d_g^*
=\Phi_{\mathsf D}b_{g,1}.
}
$$

它代表該 universe 中最強的 V-coupled directional combination。

### §87.2 Feature alignment

對 point-in-time feature $x_p$：

$$
\boxed{
A_{p,g}^{can}
=x_p^\top W_gd_g^*.
}
$$

必須像既有 target alignment 一樣做 permutation significance、temporal stability 與 OOS decay。

### §87.3 完整 downstream pipeline

canonical target 只負責回答一個額外問題：

> feature 是否能預測 universe discovery 階段找到的低維 D channel？

完整 feature analysis 仍然要計算：

$$
c_{p,g}=x_p^\top W_g\Phi_g,
$$

$$
Q_{p,g}^{D},Q_{p,g}^{V},Q_{p,g}^{L},
$$

fingerprint、stability、breadth、redundancy、coverage 與 alpha attribution。

---

## §88. Beta、Demean 與 Factor-Neutral Variant

### Baseline

既有 label 前處理為 within-date demean：

$$
Y^{XS}_{id,h}
=Y^{raw}_{id,h}-\bar Y_{d,h}.
$$

等價於投影掉

$$
\operatorname{span}\{\mathbf1\}.
$$

因此它是 unit-exposure common factor neutralization 的 rank-one 特例。

### 不會移除 heterogeneous beta dispersion

若

$$
r_i=\beta_i r_m+\epsilon_i,
$$

則 demean 後

$$
r_i-\bar r
=(\beta_i-\bar\beta)r_m
+(\epsilon_i-\bar\epsilon).
$$

所以 high-beta-dispersion subset 可能提高 $q_D$ 或產生 V–D coupling。

### Factor-neutral variant

若目標是研究 idiosyncratic geometry，先做

$$
Y_d^{idio}
=M_{F_d}Y_d,
$$

$$
M_{F_d}
=I-F_d(F_d^\top F_d)^{-1}F_d^\top,
$$

再完整重建

$$
\Phi_{\mathcal U}^{idio},
\quad
G_g^{idio},
\quad
R_{VD,g}^{idio}.
$$

禁止把 raw basis 與 idio basis 的 block 混合計算。

### 報表

若有 factor-neutral variant，同時報：

$$
q_D^{raw},Q_{VD}^{raw},\rho_1^{raw},
$$

$$
q_D^{idio},Q_{VD}^{idio},\rho_1^{idio}.
$$

raw 強、idio 弱代表 structure 主要由已知 common exposures 驅動；兩者都強表示 factor-neutral 後仍存在。

---

## §89. Cross-Sectional Breadth 與單一 Ticker 支配

### 動機

一個 universe 的高 $q_D$ 或高 $\rho_1$ 可能由極少數 names 決定。需要檢查 cross-sectional breadth。

### §89.1 D-energy contribution

對 date $d$：

$$
e_{id}^{D}
=\|\phi_{id}^{D}\|^2.
$$

normalized share：

$$
a_{id}^{D}
=\frac{e_{id}^{D}}
{\sum_{j\in\mathcal S_{g,d}}e_{jd}^{D}}.
$$

有效 contributor count：

$$
\boxed{
N_{g,d}^{eff,D}
=\frac1{\sum_i(a_{id}^{D})^2}.
}
$$

relative breadth：

$$
\boxed{
b_{g,d}^{D}
=\frac{N_{g,d}^{eff,D}}{N_{g,d}}.
}
$$

### §89.2 Canonical-channel breadth

對 frozen canonical pair $(a_g,b_g)$，每個 observation 的 coupling contribution：

$$
z_{id}^{can}
=(\phi_{id}^{V\top}a_g)
(\phi_{id}^{D\top}b_g).
$$

以 absolute contribution 定義 share：

$$
c_{id}^{can}
=\frac{|z_{id}^{can}|}
{\sum_j|z_{jd}^{can}|}.
$$

$$
\boxed{
N_{g,d}^{eff,can}
=\frac1{\sum_i(c_{id}^{can})^2}.
}
$$

高 $\rho_1$ 但極低 canonical breadth 應標記為 concentrated，而不是直接排除。

---

## §90. Dynamic Universe

### 動機

若 universe membership 隨時間更新，selection rule 必須只依賴當時以前已實現的 label history。

### §90.1 Rolling ticker operators

對 ticker $i$：

$$
\bar K_{i,t}^{D}
=
\sum_{\tau<t}\omega_{t-\tau}
K_{i\tau}^{D},
$$

$$
\bar K_{i,t}^{VD}
=
\sum_{\tau<t}\omega_{t-\tau}
K_{i\tau}^{VD}.
$$

$\omega$ 可以是 fixed rolling window 或 exponential decay，必須事前鎖定。

### §90.2 Membership update

$$
\mathcal S_{g,t}
=\mathcal A_g
(\{\bar K_{i,t}^{D},\bar K_{i,t}^{VD}\}_i).
$$

update frequency 優先使用 weekly / monthly，而不是每個 5-min bar 更新 membership。

### §90.3 Turnover

Universe turnover：

$$
\boxed{
\operatorname{Turn}_g(t)
=1-
\frac{|\mathcal S_{g,t}\cap\mathcal S_{g,t-1}|}
{|\mathcal S_{g,t}\cup\mathcal S_{g,t-1}|}.
}
$$

報告 median / p90 turnover。高 turnover 代表 universe identity 本身不穩定，後續 feature training 會受到 distribution drift 影響。

### §90.4 Freeze boundary

Dynamic rule 的 functional form、window、cluster count、size grid 與 thresholds 在 discovery 完成後全部 freeze。Evaluation 可以更新 statistics，但不能修改 rule。

---

# 第十部分：統計驗證、報表與實作

## §91. Universe Report 格式

每個 accepted / rejected candidate 都輸出一列；不得只保存 winners。

### A. Identity

| 欄位 | 內容 |
|---|---|
| `universe_id` | 唯一 ID |
| `generator` | ladder / D-shape / VD-shape / hybrid |
| `parent_cluster` | 若適用 |
| `size_median` | median $N_{g,d}$ |
| `size_p10`, `size_p90` | 動態 universe 規模 |
| `turnover_median` | dynamic membership turnover |

### B. Mechanism energy

| 欄位 | 定義 |
|---|---|
| `q_D` | $\operatorname{tr}(G_{DD})/d_D$ |
| `q_V` | $\operatorname{tr}(G_{VV})/d_V$ |
| `q_L` | $\operatorname{tr}(G_{LL})/d_L$ |
| `f_D,f_V,f_L` | dimension-adjusted fingerprint |
| `r_eff_D` | D spectrum effective rank |
| `C_D_1`, `C_D_2` | top-1 / top-2 D concentration |

### C. Canonical coupling

| 欄位 | 定義 |
|---|---|
| `rho_1...rho_r` | canonical correlations |
| `Q_D_from_V` | mean $\rho_k^2$ |
| `r_eff_can` | coupling effective rank |
| `C_can_1`, `C_can_2` | top coupling concentration |
| `can_D_loadings_1` | top D canonical coefficients |
| `can_V_loadings_1` | top V canonical coefficients |

### D. Inference

| 欄位 | 內容 |
|---|---|
| `p_qD_raw`, `p_qD_adj` | D enrichment null |
| `p_Q_raw`, `p_Q_adj` | coupling null |
| `p_rho1_raw`, `p_rho1_adj` | top canonical correlation null |
| `pi_qD` | temporal persistence |
| `pi_rho1` | coupling persistence |
| `subspace_stability` | principal-angle bootstrap stability |

### E. OOS

| 欄位 | 內容 |
|---|---|
| `rho1_train` | frozen train canonical correlation |
| `rho1_oos` | frozen OOS canonical covariance/correlation |
| `rho_decay` | $|\rho_{OOS}|/(|\rho_{train}|+\epsilon)$ |
| `rho_sign_consistency` | orientation-consistent sign |
| `qD_oos` | OOS D enrichment |

### F. Breadth / factor diagnostics

| 欄位 | 內容 |
|---|---|
| `breadth_D_median` | D contributor breadth |
| `breadth_can_median` | canonical contribution breadth |
| `qD_idio` | 若有 factor-neutral variant |
| `rho1_idio` | 若有 factor-neutral variant |

### G. Decision

| 欄位 | 內容 |
|---|---|
| `eligible_geometry` | U1–U5 |
| `significant` | permutation / max-T |
| `stable` | bootstrap + OOS |
| `family_id` | universe redundancy family |
| `accepted` | 最終是否進 local feature pipeline |
| `reject_reason` | 規模／null／不穩定／冗餘／breadth 等 |

---

## §92. 完整實作順序

以下可直接作為 pipeline checklist。

### 階段 0：設定

```text
 0. 鎖定 horizon band、母 universe、discovery/evaluation/OOS windows。
 1. 鎖定 candidate generators、size grid、clustering 方法與 cluster count range。
 2. 鎖定 canonical rank cutoff、B_U、B_time、所有 eligibility thresholds。
 3. 鎖定 dynamic universe rolling window / decay 與 update frequency（若使用）。
```

### 階段 1：Reference basis

```text
 4. 在母 universe、discovery period 對 label 做 within-date demean。
 5. fit Sigma_Y、K_eff、B_Y；transform -> Z。
 6. Phi_D := Z；建立 Phi_V、Phi_L；驗證 Phi^T W_U Phi = I。
 7. freeze reference fit objects。
```

### 階段 2：Candidate generation

```text
 8. 計算 ticker-level K_i^D、s_i^D/V/L、p_i^D、u_i^D。
 9. 計算 normalized D shape Ktilde_i^D；必要時 K_i^VD descriptor。
10. 依預先規格生成 ladder / D-shape / VD-shape / hybrid candidate family。
11. 移除低於 N_min 的 candidates；不依 evaluation performance 調整 family。
```

### 階段 3：Discovery geometry

```text
12. 對每個 candidate 建 W_g。
13. G_g = Phi^T W_g Phi；切 DD/VV/VD 等 block。
14. q_D/q_V/q_L、fingerprint、D eigenspectrum、r_eff_D。
15. 建 G_VV^{-1/2}, G_DD^{-1/2}（rank-controlled）。
16. R_VD；SVD -> rho_k, U_k, V_k。
17. H_D<-V = R^T R；Q_D<-V、r_eff_can、C_can_m。
18. 套用 geometry eligibility U1-U5。
```

### 階段 4：Evaluation inference

```text
19. freeze membership / membership rule 與 canonical rank choices。
20. transform evaluation labels 到 frozen reference Phi。
21. 重算 G_g、q_D、R_VD、rho、Q、rank。
22. size-matched stratified membership null；每次重算全部 canonical objects。
23. Westfall-Young max-T across candidate family。
24. moving-block bootstrap -> pi_qD, pi_rho1, canonical subspace stability。
25. D / canonical breadth diagnostics。
```

### 階段 5：真正 OOS freeze test

```text
26. 在 train/evaluation 結束時 freeze canonical coefficients a_gk, b_gk。
27. OOS 只 transform Phi；不得重跑 CCA。
28. 計算 frozen rho_OOS、rho_decay、sign consistency、qD_OOS。
29. 依 pre-registered stability rule 決定 stable / unstable。
```

### 階段 6：Multi-universe pruning

```text
30. 計算 pairwise Jaccard membership overlap。
31. 計算 canonical D-subspace principal-angle similarity。
32. 只有 membership 與 D mechanism 都高度重複時才合併 family。
33. family 內用 lexicographic rule 保留 representative；允許多 family 同時存在。
```

### 階段 7：Local refit 與 feature pipeline

```text
34. 對每個 accepted universe，在其 own model-training window 重建 local Phi_g。
35. 報告 reference-to-local principal angles。
36. 在 local Phi_g 上開始 feature library、C_g、Q、fingerprint、stability、selection。
37. 額外檢驗 feature 對 frozen canonical D target 的 alignment。
38. 模型 alpha 後續完全沿用原 Feature–Label–Alpha framework。
```

---

## §93. Fit / Transform 邊界總表

| 物件 | 擬合於 | 套用於 | 是否可在 OOS 重擬合 |
|---|---|---|---|
| $\Sigma_Y,B_Y,\Phi_{\mathcal U}$ | universe discovery train | evaluation / OOS transform | 否 |
| ticker $\bar K_i^D,\bar K_i^{VD}$ | discovery history | candidate generation | 否；dynamic 僅用 past rolling update |
| candidate generator | discovery | evaluation / OOS | 否 |
| static membership $\mathcal S_g$ | discovery | evaluation / OOS | 否 |
| dynamic membership rule | discovery | future times | 規則不可變，state 可更新 |
| canonical rank cutoff | discovery | evaluation / OOS | 否 |
| $U_g,V_g,\rho_g$ for selection | train/evaluation | diagnosis | OOS 不可重做後挑 top |
| frozen canonical coefficients $a_g,b_g$ | pre-OOS | OOS | 否 |
| local $\Phi_g$ | accepted-universe model train | downstream model OOS | 各 model fold 依原規格 |
| factor-neutral basis | 對應 train window | 對應 evaluation | 與 raw basis 分開 |

---

## §94. 數值恆等式測試

以下全部必須在 production implementation 中自動化。

| 測試 | 判準 |
|---|---|
| Reference orthogonality | $\|\Phi_U^TW_U\Phi_U-I\|_F<10^{-6}$ |
| Mother-universe energy | $|q_U^D-1|,|q_U^V-1|,|q_U^L-1|<10^{-8}$ |
| Mother cross-block | $\|G_{U,VD}\|_F<10^{-8}$ |
| PSD of $G_g$ | $\lambda_{min}(G_g)>-10^{-10}$ |
| Canonical bound | $\rho_k\le1+10^{-8}$ |
| $H$ identity | $\|H-R^TR\|_F<10^{-10}$ |
| Eigen/SVD match | eigenvalues$(H)$ = $\rho_k^2$ within $10^{-10}$ |
| CCA normalization | canonical V/D signals have unit $W_g$ norm |
| CCA covariance | $\langle v_k^*,d_k^*\rangle_{W_g}=\rho_k$ |
| Effective rank bound | $1\le r_{eff}\le$ retained rank |
| Top concentration monotonic | $C_1\le C_2\le\cdots\le1$ |
| Sign invariance | flipping singular vector signs leaves all scalar metrics unchanged |
| Membership identity | $\mathcal S_g=\mathcal U\Rightarrow G_g=I$ |

---

## §95. 合成資料校準測試

### T1：Random membership null

在真實 label geometry 上隨機生成與 candidate size 相同的 membership。

判準：在沒有植入 special subset 時，$p_{qD},p_Q,p_{\rho_1}$ 大致均勻；FWER 0.05 下 family false positive rate $\le5\%$。

### T2：植入 D-enrichment

生成一組 subset，其 D block covariance 乘以 $c>1$，但 V–D coupling 不變。

預期：

$$
q_D\uparrow,
\quad
Q^{D\leftarrow V}\approx\text{baseline}.
$$

這驗證 D enrichment 與 coupling 被正確分離。

### T3：植入 rank-one V–D channel

生成：

$$
d_1=\rho v_1+\sqrt{1-\rho^2}\epsilon,
$$

其餘 D directions 與 V independent。

預期：

$$
\hat\rho_1\approx\rho,
$$

$$
\rho_{k>1}\approx0,
$$

$$
r_{eff}^{can}\approx1.
$$

### T4：Diffuse coupling

植入 $r$ 個相等 canonical correlations。

預期：

$$
r_{eff}^{can}\approx r,
$$

且 $C_1^{can}\approx1/r$。

### T5：Subset scale only

只提高 subset 的 D 與 V variance，不改 normalized correlation。

預期：raw $\|G_{VD}\|$ 可上升，但 normalized $\rho_k$ 不應系統性上升。這驗證 canonical whitening 的必要性。

### T6：Search-cost calibration

用 discovery labels 生成多個 ladder / clusters，但 evaluation labels 與 discovery independent。

預期：若使用 discovery/evaluation split，adjusted p-value 保持校準；若在同期間 selection + evaluation，假陽性率明顯膨脹。

### T7：Frozen OOS channel

植入 train 與 test 同一 canonical channel。

判準：frozen $\rho_{OOS}$ 回復真值；若 test channel 旋轉 90°，frozen $\rho_{OOS}$ 應接近 0，即使 test 重新做 CCA 仍可找到高 top singular value。這個測試驗證「OOS 不得重新挑 canonical mode」。

### T8：Beta-driven subset

生成

$$
r_i=\beta_ir_m+\epsilon_i
$$

並讓 candidate 有較大 $\operatorname{Var}(\beta)$。

預期：raw $q_D$ 上升；factor-neutral variant 的 enrichment 消失或大幅下降。

### T9：Dynamic universe

讓 $N_d$ 與 membership 隨時間改變。

判準：$W_g$ 的每 date total weight 精確為 $1/T_g$；null calibration 不受 $N_d$ 變動破壞。

### T10：Single-name domination

植入一個極端 ticker 產生大部分 coupling。

預期：$\rho_1$ 可很高，但 canonical breadth 顯著降低。驗證 breadth diagnostic 能辨識 concentration。

---

## §96. 預設營運參數

| 參數 | 預設值 | 說明 |
|---|---:|---|
| discovery / evaluation | 60 / 40 | universe discovery inference |
| purge gap | $\ge H_{max}$ | label overlap |
| $N_{min}$ | $\max(100,10d_\Phi)$ | candidate minimum median size |
| ladder size grid | 100, 200, 300, 500, 800, 1200 | 依母 universe 調整，事前鎖定 |
| canonical eigen cutoff $\epsilon_{can}$ | $10^{-6}$ relative | 比純數值 orth cutoff 更保守 |
| $B_U$ membership null | 2000 | 若 candidates 很多至少 1000 |
| multiple testing | Westfall–Young max-T | family-wise control |
| $B_{time}$ | 1000 | temporal bootstrap |
| block length | $\ge H_{max}$ | overlapping label dependence |
| $q_D$ geometry eligibility | $>1$ | 非顯著性判準 |
| $f_D$ geometry eligibility | $>1/3$ | 非顯著性判準 |
| $\rho_{min}^{can}$ | 0.20 | 初始營運門檻，需 calibration |
| $C_{min}^{can}$ | 0.50 | top canonical pair 至少解釋半數 coupling energy |
| $\pi_{qD}$ | $\ge0.80$ | bootstrap persistence |
| $\pi_{\rho_1}$ | $\ge0.80$ | bootstrap persistence |
| D-subspace stability | $\ge0.70$ | bootstrap principal-angle similarity |
| Jaccard redundancy | $>0.80$ | 與 D similarity 聯合使用 |
| D-subspace redundancy | $>0.90$ | 與 Jaccard 聯合使用 |
| dynamic update | weekly / monthly | 降低 membership churn |
| numerical $\epsilon$ | $10^{-12}$ | denominator protection |

所有門檻都是營運篩選參數，不是統計判準。正式顯著性來自 §82 的 null；穩定度來自 §83–§84；真正 feature learnability 來自 §87 之後的 downstream pipeline。

---

# 第十一部分：統一觀點與方法邊界

## §97. 一個 Universe 的四個核心物件

對每個 candidate universe $g$，本章最終只需要記住四個物件。

### 1. Restriction geometry

$$
\boxed{
G_g=\Phi_U^TW_g\Phi_U.
}
$$

回答：reference label mechanisms 在 subset 上如何被重新加權。

### 2. Directional block

$$
\boxed{
G_{g,DD}.
}
$$

回答：D energy 有多強、集中在哪些 horizon combinations。

### 3. Canonical coupling

$$
\boxed{
R_{g,VD}
=G_{VV}^{-1/2}G_{VD}G_{DD}^{-1/2}.
}
$$

回答：V 與 D 的 normalized cross-geometry。

### 4. Directional explained-variance operator

$$
\boxed{
H_g^{D\leftarrow V}=R_{g,VD}^TR_{g,VD}.
}
$$

回答：哪些 D directions 在該 subset 中具有低維 activity-side representation。

四個物件的角色不可互換。

---

## §98. Universe、Feature、Alpha 的平行結構

既有 framework 對 signal $s$：

$$
c_s=s^TW\Phi.
$$

Universe 改變的是 measure，而不是提供一個一維 signal，因此自然物件是：

$$
G_g=\Phi^TW_g\Phi.
$$

兩者形成以下對照：

| Feature / Alpha Geometry | Universe Geometry |
|---|---|
| $c_s=s^TW\Phi$ | $G_g=\Phi^TW_g\Phi$ |
| $Q_s=\|c_s\|^2$ | $q_g^B=\operatorname{tr}(G_{BB})/d_B$ |
| mechanism fingerprint | universe fingerprint |
| $c_s/\|c_s\|$ | enriched D eigenvectors |
| feature mechanism rank / coverage | D effective rank |
| feature–label alignment | V–D canonical correlation |
| feature stability | universe canonical subspace stability |
| feature breadth | D / canonical contributor breadth |

因此完整研究順序為：

$$
\boxed{
\text{Universe Geometry}
\longrightarrow
\text{Local Feature Geometry}
\longrightarrow
\text{Alpha Geometry}.
}
$$

---

## §99. 本章能回答與不能回答的問題

### 可以回答

1. 某個 subset 是否有較高的 directional label energy？
2. Directional energy 是否集中在少數 modes？
3. Activity 與 direction 是否在 subset 中出現 normalized canonical coupling？
4. 這個 coupling 是否低秩？
5. canonical D subspace 是否跨時間穩定？
6. 不同 universes 是否濃縮不同的 D mechanisms？
7. structure 是否主要由 beta / common-factor dispersion 驅動？

### 不能回答

1. 是否存在任何 point-in-time feature 可以預測這些 D modes；
2. causal mechanism；
3. 交易後是否有 positive net Sharpe；
4. 某個 universe 是否「最佳」；
5. local model 的 optimal architecture；
6. canonical V–D relation 是否在 economic time ordering 上代表 V 先於 D。

最後一點特別重要：$\Phi_V$ 與 $\Phi_D$ 都由同一 future label system 建立。$V\leftrightarrow D$ coupling 是 representation relation，不是 temporal causal arrow。

---

## §100. 最小可行版本（MVP）

第一次實作只做以下十步：

1. 在母 universe discovery period 建 fixed $\Phi_U$。
2. 每 ticker 算 $s_i^D,s_i^V,s_i^L,p_i^D$。
3. 用預先指定的 top-$M$ grid 建 5–8 個 ladder universes。
4. 每個 universe 算
   $$G_g,\quad q_D,q_V,q_L,\quad f_D,\quad r_{eff}^D.$$
5. 算
   $$R_{VD}=G_{VV}^{-1/2}G_{VD}G_{DD}^{-1/2}.$$
6. 報
   $$\rho_1,\rho_2,\quad Q^{D\leftarrow V},\quad r_{eff}^{can},\quad C_1^{can}.$$
7. 在 evaluation 上做 size-matched membership null + Westfall–Young。
8. moving-block bootstrap 檢查 $q_D$、$\rho_1$ 與 top D canonical mode stability。
9. freeze top canonical coefficients，在真正 OOS 測 $\rho_1^{OOS}$。
10. 對通過的多個 universe 分別重建 local $\Phi_g$，再開始 feature engineering。

MVP 的 accept condition 不寫成單一分數，而是四道 gate：

$$
\boxed{
\text{D enrichment}
\land
\text{canonical coupling}
\land
\text{low-rank / stable representation}
\land
\text{OOS persistence}.
}
$$

只有通過這四道 gate 的 universe 才進入 downstream feature research；而「真正可學」仍必須由後續 feature–label–alpha OOS 結果確認。

---

# 附錄 U-A：核心公式速查

### Reference

$$
\Phi_U^TW_U\Phi_U=I.
$$

### Subset measure

$$
w_{id}^{(g)}
=\frac{\mathbf1\{i\in\mathcal S_{g,d}\}}{T_gN_{g,d}}.
$$

### Restriction Gram

$$
G_g=\Phi_U^TW_g\Phi_U.
$$

### D enrichment

$$
q_g^D=\frac1{d_D}\operatorname{tr}(G_{g,DD}).
$$

### D effective rank

$$
r_{eff}^D
=\exp\left(-\sum_kp_k^D\log p_k^D\right),
\qquad
p_k^D=\frac{\lambda_k^D}{\sum_j\lambda_j^D}.
$$

### Canonical V–D operator

$$
R_{g,VD}
=G_{g,VV}^{-1/2}
G_{g,VD}
G_{g,DD}^{-1/2}.
$$

### Canonical correlations

$$
R_{g,VD}=U\operatorname{diag}(\rho_k)V^T.
$$

### D explained-variance operator

$$
H_g^{D\leftarrow V}
=R_{g,VD}^TR_{g,VD}.
$$

### Coupling quality

$$
Q_g^{D\leftarrow V}
=\frac1{d_{D,g}^{can}}\sum_k\rho_k^2.
$$

### Coupling effective rank

$$
r_{eff}^{D\leftarrow V}
=\exp\left(-\sum_kp_k^{can}\log p_k^{can}\right),
$$

$$
p_k^{can}=\frac{\rho_k^2}{\sum_j\rho_j^2}.
$$

### Top-$m$ concentration

$$
C_{g,m}^{can}
=\frac{\sum_{k=1}^m\rho_k^2}{\sum_j\rho_j^2}.
$$

### Universe mechanism similarity

$$
S_{gh}^{D}
=\frac1{r_*}\|V_g^TV_h\|_F^2.
$$

### OOS frozen channel

$$
\rho_{g,k}^{OOS}
=\operatorname{Corr}_{W_g^{OOS}}
(\Phi_Va_{g,k}^{train},\Phi_Db_{g,k}^{train}).
$$

---

# 附錄 U-B：判讀速查表

| 現象 | 幾何含意 | 下一步 |
|---|---|---|
| $q_D\gg1$, $Q_{D\leftarrow V}\approx0$ | D energy 大但與 activity 無低維 coupling | 可保留作 D-enriched domain；不要稱 structural coupling |
| $q_D>1$, $Q$ 高, $r_{eff}^{can}\approx1$ | D 強且幾乎由單一 V–D channel 表徵 | 優先做 frozen OOS canonical test |
| $q_D\approx1$, $Q$ 高 | D 不特別強，但 V–D relation 特別 coherent | 可作 specialized universe，而非 D-enrichment universe |
| $q_D$ 高, $r_{eff}^D$ 高 | D 強但 diffuse | 需要多維 features / model；不是低維 representation |
| $\rho_1^{train}$ 高, $\rho_1^{OOS}\approx0$ | canonical mode selection overfit / regime-specific | reject stability |
| raw 強、idio 弱 | common-factor dispersion 驅動 | 標記 factor-driven |
| membership overlap 高、D similarity 低 | 同 names 但不同 mechanism | 可同時保留 |
| membership overlap 低、D similarity 高 | 不同 names 濃縮同一 D mode | 視研究目的可保留兩個 domain 或合併 family |
| $\rho_1$ 高但 canonical breadth 極低 | 少數 ticker 支配 | concentration warning |
| local-reference angle 大 | subset 有自己的 label basis | downstream 必須使用 local $\Phi_g$ |

---

# 附錄 U-C：一句話總結

本章的完整邏輯可以壓縮成：

$$
\boxed{
\mathcal U
\xrightarrow{\text{reference }\Phi}
\mathcal S_g
\xrightarrow{G_g}
G_{DD},G_{VV},G_{VD}
\xrightarrow{\text{CCA}}
H_g^{D\leftarrow V}
\xrightarrow{\text{null + stability + OOS}}
\mathscr S_{accepted}
\xrightarrow{\text{local refit}}
\Phi_g
\xrightarrow{X,\alpha}
\text{true predictive learnability}.
}
$$

其中 universe discovery 的核心問題不是「哪一群股票未來 return 最大」，而是：

> **哪些 market sub-domains 的 directional label structure 被濃縮成穩定、低維、可由 activity geometry 表徵的 canonical modes，並且這個結構在真正 OOS 中仍然存在？**

