# Cross-Sectional Predictive Geometry

### 一套用於 Feature 篩選、Alpha 歸因與模型診斷的統一幾何方法

---

## 本書要解決的問題

在 cross-sectional 量化研究中，一個訊號與未來報酬的關係幾乎總是被壓縮成單一個數字：rank IC。這個數字有一個根本的問題——**它把三種在交易上完全不同的預測力混在一起**：

1. 「誰會漲得比誰多」（方向）；
2. 「誰會動得比誰劇烈」（幅度、波動、活躍度）；
3. 「在大行情發生時誰的方向是對的」（條件性方向）。

一個 IC 為 0.02 的 feature，可能是一個乾淨的方向性訊號，也可能是一個波動度預測器透過報酬的偏態洩漏出來的假象。兩者需要完全不同的部位設計與風險預算，但在 IC 的尺度上無法分辨。

本方法的做法是：**不把未來報酬當成一個純量，而把它提升成一個有限維的幾何空間**，並在這個空間中為每一個訊號指派一個座標向量。有了座標向量，「這個 feature 預測什麼」「這兩個 feature 是不是同一件事」「模型是不是只是在重組已知的 feature」這些問題，都變成向量之間的幾何運算：長度、角度、投影、張成。

整套方法最終歸結為一個物件。對任何一個 cross-sectional 訊號 $s$（可以是 raw feature、模型的 alpha、或 alpha 扣除線性成分後的殘差），定義

$$c_s=s^\top W\Phi,$$

其中 $\Phi$ 是由未來報酬建構出來的標準正交 label mechanism basis，$W$ 是每個 cross-section 等權的列權重。**本書其餘所有內容，都是在說明如何正確地建構 $\Phi$、如何判斷 $c_s$ 的哪些成分不是雜訊、以及如何解讀它。**

---

## 讀者與使用方式

本書假設讀者熟悉線性代數（特別是奇異值分解與投影算子）、基礎的置換檢定與 bootstrap，以及 cross-sectional 量化研究的一般工作流程。不假設讀者熟悉任何特定的 factor model。

章節依邏輯依賴排列，建議依序閱讀。若時間有限：

| 目的 | 最少需要讀的章節 |
|---|---|
| 理解整體構想 | §1、§8、§11、§16、§67 |
| 實作 label basis | §4–§10 |
| 實作顯著性檢定 | §13–§18 |
| 實作 feature selection | §26–§27、§34–§39 |
| 只做模型診斷（feature 側已完成） | §41–§46、§54–§60 |

每一節的結構固定為：**動機**（為什麼需要這個物件）→ **定義**（精確規格）→ **推導或性質**（為什麼它是良好定義的、它的界在哪裡）→ **實例**（具體數字）→ **解讀**（報表上看到這個數字該怎麼想）。

書中所有數值例子都已校準到股票市場的真實量級。這一點需要特別強調：本方法的核心統計量 $Q$ 的量級約為「平均每期 IC 的平方」，因此一個很好的 feature 對應的 $Q$ 是 $10^{-3}$ 量級而不是 $10^{-1}$ 量級。把量級記錯一到兩個數量級，會導致所有門檻設錯，進而誤刪全部真實訊號或誤留全部雜訊。

---

## §0. 符號約定

本書採用一組經過消歧義的符號。與較早期規格的對應關係列於本節末，以便既有實作碼對照。

### §0.1 索引與維度

| 符號 | 意義 | 範圍 |
|---|---|---|
| $d$ | cross-sectional 觀測時點（「日」為慣用語，可為任何頻率） | $1,\dots,T$ |
| $i,j$ | ticker / instrument | $i\in\mathcal U_d$ |
| $p,q$ | feature 索引 | $1,\dots,P$ |
| $h$ | horizon 索引 | $1,\dots,H$ |
| $\beta$ | horizon band 索引 | $1,\dots,B_{\rm band}$ |
| $b$ | permutation / bootstrap 重抽索引 | $1,\dots,B_\bullet$ |
| $k$ | feature transform 索引（§4.3） | $1,\dots,|\mathcal T_p|$ |
| $\mathcal U_d$ | 第 $d$ 期的 universe，$N_d=|\mathcal U_d|$ | $N_d$ 可隨 $d$ 變動 |
| $T$ | 總期數（**僅**此意義） | |
| $K$ | 當前 band 保留的 label 維度（§7 之後即為 $K_{\rm eff}$） | |
| $L$ | long-format 總列數 $\sum_dN_d$ | |
| $d_{\mathsf D},d_{\mathsf V},d_{\mathsf L}$ | 三個 block 的維度 | $K,\ K(K{+}1)/2,\ K$ |
| $d_\Phi$ | $d_{\mathsf D}+d_{\mathsf V}+d_{\mathsf L}$（扣除秩截斷後的實際值） | |

三個 predictive mechanism block 一律以**無襯線體**標記，以與其他同字母符號區隔：

$$\mathsf D\ (\text{Direction}),\qquad \mathsf V\ (\text{Volatility / Activity}),\qquad \mathsf L\ (\text{Large-move directional}).$$

### §0.2 資料、基底與變換

| 符號 | 意義 |
|---|---|
| $X_d\in\mathbb R^{N_d\times P}$ | 原始 feature 面板 |
| $X^{XS}$ | within-date 標準化並全域正規化後的 feature（§4） |
| $Y_d^{raw}\in\mathbb R^{N_d\times H}$ | 原始未來報酬 |
| $Y^{XS}$ | within-date 去均值後的 label（§5） |
| $W=\operatorname{diag}(w_{id})$ | 列權重，$w_{id}=1/(TN_d)$（§3），實作上不建構 |
| $\Sigma_Y$ | label 的加權共變異數（§7） |
| $E,\Lambda$ | $\Sigma_Y=E\Lambda E^\top$ 的特徵分解 |
| $B_Y$ | whitening 矩陣（§7），$Z_d=Y_d^{(\beta)}B_Y$ |
| $Z$ | whitened label，滿足 $Z^\top WZ=I_K$（**精確**成立） |
| $G_{\mathsf D},G_{\mathsf V},G_{\mathsf L}$ | 三個 lifting map 的輸出（§8），皆已 within-date 去均值 |
| $\Theta_{\mathsf V},\Theta_{\mathsf L}$ | 正交化矩陣（§9） |
| $\Gamma_{\mathsf{DV}},\Gamma_{\mathsf{DL}},\Gamma_{\mathsf{VL}}$ | Gram–Schmidt 投影係數（§9） |
| $\Phi=[\Phi_{\mathsf D}\,|\,\Phi_{\mathsf V}\,|\,\Phi_{\mathsf L}]$ | label mechanism basis，$\Phi^\top W\Phi=I_{d_\Phi}$ |

### §0.3 核心統計量

| 符號 | 意義 |
|---|---|
| $C=X^\top W\Phi$ | predictive operator（§11） |
| $C=USV^\top$ | SVD；$V_r$ 為顯著右奇異子空間（§12、§14） |
| $P_V=V_rV_r^\top$ | 顯著機制子空間的投影算子 |
| $r_{\rm sig},r_{\rm eff}$ | 顯著秩、有效譜秩 |
| $c_s^{\rm full}=s^\top W\Phi$ | 訊號 $s$ 的完整表徵（§43） |
| $c_s=c_s^{\rm full}P_V$ | 投影到顯著子空間的表徵 |
| $Q_s=\|c_s\|^2$ | quality（能量），$Q^{\mathsf D},Q^{\mathsf V},Q^{\mathsf L}$ 為各 block 分量 |
| $f_s,\tilde f_s$ | raw fingerprint、null-adjusted enrichment fingerprint（§19） |
| $a,\ y^*$ | intended target 的 horizon 權重、目標訊號 $y_d^*=Y_d^{(\beta)}a$ |
| $r_{\mathsf D}$ | target 在 $\Phi_{\mathsf D}$ 座標中的單位方向（§20） |
| $A_s^{\rm tgt},\ \eta_s^{\mathsf D}$ | signed target alignment、directional target efficiency |
| $A^{\rm mech}_{pq},S^{\rm mech}_{pq}$ | mechanism 有號／無號相似度（§23） |
| $S^{raw}_{pq}$ | 原始 Spearman 冗餘度（§25） |
| $\pi^{\rm time},\pi^{\rm XS},M^{\rm time},M^{\rm XS}$ | 穩定度（§27、§29） |
| $N^{\rm eff}_{d,p},\ b_{d,p},\ \tilde b_p$ | 有效貢獻者數、正規化 breadth、相對 breadth（§30–§31） |
| $R^{\rm canc}_{d,p},\ \tilde R_p$ | cancellation ratio、相對 cancellation（§32） |
| $\tilde p_\bullet$ | 多重檢定調整後的 p-value（§18） |
| $\gamma_{\rm mode},\gamma_{\rm FWER},\gamma_{\rm FDR}$ | 各層檢定的顯著水準 |
| $\hat\alpha,\tilde\alpha$ | 模型原始／standardized 預測 |
| $\varepsilon_\alpha$ | OOS linear surrogate 殘差（§58） |
| $\nu_s$ | novelty fraction（§43） |
| $\epsilon$ | 純數值下限（分母保護），預設 $10^{-12}$ |

### §0.4 與舊符號的對應

| 本書 | 舊規格 | 變更理由 |
|---|---|---|
| $\mathsf L$ block | $T$ block | $T$ 已用於總期數與轉置 |
| $E$（$\Sigma_Y$ 特徵向量） | $V$ | $V$ 已用於 SVD 右奇異向量 |
| $\mathsf V$ block | $V$ block | 同上 |
| $\Theta_\bullet$ | $A_D,A_V,A_T$ | $A$ 已用於 alignment |
| $\Gamma_\bullet$ | $B_{DV},B_{DT},B_{VT}$ | $B$ 已用於重抽次數 |
| $N^{\rm eff}_{d,p}$ | $B_{d,p}$ | 同上，且語意更明確 |
| $\tilde p_p$ | $q_p$ | $q$ 保留為 feature 索引 |
| $\gamma_\bullet$ | $\alpha_{\rm mode}$ | $\alpha$ 保留給模型預測 |
| $S^{\rm mech}$ | $S^{pred}$ | 與「prediction」區隔 |

---

# 第一部分：問題設定與資料契約

## §1. 研究目標

### 動機

任何 cross-sectional alpha 研究，實際上同時在回答三個彼此獨立的問題，而這三個問題在實務上常被混為一談：

1. **哪些 feature 真的有預測力，而且預測的是什麼？**
2. **模型的輸出（alpha）到底在預測什麼？它和 feature 的預測內容是什麼關係？**
3. **模型有沒有產生任何 feature library 用線性方式看不到的東西？**

第三個問題是最難也最重要的。一個非線性模型當然會產生非線性的輸出，但「輸出是非線性的」與「預測內容是新的」是兩件完全不同的事。一個把三個動能 feature 做了複雜交互作用、最後仍然只預測同一個動能方向的模型，它的輸出高度非線性，但它沒有發現任何新東西。本方法提供把這兩者分開的精確工具。

### 三個物件

$$X_{id,p}\ \text{（feature）},\qquad Y_{id,h}\ \text{（future return）},\qquad \hat\alpha_{id}\ \text{（model prediction）}.$$

### 三個分析問題

**問題一：Feature–Label。**

$$X\longrightarrow\Phi(Y).$$

每個 feature 是否攜帶未來報酬的資訊？攜帶的是哪一類（方向／活躍度／大行情方向）？它是否穩定？它是否與其他 feature 重複？

**問題二：Label–Alpha。**

$$\hat\alpha\longrightarrow\Phi(Y).$$

模型到底預測了什麼？它的預測力是方向性的嗎？它的方向是否對準我們真正想交易的 horizon 組合？

**問題三：Feature–Alpha。**

$$X\longleftrightarrow\hat\alpha.$$

模型是重組已知 feature 的機制，還是建立了新的機制？這個問題有三個層次的答案（§51–§59），分別對應「輸出長得像不像」「機制像不像」「有沒有超出線性可重建的部分」。

### 統一座標系

三個問題必須在**同一個座標系**中回答，否則答案無法相互對照。本方法的核心設計就是：由 $Y$ 建構一組固定的標準正交基底 $\Phi$，然後把 feature、alpha、殘差全部投影到 $\Phi$ 上：

$$X_p\longmapsto c_p,\qquad \tilde\alpha\longmapsto c_\alpha,\qquad \varepsilon_\alpha\longmapsto c_\varepsilon.$$

三者是同一個空間中的三個向量，因此可以談論它們的長度比較、夾角、以及一個是否落在另一個的張成之中。

### 方法邊界

本方法有三個必須在一開始就講清楚的邊界。它們不是缺陷，而是這個工具的適用範圍。完整討論見附錄 A。

**邊界一：label 側的函數族是有限且被指定的。** $\Phi$ 由三個 lifting map 張成（§8）。未來報酬的某些函數（例如極端尾部的符號、跨日路徑形狀）不在這個空間內，任何只預測它們的訊號在本方法下會得到 $Q\approx0$。

**邊界二：訊號側的映射是線性的。** $c_s=s^\top W\Phi$ 對 $s$ 是線性的。若 $Y=X^2+\text{noise}$ 且 $X$ 對稱分布，則 $\mathbb E[XY]=\mathbb E[X^3]=0$，因此 raw $X$ 的 $Q\approx0$，儘管 $X$ 擁有完整的資訊。本方法對此的正式處理是 §4.3 的 **transform family**：每個語意 feature 展開為一組變換，逐一分析後再在 source 層級彙總。這把邊界從「$X$ 本身是否線性相關」推進到「$X$ 的指定變換族中是否有任何一個線性相關」，但仍然不是「任何 $X$ 的函數」。

**邊界三：所有結論是 association，不是 causation。** $Q_p$ 大只代表 feature 與未來報酬機制有穩定的橫截面關聯。它可能來自真實的定價機制，也可能來自共同暴險於某個未被建模的風險因子。本方法不區分這兩者。

---

## §2. 資料結構

### 定義

令 $d=1,\dots,T$ 表示 cross-sectional observation time。可以是日、5 分鐘、1 小時，只要每個 $d$ 都代表一次完整的橫截面。

第 $d$ 期的 universe 為 $\mathcal U_d$，$N_d=|\mathcal U_d|$，允許隨時間改變。每一個 observation 必須由 $(d,i)$ 唯一識別。資料為

$$X_d\in\mathbb R^{N_d\times P},\qquad Y_d^{raw}\in\mathbb R^{N_d\times H},\qquad \hat\alpha_d\in\mathbb R^{N_d}.$$

### 資料契約

> **契約 2.1（keyed alignment）**：任何矩陣運算都必須依照 `(time, ticker)` key 對齊，不得假定不同日期的 row position 對應同一 ticker。

這條契約看似瑣碎，卻是整套方法中最容易在實作中被破壞、且破壞後最難察覺的一條。典型失敗模式：universe 檔案依市值排序，而市值排序每天變動，於是 `X[d][0]` 今天是 AAPL、明天是 MSFT。所有下游統計量仍會順利算出來，數值也看似合理，只是全部沒有意義。§65 的 Dynamic universe 測試專為此而設。

> **契約 2.2（point-in-time）**：$X_{id}$ 只能使用 $d$ 時點（含）之前公開可得的資訊；$Y_{id,h}$ 只能使用 $d$ 之後的報酬。任何 restatement、任何以最終版資料回填歷史的欄位，都構成洩漏。

### 規模感

以美股 intraday 研究為例：$d$ 為 5 分鐘 bar，一天 78 個 bar，兩年約 $T\approx39{,}000$；$N_d$ 在 1200–3000 之間浮動；$P=400$；$H=10$。則 $L\approx7\times10^7$。

但本方法**從不需要 $O(L^2)$ 的運算或記憶體**。所有「全樣本」量都表達為逐期量的加權平均，這是 §3 權重設計的直接結果。實務上，整個管線的計算瓶頸只有一次面板掃描（計算 $\{C^{(d)}\}$）；此後所有分析——SVD、permutation、bootstrap、clustering、selection——都只在一個 $P\times d_\Phi$（典型為 $400\times18$）的小矩陣上進行。這是本方法在計算上極為輕量的原因，也是為什麼 §13 的 permutation 必須以 date-shift 而非 row-level permutation matrix 實作。

---

## §3. 日期權重：每個 cross-section 等權

### 動機

在計算「全樣本的 feature–label 關聯」之前，必須先回答一個問題：**基本統計觀測單位是什麼？**

有兩個自然的答案。其一是把每一個 $(d,i)$ 當作一個觀測，於是股票多的日期自然佔較大權重。其二是把「一天的 cross-section」當作一個觀測，於是每一天等權。

本方法採用後者。理由是：cross-sectional 研究中真正的獨立重複來源是**時間**，不是個股。同一天 3000 支股票的報酬共享同一個市場衝擊、同一組產業輪動，其有效獨立樣本數遠小於 3000。若採用逐列等權，則 universe 較大的時期（通常是牛市末期、IPO 潮）會主導全樣本統計量，而這與訊號品質無關，純粹是資料供給的副產品。

### 定義

$$\omega_d=\frac1T,\qquad\boxed{w_{id}=\frac1{TN_d}.}$$

因此

$$\sum_{i\in\mathcal U_d}w_{id}=\frac1T,\qquad\sum_{d,i}w_{id}=1.$$

定義 weighted inner product

$$\boxed{\langle A,B\rangle_W=A^\top WB,}\qquad W=\operatorname{diag}(w_{id}).$$

實作時不建立 $W$。所有量都以

$$\frac1T\sum_d\frac1{N_d}(\cdot)$$

計算。

### 實例

設 $T=2$，第一天 3000 檔、第二天 1200 檔。逐列等權下第一天佔 $3000/4200=71.4\%$；本方法下兩天各佔 $50\%$。若第一天恰為一次崩盤且該 feature 當天表現異常，逐列等權會把這個單一事件的影響放大 1.43 倍。

### 等權的精確範圍，以及 $T_{\rm eff}$ 診斷

上述權重保證的是**列數等權**：每一期在加總中佔的權重相同。但它不保證**影響力等權**，原因在於 label 側不做每日波動度正規化（§5 的設計選擇）。

具體地說：$\Phi_{\mathsf D}\propto Z$，而 $Z$ 的 scale 由全樣本共變異數決定，因此在高波動日 $\|z_{id}\|$ 系統性偏大。由於

$$C^{(d)}=\frac1{N_d}X_d^{\top}\Phi_d,$$

一個 cross-sectional dispersion 為平時兩倍的日子，對 $C$ 的貢獻量級也是兩倍。對 $\mathsf V$ 與 $\mathsf L$ block 情況更強：$G_{\mathsf V}\sim\|z\|^2$、$G_{\mathsf L}\sim\|z\|^2$，故貢獻量級為**四倍**。

換句話說，$Q^{\mathsf D}$ 是**以當日 dispersion 加權的平均 IC 平方**，$Q^{\mathsf V},Q^{\mathsf L}$ 則是以 dispersion 平方加權。這是刻意的設計：未來報酬的橫截面離散度本身就是我們要分析的資訊之一，把它正規化掉會讓「大行情中的方向性」這個概念失去意義。但既然如此，就必須有一個診斷量告訴我們「實際上有多少天在決定結果」。

> **定義 3.1（有效日數）**：令
> $$\theta_d^{\mathsf B}=\frac{\frac1{N_d}\|\Phi_{\mathsf B,d}\|_F^2}{\sum_{d'}\frac1{N_{d'}}\|\Phi_{\mathsf B,d'}\|_F^2},\qquad\boxed{T^{\mathsf B}_{\rm eff}=\frac1{\sum_d(\theta_d^{\mathsf B})^2}}$$
> 其中 $\mathsf B\in\{\mathsf D,\mathsf V,\mathsf L\}$。

這是 $\{\theta_d\}$ 分布的 Herfindahl 倒數，即「等效的貢獻日數」。若 $T=1000$ 而 $T^{\mathsf D}_{\rm eff}=180$，代表整個 $Q^{\mathsf D}$ 實際上只由約 180 個「有效日」決定。此時所有 permutation 與 bootstrap 的有效樣本數都必須據此重新解讀：名義上 1000 天的 moving block bootstrap，實質上只有約 180 天的獨立資訊。

$T_{\rm eff}$ 必須列入 §40 的報表。典型的健康值是 $T_{\rm eff}/T\in[0.3,0.7]$；$T^{\mathsf V}_{\rm eff}/T<0.1$ 代表 $\mathsf V$ block 的全部結論實際上由少數幾次波動度事件決定，應降級為探索性。

### 實作註記

若研究目的**不是**分析大行情，而是要求嚴格的每日等影響力，可以在 §7 之後對 $Z$ 額外做一次 per-date scaling $Z_d\leftarrow Z_d/\hat\sigma_d$，其中 $\hat\sigma_d^2=\frac1{KN_d}\|Z_d\|_F^2$。代價是 $\mathsf L$ block 的語意從「大行情中的方向性」退化為「相對大行情中的方向性」，且 $\mathsf V$ block 只剩下 horizon 之間的相對結構。兩者都是合理的設計，但必須二擇一並記錄於 §66 的設定表；本書其餘部分假設不做這個 scaling。

---
# 第二部分：Label mechanism space 的建構

這一部分建立 $\Phi$。它是整套方法的座標系，之後所有章節都在這個座標系裡工作。建構順序是：feature 前處理（§4）→ label 前處理（§5）→ 分 horizon band（§6）→ whitening 與秩控制（§7）→ 三個 lifting（§8）→ 階層正交化（§9）→ fit/transform 邊界（§10）。

## §4. Feature 前處理

### 動機

feature 側的前處理有兩個目的，必須分開理解：

1. **穩健性**：移除明顯的 bad print 與數值錯誤，使後續的矩量估計不被單一資料錯誤主導。
2. **可比性**：使 $Q_p=\|c_p\|^2$ 可以跨 feature 直接比較。這要求每個 feature 在加權內積下有相同的尺度。若缺少這一步，一個以基點計價的 feature 與一個以百分比計價的 feature，其 $Q$ 會差 $10^4$ 倍，而這與預測力毫無關係。

### §4.1 定義

對每個日期 $d$、feature $p$，取有效 observation $\{X_{idp}:i\in\mathcal U_d\}$。

**第一步：winsorize。** 預設使用 cross-sectional 1% / 99%：

$$X'_{idp}=\operatorname{clip}\big(X_{idp},\ q_{0.01,d,p},\ q_{0.99,d,p}\big).$$

**第二步：within-date 標準化。**

$$\mu_{dp}=\frac1{N_d}\sum_iX'_{idp},\qquad \sigma_{dp}=\sqrt{\frac1{N_d}\sum_i(X'_{idp}-\mu_{dp})^2},$$

$$\boxed{X^{XS}_{idp}=\frac{X'_{idp}-\mu_{dp}}{\sigma_{dp}+\epsilon}.}$$

**第三步：全域正規化。** 缺值填補與覆蓋率差異會使各 feature 的實際加權範數不同（見 §4.2），因此在進入 §11 之前做一次全域縮放：

$$\|X_p\|_W^2=\frac1T\sum_d\frac1{N_d}\sum_i\big(X^{XS}_{idp}\big)^2,\qquad\boxed{X_p\leftarrow\frac{X_p}{\|X_p\|_W}.}$$

正規化前的 $\|X_p\|_W^2$ 本身列入報表。經過第三步後

$$\boxed{X_p^\top WX_p=1\quad\text{（精確）}}$$

對每個 $p$ 成立，這是 §16 中 $Q_p\le1$ 成為嚴格不等式的前提。

### §4.2 缺值處理

記錄每日覆蓋率

$$m_{dp}=\frac{\#\{\text{valid names}\}}{N_d}.$$

規則：

- $m_{dp}<0.8$：該日該 feature 不參與貢獻，該日記入 $\mathcal D_p^c$（缺席日集合）。
- 其餘 missing value 在標準化空間填 $0$（cross-sectional neutral）。
- 報表記錄 $\text{coverage}_p=\operatorname{median}_dm_{dp}$ 與 $\text{active\_ratio}_p=T_p/T$，其中 $T_p=|\mathcal D_p|$ 為有效日數。

**在有效日上平均，而非在全期上平均。** §11 的定義為

$$C_{p,:}=\frac1{T_p}\sum_{d\in\mathcal D_p}C^{(d)}_{p,:},$$

即缺席日不以 0 計入分子。理由是：$Q_p$ 應衡量「當這個 feature 可用時它有多強」，而「它多常可用」是一個獨立的營運問題，兩者不該相乘後混為一談。代價是 $C$ 的不同列在不同日期集合上平均，因此 §12 的 SVD 嚴格來說不是單一線性算子的分解。當 feature library 的覆蓋率分布很分散時（例如一半的 feature 覆蓋率 100%、另一半 50%），應把低覆蓋 feature 分池處理，各自跑完整管線。

**覆蓋率的隱性篩選效應值得警惕。** 一個「過去 20 日分析師預估修正幅度」的 feature，只有約 60% 的 universe 有分析師覆蓋，在 $m_{dp}<0.8$ 的規則下幾乎每天都被排除，$Q_p\to0$——不是因為沒有預測力，而是因為門檻。這類 feature 應在進入管線前先歸入「限定 universe」的子研究，而非與全 universe feature 同池競爭。

### §4.3 Transform family

**動機。** $c_p=X_p^\top W\Phi$ 對 $X_p$ 是線性的。因此若 feature 與未來報酬的關係是 U 型的（例如極端的估值偏離在兩個方向上都預示反轉），則 $\mathbb E[X\cdot z]\approx0$，本方法會回報 $Q_p\approx0$，儘管這個 feature 擁有完整的預測資訊。

解決方式不是放棄線性（那會摧毀整個幾何架構），而是**把非線性搬到 feature 側**：每個語意 feature 展開為一組變換，每個變換各自作為 library 的一個條目走完整個管線。

**定義。** 對每個 source feature $p$，定義 transform family

$$\boxed{\mathcal T_p=\{t_1,t_2,\dots\},\qquad X_{p,k}=\text{標準化}\big(t_k(X^{XS}_p)\big).}$$

初始建議的最小族為三個：

$$t_1(x)=x,\qquad t_2(x)=|x|,\qquad t_3(x)=x|x|.$$

三者的分工是：$t_1$ 捕捉單調線性關係；$t_2$ 捕捉「偏離程度」（對稱、無方向），它正是 U 型關係的線性化；$t_3$ 是奇函數但凸性放大，捕捉「極端值的方向效應比中間值強」這種關係。注意 $t_2$ 的輸出必須**重新 within-date 去均值與標準化**，否則它會帶有一個很大的日間共同成分（當日橫截面離散度），該成分依引理 4.1 對 $C$ 無貢獻，但會浪費 $X$ 的變異。

需要更完整的函數族時，可改用 rank 空間的 B-spline 基底（3–5 個節點），代價是每個 source feature 產生 3–5 個條目，library 規模隨之膨脹，多重檢定的負擔加重。

**索引與彙總。** library 的條目索引為 pair $(p,k)$；$Q,A^{\rm tgt},\pi,M$ 等全部統計量都在條目層級計算。source 層級的彙總只用於報表與最終的敘事：

$$Q_p^{\rm src}=\max_kQ_{p,k},\qquad k_p^*=\arg\max_kQ_{p,k}.$$

**選擇規則。** 同一 source 的不同變換之間天然高度相關，因此：

> **規則 4.1**：同一 source feature 的所有變換條目，在 §36 的 raw redundancy 分群中**強制歸入同一 family**，除非它們落在不同的 mechanism group（§35）且各自通過 eligibility（§34）。換言之，一個 source feature 原則上只貢獻一個被選中的條目；只有當 $t_1$ 與 $t_2$ 明確預測**不同的機制**（例如 $t_1$ 預測方向、$t_2$ 預測活躍度）時，兩者才可以同時入選。

這條規則把 transform family 的成本控制住：library 的名目大小乘以 $|\mathcal T|$，但最終選出的 feature 數不會膨脹。

### §4.4 引理 4.1（每日去均值的吸收性質）

由於 $X^{XS}_d$ 每期均值為零，對任何矩陣 $\Phi_d$ 皆有

$$\frac1{N_d}X_d^{XS\top}\Phi_d=\frac1{N_d}X_d^{XS\top}\big(\Phi_d-\mathbf 1_{N_d}\bar\Phi_d\big),$$

其中 $\bar\Phi_d$ 為 $\Phi_d$ 的當期列平均。

**證明**：$X_d^{XS\top}\mathbf 1_{N_d}=0$，故 $X_d^{XS\top}\mathbf 1\bar\Phi_d=0$。∎

**推論 4.2**：$C^{(d)}$ 自動是一個**當期 cross-sectional 共變異數**，而非未中心化的交叉動差；且 $\Phi$ 中任何「當期為常數」的成分對 $C$ 的貢獻恆為零。

這個引理在 §8、§9、§13、§29 各處都會用到。它同時解釋了兩件事：為什麼 §8.2 中的 $-I_K$ 項不影響任何結果，以及為什麼 permutation 時取 ticker 交集後**必須重新中心化**（因為在交集上 $X^\top\mathbf 1\ne0$，引理的前提被破壞）。

---

## §5. Label 前處理

### 定義

對每一個 horizon：

$$\boxed{Y^{XS}_{id,h}=Y^{raw}_{id,h}-\frac1{N_d}\sum_{j\in\mathcal U_d}Y^{raw}_{jd,h}.}$$

因此 $\sum_iY^{XS}_{id,h}=0$。

### 兩個設計選擇及其後果

**選擇一：去均值，不做 beta 中性化。** 去均值移除的是當期的 equal-weighted 市場報酬。因此整個分析的語意固定為「**相對於等權市場的橫截面對比**」。這意味著：

- 一個純粹的 market timing 訊號（預測整體市場漲跌）在本方法下 $Q=0$。這是正確的——它不是 cross-sectional 訊號。
- 若 universe 的產業組成有系統性偏斜，去均值不會移除產業層級的共同變動。要分析純 idiosyncratic 的預測力，應在此步之後額外對已知因子暴險做橫截面迴歸並取殘差。這是一個可選的前處理層，不改變後續任何公式。

**選擇二：不做每日波動度正規化。** 亦即不令 $Y_d\leftarrow Y_d/\hat\sigma_d$。

後果是 $\Sigma_Y$（§7）是全樣本共變異數，而 $\|z_{id}\|$ 在高波動日系統性偏大。如 §3 所述，這使得 $Q$ 是 dispersion 加權的。這個選擇是刻意的：**未來報酬的橫截面離散度本身就是 $\mathsf V$ 與 $\mathsf L$ block 要分析的資訊**。若先把它正規化掉，$\mathsf L$ block（大行情中的方向性）就失去了定義的對象。代價是 $Q$ 不對市場 regime 免疫，必須以 $T_{\rm eff}$（定義 3.1）搭配解讀。

---

## §6. Horizon band

### 動機

不同 horizon 的預測問題在經濟上是不同的問題。5 分鐘報酬主要由微觀結構與流動性供給決定；24 小時報酬開始反映資訊擴散與隔夜風險溢酬。把它們混入同一個 whitening 與同一個 SVD，會產生一個由最高變異 horizon 主導、對其他 horizon 無意義的座標系。因此先分 band、band 內獨立分析。

### 定義

把 $H$ 個 horizon 切成 band $H_1,\dots,H_{B_{\rm band}}$。對 band $\beta$：

$$K=|H_\beta|,\qquad Y_d^{(\beta)}=Y_d^{XS}[:,H_\beta].$$

之後所有計算都在 band 內獨立執行。

### 重疊 horizon 與增量基底

這是一個必須明確選擇的規格點。若 $Y_{\cdot,h}$ 為自同一時點起算的**累積**報酬，則不同 horizon 天然巢狀：在 random walk 假設下

$$\operatorname{corr}(r_{0\to1h},r_{0\to2h})=\sqrt{1/2}\approx0.71,\qquad\operatorname{corr}(r_{0\to1h},r_{0\to4h})=\sqrt{1/4}=0.5.$$

$K=4$ 的累積 band 其相關矩陣的條件數約為 10–20，而 whitening 會把最小特徵方向放大 $\lambda_{\min}^{-1/2}$ 倍，把估計雜訊一併放大。

> **規格 6.1**：band 內優先使用**不重疊增量報酬**
> $$r_{0\to1h},\quad r_{1h\to2h},\quad r_{2h\to4h},\quad r_{4h\to8h},$$
> 其相關矩陣在 random walk 下為單位矩陣，條件數接近 1。若研究語意必須使用累積報酬（例如交易執行是持有到 $4h$ 平倉），則保留累積形式，但 §7 的秩控制成為必要而非選配。
>
> 無論選哪一種，band 的定義（累積或增量）必須寫入設定並列入報表，因為 $c_p$ 的座標語意完全依賴於此。

### 維度感

本書統一使用以下 band 切分作為範例：

$$H_1=\{5m,10m,15m,30m\},\quad H_2=\{1h,2h,4h\},\quad H_3=\{6h,12h,24h\},$$

| band | $K$ | $d_{\mathsf D}$ | $d_{\mathsf V}=K(K{+}1)/2$ | $d_{\mathsf L}$ | $d_\Phi$ |
|---|---:|---:|---:|---:|---:|
| $H_1$ | 4 | 4 | 10 | 4 | 18 |
| $H_2$ | 3 | 3 | 6 | 3 | 12 |
| $H_3$ | 3 | 3 | 6 | 3 | 12 |

請記住這張表。$d_{\mathsf V}$ 隨 $K$ 二次成長而 $d_{\mathsf D}$ 只線性成長：$K=4$ 時 Activity block 已佔全部座標的 $10/18=56\%$。這個事實在 §19（fingerprint 的維度基準）與 §38（coverage 的秩上界）兩處產生決定性後果。

---

## §7. Label whitening 與秩控制

### 動機

whitening 的目的是讓不同 horizon 的座標變成同尺度、互不相關，使「feature 預測了哪個 horizon 組合」這個問題有唯一且尺度無關的答案。若不 whitening，$\Sigma_Y$ 中變異最大的 horizon（通常是最長的）會主導所有內積。

但 whitening 同時是整套方法中**唯一會放大雜訊**的步驟：沿著 $\Sigma_Y$ 的第 $j$ 個特徵方向，資料被乘以 $\lambda_j^{-1/2}$。若 $\lambda_K/\lambda_1=10^{-2}$，則最弱方向被放大 10 倍，而那個方向上的估計誤差也被放大 10 倍。因此 whitening 必須搭配秩控制，且秩控制的門檻必須由統計穩定性決定，不能是一個純數值下限。

### §7.1 共變異數與譜

$$\boxed{\Sigma_Y=\frac1T\sum_{d=1}^T\frac1{N_d}Y_d^{(\beta)\top}Y_d^{(\beta)}\in\mathbb R^{K\times K}.}$$

（由 §5 的每日去均值，這確實是共變異數而非未中心化動差。）特徵分解

$$\Sigma_Y=E\Lambda E^\top,\qquad \lambda_1\ge\lambda_2\ge\cdots\ge\lambda_K>0.$$

### §7.2 保留維度 $K_{\rm eff}$

對日期做 moving block bootstrap（區塊長度 $\ge H_{\max}$，$B_\lambda=200$ 次），得到 $\{\lambda_j^{(b)}\}$。定義

$$\boxed{K_{\rm eff}=\max\left\{k:\ \frac{\lambda_k}{\lambda_1}\ge\frac1{\kappa_{\max}}\ \ \text{且}\ \ \frac{\operatorname{sd}_b(\lambda_k^{(b)})}{\lambda_k}\le\tau_\lambda\right\},}$$

預設 $\kappa_{\max}=100$（即最大放大倍率 10 倍）、$\tau_\lambda=0.25$。

兩個條件缺一不可：第一個控制數值放大，第二個控制統計不穩定。一個 $\lambda_k/\lambda_1=0.05$ 但 bootstrap 相對標準差高達 60% 的方向，其 whitened 座標幾乎全是估計噪音，即使它通過了條件數門檻。

**報表必須包含**：完整的 $\{\lambda_j/\lambda_1\}$、$\kappa=\lambda_1/\lambda_{K_{\rm eff}}$、$K_{\rm eff}$、以及被丟棄方向上的 horizon 載重 $E[:,K_{\rm eff}{+}1:K]$（用以判斷丟掉的是什麼）。

若 $K_{\rm eff}<K$ 且被丟棄的方向載重集中在某一兩個 horizon，正確的處理往往不是接受截斷，而是回到 §6 改用增量基底或合併 horizon。截斷是最後手段。

### §7.3 Whitening 矩陣

**情形一（無截斷，$K_{\rm eff}=K$）**：使用對稱平方根倒數

$$B_Y=E\Lambda^{-1/2}E^\top=\Sigma_Y^{-1/2}\in\mathbb R^{K\times K}.$$

選擇對稱平方根（而非 Cholesky）的好處是**不依賴 horizon 的排列順序**，因此 $Z$ 的第 $k$ 個座標仍大致對應第 $k$ 個 horizon，解釋性好。Cholesky whitening 會讓座標變成「第 1 個 horizon」「第 2 個 horizon 扣除第 1 個」這類三角結構，跨 band 不可比。

**情形二（有截斷，$K_{\rm eff}<K$）**：使用降維形式

$$B_Y=E_{1:K_{\rm eff}}\Lambda_{1:K_{\rm eff}}^{-1/2}\in\mathbb R^{K\times K_{\rm eff}}.$$

此時座標的語意變成「未來報酬的前 $K_{\rm eff}$ 個主成分方向」，不再對應個別 horizon。這個語意改變必須在報表中明說，因為 $c_p$ 的每一個分量將依 $E$ 的載重解讀。

兩種情形都定義

$$\boxed{Z_d=Y_d^{(\beta)}B_Y.}$$

### §7.4 命題 7.1（精確白化）

$$Z^\top WZ=I_{K_{\rm eff}}\quad\text{（精確成立，非近似）}.$$

**證明**：$Z^\top WZ=B_Y^\top\big(\frac1T\sum_d\frac1{N_d}Y_d^\top Y_d\big)B_Y=B_Y^\top\Sigma_YB_Y$。情形二中 $B_Y^\top\Sigma_YB_Y=\Lambda^{-1/2}E^\top E\Lambda E^\top E\Lambda^{-1/2}$（取前 $K_{\rm eff}$ 欄）$=\Lambda^{-1/2}\Lambda\Lambda^{-1/2}=I$。情形一為其特例。∎

這個命題是 §9 的基礎：正因為白化是精確的，$\Phi_{\mathsf D}$ 可以直接取為 $Z$，不需要任何額外的正交化步驟，也就不會引入任何座標歧義。

此外由 §5 的每日去均值有 $\sum_iZ_{id,:}=0$。

> **約定**：§7 之後，符號 $K$ 一律指 $K_{\rm eff}$，即實際保留的 label 維度。

---

## §8. 三個 predictive mechanism block

### 總體動機

給定 whitened 未來報酬向量 $z\in\mathbb R^K$，我們要問：一個 cross-sectional 訊號可能與未來報酬的**哪些函數**相關？

原則上答案是「任何 $z$ 的函數」，但那是無限維空間。本方法選擇三個具體的 lifting map，理由是它們對應三個在交易上意義截然不同、需要完全不同應對方式的預測型態：

| block | lifting | 交易語意 | 若模型只學到這個會怎樣 |
|---|---|---|---|
| $\mathsf D$ | $z$ | 誰會漲得比誰多 | 可直接做多空組合 |
| $\mathsf V$ | $\operatorname{svec}(zz^\top)$ | 誰會動得比誰劇烈、horizon 之間如何連動 | 只能做波動度／gamma 相關部位，不能做方向 |
| $\mathsf L$ | $\|z\|z$ | 在大行情中誰的方向對 | 方向性但高度條件化，部位與風險管理完全不同 |

這張表是整套方法存在的理由：三者在傳統 rank IC 下會互相污染，在這個座標系中被嚴格分離。

### §8.1 三個 lifting 與 within-date 中心化

對每個 observation 的 $z_{id}\in\mathbb R^K$：

$$G_{\mathsf D}(z)=z,\qquad G_{\mathsf V}(z)=\operatorname{svec}(zz^\top),\qquad G_{\mathsf L}(z)=\|z\|_2\,z.$$

維度分別為 $K$、$K(K{+}1)/2$、$K$。

隨後，**三者一律做 within-date 中心化**：

$$\boxed{G_{\mathsf B,d}\leftarrow G_{\mathsf B,d}-\mathbf 1_{N_d}\,\overline{G_{\mathsf B,d}},\qquad \overline{G_{\mathsf B,d}}=\frac1{N_d}\sum_iG_{\mathsf B}(z_{id}).}$$

$G_{\mathsf D}$ 本來就滿足（因 $\sum_iz_{id}=0$），$G_{\mathsf V}$ 與 $G_{\mathsf L}$ 不滿足，必須明確執行。

**為什麼一定要做。** 由引理 4.1，這些「當期常數」成分對 $C$ 的貢獻恆為零，因此不做中心化**不會給出錯誤的 $C$**。但它會造成兩個實質損失：

1. **浪費基底維度。** $\operatorname{orth}_W$ 作用在含有日間共同成分的矩陣上，會產生若干條主要由「市場波動度時間序列」構成的 $\Phi_{\mathsf V}$ 方向。這些方向與任何去均值的 $X$ 正交，永遠貢獻 $Q=0$，卻佔用 $d_{\mathsf V}$ 的名額，稀釋 §19 的 fingerprint 分母，並在 §12 的 SVD 中製造無意義的小奇異值。
2. **扭曲正交化的尺度。** $\Xi_R=R^\top WR$ 被日間共同變異主導，$\Lambda_R^{-1/2}$ 因而把真正的 cross-sectional 變異壓縮到不恰當的尺度，使 $\mathsf D\to\mathsf V\to\mathsf L$ 的 Gram–Schmidt 在數值上更不穩定。

**中心化讓定義變簡潔。** 注意原始的 Activity lifting 常寫成 $\operatorname{svec}(zz^\top-I_K)$。中心化之後 $-I_K$ 是一個全域常數，會被當期均值一併吸收，因此

$$\text{中心化}\big(\operatorname{svec}(zz^\top-I_K)\big)=\text{中心化}\big(\operatorname{svec}(zz^\top)\big).$$

本書因此直接寫 $G_{\mathsf V}(z)=\operatorname{svec}(zz^\top)$。

**語意宣告。** 經過這一步，整個 $\Phi$ 生活在 **within-date cross-sectional 函數空間**中。這是對的：本方法從頭到尾只關心橫截面對比，任何「所有股票共同發生的事」都不在研究對象之內。

### §8.2 svec 的定義與驗證

`svec` 把對稱矩陣攤平為向量：對角元素保持不變，非對角元素乘 $\sqrt2$，只取上三角。

**性質**：對任意對稱 $A,B$，

$$\langle\operatorname{svec}A,\operatorname{svec}B\rangle=\operatorname{tr}(A^\top B).$$

**驗證**：$\operatorname{tr}(A^\top B)=\sum_hA_{hh}B_{hh}+2\sum_{h<k}A_{hk}B_{hk}$，而 $\langle\operatorname{svec}A,\operatorname{svec}B\rangle=\sum_hA_{hh}B_{hh}+\sum_{h<k}(\sqrt2A_{hk})(\sqrt2B_{hk})$，兩者相等。∎

這個性質保證 $\mathsf V$ block 的歐氏幾何與矩陣的 Frobenius 幾何一致，因此 $\|c_p^{\mathsf V}\|^2$ 可以直接解讀為「feature 與二階結構的關聯強度」而不含任何座標人為因素。

### §8.3 $\mathsf V$ block 的內容

**實例（$K=2$，horizon 為 1h 與 4h 的增量）**：$G_{\mathsf V}$ 的三個座標為 $(z_1^2,\ \sqrt2z_1z_2,\ z_2^2)$，中心化後為各自減去當期均值。

- 第一與第三座標：該股票在 1h／4h 的**幅度**是否高於當期橫截面典型值。這是「activity」。
- 第二座標：該股票 1h 與 4h 報酬的**同向性**。若 1h 漲而 4h 也漲，$z_1z_2>0$；若 1h 漲後 4h 回吐，$z_1z_2<0$。

因此第二座標其實在測量**趨勢延續 vs 均值回歸**——一個預測它的 feature，實質上是在預測「這支股票未來的短期報酬結構會是延續型還是反轉型」。這是很有價值的資訊，且在任何單一 horizon 的 IC 中完全看不到。這一維也是最常被誤讀的：它不是波動度，而是**跨 horizon 的符號一致性**。

### §8.4 $\mathsf L$ block 的內容

$G_{\mathsf L}(z)=\|z\|z$ 與 $G_{\mathsf D}(z)=z$ 顯然高度相關（同向，只是被 $\|z\|$ 重新加權）。$\mathsf L$ block 的實際內容是 $G_{\mathsf L}$ 在投影掉 $\mathsf D$ 與 $\mathsf V$ 之後的**殘差**，直觀上約為

$$R_{\mathsf L}\approx\big(\|z\|-\mathbb E\|z\|\big)\,z\ \text{的正交化版本},$$

即「**幅度意外 × 方向**」。

一個在 $\mathsf L$ block 上載重的 feature，其含意是：它在平靜時段沒有方向性預測力，但在會出現大幅移動的名次上能正確指出方向。這類 feature 在組合中的角色是**條件性加碼訊號**，而非基礎排序訊號。實務上它的部位設計完全不同：不應該對它做等權多空，而應該把它用作既有方向性部位的規模調節器。

---

## §9. 階層正交化

### 動機

三個 lifting 的像空間彼此重疊。若直接使用，$Q^{\mathsf D}+Q^{\mathsf V}+Q^{\mathsf L}$ 會重複計算同一份資訊，fingerprint 也失去意義。階層正交化解決這個問題，代價是引入**順序依賴**。

> **原則 9.1（歸屬順序）**：固定順序 $\mathsf D\to\mathsf V\to\mathsf L$ 意味著「能被一階方向性解釋的，一律歸給 $\mathsf D$」。因此 $Q^{\mathsf V}$ 的正確讀法**不是**「這個 feature 有多少二階資訊」，而是「這個 feature 有多少**無法被一階方向性解釋的**二階資訊」。

這個區別在解讀報表時至關重要。一個同時預測方向與波動的 feature，其 $Q^{\mathsf V}$ 會偏低，不是因為它不預測波動，而是因為它的波動預測力與方向預測力共線。順序之所以這樣定，是因為方向性資訊在交易上最直接可用，把共用成分歸給它，可以確保 $\mathsf V$ 與 $\mathsf L$ 上的載重一定是「額外的、需要不同處理方式的」資訊。

### §9.1 加權正交化算子

給定矩陣 $R$（列為 observation），計算 $\Xi_R=R^\top WR$，特徵分解 $\Xi_R=U_R\Lambda_RU_R^\top$，保留 $\lambda_j>\epsilon_{\rm rank}\lambda_{\max}$（$\epsilon_{\rm rank}=10^{-10}$，純數值用途），定義

$$\Theta_R=U_R\Lambda_R^{-1/2},\qquad \operatorname{orth}_W(R)=R\,\Theta_R.$$

**驗證**：$(R\Theta_R)^\top W(R\Theta_R)=\Lambda_R^{-1/2}U_R^\top\Xi_RU_R\Lambda_R^{-1/2}=I$。∎

若有維度被 $\epsilon_{\rm rank}$ 截掉，必須記錄並回報，因為那代表 lifting 的像空間實際秩低於名目維度（例如 $K$ 很大而 $T$ 很小時可能發生）。

### §9.2 $\mathsf D$ block：直接取 $Z$

$$\boxed{\Phi_{\mathsf D}:=Z.}$$

**不做任何正交化。** 由命題 7.1，$Z^\top WZ=I_K$ 已經精確成立，因此 $\operatorname{orth}_W$ 是多餘的；而它不只是多餘，還有害。理由如下。

若對 $Z$ 呼叫 $\operatorname{orth}_W$，則 $\Xi_R=I_K$，其特徵分解為 $I=U_RIU_R^\top$，而**任何正交矩陣 $U_R$ 都是合法的特徵向量矩陣**。於是 $\Theta_R=U_R$，得到

$$\Phi_{\mathsf D}=ZU_R,$$

其中 $U_R$ 由 LAPACK 的實作細節、執行緒數、甚至記憶體對齊決定。這對 $Q^{\mathsf D}=\|c_p^{\mathsf D}\|^2$ 與任何夾角都無害（正交變換不改變範數與內積），但對 §21 的 $A^{\rm tgt}_p=(c_p^{\mathsf D})^\top r_{\mathsf D}$ 是致命的：若 $r_{\mathsf D}$ 在一個座標系中定義而 $c_p^{\mathsf D}$ 在另一個，兩者相差一個未知旋轉，alignment 的數值會隨執行環境而變。

令 $\Phi_{\mathsf D}:=Z$ 一勞永逸地消除這個歧義。$\mathsf D$ block 的座標因此有明確的經濟解釋：第 $k$ 個座標就是「未來報酬向量在第 $k$ 個白化方向上的分量」，在無截斷時大致對應第 $k$ 個 horizon。

### §9.3 $\mathsf V$ 與 $\mathsf L$ block

$$\Gamma_{\mathsf{DV}}=\Phi_{\mathsf D}^\top WG_{\mathsf V},\qquad R_{\mathsf V}=G_{\mathsf V}-\Phi_{\mathsf D}\Gamma_{\mathsf{DV}},\qquad\boxed{\Phi_{\mathsf V}=\operatorname{orth}_W(R_{\mathsf V}).}$$

$$\Gamma_{\mathsf{DL}}=\Phi_{\mathsf D}^\top WG_{\mathsf L},\quad\Gamma_{\mathsf{VL}}=\Phi_{\mathsf V}^\top WG_{\mathsf L},\quad R_{\mathsf L}=G_{\mathsf L}-\Phi_{\mathsf D}\Gamma_{\mathsf{DL}}-\Phi_{\mathsf V}\Gamma_{\mathsf{VL}},$$

$$\boxed{\Phi_{\mathsf L}=\operatorname{orth}_W(R_{\mathsf L}).}$$

最後

$$\boxed{\Phi=[\Phi_{\mathsf D}\,|\,\Phi_{\mathsf V}\,|\,\Phi_{\mathsf L}],\qquad \Phi^\top W\Phi=I_{d_\Phi}.}$$

### §9.4 命題 9.2（區塊內旋轉不變性）

$\Phi_{\mathsf V}$ 與 $\Phi_{\mathsf L}$ 內部仍有旋轉自由度（$\operatorname{orth}_W$ 的特徵向量順序與符號）。這不造成任何問題，因為：

> **命題 9.2**：令 $O=\operatorname{diag}(I_{d_{\mathsf D}},O_{\mathsf V},O_{\mathsf L})$ 為任意區塊對角正交矩陣。以 $\Phi O$ 取代 $\Phi$ 後，
> $$C\mapsto CO,\quad V_r\mapsto O^\top V_r,\quad c_s\mapsto O^\top c_s,$$
> 而下列量**完全不變**：$Q_s$、$Q_s^{\mathsf D}$、$Q_s^{\mathsf V}$、$Q_s^{\mathsf L}$、$f_s$、$\tilde f_s$、$S^{\rm mech}_{pq}$、$A^{\rm mech}_{pq}$、$A^{\rm tgt}_s$、$\eta^{\mathsf D}_s$、$N^{\rm eff}$、$R^{\rm canc}$、$r_{\rm sig}$、$r_{\rm eff}$、所有 $\pi$ 與 $M$。

**證明**：$C\mapsto CO$ 直接由定義得出。SVD 為 $CO=US(O^\top V)^\top$，故奇異值不變、$V\mapsto O^\top V$，因而 $r_{\rm sig}$ 與 $r_{\rm eff}$ 不變，$C_{\rm sig}=CV_rV_r^\top\mapsto CO\cdot O^\top V_rV_r^\top O=C_{\rm sig}O$，故 $c_s\mapsto O^\top c_s$，範數與任意兩者的內積皆不變。由於 $O$ 在 $\mathsf D$ 區塊上是單位矩陣，$c_s^{\mathsf D}$ **逐分量不變**，故 $A^{\rm tgt}$ 與 $\eta^{\mathsf D}$ 不變。各 block 的能量 $\|c_s^{\mathsf B}\|^2$ 因 $O$ 為區塊對角而各自保持。$g_{id}=X_{id}\langle\Phi_{id},\hat c\rangle$ 中 $\langle\Phi_{id}O,O^\top\hat c\rangle=\langle\Phi_{id},\hat c\rangle$，故 breadth 與 cancellation 不變。∎

這個命題是整套方法「良定義」的正式陳述：**所有被報告的量都只依賴於 $\Phi$ 的三個 block 各自張成的子空間，不依賴 block 內的基底選擇；而 $\mathsf D$ block 的基底被 $\Phi_{\mathsf D}=Z$ 固定，因此連帶把唯一有方向語意的量（target alignment）也固定了。**

實作上這給了一個免費的單元測試：以不同的 BLAS 執行緒數、不同的隨機種子重跑整個管線，上述所有量必須在數值誤差內完全相同（§65）。

---

## §10. Fit / transform 邊界

### 動機

這是把「研究」與「回測」分開的關鍵機制。若 label basis 在每個評估期重新擬合，則不同期間的 $c_p$ 落在不同座標系中，任何跨期比較、任何 feature 與 alpha 的比較都失去意義。

### 規格

在 training / research period 執行 `fit_label_basis(Y_train)`，儲存

$$B_Y;\qquad \Theta_{\mathsf V},\ \Theta_{\mathsf L};\qquad \Gamma_{\mathsf{DV}},\ \Gamma_{\mathsf{DL}},\ \Gamma_{\mathsf{VL}}.$$

新資料 $Y_{\rm new}$ 到來後：

$$Z_{\rm new}=Y_{\rm new}B_Y,\qquad \Phi_{\mathsf D}^{\rm new}=Z_{\rm new},$$

$$G_{\mathsf V}^{\rm new},G_{\mathsf L}^{\rm new}\ \text{依 §8 建構並 within-date 中心化},$$

$$R_{\mathsf V}^{\rm new}=G_{\mathsf V}^{\rm new}-\Phi_{\mathsf D}^{\rm new}\Gamma_{\mathsf{DV}},\qquad \Phi_{\mathsf V}^{\rm new}=R_{\mathsf V}^{\rm new}\Theta_{\mathsf V},$$

$$R_{\mathsf L}^{\rm new}=G_{\mathsf L}^{\rm new}-\Phi_{\mathsf D}^{\rm new}\Gamma_{\mathsf{DL}}-\Phi_{\mathsf V}^{\rm new}\Gamma_{\mathsf{VL}},\qquad \Phi_{\mathsf L}^{\rm new}=R_{\mathsf L}^{\rm new}\Theta_{\mathsf L}.$$

注意 within-date 中心化**不是擬合參數**，它在每個日期上就地執行，因此不構成洩漏。

### 關鍵性質：$\delta_{\rm drift}$

在 OOS 期間 $\Phi^{\rm new}$ 不再滿足 $\Phi^\top W\Phi=I$。這不是錯誤，而是一個免費的 regime-drift 偵測器：

$$\boxed{\delta_{\rm drift}=\big\|\Phi_{\rm OOS}^\top W_{\rm OOS}\Phi_{\rm OOS}-I\big\|_F,}$$

並額外拆成三個對角區塊的 $\delta^{\mathsf D},\delta^{\mathsf V},\delta^{\mathsf L}$。

若 $\delta_{\rm drift}$ 小，代表未來報酬的共變異結構與研究期一致，所有以 $c$ 表達的量可以直接跨期比較。若 $\delta_{\rm drift}$ 大（例如 $>0.3$），代表市場結構已變，此時 $Q_\alpha$ 的下降**可能來自座標系失效而非模型失效**——這是解讀任何 OOS 衰減時必須先排除的替代解釋。

### 實例

研究期 2022–2023，OOS 為 2024 上半年。若 2024 上半年市場波動度為研究期的一半，則 $\|z\|$ 系統性偏小；由於 $G_{\mathsf V}\sim\|z\|^2$，

$$\Phi_{\mathsf V}^{\rm OOS\top}W\Phi_{\mathsf V}^{\rm OOS}\approx0.25\,I\ \Longrightarrow\ \delta^{\mathsf V}\approx0.75\sqrt{d_{\mathsf V}}\approx2.4\ (d_{\mathsf V}=10).$$

此時 $Q_\alpha^{\mathsf V}$ 會機械性地下降約 4 倍，**不可解讀為模型的 activity 預測力衰退**。正確的處理是同時報告 $\delta^{\mathsf V}$，並在必要時以 $\delta$ 校正後的量做跨期比較。

---
# 第三部分：Feature–Label 幾何

## §11. Predictive operator

### 定義

對每一期

$$\boxed{C^{(d)}=\frac1{N_d}X_d^{XS\top}\Phi_d\ \in\mathbb R^{P\times d_\Phi},}$$

全期間（逐 feature 在其有效日集合 $\mathcal D_p$ 上平均，見 §4.2）

$$\boxed{C_{p,:}=\frac1{T_p}\sum_{d\in\mathcal D_p}C^{(d)}_{p,:}.}$$

當所有 feature 覆蓋率皆為 100% 時，這簡化為 $C=\frac1T\sum_dC^{(d)}=X^\top W\Phi$。

$C$ 是 **average daily cross-sectional predictive operator**。第 $p$ 列 $C_{p,:}$ 表示 feature $p$ 與所有 future-return mechanism 座標的平均橫截面關聯。

### 量級：$C$ 的每個元素就是一個 IC

由 §4，$X^{XS}_{\cdot p}$ 每期均值 0、變異 1；由 §9，$\Phi$ 的每一欄在加權內積下範數為 1。因此

$$C_{p,k}\approx\overline{\operatorname{corr}_d\big(X_{\cdot dp},\Phi_{\cdot dk}\big)}$$

即 feature $p$ 與第 $k$ 個機制座標的**平均每期 cross-sectional 相關係數**，也就是一個 IC。於是

$$Q_p=\|c_p\|^2\approx\sum_k\overline{IC}_{p,k}^2.$$

這給出判讀 $Q$ 的尺規：

| 平均每期 IC | 單一座標的 $Q$ 貢獻 | 在股票市場的意義 |
|---:|---:|---|
| 0.005 | $2.5\times10^{-5}$ | 幾乎不可偵測 |
| 0.01 | $1\times10^{-4}$ | 弱，但在高頻可交易 |
| 0.03 | $9\times10^{-4}$ | 良好的單一 feature |
| 0.05 | $2.5\times10^{-3}$ | 很強 |
| 0.10 | $1\times10^{-2}$ | 極為罕見，先懷疑資料洩漏 |
| 0.20 | $4\times10^{-2}$ | 幾乎必然是 look-ahead bias |

**本書所有數值例子都落在 $Q^{\mathsf D}\in[10^{-5},10^{-3}]$ 這個帶內。** 把這個尺規記錯一到兩個數量級，會導致所有門檻設錯，進而誤刪全部真實 feature 或誤留全部雜訊。

### 規模

$P=400$、$d_\Phi=18$ 時，$C$ 是 $400\times18$ 的矩陣，約 7200 個數字。整份 feature 研究的全部原始素材就是這 7200 個數字加上它們的虛無分布。

---

## §12. 奇異值分解

### 定義

$$\boxed{C=USV^\top,}\qquad S=\operatorname{diag}(s_1,\dots),\ s_1\ge s_2\ge\cdots$$

$V$ 的欄是 label-mechanism 方向；$U$ 的欄是 feature-side 組合。

### 解讀

$C$ 的 SVD 回答的問題是：**在 $d_\Phi$ 維的 label mechanism 空間中，整個 feature library 集體地「照亮」了哪些方向？**

- $v_1$：整個 library 最有共識的預測方向。在多數股票資料上，這會是一個以 $\mathsf D$ block 為主、接近「所有 horizon 同向」的方向，即動能／反轉這個古典軸。
- $u_1$：達成這個預測的 feature 線性組合，可視為自動萃取的 super-feature。
- $s_1$：該方向的強度。

**注意 $U$ 與 $V$ 的不對稱性。** $V$ 生活在 label 空間（維度 $d_\Phi$，小），$U$ 生活在 feature 空間（維度 $P$，大）。因此

$$r\le\min(P,d_\Phi)=d_\Phi.$$

**可偵測的 mode 數上限由 label 空間維度決定，與 feature 數無關。** $K=4$ 時 $r\le18$，實務上 $r_{\rm sig}$ 幾乎總是落在 2–6。這個事實在 §38 會產生決定性後果。

---

## §13. Permutation null

### 動機

$Q_p$ 需要一個參考分布。它不能由解析式給出，因為 $C$ 的虛無分布同時取決於 feature 的橫截面分布、$\Phi$ 的橫截面分布、以及兩者各自的時間自相關結構。唯一可靠的做法是建構一個保留這些結構、只破壞「哪一天配哪一天」的重抽方案。

### §13.1 虛無假設

原始配對是 $X_d\leftrightarrow\Phi_d$。虛無假設為：

> feature 的橫截面排序，與**正確時間**的 future-return mechanism 沒有預測關係。

因此破壞日期對齊：$X_{\pi(d)}\leftrightarrow\Phi_d$，但保留

1. 每一期 $X$ 的完整橫截面結構（含產業、市值等群聚）；
2. 每一期 $\Phi$ 的完整橫截面結構；
3. ticker identity；
4. $X$ 自身的時間自相關；
5. label 自身的時間自相關。

這五項保留使得 null 不會因為「打亂了橫截面結構」而人為變小——那是最常見的、會讓一切看起來都顯著的錯誤。

### §13.2 主要方案：circular date shift

給定 ordered dates，選 shift $s$，令

$$\pi_s(k)=1+\big((k+s-1)\bmod T\big).$$

合法 shift 集合為

$$\mathcal S=\{s:\Delta\le\min(s,T-s)\},\qquad\boxed{\Delta=H_{\max}+\ell_X,}$$

其中 $H_{\max}$ 為 band 內最長 horizon 的長度（以期數計），$\ell_X$ 為 feature 的記憶長度（例如 20 日移動平均則 $\ell_X=20$）。這個最小位移確保被配對的 $X$ 與 $Y$ 在時間上沒有重疊窗口。**$\Delta$ 與 §41 的 purge gap 必須是同一個數字**，否則同一份資料在兩處被以不同的洩漏標準對待。

實務上取 $|\mathcal S|$ 個 shift 中的 $B_{\rm null}$ 個（預設全取，上限 2000）。

### §13.3 Universe 不一致：matched thinning 規則

$N_d$ 隨時間變動，因此 $\mathcal U_{\pi_s(d)}\ne\mathcal U_d$。必須同時處理三件事。

**(a) 取交集。** 令

$$\mathcal I_{d,s}=\mathcal U_{\pi_s(d)}\cap\mathcal U_d.$$

若 $|\mathcal I_{d,s}|<\rho_{\min}\cdot\min(N_{\pi_s(d)},N_d)$（預設 $\rho_{\min}=0.9$），捨棄該 $(d,s)$ 配對。

**(b) 在交集上重新標準化。** 由引理 4.1，$C^{(d)}$ 是一個真正的共變異數，**前提是 $X$ 在所使用的 universe 上均值為零**。限制到 $\mathcal I_{d,s}$ 之後這個前提被破壞：若某 feature 與「是否在兩個日期都存在」相關（市值、流動性、上市年齡這類 feature 天生如此），則交集上的均值不為零，null 的 $C$ 會被一個非零的均值項系統性污染，使該類 feature 的 null 偏大或偏小。因此**必須在 $\mathcal I_{d,s}$ 上重新去均值與重新標準化 $X$**：

$$\check X_{idp}=\frac{X'_{idp}-\mu_{dp}^{\mathcal I}}{\sigma_{dp}^{\mathcal I}+\epsilon},\qquad i\in\mathcal I_{d,s}.$$

$\Phi$ 不需要重新處理：由引理 4.1，只要 $\check X$ 在 $\mathcal I$ 上均值為零，$\Phi$ 的任何當期常數成分自動被吸收。

**(c) 對觀測值施加相同的抽樣規則。** 觀測統計量若在完整的 $\mathcal U_d$ 上計算，而 null 在較小的交集上計算，則兩者的有效樣本數不同，null 的變異會被高估（偏保守）或低估（若交集有選擇性）。對稱化的規則是：

> **規格 13.1（matched thinning）**：對每個日期 $d$，獨立抽一個 $s_d\sim\text{Uniform}(\mathcal S)$，定義觀測用的分析 universe
> $$\mathcal A_d=\mathcal U_d\cap\mathcal U_{\pi_{s_d}(d)},$$
> 並在 $\mathcal A_d$ 上以 (b) 的方式重新標準化後計算觀測統計量。重複 $R_{\rm thin}=5$ 組獨立抽樣並取平均，作為最終的觀測值。

如此觀測與 null 的 universe 規模分布完全同分布。

**簡化條件**：若 $\operatorname{median}_{d,s}|\mathcal I_{d,s}|/N_d\ge0.98$，thinning 的影響小於數值誤差，可直接在 $\mathcal U_d$ 上計算觀測值。這個比率必須計算並記錄，不可默認。

### §13.4 慢變 feature 與替代虛無

circular shift 的檢定力取決於「位移 $\Delta$ 期之後，feature 的橫截面排序是否已經改變」。對一個記憶長度遠大於 $\Delta$ 的 feature（例如 250 日動能、市值、產業啞變數），

$$X_{d+\Delta,p}\approx X_{d,p},$$

於是 $C^{\rm null}\approx C^{\rm obs}$，檢定幾乎沒有檢定力。**這會使本方法系統性地偏好快速輪動的訊號**，若不處理，慢變但真實的訊號會被全數判為不顯著。

**持續性診斷。** 對每個 feature 定義

$$\boxed{\kappa_p=\operatorname{median}_d\ \operatorname{corr}_{i\in\mathcal U_d\cap\mathcal U_{d+\Delta}}\big(X_{idp},X_{i,d+\Delta,p}\big).}$$

$\kappa_p$ 必須列入報表。經驗規則：$\kappa_p>0.5$ 時 circular shift null 的檢定力已顯著下降，$\kappa_p>0.8$ 時幾乎為零。

**替代虛無（N2）：分層 within-date permutation。** 對 $\kappa_p>0.5$ 的 feature，以

$$X_{idp}\longrightarrow X_{\pi_d(i),dp}$$

取代日期位移，其中 $\pi_d$ 是**在分層內**的隨機排列：

$$\mathcal U_d=\bigsqcup_g\mathcal U_{d,g},\qquad g=\text{sector}\times\text{size bucket},$$

在每個 $g$ 內獨立排列。

**分層是必要的，不是精緻化。** 若在全 universe 內自由排列，會同時摧毀產業與市值的群聚結構，於是一個 sector-level feature 的 null 變得遠小於真實值，檢定變得嚴重 anti-conservative。分層固定了這些已知暴險的橫截面組成，使 null 只破壞「個股層級的配對」。

更乾淨的變體是先殘差化已知暴險：令 $F$ 為已知因子暴險矩陣（產業、市值、beta），

$$X^\perp_{d,p}=(I-P_{F_d})X_{d,p},\qquad P_{F_d}=F_d(F_d^\top F_d)^{-1}F_d^\top,$$

再對 $X^\perp$ 做自由 within-date permutation。此時虛無假設變成「在控制已知暴險後，個股層級的 feature 排序沒有預測力」，語意更精確，代價是結論只涵蓋殘差部分。

**報表規則**：對 $\kappa_p>0.5$ 的 feature，主要結論使用 N2，同時列出 N1 的結果與 $\kappa_p$；對 $\kappa_p\le0.5$ 的 feature，主要結論使用 N1。兩者所使用的 null 必須在報表中逐 feature 標記。

### §13.5 實作流程

```text
輸入: {X_d}, {Phi_d}, 合法 shift 集合 S, 重抽次數 B_null
預先: 對每個 d 計算 kappa_p (持續性診斷)

# ---- 觀測統計量 (matched thinning) ----
for r in 1..R_thin:
    for d in 1..T:
        s_d ~ Uniform(S)
        A_d = U_d ∩ U_{pi_{s_d}(d)}
        若 |A_d| < rho_min * min(N_d, N_{pi_{s_d}(d)}): 跳過此 d
        X_check = 在 A_d 上重新標準化的 X_d
        C_r += (1/|A_d|) * X_check^T Phi_d[A_d]
    C_obs_r = C_r / (有效日數)
C_obs = mean_r C_obs_r

# ---- N1: circular shift null ----
for b in 1..B_null:
    s = S[b]
    for d in 1..T:
        I = U_{pi_s(d)} ∩ U_d
        若 |I| < rho_min * min(...): 跳過
        X_check = 在 I 上重新標準化的 X_{pi_s(d)}
        C_b += (1/|I|) * X_check^T Phi_d[I]
    C_null[b] = C_b / (有效日數)

# ---- N2: 分層 within-date permutation (僅慢變 feature) ----
for b in 1..B_null:
    for d in 1..T:
        for g in strata(d):
            在 U_{d,g} 內隨機排列 X 的列
        C_b += (1/N_d) * X_perm_d^T Phi_d
    C_null2[b] = C_b / T
```

所有下游的虛無量（$r_{\rm sig}$、feature p-value、相似度門檻、fingerprint 基準、breadth 與 cancellation 的基準）都從同一組 $\{C^{\rm null}_b\}$ 導出，因此這個迴圈只跑一次。

---

## §14. 顯著秩 $r_{\rm sig}$

### 方法：permutation parallel analysis

對每個 $b$ 計算 null 矩陣的**完整奇異值譜** $s_1(C_b^{\rm null})\ge s_2(C_b^{\rm null})\ge\cdots$。對每個 $j$：

$$\boxed{p_j=\frac{1+\#\{b:\ s_j(C_b^{\rm null})\ge s_j(C)\}}{B_{\rm null}+1}.}$$

以 step-down 方式接受：從 $j=1$ 開始，只要 $p_j<\gamma_{\rm mode}$（預設 0.05）就接受並繼續，第一次失敗即停止。

$$\boxed{r_{\rm sig}=\#\text{accepted modes}.}$$

（分子分母各加 1 是標準的 Monte Carlo p 值修正，保證 $p\ge1/(B_{\rm null}+1)>0$，使檢定在有限 $B$ 下仍為有效水準的檢定。）

### 為什麼比較「第 $j$ 個對第 $j$ 個」

一個直覺上自然、但實際上有偏的替代做法是**序貫收縮**（sequential deflation）：接受前 $j-1$ 個 mode 後，令

$$R_j=P_U^\perp CP_V^\perp,\qquad P_U^\perp=I-U_{<j}U_{<j}^\top,\quad P_V^\perp=I-V_{<j}V_{<j}^\top,$$

以 $\|R_j\|_2$ 為觀測統計量，並以 $\|P_U^\perp C_b^{\rm null}P_V^\perp\|_2$ 為 null。

這個做法的問題在於兩側不對稱。對觀測矩陣，$P_U^\perp CP_V^\perp=\sum_{k\ge j}s_ku_kv_k^\top$，因此 $\|R_j\|_2$ **恰好等於** $C$ 的第 $j$ 個奇異值——收縮精確地移除了前 $j-1$ 個主方向。但對 null 矩陣，投影子由觀測資料決定，與 $C_b^{\rm null}$ 自己的主方向無關；從 $P$ 維空間移除 $j-1$ 個隨機取向的方向，對各向同性隨機矩陣的最大奇異值影響極小。因此

$$\|P_U^\perp C_b^{\rm null}P_V^\perp\|_2\approx s_1(C_b^{\rm null})\quad\text{對所有 }j.$$

於是檢定實際上是拿觀測的**第 $j$ 個**奇異值對抗 null 的**第 1 個**。隨機矩陣的奇異值譜相當平坦（典型比值 $s_4/s_1\approx0.8$），所以這個比較隨 $j$ 增大而越來越不公平，$r_{\rm sig}$ 被系統性低估。

具體感受：null 譜為 $(1.00,0.92,0.86,0.81,\dots)$，觀測譜為 $(3.0,1.5,1.05,0.95,\dots)$。序貫收縮給出 $r_{\rm sig}=2$（第三個 $1.05$ vs $1.00$ 邊緣、第四個 $0.95$ vs $1.00$ 失敗）；parallel analysis 則正確地把第四個 $0.95$ 與 null 的第四個 $0.81$ 比較，給出 $r_{\rm sig}=4$。

### 註記

step-down「第一次失敗即停止」不嚴格控制整體 FWER。這在實務上可接受，因為 $r_{\rm sig}$ 是後續分析的一個超參數而非最終結論；但它不應被當作嚴格的多重檢定控制來引用。若 $r_{\rm sig}$ 落在 1 或恰好等於 $d_\Phi$，都應視為警訊：前者代表整個 library 只照亮了一個方向（可能是 library 建設過於同質），後者代表 null 校準可能有誤（$d_\Phi$ 個方向全部顯著在真實資料上極為罕見）。

---

## §15. 顯著子空間與 discovery / evaluation 切分

### §15.1 為什麼需要切分

$V_r$ 是從**觀測資料**中挑出的「最強方向」。若之後又用同一份觀測資料計算

$$Q_p=\|C_{p,:}V_rV_r^\top\|^2,$$

並拿它去對抗「把 null 投影到**同一個** $V_r$」所得的分布，則兩側的處理並不對等：

$$\text{觀測}:\ \text{資料}\to\text{挑最佳子空間}\to Q,\qquad \text{null}:\ \text{null 資料}\to\text{投影到別人的最佳子空間}\to Q.$$

null 少經歷了「挑選」這一步，因此系統性偏小，p-value 偏樂觀。這是一個典型的 selective inference 問題，其後果不是理論上的潔癖：在 $d_\Phi=18$、$r_{\rm sig}=3$ 的設定下，挑選效應可以讓邊緣 feature 的 p-value 低估數倍。

有兩種正確的處理方式。本書以第一種為主要規格。

### §15.2 主要規格：discovery / evaluation 切分

把研究期依時間切成兩段，中間留 $\Delta=H_{\max}+\ell_X$ 的 purge gap：

$$\underbrace{\mathcal D}_{\text{discovery,}\ \approx60\%}\ \ \|\ \Delta\ \|\ \ \underbrace{\mathcal E}_{\text{evaluation,}\ \approx40\%}$$

**在 $\mathcal D$ 上決定座標系與子空間：**

1. 以 §7–§9 擬合 label basis（$B_Y,\Theta,\Gamma$）。
2. 計算 $C_{\mathcal D}$，以 §13 的 null 與 §14 的 parallel analysis 決定 $r_{\rm sig}$ 與

$$\boxed{V_r=[v_1,\dots,v_{r_{\rm sig}}]\quad\text{（由 }\mathcal D\text{ 決定，此後固定）}.}$$

**在 $\mathcal E$ 上做所有推論：**

3. 以 §10 的 transform 把 label basis 套用到 $\mathcal E$，計算 $C_{\mathcal E}$。
4. 定義

$$\boxed{C_{\rm sig}=C_{\mathcal E}V_rV_r^\top,\qquad c_p=C_{\rm sig}[p,:].}$$

5. §13 的 permutation **只在 $\mathcal E$ 上執行**，$V_r$ 保持固定。
6. §26–§29 的 bootstrap、§23 的相似度門檻校準、§34–§39 的所有篩選決策，全部在 $\mathcal E$ 上進行。

**為什麼這樣是有效的。** $V_r$ 只是 $\mathcal E$ 資料之外的一個固定矩陣，因此在 $\mathcal E$ 上「$X$ 與 $\Phi$ 無關聯」這個虛無假設下，$C_{\mathcal E}V_rV_r^\top$ 的 permutation 分布是精確的（至多受限於 $\mathcal D$ 與 $\mathcal E$ 之間殘存的時間相依，而這已由 purge gap 處理）。挑選效應完全被隔離在 $\mathcal D$ 之內。

**附帶好處。** 切分同時產生一個免費的診斷：$\mathcal D$ 與 $\mathcal E$ 之間的 $\delta_{\rm drift}$（§10）。若它很大，代表座標系本身在兩段之間就已經不穩定，此時任何 OOS 結論都應先降級。

**切分比例。** $\mathcal D$ 太短則 $V_r$ 估計不穩、$r_{\rm sig}$ 偏低；$\mathcal E$ 太短則所有 p-value 的解析度下降。$60/40$ 是合理起點。判斷標準是：$\mathcal D$ 必須長到讓 $\Sigma_Y$ 的譜穩定（§7.2 的 bootstrap 已提供此檢查），$\mathcal E$ 必須長到 $T^{\mathsf D}_{\rm eff}(\mathcal E)\ge200$。

### §15.3 替代規格：在每次 permutation 內重跑挑選

當樣本長度不足以切分時，改為：

> 對每一個 null 重抽 $b$，以 $C_b^{\rm null}$ **自行**計算 SVD，取其前 $r_{\rm sig}$ 個右奇異向量 $V_r^{(b)}$，再定義
> $$Q_{p,b}^{\rm null}=\big\|C^{\rm null}_{b,p,:}V_r^{(b)}V_r^{(b)\top}\big\|^2.$$

即每一個 null 都完整重現「挑選最佳子空間」這個步驟。$r_{\rm sig}$ 固定為觀測值（若連 $r_{\rm sig}$ 也重新估計，檢定會變得過度保守且難以解讀）。

計算成本極低：$C_b^{\rm null}$ 本來就要算，而對一個 $400\times18$ 矩陣做 SVD 的成本可以忽略。

兩種規格的差異在於它們校準的對象：切分法把挑選效應完全移出推論；重跑法則把挑選效應納入 null，使兩側對等。兩者都有效，但**不可混用**——同一份報表中所有 p-value 必須來自同一種規格，並在報表標頭註明。

### §15.4 最終重擬合

所有篩選決策鎖定之後，為了生產使用，可以在 $\mathcal D\cup\mathcal E$ 上重新擬合 label basis 與 $V_r$。

> **規則 15.1**：重擬合之後**不得做任何新的決策**。不得重新篩選 feature、不得重新調整門檻、不得因為某個 feature 在重擬合後變得顯著而把它加回來。重擬合的唯一目的是讓生產用的座標系使用全部資料。

違反這條規則會讓整個切分設計失效，且失效的方式無法從報表上看出來。

---

## §16. Feature Quality

### 定義

$$\boxed{Q_p=\|c_p\|^2,}\qquad Q_p^{\mathsf B}=\big\|c_p[\mathsf B]\big\|^2,\ \ \mathsf B\in\{\mathsf D,\mathsf V,\mathsf L\},$$

其中 $c_p[\mathsf B]$ 表示 $c_p$ 落在 block $\mathsf B$ 座標上的子向量。

### 命題 16.1（加總恆等式）

$$Q_p=Q_p^{\mathsf D}+Q_p^{\mathsf V}+Q_p^{\mathsf L}.$$

**證明**：$\{1,\dots,d_\Phi\}$ 被三個 block 分割為不相交的座標子集，向量的平方範數等於各子集平方和。∎

這個恆等式**恆為真**，與 $\Phi$ 是否正交無關。因此 §65 中對應的單元測試只能抓到程式的切片錯誤，不能驗證正交性；真正驗證正交性的是 $\|\Phi^\top W\Phi-I\|_F$ 測試。

### 命題 16.2（$Q$ 是一個加權 $R^2$）

令 $c_p^{\rm full}=X_p^\top W\Phi$（未投影）。則

$$\|c_p^{\rm full}\|^2=\big\|P_\Phi X_p\big\|_W^2,\qquad P_\Phi=\Phi\Phi^\top W,$$

即把 $X_p$ 對 $\Phi$ 的各欄做加權最小平方迴歸所得的 $R^2$（因 $\|X_p\|_W^2=1$）。

**證明**：$P_\Phi$ 是 $W$-內積下對 $\operatorname{span}(\Phi)$ 的正交投影（$P_\Phi^2=P_\Phi$，且 $\langle P_\Phi u,v\rangle_W=\langle u,P_\Phi v\rangle_W$，由 $\Phi^\top W\Phi=I$ 直接驗證）。故 $\|P_\Phi X_p\|_W^2=X_p^\top W\Phi\Phi^\top WX_p=\|c_p^{\rm full}\|^2$。∎

### 推論 16.3（上界）

$$\boxed{0\le Q_p\le\|c_p^{\rm full}\|^2\le\|X_p\|_W^2=1.}$$

第一個不等式因 $P_V$ 是正交投影，第二個因 $P_\Phi$ 是正交投影且 §4.1 的第三步保證 $\|X_p\|_W=1$。

### 解讀

$Q_p$ 的精確語意是：

> **feature $p$ 的橫截面變異中，有多少比例可以被「落在統計顯著機制子空間內的」未來報酬函數線性解釋。**

三個限定詞都不能省略：

- 「橫截面變異」——不含任何時間序列預測力；
- 「線性解釋」——見 §1 的邊界二與 §4.3 的 transform family；
- 「落在顯著子空間內」——$Q_p$ 已經過 $P_V$ 投影，因此一個只預測 $V_r$ 之外方向的 feature 會得到 $Q_p\approx0$。

最後一點有一個必須被意識到的後果：$V_r$ 是由**整個 library 的共識**決定的。一個真實但孤立的機制——只有一個 feature 預測它——不會出現在 $V_r$ 中，於是該 feature 被判為無效。這是本方法為了控制多重檢定而付出的代價。緩解方式是同時報告未投影的 $Q_p^{\rm full}=\|c_p^{\rm full}\|^2$ 與其 permutation p-value；若某 feature 的 $Q^{\rm full}$ 顯著而 $Q$ 不顯著，它就是這種「孤立機制」的候選，應單獨檢視而非直接丟棄。

---

## §17. $Q$ 能回答與不能回答的問題

$Q_p$ 能回答：

- 這個 feature 有多少預測資訊？（量級）
- 這些資訊屬於哪一類機制？（$Q^{\mathsf D}$ vs $Q^{\mathsf V}$ vs $Q^{\mathsf L}$）
- 與其他 feature 相比如何？（因 $\|X_p\|_W=1$，可直接比較）

$Q_p$ **不能**回答：

- 這個 feature 能賺多少錢。$Q$ 沒有考慮交易成本、容量、換手率、借券可得性。
- 這個 feature 在組合中的邊際價值。$Q$ 是單變量的；一個 $Q$ 中等但與其他 feature 正交的訊號，其組合價值可能高於一個 $Q$ 很大但與既有部位重複的訊號。這正是 §38 的 coverage 準則存在的理由。
- 方向對不對。$Q$ 是平方量，不含符號。方向由 §21 的 $A^{\rm tgt}$ 回答。

### 實例

一個 feature 的報表為

$$Q_p=4.0\times10^{-4},\quad Q_p^{\mathsf D}=3.6\times10^{-4},\quad Q_p^{\mathsf V}=3\times10^{-5},\quad Q_p^{\mathsf L}=1\times10^{-5}.$$

解讀：總 IC 量級約 $\sqrt{4\times10^{-4}}=0.02$，其中 90% 的能量在方向性上。這是一個乾淨的方向性訊號，強度中等偏上，可直接進入多空排序。

另一個 feature：

$$Q_q=9.0\times10^{-4},\quad Q_q^{\mathsf D}=5\times10^{-5},\quad Q_q^{\mathsf V}=8.0\times10^{-4},\quad Q_q^{\mathsf L}=5\times10^{-5}.$$

總量級更大（IC 約 0.03），但 89% 的能量在 $\mathsf V$ block。這**不是**一個方向性訊號，而是一個波動度／活躍度預測器。若把它放進多空排序，它會產生一個「做多即將劇烈波動的股票」的部位，其報酬特性接近做多 gamma，與方向性 alpha 完全不同。在只做方向性研究時，這個 feature 應被 §34 的 eligibility 排除（它的 $Q^{\mathsf D}$ 太小），但應保留在 library 中供波動度研究使用。

---

## §18. Feature 層級的顯著性與多重檢定

### §18.1 單一 feature 的 p-value

對每個 feature 與每個 block，permutation p-value 為

$$\boxed{p_p^{\mathsf B}=\frac{1+\#\{b:\ Q^{{\mathsf B},{\rm null}}_{p,b}\ge Q_p^{\mathsf B}\}}{B_{\rm null}+1},}$$

其中 $Q^{\rm null}_{p,b}$ 依 §15 的規格（切分法：$V_r$ 固定；重跑法：$V_r^{(b)}$ 自行重估）計算。同樣的公式套用於 $Q_p$（總量）與 §21 的 $(A^{\rm tgt}_p)^2$。

### §18.2 多重檢定：解析度問題

$P=400$ 個 feature 同時檢定，必須做多重檢定控制。這裡有一個容易被忽略的機械限制：

$$\min_b p=\frac1{B_{\rm null}+1}.$$

若採用 Benjamini–Hochberg 在 $\gamma_{\rm FDR}=0.05$ 下控制，最強的 feature 需要 $p\le\gamma_{\rm FDR}/P=1.25\times10^{-4}$，因此需要

$$B_{\rm null}\ge\frac{P}{\gamma_{\rm FDR}}-1=7999.$$

若只跑 $B_{\rm null}=1000$，則**即使是全世界最強的 feature 也不可能通過 BH**——不是因為它不顯著，而是因為 p-value 的解析度不夠。這個失敗模式在報表上表現為「所有 feature 都不顯著」，非常容易被誤診為訊號不存在。

### §18.3 主要方案：Westfall–Young max-T step-down

Westfall–Young 程序控制的是 FWER，且**不需要**解析度到 $1/P$ 的量級，因為它比較的是「觀測統計量」與「null 的跨 feature 最大值」，而後者本身就是一個極值分布。

**第一步：studentize。** 各 feature 的 null 尺度不完全相同（覆蓋率、持續性 $\kappa_p$、橫截面分布形狀都會影響）。因此先標準化：

$$t_p=\frac{Q_p^{\mathsf B}-\operatorname{mean}_b Q^{\rm null}_{p,b}}{\operatorname{sd}_bQ^{\rm null}_{p,b}+\epsilon},\qquad t_{p,b}=\frac{Q^{\rm null}_{p,b}-\operatorname{mean}_{b'}Q^{\rm null}_{p,b'}}{\operatorname{sd}_{b'}Q^{\rm null}_{p,b'}+\epsilon}.$$

**第二步：step-down 最大值。** 將觀測值排序 $t_{(1)}\ge t_{(2)}\ge\cdots\ge t_{(P)}$，對每個重抽 $b$ 計算後綴最大值

$$m_{b,k}=\max_{j\ge k}t_{(j),b}.$$

**第三步：調整後 p-value。**

$$\tilde p_{(k)}=\frac{1+\#\{b:\ m_{b,k}\ge t_{(k)}\}}{B_{\rm null}+1},\qquad \tilde p_{(k)}\leftarrow\max\big(\tilde p_{(k-1)},\tilde p_{(k)}\big).$$

最後一步強制單調性，是 step-down 程序的必要成分。

**為什麼這解決了解析度問題**：最強 feature 的 $\tilde p$ 可以低到 $1/(B_{\rm null}+1)$，而判準是 $\tilde p<\gamma_{\rm FWER}=0.05$。因此 $B_{\rm null}=2000$ 就綽綽有餘。

**為什麼 Westfall–Young 比 BH 更適合這個問題**：BH 假設檢定統計量之間為正相依或獨立。本問題中 $Q_p$ 之間高度相依（所有 feature 共用同一個 $\Phi$、同一批日期、且 feature 彼此高度相關），相依結構未知且複雜。Westfall–Young 透過 permutation **直接重現**這個相依結構，不需要任何假設。代價是它控制 FWER 而非 FDR，較為嚴格。

### §18.4 次要方案：permutation FDR

若 library 很大而希望較寬鬆的篩選，可用 permutation 估計的 FDR：在門檻 $\tau$ 下，

$$\widehat{\rm FDR}(\tau)=\frac{\frac1{B_{\rm null}}\sum_b\#\{p:t_{p,b}\ge\tau\}}{\#\{p:t_p\ge\tau\}},$$

取滿足 $\widehat{\rm FDR}(\tau)\le\gamma_{\rm FDR}$ 的最小 $\tau$。這同樣完全以 permutation 為基礎，不依賴相依性假設。

**報表規則**：同時列出 $\tilde p^{\rm WY}_p$ 與 $\widehat{\rm FDR}$ 門檻下的選取結果。主要篩選（§34）使用 Westfall–Young；FDR 版本作為敏感度分析。若兩者給出的 feature 數差距超過三倍，代表訊號集中在少數幾個極強的 feature 上，其餘皆為邊緣——這個結構本身是重要資訊。

---

## §19. Mechanism Fingerprint

### §19.1 Raw fingerprint

$$\boxed{f_p=\left(\frac{Q_p^{\mathsf D}}{Q_p},\ \frac{Q_p^{\mathsf V}}{Q_p},\ \frac{Q_p^{\mathsf L}}{Q_p}\right),\qquad f_{\mathsf D}+f_{\mathsf V}+f_{\mathsf L}=1.}$$

這是 feature 預測能量的**成分分解**，如實描述「這個 feature 的關聯落在哪些座標上」。

### §19.2 為什麼 raw fingerprint 不能直接當作機制判讀

各 block 的維度不同：$K=4$ 時 $(d_{\mathsf D},d_{\mathsf V},d_{\mathsf L})=(4,10,4)$。

考慮一個與未來報酬完全無關的 feature。其 $c_p^{\rm full}$ 在 $d_\Phi=18$ 維空間中大致各向同性，因此期望的能量分配為

$$\mathbb E[f_{\mathsf D}]\approx\frac4{18}=0.222,\qquad \mathbb E[f_{\mathsf V}]\approx\frac{10}{18}=0.556,\qquad \mathbb E[f_{\mathsf L}]\approx\frac4{18}=0.222.$$

也就是說，**一個純雜訊 feature 的 raw fingerprint 天然就是 $(0.22,0.56,0.22)$**。看到 $f_{\mathsf V}=0.55$ 就宣稱「這個 feature 主要預測活躍度」是錯誤的推論；那只是維度的算術。

更精確地說，若投影到 $V_r$ 之後，虛無期望值為

$$\mathbb E_0\big[Q_p^{\mathsf B}\big]\approx\sigma^2\sum_{j\le r_{\rm sig}}\big\|v_j[\mathsf B]\big\|^2,$$

即各 block 在 $V_r$ 中所佔的能量比例——這依 $V_r$ 的實際組成而定，通常**不等於**維度比例。因此解析近似不可靠，必須用 permutation 直接估計。

### §19.3 Null-adjusted enrichment fingerprint

$$\boxed{e_p^{\mathsf B}=\frac{Q_p^{\mathsf B}}{\operatorname{median}_b\,Q^{{\mathsf B},{\rm null}}_{p,b}+\epsilon},\qquad \tilde f_p^{\mathsf B}=\frac{e_p^{\mathsf B}}{e_p^{\mathsf D}+e_p^{\mathsf V}+e_p^{\mathsf L}}.}$$

在虛無下 $\mathbb E[\tilde f_p]\approx(1/3,1/3,1/3)$，因此 $\tilde f$ 有一個與維度無關的參考點。

注意分母使用**該 feature 自己的** null 中位數，而非全 library 的共同基準。這是必要的：覆蓋率、持續性 $\kappa_p$、橫截面分布形狀都會影響單一 feature 的 null 尺度。

**報表規則**：$f_p$ 與 $\tilde f_p$ 並列。$f_p$ 用於回答「能量在哪裡」（例如做部位歸因時），$\tilde f_p$ 用於回答「機制是什麼」（例如做分群與解讀時）。§35 的分群、§53 的判讀一律使用 $\tilde f$。

### §19.4 實例：fingerprint 如何改變決策

兩個 feature 的 $Q$ 相同，$Q=5\times10^{-4}$：

| | $Q^{\mathsf D}$ | $Q^{\mathsf V}$ | $Q^{\mathsf L}$ | $f$ | $\tilde f$ |
|---|---:|---:|---:|---|---|
| A | $1.1\times10^{-4}$ | $2.8\times10^{-4}$ | $1.1\times10^{-4}$ | $(0.22,0.56,0.22)$ | $(0.33,0.34,0.33)$ |
| B | $3.5\times10^{-4}$ | $1.0\times10^{-4}$ | $5\times10^{-5}$ | $(0.70,0.20,0.10)$ | $(0.79,0.13,0.08)$ |

A 的 raw fingerprint 看起來「以活躍度為主」，但 $\tilde f$ 顯示它在三個 block 上都只是雜訊水準——它的 $Q$ 全部來自維度效應，實際上不應通過 §34 的顯著性門檻。B 則是真正的方向性 feature。

若只看 $f$，兩者會被歸入不同的 mechanism group 並各自被保留；看 $\tilde f$ 才會發現 A 根本不該進入這個流程。

---
# 第四部分：目標方向、相似度與穩定度

## §20. Intended target direction

### 動機

$Q^{\mathsf D}_p$ 告訴我們一個 feature 有多少方向性資訊，但沒有告訴我們**那個方向是不是我們要的**。在一個 $K=3$ 的 band 中，$\mathsf D$ block 是三維的；一個 feature 可能在其中某個與我們的交易 horizon 正交的方向上非常強。這在數學上完全合法，在交易上卻毫無用處。

因此必須定義「我們要的方向」。

### 定義

使用者指定 horizon 權重 $a\in\mathbb R^K$（在**原始報酬單位**下，whitening 之前），定義目標訊號

$$\boxed{y^*_d=Y_d^{(\beta)}a\ \in\mathbb R^{N_d}.}$$

例如 $H_2=\{1h,2h,4h\}$ 的增量 band 中，$a=(1,1,1)^\top$ 代表「我關心持有到 4 小時的總報酬」；$a=(1,0,0)^\top$ 代表「我只關心第一個小時」。

$\mathsf D$ 座標中的目標方向定義為 $y^*$ 在 $\Phi_{\mathsf D}$ 上的投影並正規化：

$$\boxed{r_{\mathsf D}=\frac{\Phi_{\mathsf D}^\top Wy^*}{\big\|\Phi_{\mathsf D}^\top Wy^*\big\|},\qquad\|r_{\mathsf D}\|=1.}$$

**這是實作規範。** 它只需要一行程式，且對 §7 的秩截斷、對任何座標慣例都免疫——因為它直接量測「目標訊號本身落在 $\Phi_{\mathsf D}$ 的哪裡」，而不依賴任何閉式反推。

### 命題 20.1（閉式等價）

**無截斷情形（$K_{\rm eff}=K$）**：$\Phi_{\mathsf D}=Z=Y^{(\beta)}\Sigma_Y^{-1/2}$，故

$$\Phi_{\mathsf D}^\top Wy^*=\Sigma_Y^{-1/2}\big(Y^{\top}WY\big)a=\Sigma_Y^{-1/2}\Sigma_Ya=\Sigma_Y^{1/2}a,$$

$$\boxed{r_{\mathsf D}=\frac{\Sigma_Y^{1/2}a}{\sqrt{a^\top\Sigma_Ya}}.}$$

**有截斷情形**：$B_Y=E_{1:K_{\rm eff}}\Lambda_{1:K_{\rm eff}}^{-1/2}$，故

$$\Phi_{\mathsf D}^\top Wy^*=B_Y^\top\Sigma_Ya=\Lambda^{-1/2}E^\top\Sigma_Ya=\Lambda^{1/2}E^\top a\quad(\text{取前 }K_{\rm eff}\text{ 個分量}),$$

$$r_{\mathsf D}=\frac{\Lambda_{1:K_{\rm eff}}^{1/2}E_{1:K_{\rm eff}}^\top a}{\big\|\Lambda_{1:K_{\rm eff}}^{1/2}E_{1:K_{\rm eff}}^\top a\big\|}.$$

∎

閉式僅供理解之用：它顯示 $r_{\mathsf D}$ 的方向被 $\Sigma_Y$ 的譜「重新加權」——變異較大的報酬方向在目標中佔較重的比例。這是正確的：在白化座標中，要達成一單位的 $y^*$，需要在高變異方向上投入更多。

**實作上一律使用投影式。** 這消除了整類最隱蔽的錯誤：任何座標系的旋轉、截斷或重排，投影式都會自動跟上，閉式則不會。

### 實例

$H_2=\{1h,2h,4h\}$ 增量 band，$\Sigma_Y=\operatorname{diag}(\sigma_1^2,\sigma_2^2,\sigma_3^2)$，$\sigma=(0.004,0.005,0.007)$（每期標準差）。取 $a=(1,1,1)$：

$$\Sigma_Y^{1/2}a=(0.004,0.005,0.007),\qquad \|\cdot\|=0.00949,\qquad r_{\mathsf D}=(0.421,0.527,0.738).$$

目標方向在最長 horizon 上載重最重，因為它貢獻最多的報酬變異。若改取 $a=(1,0,0)$ 則 $r_{\mathsf D}=(1,0,0)$。

---

## §21. Signed Target Alignment

### 定義

$$\boxed{A_p^{\rm tgt}=(c_p^{\mathsf D})^\top r_{\mathsf D}.}$$

這是一個**有號**的量。它的語意是：feature $p$ 對目標訊號 $y^*$ 的平均每期橫截面 IC（在 $\Phi_{\mathsf D}$ 的度量下）。

### 解讀

- $A^{\rm tgt}_p>0$：feature 值高 $\Rightarrow$ 目標報酬高。
- $A^{\rm tgt}_p<0$：反向訊號。這**不是缺陷**——只要符號在時間上穩定（§27 的 $\pi^{\rm time,sign}$），反向訊號與正向訊號同樣可用，只需在建構部位時取負號。
- $A^{\rm tgt}_p\approx0$ 但 $Q^{\mathsf D}_p$ 很大：feature 有方向性資訊，但指向與目標正交的方向。應標記 `TARGET_MISMATCH` 並轉入其他 band 或其他 $a$ 的評估。

### 陷阱：符號必須被凍結

$A^{\rm tgt}$ 的符號在研究期一旦確定，就必須被視為 feature 定義的一部分並鎖定。若在 OOS 期間發現符號翻轉而「順手」改回來，那等於用 OOS 資料做了一次決策，整個 OOS 評估即失效。符號翻轉是一個**結果**，不是一個可以修正的 bug（§28）。

### 顯著性

$(A^{\rm tgt}_p)^2$ 以與 $Q$ 完全相同的方式做 permutation 檢定與 Westfall–Young 調整，得到 $\tilde p^{\rm tgt}_p$。這個量在 §34 的 eligibility 中扮演關鍵角色，因為它把「有多少方向性資訊」與「瞄得多準」合併成單一個可直接檢定的量——它恰好就是「對 intended target 的平均 IC 平方」。

---

## §22. Directional Target Efficiency

### 定義

$$\boxed{\eta_p^{\mathsf D}=\frac{(A_p^{\rm tgt})^2}{Q_p^{\mathsf D}+\epsilon}.}$$

### 命題 22.1（界）

$$0\le\eta^{\mathsf D}_p\le1.$$

**證明**：由 Cauchy–Schwarz，$(c_p^{\mathsf D\top}r_{\mathsf D})^2\le\|c_p^{\mathsf D}\|^2\|r_{\mathsf D}\|^2=Q_p^{\mathsf D}$。∎

### 幾何解讀

$\eta^{\mathsf D}_p=\cos^2\theta$，其中 $\theta$ 是 $c_p^{\mathsf D}$ 與 $r_{\mathsf D}$ 的夾角。

- $\eta\approx1$：feature 的方向性資訊**幾乎全部**對準目標。
- $\eta\approx0.3$：只有 30% 的方向性能量在目標方向上，其餘 70% 指向 band 內其他 horizon 組合。
- $\eta\approx0$：完全錯位。

### 三個量的分工

$$Q^{\mathsf D}_p\ \text{（有多少方向性資訊）},\quad A^{\rm tgt}_p\ \text{（指向哪邊、多準）},\quad \eta^{\mathsf D}_p\ \text{（精準度比例）}.$$

三者滿足 $(A^{\rm tgt}_p)^2=Q^{\mathsf D}_p\cdot\eta^{\mathsf D}_p$，因此只有兩個自由度。報表中三個都列，是因為它們回答不同的實務問題：$Q^{\mathsf D}$ 決定「這個 feature 值不值得保留」，$\eta^{\mathsf D}$ 決定「對當前目標值不值得」，$A^{\rm tgt}$ 決定「怎麼用」。

### 實例

$$Q^{\mathsf D}_p=4.0\times10^{-4},\quad A^{\rm tgt}_p=-1.6\times10^{-2},\quad \eta^{\mathsf D}_p=\frac{2.56\times10^{-4}}{4.0\times10^{-4}}=0.64.$$

解讀：這是一個反向的方向性 feature，對目標的平均 IC 為 $-0.016$，其方向性能量有 64% 對準目標。這是一個完全可用的訊號：以 $-X_p$ 進入多空排序。

另一個 feature：

$$Q^{\mathsf D}_q=9.0\times10^{-4},\quad A^{\rm tgt}_q=3.0\times10^{-3},\quad \eta^{\mathsf D}_q=0.01.$$

$Q^{\mathsf D}$ 更大，但 $\eta=0.01$：它的方向性資訊幾乎全部指向與目標正交的方向。若目標 $a=(1,1,1)$，這個 feature 可能是在預測「1h 漲而 4h 跌」這種 horizon 對沖結構。對當前目標無用，但它在 $a=(1,0,-1)$ 的研究中可能非常強。標記 `TARGET_MISMATCH`，保留在 library。

### $\epsilon$ 的作用

分母的 $\epsilon$ 只是數值保護。當 $Q^{\mathsf D}_p$ 本身就在雜訊水準（未通過 §18 的顯著性）時，$\eta$ 是兩個雜訊量的比值，數值上可以是 0 到 1 之間的任何值，**不具意義**。因此 $\eta$ 只應對已通過顯著性檢定的 feature 解讀；報表中對未通過者應顯示為空白而非數字。

---

## §23. Mechanism similarity

### 定義

$$\hat c_p=\frac{c_p}{\|c_p\|+\epsilon},\qquad\boxed{A^{\rm mech}_{pq}=\hat c_p^\top\hat c_q\in[-1,1],\qquad S^{\rm mech}_{pq}=|A^{\rm mech}_{pq}|\in[0,1].}$$

$A^{\rm mech}$ 保留符號，用於判斷兩個 feature 是同向還是反向預測同一個機制。$S^{\rm mech}$ 去掉符號，用於分群——因為一個 feature 與其負號版本攜帶完全相同的資訊，分群時必須視為同一件事。

### §23.1 相似度的能量解讀

這個量有一個精確的決策意義。若以 feature $p$ 作為 feature $q$ 的代表（即在選擇中保留 $p$、丟棄 $q$），則 $q$ 的機制向量被 $p$ 的方向近似，保留的能量比例為

$$\frac{\|\text{proj}_{\hat c_p}c_q\|^2}{\|c_q\|^2}=\big(S^{\rm mech}_{pq}\big)^2,$$

損失的比例為 $1-(S^{\rm mech}_{pq})^2$。

| $S^{\rm mech}$ | 保留能量 $S^2$ | 損失 |
|---:|---:|---:|
| 0.95 | 0.90 | 10% |
| 0.90 | 0.81 | 19% |
| 0.80 | 0.64 | 36% |
| 0.70 | 0.49 | 51% |

**這才是選擇門檻 $\tau$ 的正確依據**：$\tau$ 直接表達「以代表取代成員時，我願意損失多少機制能量」。$\tau=0.90$ 意味著接受最多 19% 的損失。這個判準與維度無關，也與 library 大小無關。

### §23.2 虛無校準：相似度有多少來自機率

上述決策判準回答「該不該合併」，但沒有回答「觀測到的相似結構是不是真的」。後者需要虛無校準。

**校準用的虛無與 §13 不同。** §13 使用**共同**位移（所有 feature 一起位移），以保留 feature 之間的相依結構——這是 $Q_p$ 檢定所需要的。但相似度校準要問的是「兩個**沒有共享機制**的 feature 能有多像」，因此必須使用**逐 feature 獨立位移**：

$$X_{d,p}\longrightarrow X_{\pi_{s_p}(d),p},\qquad s_p\ \text{各自獨立抽取}.$$

如此 $c_p^{\rm null}$ 與 $c_q^{\rm null}$ 獨立，而各自仍保有真實的邊際幾何（覆蓋率、持續性、橫截面分布）。

由此得到 null 的相似度分布 $\{S^{\rm null}_{pq,b}\}$，並定義**超額配對率**

$$\boxed{\zeta(\tau)=\underbrace{\frac{\#\{(p,q):S^{\rm mech}_{pq}>\tau\}}{P(P-1)/2}}_{\text{觀測}}-\underbrace{\mathbb E_b\left[\frac{\#\{(p,q):S^{\rm null}_{pq,b}>\tau\}}{P(P-1)/2}\right]}_{\text{虛無}}.}$$

### §23.3 低維空間的根本限制

以下計算必須被理解，否則分群結果會被嚴重誤讀。

若 $c_p$ 與 $c_q$ 是 $m$ 維空間中兩個獨立均勻的單位向量，則 $|\cos|$ 的尾機率為

| $m$ | $P(|\cos|>0.9)$ | $P(|\cos|>0.95)$ |
|---:|---:|---:|
| 2 | 28.7% | 20.3% |
| 3 | 10.0% | 5.0% |
| 4 | 3.74% | 1.33% |
| 5 | 1.45% | 0.36% |
| 6 | 0.58% | 0.10% |

（$m=3$ 的情形最容易記：由 Archimedes 的球帶定理，三維單位球面上均勻分布的向量，其與固定方向的餘弦在 $[-1,1]$ 上**均勻分布**，故 $P(|\cos|>\tau)=1-\tau$。）

**這意味著什麼。** $\mathsf D$ block 的相似度 $S^{\mathsf D}$（§24）生活在維度

$$m_{\mathsf D}=\dim\big(\text{span}(V_r)\cap\mathsf D\text{ 座標}\big)\le\min(K,r_{\rm sig})$$

的空間中，典型值是 2 或 3。在 $m_{\mathsf D}=3$、$P=400$（約 $8\times10^4$ 個配對）時，虛無下**預期就有約 8000 個配對的 $S>0.9$**。任何以 FWER 或 Bonferroni 方式校準的門檻都會趨近 1，因而毫無用處。

這不是校準方法的失敗，而是一個真實的結構事實：**在三維空間中，四百個方向彼此之間不可能「遠」。**

### §23.4 由此得到的設計結論

> **規則 23.1**：mechanism clustering（§35）是**描述性與組織性**的工具，用於報表分組與人工解讀；**它不是 feature selection 的決定機制**。選擇由 §38–§39 的 coverage 準則決定，該準則直接衡量「被選集合是否代表了每一個合格 feature」，因此不依賴任何配對門檻的顯著性。
>
> **規則 23.2**：分群門檻 $\tau$ 依 §23.1 的能量損失判準設定（預設 $\tau=0.90$，對應 19% 損失），並**必須同時報告** $\zeta(\tau)$ 與有效機制維度
> $$m_{\mathsf B}=\frac{\big(\sum_j\lambda_j^{\mathsf B}\big)^2}{\sum_j(\lambda_j^{\mathsf B})^2},\qquad \lambda^{\mathsf B}_j=\big\|v_j[\mathsf B]\big\|^2,$$
> 即 $V_r$ 在 block $\mathsf B$ 上能量分布的 participation ratio。
>
> 若 $\zeta(\tau)\le0$ 或 $m_{\mathsf D}<4$，報表必須明確標註「mechanism clustering 在此設定下不具鑑別力」，並且**所有以 cluster 為單位的敘述都不得作為結論**。

這是一個誠實而重要的結論：低維機制空間中，「這兩個 feature 是不是同一件事」這個問題，用夾角是答不出來的。能答出來的是「保留哪些 feature 才能代表整體」，而那是 §38 的問題。

---

## §24. Block-specific mechanism similarity

### 定義

$$\boxed{A^{\mathsf B}_{pq}=\frac{c_p^{\mathsf B\top}c_q^{\mathsf B}}{\|c_p^{\mathsf B}\|\|c_q^{\mathsf B}\|+\epsilon},\qquad S^{\mathsf B}_{pq}=\big|A^{\mathsf B}_{pq}\big|,\qquad \mathsf B\in\{\mathsf D,\mathsf V,\mathsf L\}.}$$

### 為什麼這一節不可省略

全空間的 $S^{\rm mech}_{pq}$ 可能很高，卻完全不代表兩者的方向性機制相同。考慮兩個 feature：兩者的 $c$ 都以 $\mathsf V$ block 為主（$f_{\mathsf V}\approx0.8$）且 $\mathsf V$ 成分幾乎一致，但 $\mathsf D$ 成分互相正交。則

$$S^{\rm mech}_{pq}\approx0.8,\qquad S^{\mathsf D}_{pq}\approx0.$$

在做**方向性研究**時，這兩個 feature 完全不可互相取代。

> **規則 24.1**：所有與方向性選擇有關的相似度運算（§35 的分群、§38 的 coverage、§52 的繼承判讀）一律使用 $S^{\mathsf D}$，不使用 $S^{\rm mech}$。全空間的 $S^{\rm mech}$ 只用於「這兩個 feature 整體上像不像」這個描述性問題。

同一規則的推論：若研究目標是波動度訊號，則全部改用 $S^{\mathsf V}$。block 的選擇必須與研究目標一致，這一點在 §27 的穩定度、§34 的門檻上同樣適用。

---

## §25. Raw feature redundancy

### 動機

$S^{\mathsf D}$ 衡量的是「預測機制是否相同」。但兩個 feature 可能預測同一個機制，卻攜帶不同的原始資訊（例如一個來自價格、一個來自成交量，只是恰好都預測短期反轉）。保留兩者可能有價值——它們的雜訊不同，組合起來訊噪比更好。

反過來，兩個原始上幾乎相同的 feature（例如 20 日動能與 21 日動能），保留兩者只是重複。

因此需要第二層、直接在原始資料上度量的冗餘。

### 定義

每日計算 Spearman 相關

$$\rho_{pq,d}=\operatorname{Spearman}\big(X_{d,p},X_{d,q}\big),$$

$$\boxed{S^{\rm raw}_{pq}=\operatorname{median}_d\big|\rho_{pq,d}\big|.}$$

使用 Spearman 而非 Pearson 是因為 raw feature 的分布常有厚尾；使用中位數而非平均是因為個別日期的極端值（例如財報日的資料異常）不應主導冗餘判斷。

### 與 $S^{\mathsf D}$ 的關係

兩者度量不同的東西，且**沒有必然的大小關係**：

| | $S^{\rm raw}$ 高 | $S^{\rm raw}$ 低 |
|---|---|---|
| $S^{\mathsf D}$ 高 | 真正的重複，只留一個 | 不同來源、相同機制：可能都留（雜訊不同） |
| $S^{\mathsf D}$ 低 | 原始相關但預測不同機制：都留 | 完全不同的兩件事 |

右上角的情形是 feature library 建設中最有價值的：兩個獨立的資料來源指向同一個預測機制，這是機制為真的強力佐證，且組合起來可以降低量測雜訊。

---

## §26. Temporal stability：bootstrap 設計

### 動機

$Q_p$ 是全期平均。一個在樣本前半極強、後半消失的 feature，與一個穩定中等的 feature，可能有相同的 $Q_p$，但前者不可用。因此需要量測「這個結論在不同時間子樣本上是否重現」。

### 規格

**Moving block bootstrap（對日期）**。設區塊長度 $\ell_{\rm blk}\ge H_{\max}$，每次重抽 $\lceil T/\ell_{\rm blk}\rceil$ 個起點、串接成長度 $T$ 的日期序列（允許重複），得到

$$\boxed{C^{(b)}=\frac1{T}\sum_{d\in\text{bootstrap}_b}C^{(d)},\qquad b=1,\dots,B_{\rm time}.}$$

預設 $B_{\rm time}=500$–$1000$。

**必須用 moving block 而非 iid bootstrap**：$C^{(d)}$ 在時間上有自相關（feature 持續、報酬波動度叢聚）。iid 重抽會破壞這個結構，低估變異，使所有 $\pi$ 看起來過於樂觀。區塊長度至少為 $H_{\max}$ 確保重疊 label 窗口落在同一區塊內。

**$\Phi$ 與 $V_r$ 在所有重抽中固定。** 這是刻意的：這組 bootstrap 量測的是「**給定座標系，feature 層級關聯的穩定度**」。座標系本身的不確定性由另外兩個機制量測：§7.2 的譜 bootstrap（$\Sigma_Y$ 的穩定度）與 §10 的 $\delta_{\rm drift}$（跨期漂移）。

> **規則 26.1**：$\pi$ 與 $M$ 是總不確定性的**下界**。報表在呈現它們時必須同時呈現 $\delta_{\rm drift}$ 與 $K_{\rm eff}$ 的 bootstrap 穩定度，否則讀者會把「給定座標系的穩定」誤讀為「整體穩定」。

---

## §27. Temporal stability metrics

### 定義

所有指標都是 **block-specific** 的：$\mathsf B$ 必須與研究目標一致（方向性研究一律用 $\mathsf D$）。

**存活率（強度）**：

$$\boxed{\pi^{{\rm time},\mathsf B}_p=\frac1{B_{\rm time}}\sum_b\mathbf 1\left[Q^{\mathsf B,(b)}_p>Q^{\mathsf B,{\rm null}}_{p,0.95}\right],}$$

其中 $Q^{\mathsf B,\rm null}_{p,0.95}$ 是該 feature、該 block 的虛無分布第 95 百分位。

**符號穩定率**：

$$\boxed{\pi^{\rm time,sign}_p=\frac1{B_{\rm time}}\sum_b\mathbf 1\left[\operatorname{sign}\big(A^{{\rm tgt},(b)}_p\big)=\operatorname{sign}\big(A^{\rm tgt}_p\big)\right].}$$

**機制穩定度**：

$$M^{{\rm time},\mathsf B}_{p,b}=\frac{c_p^{\mathsf B,(b)\top}c_p^{\mathsf B}}{\|c_p^{\mathsf B,(b)}\|\,\|c_p^{\mathsf B}\|+\epsilon},\qquad\boxed{M^{{\rm time},\mathsf B}_p=\operatorname{median}_bM^{{\rm time},\mathsf B}_{p,b}.}$$

**強度變異係數**：

$$\boxed{CV^{Q,\mathsf B}_p=\frac{\operatorname{sd}_b\big(Q^{\mathsf B,(b)}_p\big)}{\operatorname{mean}_b\big(Q^{\mathsf B,(b)}_p\big)+\epsilon}.}$$

### 四個量各自回答什麼

| 量 | 問題 | 失敗代表 |
|---|---|---|
| $\pi^{\rm time,\mathsf B}$ | 強度在子樣本中是否持續超過雜訊 | 只在少數時期有效 |
| $\pi^{\rm time,sign}$ | 方向是否穩定 | 方向會翻轉——最危險的失效 |
| $M^{\rm time,\mathsf B}$ | **機制**（$c$ 的方向）是否穩定 | 預測的東西本身在變 |
| $CV^{Q,\mathsf B}$ | 強度的波動程度 | 高 $CV$ 但高 $\pi$：強度忽大忽小但從未消失 |

$M$ 與 $\pi$ 的區別特別重要。$\pi^{\rm time}$ 只看「強度夠不夠大」，$M^{\rm time}$ 看「方向有沒有變」。一個 feature 可以有 $\pi=0.95$（永遠很強）但 $M=0.2$（每個子樣本預測的機制都不同）——這代表它是一個不穩定的組合訊號，$Q$ 很大只是因為它在每個時期都關聯到某個不同的東西。這種 feature 完全不可用，而只看 $\pi$ 會漏掉它。

### 與 p-value 的關係

$\pi^{\rm time,\mathsf B}$ 與 $\tilde p_p$ 看似重複，其實回答不同問題：$\tilde p_p$ 問「全樣本的點估計是否可能來自巧合」，$\pi$ 問「這個關聯在時間子樣本中出現的頻率」。一個真實但只在 30% 的時間有效的 feature，可以有極小的 $\tilde p$ 與只有 0.3 的 $\pi$。兩者都必須列入 eligibility。

### 註記：$\pi^{\rm time,\mathsf B}$ 的門檻比較

$\pi$ 的定義中，bootstrap 樣本的 $Q$ 與**全樣本 null 的**第 95 百分位比較。嚴格地說，bootstrap 子樣本的有效期數較短，其 $Q$ 的變異應略大於全樣本，因此這個比較略偏保守。這是可接受的：$\pi\ge0.8$ 本來就是一個篩選門檻而非精確檢定。若需要更精確的版本，可對每個 bootstrap 重抽同步重抽 null（成本乘以 $B_{\rm time}$），但實務上收益不足以抵銷成本。

---
## §28. 穩定度的兩種失效模式

$\pi$ 與 $M$ 的組合把不穩定分成兩類，它們的處置完全不同。

### 模式一：強度衰減（$\pi^{\rm time,\mathsf D}$ 低，$M^{\rm time,\mathsf D}$ 高）

feature 預測的**機制沒變**，只是強度在某些時期掉到雜訊以下。

典型成因：市場 regime（波動度、流動性）改變了訊號的可偵測性；或訊號被市場逐步套利掉（alpha decay）。

處置：這類 feature 仍然可用，但需要額外的資訊才能決定怎麼用。關鍵診斷是**衰減是否單調**：把 $Q^{\mathsf D,(b)}_p$ 對 bootstrap 樣本的平均日期做迴歸，若斜率顯著為負，這是 alpha decay，應降權或淘汰；若強度與某個可觀測的 regime 變數（VIX、成交量、價差）相關，則這是條件性訊號，應搭配 regime 過濾使用。

### 模式二：機制漂移（$M^{\rm time,\mathsf D}$ 低）

feature 在不同時期預測**不同的東西**。$c_p$ 的方向本身在變。

典型成因：feature 的定義隨資料供給改變（例如某個資料欄位的涵蓋範圍在中途擴大）；或 feature 是多個底層效應的混合，而各效應的相對強度隨時間變化。

處置：這類 feature **不應被使用**，即使 $Q_p$ 很大、$\tilde p_p$ 很小。原因是它在 OOS 期間預測什麼完全無法預期。正確的處理是回到 feature 建構層：把它拆解成更基本的成分，分別檢驗。

$M^{\rm time,\mathsf D}<0.3$ 應被視為硬性淘汰條件，而非一個評分項目。

### 符號翻轉（$\pi^{\rm time,sign}$ 低）

這是模式二的一個特例，但危險性最高，因此單獨標記。一個符號會翻轉的 feature，在部位上會直接把獲利變成虧損。$\pi^{\rm time,sign}<0.8$ 的 feature 一律淘汰，且**不得**以「在 OOS 中改回符號」的方式補救（§21）。

---

## §29. Cross-sectional stability

### 動機

時間穩定度回答「在不同時期是否成立」，橫截面穩定度回答「是否由少數個股驅動」。兩者是獨立的失效維度：一個 feature 可以在每一天都有效（$\pi^{\rm time}=1$），但每一天的效果都來自同樣的 20 檔小型股。

### 規格

保留所有日期。對每個日期 $d$，在 $\mathcal U_d$ 內重新 bootstrap tickers。若有 sector / cap metadata，使用**分層** bootstrap：

$$\mathcal U_d=\bigsqcup_g\mathcal U_{d,g},$$

在每個 stratum 內做 $N_{d,g}$ 次有放回抽樣，最終仍得到 $N_d$ 列。

**分層是必要的**：若不分層，某些 bootstrap 樣本會隨機過度代表某個產業，使 $C^{(d,b)}$ 的變異同時反映「產業組成的抽樣雜訊」與「feature 效力的橫截面集中度」。分層固定了產業組成，使結果單純反映後者。

### 每次重抽必須重新中心化

有放回抽樣之後，$X$ 在重抽樣本上的均值**不再為零**（某些名次被抽中多次）。由引理 4.1，$C^{(d)}$ 之所以是一個乾淨的共變異數，前提正是 $X^\top\mathbf 1=0$。前提被破壞後，

$$\frac1{N_d}X_{d,b}^\top\Phi_{d,b}=\underbrace{\operatorname{Cov}_b(X,\Phi)}_{\text{想要的}}+\underbrace{\bar X_{d,b}\,\bar\Phi_{d,b}^\top}_{\text{偽訊號}}.$$

偽訊號項對**與抽樣機率無關**的 feature 只是額外雜訊，但對與 stratum 大小相關的 feature（市值、流動性、產業啞變數）會產生系統性偏誤。

> **規格 29.1**：每次 bootstrap 重抽之後，必須在重抽樣本上重新去均值與重新標準化 $X$：
> $$\check X_{idp}=\frac{X_{idp}-\bar X^{(b)}_{dp}}{\sigma^{(b)}_{dp}+\epsilon}.$$

同樣的要求適用於 §13.3(b) 的 permutation 交集。兩處是同一個原理的兩個應用。

### 輸出

$$C^{(b)}_{\rm XS}=\frac1T\sum_d\frac1{N_d}\check X_{d,b}^\top\Phi_{d,b},$$

固定 $V_r$，以與 §27 完全相同的方式得到

$$\pi^{{\rm XS},\mathsf B}_p,\qquad \pi^{\rm XS,sign}_p,\qquad M^{{\rm XS},\mathsf B}_p.$$

預設 $B_{\rm XS}=500$。

---

## §30. Predictive breadth

### 動機

$Q_p$ 是一個平均量，它不告訴我們這個平均由多少檔股票貢獻。一個由 3000 檔股票均勻貢獻的 IC，與一個由 20 檔股票貢獻、其餘 2980 檔毫無關係的 IC，在交易上是完全不同的物件：後者的容量極小、對單一事件極度敏感、且在真實交易中極可能被成本吃掉。

### 定義

對 feature $p$，令 $\hat c_p=c_p/\|c_p\|$。第 $d$ 期第 $i$ 檔股票對預測關係的**有號貢獻**為

$$\boxed{g_{id,p}=X^{XS}_{id,p}\big\langle\Phi_{id,:},\hat c_p\big\rangle.}$$

這個量有精確意義：$\frac1{N_d}\sum_ig_{id,p}$ 恰好是 $c_p$ 沿 $\hat c_p$ 方向的當期貢獻，因此

$$\frac1{T}\sum_d\frac1{N_d}\sum_ig_{id,p}=\hat c_p^\top\,c_p^{\rm full}\approx\|c_p\|.$$

絕對貢獻權重與有效貢獻者數：

$$u_{id,p}=\frac{|g_{id,p}|}{\sum_j|g_{jd,p}|},\qquad\boxed{N^{\rm eff}_{d,p}=\frac1{\sum_iu_{id,p}^2},\qquad 1\le N^{\rm eff}_{d,p}\le N_d,}$$

$$\boxed{b_{d,p}=\frac{N^{\rm eff}_{d,p}}{N_d}\in(0,1].}$$

報表使用

$$b^{\rm med}_p=\operatorname{median}_db_{d,p},\qquad b^{\rm p10}_p=\text{第 10 百分位}_d(b_{d,p}).$$

$N^{\rm eff}$ 是 $\{|g|\}$ 分布的 Herfindahl 倒數，即「等效的貢獻檔數」。

---

## §31. Breadth 的解讀與虛無基準

### §31.1 為什麼 $b=1$ 不是參考點

直覺上 $b=1$（所有股票等量貢獻）像是「完全分散」的參考點。但這個直覺是錯的。

考慮一個**完全沒有預測力**的 feature。此時 $g_{id}=X_{id}\cdot\langle\Phi_{id},\hat c\rangle$ 是兩個近似獨立、近似標準常態的量的乘積。在 $N\to\infty$ 下，

$$b=\frac{\big(\frac1N\sum|g|\big)^2}{\frac1N\sum g^2}\longrightarrow\frac{(\mathbb E|g|)^2}{\mathbb E[g^2]}=\frac{(\mathbb E|u|\cdot\mathbb E|v|)^2}{\mathbb E[u^2]\mathbb E[v^2]}=\frac{(2/\pi)^2}{1}\approx0.405.$$

（用到 $\mathbb E|u|=\sqrt{2/\pi}$ 對標準常態。）

因此**純雜訊的 breadth 約為 0.405，不是 1**。若貢獻是單一常態量（而非兩個量的乘積），基準會是 $2/\pi\approx0.637$。看到 $b=0.35$ 就說「集中度很高」是錯的——它幾乎就是雜訊水準。

### §31.2 相對 breadth

真實資料的 $X$ 與 $\Phi$ 都不是獨立常態（$X$ 通常有厚尾，$\Phi_{\mathsf V}$ 的分量是二次型因而嚴重偏斜），所以 0.405 只是一個數量級上的定位，不應作為硬門檻。報表使用的是**相對量**：

$$\boxed{\tilde b_p=\frac{b^{\rm med}_p}{\operatorname{median}_b\big(b^{{\rm med},{\rm null}}_{p,b}\big)},}$$

其中 null 來自 §13 的同一組重抽。$\tilde b_p\approx1$ 代表集中度與雜訊無異；$\tilde b_p\ll1$ 代表真正的集中。

**為什麼相對量才有意義**：$b$ 的虛無基準依 block 而異。$\mathsf V$ block 的 $\Phi$ 分量是 $z$ 的二次型，其峰態遠高於常態，因此 $\mathsf V$ block 的虛無 breadth 顯著低於 0.405。若用同一個絕對門檻去評判 $\mathsf D$ 與 $\mathsf V$ block 的 feature，會系統性地把所有 $\mathsf V$ block feature 標記為「過度集中」。

### §31.3 集中度不是缺陷

> **規則 31.1**：breadth 是 soft diagnostic，不是 hard gate。

集中度本身不必然是缺陷——某些真實的機制（事件驅動、財報意外、指數調整）本來就只影響少數名次。硬性排除它們會系統性地偏向平庸的廣度型訊號。

標記規則：

- $\tilde b_p<0.5$ → `CONCENTRATED`
- $\tilde b^{\rm p10}_p<0.2$ → `EXTREME_CONCENTRATION`

這些標記的用途是**指示後續該做什麼檢查**，而非直接淘汰：

1. 檢查貢獻最大的名次是否有共同特徵（同一產業？同一事件類型？同一資料供應商的異常？）。
2. 檢查在剔除貢獻前 1% 的名次後，$Q_p$ 是否仍然顯著。
3. 估計實際可執行的容量。

### 實例

| feature | $b^{\rm med}$ | null $b$ | $\tilde b$ | 判讀 |
|---|---:|---:|---:|---|
| 動能 | 0.42 | 0.41 | 1.02 | 廣度型，與雜訊基準無異 |
| 財報意外 | 0.06 | 0.40 | 0.15 | 高度集中，標記並檢查容量 |
| 波動度預測 | 0.20 | 0.22 | 0.91 | 看似集中，實為 $\mathsf V$ block 的天然基準 |

第三列說明了相對化的必要性：若以 0.405 為門檻，這個 feature 會被誤標為集中。

---

## §32. Cancellation ratio

### 定義

$$\boxed{R^{\rm canc}_{d,p}=\frac{\big|\sum_ig_{id,p}\big|}{\sum_i|g_{id,p}|+\epsilon}\in[0,1].}$$

### 與 breadth 的分工

兩者都由 $\{g_{id}\}$ 導出，但問的是不同的問題：

- **breadth**：貢獻分散在多少檔股票上？（$|g|$ 的集中度）
- **cancellation**：這些貢獻是否互相抵銷？（$g$ 的符號一致性）

一個 feature 可以 breadth 很高但 cancellation 很低——3000 檔股票都有貢獻，但一半為正一半為負，淨效果來自一個很小的不平衡。這種訊號的特徵是：交易成本會吃掉絕大部分，因為你必須交易 3000 檔股票才能獲得那個小小的淨值。

### 虛無尺度

與 breadth 同理，$R^{\rm canc}=1$ 不是參考點。在虛無下 $g_{id}$ 對稱、均值為零，因此

$$\Big|\sum_ig_{id}\Big|\sim\big|\mathcal N(0,N\sigma_g^2)\big|\ \Longrightarrow\ \mathbb E\Big|\sum_ig\Big|=\sigma_g\sqrt{\frac{2N}\pi},$$

而 $\sum_i|g|\approx N\,\mathbb E|g|$。故

$$\mathbb E\big[R^{\rm canc}\big]\approx\frac{\sigma_g}{\mathbb E|g|}\sqrt{\frac2{\pi N}}=\frac\pi2\sqrt{\frac2{\pi N}}=\sqrt{\frac{\pi}{2N}},$$

（末式用到 $g=uv$ 的乘積情形：$\sigma_g=1$、$\mathbb E|g|=2/\pi$。）

$$N=500\Rightarrow0.056,\qquad N=1500\Rightarrow0.032,\qquad N=3000\Rightarrow0.023.$$

因此原始的 $R^{\rm canc}$ 不但遠小於 1，而且**隨 $N_d$ 變動**，不同日期之間、不同 universe 之間都不可比。

### 相對 cancellation

$$\boxed{\tilde R_p=\frac{\operatorname{median}_dR^{\rm canc}_{d,p}}{\operatorname{median}_b\big(\operatorname{median}_dR^{{\rm canc},{\rm null}}_{d,p,b}\big)}.}$$

$\tilde R_p\approx1$ 代表淨效果與純雜訊無異——**這是一個嚴重警訊**：它意味著 $Q_p$ 可能來自少數極端名次而非系統性的橫截面關係。$\tilde R_p\ge3$ 代表有真實的、方向一致的橫截面效應。

報表同時列出原始的 $R^{\rm canc}$（供容量與成本估算使用，因為絕對值才有交易意義）與 $\tilde R_p$（供統計判讀使用）。

---

## §33. 有效譜秩

### 定義

$$e_j=\frac{s_j^2}{\sum_ks_k^2},\qquad\boxed{r_{\rm eff}=\frac1{\sum_je_j^2}.}$$

這是奇異值能量分布的 participation ratio：若能量集中在單一 mode，$r_{\rm eff}\to1$；若 $r$ 個 mode 能量相同，$r_{\rm eff}=r$。

### 用途

$r_{\rm eff}$ 與 $r_{\rm sig}$ 回答不同問題：

- $r_{\rm sig}$：**有幾個方向在統計上超出雜訊**（檢定結果）。
- $r_{\rm eff}$：**能量實際上分散在幾個方向上**（描述量）。

兩者的比較本身就是資訊：

| 情形 | 意義 |
|---|---|
| $r_{\rm eff}\approx1\ll r_{\rm sig}$ | 整個 library 幾乎只照亮一個機制；其餘顯著 mode 很弱。library 過於同質，應補充新的資料來源。 |
| $r_{\rm eff}\approx r_{\rm sig}$ | 幾個機制強度相當，library 的多樣性良好。 |
| $r_{\rm eff}\gg r_{\rm sig}$ | 能量分散在許多微弱方向上。通常代表訊號整體很弱，$V_r$ 的估計不穩定，後續結論應降級。 |

$r_{\rm eff}$ 也提供了 feature budget 的一個下界參考：若 $r_{\rm eff}=3$，那麼選出 40 個 feature 在機制層面上必然高度重複——重複本身不是錯（雜訊分散有價值），但應該是一個有意識的決定，而非誤以為選到了 40 種不同的訊號。

---
# 第五部分：Feature Selection

這一部分把前面所有統計量轉成一個決策：**保留哪些 feature。** 設計上分成三個彼此獨立的階段，各自回答一個不同的問題：

1. **Eligibility（§34）**：這個 feature 是不是真的？——二元判定，由顯著性與穩定度決定。
2. **Grouping（§35–§36）**：哪些 feature 講的是同一件事？——描述性組織，供人工解讀。
3. **Coverage selection（§37–§39）**：保留哪個子集才能代表整體？——最佳化問題，由次模覆蓋決定。

三個階段不得混用。特別是：**分群結果不決定選擇**（§23.4 的規則 23.1），選擇由階段三的覆蓋準則決定。

## §34. Eligibility

### 規格（方向性研究）

feature 進入 selection 的條件為以下五項全部成立：

$$\boxed{
\begin{aligned}
&\text{(E1) } \tilde p^{\,\mathsf D}_p<0.05 &&\text{（Westfall–Young 調整後，§18）}\\
&\text{(E2) } \tilde p^{\,\rm tgt}_p<0.05 &&\text{（對 }(A^{\rm tgt}_p)^2\text{ 的檢定，§21）}\\
&\text{(E3) } \pi^{\rm time,\mathsf D}_p\ge0.8\ \text{且}\ \pi^{\rm time,sign}_p\ge0.8 &&\text{（§27）}\\
&\text{(E4) } \pi^{\rm XS,\mathsf D}_p\ge0.8 &&\text{（§29）}\\
&\text{(E5) } M^{\rm time,\mathsf D}_p\ge0.6\ \text{且}\ M^{\rm XS,\mathsf D}_p\ge0.6 &&\text{（§27、§29）}
\end{aligned}}$$

breadth 與 cancellation 為 soft diagnostic：標記但不排除（§31.3）。

### 每一道門檻排除什麼

| 門檻 | 排除的失效模式 |
|---|---|
| E1 | 純屬巧合的關聯 |
| E2 | 有方向性資訊，但指向與目標正交的方向 |
| E3 | 只在少數時期有效；方向會翻轉 |
| E4 | 只由少數個股驅動 |
| E5 | 機制本身在漂移（§28 模式二） |

### E2 為什麼不可省略

E1 只要求「有方向性預測力」。一個在 $\mathsf D$ block 上顯著、穩定、廣泛，但 $\eta^{\mathsf D}_p\approx0$ 的 feature，會完整通過 E1、E3、E4、E5，甚至可能成為 family representative。它預測的是 band 中另一個 horizon 組合——在數學上完全合法，在交易上對當前目標無用。

E2 直接對 $(A^{\rm tgt}_p)^2=Q^{\mathsf D}_p\cdot\eta^{\mathsf D}_p$ 做檢定，同時涵蓋「有多少方向性資訊」與「瞄得多準」。它比「$\eta^{\mathsf D}_p\ge\eta_{\min}$」這種固定門檻更好，因為後者對弱 feature 不穩定（§22 末段：$\eta$ 在 $Q^{\mathsf D}$ 接近雜訊時沒有意義），而 E2 的檢定量本身就有虛無分布。

**被 E2 排除的 feature 不應丟棄**，而應標記 `TARGET_MISMATCH` 並轉入其他 band、其他 $a$ 的評估。它們是 library 的資產，只是與當前目標不匹配。

### 門檻之間不獨立

E1 到 E5 高度相關：一個 $\tilde p$ 極小的 feature，其 $\pi$ 幾乎必然高。因此這五道門檻的**聯合**篩選力遠低於五個獨立門檻的乘積所暗示的。這不是缺陷，而是必須在解讀「通過率」時意識到的事實：若 400 個 feature 中有 60 個通過全部五項，不能推論「這 60 個各自通過了五重獨立檢驗」。真正的多重檢定控制只在 E1 與 E2（兩者都經過 Westfall–Young 調整）。

### 其他研究目標

若研究目標為波動度訊號，把全部 $\mathsf D$ 換成 $\mathsf V$，並刪除 E2（$\mathsf V$ block 沒有 target 方向的概念，因為它不是一個有號的方向性空間）。若目標為大行情方向性，換成 $\mathsf L$ 並保留 E2（$\mathsf L$ 是有號方向空間，$r_{\mathsf D}$ 的類比為 $r_{\mathsf L}=\text{normalize}(\Phi_{\mathsf L}^\top Wy^*)$）。

---

## §35. Mechanism grouping

### 規格

對通過 eligibility 的集合 $\mathcal E$，以 $1-S^{\mathsf D}_{pq}$ 為距離做 **average linkage** 階層分群，切點

$$\tau_{\mathsf D}=0.90\qquad\text{（即 }S^{\mathsf D}\ge0.90\text{ 視為同群）},$$

得到 mechanism groups $G_1,\dots,G_M$。

$\tau_{\mathsf D}=0.90$ 的依據是 §23.1 的能量損失判準：以代表取代成員時，最多損失 $1-0.9^2=19\%$ 的機制能量。

**用 average linkage 而非 complete 或 single**：single linkage 會產生鏈狀群（$p\sim q\sim r$ 但 $p\not\sim r$），在低維空間中尤其嚴重；complete linkage 對單一離群配對過度敏感。average linkage 在兩者之間，且對本問題的距離結構最穩健。

### 必須同時報告的校準量

依規則 23.2，分群結果**必須**與以下兩個量一併呈現：

$$\zeta(\tau_{\mathsf D})\ \text{（超額配對率，§23.2）},\qquad m_{\mathsf D}\ \text{（有效機制維度，§23.4）}.$$

若 $\zeta\le0$ 或 $m_{\mathsf D}<4$，報表必須標註：

> **mechanism grouping 在此設定下不具鑑別力；群的邊界是任意的，不得作為結論引用。**

這種情況在 $r_{\rm sig}\le3$ 時是常態而非例外。此時分群仍可用於報表的視覺組織，但所有實質決策交由 §38–§39。

### 群結構本身是 library 的診斷

- 群數 $M$ 遠小於 $|\mathcal E|$：library 嚴重過度取樣少數機制。這是 feature 建設的警訊——應補充新的資料來源，而非繼續增加同一類 feature 的變體。
- 存在只含一個 feature 的群：該 feature 攜帶獨特機制，在 §39 的覆蓋選擇中會自動獲得高優先序。
- 最大群佔 $|\mathcal E|$ 的一半以上：通常代表 library 由單一資料供應商或單一類型（例如全部是價格衍生量）主導。

---

## §36. Raw redundancy families

### 規格

在每個 mechanism group $G_m$ **內部**，以 $S^{\rm raw}_{pq}$（§25）做 **complete linkage** 分群，切點 $S^{\rm raw}=0.90$，產生 families $R_{m1},R_{m2},\dots$。

因此一個 family 的意義是：**預測機制接近，而且原始橫截面資訊也高度重複。**

### 兩層結構的邏輯

- **只在同一 mechanism group 內檢查 raw 冗餘。** 兩個 raw 高度相關但預測不同機制的 feature（例如原始值預測方向、絕對值預測活躍度），應該都保留。
- **不跨 group 合併。** 即使 $S^{\rm raw}=0.95$，只要機制不同就視為不同 family。

這個順序是正確的，因為最終目標是保留**預測機制的多樣性**；raw 冗餘只是在機制相同時用來挑代表的次級準則。

這裡使用 complete linkage（而 §35 用 average）是刻意的：raw 冗餘的判定應該保守，complete linkage 要求群內**所有**配對都超過門檻，避免把只是間接相關的 feature 併掉。

### Transform family 規則

依規則 4.1，同一 source feature 的所有變換條目 $(p,1),(p,2),\dots$ **強制歸入同一 raw family**，除非它們落在不同的 mechanism group。

實作上這是一個後處理步驟：完成 §36 的分群後，掃描每個 source feature，若其多個變換落在同一 group 的不同 family，則強制合併這些 family。

### 報表

記錄每個 feature 的 `mechanism_group`、`raw_family`，以及兩層的群大小。若某個 mechanism group 內含 30 個 feature 但只分出 2 個 raw family，代表這個機制被嚴重過度取樣。

---

## §37. Family representative

### 規格

每個 raw family $R_{mk}$ 選出一個代表：

$$\boxed{\text{Score}_p=\operatorname{median}_b\left[\big(A^{{\rm tgt},(b)}_p\big)^2\right],\qquad p^*=\arg\max_{p\in R_{mk}}\text{Score}_p,}$$

其中 $b$ 遍歷 §26 的 temporal bootstrap 重抽。

若研究目標未指定 $a$（探索性研究），改用 $\text{Score}_p=\operatorname{median}_b\big[Q^{\mathsf D,(b)}_p\big]$。

### 為什麼是 bootstrap 中位數，而不是點估計

直接用 $(A^{\rm tgt}_p)^2$ 的點估計會有 **winner's curse**：在一個 family 內，點估計最大的 feature 往往是「真實強度中上，且運氣特別好」的那一個。bootstrap 中位數對這個效應有內建的收縮作用——一個靠少數時期撐起來的 feature，其 bootstrap 分布會嚴重右偏，中位數遠低於點估計。

### 為什麼不用多因子乘積分數

一個直覺上自然的做法是把所有品質指標相乘：

$$Q^{\mathsf D}_p\cdot\pi^{\rm time}_p\cdot\pi^{\rm sign}_p\cdot M^{\rm time}_p\cdot\pi^{\rm XS}_p\cdot M^{\rm XS}_p\cdot h(b_p).$$

這個做法有兩個問題。

**第一，它重複計算 eligibility。** 能進入 §37 的 feature 已經滿足 $\pi\ge0.8$（三項）、$M\ge0.6$（兩項）。因此這些因子的實際變動範圍分別只有 $[0.8,1]$ 與 $[0.6,1]$，六個因子合起來的變動倍數最多約 $(1/0.8)^3(1/0.6)^2\approx5.4$ 倍。

**第二，$Q^{\mathsf D}$ 的變動範圍是 100 倍以上**（從 $10^{-5}$ 到 $10^{-3}$）。因此乘積分數在數值上幾乎完全由 $Q^{\mathsf D}$ 決定，其餘六個因子只是雜訊般的微調。既然如此，把它們寫進分數只是製造「考慮周全」的假象，實際上增加了不可追溯的自由度。

正確的設計是：**把品質門檻放在 eligibility（硬門檻，可追溯），把排序放在單一個有明確語意的量（Score）。**

### Breadth 的角色

$\tilde b_p$ 不進入 Score，而是作為**同分裁決**與**事後檢查**：若 family 內前兩名的 Score 差距在 20% 以內，選 $\tilde b$ 較高者。理由見 §31.3：集中度不是缺陷，但在其他條件相近時，廣度型訊號的容量與穩健性較好。

---

## §38. Coverage

本節定義兩個不同的覆蓋概念。**它們回答不同的問題，不可互換。**

### §38.1 Span coverage（診斷用）

$$P^{\mathsf B}_{\mathcal F}=C^{\mathsf B\top}_{\mathcal F}\big(C^{\mathsf B}_{\mathcal F}C^{\mathsf B\top}_{\mathcal F}\big)^+C^{\mathsf B}_{\mathcal F},\qquad\boxed{{\rm Cov}^{\rm span}_{\mathsf B}(\mathcal F)=\frac{\big\|C^{\mathsf B}_{\mathcal E}P^{\mathsf B}_{\mathcal F}\big\|_F^2}{\big\|C^{\mathsf B}_{\mathcal E}\big\|_F^2}\in[0,1],}$$

其中 $C^{\mathsf B}_{\mathcal F}$ 是被選集合 $\mathcal F$ 的機制向量（block $\mathsf B$ 部分）堆成的矩陣，$C^{\mathsf B}_{\mathcal E}$ 是全體合格 feature 的。

**性質 38.1（單調性）**：$\mathcal F\subseteq\mathcal F'\Rightarrow{\rm Cov}^{\rm span}_{\mathsf B}(\mathcal F)\le{\rm Cov}^{\rm span}_{\mathsf B}(\mathcal F')$。

**證明**：$\operatorname{span}(C_{\mathcal F})\subseteq\operatorname{span}(C_{\mathcal F'})$，而投影到較大子空間的範數不減。∎

**性質 38.2（飽和）**：$c_p^{\mathsf D}\in\mathbb R^{d_{\mathsf D}}=\mathbb R^K$。因此只要 $\mathcal F$ 含有 $K$ 個線性獨立的 $c^{\mathsf D}$ 向量，$P^{\mathsf D}_{\mathcal F}=I$，於是

$$\boxed{{\rm Cov}^{\rm span}_{\mathsf D}(\mathcal F)=1.}$$

$K=3$ 時，**任意三個線性獨立的 feature 就達到 100% 覆蓋**。

性質 38.2 使 span coverage 完全不適合作為選擇準則：目標 90% 在選到第三、四個 feature 時就已達成，貪婪步驟隨即停止，而選出的三個 feature 可能是任意三個。更糟的是，這個準則結構上**不可能**選出超過 $d_{\mathsf D}$ 個 feature，即使 library 中有 50 個機制各異、都值得保留的訊號。

> **規則 38.1**：${\rm Cov}^{\rm span}$ 只作為診斷使用，回答「被選集合有沒有張滿所有機制方向」。它**不進入** §39 的選擇迴圈。

它作為診斷仍有價值：若 ${\rm Cov}^{\rm span}_{\mathsf D}(\mathcal F)<1$，代表被選集合連 $\mathsf D$ 空間都沒張滿——這通常意味著所有被選 feature 都指向同一個 horizon 組合，是一個需要處理的問題。

### §38.2 Feature coverage（選擇用）

選擇要回答的真正問題是：**每一個值得保留的 feature，是否都有一個被選中的代表「夠接近」它？**

$$\boxed{{\rm Cov}^{\rm feat}_{\mathsf B}(\mathcal F)=\frac{\sum_{p\in\mathcal E}w_p\max_{j\in\mathcal F}S^{\mathsf B}_{pj}}{\sum_{p\in\mathcal E}w_p},\qquad w_p=\text{Score}_p\ \text{(§37)}.}$$

（約定 $\max_{j\in\varnothing}=0$，且 $p\in\mathcal F$ 時 $S^{\mathsf B}_{pp}=1$。）

語意：以 Score 加權的「平均代表度」。$\rm Cov^{feat}=0.9$ 表示，按重要性加權後，合格 feature 的機制平均有 0.9 的餘弦被某個被選代表捕捉到——對應 $0.9^2=81\%$ 的機制能量。

### §38.3 定理 38.3（次模性）

${\rm Cov}^{\rm feat}_{\mathsf B}$ 是單調非遞減且**次模（submodular）**的集合函數。

**證明**：對固定的 $p$，令 $h_p(\mathcal F)=\max_{j\in\mathcal F}S^{\mathsf B}_{pj}$。

*單調性*：$\mathcal F\subseteq\mathcal F'$ 時取最大值的範圍變大，故 $h_p(\mathcal F)\le h_p(\mathcal F')$。

*次模性*：設 $\mathcal A\subseteq\mathcal B$，$j\notin\mathcal B$。則

$$h_p(\mathcal A\cup\{j\})-h_p(\mathcal A)=\max\Big(0,\ S^{\mathsf B}_{pj}-\max_{k\in\mathcal A}S^{\mathsf B}_{pk}\Big)\ \ge\ \max\Big(0,\ S^{\mathsf B}_{pj}-\max_{k\in\mathcal B}S^{\mathsf B}_{pk}\Big)=h_p(\mathcal B\cup\{j\})-h_p(\mathcal B),$$

不等式成立是因為 $\max_{k\in\mathcal A}\le\max_{k\in\mathcal B}$。

最後，${\rm Cov}^{\rm feat}=\frac1{\sum w}\sum_pw_ph_p$ 是非負加權和，而單調次模函數在非負線性組合下封閉。∎

**推論 38.4（貪婪的保證）**：在基數限制 $|\mathcal F|\le n$ 下，貪婪演算法達到最佳值的 $1-1/e\approx63\%$ 以上（Nemhauser–Wolsey–Fisher）。在「達到覆蓋水準 $\theta$」的形式下，貪婪所需的元素數不超過最佳解的 $O(\log\frac1{1-\theta})$ 倍（Wolsey 的次模覆蓋定理）。

這個保證是 ${\rm Cov}^{\rm feat}$ 相對於 ${\rm Cov}^{\rm span}$ 的實質優勢：它讓「貪婪地一個一個加」這個直覺做法有了正式的品質保證，而不只是啟發式。

### §38.4 兩種 coverage 的對照

| | ${\rm Cov}^{\rm span}$ | ${\rm Cov}^{\rm feat}$ |
|---|---|---|
| 問題 | 有沒有張滿機制空間？ | 每個 feature 有沒有代表？ |
| 上限由什麼決定 | $d_{\mathsf B}$（$\mathsf D$ 為 $K$） | 無結構上限 |
| 達到 1 需要幾個 feature | $K$ 個 | 全部 |
| 次模 | 否（是 matroid rank 型，僅單調） | 是 |
| 用途 | 診斷 | 選擇 |

---

## §39. Direction-first selection

### 規格

**第一步：種子。** $\mathcal F\leftarrow\{$每個 raw family 的代表 $p^*\}$（§37）。

**第二步：方向性覆蓋。** 從剩餘的合格 feature 中貪婪加入

$$p^\dagger=\arg\max_{p\notin\mathcal F}\Big[{\rm Cov}^{\rm feat}_{\mathsf D}(\mathcal F\cup\{p\})-{\rm Cov}^{\rm feat}_{\mathsf D}(\mathcal F)\Big],$$

直到

$$\boxed{{\rm Cov}^{\rm feat}_{\mathsf D}(\mathcal F)\ge0.90}$$

或達到 feature budget $n_{\max}$。

**第三步：補 $\mathsf V$ 與 $\mathsf L$。** 若仍有 budget，對 $\mathsf V$ 與 $\mathsf L$ block 重複第二步（各自使用 $S^{\mathsf V}$、$S^{\mathsf L}$ 與對應的 Score），目標 $0.80$。

**第四步：診斷。** 計算並記錄 ${\rm Cov}^{\rm span}_{\mathsf D}(\mathcal F)$、${\rm Cov}^{\rm span}_{\mathsf V}(\mathcal F)$、${\rm Cov}^{\rm span}_{\mathsf L}(\mathcal F)$。若 ${\rm Cov}^{\rm span}_{\mathsf D}<1$，檢視未被張到的方向是什麼。

### Direction-first 的理由

「先滿足方向性、再補其他 block」這個排序原則的依據是：

1. 下游目標是 cross-sectional 相對強弱排序，直接依賴 $\mathsf D$。
2. $\mathsf V$ 與 $\mathsf L$ 的資訊在組合層級有價值（部位縮放、風險預算、條件性加碼），但它們**無法替代方向**。一個只有 $\mathsf V$ 資訊的組合無法建構多空部位。
3. 若把三個 block 放在同一個覆蓋目標中競爭，由於 $d_{\mathsf V}>d_{\mathsf D}$（§6 的維度表），$\mathsf V$ 會機械性地佔據較多名額。

### Budget 的設定

$n_{\max}$ 應由下游模型的樣本量決定，而非由本方法決定。一個實用的參考：若模型是樹模型，$n_{\max}\approx\sqrt{T^{\mathsf D}_{\rm eff}}$ 量級通常不會嚴重過擬合；若是線性或淺層模型可以更多。

$n_{\max}$ 與 $r_{\rm eff}$（§33）的比較是重要的檢查：若 $r_{\rm eff}=3$ 而 $n_{\max}=40$，那 40 個 feature 在機制層面上必然高度重複。這可能是刻意的（用多個量測降低雜訊），但必須是有意識的決定。

### 終止條件無法達成時

若貪婪跑完全部合格 feature 仍未達到 0.90，代表合格集合內部的機制多樣性很高——每個 feature 都有一部分不被其他 feature 代表。此時應記錄最終達到的覆蓋率並接受它，同時檢視是否 $\tau_{\mathsf D}$ 設得太高（§35）。

---

## §40. Feature report 格式

每個 feature（在 transform family 下，每個條目 $(p,k)$）一列：

**識別與資料品質**

| 欄位 | 說明 |
|---|---|
| `feature_id`, `source_id`, `transform` | 條目、來源 feature、變換 |
| `coverage`, `active_ratio` | $\operatorname{median}_dm_{dp}$、$T_p/T$ |
| `norm_before_renorm` | §4.1 第三步前的 $\|X_p\|_W^2$ |
| `kappa` | 持續性 $\kappa_p$（§13.4） |
| `null_used` | `N1_shift` 或 `N2_stratified` |

**強度與機制**

| 欄位 | 說明 |
|---|---|
| `Q`, `Q_D`, `Q_V`, `Q_L` | §16 |
| `Q_full` | 未投影到 $V_r$ 的 $\|c_p^{\rm full}\|^2$（§16 末段） |
| `f_D`, `f_V`, `f_L` | raw fingerprint（§19.1） |
| `ft_D`, `ft_V`, `ft_L` | enrichment fingerprint $\tilde f$（§19.3） |
| `A_target`, `eta_D` | §21、§22 |

**顯著性**

| 欄位 | 說明 |
|---|---|
| `p_D_raw`, `p_target_raw` | 原始 permutation p-value |
| `p_D_wy`, `p_target_wy` | Westfall–Young 調整後（§18.3） |
| `fdr_selected` | permutation FDR 版本是否選中（§18.4） |
| `p_full` | $Q^{\rm full}$ 的 p-value |

**穩定度**

| 欄位 | 說明 |
|---|---|
| `pi_time_D`, `pi_time_sign`, `M_time_D`, `CV_Q_D` | §27 |
| `pi_XS_D`, `M_XS_D` | §29 |

**橫截面支撐**

| 欄位 | 說明 |
|---|---|
| `b_med`, `b_p10` | §30 |
| `b_tilde` | 相對 breadth（§31.2） |
| `R_cancel_med`, `R_tilde` | §32 |
| `flags` | `CONCENTRATED` / `EXTREME_CONCENTRATION` / `TARGET_MISMATCH` |

**結構與選擇**

| 欄位 | 說明 |
|---|---|
| `mechanism_group`, `raw_family` | §35、§36 |
| `is_representative` | §37 |
| `score` | §37 的 Score |
| `selected` | §39 的最終結果 |
| `marginal_cov_gain` | 被選入時的 $\Delta{\rm Cov}^{\rm feat}_{\mathsf D}$ |

**全域欄位（每次執行一列）**

$T$、$T^{\mathsf D}_{\rm eff}$、$T^{\mathsf V}_{\rm eff}$、$K$、$K_{\rm eff}$、$\kappa(\Sigma_Y)$、$\{\lambda_j/\lambda_1\}$、band 定義（累積或增量）、$r_{\rm sig}$、$r_{\rm eff}$、$m_{\mathsf D}$、$\zeta(\tau_{\mathsf D})$、$\tau_{\mathsf D}$、$B_{\rm null}$、$B_{\rm time}$、$B_{\rm XS}$、推論規格（切分法或重跑法）、$\mathcal D$ 與 $\mathcal E$ 的日期範圍、$\delta_{\rm drift}(\mathcal D\to\mathcal E)$、thinning 比率、${\rm Cov}^{\rm feat}_{\mathsf D}$、${\rm Cov}^{\rm span}_{\mathsf D}$。

---
# 第六部分：Label–Alpha 幾何

## §41. 模型訓練與洩漏控制

feature selection 結束後，模型只使用 training 期允許的資訊訓練，產生

$$\boxed{\hat\alpha_{id}\quad\text{（必須是 strict OOS prediction）}.}$$

若使用 walk-forward，每個 test period 的預測只能來自該期之前的模型。

### 六個必須被納入 OOS 邊界的物件

「strict OOS」不只是模型參數的問題。本方法中至少有六個物件是在研究期擬合的，每一個都可能洩漏：

1. $B_Y$（label whitening）與 $K_{\rm eff}$；
2. $\Theta_{\mathsf V},\Theta_{\mathsf L},\Gamma_{\mathsf{DV}},\Gamma_{\mathsf{DL}},\Gamma_{\mathsf{VL}}$（正交化係數）；
3. $V_r$ 與 $r_{\rm sig}$；
4. feature eligibility 清單；
5. mechanism group 與 raw family 結構；
6. selected feature set。

> **規則 41.1**：即使模型本身完全 OOS，若 selected feature set 是在包含 test period 的資料上選出的，整個評估仍然無效。

這是迭代研究中最常見的洩漏途徑：研究者先用全樣本跑一次 feature selection，然後才開始做 walk-forward 建模。§15 的 discovery / evaluation 切分處理了統計推論層面的挑選效應，但**不自動**處理這一層——若最終模型的 test period 與 $\mathcal E$ 重疊，$\mathcal E$ 上的選擇決策同樣構成洩漏。

正確的做法有二，依成本選擇：

- **完整版**：把 feature selection 本身納入 walk-forward，每個 fold 各自重跑 §13–§39。成本高，但唯一嚴格正確。
- **實務版**：保留一段從未參與任何 selection 或模型決策的 final holdout，並規定該段**只能使用一次**。使用次數必須記錄在報表中。

### Purging 與 embargo

由於 label 有 $H_{\max}$ 的前視窗、feature 有 $\ell_X$ 的記憶，訓練集與測試集之間必須留出

$$\boxed{\Delta=H_{\max}+\ell_X}$$

的間隔：移除 label window 與 test period 重疊的訓練樣本（purge），並在 test period 之前額外留出 $\ell_X$ 的 embargo。

**這與 §13.2 的最小位移是同一個數字。** 實作上必須從同一個設定常數讀取，否則同一份資料會在兩處被以不同的洩漏標準對待，而這種不一致在報表上完全看不出來。

---

## §42. Alpha 前處理

### 定義

$$\boxed{\tilde\alpha_{id}=\frac{\hat\alpha_{id}-\bar\alpha_d}{s_{\alpha,d}+\epsilon},}\qquad \bar\alpha_d=\frac1{N_d}\sum_i\hat\alpha_{id},\quad s_{\alpha,d}^2=\frac1{N_d}\sum_i(\hat\alpha_{id}-\bar\alpha_d)^2.$$

隨後與 feature 同樣做一次全域正規化，使 $\|\tilde\alpha\|_W=1$。

### 這一步丟掉了什麼

within-date 標準化移除了兩個維度的資訊：

- **模型預測的絕對量級**（例如「模型認為今天整體機會很大」）。這是一個 market timing 成分，不在本方法的分析範圍內（§5）。
- **跨期的 conviction 變化**（例如「模型今天的橫截面離散度是平時的三倍」）。

第二點在實務上有價值，應單獨記錄為時間序列 $s_{\alpha,d}$ 並與已實現的報酬離散度比較：若兩者相關，模型具有條件性的 conviction 校準能力，這是一個額外的、本方法不量測的優點。

標準化是必要的，因為它使 $Q_\alpha$ 與 $Q_p$ 可以直接比較——兩者都是「單位變異訊號的預測能量」。

---

## §43. 表徵約定：full / inherited / novel

**這一節是整個第六、七部分的基礎。** 它定義了所有訊號（feature、alpha、殘差）在 $\Phi$ 中的表徵方式。

### §43.1 三個向量

對任何 within-date 標準化、且 $\|s\|_W=1$ 的訊號 $s$：

$$\boxed{c^{\rm full}_s=s^\top W\Phi\ \in\mathbb R^{d_\Phi}}$$

這是 $s$ 的**完整表徵**，不做任何投影。

令 $P_V=V_rV_r^\top$ 為 feature library 的顯著機制子空間投影（§15），定義

$$\boxed{c^{\rm inh}_s=c^{\rm full}_sP_V,\qquad c^{\rm nov}_s=c^{\rm full}_s(I-P_V).}$$

對應的能量

$$Q^{\rm full}_s=\|c^{\rm full}_s\|^2,\qquad Q^{\rm inh}_s=\|c^{\rm inh}_s\|^2,\qquad Q^{\rm nov}_s=\|c^{\rm nov}_s\|^2.$$

**命題 43.1（正交分解）**：

$$Q^{\rm full}_s=Q^{\rm inh}_s+Q^{\rm nov}_s.$$

**證明**：$P_V$ 與 $I-P_V$ 是互補的正交投影，故 $c^{\rm inh}\perp c^{\rm nov}$，由 Pythagoras 定理得證。∎

### §43.2 為什麼 alpha 不能只看投影版本

$V_r$ 的內容是「**原始 feature library 用線性方式能看到的 label mechanism**」。若一個非線性模型學到

$$\hat\alpha=f(X)$$

而這個非線性組合產生一個新的預測機制

$$v_{\rm new}\perp\operatorname{span}(V_r),$$

那麼計算 $c_\alpha^{\rm full}P_V$ 就**直接把它殺掉了**。

接著若再問「模型有沒有學到新的預測結構」，這個問題在已經把新結構投影掉的表徵上是無法回答的。這是一個邏輯上的自相矛盾。

對殘差 alpha（§58）情況更嚴重：$\varepsilon_\alpha=\tilde\alpha-\text{Linear}(X)$ 的**存在目的**就是尋找線性 feature span 之外的東西，把它投影回 feature 導出的子空間等於否定了它的定義。

### §43.3 Novelty fraction

$$\boxed{\nu_s=\frac{Q^{\rm nov}_s}{Q^{\rm full}_s}\in[0,1].}$$

$\nu_\alpha\approx0$：模型幾乎只是在重新組合 feature library 已知的預測機制。
$\nu_\alpha=0.4$：模型預測幾何中有 40% 落在線性 feature 分析原先看不見的方向。

### §43.4 Novelty 必須做虛無校準

與 fingerprint（§19.2）完全相同的維度效應在這裡再次出現，且更嚴重。

在虛無下（$\alpha$ 與 $\Phi$ 無關），$c^{\rm full}_\alpha$ 在 $d_\Phi$ 維空間中近似各向同性，因此

$$\mathbb E_0[\nu_\alpha]\approx\frac{d_\Phi-r_{\rm sig}}{d_\Phi}.$$

$d_\Phi=18$、$r_{\rm sig}=3$ 時，$\mathbb E_0[\nu_\alpha]\approx0.83$。也就是說，**一個完全沒有預測力的隨機訊號，其 novelty fraction 就有 0.83**。把 $\nu_\alpha=0.4$ 讀成「模型有 40% 的新東西」是嚴重的誤讀——0.4 其實遠低於虛無基準，代表模型的幾何**比隨機更集中在已知子空間內**。

因此定義 enrichment 版本：

$$e^{\rm inh}_s=\frac{Q^{\rm inh}_s}{\operatorname{median}_bQ^{{\rm inh},{\rm null}}_{s,b}},\qquad e^{\rm nov}_s=\frac{Q^{\rm nov}_s}{\operatorname{median}_bQ^{{\rm nov},{\rm null}}_{s,b}},$$

$$\boxed{\tilde\nu_s=\frac{e^{\rm nov}_s}{e^{\rm inh}_s+e^{\rm nov}_s}.}$$

虛無下 $\mathbb E[\tilde\nu_s]\approx0.5$。

**報表一律使用 $\tilde\nu$ 作為判讀依據**，$\nu$ 僅作為能量分解的描述量。

| $\tilde\nu_\alpha$ | 判讀 |
|---:|---|
| $<0.3$ | 模型的預測幾何集中在 feature library 已知的機制上 |
| $\approx0.5$ | 新舊成分的 enrichment 相當 |
| $>0.7$ | 模型的預測力主要來自已知機制**之外**的方向 |

高 $\tilde\nu$ 必須搭配 $Q^{\rm nov}_\alpha$ 的顯著性檢定（§46）一起讀：$\tilde\nu$ 高但 $Q^{\rm nov}$ 不顯著，代表新成分只是雜訊。

### §43.5 表徵約定總表

> **約定 43.2**：
> - **Feature 層級的選擇統計量**（§16–§39）使用**投影後**的 $c_p=c_p^{\rm full}P_V$。理由：選擇的目的是在統計上已驗證的機制方向中比較 feature，$P_V$ 是一道必要的降噪。同時報表列出 $Q_p^{\rm full}$ 供孤立機制的偵測。
> - **Alpha 與殘差**（§43–§60）一律以 $c^{\rm full}$ 為主要結果，並同時報告 inherited / novel 分解與 $\tilde\nu$。
> - **Feature 與 alpha 的交叉比較**（§52）使用 $c_p$ 與 $c^{\rm inh}_\alpha$——兩者都在 $\operatorname{span}(V_r)$ 內，夾角才有意義——並同時報告 $\tilde\nu_\alpha$ 以說明這個比較涵蓋了 alpha 的多少比例。
> - **目標對齊**（§45）使用 $c^{\rm full}_\alpha$ 的 $\mathsf D$ 分量，因為目標方向 $r_{\mathsf D}$ 定義在完整的 $\mathsf D$ 空間中，與 $V_r$ 無關。

### §43.6 這個分解的邊界

$Q^{\rm nov}$ 量測的是「在 $\operatorname{span}(\Phi)$ 之內、但在 $\operatorname{span}(V_r)$ 之外」的預測幾何。**完全落在 $\operatorname{span}(\Phi)$ 之外的東西仍然不可見**——那需要擴充 §8 的 lifting 族，不是這個分解能處理的。

因此正確的陳述是：

$$\nu_\alpha\ \text{量測的是「相對於 feature library 的新穎性」，不是「絕對的新穎性」}.$$

---

## §44. Alpha Quality 與 fingerprint

### 定義

$$\boxed{Q_\alpha=Q^{\rm full}_\alpha=\|c^{\rm full}_\alpha\|^2,}$$

各 block 分量 $Q^{\mathsf D}_\alpha,Q^{\mathsf V}_\alpha,Q^{\mathsf L}_\alpha$ 同樣取自 $c^{\rm full}_\alpha$，並滿足加總恆等式。

$$\boxed{f_\alpha=\left(\frac{Q^{\mathsf D}_\alpha}{Q_\alpha},\frac{Q^{\mathsf V}_\alpha}{Q_\alpha},\frac{Q^{\mathsf L}_\alpha}{Q_\alpha}\right),}$$

以及 enrichment 版本 $\tilde f_\alpha$（§19.3，虛無基準由 §46 的 alpha permutation 給出）。

### $\tilde f_\alpha$ 是最重要的單一診斷

$\tilde f_\alpha$ 回答一個在傳統 IC 分析中完全看不到的問題：**模型的預測力屬於哪一類？**

| $\tilde f_\alpha$ 的主導成分 | 意義 | 部位設計的後果 |
|---|---|---|
| $\mathsf D$ | 模型預測方向 | 可直接建構多空排序部位 |
| $\mathsf V$ | 模型預測活躍度 | 多空排序會變成做多 gamma，報酬特性完全不同 |
| $\mathsf L$ | 模型預測大行情中的方向 | 需要條件性部位規模；平時應降低曝險 |

一個 $\tilde f_\alpha=(0.2,0.7,0.1)$ 的模型，在傳統 IC 上可能看起來不錯（因為波動度與報酬偏態有關聯，會洩漏出正的 rank IC），但它實際上不是一個方向性模型。把它當方向性模型使用，其實現的報酬分布會與回測嚴重不符。

**這是本方法最直接的實務價值**：它讓「模型看起來有 IC 但實際部位表現不如預期」這個常見現象，在部署之前就被診斷出來。

### 量級參考

$Q_\alpha$ 與 $Q_p$ 使用同一個尺規（§11）。一個良好的組合模型，$Q^{\mathsf D}_\alpha$ 應顯著高於最佳單一 feature 的 $Q^{\mathsf D}_p$——若不然，模型並沒有從組合中獲得任何好處，這本身就是一個結論。典型值：最佳單一 feature $Q^{\mathsf D}_p\approx5\times10^{-4}$（IC 0.022），良好模型 $Q^{\mathsf D}_\alpha\approx2\times10^{-3}$（IC 0.045）。

---

## §45. Alpha target alignment

### 定義

$$\boxed{A^{\rm tgt}_\alpha=\big(c^{{\rm full},\mathsf D}_\alpha\big)^\top r_{\mathsf D},\qquad \eta^{\mathsf D}_\alpha=\frac{(A^{\rm tgt}_\alpha)^2}{Q^{\mathsf D}_\alpha+\epsilon}\in[0,1].}$$

### 三者構成的完整方向性診斷

$$Q^{\mathsf D}_\alpha\ \text{（方向性資訊有多少）},\quad A^{\rm tgt}_\alpha\ \text{（方向對不對）},\quad \eta^{\mathsf D}_\alpha\ \text{（有多少對準目標）}.$$

這組診斷能區分三種在單一 IC 下看起來相同的失效：

| 情形 | $Q^{\mathsf D}_\alpha$ | $\eta^{\mathsf D}_\alpha$ | $A^{\rm tgt}_\alpha$ | 診斷 |
|---|---:|---:|---:|---|
| 健康 | 高 | 高 | 正、大 | 模型方向性強且對準目標 |
| 錯位 | 高 | 低 | 小 | 模型有方向性，但預測的是另一組 horizon |
| 反向 | 高 | 高 | 負、大 | 模型系統性反向——檢查標籤符號或訓練目標設定 |

第二種情形在多 horizon 訓練中很常見：模型的損失函數若對所有 horizon 等權，而 $\Sigma_Y$ 的譜差異很大，模型會自動偏向變異最大的 horizon，未必是使用者想交易的那個。修正方式是調整訓練目標的 horizon 權重，使之與 $a$ 一致。

---

## §46. Alpha 的顯著性

### 規格

完全沿用 §13 的 permutation：把 $X_d$ 換成 $\tilde\alpha_d$，以 circular date shift 破壞配對，套用 §13.3 的 matched thinning 與交集重新標準化，得到 $\{c^{\rm full,null}_{\alpha,b}\}$。

由此計算全部六個 p-value：

$$p^{\rm full}_\alpha,\quad p^{\mathsf D}_\alpha,\quad p^{\mathsf V}_\alpha,\quad p^{\mathsf L}_\alpha,\quad p^{\rm tgt}_\alpha\ \big((A^{\rm tgt}_\alpha)^2\big),\quad p^{\rm nov}_\alpha\ \big(Q^{\rm nov}_\alpha\big).$$

$Q^{\rm nov}_\alpha$ 的檢定**要求 $V_r$ 與評估資料獨立**（§15 的切分規格），否則 $I-P_V$ 的方向本身就是由同一份資料挑選剩下的，null 會偏斜。這是切分規格在 alpha 分析中的直接回報。

### 模型多重性

對**單一**模型不需要多重檢定調整。但實務上研究者幾乎總是比較多個模型／超參數，此時：

> **規則 46.1**：記錄 `n_model_variants_evaluated`——在同一份評估資料上被計算過 $Q_\alpha$ 的模型變體總數（包含被放棄的）。若 $>1$，必須對模型族做 Westfall–Young 調整，或保留一段完全獨立的 final test 並只用一次。

這個欄位常被省略，而它是 alpha 研究中最大的隱性多重檢定來源。一個研究者試了 50 個超參數組合並報告最好的那個，其有效顯著水準比名目值差了一到兩個數量級。

---

## §47. Alpha signal breadth

### 定義

這一節量測的是**模型輸出本身**的集中度，與預測力無關：

$$u^{\rm sig}_{id}=\frac{|\tilde\alpha_{id}|}{\sum_j|\tilde\alpha_{jd}|},\qquad\boxed{N^{\rm eff,sig}_d=\frac1{\sum_i(u^{\rm sig}_{id})^2},\qquad b^{\rm sig}_d=\frac{N^{\rm eff,sig}_d}{N_d}.}$$

### 解讀

$b^{\rm sig}$ 低代表模型的輸出高度集中在少數名次上——多數股票的預測值接近零。這本身無所謂好壞，但它決定了部位結構：

- $b^{\rm sig}\approx0.4$（常態水準）：模型對整個 universe 都有意見，部位會廣泛分散。
- $b^{\rm sig}<0.1$：模型實質上是一個「選股器」而非「排序器」，只對少數名次有強烈意見。此時應考慮直接以 top/bottom-$k$ 建構部位，而非線性加權。

與 §31 同理，虛無基準不是 1。若 $\tilde\alpha$ 在橫截面上近似常態，$b^{\rm sig}\to2/\pi\approx0.637$（單一常態量，而非兩個量的乘積）。

---

## §48. Alpha predictive breadth

### 定義

與 §30 完全同構，但作用於 alpha：

$$g_{id,\alpha}=\tilde\alpha_{id}\big\langle\Phi_{id,:},\hat c^{\rm full}_\alpha\big\rangle,\qquad\boxed{N^{\rm eff,pred}_{d,\alpha}=\frac1{\sum_i(u^{\rm pred}_{id,\alpha})^2},\quad b^{\rm pred}_{d,\alpha}=\frac{N^{\rm eff,pred}_{d,\alpha}}{N_d}.}$$

也可以限制在 $\mathsf D$ block 上計算 $b^{{\rm pred},\mathsf D}_{d,\alpha}$。

### 關鍵對比

$b^{\rm sig}$ 與 $b^{\rm pred}$ 的差異是一個重要診斷：

| $b^{\rm sig}$ | $b^{\rm pred}$ | 意義 |
|---|---|---|
| 高 | 高 | 模型廣泛下注，預測力也廣泛分布。容量大，穩健。 |
| 高 | 低 | **模型廣泛下注，但預測力只來自少數名次。** 絕大部分部位是無效的——它們貢獻了成本與風險，卻沒有貢獻報酬。 |
| 低 | 低 | 模型集中下注且預測力集中。一致，但容量小、對單一事件敏感。 |
| 低 | 高 | 罕見。模型只對少數名次下注，但預測力分散——通常代表 $\hat c_\alpha$ 的估計不穩定。 |

第二列是最有價值的發現，因為它指出一個可以直接改善的問題：對預測力貢獻接近零的名次施加部位截斷或分數門檻，可以在幾乎不損失報酬的情況下大幅降低成本與換手。

---

## §49. 傳統 IC 必須保留

### 規格

計算目標訊號 $y^*_d=Y_d^{(\beta)}a$，每日

$$IC_d=\operatorname{Spearman}(\tilde\alpha_d,y^*_d),$$

報告

$$\overline{IC},\quad \operatorname{median}(IC_d),\quad P(IC_d>0),\quad ICIR=\frac{\operatorname{mean}(IC_d)}{\operatorname{std}(IC_d)}.$$

### 為什麼不能省略

本方法的幾何量與 IC 是互補的，不是替代的：

1. **IC 是 rank-based，對單調變換免疫。** 幾何量是線性的。兩者不一致本身就是資訊：若 $A^{\rm tgt}_\alpha$ 強而 IC 弱，代表模型的預測力集中在尾部（線性相關被少數極端值撐起來，但排序沒有改善）；反之則代表關係單調但非線性。
2. **IC 直接對應可交易的排序表現**，是與外部（其他團隊、其他方法論）溝通的共同語言。
3. **IC 有大量既有的經驗基準**，而 $Q$ 的基準需要在本方法內部建立。

$ICIR$ 還提供一個 $Q$ 無法提供的量：預測力的**時間穩定度**（以日為單位），與 §27 的 bootstrap 穩定度互補。

### 技術註記

$IC_d$ 的序列高度自相關（重疊 label 窗口），因此 $ICIR\cdot\sqrt{T}$ **不是**一個有效的 t 統計量。若需要 IC 的顯著性，使用 §13 的同一組 permutation，而非解析式。

---

## §50. Feature 與 Alpha 必須在同一個 OOS window 比較

### 規格

對 OOS period，使用**同一個** $\Phi_{\rm OOS}$（由 §10 transform 得到）與**同一個** $V_r$，重新計算

$$\boxed{c_p^{\rm OOS}=\left[\frac1{T_{\rm OOS}}\sum_{d}\frac1{N_d}X^{\rm OOS\top}_{d,p}\Phi^{\rm OOS}_d\right],\qquad c_\alpha^{\rm OOS}=\tilde\alpha^{\rm OOS\top}W_{\rm OOS}\Phi^{\rm OOS}.}$$

Feature–Alpha 的所有比較（§51–§53）必須使用這兩組 OOS 表徵。

### 為什麼這條規則不可妥協

若拿研究期的 $c_p$ 與 OOS 的 $c_\alpha$ 比較，兩者來自不同的市場期間、不同的報酬實現、不同的 regime。此時「模型繼承了哪個 feature 的機制」這個問題的答案會被兩個無關的因素污染：

1. feature 自己的機制在兩期之間的漂移；
2. label 共變異結構的漂移（$\delta_{\rm drift}$）。

兩者都與「模型學了什麼」無關，但都會改變夾角。在同一個 window 上計算，這兩個污染源被完全消除——剩下的差異只能來自模型本身。

### 附帶的 feature 衰減診斷

同一計算免費提供一個重要資訊：比較 $c_p^{\mathcal E}$ 與 $c_p^{\rm OOS}$，得到每個 selected feature 的 OOS 機制穩定度

$$M^{\rm OOS}_p=\frac{c_p^{{\rm OOS},\mathsf D\top}c_p^{\mathcal E,\mathsf D}}{\|c_p^{{\rm OOS},\mathsf D}\|\|c_p^{\mathcal E,\mathsf D}\|}.$$

若模型表現衰退而 $M^{\rm OOS}_p$ 普遍偏低，問題出在 feature 而非模型。這個區分對後續的研究方向有決定性影響。

---

## §51. Feature–Alpha 第一層：Raw output exposure

### 定義

$$\rho_{p,\alpha,d}=\operatorname{Spearman}\big(X_{d,p},\tilde\alpha_d\big),\qquad\boxed{\rho^{\rm rank}_{p,\alpha}=\operatorname{median}_d\rho_{p,\alpha,d}.}$$

同時可計算 Pearson 版本 $\rho^{\rm lin}_{p,\alpha}$。

### 這一層在問什麼

> 模型的輸出排序，有多接近某個原始 feature 的排序？

這是**輸出層面**的相似度，與預測機制無關。它的用途是：

- 偵測模型是否退化為單一 feature 的包裝（$\rho^{\rm rank}>0.9$）。
- 理解部位的實際構成：若模型輸出與市值高度相關，部位會有系統性的 size 曝險，無論模型的預測機制是什麼。

$\rho^{\rm rank}$ 與後續兩層是**不同的問題**，不應混用（§53）。

---

## §52. Feature–Alpha 第二層：Mechanism inheritance

### 定義

依約定 43.2，在 $\operatorname{span}(V_r)$ 內比較：

$$\hat c_p=\frac{c_p}{\|c_p\|},\qquad \hat c^{\rm inh}_\alpha=\frac{c^{\rm inh}_\alpha}{\|c^{\rm inh}_\alpha\|},$$

$$\boxed{A^{\rm mech}_{p,\alpha}=\hat c_p^\top\hat c^{\rm inh}_\alpha,\qquad S^{\rm mech}_{p,\alpha}=\big|A^{\rm mech}_{p,\alpha}\big|,}$$

以及 block-specific 版本

$$\boxed{A^{\mathsf D}_{p,\alpha}=\frac{c_p^{\mathsf D\top}c_\alpha^{{\rm inh},\mathsf D}}{\|c_p^{\mathsf D}\|\,\|c_\alpha^{{\rm inh},\mathsf D}\|}.}$$

方向性研究一律以 $A^{\mathsf D}_{p,\alpha}$ 為主（規則 24.1）。

### 這一層在問什麼

> 模型繼承的預測機制，與哪一個 feature 的機制最接近？

這是三層中最直接回答「模型學到了哪個 feature」的量。

### 必須同時呈現的三件事

**(a) 涵蓋比例。** $A^{\rm mech}_{p,\alpha}$ 只描述 alpha 的 inherited 成分。必須同時報告 $\tilde\nu_\alpha$，否則讀者會把「模型的機制像 feature X」誤讀為「模型就是 feature X」，而實際上可能有大半的預測幾何落在 novel 成分中。

**(b) 維度校準。** §23.3 的低維限制在這裡完全適用。$m_{\mathsf D}=3$ 時，一個隨機方向與任一 feature 的 $|\cos|>0.9$ 的機率就有 10%；$P$ 個 feature 中至少有一個超過 0.9 的機率接近 1。因此

> **$A^{\mathsf D}_{p,\alpha}=0.92$ 這個數字本身不足以宣稱「模型繼承了 feature $p$」。**

必須比較的是**相對排序與差距**：$\alpha$ 與所有 selected feature 的 $|A^{\mathsf D}|$ 排序，以及第一名與第二名的差距。若前五名都在 0.85–0.93 之間，那就是低維空間的必然結果，不是繼承關係。同時報告虛無分布（以 §23.2 的獨立位移產生 $\alpha$ 的 null）下的 $\max_p|A^{\mathsf D}_{p,\alpha}|$ 分位數。

**(c) 組合而非單一 feature。** 模型通常繼承的是 feature 的**組合**而非單一 feature。更有資訊量的量是：alpha 的 inherited 成分能被 selected feature 的機制張成多少：

$$\boxed{\chi_\alpha=\frac{\big\|c^{{\rm inh},\mathsf D}_\alpha P^{\mathsf D}_{\mathcal F}\big\|^2}{\big\|c^{{\rm inh},\mathsf D}_\alpha\big\|^2}\in[0,1],}$$

其中 $P^{\mathsf D}_{\mathcal F}$ 是 selected feature 的 $\mathsf D$ 機制向量所張成子空間的投影（§38.1）。注意由性質 38.2，當 $|\mathcal F|\ge K$ 且機制向量線性獨立時 $\chi_\alpha\equiv1$，因此 $\chi_\alpha$ 只在 $|\mathcal F|<K$ 或機制向量退化時提供資訊。這正是為什麼「模型繼承了什麼」這個問題最終必須靠 §54–§57 的線性代理來回答，而非靠幾何投影。

---

## §53. 第一層與第二層的聯合判讀

$\rho^{\rm rank}_{p,\alpha}$（輸出相似）與 $A^{\mathsf D}_{p,\alpha}$（機制相似）構成一個 $2\times2$：

| | $A^{\mathsf D}$ 高 | $A^{\mathsf D}$ 低 |
|---|---|---|
| $\rho^{\rm rank}$ 高 | **模型≈該 feature 的單調變換。** 檢查模型是否只是在包裝單一 feature；若是，其複雜度不值得。 | **輸出像但機制不同。** 通常代表 $\rho$ 被共同的橫截面結構（產業、市值）撐起來，而真正的預測內容不同。應在控制已知暴險後重算 $\rho$。 |
| $\rho^{\rm rank}$ 低 | **機制相同但輸出不同。** 這是最有價值的情形：模型以不同的方式（不同的名次權重、不同的非線性）達成同一個預測機制，可能有更好的雜訊特性或成本特性。 | **兩者無關。** 若此格對所有 feature 都成立，而 $\tilde\nu_\alpha$ 又高，則模型的預測力來自 library 之外——進入 §58 的殘差分析。 |

左下角值得特別說明。它的典型成因是：模型學到了「在某些條件下才使用 feature $p$」。輸出排序因而與 $X_p$ 不同（條件不滿足時模型不下注），但預測機制相同（條件滿足時預測的是同一件事）。這類模型往往有較好的 $b^{\rm pred}/b^{\rm sig}$ 比值（§48），是複雜度真正帶來價值的情形。

---
# 第七部分：線性代理、殘差與模型診斷

## §54. 第三層的問題

前兩層（§51 輸出相似、§52 機制相似）都是幾何的、相關性的。它們無法回答一個關鍵問題：

> **模型的複雜度有沒有換來任何東西？**

這個問題的正確對照組是：**同樣的輸入、線性的模型**。若一個 ridge 迴歸可以在 OOS 上重建 95% 的模型輸出，那麼模型的非線性只是在複製一個線性函數，複雜度沒有帶來價值——而且部署 ridge 的運維成本、可解釋性、穩定性都遠優於部署原模型。

因此第三層的設計是：建立一個 OOS linear surrogate，量測它能重建多少，並分析**重建不了的那部分**是否有預測力。

---

## §55. Walk-forward linear surrogate

### 規格

把 OOS alpha period 依時間切成連續的 folds。對第 $k$ 個 fold，在 $Train_k^{\rm sur}$ 上擬合、在 $Test_k^{\rm sur}$ 上預測：

$$\hat\beta_k=\arg\min_\beta\left[\sum_{(d,i)\in Train_k^{\rm sur}}w_{id}\big(\tilde\alpha_{id}-X_{id}\beta\big)^2+\lambda\|\beta\|_2^2\right],\qquad \hat\alpha^{\rm lin}_{id}=X_{id}\hat\beta_k\ \ \text{for }(d,i)\in Test_k^{\rm sur}.$$

所有 fold 的預測拼起來得到完整的 OOS linear surrogate $\hat\alpha^{\rm lin}$。

### 四個必須明確的規格點

**(a) 用哪一組 $X$。** 主要報告使用**模型的實際輸入集合**（§39 選出的 `selected`）。理由：問題是「模型的非線性有沒有價值」，其正確對照組是「同樣的輸入、線性的模型」。

同時附帶報告使用**全體合格 feature** 的版本，回答一個不同的問題：「這個 alpha 是否能被整個 library 線性複製？」兩個數字分別命名為

$$R^2_{\rm lin,selected},\qquad R^2_{\rm lin,library}.$$

$R^2_{\rm lin,library}\gg R^2_{\rm lin,selected}$ 代表模型使用了未被選入的 feature 所攜帶的資訊——這通常意味著 §39 的選擇過於嚴格，或模型輸入與 selected 集合不一致（一個應該被排除的實作錯誤）。

**(b) 標準化。** $X$ 使用 §4.1 處理後的 $X^{XS}$（已 within-date 標準化並全域正規化），與整個管線一致。$\tilde\alpha$ 使用 §42 的版本。**不得**在 surrogate 訓練中重新標準化，否則 $R^2$ 與 $Q$ 的尺度不再對應。

**(c) 截距。** $X^{XS}$ 與 $\tilde\alpha$ 都已 within-date 去均值，因此加權最小平方的截距恆為零。**不擬合截距**。這一點必須在程式中明確關閉（多數套件預設開啟），否則會引入一個吸收日間共同成分的自由度。

**(d) $\lambda$ 的決定。** 在 $Train_k^{\rm sur}$ **內部**以時間序 CV（依時間切，不可隨機切）決定，絕不使用 $Test_k^{\rm sur}$ 的任何資訊。每個 fold 各自決定 $\lambda_k$，並在報表中記錄 $\{\lambda_k\}$——若 $\lambda_k$ 在 fold 之間變動超過一個數量級，代表 surrogate 本身不穩定，後續的 $R^2$ 解讀應降級。

### 為什麼用 ridge 而非 OLS

selected feature 之間必然相關（§36 只移除了 $S^{\rm raw}>0.9$ 的重複，$0.6$–$0.9$ 的相關普遍存在）。OLS 的係數在這種條件下極不穩定，$\hat\beta_k$ 會在 fold 之間劇烈跳動，使 surrogate 的 OOS 表現被係數雜訊而非模型結構決定。ridge 以少量偏誤換取大幅的變異縮減，對「能不能重建」這個問題給出更穩定的答案。

注意 surrogate 的目的**不是**估計每個 feature 的貢獻（那由 §57 處理），而是估計**可重建的比例**。因此係數本身不解讀。

---

## §56. Reconstruction metrics

### 定義

$$\boxed{R^2_{\rm lin,OOS}=1-\frac{\sum_{d,i}w_{id}\big(\tilde\alpha_{id}-\hat\alpha^{\rm lin}_{id}\big)^2}{\sum_{d,i}w_{id}\tilde\alpha_{id}^2}=1-\big\|\varepsilon_\alpha\big\|_W^2,}$$

（末式用到 $\|\tilde\alpha\|_W=1$），以及排序層面的

$$\rho^{\rm rank}_{\rm lin}=\operatorname{median}_d\operatorname{Spearman}\big(\tilde\alpha_d,\hat\alpha^{\rm lin}_d\big).$$

### 兩個指標的分工

$R^2$ 量測**數值**重建，$\rho^{\rm rank}$ 量測**排序**重建。兩者可以差很多：

- $R^2=0.5$ 但 $\rho^{\rm rank}=0.9$：模型與 surrogate 的排序幾乎相同，差異在於尾部的幅度。由於部位通常由排序決定，這種情形下模型的實際價值接近 surrogate。
- $R^2=0.8$ 但 $\rho^{\rm rank}=0.6$：數值接近但排序不同——通常代表模型在中段（多數名次）與 surrogate 分歧，而中段恰好是部位權重最密集的區域。

**報表必須兩者並列**，並以 $\rho^{\rm rank}_{\rm lin}$ 作為「模型是否可被 surrogate 取代」的主要判準，因為它更接近交易語意。

### OOS $R^2$ 可以是負的

若 surrogate 在某些 fold 上嚴重失效，$R^2_{\rm lin,OOS}<0$ 是可能的，代表 surrogate 比「全部預測 0」還差。此時不應解讀為「模型非常非線性」，而應先檢查 surrogate 的實作（$\lambda$ 是否過小、fold 切分是否有洩漏或斷裂）。

---

## §57. 哪些 feature family 在重建 Alpha

### 規格

單一 ridge 係數不適合作為重要性度量（相關 feature 會使係數不穩定且可任意互換）。改用兩個互補的量，**必須同時報告**。

**(a) 條件貢獻（drop-family）**：對 redundancy family $R_g$，移除該 family 的全部 feature，以完全相同的流程重跑 surrogate，得到 $R^2_{-g}$：

$$\boxed{\Delta R^2_g=R^2_{\rm full}-R^2_{-g}.}$$

語意：**在其他所有 family 都在場的條件下**，這個 family 的邊際貢獻。

**(b) 獨立貢獻（standalone）**：只用 family $R_g$ 的 feature 擬合 surrogate：

$$\boxed{R^2_{g,{\rm own}}.}$$

語意：這個 family 單獨能重建多少。

### 為什麼兩個都需要

兩個 family 若彼此高度相關（即使 $S^{\rm raw}<0.9$ 因而未被合併），會互相遮蔽：

$$\Delta R^2_g\approx0\quad\text{且}\quad\Delta R^2_{g'}\approx0,\qquad\text{但}\quad R^2_{g,{\rm own}}=R^2_{g',{\rm own}}=0.6.$$

只看 $\Delta R^2$ 會得出「兩個 family 都不重要」的結論，而事實是「兩者擇一即可，但缺一不可」。

判讀規則：

| $\Delta R^2_g$ | $R^2_{g,\rm own}$ | 判讀 |
|---|---|---|
| 高 | 高 | 核心 family，不可替代 |
| 低 | 高 | 與其他 family 重複；可作為替代品，提供備援 |
| 高 | 低 | 只在與其他 family 組合時有用（交互作用或殘差修正） |
| 低 | 低 | 對線性重建無貢獻 |

### 更完整的版本

若 family 數不多（$\le12$），可計算 family 層級的 Shapley 值：對所有子集計算 $R^2$，取邊際貢獻的平均。成本為 $2^G$ 次 surrogate 擬合，$G=10$ 時約 1000 次，實務上可接受。Shapley 值同時解決了上述的遮蔽問題，是最乾淨的歸因，但不可擴展到大 $G$。

---

## §58. Residual alpha

### 定義

$$\boxed{\varepsilon_\alpha=\tilde\alpha-\hat\alpha^{\rm lin}.}$$

由於 $\tilde\alpha$ 與 $X^{XS}$ 都已 within-date 去均值，$\varepsilon_\alpha$ 自動滿足 $\sum_i\varepsilon_{\alpha,id}=0$。

**兩個尺度必須分開記錄。**

**絕對殘差能量**：

$$s^2_\varepsilon=\|\varepsilon_\alpha\|_W^2=1-R^2_{\rm lin,OOS}.$$

**標準化殘差訊號**：對 $\varepsilon_\alpha$ 做 within-date 標準化並全域正規化，得到 $\tilde\varepsilon_\alpha$ 滿足 $\|\tilde\varepsilon_\alpha\|_W=1$，再依 §43 計算

$$c^{\rm full}_\varepsilon=\tilde\varepsilon_\alpha^\top W\Phi,\qquad c^{\rm inh}_\varepsilon=c^{\rm full}_\varepsilon P_V,\qquad c^{\rm nov}_\varepsilon=c^{\rm full}_\varepsilon(I-P_V),$$

$$Q_\varepsilon=\|c^{\rm full}_\varepsilon\|^2,\quad Q^{\mathsf D}_\varepsilon,\ Q^{\mathsf V}_\varepsilon,\ Q^{\mathsf L}_\varepsilon,\quad A^{\rm tgt}_\varepsilon,\ \eta^{\mathsf D}_\varepsilon,\quad\tilde\nu_\varepsilon.$$

### 為什麼必須標準化後再計算 $Q_\varepsilon$

若直接用原尺度的 $\varepsilon_\alpha$，則 $\|\varepsilon\|_W^2=1-R^2$，因此

$$Q_\varepsilon^{\rm raw}\le1-R^2.$$

在 $R^2=0.9$ 時 $Q_\varepsilon^{\rm raw}$ 的上界只有 0.1 倍於 $Q_\alpha$ 的上界。把兩個數字直接比較會得出「殘差沒什麼預測力」的結論，而那只是尺度的算術，與殘差的品質無關。

標準化之後，$Q_\varepsilon$ 與 $Q_\alpha$ 在**同一個尺規**上（§11 的 IC 尺規），可以直接比較。兩者的關係為

$$\boxed{\underbrace{s^2_\varepsilon\cdot Q_\varepsilon}_{\text{殘差對總預測能量的絕對貢獻}}\quad\text{vs}\quad\underbrace{Q_\alpha}_{\text{總預測能量}}.}$$

**報表必須同時列出 $Q_\varepsilon$（品質）與 $s^2_\varepsilon Q_\varepsilon$（絕對貢獻）。** 前者回答「殘差訊號好不好」，後者回答「它在整體中佔多少」。

### 顯著性

以 §13 的同一套 permutation 對 $\tilde\varepsilon_\alpha$ 做檢定，得到 $p^{\mathsf D}_\varepsilon$、$p^{\rm nov}_\varepsilon$ 等。

值得注意的是：即使 $Q_\varepsilon$ 的**量級**因標準化而變得可比，permutation p-value 本身從來就是尺度不變的——它比較的是同一個標準化訊號的觀測值與其虛無分布。因此顯著性結論在標準化前後完全相同，標準化改變的只是「可比較的量級」。

### 這一節是整個方法的最終產物

$Q^{\mathsf D}_\varepsilon$ 與 $\tilde\nu_\varepsilon$ 的組合回答了一個在量化研究中極難回答、而本方法可以直接回答的問題：

> **模型超出線性 feature span 的那部分輸出，究竟有沒有預測價值？如果有，它預測的是已知機制的重新加權，還是一個全新的機制？**

$Q^{\mathsf D}_\varepsilon$ 顯著回答前半，$\tilde\nu_\varepsilon$ 回答後半。

---

## §59. 三層的邏輯關係

$$\underbrace{\rho^{\rm rank}_{p,\alpha}}_{\text{第一層：輸出}}\ \longrightarrow\ \underbrace{A^{\mathsf D}_{p,\alpha},\ \tilde\nu_\alpha}_{\text{第二層：機制}}\ \longrightarrow\ \underbrace{R^2_{\rm lin},\ Q_\varepsilon,\ \tilde\nu_\varepsilon}_{\text{第三層：超線性結構}}$$

三層是**逐步收緊**的：

- 第一層可以在沒有 $\Phi$ 的情況下計算，但只能看到表面相似。
- 第二層需要 $\Phi$ 與 $V_r$，能看到預測內容，但看不到「複雜度是否必要」。
- 第三層需要額外擬合 surrogate，能看到複雜度的邊際價值，但只在線性的對照下有意義。

**三層不可互相取代。** 一個常見的錯誤是只看第三層：$R^2_{\rm lin}=0.3$ 就宣稱「模型高度非線性且有價值」。但若 $Q_\varepsilon$ 不顯著，那 70% 的不可重建部分只是雜訊，模型實際上比 surrogate 更差（§60 Case C）。

另一個常見錯誤是只看第二層：$A^{\mathsf D}_{p,\alpha}=0.95$ 就宣稱「模型學的就是 feature $p$」，而忽略了 $\tilde\nu_\alpha$ 可能很高（模型大部分的幾何不在 $V_r$ 內）、以及低維空間的巧合效應（§52(b)）。

### 「非線性」的兩種完全不同的意思

本方法把一個長期含混的概念拆成兩個可量測的量：

$$\boxed{\text{輸出非線性}:\ 1-R^2_{\rm lin,OOS}\ \text{大}}\qquad\text{vs}\qquad\boxed{\text{預測上新穎的非線性}:\ Q^{\mathsf D}_\varepsilon\ \text{顯著}}$$

前者只說明「模型的輸出不是輸入的線性函數」，這對任何樹模型或神經網路都自動成立，沒有資訊量。後者才是我們真正想知道的：**那些非線性有沒有換來新的預測力。**

---

## §60. 五種診斷結果

以下五種情形涵蓋了實務上絕大多數的模型診斷結論。

### Case A：模型等價於線性代理

$$R^2_{\rm lin,OOS}=0.90,\quad \rho^{\rm rank}_{\rm lin}=0.94,\quad Q^{\mathsf D}_\varepsilon\ \text{不顯著}.$$

**結論**：alpha 幾乎完全可以由 feature 的穩定線性組合重建，超出線性 span 的部分沒有顯著的方向性預測價值。

**處置**：部署 ridge surrogate，不部署原模型。

這**不是**模型的失敗，而是一個極有價值的明確答案——它為團隊節省大量運維成本，移除一個過擬合風險來源，並大幅提升可解釋性與可監控性。能夠確信地做出這個決定，本身就是本方法的主要回報之一。

### Case B：非線性有價值

$$R^2_{\rm lin,OOS}=0.40,\quad Q^{\mathsf D}_\varepsilon\ \text{顯著},\quad A^{\rm tgt}_\varepsilon>0,\quad \tilde\nu_\varepsilon\approx0.4.$$

**結論**：模型產生了大量不能被穩定線性 feature span 重建、且對 intended target 仍具正向方向性預測價值的結構；該結構主要落在已知機制的重新加權上（$\tilde\nu_\varepsilon$ 不高）。

**處置**：部署原模型。接著用 §57 的 family 歸因與條件分析，試圖理解這個非線性結構的形式——若能理解（例如「模型在高波動日才使用動能 feature」），往往可以把它固化為新的 feature，在下一輪中變成線性可得的。這是把模型能力轉化為 library 資產的標準路徑。

### Case C：非線性是純粹的成本

$$R^2_{\rm lin,OOS}=0.30,\quad Q^{\mathsf D}_\varepsilon\ \text{不顯著}.$$

**結論**：模型輸出有大量非線性／非平穩結構，但這些結構沒有提供統計上可支持的方向性預測資訊。

**處置**：這是最糟的情形——模型引入了大量無用的變異，增加換手率、交易成本與不可預測性，而沒有補償的報酬。應加強正則化，或直接退回 surrogate。

**進一步診斷**：檢查 $Q^{\mathsf V}_\varepsilon$。若 $Q^{\mathsf V}_\varepsilon$ 顯著而 $Q^{\mathsf D}_\varepsilon$ 不顯著，代表模型的非線性部分在預測波動度。這通常是損失函數造成的：以 MSE 對原始報酬訓練的模型會自動對高波動名次投入更多容量，因為它們貢獻最多的損失。修正方式是改用波動度正規化的訓練目標，或直接以 $y^*$（§20）為訓練標籤。

### Case D：強烈的機制繼承

$$A^{\mathsf D}_{p,\alpha}=0.92,\quad\text{且顯著高於次高者},\quad \Delta R^2_g\ \text{大},\quad R^2_{g,\rm own}\ \text{大}.$$

**結論**：final alpha 強烈繼承該 feature family 的預測機制，而且該 family 對 alpha 的線性重建也有重要的獨立貢獻。這是「模型使用這個機制」的強證據。

**必要的驗證**（§52(b)）：$A^{\mathsf D}=0.92$ 本身在 $m_{\mathsf D}=3$ 時並非罕見。此結論必須滿足三個條件才成立：(i) 與次高者的差距顯著；(ii) 超過虛無分布 $\max_p|A^{\mathsf D}_{p,\alpha}|$ 的 95 百分位；(iii) $\Delta R^2_g$ 與 $R^2_{g,\rm own}$ 同時支持。單憑幾何夾角不足以下這個結論。

### Case E：新機制

$$R^2_{\rm lin,OOS}=0.50,\quad Q^{\mathsf D}_\varepsilon\ \text{顯著},\quad \tilde\nu_\varepsilon>0.7,\quad p^{\rm nov}_\varepsilon<0.05.$$

**結論**：模型不但創造了超出線性 span 的方向性預測結構，而且該結構有大部分落在現有 feature library **完全照不到**的 label 方向上。

**處置**：這是最強的正面結果，也是最需要謹慎的——因為它宣稱發現了新機制。必要的後續驗證：

1. 檢視 $c^{\rm nov}_\varepsilon$ 在 $\Phi$ 中的具體載重：這是哪個 block、哪個 horizon 組合的方向？把它翻譯成經濟語言。
2. 在一段完全獨立、從未參與任何決策的 holdout 上重複驗證（且只用一次）。
3. 檢查 $\delta_{\rm drift}$：若座標系本身在漂移，$I-P_V$ 的方向可能只是反映了 $\Phi$ 的不穩定。
4. 嘗試構造能捕捉該方向的新 feature。若成功，該機制就從「黑箱模型的專屬能力」轉化為「library 的永久資產」，下一輪它會出現在 $V_r$ 中。

第 4 步是本方法設計的最終目的：**讓模型發現的東西可以被系統性地回收。**

---

## §61. Model report

每個 model × horizon band 輸出一列。

**強度與機制**

| 欄位 | 說明 |
|---|---|
| `Q_alpha_full`, `Q_alpha_D/V/L` | §44 |
| `Q_alpha_inh`, `Q_alpha_nov` | §43.1 |
| `f_alpha_D/V/L` | raw fingerprint |
| `ft_alpha_D/V/L` | enrichment fingerprint |
| `nu_alpha`, `nu_tilde_alpha` | §43.3、§43.4 |
| `A_target_alpha`, `eta_D_alpha` | §45 |

**顯著性**

| 欄位 | 說明 |
|---|---|
| `p_alpha_full`, `p_alpha_D/V/L` | §46 |
| `p_alpha_target`, `p_alpha_nov` | §46 |
| `n_model_variants_evaluated` | 規則 46.1 |
| `holdout_used_count` | 規則 41.1 |

**橫截面支撐**

| 欄位 | 說明 |
|---|---|
| `b_signal_med` | §47 |
| `b_pred_med`, `b_pred_D_med` | §48 |
| `b_pred/b_signal` | §48 的關鍵比值 |
| `R_cancel_alpha`, `R_tilde_alpha` | §32 套用於 alpha |

**傳統指標**

| 欄位 | 說明 |
|---|---|
| `IC_mean`, `IC_median`, `IC_hit_rate`, `ICIR` | §49 |

**線性代理與殘差**

| 欄位 | 說明 |
|---|---|
| `R2_lin_selected`, `R2_lin_library` | §55(a) |
| `rho_rank_lin` | §56 |
| `lambda_by_fold` | §55(d) |
| `s2_resid` ($=1-R^2$) | §58 |
| `Q_resid`, `Q_resid_D/V/L` | §58（標準化後） |
| `Q_resid_abs` ($=s^2_\varepsilon Q_\varepsilon$) | §58 |
| `A_target_resid`, `eta_D_resid` | §58 |
| `nu_tilde_resid` | §58 |
| `p_resid_D`, `p_resid_nov` | §58 |

**座標系健康度**

| 欄位 | 說明 |
|---|---|
| `delta_drift`, `delta_D/V/L` | §10 |
| `T_eff_D`, `T_eff_V`（OOS 期） | §3 |
| `diagnosis_case` | §60 的 A–E |

---

## §62. Feature–Alpha report

每個 selected feature 一列（全部在同一個 OOS window 上計算，§50）：

| 欄位 | 說明 |
|---|---|
| `rho_rank_p_alpha`, `rho_lin_p_alpha` | §51 |
| `A_mech_p_alpha`, `S_mech_p_alpha` | §52 |
| `A_D_p_alpha` | §52（方向性研究的主要欄位） |
| `A_D_rank` | 在所有 selected feature 中的排名 |
| `A_D_gap_to_next` | 與次高者的差距 |
| `A_D_null_q95` | 虛無下 $\max_p|A^{\mathsf D}_{p,\alpha}|$ 的 95 百分位 |
| `delta_R2_family` | §57(a) |
| `R2_own_family` | §57(b) |
| `M_OOS_p` | §50 的 feature 機制 OOS 穩定度 |
| `Q_p_OOS`, `A_target_p_OOS` | OOS 期的 feature 表徵 |
| `Q_decay` ($=Q_p^{\rm OOS}/Q_p^{\mathcal E}$) | feature 衰減 |
| `cell_2x2` | §53 的四格分類 |

**全域欄位**：$\tilde\nu_\alpha$、$\chi_\alpha$、$m_{\mathsf D}$（提醒讀者維度限制）、OOS window 的起訖日期。

---
# 第八部分：執行、驗證與統一觀點

## §63. 完整執行順序

以下順序是完整的、可直接作為實作檢查表。方括號標出對應章節。

### 階段 0：設定與資料

```text
 0. 鎖定設定：band 定義(累積/增量)、目標權重 a、Delta = H_max + l_X、
    切分比例、B_null / B_time / B_XS、全部門檻              [§6,§13,§20,§66]
 1. 依 (time, ticker) key 對齊；驗證 point-in-time 契約       [§2]
 2. 切分 discovery D 與 evaluation E，中間留 Delta 的 purge   [§15.2]
 3. Feature 前處理：winsorize -> within-date 標準化
    -> 全域正規化使 ||X_p||_W = 1                            [§4.1]
 4. 展開 transform family；記錄 (source, transform) 索引      [§4.3]
 5. 缺值處理；記錄 coverage_p 與 active_ratio_p               [§4.2]
 6. Label within-date 去均值；切 horizon band                 [§5,§6]
```

### 階段 1：在 discovery 上建立座標系

```text
 7. Sigma_Y；譜 bootstrap；決定 K_eff；記錄 lambda_j/lambda_1
    與 cond(Sigma_Y)                                          [§7.1,§7.2]
 8. B_Y -> Z；驗證 Z^T W Z = I 精確成立                        [§7.3,§7.4]
 9. Phi_D := Z（不做正交化）                                   [§9.2]
10. G_V, G_L -> within-date 中心化 -> Gram-Schmidt -> Phi      [§8.1,§9.3]
11. 儲存 B_Y, Theta_V, Theta_L, Gamma_DV, Gamma_DL, Gamma_VL   [§10]
12. 計算 T_eff^D, T_eff^V                                      [§3.1]
13. C_D；permutation null（共同位移，matched thinning，
    交集重新標準化）                                           [§11,§13]
14. parallel analysis -> r_sig, V_r（此後固定）                 [§14]
15. r_eff；m_D（有效機制維度）                                  [§33,§23.4]
```

### 階段 2：在 evaluation 上做全部推論

```text
16. transform label basis 到 E；計算 delta_drift(D->E)          [§10]
17. C_E -> c_p = C_E V_r V_r^T -> Q, Q_D, Q_V, Q_L, Q_full      [§15.2,§16]
18. kappa_p；決定每個 feature 使用 N1 或 N2 虛無                [§13.4]
19. permutation on E -> 原始 p 值 -> Westfall-Young 調整         [§18]
20. fingerprint f 與 enrichment f_tilde                          [§19]
21. r_D = normalize(Phi_D^T W y*)；A_target；eta_D；
    (A_target)^2 的顯著性                                        [§20,§21,§22]
22. 獨立位移 null -> 相似度虛無分布 -> zeta(tau)；S^D             [§23,§24]
23. temporal moving-block bootstrap -> pi, M, CV（block-specific）[§26,§27]
24. cross-sectional 分層 bootstrap（每次重抽後重新中心化）
    -> pi^XS, M^XS                                               [§29]
25. breadth 與 cancellation；相對化為 b_tilde, R_tilde            [§30,§31,§32]
```

### 階段 3：選擇

```text
26. eligibility E1-E5                                            [§34]
27. mechanism grouping（S^D, average linkage, tau=0.90）；
    報告 zeta 與 m_D                                             [§35]
28. raw redundancy families（complete linkage）；
    套用 transform family 合併規則                                [§36,§4.3]
29. Score_p = median_b (A_target^(b))^2 -> family representatives  [§37]
30. 貪婪最大化 Cov^feat_D 至 0.90 或 budget；再補 V, L             [§39]
31. 記錄 Cov^span 作為診斷                                        [§38.1]
32. 輸出 feature report                                           [§40]
```

### 階段 4：模型與歸因

```text
33. 訓練模型（purge + embargo = Delta；selected set 不得含 test 資訊）[§41]
34. alpha within-date 標準化 + 全域正規化                          [§42]
35. OOS: transform Phi；記錄 delta_drift                           [§10]
36. c_alpha^full -> inherited / novel 分解；Q, f, f_tilde, nu_tilde  [§43,§44]
37. A_target_alpha, eta_D_alpha                                    [§45]
38. alpha permutation -> p_full, p_D/V/L, p_target, p_nov           [§46]
39. b_signal, b_pred；IC / ICIR                                     [§47,§48,§49]
40. 在同一 OOS window 重算 c_p^OOS；M_OOS_p；Q_decay                 [§50]
41. rho_rank；A^D_{p,alpha} + 虛無 max 分位數；2x2 判讀              [§51,§52,§53]
42. walk-forward ridge surrogate -> R2_selected, R2_library,
    rho_rank_lin                                                    [§55,§56]
43. family 歸因：delta_R2 與 R2_own                                 [§57]
44. residual：s2_resid, Q_resid（標準化後）, Q_resid_abs,
    nu_tilde_resid, p_resid                                         [§58]
45. 判定 Case A-E                                                   [§60]
46. 輸出 model report 與 feature-alpha report                       [§61,§62]
```

### 階段 5（可選）：生產重擬合

```text
47. 在 D ∪ E 上重新擬合 label basis 與 V_r，僅供生產使用。
    重擬合後不得做任何新決策。                                      [§15.4]
```

---

## §64. Fit / transform 邊界總表

下表列出每一個被擬合的物件、它的擬合資料、以及套用方式。任何一行的違反都構成洩漏。

| 物件 | 擬合於 | 套用於 | 章節 |
|---|---|---|---|
| winsorize 分位數 | 逐期計算（無跨期擬合） | 當期 | §4.1 |
| within-date $\mu,\sigma$ | 逐期計算 | 當期 | §4.1 |
| 全域正規化 $\|X_p\|_W$ | $\mathcal D$ | 全期 | §4.1 |
| $\Sigma_Y$, $K_{\rm eff}$, $B_Y$ | $\mathcal D$ | 全期 | §7 |
| $G_\bullet$ 的 within-date 中心化 | 逐期計算 | 當期 | §8.1 |
| $\Theta_{\mathsf V},\Theta_{\mathsf L},\Gamma_\bullet$ | $\mathcal D$ | 全期 | §9 |
| $r_{\rm sig}$, $V_r$ | $\mathcal D$ | 全期 | §14,§15 |
| $r_{\mathsf D}$ | $\mathcal D$（或 $\mathcal E$，須一致並記錄） | 全期 | §20 |
| eligibility 清單 | $\mathcal E$ | 模型輸入 | §34 |
| grouping / families | $\mathcal E$ | 模型輸入 | §35,§36 |
| selected set | $\mathcal E$ | 模型輸入 | §39 |
| 模型參數 | $\subseteq\mathcal D\cup\mathcal E$，且早於 OOS | OOS | §41 |
| surrogate $\hat\beta_k,\lambda_k$ | 各 fold 的 train | 該 fold 的 test | §55 |

### 跨 fold 的向量不可直接平均

若採用完整的 walk-forward feature selection（規則 41.1 的完整版），每個 fold 有自己的 $\Phi^{(k)}$ 與 $V_r^{(k)}$。此時

$$\frac1K\sum_kc_p^{(k)}\quad\text{沒有意義}$$

——不同 fold 的座標系不同，逐分量相加是在把不相干的座標混在一起。

可以做的有三種：

1. **平均純量**：$\overline{Q_p}=\frac1K\sum_kQ_p^{(k)}$、$\overline{|A^{\rm tgt}_p|}$ 等。純量是座標無關的，可以平均。
2. **平均夾角**：$\overline{M_p}=\frac1K\sum_k\cos\big(c_p^{(k)},c_p^{(\rm ref)}\big)$，其中 $c_p^{(\rm ref)}$ 來自某個指定的參考 fold。
3. **投影到共同座標系**：選定一個參考 fold 的 $\Phi^{(\rm ref)}$，用 §10 的 transform 把所有 fold 的資料投影到它上面，再平均。這是唯一能得到「平均機制向量」的正確做法，代價是參考 fold 的選擇引入了一個任意性。

---

## §65. 數值驗證

分成兩組：**恆等式測試**（驗證實作正確）與**校準測試**（驗證統計有效）。兩組都必須通過。

### §65.1 恆等式測試

| 測試 | 判準 |
|---|---|
| Feature 去均值 | $\big|\frac1{N_d}\sum_iX^{XS}_{idp}\big|<10^{-8}$ |
| Feature 正規化 | $\big|\|X_p\|_W^2-1\big|<10^{-8}$ |
| Label 去均值 | $\big|\sum_iY^{XS}_{id,h}\big|<10^{-8}\cdot$ scale |
| 精確白化 | $\|Z^\top WZ-I\|_F<10^{-10}$（**必須是精確等式**，見命題 7.1） |
| 基底正交 | $\|\Phi^\top W\Phi-I\|_F<10^{-6}$（in-sample） |
| $C$ 的兩種算法一致 | long-form 直接算 vs $\frac1T\sum_dC^{(d)}$ |
| Quality 加總 | $Q_p=Q_p^{\mathsf D}+Q_p^{\mathsf V}+Q_p^{\mathsf L}$ |
| Novelty 加總 | $Q^{\rm full}_s=Q^{\rm inh}_s+Q^{\rm nov}_s$ |
| Efficiency 界 | $0\le\eta^{\mathsf D}\le1+\epsilon$ |
| 符號反轉 | $X_q=-X_p\Rightarrow S^{\rm mech}_{pq}=1,\ A^{\rm mech}_{pq}=-1$ |
| Coverage 單調 | ${\rm Cov}^{\rm feat}(\mathcal F\cup\{p\})\ge{\rm Cov}^{\rm feat}(\mathcal F)$ |
| Coverage 次模 | 隨機抽 100 組 $\mathcal A\subset\mathcal B$，驗證邊際增益遞減（定理 38.3） |

注意 Quality 加總測試是**恆真的**（命題 16.1，來自座標分割），它只能抓到程式的切片錯誤。真正驗證正交性的是精確白化與基底正交兩項。

### §65.2 校準測試

以合成資料執行。這一組測試才是判斷「統計結論是否可信」的直接證據。

**T1 虛無校準（$Q$）。** 生成獨立的 $X$ 與 $Y$（保留真實的橫截面相關結構與時間自相關，但兩者獨立），跑完整管線。

判準：$p$ 值大致均勻；$\gamma_{\rm FWER}=0.05$ 下的假陽性率 $\le5\%$。

*這是最重要的單一測試。* 若 §15 的切分或重跑規格實作錯誤（例如 null 仍使用觀測的 $V_r$），這個測試的假陽性率會顯著高於 5%。

**T2 虛無校準（相似度）。** 同樣的合成資料下，計算 $\zeta(\tau_{\mathsf D})$。

判準：$\zeta\approx0$（在 Monte Carlo 誤差內）。

**T3 Fingerprint 校準。** 純雜訊 feature 的 $\tilde f_p$。

判準：接近 $(1/3,1/3,1/3)$；三個分量的平均偏離 $<0.05$。

（對照：同一批 feature 的 raw $f_p$ 應接近維度比例 $(4,10,4)/18$，這正好驗證了 §19.2 的分析。）

**T4 旋轉不變性。** 以 10 組不同的 BLAS 執行緒數與隨機種子重跑整個管線。

判準：命題 9.2 所列的全部量在數值誤差（$10^{-10}$ 相對）內完全相同。

*必須重複多次執行。* 這個測試的失效模式是**間歇性**的——單次通過不代表正確。

**T5 Novelty 校準。** 生成與 $\Phi$ 無關的隨機 $\alpha$。

判準：$\tilde\nu_\alpha\approx0.5$；$p^{\rm nov}_\alpha$ 均勻。

（對照：raw $\nu_\alpha$ 應接近 $(d_\Phi-r_{\rm sig})/d_\Phi$，驗證 §43.4。）

**T6 Breadth / cancellation 校準。** 純雜訊 feature。

判準：$\tilde b_p\approx1$、$\tilde R_p\approx1$。

（對照：raw $b$ 應接近 0.405，raw $R^{\rm canc}$ 應接近 $\sqrt{\pi/2N}$。）

**T7 目標還原。** 令 $\alpha=y^*$。

判準：$A^{\rm tgt}_\alpha>0$ 且 $\eta^{\mathsf D}_\alpha\approx1$；令 $\alpha=-y^*$ 則 $A^{\rm tgt}_\alpha<0$。

**T8 植入訊號還原。** 構造 $X_p=\text{IC}\cdot y^*/\|y^*\|+\sqrt{1-\text{IC}^2}\cdot\text{noise}$，取 $\text{IC}\in\{0.01,0.03,0.05\}$。

判準：估計的 $|A^{\rm tgt}_p|$ 與植入的 IC 在 $\pm20\%$ 內一致。這同時驗證了 §11 的量級尺規。

**T9 動態 universe。** 隨機讓 $N_d$ 在 $\pm40\%$ 內變動、隨機增刪 ticker。

判準：全流程不依賴等列數；T1 的校準在動態 universe 下仍然成立（這是 §13.3 matched thinning 的驗證）。

**T10 洩漏偵測。** 把 label 整體往未來平移 $2\Delta$ 期（使真實關聯被破壞但所有結構保留）。

判準：全部 $Q$ 降到虛無水準。若不然，管線中存在洩漏。

**T11 慢變 feature。** 構造 $\kappa_p=0.95$ 的持續性 feature，且**與 label 無關**。

判準：N1（circular shift）的 $p$ 值明顯偏離均勻（這正是 §13.4 所預期的檢定力損失）；N2（分層 within-date permutation）的 $p$ 值均勻。這個對照驗證了雙虛無設計的必要性。

### §65.3 驗收標準

> 一個實作在 **T1、T2、T3、T4、T5** 全部通過之前，其輸出的任何顯著性結論都不應被採信。T1 驗證推論規格、T2 驗證相似度校準、T3 與 T5 驗證維度校準、T4 驗證座標系的良定義。這五項合起來涵蓋了本方法中全部會系統性改變結論的環節。

---

## §66. 預設參數

| 參數 | 預設值 | 章節 |
|---|---:|---|
| $X$ winsorization | 1% / 99% cross-sectional | §4.1 |
| 每日最低 feature 覆蓋率 | 0.80 | §4.2 |
| transform family | $\{x,|x|,x|x|\}$ | §4.3 |
| horizon band 慣例 | 不重疊增量（優先） | §6 |
| $\kappa_{\max}$（whitening 條件數上限） | 100 | §7.2 |
| $\tau_\lambda$（特徵值 bootstrap 相對 SE 上限） | 0.25 | §7.2 |
| $\epsilon_{\rm rank}$（$\operatorname{orth}_W$ 數值截斷） | $10^{-10}$ | §9.1 |
| discovery / evaluation 比例 | 60 / 40 | §15.2 |
| $\Delta$（permutation 最小位移＝purge gap） | $H_{\max}+\ell_X$ | §13.2,§41 |
| $\rho_{\min}$（交集最低比例） | 0.90 | §13.3 |
| $R_{\rm thin}$（matched thinning 重複次數） | 5 | §13.3 |
| $\kappa_p$ 門檻（改用 N2 虛無） | 0.50 | §13.4 |
| $B_{\rm null}$ | 2000（上限為全部合法位移） | §13.2,§18.3 |
| $\gamma_{\rm mode}$（mode 顯著水準） | 0.05 | §14 |
| 多重檢定主方案 | Westfall–Young，$\gamma_{\rm FWER}=0.05$ | §18.3 |
| 多重檢定次方案 | permutation FDR，$\gamma_{\rm FDR}=0.05$ | §18.4 |
| $B_{\rm time}$（時間 bootstrap） | 1000 | §26 |
| 區塊長度 $\ell_{\rm blk}$ | $\ge H_{\max}$ | §26 |
| $B_{\rm XS}$（橫截面 bootstrap） | 500 | §29 |
| $\tau_{\mathsf D}$（mechanism grouping 切點） | 0.90（＝接受 19% 能量損失） | §23.1,§35 |
| raw redundancy 切點 | 0.90 | §36 |
| $\pi$ 門檻（三項） | 0.80 | §34 |
| $M$ 門檻（兩項） | 0.60 | §34 |
| breadth 警示 | $\tilde b<0.5$ / $\tilde b^{\rm p10}<0.2$ | §31.3 |
| ${\rm Cov}^{\rm feat}_{\mathsf D}$ 目標 | 0.90 | §39 |
| ${\rm Cov}^{\rm feat}_{\mathsf V},{\rm Cov}^{\rm feat}_{\mathsf L}$ 目標 | 0.80 | §39 |
| $m_{\mathsf D}$ 警示門檻 | $<4$ → grouping 不具鑑別力 | §23.4 |
| $\delta_{\rm drift}$ 警示門檻 | $>0.3$ | §10 |
| $T_{\rm eff}/T$ 警示門檻 | $<0.1$ | §3.1 |
| $\epsilon$（數值下限） | $10^{-12}$ | §0.3 |

**這些門檻是營運性的篩選參數，不是統計判準。** 統計顯著性一律由 permutation 決定（§14、§18、§46）。門檻的作用是把「顯著但不可用」的 feature 濾掉，其數值可以依研究目標調整，但**必須在研究開始前鎖定並記錄**——事後調整門檻直到結果好看，是本方法無法偵測的一種過擬合。

---

## §67. 統一的表徵

### 一個物件

整套方法最終歸結為：對任何 cross-sectional 訊號 $s$（已 within-date 標準化、$\|s\|_W=1$），

$$\boxed{c^{\rm full}_s=s^\top W\Phi\ \in\mathbb R^{d_\Phi},\qquad c^{\rm inh}_s=c^{\rm full}_sP_V,\qquad c^{\rm nov}_s=c^{\rm full}_s(I-P_V).}$$

從這一個物件，可以讀出訊號的全部幾何性質：

$$\begin{aligned}
&Q_s=\|c^{\rm full}_s\|^2 && \text{強度}\\
&\big(Q^{\mathsf D}_s,Q^{\mathsf V}_s,Q^{\mathsf L}_s\big) && \text{機制分解}\\
&\tilde f_s && \text{機制指紋（維度校準後）}\\
&A^{\rm tgt}_s=\big(c^{{\rm full},\mathsf D}_s\big)^\top r_{\mathsf D} && \text{方向}\\
&\eta^{\mathsf D}_s && \text{目標精準度}\\
&\tilde\nu_s && \text{相對於 feature library 的新穎性}\\
&\pi_s,\ M_s && \text{時間與橫截面穩定度}\\
&\tilde b_s,\ \tilde R_s && \text{橫截面支撐}
\end{aligned}$$

而兩個訊號之間的關係為

$$A^{\rm mech}_{s_1,s_2}=\hat c_{s_1}^\top\hat c_{s_2},\qquad S^{\rm mech}_{s_1,s_2}=\big|A^{\rm mech}_{s_1,s_2}\big|.$$

### 三個應用

**Feature 層級**：

$$X_p\ \xrightarrow{\ \Phi(Y)\ }\ c_p\quad\Longrightarrow\quad\text{篩選出真實、穩定、機制互補的 feature}.$$

**Alpha 層級**：

$$\hat\alpha\ \xrightarrow{\ \Phi(Y)\ }\ c_\alpha\quad\Longrightarrow\quad\text{模型預測的是方向、活躍度、還是大行情方向}.$$

**交叉層級**：

$$c_p\ \longleftrightarrow\ c^{\rm inh}_\alpha\quad\Longrightarrow\quad\text{模型繼承了哪些已知機制},$$

$$\varepsilon_\alpha=\tilde\alpha-\text{OOSLinear}(X)\ \xrightarrow{\ \Phi(Y)\ }\ c_\varepsilon\quad\Longrightarrow\quad\text{模型創造了什麼新東西}.$$

### 這個統一的代價

把三個問題放進同一個座標系，換來了可比較性，但也意味著**三個問題共用同一組假設**：

1. 三個 lifting 張成的空間足以描述我們關心的預測型態（§8）；
2. 訊號到 label 的映射在 feature 側是線性的，或可被 transform family 覆蓋（§4.3）；
3. $\Phi$ 在研究期與 OOS 期之間足夠穩定（以 $\delta_{\rm drift}$ 監控，§10）。

三者任何一個失效，三個問題的答案會**同時**失效，且失效方式相同。這是統一座標系的結構性風險，也是為什麼 $\delta_{\rm drift}$、$K_{\rm eff}$ 的穩定度、$T_{\rm eff}$ 這三個座標系健康度指標必須與每一份報表同行。

換來的是一個在量化研究中罕見的能力：**用同一把尺，量 feature、量模型、量模型減去 feature 之後剩下的東西。**

---

# 附錄 A：方法邊界

本附錄集中陳述本方法**不能**回答的問題。這些不是缺陷，而是適用範圍；把它們寫清楚，是為了避免把「本方法看不到」誤讀為「不存在」。

### A.1 Label 側的函數族是有限的

$\Phi$ 由三個 lifting 張成：$z$、$\operatorname{svec}(zz^\top)$、$\|z\|z$。未來報酬的其他函數不在其中，例如：

- 極端尾部的符號（$\mathbf 1[z_h>3]$ 這類指示函數）；
- 三階以上的結構（偏態、峰態的橫截面差異）；
- 跨日的路徑形狀（band 內只有同一時點起算的若干 horizon）。

一個只預測這些的訊號，本方法會給出 $Q\approx0$。

擴充方式是增加 lifting（例如加入 $\operatorname{svec}$ 的三階版本或分位數指示函數），代價是 $d_\Phi$ 增大，進而稀釋 fingerprint 的分母、降低每個方向上的檢定力。這是一個明確的取捨：**維度愈多，能看到的型態愈多，但每一種都看得愈不清楚。**

### A.2 訊號側的映射是線性的

$c_s=s^\top W\Phi$ 對 $s$ 線性。若 $Y=X^2+\text{noise}$ 且 $X$ 對稱分布，則 $\mathbb E[X\cdot z]=\mathbb E[X^3]=0$，raw $X$ 的 $Q\approx0$。

§4.3 的 transform family 把邊界從「$X$ 本身是否線性相關」推進到「$X$ 的指定變換族中是否有任何一個線性相關」。但它仍然不是「$X$ 的任何函數」。

因此本方法的正確陳述是：

> **它回答的是「這個 feature 的這些表徵是否具有預測資訊」，而不是「原始變量 $X$ 在所有可能的非線性變換下都沒有資訊」。**

這個區別在最終模型本身是非線性模型時特別重要：模型可能使用了一個 transform family 未涵蓋的變換，於是 §52 的機制繼承分析會低估繼承程度，而 §58 的殘差分析會把它歸為「新機制」。§60 Case E 的第 4 步驗證（嘗試構造捕捉該方向的新 feature）正是為了區分這兩種情形。

### A.3 Circular shift 虛無對持續性訊號的檢定力

§13.4 已完整處理：以 $\kappa_p$ 診斷、以分層 within-date permutation 作為替代虛無。但必須意識到兩個虛無檢驗的是**不同的假設**：

- N1：「feature 的橫截面排序與正確時間的報酬機制無關」；
- N2：「在控制分層（產業×市值）之後，個股層級的 feature 排序與報酬機制無關」。

N2 較弱：它無法偵測純粹的產業層級或市值層級效應（那些效應在分層內被保留）。因此對慢變 feature，本方法能給出的最強結論是「存在個股層級的預測力」，而非「存在預測力」。

### A.4 等權是列數等權，不是影響力等權

§3 已完整處理：以 $T_{\rm eff}$ 診斷。實務後果是 $Q$ 不對市場 regime 免疫，$\mathsf V$ 與 $\mathsf L$ block 尤甚（貢獻量級 $\sim\|z\|^2$）。

### A.5 低維機制空間的鑑別力

§23.3 已完整處理。當 $m_{\mathsf D}\le3$ 時，「這兩個 feature 是不是同一個機制」這個問題用夾角答不出來。本方法的設計因此把選擇建立在覆蓋準則（§38.2）而非配對門檻上，但「模型繼承了 feature $p$」這類**敘述性**結論仍然受限於此，必須搭配 §52 的虛無分位數與 §57 的線性歸因才能成立。

### A.6 關聯不是因果

$Q_p$ 大只代表 feature 與未來報酬機制有穩定的橫截面關聯。它可能來自真實的定價機制，也可能來自共同暴險於某個未被建模的風險因子。本方法不區分這兩者。

若要區分，必須在 §5 之後加入一層「對已知因子暴險做橫截面迴歸並取殘差」的前處理。這不改變後續任何公式，但改變全部結論的語意：從「這個 feature 預測報酬」變成「這個 feature 預測已知因子之外的報酬」。

### A.7 $Q$ 不是報酬

$Q$ 沒有考慮交易成本、容量、換手率、借券可得性、衝擊成本。一個 $Q$ 很高但 $\tilde b$ 極低、$\tilde R$ 接近 1 的 feature，在真實交易中可能完全不可用。§30–§32 的診斷是為此設計的，但它們只是警示，不能取代實際的成本模型與容量分析。

---

# 附錄 B：判讀速查表

| 看到 | 意義 | 該做什麼 |
|---|---|---|
| $K_{\rm eff}<K$ | whitening 遇到病態；某些 horizon 方向不穩定 | 回到 §6 改用增量基底或合併 horizon；截斷是最後手段 |
| $T_{\rm eff}/T<0.1$ | 結論由少數事件日決定 | 該 block 的全部結論降級為探索性 |
| $\delta_{\rm drift}>0.3$ | 座標系本身在漂移 | 跨期比較失效；先排除座標問題再談模型衰退 |
| $r_{\rm sig}=1$ | library 只照亮一個機制 | 補充新資料來源，而非增加同類 feature |
| $r_{\rm eff}\gg r_{\rm sig}$ | 能量分散在微弱方向 | $V_r$ 不穩，後續結論降級 |
| $m_{\mathsf D}<4$ | 機制空間太小，夾角無鑑別力 | 分群只作報表組織；選擇交由覆蓋準則 |
| $\zeta(\tau)\le0$ | 觀測相似度結構與虛無無異 | 同上；不得以 cluster 為單位下結論 |
| 所有 feature 都不顯著 | 常見原因是 $B_{\rm null}$ 太小 | 檢查 $1/(B_{\rm null}+1)$ 是否夠小（§18.2） |
| $\tilde f_p\approx(1/3,1/3,1/3)$ | 該 feature 是雜訊 | 不論 raw $f$ 看起來多有結構 |
| $\tilde b_p\approx1,\tilde R_p\approx1$ | 集中度與淨效果都與雜訊無異 | 高度懷疑 $Q_p$ 來自少數極端名次 |
| $\pi^{\rm time}$ 高但 $M^{\rm time}$ 低 | 永遠很強，但每期預測不同東西 | 硬性淘汰；回到 feature 建構層拆解 |
| $\pi^{\rm time,sign}<0.8$ | 方向會翻轉 | 淘汰；不得在 OOS 改符號 |
| $\eta^{\mathsf D}_p\approx0$ 但 $Q^{\mathsf D}_p$ 大 | 方向性資訊指向別的 horizon 組合 | 標記 `TARGET_MISMATCH`，轉入其他 band |
| ${\rm Cov}^{\rm span}_{\mathsf D}=1$ 而只選了 3 個 feature | span coverage 已飽和 | 這是 $d_{\mathsf D}=K$ 的必然；用 ${\rm Cov}^{\rm feat}$ |
| $\tilde f_\alpha$ 以 $\mathsf V$ 為主 | 模型在預測波動度而非方向 | 多空排序會變成做多 gamma；改訓練目標 |
| $\nu_\alpha=0.4$ | 低於虛無基準 $(d_\Phi-r_{\rm sig})/d_\Phi$ | 用 $\tilde\nu$ 判讀；0.4 其實代表「不新穎」 |
| $b^{\rm sig}$ 高但 $b^{\rm pred}$ 低 | 廣泛下注但預測力集中 | 對低貢獻名次施加門檻，可降成本而幾乎不損報酬 |
| $R^2_{\rm lin}$ 低但 $Q_\varepsilon$ 不顯著 | 非線性是純粹的成本（Case C） | 加強正則化或退回 surrogate |
| $R^2_{\rm lin,library}\gg R^2_{\rm lin,selected}$ | 模型用了未被選入的資訊 | 檢查模型輸入是否與 selected 一致 |
| $\Delta R^2_g\approx0$ 但 $R^2_{g,\rm own}$ 高 | 與其他 family 互相遮蔽 | 兩者擇一即可，但缺一不可 |
| $A^{\mathsf D}_{p,\alpha}=0.92$ 且前五名相近 | 低維空間的必然，不是繼承 | 需要虛無 max 分位數與 §57 共同支持 |
| $\tilde\nu_\varepsilon>0.7$ 且 $p^{\rm nov}_\varepsilon<0.05$ | 可能發現新機制（Case E） | 四步驗證；最終目標是造出捕捉它的新 feature |

---

*本書的全部定義、公式、預設值與報表欄位構成一個自洽的完整規格。§65 的恆等式測試與校準測試是驗收標準：五項核心校準測試（T1–T5）全部通過之前，任何顯著性結論都不應被採信。*
