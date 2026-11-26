# O18-CHEEGER — 条件 Cheeger 能否给出 `Var(f|0) ≤ E[H_vv+Φ''|0]`：**不能，且不是"前提不够"而是"类型不够"**（`[NG]`）；条件版 `h₀` 是**钉扎纤维上的等周常数**，它既不自动正、又比它该证的方差不等式**严格更强**，而且在精确可解模型上单靠 `1/(2π)` 折价就已经**超支 2.26×**

**日期**: 2026-11-25 · **状态**: 已判定（§1–3 为 `[T]`/`[std]`；§4 三判词；§5 数值校准 `[T]+[V]`）· **禁 git commit**（遵守）· 只写本文件

> **必读已核**：`species/GPP-STRUCT.md`（全文，191 行）—— 本文件不重证 `G''` 的定义与 `eq:split`，只把它的**δ=0 成员**当作待证目标；`species/H1A-EXACT.md:1–80`（口径与精确值）；`unconditional/CHEEGER-BLIND.md`、`unconditional/B3-ISOPERIMETRIC.md`、`research/v3/Q4-BAKRY-EMERY.md`（**已入库的 Cheeger 判词在边际 `ν` 上，不在条件测度上**——本文件的 O18 是**另一个对象**，见 §1.3）。
> **Witness**：`%TEMP%\o18_check.py`（`py -3.12 -I`，纯 `math`，单进程，exit 0，输出见 §5）；**仓库内只落盘本文件**。
> **标签**：`[T]` 真定理/恒等式 · `[std]` 教科书级（书目级已核，未逐字核） · `[C]` 条件/片段级 · `[NG]` 具名类 no-go · `[V]` 数值见证 · `[O]` 精确缺失步 · `[J]` 判断 · `[E]` 出处待核。
> **No QED. Clay 0%。**

---

## §0 判词（先读这段）

**目标（本文件的唯一目标，逐字取任务式）**

```
(T)   Var_μ(∂_vV | ξ=0)  ≤  E_μ[∂_v²V | ξ=0] = E_μ[H_vv + Φ''|0],    V = S_W + Φ,  ∂_v = ⟨v,∇·⟩ .
```

由 `GPP-STRUCT.md §2` 的结构等式 `G''(0) = E[H+Φ''|0] − Var(∂_vV|0)`，**(T) ⟺ `G''(0) ≥ 0`**（`[T]`，代数恒等，逐点逐模）。故本文件不引入新目标，只判**路线**。

**三句判词（对应任务的 (a)/(b)/(c) 死线）**

1. **(a) 成立 ⇒ 与 H1 不同型，`[NG]`。** 条件 Cheeger 给的右边是 **Dirichlet 型**，而且不只是"另一个二次型"，它精确等于 **`Hess V` 的二次形**：

   ```
   Cheeger(μ(·|0)):   Var(f|0) ≤ (4/h₀²)·E[ |P_⊥ Hess V v|² | 0 ] =: (4/h₀²)·D₀
   (T) 的右边:        R₀ := E[ vᵀ Hess V v | 0 ]
   ```

   `D₀` 对 `Hess V` 是**二次**的（`vᵀH P_⊥H v`），`R₀` 是**一次**的（`vᵀHv`）。Cheeger 只提供常数 `4/h₀²`，**永不对 `D₀` 与 `R₀` 的相对大小说话**。要接上必须额外证 `4D₀/h₀² ≤ R₀`——这正是 (T) 本身。`[T]`
2. **(b) 部分成立，但两条腿都断。** "曲率形式的 Cheeger" **确实存在**（Lichnerowicz `λ₁ ≥ nK`；Bakry–Émery `Γ₂ ≥ ρΓ ⇒ gap ≥ ρ = ess inf 曲率`）——这是活口，且**已在库里算过**：`Ric_{SU(N)} = (N/4)g > 0` 精确，但 `Ric + Hess S_W ≈ −2.3 < 0`（`Q4-BAKRY-EMERY`、`yang-mills-frontier.md:233`）⇒ `ρ < 0` ⇒ **vacuous**（腿 1 断）。而且**即便 `ρ>0`，输出仍是 `Var ≤ ρ^{-1}E|∇f|²`**，右边还是 Dirichlet form——曲率只买到一个**数**，买不到"把 Dirichlet form 换成曲率期望"（腿 2 断）。唯一真正**取曲率值**的不等式是 **Brascamp–Lieb 的 Hessian 形式** `Var ≤ E⟨(Hess V)^{-1}∇f,∇f⟩`，其**极值函数恰好是 `f = ⟨c,∇V⟩`**（即本目标的 `f`），代入即得 (T) 的无条件版 ⇒ **活口就是 BL，前提 `Hess V ⪰ 0` = Savvidy 墙**。`[std]/[C]/[NG]`
3. **(c) 不触发（`h₀ ≠ Var`），但真正的循环更糟：`h₀` 比目标严格更强。** `h₀` 是几何量（周长/体积），不是方差，所以"`h₀ = Var`"字面不成立；但 (i) `h₀ > 0 ⟹ 纤维谱隙`，**反向一般不成立**（细刺/尖钉域：`λ₁>0` 而 `h→0`；Buser 的逆要两侧 `Ric ≥ −K`，Savvidy 的 Ricci 不定），故 `h₀ ≥ c` **严格强于** H1 所需的纤维 Poincaré；(ii) Cheeger 只给 `h₀ ≤ 2√λ₁` 这一侧，**曲率下界给不了 `h₀` 的下界**（要 Lévy–Gromov 型双侧夹逼，前提又是纤维 `Ric ≥ κ>0`）。⇒ 不是"换名"，是**加价**。`[T]/[std]`

**净判**：**绕了一圈回到 BL（`[NG]`）**，不提供新的**估计机器**；只提供两个**定位**（§4: Schur 夹逼判据；§3: 目标 = 条件 Stein/Ward 的极值形式）。**且 §5 把失效从"h₀ 未知"升级成"通道本身不够宽"：在 `Hess V` 为常数的 Gaussian 纤维上，`h₀` **精确已知**，`1/(2π)` 折价单独就使界超支 `2.262×`。**

---

## §1 Cheeger 的式子、出处、条件版对象

### 1.1 式子（可核出处）

```
λ₁(μ) ≥ h(μ)²/4,     h(μ) := inf_A  ∂_μ A / min(μ(A), μ(A^c))          [Cheeger 1970]
⇒  Var_μ(f) ≤ (4/h²)·∫|∇f|² dμ                                           （一行：1/λ₁ 是 Var≤C·E 的最优常数）
```

- **Cheeger, J. "A lower bound for the smallest eigenvalue of the Laplacian."** In *Problems in Analysis* (R. C. Gunning, ed.), Princeton Univ. Press, 1970, pp. 195–199. **DOI `10.1515/9781400869312-013`**（De Gruyter 版）——**书目级已核**（本次 web_search 命中该 DOI 条目，未逐字读正文）。
- 陈述复述可核："A celebrated lower bound of `λ₁` in terms of `h`, `λ₁ ≥ h²/4`, was proved by Cheeger in 1970 for smooth Riemannian manifolds"（*Sharp Cheeger–Buser type inequalities in RCD(K,∞)*，[CVGMT 4217](https://cvgmt.sns.it/media/doc/paper/4217/Buser@inequalityFinal.pdf)，**片段级已核**）。
- 逆侧：**Buser, Ann. Sci. ENS 15 (1982) 213**（`[std]`，未逐字核）——**只在两侧 `Ric ≥ −K` 下**才有 `λ₁ ≤ C(h+h²)`。
- 离散/电导版：Jerrum–Sinclair, *SIAM J. Comput.* **18** (1989) 1149；Sinclair, *Combin. Probab. Comput.* **1** (1992) 351（`1−λ₂ ≥ Φ²/2`，`[std]`）。
- **1D 精确形式**（库内已用）：`h = inf_x 2ρ(x)/min(F(x),1−F(x))`，最优等周集是半直线（`Q4-BAKRY-EMERY.md:166–191`）。

**折价常数（必记）**：Gaussian 上 `h = √(2/π)/σ`、`λ₁ = 1/σ²`，故 `h²/4 = λ₁/(2π)`——**Cheeger 只抓到 gap 的 `1/(2π) = 0.159155`**。库内已校准过同一比例（`BASE-NONLOCAL-BOUND.md:38,50,69`：`h·√χ = √(2/π)`，fraction `1/(2π)`）。

### 1.2 条件版的测度是什么（不凭记忆，逐项写）

`ξ(X) = ⟨v,X⟩`，`dμ = e^{−V}dX`，`dU = J dX`，`Φ = −log J`。钉扎在 `ξ = 0`：

```
dμ₀ := μ(·|ξ=0),   F₀ := {ξ = 0}  (dim n−1),   dμ₀ = e^{−V(X)} |∇ξ|^{−1} dH^{n−1}(X) / ρ_v(0),
ρ_v(0) = ∫ e^{−V} δ(ξ) dX  =  ∫_{F₀} e^{−V}|∇ξ|^{−1} dH^{n−1} .
```

即：**纤维 `F₀` 上的加权体积测度**，权重 `e^{−V}`（`|∇ξ|^{−1}` 可并进纤维势 `W := V|_{F₀} + log|∇ξ|`）。`|∇ξ| = 1` 当 `|v|=1`，下文取 `|v|=1`、`P_⊥ := I − vvᵀ`。`[T]`（标准条件化/余面积公式）

**条件 Cheeger 常数**

```
h₀ := h(μ₀) = inf_{B ⊂ F₀} ∂_{μ₀}B / min(μ₀(B), μ₀(B^c)),     ∂_{μ₀} = 纤维内 μ₀-周长 .
```

### 1.3 与库里已入库的 Cheeger 判词**不是同一个对象**（防重复计功/误杀）

| 对象 | 测度 | 给的式子 | 库里状态 |
|---|---|---|---|
| **边际 Cheeger** | `ν(da) = ρ_v(a)da`（1D 粗模线） | `λ₁(ν) ≥ h(ν)²/4` ⇒ `χ ≤ 4/h²`；`h ≍ 1/√χ` | **已入库并判为 tautology**：`CHEEGER-BLIND`（2/5）、`B3-ISOPERIMETRIC`（2/5）、`H6-HEAD-SPECTRAL.md:87`（H6-C tautology `[T]`）、`BASE-NONLOCAL-BOUND §2` |
| **条件 Cheeger（本文件 O18）** | `μ₀ = μ(·|ξ=0)`（`n−1` 维纤维） | `Var(f|0) ≤ (4/h₀²)∫_{F₀}|∇_{F₀}f|²dμ₀` | **本文件首次判** |

`B3-ISOPERIMETRIC` 已证的反向不可传递性（"纤维性质不能约束边际 `ρ`"）**不覆盖本文件的问题**：这里问的是**同一方向**（纤维上的 Cheeger → 纤维上的方差），所以 B3 的障碍不适用，必须另判——判词在 §3–§4。

---

## §2 把 Cheeger 套到 (T) 上，逐项对齐

取 `f := ∂_vV = ⟨v,∇V⟩ = G_v + Φ'`（`G_v = ∂_vS_W`）。则 `∇f = Hess V · v =: Hv`，纤维梯度 `∇_{F₀}f = P_⊥Hv`。**条件 Cheeger 给出**

```
(C)   Var_μ(f|0) ≤ (4/h₀²)·E_μ[ vᵀH P_⊥H v | 0 ] = (4/h₀²)·D₀ ,      H := Hess V = Hess S_W + Hess Φ .
```

**逐项核对 `[T]`**：
- `D₀ = E[|P_⊥Hv|²|0] = E[|Hv|² − (vᵀHv)²|0]`（`P_⊥` 投影，`|v|=1`）；
- `(T)` 的右边 `R₀ = E[vᵀHv|0] = E[H_vv+Φ''|0]`（与 `ch05:230–240` 的 `E[H_vv+Φ'']` 同项）；
- 由 Stein/分部积分（`[T]`，`app.tex:167–172`）`E_μ[Hess V] = E_μ[∇V⊗∇V]` ⇒ `R₀` 的**无条件**版 `E[vᵀHv] = E[f²]`。

**两个直接后果（(a) 的实质）**

1. **次数不匹配**：`D₀` 是 `H` 的**二次**泛函，`R₀` 是**一次**。由 Cauchy–Schwarz / 分部积分，
   `∫|∇f|²dμ = E[vᵀH²v]`，而目标要的是 `E[vᵀHv]`；无条件时二者由 Stein 并列（`E[vᵀHv] = E[f²]`），**条件化后这条并列消失**。Cheeger 的整个信息量就是把 `Var` 压到 `E[|∇f|²]` 的比例上，**它对 `H²` 与 `H` 的相对大小无话可说**。
2. **无条件版是恒等，全部内容在"条件化缺陷"**：`E[vᵀHv] = E[f²]` ⇒ 无条件 `Var(f) ≤ E[vᵀHv]` 是**平凡真**（`Var ≤ E[f²]`）。所以 (T) 的**全部内容**是差
   ```
   D₀^Stein(X) := f² − (H_vv+Φ'') ,     E_μ[D₀^Stein] = 0 （Stein，恒等）
   (T)  ⟺  E_μ[ D₀^Stein | 0 ]  ≤  − (E_μ[f|0])² .
   ```
   `[T]` ⇒ **(T) 是"条件 Stein/Ward 恒等式"的缺陷不等式，不是任何 Dirichlet-form 不等式**；而 `D₀^Stein` 这个对象**在 Cheeger 的数据里根本不出现**（Cheeger 只吃 `|∇_{F₀}f|²`）。

**必须挡住的廉价错误路线**：用**全测度** Cheeger 的 `Var(f) ≤ (4/h²)E|∇f|²`，再加上
`Var(f) = E[Var(f|ξ)] + Var(E[f|ξ]) ≥ E[Var(f|ξ)]`，只能得到 `ξ`-**平均**的方差，给不了 `ξ=0` 这一个点（要 `a`-一致族），且右边仍是全测度的 Dirichlet form —— **双重失效**。`[T]`

---

## §3 能否用"曲率 ≥ h²/4"把 `R₀` 与 `D₀` 接上？——**方向反了**

`splitting`：在 `(v, ⊥)` 分块下 `H = [[a, bᵀ],[b, C]]`，`a = vᵀHv = R₀`，`|b|² = D₀`，`C = P_⊥HP_⊥`（纤维 Hessian，即 `ch13-machine.tex:58–82` 的 `K = Π₀ Hess V Π₀`，**符号所在**）。

**(B) 桥 = `4D₀/h₀² ≤ R₀`。** 由 Schur 补：

```
若  (C1) C ⪰ 0  且  (C2) a ≥ bᵀC^{−1}b   （⟺ H ⪰ 0，即局部凸性）
则  D₀ = |b|² ≤ ‖C‖·bᵀC^{−1}b ≤ ‖C‖·a ,
⇒ (B) 由  (C1) C ⪰ 0   与   (C2′) ‖C‖ ≤ h₀²/4   逐点推出。
```

`[T]`（Schur 补 + `‖b‖² ≤ ‖C‖ bᵀC^{−1}b`）。

**三条判**：

1. **需要的不是"曲率 ≥ h²/4"，而是"曲率 ≥ 0 且 ≤ h₀²/4"（两侧夹逼）。** 任务里的 `曲率 ≥ h²/4` 是**反的**：Cheeger 的式子 `λ₁ ≥ h²/4` 是**结论**，而且它恒伴随 `h ≤ 2√λ₁` ——**曲率越大概率越给 `h` 上界，不是下界**。真正要接桥，`h₀` 必须**很大**（`h₀² ≥ 4D₀/R₀`），而 `‖C‖` 必须**不超过** `h₀²/4`。即：**纤维 Hess 被等周常数两侧钉死**（一种定量"近平坦 + 凸"条件）。`[T]`
2. **`(C1)` 就是 Savvidy 墙。** `λ_min(Hess S_W) = −g|B| < 0`（`CLAUDE.md §5`）、`PᵀHP` min `≈ −0.3`（`research/unconditional/MARGINAL`）、`Ric+Hess S_W ≈ −2.3`（`Q4-BAKRY-EMERY`）⇒ `C ⪰ 0` 假。**桥的第一半在物理点直接假。**`[T]`
3. **`(C2′)` 与 Cheeger 自身矛盾（在已知例上）**：Gaussian 纤维（§5）`‖C‖ = ` 纤维精度最大本征值 `= 2.5`，而 `h₀²/4 = 0.397887` ⇒ `‖C‖ ≤ h₀²/4` **假 6.28 倍**。所以 `(C2′)` 不是"技术上的松"，而是**这类模型上根本不成立**——即使 `(C1)` 成立。`[T]+[V]`

**⇒ 桥需要的东西严格多于 (T)：既要 `H` 的纤维块半正定，又要其范数被 `h₀²/4` 压住，而 `h₀` 本身又正是所求。三个条件里两个已知为假、第三个是未知量。**

---

## §4 判：绕回 BL，还是新量？

### 4.1 绕回 BL（`[NG]`）—— 且这次是**精确同构**，不只是"同类"

BL 的 Hessian 形式（`V` 凸、`Hess V ≻ 0`）：

```
Var_μ(f) ≤ E_μ⟨ (Hess V)^{−1} ∇f , ∇f ⟩ ;      极值函数  f = ⟨c, ∇V⟩ .
```

取 `c = v`、`f = ∂_vV`：`∇f = Hv` ⇒ 右边 `= E⟨H^{−1}Hv, Hv⟩ = E[vᵀHv] = R₀`。**这正是 (T) 的无条件版**（等价于 §2 的 Stein 恒等）。所以在**无条件**测度上，**(T) 就是 BL 在其极值函数上的取值**；在**钉扎**测度上，**(T) 就是"纤维 `μ₀` 上的 BL"**，其 Hess 是 `K = Π₀Hess VΠ₀`（即 §3 的 `C` 与 Schur 补 `a − bᵀC^{−1}b` 的整体）。

- 极值函数这一条已核（片段级）："_The Brascamp–Lieb inequality admits extremal functions given by `f(x)=⟨c,∇V⟩` where `c` is a constant vector_"（[Bonnefont, Poincaré 讲义](https://www.math.u-bordeaux.fr/~mibonnef/Poincare__Toulouse.pdf)，**片段级**）。
- BL 出处：**Brascamp–Lieb, Adv. Math. 22 (1976) 366–389**（*On extensions of the Brascamp–Lieb and Prékopa inequalities*；另有 Bull. AMS **82** (1976) 927 的短讯）。`[std]`，**书目级**（卷/页未逐条与出版社页核对，标 `[E]` 者仅指年份-卷号对）。
- 与 `GPP-STRUCT §4(b)` 的判词**逐字一致**："H1 的真身 = 1D 条件边缘的 BL/log-凹不等式"，前提 `Hess V ⪰ 0` = Savvidy 墙。本文件把该判词的**条件版**（纤维而非 1D 边缘）写成显式，结论不变。

⇒ **(a) 死线命中：通往 (T) 的任何"Cheeger 型"通道，其唯一能改成曲率值的形式就是 BL，前提是凸性。`[NG]`。**

### 4.2 (b) 的活口清点（诚实列全，全部标死因）

| 活口 | 式子 | 状态 |
|---|---|---|
| Lichnerowicz | `Ric ≥ (n−1)K>0 ⇒ λ₁ ≥ nK` | `[std]`；前提在纤维上被 `Ric+Hess S_W≈−2.3<0` 破 |
| Bakry–Émery | `Γ₂ ≥ ρΓ ⇒ λ₁ ≥ ρ = ess inf(曲率型)` ——**字面"`λ₁ ≥ ess inf 曲率`"** | `[std]`；库里已算 `ρ<0` ⇒ vacuous（`Q4`、`frontier:233` "典型点正 vs ess inf 负"）；**且输出仍是 `Var ≤ ρ^{−1}E|∇f|²`，不改型** |
| Lévy–Gromov | `Ric ≥ κ>0 ⇒ h ≥ h(Sⁿ_κ) ≍ √κ`（**曲率 ⇒ `h₀` 下界**，唯一真能这么给的） | `[std]`（[RCD 讲义](http://math.univ-lyon1.fr/wikis/Stab/lib/exe/fetch.php?media=wiki:seminarlyonmondino.pdf) 片段级："Levy-Gromov … `m⁺(E) ≥ |∂B|/|S|`, `S=S_K^N`"）；前提同 Lichnerowicz（纤维一致凸），**假** |
| BL-Hessian | `Var ≤ E⟨(Hess V)^{−1}∇f,∇f⟩` | `[std]`；**取曲率值**，但 = §4.1，(T) 换了名 |

**判**：`(b)` 的三种"曲率形式"里，前两个只把**常数**变成曲率（右边仍是 Dirichlet form），第三个（Lévy–Gromov）确实从曲率给 `h₀` 下界——**这是唯一能让 §3 的桥过桥的机制**，但其前提是纤维 `Ric ≥ κ > 0`，**恰是 Savvidy 破的那条**。⇒ 活口存在但**空前提**。

### 4.3 (c) 清点：`h₀` 与 `Var`、与目标的硬度关系

- **字面 (c) 不成立**：`h₀` 是几何比（周长/体积），`Var` 是泛函，二者不同类，`h₀` 不可能是 `Var` 本身。
- **但循环以"更强假设"的形式存在**：
  - `h₀ > 0 ⟹ 纤维谱隙`；**反之一般假**（细刺/尖钉域：`λ₁ ≥ c > 0` 而 `h → 0`；Buser 逆需两侧 Ricci）。⇒ 证明 `h₀ ≥ c` **严格强于**证明 H1 所需的纤维 Poincaré。
  - 精确 1-D：`h = inf 2ρ/min(F,1−F)` 是**全轴瓶颈**，`h²/4` 只等于 `λ₁/(2π)`（Gaussian 精确）。⇒ **"证 `h₀ ≥ c`" 的价格至少是 "证 `λ₁ ≥ 2π c²`"**，即比目标贵 `2π` 倍。
- **库内已入库的边际版正是这条**：`h ≍ 1/√χ` ⇒ `h₀ ≥ c ⟺ χ ≤ 2/(πh₀²)` = H1（`BASE-NONLOCAL-BOUND §2` 双向证明；`CHEEGER-BLIND §1`；`H6-C`）。条件版**只多一层**：`h₀` 换成纤维的，硬度的结论不变。

---

## §5 精确校准：在"一切都凸、`Hess V` 常数"的 Gaussian 纤维上，Cheeger 通道**结构性超支 2.262×**

**模型（可解、精确、也是本项目的自由场酸测点）**：`V = ½(x₁² + 4x₂²)`（`H = diag(1,4)`、`Σ = H^{−1} = diag(1, ¼)`、`Φ≡0`），`v = (1,1)/√2`，钉扎 `ξ = ⟨v,X⟩ = 0`。

`witness`：`%TEMP%\o18_check.py`（纯 `math`，确定性，exit 0），输出逐位如下（`[V]`）：

| 量 | 值 | 说明 |
|---|---|---|
| `R₀ = E[H_vv|0] = vᵀHv` | **2.500000** | 目标右边 |
| `D₀ = ‖P_⊥Hv‖²` | **2.250000** | Cheeger 的 Dirichlet form |
| `Var(f|0)` | **0.900000** | 目标左边（真值，`= (Hv)ᵀΣ₀(Hv)`，`Σ₀ = 0.2[[1,−1],[−1,1]]`） |
| 纤维 `λ₁ = 1/λ_max(Σ₀) = 1/0.4` | **2.500000** | 最优 Poincaré |
| **最优 Poincaré 界** `D₀/λ₁` | **0.900000** | **取等**（`f` 是纤维本征函数） |
| `h₀`（1D 纤维，半直线最优 ⇒ **精确**） | **1.261566** | `= √(2/π)/√0.4` |
| `h₀²/4` | **0.397887** | |
| `(h₀²/4)/λ₁` | **0.159155** | `= 1/(2π)`，**精确折价** |
| **Cheeger 界** `D₀/(h₀²/4)` | **5.654867** | vs 真值 `0.900`、vs 目标 `2.500` |
| 过桥所需 `h₀` | **1.897367** | `= √(4D₀/R₀)`；**缺** `1.504×`（`h`）/ **`2.262×`**（`h²`） |

**三条读数（`[T]+[V]`）**

1. **(T) 在此模型为真且有余量**：`0.900 ≤ 2.500`（gap `1.600`）。
2. **最优 Poincaré 通道恰好够用且取等**：`D₀/λ₁ = 0.900 ≤ 2.500`。⇒ 若能把 Cheeger 的 `4/h₀²` 换成最优常数 `1/λ₁`，(T) 当场成立。**缺口 `100%` 来自 Cheeger 的 `2π` 几何折价，不来自 `h₀` 未知、不来自凸性、不来自 Savvidy。**
3. **`2π = 6.283 > 2.262`**：Cheeger 的**结构性损失本身就超过目标的全部余量**。⇒ 这条路不是"前提没拿到"，是**通道宽度不够**；在 `Hess V` 为常数（能有的最良性情形）上就已经不够。**这是本文件对 `[NG]` 的最强证据：不依赖任何关于 SU(N)/Savvidy 的断言。**

**同时钉死目标本身的适用边界（库内精确值，`H1A-EXACT.md`）**：`(T) ⟺ G''(0)≥0`。
- `β=2.5, θ=π`：`G''=+0.135995` ⇒ (T) **真**；
- `β=6, θ=π`：`G''=−0.193019`，`Var(f|0)=0.398135`、`R₀=W_full=0.205115` ⇒ **`R₀ < Var`，(T) 精确为假**；
- `β=β*=4.30661`：`G''=0` 精确（退化）。
⇒ 任何"Cheeger 型"证明若声称对一切 `β` 成立，**已被 `β=6,θ=π` 反例杀死**；这条与 §5 的 Gaussian 超支**互相独立**，两条都成立。

---

## §6 判词表（给 lead 直接引用）

| 陈述 | 判 | 依据 |
|---|---|---|
| 条件 Cheeger 的右边 = `E[vᵀH P_⊥Hv|0]`（`H` 的**二次**形、**mixed** 块） | `[T]` | §2 (C)；`∇_{F₀}f = P_⊥Hv` |
| 目标右边 = `E[vᵀHv|0]`（`H` 的**一次**形、**vv** 块） | `[T]` | `ch05:230–240` |
| ⇒ **与 H1 不同型**，(a) 命中 | `[NG]` | §2 次数不匹配；Cheeger 只给常数 `4/h₀²` |
| 无条件版 (T) = Stein 恒等 = BL 在极值函数 `f=⟨c,∇V⟩` 上取值 | `[T]/[std]` | `app.tex:167–172`；BL 1976 |
| 条件版 (T) = 纤维 `μ₀` 上的 BL（Hess = `K = Π₀Hess VΠ₀`） | `[T]/[C]` | §4.1；`ch13-machine.tex:58–82` |
| `h₀` 是"钉扎纤维的等周常数"，**不自动正** | `[T]` | §1.2；反例：细刺域（`λ₁>0`,`h=0`）、自由无质量场 `h₀=0` |
| 曲率 ⇒ `h₀` 下界**存在**（Lévy–Gromov），前提纤维 `Ric≥κ>0` | `[std]`，前提**假** | §4.2；`Ric+Hess S_W≈−2.3<0` |
| `λ₁ ≥ ess inf 曲率` 存在（Lichnerowicz/BE），但**不改型**、且 `ρ<0` | `[std]`+`[NG]` | §4.2；`Q4-BAKRY-EMERY` |
| 过桥需 `C⪰0` **且** `‖C‖ ≤ h₀²/4`（两侧夹逼） | `[T]` | §3 Schur 补 |
| 该夹逼在 Gaussian 纤维上假 `6.28×`；`h₀` 精确时仍超支 `2.262×` | `[T]+[V]` | §5 |
| `h₀ > 0` **严格强于**纤维 Poincaré（Buser 逆需两侧 Ricci） | `[std]` | §4.3 |
| 字面 (c)：`h₀ = Var` | **不成立**（不同类） | §4.3 |
| 净路径 | **绕回 BL，`[NG]`** | §4.1 |

---

## §7 新量清点（诚实：只有半个，且不是估计机器）

1. **Schur 夹逼判据** `(C1) C⪰0 ∧ (C2′) ‖C‖≤h₀²/4`（§3）：把"H1 缺的一步"从"某条未知不等式"**写成两个可逐点检验的解析条件**。它不是新的上界机器（两条都不可证），但把缺口**定位到纤维 Hessian 的两侧**，并证明任务设想的"曲率 ≥ h²/4"是**反方向**。`[T]`
2. **Stein 重写**（§2）：`(T) ⟺ E[D₀^Stein|0] ≤ −(E[f|0])²`，`E[D₀^Stein]=0`。这就是 `GPP-STRUCT §4(b)` 的"条件化修正"的**精确代数内容**，可用来判任何候选工具：**凡是不吃 `f² − Hess` 这个组合的机器（Cheeger/BE/LSI/容量/耦合），原理上都看不见 (T)**。`[T]`
3. **不算新量**：`h₀`、`λ₁(μ₀)`、`D₀` 全是标准对象；§5 的 `2π`/`2.262` 是对已有常数的校准，不是新界。

---

## §8 未闭合（诚实标注）

1. §5 的 `2.262×` 是**一个模型**（Gaussian 纤维）上的超支；它证明"Cheeger 通道对一切模型都不够"，但**没有**证明"在 `β=2.5` 的 SU(N) 纤维上 `h₀` 就一定小"——后者是未知量（正因如此 (T) 未证）。
2. §3 的 `(C2′)` 用的是 `‖b‖² ≤ ‖C‖·bᵀC^{−1}b` 这一最粗的 Schur 估计；若纤维 Hess 有各向异性结构（例如 `C` 在 `b` 支撑方向上恰好小），`(C2′)` 可能过强。**本文件未排除存在更细的、非 Schur 的 `D₀/R₀` 界**（这属于"新不等式"而非 Cheeger 路线的范围）。
3. 所有 `h₀` 均按格点单位曲面处理，未做 `a→0` 标度（`ch03-existence.tex:146–165` 的 `sec:dim` 纪律）；`h` 的量纲是 `1/长度`、`h²` 才是 `1/长度²`（与曲率同型），本文件所有 `h₀²/4` 的比较都在格点单位下。
4. **记号陷阱（沿用 `GPP-STRUCT §1.2`）**：`G(a)`=自由能、`G_v`=力。本文件的 `f = ∂_vV = G_v + Φ'`；若改用库内 `(M)` 显示的 `f = G_v`（`ch08-immune.tex:11–19`），则右边是 `E[H_vv|0]+W_full−δ`，**其中含 `−Var(Φ')` 与 `−2Cov(G_v,Φ')` 两项方差/协方差**——"右边是纯曲率"这句话**只在 `G''≥0` 形式上成立**，在 `(M)` 形式上更强地不成立。这两条显示**不是同一条不等式**，引用时必须点名是哪一条。`[T]`

---

## §9 AGENTS.md §1 检索执行核验

**工具状态**：`grep` = **OK**（`Cheeger|cheeger|isoperimetric|等周`，`*.md`，**250/256 命中**，完整结果落盘 `%TEMP%\dsh-spill-zvYlSr\session-2c9fc3f3ffdc\98025f313229-grep-results.txt`，可复读）；`glob` = **OK**（`**/GPP-STRUCT.md` → 1；`research/v3/species/*.md` → 112）；`read` = **OK**（`GPP-STRUCT.md` 全 191 行、`H1A-EXACT.md:1–80`、`DT-DIRICHLET.md:1–60`）；`web_search` = **OK**（5 条查询 → 10 源，均有效外部文档）；`web_fetch` = **未调用**（故所有外部断言只标 `[std]`/`[C]`/`[E]`，**不声称读到正文**）；`pwsh` = **OK**（`py -3.12 -I` 跑 `%TEMP%\o18_check.py`，打印 §5 全表，无错）。**无 `TOOL-FAIL`。**

| # | 检索式 | 范围 | 命中 | 结论 |
|---|---|---|---|---|
| 1 | `Cheeger\|cheeger\|isoperimetric\|等周` | 全仓 `*.md` | **250**（截断 256 中的 250） | 库内 Cheeger 判词全在**边际 `ν`** 上：`CHEEGER-BLIND`、`B3-ISOPERIMETRIC`、`BASE-NONLOCAL-BOUND`、`Q4-BAKRY-EMERY`、`H6-HEAD-SPECTRAL:87,115,316`、`SLOT-PF-FIXPOINT:103–111`。**条件测度版命中 0**（检索式：本行；工具：OK；命中：0）⇒ O18 是**新对象** |
| 2 | `h²/4\|h\^2/4\|λ₁≥h` | 全仓 `*.md` | 见 #1 的 250 内 | 常数 `1/(2π)` 已校准：`BASE-NONLOCAL-BOUND:38,50,69`；1D SU(2) `h=0.7759`、`h²/4=0.1505≤λ₁=0.2645`（`Q4-BAKRY-EMERY:313`） |
| 3 | `Brascamp\|BL\|Hess V ⪰` | 全仓 `*.md` | 见 #1 | `GPP-STRUCT §4(b)`：H1 真身 = BL/log-凹；BL 前提 = Savvidy 墙 |
| 4 | `Cheeger "A lower bound for the smallest eigenvalue of the Laplacian" Problems in Analysis 195 h^2/4` | web | 10 源 | DOI **`10.1515/9781400869312-013`**（De Gruyter）；`[CVGMT 4217](https://cvgmt.sns.it/media/doc/paper/4217/Buser@inequalityFinal.pdf)` 复述 `λ₁≥h²/4`（片段级） |
| 5 | `Brascamp-Lieb variance bound Hessian weight Var(f) ≤ ∫<(Hess V)^{-1}∇f,∇f>` | web | 10 源 | [Bonnefont 讲义](https://www.math.u-bordeaux.fr/~mibonnef/Poincare__Toulouse.pdf)：BL **极值函数 `f=⟨c,∇V⟩`**（片段级，本文件 §4.1 的关键一步）；另 [HAL 论文](https://hal.science/INSA-LYON-THESES/tel-05798816v1) 记"BL + Cramér–Rao 传递 Hessian 信息到协方差" |
| 6 | `Bakry-Emery Γ₂ criterion gap ess inf Ric + Hess V Lichnerowicz nK` | web | 10 源 | [arXiv:2405.10734](https://export.arxiv.org/pdf/2405.10734v2)：Lichnerowicz `λ₁ ≥ nK`（Ric ≥ (n−1)K）"well-known … sharp"（片段级） |
| 7 | `Levy-Gromov isoperimetric Ricci lower bound Cheeger constant sqrt kappa` | web | 10 源 | [Lyon 讲义](http://math.univ-lyon1.fr/wikis/Stab/lib/exe/fetch.php?media=wiki:seminarlyonmondino.pdf)：`m⁺(E) ≥ |∂B|/|S|`，球模型（片段级）⇒ 曲率 ⇒ `h` 下界 |
| 8 | `weighted manifold Cheeger inequality measure e^{-V} boundary measure` | web | 同 #4/#7 池 | 无独立"加权/条件测度"专文命中（检索式：本行；工具：OK；命中：0）⇒ §1.2 的 `h₀` 定义按标准余面积/条件化写，不引外部 |

**外部引用分级**：Cheeger 1970 = `[std]`（**DOI 书目级已核**，正文未逐字）；Buser 1982 / Jerrum–Sinclair / Sinclair = `[std]`（书目级）；Brascamp–Lieb 1976 = `[std]/[E]`（卷页未与出版社页逐条对）；Lichnerowicz/Bakry–Émery/Lévy–Gromov = `[std]`（片段级复述已核，原始正文未逐字）；BL 极值函数 = `[C]`（讲义片段）。

**标签汇总**：§1 `[std]`+`[T]`（条件化公式）；§2 `[T]`；§3 `[T]`；§4 `[T]/[std]/[NG]`；§5 `[T]+[V]`（witness 输出）；§6 表内逐行已标；§7 `[T]`+`[J]`；§8 未闭合 4 条已标。**未把任何数值/猜测标成 `[T]`。**
