# S08-BL — Brascamp–Lieb 能否给出 (★) 的 `−δ_W` 余量？

> 目标 (★)：`Var_μ(⟨Pᵀ∇S_j,v⟩) ≤ E_μ⟨PᵀH_jPv,v⟩ + W_j(v) − δ_W`，`μ ∝ e^{−A}`，`A = S_j + Φ`。
> **判词（一句）**：**(b)。核条件后 BL 不适用，且更强——BL 在此为「假」，不是「弱」。** BL 的假设是逐点 `Hess V ≻ 0`（严格对数凹；本 brief 的 `D²V ≥ K > 0` 是更强者，且 `K` 对精确形式并非必需），而 Savvidy `Hess A = PᵀH_jP + Φ″` 不定（`λ_min = −g|B|`）；brief 里推出 (★) 所需的那个泛函不等式**可证伪** `[T]`。非凸替代存在（精确半群恒等式 / Poincaré / 负模秩一 SCBL），但余项与 `−δ_W` **反号**，或把 (★) 换成 H1。BL 路线终点 = `W_j ≥ (Var−E[H]) + δ_W` = (◆) = **H1**。

---

## §1 · 精确陈述与条件（先核，后谈应用；`AGENTS §2.1`）

### 1.1 逐字核对：条件是「严格对数凹」，**不是**「一致强凸」

三个独立源一致：

| 源 | 逐字 | 条件 |
|---|---|---|
| Bonnefont 讲义（[Poincaré/Toulouse](https://www.math.u-bordeaux.fr/~mibonnef/Poincare__Toulouse.pdf)，读到正文 §"another proof of the Brascamp–Lieb variance inequality"） | `Let dμ = e^{−V}dx be a strictly log-concave measure on R^d ... Here strictly log-concave means Hess V(x) > 0 for all x.` | 逐点 `≻0` |
| *Ann. Probab.* **50** (2022) no.2（[DOI 10.1214/21-AOP1545](https://projecteuclid.org/journals/annals-of-probability/volume-50/issue-2/…/10.1214/21-AOP1545.pdf)，读到正文 BL 陈述句） | `Var(⟨η,X⟩) ≤ E⟨η,(Hess f)(X)^{-1}η⟩` | 严格凸 `f` |
| Brascamp–Lieb 1976（Adv. Math. **20**, 151–173，Theorem 4.2；经 [arXiv:2511.03193](https://ar5iv.labs.arxiv.org/html/2511.03193) 引作 `Brascamp and Lieb, 1976, Theorem 4.2`，`for a (strictly) log-concave measure p(x)`） | 同上 | 严格对数凹 |

**精确形式（`[T]`）**：设 `V ∈ C²`、`Hess V(x) ≻ 0` 对所有 `x` 成立（故可逆）、`Var_μ(h) < ∞`、`h` 局部 Lipschitz，`μ = e^{−V}dx/Z`，则

```
Var_μ(h)  ≤  ∫ ⟨∇h(x), (Hess V(x))^{−1} ∇h(x)⟩ dμ(x) .
```

### 1.2 条件的三层拆分（**关键更正**）

1. **本质条件**：逐点 `Hess V ≻ 0`（在每个 `{∇h ≠ 0}` 的点上）。这是 BL 的**充要级**假设——它来自一维引理 `Var ≤ ∫(h′)²/V″`，只需 `V″ > 0`。
2. **非本质条件**：`Hess V ⪰ K > 0` 的**一致** `K`。它只用于把 RHS 压成推论 `Var ≤ K^{-1}E|∇h|²`（= Poincaré，`[T]`）。**brief 的 `D²V ≥ K > 0` 是加强版，不是 BL 的原假设。**
3. **失败模式**：凸性一旦失效，BL **为假**（§2 B1 给出反例），而非"常数变差"。⇒ 「把 `K` 放松到逐点」**救不了**任何 YM 应用——因为逐点形式也假。

### 1.3 一条决定全部结论的尺度事实（`[T]`）

BL 的权重对 Hessian 是 **`−1` 次齐次**：`Hess ↦ tHess` 使 RHS `↦ t^{-1}RHS`。因此

- Hess 有**小**特征值（IR 软模、`λ↓0`）⇒ BL 的界 **爆炸**，恰在 (★) 最需要余量的方向最弱；
- Hess 有**负**特征值 ⇒ 该方向 RHS **为负**，不能上界非负的方差。

(★) 要的是 `+1` 次（`E⟨PᵀH_jPv,v⟩`）。**次数错配是结构性的，不是技术性的。**

---

## §2 · Hess 不定如何破坏 BL —— 三重破坏，逐条

### (B1) 反例：BL 的结论**为假** `[T]`

一维双井 `V(x) = x⁴/4 − x²/2`：`e^{−V} = exp(−x⁴/4 + x²/2) ∈ L¹`（可积、`Z<∞`），`V″(x) = 3x² − 1`，在 `I = (−1/√3, 1/√3)` 上 `V″ < 0`。

取光滑 `h` 使 `h′` 支撑在 `I` 的紧子区间内且 `h` 非常数（`h` 由 `h′` 积分得到，不必单调）。则

```
BL RHS = ∫ (h′)²/V″ dμ  <  0  <  Var_μ(h) .
```

**矛盾。** 无需奇点、无需发散、无需 `h′` 与 `V′` 相关。⇒ 不定 Hess 下 BL 不是"未证/常数差"，是**假命题**。

（附带，同一例的第二种失败形态：若取 `h = x`，`∫(x³−x)²/(3x²−1)dμ` 在 `x=±1/√3` 处左右异号对数发散——因为分子 `x²(x²−1)²|_{x₀}=4/27 ≠ 0`。即 `Hess` 穿零时 BL 的 RHS 甚至**未定义**。这与 Savvidy 的 `λ_min` 随背景场连续穿零完全同构。）

### (B2) 结构破坏：`Hess A` 穿零且负 ⇒ 无一致常数 `[T]` + `[C]`

`Hess A = PᵀH_jP + Φ″`。Savvidy 精确符号 `λ_min(PᵀH_jP) = −g|B|`（`PROOF-SAVVIDY-EXACT-SIGN`）。`λ_min` 随背景场/壳层指标连续变化 ⇒ **必穿过 0** ⇒ `(Hess A)^{-1}` 在塔上无界 ⇒ 任何"逐壳 BL + 输运"都不可能给出与 `a, Λ` 无关的常数。（β≈2.5 / SU(3)：`min PᵀH_jP ≈ −0.3`、21–42% 负特征值、90–97% 集中在 UV 块、`−1.12β`；`[C]`/`[V]`，窗口见 §7。）

### (B3) 形式破坏：次数错配 ⇒ 需一步 Loewner，而该步的前提正好是 Savvidy 墙 `[T]`

把 `f = ⟨Pᵀ∇S_j, v⟩`（`∇f = PᵀH_jPv`）代入 BL：

```
Var(G_v) ≤ E⟨ PᵀH_jPv , (Hess A)^{−1} PᵀH_jPv ⟩ .            (2.1)
```

要得到 (★) 的 `E⟨PᵀH_jPv,v⟩`（`Hess^{+1}`），唯一通道是 Loewner

```
(Hess A)^{−1} ⪯ H_j^{−1}   ⇔   Hess A ⪰ H_j ⪰ 0 ,           (2.2)
```

即 BL 的"等分成" `Var ≤ E[H_j]`（`S6-DOBRUSHIN §5` (5.1)→(5.2)）。`(2.2)` 的两个合取项在 head **同时为假**：`H_j` 不定（`min ≈ −0.3`），且 `Φ″` 的凸性预算不足以压过负模（`W_full ≈ 0.22–0.31` 对 `−(E[H]−Var)`）。

### (B4) brief 里"须要的泛函不等式"**可证伪** `[T]`（本轮最硬死线）

brief 的推导链要求

```
⟨∇f, (D²V)^{−1}∇f⟩  ≤  ⟨PᵀH_jPv, v⟩ + W_j − δ_W .            (2.3)
```

令 `w` 为 `Hess A` 最小特征向量的像、`λ = λ_min(Hess A) ↓ 0⁺`：LHS `= ‖w‖²/λ → +∞`，RHS `= O(1)`。**逐点假。** 若 `λ<0`，LHS `< 0` ⇒ 连"非负上界"都不是。⇒ `(2.3)` 不可能由任何"加强版 BL"获得；`(★)` 不能经 BL 得到。

---

## §3 · 酸测：BL 的**最好情形恰好是自由场** `[T]`

`Φ ≡ 0`（无 Haar）、`S = ½xᵀHx`（`H ⪰ 0` 自由 Laplacian）：`Var(G_v) = vᵀHv = E[H_vv]` **取等**，`W = 0` ⇒ `G″ = E[H] − Var + W = 0` ⇒ **(★) 失败**（`S6-DOBRUSHIN §6`）。

即：BL 的上界被**等号饱和**时，(★) 恰被否证。故

- BL 对 (★) 的**结构性贡献 = 0**；(★) 的全部内容只能是 `W_j` 的**正余量**（Prékopa/`Φ″ ⪰ N/12`，Lemma B）扣掉 Savvidy 超额：`W_j > (Var−E[H]) + δ_W` ≡ **(◆) = H1**（`S6 §5`；`PROOF-JACOBIAN-SAVES §2`）。
- **`δ_W` 不来自 BL**；`W_j` 也不来自 BL（BL 只用 `Φ″` 的凸性去**取消**饱和，从不产生 `+W`）。

---

## §4 · 非凸版的方差界：有，但都救不了 `−δ_W`

### (N1) 精确半群恒等式（**不需要任何凸性**）`[T]`

```
Var_μ f  =  2 ∫_0^∞ E_μ |∇P_t f|² dt ,      μ 可逆半群 P_t = e^{tL} ,  L = Δ − ∇V·∇ .
```

这是"非凸版方差界"的**精确**形态，且对任意 `V`（含 Savvidy 不定）成立。但闭合它需要

```
‖∇P_t f‖₂ ≤ e^{−ρ t} ‖∇f‖₂   ⇔   Γ₂ ⪰ ρ   ⇔   Hess A ⪰ ρ > 0 ,
```

即**回到被 (B1)/(B4) 否掉的假设**。⇒ (N1) 是**换名**（把 (★) 变成"半群收缩"= 粗模隙 = H1），不是替代。

### (N2) Poincaré 版 `Var ≤ C·E|∇f|²`（用户问的那条）`[T]-标准`

- 对**任意**概率测度都可定义 `λ₁ = inf ℰ/Var`，`C = 1/λ₁ < ∞ ⟺ λ₁ > 0`；非凸 `V` 也能有 `λ₁ > 0`（固定有限势垒的对称双井即是）。
- 但**常数无免费上界**：双井/亚稳 `λ₁ ≍ e^{−2h/ε}`（Helffer–Sjörstén 半经典 1984；Kramers/Eyring）。在 RG 塔上势垒 `h ∝ β` ⇒ `C` 关于 `β` **指数退化** ⇒ **无塔一致的 `δ_W`**。
- 且它 sign-blind 于 `G″`（只能给 `Var(a) ≤ 1/λ` 的**必要**型界），本仓已多次判为 H1 的必要面。

### (N3) 负模精确秩一（SCBL，本仓已有）`[T]`

`K = K₊ + λ₋ φ₋φ₋ᵀ`（`K₊ ⪰ 0`，`λ₋ < 0`，秩一）。BL 只能用在 `K₊` 上：

```
Var ≤ E⟨ K₊⁺c, c ⟩ + E|⟨φ₋, c⟩|² / |λ₋| ,        c := PᵀH_jPv .     (4.1)
```

修正项 `≥ 0`（`D = D₊ + D_Sav`，`D_Sav = E|⟨φ₋,c⟩|²/|λ₋| ≥ 0`；`SAV-CANCEL`、`SM2-RANK-ONE`），且 `λ₋ → 0⁻` 时 `→ +∞`。
**符号与 `−δ_W` 相反**：不定性只能**加**缺陷。这是"非凸替代"最深的一条，也是它**注定不能**给 (★) 的原因。

### (N4) Houdré–Kagan 不是非凸版 `[T]-摘要级`

HK = `Var G(X) ≤ E(G′(X))²`，**仅对标准 Gaussian**（= BL 在 `Hess V = I`）；其证明 `based on properties of the Gaussian density`，且被反复表述为 Gaussian 的**等号/刻画**型结果（HK, *J. Theoret. Probab.* **8** (1995) 23–30；摘要见 [NII-REO](https://reo.nii.ac.jp/oja/ART1001238138/)；KU 论文标题即 "…Houdré–Kagan's variance inequality **for Gaussian measure**"）。它**不含**任何 `Hess V` 的松弛 ⇒ **不能当非凸替代**。

### (N5) 其余候选，一句话各自判掉

| 候选 | 为什么不给 `−δ_W` | 级 |
|---|---|---|
| Efron–Stein `Var ≤ Σ influence²` | 无条件但 **sign-blind**（对 `G″` 的方向一无所知）；本仓已判 wrong-sign | `[T]` |
| LSI / Bakry–Émery | 需 `Ric + Hess ⪰ ρ`，即 `Hess A ⪰ ρ`（W6：`≈ −2.3`） | `[T]`/`[NG]` |
| Prékopa / Borell | 只给 `Φ` 的预算 `W_j ≥ N/12`（**正项，留用**），对 `H_j` 无任何界 | `[T]` |
| Simon–Lieb | 方差 = H1（sign-blind），不产生符号信息 | `[T]` |
| 加权 BL / 非对称 BL（Otto–Menz） | 仍需 `Hess ⪰ 0`；非对称版只放宽"次数配对"，不放宽符号 | `[C]` |

---

## §5 · 判定 (a) / (b) / (c)

### (a) 能否（即使 Hess 不定）给出 (★) 的条件与推导？ — **否** `[T]`

若强行写出"要成立须什么"：需同时

1. `Hess A ⪰ 0`（去掉负模，或负模被 `Φ″` 逐点压过）；
2. `Hess A ⪰ H_j`（Loewner，`(2.2)`）；
3. `λ_min(Hess A) ≥ c > 0` **且** `c` 与 `β, Λ, a` 无关（才有一致 `δ_W`）；
4. `W_j ≥ (Var−E[H]) + δ_W`（正余量，`(◆)`）。

在 Savvidy head，1–3 全假；而 4 就是 **(◆) = H1**。**⇒ 无 (a)。**（唯一例外窗口：空锥深 UV 壳，`H_j ⪰ cΛ_j² > 0`，那里 `(2.2)` 成立、(★) 为 `[C]`——条件于 `E|F_j| ≤ Cg_j`，见 `S6 §5`。这不是 BL 的功劳，是"负模被 Laplacian 地板压住"的功劳。）

### (b) 不能（因 `D²V ≥ K` / `≻0` 不成立）⇒ **「BL 用不上」** `[T]`

并且比 brief 的预期**更硬**：不是"条件不满足所以借用"，而是
**BL 在 Savvidy 型不定 Hess 上为假**（B1 反例），且 brief 推导链所需的 `(2.3)` **逐点可证伪**（B4，`λ↓0⁺ ⇒ LHS→∞`）。

### (c) 非凸替代 — **有陈述，但都不给 `−δ_W`** `[T]`/`[C]`

| # | 陈述 | 给什么 | 判 |
|---|---|---|---|
| (N1) | `Var f = 2∫₀^∞ E|∇P_t f|²dt` | 精确、无凸性 | **换名**成 H1（需 `Γ₂ ⪰ ρ`） |
| (N2) | `Var ≤ C E|∇f|²`，`C = 1/λ₁` | 任意 μ 可定义 | `C ≍ e^{2h/ε}` ⇒ **非塔一致** |
| (N3) | `Var ≤ E⟨K₊⁺c,c⟩ + E|⟨φ₋,c⟩|²/|λ₋|` | 负模精确记账 | 余项 **`≥0`，反号** |
| (N4) | HK: `Var G ≤ E(G′)²`（Gaussian） | Gaussian 等号 | **非非凸**，不可用 |

**结论：不存在"非凸 BL"能为 (★) 提供 additively-negative 的 `−δ_W`。** 所有非凸替代要么精确但把问题换成"半群收缩"，要么其不定性修正项**与所需余量反号**。

---

## §6 · 路径（本轮之后的下一靶）

1. **入库判死**：死路台账新增一条 —— `BL ⟹ (★) 的 −δ_W`：**假**（B1 反例 + B4 可证伪），原因不是 `D²V ≥ K` 太强，而是 `Hess ≻ 0` 都假。窗口：任意 `β`、任意维、任意 `N`（B1 结构性）。**不得**再以"BL 加上更强的凸性假设"复活。
2. **保留的唯一正预算**：`W_j`（Prékopa + Lemma B，`Φ″ ⪰ N/12`，点态、与 `β` 无关）。把 (★) 精确重写为
   `W_j ≥ (Var − E[H]) + δ_W` ≡ **(◆) = H1**（Savvidy 超额必须由 Haar 预算偿付）。
3. **方差不等式物种的适用域**：只在空锥深 UV 壳（`H_j ⪰ cΛ_j²`）有效 `[C]`。head 上必须换**非方差**物种（负模秩一记账 + 逆问题/约束有效势），因为 variance-bound 家族在不定 Hess 上**结构性**给不出符号。
4. **酸测门槛**：任何未来候选必须满足——(i) 在自由场（`W=0`、`G″=0`）**失效**；(ii) 在 B1 型双井上**不被反例击穿**（即不依赖 `Hess ≻ 0`）。BL 两条都不满足。

---

## §7 · 限定（`AGENTS §2.1`）

- **(B1)/(B4) 的窗口**：任意维（把 `h′` 投到负特征子空间即逐字搬用）；任意 `β`；与 `N`、`L`、格距无关。`[T]`
- **(B2) 的数值窗口**：β = 2.5、SU(3)、L = 4/6，`min PᵀH_jP ≈ −0.3`、`W_full ≈ 0.22–0.31`、21–42% 负、90–97% UV 块 `[C]/[V]`；`λ_min = −g|B|` 是精确符号 `[T]`。**不外推到 β<1/48**（那里 `H_j ⪰ cΛ²`，(★) 走 (c)/§5(a) 例外窗口）。
- **(N2) 的指数退化**是半经典标度下的结论；固定势垒的单个双井 `λ₁` 是正常数——本文件只用"塔上势垒 `∝β` ⇒ 非一致"这一点。
- **未逐字核**：Brascamp–Lieb 1976 原文正文（Adv. Math. 20, 151–173 Theorem 4.2 —— 经二手引文核对定理号，未取到原文 PDF）；Bonnefont 讲义 PDF 内容类型不支持抓取，仅逐字用了搜索返回的原句。见 §9 检索核验。

---

## §8 · 与既有内部结论的一致性（不冲突，只是收紧）

| 内部文件 | 既有说法 | 本文件 |
|---|---|---|
| `S6-DOBRUSHIN §5` | BL 给 `(5.1)`，collapse 到 `(5.2)` 需 `H_j ⪰ 0`、`Φ″ ⪰ 0` | 同判；**新增**：collapse 之前 BL 本身在不定 Hess 上已假（B1/B4） |
| `S1-CONDITIONAL-POINCARE §1.7` | `D ≤ 0` 半（BL 半边）有效，Savvidy 半边 `≥0` | 同判；本文件把"BL 半边"的适用域压到 `Hess ≻ 0` 逐点 |
| `SAV-CANCEL` / `SM2-RANK-ONE` | `D = D₊ + D_Sav`，`D_Sav ≥ 0` | 同判；本文件用它作为 (N3) 的"反号"证据 |
| `ABC-NEW-INEQ` / `ABC-FINAL-LEDGER` | "Brascamp–Lieb needs `Hess V ⪰ 0`（the Savvidy wall）" | 同判，并补上**为什么加强 `K` 无用**（§1.2）与**为什么是假而非弱**（§2 B1） |

---

## §9 · 检索核验（`AGENTS §1`：工具状态 / 检索式 / 命中数 / 失败标记）

| # | 工具 | 检索式 / 目标 | 状态 | 命中 | 备注 |
|---|---|---|---|---|---|
| 1 | `grep` | `Brascamp\|Brascamp–Lieb\|Houdr\|Poincar`，`research/**/*.md` | **OK** | **250**（共 907，已截断，全文存 spill） | 有效响应，非 provider failure |
| 2 | `pwsh` | `Get-ChildItem research\v3\species` | **OK** | 约 380 项 | 用于确认 `S08-*.md` 命名与避免覆盖 |
| 3 | `read` | `research/v3/S6-DOBRUSHIN.md` L196–265 | **OK** | 70 行 | `(5.1)/(5.2)`、`W_j`/`δ_W` 术语来源 |
| 4 | `web_search` | `Brascamp-Lieb variance inequality Hess V positive definite exact statement …` | **OK** | 10 源 | 取得 Bonnefont 讲义原句、Ann. Probab. 陈述、arXiv:2511.03193 的 Theorem 4.2 引法 |
| 5 | `web_search` | `"Brascamp-Lieb" "strictly log-concave" Var h Hess inverse theorem` | **OK** | 10 源 | 确认 `strictly log-concave means Hess V(x) > 0` |
| 6 | `web_search` | `Houdré Kagan variance inequality Gaussian characterization` | **OK** | 10 源 | HK 的 Gaussian 专属（NII-REO 摘要 + KU 论文题） |
| 7 | `web_fetch` | `sl2x.aimath.org/.../1106_0709/section1.html` | **TOOL-FAIL：HTTP 404** | 0 | 该页已迁移；未据此下任何结论 |
| 8 | `web_fetch` | `ar5iv.labs.arxiv.org/html/1404.5886`（Saumard–Wellner） | **OK**（截断） | — | 未取到 BL 定理号 ⇒ 不引号数，只引"§4 方差界一节，未逐字核定理号" |
| 9 | `web_fetch` | `math.u-bordeaux.fr/~mibonnef/Poincare__Toulouse.pdf` | **TOOL-FAIL：unsupported content type `application/pdf`** | 0 | 原文句取自检索 #4 的返回；标"讲义/摘要级" |
| 10 | `web_fetch` | `ar5iv.labs.arxiv.org/html/2511.03193` | **OK** | — | 取得 `Brascamp and Lieb, 1976, Theorem 4.2` 引法与 `(strictly) log-concave` 措辞 |

**未查到（检索式：#1 + #4；工具：OK；命中：0）**：本仓无任何文件主张"BL 在 `Hess` 不定时仍可用于 YM"。反之 #1 命中的 27 处均把 BL 与 `Hess ⪰ 0` 绑定（`ABC-FINAL-LEDGER`、`BHH-TRANSLATION`、`BLAST-YM-NONLINEAR`、`S1-CONDITIONAL-POINCARE`、`Q6-H1-STATUS` 等）——与本判词一致。

---

## §10 · 引用（`AGENTS §3`：arXiv ID / DOI + 读到哪一层）

| 对象 | 出处 | 读到 | 级 |
|---|---|---|---|
| BL 方差不等式 | Brascamp & Lieb, *Adv. Math.* **20** (1976) 151–173, **Theorem 4.2**（二手引法核对；原文未取到） | 摘要/二手 | `[C]`（号数） |
| BL 原假设逐字 | Bonnefont 讲义，<https://www.math.u-bordeaux.fr/~mibonnef/Poincare__Toulouse.pdf>：`strictly log-concave means Hess V(x) > 0` | 搜索返回原句（PDF 未抓取） | `[C]` |
| BL 陈述句 | *Ann. Probab.* **50** (2022) no.2，DOI `10.1214/21-AOP1545` | 正文句 | `[T]` |
| BL 在"严格对数凹"下的陈述 | Saumard & Wellner, arXiv:**1404.5886** §4（BL/方差界一节） | 正文（截断），定理号未核 | `[C]` |
| HK 方差不等式（Gaussian） | Houdré & Kagan, *J. Theoret. Probab.* **8** (1995) 23–30 | 摘要（NII-REO）+ 二手 | `[T]-摘要级` |
| 半群方差恒等式 `Var = 2∫E|∇P_t f|²` | Bakry–Gentil–Ledoux, *Analysis and Geometry of Markov Diffusion Operators* (2014), Γ-计算 | 标准 | `[T]` |
| 双井谱隙指数小 | Helffer & Sjörstén, *Comm. PDE* **9** (1984) 337 | 标准 | `[T]-标准` |
| Savvidy 负模精确值 | Savvidy, *Phys. Lett. B* **71** (1977) 133；本仓 `PROOF-SAVVIDY-EXACT-SIGN` | 内部 | `[T]` |
| Prékopa / `Φ″ ⪰ N/12` | Prékopa, *Acta Sci. Math. Szeged* **32** (1971) 301；`PROOF-JACOBIAN-SAVES` Lemma B | 内部 | `[T]` |

---

## §11 · 一行结论

**BL 的假设是 `Hess V ≻ 0` 逐点（不是"一致 `K`"）；Savvidy 的 `Hess A` 不定 ⇒ BL **为假**（B1 一维双井反例：`∫(h′)²/V″ < 0 < Var h`），其推导链所需的 `⟨∇f,(Hess)^{-1}∇f⟩ ≤ ⟨H_j,vv⟩ + W − δ_W` 逐点可证伪（`λ↓0⁺ ⇒ ∞`）。非凸替代三条都有，但余项全为 `+`（`D_Sav ≥ 0`）或把 (★) 换名成半群收缩 = H1。⇒ (b)「BL 用不上」，且比 (b) 更硬：**BL 在此不成立**。`δ_W` 只能来自 `W_j` 的 Prékopa 预算扣 Savvidy 超额，即 (◆)=H1。**
