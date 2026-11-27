# KKKK1-PHIDEF — `Φ_k(X)` 的逐字定义是 (5.31) 的**截断累积量级数**，域是**块的连通并**；单片 `{p}` 不在域内 ⇒ (c1)(c2)(c3) 全 `[NG]`

> 轮 KKKK1 · 子代理 `kkkk1-phidef`（委托 `lead`）· `D:\Dev\math`。只新建本文件（≤90 行）；既有文件零改动；未 `commit`；零新 MC。
> `fax:` = `research/v3/species/_ht7/faizal.txt:行`。标签只 `[T]`/`[NG]`；无百分数；每条带「下一手」。

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| **0.1** | `Φ_k(X;U\|X)` 逐字 = `Σ_{n≥1}(1/n!) Σ_{B_1..B_n⊂X 连通} ⟨V(B_1);…;V(B_n)⟩^T(U\|_{B_1∪…∪B_n})`（5.31），截断取在「限制到 `Ω_δ` 且插入 `R_k` 的归一化高斯」上 | **`[T]`**（逐字） | `fax:22002`(5.31)·`22007–22008`·`21912–21918`(5.21) | 引 `Φ_k` 必引 (5.31)＋(5.21) |
| **0.2** | 域 = 「有限连通**块并**」：`X` 是 `Λ_k` 上 `b^k` 块（4D 版 `B_k` 块）的连通并；`p`（plaquette）不是 `X` | **`[T]`** | `fax:21982`(5.29 前句)·`21988`(``vanishes on disconnected X'')·`22349–22356`(6.25) | 写 `Φ_k({p})` 前先声明域外插值 |
| **0.3** | `δ_p` 只进 `R_k` 的指数与 `Ω_δ` 判据，**不进** (5.31) 右端 | **`[T]`** | `fax:21922–21925`(5.22)·`21928–21931`(5.23)·`21945`(5.25)·`22371–22375`(6.27) | `δ_p` 记「regulator 核」列，不记「活动值」列 |
| **0.4** | (c1)(c2)(c3) **全 `[NG]`**：三者都不是 `Φ_k` 的值，而是「单片标量/其均值」 | **`[NG]`** | §2 | 三候选移出 `Φ_k` 定义栏，入「单片诊断」栏 |
| **0.5** | 前件 `‖Φ_k‖_α≤ε`：**条件式 `[T]`** —— `sup_B‖V(B)‖≤ε ⇒ ‖Φ_k‖_a≤Cε²≤ε`（`ε≤1/C`），常数与 `k` 无关 | **`[T]`**（带条件） | `fax:22015–22018`(5.32)·`22052–22059`(5.37) | 引用时并写 `sup_B‖V(B)‖≤ε` |
| **0.6** | 前件的**无条件**版 `[NG]`：种子 `sup_B‖V(B)‖≤ε` 在 `fax` 内是**假设**（`if`），`V(B)` 只有「slab 积分」定性描述，无显式式 | **`[NG]`** | `fax:21997–21998`·`22016`(`if`)| 下一手：找 `V(B)` 的显式式（本文件不造） |
| **0.7** | `U_p=−1` 见证：`δ_p=½‖1−U_p‖²_F=2N`（SU(2)：`4`）——这是 `δ_p` 的 `sup`，**不是** `‖Φ_k‖` 的见证 | **`[NG]`**（作前件反证） | `fax:21924`(`∈[0,2N]`)·`22711`(`ε∈(0,1)`）| 见证度记 `2N`，见证对象写 `δ_p` |

## §1 原文（行号 + 引文；OCR 折行已并）

**1.1 (5.29)/(5.30) `fax:21982–21996`**
```
A polymer X is a finite connected union of blocks in a fixed time layer. The effective
interaction on that layer admits the unique expansion  S_k(U) = Σ_X Φ_k(X;U|X)   (5.29)
where Φ_k(X;·) depends only on bonds in X, is gauge-invariant and reflection-covariant, and
vanishes on disconnected X.
∥Φ_k∥_a = sup_x Σ_{X∋x} ∥Φ_k(X)∥ e^{a diam(X)}                                    (5.30)
"with operator norm ∥·∥ induced by the Gaussian reference with covariance C_k^int and the
 regulator, and graph diameter diam(X)"                                          fax:21995
```

**1.2 (5.31) 逐字定义 `fax:21997–22009`**
```
Let {V(B;·)}_B be block-local base functionals obtained by integrating slab fluctuations
at fixed boundary fields. Define V(X;U|X)=Σ_{B⊂X} V(B;U|B) and set
Φ_k(X;U|X) = Σ_{n≥1} (1/n!) Σ_{B_1,…,B_n⊂X, connected} ⟨V(B_1;·);…;V(B_n;·)⟩^T (U|_{B_1∪…∪B_n})  (5.31)
with truncated expectation (5.21) taken w.r.t. the normalized Gaussian measure restricted to
Ω_δ with regulator R_k inserted.  By finite range of C_k^int and Lemma 5.2, (5.31) is well
defined and yields connected activities.
```
`(5.21)`（`fax:21912–21918`）：`⟨F_1;…;F_n⟩^T = ∂^n/∂t_1…∂t_n log E[exp Σ t_i F_i]|_{t=0}`。

**1.3 `δ_p` 的两处、也是全部去处 `fax:21922–21931`·`21945`**
```
δ_p(U) = ½‖1−U_p‖²_F ∈ [0,2N]                                                  (5.22)
R_k(U) = exp(−λ Σ_{p⊂Λ_k} φ(δ_p(U)))                                          (5.23)
Ω_δ = {U : δ_p(U) ≤ δ  for all plaquettes p ⊂ B}                              (5.25)
```
4D 版同构：`R_k(U)=exp(−λΣ_p φ(½‖1−U_p‖²_F))`，`S^F_k={U:‖1−U_p‖_F≤δ}`（6.27，`fax:22371–22375`）。

**1.4 4D 版活动与范数 `fax:22349–22370`**：`Sk(U)=Σ_{X⊂Bk, X connected} Φ_k(X;U|X)`（6.25）；`‖Φ_k‖_a = sup_{B∈Bk} Σ_{X∋B, X conn.} ess sup_{U|X}‖Φ_k(X;U|X)‖ e^{a diam(X)}`（6.26）；`‖·‖` = `C_k`＋regulator 的算子范数。

**1.5 前件的两处身份**：定理 5.3 `fax:22015–22018`「`if sup_B‖V(B)‖ ≤ ε then ‖Φ_k‖_a ≤ Cε²`」（5.32）；命题 7.2 把它当**假设**用：`fax:22708–22711`「We assume the following uniform small-activity hypothesis, which is guaranteed at each scale by the ultraviolet stability theorem … `‖Φ_k‖_a ≤ ε, ε∈(0,1)`」（7.19）。逐尺度：`fax:22440–22456`（`‖Φ_{k+1}‖_a ≤ C_1‖Φ_k‖²_a`）·`22482–22500`（`sup_k‖Φ_k‖_a ≤ 2ε`）。

## §2 三候选各判

| 候选 | 值 | 判 | 逐字依据／缺什么 |
|---|---|---|---|
| **(c1)** `Φ_k({p})=δ_p` | `sup=2N=4` | **`[NG]`** | ①域错：`X` 是块并（0.2），`{p}` 不是聚合物；②模型错：`δ_p` 只出现在 `R_k` 的指数与 `Ω_δ`（1.3），(5.31) 右端是 `V(B)` 的累积量；③`V({p})↦δ_p` 是**识别**（`HHHH4-IDENT.md:39` 转述 `BBBB1:41` 自注“假设”），非逐字。`grep Φk\(\{p\}\)`／`Φk = δp`：**0 命中** |
| **(c2)** `2I₂/I₁=1.014390` | `>1` | **`[NG]`** | 本式是 Haar Gibbs 下的 Bessel 比值读数，`fax` 内 `I₂/I₁` **0 命中**（`Bessel` 34 命中全为 `H^{−s}` 势空间一类，`I1/I2` 15 命中全在 `2027–2090` 作范数估计记号）；且 (5.31) 的测度是**截断高斯＋R_k**（1.2 末句），不是 Haar Gibbs |
| **(c3)** `2(1−I₂/I₁)=0.985610` | `<1` | **`[NG]`** | 同 (c2)：是单片标量在 Haar Gibbs 下的**均值**（`HHHH4-IDENT.md:13` 记 `⟨δ_p⟩`），不在 `Φ_k` 的定义链上；均值 ↛ 范数上界 |

**2.5 域内最接近的一支**：`X={B}`（单块）时，(5.31) 的 `n=1` 项 = `⟨V(B)⟩^T = E[V(B)]`（一阶累积量即均值），`n≥2` 项由 (5.37) 管（`fax:22052–22055`「summing over n ≥ 2 … `‖Φ_k(X)‖ ≤ C₂ε² e^{−c′diam X}`」）。故单块活动的量级是 `O(‖V(B)‖)` 与 `O(ε²)`，**既非 4、也非 1.014／0.986**。⇒ 三候选之外无第四支被 `fax` 写出：`fax` 只给 (5.31) 与 (5.32)。

## §3 `‖Φ_k‖_α ≤ ε` 真伪

**3.1 直接裁决**：作为**蕴含**它是 `[T]`（逐字 5.32：`sup_B‖V(B)‖≤ε ⇒ ‖Φ_k‖_a≤Cε²`，取 `ε≤1/C` 得 `≤ε`；`a<c′`，常数不依赖 `k`）。作为**无条件断言**它是 `[NG]`：种子 `sup_B‖V(B)‖≤ε` 在 `fax:22016` 是 `if`，`V(B)` 无显式式（`fax:21997–21998`）。**假**这一侧不成立——(c1) 的 `4` 打的是 `δ_p` 的 `sup`，打在 `R_k` 的核上，不进 (5.32) 的 `sup_B‖V(B)‖`。

**3.2 若强行按 (c1) 取见证**：`U_p=−1`（SU(2)，`U_p=−1_2`）⇒ `1−U_p=2·1`，`‖1−U_p‖²_F=tr(4·1)=8`，`δ_p=½·8=4=2N`。见证度 = `2N`（`N=2`：`4`）。该值 `>1=sup ε`，但它见证的是 `δ_p` 贴到区间上端（`fax:21924`），**不是** `‖Φ_k‖_α`；故不能作前件的反证 ⇒ `[NG]`。

**3.3 下一手（三条，互不依赖）**：①要判前件真伪，靶应改成 `sup_B sup_U‖V(B;U|B)‖≤ε`（`fax:22016` 的唯一入口）；②要到 `V(B)` 的显式式（`fax:21997` 只给定性描述，`grep "V (B;"` 6 命中全在 `21997–22030`，无展开式）；③`(7.19)` 的 `ε∈(0,1)` 与 `(5.32)` 的 `Cε²` 需并引，否则 `ε` 会被读成两个不同的门。

## §4 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 实际检索式／调用 | 命中 |
|---|---|---|---|---|
| 1 | `read` | **OK**（带行号正文全返回，32615 行文件） | `fax:21910–22019`·`22019–22110`·`22340–22384`·`22700–22729`；`HHHH4-IDENT.md:1–60` | 5/5 |
| 2 | `grep` | **OK**（ripgrep 正常，无 `provider failure`／`ENOENT`／空响应） | ①`Phi_k\|Φ_k\|\\\\Phi_k\|S_k\|R_k` ②`Φk\|base functional\|truncated` ③`V \(B\|V \(X\|V \(B;` ④`Bessel\|I1 \(\|I2 \(\|modified` ⑤`Φk \(\{p\}\)\|Φk = δp\|δp = Φk\|I2\|I1\|tr U` ⑥`hhhh4`@`*.md` | ①**0** ②**151** ③**6** ④**34** ⑤**15**（其中 `Φk({p})`／`Φk=δp`／`δp=Φk`＝**0**）⑥**1**（`HHHH4-IDENT.md:3`） |
| 3 | `glob` | **OK** | `research/v3/species/**/*hhhh4*` | **0** 文件（该名只作正文出现） |
| 4 | MC／Lean／`git`／`web_*`／`skill` | **0 调用** | 本靶＝仓内逐字核；catalog 无匹配技能 | — |

**诚实边界**：①行号对齐本仓文本副本，未核出版页码；②(5.31) 的 `Σ_{n≥1}` 与 `Σ_{B_1..B_n⊂X connected}` 在 OCR 里被折成两行、`n≥1` 与 `connected` 分置，本文件按 `22004–22006` 直录，未补字；③「(c1) 不是逐字」这一判只依赖域（0.2）与 `δ_p` 去处（0.3），不依赖 `V↦δ_p` 那条识别；④无 QED；未证之物未标 `[T]`；未改任何既有文件。

*KKKK1-PHIDEF · 只新建本文件。判：`Φ_k` 逐字 = (5.31) 截断累积量级数、域为块连通并（`[T]`）；(c1)(c2)(c3) 全 `[NG]`（域错／测度错／均值≠范数）；前件条件式 `[T]`（5.32）、无条件 `[NG]`（种子是 `if`）。*
