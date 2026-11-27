# LLLL1-VB — `V(B)` 的显式式：0 命中；只有定性句，幅度前件是假设，`ε` 无数值 ⇒ `[NG]`

> 轮 LLLL1 · 子代理 `llll1-vb`（委托 `lead`）· `D:\Dev\math`。只新建本文件（≤90 行）；既有文件零改动；未 `commit`；零新 MC。
> `fax:` = `research/v3/species/_ht7/faizal.txt:行`（本仓文本副本直读，未核出版页码）。标签只 `[T]`/`[NG]`；无百分数；每条带「下一手」。

## §0 判词

- **0.1 `[NG]`（显式式）**：`V(B)` 的唯一逐字定义是「block-local base functionals obtained by integrating slab fluctuations at fixed boundary fields」（`fax:21997–21998`）；**无方程**。检索：`V\s*\(B[^)]*\)\s*=` → **0**；`V\s*\(B[^)]*\)` → **6**（21997·21998·22002·22016·22022·22030）。下一手：写出 slab 积分式（测度/作用量/边界条件/归一化）。
- **0.2 `[T]`（身份）**：`V` **不是** (5.29) 里的 `Φ_k`（`S_k=Σ_X Φ_k`）；`V` 是 (5.31) 的**累积量种子**，`Φ_k` 由 `V` 的截断期望生成（`fax:21985`·`22002`·`21997–21998`）。下一手：引 (5.29) 时不把 `V` 当 `Φ_k`。
- **0.3 `[T]`（n=1 项）**：`V(B)` 不是 (5.31) 的 `n=1` 项（对 `B_1` 的求和）的被加项；该被加项是 `⟨V(B_1)⟩^T = E[V(B)]`（(5.21) 一阶累积量=均值；`fax:21912–21918`·`22002`）。下一手：写「`V(B)` 是 n=1 项的被平均对象」。
- **0.4 `[NG]`（居中）**：证明只算 `n≥2`（`fax:22052–22055`），(5.31) 却从 `n≥1` 起（`22002`）；若 `E[V(B)]≠0`，`n=1` 项未被 (5.37) 控制，`(5.32)` 的 `Cε²` 缺这一步。`E[V(B)]=0` / 显式减除：`fax` **未写**（`E\s*\[\s*V\s*\(B` 0 命中；`V\s*\(B[^)]*\)\s*⟩\s*=\s*0` 0 命中；`centered|zero mean` 17 命中，0 处在 21997–22030）。下一手：在 (5.31) 写 `V(B)−E[V(B)]`，或改 `n≥2`＋单块项另算。
- **0.5 `[NG]`（有界）**：`sup_B sup_U‖V(B)‖ ≤ ε`（22016 逐字只写 `sup_B`，`sup_U` 见 0.6）在 `fax:22016` 是 **`if` 假设**，不是定理；`V(B)` 无显式式 ⇒ `sup` 不可算。下一手：证 `sup_B sup_U|V(B;U|B)| ≤ ε`，或先给 `V(B)` 的显式式。
- **0.6 `[T]`（ε 约束）**：`∃ε>0`、`C<∞`、二者与 `k` 无关；**无任何数值**；`22016` 逐字**只有 `sup_B`，没有 `sup_{U|B}`**；「算子范数」句在 `21995`（不在 `22016`），`sup_U` 只在 `22361` 的 `ess sup_{U|X}` 出现（对象 `Φ_k`，`ess sup` 全篇 1 命中）。下一手：引 `ε` 写「存在、无值、与 k 无关」；`sup_U` 标为解释/`(6.26)` 类比。

## §1 原文（行号 + 引文；OCR 折行已并）

**1.1 `S_k` 与范数（`fax:21982–21996`）**
```
A polymer X is a finite connected union of blocks in a fixed time layer. The effective
interaction on that layer admits the unique expansion  S_k(U) = Σ_X Φ_k(X;U|X)   (5.29)
where Φ_k(X;·) depends only on bonds in X, is gauge-invariant and reflection-covariant, and
vanishes on disconnected X.
∥Φ_k∥_a = sup_x Σ_{X∋x} ∥Φ_k(X)∥ e^{a diam(X)}                                    (5.30)
"with operator norm ∥·∥ induced by the Gaussian reference with covariance C_k^int and the
 regulator"                                                                        fax:21995
```
**1.2 `V(B)` 的唯一逐字定义（`fax:21997–21999`；无方程）**
```
Let {V(B;·)}_B be block-local base functionals obtained by integrating slab fluctuations
at fixed boundary fields. Define V(X;U|X)=Σ_{B⊂X} V(B;U|B) and set
```
（`V(B;·)` 的积分式、测度、归一化、边界场依赖、`k`-标度：**均未给出**。）
**1.3 (5.31) 与 `n=1` 项（`fax:22002–22008`）**
```
Φ_k(X;U|X) = Σ_{n≥1} (1/n!) Σ_{B_1,…,B_n⊂X, connected} ⟨V(B_1;·);…;V(B_n;·)⟩^T (U|_{B_1∪…∪B_n})  (5.31)
with truncated expectation (5.21) taken w.r.t. the normalized Gaussian measure restricted to
Ω_δ with regulator R_k inserted.
```
`n=1`：`⟨V(B_1)⟩^T = ∂_t log E[e^{tV(B_1)}]|_{t=0} = E[V(B_1)]`（`fax:21912–21918` (5.21)）。
**1.4 前件与结论（`fax:22015–22018`）**
```
Theorem 5.3 (Scale-k polymer bound). There exist a > 0, ϵ > 0, and C < ∞, independent
of k, such that if supB ∥V (B)∥ ≤ ϵ, then  ∥Φk ∥a ≤ C ϵ2 .   (5.32)
```
**1.5 证明只算 `n≥2`（`fax:22052–22055`）**
```
... summing over n ≥ 2 and connected block families with ∥V ∥ ≤ ϵ, and choosing ϵ sufficiently small gives
∥Φk (X)∥ ≤ C2 ϵ2 exp(−c′ diam(X)).   (5.37)
```
**1.6 树界乘积（`fax:22030`）**：`⟨V(B_1);…;V(B_n)⟩^T ≤ ‖V(B_1)‖…‖V(B_n)‖ Σ_T Π α(B_i,B_j)`（(5.34)）；`‖V(B_i)‖` 只在此处与 `22052–22053` 进入不等式。

## §2 显式形：**无**（缺 5 项）

**2.1 一行**：`V(B;U|B) = ?` — `fax` 内**没有**这个方程；`V(B)` 是「块局部基泛函/块活动」，由 slab 涨落积分得到，积分式未写。

**2.2 缺什么**（逐项可 grep）：
1. **积分式**：`V(B;U|B)=∫_{slab interior}…dμ(…)` 未给出（`slab fluctuation` 仅 `21997` 1 处；`base functional` 仅 `21997` 1 处）。
2. **测度/作用量/边界条件**：`dμ_slab`、`S_W^slab`、`J^slab` 在 `21864–21873` 出现，但**未**连到 `V(B)`。
3. **归一化/居中**：`E[V(B)]=0` 未写；由 `(5.32)` 的 `ε²` 与 `22052` 的 `n≥2` **反推**需要（0.4）。
4. **范数对象**：`‖V(B)‖` 用 `21995` 的「Gaussian reference + regulator 的算子范数」，但该定义写在 `Φ_k` 的 (5.30) 里；`V(B)` 是泛函，乘子算子范数=ess sup，文本未点明。
5. **`k`-标度**：块 `B` 边长 `b^k a`（`fax:21898–21900`），`sup_B` 的 `k`-一致性未给；`ε` 只是「与 k 无关」的假设常数。

**2.3 最近的显式类比（不是 `V(B)` 的定义）**：`fax:5437–5443`（Section 9）有 `V_{k,X}`（连通聚合物、`diam=O(1)` 粗单位）满足 `‖V_{k,X}‖≤C0 e^{−σ0 diam(X)}`，且 `U_k` = 乘 `exp(Σ_X V_{k,X})` 后作 Gaussian re-centering——结构与 (5.31) 最近；`fax:11932–11933`（(F.15)）有 `‖V_{k,X}‖≤C0 e^{−μ_FRD diam(X)}`；`fax:5518`（(9.19)）有 `‖Φ_{k,X}‖_∞≤C0 e^{−σ diam(X)}`。三者都**不是** Section 5 的 `V(B)`。⇒ 下一手：核 `V(B)` 是否可**认同**为 Section 9 的 `V_{k,X}`（或 (F.15) 的 `V_{k,B}`）；若可，单块 `diam=O(1)` 给 `C0 e^{−σ0·O(1)}` 型候选 `ε`（`C0,σ0` 的 k-无关性待核）。

## §3 能否有界：**`[NG]`**

**3.1 裁决**：Section 5 的 `sup_B sup_U‖V(B)‖` 在 `fax` 内**未被证明有界**；`22016` 是 `if` 假设；`V(B)` 无显式式 ⇒ `sup` 不可算。**数**：Section 5 的 `V(B)`：显式界 **0** 条；Section 9 的 `V_{k,X}`：**1** 条 `C0 e^{−σ0 diam(X)}`（`fax:5439`）；`ε` 数值仍 **0** 个（只 `∃ε>0`、与 `k` 无关）；`C` 数值 **0** 个（只 `C<∞`）。

**3.2 三种读法都不给界**：①未归一化 slab 自由能（`−log` 型）：块体积 `~b^{3k}`，`sup_U` 可能随 `k` 增长 ⇒ `k`-一致 `ε` 不自动成立；②归一化 slab 核（正密度；`e^{−S_W}≤1`、`R_k≤1`）：点态可能 `≤1`，但 `J^slab` 未给界，文本也未把它写成核；③抽象基泛函：可任取（如 `V=0`），但那样 `Φ_k=0`，定理无内容——非空内容需显式式。⇒ 三种都不产生 `sup_B sup_U‖V(B)‖≤ε` 的证明。

**3.3 与 `(5.32)` 的自洽性**：若 `E[V(B)]≠0`，单块 `X={B}` 的 `n=1` 项给 `Φ_k(B)=E[V(B)]+O(ε²)`；`(5.37)` 只覆盖 `O(ε²)`，故证明缺 `n=1` 项。⇒ `V(B)` 的**居中**是 `(5.32)` 的隐含前提，文本未写。

**3.4 下一手（三条，互不依赖）**：①写 `V(B;U|B)` 的 slab 积分式（含 `dμ_slab/S_W/J^slab` 与归一化），再算 `sup_B sup_U`；②若走 Section 9 的 `V_{k,X}`（或 (F.15)）类比，证 `V(B)≡V_{k,B}`（或给差异界）并核 `C0,σ0` 的 `k`-无关性，把 `C0 e^{−σ0·O(1)}` 变成候选 `ε`；③在 (5.31) 明写 `V(B)−E[V(B)]`（或 `n≥2`＋单块项另算），使 `ε²` 与 `n=1` 项自洽。

## §4 台账（`AGENTS §1` 四要素）

| # | 工具 | 状态 | 实际检索式／调用 | 命中 | 失败标记 |
|---|---|---|---|---|---|
| 1 | `read` | **OK**（带行号正文全返回，32615 行） | `fax:21960–22128`·`22330–22459`·`22690–22759`·`21740–21849`·`11890–11959`·`5425–5469`·`5495–5534`·`24770–24869`；`KKKK1-PHIDEF.md`(83 行)·`HHHH4-IDENT.md`(60 行) | 全返回 | 0 |
| 2 | `grep` | **OK**（ripgrep 正常，无 `provider failure`／`ENOENT`／空响应） | ①`V\s*\(B[^)]*\)` ②`V\s*\(B[^)]*\)\s*=` ③`base functional` ④`slab fluctuation` ⑤`operator norm` ⑥`supB\|sup_B` ⑦`ess sup` ⑧`E\s*\[\s*V\s*\(B\|V\s*\(B[^)]*\)\s*⟩\s*=\s*0` ⑨`Vk,X` ⑩`5\.21\|Theorem 5\.3\|vanishes on disconnected X\|summing over n ≥ 2\|dμslab\|J slab` | ①**6**（21997·21998·22002·22016·22022·22030）②**0** ③**1**（21997）④**1**（21997）⑤**101**（21995 定义；22369/22420 为 §6）⑥**4**（仅 22016 含 `V(B)`）⑦**1**（22361）⑧**0** ⑨**8**（5437·5439·5443·5449·11933·11941·12006·12023）⑩**25**（含 21916·21989·22015·22052·21869/21873） | 0 |
| 3 | `glob` | **OK** | `research/v3/species/LLLL1*` | **0**（建前；本文件为新建） | 0 |
| 4 | MC／Lean／`git`／`web_*`／`skill` | **0 调用** | 本靶＝仓内逐字核；catalog 无匹配技能 | — | — |

**诚实边界（4 条）**：①行号对齐本仓文本副本，未核出版页码；②(5.31) 的 OCR 折行已并，`n≥1`/`connected` 分置按 `22002–22005` 直录；③「`E[V(B)]=0`」是 `(5.32)`+`(5.37)` 的**反推**，`fax` 未明文，标 `[NG]` 不标 `[T]`；④「`sup_U`」是算子范数的解释（或 `(6.26)` 类比），`22016` 逐字只有 `sup_B`；未证之物未标 `[T]`；未改任何既有文件。

*LLLL1-VB · 只新建本文件。判：`V(B)` 无显式式（`[NG]`；`V(B)=` 0 命中）；`V` 是 (5.31) 的累积量种子、非 (5.29) 的 `Φ_k`（`[T]`）；`n=1` 项是 `E[V(B)]` 不是 `V(B)`（`[T]`）；`ε²` 反推需 `E[V(B)]=0`（`[NG]`）；`sup_B sup_U‖V(B)‖≤ε` 是 `if` 假设、无数值（`[NG]`）。*
