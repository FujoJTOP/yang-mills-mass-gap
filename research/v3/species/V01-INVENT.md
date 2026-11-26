# V01-INVENT — `INVENT-TRANSPORT-INEQ.md` 核账：**把 `β_eff<1/48` 连到 `β=2.5` 的不是传输定理，是假设 `(CP)`**

> 2026-11-25 · 第 V01 轮 · 靶：`research/unconditional/INVENT-TRANSPORT-INEQ.md`（全文 **21 行**，已逐行读完）
> 交付：本文件（仅此文件；无 git commit）
> 标签：`[T]` 已证 · `[C]` 条件 · `[E]` 精确算术 · `[NG]` no-go · `[O]` 缺失步 · `[V]` 读数
> 背景引用：`research/v3/species/U02-CP.md:34-36`（A3，逐字引了 `INVENT:17`）

---

## 0 · 判词（先给结论）

| 问 | 答 | 标签 |
|---|---|---|
| 1 · 它是定理吗 | **是「条件定理」**：`(CP) ⟹ f₁²(2.5)≥Z_O²r₀>0`。其中 `(CP)` 是**假设**，且 `INVENT` 自己判它"不是一小步" | `[C]` |
| 2 · 那"搬运"解决了吗 | **没有**。搬的是**值**，不是**原子**；`(CP)` 就是缺的那一环 | `[NG]` |
| 3 · 与 `w95d` 的 `Z_O` 是同一件事吗 | **逐字同一件事**：同一公式、`W95D-TRANSPORT.md` 自己引 `INVENT` L7/L10/L13 作定义源 | `[T]` |
| 4 · 判死线 | **主判 (b)**；(c) 对**代数部分**同时成立（`0↦0`，对"原子 vs 割"盲）；**(a) 不成立** | `[NG]` |

---

## 1 · `INVENT` 的完整论证（四段，带行号逐字）

### §A 立论：不是"下界"，是"定标协变"（L3, L5）

> `L3`：`**Verdict: the missing object is a scale-*covariance*, not a positivity *bound* — recast that way it exists, and is explicit.**`

> `L5`：`Block-RG is a change of variables on *one* theory (PROOF-TRANSPORT-RG §2), so f₁² is literally the same number under blocking; transport must relate two *couplings*.`

`L5` 是要害：`f₁²` 在**同一个理论**里是**同一个数**（塔性质），所以"传输"不能靠不变性赢；它只能把**两个耦合**用算符重整化常数连起来。⇒ 这从一开始就不是"跨相搬运"的动力学陈述。

### §B 公式：迹反常 ⇒ 乘性（L7–L13）

> `L7`：`The anomaly is exact, T^μ_μ=(β/2g)O_S, T^μ_μ the RG-invariant dilatation divergence. So (β/2g)a_G is invariant, giving`

> `L10`（boxed）：`f_1^2(\mu)=Z_O(\mu,\mu_0)^2\,f_1^2(\mu_0),\quad Z_O=(\beta/2g)(\mu_0)/(\beta/2g)(\mu)>0,\quad c:=Z_O^2=[(\beta/2g)(\mu_eff)/(\beta/2g)(\mu_2.5)]^2`

> `L13`：`And c>0 **rigorously**, being a *square* — no sign input.`

`Z_O` 的合法性论证在 `L5`：`O_S=¼trF²` 的只对角项 `a_G` 的 dim-≤4 伙伴只有 恒等（对 `G≠0` 缺席）、P-奇 `F\tilde F`（与 `0⁺⁺` 解耦）、BRST-恰当导数（物理态上为零）⇒ **乘性**。

### §C 酸测（L15，三判全 PASS）

> `L15`：`(i) *0→0:* f₁²(μ₀)=0 ⇒ Z_O²·0=0. Multiplicative, so it cannot manufacture weight from zero … (ii) *Free field:* O=:φ²: has a(β)≡0 by Z₂ parity … f₁²≡0 ↦ Z_O²·0=0. **PASS.** (iii) c=Z_O²>0 explicit, V-uniform … **PASS.**`

**注意这一条就是 `w95d` 判空洞的地方**：`INVENT` 把 `0→0` 写成"通过"；同一性质换一个读法就是"对原子 vs 割盲"（见 §5）。

### §D 条件定理 + 自己点出的缺口（L17）

> `L17`：`*THEOREM.* Fix N, T=0. Assume **(CP)**: β_eff<1/48 and β=2.5 lie on one constant-physics line — a single confining phase, lightest 0⁺⁺ Kato-isolated (Λ(β)−m₁²≥κ>0) and carrying the edge throughout (no bulk transition). Then f₁²(2.5)=Z_O²f₁²(β_eff)≥Z_O²r₀(β_eff)>0 (P3 base point). Transport is *exact*, the constant *explicit*. *The gap:* (CP) is the continuum-existence / single-T=0-phase statement; through the ultra-coarse region (a·m₁≫1) the strong-coupling "glueball" is not yet a physical particle, so (CP) is no short step — it is H*(r₀)+continuity in minimal dress.`

> `L19`（Sharpest line）：`What remains is not the residue but the *phase*: that strong- and weak-coupling lattices be one T=0 theory.`

⇒ `INVENT` 自己把全部非平凡内容放在 `(CP)`，并明说它 `no short step`。**没有任何一行声称 `(CP)` 已证。**

---

## 2 · "one constant-physics line" 的精确定义

**不是**一句话口号，是 `(CP)` 这个三元合取（`INVENT:17` 逐字）：

```
(i)   单 T=0 禁闭相，无体相变                     single confining phase, no bulk transition
(ii)  最轻 0⁺⁺ 全窗 Kato-隔离                      Λ(β)−m₁² ≥ κ > 0  （κ>0，∀β∈[β_eff,2.5]）
(iii) 该边一路被携带（一致隔离，不溶解）           carrying the edge throughout
+ 前提：β_eff<1/48 与 β=2.5 属于同一条线（同一理论、同一个态）
```

**等价写法**（把 (iii) 量化）：`R1-RESIDUE-TRANSPORT.md:252-256`

> `Assume **(CP)**: (i) β_eff<1/48 and the physical point β_tgt lie on one constant-physics line (a single confining phase), and (ii) the lightest 0⁺⁺ atom is Kato-isolated and carried throughout, i.e. Ω(β_tgt,β_eff)=1 ⇔ uniform Λ_β−m₁²(β)≥δ>0 on [β_eff,β_tgt].`

⇒ **"one constant-physics line" ≡ `Ω=1`**（态重叠为 1：强耦合"胶球"就是物理胶球，无混合、跨界不溶解）。这句是**定义性等价**，不是定理：

> `LOOP-RUNBOOK.md:6341`：`Ω 大 ⟺ 同一个原子 + 一致隔离 ⟺ I2 ⟺ (CP) ⟺ 隙      [O]`

---

## 3 · 等级（`AGENTS §3`）

| 层 | 内容 | 等级 | 依据 |
|---|---|---|---|
| 代数 | `f₁²(μ)=Z_O²f₁²(μ₀)`，`Z_O=(β/2g)` 之比、`c=Z_O²>0` | **`[T]`** | `INVENT:10,13`；`R1:279-287`；`W95D:15,65` |
| 定理 | `(CP) ⟹ f₁²(2.5)≥Z_O²r₀>0` | **`[C]`**（given (CP)） | `INVENT:17`；`R1:252`「`[T]` given `(CP)`」 |
| `(CP)` 本身 | 单相 + Kato 隔离 + 一致携带 | **`[O]` 未证** | `R1:264`；`U02-CP:135` |
| "搬运已解决" | — | **`[NG]`** | `W95D:16,18,21`；`R1:243-246,266-268` |

三条独立判词，全部指向同一处：

> `R1-RESIDUE-TRANSPORT.md:264`：`**This is a [T]-given-(CP) statement; (CP) is NOT proved and is exactly the gap.**`
> `U02-CP.md:135`：`| (CP) 在 β<1/48 是否已证 | **否**：只有 I1，I2 全开 |`
> `U02-CP.md:102`：`**I2 未证** ⇒ **(CP) 在 β<1/48 也不成立**（不是"窗口内已证"）。`

⇒ 连**强耦合窗口内**都没有 `(CP)`：`P3` 只给 `I1`（残数地板），`I2`（隔离）全开。

---

## 4 · 若按 `INVENT` 的读法"它用了什么" —— 用了什么、没用什么

**用了**（`L5, L7`）：

1. **迹反常 Ward 恒等式**：`T^μ_μ=(β/2g)O_S` 精确 ⇒ `(β/2g)a_G` RG-不变 ⇒ `a_G ∝ (β/2g)^{-1}`。
2. **dim-≤4 算符混合唯一性**（Collins, *Renormalization*）：伙伴只有 恒等 / `F\tilde F` / BRST-恰当 ⇒ 乘性。
3. **`Z₂` 宇称**（酸测 (ii)：`:φ²:` 的 `a≡0`）。
4. `P3` 的强耦合地板 `r₀(β_eff)`（作为**基点**，不搬运）。

**没有用**（这是关键）：**单调性没有**、**解析延拓没有**、**RG 流没有证**。所以它**不是一台搬运机器**，而是一个**归一化协变式**。

> `R1:232-241`（N3）：`Equation (5.1) **is** the source of the multiplicative law … It fixes the *scaling* and nothing else — it contains f_G, a free matrix element. There is **no** Ward identity that lower-bounds a_G … the Ward identity is Z_O, and Z_O is free.`

### `β_eff` 是什么

`β_eff(b,β)` = **一次因子 `b` 块化后的粗格点等效耦合**，定义在 `PROOF-TRANSPORT-RG.md:7`：

> `β_eff(b,β) by matching the renormalized 2βN_p (Migdal–Kadanoff at leading order).`

**而它不是免费取到的 `<1/48`**：`PROOF-TRANSPORT-RG.md:16`

> `**One** factor-2 step gives β_eff≈β (slow UV running; PROOF-H2-GLOBAL route (i)); reaching g²~10² spans ~50 e-foldings ⇒ **tens** of steps, *through the crossover*.`

⇒ `(CP)` 的**前提**"`β_eff<1/48` 与 `β=2.5` 同线"已经**预设了迭代流穿过交叉点存在** —— 那正是 `W2` 墙。所以 `(CP)` 里其实压着**两道**墙：**(流存在) + (相/原子同一)**，`INVENT` 只写了后一道（`L19`"not the residue but the *phase*"），前一道在 `PROOF-TRANSPORT-RG` 里是明文 `[NG]`。

---

## 5 · 与 `w95d` 的 `Z_O` 是同一件事吗 —— **是，逐字同一件事**

| 项 | `INVENT` | `W95D-TRANSPORT.md` |
|---|---|---|
| 公式 | `f₁²(μ)=Z_O²f₁²(μ₀)`（`L10`） | `f₁²(μ)=Z_O²f₁²(μ₀)`（`W95D:3,90`） |
| `Z_O` 定义 | `(β/2g)(μ₀)/(β/2g)(μ)`（`L10`） | `Z_O = (β_YM/2g)(μ₀)/(β_YM/2g)(μ)`（`W95D:65`） |
| 谁引谁 | — | `W95D:169`：`定义：… research/unconditional/INVENT-TRANSPORT-INEQ.md L7, L10, L13` |

`W95D` 自己把语料的乘性搬运判成 **`Ω=1` 的特例**：

> `W95D:92`：`语料的 f₁²=Z_O²f₁²  =  特例 Ω=1`
> `W95D:21`：`f₁²(μ)=Z_O²f₁²(μ₀) 是**同一连续理论、同一个态**内的精确定标协变 [T]，跨 β 的忠实形式是 f₁²(2.5)=Z_O²·Ω·f₁²(1/48)，而 Ω ≤ 1 只有上界——**搬运搬的是*值*，不是*原子***。`

⇒ **不是"两件事碰巧同名"**：`INVENT` 的传输**就是** `w95d` 判空洞的那条 `Z_O` 传输，且 `w95d` 的两条判词各自独立地成立：

> `W95D:16`：`符号免费，"结论"不免费。等价式**不生成信息** … 0↦0 过酸测正是它对"原子 vs 割"盲的证明；有限 (L,a) 上**每个**残数都 >0 … 搬过去的"正"是空洞正` `[T]/[NG]`
> `W95D:18`：`**(b) 成立**（缺 Ω）；主判 = (b)`

**唯一的差别（且不救场）**：`INVENT` 把 `(CP)` 明确写在**假设位**（`L17` `Assume **(CP)**`），`w95d` 把同一个东西叫 `Ω=1`。**两者是同一个缺口的两个名字**，`Ω` 只有上界（Cauchy–Schwarz）：

> `R1:226-230`：`**(N2) Ω is bounded from the wrong side.** … an **upper** bound on the overlap. Every ν≥0/RP functional is Schur-concave … So Ω≥c>0 is **not** in the positivity/RP/spectral class` `[NG]`

---

## 6 · 判死线

- **(a) 不成立。** `INVENT` **没有**已证的传输定理。它有的是：`[T]` 的代数恒等式 + `[C]` 的条件定理，条件 `(CP)` 未被证，且 `INVENT:17` 自己写 `(CP) is no short step`。⇒ **本轮不产生"搬运问题已解决"这一发现**；`U02-CP:34-36` 引 `INVENT:17` 时引的是**假设的措辞**（`U02-CP:36` 逐字带 `Assume **(CP)**`），不是已证命题 —— 引用本身没错位，但**不能当作已证**。
- **(b) 主判成立。** 缺的东西：`Ω≥ω₀>0` 一致（等价 `I2`、等价 `(CP)`、等价隙）；前置还缺"迭代块化流全程存在"（`PROOF-TRANSPORT-RG:16` 的 `[NG]`）。`Z_O` 的存在性/符号/数值**都不是**缺的东西（`W95D:17` 逐字："**不需要**：`γ_O`…、a-theorem…、Ward 恒等式（**它就是** `Z_O`）"）。
- **(c) 对代数部分同时成立。** `INVENT:15` 的酸测 (i)(ii) 与 `W95D:16` 的"空洞"判词是**同一事实的两种读法**：乘性同态把 `0` 送到 `0`，故它对"原子 vs 割"**结构性盲**（`R1:284-285`：`Z_O² multiplies every spectral weight alike, so the *shape* (atom vs cut) is preserved, not created`）。⇒ `f₁²(2.5)>0` 一旦由这条搬得，就是一个在 2D、单 plaquette 模型里也成立的空洞命题（`R1:286-287`）。

---

## 7 · 工具状态（`AGENTS.md §1`）

| 工具 | 状态 | 检索式 / 对象 | 命中 |
|---|---|---|---|
| `glob` | **OK** | `**/INVENT-TRANSPORT-INEQ.md`、`**/U02-CP.md`、`**/*w95d*`、`**/PROOF-TRANSPORT-RG.md`、`**/W95D-TRANSPORT.md` | 前二各 1；`*w95d*` 文件名 **0**（真名 `W95D-TRANSPORT.md`，后经 `grep w95d` 定位） |
| `read` | **OK** | `INVENT-TRANSPORT-INEQ.md` 全文 21 行；`PROOF-TRANSPORT-RG.md` 全文 31 行；`W95D-TRANSPORT.md` 全文 178 行；`U02-CP.md` L1–80；`R1-RESIDUE-TRANSPORT.md` L225–294；`LOOP-RUNBOOK.md` L6320–6409 | 全部有效返回 |
| `grep` | **OK** | `w95d`（全仓） | **2**（`LOOP-RUNBOOK.md:6360,6361`） |
| `grep` | **OK** | `constant-physics line\|constant physics line` | **3**（`R1:254`、`U02-CP:36`、`INVENT:17`）—— 全部在**假设位**，无一条在证明位 |
| `grep` | **OK** | `\(CP\).{0,40}(已证\|proved\|proven\|定理成立)` 等 | **4**，全部指向**未证**（`R1:264`、`U02-CP:102,135`、`MS-ANALYTIC-NOZERO:14`） |

**无 TOOL-FAIL。** 未做 `git commit`。只写本文件。

---

## 8 · 一句话

`INVENT` 的传输是**真而空洞**的定标协变（`[T]`），它被写成一个**条件定理**（`[C]`），条件 `(CP)` 与 `w95d` 的 `Ω=1` **逐字同一**、且语料自认未证（`R1:264`、`U02-CP:135`）——**搬运问题没有解决，缺的是 `Ω`/`(CP)`，不是 `Z_O`。**

*状态：`(a)` 排除、`(b)` 主判、`(c)` 对代数部分成立。Clay 0%，无 QED。*
