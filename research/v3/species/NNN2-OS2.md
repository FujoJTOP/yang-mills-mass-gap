# NNN2-OS2 — OS2(RP) 弱闭：`C_b` 层闭性是 `[T]`，全锥「已 Lean 化」只到有界层

> 靶 `GGG3-OSCLOSED.md:19`（`ggg3-osclosed [T]`）：①`C_RP` 是弱闭凸锥、逐项非负取极限仍非负；②不需要一致矩界；③已 Lean 化（GL4）。
> 编号：本文件按**靶**编号（OS2=RP、OS3=对称，`GGG3:4`）；仓内 `docs/conditional-theorem.md:45`、`G1-CONTINUUM-OS.md:5` 用 OS2=`O(4)`/OS3=RP，跨文件必带号。
> 纪律：只新建本文件；不改既有文件、不 commit、零新 MC、不编译 Lean；只交 `[T]`/`[NG]`。

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| 0.1 | `C_RP:={μ:∀A∈C_b, ∫θA·A dμ≥0}` 弱闭凸锥，`μ_a∈C_RP, μ_a⇀μ ⇒ μ∈C_RP`，`C_b` 层 | **`[T]`** | §1(i)–(iii)；`REPORT.md:161–163`；Lean `GL4Transfer.lean:48–66` | 不再立项；标点修一处（概率测度集不是锥，见 §1(i)） |
| 0.2 | 「不需要一致矩界」逐字成立（`C_b` 层与全锥层都成立） | **`[T]`** | §2；`REPORT.md:150,173`；`G1-EXISTENCE-THEOREM.md:18` | 引用时把「不需要矩」写成「不需要 a-一致矩；**要**极限侧 `A∈L²(μ)`」 |
| 0.3 | 「OS2 已 Lean 化（GL4）」端到端（全锥 + 不变性） | **`[NG]`** | §3：缺 `θ_*μ=μ`、缺截断族的 μ-侧 RP 交接、`trunc` 与 `Ω→ᵇℝ` 不兼容；`lean/` 内 RP = 0 | 补 3 件中的任一件即升档（§3 末） |

**一句话。** 靶成立但要把「已 Lean 化」降口径：**有界连续层**的弱极限传递是真机检（Mathlib 真弱收敛），**全锥**的 Lean 是「一条抽象 L² 闭包引理 + 三条接过来的假设」，其交接处与 `REPORT.md:164` 的散文接缝是**同一处**。

## §1 证明（完整；所用三项逐条点名）

设 `Ω` 紧度量、`θ:Ω→Ω` 连续对合、`θ_*μ_a=μ_a`、`μ_a` 概率、`μ_a⇀μ`（弱＝对一切 `f∈C_b`：`∫f dμ_a→∫f dμ`）、每个 `μ_a` 满足 RP。

1. **凸锥**：每个 `A`，`Λ_A(μ):=∫\overline{θA}A dμ` 关于 `μ` 线性；`C_RP=⋂_{A∈C_b}{Λ_A≥0}` ＝半空间之交 ⇒ 凸锥。**标点修**：概率测度集对伸缩不闭，锥性只在 `M₊(Ω)` 上成立；`Prob` 上只剩「弱闭凸」，`GGG3:19` 的「弱闭凸锥」按此读。
2. **弱闭**：`Λ_A` 弱连续，**充要输入**是 `f_A:=\overline{θA}A∈C_b`——由 `A∈C_b`、`θ` 连续（`A∘θ` 连续，乘积有界连续）。`μ_k⇀μ ⇒ Λ_A(μ_k)→Λ_A(μ)`，`[0,∞)` 闭 ⇒ `Λ_A(μ)≥0`，`A` 任意 ⇒ `μ∈C_RP`。
3. **传递**：取 `μ_a∈C_RP`、`μ_a⇀μ`、固定 `A`：`Λ_A(μ_a)≥0 ∀a` 且 `Λ_A(μ_a)→Λ_A(μ)` ⇒ `Λ_A(μ)≥0`。∎（全程无 `a`-一致量）
4. **`θ_*μ=μ`（RP 的另一半，同样免费）**：`f∈C_b ⇒ f∘θ∈C_b`；`∫f dθ_*μ=∫f∘θ dμ=lim_a∫f∘θ dμ_a=lim_a∫f dμ_a=∫f dμ`。用 `θ` 连续 + 格点不变性。
5. **全锥（`A∈L²(μ)`）**：`A_R:=A·1_{|A|≤R}`，`‖A_R−A‖₂→0`（控制收敛，控制 `4|A|²∈L¹`）；
 `|∫(θA·A−θA_R·A_R)dμ| ≤ ‖θA‖₂‖A−A_R‖₂+‖A_R‖₂‖θA−θA_R‖₂ ≤ 2‖A‖₂‖A−A_R‖₂ → 0`，
 其中 `‖θA‖₂=‖A‖₂` 只由 `θ_*μ=μ`。⇒ 不等式过极限。∎

**用了什么**：① `θ` 连续 ⇒ `θA·A∈C_b`（这是整条的**唯一**实质结构条件，也是 `GGG3:26` 留的 F4-4 挂点：反射面/分块不对齐则 `θ_*μ_a=μ_a` 先坏）；② 弱收敛的**定义**（对一切 `C_b` 积分收敛；Lean :62 用 Mathlib `ProbabilityMeasure.tendsto_iff_forall_integral_tendsto`，非抽象 `hconv`）；③ `[0,∞)` 闭。
**没用**：tightness/Prokhorov（极限按假设已存在）、`sup_a` 的任意阶矩界、`L^p` 一致界、Minlos。

## §2 「不需要一致矩界」——逐字：对

- **为什么对**：步 3 的通篇只有「逐项 `≥0`」与「`Λ_A(μ_a)→Λ_A(μ)`」；弱收敛的定义已含「所有 `C_b` 测试的积分收敛」，故不需要 `sup_a ∫|A|^p dμ_a<∞` 之类的量来支撑取极限。全锥层同理：截断只在**极限测度**里做，控制函数是 `μ` 的 `4|A|²`，与 `a` 无关（`REPORT.md:173`）。
- **边界 1（对象在句外）**：`μ` 的**存在**要 tightness＝`UVS`/`L3`（`GGG3:17`；`CONTINUUM-TIGHTNESS.md:3`）。本句是「若有极限则传递」，不是「有极限」。
- **边界 2（无界层要极限侧 L²）**：「不需要矩界」＝不需要 **a-一致**矩界，≠ 不需要可积性；要 `A∈L²(μ)`。`REPORT.md:238` 的 F4-3（`μ_n=(1−1/n)δ₀+(1/n)δ_{n²}`、`A(x)=x`）正是把「无 L² 控制就直接对无界 `A` 取极限」打掉的反例。
- **边界 3（不变性另需输入）**：`θ_*μ=μ` 免费但非无输入——要 `θ` 连续 + `θ_*μ_a=μ_a`（F4-4，`REPORT.md:239`）。
- **分水岭成立**：OS0 的 tightness 与 OS4 的 `δ` 都要 a-一致量，RP 不要（`GGG3:11` 的一行判据）。

## §3 「已 Lean 化（GL4）」的实际范围

| `file:line` | 内容 | 覆盖 OS2 的哪一段 | 判 |
|---|---|---|---|
| `research/glue-a-limits/GL4Transfer.lean:48–66` | `rp_transfer`：`Continuous θ`(:49) + `ProbabilityMeasure` 弱拓扑真弱收敛(:52) + `A:Ω→ᵇℝ`(:53) + 逐项 `hRP`(:54) ⇒ 极限 `∫θA·A≥0`(:55) | §1 步 2–3 的 **`C_b` 层** | **`[T]`**（本轮不编译，采信依据＝陈述与证明项，非 `:97` 的 `#print axioms`） |
| `GL4L2.lean:86–93` | `rp_l2_closure`：抽象 L² 闭包 | §1 步 5 的「闭」 | **`[T]` 作为引理**：假设 `MeasurePreserving θ μ μ`(:87)、`hRP` 逐 `n` 在 **μ** 上(:92) |
| `GL4L2.lean:210–217` + `:28–29` | `rp_l2_closure_trunc`，`trunc A n := if |A ω|≤n then A ω else 0` | 同上，截断族专版 | **`[NG]` 作为端到端**：`trunc` **不连续**，且 `:25` 变量段**无 `TopologicalSpace`** ⇒ 喂不进 `rp_transfer`（要 `Ω→ᵇℝ`） |
| `θ_*μ=μ` | `grep "Measure.map\|map_eq\|map_apply\|MeasurePreserving\|eq_map" @ research/glue-a-limits/*.lean` → **3** 命中，全是 `MeasurePreserving` **假设**（`GL4L2.lean:74,87,211`）；定理 **0** | 无 | **`[NG]`**：`REPORT.md:162` 的三行散文无 Lean |
| `lean/` 交付树 | `grep "reflect\|Reflect\|rp_\|Positive\|positivity" @ lean/*.lean` → **3**，全为 `positivity` 策略（`Reduction.lean:254`、`Conditional.lean:74,95`）；RP **0** | 无 | 「Lean 化」指 `research/glue-a-limits/`，不是 `lean/` |

**裁决。** ① `rp_transfer` 真覆盖 `C_b` 层（含 `θ` 连续这个唯一结构条件）；② 全锥的 Lean 是「闭包引理 + 交接假设」，**不是**「从格点 RP 一路到极限全锥」的单一机检；③ `REPORT.md:215`／`GGG2-COMPACT.md:29` 的「唯一显式假设＝`trunc A n∈L²`」应改读为「**5 条**显式假设（`θ:Ω≃ᵐΩ`、`hθ:MeasurePreserving`、`hA:MemLp A 2`、`hAR:∀n MemLp(trunc A n) 2`、`hRP:∀n` 不等式；其中 `hAR` 可由 `hA` 导出），**载重的是 `hθ` 与 `hRP` 两条**」；④ 散文接缝与 Lean 接缝是同一处：`trunc` 不连续 ⇒ 「有界层 ⇒ 截断层」的交接在两侧都靠「支撑/代数闭包在截断下保持，故 `A_R∈C_RP`」这句，需要 `A_R` 连续或走 Lusin 密度近似，**该步未写出**。
**下一手（任一即升档）**：㈠ 把 `trunc` 换成 `clip_n∘A = max(−n,min(n,A))`（连续、`|clip|≤|A|`、逐点收敛），并给 `Ω` 加拓扑，使 `hRP` 由 `rp_transfer` 逐 `n` 供给 ⇒ 三定理并成一条；㈡ 补 `θ_*μ=μ` 的 Lean 引理（`Measure.map` 一行级）；㈢ 把「`hRP` 在 μ 上」改写成「`C_RP∩L²(μ)` 由稠密子类逼近」的显式假设，并在 `REPORT.md:164` 写下该密度步。

## §4 检索核验（AGENTS §1）

| 项 | 值 |
|---|---|
| 工具状态 | `pwsh` OK · `glob` OK · `grep` OK · `read` OK（无 `TOOL-FAIL`） |
| 检索式 | `GL4` @ 全仓；`reflect\|Reflect\|rp_\|Positive\|positivity` @ `lean/*.lean`；`Measure.map\|map_eq\|map_apply\|MeasurePreserving\|eq_map` @ `research/glue-a-limits/*.lean`；`research/v3/species/NNN*` |
| 命中 | `GL4` **250/286**（落 spill）；`lean/` 内 RP **0**（`positivity` 策略 3）；`Measure.map` **0**、`MeasurePreserving` **3**（均假设行）；`NNN*` **0** |
| 失败标记 | 无 |
| 位置图 | 本轮结论变化（0.3 由 `[T]` 降为 `[NG]`：全锥端到端）触发 AGENTS §6；受「不改既有文件」约束，`docs/position-current.html` 由 lead 回填 |
