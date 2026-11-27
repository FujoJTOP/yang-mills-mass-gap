# NNN1-OS1 — OS1（平移过弱极限）：闭成 `[T]`，但"免费"边界比原判窄一格

> 轮 NNN1 · 子代理 `nnn1-os1` · `D:\Dev\math`。写域仅本文件；既有文件一字未改；未 `commit`；零新 MC；**零 Lean 编译**。进展单位 = 轮。
> 口径：只 `[T]`/`[NG]`；无百分数；每结论带「下一手」。`AGENTS §1` 四要素见 §4。

## §0 判词

| # | 断言 | 判 | 依据 | 下一手 |
|---|---|---|---|---|
| 0.1 | **`ggg3-osclosed [T]` 成立**：`μ_a ⇒ μ` 且 `μ_a` 在 `aℤ⁴` 下不变 ⇒ `μ` 在 `ℝ⁴` 下不变。一般形式只要求「Polish 加性拓扑群 + 平移等距」 | **`[T]`** | §1 S1–S5 | 引它时把三前提（D1–D3）一起抄 |
| 0.2 | **题给的"依据"须改一格**：`f` 在 `𝒮′` 这类**非局部紧**空间上**没有**整体一致连续性，故"`f` 一致连续 ⇒ 收尾"不成立。真承重件是**紧集上一致连续（Heine–Cantor）+ 平移等距 + tightness** | **`[T]`**（更正） | §1 S3；`GGG3-OSCLOSED.md:18` 缺 `tight` 与"紧上" | 台账把"一致连续"改成"紧上一致连续 + tight" |
| 0.3 | **不需要精确不变**：只要求缺 `ε_a→0`（`\|⟨f(·+h_a)⟩_{μ_a}−⟨f⟩_{μ_a}\|≤ε_a‖f‖_∞`），结论仍 `[T]` | **`[T]`** | §1 S2 | RT/有限体积的 `o(1)` 破缺可直接套 |
| 0.4 | **仓内被引作"OS1 已证"的 `os2 §3.5 L4`（`os2_of_dense_subgroup`）走 Z-连续路，`hZ : Continuous Z` = OS0** ⇒ 那条路**不是零代价**；零代价的是 §1 的**测度路** | **`[T]`**（分账） | `os2-rotation/REPORT.md:337–357` | 引 OS1 必标走哪条路；禁把 Z 路记成"免费" |
| 0.5 | **酸测不过**：无质量自由场（Gaussian）同样满足 OS1 ⇒ OS1 不区分 YM 与自由场，**零物理内容**，是**闭条件**（只筛不造） | **`[T]`** | `GGG3-OSCLOSED.md:40`；§1 末 | 位置图记「闭条件、零代价、零内容」 |

## §1 完整证明（一般形式，五件套逐件点名）

**设定**：`E` 加性 Polish 拓扑群（`(x,y)↦x+y` 连续），度量 `d` **平移不变** `d(x+z,y+z)=d(x,y)`；`τ_t x:=x+t`。`μ_a,μ∈𝒫(E)`，`μ_a ⇒ μ`（弱收敛）。`G_a:=aℤ⁴`；**假设** `(τ_h)_*μ_a=μ_a`，`∀h∈G_a`。取定 `t∈ℝ⁴`、`f∈C_b(E)`。

**S1（稠密，显式模）** `h_a:=a⌊t/a⌋∈G_a`，逐坐标 `\|h_{a,i}−t_i\|≤a` ⇒ `\|h_a−t\|≤2a→0`。**密度只用在这里**：`cl(∪_a aℤ⁴)=ℝ⁴` 保证「每个 `t` 都有一条 `G_a∋h_a→t` 的逼近」，且给显式速率 `2a`。
**S2（不变性搬过来）** `⟨f⟩_{μ_a}=∫f(x+h_a)dμ_a(x)`，`∀a`；破缺版为 `\|⟨f⟩_{μ_a}−∫f(x+h_a)dμ_a\|≤ε_a‖f‖_∞`（⇒ 0.3）。
**S3（tightness + 平移等距 + Heine–Cantor）** 由 `E` Polish 与 `μ_a⇒μ` 得 `{μ_a}` tight（Prokhorov；应用里可由 `UVS` 直供）：`∀ε ∃K` 紧, `μ_a(K^c)≤ε ∀a`。令 `ω_a:=\sup_{x∈K}\|f(x+h_a)−f(x+t)\|`。`S:={h_a}∪{t}` 紧（收敛列 + 极限）；`F(x,y):=f(x+y)` 连续 ⇒ 在紧集 `K×S` 上**一致连续**；又因 `d` 平移不变，`d(x+h_a,x+t)=d(h_a,t)≤2a` 与 `x` **无关** ⇒ `ω_a≤η(2a)→0`。
**S4（拼）** `\|⟨f⟩_{μ_a}−⟨f⟩_{(τ_t)_*μ_a}\| ≤ ω_a+2‖f‖_∞ε`；另 `∫f(x+t)dμ_a=∫f d(τ_t)_*μ_a → ∫f d(τ_t)_*μ`（`f∈C_b` + `τ_t` 连续），`⟨f⟩_{μ_a}→⟨f⟩_μ`。故 `\|⟨f⟩_μ−⟨f⟩_{(τ_t)_*μ}\|≤2‖f‖_∞ε` 对一切 `ε` ⇒ `=0`。
**S5（分离）** `f∈C_b` 任意，`C_b` 分离 `𝒫(E)`（`E` Polish）⇒ `(τ_t)_*μ=μ`；`t` 任意 ⇒ `μ` 在 `ℝ⁴` 下不变。∎

**五件套各在哪一步**：密度 → S1；平移等距 → S3（让模与 `x` 无关）；紧上一致连续 → S3（`ω_a→0`）；tightness → S3（`K^c` 余项）；弱收敛 → S4（两端）。
**与 OS3 的分界（同一模板）**：此处 `cl(∪_a aℤ⁴)=ℝ⁴`（群随 `a` 变密）；`O(4)` 侧 `cl(∪_a W)=W⊊O(4)`（`\|W\|=384` 不随 `a` 长）⇒ 模板给不出 `O(4)`。这是 OS1 `[T]` / OS3 `[NG]` 的唯一差别。

## §2 Lean 化判

| 需要件 | Mathlib | 出处（本机 `vendor/OSforGFF/.lake/packages/mathlib`） |
|---|---|---|
| 弱收敛 ⇔ `∀f∈C_b` 积分收敛 | **有** | `MeasureTheory/Measure/ProbabilityMeasure.lean:346`；`FiniteMeasure.lean:726` |
| tight 定义 + 紧集列等价 + 相对紧 | **有** | `Measure/Tight.lean:55,61`；`Measure/Prokhorov.lean:503`（`isCompact_closure_of_isTightMeasureSet`） |
| 紧集上一致连续 | **有** | `Topology/UniformSpace/HeineCantor.lean:47`（`IsCompact.uniformContinuousOn_of_continuous`） |
| 平移是同胚（连续+可测白拿） | **有** | `Topology/Algebra/Group/Basic.lean`（`Homeomorph.addLeft`；用例 `Topology/Instances/CantorSet.lean:121`） |
| 推移复合 `map_map` | **有** | `Measure/Map.lean:203` |
| BL 版弱收敛判据（可对 Lipschitz 常数取界） | **有** | `Portmanteau.lean:688`（`tendsto_iff_forall_lipschitz_integral_tendsto`） |
| **`μ_a⇒μ ⇒ {μ_a}` tight（一般 Polish）** | **缺** | 只有 charFun（有限维，`LevyConvergence.lean:45`）与 normed/proper 版（`TightNormed.lean:57,85`）；一般 Polish **0 命中** |
| **联合推移连续 `(μ_n⇒μ,y_n→y)⇒μ_n.map(·+y_n)⇒μ.map(·+y)`** | **缺** | 全 Mathlib **0 命中**（= S3–S4 的变动被积函数步） |
| `𝒮′(ℝ⁴)` 的 TVS/Polish/Borel 实例 + 平移作用 | **缺** | 只有测试函数侧 `SchwartzMap.compSubConstCLM`（`Analysis/Distribution/SchwartzSpace/Basic.lean:1068`） |

判：**抽象骨架可 Lean**（三条显式假设 `hTendsto`/`hTight`/`hInv`；缺的两件自证 ≈40 + ≈100 行，或把 `hTight` 当假设——照 `glue-a` 的 `hMinlos` 体例）。**`𝒮′` 实体不可 Lean**（实例与格点测度都不在 Mathlib）。故 `[T]` 是**数学**的，Lean 只到骨架。

## §3 隐藏依赖与裁决

| # | 前提 | 是 OS0 吗 | 判 |
|---|---|---|---|
| D1 | **公共载体 + `ℝ⁴` 作用连续**：`μ_a` 与 `μ` 必须同在一个 `ℝ⁴`-空间上。若 `μ_a` 只活在 `𝕋⁴_{L_a}`，`ℝ⁴` 平移**没有对象**，OS1 空转（要 `L_a→∞` 与 `a→0` 同步、且嵌入与平移相容） | 不是：**建模前提** | `[T]`（前提性） |
| D2 | **极限已存在**：`μ_a⇒μ`（= GL3，消费 `UVS`/`L3`） | **是（其上游）** | OS1 **携带不了存在性**；"免费" = 条件于极限在手 |
| D3 | **`aℤ⁴` 下精确不变**（Wilson + 周期 bc） | 不是 | 可放宽到 `ε_a→0`（0.3）仍 `[T]` |
| D4 | **路选择**：`os2 §3.5 L4`（Z 连续） | **是**：`hZ` = OS0 的解析性面（`os2-rotation/REPORT.md:355`） | 该路**吃 OS0**；测度路不吃 ⇒ 二路分账 |

**裁决**：`ggg3-osclosed [T]` 成立，且**不依赖 OS0**（一般 Polish 形式；D2 除外）。"免费"的准确内容 = 「**给定极限存在 + 共同载体，平移不变性零代价**」——不含存在性、不含温性、不含 `O(4)`、**零物理内容**。下一手：`hTight` 按 `hMinlos` 体例写成显式假设；实体 Lean 化等 `𝒮′` 实例。位置图变化（OS1 由"声称"→"证明 + 边界收窄"）按 `AGENTS §6` 由 `lead` 回填（受"不改既有文件"约束）。

## §4 检索核验（`AGENTS §1`）

| 工具状态 | 检索式 | 命中 |
|---|---|---|
| `grep` **OK**（有效响应，无 `provider failure`/`ENOENT`/空响应）；`read` **OK**（4 文件带行号）；`glob` **OK** | ①`tendsto_iff_forall_integral_tendsto\|forall_integral_tendsto` @ mathlib ②`class IsTightMeasureSet\|IsTight\|Prokhorov` @ mathlib/MeasureTheory ③`addLeft` @ mathlib/Topology ④`uniformContinuousOn_of_compact\|IsCompact.add\|IsCompact.smul` @ mathlib ⑤`map_add_left\|map_sub_left` @ mathlib/MeasureTheory/Measure ⑥`CompactOpen\|locallyUniform\|tendsto_integral` @ mathlib/MeasureTheory ⑦`IsTightMeasureSet.*of_tendsto` @ mathlib ⑧`os2_of_dense_subgroup` @ research ⑨`OS0\|OS1\|平移不变` @ research/v3/*.md | ①**23** ②**225** ③**4** ④**3**（全 `smul`；`IsCompact.add` **0**）⑤**0** ⑥**74** ⑦**10**（全 charFun/normed/proper；一般 Polish **0**）⑧**3** ⑨**257**（显示 250，全表落盘） |

**失败标记**：`TOOL-FAIL` **0**。上表"缺"全部为**已执行检索 0 命中**，非工具失败。**本轮零 Lean 编译**，§2 全为静态读源码。
