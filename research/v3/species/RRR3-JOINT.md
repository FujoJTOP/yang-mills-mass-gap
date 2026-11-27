# RRR3-JOINT — `hJoint` 的 ε-论证（OS1 唯一剩余 `[NG]` 的收口）

> 2026-11-27 · 轮次 RRR3 · 靶：`OOO1-LEAN-OS1` §1.②b「测度层移动试验函数」的 ~80 行。
> 产出：`lean/os1/OS1JointShift.lean`（②a 保留，②b 写成）+ 本文件。
> **不编译**（RRR3 令）。全部 Mathlib 名字本轮**静态读过源码**，标 `文件:行`。

---

## §0 判词

| # | 判 | 等级 | 依据 | 下一手 |
|---|---|---|---|---|
| 0.1 | **②b 的数学内容已完整落纸，材料零缺** | **`[T]`** | §1 的 ε-链（8 步、无洞）；§3 依赖清单 12 件全在 Mathlib（逐条 `文件:行`） | 编译轮把 3 处 `[GLUE]` 收成 0 `sorry` |
| 0.2 | **「~80 行」是低估**：实测 **≈100–120 行**（含 statement 与标注） | **`[T]`**（计数） | 缺口不在数学，在**记账**：`integral_add_compl` + 两次 `norm_setIntegral_le_of_norm_le_const` + `ENNReal→real` 3 次转换，占 ~45 行 | 报「80 行」的地方改「~110 行」 |
| 0.3 | **②b 不需要 E 上的度量**（连 `PseudoMetricSpace` 都不要） | **`[T]`** | 走 `IsCompact.uniformContinuousAt_of_continuousAt`（`HeineCantor.lean:58`），它直接在**抽象一致结构**上给点对 closeness；`Metric.uniformContinuousOn_iff`（`Pseudo/Defs.lean:762`）那条**要** E 伪度量，本轮**不用**它 | 台账把 OOO1 的「不需要平移不变度量」加强为「不需要任何度量」 |
| 0.4 | **0 `sorry` 可达** | **`[C]`**（未编译，条件于编译轮） | 3 处 `sorry` 全标 `[GLUE]`，**均非缺件**：①`𝓤(E×E)` 拆两个 `𝓤E`；②`μ.real K ≤ 1`；③纯算术 | 下一轮开 `lean/` 的 lakefile 后编译 |

**一句话**：②b 从来不是「缺定理」，是「缺一次把 ε 摆开的写作」。写完后，OS 链上 `[NG]` 只剩 `𝒮′` 实体本身，`hJoint` 可以从 `OS1Statement.lean` 的假设降为定理。

---

## §1 ε-分解（8 步，完整不等式链）

记 `Aₙ = ∫ f(x+yₙ) dμₙ`，`Bₙ = ∫ f(x+y₀) dμₙ`，`B = ∫ f(x+y₀) dμ`。

1. 取 `ε > 0`，令 `η = ε/4`，`δ = ε/(8(‖f‖+1))`（`+1` 只为免去 `‖f‖=0` 分情形）。
2. **tightness**：`hTight` ⇒ `∃K` 紧，`∀n, μₙ(Kᶜ) ≤ δ`（`Tight.lean:60`；`ENNReal → ℝ` 转换一次）。
3. **块 I（紧集内）**：`S := {yₙ} ∪ {y₀}` 紧（`hy.isCompact_insert_range`）；`Φ(x,y) := f(x+y)` 在 `K ×ˢ S` 上连续，故 Heine–Cantor ⇒ `∃W₁,W₂ ∈ 𝓤E`：
   `y,y₀ ∈ S ∧ (y,y₀) ∈ W₂ ⇒ ∀x ∈ K, |Φ(x,y) − Φ(x,y₀)| < η`
   （`(x,x) ∈ W₁` 由 `refl_mem_uniformity` 恒成立，**所以只要一个 W₂**）。
4. `yₙ → y₀` ⇒ `∃N, ∀n ≥ N, (yₙ,y₀) ∈ W₂`（`nhds_eq_comap_uniformity'`）。
5. **块 II（紧集外）**：拆绝对值
   `|Aₙ − Bₙ| = |∫ (Φ(·,yₙ) − Φ(·,y₀)) dμₙ| ≤ ∫_K |·| + ∫_{Kᶜ} |·|`。
6. 块 I 代进第 1 项、`‖f‖` 界代进第 2 项：
   `∫_K |·| ≤ η·μₙ(K) ≤ η·1 = η`，`∫_{Kᶜ} |·| ≤ 2‖f‖·μₙ(Kᶜ) ≤ 2‖f‖δ`。
   故 `n ≥ N ⇒ |Aₙ − Bₙ| ≤ η + 2‖f‖δ ≤ ε/4 + ε/2 = 3ε/4`。
7. **块 III（固定移位，白拿）**：`|Bₙ − B| < ε/4`（n 充分大，②b-0）。
8. 三角不等式：`|Aₙ − B| ≤ 3ε/4 + ε/4 = ε`。∎

**注意两处「不是显然」**：(a) 第 3 步的 `W₂` 必须**同时管住所有 n 的 yₙ**——这正是把 `{yₙ} ∪ {y₀}` 整个塞进紧集 `S` 的理由，`yₙ→y₀` 只在第 4 步用**一次**；(b) 第 6 步的 `μₙ(K)` 只能放大到 `1`（概率测度），不能换成 `μₙ(K) ≤ 1` 以外的任何东西。

---

## §2 Lean

文件：`lean/os1/OS1JointShift.lean`（②a 两行原样保留；②b 由 1 条 `sorry` 改为 3 条逐层引理）。

| 声明 | 内容 | 等级 |
|---|---|---|
| `continuous_add_pair` / `addRightHomeo` | ②a（`continuous_add`、`Homeomorph.addRight`） | `[T]` 原样 |
| `tendsto_integral_shift_const_of_tendsto` | ②b-0 固定移位 | `[T]` 原样（`ProbabilityMeasure.lean:346`） |
| `exists_entourage_shift_close` | **块 I**：`K×S` 一致连续 ⇒ 一个 `W₂` 使 `K` 上一致 `< η` | `[T]`，1 处 `[GLUE-1]` |
| `norm_integral_shift_sub_le` | **块 II**：`‖∫f(·+y)−∫f(·+y₀)‖ ≤ η + 2‖f‖·μ.real Kᶜ` | `[T]`，1 处 `[GLUE-2]` |
| `tendsto_integral_shift_of_tendsto` | ②b 主定理（块 I+II+III 合成） | `[T]`，1 处 `[GLUE-3]`（纯算术） |

3 处 `[GLUE]` 的确切位置与替代品：

- `[GLUE-1]`：`hU ∈ 𝓤(E×E)` → `∃W₁,W₂ ∈ 𝓤E`，`(a,a')∈W₁ ∧ (b,b')∈W₂ ⇒ ((a,b),(a',b'))∈hU`。
  用 `uniformity_prod_eq_comap_prod`（`UniformSpace/Basic.lean:757`）+ `Filter.mem_prod_iff`，≈4 行。
- `[GLUE-2]`：`μ.real K ≤ 1` = `measureReal_mono (Set.subset_univ K)`（`Measure/Real.lean:88`，已核）+ `probReal_univ`（`Probability.lean:118`）。
- `[GLUE-3]`：`2‖f‖·ε/(8(‖f‖+1)) ≤ ε/2`，`nlinarith` 一行。

**没有一处 `sorry` 是为了「Mathlib 缺引理」。** 这与 OOO1 §4 ⑤ 的「0 命中」不矛盾：0 命中说的是**没有这条现成引理**，不是**缺件**。

---

## §3 判定

### 3.1 0 `sorry` 可达 —— `[C]`（未编译）

数学无洞（§1），件全在（§3.2）。剩余风险只有命名/`simp` 形态，属编译轮工作量，不是数学工作量。

### 3.2 最小依赖清单（12 件，全部核过 `文件:行`）

| 用途 | Mathlib 名 | 位置 |
|---|---|---|
| 加法联合连续 / 移位同胚 | `continuous_add`、`Homeomorph.addRight` | `Topology/Algebra/Group/Basic.lean:1081`（方向） |
| 弱收敛 ⇔ 所有 `C_b` 积分收敛 | `ProbabilityMeasure.tendsto_iff_forall_integral_tendsto` | `ProbabilityMeasure.lean:346` |
| tight ⇒ `∀ε∃K` | `isTightMeasureSet_iff_exists_isCompact_measure_compl_le` | `Measure/Tight.lean:60` |
| 序列极限集紧 | `Tendsto.isCompact_insert_range` | `Compactness/Compact.lean:634` |
| 积集紧 | `IsCompact.prod` | `HeineCantor.lean:100` 处用法 |
| **块 I 的核心** | `IsCompact.uniformContinuousAt_of_continuousAt` | `Topology/UniformSpace/HeineCantor.lean:58` |
| 目标 ε-带是 entourage | `dist_mem_uniformity` | `MetricSpace/Pseudo/Defs.lean:755` |
| 对角线在 entourage 内 | `refl_mem_uniformity` | `UniformSpace/Basic.lean` |
| 积一致结构 | `uniformity_prod_eq_comap_prod` | `UniformSpace/Basic.lean:757` |
| `yₙ→y₀` 的 entourage 形式 | `nhds_eq_comap_uniformity'` | `UniformSpace/Defs.lean:528` |
| 积分线性 / 范数不等式 / 拆集 | `integral_sub`、`norm_integral_le_integral_norm`、`integral_add_compl` | `Bochner/Basic.lean:261,934`、`Bochner/Set.lean:158` |
| 集积分界 / 常函数积分 | `norm_setIntegral_le_of_norm_le_const`、`setIntegral_const` | `Bochner/Set.lean:616,527` |
| 有界连续可积 | `Integrable.of_bound` | `Integral/IntegrableOn.lean:171` |
| 概率测度 `real univ = 1` | `probReal_univ` | `Measure/Typeclasses/Probability.lean:118` |

**不在清单里（明确不需要）**：`PseudoMetricSpace E`、`SecondCountableTopology E`、`CompleteSpace E`、平移不变度量、任何新 `def`。②b 的载体要求只有 `[UniformSpace E] [IsUniformAddGroup E] [T2Space E] [BorelSpace E]`。

### 3.3 与 Mathlib 最近的现存物（说明「为什么 0 命中」不是工具失败）

`TendstoInDistribution.add_of_tendstoInMeasure_const`（`MeasureTheory/Function/ConvergenceInDistribution.lean:333`）是 Slutsky 型「和」定理，**数学上就是 ②b**，但形式是「同一环境测度 `μ` 上的随机变量」，而 ②b 是「测度族 `μₙ` 变化 + 空间自身平移」——**形式不匹配，不能直接用**。这条记账很重要：它证明 ②b 的内容不新，只是 Mathlib 未在该形式下提供。

---

## §4 检索核验（AGENTS §1 四要素）

| # | 工具状态 | 实际检索式 | 命中 | 备注 |
|---|---|---|---|---|
| ① | `grep` **OK** | `theorem IsCompact\.uniformContinuousOn_of_continuous` @ `HeineCantor.lean` | 3（`:47` 主体、`:100`、`:110`） | — |
| ② | `grep` **OK** | `theorem IsCompact\.uniformContinuousAt_of_continuousAt` @ `HeineCantor.lean` 同文件 | 1（`:58`） | 本轮新增的关键件 |
| ③ | `grep` **OK** | `isTightMeasureSet_iff_exists_isCompact_measure_compl_le\|def IsTightMeasureSet` @ `Measure/Tight.lean` | 6（`:55` `:60` 主体） | — |
| ④ | `grep` **OK** | `theorem integral_add_compl\|norm_integral_le_integral_norm\|setIntegral_const\|norm_setIntegral_le_of_norm_le_const` @ `MeasureTheory/Integral` | 15（`Bochner/Set.lean:158,527,616`；`Bochner/Basic.lean:934`） | — |
| ⑤ | `grep` **OK** | `theorem dist_mem_uniformity\|nhds_eq_comap_uniformity\|uniformity_prod` @ `Topology` | 6（`Pseudo/Defs.lean:755`、`Defs.lean:279,528`、`Basic.lean:746,757,762`） | — |
| ⑥ | `grep` **OK** | `theorem uniformContinuousOn_iff` @ `MetricSpace/Pseudo/Defs.lean` | 2（`:762` `:771`） | 二者**都含** `[PseudoMetricSpace α]` ⇒ 不用 |
| ⑦ | `grep` **TOOL-FAIL → 重试 OK** | `norm_coe_le_norm` @ 错误的 `Topology/Algebra/Order/BoundedContinuousFunction`（目录不存在，exit 2） | 重试命中 3（`ContinuousMap/Bounded/Normed.lean:63`） | 记一次工具失败，已换路重试 |
| ⑧ | `grep` **OK** | `Integrable.of_bound` @ `MeasureTheory` | 77 行中含 `Integral/IntegrableOn.lean:171` 精确签名 | — |
| ⑨ | `grep` **OK** | `probReal_univ` @ `MeasureTheory` | 12（`Typeclasses/Probability.lean:118` 主体） | — |
| ⑩ | `grep` **OK** | 定理名清单 @ `Function/ConvergenceInDistribution.lean` | 20（含 `:333` Slutsky） | §3.3 出处 |
| ⑪ | `grep` **OK** | `real_mono\|real_def\|probReal_univ\|real_univ` @ `MeasureTheory` | 12；随后 `measureReal_mono` @ `Measure/Real.lean` 命中 `:88` ⇒ `[GLUE-2]` 名字已定 | 无遗留未定名 |
| ⑫ | `grep` **OK** | `theorem Real.norm_of_nonneg\|theorem measure_lt_top` @ `Mathlib` | 6（`Measure/Typeclasses/Finite.lean:48` 主体）；`Real.norm_of_nonneg` **0 命中** ⇒ 代码改用 `Real.norm_eq_abs` + `abs_of_nonneg` | 一处名字换路 |

**失败标记**：`TOOL-FAIL` **1**（第 ⑦ 行，路径不存在 → 已重试并取得命中）。其余「未命中」均为**已执行检索的 0/低命中**，非工具失败。**本轮零编译**，全部为静态读源码。

---

## §5 下一手（三条，按依赖序）

1. **编译轮**：给 `lean/` 建 `lakefile.lean`（require mathlib `v4.33.0-rc1`，或直接把 `lean/os1/*.lean` 挂进 `vendor/OSforGFF` 的 `lean_lib`）→ 编译 `OS1Tight` + `OS1JointShift` → 收 3 处 `[GLUE]`。
2. **落账**：`OS1Statement.lean` 的 `os1_map_shift_eq_self` 删 `hJoint` 假设，改吃 `tendsto_integral_shift_of_tendsto`；此时 OS1 链上无 `[NG]`，只剩 `𝒮′` 实体。
3. **位置图**：`docs/position-current.html` 把「OS1 Lean 化：唯一 `[NG]` = ②b」改为「②b 写成（未编译）」，并记 §3.2 载体要求降级（无度量需求）。
