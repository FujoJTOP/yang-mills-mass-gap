/-
  OS1 · 主定理骨架（[T] 组装，条件于 ②b 的 `hJoint`）

  这是「最小可编译子集」的核心：不含 `sorry`，只用 Mathlib 现成名字 + 块 ① 的结论。
  三条前提逐条点名（照 OS1 台账）：
    hConv  —— 极限已存在（= GL3 / 消费 UVS；OS1 携带不了）
    hInv   —— 格点 `G_a ∋ y_n` 下的精确不变性（D3；可放宽到 ε_n→0 仍成立）
    hJoint —— ②b（Mathlib 缺；用 `hMinlos` 体例当显式假设）
  块 ①（tight）**不再需要当假设**：由 OS1Tight.lean 的
  `isTightMeasureSet_range_of_tendsto` 从 `hConv` + Polish 直接给。

  未编译（OOO1 令）。名字待核处已标 `-- 名字待核`。
-/
import Mathlib.MeasureTheory.Measure.ProbabilityMeasure
import Mathlib.MeasureTheory.Measure.FiniteMeasure
import Mathlib.MeasureTheory.Measure.HasOuterApproxClosed
import Mathlib.MeasureTheory.Integral.Lebesgue.Basic
import Mathlib.Topology.Algebra.Group.Basic

open Filter MeasureTheory
open scoped Topology ENNReal

namespace OS1

/-! ### 分离步（OS1 §1 S5）—— [T] -/

/-- [T] `C_b` 分离概率测度：两条概率测度若对所有有界连续函数积分相同则相等。

需要什么：`[TopologicalSpace E] [HasOuterApproxClosed E] [BorelSpace E]`
  （`PseudoMetrizableSpace E` 自动给 `HasOuterApproxClosed`，HasOuterApproxClosed.lean:218）。
Mathlib 里有没有：**有** —— `FiniteMeasure.ext_of_forall_integral_eq`
  （MeasureTheory/Measure/FiniteMeasure.lean:379）。
缺的怎么补：只差一层包装（`ProbabilityMeasure` 是 `FiniteMeasure` 的子类型，
  TemperedDistribution 那种 `def` 包装用 `Subtype.ext`/`ext` 拆一层；本项目自写 ≤5 行）。 -/
theorem probabilityMeasure_ext_of_forall_integral_eq
    {E : Type*} [TopologicalSpace E] [MeasurableSpace E] [HasOuterApproxClosed E]
    [BorelSpace E] {μ ν : ProbabilityMeasure E}
    (h : ∀ f : E →ᵇ ℝ, ∫ x, f x ∂(μ : Measure E) = ∫ x, f x ∂(ν : Measure E)) :
    μ = ν := by
  -- 名字待核：`Subtype.ext` 或 `ext`（`ProbabilityMeasure` 是 def-包装的 subtype）
  refine Subtype.ext ?_
  exact FiniteMeasure.ext_of_forall_integral_eq h

/-! ### OS1 主定理 —— [T] 组装（条件于 hJoint = ②b） -/

/-- OS1（Lean 骨架）：若 `μ_n ⇒ μ`、`y_n → y₀`、`μ_n` 在每个 `y_n` 下精确不变，
且「移动试验函数」的积分收敛（②b），则 `μ` 在 `y₀` 下不变。

需要什么：
* `[AddCommGroup E] [UniformSpace E] [IsUniformAddGroup E] [T2Space E]`
  —— 让 `x ↦ x + c` 连续 + 可测（②a 用的就是这些）；
* `[PseudoMetrizableSpace E] [BorelSpace E]` —— 分离步要 `HasOuterApproxClosed`；
* 三条前提 `hConv / hInv / hJoint`（见文件头）。
Mathlib 里有没有：**组装件全有**；`hJoint` 本身缺（见 OS1JointShift.lean:②b）。
缺的怎么补：`hJoint` 按 `hMinlos` 体例留作显式假设；一旦 ②b 写成，删假设即可。 -/
theorem os1_map_shift_eq_self
    {E : Type*} [AddCommGroup E] [UniformSpace E] [IsUniformAddGroup E] [T2Space E]
    [MeasurableSpace E] [BorelSpace E] [PseudoMetrizableSpace E]
    {μs : ℕ → ProbabilityMeasure E} {μ : ProbabilityMeasure E} {y : ℕ → E} {y₀ : E}
    (hConv : Tendsto μs atTop (𝓝 μ))
    (hy : Tendsto y atTop (𝓝 y₀))
    (hInv : ∀ n, (μs n).map ((continuous_add_const (y n)).measurable.aemeasurable) = μs n)
    (hJoint : ∀ f : E →ᵇ ℝ, Tendsto (fun n => ∫ x, f (x + y n) ∂(μs n : Measure E)) atTop
      (𝓝 (∫ x, f (x + y₀) ∂(μ : Measure E)))) :
    (μ.map ((continuous_add_const y₀)).measurable.aemeasurable) = μ := by
  -- [T] 分离步：只需对所有 `f : E →ᵇ ℝ` 证 `∫ f d(μ.map τ) = ∫ f dμ`
  refine probabilityMeasure_ext_of_forall_integral_eq ?_
  intro f
  -- [T] 左端改写：`∫ f d(μ.map τ) = ∫ f(·+y₀) dμ`
  --     `MeasureTheory.Measure.integral_map`（名字待核：`integral_map` / `integral_map_of_stronglyMeasurable`）
  --     + `ProbabilityMeasure.coe_map`（`ProbabilityMeasure.map` 与 `Measure.map` 的系数桥，名字待核）
  rw [ProbabilityMeasure.coe_map, MeasureTheory.integral_map]
  -- [T] 一步到位的唯一性：`∫ f(·+y_n) dμ_n` 有两个极限
  --     (a) 由 `hJoint`：`∫ f(·+y₀) dμ`
  --     (b) 由 `hInv` + 弱收敛：`∫ f dμ`
  refine tendsto_nhds_unique (hJoint f) ?_
  -- [T] (b)：`hInv` 把 `∫ f(·+y_n) dμ_n` 逐项改回 `∫ f dμ_n`，
  --     再用 `ProbabilityMeasure.tendsto_iff_forall_integral_tendsto`（ProbabilityMeasure.lean:346）
  have hcongr : ∀ n, ∫ x, f (x + y n) ∂(μs n : Measure E) = ∫ x, f x ∂(μs n : Measure E) := by
    intro n
    have hstep : ∫ x, f (x + y n) ∂(μs n : Measure E)
        = ∫ x, f x ∂((μs n).map ((continuous_add_const (y n)).measurable.aemeasurable) :
            Measure E) := by
      rw [MeasureTheory.integral_map]  -- 名字待核（同上一处）
    simpa [hInv n] using hstep.symm
  simpa only [hcongr] using
    (ProbabilityMeasure.tendsto_iff_forall_integral_tendsto.1 hConv) f

end OS1
