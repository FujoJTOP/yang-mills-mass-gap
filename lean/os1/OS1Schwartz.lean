/-
  OS1 · 块 ③  𝒮′(ℝ⁴) 实例   —— [NG]

  判据（本节结论，非 Lean 证明）：缺的**不是一条 instance**，是**载体本身**。
  Mathlib 把 `TemperedDistribution E F` 定义成「点态收敛拓扑下的连续线性映射」
  （Analysis/Distribution/TemperedDistribution.lean:57，作者自注：因为强对偶的
  「线性映射连续性/序列收敛可由点态检验」的定理缺，所以先用点态拓扑，代价是
  「非 barrelled」，一致有界原理失效）。

  于是一个**具体阻碍**（不是「还没人写」）：
    ① 要喂块 ① 的 Prokhorov 件，需要 `PseudoMetricSpace S'` + `SecondCountableTopology S'`
       + `CompleteSpace S'`（Prokhorov.lean:545/605）；
    ② 但点态拓扑是「不可数多坐标」的积拓扑（`UniformConvergenceCLM σ F {s | Finite s}`，
       PointwiseConvergenceCLM.lean:60）。不可数积不可度量化；`E →Lₚₜ` 作为子空间要被可数坐标
       决定，需要一致有界型结构，而作者已注明 UBP 在此失效。
    ⇒ `PseudoMetricSpace S'` / `SecondCountableTopology S'` **不是缺件，是与当前定义不相容**。
       要 Polish 必须换载体：强对偶拓扑（Fréchet–Montel 型），或退到「可数稠密试验族 + 等度连续」。

  因此本文件的 [T] 只有「Borel 记账」两行；其余每条均为 [NG] 并注缺什么。

  未编译（OOO1 令）。
-/
import Mathlib.Analysis.Distribution.TemperedDistribution
import Mathlib.Analysis.Distribution.SchwartzSpace.Basic
import Mathlib.Topology.MetricSpace.Polish

open SchwartzMap MeasureTheory
open scoped SchwartzMap Topology

namespace OS1
namespace S4

/-- `E4 = ℝ⁴`（只需 `NormedAddCommGroup` + `NormedSpace ℝ`）。 -/
abbrev E4 : Type := EuclideanSpace ℝ (Fin 4)

/-- `S' = 𝒮′(ℝ⁴) = 𝒮(ℝ⁴,ℂ) →Lₚₜ[ℂ] ℂ`。 -/
abbrev S' : Type := 𝓢'(E4, ℂ)

/-! ### 已有（Mathlib 现成，import 即得）—— 供块 ② 用 -/

-- [T] 加法群：UniformConvergenceCLM.lean:167
example : AddCommGroup S' := inferInstance
-- [T] ℂ-模：UniformConvergenceCLM.lean:238
example : Module ℂ S' := inferInstance
-- [T] 一致加法群（这是「(t,x) ↦ x+t 连续」与「平移是均匀同胚」的来源）：
--     UniformConvergenceCLM.lean:197（uniform）、:203（topological）
example : IsUniformAddGroup S' := inferInstance
example : IsTopologicalAddGroup S' := inferInstance
-- [T] T2：PointwiseConvergenceCLM.lean:70
example : T2Space S' := inferInstance

/-! ### 缺件 · Borel 记账（两行，可真写成） —— [T] -/

/-- [T] `S'` 上的 Borel σ-代数。需要：`S'` 已有拓扑（有）。缺的怎么补：不缺。
注意：这两行**足够**让 `ProbabilityMeasure S'` 有类型，但**不够**让块 ① 可用（见下）。 -/
noncomputable instance instMeasurableSpaceS' : MeasurableSpace S' := borel S'

/-- [T] 上述 σ-代数就是 Borel。 -/
instance instBorelSpaceS' : BorelSpace S' := ⟨rfl⟩

/-! ### 缺件 · Polish（[NG]，且与当前定义不相容） -/

-- [NG] 目标：`instance : PseudoMetricSpace S'` —— 未写。
--   缺什么：一条把点态收敛装进单个度量的构造。点态拓扑是 `{s | Finite s}` 上的
--   `UniformConvergenceCLM`（PointwiseConvergenceCLM.lean:60），坐标集不可数 ⇒ 不存在
--   可数邻域基。要么换强拓扑（新定义 `StrongDual`），要么先证「可数子族决定拓扑」。
--   这条不是 instance 漏写，是会**假**的命题，故不落 `sorry`，只记目标。
--
-- [NG] 目标：`instance : SecondCountableTopology S'` —— 未写。同上（不可数积）。
--
-- [NG] 目标：`instance : CompleteSpace S'` —— 未写。点态拓扑下极限函数未必连续
--   （强对偶的完备性靠 Banach–Steinhaus / barrelled，作者注释已说此拓扑不 barrelled）。

/-! ### 缺件 · 平移作用 `τ_t` 与联合连续（[NG]，需新定义） -/

/-- [NG] 平移作用：`⟨τ_t T, f⟩ = ⟨T, f(·+t)⟩`。
实现：`T.comp (SchwartzMap.compSubConstCLM ℂ (-t))`，因为
`f.compSubConstCLM ℂ a x = f (x - a)`（SchwartzSpace/Basic.lean:1073）。
需要什么：无（定义本身一行）；**新定义** ⇒ 计入 [NG]。
缺的怎么补：先证 `shift 0 = id`、`shift (s+t) = shift s ∘ shift t`（后者可由
`compSubConstCLM_comp`，Basic.lean:1082，改写成加法律），再证
`Continuous fun p : E4 × S' => shift p.1 p.2`。 -/
def shift (t : E4) (T : S') : S' :=
  T.comp (SchwartzMap.compSubConstCLM ℂ (-t))

/-- [NG] 联合连续：`(t, T) ↦ τ_t T` 连续（`S'` 用点态拓扑）。
缺什么（Mathlib 0 命中）：
* 试验函数侧的联合连续 `Continuous fun p : E4 × 𝓢(E4,ℂ) => f.compSubConstCLM ℂ p.1`
  —— Mathlib 只有 `compSubConstCLM : E → (𝓢 →L 𝓢)` 的**单点** CLM
  （Basic.lean:1068），没有「对 `a` 联合连续」；
* 有了它，`PointwiseConvergenceCLM.continuous_of_continuous_eval`
  （PointwiseConvergenceCLM.lean:130）把 `S'` 侧的连续化为逐点检验；
* 另需 `IsUniformAddGroup S'` 走 `IsCompact.uniformContinuousOn_of_continuous` 收尾（块 ②b）。 -/
theorem continuous_shift_pair : Continuous fun p : E4 × S' => shift p.1 p.2 := by
  sorry

/-- [NG] OS1 需要的最后一件：`S'` 上的弱收敛 + 平移不变 ⇒ 平移不变极限（块 ①+②+上）。

这条就是 OS1 主定理在 `E = S'` 的实例化。它 [NG] 的直接原因是上文三条 Polish 缺件，
不是缺一条引理：连 `ProbabilityMeasure S'` 的弱拓扑上的 Prokhorov 件都进不去。 -/
theorem os1_at_schwartz
    {μs : ℕ → ProbabilityMeasure S'} {μ : ProbabilityMeasure S'} {y : ℕ → E4} {y₀ : E4}
    (hConv : Tendsto μs atTop (𝓝 μ)) (hy : Tendsto y atTop (𝓝 y₀))
    (hInv : ∀ n, (μs n).map ((continuous_shift_pair.continuousAt.comp
      (continuous_const.prodMk continuous_id)).measurable.aemeasurable) = μs n) :
    True := by
  -- [NG] 真结论是 `μ.map (shift y₀) = μ`，但它的前提（Polish on S'、联合连续、可测性）
  -- 都还没到手；此处只留类型骨架 `True`，避免把未到手的东西写成定理。
  trivial

end S4
end OS1
