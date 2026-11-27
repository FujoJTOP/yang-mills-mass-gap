/-
  OS1 · 块 ①  「一般 Polish：弱收敛 ⇒ tight」   —— Mathlib 组装（[T]）

  口径：本文件**未编译**（OOO1 令：只写代码 + 逐行标注依赖）。
        每条 `-- [T]` 标出所用 Mathlib 名字与 `文件:行`
        （本机 vendor/OSforGFF/.lake/packages/mathlib，toolchain v4.33.0-rc1）。
        本文件不含 `sorry`，不含新定义。

  NNN1 §2 的判「缺」需精确一格：
    Mathlib 缺的是 Prokhorov 方向（tight ⇒ 相对紧）；
    OS1 需要的方向（收敛 ⇒ tight）可走「收敛 ⇒ 相对紧」+「相对紧 ⇒ tight」两步，
    两步**都在 Mathlib**（Compact.lean:634/103、Prokhorov.lean:609）。
    故 ① 是组装，不是新数学。
-/
import Mathlib.MeasureTheory.Measure.Prokhorov
import Mathlib.MeasureTheory.Measure.ProbabilityMeasure
import Mathlib.Topology.Compactness.Compact

open Filter Set MeasureTheory
open scoped Topology ENNReal

namespace OS1

/-- ①  弱收敛序列的像集 tight（在 `Measure E` 上）。

需要什么：
* `E` Polish：`PseudoMetricSpace` + `SecondCountableTopology` + `CompleteSpace`
  （`isTightMeasureSet_of_isCompact_closure` 的变量块，Prokhorov.lean:545 与 :605）。
* `[BorelSpace E]`：`ProbabilityMeasure` 的弱拓扑实例需要它；
  并且 `PseudoMetricSpace E` 顺带给出 `HasOuterApproxClosed E`（HasOuterApproxClosed.lean:218），
  于是 `T2Space (ProbabilityMeasure E)`（ProbabilityMeasure.lean:422）到手，`IsCompact.isClosed` 可用。

Mathlib 里有没有：**有**（本文件只组装，见步 1–6）。
缺的怎么补：**不缺**。
-/
theorem isTightMeasureSet_range_of_tendsto
    {E : Type*} [PseudoMetricSpace E] [MeasurableSpace E] [BorelSpace E]
    [SecondCountableTopology E] [CompleteSpace E]
    {μs : ℕ → ProbabilityMeasure E} {μ : ProbabilityMeasure E}
    (h : Tendsto μs atTop (𝓝 μ)) :
    IsTightMeasureSet (Set.range fun n => (μs n : Measure E)) := by
  -- [T] 步 1：ℕ-序列收敛 ⇒ `insert 极限 (range)` 紧
  --     Mathlib: `Tendsto.isCompact_insert_range`  Topology/Compactness/Compact.lean:634
  have hKi : IsCompact (insert μ (Set.range μs)) := h.isCompact_insert_range
  -- [T] 步 2：紧集在 T2 里闭
  --     `IsCompact.isClosed`（Mathlib 现成，需 `[T2Space (ProbabilityMeasure E)]`：
  --      ProbabilityMeasure.lean:422，条件是 `[HasOuterApproxClosed E] [BorelSpace E]`）
  have hI : IsClosed (insert μ (Set.range μs)) := hKi.isClosed
  -- [T] 步 3：闭集 ⊇ range ⇒ ⊇ 闭包（`closure_minimal`，Mathlib 现成）
  have hsub : closure (Set.range μs) ⊆ insert μ (Set.range μs) :=
    closure_minimal (subset_insert μ (Set.range μs)) hI
  -- [T] 步 4：紧集的闭子集紧。Mathlib: `IsCompact.of_isClosed_subset` Compact.lean:103
  have hK : IsCompact (closure (Set.range μs)) :=
    hKi.of_isClosed_subset isClosed_closure hsub
  -- [T] 步 5：相对紧 ⇒ tight。Mathlib: `isTightMeasureSet_of_isCompact_closure`
  --     MeasureTheory/Measure/Prokhorov.lean:609（它是 `S : Set (ProbabilityMeasure E)` 版本，
  --     结论形如 `IsTightMeasureSet {↑μ' | μ' ∈ S}`）
  have hT := isTightMeasureSet_of_isCompact_closure (S := Set.range μs) hK
  -- [T] 步 6：集合写法对齐 `{↑μ' | μ' ∈ range μs}` ↔ `range (↑ ∘ μs)`
  --     （可能需要 `simp only [Set.range, Set.image_eq_range]`；名字待核，不影响组装性）
  simpa [Set.range] using hT

/-- ①′  同一件的「不动点」形式：把 ① 喂给 `IsTightMeasureSet` 的展开判据。
本引理是纯改写（`Tight.lean:60` 的 `isTightMeasureSet_iff_exists_isCompact_measure_compl_le`），
不需要新内容；用它是为了在 OS1 主定理里以 `∀ ε, ∃ K` 的形状消费 tightness。 -/
theorem exists_compact_measure_compl_le_of_tendsto
    {E : Type*} [PseudoMetricSpace E] [MeasurableSpace E] [BorelSpace E]
    [SecondCountableTopology E] [CompleteSpace E]
    {μs : ℕ → ProbabilityMeasure E} {μ : ProbabilityMeasure E}
    (h : Tendsto μs atTop (𝓝 μ)) :
    ∀ ε, 0 < ε → ∃ K : Set E, IsCompact K ∧ ∀ n, (μs n : Measure E) (Kᶜ) ≤ ε := by
  -- [T] 用 ① + `isTightMeasureSet_iff_exists_isCompact_measure_compl_le`（Tight.lean:60）
  --     再对 `Set.range` 做一次 `mem_range` 展开（`∀ μ' ∈ range μs` ↔ `∀ n`）。
  have hT := isTightMeasureSet_range_of_tendsto h
  rw [isTightMeasureSet_iff_exists_isCompact_measure_compl_le] at hT
  intro ε hε
  obtain ⟨K, hK, hμK⟩ := hT ε hε
  exact ⟨K, hK, fun n => hμK _ (Set.mem_range_self n)⟩

end OS1
