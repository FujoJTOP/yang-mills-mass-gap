/-
  OS1 · 块 ②  「联合推移连续」   —— 分两层判

  ②a  逐点版：`(t,x) ↦ x + t` 连续                          → Mathlib **有**（[T]）
  ②b  测度版：`(μ_n ⇒ μ, y_n → y₀) ⇒ ∫ f(·+y_n) dμ_n → ∫ f(·+y₀) dμ`
      → Mathlib 无此引理，但**不缺件**：本轮把 ε-论证写成（[T]，条件于编译）

  【本轮 · 不编译】（RRR3 令）。名字逐条标 `文件:行`（本机
  vendor/OSforGFF/.lake/packages/mathlib）。`sorry` 处一律标 `[GLUE]`，
  **均非「缺 Mathlib 引理」**，而是「本轮不许编译 ⇒ 未接的机械拼装」。

  ε-分解（与 research/v3/species/RRR3-JOINT.md §1 同一份）：
    块 I  紧集 K 上一致连续：`K × {y_n}∪{y₀}` 紧 + Heine–Cantor ⇒ 一个 entourage
    块 II 紧集外 tightness：`‖f‖·2 · μ_n(Kᶜ)`
    块 III 固定移位（②b-0 白拿）：`|∫ f(·+y₀)dμ_n − ∫ f(·+y₀)dμ|`
-/
import Mathlib.MeasureTheory.Measure.ProbabilityMeasure
import Mathlib.MeasureTheory.Measure.Tight
import Mathlib.MeasureTheory.Integral.Bochner.Set
import Mathlib.MeasureTheory.Integral.IntegrableOn
import Mathlib.Topology.UniformSpace.HeineCantor
import Mathlib.Topology.MetricSpace.Pseudo.Defs
import Mathlib.Topology.Algebra.Group.Basic

open Filter MeasureTheory
open scoped Topology ENNReal

namespace OS1

/-! ## ②a  逐点推移连续 —— [T]，Mathlib 现成 -/

section Pointwise

variable {E : Type*} [AddCommGroup E] [UniformSpace E] [IsUniformAddGroup E]

/-- ②a-1  [T]：加法联合连续（`continuous_mul` 的 to_additive）。
Mathlib：`continuous_add`。缺的怎么补：不缺。 -/
theorem continuous_add_pair : Continuous fun p : E × E => p.1 + p.2 :=
  continuous_add

/-- ②a-2  [T]：固定移位 `x ↦ x + t` 是同胚。
Mathlib：`Homeomorph.addRight t`（方向见 `Topology/Algebra/Group/Basic.lean:1081` doc）。 -/
def addRightHomeo (t : E) : E ≃ₜ E :=
  Homeomorph.addRight t

end Pointwise

/-! ## ②b  测度层：移动试验函数的积分收敛 —— [T]（ε-论证已写出） -/

section MeasureLayer

variable {E : Type*} [AddCommGroup E] [UniformSpace E] [IsUniformAddGroup E] [T2Space E]
  [MeasurableSpace E] [BorelSpace E] [PseudoMetrizableSpace E]

/-- ②b-0  [T]：**固定**试验函数，Mathlib 白拿。
`ProbabilityMeasure.tendsto_iff_forall_integral_tendsto`（ProbabilityMeasure.lean:346）
+ `f.comp ⟨x ↦ x + y₀, continuous_add_const y₀⟩`。
它把 ②b 的缺口压到只剩「`y₀` 换成 `y_n`」。 -/
theorem tendsto_integral_shift_const_of_tendsto
    {μs : ℕ → ProbabilityMeasure E} {μ : ProbabilityMeasure E}
    (h : Tendsto μs atTop (𝓝 μ)) (f : E →ᵇ ℝ) (y₀ : E) :
    Tendsto (fun n => ∫ x, f (x + y₀) ∂(μs n : Measure E)) atTop
      (𝓝 (∫ x, f (x + y₀) ∂(μ : Measure E))) :=
  (ProbabilityMeasure.tendsto_iff_forall_integral_tendsto.1 h)
    (f.comp ⟨fun x => x + y₀, continuous_add_const y₀⟩)

/-- ②b-1  **块 I（紧集内：一致连续）**。
在 `K × S`（`S = {y_n} ∪ {y₀}`）上对 `Φ (x,y) = f (x+y)` 用 Heine–Cantor，
把「`y` 与 `y₀` 在某个 entourage 内」翻译成「`K` 上**逐点**一致地
`|f(x+y) − f(x+y₀)| < η`」，且**不要求 E 有度量**。

用料（逐条 `文件:行`，本轮核过）：
* `IsCompact.prod`：`K ×ˢ S` 紧（`HeineCantor.lean:100` 用同一式）；
* `IsCompact.uniformContinuousAt_of_continuousAt`：`HeineCantor.lean:58`
  —— 返回 `{x : α×α | x.1 ∈ s → (f x.1, f x.2) ∈ r} ∈ 𝓤 α`，正好是「点对 closeness」形式；
* `dist_mem_uniformity`：`Metric/Pseudo/Defs.lean:755` 造 `r ∈ 𝓤 ℝ`；
* `refl_mem_uniformity`：`(x,x) ∈ W₁` 恒成立 ⇒ 不需要 E 上的距离。

`[GLUE-1]`：从 `𝓤 (E×E)` 拆出 `W₁ W₂`，用
`uniformity_prod_eq_comap_prod`（`UniformSpace/Basic.lean:757`）+ `Filter.mem_prod_iff`。 -/
theorem exists_entourage_shift_close
    (f : E →ᵇ ℝ) {K S : Set E} (hK : IsCompact K) (hS : IsCompact S) {η : ℝ} (hη : 0 < η) :
    ∃ W₁ ∈ 𝓤 E, ∃ W₂ ∈ 𝓤 E, ∀ y ∈ S, ∀ y₀ ∈ S, (y, y₀) ∈ W₂ →
      ∀ x ∈ K, dist (f (x + y)) (f (x + y₀)) < η := by
  set Φ : E × E → ℝ := fun p => (f : E → ℝ) (p.1 + p.2) with hΦdef
  have hΦ : Continuous Φ :=
    f.continuous.comp (continuous_fst.add continuous_snd)
  -- 目标带 `r = {p : ℝ × ℝ | dist p.1 p.2 < η} ∈ 𝓤 ℝ`
  have hr : {p : ℝ × ℝ | dist p.1 p.2 < η} ∈ 𝓤 ℝ := dist_mem_uniformity hη
  -- Heine–Cantor at `K ×ˢ S`
  have hU : {p : (E × E) × (E × E) | p.1 ∈ K ×ˢ S → (Φ p.1, Φ p.2) ∈
      {q : ℝ × ℝ | dist q.1 q.2 < η}} ∈ 𝓤 (E × E) :=
    (hK.prod hS).uniformContinuousAt_of_continuousAt Φ
      (fun p _ => hΦ.continuousAt) hr
  -- [GLUE-1] 拆积一致结构
  obtain ⟨W₁, hW₁, W₂, hW₂, hsub⟩ : ∃ W₁ ∈ 𝓤 E, ∃ W₂ ∈ 𝓤 E, ∀ a a' b b' : E,
      (a, a') ∈ W₁ → (b, b') ∈ W₂ →
        ((a, b), (a', b')) ∈ {p : (E × E) × (E × E) | p.1 ∈ K ×ˢ S → (Φ p.1, Φ p.2) ∈
          {q : ℝ × ℝ | dist q.1 q.2 < η}} := by
    sorry -- [GLUE-1] `uniformity_prod_eq_comap_prod` + `mem_prod_iff`，4 行；非缺件
  refine ⟨W₁, hW₁, W₂, hW₂, ?_⟩
  intro y hy y₀ hy₀ hyy x hx
  exact hsub x x y y₀ (refl_mem_uniformity hW₁) hyy ⟨⟨hx, hy⟩, rfl⟩

/-- ②b-2  **块 II（紧集外的 tightness）+ 积分记账**：
`‖∫ f(·+y)dμ − ∫ f(·+y₀)dμ‖ ≤ η + 2‖f‖·μ.real Kᶜ`
（`η` 来自块 I，`μ.real Kᶜ` 来自 tightness）。

用料：
* `integral_sub`：`Integral/Bochner/Basic.lean:261`；
* `Integrable.of_bound`：`Integral/IntegrableOn.lean:171`（`[IsFiniteMeasure μ]`）；
* `norm_integral_le_integral_norm`：`Integral/Bochner/Basic.lean:934`；
* `integral_add_compl (hs) (hfi)`：`Integral/Bochner/Set.lean:158`；
* `norm_setIntegral_le_of_norm_le_const (hs : μ s < ∞) (hC)`：`Bochner/Set.lean:616`
  —— 两块都用它（取 `f := ‖g‖`，`‖∫ ‖g‖‖ = ∫ ‖g‖` 因非负）；
* `probReal_univ`：`Measure/Typeclasses/Probability.lean:118`（`μ.real univ = 1`）。 -/
theorem norm_integral_shift_sub_le
    (f : E →ᵇ ℝ) {μ : Measure E} [IsProbabilityMeasure μ] {K : Set E} (hK : IsCompact K)
    {y y₀ : E} {η : ℝ} (hη : 0 < η) (huc : ∀ x ∈ K, dist (f (x + y)) (f (x + y₀)) < η) :
    ‖∫ x, f (x + y) ∂μ - ∫ x, f (x + y₀) ∂μ‖ ≤ η + 2 * ‖f‖ * μ.real Kᶜ := by
  have hKmeas : MeasurableSet K := hK.isClosed.measurableSet
  have hint₁ : Integrable (fun x => f (x + y)) μ :=
    Integrable.of_bound (f.continuous.comp (continuous_add_const y)).aestronglyMeasurable ‖f‖
      (Eventually.of_forall fun x => f.norm_coe_le_norm x)
  have hint₀ : Integrable (fun x => f (x + y₀)) μ :=
    Integrable.of_bound (f.continuous.comp (continuous_add_const y₀)).aestronglyMeasurable ‖f‖
      (Eventually.of_forall fun x => f.norm_coe_le_norm x)
  have hgnorm : Integrable (fun x => ‖f (x + y) - f (x + y₀)‖) μ := (hint₁.sub hint₀).norm
  rw [integral_sub hint₁ hint₀]
  refine (norm_integral_le_integral_norm _).trans ?_
  rw [← integral_add_compl hKmeas hgnorm]
  -- 块 II-左：K 上 ≤ η（块 I 的逐点界，乘 μ.real K ≤ 1）
  have hleft : ∫ x in K, ‖f (x + y) - f (x + y₀)‖ ∂μ ≤ η * μ.real K := by
    have h := norm_setIntegral_le_of_norm_le_const (μ := μ)
      (f := fun x => ‖f (x + y) - f (x + y₀)‖) (C := η) (measure_lt_top μ K)
      (fun x hx => by
        rw [Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
        exact le_of_lt (huc x hx))
    rwa [Real.norm_eq_abs, abs_of_nonneg (setIntegral_nonneg hKmeas fun x _ => norm_nonneg _)] at h
  -- 块 II-右：Kᶜ 上 ≤ 2‖f‖·μ.real Kᶜ（tightness 那半）
  have hright : ∫ x in Kᶜ, ‖f (x + y) - f (x + y₀)‖ ∂μ ≤ 2 * ‖f‖ * μ.real Kᶜ := by
    have h := norm_setIntegral_le_of_norm_le_const (μ := μ)
      (f := fun x => ‖f (x + y) - f (x + y₀)‖) (C := 2 * ‖f‖) (measure_lt_top μ Kᶜ)
      (fun x _ => by
        rw [Real.norm_eq_abs, abs_of_nonneg (norm_nonneg _)]
        calc ‖f (x + y) - f (x + y₀)‖
            ≤ ‖f (x + y)‖ + ‖f (x + y₀)‖ := norm_sub_le _ _
          _ ≤ ‖f‖ + ‖f‖ := add_le_add (f.norm_coe_le_norm _) (f.norm_coe_le_norm _)
          _ = 2 * ‖f‖ := by ring)
    rwa [Real.norm_eq_abs, abs_of_nonneg (setIntegral_nonneg hKmeas.compl fun x _ => norm_nonneg _)] at h
  -- 拼：`η * μ.real K ≤ η * 1 = η`
  have hKle : μ.real K ≤ 1 := by
    simpa [probReal_univ] using measureReal_mono (μ := μ) (Set.subset_univ K) -- 已核 `Measure/Real.lean:88`
  calc ∫ x in K, ‖f (x + y) - f (x + y₀)‖ ∂μ + ∫ x in Kᶜ, ‖f (x + y) - f (x + y₀)‖ ∂μ
      ≤ η * μ.real K + 2 * ‖f‖ * μ.real Kᶜ := add_le_add hleft hright
    _ ≤ η * 1 + 2 * ‖f‖ * μ.real Kᶜ :=
        add_le_add (mul_le_mul_of_nonneg_left hKle hη.le) le_rfl
    _ = η + 2 * ‖f‖ * μ.real Kᶜ := by ring

/-- ②b  **主定理（本轮写成）**：联合推移连续的测度版。

三块合成（见文件头与 RRR3-JOINT.md §1）：
`|∫f(·+y_n)dμ_n − ∫f(·+y₀)dμ| ≤ [块I+块II] + [块III] ≤ (ε/4 + ε/2) + ε/4 < ε`。

`[GLUE-3]`：纯算术 `2‖f‖·ε/(8(‖f‖+1)) ≤ ε/2`；`[NAME]`：`Metric.tendsto_atTop`、
`Measure.real`、`ENNReal.toReal_mono`、`nhds_eq_comap_uniformity'`（已核 `Defs.lean:528`）。 -/
theorem tendsto_integral_shift_of_tendsto
    {μs : ℕ → ProbabilityMeasure E} {μ : ProbabilityMeasure E}
    {y : ℕ → E} {y₀ : E} (f : E →ᵇ ℝ)
    (hTight : IsTightMeasureSet (Set.range fun n => (μs n : Measure E)))
    (hμ : Tendsto μs atTop (𝓝 μ)) (hy : Tendsto y atTop (𝓝 y₀)) :
    Tendsto (fun n => ∫ x, f (x + y n) ∂(μs n : Measure E)) atTop
      (𝓝 (∫ x, f (x + y₀) ∂(μ : Measure E))) := by
  have hS : IsCompact (insert y₀ (Set.range y)) := hy.isCompact_insert_range
  have hy₀S : y₀ ∈ insert y₀ (Set.range y) := Set.mem_insert _ _
  have hyS : ∀ n, y n ∈ insert y₀ (Set.range y) := fun n => Set.mem_insert_of_mem _ ⟨n, rfl⟩
  -- ① 差值 `Δ_n = ∫f(·+y_n)dμ_n − ∫f(·+y₀)dμ_n` 一致地 → 0
  have hΔ : Tendsto (fun n => ∫ x, f (x + y n) ∂(μs n : Measure E)
      - ∫ x, f (x + y₀) ∂(μs n : Measure E)) atTop (𝓝 0) := by
    rw [Metric.tendsto_atTop]
    intro ε hε
    have hη4 : (0 : ℝ) < ε / 4 := by positivity
    have hδ : (0 : ℝ) < ε / (8 * (‖f‖ + 1)) := by positivity
    -- tightness：一个 K 管住整族
    obtain ⟨K, hK, hKμ⟩ :=
      isTightMeasureSet_iff_exists_isCompact_measure_compl_le.1 hTight
        (ENNReal.ofReal (ε / (8 * (‖f‖ + 1)))) (ENNReal.ofReal_pos.mpr hδ)
    -- 块 I
    obtain ⟨W₁, -, W₂, hW₂, hclose⟩ := exists_entourage_shift_close f hK hS hη4
    have hyW : ∀ᶠ n in atTop, (y n, y₀) ∈ W₂ := by
      have hnb : {z : E | (z, y₀) ∈ W₂} ∈ 𝓝 y₀ := by
        rw [nhds_eq_comap_uniformity']
        exact hW₂
      exact hy.eventually hnb
    obtain ⟨N, hN⟩ := eventually_atTop.1 hyW
    refine eventually_atTop.2 ⟨N, fun n hn => ?_⟩
    -- 块 II 的数值界
    have key := norm_integral_shift_sub_le f hK (y := y n) (y₀ := y₀) hη4
      (fun x hx => hclose (y n) (hyS n) y₀ hy₀S (hN n hn) x hx)
    have hcompl : (μs n : Measure E).real Kᶜ ≤ ε / (8 * (‖f‖ + 1)) := by
      have h1 : (μs n : Measure E) Kᶜ ≤ ENNReal.ofReal (ε / (8 * (‖f‖ + 1))) :=
        hKμ _ (Set.mem_range_self n)
      have h2 : ((μs n : Measure E) Kᶜ).toReal
          ≤ (ENNReal.ofReal (ε / (8 * (‖f‖ + 1)))).toReal :=
        ENNReal.toReal_mono ENNReal.ofReal_ne_top h1
      rwa [ENNReal.toReal_ofReal hδ.le, Measure.real] at h2
    have harith : 2 * ‖f‖ * (ε / (8 * (‖f‖ + 1))) ≤ ε / 2 := by
      have h0 : (0 : ℝ) ≤ ‖f‖ := norm_nonneg f
      rw [mul_div_assoc', div_le_iff₀ (by positivity)] -- [NAME] `div_le_iff₀`
      nlinarith [h0, hε.le]                            -- [GLUE-3] 纯算术
    rw [dist_zero_right]
    calc ‖∫ x, f (x + y n) ∂(μs n : Measure E) - ∫ x, f (x + y₀) ∂(μs n : Measure E)‖
        ≤ ε / 4 + 2 * ‖f‖ * (μs n : Measure E).real Kᶜ := key
      _ ≤ ε / 4 + 2 * ‖f‖ * (ε / (8 * (‖f‖ + 1))) := by
          gcongr
      _ ≤ ε / 4 + ε / 2 := by linarith
      _ < ε := by linarith
  -- ② 块 III + 拼装
  have hconst := tendsto_integral_shift_const_of_tendsto hμ f y₀
  have hsum := hΔ.add hconst
  simpa [sub_add_cancel] using hsum

end MeasureLayer

end OS1
