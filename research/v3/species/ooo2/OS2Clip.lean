/-
OS2Clip.lean — ① 连续截断 `clipₙ`：`trunc` 的连续替代品

问题（不改既有文件，只在此重建）
--------------------------------
`research/glue-a-limits/GL4L2.lean:28–29`
    noncomputable def trunc (A : Ω → ℝ) (n : ℕ) (ω : Ω) : ℝ := if |A ω| ≤ (n : ℝ) then A ω else 0
已有 `trunc_abs_le`(:31) 给 `rp_l2_closure` 的 `hbound`、`trunc_tendsto`(:37) 给 `hpt`，
但 `trunc` **不连续**；而 `research/glue-a-limits/GL4Transfer.lean:48` 的 `rp_transfer`
要 `A : Ω →ᵇ ℝ`（有界 **连续**）。交接断在「有界连续层 ⇒ 截断层」（NNN2-OS2.md §3 第 3 行）。

本文件给 `clipₙ x = max (−n) (min n x)`：
  `continuous_clip` : `Continuous (clip n)`            ← `trunc` 缺的正是这条
  `clip_abs_le`     : `|clip n x| ≤ |x|`               ⇒ `rp_l2_closure` 的 `hbound`
  `clip_tendsto`    : `clip n x → x`                   ⇒ `rp_l2_closure` 的 `hpt`
  `clip_abs_le_n`   : `|clip n x| ≤ n`                 ⇒ 有界，造 `Ω →ᵇ ℝ` 用
  `clip_eq_self`    : `|x| ≤ n ⇒ clip n x = x`         ⇒ 截断区外与 `trunc` 同语义
  `clipBCF`         : `A` 连续（**不必有界**）⇒ `clipₙ ∘ A : Ω →ᵇ ℝ`
  `clip_memLp`      : `A ∈ L² ⇒ clipₙ ∘ A ∈ L²`

依赖核验（mathlib4 v4.33.0-rc1，逐个核到源；本文件**未编译**）
--------------------------------------------------------------
`Continuous.max` / `Continuous.min`     Topology/Order/OrderClosed.lean:695,697
`continuous_const` / `continuous_id`    Topology/Continuous.lean:161
`abs_le`（`|a| ≤ b ↔ -b ≤ a ∧ a ≤ b`）  Algebra/Order/Ring/Abs.lean:123 的 `.mpr ⟨h1,h2⟩` 用法
`min_eq_left` / `min_eq_right` / `min_le_left` / `max_eq_right` / `le_max_left` / `max_le`
                                        Order/Defs/LinearOrder.lean:155–200
`exists_nat_gt` / `Ici_mem_atTop` / `tendsto_const_nhds.congr'`   同 GL4L2.lean:39–45
`BoundedContinuousFunction.ofNormedAddCommGroup`
                                        Topology/ContinuousMap/Bounded/Normed.lean:118
`MemLp.of_le`（= `MemLp.mono`）         MeasureTheory/Function/LpSeminorm/Basic.lean:500
`Continuous.aestronglyMeasurable`       MeasureTheory/Function/StronglyMeasurable/AEStronglyMeasurable.lean:239
-/
import Mathlib.Topology.Order.OrderClosed
import Mathlib.Topology.ContinuousMap.Bounded.Normed
import Mathlib.MeasureTheory.Function.L2Space
import Mathlib.MeasureTheory.Measure.Portmanteau
import Mathlib.Algebra.Order.Ring.Abs

open MeasureTheory Filter Topology
open scoped BoundedContinuousFunction Topology

namespace GlueA.OS2

variable {Ω : Type*} [TopologicalSpace Ω] [MeasurableSpace Ω] [OpensMeasurableSpace Ω]

/-! ## 1. 定义与连续性 -/

/-- 连续截断：`clipₙ x = max (−n) (min n x)`，把 `x` 夹进 `[−n, n]`。 -/
noncomputable def clip (n : ℕ) (x : ℝ) : ℝ := max (-(n : ℝ)) (min (n : ℝ) x)

/-- **① 的核心**：`clipₙ` 连续 —— `trunc` 恰好缺这一条。 -/
lemma continuous_clip (n : ℕ) : Continuous (clip n) :=
  continuous_const.max (continuous_const.min continuous_id)

/-! ## 2. 三条给 `rp_l2_closure` 用的估计 -/

/-- `|clipₙ x| ≤ |x|`（⇒ `rp_l2_closure` 的 `hbound`）。 -/
lemma clip_abs_le (n : ℕ) (x : ℝ) : |clip n x| ≤ |x| := by
  rw [clip, abs_le]                                  -- 目标拆成 `-|x| ≤ …` 与 `… ≤ |x|`
  refine ⟨?_, ?_⟩
  · rcases le_or_gt (n : ℝ) |x| with h | h
    · exact le_trans (neg_le_neg h) (le_max_left _ _)          -- `-|x| ≤ -n ≤ max`
    · have hx : -(|x|) ≤ x := (abs_le.mp le_rfl).1             -- `|x| ≤ |x|` 取下界
      have hxn : x ≤ (n : ℝ) :=
        le_of_lt (lt_of_le_of_lt (abs_le.mp le_rfl).2 h)       -- `x ≤ |x| < n`
      rw [min_eq_right hxn]                                    -- `min n x = x`
      exact le_trans hx (le_max_right _ _)
  · apply max_le
    · exact le_trans (neg_nonpos.mpr (Nat.cast_nonneg n)) (abs_nonneg x)
    · rcases le_or_gt (n : ℝ) x with h | h
      · rw [min_eq_left h]                                     -- `min n x = n`
        exact le_trans h (abs_le.mp le_rfl).2
      · rw [min_eq_right (le_of_lt h)]                         -- `min n x = x`
        exact (abs_le.mp le_rfl).2

/-- `‖clipₙ x‖ ≤ ‖x‖`：`ℝ` 上的范数即绝对值，供 `hbound` 直接使用。 -/
lemma clip_norm_le (n : ℕ) (x : ℝ) : ‖clip n x‖ ≤ ‖x‖ := by
  simpa [Real.norm_eq_abs] using clip_abs_le n x

/-- 截断区外 `clipₙ` 就是 `x`：`trunc` 的语义（原值/夹值）被保留。 -/
lemma clip_eq_self {n : ℕ} {x : ℝ} (h : |x| ≤ (n : ℝ)) : clip n x = x := by
  rw [clip, min_eq_right (abs_le.mp h).2, max_eq_right (abs_le.mp h).1]

/-- `clipₙ` 逐点收敛到 `x`（⇒ `rp_l2_closure` 的 `hpt`）。 -/
lemma clip_tendsto (x : ℝ) : Tendsto (fun n : ℕ => clip n x) atTop (𝓝 x) := by
  obtain ⟨N, hN⟩ := exists_nat_gt |x|                            -- `|x| < N`
  have hev : ∀ᶠ n : ℕ in atTop, clip n x = x := by
    filter_upwards [Ici_mem_atTop N] with n hn
    exact clip_eq_self (le_of_lt (lt_of_lt_of_le hN (by exact_mod_cast hn)))
  exact tendsto_const_nhds.congr' (hev.mono fun n h => h.symm)

/-- `|clipₙ x| ≤ n`：有界性（`trunc` 也有，但没有连续性配它）。 -/
lemma clip_abs_le_n (n : ℕ) (x : ℝ) : |clip n x| ≤ (n : ℝ) := by
  rw [clip, abs_le]
  refine ⟨le_max_left _ _, ?_⟩
  apply max_le
  · exact le_trans (neg_nonpos.mpr (Nat.cast_nonneg n)) (Nat.cast_nonneg n)
  · exact min_le_left _ _

/-! ## 3. `clipₙ ∘ A` 作为有界连续函数（`rp_transfer` 的槽位） -/

/-- **① 的接口**：`A` 只要连续（**不必有界**），`clipₙ ∘ A` 就是 `Ω →ᵇ ℝ`。 -/
noncomputable def clipBCF (n : ℕ) (A : Ω → ℝ) (hA : Continuous A) : Ω →ᵇ ℝ :=
  BoundedContinuousFunction.ofNormedAddCommGroup (fun ω => clip n (A ω))
    ((continuous_clip n).comp hA) (n : ℝ)
    (fun ω => by simpa [Real.norm_eq_abs] using clip_abs_le_n n (A ω))

@[simp] lemma clipBCF_apply (n : ℕ) (A : Ω → ℝ) (hA : Continuous A) (ω : Ω) :
    clipBCF n A hA ω = clip n (A ω) := rfl

/-- 交接缝：`rp_transfer` 的被积函数槽位，填的正是 `clipₙ` 的那个乘积。 -/
lemma clip_integrand_eq (n : ℕ) (A : Ω → ℝ) (hA : Continuous A) (θ : Ω → Ω) (ω : Ω) :
    (clipBCF n A hA) (θ ω) * (clipBCF n A hA) ω = clip n (A (θ ω)) * clip n (A ω) := rfl

/-! ## 4. L² 层：`clipₙ ∘ A` 仍在 `L²`（`rp_l2_closure` 的 `hAR`） -/

lemma clip_memLp {μ : Measure Ω} {A : Ω → ℝ} (hA : Continuous A) (hA2 : MemLp A 2 μ)
    (n : ℕ) : MemLp (fun ω => clip n (A ω)) 2 μ :=
  MemLp.of_le hA2 ((continuous_clip n).comp hA).aestronglyMeasurable
    (Eventually.of_forall fun ω => clip_norm_le n (A ω))

end GlueA.OS2

/- 说明（残余，写进签名而不是藏起来）：
   `clip` 修的是**类型/连续性的错配**，不是 μ_a 侧的 RP 本身。μ_a 侧「`clipₙ ∘ A` 满足 RP」
   由 `OS2Wiring.rp_full_cone` 的假设 `hRP : ∀ n i, 0 ≤ ∫ clipₙ(A(θω))·clipₙ(Aω) ∂μ_a` 显式承载：
   在「RP 是对全体 `C_b(Ω)` 的不等式」这一读法下它免费（`clipₙ ∘ A ∈ C_b`）；
   在「RP 只在环多项式代数上」的读法下它是一条**额外假设**（`clipₙ ∘ A` 非线性、不在代数里）。
   两种读法下 `trunc` 都不行（不连续 ⇒ 连 `C_b` 都不是）。 -/
