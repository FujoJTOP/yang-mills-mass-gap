/-
GL4L2.lean — the L²-truncation closure completing the machine check of GL4.

`rp_transfer` in GL4Transfer.lean handles bounded continuous test functions.
This file closes the gap for the full reflection-positivity cone:

  if the truncations `A_n = A · 1_{|A| ≤ n}` satisfy the RP inequality and
  `A_n → A` in L², then `A` satisfies it.

The proof is the one in REPORT.md §3.2: triangle inequality, Cauchy–Schwarz
(Hölder with p = q = 2), and the measure-preserving isometry `f ↦ f ∘ θ`.
No `sorry`.
-/
import Mathlib.MeasureTheory.Function.L2Space
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.MeasureTheory.Measure.Typeclasses.Probability

open MeasureTheory Filter
open scoped ENNReal Topology

namespace GlueA

variable {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω} [IsProbabilityMeasure μ]

/-- The truncated test function `A · 1_{|A| ≤ n}`. -/
noncomputable def trunc (A : Ω → ℝ) (n : ℕ) (ω : Ω) : ℝ :=
  if |A ω| ≤ (n : ℝ) then A ω else 0

lemma trunc_abs_le (A : Ω → ℝ) (n : ℕ) (ω : Ω) : ‖trunc A n ω‖ ≤ ‖A ω‖ := by
  unfold trunc
  split_ifs with h
  · exact le_rfl
  · simpa using norm_nonneg (A ω)

lemma trunc_tendsto (A : Ω → ℝ) (ω : Ω) :
    Tendsto (fun n : ℕ => trunc A n ω) atTop (𝓝 (A ω)) := by
  obtain ⟨N, hN⟩ := exists_nat_gt (|A ω|)
  have hev : ∀ᶠ n : ℕ in atTop, trunc A n ω = A ω := by
    filter_upwards [Ici_mem_atTop N] with m hm
    unfold trunc
    rw [if_pos]
    exact le_trans (le_of_lt hN) (by exact_mod_cast hm)
  exact tendsto_const_nhds.congr' (hev.mono fun n h => h.symm)

/-- `‖f‖₂ = sqrt ∫ f²` (with respect to `μ`). -/
noncomputable def nrm (μ : Measure Ω) (f : Ω → ℝ) : ℝ :=
  Real.sqrt (∫ ω, f ω ^ 2 ∂μ)

lemma nrm_nonneg (f : Ω → ℝ) : 0 ≤ nrm μ f := Real.sqrt_nonneg _

/-- Cauchy–Schwarz in the form used by GL4: `∫ ‖f‖‖g‖ ≤ ‖f‖₂‖g‖₂`. -/
lemma holder2 {f g : Ω → ℝ} (hf : MemLp f 2 μ) (hg : MemLp g 2 μ) :
    ∫ ω, ‖f ω‖ * ‖g ω‖ ∂μ ≤ nrm μ f * nrm μ g := by
  have hf2 : MemLp f (ENNReal.ofReal (2 : ℝ)) μ := by
    rwa [show ENNReal.ofReal (2 : ℝ) = 2 by norm_num]
  have hg2 : MemLp g (ENNReal.ofReal (2 : ℝ)) μ := by
    rwa [show ENNReal.ofReal (2 : ℝ) = 2 by norm_num]
  have h := integral_mul_norm_le_Lp_mul_Lq (p := (2 : ℝ)) (q := (2 : ℝ))
    Real.HolderConjugate.two_two hf2 hg2
  have hf' : (∫ ω, ‖f ω‖ ^ (2 : ℝ) ∂μ) = ∫ ω, f ω ^ 2 ∂μ := by
    apply integral_congr_ae
    filter_upwards with ω
    rw [Real.norm_eq_abs, Real.rpow_two, sq_abs]
  have hg' : (∫ ω, ‖g ω‖ ^ (2 : ℝ) ∂μ) = ∫ ω, g ω ^ 2 ∂μ := by
    apply integral_congr_ae
    filter_upwards with ω
    rw [Real.norm_eq_abs, Real.rpow_two, sq_abs]
  rw [hf', hg'] at h
  simpa only [nrm, Real.sqrt_eq_rpow] using h

/-- `nrm` is invariant under the measure-preserving involution `θ`. -/
lemma nrm_comp (θ : Ω ≃ᵐ Ω) (hθ : MeasurePreserving (⇑θ) μ μ) (f : Ω → ℝ) :
    nrm μ (fun ω => f (θ ω)) = nrm μ f := by
  unfold nrm
  congr 1
  exact hθ.integral_comp' (fun y => f y ^ 2)

/-- **GL4, full cone (L² closure).**
If `A_n → A` pointwise, `|A_n| ≤ |A|`, each `A_n ∈ L²`, and every `A_n`
satisfies the reflection-positivity inequality, then so does `A`.

This is the step that `rp_transfer` (bounded continuous) leaves open; with it
GL4 is machine-checked for the whole cone, not only the bounded layer. -/
theorem rp_l2_closure
    (θ : Ω ≃ᵐ Ω) (hθ : MeasurePreserving (⇑θ) μ μ)
    (A : Ω → ℝ) (hA : MemLp A 2 μ)
    (AR : ℕ → Ω → ℝ) (hAR : ∀ n, MemLp (AR n) 2 μ)
    (hpt : ∀ ω, Tendsto (fun n => AR n ω) atTop (𝓝 (A ω)))
    (hbound : ∀ n ω, ‖AR n ω‖ ≤ ‖A ω‖)
    (hRP : ∀ n, 0 ≤ ∫ ω, AR n (θ ω) * AR n ω ∂μ) :
    0 ≤ ∫ ω, A (θ ω) * A ω ∂μ := by
  haveI : ENNReal.HolderTriple (2 : ℝ≥0∞) 2 1 :=
    ⟨by simpa using ENNReal.inv_two_add_inv_two⟩
  -- (1) L² convergence of the differences
  have hL2 : Tendsto (fun n => ∫ ω, (AR n ω - A ω) ^ 2 ∂μ) atTop (𝓝 0) := by
    have hmeas : ∀ n, AEStronglyMeasurable (fun ω => (AR n ω - A ω) ^ 2) μ :=
      fun n => ((hAR n).1.sub hA.1).pow 2
    have hint : Integrable (fun ω => 4 * A ω ^ 2) μ := by
      have h := hA.integrable_sq
      simpa [mul_comm] using h.const_mul 4
    have hdom : ∀ n, ∀ᵐ ω ∂μ, ‖(AR n ω - A ω) ^ 2‖ ≤ 4 * A ω ^ 2 := by
      intro n
      filter_upwards with ω
      have h1 : ‖AR n ω - A ω‖ ≤ 2 * ‖A ω‖ := by
        calc ‖AR n ω - A ω‖ ≤ ‖AR n ω‖ + ‖A ω‖ := norm_sub_le _ _
          _ ≤ ‖A ω‖ + ‖A ω‖ := by gcongr; exact hbound n ω
          _ = 2 * ‖A ω‖ := by ring
      have h2 : ‖(AR n ω - A ω) ^ 2‖ = ‖AR n ω - A ω‖ ^ 2 := by
        rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _), Real.norm_eq_abs, sq_abs]
      have hA2 : ‖A ω‖ ^ 2 = A ω ^ 2 := by rw [Real.norm_eq_abs, sq_abs]
      rw [h2, ← hA2]
      nlinarith [h1, norm_nonneg (AR n ω - A ω), norm_nonneg (A ω)]
    have hlim : ∀ᵐ ω ∂μ, Tendsto (fun n => (AR n ω - A ω) ^ 2) atTop (𝓝 0) := by
      filter_upwards with ω
      have hAconst : Tendsto (fun _ : ℕ => A ω) atTop (𝓝 (A ω)) := tendsto_const_nhds
      simpa using (hpt ω).sub hAconst |>.pow 2
    have h := tendsto_integral_of_dominated_convergence (μ := μ)
      (F := fun n ω => (AR n ω - A ω) ^ 2) (f := fun _ => (0 : ℝ))
      (bound := fun ω => 4 * A ω ^ 2) hmeas hint hdom hlim
    simpa using h
  -- (2) the sqrt norms tend to zero
  have hnrm : Tendsto (fun n => nrm μ (fun ω => AR n ω - A ω)) atTop (𝓝 0) := by
    have hcomp : Tendsto (fun n => Real.sqrt (∫ ω, (AR n ω - A ω) ^ 2 ∂μ)) atTop
        (𝓝 (Real.sqrt 0)) := (Real.continuous_sqrt.tendsto 0).comp hL2
    simpa [nrm] using hcomp
  -- (3) uniform bound on the truncation norms
  have hnrm_le : ∀ n, nrm μ (AR n) ≤ nrm μ A := by
    intro n
    apply Real.sqrt_le_sqrt
    apply integral_mono
    · exact (hAR n).integrable_sq
    · exact hA.integrable_sq
    · intro ω
      have h := hbound n ω
      rw [Real.norm_eq_abs, Real.norm_eq_abs] at h
      have h1 : |AR n ω| ^ 2 = AR n ω ^ 2 := sq_abs _
      have h2 : |A ω| ^ 2 = A ω ^ 2 := sq_abs _
      nlinarith [h, abs_nonneg (AR n ω), abs_nonneg (A ω), h1, h2]
  -- (4) the key estimate
  have heach : ∀ n, |∫ ω, AR n (θ ω) * AR n ω ∂μ - ∫ ω, A (θ ω) * A ω ∂μ|
      ≤ 2 * nrm μ A * nrm μ (fun ω => AR n ω - A ω) := by
    intro n
    have hBθ : MemLp (fun ω => AR n (θ ω)) 2 μ := (hAR n).comp_measurePreserving hθ
    have hB : MemLp (AR n) 2 μ := hAR n
    have hAθ : MemLp (fun ω => A (θ ω)) 2 μ := hA.comp_measurePreserving hθ
    have hPn : Integrable (fun ω => AR n (θ ω) * AR n ω) μ := hBθ.integrable_mul hB
    have hP : Integrable (fun ω => A (θ ω) * A ω) μ := hAθ.integrable_mul hA
    have hptw : ∀ ω, ‖AR n (θ ω) * AR n ω - A (θ ω) * A ω‖
        ≤ ‖AR n (θ ω)‖ * ‖AR n ω - A ω‖ + ‖AR n (θ ω) - A (θ ω)‖ * ‖A ω‖ := by
      intro ω
      have hdec : AR n (θ ω) * AR n ω - A (θ ω) * A ω
          = AR n (θ ω) * (AR n ω - A ω) + (AR n (θ ω) - A (θ ω)) * A ω := by ring
      calc ‖AR n (θ ω) * AR n ω - A (θ ω) * A ω‖
          = |AR n (θ ω) * AR n ω - A (θ ω) * A ω| := Real.norm_eq_abs _
        _ = |AR n (θ ω) * (AR n ω - A ω) + (AR n (θ ω) - A (θ ω)) * A ω| := by rw [hdec]
        _ ≤ |AR n (θ ω) * (AR n ω - A ω)| + |(AR n (θ ω) - A (θ ω)) * A ω| := abs_add_le _ _
        _ = |AR n (θ ω)| * |AR n ω - A ω| + |AR n (θ ω) - A (θ ω)| * |A ω| := by
              rw [abs_mul, abs_mul]
        _ = ‖AR n (θ ω)‖ * ‖AR n ω - A ω‖ + ‖AR n (θ ω) - A (θ ω)‖ * ‖A ω‖ := by
              rw [Real.norm_eq_abs, Real.norm_eq_abs, Real.norm_eq_abs, Real.norm_eq_abs]
    calc |∫ ω, AR n (θ ω) * AR n ω ∂μ - ∫ ω, A (θ ω) * A ω ∂μ|
        = |∫ ω, (AR n (θ ω) * AR n ω - A (θ ω) * A ω) ∂μ| := by
            rw [integral_sub hPn hP]
      _ ≤ ∫ ω, ‖AR n (θ ω) * AR n ω - A (θ ω) * A ω‖ ∂μ :=
            abs_integral_le_integral_abs
      _ ≤ ∫ ω, (‖AR n (θ ω)‖ * ‖AR n ω - A ω‖ + ‖AR n (θ ω) - A (θ ω)‖ * ‖A ω‖) ∂μ := by
            apply integral_mono
            · exact (hPn.sub hP).norm
            · exact ((hBθ.norm).integrable_mul ((hB.sub hA).norm)).add
                (((hBθ.sub hAθ).norm).integrable_mul hA.norm)
            · intro ω; exact hptw ω
      _ = (∫ ω, ‖AR n (θ ω)‖ * ‖AR n ω - A ω‖ ∂μ)
            + (∫ ω, ‖AR n (θ ω) - A (θ ω)‖ * ‖A ω‖ ∂μ) := by
            rw [integral_add]
            · exact (hBθ.norm).integrable_mul ((hB.sub hA).norm)
            · exact ((hBθ.sub hAθ).norm).integrable_mul hA.norm
      _ ≤ nrm μ (fun ω => AR n (θ ω)) * nrm μ (fun ω => AR n ω - A ω)
            + nrm μ (fun ω => AR n (θ ω) - A (θ ω)) * nrm μ A := by
            gcongr
            · exact holder2 hBθ (hB.sub hA)
            · exact holder2 (hBθ.sub hAθ) hA
      _ = nrm μ (AR n) * nrm μ (fun ω => AR n ω - A ω)
            + nrm μ (fun ω => AR n ω - A ω) * nrm μ A := by
            rw [nrm_comp θ hθ, nrm_comp (f := fun ω => AR n ω - A ω) θ hθ]
      _ ≤ nrm μ A * nrm μ (fun ω => AR n ω - A ω)
            + nrm μ (fun ω => AR n ω - A ω) * nrm μ A := by
            exact add_le_add
              (mul_le_mul_of_nonneg_right (hnrm_le n)
                (nrm_nonneg (μ := μ) (f := fun ω => AR n ω - A ω))) le_rfl
      _ = 2 * nrm μ A * nrm μ (fun ω => AR n ω - A ω) := by ring
  -- (5) squeeze and conclude
  have hconv : Tendsto (fun n => ∫ ω, AR n (θ ω) * AR n ω ∂μ) atTop
      (𝓝 (∫ ω, A (θ ω) * A ω ∂μ)) := by
    rw [tendsto_iff_dist_tendsto_zero]
    refine squeeze_zero
      (f := fun n => dist (∫ ω, AR n (θ ω) * AR n ω ∂μ) (∫ ω, A (θ ω) * A ω ∂μ))
      (g := fun n => 2 * nrm μ A * nrm μ (fun ω => AR n ω - A ω))
      (fun n => dist_nonneg) (fun n => ?_) ?_
    · rw [Real.dist_eq]; exact heach n
    · have h2 : Tendsto (fun n => 2 * nrm μ A * nrm μ (fun ω => AR n ω - A ω)) atTop
          (𝓝 (2 * nrm μ A * 0)) := tendsto_const_nhds.mul hnrm
      simpa using h2
  exact ge_of_tendsto hconv (Eventually.of_forall hRP)

/-- **GL4, full cone, concrete form.**  The hypotheses are exactly the ones the
truncation `A_n = A · 1_{|A| ≤ n}` satisfies: `A_n → A` pointwise, `|A_n| ≤ |A|`,
and `A_n ∈ L²`.  So reflection positivity for all truncations gives it for `A`. -/
theorem rp_l2_closure_trunc
    (θ : Ω ≃ᵐ Ω) (hθ : MeasurePreserving (⇑θ) μ μ)
    (A : Ω → ℝ) (hA : MemLp A 2 μ)
    (hAR : ∀ n, MemLp (trunc A n) 2 μ)
    (hRP : ∀ n, 0 ≤ ∫ ω, trunc A n (θ ω) * trunc A n ω ∂μ) :
    0 ≤ ∫ ω, A (θ ω) * A ω ∂μ :=
  rp_l2_closure θ hθ A hA (fun n => trunc A n)
    hAR (trunc_tendsto A) (trunc_abs_le A) hRP

end GlueA

/-- info: 'GlueA.rp_l2_closure_trunc' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueA.rp_l2_closure_trunc

/-- info: 'GlueA.rp_l2_closure' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueA.rp_l2_closure

/-- info: 'GlueA.holder2' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueA.holder2
