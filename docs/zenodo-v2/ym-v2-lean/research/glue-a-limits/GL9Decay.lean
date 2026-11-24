/-
GL9Decay.lean — the spectral-decay core of glue lemma GL9.

GL9 (see REPORT.md §4): a uniform transfer-matrix gap implies volume-uniform
exponential clustering.  The only non-elementary step is the following purely
algebraic lemma, valid in any real inner-product space (finite- or
infinite-dimensional):

  if  T Ω = Ω,  ‖Ω‖ = 1,  T preserves Ω^⊥,  and  ‖T x‖ ≤ λ‖x‖ on Ω^⊥,
  then  |⟪A, Tⁿ B⟫ − ⟪Ω,B⟫⟪A,Ω⟫| ≤ λⁿ ‖A‖ ‖B‖.

This is the spectral decomposition with the vacuum removed (the operator
`offVac` is the projection off the Ω direction).  In GL9 one substitutes
λ = exp(−aδ), n a = t, and reads off |c(t)| ≤ exp(−δ t) ‖A‖‖B‖.

No `sorry`, no extra axioms: the file compiles against the pinned Mathlib and
`#print axioms` reports only [propext, Classical.choice, Quot.sound].
-/
import Mathlib.Analysis.InnerProductSpace.Basic

namespace GlueA

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- `x` minus its component along the vacuum direction `Ω`. -/
noncomputable def offVac (Ω x : E) : E := x - (inner ℝ Ω x) • Ω

/-- `T` preserves the orthogonal complement of `Ω`
(the consequence of self-adjointness used by the spectral argument). -/
def PreservesPerp (T : E →ₗ[ℝ] E) (Ω : E) : Prop :=
  ∀ x, inner ℝ Ω x = 0 → inner ℝ Ω (T x) = 0

/-- The transfer-matrix gap hypothesis, stated on `Ω^⊥`. -/
def GapOnPerp (T : E →ₗ[ℝ] E) (Ω : E) (lam : ℝ) : Prop :=
  ∀ x, inner ℝ Ω x = 0 → ‖T x‖ ≤ lam * ‖x‖

/-- `offVac Ω x` is orthogonal to `Ω` when `Ω` is a unit vector. -/
lemma inner_offVac (hΩ : ‖Ω‖ = 1) (x : E) : inner ℝ Ω (offVac Ω x) = 0 := by
  simp only [offVac, inner_sub_right, inner_smul_right, real_inner_self_eq_norm_sq, hΩ,
    one_pow, mul_one, sub_self]

/-- The projection off `Ω` does not increase the norm. -/
lemma norm_offVac_le (hΩ : ‖Ω‖ = 1) (x : E) : ‖offVac Ω x‖ ≤ ‖x‖ := by
  have hsq : ‖offVac Ω x‖ ^ 2 = ‖x‖ ^ 2 - (inner ℝ Ω x) ^ 2 := by
    rw [offVac, norm_sub_sq_real, norm_smul, Real.norm_eq_abs, hΩ, mul_one,
      inner_smul_right, sq_abs, real_inner_comm x Ω]
    ring
  have hnn : 0 ≤ ‖offVac Ω x‖ := norm_nonneg _
  nlinarith [hsq, hnn, norm_nonneg x, sq_nonneg (inner ℝ Ω x)]

/-- Decompose `Tⁿ B` into its vacuum part plus the recursively shifted remainder. -/
lemma iterate_decomp (T : E →ₗ[ℝ] E) (hTΩ : T Ω = Ω) (B : E) (n : ℕ) :
    (T : E → E)^[n] B = (inner ℝ Ω B) • Ω + (T : E → E)^[n] (offVac Ω B) := by
  induction n with
  | zero => simp [offVac]
  | succ n ih =>
      rw [Function.iterate_succ_apply', ih]
      simp only [map_add, map_smul, hTΩ]
      rw [Function.iterate_succ_apply']

/-- Membership in `Ω^⊥` is preserved under all iterates. -/
lemma perp_iterate (hperp : PreservesPerp T Ω) {y : E} (hy : inner ℝ Ω y = 0) (n : ℕ) :
    inner ℝ Ω ((T : E → E)^[n] y) = 0 := by
  induction n with
  | zero => simpa using hy
  | succ n ih =>
      rw [Function.iterate_succ_apply']
      exact hperp _ ih

/-- Geometric decay of the iterates on `Ω^⊥`: `‖Tⁿ y‖ ≤ λⁿ‖y‖`. -/
lemma iterate_norm_le (hperp : PreservesPerp T Ω) (hgap : GapOnPerp T Ω lam)
    (hlam : 0 ≤ lam) {y : E} (hy : inner ℝ Ω y = 0) (n : ℕ) :
    ‖(T : E → E)^[n] y‖ ≤ lam ^ n * ‖y‖ := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [Function.iterate_succ_apply']
      calc ‖T ((T : E → E)^[n] y)‖ ≤ lam * ‖(T : E → E)^[n] y‖ :=
            hgap _ (perp_iterate hperp hy n)
        _ ≤ lam * (lam ^ n * ‖y‖) := mul_le_mul_of_nonneg_left ih hlam
        _ = lam ^ (n + 1) * ‖y‖ := by rw [pow_succ]; ring

/-- **GL9 core.**  With the vacuum removed, the `n`-step correlation is bounded
by `λⁿ‖A‖‖B‖`. -/
theorem spectral_decay (hΩ : ‖Ω‖ = 1) (hTΩ : T Ω = Ω) (hperp : PreservesPerp T Ω)
    (hgap : GapOnPerp T Ω lam) (hlam : 0 ≤ lam) (A B : E) (n : ℕ) :
    |inner ℝ A ((T : E → E)^[n] B) - inner ℝ Ω B * inner ℝ A Ω|
      ≤ lam ^ n * ‖A‖ * ‖B‖ := by
  have hkey : inner ℝ A ((T : E → E)^[n] B) - inner ℝ Ω B * inner ℝ A Ω
      = inner ℝ A ((T : E → E)^[n] (offVac Ω B)) := by
    rw [iterate_decomp (T := T) hTΩ B n, inner_add_right, inner_smul_right]
    ring
  rw [hkey]
  calc |inner ℝ A ((T : E → E)^[n] (offVac Ω B))|
        ≤ ‖A‖ * ‖(T : E → E)^[n] (offVac Ω B)‖ := abs_real_inner_le_norm A _
    _ ≤ ‖A‖ * (lam ^ n * ‖offVac Ω B‖) := by
          gcongr
          exact iterate_norm_le hperp hgap hlam (inner_offVac hΩ B) n
    _ ≤ ‖A‖ * (lam ^ n * ‖B‖) := by
          gcongr
          exact norm_offVac_le hΩ B
    _ = lam ^ n * ‖A‖ * ‖B‖ := by ring

end GlueA

-- Axiom audit: only the three standard axioms, no `sorryAx`.
#print axioms GlueA.spectral_decay
