/-
# Transfer-matrix gap ⇒ exponential decay

This module isolates the *easy half* of the mass-gap equivalence: if the transfer
matrix `T t` decays like `e^{-Δt}` on the orthogonal complement of the vacuum, then
matrix elements decay at the same rate.  The hard half — deriving that decay from a
lower bound on the Hamiltonian — needs the spectral theorem and lives behind the
`SpectralTheoremInput` interface.

Nothing here uses unbounded operators, so everything compiles against plain mathlib.

Main results:
* `SemigroupGap` — the decay hypothesis, bundled;
* `twoPoint_bound` — matrix elements decay with the same rate (Cauchy–Schwarz);
* `decay_at_nat_multiples` — decay survives iteration along the semigroup.
-/
import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.Normed.Operator.Basic
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Pow.Real

noncomputable section

open scoped InnerProductSpace
open RCLike

namespace YangMillsInfra
namespace Spectral

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H]

/-- Transfer-matrix decay on the orthogonal complement of a subspace (the vacuum).

This is the *Euclidean* form of a mass gap and is what the clustering statements use.
It is stated for a one-parameter family `T : ℝ → H →L[ℂ] H` that need not literally be
a semigroup; `decay_at_nat_multiples` below adds the semigroup law. -/
structure SemigroupGap (T : ℝ → H →L[ℂ] H) (V : Submodule ℂ H) (Δ : ℝ) : Prop where
  /-- The rate is positive: this is the mass gap. -/
  Δ_pos : 0 < Δ
  /-- Exponential decay on the vacuum complement. -/
  decay : ∀ t : ℝ, 0 ≤ t → ∀ ψ : H,
    (∀ v ∈ V, inner ℂ ψ v = 0) → ‖T t ψ‖ ≤ Real.exp (-Δ * t) * ‖ψ‖

namespace SemigroupGap

variable {T : ℝ → H →L[ℂ] H} {V : Submodule ℂ H} {Δ : ℝ}

/-- **Matrix elements decay at the rate `Δ`.**  This is the form in which clustering
of the two-point function is usually stated; it follows from the operator-norm bound by
Cauchy–Schwarz. -/
theorem twoPoint_bound (h : SemigroupGap T V Δ) {t : ℝ} (ht : 0 ≤ t)
    (φ ψ : H) (hψ : ∀ v ∈ V, inner ℂ ψ v = 0) :
    ‖inner ℂ φ (T t ψ)‖ ≤ ‖φ‖ * Real.exp (-Δ * t) * ‖ψ‖ := by
  have hdecay : ‖T t ψ‖ ≤ Real.exp (-Δ * t) * ‖ψ‖ := h.decay t ht ψ hψ
  have hexp : 0 ≤ Real.exp (-Δ * t) := Real.exp_nonneg _
  calc ‖inner ℂ φ (T t ψ)‖
      ≤ ‖φ‖ * ‖T t ψ‖ := norm_inner_le_norm _ _
    _ ≤ ‖φ‖ * (Real.exp (-Δ * t) * ‖ψ‖) := by
        exact mul_le_mul_of_nonneg_left hdecay (norm_nonneg _)
    _ = ‖φ‖ * Real.exp (-Δ * t) * ‖ψ‖ := by ring

/-- **Iterated decay.**  If the family is a semigroup and it preserves the vacuum
complement, the decay exponent adds up.  This is what turns a per-step estimate in an
RG argument into a uniform bound over all distances. -/
theorem decay_at_nat_multiples (h : SemigroupGap T V Δ)
    (hmul : ∀ s t : ℝ, 0 ≤ s → 0 ≤ t → T (s + t) = T s * T t)
    (hpres : ∀ t : ℝ, 0 ≤ t → ∀ ψ : H, (∀ v ∈ V, inner ℂ ψ v = 0) →
      ∀ v ∈ V, inner ℂ (T t ψ) v = 0)
    {t : ℝ} (ht : 0 ≤ t) :
    ∀ n : ℕ, ∀ ψ : H, (∀ v ∈ V, inner ℂ ψ v = 0) →
      ‖T ((n : ℝ) * t) ψ‖ ≤ Real.exp (-Δ * ((n : ℝ) * t)) * ‖ψ‖ := by
  intro n
  induction n with
  | zero =>
      intro ψ _
      simp only [Nat.cast_zero, zero_mul, mul_zero]
      rw [Real.exp_zero]
      -- `T 0 = T (0+0) = T 0 * T 0`; we only know the semigroup law for `0+0`.
      -- `T 0` is a contraction by the decay estimate with `t = 0`.
      have hT0 : ‖T 0 ψ‖ ≤ ‖ψ‖ := by
        have := h.decay 0 le_rfl ψ (by assumption)
        simpa using this
      simpa using hT0
  | succ n _ih =>
      intro ψ hψ
      have hn : (0 : ℝ) ≤ (n : ℝ) * t := mul_nonneg (Nat.cast_nonneg n) ht
      have hTpres : ∀ v ∈ V, inner ℂ (T t ψ) v = 0 := hpres t ht ψ hψ
      have hstep : ‖T t ψ‖ ≤ Real.exp (-Δ * t) * ‖ψ‖ := h.decay t ht ψ hψ
      have hcast : (↑(n + 1) : ℝ) * t = (n : ℝ) * t + t := by push_cast; ring
      have hmul' : T ((n : ℝ) * t + t) = T ((n : ℝ) * t) * T t := hmul _ _ hn ht
      rw [hcast, hmul']
      have hnorm : ‖(T ((n : ℝ) * t) * T t) ψ‖ = ‖T ((n : ℝ) * t) (T t ψ)‖ := rfl
      rw [hnorm]
      calc ‖T ((n : ℝ) * t) (T t ψ)‖
          ≤ Real.exp (-Δ * ((n : ℝ) * t)) * ‖T t ψ‖ :=
            h.decay ((n : ℝ) * t) hn (T t ψ) hTpres
        _ ≤ Real.exp (-Δ * ((n : ℝ) * t)) * (Real.exp (-Δ * t) * ‖ψ‖) := by
            gcongr
        _ = Real.exp (-Δ * ((n : ℝ) * t + t)) * ‖ψ‖ := by
            rw [← mul_assoc, ← Real.exp_add]
            congr 1
            ring

end SemigroupGap

end Spectral
end YangMillsInfra
