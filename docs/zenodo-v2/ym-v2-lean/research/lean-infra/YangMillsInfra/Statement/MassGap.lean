/-
# The mass-gap statement layer

The statement of the mass gap that the rest of the program targets is *not* the bare
clustering axiom OS4.  As hexagonmath's reconstruction (arXiv 2611.00492) shows, there
is an OS model satisfying all 16 certified inputs and OS4 clustering whose Hamiltonian
has no spectral gap.  What is needed is clustering **with a rate**, uniformly in the
volume:

    ∃ m > 0, ∀ f g, ∃ C, ∀ a, |Z[f + T_a g] − Z[f] Z[g]| ≤ C e^{−m‖a‖}.

This file:
* bundles that statement (`MassGapStatement`);
* proves that the rate form implies bare clustering (`exponentialClustering_implies_clustering`),
  so the rate form really is a strengthening of OS4 and not an unrelated condition.

The converse (bare clustering does **not** imply the rate form) is a counterexample
result, not formalized here; see REPORT.md §"What is deliberately not formalized".
-/
import YangMillsInfra.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Exp

noncomputable section

namespace YangMillsInfra

variable {d : ℕ}

/-- The Euclidean mass-gap statement for a theory `T`, relative to a supplied
non-triviality predicate.

`NonTrivial` is an explicit parameter rather than a definition here because the
nontriviality line owns that notion (free-field detection, connected four-point
witness, ...).  The statement layer only needs to know that the theory is not
Gaussian. -/
structure MassGapStatement (T : EuclideanTheory d) (τ : TranslationAction d)
    (NonTrivial : EuclideanTheory d → Prop) : Prop where
  nonTrivial : NonTrivial T
  gap : ∃ m : ℝ, ExponentialClustering T τ m

/-- **The rate form refines OS4.**  Exponential clustering with rate `m > 0` implies
bare clustering: for every `ε > 0` there is a radius beyond which the connected
generating functional is smaller than `ε`.

This is the formal content of the claim that the mass gap implies the OS4 clustering
axiom in `OSforGFF.OS.Axioms.OS4_Clustering`; it also shows the rate form is not
vacuous in the `ε`-direction. -/
theorem exponentialClustering_implies_clustering {T : EuclideanTheory d}
    {τ : TranslationAction d} {m : ℝ} (h : ExponentialClustering T τ m) :
    Clustering T τ := by
  obtain ⟨hm, hb⟩ := h
  intro f g ε hε
  obtain ⟨C, hC, hbound⟩ := hb f g
  rcases eq_or_lt_of_le hC with hC0 | hCpos
  · -- `C = 0`: the difference vanishes identically, so any radius works.
    refine ⟨1, one_pos, fun a _ => ?_⟩
    have := hbound a
    rw [← hC0] at this
    simp only [zero_mul] at this
    linarith
  · -- `C > 0`: choose `R` large enough that `C e^{−mR} < ε`.
    refine ⟨max 1 (Real.log (C / ε) / m), lt_of_lt_of_le one_pos (le_max_left _ _),
      fun a ha => ?_⟩
    have hR : Real.log (C / ε) / m ≤ ‖a‖ := by
      have h1 : Real.log (C / ε) / m ≤ max 1 (Real.log (C / ε) / m) := le_max_right _ _
      exact le_trans h1 (le_of_lt ha)
    have hexp : Real.exp (-m * (Real.log (C / ε) / m)) = ε / C := by
      have hsimp : -m * (Real.log (C / ε) / m) = -Real.log (C / ε) := by
        field_simp
      rw [hsimp, Real.exp_neg, Real.exp_log (div_pos hCpos hε)]
      field_simp
    have hdecay : Real.exp (-m * ‖a‖) < ε / C := by
      have hlt : Real.log (C / ε) / m < ‖a‖ := lt_of_le_of_lt (le_max_right _ _) ha
      have hneg : -m * ‖a‖ < -m * (Real.log (C / ε) / m) := by nlinarith
      calc Real.exp (-m * ‖a‖) < Real.exp (-m * (Real.log (C / ε) / m)) :=
            Real.exp_lt_exp.mpr hneg
        _ = ε / C := hexp
    calc ‖T.Z (f + τ.translate a g) - T.Z f * T.Z g‖
        ≤ C * Real.exp (-m * ‖a‖) := hbound a
      _ < C * (ε / C) := mul_lt_mul_of_pos_left hdecay hCpos
      _ = ε := by
          have hCne : C ≠ 0 := ne_of_gt hCpos
          field_simp

end YangMillsInfra
