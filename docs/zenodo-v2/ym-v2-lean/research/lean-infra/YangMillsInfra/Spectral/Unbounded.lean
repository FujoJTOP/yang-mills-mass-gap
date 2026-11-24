/-
# Unbounded self-adjoint operators and the quadratic-form gap

mathlib (v4.33.0-rc1) has `LinearPMap` — partially defined linear maps — together with
`LinearPMap.adjoint`, the `Star` instance, `IsSelfAdjoint`, and
`IsSelfAdjoint.isClosed`.  What it does **not** have is any of:

* the spectral theorem for unbounded self-adjoint operators;
* projection-valued measures / spectral families;
* the measurable functional calculus for unbounded operators;
* Stone's theorem / the generator of a strongly continuous unitary group.

This file therefore takes the following design decision, which is the one the rest of
the mass-gap program should use:

> The **mass gap** is stated as a lower bound on the quadratic form (numerical range)
> of the Hamiltonian on the orthogonal complement of the vacuum.  This is equivalent to
> `spec (H) ∩ (0, Δ) = ∅` in the presence of the spectral theorem, but it is a
> first-order statement about pairings and needs no spectral measure.

The equivalence with the spectral-theoretic gap is isolated as the hypothesis bundle
`SpectralTheoremInput` (see `Spectral/GapEquivalence.lean` in the design), so downstream
work can be conditional on it explicitly.

Main definitions:
* `SelfAdjointOperator` — a bundled self-adjoint `LinearPMap`;
* `HasGapAbove` — the quadratic-form gap;
* `RayleighQuotient` and its basic properties.
-/
import Mathlib.Analysis.InnerProductSpace.LinearPMap
import Mathlib.Analysis.InnerProductSpace.Adjoint
import Mathlib.Data.Real.Basic

noncomputable section

open scoped InnerProductSpace
open RCLike

namespace YangMillsInfra
namespace Spectral

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℂ H] [CompleteSpace H]

/-- A self-adjoint (possibly unbounded) operator on a complex Hilbert space, bundled with
the density of its domain (which mathlib derives from self-adjointness, but which we
record for convenience). -/
structure SelfAdjointOperator (H : Type*) [NormedAddCommGroup H] [InnerProductSpace ℂ H]
    [CompleteSpace H] where
  /-- The underlying partially defined operator. -/
  op : H →ₗ.[ℂ] H
  /-- Self-adjointness, in the sense of `LinearPMap.adjoint`. -/
  selfAdjoint : IsSelfAdjoint op

namespace SelfAdjointOperator

variable (A : SelfAdjointOperator H)

/-- The domain of `A`, as a submodule of `H`. -/
def domain : Submodule ℂ H := A.op.domain

/-- The domain is dense. -/
theorem dense_domain : Dense (A.domain : Set H) :=
  A.selfAdjoint.dense_domain

/-- Application of `A` to a vector in its domain. -/
def apply (ψ : A.domain) : H := A.op ψ

/-- The graph of `A` is closed. -/
theorem isClosed : A.op.IsClosed := A.selfAdjoint.isClosed

/-- The Rayleigh quotient of `A` at a nonzero vector of its domain.  For a self-adjoint
operator the numerator is real, so this is a real number; it is the quantity whose
infimum over the vacuum complement is the bottom of the spectrum there. -/
def rayleigh (ψ : A.domain) : ℝ :=
  (inner ℂ (ψ : H) (A.op ψ)).re / ‖(ψ : H)‖ ^ 2

/-- The quadratic-form gap above a subspace `V` (intended: the span of the vacuum).

`HasGapAbove A V Δ` says `A ≥ Δ` as a quadratic form on `Vᗮ`, restricted to the domain.
Non-vacuity: for the free massive field, `A = √(-Δ + m²)` and `V = ℂ·Ω`, the statement
holds with `Δ = m`; for the massless field it fails for every `Δ > 0`. -/
def HasGapAbove (V : Submodule ℂ H) (Δ : ℝ) : Prop :=
  0 < Δ ∧
    ∀ (ψ : A.domain), (ψ : H) ≠ 0 → (∀ v ∈ V, inner ℂ (ψ : H) v = 0) →
      Δ * ‖(ψ : H)‖ ^ 2 ≤ (inner ℂ (ψ : H) (A.op ψ)).re

/-- `HasGapAbove` restricted to a single vacuum vector. -/
def HasGapAboveVacuum (Ω : H) (Δ : ℝ) : Prop :=
  A.HasGapAbove (Submodule.span ℂ {Ω}) Δ

/-- The quadratic form of a self-adjoint operator is real on its domain. -/
theorem inner_apply_im_zero (ψ : A.domain) : (inner ℂ (ψ : H) (A.op ψ)).im = 0 := by
  -- `A` is self-adjoint, so on the domain `⟪Aψ, ψ⟫ = ⟪ψ, Aψ⟫`.
  have hstar : LinearPMap.adjoint A.op = A.op :=
    LinearPMap.isSelfAdjoint_def.mp A.selfAdjoint
  have hsym : inner ℂ (A.op ψ) (ψ : H) = inner ℂ (ψ : H) (A.op ψ) := by
    have hform := LinearPMap.adjoint_isFormalAdjoint (T := A.op) A.dense_domain
    rw [hstar] at hform
    exact hform ψ ψ
  have hconj : inner ℂ (ψ : H) (A.op ψ)
      = starRingEnd ℂ (inner ℂ (ψ : H) (A.op ψ)) := by
    calc inner ℂ (ψ : H) (A.op ψ)
        = starRingEnd ℂ (inner ℂ (A.op ψ) (ψ : H)) :=
          (inner_conj_symm (𝕜 := ℂ) (ψ : H) (A.op ψ)).symm
      _ = starRingEnd ℂ (inner ℂ (ψ : H) (A.op ψ)) := by rw [hsym]
  exact Complex.conj_eq_iff_im.mp hconj.symm

/-- The gap is monotone in `Δ` in the expected way: a larger `Δ` is a stronger hypothesis. -/
theorem hasGapAbove_mono {V : Submodule ℂ H} {Δ Δ' : ℝ} (h : A.HasGapAbove V Δ)
    (hΔ' : 0 < Δ') (hΔ : Δ' ≤ Δ) : A.HasGapAbove V Δ' := by
  refine ⟨hΔ', ?_⟩
  intro ψ hψ horth
  have := h.2 ψ hψ horth
  have hnorm : 0 ≤ ‖(ψ : H)‖ ^ 2 := sq_nonneg _
  nlinarith

end SelfAdjointOperator

end Spectral
end YangMillsInfra
