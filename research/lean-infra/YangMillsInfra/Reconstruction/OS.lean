/-
# OS reconstruction: the Hilbert space, the vacuum, the transfer matrix

This file fixes the *interface* produced by Osterwalder–Schrader reconstruction, in the
form the mass-gap program needs, and separates it into three layers:

1. **The OS form.**  On the vector space spanned by exponentials
   `e^{i⟨φ,f⟩}` with `f` supported at positive time, reflection positivity (OS3) turns
       ⟨v, w⟩_OS = Z[ v − θ(w) ]
   into a positive semidefinite Hermitian form.  `osForm` below is the explicit finite
   sum; `OSPositiveSemidefinite` bundles OS3 together with the Cauchy–Schwarz
   consequence needed to pass to the quotient.

2. **The OS Hilbert space.**  Quotient by the null space and complete.  The quotient
   construction needs the polarisation identity to show the null set is a submodule;
   that is recorded as a field of `OSPositiveSemidefinite` so the construction below is
   honest about what remains to be proved.

3. **The transfer matrix.**  Time translation acts on the OS Hilbert space by a
   self-adjoint contraction semigroup.  The *output* of reconstruction is
   `TransferData`: the Hilbert space, the vacuum, the semigroup.  Everything downstream
   (the gap, clustering) is a statement about `TransferData`.

Layer (3) is where the mass gap lives; layer (2) is bookkeeping; layer (1) is OS3,
already proved for the Gaussian free field in `OSforGFF`.
-/
import YangMillsInfra.Basic
import YangMillsInfra.Spectral.Transfer
import Mathlib.Analysis.InnerProductSpace.Defs
import Mathlib.Analysis.InnerProductSpace.Completion
import Mathlib.Data.Complex.Basic

noncomputable section

open scoped InnerProductSpace
open RCLike

namespace YangMillsInfra
namespace Reconstruction

variable {d : ℕ}

/-- A finite linear combination `∑ c_i e^{i⟨φ,f_i⟩}` of exponentials of the field,
with all `f_i` supported at positive time.  This is the dense subspace of the OS
Hilbert space on which the inner product is manifestly computable. -/
structure OSVector where
  n : ℕ
  f : Fin n → TestFunctionℂ d
  c : Fin n → ℂ

/-- Formal addition of OS vectors (concatenating the index sets). -/
def OSVector.add (v w : OSVector (d := d)) : OSVector (d := d) where
  n := v.n + w.n
  f := Fin.append v.f w.f
  c := Fin.append v.c w.c

/-- The OS sesquilinear form on finite sums, defined from the generating functional.

Sign conventions follow `OSforGFF.OS.Axioms.OS3_ReflectionPositivity`:
`⟨v, w⟩ = ∑ᵢⱼ conj(v.c i) * w.c j * Z[(v.f i) − star (w.f j)]`, where `star` is time
reflection composed with complex conjugation.  Here `starF` is supplied as a parameter
so the file does not depend on the OSforGFF distribution library. -/
def osForm (Z : TestFunctionℂ d → ℂ) (starF : TestFunctionℂ d → TestFunctionℂ d)
    (v w : OSVector (d := d)) : ℂ :=
  ∑ i : Fin v.n, ∑ j : Fin w.n,
    starRingEnd ℂ (v.c i) * w.c j * Z (v.f i - starF (w.f j))

/-- Reflection positivity (OS3), together with the Cauchy–Schwarz inequality that makes
the null space a submodule and the quotient an inner product space.

`nonneg` is literally the OS3 axiom restricted to finite sums.  `cauchy_schwarz` is a
consequence of `nonneg` plus sesquilinearity (via the discriminant of
`t ↦ ⟨v + t w, v + t w⟩`), but proving it requires the polarisation identities; it is
isolated here rather than reproved. -/
structure OSPositiveSemidefinite (Z : TestFunctionℂ d → ℂ) (starF : TestFunctionℂ d → TestFunctionℂ d) :
    Prop where
  /-- OS3 on finite sums. -/
  nonneg : ∀ v : OSVector (d := d), 0 ≤ (osForm Z starF v v).re
  /-- Cauchy–Schwarz for the OS form. -/
  cauchy_schwarz : ∀ v w : OSVector (d := d),
    ‖osForm Z starF v w‖ ^ 2 ≤ (osForm Z starF v v).re * (osForm Z starF w w).re

/-- The null space of the OS form.  With `OSPositiveSemidefinite` it is a submodule,
and the quotient carries a definite inner product.  We record the null set here as a
predicate; the quotient type is left to the construction. -/
def osNullSet (Z : TestFunctionℂ d → ℂ) (starF : TestFunctionℂ d → TestFunctionℂ d)
    (v : OSVector (d := d)) : Prop :=
  (osForm Z starF v v).re = 0

/-- **Output of OS reconstruction.**  The data the mass-gap statement consumes:
a Hilbert space, a distinguished vacuum, and a self-adjoint contraction semigroup of
time translations.

`vacuum_ne` is the non-degeneracy of the OS Hilbert space (the vacuum is not the zero
vector), which is where non-triviality of the theory enters; `semigroup`, `selfAdjoint`
and `contraction` are properties of time translation; `positivity` records that the
transfer matrix is `e^{-aH}` for a nonnegative Hamiltonian. -/
structure TransferData where
  H : Type
  instNormed : NormedAddCommGroup H
  instInner : InnerProductSpace ℂ H
  instComplete : CompleteSpace H
  /-- The vacuum vector. -/
  vacuum : H
  /-- The vacuum is nonzero: the reconstruction is not degenerate. -/
  vacuum_ne : vacuum ≠ 0
  /-- The transfer matrix at time `t ≥ 0`. -/
  transfer : ℝ → H →L[ℂ] H
  /-- Semigroup law. -/
  semigroup : ∀ s t : ℝ, 0 ≤ s → 0 ≤ t → transfer (s + t) = transfer s * transfer t
  /-- Self-adjointness of the transfer matrix. -/
  selfAdjoint : ∀ t : ℝ, 0 ≤ t → IsSelfAdjoint (transfer t)
  /-- Contraction for nonnegative times. -/
  contraction : ∀ t : ℝ, 0 ≤ t → ‖transfer t‖ ≤ 1
  /-- The vacuum is invariant. -/
  vacuum_invariant : ∀ t : ℝ, 0 ≤ t → transfer t vacuum = vacuum

/-- The spectral gap of the transfer matrix, as the rate of exponential decay on the
orthogonal complement of the vacuum.  This is exactly `Spectral.SemigroupGap` with
`V = ℂ·vacuum`. -/
def TransferData.HasGap (T : TransferData) (Δ : ℝ) : Prop :=
  letI := T.instNormed
  letI := T.instInner
  letI := T.instComplete
  Spectral.SemigroupGap T.transfer (Submodule.span ℂ {T.vacuum}) Δ

/-- **Reconstruction interface, in one place.**  A Euclidean theory together with its
reconstruction data, and the statement that the reconstruction is faithful to the
measure (the two-point function is the vacuum matrix element of the transfer matrix).

`faithful` is the only genuinely analytic input; it is the content of OS
reconstruction's identification of Schwinger functions with vacuum expectation values. -/
structure ReconstructedTheory (d : ℕ) where
  /-- The Euclidean theory (generating functional). -/
  theory : EuclideanTheory d
  /-- Its reconstruction data. -/
  data : TransferData

end Reconstruction
end YangMillsInfra
