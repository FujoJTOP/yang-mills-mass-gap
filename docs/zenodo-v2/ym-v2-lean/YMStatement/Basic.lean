import OSforGFF.OS.Axioms
import OSforGFF.Schwinger.Defs

/-!
# Yang--Mills mass-gap statement layer

This file fixes the *statement* we intend to prove, before any attempt at a proof.
The point is to make every field non-vacuous:

* `SatisfiesAllOS` is OSforGFF's bundled OS0--OS4 predicate;
* `NonTrivial` rules out the free (Gaussian) field, which would otherwise satisfy
  the axioms and the gap trivially;
* `ExponentialClustering` is the *rate* form of clustering, not the bare
  "correlations factorize" form (the latter does not imply a spectral gap).

The final Clay statement will quantify over a compact simple gauge group and
assert the existence of a measure satisfying `MassGapStatement`.
-/

noncomputable section

namespace YMStatement

open MeasureTheory OSforGFF

variable {d : ℕ} [Fact (2 ≤ d)]

/-- A Euclidean measure is Gaussian if its generating functional has the
characteristic form of a centred Gaussian.  This is the negation target for
`NonTrivial`. -/
def IsGaussianMeasure (μ : ProbabilityMeasure (FieldConfiguration d)) : Prop :=
  ∃ B : SchwartzTestFunction d →ₗ[ℝ] SchwartzTestFunction d →ₗ[ℝ] ℝ,
    ∀ J : SchwartzTestFunction d,
      GJGeneratingFunctional μ J =
        Complex.exp (-(1 / 2 : ℂ) * ((B J J : ℝ) : ℂ))

/-- The connected four-point function.

**Scalar theories only.**  For a pure gauge theory this is *not* a valid
interaction witness: every gauge-invariant local operator is at least quadratic
in the field strength, so `tr F²` is a Wick square of a Gaussian field and a
free gauge field already has nonzero connected four-point functions.  Use
`RGCharge` below for the gauge case. -/
def ConnectedFourPoint (μ : ProbabilityMeasure (FieldConfiguration d))
    (f₁ f₂ f₃ f₄ : SchwartzTestFunction d) : ℝ :=
  SchwingerFunction μ 4 ![f₁, f₂, f₃, f₄]
    - SchwingerFunction₂ μ f₁ f₂ * SchwingerFunction₂ μ f₃ f₄
    - SchwingerFunction₂ μ f₁ f₃ * SchwingerFunction₂ μ f₂ f₄
    - SchwingerFunction₂ μ f₁ f₄ * SchwingerFunction₂ μ f₂ f₃

/-- Non-triviality: the theory is not a free (Gaussian) field. -/
def NonTrivial (μ : ProbabilityMeasure (FieldConfiguration d)) : Prop :=
  ¬ IsGaussianMeasure μ

/-- A gauge-invariant renormalized-coupling witness.

For a pure gauge theory the correct non-triviality witness is the renormalized
coupling — a ratio of a gauge-invariant vertex to a propagator — not a connected
four-point function (see `ConnectedFourPoint`).  This structure records the
interface without committing to a construction: the `value` is a real number,
and `free_vanishes` is the property that makes the witness non-vacuous. -/
structure RGCharge (μ : ProbabilityMeasure (FieldConfiguration d)) where
  value : ℝ
  free_vanishes : IsGaussianMeasure μ → value = 0

/-- Non-triviality in the gauge-theory sense: some renormalized-coupling witness
is nonzero.  This is the witness form of `NonTrivial`; the structure exists so
that the witness cannot be a vacuous connected four-point function. -/
def GaugeNonTrivial (μ : ProbabilityMeasure (FieldConfiguration d)) : Prop :=
  ∃ c : RGCharge μ, c.value ≠ 0

/-- Uniform exponential clustering of the generating functional.
The rate `m` is what the mass gap is equivalent to; bare clustering (without a
rate) is not enough. -/
def ExponentialClustering (μ : ProbabilityMeasure (FieldConfiguration d)) (m : ℝ) : Prop :=
  m > 0 ∧ ∀ f g : SchwartzTestFunction d, ∃ C : ℝ, 0 ≤ C ∧
    ∀ a : SpaceTime d,
      ‖GJGeneratingFunctional μ (f + g.translate a) -
        GJGeneratingFunctional μ f * GJGeneratingFunctional μ g‖
        ≤ C * Real.exp (-m * ‖a‖)

/-- The Euclidean form of the Clay statement for a fixed gauge group: a
non-Gaussian OS measure with a uniform exponential clustering rate. -/
structure MassGapStatement (μ : ProbabilityMeasure (FieldConfiguration d)) : Prop where
  os : SatisfiesAllOS μ
  nonTrivial : NonTrivial μ
  gap : ∃ m : ℝ, ExponentialClustering μ m

end YMStatement
