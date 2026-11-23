/-
# The Wilson lattice measure

For a finite lattice `Λ` and a compact gauge group `G`, the Wilson measure at inverse
coupling `β` is the probability measure on link configurations

    dμ_β(U) = Z_β⁻¹ exp(−S_β(U)) ∏_{l ∈ Λ.L} dg_l ,

where `dg` is normalized Haar measure on `G` and `S_β` is the Wilson action of
`Gauge/Lattice.lean`.

The finite-volume measure is where the whole construction starts, so it is worth having
its exact type.  The only mathematically substantive input is bi-invariance of Haar
measure on a compact group (mathlib has left-invariance; compact groups are
unimodular, but the right-invariance instance needed for the gauge action is not
obviously available) — gap G-3, recorded as a hypothesis field below.

Nothing here assumes the measure exists: `wilsonMeasure` is *defined* as a normalized
`withDensity`, and its total mass is computed from the normalization hypotheses.
-/
import YangMillsInfra.Gauge.Lattice
import Mathlib.MeasureTheory.Measure.Haar.Basic
import Mathlib.Topology.Sets.Compacts
import Mathlib.MeasureTheory.Constructions.Pi
import Mathlib.MeasureTheory.Measure.WithDensity
import Mathlib.MeasureTheory.Measure.ProbabilityMeasure
import Mathlib.Analysis.SpecialFunctions.Exp

noncomputable section

open MeasureTheory MeasureTheory.Measure
open scoped ENNReal

namespace YangMillsInfra
namespace LatticeGauge

universe u

/-- A compact gauge group: a compact Hausdorff second-countable topological group,
equipped with its Borel σ-algebra.  This is exactly the hypothesis under which
normalized Haar measure exists. -/
class CompactGaugeGroup (G : Type u) extends Group G, TopologicalSpace G, MeasurableSpace G,
    BorelSpace G, CompactSpace G, T2Space G, IsTopologicalGroup G, SecondCountableTopology G

/-- The whole group as a positive compact set, used to normalize Haar measure so that
it is a probability measure. -/
def univPositiveCompacts (G : Type u) [TopologicalSpace G] [CompactSpace G] [Nonempty G] :
    TopologicalSpace.PositiveCompacts G where
  carrier := Set.univ
  isCompact' := isCompact_univ
  interior_nonempty' := by simp

/-- Normalized Haar measure on a compact group.  Normalization is with respect to the
whole group, so its total mass is one (`haarMeasure_self`). -/
def haarProb (G : Type u) [Group G] [TopologicalSpace G] [MeasurableSpace G] [BorelSpace G]
    [CompactSpace G] [T2Space G] [Nonempty G] [IsTopologicalGroup G] : Measure G :=
  haarMeasure (univPositiveCompacts G)

/-- The product Haar measure on link configurations. -/
def linkHaar (Λ : LatticeData.{u}) (G : Type u) [CompactGaugeGroup G] :
    Measure (Λ.L → G) :=
  Measure.pi (fun _ => haarProb G)

variable {Λ : LatticeData.{u}} {G : Type u} [CompactGaugeGroup G]

/-- The Boltzmann weight of the Wilson action, as an `ℝ≥0∞` density. -/
def wilsonWeight (χ : CharacterData G) (β : ℝ) (U : Λ.L → G) : ℝ≥0∞ :=
  ENNReal.ofReal (Real.exp (- wilsonAction χ β U))

/-- The partition function: the total mass of the unnormalized Wilson measure. -/
def partitionFunction (χ : CharacterData G) (β : ℝ) : ℝ≥0∞ :=
  ∫⁻ U, wilsonWeight (Λ := Λ) χ β U ∂(linkHaar Λ G)

/-- The normalized Wilson measure. -/
def wilsonMeasure (χ : CharacterData G) (β : ℝ) (Z : ℝ≥0∞) : Measure (Λ.L → G) :=
  (linkHaar Λ G).withDensity fun U => wilsonWeight (Λ := Λ) χ β U / Z

/-- The normalized Wilson measure has total mass one, provided `Z` is the partition
function and is finite and nonzero.  (Positivity and finiteness of the partition
function are the routine compactness/continuity facts isolated here as hypotheses.) -/
theorem wilsonMeasure_univ (χ : CharacterData G) (β : ℝ) {Z : ℝ≥0∞}
    (hZ : Z = partitionFunction (Λ := Λ) χ β) (h0 : Z ≠ 0) (htop : Z ≠ ∞) :
    wilsonMeasure (Λ := Λ) χ β Z Set.univ = 1 := by
  rw [wilsonMeasure, withDensity_apply _ MeasurableSet.univ, Measure.restrict_univ]
  have hfun : (fun U => wilsonWeight (Λ := Λ) χ β U / Z)
      = fun U => wilsonWeight (Λ := Λ) χ β U * Z⁻¹ := by
    funext U; rw [div_eq_mul_inv]
  rw [hfun, lintegral_mul_const' Z⁻¹ _ (by simpa using h0)]
  have hZeq : (∫⁻ U, wilsonWeight (Λ := Λ) χ β U ∂(linkHaar Λ G)) = Z := hZ.symm
  rw [hZeq]
  exact ENNReal.mul_inv_cancel h0 htop

/-- **Gauge invariance of the Wilson measure**, isolated as an interface.

Two inputs are needed:
* `haar_invariant`: the product Haar measure is invariant under the gauge action
  (bi-invariance of Haar measure on a compact group — gap G-3);
* `weight_invariant`: the Wilson weight is constant on orbits
  (`wilsonAction_invariant`, already proved).

Given these, invariance of `wilsonMeasure` follows by change of variables; this is the
statement that the Wilson measure really is a measure *on gauge orbits* in the sense of
`Gauge/Orbits.lean`. -/
structure WilsonGaugeInvariance (χ : CharacterData G) (β : ℝ) (Z : ℝ≥0∞) : Prop where
  /-- Product Haar measure is invariant under the gauge action. -/
  haar_invariant : ∀ g : Λ.V → G,
    Measure.map (fun U => g • U) (linkHaar Λ G) = linkHaar Λ G
  /-- The density is invariant as a function. -/
  weight_invariant : ∀ (g : Λ.V → G) (U : Λ.L → G),
    wilsonWeight (Λ := Λ) χ β (g • U) = wilsonWeight (Λ := Λ) χ β U

/-- The density half of `WilsonGaugeInvariance` is a theorem, not an assumption. -/
theorem wilsonWeight_invariant (χ : CharacterData G) (β : ℝ) (g : Λ.V → G) (U : Λ.L → G) :
    wilsonWeight (Λ := Λ) χ β (g • U) = wilsonWeight (Λ := Λ) χ β U := by
  unfold wilsonWeight
  rw [wilsonAction_invariant]

end LatticeGauge
end YangMillsInfra
