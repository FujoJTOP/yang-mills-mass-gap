/-
# Measures on gauge orbits

Two facts, of very different difficulty, are conflated under "the measure on the space
of gauge orbits":

1. **Pushing forward.**  A gauge-invariant weight defines a measure on the orbit space
   by pushforward along the quotient map, and gauge-invariant observables have the
   same expectation computed upstairs or downstairs.  This is elementary
   (`integral_map`) and is what every lattice computation actually uses.

2. **Disintegrating.**  Writing the upstairs measure as an integral of measures on
   orbits requires a measurable selection of orbit measures; equivalently a
   fundamental domain and a quotient-measure theorem.  mathlib has this only for
   *countable* group actions (`MeasureTheory.QuotientMeasureEqMeasurePreimage`,
   `IsFundamentalDomain`), while the lattice gauge group `Λ.V → G` is a compact Lie
   group and uncountable.  This is gap **G-2** in REPORT.md; here it appears as the
   `GaugeFixing` structure, whose `integral_eq` field is exactly the missing theorem.

The separation matters: nothing in the mass-gap program needs (2) — it needs (1) plus
the fact that the Wilson weight is constant on orbits, which is
`LatticeGauge.wilsonAction_invariant`.
-/
import YangMillsInfra.Gauge.Lattice
import Mathlib.MeasureTheory.Measure.Map
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.MeasureTheory.Group.FundamentalDomain
import Mathlib.GroupTheory.GroupAction.Quotient

noncomputable section

open MeasureTheory

namespace YangMillsInfra
namespace LatticeGauge

variable {Λ : LatticeData.{u}} {G : Type u} [Group G] [MeasurableSpace G]

/-- The orbit equivalence relation of the gauge action on link configurations. -/
abbrev gaugeOrbitRel : Setoid (Λ.L → G) := MulAction.orbitRel (Λ.V → G) (Λ.L → G)

/-- The space of gauge orbits. -/
abbrev OrbitSpace (Λ : LatticeData.{u}) (G : Type u) [Group G] : Type u :=
  Quotient (gaugeOrbitRel (Λ := Λ) (G := G))

/-- A gauge-invariant observable descends to the orbit space. -/
def orbitLift {E : Type*} (O : (Λ.L → G) → E)
    (hO : ∀ (g : Λ.V → G) (U : Λ.L → G), O (g • U) = O U) :
    OrbitSpace Λ G → E :=
  Quotient.lift O fun U V h => by
    obtain ⟨g, hg⟩ := (MulAction.orbitRel_apply).mp h
    rw [← hg]
    exact hO g V

omit [MeasurableSpace G] in
@[simp] theorem orbitLift_mk {E : Type*} (O : (Λ.L → G) → E)
    (hO : ∀ (g : Λ.V → G) (U : Λ.L → G), O (g • U) = O U) (U : Λ.L → G) :
    orbitLift O hO (Quotient.mk'' U) = O U := rfl

/-- The pushforward of a configuration measure to the orbit space. -/
def orbitPushforward (μ : Measure (Λ.L → G)) : Measure (OrbitSpace Λ G) :=
  μ.map (Quotient.mk'')

/-- **Expectations of gauge-invariant observables do not depend on whether one works
on configurations or on orbits.**  Stated for a general measure; the gauge invariance
of the weight enters only through `hO`. -/
theorem integral_orbitPushforward
    (μ : Measure (Λ.L → G)) {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (O : (Λ.L → G) → E) (hO : ∀ (g : Λ.V → G) (U : Λ.L → G), O (g • U) = O U)
    (hπ : AEMeasurable (Quotient.mk'' : (Λ.L → G) → OrbitSpace Λ G) μ)
    (hOm : AEStronglyMeasurable (orbitLift O hO) (orbitPushforward μ)) :
    ∫ ω : OrbitSpace Λ G, orbitLift O hO ω ∂(orbitPushforward μ) = ∫ U, O U ∂μ := by
  unfold orbitPushforward at hOm ⊢
  simpa using integral_map hπ hOm

/-- **Gauge fixing, as an interface.**

`slice` is a set of configurations meeting every orbit exactly once (a.e. — mathlib's
`IsFundamentalDomain`), and `integral_eq` says that integrating a gauge-invariant
observable over configurations equals integrating its descent over the orbit space.

For a countable gauge group, `IsFundamentalDomain` together with
`MeasureTheory.QuotientMeasureEqMeasurePreimage` supplies `integral_eq` in mathlib.
For the compact Lie groups used in lattice gauge theory (uncountable),
`integral_eq` is *not* available from mathlib: this is gap G-2.  Recording it as a
field keeps the gap visible and local. -/
structure GaugeFixing (μ : Measure (Λ.L → G)) where
  /-- A gauge slice. -/
  slice : Set (Λ.L → G)
  /-- The slice meets every orbit in exactly one point, almost everywhere. -/
  isFundamentalDomain : IsFundamentalDomain (Λ.V → G) slice μ
  /-- The missing theorem, for uncountable gauge groups: integration over
  configurations agrees with integration of the descended observable over the orbit
  space. -/
  integral_eq : ∀ {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (O : (Λ.L → G) → E) (hO : ∀ (g : Λ.V → G) (U : Λ.L → G), O (g • U) = O U),
    ∫ U, O U ∂μ = ∫ ω : OrbitSpace Λ G, orbitLift O hO ω ∂(orbitPushforward μ)

/-- **Tree gauge fixing is the intended instance.**  On a lattice with a maximal tree,
every configuration is gauge-equivalent to a unique one that is the identity on the
tree.  Formalizing this requires the tree combinatorics and an induction along the
tree; it is the concrete instance of `GaugeFixing` that a lattice-YM development
should build first.

This is recorded as a definition of the slice, not as a theorem, so that the statement
of the intended instance is visible without pretending it is proved here. -/
def treeSlice (T : Finset Λ.L) : Set (Λ.L → G) :=
  {U | ∀ l ∈ T, U l = 1}

end LatticeGauge
end YangMillsInfra
