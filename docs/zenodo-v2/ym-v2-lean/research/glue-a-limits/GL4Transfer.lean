/-
GL4Transfer.lean — machine-checked core of glue lemma GL4 (weak-limit transfer).

Two statements, matching REPORT.md §3:

* `rp_transfer` — reflection positivity is preserved by a weak limit.
  The weak-convergence hypothesis is Mathlib's `Tendsto μs F (𝓝 μ)` in the
  weak topology on `ProbabilityMeasure Ω`; the continuity hypothesis is
  literally `Continuous θ`, and the test function is the bounded continuous
  `ω ↦ A (θ ω) * A ω` built from `A : Ω →ᵇ ℝ`.

  Scope: this is the bounded-continuous form.  The extension to unbounded
  `A ∈ L²(μ)` is the L²-truncation closure proved on paper in REPORT.md §3.2
  (truncate `A_R = A · 1_{|A| ≤ R}`, use `A_R → A` in `L²(μ)` and Cauchy–
  Schwarz).  It is *not* formalised here.

* `clustering_transfer` — an exponential clustering bound passes to the pointwise
  limit.  The uniform-in-`a` bound `(U_a)` is an **explicit hypothesis**
  (`hU`), exactly as demanded: without it the statement is false (massive
  Gaussian with `m → 0`; see `gl4_gaussian_counterexample.py`).

* `exists_convergent_subseq_of_abs_le_one` — the one-dimensional Bolzano–
  Weierstrass step used by GL3-B (uniformly bounded moments have convergent
  subsequences).  The finite/countable diagonal over several moments is prose
  in REPORT.md §2.2.

No `sorry`; `#print axioms` below reports only the three standard axioms.
-/
import Mathlib.MeasureTheory.Measure.Portmanteau
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Topology.MetricSpace.Sequences
import Mathlib.Topology.Order.OrderClosed

open MeasureTheory Filter Topology
open scoped BoundedContinuousFunction

namespace GlueA

variable {Ω : Type*} [TopologicalSpace Ω] [MeasurableSpace Ω] [OpensMeasurableSpace Ω]

/-- **GL4-RP, bounded-continuous form.**
If `μs → μ` weakly (in the weak topology of probability measures), each `μs`
satisfies reflection positivity for the involution `θ`, and `A` is bounded
continuous, then `μ` satisfies the same reflection-positivity inequality.

The only continuity input is `Continuous θ`; `OpensMeasurableSpace Ω` makes
bounded continuous functions measurable. -/
theorem rp_transfer
    (θ : Ω → Ω) (hθ : Continuous θ)
    {ι : Type*} {F : Filter ι} [NeBot F] (μs : ι → ProbabilityMeasure Ω)
    (μ : ProbabilityMeasure Ω)
    (hweak : Tendsto μs F (𝓝 μ))
    (A : Ω →ᵇ ℝ)
    (hRP : ∀ i, 0 ≤ ∫ ω, A (θ ω) * A ω ∂(μs i : Measure Ω)) :
    0 ≤ ∫ ω, A (θ ω) * A ω ∂(μ : Measure Ω) := by
  let g : Ω →ᵇ ℝ := A.compContinuous ⟨θ, hθ⟩
  let f : Ω →ᵇ ℝ := g * A
  have hf : (fun ω => f ω) = fun ω => A (θ ω) * A ω := by
    funext ω
    simp only [f, g, BoundedContinuousFunction.mul_apply,
      BoundedContinuousFunction.compContinuous_apply, ContinuousMap.coe_mk]
  have hlim := (ProbabilityMeasure.tendsto_iff_forall_integral_tendsto.mp hweak) f
  have hRP' : ∀ᶠ i in F, 0 ≤ ∫ ω, f ω ∂(μs i : Measure Ω) :=
    Eventually.of_forall (fun i => by simpa only [← hf] using hRP i)
  have key : 0 ≤ ∫ ω, f ω ∂(μ : Measure Ω) := ge_of_tendsto hlim hRP'
  rwa [hf] at key

/-- **GL4-C, machine-checked.**
A uniform (in the index `i`, i.e. in `a`) exponential bound on the connected
correlations passes to the pointwise limit.  `(U_a)` is the hypothesis `hU`.

Note the bound function `C * exp (−δ * x)` is not assumed nonnegative; the
argument is the elementary "pointwise inequality survives a limit". -/
theorem clustering_transfer
    {ι : Type*} {F : Filter ι} [NeBot F] (c : ι → ℝ → ℝ) (clim : ℝ → ℝ)
    (hconv : ∀ x, Tendsto (fun i => c i x) F (𝓝 (clim x)))
    (C δ : ℝ)
    (hU : ∀ i x, |c i x| ≤ C * Real.exp (-δ * x)) :
    ∀ x, |clim x| ≤ C * Real.exp (-δ * x) := by
  intro x
  have habs : Tendsto (fun i => |c i x|) F (𝓝 |clim x|) :=
    (continuous_abs.tendsto _).comp (hconv x)
  exact le_of_tendsto habs (Eventually.of_forall (fun i => hU i x))

/-- One-dimensional Bolzano–Weierstrass: a sequence in `[-1,1]` has a convergent
subsequence.  This is the `n = 1` compactness step of GL3-B. -/
theorem exists_convergent_subseq_of_abs_le_one
    (S : ℕ → ℝ) (hS : ∀ n, |S n| ≤ 1) :
    ∃ (L : ℝ) (φ : ℕ → ℕ), StrictMono φ ∧ Tendsto (fun k => S (φ k)) atTop (𝓝 L) := by
  obtain ⟨a, _ha, φ, hφ, hlim⟩ := tendsto_subseq_of_bounded (s := Set.Icc (-1 : ℝ) 1)
    (Metric.isBounded_Icc (-1 : ℝ) 1) (x := S) (fun n => abs_le.mp (hS n))
  exact ⟨a, φ, hφ, hlim⟩

end GlueA

-- Axiom audit: standard three only, no `sorryAx`.
#print axioms GlueA.rp_transfer
#print axioms GlueA.clustering_transfer
#print axioms GlueA.exists_convergent_subseq_of_abs_le_one
