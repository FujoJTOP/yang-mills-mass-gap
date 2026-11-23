/-
GL3Limits.lean — machine-checked skeleton of glue lemma GL3 (uniform-in-`a`
compactness), with Minlos as an explicit external hypothesis.

Modelling.  Fix L and enumerate the countably many moments
`(n, f₁, …, f_n)` by `j : ℕ`; `k : ℕ` indexes the lattice spacing `a_k → 0`.
Thus a moment family is `S : ℕ → (ℕ → ℝ)`, with `S k j` the `j`-th moment at
spacing `a_k`, and `|S k j| ≤ C j` is the uniform bound (UVS in REPORT.md §2).

What is proved here (steps 1–2):

* `banach_alaoglu_compact` — Mathlib's Banach–Alaoglu
  (`WeakDual.isCompact_closedBall`): weak-* closed balls are compact.
* `exists_subseq_tendsto_pi` — uniform coordinate bounds ⇒ the moment family
  lies in a compact product of intervals (`isCompact_univ_pi`, Tychonoff);
  compact + first-countable ⇒ a subsequence converges simultaneously in every
  coordinate (`IsCompact.tendsto_subseq` + `tendsto_pi_nhds`).

What is **not** proved here (step 3, an explicit hypothesis):

* `hMinlos` — Minlos' theorem (equivalently, the kernel form of Prokhorov):
  a positive, normalized, uniformly bounded moment functional on the cylinder
  algebra is the moment functional of a probability measure.  Mathlib has no
  such theorem, so it appears verbatim as a hypothesis in the signature.
  The predicate `P` abstracts the positivity/normalization of the limit; in the
  application it is supplied by GL4 (reflection positivity at the limit).

No `sorry`; the axioms of the two theorems are frozen by `#guard_msgs` below.
-/
import Mathlib.Analysis.Normed.Module.WeakDual
import Mathlib.MeasureTheory.Measure.ProbabilityMeasure

open MeasureTheory Filter Topology
open scoped Topology

namespace GlueA

/-- **Banach–Alaoglu (Mathlib, `WeakDual.isCompact_closedBall`).**
Closed balls of the dual of a normed space are compact in the weak-* topology.
This is the abstract form of the compactness step; the concrete moment encoding
below uses Tychonoff (`isCompact_univ_pi`), which is the same statement in
coordinates. -/
theorem banach_alaoglu_compact
    {𝕜 E : Type*} [NontriviallyNormedField 𝕜] [ProperSpace 𝕜]
    [SeminormedAddCommGroup E] [NormedSpace 𝕜 E]
    (x' : StrongDual 𝕜 E) (r : ℝ) :
    IsCompact (WeakDual.toStrongDual ⁻¹' Metric.closedBall x' r : Set (WeakDual 𝕜 E)) :=
  WeakDual.isCompact_closedBall x' r

/-- **Steps 1–2: compactness and the countable diagonal.**
A moment family with uniform coordinate bounds has a subsequence converging
simultaneously in every coordinate.

* The compact set is `∏_j [-C_j, C_j]` (`isCompact_univ_pi`, Tychonoff).
* `ℕ → ℝ` is first-countable in the product topology, so
  `IsCompact.tendsto_subseq` extracts a converging subsequence.
* `tendsto_pi_nhds` turns product convergence into coordinatewise convergence. -/
theorem exists_subseq_tendsto_pi (S : ℕ → (ℕ → ℝ)) (C : ℕ → ℝ)
    (hS : ∀ k j, |S k j| ≤ C j) :
    ∃ (L : ℕ → ℝ) (φ : ℕ → ℕ), StrictMono φ ∧
      ∀ j, Tendsto (fun k => S (φ k) j) atTop (𝓝 (L j)) := by
  let K : Set (ℕ → ℝ) := Set.univ.pi (fun j => Set.Icc (-(C j)) (C j))
  have hK : IsCompact K := isCompact_univ_pi (fun j => isCompact_Icc)
  have hmem : ∀ k, S k ∈ K := by
    intro k
    simp only [K, Set.mem_pi, Set.mem_univ, true_implies]
    exact fun j => abs_le.mp (hS k j)
  obtain ⟨L, _hL, φ, hφ, hlim⟩ := hK.tendsto_subseq hmem
  exact ⟨L, φ, hφ, fun j => tendsto_pi_nhds.mp hlim j⟩

/-- A bounded moment sequence has a probability-measure realization.  This is
the conclusion of Minlos, not a definition we prove. -/
def HasRealization (L : ℕ → ℝ) : Prop :=
  ∃ μ : Measure (ℕ → ℝ), IsProbabilityMeasure μ ∧ ∀ j, ∫ x, x j ∂μ = L j

/-- **GL3, machine-checked skeleton, Minlos as an explicit hypothesis.**

`P` abstracts positivity/normalization of the limit moment functional (in the
application: reflection positivity at the limit, GL4-RP).  `hMinlos` is the
external theorem (Minlos / Prokhorov kernel form), **not** proved in Mathlib
and **not** proved here; it is written verbatim as a hypothesis. -/
theorem gl3_imp_with_minlos
    (S : ℕ → (ℕ → ℝ)) (C : ℕ → ℝ) (hS : ∀ k j, |S k j| ≤ C j)
    (P : (ℕ → ℝ) → Prop)
    (hMinlos : ∀ L, (∀ j, |L j| ≤ C j) → P L → HasRealization L) :
    ∃ (L : ℕ → ℝ) (φ : ℕ → ℕ), StrictMono φ ∧
      (∀ j, Tendsto (fun k => S (φ k) j) atTop (𝓝 (L j))) ∧
      (∀ _hP : P L, HasRealization L) := by
  obtain ⟨L, φ, hφ, hlim⟩ := exists_subseq_tendsto_pi S C hS
  refine ⟨L, φ, hφ, hlim, fun hP => hMinlos L ?_ hP⟩
  intro j
  have habs : Tendsto (fun k => |S (φ k) j|) atTop (𝓝 |L j|) :=
    (continuous_abs.tendsto _).comp (hlim j)
  exact le_of_tendsto habs (Eventually.of_forall (fun k => hS (φ k) j))

end GlueA

/-- info: 'GlueA.exists_subseq_tendsto_pi' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueA.exists_subseq_tendsto_pi

/-- info: 'GlueA.gl3_imp_with_minlos' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueA.gl3_imp_with_minlos
