/-
  GL2a — reflection positivity descends to the blocked measure.

  Statement (REPORT.md §2.1): let `μ` be a probability measure with the
  Osterwalder–Seiler reflection positivity (RP) property w.r.t. a configuration
  reflection `θf`, let `F : Ω → Ωc` be a deterministic block map satisfying

    (B1) reflection covariance   `θc (F x) = F (θf x)`  for all `x`,
    (B2) half-space separation   `A ∘ F` is RP-admissible whenever `A` is,

  and let `μc := μ.map F` be the pushforward.  Then for every `A` that is
  measurable on the *coarse* configuration space,

    `0 ≤ ∫ y, A (θc y) * A y ∂μc`.

  The proof is the five-liner of REPORT.md §2.2: block-measurable functions have
  the same integral under `μ` and under `μc` (`integral_map`), and (B1) turns the
  coarse RP form into the fine one.

  Design notes (for `lean-infra`):

  * `MeasureTheory.condExp` is *not* needed.  The conditional expectation onto
    `σ(F)` is by definition the orthogonal projection whose range is the
    block-measurable functions, and the only property used below is that such a
    function has equal integrals under `μ` and `μc` — which is exactly
    `MeasureTheory.integral_map`.  Keeping `condExp` out avoids a heavyweight
    import for zero gain.
  * `θE = Eθ` (θ-covariance of the conditional expectation) is *not* an
    assumption: it is a consequence of (B1).  We assume the checkable primitive
    (B1) instead.  This is the difference flagged in REPORT.md §2.2: the
    "θ-covariant conditional expectation" framing is automatic; the geometry is
    not.
  * `P` is an abstract predicate standing for "supported in the `t > 0`
    half-space".  Take `P := C_RP` at instantiation.  (B2) is then the single
    geometric input.

  Real-valued form; the complex form replaces `*` by `star` and is otherwise
  identical.  No `sorry`, standard axioms only (see `check.sh`).
-/
import Mathlib

open MeasureTheory

namespace YMGlue

variable {Ω Ωc : Type*} [MeasurableSpace Ω] [MeasurableSpace Ωc]

/-- **GL2a.** Reflection positivity descends to the pushforward (blocked) measure.

`F` is the deterministic block map, `θf`/`θc` the fine/coarse reflections,
`P` the fine RP cone (`t > 0` half-space), and `hRP` the fine reflection
positivity hypothesis. -/
theorem rp_descends
    (μ : Measure Ω) [IsProbabilityMeasure μ]
    (F : Ω → Ωc) (hF : Measurable F)
    (θf : Ω → Ω) (θc : Ωc → Ωc)
    -- (B1) reflection covariance of the block map
    (hB1 : ∀ x, θc (F x) = F (θf x))
    -- fine reflection positivity
    (P : (Ω → ℝ) → Prop)
    (hRP : ∀ φ : Ω → ℝ, P φ → 0 ≤ ∫ x, φ (θf x) * φ x ∂μ)
    -- the coarse observable, and (B2): its pullback is fine-RP-admissible
    (A : Ωc → ℝ) (hA : Measurable A) (hθc : Measurable θc)
    (hB2 : P (fun x => A (F x))) :
    0 ≤ ∫ y, A (θc y) * A y ∂(μ.map F) := by
  -- Push the integral to `μ` along `F`.
  have hmeas : AEStronglyMeasurable (fun y => A (θc y) * A y) (μ.map F) :=
    ((hA.comp hθc).mul hA).aestronglyMeasurable
  have hmap : ∫ y, A (θc y) * A y ∂(μ.map F)
      = ∫ x, A (θc (F x)) * A (F x) ∂μ := by
    rw [integral_map hF.aemeasurable hmeas]
  rw [hmap]
  -- (B1) identifies the coarse-composed integrand with the fine RP form of `A ∘ F`.
  have hshape : (fun x => A (θc (F x)) * A (F x))
      = fun x => (fun y => A (F y)) (θf x) * (fun y => A (F y)) x := by
    funext x
    rw [hB1 x]
  rw [hshape]
  exact hRP (fun y => A (F y)) hB2

/-- **GL2a, second half: the compression is positive.**

If `T` has nonnegative form then so does `P T P`, for *any* bounded `P`.
This is the "RP is preserved at every scale ⇒ `T_c = E T_f E ≥ 0`" step, and the
point of stating it is that it is one substitution: no projection property of
`E` is used.  The content of GL2a is therefore entirely in (B1)+(B2), not here. -/
theorem compression_form_nonneg
    {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (T P : E →L[ℝ] E) (hT : ∀ φ, 0 ≤ inner ℝ φ (T φ)) (ψ : E) :
    0 ≤ inner ℝ (P ψ) (T (P ψ)) :=
  hT (P ψ)

end YMGlue

/-- info: 'YMGlue.rp_descends' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMGlue.rp_descends

/-- info: 'YMGlue.compression_form_nonneg' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMGlue.compression_form_nonneg
