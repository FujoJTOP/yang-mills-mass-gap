/-
# YangMillsInfra — shared vocabulary

This file fixes conventions shared by every module of the infrastructure package:

* `SpaceTime d` is `ℝ^d` with the Euclidean inner product (same abbreviation as
  `OSforGFF.Spacetime.Basic`, so the adapter layer can identify the two).
* `TestFunction d` is real Schwartz space; field configurations live in its weak dual.
* `EuclideanTheory` is the *interface* a Euclidean field theory must expose to the
  mass-gap statement layer: a generating functional `Z` on test functions.  Everything
  else (measure, OS axioms, moments) is a property of `Z`, not part of the type.

Deliberately, this file does **not** import OSforGFF: the infrastructure is meant to be
usable by any of the parallel development lines, and the adapter
`YangMillsInfra/Bridge/OSforGFF.lean` is the single place where the two vocabularies meet.
-/
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.Distribution.SchwartzSpace.Basic
import Mathlib.Data.Complex.Basic

noncomputable section

open MeasureTheory

namespace YangMillsInfra

/-- Euclidean `d`-space, `ℝ^d` with its Euclidean inner product. -/
abbrev SpaceTime (d : ℕ) : Type := EuclideanSpace ℝ (Fin d)

/-- Real-valued Schwartz test functions on Euclidean space. -/
abbrev TestFunction (d : ℕ) : Type := SchwartzMap (SpaceTime d) ℝ

/-- Complex-valued Schwartz test functions. -/
abbrev TestFunctionℂ (d : ℕ) : Type := SchwartzMap (SpaceTime d) ℂ

/-- Field configurations: tempered distributions, i.e. the weak dual of Schwartz space.

    We only need the *type* here; the measure and its σ-algebra come from
    `OSforGFF`/`bochner` in the adapter layer. -/
abbrev Configuration (d : ℕ) : Type := WeakDual ℝ (TestFunction d)

/-- A Euclidean field theory, as seen by the mass-gap statement layer.

    `Z` is the Glimm–Jaffe generating functional `Z[J] = ∫ exp(i⟨ω,J⟩) dμ(ω)`.
    Bundling only `Z` keeps the statement layer independent of the choice of
    configuration space and of how the measure is constructed; the map from a
    measure to `Z` is `EuclideanTheory.ofMeasure` in the adapter.

    Note that the axioms OS0–OS4 are predicates on `Z` (or on the underlying measure);
    they are intentionally *not* fields, so that the same type can be used for
    theories for which only some axioms are known. -/
structure EuclideanTheory (d : ℕ) where
  /-- Generating functional `Z[J]`. -/
  Z : TestFunction d → ℂ

variable {d : ℕ}

/-- A *translation action* on test functions, abstracted as a monoid action.

    `OSforGFF` instantiates this with `SchwartzMap.translate` (which mathlib does not
    have).  Keeping it abstract lets the statement layer state clustering without
    importing the distribution library. -/
structure TranslationAction (d : ℕ) where
  /-- Translate by a spacetime vector. -/
  translate : SpaceTime d → TestFunction d → TestFunction d
  translate_zero : ∀ f, translate 0 f = f
  translate_add : ∀ a b f, translate (a + b) f = translate a (translate b f)

/-- Two-point clustering of the generating functional, *with a rate* `m`.

    This is the Euclidean form the mass gap is equivalent to (see
    `Statement/MassGap.lean` and `RG/Handoff.lean`); the bare clustering axiom OS4 is
    this statement with `Real.exp (-m * ‖a‖)` replaced by an arbitrary null sequence,
    and does **not** imply a spectral gap (hexagonmath 2611.00492). -/
def ExponentialClustering (T : EuclideanTheory d) (τ : TranslationAction d) (m : ℝ) : Prop :=
  0 < m ∧
    ∀ f g : TestFunction d, ∃ C : ℝ, 0 ≤ C ∧
      ∀ a : SpaceTime d,
        ‖T.Z (f + τ.translate a g) - T.Z f * T.Z g‖ ≤ C * Real.exp (-m * ‖a‖)

/-- Bare clustering (OS4 as written in `OSforGFF.OS.Axioms.OS4_Clustering`), for comparison:
the same difference tends to zero along `‖a‖ → ∞`, with no rate. -/
def Clustering (T : EuclideanTheory d) (τ : TranslationAction d) : Prop :=
  ∀ f g : TestFunction d, ∀ ε : ℝ, 0 < ε → ∃ R : ℝ, 0 < R ∧
    ∀ a : SpaceTime d, R < ‖a‖ →
      ‖T.Z (f + τ.translate a g) - T.Z f * T.Z g‖ < ε

end YangMillsInfra
