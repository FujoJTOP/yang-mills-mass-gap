/-
GL8: the handoff assembly lemma is group-blind.

The gauge group enters only through four real constants per group
`q0, q1, δW, δS` and two antitone functions `q`, `gsq`.  The conclusion
`min δW δS ≤ cert μ` never mentions the group.  Hence "group uniformity"
of the assembled mass-gap statement reduces to uniformity of these four
scalars -- it is not a new analytic statement about Lie groups.

This restates `YMStatement.handoff_gap` (lean/YMStatement/Handoff.lean)
and adds the uniform-over-a-class form and the positivity step.  It does
NOT prove `q0 ≥ q1` for any group; that is the open handoff hypothesis
(gluing R25), carried as the field `overlap`.

Run: bash check-glue-c.sh   (uses the pinned leanprover/lean4:v4.33.0-rc1)
-/

import Mathlib.Data.Real.Basic

namespace GlueC

/-- Per-group handoff data. -/
structure HandoffData where
  q : ℝ → ℝ
  gsq : ℝ → ℝ
  μup : ℝ
  q0 : ℝ
  q1 : ℝ
  δW : ℝ
  δS : ℝ
  μup_pos : 0 < μup
  q_up : q μup = q0
  overlap : q1 ≤ q0
  q_antitone : ∀ a b : ℝ, 0 < a → a ≤ b → q b ≤ q a
  gsq_antitone : ∀ a b : ℝ, 0 < a → a ≤ b → gsq b ≤ gsq a

/-- Covering lemma: weak domain `[μup,∞)`, strong domain `{q ≥ q1}`; overlap
`q1 ≤ q0` makes them cover `(0,∞)`. -/
theorem handoff_gap (D : HandoffData) (cert : ℝ → ℝ)
    (hW : ∀ μ : ℝ, 0 < μ → D.gsq μ ≤ D.gsq D.μup → D.δW ≤ cert μ)
    (hS : ∀ μ : ℝ, 0 < μ → D.q1 ≤ D.q μ → D.δS ≤ cert μ) :
    ∀ μ : ℝ, 0 < μ → min D.δW D.δS ≤ cert μ := by
  intro μ hμ
  by_cases hle : μ ≤ D.μup
  · have hq : D.q1 ≤ D.q μ := by
      calc
        D.q1 ≤ D.q0 := D.overlap
        _ = D.q D.μup := D.q_up.symm
        _ ≤ D.q μ := D.q_antitone μ D.μup hμ hle
    exact le_trans (min_le_right _ _) (hS μ hμ hq)
  · have hlt : D.μup < μ := lt_of_not_ge hle
    have hg : D.gsq μ ≤ D.gsq D.μup :=
      D.gsq_antitone D.μup μ D.μup_pos (le_of_lt hlt)
    exact le_trans (min_le_left _ _) (hW μ hμ hg)

/-- GL8 reduction.  The index type `ι` is an arbitrary class of groups and is
absent from the conclusion: the assembled bound is pointwise in `ι`. -/
theorem uniform_over_groups {ι : Type*} (D : ι → HandoffData) (cert : ι → ℝ → ℝ)
    (hW : ∀ (i : ι) (μ : ℝ), 0 < μ →
      (D i).gsq μ ≤ (D i).gsq (D i).μup → (D i).δW ≤ cert i μ)
    (hS : ∀ (i : ι) (μ : ℝ), 0 < μ →
      (D i).q1 ≤ (D i).q μ → (D i).δS ≤ cert i μ) :
    ∀ (i : ι) (μ : ℝ), 0 < μ → min (D i).δW (D i).δS ≤ cert i μ :=
  fun i => handoff_gap (D i) (cert i) (hW i) (hS i)

/-- Per-group positivity of the assembled constant. -/
theorem min_pos {a b : ℝ} (ha : 0 < a) (hb : 0 < b) : 0 < min a b :=
  lt_min ha hb

/-- A *class-uniform* lower bound is an extra hypothesis `δ`, not a consequence
of pointwise positivity.  (Concretely, `δW = δS = 1/(n+1)` is positive for
every `n` but has infimum `0`; see `gl8_check.py` (G4).) -/
theorem class_uniform {ι : Type*} (D : ι → HandoffData) (δ : ℝ)
    (hδ : 0 < δ) (h : ∀ i, δ ≤ min (D i).δW (D i).δS) :
    ∀ i, 0 < min (D i).δW (D i).δS :=
  fun i => lt_of_lt_of_le hδ (h i)

/-! ## GL8-S: constant supply (✅-implication)

Following the GL9 pattern (`✅ as an implication`), the analytic content is a
**named input**, not proved here.  `IN-3` (weak side, `CPB-W`) supplies
`(δW, λw)` with its certificate bound; `IN-4` (strong side, `CPB-S`) supplies
`(δS, q1)`.  GL8 only reads their positivity.

The `bound` fields have exactly the shape of `HandoffData.hW` / `hS`, so
supplying them instantiates `uniform_over_groups`; the remaining assembly
hypothesis `overlap : q1 ≤ q0` is `IN-5` (the handoff overlap). -/

/-- `IN-3`: weak-side certificate lower bound (named analytic input). -/
structure WeakSideInput where
  cert : ℝ → ℝ
  lam : ℝ → ℝ
  δW : ℝ
  lamw : ℝ
  δW_pos : 0 < δW
  lamw_pos : 0 < lamw
  bound : ∀ μ : ℝ, 0 < μ → lam μ ≤ lamw → δW ≤ cert μ

/-- `IN-4`: strong-side certificate lower bound (named analytic input). -/
structure StrongSideInput where
  cert : ℝ → ℝ
  q : ℝ → ℝ
  δS : ℝ
  q1 : ℝ
  δS_pos : 0 < δS
  q1_pos : 0 < q1
  bound : ∀ μ : ℝ, 0 < μ → q1 ≤ q μ → δS ≤ cert μ

/-- **GL8-S (✅-implication).**  If `IN-3` and `IN-4` supply `(δW, δS)` for
every group, then `min δW δS > 0` for every group.  The antecedent is the
named input; this theorem proves only the implication. -/
theorem supply_min_pos {ι : Type*}
    (W : ι → WeakSideInput) (S : ι → StrongSideInput) :
    ∀ i, 0 < min (W i).δW (S i).δS :=
  fun i => min_pos (W i).δW_pos (S i).δS_pos

end GlueC

/-- info: 'GlueC.uniform_over_groups' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueC.uniform_over_groups

/-- info: 'GlueC.min_pos' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueC.min_pos

/-- info: 'GlueC.class_uniform' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueC.class_uniform

/-- info: 'GlueC.supply_min_pos' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueC.supply_min_pos
