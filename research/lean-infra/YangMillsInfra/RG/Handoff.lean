/-
# The weak → strong coupling handoff, as a machine-checkable inequality

This file formalizes §3 of `research/physics_handoff.md`: the statement that the
control parameter across the two regimes is the RG-invariant ratio

    q(μ) = m(μ) / μ  =  a / ξ_phys(μ),

not `g²(μ)`.  The mathematical content is small and entirely rigorous:

* weak side:  `g²(μ) ≤ g*²  →  gap(μ) ≥ δ_W`   (valid for large μ);
* strong side:`q(μ) ≥ q₁   →  gap(μ) ≥ δ_S`   (valid all the way down to μ → 0);
* overlap:    `q(μ_up) ≥ q₁`, where `μ_up` is defined by `g²(μ_up) = g*²`.

If both `g²` and `q` are strictly decreasing in `μ`, then the overlap hypothesis
`q(μ_up) ≥ q₁` is *equivalent* to `μ_up ≤ μ*` where `q(μ*) = q₁`, and implies that
every scale `μ > 0` lies in one of the two regimes.  Hence `gap ≥ min δ_W δ_S`
everywhere — which is the "step theorem is uniform in the scale" conclusion that the
RG program needs.

The two substantive inputs (weak-side and strong-side bounds) are *hypotheses* here,
because they are exactly the open mathematical content.  This file proves that the
handoff itself is not where the difficulty lies.

There is no `sorry` in this file and no appeal to unproved analysis: the only
order-theoretic input is `StrictAntiOn`.
-/
import Mathlib.Order.Monotone.Basic
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Data.Real.Basic

noncomputable section

namespace YangMillsInfra
namespace RG

/-- Data for the handoff: the running coupling, the invariant ratio, the effective gap,
and the scale `μ_up` at which the weak-coupling analysis stops being valid. -/
structure HandoffData where
  /-- Effective mass over scale: `q(μ) = m(μ)/μ`. -/
  q : ℝ → ℝ
  /-- Running coupling squared. -/
  g2 : ℝ → ℝ
  /-- Effective transfer-matrix gap at scale `μ`. -/
  gap : ℝ → ℝ
  /-- Boundary of the weak-coupling region: `g2 μ_up = gstar^2`. -/
  μup : ℝ
  gstar : ℝ
  q1 : ℝ
  δW : ℝ
  δS : ℝ
  /-- `q` is strictly decreasing on the positive half-line. -/
  q_strictAnti : StrictAntiOn q (Set.Ioi 0)
  /-- `g²` is strictly decreasing on the positive half-line (asymptotic freedom). -/
  g2_strictAnti : StrictAntiOn g2 (Set.Ioi 0)
  /-- `μ_up` is a positive scale. -/
  μup_pos : 0 < μup
  /-- `μ_up` really is the weak-coupling boundary. -/
  g2_μup : g2 μup = gstar ^ 2
  /-- `g*` is a positive coupling. -/
  gstar_pos : 0 < gstar

/-- The weak-side regime: below the coupling threshold, the gap is at least `δ_W`. -/
def WeakSide (D : HandoffData) : Prop :=
  ∀ μ : ℝ, 0 < μ → D.g2 μ ≤ D.gstar ^ 2 → D.δW ≤ D.gap μ

/-- The strong-side regime: once `q ≥ q₁`, the gap is at least `δ_S`.  The statement is
deliberately not restricted to small `μ`: this is what "the strong-coupling analysis
runs all the way to the infrared" means. -/
def StrongSide (D : HandoffData) : Prop :=
  ∀ μ : ℝ, 0 < μ → D.q1 ≤ D.q μ → D.δS ≤ D.gap μ

/-- The overlap hypothesis: at the weak-coupling boundary the invariant ratio is still
in the strong-coupling regime.  This is the one hypothesis that is neither weak-side
analysis nor strong-side analysis; it is the open content. -/
def Overlap (D : HandoffData) : Prop :=
  D.q1 ≤ D.q D.μup

/-- The set of scales covered by the weak analysis. -/
def weakRegion (D : HandoffData) : Set ℝ :=
  {μ | 0 < μ ∧ D.g2 μ ≤ D.gstar ^ 2}

/-- The set of scales covered by the strong analysis. -/
def strongRegion (D : HandoffData) : Set ℝ :=
  {μ | 0 < μ ∧ D.q1 ≤ D.q μ}

/-- **The handoff lemma.**  Under weak-side and strong-side bounds and the overlap
hypothesis, the effective gap is bounded below uniformly in the scale. -/
theorem gap_uniform (D : HandoffData) (hW : WeakSide D) (hS : StrongSide D)
    (hOv : Overlap D) :
    ∀ μ : ℝ, 0 < μ → min D.δW D.δS ≤ D.gap μ := by
  intro μ hμ
  by_cases h : D.g2 μ ≤ D.gstar ^ 2
  · exact le_trans (min_le_left _ _) (hW μ hμ h)
  · push Not at h
    -- `g²` strictly decreasing gives `μ < μ_up`
    have hlt : μ < D.μup := by
      by_contra hge
      push Not at hge
      rcases lt_or_eq_of_le hge with hlt' | heq
      · -- μ_up < μ contradicts `g2 μ > g2 μ_up`
        have hmono : D.g2 μ ≤ D.g2 D.μup :=
          (D.g2_strictAnti.antitoneOn D.μup_pos hμ hlt'.le)
        rw [D.g2_μup] at hmono
        linarith
      · subst heq
        rw [D.g2_μup] at h
        linarith
    -- `q` strictly decreasing gives `q μ > q μ_up ≥ q₁`
    have hq : D.q1 ≤ D.q μ := le_trans hOv (le_of_lt (D.q_strictAnti hμ D.μup_pos hlt))
    exact le_trans (min_le_right _ _) (hS μ hμ hq)

/-- **Covering.**  The same case split shows that the weak and strong regions cover all
positive scales.  (This is the content that the overlap hypothesis buys: without it,
there is a gap of scales between `μ*` and `μ_up`.) -/
theorem regions_cover (D : HandoffData) (hOv : Overlap D) :
    Set.Ioi (0 : ℝ) ⊆ weakRegion D ∪ strongRegion D := by
  intro μ hμ
  by_cases h : D.g2 μ ≤ D.gstar ^ 2
  · exact Or.inl ⟨hμ, h⟩
  · push Not at h
    have hlt : μ < D.μup := by
      by_contra hge
      push Not at hge
      rcases lt_or_eq_of_le hge with hlt' | heq
      · have hmono : D.g2 μ ≤ D.g2 D.μup :=
          (D.g2_strictAnti.antitoneOn D.μup_pos hμ hlt'.le)
        rw [D.g2_μup] at hmono
        linarith
      · subst heq
        rw [D.g2_μup] at h
        linarith
    exact Or.inr ⟨hμ, le_trans hOv (le_of_lt (D.q_strictAnti hμ D.μup_pos hlt))⟩

/-- The overlap hypothesis is exactly non-emptiness of the overlap interval
`[μ_up, μ*]` where `q(μ*) = q₁`: if `μ*` is defined as a preimage of `q₁`, then
`q(μ_up) ≥ q₁ ↔ μ_up ≤ μ*`.  This records that `q` -- not `g²` -- is the parameter in
which the criterion is stated. -/
theorem overlap_iff_order (D : HandoffData) {μstar : ℝ} (hμstar : 0 < μstar)
    (hqstar : D.q μstar = D.q1) :
    Overlap D ↔ D.μup ≤ μstar := by
  constructor
  · intro hOv
    by_contra hlt
    push Not at hlt
    -- μstar < μup, and q is strictly decreasing, so q(μup) < q(μstar) = q₁
    have hlt' : D.q D.μup < D.q μstar := D.q_strictAnti hμstar D.μup_pos hlt
    rw [hqstar] at hlt'
    exact absurd hOv (not_le.mpr hlt')
  · intro hle
    rcases lt_or_eq_of_le hle with hlt | heq
    · have : D.q μstar < D.q D.μup := D.q_strictAnti D.μup_pos hμstar hlt
      rw [hqstar] at this
      exact this.le
    · subst heq
      unfold Overlap
      exact hqstar.ge

end RG
end YangMillsInfra
