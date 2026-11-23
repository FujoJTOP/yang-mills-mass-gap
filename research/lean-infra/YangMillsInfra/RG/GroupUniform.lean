/-
# Group-uniformity: the handoff lemma with the gauge group as a parameter

The group-uniformity analysis (`research/group-uniformity/REPORT.md` §6, §9.1–9.2)
asks for the weak→strong handoff to be stated **with the gauge group as a parameter**,
so that "uniformity in `G`" becomes the statement that four constants have uniform
bounds and nothing else depends on `G`.

This file does exactly that:

* `GroupHandoffData Grp` carries a group label and the four group-dependent constants
  `q0`, `q1`, `δW`, `δS` (plus the monotone data);
* `GroupHandoffData.gap_uniform` is the interval-cover assembly, proved;
* `uniform_over_groups` is the **no-hidden-dependence certificate**: if the constants
  are uniform across a family of groups, the gap bound is uniform, and the proof never
  inspects the group.

The four constants are the ones identified in the group-uniformity report:
`q0(G) = q(μ_up)` is the invariant ratio at the weak-coupling boundary,
`q1(G)` the strong-coupling threshold, and `δ_W(G)`, `δ_S(G)` the two gap bounds.
The conclusion `min δ_W δ_S ≤ gap μ` depends only on `q1`, `δW`, `δS` and the
monotonicity of `q`, `g²`; `q0` enters only through the overlap hypothesis.
-/
import YangMillsInfra.RG.Handoff
import Mathlib.Order.Monotone.Basic

noncomputable section

namespace YangMillsInfra
namespace RG

universe u

/-- Handoff data with the gauge group as an explicit parameter.

`Grp` is a label: no field may depend on it except through the four constants
`q0`, `q1`, `δW`, `δS` and the monotone functions `q`, `g2`, `gap`.  This is a
typing discipline, not a mathematical restriction — it is what makes the
no-hidden-dependence theorem below meaningful. -/
structure GroupHandoffData (Grp : Type u) where
  /-- Effective mass over scale, `q(μ) = m(μ)/μ`. -/
  q : ℝ → ℝ
  /-- Running coupling squared. -/
  g2 : ℝ → ℝ
  /-- Effective transfer-matrix gap. -/
  gap : ℝ → ℝ
  /-- `q` at the weak-coupling boundary: `q0 = q μ_up`. -/
  q0 : ℝ
  /-- Strong-coupling threshold for `q`. -/
  q1 : ℝ
  /-- Weak-side gap bound. -/
  δW : ℝ
  /-- Strong-side gap bound. -/
  δS : ℝ
  q_strictAnti : StrictAntiOn q (Set.Ioi 0)
  g2_strictAnti : StrictAntiOn g2 (Set.Ioi 0)
  /-- The weak-coupling boundary scale. -/
  μup : ℝ
  μup_pos : 0 < μup
  gstar : ℝ
  gstar_pos : 0 < gstar
  g2_μup : g2 μup = gstar ^ 2
  /-- The stated value of `q0` really is the value at the boundary. -/
  q0_eq : q0 = q μup
  /-- The overlap hypothesis, in the group-parameterized form `q1 ≤ q0`. -/
  overlap : q1 ≤ q0

namespace GroupHandoffData

variable {Grp : Type u} (D : GroupHandoffData Grp)

/-- Weak-side input: below the coupling threshold, the gap is at least `δ_W`. -/
def WeakSide : Prop :=
  ∀ μ : ℝ, 0 < μ → D.g2 μ ≤ D.gstar ^ 2 → D.δW ≤ D.gap μ

/-- Strong-side input: once `q ≥ q₁`, the gap is at least `δ_S`, all the way to the IR. -/
def StrongSide : Prop :=
  ∀ μ : ℝ, 0 < μ → D.q1 ≤ D.q μ → D.δS ≤ D.gap μ

/-- **Group-parameterized handoff lemma.**  The same interval-cover argument as
`HandoffData.gap_uniform`, with the overlap hypothesis stated as `q1 ≤ q0` and `q0`
tied to `q μ_up` by a field. -/
theorem gap_uniform (hW : D.WeakSide) (hS : D.StrongSide) :
    ∀ μ : ℝ, 0 < μ → min D.δW D.δS ≤ D.gap μ := by
  intro μ hμ
  by_cases h : D.g2 μ ≤ D.gstar ^ 2
  · exact le_trans (min_le_left _ _) (hW μ hμ h)
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
    have hov : D.q1 ≤ D.q D.μup := by
      rw [← D.q0_eq]
      exact D.overlap
    have hq : D.q1 ≤ D.q μ :=
      le_trans hov (le_of_lt (D.q_strictAnti hμ D.μup_pos hlt))
    exact le_trans (min_le_right _ _) (hS μ hμ hq)

end GroupHandoffData

/-- **No-hidden-dependence certificate.**

Given a family of groups `Grp i` with handoff data `D i`, if the three constants that
appear in the conclusion (`q1`, `δW`, `δS`) are the same for all `i`, then the gap
bound is uniform in `i` — and the proof never mentions `Grp`.  This is the formal
content of "group-uniformity is localised in four numbers": the group enters the
argument only through the constants, never through the analysis.

(`q0` does not appear in the conclusion; it enters only via the overlap hypothesis,
which is why it is not among the equality hypotheses.) -/
theorem uniform_over_groups {ι : Type u} {Grp : ι → Type u}
    (D : ∀ i, GroupHandoffData (Grp i))
    (hW : ∀ i, (D i).WeakSide) (hS : ∀ i, (D i).StrongSide)
    {q1 δW δS : ℝ}
    (_hq1 : ∀ i, (D i).q1 = q1) (hWc : ∀ i, (D i).δW = δW)
    (hSc : ∀ i, (D i).δS = δS) :
    ∀ i (μ : ℝ), 0 < μ → min δW δS ≤ (D i).gap μ := by
  intro i μ hμ
  have h := (D i).gap_uniform (hW i) (hS i) μ hμ
  rw [hWc i, hSc i] at h
  exact h

/-- The constant ledger `(α, β, γ)` of the group-uniformity report §6: powers of
`dim G`, `h^∨` (dual Coxeter number) and `dim f` (dimension of the defining
representation) that a quantity is allowed to depend on.

It is a bookkeeping device, not a mathematical object; its purpose is to make
"this bound must not grow like `dim G`" a checkable annotation on statements. -/
structure Ledger where
  /-- Power of `dim G`. -/
  alpha : ℕ
  /-- Power of `h^∨`. -/
  beta : ℕ
  /-- Power of `dim f`. -/
  gamma : ℕ
deriving DecidableEq, Repr

/-- The ledger of an intensive quantity: allowed to depend only on `h^∨`. -/
def Ledger.intensive : Ledger := ⟨0, 1, 0⟩

/-- The ledger of a representation-dependent quantity. -/
def Ledger.representation : Ledger := ⟨0, 0, 1⟩

/-- The ledger of an extensive quantity: linear in `dim G`. -/
def Ledger.extensive : Ledger := ⟨1, 0, 0⟩

end RG
end YangMillsInfra
