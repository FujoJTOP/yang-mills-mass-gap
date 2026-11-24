import Mathlib

/-!
# Parameterized RG handoff lemma

This file isolates the purely real-analytic core of the weak→strong coupling
handoff.  The gauge group enters only through the four constants
`q0, q1, deltaW, deltaS` and the two monotone functions `q` and `gsq`; the
conclusion contains no group data.  This is the "no hidden G dependence"
certificate described in `research/group-uniformity/REPORT.md` §9.2.

`q μ = m(μ)/μ` is the scheme-independent correlation-length ratio and
`gsq μ` is the running coupling; both are antitone on `(0, ∞)`.  `μ_up` is the
weak-side boundary where `gsq` reaches its small-coupling radius, and
`q0 = q μ_up`.  `q1` is the strong-side threshold.  The overlap condition
`q1 ≤ q0` is exactly the open handoff hypothesis.
-/

namespace YMStatement

/-- Data for the handoff. -/
structure HandoffData where
  q : ℝ → ℝ
  gsq : ℝ → ℝ
  μ_up : ℝ
  q0 : ℝ
  q1 : ℝ
  deltaW : ℝ
  deltaS : ℝ
  μ_up_pos : 0 < μ_up
  q_up : q μ_up = q0
  overlap : q1 ≤ q0
  q_antitone : ∀ a b : ℝ, 0 < a → a ≤ b → q b ≤ q a
  gsq_antitone : ∀ a b : ℝ, 0 < a → a ≤ b → gsq b ≤ gsq a

/-- The handoff lemma.

`cert μ` is the *certificate*: a lower bound on the **same fixed** dimensionless
number `Δ/Λ` (or `√σ/Λ`), obtained by whichever analysis is available at scale
`μ`.  It is **not** the physical gap at scale `μ` — that would be `m/μ`, and
demanding `δ_W ≤ m/μ` as `μ → ∞` would be unsatisfiable.  The semantic link
`cert μ ≤ Δ/Λ` is the RG gap-transfer lemma and belongs to the instantiation
(see `mass_gap_of_certificate`), not to this skeleton.

If the weak analysis covers `[μ_up, ∞)`, the strong analysis covers
`{μ | q μ ≥ q1}`, and the domains overlap (`q1 ≤ q0`), then the certificate is
bounded below by `min deltaW deltaS` at every scale. -/
theorem handoff_gap (D : HandoffData) (cert : ℝ → ℝ)
    (hW : ∀ μ : ℝ, 0 < μ → D.gsq μ ≤ D.gsq D.μ_up → D.deltaW ≤ cert μ)
    (hS : ∀ μ : ℝ, 0 < μ → D.q1 ≤ D.q μ → D.deltaS ≤ cert μ) :
    ∀ μ : ℝ, 0 < μ → min D.deltaW D.deltaS ≤ cert μ := by
  intro μ hμ
  by_cases hle : μ ≤ D.μ_up
  · -- Strong side: `μ ≤ μ_up`, so `q μ ≥ q μ_up = q0 ≥ q1`.
    have hq : D.q1 ≤ D.q μ := by
      calc
        D.q1 ≤ D.q0 := D.overlap
        _ = D.q D.μ_up := D.q_up.symm
        _ ≤ D.q μ := D.q_antitone μ D.μ_up hμ hle
    exact le_trans (min_le_right _ _) (hS μ hμ hq)
  · -- Weak side: `μ > μ_up`, so `gsq μ ≤ gsq μ_up`.
    have hlt : D.μ_up < μ := lt_of_not_ge hle
    have hg : D.gsq μ ≤ D.gsq D.μ_up :=
      D.gsq_antitone D.μ_up μ D.μ_up_pos (le_of_lt hlt)
    exact le_trans (min_le_left _ _) (hW μ hμ hg)

/-- The bridge from the certificate to the physical dimensionless gap.

If every certificate value is a lower bound on the fixed number `Δ_over_Λ`
(this is the RG gap-transfer statement: the gap is invariant along an exact
renormalized trajectory, and controlled along a truncated one), then the
handoff lemma yields `min deltaW deltaS ≤ Δ_over_Λ`.  The hypothesis `hcert` is
exactly the part that is *not* free. -/
theorem mass_gap_of_certificate (D : HandoffData) (cert : ℝ → ℝ) (Δ_over_Λ : ℝ)
    (hW : ∀ μ : ℝ, 0 < μ → D.gsq μ ≤ D.gsq D.μ_up → D.deltaW ≤ cert μ)
    (hS : ∀ μ : ℝ, 0 < μ → D.q1 ≤ D.q μ → D.deltaS ≤ cert μ)
    (hcert : ∀ μ : ℝ, 0 < μ → cert μ ≤ Δ_over_Λ) :
    min D.deltaW D.deltaS ≤ Δ_over_Λ := by
  have h1 : min D.deltaW D.deltaS ≤ cert 1 := handoff_gap D cert hW hS 1 zero_lt_one
  exact le_trans h1 (hcert 1 zero_lt_one)

end YMStatement
