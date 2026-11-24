/-
GL7Dual.lean — the real-analysis core of glue lemma GL7 (dual handoff, one RT).

GL7 evaluates `I-HANDOFF[m]` and `I-HANDOFF[σ]` on the same renormalization
trajectory.  With a shared weak boundary `μ_up` (hence a shared RG time
`t = t_up`), each certificate has the numerical-overlap condition

    H3_i : q1_i ≤ q0_i,     q0_i = (invariant_i/Λ) · exp(−t_up).

`dual_overlap_iff` proves the collapse: both overlaps hold iff a single
inequality holds, with the maximum of the two thresholds,

    R_m ≥ exp(t) · max(q1^m, (m/√σ)·q1^σ),

where `R_m = m/Λ`, `R_σ = √σ/Λ` and `c := R_m/R_σ = m/√σ`.

`dual_overlap_iff_muStar` gives the falsification form: `H3_i ↔ μ_up ≤ μ*_i`
with `μ*_i = invariant_i/q1_i`, so when the two `μ*` differ the smaller one
decides — coverage is `μ_up ≤ min(μ*_m, μ*_σ)`.

Inputs consumed (all as hypotheses, none proved here):
  * the shared RT and shared `t_up` (GL7-a, an interface obligation);
  * `CPB-W` (weak side, ≥2-loop `t_up`), `CPB-S^m`, `CPB-S^σ` (strong sides);
  * `C = C_RG` (the shared scheme constant).
The open estimates themselves are NOT formalized here.

No `sorry`; `#print axioms` must stay `[propext, Classical.choice, Quot.sound]`.
-/
import Mathlib

namespace GlueD

/-- Handoff coordinate at the weak boundary: `q₀ = R · e^{−t}`. -/
noncomputable def q0 (R t : ℝ) : ℝ := R * Real.exp (-t)

/-- Numerical overlap `H3` of one certificate. -/
def H3 (R t q1 : ℝ) : Prop := q1 ≤ q0 R t

/-- `q₀ · e^{t} = R`. -/
lemma q0_mul_exp (R t : ℝ) : q0 R t * Real.exp t = R := by
  unfold q0
  have h : Real.exp (-t) * Real.exp t = 1 := by
    rw [← Real.exp_add, neg_add_cancel, Real.exp_zero]
  calc R * Real.exp (-t) * Real.exp t
      = R * (Real.exp (-t) * Real.exp t) := by ring
    _ = R := by rw [h, mul_one]

/-- `x ≤ q₀` iff `x · e^{t} ≤ R`. -/
lemma le_q0_iff (x R t : ℝ) (h : 0 < Real.exp t) :
    x ≤ q0 R t ↔ x * Real.exp t ≤ R := by
  constructor
  · intro hx
    calc x * Real.exp t ≤ q0 R t * Real.exp t :=
          mul_le_mul_of_nonneg_right hx (le_of_lt h)
      _ = R := q0_mul_exp R t
  · intro hx
    have hexp : 0 < Real.exp (-t) := Real.exp_pos _
    have h1 : x * Real.exp t * Real.exp (-t) ≤ R * Real.exp (-t) :=
      mul_le_mul_of_nonneg_right hx (le_of_lt hexp)
    have h2 : x * Real.exp t * Real.exp (-t) = x := by
      rw [mul_assoc, ← Real.exp_add, add_neg_cancel, Real.exp_zero, mul_one]
    unfold q0
    rw [← h2]
    exact h1

/-- **Collapse.** With one shared RG time, the two overlaps are equivalent to
one inequality whose threshold is the maximum of the two.  `c = m/√σ` is the
ratio of the invariants, expressed by `hRm : R_m = c · R_σ`. -/
theorem dual_overlap_iff (Rm Rs c t q1m q1s : ℝ)
    (hRm : Rm = c * Rs) (hc : 0 < c) :
    (H3 Rm t q1m ∧ H3 Rs t q1s) ↔
      Real.exp t * max q1m (c * q1s) ≤ Rm := by
  have he : 0 < Real.exp t := Real.exp_pos t
  have hm := le_q0_iff q1m Rm t he
  have hs := le_q0_iff q1s Rs t he
  unfold H3
  rw [hm, hs]
  have hassoc : c * (q1s * Real.exp t) = (c * q1s) * Real.exp t := by ring
  have hcs : q1s * Real.exp t ≤ Rs ↔ (c * q1s) * Real.exp t ≤ Rm := by
    constructor
    · intro h
      have := mul_le_mul_of_nonneg_left h (le_of_lt hc)
      simpa [hRm, hassoc] using this
    · intro h
      have h' : c * (q1s * Real.exp t) ≤ c * Rs := by
        rw [hassoc]; simpa [hRm] using h
      exact le_of_mul_le_mul_left h' hc
  rw [hcs]
  have hmax : Real.exp t * max q1m (c * q1s)
      = max (q1m * Real.exp t) ((c * q1s) * Real.exp t) := by
    rcases le_total q1m (c * q1s) with hle | hle
    · rw [max_eq_right hle,
          max_eq_right (mul_le_mul_of_nonneg_right hle (le_of_lt he))]
      ring
    · rw [max_eq_left hle,
          max_eq_left (mul_le_mul_of_nonneg_right hle (le_of_lt he))]
      ring
  rw [hmax]
  exact max_le_iff.symm

/-- `H3` is exactly `μ_up ≤ μ*`, with `μ* = R/q1` in units of `Λ`. -/
theorem overlap_iff_muStar (R t q1 : ℝ) (hq : 0 < q1) :
    H3 R t q1 ↔ Real.exp t ≤ R / q1 := by
  have he : 0 < Real.exp t := Real.exp_pos t
  unfold H3
  rw [le_q0_iff q1 R t he, le_div_iff₀ hq]
  constructor
  · intro h
    linarith [mul_comm q1 (Real.exp t)]
  · intro h
    linarith [mul_comm q1 (Real.exp t)]

/-- **Falsification form.** The dual overlap holds iff the covered strong
region reaches the weak boundary for the smaller of the two `μ*`.  When the two
`μ*` differ, the smaller decides: `H3_m ∧ H3_σ ↔ μ_up ≤ min(μ*_m, μ*_σ)`. -/
theorem dual_overlap_iff_muStar (Rm Rs t q1m q1s : ℝ)
    (hq1m : 0 < q1m) (hq1s : 0 < q1s) :
    (H3 Rm t q1m ∧ H3 Rs t q1s) ↔
      Real.exp t ≤ min (Rm / q1m) (Rs / q1s) := by
  rw [overlap_iff_muStar Rm t q1m hq1m,
      overlap_iff_muStar Rs t q1s hq1s, le_min_iff]

end GlueD

/-! ## Axiom audit (build-enforced) -/

/-- info: 'GlueD.dual_overlap_iff' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueD.dual_overlap_iff

/-- info: 'GlueD.dual_overlap_iff_muStar' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms GlueD.dual_overlap_iff_muStar
