/-
# Minimal handoff: the RG-time comparison (R3)

This file formalizes only the *one genuine reduction* of `MINIMAL.md` §3:
if `bet ≤ bet2 < 0` pointwise on `[lam, lamc]`, then the true RG time

    t_up(lam) = ∫_{lam}^{lamc} dx / (−bet(x))

is at most the reference (2-loop) RG time `∫_{lam}^{lamc} dx / (−bet2(x))`.

The other reductions of `MINIMAL.md` (R1 drop sigma, R2 uniqueness, R4 spectral
decay) are standard logic / spectral theory and are proved in the prose there;
they are not repeated here.

The point of the lemma: the *quantitative* remainder bound
`|bet − bet2| ≤ C3 * x^4` (NPB) is not needed — a one-sided sign condition
suffices.
-/

import Mathlib

noncomputable section
open MeasureTheory intervalIntegral

namespace YMStatement.MinimalHandoff

/-- **R3.** A one-sided sign condition on `bet` bounds the true RG time by the
reference RG time.  Integrability is an explicit hypothesis (`hint1`, `hint2`);
in the application `bet` and `bet2` are smooth beta functions on `(0, lamc]`,
so these hold automatically. -/
theorem rgTime_le_of_beta_le {lamc lam : ℝ} {bet bet2 : ℝ → ℝ}
    (_hlam : 0 < lam) (hllam : lam ≤ lamc)
    (hle : ∀ x ∈ Set.Icc lam lamc, bet x ≤ bet2 x)
    (hbet2 : ∀ x ∈ Set.Icc lam lamc, bet2 x < 0)
    (hint1 : IntervalIntegrable (fun x => 1 / (-bet x)) volume lam lamc)
    (hint2 : IntervalIntegrable (fun x => 1 / (-bet2 x)) volume lam lamc) :
    ∫ x in lam..lamc, (1 / (-bet x)) ≤ ∫ x in lam..lamc, (1 / (-bet2 x)) := by
  refine intervalIntegral.integral_mono_on hllam hint1 hint2 ?_
  intro x hx
  -- Both `−bet2 x` and `−bet x` are positive, with `−bet2 x ≤ −bet x`.
  have h2 : 0 < -bet2 x := by linarith [hbet2 x hx]
  have h1 : -bet2 x ≤ -bet x := by linarith [hle x hx]
  -- Hence `1 / (−bet x) ≤ 1 / (−bet2 x)`.
  exact one_div_le_one_div_of_le h2 h1

/-! ## One scale, one side

The handoff proof only ever reads the certificate at a **single** scale
(`mass_gap_of_certificate` evaluates at `μ = 1`), and at that single scale only
**one** of the weak / strong hypotheses is used (the `by_cases` in
`handoff_gap`).  So neither *uniformity in the scale* nor *both sides* is
needed: one certificate at one scale suffices.  These two lemmas record that
logical core; the source of the disjunct (W or S at that scale) is supplied by
whichever analysis is available there. -/

/-- **One side, one scale.** A single lower bound at a single scale suffices. -/
theorem min_le_of_one_certificate {dW dS c : ℝ} (h : dW ≤ c ∨ dS ≤ c) :
    min dW dS ≤ c := by
  rcases h with h | h
  · exact le_trans (min_le_left _ _) h
  · exact le_trans (min_le_right _ _) h

/-- **One scale + bridge.** With the certificate bounded above by `Delta` at the
same single scale, the gap follows. -/
theorem gap_of_one_certificate {dW dS c Delta : ℝ}
    (h : dW ≤ c ∨ dS ≤ c) (hc : c ≤ Delta) :
    min dW dS ≤ Delta :=
  le_trans (min_le_of_one_certificate h) hc

end YMStatement.MinimalHandoff

/-- info: 'YMStatement.MinimalHandoff.rgTime_le_of_beta_le' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMStatement.MinimalHandoff.rgTime_le_of_beta_le

/-- info: 'YMStatement.MinimalHandoff.min_le_of_one_certificate' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMStatement.MinimalHandoff.min_le_of_one_certificate

/-- info: 'YMStatement.MinimalHandoff.gap_of_one_certificate' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMStatement.MinimalHandoff.gap_of_one_certificate
