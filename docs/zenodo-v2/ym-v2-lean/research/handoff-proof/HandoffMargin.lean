/-
# Handoff margin: the quantitative form of `q₀ ≥ q₁`

Companion to `REPORT.md` §1.2 / §7 and to
`research/lean-infra/YangMillsInfra/RG/Handoff.lean`.

The existing `Handoff.lean` proves the real-analytic covering lemma from the
abstract hypothesis `Overlap := q1 ≤ q μup`.  This file exposes the *constants*
of the one-loop handoff, so that the hypothesis becomes a single numerical
inequality.

## Convention (important)

This file uses the **physical normalization `T(f) = 1/2`**:

    b₀ = 11 h^∨ / (48 π²),    B = 2 b₀ = 11 h^∨ / (24 π²),
    d(1/g²)/d ln μ = B,       1/g²(μ) = B ln(μ/Λ),
    t_up = 1/(B g*²) = 24 π² / (11 λ*),     λ* := h^∨ g*².

SU(3) `g² = 1` gives `t_up = 7.178`, matching `numerics/REPORT.md`.

`research/group-uniformity/REPORT.md` uses `T(adj) = 2h^∨`, where the coupling
is **half** the physical one; its `λ*` is half and its exponent reads
`12 π² / (11 λ*)` in those units.  After conversion the two formulas agree.
Both are correct in their own normalization; the only mistake is mixing them.
Do not compare a `λ*` from one convention with a `t_up` from the other.
-/

import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Real.Basic

noncomputable section

namespace YangMillsInfra
namespace RG

/-- One-loop RG time of the weak-coupling boundary, `t_up = C_RG / λ*`.
The constant is a parameter so that a different scheme normalization can be
threaded through without changing the statements. -/
def tUp (C_RG lamStar : ℝ) : ℝ := C_RG / lamStar

/-- The invariant ratio at the weak boundary: `q₀ = (m/Λ) exp(−t_up)`. -/
def q0 (C_RG m Λ lamStar : ℝ) : ℝ :=
  (m / Λ) * Real.exp (-(tUp C_RG lamStar))

/-- The handoff margin `M = ln(m/Λ) − t_up − ln q₁`. -/
def margin (C_RG m Λ lamStar q1 : ℝ) : ℝ :=
  Real.log (m / Λ) - tUp C_RG lamStar - Real.log q1

/-- **`H3` in pure-number form.**  For positive data,
`q₁ ≤ q₀ ↔ 0 ≤ margin`. -/
theorem handoff_iff_margin (C_RG m Λ lamStar q1 : ℝ)
    (hm : 0 < m) (hΛ : 0 < Λ) (_hlam : 0 < lamStar) (hq : 0 < q1) :
    q1 ≤ q0 C_RG m Λ lamStar ↔ 0 ≤ margin C_RG m Λ lamStar q1 := by
  unfold q0 margin tUp
  have hmL : 0 < m / Λ := div_pos hm hΛ
  have hexp : 0 < Real.exp (-(C_RG / lamStar)) := Real.exp_pos _
  -- `q₁ ≤ (m/Λ) e^{−t}`  ↔  `log q₁ ≤ log(m/Λ) − t`  ↔  `M ≥ 0`
  rw [← Real.log_le_log_iff hq (mul_pos hmL hexp),
      Real.log_mul (ne_of_gt hmL) (ne_of_gt hexp),
      Real.log_exp]
  constructor <;> intro h <;> linarith

/-- **Inverse form.**  Given `q₁ < m/Λ`, the weak reach needed for `q₀ = q₁` is
`λ* = C_RG / log((m/Λ)/q₁)`. -/
theorem required_lambda (C_RG m Λ q1 : ℝ)
    (hC : 0 < C_RG) (hm : 0 < m) (hΛ : 0 < Λ) (hq : 0 < q1) (h : q1 < m / Λ) :
    margin C_RG m Λ (C_RG / Real.log ((m / Λ) / q1)) q1 = 0 := by
  unfold margin tUp
  have hmL : 0 < m / Λ := div_pos hm hΛ
  have hratio : 1 < (m / Λ) / q1 := by
    rw [lt_div_iff₀ hq]   -- `lt_div_iff` in older Mathlib
    linarith
  have hlogpos : 0 < Real.log ((m / Λ) / q1) := Real.log_pos hratio
  have hne : Real.log ((m / Λ) / q1) ≠ 0 := ne_of_gt hlogpos
  have hCne : C_RG ≠ 0 := ne_of_gt hC
  have hdiv : C_RG / (C_RG / Real.log ((m / Λ) / q1))
      = Real.log ((m / Λ) / q1) := by
    field_simp [hCne, hne]
  rw [hdiv, Real.log_div (ne_of_gt hmL) (ne_of_gt hq)]
  ring

/-! ## Certificate type (do not collapse)

The handoff interface must keep the certificate a **function of the scale**:

    (cert : ℝ → ℝ)
    (hW : ∀ μ, 0 < μ → g2 μ ≤ gstar^2 → δW ≤ cert μ)
    (hS : ∀ μ, 0 < μ → q1 ≤ q μ     → δS ≤ cert μ)
    (hOv : q1 ≤ q0)
    :  ∀ μ, 0 < μ → min δW δS ≤ cert μ

That is `Handoff.lean`'s `gap_uniform` (with `gap : ℝ → ℝ`).  `gluing` §2.7
has been corrected to this `cert : ℝ → ℝ` shape (`mass_gap_of_certificate`);
the constant-`gapBar` caveat is historical, kept only as an interface rule.
-/

/-- The canonical constant in the physical normalization: `C_RG = 24π²/11`. -/
def C_RG : ℝ := 24 * Real.pi ^ 2 / 11

/-- Sanity: the SU(3) one-loop value at `g² = 1` is `t_up = 7.178…` -/
example : tUp C_RG (3 * 1) = 24 * Real.pi ^ 2 / 33 := by
  unfold tUp C_RG
  ring

end RG
end YangMillsInfra

/-! ## Axiom audit (build-enforced)

These guards make `bash lean/check-handoff.sh` fail if the axiom footprint of the
two headline theorems changes (e.g. a `sorry` leaks in). -/

/-- info: 'YangMillsInfra.RG.handoff_iff_margin' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YangMillsInfra.RG.handoff_iff_margin

/-- info: 'YangMillsInfra.RG.required_lambda' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YangMillsInfra.RG.required_lambda
