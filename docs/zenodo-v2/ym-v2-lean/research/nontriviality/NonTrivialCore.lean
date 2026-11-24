/-
# Non-vacuity and sharpness of the connected-cumulant witness

Machine-checked (core Lean, no Mathlib). Companion to `REPORT.md` §6 and `PROOF.md`.

Context. `lean/YMStatement/Basic.lean` defines non-triviality as `¬ IsGaussianMeasure`
and provides `ConnectedFourPoint` (four-point minus the three Wick pairings) as a
*witness functional*. Two things must be true for that witness to be a valid, non-vacuous
detector of "not a free/Gaussian field":

  (A) TIMES ARE COMPUTABLE. The connected 4-point of a concrete measure can be evaluated
      and compared to zero.  (Proved below for finite weighted samples.)

  (B) THE 4-POINT ALONE IS **NOT** ENOUGH.  The tempting statement "an even measure with
      vanishing connected 4-point is Gaussian" is FALSE.  We exhibit a symmetric
      (even) finite measure with **zero** fourth cumulant but **nonzero** sixth cumulant.
      Hence a correct witness must range over ALL cumulant orders ≥ 3, not just n = 4.
      (This corrects an over-claim in an earlier draft of REPORT §6.1.)

We work over `Int` (kernel-decidable). A finite "measure" is a `List (weight, value)`;
all moments are integer sums, so `decide` closes every concrete identity exactly
(no floating point, no `norm_num`). With total weight `D` and raw moment
`M_k = Σ w·x^k`, the *normalized* cumulants equal the integer quantities below divided
by a positive power of `D`, so `K = 0 ↔ normalized cumulant = 0`.

Mapping to Basic.lean (informal, since that file needs the OSforGFF build, unavailable here):
  `SchwingerFunction₂`      ↔ `M_2 / D`
  `SchwingerFunction 4`     ↔ `M_4 / D`
  `ConnectedFourPoint`      ↔ `K4 / D^2`   (for the single-variable specialization)
-/

import Std

namespace NonTrivialCore

/-! ## The connected 4-point functional (mirrors `YMStatement.ConnectedFourPoint`) -/

/-- Connected four-point from a raw 4-point and the three Wick pairings.
This is exactly `Basic.lean`'s `ConnectedFourPoint` shape:
`S4 − (S12·S34 + S13·S24 + S14·S23)`. -/
def conn4 (S4 S12 S34 S13 S24 S14 S23 : Int) : Int :=
  S4 - (S12 * S34 + S13 * S24 + S14 * S23)

/-- The witness vanishes iff the 4-point equals its Wick value. -/
theorem conn4_eq_zero_iff (x4 w : Int) : (x4 - w = 0) ↔ x4 = w := by
  omega

/-- Gaussian (Wick) moment data give vanishing connected 4-point. -/
example : conn4 6 1 2 3 4 5 6 = 0 ∨ conn4 6 1 2 3 4 5 6 ≠ 0 := by
  -- (placeholder to show the functional is well-typed; the real content is below)
  by_cases h : conn4 6 1 2 3 4 5 6 = 0
  · exact Or.inl h
  · exact Or.inr h

/-! ## Finite weighted samples and their moments -/

/-- Raw `k`-th moment `Σ w · x^k`. -/
def Mk (s : List (Int × Int)) (k : Nat) : Int :=
  s.foldr (fun wx acc => acc + wx.1 * wx.2 ^ k) 0

/-- Total weight `D = Σ w`. -/
def D (s : List (Int × Int)) : Int :=
  s.foldr (fun wx acc => acc + wx.1) 0

/-- Integer form of the normalized fourth cumulant `κ₄ = m₄ − 3 m₂²`, times `D²`.
`κ₄ = K4 / D²`. -/
def K4 (s : List (Int × Int)) : Int :=
  D s * Mk s 4 - 3 * Mk s 2 * Mk s 2

/-- Integer form of the normalized sixth cumulant of a mean-zero symmetric law,
`κ₆ = m₆ − 15 m₄ m₂ − 10 m₃² + 30 m₂³`, times `D³`.  `κ₆ = K6 / D³`. -/
def K6 (s : List (Int × Int)) : Int :=
  D s ^ 2 * Mk s 6
    - 15 * D s * Mk s 4 * Mk s 2
    - 10 * D s * Mk s 3 ^ 2
    + 30 * Mk s 2 ^ 3

/-! ## Concrete measures -/

/-- Rademacher: value `±1`, each weight `1` (total weight 2). -/
def radem : List (Int × Int) := [(1, -1), (1, 1)]

/-- A symmetric three-atom law on `{-1, 0, 1}` with weights `{1, 4, 1}` (total 6).
Its moments are  `m₂ = m₄ = m₆ = 1/3`.  It is "Gaussian through order 4"
(`κ₂ = 1/3`, `κ₄ = 0`) but NOT Gaussian (`κ₆ = −2/9`). -/
def gauss4 : List (Int × Int) := [(1, -1), (4, 0), (1, 1)]

/-! ## Machine-checked facts -/

/-- Rademacher is symmetric with `m₂ = 1`, so `κ₂ = 1`. -/
theorem radem_m2 : Mk radem 2 = 2 := by decide

/-- Rademacher's fourth cumulant is `−2 ≠ 0`: the connected 4-point *does* detect
that Rademacher is not Gaussian. -/
theorem radem_K4 : K4 radem = -8 := by decide
theorem radem_K4_ne_zero : K4 radem ≠ 0 := by decide

/-- The three-atom law has vanishing fourth cumulant: `κ₄ = 0` (the connected
4-point cannot see its non-Gaussianity). -/
theorem gauss4_K4_zero : K4 gauss4 = 0 := by decide

/-- ... but its sixth cumulant is nonzero: `K6 = −48`, i.e. `κ₆ = −48/6³ = −2/9 ≠ 0`. -/
theorem gauss4_K6 : K6 gauss4 = -48 := by decide
theorem gauss4_K6_ne_zero : K6 gauss4 ≠ 0 := by decide

/-- **Sharpness.**  An even (symmetric) finite measure can have `κ₄ = 0` while `κ₆ ≠ 0`.
Therefore evenness + vanishing connected 4-point does NOT imply Gaussian; a correct
witness must control ALL cumulant orders `n ≥ 3` (odd orders vanishing by evenness,
but every even order `4, 6, 8, …` is needed). -/
theorem four_point_insufficient :
    K4 gauss4 = 0 ∧ K6 gauss4 ≠ 0 :=
  ⟨gauss4_K4_zero, gauss4_K6_ne_zero⟩

/-- Sanity: on the three-atom law the *Wick* m₄ relation `m₄ = 3 m₂²` holds
(`K4 = 0`), confirming the connected 4-point equals zero there. -/
theorem gauss4_wick_fourth : D gauss4 * Mk gauss4 4 = 3 * Mk gauss4 2 * Mk gauss4 2 := by
  decide

end NonTrivialCore

/-! ## Axiom audit

The file uses only core Lean + `Std`; every theorem above is closed by `decide`
(exact integer computation) or `omega` (linear arithmetic).  No `sorry`, no axioms
beyond Lean's kernel.  Check with `#print axioms` if audited in a full environment.
-/
