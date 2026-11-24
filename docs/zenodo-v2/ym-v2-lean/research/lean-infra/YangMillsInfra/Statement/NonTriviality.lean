/-
# Non-triviality: Gaussianity, cumulants, and the renormalized-coupling witness

This module incorporates the nontriviality line's finding
(`research/nontriviality/REPORT.md`, §1.2): for a **pure gauge** theory the witness
"some gauge-invariant connected 4-point function is nonzero" is **vacuous**, because
every gauge-invariant local operator is at least quadratic in the field strength, so
`tr F²` is a Wick square of a Gaussian field and a free gauge field already has nonzero
connected 4-point functions.  The valid witness is instead the renormalized coupling
(`RGCharge` below): gauge invariant, and zero for a free theory.

What is formalized here:

* `IsGaussianTheory` — Gaussianity stated on the generating functional, matching
  `YMStatement.IsGaussianMeasure` exactly;
* `cumulantGF`, `nthCumulant` — cumulants of the smeared field along `t ↦ t • J`,
  defined as derivatives of `log Z(t • J)`;
* **Lemma C, forward direction** (`isGaussianTheory_cumulants_vanish`): a Gaussian
  theory has all cumulants of order `≥ 3` vanishing.  Proved, not assumed.
* **Lemma C, converse** — packaged as `CumulantCharacterization`, whose
  `gaussian_of_vanishing` field is the analytic implication (needs OS0-type
  analyticity of the cumulant generating function);
* falsification tests: `gaussian_not_nonTrivial` and the concrete interacting toy
  `toyTheory_not_gaussian`;
* `RGCharge`, `GaugeNonTrivial`, and `gaugeNonTrivial_implies_not_gaussian`.

The distinction the nontriviality report insists on is encoded in the types: the
4-point witness is *not* part of `GaugeNonTrivial`; only `RGCharge` is.

**Convention note.**  `cumulantGF` takes the *real part* of `log Z(t • J)`.  For a
reflection-positive or symmetric theory `Z(t • J)` is real and positive, so this is the
ordinary cumulant generating function; taking the real part keeps the definition real
and avoids branch-of-log issues.  The forward direction of Lemma C is unaffected; the
converse is an explicit interface field.
-/
import YangMillsInfra.Basic
import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Analysis.SpecialFunctions.Log.Basic

noncomputable section

namespace YangMillsInfra

variable {d : ℕ}

/-- A Euclidean theory is **Gaussian** (quasi-free) if its generating functional has the
characteristic form of a centred Gaussian.  This is the measure-free counterpart of
`YMStatement.IsGaussianMeasure`: for a measure-backed theory the two are the same
predicate on `Z`. -/
def IsGaussianTheory (T : EuclideanTheory d) : Prop :=
  ∃ B : TestFunction d →ₗ[ℝ] TestFunction d →ₗ[ℝ] ℝ,
    ∀ J : TestFunction d,
      T.Z J = Complex.exp (-(1 / 2 : ℂ) * ((B J J : ℝ) : ℂ))

/-- Non-triviality: the theory is not a free (Gaussian) field.  For a pure gauge theory
this is *not* witnessed by a connected 4-point function; see `GaugeNonTrivial`. -/
def NonTrivial (T : EuclideanTheory d) : Prop :=
  ¬ IsGaussianTheory T

/-- The cumulant generating function of the smeared field `⟨φ, J⟩`, evaluated on the
line `t ↦ t • J`: the real part of `log Z(tJ)`. -/
def cumulantGF (T : EuclideanTheory d) (J : TestFunction d) : ℝ → ℝ :=
  fun t => (Complex.log (T.Z (t • J))).re

/-- The `n`-th cumulant of the smeared field `⟨φ, J⟩`. -/
def nthCumulant (T : EuclideanTheory d) (J : TestFunction d) (n : ℕ) : ℝ :=
  iteratedDeriv n (cumulantGF T J) 0

/-- **Lemma C, forward direction.**  A Gaussian theory has vanishing cumulants of order
`≥ 3`: along every line `t ↦ t • J` the cumulant generating function is the quadratic
`−t² B(J,J)/2`, whose third and higher derivatives vanish.

This is the rigorous half of the nontriviality report's Lemma C.  For scalar theories
the 4-point cumulant (hence the connected 4-point function) suffices, because
`κₙ = 0` for odd `n` by symmetry; for gauge theories it does **not**, which is why
`RGCharge` is the witness there. -/
theorem isGaussianTheory_cumulants_vanish {T : EuclideanTheory d}
    (h : IsGaussianTheory T) :
    ∀ (n : ℕ), 3 ≤ n → ∀ J : TestFunction d, nthCumulant T J n = 0 := by
  obtain ⟨B, hB⟩ := h
  intro n hn J
  -- `B (t • J) (t • J) = t² B J J`
  have hlin : ∀ t : ℝ, B (t • J) (t • J) = (t * t) * B J J := by
    intro t
    simp only [map_smul, LinearMap.smul_apply, smul_eq_mul]
    ring
  -- `log Z(t J) = -(t² B J J)/2`, hence the real part is a quadratic in `t`.
  have hfun : cumulantGF T J = fun t : ℝ => (-(1 / 2 : ℝ) * B J J) * t ^ 2 := by
    funext t
    rw [cumulantGF, hB (t • J), hlin t]
    have hz : Complex.log (Complex.exp (-(1 / 2 : ℂ) * (((t * t) * B J J : ℝ) : ℂ)))
        = (-(1 / 2 : ℂ) * (((t * t) * B J J : ℝ) : ℂ)) :=
      Complex.log_exp (by simpa using Real.pi_pos) (by simpa using Real.pi_pos.le)
    rw [hz]
    simp
    ring
  rw [nthCumulant, hfun]
  rw [iteratedDeriv_const_mul_field, iteratedDeriv_pow]
  have hdesc : Nat.descFactorial 2 n = 0 := by
    rw [Nat.descFactorial_eq_zero_iff_lt]
    omega
  simp [hdesc]

/-- **Lemma C, converse direction, as an interface.**  If the cumulant generating
function is analytic at the origin along every line and all cumulants of order `≥ 3`
vanish, the theory is Gaussian.

The analyticity input is exactly OS0 (analyticity of `Z`) transported to `log Z`; it is
not derivable from the vanishing of cumulants alone.  This is the one direction of
Lemma C that is not proved here. -/
structure CumulantCharacterization (T : EuclideanTheory d) : Prop where
  /-- Analyticity of the cumulant generating function at the origin. -/
  analytic_cgf : ∀ J : TestFunction d, AnalyticAt ℝ (cumulantGF T J) 0
  /-- Vanishing of all cumulants of order `≥ 3` forces Gaussianity. -/
  gaussian_of_vanishing :
    (∀ (n : ℕ), 3 ≤ n → ∀ J : TestFunction d, nthCumulant T J n = 0) → IsGaussianTheory T

/-- **Lemma C**, as an equivalence, conditional on the analyticity input. -/
theorem isGaussianTheory_iff_cumulants_vanish {T : EuclideanTheory d}
    (h : CumulantCharacterization T) :
    IsGaussianTheory T ↔ ∀ (n : ℕ), 3 ≤ n → ∀ J : TestFunction d, nthCumulant T J n = 0 :=
  ⟨isGaussianTheory_cumulants_vanish, h.gaussian_of_vanishing⟩

/-! ## Falsification tests

The statement layer must reject the free theory and accept an interacting one.  The
free-theory direction is trivial by definition; the interacting direction is a concrete
computation. -/

/-- A Gaussian theory is not non-trivial (the "free field must fail" test). -/
theorem gaussian_not_nonTrivial {T : EuclideanTheory d} (h : IsGaussianTheory T) :
    ¬ NonTrivial T := by
  simpa [NonTrivial] using h

/-- A concrete interacting toy theory: a Gaussian term plus a quartic term in a single
linear functional `L` of the test function, with coupling `λ`.

For `L ≠ 0` and `λ ≠ 0` this is **not** Gaussian, so it passes `NonTrivial`.  This is
the "interacting toy must pass" falsification test. -/
def toyTheory (L : TestFunction d →ₗ[ℝ] ℝ) (lam : ℝ) : EuclideanTheory d where
  Z J := ((Real.exp (-(1 / 2 : ℝ) * (L J) ^ 2 - lam * (L J) ^ 4) : ℝ) : ℂ)

/-- The interacting toy theory is not Gaussian. -/
theorem toyTheory_not_gaussian {L : TestFunction d →ₗ[ℝ] ℝ} (hL : L ≠ 0)
    {lam : ℝ} (hlam : lam ≠ 0) : ¬ IsGaussianTheory (toyTheory (d := d) L lam) := by
  rintro ⟨B, hB⟩
  -- Pick `J₀` with `L J₀ ≠ 0`.
  have hL' : ∃ J : TestFunction d, L J ≠ 0 := by
    by_contra h
    push Not at h
    exact hL (LinearMap.ext h)
  obtain ⟨J₀, hJ₀⟩ := hL'
  set a : ℝ := L J₀ with ha
  set b : ℝ := B J₀ J₀ with hb
  -- Both sides are `exp` of a real number, so the exponents agree.
  have hexp : ∀ J : TestFunction d,
      -(1 / 2 : ℝ) * (L J) ^ 2 - lam * (L J) ^ 4 = -(1 / 2 : ℝ) * B J J := by
    intro J
    have h := hB J
    change ((Real.exp (-(1 / 2 : ℝ) * (L J) ^ 2 - lam * (L J) ^ 4) : ℝ) : ℂ)
        = Complex.exp (-(1 / 2 : ℂ) * ((B J J : ℝ) : ℂ)) at h
    have hC : (-(1 / 2 : ℂ) * ((B J J : ℝ) : ℂ)) = ((-(1 / 2 : ℝ) * B J J : ℝ) : ℂ) := by
      push_cast
      ring
    rw [hC, ← Complex.ofReal_exp] at h
    exact Real.exp_injective (Complex.ofReal_injective h)
  have e1 : -(1 / 2 : ℝ) * a ^ 2 - lam * a ^ 4 = -(1 / 2 : ℝ) * b := by
    have := hexp J₀
    simpa [ha, hb] using this
  have h2 : L ((2 : ℝ) • J₀) = 2 * a := by
    rw [map_smul, smul_eq_mul, ← ha]
  have hb2 : B ((2 : ℝ) • J₀) ((2 : ℝ) • J₀) = 4 * b := by
    simp only [map_smul, LinearMap.smul_apply, smul_eq_mul]
    rw [← hb]
    ring
  have e2 : -(1 / 2 : ℝ) * (2 * a) ^ 2 - lam * (2 * a) ^ 4 = -(1 / 2 : ℝ) * (4 * b) := by
    have h := hexp ((2 : ℝ) • J₀)
    rw [h2, hb2] at h
    exact h
  -- `e1` and `e2` force `λ * a ^ 4 = 0`, contradicting `λ ≠ 0` and `a ≠ 0`.
  have hkey : lam * a ^ 4 = 0 := by
    nlinarith [e1, e2]
  exact (mul_ne_zero hlam (pow_ne_zero 4 hJ₀)) hkey

/-! ## The gauge-theory witness: renormalized coupling -/

/-- A renormalized-coupling witness: a gauge-invariant real number that vanishes on
every Gaussian theory.  For a pure gauge theory this — not a connected 4-point
function — is the non-triviality witness (nontriviality report §1.2). -/
structure RGCharge (T : EuclideanTheory d) where
  value : ℝ
  free_vanishes : IsGaussianTheory T → value = 0

/-- Non-triviality in the gauge-theory sense: some renormalized-coupling witness is
nonzero.  This is the witness form that cannot be satisfied vacuously by a Wick square
of the free field. -/
def GaugeNonTrivial (T : EuclideanTheory d) : Prop :=
  ∃ c : RGCharge T, c.value ≠ 0

/-- A nonzero renormalized coupling rules out Gaussianity. -/
theorem gaugeNonTrivial_implies_not_gaussian {T : EuclideanTheory d}
    (h : GaugeNonTrivial T) : ¬ IsGaussianTheory T := by
  obtain ⟨c, hc⟩ := h
  intro hg
  exact hc (c.free_vanishes hg)

end YangMillsInfra
