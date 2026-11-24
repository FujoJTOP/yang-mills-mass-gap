/-
  GL2b' — gap continuity across the UV limit (REPORT.md §3.8).

  Paper statement:

      SPEC-RES  ⟹  Δ_cont ≥ liminf_{a→0} Δ_a .

  `SPEC-RES` here is the strongest (and most classical) form of the spectral
  convergence hypothesis, stated explicitly as a hypothesis:

      SPEC-RES-norm :  ‖H_a − H‖ → 0      (norm convergence of the rescaled
                                          Hamiltonians `H_a = −(1/a) log T_a`)

  Norm convergence ⇒ Hausdorff convergence of spectra ⇒ strong resolvent
  convergence, so `SPEC-RES-norm` implies the hypothesis used in REPORT.md §3.8.
  It is **not** the conclusion: it constrains the operators `H_a`, not `gap` or
  `liminf`.

  ## What is proved, what is assumed

  * `gap` is the **variational** gap (bottom of the Rayleigh quotient on the
    vacuum-orthogonal unit sphere).  Identifying it with `inf (σ H \ {0})`
    needs `ker H = ℝ Ω` (vacuum simplicity); that identification is standard,
    **not formalised here**, and *not needed*: the analytic inequality is what
    downstream consumes.  Declared gap, not a `sorry`.
  * `gap_lipschitz`: for bounded `H, K` with nonnegative form on `Ω^⊥`, the two
    variational bottoms differ by at most `‖H − K‖`.  The whole content is one
    observation — the *same* test vector is admissible for both operators — plus
    `|⟪ψ, Tψ⟫| ≤ ‖T‖` on the unit sphere.
  * `gap_tendsto` / `gap_ge_liminf`: `gap (H k) → gap H`, hence the `≥ liminf`
    form that GL2b' consumes.

  ## Declared Mathlib gaps

  * No `σ(H) ⊆ limsup σ(H_a)` for strong-resolvent convergence (Reed–Simon I,
    Thm VIII.24(b)), and no spectral perturbation theory for general
    self-adjoint operators.  We therefore formalise the **norm-convergence**
    route, where everything is elementary and complete.  The weaker
    strong-resolvent route is the paper remark in REPORT.md §3.8 and is *not*
    claimed here.
  * v4.33.0-rc1 has no `⟪_,_⟫_ℝ` notation in scope; `inner ℝ` is used.
  * No `sorry`; standard axioms only (`check.sh`).
-/
import Mathlib

namespace YMGlue

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- Unit vectors orthogonal to the vacuum: the variational domain of the gap. -/
def VacOrth (Ω : E) : Set E := {ψ | ‖ψ‖ = 1 ∧ inner ℝ ψ Ω = 0}

/-- Variational gap: bottom of the Rayleigh quotient on `Ω^⊥`. -/
noncomputable def gap (H : E →L[ℝ] E) (Ω : E) : ℝ :=
  sInf {r : ℝ | ∃ ψ ∈ VacOrth Ω, r = inner ℝ ψ (H ψ)}

/-- **GL2b' core.**  The variational gap is 1-Lipschitz in the operator norm. -/
theorem gap_lipschitz (H K : E →L[ℝ] E) (Ω : E)
    (hne : (VacOrth Ω).Nonempty)
    (hHnn : ∀ ψ ∈ VacOrth Ω, 0 ≤ inner ℝ ψ (H ψ))
    (hKnn : ∀ ψ ∈ VacOrth Ω, 0 ≤ inner ℝ ψ (K ψ)) :
    |gap H Ω - gap K Ω| ≤ ‖H - K‖ := by
  have hbdH : BddBelow {r : ℝ | ∃ ψ ∈ VacOrth Ω, r = inner ℝ ψ (H ψ)} :=
    ⟨0, by rintro r ⟨ψ, hψ, rfl⟩; exact hHnn ψ hψ⟩
  have hbdK : BddBelow {r : ℝ | ∃ ψ ∈ VacOrth Ω, r = inner ℝ ψ (K ψ)} :=
    ⟨0, by rintro r ⟨ψ, hψ, rfl⟩; exact hKnn ψ hψ⟩
  -- A single Rayleigh quotient moves by at most `‖A − B‖`.
  have hkey : ∀ (A B : E →L[ℝ] E) (ψ : E), ψ ∈ VacOrth Ω →
      inner ℝ ψ (A ψ) ≤ inner ℝ ψ (B ψ) + ‖A - B‖ := by
    intro A B ψ hψ
    have hsplit : inner ℝ ψ (A ψ) = inner ℝ ψ (B ψ) + inner ℝ ψ ((A - B) ψ) := by
      have hA : A ψ = B ψ + (A - B) ψ := by simp [sub_apply]
      rw [hA, inner_add_right]
    have hsmall : inner ℝ ψ ((A - B) ψ) ≤ ‖A - B‖ := by
      have h1' : |inner ℝ ψ ((A - B) ψ)| ≤ ‖ψ‖ * ‖(A - B) ψ‖ :=
        abs_real_inner_le_norm ψ ((A - B) ψ)
      have h2' : ‖(A - B) ψ‖ ≤ ‖A - B‖ * ‖ψ‖ := (A - B).le_opNorm ψ
      simp only [hψ.1, one_mul, mul_one] at h1' h2'
      have h1'' : inner ℝ ψ ((A - B) ψ) ≤ ‖(A - B) ψ‖ := (abs_le.mp h1').2
      linarith
    rw [hsplit]; linarith
  have h1 : gap H Ω ≤ gap K Ω + ‖H - K‖ := by
    have hle : gap H Ω - ‖H - K‖ ≤ gap K Ω := by
      obtain ⟨ψ₀, hψ₀⟩ := hne
      refine le_csInf ⟨inner ℝ ψ₀ (K ψ₀), ψ₀, hψ₀, rfl⟩ ?_
      rintro r ⟨ψ, hψ, rfl⟩
      have hmem : inner ℝ ψ (H ψ) ∈ {r : ℝ | ∃ ψ ∈ VacOrth Ω, r = inner ℝ ψ (H ψ)} :=
        ⟨ψ, hψ, rfl⟩
      have hc := csInf_le hbdH hmem
      have hk := hkey H K ψ hψ
      simp only [gap]
      linarith
    linarith
  have h2 : gap K Ω ≤ gap H Ω + ‖H - K‖ := by
    have hle : gap K Ω - ‖H - K‖ ≤ gap H Ω := by
      obtain ⟨ψ₀, hψ₀⟩ := hne
      refine le_csInf ⟨inner ℝ ψ₀ (H ψ₀), ψ₀, hψ₀, rfl⟩ ?_
      rintro r ⟨ψ, hψ, rfl⟩
      have hmem : inner ℝ ψ (K ψ) ∈ {r : ℝ | ∃ ψ ∈ VacOrth Ω, r = inner ℝ ψ (K ψ)} :=
        ⟨ψ, hψ, rfl⟩
      have hc := csInf_le hbdK hmem
      have hk := hkey K H ψ hψ
      have hnorm : ‖K - H‖ = ‖H - K‖ := norm_sub_rev K H
      simp only [gap]
      linarith
    linarith
  exact abs_sub_le_iff.mpr ⟨by linarith, by linarith⟩

/-- **GL2b'.**  Norm convergence of the rescaled Hamiltonians gives convergence
of the variational gaps. -/
theorem gap_tendsto (H : ℕ → E →L[ℝ] E) (Hlim : E →L[ℝ] E) (Ω : E)
    (hne : (VacOrth Ω).Nonempty)
    (hnn : ∀ k, ∀ ψ ∈ VacOrth Ω, 0 ≤ inner ℝ ψ (H k ψ))
    (hnnlim : ∀ ψ ∈ VacOrth Ω, 0 ≤ inner ℝ ψ (Hlim ψ))
    (hconv : Filter.Tendsto (fun k => ‖H k - Hlim‖) Filter.atTop (nhds 0)) :
    Filter.Tendsto (fun k => gap (H k) Ω) Filter.atTop (nhds (gap Hlim Ω)) := by
  have hbound : ∀ k, ‖gap (H k) Ω - gap Hlim Ω‖ ≤ ‖H k - Hlim‖ := by
    intro k
    have := gap_lipschitz (H k) Hlim Ω hne (hnn k) hnnlim
    simpa [Real.norm_eq_abs] using this
  have hzero : Filter.Tendsto (fun k => gap (H k) Ω - gap Hlim Ω) Filter.atTop (nhds 0) :=
    squeeze_zero_norm (f := fun k => gap (H k) Ω - gap Hlim Ω)
      (a := fun k => ‖H k - Hlim‖) hbound hconv
  have h := hzero.add_const (gap Hlim Ω)
  simpa using h

/-- The form consumed by the handoff: `Δ_cont ≥ liminf_{a→0} Δ_a`. -/
theorem gap_ge_liminf (H : ℕ → E →L[ℝ] E) (Hlim : E →L[ℝ] E) (Ω : E)
    (hne : (VacOrth Ω).Nonempty)
    (hnn : ∀ k, ∀ ψ ∈ VacOrth Ω, 0 ≤ inner ℝ ψ (H k ψ))
    (hnnlim : ∀ ψ ∈ VacOrth Ω, 0 ≤ inner ℝ ψ (Hlim ψ))
    (hconv : Filter.Tendsto (fun k => ‖H k - Hlim‖) Filter.atTop (nhds 0)) :
    gap Hlim Ω ≥ Filter.liminf (fun k => gap (H k) Ω) Filter.atTop := by
  have h := gap_tendsto H Hlim Ω hne hnn hnnlim hconv
  rw [h.liminf_eq]

end YMGlue

/-- info: 'YMGlue.gap_lipschitz' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMGlue.gap_lipschitz

/-- info: 'YMGlue.gap_tendsto' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMGlue.gap_tendsto

/-- info: 'YMGlue.gap_ge_liminf' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms YMGlue.gap_ge_liminf
