/-
# H1 sharpening — machine-checked

This file machine-checks the sharpened statements about condition `H1`
(the susceptibility / coarse-mode Poincaré ceiling `χ ≤ 1/λ`), complementing
`YMStatement/Reduction.lean`.

All four theorems below are **0 sorry** (ordered-field algebra only).

* `lemmaG_twoSided`  — Lemma G's two-sided squeeze `f/m ≤ χ ≤ M/m`.
* `cert_necessity`   — the certificate `χ ≤ 1/λ` forces a gap floor `λ·f ≤ m`
                       (given `f ≤ χ·m`).
* `cert_sufficiency` — a gap floor `λ·M ≤ m` together with `χ·m ≤ M` gives the
                       certificate `χ ≤ 1/λ`; i.e. the certificate is *not
                       weaker* than a gap floor.
* `H1_blind_to_thin_atom` / `H1_does_not_imply_atom`
                     — an explicit model with `H1` and `m > 0` but atom weight
                       `f = 0`: the machine-checked statement that `H1` alone
                       does NOT force a positive atom weight — "H1 is necessary
                       but not sufficient / blind to a thin atom".

Honest note: this machine-checks a *reduction* claim about `H1`, not the
Yang–Mills mass gap. No QED.
-/

import Mathlib

set_option linter.unusedVariables false

namespace YMStatement.H1Sharpening

variable {f m M χ lam : ℝ}

/-- Lemma G, two-sided squeeze: `f ≤ χ·m` and `χ·m ≤ M` (with `0 < m`, `0 < χ`)
give `f/m ≤ χ ≤ M/m`. -/
theorem lemmaG_twoSided (hm : 0 < m) (hχ : 0 < χ) (hlo : f ≤ χ * m)
    (hup : χ * m ≤ M) : f / m ≤ χ ∧ χ ≤ M / m := by
  constructor
  · rw [div_le_iff₀ hm]; exact hlo
  · rw [le_div_iff₀ hm]; exact hup

/-- Necessity: the certificate `χ ≤ 1/λ` (with `f ≤ χ·m`) forces a gap floor
`λ·f ≤ m`. -/
theorem cert_necessity (hf : 0 < f) (hχ : 0 < χ) (hlam : 0 < lam)
    (hχle : χ ≤ 1 / lam) (hlo : f ≤ χ * m) (hm : 0 < m) : lam * f ≤ m := by
  have hx : χ * lam ≤ 1 := by
    have h := mul_le_mul_of_nonneg_right hχle (le_of_lt hlam)
    rwa [one_div, inv_mul_cancel₀ (ne_of_gt hlam)] at h
  have h1 : lam * f ≤ lam * (χ * m) := mul_le_mul_of_nonneg_left hlo (le_of_lt hlam)
  have h2 : lam * (χ * m) ≤ m := by
    have e : lam * (χ * m) = (χ * lam) * m := by ring
    rw [e]
    calc (χ * lam) * m ≤ 1 * m := mul_le_mul_of_nonneg_right hx (le_of_lt hm)
      _ = m := one_mul m
  linarith [h1, h2]

/-- Sufficiency: a gap floor `λ·M ≤ m` together with `χ·m ≤ M` (and `0 < M`,
`0 < λ`) gives the certificate `χ ≤ 1/λ`; the certificate is not weaker than a
gap floor. -/
theorem cert_sufficiency (hM : 0 < M) (hχ : 0 < χ) (hlam : 0 < lam)
    (hup : χ * m ≤ M) (hgap : lam * M ≤ m) : χ ≤ 1 / lam := by
  have hm : 0 < m := lt_of_lt_of_le (mul_pos hlam hM) hgap
  have h1 : (χ * lam) * m ≤ lam * M := by
    have e : (χ * lam) * m = lam * (χ * m) := by ring
    rw [e]; exact mul_le_mul_of_nonneg_left hup (le_of_lt hlam)
  have h2 : (χ * lam) * m ≤ m := le_trans h1 hgap
  have h4 : (χ * lam) * m ≤ 1 * m := by rw [one_mul]; exact h2
  have h3 : χ * lam ≤ 1 := le_of_mul_le_mul_right h4 hm
  rw [le_div_iff₀ hlam]
  exact h3

/-- An explicit model satisfying `H1` and `m > 0` but with atom weight `f = 0`:
`(f, m, M, χ, λ) = (0, 1, 1, 1, 1)`. -/
theorem H1_blind_to_thin_atom :
    ∃ f m M χ lam : ℝ,
      0 < lam ∧ 0 < χ ∧ χ ≤ 1 / lam ∧ 0 < m ∧ f = 0 ∧ f ≤ χ * m ∧ χ * m ≤ M :=
  ⟨0, 1, 1, 1, 1, by norm_num⟩

/-- **H1 is necessary but not sufficient — blind to a thin atom.** There is no
way to deduce a positive atom weight `0 < f` from `H1` (`χ ≤ 1/λ`) together with
`m > 0`, `f ≤ χ·m`, `χ·m ≤ M`: the tuple `(0, 1, 1, 1, 1)` is a countermodel. -/
theorem H1_does_not_imply_atom :
    ¬ (∀ f m M χ lam : ℝ,
        0 < lam → 0 < χ → χ ≤ 1 / lam → 0 < m →
          f ≤ χ * m → χ * m ≤ M → 0 < f) := by
  intro h
  obtain ⟨f, m, M, χ, lam, hlam, hχ, hχle, hm, hf, hlo, hup⟩ := H1_blind_to_thin_atom
  have hpos : 0 < f := h f m M χ lam hlam hχ hχle hm hlo hup
  rw [hf] at hpos
  exact lt_irrefl 0 hpos

end YMStatement.H1Sharpening

#print axioms YMStatement.H1Sharpening.lemmaG_twoSided
#print axioms YMStatement.H1Sharpening.cert_necessity
#print axioms YMStatement.H1Sharpening.cert_sufficiency
#print axioms YMStatement.H1Sharpening.H1_blind_to_thin_atom
#print axioms YMStatement.H1Sharpening.H1_does_not_imply_atom
