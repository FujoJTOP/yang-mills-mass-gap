# PROOF-EDGE-E2 — lower-bounding the second level `Λ=E₂`: does the block-local gap break the Temple regress?

> 2026-11-23 · The residue–level-gap bound `f₁² ≥ (ΛM₀−M₁)/(Λ−m₁²)` is empty without a lower bound on the second threshold `Λ=E₂`. Three routes attacked — Temple/Kato–Weinstein, OPE, two-level GEVP. **Verdict: not closed at any coupling; the failure is a theorem about the *tools*: every available spectral instrument is a bottom-Dirichlet functional, blind to the second gap `λ₁−λ₂`; the proven block-local gap does not break the regress.** The datum is exactly the "upper gap" `I2` of `mass-spectrum` (there open); the exact missing estimate is a two-mass bound. No fabricated proof.

## Target, reframed
`E₂≥E₁>0` is free (eigenvalue ordering + the proven first gap) but *empty* here: `M₁/M₀>m₁²`, so `Λ=E₂=E₁` gives `ΛM₀−M₁≤0`. Nonvacuity needs isolation, `Λ>~M₁/M₀` (equivalently `∫(s−Λ)dρ < f₁²(Λ−m₁²)`), i.e. `I2`: `spec(H)|_{0⁺⁺} ∩ (m₁,(1+δ)m₁)=∅`. By **Theorem A** of `mass-spectrum` (proved) `atom(I1) ∧ I2 ⟺ sharp two-exponential decay Ze^{−m₁t}+O(e^{−m₂t})`. `I2` is open.

## (i) Temple/Weinstein — the regress is structural
Temple `E₀ ≥ R − Var_ψ(H)/(ρ−R)` is valid **only** for `ρ ≤ E₁` (a lower bound on the *second* level; proof `Var ≥ (ρ−R)(R−E₀)` given `E₀≤ρ≤E₁, R<ρ`, equality iff a 2-level trial with `ρ=E₁`). To bound `E₂` apply it to `H|_{(vac⊕glueball)^⊥}`, whose bottom is `E₂` and second is `E₃` ⇒ needs `ρ≤E₃`: one rung up (LEVER §1). The **block-local gap** cannot supply it: it is a Poincaré inequality, and a Poincaré/spectral-gap datum encodes **only** `μ₁=−log λ₁=a m₁`; `λ₂` never enters. The spectra `{1,λ₁,λ₁,…}` and `{1,λ₁,λ₂,…}` share it and differ in `λ₁−λ₂`; no bottom functional sees `λ₁−λ₂`, which is the bottom of `H` on the *first-eigenspace complement* — not a quadratic-form quantity. (At β=2.5 it does not even globalize: `C≈15`, PROOF-H2-GLOBAL.) **Dead, precisely.**

## (ii) OPE — the approach is *from above*
For `ν≥0`, `m_eff(t)=−∂_t log C(t)` is a mean of `√s` under the `e^{−√s t}`-reweight, hence `≥m₁ ∀t`, `↑m₁`. Its first correction is governed by `Λ−m₁`, so lower-bounding the rate — the only way to lower-bound `Λ` — is the same second-level datum (moment-ratio lemma, GAP-G2-edge). Individual glueball coefficients are unprotected and the sign Schur-concave (EDGE-f1²; INEQ-T2: the residue-dual is an infimum). *Symmetry is blind too*: in `0⁺⁺`, C and P cannot separate `gg` from `(gg)(gg)` (both C-even, P-even), so no Z₂ protects the upper gap — unlike O(3)'s spin-flip S2. **Dead.**

## (iii) GEVP `R_opt` — variational ⇒ upper only
Courant–Fischer: every GEVP eigenvalue is `λ_k^{(n)}≥E_k`, an **upper** bound; `R_opt` upper-bounds `m₁²`. A difference of upper bounds bounds neither `E₂` nor `Λ−m₁`. No rigorous second-level lower bound is extractable; at β=2.5 `R_opt` is cutoff-dominated anyway (ATOM). **Dead.**

## The exact missing estimate
A two-mass asymptotic in the `tr F²` channel,
`C(t)=A e^{−m₁t}(1+η(t)), |η(t)|≤B e^{−δt}, A≥a₀>0, δ>0,`
whose proof gives `supp ν ∩ (m₁²,(m₁+δ)²)=∅`, i.e. `Λ≥(m₁+δ)²` — exactly `I2`+`I1`, hence (Theorem A) sharp two-exponential decay. It dominates the one-gap cluster/Poincaré theorem (it pins the subleading exponent) and, by the residue–level-gap lemma, is *equivalent* to `f₁²≥a₀`. Strong coupling supplies only `m₁>0`; the two-mass bound is no known theorem there, so `E₂−E₁>0` relocates to the same `H*`.

**Sharpest line.** The second-level bound is missing not for a weak constant but because Poincaré/SZZ, the block-local gap, cluster one-gap, OPE-approach and GEVP are all *bottom* functionals — identical under multiplicities of the lightest level, blind to `λ₁−λ₂`, and even C/P cannot separate `gg` from `(gg)(gg)` in `0⁺⁺`. Residue–level-gap converts `λ₁−λ₂` into `f₁²` and back: closing `I2` *is* `f₁²>0` *is* `H*`; the concrete open theorem is the two-mass isolation bound.

*Refs.* GAP-G2-edge (residue–level-gap lemma); RESIDUE-CONDITIONAL (`H*`); LEVER-nonconvex-edge §1; INEQ-T2-variational-residue (E₂-collapse); EDGE-f1sq-attack (sign); PROVE-H2-uvgap (block gap); PROOF-H2-GLOBAL (`C≈15`); ATOM-trial-state; mass-spectrum/STATUS (Theorem A, `I1`/`I2`, `O2`); Temple 1928; Courant–Fischer; SZZ arXiv:2204.12737.
