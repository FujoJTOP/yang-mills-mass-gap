# SCORE-HYPOTHESIS — is `λ > 2 sup Var_μ(score)` weaker than the mass gap?

> 2026-11-23 · The decisive question left by PROVE-H3(b). **Verdict: EQUIVALENT — in fact ≥ gap. The repair is sound (H3 was false), but it is not a partial reduction.** `Var_μ(score)` itself is a genuine *new* object and plausibly a *provable UV* constant; the absorption hypothesis is nonetheless gap-strength, because the `λ` on the other side of the inequality is H1's *global* coarse Poincaré/SLC constant, whose uniform infimum is the `k→0` transverse mode. No fabricated proof.

## The object, and which `λ`
H3(b)'s contraction fails; the chain-rule repair (PROVE-H3 line 10) leaves `∇φ = E_μ[∇_∥F|U'] + Cov_μ(F,s)`, `s = ∂_{U'}\log μ(·|U')`. Absorption gives `Var_μ(F)(1 − 2V_s/λ) ≤ …`, needing
`λ > 2V_s`,  `V_s := sup_{|v|=1} Var_μ(⟨s,v⟩)`.
Here `λ` is *not* free: it is the constant in Step 1, `Var_ν(φ) ≤ (1/λ)E_ν|∇φ|²`, i.e. H1's coarse curvature `Hess S_eff⪰λ` — a bound required for the **whole** coarse marginal, whose uniform infimum is taken at `k→0`.

## Side (a): `Var(score)` is UV — and likely provable
Per coarse link `U'_{(B)}`, the score touches only the `O(1)` block-boundary plaquettes (factor-2 block), each with `|∂_{U'}Re tr U_p/N|≤1`, uniformly in `ξ,U',L`. The fiber `ξ=ker P` is *band-limited* (`|k|≥π/2a`), so its correlations — hence `Cov(s(B),s(B'))` — have range `O(a)`. Then `V_s ≤ C_UV(β,b,N)·(Σ_{B'}|Cov|) < ∞`, uniform in `V,a,U'`. This is a UV Fisher-information bound, the same *caliber* as H2-b-uniform (it needs the fiber to be a massive short-range shell, ~`π/2a`), **not** the IR gap.

## Side (b): `λ` is IR — so the hypothesis is the gap
Grant even a *perfect* bound `V_s ≤ C_UV`. The absorption becomes `λ > 2C_UV`: the **global** coarse Poincaré constant must exceed a *known explicit* number. Via the soft mode (PROVE-Gpp §5: `H1 ⟺ Var(soft)≤1/λ`; W6) this says
`m_eff²(0) ≥ λ > 2C_UV`,
a **quantitative** mass gap. One cannot dodge to the block scale: a block-local `λ` needs localization back to the global fiber, and that is exactly H2-globalization, where the Zegarlinski/Nachtergaele contraction fails at `β=2.5` (open). So the global `λ` — the gap — is unavoidable.

**Strength ordering.** `λ>2C_UV` implies `λ>0` (bare H1 ⇒ gap), but the bare gap (`m_eff²(0)>0`) does **not** imply `m_eff²(0)≥2C_UV` (a gap permits arbitrarily small mass). Hence the repaired hypothesis is **strictly stronger** than the gap — certainly not weaker.

## Relation to the equivalence circle
`{m_eff²(0), G″, f₁², γ_LGP, Polchinski c}` are **one IR object**: the `k→0` transverse mode in five vocabularies (CONVERGENCE-redteam §c). `Var(score)` is **not** in this circle — it is a UV, finite-block Fisher information, a different *kind* of quantity. So it is a "new member" only in the trivial (computable) sense; as a hypothesis it is IR-hard because it is paired with `λ∈` the circle. The score supplies no weaker coordinate.

## Bottom line for the conditional theorem
The repair (H3(b) → `λ>2V_s`) makes the theorem *correct* (the old H3 was Jensen-only, false for a `U'`-dependent fiber) and cleanly **localizes all irreducibility in `λ`**: `Var(score)` is UV and (plausibly) provable, so nothing IR hides in the mixed-Hessian term itself. But `gap ⟸ (H1 with λ>2C_UV)` assumes gap-strength, so it is **not a partial reduction** — merely a sound sufficient condition. Honest deliverable: a UV lemma worth proving (`V_s ≤ C_UV`, ~H2-b caliber) that *confirms* the gap sits entirely in the coarse curvature `λ`, and **no** weakening of it.

*Refs.* PROVE-H3 (b); CONDITIONAL-LATTICE-GAP H1/Step1; PROVE-Gpp-positive §5; W6-induced-mass; W6-BE-at-2.5; H2-b-uniform; H2-globalization; VERIFY-LOAD-BEARING (ii); CONVERGENCE-redteam (c).
