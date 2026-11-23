# k0-TRANSLATION-RIGIDITY — does translation-rigidity make the zero-momentum mode stiff?

> 2026-11-23 · Forge of the "translation-rigidity" species for the `k→0` coarse mode `a_0`, `Var_μ(a_0)≤1/λ` (`⟺ V_eff″(a_0)≥λ ⟺ m_eff²(0)≥λ` = H1). **Verdict (判词): NO — a no-go, and a sharp one.** The intended lever is self-annihilating: `a_0` is the zero-momentum (constant) mode, hence **translation-invariant** — every translation generator `v_μ=∂_μ·` has *exactly zero* overlap with `a_0`. A "Poincaré in the translation direction" therefore reads `Var(a_0)≤C·E|∇_trans a_0|²=C·0=0` — false. Translation is the one direction where the restoring curvature on `a_0` is *provably zero*: the wrong sign, not a lower bound. The other reading ("free energy convex in `a_0`") is Cramér/Legendre convexity — **automatic and content-free**, since the rate-function curvature at its minimizer *is* `1/Var(a_0)`. So "strict convexity with uniform `λ`" is literally the target restated; nothing sits between "convex" and "convex with curvature `λ`." The acid test kills it cold: translation invariance is shared by the **massless** free field, so a translation-rigidity claim would wrongly gap it (`Var=∞` there) — vacuous.

## 1. The object is fixed by translation — the core failure
`a_0=V^{-1}∫ℓ(U)d⁴x`. Global translation `(T_tA)(x)=A(x+t)` gives `a_0(T_tA)=a_0(A)` **for all `t`**; the generator `v_μ=∂_μU` satisfies `⟨∂a_0,v_μ⟩=V^{-1}∫∂_μ(…)dx=0` (total derivative on the torus) — indeed `|∇_trans a_0|²≡0`, *exactly*. The translation Dirichlet form is the **zero form** on `a_0`; a Poincaré along it would force `Var(a_0)=0`. The stiffness lives in the **constant-shift** direction `A→A+c` (a symmetry of the *massless free* field, absent in YM), **not** in translation. Same anatomy as PROOF-IR-H1 (`∇_⊥` annihilates coarse functions): the named form is kernel on the target mode.

## 2. Free-energy convexity is automatic — circular
`Ψ(h)=log E_μ e^{h a_0}` is convex (Hölder); the rate function `I(x)=sup_h(hx−Ψ)` is *always* convex (Legendre), and at its minimizer `I″(x*)=1/Ψ″(0)=1/Var(a_0)`. Thus "`V_eff″≥λ`" **is** "`Var(a_0)≤1/λ`." Convexity `I″≥0` is free and uniform in every theory; the strict, uniform `λ>0` is exactly the gap. Massless free field: `I″(min)=0` (flat); massive: `I″=m²`. Free-energy rigidity does not *supply* `λ` — it *is* `λ`.

## 3. Acid test — vacuous on the free field
Translation invariance holds at `m=0` and `m>0` alike. "Translation invariance ⇒ `a_0` stiff" predicts `G″>0` also at `m=0` (where `Var=∞`) — plainly false. If rephrased as "rigid under volume-averaging ⇒ `Var` small," that needs decay of connected correlations — the gap itself (PROOF-CONCENTRATION: concentration is downstream).

## 4. Five-point verdict
1. **Nonperturbative** — yes in form (a global measure statement), but this buys nothing.
2. **Sign-controlling** — **NO.** `|∇_trans a_0|²≡0`; the form certifies zero dissipation (would force `Var=0`), the opposite of a positive lower bound.
3. **Non-circular** — **NO.** Strict convexity with curvature `λ` ⟺ `Var(a_0)≤1/λ` = H1.
4. **Fails free field** — **NO.** Translation invariance is a *massless*-free-field property; the mechanism would gap it, hence vacuous. (Not SUSY-related; it does not exploit non-SUSY dynamics.)
5. **Uniform lower bound** — **NO.** Cramér gives `λ` possibly `0`; Brascamp–Lieb would need an input Hessian `⪰λ` = H1.

## 5. Exact failure / what is true
`a_0` *can* be stiff precisely because **no symmetry fixes it** — and translation *does* fix it, so translation cannot move it, let alone lower-bound its variance. The honest content: `a_0` is unprotected (no shift symmetry in YM, unlike the free field), but its stiffness is dynamical — `V_eff″=Δ_dyn+W_geo+X` (PROOF-JACOBIAN-SAVES), whose positivity is `(◆)`=H1=the gap. **Missing step, unchanged:** a non-symmetry, dynamical `V_eff″(a_0)≥λ` uniform in `(V,a,boundary)` — `(CP)`/H1. (Distinct from the *index*-rigidity no-go, INVENT-NEW-INDEX; here it is *thermodynamic/translation* rigidity: different object, same conclusion.)

*Refs.* PROOF-COARSE-POINCARE (H1 at k→0); PROOF-IR-H1 (orthogonal split, `∇_⊥` annihilates coarse fns); PROOF-DIRICHLET-COMPARISON; PROOF-JACOBIAN-SAVES (`G″≈1/Var`, `G″=Δ_dyn+W_geo+X`); PROOF-CONCENTRATION (concentration downstream); PROOF-THERMODYNAMIC-FORMALISM (pressure sees `λ₀`, not `λ₁`); INVENT-NEW-INDEX (index-rigidity no-go, distinguished); Brascamp–Lieb (1976); Cramér/Varadhan; Holley–Stroock (1987); Shen–Zhu–Zhu, CMP **400** (2023).
