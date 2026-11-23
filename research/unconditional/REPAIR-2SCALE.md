# REPAIR-2SCALE — fixing H3(b): absorb the score on the fiber, not against the total variance

**Verdict.** A valid repair exists (options **(a)+(c)**, stated as **(d)**). The false H3(b) — the contraction `E_ν|∇φ|² ≤ E_μ|∇_∥F|²` — is replaced by the *exact* chain rule plus a Stein / Fokker–Planck absorption of the score on the fiber. H3 collapses to a structural identity plus a **UV** condition (`S<∞`, automatic by compactness on `G^{n_B}`); it does **not** need the gap-strength `λ>2 sup Var(score)` of SCORE-HYPOTHESIS. All IR content stays in **H1 = the gap in curvature language**. The theorem is valid again as a conditional theorem; its one open obligation is the pre-existing global-fiber Poincaré (H2-globalization), untouched by the score. **(b) BL-on-the-true-split does not help on its own** — the coarse marginal is exactly where the gap lives, so BL there re-assumes H1; the non-orthogonality only adds harmless `O(1)` metric factors.

## Exact chain rule (route a)
`φ(U′)=E_μ[F|U′]`, `s_B=∂_{U′_B}\log μ(·|U′)=−(∂_{U′_B}S_W−E[·])` (centred coarse force; `Var_μ(s_B)=Var(∂_{U′_B}S_W|U′)`). Then
`∂_{U′_B}φ = E_μ[∂_{U′_B}F|U′] + Cov_μ(F,s_B).`
The first term contracts by conditional Jensen (`E_ν|E[X|U′]|²≤E_μ|X|²`). Only `Cov` is the gap.

## Stein / fiber-generator absorption (route c)
Fix `U′`; let `L_⊥` be the fiber Langevin generator (invariant `μ(·|U′)`, Dirichlet form `E_⊥(f)=E|∇_⊥f|²`, **gap `Λ_UV`** by H2). Since `s_B` is mean-zero and smooth, solve `L_⊥ g_B=s_B` on the mean-zero fiber. Reversibility gives
`Cov_μ(F,s_B) = E_μ[F L_⊥g_B] = −E_μ⟨∇_⊥F,∇_⊥g_B⟩`.
Cauchy–Schwarz, and `‖∇_⊥g_B‖²=E_⊥(g_B)=−⟨g_B,s_B⟩≤(1/Λ_UV)Var_μ(s_B)`:
```
E_ν|Cov(F,s_B|U′)|² ≤ (S/Λ_UV)·E_μ|∇_⊥F|² ,   S := sup_{U′,B}Var_μ(s_B|U′) ≤ 4‖∂_{U′}S_W‖_∞² < ∞.
```
This routes the score into the **fiber Dirichlet budget** using the (proved, block-local) fiber gap `Λ_UV` — *not* into `Var_μ(F)`. No bootstrap on `Var_μ(F)`, hence **no `λ>2S`**. *(Numerically verified on the PROVE-H3 toy `½a²+½b²+ga²b²`, `F=b²`: `Cov=−0.891≈∂_aφ`, and `‖∇g‖²=0.297≤Var(s)/Λ_UV=0.595`.)*

## Repaired theorem (d), exact
Assume **(H1)** `Hess S_eff⪰λ>0` (IR); **(H2)** fiber gap `Λ_UV>0` (UV — its *global* form is the open H2-globalization); **(H3′)** exact `(U′,ξ)`-measurability + `S<∞` (structural/UV, automatic). Then with
```
1/m = max( 2/λ , (1/Λ_UV)(1+2S/λ) ),   i.e.  m = min( λ/2 , Λ_UV/(1+2S/λ) ) > 0 ,
```
`Var_μ(F) ≤ (1/m)Σ_x E_μ|∇_xF|²`, hence `Δ_Λ≥m`. Proof: `Var_μ(F)=E_νVar_fib(F)+Var_ν(φ) ≤ (1/Λ_UV)E|∇_⊥F|² + (1/λ)(2E|∇_∥F|² + 2(S/Λ_UV)E|∇_⊥F|²)`. Under H2 stated globally, `L_⊥^{-1}` is global and this is unconditional; under block-local H2 alone the score shares the same open localization the fiber term already needs (UV, per H2-b-uniform).

## Strength vs the bare gap
- **H3′ is strictly weaker than the gap**: a UV compactness bound, satisfied by free/massless theory (finite force variance). It injects no IR.
- **The total IR burden is H1 alone — which *is* the bare gap** (`G″⪰λ ⟺ m_eff²≥λ`, PROVE-Gpp §5). So the repaired theorem is **equal to the bare gap in IR content**: neither weaker nor stronger. The repair buys *correctness*, not leverage.
- **Versus the crude `λ>2S`** (SCORE-HYPOTHESIS): that closes a bootstrap on `Var_μ(F)` and is **stronger** than the bare gap (a quantitative floor `>2S`). The generator route **avoids it** — `λ` stays the global H1 constant (so SCORE-HYPOTHESIS's "cannot dodge λ to the block scale" does not bite); the price is a factor `(1+2S/λ)` on the `∇_⊥` coefficient, large precisely when `λ` is small (the physical soft regime, `λ≈G″≈0.27`, `S≈Var(force)≈2.4` ⇒ factor `≈19`).

**Net.** CONDITIONAL-LATTICE-GAP is valid again as a conditional theorem; its obligations shrink `{H1,H2,H3}→{H1=gap (IR), H2-global (UV, open)}`. No score-generated IR; no contraction needed; no leverage on the Clay gap.

*Refs.* PROVE-H3(b); SCORE-HYPOTHESIS (crude route, corrected here); H2-TIEBREAKER; H2-globalization; PROVE-H2 (block-local `Λ_UV`); PROVE-Gpp §5; W7; Bakry–Émery/Stein; Holley–Stroock.
