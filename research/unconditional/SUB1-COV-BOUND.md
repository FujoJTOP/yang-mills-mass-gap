# SUB1-COV-BOUND — the SD identity + conditioning correction does not bound `Cov_μ(a_{k→0})`: exact but gap-blind, and `D ≤ c·W_geo` is false on the Savvidy cone

> **2026-11-23 · Target:** `‖Cov_μ(a_{k→0})‖ ≤ 1/λ` uniformly in `(L,a)`, the QV martingale having collapsed the backward to this one covariance. **Tested:** the exact (unconditioned) Schwinger–Dyson/Stein identity `E_μ[H_n]=E_μ[f²]`, plus a bound on the *conditioning correction* `D=Var_μ(f|a{=}0)−E_μ[⟨K⁻¹c,c⟩|0]` (the fibre Brascamp–Lieb defect), by either `c·W_geo` (`c<1`) or the exact Savvidy cone structure. **Verdict: does NOT reach 5/5 — 2 passes, 3 fails, and the three failures are one object.** The SD identity is exact and acid-clean but *unconditioned*: at the reflection-symmetric vacuum `E f=0` it reads `Var(f)=E[H]` (equipartition), so it says nothing about the *marginal* variance `Var_μ(a_{k→0})`. The whole gap is the conditioning step, and `D` cannot be bounded without assuming it.

## The identity that nearly closes — and why conditioning is the whole gap
Flat chart `X=log U`, `μ∝e^{−V}dX` (`V=S_W+Φ`), mode `a=⟨v,X⟩`, force `f=∂_nV`. Integration by parts gives exactly, for **every** measure, any `β`, no convexity/positivity, `E_μ[H_n]=E_μ[f²]`, hence `Var_μ(f) ≤ E_μ[H_n]` (the measured `E[H]−Var≈0` at β=2.5 *is* this identity at `E f=0`). But the target is **conditioned**: `G″(0)=E_μ[H_n|a{=}0]−Var_μ(f|a{=}0)+W`. Conditioning on the coarse slice `a=0` **replaces the slack `(E f)²` by the coarse-slice variance** — and that slice term *is* the entire `W`+coarse content.

**Acid test (decisive).** Free massless `k=0`: `H_n=A_vv=0`, so `E[H]=E[f²]=0` — finite — while `Var(a)=1/(A^{-1})_{vv}=∞`. The identity is satisfied while the target diverges: it is **gap-blind**, giving `G″=0`, `f₁²=0` vacuously (`FORGE5-CONE-VARIANCE`, `k0-CONSTRAINT-EFFPOT`).

## Can `D` be bounded by `c·W_geo` (`c<1`)? No — fails on the cone, wrong sign
`D = Var(f|0) − E[K⁻¹c,c]` is controlled by `Var(f|0) ⊇ Var(G|0)`, the conditioned *dynamical* force variance. On the Savvidy cone `{2g|F|>μ²}` the compactness–sensitivity (Efron–Stein) bound gives `Var(G|0) ≤ diam²·E‖Hess V·v‖² ≥ (g|F|)²` — **largest exactly on the cone** — while `E[H|0] ≥ −Cg|F|` (Savvidy floor), `W_geo ≤ N/12`. Hence `G″ ≥ E[H] − diam²‖HessV·v‖² + W ≲ −diam²(g|F|)² < 0`: a *magnitude* bound of the wrong sign. The shortfall `~g|F|/(N/12)`, with `g|B|=√(2/β)|B|`, is ≈3.6, 7.2, 10.7 at `|B|=1,2,3`, β=2.5 — the recorded "≳10×" (`PROOF-SAVVIDY-EXACT-SIGN` iv). So `D ≤ c·W_geo` is uniformly **false**: `W_geo ⪰ N/12` is too small by ~10× and the wrong object (Jacobian curvature, not dynamical force). `Var(G|0)` is the induced-mass/H1 quantity (`PROOF-JACOBIAN-SAVES` (◆)).

## Can `D` be bounded by the exact Savvidy cone structure, gap-free? No — inverse-Hessian singularity = H1
The cone fixes the *forward* negative spectrum: count `N_neg/V=c_N(gB)²/8π² ∝ ∫|F|²`, depth `λ_min=−g|B|`. But `D` depends on `K⁻¹` (the *inverse* fibre Hessian). `K` is singular when a negative eigenvalue crosses zero (a transverse massless mode, `det K=0`), where `‖K⁻¹‖=1/|λ_min|→∞`. Count and depth control `K` *deep* in the cone (`K⁻¹` small) but say nothing about `K⁻¹` **near the edge**, where the measure still has weight; the uniform bound on `E‖K⁻¹‖` that would make `D` finite *is* the fibre spectral gap. The exact-RHS form of any such bound is `Var(G|0) ≤ E[H|0]+W−δ ⟺ G″≥δ ⟺ m_eff²(k→0)≥δ` (`FORGE5-CONE-VARIANCE` #3) — literally H1.

## Honest five-point scorecard
1. **Nonperturbative — PASS.** IBP identity exact, all `β`, no expansion.
2. **Sign-controlling — FAIL.** Needs `D<S`; unproven, and the cone bound is negative-signed.
3. **Non-circular — FAIL at closing.** Bounding `D` = bounding `Var(G|0)` = H1.
4. **Non-SUSY / free field — vacuous.** Returns `G″=0` for free massless — but only because it is gap-blind there; no positivity extracted for YM either.
5. **Uniform `≥λ` — FAIL.** `K⁻¹` singular at the cone edge; no `λ` independent of `(L,a)`.

## Exact failing step
Bounding `D` **from above** requires the conditioned dynamical-force variance `Var(G|0)` (equivalently `E_μ[K⁻¹c,c]` with `K` uniformly invertible). The SD identity is *forward* and cannot supply it — the free `k=0` mode saturates it with `Var(a)=∞`. On the cone the only gap-free bound on `Var(G|0)` is a *magnitude* bound largest where the sign must flip, missing `W_geo` by ~10×; the cone-structure route needs `K⁻¹` near the zero-crossing, i.e. the fibre gap. This is **H1 renamed**, not evaded. No 5/5.

*Refs.* FORGE5-CONE-VARIANCE (#2/#3/#5); INVENT-VARIANCE-INEQ (`G″=S−D`); PROOF-JACOBIAN-SAVES (◆); PROOF-SAVVIDY-EXACT-SIGN (count ∝∫|F|², depth `−g|B|`); FORGE5-SL-MARTINGALE (BT′); k0-CONSTRAINT-EFFPOT; PROOF-COARSE-POINCARE (H1); PROOF-GAUSS-INVARIANT-PF; Efron–Stein 1981; Stein/Schwinger–Dyson; Brascamp–Lieb 1976; Prékopa 1971.
