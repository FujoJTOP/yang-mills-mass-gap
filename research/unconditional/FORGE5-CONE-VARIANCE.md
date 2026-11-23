# FORGE5-CONE-VARIANCE — the Efron–Stein / compactness cone-variance bound: nonperturbative, but dead on sign

> 2026-11-23 · FORGE5 target: a *genuinely new* upper bound `Var(G|0) ≤ E[H]+W−δ` on the Savvidy cone `{2g|F|>μ²}`, not positivity, not convexity, acid-test-clean. **VERDICT: cannot satisfy all 5.** I construct the **compactness–sensitivity bound** `Var(G|0) ≤ diam²·E‖Hess V·v‖²` (Efron–Stein + SU(N) compactness) — a sixth species (sensitivity, not BL/Prékopa). It **PASSES #1** (nonperturbative, all β), **FAILS #2** (it is a *magnitude* bound, `≥(g|F|)²` on the cone, so its lower bound on `G″` is *negative* — the wrong sign), **FAILS #5** (no `λ` follows), and in its exact-target form is literally **#3** (H1, the gap). **#4 is nominally met but vacuous** (the noncompact free field escapes via unbounded fiber fluctuation, extracting no positivity — but neither does YM).

## The exact identity that nearly closes it — and why conditioning is the whole gap
Flat chart, `μ∝e^{−V}dX` (`V=S_W+Φ`), mode `a=⟨v,log U⟩`, force `f=∂_nV`. Integration by parts (`∂_ne^{−V}=−fe^{−V}`), taken with `g=f`, gives the Schwinger–Dyson / Stein identity
`E_μ[H_n]=E_μ[f²]`,  hence  `Var_μ(f)=E_μ[H_n]−(E_μf)² ≤ E_μ[H_n]`.
This is a genuine variance upper bound of the requested *shape*, valid for **any** measure, any β, **no convexity, no positivity** (numerically exact on an anharmonic 1D measure: `E[H]=3.834`, `E[f²]=3.833`). But it is the **unconditioned** variance. The gap lives in `G″(0)=E_μ[H_n|a=0]−Var_μ(f|a=0)+W` — the **conditioned** object; conditioning on `a=0` replaces the slack `(E_μf)²` by a coarse-slice variance, and *that slice term is the entire `W`+coarse content*. At the reflection-symmetric vacuum `E f=0`, so unconditioned `Var=E[H]` **exactly** (equipartition — the repeatedly measured `E[H]−Var≈0`); the gap is 100% a conditioning effect, invisible here. The Gaussian free field saturates it: `G″=E[H]−Var+W=0`, massless mode `δ=0`.

## The constructed inequality (boxed)
Efron–Stein, with compactness `|δX_i|²≤diam²` on each SU(N) link fiber:
```
   SCV (Savvidy cone-variance, compactness–sensitivity):
   Var_μ(∂_nV | a=0) ≤ diam_SU(N)² · E_μ[ ‖Hess V·v‖² | a=0 ],   all β,
```
nonperturbative, non-convexity, non-positivity.

## Exact failure
- **(#2 SIGN — fails.)** On the cone `‖HessV·v‖²≥(g|F|)²≥(μ²/2)²`, while `E[H|0]≥−Cg|F|` (Savvidy floor) and `W≤N/12`. Substituting, `G″≥E[H]−diam²‖HessV·v‖²+W ≲ −diam²(g|F|)² <0`: a **magnitude** bound, largest *exactly* on the cone, yielding a **negative** lower bound — no `δ`. It is off by a factor `~g|F|`, the recorded "Cauchy–Schwarz `Var≤E|G|²∼β²` misses `E[H]∼β` by β" (PROOF-MULTISCALE-4D §2).
- **(#3 CIRCULAR — fails as sharp form.)** With the *exact* RHS, `Var(G|0)≤E[H|0]+W−δ ⟺ G″≥δ ⟺ m_eff²(k→0)≥δ` (PROVE-Gpp §5). Proving it **is** the gap. Escaping #3 via the loose §2 bound is precisely what forces #2.
- **(#5 UNIFORM — fails.)** No `λ` follows, since #2 gives no `≥δ`.
- **(#1 PASS; #4 nominal.)** Nonperturbative ✓, non-SUSY ✓; the noncompact free field violates compactness (unbounded fiber fluctuation) so no positivity is extracted — but none is extracted for YM either → the acid-test pass is vacuous.

## Sharpest line
`E[H]=E[f²]` is a theorem for *every* measure — YM and free alike, unconditioned. The conditioned `Var(G|0)` that signs `G″` differs from it by exactly the coarse-slice (W) term; every non-positivity, non-convexity magnitude bound (Efron–Stein, `Var≤E|G|²∼β²`) is **largest on the Savvidy cone, the wrong sign**. To get `−δ` with the right sign requires `W>Var−E[H]` on the cone — which is H1. **Fails #2, #3(sharp), #5; #4 vacuous; #1 only.**

*Refs.* PROVE-Gpp-positive §5 (`G″≈1/Var`); PROOF-MULTISCALE-4D §2 (`Var≤E|G|²` misses by β_j); PROOF-COARSE-POINCARE (H1); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`); INVENT-VARIANCE-INEQ (`G″=S−D`, `D>0`); PROOF-JACOBIAN-SAVES (`W_geo≥N/12`); Efron–Stein (1981); Stein's lemma / Schwinger–Dyson identity.
