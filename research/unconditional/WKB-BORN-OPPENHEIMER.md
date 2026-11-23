# WKB-BORN-OPPENHEIMER — BO with a rigorous fast-sublevel lower bound: the indefinite Savvidy Hessian is *handleable for BO's validity, not for the gap*; the slow stiffness `Hess E_0=G″` is H1, and the cone puts the side BO needs gapped exactly on the instability

> **判词 (verdict).** BO splits the coarse mode `R=a₀` (slow) from the Savvidy fiber `ξ` (fast). Plain BO renames H1 (PROOF-BORN-OPPENHEIMER); the *new* question is whether a rigorous **WKB/adiabatic lower bound on the fast sublevels** — a non-Gaussian min–max bound surviving an indefinite Hessian — becomes a sign-controlling bound on `G″(k→0)`. **It separates two objects, then cannot manufacture the gap.** The fast gap `ε_1(R)` is positive on the compact fiber and certifiable by min–max even where `H_ff` is indefinite (a genuine new bound, `ε_1=Λ_UV`), but `ε_1` enters only the adiabatic correction `A(R)=Σ|A_{n0}|²/(ε_n−ε_0)` and the `1/ε_1` memory scale, never the static stiffness `Hess E_0=G″`. So `Hess H_slow=G″+Hess A`; the gap is the leading `G″=Δ_dyn+W_geo+X` = **H1**, and WKB bounds `A`'s size, not its sign (`X` sign-free, C8). **BO's premise fails where it matters:** negatives are 90–97% *fiber* (W6) — the side BO needs gapped carries the instability, `ε_1` only marginal (`~e^{−2βN_p}`) with a Landau–Zener crossing at the cone. **The indefinite Hessian breaks the Gaussian/Schur potential and the WKB fast gap there; it does not break the exact compact-fiber integral `E_0`.**

## §1 · The two objects BO separates
`μ∝e^{−S_W(R,ξ)}J\,dR\,dξ`, exact fiber free energy `E_0(R)=−log∫dξ\,e^{−S_W}J`. Klein–Martinez–Seiler/Teufel: `H_slow=−(1/2M)Δ_R+E_0(R)+A(R)`, `A(R)≥0`, error `O(1/ε_1)`; slow gap `=Hess E_0(R*)+Hess A(R*)`. The k0/W7 identity fixes the leading term: `Hess E_0(R*)=G″=Δ_dyn+W_geo+X`, the static fiber-free-energy curvature, **independent of `ε_1`**. Fast gap and slow stiffness are distinct; the WKB bound is on the first.

## §2 · What the WKB fast-sublevel bound does and does not buy
Fast generator = the fiber part of E1-FRG's exact FP operator at frozen `R`; spectrum `0=ε_0<ε_1(R)≤…`. Fiber compact, `e^{−S_W}J>0` ⇒ `ε_1(R)>0`; **min–max / Birman–Schwinger** certifies `ε_1≥γ` with no Gaussian — valid even where `H_ff=−D²+2adF` is indefinite (`λ_min=−g|B|`). This is the genuine new object. But `ε_1=Λ_UV` is a **relaxation rate** bounding `|A|`'s smallness, not `Hess E_0`. A lower bound on fast *excitation values* is not a lower bound on the slow *curvature*; strong convexity `V_eff(R)−V_eff(R*)≥(λ/2)|R−R*|²` — the only shape giving `G″≥λ` — is H1.

**Berry–Landau sub-idea (creative; it dies).** `A`'s Berry curvature `F` could gap a flat `E_0` magnetically. But the goal mode is one direction (`F_{μν}` antisymmetric ⇒ `F≡0` in 1D), and at the symmetric minimum `R*=0`, `F=0` by symmetry ⇒ gap `=Hess E_0=G″`. Discrete breaking gives only a tunneling splitting `e^{−cV}→0` in the continuum. No uniform `λ`.

## §3 · Does the indefinite Hessian break BO?
**Validity — no.** Exact `E_0` is finite and `C^∞` (compact fiber, bounded `S_W`); a non-Gaussian min–max fast-gap bound handles `H_ff`'s negative mode. **Gap — yes.** The Gaussian/Schur potential `½Tr log H_ff` is complex (PROOF-BORN-OPPENHEIMER §2; COARSE-POINCARE (a)(iii)); the inverted-oscillator WKB at the cone has no bound-state gap; negatives are 90–97% *fiber* (W6; FORGE5-SL-UV-CONE-EMPTY K2), violating the fast-gapped premise on the instability side; the Landau–Zener crossing makes `1/ε_1` and `sign Hess A` uncontrolled on the `k→0` goal mode.

## §4 · Exact residual
```
∃λ>0:∀R*, Hess E_0(R*) = G″(R*)|_{k→0} ≥ λ   (≡ H1),   AND
        inf_α ε_1(α) ≥ γ>0 through the cone       (now ~e^{−2βN_p}, W6/BE).
```
First line = H1; second = the marginal-gap/BE wall.

## Five-point (honest)
1. **Nonperturbative — PASS.** Exact compact-fiber integral + min–max/WKB fast-sublevel bounds; indefinite Hessian handled non-Gaussianly, no loop expansion.
2. **Sign lower bound — PARTIAL.** Genuine rigorous `ε_1(α)≥γ` where convex (new); **FAILS at the cone** (marginal + Landau–Zener) and gives no sign on the *gap* object `Hess E_0`.
3. **Non-circular — FAIL at closing.** `Hess E_0=G″`; `G″≥λ` is H1. WKB bounds the error `A`, not the gap.
4. **Non-SUSY + fails free field — PASS on value, with a sharp nuance.** Free: `S_W≡0`, fiber decouples, `E_0=const`, `Hess E_0=A=F=0` ⇒ `G″=0`. But the species's own fast sublevels `ε_1` are *nonzero for free too* (Laplacian gap `→0` continuum), so they **do not discriminate**; only `Hess E_0=G″` does (free 0, YM +0.27/+0.49). The discriminator is the one term WKB cannot reach.
5. **Uniform `≥λ` — FAIL.** No uniform `λ`; cone marginal + `Hess E_0`=H1.

**Sharpest line.** BO-with-WKB cleanly separates two objects the corpus conflates: the fast sublevel gap `ε_1` (=`Λ_UV`, relaxation) and the slow static stiffness `Hess E_0` (=`G″`). WKB genuinely certifies `ε_1` under an indefinite Hessian — a real new bound — but `ε_1` controls the BO *error*, never the *gap*; the gap is `Hess E_0`, untouched. The indefinite Savvidy Hessian, being 90–97% fiber, breaks the very *fast-gapped* premise BO needs at the `k→0` goal mode. Net: an exact separation and one certified new object, but `Hess E_0=G″≥λ` is H1 verbatim, and the cone's marginal fast gap is the BE wall (W6) in adiabatic clothes.

*Refs.* PROOF-BORN-OPPENHEIMER (BO reroutes to H1; Schur complex); PROOF-COARSE-POINCARE (a)(iii), H1; k0-CONSTRAINT-EFFPOT (`V_eff″=G″`); E1-FRG (exact FP/fiber supertrace); W6-BE-at-2.5 (negatives 90–97% fiber; BE wall); FORGE5-SL-UV-CONE-EMPTY (K2 cone high-momentum); C8-SPECIES-NOGO (`X` sign-varying); PROOF-SAVVIDY-EXACT-SIGN (`λ_min=−g|B|`); PROOF-JACOBIAN-SAVES (`G″=Δ+W+X`); W7-Gpp0 (+0.27/+0.49). External: Combes–Seiler, CMP **46** (1976) 1; Klein–Martinez–Seiler, Ann. IHP **49** (1985); Teufel, *Adiabatic Perturbation Theory in Quantum Dynamics* (2003); Berry, Proc. R. Soc. A **392** (1984) 45; Simon, *Trace Ideals & Their Applications* (min–max/Birman–Schwinger); Landau–Zener; Kato.
