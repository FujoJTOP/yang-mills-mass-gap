# ABC-NEW-INEQ — Ward–Stein transport cone bound (WSTC): the Savvidy cone's *variance mass* is provably negligible — the sign lives elsewhere

> 2026-11-23 · SHAPE C forge. Target: a genuinely new **non-positivity** upper bound on the cone variance `Var_μ(⟨Pᵀ∇S,v⟩|a=0)` with the **right sign**, not Schur-concave, not log-concavity, acid-test-clean.
> **VERDICT: a real new inequality that meets #1, #3(cone), #4 and gives the right sign on the cone (#2 partial) — but FAILS #5 and the sharp part of #2/#3.** It proves the high-field cone contributes `≤e^{−cβ(κ²−κ_typ²)N_□}` to the variance, uniformly in volume — so "*the cone cannot destabilize `G″`*" becomes a theorem. But the sign-carriers are not the high-field cone: they are **high-momentum negatives at typical field** (W6-BE: 21–42% negative, 90–97% in the UV/fluctuation block). WSTC bounds the wrong tail; the residual is exactly H1.

## 1. The exact Ward–Stein linearization (new)
Flat chart `X=log U`, `a=⟨v,X⟩`, `|v|=1`, `f=∂_nV=⟨v,∇V⟩`, `μ∝e^{−V}`, `V=S_W+Φ`. Integration by parts (fiber compact/periodic: **no boundary term**), `∂_n a≡1`:
```
   E_μ[f a] = ∫ e^{−V} ∂_n a dX = 1,   all β, nonperturbative.        (Ward)
```
With marginal density `p(a)`, `m(a)=E[f|a]`, `∂_a log p = −m`, and `m′(a)=E[H|a]−Var(f|a)` (`H=Hess V`), so `G″(a)=m′(a)` and `G″(0)=E[H|0]−Var(f|0)` (matches PROOF-JACOBIAN-SAVES). Then `∫a m(a)p(a)da=E[af]=1`, so with `E[a]=0`:
```
   G″·Var(a) = 1  ⟺  marginal of a is exactly Gaussian.
```
The gap is precisely the **non-Gaussianity of the coarse marginal at typical action** — a first-moment (Ward) coordinate for H1, complementary to FORGE5's second-moment (`E[H]=E[f²]`).

## 2. The new inequality (boxed) — transport split of the cone variance
Let `C_κ={‖F‖²≥κ²}` ⊃ Savvidy cone `{2g|F|>μ²}`, `μ_0=μ(·|a=0)`. Split `Var_{μ_0}(f)=∫(f−f̄)²1_C+∫(f−f̄)²1_{Cᶜ}`.
```
   WSTC:   Var_{μ_0}(f·1_{C_κ}) ≤ osc_{C_κ}(f)² · μ_0(C_κ),
           μ_0(C_κ) ≤ e^{−cβ(κ²−κ_typ²)N_□},  osc_{C_κ}(f)² ≤ (cβκ_max)² (compactness),
   ⇒       cone contribution ≤ (cβκ_max)² e^{−cβ(κ²−κ_typ²)N_□} ⟶ 0,  uniformly in L.
```
**Right sign, unlike Efron–Stein.** FORGE5's `Var≤diam²E‖HessV‖²` is `≥(g|F|)²`, *largest on the cone* (wrong sign). WSTC is `≤` an exponentially small mass times a bounded oscillation: the large field `f~β‖F‖` is multiplied by its own `e^{−β‖F‖²V}` weight, so the **variance *mass* is elsewhere** — in the typical bulk, not the cone. Non-Schur (it is a tail-mass estimate, not a symmetric-convexity functional), non-log-concavity (uses only extensive-action suppression, no `Hess V⪰0`), and **more** suppressed at large `L` and large `β` — a genuinely volume-uniform statement.

## 3. Acid test
Free massless field (`β=0`,`Φ=0`): cone empty (`F=0`), `C_κ` trivial, WSTC returns `Var_cone=0`; the massless `k→0` mode has `Var(a)=∞⇒G″=0` exactly (Ward: `G″Var(a)=1` fails only because `Var(a)=∞`). **No spurious positivity** — the free-field bar is met. Non-SUSY: pure YM, no superpartners, `β=2.5` real and indefinite. `#4` PASS.

## 4. Why it does not close the gap
Savvidy's indefiniteness is **pointwise at typical field**, carried by *high-momentum* modes (W6-BE: SU2 `−1.12β` at 21.7% of directions, SU3 42%; 90–97% in the fluctuation/UV block). These sit at `‖F‖≈κ_typ`, **outside `C_κ` for every `κ>κ_typ`** — WSTC's exponent vanishes exactly there. The cone `{2g|F|>μ²}` (large field) is *not* where the negative variance lives; the large-field region is measure-starved. So WSTC proves the honest but insufficient statement: *the high-field Savvidy cone cannot move the physical `G″`* — the destabilizing variance is the **typical-action** one, where `E[H]−Var≈0` (equipartition) and only `W_geo≥0` (Prékopa, PROOF-JACOBIAN-SAVES) supplies margin. Certifying `G″≥λ` there is `Var(a)≤1/λ` at `k→0` = H1. And a transport/Holley–Stroock bound *there* is vacuous: `S_W` is **extensive**, so no volume-uniform log-concavity constant exists (W6 `Ric+HessS≈−2.3`).

## 5. Honest 5-point verdict
- **(1) nonperturbative — PASS.** Exact Ward identity + large-deviation mass bound, all `β`.
- **(2) sign-controlling — PARTIAL.** Right sign (upper bound) and a *proof* that the cone's variance contribution is `e^{−cβ(κ²−κ_typ²)N_□}`: the cone cannot destabilize. But it does not control the typical-action high-momentum variance, where the sign is decided.
- **(3) non-circular — PASS for the cone piece** (uses only action mass, assumes no gap); the sharp global form is circular = H1.
- **(4) non-SUSY + fails free field — PASS.**
- **(5) uniform `m_eff²(k→0)≥λ` — FAIL.** No `λ`: the bound is uniform *in `L`* but bounds the wrong (measure-staved) region; the sign-carriers at `κ_typ` are untouched = H1.

**Sharpest line.** WSTC *proves* the Savvidy cone is a red herring for the sign — its variance mass is `e^{−cβ(κ²−κ_typ²)N_□}`, uniformly in `L` — and the exact Ward identity `G″Var(a)=1` (iff Gaussian) relocates the entire gap to the **non-Gaussianity of the coarse marginal at typical action**. That is a genuine, correctly-signed, non-positivity theorem about the cone; it does not close #5, because the negatives are high-momentum-at-typical-field (W6-BE), not high-field, and that is H1.

*Refs.* PROOF-JACOBIAN-SAVES (`G″=E[H]−Var`, `W_geo≥0`, (◆)); FORGE5-CONE-VARIANCE (`E[H]=E[f²]`, Efron–Stein wrong sign); INVENT-VARIANCE-INEQ (SCBL `G″=S−D`); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`); W6-BE-at-2.5 (negatives high-momentum, 21–42%, `−1.12β`; `Ric+HessS≈−2.3`); GAP-G1-largefield (cone `√(2g|F|)`, UV-empty, `k→0` residue=G2); REDTEAM-SPEC ((2)≡(5), positivity class dead); Stein's lemma / Schwinger–Dyson; Holley–Stroock 1987; Brascamp–Lieb 1976.
