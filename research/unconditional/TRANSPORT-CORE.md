# TRANSPORT-CORE — the Bakry–Émery Γ₂ / Otto–Wasserstein species: its native output *is* a gap lower bound (the right anatomy), but the integrated form degenerates on the affine family to `G″⪰ρ` exactly, and the composition `W_geo≥0 + WSTC` covers the complement of the gap's cell — the eighth species is a faithful dialect of H1, score 2/5

> 2026-11-23 · CORE attack on the candidate **8th species**: Bakry–Émery `Γ₂`-calculus / optimal-transport (Otto–Wasserstein) geometry of `μ = e^{−S}dU`. This is the one positivity-type machinery whose **native output is a lower bound on a spectral gap** (the carré-du-champ `Γ₂≥ρΓ ⇒ gap ≥ ρ`), i.e. the anatomically correct object the Clay problem needs. **判词 (verdict, up front): NO NEW LEVER — score 2/5.** It passes (1) nonperturbative and (4) free-field, but fails (2) sign-control, (3) non-circularity, (5) uniformity. The hoped-for composition ("`W_geo≥0` + `WSTC` starvation combine via the *integrated* `Γ₂` into `E[Γ₂]≥ρ`") is **false**: on the coarse marginal the integrated criterion is **equivalent** (no slack) to `Hess V_eff⪰ρ ⟺ G″⪰ρ ⟺ H1`, because the affine coarse family kills the one positive term `‖Hess f‖²`; and the two proven pieces live on the **complement** of the sign-carrying cell (typical field × `k→0`). The single obstruction — *anti-localization of every positive compensation* (`‖Hess f‖²` at high `k`, `Φ″` surplus at the log-cut `d→2π`, WSTC at high field) — is one obstruction in three guises. Exact stuck step = one conditioned-variance bound = `(◆)` = H1.

---

## §1 · The `Γ₂` decomposition, written out

**Setting (真定理 / exact).** Lattice connection space `X = G^E/G_site`, flat chart `X_ℓ = log U_ℓ ∈ su(N)`, `dU = J dX`, `Φ = −log J = Σ_{i<j} g(λ_i−λ_j)`, `g(d) = −2log(sin(d/2)/(d/2))`. Wilson measure `μ = Z⁻¹ e^{−S_W} dU = Z⁻¹ e^{−V} dX` with the **effective flat potential**

```
V = S_W + Φ,    S_W = (1/2g²)∫tr F²  = β Σ_□ (1 − Re tr U_□/N).
```

The reversible Langevin/OU generator on functions is `L = Δ_gauge − ∇V·∇`, and the carré du champ is

```
Γ(f) := ½ L(f²) − f Lf = |∇f|²_{G}      (G = DeWitt / bi-invariant metric on X).
```

The **Bakry–Émery `Γ₂` operator** is `Γ₂(f) := ½ L Γ(f) − Γ(f, Lf)`, and the Bochner/Pratt identity gives the **decomposition**

```
┌───────────────────────────────────────────────────────────────────────────┐
│  Γ₂(f) = ‖ Hess f ‖²_{G}  +  Ricₕₐₐᵣ(∇f,∇f)  +  Hess_S(∇f,∇f)               │   (Γ₂)
│          ╰ always ≥0 ╯      ╰ geometric, ≥0 ╯   ╰ action, INDEFINITE ╯         │
└───────────────────────────────────────────────────────────────────────────┘
```

where `Ricₕₐₐᵣ` is the Ricci curvature of the configuration manifold and `Hess_S` is the action Hessian `= β·Hess(−Re tr U_□/N)` per plaquette.

**The flat-chart reconciliation (真定理, key for §2).** In the `X = log U` chart the group-manifold Ricci `Ricₕₐₐᵣ` (the W6/W6-BE object) and the Jacobian Hessian `Hess Φ` (the PROOF-JACOBIAN-SAVES object) are **the same geometric content**, related by the Riemann-volume expansion `Hess(−log J)|_{0} = −⅓ Ricₕₐₐᵣ` (the standard `J = 1 + ⅙Ric_{ij}x^ix^j + …`). So `Ricₕₐₐᵣ ↔ Hess Φ` up to the `O(1)` constant, and `Hess V = Hess S_W + Hess Φ`, so

```
Γ₂(f) = ‖ Hess f ‖²_G  +  Hess S_W(∇f,∇f)  +  Hess Φ(∇f,∇f).
                                    ╰ indefinite, Savvidy ╯   ╰ ≥0, Prékopa "W_geo" ╯
```

Two facts about the third term are theorems and are the corpus's load-bearing pieces:
- **Pointwise floor (Lemma B).** `Φ″ ⪰ (N/12)` per link (`g″ ≥ 1/6`, SU2 exact; measured `E[Φ″] = 0.34` SU2 / `0.46` SU3).
- **Prékopa (真定理).** `Φ` convex ⇒ `e^{−Φ}dX` log-concave ⇒ pushed coarse free energy convex ⇒ `W_geo := E[Φ″] − Var(Φ′) ≥ 0`, **every mode incl `k→0`, every `L`** (PROOF-JACOBIAN-SAVES).

**The exact marginal identity (真定理; PROOF-JACOBIAN-SAVES, ABC-NEW-INEQ).** Push `μ` to the coarse coordinate `a = ⟨v, log U⟩` (`|v|=1`, `k→0` mode); let `m(a) = E[⟨v,∇S⟩|a]`. Then `m′(a) = E[H|a] − Var(⟨v,∇S⟩|a)` with `H = Hess S_W`, and `G″ := V_eff″ = m′(0)`:

```
G″(0)[v,v] = E[H + Φ″|a=0][v,v] − Var(G + Φ′|a=0)[v,v]
           = Δ_dyn  +  W_geo  +  X,        Δ_dyn = E[H]−Var(G),  X = −2Cov(G,Φ′).   (G″-split)
```

Unconditioned at the `P`-invariant vacuum `E[H] = E[⟨v,∇S⟩²]` (equipartition — ABC-NEW-INEQ), so `G″` **is** the conditioning defect; the Ward–Stein identity reads `G″·Var_μ(a) = 1` **iff** the coarse marginal is exactly Gaussian (measured to 2–5%).

This is the whole machine. Everything below tests whether `Γ₂` can turn (Γ₂) into `≥ρ>0` on YM.

---

## §2 · The integrated-`Γ₂` question: can `W_geo≥0` + `WSTC` combine to `E[Γ₂] ≥ ρ`?

The **hypothesis under test** (stating it at full strength): pointwise `Γ₂ ≥ ρΓ` fails at the Savvidy cone (W6-BE: `min ≈ −2.0`), but perhaps the *integrated / measure-weighted* form `E_μ[Γ₂(f)] ≥ ρ E_μ[Γ(f)]` survives, because (a) the geometric part `Hess Φ/W_geo ≥ 0` is proven, and (b) the negative action part is *measure-starved* on the cone (`WSTC`), so averaging kills it. If true, this composes two proven theorems into the first genuine gap lower bound.

### §2.1 The integrated criterion is sufficient (真定理 — this is the honest reason the hope is tempting)

**Lemma (energy decay).** Fix `ρ>0`. Suppose, for every `f ∈ C^∞`,

```
   (ICD)      E_μ[Γ₂(f)]  ≥  ρ · E_μ[Γ(f)].
```

Then `μ` satisfies `Poincaré(ρ)`: `Var_μ(f) ≤ (1/ρ) E_μ[Γ(f)]` for all `f`, i.e. the spectral gap of `L` is `≥ ρ`.

*Proof.* `μ` invariant ⇒ `∫ L g dμ = 0` for all `g`. With `f_t = P_t f`, `E(t) := E_μ[Γ(f_t)]`, use `Γ₂(f_t) = ½ LΓ(f_t) − Γ(f_t, Lf_t)`:
`E_μ[Γ₂(f_t)] = ½∫LΓ(f_t)dμ − ∫Γ(f_t,Lf_t)dμ = −½E′(t)`. So `E′(t) = −2E_μ[Γ₂(f_t)] ≤ −2ρ E(t)` by (ICD) applied to `f_t`; hence `E(t) ≤ e^{−2ρt}E(0)`. Finally `−d/dt Var(f_t) = 2E(t)`, so `Var_μ(f) = ∫₀^∞ 2E(t)dt ≤ E(0)/ρ = (1/ρ)E_μ[Γ(f)]`. ∎

So an *integrated* `Γ₂` bound is a **bona fide, non-circular gap-lowering criterion** — strictly weaker than pointwise `Γ₂≥ρΓ` (the `‖Hess f‖²` term is allowed to contribute only on average). This is exactly the mechanism `W_geo≥0 + WSTC` would exploit. The question is whether (ICD) is actually reachable.

### §2.2 No escape: on the coarse marginal (ICD) is *equivalent* to `H1` (真定理)

**Lemma (affine kills the only slack).** Let `ν = P_*μ` be the coarse marginal with effective potential `V_eff` (`Hess V_eff = G″`) on **flat** coarse space. The coarse generator `L_ν = Δ_c − ∇V_eff·∇_c` has `Γ^c(φ)=|∇φ|²`, `Γ₂^c(φ)=‖Hess φ‖²_{flat}+Hess V_eff(∇φ,∇φ)`. Then

```
   (ICD)_c   ⟺   Hess V_eff ⪰ ρ      ( = G″ ⪰ ρ = m_eff²(k→0) ≥ ρ = H1 ).
```

*Proof.* (⇐) If `G″⪰ρ` then `Γ₂^c(φ) = ‖Hess φ‖² + Hess V_eff(∇φ,∇φ) ≥ ρ|∇φ|² = ρΓ^c(φ)` pointwise, so (ICD)_c holds. (⇒) Conversely, take `φ(a) = ⟨v,a⟩` **affine**: `Hess φ ≡ 0`, so `Γ₂^c(φ) = Hess V_eff(v,v)` *deterministically* and `Γ^c(φ)=|v|²`; (ICD)_c reads `Hess V_eff(v,v) ≥ ρ|v|²` for every `v`, i.e. `G″⪰ρ`. ∎

The point: the always-`≥0` term `‖Hess φ‖²`, the *only* thing that could let an average beat a pointwise negative, **vanishes on the affine family** — and the affine family is in the test set (the criterion must hold for all `f`, in particular linear ones). So on the coarse marginal the integrated criterion buys **zero** slack. It is H1, lossless.

**And the full-space criterion cannot be easier either (pushforward sandwich).** For any coarse function `ψ(U′)`, `Var_μ(ψ∘π) ≤ (1/ρ)E_μ(|∇(ψ∘π)|²) = (1/ρ)E_ν(|∇_cψ|²)` under (ICD) on the full generator. Hence **full `Γ₂`-gap `≥ρ` ⇒ coarse Poincaré(ρ) = H1**. Combined with §2.2: a successful integrated-`Γ₂` proof on the full YM measure would *itself prove H1*. There is no reading — coarse or full, pointwise or averaged — in which `Γ₂` is weaker than H1 where the gap lives. (This sharpens B3-FUNCTIONAL-BE's "lossless on flat space" from a remark to a two-sided equivalence.)

### §2.3 The composition `W_geo≥0 + WSTC` does *not* reach (ICD) — and exactly why

By §2.2 the only question is whether the composition supplies `G″⪰ρ` on the coarse marginal. Use (G″-split) on the `k→0` mode and insert the two proven pieces:

```
   G″ = E[Φ″|0]  +  E[H|0]  −  Var(G+Φ′|0)
      ≥  N/12     +  E[H|0]  −  Var_cone(G|0)  −  Var_typ(G|0) .
        ╰ ≥0 ╯     ╰ sign? ╯   ╰ WSTC: ≤ e^{−cβ(κ²−κ_typ²)N_□} → 0 ╯   ╰ UNCONTROLLED ╯
```

- `E[Φ″|0] ≥ N/12` is **true** (Lemma B) — the geometric positivity the hypothesis relies on. ✔
- **WSTC** (真定理, ABC-NEW-INEQ) bounds the **cone** part of the force variance: `Var_μ(f·1_{C_κ}) ≤ osc² · e^{−cβ(κ²−κ_typ²)N_□}`, uniformly in `L` — so `Var_cone → 0`. ✔
- **But** the sign-carriers are **not** on the cone. W6-BE: the Savvidy negatives are **high-momentum at *typical* field** — 21.7% of directions (SU2), 42% (SU3), `−1.12β`, **90–97% in the fluctuation/UV block**, at `‖F‖ ≈ κ_typ`, **outside `C_κ` for every `κ > κ_typ`**. WSTC's exponent vanishes exactly there. So `Var_typ(G|0)` — the typical-field force variance — is **untouched** by WSTC. ✘
- `E[H|0]` (the conditioned mean of the *bare* action Hessian) has **uncontrolled sign**: pointwise bare coarse `PᵀHP` min `≈ −0.3` (MARGINAL); the measured `+0.22` is a *measure* fact, not a theorem. ✘

Hence the composition yields `G″ ≥ N/12 − Var_typ(G|0) − (sign-uncontrolled E[H|0])`, and since `Var_typ` is unbounded without H1 and the `E[H|0]` sign is open, **no `ρ>0` follows**. The residual line is literally

```
   Var_typ(G|0)  ≤  E[Φ″|0] + E[H|0] + (cone, →0) − ρ      ⟺   Var_μ(a_{k→0}) ≤ 1/ρ   ⟺ H1.
```

**Verdict of §2.** The integrated form is sufficient (Lemma 2.1) — tempting — but it is *exactly* H1 on the coarse marginal (§2.2), and the two proven pieces `W_geo≥0` + `WSTC` cover the **complement of the sign-carrying cell**: `W_geo` owns the `|F|`-free Haar floor `N/12`, `WSTC` owns the high-field cone, `‖Hess f‖²` owns high `k` — and the gap lives at **typical field × `k→0`**, where all three vanish. The pieces do **not** compose into `E[Γ₂]≥ρ`; they compose into `G″ ≥ N/12 − (uncontrolled typical-field variance)`.

---

## §3 · The two acid tests, decisively

### §3.1 Acid test (1) — does it escape the Savvidy-sign / positivity-circularity wall?

**Outcome: FAIL as a new lever; but it is the *first* species whose output-type is correct.** Two separate questions:

- *(1a) Direction.* `Γ₂` is built *on* positivity (`Γ` is a variance density, `Γ₂ = ½LΓ−Γ(f,Lf)`), and `L` is elliptic, so the native inequality is a **lower** bound on the gap — the **correct direction**, unlike all six dead positivity species (which give variance **lower** bounds / decay **upper** bounds and are blind to the atom via the modulus-square sign). **This is the genuine qualitative advance of the species and the reason it earns 2 rather than 0–1:** its output *is* the target-type object.
- *(1b) Escape.* It does **not** escape. By §2.2 the criterion is *equivalent* to `G″⪰ρ = H1` on the coarse marginal — a faithful **rename**, zero slack (matches B3-FUNCTIONAL-BE §5(3): "faithful dialect, not a new functional"). Physically it degenerates *to* the Savvidy-sign wall, not past it: the failure is at `v₀` = the `k→0` Savvidy-negative coarse mode (`PᵀHP` min `≈−0.3`), where `E[Φ″|0]≈N/12` does not beat `−E[H|0]+Var(V′|0)+ρ` (B3 §4; coarse slack `−0.08`).

**Result: (1) FAIL for the *composition* hypothesis, PASS for *direction*.** The integrated/measure-weighted form is the strongest possible version of the hope, and Lemma 2.1 shows it is sufficient — but Lemma 2.2 shows it is *lossless*: no slack, hence no escape.

### §3.2 Acid test (2) — free MASSLESS gauge field (gapless, `f₁²=0`)

**Outcome: PASS in the continuum, VACUOUS on the lattice — and that split is itself the decisive finding.**

- **Continuum free massless field.** `S_0 = (1/2g²)∫tr F² = ½∫A(−Δ)A` on transverse modes; `Hess S_0 = −Δ ⪰ 0` (PSD, **zero mode at `k=0`**). Continuum `R⁴` has Lebesgue measure ⇒ **no compact Jacobian** (`Φ ≡ 0`). So `Γ₂ = ‖Hess f‖² + 0 + (−Δ)(∇f,∇f) ⪰ 0`, and its infimum over modes is at the zero mode ⇒ `ρ = 0`. **Returns exactly `0`, not `>0`.** ✔ This is the desired free-field failure: `BE` correctly distinguishes massless (`ρ=0`) from massive Proca (`Hess = −Δ+m² ⇒ ρ=m²>0`).

- **Lattice free field — VACUOUS.** On the compact group the Haar Jacobian is present *at every β*: `Φ″(λ→0) → 1/6 ≥ N/12 > 0` (the eigenvalue-collision floor, Lemma B — it does **not** vanish for small fluctuations). So a lattice `Γ₂` computation returns `ρ ≥ N/12 ≈ 0.25 (SU3) / 0.17 (SU2)` **on the free field**. By the charter's acid rule this is **vacuous**: the tool certifies the **compactness / kinematic floor**, not a particle. The corpus already flags precisely this (registry A5 "只认证自由格点隙（非空）"; INVENT-NEW-MONOTONE §4: "forces only `λ_min(a,L)→0`… no thermodynamic gap"). And it is dimensionally wrong: `G″_a = a²m_phys²`, so a *constant* `G″ ≥ N/12` at fixed `a` means `m_phys = √λ/a → ∞` (a UV artifact, PROOF-SKELETON §0 量纲更正).

**The decisive paradox (the sharpest acid-test content).** The species' **only positivity** — `W_geo≥0`, the Haar floor `N/12` — is **nonzero on the free field** (on the lattice, at all β). So the positivity that would "save" the YM sign and the failure-mode that makes the test vacuous are **the same object**: the tool cannot distinguish the **physical gap** from the **compact-group kinematic mass**. To dodge this one must *relocate* the certified positivity to the **conditioning** `E[H|0] = +0.22` (B3 §5(4)), i.e. to the *action* — which is Savvidy-indefinite and controlled only by `Var(G|0)` = H1. This is the **Schur-convex signature in Riemannian dress**: the only positivity available is one that is nonzero at `g=0`, hence cannot be the atom.

**Result: (2) PASS in the continuum, VACUOUS in the lattice regularization actually used.** The species is "half-vacuous": its action part passes the free test, its Haar part *is* the free-field failure.

---

## §4 · The Savvidy cone in `Γ₂` language — is it fatal?

**The cone Hessian (真定理 / exact).** `Hess S_W = −D² + 2 ad F_B` on a background `F_B`; the softest unstable direction is transverse-magnetic at `k_∥ = 0`, with `λ_min = −g|B| → −g²(μ)|F|` (PROOF-SAVVIDY-EXACT-SIGN; Savvidy 1977, Nielsen–Olesen 1978). In `Γ₂`:

```
   Γ₂(f) = ‖ Hess f ‖²_G  +  Hess Φ(∇f,∇f)  +  ⟨∇f, (−D² + 2 ad F_B) ∇f⟩ .
```

On a plane-wave test `f = e^{ikx}`, `‖Hess f‖² ∼ k⁴` while the cone term `∼ −g|B|·k²`. Therefore:

| regime | dominant term | `Γ₂` sign | gap? |
|---|---|---|---|
| high `k` (`k² ≫ g|B|`) | `‖Hess f‖² ∼ k⁴` | **`>0`** (compensated) | irrelevant — gap is IR |
| `k → 0` (the marginal mode) | `‖Hess f‖² → 0`; `Γ₂ ≈ Hess Φ − g|B|` | **can be `<0`** | **this is where the gap lives** |

**Is it fatal? YES — and for a structural reason, not a missing estimate.** The one positive term that could beat the negative cone, `‖Hess f‖²`, is **anti-localized to high `k`** — largest exactly where the mass gap is *not* decided. At `k→0` (the lightest mode, where `ρ` is set) it vanishes, leaving `Γ₂ ≈ Hess Φ − g|B|`. And `Hess Φ`'s invariant content is the `|F|`-free floor `N/12`; its **surplus** convexity (`g″→∞`) sits only at the **log-cut** `d→2π` (links near the center element `zˣ`, where `U_□ ≈ 1` hence `F≈0`) — **anti-correlated with the driver `∝|F|`** (PROOF-SAVVIDY-WGEO-DOMINATION: no gauge-invariant norm makes `W_geo` dominate with `c<1`; spectral ratio `ρ≈3–16`; only `β < 0.15` works). So at the gap mode: `Hess Φ` sees only `N/12`, the cone sees `−g|B|`, and the `−Var(G|0)` conditioning term subtracts further.

**The unifying obstruction (this file's synthesis).** Every positive compensation in the `Γ₂` decomposition is **anti-localized to the complement of the gap's cell** (typical field × `k→0`):

| positive term | large where | gap cell where it vanishes |
|---|---|---|
| `‖Hess f‖²` | high momentum `k ≫ √(g/β)` | `k→0` |
| `Φ″` surplus (`→``∞`) | log-cut `d→2π` (`F≈0`) | per-link phase `O(π)` (Savvidy on) |
| `WSTC` (only *suppresses*) | high field `‖F‖ > κ_typ` | typical field `‖F‖ ≈ κ_typ` |

All three are largest where the gap is *not*, and all three vanish on the *one* cell (typical field, `k→0`, `O(π)` phase) where the sign is decided. **Anti-localization of positivity is the single wall, seen in three costumes.** This is why the composition `W_geo + WSTC` cannot close: not a magnitude shortfall (`D/S≈2%` is *large* margin on average), but a **support mismatch** — the pieces and the sign live on disjoint regions.

**Is the Savvidy cone fatal *as a no-go*?** It is **physical, not gauge** (transverse-magnetic `δF ≠ 0`, survives the gauge quotient, W4 measures 7.1% negative after projecting) — so it cannot be gauged away. But it is **not fatal to the *sign*** in the corpus's measured reality: the exact marginal `G″ = +0.51` (SU3 β=2.5, W7 Part B), `1/Var` volume-independent — the sign is **carried by the conditioning**, `E[H|0]=+0.22` + `W_full≈+0.54`, i.e. by `Var(G|0)`. The cone is fatal to **certifying** the sign in `Γ₂` language because the certification needs `Var(G|0)`, which is exactly the piece `Γ₂` rewrites rather than bounds.

---

## §5 · Scorecard, sharpest line, exact stuck step

### Five-condition grade (0–5): **2/5**

| # | condition | grade | one-line basis |
|---|---|---|---|
| 1 | **nonperturbative** | **PASS** | exact `Γ₂` identity (Γ₂), exact G′-split, Prékopa/BL, Savvidy spectrum; no loop, no `g`-expansion, all β |
| 2 | **sign-controlling LOWER bound** | **FAIL** | native output is a lower bound (right anatomy), but at β=2.5 the provable `ρ ≤ 0` (pointwise `−2.0`, coarse `−0.3`); `ρ≥N/12` only as a lattice artifact — no `λ>0` |
| 3 | **non-circular** | **FAIL** | on coarse marginal `Γ₂`-criterion ≡ `G″⪰ρ` ≡ H1 (Lemma 2.2, lossless); full form ⇒ H1 (pushforward) — a faithful rename |
| 4 | **non-SUSY + fails free massless field** | **PASS (caveat)** | continuum free field ⇒ `ρ=0` ✔; lattice free field ⇒ `ρ≥N/12` (vacuous) — half-vacuous |
| 5 | **`(L,a)`-uniform `m_eff²(k→0)≥λ`** | **FAIL** | dimension-correct uniform bound is H1; constant `G″` is dimensionally a UV artifact |

**Score = 2/5** (conditions 1, 4 pass; 2, 3, 5 fail). It is the most *structurally honest* 2 in the corpus: the **right anatomical type** (a genuine spectral-gap lower bound, the only species with that output-type) whose **value degrades exactly to H1**.

### Sharpest line
`Γ₂` is not the eighth species — it is the **genotype** of H1. On flat coarse space it is *lossless* (`Γ₂^ν = G″ = H1`), the integrated form degenerates on the affine family (Lemma 2.2), and the composition `W_geo≥0 + WSTC` covers the **complement** of the gap's cell; worse, the one positivity it owns (the Haar floor `N/12 ≈ 0.25`) is **nonzero on the free field**, so the tool's positivity and its free-field failure are the *same object* — the Schur-convex signature in Riemannian dress. The only species whose output *is* a gap lower bound turns out to be the one whose lower bound *is* the gap: `ρ` collapses to H1 not by a missing estimate but by an identity.

### EXACT stuck step (缺失步)
> **Control the conditioned coarse force variance on the `k→0` typical-field mode:**
> ```
>   Var_μ( ⟨v,∇S⟩ | a=0 )  ≤  E_μ[H_vv + Φ″|a=0]  −  ρ ,     uniform (L, a),  v = v₀ (softest),
> ```
> equivalently (Ward–Stein) `Var_μ(a_{k→0}) ≤ 1/ρ` with the correct asymptotics `ρ_a ∼ a²Δ₀²`. This is **exactly `(◆)` of PROOF-JACOBIAN-SAVES = H1 = the gap itself.**
>
> Every component of `Γ₂` is a faithful rewriting of this one inequality: `Hess Φ/W_geo` supplies its `|F|`-free floor `N/12`; `WSTC` supplies its high-field-cone part; `‖Hess f‖²` supplies its high-`k` compensation; and the residual cell — **typical field × `k→0`** — is the conditioned variance itself. The species does not cross the wall; it *is* the wall, viewed through a Riemannian lens.

### What would change the verdict (honest)
Not a sharper estimate — a **support relocation**: (i) a proof that the `k→0` negative mode is *cone*-localized (contra W6-BE, which puts it at typical field — this is the only lever that would let `WSTC` reach the sign); or (ii) an *upper* bound `Var_typ(G|0) ≤ N/12 + E[H|0] − ρ` that does **not** presuppose `Var(a)≤1/ρ` (this is H1; C5/C6 prove the nine dependency-inequality classes all fail to supply it for frustated, non-attractive SU(N)). Absent one of these, the eighth species is closed as a **faithful dialect of H1** — a real theorem-shaped object (`Γ₂` is a correct gap-lowering criterion) whose value is the open problem.

---

*Refs.* **Internal:** B3-FUNCTIONAL-BE (`Γ₂` on `ν` = lossless `G″=H1`; anti-localization; coarse slack `−0.08`); W6-BE-at-2.5 (pointwise `min −2.0`, 21–42% neg, 90–97% fiber); OTHER-H1-POINCARE (BE≡H1 two-route no-go; `PᵀHP` min `−0.3`); INVENT-NEW-MONOTONE (`Γ₂` sign inversion; free forces `λ_min→0`); TRANSPORT-MACHINE (Stage-3 = H1/I2; factorization); P3-wasserstein-ym (CD open, per-link negatives UV-localized); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `Φ″≥N/12`, Prékopa `W_geo≥0`, `(◆)`); ABC-NEW-INEQ (Ward `G″Var(a)=1`; WSTC `≤ e^{−cβ(κ²−κ_typ²)N_□}`; negatives are typical-field high-momentum, *not* high-field); PROOF-SAVVIDY-WGEO-DOMINATION (`c<1` false, `ρ≈3–16`, log-cut anti-correlation); OTHER-SAVVIDY-SIGN (coarse = susceptibility ceiling, not magnitude); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF_B`, `λ_min=−g|B|`); REDUCTION-AUDIT (gap⇒H1 one-way; `H*=I1∧I2` strictly stronger); PROOF-SKELETON (`G″_a=a²m²` dimension correction); new-tools-registry A5 (W-entropy, free-field flag). **External (each verified to the extent cited in the corpus):** Bakry–Émery, *Diffusions hypercontractives*, LNM 1123 (1985); Bakry–Gentil–Ledoux, *Analysis and Geometry of Markov Diffusion Operators*, Grundlehren 348 (2014) — `Γ₂`, Bochner/Pratt, the integrated energy-decay criterion; Otto, Comm. PDE 26 (2001) 101; Sturm, Acta Math. 196 (2006) 65; Lott–Villani, Ann. Math. 169 (2009) 903; von Renesse–Sturm, CPAM 58 (2005) 1014; Prékopa, Acta Sci. Math. 34 (1973) 335 (log-concavity); Brascamp–Lieb, J. Funct. Anal. 22 (1976) 366; Holley–Stroock, *Logarithmic Sobolev inequalities and stochastic Ising models*, J. Stat. Phys. 46 (1987) 1159; Savvidy, Phys. Lett. B 71 (1977) 133; Nielsen–Olesen, Nucl. Phys. B 144 (1978) 376; Shen–Zhu–Zhu, CMP 400 (2023) 805, arXiv:2204.12737 (volume-uniform LSI at `β<1/48`).
