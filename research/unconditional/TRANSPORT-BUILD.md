# TRANSPORT-BUILD — Bakry–Émery Γ₂ composition for 4D pure SU(N) YM: the right-signed tool, and the exact place it stops

> 2026-11-23 · CONSTRUCTOR build of the **Bakry–Émery Γ₂-calculus / optimal-transport** species (§1–2 formulas), the **key composition** `Γ₂ = Ric_Haar + Hess_S` against `W_geo≥0` and WSTC starvation (§2), three honest regime tests (§3), a 5-condition scorecard (§4), sharpest line and the exact stuck step (§5).
>
> **VERDICT.** Bakry–Émery is the **directionally-correct** tool (a lower bound on the gap — the only species here that is not positivity-closed). Its only *usable* form is the **pointwise** criterion `Γ₂(f) ≥ ρΓ(f) ∀f`, because the **measure-weighted** form is not a weakening of the gap — it is **an identity for it**: `∫Γ₂(f)dμ = ‖Lf‖²_μ` and `inf_{f⊥1} ‖Lf‖²/‖∇f‖² = λ₁` (both proved below and numerically checked). The pointwise form holds rigorously in a strong-coupling window (Ric dominates, `ρ=κ_H−cβ>0`, volume-uniform — consistent with P3/Lemma SC) and returns exactly `ρ=0` on the free massless field (acid test met). At the **crossover β≈2.5 it dies**, and it cannot be saved: the Haar–Prékopa budget `W_geo≥0` is real but is a *measure-average*, the Savvidy negatives are *pointwise and physical* (`−g|B|`, 22–49% of directions), and **WSTC starvation does not cover them** because the negative variance mass sits at **typical field**, where the starvation factor is exactly `1`. The composition therefore lands, exactly, on the existing wall: `G″ = E[H] − Var(G) + W_geo + X ≥ ρ`, i.e. **`E[Hess_S] ≥ ρ − W_geo` — the Savvidy-sign wall with `W_geo` as a (real but insufficient) budget — which is H1.** Score **2.5 / 5**. Not a proof; not vacuous; a clean no-go and a precise reduction.

---

## 1. The generator, the carré du champ, Γ₂

**Setting.** Configuration space `M = G^E / G_site` (gauge-fixed `G^E` equivalent), product of `|E|` copies of the compact group with the **bi-invariant metric** `g` (scaled `−Tr`, so `Ric>0`). The YM measure is `μ = Z^{-1}e^{−S_W}dU`, `dU =` product Haar, `S_W = (1/2g²)∫tr F² = β Σ_□ (1 − (1/N)Re tr U_□)`.

**Langevin generator** (μ-symmetric, negative semidefinite):
```
L f = Δ_gauge f − ⟨∇S_W, ∇f⟩,          L = e^{S_W}∇·(e^{−S_W}∇·).            (1.1)
```
**Carré du champ** (Leibniz, μ-symmetric):
```
Γ(f,g) = ½[ L(fg) − f Lg − g Lf ] = ⟨∇f, ∇g⟩,      Γ(f) := Γ(f,f) = |∇f|².   (1.2)
```
**Iterated carré du champ / Γ₂** (Bochner–Lichnerowicz):
```
Γ₂(f) = ½ L Γ(f) − Γ(f, Lf)
      = |Hess_g f|² + Ric_M(∇f,∇f) + Hess S_W(∇f,∇f)                          (1.3)
      = |Hess_g f|² + [Ric_Haar](∇f,∇f) + [Hess S_W](∇f,∇f).
```
Two charts — and this is where the species splits, so state it once:
- **Intrinsic** (above): geometric term = `Ric_Haar` (bi-invariant `SU(N)`: `Ric = ½g` (SU2), `¾g` (SU3) per unit; a β-independent **positive** constant `κ_H`), action term = `Hess S_W`.
- **Flat chart** `X = log U`, `dU = J dX`, `Φ = −log J`: `μ ∝ e^{−(S_W+Φ)}dX`, so the generator is `L = Δ_flat − ⟨∇V,∇·⟩` with `V = S_W + Φ` and `Ric_flat = 0`; the geometry is now `Hess Φ`. The intrinsic `Ric_Haar` and the flat `Hess Φ` are **the same geometric budget in two charts** (the difference is the connection term, lower order and integrating away). The numerics (W6, W7) and Prékopa live in the flat chart; §2 uses both.

**Bakry–Émery criterion (explicit).**
```
   (BE)   ∃ρ>0 :  Γ₂(f)(x) ≥ ρ Γ(f)(x)   ∀f∈C^∞_c(M), ∀x∈M
   ⟹  spectral gap of −L on L²(μ)  ≥ ρ,
        i.e.  Var_μ(f) ≤ (1/ρ)∫|∇f|²dμ  and  ‖P_t f − μ(f)‖ ≤ e^{−ρt}‖f‖.
```
Equivalently `(BE) ⟺ Ent_μ(f²) ≤ (2/ρ)∫|∇f|²dμ` (LSI(ρ)) ⟺ entropy `ρ`-convex along `W₂` geodesics (CD(ρ,∞)). Proof is the Bochner identity + Gronwall on `Γ(P_t f)`, Bakry–Émery [LNM 1123; Otto 2001; Sturm 2006; Lott–Villani 2009; von Renesse–Sturm 2005].

**Direction (the whole point of the species).** `(BE)` yields `m_eff² ≥ ρ`: a **LOWER bound on the gap**, `ρ` enters on the *good* side. This is why the species is not already in the "positivity ⇒ upper bound / vacuum" immune class — it is the one correct-signed lever. The price is that `(BE)` is only *sufficient*: `Γ₂ ≥ ρΓ` is strictly stronger than `λ₁ ≥ ρ` (a measure can have a gap and a non-convex `Γ₂` somewhere). That strength is exactly what §2 is forced to give up — and giving it up costs everything.

---

## 2. The key composition attempt: `Γ₂ = Ric_Haar + Hess_S`

### 2.1 The two proven halves, and the indefinite core
- **Haar geometry `≥0` (theorem, Prékopa).** `Φ = Σ_{i<j} g(λ_i−λ_j)`, `g(d) = −2log(sin(d/2)/(d/2))`, `g″ ≥ 1/6` on `|d|<2π` ⇒ `Φ` convex ⇒ `e^{−Φ}dX` log-concave ⇒ (Prékopa) every linear marginal has convex free energy:
  ```
  W_geo := E_μ[Φ″] − Var_μ(Φ′) ≥ 0,   all modes, all k (incl. k→0), all β.      (2.1)
  ```
  [PROOF-JACOBIAN-SAVES; Brascamp–Lieb; Prékopa.] In intrinsic language this is the measure-average of `Ric_Haar`; in flat language it is the Haar-Jacobian curvature. Either way it is a **non-negative budget**, and it is a *measure-average*, not a pointwise operator bound (PROOF-SAVVIDY-WGEO-DOMINATION: pushing Prékopa to a per-config Loewner bound is a category error).
- **Action `Hess S_W` indefinite at the Savvidy cone.** Exact transverse operator `H = −D² + 2ad F`, constant abelian background `F₁₂=B`: Landau ladder `λ_{n,s} = k_∥² + (2n+1)g|B| + 2s gB`, bottom `λ_min = −g|B|`. **Physical, not gauge:** every negative mode is transverse-magnetic with `δF = d_AδA ≠ 0`, survives the gauge quotient ([PROOF-SAVVIDY-EXACT-SIGN; PROOF-SAVVIDY-WGEO-DOMINATION]). Measured at β=2.5: `min(Ric+Hess S) ≈ −2.3` (SU2) / `−2.0` (SU3), **21–49% of directions negative**, 90–97% inside the high-momentum fiber/fluctuation block (W6-BE; MARGINAL).
- **The core is indefinite.** On the Savvidy-negative directions `Hess S ≤ −g|B| < −κ_H`, so `Γ₂ < 0` pointwise; against the positive `Ric = κ_H = O(1)` and `Hess S ∝ β`, this fails for every `m>0` once `β ≳ κ_H/1.1`. **Pointwise `(BE)` is FALSE at β=2.5** — this is W6, not new here; the new content is what follows.

### 2.2 Three candidate relaxations, and the exact verdict on each
The task's question — *measure-weighted or truncated `E_μ[Γ₂] ≥ ρ`* — has a sharp answer. Define
```
 (C_pt)   Γ₂(f)(x) ≥ ρΓ(f)(x)  ∀f,∀x                         (pointwise BE)
 (C_avg)  ∃ρ>0, ∀f mean-zero:   ∫Γ₂(f) dμ ≥ ρ ∫Γ(f) dμ        (measure-weighted)
 (C_tr)   ∃ρ>0:  ∫ Γ₂(f)·1_{Hess S≥0 or C_κ^c} dμ ( + Haar budget ) ≥ ρ ∫Γ dμ  (WSTC-truncated)
```

**(C_pt) — FALSE at β=2.5.** Established (W6; §2.1). No.

**(C_avg) — not a relaxation; it is the gap. (NEW, central.)** Two exact facts:
```
 (i)  ∫_M Γ₂(f) dμ  =  ‖Lf‖²_μ.                                             (2.2)
      [½∫LΓ dμ = 0 (∫Lg dμ=0); −∫Γ(f,Lf)dμ = ∫ f L²f dμ = ∫(Lf)²dμ, L self-adjoint.]
 (ii) inf_{f⊥1} ‖Lf‖²_μ / ‖∇f‖²_μ  =  λ₁  (the spectral gap of −L).          (2.3)
      [f = Σ a_k φ_k, −Lφ_k=λ_kφ_k:  ratio = Σλ_k²a_k / Σλ_k a_k = weighted avg
       of λ_k with weights λ_k a_k ≥ 0, hence ≥ λ₁, equality at f=φ₁.]
```
Therefore
```
   (C_avg)  ⟺  λ₁ ≥ ρ  ⟺  the spectral gap is ≥ ρ.
```
`(C_avg)` is **exactly as hard as the target.** It is **circular** — not a useful relaxation. Numerically checked on a random SPD `6×6` generator (λ₁=6.856; a random all-mode `f` gives ratio 15.52 ≥ λ₁; the pure lowest eigenmode gives 6.856=λ₁) and on 1D OU (`∫Γ₂dμ = ‖Lf‖² = 8` for `f=x²`, `=672` for `f=x⁴`, quadrature 80-pt Gauss–Hermite). **This is the species' sharp no-go: measure-averaged Bakry–Émery is the gap, not weaker than it.** The power of BE lies *entirely* in pointwise-ness (`(C_pt) ⟹ (C_avg) ⟹ λ₁≥ρ`); averaging destroys the only thing that made it a tool.

Along the coarse atom mode `a = ⟨v,X⟩` (linear, `Γ(a)=1`) the averaged quantity is the **positive identity** `E_μ[Γ₂(a)] = κ_H + E_μ[H] = κ_H + Var_μ(G)` (`G=⟨∇S,v⟩`, Stein `E[H]=E[G²]=Var G`). It is `>0` (≈`κ_H+0.2`), but it is the trivial Stein identity, carries **no** information about `Var(a)`, and it is exactly `0` on the free field (`G=0, κ_H^{spatial}=0`) — i.e. it is free-field-blind in the averaged sense and again confirms `(C_avg)` can't be the lever.

**(C_tr) — the one non-circular hope: WSTC truncation. It fails, and here is exactly where.** WSTC [ABC-NEW-INEQ] is a real, correctly-signed, volume-uniform *upper* bound on the cone variance:
```
   Var_{μ₀}(f·1_{C_κ}) ≤ osc_{C_κ}(f)²·μ₀(C_κ),   μ₀(C_κ) ≤ e^{−cβ(κ²−κ_typ²)N_□} → 0,
   C_κ = {‖F‖² ≥ κ²} ⊃ Savvidy cone {2g|F| > μ²},   uniform in L.              (2.4)
```
The idea: excise the high-field cone where `Hess S` is deepest-negative, keep the negative budget outside the truncated region. **It does not compose.**
- The starvation factor is `e^{−cβ(κ²−κ_typ²)N_□}`, i.e. `= 1` at `κ = κ_typ`. It suppresses only **high field**, `‖F‖ ≫ κ_typ`.
- The Savvidy bottom `λ_min = −g|B|` at longitudinal `k_∥ = 0` exists **for every background strength `B`**, including **typical** `|B| ≈ κ_typ`. That soft negative is *not* in `C_κ` for any `κ > κ_typ`, and for `κ = κ_typ` starvation is absent.
- W6 locates the actual negative variance mass at **typical field, high momentum** (90–97% in the fluctuation block), i.e. **outside** `C_κ` for every useful `κ`. WSTC bounds the *wrong tail* (ABC §4: the cone is a red herring; the sign is decided at typical action).
So excising `C_κ` removes a region carrying a vanishing fraction of the negative `Γ₂`, and **leaves untouched exactly the cell that decides the sign**. `(C_tr)` reduces to controlling `Γ₂` on `{typical field} ∩ {high momentum}` — which is the conditional variance bound, i.e. H1. **Answer to the task's sub-question: NO — WSTC starvation does not cover the soft `k_∥=0` negatives, because they sit at typical field where the exponent is 1.**

### 2.3 The honest reduction (the valuable output)
Put the three facts together via the exact marginal identity `G″(0) = Δ_dyn + W_geo + X`, `Δ_dyn = E[H]−Var(G)`, `X = −2Cov(G,Φ′)` [PROOF-JACOBIAN-SAVES]. The composition's `ρ`, on the mode where the sign is decided, is
```
   ρ  =  κ_H  +  ( E_μ[H|0] − Var_μ(G|0) )  +  W_geo  +  X.
```
For `ρ>0` one needs the **Savvidy-sign deficit** `E[H|0] − Var(G|0)` to be covered by the **Haar–Prékopa budget** `κ_H + W_geo` (plus `X`). Since `W_geo ≥ 0` is proven *as a measure-average* and the deficit is *pointwise/conditional*, the composition lands on
```
   ┌──────────────────────────────────────────────────────────────────────┐
   │   (◆_Γ₂)   E_μ[Hess S_W | 0]  ≥  ρ − W_geo          (≡ H1)            │
   │            ↔  Var_μ(⟨v,∇S⟩|a=0)  ≤  E_μ[H|a=0] − ρ + W_geo,           │
   │            uniform in (L,a),  at β=2.5.                               │
   └──────────────────────────────────────────────────────────────────────┘
```
This is **exactly the Savvidy-sign wall with `W_geo` as a budget.** The Haar convexity is a real credit (measured `W_full≈+0.22–0.31` beats the deficit `−Δ_dyn≈+0.02`), but certifying `(◆_Γ₂)` uniformly needs an **upper bound on `Var(G|0)`** — which *is* H1 in curvature language (`G″≈1/Var(a)` to 2–5%; `(◆_Γ₂) ⟺ Var(a)≤1/δ ⟺ m_eff²(k→0)≥δ`). The Bakry–Émery species does not add a new handle; it re-derives H1 with `W_geo` named as the credit.

---

## 3. Three-regime test (honest)

| regime | pointwise `Γ₂` | measure-avg `(C_avg)` | `(C_tr)` WSTC | acid test / atom | verdict |
|---|---|---|---|---|---|
| **Strong coupling `β<1/48`** | `Γ₂ ≥ κ_H − cβ > 0` (Ric dominates; `Hess S∝β`) | `(C_avg)⟺λ₁`; and `λ₁≥ρ(β)>0` genuine | not needed | `f₁² ≥ r₀(β)>0` | **PASS / consistent** |
| **Free massless `(β=0, S=0, flat)`** | `Γ₂(f)=\|Hess f\|² ≥ 0`, ratio `Γ₂/Γ = k² → 0` at `k→0` | red herring (gap=0) | cone empty `C_κ=∅` | returns `ρ=0`, `f₁²=0` | **PASS (non-vacuous)** |
| **Crossover `β≈2.5` (Savvidy)** | `min(Γ₂) ≈ −2.3 < 0`; 22–49% dirs neg, 90–97% fiber | `(C_avg)=`gap (circular) | negatives at typical field, `κ=κ_typ` ⇒ **no starvation** | `G″` sign = H1, uncertified | **FAIL (exact break)** |

**Strong coupling.** Pointwise `(BE)` holds: `Ric=κ_H=O(1)` (β-independent) vs `Hess S_W ≥ −cβ` (the plaquette Hessian is `β×` a bounded form, `|ad F|≤2‖F‖`, `‖F‖` bounded on compact links), so `ρ(β)=κ_H−cβ>0`, **uniform in `(L,a)`**. This yields a volume-uniform lattice spectral gap in a strong-coupling window — the same qualitative content as SZZ's rigorous LSI at `β<1/48` and P3's residue `f₁²≥r₀>0`; **no contradiction, mutually consistent.** Honest scope: the rigorous window is `β<β₀=κ_H/c` with `c` explicit but `β₀` modest; the measured `c≈1.12` is numerical, so BE does not by itself *extend* SZZ's `1/48` without a proven per-link `Hess S_W` bound. (BE here is a *worse-constant* reproduction of the known strong-coupling gap, not a new regime.)

**Free massless — acid test.** `S=0`, `μ=`Haar, flat chart `μ=`:Lebesgue on the torus, `L=Δ`, `Ric_flat=0`, `Γ₂(f)=|Hess f|²`. A purely spatial plane wave `f=e^{ikx}` has `|Hess f|²=k⁴`, `Γ(f)=k²`, so `Γ₂/Γ=k²→0` as `k→0`: the criterion returns **`ρ=0` exactly**, and `Var(a)=∞⇒G″=0⇒f₁²=0`. **Not vacuous** — the candidate goes to `0`, it does not manufacture a spurious positive. (The β=0 *Haar* group-fiber gap `κ_H` is a different, real lattice gap — the local `Ad` zero mode — not the physical `k→0` mode; on the physical spatial `k→0` direction `Ric` contributes `0`.)

**Crossover — exact break.** Pointwise `Γ₂` is negative on the physical Savvidy directions (W6). The measure-weighted escape is circular (§2.2). The truncated escape is defeated by location: whose question — *does WSTC cover the soft `k_∥=0` negatives?* — has answer **NO**, because the soft bottom `−g|B|` at `k_∥=0` lives at **typical** `|B|≈κ_typ` where `μ₀(C_κ)` starvation `= e^{0}=1`; WSTC suppresses only the measure-negligible high-field tail, while the sign-carrying negatives are at typical field, high momentum (W6: `−1.12β`, 90–97% fluctuation block). The residue is `Var(G|0)` = H1. This is where the composition first stops in the three-regime test (and, as β is raised from 0, the *pointwise* form first stops at `β≈κ_H/1.12≈0.45` (SU2), long before 2.5).

---

## 4. Scorecard (0–5)

| # | condition | verdict | why |
|---|---|---|---|
| 1 | **nonperturbative** | **PASS** | generator, Prékopa `W_geo≥0`, Ward–Stein, WSTC are exact/all-β; no weak-coupling expansion |
| 2 | **sign-controlling LOWER bound** | **PARTIAL** | `(BE)` gives `m_eff²≥ρ` (correct direction) — but only pointwise; the pointwise form is false at 2.5 and the averaged form is circular |
| 3 | **non-circular** | **FAIL (at the decisive form)** | `(C_avg)` *is* `λ₁≥ρ` (identity `∫Γ₂=‖Lf‖²`, `inf=λ₁`); `(C_tr)` non-circular but bounds the wrong region |
| 4 | **non-SUSY + fails free massless** | **PASS** | pure YM, real indefinite β=2.5; free field returns `ρ=0`, `f₁²=0` (acid test met) |
| 5 | **`(L,a)`-uniform `m_eff²(k→0)≥λ`** | **FAIL** | pointwise `ρ>0` only for `β<β₀`; at 2.5 no `L`-uniform `ρ`, WSTC is `L`-uniform but maps to H1 |

**Score 2.5 / 5.** Right-signed, nonperturbative, acid-test-clean *as a criterion* — but it does not close the gap; it reduces to it. The 0.5 credit is for condition (2): the species alone among the ~60 routes points in the gap-increasing direction, and honestly names its own failure.

---

## 5. Sharpest line, exact stuck step, refs

**Sharpest line.** *Bakry–Émery is the one correctly-signed tool, but its only usable form is pointwise; its measure-weighted form is not a relaxation of the gap but an identity for it (`∫Γ₂dμ=‖Lf‖²`, `inf_f=λ₁`). The pointwise form dies at the Savvidy bottom `−g|B|` (physical, 22–49% of directions), and it cannot be rescued: the Haar–Prékopa budget `W_geo≥0` is a measure-average against a pointwise negativity (category mismatch), and WSTC starvation does not reach the soft `k_∥=0` negatives because they sit at typical field where its exponent is exactly `1`. So the composition lands, word for word, on `G″ = E[H]−Var(G)+W_geo+X ≥ ρ` — `E[Hess_S] ≥ ρ−W_geo`, the Savvidy-sign wall with `W_geo` as a real but insufficient credit. That is H1.*

**Exact stuck step (one inequality).**
```
   ∃ρ>0, ∀f mean-zero, uniform in (L,a):
        ∫_M Γ₂(f) dμ  =  ‖L f‖²_μ  ≥  ρ ‖∇f‖²_μ = ρ∫Γ(f)dμ.
```
By (2.2)–(2.3) this is `λ₁≥ρ` — the gap. The *pointwise sufficient* version `Γ₂(f)(x)≥ρΓ(f)(x)` fails along the typical-field transverse-magnetic Savvidy negatives (where `Hess S≤−g|B|<−κ_H`, `δF≠0`, physical); the *measure-averaged* version is identical to `λ₁` (circular); the *WSTC-truncated* version leaves those negatives untouched. The single missing estimate is the **conditional variance bound**
```
   Var_μ(⟨v,∇S⟩ | a=0)  ≤  E_μ[H|a=0] − ρ + W_geo,   uniform (L,a),  at β=2.5   (= (◆_Γ₂) = H1),
```
equivalently `Var_ν(φ) ≤ (1/λ)∫|∇φ|²dν` with `λ=ρ` on the coarse mode. Nothing in the transport species supplies it; that is the wall.

**Novel content of this build.** (i) The exact identity `∫Γ₂dμ = ‖Lf‖²` and the spectral-min `inf_f ‖Lf‖²/‖∇f‖² = λ₁`, hence the **no-go**: *measure-averaged Bakry–Émery ⟺ the spectral gap* (removes the naive "weighted BE is weaker" hope, in-species). (ii) The exact composition statement `E[Hess_S] ≥ ρ − W_geo` identifying `W_geo` as the Haar budget. (iii) The precise disjunction "starvation region ≠ negative-`Γ₂` region," pinning why WSTC truncation cannot rescue the pointwise form. All checked, none overclaimed: no QED.

**Refs.** Bakry–Émery, *Diffusions hypercontractives*, LNM **1123** (1985); Bakry–Ledoux–Gentil *Markov Semigroups*; Otto, Comm. PDE **26** (2001) 101; Sturm, Acta Math. **196** (2006) 65; Lott–Villani, Ann. Math. **169** (2009) 903; von Renesse–Sturm, CPAM **58** (2005) 923; Holley–Stroock, Ann. Probab. **15** (1987) 1151; Cattiaux–Guillin–Wu, Ann. Probab. **42** (2014) 2161; Prékopa, Acta Sci. Math. **32** (1971) 301; Brascamp–Lieb, J. Funct. Anal. **22** (1976) 366; Savvidy, Phys. Lett. B **71** (1977) 133; Nielsen–Olesen, Nucl. Phys. B **144** (1978) 376; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737.
**Internal.** W6-BE-at-2.5 (pointwise `(BE)` false: `min −1.12β`, 22–49%); PROOF-JACOBIAN-SAVES (`G″=E[H]−Var+W_geo+X`, `W_geo≥0` Prékopa, `(◆)`); ABC-NEW-INEQ (WSTC `(2.4)`, "cone is a red herring"); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2ad F`, `λ_min=−g|B|`, physical); PROOF-SAVVIDY-WGEO-DOMINATION (per-config domination false; `W_geo` is measure-level); PROOF-COARSE-POINCARE (H1, `β<1/48` anchor); SCORE-HYPOTHESIS (`λ` = H1's coarse constant); P3-wasserstein-ym; new-tools-registry A5 (gauge `W`-entropy, same `Γ₂`, "correct gapless-side / blocked gapped-side"); INVENT-NEW-MONOTONE.
