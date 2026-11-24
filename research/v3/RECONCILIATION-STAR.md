# RECONCILIATION-STAR — which `G″` is the object of (★), what G1/G3 actually computed, and the resulting status of (★)

> 2026-11-24 · Adjudication of the definitional conflict **G1/G3 (`G″=2W_geo`) vs G4 (`G″=E[Σ]−D`, `2W_geo` invalid)**.
> Reads: `docs/paper-v2/ch05-theorems.tex` (`eq:split`, `eq:vardecomp`, `eq:wgeo`, `sec:wgeo`);
> `research/v3/C2-BALABAN-UNIFORM.md` §2.4; `research/v3/S2-BOUNDARY-COUNTEREXAMPLE.md`,
> `G1-BLOCK-IDENTITY.md`, `G3-PROJECTION-TERM.md`, `G4-GENERAL-BOUND.md`, `S3-UV-FESHBACH.md`;
> `research/unconditional/PROOF-JACOBIAN-SAVES.md`, `INVENT-VARIANCE-INEQ.md`.
> Witness: `research/v3/w7/r1_model_curvature.py` (`py -3.12`; deterministic Gauss–Legendre, output quoted).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[V]** numerical (deterministic) ·
> **[NG]** no-go / disproof for a named reading · **[O]** exact missing step. **Clay 0%. No QED.**

---

## 0. Verdict (one paragraph)

**The right object for (★) is the coarse-field-conditioned marginal curvature (C):**
`G″_j(v;a₀)=E_{μ(·|U′)}[∂_v²V|a=a₀]−Var_{μ(·|U′)}(∂_vV|a=a₀)`, with `V=S_W+Φ`, conditioning on the
coarse coordinate `a=⟨PᵀX,v⟩=a₀`. **It includes `X`** (via `W_j=W_geo+X`) **and** the conditioning `|a=a₀`.
It equals `E[Σ]−D` (S3). **G1/G3 computed a different object** — the *unconditioned, X-less*
`G″_0:=E_μ[H]−Var_μ(G)+W_geo,μ`, which by the unconditional Stein identity is exactly `2W_geo,μ`. That is a
true algebraic identity (`[T]`) but **not the `G″` of (★)**. Consequently: (i) **the S2/G1/G3 reduction
`(★)⇔2W_geo≥δ_W` is invalid `[NG]`** (G4 is right); (ii) reading C2 literally with `W_j=W_geo+X` and no
conditioning gives `G″≡0` (vacuous, `[NG]`); (iii) **the corrected (★) at the RG-relevant coarse field is
`[O]`** — it reduces to a conditioned Schur/coarse-gap bound, not to `W_geo`; (iv) **the "uniform in every
coarse field `a₀`" reading is `[NG]`** — G4's barrier counterexample displaces it (at `a₀=0`, `B→−I`,
longitudinal `e₃`, `β>β*≈4.31`), but that point is a suppressed *barrier*, not the effective-action minimum.

---

## 1. The definition: does `G″` for (★) include `X`? — **Yes.** (exact quotes + identity)

### 1.1 The campaign's exact split (`ch05`, eq. `eq:split`) — the definition of `G″(0)`

```
G''(0) = E[ H_vv + Φ'' | ξ=0 ] − Var( G_v + Φ' | ξ=0 ) = Δ_dyn + W_geo + X ,
Δ_dyn = E[H] − Var(G),          X = −2 Cov(G, Φ') .
```
Every `E`, `Var` carries the conditioning `| ξ=0` (the coarse field pinned to `0`; `ch05` eq. `eq:vardecomp`
writes the same as `G″(0)=S(v)−D(v)`, `S=E[h_n−⟨K⁻¹c,c⟩|0]`, `D=Var(∂_nV|0)−E[⟨K⁻¹c,c⟩|0]`). This is the
marginal curvature `−d²/da² log ρ_v(a)|_{a₀}` of the coarse coordinate (`INVENT-VARIANCE-INEQ` §1, S3 Lemma S,
G4 Lemma G1: `G″(v;a₀)=E[∂_v²V|a₀]−Var(∂_vV|a₀)`, pure differentiation under the integral, no IBP).

### 1.2 The C2 statement, with the load-bearing `W_j` (§2.4)

```
(★)  Var_{μ(·|U′)}(⟨Pᵀ∇S_j, v⟩) ≤ E_{μ(·|U′)}⟨PᵀH_jP v, v⟩ + W_j(v) − δ_W   (δ_W>0),
     equivalently  G″_j(v) = E[H_j] − Var_j + W_j ≥ δ_W .
```
**`W_j` means `W_j = W_geo + X`, not `W_geo`.** S3 fixes the convention explicitly (§0.3, §2.3):
> "In the S-split, `W_j = W_geo + X`, `W_geo = E[Φ″] − Var(Φ′)`, `X = −2Cov(G,Φ′)`, `G = ⟨Pᵀ∇S,v⟩`."
> "`W_j := W_geo + X = E[Φ″] − Var(Φ′) − 2Cov(G,Φ′)`. … `X` is indefinite. `W_j` **cannot absorb `X`**
> without a bound on the force–Jacobian covariance."

`PROOF-JACOBIAN-SAVES` §1 agrees ("W7's `W` is `W_full = W_geo+X`"). **So the C2 inequality, read at face
value, is the object (C) with `X` included** — *provided* the moments are the conditioned ones and `W_j` is the
full `W_geo+X`. The C2 shorthand suppresses two load-bearing facts: the conditioning `|a=a₀` and the `X` in `W_j`.

### 1.3 What G1/G3 actually computed — call it `G″_0`

G1 §1 and G3 §1 define (identical):
```
G_v := ∂_v S,  H_vv := ⟨v, Hess S·v⟩,  φ := ∂_v Φ,
W_geo(v) := E[Φ″_vv] − Var(φ),   G″_0(v) := E[H_vv] − Var(G_v) + W_geo(v) ,
```
with `E`, `Var` over the **full fibre measure `μ(·|U′)`** — **no coarse-field conditioning, and `W_geo` only
(no `X`).** The Stein identities `(S0)–(S3)` (G1 §3) are then **unconditioned**, and give (G1 Theorem, G4 Lemma G2)
```
Δ_dyn,μ := E_μ[H]−Var_μ(G) = E_μ[Φ″]−Var_μ(Φ′) = Cov_μ(G,Φ′) = W_geo,μ ,
so   G″_0 = Δ_dyn,μ + W_geo,μ = 2 W_geo,μ ,      and   X_μ = −2W_geo,μ  ⇒  (literal full sum) = 0 .
```

### 1.4 The difference `G″ − G″_0`

Write `a` for the coarse coordinate and subscript `a`/`μ` for conditioned/unconditioned moments.

| term | G1/G3/S2 object `G″_0` (U) | true (★) object `G″` (C) |
|---|---|---|
| definition | `E_μ[H]−Var_μ(G)+W_geo,μ` | `E_a[H+Φ″]−Var_a(G+Φ′) = Δ_dyn,a+W_geo,a+X_a` |
| includes `X`? | **no** | **yes** (`W_j=W_geo+X`) |
| coarse-conditioned? | **no** | **yes** (`|a=a₀`) |
| equals | `2W_geo,μ` (Stein, `[T]`) | `E_a[Σ]−D` (S3, `[T]`) |

The exact difference is
```
G″ − G″_0 = X_a  +  (Δ_dyn,a − Δ_dyn,μ)  +  (W_geo,a − W_geo,μ) .
```
- Within a **fixed** measure, `G″ − (E[H]−Var(G)+W_geo) = X` exactly (the `ch05` statement). G1/G3 dropped `X`.
- But G1/G3 *also* used the unconditioned measure; the full gap between (C) and `2W_geo,μ` is `X_a` **plus**
  the conditioning shifts of `Δ_dyn` and `W_geo`. **The conditioning is the load-bearing piece, not `X` alone**
  (see §2: on G4's own counterexample slice `X_a=0`, and the negativity is entirely the conditioning shift of
  `Δ_dyn`). This is a sharper, honest correction to G4's "it omits the cross term `X`": the cross term is
  omitted too, but at the specific counterexample it vanishes.

---

## 2. Exact recomputation of the one-plaquette S2 model for the true `G″` (`r1_model_curvature.py`)

Model (S2, unchanged): SU(2) one plaquette, one free link `U=e^ξ`, `ξ∈ℝ³`, `ρ∈[0,2π]`, boundary
`B=e^{iθσ₃/2}`, `S=β(1−½tr UB)`, `Φ=−log J`, `dμ∝e^{−S−Φ}d³ξ`. Coarse mode `v=e₃`, condition on `a=ξ₃=a₀`.

### 2.1 The exact slice decomposition at `θ=π`, `v=e₃`, `a₀=0` — **[T]**

At `θ=π` and on the slice `ξ₃=0`: `S(r,0)=β` (a **constant** in `r`), so `e^{−S}` factors out and the slice
measure is **β-independent Haar** `∝Φ`-weight. Directly (radial `Φ`): `Φ′_v=0`, `Φ″_vv=Φ′(r)/r`,
`G_v=β sin(r/2)/r`, `H_vv=0`. Hence **exactly**
```
G″(e₃; π, 0) = E_Haar[Φ′(r)/r] − β² Var_Haar[ sin(r/2)/r ]
             =   W_geo,a          +   Δ_dyn,a        ( X_a = 0 ) .
```
So on this slice **`X_a=0`** (radial symmetry kills `Cov(G,Φ′)`), **`W_geo,a = E_Haar[Φ′/r] > 0` (β-independent)**,
and the **only** negative ingredient is `Δ_dyn,a = −β² Var_Haar[sin(r/2)/r]` — the Wilson **force variance**
(the Savvidy double well in the marginal). Witness (`r1`, Gauss–Legendre, converged to 8 digits):

| `β` | `G″` (exact slice) | `G″` (FD of `−d²/da²logρ₃`) | `Δ_dyn,a` | `W_geo,a` | `X_a` | `2W_geo,μ` (G1/G3) |
|---|---|---|---|---|---|---|
| 0 | **0.205115** | 0.205115 | −0.000000 | 0.205115 | 0 | **0.000000** (6.97e−8) |
| 1.0 | 0.194056 | 0.194056 | −0.011059 | 0.205115 | 0 | 0.06390 |
| 2.5 | 0.135995 | 0.135995 | −0.069121 | 0.205115 | 0 | 0.26182 |
| 4.0 | 0.028167 | 0.028167 | −0.176949 | 0.205115 | 0 | 0.39792 |
| 5.0 | **−0.071367** | −0.071367 | −0.276482 | 0.205115 | 0 | 0.44775 |
| 10.0 | −0.900814 | −0.900814 | −1.105929 | 0.205115 | 0 | 0.53156 |
| 20.0 | −4.218602 | −4.218602 | −4.423718 | 0.205115 | 0 | 0.56462 |

**Exact β-law:** `G″(e₃;π,0) = A − Bβ²`, with `A = E_Haar[Φ′/r] = 0.205115` and
`B = Var_Haar[sin(r/2)/r] = 0.0110593` (both converged and `n`-independent for `n=800…12800`);
zero crossing `β* = √(A/B) = 4.3066`. *(This is a small numerical correction to G4's quoted
`0.2057−0.0123β²`, `β*=4.09`: G4's own finite-difference column already agrees with the converged
values here; only its analytic-constant column was slightly off. The qualitative verdict — `G″` turns
negative for `β` above `≈4.3` — is unchanged.)*

### 2.2 The β=0 check against the campaign theorem `G″ ≥ N/12` — **consistent, and resolves S2's `W_geo=0`**

- `N=2`: `N/12 = 0.166667`. True (C) at `β=0`: **`G″ = 0.205115 ≥ N/12`** ✓ (`Δ_dyn,a=0`, `W_geo,a=0.2051`, `X_a=0`).
  This is the conditioned `Σ`-form of the `β=0` theorem (`V=Φ` convex ⇒ `K⪰c_Φ`, `D≤0`, `G″≥E[Σ]≥c_Φ=N/12`;
  S3 Cor. S, `PROOF-JACOBIAN-SAVES` Lemma B).
- S2's "`W_geo=0` at `β=0`" is the **unconditioned** `W_geo,μ = E_μ[Φ″]−Var_μ(Φ′) = 0` (by `ξ→−ξ` symmetry
  `E_μ[Φ′]=0` and Stein `Var_μ(Φ′)=E_μ[Φ″]`); its unconditioned `G″_0=2W_geo,μ≈7e−8`.
- **Resolution.** Both are correct and different objects. `W_geo,μ(β=0)=0` (unconditioned) is consistent with
  the `β=0` theorem because the theorem is about the **conditioned** curvature, whose value is the
  **conditioned** `W_geo,a = E_Haar[Φ′/r] = 0.2051 > N/12`, **not** the unconditioned `W_geo,μ=0`. S2 is
  internally consistent but was computing the (U) object, not (C).

**Does the single-link model differ from the lattice coarse mode?** Not in the mechanism: both give the
`β=0` conditioned curvature `= E_a[Φ″|a₀]−Var_a(Φ'|a₀) ≥ N/12` by the same Prékopa/`Φ″⪰N/12` argument. It
differs in one quantitative respect that matters at `β>0`: the 0-dimensional toy has **no spatial free gauge
gap `c_L=4π²/L²`** on the coarse mode and **no gapped transverse fibre** `K⪰γ_j`. The real lattice UV block has
both, so `E[Σ]≥c_L+c_Φ−`Schur-loss `>0` at the RG-relevant coarse field (S3 §2.2); the toy's only floor is
`c_Φ≈0.205`, which `Δ_dyn,a=−β²Var` overwhelms at large `β` (this is the origin of the `β>4.3` negativity).

### 2.3 The unconditioned identity (cross-check, `[T]`+`[V]`)

`Δ_dyn,μ=W_geo,μ`, `X_μ=−2W_geo,μ`, `Δ_dyn,μ+W_geo,μ+X_μ=0` — verified to `≤4.5×10⁻⁷` at all 9 `(β,θ)`
(G4 Lemma G2). `2W_geo,μ>0` everywhere for `β>0` (`β=2.5: 0.262`, `β=5: 0.448`) — yet the true (C) turns
negative at the barrier. This is the quantitive heart of the conflict.

---

## 3. Scope of G4's counterexample — is it inside (★)?

G4's counterexample: one plaquette, boundary `B→−I` (`θ=π`), **longitudinal** `e₃`, coarse field `a₀=0`,
`β>β*≈4.31`: `G″=A−Bβ²<0` (`β=5: −0.0714`). Three scope questions:

1. **Is the model a legitimate block?** Yes — the smallest legitimate block (one plaquette, one free link,
   fixed boundary). It is exactly S2's model, and its `2W_geo` reproduces S2. `[T]` as a computation of (C).
2. **Is `e₃` admitted by the projection `P`?** As a **direction**, yes: `e₃` is a tangent direction of the
   fibre, admissible for a constant-coefficient `P` (S2/G3 Reading 1) and representable under a covariant `P`.
   But G4's object further **conditions on the value `a=⟨e₃,ξ⟩=a₀`** — an **interior coarse-field (block-spin)
   constraint**, which is reading (C) (`ch05`'s `|0`), *not* C2's literal boundary-only `μ(·|U′)`. That
   interior constraint is precisely the one door S2 §6 and G1/G3 §5(c) left `[O]`.
3. **Is the barrier configuration in the support at `B→−I`?** Yes — `ρ₃(0|θ=π)>0` (strictly positive marginal
   density on the whole disk), so `(B→−I, a₀=0)` is in `supp μ(·|B)`. **But it is the barrier, exponentially
   suppressed relative to the effective-action minimum.** Witness (`r1`):

| `β` | argmax `a` of `ρ₃` (well) | `G″` at well | `ρ₃(0)/ρ₃(peak)` | `G″` at `a₀=0` | `G″` for `a≳0` (far barrier) |
|---|---|---|---|---|---|
| 5.0, θ=π | −2.55 | **+1.44** (`a=−3`) / +1.11 (`a=−2`) | **0.045** | −0.071 | −0.40 (`a=+1`), −0.34 (`a=+2`) |
| 2.5, θ=π | −2.10 | +0.66 (`a=−2`) | 0.31 | **+0.136** | −0.037 (`a=+1`), −0.077 (`a=+2`) |

So `a₀=0` is the **top of the double-well barrier**, carrying `≈4.5%` (β=5) of the peak density; the
**effective-action minimum** (the mode of the coarse marginal = the RG-relevant coarse spin) is at
`a≈−2.5`, where **`G″>0`** (and `≫N/12`). The far side of the barrier (`a>0`) is even more suppressed and also
negative.

**Precise scope statement.** G4's counterexample is a **legitimate, in-support computation of the true (C)**
object. It therefore **disproves** (★) under the reading:
> **(★_unif-a₀)** `G″(v;a₀) ≥ δ_W` for **every** coarse field `a₀∈supp ρ_v(·|U′)`, uniformly in `v,U′`.

That is the *global convexity* of the coarse effective action `S_eff(a)=−log ρ_v(a)` — **false** `[NG]`
(the double well at strong background; consistent with S1 `Prop U′-NG`). It does **not** disprove (★) under
the reading the campaign needs:
> **(★_RG)** `G″(v;a₀(U′)) ≥ δ_W` at the **RG-relevant** coarse field `a₀(U′)` = the effective-action minimum
> (mode of the coarse marginal / the coarse spin actually carried into the next scale).

At the barrier the toy also lacks the spatial `c_L` that dominates the real lattice UV block (§2.2), so the
barrier negativity does **not** transfer to the lattice RG-relevant point. **Verdict of §3: (★_unif-a₀) `[NG]`;
(★_RG) untouched by this counterexample and `[O]` (§4).**

---

## 4. Corrected statement of (★) and verdict

### 4.1 Corrected (★) — all terms carried (conditioning, `W_j`, `X`, transport defect) — **[definition, T]**

For shell `j`, environment `U′`, cutoff `a`, coarse mode `v=Pw`, at the RG-relevant coarse field `a₀(U′)`:
```
(★)   G″_j(v; a₀) := E_{μ(·|U′)}[ ⟨Pᵀ Hess V_j P v, v⟩ | a=a₀ ]
                   − Var_{μ(·|U′)}( ⟨Pᵀ ∇ V_j, v⟩ | a=a₀ )   ≥  δ_W > 0 ,
      V_j = S_j + Φ ,   a := ⟨Pᵀ X, v⟩ ,   a₀ = a₀(U′) (effective-action minimum) ,
```
with the **exact** decomposition (all moments conditional on `a=a₀`; identities `[T]`):
```
G″_j = Δ_dyn + W_j = E[Σ] − D ,
Δ_dyn = E[H]−Var(G) ,   W_j = W_geo + X ,   W_geo = E[Φ″]−Var(Φ′) ,   X = −2Cov(G,Φ′) ,
Σ = vᵀPᵀ Hess V_j Pv − ⟨K⁻¹c,c⟩ ,   D = Var(∂_vV) − E[⟨K⁻¹c,c⟩] ,
```
and for a **covariant** projection `P` the transport (divergence / anholonomy) defect is carried explicitly
(G3 §4 / G1 §4): in the unconditioned picture `G″_P = 2W_geo + R_proj`, with
```
R_proj = E[⟨(Dw)w,∇Φ⟩ − ⟨(Dw)w,∇S⟩] + Cov(⟨w,∇Φ⟩, div w) − Cov(⟨w,∇S⟩, div w) ,
```
which is `O(g_j)`-small in the UV tower `[C]` and `O(W_geo)` sign-indefinite at the head `[NG] for "controlled
by W_geo"` (G3 §5–6). This `R_proj` is the "covariant divergence" term; in the conditioned `Σ/D` picture it is
absorbed into `Σ` (via `Φ″`, `Π₀ Hess Φ`) and `D` (via `Var(Φ′)`, `Cov`) — S3 §2.3.

### 4.2 Verdict ledger (no QED)

| statement | label |
|---|---|
| Correct object for (★) = conditioned marginal curvature `G″(v;a₀)=E_a[∂_v²V]−Var_a(∂_vV)=E[Σ]−D`; includes `X` and `|a=a₀` | **[T]** (identity: `ch05` `eq:split`/`eq:vardecomp`, S3 Lemma S, G4 Lemma G1) |
| `W_j = W_geo + X` (not `W_geo`); `W_j` cannot absorb `X` | **[T]** (S3 §0.3, §2.3) |
| G1/G3 computed `G″_0=E_μ[H]−Var_μ(G)+W_geo,μ = 2W_geo,μ` (unconditioned, no `X`) — a true identity but **not** (★) | **[T]** + **[V]** (`r1`, §2.3) |
| S2/G1/G3 reduction `(★)⇔2W_geo≥δ_W` as a certificate of (C) | **[NG]** (compares the wrong object; G4 §0.2/§1.2 correct) |
| Literal C2 inequality `E[H]−Var(G)+W_j`, `W_j=W_geo+X`, no conditioning (`G″≡0`) | **[NG]** (vacuous `0≥δ_W`; G4 Lemma G2) |
| `β=0`: true (C) `=0.2051 ≥ N/12=0.1667` while unconditioned `W_geo,μ=0` — S2 consistent, reduction mis-predicts | **[T]**+**[V]** (`r1`, §2.2) |
| (★_unif-a₀): `G″≥δ_W` for every `a₀∈supp ρ(·|U′)`, uniformly in `U′` (global convexity) | **[NG]** (G4 counterexample: `e₃`, `B→−I`, `a₀=0`, `β>4.31`; independently `S1 Prop U′-NG`) |
| **Corrected (★_RG)** at the effective-action minimum, uniform in `(j,U′,a)` | **[O]** — reduces to `(H-C) ∧ (H-SL) ∧ (★dec)`: fibre gap `K⪰γ_j I`; conditioned Schur bound `E_{μ(·|U′)}[Σ|a₀]≥δ_W`; head clustering. **Named missing lemma: a uniform lower bound on the conditioned expected Schur complement `E_{μ(·|U′)}[Σ|a₀(U′)] ≥ δ_W` at the RG-relevant coarse field** (`= ` the coarse spectral gap H1 at `a₀(U′)`, plus the `O(g_j)` transport-defect bound for covariant `P`) |
| `W_geo>0` at `β>0` is the gate | **no** — `(★)≠2W_geo`; the gate is `E[Σ|a₀]≥δ_W` = coarse gap H1 (**[T]** that it is not the gate; **[O]** the gate) |
| free/flat field acid test (`Φ≡0 ⇒ G″=0`) and `β=0` floor | **[T]**+**[V]** |
| Clay mass gap | still open, **0%** |

### 4.3 Bottom line

- **Which object is right:** the **coarse-field-conditioned marginal curvature (C)**, `G″ = Δ_dyn + W_j =
  E[Σ]−D`, with `W_j=W_geo+X` and the conditioning `|a=a₀` — i.e. `G″` **does include `X`**, and (more
  importantly) **must be conditioned on the coarse field**. G1/G3's `G″_0 = 2W_geo,μ` is a true identity for
  the *unconditioned, X-less* object and is **not** the `G″` of (★).
- **Status of (★):** the candidate reduction `(★)⇔2W_geo≥δ_W` is **disproved `[NG]`** (G4 correct). The
  literal no-conditioning C2 form is `[NG]` (identically zero). The **uniform-in-every-coarse-field** reading
  is `[NG]` (G4's barrier counterexample, `β>4.31`). The **corrected (★)** at the **RG-relevant** coarse field
  (effective-action minimum) is **`[O]`**: it is exactly the conditioned Schur/coarse-gap bound
  `E[Σ|a₀(U′)]≥δ_W` — the coarse mass gap H1 — and is **not** upgraded to a theorem. **No QED. Clay 0%.**

---

## 5. Verification ledger

| item | check | status |
|---|---|---|
| `eq:split` conditioned reading, `W_j=W_geo+X` | read from `ch05-theorems.tex` §`sec:wgeo`; S3 §0.3/§2.3 | **[T]** (text) |
| G1/G3 object = unconditioned `2W_geo,μ` | G1 §1/§3, G3 §1/§2; G4 Lemma G2; `r1` §(C) residual `≤4.5e-7` | **[T]**+**[V]** |
| exact slice decomposition `G″=Δ_dyn,a+W_geo,a+X_a`, `X_a=0` | `r1` slice formula = FD of `−d²/da²logρ₃` to 6 digits, all β | **[T]**+**[V]** |
| `G″=A−Bβ²`, `A=0.205115`, `B=0.0110593`, `β*=4.3066` | converged `n=800…12800`; matches G4 FD column | **[V]** |
| `β=0`: `G″=0.2051≥N/12=0.1667`, `W_geo,μ≈0` | `r1` §(A) | **[T]**+**[V]** |
| barrier vs well: `ρ₃(0)/ρ₃(peak)=0.045` (β=5), `G″`(well)>0 | `r1` marginal scan | **[V]** |
| covariant transport defect `R_proj` formula | G3 §4 (residual `≤1.1e-7`) | **[T]**+**[V]** (carried) |

**Corrections to the record (honest).** (i) G4's analytic constant `0.2057−0.0123β²` (β*=4.09) is slightly
off; the converged value is `0.205115−0.0110593β²` (β*=4.31), and **G4's own FD column already agrees with
this** — the sign/verdict is unaffected. (ii) G4's framing "the missing term is the cross term `X`" is
incomplete for the specific counterexample: on the `θ=π`, `e₃`, `a₀=0` slice **`X_a=0`** and the negativity is
the *conditioning shift of `Δ_dyn`* (the Wilson force variance `−β²Var`). Both `X` and the conditioning are
missing from G1/G3; the conditioning is the load-bearing one.

**Witness.** `research/v3/w7/r1_model_curvature.py` (`py -3.12`; deterministic Gauss–Legendre, output quoted
in §2–§3). No third-party sources downloaded.

*Refs.* `ch05-theorems.tex` `eq:split`, `eq:vardecomp`, `eq:wgeo`, `sec:wgeo`; `C2-BALABAN-UNIFORM.md` §2.4;
`S2-BOUNDARY-COUNTEREXAMPLE.md` §1–2, §6; `G1-BLOCK-IDENTITY.md` §1–4; `G3-PROJECTION-TERM.md` §1–8;
`G4-GENERAL-BOUND.md` §1–5; `S3-UV-FESHBACH.md` §0–2; `PROOF-JACOBIAN-SAVES.md` §1–2;
`INVENT-VARIANCE-INEQ.md` §1–3; `S1-CONDITIONAL-POINCARÉ.md` (`Prop U′-NG`).
