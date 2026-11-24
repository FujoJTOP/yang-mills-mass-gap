# V1-UV-NORMALIZATION — the UV-end normalization of (★): the free-block limit, the δ_W scaling, and the renormalized statement L3 actually needs

> 2026-11-24 · Attack angle V1. **The last open inequality (★) of the campaign, read at the *UV end*
> of the RG tower** (`β_j → ∞`, weak coupling), where L3 lives — as opposed to the IR end (`β≈2.5`),
> where the gap / H1 lives. The task: settle whether the **bare** inequality survives the UV limit,
> compute the free block exactly, derive the scaling of δ_W, and write the **renormalized** `(★_ren)`
> that L3 needs with `Z_j`.
> Reads: `docs/paper-v2/ch05-theorems.tex` (`eq:wgeo`, `eq:split`, `eq:vardecomp`),
> `ch03-existence.tex` (§L3, `eq:l3`, §dim), `C2-BALABAN-UNIFORM.md` §2.4, `S5-N-UNIFORM.md`,
> `G1-BLOCK-IDENTITY.md`, `G4-GENERAL-BOUND.md`, `S3-UV-FESHBACH.md` (`G″=E[Σ]−D`),
> `G1-L3.md`, `G1-UNIFORM-ZA.md`, `G1-L3-TOWER-POINCARE.md`, `G1-L3-UV-FORCE-BOUND.md`,
> `G3-PROJECTION-TERM.md` (the `g_j`-scaling table).
> Witness: `research/v3/w7/v1_uv_normalization.py` (`py -3.12`, deterministic, output quoted §6).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional / proof-sketch · **[V]** numerical
> (deterministic; residual quoted) · **[NG]** no-go for a named class · **[O]** exact missing step.
> No QED. Clay 0%.

---

## 0 · Verdict (no QED)

**A constant bare `δ_W` does *not* survive as a renormalized, `j`-uniform, positive-margin floor.** The
bare objects `H_j, E[H], Var` all scale like `β_j = 2N/g_j² → ∞`; the physically correct coarse
curvature `G″_j` scales like `β_j Q_c` with it; a *constant* `δ_W` is therefore a fraction
`O(1/β_j)=O(g_j²)` of the scale it sits under, and its renormalized value `Z_j^{-1}δ_W → 0`. The bare
floor must track the renormalized Hessian:

```
   δ_W  =  Θ(Z_j)  =  Θ(β_j)  =  Θ(1/g_j²)  ∝  log(1/a_j)      (asymptotic freedom → LOG, not a power of a_j),
```

so that the **dimensionless, j-uniform floor** is `δ̂_W := Z_j^{-1}δ_W = Θ(1) > 0`, and the correct
renormalized statement is

```
 (★_ren)   Z_j^{-1} G″_j(v)  ≥  δ̂_W > 0 ,     G″_j := −d²/da² log ρ^{(j)}_v(a)|_{a₀} = E_{μ(·|U′,a=a₀)}[Σ_j] − D_j ,
```

equivalently (coarse-pinned variance form) `Z_j^{-1}[Var − E[H] − W_j] ≤ −δ̂_W` at `a=a₀`; here `G″_j` is
the **coarse-pinned marginal curvature**
(G4's (C) = S3's `E[Σ]−D`), not the boundary-conditioned `E[H]−Var+W`. The `j`-uniform floor is the
**free-box gap** `δ̂_W = Q_c = c_L = 4π²/L²` **[T for the free block]** (interacting corrections are the
open `[O]`).

Three finite-size bits of bookkeeping, all verified exactly (§2, §6):

1. **Free block, exact.** In the quadratic-Wilson + Haar block, `E[H] = Var(G) = β_j Q_c + o(β_j)` (the
   `O(β_j)` parts cancel — Gaussian equipartition), the residual `E[H]−Var = W_geo = m β_jQ_c/(β_jQ_c+m)`
   `→ m = N/12` (coupling-free Haar floor, `j`-uniform in lattice units). **`G″_lit := E[H]−Var+W_geo+X
   = 0` identically** (`X=−2Cov`, `Cov=W_geo`, G4 Lemma G2); `G″_model := E[H]−Var+W_geo = 2W_geo → N/6`
   (constant, `j`-uniform). Both are `O(N)`; neither carries the `β_j` scale. **[T]**+**[V]**.
2. **The `β_j` scale lives in the coarse-pinned curvature.** `G″_cond = −d²/da² log ρ_v(a)|_0 = 1/(vᵀCv)
   = β_j Q_c + m → β_j Q_c`. This is the only object of the three whose normalization matches the
   Hessian; `Z_j^{-1}G″_cond → Q_c = c_L` is the `j`-uniform positive floor. **[T]**+**[V]**.
3. **Acid test, two faces.** The **IR** coarse mode (`k→0`) of the free massless field fails `(★)`
   (`Q(0)=0 ⇒ Var(ξ)=∞ ⇒ G″=0`), *as required for H1*. The **UV** block must **pass** on the free field
   (`Q_c=4π²/L²>0 ⇒ G″=β_jQ_c>0`): the free theory is UV-stable, so the opposite behaviour is correct.
   Conflating the two faces is the exact trap. **[T]**.

No new tool is manufactured here; this is bookkeeping of the *existing* `(★)`, but it removes a real
error risk (a constant `δ_W` silently read as renormalized) and pins the state L3 is in.

---

## 1 · The tower conventions (fixed once)

Wilson action `S_j = β_j Σ_p (1 − N^{-1}Re tr U_p)`, `β_j = 2N/g_j²`. Tower
(`C2` §2.1, `S5` §1): `g_j² = 1/(2β₀j)`, `β₀ = 11N/(48π²)`, so

```
   β_j = 4Nβ₀ j = 2N/g_j²  → ∞     (deep UV = large j),
   a_j = a₀ L^{-j}        (deep UV = small a_j),      j = log_L(a₀/a_j).
```

**Asymptotic freedom makes `β_j` logarithmic in `a_j`, not a power:** `β_j ∝ j ∝ log(1/a_j)`, and
`g_j² ∝ 1/log(1/a_j)`. The wave-function / field-strength renormalization `Z_j` is the log-divergent
object that absorbs this; in the free block it is exactly `Z_j = β_j` (§3).

Flat chart `X_ℓ = log U_ℓ ∈ su(N)`, Haar `dU = J dX`, `Φ = −log J = Σ_ℓ Φ_0(X_ℓ)`. Near the identity
(`SU(2)`: `Φ_0 ≃ ρ²/12`),

```
   Φ_0(X) ≃ (1/2) m |X|² + O(X⁴),     m := N/12,     Hess Φ_0 = m I = (N/12) I ,
```

which is `ch05` Lemma B (`Φ″ ⪰ (N/12)I`) read at the flat point.

---

## 2 · The free block, computed exactly

**The block.** A block of free links at shell `j`, flat boundary `U′ = 1`, small fluctuations `ξ`. The
quadratic Wilson action is

```
   S_j ≃ (1/2) ξᵀ H_j ξ ,      H_j = β_j Q ,       Q = lattice Maxwell form,  Q(k) = 2 Σ_μ (1 − cos k_μ) ,
```

with `Q` dimensionless and `O(1)` per link (bounded by `8`). The measure is exactly Gaussian in the flat
chart, `e^{−S_j−Φ}dξ`, with covariance

```
   C = (β_j Q + m I)^{-1},     m = N/12 .
```

Take a **unit coarse eigenmode** `v` of `Q` (eigenvalue `Q_c`): for a block of side `L` lattice sites,
the lowest nonzero mode has `Q_c = 2(1 − cos(2π/L)) → 4π²/L² =: c_L` (the free-box gap). All moments are
closed-form:

| object | exact free-block value | `β_j→∞` (deep UV) | scale |
|---|---|---|---|
| `H_vv = vᵀH_jv` | `β_j Q_c` | `→∞` | `Θ(β_j)` |
| `E[H_vv]` | `β_j Q_c` | `→∞` | `Θ(β_j)` |
| `Var(G_v)`, `G_v=⟨v,∇S_j⟩` | `β_j²Q_c² /(β_jQ_c+m)` | `β_jQ_c` | `Θ(β_j)` |
| defect `E[H]−Var` | `m β_jQ_c/(β_jQ_c+m)` | `m = N/12` | `Θ(1)` |
| `E[Φ″_vv]=m`; `Var(Φ′_v)=m²/(β_jQ_c+m)` | `m`; `→0` | `m`; `0` | `Θ(1)`; `o(1)` |
| **`W_geo = E[Φ″]−Var(Φ′)`** | `m β_jQ_c/(β_jQ_c+m)` | **`m = N/12`** | **`Θ(1)`, coupling-free** |
| `G″_model = E[H]−Var+W_geo` | `2m β_jQ_c/(β_jQ_c+m)` | `2m = N/6` | `Θ(1)` |
| `Cov(G,Φ′)` | `m β_jQ_c/(β_jQ_c+m)` | `m = N/12` | `Θ(1)` |
| `G″_lit = Δ_dyn+W_geo+X`, `X=−2Cov` | `0` (exact) | `0` | `0` |
| **`G″_cond = 1/(vᵀCv)`** | `β_jQ_c+m` | **`β_jQ_c`** | **`Θ(β_j)`** |

**Exact identities (machine-verified, §6).** The witness confirms, to `≲10⁻¹⁵`:
`Δ_dyn = E[H]−Var = W_geo` (G4 Lemma G2 = `Cov(G,Φ′)`), `G″_model = 2W_geo`, and `G″_lit = 0`. This is
the `S2`/`G4` algebraic structure reproduced on the free block: the *unconditioned / boundary-conditioned*
moments give `G″_lit=0` and `G″_model=2W_geo=O(N)`.

**The three objects, and which one is `(★)`.** Following G4 §1:
- `G″_lit` (boundary-conditioned, with `X`) `= 0` — the coarse-mode conditioning is load-bearing; the
  literal unconditioned `(★)` is `0≥δ_W`, false for any `δ_W>0`. **[T]**
- `G″_model = 2W_geo` (`S2`'s object) `→ N/6` — a true identity, but for the wrong functional (no coarse
  pinning); it is `O(N)`, **does not carry the `β_j` Hessian scale**, and misses the Savvidy negativity
  at the head (G4 §3). **[T]**
- `G″_cond` (`G4`'s (C) = `S3`'s `E[Σ]−D` at the coarse slice `a=a₀`) `→ β_j Q_c` — the **only** object
  with the Hessian normalization. This is the object L3/`(★)` must use. **[T]**

**Pure free field (acid test).** Set `Φ≡0` (abelian): `m=0`, so `defect=W_geo=0` and
`G″_lit=G″_model=0`; `G″_cond = β_jQ_c` (pure Gaussian stiffness). All three finite-`v` curvatures of the
*effective action* vanish with `Φ`; no false positivity is manufactured. **[T]**+**[V]**.

---

## 3 · Does the bare δ_W survive? The scaling

The question is a **normalization** question, and it has two inequivalent answers depending on which
`G″` one means — which is exactly why the campaign must be explicit.

### 3.1 Scale mismatch — a constant bare δ_W does not survive **[T]**

All terms that `(★)` compares carry the **bare lattice Hessian scale** `β_j = 2N/g_j²`:

```
   E[H] ~ Var ~ β_j  ~  g_j^{-2}  ~  log(1/a_j)      (Gaussian equipartition; §2),
   H_j  ~  β_j = 2N/g_j²  → ∞ .
```

A **constant** `δ_W = O(N)` is a fraction

```
   δ_W / E[H]  ~  O(1)/β_j  =  O(g_j²)  →  0
```

of the scale it sits under. After the field-strength renormalization `Z_j` that L3 applies (the
constant making `Z_j^{-1}H_j` finite on the core `D`; free-block value `Z_j = β_j`, §1), the **renormalized
margin** is

```
   δ̂_W := Z_j^{-1} δ_W  ~  O(1)/β_j  =  O(g_j²)  →  0 .
```

So a constant bare `δ_W` gives a renormalized floor that **vanishes like `g_j²`**. It survives as a
literal inequality (there is slack — `δ_W ≤ β_jQ_c` for large `j`) but **not** as a `j`-uniform positive
renormalized margin. This is the precise content of "the bare `δ_W` must scale."

### 3.2 The correct scaling — `δ_W = Θ(β_j) = Θ(g_j^{-2})` ∝ `log(1/a_j)` **[T]**

The object that tracks the Hessian is the coarse-pinned curvature `G″_cond = β_jQ_c + m → β_jQ_c`. A
floor taken as a fixed fraction of it,

```
   δ_W  =  c · G″_cond  =  c (β_j Q_c + m) ,     c ∈ (0,1) ,
```

scales as

```
   δ_W  =  Θ(β_j)  =  Θ(1/g_j²)  ∝  log(1/a_j)        (NOT a power of a_j — asymptotic freedom).
```

In terms of the task's two suggestions: the bare `δ_W` scales like **`∝ g_j^{-2}`** (a power of the
*coupling*), and only *logarithmically* like `a_j`: `δ_W ∝ log(1/a_j)`. It is **not** `∝ a_j^{-p}` for any
`p>0`; any power-law `a_j` piece is precisely what the renormalization `Z_j` subtracts (the `a⁻⁴`
triangle of `L2`/BPHZ, `G1-UNIFORM-ZA` §2), and belongs to the **subtracted** constant, never to `δ̂_W`.

Equivalently, the natural **dimensionless** bare floor is

```
   δ_W  =  Z_j · δ̂_W ,      δ̂_W = c · Q_c = c · c_L = c · 4π²/L²  =  Θ(1) ,
```

with `δ̂_W` the free-box gap — `j`-uniform, `a`-independent at fixed `L`. This is `S3`'s
`δ_W ≳ c_L + c_Φ − o(1)` and `TOWER-POINCARE`'s `c_L + δ_W` read correctly: the **renormalized** floor is
`Θ(1)`; the **bare** floor is `Θ(β_j)`. **[T]**

> **One-line answer.** The bare `δ_W` does **not** survive as a constant once `Z_j` is applied; it must
> scale like the renormalized Hessian, `δ_W ∝ β_j ∝ g_j^{-2} ∝ log(1/a_j)`, leaving the dimensionless
> `j`-uniform floor `δ̂_W = Z_j^{-1}δ_W = c_L = 4π²/L²`.

---

## 4 · The renormalized statement `(★_ren)` L3 actually needs

L3 is `sup_a Z_a^{-1}‖H_a‖_D < ∞` (`ch03` `eq:l3`, `G1-UNIFORM-ZA` §1). `(★)` is its per-shell
force-variance engine. Written on the **coarse-pinned** curvature (`G4` (C) = `S3` `E[Σ]−D`), and
normalized by `Z_j`, the correct statement is

```
 (★_ren)   Z_j^{-1} G″_j(v)  ≥  δ̂_W > 0 ,      ∀ v, j ≤ J, U′ ,
           G″_j(v) := −d²/da² log ρ^{(j)}_v(a)|_{a₀}
                    =  E_{μ(·|U′,a=a₀)}[∂²_v V_j] − Var_{μ(·|U′,a=a₀)}(∂_v V_j)
                    =  E_{μ(·|U′,a=a₀)}[Σ_j] − D_j        (S3 form, coarse slice a=a₀),
```
i.e. the **coarse-pinned marginal curvature** (`G4` (C)), equivalently the variance form
`Var_{μ(·|U′,a=a₀)}(⟨Pᵀ∇S_j,v⟩) ≤ E_{μ(·|U′,a=a₀)}⟨PᵀH_jP v,v⟩ + W_j(v) − δ_W` with `δ_W = Z_j δ̂_W`.

(For the *literal* boundary-conditioned object `E[H]−Var+W` — averaging over `a`, not pinning it — the
LHS is `Z_j^{-1}·2W_geo = O(N)/β_j → 0`: the **wrong** object. This is exactly the G4 correction that makes
the coarse-pinning in `G″_j` load-bearing; `(★_ren)` must be stated on the pinned curvature, which is
what `S3`'s `E[Σ]−D` already is.)

**Dimensionless.** In the flat chart `ξ` is dimensionless (lattice link log), so `S_j` is dimensionless,
`H_j=∂²S_j/∂ξ²`, `G_v=∂_vS_j` and `W_geo = E[Φ″]−Var(Φ′)` are all dimensionless numbers (lattice units);
`Z_j` is a pure ratio (renormalization constant). Hence both sides of `(★_ren)` are dimensionless. **[T]**

**`j`-uniform.** Two facts combine:
1. `Z_j^{-1}H_j = O(Q)` — bounded by construction of `Z_j` (this *is* L3's content at the operator level).
2. `Z_j^{-1}G″_cond = Q_c + m/β_j → Q_c = c_L`, a positive `j`-independent limit (the free-box gap).

So with `δ̂_W = c·c_L` `(0<c<1)`, `(★_ren)` holds **uniformly in `j` (and in `U′` off the
`e^{−cNj}`-small large-field set, `S5` §5; and in `a` at fixed `L`)** — for the free block. The
interacting content (the equipartition defect `Var−E[H]` at higher loops, `(H-SL)` of `S3`, the
`a⁻⁴` renormalon) is the open `[O]`. **[T free-block / C interactive / O remainder]**

**Why this is the form L3 needs.** `(★_ren)` is exactly `S3`'s `E[Σ] ≥ δ̂_W` (with `D≤0` from the fibre
gap `(H-C)`), i.e. `Z_j^{-1}E[Σ_j] ≥ δ̂_W`: the renormalized expected Schur complement stays above the
`j`-uniform free-box floor. Dividing by `Z_j` is what turns the `β_j`-stiff free-block Gaussian into the
`O(1)` free-box gap and makes the tower's `min_j` constant meaningful (= `(iii)` of L3,
`G1-UNIFORM-ZA` §1).

---

## 5 · The acid test — the two faces must not be confused

The charter's acid test (`f₁²=0`, `G″=0` on the free massless field) is calibrated for the **gap / H1**
statement. The **UV** statement has the *opposite* requirement. Precisely:

| | **IR face** (H1 / gap) | **UV face** (L3) |
|---|---|---|
| mode | coarse `k→0` of the physical measure | block-local mode at shell `j` |
| free-field `Q`-eigenvalue | `Q(0)=0` (zero mode) | `Q_c = 4π²/L² > 0` (lowest box mode) |
| free-field variance `vᵀCv` | `∞` (`Var(ξ)=∞`) | finite `≈1/(β_jQ_c)` |
| free-field `G″` | `1/∞ = 0` | `β_jQ_c > 0` |
| free-field `(★)` | **must FAIL** (`δ_W=0`) — required | **must PASS** — required |
| correct `δ_W` | `→0` on the free field (acid test) | `Z_j·c_L`, `j`-uniform (free UV-stable) |

- **IR face.** `G″_cond(k=0) = 1/(vᵀCv) = 0` because the massless zero mode has infinite variance: the
  free field manufactures no atom, `f₁²=0`, `G″=0`. This is the **correct** acid test for H1, and the one
  the charter mandates. **[T]**
- **UV face.** `G″_cond = β_jQ_c + m > 0`, and renormalized `Z_j^{-1}G″_cond → Q_c = 4π²/L² > 0`: the free
  UV block **passes** `(★)`, as it must (`C2` §5.1, `S5` §5.1, `S3` §4: "the free field is trivially
  UV-stable — L3 *passing* on the free field is correct"). **[T]**
- **The trap.** Applying the H1 acid test ("free field must fail") to the UV statement would demand
  `G″_UV=0` on the free field, i.e. would declare the *correct* UV stability a defect. The two faces
  differ in **scale** (IR `Q=0` vs UV `Q=4π²/L²`), not in algebraic form (`S3`: "same algebra, dual
  small blocks"). Do **not** carry the `f₁²=0 ⇒ vacuous` test across. **[T]** (+**[NG]** for the conflation)

---

## 6 · Witness — `research/v3/w7/v1_uv_normalization.py` (`py -3.12`, deterministic)

Closed-form free block (`N=2`, `4D` periodic `L=16`, `Q_c=0.152241`, `m=N/12=0.166667`), exact identities
to machine precision, and the scaling:

```
   j   beta_j    E[H]    Var(G)   defect=E[H]-Var   W_geo    Gpp_model   Gpp_lit    Gpp_cond  Gpp_cond/beta
   1   0.3715   0.0566   0.0143      0.042229      0.042229   0.084457   1.4e-17     0.2232    0.600859
  10   3.7151   0.5656   0.4369      0.128732      0.128732   0.257465  -1.7e-16     0.7323    0.197103
  40  14.8604   2.2624   2.1071      0.155231      0.155231   0.310462  -3.9e-16     2.4290    0.163456
                      ->              N/12=.16667   N/12       N/6=.3333    0          ~beta_j*Q_c  -> Q_c=.152241
```

- `defect = W_geo` exactly (G4 Lemma G2), `G″_model = 2W_geo` exactly, `G″_lit = 0` exactly — all three
  confirmed to `≲10⁻¹⁵` on the closed forms **and** on an explicit `d=6` dense Gaussian block
  (`|lit|=1.1e-15`, `|model−2W|=1.1e-15`, `|defect−W|=1.0e-15`; `G″_cond=3.19, 47.5, 472.8` at
  `β=2.5,40,400`).
- `G″_cond` grows like `β_j`; `G″_cond/β_j = Q_c + m/β_j → Q_c` (j-uniform, positive). **[V]**
- Acid test `Φ=0`: `W_geo = G″_model = G″_lit = 0`; `defect = 0`. **[V]**

The witness reproduces the `S2`/`S3`/`G4` algebra on the free block and exhibits the exact scaling
`δ_W ∝ β_j`, `δ̂_W = Z_j^{-1}δ_W → c_L`.

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| free-block Hessian `H_j=β_j Q`, `E[H]=Var=β_jQ_c+o(β_j)` (equipartition) | exact closed form | **[T]** (+**[V]**) |
| free-block Haar floor `W_geo = E[H]−Var → N/12` (coupling-free) | exact; `→N/12` | **[T]** (+**[V]**) |
| `G″_lit = Δ_dyn+W_geo+X ≡ 0`; `G″_model = 2W_geo → N/6` (wrong object, `O(N)`) | G4 G2 on free block | **[T]** (+**[V]** 1e-15) |
| coarse-pinned curvature `G″_cond = 1/(vᵀCv) = β_jQ_c+m → β_jQ_c` | exact; carries the Hessian scale | **[T]** (+**[V]**) |
| **bare constant `δ_W` does NOT survive**: renormalized margin `Z_j^{-1}δ_W ~ g_j² → 0` | §3.1 | **[T]** |
| **correct scaling `δ_W = Θ(β_j) = Θ(g_j^{-2}) ∝ log(1/a_j)`** (not a power of `a_j`) | §3.2 | **[T]** |
| `(★_ren)` dimensionless, `j`-uniform with floor `δ̂_W = Z_j^{-1}δ_W = c_L = 4π²/L²` | §4 (free block) | **[T]** |
| `(★_ren)` for the *literal* boundary-conditioned object | margin `→0` (wrong object) | **[NG]** (of that reading) |
| `(★_ren)` interacting content (equipartition defect, `a⁻⁴` renormalon, `(H-SL)`) | the open remainder | **[O]** |
| acid test: IR coarse mode `G″=0` (free, required fail); UV block `G″=β_jQ_c>0` (free, required pass) | §5 | **[T]** (+**[V]**) |
| applying the H1 acid test to the UV statement | category error | **[NG]** |
| `Z_j = Θ(β_j)` so `Z_j^{-1}H_j` bounded; interaction adds `(log)^γ` only | §1, §3 | **[T]** free / **[C]** interacting |
| **(★) itself / Clay** | unbuilt | **0%** |

---

## 8 · Sharpest line

**At the UV end the bare `(★)` is an inequality between quantities that all scale like
`β_j = 2N/g_j² = Θ(g_j^{-2}) = Θ(log(1/a_j))`; the only `j`-uniform object is the dimensionless ratio,
and it forces the floor to scale with the renormalized Hessian: `δ_W = Z_j δ̂_W` with
`δ̂_W = c_L = 4π²/L²` the free-box gap.** A constant bare `δ_W` is a `O(g_j²)`-fraction of the scale it
sits under and collapses to `0` after `Z_j` — so `(★_ren)`, not the bare `(★)`, is the form L3 needs:
`Z_j^{-1}G″_j(v) ≥ δ̂_W > 0` on the *coarse-pinned* curvature. The free block computes this exactly
(`G″_cond = β_jQ_c + m → β_jQ_c`, `Z_j^{-1}G″_cond → c_L`), and the acid test splits cleanly: the **IR**
zero mode fails `(★)` (`Q(0)=0`, `G″=0`, the H1 test) while the **UV** block passes it (`Q_c>0`,
`G″>0`, UV stability) — the two faces differ in scale, not form. No new tool, no QED. Clay 0%.

---

## 9 · Verification ledger (this note)

| item | check | status |
|---|---|---|
| tower `g_j²=1/(2β₀j)`, `β_j=2N/g_j²=4Nβ₀j`, `a_j=a₀L^{-j}` | `C2` §2.1, `S5` §1 | **[V]** (in-repo) |
| `β_j ∝ log(1/a_j)` (AF), not a power of `a_j` | `β_j ∝ j = log_L(a₀/a_j)` | **[T]** |
| free-block closed forms (H, Var, defect, W_geo, G″, G″_cond) and identities | `v1_uv_normalization.py` §A (j=1..40), §B dense `d=6` | **[V]** numerical |
| `Δ_dyn=W_geo`, `G″_model=2W_geo`, `G″_lit=0` | `G4-GENERAL-BOUND.md` §2 Lemma G2; witness residual `≲10⁻¹⁵` | **[T]**+**[V]** |
| `G″_cond=1/(vᵀCv)=β_jQ_c+m`; `Z_j^{-1}G″_cond→Q_c` | witness §A/§D | **[V]** numerical |
| acid test `Φ=0 ⇒ G″=0` | witness §C | **[T]**+**[V]** |
| UV free field must PASS (opposite of H1) | `C2` §5.1, `S5` §5.1, `S3` §4 | **[T]** (reading) |
| `(★)`/L3 open remainder (defect, renormalon, `(H-SL)`) | `G1-L3` §4, `S3` §4 | **[O]** |

*Refs.* `ch05-theorems.tex` (`eq:wgeo`, `eq:split`, `eq:vardecomp`, Lemma B `Φ″⪰N/12`); `ch03-existence.tex`
(§L3, `eq:l3`, §dim `G″_a=a²m²`); `C2-BALABAN-UNIFORM.md` §2.1–2.4, §5; `S5-N-UNIFORM.md` §1–5;
`G1-BLOCK-IDENTITY.md` §2–3, §6; `G4-GENERAL-BOUND.md` §1–3 (Lemma G1/G2/G3, (C) vs (U));
`S3-UV-FESHBACH.md` §1–2 (`G″=E[Σ]−D`, `(H-C)∧(H-SL)`); `G1-L3.md` §2–4; `G1-UNIFORM-ZA.md` §1–3;
`G1-L3-TOWER-POINCARE.md` §2–5; `G1-L3-UV-FORCE-BOUND.md` §1,§3,§6; `G3-PROJECTION-TERM.md` §6
(`grad S~β_jg_j`, `W_geo=O(1)`); witness `research/v3/w7/v1_uv_normalization.py`.
No third-party sources downloaded; all external inputs are standard and cited from the corpus.
