# G1-BLOCK-IDENTITY — the exact one-plaquette identity `G″=2W_geo` generalizes to a general block **iff** the direction is constant in the chart; a covariant/parallel-transported direction adds an explicit divergence defect

> 2026-11-24 · Generalization angle G1 on the single live inequality `(★)`.
> Reads: `S2-BOUNDARY-COUNTEREXAMPLE.md` §1–§2, `STAR-WAVE-SYNTHESIS.md`,
> `docs/paper-v2/ch05-theorems.tex` (eqs. `eq:wgeo`, `eq:split`, `eq:vardecomp`),
> `research/unconditional/PROOF-JACOBIAN-SAVES.md`, `research/unconditional/G1-L3-UV-FORCE-BOUND.md` §1.
> Witness: `research/v3/w7/g1_block_identity.py` (`py -3.12`; Gaussian block closed-form + exact SU(2) Haar quadrature).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional / proof-sketch · **[V]** numerical
> (deterministic; residual quoted) · **[NG]** no-go for a named class · **[O]** exact missing step.

---

## 0. Verdict (no QED)

**The identity generalizes — exactly and without modification — to a general block of many free links
with a fixed boundary `U′`, provided the direction `v` is constant in the flat chart (projection class
(a)).** Block size, the number of links, and the boundary `U′` are irrelevant to the identity: they enter
only through the action `S`, and the Haar potential `Φ` is a *per-link sum* because Haar measure factorises.
Precisely, for every constant `v` in the block tangent space,

```
             G″(v) := E[H_vv] − Var(G_v) + W_geo(v)  =  2 W_geo(v)      EXACTLY.        [T]
```

**It does *not* generalize to a configuration-dependent (covariant / parallel-transported) direction
`v = v(ξ)` (class (b)).** Then the exact identity acquires one extra, explicitly computable term:

```
             G″(v) = 2 W_geo(v) + Δ(v),
             Δ(v)  = E[g]² − E[φ]² + E[(∂_v v)·∇(Φ−S)] + E[(φ−g)·div v],
             g := ∂_v S,  φ := ∂_v Φ.                                                [T]
```

`Δ ≡ 0` for a constant field; it is generically non-zero otherwise (verified to `1.2×10⁻¹³` in a Gaussian
block and `≤1.6×10⁻⁷` on the real SU(2) Haar model). **`Δ` is exactly the "divergence of the collective
(gauge-projected) mode's vector field" named in `G1-L3-UV-FORCE-BOUND.md` §1** — the equipartition defect
`E[H]−Var` that would vanish identically in flat space for a constant direction. It can have either sign,
so it is *not* controlled by `W_geo` alone.

Consequences, stated with their labels:

1. **[T] Class (a) — constant projection.** `G″(v)=2W_geo(v)` exactly, for any block, any `U′`, any fixed
   mode `w = Pv` with `P` constant in the chart. Corollary: for the block,
   `(★) ⟺ 2W_geo(Pv) ≥ δ_W` — verbatim the S2 corollary, now block-level. **[T]**
2. **[T] Class (b) — covariant projection.** The identity *fails*; `G″=2W_geo+Δ`. This is a sharpening of
   S2's "survives any projection `P`": S2's `P` was implicitly constant (`w=Pv` fixed). A covariant `P` is
   not covered and breaks the clean form. **[T]** (+ **[V]** on the witness).
3. **[T] Chart edge & fixed boundary.** The Haar density vanishes on the boundary of every free link's
   chart domain, uniformly in `U′` (`S` bounded); **fixed boundary links produce no surface term** — they
   are not integrated. The only surviving surface term would come from an *interior* coarse constraint
   (pinning the block average), S2 §6 — still **[O]**.
4. **[T] Acid test.** Free field (`U(1)`, `Φ≡0`) gives `W_geo=0`; class (a) then has `G″=0` and class (b)
   with the Euler field `v=ξ` has `G″=−D<0`. No positive floor is manufactured. **[T]**+**[V]**
5. **[C]/[V] Positivity of `W_geo` at `β>0`** — the value/sign input — is *not* part of the identity and is
   unchanged from S2: `W_geo≥0` at `β=0` (Prékopa) and `W_geo>0` numerically for `β>0`, worst case at the
   Savvidy-aligned boundary. Same status as S2. **[C]**+**[V]**

Clay 0%. No QED.

---

## 1. General setup — the block measure

Let `𝔅` be a set of **free** links (the block interior), `U′` the **fixed** boundary (all links of the
lattice incident to `𝔅` but not in it). Write `φ = (U_ℓ)_{ℓ∈𝔅}` for the free links, `U_ℓ ∈ SU(N)`. Let
`d = N²−1` and `D = d·|𝔅|`.

**Flat (exponential) chart.** `U_ℓ = exp(ξ_ℓ)`, `ξ_ℓ ∈ su(N) ≅ ℝ^d`, in the fundamental domain
`𝒟_ℓ ⊂ ℝ^d` on which `exp` is a diffeomorphism onto `SU(N)` (the Weyl-chamber × maximal-torus cell). Let
`ξ = (ξ_ℓ) ∈ Ω := ∏_{ℓ∈𝔅} 𝒟_ℓ ⊂ ℝ^D`. Per-link Haar is `dU_ℓ = J(ξ_ℓ) d^dξ_ℓ`, and
`Φ(ξ) := −log J(ξ_ℓ)` summed, so

```
   ∏_{ℓ∈𝔅} dU_ℓ = e^{−Φ(ξ)} d^Dξ,        Φ(ξ) = Σ_{ℓ∈𝔅} Φ_ℓ(ξ_ℓ),    Φ_ℓ := −log J(ξ_ℓ).
```

**The key structural fact:** `Φ` is a **sum of single-link functions** — because Haar measure factorises
over links. Its Hessian is block-diagonal, `Hess Φ = ⊕_ℓ Hess Φ_ℓ`.

The conditional measure at fixed boundary is

```
   dμ(φ | U′) = Z(U′)^{−1} e^{−S(φ;U′)} e^{−Φ(φ)} d^Dξ ,     Z(U′) = ∫_Ω e^{−S−Φ} d^Dξ .
```

`S(φ;U′)` is the Wilson lattice action restricted to the block; it is a bounded trigonometric polynomial on
the compact group `SU(N)^{𝔅}`, for every fixed `U′`. The boundary enters **only** as a parameter of `S`.

The quantities in `(★)` (all directional, along a field `v = v(ξ)` in the tangent space `TΩ`):

```
   G_v := ∂_v S ,        H_vv := ⟨v , Hess S · v⟩ ,          (Hessian bilinear form)
   φ   := ∂_v Φ ,        Φ″_vv := ⟨v , Hess Φ · v⟩ ,
   W_geo(v) := E[Φ″_vv] − Var(φ) ,     G″(v) := E[H_vv] − Var(G_v) + W_geo(v) .
```

Here `∂_v = v(ξ)·∇_ξ`, while `H_vv, Φ″_vv` use the genuine Hessian form `vᵀ(Hess)v` (not `∂_v²`). This is
the definition of `G″` in `(★)` and in `ch05` eq. `eq:split`.

**Block-level `W_geo`.** Since `Φ = Σ_ℓ Φ_ℓ`,

```
   W_geo(v) = Σ_ℓ E[ v_ℓᵀ (Hess Φ_ℓ) v_ℓ ] − Var( Σ_ℓ v_ℓ·∇Φ_ℓ ).
```

The cross-link covariances `−Σ_{ℓ≠ℓ′} Cov(v_ℓ·∇Φ_ℓ , v_{ℓ′}·∇Φ_{ℓ′})` are physical (the blocks are coupled
through `S`) and are *not* dropped — they do not affect the identity.

---

## 2. The master integration by parts (arbitrary `v`)

Let `𝒱 := S + Φ`. Since `S` is bounded and `Φ_ℓ → +∞` on `∂𝒟_ℓ` (Vandermonde), `e^{−𝒱}` **vanishes on
`∂Ω`**. Then for any smooth `f` and any smooth vector field `v`,

```
   ∫_Ω (∂_v f) e^{−𝒱} d^Dξ  =  ∫_Ω f · ( ∂_v𝒱 − div v ) e^{−𝒱} d^Dξ ,      div v := ∇·v .
   (M)
```

*Proof.* `∫ v·∇f e^{−𝒱} = −∫ f ∇·(v e^{−𝒱})` (the boundary term is `∫_{∂Ω} f v e^{−𝒱}=0` because
`e^{−𝒱}|∂Ω=0`) `= ∫ f( v·∇𝒱 − div v) e^{−𝒱}`. ∎

For a **constant** `v`, `div v = 0` and `(M)` is the ordinary Stein identity; `∂_v² = ⟨v,Hess·v⟩` as well.
For a **varying** `v`, the divergence `div v` survives, and `∂_v² = ⟨v,Hess·v⟩ + (∂_v v)·∇`.

Everything below is a specialisation of `(M)` with `f ∈ {1, G_v, φ, ∂_v𝒱}`.

---

## 3. Case (a) — constant direction: the exact identity — **[T]**

Suppose `v ∈ ℝ^D` is constant in the chart (a constant-coefficient linear combination of single-link
directions; equivalently a fixed projection `w=Pv` with `P` a constant matrix). Then `div v=0`, `∂_v v=0`,
`∂_v² = ⟨v,Hess·v⟩`, and `(M)` gives the four Stein relations (identical in form to S2 Lemmas 1–2):

```
 (S0)  E[∂_v𝒱] = 0                     (f=1)
 (S1)  E[H_vv + Φ″_vv] = E[(∂_v𝒱)²]    (f=∂_v𝒱)          [ = Var(∂_v𝒱), by (S0) ]
 (S2)  E[Φ″_vv] = E[∂_v𝒱 · φ]          (f=φ)
 (S3)  E[H_vv]  = E[∂_v𝒱 · G_v]        (f=G_v)            [ = E[G_v²] + E[G_v φ] ]
```

**Theorem (block Stein identity).**
```
   Var(G_v + φ) = E[H_vv + Φ″_vv]                                   [T]
   E[G_v]       = −E[φ]                                             [T]
   Cov(G_v, φ)  = W_geo(v)                                          [T]
   E[H_vv] − Var(G_v) = W_geo(v)   ⇒    G″(v) = 2 W_geo(v)          [T]
```

*Proof.* `(S0)` is (S0). `(S1)` with (S0) is the first line. `(S0)` splits as `E[G_v]=−E[φ]`. For the
covariance: `(S2)` gives `E[G_v φ]=E[φ∂_v𝒱]−E[φ²]=E[Φ″_vv]−E[φ²]`, so
`Cov(G_v,φ)=E[Φ″_vv]−E[φ²]+E[φ]²=E[Φ″_vv]−Var(φ)=W_geo(v)`. For the Hessian: `(S3)` gives
`E[G_v²]=E[H_vv]−E[G_vφ]`, hence `Var(G_v)=E[H_vv]−E[Φ″_vv]+Var(φ)`, i.e.
`E[H_vv]−Var(G_v)=E[Φ″_vv]−Var(φ)=W_geo(v)`. Adding `W_geo` once more gives `G″=2W_geo`. ∎

**Why the block is irrelevant.** `(S0)–(S3)` are *pointwise algebraic* consequences of `(M)`; they never use
the number of links, the coupling structure, or the boundary `U′`. The only inputs are: (i) `e^{−𝒱}` vanishes
on `∂Ω` (chart edge, §5); (ii) `v` is constant. Since `Φ` is a per-link sum, `Hess Φ` is block-diagonal and
`Φ″_vv, φ` are the correct per-link sums. The theorem holds for **every** block, **every** `U′`, **every**
constant `v`. This is a direct, non-perturbative generalization of S2 §2.

**Why a constant projection `P` is harmless.** `(★)` uses `⟨Pᵀ∇S, v⟩ = ⟨∇S, Pv⟩ = ∂_{Pv}S` with
`w := Pv`. If `P` is constant in the chart, `w` is a constant field and §3 applies to `w` verbatim:
`G″(w)=2W_geo(w)`. This recovers and block-generalizes S2's "survives any projection `P`" — *for constant `P`*.

---

## 4. Case (b) — configuration-dependent direction: the defect — **[T]**

Now let `v = v(ξ)` depend on the field (the covariant / parallel-transported situation: a direction carried
across the block by the links themselves, `v_ℓ(ξ) = Ad_{g_ℓ(ξ)} v⁰_ℓ`, or any gauge-covariantly defined
mode). Then `div v ≢ 0` and `∂_v v ≢ 0`. Apply `(M)` with `f = ∂_v𝒱`, `φ`, `G_v`; keep the divergence and
the acceleration. Writing `g := ∂_v S = G_v`, `φ := ∂_v Φ`, `u := ∂_v𝒱 = g+φ`, `D := div v`,
`acc := ∂_v v`:

```
 (S0')  E[u] = E[D]                                    (f=1)
 (S1')  E[⟨v,Hess𝒱 v⟩] + E[acc·∇𝒱] = E[u²] − E[uD]    (f=u)
 (S2')  E[⟨v,HessΦ v⟩] + E[acc·∇Φ] = E[uφ] − E[φD]    (f=φ)
 (S3')  E[⟨v,Hess S v⟩] + E[acc·∇S] = E[u g] − E[gD]   (f=g)
```

Eliminating `E[u²], E[uφ], E[ug]` between these (proof in §4.1) yields the **corrected identity**

```
        G″(v) = 2 W_geo(v) + Δ(v),                                                   (★Δ)
        Δ(v)  = E[g]² − E[φ]²  +  E[acc·∇(Φ−S)]  +  E[(φ−g) D],        [T]
        acc = ∂_v v = v·∇v,   D = div v,   g = ∂_vS,  φ = ∂_vΦ.
```

`Δ` is a **divergence/acceleration defect**: it measures the failure of the direction field to be translation
invariant. For constant `v`, `acc=0`, `D=0`, and `E[g]=−E[φ]` (by S0'), so `Δ=0` and `(★Δ)` collapses to
§3 — the two cases are one formula. For a varying field `Δ` is generically non-zero and of either sign.

### 4.1 Algebra

From (S1'): `E[H_vv] = E[u²] − E[uD] − E[Φ″_vv] − E[acc·∇𝒱]`. Hence
`E[H_vv] − Var(g) = E[u²] − E[g²] − E[uD] − E[Φ″_vv] − E[acc·∇𝒱] + E[g]²`.
Since `u=g+φ`: `E[u²]−E[g²] = 2E[gφ]+E[φ²]`. So
`E[H_vv]−Var(g) = 2E[gφ] + E[φ²] − E[uD] − E[Φ″_vv] − E[acc·∇𝒱] + E[g]²`.
Adding `W_geo = E[Φ″_vv] − E[φ²] + E[φ]²`, the `E[Φ″_vv]` and `E[φ²]` terms cancel:
`G″ = 2E[gφ] − E[uD] − E[acc·∇𝒱] + E[g]² + E[φ]²`.
Now `2W_geo = 2E[Φ″_vv]−2E[φ²]+2E[φ]²`, and from (S2')
`E[Φ″_vv] = E[uφ]−E[φD]−E[acc·∇Φ] = E[gφ]+E[φ²]−E[φD]−E[acc·∇Φ]`. Substituting into `G″−2W_geo`, the
`E[gφ]` and `E[φ²]` terms cancel, and `2E[φD]−E[uD]=E[φD]−E[gD]`, `2E[acc·∇Φ]−E[acc·∇𝒱]=E[acc·∇(Φ−S)]`;
this leaves exactly `(★Δ)`. ∎

### 4.2 Witness confirmation

`research/v3/w7/g1_block_identity.py`:

* **Gaussian block, `D=6` coordinates, random non-commuting `A` (action), `K` (Haar), fixed boundary offset
  `μ`** — closed-form Gaussian moments. Constant `v`: `G″−2W_geo = 1.4×10⁻¹⁴`, independent of `μ`
  (`1.4×10⁻¹⁴` over four boundary offsets). Euler field `v=ξ`: `G″−2W_geo = 17.71… = Δ` to `1.2×10⁻¹³`.
* **Exact SU(2) one-link Haar quadrature, Euler field `v=ξ`** (`∂_v=ρ∂_ρ`, `div ξ=3`, `acc=ξ`),
  `Δ = ⟨g⟩²−⟨φ⟩²+4⟨ρ(Φ_ρ−S_ρ)⟩`:

```
  β     θ        G″            2W_geo        Δ            |G″−2W−Δ|
  0.5   0.000   −2.42031506   −4.18776828   1.76745327    4.6e-08
  1.0   π/2     −1.82251668   −3.12781075   1.30529413    5.1e-08
  2.5   π       +0.96843686   −1.48732622   2.45576309    1.6e-08
```

  The corrected identity `G″=2W_geo+Δ` holds on the real Haar model; `Δ≠0`; and crucially **`G″` is not
  sign-locked to `W_geo`** — the Euler defect makes it negative at small `β` (`G″=−2.42` at `β=0.5`,
  `θ=0`) even though the constant-directions of S2 are protected. A varying direction is *not* protected by
  the S2 mechanism.

**Reading.** S2's positivity result is a statement about **constant** directions (it says so: "take the
fixed direction `w=Pv`"). A genuinely covariant block mode is a *varying* field, and then the clean identity
is spoiled by `Δ` — the same object `G1-L3` §1 calls the divergence of the collective mode. This is a real,
labelled correction to the scope of S2, not a numerical artifact.

---

## 5. Chart edge and fixed-boundary surface terms — **[T]**

**(a) Chart-edge vanishing, every free link.** For `ℓ∈𝔅`, `J(ξ_ℓ) = |Δ_N(ξ_ℓ)|²·(regular)`, where the
Vandermonde `Δ_N ∝ ∏_{i<j} 2 sin((λ_i−λ_j)/2)` vanishes on the cut locus `λ_i−λ_j = ±2π` and on the Weyl
walls `λ_i=λ_j`. For SU(2), `J ∝ sin²(ρ/2)/ρ² → 0` at `ρ=2π` (the only singular edge) and `J=1` at the
regular interior point `ρ=0`. Thus `e^{−Φ_ℓ}→0` on **all** of `∂𝒟_ℓ`. Since `S` is a bounded trig
polynomial on the compact group for every fixed `U′`, `e^{−S}` is bounded, hence `e^{−𝒱}=e^{−S}e^{−Φ}→0` on
`∂Ω = ⋃_ℓ (∂𝒟_ℓ × ∏_{ℓ′≠ℓ}𝒟_{ℓ′})`. The IBP boundary term in `(M)` vanishes **uniformly in `U′` and in `v`.**
The interior point `ρ_ℓ=0` is regular and needs no special treatment.

**(b) No surface term from fixed boundary links.** The boundary links are **not integrated**; `U′` appears
only inside `S(φ;U′)` as a fixed, bounded parameter. The integration domain `Ω` and the Haar potential `Φ`
are *independent of `U′`*. Hence differentiation/IBP over the free links never meets `U′`; there is no
"surface at the boundary links." The identity is therefore uniform in `U′` **at the level of the
algebra** (its value depends on `U′` only through `S`, hence through `W_geo` and `Δ`).

**(c) The one caveat (unchanged).** A genuine surface term appears only if one restricts the free integral
to a *sub-domain* of `Ω` — an **interior coarse-projection constraint** (pinning the block average / a
coarse block-spin to a saddle), which truncates the integral with a non-vanishing edge `e^{−𝒱}|_{∂(sub)}≠0`.
That is exactly the S2 §6 hiding place, and remains **[O]**. §5(a)–(b) say it is the *only* one: neither
boundary holonomy nor the chart edge can generate a surface defect.

---

## 6. Corollary and acid test

**Corollary (block `(★)`, class (a)).** Let `P` be constant in the chart and `w=Pv`. Then for the block,
```
   (★):  Var(G_w) ≤ E[H_ww] + W_geo(w) − δ_W
      ⟺  E[H_ww] − Var(G_w) + W_geo(w) ≥ δ_W
      ⟺  G″(w) ≥ δ_W
      ⟺  2 W_geo(w) ≥ δ_W .                                                          [T]
```
So, exactly as in the one-plaquette model, in projection class (a) the whole UV inequality reduces to the
**strict positivity and uniform lower bound of the interacting Haar curvature along the blocked mode `w`.**
For class (b) the reduction is `2W_geo(w)+Δ(w) ≥ δ_W`, and `Δ` must be carried (it is not sign-definite).

**Acid test (free field must fail).** Take the group `U(1)^{𝔅}` (or `ℝ^D`): `J≡1`, `Φ≡0`, so `W_geo≡0`.
* Class (a), Gaussian `S`: `S` quadratic ⇒ `Var(G_w)=E[H_ww]` exactly ⇒ `G″=0`. `(★)` needs `0≥δ_W`,
  false for any `δ_W>0`. **[T]** (witness: `2.8×10⁻¹⁴`).
* Class (b), Euler `v=ξ`, isotropic Gaussian: `G″ = D − 2D = −D < 0` (witness: `−6.000000` at `D=6`).
  **[T]**+**[V]**.
* SU(2) pure Haar (`β=0`): `W_geo=0` by `ξ→−ξ` symmetry, `G″→0` as `β↓0` (`∝β²`, S2 §3). The positive
  floor is generated **entirely by the Wilson interaction** `β>0`.

Thus the generalization manufactures no positivity for the free theory; the gap-causing ingredient is still
the Haar curvature `W_geo` (and, for covariant modes, the defect `Δ`), exactly as the campaign identifies.

---

## 7. Status ledger

| item | statement | label |
|---|---|---|
| master IBP `(M)` for arbitrary `v`, boundary `e^{−𝒱}\|∂Ω=0` | §2 | **[T]** |
| block Stein identity `Var(G_v+φ)=E[H_vv+Φ″_vv]`, `Cov=W_geo`, `G″=2W_geo` for **constant** `v`, any block, any `U′` | §3 | **[T]** (+**[V]** 1.4e-14) |
| generalization is **block- and boundary-independent** (only `S` carries them; `Φ` per-link sum) | §1, §3 | **[T]** |
| corrected identity `G″=2W_geo+Δ` for **varying** `v`; `Δ=E[g]²−E[φ]²+E[acc·∇(Φ−S)]+E[(φ−g)D]` | §4 | **[T]** (+**[V]** 1.2e-13 / 1.6e-7) |
| chart-edge vanishing per free link; **no** surface term from fixed `U′`; only interior constraint gives one | §5 | **[T]** |
| corollary `(★) ⟺ 2W_geo(Pv)≥δ_W` for the block, constant `P` | §6 | **[T]** |
| free-field failure (`W_geo=0 ⇒ G″=0` (a), `=−D` (b)) | §6 | **[T]**+**[V]** |
| `W_geo>0` for `β>0`, worst case Savvidy-aligned boundary | S2 §3, unchanged | **[C]**+**[V]** |
| identity holds *unchanged for every projection class* (incl. covariant `P`) | refuted: class (b) has `Δ≠0` | **[NG]** (of that reading) |
| interior coarse-projection surface term / `(★dec)` | untested; only remaining hiding place | **[O]** |

**Sharpest line.** `G″=2W_geo` is a theorem of the **constant-direction** class — and it is exactly that:
the flat-chart Stein relations are algebraic and know nothing of the block's size, its boundary, or its
couplings, so the one-plaquette identity lifts verbatim to any block. The moment the direction is carried
by the field (the physically covariant case), the identity picks up the divergence-and-acceleration defect
`Δ` — the same "divergence of the collective mode" that `G1-L3` §1 isolates. So the S2 reduction survives
blocking **iff** the coarse projection is constant; a covariant projection is a *strictly larger* object whose
extra term is now explicit and labelled. The only uncomputed hidden place remains the interior
coarse-projection boundary term. Clay 0%.

---

## 8. Verification ledger

| item | checked | how |
|---|---|---|
| block measure, `Φ=Σ_ℓΦ_ℓ`, `Hess Φ` block-diagonal | exact | Haar factorisation; `Φ_ℓ=−log J_ℓ` |
| `(M)` master IBP, boundary `e^{−𝒱}\|∂Ω=0` | exact | §2; chart-edge vanishing §5(a) |
| class (a) `G″=2W_geo`, boundary-independence | `1.4×10⁻¹⁴`, `D=6`, non-commuting `A,K`, 4 boundary offsets | `g1_block_identity.py` Part A |
| class (b) `G″=2W_geo+Δ`, `Δ≠0` (Gaussian) | residual `1.2×10⁻¹³`, `Δ=17.71` | Part A |
| class (b) `G″=2W_geo+Δ`, `Δ≠0` (real SU(2) Haar, Euler `v=ξ`) | residual `≤1.6×10⁻⁷`, `Δ` listed | Part B |
| acid test: free field `G″=0` (a) / `−D` (b) | `2.8×10⁻¹⁴` / `−6.000000` | Part A |
| chart-edge `J→0` at cut locus (SU(2): `ρ=2π`), `S` bounded | analytic | §5(a) |
| `W_geo>0` for `β>0` (value input) | inherited from S2 (worst case `θ=π`, transverse, `δ_W≈0.13` at `β=2.5`) | S2 §3 |

*Refs.* `S2-BOUNDARY-COUNTEREXAMPLE.md` §1–2, §6; `STAR-WAVE-SYNTHESIS.md` §1–2; `docs/paper-v2/ch05-theorems.tex`
eqs. `eq:wgeo`, `eq:split`, `eq:vardecomp`; `research/unconditional/PROOF-JACOBIAN-SAVES.md`;
`research/unconditional/G1-L3-UV-FORCE-BOUND.md` §1; witness `research/v3/w7/g1_block_identity.py`.
