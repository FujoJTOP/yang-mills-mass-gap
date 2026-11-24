# G2-TWO-PLAQUETTE — `R = G″ − 2W_geo = 0` exactly for every linear direction in a two-plaquette / two-free-link block; nonzero (and the β=0 acid test fails) for **field-dependent** directions

> 2026-11-24 · Follow-up to `S2-BOUNDARY-COUNTEREXAMPLE.md` (the one-plaquette/one-link exact model).
> Reads: `S2-BOUNDARY-COUNTEREXAMPLE.md`, `docs/paper-v2/ch05-theorems.tex` (§`sec:wgeo`, `eq:wgeo`,
> `eq:split`, `eq:vardecomp`), `STAR-WAVE-SYNTHESIS.md` §2.
> Witness: `research/v3/w7/g2_two_plaquette.py` (`py -3.12`; uniform-chart sampling + Haar sampling,
> validated against the S2 exact quadrature).
>
> **Labels (charter §2).** **[T]** genuine theorem · **[C]** conditional · **[V]** numerical ·
> **[NG]** no-go · **[O]** exact missing step.

---

## 0. Verdict (no QED)

The one-plaquette identity `G″ = 2W_geo` (S2) is **not** a one-link accident. It survives the step to a
genuine **two-free-link** block, and the reason is structural. Precisely:

1. **[T] `R=0` for every linear direction, any block.** In the fixed-boundary two-plaquette block of
   §1 the generalized flat-chart Stein identity (§2.2) gives, for **every constant** direction
   `v=(v_A,v_V)∈ℝ⁶` — in particular `D3=(e₃,e₃)/√2`, a constant linear combination of **both** free
   links (the *linear-P* case) —
   ```
     Cov_μ(G_v,Φ′_v) = W_geo(v),   E_μ[H_vv] − Var_μ(G_v) = W_geo(v),
     ⇒  G″(v) := E[H]−Var(G)+W_geo = 2W_geo(v),   R(v) := G″−2W_geo = 0.
   ```
   This holds for **any** number of free links, **any** fixed boundary `B`, **any** `β`; the only
   hypothesis is that the fibre density `e^{−S}J` vanishes on the chart edge (§2.1), which is
   dimension-independent and does **not** use that the measure factorizes — so it survives the `A–V`
   coupling. **Hence `(★)` in the boundary-conditioned fibre is again exactly `2W_geo(v|U′) ≥ δ_W`,
   with the S2 `δ_W` structure; the block size and the extra free link cost nothing.** `[T]`
2. **[V] Linear scan.** On the coupled two-link model, `R` is zero to Monte-Carlo noise `~10⁻³`
   (`~10⁻⁴` where the heavy tail is quiet) over `β∈{0.5,1,2.5}`, `θ∈{0,π/2,π}` and six directions,
   while `W_geo>0` and `G″=2W_geo`. The one-link code reproduces the S2 **exact** quadrature to
   `≤5×10⁻⁴` (§3.1). (§4.)
3. **[T] A field-dependent (“covariant”) direction is *forced* out by the acid test.** The two exact
   IBP identities read (§2.3)
   ```
     (IBP-1)  E_μ[∂_v𝒱]   = E_μ[∇·v],
     (IBP-2)  E_μ[∂_v²𝒱]  = E_μ[(∂_v𝒱)²] − E_μ[(∂_v𝒱)(∇·v)] ,
   ```
   where `𝒱=S+Φ`, and the correction terms carry the **divergence of the direction field** `∇·v`.
   They vanish **iff `v` is constant** (`∇·v=0`). For any **field-dependent** `v` the β=0 acid test
   `G″(β=0)=0` **fails**: on pure Haar with `S=0` one has `G″(β=0)=W_geo(v|β=0)`, which is `0` for
   constant `v` but **nontrivial** otherwise —
   * bounded field-dependent (radial `v=t̂`, `∇·v=2/ρ>0`): `W_geo(β=0)=+0.151>0` ⇒ **false positivity**
     on the free theory;
   * unbounded field-dependent (right-invariant `U→Ue^{εw}`, `|v|~1/ε` at the chart antipode):
     `W_geo(β=0)=−0.665<0` ⇒ **wrong sign**, and `R=+0.665`.
   So the exact reduction and the acid test **single out the linear (constant) direction class**:
   every constant `v` gives `R=0` and `G″(β=0)=0` `[T]`; every field-dependent `v` examined breaks
   both, and the acid test rules field-dependence out generically (it needs `E_0[Φ′_v]=0` and
   `∇·v=0` simultaneously, which constant `v` supplies and a generic field does not). `[T]`
4. **[V] Covariant scan.** For the right-invariant direction `R≠0`, `O(0.1–0.6)`, sign varying with
   `β,θ`; `G″≈W_geo<0` (so `(★)` fails badly — the covariant direction is **not** a valid tool).
   For the bounded radial direction `R` is `O(0.03–0.15)`. (§5.)
5. **[O] Reading for the campaign.** The remainder `R` is nonzero precisely through the `∇·v`
   corrections to IBP-1/IBP-2 (equivalently the cross term `X=−2Cov(G,Φ′)` of S2), the one surviving
   “interior coarse-projection” term S2 §6 flagged. The UV
   shell uses `⟨Pᵀ∇S,v⟩` with a **constant coarse direction** `v`, i.e. exactly the linear-P class
   where `R=0`; a **field-dependent** coarse map (a covariant block-spin) is the one that revives
   `R≠0` — and it is automatically excluded by the acid test. This sharpens, rather than closes, the
   gap: the `(★)` input is *not* obstructed by block size, but it is a theorem about **linear**
   projections only. Clay 0%.

---

## 1. The model (two free links, fixed boundary)

SU(2) lattice gauge theory on a `1×2` strip. Two plaquettes share the vertical link `V`; the shared
link `V` and a horizontal link `A` are **free**, all other links are **fixed boundary**. Choosing the
fixed boundary products so the plaquette holonomies read
```
   P₁ = A V ,            P₂ = B V⁻¹ ,        A, V ∈ SU(2) free,   B ∈ SU(2) fixed boundary,
```
the Wilson action and the conditioned measure are
```
   S(A,V) = β(1 − ½tr P₁) + β(1 − ½tr P₂) = β( 2 − ½tr(AV) − ½tr(BV⁻¹) ),
   dμ(·|B) = Z_B⁻¹ e^{−S(A,V)} dA dV ,        dA, dV Haar.
```
This **is** `μ(·|U′)` with `U′=B`: fixed boundary = fixed `U′`, integrate the two interior links. It
is the minimal block with **two** free links (the S2 model had one).

**Structural fact — the model factorizes in the *right* variables [T].** Put `W:=AV`. Haar
translation invariance gives `dA=dW`, so in `(W,V)` the measure is a product of two one-plaquette
measures,
```
   dμ = [ Z_I⁻¹ e^{−β(1−½tr W)} dW ] ⊗ [ Z_B⁻¹ e^{−β(1−½tr(BV⁻¹))} dV ] ,   W ⊥ V ,
```
`W` with trivial boundary `I`, `V` with boundary `B`. (Verified: `E[½tr AV]=E[½tr BV⁻¹]` equals the
one-plaquette value at the same `β,θ`, §3.1.) **But** the campaign’s free links are `A` and `V`, and
the flat chart is `ξ_A=log A`, `ξ_V=log V` — **not** `log W`. Derivatives `∂_{v_A}` at fixed `V`
couple `A,V` through `W=AV`, so the two-link chart is genuinely non-factorized even though the
measure factors in `(W,V)`. This is exactly why the block is a real test and not two copies of S2.

**Two direction types.**
* **Linear (flat) / “linear P” [T].** A **constant** `v=(v_A,v_V)∈ℝ³⊕ℝ³`: `ξ(ε)=ξ+εv`. This is
  `⟨Pᵀ∇S,v⟩` with `P` linear; the coarse direction `w=Pv` is a constant linear combination of the
  free links. Includes `A3,V3,A1,V1` (single-link) and `D3,D1=(e,e)/√2` (diagonal, **both** links).
* **Field-dependent / “covariant” [T]+[V].** e.g. radial `v=t̂` (bounded), or the right-invariant
  (gauge-covariant) generator `U(ε)=U·e^{εw}` (bounded in the group, unbounded in the chart).

---

## 2. The exact identities — **[T]**

Per free link, flat chart `ξ_i∈B(2π)⊂ℝ³`, `Haar = J(ξ_i)d³ξ_i`,
`Φ(ξ_i)=−log J(ξ_i)=2log(ρ_i/2)−2log sin(ρ_i/2)`, `J(ρ)=sin²(ρ/2)/ρ²` (`J(0)=¼`, `J(2π)=0`). Put
`ξ=(ξ_A,ξ_V)∈D=∏_i B(2π)`, `𝒱=S+Φ`, `Φ=Σ_iΦ(ξ_i)`, `dμ=Z⁻¹e^{−𝒱}dξ`.

### 2.1 Edge-vanishing (§2.1 of S2, dimension-independent)
On a chart face `|ξ_i|=2π` (`ε=2π−ρ_i→0`): `e^{−Φ}=J~ε²`, `|∂Φ|~1/ε`, `|∂²Φ|~1/ε²`, so
`e^{−Φ}|∂Φ|~ε→0`, `e^{−Φ}|∂²Φ|=O(1)`. Since `S` is bounded, the edge flux `e^{−𝒱}(∂_v𝒱)(v·n)` and
`e^{−𝒱}` vanish on `∂D` for any **bounded** direction field `v` (`|v|=O(1)`) and any action —
`|flux| ≲ J|v||∇Φ| ~ |v|·ε → 0`. This is the only structural hypothesis; it holds for every block
size and every fixed boundary. **An unbounded field** (the right-invariant generator has `|v|~1/ε` at
the antipode) has `|flux| ~ O(1)` and violates it — see §5.1, where this is exactly what breaks
`R=0`.

### 2.2 IBP identities (valid for any bounded direction field with vanishing edge flux)
```
   (IBP-1)  E_μ[∂_v𝒱]  = E_μ[∇·v],
   (IBP-2)  E_μ[∂_v²𝒱] = E_μ[(∂_v𝒱)²] − E_μ[(∂_v𝒱)(∇·v)] .
```
*Proof.* `(v·∇𝒱)e^{−𝒱} = −∇·(v e^{−𝒱}) + (∇·v)e^{−𝒱}` ⇒ IBP-1. And
`v·∇( (v·∇𝒱)e^{−𝒱} ) = (∂_v²𝒱)e^{−𝒱} − (∂_v𝒱)²e^{−𝒱}`; integrating the divergence
`∇·( (∂_v𝒱)e^{−𝒱} v )` and using `∇·v` gives IBP-2. ∎
**For constant `v`, `∇·v=0`, so both corrections vanish:** `E[∂_v𝒱]=0`, `E[∂_v²𝒱]=Var(∂_v𝒱)`.

### 2.3 Reduction and `R=0` for **linear** directions
Constant `v`: `E[∂_v𝒱]=0`, `E[∂_v²𝒱]=Var(∂_v𝒱)`, and the cross-Stein `E[∂_v𝒱·Φ′_v]=E[Φ″_vv]`.
With `G_v=∂_v𝒱−Φ′_v`, `H_vv=∂_v²𝒱−Φ″_vv`, this reproduces the S2 corollary verbatim:
`Cov_μ(G_v,Φ′_v)=W_geo(v)`, `E[H]−Var(G)=W_geo(v)`, `G″=2W_geo(v)`, **`R=0`**. ∎
`W_geo(v)=E_μ[Φ″_vv]−Var_μ(Φ′_v)` (full conditioned measure `μ`, S2 convention). The derivation uses
**no** property of `S` beyond boundedness and **no** count of free links ⇒ holds for the two-link
block and for any `m`-link block. Projection `P` is harmless per-direction (`w=Pv` fixed).

### 2.4 The acid test localizes the direction class — **[T]**
At `β=0`: `S=0 ⇒ G=H=0 ⇒ G″(β=0)=W_geo(v|β=0)`.
* Constant `v`: `E_0[Φ′_v]=v·E_0[∇Φ]=0` (isotropy) and `∇·v=0`, so by IBP-2
  `E_0[Φ″]=E_0[Φ′²]` and `W_geo(β=0)=E_0[Φ″]−Var_0(Φ′)=E_0[Φ′]²=0`. **Acid test passes.**
* Field-dependent `v`: `E_0[Φ′_v]≠0` in general. If the edge flux still vanishes (bounded `v`, e.g.
  radial), IBP-2 gives `E_0[Φ″]=E_0[Φ′²]−E_0[Φ′(∇·v)]`, hence
  `W_geo(β=0)=E_0[Φ′]²−E_0[Φ′(∇·v)]≠0` in general (radial: `+0.151`) — a **false positive on the
  free theory**. If the direction blows up at the antipode (right-invariant `v`), even the edge flux
  fails to vanish and `W_geo(β=0)` is `O(1)` and **negative** (`−0.66`). **Acid test fails in both
  cases.**
Hence **`R=0` and `G″(β=0)=0` hold simultaneously for every constant (`∇·v=0`) direction** (the
linear-P class), and fail generically for a field-dependent one. This is the sharp structural
statement of G2.

---

## 3. Numerical method and validation

Witness `g2_two_plaquette.py`. All expectations are Monte-Carlo over the two free links, with two
sampling modes (both unbiased self-normalized importance sampling):
* **uniform chart**: `ξ_i` uniform in `B(2π)`, weight `e^{−S}J(ρ_A)J(ρ_V)∝e^{−𝒱}` (finite variance
  for `Φ″`,`Φ′²` because the `J~ε²` weight beats the `1/ε²` of `Φ″` — good at low `β`);
* **Haar**: `ξ_i` from Haar (density `J`), weight `e^{−S}` (the antipode is rare in Haar — good at
  high `β`).
`G_v,H_vv` by central differences in the flat flow (a `β>0` action is smooth); `Φ′_v,Φ″_vv` analytic
(`phi_grad_hess`) for **constant** `v`, and by central differences along the true flow for a
field-dependent `v`. Independent seeds are used for `A` and `V`.

### 3.1 One-link validation [V]
The same code, one free link, reproduces the S2 **exact** deterministic quadrature (S2 §3):
```
   b=0.5 th=0   e3: W=0.04513 (exact 0.04552)  G″=0.09051 (exact 0.09104)  R=+0.00024
   b=2.5 th=0   e3: W=0.13183 (exact 0.13196)  G″=0.26318 (exact 0.26391)  R=−0.00047
   b=2.5 th=π   e1: W=0.06566 (exact 0.06587)  G″=0.13145               R=+0.00014
```
Agreement to `≤5×10⁻⁴`: the machinery is calibrated.

### 3.2 Factorization check [V]
`E[½tr AV]=0.1245`, `E[½tr BV⁻¹]=0.1229` at `β=0.5, θ=0` (both = one-plaquette value); and
`0.5110 / 0.5120` at `β=2.5, θ=0`. `W=AV` and `V` are independent, as §1 predicts.

---

## 4. Linear directions — `R=0` — **[T]**+**[V]**

Representative rows (uniform sampling, `n=2×10⁶`; `E[H]−VarG = E_μ[H_vv]−Var_μ(G_v)`). The theorem
says `E[H]−VarG = W_geo` and `G″=2W_geo`, i.e. `R=0`.

| β | θ | dir | `E[H]−VarG` | `W_geo` | `G″` | `R=G″−2W` |
|---|---|---|---|---|---|---|
| 0.5 | 0 | A3 | 0.00609 | 0.00621 | 0.01230 | −0.0001 |
| 0.5 | 0 | V3 | 0.04550 | 0.04538 | 0.09088 | +0.0001 |
| 0.5 | 0 | **D3** | 0.04455 | 0.04454 | 0.08910 | **+0.00001** |
| 0.5 | π/2 | D3 | 0.03984 | 0.03983 | 0.07966 | +0.00001 |
| 0.5 | π | V3 | 0.00862 | 0.00837 | 0.01699 | +0.00025 |
| 0.5 | π | **D3** | 0.02505 | 0.02503 | 0.05008 | **+0.00001** |
| 1.0 | 0 | A3 | 0.02236 | 0.02281 | 0.04517 | −0.0005 |
| 1.0 | 0 | V3 | 0.07952 | 0.07966 | 0.15918 | −0.0001 |
| 1.0 | 0 | **D3** | 0.08237 | 0.08283 | 0.16520 | **−0.0005** |
| 1.0 | π/2 | V3 | 0.07136 | 0.07124 | 0.14260 | +0.0001 |
| 1.0 | π | A3 | 0.00167 | 0.00186 | 0.00353 | −0.0002 |
| 1.0 | π | **D3** | 0.05437 | 0.05441 | 0.10878 | **−0.00003** |
| 2.5 | π | V3 | 0.13200 | 0.13065 | 0.26265 | +0.0014 |
| 2.5 | π | **D3** | 0.13752 | 0.13661 | 0.27413 | **+0.0009** |

At `β=0.5,1.0` the residual is `≤10⁻³` on every row; at `β=2.5, θ∈{0,π/2}` the heavy `Φ″` tail at
the antipode inflates the uniform-sampling noise to `~10⁻²` on the largest-`W` directions (`V1,D3,D1`)
— a sampling artefact, removed by Haar sampling at high `β` (see note). `W_geo>0` throughout.
The factorization makes `E[½tr AV]` and `E[½tr BV⁻¹]` one-plaquette values (checked §3.2).

**Acid test, linear (β=0, pure Haar)** — must give `G″=0`:
```
   β=0 A3: E[H]−VarG=0.000000  W_geo=+0.000012  G″=+0.000012
   β=0 D3: E[H]−VarG=0.000000  W_geo=−0.000138  G″=−0.000138
```
`G″=0` to `~10⁻⁴`. **Pass.** The free/flat `U(1)` group (`ad=0, Φ=0`, Gaussian `S ⇒ Var(G)=E[H]`)
gives `G″=0` identically — the further required failure. `[T]`

*Note on `β=2.5` sampling:* the flat-chart `Φ″` has a genuine `~1/ε²` singularity at the antipode
(`ε=2π−ρ`); the mean is finite but the naive estimator is heavy-tailed. Haar sampling (antipode rare)
gives at `β=2.5` (`n=2.5×10⁶`): `θ=0` `R(A3,V3,D3)=−0.0001,−0.0004,−0.0004`; `θ=π/2`
`≈−0.0003…−0.0017`; `θ=π` `≈−0.0001…−0.0037`. So `R=0` to `~10⁻⁴` at `θ=0`, `~10⁻³` at `θ=π`. The
uniform-sampling `~10⁻²` at `β=2.5, θ∈{0,π/2}` is purely this tail (and is removed by Haar sampling —
note the sign flips with sampling, confirming it is noise around `0`, not a systematic `R`).

---

## 5. Field-dependent (“covariant”) directions — `R≠0`, acid test FAILS — **[T]**+**[V]**

### 5.1 right-invariant `U(ε)=U·e^{εw}` (gauge-covariant; chart-norm `|v|~1/ε` at the antipode)
```
   β     θ      dir   E[H]−VarG     W_geo       G″          R
   0.5   0      C3    −0.00010     −0.47300   −0.47310   +0.47290
   0.5   0      C30    0.00001     −0.31384   −0.31383   +0.31385
   0.5   π/2    C3    −0.00049     −0.50808   −0.50857   +0.50759
   0.5   π      C3    −0.00022     −0.60882   −0.60904   +0.60859
   2.5   0      C3    −0.00732     +0.06823   +0.06091   −0.07555
   2.5   π/2    C30   −0.00400     −0.11814   −0.12214   +0.11414
   2.5   π      C3    −0.00094     −0.26043   −0.26137   +0.25949
```
`C3` = right-fluctuate **both** links by `w=e₃`; `C30` = only link `A`. `E[H]−VarG≈0` (the action
contributes almost no dynamic curvature along this covariant direction), so
`G″≈W_geo` and `R=−W_geo` to leading order. `R` is `O(0.1–0.6)`, sign following `−sign(W_geo)`.
**Acid test, β=0 (pure Haar):**
```
   β=0 C3 : W_geo=−0.664914  G″=−0.664914   (E[H]−VarG=0)
   β=0 C30: W_geo=−0.332181  G″=−0.332181
```
`G″(β=0)=−0.66 ≠ 0`: **FAIL** (and negative). Diagnostic: along this direction `E_0[Φ″]≈0.0001` while
`E_0[Φ′²]≈0.33`, so `E_0[Φ″]≠E_0[Φ′²]` — the IBP boundary flux is `~|v|·ε ~ O(1)` and does **not**
vanish, because `|v|~c/ε` at the antipode. The negative `W_geo` is this chart pathology, *not* the
`W_geo≥0` theorem (which is a statement about **linear** marginals, Prékopa; §2.3 of ch05).

### 5.2 bounded radial `v=t̂` (one-link, to isolate “field-dependence” from “chart blow-up”)
```
   β=0    θ=0 : E[H]−VarG=0.00000  W_geo=+0.15096  G″=+0.15096  R=−0.15096   (E_0[Φ′]=+0.777)
   β=0.5  θ=0 : E[H]−VarG=0.01337  W_geo=+0.15167  G″=+0.16504  R=−0.13830
   β=2.5  θ=0 : E[H]−VarG=0.24651  W_geo=+0.17279  G″=+0.41930  R=+0.07373
   β=2.5  θ=π : E[H]−VarG=0.17662  W_geo=+0.20591  G″=+0.38253  R=−0.02929
```
Here the flux **does** vanish (`|v|=1`), IBP-2 holds, and
`W_geo(β=0)=E_0[Φ′]²−E_0[Φ′(∇·v)]=+0.151>0` (with `E_0[Φ′]=+0.777`, `∇·v=2/ρ`) exactly as §2.4
predicts — and `G″(β=0)=+0.151>0`: a **false positivity on the free theory**. `R≠0` (`O(0.03–0.15)`)
even though `|v|` is bounded, caused purely by `∇·v=2/ρ≠0` in IBP-1/2. This cleanly separates the two
failure routes: field-dependence alone (`∇·v≠0`) already breaks `R=0` and the acid test; the extra
chart blow-up of the right-invariant field only flips the sign of `W_geo`.

### 5.3 Summary of the field-dependent case [T]
```
   direction            ∇·v      W_geo(β=0)       G″(β=0)     acid test    R(β>0)
   constant (linear)    0        0                0           PASS         0  (exactly)
   radial t̂ (bounded)   2/ρ>0    +0.15 (>0)       +0.15       FAIL(>0)     O(0.03–0.15)
   right-invariant      −Φ′      −0.66 (chart)    −0.66       FAIL(<0)     O(0.1–0.6)
```
**One-plaquette limit.** For a **linear** direction `R=0` at every block size, so it “vanishes” in the
one-plaquette limit only in the trivial sense that it is already identically zero. For a
**field-dependent** direction `R` is **not** a block-size effect: §5.2 is already the **one-link**
(S2) model and gives `R=−0.15` (`β=0`, radial). So `R` is governed by the **direction type**
(`∇·v`), not by the number of free links — the two-link block merely adds the coupled `A–V` chart,
which (for linear directions) costs nothing.

---

## 6. Honest notes

* **What is theorem vs. numerics.** `R=0` for linear `v` [T] and the acid-test characterization
  [T] are proved in §2; the tables are [V] corroboration. The `β=2.5` uniform-sampling rows carry a
  heavy-tail artefact (§4 note) that Haar sampling removes; the exact theorem does not need either.
* **Two charts, two conventions.** The absolute value of `W_geo` depends on the flat-chart / mode
  normalization (the S2 §4 caveat). What is convention-independent here is (i) the **identity**
  `G″=2W_geo` (`R=0`) for linear directions, and (ii) the **sign/vanishing** of `G″(β=0)` — the acid
  test. Only (i)–(ii) are used.
* **Why `W_geo` can be negative for the covariant direction without contradicting `W_geo≥0`.** The
  theorem `W_geo≥0` (ch05 `sec:wgeo`) is the Prékopa statement for a **linear** marginal
  `ξ=⟨v,log U⟩`, `v` constant. A field-dependent `v` is outside its hypothesis; `W_geo(v)<0` there is
  consistent (and is what fails the acid test).
* **The covariant tie to S2 §6.** The failure of the clean reduction is governed by the `∇·v` terms
  `E[(∂_v𝒱)(∇·v)]` (equivalently the cross term `X=−2Cov(G,Φ′)` of S2) — the “interior
  coarse-projection boundary term” S2 §6 left open. A **covariant block-spin** (`P` field-dependent)
  is exactly a `∇·v≠0` direction; §5 models it minimally.
* **Scope.** SU(2), fixed boundary, full fibre integrated (no interior coarse constraint beyond the
  field-dependent direction). SU(N) generalizes by the same edge-vanishing argument (Weyl-chamber
  flat chart, Vandermonde Jacobian → 0 on the walls). Not attempted here.

---

## 7. Status ledger

| item | statement | label |
|---|---|---|
| two-link block, coupled `A–V`, measure factorizes in `(W=AV,V)`; one-plaquette values reproduced | §1, §3.2 | **[T]**+**[V]** |
| generalized flat-chart Stein IBP-1/IBP-2 (§2.2) for any bounded field (vanishing edge flux), dimension-independent | §2.2 | **[T]** |
| `R(v)=G″−2W_geo=0` for every **constant** linear direction, any block, any boundary, any β | §2.3, §4 | **[T]** (+**[V]**) |
| `(★)` in boundary-conditioned two-link fibre ⟺ `2W_geo(v|U′)≥δ_W`, same `δ_W` as S2 | §2.3 | **[T]** |
| acid test `G″(β=0)=0` holds for constant `v` (`∇·v=0`, `E_0[Φ′_v]=0`); fails generically for field-dependent `v` | §2.4, §4 | **[T]** (+**[V]**) |
| field-dependent direction: `R≠0`, `O(0.03–0.6)`; bounded radial `W_geo(β=0)>0` (false positive), right-invariant `W_geo(β=0)<0` | §5 | **[T]**+**[V]** |
| `R=0` controlled by block size / boundary (linear) | **no**: identically `0`, independent of both | **[T]** |
| covariant / field-dependent `G″` as a mass-gap tool | **fails the acid test** ⇒ not a valid tool | **[NG]** |
| SU(N) two-link | same edge-vanishing argument; not computed | **[O]** |

**Sharpest line.** The exact reduction `G″=2W_geo` is not a one-plaquette accident and does not degrade
with block size: for every **constant** (linear-P) direction in a fixed-boundary two-free-link block
it holds exactly, `R=0`, by the same IBP that proved S2 — the block’s extra free link and the `A–V`
coupling cost nothing because the edge-vanishing lemma is dimension-independent. The remainder `R`
is nonzero (generically) for a **field-dependent** (covariant) direction, through `∇·v`; and the β=0
acid test `G″=0` forces the direction into the bounded-flux, `E_0[Φ′_v]=0` class, which every
constant `v` supplies and a generic field does not. So the linear-P class is the one that is
simultaneously exact and free-field-safe for every `β,θ`. The interior coarse-projection term S2 §6
flagged is thus the covariant (`∇·v≠0`) direction, and it is automatically excluded by the acid test.
The UV `(★)` input is obstructed by neither boundary holonomy nor block size — only by the
requirement that the coarse projection be **linear**, which the campaign’s `⟨Pᵀ∇S,v⟩` already is.
Clay 0%.

---

*Refs.* `S2-BOUNDARY-COUNTEREXAMPLE.md` (one-plaquette model, Lemmas 1–2, corollary); `ch05-theorems.tex`
`sec:wgeo` (`W_geo≥0`, Prékopa/linear-marginal), `eq:split`, `eq:vardecomp`; `STAR-WAVE-SYNTHESIS.md` §2;
witness `research/v3/w7/g2_two_plaquette.py`.
