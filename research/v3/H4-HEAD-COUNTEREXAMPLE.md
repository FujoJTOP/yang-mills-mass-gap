# H4 — the corrected head inequality at the coarse minimiser: **no typical-field counterexample exists**, and none can (the sign is an unconditional theorem `[T]`); the head's real content is the uniform margin `δ_W` = H1 `[O]`

> 2026-11-24 · Attack H4 on the **head** of the UV tower — the last open piece of `(★)` after
> `STAR2-CLASSIFICATION.md`. The V3/G4 counterexample (`G″=A−Bβ²<0` for `β>4.307`) sits at the **barrier**
> `a₀=0`; the corrected `(★_loc)` evaluates at the **coarse minimiser** `a₀*`. This note asks whether a
> counterexample survives there, in a **larger block** (SU(2) 1×2 and 1×3, coupled `A–V` chart, `X≠0`).
> Reads: `STAR2-CLASSIFICATION.md`, `V3-COUNTEREXAMPLE-SCOPE.md`, `G2-TWO-PLAQUETTE.md`, `G4-GENERAL-BOUND.md`,
> `docs/paper-v2/ch04-gap.tex`, `ch05-theorems.tex` (`eq:split`, `eq:vardecomp`).
> Witness: `research/v3/w7/h4_head_counterexample.py` (`py -3.12 [quick|full]`; deterministic 1×1 GL
> quadrature + Haar-proposal importance sampling for the blocks; output quoted below).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[V]** numerical (deterministic
> quadrature / importance sampling; validations and residuals quoted) · **[NG]** no-go for a named class ·
> **[O]** exact missing step.

---

## 0 · Verdict (no QED; Clay 0%)

**No typical-field counterexample at the coarse minimiser — and, structurally, none can exist.**

1. **[T] The sign of the corrected `(★_loc)` is an unconditional theorem.** For the Haar measure on any
   block, the coarse minimiser `a₀*=argmax ρ_v` is **interior** (the marginal vanishes at the chart edge
   `|ξ_i|=2π` because `J→0`; §1, verified numerically `ρ_v(2π)~5·10⁻¹⁴`). At an interior maximum
   `ρ_v'(a₀*)=0`, `ρ_v''(a₀*)≤0`, hence
   ```
     G″(v;a₀*) = −d²/da² log ρ_v|_{a₀*} = (ρ_v'/ρ_v)² − ρ_v''/ρ_v = −ρ_v''/ρ_v ≥ 0 .            (0.1)
   ```
   for every block, boundary, `β` and linear direction `v`, with `X` included automatically. `G″>0`
   **strictly** iff the maximiser is non-degenerate (`ρ_v''(a₀*)<0`); a quartic-flat top (`ρ_v''=0`) would
   give `G″(a₀*)=0`. So the sign is `≥0` **unconditionally** and `>0` in every non-degenerate case (all
   scanned configs). The degenerate `=0` case is not a harmful counterexample — it is exactly the `δ→0`
   boundary of the uniform bound, i.e. H1. Thus the corrected `(★_loc)` **cannot be refuted with a negative
   sign**; the only open content is the *uniform lower bound* `δ_W`, which is the pre-existing wall
   **H1 `[O]`**.
2. **[V] The value of the minimum, and its scaling.** In blocks 1×1 / 1×2 / 1×3 (SU(2), `β=2.5`,
   boundaries `θ∈[0,π]`, longitudinal/`k=0` and single-link modes), the measured `G″(a₀*)` is **O(1),
   positive, and does not soften with block size**; its **global minimum over the scan is at small `β`**,
   where it stays just above the Haar floor `N/12=0.1667`. Headline numbers (`[V]`, witness):
   - `β=2.5, θ=π` k=0 mode `D3=(e₃,e₃)/√2`: `G″(a₀*)=+1.022` (vs one-plaquette barrier value `+0.136` at
     `a₀=0`, and the V3 hairline `2W_geo`);
   - modes `A3, V3, Diff3` at the same point: `+0.350, +0.692, +0.335`, all `> N/12`;
   - boundary scan (k=0, β=2.5): `G″(a₀*)=+1.170,+1.151,+1.142,+1.054,+0.981` at `θ=0,π/4,π/2,3π/4,π`
     (worst at maximal holonomy `θ=π`), non-abelian boundary `+1.058`;
   - block size (k=0, β=2.5, θ=π): `L=1:+0.696 (1×1 exact)`, `L=2:+1.012`, `L=3:+0.628`;
   - `β`-scan (k=0, 1×2, θ=π): `β=0:+0.240, 1:+0.484, 2.5:+0.989, 4:+1.568, 5:+2.014` — increasing; the
     **minimum is at `β=0`** and `≈0.24 > N/12=0.167` (one-plaquette `β=0`: `0.205 > N/12`).
3. **[V] Negative-curvature shoulders exist but are never at the minimiser.** The marginal-curvature
   `G″(a)` has genuine negative shoulders at `β=2.5` already in the **one-plaquette** model (`a∈[0.7,3.0]`,
   min `−0.077` at `a=2.0`, exact), and the 1×2 coupling deepens them (bulk min `≈−0.33` at `a≈−3.1` for
   `D3`; §3.5). They sit on the **flanks / near the barrier**, never on the marginal's maximum `a₀*`. This is
   exactly the V3 mechanism (`−Var(G|a₀)` exposed off the well); the corrected `(★_loc)` avoids it by
   construction.
4. **Acid tests pass.** Free massless Gaussian (`Φ≡0`, Gaussian `S`): `G″=κ→0` at every `a₀` (no
   manufactured floor). Pure Haar `β=0`: `G″(a₀*)≈0.24 > N/12` (the Haar floor, *not* the free field).

**Bottom line.** The "falsify the head at the minimiser" programme returns a **negative result**: the sign
is a theorem, the value stays `O(1)>N/12` across the scanned `(block, boundary, mode, β)`, and the only
surviving inequality at the head is the **uniform** `δ_W>0`, which is H1. The wall did not move; H4 closes
the question "is there a typical-field counterexample at `a₀*`?" with a rigorous **no** and a numeric map.

---

## 1 · The object, the model, and why the minimiser is interior

**Object.** For a unit **linear** (constant, `∇·v=0`) direction `v` in the fibre tangent space and coarse
coordinate `a=⟨v,ξ⟩`, the true coarse curvature (ch05 `eq:split`/`eq:vardecomp`, G4 Lemma G1) is
```
  G″(v;a₀) := −d²/da² log ρ_v(a)|_{a₀} = E_{μ}[∂_v²V | a₀] − Var_{μ}(∂_vV | a₀),      V = S + Φ,
  = Δ_dyn + W_geo + X,   Δ_dyn=E[H]−Var(G),  W_geo=E[Φ″]−Var(Φ′),  X=−2Cov(G,Φ′) ,      (1.1)
```
all expectations **conditioned at `a=a₀`**; `G=∂_vS`, `H=∂_v²S`, `Φ′=∂_vΦ`, `Φ″=∂_v²Φ`. The identity (1.1)
is pure differentiation (Lemma G1, `[T]`), model-independent, and includes `X`. The **coarse minimiser** is
`a₀*=argmax_a ρ_v(a)` (the mode of the marginal = minimiser of the coarse effective action `S_eff=−log ρ_v`).

**Model (SU(2), flat exponential chart `ξ∈ball(2π)` per link, Haar `dU=J(ρ)d³ξ`, `J=sin²(ρ/2)/ρ²`,
`Φ=−log J`).** Free interior links; fixed boundary links enter as fixed SU(2) matrices.
```
  1×1:  free U ;  P = U B .                                                   (reproduces V3/G4 exactly)
  1×2:  free A,V ;  P₁ = A V ,  P₂ = B V⁻¹ .                                   (G2 two-plaquette block)
  1×3:  free A,V,W ; P₁ = A V , P₂ = B₁ V⁻¹ W , P₃ = B₂ W⁻¹ .                    (3-plaquette chain)
  S = β Σ_□ (1 − ½tr P_□) ,   dμ ∝ e^{−S} Π dU = e^{−(S+Φ)} dξ ,   B=e^{i(θ·σ)/2}.
```
The **k=0 coarse mode** is the block-average along the longitudinal axis: `v=(e₃,…,e₃)/√L` over the `L`
free links (`L=1,2,3`). Longitudinal (`e₃`) single-link modes `A3,V3` and the difference `Diff3` are also
scanned. The two-plaquette chart is genuinely coupled (`ξ_A,ξ_V` enter `W=AV` and `BV⁻¹`), so `X≠0` at
`a₀*`; we verify the decomposition (1.1) term-by-term.

**Why `a₀*` is interior `[T]`.** `a=⟨v,ξ⟩` ranges over `[−Σ_i2π|v_i|, +Σ_i2π|v_i|]`. The endpoint `a_max`
is attained only at `ξ_i=2π v̂_i` for all `i`, where `J(ρ_i)=sin²(ρ_i/2)/ρ_i²→0`; hence `ρ_v(a_max)=0`
(and `=0` at `a_min`). `ρ_v` is continuous, positive on the open interval, and real-analytic (`V` analytic,
pushforward of an analytic density). A positive density vanishing at both endpoints attains its **maximum
strictly inside**, and a real-analytic non-constant density has a non-degenerate isolated max, so
`ρ_v''(a₀*)<0` unless the top is degenerate to all orders (a quartic-flat plateau, non-generic). Numerically
`ρ_v(2π)≈5·10⁻¹⁴` (one-link, GL quadrature), confirming the edge vanishing. ∎

Therefore (0.1) holds and `G″(a₀*)≥0` **unconditionally**, with `G″(a₀*)>0` in every non-degenerate
maximiser (a quartic-flat top would give `0`, the `δ→0` edge of H1). The *edge loophole* (a boundary
minimiser, where `ρ_v'≠0` and the one-sided curvature could be either sign) is **closed** by the Haar edge
`J→0`.

---

## 2 · Validations of the numerical machinery — `[V]`

Two independent exact references pin the estimators before any `β>0` block result is trusted.

**(a) One-plaquette exact reproduction (V3/G4).** Deterministic Gauss–Legendre on the exact slice marginal
reproduces V3 to all quoted digits:
```
  A = E_Haar[Φ'(r)/r] = 0.2051152978 ,  B = Var_Haar[sin(r/2)/r] = 0.0110592938 ,  β* = √(A/B) = 4.30661 .
  G''(e₃;θ=π,a₀=0) = A − B β² :   β=0  0.205115 | 2.5  0.135995 | 5  −0.071367 | 10  −0.900814 .
  coarse minimiser a₀* & G''(a₀*): β=2.5 a₀*=−2.15 G''=+0.696 | 5 a₀*=−2.575 +1.326 | 10 a₀*=−2.850 +2.597 .
```
**(b) End-to-end MC at `β=0` (pure Haar) vs the exact convolution.** At `β=0` the two free links are
independent Haar and `a` is a sum of two scaled one-link marginals, with the one-link marginal in closed
form `ρ₁(a) ∝ 2π∫_{|a|}^{2π} sin²(ρ/2)/ρ dρ` (the normalization constant cancels in `G″`). The
importance-sampling scanner (sampling, `a`-binning, exact pointwise Hessian moments) reproduces the exact
convolution:
```
  D3 :   MC a₀*=+0.163  G''=+0.2389   |   EXACT a₀*=0.000  G''=+0.2396     (0.3%)
  Diff3: MC a₀*=−0.015  G''=+0.2392   |   EXACT a₀*=0.000  G''=+0.2396
```
**(c) End-to-end MC at `β=2.5` vs the exact one-plaquette** (validates the `e^{−S}` reweighting at
interacting coupling): `G″(a₀*)_{MC}=+0.713` vs exact `+0.696` — the `0.017` gap is entirely the `a₀*` bin
location (`−2.231` vs `−2.150`); at `a=−2.231` the exact curve is `≈0.71`, so the MC matches the exact to
`<0.5%` at its own `a₀*`. `G″(0)_{MC}=+0.133` vs exact `+0.136`. Machinery validated.

---

## 3 · The scan — `[V]`

Per-`a` histogram of the **exact pointwise** Hessian moments; `G″` is read from the bin at the marginal
maximum (`Nb=600` bins; ESS per bin quoted; `Hb=250000`-sample streaming chunks). `N=6·10⁶` (S2),
`2·10⁶` (scans). **`G″(a)=−d²logρ/da²` is read directly at the dense max bin (no smoothing window), so the
reported `G″(a₀*)` is the bin-averaged true curvature at the coarse minimiser.**

### 3.1 Modes at the worst point `(β=2.5, θ=π)`, 1×2 block (`N=6·10⁶`)
| mode `v` | `a₀*` | `G″(a₀*)` | `Δ_dyn` | `W_geo` | `X` | bulk shoulder `G″` | `@a` |
|---|---|---|---|---|---|---|---|
| `A3`  (A only)      | −1.330 | **+0.350** | +0.139 | +0.221 | −0.010 | −0.21 | +3.2 |
| `V3`  (V only)      | +2.126 | **+0.692** | +0.463 | +0.249 | −0.019 | −0.32 | −2.5 |
| `D3`  (k=0, both)   | +0.311 | **+1.022** | +0.797 | +0.261 | −0.036 | −0.33 | −3.1 |
| `Diff3` (gradient)  | −2.473 | **+0.335** | +0.110 | +0.251 | −0.026 | −0.14 | +3.3 |
`G″(a₀*)>0` for every mode, all `> N/12=0.1667`; the decomposition `Δ_dyn+W_geo+X` matches `G″(a₀*)` to
the digits shown (identity (1.1) is exact pointwise). `X≠0` (the coupled `A–V` chart) but small
(`|X|≤0.04`). The negative shoulders are on the flanks (`|a−a₀*|≥2.4`), never at the minimiser.

### 3.2 Boundary scan, k=0 mode `D3`, `β=2.5` (`N=2·10⁶`)
| `θ` (longitudinal boundary) | `a₀*` | `G″(a₀*)` | `Δ_dyn` | `W_geo` | `X` | bulk shoulder `G″` | `@a` |
|---|---|---|---|---|---|---|---|
| 0        | +0.074 | **+1.170** | 0.971 | 0.198 | +0.001 | −0.55 | −4.3 |
| π/4      | +0.222 | **+1.151** | 0.951 | 0.202 | −0.002 | −0.47 | −4.5 |
| π/2      | +0.222 | **+1.142** | 0.938 | 0.212 | −0.008 | −0.51 | −3.3 |
| 3π/4     | +0.429 | **+1.054** | 0.844 | 0.233 | −0.024 | −0.28 | −2.9 |
| π        | +0.429 | **+0.981** | 0.760 | 0.261 | −0.040 | −0.33 | −3.1 |
| non-abelian `b=(π/2,0,π/2)` | +0.133 | **+1.058** | — | — | — | −0.30 | — |
`G″(a₀*)≥0.98` across the entire boundary range (worst at maximal holonomy `θ=π`); `X→0` at `θ=0` by
isotropy and grows with `θ`. All above the β=2.5 good-set modulus `ρ₀≈0.217` (`STAR2` §2).

### 3.3 Block-size scaling, k=0 mode, `β=2.5, θ=π` (`N=2·10⁶` for MC)
| block `L` | `a₀*` | `G″(a₀*)` | `X` |
|---|---|---|---|
| 1×1 (exact)  | −2.150 | **+0.696** | 0 (e₃ slice) |
| 1×2          | +0.281 | **+1.012** | −0.038 |
| 1×3          | +0.998 | **+0.628** | −0.044 |
`G″(a₀*)` is `O(1)` and **does not soften with block size** (0.63–1.01); no `1/L²` decay. (The `a₀*` is
block-dependent — the well shifts — but the curvature at it stays `O(1)`.)

### 3.4 `β`-scan, k=0 mode `D3`, `θ=π` (`N=2·10⁶`)
| `β` | `a₀*` | `G″(a₀*)` | bulk shoulder `G″` | `@a` |
|---|---|---|---|---|
| 0   | +0.133 | **+0.240** | +0.01 (none) | — |
| 1   | +0.400 | **+0.484** | −0.28 | −2.4 |
| 2.5 | +0.459 | **+0.989** | −0.33 | −3.1 |
| 4   | +0.429 | **+1.568** | −0.25 | −1.9 |
| 5   | +0.311 | **+2.014** | −0.22 | −1.9 |
`G″(a₀*)` increases monotonically with `β`; the **minimum is at `β=0`** (`≈0.24 > N/12`), where the pure
Haar measure is log-concave and there is **no** negative shoulder at all. The `β=0` value `0.240` matches
the exact convolution `0.2396` (§2b). (`β=5` resolves with the adaptive ESS floor; it agrees with the exact
one-plaquette growth `+1.326,+2.597` at `β=5,10`, §2a.)

### 3.5 The `G″(a)` curve and the negative shoulders (D3, β=2.5, θ=π, `N=6·10⁶`)
The marginal-curvature curve `G″(a)` rises to a **central well maximum `≈+1.05` near `a=0`** and has
**two negative shoulders** on the flanks, both far from the mode `a₀*=+0.31`; the well value is
`G″(a₀*)=+1.02`:
```
  a=−3.15 G″=−0.33 | a=−2.27 −0.17 | a=−1.38 +0.44 | a=−0.49 +0.96 | a=−0.04 +1.05 | a=+0.40 +0.99
  a=+0.84 +0.83 | a=+1.73 +0.34 | a=+2.18 +0.11 | a=+2.62 −0.005 | a=+3.07 −0.05 | a=+3.51 −0.05
```
*(The single deepest raw bin, `−4.7` at `a=−4.2`, is a high-`|a|` tail artifact — neighbouring curve
points there are `−0.04` and `+0.18`; the honest **bulk** shoulder depth is `≈−0.33` at `a≈−3.15`, quoted
above.)* The **exact one-plaquette** at the same `β=2.5` already shows a (shallower) negative shoulder
`a∈[0.7,3.0]`, min `−0.077` at `a=2.0` (§2a), confirming the mechanism is a real one-plaquette feature
that the block coupling deepens. The minimiser never sits on a shoulder.

---

## 4 · Acid tests and the free-field limit — `[T]+[V]`

1. **Free massless field** (`Φ≡0`, `V=κ|ξ|²/2`): the marginal of `a=⟨v,ξ⟩` is Gaussian `∝e^{−κa²/2}`, so
   `G″(a₀)=κ` at *every* `a₀`, and `κ→0` as the mass →0. **The tool manufactures no floor** — a positive
   `δ_W` requires the Wilson interaction `β>0`. `[T]`
2. **Pure Haar `β=0`**: `G″(a₀*)≈0.24` (D3, 1×2) and `0.205` (1×1), both `> N/12=0.1667`, matching the
   `β=0` theorem `G″≥N/12` (`ch05`). Note this is the **Haar floor**, not the free field — the free-field
   failure is item 1. `[T]+[V]`
3. **Continuity / monotonicity in `β`.** `G″(a₀*)` rises monotonically with `β`
   (`0.24→0.48→0.99→1.57→2.01` at `β=0,1,2.5,4,5`), so the head margin is *weakest at small `β`*, exactly
   where it approaches (but stays above) the Haar floor; at `β=0` there is no negative shoulder at all
   (`#neg=0`). `[V]`

---

## 5 · Verdict and status ledger

| item | statement | label |
|---|---|---|
| marginal-curvature identity `G″(v;a₀)=−d²/da²logρ_v\|_{a₀}=E[∂_v²V\|a₀]−Var(∂_vV\|a₀)`, any block, no IBP | Lemma G1, exact pointwise | **[T]** |
| coarse minimiser `a₀*=argmax ρ_v` is **interior** for the Haar measure (`ρ_v→0` at `\|ξ_i\|=2π`, `J→0`) | §1; `ρ_v(2π)≈5·10⁻¹⁴` | **[T]** |
| `G″(v;a₀*)=−ρ_v''/ρ_v ≥ 0` at every block/boundary/β/linear-`v`, `X` included | (0.1) | **[T]** |
| strict `G″(a₀*)>0` iff non-degenerate max; `=0` only for a quartic-flat top (= `δ→0` edge of H1) | §1 | **[T]** |
| **a typical-field counterexample `G″(a₀*)<0` at the coarse minimiser** | **DOES NOT EXIST** (edge loophole closed by `J→0`) | **[T]** (+**[V]**) |
| exact one-plaquette `A=0.20511530, B=0.01105929, β*=4.30661`; `G″(a₀*)=+0.696,1.326,2.597` at β=2.5,5,10 | §2a, reproduces V3/G4 | **[T]**+**[V]** |
| `β=0` MC vs exact convolution (`0.2389` vs `0.2396`); 1×1 MC vs exact (`0.713` vs `0.696`@its `a₀*`, `<0.5%`; `0.133` vs `0.136`) | §2b,c machinery validations | **[V]** |
| 1×2 modes @β=2.5,θ=π: `G″(a₀*)={0.350,0.692,1.022,0.335}`, all `> N/12`; `X≠0` but small | §3.1 | **[V]** |
| boundary scan: `G″(a₀*)∈[0.98,1.17]`, worst at maximal holonomy `θ=π`, non-abelian `+1.06` | §3.2 | **[V]** |
| block size `L=1,2,3`: `G″(a₀*)=0.696, 1.012, 0.628` — no softening, no `1/L²` | §3.3 | **[V]** |
| `β`-scan: `G″(a₀*)=0.240→2.014` (β=0→5); minimum at β=0 (`>N/12`); `β=0` no negative shoulder | §3.4 | **[V]** |
| negative shoulders on the **flanks**, never at `a₀*` (exact 1×1 min `−0.077`; coupled bulk `≈−0.33`) | §3.5, §2a | **[V]** |
| acid: free Gaussian `G″=κ→0`; pure Haar β=0 `G″≈0.24>N/12` | §4 | **[T]**+**[V]** |
| **uniform** `δ_W>0` at `a₀*` (uniform in block, boundary, β, mode) | **H1 / coarse spectral gap** — the surviving open step | **[O]** |

**Reading.** The corrected `(★_loc)` at the coarse minimiser is **true and unconditionally signed**: a
negative counterexample cannot exist there. The numeric scan shows the head value is `O(1)` (minimum
`≈0.20–0.24 > N/12`) and stable across block size and boundary, with negative curvature confined to the
flanks. What remains open is *only* the **uniform lower bound** `δ_W` — the same coarse spectral gap H1 the
campaign already identified. The "falsify the head at the minimiser" programme therefore **fails to
produce a counterexample, and is closed by a structural argument**, leaving the wall exactly where it was.

---

## 6 · Witness — `research/v3/w7/h4_head_counterexample.py`

`py -3.12 h4_head_counterexample.py full` (exit 0). Deterministic Gauss–Legendre for the exact 1×1
references (`n=6000` radial nodes); Haar-proposal importance sampling with streaming per-`a` accumulation
(`Nb=600` bins; ESS reported per bin) for the 1×2/1×3 blocks. The pointwise Hessian `∂_v²V` is built from
the analytic `Φ',Φ''` (stable `1/x−cot x`, `½(−x⁻²+csc²x)`) plus a central finite-difference of the smooth
`S` along the constant flow (`δ=10⁻⁴`). Validations §2. No third-party sources downloaded.

---

## 7 · Sharpest line

**At the coarse minimiser the corrected head inequality `(★_loc)` is not merely surviving — it is an
unconditional theorem: the coarse minimiser is interior (the Haar marginal vanishes at the chart edge), so
`G″(a₀*)=−ρ_v''/ρ_v≥0` (strictly `>0` for a non-degenerate max) for every block, boundary, `β` and linear
mode, with `X` included. The larger
1×2/1×3 blocks change nothing qualitatively: `G″(a₀*)` stays `O(1)>N/12` (minimum `≈0.20–0.24` at `β=0`;
`+1.02` for the k=0 mode at `β=2.5,θ=π`), while genuine negative shoulders (`−0.08` exact one-plaquette,
`≈−0.33` bulk in the coupled 2-plaquette block) sit on the flanks, never on the minimiser.** So there is
**no typical-field counterexample at `a₀*`** `[T]+[V]`; the head's only open inequality is the **uniform**
margin `δ_W>0`, which is H1 `[O]`. **No QED. Clay 0%.**

---

*Refs.* `V3-COUNTEREXAMPLE-SCOPE.md` §0,§2.4,§4.2 (`(★_loc)`, the barrier counterexample, the minimiser);
`G2-TWO-PLAQUETTE.md` §1 (the 1×2 block, `A–V` coupling); `G4-GENERAL-BOUND.md` §1–3 (objects (L)/(U)/(C),
Lemma G1); `STAR2-CLASSIFICATION.md` §2–3 (good-set theorem, head=H1); `docs/paper-v2/ch04-gap.tex`
(G2.5 open cell = H1), `ch05-theorems.tex` `eq:split`/`eq:vardecomp`. Witness
`research/v3/w7/h4_head_counterexample.py`. No third-party sources downloaded.
