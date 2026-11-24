# D6-DECISION — decompactification of the Yang–Mills constant (soft) mode: **decay to 0 (Route A) vs saturation at an L-independent floor (Route B)** — adversarial decision attempt

> **2026-11-24 · Wave D6**, the decision note the UV-face campaign was waiting for.
> Question: at `β=2.5`, does the constant-mode curvature `G″(a₀*)` **decay to 0** as `L→∞`
> (then the uniform floor is `H1`, and any small-volume positivitiy is a finite-size mirage), or
> **saturate at a positive, L-independent value** (the Haar / Prékopa floor genuinely survives)?
> Reads: `SM1-SOFT-EXACT.md`, `SM4-CONSTANT-MODE.md`, `SM3-ANISOTROPIC-RG.md`, `SM2-RANK-ONE.md`,
> `H2-DEFECT-BOUND.md`, `docs/paper-v3/main.tex` §6.7.
> Witness: `research/v3/w7/d6_decision.py` (`py -3.12 -I`, numpy 2.4.6 + scipy 1.18.1; exact
> Gauss–Legendre quadrature + FFT convolution; output quoted §7). **No third-party sources.**
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

**Route B (saturation) survives; Route A (decay) is dead as an asymptotic of the models the data come
from, and the observed `≈1/L` is a finite-size transient of a quantity that saturates.** Three
independent legs:

1. **A rigorous `[T]` identity turns the question into a susceptibility.** For a Gaussian-limit
   marginal, `G″(a₀*) = 1/Var(a)`, and with the unit-norm constant mode
   `Var(a) = χ_L = Σ_d (1−|d|/L) C(d) → χ := Σ_d C(d)` (the zero-momentum susceptibility). Hence
   ```
   χ < ∞  (finite constant-mode correlation length)  ⟹  G″ → 1/χ > 0   [SATURATION, Route B]
   χ = ∞  (gapless / critical)                       ⟹  G″ → 0          [DECAY, Route A]
   ```
   with the exact decay exponent set by the growth of `χ_L`: `χ_L ≍ L^p ⇒ G″ ≍ L^{−p}`
   (`p=2` free box-gap `4π²/L²`; `p=1` marginal/`C(d)~1/d`; `p=0` saturation). **This is the whole
   dichotomy in one line.**
2. **At `β=0` it is saturation, exactly `[T]+[V]`.** The pure-Haar marginal is the `L`-fold convolution
   of the one-link marginal; by the CLT it converges to a Gaussian of variance `χ=Var_Haar(ξ₃)=3.7198`,
   so `G″ → 1/χ = 0.2688`, **L-independent**. Exact convolution (`L=1…16`): `0.205 → 0.265 → 0.2688`;
   it does **not** decay (`[D1,D2]`). The Haar floor is the L-independent band `[A, 1/χ]=[0.205,0.269]`.
3. **For the `1×L` chain the data live in, Route A is a `[NG]`.** A 1D finite-range Gibbs measure has
   exponential clustering (unique transfer-matrix state, no transition), so `Σ_d C(d) < ∞`, hence
   `χ_L` bounded and `G″ → 1/χ > 0`. **The `1×L` chain cannot decay as `1/L→0`** — the measured
   `1.01, 0.65, 0.42` (`L=2,3,4`) is the finite-`L` transient; a `b + a/L²` fit gives `b=+0.266 ≈` the
   Haar value `0.269`, while a pure `b + a/L` decay forces a **negative** (unphysical) intercept `−0.15`
   (`[D3]`).

**What is `[O]` (and why this is not a free lunch).** The constant-mode curvature is the inverse
zero-momentum susceptibility `1/χ`, **not** the pole mass. A mass gap `m>0` *implies* `χ<∞` *implies*
saturation — so **saturation is a necessary condition the program must pass** — but `1/χ` is not by
itself `m²`. For the **`1×L` chain** the saturation is a theorem `[T]`; for **4D at `β=2.5`** the
statement `χ<∞` **is** `H1`/the mass gap, unchanged `[O]`. The one inequality Route B needs is *not*
global log-concavity — that is **false** pointwise (`[NG]`, §4, `[D4]`: the Savvidy-indefinite flanks
locally dip below the Haar curve) — but the weaker **mode-value / finite-susceptibility** bound (§5).

**One line.** *Decay (Route A) is falsified for the `1×L` chain by 1D exponential mixing and disfavoured
by the `β=2.5` finite-size fit (`1/L² → 0.27`, not `1/L → 0`); saturation (Route B) is exact at
`β=0` (`1/χ=0.269`, L-independent) and is the only behaviour consistent with a finite constant-mode
susceptibility `χ`. The deciding hypothesis is `χ<∞` — a theorem in 1D, and in 4D exactly `H1`. The
`≈1/L` trend is the transient, not the asymptote.* **No QED. Clay 0%.**

---

## 1 · The object, and why the dichotomy is exact

**The curvature.** `G″(v;a₀) = −d²/da² log ρ_v(a)|_{a₀} = E[∂_v²V|a₀] − Var(∂_vV|a₀)` (`V=S_W+Φ`,
`a=⟨v,ξ⟩`, flat exponential chart, `Φ=−log J` the Haar Jacobian `[T]`, SM1 (1.1)/H4). `a₀*` is the
coarse minimiser, so `G″(a₀*) = −ρ_v''/ρ_v ≥ 0` `[T]`.

**The Gaussian/CLT marginal.** For a *delocalised* mode (the `k=0` mode `v=(e₃,…,e₃)/√L` is a spread
weighted sum), `ρ_v(a)` is asymptotically Gaussian by the CLT, and then
```
   G″(a₀*) ≃ 1 / Var(a) ,   Var(a) = Σ_{l,l'} v_l v_{l'} C(l,l') = Σ_d (1−|d|/L) C(d) =: χ_L ,   (1.1)
```
where `C(l,l') = ⟨ξ₃,l ξ₃,l'⟩_c` is the connected link-field correlation. `[T]` (Gaussian marginal) /
`[C]` (CLT needs enough delocalisation: equal weights `v_l=1/√L` satisfy Lindeberg for bounded support).

**The dichotomy.** `χ_L = Σ_d(1−|d|/L)C(d)`:
- `C(d)` **integrable** (`ξ_corr<∞`) ⟹ `χ_L = χ + O(1/L)`, `χ=Σ_d C(d)<∞` ⟹ `G″→1/χ>0`: **saturation**.
- `C(d)~|d|` (free zero mode) ⟹ `χ_L≍L²` ⟹ `G″≍4π²/L²→0`: **decay, box-gap exponent 2**.
- `C(d)~1/|d|` (marginal / critical) ⟹ `χ_L≍L` ⟹ `G″≍c/L→0`: **decay, exponent 1** (the `≈1/L`).

So "decay vs saturation" **is** "is the constant-mode zero-momentum susceptibility `χ` finite or
infinite" **is** "is the constant-mode channel gapped." There is no other content — and that is exactly
why the question is hard (it touches the mass gap) and exactly why it is decidable in the models.

---

## 2 · Route A (decay): statement, exact scaling, hypothesis `[C]/[NG]`

**Statement A.** If the constant mode is gapless, `Var(a)=χ_L→∞` as `L→∞`, then `G″(a₀*)→0`.

**Exact scaling `[T]` (given a Gaussian-limit marginal).** `G″(a₀*) = 1/χ_L` and
```
   χ_L = Σ_{d} (1−|d|/L) C(d) .
```
Hence the decay rate is the growth rate of `χ_L`:
| mechanism | `C(d)` | `χ_L` | `G″` | exponent |
|---|---|---|---|---|
| free massless zero mode (box gap) | `~|d|` | `≍L²` | `4π²β/L²` | **2** |
| critical / marginal (1D at `T_c`) | `~1/|d|` | `≍L` | `c/L` | **1** (`≈1/L`) |
| any finite `ξ_corr` | decays exp. | `→χ<∞` | `→1/χ` | **0** (saturation) |

**Hypothesis (A):** `(A)` *`χ_L→∞`* — equivalently the constant mode's correlation length diverges
(`ξ_corr∝L`) or the theory is gapless in this channel. **Label:** conditional `[C]`; and for the
`1×L` chain it is **`[NG]`** (§4.1: 1D exponential mixing forbids it). For the free massless field it
is realised `[T]+[V]` (box-gap, exponent 2), which is the **acid test**, not the interacting head.

**Why Route A is wrong for the actual data — two precise reasons.**
- **Wrong exponent.** The free/box-gap `4π²/L² = 9.87, 4.39, 2.47` is ~10× the measured
  `1.01, 0.65, 0.42` and drops ×4 per doubling; the data drop ≈×1.6 (`[D5]`). The constant mode is
  **not** the free box-gap mode.
- **Wrong extrapolation.** A pure `1/L→0` fit to the three points needs an intercept
  `b=−0.149` (`[D3]`): i.e. `1/L` decay does **not** head to `0` monotonically — it over-shoots below
  zero, the signature that the true asymptote is a *positive* floor with a faster transient, not `0`.

---

## 3 · Route B (saturation): statement, the `β=0` floor `[T]`, the needed inequality

**Statement B.** If `χ_L` is bounded uniformly in `L` (finite constant-mode susceptibility /
correlation length), then `G″(a₀*) ≥ 1/χ_max > 0` for all `L` — **saturation at an L-independent floor.**

### 3.1 The `β=0` Haar floor is rigorous and L-independent `[T]+[V]`

At `β=0` the links are i.i.d. Haar, so `C(l,l')=δ_{l,l'}·Var_Haar(ξ₃)` and the marginal of the unit-norm
`k=0` mode is the `L`-fold (weighted) convolution of the one-link marginal. By the CLT it converges to
`N(0,χ)` with the **exact closed form**
```
   χ = Var_Haar(ξ₃) = (1/3)E_Haar[ρ²] = (4π²−6)/9 = 3.71981 ,      G″ → 1/χ = 9/(4π²−6) = 0.26883 .
```
(`E_Haar[ρ²]=2(2π²/3−1)` for the SU(2) exponential chart; this equals `9/(4π²−6)` — the same value the
companion witnesses `d2_decompact.py`/`d4_haar_wilson.py` quote, an independent cross-check.)
Exact FFT convolution (`L=1…16`, `[D1,D2]`):
```
   L:        1        2        3        4        6        8       12       16
   G″(a₀*): 0.210    0.240    0.249    0.254    0.259    0.262    0.264    0.265   →  0.2688
   kurtosis: 2.57 → 2.78 → 2.86 → 2.89 → 2.93 → 2.95 → 2.96 → 2.97    →  3   (Gaussian)
   Var(a):  3.7198   3.7198   3.7198   3.7198   3.7198   3.7198   3.7198   3.7196  (L-independent)
```
So the pure-Haar curvature is the **L-independent band `[A, 1/χ] = [0.205, 0.269]`** (the `L=1` exact
value is the closed form `A=0.20511530`; the grid reading `0.210` sits `0.005` above it), rising
monotonically to `0.2688` and **not decaying**. This is
the `β=0` Haar floor, `[T]` (CLT) + `[V]` (exact). It is **not** a free field: the free massless field
has no `Φ` and gives `4π²/L²→0` (`[D5]`); the discriminator is precisely the Haar Jacobian `Φ`.

### 3.2 The needed inequality for `β>0`

The `β>0` marginal is `ρ_v^{(β)}(a) = ρ_v^{(0)}(a)·w(a)/Z` with `w(a)=E_Haar[e^{−S_W}|⟨v,ξ⟩=a]`. Then
```
   G″^{(β)}(a) − G″^{(0)}(a) = −(log w)''(a) = −[ (log ρ^{(β)})'' − (log ρ^{(0)})'' ] ,
```
so the pointwise statement **`G″^{(β)}≥G″^{(0)}` ⟺ `w` log-concave**. Two forms of the needed inequality:

- **(LC-W) [strong, pointwise]:** `w(a)` is log-concave (`(log w)''≤0`). Then `G″^{(β)}(a)≥G″^{(0)}(a)≥A`
  for **every** `a`, and Route B is immediate. **This is FALSE** at `β=2.5` (§4.2, `[D4]`): the
  Savvidy-indefinite fibre makes `w` locally *log-convex* on the flanks.
- **(B-ineq) [weak, mode-value / susceptibility] — the one Route B actually needs:**
  ```
     Var(a) = χ_L ≤ χ_max < ∞  uniformly in L      (equivalently  G″(a₀*) ≥ 1/χ_max > 0),
     sufficient if |C(l,l')| ≤ C₀ e^{−|l−l'|/ξ}  with ξ<∞  (exponential clustering).
  ```
  **Label:** `[C]` in 4D (`χ<∞` = `H1`), **`[T]`** in the `1×L` chain and at `β=0`.

The honest lesson of `[D4]`: Route B **cannot** be run through global Prékopa/Brascamp–Lieb, because the
Wilson action is not convex (Savvidy). It must be run through the **mode-value** bound `(B-ineq)`, which
is exactly the corpus's "measure-carried, `D`-defect" positivity — the `D`-defect is what breaks (LC-W)
pointwise while `(B-ineq)` at the mode still holds.

---

## 4 · Falsification attempt against the small-volume data (SM1/SM4)

### 4.1 `1×L` chain: Route A is a no-go `[NG]`

A 1D finite-range Gibbs measure at finite `β` has a **unique** transfer-matrix ground state and
exponential clustering `|C(l,l')|≤C₀e^{−|l−l'|/ξ}`, `ξ<∞` (standard; there is no phase transition in
1D with finite-range interactions). Then `Σ_d C(d)<∞`, `χ_L` bounded, and `G″→1/χ>0`. **Route A
(`1/L→0`) is therefore impossible for the `1×L` chain** — the model that produced the `1.01,0.65,0.42`
numbers. The observed decrease is the finite-`L` crossover from the UV `O(1)` regime into the plateau.
This is the strongest falsification and it is `[T]` (given the standard 1D mixing).

### 4.2 (LC-W) is false pointwise, true at the mode `[V]` (`1×1` exact, `β=2.5`, `θ=π`) `[D4]`

Exact `1×1` marginal, both constant-mode classes:
```
   class            -(log w)'' min      median      frac<-0.02     at β-mode a₀*:  G″^β   G″^0   Δ
   Cartan  e₃        −0.369            −0.069       0.53          Δ=+0.431 (a₀*=−2.20)
   Savvidy e₁        −0.113            −0.017       0.49          Δ=+0.157 (a₀*= 0.00)
```
So **(LC-W) fails on ~half the support** (the marginal is locally *less* curved than Haar at
`a=0`, `Δ=−0.069` for `e₃`) — the Savvidy-indefinite flanks produce genuine log-convexity. **But at the
β-mode `a₀*` the curvature is *above* the Haar value by `+0.43` (Cartan) / `+0.16` (Savvidy)** — the
floor holds *where it is evaluated*. This is the precise, non-obvious content: **positivity is a
mode-value fact, not a pointwise one** — matching H2/SM2's "BL-carried, `D`-defect" reading. It also
tells you *where* Route B is fragile: the Savvidy flank `Δ=−0.11` is the first place a large-`L`
reversal would appear.

### 4.3 The finite-size fit against the data `[V]` (`[D3]`)

`β=2.5`, `k=0`, `θ=π`: `(L,G″) = (2,1.01),(3,0.65),(4,0.42)` (MC scatter `0.05–0.15`).
```
   model                 fit                 pred L=2,3,4            rms     asymptote
   Route A   a/L         a=1.952             [0.976, 0.650, 0.488]   0.044   → 0   (intercept 0)
   Route A'  b+a/L       b=−0.149, a=2.331   [1.017, 0.629, 0.434]   0.015   → b=−0.149  (NEGATIVE, unphysical)
   Route B   b+a/L²      b=+0.266, a=3.024   [1.022, 0.602, 0.455]   0.035   → b=+0.266  (≈ Haar 0.269)
   [free]    4π²/L²      —                   [9.87, 4.39, 2.47]      —        → 0   (excluded: ~10× too big)
```
**Read honestly:** all rms are below the MC scatter, so **three points cannot decide**. But the pattern
is unambiguous: `b+a/L` (a *pure* `1/L` decay toward a floor) needs a **negative floor**, which is
unphysical, while `b+a/L²` lands on `b=+0.266`, within noise of the exact `β=0` Haar value `0.269`.
So the data are **consistent with saturation at the Haar scale** and **disfavour `1/L→0`**. The
box-gap is excluded outright.

### 4.4 What additional data decides

1. **`G″` at `β=2.5` for `L=5,6,7,8`** (MC, tight bars, `Landau`/block-Hessian): does it **flatten at
   `≈0.27`** (Route B) or **continue down through `0.27` toward `0`** (Route A)? This is the decisive
   curve.
2. **The susceptibility itself `χ_L = Var(a)` at `β=2.5`** (cleaner than the curvature: it is a first
   moment of the exact marginal / a direct MC average). Route A ⟺ `χ_L` grows with `L`; Route B ⟺
   `χ_L` plateaus. One plot settles it.
3. **Exact `β=0` reference curve** `L=5…16` (already computed here `[D1]`: `→0.2688`) to fix the
   expected plateau value.
4. **4D** constant-mode `χ_L` (the real target): a finite plateau of `Var(A₀)` in the 4D `L=8,12` model
   would be the `H1`-consistent signal; growth would be a no-go.

---

## 5 · Cleanest conditional statements (the hypothesis, named)

> **[C–B · saturation]** *If the constant-mode zero-momentum susceptibility is bounded,
> `(B-ineq)  χ_L = Σ_d(1−|d|/L)C(d) ≤ χ_max < ∞` uniformly in `L` (equivalently the constant-mode
> correlation length is finite, `|C(l,l')|≤C₀e^{−|l−l'|/ξ}`, `ξ<∞`), then `G″(a₀*) ≥ 1/χ_max > 0` for
> all `L` — the soft-mode curvature **saturates** at a positive L-independent floor, and no uniform
> gap is lost to decompactification.* **Label:** the implication `[T]` (Gaussian marginal (1.1));
> the hypothesis is **`[T]`** in the `1×L` chain (1D exponential mixing) and at `β=0` (i.i.d. Haar →
> `χ=3.7198`), and **`[C]`=`H1`** in 4D.

> **[C–A · decay]** *If instead `χ_L→∞` (zero-momentum susceptibility diverges / gapless constant
> mode), then `G″(a₀*) → 0`, with the exact rate `G″ ≍ L^{−p}` where `χ_L ≍ L^p` (`p=2` free box-gap
> `4π²/L²`; `p=1` critical `c/L`).* **Label:** the implication `[T]`; the hypothesis is **`[NG]`** for
> the `1×L` chain and is the *definition of no mass gap* in 4D. Realised only on the free field (acid
> test).

**The single hypothesis that decides it: `χ = Σ_d C(d) < ∞` — is the constant-mode zero-momentum
susceptibility finite?** This is the mass gap read in the constant-mode channel: a gap `m>0` forces
`χ<∞` (saturation) and `χ=∞` means no gap (decay). The note's real deliverable is the **exact
equivalence** `(decay vs saturation) ⟺ (χ=∞ vs χ<∞) ⟺ (no gap vs gap)`, plus the two hard facts that
move it: `β=0` gives `χ=3.7198` exactly `[T]`, and the `1×L` chain has `χ<∞` by 1D mixing `[T]`.

**Logical status, stated plainly.** Saturation is a **necessary** condition for the mass gap (since
`m>0 ⇒ χ<∞`), but `1/χ` is the zero-momentum susceptibility, **not** the pole mass — proving
saturation does **not** prove the gap, and proving `1/χ` is the physical mass is the further step. What
this note *does* settle: **decay would be a no-go** (it would kill the constant-mode route outright),
and in the models where the data live, decay does not happen.

---

## 6 · Acid tests `[T]+[V]`

1. **Free massless field.** `Φ≡0` (no Haar Jacobian), the constant mode is the free zero mode:
   `χ_L ≍ L²`, `G″ = β·4π²/L² → 0` (exponent 2). The tools return **0** — no manufactured floor
   `[D5]`: `4π²/L² = 9.87, 2.47, 0.62, 0.15, 0.039` for `L=2,4,8,16,32`, `→0`. **Passes.**
2. **Pure Haar `β=0`.** `G″ ∈ [0.205, 0.269]`, L-independent, `→1/χ=0.2688` (`[D1,D2]`). **Not** the
   free field (free is 0); the discriminator is `Φ`. **Passes discriminatingly.**
3. **The `(LC-W)` test is the acid test of the *sharp* route.** It fails pointwise (`[D4]`) — the tool
   honestly reports the Savvidy log-convexity rather than manufacturing log-concavity. **Passes (by
   failing where it should).**
4. **The `β=0` band `[A,1/χ]`.** `A=0.20511530` is reproduced exactly from the closed form
   (`G″=A−Bβ²` at `β=0`) — the exactness check on the convolution solver. `[T]`.

---

## 7 · Witness — `research/v3/w7/d6_decision.py`

`py -3.12 -I d6_decision.py` (numpy 2.4.6, scipy 1.18.1; exit 0). Exact `1×1` marginal by 6000-pt
Gauss–Legendre; `β=0` by FFT convolution of the one-link marginal; the `β=2.5` data by least squares.
Quoted output (abridged):
```
[D1] beta=0 L-fold convolution: G''(a₀*) L=1..16: 0.210,0.240,0.249,0.254,0.259,0.262,0.264,0.265
     -> 1/Var(rho1)=0.26883 (CLT), L-INDEPENDENT; kurtosis 2.57->2.97->3.   [T]+[V]
[D2] Var(a)=3.71982 (L-independent) => G'' -> 1/Var = 0.2688; CLT converges the marginal. [T]
[D3] beta=2.5 {1.01,0.65,0.42}: a/L: a=1.952 rms0.044 (->0); b+a/L: b=-0.149 rms0.015 (NEG floor);
     b+a/L^2: b=+0.266 rms0.035 (≈ Haar 0.269); box-gap 4pi^2/L^2 ~10x too big (excluded). [V]
[D4] (LC-W) exact 1x1: -(log w)'' min -0.37 (e3),-0.11 (e1): VIOLATED on ~half support (flanks);
     at beta-mode: G''^b-G''^0 = +0.431 (e3), +0.157 (e1): floor holds AT THE MODE. [V]
[D5] free box-gap 4pi^2/L^2 -> 0; Haar floor 0.269 L-independent. Tool clean on free field. [T]+[V]
```
No third-party sources downloaded; standard inputs (1D Gibbs exponential mixing / transfer matrix,
CLT, Prékopa, Brascamp–Lieb, Schur) cited from the corpus.

**Companion witnesses (same Wave D, same question, independent methods).** `research/v3/w7/d1_exact_L.py`
(exact transfer-matrix / FFT marginal for the `1×L` chain), `d2_decompact.py` + `d2_quick.py` (exact-target
Metropolis on `Var(a)=⟨a²⟩_c` directly — the criterion "`Var(a)` bounded ⟺ `G″→positive`; `∝L` ⟺ `G″→0`"
is identical to (1.1)), `d4_haar_wilson.py` (exact split `G″=G″_Haar+G″_Wilson` via conditional cumulants
of `S_W` — the same `w(a)=E_Haar[e^{−S_W}|a]` as §3.2). Both independent routes report the SAME Haar
floor `1/Var(ξ₃)=9/(4π²−6)=0.26883`, corroborating §3.1. The decisive measurement this note calls for
(§4.4 item 2) is exactly the `Var(a)`-vs-`L` plateau test of `d2_*.py`.

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| dichotomy `G″(a₀*)→1/χ`, `χ=Σ_d C(d)`; `χ<∞⇒saturation`, `χ=∞⇒decay`, rate `L^{−p}` | §1 (1.1) | **[T]** (Gaussian marginal) / **[C]** (CLT) |
| `β=0` Haar marginal → Gaussian, `χ=Var_Haar(ξ₃)=3.7198`, `G″→0.2688`, L-independent band `[0.205,0.269]` | §3.1, `[D1,D2]` | **[T]**+**[V]** |
| `1×L` chain: 1D exponential mixing ⇒ `χ<∞` ⇒ saturation; Route A `1/L→0` impossible | §4.1 | **[T]** (standard 1D) → **[NG]** for Route A on the chain |
| `(LC-W)` global log-concavity of `w` ⇒ `G″^β≥G″^0` | §3.2 | **[T]** (implication); **(LC-W) false pointwise** `[NG]` |
| at the β-mode `G″^β(a₀*) ≥ G″^0(a₀*)` (1×1 exact, Cartan+Savvidy, β=2.5) | §4.2, `[D4]` | **[V]** |
| `β=2.5` data consistent with `b+a/L²` to `b=+0.266≈0.269`; pure `1/L→0` needs negative intercept | §4.3, `[D3]` | **[V]** (3 points, undecided) |
| box-gap `4π²/L²` excluded (~10× data) | §4.3, `[D5]` | **[V]** |
| `(B-ineq)` finite `χ` ⇒ saturation (the needed inequality) | §5 | **[C]** in 4D (`=H1`), **[T]** in 1D/β=0 |
| **`χ<∞` in 4D at `β=2.5`** | — | **[O] = H1** |
| **uniform-in-`L` soft-mode gap `δ_W>0`; Clay** | — | **[O]** / **0%** |

---

## 9 · Sharpest line

**The decompactification question is exactly "is the constant-mode zero-momentum susceptibility
`χ=Σ_d C(d)` finite", i.e. the mass gap in the `k=0` channel: `G″(a₀*)→1/χ` saturates iff `χ<∞`,
decays as `L^{−p}` iff `χ_L≍L^p`.** At `β=0` the answer is a rigorous `χ=3.7198`, `G″→0.2688`,
L-independent `[T]+[V]`; in the `1×L` chain (where the `1/L` data live) 1D exponential mixing gives
`χ<∞`, so **Route A is a no-go and the `≈1/L` is a finite-size transient** `[T]/[NG]`. The `β=2.5`
data agree (`1/L²` fit → `b=+0.266≈0.269`; a pure `1/L→0` fit needs a negative, unphysical intercept;
the free box-gap is ~10× too big) `[V]`. The sharp Prékopa/log-concavity route is **false pointwise**
(the Savvidy flanks dip below Haar) — the floor is the *mode-value* / susceptibility bound `(B-ineq)`,
the corpus's measure-carried `D`-defect. So **Route B survives and Route A dies**; the residue for 4D is
the named hypothesis `χ<∞` = `H1` = the mass gap, now shown to be exactly what the decompactification
of the soft mode *is*. Saturation is **necessary** for the gap, not sufficient (`1/χ` is the
susceptibility, not the pole mass). **No QED. Clay 0%.**

---

*Refs.* `SM1-SOFT-EXACT.md` §0,§4,§6 (`G″(a₀*)` data, β=0 `0.23–0.25` L-independent, closed form `A`,
Fisher); `SM4-CONSTANT-MODE.md` §1–§3 (exact constant-mode marginal, `L`-dependence, twist, β-scan);
`SM3-ANISOTROPIC-RG.md` §3–§5 (fixed-point `[NG]` "no blocking lifts `G″_soft`"; good-set floor = Haar
`N/12`); `SM2-RANK-ONE.md` §3 (ceiling `m_soft≤G″`, BL-carried positivity); `H2-DEFECT-BOUND.md`
(`G″=S−D`, measure-carried defect, Fisher identity); `docs/paper-v3/main.tex` §6.7 (the deferred
question). External, standard, from the corpus: 1D Gibbs measure exponential mixing / transfer matrix;
central limit theorem; Prékopa (1971); Brascamp–Lieb (1976); Schur/Haynsworth. Witness
`research/v3/w7/d6_decision.py`. **No third-party sources downloaded.**
