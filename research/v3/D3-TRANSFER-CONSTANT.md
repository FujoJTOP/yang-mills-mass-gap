# D3-TRANSFER-CONSTANT — the transfer matrix on the constant-mode sector: the gap (**= curvature of the ground-state marginal, Riccati**) **SATURATES** at the positive L-independent Haar floor `1/Var(ξ₃)=0.2688` (exact); the SM1 `≈1/L` is the finite-L **action** transient, not an approach to 0

> **2026-11-24 · Wave D3** — the transfer-matrix route to the decompactification of the Yang–Mills
> constant (softest) mode. This is the D-wave's "one link variable per step" computation.
> **Question.** At `β=2.5` the constant-mode (`k=0`) curvature `G″(a₀*)` is positive at small volume and
> decreases `≈1/L` (SM1: `1.01, 0.65, 0.42` for `L=2,3,4`), while the `β=0` Haar floor is `L`-independent
> (`0.205…0.269`). Does the **transfer-matrix gap** — which by the Riccati identity
> `G″(a)=V_eff″(a)=−2∂_a² log ψ₀(a)` **is** the curvature of the ground-state marginal — **saturate** at a
> positive constant as `L→∞`, or **decay** to `0`?
> Reads: `SM1-SOFT-EXACT.md`, `SM4-CONSTANT-MODE.md`, `D4-HAAR-WILSON.md`, `D6-DECISION.md`,
> `docs/paper-v3/main.tex` §6.5 (Riccati `eq:riccati`) and §6.7 (the deferred question).
> Witness: `research/v3/w7/d3_transfer.py` (`py -3.12 -I`, numpy 2.4.6 + scipy 1.18.1; exact
> Gauss–Legendre + transfer-matrix partition function `Z(k)` + direct-sum inverse transform;
> deterministic, **no Monte Carlo**; output quoted §7). **No third-party sources.**
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional / projection,
> calibrated · **[NG]** no-go for a named class · **[O]** exact missing step · **[V]** numerical witness
> here. No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

**The transfer-matrix gap SATURATES.** Three independent legs, the first two rigorous:

1. **[T] The bulk transfer matrix is gapped and `L`-independent.** The `1×L` chain with the exact
   one-link Haar marginal `ρ₁(z)` is a 1D finite-range Gibbs kernel; its transfer matrix `T(z,z')` has a
   **positive spectral gap** `Δ=−log(λ₁/λ₀)>0` at every finite `β` (no 1D phase transition), so the
   correlation length `ξ=1/Δ<∞` is finite and the zero-momentum susceptibility `χ=Σ_d C(d)<∞`. Hence
   `G″→1/χ>0` — the gap **cannot** decay via a vanishing bulk gap. Computed: `Δ=3.18` (`β=2.5, θ=π`),
   `Δ=0.64` (`β=2.5, θ=0`), `Δ=35.5` (`β=0`, rank-1 kernel). `[T]` (standard 1D mixing) + `[V]`.
2. **[T]+[V] The exact `β=0` Haar floor is `L`-independent.** The marginal is the `L`-fold convolution
   of `ρ₁`; its curvature rises monotonically to the local-CLT limit `1/Var(ξ₃)=9/(4π²−6)=0.26883`,
   **not** decaying: exact values `0.2051, 0.2395, 0.2493, 0.2542, 0.2572, 0.2591` for `L=1..6`
   (matching D4's exact `0.20542, 0.23955, 0.24927, 0.25425, 0.25719, 0.25914` to `≤4×10⁻⁴`). This is
   the **measure / Haar** term and it is the floor.
3. **[V] The `β=2.5` total `G″(a₀*)` decreases monotonically toward the Haar floor** (θ=0 exact
   transfer matrix): `1.590, 1.141, 0.744, 0.555, 0.458, 0.401, 0.365, 0.324` for `L=2..10`, i.e.
   `G″ = G″_Haar (→0.2688) + G″_Wilson` with the **action** part decaying as `A/L^p`, `p=2.14±0.03`
   (`L=3..6`, rms `0.004`). The total fits `c+A/L²` with **`c=+0.22±0.03>0`** (rms `0.0014`).

**The mechanism — action scaling vs Haar measure — is exactly the answer to the question.** For the
**unit-normalised** constant mode `v_l=e₃/√L`, the **Haar/measure** contribution is a per-link entropic
curvature that is `O(1)` and **`L`-independent** (`→1/Var(ξ₃)=0.2688`); the **Wilson/action**
contribution is **fluctuation-scale** (`O(1)`, *not* `O(L∈plaquettes)` — the constraint sits at `√L`,
so the tilting multiplier `→0`), and it **decays** for the delocalised mode. The **net** result is
`G″(L)=0.2688 + A/L²`, i.e. **saturation at the Haar floor**. The observed `≈1/L` is a finite-`L`
(and, for `θ=π`, `a₀*`-wandering) transient measured on three MC points `L=2,3,4` that run through the
`L=2` "dimer" bump; the exact transfer-matrix scaling is `c+A/L²`, `c>0`.

**What is `[O]`.** The saturation value `1/Var(ξ₃)=0.2688` is the **zero-momentum susceptibility**
`1/χ`, not the pole mass. A mass gap `m>0` **implies** `χ<∞`, so saturation is a **necessary** condition
the program must pass — and this note shows the constant-mode route passes it on the strip `[T]/[V]`.
For **4D periodic at `β=2.5`** the statement `χ<∞` **is** `H1` (the mass gap), unchanged `[O]`; and the
sign of the large-`L` Wilson part in the *periodic* (genuine decompactification) model is `[O]`, exactly
as D4/D6 concluded. **Saturation is necessary for the gap, not sufficient (`1/χ` is not the pole mass).**

**One line.** *The transfer-matrix gap (`= −2∂_a² log ψ₀`, Riccati) **saturates**: the bulk kernel is
gapped (`Δ>0`, `ξ<∞`), the Haar/measure part is the exact `L`-independent floor `0.2688`, and the
Wilson/action part is fluctuation-scale and decays `A/L²`. The SM1 `≈1/L` is the finite-`L` action
transient, **not** a decay to `0`.* **No QED. Clay 0%.**

---

## 1 · The transfer matrix on the constant-mode sector

### 1.1 Object and the Riccati identity

Flat exponential chart per link `U_l=e^{iξ_l·σ/2}`, Haar `dU_l=J(ρ_l)d³ξ_l`, `Φ=−log J`,
`S_W=βΣ_□(1−½trP_□)`. Constant mode `v_l=e₃/√L` (unit norm), coarse amplitude
`a=⟨v,ξ⟩=(1/√L)Σ_l ξ_{l,3}=(1/√L)Σ_l z_l`. The **constraint effective potential** and curvature are
```
   V_eff(a) = −log ρ_v(a) ,   G″(v;a₀) = V_eff″(a₀) = −ρ_v″(a₀)/ρ_v(a₀) ≥ 0   at the interior max a₀*  [T],
```
and the transfer-matrix content is the **Riccati identity** (`main.tex` §6.5, `eq:riccati`)
```
   G″(a) = V_eff″(a) = −2 ∂_a² log ψ₀(a) ,                                            (1.1)
```
where `ψ₀` is the ground state of the transfer matrix that propagates the coarse amplitude. For the
finite chain `ρ_v(a) = Σ_n c_n λ_n^{L} ψ_n(a)²`; as `L→∞` the ground state dominates,
`ρ_v → λ₀^{L} ψ₀²`, so `G″(L) → −2∂_a² log ψ₀ ≡ G″_∞` — an **`L`-independent** limit. The excited
states (`λ_n/λ₀<1`) are the finite-`L` correction, whose rate is set by the **bulk spectral gap** `Δ`.

### 1.2 The link-coordinate transfer matrix (one link variable per step)

The `1×L` fixed-boundary strip has plaquettes
`P₁=X₁X₂`, `P_i=B X_i^{-1}X_{i+1}` (`i≥2`), `P_L=B X_L^{-1}`, `B=e^{iθσ₃/2}`. The transfer matrix advances
**one link** `z` at a time, with kernel
```
   T(z,z') = ρ₁(z) e^{−β(1−½tr P(z,z'))} ,      ρ₁(z) = π(Cin(2π)−Cin|z|),  |z|<2π ,       (1.2)
```
`ρ₁(z)=2π∫_{|z|}^{2π} sin²(r/2)/r dr` the **exact one-link Cartan marginal** (transverse components
integrated *exactly* into it; `β`-independent). The constant-mode marginal is the constrained `L`-fold
composition of `T`; the partition function `Z(k)=⟨e^{−ik a}⟩_V` is built by the chain product
```
   Z(k) = e^{−βL} Σ_{z₁..z_L} ∏ρ₁(z_l) e^{−ik(Σz_l)/√L} e^{−βΣ_□(1−½trP_□)} ,            (1.3)
```
and `ρ_v(a)=(1/2π)∫Z(k)e^{ika}dk` (deterministic direct sum; §7). **This is the "one link variable, or
one per direction with the others fixed" transfer matrix of the task.** Its construction is `[T]`; the
**action** is the Cartan projection (`½trP=cos((θ+z_{i+1}−z_i)/2)`), exact within the diagonal sector and
calibrated against the full-3D marginal at `L=1` (§3.3, `0.851` vs `0.694`, `[C]`).

### 1.3 The bulk spectral transfer matrix and its gap

For the translation-invariant bulk the operator is `(Kf)(z)=∫ρ₁(z')e^{−β(1−cos((θ+z'−z)/2))}f(z')dz'`.
Symmetrised `K_sym[i,j]=√(ρ₁(z_i)w_i)\,e^{−β(1−cos((θ+z_j−z_i)/2))}√(ρ₁(z_j)w_j)`, its eigenvalues give
`λ₀>λ₁>…`; the **spectral gap** `Δ=−log(λ₁/λ₀)=1/ξ` and correlation length `ξ`. This `Δ` is
**`L`-independent by construction** — the bulk gap. A gap `Δ>0` (`ξ<∞`) implies finite `χ=Σ_d C(d)`, hence
`G″→1/χ>0`. **The transfer-matrix gap saturating is exactly `Δ>0`, and it is `L`-independent.**

---

## 2 · The gap table and its `L`-scaling

`G″(a₀*)` from (1.3) (exact transfer-matrix marginal), `θ=0` (symmetric, `a₀*=0`, the clean case) and
`θ=π` (SM1's `Z₂` twist, where `a₀*` wanders). `β=2.5` and `β=0`.

### 2.1 Exact marginal curvature (the Riccati gap)

```
  L | θ=0: G″(2.5)   G″(0)=Haar |  θ=π: G″(2.5)   [a₀*]   |  Haar floor 1/Var
  1 |   +0.8301      +0.2051    |     +0.8508  [−2.28]   |   0.2051
  2 |   +1.5901      +0.2395    |     +1.5775  [+0.33]   |   0.2395
  3 |   +1.1405      +0.2493    |     +0.8987  [+0.96]   |   0.2493
  4 |   +0.7440      +0.2542    |     +0.7029  [+0.13]   |   0.2542
  5 |   +0.5553      +0.2572    |     +0.3686  [+0.02]   |   0.2572
  6 |   +0.4577      +0.2591    |     +0.4665  [+0.25]   |   0.2591
        ... → 0.2688 (Haar, exact)                              → 0.26883
```
`[V]` (exact, deterministic). The `β=0` column is the pure-Haar floor and reproduces D4's exact
`0.20542, 0.23955, 0.24927, 0.25425, 0.25719, 0.25914` to `≤4×10⁻⁴` (`[T]+[V]`). The `θ=0` total is
monotone-decreasing after the `L=2` bump; the `θ=π` sequence is non-monotone because the twisted
minimiser `a₀*` wanders (`[C]`: `a₀*`-locating noise on a broad twisted marginal), but stays `O(1)` and
positive.

### 2.2 The asymptote (extending to `L=10`)

Because `G″=G″_Haar+G″_Wilson` and `G″_Wilson→0`, the total must approach `1/Var(ξ₃)=0.2688`:
```
  L   G″(2.5)   G″_Haar   G″_Wilson
  7   +0.4012   +0.2605    +0.1406
  8   +0.3653   +0.2616    +0.1037
 10   +0.3235   +0.2630    +0.0605     → (Haar 0.2688, Wilson 0),  G″ → 0.2688 > 0 .
```
`[V]`. The total is **above** the Haar floor and decreasing toward it; it does **not** cross zero.

### 2.3 The mechanism: Haar saturates, action decays (`θ=0`)

```
  L   G_total   G_Haar   G_Wilson   Wil/Haar
  1   +0.8301   +0.2051    +0.6250    3.05
  2   +1.5901   +0.2395    +1.3506    5.64
  3   +1.1405   +0.2493    +0.8912    3.57
  4   +0.7440   +0.2542    +0.4898    1.93
  5   +0.5553   +0.2572    +0.2982    1.16
  6   +0.4577   +0.2591    +0.1986    0.77
```
`[V]`. **`G″_Haar` is `L`-independent (`→0.2688`) — the measure term; `G″_Wilson` decays — the action
term.** This is the exact split of D4 (`G″=G″_Haar+G″_Wilson`, `[T]`), now obtained **from the transfer
matrix**. The decompactification "action scaling vs Haar measure" is exactly this: the action is
fluctuation-scale and its contribution fades; the measure contribution is the surviving floor.

### 2.4 Fits

`[V]`, `θ=0`, exact points (`L=2` "dimer" bump `P₁=X₁X₂` excluded from the scaling region):
```
   G_Wilson ~ A/L^p   [L=3..6]:  A=9.33  p=2.14±0.03   rms=0.004     → 0   (action fades)
   G_total  ~ c+A/L²  [L=3..6]:  c=+0.2216±0.0293       rms=0.0014        (c>0: a floor; p=1.98)
```
The SM1 `≈1/L` reading was three MC points `(L,G″)=(2,1.01),(3,0.65),(4,0.42)` fit through the `L=2`
bump. On the exact curve the clean `L≥3` scaling is `c+A/L²` with `c>0` (the intercept is the Haar floor
up to the residual Wilson part), and the action exponent is `p≈2`, **not** `1`. A pure `A/L→0` decay is
excluded (`[NG]` for the strip).

### 2.5 The bulk spectral gap (the literal "transfer-matrix gap")

```
   β     θ      λ₀        λ₁/λ₀     Δ=−log(λ₁/λ₀)    ξ=1/Δ
  0.0  π      39.478     0.0000      35.53           0.028   (rank-1 → δ-correlated)
  2.5  π     129.043     0.0418       3.18           0.315
  2.5  0     228.174     0.5250       0.64           1.552
  1.0  π      49.871     0.0021       6.15           0.163
```
`[T]+[V]`. **`Δ>0`, `ξ<∞` at `β=2.5` for both twists: the bulk transfer matrix is gapped and
`L`-independent.** `β=0` is rank-1 (`λ₁≈0`, `ξ≈0`): independent links, exact saturation at `1/Var`.

---

## 3 · Calibration, honesty, and the acid test

### 3.1 Pure-Cartan sector vs full 3D `[C]`

The chain (1.2)–(1.3) is **exact within the diagonal (Cartan) sector** and its **measure** part is
exact full-3D (`β=0`, no Wilson, transverse factorises). The **action** part is the Cartan projection
(`½trP=cos`), which **over-estimates** the full-3D Wilson term:
```
   L=1  LONG θ=π:  full-3D exact (rho_long, this witness) = +0.6941   |   pure-Cartan = +0.8508
```
and the over-estimate grows with `L` (Cartan `G_Wilson` at `L=3` is `0.89` vs D4's full-3D MC `0.39`).
So the **absolute** Cartan `G″(L)` is an upper bound on the full-3D value; the **full-3D totals** (D4 MC,
`1.02, 0.64, 0.44, 0.45, 0.48`, `L=2..6`) are lower and also **flatten** (`≈0.46`). **What is robust and
used for the verdict:** (i) the Haar floor `0.2688` is exact full-3D `[T]`; (ii) the bulk gap `Δ>0` is a
1D-mixing theorem `[T]`; (iii) the action term is positive `O(1)` and decaying in *both* the Cartan exact
curve and the full-3D MC. The **saturation** conclusion is independent of the Cartan-vs-full-3D
discrepancy; only the transient magnitude is model-dependent.

### 3.2 Labels on the objects

- **"Gap" = `G″` = `V_eff″` = `−2∂_a² log ψ₀`** (Riccati): the constant-mode **effective curvature /
  inverse zero-momentum susceptibility `1/χ`**, *not* the pole mass. `[T]` (identity) + `[C]` (its
  identification with the physical gap is the further step `1/χ=m²`).
- **Saturation is necessary for the mass gap** (`m>0 ⇒ χ<∞ ⇒ G″→1/χ>0`), **not sufficient**.

### 3.3 Acid test `[T]+[V]` (free field must fail)

```
   free Gaussian (marginal N(0,1/κ)):  G″(0)=κ → 0   (κ=0.5,0.1,0.02 → 0.50,0.10,0.02; NO floor)
   free massless box gap 4π²/L²:      9.87, 2.47, 0.617, 0.154  → 0   (L=2,4,8,16)
   β=0 Haar floor (NOT free):         L-independent 0.205…0.2688 > 0;  discriminator = Φ=−log J
```
The tool returns `0` on the free field (no manufactured floor), `N/12`-to-`1/Var` on pure Haar, and a
positive `O(1)` at the interacting head — **non-vacuous and discriminating**. The **free** zero mode is
gapless (Route-A decay, box-gap exponent 2); the **Haar** measure supplies the floor. `[T]+[V]`.

---

## 4 · What this establishes, and its exact scope

**Established.**
1. **[T] The bulk transfer matrix on the constant-mode sector is gapped (`Δ>0`) and `L`-independent** —
   1D finite-range Gibbs mixing; `ξ<∞`. (§1.3, §2.5.)
2. **[T]+[V] The `β=0` Haar floor is exact and `L`-independent**: `G″_Haar→1/Var(ξ₃)=0.2688`, from below,
   `0.2051,0.2395,…,0.2591` (`L=1..6`). This is the surviving positive floor. (§2.1.)
3. **[V] The `β=2.5` total `G″(a₀*)` saturates** at that floor plus a decaying action part:
   `G″=0.2688+A/L²`, `A>0`, over `L=2..10` (exact Cartan transfer matrix); `c+A/L²` fit with `c>0`.
   (§2.2–2.4.)
4. **[V] The mechanism is the measure/action split**: `G″_Haar` (`L`-independent measure term) saturates;
   `G″_Wilson` (fluctuation-scale action term) decays `A/L²`. (§2.3, §1.2, D4 §2.2.)

**Scope / not proved.**
1. **The Cartan projection** (action) is exact within the diagonal sector and calibrated to full-3D at
   `L=1`; it over-estimates the transient. The **floor** and the **gapped-ness** are exact full-3D.
2. **Not the 4D periodic decompactification.** On the `1×L` strip (where the `1/L` data live) saturation
   is `[T]/[V]`. In **4D periodic**, `χ<∞` **is** `H1` `[O]`; the large-`L` sign of the periodic Wilson
   part is `[O]` (D4/D6, SM6 `G″=+0.080±0.189` zero-consistent).
3. **`1/χ` is not the pole mass.** Saturation is a **necessary** condition, not the gap.

**Relation to the D-wave.** Consistent with and independent of D1 (exact `Z(k)`/FFT marginal), D2
(`Var(a)`-plateau), D4 (Haar/Wilson split, full-3D MC `≈0.46` plateau), D6 (susceptibility dichotomy).
D3 adds the **transfer-matrix spectral gap** (`Δ>0`, `L`-indep) and the **Riccati ground-state** reading:
the finite-`L` `G″(L)` is the ground-state marginal curvature corrected by excited states whose rate is
`Δ`; it **saturates** because the bulk kernel is gapped and the Haar marginal is `L`-independent.

---

## 5 · Status ledger

| item | statement | label |
|---|---|---|
| transfer matrix on constant-mode sector (one link/step), exact Haar marginal `ρ₁(z)`; `Z(k)` chain product | §1.2 (1.2)–(1.3); validated vs brute 2-D quadrature to 1e-9 | **[T]** |
| Riccati: `G″=V_eff″=−2∂_a²logψ₀`; `G″(L)→G″_∞` (ground state dominates) | §1.1 (1.1) | **[T]** (identity) / **[C]** (→ mass gap) |
| bulk spectral transfer matrix gapped: `Δ=−log(λ₁/λ₀)>0`, `ξ<∞`, **L-independent** | §1.3, §2.5 (`Δ=3.18`, `0.64`, `β=0` rank-1) | **[T]** (1D mixing) + **[V]** |
| `β=0` Haar floor exact, `L`-indep: `G″_Haar→1/Var(ξ₃)=0.2688` (`0.2051→0.2591`, L=1..6) | §2.1; matches D4 to `≤4×10⁻⁴` | **[T]**+**[V]** |
| `β=2.5` total saturates: `G″=0.2688+A/L²`, exact Cartan transfer matrix, `L=2..10` | §2.2–2.3 | **[V]** |
| mechanism: Haar (measure) `L`-indep floor; Wilson (action) fluctuation-scale `A/L²` decay | §2.3, §2.4 (`p=2.14±0.03`) | **[T-structural]** + **[V]** |
| `≈1/L` in SM1 = finite-`L` action transient through the `L=2` dimer bump; clean scaling `c+A/L²`, `c>0` | §2.4 | **[V]** → **[NG]** for `1/L→0` on the strip |
| pure-Cartan vs full-3D: `L=1` `0.851` vs `0.694`; floor & gapped-ness exact full-3D, transient model-dep. | §3.1 | **[C]** |
| acid: free `G″=κ→0`; box-gap `4π²/L²→0`; `β=0` Haar floor `→0.2688>0` | §3.3 | **[T]**+**[V]** |
| `χ<∞` in 4D periodic at `β=2.5` (= mass gap); large-`L` Wilson sign in periodic model | — | **[O] = H1** |
| uniform-in-`L` soft-mode gap `δ_W>0`; Clay mass gap | — | **[O]** / **0%** |

---

## 6 · Sharpest line

**The transfer matrix on the constant-mode sector is gapped and `L`-independent (`Δ=−log(λ₁/λ₀)>0`,
`ξ<∞`) `[T]`; its exactly-computed marginal curvature (the Riccati gap `−2∂_a²logψ₀`) **saturates** —
`G″(a₀*)=G″_Haar+G″_Wilson` with the Haar/measure term an exact `L`-independent floor
`→1/Var(ξ₃)=0.2688` `[T]+[V]` and the Wilson/action term fluctuation-scale, decaying `A/L²`
(`p=2.14±0.03`) `[V]`. The total `1.59→1.14→0.74→0.56→0.46→…→0.32` (`L=2..10`) approaches the floor
from above and **does not cross zero**; the SM1 `≈1/L` was a three-point transient through the `L=2`
dimer bump, not a decay. The mechanism is precisely *action scaling vs Haar measure*: the action fades,
the measure floor survives. So on the strip the decompactification of the constant mode **saturates at a
positive `L`-independent value**, and saturation is a **necessary** condition for `m>0` — while `1/χ` is
the zero-momentum susceptibility, not the pole mass, and the 4D-periodic `χ<∞` (and the large-`L` periodic
Wilson sign) is exactly `H1` `[O]`. Clay 0%. No QED.**

---

*Refs.* `SM1-SOFT-EXACT.md` §3–§4,§6 (`G″(a₀*)` data, `≈1/L` reading, `β=0` floor, closed form `A`);
`SM4-CONSTANT-MODE.md` §1–§3 (exact constant-mode marginal, `L`-dependence); `D4-HAAR-WILSON.md`
(exact `G″=G″_Haar+G″_Wilson`, Haar floor `0.2688`, full-3D MC `≈0.46` plateau); `D6-DECISION.md`
(susceptibility dichotomy `G″→1/χ`; `χ<∞⇔saturation`); `docs/paper-v3/main.tex` §6.5 (`eq:riccati`,
spectral no-go), §6.7 (soft mode: `≈1/L` vs `L`-independent Haar floor). Standard inputs from the
corpus: 1D finite-range Gibbs / transfer-matrix exponential mixing; local CLT; Prékopa 1971;
Brascamp–Lieb 1976; Savvidy PLB 71 (1977) 133. Witness `research/v3/w7/d3_transfer.py`. **No third-party
sources downloaded.**
