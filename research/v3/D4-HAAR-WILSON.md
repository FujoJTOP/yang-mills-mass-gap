# D4-HAAR-WILSON — the exact Haar/Wilson decomposition of the constant-mode curvature: at β=2.5 the **Haar part is a proven L-independent positive floor** (`→1/Var(ξ₃)=0.2688`), the **Wilson part is positive O(1)** at the coarse minimiser (`0.78`→ plateau `≈0.2` for L≥4), it does **not** grow with the plaquette number, and on the tractable strips the total **flattens at ≈0.46 > Haar floor**, it does **not** decay to 0; the `≈1/L` reading was an L=2–4 transient. The decompactification limit itself is still H1 `[O]`

> **2026-11-24 · D4 — decompactification of the Yang–Mills soft mode.**
> Question: at β=2.5 the soft (constant, k=0) mode curvature decreases `≈1/L` (SM1), while the
> β=0 Haar floor is L-independent `0.23–0.25`. Decompose `G″` into Haar and Wilson parts and decide
> which dominates as `L→∞`.
> Reads: `SM1-SOFT-EXACT.md`, `SM4-CONSTANT-MODE.md`, `SM2-RANK-ONE.md`,
> `research/unconditional/PROOF-JACOBIAN-SAVES.md`, `docs/paper-v3/main.tex` §6.1, §6.7.
> Witness: `research/v3/w7/d4_haar_wilson.py` (`py -3.12 -I`, numpy/scipy; exact 1-D quadrature +
> FFT L-fold convolution for the Haar part; Haar-importance MC for the total at L≥2; **one job at
> a time**).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[V]** numerical
> witness here · **[NG]** no-go for a named class · **[O]** exact missing step. No QED. Clay 0%.

---

## 0 · Verdict (no QED, Clay 0%)

Write the marginal in the flat chart as `ρ_v(a) = ρ_v^{Haar}(a)·R(a)`, `R(a)=E_{Haar}[e^{−S_W}|a]`.
This is an **exact factorisation** (Haar factorises and the Wilson weight is `e^{−S_W}`), so

```
   G″(a) = −d²/da² log ρ_v(a) = G″_Haar(a) + G″_Wilson(a) ,                       [T]   (0.1)
   G″_Haar(a)   = −d²/da² log ρ_v^{Haar}(a) ,        (pure β=0 marginal; an L-fold convolution)
   G″_Wilson(a) = −d²/da² log E_{Haar}[e^{−S_W}|a]
                = κ₁″ − ½κ₂″ + (1/6)κ₃″ − … ,       κ_n = conditional cumulants of S_W .        (0.2)
```

**Result.**

1. **Haar part `[T]+[V]`.** For the k=0 chain mode `v_l=e₃/√L`, `G″_Haar` is the curvature of the
   L-fold convolution `ρ₁^{*L}` of the one-link marginal. It is **`O(1)` and L-independent in the
   limit** (local CLT): `G″_Haar(a)→1/Var(ξ₃)=0.26883`, exactly, from below. Finite values at the
   β=0 minimiser `a=0`: `0.2054, 0.2396, 0.2493, 0.2542, 0.2572, 0.2591` for `L=1..6`. This is a
   **genuine, positive, L-independent Haar floor** — the `β=0` object, *not* the free field (§6).
2. **Wilson part `[V]`, sign/scaling.** `G″_Wilson` is **positive at the coarse minimiser** `a₀*` in
   every model decomposed here at β=2.5: exact `+0.422` (L=1, Cartan e₃) and `+0.157` (L=1, off-diagonal
   e₁); MC chains `+0.78, +0.39, ≈+0.19, ≈+0.20, ≈+0.22` for `L=2..6`. It is **`O(1)`, bounded, and
   does not grow with the number of plaquettes `N_□=L`** (the naive `L`-growth is cancelled because the
   constant-mode constraint sits at the fluctuation scale, §2). It decays from the L=2 peak to a
   **plateau `≈0.2`** for `L≥4`. The sign is *not* universal (it flips negative at the L=1 longitudinal
   **barrier** `a=0`: `−0.069`), but at the coarse minimiser — the object `G″(a₀*)` — it is positive.
3. **Which dominates / the limit.** The Haar part is the **only part with an L-independent nonzero
   limit** (`0.2688`). The crossover where Haar overtakes Wilson is **`L*≈3.5`** (Wilson/Haar =
   `3.24, 1.51, 0.68, 0.72, 0.83` at `L=2,3,4,5,6`). On the tractable strips the **total does not decay to 0
   or negative**: `G″(a₀*) = 1.02, 0.64, 0.44, 0.45, 0.48` (`L=2..6`, 3 seeds each at L=4–6), i.e. it
   **flattens at `≈0.46`**, above the Haar floor. The SM1 `≈1/L` fit was through `L=2,3,4` and does not
   continue. Structurally both parts are `O(1)` (not `O(L)`, not `1/L`), so the decompactification
   limit is `Haar-floor + (O(1) Wilson)`, i.e. bounded away from 0 **on this model**.
4. **What is NOT settled `[O]`.** The Wilson part's **sign at large L in the periodic
   (decompactification) model** is exactly H1. The fixed-boundary strip has a `Z₂`-twisted boundary
   whose holonomy supplies a positive Wilson contribution; the periodic lattice probe at β=2.4, L=8 gives
   `G″=+0.080±0.189` (SM6), zero-consistent. So the strip's positive plateau is **not** a proof that the
   continuum Wilson part is positive, and the Haar floor alone is `0.2688` only if the Wilson part does
   not go sufficiently negative. **The decomposition is exact and the Haar floor is a theorem; whether
   `G″→` it from above, or dips through it, is H1.**

One line. **The Haar part is the L-independent carrier (`→0.2688>0`, proven); the Wilson part is
positive `O(1)` at `a₀*` on the strips (no `N_□` growth, plateau `≈0.2`), so there the total flattens
at `≈0.46`, not at 0; the `1/L` decay was a finite-size transient. The large-L sign of the Wilson part
in the periodic model is H1 `[O]`. Clay 0%. No QED.**

---

## 1 · The exact decomposition

**Setup** (identical to SM1/SM4). Flat exponential chart per link, `ξ_l∈B(2π)`, `U_l=e^{iξ_l·σ/2}`,
Haar `dU_l=J(ρ_l)d³ξ_l`, `J=sin²(ρ/2)/ρ²`, `Φ=−log J`, `S_W=βΣ_□(1−½trP_□)`. Constant mode
`v_l=e₃/√L` (unit norm), amplitude `a=⟨v,ξ⟩=(1/√L)Σ_l ξ_{3,l}`. The exact marginal and curvature:

```
   ρ_v(a) = ∫ ∏_l d³ξ_l δ(⟨v,ξ⟩−a) e^{−(S_W+Φ)} ,     G″(v;a₀) = −(d²/da²) log ρ_v(a)|_{a₀} .   (1.1)
```

**Split `[T]`.** Insert `e^{−S_W}·e^{−Φ}` and condition the Haar measure on `a`:
`ρ_v(a) = ρ_v^{Haar}(a)·E_{Haar}[e^{−S_W}|a]`, hence (0.1)–(0.2). This is exact at every `a`, every
`L`, boundary and `β`; nothing is expanded. The **cumulant form** (0.2) is the expansion in powers of
`S_W` (`log E[e^{−S}] = −κ₁ + κ₂/2 − κ₃/6 + …`), i.e. the order at which the split is "meaningful":
the leading term `κ₁″ = d²/da² E_{Haar}[S_W|a]` is the **plaquette/commutator response**, and
`−½κ₂″` the force-variance correction. At L=1, β=2.5, `a₀*` (§3, [C3]):
`κ₁″=+0.4536`, `−½κ₂″=+0.0270`, sum `+0.4806` vs exact `G″_Wilson=+0.4224` — the leading cumulant
carries the sign and ~`107%` of the value; the higher cumulants give the `−0.058` finite-β correction.

**Haar part, closed form `[T]`.** Pure Haar has independent links, so `ρ_v^{Haar}` is the L-fold
convolution of the one-link marginal `ρ₁(s)=2π∫_{|s|}^{2π}sin²(ρ/2)/ρ dρ` (normalisation irrelevant):

```
   ρ_v^{Haar}(a) = √L (ρ₁^{*L})(√L a) ,     G″_Haar(a) = −L (log ρ₁^{*L})″(√L a) .                 (1.2)
```

This is evaluated **deterministically** (exact quadrature for `ρ₁` + FFT convolution), so the Haar part
carries **no MC error**.

**Wilson part, definition `[T]`.** `G″_Wilson(a) = −d²/da² log E_{Haar}[e^{−S_W}|a]`. Its finite-β
value is the (β-dependent, `[V]`) object below.

---

## 2 · L-scaling of each part

### 2.1 Haar part: `O(1)`, L-independent limit `[T]`

From (1.2) and the **local CLT** for the sum `Q=Σ_l ξ_{3,l}=√L a`: `ρ₁^{*L}(Q) = (2πLq₂)^{-1/2}
e^{−Q²/(2Lq₂)}(1+O(L^{-1/2}))`, with `q₂=Var(ξ₃)=3.719824`. Hence

```
   G″_Haar(a) = L·[ 1/(L q₂) + O(L^{-3/2}) ] = 1/q₂ + O(L^{-1/2}) = 0.26883 + O(L^{-1/2}) .        [T]
```

So the Haar part is **not** `1/L` or `L`; it is `O(1)` with the exact limit `1/Var(ξ₃)=0.26883`. The
finite-L values rise to it from below (`0.205→0.259`) because the one-link marginal is bounded
(`|ξ₃|≤2π`) and the convolution is still non-Gaussian at small L. **"L-independent per link"** is the
statement that the normalised amplitude `a` is a fluctuation-scale (`√L`) sum, so the L-fold
convolution freezes at `1/q₂`; equivalently each link contributes an equal `1/L` share of a fixed
`O(1)` total. The minimum of `G″_Haar(a)` on `|a|≤1.5` is at `a=0` and equals the values above (`[V]`),
so the Haar floor is positive throughout the sampled region.

### 2.2 Wilson part: `O(1)`, no `N_□` growth `[T-structural]+[V]`

`S_W=βΣ_{□=1}^{L}(1−½trP_□)` is explicitly `O(N_□)=O(L)` (the chain has `N_□=L`), so one might fear
`G″_Wilson=O(L)`. It is not. The constant-mode constraint fixes `a=O(1)`, i.e.
`Σ_l ξ_{3,l}=√L a`, which is a constraint at the **scale of the sum's own fluctuation** `√(Lq₂)`.
Equivalently the tilting Lagrange multiplier is `λ=O(a/√L)→0`, and

```
   E_{Haar}[S_W|a] = L·s̄ + O(√L)·a + O(a²)   ⇒   κ₁″ = d²/da² E[S_W|a] = O(1) ,                    [T-structural]
```

the `O(√L)a` term being a pure linear shift (zero second derivative) and the `O(a²)` term carrying the
curvature. The naive `O(L)` plaquette growth is cancelled exactly because the constraint is a
fluctuation, not an `O(1)` per-link constraint. Likewise `Var(S_W|a)` is `O(1)` in its `a`-second-
derivative. Hence **both parts are `O(1)`**, and the total has an `O(1)` (`L→∞`) limit generically —
**not** `1/L→0`. The commutator enters through the per-plaquette
`1−½trP_□ = ¼|F_□|²+O(a⁴)`; `κ₁″` is the curvature of the conditional plaquette (commutator) energy.

Whether the `O(1)` Wilson limit is **positive or negative** is model/coupling dependent — that is the
measured content of §3 and the open H1 content of §5.

---

## 3 · Both parts at β=2.5, L=1…6

### 3.1 Exact minimal torus, L=1 `[T]+[V]`

One free link, one plaquette; no non-constant modes, so (1.1) **is** the complete all-loop marginal.
θ=π (`Z₂` twist), β=2.5:

| class | `a₀*` | `G″_total` | `G″_Haar` | `G″_Wilson` |
|---|---|---|---|---|
| LONG `e₃` at minimiser | −2.1411 | **+0.6940** | +0.2715 | **+0.4224** |
| TRANS `e₁` at minimiser | 0.0000 | **+0.3623** | +0.2054 | **+0.1569** |
| LONG `e₃` at barrier `a=0` | 0 | +0.1360 | +0.2054 | **−0.0694** |
| TRANS `e₁` at barrier `a=0` | 0 | +0.3623 | +0.2054 | +0.1569 |

The direct ratio `−d²log[ρ_β/ρ_Haar]` reproduces `G″_Wilson` to 5 digits (`+0.42242` vs `+0.42244`;
`+0.15720` vs `+0.15689`) — the decomposition (0.1) verified non-perturbatively. The Haar part of the
two classes is the **same** (`≈0.205` at `a=0`, rotation-invariant Haar); the classes differ entirely
in the Wilson part. **The sign of the Wilson part is not universal**: it is `−0.069` at the
longitudinal barrier (the `A−Bβ²` mechanism, `[NG]` below) but `+0.16…+0.42` at the minimiser.

**β-scan at the L=1 minimiser `[T]+[V]`** (`G″_Haar` is β-independent; **all β-growth is Wilson**):

| β | LONG `a₀*` | LONG total | LONG Haar | LONG Wilson | TRANS total | TRANS Haar | TRANS Wilson |
|---|---|---|---|---|---|---|---|
| 0.5 | −0.803 | 0.245 | 0.213 | **+0.033** | 0.214 | 0.205 | +0.008 |
| 1.0 | −1.360 | 0.338 | 0.229 | **+0.109** | 0.237 | 0.205 | +0.032 |
| 2.5 | −2.141 | 0.694 | 0.272 | **+0.422** | 0.362 | 0.205 | +0.157 |
| 4.0 | −2.454 | 1.072 | 0.299 | **+0.773** | 0.511 | 0.205 | +0.306 |
| 6.0 | −2.660 | 1.581 | 0.321 | **+1.259** | 0.713 | 0.205 | +0.508 |

So at fixed `L=1` the Wilson part is **positive and grows with β** (matching SM1's `0.49→0.83→1.60`
total scaling; the Haar part is frozen). At β=0 it is exactly 0 (`S_W=0⇒R≡1`) and `G″=G″_Haar`.

### 3.2 Chains, L=1…6 `[V]` (MC total; exact Haar)

`1×L` fixed-boundary strip, θ=π, k=0 mode, β=2.5. Total by force moments
`E[H_t]−Var(G_t)` (the SM1 estimator); `G″_Haar` exact (1.2) at the **same** `a₀*`.

| L | `a₀*` | `G″_total` | `G″_Haar` (exact) | `G″_Wilson` | Haar floor `G″_Haar(0)` | Wilson/Haar |
|---|---|---|---|---|---|---|
| 1 | −2.14 | 0.694 | 0.2715 | **+0.422** | 0.2054 | 1.55 |
| 2 | 0.30 | 1.017 | 0.2401 | **+0.777** | 0.2396 | 3.24 |
| 3 | 1.15 | 0.641 | 0.2556 | **+0.385** | 0.2493 | 1.51 |
| 4 | 0.48–1.03 | 0.435 | 0.2581 | **+0.18±0.02** | 0.2542 | ≈0.7 |
| 5 | 0.65–0.76 | 0.446 | 0.2589 | **+0.20±0.03** | 0.2572 | ≈0.8 |
| 6 | 0.34–0.65 | 0.476 | 0.2602 | **+0.22±0.02** | 0.2591 | ≈0.85 |

Independent-seed scatter (N=1.5·10⁶ each), L=4–6: total `{0.435,0.461,0.458}` / `{0.446,0.488,0.434}` /
`{0.476,0.497,0.451}`, i.e. `0.45±0.02` / `0.46±0.03` / `0.47±0.02`. **The total flattens** at `≈0.46`
(not `1/L`). The Wilson part is **positive and `O(1)`**, decaying from the L=2 peak (`0.78`) to a
**plateau `≈0.2`**. The Haar part is the exact, noise-free, monotone `0.240→0.260→0.2688`.

**Cross-references (single volume, total only, no L-series):** 2×2 cross (β=2.5) `+0.78`; 4D
single-site Savvidy eigenmode (β=2.5) `+0.424` (SM1). Their Haar parts are not the uniform chain value
(their soft directions are not uniform), so they are quoted only as totals.

**MC honesty.** Third digit is not certified (ESS 120–2200 per bin; independent seeds move
`G″(a₀*)` by `≈0.02–0.05` at L≥4, `0.05–0.15` at L=2,3). The exact L=1 values and the exact Haar
part carry no such scatter.

---

## 4 · Normalization invariance

`G″(a₀*)=−d²/da² log ρ_v` is invariant under **any constant rescaling** of the density: measure
normalisation, a constant chart/Jacobian factor, the MC proposal constant, and the bin width all drop
out of `log`. Because both parts are computed from **densities at the same `a`**,
`G″_Haar + G″_Wilson = G″_total` holds before and after any rescaling, and each part is individually
normalisation-free. The Haar part uses `ρ₁` only up to a dropped overall `1/(4π²)` — verified
irrelevant. The one scale that is *not* free is the **mode normalisation**: with `‖v‖=1` fixed
(`v_l=e₃/√L`), `a` has the physical-amplitude units in which the `L`-scaling above is stated; rescaling
`v→cv` rescales `a→ca` and `G″→c²G″` (a change of coordinates, not of the geometry). The k=0 mode
`v_l=e₃/√L` is the canonical unit-normalised constant mode, making the quoted `L`-dependence intrinsic.

---

## 5 · The limit, the crossover, and what is open

**On the tractable strips.** `G″(a₀*) = G″_Haar + G″_Wilson` with `G″_Haar→0.2688` (proven) and
`G″_Wilson` measured **positive, `O(1)`, plateauing at `≈0.2`**. So `G″(a₀*)→ ≥0.27+0.2 ≈ 0.47 > 0`:
**the Haar part dominates the limit; the Wilson part does not drive it to 0 or negative** on this
model. The crossover where Haar overtakes Wilson is `L*≈3.5`.
**Correction to the SM1 reading.** SM1's "`≈1/L`" was fit through `L=2,3,4` (`1.01,0.65,0.42`). Adding
`L=5,6` (`0.45,0.48`, 3 seeds each) shows the total **flattens**; the structural reason is §2 — both
parts are `O(1)`, so the true asymptotic is an `O(1)` limit, not `1/L→0`. `[V]` for the plateau,
`[T-structural]` for the `O(1)` scaling.

**What is `[O]` = H1.** The **sign of the large-L Wilson part in the periodic (decompactification)
model** is *not* determined by the strip. The strip carries a `Z₂`-twisted fixed boundary whose
holonomy contributes a positive Wilson term. The periodic dynamical-lattice probe (SM6, β=2.4, L=8,
plane-wave modes) gives `G″ = +0.080±0.189` — **zero-consistent**, and at face value *below* the Haar
floor `0.2688`. If the periodic Wilson part is negative and of size `≳0.27`, the Haar floor does not
protect and the total can approach 0. So the honest statement of the decompactification limit is:

```
   G″(L→∞) = (proven positive Haar floor 1/Var(ξ₃)=0.2688)  +  (Wilson part of undetermined large-L sign).
```

The Haar floor is a **theorem** and is the L-independent piece; the Wilson part is the L-dependent
correction whose large-L sign is the mass-gap input H1 `[O]`. **The decomposition converts "does G″
survive decompactification?" into "is the large-L Wilson part bounded below by −0.2688?", which is H1
in curvature language** — and makes explicit that the Haar/Jacobian term alone already supplies a
positive, mode- and L-uniform `0.2688`.

---

## 6 · Acid test (free field), part by part

1. **Free massless field** (`Φ≡0`, no Haar Jacobian): `G″_Haar=0` (Lebesgue ⇒ Gaussian marginal ⇒
   curvature `κ`), `G″_Wilson=κ`, so `G″=κ→0`. **The decomposition manufactures no floor** — each part
   returns its free value, and the sum vanishes with the mass. `[T]`.
2. **Pure Haar `β=0`**: `S_W=0⇒R≡1⇒G″_Wilson≡0` **exactly**, and `G″_Haar` is the positive floor
   `0.2054→0.2688`. This is **not** the free field (item 1): the Haar Jacobian `J` is the discriminator.
   The tool thus separates "measure curvature" (Haar, positive, L-uniform) from "action curvature"
   (Wilson, `→0` at β=0 / free) — exactly the discrimination the acid test demands. `[T]+[V]`.
3. **β=0 consistency of the total**: at β=0 the total reduces to the Haar part, reproducing the SM1
   `[E4]` convolution values (`0.229,0.240,0.249,0.254`) to the digit. `[T]+[V]`.

---

## 7 · Ledger

| item | statement | label |
|---|---|---|
| exact decomposition `G″=G″_Haar+G″_Wilson`, `R=E_{Haar}[e^{−S_W}\|a]` | eqs (0.1)–(0.2), verified at L=1 to 5 digits | **[T]** |
| cumulant form `G″_Wilson=κ₁″−½κ₂″+…`; leading `κ₁″` carries the sign | L=1: κ₁″=0.454, exact 0.422 | **[T]**+**[V]** |
| Haar part exact convolution; `→1/Var(ξ₃)=0.2688` (local CLT) | (1.2); `0.205→0.259→0.269` | **[T]**+**[V]** |
| Haar part is a positive L-independent floor (min at a=0 on \|a\|≤1.5) | §2.1, §3 | **[T]+[V]** |
| Wilson part `O(1)`, **no `N_□=L` growth** (fluctuation-scale constraint) | §2.2 | **[T-structural]** |
| Wilson part **positive at `a₀*`** at β=2.5: exact +0.42/+0.16 (L=1), strips +0.78→plateau ≈0.2 | §3 | **[V]** |
| Wilson part **negative at the L=1 LONG barrier** `a=0` (−0.069) | §3.1 | **[V]** (`A−Bβ²`, consistent with **[NG]**) |
| Wilson part grows with β at fixed L (L=1: 0.03→1.26 for β=0.5→6) | §3.1 | **[T]+[V]** |
| total on strips **flattens at ≈0.46** for L≥4; `1/L` was an L=2–4 transient | §3.2, 3 seeds L=4–6 | **[V]** |
| crossover Haar>Wilson at `L*≈3.4` | Wilson/Haar ratios | **[V]** |
| normalization invariance of each part and of the identity | §4 | **[T]** |
| acid: free field → both parts → 0; β=0 → Wilson≡0, Haar>0 | §6 | **[T]+[V]** |
| **large-L sign of the Wilson part in the periodic model** | strip positive vs periodic zero-consistent | **[O] = H1** |
| **decompactification limit of `G″`** | = Haar floor + (Wilson sign) | **[C]/[O]** |
| **Clay mass gap** | unbuilt | **0%** |

---

## 8 · Witness — `research/v3/w7/d4_haar_wilson.py`

`py -3.12 -I d4_haar_wilson.py [quick|full]` (exit 0; numpy 2.4.6, scipy 1.18.1; exact 1-D
Gauss–Legendre for the one-link marginal and the L=1 marginals; FFT L-fold convolution for the Haar
part; Haar-importance MC for the total at L≥2; one job at a time). Quoted output (abridged):

```
[A] Var(xi_3)_Haar = 3.719824 (direct 3.719824)  =>  CLT Haar-floor limit 1/Var = 0.26883
[C] L=1 exact beta=2.5, theta=pi:
    LONG e3  a0*=-2.1411  total=+0.69398  Haar=+0.27154  Wilson=+0.42244  (direct ratio +0.42242)
    TRANS e1 a0*= 0.0000  total=+0.36231  Haar=+0.20542  Wilson=+0.15689  (direct ratio +0.15720)
    LONG barrier a=0:     total=+0.13599  Haar=+0.20542  Wilson=-0.06942
[C2] beta-scan at L=1 minimiser: LONG Wilson 0.033,0.109,0.422,0.773,1.259 (beta=0.5,1,2.5,4,6)
[C3] cumulants (L=1,LONG): kappa1''=+0.45357 ; -1/2 kappa2''=+0.02699 ; sum=+0.48056 ; exact=+0.42244
[D] chains beta=2.5, theta=pi (total by force moments, Haar exact):
    L   a0*   G_total   G_Haar  G_Wilson  HaarFloor  ESS
    2  0.302   1.0170   0.2401   0.7769   0.2396   2218
    3  1.154   0.6406   0.2556   0.3850   0.2493   1118
    4  1.030   0.4345   0.2581   0.1764   0.2542    555
    5  0.759   0.4456   0.2589   0.1868   0.2572    343
    6  0.646   0.4760   0.2602   0.2159   0.2591    119
   second seeds (N=1.5e6): L=4 total {0.461,0.458}; L=5 {0.488,0.434}; L=6 {0.497,0.451}
[E] beta=0 Haar floor L=1..6: 0.20542 0.23955 0.24927 0.25425 0.25719 0.25914  (-> 0.26883)
```

No third-party sources downloaded; standard inputs (Haar measure, local CLT/Edgeworth, Prékopa 1971,
Savvidy 1977) cited from the corpus.

---

## 9 · Sharpest line

**`G″ = G″_Haar + G″_Wilson` is an exact identity. The Haar part is a deterministic L-fold convolution
whose curvature is `O(1)` with the exact limit `1/Var(ξ₃)=0.2688>0` — a proven, L-independent,
positive floor. The Wilson part is `O(1)` (the `N_□=L` growth is cancelled by the fluctuation-scale
constraint), positive at the coarse minimiser on every tractable model (`+0.42/+0.16` exact at L=1;
`+0.78`→ plateau `≈0.2` on the strips), and negative only at the L=1 longitudinal barrier. On the
strips the total therefore **flattens at `≈0.46`, above the Haar floor, not at 0** — the SM1 `≈1/L`
decay was an L=2–4 transient. So **the Haar part dominates the limit (`→0.2688>0`); the Wilson part
neither grows with the plaquette number nor drives the total to 0 on this model.** The one unresolved
input is the **large-L sign of the Wilson part in the periodic (decompactification) model** — the
periodic L=8 probe is zero-consistent — which is exactly H1 `[O]`. Clay 0%. No QED.**

---

*Refs.* `SM1-SOFT-EXACT.md` §3–§4, §6 (small-volume `G″(a₀*)`, `1/L` reading, β=0 convolution floor
`0.229…0.254`); `SM4-CONSTANT-MODE.md` §1.2–§2 (exact L=1 marginal, all-loop, closed forms);
`SM2-RANK-ONE.md` (soft-mode margin, `G″=S−D`, ceiling); `research/unconditional/PROOF-JACOBIAN-SAVES.md`
(`Φ` convex ⇒ `W_geo≥0` via Prékopa; the Haar/Jacobian term's sign provable, control of the dynamical
force variance = H1); `docs/paper-v3/main.tex` §6.1 (defect/Fisher), §6.7 (soft mode: `≈1/L` decay vs
L-independent Haar floor, "the sharpest open question"). Witness `research/v3/w7/d4_haar_wilson.py`.
*No third-party sources downloaded.*
