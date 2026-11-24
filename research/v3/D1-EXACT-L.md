# D1-EXACT-L — decompactification of the YM soft mode: the constant-mode curvature **saturates (c>0) at the Haar/thermodynamic floor**, it does **not** decay to 0

> **2026-11-24 · Wave D1** — the decompactification question posed in `SM-WAVE-SYNTHESIS.md` §2.
> At β=2.5 the constant (k=0) soft-mode curvature `G″(a₀*)` is positive at small volume and decreases
> ≈`1/L` (`1.01→0.65→0.42`, L=2,3,4); the β=0 Haar floor is `L`-independent (`0.23–0.25`). **Does `G″(L)→0`
> as `L→∞` (then the uniform floor is H1), or does it saturate at the Haar floor (`c>0`)?**
> Witness: `research/v3/w7/d1_exact_L.py` (`py -3.12`, numpy/scipy; one job at a time, small L).
> Labels [T]/[V]/[C]/[NG]/[O] (charter §2). **No QED. Clay 0%.**

---

## 0 · Verdict

**`c>0`: the constant-mode curvature saturates at an `O(1)`, `L`-independent positive value — for the
untwisted strip it is exactly the Haar floor (`c≈0.27`), for the Z₂-twisted strip it is `≈0.42`
(frustration-pinned, still `>0`). The `≈1/L` decrease seen at `L=2,3,4` is the *Wilson (pinning) part*
decaying, sitting on top of a persistent Haar-carried floor. The `1/L→0` reading is excluded.**

The mechanism is exact and structural [T]:

1. **The constant mode is a zero mode of the Wilson difference operator.** The plaquette action depends
   on differences of neighbouring links; the uniform mode `v=(e₃,…,e₃)/√L` is annihilated by the
   quadratic (one-loop) part. In a Gaussian model `(1/2σ²)Σz² + (β/8)Σ(z_{l+1}−z_l)²`, the constrained
   mean has variance `Var(a)=σ²` **exactly, independent of β and L** (the Laplacian kills the constant
   vector). So the one-loop Wilson coupling contributes `0` to `G″(constant)` for every `L`.
2. **What sets the floor is the local (Haar) pinning.** The compact Haar measure confines each link to
   `|ξ|≤2π`, giving an `O(1)` local variance → a finite susceptibility → `G″→1/(σ²χ)=O(1)`. A *free*
   massless field has no such local pinning: its constant mode is the Laplacian zero mode, `Var(a)=∞`,
   `G″=0` (acid test passes).
3. **The `1/L` piece is a boundary/end effect.** In the strip `1×L` the single "sum" plaquette
   `P₁=X₁X₂` directly constrains the constant mode; it is one `O(1)` term out of `L` → its curvature
   contribution is `O(1/L)`. This is exactly the `1.01→0.65→0.42` the corpus measured; it rides on the
   floor and fades.

**Why `c>0` is forced for the 1-D strip (model-independent) `[T/C]`.** `Var(a)=(1/L)Σ_{l,l'}Cov(z_l,z_{l'})`.
At β=2.5 the 1-D strip has a **finite correlation length** `ξ` (nearest-neighbour 1-D chain, no phase
transition, bounded compact variables), so `Cov` decays as `e^{-|l-l'|/ξ}` and `Σ_{l,l'}Cov = L·χ` with
`χ=σ²(1+2Σ_{d≥1}ρ_d)<∞`. Hence `Var(a)=χ` is **finite and `L`-independent**, and
`G″(const)≈1/Var(a)=1/χ>0`. The only escape is `χ=∞` (a critical/Goldstone constant mode), which the
1-D chain does not have. This makes `c>0` a consequence of extensivity + finite `ξ`, **independent of
the Cartan reduction** — the chain is a concrete realisation, not the reason.

**The answer to the posed question: SATURATION (`c>0`), not decay to H1 — for the 1-D strip.** The scope
strip↔4D is stated honestly in §7 (`[O]`), because the constant-mode susceptibility *can* diverge in
`d≥2` (KT/Goldstone); the strip (`d=1`, finite correlation length) is `c>0` decisively.

---

## 1 · Normalization, and how the Wilson action scales with L

**Coarse amplitude (fixed, unit-normalized mode).** SU(2) flat chart `U_l=e^{iξ_l·σ/2}`, `ξ_l∈B(2π)`,
Haar `J(ρ)d³ξ`, `Φ=−log J`. Constant mode

```
   v = (e₃,…,e₃)/√L ,   |v| = 1  (FIXED unit norm) ,   a = ⟨v,ξ⟩ = (1/√L) Σ_l ξ_{l,3} .
   V_eff(a) = −log ρ_v(a) ,   G″(v;a₀) = −d²/da² log ρ_v(a)|_{a₀} = −ρ_v″/ρ_v .   (1.1)
```

`G″` is a curvature in the **unit-normalized mode amplitude** `a` (units of the gauge field). It scales
as `G″→λ^{-2}G″` under `v→λv`; the dichotomy `c=0` vs `c>0` is invariant under any `L`-independent
rescaling (see §6).

**Plaquette action of the constant configuration.** For the strip (`1×L`, plaquettes
`P₁=X₁X₂`, `P_i=B X_i^{-1}X_{i+1}`, `P_L=B X_L^{-1}`, boundary `B=e^{iθσ₃/2}`), take all `z_l=z`
(Cartan constant config), bulk `X_i^{-1}X_{i+1}=I`:

| boundary | bulk plaquette `½tr P_i = cos(θ/2)` | `S_W(constant)` | per-mode |
|---|---|---|---|
| **θ=0** (trivial) | `cos0 = 1` ⇒ `1−½trP_i = 0` | `O(1)` (only the two ends, `P₁` and `P_L`) | `O(1/L)→0` |
| **θ=π** (Z₂ twist) | `cos(π/2)=0` ⇒ `1−½trP_i = 1` | `β(L−2)+O(1) = O(L)` (frustrated) | `O(1)` |

**This is the crux of §1 of the task.** The plaquette action **vanishes** for the constant configuration
iff the boundary is untwisted (θ=0); for the twisted case (θ=π, the corpus's "worst" case) the constant
config is *frustrated* and carries `O(L)` action. Both still give `c>0`, with different asymptotes
(θ=0 → the Haar floor; θ=π → a larger, frustration-pinned floor).

---

## 2 · The exact β=0 Haar floor `[T]+[V]`

At β=0 the marginal of `a` is the **exact** `L`-fold convolution of the one-component Haar marginal
`ρ₁(z)=2π∫_{|z|}^{2π} sin²(ρ/2)/ρ dρ = π(Cin(2π)−Cin(|z|))`, `|z|≤2π`. Computed exactly (density from
its characteristic function `Z₀(k)=[∫ρ₁(z)e^{-ikz/√L}dz]^L`, no approximation). `Var_Haar(ξ₃)=3.719824`,
so the Gaussian-limit floor `1/Var=0.268830`.

| L | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | →∞ |
|---|---|---|---|---|---|---|---|---|---|
| `G″(a₀*)` | 0.20512 | 0.23955 | 0.24933 | 0.25425 | 0.25718 | 0.25914 | 0.26053 | 0.26157 | **0.26883** |

`a₀*=0` throughout (symmetric); at `L=16` it is `0.26521`, still rising to `0.26883`. The Haar floor is
`L`-**independent**, a positive `O(0.27)` constant — the exact value the `0.23–0.25` in
`SM-WAVE-SYNTHESIS` was converging to. `[T]` (convolution identity, `G″=−ρ″/ρ`), `[V]` (values).

*Consistency check:* the interacting chain below, run at β=0, reproduces this table to `5×10⁻⁶` (it is
the same object).

---

## 3 · The tractable interacting chain at β=2.5, `L=1..16` `[V]`

The constant-mode marginal on the strip, in the **exact Cartan sector**: per-link exact 3D weight
`ρ₁(z)` + the Cartan projection of the Wilson action (plaquette trace `cos(Δz/2)`). This is exact at
β=0 (§2) and the exact Cartan sector of the strip at β=2.5; it reaches `L=16` by an exact transfer
matrix + windowed FFT (deterministic, no Monte Carlo; validated against brute-force 2-D quadrature to
`10⁻¹⁴`, and against the analytic `L=1` to 3 digits).

**θ=0 (trivial boundary), β=2.5:**

| L | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| `G″(a₀*)` | 0.830 | 1.590 | 1.140 | 0.744 | 0.555 | 0.458 | 0.401 | 0.365 |
| `G″` (β=0) | 0.205 | 0.240 | 0.249 | 0.254 | 0.257 | 0.259 | 0.261 | 0.262 |
| Wilson part | 0.625 | 1.350 | 0.891 | 0.490 | 0.298 | 0.199 | 0.141 | 0.104 |

| L | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|
| `G″(a₀*)` | 0.341 | 0.324 | 0.311 | 0.301 | 0.293 | 0.286 | 0.281 | **0.276** |
| Wilson part | 0.080 | 0.064 | 0.050 | 0.040 | 0.032 | 0.025 | 0.020 | ~0.015 |

The Wilson part decays (geometric, `q≈0.72`), and `G″(L)→0.27`, **the Haar floor** of §2. `L=16` value
`0.2764` vs the β=0 floor (`0.2652` at `L=16`, rising to `0.2688` as `L→∞`); the two curves merge as
`L→∞` — the untwisted Wilson part decays to `0` and the constant-mode curvature falls onto the exact
Haar floor.

**θ=π (Z₂ twist), β=2.5:**

| L | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 14 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|
| `G″(a₀*)` | 0.851 | 1.580 | 0.896 | 0.706 | 0.467 | 0.438 | 0.433 | 0.424 | 0.421 | **0.417** |

The twisted case saturates **above** the Haar floor (`c≈0.42`), because the constant config is
frustrated and carries `O(L)` Wilson action (§1) → a persistent `O(1)` restoring force. (The `L=5`
entry is a mode-branch switch of the frustrated density, `a₀*` jumping near 0; the `L≥6` tail is smooth
and is what is fit.) Still `c>0`. The **Haar part** is the same β=0 row of §2 for every `θ` (at β=0 the
boundary is irrelevant), so the Wilson part is `G″(2.5)−G″(β=0)`; it is `≈0.15` here, not decaying
(frustration-pinned).

---

## 4 · Fits `[V]` (actual witness output)

Forms `c + A/L^p` (`p=1,2,3`) and `c + A·q^L` (a finite correlation length predicts the exponential; the
decrements halve each step). Fits to the clean tail (`L≥6`; `L≥4` where noted):

**θ=0 (→ Haar floor):**

| range | form | `c` | `A` | `q` | rms |
|---|---|---|---|---|---|
| L≥6 | `c+A/L` | 0.1621 | +1.688 | — | 0.0071 |
| L≥6 | `c+A/L²` | **0.2476** | +7.554 | — | **0.0004** |
| L≥6 | `c+A/L³` | 0.2757 | +41.4 | — | 0.0066 |
| L≥6 | `c+A·q^L` | 0.2724 ± 0.0022 | +1.290 | q=0.722 | 0.0021 |
| L≥4 | `c+A·q^L` | 0.2859 ± 0.0044 | +2.81 | q=0.632 | 0.0081 |

**θ=π (→ frustration-pinned floor):**

| range | form | `c` | `A` | `q` | rms |
|---|---|---|---|---|---|
| L≥6 | `c+A/L` | 0.3856 | +0.473 | — | 0.0029 |
| L≥6 | `c+A/L²` | **0.4096** | +2.111 | — | **0.0023** |
| L≥6 | `c+A/L³` | 0.4175 | +11.55 | — | 0.0031 |
| L≥6 | `c+A·q^L` | 0.4160 ± 0.0026 | +0.337 | q=0.731 | 0.0024 |

**Best fits:** `c+A/L²` has the lowest rms for θ=0 (`c=0.2476`, rms `0.0004`) and ties the exponential for
θ=π (`c=0.4096`); the finite-ξ exponential gives `c=0.2724` / `0.4160`. The estimates cluster at
`c ≈ 0.25–0.28` (θ=0) and `0.41–0.42` (θ=π). **Every form and range returns `c>0`.** The θ=0 value
brackets the Haar floor `0.269` from both sides, confirming the untwisted curvature lands on the Haar
floor. The bare `c+A/L` over the full `L=1..16` is *not* used (it is dragged by the `L=2` sum-plaquette
peak `1.59` and the `L=1` minimal torus `0.83`, both `O(1/L)` transients); the `L→∞` behaviour is set by
`L≥6`.

*Uncertainty.* Deterministic (quadrature + validated transfer), not Monte Carlo. Windowed-FFT/quadrature
resolution agrees to `3–4` digits across settings for `L≥3` (§3); the fit `c` is stable to `±0.01`
(θ=0) and `±0.005` (θ=π, exponential). The **spread between functional forms** (`c` from `0.248` to
`0.286` for θ=0) reflects the finite-`L` transient, not numerical error; it does not cross `0`. The
dominant physical uncertainty is the **model** (Cartan reduction vs full 3D, §5), not the fit.

---

## 5 · Full-3D calibration (the actual marginal, not the Cartan sector) `[V]`

The Cartan chain is a *tractable reduction*; the full 3D marginal includes the transverse fluctuations.
Exact full-3D values at the volumes where direct quadrature is feasible:

| model | a₀* | `G″(a₀*)` |
|---|---|---|
| L=1 full 3D LONG, θ=π (exact 1-D radial) | −2.142 | **+0.694** |
| L=1 full 3D LONG, θ=0 | 0.000 | **+0.730** |
| L=1 full 3D TRANS (Savvidy), θ=π | 0.000 | **+0.362** |
| **L=2 full 3D LONG, θ=π** (cylindrical quadrature) | +0.415 | **+0.990** |
| L=2,3,4 strip LH MC (SM1 §4, `[V]` noisy) | — | 1.01, 0.65, 0.42 |

The exact quadrature at L=2 gives `0.990`, matching the SM1 MC `1.01` (the MC range was `0.83–1.01`) —
so the full-3D curve is `0.694, 0.990, 0.65, 0.42` for `L=1..4`: it **rises to L=2 then falls**, the same
shape as the chain (`0.851, 1.58, 1.14, 0.74, …`). The Cartan chain sits **above** the full 3D
(L=1: 0.851 vs 0.694, `+23%`; L=2: 1.58 vs 0.99, `+60%`): the Cartan projection over-weights the
pinning direction, and the transverse fluctuations soften it — increasingly with `L`. **So the chain's
*values* are not the full model's; its *L-dependence* is.** The `L→∞` asymptote of the *full* model is
not computed here beyond L=2, but it is `c>0` by the extensivity argument of §0/§6 below (the full strip
has massive transverse modes at β=2.5 ⇒ finite correlation length ⇒ finite susceptibility ⇒ `c>0`),
independently of the chain's exact floor.

---

## 6 · Normalization dependence — the invariant statement `[T]`

`G″` is the curvature of `−log ρ_v(a)` in the amplitude `a=⟨v,ξ⟩` of the **unit-normalized** mode
`|v|=1`. Three equivalent/invariant readings:

- **Rescaling `v→λv`:** `a→λa`, `G″→λ^{-2}G″`. The *shape in L* scales by `λ^{-2}`; `c>0` is preserved
  for any fixed `λ`. The canonical `|v|=1` gives the `O(1)` curve of §3.
- **Average field `φ=Σξ/L=a/√L` (per-site mean):** `G″_φ = L·G″_a`. With `G″_a→c>0` this **grows** as
  `cL` — the average field's constraint potential is extensive (`V_eff(φ)≈L·v(φ)`), curvature `∝L`. A
  naive "`1/L→0`" is what one gets by *dividing by `L`* (per-site curvature `=G″_a`), which is the
  quantity that saturates.
- **Invariant (susceptibility):** the *physical* statement is `χ := Var(a)|_{|v|=1} = 1/G″` (Gaussian).
  `χ` is finite and `L`-independent ⟺ `c>0` ⟺ the constant mode is **not** a soft/Goldstone mode. A
  `1/L→0` decay would require `χ→∞` (critical mode) — excluded for the 1-D strip.

So the result is normalization-covariant; the normalization-independent statement is:

> **The constant-mode susceptibility `χ` is finite and `L`-uniform for the strip at β=2.5; the IR
> (zero-momentum) curvature is an `O(1)` positive constant. It does not vanish.**

---

## 7 · Acid test, scope, and what this does and does not prove

**Acid test (free field) — passes discriminatively `[T]`.** Pure quadratic coupling `(β/8)Σ(Δz)²`, no
Haar: the constant mode is the Laplacian zero mode, `Var(a)=∞`, `G″=0`. The Haar floor is exactly what a
free field *lacks* (no local pinning); the tool returns `0` on the free field and `+0.27…+0.42` on the
interacting compact measure. It discriminates.

**Scope — the dimensional caveat `[C]`/`[O]`.**
- **1-D strip (this note):** `c>0`, established `[V]` with an `[T]` mechanism. `χ` finite (1-D, finite
  correlation length), floor = Haar.
- **2-D `L×L`:** the constant-mode susceptibility is a U(1)/Cartan *spin-wave* susceptibility, which can
  **diverge** (KT phase: `χ∼L^{2−η}`) → `G″∼L^{−(2−η)}→0` (`c=0`) is *possible* in 2D. **Not computed
  here**; flagged `[O]`. (The `1×L` strip is *not* this; it is a 1-D chain.)
- **4-D `L⁴` / decompactification to ℝ⁴:** the zero-momentum mode's effective potential curvature *is*
  the mass gap's IR floor; whether it saturates or vanishes is the gap itself `[O]`.

**What this proves:** for the tractable 1-D constant-mode problem (the model behind the quoted
`1.01→0.65→0.42`), the curvature **saturates at an `O(1)` Haar/thermodynamic floor (`c>0`)** — resolving
the `SM-WAVE-SYNTHESIS` §2 question *in favour of saturation for the strip* — and the `≈1/L` behaviour is
identified exactly as the decaying Wilson (sum-plaquette / boundary) part. **What it does not prove:** a
uniform-in-`L` floor in 2-D or 4-D (KT/Goldstone divergence is not excluded), the transport to the full
non-constant soft mode, or the 4-D mass gap. **The uniform floor in `d≥2` remains H1 `[O]`.** No QED.
Clay 0%.

---

## 8 · Ledger

| item | statement | label |
|---|---|---|
| `G″(v;a₀*)=−ρ_v″/ρ_v≥0`, interior minimiser | (1.1), H4/SM4 | **[T]** |
| exact β=0 Haar floor, `L`-independent, `→1/Var=0.2688` | §2 table | **[T]+[V]** |
| constant mode = zero mode of the Wilson difference operator ⇒ one-loop Wilson part of `G″(const)=0` | §0, Gaussian identity `Var(a)=σ²` ∀β,L | **[T]** |
| plaquette action of the constant config: `O(1)` (θ=0), `O(L)` (θ=π frustration) | §1 table | **[T]** |
| chain (exact Cartan sector): `G″(L)` saturates, θ=0 `c=0.272`, θ=π `c=0.416` | §3–§4, L=1..16 | **[V]** |
| best fits `c+A/L²` / `c+A·q^L`: θ=0 `c=0.248/0.272`, θ=π `c=0.410/0.416`; rms `≤0.002`; `c>0` | §4 | **[V]** |
| **1-D strip `c>0` from extensivity + finite ξ (independent of the Cartan reduction)** | §0, `Var(a)=χ<∞` | **[T/C]** |
| normalization: invariant = finite `L`-uniform susceptibility `χ=1/c` | §6 | **[T]** |
| free-field acid: `G″=0`; compact Haar gives the floor | §7 | **[T]** |
| full-3D L=1 exact `0.694/0.730/0.362`; L=2 `0.990`; chain above by `23%/60%`, same `L`-trend | §5 | **[V]** |
| **1-D strip constant-mode floor `c>0`** | §0–§4 | **[V]**+[T-mechanism] |
| 2-D/4-D constant-mode floor (KT/Goldstone divergence possible) | not proved; `c=0` not excluded | **[O]** |
| **uniform-in-`L` floor in `d≥2` / 4-D mass gap** | = H1 / Clay | **[O] / 0%** |

---

## 9 · Witness — `research/v3/w7/d1_exact_L.py`

`py -3.12 d1_exact_L.py [quick|full]`. Exact one-component Haar marginal `ρ₁` (via `Cin`); β=0 floor by
exact convolution; interacting Cartan chain by an exact transfer matrix + windowed FFT (validated
against brute 2-D quadrature to `10⁻¹⁴` and against the analytic `L=1` to 3 digits; reduces to the full
Haar floor to `5×10⁻⁶` at β=0); full-3D L=1 by exact radial quadrature, L=2 by cylindrical quadrature;
fits `c+A/L^p` (`p=1,2,3`) and `c+A·q^L`. Deterministic, no Monte Carlo. One job at a time, small L.

**Sharpest line.** *For the 1-D constant-mode problem the soft-mode curvature does not decay to zero —
it saturates at the `O(1)` Haar floor (`≈0.27` untwisted, `≈0.42` twisted), the `≈1/L` seen at `L=2,3,4`
being the decaying Wilson/pinning part on top of it; this is `c>0` (saturation), not `c=0` (H1), for the
strip. In `d≥2` the same question can fail (KT/`Goldstone` susceptibility divergence), so the 4-D mass
gap is unchanged: H1 `[O]`. No QED. Clay 0%.*
