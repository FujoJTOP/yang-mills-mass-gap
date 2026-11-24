# SM4-CONSTANT-MODE — the exact constant (softest) mode effective potential at β=2.5: a genuine positive floor for the constant mode `[T]+[V]`, in **both** the Cartan and the off-diagonal (Savvidy) class; the one-loop determinant is a free-gluon `1/L²` object and is **not** the carrier; the uniform-in-`L` bound is still H1 `[O]`

> **2026-11-24 · SM4 — the last open mode of the Yang–Mills UV face.**
> The soft mode on a small torus *is* the **constant** mode. Its effective potential is the
> constraint effective potential `V_eff(a)=−log ρ_v(a)` of the coarse amplitude `a=⟨v,ξ⟩`,
> and on the **minimal torus** (one plaquette / one free link) it is **exactly** the
> finite-dimensional marginal — Haar Jacobian + Wilson action, no non-constant modes — so the
> finite-dimensional integral *is* the complete (all-loop) answer. This note derives that exact
> marginal in closed form for the two **distinct** constant-mode classes (Cartan/holonomy `v=e₃`
> and off-diagonal/Savvidy-tachyon `v=e₁`), computes `G″(a₀*)` at `β=2.5`, maps its `L`,
> boundary and twist dependence, and decides the soft-mode margin.
> Reads: `research/v3/H4-HEAD-COUNTEREXAMPLE.md`, `H2-DEFECT-BOUND.md`,
> `research/unconditional/K0-BLIND.md`, `K0-CONSTANT-MODE-PROTECTION.md`, `k0-EXACT-ONEMODE.md`,
> `k0-HOLONOMY-EFFPOT.md`, `k0-CONSTRAINT-EFFPOT.md`, `docs/paper-v2/ch04-gap.tex`.
> Witness: `research/v3/w7/sm4_constant_mode.py` (`py -3.12 [quick|full]`, exit 0; output quoted below;
> exact 1-D Gauss–Legendre `n=4000`, `L≥2` via the H4-validated MC, imported).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[V]** numerical
> (deterministic quadrature; exact at `L=1`; MC at `L≥2`) · **[NG]** no-go for a named class ·
> **[O]** exact missing step.

---

## 0 · Verdict (no QED; Clay 0%)

**The constant mode has a genuine, strictly positive effective-potential curvature at its coarse
minimiser, in every class and at every `L` examined — and it is largest exactly where the Savvidy
tachyon lives.**

1. **[T] Sign at the minimiser is an unconditional theorem.** For the constant mode, the flat-chart
   Haar edge `J(2π)=0` forces the marginal maximum `a₀*=argmax ρ_v` to be **interior**, so
   `ρ_v'(a₀*)=0`, `ρ_v''(a₀*)<0`, hence
   ```
      G″(v;a₀*) = −d²/da² log ρ_v|ₐ₀* = −ρ_v''/ρ_v ≥ 0 ,   strict >0 for a non-degenerate max,  (0.1)
   ```
   for every `L`, boundary/twist, and every **constant** linear direction `v` (Cartan **and**
   off-diagonal). This is `H4-HEAD-COUNTEREXAMPLE` §1 restricted to the constant mode; SM4 supplies the
   exact closed-form density and extends the *measurement* to the off-diagonal class (`H4` scanned Cartan
   `e₃` only).
2. **[V] Value at `β=2.5` — positive, `O(0.4–1.0)`, no softening.** Exact minimal-torus `(L=1)`,
   worst twist `θ=π`: `G″(a₀*) = +0.694` (Cartan), `+0.362` (off-diagonal). Blocks `L=1,2,3`
   (Cartan / off-diagonal): `{0.694, 1.012, 0.646}` / `{0.362, 0.867, 0.427}`. **The softest
   constant-mode margin at `β=2.5` is `≈+0.36`** (off-diagonal, `L=1`); all values `> 0`, all above the
   `β=0` Haar floor `N/12=0.167`.
3. **[NG] The one-loop determinant of non-constant modes is *not* the carrier.** On the minimal torus
   there are no non-constant modes (the integral is exact). On a larger torus the non-constant Gaussian
   determinant is the GPY/holonomy **free-gluon** factor: nonzero for the *free* theory, of
   compactification scale (`m²∼c/L²→0` in the decompactification sense), and in the electric/center
   channel the neutral `0⁺⁺` flux is blind to. It does **not** supply the physical floor (banked
   `k0-HOLONOMY-EFFPOT`). The floor is carried by the **Haar+Wilson constraint potential** (`W_geo≥0`
   plus the dynamical part), which is `O(1)` and does not decay with `L`.
4. **What this proves, and only this.** A positive floor for the **constant mode on a finite small
   torus** (all `L`, all boundary/twist, both classes) — the **exact** face of the soft mode. It does
   **not** cover: the full soft mode (the non-constant `k→0` modes), a **uniform-in-`L`** lower bound
   `δ_W` (that is **H1 `[O]`**), or the decompactification `small-L→ℝ⁴` (adiabatic continuity = the
   gap itself, `[C]`/circular).

**Clay 0%. No QED. The wall did not move; the constant-mode face of it is now computed exactly.**

---

## 1 · The exact constant-mode marginal `ρ_v(a)` — Haar Jacobian + Wilson, and the one-loop determinant's true role

### 1.1 Setup (flat exponential chart, SU(2))

Per link `l`: `ξ_l∈B(2π)⊂ℝ³`, `U_l=exp(iξ_l·σ/2)`, Haar `dU_l=J(ρ_l)d³ξ_l`,
```
   J(ρ)=sin²(ρ/2)/ρ² ,   J(0)=1/4 ,   J(2π)=0 ,   Φ₀(ρ)=−log J(ρ)=2log(ρ/2)−2log sin(ρ/2) ,   Φ₀″(0)=N/12=1/6 .
```
Wilson `S_W=βΣ_□(1−½trP_□)`. The **constant mode** amplitude along a uniform unit direction `v` is
`a=⟨v,ξ⟩`; the exact marginal and the coarse curvature are
```
   ρ_v(a) = ∫ ∏_l J(ρ_l)d³ξ_l  δ(⟨v,ξ⟩−a)  e^{−S_W(ξ)} = ∫ δ(⟨v,ξ⟩−a) e^{−(S_W+Φ)} ∏d³ξ_l ,      (1.1)
   G″(v;a₀) = −d²/da² log ρ_v(a)|_{a₀}   (ch05 eq:split / G4 object (C)) .                          (1.2)
```

### 1.2 Minimal torus: one plaquette, one free link — the *complete* one-loop-and-all answer

Free `U`, boundary holonomy `B=e^{iθσ₃/2}`, `P=UB`:
```
   ½tr P = cos(ρ/2)cos(θ/2) − sin(ρ/2)sin(θ/2)(n·b̂) ,   n=ξ/ρ, b̂=e₃ ,
   S_W = β[1 − cos(ρ/2)cos(θ/2) + sin(ρ/2)sin(θ/2) n₃] .                                           (1.3)
```
On this torus every mode is constant, so (1.1) **is** the full determinant; there is nothing to
"loop-expand". The two constant-mode **classes** differ by which color component is frozen:

**(L) Cartan / holonomy direction `v=e₃`.** Constrain `ξ₃=a`, so `n₃=a/ρ` is φ-independent and the
azimuth about `e₃` gives `2π`:
```
   ρ_long(a) = 2π ∫_0^{√(4π²−a²)} r dr e^{−V} ,
   V = β[1 − cos(ρ/2)cos(θ/2) + sin(ρ/2)sin(θ/2)(a/ρ)] + Φ₀(ρ) ,   ρ=√(a²+r²) .                    (1.4)
```
(This reproduces H4's exact `rho_v_1x1`, `vdir='e3'`.) At `a₀=0`, `G″(0)=A−Bβ²` with
`A=E_Haar[Φ₀′(ρ)/ρ]=0.20511530`, `B=Var_Haar[sin(ρ/2)/ρ]=0.01105929`, `β*=4.30661` — the V3/G4/STAR2
counterexample ([NG] for `β>β*` **at the barrier**).

**(T) Off-diagonal / Savvidy-tachyon direction `v=e₁`.** Constrain `ξ₁=a`; now `n₃=ξ₃/ρ=(r sinψ)/ρ`
depends on the azimuth `ψ`. Because `ρ=√(a²+r²)` is `ψ`-independent, the inner `ψ`-integral is a
**modified Bessel function** in closed form:
```
   ρ_trans(a) = 2π ∫_0^{√(4π²−a²)} r dr J(ρ) e^{−β[1−cos(ρ/2)cos(θ/2)]} I₀( β sin(ρ/2)sin(θ/2) r/ρ ) .  (1.5)
```
This is the exact constant-mode marginal in the off-diagonal color direction — the **lattice
constant-mode home of the homogeneous Savvidy tachyon** (the continuum tachyon being the `F_B`-coupled
transverse-magnetic lowest Landau level; `e₁` is its `k=0` block analog). It is the physically
dangerous face, and it is *not* the Cartan `e₃` slice H4 scanned. (Its precise continuum identification
is `[C]`; the curvature computation below is `[T]+[V]` for the lattice constant mode regardless.)

### 1.3 Why `G″` is computable in closed form

Because `J(2π)=0`, the integrand of (1.4)/(1.5) vanishes at the moving endpoint `r=√(4π²−a²)` (where
`ρ=2π`). So `d/da∫_0^{R(a)}` carries **no boundary term**, and `G″=−d²log ρ/da²` is an ordinary
derivative of the smooth integral — evaluated by exact Gauss–Legendre quadrature + 5-point second
difference (error `<10⁻⁶`, checked against the analytic `A−Bβ²` at `a₀=0`, §6). The edge vanishing is
also *why the minimiser is interior*: a positive analytic density vanishing at both chart ends attains
its maximum strictly inside, giving (0.1).

### 1.4 The one-loop determinant's true role

The phrase "one-loop determinant" splits into two different objects, and only one is a floor:

| object | what it is | status for the gap |
|---|---|---|
| **(i) minimal-torus exact marginal** (1.4)/(1.5) | full finite-dim. integral = all-loop on `L=1` | **[T]+[V]** — the floor |
| **(ii) non-constant Gaussian det. on `L>1`** | GPY/holonomy **free-gluon** determinant around a fixed background | **[NG]** as a gap tool |

Object (ii) is nonzero for the **free** theory (`β=∞`: it is the free determinant), of compactification
scale (`V_1loop(Ω)=(1/L⁴)f(Ω)`, `Ω=aL` dimensionless ⇒ `∂²V/∂a²=(1/L²)f″(Ω)`), and lives in the
electric / center channel (the neutral `0⁺⁺` magnetic flux is blind to it — `PROOF-1FORM-ANOMALY`). It
therefore **cannot** manufacture the physical mass floor and it *vanishes* in the decompactification
sense. This is exactly the banked `k0-HOLONOMY-EFFPOT` verdict (`V″=c/L²`, free-gluon, wrong operator),
and it is why SM4 attributes the floor to the **Haar+Wilson** constraint potential (object (i)), whose
`W_geo≥0` is `O(1)` and `L`-uniform, not to the one-loop holonomy piece. (It is also `k0-EXACT-ONEMODE`'s
`G″=S−D`: the Gaussian/one-loop piece is the `D=0` slice; the exact finite-dim. marginal carries the
non-Gaussian `D`, which is the sign-relevant part.)

---

## 2 · `G″(a₀*)`: the sign theorem and its value at `β=2.5`

### 2.1 The sign is `[T]`

`a₀*` is interior (edge vanishing, §1.3), so at an interior maximum `ρ_v'(a₀*)=0`, `ρ_v''(a₀*)≤0`, and
```
   G″(v;a₀*) = (ρ_v'/ρ_v)² − ρ_v''/ρ_v = −ρ_v''/ρ_v ≥ 0 ,
```
strictly `>0` unless the top is quartic-flat (non-generic; the `=0` edge is the `δ→0` boundary of H1).
This holds for every constant `v`, every `L`, every boundary/twist, `β` — **[T]**, the H4 local sign
theorem specialised to the constant mode. Nothing here is numerical.

### 2.2 The values are `[V]` (exact at `L=1`)

Exact minimal-torus marginal, `β=2.5` (`N/12=0.1667`), twist scan `θ∈{0,π/4,π/2,3π/4,π}`:

| `θ` | LONG `a₀*` | LONG `G″(a₀*)` | LONG `G″(0)` | TRANS `a₀*` | TRANS `G″(a₀*)` | TRANS `G″(0)` |
|---|---|---|---|---|---|---|
| 0 | 0.000 | **0.730** | 0.730 | 0.000 | **0.730** | 0.730 |
| π/4 | −0.592 | **0.730** | 0.681 | 0.000 | **0.701** | 0.701 |
| π/2 | −1.162 | **0.730** | 0.543 | 0.000 | **0.618** | 0.618 |
| 3π/4 | −1.691 | **0.723** | 0.341 | 0.000 | **0.499** | 0.499 |
| π | −2.142 | **0.694** | 0.136 | 0.000 | **0.362** | 0.362 |

Every entry is `>0`. The **softest** constant-mode value at `β=2.5` is the off-diagonal `L=1`
`θ=π` entry, `G″(a₀*)=+0.362` — still `> N/12` and `>0` by a factor `>2`. The `β`-scan (`θ=π`) shows the
minimiser value **growing** with `β` and only the *barrier* `G″(0)` crossing zero:

| `β` | LONG `G″(a₀*)` | LONG `G″(0)` | TRANS `G″(a₀*)` |
|---|---|---|---|
| 0 | 0.205 | 0.205 | 0.205 |
| 2.5 | **0.694** | **+0.136** | **0.362** |
| 4.31 | 1.150 | **−0.0003** (`=β*`) | 0.543 |
| 5 | 1.326 | −0.071 | 0.612 |
| 10 | 2.595 | −0.901 | 1.118 |

So at `β=2.5` **both** the barrier and the minimiser are positive in **both** classes. The [NG]
strong-form violation (`G″(a₀=0)<0`) appears only for `β>β*=4.307` at the **longitudinal barrier**; at
the **minimiser** `G″(a₀*)` never goes negative (`[T]`), and the **off-diagonal** `G″` is positive at
the barrier too, for every `β` in the scan (`0–10`).

### 2.3 Comparison with the H4 local sign theorem

- **Agreement on the theorem.** H4 §1's `G″(a₀*)=−ρ_v''/ρ_v≥0` is `[T]` for *any* linear `v`; SM4's
  exact densities (1.4)/(1.5) are two instances, and the computed values obey it exactly.
- **New measurement — the dangerous class.** H4 scanned Cartan directions (`e₃`) only. SM4 adds the
  **off-diagonal** (`e₁`, Savvidy-tachyon) class, which by K0-CONSTANT-MODE-PROTECTION is "the tachyon's
  home": its constant-mode curvature is `+0.362` at `β=2.5` — **positive but the smallest of all
  classes** (`≈½` the Cartan value). That it is positive here is the sharp statement: even the
  homogeneous tachyon direction has a positive *coarse* curvature at its minimiser, because the
  coarse curvature is measure-carried (`D`-defect, `H2`), not a bare Hessian eigenvalue.
- **Scope unchanged.** The only surviving open content is the **uniform** `δ_W>0` (`H1`), exactly as H4
  concluded.

---

## 3 · `L`-, boundary-, and twist-dependence

**(i) `L`-dependence — no `1/L²` softening.** Constant (`k=0`) mode, `β=2.5`, `θ=π`; `L=1` exact,
`L=2,3` H4-validated MC:

| `L` | LONG `G″(a₀*)` | TRANS `G″(a₀*)` |
|---|---|---|
| 1 | **0.694** (exact) | **0.362** (exact) |
| 2 | **1.012** (MC) | **0.867** (MC) |
| 3 | **0.646** (MC) | **0.427** (MC) |

The constant-mode curvature is `O(0.4–1.0)` and **does not decay as `1/L²`** — the hallmark of a finite
constraint-potential quantity, *not* of the holonomy one-loop piece. (Contrast the one-loop free-gluon
`V″=c/L²`, which would fall by `×9` from `L=1→3`; the observed change is `≈×1`, LONG `0.694→0.646`,
TRANS `0.362→0.427`.)

**(ii) Boundary / twist.** `θ` is the boundary holonomy angle; `θ=π` (`B=−I∈Z₂` center) **is the `Z₂`
twist**, `θ=0` is trivial boundary. At `β=2.5` the minimiser value decreases monotonically
`θ:0→π` — LONG `0.730→0.694`, TRANS `0.730→0.362` — and the barrier LONG `G″(0) 0.730→0.136`. The
twist **reduces** the margin but never removes it; the worst case (off-diagonal, `θ=π`, `L=1`) is
`+0.362`. A non-abelian boundary (`b=(π/2,0,π/2)`) gives `G″(a₀*)≈+1.06` in the `L=2` block (H4 §3.2),
above the abelian worst case. So **twist is the minimizing boundary**, and even there the floor holds.

**(iii) `β`.** Reported in §2.2: the minimiser value increases with `β`; the floor is weakest at small
`β`, approaching (from above) the `β=0` Haar floor `N/12`.

**Reading.** The constant-mode floor is `L`-stable, twist-robust, and grows with coupling. It behaves
like a genuine constraint-potential (measure/entropic) quantity, not like a Casimir/one-loop artifact.

---

## 4 · What the constant-mode potential proves — the theorem and its exact scope

**Theorem SM4 (constant-mode positive floor).** *For SU(2) lattice gauge theory on a finite torus with
the constant-mode direction `v` (uniform in `x`), in the flat exponential chart with Haar measure, for
any block size `L`, any fixed boundary/twist `B`, any `β`, and any constant linear `v` (Cartan or
off-diagonal): let `a₀*` be the mode of the exact constant-mode marginal `ρ_v`. Then*
```
   G″(v;a₀*) = −d²/da² log ρ_v|_{a₀*} = −ρ_v''(a₀*)/ρ_v(a₀*) ≥ 0 ,   strictly >0 for a non-degenerate max.  [T]
```
*At `β=2.5` the value is `+0.694` (Cartan) and `+0.362` (off-diagonal, Savvidy class) on `L=1`, and
`O(0.4–1.0)` for `L=1,2,3`; the soft constant-mode margin is `≈+0.36`. `[V]`*

**What it covers:**
- The **constant mode only** — the `x`-independent collective coordinate, the softest mode on a small
  torus. Both its physically distinct classes: the **Cartan/holonomy** direction and the
  **off-diagonal / homogeneous-Savvidy-tachyon** direction.
- A finite **small torus**, all `L` examined, all boundary/twist, all `β`; the sign is `[T]`, the values
  `[V]` (exact at `L=1`).
- It is the exact **constraint effective potential** face — the right object (`k0-CONSTRAINT-EFFPOT`):
  `V_eff″(a₀*)=G″` is an identity, and the positive floor is the measure-carried (`W_geo` + dynamical)
  part, not a bare Hessian eigenvalue.

**What it does NOT prove (honest scope):**
1. **Not the full soft mode.** The non-constant `k→0` modes (the coarse mode with a long-wavelength
   profile) are *not* the constant mode; the constant-mode floor does not transport to them. That
   transport is exactly the `D`-defect / uniform Fisher bound — **H1 `[O]`**.
2. **Not a uniform-in-`L` floor.** The values `[V]` are positive but their *uniform* lower bound
   `δ_W=inf_{L,boundary,v}G″(a₀*)>0` is unproved — **H1 `[O]`**. (The `1/L²` free-gluon piece is the
   wrong object and gives no help; §1.4.)
3. **Not decompactification.** Reading small-`L` as a proxy for `ℝ⁴` is **adiabatic continuity**, which
   *is* the gap (`k0-HOLONOMY-EFFPOT`; `[C]` at best, circular at worst).
4. **[NG] the strong form is false at the barrier.** For `β>β*=4.307`, `G″(a₀=0)<0` in the longitudinal
   class — a genuine `[NG]` counterexample (V3/G4/STAR2). At `β=2.5` it is not in this regime
   (`G″(0)=+0.136`), so the *minimiser* theorem is not threatened; only the "uniform in `a₀`" form dies.

So SM4 delivers the **exact** face of the soft mode — a positive floor, computed, in both classes, with
its `L`/twist/`β` map — while pinning the residue to the pre-existing **H1**: the uniform-in-`L`
constant-mode (and full soft-mode) margin.

---

## 5 · Acid test and ledger

**Acid test (free field) — passes discriminatively.**
- **Free massless Gaussian** (`Φ≡0`, `S=κ|ξ|²/2`): the marginal of `a=⟨v,ξ⟩` is `N(0,1/κ)`, so
  `G″(a₀)=κ` at *every* `a₀`, and `κ→0` as the mass →0. **No manufactured floor** `[T]`. (The flat
  chart's Haar Jacobian, absent for the free Gaussian, is the discriminator.)
- **Pure Haar `β=0`:** the Haar marginal curvature is `+0.2051 > N/12=0.1667` — the Haar floor, *not*
  the free field `[T]+[V]`.

| item | statement | label |
|---|---|---|
| exact constant-mode marginal `ρ_v(a)`, Haar Jacobian + Wilson; Bessel reduction (1.5) for off-diagonal | eqs (1.4)–(1.5); one 1-D quadrature, `n=4000` | **[T]** (form) + **[V]** (values) |
| `G″(v;a₀*) = −ρ_v''/ρ_v ≥ 0` at the interior minimiser, all `L`/boundary/`β`/constant-`v` | §2.1, H4 §1 specialised | **[T]** |
| `β=2.5` values: LONG `+0.694`, TRANS `+0.362` (`L=1` exact); `L=2,3` ∈ `[0.36,1.01]` | §2.2, §3 | **[V]** |
| `L`-dependence: no `1/L²` softening (constraint-potential quantity) | `{0.694,1.012,0.646}` / `{0.362,0.867,0.427}` | **[V]** |
| twist `θ=π` (`Z₂`) is the minimizing boundary; floor survives (`+0.362`) | §3(ii) | **[V]** |
| one-loop non-constant determinant = GPY free-gluon, `c/L²`, electric/center channel | §1.4; banked `k0-HOLONOMY-EFFPOT` | **[NG]** as gap tool |
| `[NG]` strong form false: `G″(a₀=0)<0` for `β>β*=4.307` (longitudinal barrier) | §2.2 (V3/G4/STAR2) | **[NG]** |
| **uniform** `δ_W=inf G″(a₀*)>0`; full soft mode; decompactification | not proved | **[O] = H1** |
| free Gaussian `G″=κ→0`; β=0 Haar `G″=0.205>N/12` | §5 | **[T]+[V]** |

---

## 6 · Witness — `research/v3/w7/sm4_constant_mode.py`

`py -3.12 sm4_constant_mode.py [quick|full]` (exit 0). Exact 1-D Gauss–Legendre (`n=4000`) for the
minimal-torus marginal (1.4)/(1.5) with the closed-form inner Bessel integral; `G″` by 5-point second
difference of `log ρ` (validated `≤10⁻⁶` vs the analytic `A−Bβ²` at `a₀=0`, reproducing `β*=4.30661`).
`L=2,3` values by importing the H4-validated Haar-proposal MC (`scan_block`/`summarize`). Quoted output
(`quick`, `N=1.2×10⁶` for `L≥2`):
```
[E1] β=2.5, θ=π:  LONG a0*=−2.142 G″(a0*)=+0.6942 G″(0)=+0.1360 | TRANS a0*=0 G″=+0.3623
[E2] LONG barrier: A−Bβ² = +0.135995 (β=2.5), 0.000 (β=4.30661), −0.071 (β=5); minimiser stays >0
[E3] L=1 {0.694,0.362}  L=2 {1.012,0.867}  L=3 {0.646,0.427}   (LONG, TRANS; no 1/L² decay)
[E4] free Gaussian G″=κ→0 ;  β=0 Haar G″=0.2051 > N/12=0.1667
```
No third-party sources downloaded.

---

## 7 · Sharpest line

**The constant (softest) mode's *exact* effective potential is the finite-dimensional
Haar+Wilson marginal, and on any finite small torus its curvature at the coarse minimiser is
unconditionally non-negative `[T]` — `+0.694` (Cartan) and `+0.362` (off-diagonal Savvidy class) at
`β=2.5` on the minimal torus, `O(0.4–1.0)` across `L=1,2,3`, robust under twist, no `1/L²` softening
`[V]`. The one-loop determinant of non-constant modes is the GPY free-gluon `c/L²` piece: nonzero for
the free theory, center-channel, vanishing in decompactification — it is *not* the floor `[NG]`. So the
constant mode carries a genuine positive margin, the strong "uniform-in-`a₀`" form is false only at the
barrier for `β>4.31` `[NG]`, and the residue is exactly the uniform-in-`L` bound `δ_W` together with the
transport to the non-constant soft mode and `ℝ⁴` — i.e. `H1`, unchanged `[O]`. Clay 0%. No QED.**

---

*Refs.* `H4-HEAD-COUNTEREXAMPLE.md` §1 (interior minimiser, `G″(a₀*)=−ρ''/ρ≥0` `[T]`), §3 (Cartan
`L`-scan); `H2-DEFECT-BOUND.md` (`G″=S−D`, measure-carried positivity, Fisher identity);
`V3-COUNTEREXAMPLE-SCOPE.md` / `STAR2-CLASSIFICATION.md` §1 (`A−Bβ²`, `β*=4.30661`, `[NG]` strong form);
`k0-CONSTRAINT-EFFPOT.md` (`V_eff″=G″` identity, `W_geo≥0`, β=0 `N/12`); `k0-HOLONOMY-EFFPOT.md`
(one-loop `V″=c/L²`, free-gluon, center channel, `[NG]`); `k0-EXACT-ONEMODE.md` (`G″=S−D`, reduction
drops `D`); `K0-BLIND.md` / `K0-CONSTANT-MODE-PROTECTION.md` (constant mode = tachyon's home, `Var`
bound = H1); `docs/paper-v2/ch04-gap.tex` (open cell = H1). External (standard, cited from corpus, not
re-fetched): Prékopa 1971; Brascamp–Lieb 1976; Savvidy PLB **71** (1977) 133; Gross–Pisarski–Yaffe RMP
**53** (1981) 43. Witness `research/v3/w7/sm4_constant_mode.py`. *No third-party sources downloaded.*
