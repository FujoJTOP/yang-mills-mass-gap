# CONSTRUCT-LUSCHER — effective-string prediction: the 0⁺⁺ as the closed-string worldsheet scalar, the finite-V law is `σL − π/(3L)`, but the absolute mass needs the rigidity `κ` (= `H1`)

> **VERDICT (no QED).** I construct the Lüscher / Nambu–Goto effective-string prediction for the
> 0⁺⁺ glueball. Three results. **(i)** For a closed flux tube (the glueball model), the spectrum is
> `E_N(ℓ) = σℓ + (2π/ℓ)(N+N̄ − (D−2)/12)`, level-matched `N=N̄`: **`E₀(ℓ) = σℓ − π/(3ℓ)`** (ground state)
> and `E₁(ℓ) = σℓ + 11π/(3ℓ)`, gap `4π/ℓ`. The **closed Lüscher term is `−π(D−2)/(6ℓ) = −π/(3ℓ)`, four
> times the famous open-string coefficient `−π(D−2)/(24ℓ) = −π/12ℓ`**. **(ii)** The 0⁺⁺ **is** the
> closed-string ground state — the worldsheet vacuum (no transverse excitation) is a scalar, spin 0, even
> parity/charge. **(iii)** The finite-volume scaling law is the torelon formula
> **`M(L) = σL − π/(3L) + O(1/L³)`** (approach from below, relative correction `−π/(3σL²)`); for the
> *stable* glueball the shift is Lüscher-exponential, `Δm(L) ∝ e^{−√3 mL}/(mL)` (two-glueball threshold,
> no pions in pure YM). **But pure NG has no interior minimum** (`E₀′(ℓ)=σ+π/(3ℓ²)>0` ⇒ the closed string
> shrinks), so the **absolute** mass `m_{0⁺⁺}/√σ≈3.9` requires the rigidity `κ` — exactly `MASTERFIELD`'s
> `m/√σ = 2π√(2κ)`, fitted `κ∈[0.138,0.193]`. **The effective string fixes the scaling, not the mass.** `κ=H1`.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · Setup: the effective string in D=4

A confining flux tube of length `ℓ` is modelled by a Nambu–Goto string in `D=4` spacetime; the transverse
fluctuations are `D−2 = 2` free massless 2D bosons `X^i(τ,σ)`, `i=1,2`. The glueball is modelled as a
**closed** flux tube of circumference `ℓ` (the Isgur–Paton / Lund "glueball as a closed string" picture,
`[C]` — a model, not a theorem).

**Closed-string spectrum** (light-cone, level matching `N = N̄`):

```
E(N,N̄; ℓ) = σℓ + (2π/ℓ)(N + N̄ − (D−2)/12),   N = N̄ = 0,1,2,…
```

- **Ground state** `N=N̄=0`:  `E₀(ℓ) = σℓ − π(D−2)/(6ℓ) = σℓ − π/(3ℓ)`  in `D=4`.
- **First excited** `N=N̄=1`:  `E₁(ℓ) = σℓ + 11π/(3ℓ)`;  gap `E₁−E₀ = 4π/ℓ`.

(The zero-point sum `Σ_{m≠0}(2π|m|/ℓ)` regularized gives `−(D−2)·(2π/ℓ)·(1/12) = −π(D−2)/(6ℓ)`, hence the
closed-string Casimir.)

---

## §2 · (ii) The 0⁺⁺ is the worldsheet scalar  `[T]`

Under the transverse rotation group `SO(D−2)=SO(2)`:

- the **ground state** (no oscillators) is invariant ⇒ **spin 0**, parity `+`, charge-conjugation `+` ⇒
  **`J^{PC}=0^{++}`**. So the *lightest* closed-string state is a `0⁺⁺` scalar.
- the **first excited** level `α^i_{−1}ᾱ^j_{−1}|0⟩` decomposes into symmetric-traceless (`2⁺⁺`) **and** the
  trace `δ_{ij}α^i_{−1}ᾱ^j_{−1}|0⟩` (`0⁺⁺`). So a second `0⁺⁺` sits at the first excited level.

⇒ In the effective-string picture the `0⁺⁺` **is** the worldsheet scalar (ground state), consistent with
the flux-tube expectation that the lightest glueball is the "dilaton"/breathing mode.

---

## §3 · (i) The pure-NG obstruction: no interior minimum  `[NG]`

`E₀′(ℓ) = σ + π/(3ℓ²) > 0` for **all** `ℓ>0`. The closed NG string energy is **strictly increasing** with
circumference, and `E₀ → −∞` as `ℓ→0`: the string **shrinks**. There is **no stable finite-size closed
string**, hence no glueball mass — the same obstruction `CONSTRUCT-MASTERFIELD` found (`TRANSMUTE-STRING-SPECTRUM`
killed pure NG earlier). The absolute `0⁺⁺` mass cannot come from pure NG. `[NG]` for pure NG as the mass
source.

---

## §4 · (i′) The rigid completion gives the absolute mass  `[C]`

Adding the rigidity (extrinsic-curvature) term `πκ/ℓ` (the `∮K²ds` term, `MASTERFIELD`):

```
E(ℓ) = 2πσℓ + πκ/ℓ,   ℓ* = √(κ/2σ),   m = E(ℓ*) = 2π√(2σκ),   m/√σ = 2π√(2κ).
```

This has an **interior minimum** (unlike pure NG). Fitting the lattice band `m_{0⁺⁺}/√σ ∈ [3.3, 3.9]`:
`κ ∈ [0.138, 0.193]` (`[V-num]`). **The rigidity `κ` is not predicted by the effective string** — it is the
flux-tube core/rigidity, exactly `H1` in string coordinates (`[O]`).

The Lüscher term alone gives no absolute mass; the rigidity does, at the price of the one input `κ`.

---

## §5 · (iii) The finite-volume scaling law (concrete, checkable)

**Torelon** (a flux tube wrapping an `L`-torus) — the clean effective-string finite-V datum `[T]`:

```
M_tor(L) = σL − π(D−2)/(6L) = σL − π/(3L) + O(1/L³),   relative correction −π/(3σL²).
```

Concrete values (σ=1): `L=2 → M=1.476 (rel −0.262)`, `L=4 → 3.738 (−0.065)`, `L=10 → 9.895 (−0.010)` —
the torelon approaches `σL` **from below** with the `−π/(3L)` Lüscher correction. This is a sharp,
finite-volume-checkable prediction.

**Stable glueball** on a finite volume `L³`: for a particle of mass `m` below its two-body threshold with no
light hadrons (pure YM has no pions), Lüscher's finite-volume formula gives an **exponentially small** shift
`[T]`:

```
Δm(L) ∝ − (g_{0⁺⁺0⁺⁺0⁺⁺}²) · e^{−√3 mL} / (mL)   (threshold √(4m²−m²) = √3 m),
```

where `g` is the `0⁺⁺–0⁺⁺–0⁺⁺` coupling (an external input). The scaling `Δm ∝ e^{−√3 mL}/(mL)` is
structure-independent; only the coefficient needs `g` `[O]`.

---

## §6 · Comparison to the lattice and to the Lüscher coefficient

- **Lüscher coefficient `[T]`:** the famous **`−π/12`** is the **open-string** `D=4` term
  `−π(D−2)/(24R) = −π/(12R)` (static potential). The **closed-string** (glueball) analog is
  `−π(D−2)/(6ℓ) = −π/(3ℓ)` — **exactly 4× the open coefficient** (periodic vs Dirichlet zero-point sum).
  Both are universal (regulator/action-independent for the NG class).
- **Lattice `m_{0⁺⁺}/√σ ≈ 3.9`:** the effective string alone **does not** reach 3.9 (pure NG has no
  minimum). With the rigid completion it does, **given** `κ=0.193`. So the lattice value **selects** `κ`
  but does not **derive** it.
- **Torelons at large `L`:** the `−π/(3L)` correction is the sharp check; agreement with lattice torelon
  data would confirm the closed-string Casimir `−π/3`.

---

## §7 · Honest bottom line

- **What the effective string DOES give `[T]`:** the closed-string spectrum (`E₀=σℓ−π/(3ℓ)`), the `0⁺⁺`
  = worldsheet scalar identification, the `−π/3` (closed) Lüscher term (4× the open `−π/12`), and the
  finite-volume scaling laws (torelon `σL−π/(3L)`; glueball Lüscher-exponential `~e^{−√3 mL}`).
- **What it does NOT give `[O]`:** the **absolute** `m_{0⁺⁺}` — pure NG has no interior minimum, so the
  mass needs the rigidity `κ` (fitted `0.138–0.193`), which the effective string cannot predict. `κ=H1`.
- **The division in one line:** effective string = the **scaling and the quantum numbers**; `κ` (rigidity,
  core thickness) = the **absolute mass** = the one missing input = `H1`.

**Clay 0%; no QED.**
