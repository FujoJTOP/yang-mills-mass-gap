# CONSTRUCT-WORLDSHEET-CFT — the flux-tube worldsheet is a 2D CFT whose consistency fixes `c=D−2` and the universal Lüscher term, but the rigidity `κ` is an **irrelevant dimension-4 low-energy constant**; the worldsheet bootstrap therefore does **not** pin `κ`, hence not `m_{0⁺⁺}`

> **VERDICT (construction, no QED).** The confining flux tube is a 2D worldsheet CFT of `c=D−2` transverse
> bosons. **2D CFT consistency** — Virasoro (`c`), modular invariance, reflection positivity — **fixes the
> universal spectrum** (the Lüscher finite-size term `−π(D−2)/(6L)` and the level spacing `2π/L`) but
> **cannot fix the rigidity `κ`**, because `κ` is the coefficient of the worldsheet operator `(∂²X)²`, of
> scaling dimension **4 > 2** ⇒ **irrelevant** ⇒ a **free low-energy constant (LEC)**, invisible to the CFT
> data. Since the rigid-ring `0⁺⁺` mass is `m_{0⁺⁺}/√σ = 2π√(2κ)` (`CONSTRUCT-MASTERFIELD`), the worldsheet
> bootstrap pins `m_{0⁺⁺}` **iff** it pins `κ` — which it does **not**. The genuine positive content: the CFT
> *does* prove **pure Nambu–Goto is unstable** (the closed string collapses; no glueball) and *does* fix the
> universal Lüscher term, isolating `κ` as an **irrelevant LEC = `H1` in worldsheet coordinates.**
> **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric witness. Units `σ=1`.

---

## §1 · The worldsheet CFT

A confining flux tube of length `L` is an effective 2D surface `X^μ(τ,σ)`, `μ=0,…,3`, gauge-fixed to the
static/conformal gauge. The transverse fluctuations `X^i` (`i=1,2`, `D=4`) are **2 free bosons**:
```
c = D − 2 = 2   (transverse central charge),
```
and the low-energy action is Nambu–Goto `σ∫√det` plus, in general, the Polyakov–Kleinert **rigidity**
`(κ/2)∫(∂²X)²`. On a spatial circle of circumference `L`, this is a 2D CFT on a cylinder.

**Spectrum `[T]`.** For the closed flux tube,
```
E_n(L) = σL + (2π/L)(n − c/24) + …     (Nambu–Goto + Lüscher)
       = σL − π c /(6L) + 2π n /L + … ,  c = D−2 .
```
The `−π(D−2)/(6L)` piece is the **universal Lüscher term**, fixed by `c=D−2` alone (the Casimir of the
transverse bosons). For `D=4`: `−π/(3L)·… `, coefficient `−π(D−2)/6 = −π/3`.

---

## §2 · What 2D CFT consistency does and does not fix

**Fixes `[T]`.** (a) **Virasoro / central charge**: `c=D−2` is fixed by the number of transverse directions.
(b) **Modular invariance**: the torus partition function `Z(τ)=∏_{i=1}^{D−2}|η(τ)|^{−2}·(q q̄)^{−c/24}…` is
modular-invariant, fixing the vacuum energy `−π c/(6L)` (the Lüscher term). (c) **Reflection positivity**:
unitarity ⇒ the spectrum is real, the tower `n=0,1,2,…` is gapped above the ground state with spacing
`2π/L`, and the `0⁺⁺` (scalar, zero transverse momentum) is the **closed-string ground state** or its radial
("breathing") excitation.

**Does NOT fix `[T]`/`[NG]`.** The rigidity `κ` multiplies `(∂²X)²`. In worldsheet counting `X` has
dimension `0`, `∂X` dimension `1`, so `(∂²X)²` has **dimension `2×2 = 4 > 2`** — an **irrelevant operator**.
Its coupling `κ` is a **low-energy constant**: CFT data (Virasoro, modular, RP) constrain only *marginal*
(dim 2) and *relevant* (dim < 2) couplings; an irrelevant coupling is **not** fixed by the CFT. This is the
standard EFT statement — `κ` is in the same class as every Wilson coefficient of a higher-derivative term.

---

## §3 · Pure Nambu–Goto is unstable ⇒ no glueball (matches the corpus)

The closed NG string energy `E(L) = σL − π c/(6L)` has derivative `E'(L) = σ + π c/(6L²) > 0` for all
`L>0` ⇒ **monotonically increasing** ⇒ minimum at `L→0`, `E→−∞`: the **shrink modulus** `[T]/[V-num]`.
This reproduces `PROOF-EFFECTIVE-STRING` / `TRANSMUTE-STRING-SPECTRUM`: **pure NG has no stable closed
string, hence no neutral `0⁺⁺`.** Stabilisation needs a **positive `1/L` coefficient**, i.e. a rigidity that
overcomes the negative Lüscher Casimir — `κ > 0` (with the sign of the bending term chosen so `+πκ/ρ` beats
`−πc/(6ρ)`; the pure-NG point `κ=0` is the unstable boundary).

---

## §4 · The `0⁺⁺` as a worldsheet scalar, and its mass

With the rigidity, `CONSTRUCT-MASTERFIELD`'s neutral circular loop has
```
E(ρ) = 2πσρ + πκ/ρ ,  interior minimum ρ* = √(κ/2σ) ,
m_{0⁺⁺} = E(ρ*) = 2π√(2σκ)  ⇒  m_{0⁺⁺}/√σ = 2π√(2κ) .
```
The `0⁺⁺` is the **radial (breathing) mode** of the worldsheet — a worldsheet scalar (the loop-radius
fluctuation), i.e. the lowest `J^{PC}=0^{++}` closed-string scalar. `[C]` (model identification).

---

## §5 · Verdict: the bootstrap pins the universal part, not `κ`

- **What is pinned `[T]`:** `c=D−2`; the Lüscher term `−π(D−2)/(6L)`; the gapped tower `2πn/L`; pure NG
  instability (shrink modulus). These are genuine 2D-CFT facts and pass the acid test (`m→0` as the string
  disappears).
- **What is NOT pinned `[NG]`:** the rigidity `κ` — a dimension-4 **irrelevant LEC**, outside the CFT data.
  Hence `m_{0⁺⁺}/√σ = 2π√(2κ)` is **not** fixed by the worldsheet bootstrap. **`κ` = `H1` in worldsheet
  coordinates** (the flux-tube rigidity / core thickness — the same object `TRANSMUTE-STRING-SPECTRUM` and
  `CONSTRUCT-MASTERFIELD` named).
- **The one loophole, and why it does not help `[C]/[O]`:** if the flux tube were the **exact** NG string,
  the `(∂²X)²` term is *generated* by the NG determinant and `κ` would be calculable — but that is the
  **unstable** `κ=0` boundary (no glueball). A physical `κ>0` (the rigid string) is an **extra** LEC; the
  CFT does not distinguish "generated" from "free", and that distinction is non-perturbative = `H1`.

---

## §6 · Witness (`w7/construct_worldsheet_cft.py`)

Exact arithmetic: confirms `c=2`, the Lüscher coefficient `−π/3`; the NG monotone-collapse; the rigid-ring
`m/√σ=2π√(2κ)` with `κ∈[0.138,0.193]` ⇔ `m/√σ∈[3.3,3.9]`; and that `m/√σ` varies with `κ` at fixed `c` (so
`κ` is not fixed by `c`). No QED; Clay 0%.
