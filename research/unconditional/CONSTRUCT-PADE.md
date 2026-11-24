# CONSTRUCT-PADE — Padé / Borel–Padé continuation of the `0⁺⁺` strong-coupling series: the METHOD works, but the physical input (only the leading term) does **not** pin β=2.5 — estimate `m/√σ ∈ [3.3, 3.9]`, ~15% systematic

> **VERDICT (no QED).** The Padé / Borel–Padé machinery **converges** on Borel-summable asymptotic series
> (calibrated on two models with known exact values), so a *controlled* crossover estimate is possible **in
> principle given enough genuine coefficients**. But the available physical input is only the **leading**
> strong-coupling term `m/√σ = 4√ln(1/u)`; at `u(β=2.5)=0.507` the `O(u)` correction is **not small**, the
> leading continuation **undershoots** the lattice (`m/√σ≈3.9`) by **~15%**, and the unknown correction
> coefficient **spans the whole band** — the continuation is **underdetermined**. Honest estimate:
> `m_{0⁺⁺}/√σ(β=2.5) ∈ [3.3, 3.9]`, best leading point **3.30**, systematic ~15%. **The exact missing input
> is the higher SU(2) strong-coupling coefficients (OS'78 / Münster / SZZ beyond leading).** `[V-num]/[C]/[O]`

---

## §1 · Exact SU(2) input and the leading anchor

`u(β)=I₂(β)/I₁(β)` (fundamental character ratio, exact via `ive`), leading strong-coupling
`m_{0⁺⁺}/√σ = 4√ln(1/u)` (the four-wall tube / area-law ratio, `ROOF-CONFINEMENT-GAP`):

| β | u=I₂/I₁ | ln(1/u) | `4√ln(1/u)` |
|---|---|---|---|
| 0.0208 | 0.00521 | 5.258 | 9.17 |
| 0.5 | 0.1237 | 2.090 | 5.78 |
| 1.0 | 0.2402 | 1.426 | 4.78 |
| 1.625 | 0.3678 | 1.000 | 4.00 |
| 2.2 | 0.4645 | 0.767 | 3.50 |
| **2.5** | **0.5072** | **0.679** | **3.30** |
| 3.0 | 0.5679 | 0.566 | 3.01 |

Lattice `m_{0⁺⁺}/√σ ≈ 3.9` (Morningstar–Peardon, hep-lat/9901004). ⇒ the leading continuation at β=2.5 is
**~15.5% below lattice**. The `O(u)` correction at `u=0.507` is not a small parameter — the leading formula
is the whole problem's difficulty in one number.

---

## §2 · Method calibration (the machine works)

Two Borel-summable asymptotic series with **known exact** values at `u=0.5`:

- **M1** Stieltjes `f₁(u)=∫₀^∞ e^{−t}/(1+ut)dt`, `a_n=(−1)ⁿn!`, exact `f₁(0.5)=0.72265723`;
- **M2** sqrt-branch `f₂(u)=∫₀^∞ e^{−t}/√(1+ut)dt`, `a_n=(−1)ⁿ(2n)!/(4ⁿn!)`, exact `f₂(0.5)=0.84273846`.

**Padé [L/M] of the asymptotic series at `u=0.5`** (error vs exact):

| N=L+M | M1 err | M2 err |
|---|---|---|
| 4 | 6.4e−3 | 2.8e−3 |
| 8 | 3.8e−4 | 1.6e−4 |
| 12 | 3.8e−5 | 1.5e−5 |
| 16 | 5.1e−6 | 2.0e−6 |

**Borel–Padé [L/M]** (`b_n=a_n/n!`): M1 exact from `N=2` (its Borel transform `1/(1+t)` is rational); M2
reaches **1.7e−10** at `N=16`. ⇒ **Borel–Padé is the right engine**: it beats plain Padé by ~4 orders at
`N=16`, and converges monotonically in `N` for a Borel-summable asymptotic series. **A crossover estimate is
deliverable — provided the genuine coefficients are in hand.**

---

## §3 · Physical underdetermination (why β=2.5 is not pinned)

Write `m/√σ(u) = 4√ln(1/u) · g(u)` with `g(0)=1` and `g`'s corrections unknown. The leading data fixes only
`g(0)=1`. A minimal one-parameter family `g_c(u)=1/(1−cu)` matches the leading term exactly and, at `u=0.5`:

| c | g_c(0.5) | m/√σ(2.5) |
|---|---|---|
| 0.00 | 1.000 | **3.30** |
| 0.10 | 1.053 | 3.47 |
| 0.20 | 1.111 | 3.66 |
| 0.30 | 1.176 | **3.88** |
| 0.40 | 1.250 | 4.12 |

The correction coefficient `c` is **unconstrained** by the leading strong-coupling data: **any** `c∈[0,0.5]`
gives `g(0)=1` yet **spans the entire 3.3–3.9 band**. This is the precise sense in which the continuation is
underdetermined — a Padé of a *one-term* series is the one-term series.

---

## §4 · Honest verdict

- **Method `[C]`:** Padé / Borel–Padé converge on Borel-summable asymptotic series (Table §2); a controlled
  crossover estimate is achievable **in principle**.
- **Physical result `[V-num]`:** leading anchor `m/√σ(2.5)=3.30`; lattice `3.9`; leading undershoots by
  **~15.5%**; estimate **`∈ [3.3, 3.9]`**, ~15% systematic. **Not a robust continuation.**
- **Exact missing input `[O]`:** the **higher SU(2) strong-coupling coefficients** of the `0⁺⁺` glueball mass
  and string tension (OS'78 / Münster / SZZ beyond `a_0`). With, say, `N≳8–12` genuine coefficients, the
  calibration shows Borel–Padé would give the crossover value to `≲1%` *if* the series is Borel-summable
  with a single dominant renormalon (not guaranteed; `ROOF-RESURGENCE` shows the IR renormalon has no saddle,
  so the Borel representation is itself the open object). So the numerical continuation is a **useful
  diagnostic**, not a certificate, and its quality is bounded by coefficient availability.

*No QED; Clay 0%.* Witness: `w7/construct_pade.py` (numpy/scipy, single run).
Refs: `PNT-STRONG-BRIDGE-EXTEND.md`, `ROOF-RESURGENCE.md`, `ROOF-CONFINEMENT-GAP.md`.
