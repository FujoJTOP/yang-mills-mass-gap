# CONSTRUCT-WIDTH-KAPPA-LATTICE — the flux-tube width fixes `α·r₀²`, not the rigidity: `κ` is NOT extractable from width data (UV/core degeneracy = the same UV/IR separation), and PK `α` ≠ MASTERFIELD `κ`

> **VERDICT (no QED).** The confining-string width `w²(R)=(D−2)/(2πσ)ln(R/r₀)+w₀²` has a **universal (rigidity-free) logarithmic coefficient** `[T-mod. EFT]` and a constant term `w₀²` that, in the Polyakov–Kleinert rigid string, equals `(D−2)/(4πσ)·ln(ασr₀²)`. This is **invertible for the rigidity `α` only up to the UV core scale `r₀`**: width data fixes the combination **`α·r₀²`**, not `α` (`[O]`, a UV/IR-separation degeneracy of the same type as `H1`). Worse, the MASTERFIELD object `κ` (a **dimensionless** boundary-curvature stiffness, `m_{0⁺⁺}/√σ=2π√(2κ)`) is a **different parameter** from the **dimensionful** worldsheet rigidity `α` — there is **no clean map**. So `κ` is **not extractable** from width data; it stays a **fit** to the spectrum (`κ∈[0.138,0.193]`), not a derivation. **Clay 0%.**

---

## §0 · Labels
**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric · **[T-mod. EFT]** true modulo the effective-string (Polyakov–Kleinert) assumption. Acid: `≤0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The width and its universality

For a confining flux tube of length `R` in `D=4`, the mean-square transverse width is (Lüscher 1981; lattice: Caselle–Necco; Bali–Schlichter):

    w²(R) = (D−2)/(2πσ) · ln(R/r₀) + w₀²           (per-direction: 1/(2πσ))

The **logarithmic coefficient is universal** and **rigidity-free**: it is fixed by the Nambu–Goto measure alone (it does not renormalize under the rigidity term). So the *leading* width gives **no** information about `κ`/`α`.

---

## §2 · The Polyakov–Kleinert constant term

PK action `S = σ∫√g + (1/2α)∫(ΔX)²`. The transverse propagator is `1/(σk² + k⁴/α)`, and

    w²(R) = (D−2)/(2π)∫_{1/R}^{1/r₀} dk /(k(σ + k²/α)).

Evaluating (UV cutoff `1/r₀`): the `ln R` part is universal, and the **constant** is

    w₀² = (D−2)/(4πσ) · ln(α σ r₀²) (+ UV const independent of α).

This is the clean constructive relation `w₀² = w₀²(α, σ, r₀)`, **invertible for `α` given `r₀`:**

    α = (1/(σ r₀²)) · exp( 4πσ w₀²/(D−2) ).

---

## §3 · The `r₀` degeneracy `[O]`

`α` enters only through the product **`α·r₀²`** (`r₀` = UV/core cutoff). The width measures `α·r₀²`; the core scale `r₀` is **not independently known**, so `α` is **not determined**:

    α ∝ 1/r₀²      (fixed `w₀²`, σ).

This is the **same UV/IR-separation obstruction** that `PROMISING-SYNTHESIS` proved is `H1`: you cannot separate the IR rigidity scale from the UV core thickness using width data.

---

## §4 · Numbers (`[V-num]`, σ = (440 MeV)² ⇒ √σ = 2.23 fm⁻¹, σ = 4.97 fm⁻²)

| quantity | value |
|---|---|
| Lüscher per-direction log coeff `1/(2πσ)` | 0.0320 fm² |
| summed `(D−2)/(2πσ)` | 0.0640 fm² |
| const-term coeff `(D−2)/(4πσ)` | 0.0320 fm² |
| `w₀²=0.15 fm²` ⇒ `α·r₀²` | 21.8 fm² |
| then `α` at `r₀=0.15 / 0.20 / 0.30 fm` | **970 / 545 / 242** (∝1/r₀²) |

Representative `w₀²` values (Bali–Schlichter/Caselle–Necco order; **numbers illustrative, verify against the source table**). The spread 970→242 is **pure `r₀`-ambiguity**, not a measurement.

---

## §5 · MASTERFIELD `κ` is a fit, and a different object

`m_{0⁺⁺}/√σ = 2π√(2κ) ⇒ κ = (m/√σ)²/(8π²)`:

| `m/√σ` | `κ` |
|---|---|
| 3.3 | 0.138 |
| 3.7 | 0.173 |
| 3.9 | 0.193 |

So `κ∈[0.138,0.193]` is recovered **by fitting the spectrum** — exactly the input MASTERFIELD left open. Crucially **PK `α` ≠ MASTERFIELD `κ`**: `α` is dimensionful (fm²) and `r₀`-degenerate; `κ` is dimensionless. No established dictionary connects them; asserting `κ=α` would be a model choice, not a derivation.

---

## §6 · Honest bottom line

- **Width data does NOT determine κ.** The log is universal (rigidity-free); the constant fixes `α·r₀²`, degenerate with the UV core `r₀` (spread ~4× over physical `r₀`).
- **The degeneracy is `H1` in width coordinates**: separating the IR rigidity from the UV core is the UV/IR separation.
- **Type mismatch**: PK worldsheet rigidity `α` (dimensionful) vs MASTERFIELD boundary stiffness `κ` (dimensionless) — no clean map.
- **What the construction yields**: a clean `w₀² = (D−2)/(4πσ)ln(ασr₀²)` relation (useful diagnostic) and the numerical band `κ∈[0.138,0.193]` as a *fit*.
- **κ remains underived** = `H1` in string coordinates. **Clay 0%; no QED.**

*Witness:* `research/unconditional/w7/construct_width_kappa_lattice.py` (numpy, runs clean). *Refs (verify):* Lüscher, Nucl. Phys. B180 (1981) 317; Caselle–Necco, flux-tube width; Bali–Schlichter, flux-tube profile.
