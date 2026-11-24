# KAPPA-EFT-MATCH — microscope-matching the rigid-string rigidity κ: no lower-dimension coefficient fixes it; κ = (m/√σ)²/(8π²) = the gap in EFT coordinates

> **VERDICT ([O]/[C-model]; no QED).** The rigidity κ of the rigid-string 0⁺⁺ (model
> `E(R)=2πσR+πκ/R`, `m₀₊₊/√σ=2π√(2κ)`, `κ∈[0.138,0.193]`) is **not pinned by any
> lower-dimension coefficient**. Geometrically κ is the tube's **bending stiffness**,
> `κ = O(1)·σR_core²` — a pure **dimensional** statement (`R_core ~ c/√σ`, `κ = 2c²`,
> band ⇔ `c∈[0.263,0.311]`), no number. The one place a microscopic input enters — the
> core radius set by the gluon condensate — brings only a **logarithm `ln(σ/Λ²)`** and
> **drifts with `Λ` (= the condensate scale = H1)**. For the pure NG string the rigidity
> is **redundant** (Gauss–Codazzi) and the string is **unstable**, so NG fixes no physical
> κ. In the model the inversion is **exact**: `κ = (m/√σ)²/(8π²)` — **matching κ is
> determining the gap.** `H1` in EFT coordinates, confirmed by a fifth route.

---

## §0 · Labels / acid

**[T]** 真定理 · **[C-model]** conditional on the rigid-string model · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric. Acid: must vanish on free massless and fail on `:φ²:` — the rigid-string `0⁺⁺` is a gapped bound state, and `κ→0` collapses it (no atom), so κ carries the discriminating content (non-vacuous).

---

## §1 · The model and the exact inversion `[T-model]`

Closed flux ring, rigid string, `E(R)=2πσR+πκ/R`. Interior minimum `R*=√(κ/2σ)`, ground
`m₀₊₊=2π√(2σκ)`. Invert:

```
m₀₊₊/√σ = 2π√(2κ)   ⟺   κ = (m₀₊₊/√σ)² / (8π²).
```
`m/√σ = 3.3, 3.5, 3.7, 3.9 ⇒ κ = 0.138, 0.155, 0.173, 0.193` exactly (witness §A). **κ and
the gap are the same information** — this alone makes "match κ" ≡ "compute the gap".

## §2 · Geometric matching: κ = bending stiffness `[C-model]`

For a tube of tension σ and core radius `R_core`, the bending energy per unit length is
`(1/2)B K²` with stiffness `B = O(1)·σR_core²`; the ring bending energy `πB/R` matches the
model's `πκ/R` at `κ = B`. With `R_core = c/√σ`: **`κ = 2c²`** (witness §B). The band
`[0.138,0.193]` ⇔ `c ∈ [0.263,0.311]`, i.e. `R_core ~ 0.3/√σ` — a confinement-length core.
**This is dimensional analysis: O(1) magnitude right, no number.**

## §3 · Pure NG fixes nothing physical `[T]`

The extrinsic curvature squared `K² = (ΔX)²` is **redundant** modulo the leading NG equations
of motion (Gauss–Codazzi identity) — it is not an independent coupling for the pure NG string.
The one-loop NG determinant fixes all coefficients, but at the **unstable** `κ=0` boundary:
`E(R)=2πσR` is monotone increasing ⇒ the closed NG string shrinks (`R→0`, shrink modulus),
there is **no stable closed string and no glueball**. So NG supplies **no physical κ>0**.

## §4 · Condensate / (F²)² enters only logarithmically `[O]`

The tube core radius in the dual/bag picture is `R_core² ~ (1/2πσ)·ln(σ/Λ²)` with `Λ` the
condensate/IR scale. Then `κ = 2σR_core² ~ (1/π)·ln(σ/Λ²)` (witness §C): crude prefactor
**overshoots** (0.36–0.69 for `Λ=0.25–0.15 GeV`), and — decisively — κ **drifts with `Λ`**,
i.e. with the gluon condensate ⟨F²⟩ = `H1`. Any `(F²)²` Wilson coefficient enters only inside
this log; **no lower-dimension coefficient pins κ**.

## §5 · Worldsheet RG: κ is a boundary condition `[T]`

`κ(∂²X)²` is a **dimension-4 (irrelevant)** worldsheet operator; its IR value is a UV
**boundary condition** from the core, not a flow output (WORLDSHEET-CFT). RG cannot
manufacture κ.

---

## §6 · Verdict

Every microscopic route leaves κ undetermined by lower-dimension data:
- pure NG: redundant/unstable `[T]`;
- geometric: `κ=2c²`, dimensional only `[C-model]`;
- condensate: only `ln(σ/Λ²)`, drifts with `H1` `[O]`;
- RG: irrelevant, boundary condition `[T]`.

And the model's own inversion is exact: `κ = (m/√σ)²/(8π²)`. **Matching κ microscopically
is determining the `0⁺⁺` gap** — `κ = H1` in EFT coordinates, a fifth independent route
(joining KAPPA-DERIVE, WIDTH-KAPPA, DUAL-CONDENSATE, WORLDSHEET-CFT).

*Files:* `w7/kappa_eft_match.py` (runs `py -3.12 -I`, clean). Clay 0%; no QED.
