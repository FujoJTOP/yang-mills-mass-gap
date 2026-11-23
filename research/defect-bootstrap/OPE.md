# Small-ξ defect-OPE expansion

> CPB-attack probe (sub-task 5). Fixes only the spin-0 power law; the `SO(3)` transverse-spin
> structure is the missing ingredient.

## Correction

A 1D defect in 4D preserves `SO(2,1)×SO(3)` (transverse `R³`), **not** `SO(2,1)×SO(2)` (that
is codim-2). See `KINEMATICS.md`.

## Cross-ratio

`ξ = |y|²/(t²+|y|²)`; `ξ→0` is the defect-OPE limit.

## Expansion

The defect OPE `φ(x) = Σ_O c_O |y|^{Δ_O−Δ_φ} O(t)` gives

```
⟨φφ⟩_def = t^{−2Δ_φ} Σ_O c_O² ξ^{Δ_O−Δ_φ}(1−ξ)^{Δ_φ−Δ_O} .
```

So the spin-0 block is `𝔽_O(ξ) ≃ ξ^{Δ_O−Δ_φ}` as `ξ→0`; the leading term is the identity
(`Δ_O=0`): `ξ^{−Δ_φ}`.

## Consistency

The Casimir ODE's indicial exponents at `ξ=0` must be `{Δ_O−Δ_φ}`; at `ξ=1` the bulk-channel
exponents `{Δ_b−Δ_φ}`. The check: does the Casimir solution reproduce both?

## Honest boundary

Fixes only the spin-0 power law. The full 4D block's `SO(3)` transverse-spin structure — the
subleading, spin-dependent part — is exactly the missing ingredient; this is not a full block.
