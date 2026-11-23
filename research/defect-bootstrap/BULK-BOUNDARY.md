# `ξ→1` bulk-channel boundary condition

> CPB-attack probe (round 3). Verdict: the `ξ→1` data list is fixed; with `ξ→0` it makes the PDE
> a **two-point boundary-value problem**; the explicit boundary function is not closed-form.

## Convention

Take `ξ∈[0,1]`: `ξ→0` the defect-OPE (`ρ_i→0`), `ξ→1` the bulk-OPE (`x₁→x₂`), with
`(1−ξ) ∝ |x₁−x₂|²` near the corner. (The parent's `ξ=ρ₁ρ₂/ẑ²` is the reciprocal; the physics is
convention-independent: bulk-OPE = `x₁→x₂`.)

## Bulk-OPE data fixing `ξ→1`

The operators `O_b` in the `φφ` bulk OPE, with dimensions/spins and coefficients `C_b` weighted
by the defect-background one-point functions `a_b`:

- identity (`Δ=0`);
- `T_{μν}` (`Δ=4`, `ℓ=2`);
- **`tr F²` (`Δ=4`, `ℓ=0`)** — lightest scalar, the `0^{++}` creator; its `Δ` **is** the gap;
- `tr F⁴` (`Δ=8`, `ℓ=0`); double-traces `[φφ]_{n,ℓ}` (`Δ=2Δ_φ+2n+ℓ`).

The block must reproduce, as `ξ→1`,

```
𝔽_j(ξ,θ) ~ Σ_b (data) × (1−ξ)^{(Δ_b−2Δ_φ)/2+…} × C_{j_b}^{1/2}(cos θ) ,
```

i.e. bulk-channel exchange of each `O_b`; the leading scalar power is set by `tr F²` (`Δ=4`), the
tensor power by `T_{μν}` (`Δ=4`, `ℓ=2`).

## How it enters

With the defect-OPE data at `ξ→0` (`𝔽_j ~ ξ^{Δ_j} P_j(cos θ)`), this `ξ→1` matching makes the
PDE a **two-point boundary-value problem**.

## Honest boundary

Operator content is solid; exact exponents/normalisations depend on the cross-ratio convention
and the one-point functions `a_b`, and the bulk-channel block is itself not closed-form. This
fixes the structure and data list, not the explicit boundary function.

Sources: Billò–Gonçalves–Lauria–Meineri 1601.02883 (bulk–defect OPE); arXiv:2501.06900;
Dolan–Osborn hep-th/0011040 (bulk block / Casimir).
