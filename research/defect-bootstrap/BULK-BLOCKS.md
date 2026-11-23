# 4D bulk conformal blocks (Dolan–Osborn)

> CPB-attack probe (sub-task 4). Verdict: the **bulk half is known**; the defect half is the
> missing ingredient.

## Cross-ratios

`u = z z̄`, `v = (1−z)(1−z̄)`.

## 4D block (Dolan–Osborn)

For external scalars of equal dimension `Δ_φ`, exchanged primary of dimension `Δ`, spin `ℓ`:

```
g_{Δ,ℓ}(z,z̄) = (z z̄/(z−z̄)) · [ k_{Δ+ℓ}(z) k_{Δ−ℓ−2}(z̄) − k_{Δ−ℓ−2}(z) k_{Δ+ℓ}(z̄) ] ,
k_β(z) = z^{β/2} ₂F₁(β/2, β/2; β; z) .
```

For unequal external dimensions `Δ_i`, replace the `₂F₁` arguments by `(β+Δ₂−Δ₁)/2` and
`(β+Δ₃−Δ₄)/2`. (Normalisation: the conformal block up to the standard OPE prefactor; the
"conformal partial wave" differs only by a `Δ`-dependent normalisation.)

## Defect crossing conventions

Bulk-OPE channel `⟨φφ⟩_def = Σ_O C^{φφ}_O 𝒢_O(ξ)`, with `𝒢_O` the above block at the defect
cross-ratio `ξ`; the defect-OPE channel `Σ_i a_i² 𝔽_i(ξ)` is the (missing) other side.

## Free benchmark

Free `g=0` bulk scalar decomposes into **double-trace blocks** (`Δ = 2Δ_φ + 2n + ℓ`, `ℓ`
even) — the exact bulk-side solution to beat.

## Honest boundary

Bulk blocks are known (Dolan–Osborn hep-th/0011040, hep-th/0309180; Simmons-Duffin TASI,
arXiv:1602.07982, eq. 3.16). The **defect** blocks are the missing ingredient; this covers
only the known bulk half.
