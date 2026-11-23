# Correct codim-3 free benchmark

> CPB-attack probe (round 5). Verdict: the free block is the **Legendre function of the second
> kind** `Q_j(cosh η)` — the exact anchor for the `u`-ODE.

## Setup

Free scalar (`Δ_φ=1`), line defect **codim-3** (`SO(2,1)×SO(3)`, transverse `R³`,
`cos θ = ŷ₁·ŷ₂`). **Not** the codim-2 monodromy.

## Cross-ratios

`Z=z₁−z₂`, `d²=Z²+(ρ₁−ρ₂)²`, `u = tanh²(η/2)` (H² cross-ratio), `cosh η = (1+u)/(1−u)`.

## Free 2-point

```
|x₁−x₂|² = (d²/(2u))·[(1+u) − (1−u)cos θ] ,
G_free = 1/(4π²|x₁−x₂|²)
       = u/(2π² d² (1−u)) · Σ_{l≥0} (2l+1) Q_l(cosh η) P_l(cos θ) .
```

(Heine expansion `1/(z−x) = Σ (2l+1) Q_l(z) P_l(x)`; verified at `θ=0` ⟹ `1/(4π²d²)` exactly;
`l=0` gives `Q₀ = −½ ln u`, matching the 1D `Δ=1` log.)

## Free block

```
f_l(u) ∝ (2l+1) Q_l(cosh η)
      = N_l·((1−u)/(1+u))^{l+1}·₂F₁((l+1)/2, (l+2)/2; l+3/2; ((1−u)/(1+u))²) .
```

So the free block is the **Legendre function of the second kind**, an **H²-hypergeometric** in
`u`. At the free point the parent's `f_{Δ_j}(u)` must reproduce `Q_j(cosh η)` (up to
normalisation); the free OPE ties the radial label to the spin `j` (a single `l` labels both).

## Honest boundary

This is the **trivial / identity** line defect (= the free propagator). A genuinely nontrivial
conformal line defect (a GFF defect) gives a *different* 2-point, not computed here; `u→0` is
`x₁→x₂` (OPE), `u→1` is a point on the defect — **opposite** to the parent's `ξ` convention.
Codim-3 throughout (`P_l(cos θ)`, transverse `R³`); it is **not** the codim-2 monodromy.

Sources: standard Heine Legendre expansion; DLMF 14.18 / `Q_l`; no codim-2 monodromy used.
