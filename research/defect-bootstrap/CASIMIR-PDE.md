# 3 cross-ratios + Casimir PDE system

> CPB-attack probe (sub-task 4 of round 3). The well-posed PDE promised in
> `BLOCK-DERIVATION.md` is now explicit; the last piece is to *solve* it.

## Setup

`R⁴ = (z∈R along line) × (y∈R³ transverse)`, `ρ=|y|`; bulk scalars at `x_i=(z_i,y_i)`, `i=1,2`.

## Count

`8` coordinates − generic orbit `5` (`SO(2,1)` on the `z`-pair: orbit `2`, stabiliser `1`;
`SO(3)` on the `y`-pair: `3`) = **3 invariants**.

## Explicit cross-ratios (candidate)

With `ẑ=z₁−z₂`, `ρ_i=|y_i|`, `c=(y₁·y₂)/(ρ₁ρ₂)`:

```
r = ρ₁ρ₂/ẑ² ,    s = ρ₁/ρ₂ ,    c = cos θ .
```

Dimensionless and invariant under `P₃, D, SO(3)`. **Caveat:** `r` is **not manifestly
`K₃`-invariant** (special conformal along `z` scales `ρ` and `ẑ` by different powers). The
`K₃`-invariant normalisation is the inversion-built one; take the exact convention from
2501.06900.

## Casimir PDE system

Exchange a defect primary `O_j` of dimension `Δ_j`, `SO(3)` spin `j`. Two quadratic Casimirs
act:

```
D_{SO(2,1)} G = Δ_j(Δ_j−1) G ,     D_{SO(3)} G = j(j+1) G ,
```

with `D_{SO(2,1)}` the sum-of-squares of the conformal Killing vectors `{P₃,D,K₃}` (in `z,ρ`)
and `D_{SO(3)}` the transverse rotation Casimir (in `ρ_i, c`). The block is the **simultaneous
eigenfunction in the 3 cross-ratios**; OPE limits fix the boundary conditions (`r→0`:
defect-OPE power; `ξ→1`: bulk channel).

## Derived vs conjectured

**Derived:** the count (3), the eigenvalues `Δ_j(Δ_j−1)` and `j(j+1)`, the PDE structure.
**To-verify:** the explicit `K₃`-invariant normalisation of `(r,s,c)` and the resulting
coefficient functions — match to 2501.06900 App. B (which gives only the *defect-4-point*
analogue `(t/(1+t))^Δ ₂F₁`). The 3-variable bulk–defect block itself is not in that paper.
