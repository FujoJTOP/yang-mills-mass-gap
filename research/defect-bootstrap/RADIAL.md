# 2-variable radial block `R_{(Δ,j)}(u,v)`

> CPB-attack probe (round 4). Verdict: the full radial closed form is **not obtained** (Heun-type);
> the crux is whether it is **1- or 2-variable**, set by the `K₃`-invariant cross-ratios.

## PDE

With `𝔽_j = Σ_j R_{(Δ,j)}(u,v)·C_j^{(1/2)}(cos θ)`, the `SO(3)` Casimir contributes `j(j+1)` and
the radial block is fixed by the `SO(2,1)` Casimir on the two radial invariants:

```
𝕃²_{SO(2,1)}(u,v) R_{(Δ,j)} = Δ_j(Δ_j−1) R_{(Δ,j)} .
```

The radial pair: the two `SO(2,1)`-invariants of `(z₁,ρ₁,z₂,ρ₂)`. A `K₃`-invariant normalisation
is **not** the parent's `(ξ, ρ₁/ρ₂)`: `ρ₁/ρ₂` is dilation- but not `K₃`-invariant. The exact
`K₃`-invariant pair must be taken from 2501.06900.

## Solution (tractable limits)

- **OPE (`ξ→0`):** `R ≃ ξ^{Δ_j}`, independent of `v` at leading order — matches `PDE-SOLUTION.md`.
- **Collinear / planar limit** (one radial invariant frozen): the radial equation reduces to a
  boundary-CFT-type hypergeometric ODE; the solution is `ξ^{Δ_j} ₂F₁(·;ξ)`, the same 1D-block
  class computed to `1e-15` (`NUMERICS.md`).
- Representative: `(2,0)`: `ξ²(1+O(ξ))`; `(2,1)`: `ξ²(1+O(ξ))·(angular cos θ handled by
  C_j^{(1/2)})`.

## Free limit

At `g=0` the Wilson line decouples (path ordering trivialises), so the free bulk–defect block is
the free bulk propagator — a **trivial** benchmark; the informative comparison is interacting.

## Honest boundary

The full 2-variable closed form is **not obtained** (coupled Heun-type PDE), and there is an
**open tension**: `PDE-SOLUTION.md` says `r` enters only in the normalisation (→ effectively
1-variable), while the invariant count gives **2** radial variables. Resolving this requires the
explicit `K₃`-invariant cross-ratios. A numerical spectral solve of `𝕃²_{SO(2,1)}` on the correct
`(u,v)` is the next concrete step; not done here.
