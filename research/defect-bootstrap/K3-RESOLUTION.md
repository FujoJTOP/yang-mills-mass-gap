# The block is a product: one radial variable (K₃ resolution)

> CPB-attack probe (the crux). Verdict: the **3-variable object does not exist**; the block is
> `f_{Δ_j}(u)·P_j(cos θ)`, **one** radial variable `u` (H² cross-ratio), a **single ODE**.

## The action

The line defect along `z` preserves `{P_z, D, K_z} ≅ SO(2,1) × SO(3)`. Critically, `K_z` acts
on the transverse radius `ρ`, not just `z`: to `O(b)`, `ρ → ρ−2bzρ`, `z → z−b(z²−ρ²)`. This is
the Möbius action `w → w/(1+bw)` on `w = z+iρ`, i.e. `SO(2,1) = PSL(2,R)` = the **H² isometry
group** on the `(z,ρ)` half-plane.

## (a) Two `K₃`-invariant cross-ratios, not three

- radial: `u = |w₁−w₂|²/|w₁−w̄₂|² = [(z₁−z₂)²+(ρ₁−ρ₂)²]/[(z₁−z₂)²+(ρ₁+ρ₂)²] ∈ [0,1]` (the H²
  cross-ratio; `= (cosh d − 1)/(cosh d + 1)`);
- angular: `cos θ = ŷ₁·ŷ₂`.

The parent's `ξ = ρ₁ρ₂/ẑ²` is **not** `K₃`-invariant — explicit `O(b)` variation
`2b(ρ₂²−ρ₁²)/(z₁−z₂)² ≠ 0`. The earlier "3-count" (`KINEMATICS.md`) treated `SO(2,1)` as acting
on the `z`-pair only (orbit 2); with the correct `ρ`-action the `SO(2,1)` orbit on
`(z₁,ρ₁,z₂,ρ₂)` is **3**-dimensional, leaving **one** radial invariant. `8 − 6 = 2`.

## (b) Product

`𝕃²_{SO(2,1)}` acts on the single variable `u`; `𝕃²_{SO(3)}` forces `P_j(cos θ)`. The block is a
**product** `f_{Δ_j}(u)·P_j(cos θ)`.

## (c) Resolution

`PDE-SOLUTION.md` is right; `RADIAL.md` is wrong. "`r` only in the normalisation" means **one**
radial invariant; `RADIAL.md`'s "2 radial variables" came from the non-invariant `(ξ, ρ₁/ρ₂)`.
In `K₃`-invariant coordinates there is **one** radial variable. **The 3-variable object does not
exist.**

## Honest boundary

The invariants are derived exactly (infinitesimal `K_z` check). The radial ODE in `u` for
`f_{Δ_j}` — hypergeometric-type, a **single** ODE — is the object to solve next, not done here.
