# 3-variable Casimir PDE: leading term + the reduced problem

> CPB-attack probe (sub-task 5 of round 3). The "3-variable" block reduces to a **coupled
> 2nd-order PDE in `(ξ,θ)`** with `r` only in the normalisation; leading term in hand.

## Kinematics

Bulk points `(z_i, y_i)`, `y_i∈R³`. Three cross-ratios:

```
ξ = |y₁||y₂|/|z₁−z₂|²  (dilation-invariant; convention-dependent),
cos θ = (y₁·y₂)/(|y₁||y₂|),    r = |y₁|/|y₂| .
```

## Casimir PDE

The quadratic `SO(2,1)×SO(3)` Casimir acting on the bulk 2-point function has eigenvalue
`Δ_j(Δ_j−1) + j(j+1)` on the exchanged defect primary:

```
[ 𝕃²_{SO(2,1)}(ξ) + 𝕃²_{SO(3)}(θ) ] 𝔽_j = [ Δ_j(Δ_j−1) + j(j+1) ] 𝔽_j ,
```

a **coupled 2nd-order PDE in `(ξ, θ)`**; `r` enters only through the normalisation.

## Leading term (defect-OPE, `ξ→0`)

Each bulk field couples as `φ ~ |y|^{Δ_j−Δ_φ} O_j(z)`; the 2-point function gives

```
𝔽_j(ξ, θ, r) ≃ ξ^{Δ_j} P_j(cos θ) · [1 + O(ξ)] ,   independent of r at leading order,
```

with `P_j` the Legendre polynomial. Subleading terms are the 1D conformal `₂F₁` in `ξ`;
`r`-dependence begins at higher order.

## Table (leading `ξ^{Δ_j}·P_j`)

| `(Δ_j, j)` | leading |
|---|---|
| (1, 0) | `ξ` |
| (2, 0) | `ξ²` |
| (2, 1) | `ξ² cos θ` |
| (2, 2) | `ξ² (3cos²θ−1)/2` |
| (3, 0) | `ξ³` |

`P₂(cos θ)` at `θ = 0°, 60°, 90°`: `1, −0.125, −0.5`.

## Honest boundary

This is the **leading** (defect-OPE) term, **not** the full 3-variable block. The full solution
needs the coupled PDE beyond leading order (the `ξ`- and `r`-subleading `₂F₁` structure) —
exactly the missing object. A numerical check is possible via the free-scalar
monodromy-defect 2-point function (known), not run here.
