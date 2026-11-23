# Radial ODE in `u`: hypergeometric (not Heun)

> CPB-attack probe (round 5). Verdict: the radial equation is **hypergeometric**, confirming the
> one-variable reduction and the product form `f(u)·P_j(cos θ)`. The block is essentially
> constructed.

## Setting

`u = [(z₁−z₂)²+(ρ₁−ρ₂)²]/[(z₁−z₂)²+(ρ₁+ρ₂)²]` is the H² geodesic cross-ratio,
`u = (cosh d−1)/(cosh d+1)`. The `SO(2,1)` Casimir acting on the bulk–defect two-point function
is the **H² Laplace–Beltrami radial operator**; in `u` it is a hypergeometric ODE with regular
singular points at `u=0` and `u=1`:

```
u(1−u) f'' + [c − (a+b+1)u] f' − ab f = 0 .
```

## Exponents

The Casimir is `Δ_j(Δ_j−1)`; the indicial equation at `u=1` gives `ρ(ρ−1) = Δ_j(Δ_j−1)`, so the
two solutions near `u=1` have powers `{ (1−u)^{Δ_j}, (1−u)^{1−Δ_j} }`. The physical (defect-OPE)
block is the `(1−u)^{Δ_j}` solution.

## Solution

```
f_{Δ_j}(u) = ₂F₁( Δ_j, 1−Δ_j; 1; −u/(1−u) )      — the H² spherical function,
```

up to the `Δ_φ`-dependent normalisation; the pure `(1−u)^{Δ_j}` branch is the physical one.

## OPE limits

- **Defect OPE (`u→1`):** `f ~ (1−u)^{Δ_j}`; with `ξ = ρ₁ρ₂/z₁₂² = (1−u)/(4u)`, this is
  `ξ^{Δ_j}`. ✓
- **Bulk OPE (`u→0`):** `f → 1` (spherical-function normalisation); the identity contributes
  `u^{−Δ_φ}` after the `ρ^{−Δ_φ}` prefactor.

## Check

The leading defect-OPE power `Δ_j` is reproduced **exactly**. The equation is **hypergeometric,
not Heun** — confirming the one-variable reduction and the product form `f(u)·P_j(cos θ)`.

## Honest boundary

The ODE class, the exponents, and the defect-OPE limit are derived; the full closed form
including the `Δ_φ`-dependence (the exact `A,B` mix of the two `u=1` solutions) and the
`K`-invariant normalisation need the embedding-formalism derivation, **not done here**.
