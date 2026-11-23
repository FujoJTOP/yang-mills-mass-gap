# `SO(3)` transverse-spin structure of the bulk–defect block

> CPB-attack probe (sub-task 6). Verdict: the angular factor is **standard (Gegenbauer)**; the
> missing piece is the **2-variable radial block**.

## `SO(3)` structure

The external bulk scalars are `SO(3)`-neutral, so the only transverse-angle data is
`cos θ = ŷ₁·ŷ₂`. An exchanged defect operator of `SO(3)` spin `j` contributes, by the `SO(3)`
addition theorem, the Gegenbauer factor

```
Σ_m Y_{jm}(Ω₁) Y*_{jm}(Ω₂) = ((2j+1)/4π) · C_j^{(1/2)}(cos θ) ,   C_j^{(1/2)} ∝ P_j .
```

So the block, as a function of the transverse angle, is a **Gegenbauer expansion**

```
𝔽(u, v, θ) = Σ_j R_{(Δ_j, j)}(u, v) · C_j^{(1/2)}(cos θ) ,
```

where `(u, v, cos θ)` are the 3 cross-ratios; the angle carries the spin decomposition.

## Spin dependence

`j` enters the radial block via the `SO(2,1)×SO(3)` quadratic Casimir
`C = Δ_j(Δ_j−1) + j(j+1)`.

## Honest boundary

The Gegenbauer angular factor is standard and cheap; the genuinely missing piece is the
**2-variable radial block** `R_{(Δ_j,j)}(u,v)`. In the free theory it is a Gaussian integral (a
benchmark); the interacting `R_j` is computed nowhere — but see `BLOCK-DERIVATION.md`: it is a
well-posed PDE problem.

Sources: Billò–Gonçalves–Lauria–Meineri 1601.02883 (transverse harmonic decomposition);
arXiv:2501.06900; SO(3) addition theorem / Gegenbauer.
