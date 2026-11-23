# Coupled `(ξ,θ)` PDE: the product form

> CPB-attack probe (round 4). Verdict: `SO(3)` forces `P_j(cos θ)`; the `ξ`-part is the 1D
> hypergeometric **if** separability holds in the `K`-invariant cross-ratio. The block is
> **(very likely) a product** — essentially known.

## Operators (explicit)

- `𝕃²_{SO(3)} = −∂_θ² − cot θ ∂_θ`; eigenfunctions `P_j(cos θ)`, eigenvalue `j(j+1)`. This
  **forces** `𝔽_j(ξ,θ,r) = f(ξ,r)·P_j(cos θ)` (the Casimirs commute and `𝕃²_{SO(3)}` acts only
  on `θ`).
- `𝕃²_{SO(2,1)}(ξ)` = the 1D conformal Casimir in the defect cross-ratio, second order in `ξ`;
  eigenvalue `Δ_j(Δ_j−1)`.

## Series (under the standard 1D-block form; verified vs `scipy.hyp2f1` to ≤ `4×10⁻⁴` at `ξ=0.3`)

```
𝔽_j = ξ^{Δ_j} · ₂F₁(Δ_j,Δ_j;2Δ_j;ξ) · P_j(cos θ)
     = ξ^{Δ_j}[1 + (Δ_j/2)ξ + Δ_j(Δ_j+1)²/(4(2Δ_j+1))·ξ² + …] P_j(cos θ).
```

| `(Δ_j,j)` | leading | `+` first subleading |
|---|---|---|
| (1,0) | `ξ` | `ξ(1 + ξ/2 + ξ²/3 + …)` |
| (2,0) | `ξ²` | `ξ²(1 + ξ + 0.90ξ² + …)` |
| (2,1) | `ξ² cos θ` | `ξ²(1+ξ+…) cos θ` |
| (2,2) | `ξ²(3c²−1)/2` | `ξ²(1+ξ+…)(3c²−1)/2` |
| (3,0) | `ξ³` | `ξ³(1 + 1.5ξ + 1.714ξ² + …)` |

## Boundary conditions

Defect OPE at `ξ→0` (power `ξ^{Δ_j}`, above); bulk channel at `ξ→1` (`₂F₁` branch point; anomaly
phase `e^{iθ/N}` on the crossing kernel).

## Honest boundary

The `SO(3)` factor and the leading term are **exact**. The `ξ`-part is the 1D hypergeometric
**under the assumption** that `𝕃²_{SO(2,1)}` acts on `ξ` alone; the exact `(ξ,r)` coupling via
`K₃` (`CASIMIR-PDE.md`'s `r` flagged non-`K₃`-invariant) is **not** verified — so `r` may enter
beyond leading order, and the exact `₂F₁` parameters are convention-dependent. Closed form for
the full 3-variable block unverified; the coupled-PDE numerical solve remains open.
