# Conformal-block recursion for the defect `(ξ,θ)` block

> CPB-attack probe (round 4). Verdict: the recursion is **trivial** *if* the block is separable
> (product form); the crux is the `K`-invariance of `ξ`.

## Recursion (if separable)

The Casimir `𝕃²_{SO(2,1)}(ξ) + 𝕃²_{SO(3)}(θ)` acts on separate variables, so the block is a
**product** `𝔽_{Δ,j} = f_Δ(ξ)·P_j(cos θ)`. Then the recursion collapses to the 1D coefficient
recursion

```
h_n = h_{n−1}·[(n−1)(n+2Δ−1)+Δ²]/[n(n+2Δ−1)]      (f_Δ = ξ^Δ ₂F₁(Δ,Δ;2Δ;ξ)),
```

already validated to `1e-15` (`NUMERICS.md`), times the Legendre `P_j`. **Cost `O(N·j)`,
trivial.** This *is* the Zamolodchikov recursion for the 1D factor; no new machinery.

## Applicability — conditional

Two caveats. (i) The naive `ξ = ρ₁ρ₂/z₁₂²` is **not `K₃`-invariant**: under special conformal
along `z` it scales by `(1+cz₁)(1+cz₂)`. The true conformal invariant mixes `ξ` and `r`. So the
product form is exact only in the **`K`-invariant** cross-ratio, whose convention must come from
2501.06900 and be checked. (ii) If it fails at subleading order, the block is genuinely
3-variable and **no standard 3-variable Zamolodchikov recursion exists** — fallback is the
numerical spectral solve of the coupled PDE (`BLOCK-DERIVATION.md`).

## Honest boundary

The cheap recursion is conditional on separability; the crux is the `K`-invariance of `ξ`.

Sources: Zamolodchikov 1984 (block recursion); 1D conformal-block literature; arXiv:2501.06900
(defect conventions).
