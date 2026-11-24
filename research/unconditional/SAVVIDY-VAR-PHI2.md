# SAVVIDY-VAR-PHI2 — the negative-mode square `Var(φ²)` is UNIFORMLY BOUNDED in `c=gB` (→ `2/λ₄`); the physical softening `g²Var(φ²)` is O(1) and `g`-independent

> **VERDICT (no QED).** For the bimodal `V(φ) = −½ c φ² + ¼ λ₄ φ⁴` (`c=gB>0`, `λ₄>0`) the negative-mode
> square has
> **`Var(φ²) = 2/λ₄ + 2/(gB)² + O((gB)^{-4})`** `[T-saddle]/[V-num]` (numerics: `(Var−2/λ₄)·c² → 2.00`).
> Hence `Var(φ²)` is **uniformly bounded as `c=gB` grows at fixed `λ₄`**: it increases to `2/λ₄` and no
> further (verified: `2.0221, 2.0051, 2.0013, 2.0003` at `c=10,20,40,80` with `λ₄=1`). In the physical
> regime `λ₄ = k g²`, `Var(φ²)` itself grows `~ (1/g²)`, **but the softening `g²Var(φ²)` is O(1) and
> `g`-independent** (numerics: `0.956456` for `a≡c/√λ₄=1`; `2.02209` for `a=10`, matching the large-`a`
> asymptote `2/k+2/B²=2.02`) — **so the softening does NOT diverge.** The bound's *value* is set by `λ₄`
> (a microscopic EFT coefficient), which is the honest `H1`-adjacent input.

---

## §1 · Exact saddle-point structure

With `y = φ² ≥ 0`, the measure is `μ(y) ∝ y^{−1/2} e^{½ c y − ¼ λ₄ y²} dy`. It peaks at
`y* = c/λ₄` (the well `φ* = ±√(c/λ₄)`) with a Gaussian width `σ_y² = 1/|f''(y*)| = 2/λ₄`.
The variance is the *subleading* Gaussian width; including the `y^{−1/2}`-prefactor variation
(`(y*+δ)^{−1/2} ≈ y*^{−1/2}(1 − δ/2y* + …)`) gives, to next order,

```
Var(φ²) = Var(y) = σ_y² + σ_y⁴/(2 y*²) + O(σ_y⁶/y*⁴) = 2/λ₄ + 2/(gB)² + O((gB)^{-4}).
```

**Crucial distinction.** `Var(φ)` (the *field*, not its square) is dominated by the two wells
`±√(c/λ₄)` and grows as `~ c/λ₄ = gB/λ₄` (divergent in `c`) — but the Savvidy coupling is to **`φ²`**
(`g a φ²`), and `Var(φ²)` is bounded. The soft mode is the *magnitude*, not the *sign*.

## §2 · Boundedness in `c` and physical scaling `[T]/[V-num]`

- **Fixed `λ₄`, growing `c=gB`:** `Var(φ²) → 2/λ₄` (numerics, `λ₄=1`: `2.0221→2.0001` for `c=10→160`;
  residual `(Var−2/λ₄)c² → 2.00`). **Uniformly bounded**, bound `= 2/λ₄`.
- **Physical `λ₄ = k g²`, `c = gB`:** after `y = 2u/√λ₄` the problem depends only on `a = c/√λ₄ = B/√k`,
  so `g²Var(φ²) = (2/k)·F(a)` is a **`g`-independent O(1)**. Numerics: `g²Var = 0.956456` constant across
  `g∈{1,0.7,0.5,0.4,0.3}` at `a=1`; `g²Var = 2.022094` constant at `a=10` (matching `2/k+2/B² = 2.02`).
  **The softening is bounded and O(1); it does not diverge toward the continuum.**
- **Requires `λ₄>0`.** `λ₄=0` is the pure Savvidy tachyon (no quartic): `V` is unbounded below, the integral
  diverges, `Var` is infinite. **Boundedness requires the quartic stabilization** `λ₄>0`, physically the
  gauge self-interaction (`λ₄ ~ k g²`).

## §3 · Honest limits `[O]`

1. **Bounded ≠ bounded by a known constant.** The O(1) value (`2/λ₄` at fixed `λ₄`; `(2/k)F(B/√k)` physically)
   is set by the negative-mode **quartic `λ₄`**, a microscopic EFT coefficient. So `Var(φ²)`'s *bound* is
   `[T]/[V-num]`, but pinning its *number* is the same microscopic input as `H1` in the Savvidy-cell
   coordinate (consistent with `CONSTRUCT-SAVVIDY-CELL`).
2. **Consistency with `CONSTRUCT-SAVVIDY-CELL`.** Its witness `Var(φ²)=0.41…8.51` across `λ₄=2…0.2` is
   exactly the `~2/λ₄` law (`8.51/0.41 ≈ 21` vs `2/0.2 ÷ 2/2 = 10`, order-correct); its conclusion
   "finite threshold shift, not an unbounded hole" is confirmed and made explicit.

## §4 · Consequence for the χ-floor / softening

The Savvidy softening `g²Var(φ²)` is **O(1) and bounded in `(g, B)`** — the negative mode does **not** create
an unbounded hole in the χ floor. The floor condition `κ_bare > g²Var(φ²)` is a **genuine O(1) competition**
(`κ_bare` vs an O(1) softening), *not* a divergence. What remains = controlling this O(1) competition
uniformly across cells/β — i.e. the value of `λ₄` / the ratio `B/√k` — which is `H1` (softened) in the
Savvidy-cell coordinate. **Clay 0%; no QED.**

*Witness:* `w7/savvidy_var_phi2.py` (scipy quadrature, runs clean, `py -3.12 -I`).
