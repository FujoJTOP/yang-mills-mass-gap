# Block: two agents disagreed — settled by direct numerics

> The two probes contradicted each other (`U-ODE.md` says `f` is the block; the normalization probe
> says `f` fails and `Q` is the block). **Neither raw function is the block.** Settled by evaluating
> both: each is right at *one* OPE end and wrong at the other.

## What both agree on (settled)

- **Shape:** `𝔽_j(u,θ) = f(u)·P_j(cos θ)`, with `u` the **one** `K_z`-invariant H² cross-ratio and
  `SO(3)` `P_j(cos θ)`. The invariants of `(z₁,ρ₁,z₂,ρ₂)` are `u` + `cos θ`. ✓ (`RESOLUTION.md`)
- **`Δ_φ`:** enters **only** through the `u→0` (bulk-OPE) normalization and the dynamical `a_O²`;
  the defect-channel shape is `Δ_φ`-independent. ✓
- **Conformal Casimir:** `SO(2,1)` quadratic Casimir = hypergeometric-type operator on `u`
  (regular singular at `u=0,1`), **not** the naive H² Laplacian `(z²−1)f''+zf'`.
- **`ξ = ρ₁ρ₂/z₁₂²` is not `K_z`-invariant** — it is the *defect-channel* variable (`Σc_i²ξ^{Δ_i}`),
  **not** the full block's variable. (`RADIAL-SOLVE.md`'s `G(ξ)` is a *different* correlator.)

## The disagreement, and why neither raw form works

Two OPE ends must both hold:
- **defect OPE (`u→1`):** block decays as `(1−u)^{Δ_j}`.
- **bulk OPE (`u→0`):** block has the standard bulk power (the `ρ^{−Δ_φ}`-reduced block `→ 1`).

Direct evaluation (`scipy.hyp2f1`, `D=Δ_j`, `u=(cosh η−1)/(cosh η+1)`):

| `Δ_j` | `u` | `f = ₂F₁(Δ,1−Δ;1;−u/(1−u))` | `Q_{Δ−1}(cosh η)` | ratio |
|---|---|---|---|---|
| 2 | 0.5 | 3.00 | 3.97×10⁻² | 75 |
| 2 | 0.95 | 39.0 | 2.19×10⁻⁴ | 1.8×10⁵ |
| 3 | 0.8 | 121 | 1.85×10⁻⁴ | 6.5×10⁵ |
| 3 | 0.95 | 2281 | 2.25×10⁻⁶ | 1.0×10⁹ |

- **`f` (U-ODE):** ✓ `f→1` at `u→0` (bulk OPE right); ✗ **grows** at `u→1` (dominated by the
  `(1−u)^{1−Δ_j}` branch) — **wrong defect OPE**.
- **`Q_{Δ−1}(cosh η)` (normalization probe):** ✓ decays `~(1−u)^{Δ_j}` at `u→1` (defect OPE right);
  ✗ **log-diverges** at `u→0` — **wrong bulk OPE**.
- The ratio is `u`- and `Δ`-dependent (not a constant), so `f ≠ const·Q`: **genuinely different.**

## The correct block (up to the open derivation)

The physical block is the **second-kind (`Q`-type) solution carrying the `u→1` decay AND a fixed
`u→0` normalization** — i.e. `(1−u)^{Δ_j} × ₂F₁(…; u)` (or the standard defect block
`ξ^{Δ_j/2}₂F₁((Δ_j+Δ_φ)/2,(Δ_j−Δ_φ)/2;Δ_j;ξ)` written in `u`), **not** either raw function above. The
`u→0` boundary condition that picks the `Q`-branch coefficients is exactly the bulk data
`(Δ_φ, C^{φφ}_O)` — this is where `Δ_φ` lives.

## Honest boundary

The **shape is settled** (`f(u)P_j`, one invariant `u`, hypergeometric conformal Casimir,
`Δ_φ` only in the bulk-OPE normalization). The **explicit `u`-form with both OPE limits and the
`Δ_φ` normalization is open** — `U-ODE.md`'s `f` is wrong at `u→1`, `G(ξ)` tests the wrong
correlator, and `Q` alone is wrong at `u→0`.
