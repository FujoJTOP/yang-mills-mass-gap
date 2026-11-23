# 1D (codim-2) block: numerical computation

> CPB-attack probe (sub-task 11). The 1D block is computed to `1e-15`; it is a **building
> block**, not the full 4D (codim-3, 3-variable) block.

## Correction (matches corrected `SETUP.md`)

A line in 4D is **codim-3** (transverse `R³`, group `SO(2,1)×SO(3)`, **three invariants**), so
the true 4D non-SUSY block is **multi-variable** — the stated reason it is missing. What is
computed here is the **1D (codim-2) block**, single `ξ` — a building block, not the full 4D
object.

## Method

Exchanged-dimension `Δ_e` 1D block `G = z^{Δ_e} H`, with `H = ₂F₁(Δ_e,Δ_e;2Δ_e;z)` solving
`z(1−z)H'' + [2Δ_e − (2Δ_e+1)z]H' − Δ_e² H = 0`. Series recursion
`h_n = h_{n−1}·[(n−1)(n+2Δ_e−1)+Δ_e²]/[n(n+2Δ_e−1)]`; series vs `hyp2f1` agrees to ≤ `7×10⁻¹⁵`.

## Table `G(ξ)` at `ξ = 0.1, 0.3, 0.5, 0.7, 0.9`

| `Δ_e` | 0.1 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|
| 1 | 0.10536 | 0.35667 | 0.69315 | 1.20397 | 2.30259 (= `−ln(1−ξ)`) |
| 2 | 0.01110 | 0.12695 | 0.47665 | 1.41570 | 4.88562 |
| 3 | 0.00117 | 0.04517 | 0.32740 | 1.65904 | 10.24611 |
| 4 | 0.00012 | 0.01607 | 0.22480 | 1.94202 | 21.40270 |

## Limits

`ξ→0`: `G/ξ^{Δ_e} → 1.000` (leading twist). `ξ→1`: `c = a+b` ⟹ logarithmic divergence
(`Δ_e=1` is exactly `−ln(1−ξ)`).

## Free benchmark

Free massless scalar, straight line, collinear: `f_free(ξ) = 1/(16π²ξ)` (single
identity-defect exchange).

## Honest boundary

This is the **1D/codim-2** block (single `ξ`). The true 4D line (codim-3) block is
**3-variable** and remains the missing research object; the `SO(3)` transverse-spin
generalisation is not derived here.
