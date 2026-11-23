# 4D line-defect CFT kinematics (corrected)

> CPB-attack probe (sub-task 2). **Corrects the group and flags the multi-variable block.**

## Group correction (firm)

A line in 4D is **codimension 3**: the transverse space is `R³` and the preserved group is
**`SO(2,1)×SO(3)`** (1D conformal on the line × transverse rotations). The `SO(2,1)×SO(2)`
in `ROUTE.md`/`SETUP.md` is the **codim-2** case (a line in 3D, or a surface in 4D). Use
`SO(3)`.

## Conventions

Euclidean `R⁴`, `x=(z,y)`, `z∈R` along the line, `y∈R³` transverse, `ρ=|y|`. Defect
operators `O_i(τ)` carry `(Δ_i, j_i)`, `j_i` an `SO(3)` spin. Bulk scalar `φ`, dimension
`Δ_φ`.

## Two-point functions

- Defect: `⟨O_i(τ)O_j(0)⟩ = δ_{ij}δ_{R_iR_j}/|τ|^{2Δ_i}`.
- Bulk one-point (`z`-translation ⟹ no `z`-dependence): `⟨φ(x)⟩ = a_φ/ρ^{Δ_φ}`.
- Bulk–defect: `⟨φ(x)O_i(τ)⟩ = a_i (ρ²+(z−τ)²)^{−Δ_φ/2} × (SO(3) structure)`.

## Crossing

`⟨φ(x₁)φ(x₂)⟩_def = ρ₁^{−Δ_φ}ρ₂^{−Δ_φ} G(ξ) = Σ_O C^{φφ}_O 𝒢_O(ξ)` (bulk)
`= Σ_i a_i² 𝔽_i(ξ)` (defect).

## Cross-ratios

A 4D line has **three** invariants (8 coordinates − 5 symmetries): transverse-distance
ratio, along-line ratio, transverse angle — **not** the single `ξ` of codim-2 / 3D. The
**multi-variable blocks** (3 cross-ratios) are the main reason the 4D non-SUSY blocks are
missing.

## Honest boundary

Group correction is firm; verify the explicit 3-cross-ratio parameterisation against
Billò–Gonçalves–Lauria–Meineri, *Defects in CFT* (1601.02883) before using verbatim.
