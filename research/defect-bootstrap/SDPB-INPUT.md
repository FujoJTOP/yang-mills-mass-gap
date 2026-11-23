# SDPB input assembly (contingent on the block)

> CPB-attack probe (round 5). Verdict: the SDPB input is assembled concretely, but **depends on
> the unsettled block** (`U-ODE.md` vs `RADIAL-SOLVE.md`); untested.

## Grid

`SO(3)` is exact (`P_j`), so project the crossing onto each transverse partial wave `j` (multiply
by `P_j(cos θ)`, integrate over `cos θ`): the problem splits into **independent 1D problems per
`(j,q)`**. Nodes: `N_u = 60`, mapped-Chebyshev in `w = log((1−u)/u)` over `u ∈ [10⁻⁴, 1−10⁻⁴]`
(dense at both OPE ends); `j = 0…J` (`J≈10`); `q ∈ Z_N` (SU(3): 3 sectors).

## Block vectors

- **Defect:** `F_{Δ,j}(u_k) = f_Δ(u_k)`, `f_Δ = ₂F₁(Δ,1−Δ;1;−u/(1−u))` (`U-ODE.md`), one column
  per allowed `(Δ,j)`.
- **Bulk:** `G_{O,j}(u_k)` = `j`-th partial wave of the Dolan–Osborn block (in `z,z̄` restricted to
  the defect kinematics), one per `O` = identity, `T_{μν}` (`ℓ=2`), `tr F²` (`ℓ=0`).

## Crossing

Per `(q,j)`: `Σ_{i:q_i=q} a_i² F_i(u_k) = Σ_O C_O² G_{O,j}(u_k)` at each `u_k`. Junction phases
multiply `a_i² → e^{2πi k q_i q_j/N} a_i²` **before** projection.

## Constraints / objective

`Re(e^{−2iφ_i} Λ(F_i)) ≥ 0` (`φ_i` from `ω`); `Λ(G_O) ≥ 0`; `Λ(F_id) = 1`; `Im`-part `= 0` per
`q`. Objective: `max Λ(F_test)` for a fixed trial block (a chosen `Δ_test`) — the standard
exclusion test.

## Honest boundary

The bulk partial waves `G_{O,j}` depend on the defect-kinematics convention (**unverified** — the
`z,z̄ ↔ u` map is the weak point); `φ_i` depend on level `k` / BF normalisation; `u`-nodes must
resolve `(1−u)^{Δ_j}` near `u=1`. **Untested; no code run.** And the defect block `f_Δ` is itself
in conflict (`U-ODE.md` vs `RADIAL-SOLVE.md`) — the resolution agent must settle it first.
