# 1-form-symmetry SDPB layer

> CPB-attack probe (round 5). Verdict: the "needs-new-solver" is **real but bounded** — sector
> refinement is standard SDPB; the projective phase needs a **modified linear-constraint layer**;
> the **twisted** sector is outside SDPB entirely.

## (1) Constraints

The `Z_N` 1-form acts on lines; `H_def` graded by `q∈Z_N`, defect OPE additive mod `N`. The
`Z_N`–θ anomaly is a 5d BF symmetry TQFT; anomaly inflow gives a **1d 't Hooft anomaly** on the
worldline, classified by `H²(Z_N,U(1)) = Z_N` at level `k`. It appears as a **projective
2-cocycle** `ω(q_i,q_j) = e^{2πi k q_i q_j/N}` in defect fusion/crossing — non-trivial class,
**irremovable** by rephasing. Crossing splits into `Z_N` sectors with a **fixed** inter-sector
phase (a selection rule, not a parameter).

## (2) Encoding

Standard SDPB: real `Σ a_i²F_i = Σ C_O²G_O`, `a_i², C_O² ≥ 0`. Here the crossing is **complex**:
one equation per total charge `q`, `ω` phase fixed; coefficients live in a projective rep, so
their real/imag parts are linked by a **linear** relation, not independent.

## (3) Genuinely new

- **(a) Sector refinement** — more equations, standard SDPB-solvable.
- **(b) Projective phase** shrinks the positivity cone from "independent `a_i²≥0`" to
  "`|a_i|²≥0` with fixed relative phases"; SDPB's semidefinite ansatz may not represent this cone
  efficiently — a **modified linear-constraint layer** is needed.
- **(c) Twisted (line-changing) sector** — a defect-of-defect bootstrap, **outside standard
  SDPB** entirely.

## Honest boundary

Level/normalisation convention-dependent; no built 4D pure-gauge twisted solver.

Sources: GKSW arXiv:1412.5148; Cordova–Tachikawa (symmetry TQFT); Simmons-Duffin, SDPB
arXiv:1502.02033; arXiv:2501.06900.
