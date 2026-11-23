# Can the 3-cross-ratio block be derived from the 1D blocks? — decisive reframing

> CPB-attack probe (sub-task 7, the lever). Verdict: **no algebraic shortcut, but the missing
> block is a well-posed, numerically solvable PDE problem** — missing in *closed form*, not in
> *computability*.

## No algebraic shortcut

The 3-variable block solves the combined `SO(2,1)×SO(3)` Casimir PDE system — a multi-variable
PDE, **not** a product or finite combination of the 1D `SO(2,1)` block and the bulk `SO(2,2)`
block. The cross-ratio `ξ = |y|²/(t²+|y|²)` mixes longitudinal and transverse, so the Casimir
does **not** separate into (1D block)×(bulk block). No defect partial-wave formula reduces it to
lower-dimensional blocks.

## But the failure is algebraic, not computational

The 3-variable block is **determined** by (i) the Casimir PDE system and (ii) the OPE limits —
and those limits **are** built from the known blocks:

- the bulk→defect limit gives `1D block × bulk–defect OPE coefficients`;
- the far-defect limit gives the bulk block.

So it is a **well-posed PDE + boundary-data problem**, solvable numerically (spectral method /
conformal-block recursion) — exactly what SDPB needs (sampled vectors, not closed forms).

## Decisive

**The "missing block" is missing in *closed form*, not in *computability*.** The real task is
to set up and solve the Casimir PDE with the known OPE boundary data — **feasible and bounded**,
not an open-ended hunt. A closed form may not exist; a numerical block almost certainly does.

## Honest boundary

"Solvable numerically" is a claim about well-posedness, not a verified computation; the
transverse `SO(3)` spin structure and the 3-variable boundary conditions still need explicit
derivation.
