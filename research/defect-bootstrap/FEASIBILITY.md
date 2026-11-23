# Defect-bootstrap feasibility (4D line defects, 1-form symmetry)

> CPB-attack probe. Verdict: **needs-new-solver (full 1-form route)**; a **simplified
> defect-crossing test is feasible now**.

## Equations / blocks

A line-defect bootstrap uses the crossing of a bulk 4-point function **with a defect**
(one bulk scalar, defect conformal group `SO(2,2)×SO(2)`), expanded in **defect conformal
blocks**, plus the defect OPE and the bulk–defect OPE. The 1-form symmetry does **not** act
on local operators, so the *local* defect crossing does not see it directly; it constrains
the **line-operator** spectrum (line OPE / fusion), a different, less-standard set
(fusion-category / anomaly-inflow constraints).

## Solver

SDPB (Simmons-Duffin 2015) handles the **standard** defect crossing — it is the usual
positive-semidefinite form, defect blocks in 3D/4D are computable, and line-defect
bootstraps have been run since ~2018. The **1-form-symmetry constraint** (line-OPE / fusion
+ anomaly inflow) is a 2023+ frontier and needs new machinery (symmetry-TQFT /
generalized-symmetry constraints) — **not** out-of-the-box SDPB.

## Proof

The bootstrap yields **bounds / islands** (rigorous up to interval-arithmetic precision). A
*proof* of "the only solution is free" needs exclusion of **all** interacting solutions —
stronger than an island. Possible in principle (rigorous bootstrap), but typically hard.

## Cost / maturity

4D defect bootstrap: current frontier; SDPB exists; the full 1-form version needs a new
solver / constraint layer.

## Verdict

**needs-new-solver** for the full 1-form-symmetry route; a **simplified** defect-crossing
test (bulk 4-pt with a `Z_N`-charged line defect, no fusion layer) is **feasible now**.

## Honest boundary

Assessment only, no proof; the generalized-symmetry bootstrap is immature and may give
bounds, not uniqueness.

Sources: Billò–Gonçalves–Lauria–Meineri, *Defects in CFT* (arXiv:1601.02883);
Simmons-Duffin, *SDPB* (arXiv:1502.02033); Poland–Rychkov, *The Conformal Bootstrap*, Rev.
Mod. Phys. 91 (2019).
