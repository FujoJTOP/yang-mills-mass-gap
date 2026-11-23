# Defect-bootstrap setup: the crossing equations

> CPB-attack probe. Verdict: the crossing equation is written; **the first test is not
> runnable yet** — the 4D non-SUSY line-defect blocks must be built first.

## Setup

Line along `x³`; preserved symmetry **`SO(2,1)×SO(3)`** (1D conformal on the line ×
transverse `R³` rotations). **CORRECTION (2026-11-23):** a line in 4D is codim-3, so the
transverse space is `R³` and the group is `SO(3)`; the earlier `SO(2,1)×SO(2)` is the
**codim-2** case (a line in 3D, or a surface in 4D). This is the correct group.

**Multi-variable:** a 4D line has **three** invariants (8 coordinates − 5 symmetries) —
transverse-distance ratio, along-line ratio, transverse angle — not the single `ξ` of the
codim-2 / 3D case. So the 4D non-SUSY blocks are **multi-variable** (3 cross-ratios); this
is the main reason they are missing. Defect operators `O_i(τ)` carry dimension `Δ_i` and
`SO(3)` spin `j_i`. The `Z_N` 1-form symmetry acts on the **line** (charge `q∈Z_N`, `Z_N`
fusion of lines) and **trivially on local (bulk) operators**.

## Crossing (bulk–defect, the well-defined object)

For a neutral bulk scalar `φ`:

```
⟨φ(x₁)φ(x₂)⟩_def = Σ_O C^{φφ}_O 𝒢_O(ξ)     (bulk OPE channel)
                  = Σ_i a_i² 𝔽_i(ξ)          (defect OPE channel),
```

`ξ` = defect cross-ratio; `𝒢_O` bulk block, `𝔽_i` defect block. Crossing: the two
expansions are equal. The `Z_N` charge `q` is a superselection label; the mixed `Z_N–θ`
anomaly (5D symmetry TQFT) adds a **1D 't Hooft anomaly constraint (a "twist")** on the
charged line's defect spectrum.

## Free benchmark

Free gluons (`g=0`); the Wilson line's defect CFT is known perturbatively — its primaries
are built from the gluon modes restricted to the line, `q` fixed by the representation.
This is the exact crossing solution to beat.

## What exclusion requires

The bootstrap must show that, with the `Z_N` line and its anomaly, the **only** crossing
solution is the free one — a **uniqueness** statement, strictly stronger than the usual
bounds.

## Honest boundary

**Closed-form 4D line-defect conformal blocks are NOT available.** The 2023–2025
line-defect bootstrap is developed mainly for SUSY (`N=4 SYM`) lines and for boundary
bootstrap in 3D; the **non-SUSY 4D defect blocks are the missing ingredient**. So the first
test is not yet runnable — the blocks must be built first.

Sources: line-defect bootstrap (N=4, ~2020–2025); boundary/defect bootstrap; GKSW anomaly;
symmetry TQFT.
