# Cross-checks for a candidate 4D non-SUSY line-defect block

> CPB-attack probe (sub-task 12). Verdict: a list of *necessary* checks; the most
> diagnostic is the small-ξ exponent; the free-theory benchmark is the only exact known value.

## Checks (a candidate must satisfy all)

1. **Small-ξ (defect OPE):** `𝔽_i(ξ) ~ ξ^{Δ_i/2}·(1+…)` — the leading power is set by the
   defect-operator dimension `Δ_i`; a wrong prefactor/spin structure shows up here first.
2. **ξ→1 (bulk OPE):** the block must reproduce the bulk-OPE singularity controlled by the
   exchanged bulk dimension `Δ_O`; a mismatch falsifies the bulk–defect normalisation.
3. **Large-Δ:** `𝔽_i ~ exp(−Δ·σ(ξ))` with the defect geodesic exponent `σ(ξ)`; wrong
   asymptotic decay ⟹ wrong block.
4. **Free theory:** the `g=0` Wilson-line two-point function is the **exact** benchmark; a
   correct block reduces to it.

## Which pass

Without the explicit candidate block these are *tests to apply*, not verdicts. The one a
candidate should pass **first** is the free-theory benchmark — it is the only exact known
value.

## Most diagnostic

The **small-ξ exponent**: fixed by the defect spectrum, cheap to compute; a wrong block
almost always fails it before the free-theory normalisation does.

## Honest boundary

Necessary, not sufficient; passing them does not prove the block. Closed-form 4D non-SUSY
defect blocks may not exist (only integral/numerical representations), so the free-theory
benchmark itself may need numerical evaluation.
