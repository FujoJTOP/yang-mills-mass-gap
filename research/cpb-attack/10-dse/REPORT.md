# Method 10 — Ward / Schwinger–Dyson (loop) equations + reflection positivity

> CPB-attack probe. Verdict: **DEAD as a standalone route** (reach assessment, not a disproof).

## The equation

The lattice Makeenko–Migdal (loop) equation at finite `N`:

```
∂_μ δ/δσ_{μν}(x) ⟨W(C)⟩  ∝  g² [ ⟨W(C₁)W(C₂)⟩_c + O(1/N²) ] .
```

It ties the loop-derivative of the 1-loop expectation to a **connected 2-loop**
correlator; the 2-point decay rate is tied to a 3-point object, and so on. **The
hierarchy never closes at finite `N`.** Only the `1/N²`-suppressed connected part
vanishes at large `N`, where the planar master-field equation closes
(Makeenko–Migdal 1979/81; Migdal 1983).

## Calibration

- **2D YM:** equations close; exact `W(C)=exp(−g²A/2)` (area law) — but 2D has no
  propagating degrees of freedom; the "gap" is topological, **not informative for 4D**.
- **Large `N`:** the planar equation closes but is not solved explicitly in 4D; no
  rigorous gap extraction (Chatterjee, *Yang–Mills for probabilists*).
- **Finite `N`:** open hierarchy.

## Positivity

Reflection positivity makes the loop functional positive-definite, but **gapped and
gapless functionals can both be positive-definite** — positivity constrains but does
not determine the decay rate. Osterwalder–Seiler reconstruction recovers the measure
from the moments, not the spectral gap.

## Verdict

**DEAD as a standalone route.** At `β≈2.5` the DSE + positivity do not force
exponential decay. The only conditional salvage — proving a decay/cluster input from
the DSE — *is* the crossover certificate itself.

## Honest boundary

A reach assessment, not a disproof; a new positivity–DSE theorem could change it.
Overlaps method 9 (area law) through the random-surface representation (CNS25a).
