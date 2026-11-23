# Defect-of-defect (twisted) bootstrap: concrete setup

> CPB-attack probe (round 5, reuse). Verdict: the **line-changing-local-operator** reading is
> standard SDPB with a phase-twisted kernel — **not** outside SDPB; the **point-defect-on-line**
> reading is the genuinely harder object.

## Two readings

- **(i) Line-changing local operator** `O_{q→q'}`: a defect operator on the line jumping the
  `Z_N` charge by `Δq=q'−q`, living in `Hom(H^{(q)}, H^{(q')})`. On a closed line insertions must
  sum to zero, so they come in **pairs**; the simplest nontrivial object is the 4-point of two
  junctions, `q₁→q₂→q₃→q₄→q₁`.
- **(ii) Genuine point defect on the line** (a 0d defect with its own defect data): this breaks
  the 1D conformal group — a separate, harder object.

## Crossing (reading i)

Four points on the line, `χ = τ₁₂τ₃₄/τ₁₃τ₂₄`. Crossing is a **braid** of junctions; the anomaly
2-cocycle `ω(q_i,q_j) = e^{2πi k q_i q_j/N}` multiplies OPE coefficients:

```
G(χ) = Σ_k C₁₂^k C₃₄^k G_k(χ) = Σ_k ω·C₁₄^k C₂₃^k G_k(1−χ) .
```

Unlike neutral externals (`q_i=0 ⟹ ω=1`), here `ω≠1`.

## Blocks

For the line-changing-local reading: still **one cross-ratio**, and the blocks are the **same 1D
hypergeometric** `χ^Δ ₂F₁(Δ,Δ;2Δ;χ)` (already `1e-15`). The twist is entirely in the crossing
matrix / OPE phases, **not** in the blocks — **no new block** needed.

## Machinery

Standard 1D defect bootstrap in **SDPB with a phase-twisted crossing kernel** — inserting `ω` is
free. **Not** outside SDPB. New work: pin the `k`/θ convention and the `Δ_{q→q'}`.

## Well-posed?

**Yes** for reading (i): one cross-ratio, known blocks, standard crossing with a fixed phase;
needs only charge-closure (pairs) and the BF convention. The point-defect-on-line (reading ii)
breaks the 1D conformal group — the "problematic" piece if meant literally.

## Honest boundary

Makes the `Z_N` grading/anomaly **visible** (neutral was blind) and is runnable now, but it
constrains the **defect** spectrum/phases (a selection rule), not the **bulk** gap directly; a
proof (island vs exclusion) still needs the 1-form layer.

Sources: GKSW arXiv:1412.5148; Cordova–Tachikawa; arXiv:2501.06900; arXiv:2208.11715.
