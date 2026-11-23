# Free (`g=0`) Wilson-line defect 4-point

> CPB-attack probe (sub-task 2 of round 3). Verdict: the Gaussian structure is explicit; the
> exact free normalisations need the bulk–defect mode expansion.

## Free defect 4-point

Line along `τ`, defect operators in rep `R'`, ordered `τ₁>…>τ₄`, cross-ratio
`χ = (τ₁₂τ₃₄)/(τ₁₃τ₂₄)`. At `g=0` the line is Gaussian (generalized-free), and Wick gives,
for dimension `Δ_D`,

```
⟨D₁D₂D₃D₄⟩ = (τ₁₃τ₂₄)^{−2Δ_D} · g(χ) ,
g(χ) = 1 + χ^{2Δ_D} + (1−χ)^{2Δ_D}          (Δ_D = 2 for the displacement D^a).
```

Crossing-symmetric under `χ → 1−χ`; the only free data are `Δ_D` and the normalisation.

## OPE content

- **Defect channel:** identity + double-trace `[DD]_n` (`Δ = 2Δ_D+2n = 4+2n`); **no**
  single-trace defect primaries in the Gaussian limit.
- **Bulk channel:** the free 4-point of `tr F²` (`Δ=4`) and descendants, entering via the
  bulk–defect OPE.

## `Z_N` charge

`q = N-ality(R) mod N` is a superselection label of the **line**; `D^a` (adjoint) is `q'=0`.
At `g=0` path ordering trivialises, so the `Z_N` 't Hooft-anomaly phase is **absent** in the
benchmark — it appears only at `g≠0`.

## Honest boundary

This is the Gaussian structure. The exact free normalisations (the `D`–gluon mode map, the
`tr F²` coefficient) need the explicit bulk–defect mode expansion; the full free spectrum (all
primaries + coefficients) is not given here.

Sources: arXiv:2501.06900, 2208.11715, 1806.01862.
