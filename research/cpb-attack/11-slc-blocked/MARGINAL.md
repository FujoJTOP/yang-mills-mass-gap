# Marginal blocked Hessian — exact fine Hessian + coarse restriction

> CPB-attack probe. Verdict: **NEEDS-X (negative lean)** — the Gaussian/one-loop marginal
> is ill-defined; the exact marginal is uncomputed.

## Method

Built the **exact** fine Wilson Hessian (8192-dim, L=4) by per-plaquette finite
differences of `−(β/3)Re tr U_p` in left variations.

## Spectrum (β=2.5)

- Full fine Hessian: **min ≈ −2.7**, 20/20 smallest negative, 100% of configs
  (L=2 control ≈ −2.6). **Indefinite.**
- Coarse fraction of the lowest eigenvectors (b=2 block-uniform projector): **0.05–0.10**
  → **~90–95% fluctuation / ~5–10% coarse**.
- Coarse restriction `PᵀHP` (512-dim): **min ≈ −0.3** (configs −0.33, −0.28), 20/20
  negative — ~9× smaller than the full min, but **still negative**.

## Implications

1. The Gaussian/one-loop marginal is **ill-defined**: the negative modes are ~95% in the
   fluctuation block `QᵀHQ`, so the Gaussian fluctuation integral diverges; the Schur
   complement is not a meaningful Gaussian marginal.
2. Blocking **reduces but does not remove** the Savvidy negativity: **−2.7 → −0.3**.
3. The marginal is **not computed**: it needs the conditional average over `A⊥`
   (constrained sampling). `PᵀHP` is pointwise, not the conditional average, and
   `Var ≥ 0` only subtracts — so the marginal can only be *more* negative than `E[PᵀHP]`.

## Honest boundary

L=4 → coarse `L'=2` is small (finite-size); pointwise ≠ conditional average; the Gaussian
proxy is invalid here; the exact marginal needs constrained sampling or the exact blocked
second variation. **Do not read this as closing or opening the entry.**

Scripts: `full_hessian.py`, `coarse_fraction.py`, `coarse_block_hessian.py`.
