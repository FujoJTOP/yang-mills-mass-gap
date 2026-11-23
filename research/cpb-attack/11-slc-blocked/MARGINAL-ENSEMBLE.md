# Ensemble-average coarse-projected Hessian (the E-part of the marginal)

> L=4, β=2.5, b=2. Script `marginal_ensemble.py`. This tests the *ensemble
> average* `Ĥ = E[PᵀHP]`, not the pointwise min.

## Result

- **Per-config** `PᵀHP` min: −0.310, −0.358, −0.299, −0.278 (4/4 negative).
- **Ensemble average** `Ĥ = E[PᵀHP]` (dim 512): **min = +0.043, 0/512 negative**;
  lowest 10: 0.043, 0.049, 0.056, 0.060, 0.062, 0.065, 0.067, 0.071, 0.074, 0.077.
- Cold `U=1`: min = +0.0000 (gauge zero mode).

## Reading

The pointwise coarse Hessian is negative (−0.3) at **every** config, but the **ensemble
average is positive definite** (min +0.043) — the negative modes are **config-dependent
and average out**. This **corrects the previous "negative lean"**, which used the
pointwise min rather than the averaged operator. The marginal's E-part is the
*conditional* average `E[PᵀHP | A']`, whose ensemble mean is `Ĥ` (positive).

## Caveat

The `+0.043` min may be the **gauge zero mode** lifted by the finite-difference
(`EPS=1e-3`) error at thermal configs (at `U=1` the min is exactly 0). Excluding the
gauge modes could change the physical min.

## Verdict

**INCONCLUSIVE (not dead).** The E-part average is positive; the exact marginal still
needs the conditional `Var(∂S | A')`. The earlier "negative lean" is **not supported** by
the ensemble-averaged operator.

## Boundary

L=4 (`L'=2`, small); 4 configs; pointwise ≠ conditional average; gauge modes not removed.
