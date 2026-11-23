# Method 3b/4b — blocked-scale SLC: collective plane-wave Hessian

> CPB-attack probe. Verdict: **needs-X (positive proxy)** — the first positive signal.

## Instrument

Collective plane-wave Hessian of the Wilson action: perturb **all** links of direction
`ν` by `exp(iε·cos(p·x)·T)` (`T` = random algebra direction); symmetric finite
difference. Sanity: at `U=1` the low-`p` Hessian is positive (mean +202, min ≈ 0).
Batched over L=6; 25 configs × (4 low-`p` `n=1` + 4 high-`p` `n=3`) × 4 `ν` × 4 random
color vectors.

## Results (L=6, fraction of modes with `H<0` = **0.000** in every cell)

| β | LOW-p mean / min | HIGH-p mean / min |
|---|---|---|
| 2.5 | +271.3 / +234.6 | +541.2 / +465.5 |
| 2.6 | +292.0 / +254.6 | +583.7 / +511.6 |
| 3.0 | +397.7 / +345.2 | +793.8 / +708.0 |

## Read against Method 4

The **per-link** Hessian has 49% negative blocks at `β=2.5` (Method 4), but **every
collective plane-wave mode is strictly positive**. So the Savvidy negatives are
**localized / high-momentum**, while the **low-momentum (blocking-surviving) collective
modes are stable** — consistent with the blocked action being convex because the unstable
directions are UV and get integrated out.

## Verdict

**needs-X (positive proxy).** The plane-wave collective Hessian is strictly positive at
all three β, which *supports* blocked-SLC being alive; but plane waves do not span
localized modes, and the **marginal** blocked Hessian (with the conditional-variance
subtraction `E[∂²F] − Var(∂F)`) is not computed. Not a proof.

## Honest boundary

Fine collective Hessian ≠ marginal blocked Hessian; the proper object needs constrained
sampling or the exact second variation of the blocked plaquette action.
