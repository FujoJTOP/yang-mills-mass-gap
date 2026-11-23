# Log-critical tree-weight accounting for 4D YM (SPDE)

> CPB-attack probe. Verdict: **needs-further** (plausible at one loop; unproved beyond).

## Scaling

Parabolic metric (`[x]=1`, `[t]=2`); space-time white noise; Gaussian-free-field regularity
`α = (3−d)/2 − κ`, so `α = −1/2 − κ` at `d=4` (one unit rougher than the 3D CCHS field,
`α = −κ`). Equilibrium marginality: `[g]=0` at `d=4`.

## Tree accounting

Naive power counting for the YM nonlinearities (quadratic `A∧∂A`, `j=1`; cubic
`A∧A∧A`, `j=0`) gives margin `(k−1)α − j + 2 = +0.5` (quadratic) and `+1` (cubic) at
`d=4` — **positive**, suggesting marginal-subcriticality.

**But this counting is unreliable for a gauge theory:** the cubic term is part of the
marginal `F²`, and gauge invariance reorganises the homogeneity, so the naive `A³` weight
is not the true one. The parabolic status is therefore **not settled by power counting**.

At equilibrium marginality the weight-0 trees are `F²`-type: coupling log (`g log F²`),
field-strength log (`log F²`), plus gauge-fixing/ghost logs if the gauge is fixed.

## Gauge covariance

`F²`-type logs are gauge-covariant by the Ward identity and close at one loop. The
obstacle is not the individual counterterm but:

1. survival of the gauge-covariant (unfixed, parallel-transport) frame under
   renormalisation;
2. the **invariant measure**: the 4D Φ⁴ SPDE is well-posed yet its invariant measure is
   **trivial** (Aizenman–Duminil-Copin) — so local well-posedness ≠ a non-trivial measure.

## Verdict

**needs-further.** The log-critical structure is plausible at one loop, unproved beyond;
the deeper obstruction is the **invariant measure (UV)**, not the local log-counterterms.

## Boundary

Pen-and-paper; the gauge structure makes the naive margin unreliable; no 4D theorem.

**Correction to Method 5:** its "4D marginal ⟹ regularity structures fail" conflates
**equilibrium** marginality (the beta-function dimension) with **parabolic** criticality
(the SPDE scaling). These are different; the failure must be re-derived from the parabolic
counting, not asserted from `[g]=0`. Re-check.
