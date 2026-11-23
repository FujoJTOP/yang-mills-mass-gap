# W4 — Blocked-positivity numerical test

**Verdict: the positivity seed is not dead — but the blocked Hessian is positive only near the crossover, only in a factor-2 low-pass band, and this is not a proof of SLC.**

**Method.** SU(2)/SU(3) Wilson action, heat-bath thermalized. For each config: intrinsic Hessian H=d²S on all links (left-geodesic FD; the antisymmetric part is per-link only — cross-link max 5e-11 — so symmetrisation is exact); project out the gauge orbit (Q⊥); "blocked" Hessian = H restricted to the factor-2 low-pass (|k_μ|≤π/2, Nyquist plane dropped) ∩ Q⊥. Checks passed: identity-config transverse Hessian positive-semidefinite (min≈0, free-gluon max=3β); gauge·gradient≈1e-15; analytic vs FD gradient≈1e-7.

**Numbers** (eigenvalues in units of β; 7 SU(2), 2+2 SU(3) configs, L=4):

| case | per-link neg | full Q⊥ trans neg (min) | blocked neg (min) | blocked dim |
|---|---|---|---|---|
| SU(2) β=2.5 (near β_c) | ~0% | 7.1% (−0.82) | 0% (**+0.81**) | 204 |
| SU(3) β=2.5 (deep SC, 04 regime) | 48% | — | 12% (**−0.14**) | 544 |
| SU(3) β=5.7 (near β_c) | 0% | — | 0% (**+0.44**) | 544 |

One APE smoothing step at SU(2) β=2.5 cuts full transverse negativity 7.1%→~1% (min −0.82→−0.30β). Only ~8% of the negative-eigenvalue density sits in the blocked band.

**Findings.**
- The 49% per-link negativity is an **SU(3) strong-coupling** effect; SU(2) at its crossover has ~0% per-link negative (04's number does not describe the SU(2) crossover).
- The negative modes are **not pure gauge**: the gauge-projected transverse Hessian is genuinely indefinite (7% negative at the SU(2) crossover).
- But they are largely **high-momentum**: one factor-2 block leaves a positive gauge-invariant Hessian (min +0.8β) at the crossover (SU(2) β=2.5; SU(3) β=5.7, min +0.44β). At SU(3) deep strong coupling one block is **insufficient** (12% still negative).

**Boundary.** Pointwise over ~10 L=4 configs, not a uniform SLC bound; blocked band is tiny (204/2304 modes); positivity of a sampled restriction is necessary, not sufficient, for SLC. Net: the W2-1 assumption "blocking must re-import Savvidy negativity" is **falsified at the physical crossover**, but blocked-scale SLC remains unproved.
