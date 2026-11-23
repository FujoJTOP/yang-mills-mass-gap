# T5 — Rigorous methods for signed / complex weights: any route past the Savvidy IR?

> G2 tool probe. Verdict: **none is rigorous at β≈2.47**; the only theorem-bearing tool (cluster) closes at β<β₀.

## Reframing

The **bare** measure at real β>0 is *positive*: `exp((β/N)Σ Re tr U_p)` is real and positive, and the Wilson character expansion has all-positive coefficients. Strong-coupling convergence comes from **small activity**, not sign cancellation. Our "signed gas" is thus an artifact of the **blocked effective action** (Savvidy-negative Hessian), not the path integral. Thimbles/CL target complex measures bare YM lacks; useful only if effective = physical, unproven. No theorem controls a "phase" of the SU(N) Wilson measure — at real β there is none; β-analyticity (positive coefficients) holds only for |β|<β₀.

## Per method

| Method | Rigorously proves | Death point |
|---|---|---|
| **Lefschetz thimble / PL** | Finite-dim PL is a theorem; path-integral version heuristic (Witten, arXiv:1001.2933) | ∞-dim interchange, gauge invariance (Gribov), residual Jacobian phase not bounded (arXiv:1205.3996, 1303.7204). |
| **Complex Langevin** | Sufficient criterion: holomorphic drift + localized distribution ⇒ boundary terms vanish (arXiv:0912.3360, 1101.3270, 2304.00563) | Sufficient, not theorem; singular drift / **wrong convergence** known. Real-β pure YM has nothing to complexify. |
| **Resurgence/trans-series** | Formal/asymptotic; resums 2D CP^{N−1} (arXiv:1210.2423) | ℝ⁴ renormalons defeat Borel–Laplace; bions only on ℝ³×S¹ (arXiv:1206.1890). No ℝ⁴ definition. |
| **Complex-β / signed cluster** | Kotecký–Preiss: small complex activity ⇒ absolute convergence, analyticity (CMP 103, 1986) | Only |β|<β₀; SZZ β<1/48 (arXiv:2202.10375) — crossover 15–120× beyond. |
| **Signed gas / cancellation** | Strong-coupling surface rep. non-negative; SU(2) κ=1 (arXiv:2307.06790) | No mechanism forces κ_c<1 (BL-6); cancellation uncontrolled without small activity. |
| **SU(N) Wilson phase** | Only positivity ⇒ strong-coupling gap (Osterwalder–Seiler, Ann. Phys. 110, 440; Seiler LNP 159) | Blocked effective weight non-positive; no theorem controls it. |

## Verdict

None crosses the wall: provable results live at β<β₀; thimbles/CL are theorem-free heuristics, resurgence is formal. **Needs-X:** extend complex-activity convergence to Re β~2.47, or prove the thimble/CL representation of the *blocked* theory equals the physical measure with controlled phase.
