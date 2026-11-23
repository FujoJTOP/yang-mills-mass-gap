# P2 — Porting the BBS d=4 marginal trap to gauge theory

**Verdict.** BBS's trap is a theorem about a *Gaussian-anchored* flow. Its gauge-covariant analogue states cleanly (W4's T1–T3), but one hypothesis has no gauge transplant: a **positive-definite (Savvidy-safe) large-field bound**. No current tool removes it: both bypasses fail, and the fallback is the marginal log-sum.

**(1) What BBS proves and requires.** On Λ⊂Z⁴, φ∈Rⁿ, reference = Gaussian covariance C (free massive propagator). Interaction K_j=λ_j∫φ⁴+R_j. The *trap* is 0<λ_{j+1}≤λ_j≤λ*, ‖R_j‖_{T_φ}≤O(λ_j²); the theorem: the RG map maps the trap into itself for all j≤J~log(1/a), uniformly in cutoff (BBS, LNM 2242 (2019); JSP 157 (2014), 159 (2015); AHP 16 (2015) for the non-hyperbolic flow). Ingredients:
- **running coupling** λ_{j+1}=λ_j−β₀λ_j²+O(λ_j³), λ_j≈1/(β₀j);
- **covariance** C positive-definite, translation-invariant;
- **large field** SF/LF split; LF bounded by positivity λ∫φ⁴≥0 plus Gaussian damping;
- **log resummation** Σλ_j²~log J diverges, so track the exact marginal trajectory; the deviation is forced at O(λ_j²)=O(j⁻²) (summable) and contracts at e^{−c}.

**(2) Gauge-covariant analogue.** On P(𝒜/𝒢) (CCHS, arXiv:2201.03487), K_j=g_j^{−2}∫Tr F²+R_j, norm on gauge-invariant dim>4 functionals; trap: 0<g_j≤g*, ‖R_j‖_X≤O(g_j²), g_{j+1}=g_j−β₀g_j³+O(g_j⁵), uniform in a. AF gives β₀>0 (Gross–Wilczek/Politzer 1973; 4-loop van Ritbergen–Vermaseren–Larin 1997), so the UV point is Gaussian—the *sign* matches BBS.

**(3) Breaks.**
(a) *Gauge.* Ward/BRST is benign: the invariant subspace is closed, the anomaly slot is empty (W3); Gribov/Singer (NPB 139 (1978); CMP 60 (1978)) kills only the gauge-*fixed* route. The break is the **orbit volume / zero mode**: C is degenerate along gauge directions, so the BBS norm does not descend to P(𝒜/𝒢).
(b) *Flowing marginal.* Tr F² is not a scalar monomial: it is kinetic+cubic+quartic, with wave-function renormalization Z_j. The trap must track (g_j,Z_j), not a single λ_j—the "non-Gaussian marginal coupling."

**(4) The irreducible hypothesis.** The **Savvidy-safe large-field bound**. BBS's LF estimate rests on a positive-definite quadratic form controlling fluctuations. For gauge, expanding about a background, the Hessian −D²+2 ad F₀ is indefinite (Savvidy 1977; Nielsen–Olesen 1978): constant chromomagnetic fields carry negative modes, so the large-field region is not damped and the trap leaks along Savvidy directions.

**(5) Removable?** No cleanly. Gauge-fix → FP determinant sign/zeros, non-global (Gribov). Covariant → no Gaussian anchor at d=4 (Hairer's degree condition fails; CCHS is d=3), and the fallback contraction *is* the marginal log-sum. The obstruction is dual to the wall: keep the Gaussian and lose covariance, or keep covariance and lose log resummation. The port stalls exactly at W4's point—G1 (UV) plausible, G2 (IR gap) unreached.
