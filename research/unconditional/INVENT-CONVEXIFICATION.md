# INVENT-CONVEXIFICATION — a log-concave surrogate, and why it cannot carry the YM gap back to μ

> **Verdict.** No convexification of the Wilson measure yields a volume-uniform gap at β=2.5. The obstruction is topological: μ's density is **periodic** (it lives on the compact group), and **no nonconstant periodic concave function exists** (superharmonic on a closed manifold ⟹ constant). Every construction thus fails where periodicity bites — at U_p≈−1, the Savvidy region. **(a)** Heat smoothing gives *semi*-concavity ∇²log u ⪰ −σ⁻²I, the **wrong sign** for Brascamp–Lieb, with unbounded deconvolution → no transfer inequality. **(b)** The reweight e^{−S_W+c|F|²} gives a **genuine** transfer Var_μ ≤ (b/a)C_μ̃ E_μ, but b/a = e^{c·osc|F|²} is exponential in N_□ and its reference gap is present at β=0 — **acid test #1 fires**. **(c)** A Prékopa lift is *forbidden* by Prékopa–Leindler: H_ff is indefinite. All respect PROOF-DIRICHLET: the Dirichlet form is positive throughout; only curvature/log-concavity dies.

## The periodic no-go — why "log-concave surrogate" is the wrong target
Flat chart X=log U∈g^E: μ ∝ e^{−S_W−Φ}dX with Φ=−log J convex (PROOF-JACOBIAN). e^{−Φ} is log-concave, yet the physical density is **periodic** — and a periodic concave function is constant. Hence −(S_W+Φ) is never concave on all of g^E; Savvidy's Hess ⪰ −1.1β (W6-BE) is exactly periodicity's cost at the wrap |x|>π/2. Any μ̃ that *is* log-concave therefore differs from μ in the wrapped region — precisely where |F| is largest.

## (a) Heat smoothing — wrong-sign semi-concavity
u = e^{−S_W}*γ_σ. **Tweedie's identity:** ∇²log u = σ⁻⁴(Cov_{ν_x}(Y) − σ²I), ν_x the posterior law of Y given Y+√σ²Z=x. As Cov ⪰ 0, **∇²log u ⪰ −σ⁻²I**: smoothing gives *semi*-concavity, never concavity (which needs Cov_{ν_x} ⪯ σ²I for all x, false once f spreads beyond σ). BL needs Hess ⪰ +λ; smoothing supplies a negative floor. Forcing λ>0 takes σ→∞, i.e. μ→Haar, gapped already at β=0 (test #1 fails). Transfer needs the symbol 1/ĝ_σ, unbounded: no bounded density ratio, no inequality. **Fails at: direction, and transfer.**

## (b) Multiplicative reweight — a true but vacuous transfer
Take μ̃ ∝ e^{−S_W+c|F|²}dU, so w=dμ/dμ̃=(Z̃/Z)e^{−c|F|²}. The transfer is a real theorem.

> **Construction + Transfer (genuine).** Let dμ = w dμ̃ with 0<a≤w≤b and μ̃ Poincaré with constant C_μ̃. Then
> **Var_μ(f) ≤ (b/a)·C_μ̃·E_μ(f,f),  E_μ(f,f)=∫|∇f|²dμ.**
> *Proof.* Var_μ f = inf_c∫(f−c)²w dμ̃ ≤ b·Var_μ̃ f ≤ bC_μ̃∫|∇f|²dμ̃ ≤ (b/a)C_μ̃ E_μ f. ∎

Since w=e^{−c|F|²}, b/a = e^{c·osc|F|²} — exponential in N_□, the PROOF-DIRICHLET loss verbatim. One c cannot serve both ends: the c fixing the high-k Savvidy floor (≈−1.1β) exceeds the O(β) free-field coefficient of |F|², making the k→0 log-Hessian *positive* — IR log-convexity, not log-concavity. **Fails #1 both ways; passes #2** (an honest inequality).

## (c) Prékopa lift — forbidden, not merely hard
e^{−S_eff(U′)} = ∫ e^{−S_W−Φ}dξ. Prékopa–Leindler: a marginal of a jointly log-concave density is log-concave. Contrapositive — G″=E[PᵀHP]−Var(Pᵀ∇S) is Savvidy-indefinite (PᵀHP min ≈ −0.3, MARGINAL) — so **no jointly-log-concave lift exists**; concretely the Schur complement is undefined because H_ff is indefinite (90–95% of negative directions in the fluctuation block). Inf-convolution does not help: the Moreau envelope of a non-convex −log μ is not convex. **Fails at: the lift itself.**

## Acid-test result
Only (b) passes test #2; all three fail test #1. With its best reference (Haar, C_G) the boxed inequality reads Var_μ ≤ e^{2βN_□}C_G E_μ — the massive-Gaussian loss of PROOF-DIRICHLET. **Where it dies:** log-concavity dies at the periodicity wrap |x|>π/2, and transfer dies because that is exactly where the density ratio blows up — one convexifies where μ is small and w is huge.

**Sharpest line.** The two-sided-density transfer is genuine but its constant *is* the gap: a uniformly gapped log-concave reference (Haar, C_G) gives ratio e^{2βN_□}; a ratio-bounded reference carries only the free-field gap (absent uniform-in-L). No construction threads both acid tests, and the residual lands on H1: a uniform C_μ equals m_eff²(k→0)⁻¹.

*Refs.* PROOF-JACOBIAN-SAVES (Φ convex, W_geo≥0); W6-BE-at-2.5 (Hess≈−1.1β, 92–97% fiber); MARGINAL (PᵀHP min −0.3); PROOF-DIRICHLET-COMPARISON (massive-Gaussian ratio e^{−2βN_□}); PROVE-H2-uvgap (H_ff indefinite); PROOF-IR-H1 + GAP-G2-mass (H1 = m_eff²(k→0)); W7-Gpp0 (G″≈+0.27≈1/Var); Brascamp–Lieb 1976; Prékopa 1971, Prékopa–Leindler; Holley–Stroock 1987; Tweedie 1947; Otto 2001; Ledoux, *Concentration of Measure* (2001).
