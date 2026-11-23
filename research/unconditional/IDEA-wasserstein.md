# IDEA — Wasserstein–Otto / CD(m,∞) for the 4D YM Gibbs measure

**Verdict: REFUTED as an independent sufficient condition for the microscopic Gibbs measure — displacement convexity is FALSE for m>0 at physical coupling. The only surviving slice (collective/marginal) does not give a new route: it reduces exactly to the induced gluon mass, i.e. to the gap.**

**1. Precise statement — Gibbs (not fine-blocked).** Let M=G^E/G_site with the L² bi-invariant metric and μ_β=Z⁻¹e^{−S_W}dVol, S_W=βΣ_p(1−Re Tr U_p/N). S_W is gauge-invariant ⇒ its orbit Hessian vanishes and the BE condition descends to the quotient [Sturm 2006]. CD(m,∞) [Otto 2001; Sturm 2006; Lott–Villani 2009] is the pointwise bound
  Ric_M(v,v)+Hess S_W(v,v) ≥ m|v|²,  ∀x∈M, ∀v∈T_xM,  (★)
⇔ m-convexity of Ent(·|μ_β) along W₂-geodesics ⇔ W₂(P_tμ,P_tν)≤e^{−mt}W₂ [von Renesse–Sturm 2005] ⇒ LSI 2/m, Poincaré 1/m, gap≥m. Since S_W is bounded with e^{−S_W}>0 everywhere, μ_β has full support, so (★) is pointwise on all of M — measure-zero carve-outs are impossible.

**2. Full measure: REFUTED (proven false).** `witness`. Take a full plaquette with holonomy U_p=−I∈SU(2) (one link −I, the other three =1); it is horizontal (changes field strength). Restricted to the maximal torus the potential is β(1−cos θ)-type; every plaquette meeting that link contributes Hess βcos θ, net ≤−cβ with c=O(1) at θ=π, while Ric_M=O(1) (SU(2) ½). Hence (★) forces m ≤ c′−cβ<0 whenever β>β*(N)=O(1). Numerically at SU(2) β=2.5 the gauge-projected transverse Hessian reaches −0.82β [W4], the BE tensor ≈ Ric−2.8≈−2.3 [W6-BE]. So CD(m,∞) fails for every m>0; the Savvidy modes are genuine (not pure gauge) and 90–95% orthogonal to any collective band.

**3. Collective slice: only averaged, blocked-scale positivity — NOT CD.** On |k_μ|≤π/2∩Q⊥ the ensemble-averaged coarse Hessian E[PᵀHP] is positive definite (min +0.043 SU(2) β=2.5; +0.44β SU(3) β=5.7) [W4,P3]. But CD is pointwise, not averaged; and the exact marginal curvature is G″=E[PᵀHP|A′]−Var(Pᵀ∇S|A′), whose softest proxy is −0.165<0 [VAR-FINAL]. Averaged subspace positivity is strictly weaker than displacement convexity.

**4. Exact reduction to the gap.** Free curvature on momentum k is λ(k)=(4/a²)sin²(ka/2)→0 as k→0 [P1]. On the k→0 transverse mode (★) reads Ric+interaction-curvature≥m, and the interaction curvature is exactly m_eff²(0), the induced gluon mass (P1's c=m²; P4's L*). Therefore
  CD(m,∞) on the k→0 transverse mode  ⇔  m_eff²(0)≥m>0 uniformly  ⇔  volume-uniform mass gap [W6-induced-mass].
Full-measure CD is false (§2); the only true slice is the marginal one, where CD ⇔ induced mass ⇔ the target. CD is also ≠ RG-contraction (P3 §3).

**Sharpest honest line.** Proven FALSE on the Gibbs measure for m>0 (Savvidy/holonomy witness, β>β*≈1); on the collective-marginal slice it is *equivalent* to the induced-mass/gap. No independent Wasserstein route; the open content is again m_eff²(0)>0.
