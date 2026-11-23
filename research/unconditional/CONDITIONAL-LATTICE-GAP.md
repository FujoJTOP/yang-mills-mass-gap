# CONDITIONAL-LATTICE-GAP — blocked SLC ⇒ fixed-β lattice mass gap

**Verdict.** A precise, self-contained conditional theorem for 4D pure SU(N) Yang–Mills. Fix β in the crossover window (SU(2) β≈2.5; SU(3) β≈5.7). If **(H1)** the blocked marginal is uniformly strongly log-concave with constant λ; **(H2)** the UV fluctuation block has a uniform local Poincaré gap Λ_UV; **(H3)** the factor-2 blocking is faithful — **then** the Wilson measure obeys a volume-uniform local Poincaré inequality with constant 1/m, m=λΛ_UV/(λ+Λ_UV), equivalently the transfer matrix has gap Δ_Λ≥m. This is a **fixed-β lattice** gap, **not** a continuum mass gap. It is a sufficient route, not a restatement. Every hypothesis is open; the decisive one (H1) is currently disfavoured by the softest mode. No step below is fabricated — each is a standard implication, stated with its constant.

## Setup
Λ=Λ_{L,a}: torus T⁴_L of side L, spacing a, N_s=L/a sites/dir, gauge group G=SU(N). Links U∈G^E, Wilson measure μ_Λ∝e^{−S_W(U)}dU (Haar), S_W=βΣ_p(1−(1/N)Re tr U_p). Let P be a gauge-covariant factor-2 low-pass; coarse field U′=PU; marginal ν_Λ=P_*μ_Λ; S_eff(U′)=−log(dν_Λ/dU′) on coarse Haar. Exact marginal-curvature identity (W4; MARGINAL/VAR-FINAL):
$$G''(U')[v]=\mathbb E_\mu[P^{\!\top}HP\mid U']-\operatorname{Var}_\mu(P^{\!\top}\nabla S\mid U'),\quad H=\operatorname{Hess}S_W,$$
the coarse Hessian = W4's positive object **minus** the entropic variance penalty. Write μ(·|U′) for the regular conditional law of the fluctuation fiber, so μ_Λ=∫μ(·|U′)ν_Λ(dU′). The tangent space splits **orthogonally** as coarse⊕fiber with gradients ∇_∥, ∇_⊥, hence |∇F|²=|∇_∥F|²+|∇_⊥F|². On G^{E′}, Ric≥ρ_G g, ρ_G>0 (SU(2): ½g; SU(3): ¾g — W6); this only helps.

## Hypotheses (uniform in L, a, block scale)
**(H1) Uniform blocked SLC.** ∃λ>0: for every coarse field U′ and every coarse tangent v⊥ gauge orbit, ⟨G″(U′)v,v⟩≥λ‖v‖², i.e. Hess S_eff⪰λg on coarse links.

**(H2) UV fluctuation gap.** ∃Λ_UV>0: for every U′ and every local f on the block,
$$\operatorname{Var}_{\mu(\cdot|U')}(f)\le\tfrac1{\Lambda_{\rm UV}}\mathbb E_{\mu(\cdot|U')}\bigl[\textstyle\sum|\nabla_\perp f|^2\bigr].$$

**(H3) Faithful blocking.** μ_Λ=∫μ(·|U′)ν_Λ(dU′) is an **exact** regular-conditional decomposition; every gauge-invariant local F is (U′,ξ)-measurable and, for φ(U′)=E_μ[F|U′], the gradient contraction E_ν[|∇φ|²]≤E_μ[|∇_∥F|²] holds (no IR content hides in the fiber ξ).

## Theorem
Under (H1)–(H3), with $$\boxed{m=\lambda\Lambda_{\rm UV}/(\lambda+\Lambda_{\rm UV})>0},$$ for every local gauge-invariant F and every Λ,
$$\operatorname{Var}_{\mu_\Lambda}(F)\le\tfrac1m\sum_x\mu_\Lambda(|\nabla_xF|^2).$$
Hence the OS transfer matrix T=e^{−aH} (0≤T≤1) has Δ_Λ=−(1/a)lnλ₁(T)≥m, uniformly in L,a; by OS reconstruction this is exponential clustering at rate m.

## Proof (chain and constants)
1. **SLC ⇒ coarse Poincaré (1/λ).** Bakry–Émery/Γ₂: Ric+Hess S_eff⪰(λ+ρ_G)g ⇒ Var_ν(F′)≤(1/(λ+ρ_G))∫|∇F′|²dν ≤(1/λ)∫|∇F′|²dν. (Brascamp–Lieb is the flat special case.)
2. **Two-scale law of total variance (exact).** $\operatorname{Var}_\mu(F)=\mathbb E_\nu\operatorname{Var}_{\mu(\cdot|U')}(F)+\operatorname{Var}_\nu(\mathbb E_\mu[F|U'])$. Fiber term ≤(1/Λ_UV)E_μ[Σ|∇_⊥F|²] by (H2). Coarse term ≤(1/λ)E_μ[Σ|∇_∥F|²] by Step 1 with (H3)'s contraction. Since 1/m=1/λ+1/Λ_UV, summing gives the boxed constant.
3. **Poincaré ⇒ gap.** For a Markov generator gap = inf_f Dirichlet/Var; a volume-uniform **local** Poincaré constant 1/m is equivalent to Δ_Λ≥m (finite-range localization; P4 §4; Helffer–Sjöstrand; SZZ).

## What discharges each hypothesis
- **H1 ⟸ the G″>0 sign test.** Need a uniform-in-U′ lower bound on E[PᵀHP|U′] and upper bound on Var(Pᵀ∇S|U′) with positive gap, via the exact G″ identity. W4 samples ~10 L=4 configs (SU(2) β=2.5 blocked min +0.81β; necessary, not uniform); VAR-FINAL's softest k→0 mode gives +0.043−0.208≈**−0.165<0**, so H1 is **currently disfavoured**. The binary experiment is the sign of G″(0) by exact constrained sampling (W6 foothold).
- **H2 ⟸ a provable UV block gap.** A fixed-scale (b≈2a), finite-range block Poincaré gap independent of U′,L,a. One block suffices for SU(2) at β=2.5; SU(3) deep-SC retains 12% negative modes (W4), so H2 is open there.
- **H3 ⟸ blocking faithfulness.** Exactness of the decomposition plus sufficiency of (U′,ξ); a structural property of P (checkable for invertible covariant blocks), not implied by H1–H2.

## Honest status
Sufficient, not equivalent: Poincaré does not imply SLC, so proving (H1)–(H3) would prove a genuine gap, not assume it. **But it does not evade the Clay problem.** H1 at the softest coarse (k→0) mode is the mass gap in curvature language, and positivity at b≈2a says nothing at k=0 without an RG step (W6 iii.b); H2 is unproved for SU(3); both are open. The conclusion is a fixed-β **lattice** gap; the continuum a→0 along the RG trajectory (IN-9) needs separate uniform control and is not addressed. **Nothing here is proved.** Value: the hard problem is isolated into three concrete, individually testable obligations — none of them "prove the mass gap."

*Refs.* Osterwalder–Seiler, Ann. Phys. 110 (1978) 440; Bakry–Émery / Brascamp–Lieb; Helffer–Sjöstrand; Shen–Zhu–Zhu arXiv:2204.12737; Nachtergaele–Sims, CMP 265 (2006) 119.
