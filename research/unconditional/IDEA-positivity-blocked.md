# IDEA — Blocked-measure positivity: a conditional SLC ⇒ Poincaré ⇒ gap theorem

**Verdict.** *If* the blocked (marginal) Yang–Mills measure is **uniformly strongly log-concave**, then SLC ⇒ Poincaré ⇒ spectral gap ⇒ a rigorous **lattice mass gap at β≈2.5**, far beyond SZZ's β<1/48. The implication chain is standard and sound. The hypothesis is open, is **strictly stronger than the conclusion** (sufficient, not equivalent), and the one decisive term in it currently leans negative (VAR-FINAL). A route, not a proof.

## Setup
Let Λ=Λ_{L,a}, links U∈G^E, μ_Λ∝e^{−S_W}dU. Block by a gauge-covariant factor-2 low-pass P, marginal ν_Λ=P_*μ_Λ on coarse links U′, effective action S_eff=−log dν. For a marginal, the exact identity is
$$G''(U')[v]=\operatorname{Hess}S_{\rm eff}=\mathbb E_\mu[P^{\!\top}HP\mid U']-\operatorname{Var}_\mu(P^{\!\top}\nabla S\mid U'),$$
the curvature of the coarse log-density. This is W4's positive object **minus** the entropic variance penalty (MARGINAL/VAR-FINAL). Compactness is free help: on coarse links the group Ricci (SU(2)=½g, SU(3)=¾g; W6) adds a positive constant.

## Conditional theorem
Assume, uniformly in L,a and block scale:
- **(H1) SLC.** ∃λ>0 with G″(U′)≥λ·Id for **every** coarse field U′ and every coarse tangent vector (the blocked Hessian is uniformly positive).
- **(H2) UV regularity.** The conditional law of the fluctuation block given U′ has a local Poincaré/spectral gap ≥Λ_UV>0 (a finite-range, block-level statement).
- **(H3) Faithful blocking.** Every gauge-invariant local observable is a function of U′ plus a fluctuation H2 controls (the coarse field carries all IR gauge-invariant content).

**Then** ∃m=λΛ_UV/(λ+Λ_UV)>0, uniform in L (and a along the block), with
$$\operatorname{Var}_{\mu_\Lambda}(F)\le \tfrac1m\sum_{x}\mu_\Lambda(|\nabla_xF|^2)$$
for all local gauge-invariant F; equivalently Δ_Λ≥m: a positive β=2.5 **lattice** mass gap, volume-uniform.

## Implication chain
1. **SLC ⇒ Poincaré.** H1 ⇒ Brascamp–Lieb/Bakry–Émery on ν: Var_ν(F′)≤(1/λ)ν(|∇_{U′}F′|²) — coarse constant 1/λ.
2. **Two-scale (law of total variance).** Var_μ(F)=E_ν Var_{μ(·|U′)}(F)+Var_ν(E_μ[F|U′]). First term ≤(1/Λ_UV)E[Σ|∇_⊥F|²] by H2; second ≤(1/λ)E[Σ|∇F|²] by H1 + Jensen on conditional expectation. Sum ⇒ fine local Poincaré, constant 1/m.
3. **Poincaré ⇒ gap.** A volume-uniform local Poincaré constant is the transfer-matrix gap (OS / Helffer–Sjöstrand; P4 §4; SZZ); via OS reconstruction it is exponential clustering at rate m. ∎

## What G″(0)>0 buys — and its obligations
- **Buys:** the *coarse/IR* Dirichlet-form factor 1/λ, the uniform constant exactly where the free propagator is massless (P1, λ_min→0). It replaces a global spectral gap by a **pointwise curvature sign** — binary and computable (W6 foothold).
- **(a)** Uniform lower bound on E[PᵀHP|U′] for **all** U′; W4 is pointwise over ~10 L=4 configs — not uniform.
- **(b), the wall:** Uniform upper bound on Var(Pᵀ∇S|U′). Variance only *subtracts*, so ensemble-mean positivity (+0.043) is compatible with G″<0. VAR-FINAL's softest k→0 mode gives **+0.043−0.208≈−0.165**: the decisive mode fails.
- **(c)** H2: the fluctuation block's own uniform Poincaré at β=2.5 (Savvidy/dilute-gas at scale a; one block suffices for SU(2), not deep-SC SU(3)).
- **(d)** H3 and persistence of positivity through the RG step to k→0 (W6 iii.b: positivity at b~2a says nothing at k=0).

## Honest status
SLC is strictly stronger than a gap (Poincaré does **not** imply SLC), so this is a genuine sufficient route, not a restatement. The conclusion is a fixed-β **lattice** gap; the continuum still needs a→0 along the RT (IN-9). Current numerics on the softest mode **disfavour** H1; the route survives only if exact constrained-sampling flips the G″ sign. Nothing here is proved.
