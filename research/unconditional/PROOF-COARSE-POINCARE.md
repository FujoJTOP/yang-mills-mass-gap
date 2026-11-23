# PROOF-COARSE-POINCARE — attacking the coarse Poincaré (H1) directly: effective action, induced mass, strong coupling

> **Verdict.** The coarse Poincaré `Var_ν(φ)≤(1/λ)∫|∇φ|²dν` for block-averaged `φ` is **not proved at β=2.5**, and the three direct routes fail for *independent, nameable* reasons. **(a)** The effective action is **Savvidy-indefinite**, not provably convex: its exact Hessian `G″=E_μ[PᵀHP|U′]−Var_μ(Pᵀ∇S|U′)` has a bare part still negative on the coarse subspace (min ≈−0.3) and a one-loop (Schur) form **undefined** because the fiber Hessian is indefinite. **(b)** The block gap is a **fast relaxation rate**, not the slow stiffness; adiabatic elimination controls the memory kernel, never `Hess S_eff`, and H2's variance bound runs the **wrong way** (`1/Λ_UV≳e^{5N_p}`). **(c)** At **β<1/48 the coarse Poincaré is provable** (push SZZ's LSI through the Lipschitz block map) — a genuine anchor that does not transport to 2.5. The residual is `G″≥λ` uniformly ⟺ `m_eff²(k→0)≥λ` — H1, the gap in curvature language.

## The object, exactly
`ν=P_*μ`, `S_eff=−log(dν/dU′)`. For block-averaged `φ`, `∫|∇φ|²dν=∫|∇_∥(φ∘P)|²dμ` (`∇_⊥` annihilates `φ`). Need `Hess S_eff⪰λg`, uniform in `(L,a)`. The exact identity `G″(U′)[v]=E_μ[PᵀHP|U′][v,v]−Var_μ(⟨Pᵀ∇S,v⟩|U′)` is the whole content.

## (a) Effective action — Savvidy-indefinite at every order
**(i) Pointwise:** the coarse-restricted bare Hessian `PᵀHP` has min ≈−0.3 at β=2.5 (MARGINAL) — the Savvidy negative mode survives blocking, being the k→0 mode of a slow background. **(ii) Entropic:** `Var≥0` subtracts, so `G″≤E[PᵀHP]`, and the action margin is saturated by force variance (W7: bare `E[H]≈99%` cancelled) — no convexity margin. **(iii) One-loop:** the Schur complement needs `H_ff⪰0`, but 90–97% of negative BE directions lie in the fiber (W6-BE; MARGINAL); `H_ff` is indefinite, so the Gaussian effective Hessian is **undefined**. **(iv) Prékopa:** marginalization preserves log-concavity only for jointly log-concave `S_W`, which fails. So the effective action is not convex; its Hessian is `G″`, and `G″⪰λ` is H1.

## (b) Induced mass — fast gap ≠ slow stiffness; wrong direction
Mori–Zwanzig elimination yields a slow drift plus a memory kernel; `Λ_UV` sets the **memory decay rate**, not the static curvature. Independence is exact: `ξ~N(0,1/Λ_UV)` independent of a non-log-concave slow `X` gives fiber gap `Λ_UV→∞`, marginal gap `0`. Quantitatively worse: H2 gives only `Var_μ(⟨Pᵀ∇S,v⟩|U′)≤(1/Λ_UV)E|∇_⊥·|²`, yet `1/Λ_UV=e^{2βN_p}/C_G≳e^{5N_p}` while the variance is `O(0.2)` — vacuous by a factor `~e^{5N_p}` (`N_p` a block's plaquettes, O(10²)). H2 supplies only a *tiny lower* bound on the fiber gap; to force `Var` below `E[H]−λ` one needs an *upper* bound on the LSI constant `1/Λ_UV` (a usable gap bound) or a direct variance estimate — neither exists. The induced mass `m_eff²(0)=m_0²+Π(0)`, `m_0²=0`, is `Π(0)`, the Schwinger self-energy, signed by the **dynamical input** and fed by the fiber's negative Savvidy directions. RPA reproduces the gap equation; its nonzero branch *is* the gap.

## (c) Strong coupling — provable, does not transport
At β<1/48 SZZ give `Var_μ(F)≤(1/λ_μ)∫|∇F|²dμ` uniformly in `V` (CMP 400 (2023), arXiv:2204.12737). With `Lip(P)=C_blk=O(1)` (averaging is a contraction, volume-independent), `F=φ∘P` gives `Var_ν(φ)≤(C_blk²/λ_μ)∫|∇φ|²dν`. **Lemma SC: at β<1/48 the coarse Poincaré holds with `λ=C_blk^{-2}λ_μ`, uniform in `(L,a)`.** P3's residue `f₁²` is the spectral atom, not the coarse marginal, so the anchor is SZZ's. No leverage at 2.5: the pulled-back LSI is the full-measure gap unproved there, and β=2.5≫1/48 lies outside every high-temperature regime.

## Exact missing estimate
```
∃λ=λ(β,N)>0: ∀U′, ∀ coarse v⊥orbit,
  E_μ[PᵀHP|U′][v,v]−Var_μ(⟨Pᵀ∇S,v⟩|U′) ≥ λ‖v‖²,  uniform in (L,a);
⟺ m_eff²(k→0)≥λ  ⟺  Var_ν(φ)≤(1/λ)∫|∇φ|²dν.
```
Every banked piece lies on the other side: H2/H2-global bound `∇_⊥`, in the kernel of the coarse form; block positivity at `b≈2a` is silent at k=0; the softest `G″` sign is undetermined (W7 +0.27 at k=π/2, action margin ≈0; VAR-FINAL −0.165), and a *uniform* sign **is** G2.

**Sharpest line.** No banked lemma yields the coarse Poincaré at β=2.5: the effective action is Savvidy-indefinite so not convex, the fiber gap is a fast relaxation rate whose only variance bound runs the wrong way by `e^{5N_p}`, and the one regime where it *is* provable (β<1/48, pushing SZZ through a 1-Lipschitz block map) is exactly where the full LSI already holds. Direct attack confirms PROOF-IR-H1: H1 is the bare gap, and no coarse-grained Dirichlet form supplies it.

*Refs.* W4; CONDITIONAL-LATTICE-GAP (H1, exact `G″`); MARGINAL (`PᵀHP` min −0.3); W7-Gpp0 (`G″` +0.27, action margin ≈0); VAR-FINAL (−0.165); W6-BE-at-2.5 (fiber indefinite); PROVE-H2-uvgap + H2-b-uniform (`Λ_UV=e^{−2βN_p}`); PROOF-DIRICHLET-COMPARISON (fiber form kernel = coarse functions); PROOF-IR-H1 (orthogonal split, H1=gap); PROOF-GRADIENT-FLOW (LSI⇒pushforward); PROVE-H3 (mixed-Hessian score term); W6-induced-mass + ATOM-FRG-mass + GAP-G2-mass (Schwinger `Π(0)`, sign dynamical); Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737 (β<1/48); Mori–Zwanzig; Bakry–Émery; Prékopa; Holley–Stroock; Ledoux (LSI pushforward); Cornwall, PRD **26** (1982) 1453.
