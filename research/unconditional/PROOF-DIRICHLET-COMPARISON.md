# PROOF-DIRICHLET-COMPARISON — can the YM Dirichlet form inherit a gap from a gapped reference?

> 2026-11-23 · Deuschel–Stroock route to the 4D fixed-β lattice gap; *comparison*, not LSI/BE construction.

> **Verdict. NO gap by comparison — and not at Savvidy.** The YM Dirichlet form `E_YM(f)=μ_Λ(|∇f|²)` is positive for every field and every test function; Savvidy's negative Hessian modes never make it indefinite (a Dirichlet form is a quadratic form in `f`, not in `U`). The comparison fails structurally: every form on `L²(μ_Λ)` that is *actually* gapped is `E_YM` itself or already assumes the coarse Poincaré. The proven gapped reference (H2 fiber form `E_⊥`, constant `Λ_UV`) has every coarse function `φ(U′)` in its kernel, so its gap on `L²(μ_Λ)` is **zero**; a massive-Gaussian field reference is a *different measure* whose transfer constant obeys `c≲e^{−2βN_□}`, non-uniform in `(L,a)`. The exact missing inequality is the coarse Poincaré (H1) — comparison renames the gap by zero distance, exactly as PROOF-IR-H1.

## 1. Deuschel–Stroock, precisely
Let μ be a probability measure and `E₁,E₂` symmetric Markovian (Dirichlet) forms on `L²(μ)`, **same measure**, `D(E₁)⊂D(E₂)`. If `E₁(f,f) ≥ c E₂(f,f)` for all `f` and `E₂` has Poincaré constant `C₂` (`Var_μ f ≤ C₂E₂(f,f)`), then `Var_μ f ≤ (C₂/c)E₁(f,f)`: the dominating form inherits the smaller one's gap. The load-bearing hypothesis is the **shared measure**. For forms of `μ_A≠μ_B` one needs two-sided bounded densities; but even then the transfer gives a gap only for a measure in that equivalence class, never for `μ_Λ` unless `μ_Λ` itself is in it. (Holley–Stroock is the Gibbs-measure instance.)

## 2. The two references, and why each degenerates
Split `T(G^E)=coarse⊕fiber` orthogonally; then `E_YM=E_∥+E_⊥`, `E_∥=μ(|∇_∥f|²)`, `E_⊥=μ(|∇_⊥f|²)`.
- **Fiber / H2 (PROVE-H2, proved).** `E_⊥` is gapped *per block*, `Λ_UV=e^{−2βN_p}/C_G`, uniform in `U′,L,a` at `b=2a`. But every `f=φ(U′)` has `∇_⊥φ=0`, so `E_⊥(φ)=0` while `Var_{μ_Λ}(φ)>0`: on `L²(μ_Λ)` the gap is **0**. The trivial `E_YM≥E_⊥` transfers nothing. (That the block-local gap does not globalize is the *block coupling*, `C(β)>1` — explicitly **not** Savvidy, H2-GLOBALIZATION.)
- **Coarse (H1).** `E_∥` annihilates fiber-only `f`; its gap on `L²(μ_Λ)` is likewise 0. Its gap *on coarse functions* is `λ` — but that **is** H1, the target.
- **Massive Gaussian on the field.** Gapped (`m²` on mean-zero), but `E_G≤E_YM/c` forces `dμ_G≤(1/c)dμ_Λ`, i.e. `c ≤ inf (dμ_Λ/dμ_G) ≲ e^{−2βN_□}`, `N_□∼L⁴/a⁴`: no uniform `c`; in the continuum the measures are mutually singular (W2).

So the only form gapped on all of `L²(μ_Λ)` is `E_YM=E_∥⊕E_⊥` itself — comparison is either tautological or assumes H1. Along the `a→0` tower the fiber constant also degrades, `Λ_UV=e^{−2β_jN_p}` (PROOF-L3 §3).

## 3. Savvidy: where it bites, and why quotienting fails
`E_YM(f,f)=∫|∇f|²dμ≥0` always. Savvidy's `Hess S_W<0` defeats the **Γ₂/BE construction** (`Γ₂` carries `+Hess S_W(∇f,∇f)`), not the positivity of `E_YM`. In the Gaussian comparison it marks the large-field region where `e^{S_W}` is maximal (`S_W∈[0,2βN_□]`), i.e. where the constant `c` is worst — but the failure is *non-uniformity*, not sign. Quotienting cannot repair it: the unstable directions carry positive μ-mass, and marginalizing them injects a Jacobian reproducing exactly the `−Var(P^⊤∇S)` term that makes `G″` indefinite (GAP-G1: `G″≤E[H]` because `Var` subtracts).

## Exact missing inequality
A comparison gap needs a symmetric form `E_ref` on `L²(μ_Λ)` with (i) `E_ref≤E_YM/c`, and (ii) Poincaré on **all** of `L²(μ_Λ)`, coarse functions included. Condition (ii) forces `E_ref⊇E_∥` gapped, i.e.
$$\operatorname{Var}_\nu(\varphi)\le \tfrac1\lambda\int|\nabla\varphi|^2\,d\nu,\qquad \lambda>0\ \text{uniform in }(L,a),$$
which is **H1 = the mass gap in curvature language** (`H1 ⟺ G″⪰λ ⟺ m_eff²(0)≥λ`, PROOF-Gpp §5). The Gaussian variant additionally needs `‖dμ_Λ/dμ_G‖_∞<∞`, false in the continuum. Neither is supplied by comparison.

**Sharpest line.** There are no negative directions in the Dirichlet form to quotient; the comparison degenerates because each rigorously gapped reference (fiber H2, coarse H1) has the *opposite* scale in its kernel, and the unique measure gapped on both scales is `μ_Λ` itself — whose gap is the target. Comparison relocates H1 by zero distance; it does not bound it.

*Refs.* Deuschel–Stroock, *Large Deviations* (Academic Press, 1989), Ch. 5; Holley–Stroock 1987; Bakry–Gentil–Ledoux; PROVE-H2 (`Λ_UV`); PROOF-IR-H1 (orthogonal split; H1=gap); H2-GLOBALIZATION (block coupling, not Savvidy); CONDITIONAL-LATTICE-GAP (H1–H3); GAP-G1 (Savvidy cone; `G″≤E[H]`); PROOF-Gpp §5; PROOF-L3 §3; REPAIR-2SCALE; W2 (4D marginal/measure).
