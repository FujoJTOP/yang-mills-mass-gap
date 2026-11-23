# REDTEAM-FORGE — adversarial audit of every forged tool (2026-11-23)

> **Verdict: no candidate satisfies all five hard constraints.** The library is dominated by honest no-gos (which yield no lower bound) plus one genuine but regime-locked residue theorem. Constraint key: **(1)** nonperturbative **(2)** sign on Savvidy cone **(3)** non-circular **(4)** non-SUSY & fails for free field **(5)** uniform lower bound.

| Candidate | Breaks | Exact flaw |
|---|---|---|
| INVENT-GAUGE-INVARIANT-MASS | 3,5 | `m_GI`=edge of the `trF²` correlator is a *definition* of the 0⁺⁺ mass; "`m_GI>0`" restates the gap. Bounds only the 0⁺⁺ channel (needs the unproven "lightest state is 0⁺⁺"). Free Maxwell passes (`=0`) — content is the gap renamed. |
| INVENT-TRANSPORT-INEQ | 3 | (CP) controls the *energy* (Kato), not the *residue*: the operator can decouple (residue→0) at a fixed below-threshold pole — PROOF-DISSOLUTION's own wall. `f₁²=Z_O²f₁²(β_eff)` is not implied by (CP); and `Z_O` is one-loop, inapplicable at β<1/48 (g→∞) and at 2.5. Circular. |
| INVENT-CONVEXIFICATION | 1,5 | Only honest inequality (b): `Var ≤ (b/a)C E`, with `b/a=e^{c·osc|F|²}=e^{cV}` — volume-exponential, not uniform. (a),(c) fail sign/lift. Periodic-concave no-go is real; yields no tool. |
| INVENT-SECOND-CONE | 4 | Karlin cones *nested* `K₀⊂K₁⊂…`, so the strictly-positive free field (`f₁²=0`) lies interior to every level — vacuous. N=2 only (no Sturm flag for N≥3). |
| INVENT-VARIANCE-INEQ | 3 | `G″=S−D` exact, but `K=Π₀HessVΠ₀` is indefinite at 2.5 ⇒ `⟨K⁻¹c,c⟩` sign-meaningless and `D<S` is verbatim `G″>0`=the gap. `D/S≈2%` is numeric. |
| INVENT-NEW-INDEX | 4,5 | Rigidity no-go: any locally-constant index = free value = 0. Correct, but no lower bound. |
| INVENT-NEW-IDENTITY | 4 | "Pole below cut" is trivial (`m₁<2m₁`); at finite L the *free* field also has an isolated zero+residue. The distinction (residue→0 free vs `≥r₀` YM) is *only* the `L→∞` limit = the theorem. Identity exact; acid test deferred to the limit. |
| INVENT-NEW-MONOTONE | 3,5 | The monotone that holds (Fisher/de Bruijn) is relative entropy to μ — circular; the gap-certifying `Γ₂` is sign-indefinite at 2.5. |
| FORGE5-DUALITY | 1,2,3,5 | No exact 4D non-SUSY duality exists; exact duality = isomorphism ⇒ preserves the spectrum (circular); wrong channel (Proca `√σ` ≠ 0⁺⁺). #4 only "vacuous." Self-adjudicated fail. |
| PROOF-JACOBIAN-SAVES | 2,5 | `W_geo≥0` is a theorem but at pure-Haar/β<1/48; at 2.5 swamped ≈10× *pointwise*. "+0.51" is L=4,6 numerics at k=π/2,π/3 — never k→0, where the gap lives. |
| PROOF-P3-UNCONDITIONAL | 2,5 | Genuine `f₁²≥r₀>0`, V,L-uniform — but β<β₁=O(1/N)→0, **fixed a** (not continuum), not on Savvidy. Its "grading yields isolation (I2)" contradicts PROOF-ISOLATION ("two-exponential bound not closed"). Strong-coupling gap itself is classical OS78. |
| PROOF-GMT-CURRENTS | 4 | Every GMT bound is a metric filling volume >0 already at σ=0 ⇒ survives the free field; thin-string regime ≠ fat glueball. |
| PROOF-MORSE-FLOER | 4 | Hess at vacuum `=−Δ`=free Maxwell for all g; Savvidy negatives lie *off* the critical locus. Coupling-free topology vs coupling-dependent gap. |
| PROOF-OPTIMAL-CONTROL-HJB | 3,5 | A lower bound needs the inequality uniform over controls = LSI/H1; the optimizer is ψ₁, the very eigenvector sought — circular. Quadratic-cost HJB `=+∞` at 2.5. |
| PROOF-QUANTUM-GROUP-DEFORM | 4 | Positivity is aggregate (`TV=|WRT|²`), truncation-bought; q→1 kills it (`[k+N]_q=0`) and the Weingarten sign returns; only d≤3 topological models — wrong dimension. Both horns fail. |
| PROOF-RKHS-GAUSSIAN-PROCESS | 2,5 | Mercer eigenvalues are positive functionals of ν ⇒ `f₁²≤λ_k/w_k` — an **upper** bound (wrong direction). BBP spike needs bounding the continuum = the gap. |
| PROOF-ROUGH-PATH-SIGNATURE | 2,4 | The signature records increments, *discards the constant k→0 mode* — exactly the gap's mode. The Gaussian free field is the prototype rough path; every bound holds at g=0. |
| PROOF-SAVVIDY-EXACT-SIGN | 4 | Correct no-go: `λ_min=−g|B|`, `N_neg∝(gB)²V`, weight `∝∫|F|³` all grow with the field; only `Tr H=4Tr(−D²)≥0` ⇒ `Tr(H₋)≤Tr(H₊)`, saturated (0.94), not a gap. No positivity. |
| PROOF-STOCHASTIC-LOCALIZATION | 3,5 | Convexifies only at t>2.8 (`μ_t→δ`, trivially log-concave); the gap is the t→0 slice, whose energy `=G″=S−D=H1`; backward transfer costs `e^{ct/a²}`. |
| PROOF-SYMPLECTIC-REDUCTION | 4 | DH pushforward is absolutely continuous (atom-free) ⇒ structurally cannot carry `f₁²`; the free Gaussian case is exactly localizable and gives 0 everywhere. |
| PROOF-THERMODYNAMIC-FORMALISM | 3,4 | Pressure `=log λ₀`; the gap is `log(λ₀/λ₁)` and λ₁ has no bulk pressure. Analyticity domain = Fisher zero-free = β<1/48 (circular with P3); free field: analytic p, Δ=0. |
| PROOF-TOPOLOGICAL-DATA | 4 | `σ(H)` 0-dim ⇒ all bars infinite (vacuous); **and** at finite L the *free* field has `E₁(L)=ω(2π/L)>0` — a positive "bar" with `f₁²=0`. A uniform-in-L persistence bound would falsely certify the free field. |

**Genuine survivors (exact content only).** (i) PROOF-P3/P3-UNCONDITIONAL: a real, **V,L-uniform residue bound** `f₁²≥r₀(β)>0` for β<β₁=O(1/N) at *fixed* a — a regime-locked residue theorem, not a continuum or Savvidy-cone gap. (ii) The no-go theorems are genuine and load-bearing — Schur-convex master barrier (no protected g-continuous inequality forces `f₁²>0`; the free field is the witness), periodic-concave, invariant-inertia (Sylvester) and index-rigidity. They **close routes**; none produces a lower bound, so none is a *tool*.

**Cross-cutting flaw (defeats several claimed "passes").** The acid test `f₁²(free)=0` is true only in the `L→∞, a→0` limit; at any finite regulator the free lattice field has a positive lowest mode and nonzero residue. Tools claiming to pass (NEW-IDENTITY, TOPOLOGICAL-DATA) do so by evaluating the free field *in the very limit that is the theorem* — a disguised circularity (#3). Likewise "H1 = the gap" is an iff only as a slogan; the valid content is "H1 ⇒ gap, H1 unreduced," a rename, not a closure. Every "sharpened criterion / exact identity / same wall" verdict is honest self-description of a **failure**, not a survivor.
