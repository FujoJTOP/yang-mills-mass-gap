# HODGE-ORBIT-LAPLACIAN — the de Rham–Hodge complex on `A/G` is *kinematic*: `d` and `δ` intertwine the form-Laplacians, so the "exact-sector gap" is literally the 0-form Poincaré (= H1), the co-exact sector is its Hodge dual (same problem), and the harmonic sector is β-**rigid** cohomology; the one geometric lever (Bochner–Weitzenböck `Ric_μ⪰ρ`) is the already-banked BE-Ricci with `K_*<0`

> **2026-11-23 · NEW SPECIES — de Rham–Hodge Laplacian on the gauge-orbit space `A/G`.** 判词 (verdict). Forge the weighted Hodge Laplacian `Δ_k=dδ_k+δ_{k+1}d` on `Ω^k(A/G)` (`A=su(N)^V`, `μ=Haar·e^{−S_W}`, orbital quotient; coarse amplitude `a=⟨v,X⟩`) and seek a **Hodge gap** in the exact/co-exact decomposition that lower-bounds `G″(k→0)`. **Verdict: the complex is real but kinematic.** `Δ_0^{(μ)}=δ_μ d` **is** the YM Dirichlet-form generator (PROOF-DIRICHLET-COMPARISON), so its gap **is** the Poincaré constant `=H1` as an *operator identity* (sharper than an equivalence). The hoped separation "coarse mode = 0-form / dangerous direction = divergence" cannot exist: `dΔ_0=Δ_1d` makes the exact-1-form spectrum `=` the nonzero 0-form spectrum (offset 0), `δΔ_1=Δ_2δ` makes the co-exact sector its Hodge dual (same spectral problem), and the harmonic sector `kerΔ=H^*_{dR}(A/G)` is finite-dim, topological, **β-rigid** (RIGIDITY no-go) — it cannot jump. The only geometric lever, Bochner–Weitzenböck `Δ_1=∇*∇+Ric_μ` with Lichnerowicz `m²≥c·Ric_min`, is the BE-Ricci/O'Neill curvature already adjudicated (`RCD-SYNTHETIC-RICCI`, `B3-FUNCTIONAL-BE`, both 2/5). Free massless `[0,∞)` → `λ=0` cleanly. **#1,#4 PASS; #2,#3,#5 FAIL. NOT 5/5.**

## §1 The forge and the exact dictionary
Weighted Hodge Laplacian on `k`-forms, `δ_μ=d^*` in `L²(Ω^k,μ)`. The **0-form** operator is the smooth object: `Δ_0^{(μ)}f=δ_μ df=−Δf+∇S_W·∇f`, the generator of `E_YM(f)=μ(|∇f|²)` on `L²(A/G,μ)` (PROOF-DIRICHLET-COMPARISON §2: "the only form gapped on all of `L²` is `E_YM` itself"). Hence `gap(Δ_0^{(μ)}|_{1^⊥})` `=` `Var_ν(φ)≤(1/λ)∫|∇φ|²` `=` `G″(k→0)≥λ` `=` **H1**, as an identity of operators (`G″=1/Var(a)`, D2).

## §2 The exact/co-exact split is kinematic — duality relabels, never separates
Hodge decomposition `Ω^1=dΩ^0⊕δΩ^2⊕Harm` holds on *any* Riemannian manifold, **independent of dynamics**. Crucially `d` and `δ` *intertwine* the Hodge Laplacians: `Δ_1d=dΔ_0` and `Δ_1δ=δΔ_2`. So the **exact-1-form spectrum equals the nonzero 0-form spectrum** (the coarse mode is a 0-form — its own gradient's gap *is* the 0-form gap), the **co-exact sector equals the 2-form spectrum** (Hodge-dual, `specΔ_k=specΔ_{n−k}`), and there is **no independent "exact vs divergence" gap** to read off. The phrase "the coarse mode is a 0-form vs a divergence" is the whole problem restated: the divergence sector is not a separable safe harbour — under Hodge duality it *is* the field-strength (`2`-form) sector, fed by the same indefinite curvature.

## §3 The harmonic sector is β-rigid — it cannot carry the dynamic gap
`kerΔ_k=H^k_{dR}(A/G)=` gauge-invariant (=BRST/equivariant cohomology) Betti space, finite-dimensional. By the **RIGIDITY no-go** (INVENT-NEW-INDEX), any discrete invariant locally constant under symmetry-preserving deformation is constant along the connected `β`-family and equals its free value; the gap must be a β-**non-analytic** jump, which no cohomology dimension can supply. The harmonic forms give only the 0-eigenvalue (vacuum); the gap lives in the *non-harmonic continuous* spectrum — i.e. H1.

## §4 The geometric face is the already-adjudicated BE-Ricci
Bochner–Weitzenböck `Δ_1=∇*∇+Ric_μ`; a Lichnerowicz gap `m²≥c·Ric_min` needs `Ric_μ⪰ρ>0` on horizontal forms. This `Ric_μ` is exactly the Bakry–Émery `Γ₂`/O'Neill curvature of the submersion `A→A/G`: `K_*=infl(Hess S_W)+N/12=−g|B|+N/12<0` at the Savvidy negative mode, the O'Neill `|F|²` averaging to the gluon condensate. This is `RCD-SYNTHETIC-RICCI` and `B3-FUNCTIONAL-BE` verbatim (both 2/5); Hodge vocabulary supplies no upgrade. The singular `A/G` (Gribov, Singer) takes intersection-cohomology Hodge theory (Cheeger), but its L²-cohomology is β-rigid (§3) and the pointwise curvature obstruction survives on the smooth, μ-typical Savvidy point (C8: per-mode sign flips; `PᵀHP` min `≈−0.3`, PROOF-COARSE-POINCARE).

## §5 Acid test (clean)
Free massless: flat orbit space `⇒` Hodge spectrum `[0,∞)` continuous `⇒λ=0`, `G″=0`. Massive free: gap `m²` correctly returned — **not refuted** (unlike Jones, which returned 0 for a gapped field). Non-SUSY. PASS.

## Five-point scorecard
1. **Nonperturbative — PASS.** Exact operator identity, all β, no loops / `1/N`.
2. **Sign-controlling LOWER bound `G″≥λ` — FAIL.** Hodge gap `=G″=H1` (identity); no independent `λ`.
3. **Non-circular — FAIL (new grounds).** `Δ_0^{(μ)}=` Dirichlet form `⇒` gap `=H1`; exact/co-exact tied by `d,δ` and Hodge duality; harmonic sector β-rigid.
4. **Non-SUSY + fails free field — PASS.** Free massless `[0,∞)`, massive `m²`; exact.
5. **Uniform `≥λ` — FAIL.** The non-harmonic bottom is the unproven H1 object; cohomology is β-rigid.

## Exact failing estimate
```
target:  G″(k→0) ≥ λ  ⟺  Var(a) ≤ 1/λ                          [D2, H1]
Hodge:   Δ_k = dδ_k + δ_{k+1}d ;  0-form Δ_0^{(μ)} = δ_μ d = generator(E_YM)  [identity]
  exact:    dΔ_0 = Δ_1 d ⇒ spec(dΩ^0) = spec(Δ_0)|_{1^⊥} = Poincaré = H1
  co-exact: δΔ_1 = Δ_2 δ ⇒ spec(δΩ^2) = spec(Δ_2)            (Hodge-dual — same problem)
  harmonic: ker Δ_k = H^*_dR(A/G) = BRST cohomology, finite-dim, β-RIGID  [RIGIDITY no-go]
  ⇒ no exact/co-exact "gap separation": each bottom is a form-Laplacian gap = H1
geometric: Δ_1 = ∇*∇ + Ric_μ ;  Lichnerowicz m² ≥ c·Ric_min needs Ric_μ ⪰ ρ>0
  = BE-Γ₂ / O'Neill ;  K_* = −g|B| + N/12 < 0   [RCD-SYNTHETIC-RICCI, B3 — banked 2/5]
free massless: flat ⇒ spec Δ = [0,∞), λ=0, G″=0                [acid PASS]
missing: an independent, uniform, β-jumping lower bound in the non-harmonic spectrum — that IS H1.
```

**Sharpest line.** Hodge theory supplies the one thing no prior species owned — a *complex* (`d,δ`, exact/co-exact/harmonic) — and it is exactly what exposes the illusion: `d` and `δ` **intertwine** the form-Laplacians, so the exact-sector gap is *literally* the 0-form Poincaré constant `=H1`, the co-exact sector is its Hodge dual (`2`-form, same spectral problem), and the only finite-dimensional piece — the harmonic/cohomological sector — is β-rigid topology that cannot jump. So "separate the coarse 0-form from the divergence" is unavailable: the decomposition is orthogonal but **energy-blind**. The single geometric lever, `Bochner–Weitzenböck`/Lichnerowicz `Ric_μ⪰ρ`, is the already-refuted BE-Ricci (`K_*<0` at Savvidy) in form dress. Free massless is correctly `[0,∞)`. **#1,#4 pass; #2,#3,#5 fail — on new (duality + rigidity) and old (BE-Ricci) grounds. NOT 5/5; no `λ` fabricated.**

*Refs.* (external) Hodge, *The Theory and Applications of Harmonic Integrals* (Cambridge, 1941); de Rham, *Variétés différentiables* (1955); Weitzenböck (1923); Bochner, *Bull. AMS* **52** (1946) 776; Lichnerowicz, *Géométrie des groupes de transformations* (Dunod, 1958); O'Neill, *Michigan Math. J.* **13** (1966) 459 (submersion A-tensor); Singer, *CMP* **60** (1978) 7 (Gribov: `A/G` not a manifold); Cheeger, *PNAS* **76** (1979) 2103 (Hodge theory of cone-like singularities); Atiyah–Bott, *Phil. Trans. R. Soc. A* **308** (1982) 523 (Hodge-theoretic YM); Bakry–Émery, *LNM* **1123** (1985) 177. (internal) PROOF-DIRICHLET-COMPARISON (§2: 0-form Hodge form = `E_YM`, gapped only on all `L²`); PROOF-IR-H1 (fiber/coarse orthogonal split; H1 ⟺ `m_eff²(k→0)≥λ`); PROOF-COARSE-POINCARE (`PᵀHP≈−0.3`; `G″=S−D`); C8-SPECIES-NOGO (per-mode sign flip); RCD-SYNTHETIC-RICCI (O'Neill `|F|²`, `K_*<0`, 2/5); B3-FUNCTIONAL-BE (`Γ₂^ν=G″=H1`, 2/5); D2-COUPLING-DOEBLIN (`G″=−b′(0)=1/Var(a)`); INVENT-VARIANCE-INEQ (`G″=S−D`); INVENT-NEW-INDEX (RIGIDITY no-go); INVENT-NEW-MONOTONE (`Γ₂=Ric+Hess S_W`); W7-Gpp0 (`G″=+0.51`). No fabricated bound.
