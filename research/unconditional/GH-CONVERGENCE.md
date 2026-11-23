# GH-CONVERGENCE — measured-GH / spectral convergence is **field-blind at the base** (lattice metric-measure limit is `R⁴` for free and YM alike, gapless) and, on the **field-sensitive Wasserstein/LSI structure**, its spectral gap is *literally* `G″=1/Var(a)`; lower-semicontinuity gives the right **sign** but the value is `0` at the UV and the unknown at the IR, and the `(L,a)`-uniform floor needs the same **quantitative RG-stability modulus** as L3. **2/5 — not 5/5.**

> **2026-11-23 · NEW SPECIES — Gromov–Hausdorff / spectral (measured-GH) convergence.** 判词 (verdict). Forge RG/coarse-graining as metric-measure convergence and ask whether spectral convergence transfers a gap from a reference (or gapped limit) to the finite-`(L,a)` theory **uniformly**. The metric splits into two faces, and both fail by nameable mechanisms. **(F1) Base — field-blind.** `(T⁴_L, a·d_graph, a⁴#) →mGH,spectral (T⁴_L,Eucl,Leb) → (R⁴,Leb)`, a limit **independent of the gauge field** — identical for free Maxwell and YM at every `β`; spectrum is the flat Laplacian, `λ₁=4π²/L²→0`, `[0,∞)` on `R⁴`, **gapless**. **(F2) Field space — the Wasserstein/LSI structure.** The only field-sensitive metric-measure space `(SU(N)^{links}, W₂)` has spectral gap = LSI constant = `G″=1/Var(a)` = H1 (D2, von Renesse–Sturm). Its Mosco/spectral convergence is the limit theory of the gap itself: right-**signed** (`λ₁` lower-semicontinuous) but value `0` at the UV (Gaussian fixed point by asymptotic freedom) or the unknown at the IR (`=Δ`). **#1 PASS; #2 FAIL (sign yes, no floor); #3 FAIL (IR limit *is* the gap); #4 literal PASS but field-blind; #5 FAIL — it is the target. NOT 5/5.**

## §1 · The forge: two metric-measure spaces, only one carries the field
The base construction is *trivial and field-blind*. `(Z⁴_L, a·d_1, a⁴#) → (T⁴_L, d_Eucl, dx)` in mGH as `a→0`, then `T⁴_L→R⁴` as `L→∞`; Fukaya–Cheeger–Colding give **spectral** convergence of `Δ`. But the Wilson field does **not** backreact on the metric (no dynamical geometry), so free and interacting, weak and strong coupling, all share the *same* mGH limit `(R⁴, dx)` with spectrum `[0,∞)` and `λ₁=0`: zero gauge-field information.

## §2 · F2 — the field-sensitive structure *is* H1 (via D2)
The content must live in the **field** measure. The correct metric-measure object is the equilibrium `μ` with the `W₂` transport metric and its Dirichlet form; its Bakry–Émery spectral gap / LSI constant is, by von Renesse–Sturm (heat-flow contraction ⇔ LSI), the relaxation rate of the block-spin mode, i.e. **`G″(k→0)=1/Var(a)` exactly** (D2; Ward `G″·Var(a)=1`). Spectral convergence of this family is the Mosco/Γ-limit theory of the very constant we want: **`λ_{(L,a)}→λ_{limit}`**. Mosco convergence gives `liminf λ_{(L,a)} ≥ λ_{limit}` — the **correct sign** — but the delivered *value* is degenerate: at the UV fixed point `λ_{limit}=0` (asymptotically-free Gaussian), at the IR `λ_{limit}=Δ` (the gap). Either the bound is `0` or it is tautological.

## §3 · Why it cannot be made uniform — the missing modulus
Turning `λ_{(L,a)}→λ_{limit}` into `λ_{(L,a)} ≥ λ*` **uniformly in (L,a)** requires a **quantitative** spectral-stability modulus — the RG step's Lipschitz action on the Dirichlet form in Mosco/spectral-propinquity distance,
```
|G″_{L,a} − G″_{L',a'}| ≤ ω(|log a−log a′| + |L−L′|),   Σ_shells ω < ∞.
```
Cheeger–Colding quantitative stratification supplies such rates **only under a two-sided Ricci bound**, which the Savvidy-**indefinite** action-Hessian (`PᵀHP` min `≈−0.3`) does not meet. The modulus `ω` is *exactly* the RG-stability modulus O-MINIMAL §3 named and the tightness `(★)+ (★dec)` of L3 (G1-ASSEMBLE-100 §5). Measured-GH needs it; it does not supply it.

## §4 · Free-field acid test — literal PASS, but **field-blind**
Free massless: base limit `R⁴`, `λ₁=0`, no gap transferred — the tool correctly returns **gapless** (literal acid **PASS**). But the pass is *uninformative*: the base limit is `R⁴` for YM *and* free alike, so "failing the free field" and "failing YM at the UV" are the **same event** (asymptotic freedom). The only discriminator — the trace-anomaly non-conservation (`θ^μ_μ≠0`; HYDRODYNAMIC) or the field-space LSI (D2) — is precisely where the IR circularity bites. Geometry-detector, not gap-detector.

## Five-point scorecard
1. **Nonperturbative — PASS.** mGH/Mosco/spectral convergence uses no loop or `1/N` expansion, all `β`.
2. **Sign-controlling LOWER bound — FAIL.** Mechanism right (λ₁ lower-semicontinuous), but the emitted floor is `0` (UV) or the unknown (IR): no `G″≥λ`.
3. **Non-circular — FAIL.** IR limit gap `=Δ_phys` = the gap; the only *provable* gapped references are free-**massive** (mass inserted ⇒ circular) or `β<1/48` (P3/SZZ ⇒ redundant, non-transporting — PROOF-COARSE-POINCARE §c, Lemma SC).
4. **Non-SUSY + fails free field — PASS (literal), but field-blind.** Free→`R⁴`→gapless, no false positivity; the pass reflects *blindness*, not discrimination.
5. **Uniform `≥λ` in `(L,a)` — FAIL.** This *is* the target; requires the quantitative spectral-stability modulus = L3-tightness.

## Exact failing estimate
```
base (F1, field-blind):  (T⁴_L, a·d_graph, a⁴#) →mGH,spectral (R⁴,Leb),  SAME for free & YM;
   λ₁(flat) = 4π²/L²→0 ;  σ(−Δ_{R⁴})=[0,∞) ⇒ λ₁=0  (gapless; node-6 box gap only, ≠ Δ)
field (F2, Wasserstein/LSI):  spectral gap = LSI const = G″(k→0) = 1/Var(a) = H1   [D2; vRS]
   Mosco:  liminf λ_{(L,a)} ≥ λ_limit   (RIGHT sign)   but
   UV limit  : a→0 ⇒ Gaussian fixed pt (asymptotic freedom) ⇒ λ_limit = 0
   IR limit  : λ_limit = Δ_phys  = the gap
   uniform λ*: needs  |G″_{L,a}−G″_{L',a'}| ≤ ω ,  Σω<∞   = L3 tightness / RG-stability modulus (O-MIN §3);
   Cheeger–Colding rate needs  Ric ≥ −K two-sided ;  Savvidy H ⪰̸ 0 (PᵀHP min ≈ −0.3) ⇒ NO rate
free: R⁴ ⇒ λ₁=0 (acid PASS, but same limit as YM ⇒ no discrimination)
⇒ measured-GH/spectral = limit-theory face of the variance-upper-bound species: language, not lever.
```

**Sharpest line.** Measured-GH / spectral convergence splits cleanly and both halves fail by nameable mechanisms. The **base** limit is `R⁴` for **every** theory — field-blind, so its "free-field acid pass" reflects blindness, and it would equally "fail" interacting YM at the UV (asymptotic freedom *is* that limit). The **field**-sensitive structure is the Wasserstein/LSI object of D2/P3, whose spectral gap is *identically* `G″=1/Var(a)` = H1: Mosco gives the **right-signed** `liminf λ ≥ λ_limit`, but `λ_limit` is `0` (UV) or unknown (IR) — sign on an empty value. A **uniform `(L,a)`** floor needs a quantitative spectral-stability (propinquity) modulus of the block-spin map — exactly the L3-tightness/RG-stability modulus that O-MINIMAL §3 and G1-ASSEMBLE-100 §5 already named as the missing machine. The species *sharpens the name* ("uniform spectral gap of the RG-flowing Dirichlet form, stable in Mosco/propinquity distance") but adds no lever: it is the limit-theory face of the one variance-upper-bound species. **#1, #4(literal) pass; #2, #3, #5 fail on the identical object.** NOT 5/5.

*Refs.* (external) M. Gromov, *Metric Structures for Riemannian and Non-Riemannian Spaces* (Birkhäuser 1999); K. Fukaya, Invent. Math. **87** (1987) 517 (spectral convergence under GH); J. Cheeger, T. Colding, J. Differential Geom. **46** (1997) 406; **54** (2000) 13; **56** (2000) 189 (spectral convergence + quantitative stratification); K.-T. Sturm, Acta Math. **196** (2006) 65/133 (measured GH / RCD); J. Lott, C. Villani, Ann. of Math. **169** (2009) 903; U. Mosco, J. Funct. Anal. **123** (1994) 368 (Mosco/Γ-convergence); K. Kuwae, T. Shioya, Comm. Anal. Geom. **11** (2003) 599 (convergence of spectral structures; λ₁ lower-semicontinuity); F. Latrémolière, J. Math. Pures Appl. **103** (2015) 303 (spectral propinquity); A. Connes, *Noncommutative Geometry* (1994); M.-K. von Renesse, K.-T. Sturm, Comm. Pure Appl. Math. **58** (2005) 923 (heat-flow ↔ LSI); D. Bakry, M. Émery, *Diffusions hypercontractives*, Sém. Probab. XIX (1985) 177; Gross–Wilczek, PRL **30** (1973) 1343; Politzer, PRL **30** (1973) 1346 (asymptotic freedom ⇒ Gaussian UV fixed point). (internal) O-MINIMAL-UNIFORMITY §2–3 (same compactness-without-quantitative-modulus; `(L,a)` noncompact); G1-ASSEMBLE-100 §5–6 (L3 = `sup_a Z⁻¹‖H_a‖<∞` = tightness = H1); PROOF-COARSE-POINCARE §c (β<1/48 anchor does not transport; Lemma SC); H1-T5-NCG (spectral-triple cousin — same circular gap); P3-WASSERSTEIN-ym + D2-COUPLING-DOEBLIN (`G″=1/Var`, von Renesse–Sturm — the field-sensitive metric-measure structure *is* H1); HYDRODYNAMIC-GAP (detector-not-gap pattern); PROPAGATION-OF-CHAOS + C8-SPECIES-NOGO (same-wave format; empty-referent/false-implication). No fabricated bound.
