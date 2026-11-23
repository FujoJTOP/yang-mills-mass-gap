# D1-BLOCK-MARTINGALE — the amplitude-variance (susceptibility) block decomposition: right object, right sign, exact identity — but the V-uniform bound on `Σ_i E[(Δ_i a)²]` *is* short-range correlation, and the certificate normalization `1/√V` is exactly the one that destroys the compactness bound

> **判词 (verdict).** D1 correctly retargets Efron–Stein from the **force** (FORGE5-CONE-VARIANCE: `Var(f|0)≤diam²E‖HessV·v‖²`, *largest on the Savvidy cone*, wrong sign) to the **amplitude** `a`: since `Var(a)=χ=G_T(0)` is the zero-momentum connected propagator, an upper bound `Var(a)≤1/λ` **is** the lower bound `m_eff²:=Γ_T(0)=1/χ≥λ` — the right object and the right direction the force route got backwards. But it cannot close: the **martingale identity** `Var(a)=Σ_i E[(Δ_i a)²]` is an equality (sign-blind, bounds nothing), and the **bound form** `≤C` is V-uniform **iff** correlations decay (`ξ<∞`) — clustering — H1. **The central no-go — a normalization conflict:** the only normalization that makes it a certificate (`|v|₂=1`, `a=(1/√V)ΣX`, `Var(a)=χ`) is the one where `a` is delocalized (range `√V`), so SU(N) compactness gives only `Var(a)≤π²dV` — V-**growing**; the normalization where compactness *does* give a V-uniform bound (`a=(1/V)ΣX`, range `O(1)`) is **vacuous** (free also `→0`). No normalization is both certificate-correct and compactness-bounded. **Acid test:** identity gives `χ_free≈L²/4π²→∞` and `χ_YM=O(1)` — a real distinction — but *proving* V-uniformity of the RHS **is** H1. 5-pt: **#1 PASS; #2 FAIL (right direction, undelivered); #3 FAIL; #4 vacuous-pass; #5 FAIL.**

## §0 The object and the correct normalization
`a=⟨v,X⟩`, `X=log U`, flat chart `μ∝e^{−V}dX`. Two normalizations.
- **Certificate (`1/√V`):** `|v|₂=1`, `v_ℓ=e^{ikx_ℓ}/√V`, so `a=(1/√V)Σ_ℓ X_ℓ`. Then `Var(a)=(1/V)Σ_{ℓℓ'}Cov(ℓ,ℓ')=Σ_r C(r)=χ`, and `χ=G_T(0)=1/Γ_T(0)` **exactly** (susceptibility = zero-momentum propagator = inverse inverse-propagator). Hence `Var(a)≤1/λ ⟺ m_eff²:=Γ_T(0)≥λ`. Free massless: `χ=1/λ₁≈L²/4π²→∞`. YM: `Var(a)≈2.2–3.7`, V-stable, `1/Var≈0.27–0.51=G″` (W7). **Right object, right direction.**
- **Magnetization (`1/V`):** `a=(1/V)ΣX=\bar X`, `Var(\bar X)=χ/V→0` for **both** free (`~1/L²`) and YM (`1/(m²V)`). No distinction ⇒ vacuous. **REJECT.**

So `1/√V` is mandatory for a certificate.

## §1 Exact identity — and why it is not a bound
Block partition `B_1…B_m`, `m=V/b`, filtration `𝔉_i`. Martingale increments `Δ_i a=E[a|𝔉_i]−E[a|𝔉_{i−1}]` obey, **exactly and for every measure**, `Var(a)=Σ_{i=1}^m E[(Δ_i a)²]`. Free: RHS `=χ_free≈L²/4π²→∞`; gapped YM: `=χ_YM=O(1)`. The V-growing-vs-V-uniform distinction is real — but this is an **equality**: it computes `χ`, it does not bound it. A certificate needs the strict Efron–Stein form with a **V-uniform upper bound** on the RHS, and that is where it dies.

## §2 The compactness conflict (the central no-go)
The per-link bound `|X_ℓ|≤π` does **not** bound `Var(a)`: `a` is delocalized, `|a|≤π√d·√V`, so range alone gives `Var(a)≤π²dV` — V-**growing**. Equivalently the martingale increment is **global**: for a long-range-correlated field, revealing block `i` shifts the conditional mean of *every* `X_ℓ` (perfectly-correlated "all-`Y`" field: one block reveals `Y`, `Δa=O(√V)`), so `(Δ_i a)²` is not `O(1/V)`. The Gibbs/conditional-resampling bound is V-uniform **iff** the block-`i` influence is confined within `ξ`: `Σ_i E[(Δ_i a)²]≲ξ⁸/b<∞` **iff `ξ<∞`**. Compactness does not give `ξ<∞`; a compact *massless* theory (4D compact U(1) Coulomb phase, `1/k²` pole) has `χ~L²→∞` with bounded links. **A V-uniform bound on `Σ_i E[(Δ_i a)²]` is precisely short-range order = clustering = the mass gap.**
```
certificate  a=(1/√V)ΣX, |v|₂=1 ⇒ range(a)=√V ⇒ compactness ⇒ Var≤π²dV  (V-GROWING, useless)
compactness  a=(1/V)ΣX,  range O(1) ⇒ Var≤π²d    (V-uniform)  ⇒ free also →0  (VACUOUS)
   ⇒ NO normalization is both certificate-correct and compactness-bounded.
```

## §3 Gauge dependence makes it worse, not better
Raw `a=⟨v,log U⟩` is gauge-variant; the longitudinal (pure-gauge) component shifts coherently under a global gauge parameter, giving `Var(a_gauge)~V` even for gapped YM (B5: longitudinal dominates `Var(a)`). So the **raw** amplitude's decomposition blows up for **both** free and YM (gauge artifact), losing the acid test's *discriminatory* content. Transverse projection recovers the physical `χ`=H1, and its gauge-fixing (Coulomb) is itself long-range: the "long-range gauge dependence" is an artifact to strip, not a tool.

## §4 Acid test (decisive)
Free massless: long-range Gaussian correlations make `Δ_i a` global ⇒ `Σ_i E[(Δ_i a)²]=χ_free≈L²/4π²→∞`; the compactness bound fails (noncompact, unbounded range). The identity correctly blows up. **But the pass is vacuous** (FORGE5-CONE-VARIANCE caveat): the finiteness invoked is range-boundedness, shared by any compact *massless* theory with the same `χ~L²→∞`; no YM positivity is extracted that a compact gapless theory lacks.

## §5 Five-point verdict
1. **Nonperturbative — PASS.** Efron–Stein/martingale exact, all β, any measure.
2. **Sign-controlling lower bound — FAIL (undelivered).** Direction right (`Var(a)≤1/λ ⟹ m²≥λ`, §0) — a genuine structural improvement over the force route — but the RHS bound is circular (#3), so no `λ`.
3. **Non-circular — FAIL.** `Σ_i E[(Δ_i a)²]≤C` uniformly ⟺ `ξ<∞` ⟺ clustering ⟺ H1; compactness only yields `π²dV`.
4. **Non-SUSY + fails free — vacuous PASS.** Free `→∞` via noncompactness, not via gap (compact massless also `~L²→∞`); no extracted positivity.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** No `λ`; the V-uniform bound is H1.

**Sharpest line.** D1 fixes the *object*: `Var(a)=χ=1/Γ_T(0)` is the susceptibility, and its upper bound is exactly a mass-gap lower bound — the correct target the force-Efron–Stein missed. But the certificate normalization `1/√V` that makes `Var(a)=χ` is the one that delocalizes `a` (`√V` range), so SU(N) compactness yields only the V-growing `π²dV`; the `1/V` normalization that compactness *can* bound is vacuous for the free field. A V-uniform bound on `Σ_i E[(Δ_i a)²]` is not a consequence of compactness — it **is** `ξ<∞`, clustering, the gap. Right object, right sign, exact identity — and the bound is H1. **No 5/5.**

*Refs.* FORGE5-CONE-VARIANCE (force-Efron–Stein, wrong sign on cone); FORGE5-SL-MARTINGALE (QV martingale, equality tight for Gaussian); SUB1-COV-BOUND (SD identity gap-blind, free `Var(a)=∞`); CONTINUUM-VARIANCE-DEFECT (`S=1/Var` is H1); PROOF-IR-H1 (coarse orthogonal term = H1); PROVE-Gpp-positive §5 (H1 ⟺ `Var(a_{k→0})≤1/λ`); INVENT-VARIANCE-INEQ (`S=1/Var(a)` Schur); B5-CORRELATION-CONTROL (longitudinal dominates `Var(a)`); OTHER-SAVVIDY-SIGN (susceptibility ceiling, `1/Var` V-stable); G2-CORRELATION-GEOMETRY (soft mode `O(1)`, `G″·Var(a)=1`); W7 (`G″≈1/Var`, L4=L6); Efron–Stein, Ann. Statist. 9 (1981) 586; Steele, Ann. Probab. 14 (1986) 1182; compact QED Coulomb-phase `1/k²` propagator.
