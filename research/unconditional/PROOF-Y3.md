# PROOF-Y3 — IR-limit uniqueness (Y3): strictly sub-gap, not closable, exact missing lemma

> **Verdict.** Y3 is **not proved and not closable from the proven pieces.** Split as in IN-6: the **VAC** half (each limit point has a unique vacuum) is free from OS4 and strictly weaker than the gap; the **limit-uniqueness** half is open. Y3 is in fact *narrower* than UNI-lim: `g_W(r)=c_r^{−1/2}` is a **neutral, local, blocked** datum, so Y3 holds iff the scalar family `{c_r}_{r>0}` is constant on the set of limit points — which follows from a **unique neutral restriction** of the infinite-volume state, i.e. from **W1₀ alone**. That is strictly weaker than UNI-lim (which also needs the string sector W1₁) and strictly below the *bc-uniform* gap; it is **incomparable** to the single-state Clay gap. No proven lemma delivers W1₀ and no 4D clustering bound exists at the relevant coupling. Exact missing lemma at the end.

## What Y3 requires
`g_W(r)=c_r(μ)`, the quadratic coefficient of the scale-`r` blocked effective action (Def. 9.17); it is built from finitely many block-local observables `tr F²`, `Re tr U_p` — all in the **N-ality-0** algebra `𝒜₀`. Hence `c_r` is a continuous functional of the *neutral* local marginals of `μ`. Consequences:
- **Y3 ⟺ `{c_r}` const on the limit set** (BREAK-IN-0 clean lemma: `c_r` is a function of `μ`).
- **Y3 ⟸ unique neutral restriction ⟸ W1₀** (L2 applied on `𝒜₀`). The string sector W1₁ is **not** needed: two phases differing only in the area law share the local neutral blocked action, so their `c_r` agree and `δ(ν)=0`. Thus **Y3 < UNI-lim strictly.**

Not automatic, because `g_W` is nonlinear: blocking is deterministic, so the mixture is the log-sum-exp `c_r(mix)=c_r(ω)+δ(ν)`, `δ` the Jensen variance correction (mixing *lowers* curvature — the mixture is wider). VAC leaves distinct extreme states, and `δ(ν)≠0` across them splits `g_W`. VAC is necessary, not sufficient.

## The three routes
**(i) Dobrushin/decoupling.** Uniqueness via the influence matrix needs `β·`(local plaquette degree)`<1` — strong coupling. At the physical `β` (SU2≈2.5, SU3≈5.7) the per-shared-plaquette influence is `2β=5≫1` (H2-globalization), so the criterion fails; multiscale/blocked Dobrushin needs `β_eff<1` under RG = IN-9 = the gap wall. **No 4D YM clustering bound exists at the relevant coupling**; only strong-coupling clustering is proved (Osterwalder–Seiler '78; SZZ arXiv:2204.12737). Route needs an input that does not exist.

**(ii) OS4 + RP + block-local H2.** OS4 ⇒ MIX ⇒ ergodic ⇒ VAC (L4): the invariant subspace of `L²(μ)` is `ℂ·1`. RP gives `T=e^{−aH}≥0`, `TΩ=Ω`, but only a *sign*; finite-`V` Perron–Frobenius isolates `λ₀` with an `e^{−cV}` constant, never volume-uniform. Block-local H2 (`Λ_UV=e^{−2βN_p(B)}/C_G`) and Lemma B (Haar log-concavity, `W≥δ_W` on `B_{π/2}`) are **intra-block/UV**; H2-globalization shows they do not globalize (inter-block `C(β)≫1`). Even globalized they give at most a neutral gap (⇒W1₀), never UNI-lim. **Route yields VAC only.**

**(iii) Transport/heat flow.** The Langevin generator's gap is fictitious-time, ≠ transfer-matrix gap (NA-1 v2: finite-`V` mixing `τ_P≈0.5` with the gap open). LSI⇒clustering needs Bakry–Émery, false at `β=2.5` (W6; `CD(m,∞)` refuted on the Gibbs measure). Its surviving slice is induced-mass = the gap. Infinite-volume uniqueness via the flow needs a well-posed 4D Langevin with a unique invariant measure = the missing marginal BPHZ theorem (W2). **No independent closure.**

## Exact missing lemma
**Lemma W1₀ (neutral boundary screening).** There is `ε₀(d)→0` (any rate) such that for all `Λ₀⋐Δ`, all `A∈𝒜₀` supported in `Λ₀`, `‖A‖≤1`, and all boundary data `ξ,ξ′` agreeing off distance `d`, `|γ_Δ(A|ξ)−γ_Δ(A|ξ′)|≤ε₀(d)`. By L2 this gives a unique neutral Gibbs state, hence bc-independence of every `c_r`, hence Y3. It is implied by (neutral) bc-uniform MG-spec, i.e. by the gap; conversely no rate-free proof of W1₀ is known (no FKG/positive association in non-abelian YM). **This is the exact missing input** — strictly weaker than the gap and than UNI-lim, but unreachable from the proven pieces.

**Sharpest line.** Y3 needs only the *neutral* limit to be unique — strictly less than IN-6's UNI-lim — but even that is W1₀, which the proven UV/block-local structure does not supply and which no 4D tool produces below strong coupling.

*Refs.* IN-0-attack; IN-6-attack; BREAK-IN-0; uniqueness/REPORT.md L2,L4,§5,§5.5; PROVE-H2-uvgap; H2-globalization; P4-RP-uniform-gap; H1-T1-stochastic-quant; IDEA-wasserstein; W6-BE; SZZ arXiv:2204.12737; Osterwalder–Seiler, Ann. Phys. 110 (1978) 440; Georgii, *Gibbs Measures and Phase Transitions* (2011) Ch. 7–8.
