# B3-ISOPERIMETRIC — SHAPE B round 3: building the coarse Poincaré constant ρ **directly** via Cheeger/isoperimetry on ν. The fiber's own LSI is provably *silent* on exactly the sets that define h(ν) (fiber-saturated ones carry zero fiber boundary); diameter/spectral-volume routes are the variance itself (circular) or need a Ricci = Hessian lower bound (=H1); two-scale collapses to the soft shell (=H1). ρ>0 uniform in L ⟺ H1.

> **2026-11-23 · 判词 (verdict). No constructive ρ exists.** `dν(a)=∫δ(a−P̄A)e^{−S(A)}dA_fiber` is a **marginal**; its Poincaré constant is a functional of `ν` alone. I prove the obstruction: the fiber conditionals `{μ(·|a)}` and the marginal `ν` are two **independent coordinates** of the joint, so *no* fiber property — LSI, isoperimetry, Prékopa convexity — can constrain `ρ(ν)`. Concretely, `ν`-isoperimetric competitors are **fiber-saturated** sets `A=P̄^{−1}(B)`, on which the fiber-tangent boundary measure is *identically zero*: the fiber LSI lower-bounds a boundary that vanishes there. Cheeger gives `ρ≥h(ν)²/4`; bounding `h(ν)>0` is bounding the transverse bottleneck of `ν` = the conditioning defect = **H1**. **Passes #1, #4(vacuous); fails #2, #3, #5.**

## The two-coordinate obstruction (new no-go)
The joint splits exactly as `μ = ν(da)·μ(·|a)`; given `{μ(·|a)}` the marginal `ν` is a **free** choice (and conversely). Every fiber property is a functional of the conditionals at fixed `a`, hence *independent* of `ρ(ν)`: no fiber bound upper-bounds `1/ρ=Var_ν(a_{k→0})`. This is the isoperimetric mirror of PROOF-IR-H1 — "no `∇_⊥` bound gives a `∇_∥` bound," restated as "conditionals do not constrain the marginal." The only tie of `ν` to the fiber is the *specific* coupling `e^{−S}`, and using it to bound `ν` **is** H1.

## Route 1 — Cheeger with the fiber's own LSI: silent on fiber-saturated sets
`ρ(ν) ≥ h(ν)²/4`, `h(ν)=inf_B ∂_ν B / min(ν(B),ν(B^c))`. Coarea: for a fiber-saturated `A=P̄^{−1}(B)` every fiber slice `A_a` is the whole fiber or empty, so its fiber-perimeter `P_{μ(·|a)}(A_a)=0`. The fiber LSI (an isoperimetric bound *within* each fiber) lower-bounds `P_{μ(·|a)}` only on proper subsets — never on `P̄^{−1}(B)`. So `h(ν)` sees **only** the transverse boundary of saturated sets, governed by `dν`, nothing in the fiber; slicing lower-bounds `μ`'s isoperimetry (a `∇_⊥` statement), never `ν`'s. **The fiber LSI is provably silent on the Cheeger competitors for `ν`.** The residue — `h(ν)>0` uniformly — is the transverse bottleneck of the marginal, `Var_ν(a_{k→0})≤1/λ` = H1 (Cheeger/Buser: `h(ν)≍√ρ`).

## Route 2 — diameter / spectral-volume: circular or Ricci=H1
Lévy–Gromov / spectral-volume isoperimetry need `Ric≥(n−1)K` and dimension — Hessian data; on `ν` the Hessian is `Hess V_eff=G″` (min `−0.3` on the k→0 slice), false at 2.5. A bare "diameter" bound is `diam²≈Var(a)=1/ρ` — the object itself, circular by one symbol. Prékopa/Borell give `ν` log-concavity ⇒ KLS Poincaré only if `S_W` were convex; it is Savvidy-indefinite, and only the Haar sub-factor `Φ` is log-concave, its isoperimetry at the log-cut `F≈0` anti-correlated with the `∝|F|` driver (OTHER-H1-POINCARE).

## Route 3 — two-scale: `ρ_stiff` gapped, `ρ_soft`=H1
Split `ℰ_ν=ℰ_soft⊕ℰ_stiff` by momentum shell (orthogonal under the *free* Gaussian only; interactions couple). `ρ=min(ρ_stiff,ρ_soft)`. Stiff (high-k, fiber-shadowed) shells are gapped (block positivity/SZZ); the soft `k→0` shell is where Savvidy dominates (`k²−2g|F|` indefinite), `ρ_soft` *is* H1, and it is the smaller constant (hence the target). Cross-shell coupling is the running interaction, unsigned. A reorganization, not a positivity.

## Exact failure
```
ρ ≥ h(ν)²/4,   h(ν) = inf_B ∂_ν B / min(ν(B),ν(B^c)),  B fiber-saturated;
∃ρ>0 uniform in L  ⟺  h(ν)≥2√ρ  ⟺  transverse bottleneck of ν bounded below
                   ⟺  Var_ν(a_{k→0}) ≤ 1/λ  ⟺  m_eff²(k→0) ≥ λ  (H1).
Fiber LSI active on proper subsets of fibers (P_{μ(·|a)}(E)>0) — NOT on P̄^{−1}(B), where it is 0.
```

## Five-point verdict (honest)
1. **Nonperturbative — PASS.** Cheeger + coarea + Prékopa; no loop expansion.
2. **Sign-controlling LOWER bound — FAIL at the bottleneck.** Fiber LSI cannot sign `h(ν)`; the transverse bottleneck is indefinite (`PᵀHP` min −0.3), its lower bound is H1.
3. **Non-circular — FAIL at closing.** R1 needs `h(ν)>0` (H1); R2 is `diam²=1/ρ` or `Ric=G″⪰0` (false); R3 is `ρ_soft`=H1.
4. **Non-SUSY + fails free field — PASS (vacuous).** Free massless: `ν` Gaussian, `h(ν)=0`, `ρ=0`, `f₁²=0` — correctly no positivity; but the object is *acid-test-blind* (free transverse modes also carry a fiber LSI), so it never discriminates YM from free.
5. **Uniform LOWER `m_eff²(k→0)≥λ` — FAIL.** Only genuine constructive `ρ` is Lemma SC at β<1/48 (`ρ=C_blk^{−2}λ_μ`) — exactly where the full LSI already holds; nothing transports to 2.5.

**Sharpest line.** Every fiber tool is a functional of the *conditionals*; `ρ` is a functional of the *marginal*; these are independent coordinates of `μ=ν·μ(·|a)`, so the fiber cannot bound `ρ`. Literally: the `ν`-Cheeger competitors `P̄^{−1}(B)` carry **zero** fiber boundary, so "Cheeger with the fiber's own LSI" is empty by construction. Diameter/spectral-volume return `1/ρ` (circular) or demand `Ric=G″⪰0` (H1, false); two-scale names `ρ_soft`=H1. The isoperimetric route does not cross the wall — it terminates at the bottleneck, which is the conditioning defect = H1.

*Refs.* PROOF-IR-H1 (orthogonal split; `∇_⊥ ≠ ∇_∥`); PROOF-DIRICHLET-COMPARISON (fiber form kernel = coarse functions; only gapped form is H1); PROOF-COARSE-POINCARE (Lemma SC β<1/48, `H_ff` indefinite, `G″` exact); OTHER-H1-POINCARE (Haar half, log-cut anti-correlation); FORGE5-CONE-VARIANCE (#4-vacuous precedent); k0-CONSTRAINT-EFFPOT (`Hess V_eff=G″`); MARGINAL (`PᵀHP` min −0.3); Shen–Zhu–Zhu CMP 400 (2023) 805, arXiv:2204.12737; Cheeger 1970; Buser 1982; Gromov–Lévy; Prékopa 1971; Borell 1975; Bobkov–Ledoux / KLS; Ledoux (LSI pushforward); Bakry–Émery.
