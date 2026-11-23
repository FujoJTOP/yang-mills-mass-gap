# BREAK-IN-7 — attempting `Λ>0` for the continuum limit

**Verdict: does not close.** The chain `b₀>0` → two-loop flow → finite `λ_w` (IN-3) → `Λ>0` fails structurally: **`b₀` is the linearization of the RG flow at the Gaussian fixed point, while `Λ>0` says the limit is *not* the fixed-point trajectory. A coefficient of a series around a point cannot select which trajectory passes through it.**

**Steelman.** (1) `b₀=11N/(48π²)>0` (IN-0). (2) Lemma DT: along the RT `Λ:=μ₀e^{−1/(2b₀g²(μ₀))}`, with `Λ>0 ⟺ g²(μ₀)>0` at some finite `μ₀`; `Λ=0 ⟺ g₀≡0`. (3) IN-3 two-loop: `NPB: |β−β_{≤2}|≤C₃λ⁴` on `0<λ≤7.3`, closed form `t_up^{2L}=C_RG/λ−(51/121)ln(1+121C_RG/(51λ))`; reaching `λ_w≈7.3` at `μ_phys` gives `g²(μ_phys)>0`, hence `Λ>0`.

**Why (3) fails.** It is circular: finite `λ_w` at the physical scale *presupposes* nonzero coupling on the trajectory — `N2` restated. The two-loop flow says how a `g` runs *given* a `g`; it does not certify the limit *measure* carries one.

**The a⁻⁴ / renormalon / triviality triad (what blocks `b₀→Λ`):**
- **(a) a⁻⁴ subtraction.** The only gauge-invariant scalar is `tr F²`, with `⟨tr F²⟩=c/a⁴+[finite]`. The finite part carrying `g_ren`/`Λ` exists only after an RT-consistent (`C`) subtraction whose `a`-uniformity and **limit-independence** are unproven (IN-0 Y3; IN-6). `b₀` is a UV log-coefficient, constraining the divergent piece, not the finite sign.
- **(b) Renormalon.** `b₀,b₁` lead a divergent β-series with IR renormalons on the positive Borel axis; the 4D transseries parameter is **undetermined** (W2(ii); 07-resurgence). So `Λ` is ambiguous at `O((Λ/μ)²)`, which `b₀` cannot resolve.
- **(c) Triviality / limit-selection.** Continuum limit is `β(a)→∞` (`g₀(a)→0`); free (`Λ=0`) and nontrivial AF trajectories share `g₀(a)→0` and the *same* leading log. Distinguishing them is proving `ξ/a→∞` — the marginal-AF counterpart of φ⁴₄ triviality, with no rigorous-RG/SPDE theorem (05-/12-spde).

**The one escape is the same wall:** fixing a physical scale by construction, e.g. `a²σ_lat(β)→c_σΛ²>0` along the RT — but that is confinement L7 (`σ`-scaling to `β→∞`), i.e. IN-3 again.

**Conditional lemma (its hypothesis *is* N2).** If (i) `λ(μ_phys)=λ_w∈(0,7.3]` with NPB, and (ii) the renormalized data are `a`-uniform and limit-independent (Y3), then `Λ>0` and Theorem N gives `¬IsGaussian`. But (i)∧(ii)∧`g²>0` = `N2` verbatim.

**EXACT obstruction.** `b₀>0` fixes the tangent to the flow at the Gaussian FP; `Λ>0`-for-the-limit is the selection of a non-FP trajectory surviving `a→0`. The bridge is an IR/nonperturbative input — the a⁻⁴/renormalon/triviality triad — none available. IN-7 reduces not to `b₀` but to `N2` = "the limit is not free," i.e. IN-3 (intermediate coupling) ∧ IN-6/Y3 (limit independence), both open.
