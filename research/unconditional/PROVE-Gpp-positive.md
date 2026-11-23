# PROVE G″(0)-positive — what is provable now, and the one failing estimate

**Verdict.** No uniform-in-(V,a) proof follows from the presently available estimates. Of the two halves of H1, the **variance term** now has a rigorous upper bound (§1) and the **Jacobian curvature** a rigorous lower bound (§2); the **coarse-curvature term** admits only its *negative* Savvidy floor (§3). The bottleneck is one inequality, (★), which is numerically true (W4: blocked Q⊥ min +0.81β) and is logically equivalent to induced mass generation (W6).

## The object, and a near-collapse that fixes the strategy
On the flat link-algebra chart the exact marginal curvature is
`G″(A′)[v] = E_μ[H_vv|a=0] − Var_μ(G_v|a=0) + W(v)`,  `H=Hess S_W`, `G=∇S_W`, `a=⟨v,log U⟩`,
where `W` is the Weyl/Jacobian curvature of Haar (W7's "0.25"). Numerically **G″(0)=1/Var(a) to ≈2%** (SU2: 0.274 vs 0.268; SU3: 0.485 vs 0.45; quartic coeff ≈0), while the literal `E[H]−Var` ≈ 0 (+0.06, −0.01). So the marginal is essentially Gaussian and the surviving positive curvature is *not* `E[H]−Var` but the geometric/entropic term. This is why even a perfect bound on `E[H]−Var` would not settle the sign.

## §1 Variance term — provable upper bound
**Lemma A.** If `μ(·|A′)` satisfies a Poincaré inequality with spectral gap `λ_f` on the fibre, uniformly in `A′,L,a`, then
`Var_μ(Pᵀ∇S|A′) ≤ (1/λ_f)·E_μ[‖Pᵀ∇_⊥∇S‖²|A′] ≤ (β² N_□(v))/λ_f`,
`N_□(v)` = plaquettes in the coarse mode's support. *(Poincaré = orthogonal Brascamp–Lieb; `∇_⊥` is the fibre gradient so `Pᵀ∇_⊥∇S` is the **mixed** coarse–fibre Hessian, bounded by a per-plaquette operator norm.)* The bound uses the small off-diagonal object, not the full 2.44 force variance — that is what makes it usable. Rigorous **given** `λ_f>0`.

## §2 Jacobian curvature — provable positive lower bound
**Lemma B.** On bi-invariant SU(N) the Haar volume element in exponential coordinates is `∏_α (sin x_α/x_α)`, a product of functions log-concave for `|x_α|≤π/2`; hence `−log J` is convex and `W(v)≥0` **on the geodesic ball B_{π/2}** (SU2: ≥½g, SU3: ≥¾g — W6). Uniform in V,a. Across the cut (`|x|>π/2`) `W` is uncontrolled — W7's 10–20% SU(3) systematics.

## §3 Curvature term — only a negative floor
`H_vv = β Σ_□ (−∂²(Re tr U_□/N))[v]`. The per-plaquette operator-norm bound gives the *pointwise* floor
`E_μ[H_vv|A′] ≥ −β·C_N·(#□)`  (`C_N=O(N)`; ≈ −1.1β, W6 — the Savvidy floor).
This is negative and proves nothing about positivity. The identity-background value `g_free(v)=4Σ_μ sin²(k_μ/2) ≥ 4 sin²(π/L)` (=2.0 at L=4,k=π/2) is *not* a lower bound on `E[H]` — it is cancelled to ≈0 by the force variance.

## §4 The single failing estimate
Positivity needs, from §1–§3,
`g_free + E(H−H_free|a=0) − β²N_□/λ_f + W ≥ c₀`.
Lemmas A,B supply the last two terms; `g_free` is explicit. The one unproven inequality is the middle:
**(★)** `E_μ[PᵀHP − H_free | a=0] ≥ −g_free + β²N_□/λ_f + O(c₀)`,
i.e. the conditional mean of the Savvidy-negative coarse Hessian must stay above `−g_free`+variance-bound. This is *exactly* W4's "collective Hessian positive" — verified on 10 sampled configs, but pointwise-over-samples, so it cannot be used as a uniform black box. Its uniform version is H1.

## §5 Reformulation that pins the wall
Because G″≈1/Var(a): **H1 ⟺ ∃λ: Var_μ(⟨v,log U⟩|A′) ≤ 1/λ for every coarse transverse mode, uniformly in V,a.** Free theory gives Var ~ 1/g_free(k) → ∞ as k→0, so a uniform *upper* bound on the low-mode amplitude **is** a transverse induced mass `m_eff²≥λ`. It cannot come from Lemmas A/B (fibre and Jacobian only). (★) needs coarse log-concavity `S_eff⪰λ`, which is H1 — circular (Brascamp–Lieb on the coarse field re-assumes the conclusion).

**Ordering of conclusions.** Lemma A (variance) and Lemma B (Jacobian positivity) are rigorous now, at fixed block scale / group-ball. The coarse positivity (★) is not: its only rigorous lower bound is the negative Savvidy floor, and upgrading it is induced mass. There is no uniform-in-(V,a) theorem here; one estimate is done, the other is the whole gap.
