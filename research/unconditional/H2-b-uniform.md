# H2-b-uniform — is the UV fiber gap uniform in block size b?

**Verdict.** H2-b-uniform (`Λ_UV ≥ λ₀>0` independent of the block side b, at β≈2.5) is **not proved and not disproved**. The Holley–Stroock loss `e^{−2β·O(b⁴)}` is a **method artifact**, not a physical decay: the fiber is band-limited, so its *free* gap and correlation length are b-uniform and the fiber is physically a massive UV shell. But no available technique upgrades HS to a b-uniform constant, and **all four prescribed routes die for one shared reason**: at β=2.5 the UV-band Hessian is indefinite, its negative eigenvalues matching the free gap, so there is no "free gap dominates interaction" margin to localize on. H2-b is *not* the IR gap (the k→0 mode lives in U′, excluded from the fiber), yet proving it is blocked by the same β=2.5 strong-coupling wall. It hands no leverage on the gap.

## 1. The fiber is band-limited — free gap and ξ are b-uniform
`P` is a fixed factor-2 low-pass; the fiber `ξ = ker P` is the high-frequency half, wavevectors `π/2a ≤ |k| ≤ π/a`. Interval-arithmetic check: the band-bottom free eigenvalue `min 4Σ_μ sin²(k_μ/2) = 2.0` (a=1), **independent of b**. Growing the block to side b adds modes but keeps the floor at `2/a²`; the fiber's correlation length is `ξ_fib ≤ O(a)`, fixed. A field with fixed frequency floor has no soft mode at any b. So the physical fiber gap is expected **b-uniform**, equal to that of an O(a)-patch.

## 2. Why HS loses e^{−2βb⁴}
HS charges `osc(βH_fib) ≤ 2β N_p(B)`, and `N_p(B) ~ 3(b/a)⁴ → ∞`. That is the *total* oscillation sum over an extensive block. But the gap of a sum of gapped, short-range-correlated patches of size `ξ_fib ≤ O(a)` is the **min** over patches, not degraded by their number. The volume factor is therefore bookkeeping, not a physical `b⁴` decay.

## 3. The four routes, one death
The free fiber gap is `2/a²`; the Savvidy/interaction scale is `≈ −1.1β = −2.75` at β=2.5 (W6). Negative exceeds positive: **the band Hessian is indefinite** (MARGINAL: coarse restriction min ≈ −0.3; 90–97% of worst negative BE directions sit *in the fiber*, W6-BE). Every route needs a *local* positivity/contraction that this kills:
- **(a) Renormalized Dobrushin** — the fiber is already the smallest scale (a single fixed band, no sub-scales to renormalize); a single-scale local Dobrushin needs `β·deg<1`, false at β=2.5.
- **(b) Block-level BE/LSI after coarse-graining** — contradicted by data: the fiber is the *least* convex part, so "genuinely convex at the block level" is empirically false.
- **(c) Reflection positivity** — a sign inequality only; it gives no size of `1−λ₁` (P4 §1).
- **(d) Dimension count** — `N_fib ~ (b/a)⁴` modes each with free gap 2/a² do argue for b-uniformity (sum of gapped patches), but the within-patch local bound again needs free-gap > interaction, absent here.

## 4. Genuinely gapless as b→∞? No.
Fixed frequency floor ⇒ no expected gaplessness; β=2.5 is confined, UV modes massive in lattice units. **But a proof needs gauge-invariant clustering at rate `1/ξ_fib`; plaquette–plaquette connected functions are glueball correlators decaying at the *physical* mass m.** On measure-typical U′ that demands `m>0` — i.e., the gap — while the only fiber-localized soft content is the Savvidy mode on atypical strong-field backgrounds. So H2-b is a genuinely *different, UV* statement that nonetheless requires the same β=2.5 analyticity/contraction that fails, **not** the k→0 marginal mode (CONVERGENCE §c).

## 5. The honest decider (finite, not the wall)
The fiber is band-limited and β-bounded, so `Λ_UV(b)` is a **certifiable finite object** for `b=2a,4a,8a` (block-measure spectrum / interval arithmetic on the band) — unlike the IR gap. Plot `Λ_UV(b)` vs b: a **flat nonzero plateau** ⟹ H2-b TRUE and the obstruction is only HS; a **decay** ⟹ FALSE, rate to be compared against `e^{−2βb⁴}`. No fabricated proof; this is the test that settles it.

*Refs.* PROVE-H2-uvgap; W6-BE-at-2.5; PROVE-Gpp-positive (Lemma A); MARGINAL; GAP-G2-mass(c); P4 §1; CONVERGENCE-redteam(c).
