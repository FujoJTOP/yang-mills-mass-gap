# G1 — Theorem (fixed-volume existence)

**Scope (read first).** This is existence of the UV (`a→0`) limit of pure SU(N) Yang–Mills at **fixed spatial volume** `T⁴_L`, along the two-loop renormalised trajectory, in a fixed scheme C. It is **not** the mass gap, **not** the continuum limit `L→∞`, **not** uniqueness (the limit is subsequential), **not** clustering. The gap never enters: the single open input (`L3`) is a UV, one-derivative bound, `⊥` gap.

## Setting

Lattice `T⁴_L`, spacing `a`; Wilson measure `μ_{L,a,β(a)}=Z⁻¹e^{−S_W}∏dU`; RT `β(a)=2N/g²(a)`, `g²(a)∼1/(2b₀log(1/aΛ))`; scheme C fixed. `Obs` = gauge-invariant local `*`-algebra; `S^{(n)}_{L,a}` its Schwinger moments; `θ` = time reflection; `C_RP` = forward-half-space cone.

## Sub-lemma ledger

**PROVEN.**
- **L1 (scheme C).** Fixing C is a finite reparametrisation `g'²=g²+cg⁴`, `Λ'=Λe^{c/2b₀}`; leaves `m/Λ,√σ/Λ,a/ξ` invariant (A2 §L1). Bookkeeping, not an analytic wall.
- **L2 (perturbative all-orders lattice↔RT match).** BRST cohomology (ghost # 0, dim 4) leaves only `Tr F²` and `Tr F F̃`; Reisz lattice BPHZ + Symanzik `a`-expansion give order-by-order convergence of renormalised Green's functions with the *same* counterterm structure as the continuum, up to C's finite constants (A2 §L2). A genuine theorem.
- **Lemma B.** On bi-invariant SU(N) the Haar volume `∏_α(sin x_α/x_α)` is log-concave on `B_{π/2}`, so `W≥δ_W>0` there (SU2 ≥½g, SU3 ≥¾g), uniform in `V,a` (PROVE-Gpp §2).
- **Free-box-gap a-uniformity (P1-L0).** `λ_a(p)≥c₀|p|²`, and `c_L=λₐ(2πe₁/L)=4π²/L²(1+O(a²/L²))>0` is **a-uniform at fixed L**; zero modes are pure gauge and decouple. Sharp — and *not* `L`-uniform (that failure is G2).
- **Single-block retargeted `G″≥δ_W`.** Exact identity `G″=E[H]−Var+W`; one-loop equipartition gives `E[H]−Var=0` in the UV; the Savvidy cone is tail-suppressed (`O(e^{−c/g_j²})`); Lemma B supplies `W≥δ_W·μ(typ)`; the free-box gap supplies the a-uniform positive constant and `g_j≤g*` the small correction (A1 §4–5). This is "(c), repaired" — the δ is `W`, not `E[H]−Var`.
- **GL3.** `UVS(L)` (a-uniform renormalised moment bounds + Minlos equicontinuity) `⇒` relative compactness, diagonal subsequence, weak-`*` moments, and Minlos realisation `μ_L`. Lean: `exists_subseq_tendsto_pi`, `gl3_imp_with_minlos` (Minlos an explicit hypothesis).
- **GL4.** Reflection positivity is weakly closed (`rp_transfer`, `rp_l2_closure_trunc`); holds for *bounded* and `L²` test observables, no uniform-moment input needed. Unconditional as a lemma.

**OPEN (named).**
- **L3 (scale-uniformity across `J∼log(1/a)`).** The effective-Hessian lower bound `Λ_UV>0` must be uniform across the whole tower of `J∼log(1/a)` RG steps: `sup_a Z_a⁻¹‖H_a‖_D<∞` on a common core (A2 §L3). This is the **one-derivative face** of the marginal/Savvidy wall, strictly milder than the two-sided Hessian sign of (c), and `⊥` gap. It replaces the *a⁻⁴ resummation*: the `a⁻⁴` coefficient of `tr F²` is an asymptotic renormalon series, non-resummable; L2 supplies only the finite per-order scheme constants, and the nonperturbative a-uniform bound is exactly L3.

**NOT needed.** Uniqueness (GL3 is subsequential); `UVS-L` (stage 2 / `L→∞`); the literal `E[H]−Var≥δ` (false — equipartition, A1).

## Theorem (fixed-volume existence)

Fix `L` and compact simple simply-connected `G=SU(N)`. Along the RT in scheme C, assume only **L3** (all other sub-lemmas above are proven). Then there exist a subsequence `a_n→0` and a probability measure `μ_L` on `Obs` such that, for every `n` and every `f_i∈C_c^k(B(0,L/4))`,
```
S^{(n)}_{L,a_n}(f_1,…,f_n) → S^{(n)}_L(f_1,…,f_n)   (weak-*),
```
and `μ_L` satisfies **OS0** (regularity/moment bounds), **OS1** (Euclidean invariance on the torus — automatic from lattice symmetry), and **OS3** (reflection positivity): `∫\overline{θA}A dμ_L ≥ 0` on `C_RP`.

**Proof (assembly).** L1+L2 give the `a⁻⁴` subtraction in scheme C; the single-block `G″≥δ_W` (Lemma B + free-box gap + weak coupling) plus L3 give `UVS(L)` with a-uniform constants; GL3 turns these bounds into a subsequential weak-`*` limit and Minlos realisation `μ_L`; GL4 passes reflection positivity to the weak limit. OS0 from the a-uniform moment bounds; OS1 from lattice translation/rotation invariance. No step uses the gap.

## Honest residue

Exactly one open sub-lemma: **L3**, a UV one-derivative scale-uniformity, `⊥` gap. Everything else is proven. This theorem asserts **existence at fixed volume** (the `a→0` limit); it does **not** assert a mass gap, does **not** construct the continuum `L→∞` theory, and does **not** claim uniqueness. The route to it required *retargeting* A1's (c) from the false `E[H]−Var` to the true `G″=E[H]−Var+W`; the delta is the geometric Haar curvature `W`, not the action Hessian.

*Refs.* A1-twosided-UV, A2-scheme-match, G1-fixed-volume, P1-L0-proof, IN-1-attack, PROVE-Gpp-positive, W7-Gpp0, glue-a REPORT (GL3/GL4).
