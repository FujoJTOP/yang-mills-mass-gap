# H2-TIEBREAKER — does the fiber factorize over blocks given U′?

**Verdict.** No — not for the P the corpus actually defines. PROVE-H2 and H2-globalization are computing in two different microscopic blockings; the disagreement is that modeling choice, not the k→0 physics. Corrected: the fiber does **not** factorize (H2-globalization's premise holds), so CONDITIONAL Step 2 is **not** discharged by block-local H2 — but the globalization is an **open UV lemma, not the IR mass gap**, so H2-globalization's attribution is wrong and P4 is not contradicted. Net: CONDITIONAL-LATTICE-GAP is **not reducible to H1 alone**, yet not entangled by extra IR either — it is **H1 (the IR wall) + one open UV lemma**.

## The decisive fact: what P is
PROVE-H3 §P fixes it: `U′_μ(B)=proj_G exp(Σ_e w_e Ad(T_e) log U_e)`, Gaussian kernel `w` at block scale; `dP` is convolution with low-pass symbol `p(k)`, `p(0)≠0`, `p(π/a)→0`. This is a **nonlocal, overlapping average**, *not* a decimation pinning the block-face links. Under it `(U′,ξ)↦U` rebuilds each fine link from a whole neighborhood, so a straddling plaquette's fiber support spans ≥2 blocks:
`H_fib(ξ;U′) = Σ_p β(1−1/N Re tr U_p)` is **not** `Σ_i h_i(ξ_i;U′)`.
Factorization ⟺ every plaquette's fiber support lies in one block ⟺ the coarse field pins *all* inter-block links. That happens only for a checkerboard decimation (odd links coarse, even fiber) — which is exactly PROVE-H2's "each boundary plaquette sees U′ as a fixed link," and is *not* the defined low-pass. (The numeric `decimate2` is a third scheme, a sharp 2-link product.) **H2-globalization is factually right: straddling plaquettes couple adjacent blocks.**

## Two obstructions, conflated under "block coupling"
- **(i) U′-dependence of the fiber conditional** (mixed Hessian; score `∂_{U'}log μ(ξ|U′)`), breaking H3(b)'s Jensen step. Its toy `V=½a²+½b²+ga²b²` has **one** block, no spatial coupling: (i) is a *single-site* effect, governed by the marginal curvature at the **soft** coarse mode — i.e. it is H1's content (`G″(0)≈+0.043−0.208=−0.165`), IR only because that mode is k→0, benign at k>0.
- **(ii) Inter-block fiber coupling at fixed U′** (straddling plaquettes), which blocks tensorization of the **fiber** term.

H2-globalization's `C(β)≫1`, HS `e^{cV}`, and "the gap between block-local and full-fiber *is* the mass-gap content" all concern (ii). But (ii) lives in the **high band**: free gap floor 2.0, against the Savvidy −2.75 curvature that can make the band Hessian indefinite — the H2-b-uniform wall. A volume-uniform Poincaré for a band-limited shell **does not imply a system gap**, so it cannot contradict P4. `C(β)≫1` defeats only Dobrushin/Zegarlinski, and `e^{cV}` is a *technique artifact* — the same status the corpus itself assigns `e^{−2βb⁴}` in H2-b-uniform §2. So globalization (ii) is **open and UV**, expected true, not IR.

## Adjudication
- **PROVE-H2** right: block-local Λ_UV proved; k→0 frozen (H3(c)); difficulty not the soft mode. **Wrong** to call globalization closed — block-local H2 does not discharge Step 2.
- **H2-globalization** right: no tensorization; Step 2 needs an unproved globalization. **Wrong** that this globalization *is* the IR/mass-gap content or violates P4.
- **CONDITIONAL-LATTICE-GAP**: not reducible to H1 alone (Step 2 needs a *global* fiber Poincaré; H3(b) needs the score repair, `λ>2 sup Var(score)`). But not entangled by IR in the fiber — that term is UV. Net obligation = **H1 (IR soft-mode wall, likely to fail) + global (V,b)-uniform high-band Poincaré (UV, open, expected true)** — asymmetric, neither installed.

*Refs.* PROVE-H3 §P,(b),(c); PROVE-H2; H2-globalization; H2-b-uniform §2; W6-BE-at-2.5; VAR-FINAL; `w4_blocked.py:decimate2`.
