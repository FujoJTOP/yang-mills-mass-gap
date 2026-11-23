# PADE-ROUCHE-DEEP — the `[M−1/M]` Padé of the KL Stieltjes function is a correct, acid-clean band detector, but it needs `M≳2/(am₁)` exact moments — bounded (`~O(10)`) at one β=2.5 lattice, divergent as `a→0` — and finitely many moments + positivity cannot *certify* a gap; **RELOCATES** to marginal BPHZ/transport. **2/5.**

> 2026-11-23 · Target: derive `I2` (no `ν`-weight in `(m₁,Λ)`, `Λ=(1.5m₁)²`) at β=2.5 from the `[M−1/M]` Padé of `G(z)=∫dν(λ)/(1−zλ)`, using only finitely many exact moments `{C(k)}_{k<2M}` and `ν≥0`. **判词 (verdict, up front): the lever is real and free-field-failing — it rejects the two-mass GFF exactly — but it fails the two questions that decide Clay. (i) The isolation gap *collapses in the only variable where `C(k)` are power moments*: `λ=e^{−aE}`, relative `z`-gap `z_Λ/z₁−1=e^{a(√Λ−m₁)}−1≈½am₁→0`, so the signal-to-noise count is `M≳1/(½am₁)=2/(am₁)`, finite (`≈10`) at fixed β=2.5 but `→∞` as `a→0`. (ii) Finitely many moments + positivity cannot *certify* absence of an atom (truncated Stieltjes indeterminacy), so a finite Padé certifies the approximant, not `G`. Net: it consumes `O(1/(am₁))` uniform exact moments = W3's non-summable `Σ_jg_j²`. Scorecard 2/5 — NOT 5/5.**

## 1 · The clean object
`G(z)=Σ_k C(k)z^k=∫_{(0,1]}dν(λ)/(1−zλ)`, `λ=e^{−aE}`; `ν≥0` makes `G` a Stieltjes/Herglotz function, poles at `z=1/λ=e^{aE}`. The `[M−1/M]` Padé `P/Q` is fixed by `m_0…m_{2M−1}=C(0)…C(2M−1)`; its denominator roots are the **`M`-point Gauss nodes** of `ν` (Stieltjes' S-fraction), and `P/Q→G` on compacta off `supp`. Atom at `m₁` ⇒ pole `z₁=e^{am₁}`, residue `f₁²`; isolation `supp ν_rest⊂[Λ,∞)` ⇒ *one* pole below `z_Λ=e^{a√Λ}`, `√Λ=1.5m₁`, none in `(z₁,z_Λ)`.

## 2 · The gap collapses — required `M`
The correlators are power moments of `ν` **only in the transfer variable** `λ`. There the annulus `(z₁,z_Λ)=(e^{am₁},e^{1.5am₁})` has relative width `δ_λ:=z_Λ/z₁−1≈½am₁`. The atom's moment signature vs the continuum is the ratio `(λ_Λ/λ₁)^k=(1−δ_λ)^k≈e^{−kδ_λ}`; the atom resolves once `kδ_λ≳1` — the corpus template `k\ln(1/r)≳1` with `r=λ_Λ/λ₁≈1−δ_λ`. Hence
**`M ≳ 1/δ_λ = 2/(am₁)`.** (Estimate: am₁≈0.2 on a β=2.5 SU(3) lattice ⇒ `M≈10`, `2M≈20` moments `C(0..19)`.) **No reparametrisation opens the gap:** in `s=aE` the `C(k)` are *exponential*, not power, moments, so the Stieltjes–Padé machinery does not apply. The gap is fixed in *mass*, zero in the *eigenvalue*; `am₁→0` forces `M→∞` — the exact unboundedness Clay forbids.

## 3 · Finitely many moments + positivity cannot certify absence
Even granting exact `m_0…m_{2M−1}`: the truncated Stieltjes moment problem has a Nevanlinna family of representing measures; for *any* open gap interval some member carries an atom there while matching all `2M` moments. So "no pole in `(z₁,z_Λ)`" for the *Padé* is not a certificate for `G`. Certification needs either `M→∞` or an a-priori state-density bound — verbatim the convergence `I2`. The pure two-atom GFF is the degenerate exception (`[1/2]` exact, hence non-vacuous), but real YM carries a continuum, not two atoms.

## 4 · Acid test — PASS (the one clean win)
*Two-mass GFF* `wδ_ε+(1−w)δ_S`: exact moments recover the rational `G`, poles at `1/ε,1/S`; Rouché rejects the fake gap once a node resolves `ε`. *Free massless `trF²`*: pure cut, no pole ⇒ "no atom" ✓. *Free massive `:φ²:`*: cut on `[4m²,∞)`, no pole ⇒ `f₁²=0` ✓. Correct *and* free-field-failing — a true new detector, not a repackaged inequality.

## 5 · Cross or relocate? **RELOCATE.**
The `2M=O(1/(am₁))` *uniform exact moments* are precisely what W3 cannot supply: KP/P3 moments are rigorous only at β<1/48; at β=2.5 the marginal scale-sum `Σ_jg_j²∼log log(1/a)` is non-summable (no cutoff-uniform BPHZ bound). Renormalising the high moments *is* the marginal 4-D gauge BPHZ problem (CONTINUUM-FREDHOLM, W3, G1). Padé–Rouché is a downstream consumer of the transport wall, not a crossing of it.

## 6 · Five-point scorecard (target = uniform `I2`; NOT 5/5)
1. **Nonperturbative — PASS in form** (moments from KP/H; no loop expansion).
2. **Sign / lower-bounds edge — FAIL.** Rouché gives an *absence* (upper) certificate; residue only `f₁²≥0` (positivity is Schur-concave).
3. **Non-circular — FAIL.** Uniform certificate ⇔ uniform `M` ⇔ the state density = `I2` itself.
4. **Fails free field — PASS.** §4.
5. **Uniform `≥λ` to `2.5`/`a→0` — FAIL.** `M∼2/(am₁)→∞`; exact moments absent at β=2.5.
**2/5.**

**Sharpest line.** Padé–Rouché is the sharpest *band* probe writable from moments — Stieltjes-Padé poles are Gauss nodes and Rouché counts them, acid-clean. But the isolation gap is fixed in *mass* and shrinking in the *transfer eigenvalue*, where alone the correlators are power moments: resolving it demands `M∼1/(am₁)→∞`, and finitely many moments never certify a gap they cannot uniquely represent. It localises `I2` to "transport `O(1/(am₁))` exact moments uniformly" — the marginal-BPHZ log-sum wall. Clay stays 0%.

*Refs.* LOC-DETERMINANT (§4 lever); INVENT-NEW-IDENTITY (determinant/Stieltjes identity); G2F-ATOM-PLUS-ISOLATION; G2D-ISOLATION-I2 (`λ₂/λ₁` separation; 5-pt format); G2B (§1 positivity upper-only); CONTINUUM-FREDHOLM (`z₂−z₁≈a(E₂−m₁)`; renorm `D_O` = marginal BPHZ); W3-marginal-BPHZ (log-sum, no cutoff-uniform bound); PROOF-ISOLATION §4 (`Λ=1.5m₁`); SUB2-I2-STRONG. External: Stieltjes (1894); Baker–Graves-Morris, *Padé Approximants*; de Montessus de Ballore; Nevanlinna (truncated moment problem, indeterminacy); Gautschi (Gauss quadrature / orthogonal polynomials, nodes = Padé poles); Kotecký–Preiss CMP **103** (1986) 19; Shen–Zhu–Zhu CMP **400** (2023) 805.
