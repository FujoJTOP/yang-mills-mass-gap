# F2-KREIN-RUTMAN-QUANT — quantitative Birkhoff–Bushell cone-contraction on `T_{0⁺⁺}`: a genuinely sign-correct gap bound that pins the *bottom* `m₁`, not the atom `f₁²`, with uniform rate `= G″ = H1`

> **2026-11-23 · NEW SPECIES F2 (quantitative Krein–Rutman).** 判词 (verdict). The qualitative Perron–Frobenius locks the vacuum; the glueball is the *second* `0⁺⁺` eigenstate in the same sector. Forge a **quantitative** Hilbert-metric (Birkhoff) / Bushell cone-contraction gap on an RP-adapted cone that separates the node-0 vacuum from the node-1 glueball. **Verdict: F2 is real and is the one species with the right *sign* — Bushell genuinely lower-bounds a spectral gap (no variance/moment blindness).** But the eigenvector it pins is the *bottom* of `Ω^⊥∩(0⁺⁺)` `= e^{−a m₁}`, i.e. **`m₁`, the gap** — a `bottom` functional, exactly the class PROOF-EDGE-E2 proves blind to `λ₁−λ₂` (I2). Its uniform contraction rate is the marginal curvature `G″(0)=S−D`, **verbatim H1**; the only cone that would separate the glueball-*atom* from the `0⁺⁺` continuum is sign-indefinite and **presupposes `f₁²>0`** (INVENT-SECOND-CONE). The free-field acid test is *passed* (`k→0` flat ⇒ `G″=0` ⇒ bound `0`), but only *vacuously* — without H1 the bound is `0` for YM too. **#1,2,4 PASS; #3,5 FAIL = the entire Clay content. No closure.**

## §1 The cone and the forged bound
On `H_phys=L²(𝒜/G,μ̄)` restrict to `0⁺⁺`; `T=e^{−aH}`, `1=λ₀>λ₁≥…`, `λ_n=e^{−aE_n}`. Pick the `k→0` block-uniform mode `a=⟨v,X⟩` as the *coarse/radial* coordinate and the RP-adapted cone
```
CONE C_rad := { f : f = φ(a)⊗1_fast + … , φ ≥ 0, φ single-signed (≤0 sign changes) in a }.
```
`Ω≡1` (`φ` constant) is **interior**; the radial excitation `ψ₁^{0⁺⁺}` (node-1 in `a`) lies on `∂C_rad` — the requested vacuum/glueball separation. **Birkhoff:** a positive `T` with finite projective diameter `Δ(T)=sup_{u,v∈C}d_H(Tu,Tv)` contracts the Hilbert metric, `c(T)=tanh(Δ/4)`. **Bushell:** if `c<1` then `λ₀` is simple and **`|λ|≤c·λ₀` for every `λ≠λ₀`** — a *concrete* spectral-gap floor `a·m_eff ≥ −ln c(Δ)>0` read straight off the cone-contraction ratio. This is the forged object: a lower bound on the second eigenvalue from projective contraction, no loop expansion.

## §2 What the cone's "second eigenvalue" *is* — bottom, not atom
`T` restricted to `0⁺⁺` has spectrum `{e^{−a s}: s∈spec(H|_{0⁺⁺})}`. Bushell bounds `λ₁`, the **bottom of `Ω^⊥∩(0⁺⁺)`**, `=e^{−a·inf spec(H|_{0⁺⁺})}`. That bottom is a *continuum edge whenever the lightest `0⁺⁺` is two-particle* — in the free field it is the two-photon cut at `s=0`. Bounding `λ₁<1` certifies a **gap** (`m₁>0`); it does **not** certify the bottom is an **isolated atom** (`f₁²=ν({m₁²})>0`) separated from the `0⁺⁺` continuum above. This is EDGE-E2 exactly: Poincaré/SZZ/GEVP/Bushell all encode only `μ₁=−ln λ₁=a m₁` and are **blind to `λ₁−λ₂`**. F2 changes the *constant*, not the *functional*: it re-quantifies the gap the qualitative PF already gives, and never touches `f₁²` or I2.

## §3 The rate is `G″ = H1`; the sharp cone presupposes the atom
Two exact obstructions, matching INVENT-SECOND-CONE. **(i) 4D operator, orderless.** On the true `G^E/G` action there is **no order** — the plaquette product `U_μU_νU_μ^{-1}U_ν^{-1}` is non-monotone and no section descends (PROOF-F1SQ iii) — so no Sturm/Gantmacher–Krein sign cone exists; the only invariant cone is the **full positive one** (positivity-improving, PROOF-GAUSS-INVARIANT-PF), which contains `ψ₁` on its boundary and separates nothing. **(ii) 1-D block diffusion, nested.** The one coarse generator that *does* act on the radial variable is the block-spin diffusion with drift `b(a)=G″\,a` (von Renesse–Sturm/Bakry–Émery; D2-COUPLING-DOEBLIN); its Gaussian kernel *is* TP₂, so the Sturm sign-change cone **is** invariant — but it is **nested above positivity**, so the strictly-positive free field lies in the interior of *every* level: the cone is **vacuous exactly where the acid test bites** (it certifies a finite-`V` gap for the free field too). On that generator the contraction rate is `−b′(0)=G″=S−D`, so uniform `c<1` `⟺ G″(0)≥λ>0` `⟺` **H1**. The *sharp* cone that would split the atom from the continuum is sign-indefinite (below positivity, INEQ-T3), where Schur-convexity/isolation is already assumed. **Vacuous (i/ii) or circular (the rate); no non-vacuous sharp cone.**

## §4 Free-field acid test
`g=0`: `k→0` is a flat direction, `b≡0`, `G″=0` ⇒ the radial generator is pure diffusion ⇒ `c=1`, bound `0`; equivalently `Δ_V→∞` at `∞`-volume so `c_V→1`. **Acid test PASSED** — the bound correctly gives no gap for the gapless free field. *But the degeneration is universal*: the YM positive-Haar kernel is equally spread at `∞ V`, so `c→1` there too **unless H1** restores a uniform drift. The pass is therefore *vacuous* (0 for free **and**, absent H1, 0 for YM); it discriminates only via `G″>0` — the unproven wall (the FORGE5-CONE-VARIANCE "nominal #4").

## Five-point scorecard
1. **Nonperturbative — PASS.** Exact transfer kernel, all β; Birkhoff/Bushell on the cone, no `1/N`/`g`.
2. **Sign-controlling `m_eff²(k→0)≥λ` — PASS (shape/sign).** Genuine lower bound, the one species not sign-blind — **but it bounds the bottom gap `m₁`, not the atom `f₁²`** (the pass is the sought shape on the wrong functional, §2).
3. **Non-circular — FAIL.** Uniform rate `=G″=H1`; sharp sign-separating cone presupposes `f₁²>0`.
4. **Non-SUSY + fails free field — PASS (clean).** `G″=0` ⇒ `c=1` ⇒ `0`; non-SUSY. Vacuous-uniform caveat (§4).
5. **Uniform `≥λ` in `(V,a)` — FAIL.** `c_V→1`; `c₀<1 ⟺ G″(0)≥λ`.

## Exact failing estimate
```
Bushell: c(T)=tanh(Δ/4)<1 ⇒ |λ|≤cλ₀ ⇒ a m_eff ≥ −ln c(Δ)      [a lower bound — right sign]
but  λ₁ e^{-a m₁} = bottom of spec(H|0⁺⁺) = GAP, blind to λ₁-λ₂   [EDGE-E2 bottom functional]
4D G^E/G orderless ⇒ no Sturm cone; only full positive cone (holds ψ₁) [INVENT-SC i; F1SQ iii]
1-D block diffusion: Gaussian kernel TP₂ ⇒ Sturm cone invariant BUT nested ⊃ positivity
  ⇒ free field interior to every level ⇒ vacuous at the acid test      [INVENT-SC ii]
that generator: drift b(a)=G″a ⇒ rate = G″ = S−D  [von Renesse–Sturm; D2]
  ⇒ uniform c<1 ⟺ G″(0)=S−D ≥ λ > 0 ⟺ H1 .
sharp cone (atom vs continuum) = sign-indefinite ⇒ presupposes f₁²>0  [INEQ-T3]
free field: k→0 flat ⇒ G″=0 ⇒ c=1, bound 0  ⇒ acid test PASSED (but vacuously: 0 for YM too).
```

**Sharpest line.** F2 is the only positive species whose native inequality points the *right way*: Bushell gives a genuine *lower* bound on a spectral gap from a projective contraction ratio — exactly what the variance/moment routes (sign-blind, upper-only) cannot do. The catch is twofold and is the whole problem: (i) the "second eigenvalue" the cone pins is the **bottom** of `Ω^⊥∩(0⁺⁺)`, i.e. `m₁` — the gap, a `bottom` functional that PROOF-EDGE-E2 proves encodes only `λ₁` and is blind to `λ₁−λ₂`, so F2 re-quantifies SZZ's gap and never reaches the atom `f₁²` or the isolation I2 (`SUB2-I2-STRONG`'s `Be^{−δt}`, `PROOF-ISOLATION`'s `(P3)`); and (ii) any *uniform* contraction is the marginal curvature `G″(0)=S−D` (`D2`/H1), while any cone *sharp* enough to split the atom from the `0⁺⁺` continuum is sign-indefinite and already assumes `f₁²>0`. The airy hope — "sign structure replaces dynamics" — dies because the sign structure alone gives no invariant cone, and the invariant one is the diffusion whose rate *is* the dynamics. **The cone-contraction ratio is a fourth face of H1, and a bottom face at that; #1,2,4 pass, the two fails (#3,#5) are the Clay content. No fabrication.**

*Refs.* Birkhoff, *Trans. AMS* **85** (1957) 219; Bushell, *Arch. Rational Mech. Anal.* **52** (1973) 330; Krasnosel'skii, *Positive Solutions of Operator Equations* (1964); Krein–Rutman 1948; Lemmens–Nussbaum, *Nonlinear Perron–Frobenius Theory* (CUP 2012); von Renesse–Sturm, *CPAM* **58** (2005) 923; Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737. Internal: PROOF-GAUSS-INVARIANT-PF (PF locks vacuum; positive kernel); PROOF-EDGE-E2 (bottom-functional blindness, I2); SUB2-I2-STRONG (`Be^{−δt}` vs `O(u)`); PROOF-ISOLATION (`(P3)`, operative `Λ≈1.5m_g`); INEQ-R1-perron-positive (Perron vector = vacuum; second cone absent); INVENT-SECOND-CONE (Karlin flag nested ⇒ free field interior); D2-COUPLING-DOEBLIN (`b′=−G″`, rate = H1); FORGE5-CONE-VARIANCE (nominal/vacuous #4); PROOF-F1SQ (gap ⇏ residue). No fabricated bound.
