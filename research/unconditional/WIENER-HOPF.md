# WIENER-HOPF — the gap as a finite-width analyticity strip of the 0⁺⁺ resolvent: Wiener–Hopf controls the *index* of a symbol, not the *width* of the strip; in ≥2D there is no symbol to factorize, and scattering theory presupposes the mass shell it would produce

> **2026-11-23 · NEW SPECIES — Wiener–Hopf / scattering.** 判词 (verdict). Forge: `G″(k→0)≥λ` ⟺ a `0⁺⁺` resolvent analytic in a strip of width `√λ` ⟺ no weight in `(0,λ)`; factorize the collective-mode symbol `S=S_+S_−` and read the strip off it. **Verdict: shape-correct, acid-test-clean, content-empty.** WH controls the *index/topology* of a symbol on a 1-D contour; the gap is the *width* of the analyticity strip — real, non-topological, never computed by WH. Worse: in ≥2D there is no symbol factorization (the quarter-plane problem reduces to a Fredholm equation), and in 3+1D **Haag–Ruelle scattering assumes the isolated mass shell — the gap — as a hypothesis**, so "scattering states ⇒ gap" runs backward. The one real theorem is 1+1D integrable TBA, where integrability *furnishes* the `S`-matrix; 4D SU(N) has none. **#1 PASS (character); #2,#3,#5 FAIL; #4 PASS (vacuous). No `λ`.**

## §1 The exact dictionary — target `=` strip width `= H1`
KL: `R(z)=∫dν(s)/(s−z)`, `ν≥0`, `supp ν⊂[m₁²,∞)`, atom `f₁²=ν({m₁²})` (PROOF-DISSOLUTION). `R` is analytic on `ℂ∖supp ν`; analytic in a strip *just below* `m₁²` **iff** `m₁²=inf supp ν`. So "analytic in a strip" = the gap, an *identity*, not an estimate. For the target: `G″(k→0)` is the collective-mode stiffness conjugate to `Var(a)` (`G″=1/Var`, D2/F2); `C_a(t)` analytic in `|Im t|<1/m_a` makes `G″≥λ` exactly the width statement `m_a≥√λ` — verbatim H1. WH factorization `S=S_+S_−` yields the *index* (zeros on the contour; a zero *on* it is the threshold). Index is **integer, protected, topological**; gap is **real, unprotected, analytic** — wrong resolution.

## §2 Haag–Ruelle circularity — scattering presupposes the shell
Asymptotic scattering states (LSZ; Haag–Ruelle) require the spectrum condition **plus an isolated one-particle mass shell above the vacuum** — i.e. the gap — as a *hypothesis*. Free massless: no isolated shell, no such states, "scattering states absent" — but as a *consequence* of no gap, not its cause. Below `s_th=4m₁²` there are by definition no multiparticle scattering states, only the atom (PROOF-DISSOLUTION: the lightest `0⁺⁺` sits at `m₁²`, a quarter of `s_th`). "Absence of scattering states below threshold" is a kinematic tautology; the gap is the *atom*, which scattering takes as input.

## §3 No symbol in ≥2D — the quarter-plane obstruction
WH is 1+1D by construction (one contour). Lifting to a 3D-momentum symbol is the multidimensional / quarter-plane problem: no closed-form factorization, reduces to a Fredholm equation; matrix symbols (color×spin) have unstable partial indices (Gohberg–Krein canonical factorization generically fails). Granting a factorization, the 3+1D object is a 4D eigenvalue problem — `inf_{T³} log λ_max` — whose answer is the gap; WH adds nothing.

## §4 Free-field acid test — passes, vacuously
Free massless: `ν` continuous from `0`, branch point on the contour, strip width `0`, factorization degenerate, gap `0`. The species **correctly fails the free field** (claims no spurious gap); gapped free field gives strip `=m`. It discriminates — but only by *reading the strip off a spectrum it cannot obtain without assuming it* (F2's vacuous pass, AQFT-CAUSAL's mass-blindness).

## §5 What is genuinely true
1+1D integrable (Sinh–Gordon, O(3) NLσM): integrability gives the exact `S`-matrix; TBA/WH *derives* the gap nonperturbatively (Hasenfratz–Niedermayer; Yang–Yang). 2D Ising: Toeplitz determinant / Szegő gives the exact edge from a symbol. Real — and dimension-/integrability-locked; neither available for 4D SU(N).

## Five-point scorecard
1. **Nonperturbative — PASS (character).** No loop/`1/N`; but concretely only 1+1D integrable; 4D `S`-matrix nonperturbatively unknown.
2. **Sign-controlling LOWER bound — FAIL.** Right *shape* (lower, unlike moment/C-S), but in 4D no factorization (quarter-plane) and "analytic in a strip" is the hypothesis restated; no `δ` produced.
3. **Non-circular — FAIL.** Haag–Ruelle *assumes* the isolated shell; strip width `=inf supp ν=` gap; WH of the `a`-mode symbol `=` H1.
4. **Non-SUSY + fails free field — PASS (vacuous).** Free massless: zero strip ⇒ gap `0`, correctly; runs only after the shell is assumed.
5. **Uniform `≥λ` — FAIL.** Index integer, width `=` gap, quarter-plane open; no uniform `(V,a,β)` `λ`.

## Exact failing estimate
```
Target: G″(k→0)≥λ  ⟺  R|0⁺⁺ analytic in strip of width √λ  ⟺  Var(a)≤1/λ   [D2; F2; H1]
spectral thm: "analytic in strip below m₁²" ⟺ m₁²=inf supp ν       [identity, not estimate]
WH: S(θ)=S₊S₋ ⇒ INDEX (integer, protected); zero on contour = threshold  [topology]
     gap = real WIDTH of strip (unprotected)  ⇒ wrong resolution, no width computed
Haag–Ruelle/LSZ: isolated mass shell (the gap) assumed to build scattering states ⇒ circular
≥2D: no symbol factorization (quarter-plane ⇒ Fredholm; matrix partial-index unstable)
4D: inf_{T³} log λ_max = the gap itself; WH adds nothing.
free massless: ν from 0 ⇒ strip 0 ⇒ gap 0 (acid test PASS) — but vacuous (reads strip off spectrum).
```

**Sharpest line.** Wiener–Hopf and scattering are the sharpest *qualitative* tools — poles, thresholds, bound-state counting — and the weakest *quantitative*: WH computes the *index* of a symbol while the gap is the *width* of the analyticity strip, a real number the topological machinery never touches. The engine is circular from the first line: Haag–Ruelle scattering states *require* the isolated mass shell (the gap) before they exist, so they cannot be its source. The one honest theorem — 1+1D integrable TBA — works precisely because integrability supplies the `S`-matrix 4D SU(N) withholds, deriving its gap *from a datum that does not exist here*. Same wall, new alphabet: the missing object is a nonperturbative `0⁺⁺` `S`-matrix (the uniform strip width), which *is* `G″(k→0)≥λ` `=` H1.

*Refs.* (external) Wiener–Hopf 1931; Noble, *Methods Based on the Wiener–Hopf Technique* (1958); Böttcher–Silbermann, *Analysis of Toeplitz Operators* (2nd ed., Springer 2006); Gohberg–Krein, *Systems of Integral Equations…* (1960, canonical matrix factorization); Lehmann–Symanzik–Zimmermann, Nuovo Cim. **1** (1955) 205; Haag–Ruelle, Helv. Phys. Acta **35** (1962) 290; Haag–Swieca, CMP **1** (1965) 308; Yang–Yang, Phys. Rev. **170** (1968) 1608; Al. Zamolodchikov, *Thermodynamic Bethe Ansatz* (Nucl. Phys. B **342** (1990) 695); Hasenfratz–Niedermayer, Nucl. Phys. B **337** (1990) 233; Szegő 1952; McCoy–Wu, *The Two-Dimensional Ising Model* (1973). (internal) PROOF-DISSOLUTION-SPECTRAL-DENSITY (`ν`, `s_th=4m₁²`, atom `f₁²`); GAP-G2-edge (moment blindness, upper-only); F2-KREIN-RUTMAN-QUANT (bottom functional, rate `=G″=H1`, vacuous acid-test pass); AQFT-CAUSAL (mass-blindness, clustering rate `=` H1); C8-SPECIES-NOGO (empty-referent / false-implication pattern); D2-COUPLING-DOEBLIN (`b′=G″`, `G″=1/Var`); AREALAW-GAP (shape-pass, channel-fail template); PROOF-P3-UNCONDITIONAL (`β<1/48`). No fabricated bound.
