# G2G-TRANSFER-SPECTRAL — rigorous `H` makes `Δ>0` well-posed and RP+PF(Bushell) gives the first right-sign *lower* bound; uniform `Δ>0` and second-level isolation both reduce to `H1`

> **2026-11-23 · G2G (transfer-spectral).** 判词 (verdict). Grant the OpenAI/Lean-2026 premise: `H≥0` self-adjoint on `H_phys=ker Gauss`, `T=e^{−aH}`. This buys exactly three things over the lattice proxies: **(i)** "the mass gap of the 0⁺⁺ sector" is now a statement about a genuine self-adjoint operator, not a finite-`V` proxy; **(ii)** reflection positivity (`T≥0`) + Perron–Frobenius (positivity-improving) give `λ₀=1` simple with strictly-positive vacuum, and the *quantitative* Bushell/Krein–Rutman cone-contraction yields the project's first genuinely **lower-bounding, right-sign** spectral instrument — a lower bound on the gap `E₁−E₀`, not a susceptibility; **(iii)** the Källén–Lehmann measure `ν` exists, so `f₁²` is a well-defined real number. **But (i)–(iii) are the OBJECT, not the gap** (G1-CONTINUUM-OS: existence of `ν` is a tautology of the spectral theorem once `H` exists). The *uniform* sufficient condition `Δ:=inf_{V,a}(E₁−E₀)>0` is equivalent to the Bushell contraction staying `c_V<1` uniformly, and `c_V→1` unless **`G″(0)=S−D≥λ` = H1**; the *second-level* isolation (`E₂−E₁`, the atom `f₁²`) is provably not a bottom functional (PROOF-EDGE-E2) and is unreachable. **#1,2,4 PASS; #3,5 FAIL = the Clay content. No closure; 0%.**

## §1 · What rigorous `H` buys that lattice proxies did not
Lattice proxies bound the gap from *above* (GEVP/Courant–Fischer `≥E_k`), measure susceptibilities `χ=1/G″`, or fit `m_eff(t)` (an approach from above). None is a *lower* bound with the right sign, and none can cleanly *state* `Δ>0`: finite-`V` gaps are upper bounds on the mass gap and the `V→∞,a→0` limits were heuristic. The rigorous `H` converts this into a well-posed spectral problem — `spec(H|_{0⁺⁺})`, `λ_n=e^{−aE_n}` — and separates the walls: **L3** (UV existence) from **H1** (IR atom) (G1-CONTINUUM-OS). That is the buy — **well-posedness plus a legitimate RP+PF operator**, not a gap.

## §2 · Two readings of "second eigenvalue bounded away from the first"
**(A) The gap.** `λ₀−λ₁ = 1−e^{−aE₁} ≥ 1−e^{−aΔ}` uniformly `⇒ E₁≥Δ>0`: no massless `0⁺⁺`. **(B) The isolation.** `λ₁−λ₂ = e^{−aE₁}−e^{−aE₂} ≥` … uniformly `⇒` the lightest `0⁺⁺` is an isolated atom above the two-particle continuum, i.e. `f₁²>0` (`I2`). (A) is the sufficient condition for "gapped"; (B) is the sufficient condition for "a particle." **Both are attacked below; both land on `H1`.**

## §3 · RP + PF bottom, Bushell the gap — the genuine lower bound
`T̄` on `𝒜/G` is positivity-improving (positive Haar pushforward, no Gribov/FP sign — PROOF-GAUSS-INVARIANT-PF §1). Jentzsch/Krein–Rutman: `λ₀=1` simple, `Ω≡1` strictly positive, `0≤λ₁<λ₀`. **Quantitatively** (F2 §1): Birkhoff `c(T)=tanh(Δ(T)/4)` and Bushell `|λ|≤cλ₀` for `λ≠λ₀` give, at each finite `V`, `λ₁/λ₀ ≤ c(T)<1 ⇒ Δ_V := E₁−E₀ ≥ −(1/a)·ln c(T)>0`. A real *lower* bound on the `0⁺⁺` gap — the one species not sign-blind.

## §4 · The rate is `G″`; isolation is second-level, not bottom
**(a)** The uniform contraction rate is the marginal curvature: on the coarse block generator `b(a)=G″a`, `−b′(0)=G″=S−D`, so `sup_V c_V<1 ⟺ G″(0)≥λ ⟺ H1` (F2 §3, D2; Ward `G″·Var(a)=1` ties it to `χ` — so bounding the gap "not just `χ`" lands on the *same* number in another vocabulary). **(b)** Reading (B): the physical glueball needs `E₂−E₁=Λ−E₁≥δ`. Every instrument — Temple/Kato–Weinstein, OPE-from-above, GEVP, the block-local Poincaré gap, Bushell itself — is a bottom-Dirichlet functional encoding only `λ₁`, **blind to `λ₁−λ₂`** (PROOF-EDGE-E2): `{1,λ₁,λ₁,…}` and `{1,λ₁,λ₂,…}` share every bottom datum. Temple on `(vac⊕glueball)^⊥` needs an `E₃` bound — one rung up, a structural regress. "Rigorous isolation of the second eigenvalue" is **not a tool we possess**; the missing estimate is a two-mass bound, which is `f₁²>0`, which is `H*=H1`.

## §5 · Free-field acid test
Free massless Maxwell: `H` has continuous spectrum from 0 (`k→0` flat ⇒ `G″=0` ⇒ `b≡0` ⇒ `c=1` ⇒ `Δ_V→0`), and `V→∞` gives `c_V→1`. **Acid test PASSED** — the bound correctly gives *no* gap for the gapless free field, and is non-SUSY. But the degeneration is universal: the YM kernel is equally spread at `∞V`, so `c_V→1` there too **unless H1 restores a uniform drift**. The pass is exact but *vacuous for YM* absent `G″>0`.

## Five-point scorecard
1. **Nonperturbative — PASS.** Exact self-adjoint `H`, all β; Bushell on the cone, no `1/N`/`g`.
2. **Sign-controlling lower bound on the gap `E₁−E₀` — PASS (shape/sign).** Genuine `Δ_V>0`, right sign — the gain of rigorous `H`; bounds the *gap*, not `χ`. (One functional short of the atom, §4.)
3. **Non-circular — FAIL.** Uniform rate `=G″=H1`; second-level isolation presupposes `f₁²>0`.
4. **Non-SUSY + fails free field — PASS (clean).** `G″=0` ⇒ `c=1` ⇒ `0`; continuous spectrum from 0 correctly unbounded.
5. **Uniform `Δ≥λ` in `(V,a)` at β≈2.5 — FAIL.** `c_V→1`; `c₀<1 ⟺ G″(0)≥λ`.

## Exact failing estimate
```
RP: T≥0 ; PF: λ₀=1 simple, Ω≡1>0 ; KR: λ₁<λ₀             [bottom — proven]
Bushell: λ₁/λ₀ ≤ tanh(Δ/4)<1 ⇒ E₁−E₀ ≥ −(1/a)ln c         [lower bound, right sign]
BUT  coarse generator drift b(a)=G″a ⇒ rate = G″ = S−D    [D2]
  ⇒ sup_V c_V<1 ⟺ G″(0)=S−D ≥ λ > 0 ⟺ H1 .               [the wall]
AND  atom needs E₂−E₁=Λ−E₁≥δ, i.e. λ₁−λ₂>0 :
  Temple/OPE/GEVP/Poincaré/Bushell all bottom ⇒ blind to λ₁−λ₂   [EDGE-E2]
  ⇒ missing estimate = two-mass isolation I2 (f₁²≥r₀) = H* = H1 .
free field: G″=0 ⇒ c=1, Δ_V→0 ; acid PASSED (vacuous for YM without H1).
```

**Sharpest line.** The rigorous existence of `H` is the precondition, not the content: it makes `Δ>0` a well-posed question about a real self-adjoint operator and lets RP+PF+Bushell deliver the first *lower*-bounding, right-sign instrument the project has — a genuine advance over sign-blind susceptibilities and above-only GEVP. But the Bushell contraction rate *is* the marginal curvature `G″=S−D`, so "bound the gap, not the susceptibility" lands on the same unproven number (Ward); and the *isolation* of the second eigenvalue is a two-level (`λ₁−λ₂`) datum no bottom functional can see, so "rigorous isolation of the second eigenvalue" is precisely the tool we lack. **Reading (A) uniform `Δ>0` ⟺ H1; Reading (B) isolation ⟺ `H*`=H1. 3/5, not 5/5; Clay 0%.**

*Refs.* Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Lüscher, *CMP* **54** (1977) 283; Birkhoff, *Trans. AMS* **85** (1957) 219; Bushell, *Arch. Rational Mech. Anal.* **52** (1973) 330; Krein–Rutman 1948; von Renesse–Sturm, *CPAM* **58** (2005) 923. Internal: G1-CONTINUUM-OS (object ≠ gap; L3/H1 split); PROOF-GAUSS-INVARIANT-PF (positivity-improving, RK vacuum+gap, edge obstruction); F2-KREIN-RUTMAN-QUANT (Bushell lower bound = `G″` = H1); PROOF-EDGE-E2 (bottom-functional blindness, Temple regress, two-mass I2); G2-ASSEMBLE-100 (G2 modulo `{H1}`, Clay 0%); D2-COUPLING-DOEBLIN; GAP-G2-edge (residue–level-gap lemma); mass-spectrum (Thm A/B, I1/I2). No fabricated bound.
