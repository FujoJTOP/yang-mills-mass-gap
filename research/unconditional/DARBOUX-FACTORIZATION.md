# DARBOUX-FACTORIZATION — the Darboux/Infeld–Hull chain around the exact Perron vacuum is a *supersymmetric isospectral* rearrangement: it yields the gap only when shape invariance closes, and every available closure is either blocked by the `0⁺⁺` continuum (global) or is the Gaussian block whose shape-invariance constant *is* `G″` (local) = H1. Sign-blind, exactly as BHH. **#1,#4 PASS; #2,#3,#5 = the Clay content.**

> **2026-11-23 · NEW SPECIES — Darboux transform / Infeld–Hull shape invariance (pure spectral device; the physical theory stays non-SUSY YM).** 判词 (verdict). Factor the `0⁺⁺` transfer generator as `H−E₀=A†A`, `A=∂_a+W`, `W=−∂_aφ₀/φ₀`. Because the vacuum `φ₀` is the *exact nodeless Perron vector* (`PROOF-GAUSS-INVARIANT-PF`), `W` is real globally on finite volume — the factorization exists, unlike in the massless continuum. But the Darboux transform is an **isospectral identity**: it says `spec(AA†+E₀)=spec(H)\{E₀}`, i.e. the gap `m₁` *is* the ground energy of the partner `H₂`. That is an equality, sign-blind — it gives `m₁` only when **shape invariance** `V₂(a;α)=V₁(a;α₁)+R(α)` closes the chain, and then `m₁=R(α)`. The two available closures both collapse: the *global* chain is impossible (the `0⁺⁺` spectrum contains a two-gluon continuum above threshold, §3i), and the *local/coarse block* is the Ornstein–Uhlenbeck harmonic oscillator, whose shape-invariance constant is the drift slope `R=−b′(0)=G″(0)=S−D` — **verbatim H1** (§3ii). **#1,#4 pass; #2,#3,#5 are the entire Clay content. No closure.**

## §1 What the Darboux chain gives — an isospectral equality, not a bound
`H=A†A+E₀`, `Aφ₀=0`. Then `Aφ_n` are eigenvectors of `H₂=AA†+E₀` with the *same* eigenvalues `E_n`, `n≥1`; so `spec(H₂)={E_n:n≥1}` and `m₁=E₁−E₀=inf spec(H₂)−E₀`. The whole gap is the **ground energy of the partner potential** `V₂=W²+W′+E₀`. This is a re-arrangement (Witten SUSY-QM), not an estimate: `A†A≥0` already just restates `E₀=inf spec`, and bounding `inf spec(H₂)` to get `m₁` from below needs a *lower* bound on a *new* bottom — the Temple regress one rung up (`TEMPLE-WEINSTEIN` §2). The sign of `R` is uncontrolled by the transform itself.

## §2 Shape invariance is the only algebraic handle — and `m₁=R`
If `V₂(a;α)=V₁(a;α₁)+R(α)`, Gendenshtein's chain gives `E_n=E₀+Σ_{k<n}R(α_k)` and `m₁=R(α₀)`. Solvability is thus *definitionally* the gap written as a potential-difference constant. Infeld–Hull shape invariance is a measure-zero structure on potentials (Pöschl–Teller, Morse, Coulomb, oscillator families); nothing in the gauge action selects it. Asserting it = asserting the transfer operator is exactly solvable = asserting the spectrum, i.e. assuming the answer.

## §3 The two closures both die
**(i) Global — blocked by the continuum.** An exact shape-invariant chain yields a pure discrete tower. The `0⁺⁺` gauge spectrum carries a genuine **two-gluon cut** above the lightest glueball — the same bottom/continuum obstruction as `PROOF-EDGE-E2`/`F2` §2. So no exact global Darboux chain exists; shape invariance can at most be an *effective below-threshold* statement.
**(ii) Local block — it is the oscillator, and `R=G″`.** The one coarse generator acting on the radial mode is the block diffusion with drift `b(a)=G″a` (`F2` §3; `D2-COUPLING-DOEBLIN`). Its Schrödinger form is the harmonic oscillator — the *archetype* of shape invariance — with gap `R=−b′(0)=G″(0)=S−D`. So "the gauge Hamiltonian is block-factorizable ⇒ the gap follows algebraically" reads, concretely, `m₁=G″(0)`, which is **H1**. And like every SUSY/Ward structure this factorization is **sign-blind** (`BHH-TRANSLATION`, `C8-SPECIES-NOGO`): the only surviving exact content is the divergence/Stein identity `E[G′]−Var=G″`.

## §4 Free-field acid test + uniformity
Free massless: `W=0`, `V₂=V₁+0`, trivially shape-invariant with `R=0` ⇒ gap `0`; the two-photon cut at `s=0` is the trivial fixed point. **Acid test PASSED** (the method correctly reports no gap; non-SUSY physical theory, `#4` clean). *But vacuously* (as F2/BHH): `R=0` for free and, absent `G″>0`, `R=0` for YM too — the discriminator is exactly `G″`. Uniform `(V,a)` FAILS: the block `R=G″(0)` is uniform `⟺ G″(0)≥λ`, i.e. H1 unproven.

## Five-point scorecard
1. **Nonperturbative — PASS.** Exact Darboux around the nodeless Perron vacuum, all β; no `1/N`.
2. **Sign-controlling `m₁≥λ` — FAIL.** Darboux is an isospectral *equality*, sign-blind; `R`'s sign uncontrolled (the BHH face).
3. **Non-circular — FAIL.** Algebraic gap needs shape invariance = exact solvability = the spectrum assumed; block closure `R=G″=H1`.
4. **Non-SUSY + fails free field — PASS.** Physical theory non-SUSY; free shape-invariant `R=0` ⇒ no gap (vacuous-uniform caveat).
5. **Uniform `≥λ` — FAIL.** `R=G″(0)`; uniform `⟺ G″≥λ` = H1.

## Exact failing estimate
```
Darboux:  H=E₀+A†A, A=∂_a+W ⇒ spec(AA†+E₀)=spec(H)\{E₀} ⇒ m₁ = inf spec(H₂)          [isospectral — equality]
shape inv.: V₂(α)=V₁(α₁)+R(α) ⇒ E_n=E₀+Σ R(α_k), m₁=R(α₀)                              [only algebraic handle]
global:  exact chain ⇒ pure discrete tower ; 0⁺⁺ has 2-gluon cut ⇒ NO exact chain      [EDGE-E2 continuum]
block:   b(a)=G″a ⇒ Schrödinger form = oscillator (shape-invariant) ⇒ R=−b′(0)=G″=S−D  [F2 §3; D2]
  ⇒ "block-factorizable gap" ≡ m₁=G″(0) ≡ H1 ; sign-blind (Stein E[G′]−Var=G″)          [BHH; C8]
free:    W=0 ⇒ R=0 ⇒ gap 0                                                              [trivial fixed point]
exact missing:  shape-invariant closure with independent R≥λ ; block closure R=G″ ⇒ G″(0)≥λ (H1)
```

**Sharpest line.** Darboux/Infeld–Hull is the *supersymmetric face* of the gap: factor `H−E₀=A†A` around the exact Perron vacuum and the gap is the partner's ground — an isospectral equality, structurally the same sign-blindness `BHH-TRANSLATION` already exposed. "Algebraic gap from shape invariance" is honest only if the chain closes; but a closed chain forces a discrete tower the `0⁺⁺` continuum forbids (global), and the only closure the coarse slice offers is the harmonic oscillator whose shape-invariance constant `R` *is* the marginal curvature `G″(0)=S−D` (local). So the factorization genuinely reorganizes the `0⁺⁺` spectrum around the vacuum — the one new asset — and then hands back the same number H1 has been carrying all along. `#1,#4` pass; `#2,#3,#5` are H1. **No fabrication; the wall is relabeled, not crossed.**

*Refs.* G. Darboux, *C. R. Acad. Sci. Paris* **94** (1882) 1343; E. Schrödinger, *Proc. Roy. Irish Acad.* **A46** (1940) 9; L. Infeld, T. E. Hull, *Rev. Mod. Phys.* **23** (1951) 21; L. Gendenshtein, *JETP Lett.* **38** (1983) 356; E. Witten, *Nucl. Phys. B* **185** (1981) 513; F. Cooper, A. Khare, U. Sukhatme, *Phys. Rep.* **251** (1995) 267. Internal: PROOF-GAUSS-INVARIANT-PF (Perron vacuum nodeless ⇒ real `W`); TEMPLE-WEINSTEIN (§2 partner/`E₁` regress); F2-KREIN-RUTMAN-QUANT (§3 drift `b=G″a`); D2-COUPLING-DOEBLIN (`b′=−G″`); BHH-TRANSLATION (SUSY/Ward sign-blind); C8-SPECIES-NOGO (only `E[G′]−Var=G″` survives); PROOF-EDGE-E2 (`0⁺⁺` two-gluon cut; bottom-functional blindness); PROOF-SAVVIDY-EXACT-SIGN (§3ii block curvature). No fabricated bound.
