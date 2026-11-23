# JONES-SUBFACTOR — the local inclusion A(O₁)⊂A(O₂) is type III₁ with **infinite, mass-independent** Jones index (free massless *and* massive alike); finite index is the DHR/chiral-sector datum that thrives in *massless* rational CFT, so the index is mass-blind and anti-correlated with the gap; pure YM has no charged sector to index

> **判词 (verdict).** Forge: treat `Ω⊂A(O)` (equivalently the nested inclusion `A(O₁)⊂A(O₂)`) as a subfactor inclusion and read `λ` off the **Jones index** `κ=[A(O₂):A(O₁)]` — finite index ⇒ gap, as for DHR particle sectors — i.e. `G″(k→0)=1/Var_ω(a₀) ≥ c/κ`. **Verdict: dead — a no-go, not a missing lemma, and sharper than the acid test requires.** The free *massless* acid test passes (infinite index ⇒ returns 0), but the tool **also returns 0 for the massive free field, which HAS a gap** — so it cannot certify even a *known* gap. Three exact walls. **(i) Mass-blind (§2):** `A(O)` is type III₁ at every mass (Fredenhagen), and `κ=∞` for the free field at `m=0` *and* `m>0` (relative commutant type III); no κ-functional sees `m`. **(ii) Anti-correlated (§3):** the sharpest finite-index objects — rational CFT / ADE-classified chiral nets — are exactly **gapless**, so "finite index ⇒ gap" is **false** (rational CFT refutes it, as the free massless net refutes AQFT-CAUSAL). **(iii) Empty referent (§4):** finite index measures **DHR/chiral sector content**; pure YM is superselection-complete (confinement, all glueballs in the vacuum sector), so there is no charged endomorphism to index and the only available inclusion is type III₁/∞. **Score 1/5 — no `λ`, no fabrication.**

## §1 The forge and the exact dictionary
`G″(k→0)=1/Var_ω(a₀)` (k0-CONSTRAINT-EFFPOT; MODULAR-TAKESAKI §1). Candidate: `G″≥c/κ`, `κ=[A(O₂):A(O₁)]`. Jones: `κ∈{4cos²(π/n):n≥3}∪[4,∞]`; finite κ carries the standard invariant / principal graph. The hope (imported from DHR–Buchholz–Fredenhagen): finite κ ↔ a superselection sector ↔ a particle ↔ a mass shell. The `k→0` mode `a₀` is a vacuum-sector object; the question is whether its variance is bounded by κ.

## §2 Mass-blind inclusion lemma (rigorous)
Fredenhagen: `A(O)` is type III₁ for the free field at *every* mass and for interacting nets. For `O₁⋐O₂` the relative commutant `A(O₁)′∩A(O₂)` is type III (infinitely many collar modes), so `κ=∞` — at `m=0` **and** `m>0`. The standard/split refinement (Doplicher–Longo) is the **split property**, a UV/nuclearity condition (AQFT-CAUSAL §4), mass-independent. Every κ-functional is thus invariant under `m→0`: it returns `0` for the free massless field (acid test met) **and equally `0` for the massive free field, which has a genuine gap** — it cannot *detect* a gap, hence cannot *certify* YM's.

## §3 Finite index ≠ gap — the rational-CFT refutation
The canonical finite-index structures — Jones' ADE classification (Cappelli–Itzykson–Zuber), Kawahigashi–Longo rational chiral nets — are **massless**. Chiral/CFT nets have finite κ yet `G″(k→0)=0` (no mass shell). So κ can be finite with zero gap: "κ<∞ ⇒ gap" is **false**, not unproven. The index counts chiral primaries / quantum dimensions — a *species* count, not a mass scale. It is tied to conformal covariance (massless regime), i.e. anti-correlated with the very gap we want.

## §4 Empty referent in pure YM
`κ<∞` is the DHR statistic index [Longo 1989/90] = square of the statistical dimension of a superselection endomorphism. Pure YM has **trivial superselection structure** (confinement): no localizable charged sectors; all asymptotic states are vacuum-sector glueballs. So there is no finite-index endomorphism to appeal to, and the vacuum-sector `A(O₁)⊂A(O₂)` is type III₁ / ∞ (§2). The one place finite index *does* imply a particle (DHR) presupposes the gap (particles ⇒ mass shell), so invoking it is H1 in subfactor dress.

## Five-point scorecard
1. **Nonperturbative — PASS.** Jones theory is exact, all β, no loops / `1/N`.
2. **Sign lower bound `m_eff²≥λ` — FAIL.** Mass-blind: `κ=∞` ⇒ `0` for free massless *and* massive; cannot lower-bound a gap it cannot see.
3. **Non-circular — FAIL.** The only finite-index datum is the DHR/chiral sector = the particle = H1; the vacuum-sector index is mass-free.
4. **Non-SUSY + fails free field — FAIL (refuted).** Free massless → 0 (good), but finite-index **massless** rational CFT refutes "κ<∞⇒gap"; no YM positivity.
5. **Uniform `≥λ` — FAIL.** No `λ`; `κ` is mass-independent.

## Exact failing estimate
```
Target:   G″(k→0)=1/Var_ω(a₀) ≥ λ        (k0-CONSTRAINT-EFFPOT; MODULAR-TAKESAKI §1)
Forge:    λ = c/κ,   κ=[A(O₂):A(O₁)]
(i)  A(O) type III₁ ∀m;  relative commutant type III ⇒ κ=∞  for m=0 AND m>0  ⇒ mass-blind ⇒ λ=0
(ii) rational/CFT nets: κ<∞  yet G″(k→0)=0  ⇒ "finite index ⇒ gap" FALSE  (counterexample)
(iii)pure YM: no DHR sectors (confinement) ⇒ no finite-index endomorphism ⇒ empty referent
missing:  a mass-sensitive functional of the standard invariant bounding 1/Var_ω(a₀).  Does not exist.
```

**Sharpest line.** The species is a no-go sharper than its own acid test: the Jones index returns `0` not only for the gapless free massless field (as required) but **also for the gapped massive free field** — it is blind to `m` — and the only finite-index theories we possess (rational CFT) are **massless**, so finite index is *anti-correlated* with the gap, not a certificate of it. In pure YM even the DHR harbour is empty: confinement erases the charged sectors whose statistical dimension the index would count. The index measures the *species* of the vacuum's superselection content, not the *stiffness* of its softest mode; the missing object is a mass-sensitive functional of the standard invariant — every known one is either infinite (typed out) or the particle it presupposes.

*Refs.* (internal) AQFT-CAUSAL (split/nuclearity = UV; type III₁; mass-blindness), MODULAR-TAKESAKI (§1 dictionary; BW), C8-SPECIES-NOGO (empty-referent pattern), k0-CONSTRAINT-EFFPOT, PROOF-IR-H1 (gap = clustering rate), COMPUTE-GPP-SU3 (`G″≈1/Var`), PROOF-GAUSS-INVARIANT-PF. External: Jones, Invent. Math. **72** (1983) 1; Doplicher–Longo, Invent. Math. **75** (1984) 493 (standard/split, type III index); Longo, CMP **126** (1989) 217 / **130** (1990) 285 (index = statistics²); Doplicher–Haag–Roberts, CMP **23** (1971) 199 / **35** (1974) 49; Fredenhagen, CMP **97** (1985) 79 (type III₁); Buchholz–Fredenhagen, CMP **84** (1982) 1 (particles ⇐ gap); Cappelli–Itzykson–Zuber, Nucl. Phys. B **280** (1987) 445 (ADE); Kawahigashi–Longo, Ann. Math. **160** (2004) 493 (rational nets); Haag, *Local Quantum Physics* (1996).
