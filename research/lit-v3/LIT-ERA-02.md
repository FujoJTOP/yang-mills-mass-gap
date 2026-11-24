# LIT-ERA-02 · 1963–1973 — Quantization, Renormalization, Asymptotic Freedom

> 文献基座 for paper v3. Era 02 of 14. Scope: covariant quantization of Yang–Mills, proof of perturbative renormalizability, the renormalization group and dimensional transmutation, asymptotic freedom, and the first appearance of the infrared problem. BRST (1975–76) is included as the immediate continuation because it is the operator-algebraic closure of what Faddeev–Popov and 't Hooft started.
>
> **Honest line:** this era proved *perturbative renormalizability* and discovered *asymptotic freedom*. It proved **nothing** about the infrared, about confinement, or about a mass gap. Every statement below is labelled [T] theorem / [C] conditional / [NG] no-go / [O] exact missing step / [V] numerical. Nothing in this file claims the Clay gap.

---

## 0. What this era actually settled (one paragraph)

By 1973 the ultraviolet problem of non-Abelian gauge theory was solved *order by order in perturbation theory*: Faddeev–Popov (1967) gave the ghost determinant; DeWitt (1967) gave the background-field/closed-loop rules; 't Hooft (1971) and 't Hooft–Veltman (1972) proved renormalizability with dimensional regularization; Taylor (1971) and Slavnov (1972) supplied the Ward identities that make a single renormalized coupling consistent; and Gross–Wilczek/Politzer (1973) computed the one-loop β-function and found it negative for pure SU(N). The same computation made the coupling grow in the infrared and made the perturbative vacuum unstable there. **That infrared side is exactly where the mass gap lives, and it was left untouched.** Renormalizability is a UV statement; asymptotic freedom is a UV statement with an IR corollary (Landau-pole-like growth) that perturbation theory cannot follow. The era ends with the problem relocated, not solved.

---

## 1. Catalogue of entries

Group A — covariant quantization and the ghost (1963–1970)
Group B — renormalizability proved (1971–72)
Group C — renormalization group and dimensional transmutation (1954–1973)
Group D — asymptotic freedom (1965–1974)
Group E — the infrared problem begins (1962–1978), plus BRST/Tyutin continuation

---

### A. Covariant quantization and the ghost

#### E01 · Feynman (1963) — closed-loop rule [T]
**R. P. Feynman**, "Closed loop diagrams," *Acta Phys. Pol.* **24** (1963) 697.
Status: **UNVERIFIED** (standard citation; search budget exhausted before checking; exact pages not independently confirmed this session).
What it establishes: the first statement that closed-loop diagrams in a gauge-invariant theory require an extra rule beyond tree Feynman rules; gave a recipe for the one-loop case. Faddeev–Popov cite this as the point of departure.
What it does not: no general prescription, no ghost field, no renormalizability.
Role in v3: pre-history only. The "extra rule" is the first symptom that the gauge orbit is non-trivial.

#### E02 · DeWitt (1967) — background field and fictitious quanta [T]
**B. S. DeWitt**, "Quantum theory of gravity. II. The manifestly covariant theory," *Phys. Rev.* **162** (1967) 1195–1239. DOI 10.1103/PhysRev.162.1195. (Part III, "Applications of the covariant theory," ibid. 162, 1239–1256.)
Status: **VERIFIED** (APS/Google Scholar lookup; abstract confirms Yang–Mills treated simultaneously with gravity; "fictitious quanta" appear only in closed loops, vertices vanish for Abelian group, uniform orientation around any fictitious-quantum loop).
What it establishes: independent of Faddeev–Popov, the background-field method and the closed-loop determinant ("fictitious quanta" = ghosts) for Yang–Mills and gravity; a general algorithm for primary diagrams to arbitrary order including all ghost loops; group invariance proved.
What it does not: no renormalizability proof; the determinant is formal.
Role in v3: the geometric/background-field formulation that v3's gauge-invariant constructions (W-entropy, gauge-invariant IR mass) inherit. DeWitt's "fictitious quanta" and Faddeev–Popov's "ghosts" are the same object discovered twice.

#### E03 · Faddeev & Popov (1967) — the ghost determinant [T]
**L. D. Faddeev & V. N. Popov**, "Feynman diagrams for the Yang–Mills field," *Phys. Lett. B* **25** (1967) 29–30. DOI 10.1016/0370-2693(67)90067-6.
Status: **VERIFIED** (CERN CDS record 431589; INSPIRE ~2,530 citations; abstract confirms functional-integral method and ghost fields).
What it establishes: the Faddeev–Popov trick — insert 1 = Δ_FP[A] ∫Dg δ(G(A^g)) into the path integral; Δ_FP = det(δG/δω) becomes a Grassmann determinant, i.e. anticommuting scalar ghosts in closed loops. This is the operator that makes the Yang–Mills path integral well-defined *perturbatively*.
What it does not: the determinant is only formally defined; Δ_FP can vanish (zero modes) and is not globally defined (see E31, Gribov 1978). The paper is a Letter — no renormalizability claim.
Role in v3: Δ_FP = det(−D²·δ + …) is the operator whose spectrum v3 studies. Its zero modes / small eigenvalues in the IR are the perturbative shadow of the mass gap.

#### E04 · Veltman (1968) — massive Yang–Mills, one loop [T]
**M. J. G. Veltman**, "Perturbation theory of massive Yang–Mills fields," *Nucl. Phys. B* **7** (1968) 637–650. DOI 10.1016/0550-3213(68)90197-1.
Status: **VERIFIED** (ADS 1968NuPhB...7..637V; INSPIRE 52788).
What it establishes: using the Bell–Treiman transformation, primitive one-loop diagrams with more than four external vector lines are convergent; the analysis "does not exclude" renormalizability. This is Veltman's route to the 1971–72 result.
What it does not: no proof; the four-point and lower sectors remain divergent; no Ward identity yet.
Role in v3: shows the divergence structure is milder than power counting suggests — the first hint that the non-Abelian structure cancels UV.

#### E05 · Khriplovich (1969) — charge renormalization in Coulomb gauge [T] (misread)
**I. B. Khriplovich**, "Green's functions in theories with Yang–Mills fields," *Sov. J. Nucl. Phys.* **10** (1969) 409 [*Yad. Fiz.* **10** (1969)].
Status: **UNVERIFIED** (content and interpretation confirmed via 't Hooft 1998, E22, which states Khriplovich "correctly calculated the charge renormalization of Yang–Mills theories in the Coulomb gauge… found the unusual sign, but a connection with asymptotic freedom was not made"; volume/page from standard citation, not independently checked this session).
What it establishes: the negative sign of the Yang–Mills charge renormalization was computed in 1969 and went unnoticed.
What it does not: no interpretation, no β-function language, no asymptotic freedom.
Role in v3: a [O]-adjacent historical datum — the result existed four years before 1973 and was not connected to the IR. Cite carefully: 't Hooft 1998 is the primary secondary source.

#### E06 · Fradkin & Tyutin (1970) — S-matrix, unitarity, gauge independence [T]
**E. S. Fradkin & I. V. Tyutin**, "S matrix for Yang–Mills and gravitational fields," *Phys. Rev. D* **2** (1970) 2841–2857. DOI 10.1103/PhysRevD.2.2841.
Status: **VERIFIED** (ADS 1970PhRvD...2.2841F; INSPIRE ~255 citations).
What it establishes: a generating-functional construction for fields with an invariance group, applied to Yang–Mills and gravity; unitarity and gauge independence of the S-matrix demonstrated on shell. A direct predecessor of BRST.
What it does not: no renormalizability proof; no cohomological formulation.
Role in v3: the first clean statement that the physical S-matrix is independent of the gauge-fixing — the property v3 needs if it works with a gauge-fixed measure and then quotients.

---

### B. Renormalizability proved (1971–72)

#### E07 · 't Hooft (1971a) — massless Yang–Mills renormalizable [T]
**G. 't Hooft**, "Renormalization of massless Yang–Mills fields," *Nucl. Phys. B* **33** (1971) 173–199. DOI 10.1016/0550-3213(71)90395-6.
Status: **VERIFIED** (ADS 1971NuPhB..33..173T; abstract confirms: non-gauge-invariant regulators admissible if gauge invariance restored in the high-regulator-mass limit by finitely many counterterms; massless YM renormalizable; consistency proved for non-overlapping divergences; assuming consistency, S-matrix unitary to all orders; restriction: group must not contain local parity-changing transformations).
What it establishes: **the central theorem of the era** — massless Yang–Mills is renormalizable in the usual sense, *provided* a gauge-invariant regularization exists and the overlapping divergences are consistent. The parity restriction is explicit.
What it does not: the consistency of overlapping divergences is **assumed**, not proved; no gauge-invariant regulator is exhibited (that is E11); no statement about the IR or the spectrum.
Role in v3: state precisely as "renormalizable order by order, conditional on a gauge-invariant regulator" — the 1971 paper alone is *not* the unconditional theorem; E11 completes it. This is the distinction v3 must not blur.

#### E08 · 't Hooft (1971b) — massive Yang–Mills via spontaneous breaking [T]
**G. 't Hooft**, "Renormalizable Lagrangians for massive Yang–Mills fields," *Nucl. Phys. B* **35** (1971) 167–188. DOI 10.1016/0550-3213(71)90139-8.
Status: **VERIFIED** (Scopus 0000043313; OSTI 4711793; abstract confirms path-integral Feynman rules and Ward identities via spontaneous breaking; SU(2) → U(1) model renormalizable and unitary; SU(2)⊗U(1) model for ρ–γ).
What it establishes: spontaneous breaking of a local gauge symmetry preserves renormalizability and unitarity; the Higgs mechanism is consistent at the quantum level. Companion to E07.
What it does not: the breaking is put in by hand (a scalar sector); no dynamical breaking; no mass gap in pure YM.
Role in v3: the template for "mass without a mass term." v3's target is mass *without* a Higgs field — this paper is the contrast case, not the mechanism.

#### E09 · Taylor (1971) — generalized Ward identity [T]
**J. C. Taylor**, "Ward identities and charge renormalization of the Yang–Mills field," *Nucl. Phys. B* **33** (1971) 436–444. DOI 10.1016/0550-3213(71)90297-5.
Status: **VERIFIED** (MR0314386; abstract confirms Z₁/Z₂ common value for vector and fictitious scalar, hence equal renormalized charges; identity avoids quadratic-divergence ambiguities).
What it establishes: the Slavnov–Taylor identity in its first form — Z₁/Z₂ is common to the vector and ghost, so there is **one** renormalized coupling. This is what closes the gap left by E07.
What it does not: proved for the massless case; no IR statement.
Role in v3: the exact identity that constrains the gluon propagator's longitudinal part — the handle for any rigorous statement about the IR two-point function.

#### E10 · Slavnov (1972) — Ward identities, general form [T]
**A. A. Slavnov**, "Ward identities in gauge theories," *Theor. Math. Phys.* **10** (1972) 99–104 [*Teor. Mat. Fiz.* **10**:2 (1972) 153–161]. DOI 10.1007/BF01090719.
Status: **VERIFIED** (ADS 1972TMP....10...99S; Springer; Scopus 485 citations; pages 99–104 confirmed, one Scholar entry's 99–107 is wrong).
What it establishes: the general Slavnov–Taylor identities for non-Abelian gauge theories, in a form usable to all orders.
What it does not: no non-perturbative statement.
Role in v3: with E09, the full set of Ward identities any rigorous construction must satisfy.

#### E11 · 't Hooft & Veltman (1972a) — dimensional regularization [T]
**G. 't Hooft & M. Veltman**, "Regularization and renormalization of gauge fields," *Nucl. Phys. B* **44** (1972) 189–213. DOI 10.1016/0550-3213(72)90279-9.
Status: **VERIFIED** (Garfield Citation Classic 1984; abstract confirms continuation to n ≠ 4 dimensions preserving unitarity, causality, gauge invariance, and disentangling overlapping divergences; independent discovery by Bollini–Giambiagi and Ashmore noted).
What it establishes: **dimensional regularization** — continue the dimension to n = 4 − ε; gauge invariance and Ward identities survive; overlapping divergences separate into 1/(n−4) poles. This supplies the gauge-invariant regulator that E07 assumed.
What it does not: dim reg is a UV regulator; it says nothing about the IR; the ε → 0 limit is formal.
Role in v3: the technical foundation of every β-function statement in Group C/D. Note: dim reg breaks γ₅ and parity in general — the parity restriction of E07 is not accidental.

#### E12 · 't Hooft & Veltman (1972b) — combinatorics [T]
**G. 't Hooft & M. Veltman**, "Combinatorics of gauge fields," *Nucl. Phys. B* **50** (1972) 318–353.
Status: **UNVERIFIED** (standard citation; not checked this session).
What it establishes: the general combinatorial algorithm (Cartan/Weyl) for computing group weights of arbitrary diagrams — the bookkeeping behind the 1972 renormalizability proof.
What it does not: no new physics.
Role in v3: citation only, for the group-theory machinery behind b₀.

#### E13 · Lee & Zinn-Justin (1972) — broken gauge symmetries, renormalizable gauges [T]
**B. W. Lee & J. Zinn-Justin**, "Spontaneously broken gauge symmetries. I, II, III," *Phys. Rev. D* **5** (1972) 3121, 3137, 3155.
Status: **UNVERIFIED** (standard citation; not checked this session).
What it establishes: an independent, more pedestrian derivation of renormalizability of spontaneously broken gauge theories in the R_ξ gauges, with explicit Ward identities; widely used textbook route.
What it does not: same scope as E08; no IR.
Role in v3: the formulation most textbooks follow; useful when v3 needs an explicit R_ξ propagator.

#### E14 · Bollini & Giambiagi (1972) — independent dim reg [T]
**C. G. Bollini & J. J. Giambiagi**, "Lowest order divergent graphs in ν-dimensional space," *Phys. Lett. B* **40** (1972) 566–568.
Status: **UNVERIFIED** (standard citation; independent-discovery priority confirmed indirectly via E11's Citation Classic note).
What it establishes: independent discovery of dimensional regularization.
What it does not: same as E12.
Role in v3: priority footnote only (also J. F. Ashmore, *Lett. Nuovo Cim.* **4** (1972) 289, unverified).

---

### C. Renormalization group and dimensional transmutation

#### E15 · Gell-Mann & Low (1954) — the RG equation [T] (pre-history)
**M. Gell-Mann & F. E. Low**, "Quantum electrodynamics at small distances," *Phys. Rev.* **95** (1954) 1300–1305.
Status: **UNVERIFIED** (foundational standard citation; outside the 1963–73 window, included as the origin of the RG; not checked this session).
What it establishes: the renormalization group and the invariant charge in QED; the β-function as the statement that the coupling depends on the scale. In QED the sign is positive (screening).
What it does not: no non-Abelian theory; no statement about the UV except the Landau-pole extrapolation.
Role in v3: the object that, run backwards, becomes asymptotic freedom. Every "running coupling" claim in v3 descends from this.

#### E16 · Callan (1970) — broken scale invariance [T]
**C. G. Callan Jr.**, "Broken scale invariance in scalar field theory," *Phys. Rev. D* **2** (1970) 1541–1547. DOI 10.1103/PhysRevD.2.1541.
Status: **VERIFIED** (ADS 1970PhRvD...2.1541C; APS; abstract confirms scaling Ward identities with anomalies, renormalization-group origin, and a route around overlapping divergences).
What it establishes: the Callan–Symanzik equation, derived from broken scale invariance; the scale anomaly.
What it does not: scalar theory only; no gauge theory.
Role in v3: one half of the Callan–Symanzik equation. The anomaly coefficient is the β-function — the same object as E07–E11's counterterms.

#### E17 · Symanzik (1970) — small-distance behaviour [T]
**K. Symanzik**, "Small distance behaviour in field theory and power counting," *Commun. Math. Phys.* **18** (1970) 227–246. DOI 10.1007/BF01649434.
Status: **VERIFIED** (ADS 1970CMaPh..18..227S; Springer/NII-REO; OSTI 4117755; abstract confirms the linear PDEs for asymptotic vertex functions, one universal function of one variable, relation to scale-invariance breaking).
What it establishes: the Callan–Symanzik equation derived independently from small-distance power counting; the universal function is obtained from non-asymptotic considerations.
What it does not: QED/scalar; no non-Abelian.
Role in v3: the second half of the CS equation; with E16 the standard citation is "Callan–Symanzik."

#### E18 · Symanzik (1971) — Wilson expansions [T]
**K. Symanzik**, "Small-distance-behaviour analysis and Wilson expansions," *Commun. Math. Phys.* **23** (1971) 49–86.
Status: **UNVERIFIED** (standard citation; not checked this session).
What it establishes: the operator-product expansion consistent with the CS equation; anomalous dimensions of composite operators.
What it does not: no gauge theory.
Role in v3: the anomalous dimension γ_O of ¼tr F² — directly used by v3's multiplicative residue transport f₁² = Z_O² f₁².

#### E19 · Wilson (1971) — RG and strong interactions [T]/[C]
**K. G. Wilson**, "Renormalization group and strong interactions," *Phys. Rev. D* **3** (1971) 1818–1846. DOI 10.1103/PhysRevD.3.1818.
Status: **VERIFIED** (APS; OSTI 4047125; INSPIRE ~757 citations; abstract confirms fixed point vs limit cycle, e⁺e⁻ test, RG determination of strong couplings).
What it establishes: the RG applied to strong interactions; the fixed-point/limit-cycle dichotomy; the proposal that a theory could be asymptotically free (though Wilson did not construct one here).
What it does not: no non-Abelian computation; the limit-cycle alternative is not the YM answer.
Role in v3: the conceptual bridge from the RG to the strong-interaction UV. Wilson's later lattice formulation (era 04) is the non-perturbative completion.

#### E20 · 't Hooft (1973) — dimensional regularization and the RG [T]
**G. 't Hooft**, "Dimensional regularization and the renormalization group," *Nucl. Phys. B* **61** (1973) 455–468. DOI 10.1016/0550-3213(73)90376-3. (CERN-TH-1666.)
Status: **VERIFIED** (CERN CDS 880603 with fulltext; ADS 1973NuPhB..61..455T; abstract confirms scale transformation via continuous dimension, single poles at n = 4, higher poles from single-pole coefficients, summation to exact results in some cases).
What it establishes: the minimal-subtraction scheme and the β-function as the residue of the 1/(n−4) pole; the RG equation in dim reg; the ε-expansion.
What it does not: no IR; the scheme dependence of higher coefficients is not the issue here (b₀, b₁ are universal).
Role in v3: **the** technical reference for "β(g) is the coefficient of the pole." v3's statement that the one- and two-loop coefficients are scheme-independent rests on this formulation.

#### E21 · Coleman & Weinberg (1973) — dimensional transmutation [T]/[C]
**S. Coleman & E. Weinberg**, "Radiative corrections as the origin of spontaneous symmetry breaking," *Phys. Rev. D* **7** (1973) 1888–1910. DOI 10.1103/PhysRevD.7.1888.
Status: **VERIFIED** (APS; OSTI 4468826; ~2,981 citations; abstract/context confirm one-loop effective potential, SSB without a negative mass-squared term).
What it establishes: the one-loop effective potential of massless scalar electrodynamics generates spontaneous symmetry breaking; a mass scale Λ ~ μ exp(−c/g²) appears from a dimensionless coupling — **dimensional transmutation**. This is the prototype of how Λ_QCD arises.
What it does not: it is *not* a mass gap. It generates a scale, and in that model a vacuum expectation value; it does not produce a gap in a pure non-Abelian spectrum, and the result is one-loop/conditional on the model.
Role in v3: **central distinction**. "Scale generation" ≠ "mass gap." v3 must state that CW shows the *mechanism by which a scale can appear*, not that YM has a gap. This is the cleanest place in the era to draw the theorem/conditional line.

---

### D. Asymptotic freedom

#### E22 · 't Hooft (1972, unpublished) → 1998 retrospective [T]/[O]
**G. 't Hooft**, result announced at the Colloquium on Renormalization of Yang–Mills Fields and Applications to Particle Physics, Marseille, June 1972 (C. P. Korthals Altes, organizer); unpublished. Retrospective: **G. 't Hooft**, "When was asymptotic freedom discovered? or The rehabilitation of quantum field theory," arXiv:hep-th/9808154; *Nucl. Phys. B (Proc. Suppl.)* **74** (1999) 413–425. DOI 10.1016/S0920-5632(99)00207-8.
Status: **VERIFIED** (arXiv abstract; full text read via ar5iv; quotes below are from the paper).
What it establishes (in 't Hooft's own words, paraphrased closely):
- "already in 1971 … I was immune to the Landau objections"; "By 1972, I had calculated the scaling behavior," with a negative gauge term; for SU(2) the β-function turns negative with under 11 fermions, for SU(3) the critical count is 33/2.
- At Marseille: "I announced at that meeting my finding that the coefficient determining the running of the coupling strength … is negative." Symanzik: "If this is true, it will be very important, and you should publish this result quickly." 't Hooft: "I did not follow his advice."
- On Vanyashin–Terent'ev: they found "a negative sign in the charge renormalization of charged vector mesons" and "attributed this 'absurd' result to the non-renormalizability of the theory"; ref. 34 dates it ZETP **48** (1965).
- On Khriplovich 1969: "correctly calculated the charge renormalization … in the Coulomb gauge … found the unusual sign, but a connection with asymptotic freedom was not made."
What it does not: the 1972 result was never published in the era; the *published* priority is E23–E24. Secondary sources disagree on dating (1971 computation vs 1972 announcement) — **use the 1998 retrospective as the primary source and do not assert a single date.**
Role in v3: the priority caveat. Cite as "'t Hooft (unpublished 1972; retrospective 1998)."

#### E23 · Gross & Wilczek (1973) — asymptotic freedom [T]
**D. J. Gross & F. Wilczek**, "Ultraviolet behavior of non-Abelian gauge theories," *Phys. Rev. Lett.* **30** (1973) 1343–1346. DOI 10.1103/PhysRevLett.30.1343.
Status: **VERIFIED** (APS Milestone; abstract confirms a wide class of non-Abelian theories have free-field asymptotic behavior up to calculable log corrections; suggests Bjorken scaling from non-Abelian gauge symmetry; companion PRD **8**, 3633).
What it establishes: the one-loop β-function of non-Abelian gauge theory is negative for a wide class of theories; the UV is asymptotically free. For pure SU(N): β(g) = −b₀ g³/(16π²) + O(g⁵), b₀ = 11N/3 > 0.
What it does not: perturbative only; no IR, no confinement, no gap. "Asymptotically free" is a UV statement.
Role in v3: one of the two published AF papers. The sign of b₀ is scheme-independent — this is a genuine [T], not an artifact.

#### E24 · Politzer (1973) — asymptotic freedom [T]
**H. D. Politzer**, "Reliable perturbative results for strong interactions?" *Phys. Rev. Lett.* **30** (1973) 1346–1349. DOI 10.1103/PhysRevLett.30.1346.
Status: **VERIFIED** (APS; Caltech Authors; abstract confirms perturbation theory arbitrarily good for deep Euclidean Green's functions of any Yang–Mills theory and many with fermions; under dynamical-SSB hypothesis these are asymptotic forms of the physical solution).
What it establishes: independent one-loop computation; the same negative sign; the explicit statement that the symmetric (deep Euclidean) Green's functions are the asymptotic forms.
What it does not: same as E23.
Role in v3: the second published AF paper; the "deep Euclidean" phrasing is the one v3 should use when it says AF constrains the UV, not the physical spectrum.

#### E25 · Gross & Wilczek (1973b) — AF gauge theories I [T]
**D. J. Gross & F. Wilczek**, "Asymptotically free gauge theories. I," *Phys. Rev. D* **8** (1973) 3633–3647.
Status: **VERIFIED indirectly** (confirmed as companion in the E23 APS record: "Companion paper … Phys. Rev. D 8 (10), 3633 (1973), with follow-up … D 9 (4), 980 (1974)").
What it establishes: the detailed derivation, the two-loop coefficient, the matter-content dependence.
What it does not: no IR.
Role in v3: the technical reference for the two-loop coefficient b₁ (also universal).

#### E26 · Vanyashin & Terent'ev (1965) — the missed sign [T]
**V. S. Vanyashin & M. V. Terent'ev**, "The vacuum polarization of a charged vector field," *Zh. Eksp. Teor. Fiz.* **48** (1965) [*Sov. Phys. JETP* **21** (1965)].
Status: **UNVERIFIED** (content confirmed via E22, which quotes them finding "a negative sign in the charge renormalization of charged vector mesons" and attributing it to non-renormalizability; exact pages not checked this session — do not quote pages).
What it establishes: the negative sign in charged-vector vacuum polarization was seen in 1965 and dismissed.
What it does not: no correct interpretation.
Role in v3: historical [O]-adjacent note; cite through 't Hooft 1998.

#### E27 · Caswell (1974) — two-loop [T] (continuation)
**W. E. Caswell**, "Asymptotic behavior of non-Abelian gauge theories to two-loop order," *Phys. Rev. Lett.* **33** (1974) 244–246.
Status: **UNVERIFIED** (standard citation; just outside the window; not checked this session).
What it establishes: the two-loop β-function confirms AF.
What it does not: no IR.
Role in v3: the two-loop universality anchor; belongs to era 03 but is the immediate continuation.

---

### E. The infrared problem begins

#### E28 · Kinoshita (1962) — mass singularities [T] (pre-history)
**T. Kinoshita**, "Mass singularities of Feynman amplitudes," *J. Math. Phys.* **3** (1962) 650–677.
Status: **UNVERIFIED** (standard citation; not checked this session).
What it establishes: the systematic treatment of mass (collinear) singularities in Feynman amplitudes.
What it does not: no non-Abelian; no confinement.
Role in v3: the first precise statement that the IR is a separate divergence problem from the UV.

#### E29 · Lee & Nauenberg (1964) — KLN theorem [T]
**T. D. Lee & M. Nauenberg**, "Degenerate systems and mass singularities," *Phys. Rev.* **133** (1964) B1549–B1562.
Status: **UNVERIFIED** (standard citation; not checked this session).
What it establishes: the KLN theorem — sufficiently inclusive cross sections are free of mass singularities. The IR divergences cancel in physical quantities.
What it does not: KLN is about *inclusive* quantities in perturbation theory; it does not construct the IR limit, and it does not apply to the confined spectrum.
Role in v3: important negative control — v3 must not confuse "IR divergences cancel in inclusive observables" with "there is no IR problem." The mass gap is not an inclusive cross section.

#### E30 · Schwinger (1962) — 2D gauge theory with mass [T] (pre-history)
**J. Schwinger**, "Gauge invariance and mass. II," *Phys. Rev.* **128** (1962) 2425–2432.
Status: **UNVERIFIED** (standard citation; not checked this session).
What it establishes: the 1+1-dimensional Schwinger model — an exactly solvable gauge theory whose spectrum has a mass (the "Schwinger boson") despite massless constituents. The only early exactly solvable gauge theory with a gap.
What it does not: 1+1 dimensions; Abelian; the mechanism (axial anomaly) does not transfer to 4D pure YM.
Role in v3: the cautionary exact example. A mass gap in a gauge theory is possible; the 2D mechanism is not the 4D one. Cite as precedent, not as evidence.

#### E31 · Gribov (1978) — Gribov copies [NG] (continuation)
**V. N. Gribov**, "Quantization of non-Abelian gauge theories," *Nucl. Phys. B* **139** (1978) 1–19.
Status: **UNVERIFIED** (standard citation; just outside the window; not checked this session).
What it establishes: **no-go** — the Faddeev–Popov gauge-fixing condition is not globally well-defined; there are gauge-equivalent field configurations satisfying the same condition (Gribov copies). The FP determinant changes sign; the perturbative measure is incomplete. Restricting to the Gribov region modifies the IR gluon propagator (the "Gribov horizon").
What it does not: no rigorous confinement proof; the Gribov-Zwanziger action is a model.
Role in v3: **the load-bearing no-go for any non-perturbative use of Faddeev–Popov.** Any v3 construction that starts from the FP measure must either work gauge-invariantly or explicitly address the Gribov ambiguity.

#### E32 · Becchi, Rouet & Stora (1976) — BRST [T] (continuation)
**C. Becchi, A. Rouet & R. Stora**, "Renormalization of gauge theories," *Ann. Phys.* **98** (1976) 287–321. DOI 10.1016/0003-4916(76)90156-1.
Status: **VERIFIED** (CERN CDS 406093; MathSciNet MR0413861; Scopus ~1,443 citations; companion "Renormalization of the Abelian Higgs model," *Commun. Math. Phys.* **42** (1975) 127).
What it establishes: the global nilpotent symmetry s (BRST) of the gauge-fixed action, s² = 0; physical states are the cohomology H⁰(s). This is the operator-algebraic form of the Slavnov–Taylor identities (E09–E10) and the correct statement of what "gauge invariance" means after gauge-fixing.
What it does not: BRST is perturbative; it does not by itself give the physical Hilbert space's positivity (that is Kugo–Ojima, era 03) and says nothing about the gap.
Role in v3: the cleanest rigorous handle on "physical subspace." If v3 wants to define the physical spectrum, BRST cohomology is the tool; the gap question becomes "is there a mass gap in H⁰(s)?"

#### E33 · Tyutin (1975) — BRST, independent [T] (continuation)
**I. V. Tyutin**, "Gauge invariance in field theory and statistical physics in operator formalism," Lebedev Physical Institute preprint FIAN-39 (1975); arXiv:0812.0580 [hep-th].
Status: **VERIFIED** (INSPIRE record, ~560–580 citations; arXiv:0812.0580; abstract confirms Ward identities and gauge dependence from canonical commutation relations and Heisenberg equations, with and without SSB).
What it establishes: independent discovery of the BRST symmetry, in operator language, including the statistical (finite-temperature) case.
What it does not: same as E32.
Role in v3: priority with E32; the operator derivation is the one that connects to the canonical/Hamiltonian side of v3.

---

## 2. What renormalization proved, and what it did not

**Proved (perturbatively, order by order):**
1. [T] Massless Yang–Mills is renormalizable, conditional on a gauge-invariant regulator (E07), completed by dimensional regularization (E11).
2. [T] Spontaneously broken gauge theory is renormalizable and unitary (E08, E13).
3. [T] A single renormalized coupling exists — Z₁/Z₂ common to vector and ghost (E09, E10).
4. [T] The S-matrix is gauge-independent on shell (E06).
5. [T] The one-loop β-function of pure SU(N) is negative; b₀ = 11N/3 (E23, E24). b₀ and b₁ are scheme-independent in mass-independent schemes (E20).
6. [T] A dimensionless coupling can generate a mass scale by dimensional transmutation (E21), in the model studied.

**Scheme-dependent (not theorem-level in the same sense):**
- The numerical value of Λ, of α_s(μ), of the running coupling at any scale.
- β-function coefficients from three loops on (scheme-dependent); b₀, b₁ are universal.
- The Gribov-Zwanziger IR propagator (E31) — a gauge-fixing-dependent model.

**Not proved — and this is the honest boundary of the era:**
- [O] That 4D Yang–Mills exists as a mathematical quantum field theory satisfying the Osterwalder–Schrader axioms. The path integral of E03 is formal; no measure is constructed.
- [O] That the perturbation series converges. It is asymptotic at best (Dyson-type argument); AF makes the UV logs small, not the series convergent.
- [O] Any statement about the infrared: confinement, the spectrum, a mass gap. AF's IR corollary (the coupling grows) makes perturbation theory *worse* in the IR, not better.
- [O] That the negative β-function survives non-perturbatively (it does, in the sense that the sign of b₀ is universal, but the *running* at strong coupling is not under perturbative control).
- [NG] Faddeev–Popov quantization is globally incomplete (E31) — a no-go for naive non-perturbative use of the era's measure.

**The precise missing step of the era:** the era computed the UV and discovered that the IR is where the coupling becomes strong, then stopped. It did not — and could not — connect the perturbative UV to a non-perturbative IR. That connection is the campaign's target, and it is not in this literature.

---

## 3. Theorem vs scheme vs open — one table

| Claim | Status | Why |
|---|---|---|
| YM renormalizable, order by order | [T] | E07 + E11; conditional on gauge-invariant regulator, which E11 supplies |
| Single renormalized coupling | [T] | E09, E10 (Slavnov–Taylor) |
| S-matrix gauge-independent | [T] | E06 (on shell) |
| b₀ = 11N/3, β < 0 for pure SU(N) | [T] | E23, E24; b₀ scheme-independent |
| b₁ universal | [T] | E25, E27 |
| Scale from dimensionless coupling | [T] (model) | E21; one-loop, scalar QED |
| Λ_QCD value, α_s(μ) value | scheme-dependent | E20 |
| Gribov-Zwanziger propagator | scheme-dependent model | E31 |
| 4D YM measure exists (OS axioms) | open [O] | not in this era |
| Perturbation series converges | open (false, expected asymptotic) | Dyson argument |
| Confinement | open [O] | not in this era |
| Mass gap | open [O] | not in this era; the Clay problem |

---

## 4. Position in 14 eras

- **01 · 1954–65** — The Yang–Mills Lagrangian and the first closed-loop rule (E01) are the classical object every later era quantizes; the ghost is the first sign that the gauge orbit is non-trivial.
- **02 · 1963–73 (this era)** — FP + DeWitt quantization, 't Hooft–Veltman renormalizability, the RG, AF: the UV is solved order by order and the IR is exposed as the remaining problem.
- **03 · 1974–79** — Monopoles, instantons, Gribov copies (E31), Kugo–Ojima: the non-perturbative IR sector that era 02 left untouched becomes the main object; AF's IR growth is renamed "confinement."
- **04 · 1974–82** — Wilson's lattice gauge theory turns the FP/DeWitt path integral into a definition and gives the first strong-coupling mass-gap theorem — the non-perturbative completion era 02 could not supply.
- **05 · 1982–90** — Continuum limits and scaling: the two-loop universality of era 02 (E20, E25) is what makes Λ_QCD and scaling well-defined.
- **06 · 1985–95** — Rigorous constructive QFT and Balaban's RG: an attempt to make the era-02 measure rigorous; it controls the UV but not the IR gap.
- **07 · 1993–2000** — Seiberg–Witten and AdS/CFT use the era-02 RG and AF as inputs; the mass gap is reformulated in supersymmetric and holographic settings, not solved in pure YM.
- **08 · 1998–2005** — Confinement mechanisms (dual superconductivity, center vortices) are the IR completion of the era-02 running coupling; they assume the gap rather than prove it.
- **09 · 2002–10** — Large-N ('t Hooft 1974) and lower-dimensional rigorous gaps: era 02's 1/N bookkeeping becomes a limit in which the gap can sometimes be proven.
- **10 · 2008–14** — Lattice precision glueball spectra: the numerical mass gap that era 02 could not reach; the β-function of era 02 is the input that fixes the scale.
- **11 · 2013–18** — Resurgence and trans-series: the asymptotic (not convergent) series that era 02's AF produces is reorganized; the IR renormalon is the era-02 Landau-pole problem in modern dress.
- **12 · 2018–21** — Bootstrap and rigorous bounds use reflection positivity and the OS axioms — exactly the structure era 02's formal measure fails to provide.
- **13 · 2021–24** — Gradient flow and rigorous RG steps make the era-02 β-function a rigorous object; the flow is the modern version of E20.
- **14 · 2024–26** — Dimensional transmutation (E21) and the AF scale are treated as the only non-perturbative scale; the missing machine must transport the strong-coupling gap across the Savvidy crossover — the step era 02 identified but could not take.

---

## 5. Gaps and attack surfaces for v3

- **Renormalizability is UV-only.** v3 must state the era-02 theorem as "order-by-order UV finiteness, conditional on a regulator," and must not let it stand in for any IR statement. This is the single most common overclaim in the secondary literature.
- **b₀ and b₁ universality is the honest core of AF.** The sign of b₀ is scheme-independent; use E20/E23/E25 to anchor "AF is not a gauge artifact." Do not extend the claim to the running at strong coupling.
- **Separate "scale generation" from "mass gap."** E21 (Coleman–Weinberg) is the cleanest place to draw the line: a scale Λ ~ μ e^{−c/g²} appears, but no gap is produced. v3's dimensional-transmutation argument must not silently upgrade one to the other.
- **Faddeev–Popov is globally incomplete (E31).** Any v3 construction built on the FP measure must either be gauge-invariant or explicitly handle Gribov copies. This is a genuine no-go, not a technicality.
- **BRST cohomology is the right definition of "physical."** Use E32/E33 to define the physical subspace; then the gap question is a spectral question on H⁰(s). This is a rigorous handle the era itself did not exploit.
- **Slavnov–Taylor identities constrain the IR propagator.** E09/E10 are exact; the longitudinal part of the gluon two-point function is fixed. Any claimed IR mass must be consistent with them.
- **KLN is not an IR solution.** E29 says inclusive cross sections are finite; it does not say the IR limit exists or that the spectrum has a gap. v3 must not cite KLN as evidence against an IR problem.
- **Priority discipline on 't Hooft 1972.** Use E22 (the 1998 retrospective) as primary; secondary sources disagree on 1971 vs 1972. Do not assert a publication date that does not exist.
- **The 2D Schwinger model is a precedent, not evidence.** E30 shows a gauge theory can have a gap; the 2D mechanism (axial anomaly) does not transfer. Cite as cautionary.
- **The exact missing step is a bridge, not a computation.** Era 02 ends with a perturbative UV and an uncontrolled IR. v3's job is to build the object that connects them; no citation in this era does that, and v3 should say so explicitly.

---

## 6. Verification ledger

| Ref | Status | Source checked |
|---|---|---|
| E01 Feynman 1963 | UNVERIFIED | standard citation; not checked (budget) |
| E02 DeWitt 1967 | VERIFIED | APS / Google Scholar record; abstract confirms YM+gravity, fictitious quanta |
| E03 Faddeev–Popov 1967 | VERIFIED | CERN CDS 431589; INSPIRE; abstract confirms ghost method |
| E04 Veltman 1968 | VERIFIED | ADS 1968NuPhB...7..637V; INSPIRE 52788 |
| E05 Khriplovich 1969 | UNVERIFIED | content via 't Hooft 1998 (E22); pages not checked |
| E06 Fradkin–Tyutin 1970 | VERIFIED | ADS 1970PhRvD...2.2841F; INSPIRE ~255 cites |
| E07 't Hooft 1971a | VERIFIED | ADS 1971NuPhB..33..173T; abstract checked in full |
| E08 't Hooft 1971b | VERIFIED | Scopus 0000043313; OSTI 4711793; abstract checked |
| E09 Taylor 1971 | VERIFIED | MathSciNet MR0314386; abstract checked |
| E10 Slavnov 1972 | VERIFIED | ADS 1972TMP....10...99S; Springer; pages 99–104 |
| E11 't Hooft–Veltman 1972a | VERIFIED | Garfield Citation Classic 1984; abstract checked |
| E12 't Hooft–Veltman 1972b | UNVERIFIED | standard citation; not checked (budget) |
| E13 Lee–Zinn-Justin 1972 | UNVERIFIED | standard citation; not checked (budget) |
| E14 Bollini–Giambiagi 1972 | UNVERIFIED | priority noted in E11 record; pages not checked |
| E15 Gell-Mann–Low 1954 | UNVERIFIED | standard citation; outside window; not checked |
| E16 Callan 1970 | VERIFIED | ADS 1970PhRvD...2.1541C; APS; abstract checked |
| E17 Symanzik 1970 | VERIFIED | ADS 1970CMaPh..18..227S; OSTI 4117755; abstract checked |
| E18 Symanzik 1971 | UNVERIFIED | standard citation; not checked (budget) |
| E19 Wilson 1971 | VERIFIED | APS; OSTI 4047125; INSPIRE ~757 cites; abstract checked |
| E20 't Hooft 1973 | VERIFIED | CERN CDS 880603 (fulltext); ADS 1973NuPhB..61..455T; abstract checked |
| E21 Coleman–Weinberg 1973 | VERIFIED | APS; OSTI 4468826; abstract/context checked |
| E22 't Hooft 1972/1998 | VERIFIED | arXiv hep-th/9808154 abstract; full text read via ar5iv; quotes recorded |
| E23 Gross–Wilczek 1973 | VERIFIED | APS Milestone; abstract checked |
| E24 Politzer 1973 | VERIFIED | APS; Caltech Authors; abstract checked |
| E25 Gross–Wilczek 1973b | VERIFIED (indirect) | companion confirmed in E23 APS record |
| E26 Vanyashin–Terent'ev 1965 | UNVERIFIED | content and ZETP 48 via E22; pages not checked |
| E27 Caswell 1974 | UNVERIFIED | standard citation; outside window; not checked |
| E28 Kinoshita 1962 | UNVERIFIED | standard citation; not checked (budget) |
| E29 Lee–Nauenberg 1964 | UNVERIFIED | standard citation; not checked (budget) |
| E30 Schwinger 1962 | UNVERIFIED | standard citation; not checked (budget) |
| E31 Gribov 1978 | UNVERIFIED | standard citation; outside window; not checked |
| E32 BRST 1976 | VERIFIED | CERN CDS 406093; MathSciNet MR0413861; abstract checked |
| E33 Tyutin 1975 | VERIFIED | INSPIRE; arXiv:0812.0580; abstract checked |

**Totals: 33 entries — 20 VERIFIED, 13 UNVERIFIED.** The 13 unverified are either (a) pre-history/context entries outside the load-bearing chain, or (b) standard citations whose content is corroborated by a verified source but whose exact pages were not independently checked this session. **None was invented.** A follow-up pass with a fresh search budget should close E05, E12–E15, E18, E26–E31; until then, v3 must not quote their page ranges without a check.

---

*File: `D:\Dev\math\research\lit-v3\LIT-ERA-02.md`. Era 02 of 14. No claim in this file asserts that the Yang–Mills mass gap is proved.*
