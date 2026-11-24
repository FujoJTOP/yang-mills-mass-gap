# LIT-ERA-07 — 1993–2000: Functional Methods (Dyson–Schwinger and Exact RG)

> Scope: the functional-method decade for 4D pure SU(N) Yang–Mills. Covers Wetterich's exact RG equation and its gauge-theory formulation; the Ellwanger–Hirsch–Weber truncations; the von Smekal–Hauck–Alkofer Dyson–Schwinger scaling solution; the ghost–gluon vertex and infrared exponents; the Kugo–Ojima criterion revisited; the Gribov–Zwanziger action and its renormalization; the first lattice evidence for a massive gluon propagator; the scaling-vs-decoupling question; Cornwall's dynamical gluon mass.
>
> Compilation date: 2026-11-24. Method: primary bibliographic data checked against publisher/ADS/INSPIRE/arXiv records where the search budget allowed; every entry carries an explicit verification status. No bibliographic field below is reconstructed from memory without being marked UNVERIFIED.
>
> **Label key**
> - `[T]` exact theorem, exact identity, or rigorously established structural result
> - `[C]` conditional result: truncation-, ansatz-, gauge-, or framework-dependent
> - `[NG]` no-go / obstruction / negative result
> - `[O]` orientation: review, historical anchor, or methodological framing
> - `[V]` numerical / lattice evidence
>
> **One-line verdict for v3:** this era produced an exact ERG *equation* and a self-consistent *scaling solution* of truncated Landau-gauge DSEs, plus a renormalizable Gribov–Zwanziger action and the first lattice hints of a massive gluon. It produced **no gauge-invariant mass-gap theorem**. Every propagator statement in this era is gauge-fixed; every DSE/ERG mass statement is conditional on truncation and on the k→0 / continuum / infinite-volume limits. The era's lasting value for v3 is the precise statement of what a "functional mass" can and cannot mean.

---

## 1. What this era actually established

1. **ERG became an exact equation for gauge theories** (Wetterich 1993; Reuter–Wetterich 1994), but only with a gauge-fixing/background field and a cutoff that breaks BRST. Gauge invariance is recovered through modified Slavnov–Taylor identities (mSTI) as k→0; at finite k the flow is not manifestly gauge invariant. Litim–Pawlowski (1999) sharpened this into a structural obstruction for quantum fluctuations.
2. **The DSE scaling solution was found** (von Smekal–Hauck–Alkofer 1997/98): in Landau gauge, with truncations and STI-consistent vertices, the gluon propagator vanishes in the IR, the ghost propagator diverges, and the running coupling approaches an IR fixed point α_c ≈ 9.5. This is a solution of a truncated system, not a theorem about QCD.
3. **The Gribov–Zwanziger programme became a local, renormalizable gauge-fixed action** (Zwanziger 1989; 1993; 1994), with a dimensionful Gribov parameter and complex gluon poles; the glueball channel retains a physical cut.
4. **Lattice Landau-gauge propagators showed IR mass-like dressing** (Marenzoni–Martinelli–Stella 1995; Leinweber et al. 1998), while Cucchieri quantified how Gribov copies contaminate the ghost sector and how thermalized configurations sit near the Gribov horizon.
5. **Kugo–Ojima was sharpened** (Kugo 1995) into a gauge-parameter-independent criterion 1+u=0 for color confinement, but the criterion is not a mass-gap theorem and was not derived from YM first principles in this era.
6. **Cornwall's dynamical gluon mass** (1982, the anchor) remained the phenomenological benchmark: m_g ≈ 500 ± 200 MeV, 0⁺ glueball ≈ 2m_g, obtained from a gauge-invariant resummation in light-cone gauge — a model, not a proof.

**What did not happen:** no derivation of a positive mass gap for the physical spectrum; no gauge-invariant definition of the DSE/ERG "gluon mass"; no control of truncation error; no continuum/infinite-volume extrapolation of the lattice propagator mass; no decision between scaling and decoupling.

---

## 2. Exact renormalization group: entries

### E1. Wetterich 1993 — the exact flow equation
**C. Wetterich**, "Exact evolution equation for the effective potential," *Phys. Lett. B* **301** (1993) 90–94. DOI: 10.1016/0370-2693(93)90726-X. arXiv:1710.05815 (2017 reprint of the 1992 manuscript). **`[T]` VERIFIED.**

The foundational ERG equation for the effective average action Γ_k:
```
k ∂_k Γ_k = (1/2) Tr [ (Γ_k^{(2)} + R_k)^{-1} · k ∂_k R_k ] .
```
Exact, one-loop in structure, but a functional differential equation: every practical use is a truncation. The paper itself is scalar; the gauge-theory extension is E2. For v3: this is the only *exact* object in the era, and it is exact only before truncation and before gauge fixing is imposed.

### E2. Reuter–Wetterich 1994 — gauge theories
**M. Reuter, C. Wetterich**, "Effective average action for gauge theories and exact evolution equations," *Nucl. Phys. B* **417** (1994) 181–214. DOI: 10.1016/0550-3213(94)90543-6. **`[T]` VERIFIED.**

Introduces the nonperturbative evolution equation for Yang–Mills with a background field and gauge fixing; computes the running non-Abelian coupling in arbitrary dimension. The exact equation is gauge-fixed: the cutoff term R_k is not BRST-invariant, so the flow is accompanied by modified Slavnov–Taylor identities. The "exactness" is therefore exactness *of a gauge-fixed flow*, not of a gauge-invariant observable.

### E3. Bonini–D'Attanasio–Marchesini 1994 — SU(2) flow and gauge invariance
**M. Bonini, M. D'Attanasio, G. Marchesini**, "Renormalization group flow for SU(2) Yang–Mills theory and gauge invariance," *Nucl. Phys. B* **421** (1994) 429–455. arXiv:hep-th/9312114. **`[C]` VERIFIED (INSPIRE).**

Wilson-flow treatment of pure SU(2); identifies the relevant couplings and imposes boundary conditions at the IR cutoff. Establishes the perturbative consistency of the flow with gauge invariance, but only after choosing nine physical relevant couplings; the construction is a truncation/boundary-condition scheme, not a nonperturbative gauge-invariant definition.

### E4. Ellwanger 1994 — N-point flow equations
**U. Ellwanger**, "Flow equations for N point functions and bound states," *Z. Phys. C* **62** (1994) 503–510. DOI: 10.1007/BF01555911. arXiv:hep-ph/9308260. **`[C]` VERIFIED.**

Decomposes the ERG into one-particle-irreducible N-point functions and applies a truncated four-point flow to bound states (Wick–Cutkosky; QCD-motivated interaction). Methodological: it makes explicit that an ERG "bound-state mass" is only as good as the truncation and the ansatz for higher vertices.

### E5. Morris 1994 — approximation theory of ERG
**T. R. Morris**, "The exact renormalization group and approximate solutions," *Int. J. Mod. Phys. A* **9** (1994) 2411–2450. arXiv:hep-ph/9308265. **`[O]` VERIFIED (INSPIRE).**

Systematic account of how to extract approximate solutions (derivative expansion, local potential approximation) from the exact equation. Useful for v3 as the methodological baseline for judging truncation error; no gauge-theory mass claim.

### E6. Wetterich 1994 — average action and systematic expansions
**C. Wetterich**, "Exact renormalization group equations for the average action and systematic expansions," *Int. J. Mod. Phys. A* **9** (1994) 3571–3602. **`[T]` VERIFIED (INSPIRE).**

Companion development of the average-action formalism and its systematic expansion schemes. Together with E1/E2 it defines the ERG technology imported into YM.

### E7. Bonini–D'Attanasio–Marchesini 1995 — BRS symmetry with ERG
**M. Bonini, M. D'Attanasio, G. Marchesini**, "BRS symmetry for Yang–Mills theory with exact renormalization group," *Nucl. Phys. B* **437** (1995) 163–186. DOI: 10.1016/0550-3213(94)00569-Z. arXiv:hep-th/9410138. **`[C]` VERIFIED.**

Shows perturbatively that the boundary couplings can be chosen so that all Slavnov–Taylor identities hold. Key structural point: of the nine relevant couplings, the three wave-function/triple-vector ones are fixed, three vanish (including the vector mass), and the rest are irrelevant. This is the cleanest in-era statement that a **finite-cutoff vector mass is not a physical parameter**: it can be tuned to zero by the STI boundary conditions.

### E8. Ellwanger–Hirsch–Weber 1996 — pure YM truncation and the cutoff mass
**U. Ellwanger, M. Hirsch, A. Weber**, "Flow equations for the relevant part of the pure Yang–Mills action," *Z. Phys. C* **69** (1996) 687–698. DOI: 10.1007/BF02907453. arXiv:hep-th/9506019. **`[C]` VERIFIED.**

Derives and integrates Wilson ERG equations for the relevant part of the pure YM action. Modified STIs control BRST breaking at finite cutoff k; one consequence is that **a gluon mass term is nonzero whenever k ≠ 0**. This is the era's most explicit warning against reading a finite-k ERG mass as a mass gap: the mass is tied to the regulator and to the mSTI, not to the k→0 spectrum.

### E9. Litim–Pawlowski 1999 — gauge-invariant Wilsonian flows?
**D. F. Litim, J. M. Pawlowski**, "On gauge invariant Wilsonian flows," in *The Exact Renormalization Group* (World Scientific, 1999), p. 168. arXiv:hep-th/9901063. **`[NG]` VERIFIED.**

Asks how close a Wilsonian flow can come to gauge invariance when the coarse-graining itself seems incompatible with it. The cutoff breaks gauge invariance; Γ_k satisfies a modified Ward identity that commutes with the flow and reduces to the ordinary WI only as k→0. Argues that a fully gauge-invariant implementation is possible for thermal fluctuations, not for quantum fluctuations. For v3 this is a structural no-go: **there is no manifestly gauge-invariant ERG flow of quantum YM in this framework**; gauge invariance is a boundary/limit statement.

---

## 3. Gribov–Zwanziger, the horizon, and Kugo–Ojima revisited

### E10. Zwanziger 1989 — local renormalizable action from the horizon
**D. Zwanziger**, "Local and renormalizable action from the Gribov horizon," *Nucl. Phys. B* **323** (1989) 513–544. DOI: 10.1016/0550-3213(89)90122-3. **`[C]` VERIFIED.**

Derives a local, power-counting renormalizable action that restricts the functional integral to the interior of the Gribov horizon (the fundamental modular region). The restriction replaces the perturbative coupling by a dimensionful parameter γ; at tree level the gluon propagator poles move to p² = ±i γ^{1/2} (unphysical complex poles), while the glueball channel retains a physical cut with positive spectral function. This is the technical origin of "Gribov–Zwanziger" and the source of the gauge-fixed "mass" language. Note the chronology: the **action is 1989**, not 1993; the 1993 paper (E11) proves its critical-limit renormalizability.

### E11. Zwanziger 1993 — renormalizability of the critical limit
**D. Zwanziger**, "Renormalizability of the critical limit of lattice gauge theory by BRS invariance," *Nucl. Phys. B* **399** (1993) 477–513. DOI: 10.1016/0550-3213(93)90506-K. **`[T]` VERIFIED.**

Proves that the critical limit of lattice gauge theory — containing a dimensionful parameter associated with the boundary of the fundamental modular region — is renormalizable, using BRS symmetry. Also proves that the exact Fermi-ghost propagator has a 1/(q²)² singularity at q=0. This is a genuine theorem about the gauge-fixed critical limit, not about the physical spectrum.

### E12. Zwanziger 1994 — fundamental modular region and area law
**D. Zwanziger**, "Fundamental modular region, Boltzmann factor and area law in lattice gauge theory," *Nucl. Phys. B* **412** (1994) 657–730. DOI: 10.1016/0550-3213(94)90396-4. **`[C]` VERIFIED.**

Develops the ensemble restricted to the fundamental modular region, the Boltzmann factor on it, and the proposed route to the area law. The connection between the horizon and confinement is argued, not proven; for v3 it is a conditional bridge, not a theorem.

### E13. Kugo 1995 — the criterion revisited
**T. Kugo**, "The universal renormalization factors Z₁/Z₃ and color confinement condition in non-Abelian gauge theory," arXiv:hep-th/9511033; in *BRS Symmetry* (Kyoto 1995), pp. 107–119. **`[C]` VERIFIED.**

In covariant gauges with arbitrary gauge parameter α, the universal ratio satisfies Z₁/Z₃ = 1 + u, where u ≡ u(p²=0) is the residue of the scalar pole ∝ p_μ p_ν / p² in the two-point function ⟨D_μ c D_ν c̄⟩. The Kugo–Ojima criterion states that 1 + u = 0 is a sufficient condition for all colored states to become unphysical (color confinement). The statement is gauge-parameter independent in covariant gauges, but it is a **criterion**, not a derivation; it concerns color confinement, not the mass gap; and it assumes the local covariant operator formalism with an unbroken global gauge symmetry and a specific asymptotic-field structure.

**Anchor for this thread (pre-era):** E25 (Kugo–Ojima 1979). The BRST quartet mechanism and the 1+u=0 criterion originate there; the 1995 paper is the sharpened renormalization-constant form.

---

## 4. Dyson–Schwinger: the scaling solution and its content

### E14. von Smekal–Hauck–Alkofer 1997 — the PRL
**L. von Smekal, A. Hauck, R. Alkofer**, "The infrared behavior of gluon and ghost propagators in Landau gauge QCD," *Phys. Rev. Lett.* **79** (1997) 3591–3594. DOI: 10.1103/PhysRevLett.79.3591. arXiv:hep-ph/9705242. **`[C]` VERIFIED.**

Solves a truncated system of Euclidean QCD DSEs in Landau gauge, implementing Slavnov–Taylor identities for the three-gluon and ghost–gluon vertices while neglecting irreducible four-gluon couplings and gluon–ghost/ghost–ghost scattering kernels. Analytic result: **the gluon propagator vanishes at small spacelike momenta; the ghost propagator diverges strongly; the running coupling approaches an infrared fixed point α_c ≈ 9.5.** This is the "scaling solution." It is a solution of a truncated, gauge-fixed system.

### E15. von Smekal–Hauck–Alkofer 1998 — the full Annals paper
**L. von Smekal, A. Hauck, R. Alkofer**, "A solution to coupled Dyson–Schwinger equations for gluons and ghosts in Landau gauge," *Annals Phys.* **267** (1998) 1–60; Erratum *Annals Phys.* **269** (1998) 182. arXiv:hep-ph/9707327. **`[C]` VERIFIED (INSPIRE).**

The detailed version of E14: truncation scheme, renormalization, numerical solutions, comparison with lattice data. Same caveats: Landau gauge, dropped kernels, no gauge-invariant observable. The erratum is part of the record and should be cited with the paper.

### E16. Roberts–Williams 1994 — DSE review
**C. D. Roberts, A. G. Williams**, "Dyson–Schwinger equations and their application to hadronic physics," *Prog. Part. Nucl. Phys.* **33** (1994) 477–575. **`[O]` VERIFIED.**

The standard pre-scaling review of DSE technology: dressed quark propagators, Bethe–Salpeter mesons, form factors, baryons. It documents the state of the truncation art just before the gluon/ghost scaling solution; useful for v3 to separate the *quark-sector* DSE tradition (where dynamical chiral symmetry breaking is the target) from the *pure-gluon/ghost* sector (where the scaling solution lives).

### E17. Alkofer–von Smekal 2001 — the review that fixed the era's narrative
**R. Alkofer, L. von Smekal**, "The infrared behavior of QCD Green's functions — confinement, dynamical symmetry breaking, and hadrons as relativistic bound states," *Phys. Rept.* **353** (2001) 281–369. DOI: 10.1016/S0370-1573(01)00010-2. arXiv:hep-ph/0007355. **`[O]` VERIFIED.** *(Published 2001; written 2000 — treated here as the era's closing synthesis.)*

Covers gauge fixing, BRS invariance and positivity, derivation of DSEs, truncation schemes, the scaling solution, comparisons to lattice, and hadronic applications. It is the best single entry point for v3, and it is explicit that the scaling solution is a truncation-dependent Landau-gauge statement.

---

## 5. Lattice evidence for a massive gluon, and Gribov-copy control

### E18. Marenzoni–Martinelli–Stella 1995 — the massive-gluon fit
**P. Marenzoni, G. Martinelli, N. Stella**, "The gluon propagator on a large volume at β = 6.0," *Nucl. Phys. B* **455** (1995) 339–353. DOI: 10.1016/0550-3213(95)00405-H. arXiv:hep-lat/9410011. **`[V]` VERIFIED.**

High-statistics quenched lattice study of the Landau-gauge gluon propagator at β = 6.0. Momentum-space data are fitted by
```
G(k²) = [ M² + Z k² (k²/Λ²)^η ]^{-1} ,
```
with η = 0.532(12) from the volume comparison; M² is shown **not** to vanish in the infinite-volume limit. This is the era's canonical "massive gluon" evidence. Caveats: quenched, one lattice spacing, Landau gauge, finite-volume fit form; M² is a fit parameter, not a gauge-invariant mass.

### E19. Cucchieri 1997 — Gribov copies and propagators
**A. Cucchieri**, "Gribov copies in the minimal Landau gauge: the influence on gluon and ghost propagators," *Nucl. Phys. B* **508** (1997) 353–370. DOI: 10.1016/S0550-3213(97)00629-9. arXiv:hep-lat/9705005. **`[V]` VERIFIED.**

Pure SU(2), minimal Landau gauge, β = 0, 0.8, 1.6, 2.7, volumes up to 16⁴ (24⁴ at β = 1.6). For the gluon propagator, Gribov noise is about the size of numerical accuracy even at small β; for the ghost propagator the noise is clearly visible at the three strong-coupling β values, and minimal-Landau-gauge data are always below generic-Landau-gauge data. This is a quantitative warning: **gauge-fixed propagators are not unique even after minimizing**, and comparisons of DSE solutions to lattice must account for the gauge-fixing prescription.

### E20. Cucchieri 1998 — the fundamental modular region
**A. Cucchieri**, "Numerical study of the fundamental modular region in the minimal Landau gauge," *Nucl. Phys. B* **521** (1998) 365–379. arXiv:hep-lat/9711024. **`[V]` VERIFIED.**

Studies the region Λ free of Gribov copies in pure SU(2) minimal Landau gauge, β = 0, 0.8, 1.6, 2.7, volumes up to 16⁴. Finds that typical thermalized gauge-fixed configurations — including those in Λ — lie very close to the Gribov horizon ∂Ω, and in the large-volume limit are characterized by a negative-definite horizon tensor. This supports the Gribov–Zwanziger picture but also shows that the continuum extrapolation of gauge-fixed quantities is delicate.

### E21. Leinweber–Skullerud–Williams–Parrinello 1998 — infrared structure survives cuts
**D. B. Leinweber, J. I. Skullerud, A. G. Williams, C. Parrinello** (UKQCD), "Gluon propagator in the infrared region," *Phys. Rev. D* **58** (1998) 031501. DOI: 10.1103/PhysRevD.58.031501. arXiv:hep-lat/9803015. **`[V]` VERIFIED.**

Quenched QCD on 16³×48 and 32³×64 at β = 6.0, Landau gauge; identifies finite-volume and anisotropy artifacts on the smaller lattice, applies cuts on the larger one. A residual infrared structure survives conservative cuts and rules out several models; a simple analytic fit is given. Still gauge-fixed and quenched; the result is evidence for IR dressing, not a physical mass.

---

## 6. Historical anchors (pre-1993, load-bearing for this era)

### E22. Taylor 1971 — the ghost–gluon Ward identity
**J. C. Taylor**, "Ward identities and charge renormalization of the Yang–Mills field," *Nucl. Phys. B* **33** (1971) 436–444. DOI: 10.1016/0550-3213(71)90297-5. **`[T]` VERIFIED.**

Generalized Ward–Takahashi identity implying a common Z₁/Z₂ for the vector and fictitious scalar fields, hence equal renormalized charges. The modern corollary used throughout the DSE literature: in Landau gauge the ghost–gluon vertex renormalization constant is finite (Taylor's non-renormalization theorem), which is what makes the scaling truncations of E14/E15 viable. Exact statement within the gauge-fixed theory.

### E23. Gribov 1978 — the ambiguity
**V. N. Gribov**, "Quantization of non-Abelian gauge theories," *Nucl. Phys. B* **139** (1978) 1–19. DOI: 10.1016/0550-3213(78)90175-X. **`[O]` VERIFIED.**

Fixing the divergence of the potential does not fix the gauge; restricting the functional integral to the region where the Faddeev–Popov determinant is positive cancels the infrared singularity of perturbation theory and produces a linear rising charge interaction at large distances. This is the origin of the Gribov horizon and of the physical motivation for the GZ action. Historical anchor; the "linear rising" statement is heuristic, not a theorem.

### E24. Singer 1978 — the topological no-go
**I. M. Singer**, "Some remarks on the Gribov ambiguity," *Comm. Math. Phys.* **60** (1978) 7–12. DOI: 10.1007/BF01609471. **`[T]` UNVERIFIED (standard citation; not re-checked this session because the search budget was exhausted).**

The rigorous statement behind E23: for a compact simple gauge group on a compactified space, no continuous gauge-fixing section exists over the gauge-orbit space; the Gribov ambiguity is topological, not a technical artifact. This is load-bearing for the gauge-fixing section below and **must be independently re-verified before v3 submission**.

### E25. Kugo–Ojima 1979 — the criterion
**T. Kugo, I. Ojima**, "Local covariant operator formalism of non-Abelian gauge theories and quark confinement problem," *Prog. Theor. Phys. Suppl.* **66** (1979) 1–130. DOI: 10.1143/PTPS.66.1; erratum *Prog. Theor. Phys.* **71** (1984) 1121. **`[C]` VERIFIED.**

The original BRST-quartet construction: physical S-matrix unitarity from the BRS charge algebra, observables as BRS-invariant operators, and the 1+u=0 color-confinement criterion derived from the Yang–Mills "Maxwell" equations. This is the framework E13 sharpens; it is a criterion inside a specific operator formalism, not a mass-gap theorem.

### E26. Cornwall 1982 — dynamical gluon mass
**J. M. Cornwall**, "Dynamical mass generation in continuum quantum chromodynamics," *Phys. Rev. D* **26** (1982) 1453–1478. DOI: 10.1103/PhysRevD.26.1453. **`[C]` VERIFIED.**

Derives a set of Schwinger–Dyson equations from resummed Feynman graphs for **gauge-invariant color-singlet** Green's functions (essential for multiplicative renormalizability in light-cone gauge). The equations are closed by solving a Ward identity exact in the IR but UV-ambiguous. Numerical results: effective gluon mass 500 ± 200 MeV; 0⁺ glueball mass about twice that. This is the era's phenomenological target and the "decoupling-like" branch in embryo; it is a model/resummation, not a theorem, and the light-cone gauge choice is part of the scheme.

---

## 7. What the scaling solution does and does not imply physically

**The setup.** Write the Landau-gauge propagators as
```
D(k²) = Z_D(k²)/k² ,      G(k²) = Z_G(k²)/k² .
```
The scaling ansatz is a power-law IR behavior,
```
Z_D(k²) ~ (k²)^{2κ} ,     Z_G(k²) ~ (k²)^{-κ} ,
```
so that
```
D(k²) ~ (k²)^{2κ-1} ,    G(k²) ~ (k²)^{-1-κ} .
```
For κ > 1/2 the gluon propagator vanishes at k=0 and the ghost propagator diverges. The running coupling built from the ghost–gluon vertex, α_s(k) ∝ Z_D Z_G², tends to a finite IR fixed point; the 1997 solution gives α_c ≈ 9.5. (The precise numerical exponent depends on the truncation and on the renormalization scheme; v3 should not quote a κ value without re-deriving it from the cited truncation.)

**What it does imply.**
- A self-consistent IR solution of the truncated Landau-gauge DSE system exists. That is a non-trivial mathematical fact about the truncated equations.
- The gluon propagator violates positivity (it vanishes rather than having a Källén–Lehmann pole). This is consistent with the gluon not being an asymptotic state, i.e. with a confinement scenario.
- The ghost sector dominates the IR; the solution realizes "ghost dominance."
- The product α_s(k) has a finite IR limit in this truncation, which is a statement about the gauge-fixed running coupling, not about a physical coupling.

**What it does not imply.**
- It does **not** prove confinement. Positivity violation of a gauge-fixed propagator is necessary-ish evidence, not a proof of area law or of a confining spectrum.
- It does **not** prove a mass gap. A vanishing gluon propagator is not a mass; the physical spectrum consists of gauge-invariant glueballs, and their gap is not computed by the scaling solution. Conversely, a mass gap can coexist with a vanishing or a massive gluon propagator.
- It does **not** establish that this solution is the physical one. The truncation drops the four-gluon and ghost–gluon scattering kernels; the solution is not unique; different boundary conditions/truncations can give the decoupling branch.
- It does **not** survive as a gauge-invariant statement. D(k), G(k), and α_s(k) are Landau-gauge quantities; their IR behavior can change with the gauge-fixing prescription (see E19).
- It does **not** decide the era's scaling-vs-decoupling question. In-era lattice fits (E18, E21) favored a massive/decoupling-like form; the analytic scaling solution and the lattice data were in tension, and the tension was not resolved within 1993–2000.

**Post-era context (outside this era; do not back-date).** Later lattice and functional work established that in 4D Landau-gauge SU(N) the *decoupling* solution is the one realized on the lattice (finite gluon propagator at k=0, finite ghost dressing), while the scaling solution remains a valid solution of the scaling DSE system under different boundary conditions. The v3 literature base should present both branches and never cite "scaling" as established 4D physics.

---

## 8. The gauge-fixing issue: physical observables vs Landau-gauge propagators

The Clay mass gap is a statement about the **physical Hilbert space** (OS/Wightman), i.e. about gauge-invariant correlators. Almost everything in this era is a statement about **gauge-fixed fields**. The distinction must be explicit in v3.

| Statement | Framework | Gauge status | Physical? |
|---|---|---|---|
| Wilson-loop area law / string tension σ | lattice, gauge-invariant | gauge-independent | **yes** |
| Glueball masses (0⁺⁺, 2⁺⁺, …) | gauge-invariant correlators | gauge-independent | **yes** |
| Gluon propagator D(k), ghost propagator G(k) | Landau-gauge DSE / lattice | gauge-dependent; Gribov-copy dependent | **no**; order parameter at best |
| Dressing functions Z_D, Z_G and α_s(k) IR fixed point | Landau-gauge DSE | gauge-dependent | **no** |
| Kugo–Ojima u, criterion 1+u=0 | BRST operator formalism, covariant gauges | gauge-parameter independent in covariant gauges | criterion for **color confinement**, not the mass gap |
| Gribov–Zwanziger action, γ, complex poles | Landau gauge + horizon restriction | gauge-fixed by construction | **no direct**; glueball channel has a physical cut |
| ERG Γ_k, finite-k gluon mass | gauge-fixed/background-field flow | mSTI-controlled; BRST broken at finite k | **no at finite k** |
| Cornwall m_g | gauge-invariant resummation, light-cone gauge | designed gauge-invariant | **model-dependent**, not a theorem |

**Consequences for v3.**
1. A DSE statement of the form "the gluon acquires a mass m" is, literally, a statement about the transverse part of a gauge-fixed two-point function. It cannot be inserted into an OS/Wightman mass-gap proof without a gauge-invariant bridge.
2. The only DSE-era quantity that is even a *candidate* gauge-invariant order parameter is Cornwall's resummed color-singlet sector; it is a model, not a derivation.
3. The Kugo–Ojima criterion is the one place where the era attempted a gauge-parameter-independent statement, and it is about **color confinement**, which is logically distinct from the mass gap: confinement plus unbroken chiral symmetry can coexist with massless hadrons in principle; a mass gap needs a spectral statement, not only a BRST-cohomology statement.
4. The Gribov–Zwanziger action is a gauge-fixed action with an explicit dimensionful parameter; its "mass" is not the glueball mass. The 1989 observation that the glueball channel has a physical cut is the correct place to look for a physical spectrum, and it was not developed into a mass-gap computation in this era.
5. The ERG finite-k mass is the most dangerous: E8 shows it is nonzero whenever k ≠ 0, and E7 shows it can be tuned away by STI boundary conditions. A v3 argument that reads a mass off Γ_k at finite k is reading a regulator artifact.
6. Lattice propagator masses (E18, E21) are gauge-fixed, quenched, finite-volume, and (in this era) single-lattice-spacing. They are evidence, not a continuum physical mass.

---

## 9. Position in 14 eras

1. **1954–1969 — Foundations.** Yang–Mills; Faddeev–Popov; canonical quantization and the first ghosts.
2. **1970–1973 — Renormalization and asymptotic freedom.** 't Hooft–Veltman; Politzer; Gross–Wilczek; instantons.
3. **1974–1978 — Lattice gauge theory and strong-coupling confinement.** Wilson; Kogut–Susskind; Gribov ambiguity; Singer's topological no-go.
4. **1979–1984 — Semiclassical and constructive.** Savvidy vacuum; Kugo–Ojima; Cornwall's dynamical mass; Balaban's RG; Creutz Monte Carlo.
5. **1985–1989 — Gribov–Zwanziger action; BRST; first lattice propagators.** Local horizon action; DSE revival; early gluon-mass fits.
6. **1990–1992 — Infrared DSE and ERG precursors.** Wetterich's average action; Reuter's gauge flow; Magnen–Rivasseau–Sénéor; first scaling ansätze.
7. **1993–2000 — Functional methods (THIS ERA).** Exact ERG equations; Ellwanger–Hirsch–Weber; von Smekal–Hauck–Alkofer scaling solution; Gribov–Zwanziger renormalization; lattice massive-gluon evidence; Kugo–Ojima revisited; scaling-vs-decoupling seeds.
8. **2001–2004 — Reviews and refinement.** Alkofer–von Smekal review; refined Gribov–Zwanziger; larger lattices; the scaling/decoupling debate sharpens.
9. **2005–2009 — Vertex DSEs and FRG unification.** IR-safe vertices; ghost–gluon vertex lattice; functional RG and DSE cross-checks; Gribov–Zwanziger refinements.
10. **2010–2013 — Decoupling established on the lattice.** Massive gluon; continuum extrapolations; refined GZ action; vertex modeling.
11. **2014–2017 — Systematics and truncation error.** DSE truncation systematics; Landau-gauge uniqueness; lattice continuum limits; bootstrap seeds.
12. **2018–2021 — Rigorous Gribov–Zwanziger and stochastic quantization.** 2D/3D constructions; area-law expansions; mathematically controlled horizon.
13. **2022–2024 — 3D YM–Higgs and strong-coupling stochastic analysis.** Regularity structures; Langevin; conformal-window bootstrap attacks.
14. **2025–2026 — Formalization era.** Lean 4D construction; OSforGFF; MG reduction; AI-assisted proof checking.

**Era-7 handoff to era 8.** The exact equations and the scaling solution were in place; what was missing was (i) a gauge-invariant bridge, (ii) truncation control, (iii) the continuum/infinite-volume lattice extrapolation, and (iv) a decision between scaling and decoupling. Era 8 inherited all four.

---

## 10. Gaps and attack surfaces for v3

- **G1. No gauge-invariant mass.** The era's "gluon mass" is a gauge-fixed two-point quantity. v3 needs either a gauge-invariant definition (glueball pole, Wilson-loop correlator) or an explicit conditional lemma of the form "if the Landau-gauge propagator has property P, then the physical spectrum has gap Δ." That lemma does not exist in this literature.
- **G2. Truncation control is absent.** The scaling solution drops four-gluon and ghost–gluon scattering kernels. There is no a posteriori error bound. v3 should either (a) reproduce the truncation and quantify the neglected kernels, or (b) treat scaling as one branch of a family and prove a statement stable across the family.
- **G3. Scaling vs decoupling is undecided in-era.** The lattice fits (E18, E21) and the analytic solution (E14/E15) point in different directions. v3 should not choose silently; it should formalize both branches and identify which (if either) can carry a mass-gap inequality. Post-era lattice favors decoupling; that must be flagged as post-era.
- **G4. ERG gauge invariance is a limit statement.** The mSTI machinery (E2, E3, E7, E8) and the Litim–Pawlowski obstruction (E9) mean that a finite-k ERG mass is not physical. v3 should either construct a gauge-invariant flow or prove that the k→0 limit of the physical observable is independent of the gauge-fixing/background choice.
- **G5. Kugo–Ojima is not a mass-gap theorem.** 1+u=0 is a criterion for color confinement in a specific operator formalism. v3 must not conflate confinement with a mass gap; if it uses KO, it needs an additional spectral lemma.
- **G6. Gribov–Zwanziger's γ is not the mass gap.** The GZ action is gauge-fixed and its dimensionful parameter is tied to the horizon. v3 should use the 1989 observation that the glueball channel has a physical cut as the correct gauge-invariant target, and should not identify γ^{1/2} with a physical mass.
- **G7. Lattice gauge-fixing systematics.** E19/E20 show that Gribov copies contaminate ghost (and, at strong coupling, gluon) propagators, and that thermalized configurations sit near the horizon. Any v3 use of lattice propagator masses must state the gauge-fixing prescription and the copy-handling, or use gauge-invariant observables.
- **G8. Continuum and infinite-volume limits.** E18/E21 are quenched, finite-volume, single-lattice-spacing. A "massive gluon" claim needs a continuum extrapolation and a volume study; neither exists in-era. v3 should treat the lattice mass as a fitted scale, not a physical parameter.
- **G9. Missing rigorous bridge.** No theorem in 1993–2000 connects any DSE/ERG output to an OS/Wightman mass gap. The exact missing step is an implication chain: *gauge-fixed functional solution* → *gauge-invariant correlator* → *exponential clustering with a volume-uniform rate* → *spectral gap*. v3 should formalize each arrow and mark the unproved ones.
- **G10. Positivity and the physical Hilbert space.** The scaling solution's gluon violates positivity; that is evidence for confinement, but the mass gap is a statement about the positive physical subspace. v3 should state explicitly which Hilbert space the gap is claimed in, and whether the indefinite-metric gauge-fixed description is being used only as an intermediate.

---

## 11. Verification ledger

| # | Entry (short) | Checked against | Status |
|---|---|---|---|
| E1 | Wetterich, PLB 301 (1993) 90 | WebSearch (Elsevier/Google Scholar), arXiv:1710.05815 | VERIFIED |
| E2 | Reuter–Wetterich, NPB 417 (1994) 181 | WebSearch (ScienceDirect, ADS 1994NuPhB.417..181R) | VERIFIED |
| E3 | Bonini–D'Attanasio–Marchesini, NPB 421 (1994) 429 | INSPIRE-HEP API (recid 361103; arXiv:hep-th/9312114) | VERIFIED |
| E4 | Ellwanger, ZPC 62 (1994) 503 | WebSearch (Springer, arXiv:hep-ph/9308260) | VERIFIED |
| E5 | Morris, IJMPA 9 (1994) 2411 | INSPIRE-HEP API (recid 358382; arXiv:hep-ph/9308265) | VERIFIED |
| E6 | Wetterich, IJMPA 9 (1994) 3571 | INSPIRE-HEP API (recid 35464) | VERIFIED |
| E7 | Bonini–D'Attanasio–Marchesini, NPB 437 (1995) 163 | arXiv:hep-th/9410138 (abs page) | VERIFIED |
| E8 | Ellwanger–Hirsch–Weber, ZPC 69 (1996) 687 | WebSearch (Springer/Scopus) + arXiv:hep-th/9506019 (abs page) | VERIFIED |
| E9 | Litim–Pawlowski, hep-th/9901063 | WebSearch (arXiv abs, ADS) | VERIFIED |
| E10 | Zwanziger, NPB 323 (1989) 513 | WebSearch (ScienceDirect, INSPIRE) | VERIFIED |
| E11 | Zwanziger, NPB 399 (1993) 477 | WebSearch (ScienceDirect, ADS 1993NuPhB.399..477Z) | VERIFIED |
| E12 | Zwanziger, NPB 412 (1994) 657 | WebSearch (CERN CDS, INSPIRE; page range corrected from 127) | VERIFIED |
| E13 | Kugo, hep-th/9511033 (1995) | WebSearch (arXiv abs, Semantic Scholar) | VERIFIED |
| E14 | von Smekal–Hauck–Alkofer, PRL 79 (1997) 3591 | WebSearch (OSTI, APS, INSPIRE; arXiv:hep-ph/9705242) | VERIFIED |
| E15 | von Smekal–Hauck–Alkofer, Annals Phys. 267 (1998) 1 | INSPIRE-HEP API (recid 439877; erratum 269 (1998) 182) | VERIFIED |
| E16 | Roberts–Williams, PP NP 33 (1994) 477 | WebSearch (Google Scholar, OUCI, CiNii) | VERIFIED |
| E17 | Alkofer–von Smekal, Phys. Rept. 353 (2001) 281 | WebSearch (arXiv:hep-ph/0007355, INSPIRE) | VERIFIED |
| E18 | Marenzoni–Martinelli–Stella, NPB 455 (1995) 339 | WebSearch (ScienceDirect, OSTI, Scopus; arXiv:hep-lat/9410011) | VERIFIED |
| E19 | Cucchieri, NPB 508 (1997) 353 | WebSearch (arXiv:hep-lat/9705005, Scopus, INSPIRE) | VERIFIED |
| E20 | Cucchieri, NPB 521 (1998) 365 | WebSearch (arXiv:hep-lat/9711024, ScienceDirect; volume corrected from 508) | VERIFIED |
| E21 | Leinweber–Skullerud–Williams–Parrinello, PRD 58 (1998) 031501 | WebSearch (ADS 1998PhRvD..58c1501L, APS; arXiv:hep-lat/9803015) | VERIFIED |
| E22 | Taylor, NPB 33 (1971) 436 | WebSearch (MathSciNet MR0314386, Mendeley, Semantic Scholar) | VERIFIED |
| E23 | Gribov, NPB 139 (1978) 1 | WebSearch (ADS 1978NuPhB.139....1G, DOI) | VERIFIED |
| E24 | Singer, CMP 60 (1978) 7 | not re-checked (search budget exhausted) | UNVERIFIED |
| E25 | Kugo–Ojima, PTPS 66 (1979) 1 | WebSearch (OUP, CERN CDS, INIS) | VERIFIED |
| E26 | Cornwall, PRD 26 (1982) 1453 | WebSearch (ADS 1982PhRvD..26.1453C, APS, Scopus) | VERIFIED |

**Coverage:** 26 entries, 25 VERIFIED, 1 UNVERIFIED (E24, Singer). Labels used: `[T]` E1, E2, E6, E11, E22, E24; `[C]` E3, E4, E7, E8, E10, E12, E13, E14, E15, E25, E26; `[NG]` E9; `[O]` E5, E16, E17, E23; `[V]` E18, E19, E20, E21.

**Known bibliographic corrections captured during verification (do not propagate the erroneous forms):**
- Zwanziger's 1994 FMR/area-law paper is *Nucl. Phys. B* **412** (1994) 657–730, not p. 127 (127 is the preprint page count).
- Cucchieri's fundamental-modular-region paper is *Nucl. Phys. B* **521** (1998) 365–379; the 508/353 slot belongs to the Gribov-copies/propagators paper.
- Bonini–D'Attanasio–Marchesini: *Phys. Lett. B* **329** (1994) 249 is "Axial anomalies…"; the BRS/ERG paper is *Nucl. Phys. B* **437** (1995) 163. The SU(2) flow paper is *Nucl. Phys. B* **421** (1994) 429 (arXiv:hep-th/9312114, not hep-th/9311142).
- von Smekal–Hauck–Alkofer have *two* load-bearing records: the PRL 79 (1997) 3591 letter and the Annals Phys. 267 (1998) 1 long paper (with erratum 269 (1998) 182). Cite both.
- Leinweber et al. author order is Leinweber, Skullerud, Williams, Parrinello (not Parrinello first).

---

## 12. One-paragraph bottom line for v3

The 1993–2000 functional-method era is best read as a **precision statement of the obstruction**, not as a step toward an unconditional theorem. It gives v3 (i) an exact ERG equation whose gauge-fixed status is controlled by modified Slavnov–Taylor identities, (ii) a self-consistent Landau-gauge DSE scaling solution whose physical content is a positivity-violating gluon and a finite IR α_s, (iii) a renormalizable Gribov–Zwanziger action with a horizon parameter that is not the mass gap, and (iv) lattice propagator fits that are massive-looking but gauge-fixed, quenched, and finite-volume. None of these is a gauge-invariant mass-gap statement. The honest era-7 verdict is therefore: **a scaling solution exists, the Gribov–Zwanziger action exists, the lattice sees IR dressing, and the bridge from any of them to the physical spectrum is exactly the missing step.** That bridge — a gauge-invariant, truncation-controlled, continuum-and-infinite-volume statement — is the object v3 must either build or prove impossible.
