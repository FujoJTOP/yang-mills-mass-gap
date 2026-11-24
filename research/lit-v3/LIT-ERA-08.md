# LIT-ERA-08 · 1998–2005 — Holography, large-N, and AdS/QCD

> Literature base for paper v3. Era 8 of 14. Compiled 2026-11-24.
> Scope: Maldacena's AdS/CFT; Gubser–Klebanov–Polyakov; Witten; thermal AdS and confinement; the holographic Wilson loop and the area law; Polchinski–Strassler and the soft wall; Sakai–Sugimoto and AdS/QCD; hard-wall models and the glueball spectrum; the large-N master field and loop equations (Makeenko–Migdal, Eguchi–Kawai); the status of a 4D pure-YM string dual.
> Method: primary sources checked against arXiv abstract pages and INSPIRE-HEP records. **VERIFIED** below means title/authors/venue/DOI were checked against those records; it does **not** certify the scientific content. No bibliographic datum is invented; where a publisher DOI could not be confirmed it is omitted.

---

## 0. Verdict in one paragraph

Holography does **not** prove a mass gap for 4D pure SU(N) Yang–Mills, and the 1998–2005 corpus does not contain the ingredients for such a proof. What the era contains is (i) a precise, extensively checked **conjecture** (AdS/CFT) for supersymmetric conformal theories; (ii) a **mechanism** by which a confining background produces an area law and a gapped fluctuation spectrum; (iii) a family of **models** (hard wall, soft wall, Sakai–Sugimoto) that reproduce glueball and meson spectra by inserting an IR scale by hand or by compactification. None is a dual of 4D pure YM, and none supplies a controlled derivation of Δ > 0 for SU(3). The exact large-N structures of the era (planar expansion, loop equations, Eguchi–Kawai reduction) are theorem-grade in places but do not close in four dimensions. The honest label for the era is: **one conjecture, one mechanism, many models; no theorem about the Clay gap.**

---

## 1. Labels

| Label | Meaning |
|---|---|
| **[T]** | Theorem / exact identity: mathematical proof, or an exact equation/identity. |
| **[C]** | Conjecture: well-posed but unproven. |
| **[NG]** | No-go / obstruction: a result or structural argument that a class of approaches cannot work as stated. |
| **[O]** | Observation / model / phenomenological computation: standard physics derivation, model fit, or strong-coupling computation inside a conjectured duality; not a proof. |
| **[V]** | Verified numerical / lattice result. |

A single entry carries one primary label; secondary aspects are noted in its annotation.

---

## 2. Annotated entries (35)

### A. The duality and its consistency (5)

**A1. [C]** J. M. Maldacena, "The Large N Limit of Superconformal Field Theories and Supergravity," *Adv. Theor. Math. Phys.* **2** (1998) 231–252; arXiv:hep-th/9711200; DOI 10.1023/A:1026654312961. **VERIFIED.**
The founding statement: N=4 SU(N) SYM at large N and large λ is dual to type IIB string theory on AdS₅×S⁵. It is a conjecture; the evidence is symmetry matching, the supergravity spectrum, and later precision checks. It is **not** a statement about pure YM.

**A2. [C]** S. S. Gubser, I. R. Klebanov, A. M. Polyakov, "Gauge Theory Correlators from Non-Critical String Theory," *Phys. Lett. B* **428** (1998) 105–114; arXiv:hep-th/9802109; DOI 10.1016/S0370-2693(98)00377-3. **VERIFIED.**
The GKPW prescription: the generating functional of boundary correlators equals the string partition function with boundary conditions, `Z_CFT[J] = Z_string[φ→J]`. This is the operational content of the duality; a conjecture, checked in the supergravity limit.

**A3. [C]** E. Witten, "Anti de Sitter Space and Holography," *Adv. Theor. Math. Phys.* **2** (1998) 253–291; arXiv:hep-th/9802150. **VERIFIED.**
Sharpens the dictionary, discusses the holographic bound, and gives the boundary-to-bulk map for fields of arbitrary spin. Foundational conjecture.

**A4. [O]** O. Aharony, S. S. Gubser, J. Maldacena, H. Ooguri, Y. Oz, "Large N Field Theories, String Theory and Gravity," *Phys. Rept.* **323** (2000) 183–386; arXiv:hep-th/9905111; DOI 10.1016/S0370-1573(99)00083-6. **VERIFIED.**
The standard review of the conjecture and its checks; the era's map. A review, not a proof.

**A5. [T]** P. Breitenlohner, D. Z. Freedman, "Positive Energy in anti-De Sitter Backgrounds and Gauged Extended Supergravity," *Phys. Lett. B* **115** (1982) 197–201; DOI 10.1016/0370-2693(82)90643-8. **VERIFIED.**
The BF bound: a scalar of mass m in AdS_{d+1} is stable if m² ≥ −d²/4. Theorem-grade stability fact that makes AdS/CFT with relevant operators possible; without it the "conformal boundary + normalizable modes" picture fails. Pre-era but load-bearing.

### B. Thermal AdS, Wilson loops, and the dictionary (7)

**B6. [C]** E. Witten, "Anti-de Sitter Space, Thermal Phase Transition, and Confinement in Gauge Theories," *Adv. Theor. Math. Phys.* **2** (1998) 505–532; arXiv:hep-th/9803131. **VERIFIED.**
The Hawking–Page transition in AdS is identified with confinement/deconfinement in the dual gauge theory: the low-temperature thermal-AdS phase has a mass gap and an area law, the high-temperature black-hole phase does not. The same paper proposes a non-supersymmetric 4D gauge theory from D4-branes on a circle with antiperiodic fermions (the "Witten model"). The gravity-side transition is a classical GR computation [T]; the identification with a **4D pure-YM** gap is [C]/[O], because the model is a 5D theory on a circle.

**B7. [C]** J. M. Maldacena, "Wilson loops in large N field theories," *Phys. Rev. Lett.* **80** (1998) 4859–4862; arXiv:hep-th/9803002; DOI 10.1103/PhysRevLett.80.4859. **VERIFIED.**
The holographic Wilson loop: the quark–antiquark potential is the minimal string worldsheet ending on the boundary. At zero T in N=4 the potential is Coulombic (conformal), not confining; at finite T the black-hole geometry screens. The "area law" is a strong-coupling computation in a specific background, not a proof of confinement for YM.

**B8. [C]** S.-J. Rey, J.-T. Yee, "Macroscopic strings as heavy quarks: Large-N gauge theory and anti-de Sitter supergravity," *Eur. Phys. J. C* **22** (2001) 379–394; arXiv:hep-th/9803001; DOI 10.1007/s100520100799. **VERIFIED.**
Independent derivation of the Wilson-loop/string correspondence, with the same caveat.

**B9. [C]** A. Brandhuber, N. Itzhaki, J. Sonnenschein, S. Yankielowicz, "Wilson Loops in the Large N Limit at Finite Temperature," *Phys. Lett. B* **434** (1998) 36–40; arXiv:hep-th/9803137; DOI 10.1016/S0370-2693(98)00730-8. **VERIFIED.**
Finite-temperature Wilson loops: the deconfined phase screens, the confined phase area-laws. Model computation.

**B10. [C]** A. M. Polyakov, "String Theory and Quark Confinement," *Nucl. Phys. B Proc. Suppl.* **68** (1998) 1–8; arXiv:hep-th/9711002; DOI 10.1016/S0920-5632(98)00135-2. **VERIFIED.**
Polyakov's proposal that confinement is a string phenomenon and that the QCD string is a dual string; no explicit 4D pure-YM dual is constructed. The conceptual bridge from 't Hooft's planar expansion to holography.

**B11. [C]** E. Witten, "Baryons and Branes in Anti de Sitter Space," *JHEP* **9807** (1998) 006; arXiv:hep-th/9805112; DOI 10.1088/1126-6708/1998/07/006. **VERIFIED.**
Baryons as wrapped D-branes; the dictionary for baryon operators. Relevant because the Sakai–Sugimoto construction uses wrapped D-branes and because it shows the dictionary extends beyond local operators.

**B12. [C]** D. J. Gross, H. Ooguri, "Aspects of Large N Gauge Theory Dynamics as Seen by String Theory," *Phys. Rev. D* **58** (1998) 106002; arXiv:hep-th/9805129; DOI 10.1103/PhysRevD.58.106002. **VERIFIED.**
Large-N gauge dynamics (confinement, chiral symmetry, the string tension) read through the duality; a bridge between the CFT and confining regimes. Conjectural.

### C. Glueballs and confining geometries (7)

**C13. [O]** C. Csáki, H. Ooguri, Y. Oz, J. Terning, "Glueball Mass Spectrum From Supergravity," *JHEP* **9901** (1999) 017; arXiv:hep-th/9806021; DOI 10.1088/1126-6708/1999/01/017. **VERIFIED.**
Computes 0⁺⁺, 2⁺⁺, … glueball masses from the supergravity fluctuation spectrum of a confining background; the prototype "glueball masses from holography" computation. Supergravity-level [O], inside the conjecture.

**C14. [O]** R. de Mello Koch, A. Jevicki, M. Mihailescu, J. P. Nunes, "Evaluation Of Glueball Masses From Supergravity," *Phys. Rev. D* **58** (1998) 105009; arXiv:hep-th/9806125; DOI 10.1103/PhysRevD.58.105009. **VERIFIED.**
Independent evaluation, same status.

**C15. [O]** R. C. Brower, S. D. Mathur, C.-I Tan, "Glueball Spectrum for QCD from AdS Supergravity Duality," *Nucl. Phys. B* **587** (2000) 249–276; arXiv:hep-th/0003115; DOI 10.1016/S0550-3213(00)00435-1. **VERIFIED.**
A systematic hard-wall-type computation of the glueball spectrum, compared to lattice ratios. Model.

**C16. [C]** I. R. Klebanov, E. Witten, "Superconformal Field Theory on Threebranes at a Calabi-Yau Singularity," *Nucl. Phys. B* **536** (1998) 199–218; arXiv:hep-th/9807080; DOI 10.1016/S0550-3213(98)00654-3. **VERIFIED.**
The conifold theory: N=1 SUSY, a conformal fixed point dual to AdS₅×T^{1,1}. The seed of the cascading construction.

**C17. [C]** I. R. Klebanov, A. A. Tseytlin, "Gravity Duals of Supersymmetric SU(N)×SU(N+M) Gauge Theories," *Nucl. Phys. B* **578** (2000) 123–138; arXiv:hep-th/0002159; DOI 10.1016/S0550-3213(00)00206-6. **VERIFIED.**
The KT solution with 3-form flux; the geometry that cascades.

**C18. [C]** I. R. Klebanov, M. J. Strassler, "Supergravity and a Confining Gauge Theory: Duality Cascades and χSB-Resolution of Naked Singularities," *JHEP* **0008** (2000) 052; arXiv:hep-th/0007191; DOI 10.1088/1126-6708/2000/08/052. **VERIFIED.**
The KS cascade: a confining N=1 SUSY gauge theory with a mass gap and chiral symmetry breaking, dual to a warped deformed conifold with a smooth IR cap (the naked singularity is resolved). The closest thing to a "derived" confining mass gap in the era — but it is N=1 SUSY with a cascade, not pure YM.

**C19. [C]** J. M. Maldacena, C. Nuñez, "Towards the large N limit of pure N=1 super Yang Mills," *Phys. Rev. Lett.* **86** (2001) 588–591; arXiv:hep-th/0008001; DOI 10.1103/PhysRevLett.86.588. **VERIFIED.**
A supergravity dual of N=1 SYM (with KK modes); the theory confines and has a mass gap. Still not pure YM, and the dual has a singularity/compactification.

### D. AdS/QCD models (6)

**D20. [O]** J. Polchinski, M. J. Strassler, "Hard Scattering and Gauge/String Duality," *Phys. Rev. Lett.* **88** (2002) 031601; arXiv:hep-th/0109174; DOI 10.1103/PhysRevLett.88.031601. **VERIFIED.**
Introduces a hard-wall IR cutoff in AdS to break conformal invariance and produce scaling; the origin of the hard-wall AdS/QCD model. Model, not a derivation.

**D21. [O]** T. Sakai, S. Sugimoto, "Low energy hadron physics in holographic QCD," *Prog. Theor. Phys.* **113** (2005) 843–882; arXiv:hep-th/0412141; DOI 10.1143/PTP.113.843. **VERIFIED.**
The D4/D8 top-down model (based on Witten's D4 background): a 5D gauge theory on a circle with flavor branes; reproduces mesons, baryons, and chiral symmetry breaking. The closest top-down "QCD-like" construction; it is a 5D theory on a circle, not 4D pure YM, and the duality is a conjecture.

**D22. [O]** T. Sakai, S. Sugimoto, "More on a holographic dual of QCD," *Prog. Theor. Phys.* **114** (2005) 1083–1118; arXiv:hep-th/0507073; DOI 10.1143/PTP.114.1083. **VERIFIED.**
Extends the model (vector/axial-vector mesons, baryons, etc.). Model.

**D23. [O]** J. Erlich, E. Katz, D. T. Son, M. A. Stephanov, "QCD and a Holographic Model of Hadrons," *Phys. Rev. Lett.* **95** (2005) 261602; arXiv:hep-ph/0501128; DOI 10.1103/PhysRevLett.95.261602. **VERIFIED.**
The bottom-up hard-wall model: AdS₅ cut at z₀, with a bulk scalar for chiral symmetry breaking; fits meson masses. Parameters are fit, so the "gap" is an input.

**D24. [O]** L. Da Rold, A. Pomarol, "Chiral symmetry breaking from five dimensional spaces," *Nucl. Phys. B* **721** (2005) 79–97; arXiv:hep-ph/0501218; DOI 10.1016/j.nuclphysb.2005.05.009. **VERIFIED.**
A related hard-wall/warped model with chiral symmetry breaking; model.

**D25. [O]** A. Karch, E. Katz, D. T. Son, M. A. Stephanov, "Linear Confinement and AdS/QCD," *Phys. Rev. D* **74** (2006) 015005; arXiv:hep-ph/0602229; DOI 10.1103/PhysRevD.74.015005. **VERIFIED.**
The soft-wall model: a dilaton-like background produces m² ∝ n (linear Regge trajectories). The smooth IR cutoff is the key feature; see the NEC discussion in §5. Published 2006 (preprint Feb 2006): the bridge to era 9.

### E. Large-N exact structures (5)

**E26. [C]** G. 't Hooft, "A Planar Diagram Theory for Strong Interactions," *Nucl. Phys. B* **72** (1974) 461–473; DOI 10.1016/0550-3213(74)90154-0. **VERIFIED.**
The 1/N expansion: planar diagrams dominate at fixed λ = g²N, and the topological expansion resembles a string coupling expansion. The combinatorial planar-dominance statement is theorem-grade; the identification with a string dual is a conjecture. The seed of the holographic idea.

**E27. [T]** Yu. M. Makeenko, A. A. Migdal, "Exact equation for the loop average in multicolor QCD," *Phys. Lett. B* **88** (1979) 135–137; DOI 10.1016/0370-2693(79)90131-X; Erratum *Phys. Lett. B* **89** (1980) 437. **VERIFIED.**
The Makeenko–Migdal loop equation: a closed equation for the Wilson-loop average at N=∞. Exact identity; solvable in 2D, not solved in 4D.

**E28. [T]** Yu. M. Makeenko, A. A. Migdal, "Quantum Chromodynamics as Dynamics of Loops," *Nucl. Phys. B* **188** (1981) 269; DOI 10.1016/0550-3213(81)90258-3. **VERIFIED (page end not in INSPIRE record).**
The loop-space formulation; the master-field/loop-equation programme. Exact but non-closing in 4D.

**E29. [T]** T. Eguchi, H. Kawai, "Reduction of Dynamical Degrees of Freedom in the Large-N Gauge Theory," *Phys. Rev. Lett.* **48** (1982) 1063–1066; DOI 10.1103/PhysRevLett.48.1063. **VERIFIED.**
Large-N volume reduction: at N=∞ the lattice gauge theory on a finite volume is equivalent to a single-site matrix model, **provided the center symmetry is unbroken**. A theorem under stated conditions; the reduced model is still a matrix integral, not a solution.

**E30. [NG]** G. Bhanot, U. M. Heller, H. Neuberger, "The Quenched Eguchi-Kawai Model," *Phys. Lett. B* **113** (1982) 47–50; DOI 10.1016/0370-2693(82)90106-X. **VERIFIED.**
Introduces the quenched modification after the naive Eguchi–Kawai reduction is found to fail for D>2 because the center symmetry breaks; the quenched version restores it. A structural obstruction to using reduction as a shortcut.

### F. Holographic renormalization (3)

**F31. [O]** M. Henningson, K. Skenderis, "The Holographic Weyl Anomaly," *JHEP* **9807** (1998) 023; arXiv:hep-th/9806087; DOI 10.1088/1126-6708/1998/07/023. **VERIFIED.**
Holographic renormalization: the boundary Weyl anomaly is computed from the bulk; a consistency condition on candidate duals. Standard computation.

**F32. [O]** V. Balasubramanian, P. Kraus, "A Stress Tensor for Anti-de Sitter Gravity," *Commun. Math. Phys.* **208** (1999) 413–428; arXiv:hep-th/9902121; DOI 10.1007/s002200050764. **VERIFIED.**
The boundary stress tensor and conserved charges; a consistency tool.

**F33. [T]** E. Witten, S.-T. Yau, "Connectedness of the Boundary in the AdS/CFT Correspondence," *Adv. Theor. Math. Phys.* **3** (1999) 1635–1655; arXiv:hep-th/9910245. **VERIFIED.**
Proves, under stated hypotheses, that the conformal boundary of a conformally compact Einstein manifold is connected. Rules out candidate duals with disconnected boundaries; a rigorous consistency result.

### G. Benchmarks and IR consistency (2)

**G34. [V]** C. J. Morningstar, M. Peardon, "The glueball spectrum from an anisotropic lattice study," *Phys. Rev. D* **60** (1999) 034509; arXiv:hep-lat/9901004; DOI 10.1103/PhysRevD.60.034509. **VERIFIED.**
The lattice glueball spectrum (0⁺⁺, 2⁺⁺, …) used as the benchmark against which AdS/QCD model spectra are compared. Numerical, not a proof.

**G35. [O]** U. Gürsoy, E. Kiritsis, F. Nitti, "Exploring improved holographic theories for QCD: Part I," *JHEP* **0802** (2008) 032; arXiv:0707.1324; DOI 10.1088/1126-6708/2008/02/032; "Part II," *JHEP* **0802** (2008) 019; arXiv:0707.1349; DOI 10.1088/1126-6708/2008/02/019. **VERIFIED.**
Systematic analysis of the IR asymptotics of holographic QCD and the constraints (including the NEC) on backgrounds that produce linear confinement and a mass gap. Outside the era, included as the cleanest statement of the NEC/IR-throat obstruction.

---

## 3. Does holography prove a mass gap for 4D pure SU(N)?

### 3.1 Short answer

**No.** In 1998–2005 holography provides a **mechanism** and a family of **models**, not a proof. The only exact dualities are for supersymmetric theories (N=4 SYM; N=1 cascades). Pure SU(N) YM is asymptotically free, non-supersymmetric, and confining; no string dual of it is known. The AdS/QCD constructions (hard wall, soft wall, Sakai–Sugimoto) reproduce spectra but insert the IR scale by hand or by compactification. The large-N structures (planar expansion, loop equations, Eguchi–Kawai) are exact in places but do not close in 4D.

### 3.2 What "prove a mass gap" means

Precisely: construct the Euclidean measure dμ on the space of SU(N) connections, show it satisfies the OS axioms, and show that the transfer matrix H on the physical Hilbert space has `spec(H) ⊂ {0} ∪ [Δ, ∞)` with `Δ > 0` uniform in the infinite-volume limit. Equivalently, every connected gauge-invariant correlator decays exponentially at a uniform rate. This is the campaign's MG statement; it is a statement about **fixed N** (e.g. SU(3)), not about a large-N limit.

### 3.3 What a holographic proof would need

A holographic proof would need six things:

| # | Requirement | Content |
|---|---|---|
| **H1** | Duality theorem | An exact equality `Z_YM[J] = Z_string[boundary J]` for 4D pure SU(N) YM, not for N=4. At minimum a controlled large-N, large-λ expansion with proven error control. |
| **H2** | Legitimate background | A weakly curved (`R/ℓ_s ≫ 1`, `g_s ≪ 1`), asymptotically AdS/holographic solution of the full string theory, compact internal space, normalizable spectrum, no IR pathology. |
| **H3** | NEC / IR completion | The background satisfies the null energy condition, or is a known consistent string background whose NEC violation is sourced by a controlled object (e.g. a tachyon condensate with a known completion). A wall/throat that caps the IR must be smooth or resolved — not a delta-function brane or a naked singularity. |
| **H4** | Gravity-side gap | The normalizable fluctuation spectrum has a positive lowest eigenvalue `Δ_grav > 0`, uniform in the IR cutoff and surviving decompactification to R⁴. |
| **H5** | Stability of the gap | The gap survives α′ and g_s corrections, 1/N corrections, and the extrapolation from large N, large λ to N=3, finite λ. Pure YM has **no SUSY nonrenormalization theorem** to protect it; this is the hardest step. |
| **H6** | Correct theory | The dual is of 4D pure SU(N) — not a 5D theory on a circle (Witten/Sakai–Sugimoto), not an N=1 cascade (KS/MN), not a model with a fitted IR scale (hard/soft wall). |

### 3.4 Which of H1–H6 the corpus supplies

| Req. | Status in 1998–2005 | Why |
|---|---|---|
| H1 | **No** | AdS/CFT is a conjecture; exact duals are N=4 and N=1 SUSY. No pure-YM dual. |
| H2 | **Partial** | D4-brane and KS backgrounds are legitimate solutions; the hard wall is not a solution; the soft wall is not a standard solution. |
| H3 | **No for soft wall; yes for D4/KS** | Soft wall violates the NEC or needs an uncontrolled tachyon condensate; D4/KS have smooth caps but are not pure YM. |
| H4 | **Supergravity level only** | The glueball computations [O] find positive masses; not uniform under decompactification in the required sense. |
| H5 | **No** | No protection; the large-N/large-λ extrapolation is uncontrolled. |
| H6 | **No** | Every construction is higher-dimensional, SUSY, or a fitted model. |

### 3.5 The status of a 4D pure-YM string dual

- **'t Hooft 1974:** planar diagrams suggest a string dual; no construction.
- **Polyakov 1998:** confinement as a string phenomenon; no explicit dual.
- **Maldacena 1998:** the dual is N=4 SYM; pure YM is not covered.
- **Witten 1998:** the D4-brane model gives a 4D gauge theory below the KK scale; it is a 5D theory on a circle and the duality is a conjecture.
- **Klebanov–Strassler / Maldacena–Nuñez:** confining N=1 SUSY duals; not pure YM.
- **AdS/QCD:** models, not duals.

**Conclusion.** As of 2005 there is **no known string dual of 4D pure SU(N) YM.** The existence of such a dual is a conjecture; the models are phenomenological. The precise missing step is **H1 + H5**: a controlled duality for the non-supersymmetric theory and a protection/continuation argument to N=3, finite λ.

### 3.6 The gap is not the same as a compactification gap

Two subtleties that are routinely blurred:

1. In Witten's thermal AdS and in the D4 model, the "mass gap" is often the **KK/compactification scale** 1/R. The 4D pure-YM limit requires R → 0 (or the KK scale far above Λ_QCD), a strong-coupling/large-N extrapolation; the gap at finite R is not the dynamically generated 4D gap.
2. In N=4 on S³ the "gap" is 1/R and **vanishes as R → ∞**; the infinite-volume CFT has no gap. So even the cleanest holographic "gap" is not the Clay gap.

A genuine holographic proof would have to show that the gap survives the decompactification and is generated by the dynamics, not by the compactification.

---

## 4. Position in 14 eras

Working periodisation of the mass-gap literature; boundaries are approximate and several programmes straddle them. Era 8 is this file.

| Era | Dates | One-line description |
|---|---|---|
| **E01** | 1954–1959 | Yang–Mills non-abelian gauge theory; quantisation and the mass problem. |
| **E02** | 1960–1972 | Quarks, current algebra, scaling; pre-QCD phenomenology. |
| **E03** | 1973–1975 | Asymptotic freedom; QCD formulated; confinement and mass gap posed. |
| **E04** | 1974–1979 | Lattice gauge theory; strong-coupling area law; 't Hooft loop; monopoles. |
| **E05** | 1979–1983 | Instantons, Savvidy, Gribov copies; large-N expansion; first lattice numerics. |
| **E06** | 1983–1990 | Rigorous RG and constructive 4D (Balaban); strong-coupling mass gap; cluster expansions. |
| **E07** | 1990–1998 | Lattice glueball spectroscopy; continuum extrapolation; dual superconductivity; Seiberg–Witten exact results. |
| **E08** | **1998–2005** | **Holography, large-N, and AdS/QCD — this file.** |
| **E09** | 2005–2010 | Soft wall, improved holographic QCD; gauge–gravity model building; lattice glueball precision. |
| **E10** | 2010–2014 | Gluon/ghost propagators, Gribov–Zwanziger, functional RG; refined lattice. |
| **E11** | 2014–2017 | Stochastic quantisation and regularity structures (Hairer); rigorous 3D; Balaban revival. |
| **E12** | 2017–2020 | Constructive 4D gauge attempts; rigorous large-N (Chatterjee); bootstrap/resurgence seeds. |
| **E13** | 2020–2023 | Bootstrap/spectral bounds; resurgence; rigorous lower bounds; volume reduction. |
| **E14** | 2023–2026 | Formalisation (Lean 4), AI-assisted proofs, Clay attack; unconditional programmes. |

---

## 5. Gaps and attack surfaces for v3

1. **No dual for pure YM (the central gap).** Every 1998–2005 duality is for a SUSY CFT (N=4) or a SUSY confining theory (KS, MN); AdS/QCD are models. *Attack:* turn the holographic dictionary into a consistency sieve — state necessary conditions a pure-YM dual must satisfy (H1–H6) and prove that the known candidates fail at least one. This is a publishable conditional/no-go result and fits the campaign's conditional-theorem strategy.

2. **NEC / IR-throat issue.** A confining holographic background needs an IR cap or throat. Hard wall = not a solution (delta-function source; the gap is `1/z₀` by hand). Soft wall = the smooth cutoff that gives `m² ∝ n` requires a dilaton stress tensor violating the null energy condition (or an uncontrolled tachyon condensate). Top-down D4 cigar = a legitimate solution with a smooth cap, but it is a 5D theory on a circle, not 4D pure YM. KS throat = a legitimate solution, but N=1 SUSY with a cascade. *Attack:* prove a no-go — any smooth, NEC-satisfying, asymptotically AdS background with a mass gap and a 4D pure-YM boundary either carries extra degrees of freedom/compactification or is not weakly curved. This would formalise why holography cannot shortcut the Clay problem. GKN 2008 is the cleanest existing analysis of the NEC constraint on the IR asymptotics.

3. **The gap is inserted by hand.** In hard/soft wall the IR scale (z₀, c) is a fit; in Sakai–Sugimoto it is 1/R. The "prediction" of the glueball spectrum is a fit to `m/√σ`, not a derivation of Δ from the UV. *Attack:* identify a UV-complete quantity that fixes the IR scale (e.g. the 't Hooft coupling at the IR scale, or the campaign's residue `f₁²`) and check whether the models can predict it rather than fit it.

4. **Large-N is not the Clay problem.** The Clay statement is fixed N (SU(3)); holography is N→∞, λ→∞. There is no controlled interpolation to N=3, finite λ, and pure YM has no SUSY protection. *Attack:* prove a gap-protection/nonrenormalization theorem, or a rigorous large-N → finite-N stability bound. Without one, a holographic gap at N=∞ says nothing about SU(3).

5. **Gravity-side gap ≠ field-theory gap.** Even with a dual, a classical supergravity gap must survive α′, g_s, 1/N, and decompactification. *Attack:* formulate the gap as a spectral lower bound for the dual Sturm–Liouville operator and test its stability under the known corrections; compare with the campaign's "uniform in volume" requirement. A gap that is only uniform at fixed compactification radius is not the Clay gap (§3.6).

6. **Wilson-loop area law is not confinement.** The Maldacena/Rey–Yee computation gives the quark potential in a specific background at strong coupling; N=4 at zero T is Coulombic, and the finite-T area law is a screening/compactification effect. *Attack:* separate "area law in a model" from "area law for pure YM"; use the campaign's gauge-invariant IR mass `m_GI = inf supp ρ_s` as the bridge and check which holographic quantity maps to it.

7. **Loop equations are exact but unsolved.** Makeenko–Migdal gives a closed equation in 4D at N=∞, but no solution; the master field is known in 2D (and on the lattice in the large-N limit). *Attack:* look for a 4D solution restricted to a truncation (e.g. the campaign's graded-strip / residue sector); test whether `f₁²` appears as an integration constant of the loop equation.

8. **Eguchi–Kawai reduction is conditional.** Naive reduction fails for D>2 because the center symmetry breaks; quenched/twisted versions restore it but change the model. *Attack:* use reduction only as a large-N consistency check, not as a solution; identify the center-symmetry order parameter with the campaign's confinement diagnostics. A reduction that requires altering the model cannot prove a gap for the original theory.

9. **Holographic renormalization is a consistency tool, not a gap proof.** Witten–Yau, Balasubramanian–Kraus, Henningson–Skenderis constrain the boundary theory (connectedness, stress tensor, Weyl anomaly) but do not produce a mass gap. *Attack:* use them to rule out candidate duals (disconnected boundary, anomaly mismatch, wrong central charges) — a cheap filter on any proposed pure-YM dual.

10. **Coupling-regime mismatch with the campaign's wall.** The holographic models live at large λ (strong coupling); the Savvidy wall (`Hess −D² + 2adF_B` indefinite; negative modes) is a weak-coupling/large-field phenomenon. No 1998–2005 construction resolves the crossover. *Attack:* this is the same missing "machine" the campaign identifies — a nonperturbative object that transports the strong-coupling gap/residue across the Savvidy crossing. Holography suggests what the machine must output (a gapped normalizable spectrum) but does not build it.

---

## 6. Verification ledger

`VERIFIED` = title/authors/venue/DOI checked against arXiv abstract page or INSPIRE-HEP record. Content is not independently re-derived. All 35 entries verified bibliographically; the two explicitly out-of-era entries are marked.

| # | Entry (short) | Checked against | Status |
|---|---|---|---|
| A1 | Maldacena, AdS/CFT | arXiv:hep-th/9711200 | VERIFIED |
| A2 | Gubser–Klebanov–Polyakov | arXiv:hep-th/9802109 | VERIFIED |
| A3 | Witten, AdS and holography | arXiv:hep-th/9802150 | VERIFIED |
| A4 | Aharony et al., review | arXiv:hep-th/9905111 | VERIFIED |
| A5 | Breitenlohner–Freedman | INSPIRE | VERIFIED |
| B6 | Witten, thermal AdS/confinement | arXiv:hep-th/9803131 | VERIFIED |
| B7 | Maldacena, Wilson loops | arXiv:hep-th/9803002 | VERIFIED |
| B8 | Rey–Yee | arXiv:hep-th/9803001 | VERIFIED |
| B9 | Brandhuber et al., finite-T loops | arXiv:hep-th/9803137 | VERIFIED |
| B10 | Polyakov, string & confinement | arXiv:hep-th/9711002 | VERIFIED |
| B11 | Witten, baryons and branes | arXiv:hep-th/9805112 | VERIFIED |
| B12 | Gross–Ooguri | arXiv:hep-th/9805129 | VERIFIED |
| C13 | Csáki–Ooguri–Oz–Terning | arXiv:hep-th/9806021 | VERIFIED |
| C14 | de Mello Koch et al. | arXiv:hep-th/9806125 | VERIFIED |
| C15 | Brower–Mathur–Tan | arXiv:hep-th/0003115 | VERIFIED |
| C16 | Klebanov–Witten, conifold | arXiv:hep-th/9807080 | VERIFIED |
| C17 | Klebanov–Tseytlin | arXiv:hep-th/0002159 | VERIFIED |
| C18 | Klebanov–Strassler | arXiv:hep-th/0007191 | VERIFIED |
| C19 | Maldacena–Nuñez | arXiv:hep-th/0008001 | VERIFIED |
| D20 | Polchinski–Strassler, hard wall | arXiv:hep-th/0109174 | VERIFIED |
| D21 | Sakai–Sugimoto I | arXiv:hep-th/0412141 | VERIFIED |
| D22 | Sakai–Sugimoto II | arXiv:hep-th/0507073 | VERIFIED |
| D23 | Erlich–Katz–Son–Stephanov | arXiv:hep-ph/0501128 | VERIFIED |
| D24 | Da Rold–Pomarol | arXiv:hep-ph/0501218 | VERIFIED |
| D25 | Karch–Katz–Son–Stephanov, soft wall | arXiv:hep-ph/0602229 | VERIFIED (2006; bridge to E09) |
| E26 | 't Hooft, planar diagrams | ADS/INSPIRE | VERIFIED |
| E27 | Makeenko–Migdal 1979 | ADS/INSPIRE | VERIFIED |
| E28 | Makeenko–Migdal 1981 | INSPIRE | VERIFIED (page end not in record) |
| E29 | Eguchi–Kawai | ADS/INSPIRE | VERIFIED |
| E30 | Bhanot–Heller–Neuberger | INSPIRE | VERIFIED |
| F31 | Henningson–Skenderis | arXiv:hep-th/9806087 | VERIFIED |
| F32 | Balasubramanian–Kraus | arXiv:hep-th/9902121 | VERIFIED |
| F33 | Witten–Yau | arXiv:hep-th/9910245 | VERIFIED |
| G34 | Morningstar–Peardon | arXiv:hep-lat/9901004 | VERIFIED |
| G35 | Gürsoy–Kiritsis–Nitti I–II | arXiv:0707.1324, 0707.1349 | VERIFIED (2008; context) |

**Counts:** 35 entries; 35 bibliographically verified; 0 unverified. Out-of-era by design: D25 (2006), G35 (2008); both flagged in place.

---

## 7. One-line handoff to v3

The era's usable residue for v3 is negative and precise: **holography supplies a dictionary and a mechanism, but not a dual of pure YM and not a gap theorem.** The three exploitable assets are (i) the H1–H6 checklist as a no-go/conditional-theorem template; (ii) the NEC/IR-throat obstruction as a hard constraint on any proposed dual; (iii) the mapping `m_GI ↔ lowest normalizable fluctuation mass` and `f₁² ↔ boundary-mode overlap` as the bridge between the campaign's gauge-invariant IR mass and any future holographic construction. None of the 1998–2005 models crosses the Savvidy wall; that remains the missing machine.
