# LIT-ERA-03 — 1974–1979: instantons, the Savvidy vacuum, and confinement pictures

> Literature base for paper v3 (4D pure SU(N) Yang–Mills mass gap).
> **Scope.** BPST instanton; 't Hooft's instanton calculus, θ-vacua, axial anomaly, U(1) problem; the Savvidy vacuum and the negative mode (Savvidy 1977, Nielsen–Olesen 1978, Leutwyler 1981); dual-superconductivity confinement ('t Hooft 1978, Mandelstam 1976); the Wilson lattice formulation and strong-coupling area law; the first glueball / mass-gap expectations.
> **Method.** Primary records checked against publisher / NASA ADS / INSPIRE-HEP / CERN Document Server / MathSciNet via WebSearch and WebFetch (INSPIRE API). 22 of 31 entries **VERIFIED**; 9 **UNVERIFIED** (standard attribution, not re-checked). No bibliographic datum is invented; corrections caught are listed in §5.3.
> **Era boundaries.** Leutwyler 1981 and 't Hooft 1981 lie just outside 1979 but are load-bearing for the Savvidy and dual-superconductor cells. They are included and flagged **post-era anchor**.

## 0. Bottom line

- This era supplies the **objects** (instantons, θ-sectors, lattice link variables, monopoles, the constant-field Hessian) and the **pictures** (dual superconductor, strong-coupling string, chromomagnetic condensate). It supplies **no theorem about the 4D continuum mass gap**.
- The only theorem-level statements are at the level of a fixed background, a model, or the lattice strong-coupling regime: BPST is an exact solution; the instanton determinant is exact in a fixed background; the Nielsen–Olesen negative eigenvalue is exact in a fixed constant background; Leutwyler's classification of stable constant fields is exact; Osterwalder–Seiler give a rigorous strong-coupling lattice gap.
- The **Savvidy "negative mode" is real and is not a gauge artifact**, but it is the curvature of a **saddle** of the one-loop effective potential, not of the true vacuum. The 1970s papers did not have the "barrier-top" reading; that is a later re-interpretation (§2.4).
- **Instanton-based mass estimates are not a mass-gap proof.** The single-instanton calculus is controlled; the instanton gas/liquid in 4D is not (large-ρ IR divergence; Witten 1979 argues the gas disappears). No 1970s paper proves Δ > 0.

---

## 1. Entries

### A. Lattice formulation and strong-coupling area law

**E01. K. G. Wilson, "Confinement of quarks," *Phys. Rev. D* **10** (1974) 2445–2459. DOI 10.1103/PhysRevD.10.2445.** `[T]` **VERIFIED**
Introduces lattice gauge theory: link variables U_μ(n) ∈ G, Wilson action S = (1/g²) Σ_p Re tr(1 − U_p), exact gauge invariance with no gauge fixing. The strong-coupling expansion is a convergent polymer expansion. Wilson argues the Wilson loop obeys an area law at small β (large bare coupling), σ > 0, with a confining flux-tube / string picture. The area law is *argued*, not proved here; the rigorous proof is E04 / Seiler's monograph. This is the origin of the lattice route to the gap: transfer matrix T = e^{−aH}, gap = −(1/a) ln(λ₂/λ₁).

**E02. R. Balian, J. M. Drouffe, C. Itzykson, "Gauge fields on a lattice. I. General outlook," *Phys. Rev. D* **10** (1974) 3376–3395. DOI 10.1103/PhysRevD.10.3376.** `[O]` **VERIFIED**
Systematic study of Wilson's lattice model with matter; algebraic structure for abelian and non-abelian groups; mean-field phase diagram; the gauge transition is first order, the matter transition second order. Establishes the model as a statistical-mechanical system. No rigorous area law. (Series: II, D11 (1975) 2098; III, D11 (1975) 2104, erratum D19 (1979) 2514.)

**E03. J. B. Kogut, L. Susskind, "Hamiltonian formulation of Wilson's lattice gauge theories," *Phys. Rev. D* **11** (1975) 395–408. DOI 10.1103/PhysRevD.11.395.** `[T]` **VERIFIED**
Canonical (Hamiltonian) formulation; gauge fields as angular variables / rigid rotators; electric-flux basis; strong-coupling expansion in the Hamiltonian language; confinement in the strong-coupling limit. Supplies the transfer-matrix / Hamiltonian side of the lattice gap.

**E04. K. Osterwalder, E. Seiler, "Gauge field theories on a lattice," *Ann. Phys.* **110** (1978) 440–471. DOI 10.1016/0003-4916(78)90039-8.** `[T]` **VERIFIED**
Rigorous Euclidean lattice gauge theory: proves the Osterwalder–Schrader axioms for the lattice measure; the strong-coupling cluster expansion gives the area law and exponential clustering with a positive mass gap at sufficiently small β, for any compact gauge group. This is the first theorem-level strong-coupling lattice mass gap. It does **not** reach the continuum limit β → ∞.

**E05. K. G. Wilson, "Quarks and strings on a lattice," in *New Phenomena in Subnuclear Physics* (Erice, 1975), ed. A. Zichichi (Plenum, New York, 1977).** `[O]` **UNVERIFIED**
Expository lecture; quarks and strings on a lattice; the strong-coupling string picture and the area law. Standard attribution; not re-checked against a record.

### B. Instantons, θ-vacua, axial anomaly, U(1)

**E06. A. A. Belavin, A. M. Polyakov, A. S. Schwarz, Yu. S. Tyupkin, "Pseudoparticle solutions of the Yang–Mills equations," *Phys. Lett. B* **59** (1975) 85–87.** `[T]` **VERIFIED**
Explicit finite-action Euclidean SU(2) solution, Q = 1, action 8π²/g², self-dual F = *F; parameters: position (4), size ρ (1), global gauge orientation (3). Topological charge Q = (g²/32π²) ∫ tr F*F. Generalizes to SU(N). The solution has **8 zero modes and no negative mode**; it is a local minimum of the action in the Q = 1 sector. This is the object that generates the instanton calculus.

**E07. G. 't Hooft, "Symmetry breaking through Bell–Jackiw anomalies," *Phys. Rev. Lett.* **37** (1976) 8–11. DOI 10.1103/PhysRevLett.37.8.** `[T]` **VERIFIED**
Anomaly-induced nonperturbative violation of U(1)_A by instantons; the axial charge is not conserved; consequences: η′ mass, baryon/lepton number violation in V−A theories with charm. The anomaly statement is theorem-level; the η′ mass connection is a dynamical inference.

**E08. G. 't Hooft, "Computation of the quantum effects due to a four-dimensional pseudoparticle," *Phys. Rev. D* **14** (1976) 3432–3450. DOI 10.1103/PhysRevD.14.3432; Erratum *Phys. Rev. D* **18** (1978) 2199. DOI 10.1103/PhysRevD.18.2199.3.** `[T]` **VERIFIED**
Instanton calculus: one-loop determinant in the single-instanton background; zero-mode Jacobian; the 't Hooft vertex; effective multi-fermion interactions; dilute instanton gas. The determinant is exact in the given background. The gas is **not** controlled in 4D. The erratum must be cited with the paper.

**E09. R. Jackiw, C. Rebbi, "Vacuum periodicity in a Yang–Mills quantum theory," *Phys. Rev. Lett.* **37** (1976) 172–175. DOI 10.1103/PhysRevLett.37.172.** `[T]` **VERIFIED**
θ-vacua: the vacuum is a superposition Σ_n e^{inθ} |n⟩ over winding sectors; θ is a new parameter; stability under large gauge transformations forces chiral non-invariance. Exact structural statement about the sector decomposition.

**E10. C. G. Callan Jr., R. F. Dashen, D. J. Gross, "The structure of the gauge theory vacuum," *Phys. Lett. B* **63** (1976) 334–340. DOI 10.1016/0370-2693(76)90277-X.** `[C]` **VERIFIED**
Instantons tunnel between topologically distinct vacua; diagonalizing the Hamiltonian gives a continuum of θ-vacua; implications for the strong interactions. Structural / heuristic, not a construction of the 4D measure.

**E11. A. M. Polyakov, "Quark confinement and topology of gauge theories," *Nucl. Phys. B* **120** (1977) 429–458. DOI 10.1016/0550-3213(77)90086-4.** `[C]` **VERIFIED**
Pseudoparticle contributions to correlation functions; argues quantum fluctuations are inessential and the one-loop approximation is effectively correct even for large pseudoparticles; reduces YM to a gas of interacting pseudoparticles. This is a conjecture; it is challenged by E13 and by the large-instanton IR divergence.

**E12. C. G. Callan Jr., R. F. Dashen, D. J. Gross, "Toward a theory of the strong interactions," *Phys. Rev. D* **17** (1978) 2717–2763. DOI 10.1103/PhysRevD.17.2717.** `[C]` **VERIFIED**
Systematic instanton-based theory of the strong interactions; vacuum structure from tunneling; the instanton size distribution and the failure of the dilute gas at large ρ. Contains the seed of the IR problem later sharpened by Witten.

**E13. E. Witten, "Instantons, the quark model, and the 1/N expansion," *Nucl. Phys. B* **149** (1979) 285–320. DOI 10.1016/0550-3213(79)90243-8.** `[NG]` **VERIFIED**
Argues that quantum corrections make the instanton gas disappear in QCD; a 2D model (CP^{N−1}, independently studied by D'Adda–Di Vecchia–Lüscher) where this is explicit. A no-go for the **naive instanton-gas picture**; not a no-go for instantons as such.

**E14. C. G. Callan Jr., R. F. Dashen, D. J. Gross, "A theory of hadronic structure," *Phys. Rev. D* **19** (1979) 1826.** `[C]` **UNVERIFIED**
Instanton-based hadronic structure; companion to E12. Cited as companion in the CERN record of E12; full record not re-checked.

### C. Confinement pictures: monopoles, vortices, dual superconductivity

**E15. A. M. Polyakov, "Particle spectrum in quantum field theory," *JETP Lett.* **20** (1974) 194–195 [*Pis'ma Zh. Eksp. Teor. Fiz.* **20** (1974) 430].** `[T]` **UNVERIFIED**
Compact 3D U(1) gauge theory (Polyakov model): monopole instantons produce a linear potential and confinement; the photon acquires a mass. Theorem-level within the model. The 3D mechanism does not transfer to 4D.

**E16. G. 't Hooft, "Magnetic monopoles in unified gauge theories," *Nucl. Phys. B* **79** (1974) 276–284.** `[T]` **UNVERIFIED**
't Hooft–Polyakov monopole in a non-abelian Higgs model; exact classical solution with quantized magnetic charge; smooth core, no Dirac string. Classical theorem.

**E17. A. M. Polyakov, "Compact gauge fields and the infrared catastrophe," *Phys. Lett. B* **59** (1975) 82–84. DOI 10.1016/0370-2693(75)90162-8.** `[C]` **VERIFIED**
Classical solutions (monopoles) govern the infrared; in compact abelian 3D the correlation length and confinement are solved. The 4D non-abelian case is left open.

**E18. S. Mandelstam, "Vortices and quark confinement in non-Abelian gauge theories," *Phys. Rep.* **23** (1976) 245–249. DOI 10.1016/0370-1573(76)90043-0.** `[V]` **VERIFIED**
Non-abelian vortices of Nielsen–Olesen type carry one unit of quantized flux, absorbed by a Dirac monopole at each end; monopole condensation confines; colour = isospin. A **picture**, not a construction; the condensate is assumed. (Proceedings of "Extended Systems in Field Theory," Paris, June 1975.)

**E19. G. 't Hooft, "On the phase transition towards permanent quark confinement," *Nucl. Phys. B* **138** (1978) 1–25. DOI 10.1016/0550-3213(78)90153-0.** `[V]` **VERIFIED**
Topological operators that modify the gauge-topological structure with local physical effects. In 2+1D, monopole operators develop a VEV when the gauge symmetry is unbroken ⇒ permanent confinement. In 3+1D, four phases: (i) Higgs, (ii) confinement, (iii) Higgs + confinement, (iv) intermediate. The dual-superconductor **phase classification**; not a construction of the 4D vacuum.

**E20. G. 't Hooft, "A property of electric and magnetic flux in non-Abelian gauge theories," *Nucl. Phys. B* **153** (1979) 141–160. DOI 10.1016/0550-3213(79)90595-9.** `[T]` **VERIFIED**
Electric and magnetic flux defined gauge-invariantly in a periodic box ('t Hooft flux / twisted boundary conditions); exact duality equation for the free energy from Euclidean invariance; for N ≤ 3 one cannot have both electric and magnetic confinement; the Georgi–Glashow IR is self-dual. Exact structural result.

**E21. V. N. Gribov, "Quantization of non-Abelian gauge theories," *Nucl. Phys. B* **139** (1978) 1–19. DOI 10.1016/0550-3213(78)90175-X.** `[O]` **VERIFIED**
Gauge-fixing ambiguity (Gribov copies); the Faddeev–Popov determinant changes sign; restricting the functional integral to the first Gribov region cancels the IR singularity and gives a linear rising potential. The copy/ambiguity is a fact; the confinement claim is heuristic. (Page range is 1–19; some citing works misprint 1–23.)

**E22. G. 't Hooft, "Topology of the gauge condition and new confinement phases in non-Abelian gauge theories," *Nucl. Phys. B* **190** [FS3] (1981) 455–478.** `[V]` **UNVERIFIED** *(post-era anchor)*
Classification of confinement phases via the 't Hooft loop; electric vs magnetic confinement; the 't Hooft algebra. The loop operator is the exact order parameter for the dual-superconductor scenario.

**E23. A. M. Polyakov, "Thermal properties of gauge fields and quark liberation," *Phys. Lett. B* **72** (1978) 477–480.** `[C]` **UNVERIFIED**
Finite-temperature gauge theory; quark liberation at high T; the deconfinement transition; the Polyakov loop. A picture / calculation, not a proof.

### D. The Savvidy vacuum and the negative mode

**E24. G. K. Savvidy, "Infrared instability of the vacuum state of gauge theories and asymptotic freedom," *Phys. Lett. B* **71** (1977) 133–134. DOI 10.1016/0370-2693(77)90759-6.** `[O]` **VERIFIED**
One-loop effective potential for a constant chromomagnetic background in pure YM. The radiative correction is negative (fixed by asymptotic freedom), so Re V_eff(B) < V_eff(0) for B ≠ 0; the perturbative vacuum is unstable; the effect is attributed to infrared singularities. The **computation** is one-loop; the existence of a chromomagnetic condensate is a **conjecture**.

**E25. I. A. Batalin, S. G. Matinyan, G. K. Savvidy, "Vacuum polarization by a constant chromomagnetic field," *Yad. Fiz.* **26** (1977) 407 [*Sov. J. Nucl. Phys.* **26** (1977) 214].** `[O]` **UNVERIFIED**
Companion one-loop calculation of the effective Lagrangian for a constant chromomagnetic field. Standard attribution; page range not re-checked.

**E26. S. G. Matinyan, G. K. Savvidy, "Vacuum polarization induced by the intense gauge field," *Nucl. Phys. B* **134** (1978) 539–545. DOI 10.1016/0550-3213(78)90463-7.** `[O]` **VERIFIED**
RG-improved effective Lagrangian in strong fields; relates strong-field asymptotics to short-distance behaviour; confirms the negative one-loop correction and the instability. **Title correction:** the paper is *"Vacuum polarization induced by the intense gauge field,"* not the commonly paraphrased *"…induces instability of gauge fields."*

**E27. N. K. Nielsen, P. Olesen, "An unstable Yang–Mills field mode," *Nucl. Phys. B* **144** (1978) 376–396.** `[T]` **VERIFIED**
Explicit diagonalization of the charged-vector fluctuation operator in a homogeneous magnetic field. The operator is (M²)_{μν} = −D² δ_{μν} − 2g F_{μν}; the spectrum is E² = k₃² + (2n+1)gB + 2gS·B, S = ±1. The lowest Landau level with S = −1 has **E² = k₃² − gB**, a negative eigenvalue **−gB** at k₃ = 0 (normalization-dependent; −2gB in conventions where the Landau spacing is 2gB). The one-loop effective action acquires an **imaginary part**; the unstable mode's vacuum polarization is a tachyon self-energy in 1+1D. Exact operator statement in the given background.

**E28. H. Leutwyler, "Constant gauge fields and their quantum fluctuations," *Nucl. Phys. B* **179** (1981) 129–170. DOI 10.1016/0550-3213(81)90252-2.** `[T]` **VERIFIED** *(post-era anchor)*
Constant gauge fields are **unstable unless (anti-)self-dual and abelian**; self-dual constant fields support infinitely many zero modes (instanton-like); one-loop corrections to the classical action are computed. The definitive one-loop classification. **Journal correction:** *Nuclear Physics B* 179, not *Physics Letters B* (a misprint in several citing works).

### E. First glueball / mass-gap expectations

**E29. H. Fritzsch, P. Minkowski, "Psi-resonances, gluons and the Zweig rule," *Nuovo Cimento A* **30** (1975) 393–420.** `[C]` **UNVERIFIED**
Gluons as constituents; gluonium states; the ψ/J resonances and the Zweig rule. A model proposal, not a mass bound.

**E30. R. L. Jaffe, K. Johnson, "Unconventional states of confined quarks and gluons," *Phys. Lett. B* **60** (1976) 201–204. DOI 10.1016/0370-2693(76)90423-8.** `[C]` **VERIFIED**
Bag-model spectrum of glueballs, Q Q̄ Q Q̄ exotics, and unconventional states; glueball masses in the range **600–1600 MeV**; the 0⁺⁺ mesons may be a four-quark nonet. Model-dependent; **no rigorous lower bound**. This is the first systematic bag-model glueball estimate.

**E31. J. B. Kogut, "An introduction to lattice gauge theory and spin systems," *Rev. Mod. Phys.* **51** (1979) 659–713.** `[O]` **UNVERIFIED**
Review of lattice gauge theory and spin systems: strong-coupling expansion, duality, phase structure. Standard attribution; not re-checked.

---

## 2. The Savvidy cluster: what is computed, what is a theorem, what is not

### 2.1 The computed object

The 1977–1981 cluster computes a **one-loop effective potential for a constant background** and the **spectrum of the quadratic fluctuation operator** about that background.

- **Background.** F_{μν} = constant in a fixed colour direction (a Cartan generator), so the field is abelian-like and the fluctuation operator is exactly solvable. This is a **variational ansatz**, not a generic field.
- **Hessian.** (M²)_{μν}^{ab} = −(D²)^{ab} δ_{μν} − 2g f^{acb} F^c_{μν}. The charged components see Landau levels; the neutral components see zero field.
- **Spectrum.** E² = k₃² + (2n+1)gB + 2gS·B, S = ±1. At n = 0 and S = −1: **E² = k₃² − gB**, i.e. a negative eigenvalue **−gB** at k₃ = 0 (normalization-dependent: −2gB in some conventions). For SU(2) the negative mode is unique; for SU(N) the negative modes sit in those root directions whose colour-magnetic coupling has the unstable sign.
- **Effective potential.** Re V_eff(B) = ½B² + c_N g²B² [ln(gB/μ²) − ½] + …, with c_N > 0 proportional to the one-loop β-function coefficient b₀ = 11N/3. For gB below the renormalization scale the bracket is negative, so the one-loop correction lowers V_eff and B = 0 is not a minimum; there is a stationary point at B_c ≠ 0 with V_eff(B_c) < 0. The imaginary part from the negative mode makes the constant-field configuration a **saddle**.

**Important distinction.** The negative mode is a property of the **constant-field background**, not of the instanton. The BPST instanton has 8 zero modes and **no** negative mode; it is a local minimum of the action in its sector. Conflating the two is a common error.

### 2.2 What is a theorem

1. **The operator spectrum in a fixed constant background is exact** (Nielsen–Olesen 1978). The negative eigenvalue is not a gauge artifact: the operator is gauge-covariant, and in a fixed background the spectrum is gauge-invariant. (The background itself is not gauge-invariant; only the gauge-invariant content of the ansatz is physical.)
2. **Leutwyler's classification is exact at one loop:** a constant gauge field is unstable unless it is (anti-)self-dual and abelian. This implies the Savvidy chromomagnetic background (non-self-dual) is unstable.
3. **The sign of the one-loop correction is fixed by the β-function**, a renormalization theorem.

### 2.3 What is *not* a theorem

- **That the true Yang–Mills vacuum is a Savvidy chromomagnetic condensate.** The one-loop potential has a stationary point at B_c ≠ 0, but it is complex; the constant-field ansatz is a saddle, not a minimum.
- **That the instability survives nonperturbatively**, or that it implies the physical vacuum is unstable / the mass gap vanishes. The negative mode is a property of a specific saddle; it says nothing about the measure at the true vacuum.
- **That the "magnetic mass" m² ~ gB is a physical pole mass.** It is a one-loop screening mass in the background.
- **That the Savvidy crossing is a genuine phase transition of the full theory.** It is an instability of the one-loop effective potential; its endpoint is not established.
- **Gauge invariance of the constant-field ansatz.** The ansatz breaks the gauge group to the stabilizer of the background; only gauge-invariant combinations are physical.

### 2.4 The later re-reading: barrier-top vs vacuum

- **1970s reading.** Savvidy proposed the chromomagnetic condensate as the true vacuum; Nielsen–Olesen exhibited the imaginary part; the perturbative vacuum was thought to be unstable, and the endpoint was taken to be a condensate.
- **1980s–2000s.** The constant-field picture became a mean-field caricature: the "Copenhagen vacuum" (Nielsen–Olesen 1980; Ambjørn–Nielsen–Olesen 1980) replaced it by a vortex/condensate structure; SVZ sum rules (1979) extracted ⟨(g²/4π²)G²⟩ ≠ 0; lattice confirmed a gluon condensate and a mass gap but **not** a constant chromomagnetic field. The QCD vacuum is a fluctuating condensate.
- **Modern effective-potential re-reading (campaign's SAVVIDY-* series).** In a model where the order parameter φ ~ gB has a double-well potential, the negative mode sits at the **barrier top** φ = 0, while the true vacua at φ = ±v have V″ = 2c > 0. The indefinite Hessian is the curvature of a **saddle**, not of the measure. This is a model-level statement; it is *not* a derivation from the 1970s papers and should be cited as reinterpretation.
- **Campaign's sharper statement.** The negative mode is rank-one in the SU(2) cell; the open question is whether that one direction can be dominated uniformly across the crossover ("Savvidy cell" = H1 / (CP)). This is a modern reduction, not a 1970s result.

### 2.5 The open cell

What survives is an **exact rank-one negative direction** plus a one-loop inequality. What is missing is a **nonperturbative, volume-uniform domination** of that direction. State this as an attack surface, not as a theorem. The 1970s literature gives the operator and the one-loop potential; it does not give the domination.

---

## 3. Position in 14 eras

| Era | Span | One line |
|---|---|---|
| 01 | 1954–65 | Yang–Mills theory formulated; non-abelian quantization obstructed (ghosts, unitarity); no confinement or gap question yet. |
| 02 | 1963–73 | Higgs mechanism and EW unification; 't Hooft (1971) renormalizability; asymptotic freedom (Gross–Wilczek, Politzer 1973) makes YM the strong-interaction candidate; confinement becomes a dynamical problem. |
| **03** | **1974–79** | **This era: lattice formulation and strong-coupling area law; instantons, θ-vacua, axial anomaly/U(1); Savvidy instability and negative mode; dual-superconductor pictures; first glueball estimates. Mechanisms, no 4D continuum theorem.** |
| 04 | 1974–82 | Rigorous lattice / constructive side: Osterwalder–Seiler axioms, Seiler's monograph, cluster expansions; area law and mass gap proved only at strong coupling. |
| 05 | 1982–90 | Balaban RG (UV stability), SVZ sum rules / gluon condensate, Gribov–Zwanziger, Creutz Monte Carlo; instanton liquid; continuum limit not under control. |
| 06 | 1985–95 | Lattice QCD matures; glueball spectrum; chiral perturbation theory; QCD sum rules; rigorous 4D construction still absent. |
| 07 | 1993–2000 | Seiberg–Witten N = 2 and Seiberg duality; AdS/CFT; MRG YM₄ with IR cutoff; confinement as dual Meissner in SUSY. |
| 08 | 1998–2005 | Holographic confinement models (Witten, Polchinski–Strassler, Sakai–Sugimoto); precision lattice glueball; domain-wall / overlap fermions. |
| 09 | 2002–10 | Continuum extrapolation, chiral fermions, Lüscher flow; FRG / QCD propagators; Gribov–Zwanziger scenarios; no proof. |
| 10 | 2008–14 | Conformal-window lattice; large-N lattice; tensor networks; resurgence; exact 2D/3D results. |
| 11 | 2013–18 | Rigorous stochastic quantization in 2D/3D; Shen–Zhu–Zhu strong-coupling gap; conformal bootstrap; 4D existence open. |
| 12 | 2018–21 | Random-surface area-law program (Cao–Sheffield); rigorous lattice YM at strong coupling; gradient-flow rigor; mass gap still the wall. |
| 13 | 2021–24 | CPS/CNS area law; stochastic 3D YM–Higgs; rigorous 2D YM; 4D mass gap open. |
| 14 | 2024–26 | Lean / formalization wave; OpenAI 4D existence claim and Gross committee; hexagonmath 16-input reduction to MG; CPS/CNS λ < 1/24; existence conditional, gap open. |

---

## 4. Gaps and attack surfaces for v3

- **Savvidy cell — separate the three layers.** (i) exact operator spectrum in a fixed background (theorem); (ii) one-loop V_eff (computation); (iii) stability of the true vacuum (open). Do **not** cite Savvidy 1977 / Nielsen–Olesen 1978 as evidence that the QCD vacuum is unstable. Cite them for the operator and the one-loop potential, and mark the vacuum reading as a conjecture.
- **Negative-mode counting for SU(N).** The 1970s literature is SU(2)-centric. v3 needs the exact root-space count (one per unstable root direction) and its gauge-invariant status. The campaign's "rank-one" reading is an SU(2)-cell statement, not a general theorem.
- **Barrier-top vs vacuum.** The double-well re-reading is a model-level statement about the effective potential; it is not in the 1970s papers. Cite it as reinterpretation, and keep the "true vacuum stability" as an open cell.
- **Instanton-based mass estimates.** 't Hooft's determinant and zero modes are exact in the single-instanton background, but the 4D gas/liquid is uncontrolled (large-ρ IR divergence; Witten 1979 argues the gas disappears). No 1970s paper proves Δ > 0. v3 must not present the instanton liquid as a mass-gap derivation.
- **U(1)_A vs mass gap.** The anomaly resolution is a theorem about the axial current; it does not give a pure-YM mass gap. Keep the η′ / U(1) problem separate from the 0⁺⁺ glueball gap.
- **Dual superconductivity is a picture.** 't Hooft 1978 / Mandelstam 1976 give a phase classification and a topological algebra, not an order-parameter theorem. The monopole condensate is an assumption. v3 should state the missing order parameter explicitly.
- **Area law vs mass gap.** The Wilson-loop area law and the transfer-matrix gap are different statements; the 1970s literature often conflates them. This is exactly the campaign's G1/G2 split, and v3 should state it against the historical conflation.
- **Strong-coupling regime.** Wilson 1974 / Osterwalder–Seiler 1978 give a genuine gap only at small β. The continuum limit β → ∞ is untouched. Quote the exact regime and what fails at large β.
- **Glueball masses are model numbers.** Bag-model values (600–1600 MeV) are model-dependent; there is no rigorous lower bound in this era. Do not use them as evidence for Δ > 0.
- **θ-vacua are a sector label, not a gap tool.** θ labels the vacuum family; its physical effect in pure YM is a vacuum-energy question, not a mass-gap mechanism. Keep it out of the gap argument.

---

## 5. Verification ledger

### 5.1 Method

- **WebSearch** against NASA ADS, INSPIRE-HEP, CERN Document Server, publisher pages (APS, Elsevier/ScienceDirect), MathSciNet, and Google Scholar records.
- **WebFetch** against the INSPIRE-HEP REST API for one record (E08), which returned the journal, volume, pages, DOI, and erratum.
- DOIs are listed only where returned by a record, or where given by the deterministic APS scheme 10.1103/<journal>.<volume>.<page> (E01). No DOI is guessed for an UNVERIFIED entry.
- Search budget for this turn was exhausted after 22 entries; the remaining 9 were not re-checked and are marked **UNVERIFIED**. They are standard attributions, not inventions; a follow-up pass can upgrade them.

### 5.2 Ledger

| # | Entry (short) | Check performed | Status |
|---|---|---|---|
| E01 | Wilson 1974 | ADS/publisher: title, PRD 10, 2445–2459 | **VERIFIED** |
| E02 | Balian–Drouffe–Itzykson 1974 | ADS/APS: title, PRD 10, 3376–3395 | **VERIFIED** |
| E03 | Kogut–Susskind 1975 | APS: title, PRD 11, 395–408 | **VERIFIED** |
| E04 | Osterwalder–Seiler 1978 | Scopus/MathSciNet: Ann. Phys. 110, 440–471 | **VERIFIED** |
| E05 | Wilson 1975 (Erice) | not checked | **UNVERIFIED** |
| E06 | BPST 1975 | ADS/CERN: PLB 59, 85–87 | **VERIFIED** |
| E07 | 't Hooft PRL 37 1976 | ADS/APS: PRL 37, 8–11 | **VERIFIED** |
| E08 | 't Hooft PRD 14 1976 | INSPIRE API: PRD 14, 3432–3450 + erratum PRD 18, 2199 | **VERIFIED** |
| E09 | Jackiw–Rebbi 1976 | ADS/APS: PRL 37, 172–175 | **VERIFIED** |
| E10 | CDG PLB 63 1976 | CERN: PLB 63, 334–340 | **VERIFIED** |
| E11 | Polyakov NPB 120 1977 | ADS/Scopus: NPB 120, 429–458 | **VERIFIED** |
| E12 | CDG PRD 17 1978 | CERN/APS: PRD 17, 2717–2763 | **VERIFIED** |
| E13 | Witten NPB 149 1979 | INSPIRE/Semantic Scholar: NPB 149, 285–320 | **VERIFIED** |
| E14 | CDG PRD 19 1979 | indirect (CERN record of E12); full record not checked | **UNVERIFIED** |
| E15 | Polyakov 1974 | not checked | **UNVERIFIED** |
| E16 | 't Hooft NPB 79 1974 | not checked | **UNVERIFIED** |
| E17 | Polyakov PLB 59 1975 | ADS/OSTI: PLB 59, 82–84 | **VERIFIED** |
| E18 | Mandelstam Phys. Rep. 23 1976 | ScienceDirect/CERN: 23, 245–249 | **VERIFIED** |
| E19 | 't Hooft NPB 138 1978 | INSPIRE/ADS: NPB 138, 1–25 | **VERIFIED** |
| E20 | 't Hooft NPB 153 1979 | ADS: NPB 153, 141–160 | **VERIFIED** |
| E21 | Gribov NPB 139 1978 | ADS/MathSciNet: NPB 139, 1–19 | **VERIFIED** |
| E22 | 't Hooft NPB 190 1981 | not checked | **UNVERIFIED** |
| E23 | Polyakov PLB 72 1978 | not checked | **UNVERIFIED** |
| E24 | Savvidy PLB 71 1977 | ADS/SLUB: PLB 71, 133–134 | **VERIFIED** |
| E25 | Batalin–Matinyan–Savvidy 1977 | not checked | **UNVERIFIED** |
| E26 | Matinyan–Savvidy NPB 134 1978 | ADS/Semantic Scholar: NPB 134, 539–545; title corrected | **VERIFIED** |
| E27 | Nielsen–Olesen NPB 144 1978 | ADS: NPB 144, 376–396 | **VERIFIED** |
| E28 | Leutwyler NPB 179 1981 | INSPIRE/ScienceDirect: NPB 179, 129–170; journal corrected | **VERIFIED** |
| E29 | Fritzsch–Minkowski 1975 | not checked | **UNVERIFIED** |
| E30 | Jaffe–Johnson PLB 60 1976 | ADS/OSTI: PLB 60, 201–204 | **VERIFIED** |
| E31 | Kogut RMP 51 1979 | not checked | **UNVERIFIED** |

**Counts:** 31 entries; 22 **VERIFIED**, 9 **UNVERIFIED**.

### 5.3 Corrections caught

1. **Matinyan–Savvidy title.** The paper is *"Vacuum polarization induced by the intense gauge field,"* Nucl. Phys. B 134 (1978) 539–545 — not the commonly paraphrased *"Vacuum polarization induces instability of gauge fields."* Correct before v3 cites it.
2. **Leutwyler journal.** *Nuclear Physics B* 179 (1981) 129–170, **not** *Physics Letters B* 179. Several citing works misprint this.
3. **Gribov page range.** Nucl. Phys. B 139 (1978) **1–19**, not 1–23. Some citing works misprint.
4. **'t Hooft 1976 erratum.** Phys. Rev. D 14 (1976) 3432 must be cited together with the erratum Phys. Rev. D 18 (1978) 2199 (DOI 10.1103/PhysRevD.18.2199.3). The erratum corrects a factor in the determinant.
5. **Nielsen–Olesen page range.** Nucl. Phys. B 144 (1978) **376–396**, not merely 376. The negative-mode analysis is in this paper, not in the 1973 Nielsen–Olesen vortex paper (Nucl. Phys. B 61 (1973) 45), which is a different object.
6. **BPST vs Savvidy.** The BPST instanton has zero modes and no negative mode; the Savvidy constant-field background has the negative mode. Do not transfer the "negative mode" language to the instanton.

### 5.4 What is deliberately absent

- No claim that any 1970s paper proves the 4D mass gap. The era supplies mechanisms and exact local computations only.
- No claim that the Savvidy condensate is the true vacuum. The one-loop potential is a saddle; the vacuum reading is a conjecture.
- No claim that the instanton gas gives a controlled mass estimate in 4D. The single-instanton calculus is controlled; the gas is not.
- No invented page ranges, DOIs, or author orders. Where a datum was not re-checked, the entry is marked **UNVERIFIED**.
