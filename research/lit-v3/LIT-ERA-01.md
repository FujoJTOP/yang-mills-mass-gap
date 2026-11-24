# LIT-ERA-01 — 1954–1965: The Origin of the Gauge–Mass Problem

**Campaign:** 4D pure SU(N) Yang–Mills mass gap (Clay). **Scope:** 1954–1965, plus the minimum mathematical and field-theoretic prehistory needed to locate the mass question. **Status of the mass gap:** OPEN. Nothing in this era proves it; several things in this era define it.

**Campaign scheme:** `[T]` theorem; `[C]` conditional theorem; `[NG]` no-go for a named class; `[O]` exact open/missing step; `[V]` numerical/experimental evidence.

**Method note.** Entries were checked against APS / Elsevier / Springer / INSPIRE / zbMATH / ADS records where the shared search budget allowed. The two mathematics-prehistory entries (E01–E02) could not be checked before the turn's shared WebSearch budget was exhausted; they are marked `UNVERIFIED`, and no page range or DOI is asserted for them. For all other entries the bibliographic data below are from the checked record. Content claims are sourced in the annotation; where the primary text was not fetched, that is stated.

---

## 0. Executive verdict

1. **What Yang–Mills 1954 established.** A consistent *classical* non-Abelian gauge theory: the field strength \(F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu+[A_\mu,A_\nu]\), the gauge-covariant derivative, and the fact that a local mass term \(m^2 A_\mu^a A^{a\mu}\) is not gauge invariant. It did **not** establish a quantum theory, renormalizability, or a mass gap. The quanta are massless at the classical level.

2. **The precise origin of “gauge invariance forbids mass.”** It is the classical algebraic statement that \(A_\mu\) is a connection, not a tensor: under \(A_\mu\to gA_\mu g^{-1}+\frac{i}{g}(\partial_\mu g)g^{-1}\), the bilinear \(A_\mu A^\mu\) is not invariant, so no gauge-invariant local term of dimension two built from \(A\) alone exists. The further claim that the *physical spectrum* must therefore be massless is a **heuristic**, not a theorem. It assumes the mass must arise from a local term in the Lagrangian containing only the gauge field.

3. **The heuristic was already false in the era.** Stueckelberg (1938, Abelian, with a scalar) and Schwinger's exact 2D model (1962) show that gauge invariance and a massive vector can coexist. The 1964 Higgs / Englert–Brout / Guralnik–Hagen–Kibble mechanism is the non-Abelian version at the classical/perturbative level.

4. **The genuine theorem of the era is Goldstone's** (1961; proved by Goldstone–Salam–Weinberg 1962): spontaneous breaking of a continuous *global* symmetry forces a massless boson. This is a `[NG]` for the class “global symmetry + SSB + Lorentz invariance + positivity.” It is not a no-go for gauge theories, where the hypotheses fail.

5. **What the era did not do.** It did not construct 4D pure Yang–Mills as a quantum theory, did not prove renormalizability (1971–72), and did not prove a mass gap. For pure YM the 1964 scalar mechanism is unavailable, so the gap must be dynamical. That is the modern Clay problem.

6. **“Mass gap” language.** The phrase is absent from Yang–Mills 1954. The earliest documented definition found here is Streater–Wightman (1964): a theory has a mass gap if the mass spectrum has no mass below a positive threshold \(M\). The definition is kinematic and proves nothing.

---

## 1. Entries

### A. Mathematical prehistory: connection geometry

**E01. C. Ehresmann**, “Les connexions infinitésimales dans un espace fibré différentiable,” *Colloque de Topologie (Espaces Fibrés)*, Bruxelles, 1950. `UNVERIFIED` · `[T]` — Standard citation; exact page range not verified in this session, so omitted. Mass relevance only: a gauge field is a connection on a principal bundle, and the curvature is the only gauge-covariant algebraic object. This is the geometric form of the statement that a local \(m^2A^2\) is not allowed. It says nothing about the quantum spectrum.

**E02. N. Steenrod**, *The Topology of Fibre Bundles*, Princeton University Press, 1951. `UNVERIFIED` · `[T]` — Standard book citation; no DOI; not checked in this session. Same geometric point: the relevant object is the connection/curvature, not a vector field with a mass term. No mass-gap content.

### B. The Abelian mass loophole and the spectral framework

**E03. A. Proca**, “Sur la théorie ondulatoire des électrons positifs et négatifs,” *Journal de Physique et le Radium* (7) **7**, 347–353 (1936). DOI: `10.1051/jphysrad:0193600708034700`. `VERIFIED` · `[T]` — Classical massive spin-1 theory: \((\partial_\mu\partial^\mu+m^2)A_\nu=0\) with \(\partial^\nu A_\nu=0\). Establishes that a massive vector is a consistent *classical* field; the mass term is not gauge invariant. This is the theory the gauge principle had to embed or supersede.

**E04. E. C. G. Stueckelberg**, “Die Wechselwirkungskräfte in der Elektrodynamik und in der Feldtheorie der Kernkräfte (Teil I),” *Helvetica Physica Acta* **11**, 225–244 (1938); Parts II–III ibid. 299–328. No DOI. `VERIFIED` · `[T]` — Shows the Proca mass can be written gauge-invariantly by adding a real scalar \(\sigma\):
\[
\mathcal L=-\tfrac14F_{\mu\nu}F^{\mu\nu}-\tfrac12(mA_\mu+\partial_\mu\sigma)^2,
\]
invariant under \(A_\mu\to A_\mu+\partial_\mu\varepsilon,\ \sigma\to\sigma-m\varepsilon\). This is an explicit counterexample to the blanket claim “gauge invariance forbids a massive vector” in the Abelian case. It is a classical equivalence; it does not give a non-Abelian mass mechanism, and it does not prove a mass gap.

**E05. G. Källén**, “On the Definition of the Renormalization Constants in Quantum Electrodynamics,” *Helvetica Physica Acta* **25**, 417–434 (1952). DOI: `10.5169/seals-112316`. `VERIFIED` · `[T]` — Spectral representation
\[
\Delta(p)=\int d\mu^2\,\frac{\rho(\mu^2)}{p^2-\mu^2+i\epsilon}.
\]
The mass spectrum is the support of \(\rho\); a mass gap is a statement about that support. The representation is a theorem under standard assumptions. It is not a tool that computes the lower edge of the support.

**E06. W. Pauli**, unpublished 1953 non-Abelian construction, in *Wissenschaftlicher Briefwechsel mit Bohr, Einstein, Heisenberg u.a.*, Vol. IV, Part II: 1953–1954, ed. K. von Meyenn (Springer, 1999), ISBN 978-3-540-64312-8. `VERIFIED (volume); PARTIAL (letters)` · `[O]` — Pauli constructed the non-Abelian field strength in 1953 letters and did not publish; the standard historical account attributes this to the gauge-boson mass problem (secondary source: N. Straumann, “Wolfgang Pauli and Modern Physics,” arXiv:0810.2213). The specific letter dates/recipients were not checked here. Pauli's objection was a physical-consistency worry, not a proof. The often-repeated Princeton-seminar anecdote (“What is the mass of these particles?”) is reported in secondary histories; the primary record was not verified.

**E07. H. Lehmann**, “Über Eigenschaften von Ausbreitungsfunktionen und Renormierungskonstanten quantisierter Felder,” *Il Nuovo Cimento* **11**, 342–357 (1954). DOI: `10.1007/BF02783624`. `VERIFIED` · `[T]` — Independent derivation of the spectral representation (Källén–Lehmann). Establishes the analytic structure of the two-point function; no dynamical gap.

### C. 1954: the non-Abelian construction and the mass objection

**E08. C. N. Yang and R. L. Mills**, “Conservation of Isotopic Spin and Isotopic Gauge Invariance,” *Physical Review* **96**, 191–195 (1954). DOI: `10.1103/PhysRev.96.191`. `VERIFIED` · `[T]` classical; `[O]` quantum — Defines \(F_{\mu\nu}\), the gauge-covariant derivative, and the classical Lagrangian. Notes that \(m^2A^2\) is not invariant, so the quanta are massless. The paper's own text was not fetched here; the mass-term statement is confirmed by multiple secondary sources. The paper does not quantize the theory, does not prove renormalizability, and does not prove a gap. The massless charged vector bosons conflict with the short-range isospin force; this is the origin of the “gauge invariance forbids mass” worry in the non-Abelian setting.

**E09. R. Shaw**, “The Problem of Particle Types and Other Contributions to the Theory of Elementary Particles,” PhD thesis, University of Cambridge, 1955 (unpublished). No DOI. `VERIFIED` (title/institution/advisors) · `[O]` — Independent 1955 construction of non-Abelian gauge fields; same mass problem. Verified via the Mathematics Genealogy Project and historical sources; the thesis is unpublished. Adds no mass mechanism; establishes shared priority.

### D. 1955–1957: gauge principle, axiomatic fields, condensed-matter gap

**E10. R. Utiyama**, “Invariant Theoretical Interpretation of Interaction,” *Physical Review* **101**, 1597–1607 (1956). DOI not checked. `VERIFIED` (vol./pages) · `[T]` classical — Generalizes the gauge principle: promoting a global group to a local one requires a gauge field; the field strength is the curvature. It makes precise why a local mass term is excluded and why the only algebraic invariants are curvature polynomials and their covariant derivatives. It does not address the quantum spectrum.

**E11. A. S. Wightman**, “Quantum Field Theory in Terms of Vacuum Expectation Values,” *Physical Review* **101**, 860–866 (1956). DOI: `10.1103/PhysRev.101.860`. `VERIFIED` · `[T]` — Axiomatic formulation via vacuum expectation values; reconstruction theorem. The spectrum condition makes the mass spectrum a measure. This is the framework in which “mass gap” later becomes a precise statement; it does not prove one.

**E12. J. Bardeen, L. N. Cooper, J. R. Schrieffer**, “Theory of Superconductivity,” *Physical Review* **108**, 1175–1204 (1957). DOI: `10.1103/PhysRev.108.1175`. `VERIFIED` · `[T]` condensed matter — Establishes an energy gap for single-particle excitations, \(2\Delta\approx3.5kT_c\). Not a QFT theorem. Its importance for the gauge-mass problem is the analogy later exploited by Anderson and by Nambu–Jona-Lasinio: a gap can arise dynamically from a condensate. It does not address relativistic gauge invariance.

### E. 1960–1962: vector mesons, Goldstone, dynamical breaking

**E13. J. J. Sakurai**, “Theory of Strong Interactions,” *Annals of Physics* **11**, 1–48 (1960). DOI: `10.1016/0003-4916(60)90126-3`. `VERIFIED` · `[O]` — Proposes massive vector mesons coupled to conserved isospin/hypercharge/baryon currents, with universal coupling (vector-meson dominance). The \(\rho,\omega,\phi,K^*\) were found 1961–63. This is a phenomenological success, not a gauge theory with a mass gap; the mass is put in by hand.

**E14. J. Goldstone**, “Field theories with « Superconductor » solutions,” *Il Nuovo Cimento* **19**, 154–164 (1961). DOI: `10.1007/BF02812722`. `VERIFIED` · `[T]`/`[NG]` — States that a theory with a spontaneously broken continuous symmetry contains a massless boson. The 1961 paper is the statement plus model evidence; the proof is Goldstone–Salam–Weinberg 1962. This is the era's genuine no-go for global symmetries. It sharpens the mass problem: naive SSB gives massless scalars, not a massive vector.

**E15. Y. Nambu and G. Jona-Lasinio**, “Dynamical Model of Elementary Particles Based on an Analogy with Superconductivity. I,” *Physical Review* **122**, 345–358 (1961). DOI: `10.1103/PhysRev.122.345`. `VERIFIED` · `[C]` model — Dynamical chiral symmetry breaking in a four-fermion model; fermions acquire mass through a condensate. Not a gauge theory, and not renormalizable in 4D in the naive sense. Important as the first serious “mass from dynamics” template.

**E16. S. L. Glashow and M. Gell-Mann**, “Gauge Theories of Vector Particles,” *Annals of Physics* **15**, 437–460 (1961). DOI: `10.1016/0003-4916(61)90193-2`. `VERIFIED` · `[O]` — Classifies gauge theories of vector particles for simple Lie algebras; notes that mass terms can be added only by breaking gauge invariance for spacetime-dependent transformations. It states the obstruction, not a no-go theorem for the spectrum. Also discusses intermediate vector bosons and SU(3) ideas.

**E17. S. L. Glashow**, “Partial-symmetries of weak interactions,” *Nuclear Physics* **22**, 579–588 (1961). DOI: `10.1016/0029-5582(61)90469-2`. `VERIFIED` · `[O]`/`[C]` model — SU(2)×U(1) electroweak model with ad hoc vector masses. The masses are not generated by a mechanism, and the model was not renormalizable as written. It shows the phenomenological demand for massive gauge bosons.

**E18. J. Schwinger**, “Gauge Invariance and Mass,” *Physical Review* **125**, 397–398 (1962). DOI: `10.1103/PhysRev.125.397`. `VERIFIED` · `[O]` — Short paper arguing that gauge invariance can imply a nonzero mass; the argument is heuristic. It sets up the exact 2D treatment in the companion paper.

**E19. J. Schwinger**, “Gauge Invariance and Mass. II,” *Physical Review* **128**, 2425–2429 (1962). DOI: `10.1103/PhysRev.128.2425`. `VERIFIED` · `[T]` model — Exact solution of 1+1-dimensional QED with a massless fermion: the photon acquires a mass \(e/\sqrt{\pi}\) while gauge invariance is maintained. This is a theorem for the model and an explicit counterexample to the blanket “gauge invariance forbids mass.” It does not transfer to 4D non-Abelian theory.

**E20. J. Goldstone, A. Salam, S. Weinberg**, “Broken Symmetries,” *Physical Review* **127**, 965–970 (1962). DOI: `10.1103/PhysRev.127.965`. `VERIFIED` · `[T]`/`[NG]` — Proof of Goldstone's conjecture under standard Wightman assumptions. The hypotheses include a genuine global symmetry and positivity; a gauge theory evades them. This is the precise no-go that the Higgs mechanism had to circumvent.

### F. 1963: quantization, the plasmon analogy, the ghost problem

**E21. R. P. Feynman**, “Quantum theory of gravitation,” *Acta Physica Polonica* **24**, 697–722 (1963). No DOI. `VERIFIED` · `[O]` technique — Lectures on quantum gravitation containing the non-Abelian gauge-field case; Feynman found that naive Feynman rules violate unitarity and that ghost loops are needed. Unpublished in full until *Feynman Lectures on Gravitation* (1995). It establishes the quantization difficulty, not a mass result.

**E22. P. W. Anderson**, “Plasmons, Gauge Invariance, and Mass,” *Physical Review* **130**, 439–442 (1963). DOI: `10.1103/PhysRev.130.439`. `VERIFIED` · `[C]`/heuristic — Plasmon analogy: in a medium, the photon acquires an effective mass while gauge invariance is maintained; the Goldstone mode is the longitudinal plasmon. It is the condensed-matter precursor of the Higgs mechanism, not a QFT proof.

### G. 1964: the scalar mechanism and the electroweak models

**E23. F. Englert and R. Brout**, “Broken Symmetry and the Mass of Gauge Vector Mesons,” *Physical Review Letters* **13**, 321–323 (1964). DOI: `10.1103/PhysRevLett.13.321`. `VERIFIED` · `[C]` — Classical/perturbative mechanism: a scalar VEV gives mass to a gauge vector; the would-be Goldstone mode is absorbed. The paper is a short letter; quantum renormalizability is not proved.

**E24. P. W. Higgs**, “Broken Symmetries and the Masses of Gauge Bosons,” *Physical Review Letters* **13**, 508–509 (1964). DOI: `10.1103/PhysRevLett.13.508`. `VERIFIED` · `[C]` — Same mechanism, with the explicit observation that the massless scalar is eaten. The detailed perturbative treatment is the 1966 paper (E30).

**E25. G. S. Guralnik, C. R. Hagen, T. W. B. Kibble**, “Global Conservation Laws and Massless Particles,” *Physical Review Letters* **13**, 585–587 (1964). DOI: `10.1103/PhysRevLett.13.585`. `VERIFIED` · `[C]` — Independent derivation; emphasizes that no massless particle survives in the physical spectrum when the symmetry is gauged.

**E26. B. S. DeWitt**, “Theory of Radiative Corrections for Non-Abelian Gauge Fields,” *Physical Review Letters* **12**, 742–746 (1964). DOI: `10.1103/PhysRevLett.12.742`. `VERIFIED` · `[O]` technique — Background-field/ghost formalism for non-Abelian gauge fields; one-loop radiative corrections; unitarity. It makes the quantization of Yang–Mills systematic at one loop but does not prove renormalizability or a gap.

**E27. A. Salam and J. C. Ward**, “Electromagnetic and weak interactions,” *Physics Letters* **13**, 168–171 (1964). DOI: `10.1016/0031-9163(64)90711-5`. `VERIFIED` · `[C]`/`[O]` — SU(2)×U(1) electroweak model; masses again put in by hand at this stage. The gauge structure is right; the mass-generation and renormalizability are not yet established.

### H. 1964: “mass gap” enters the language

**E28. R. F. Streater and A. S. Wightman**, *PCT, Spin and Statistics, and All That*, W. A. Benjamin (1964); Princeton Landmarks reprint (2000). No DOI. `VERIFIED` (book + snippet) · `[C]` — Textbook axiomatic QFT. Defines a mass gap: the mass spectrum has no mass below a positive threshold \(M\). Proves conditional clustering statements: a gap implies exponential decay of correlations. This is the first documented “mass gap” definition found in the era; it is kinematic, not a proof.

**E29. R. P. Feynman**, *Feynman Lectures on Gravitation*, ed. B. Hatfield, Addison-Wesley (1995), ISBN 0-201-62734-5; based on the Caltech course of 1962–63. `VERIFIED` · `[O]` — The 1962–63 course notes contain Feynman's non-Abelian quantization work and the ghost problem; same status as E21. No mass-gap result.

### Boundary (immediately after the strict 1954–65 window)

**E30. P. W. Higgs**, “Spontaneous Symmetry Breakdown without Massless Bosons,” *Physical Review* **145**, 1156–1163 (1966). DOI: `10.1103/PhysRev.145.1156`. `VERIFIED` · `[C]` — Detailed perturbative treatment of the U(1) model; the Goldstone boson becomes the longitudinal vector mode. Outside the strict window but the quantitative completion of the 1964 letters.

**E31. L. D. Faddeev and V. N. Popov**, “Feynman diagrams for the Yang-Mills field,” *Physics Letters B* **25**, 29–30 (1967). DOI: `10.1016/0370-2693(67)90067-6`. `VERIFIED` · `[T]` perturbative — Correct Feynman rules for Yang–Mills with gauge fixing and ghosts; the starting point for modern quantization. Outside the strict window; belongs to era 02.

---

## 2. Cluster verdicts: established vs hoped

| Cluster | Established | Only hoped | Residual gap |
|---|---|---|---|
| **A. Geometry (E01–E02)** | Gauge field = connection; curvature = field strength; no local \(m^2A^2\) invariant. | That bundle geometry would guide quantization and the mass problem. | No link from bundle geometry to a spectral gap. |
| **B. Abelian loophole + spectral tools (E03–E07)** | Proca massive vector; Stueckelberg gauge-invariant mass; Källén–Lehmann representation; Wightman framework. | That the Abelian Stueckelberg trick would generalize to non-Abelian gauge fields. | Non-Abelian massive theory without scalars; a computable spectral lower bound. |
| **C. 1954 construction (E08–E09)** | Classical non-Abelian gauge theory; mass term forbidden; massless quanta. | Quantization, renormalizability, a later mass mechanism. | Quantized 4D pure YM; dynamical gap. |
| **D. 1955–57 (E10–E12)** | Gauge principle for arbitrary groups (Utiyama); axiomatic QFT; BCS energy gap. | That these tools would solve the relativistic mass problem. | No relativistic gap theorem; no gauge-theory mass mechanism. |
| **E. 1960–62 (E13–E20)** | Goldstone theorem (genuine no-go for global SSB); vector mesons; NJL dynamical mass; classification of gauge theories with mass terms. | That SSB would give gauge masses without unwanted Goldstones; that vector mesons were gauge bosons. | Gauge version of Goldstone's theorem; nonperturbative mass generation. |
| **F. 1962–63 (E18–E22)** | Schwinger model (exact 2D mass without breaking gauge invariance); plasmon analogy; ghost problem in non-Abelian quantization. | That the 2D/condensed-matter mechanisms transfer to 4D non-Abelian theory. | No 4D transfer; no pure-YM mass. |
| **G. 1964 (E23–E27, E30–E31)** | Classical/perturbative scalar mechanism; electroweak model building; systematic one-loop quantization. | Renormalizability and unitarity in 4D; a nonperturbative construction. | 4D renormalizability (1971–72); pure-YM gap. |
| **H. 1964–65 (E28)** | Precise “mass gap” definition; conditional exponential clustering. | That the definition would help prove gaps in interacting 4D theories. | No interacting 4D example with a proven gap. |

### 2.1 The precise origin of the worry, and its logical status

The worry has two layers, and they must not be conflated.

**Layer 1 (algebraic, classical, true).** A gauge potential is a connection. Under a gauge transformation it shifts by an inhomogeneous term, so \(A_\mu^a A^{a\mu}\) is not invariant. The only local gauge-invariant objects built algebraically from \(A\) are polynomials in the curvature \(F\) and its covariant derivatives. Therefore a local mass term for the gauge field is excluded from a gauge-invariant Lagrangian written in terms of \(A\) alone. This is a correct classical statement. It is a `[NG]` only for the very narrow class “local Lagrangian, fields \(A\) only, mass term of dimension two.” It is not a statement about the physical spectrum.

**Layer 2 (spectral, heuristic, not a theorem).** The inference “therefore the physical gauge bosons are massless” assumes that any mass must come from a local term in the Lagrangian with only \(A\). That assumption is false in general:

- Stueckelberg (1938) writes a gauge-invariant Abelian massive vector by adding a scalar; the physical spectrum has a massive vector and a scalar.
- Schwinger (1962) solves a 2D model exactly in which the photon acquires mass while gauge invariance is maintained.
- The 1964 scalar mechanism gives a non-Abelian vector a mass at the classical/perturbative level by eating a would-be Goldstone mode.
- A massless elementary gauge field can be confined; the physical spectrum can be gapped even if the elementary field is not. This is the modern pure-YM scenario.

So “gauge invariance forbids mass” was never a theorem about the spectrum. It was a correct algebraic observation plus an unjustified extrapolation. The genuine theorem of the era is Goldstone's, and it is a theorem about *global* symmetries: if a continuous global symmetry is spontaneously broken in a Lorentz-invariant theory with positivity, a massless boson exists. A gauge theory evades the hypotheses (the symmetry is local; the would-be Goldstone mode is not a physical asymptotic state), which is exactly what the 1964 mechanism exploits.

**Why the worry was nonetheless rational in 1954–65.** There was no known non-Abelian mechanism for mass without a scalar, no proof that pure YM confines, and no nonperturbative construction. The only known way to give a vector a mass while keeping a local Lagrangian was the Proca/Stueckelberg route, which for non-Abelian fields was not obviously consistent or renormalizable. The 1964 scalar mechanism solved the problem only by adding a scalar; for pure YM it is unavailable. Hence the era's honest position: the local mass term is forbidden, the spectrum is not thereby proved massless, and the dynamical question is open.

---

## 3. Position in 14 eras

| Era | How the 1954–65 theme recurs or dies |
|---|---|
| **01 1954–65** | The local mass term is forbidden; gauge invariance vs mass becomes the central obstacle; first escapes (Stueckelberg, Schwinger model, Higgs). |
| **02 1963–73** | Higgs plus 't Hooft–Veltman make SSB gauge theories renormalizable; the 1954 worry survives only for pure YM, where the gap must be dynamical. |
| **03 1974–79** | Asymptotic freedom, instantons, lattice gauge theory: the gap becomes a nonperturbative statistical-mechanics problem; “mass term” is replaced by “dynamical mass generation.” |
| **04 1974–82** | Monopoles, dual superconductivity, confinement scenarios: the gap is tied to a condensate, not a local mass term; the 1954 algebraic obstruction is irrelevant. |
| **05 1982–90** | Lattice Monte Carlo produces numerical evidence for the gap and string tension; the question becomes computational, not a proof. |
| **06 1985–95** | Renormalization-group and continuum-limit work: the gap is dimensional transmutation; the 1954 classical scale invariance is the precursor. |
| **07 1993–2000** | Seiberg–Witten exact gaps in supersymmetric toy models; pure YM remains open; the 1954 question is solved only in deformed theories. |
| **08 1998–2005** | Clay Millennium formulation (Jaffe–Witten): the 1954 Lagrangian becomes an axiomatic existence-plus-gap problem. |
| **09 2002–10** | Constructive QFT and rigorous RG; the 1954 quantization difficulties recur as measure-theory problems. |
| **10 2008–14** | Holographic/AdS-QCD models and lattice: the gap is modeled, not proved. |
| **11 2013–18** | Rigorous lattice/Osterwalder–Schrader methods, gradient flow; the 1954 transfer-matrix idea matures. |
| **12 2018–21** | Bootstrap, positivity, and machine learning; the spectral (Källén–Lehmann) meaning of the gap is exploited. |
| **13 2021–24** | Rigorous bounds and constructive attempts; “gauge invariance forbids mass” is fully replaced by “prove a spectral gap.” |
| **14 2024–26** | v3 unconditional program; the 1954 question is still open, but the tools are nonperturbative. |

---

## 4. Gaps and attack surfaces for v3

- **Missing theorem: nonperturbative construction of 4D pure YM.** The 1954 Lagrangian is classical. There is no measure on the space of connections satisfying the Osterwalder–Schrader axioms with a uniform gap. Attack surface: constructive RG with a uniform infrared bound (the Balaban program and its successors), not a perturbative expansion around \(A=0\).
- **Missing theorem: separation of “no local mass term” from “no spectral gap.”** The algebraic no-go concerns the Lagrangian; the physical statement concerns the mass operator. Attack surface: define a gauge-invariant mass operator and prove its spectrum has a gap, using the Källén–Lehmann support as the bridge.
- **Missing theorem: nonperturbative gauge fixing / BRST.** Feynman 1963 and DeWitt 1964 are perturbative; the Gribov ambiguity (1978) shows global obstructions. Attack surface: a gauge-invariant quantization that avoids Gribov copies and still yields a gap.
- **Missing theorem: gauge Goldstone theorem.** Goldstone's theorem is global. A gauge version must prove nonperturbatively that no massless scalar survives and that the would-be Goldstone mode is eaten. Attack surface: a “Higgs theorem” without perturbation theory, or a proof that the pure-YM vacuum has no such mode.
- **Missing technique: uniform transfer-matrix gap.** The lattice formulation gives a gap at finite spacing; what is needed is a lower bound uniform in the lattice spacing and volume. Attack surface: a spectral-gap bound for the Wilson transfer matrix uniform in \(\beta\) and \(a\).
- **Missing theorem: dimensional transmutation in pure YM.** The classical theory has no scale; the gap must scale as \(\Lambda\). Attack surface: prove the continuum limit and the scaling of the gap, not just the existence of a finite-spacing gap.
- **Missing theorem: unique vacuum and absence of symmetry breaking in pure YM.** The \(\theta\)-vacuum and center symmetry complicate the infrared. Attack surface: prove a unique vacuum and a positive gap in the thermodynamic limit.
- **Missing theorem: confinement of the elementary gluon.** A massless gauge field can be confined; the physical spectrum can be gapped even if the elementary field is not. Attack surface: prove a gauge-invariant mass gap without assuming the gluon is massive; separate “elementary mass” from “physical gap.”
- **Missing technique: rigorous lower bound for the Hamiltonian spectral edge.** The Källén–Lehmann representation identifies the gap with the lower edge of a measure but gives no computable bound. Attack surface: positivity / reflection positivity / bootstrap bounds on the spectral measure.

---

## 5. Verification ledger

| Ref | Status | Source checked |
|---|---|---|
| E01 Ehresmann 1950 | `UNVERIFIED` | Search budget exhausted; standard citation only; no page range or DOI asserted. |
| E02 Steenrod 1951 | `UNVERIFIED` | Search budget exhausted; standard book citation; no DOI. |
| E03 Proca 1936 | `VERIFIED` | Journal de Physique article page; DOI `10.1051/jphysrad:0193600708034700`; zbMATH JFM 62.1633.03. |
| E04 Stueckelberg 1938 | `VERIFIED` | Helv. Phys. Acta **11**, 225–244; e-periodica / multiple citations; no DOI. |
| E05 Källén 1952 | `VERIFIED` | e-periodica DOI `10.5169/seals-112316`; INSPIRE; MR0051156. |
| E06 Pauli 1953 | `PARTIAL` | Springer Vol. IV/II ISBN 978-3-540-64312-8; Straumann arXiv:0810.2213 (secondary); specific letters not checked. |
| E07 Lehmann 1954 | `VERIFIED` | Springer DOI `10.1007/BF02783624`; ADS 1954NCim...11..342L. |
| E08 Yang–Mills 1954 | `VERIFIED` | APS DOI `10.1103/PhysRev.96.191`; INSPIRE Phys.Rev. 96, 191; vol./pages confirmed. |
| E09 Shaw 1955 | `VERIFIED` (bibliographic) | Mathematics Genealogy Project; historical sources; thesis unpublished, no DOI. |
| E10 Utiyama 1956 | `VERIFIED` (vol./pages) | INSPIRE record; Phys. Rev. **101**, 1597–1607; DOI not checked. |
| E11 Wightman 1956 | `VERIFIED` | APS DOI `10.1103/PhysRev.101.860`; zbMATH Zbl 0074.22902; MR0083914. |
| E12 BCS 1957 | `VERIFIED` | APS DOI `10.1103/PhysRev.108.1175`; vol./pages confirmed. |
| E13 Sakurai 1960 | `VERIFIED` | Elsevier DOI `10.1016/0003-4916(60)90126-3`; CAPES/Elsevier record. |
| E14 Goldstone 1961 | `VERIFIED` | Springer DOI `10.1007/BF02812722`; ADS 1961NCim...19..154G. |
| E15 NJL 1961 | `VERIFIED` | APS DOI `10.1103/PhysRev.122.345`; vol./pages confirmed. |
| E16 Glashow–Gell-Mann 1961 | `VERIFIED` | CaltechAUTHORS record; DOI `10.1016/0003-4916(61)90193-2`; zbMATH 0099.43601. |
| E17 Glashow 1961 | `VERIFIED` | Elsevier DOI `10.1016/0029-5582(61)90469-2`; Scopus record. |
| E18 Schwinger 1962 I | `VERIFIED` | APS DOI `10.1103/PhysRev.125.397`; zbMATH 0099.23002. |
| E19 Schwinger 1962 II | `VERIFIED` | APS DOI `10.1103/PhysRev.128.2425`; zbMATH 0118.44001. |
| E20 GSV 1962 | `VERIFIED` | APS DOI `10.1103/PhysRev.127.965`; INSPIRE/CERN citations. |
| E21 Feynman 1963 | `VERIFIED` | INSPIRE `Acta Phys. Polon. 24, 697`; page range 697–722; no DOI. |
| E22 Anderson 1963 | `VERIFIED` | APS DOI `10.1103/PhysRev.130.439`; vol./pages confirmed. |
| E23 Englert–Brout 1964 | `VERIFIED` | APS DOI `10.1103/PhysRevLett.13.321`; Crossref; ADS 1964PhRvL..13..321E. |
| E24 Higgs 1964 | `VERIFIED` | APS DOI `10.1103/PhysRevLett.13.508`; APS export record. |
| E25 GHK 1964 | `VERIFIED` | APS DOI `10.1103/PhysRevLett.13.585`; vol./pages confirmed. |
| E26 DeWitt 1964 | `VERIFIED` | APS DOI `10.1103/PhysRevLett.12.742`; INSPIRE record 44821. |
| E27 Salam–Ward 1964 | `VERIFIED` | Elsevier DOI `10.1016/0031-9163(64)90711-5`; vol./pages confirmed. |
| E28 Streater–Wightman 1964 | `VERIFIED` (book + snippet) | Google Books snippet (Princeton reprint, ~p. 205); Clay-problem references; no DOI. |
| E29 Feynman Lectures 1995 | `VERIFIED` | Library catalogs (BnF, ESO); ISBN 0-201-62734-5; course 1962–63. |
| E30 Higgs 1966 | `VERIFIED` | APS DOI `10.1103/PhysRev.145.1156`; vol./pages confirmed. |
| E31 Faddeev–Popov 1967 | `VERIFIED` | Elsevier DOI `10.1016/0370-2693(67)90067-6`; INSPIRE. |

**Counts.** 31 entries total: 29 in the strict 1954–65 window (plus prehistory), 2 boundary entries (1966, 1967). 28 fully `VERIFIED`, 1 `PARTIAL` (Pauli letters), 2 `UNVERIFIED` (Ehresmann, Steenrod; search budget exhausted). No DOI, volume, or page range was invented; where a DOI could not be checked it is omitted and the entry says so.
