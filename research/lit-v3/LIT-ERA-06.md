# LIT-ERA-06 — 1985–1995: Confinement mechanisms and lattice phenomenology

**Window.** Core window 1985–1995. Pre-1985 foundational anchors are included where the era's mechanisms are unintelligible without them ('t Hooft 1978/1979/1981, Gribov 1978, Kugo–Ojima 1979, Mack–Petkova 1980, Wilson 1974, Polyakov 1978). Two immediately-post-window endpoints named in the brief are included and flagged: Del Debbio–Faber–Greensite–Olejník 1997 and Morningstar–Peardon 1999.

**Purpose.** This is the era in which the *phenomenology* of the mass gap was fixed (glueball spectrum, string tension ratio, deconfinement temperature, running coupling) while the *mechanisms* remained either gauge-dependent models or strong-coupling/hierarchical theorems. The file separates these three registers explicitly.

**Label key.**
- **[T]** rigorous theorem/result, unconditional within its stated model or regime.
- **[C]** conditional result: theorem modulo an explicit unproven hypothesis, or a criterion/equivalence whose premise is open.
- **[NG]** no-go / obstruction.
- **[O]** numerical or lattice observation/measurement; no rigorous continuum error control.
- **[V]** model / view / scenario: a mechanism proposal not derived from the Yang–Mills measure.

**Verification key.** VERIFIED = bibliographic record checked this session against a publisher or indexing record via WebSearch (including incidental confirmation). UNVERIFIED = standard citation not re-checked this session (the session WebSearch budget was exhausted after the eleven checks below); it is flagged, not asserted as checked. No citation below is fabricated; uncertain details are marked in the ledger notes.

**Numbering is thematic, not chronological.** E06-01…E06-35 run through seven groups.

---

## 1. The logical target: what a "confinement mechanism" would have to prove

### 1.1 Definitions

Work on a finite hypercubic lattice Λ with gauge group G (here SU(N)), Wilson action, and a positive transfer matrix T = e^{-aH} on the gauge-invariant Hilbert space H_phys (or on BRST cohomology in the continuum). The **mass gap** is the statement

  spec(H_phys) = {0} ∪ [Δ, ∞),  Δ > 0,  vacuum simple,

equivalently: for every local gauge-invariant operator O,

  |⟨O(x)O(0)⟩_c| ≤ C e^{-Δ|x|}   (Euclidean, connected),

with C, Δ **uniform in the volume** and surviving the continuum limit a → 0. On a finite lattice the transfer matrix is finite-dimensional, so a gap always exists trivially; the content is *uniformity* in volume and in the cutoff, plus the existence of the continuum limit. This is the precise sense in which "the mass gap is a lattice observable": one measures the exponential decay rate of the lightest gauge-invariant correlator and must show the rate has a positive limit.

A **confinement mechanism** is a dynamical explanation of the area law

  ⟨W(C)⟩ ~ exp(−σ · Area(C)),  σ > 0,

for large fundamental-representation Wilson loops, equivalently a static potential V(R) ~ σR.

### 1.2 What is *not* equivalent

- **Area law ⇏ mass gap.** The area law constrains one family of observables (static sources). It does not bound the spectral density of glueballs. A confining sector can in principle coexist with massless gauge-invariant excitations; only a separate spectral statement excludes them.
- **Mass gap ⇏ area law.** A Higgs phase has a gap and no confinement. Conversely the deconfined high-temperature phase has a gap (electric screening mass) and a perimeter law.
- **Confinement of static quarks ⇏ mass gap.** Confinement is a statement about the potential between external probes; the mass gap is a statement about the autonomous glue sector.

So the slogan "confinement ⇒ mass gap" is not a theorem. It becomes one only after adding a clustering/spectral input.

### 1.3 Sufficient routes

A mechanism yields a mass gap if it proves at least one of:

- **(S1) Spectral route.** A uniform positive lower bound on the first excited eigenvalue of T (or on the propagator of every gauge-invariant operator), e.g. via Perron–Frobenius plus a variational estimate.
- **(S2) Clustering route.** Exponential decay of *all* connected gauge-invariant correlators with a common rate Δ > 0, uniform in volume and cutoff; Osterwalder–Schrader reconstruction then gives a Wightman theory with a mass gap.
- **(S3) Defect/Peierls route.** A free-energy lower bound for the relevant defect (vortex, monopole, Dirac string) that dominates every Wilson loop and every gauge-invariant observable; combined with a cluster expansion this gives area law *and* exponential clustering.

Each route must additionally survive the continuum limit: the bound must be uniform along a renormalization-group trajectory a → 0. None of the 1985–1995 proposals does all of this in 4D SU(N).

### 1.4 Verdict per mechanism (details in §2 and the scorecard in §3)

| Mechanism | What it actually proves/establishes | Does it yield a mass gap? |
|---|---|---|
| Center vortices ('t Hooft 1978/1979; Mack–Petkova 1980; Del Debbio et al. 1997) | Center flux algebra [T]; a sufficient condition for area law from vortex condensation [C]; lattice center dominance [O] | **No.** Kinematics + a conditional Peierls bound + numerics. No 4D construction, no spectral bound. |
| Dual superconductivity / monopole condensation ('t Hooft 1981; KSW 1987; Suzuki 1988; Suganuma et al. 1995) | Gauge-dependent monopole condensation; dual Ginzburg–Landau linear potential | **No.** Model physics [V] plus gauge-fixed numerics [O]; no gauge-invariant derivation. |
| Gribov–Zwanziger (1989–1994) + Dell'Antonio–Zwanziger | Horizon geometry [T]/[NG]; local action [C]; gluon IR suppression [C] | **No.** Statements concern gauge-fixed colored correlators, not the gauge-invariant spectrum. |
| Kugo–Ojima (1979) | Criterion u(0) = −1 for color confinement | **No.** A conditional equivalence [C] whose premise (no Gribov copies) is false in general; silent on the glueball mass. |
| Strong-coupling / hierarchical / 3D rigor (Wilson 1974; Seiler 1982; Ito 1985; Balaban 1985) | Area law and clustering in restricted regimes | **No.** Genuine [T], but strong-coupling, hierarchical, or 3D — not the 4D continuum. |
| Lattice spectrum (Michael–Tickle 1988; Bali et al. 1993; Morningstar–Peardon 1999) | Δ ≈ m(0++) ≈ 1.7 GeV, m/√σ ≈ 3.5–4 | Gives the **number**, not a proof or a mechanism. |

---

## 2. Entries

### A. Center vortices and the vortex picture

**E06-01. [V] G. 't Hooft, "On the phase transition towards permanent quark confinement," Nucl. Phys. B 138, 1–25 (1978).**
Proposes the center-vortex condensation picture: confinement is a phase in which Z(N) magnetic flux vortices percolate through the vacuum. This is a scenario, not a derivation. It supplies the target later made precise by Mack–Petkova and tested on the lattice. No statement about the glueball mass.

**E06-02. [T] G. 't Hooft, "A property of electric and magnetic flux in non-Abelian gauge theories," Nucl. Phys. B 153, 141–160 (1979).**
Defines 't Hooft loops (magnetic flux) and proves the commutation algebra between Wilson (electric) and 't Hooft (magnetic) loop operators, with Z(N) phases when the loops link. This is a rigorous *kinematic* result about the center. It constrains which phases are possible; it does not decide which phase 4D SU(N) realizes and gives no mass gap.

**E06-03. [C] G. Mack and V. B. Petkova, "Sufficient condition for confinement of static quarks by a vortex condensation mechanism," Ann. Phys. 125, 117–134 (1980).**
Proves a sufficient condition: if center vortices of a compact gauge group with nontrivial center condense with a free-energy cost scaling with the enclosed volume (a Peierls-type bound), then static quarks confine with a linear potential. This is the cleanest rigorous-style route from a vortex picture to an area law. It is a conditional theorem in a model setting: it does not construct the 4D continuum theory and does not prove a mass gap. This is the single most load-bearing entry in group A for v3.

**E06-04. [O] L. Del Debbio, M. Faber, J. Greensite, Š. Olejník, "Center dominance and Z2 vortices in SU(2) lattice gauge theory," Phys. Rev. D 55, 2298–2305 (1997); arXiv:hep-lat/9610005.**
Introduces center projection and shows that center-projected (vortex-only) configurations reproduce the full SU(2) string tension — "center dominance." Strong numerical evidence for the vortex picture, but the projection is gauge-dependent and there is no continuum control. Immediately post-window; included because the brief names it as the vortex endpoint.

### B. Dual Ginzburg–Landau and monopole condensation

**E06-05. [V] G. 't Hooft, "Topology of the gauge condition and new confinement phases in non-Abelian gauge theories," Nucl. Phys. B 190 [FS3], 455–478 (1981).**
Introduces Abelian projection: after a partial gauge fixing leaving the maximal torus, SU(N) Yang–Mills looks like an Abelian gauge theory with magnetic monopoles. Proposes monopole condensation as the confinement mechanism. The projection is gauge-dependent, and monopole condensation is not manifestly gauge-invariant — the central weakness of the entire programme.

**E06-06. [O] A. S. Kronfeld, G. Schierholz, U.-J. Wiese, "Topology and dynamics of the confinement mechanism," Nucl. Phys. B 293, 461–478 (1987).**
Performs the Abelian projection of SU(N) in continuum and on the lattice and measures magnetic-monopole currents, finding monopole condensation signals correlated with confinement. Numerical, in a fixed (maximal Abelian) gauge; not a gauge-invariant proof.

**E06-07. [O] A. S. Kronfeld, M. L. Laursen, G. Schierholz, U.-J. Wiese, "Monopole condensation and color confinement," Phys. Lett. B 198, 516 (1987).**
Follow-up lattice study reporting monopole condensation and its relation to color confinement in SU(2). Numerical evidence within the same gauge-dependent framework.

**E06-08. [V] T. Suzuki, "A Ginzburg–Landau type theory of quark confinement," Prog. Theor. Phys. 80, 929–934 (1988).**
Proposes an effective dual Ginzburg–Landau theory in which condensed monopoles produce a dual Meissner effect and a linear potential. The string tension is an input of the dual Higgs scale, not derived from the SU(N) measure.

**E06-09. [V] H. Suganuma, S. Sasaki, H. Toki, "Color confinement, quark pair creation and dynamical chiral-symmetry breaking in the dual Ginzburg–Landau theory," Nucl. Phys. B 435, 207–240 (1995).**
Develops the dual Ginzburg–Landau theory with QCD-monopole condensation, deriving a linear potential, a vortex string-tension formula, quark-pair creation, and dynamical chiral symmetry breaking. A model: the monopole condensate is an input. It is the most complete dual-superconductor phenomenology of the era but proves nothing about the SU(3) measure.

### C. Gribov–Zwanziger framework and the Kugo–Ojima criterion

**E06-10. [NG] V. N. Gribov, "Quantization of non-Abelian gauge theories," Nucl. Phys. B 139, 1–10 (1978).**
Shows that the Faddeev–Popov determinant vanishes on the Gribov horizon, so no single continuous gauge slice exists and perturbative gauge fixing is incomplete. An obstruction to naive BRST/FP quantization and the origin of the Gribov ambiguity. It does not by itself give confinement or a mass gap; it is the reason the Kugo–Ojima premise is not free.

**E06-11. [C] T. Kugo and I. Ojima, "Local covariant operator formalism of non-Abelian gauge theories and quark confinement problem," Prog. Theor. Phys. Suppl. 66, 1–130 (1979).**
Derives the BRST quartet mechanism and the confinement criterion u(0) = −1 (equivalently a 1/p⁴ ghost propagator), under the assumptions that the FP operator has no zero modes (no Gribov copies) and that the physical space is the BRST cohomology. It is a criterion/equivalence, not a proof: the premise is exactly what Gribov showed can fail, and the criterion does not by itself give a glueball mass gap. (Erratum: Prog. Theor. Phys. 71 (1984) 1121.)

**E06-12. [C] D. Zwanziger, "Local and renormalizable action from the Gribov horizon," Nucl. Phys. B 323, 513–544 (1989).**
Derives the Gribov–Zwanziger local, renormalizable action by restricting the functional integral to the interior of the Gribov horizon, introducing a dimensionful parameter γ. The action is local and perturbatively renormalizable, but the restriction is a gauge-fixing prescription; the physical spectrum and any mass gap are not derived. Gluon propagator poles shift to unphysical p² = ±iγ^{1/2}; the glueball channel is argued to retain a physical cut with positive spectral function.

**E06-13. [NG] G. Dell'Antonio and D. Zwanziger, "Ellipsoidal bound on the Gribov horizon contradicts the perturbative renormalization group," Nucl. Phys. B 326, 333–350 (1989).**
Proves an ellipsoidal bound on the Gribov horizon and shows it is incompatible with the perturbative renormalization-group asymptotics of the gluon two-point function in D = 4. A no-go: the Gribov region cannot be treated as a small perturbation of the perturbative vacuum.

**E06-14. [T] G. Dell'Antonio and D. Zwanziger, "Every gauge orbit passes inside the Gribov horizon," Commun. Math. Phys. 138, 291–299 (1991).**
Proves, using L² topology on connections and a natural topology on local gauge transformations, that every gauge orbit intersects the Gribov region. A rigorous geometric result. It justifies restricting to the Gribov region but neither selects a unique representative nor proves confinement.

**E06-15. [C] D. Zwanziger, "Vanishing of zero-momentum lattice gluon propagator and color confinement," Nucl. Phys. B 364, 127–161 (1991).**
Defines minimal Landau gauge on the lattice and argues that the zero-momentum gluon propagator vanishes in the infinite-volume limit, implying maximal violation of reflection positivity for colored gluons. An analytic argument resting on lattice-motivated hypotheses; it concerns gauge-fixed colored correlators, not the gauge-invariant mass gap.

**E06-16. [C] D. Zwanziger, "Fundamental modular region, Boltzmann factor and area law in lattice gauge theory," Nucl. Phys. B 412, 127–142 (1994).**
Argues that the fundamental modular region plus a Boltzmann-factor suppression of the horizon yields an area law for Wilson loops. A conditional analytic argument in lattice gauge theory, not a construction of the 4D continuum and not a mass-gap proof.

### D. Lattice gluon propagator and running coupling

**E06-17. [O] P. Marenzoni, G. Martinelli, N. Stella, "The gluon propagator in momentum space," Phys. Lett. B 318, 511–516 (1993).**
Measures the gluon propagator in momentum space in lattice Landau gauge and finds IR suppression relative to the perturbative form. A gauge-fixed numerical result; the IR extrapolation and continuum limit are not controlled. Author list/venue should be re-checked before quoting (see ledger).

**E06-18. [O] G. P. Lepage and P. B. Mackenzie, "On the viability of lattice perturbation theory," Phys. Rev. D 48, 2250–2264 (1993).**
Analyzes the failure of naive lattice perturbation theory and introduces tadpole improvement, matching lattice operators to continuum schemes. A method validated numerically; essential for continuum coupling extraction, but not a dynamical result.

**E06-19. [T] M. Lüscher, R. Narayanan, P. Weisz, U. Wolff, "The Schrödinger functional — a renormalizable probe for non-Abelian gauge theories," Nucl. Phys. B 384, 168–228 (1992).**
Introduces the Schrödinger functional as a gauge-invariant, renormalizable probe of non-Abelian gauge theories with Dirichlet boundaries and proves its perturbative renormalizability. The basis of the precision running-coupling measurements that follow.

**E06-20. [O] R. Sommer, "A new way to set the energy scale in lattice gauge theories and its application to the static force and α_s in SU(2) Yang–Mills theory," Nucl. Phys. B 411, 839–854 (1994).**
Introduces the r₀ scale from the static quark potential and uses it to set the scale and extract α_s in SU(2) Yang–Mills. A precise, widely used scale-setting method; not a mass-gap measurement.

**E06-21. [O] M. Lüscher, R. Sommer, P. Weisz, U. Wolff, "A precise determination of the running coupling in the SU(3) Yang–Mills theory," Nucl. Phys. B 413, 481–502 (1994).**
Determines the running coupling in SU(3) Yang–Mills non-perturbatively with the Schrödinger functional and a controlled continuum extrapolation. Demonstrates the SF scheme works; does not address the IR mass gap.

### E. Finite-temperature deconfinement

**E06-22. [V] A. M. Polyakov, "Thermal properties of gauge fields and quark liberation," Phys. Lett. B 72, 477–480 (1978).**
Introduces the Polyakov loop as the order parameter for finite-temperature deconfinement and argues that the center symmetry Z(N) is spontaneously broken in the deconfined phase. Foundational scenario; the argument is strong-coupling/mean-field, not a proof for 4D SU(N).

**E06-23. [O] L. D. McLerran and B. Svetitsky, "Quark liberation at high temperature: A Monte Carlo study of SU(2) gauge theory," Phys. Rev. D 24, 450–460 (1981).**
Monte Carlo study of SU(2) lattice gauge theory at finite temperature; locates the deconfinement transition and measures the Polyakov-loop susceptibility. Numerical evidence for the Z(2) transition; no continuum extrapolation.

**E06-24. [C] B. Svetitsky and L. G. Yaffe, "Critical behavior at finite-temperature confinement transitions," Nucl. Phys. B 210 [FS6], 423–447 (1982).**
Argues by universality that the finite-temperature deconfinement transition of a gauge theory with center Z(N) is in the universality class of a Z(N) spin model in d−1 dimensions (for SU(2): 3D Ising). A renormalization-group universality hypothesis supported by the loop algebra, not a proof for every N_t.

**E06-25. [O] J. Fingberg, U. Heller, F. Karsch, "Scaling and asymptotic scaling in the SU(2) gauge theory," Nucl. Phys. B 392, 493–517 (1993).**
Studies the scaling of the SU(2) deconfinement temperature T_c across N_t = 4, 6, 8, … and extracts β_c(N_t); tests asymptotic scaling. Numerical; establishes the N_t-dependent crossover and its scaling behavior. The SU(3) analogue gives T_c ≈ 270 MeV in the quenched theory (scale-dependent).

### F. Glueball spectrum and the mass gap as a lattice observable

**E06-26. [O] C. Michael and G. A. Tickle, "The glueball spectrum in SU(3) lattice gauge theory," Phys. Lett. B 207, 313–317 (1988).**
Lattice computation of the low-lying SU(3) glueball spectrum (0++, 2++, …) in the quenched approximation; establishes the 0++ as the lightest glueball and gives m(0++) in units of √σ. Sets the phenomenological target for any mass-gap mechanism.

**E06-27. [O] G. S. Bali, K. Schilling, A. Hulsebos, A. C. Irving, C. Michael, P. W. Stephenson, "A comprehensive lattice study of SU(3) glueballs," Phys. Lett. B 309, 378–384 (1993).**
Comprehensive quenched SU(3) glueball study on large lattices; refines the 0++, 2++, and excited glueball masses and the m/√σ ratios. Precision phenomenology, still with systematic (continuum, quenching) uncertainties.

**E06-28. [O] C. J. Morningstar and M. Peardon, "Glueball spectrum from an anisotropic lattice study," Phys. Rev. D 60, 034509 (1999); arXiv:hep-lat/9901004.**
Precision anisotropic-lattice computation of the quenched glueball spectrum: m(0++) ≈ 1.7 GeV, m(0++)/√σ ≈ 3.5–4. The modern benchmark; still quenched and not a mechanism. Immediately post-window; included as the named endpoint.

**E06-29. [T] M. Lüscher, "Some analytic results concerning the mass spectrum of Yang–Mills gauge theories on a torus," Nucl. Phys. B 219, 233–261 (1983).**
Analytic results on the mass spectrum on a torus: the transfer matrix has a discrete spectrum and masses are extracted from exponential decay of correlators. A finite-volume spectral framework that *assumes* the theory is in the massive phase; it does not prove the gap. Pre-window anchor for the "mass gap as observable" line.

**E06-30. [T] M. Lüscher, "Volume dependence of the energy spectrum in massive quantum field theories. I. Stable particle states," Commun. Math. Phys. 104, 177–206 (1986).**
Derives the finite-volume mass formula relating the exponential volume dependence of a stable particle's energy to its infinite-volume mass and forward scattering amplitude. Conditional on the existence of the massive particle; the standard tool for extracting masses, not a proof of the gap.

**E06-31. [T] M. Lüscher, "Volume dependence of the energy spectrum in massive quantum field theories. II. Scattering states," Commun. Math. Phys. 105, 153–188 (1986).**
Extends the finite-volume analysis to multi-particle/scattering states, giving the quantization condition for two-particle levels. Rigorous in the same conditional sense; essential for the modern precision spectrum but not a gap proof.

### G. Rigorous strong-coupling, hierarchical, and 3D results

**E06-32. [T] K. G. Wilson, "Confinement of quarks," Phys. Rev. D 10, 2445–2459 (1974).**
Defines lattice gauge theory with the Wilson action and proves, by a convergent strong-coupling (character) expansion for sufficiently small β = 2N/g², an area law for fundamental Wilson loops. A genuine theorem in the strong-coupling regime and the origin of the "confinement = area law" identification. It does not address the continuum limit or a mass gap; the strong-coupling expansion has finite radius and does not reach the scaling regime.

**E06-33. [T] E. Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics*, Lecture Notes in Physics 159 (Springer, Berlin, 1982).**
Develops lattice gauge theory as constructive QFT: reflection positivity, transfer matrix, cluster expansions, and strong-coupling results including area law and exponential clustering. Rigorous for the strong-coupling regime; the weak-coupling/continuum limit is not covered. The canonical entry point for what "rigorous" meant in this era.

**E06-34. [T] K. R. Ito, "Permanent quark confinement in four-dimensional hierarchical lattice gauge theories of Migdal–Kadanoff type," Phys. Rev. Lett. 55, 558–561 (1985).**
Proves permanent quark confinement for 4D *hierarchical* lattice gauge theories of Migdal–Kadanoff type, for SU(N) and U(N), using the exact MK recursion. A theorem in the hierarchical model, **not** the true 4D Wilson lattice. This distinction is frequently lost in citations: the result does not transfer to the Clay model.

**E06-35. [T] T. Balaban, "Ultraviolet stability of three-dimensional lattice pure gauge field theories," Commun. Math. Phys. 102, 255–275 (1985).**
Rigorous renormalization-group analysis of 3D lattice pure gauge theories: ultraviolet stability of the effective actions under block-spin transformations. A theorem in 3D, part of the program toward a rigorous continuum limit; not a 4D mass-gap proof.

---

## 3. Scorecard: rigorous vs model vs numerical

| Proposal / line | What is actually established | Register | Yields mass gap? |
|---|---|---|---|
| Center vortex picture (E06-01, 02, 03, 04) | Center algebra [T]; sufficient condition for area law [C]; lattice center dominance [O] | theorem + conditional + numerical | No: area law only; no 4D construction or spectral bound |
| Dual superconductivity (E06-05…09) | Gauge-dependent monopole condensation; dual GL linear potential | model + numerical | No: model physics; not gauge-invariant |
| Gribov–Zwanziger (E06-10…16) | Horizon geometry [T]/[NG]; local action [C]; gluon suppression [C] | theorem (geometry) + conditional (dynamics) | No: colored correlators only |
| Kugo–Ojima (E06-11) | Criterion u(0) = −1 | conditional equivalence | No: criterion, not proof |
| Lattice propagator/coupling (E06-17…21) | Scheme-defined running coupling, r₀, IR suppression | numerical (+ one perturbative theorem) | No: not the gauge-invariant spectrum |
| Finite T (E06-22…25) | Polyakov-loop order parameter; β_c(N_t); universality class | model + numerical | No: deconfinement scale, not Δ |
| Glueball spectrum (E06-26…31) | m(0++) ≈ 1.7 GeV; m/√σ ≈ 3.5–4 | numerical | Gives the number, not a proof |
| Strong-coupling / hierarchical / 3D (E06-32…35) | Area law, clustering, confinement in restricted models | rigorous | No: restricted regime only |

**Era synthesis.** The era fixed the target quantitatively (Δ ≈ 1.7 GeV, Δ/√σ ≈ 3.5–4, T_c ≈ 270 MeV, a precise running coupling) and produced the two enduring mechanism *languages* (center vortices; dual superconductivity) plus two *frameworks* (Gribov–Zwanziger; Kugo–Ojima). It produced **no** rigorous 4D continuum mass gap. Every rigorous theorem is confined to strong coupling, a hierarchical model, 3D, or finite-volume spectral technology. Every 4D mechanism is either a gauge-dependent model or a conditional Peierls-type bound.

---

## 4. Position in 14 eras

The campaign's 14-era division (working periodization, not a claim about sharp historical breaks). Era 06 is this file.

1. **1954–1972 — Non-Abelian gauge theory founded.** Yang–Mills, quantization, the renormalizability question; no confinement problem yet.
2. **1973–1975 — Asymptotic freedom and dimensional transmutation.** The infrared problem and the mass-gap question become sharp.
3. **1976–1978 — Instantons, θ-vacua, 't Hooft–Polyakov monopoles.** Early dual-superconductor and vortex ideas appear.
4. **1978–1980 — Gribov copies; center algebra; lattice gauge theory.** Wilson's strong-coupling confinement; the FP obstruction identified.
5. **1981–1984 — Monte Carlo era.** String tension, scaling, Savvidy instability; first continuum-limit attempts.
6. **1985–1995 — Confinement mechanisms and lattice phenomenology.** ← **this file.** Vortices, dual SC, GZ/KO, propagator, finite T, glueballs; no rigorous 4D gap.
7. **1995–2000 — Abelian dominance and vortex detection.** Dual-superconductivity precision; first Gribov–Zwanziger propagator studies; center dominance.
8. **2000–2004 — Refined GZ and functional methods.** IR propagators, Kugo–Ojima tests; DSE/FRG mature.
9. **2005–2008 — Scaling vs decoupling.** Ghost dominance debate; gauge-invariant order parameters.
10. **2008–2011 — Precision spectrum and finite-T continuum extrapolations.** Large-N lattice; controlled systematics.
11. **2011–2014 — Rigorous lattice RG and gauge-fixing rigor.** Constructive progress; OS reconstruction advances.
12. **2014–2018 — Gradient flow and rigorous continuum limits.** 2D/3D constructions; partial 4D results; stochastic quantization.
13. **2018–2022 — Constructive and weak-coupling partial theorems.** 't Hooft-limit and related mass-gap results.
14. **2022–2026 — This campaign.** Conditional theorems, H1/f₁² residue transport, numerical + Lean.

---

## 5. Gaps and attack surfaces for v3

- **No era mechanism produces a gauge-invariant spectral bound.** The missing step is (S1)/(S2), not another area-law derivation. Attack: use the project's exact variance decomposition G″ = S − D to turn a vortex/monopole free-energy bound into a lower bound on a gauge-invariant connected correlator; the defect cost must dominate D, not merely S.
- **Mack–Petkova is the only clean sufficient condition, and it is a Peierls bound in a model.** Attack: formulate a continuum-compatible vortex free energy using the center-probe/'t Hooft-loop algebra (E06-02), and prove the bound is uniform as a → 0. That uniformity is precisely the missing step.
- **Abelian projection is gauge-dependent (E06-05…07).** Attack: build a gauge-invariant monopole order parameter (a 't Hooft loop, or a monopole dressed by a Dirac string to infinity) and test condensation with the project's gauge-invariant IR mass. If the condensate survives only in maximal Abelian gauge, it is an artifact.
- **Kugo–Ojima's premise is false in the presence of Gribov copies (E06-10, 11).** Attack: derive a modified criterion inside the GZ ensemble and test whether u(0) = −1 is compatible with the horizon restriction and with a positive glueball mass; a no-go here is publishable.
- **GZ predicts colored positivity violation but is silent on the gauge-invariant spectrum (E06-12, 15, 16).** Attack: prove either (a) the GZ spectral function has a positive-mass pole in the glueball channel, or (b) GZ cannot produce a physical mass gap. Either closes a loophole.
- **The lattice Landau-gauge propagator is scheme- and gauge-dependent, and its IR extrapolation was uncontrolled in this era (E06-17).** Attack: replace it with the project's gauge-invariant IR mass and cutoff-consistent G″(0); match to the Schrödinger-functional running coupling (E06-21) at a scale where both are defined.
- **Finite-T β_c(N_t) gives a deconfinement scale, not Δ (E06-22…25).** Attack: use Svetitsky–Yaffe universality as a consistency check; any candidate order parameter must reproduce the Z(N) universality class, which is a sharp falsification test.
- **The precision ratio m(0++)/√σ ≈ 3.5–4 is the target (E06-26…28).** Attack: a mechanism that only gives σ > 0 cannot explain the ratio. Derive it from the f₁² residue transport / H1; this is the concrete quantitative discriminator between mechanisms.
- **The strong-coupling cluster expansion is rigorous but stops at β < β₀ (E06-32, 33); the 4D RG completion is the actual missing "machine" (E06-34, 35).** Attack: this is exactly the charter's demanded object — a cutoff-consistent, non-perturbative transport of the strong-coupling gap across the Savvidy crossover. Balaban's RG is the closest rigorous precedent; its 4D completion is the step to isolate.
- **Bibliographic honesty trap: Ito 1985 is hierarchical, not 4D lattice SU(2).** Do not cite it as a 4D theorem in v3. Use it only as a model check of the MK recursion, and state the scope explicitly.

---

## 6. Verification ledger

| # | Entry (short) | Checked via | Status | Note |
|---|---|---|---|---|
| E06-01 | 't Hooft, NPB 138 (1978) 1 | not re-checked | UNVERIFIED | standard citation |
| E06-02 | 't Hooft, NPB 153 (1979) 141 | not re-checked | UNVERIFIED | standard citation |
| E06-03 | Mack–Petkova, Ann. Phys. 125 (1980) 117 | WebSearch (ScienceDirect/ADS) | VERIFIED | pages 117–134; DOI 10.1016/0003-4916(80)90121-9 |
| E06-04 | Del Debbio et al., PRD 55 (1997) 2298 | not re-checked | UNVERIFIED | arXiv:hep-lat/9610005 |
| E06-05 | 't Hooft, NPB 190 (1981) 455 | not re-checked | UNVERIFIED | standard citation |
| E06-06 | KSW, NPB 293 (1987) 461 | WebSearch (Elsevier/Google Scholar) | VERIFIED | pages 461–478; DOI 10.1016/0550-3213(87)90080-0 |
| E06-07 | KLSW, PLB 198 (1987) 516 | WebSearch (follow-up confirmed) | VERIFIED | start page 516 confirmed; end page not re-checked |
| E06-08 | Suzuki, PTP 80 (1988) 929 | not re-checked | UNVERIFIED | standard citation |
| E06-09 | Suganuma–Sasaki–Toki, NPB 435 (1995) 207 | WebSearch (ADS/Elsevier) | VERIFIED | pages 207–240; DOI 10.1016/0550-3213(94)00392-R; exact title differs from shorthand |
| E06-10 | Gribov, NPB 139 (1978) 1 | not re-checked | UNVERIFIED | standard citation |
| E06-11 | Kugo–Ojima, PTPS 66 (1979) 1 | WebSearch (OUP/CERN) | VERIFIED | pp. 1–130; DOI 10.1143/PTPS.66.1; erratum PTP 71 (1984) 1121 |
| E06-12 | Zwanziger, NPB 323 (1989) 513 | WebSearch (Elsevier) | VERIFIED | pp. 513–544; DOI 10.1016/0550-3213(89)90122-3 |
| E06-13 | Dell'Antonio–Zwanziger, NPB 326 (1989) 333 | WebSearch (ADS/Semantic Scholar) | VERIFIED | pp. 333–350; the ellipsoidal-bound paper is NPB 326, not CMP 138 |
| E06-14 | Dell'Antonio–Zwanziger, CMP 138 (1991) 291 | WebSearch (Springer) | VERIFIED | pp. 291–299; DOI 10.1007/BF02099494 |
| E06-15 | Zwanziger, NPB 364 (1991) 127 | WebSearch (Elsevier/INSPIRE) | VERIFIED | pp. 127–161; DOI 10.1016/0550-3213(91)90581-H |
| E06-16 | Zwanziger, NPB 412 (1994) 127 | budget exhausted | UNVERIFIED | standard citation; check before quoting |
| E06-17 | Marenzoni–Martinelli–Stella, PLB 318 (1993) 511 | budget exhausted | UNVERIFIED | author list/venue uncertain; check before quoting |
| E06-18 | Lepage–Mackenzie, PRD 48 (1993) 2250 | not re-checked | UNVERIFIED | standard citation |
| E06-19 | LNWW, NPB 384 (1992) 168 | not re-checked | UNVERIFIED | standard citation |
| E06-20 | Sommer, NPB 411 (1994) 839 | not re-checked | UNVERIFIED | standard citation |
| E06-21 | LSWW, NPB 413 (1994) 481 | not re-checked | UNVERIFIED | standard citation |
| E06-22 | Polyakov, PLB 72 (1978) 477 | not re-checked | UNVERIFIED | standard citation |
| E06-23 | McLerran–Svetitsky, PRD 24 (1981) 450 | not re-checked | UNVERIFIED | standard citation |
| E06-24 | Svetitsky–Yaffe, NPB 210 (1982) 423 | not re-checked | UNVERIFIED | standard citation |
| E06-25 | Fingberg–Heller–Karsch, NPB 392 (1993) 493 | budget exhausted | UNVERIFIED | standard citation |
| E06-26 | Michael–Tickle, PLB 207 (1988) 313 | not re-checked | UNVERIFIED | standard citation |
| E06-27 | Bali et al., PLB 309 (1993) 378 | not re-checked | UNVERIFIED | standard citation |
| E06-28 | Morningstar–Peardon, PRD 60 (1999) 034509 | not re-checked | UNVERIFIED | arXiv:hep-lat/9901004 |
| E06-29 | Lüscher, NPB 219 (1983) 233 | not re-checked | UNVERIFIED | standard citation |
| E06-30 | Lüscher, CMP 104 (1986) 177 | not re-checked | UNVERIFIED | standard citation |
| E06-31 | Lüscher, CMP 105 (1986) 153 | not re-checked | UNVERIFIED | standard citation |
| E06-32 | Wilson, PRD 10 (1974) 2445 | not re-checked | UNVERIFIED | standard citation |
| E06-33 | Seiler, LNP 159 (1982) | not re-checked | UNVERIFIED | book; Springer |
| E06-34 | Ito, PRL 55 (1985) 558 | WebSearch (APS/INSPIRE) | VERIFIED | pp. 558–561; **hierarchical MK model, not 4D Wilson lattice** |
| E06-35 | Balaban, CMP 102 (1985) 255 | WebSearch (incidental confirmation) | VERIFIED | 3D UV stability; not a 4D mass-gap proof |

**Counts.** 35 entries. 11 VERIFIED, 24 UNVERIFIED. Of the 24 UNVERIFIED, 2 (E06-16, E06-17) carry explicit field-level uncertainty and should be checked before v3 submission; the rest are standard citations.

---

*End of LIT-ERA-06.*
