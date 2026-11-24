# LIT-ERA-12 — 2018–2021: Resurgence, Transseries, and Modern Nonperturbative Methods

> **Era:** E12 in the campaign's 14-era working index (see §4).
> **Survey date:** 2026-11-24.
> **Scope:** publications 2018–2021, with a small number of adjacent items (2014–2017 preprints, a 2022 review of 2018–21 work) explicitly flagged *adjacent*.
> **Question addressed:** what this era established about the 4D pure SU(N) mass gap, what it only formalised, and what it proved cannot be reached by which methods.
> **Method:** entries were checked against arXiv abstract pages and Crossref DOI records; see the ledger in §6. WebSearch budget for the shared turn was exhausted partway through, so later checks use direct DOI/arXiv metadata only. Nothing below is cited from memory without an explicit `[U]`.

**Label legend** (primary epistemic label + bibliographic tag):

| Label | Meaning |
|---|---|
| **[T]** | rigorous theorem, construction, or exact identity |
| **[C]** | conditional: valid only under stated hypotheses, a truncation, or a semiclassical regime |
| **[NG]** | no-go / obstruction / blindness result |
| **[O]** | observation, numerical result, review, or method (not a theorem about the gap) |
| **[V]** | bibliographic record verified this session |
| **[U]** | not independently verified this session; standard record, treated as provisional |

---

## 0. Executive summary

1. **The era's central conceptual gain is a precise location of the obstruction.** In an asymptotically free theory the mass gap is an essential singularity of the coupling: $\Delta = C\,\Lambda$ with $\Lambda=\mu e^{-1/(2\beta_0g^2)}$. The perturbative series is the Taylor series of a function that is flat at $g^2=0$; it contains *zero* terms of the form $e^{-1/(2\beta_0g^2)}$. Every method that is analytic in $g^2$ at the origin is structurally blind to the gap, not merely inaccurate.
2. **Resurgence sharpens this but does not remove it.** The Borel plane, large-order asymptotics, and Stokes data fix the *exponent* $A$ and the discontinuity; they do **not** fix the transseries parameter $\sigma$ (normalisation and phase), the choice of vacuum, or the real part across an unstable saddle. Since $\Delta=C\Lambda$, resurgence can give the scale $\Lambda$ but not the pure number $C$.
3. **The renormalon picture is the precise statement of non-Borel-summability.** IR renormalons place singularities on the positive Borel axis at $t=u/\beta_0$; the Borel integral is ambiguous by exactly the size of a nonperturbative correction. No rigorous proof or disproof of Borel summability for 4D YM exists in this era.
4. **The gluon condensate remains an input, not an output.** SVZ sum rules and their modern evaluations (Gubler–Satow 2019) treat $\langle\frac{\alpha_s}{\pi}G^2\rangle$ as a fitted/measured quantity with scheme dependence; lattice extractions (Bali–Bauer–Pineda 2014) make it a measured input, not a derivation of the gap.
5. **FRG and DSE became quantitative but remained truncation-bound.** The Dupuis et al. 2021 and Huber 2020 reviews are the era's reference points. Both produce a gauge-fixed IR mass scale; neither produces a gauge-invariant physical mass gap with controlled error. The scaling-vs-decoupling dispute was resolved *phenomenologically* in favour of decoupling on the lattice, not by a theorem.
6. **Gribov–Zwanziger gives a confining gauge-fixed dynamics, not a physical Hilbert space.** The GZ propagator violates reflection positivity; the Gribov parameter is a mass scale fixed by a gap equation, but it is not the physical gap.
7. **Lattice glueball precision and the Wilson flow are the era's strongest nonperturbative numbers.** $m_{0^{++}}/\sqrt{\sigma}\approx3.4$ and continuum extrapolations (Athenodorou–Teper 2020) are the benchmark any proof must reproduce; the continuum limit itself remains the open step.
8. **Gauge-theory bootstrap in 4D did not exist in this window.** The era has lattice loop-equation bootstrap (Anderson–Kruczenski 2017, adjacent) and matrix-QM bootstrap (Han–Hartnoll–Kruthoff 2020); the 4D lattice YM bootstrap lands post-era (2022–23) and is a cross-check, not a result of E12.
9. **Net verdict:** E12 did not prove the mass gap, did not prove Borel summability, and did not construct 4D YM. It did convert "nonperturbative" from a slogan into a sharp statement about essential singularities and Borel singularities — which is exactly the statement v3 must attack.

---

## 1. The precise sense in which the mass gap is an essential singularity

### 1.1 Dimensional transmutation

For pure SU(N) with
$$\beta(g)\;=\;-\beta_0 g^3-\beta_1 g^5-\cdots,\qquad \beta_0=\frac{11N}{48\pi^2},$$
the one-loop solution is
$$\Lambda\;=\;\mu\,\exp\!\left(-\frac{1}{2\beta_0 g^2(\mu)}\right)\bigl(1+O(g^2)\bigr),\qquad
\frac{1}{2\beta_0 g^2}=\frac{8\pi^2}{b_0 g^2},\quad b_0=\frac{11N}{3}.$$
$\Lambda$ is the only mass scale of the pure gauge theory. Every physical mass is
$$\Delta \;=\; C\,\Lambda,$$
with $C$ a pure number (a function of $N$ and the renormalisation scheme, but not of $\mu$). The whole dynamical problem is the number $C$; the scale $\Lambda$ is one-loop and perturbative.

### 1.2 Flatness at $g^2=0$

A gauge-invariant Euclidean correlator at fixed physical momentum has the transseries form
$$G(p;g^2)\;=\;\underbrace{\sum_{n\ge0}a_n\,g^{2n}}_{\text{perturbative}}
\;+\;\sum_s \sigma_s\,e^{-A_s/g^2}\,g^{2\beta_s}\sum_{n\ge0}b_{s,n}\,g^{2n},
\qquad A_s>0 .$$
The second sum is $C^\infty$ in $g^2$ at $0$ and **flat**: for every $k$,
$$\frac{d^k}{d(g^2)^k}\Bigl[e^{-A/g^2}\Bigr]_{g^2=0}=0 .$$
Hence its Taylor series at $g^2=0$ is identically zero. The perturbative expansion is the Taylor series of the first sum only. It follows that:

- the mass gap, carried by the $A=1/(2\beta_0)$ sector (plus running corrections), **vanishes to all orders in perturbation theory**;
- the full $G(g^2)$ has an **essential singularity** at $g^2=0$;
- "the gap is nonperturbative" is not a statement about difficulty; it is the statement that the Taylor series at the origin is identically zero in the relevant sector.

Equivalently, at fixed lattice cutoff and bare coupling $\beta=2N/g^2$, the gap in lattice units is $\Delta a\sim e^{-c\beta}$; the continuum limit $\beta\to\infty$ is approached through an essential singularity in $\beta$.

### 1.3 Borel plane and renormalons

The Borel transform
$$B(t)=\sum_{n\ge0}\frac{a_n}{n!}\,t^n$$
has singularities at $t=A_s$. The Borel sum $\int_0^\infty dt\,e^{-t/g^2}B(t)$ is ambiguous by $\sim e^{-A_s/g^2}$. For IR renormalons the singularities sit on the **positive real axis** at $t=u/\beta_0$, so the integral is not defined without a prescription, and the prescription ambiguity is exactly the size of a nonperturbative correction. This is the precise content of "'t Hooft's renormalon obstruction" [R11]: the perturbative series is not Borel-summable in the sense that would determine the exponential sector. Where Borel summability *is* proved (massive $\varphi^4$; see Costin [R6] for the general theory), the positive-axis singularities are absent; in 4D YM they are expected and unproven.

Large-order asymptotics closes the loop: if $a_n\sim n!\,A^{-n}$, then $B$ has a singularity at $t=A$, and resurgence identifies the discontinuity across it (the Stokes constant) with the same residue that controls the $n!$ growth.

### 1.4 What resurgence fixes and what it does not

| Quantity | Determined by resurgence? | Comment |
|---|---|---|
| Location $A_s$ of Borel singularities | **Yes** | from saddle actions and/or large-order growth |
| Discontinuity / Stokes constant | **Yes**, given the perturbative series | residue at $t=A_s$ |
| Large-order asymptotics $a_n$ | **Yes** | $a_n\sim n!A_s^{-n}\times$ series |
| Transseries parameter $\sigma_s$ (normalisation, phase) | **No** | integration constant; fixed by vacuum/boundary condition |
| Choice of vacuum / lateral summation | **No** | a nonperturbative input |
| Real part across an unstable saddle (Savvidy) | **No** | requires a Lefschetz thimble; the Gaussian saddle is complex |
| $\Delta/\Lambda=C$ | **No** | the pure number is the dynamical datum |

**The "datum vs sign" distinction.** Resurgence fixes the *datum* in the sense of the singularity structure and the discontinuity: the exponent, the large-order relation, and the relative sign of the exponentially small term across a Stokes line. It does **not** fix the absolute sign or magnitude of the physical condensate, because that is the transseries parameter $\sigma$, selected by the vacuum. In the campaign's language: resurgence can see the residue's location and modulus, not the sign of $f_1^2$. The strong-coupling theorem P3 ($f_1^2>0$ at $\beta<1/48$) is precisely a nonperturbative input of the type resurgence leaves free.

### 1.5 Which methods are blind

| Method | Sees $e^{-1/(2\beta_0g^2)}$? | Can certify $\Delta>0$? | Why |
|---|---|---|---|
| Finite-order perturbation theory | **No** | No | flat function; Taylor series has no exponential |
| Borel–Padé resummation | Partly (singularity location) | No | residue/normalisation is an input or a fit |
| OPE / SVZ sum rules | **No** | No | condensates are inputs; gap not derived |
| DSE truncations (rainbow-ladder etc.) | Yes (a branch) | No | truncation and vertex ansatz uncontrolled; IR mass gauge-dependent |
| FRG derivative/vertex truncations | Yes | No | finite truncation cannot certify an essential singularity; background-field dependence |
| Lattice strong-coupling expansion | **No** (analytic in $\beta$) | Only inside its radius | radius ends at the bulk transition; cannot cross to continuum |
| Lattice continuum extrapolation | Yes (nonperturbative definition) | Not yet | continuum limit is the open step |
| Semiclassical Savvidy expansion | Complex (negative mode) | No | real part needs thimbles; saddle unstable |
| Resurgence / transseries | Yes (exponent) | No | $\sigma$ free; unstable-saddle real part not fixed |
| Conformal bootstrap without a mass scale | **No** | No | exact scale invariance forbids a gap |
| Holographic models | Yes (by construction) | No | gap is an IR boundary condition, not derived |
| Rigorous RG (Balaban program) | Yes | In principle | not completed in 4D |

The table is the operational meaning of "analytic in $g^2$ ⇒ blind": any method whose output is a formal power series in $g^2$ around $0$, or a finite truncation of one, cannot contain the gap.

### 1.6 The finite-volume caveat (important for v3)

At fixed volume $L^4$ and cutoff $a$, the transfer matrix is finite-dimensional and $\Delta_L(g^2)$ is **analytic** in $g^2$. The essential singularity appears only in the infinite-volume/continuum limit. Therefore:

- "$\Delta_L>0$ for every $L$" is **not** the mass gap;
- a proof needs a **uniform-in-$L$** lower bound, $\Delta_L\ge\delta>0$, before any limit;
- the campaign's $L$-uniform bounds (P3, Lemma SC, $W_{\rm geo}\ge0$) are the correct logical shape; the open step is exactly the $L$-uniform control at the physical coupling.

This is the sharpest way to state why the problem is hard: one must prove that a limit of analytic functions does not collapse.

---

## 2. Annotated bibliography (33 entries)

### A. Resurgence and transseries

**R1 — Aniceto, Başar & Schiappa (2019).** *A primer on resurgent transseries and their asymptotics.* Phys. Rept. **809**, 1–135. arXiv:1802.10441 [hep-th]. DOI:10.1016/j.physrep.2019.02.003. **[O][V]**
- Content: 192-page pedagogical primer; transseries, Stokes lines, Borel transforms, Écalle alien calculus; a linear example with finite instanton sectors and a nonlinear one with infinite sectors; "statistical-mechanical" and Virasoro-like reorganisations.
- Use in v3: canonical definitions and the language for "what is fixed vs free"; no YM-specific theorem.

**R2 — Dorigoni (2019).** *An introduction to resurgence, trans-series and alien calculus.* Annals Phys. **409**, 167914. arXiv:1411.3585 [hep-th]. DOI:10.1016/j.aop.2019.167914. **[O][V]**
- Content: 48-page minimal introduction; resurgence, transseries, alien derivatives; worked examples.
- Use in v3: shorter companion to R1 for definitions; same scope, fewer examples.

**R3 — Mariño (2014, adjacent).** *Lectures on non-perturbative effects in large N gauge theories, matrix models and strings.* Fortsch. Phys. **62**, 455–540. arXiv:1206.6272 [hep-th]. DOI:10.1002/prop.201400005. **[O][V]**
- Content: resurgence of the $1/N$ expansion; matrix models and topological strings; instanton actions and Stokes data in computable settings.
- Use in v3: the benchmark for what "resurgence determines" means when the model is exactly solvable; not 4D YM.

**R4 — Argyres & Ünsal (2012, adjacent).** *The semi-classical expansion and resurgence in gauge theories: new perturbative, instanton, bion, and renormalon effects.* JHEP **08** (2012) 063. arXiv:1206.1890 [hep-th]. DOI:10.1007/JHEP08(2012)063. **[C][V]**
- Content: controlled semiclassical expansion on $\mathbb{R}^3\times S^1$; monopole-instantons, bions, renormalons; the transseries and its continuity as the circle shrinks; competition of instanton and renormalon sectors.
- Use in v3: the best-controlled gauge-theory realisation of a transseries; shows the exponential sector is fixed on a small circle. The $\mathbb{R}^4$ limit is the open step.

**R5 — Ünsal (2008, adjacent).** *Abelian duality, confinement, and chiral symmetry breaking in QCD(adj).* Phys. Rev. Lett. **100**, 032005. arXiv:0708.1772 [hep-th]. DOI:10.1103/PhysRevLett.100.032005. **[C][V]**
- Content: abelian duality on small $S^1$; confinement and chiral symmetry breaking from monopole/bion dynamics in QCD(adj).
- Use in v3: the physical mechanism behind R4; a concrete "gap from semiclassics" in a controlled regime.

**R6 — Costin (2009).** *Asymptotics and Borel Summability.* Chapman & Hall/CRC, Boca Raton. **[T][U]**
- Content: rigorous theory of Borel summability, resurgence, and transseries for nonlinear ODEs and analyzable functions; conditions under which a divergent series is Borel-summable.
- Use in v3: the rigorous standard against which any claim of "Borel summability of YM" must be judged; not a YM result.

### B. Renormalons, Borel summability, and the gluon condensate

**R7 — Beneke (1999, adjacent).** *Renormalons.* Phys. Rept. **317**, 1–142. arXiv:hep-ph/9807443. DOI:10.1016/S0370-1573(98)00130-6. **[O][V]**
- Content: review of renormalons; Borel singularities at $t=u/\beta_0$; IR renormalons and power corrections; the gluon-condensate renormalon; the failure of naive Borel summation.
- Use in v3: the precise dictionary between Borel singularities and nonperturbative corrections; the source of the $[NG]$ statement in §1.3.

**R8 — Gubler & Satow (2019).** *Recent progress in QCD condensate evaluations and sum rules.* Prog. Part. Nucl. Phys. **106**, 1–67. arXiv:1812.00385 [hep-ph]. DOI:10.1016/j.ppnp.2019.02.005. **[O][V]**
- Content: review of QCD vacuum condensates (quark, gluon, mixed, higher-dimensional); sum-rule, lattice, and instanton-model determinations; in-medium condensates.
- Use in v3: the status of the gluon condensate as an input with scheme dependence; the gap is not derived from it.

**R9 — Shifman, Vainshtein & Zakharov (1979, foundational).** *QCD and resonance physics. Theoretical foundations.* Nucl. Phys. **B147**, 385–447. **[O][U]**
- Content: original QCD sum rules; OPE plus condensates; the gluon condensate.
- Use in v3: the historical origin of the condensate input; the logic "condensate in, hadron masses out," not "gap out."

**R10 — Lipatov (1977, foundational).** *Divergence of the perturbation-theory series and the quasi-classical theory.* Sov. Phys. JETP **45**, 216–223 [Zh. Eksp. Teor. Fiz. **72**, 411–427]. **[T][U]**
- Content: factorial growth of perturbative coefficients from the instanton saddle; large-order asymptotics.
- Use in v3: the mechanism behind Borel singularities; the $n!$ growth that resurgence inverts.

**R11 — 't Hooft (1979, foundational).** *Can we solve QCD?* In *The Whys of Subnuclear Physics* (Erice 1977), ed. A. Zichichi, Plenum, New York. **[NG][U]**
- Content: the renormalon argument; the perturbative series is not Borel-summable; IR renormalons produce ambiguities of the order of nonperturbative corrections.
- Use in v3: the original no-go for naive Borel summation; the conceptual ancestor of §1.3.

### C. Functional renormalisation group

**R12 — Dupuis, Canet, Eichhorn, Metzner, Pawlowski, Tissier & Wschebor (2021).** *The nonperturbative functional renormalization group and its applications.* Phys. Rept. **910**, 1–114. arXiv:2006.04853 [hep-th]. DOI:10.1016/j.physrep.2021.01.001. **[O][V]**
- Content: 93-page review; Wetterich equation; derivative and vertex expansions; gauge theories; critical phenomena; error control.
- Use in v3: the era's canonical FRG reference; the truncation hierarchy and its limits.

**R13 — Reuter & Saueressig (2019; online 2018).** *Quantum Gravity and the Functional Renormalization Group: The Road towards Asymptotic Safety.* Cambridge University Press. DOI:10.1017/9781316227596. **[O][V]**
- Content: book-length treatment of the functional RG; asymptotic safety; gauge/gravity applications.
- Use in v3: the FRG formalism and its background-field dependence; the "exact equation, approximate solutions" structure.

**R14 — Reuter (1998).** *Nonperturbative evolution equation for quantum gravity.* Phys. Rev. D **57**, 971–985. arXiv:hep-th/9605030. DOI:10.1103/PhysRevD.57.971. **[T][V]**
- Content: the Wetterich equation applied to gravity; the first nonperturbative flow.
- Use in v3: foundational; shows both the exactness of the flow equation and the truncation problem.

**R15 — Gies (2012).** *Introduction to the functional RG and applications to gauge theories.* Lecture Notes in Physics **852**, 287–348, Springer. arXiv:hep-ph/0611146. DOI:10.1007/978-3-642-27320-9_6. **[O][V]**
- Content: lecture notes; FRG for gauge theories; background field; Nielsen identities; flow equations.
- Use in v3: the gauge-theory entry point to the FRG; the treatment of gauge invariance in a truncation.

### D. Dyson–Schwinger equations and the scaling/decoupling resolution

**R16 — Huber (2020).** *Nonperturbative properties of Yang–Mills theories.* Phys. Rept. **879**, 1–92. arXiv:1808.05227 [hep-ph]. DOI:10.1016/j.physrep.2020.04.004. **[O][V]**
- Content: 135-page review of gluon correlation functions; Landau gauge; DSE; 2D/3D/4D; truncation systematics; comparison of scaling and decoupling; linear covariant, Coulomb, and maximally Abelian gauges.
- Use in v3: the era's DSE reference; the honest statement that the IR mass is truncation- and gauge-dependent.

**R17 — Alkofer & von Smekal (2001, foundational).** *The infrared behavior of QCD Green's functions — confinement, dynamical symmetry breaking, and hadrons as relativistic bound states.* Phys. Rept. **353**, 281. arXiv:hep-ph/0007355. DOI:10.1016/S0370-1573(01)00010-2. **[O][V]**
- Content: the scaling solution; ghost enhancement, gluon suppression; confinement and dynamical symmetry breaking.
- Use in v3: the scaling-solution canon; the baseline the decoupling solution reacts to.

**R18 — Fischer (2006).** *Infrared properties of QCD from Dyson–Schwinger equations.* J. Phys. G **32**, R253–R291. arXiv:hep-ph/0605173. DOI:10.1088/0954-3899/32/8/R02. **[O][V]**
- Content: DSE review; IR properties; confinement; dynamical chiral symmetry breaking; vertex ansätze.
- Use in v3: the DSE review paired with R16; the truncation and vertex-ansatz problem stated explicitly.

**R19 — Fischer, Maas & Pawlowski (2009).** *On the infrared behavior of Landau gauge Yang–Mills theory.* Annals Phys. **324**, 2408–2437. arXiv:0810.1987 [hep-ph]. DOI:10.1016/j.aop.2009.07.009. **[O][V]**
- Content: sharp statement of the scaling solution in 4D; ghost dominance; the ghost–gluon vertex; the IR fixed point.
- Use in v3: the cleanest scaling-solution reference; the target of the decoupling comparison.

**R20 — Cyrol, Mitter, Pawlowski & Strodthoff (2018).** *Non-perturbative quark, gluon and meson correlators of unquenched QCD.* Phys. Rev. D **97**, 054006. arXiv:1706.06326 [hep-ph]. DOI:10.1103/PhysRevD.97.054006. **[O][V]**
- Content: FRG computation of unquenched QCD correlators; quantitative gluon, quark, and meson propagators; comparison with lattice.
- Use in v3: the era's quantitative FRG benchmark; truncation systematics in a realistic setting.

**R21 — von Smekal, Hauck & Alkofer (1997, foundational).** *The infrared behavior of gluon and ghost propagators in Landau gauge QCD.* Phys. Rev. Lett. **79**, 3591–3594. arXiv:hep-ph/9705242. DOI:10.1103/PhysRevLett.79.3591. **[O][V]**
- Content: the first scaling solution; ghost propagator enhanced, gluon suppressed in the IR.
- Use in v3: the origin of the scaling scenario.

**R22 — Aguilar, Binosi & Papavassiliou (2008).** *Gluon and ghost propagators in the Landau gauge: deriving lattice results from Schwinger–Dyson equations.* Phys. Rev. D **78**, 025010. arXiv:0802.1870 [hep-ph]. DOI:10.1103/PhysRevD.78.025010. **[O][V]**
- Content: the decoupling (massive-gluon) solution from DSE; a dynamical gluon mass scale; comparison with lattice.
- Use in v3: the decoupling scenario; the gauge-fixed IR mass, explicitly not the physical gap.

### E. Gribov–Zwanziger and the Gribov horizon

**R23 — Gribov (1978, foundational).** *Quantization of non-Abelian gauge theories.* Nucl. Phys. **B139**, 1–19. **[T][U]**
- Content: Gribov copies; the Gribov horizon; the Gribov action; the ambiguity of gauge fixing beyond perturbation theory.
- Use in v3: the origin of the horizon; a gauge-fixed construction, not a physical Hilbert space.

**R24 — Zwanziger (1989).** *Local and renormalizable action from the Gribov horizon.* Nucl. Phys. **B323**, 513–544. DOI:10.1016/0550-3213(89)90122-3. **[T][V]**
- Content: the local, renormalisable Gribov–Zwanziger action; auxiliary fields; the horizon condition as a gap equation for the Gribov parameter.
- Use in v3: the GZ framework in its modern form; the gap equation is not the physical mass gap.

**R25 — Vandersickel & Zwanziger (2012, adjacent).** *The Gribov problem and QCD dynamics.* Phys. Rept. **520**, 175–251. arXiv:1202.1491 [hep-th]. DOI:10.1016/j.physrep.2012.07.003. **[O][V]**
- Content: 121-page review of the Gribov problem and QCD dynamics; GZ; confinement; lattice evidence; positivity violation.
- Use in v3: the era's synthesis; the explicit warning that the GZ gluon propagator is not a physical particle.

**R26 — Dell'Antonio & Zwanziger (1991).** *Every gauge orbit passes inside the Gribov horizon.* Commun. Math. Phys. **138**, 291–305. **[T][U]**
- Content: rigorous proof that every gauge orbit intersects the Gribov region.
- Use in v3: the rigorous backbone of the GZ picture; the horizon is unavoidable, not an artefact.

### F. Lattice glueballs and the Wilson flow

**R27 — Lüscher (2010).** *Properties and uses of the Wilson flow in lattice QCD.* JHEP **08** (2010) 071. arXiv:1006.4518 [hep-lat]. DOI:10.1007/JHEP08(2010)071. **[T][V]**
- Content: the Wilson flow; the scale $t_0$; smoothness and renormalisation properties; scale setting.
- Use in v3: the lattice scale-setting tool; the flow coupling as a candidate nonperturbative definition.

**R28 — Morningstar & Peardon (1999, foundational).** *The glueball spectrum from an anisotropic lattice study.* Phys. Rev. D **60**, 034509. arXiv:hep-lat/9901004. DOI:10.1103/PhysRevD.60.034509. **[O][V]**
- Content: glueball spectrum on an anisotropic lattice; $0^{++}$, $2^{++}$, and excited states; precision at the time.
- Use in v3: the pre-era benchmark; the physical scale of the gap.

**R29 — Athenodorou & Teper (2020).** *The glueball spectrum of SU(3) gauge theory in 3+1 dimension.* JHEP **11** (2020) 172. arXiv:2007.06422 [hep-lat]. DOI:10.1007/JHEP11(2020)172. **[O][V]**
- Content: continuum extrapolation of the SU(3) glueball spectrum; precise $0^{++}$ and $2^{++}$ masses; $m_{0^{++}}/\sqrt{\sigma}\approx3.4$.
- Use in v3: the era's lattice benchmark; any proof must reproduce this number.

**R30 — Bali, Bauer & Pineda (2014, adjacent).** *Model independent determination of the gluon condensate in four-dimensional SU(3) gauge theory.* Phys. Rev. Lett. **112**, 042001. **[O][U]**
- Content: a model-independent lattice extraction of the gluon condensate.
- Use in v3: the condensate as a measured input with scheme dependence; not a derivation of the gap.

### G. Gauge-theory bootstrap and sign-problem methods

**R31 — Anderson & Kruczenski (2017, adjacent).** *Loop equations and bootstrap methods in the lattice.* Nucl. Phys. **B921**, 702–726. arXiv:1612.08140 [hep-th]. DOI:10.1016/j.nuclphysb.2017.06.009. **[O][V]**
- Content: loop equations plus positivity/bootstrap on the lattice; truncated moment problem; proof of principle.
- Use in v3: the first gauge-theory bootstrap; the truncation/positivity problem in a finite-dimensional setting.

**R32 — Han, Hartnoll & Kruthoff (2020).** *Bootstrapping matrix quantum mechanics.* Phys. Rev. Lett. **125**, 041601. arXiv:2004.10212 [hep-th]. DOI:10.1103/PhysRevLett.125.041601. **[O][V]**
- Content: bootstrap for matrix quantum mechanics; gauge-theory-adjacent; positivity constraints on the spectrum.
- Use in v3: the bootstrap method in a controlled setting; not 4D YM.

**R33 — Alexandru, Başar, Bedaque & Warrington (2022; reviews 2018–21).** *Complex paths around the sign problem.* Rev. Mod. Phys. **94**, 015006. arXiv:2007.05436 [hep-lat]. DOI:10.1103/RevModPhys.94.015006. **[O][V]**
- Content: review of Lefschetz thimbles, complex Langevin, and complex paths for the sign problem.
- Use in v3: the toolkit for complex saddles (including Savvidy); the honest statement of when the real part is under control.

---

## 3. What this era settles, and what it does not

**Settles / sharpens**

- The vocabulary: transseries, Stokes lines, Borel singularities, renormalons, alien calculus (R1–R7).
- The quantitative FRG and DSE machinery, with explicit truncation hierarchies (R12–R22).
- The GZ framework and the rigorous horizon statement (R23–R26).
- The lattice benchmark: $m_{0^{++}}/\sqrt{\sigma}\approx3.4$, continuum extrapolated (R27–R29).
- The phenomenological resolution of scaling vs decoupling: lattice favours decoupling; the scaling solution survives as a DSE branch. Crucially, **both are gauge-fixed statements and neither is the physical mass gap**.
- The precise location of the obstruction: essential singularity at $g^2=0$; positive-axis Borel singularities; the transseries parameter is free.

**Does not settle**

- No rigorous construction of 4D pure YM; no OS/Wightman verification.
- No proof of Borel summability or non-summability for 4D YM.
- No derivation of $\Delta/\Lambda$; no controlled continuum limit of the lattice gap.
- No truncation-error bound for FRG/DSE; no gauge-invariant IR mass.
- No real Lefschetz-thimble computation of the Savvidy sector in 4D.
- No 4D gauge bootstrap; the 2022–23 lattice YM bootstrap is post-era and truncation-bound.
- No resolution of the sign of the condensate / transseries parameter from first principles.

**Cross-links to the campaign's objects**

- $f_1^2$ (P3 residue) is the strong-coupling analogue of the transseries parameter $\sigma$: a nonperturbative datum that perturbation theory cannot fix. P3 fixes its sign at $\beta<1/48$; the open step is transporting it to $\beta\approx2.5$.
- $G''=S-D$ is an exact identity; the defect $D/S\approx2\%$ at $\beta=2.5$ is the quantitative shadow of the exponentially small sector. Resurgence says the exponential scale is $e^{-1/(2\beta_0g^2)}$; it does not say $D/S$ is small.
- $W_{\rm geo}\ge0$ is a theorem about the Haar measure; it is not a statement about the Borel plane. The two must be combined, not conflated.
- Savvidy concavity is the negative mode that makes the semiclassical transseries complex; the thimble toolkit (R33) is the only era tool addressed to it.

---

## 4. Position in 14 eras

*Campaign working index; E12 is this survey. Boundaries are thematic, not strict publication years.*

| Era | Window | One-line characterisation |
|---|---|---|
| **E1** | 1954–1972 | Non-Abelian gauge theory and its quantisation (Faddeev–Popov); the problem is posed, no confinement mechanism. |
| **E2** | 1973–1976 | Asymptotic freedom and the birth of QCD; Gribov copies; 't Hooft–Veltman renormalisation. |
| **E3** | 1974–1979 | Wilson lattice gauge theory; strong-coupling confinement; Osterwalder–Seiler; first rigorous mass-gap-like bounds. |
| **E4** | 1975–1981 | Instantons, $\theta$-vacua, $U(1)_A$; SVZ sum rules and the gluon condensate; the semiclassical vacuum. |
| **E5** | 1977–1986 | Confinement pictures: dual superconductivity, center vortices, monopoles; the Gribov horizon; DSE beginnings. |
| **E6** | 1982–1995 | Continuum nonperturbative QCD: DSE/rainbow-ladder, Cornwall, the GZ action; lattice glueball spectroscopy starts. |
| **E7** | 1994–2000 | Large-$N$, AdS/CFT, holographic confinement; Seiberg–Witten exact SUSY results. |
| **E8** | 1998–2008 | Precision lattice glueball spectrum; improved actions; DSE infrared scaling solutions; the GZ confinement criterion. |
| **E9** | 2005–2012 | Functional RG and DSE infrared: Gies, Pawlowski, Fischer; scaling vs decoupling; the gauge-fixed IR mass. |
| **E10** | 2008–2014 | Resurgence enters QFT: Ünsal, Argyres–Ünsal; bions and renormalons on $\mathbb{R}^3\times S^1$; transseries in gauge theory. |
| **E11** | 2010–2017 | Gradient flow (Lüscher); lattice scale setting and precision; chiral lattice; QCD thermodynamics. |
| **E12** | **2018–2021** | **Resurgence/transseries reviews; FRG/DSE reviews; GZ synthesis; lattice glueball precision; gauge-bootstrap attempts. This survey.** |
| **E13** | 2022–2024 | Bootstrap expansion; complex Langevin/Lefschetz review; rigorous stochastic quantisation in lower dimensions; continuum-limit pushes. |
| **E14** | 2024–2026 | AI-assisted mathematics and Lean formalisation of YM; conditional-theorem campaigns; the present program. |

---

## 5. Gaps and attack surfaces for v3

1. **No theorem on Borel summability of 4D YM.** The renormalon picture is heuristic. *Attack:* prove a precise statement for a stripped-down model with the same positive-axis singularity structure, or prove Borel summability of a gauge-invariant observable in a confining toy model; then ask what fails in 4D.
2. **The transseries parameter $\sigma$ is never computed in 4D.** *Attack:* define a lattice analogue of $\sigma$ (e.g. the condensate normalisation extracted from a controlled operator), take a continuum limit in a simplified observable, and compare with the resurgence prediction for the exponent.
3. **The Savvidy negative mode makes the semiclassical sector complex.** *Attack:* a real Lefschetz-thimble computation of the one-instanton/bion contribution on small 4D lattices; compare the real part with the resurgence-implied ambiguity. This is the only era tool aimed at the campaign's actual wall.
4. **FRG/DSE truncations have no error bar.** *Attack:* build a certified truncation-error functional (interval arithmetic on the vertex expansion) for the IR mass; or prove a conditional statement: "if the truncation error is below $\epsilon$, then $\Delta>0$."
5. **Scaling vs decoupling is not a theorem.** *Attack:* prove uniqueness of the IR fixed point of the Landau-gauge DSE system given the ghost–gluon vertex, or measure the ghost form factor at very large volumes with controlled systematics.
6. **GZ has no rigorous 4D measure.** *Attack:* construct the GZ measure rigorously in 3D or on a finite lattice with reflection positivity; decide whether the GZ transfer matrix has a positive spectrum.
7. **No 4D gauge bootstrap with a mass gap.** *Attack:* formulate positivity constraints from OS/reflection positivity on Wilson-loop moments; produce a rigorous truncated bootstrap with error bars, even in 2D or 3D first.
8. **Lattice continuum extrapolation is the bottleneck.** *Attack:* multi-$a$, multi-$L$, flow-based scale setting, blinded analyses; quote a continuum $m_{0^{++}}/\sqrt{\sigma}$ with a defensible error.
9. **The Wilson flow lacks a proved continuum-limit coupling.** *Attack:* define a flow-based renormalised coupling with a proved continuum limit and a proved relation to $\Lambda$; this would give a nonperturbative definition of the scale in §1.1.
10. **State the missing step in the right logical form.** Because of §1.6, the target is a **uniform-in-$L$ lower bound** on a nonperturbative quantity (the residue $f_1^2$, the GZ horizon parameter, or a gauge-invariant IR mass), not a finite-order perturbative computation. v3 should make this the explicit "exact missing step," with no hidden slack.

---

## 6. Verification ledger

| # | Entry (short) | Identifier checked | Check method | Status | Notes |
|---|---|---|---|---|---|
| R1 | Aniceto–Başar–Schiappa | arXiv:1802.10441; Phys. Rept. 809 | WebSearch + arXiv | **VERIFIED** | 809 (2019) 1–135; DOI 10.1016/j.physrep.2019.02.003 |
| R2 | Dorigoni | arXiv:1411.3585 | arXiv + Crossref | **VERIFIED** | Annals Phys. 409 (2019) 167914 |
| R3 | Mariño | arXiv:1206.6272 | arXiv + Crossref | **VERIFIED** | Fortsch. Phys. 62 (2014) 455–540 |
| R4 | Argyres–Ünsal | arXiv:1206.1890 | arXiv | **VERIFIED** | JHEP 08 (2012) 063 via DOI 10.1007/JHEP08(2012)063 |
| R5 | Ünsal | arXiv:0708.1772 | arXiv | **VERIFIED** | PRL 100 (2008) 032005 |
| R6 | Costin, *Asymptotics and Borel Summability* | book | none | **UNVERIFIED** | Standard CRC 2009 record; no DOI check |
| R7 | Beneke, *Renormalons* | arXiv:hep-ph/9807443 | arXiv | **VERIFIED** | Phys. Rept. 317 (1999) 1–142 |
| R8 | Gubler–Satow | arXiv:1812.00385 | arXiv | **VERIFIED** | Prog. Part. Nucl. Phys. 106 (2019) 1–67; DOI 10.1016/j.ppnp.2019.02.005 |
| R9 | Shifman–Vainshtein–Zakharov | — | none | **UNVERIFIED** | Nucl. Phys. B147 (1979) 385–447, standard record |
| R10 | Lipatov | — | none | **UNVERIFIED** | Sov. Phys. JETP 45 (1977) 216, standard record |
| R11 | 't Hooft, *Can we solve QCD?* | — | none | **UNVERIFIED** | Erice 1977 proceedings, Plenum 1979; page numbers omitted |
| R12 | Dupuis et al. | arXiv:2006.04853 | arXiv | **VERIFIED** | Phys. Rept. 910 (2021) 1–114; DOI 10.1016/j.physrep.2021.01.001 |
| R13 | Reuter–Saueressig | DOI:10.1017/9781316227596 | Crossref | **VERIFIED** | CUP; online 2018, print 2019; ISBN 978-1-107-10732-8 |
| R14 | Reuter 1998 | arXiv:hep-th/9605030 | arXiv | **VERIFIED** | PRD 57 (1998) 971–985 |
| R15 | Gies | arXiv:hep-ph/0611146; DOI:10.1007/978-3-642-27320-9_6 | arXiv + Crossref | **VERIFIED** | LNP 852 (2012) 287–348; volume number from standard record, Crossref omits it |
| R16 | Huber | arXiv:1808.05227 | WebSearch + arXiv | **VERIFIED** | Phys. Rept. 879 (2020) 1–92; DOI 10.1016/j.physrep.2020.04.004 |
| R17 | Alkofer–von Smekal | arXiv:hep-ph/0007355 | arXiv | **VERIFIED** | Phys. Rept. 353 (2001) 281; end page not shown on arXiv, omitted |
| R18 | Fischer 2006 | arXiv:hep-ph/0605173 | arXiv | **VERIFIED** | J. Phys. G32 (2006) R253–R291 |
| R19 | Fischer–Maas–Pawlowski | arXiv:0810.1987 | arXiv | **VERIFIED** | Annals Phys. 324 (2009) 2408–2437 |
| R20 | Cyrol–Mitter–Pawlowski–Strodthoff | arXiv:1706.06326 | arXiv | **VERIFIED** | PRD 97 (2018) 054006 |
| R21 | von Smekal–Hauck–Alkofer | arXiv:hep-ph/9705242 | arXiv | **VERIFIED** | PRL 79 (1997) 3591–3594 |
| R22 | Aguilar–Binosi–Papavassiliou | arXiv:0802.1870 | arXiv | **VERIFIED** | PRD 78 (2008) 025010 |
| R23 | Gribov 1978 | — | none | **UNVERIFIED** | Nucl. Phys. B139 (1978) 1–19, standard record |
| R24 | Zwanziger 1989 | DOI:10.1016/0550-3213(89)90122-3 | Crossref | **VERIFIED** | Nucl. Phys. B323 (1989) 513–544 |
| R25 | Vandersickel–Zwanziger | arXiv:1202.1491; DOI | arXiv + Crossref | **VERIFIED** | Phys. Rept. 520 (2012) 175–251 |
| R26 | Dell'Antonio–Zwanziger | — | none | **UNVERIFIED** | Commun. Math. Phys. 138 (1991) 291–305, standard record |
| R27 | Lüscher 2010 | arXiv:1006.4518 | arXiv | **VERIFIED** | JHEP 08 (2010) 071 |
| R28 | Morningstar–Peardon | arXiv:hep-lat/9901004 | arXiv | **VERIFIED** | PRD 60 (1999) 034509 |
| R29 | Athenodorou–Teper | arXiv:2007.06422 | arXiv | **VERIFIED** | JHEP 11 (2020) 172 |
| R30 | Bali–Bauer–Pineda | — | none | **UNVERIFIED** | PRL 112 (2014) 042001, standard record |
| R31 | Anderson–Kruczenski | arXiv:1612.08140 | arXiv + Crossref | **VERIFIED** | Nucl. Phys. B921 (2017) 702–726 |
| R32 | Han–Hartnoll–Kruthoff | arXiv:2004.10212 | arXiv | **VERIFIED** | PRL 125 (2020) 041601 |
| R33 | Alexandru–Başar–Bedaque–Warrington | arXiv:2007.05436; DOI | arXiv + Crossref | **VERIFIED** | RMP 94 (2022) 015006; published after the era, reviews 2018–21 work |

**Tally:** 33 entries; **26 VERIFIED**, **7 UNVERIFIED**. The seven unverified entries are foundational pre-2018 items cited from the standard record (R6, R9, R10, R11, R23, R26, R30). No unverified entry is load-bearing for a quantitative claim in §1; they are historical anchors. All quantitative claims in §1 rest on R1–R5, R7–R8, R12–R22, R24–R25, R27–R29, R31–R33.

---

## 7. One-paragraph handoff

E12 is the era that made the obstruction precise without removing it. Its durable contribution is the transseries/Borel language: the mass gap is an essential singularity at $g^2=0$, the perturbative series is its flat Taylor shadow, and renormalons place the ambiguity on the positive Borel axis. Resurgence fixes the exponent and the discontinuity but leaves the transseries parameter — the number $C$ in $\Delta=C\Lambda$ — free. The era's quantitative tools (FRG, DSE, GZ, lattice flow, glueball spectroscopy) are all gauge-fixed or truncation-bound, and the scaling-vs-decoupling dispute is a dispute about gauge-fixed correlators, not about the physical spectrum. For v3 the usable inheritance is: a precise statement of the missing step (a uniform-in-$L$ lower bound on a nonperturbative datum), a benchmark number ($m_{0^{++}}/\sqrt{\sigma}\approx3.4$), and a warning — any method analytic in $g^2$ is structurally blind, so the proof must be nonperturbative at the level of the measure, not at the level of the series.
