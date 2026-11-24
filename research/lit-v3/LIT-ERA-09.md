# LIT-ERA-09 — 2002–2010: lattice precision, step scaling, and the conformal window

> V3 literature wave, era 09. Scope: the assigned 2002–2010 window plus the precursor papers the
> assignment names (Schrödinger functional 1992, step scaling 1994, $r_0$ 1994, two-loop matching
> 1995, glueballs 1999) and the two 2013 LatKMI papers needed to close the many-flavour story.
> Labels: `[T]` rigorous theorem (inside explicit hypotheses); `[C]` conditional theorem /
> controlled approximation; `[NG]` no-go for a named class; `[O]` exact open/missing step, or a
> method that defines the object without closing it; `[V]` numerical evidence (lattice), with
> continuum/artifact status stated. The label is the strongest statement the entry contributes to
> the mass-gap problem, not the paper's genre.
>
> Verification method: `VERIFIED` means the bibliographic data was checked during this pass against
> the arXiv abstract page or a publisher/ADS/CERN record. The shared WebSearch budget was exhausted
> mid-task; remaining items were either checked by direct arXiv fetch or left `UNVERIFIED` rather
> than guessed. No bibliographic field below is invented.

---

## 0. Verdict

This era turned the lattice from a qualitative confinement demonstration into a **continuum-controlled
quantitative instrument** — but only in the ultraviolet. Step scaling in the Schrödinger functional
(SF) gives a non-perturbative running coupling whose continuum limit is controlled to a few percent
in pure SU(3) YM; the static potential and $r_0$ are continuum-extrapolated; the lightest glueball
mass is continuum-extrapolated to $\approx 1.7$ GeV; the Wilson flow gives a rigorous smoothing
map and a renormalized energy density. None of this bounds the mass gap. In parallel, the
Banks–Zaks conformal window was revived as the correct language for the **gapless alternative**:
a gauge theory can flow to an interacting IR CFT instead of confining. The lattice many-flavour
programme produced strong but finite-$(a,L,m)$ evidence that $N_f=12$ is (near-)conformal and
$N_f=8$ is walking/near-conformal; the $N_f=12$ classification remains contested. The logically
decisive statement — pure YM ($N_f=0$) is outside the window — is **perturbatively** rigorous only:
the two-loop $\beta$-function has no positive IR zero. A non-perturbative proof would require either
a construction plus a gap proof, or a rigorous no-IRFP theorem, or a rigorous monotonicity in $N_f$.
That is precisely the missing step.

---

## 1. Cluster A — Schrödinger functional, step scaling, running coupling (8 entries)

**[E09-01] [O] The Schrödinger functional is defined.**
M. Lüscher, R. Narayanan, P. Weisz, U. Wolff, “The Schrödinger functional — a renormalizable probe
for non-Abelian gauge theories,” *Nucl. Phys.* **B384** (1992) 168–228.
DOI: 10.1016/0550-3213(92)90466-O. arXiv:hep-lat/9207009. **VERIFIED** (arXiv abs).
Establishes: a finite-volume gauge theory with Dirichlet boundary conditions in time and periodic
space; the boundary values act as external sources; the effective action is renormalizable and its
derivatives define a running coupling at scale $\mu=1/L$. Continuum status: the SF is a
continuum-defined scheme; the lattice is an approximation. It does **not** address the mass gap;
its finite-volume spectrum has a gap $\sim 1/L$ by construction.

**[E09-02] [V] Step scaling determines the pure-YM running coupling.**
M. Lüscher, R. Sommer, P. Weisz, U. Wolff, “A precise determination of the running coupling in the
SU(3) Yang–Mills theory,” *Nucl. Phys.* **B413** (1994) 481–502.
DOI: 10.1016/0550-3213(94)90629-7. arXiv:hep-lat/9309005. **VERIFIED** (arXiv abs + publisher).
Establishes: the step-scaling method — evolve the SF coupling by a fixed factor $s$ in finite volume,
take the continuum limit at each step, and iterate from small to large scales; the low-energy result
is connected to $\sqrt{K}$ (string tension) and the high-energy result to perturbation theory; the
$\overline{\rm MS}$ coupling is obtained with a few-percent total error. Continuum status:
continuum-controlled in the UV; $O(a^2)$ artefacts removed by extrapolation. This is the canonical
“running coupling in pure YM” result. It says nothing about the IR residue $f_1^2$.

**[E09-03] [V] The $r_0$ scale from the static force.**
R. Sommer, “A new way to set the energy scale in lattice gauge theories and its application to the
static force and $\alpha_s$ in SU(2) Yang–Mills theory,” *Nucl. Phys.* **B411** (1994) 839–854.
DOI: 10.1016/0550-3213(94)90473-1. arXiv:hep-lat/9310022. **VERIFIED** (arXiv abs).
Establishes: $r_0$ from the dimensionless combination $r^2 F(r)$ at the inflection point of the
static force; $r_0\approx 0.5$ fm; a hadronic scale that is continuum-extrapolable and insensitive
to the far IR. Continuum status: continuum-controlled; used throughout this era as the scale that
converts lattice numbers to physical units. It is a scale-setting device, not a gap statement.

**[E09-04] [T] Two-loop matching between the bare lattice coupling and $\overline{\rm MS}$.**
M. Lüscher, P. Weisz, “Computation of the relation between the bare lattice coupling and the
$\overline{\rm MS}$ coupling in SU(N) gauge theories to two loops,” *Nucl. Phys.* **B452** (1995)
234–260. DOI: 10.1016/0550-3213(95)00338-S. arXiv:hep-lat/9505011. **VERIFIED** (arXiv abs).
Establishes: the two-loop conversion factor, which is what turns the step-scaling result into
$\Lambda_{\overline{\rm MS}}$. Label `[T]` is for the perturbative computation inside BPHZ; it is
not a non-perturbative theorem. Continuum status: continuum perturbation theory; lattice artefacts
enter only through the bare-coupling definition.

**[E09-05] [T] One-loop SF running coupling.**
S. Sint, R. Sommer, “The running coupling from the QCD Schrödinger functional: a one-loop analysis,”
*Nucl. Phys.* **B465** (1996) 71–98. DOI: 10.1016/0550-3213(96)00020-X. arXiv:hep-lat/9508012.
**VERIFIED** (arXiv abs).
Establishes: the perturbative definition and one-loop $\beta$-function of the SF coupling, including
the boundary-counterterm tuning. Continuum status: continuum perturbation theory. It fixes the
scheme used by E09-02; without it the step-scaling result has no perturbative anchor.

**[E09-06] [V] Non-perturbative quark-mass renormalization in the quenched theory.**
S. Capitani, M. Lüscher, R. Sommer, H. Wittig, “Non-perturbative quark mass renormalization in
quenched lattice QCD,” *Nucl. Phys.* **B544** (1999) 669–698.
DOI: 10.1016/S0550-3213(98)00857-8. arXiv:hep-lat/9810063. **VERIFIED** (arXiv abs).
Establishes: SF-based non-perturbative renormalization of the quark mass in the quenched
approximation ($N_f=0$ with valence fermions). Continuum status: continuum-controlled in the
quenched theory. Relevance: it shows the SF machinery can renormalize composite operators, which is
the prerequisite for any SF determination of $Z_O$ for $O=\tfrac14{\rm tr}F^2$ — not done here.

**[E09-07] [V] Continuum static potential in pure YM.**
S. Necco, R. Sommer, “The $N_f=0$ heavy quark potential from short to intermediate distances,”
*Nucl. Phys.* **B622** (2002) 328–346. DOI: 10.1016/S0550-3213(01)00582-X. arXiv:hep-lat/0108008.
**VERIFIED** (arXiv abs + publisher).
Establishes: the static potential for $0.05\,{\rm fm}\le r\le 0.8\,{\rm fm}$ on lattices up to
$64^4$; lattice artefacts in the potential and force analysed; continuum quantities by $a\to0$
extrapolation; scaling violations consistent with Symanzik; $r_0/a$ determined up to $\beta=6.92$.
Continuum status: continuum-extrapolated. This is the cleanest pure-YM continuum observable of the
era below the glueball scale, and a template for how to control $a\to0$.

**[E09-08] [V] Step scaling with dynamical fermions.**
M. Della Morte, R. Frezzotti, J. Heitger, J. Rolf, R. Sommer, U. Wolff (ALPHA Collaboration),
“Computation of the strong coupling in QCD with two dynamical flavours,” *Nucl. Phys.* **B713**
(2005) 378–406. DOI: 10.1016/j.nuclphysb.2005.02.013. arXiv:hep-lat/0411025. **VERIFIED** (arXiv abs).
Establishes: the SF/step-scaling method survives dynamical fermions; a continuum-controlled
$\alpha_s$ in two-flavour QCD. Continuum status: continuum-controlled. Relevance to pure YM: it
demonstrates that the method is not an artefact of the quenched approximation, but it also shows how
much harder the $N_f\neq0$ continuum limit is — relevant to the window programme.

**Cluster A: established vs hoped.** Established: a non-perturbative, continuum-controlled running
coupling and $\Lambda_{\overline{\rm MS}}$ in pure YM; $r_0$; the SF as a renormalization tool.
Hoped: that the same machinery would reach the IR residue. It does not: the SF is a UV probe, and
its finite-volume gap is $\sim 1/L$, not the physical gap.

---

## 2. Cluster B — precision glueball masses (2 entries)

**[E09-09] [V] The anisotropic-lattice glueball spectrum.**
C. J. Morningstar, M. J. Peardon, “The glueball spectrum from an anisotropic lattice study,”
*Phys. Rev. D* **60** (1999) 034509. DOI: 10.1103/PhysRevD.60.034509. arXiv:hep-lat/9901004.
**VERIFIED** (arXiv abs + publisher).
Establishes: the SU(3) pure-gauge glueball spectrum below 4 GeV on several anisotropic lattices with
spatial spacings $0.1$–$0.4$ fm; continuum spin quantum numbers identified; single glueball states
separated from two-glueball and torelon-pair states; the $0^{++}$ mass comes out near
$1.7$ GeV, with the $2^{++}$ and other states above it. Continuum status: continuum-extrapolated in
the spatial lattice spacing, with quoted systematic errors from anisotropy and volume; this is
numerical evidence, not a theorem. The paper measures masses, **not** the residue
$\langle0|\tfrac14{\rm tr}F^2|G\rangle$.

**[E09-10] [V] Glueball spectrum and matrix elements, continuum extrapolated.**
Y. Chen, A. Alexandru, S. J. Dong, T. Draper, I. Horvath, F. X. Lee, K. F. Liu, N. Mathur,
C. Morningstar, M. Peardon, S. Tamhankar, B. L. Young, J. B. Zhang, “Glueball spectrum and matrix
elements on anisotropic lattices,” *Phys. Rev. D* **73** (2006) 014516.
DOI: 10.1103/PhysRevD.73.014516. arXiv:hep-lat/0510074. **VERIFIED** (arXiv abs).
Establishes: an independent anisotropic-lattice determination of the low-lying glueball spectrum
with matrix elements, multiple lattice spacings, and a continuum extrapolation; confirms the
$0^{++}$ near $1.7$ GeV and gives the $2^{++}$, $0^{-+}$, etc. Continuum status:
continuum-extrapolated, with the usual caveats on topology freezing and two-glueball contamination.
Relevance: matrix elements are the closest this era comes to a residue, but the operators are not
the campaign's $\tfrac14{\rm tr}F^2$ spectral weight.

**Cluster B: established vs hoped.** Established: $m_{0^{++}}\approx1.7$ GeV and a low-lying
spectrum, continuum-extrapolated at the few-percent level. Hoped: that the mass would pin the
residue or the IR fixed point. It does not; the residue is a separate observable.

---

## 3. Cluster C — improvement, chirality, and the lattice-fermion toolbox (7 entries)

**[E09-11] [T] Symanzik improvement.**
K. Symanzik, “Continuum limit and improved action in lattice theories. (I). Principles and $\phi^4$
theory,” *Nucl. Phys.* **B226** (1983) 187; “(II). $O(N)$ nonlinear sigma model in perturbation
theory,” *Nucl. Phys.* **B226** (1983) 205. **UNVERIFIED** (pre-arXiv; standard citation, not
re-checked this pass).
Establishes: the effective-field-theory statement that lattice artefacts can be removed
order-by-order by local counterterms; the basis of all $O(a)$ improvement. Continuum status: the
theorem is continuum EFT; the coefficients are tuned non-perturbatively. Relevance: it is what makes
the continuum extrapolations of Clusters A and B controlled rather than ad hoc.

**[E09-12] [O] The clover improvement term.**
B. Sheikholeslami, R. Wohlert, “Improved continuum limit lattice action for QCD with Wilson
fermions,” *Nucl. Phys.* **B259** (1985) 572. **UNVERIFIED** (pre-arXiv; standard citation, not
re-checked this pass).
Establishes: the Sheikholeslami–Wohlert (clover) $O(a)$ improvement of Wilson fermions. Continuum
status: the improvement is a lattice construction; its coefficient is tuned. Relevance: the
standard tool for reducing $O(a)$ artefacts in dynamical-fermion window studies; not directly a
pure-YM gap statement.

**[E09-13] [T] Ginsparg–Wilson relation.**
P. H. Ginsparg, K. G. Wilson, “A remnant of chiral symmetry on the lattice,” *Phys. Rev. D* **25**
(1982) 2649. **UNVERIFIED** (pre-arXiv; standard citation, not re-checked this pass).
Establishes: the algebraic relation $\{D,\gamma_5\}=2aD\gamma_5D$ that defines a lattice Dirac
operator with an exact remnant of chiral symmetry at finite $a$. Continuum status: the relation is
exact on the lattice and has the correct continuum limit; it is the foundation of the overlap/domain
wall programme. Relevance: exact chiral symmetry at finite $a$ is what makes chiral-limit statements
in window studies continuum-controlled rather than extrapolated.

**[E09-14] [NG] Nielsen–Ninomiya no-go.**
H. B. Nielsen, M. Ninomiya, “No-go theorem for regularizing chiral fermions,” *Phys. Lett.* **B105**
(1981) 219; see also *Nucl. Phys.* **B185** (1981) 20 and **B193** (1981) 173. **UNVERIFIED**
(pre-arXiv; standard citation, not re-checked this pass).
Establishes: under locality, correct continuum limit, and exact chiral symmetry, a lattice fermion
spectrum must double. Label `[NG]`: it is a no-go for a named class (naive chiral lattice fermions).
Continuum status: a theorem about lattice regularizations, not about the continuum theory. Relevance:
it explains why exact chiral symmetry requires the Ginsparg–Wilson route, and why domain-wall/overlap
fermions are the clean window tools.

**[E09-15] [O] Domain-wall fermions.**
D. B. Kaplan, “A method for simulating chiral fermions on the lattice,” *Phys. Lett.* **B288** (1992)
342–347. DOI: 10.1016/0370-2693(92)91112-M. arXiv:hep-lat/9206013. **VERIFIED** (arXiv abs).
Establishes: domain-wall fermions — an extra dimension with a domain wall localizes a chiral zero
mode; the residual chiral symmetry breaking is exponentially small in the wall separation.
Continuum status: the construction has a continuum limit with chiral fermions; the residual breaking
is a lattice artefact controlled by $L_s$. Relevance: the fermion formulation used by the LSD
collaboration in E09-29.

**[E09-16] [T] Overlap fermions are exactly massless.**
H. Neuberger, “Exactly massless quarks on the lattice,” *Phys. Lett.* **B417** (1998) 141–144.
DOI: 10.1016/S0370-2693(97)01368-3. arXiv:hep-lat/9707022. **VERIFIED** (arXiv abs).
Establishes: the overlap Dirac operator, which satisfies the Ginsparg–Wilson relation and has exact
zero modes; chiral symmetry is exact at finite $a$. Continuum status: exact lattice chiral symmetry
with a controlled continuum limit. Relevance: the cleanest route to a controlled chiral limit in
many-flavour studies.

**[E09-17] [T] Exact chiral symmetry from the Ginsparg–Wilson relation.**
M. Lüscher, “Exact chiral symmetry on the lattice and the Ginsparg–Wilson relation,” *Phys. Lett.*
**B428** (1998) 342–345. DOI: 10.1016/S0370-2693(98)00423-7. arXiv:hep-lat/9802011. **VERIFIED**
(arXiv abs).
Establishes: the precise sense in which Ginsparg–Wilson fermions have exact chiral symmetry —
a modified chiral transformation under which the action and measure are invariant; the anomaly is
reproduced. Label `[T]`: a theorem inside the GW framework. Continuum status: exact on the lattice,
correct continuum limit. Relevance: it licenses treating the chiral limit as controlled in window
studies, unlike Wilson/staggered fermions.

**Cluster C: established vs hoped.** Established: a complete toolbox (Symanzik improvement, clover,
GW/overlap/domain-wall) that makes continuum extrapolation and the chiral limit controlled.
Hoped: that exact chiral symmetry would settle the window boundary. It does not; it removes one
systematic error, not the finite-size and continuum ones.

---

## 4. Cluster D — Wilson/gradient flow (3 entries)

**[E09-18] [V] Smoothed Wilson loops and the gradient-flow precursor.**
R. Narayanan, H. Neuberger, “Infinite $N$ phase transitions in continuum Wilson loop operators,”
*JHEP* **03** (2006) 064. DOI: 10.1088/1126-6708/2006/03/064. arXiv:hep-th/0601210. **VERIFIED**
(arXiv abs).
Establishes: smoothed Wilson-loop operators on the lattice have a finite non-trivial continuum
limit; their eigenvalue distribution undergoes a phase transition at infinite $N$ as the loop is
dilated. Continuum status: numerical evidence for a continuum limit of the smoothed operators.
Relevance: the conceptual precursor of the Wilson flow; it shows that smoothing can produce
renormalized operators, not that smoothing sees the IR.

**[E09-19] [T] Trivializing maps and the Wilson flow.**
M. Lüscher, “Trivializing maps, the Wilson flow and the HMC algorithm,” *Commun. Math. Phys.* **293**
(2010) 899–919. DOI: 10.1007/s00220-009-0953-7. arXiv:0907.5491. **VERIFIED** (arXiv abs).
Establishes: a rigorous construction of trivializing maps for lattice gauge fields; the Wilson flow
is the gradient flow of the gauge action and provides a trivializing map under suitable conditions;
the Jacobian and the HMC application are derived. Label `[T]`: rigorous inside the stated
hypotheses. Continuum status: the flow is defined on the lattice and has a continuum limit; the
theorem is about the lattice map. Relevance: it gives the flow a mathematical status, but it is a
smoothing/UV statement, not a mass-gap statement.

**[E09-20] [V] The Wilson flow as a renormalized probe.**
M. Lüscher, “Properties and uses of the Wilson flow in lattice QCD,” *JHEP* **08** (2010) 071.
DOI: 10.1007/JHEP08(2010)071. arXiv:1006.4518. **VERIFIED** (arXiv abs). Erratum: *JHEP* **03**
(2014) 092, DOI: 10.1007/JHEP03(2014)092.
Establishes: the flowed gauge field at $t>0$ is smooth; local gauge-invariant expectation values at
flow time $t$ are renormalized and probe scales $\sim\sqrt t$; a one-loop continuum calculation of
$\langle E\rangle$; the $t_0$ scale; the emergence of topological sectors as $a$ is reduced below
$0.1$ fm. Continuum status: $\langle E\rangle$ has a continuum limit after renormalization; $t_0$ is
a lattice-defined, scheme-dependent scale. **Caution:** the 2014 erratum corrected Figs. 2–3 and the
$t_0/a^2$ values in Table 1 after an analysis bug — a concrete reminder that flow quantities are
easy to over-interpret. Relevance: the flow is a UV smoother; it cannot detect the IR mass gap.

**Cluster D: established vs hoped.** Established: a rigorous smoothing map and a renormalized
energy density; the $t_0$ scale. Hoped: that the flow would provide a clean non-perturbative
definition of the coupling or of the gap. It provides a coupling; it does not provide the gap, and
large flow times destroy the IR modes one would need.

---

## 5. Cluster E — the conformal window (13 entries)

**[E09-21] [C] Banks–Zaks: a perturbative IR fixed point.**
T. Banks, A. Zaks, “On the phase structure of vector-like gauge theories with massless fermions,”
*Nucl. Phys.* **B196** (1982) 189–204. DOI: 10.1016/0550-3213(82)90035-9. **VERIFIED** (ADS +
publisher).
Establishes: for $N_f$ just below $11N/2$, the two-loop $\beta$-function has an IR-stable zero at
$g_*^2\propto \epsilon\equiv 11N/2-N_f$; an expansion in $\epsilon$ controls the IR; a large class of
models with fermions in real representations do not break chiral symmetry; the transition between
chiral phases is argued to be first order. Label `[C]`: the IRFP is established inside the
two-loop/$\epsilon$ expansion, not non-perturbatively. Continuum status: continuum perturbation
theory. This is the origin of the “conformal window.”

**[E09-22] [C] Walking technicolor.**
T. Appelquist, J. Terning, L. C. R. Wijewardhana, “Postmodern technicolor,” *Phys. Rev. Lett.* **79**
(1997) 2767–2770. DOI: 10.1103/PhysRevLett.79.2767. arXiv:hep-ph/9706238. **VERIFIED** (arXiv abs).
Establishes: a technicolor framework built on a nontrivial IR fixed point, with a relation between
the electroweak and QCD confinement scales. Label `[C]`: model-building conditional on the IRFP.
Continuum status: continuum model, not lattice. Relevance: it is the phenomenological reason the
window was revived; it does not prove anything about pure YM.

**[E09-23] [C] Phase structure of SU(N) with $N_f$ flavours.**
T. Appelquist, A. Ratnaweera, J. Terning, L. C. R. Wijewardhana, “The phase structure of an SU(N)
gauge theory with $N_f$ flavors,” *Phys. Rev. D* **58** (1998) 105017.
DOI: 10.1103/PhysRevD.58.105017. arXiv:hep-ph/9806472. **VERIFIED** (arXiv abs).
Establishes: a critical $N_f$ separating confining/chirally-broken from conformal behaviour, using
the ladder Schwinger–Dyson approximation; the phase diagram as a function of $N_f$ and $N$.
Label `[C]`: conditional on the truncation. Continuum status: continuum approximation. Relevance:
it frames the lattice target — find the critical $N_f$ non-perturbatively.

**[E09-24] [C] Conformal phase transition.**
V. A. Miransky, K. Yamawaki, “Conformal phase transition in gauge theories,” *Phys. Rev. D* **55**
(1997) 5051; Erratum *Phys. Rev. D* **56** (1997) 3768. **UNVERIFIED** (arXiv ID not confirmed this
pass; standard citation).
Establishes: the scenario in which the order parameter goes to zero with an essential singularity
(Miransky scaling) at the window edge, rather than a mean-field transition; a “walking” regime with
a large anomalous dimension. Label `[C]`: a controlled scenario, not a theorem. Continuum status:
continuum effective description. Relevance: it is the language used to interpret the $N_f=8$
walking signals in E09-30.

**[E09-25] [C] The window for higher representations.**
D. D. Dietrich, F. Sannino, “Walking in the SU(N),” *Phys. Rev. D* **75** (2007) 085018.
DOI: 10.1103/PhysRevD.75.085018. arXiv:hep-ph/0611341. **VERIFIED** (arXiv abs).
Establishes: using an all-orders $\beta$-function, maps the conformal window for SU(N) with fermions
in higher-dimensional representations; fewer flavours are needed to reach the walking regime;
discusses the spectrum near the fixed point. Label `[C]`: the all-orders $\beta$-function is a
conjecture, not derived from first principles. Continuum status: continuum gap-equation/all-orders
approximation. Relevance: it broadens the window beyond the fundamental representation, which
matters for model-building but not for pure YM.

**[E09-26] [C] A supersymmetry-inspired all-orders $\beta$-function.**
T. A. Ryttov, F. Sannino, “Supersymmetry inspired QCD beta function,” *Phys. Rev. D* **78** (2008)
065001. DOI: 10.1103/PhysRevD.78.065001. arXiv:0711.3745. **VERIFIED** (arXiv abs).
Establishes: an all-orders $\beta$-function inspired by the NSVZ formula; it bounds the conformal
window; for one adjoint Weyl fermion it matches the supersymmetric YM $\beta$-function; it computes
the running of the pure-YM coupling and compares the deviation from two loops with lattice data.
Label `[C]`: the $\beta$-function is conjectural outside SUSY. Continuum status: continuum
all-orders ansatz. Relevance: it is the sharpest perturbative-inspired statement about the window,
and it explicitly discusses the pure-YM ($N_f=0$) running — but still not a non-perturbative proof.

**[E09-27] [V] Lattice evidence for the window: $N_f=12$ conformal, $N_f\le8$ confining.**
T. Appelquist, G. T. Fleming, E. T. Neil, “Lattice study of the conformal window in QCD-like
theories,” *Phys. Rev. Lett.* **100** (2008) 171607. DOI: 10.1103/PhysRevLett.100.171607.
arXiv:0712.0609. **VERIFIED** (arXiv abs + INSPIRE). Erratum: *Phys. Rev. Lett.* **102** (2009)
149902.
Establishes: SU(3) gauge theory with $N_f$ fundamental Dirac fermions; evidence that the IR is
conformal for $12\le N_f\le16$, while confinement and chiral symmetry breaking persist for
$N_f\le8$. Continuum status: lattice evidence at finite $a,L,m$; the continuum limit is not taken.
This is the paper that set the modern $N_f=12$ target — and its erratum is a reminder that the
classification is delicate.

**[E09-28] [V] MCRG: $N_f=8$ QCD-like, $N_f=12$ IRFP or slow walking.**
A. Hasenfratz, “Conformal or walking? Monte Carlo renormalization group studies of SU(3) gauge
models with fundamental fermions,” *Phys. Rev. D* **82** (2010) 014506.
DOI: 10.1103/PhysRevD.82.014506. arXiv:1004.1004. **VERIFIED** (arXiv abs).
Establishes: a two-lattice matching MCRG study of $N_f=8$ and $N_f=12$; the $N_f=8$ perturbative and
confining regimes can be connected, verifying its QCD-like nature; for $N_f=12$ the data favour an
IRFP and a conformal phase, but are also consistent with very slow walking; $\gamma_m$ is measured
and found close to the free-field value. Continuum status: lattice evidence at finite $a$; the MCRG
fixed-point search is a finite-lattice construction. Relevance: it is the strongest
“finite-size-scaling” counterpoint to the simple $N_f=12$ conformal claim.

**[E09-29] [V] LSD: chiral condensate enhancement toward the window.**
T. Appelquist, A. Avakian, R. Babich, R. C. Brower, M. Cheng, M. A. Clark, S. D. Cohen,
G. T. Fleming, J. Kiskis, E. T. Neil, J. C. Osborn, C. Rebbi, D. Schaich, P. Vranas (LSD
Collaboration), “Toward TeV conformality,” *Phys. Rev. Lett.* **104** (2010) 071601.
DOI: 10.1103/PhysRevLett.104.071601. arXiv:0910.2224. **VERIFIED** (arXiv abs + INSPIRE).
Establishes: domain-wall fermions, SU(3) with $N_f=2$–$6$ massless Dirac fermions; the chiral
condensate $\langle\bar\psi\psi\rangle/F^3$ agrees with QCD at $N_f=2$ and is significantly enhanced
at $N_f=6$, suggesting a larger enhancement toward the critical $N_f$ for the
confinement-to-conformality transition. Continuum status: lattice evidence with domain-wall
fermions at finite $a,L$; chiral symmetry is exact at finite $a$, but the continuum limit is not
taken. Relevance: it brackets the window from below ($N_f=6$ still confining) and sets up the
$N_f=8$–$12$ controversy.

**[E09-30] [V] LatKMI: walking signals at $N_f=8$.**
Y. Aoki, T. Aoyama, M. Kurachi, T. Maskawa, K.-I. Nagai, H. Ohki, A. Shibata, K. Yamawaki,
T. Yamazaki (LatKMI Collaboration), “Walking signals in $N_f=8$ QCD on the lattice,” *Phys. Rev. D*
**87** (2013) 094511. DOI: 10.1103/PhysRevD.87.094511. arXiv:1302.6859. **VERIFIED** (arXiv abs).
Establishes: HISQ simulations of eight-flavour QCD; data consistent with spontaneous chiral symmetry
breaking near the chiral limit, but also finite-size hyperscaling at larger fermion masses suggesting
a large anomalous dimension $\gamma_m\sim1$ — approximate IR conformality, as in walking technicolor.
Continuum status: lattice evidence at finite $a,L,m$; staggered taste-breaking is a systematic.
Relevance: $N_f=8$ is the closest lattice proxy to the window edge; the “walking” interpretation is
not a continuum theorem.

**[E09-31] [V] LatKMI: a light composite scalar at $N_f=12$.**
Y. Aoki, T. Aoyama, M. Kurachi, T. Maskawa, K.-I. Nagai, H. Ohki, E. Rinaldi, A. Shibata,
K. Yamawaki, T. Yamazaki (LatKMI Collaboration), “Light composite scalar in twelve-flavor QCD on the
lattice,” *Phys. Rev. Lett.* **111** (2013) 162001. DOI: 10.1103/PhysRevLett.111.162001.
arXiv:1305.6006. **VERIFIED** (arXiv abs).
Establishes: improved-staggered simulations of twelve-flavour QCD at several fermion masses; a
flavour-singlet scalar state lighter than the pion; interpreted as evidence for IR conformality and a
possible light composite Higgs/technidilaton. Continuum status: lattice evidence at finite $a,L,m$;
the light scalar could be a finite-volume or taste artefact, and the $N_f=12$ classification is
contested (see E09-32). Relevance: it is the phenomenological high point of the window programme,
not a rigorous statement about pure YM.

**[E09-32] [V] The counter-claim: $N_f=12$ below the window.**
Z. Fodor, K. Holland, J. Kuti, D. Nógrádi, C. Schröder, “Twelve massless flavors and three colors
below the conformal window,” *Phys. Lett.* **B703** (2011) 348. **UNVERIFIED** (journal citation
standard; not re-checked this pass; arXiv ID not confirmed).
Establishes: the opposite classification — twelve massless flavours in SU(3) are below the conformal
window (confining/chirally broken), based on finite-temperature and chiral-condensate observables.
Continuum status: lattice evidence at finite $a,L$; the tension with E09-27/E09-31 is the central
unresolved systematic of the window programme. Relevance: it shows that the “$N_f=12$ is conformal”
claim is not settled; any v3 use of the window must cite both sides.

**[E09-33] [C] Review of conformal dynamics.**
F. Sannino, “Conformal dynamics for TeV physics and cosmology,” *Acta Phys. Polon.* **B40** (2009)
3533. arXiv:0911.0931. **UNVERIFIED** (arXiv preprint verified; journal volume/page not confirmed
this pass).
Establishes: a 207-page review of the conformal window, walking dynamics, $\gamma_m$, and the
lattice/model-building status as of 2009. Label `[C]`: a review of conditional statements. Continuum
status: mixed; it clearly separates perturbative window estimates from lattice evidence. Relevance:
the best single entry point to the era's window literature.

**Cluster E: established vs hoped.** Established: a perturbative IRFP near $N_f=11N/2$
(Banks–Zaks); a set of conditional all-orders $\beta$-functions; lattice evidence that $N_f=12$ is
near-conformal (contested) and $N_f=8$ is walking; exact chiral symmetry tools to control the
chiral limit. Hoped: a continuum-controlled critical $N_f$ and a proof that $N_f=0$ is on the
confining side. Neither exists.

---

## 6. Continuum-controlled vs lattice artefact

| Claim from this era | Status | Why |
|---|---|---|
| $\alpha_{\rm SF}(\mu)$ in pure YM at high $\mu$ | Continuum-controlled | SF scheme defined in the continuum; $a\to0$ extrapolation at each step; $O(a^2)$ artefacts removed; two-loop matching to $\overline{\rm MS}$ (E09-02, E09-04). |
| $\Lambda_{\overline{\rm MS}}$ in pure YM (few %) | Continuum-controlled (few %) | Step scaling + perturbative matching; dominant errors are statistics and truncation, not $a$. |
| $r_0/a$ up to $\beta=6.92$; static potential $0.05$–$0.8$ fm | Continuum-extrapolated | Necco–Sommer; Symanzik scaling violations fitted (E09-07). |
| $m_{0^{++}}\approx1.7$ GeV; spectrum below 4 GeV | Continuum-extrapolated, few-% | Anisotropic lattices, several $a_s$, volume and two-glueball systematics studied (E09-09, E09-10). Not a theorem. |
| Glueball matrix elements | Continuum-extrapolated with caveats | Same ensembles; topology freezing and operator mixing are the limiting systematics (E09-10). |
| $t_0$ scale | Lattice-defined; continuum limit exists | $t_0$ is a flow-time scale, scheme-dependent, not a physical mass; erratum corrected values (E09-20). |
| $\langle E\rangle$ at flow time $t$ | Continuum-controlled after renormalization | One-loop coefficient computed in continuum DR; flow smooths UV; finite continuum limit (E09-19, E09-20). |
| Existence of an IRFP at $N_f=12$ | Lattice evidence, **not** continuum-controlled | Finite $a,L,m$; AFN/LSD vs Fodor conflict; hyperscaling vs $\chi$SB; no $a\to0$ extrapolation of the fixed-point existence (E09-27, E09-28, E09-31, E09-32). |
| Critical $N_f$ (window boundary) | Not a rigorous continuum quantity | Depends on the definition of “window”; lattice gives brackets, not a continuum critical point. |
| Chiral symmetry at finite $a$ (GW/overlap/domain wall) | Exact on the lattice, correct continuum limit | Ginsparg–Wilson relation; no tuning; residual breaking controlled by $L_s$ (E09-13, E09-15, E09-16, E09-17). |
| Chiral symmetry with Wilson/staggered fermions | Lattice artefact | Explicit breaking; chiral limit is an extrapolation; staggered taste-breaking. |
| $O(a)$ improvement | Continuum EFT statement | Symanzik; local counterterms; coefficients tuned non-perturbatively (E09-11, E09-12). |
| Strong-coupling mass gap at small $\beta$ | Rigorous on the lattice | Osterwalder–Seiler / SZZ; **not** continued to the continuum. |
| Mass gap at physical $\beta$ | **Open** | Clay. Nothing in this era closes it. |

**Reading the table.** The left column of the era — step scaling, $r_0$, glueball masses, flow —
is genuinely continuum-controlled (with quoted errors). The right column — the window boundary,
the $N_f=12$ classification, the chiral limit with non-GW fermions — is lattice evidence at finite
$(a,L,m)$ and must not be promoted to a continuum statement. The single most important distinction
for v3: **a continuum-controlled running coupling is not a continuum-controlled mass gap.**

---

## 7. Why pure YM ($N_f=0$) being outside the window is only *perturbatively* rigorous

### 7.1 The perturbative statement, exactly

In a mass-independent scheme,
`β(g) = -β₀ g³/(16π²) - β₁ g⁵/(16π²)² - …`,
with
`β₀ = 11N/3 - 2N_f/3`,
`β₁ = 34N²/3 - (13N/3 - 1/N) N_f`.

For $N_f=0$ and any $N$: `β₀ = 11N/3 > 0`, `β₁ = 34N²/3 > 0`. The two-loop $\beta$-function is
negative for all $g>0$ and has **no positive zero**: the coupling grows into the IR, and the
two-loop truncation predicts no IR fixed point.

For $N=3$, the two-loop window is `8.05 < N_f < 16.5`: at `N_f=16.5`, `β₀=0`; at `N_f≈8.05`,
`β₁=0`; between them `β₀>0`, `β₁<0`, so `g_*² = -16π² β₀/β₁ > 0`. This is the Banks–Zaks
IRFP, controlled by `ε = 11N/2 - N_f`. The lower edge `N_f≈8.05` is where the two-loop coupling
becomes `O(1)`; perturbation theory is already failing there. For $N_f=0$ one is infinitely far
from the upper edge, and the two-loop “no zero” is not a statement about the exact $\beta$ at
strong coupling.

### 7.2 What “perturbatively rigorous” does and does not mean

- **What is rigorous:** `β₀` and `β₁` are universal (scheme-independent) coefficients of the
  mass-independent $\beta$-function, computable in BPHZ; the statement “the two-loop polynomial has
  no positive zero for $N_f=0$” is a theorem about that polynomial.
- **What is not rigorous:** the exact $\beta$-function is not the polynomial. The physical IR is at
  `g~1`, where the truncation is uncontrolled. The exact $\beta$-function could, in principle, have
  a strong-coupling zero. Lattice evidence says it does not (confinement), but that is evidence,
  not proof.
- **A prior problem:** even the object “the $\beta$-function of the continuum pure-YM theory” is not
  rigorously defined non-perturbatively until the continuum measure is constructed (Clay G1). A
  statement about the IR phase presupposes the theory exists.

### 7.3 What would make it non-perturbatively rigorous

Five sufficient routes, roughly from strongest to weakest:

1. **Construct the 4D continuum SU(N) YM measure with Osterwalder–Schrader axioms (G1), then prove
   exponential clustering / a positive mass gap (G2).** A gapped theory is not a gapless CFT, so
   this implies “outside the window.” This is the Clay problem itself; it is sufficient but not
   necessary for the window statement.
2. **Prove a non-perturbative no-IRFP theorem:** a rigorous RG map or exact $\beta$-function for
   $N_f=0$ has no zero at any coupling. This is weaker than a gap proof and is the natural target
   for a “new monotone” (the campaign's sixth-species search).
3. **Prove monotonicity in $N_f$:** e.g., a universal IR observable (mass gap, $\gamma_m$, or a
   spectral weight) is monotone in $N_f$, and at some $N_f^{+}$ the Banks–Zaks IRFP is rigorously
   established (large $N_f$, perturbative). Then $N_f=0$ lies on the gapped side. **No such
   monotonicity theorem is known.** This is a concrete, unoccupied attack surface.
4. **Exclude the gapless CFT directly:** a rigorous $a$-theorem plus anomaly matching, or a
   bootstrap bound, that no interacting 4D CFT with the symmetries and anomaly content of pure YM
   exists. Era 10 is the relevant literature; the known constraints do not exclude the candidate.
5. **Give a rigorous definition of the conformal window** (the set of $N_f$ for which the continuum
   theory has an IRFP of a specified universality class), then prove $N_f=0$ is not in it. Today the
   window is a physicists' notion; the critical $N_f$ is not a rigorously defined continuum
   quantity. Route 5 is logically prior to a sharp statement of the problem.

### 7.4 The logical gap between “outside the window” and “mass gap”

“Outside the conformal window” means “no IRFP.” That is not, by itself, the same as “gapped” unless
one also classifies the possible gapless phases. A gapless interacting 4D theory is *expected* to be
a CFT, but this is not proven. The campaign's precise alternative is exactly (iii): the YM measure is
a gapless interacting 4D CFT. So the conformal-window literature is the map of the **alternative**,
not a shortcut to the gap. Even a proof of “no IRFP at $N_f=0$” would still leave the possibility of
a non-conformal gapless phase; closing that requires the same non-perturbative control.

### 7.5 What the era actually established

- **Lattice:** pure YM confines; the $0^{++}$ glueball is at $\approx1.7$ GeV; the static potential
  is continuum-extrapolated; the running coupling is continuum-controlled in the UV. All at finite
  $(a,L)$, with quantified systematics — **evidence, not theorem.**
- **Strong-coupling lattice:** a rigorous mass gap at small $\beta$ (Osterwalder–Seiler / SZZ), but
  only in the strong-coupling region, not the continuum.
- **Perturbation theory:** no positive two-loop IR zero at $N_f=0$; a perturbative IRFP near
  $N_f=11N/2$. **This is the only rigorous “outside the window” statement.**

**Sharp summary.**
`Theorem (perturbative):` for $N_f=0$, `β₀,β₁>0`; the two-loop $\beta$ has no positive zero.
`Conditional:` if the exact $\beta$ is approximated by the two-loop truncation down to the IR, there
is no IRFP.
`Open:` the exact $\beta$ has no IRFP; equivalently, pure YM is not a gapless CFT. This is Clay.

---

## 8. Position in 14 eras

| # | Window | Theme | What this era's angle says |
|---|---|---|---|
| 01 | 1954–1965 | origin of the gauge-mass problem | The question is posed kinematically: can a gauge theory have a mass without breaking the symmetry? No tool yet. |
| 02 | 1963–1973 | quantization, renormalization, asymptotic freedom | AF makes the UV safe; the two-loop $\beta$ computed here is the seed of Banks–Zaks. The IR is still uncontrolled. |
| 03 | 1974–1979 | instantons, Savvidy vacuum, confinement pictures | First evidence that the IR is strong and non-perturbative; no controlled calculational scheme. |
| 04 | 1974–1982 | lattice gauge theory and rigorous strong coupling | Wilson's formulation plus a rigorous strong-coupling mass gap; the continuum limit is open. |
| 05 | 1982–1990 | constructive RG and the gluon condensate | OPE and condensates; a language for non-perturbative matrix elements, but no gap proof. |
| 06 | 1985–1995 | confinement mechanisms and lattice phenomenology | Flux tubes and glueballs become numerical; still no continuum control. |
| 07 | 1993–2000 | Dyson–Schwinger and exact RG | Continuum truncations and scheme dependence; no controlled continuum limit. |
| 08 | 1998–2005 | holography, large-N, AdS/QCD | Qualitative gap models; not pure YM, not rigorous. |
| **09** | **2002–2010** | **lattice precision, step scaling, conformal window** | **This report: continuum-controlled running coupling and glueball masses; Wilson flow; the window revived as the gapless alternative. Quantitative, not rigorous; $N_f=0$ outside the window only perturbatively.** |
| 10 | 2008–2014 | CFT constraints: collider bounds, $a$-theorem, bootstrap | Constrains the gapless alternative (iii); no direct YM gap. |
| 11 | 2013–2018 | rigorous probabilistic and constructive QFT | 2D/3D constructions; 4D YM still open. |
| 12 | 2018–2021 | resurgence, transseries, modern non-perturbative methods | Asymptotic series and renormalons; no gap proof. |
| 13 | 2021–2024 | rigorous frontier and large-scale lattice | Bootstrap and exascale lattice; still no 4D YM gap. |
| 14 | 2024–2026 | formalization, AI-assisted mathematics, current frontier | Lean/formalization and the present campaign; the gap remains open. |

---

## 9. Gaps and attack surfaces for v3

- **Step scaling is UV-only.** E09-02/E09-07 give a continuum-controlled $\alpha_{\rm SF}$ and
  $r_0$, but the SF is defined at $\mu=1/L$ with a finite-volume gap $\sim1/L$. It cannot bound the
  infinite-volume $0^{++}$ residue. *Attack:* extend step scaling to composite-operator matching for
  $O=\tfrac14{\rm tr}F^2$ — compute $Z_O(\mu)$ non-perturbatively and transport $f_1^2$ from the
  strong-coupling end. The campaign's multiplicative transport already assumes this; the missing
  input is a continuum-controlled $Z_O$ across the Savvidy crossing.
- **The glueball continuum limit is few-percent, not exact.** E09-09/E09-10 give
  $m_{0^{++}}\approx1.7$ GeV with $a_s\to0$ extrapolations, but the finest $a_s\sim0.1$ fm and the
  quoted systematics are 3–5%. *Attack:* a dedicated fine-lattice run at $\beta\in[2.4,2.6]$ with
  the $t_0$ scale, controlled topology, and a simultaneous extraction of $m_{0^{++}}$ and the
  residue $\langle0|\tfrac14{\rm tr}F^2|G\rangle$.
- **The residue $f_1^2$ is essentially unmeasured.** None of the era's glueball papers extracts the
  $\tfrac14{\rm tr}F^2$ spectral weight; E09-10's matrix elements are different operators.
  *Attack:* measure the overlap on the same ensembles used for the mass; a nonzero continuum limit
  is direct evidence for the residue, a vanishing one is a no-go signal.
- **The Wilson flow cannot see the IR.** E09-19/E09-20 smooth modes with $p^2\gtrsim1/t$; $t_0$ is a
  scheme-dependent scale, and the 2014 erratum shows how easy it is to over-read flow quantities.
  *Attack:* use the flow only for UV renormalization/$Z_O$ and for topology, never as a gap detector;
  cross-check $t_0$ against $r_0$ and $\sqrt\sigma$.
- **The window boundary is not a rigorous object, and $N_f=12$ is contested.** E09-27/E09-31
  (conformal) vs E09-32 (confining) differ at the few-percent level; E09-28 finds IRFP-or-slow-walking.
  *Attack:* define the window via universal IR observables (existence of an IRFP, $\gamma_m$,
  spectrum) and a finite-size-scaling protocol, then scan $N_f$ with the same action and scale
  setting. Do not cite one side of the $N_f=12$ controversy as settled.
- **No $N_f$-monotonicity theorem.** A proof that a universal IR observable is monotone in $N_f$,
  plus a rigorous Banks–Zaks IRFP at large $N_f$, would put $N_f=0$ on the gapped side. This is the
  cleanest “new tool” target: a rigorous inequality in $N_f$, not another $\beta$-function fit.
- **SF boundary conditions break translation invariance.** The SF's finite-volume gap is not the
  physical gap, and the $L\to\infty$ limit of the SF is not the periodic torus. *Attack:* use the SF
  to define the renormalized coupling and $Z_O$, then match to a periodic-volume calculation before
  taking $L\to\infty$; do not conflate the two limits.
- **Perturbative rigor stops at $g\sim1$.** The two-loop no-zero statement for $N_f=0$ is a theorem
  about a truncated series; the exact strong-coupling $\beta$ is open. *Attack:* a rigorous
  non-perturbative $\beta$-function or a no-IRFP theorem at strong coupling — the sixth-species
  target. This is the precise missing step between “perturbatively outside the window” and “gapped.”
- **Topology and the $\theta$-angle are uncontrolled.** Flow/topology measurements freeze at
  $a\gtrsim0.1$ fm; the window literature ignores $\theta$. *Attack:* open boundaries or GW
  fermions to control topology; verify that the $0^{++}$ residue is $\theta$-independent at
  $\theta=0$ before using it as a gap order parameter.

---

## 10. Verification ledger

`VERIFIED` = bibliographic data checked this pass against the arXiv abstract page or a
publisher/ADS/CERN record. `UNVERIFIED` = not checked this pass; the citation is the standard one,
but no field was invented. The shared WebSearch budget was exhausted mid-task; direct arXiv fetch
was used for the core.

| # | Short cite | Checked against | Status | Note |
|---|---|---|---|---|
| E09-01 | Lüscher–Narayanan–Weisz–Wolff, SF, NPB 384 (1992) 168 | arXiv:hep-lat/9207009 abs | VERIFIED | DOI confirmed |
| E09-02 | Lüscher–Sommer–Weisz–Wolff, running coupling, NPB 413 (1994) 481 | arXiv:hep-lat/9309005 abs + publisher | VERIFIED | DOI confirmed |
| E09-03 | Sommer, $r_0$, NPB 411 (1994) 839 | arXiv:hep-lat/9310022 abs | VERIFIED | DOI confirmed |
| E09-04 | Lüscher–Weisz, two-loop matching, NPB 452 (1995) 234 | arXiv:hep-lat/9505011 abs | VERIFIED | DOI confirmed |
| E09-05 | Sint–Sommer, SF one-loop, NPB 465 (1996) 71 | arXiv:hep-lat/9508012 abs | VERIFIED | DOI confirmed |
| E09-06 | Capitani–Lüscher–Sommer–Wittig, quark mass, NPB 544 (1999) 669 | arXiv:hep-lat/9810063 abs | VERIFIED | DOI confirmed |
| E09-07 | Necco–Sommer, static potential, NPB 622 (2002) 328 | arXiv:hep-lat/0108008 abs + publisher | VERIFIED | DOI confirmed |
| E09-08 | Della Morte et al., dynamical $\alpha_s$, NPB 713 (2005) 378 | arXiv:hep-lat/0411025 abs | VERIFIED | DOI confirmed |
| E09-09 | Morningstar–Peardon, glueballs, PRD 60 (1999) 034509 | arXiv:hep-lat/9901004 abs + publisher | VERIFIED | DOI confirmed |
| E09-10 | Chen et al., glueball matrix elements, PRD 73 (2006) 014516 | arXiv:hep-lat/0510074 abs | VERIFIED | DOI confirmed |
| E09-11 | Symanzik, improvement, NPB 226 (1983) 187/205 | not checked this pass | UNVERIFIED | Pre-arXiv; standard citation |
| E09-12 | Sheikholeslami–Wohlert, clover, NPB 259 (1985) 572 | not checked this pass | UNVERIFIED | Pre-arXiv; standard citation |
| E09-13 | Ginsparg–Wilson, PRD 25 (1982) 2649 | not checked this pass | UNVERIFIED | Pre-arXiv; standard citation |
| E09-14 | Nielsen–Ninomiya, PLB 105 (1981) 219 | not checked this pass | UNVERIFIED | Pre-arXiv; standard citation |
| E09-15 | Kaplan, domain wall, PLB 288 (1992) 342 | arXiv:hep-lat/9206013 abs | VERIFIED | DOI confirmed |
| E09-16 | Neuberger, overlap, PLB 417 (1998) 141 | arXiv:hep-lat/9707022 abs | VERIFIED | DOI confirmed |
| E09-17 | Lüscher, GW chiral symmetry, PLB 428 (1998) 342 | arXiv:hep-lat/9802011 abs | VERIFIED | DOI confirmed |
| E09-18 | Narayanan–Neuberger, smoothed loops, JHEP 0603:064 | arXiv:hep-th/0601210 abs | VERIFIED | DOI confirmed |
| E09-19 | Lüscher, trivializing maps, CMP 293 (2010) 899 | arXiv:0907.5491 abs | VERIFIED | DOI confirmed |
| E09-20 | Lüscher, Wilson flow, JHEP 1008:071 | arXiv:1006.4518 abs | VERIFIED | Erratum JHEP 03 (2014) 092 confirmed |
| E09-21 | Banks–Zaks, NPB 196 (1982) 189 | ADS + publisher | VERIFIED | DOI confirmed |
| E09-22 | Appelquist–Terning–Wijewardhana, PRL 79 (1997) 2767 | arXiv:hep-ph/9706238 abs | VERIFIED | DOI confirmed |
| E09-23 | Appelquist–Ratnaweera–Terning–Wijewardhana, PRD 58 (1998) 105017 | arXiv:hep-ph/9806472 abs | VERIFIED | DOI confirmed |
| E09-24 | Miransky–Yamawaki, PRD 55 (1997) 5051 | not checked this pass | UNVERIFIED | Standard citation; arXiv ID not confirmed |
| E09-25 | Dietrich–Sannino, walking SU(N), PRD 75 (2007) 085018 | arXiv:hep-ph/0611341 abs | VERIFIED | DOI confirmed |
| E09-26 | Ryttov–Sannino, all-orders $\beta$, PRD 78 (2008) 065001 | arXiv:0711.3745 abs | VERIFIED | DOI confirmed |
| E09-27 | Appelquist–Fleming–Neil, PRL 100 (2008) 171607 | arXiv:0712.0609 abs + INSPIRE | VERIFIED | Erratum PRL 102 (2009) 149902 confirmed |
| E09-28 | Hasenfratz, MCRG, PRD 82 (2010) 014506 | arXiv:1004.1004 abs | VERIFIED | DOI confirmed |
| E09-29 | LSD, toward TeV conformality, PRL 104 (2010) 071601 | arXiv:0910.2224 abs + INSPIRE | VERIFIED | Correct arXiv ID; 0910.2226 is a solar paper |
| E09-30 | Aoki et al. (LatKMI), $N_f=8$ walking, PRD 87 (2013) 094511 | arXiv:1302.6859 abs | VERIFIED | DOI confirmed |
| E09-31 | Aoki et al. (LatKMI), light scalar $N_f=12$, PRL 111 (2013) 162001 | arXiv:1305.6006 abs | VERIFIED | DOI confirmed |
| E09-32 | Fodor et al., 12 flavours below window, PLB 703 (2011) 348 | not checked this pass | UNVERIFIED | Standard citation; arXiv ID not confirmed |
| E09-33 | Sannino, conformal dynamics review, Acta Phys. Polon. B40 (2009) 3533 | arXiv:0911.0931 abs (preprint only) | UNVERIFIED | arXiv verified; journal volume/page not confirmed |

**Counts:** 33 entries; 26 VERIFIED; 7 UNVERIFIED. The 7 unverified are either pre-arXiv
foundational papers (E09-11–E09-14) or citations whose journal record could not be re-checked after
the search budget was exhausted (E09-24, E09-32, E09-33). None of the unverified entries is
load-bearing for the era's central claims; E09-32 is load-bearing only as the **counter-claim** to
the $N_f=12$ conformal interpretation, and is labelled as such.
