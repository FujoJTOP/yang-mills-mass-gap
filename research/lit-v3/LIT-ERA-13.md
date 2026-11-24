# LIT-ERA-13 — 2021–2024: the rigorous frontier and large-scale lattice

Verification date: 2026-11-24. WebSearch quota was exhausted partway through this survey; all
remaining checks were done against the arXiv API / arXiv abstract pages (publisher-side metadata)
and, where noted, zbMATH/MR snippets already retrieved. Every entry below is either `VERIFIED`
(identifier and metadata checked) or `UNVERIFIED` (flagged, not guessed). No bibliographic field
in this file is invented.

Labels: **[T]** theorem, **[C]** conditional theorem, **[NG]** no-go for a named class,
**[O]** exact open/missing step or proposal, **[V]** numerical evidence.

---

## 0. Headline assessment

1. The era's one unconditional 4D mass-gap theorem is **Shen–Zhu–Zhu (CMP 400, 2023)**. It is a
   *strong-coupling* theorem: for SU(N) in d=4 it requires |β| < 1/48 in the paper's 't Hooft-scaled
   normalization. It proves a strictly positive, volume-uniform mass gap, but the constant is not
   exhibited in closed form and the result never approaches the physical coupling.
2. The only rigorous *weak-coupling* correlation decay of the era is **Adhikari–Cao** — and it is
   for **finite gauge groups**, not continuous SU(N). There is therefore a theorem at strong
   coupling (continuous groups) and a theorem at weak coupling (finite groups), and **nothing for
   continuous SU(N) in between or at weak coupling in 4D**. That split is the structural fact of
   the era.
3. The "CPS/CNS follow-up" the campaign asked about is **not in this era**. Cao–Nissim–Sheffield
   (area law, 't Hooft regime) is arXiv:2505.16585 (May 2025) and arXiv:2509.04688 (Sep 2025);
   Nissim's U(N) 't Hooft-regime mass gap is arXiv:2510.22788 (Oct 2025). Within 2021–2024 the
   rigorous follow-ups are Adhikari–Cao (weak coupling, finite groups), Cao–Park–Sheffield and
   Borga–Cao–Shogren-Knaak (exact surface expansions, no gap), and the 3D constructions
   Cao–Chatterjee / Chandra–Chevyrev–Hairer–Shen. The CPS 2025 papers belong to Era 14.
4. The lattice side consolidated rather than surprised: **Athenodorou–Teper (JHEP 12 (2021) 082)**
   is the era's reference continuum extrapolation of the pure-YM glueball spectrum and string
   tension for N=2…12, with the 't Hooft scaling g²(a) ∝ 1/N verified and Λ_MSbar/√σ extracted.
5. The 12-flavor conformal window is **still contested numerically, not settled**: Fodor et al.
   (LSD, PRD 94 (2016) 091501) eliminated the earlier IRFP evidence in the 6.0<g²<6.4 window in
   their scheme, and Peterson–Hasenfratz (PRD 109 (2024) 114507) report an IRFP at
   g²_GF⋆ = 6.60(62) in a different (gradient-flow, infinite-volume) scheme. Neither is a theorem.
   Quantum-simulation, tensor-network and ML work is *proposal or numerical method*, not evidence
   for the 4D mass gap.

---

## 1. Rigorous strong coupling and the mass gap (11 entries)

**13.01 [T] VERIFIED** — H. Shen, R. Zhu, X. Zhu, "A stochastic analysis approach to lattice
Yang–Mills at strong coupling," *Comm. Math. Phys.* **400** (2023), no. 2, 805–851.
DOI 10.1007/s00220-022-04609-1; arXiv:2204.12737 [math.PR].
Proves, for any d>1 and Wilson-type lattice YM with 't Hooft scaling βN: under
|β| < (N−2)/(32(d−1)N) for SO(N), or |β| < 1/(16(d−1)) for SU(N), (i) uniqueness of the invariant
measure of the Langevin dynamic on the whole lattice; (ii) exponential ergodicity in Wasserstein
distance; (iii) convergence of finite-volume measures to the infinite-volume one, with
Log-Sobolev and Poincaré inequalities; (iv) exponential decay of correlations for a large class of
observables, i.e. a strictly positive mass gap; (v) factorized Wilson-loop correlations and
deterministic large-N limits. In d=4 the SU(N) hypothesis is **β < 1/48**. The bound is uniform in
the lattice volume. It does **not** give a numerical lower bound on the gap, does not treat the
continuum limit, and does not reach weak coupling.

**13.02 [T] VERIFIED** — A. Adhikari, S. Cao, "Correlation decay for finite lattice gauge theories
at weak coupling," *Ann. Probab.* **53** (2025), no. 1, 140–174; arXiv:2202.10375 [math.PR]
(preprint Feb 2022; journal 2025).
At **weak coupling** (large β) for lattice gauge theories with **finite (possibly non-Abelian)
gauge groups**, proves exponential decay of correlations for a wide class of gauge-invariant
functions, including arbitrary functions of Wilson-loop observables. This is the era's only
weak-coupling correlation-decay theorem. It does **not** cover continuous SU(N); finite groups have
no asymptotic freedom and no continuum limit, so it cannot be transported to the physical theory.
The dimension is not stated on the abstract page and was not independently checked — flagged.

**13.03 [C] VERIFIED** — S. Cao, S. Chatterjee, "A state space for 3D Euclidean Yang–Mills
theories," *Comm. Math. Phys.* **405** (2024), no. 1, Paper No. 3, 69 pp.; MR4689958;
arXiv:2111.12813 [math.PR].
Conditional on GFF-like short-distance behaviour, constructs 3D Euclidean YM on the unit torus as
random distributional gauge orbits and defines regularized Wilson loops. The key tool (YM heat
flow from GFF-like initial data) is in the companion paper, *Comm. Partial Differential Equations*
**48** (2023) 209–251, arXiv:2111.10652. There is **no mass gap** here and no 4D statement; the
condition is the whole content. Useful to v3 as a precise formulation of "what a 4D construction
would have to assume."

**13.04 [T] VERIFIED** — A. Chandra, I. Chevyrev, M. Hairer, H. Shen, "Stochastic quantisation of
Yang–Mills–Higgs in 3D," *Invent. Math.* **237** (2024), 541–696.
DOI 10.1007/s00222-024-01264-2; arXiv:2201.03487 [math.PR].
Well-posedness and regularity of the stochastic-quantisation dynamic for 3D YM–Higgs. A local/
global construction theorem, not a mass-gap result, and strictly 3D. Its relevance to v3 is
methodological: it is the current high-water mark for *constructing* a non-Abelian gauge dynamic
in a physically relevant dimension, and 4D remains out of reach.

**13.05 [T] VERIFIED** — S. Cao, M. Park, S. Sheffield, "Random surfaces and lattice Yang–Mills,"
arXiv:2307.06790 [math.PR] (Jul 2023); to appear in *Comm. Amer. Math. Soc.*
Exact rewriting of Wilson-loop expectations for compact G as weighted sums over embedded maps
(surfaces), for **all N ≥ 1, all β > 0, all d ≥ 2**, with Weingarten edge weights and weight
∝ N^χ. Covers U(N), SU(N), O(N), SO(N), Sp(N/2). This is an identity, not an estimate: **no mass
gap and no area law** are claimed. It is the combinatorial substrate that the 2025 CPS area-law
work builds on.

**13.06 [T] VERIFIED** — J. Borga, S. Cao, J. Shogren-Knaak, "Surface sums for lattice Yang–Mills
in the large-N limit," arXiv:2411.11676 [math.PR] (Nov 2024).
In the **large-N, strongly coupled** limit, each Wilson-loop expectation is a sum over planar
surfaces weighted by products of **signed Catalan numbers**, obtained by a peeling exploration of
the master loop equation, in any dimension. Explicit computations for general loops. Again a
combinatorial theorem, **not** a gap or area-law result. Confirms that the surface program is
still at the level of exact bookkeeping in this era.

**13.07 [O] VERIFIED** — S. Cao, S. Sheffield, "Fractional Gaussian forms and gauge theory: an
overview," arXiv:2406.19321 [math.PR] (Jun 2024).
Programmatic survey: short-distance behaviour of gauge theories modelled by fractional Gaussian
forms, with the 3D construction of 13.03 as the worked example. Best single statement of the
"GFF-input" research programme and its open endpoints. No new theorem.

**13.08 [T] VERIFIED** — R. Bauerschmidt, B. Dagallier, "Log-Sobolev inequality for the φ⁴₂ and
φ⁴₃ measures," *Comm. Pure Appl. Math.* **77** (2024), 2579–2612; arXiv:2202.02295.
Not a gauge paper, but load-bearing method: a multiscale (Polchinski) proof of a volume-uniform
log-Sobolev inequality for a non-convex, non-gradient model in d=2,3. This is the technique any
d=4 gauge LSI would have to imitate; the era did not get there for SU(N).

**13.09 [T] VERIFIED** — R. Bauerschmidt, T. Bodineau, B. Dagallier, "Stochastic dynamics and the
Polchinski equation: an introduction," *Probab. Surv.* **21** (2024), 200–290; arXiv:2307.07619.
Survey of the multiscale Bakry–Émery / Polchinski-flow machinery for volume-independent LSI,
including the renormalized potential and its Hamilton–Jacobi structure. The reference to cite for
"what a cutoff-consistent LSI requires." Explicitly a framework paper; no gauge application.

**13.10 [O] UNVERIFIED** — S. Chatterjee, *Gauge Theory*, World Scientific, 2022.
Widely cited monograph on rigorous lattice gauge theory (loop equations, 1/N expansion, confinement
mechanisms). Publisher record was **not** fetched; title/author/publisher/year only. Do not quote
page numbers or ISBN from this entry until the publisher page is checked.

**13.11 [O] UNVERIFIED** — R. Bauerschmidt, B. Dagallier, "Log-Sobolev inequality for lattice gauge
theory with U(1) gauge group in d=4," reported 2023–2024 (Oberwolfach Workshop 2351 programme;
secondary citations describe it as the first d=4 gauge LSI).
**Caution:** the complete arXiv author listing for B. Dagallier (checked 2026-11-24, 16 records
through 2026) contains **no paper with this title**. It appears to be talk/Oberwolfach-abstract
level, or an unpublished manuscript, not an arXiv/journal paper. Treat as a claimed-but-unlocated
result. This is exactly the kind of item v3 must not cite as established.

---

## 2. Lattice determinations: glueball spectrum, string tension, large N, step scaling (6 entries)

**13.12 [V] VERIFIED** — A. Athenodorou, M. Teper, "SU(N) gauge theories in 3+1 dimensions:
glueball spectrum, string tensions and topology," *JHEP* **12** (2021) 082.
DOI 10.1007/JHEP12(2021)082; arXiv:2106.00364 [hep-lat].
The era's reference quenched calculation: Wilson plaquette action, **N = 2,3,…,12**, low-lying
glueball spectrum in all cubic-group representations with both parity and charge conjugation,
continuum extrapolation for each N and then to N=∞; continuum spins identified for many states;
fundamental and k=2 string tensions and their N-dependence; running of the lattice coupling
confirmed as **g²(a) ∝ 1/N** at fixed physics; Λ_MSbar/√σ from 3-loop β-function fits, with
analytic formulas for √σ at given coupling; topological charge for N=2…6 and quantitative
"freezing" of topology as a→0 and N→∞. Standard SU(3) continuum anchors it reproduces are
m_{0++} ≈ 1.7 GeV and √σ ≈ 0.42–0.44 GeV. Systematic caveats: quenched; topology freezing limits
the smallest a at large N; continuum extrapolations are per-N, so the N→∞ limit carries the
combined systematic.

**13.13 [V] VERIFIED** — B. Lucini, M. Panero, "SU(N) gauge theories at large N," *Phys. Rept.*
**526** (2013) 93–163; arXiv:1210.4997 [hep-th].
Pre-era baseline that defines the large-N lattice programme (glueball and string-tension scaling,
topological susceptibility, 1/N corrections) which 13.12 completes at N=12. Included because v3's
large-N statements rest on this review plus 13.12; it contains no mass-gap theorem.

**13.14 [V] VERIFIED** — Y. Aoki et al. (FLAG), "FLAG Review 2021," *Eur. Phys. J. C* **82** (2022)
869; arXiv:2111.09849 [hep-lat].
The lattice-QCD precision/error-budget standard (quark masses, α_s, decay constants). It does
**not** cover pure-YM glueballs, and must not be cited as a pure-YM spectrum source; its role here
is methodological — the continuum-extrapolation and error-correlation standards that a v3 lattice
claim should meet.

**13.15 [V] VERIFIED** — M. Dalla Brida, "Past, present, and future of precision determinations of
the QCD coupling from lattice QCD," *Eur. Phys. J. A* **57** (2021) 66; arXiv:2012.01232 [hep-lat].
Review of step-scaling / Schrödinger-functional determinations of α_s and Λ. The pure-gauge sector
is where the continuum Λ is best controlled, so this is the correct citation for "step scaling"
in the era. Review, no new theorem.

**13.16 [V] VERIFIED** — C. T. Peterson, A. Hasenfratz, "Infrared fixed point in the massless
twelve-flavor SU(3) gauge-fermion system," *Phys. Rev. D* **109** (2024) 114507;
arXiv:2402.18038 [hep-lat]. (Author order as on arXiv: Peterson, Hasenfratz.)
Nf=12 in the fundamental representation, nHYP-smeared staggered fermions with Pauli–Villars
improvement; no first-order bulk phase transition in the explored region; infinite-volume
gradient-flow renormalization scheme. The GF β-function has a zero at **g²_GF⋆ = 6.60(62)** and
leading irrelevant exponent **γ_g⋆ = 0.199(32)**, consistent with earlier literature at 1–2σ.
Strong numerical evidence for IR conformality — **numerical, not a theorem**, and scheme-dependent.

**13.17 [V] VERIFIED** — Z. Fodor, K. Holland, J. Kuti, S. Mondal, D. Nógrádi, C. H. Wong, "Fate of
the conformal fixed point with twelve massless fermions and SU(3) gauge group," *Phys. Rev. D*
**94** (2016) 091501(R); arXiv:1607.06121 [hep-lat].
Pre-era but essential contrast: using the **same renormalization scheme** as the earlier positive
study, finds **no β-function fixed point in the reported interval 6.0 < g² < 6.4**, eliminating
what had been the only seemingly credible IRFP evidence; the abstract states the IR properties
"remain unresolved," not that the theory is chirally broken. Read together with 13.16, the honest
status is: **two continuum-oriented lattice programmes disagree, and the discrepancy tracks
fermion formulation and scheme.** No conformal-window claim from this era should be treated as
settled.

---

## 3. Quantum simulation and tensor networks (9 entries — proposals and benchmarks)

**13.18 [O] VERIFIED** — C. W. Bauer et al., "Quantum simulation for high-energy physics," *PRX
Quantum* **4** (2023) 027001; arXiv:2204.03381 [quant-ph].
Community roadmap (69 authors) covering lattice gauge theory, parton physics, and resource
estimates. A **proposal/strategy document**, not a result. Cite for the field's consensus targets,
never as evidence about 4D SU(N).

**13.19 [O] VERIFIED (arXiv)** — M. C. Bañuls et al., "Simulating lattice gauge theories within
quantum technologies," *Eur. Phys. J. D* **74** (2020) 165; arXiv:1911.00003 [quant-ph].
Journal volume/page from standard citation, not independently fetched. The founding review of the
quantum-tech LGT programme (trapped ions, cold atoms, tensor networks). Pre-era boundary (2019/20)
but the reference every 2021–2024 proposal cites.

**13.20 [O] VERIFIED (arXiv)** — E. Zohar, "Quantum simulation of lattice gauge theories in more
than one space dimension — requirements, challenges, methods," *Phil. Trans. R. Soc. A* **380**
(2021) 20210069; arXiv:2106.04609 [quant-ph]. Journal ref not independently fetched.
Sharp statement of why 2+1D and 3+1D non-Abelian simulation is hard (gauge-invariance protection,
dynamical fermions, locality). Proposal-level.

**13.21 [O] VERIFIED** — A. Ciavarella, N. Klco, M. J. Savage, "A trailhead for quantum simulation
of SU(3) Yang–Mills lattice gauge theory in the local multiplet basis," *Phys. Rev. D* **103**
(2021) 094501; arXiv:2101.10227 [quant-ph].
Concrete SU(3) digitization and resource estimate in the local multiplet basis. Proposal; no
hardware, no continuum.

**13.22 [O] VERIFIED** — Z. Davoudi, M. Hafezi, C. Monroe, G. Pagano, A. Seif, A. Shaw, "Towards
analog quantum simulations of lattice gauge theories with trapped ions," *Phys. Rev. Research* **2**
(2020) 023015; arXiv:1908.03210 [quant-ph].
Analog trapped-ion scheme (U(1)-type). Pre-era boundary; proposal, with explicit experimental
error budgets. Not non-Abelian SU(N).

**13.23 [V] VERIFIED** — R. C. Farrell, M. Illa, A. N. Ciavarella, M. J. Savage, "Quantum
simulations of hadron dynamics in the Schwinger model using 112 qubits," *Phys. Rev. D* **109**
(2024) 114510; arXiv:2401.08044 [quant-ph].
Circuit-level simulation of hadron dynamics in 1+1D Schwinger model at 112 qubits. Numerical
(classical emulation of the quantum circuit), 1+1D, Abelian: a capability benchmark, **not**
evidence for 3+1D non-Abelian physics.

**13.24 [O] VERIFIED (arXiv)** — G. Magnifico, G. Cataldi, M. Rigobello, P. Majcen, D. Jaschke,
P. Silvi, S. Montangero, "Tensor networks for lattice gauge theories beyond one dimension: a
roadmap," arXiv:2407.03058 [quant-ph] (Jul 2024); *Commun. Phys.* **8** (2025) 322.
Roadmap for 2D/3D tensor-network LGT: gauge-invariant PEPS, entanglement structure, truncation
strategies. Proposal/roadmap, no 3+1D non-Abelian result.

**13.25 [O] VERIFIED** — A. Kelman, U. Borla, I. Gomelski, J. Elyovich, G. Roose, P. Emonts,
E. Zohar, "Gauged Gaussian PEPS — a high dimensional tensor network formulation for lattice gauge
theories," *Phys. Rev. D* **110** (2024) 054511; arXiv:2404.13123 [quant-ph].
Construction of a Gaussian-PEPS ansatz with gauging in higher dimensions. Ansatz/method paper;
benchmarks are small systems.

**13.26 [V] VERIFIED** — J. Knaute, M. Feuerstein, E. Zohar, "Entanglement and confinement in
lattice gauge theory tensor networks," *JHEP* **2024** (2024) 174; arXiv:2401.01930 [quant-ph].
Numerical tensor-network study of entanglement/confinement in low-dimensional LGT. A genuine
numerical result, but in 1+1D/2+1D and not a continuum mass-gap determination.

---

## 4. Machine learning and flow-based sampling (5 entries)

**13.27 [V] VERIFIED** — M. S. Albergo, G. Kanwar, P. E. Shanahan, "Flow-based generative models
for Markov chain Monte Carlo in lattice field theory," *Phys. Rev. D* **100** (2019) 034515;
arXiv:1904.12072 [hep-lat].
Proof of principle that normalizing flows can sample scalar lattice field theories. Pre-era
boundary; the seed of the programme. Numerical.

**13.28 [V] VERIFIED** — D. Boyda, G. Kanwar, S. Racanière, D. J. Rezende, M. S. Albergo,
K. Cranmer, D. C. Hackett, P. E. Shanahan, "Sampling using SU(N) gauge equivariant flows," *Phys.
Rev. D* **103** (2021) 074504; arXiv:2008.05456 [hep-lat].
Gauge-equivariant continuous flows for SU(N) lattice gauge theory. Establishes that exact gauge
symmetry can be built into the sampler; demonstrated at small volume / coarse coupling. Numerical.

**13.29 [V] VERIFIED** — M. S. Albergo, G. Kanwar, S. Racanière, D. J. Rezende, J. M. Urban,
D. Boyda, K. Cranmer, D. C. Hackett, P. E. Shanahan, "Flow-based sampling for fermionic lattice
field theories," *Phys. Rev. D* **104** (2021) 114507; arXiv:2106.05934 [hep-lat].
Extends flows to fermionic (pseudofermion) path integrals. Numerical; small systems.

**13.30 [V] VERIFIED (arXiv)** — R. Abbott et al., "Normalizing flows for lattice gauge theory in
arbitrary space-time dimension," arXiv:2305.02402 [hep-lat] (May 2023). No journal ref in arXiv
metadata at check time.
Scaling flow architectures to arbitrary d with a proof-of-principle **4D SU(3)** application.
Numerical; no continuum extrapolation and no error-controlled comparison to a known answer at
physical coupling. The honest description is "sampler technology," not "lattice solution."

**13.31 [V] VERIFIED** — D. Bachtis, G. Aarts, B. Lucini, "Mapping distinct phase transitions to a
neural network," *Phys. Rev. E* **102** (2020) 053306; arXiv:2007.00355 [cond-mat.stat-mech].
Demonstrates NN classification of phase transitions (including Berezinskii–Kosterlitz–Thouless)
from configurations. Numerical methodology; not gauge-specific and not a mass-gap tool.

---

## 5. Coupling range, precisely, and the distance to the physical coupling

Normalization warning: "β" is not one object across this literature. SZZ write the lattice action
with 't Hooft scaling βN; their β equals 1/λ with λ = g²N the 't Hooft coupling. The standard
Wilson normalization is β_W = 2N/g². The campaign's own tables use β = 1/g² with the physical point
near β ≈ 2.5. **These three must be reconciled once in v3**; the numbers below are given in each
paper's own convention.

| Result | Group / dimension | Coupling reached | Volume-uniform? | Gap? |
|---|---|---|---|---|
| SZZ 2023 (13.01) | SU(N), SO(N); any d>1 | SU(N), d=4: β < 1/48, i.e. λ > 48 | yes | yes, strictly positive |
| Adhikari–Cao (13.02) | finite groups; d not checked | weak coupling (β→∞) | infinite volume | exponential decay |
| Athenodorou–Teper (13.12) | SU(N), N=2…12; d=4 | continuum extrapolation per N (β_W ≈ 5.7–6.9 at N=3, scaling ∝ N² at fixed λ) | finite volume | numerical m/√σ |
| Peterson–Hasenfratz (13.16) | SU(3), Nf=12; d=4 | g²_GF up to ≈7, zero at 6.60(62) | infinite-volume scheme | n/a (conformal) |

**Distance to the physical coupling.** The physical SU(3) lattice coupling used in the continuum
extrapolations is λ = g²N ≈ 1–3 (β_W ≈ 5.8–6.5). SZZ require λ > 48. The rigorous strong-coupling
window therefore stops a factor **≈16–48 in λ** — one to two orders of magnitude — short of the
physical point. In the campaign's β = 1/g² normalization, with the physical point at β ≈ 2.5, the
shortfall is of the same order (≈10²); the exact conversion depends on which action normalization
is matched and should be fixed once. The lower endpoint of the campaign's interval [1/48, 2.5] is
literally SZZ's d=4 SU(N) bound; the literature supplies that endpoint and **nothing between it
and 2.5**.

**What is rigorous vs numerical vs proposal.**
- Rigorous, unconditional, 4D: only 13.01, and only at strong coupling.
- Rigorous, conditional: 13.03 (3D, GFF hypothesis).
- Rigorous, exact identities: 13.05, 13.06 (surface expansions; no gap).
- Rigorous, weak coupling: 13.02, finite groups only.
- Numerical: 13.12, 13.13, 13.14, 13.15, 13.16, 13.17, 13.26, 13.27–13.31.
- Proposal/roadmap: 13.18–13.22, 13.24, 13.25.
- Unlocated/claimed: 13.11.

---

## 6. Position in 14 eras

| # | Window | What this era owes it / how it stands |
|---|---|---|
| 01 | 1954–1965 | Inherits the gauge-mass question; contributes no 4D answer, only sharper formulations. |
| 02 | 1963–1973 | Asymptotic freedom explains why the physical coupling is weak — and why strong-coupling theorems (13.01) cannot reach it. |
| 03 | 1974–1979 | Instantons/Savvidy: lattice topology data (13.12) constrain the topological sector, but the Savvidy instability is untouched rigorously. |
| 04 | 1974–1982 | SZZ 2023 is the modern completion of Osterwalder–Seiler 1978 — same β = O(1) strong-coupling regime, now with a full stochastic-analysis proof. |
| 05 | 1982–1990 | Constructive RG becomes the Polchinski/multiscale-LSI framework (13.08–13.09); still no 4D gauge LSI. |
| 06 | 1985–1995 | Confinement pictures: the era's quantitative monument is the 2021 string-tension scan (13.12); no rigorous mechanism. |
| 07 | 1993–2000 | DSE/ERG largely superseded as evidence by lattice; ERG survives only inside 13.09. |
| 08 | 1998–2005 | Holography/large-N: 13.12–13.13 supply the lattice anchor (N≤12, g²∝1/N) that any holographic gap claim must match. |
| 09 | 2002–2010 | Precision/step-scaling/conformal window: continued by 13.15; the 12-flavor question is still numerically contested (13.16 vs 13.17). |
| 10 | 2008–2014 | CFT constraints: no contact — bootstrap does not constrain 3+1D pure YM. |
| 11 | 2013–2018 | Direct ancestor of 13.02–13.04 and 13.01 (stochastic analysis, regularity structures, Langevin). |
| 12 | 2018–2021 | Resurgence/transseries: no contact; remains formal, no lattice or rigorous bridge in this era. |
| 13 | 2021–2024 | This era: one unconditional strong-coupling gap theorem, one weak-coupling finite-group theorem, no continuous-SU(N) weak-coupling result, consolidated lattice numbers, and a wave of proposals. |
| 14 | 2024–2026 | Hands over SZZ as the theorem to formalize, and the 2025 CPS area-law / Nissim 't Hooft-regime programme as the immediate successor — *outside* this era. |

---

## 7. Gaps and attack surfaces for v3

1. **The coupling gap is the whole problem.** SZZ stops at β < 1/48 (λ > 48); the physical point is
   λ ≈ 1–3. No theorem in 2021–2024 crosses it for continuous SU(N). The campaign's
   cutoff-consistent G″(0)>0 is precisely the missing bridge; the literature confirms the bridge
   does not exist, it does not suggest it is easy.
2. **Strong/weak split by group type.** The only weak-coupling decay theorem (13.02) is for finite
   groups. Extending it to SU(N) would require controlling the non-compact direction of the gauge
   group at weak coupling — nothing in this era does. This is a concrete, citable open step.
3. **No N-uniform gap constant.** SZZ's hypothesis is N-uniform for SU(N), but no explicit
   N-dependence of the gap lower bound is stated. If the constant degrades exponentially in N, the
   large-N route to the continuum is blocked; if it does not, that is a publishable lemma. Check the
   LSI constant's N-scaling in the proof.
4. **d=4 LSI status is unresolved.** The alleged Bauerschmidt–Dagallier U(1) d=4 LSI (13.11) has no
   locatable paper; the verified d=4-adjacent LSI results are φ⁴₂/φ⁴₃ (13.08) and the strong-coupling
   gauge LSI inside 13.01. v3 should either locate the U(1) paper or state the gap as open.
5. **Lattice systematics, not statistics, are the limit.** Athenodorou–Teper's dominant uncertainties
   are the a→0 range and topological freezing at large N — precisely the regime where a continuum
   mass-gap-in-units-of-Λ claim would be tested. A v3 lattice claim needs an explicit freezing
   budget.
6. **The 12-flavor contradiction is a warning, not a tool.** LSD finds no IRFP in 6.0<g²<6.4 in one
   scheme; Hasenfratz–Peterson find one at 6.60(62) in another. Any v3 argument that leans on
   "conformal window" or "walking" must state which scheme and why the other is excluded.
7. **Quantum simulation is pre-evidential.** All 2021–2024 entries are proposals or small-system
   circuit emulations (1+1D Abelian for the largest, 13.23). None constrains 3+1D non-Abelian
   physics. Cite as motivation only.
8. **Tensor networks have not reached 3+1D non-Abelian.** The 2024 roadmap (13.24) and PEPS papers
   (13.25–13.26) are 1+1D/2+1D. A 3+1D continuum extrapolation is not on the horizon in this
   literature.
9. **Flow-based sampling is a sampler, not a solution.** 4D SU(3) demonstrations (13.28, 13.30) are
   small-volume and coarse-coupling, with statistical errors only. No flow paper in the era reports
   a continuum extrapolated observable at physical coupling.
10. **The successor programme is the attack surface.** CPS 2025 (area law in the 't Hooft regime,
    arXiv:2505.16585, 2509.04688) and Nissim 2025 (U(N) 't Hooft-regime mass gap, arXiv:2510.22788)
    are the first results that combine the SZZ dynamic with the large-N/t Hooft regime. v3's
    "machine" should be positioned relative to them, not to SZZ alone — and must note they are
    2025, hence Era 14, not this era.

---

## 8. Verification ledger

| # | Entry (short) | Identifier checked | Checked against | Status | Note |
|---|---|---|---|---|---|
| 13.01 | Shen–Zhu–Zhu | arXiv:2204.12737; CMP 400 (2023) 805–851 | arXiv abs page + zbMATH Zbl 1522.81127 snippet | VERIFIED | abstract quoted; β<1/48 confirmed |
| 13.02 | Adhikari–Cao | arXiv:2202.10375; Ann. Probab. 53 (2025) 140–174 | arXiv API + search snippet | VERIFIED | dimension not checked |
| 13.03 | Cao–Chatterjee | arXiv:2111.12813; CMP 405 (2024) Paper 3 | arXiv API + MR4689958 snippet | VERIFIED | conditional result |
| 13.04 | Chandra–Chevyrev–Hairer–Shen | arXiv:2201.03487; Invent. Math. 237 (2024) 541–696 | arXiv API (journal ref in feed) | VERIFIED | 3D only |
| 13.05 | Cao–Park–Sheffield | arXiv:2307.06790 | arXiv abs page | VERIFIED | Comm. AMS "to appear" |
| 13.06 | Borga–Cao–Shogren-Knaak | arXiv:2411.11676 | arXiv abs page | VERIFIED | no venue yet |
| 13.07 | Cao–Sheffield overview | arXiv:2406.19321 | arXiv API | VERIFIED | survey |
| 13.08 | Bauerschmidt–Dagallier φ⁴ | arXiv:2202.02295; CPAM 77 (2024) 2579–2612 | arXiv API | VERIFIED | method paper |
| 13.09 | Bauerschmidt–Bodineau–Dagallier | arXiv:2307.07619; Probab. Surv. 21 (2024) 200–290 | arXiv API | VERIFIED | method paper |
| 13.10 | Chatterjee, *Gauge Theory* | World Scientific 2022 | none | UNVERIFIED | publisher not fetched |
| 13.11 | Bauerschmidt–Dagallier U(1) d=4 LSI | no arXiv id | full arXiv author listing (16 records, through 2026) | UNVERIFIED | no such paper located; talk-level |
| 13.12 | Athenodorou–Teper | arXiv:2106.00364; JHEP 12 (2021) 082 | arXiv API + abs page + search snippet | VERIFIED | N=2…12 confirmed |
| 13.13 | Lucini–Panero | arXiv:1210.4997; Phys. Rept. 526 (2013) 93–163 | arXiv API | VERIFIED | pre-era baseline |
| 13.14 | FLAG Review 2021 | arXiv:2111.09849; EPJC 82 (2022) 869 | arXiv API | VERIFIED | QCD, not pure YM |
| 13.15 | Dalla Brida review | arXiv:2012.01232; EPJA 57 (2021) 66 | arXiv API | VERIFIED | review |
| 13.16 | Peterson–Hasenfratz | arXiv:2402.18038; PRD 109 (2024) 114507 | arXiv API + abs page + search snippet | VERIFIED | author order per arXiv |
| 13.17 | Fodor et al. | arXiv:1607.06121; PRD 94 (2016) 091501 | arXiv API | VERIFIED | corrected from 2017/96 |
| 13.18 | Bauer et al. | arXiv:2204.03381; PRX Quantum 4 (2023) 027001 | arXiv API | VERIFIED | roadmap |
| 13.19 | Bañuls et al. | arXiv:1911.00003; EPJD 74 (2020) 165 | arXiv API | VERIFIED (arXiv) | journal ref not fetched |
| 13.20 | Zohar | arXiv:2106.04609; Phil. Trans. R. Soc. A 380 (2021) | arXiv API | VERIFIED (arXiv) | journal ref not fetched |
| 13.21 | Ciavarella–Klco–Savage | arXiv:2101.10227; PRD 103 (2021) 094501 | arXiv API | VERIFIED | proposal |
| 13.22 | Davoudi et al. | arXiv:1908.03210; PRResearch 2 (2020) 023015 | arXiv API | VERIFIED | proposal |
| 13.23 | Farrell et al. | arXiv:2401.08044; PRD 109 (2024) 114510 | arXiv API | VERIFIED | 1+1D benchmark |
| 13.24 | Magnifico et al. | arXiv:2407.03058; Commun. Phys. 8 (2025) 322 | arXiv API | VERIFIED (arXiv) | journal 2025 |
| 13.25 | Kelman et al. | arXiv:2404.13123; PRD 110 (2024) 054511 | arXiv API | VERIFIED | ansatz paper |
| 13.26 | Knaute–Feuerstein–Zohar | arXiv:2401.01930; JHEP 2024 (2024) 174 | arXiv API | VERIFIED | numerical |
| 13.27 | Albergo–Kanwar–Shanahan | arXiv:1904.12072; PRD 100 (2019) 034515 | arXiv API | VERIFIED | pre-era seed |
| 13.28 | Boyda et al. | arXiv:2008.05456; PRD 103 (2021) 074504 | arXiv API | VERIFIED | numerical |
| 13.29 | Albergo et al. (fermionic) | arXiv:2106.05934; PRD 104 (2021) 114507 | arXiv API | VERIFIED | numerical |
| 13.30 | Abbott et al. | arXiv:2305.02402 | arXiv API + title query | VERIFIED (arXiv) | no journal ref |
| 13.31 | Bachtis–Aarts–Lucini | arXiv:2007.00355; PRE 102 (2020) 053306 | arXiv API | VERIFIED | corrected venue from PRL guess |

**Counts:** 31 entries; 29 VERIFIED (of which 6 are arXiv-verified with journal field not
independently fetched), 2 UNVERIFIED (13.10 book, 13.11 unlocated LSI claim).

**Corrections made during checking (record these, do not propagate the wrong forms):**
- Fodor et al. conformal-fixed-point paper is **PRD 94 (2016) 091501**, arXiv:1607.06121 — not
  PRD 96 (2017).
- The large-N review "SU(N) gauge theories at large N" is **Lucini–Panero**, not Lucini–Teper.
- The NN phase-classification paper is **Phys. Rev. E 102 (2020) 053306**, not PRL.
- "Normalizing flows for lattice gauge theory in arbitrary space-time dimension" is
  **Abbott et al. 2023 (arXiv:2305.02402)**, not the 2020 Kanwar et al. PRL.
- The Bauerschmidt–Dagallier U(1) d=4 lattice-gauge LSI **could not be located** in the authors'
  complete arXiv listing; cite only with the caveat in 13.11.
