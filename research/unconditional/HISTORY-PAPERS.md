# HISTORY-PAPERS — verified history and key papers of the Clay Yang–Mills mass gap

> 2026-11-23 · claims checked against primary sources (CMI page, journal/DOI pages, arXiv abstracts). Scope = the six items + timeline + citation audit.

> **判词 (verdict).** Still **open**; CMI labels it "Unsolved" ("no proof of this property is known"). No claimed proof survived peer review. The frontier is asymmetric: the **strong-coupling *lattice* gap and area law are theorems** (Osterwalder–Seiler 1978; Shen–Zhu–Zhu 2023; Cao–Nissim–Sheffield 2025), while the **continuum 4D construction** (OS axioms + infinite volume + Δ>0) is unbuilt. Balaban's UV stability is *not* a continuum limit; 2025–26 "completions" are unverified preprints.

## §1 · Clay statement
A. Jaffe & E. Witten, "Quantum Yang–Mills Theory," official CMI Millennium problem description, **2000** (reprinted in *The Millennium Prize Problems*, eds. Carlson–Jaffe–Wiles, AMS/CMI 2006). Ask: for every compact simple gauge group G, construct a non-trivial quantum YM theory on ℝ⁴ whose existence meets axioms "at least as strong as" Streater–Wightman (1964) and Osterwalder–Schrader (1973, 1975), **and** prove a mass gap Δ>0. CMI page still "Unsolved."

## §2 · Foundations of rigorous Euclidean QFT
- **Osterwalder–Schrader axioms:** K. Osterwalder & R. Schrader, "Axioms for Euclidean Green's Functions," *CMP* **31** (1973) 83–112; II, *CMP* **42** (1975) 281–305 — analyticity, symmetry, reflection positivity, cluster; OS reconstruction ⇒ Wightman theory.
- **Glimm–Jaffe:** *Quantum Physics: A Functional Integral Point of View*, Springer **1981** (2nd ed. 1987) — standard CQFT text; 2D/3D existence theorems and OS reconstruction.
- **Guerra–Rosen–Simon:** "Nelson's symmetry and the infinite-volume behavior of the vacuum in P(φ)₂," *CMP* **27** (1972) 10–22; "The P(φ)₂ Euclidean QFT as classical statistical mechanics," *Ann. of Math.* **101** (1975) 111–259.

## §3 · 4D constructive programs
- **Balaban (1980s–90s):** "UV stability of three-dimensional lattice pure gauge field theories," *CMP* **102** (1985) 255–275; RG series *CMP* **109** (1987) 249–301, **119** (1989) 243; "Large field renormalization I/II," *CMP* **122** (1989) 175–202 & 355–392. Established **UV stability of the finite-volume lattice effective action**; did **not** take the continuum limit, verify OS data, or prove a gap.
- **Magnen–Rivasseau–Sénéor:** "Construction of YM₄ with an infrared cutoff," *CMP* **155** (1993) 325–383 — pure SU(2), trivial topological sector, axial gauge, IR cutoff, no UV cutoff; verified Slavnov identities non-perturbatively. Foundations only; no gap.
- **Feldman–Osterwalder:** the landmark is "The Wightman axioms and the mass gap for weakly coupled (φ⁴)₃ QFT," *Ann. Physics* **97** (1976) 80–135 — **3D φ⁴, not 4D gauge**. 4D φ⁴ is Feldman–Magnen–Rivasseau–Sénéor, *CMP* **109** (1987) 437–480.

## §4 · Strong-coupling lattice theorems
- **Osterwalder–Seiler:** "Gauge field theories on a lattice," *Ann. Physics* **110** (1978) 440–471 — reflection positivity, positive transfer matrix, infinite-volume strong-coupling limit, Wilson confinement bound.
- **Brydges–Fröhlich–Seiler:** I. *Ann. Physics* **121** (1979) 227–284; II. *CMP* **71** (1980) 159–205; III. *CMP* **79** (1981) 353–399 — diamagnetic inequalities, mass generation, 2D Abelian–Higgs without cutoffs, OS axioms for gauge-invariant fields.

## §5 · Modern (2015–2025)
- **Chatterjee large-N:** "Rigorous solution of strongly coupled SO(N) lattice gauge theory in the large-N limit," *CMP* **366** (2019) 203–268 (arXiv:1502.07719, 2015) — gauge-string duality. Survey "Yang–Mills for probabilists" (arXiv:1803.01950, 2019); "A probabilistic mechanism for quark confinement," *CMP* **385** (2021) 1007–1039.
- **Shen–Zhu–Zhu:** "A stochastic analysis approach to lattice Yang–Mills at strong coupling," *CMP* **400** (2023) 805–851 (arXiv:2204.12737) — strong coupling |β| < 1/(16(d−1)) for SU(N), **uniform in N**; unique invariant measure, exponential ergodicity, LSI/Poincaré, correlation decay ⇒ **lattice mass gap**.
- **Cao–Nissim–Sheffield ("CNS"):** "Expanded regimes of area law for lattice Yang–Mills," arXiv:2505.16585 (May 2025; improves Osterwalder–Seiler 1978); "Dynamical approach to area law," arXiv:2509.04688 (Sep 2025; 't Hooft regime, U(N)/SU(N)/SO(2N)).
- **Nissim:** "U(N) lattice Yang–Mills in the 't Hooft regime," arXiv:2510.22788 (Oct 2025) — mass gap, unique infinite-volume limit, new large-N proof for U(N). **Lattice only.**
- **Context:** Chandra–Chevyrev–Hairer–Shen, "Stochastic quantisation of Yang–Mills–Higgs in 3D," *Invent. Math.* **237** (2024) 541–696 (arXiv:2201.03487); Bauerschmidt–Bodineau–Dagallier, "Stochastic dynamics and the Polchinski equation," *Probab. Surveys* (2024).

## §6 · Claimed / rumored proofs — status
All **unaccepted**: (i) Zenodo/OSF/philarchive 2025–26 "proofs" (e.g. *Spectral Operator with a Positive Gap*; *Balaban–Tamesis closure*; *Formal Proof … Synthesizing Balaban's UV Stability*) are non-peer-reviewed and often misattribute Balaban's UV-stability as a gap proof; (ii) arXiv:2506.00284 (Jacobson, 5D-orbifold SU(3)) was **withdrawn**; (iii) the 2025 OpenAI/Lean "Humanise Yang–Mills" formalization explicitly places the mass gap **out of scope**. None is a Clay solution.

## §7 · Timeline (what each established)
- **1964** Streater–Wightman axioms · **1972–75** GRS: Nelson symmetry, P(φ)₂ · **1973/75** OS axioms · **1976** Feldman–Osterwalder: φ⁴₃ mass gap (3D) · **1978** Osterwalder–Seiler: lattice strong-coupling confinement · **1979–81** BFS: gauge-field construction · **1981/87** Glimm–Jaffe book · **1985–89** Balaban: 3D/4D UV stability (no continuum limit) · **1987** FMRS: infrared φ⁴₄ · **1993** MRS: YM₄ with IR cutoff · **2000** Clay problem posed · **2015/19** Chatterjee: large-N SO(N) · **2021** Chatterjee: confinement mechanism · **2023** SZZ: lattice strong-coupling gap · **2024** CCHS 3D YMH; BBD Polchinski · **2025** CNS area law; Nissim U(N) 't Hooft gap · **2026** still open.

## §8 · Citation errors / overclaims flagged
1. **"Osterwalder–Seiler construction axioms"** (item 2) is a misnomer: the axioms are **Osterwalder–Schrader** (1973/75); Osterwalder–Seiler (1978) is a lattice *construction*, listed correctly in item 4.
2. **Feldman–Osterwalder** grouped under "4D constructive" is wrong: their 1976 mass-gap theorem is **(φ⁴)₃, not 4D gauge**. 4D gauge = Balaban/MRS; 4D φ⁴ = FMRS.
3. **Shen–Zhu–Zhu** *CMP* **400** (2023) 805–851 is correct, but it is a **strong-coupling lattice** gap, not the continuum/large-N Clay statement.
4. **Balaban "4D UV stability" ≠ continuum limit or mass gap** — treat 2025–26 "completion" claims as unverified.
5. **"CPS"** could not be verified; the real 2025 trio is **CNS = Cao–Nissim–Sheffield**. The SZZ threshold |β|<1/(16(d−1)) is the concrete strong-coupling bound (a 2026 preprint's β<195/7168 is unrefereed).
6. Every Zenodo/OSF/philarchive "solution" is **not** a proof.

*Sources: CMI Yang–Mills page; Springer/DOI and Project Euclid journal pages; arXiv abstracts 2204.12737, 2505.16585, 2509.04688, 2510.22788, 2201.03487, 1502.07719, 1803.01950, 2506.00284.*
