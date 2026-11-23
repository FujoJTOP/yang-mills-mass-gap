# CONSTRUCTION-METHODS — rigorous routes to 4D pure SU(N) Yang–Mills + mass gap

> Each method: mechanism / exact theorem / 4D stall / concrete extension. Proved, stalled, and missing are kept separate.

## (a) Bałaban multiscale large-field/small-field RG (CMP 1984–89)
**How.** Block-spin RG. Each scale: small-field Gaussian expansion about a background (effective action + coupling renormalization) and large-field R-operation (localization, exponentiation, exponential suppression). Finite volume, background/axial gauge, SU(2).
**Proves.** Convergence of the RG transformations and a cutoff-uniform UV-stability bound, with coupling renormalization and bounded large fields (Bałaban, *CMP* 95:17; 96:223; 98:17; 109:249; 116:1; 119:243; 122:175,355). Dimock (*Rev. Math. Phys.* 25:1330010; *J. Math. Phys.* 54:092301; *AHP* 15:2133) proves the φ⁴₃ analogue.
**Stalls.** No gauge-invariant continuum measure, Schwinger functions, OS positivity, or gap. Large-field control is UV-adapted; in the IR the coupling is O(1) and large fields are *typical*. The bounds control the action's value/norm, not the sign of its second variation.
**Extension.** Uniform all-scale d=4 large-field estimate; gauge-invariant OS reconstruction of the J→∞ limit; transfer-matrix gap bound.

## (b) Magnen–Rivasseau–Sénéor (*CMP* 155, 1993)
**How.** Constructive cluster expansion with the UV cutoff removed and a fixed IR cutoff; regularized axial gauge (positive large-field weight) plus a perturbative gauge in a lower-momentum background; nonperturbative Slavnov identities restore the broken gauge invariance. Trivial sector, pure SU(2).
**Proves.** UV-renormalized YM₄ Schwinger functions with an IR cutoff, in regularized axial gauge, obeying nonperturbative Slavnov–Taylor identities (*CMP* 155:325–383).
**Stalls.** IR cutoff never removed; trivial sector only; no OS positivity in the limit; no gap.
**Extension.** Infinite-volume/IR limit, all topological sectors, gauge-invariant OS positivity, clustering.

## (c) Rigorous RG / Polchinski–Wetterich
**How.** Exact cutoff flow for the effective action: Polchinski's diffusion equation (*Nucl. Phys.* B231:269) and Wetterich's 1PI form (*Phys. Lett.* B301:90); gauge versions use modified Ward identities or a gauge-invariant regulator (Arnone–Morris–Rosten, arXiv:hep-th/0506162).
**Proves.** A rigorous nonperturbative Wilsonian construction for φ⁴₄ — a trapping region about the Gaussian fixed point with marginal-log resummation (Bauerschmidt–Brydges–Slade, *LNM* 2242, 2019). For gauge: perturbative renormalizability and β-functions only; no nonperturbative 4D construction; Wetterich truncations carry no convergence control.
**Stalls.** Gauge invariance/Gribov; the non-Gaussian *marginal* large-field term; the Savvidy sign.
**Extension.** Gauge-covariant all-orders flow with a trapping region uniform in cutoff and volume, then OS and gap.

## (d) Stochastic quantization + regularity structures
**How.** Parisi–Wu Langevin flow ∂ₜA = −δS/δA + ξ; super-renormalizable in low d and renormalized by regularity structures/BPHZ. Gauge covariance of the deterministic flow gives a quotient by gauge orbits; symmetry fixes the counterterms.
**Proves.** Renormalized *local* solutions and gauge-covariant Markov processes on gauge orbits in 2D (Chandra–Chevyrev–Hairer–Shen, *Publ. Math. IHÉS* 136:1, arXiv:2006.04987) and 3D (CCHS, *Invent. Math.* 237:541, arXiv:2201.03487); the 2D invariant measure (Chevyrev–Shen, arXiv:2302.12160); unique 3D gauge-covariant renormalization (Chevyrev–Shen, *ARMA* 250, arXiv:2503.03060). Framework: Hairer, *Invent. Math.* 198:269.
**Stalls.** 4D is critical/marginal ([g]=0); the subcritical degree condition fails and there is no gauge-covariant BPHZ. The analogous 4D φ⁴ is solved only as *trivial* (Aizenman–Duminil-Copin, arXiv:1912.07973), so well-posedness would not yield a nontrivial measure. 3D global existence/invariant measure is already open; the Langevin gap is fictitious-time.
**Extension.** A critical 4D regularity structure, global existence, invariant measure, physical gap — beyond current technology.

## (e) Shen–Zhu–Zhu (*CMP* 400, 2023)
**How.** Langevin dynamics for lattice YM under 't Hooft scaling βN; unique invariant measure, Wasserstein ergodicity, volume convergence, log-Sobolev/Poincaré; these give exponential correlation decay (= gap); Wilson loops factorize at large N.
**Proves.** An unconditional, volume-uniform mass gap for the infinite-volume lattice measure at strong coupling |β| < 1/(16(d−1)) for SU(N) — i.e. β < 1/48 in d=4 — under 't Hooft scaling, with factorized large-N Wilson loops (*CMP* 400:805–851, arXiv:2204.12737). The only proven 4D gauge uniform gap.
**Stalls.** Strong coupling only; the continuum needs β→∞; the constants do not reach the physical crossover β≈2.5, where Bakry–Émery fails on the Savvidy region.
**Extension.** A log-Sobolev/Poincaré bound across the crossover, uniform in a and volume.

## (f) Chatterjee rigorous large-N lattice
**How.** Master loop (Makeenko–Migdal) equations valid at all β; solve by a convergent lattice-string expansion at strong coupling; derive the 1/N expansion, factorization, and gauge–string duality.
**Proves.** A rigorous solution of strongly coupled SO(N) (and SU(N), Jafarov arXiv:1610.03821) lattice gauge theory in any d, large N: Wilson loops are an absolutely convergent sum over lattice strings (Chatterjee, *CMP* 366:203–268, arXiv:1502.07719; Chatterjee–Jafarov, arXiv:1604.04777). The scaling limit it constructs (any d≥2, including 4D) is YM–Higgs → a massive Gaussian Proca field (Chatterjee, arXiv:2401.10507); the non-Gaussian limit is open.
**Stalls.** Small fixed β only; the weak-coupling/continuum master field is unbuilt; loop equations plus reflection positivity do not force decay; the asymptotic 1/N series gives no N=3 bound; the YM–Higgs limit is Gaussian. No gap is claimed.
**Extension.** A weak-coupling large-N construction with exponential decay, uniform in a,L and stable at N=3.

## Verdict — the most promising route and its precise blocker
The most promising *construction* route is **(a) completed by (c)**: multiscale RG with a Polchinski/BBS trapping region. It alone is built for the actual 4D asymptotically-free theory, has a rigorous 4D *marginal*-coupling template (BBS φ⁴), and can in principle reach physical coupling rather than stopping at strong coupling. The strong-coupling gap (e) is a genuine theorem but confined to β<1/48, with the same obstruction at the crossover; (d) fails structurally at d=4 criticality; (b) is UV-only; (f) is strong-coupling and Gaussian in its one continuum limit.

**Precise blocker:** a uniform (all scales, volumes, cutoffs) large-field estimate in the Savvidy background — a lower bound on the second variation of the effective action,
`G″_j = E[H] − Var(G) + W ≥ δ > 0`, uniform in j,a,L.
Value/coercivity/norm bounds cannot supply it: the Savvidy Hessian `−D² + 2 ad F_B` is indefinite (lowest mode `−2g|F_B|` on `|k| < √(2g|F_B|)`). This is exactly where Bałaban's R-operation (value bounds), BBS trapping (works for φ⁴ because `+12λφ²` makes the large field *more* convex), and regularity structures (critical degree) stop. The mass gap is the IR face of this one non-abelian sign — the same object behind the Weingarten/`6j` sign and the `⟨¼trF²⟩` condensate.
