# SURVEY-CONSTRUCTIVE — Can any modern constructive / rigorous-QFT method cross the 4D Yang–Mills crossover?

> **2026-11-24.** Companion to `SURVEY-RECENT-ARXIV.md` (claim-level 2025–26 sweep), `SURVEY-AI-FORMAL.md` (Lean / AI / verified numerics), `SURVEY-ADJACENT-TOOLS.md` (other fields), `SURVEY-WILDCARD-DOMAINS.md` (algebra/topology/arithmetic). This document maps the **method classes** of constructive / rigorous QFT and answers (a)–(e) of the brief. It does **not** re-report the excluded items (Balaban's multiscale; Magnen–Rivasseau–Sénéor; Chandra–Chevyrev–Hairer–Shen 3D YM–Higgs; Bringmann–Cao–Hairer–Zhao 2D YM–Higgs; CNS area law; the 2D stochastic-quantisation line); those appear only as status markers.
>
> **Labels (charter §2).** **[T]** true theorem (complete proof; "(preprint)" marks a complete proof not yet peer-reviewed) · **[C]** conditional · **[O]** exact missing step / open · **[NG]** no-go for a named class · **[R]** reported / unverified / unrefereed. **[V]** = citation and content verified this session against arXiv / DOI / INSPIRE.
>
> **Coupling normalisation.** β is the project's inverse coupling: the SZZ SU(N) d=4 strong-coupling radius is β<1/48; the physical crossover is β≈2.5 (ratio ≈120). Under a common rescaling of β the ratio is invariant. The standard Wilson β_W = 2N/g² differs by an O(N) factor; no known result moves the threshold by more than an O(1) factor.

---

## 0. TOP-LINE (read this first)

**No.** No constructive or rigorous-QFT method in the literature reaches G2 — `spec(H) ⊆ {0} ∪ [Δ,∞)`, Δ>0, with `f₁² = |⟨0|trF²|G⟩|² > 0` — for 4D pure gauge across the crossover (`α_s ~ O(1)`, project β≈2.5), with an `(L,a)`-uniform lower bound.

Precisely:

1. **No 4D nonabelian continuum YM measure has been constructed.** G1 itself is open. Beyond 2D, even the measure is open: 2D has a global invariant measure; 3D has only local-in-time solutions / a Markov process up to possible blow-up; 4D has nothing. Chevyrev's review states verbatim: *"The cases d = 2, 3 are considered substantially simpler than d = 4 as they correspond to super-renormalisable theories … (vs. renormalisable for d = 4 …)"*, and *"a construction of the 4D YM measure … remains open."* **[V]** (Chevyrev, *J. Math. Phys.* **63** (2022), arXiv:2202.13359).

2. **Every 4D G2 theorem is regime-limited.** There are exactly three classes in the literature:
   - **Strong coupling, continuous groups, all volumes:** Osterwalder–Seiler 1978; Shen–Zhu–Zhu (SZZ) 2023 — SU(N), any d>1, `|β| < 1/(16(d−1))`, d=4 ⇒ **β<1/48**; Nissim 2025 — U(N), 't Hooft / small β, all d≥2. **[T]**
   - **Weak coupling, 4D, but only finite gauge groups:** Adhikari–Cao 2025 (*Ann. Probab.*) — volume-uniform exponential clustering for `β ≥ (11/4 + 4 log|G|)/Δ_G`. **[T]** Closest in coupling; unusable for SU(N) because `Δ_G → 0`.
   - **All couplings, but only 2D:** Piazza–Shen 2026 — log-Sobolev for 2D SU(N) at every β, volume-uniform. **[T]** The only all-coupling LSI for SU(N) gauge found in this sweep.

3. **The coupling gap is ≥ O(10²)** (1/48 vs 2.47 in the project's normalisation). Unrefereed "improvements" move the threshold by O(1) factors in their own normalisation (e.g. one claims ×1.31 over 1/48); none reaches the physical crossover. **[R]**

4. **The barrier is structural, not quantitative.** Three independent method classes fail for three precise reasons at the same physical point:
   - *Strong-coupling / high-temperature expansions*: finite radius (Pringsheim; 1/48), no known resummation.
   - *Weak-coupling / low-temperature (Pirogov–Sinai) expansions*: need an isolated ordered ground state; exist for finite groups (`Δ_G>0`), fail for continuous SU(N) because `Δ_G=0` (massless gluons at tree level).
   - *SPDE / regularity structures*: need subcriticality; 4D YM is renormalisable / marginal (the same level as 4D φ⁴), not super-renormalisable. **[V]**

5. **The abelian no-go kills any general "large β ⇒ gap" theorem.** 4D U(1) at weak coupling has a nonconfining, massless Coulomb phase (Guth 1980; Fröhlich–Spencer 1982). Hence the nonabelian gap at weak coupling is a nonperturbative / confinement phenomenon, not a consequence of small coupling. **[T]**, and a **[NG]** for the whole "large β ⇒ gap" method class.

6. **Acid test.** None of the successful methods is a candidate "crossing tool": either it does not apply to the free massless field (finite groups), or it is regime-limited, or its weak-coupling continuum limit is itself a *free massless* field (3D abelian). No method gives a sign-controlling lower bound at β≈2.5 that fails on the free field. **[O]**

**Verdict in one line:** the constructive literature has sharpened both ends of the interval `[1/48, 2.5]` but has not built the bridge; every method stops at its natural radius, and the missing object is exactly the project's H1 / `f₁²>0` transport.

---

## 1. Direct answers to (a)–(e)

### (a) 4D stochastic quantization — has anyone globally constructed 4D YM?

**No. Only 2D globally; 3D locally; 4D not at all.**

- **2D:** global invariant measure; the YM measure is invariant for the renormalised Langevin dynamic and universal (Chevyrev–Shen, *Comm. Pure Appl. Math.* (2026), DOI 10.1002/cpa.70043). **[T]**
- **3D:** local-in-time solutions of the renormalised stochastic YM–Higgs flow, a unique gauge-covariant choice of counterterms, and a canonical Markov process on gauge orbits **"up to a potential finite time blow-up"** (Chandra–Chevyrev–Hairer–Shen, *Invent. Math.* **237** (2024) 541–696). No invariant measure, no global solution. **[T]** (excluded from re-reporting; status only).
- **4D:** nothing. The obstruction is not technical: d=4 is the **renormalisable / marginal** dimension, so the small parameter (subcriticality) that drives the regularity-structure/paracontrolled machinery is absent. The same is true of 4D φ⁴; there is no global stochastic-quantisation construction at this critical level. **[V]** The gauge SPDE line (Hairer, Chandra, Chevyrev, Mattingly, Bringmann) has no 4D gauge result.

### (b) Regularity structures / paracontrolled / SPDE in the "subcritical regime"

**There is no such regime for 4D gauge.** 4D YM is not subcritical; the phrase in the brief has no referent in the literature. What exists:

- **2D:** paracontrolled/regularity-structure treatments (global; invariant measure).
- **3D:** local solutions and renormalisation; uniqueness of the gauge-covariant mass renormalisation (Chevyrev–Shen, *Arch. Ration. Mech. Anal.* **250** (2026), Paper No. 11, arXiv:2503.03060); a tree-free regularity-structure model (Sevostyanov, arXiv:2605.14616 (2026), preprint). No measure, no global solution.
- **3D abelian at weak coupling:** the only **weak-coupling continuum gauge limit** in this line. Chevyrev–Qu–Shen (arXiv:2608.27828 (2026)) prove the rescaled U(1) phase field converges in probability to the solution of the **stochastic heat equation** — a *free, massless* Gaussian limit, universal in the plaquette action. **[T]**

This is the acid test in action: the one weak-coupling continuum gauge construction in the SPDE line produces the free field, i.e. exactly the object on which the project's tools must return zero.

### (c) Variational / spectral (log-Sobolev, Bakry–Émery, spectral gap)

- **4D SU(N), strong coupling:** SZZ (*Comm. Math. Phys.* **400** (2023) 805–851; arXiv:2204.12737) prove, for any d>1, `|β| < 1/(16(d−1))` for SU(N) [d=4 ⇒ 1/48] and `|β| < (N−2)/(32(d−1)N)` for SO(N): a unique infinite-volume invariant measure, exponential ergodicity in a Wasserstein distance, **log-Sobolev and Poincaré inequalities uniform in volume**, and exponential decay of correlations — "a strictly positive mass gap". **[T]**
- **4D U(N), 't Hooft / small β:** Nissim (arXiv:2510.22788 (2025)) proves a mass gap, a unique infinite-volume limit, and a new proof of the large-N limit for all d≥2 (lattice, small β, constants independent of L; intermediate estimates use N>8π, so read it as a large-N/'t Hooft result). **[T]**
- **4D YM–Higgs, small (β,κ):** SZZ (arXiv:2401.13299 (2024; v3 2025)) prove exponential ergodicity and a positive mass gap, but not pure YM and not weak coupling. **[T]**
- **All couplings, 2D SU(N):** Piazza–Shen (arXiv:2610.35286 (2026)) prove the LSI at **arbitrary fixed β**, with a constant uniform in the lattice volume, via a two-scale disintegration argument. **[T]** No 3D or 4D all-coupling LSI was found in this sweep; the only >2D results at large β are for finite groups (Adhikari–Cao, §3.2) or fixed-cutoff / conditional criteria (Mihajlovic, Gauvin — see `SURVEY-RECENT-ARXIV`).
- **Bakry–Émery in 4D at β≈2.5:** the Γ₂ curvature input is positive only while β is small; at β≈2.5 the Savvidy Hessian `−D² + 2 ad F_B` is indefinite, so Γ₂ goes negative and the method returns nothing. The "Dobrushin ceiling" is the same statement in single-link form. This is a **method** failure, not a proof that the LSI fails. **[O]/[NG for the method class]**
- **No spectral-gap lower bound for continuous SU(N) at weak coupling was found in any dimension >2.**

### (d) Cluster / polymer / lace expansion for gauge

- **Strong-coupling cluster/polymer expansions:** 4D continuous groups, small β, gap. Osterwalder–Seiler 1978 (reflection positivity, positive transfer matrix, existence/analyticity of the infinite-volume limit, Wilson's area law); Seiler's monograph (Lecture Notes in Physics 159, Springer 1982) for the cluster-expansion form; SZZ 2023 as the modern LSI version. **[T]**
- **Weak-coupling (low-temperature) polymer / Pirogov–Sinai:** 4D, volume-uniform, but **finite groups only** (Adhikari–Cao 2025). The method needs a unique ordered ground state and a positive defect cost `Δ_G`. **[T]**
- **Lace expansion:** **no gauge application at weak coupling.** Its successes are scalar φ⁴, self-avoiding walks, percolation; 4D gauge at weak coupling would need RG (Balaban's programme), which is incomplete. **[R]/[O]**
- **Area law at large N:** covered in `SURVEY-RECENT-ARXIV` (Cao–Nissim–Sheffield; not re-reported). It is confinement, not a gap, and not the crossover.
- **No 4D weak-coupling gap for a continuous gauge group by any cluster/polymer/lace method exists.**

### (e) Any rigorous G2 for 4D pure gauge at any coupling / volume?

- **Yes, at strong coupling, all volumes, continuous SU(N):** Osterwalder–Seiler 1978; SZZ 2023 (β<1/48, d=4); Nissim 2025 (U(N), 't Hooft). **[T]** These are lattice statements at fixed a; none takes a→0.
- **Yes, at weak coupling, 4D, but only finite groups (and a finite lattice):** Adhikari–Cao 2025. **[T]**
- **No for continuous SU(N) at weak coupling, any volume.** No for the continuum, any coupling.
- **Fixed finite volume, any β:** the transfer matrix is a finite matrix, so positivity is trivial; the content is *uniformity* in volume and a, which is exactly what is missing. **[O]**

---

## 2. Master table

| # | Method / line | Exact citation | What is actually proved | Regime (d, group, coupling, volume) | G1 / G2 | Crosses crossover? | Label |
|---|---|---|---|---|---|---|---|
| 1 | Strong-coupling cluster / transfer matrix | Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440–471, DOI 10.1016/0003-4916(78)90039-8 | Reflection positivity, positive transfer matrix, infinite-volume limit (existence + analyticity), Wilson area law; exponential clustering = lattice mass gap | 4D, SU(N)/SO(N)/U(N), **β small (strong coupling)**, all volumes | G1 lattice (fixed a) + G2 lattice | **No** — strong coupling | [T][V] |
| 2 | Log-Sobolev / Bakry–Émery | Shen–Zhu–Zhu, *Comm. Math. Phys.* **400** (2023) 805–851, DOI 10.1007/s00220-022-04609-1, arXiv:2204.12737 | Unique infinite-volume measure, LSI + Poincaré uniform in volume, exponential decay ⇒ positive mass gap | any d>1, SU(N)/SO(N), **SU(N): \|β\|<1/(16(d−1))**; d=4 ⇒ **β<1/48**, all volumes | G1 lattice + G2 lattice | **No** — strong coupling | [T][V] |
| 3 | U(N) random-environment + cluster + Langevin | Nissim, arXiv:2510.22788 (2025), preprint | Mass gap, unique infinite-volume limit, new large-N limit proof | all d≥2, U(N), **'t Hooft / small β** (intermediate steps N>8π), all volumes | G1 lattice + G2 lattice | **No** — 't Hooft / strong coupling | [T][V] (preprint) |
| 4 | YM–Higgs Langevin | Shen–Zhu–Zhu, arXiv:2401.13299 (2024; v3 2025) | Exponential ergodicity, positive mass gap | small (β,κ), compact/unbounded Higgs, any d | G1 lattice + G2 lattice | **No** — not pure YM, strong coupling | [T][V] |
| 5 | Weak-coupling Pirogov–Sinai / swapping | Adhikari–Cao, *Ann. Probab.* **53** (2025) 140–174, DOI 10.1214/24-AOP1702, arXiv:2202.10375 | Exponential decay of correlations for all Wilson-loop functions; `Cov ≤ 4⁴·10²⁴\|G\|^{2(\|B₁\|+\|B₂\|)}‖f₁‖‖f₂‖ e^{−(β/2)Δ_G(L−1)}` (constants as in Thm 1.1) | **4D, finite group G, β ≥ (11/4+4log\|G\|)/Δ_G (weak coupling), finite lattice; bound volume-independent** | G2 only (clustering) | **No** — finite group only; `Δ_G→0` for SU(N) | [T][V] |
| 6 | Two-scale disintegration | Piazza–Shen, arXiv:2610.35286 (2026), preprint | LSI at arbitrary fixed β, constant uniform in volume | **2D, SU(N)**, all couplings | G2-type (LSI ⇒ Poincaré) in 2D | **No** — 2D only | [T][V] (preprint) |
| 7 | 2D stochastic quantisation | Chevyrev–Shen, *Comm. Pure Appl. Math.* (2026), DOI 10.1002/cpa.70043 | YM measure invariant for renormalised dynamic; universality | 2D, compact G, finite volume | 2D G1 (invariance/universality) | **No** — 2D | [T][V] (status only) |
| 8 | 3D stochastic quantisation | Chandra–Chevyrev–Hairer–Shen, *Invent. Math.* **237** (2024) 541–696 | Local-in-time solutions; unique gauge-covariant counterterms; Markov process **up to blow-up** | 3D, finite volume | partial G1; no G2 | **No** — 3D, local | [T][V] (status only) |
| 9 | 3D abelian weak-coupling scaling limit | Chevyrev–Qu–Shen, arXiv:2608.27828 (2026), preprint | Rescaled U(1) phase → stochastic heat equation (free Gaussian), universal | 3D, U(1), weak-coupling scaling β=β̄ε⁻¹ | G1 (free limit); **no G2** | **No** — 3D abelian, and the limit is massless | [T]/[R] (preprint) |
| 10 | 4D abelian weak coupling (no-go) | Guth, *Phys. Rev. D* **21** (1980) 2291, DOI 10.1103/PhysRevD.21.2291; Fröhlich–Spencer, *Comm. Math. Phys.* **83** (1982) 411–454, DOI 10.1007/BF01213610 | Existence of a nonconfining / massless Coulomb phase at weak coupling | 4D, U(1), large β | **No G2** | **NG** for any general "large β ⇒ gap" theorem | [T][NG][V] |
| 11 | 3D U(1) all-coupling confinement (contrast) | Göpfert–Mack, *Comm. Math. Phys.* **82** (1981) 545, DOI 10.1007/BF01961240 | Confinement of static quarks for all couplings | 3D, U(1), all β | confinement, not 4D G2 | **No** — 3D abelian | [T][V] |
| 12 | Large-N area law | Cao–Nissim–Sheffield, arXiv:2505.16585 / 2509.04688 (2025) | Area law in expanded / 't Hooft regimes | 4D, U(N)/SU(N)/SO(2N), large N | confinement | **No** — large N; not a gap | [T] (see `SURVEY-RECENT-ARXIV`) |
| 13 | Fixed-cutoff / conditional criteria | Mihajlovic, arXiv:2611.11534 (2026); Gauvin, arXiv:2503.15539 (2025/26) | Spectral bounds / gap certificates at fixed cutoff; explicit criteria for a continuum limit; **no continuum gap claimed** | 4D regulated Hamiltonians, SU(2)/SU(3) | conditional G2 only | **No** — fixed cutoff | [C]/[R] (see `SURVEY-RECENT-ARXIV`) |

---

## 3. What is actually proved, line by line

### 3.1 Strong coupling, continuous groups, 4D, all volumes

**Osterwalder–Seiler 1978** is the anchor. It establishes the Euclidean lattice theory at strong coupling: reflection positivity (hence a positive self-adjoint transfer matrix), existence and analyticity of the infinite-volume limit, and Wilson's confinement bound. The same cluster expansion yields exponential clustering, i.e. a mass gap in the lattice transfer-matrix sense, at fixed lattice spacing. **[T]** It does **not** take a→0 and does **not** reach weak coupling.

**SZZ 2023** is the modern functional-inequality version. Their exact SU(N) condition is `|β| < 1/(16(d−1))`; in d=4 this is **β<1/48**. They prove the infinite-volume measure exists and is unique, the LSI and Poincaré inequality hold uniformly in volume, and correlations decay exponentially — "a strictly positive mass gap". **[T]** The result is independent of N, so it covers N=3. It is a genuine lattice G2 theorem at fixed a and strong coupling.

**Nissim 2025** covers U(N) in the 't Hooft regime, all d≥2, with a mass gap and a unique infinite-volume limit; the proof recasts U(N) as a random-environment SU(N) model and combines cluster and Langevin techniques. The intermediate estimates use N>8π; read it as a large-N / 't Hooft strong-coupling result. **[T]**

All three are **lattice, fixed a, strong coupling**. None is a continuum statement. None crosses β≈2.5.

### 3.2 Weak coupling, finite groups, 4D, volume-uniform — Adhikari–Cao 2025

This is the **closest approach in coupling** in the entire literature, and the most important new item for this survey.

The setting is a 4D finite lattice `Λ = [a₁,b₁]×…×[a₄,b₄] ⊂ ℤ⁴`, a **finite** gauge group `G` (possibly non-Abelian), unitary representation χ, and the Wilson measure at inverse coupling β, in the **large-β (weak-coupling)** regime. Their Theorem 1.1 states: for `β ≥ (11/4 + 4 log|G|)/Δ_G`, for any two rectangles `B₁,B₂` at distance ≥L, and any conjugacy-invariant functions `f₁,f₂` of loops in `B₁,B₂` (including arbitrary functions of Wilson loops),

```
|Cov(f₁,f₂)| ≤ 4·4·10²⁴ |G|^{2(|B₁|+|B₂|)} ‖f₁‖ ‖f₂‖ e^{−(β/2)Δ_G(L−1)}.
```

So: **exponential clustering at weak coupling, in 4D, with a bound independent of the lattice volume**, for finite non-Abelian groups. They explicitly interpret this as a mass gap in their setting. The theorem is stated on a finite lattice; the volume-independence of the bound makes an infinite-volume limit plausible but no such limit theorem is stated. **[T][V]**

**Why it does not reach SU(N) — the exact barrier.** There are four independent obstructions, each fatal:

1. **G is finite.** The authors say it plainly: *"we will make a mathematical simplification and take G to be a finite group. Of course, this lessens the direct relevance to physics."* SU(N) is never treated.
2. **The threshold diverges with |G|.** `β ≥ (11/4 + 4 log|G|)/Δ_G` grows like `log|G|`; any sequence of finite groups approximating SU(N) has `|G|→∞`, so the theorem gives no uniform statement in the approximation.
3. **The rate vanishes for continuous groups.** `Δ_G = min_{g≠1} Re[χ(1)−χ(g)]`. For a finite group this is positive. For any continuous group, `g→1` gives `Δ_G = 0`, so the decay rate `(β/2)Δ_G` is zero: the method's mass gap closes exactly in the SU(N) limit. This is the quantitative form of "the gluon is massless at tree level."
4. **The weak-coupling phase is frozen, not a continuum QFT.** The correlation length is `ξ = 2/(βΔ_G) → 0` as β→∞. The finite-group large-β model is an ordered/frozen phase with a gap that *grows* with β — the opposite of an asymptotically free continuum theory.

The authors themselves frame the open problem: *"The problem of showing exponential decay of correlations for a 4D lattice gauge theory (with gauge group a non-Abelian Lie group, say SU(3)) is one of the central questions of lattice gauge theory, as it is intimately related to showing the existence of a mass gap."* **[V]** This is the cleanest published statement that the continuous-group weak-coupling case is open.

**Acid test.** The theorem is not vacuous: the finite-group model really does have a gap. But it is the wrong object. The free massless field is not in the hypothesis class (finite group), and the mechanism is a low-temperature/Pirogov–Sinai expansion around a unique frozen ground state. It cannot manufacture positivity for SU(N), because for SU(N) the "defect cost" `Δ_G` is zero.

### 3.3 All couplings, 2D only — Piazza–Shen 2026

Piazza–Shen prove the LSI for the 2D Wilson SU(N) measure at **arbitrary fixed β**, with a constant depending on β and N but **uniform in the lattice volume**, via a two-scale disintegration argument. **[T][V]** This is the only all-coupling LSI for SU(N) gauge found in this sweep.

Why it does not lift to 4D:
- 2D YM is exactly solvable (Migdal/Makeenko); the all-coupling statement exploits this solvability and a two-scale disintegration.
- The coarse-scale step requires **uniformly positive Ricci curvature** of the renormalised block measure. In 2D this is available. In 4D the coarse Hessian is the Savvidy-indefinite object; there is no positive curvature input.
- The related 2D U(1)/U(N) extensions (reported in the search) also need positive Ricci curvature; groups with a nontrivial Lie-algebra centre obstruct it.

So 2D "all couplings" is a consequence of solvability plus a favourable curvature step, not a template that can be transported to 4D.

### 3.4 Stochastic quantisation and regularity structures

See §1(a)–(b). The status is:
- 2D: global, invariant measure.
- 3D: local-in-time, Markov process up to blow-up; no invariant measure/global solution; renormalisation uniqueness.
- 4D: nothing; d=4 is renormalisable/marginal, so the subcriticality that the machinery needs is absent.
- The one weak-coupling continuum gauge limit is 3D **abelian**, and it is the free stochastic heat equation (gapless).

**This is the strongest structural statement in the SPDE line:** the method's only weak-coupling success produces the free field. That is exactly the object the project's acid test says a valid gap tool must fail on.

### 3.5 Cluster / polymer / lace

The cluster/polymer line splits cleanly by coupling:
- **Small β (strong coupling):** works for continuous groups, gives a gap (Osterwalder–Seiler; Seiler; SZZ).
- **Large β (weak coupling):** works only for finite groups (Adhikari–Cao), because the low-temperature expansion needs a unique ordered ground state and `Δ_G>0`.

There is **no lace-expansion result for gauge theory at weak coupling**. Lace expansion is a technique for scalar/branching systems; 4D gauge at weak coupling needs a multiscale RG (Balaban), and that programme is incomplete. **[R]/[O]**

### 3.6 Large N and area law

Cao–Nissim–Sheffield and related work prove area laws in expanded / 't Hooft regimes. These are **confinement** statements, not gap statements, and they are large-N / strong-coupling. They are covered in `SURVEY-RECENT-ARXIV` and are not re-reported here.

### 3.7 Fixed-cutoff criteria / conditional results

Mihajlovic (arXiv:2611.11534 (2026)) and Gauvin (arXiv:2503.15539 (2025/26)) prove spectral bounds / gap certificates for **regulated** Hamiltonians and isolate explicit criteria for a continuum limit. Both explicitly state that the interacting 4D construction and a uniform physical mass gap remain unproved. **[C]/[R]** See `SURVEY-RECENT-ARXIV` for details.

### 3.8 The abelian no-go

**Guth 1980** and **Fröhlich–Spencer 1982** establish a nonconfining, massless Coulomb phase in 4D U(1) lattice gauge theory at weak coupling. **[T][V]** (Contrast: 3D U(1) is confined for all couplings — Göpfert–Mack 1981. **[T][V]**)

This is a decisive structural fact for the constructive programme: **there is no theorem of the form "compact group + large β ⇒ mass gap."** Any such theorem is false, with U(1)₄ as the counterexample. The nonabelian gap at weak coupling is therefore a genuinely nonperturbative/confinement phenomenon. This is the constructive-literature version of the project's thesis "transport ≡ confinement": crossing the crossover is not a matter of a better expansion at large β; it is the same difficulty as proving confinement.

---

## 4. Deepest analysis: the closest approaches and their exact barriers

### 4.1 Adhikari–Cao 2025 — closest in coupling

**What it achieves.** A genuine 4D, weak-coupling, volume-uniform exponential-clustering theorem for Wilson loops, for finite (including non-Abelian) groups. Nothing else in the literature reaches large β in 4D with a gap, even for a toy group. It is the strongest evidence that the *low-temperature/Pirogov–Sinai* route can be made to work in 4D when the ground state is isolated.

**Exact barrier.** The rate is `(β/2)Δ_G` and the threshold is `(11/4+4log|G|)/Δ_G`. For SU(N), `Δ_G=0`: the identity is not isolated, there are arbitrarily small group elements, and the "defect cost" vanishes. Equivalently, the perturbative SU(N) vacuum is critical, with massless gluons; there is no unique gapped ground state around which to expand. The theorem's entire mechanism disappears in the SU(N) limit. This is not a technical gap to be closed by a better estimate; it is the wrong phase.

**Acid test.** The finite-group model has a gap that grows with β; the free massless field does not. The theorem is therefore not vacuous, but it certifies a frozen phase, not a continuum gauge theory. It passes "has a gap" and fails "is the right object."

### 4.2 SZZ 2023 / Nissim 2025 — closest in physicality

**What they achieve.** For actual SU(N) (SZZ) or U(N) (Nissim), in 4D, infinite volume, they prove a unique measure, an LSI/Poincaré inequality uniform in volume, and exponential clustering = a positive lattice mass gap. This is the closest the literature comes to a *physical* G2 theorem.

**Exact barrier.** The SZZ proof is Bakry–Émery: the LSI constant is controlled by a curvature lower bound `Γ₂ ≥ K > 0`. On the lattice gauge configuration space, the curvature is a sum of a group-geometric term and the Hessian of the Wilson action. At small β the geometric term dominates and K>0; at β≈2.5 the action Hessian `−D² + 2 ad F_B` is indefinite (the Savvidy negative mode), so `Γ₂` can be negative and the method gives nothing. The project has measured this negative mode (7.1% negative after the gauge quotient). The "Dobrushin ceiling" is the single-link version of the same failure: the Dobrushin coefficient exceeds 1 beyond a threshold. **[O]/[NG for the method class]**

This is a **method** barrier, not a proof that the LSI constant actually becomes negative. The true LSI constant could remain positive at β≈2.5; the point is that no known technique proves it. The project's H1 is precisely the statement that the SZZ-type variance upper bound survives to β≈2.5; the constructive literature proves it at β<1/48 and nowhere else in 4D.

### 4.3 Piazza–Shen 2026 — the only all-coupling LSI

**What it achieves.** LSI for 2D SU(N) at every β, volume-uniform. This is a genuine all-coupling result and shows that "all couplings" is not logically impossible.

**Exact barrier.** 2D solvability plus a two-scale disintegration whose coarse step needs uniformly positive Ricci curvature. In 4D the coarse Hessian is Savvidy-indefinite; the curvature input is absent. The 2D result is a template for *how* an all-coupling LSI could look, not a method that transports.

### 4.4 Why the barrier is structural: three routes, one wall

The constructive literature contains exactly three routes to a 4D G2 theorem, and each dies for a precise reason at the same physical point:

| Route | Natural radius | Why it stops | What would be needed |
|---|---|---|---|
| Strong-coupling / high-temperature expansion (OS, Seiler, SZZ, Nissim) | β < 1/48 (d=4) | Finite radius of convergence (Pringsheim); no resummation known; Bakry–Émery curvature positive only at small β | A resummation / analytic continuation past the radius that preserves a sign-controlling lower bound |
| Weak-coupling / low-temperature Pirogov–Sinai (Adhikari–Cao) | finite G, β large | Needs an isolated ordered ground state (`Δ_G>0`); `Δ_G=0` for SU(N) | A low-temperature expansion for a *critical* vacuum with massless gluons — i.e. a proof of confinement/nonperturbative mass generation |
| SPDE / regularity structures (CCHS, Chevyrev–Shen, Chevyrev–Qu–Shen) | d=2,3 | d=4 is marginal; subcriticality fails; 3D global solution/invariant measure still open | A critical SPDE construction at the level of 4D φ⁴, which does not exist |

All three converge on the same missing object: a **sign-controlling, nonperturbative, `(L,a)`-uniform lower bound at β≈2.5** that fails on the free massless field. That is the project's H1 / `f₁²>0` transport. The constructive literature has not produced it, and its three natural methods each fail at their own boundary before reaching it.

### 4.5 Acid test

For each candidate:

- **OS / SZZ / Nissim:** genuine lattice gaps, strong coupling. They do not apply to the free massless field (the model is a compact-group lattice gauge theory at small β, far from Gaussian), and they give no false positivity. They are not vacuous, but they are regime-limited.
- **Adhikari–Cao:** finite groups; the free field is not in the class; the gap grows with β (frozen phase). Not vacuous, but the wrong object.
- **Piazza–Shen:** 2D; the 2D theory is solvable and gapped. Not vacuous, but 2D.
- **Chevyrev–Qu–Shen:** the weak-coupling limit *is* the free massless field. This is the acid test failing in the clearest possible way: the only weak-coupling continuum gauge construction returns zero gap.
- **Bakry–Émery / Dobrushin:** at β≈2.5 the curvature input is negative; the method returns nothing precisely on the gapped side, while it is correct (gives no false positivity) on the free side. It cannot certify the gap.

**Conclusion of the acid test:** no method in the literature is a candidate crossing tool. The successful ones are regime-limited; the weak-coupling one is free; the curvature one fails at the physical point.

---

## 5. The precise missing step (honest)

A method that could cross the crossover would need all of:

1. **Continuous SU(N)** (not a finite group, not U(1), not 2D).
2. **β≈2.5** (α_s~O(1)), with **no small parameter**.
3. **A sign-controlling lower bound** on the `0⁺⁺` spectral edge `m₁` and/or the residue `f₁²`, not an upper bound and not a value.
4. **`(L,a)`-uniform** (infinite volume and continuum limit).
5. **Nonperturbative and non-circular** (does not assume the gap/condensate it is trying to prove).
6. **Fails on the free massless field** (`f₁²=0`, `G″=0`).

No method in the constructive literature satisfies 1–6. The closest structural ideas each satisfy a subset:
- SZZ/Nissim: 1, 3, 4, 5, 6 — but only at β<1/48 (fail 2).
- Adhikari–Cao: 2, 3, 4, 6 — but only finite groups (fail 1), and `Δ_G=0` for SU(N).
- Piazza–Shen: 1, 2, 3, 4, 5 — but only 2D (fail 4 in 4D).
- Chevyrev–Qu–Shen: 1 (abelian), 2, 4 — but the limit is free (fail 3, 6).
- Bakry–Émery / Dobrushin: 1, 5, 6 — but fails 2, 3 at β≈2.5.

This is a precise, citable statement of the missing step. It is the constructive-literature shadow of the project's H1.

---

## 6. Claims not counted (unrefereed / already audited)

These are **not** counted as [T] and do not change the top-line. They are listed for completeness; see `SURVEY-RECENT-ARXIV` and `SURVEY-AI-FORMAL` for fuller audits.

- **Faizal–Shabir**, *Reflection-Positive Construction of a Four-Dimensional SU(N) Yang–Mills Theory with Mass Gap and Confinement*, arXiv:2606.19362, *Fortschr. Phys.* **74** (2026) e70097. Already audited in `BRIEF-FOR-MODEL §5` (T3): the local-to-global step, the weak-coupling identification of `g_bare`, the common-action hypothesis, and the harmonic divergence at Eq. 10.72 do not close. Not counted.
- **Bhattacharjee–Bhattacharya–Lockwood** Zenodo/Lean series ("Dobrushin ceiling", "one cube estimate", "complete construction"). Non-refereed; the claimed thresholds are still strong-coupling (`β_W ≲ 0.35–0.76`), and the papers themselves state the continuum/crossover remains open. **[R]** Not counted.
- **Sakib**, Cambridge Engage working paper DOI 10.33774/coe-2026-4q17x (2026). Non-refereed; claims `β<195/7168` (×1.31 over 1/48) in 4D from a gauge-orbit curvature improvement; the author concedes production couplings are ~24× larger and the continuum is out of reach. **[R]** Not counted.
- **Eriksson** viXra series 2602.0051–2602.0054 (2026). Conditional/windowed weak-coupling LSI claims; later versions withdraw the unconditional claims. **[R]** Not counted.
- **Mihajlovic** arXiv:2611.11534 (2026), **Gauvin** arXiv:2503.15539 (2025/26). Fixed-cutoff / conditional criteria; both explicitly say no continuum gap is proved. **[C]/[R]** Not counted as G2.
- **Rumoured "OpenAI / Sussman–Alpöge" Balaban-route constructions** (Oct–Nov 2026, non-arXiv, no public manuscript at the time of writing). Unverifiable from the literature; not counted. (The OpenAI Lean G1-only result is discussed in `SURVEY-AI-FORMAL`; it explicitly excludes the mass gap.)

---

## 7. Verified citation list

**Primary theorems [T]**

1. K. Osterwalder, E. Seiler, "Gauge field theories on a lattice," *Ann. Phys.* **110** (1978) 440–471. DOI 10.1016/0003-4916(78)90039-8. **[V]**
2. E. Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics*, Lecture Notes in Physics **159**, Springer (1982). **[V]**
3. H. Shen, R. Zhu, X. Zhu, "A stochastic analysis approach to lattice Yang–Mills at strong coupling," *Comm. Math. Phys.* **400** (2023) 805–851. DOI 10.1007/s00220-022-04609-1; arXiv:2204.12737. **[V]**
4. R. Nissim, "U(N) lattice Yang–Mills in the 't Hooft regime," arXiv:2510.22788 (2025), preprint. **[V]**
5. H. Shen, R. Zhu, X. Zhu, "Langevin dynamics of lattice Yang–Mills–Higgs and applications," arXiv:2401.13299 (2024; v3 2025); reported *Comm. Math. Phys.* (2026). **[V]** (arXiv)
6. A. Adhikari, S. Cao, "Correlation decay for finite lattice gauge theories at weak coupling," *Ann. Probab.* **53** (2025) 140–174. DOI 10.1214/24-AOP1702; arXiv:2202.10375. **[V]**
7. J. Piazza, H. Shen, "Log-Sobolev inequality for 2D SU(N) lattice gauge theory at all couplings," arXiv:2610.35286 (2026), preprint. **[V]**
8. A. Guth, "Existence proof of a nonconfining phase in four-dimensional U(1) lattice gauge theory," *Phys. Rev. D* **21** (1980) 2291. DOI 10.1103/PhysRevD.21.2291. **[V]**
9. J. Fröhlich, T. Spencer, "Massless phases and symmetry restoration in abelian gauge theories and spin systems," *Comm. Math. Phys.* **83** (1982) 411–454. DOI 10.1007/BF01213610. **[V]**
10. M. Göpfert, G. Mack, "Proof of confinement of static quarks in three-dimensional U(1) lattice gauge theory for all values of the coupling constant," *Comm. Math. Phys.* **82** (1981) 545. DOI 10.1007/BF01961240. **[V]**

**Status / structural [T]/[R]**

11. I. Chevyrev, "Stochastic quantisation of Yang–Mills," *J. Math. Phys.* **63** (2022). DOI 10.1063/5.0089431; arXiv:2202.13359. **[V]** (quotes in §0–§1).
12. A. Chandra, I. Chevyrev, M. Hairer, H. Shen, "Stochastic quantisation of Yang–Mills–Higgs in 3D," *Invent. Math.* **237** (2024) 541–696. DOI 10.1007/s00222-024-01264-2; arXiv:2201.03487. **[V]** (status only; excluded from re-reporting).
13. I. Chevyrev, H. Shen, "Invariant measure and universality of the 2D Yang–Mills Langevin dynamic," *Comm. Pure Appl. Math.* (2026). DOI 10.1002/cpa.70043; arXiv:2302.12160. **[V]**
14. I. Chevyrev, Y. Qu, H. Shen, "Scaling limit of the 3D abelian Yang–Mills Langevin dynamics," arXiv:2608.27828 (2026), preprint. **[V]**
15. I. Chevyrev, H. Shen, "Uniqueness of gauge covariant renormalisation of stochastic 3D Yang–Mills–Higgs," *Arch. Ration. Mech. Anal.* **250** (2026) Paper No. 11. DOI 10.1007/s00205-025-02163-3; arXiv:2503.03060. **[V]**
16. A. Sevostyanov, "A tree-free approach to 3D Yang–Mills Langevin dynamic," arXiv:2605.14616 (2026). **[V]** (preprint).
17. S. Cao, R. Nissim, S. Sheffield, "Expanded regimes of area law for lattice Yang–Mills theories," arXiv:2505.16585 (2025); "Dynamical approach to area law for lattice Yang–Mills," arXiv:2509.04688 (2025). **[V]** (pointer only; see `SURVEY-RECENT-ARXIV`).
18. S. Chatterjee, "A probabilistic mechanism for quark confinement," *Comm. Math. Phys.* **385** (2021) 1007–1039; arXiv:2006.16229. **[R]** (strong-coupling/confinement, not a weak-coupling gap).
19. S. Cao, S. Chatterjee, "A state space for 3D Euclidean Yang–Mills theories," *Comm. Math. Phys.* **405** (2024) art. 3. **[R]** (3D).
20. A. Mihajlovic, "SU(2) Yang–Mills: Spectral bounds, vacuum concentration, and criteria for the continuum limit," arXiv:2611.11534 (2026). **[V]** (preprint; conditional/fixed-cutoff).
21. S. Gauvin, "Pinched multi affine geometry and confinement: describing the Yang–Mills mass gap," arXiv:2503.15539 (2025; v3 2026). **[V]** (preprint; conditional/fixed-cutoff).

**Unrefereed / not counted [R]**

22. M. Faizal, A. Shabir, arXiv:2606.19362, *Fortschr. Phys.* **74** (2026) e70097. (Already audited in `BRIEF-FOR-MODEL §5`.)
23. D. Bhattacharjee, U. Bhattacharya, et al., Zenodo records 23484043 / 23788087 (2026).
24. S. M. N. Sakib, Cambridge Engage working paper, DOI 10.33774/coe-2026-4q17x (2026).
25. L. Eriksson, viXra 2602.0051–2602.0054 (2026).
26. M. Douglas, *Nature Rev. Phys.* **8** (2026) 86–97, DOI 10.1038/s42254-025-00909-2. (Review; "no clear path known".)

---

## 8. Verdict

**No constructive or rigorous-QFT method crosses the 4D Yang–Mills crossover.** The literature contains genuine 4D G2 theorems at strong coupling for continuous SU(N), all volumes (Osterwalder–Seiler 1978; SZZ 2023; Nissim 2025), one genuine 4D weak-coupling G2 theorem for finite groups (Adhikari–Cao 2025), and one all-coupling LSI for SU(N) gauge in 2D (Piazza–Shen 2026). It contains no 4D weak-coupling G2 theorem for a continuous gauge group, no 4D nonabelian continuum measure, and no method that produces a sign-controlling `(L,a)`-uniform lower bound at β≈2.5.

The three natural method classes fail for three precise, independent reasons at the same physical point: the strong-coupling radius (1/48), the vanishing defect cost `Δ_G=0` for continuous groups, and the failure of subcriticality in 4D. The abelian no-go (4D U(1) is massless at weak coupling) shows that "large β ⇒ gap" is false as a general principle. The only weak-coupling continuum gauge construction in the SPDE line is 3D abelian and its limit is the free massless field.

The missing object is exactly the project's H1: a nonperturbative, non-circular, `(L,a)`-uniform, sign-controlling lower bound on the `0⁺⁺` edge/residue at β≈2.5 that fails on the free massless field. No such method exists in the constructive literature. **Clay remains open; no false hope is warranted from this sweep.**
