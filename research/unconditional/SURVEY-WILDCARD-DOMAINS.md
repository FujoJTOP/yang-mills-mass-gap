# SURVEY — Wildcard Domains Sweep: Is There a Genuinely New Handle on the 4D Gap?

> 2026-11-24 · Scope: algebraic topology / K-theory / homotopy / equivariant; representation theory / tensor models / spin foam / recoupling; noncommutative geometry (Connes, spectral triples); arithmetic / number theory (incl. the alleged Douglas program, zeta, p-adic); category theory / higher algebra.
> Question for every candidate: exact citation, actual mechanism, and **does it give a new handle on (`m₁²`, `f₁²>0` isolated) at the crossover `α_s~O(1)`, or does it reduce to a known wall** (atom-blind positivity / rigidity / wrong channel / crossover-no-small-parameter)? Acid test every proposed gap tool on the **free massless field** (`f₁²=0`, `G″=0`): a legitimate tool must **fail** there.
> Labels (project charter §2): **[T]** true theorem · **[C]** conditional theorem · **[O]** exact missing step · **[NG]** no-go for a named class · **[R]** reported/verified external result · **[J]** judgment. **[V]** = citation verified this session.

---

## 0. TOP-LINE

**No wildcard domain yields a genuinely new, non-circular, acid-test-passing lower bound on the `0⁺⁺` edge atom `f₁²` across the crossover.** [J] Every candidate either (i) is **rigid** (discrete deformation-invariant ⇒ frozen at the free value 0, registry B3), (ii) is **atom-blind positivity** (Schur-concave ⇒ only upper-bounds `m₁²`, registry §3), (iii) reads the **wrong channel / wrong sector** (instantons→`0⁻⁺`, spin-foam/tensor→gravity or large-N), or (iv) has the **correct logical shape** but is **circular/equivalent-hard**.

**One candidate is genuinely serious and worth the deepest analysis** (below, §3): the **arithmetic / trace-formula cluster** — Weil positivity & the explicit formula (Weil 1952, Burnol 1998), Connes' noncommutative trace formula with zeros-as-absorption-spectrum (1999), and Selberg/Gutzwiller. Reason it is serious: it is the **only** family in the list with the *correct logical direction* — a **positivity ⟺ exact spectral location** (zeros on a line / atom at `m₁`) rather than "positivity ⟹ upper bound." It is atom-seeing and free-field-failing (passes the acid test in principle). **Verdict [O]:** it reduces **exactly** to the project's own **A3 (Fredholm determinant identity)** plus **B6 (reflection positivity = `ν≥0`)**, and its residual — locating the determinant zero in the continuum — is the marginal BPHZ / H1 wall. The arithmetic dressing adds **no new material**; it *independently re-derives A3 as the unique atom-seeing structure*. [J]

**Important correction to the task premise [V]:** there is **no "M. Douglas arithmetic Yang–Mills" program.** Douglas's 2026 review (*Nature Rev. Phys.* **8**, 86–97) advocates **rigorous stochastic quantization + strong-coupling expansion** and states plainly *"no clear path is known"* to the mass gap; it contains **no arithmetic / number-theory route**. "Arithmetic gauge theory" (Minhyong Kim) is real mathematics but is arithmetic geometry (primes as knots, arithmetic Chern–Simons/BF), **not** a 4D mass-gap handle. The Douglas-arithmetic lead is a **false premise**; the stochastic-quantization lead is already the project's G1/L3 line. [J][V]

---

## 1. TABLE — candidates, mechanism, reduction, verdict

### A. Algebraic topology / equivariant / K-theory / homotopy

| domain | candidate | citation | mechanism | reduces to? | verdict |
|---|---|---|---|---|---|
| A-top | Atiyah–Singer index / APS **eta invariant** of the Dirac coupled to `A` | Atiyah–Patodi–Singer, *Math. Proc. Camb. Phil. Soc.* **77** (1975) 43; standard [R] | integer (index) or signed spectral asymmetry `η`; detects instanton number / spectral flow | **Rigidity (B3):** integer & locally constant along the connected β-family ⇒ equals free value; `η` only on an odd-dim (5D) boundary, wrong channel | [NG] |
| A-top | **Analytic (Ray–Singer) torsion** | Ray–Singer 1971; **Cheeger 1979 / Müller 1978** (Cheeger–Müller thm), Müller *JAMS* **6** (1993) 721 [V] | zeta-regularized `det` of Hodge Laplacian; a *secondary* (non-integer) spectral invariant | **normalized torsion = Reidemeister torsion (topological/rigid) ⇒ B3**; raw torsion = **determinant = A3** | [NG] |
| A-top | **Equivariant localization / Nekrasov partition function** on instanton moduli | Nekrasov, *Adv. Theor. Math. Phys.* **7** (2003) 831 [R] | equivariant index with `ε₁,ε₂`; exact in SUSY | requires **SUSY**; rigid; non-SUSY YM has no protected localization ⇒ condensation | [NG] |
| A-top | **Homotopy of `𝒜/𝒢`** (`π₀=ℤ` instanton, `π₁`Gribov) ⇒ θ-vacua; **Gribov horizon** | θ-vacua standard (e.g. KEK 2000 preprint); Gribov 1978; Zwanziger horizon cond. [R] | instanton sectors → topological susceptibility → Witten–Veneziano `η′` mass; horizon condition fixes a scale by a *gap equation* | **wrong channel `0⁻⁺`** (registry D-dimension: instanton native `0⁻⁺`); Gribov scale **circular** (assumes scenario) ⇒ H1 | [NG] |
| A-top | **Twisted equivariant K-theory / Freed–Hopkins–Teleman Verlinde** | FHT, *Ann. Math.* **172** (2010) 1875 [R] | `K_G^τ(G)` = Verlinde algebra for 2D/3D (Chern–Simons/WZW) | **2D/3D only** (loop-group level); no 4D transfer-operator statement; K-theory discrete ⇒ rigid | [NG] |
| A-top | **Cobordism / generalized-symmetry anomaly** (center 1-form `ℤ_N`, `Ω_5^Spin(Bℤ_N)`) | Gaiotto–Kapustin–Komargodski–Seiberg, *JHEP* **05** (2017) 091, arXiv:1703.00501 [V] | discrete 't Hooft anomaly: θ=π vacuum cannot be trivially gapped; constrains confinement / center breaking | **rigid (discrete)** + **atom-blind to `0⁺⁺`** (the glueball is *neutral* under center 1-form) ⇒ sees confinement, not the atom | [NG] |

### B. Representation theory / tensor models / spin foam / recoupling

| domain | candidate | citation | mechanism | reduces to? | verdict |
|---|---|---|---|---|---|
| B-rep | **Character / recoupling (6j,9j) expansion** | lattice character expansion; Creutz; recoupling for SU(N) [R] | sum over irreps → the **P3 strong-coupling residue theorem** `f₁²≥c₁u^p>0` at `β<1/48` | **already the engine of the one unconditional result (P3)**; series radius `β<1/48`; **crossover = the wall** (no analytic continuation) | [T] strong-coupling / [O] crossover |
| B-rep | **Spin foam** (Barrett–Crane 4D, Ponzano–Regge 3D) | Barrett–Crane, *J. Math. Phys.* **39** (1998) 3296 [R] | recoupling amplitudes = BF/gravity state sums | **wrong sector**: BF/gravity, not YM; no dynamical mass scale | [NG] |
| B-rep | **Random tensor models** (Gurau melonic; SYK-like) | Gurau, *CMP* **304** (2011) 69; Klebanov–Tarnopolsky; Witten (SYK) [R][V] | `1/N` expansion, melonic dominance ⇒ conformal `0+1`-D / emergent-gravity universality class | **large-N melonic / wrong dimension & universality class** (= large-N line already dead) | [NG] |
| B-rep | **Representation-ring / Bochner heat coefficients** | (already in registry §3) | heat-kernel coefficients all strictly positive | **free-field vacuous**: free massless field has positive coefficients yet `f₁²=0` ⇒ tool satisfied by no-gap field | [NG] |

### C. Noncommutative geometry (Connes, spectral triples)

| domain | candidate | citation | mechanism | reduces to? | verdict |
|---|---|---|---|---|---|
| C-NCG | **Spectral action** `Tr f(D/Λ)` on almost-commutative triple | Chamseddine–Connes, *CMP* **186** (1997) 731; Connes–Marcolli, *NCG, Quantum Fields and Motives* (2008) [R] | seeley–deWitt heat expansion yields YM+Higgs+gravity at fixed `Λ` | **heat kernel = representation ring ⇒ free-field vacuous** (positive coefficients, `f₁²=0`); "spectral gap" of internal `D_F` is put in by hand (Higgs mass) | [NG] |
| C-NCG | **Spectral-triple reconstruction / Dirac index** | Connes, *J. Noncomm. Geom.* (2013); index ⇒ Atiyah–Singer | reconstructs geometry from `(A,H,D)`; index is topological | **rigidity (B3)**; no dynamical atom | [NG] |
| C-NCG | **Connes–Marcolli "motives" / spectral zeros** | Connes, *Selecta Math.* **5** (1999) 29–106 [V] | zeros as **absorption spectrum** on adele classes | **= A3 + arithmetic (§3)** | [O] |

### D. Arithmetic / number theory

| domain | candidate | citation | mechanism | reduces to? | verdict |
|---|---|---|---|---|---|
| D-arith | **"Douglas arithmetic YM"** | — (**does not exist**) | — | Douglas 2026 advocates stochastic quantization, "no clear path" [V] | [J] false premise |
| D-arith | **Weil positivity / explicit formula** | Weil 1952; **Burnol** arXiv:math/9809119, math/9810169 [V] | `RH ⟺ Q_W(g)≥0`; spectral side ⟷ geometric side (periodic orbits); `−log|x|` = 2D free-boson propagator | for ζ **= RH itself (circular)**; YM analog **= A3 + B6** (positivity=`ν≥0`, upper bounds) | [O] |
| D-arith | **Selberg / Gutzwiller trace formula** | Selberg 1956; Hejhal; Gutzwiller [R] | `Σ eigenvals = Σ periodic orbits`; `Tr T_O^k = Z_{S¹_k×T³}` (A3) | lowest eigenvalue needs **prime-geodesic theorem = short flux-tube count at crossover = H1/confinement** | [O]/[NG] |
| D-arith | **p-adic / adelic QFT** | Vladimirov–Volovich–Zelenov (p-adic QFT) [R] | product formula over places | **different theory**; real place carries the gap and is decoupled; no 4D statement | [NG] |
| D-arith | **Arithmetic gauge theory** (primes = knots) | Minhyong Kim, *Mod. Phys. Lett. A* **33** (2018) 1830012, arXiv:1712.07602 [V] | arithmetic Chern–Simons/BF; Galois analog of YM moduli | **analogy only** (arithmetic schemes); no dynamical mass gap | [NG] |

### E. Category theory / higher algebra

| domain | candidate | citation | mechanism | reduces to? | verdict |
|---|---|---|---|---|---|
| E-cat | **Factorization algebras** (quantum observables `E_n`) | Costello–Gwilliam, *Factorization Algebras in QFT* I (2017), II (2021), CUP [R][V] | BV quantization; mass gap = exponential clustering rate of connected correlators | **re-describes** gap as decay rate; gives **no lower-bound tool** on the rate ⇒ H1 | [O]/[NG] |
| E-cat | **Generalized / categorical symmetry** (1-form center, non-invertible) | Gaiotto–Kapustin–Komargodski–Seiberg 2017 [V]; non-invertible (Córdova–Ohmori-type) [R] | symmetry defects, anomaly, confinement order parameter | **atom-blind**: `0⁺⁺` glueball neutral under center 1-form; rigid anomaly ⇒ B3 | [NG] |
| E-cat | **Higher deformation quantization** (twisted `N=4`, `E₄`-algebras) | Costello; Kapustin–Witten (2007); Gwilliam–... [R] | `E₄`-algebra of local operators in topological twist | **SUSY/TQFT only**; wrong sector | [NG] |
| E-cat | **Koszul duality / little-disks operad** | standard (Ayala–Francis; Lurie) [R] | algebraic constraints on `E_n`-observables | constrains **algebra of observables**, not Hamiltonian spectrum lower bound | [NG] |

---

## 2. What the sweep establishes (structural, not new claims)

1. **The death map is closed under these five domains.** [J] Every row lands in registry §3 (positivity / rigidity / wrong channel) or in the A3/H1 residual. No row is class-(c) "genuinely independent open."
2. **The domains split cleanly into the two known faces** (charter §4, "两坐标墙"):
   - **Positivity face** (Schur-concave, wrong direction, atom-blind): spectral action / heat-kernel / Bochner, Weil/asymptotic positivity, factorization-algebra clustering, Gribov variational scale.
   - **Rigidity face** (discrete, right channel impossible, frozen at free value): index, eta, analytic torsion (Cheeger–Müller), Verlinde K-theory, discrete 't Hooft anomaly, Nekrasov-localization, spin-foam BPS.
   - **Wrong-channel / wrong-sector** (neither face applies cleanly): instanton→`0⁻⁺`, spin foam→gravity/BF, tensor→melonic large-N, twisted N=4→SUSY.
3. **No third face appears.** [J] The one candidate that *looked* like it might be a third face — the arithmetic **"absorption spectrum"** (an atom = a *deficit* of continuous spectral weight, a non-positivity) — is on inspection the **thimble/Schur-convex face in arithmetic clothing**: computing the deficit needs the continuum *exactly*, which is A3's zero-finding = circular.

---

## 3. DEEPEST ANALYSIS — the arithmetic / trace-formula cluster (the single serious candidate)

**Why it is the most serious.** [J] The project's entire death map rests on one asymmetry: *positivity tools are Schur-concave, so they bound `m₁²` from above but never certify the edge atom `f₁²>0`.* The arithmetic trace-formula family is the **only** item in the wildcard list that does not have this shape. It is built on an **equivalence**:

> **Weil criterion [R]:** the Riemann Hypothesis ⟺ `Q_W(g) ≥ 0` for all test functions `g` (positivity of the Weil functional, `Q_W(g)=C(g∗g^τ)`).
> **Connes [V]:** the nontrivial zeros are the **absorption spectrum** of a trace formula on the noncommutative space of adele classes; any off-line zero is a *resonance*.

The logical form is `positivity ⟺ (all zeros on a line)`, i.e. `positivity ⟺ exact spectral location` — **not** `positivity ⟹ upper bound`. If a YM analog existed, it would be exactly the missing object: a positivity certifying that the spectral edge (the `0⁺⁺` atom) sits at `m₁>0`, isolated. **Acid test passes in principle:** the free massless field has a continuous spectrum from `m=0` (no isolated edge atom, `f₁²=0`); any "zeros lie outside a disk / on a line at `m≥Δ`" criterion is **automatically violated** by the free field (its support touches `m=0`). So the criterion is non-vacuous and fails correctly on the free field. [J]

**Why it nonetheless reduces (the exact circularity).** [O]

1. **The YM "zeta function" is the dynamical (Ruelle/Fredholm) determinant `D_O(z)=det(1−zT_O)`, and the atom is its smallest zero** — this is **precisely registry A3**. [T-identity, already proved in project] The "spectral side" of the explicit formula is the mass/atom sum; the "geometric side" is the sum over closed flux tubes `Tr T_O^k = Z_{S¹_k×T³}` — again literally A3's Selberg–Gutzwiller form.

2. **The "positivity of the Weil functional" is reflection positivity `ν≥0`** — this is **registry B6**: complete monotonicity / `C(t)=∫e^{−st}dν(s)≥0` is *the same statement*. B6 proved this yields **only decreasing upper bounds** `r(t):=−(1/t)log(C(t)/M₀) ↓ m₁²`. So the geometric-side positivity, transported back to the spectral side, is exactly the positivity that the death map already classified as Schur-concave.

3. **Therefore the arithmetic criterion is logically equivalent to the gap, with no slack:** proving `Q_W ≥ 0` over the full infinite test space is, for `ζ`, **the RH itself** (known only for restricted-support test functions — Burnol's invariant analysis, `[1/c,c]`). For YM it is the statement that `D_O` has no zero inside the disk of radius `e^{aΔ}` — i.e. **evaluating the Fredholm determinant zeros in the continuum limit = A3's residual = the marginal BPHZ / H1 wall.** Burnol himself judged the naive "temporal/propagator" interpretation "misguided." [V]

**Net.** The arithmetic cluster is a **faithful re-encoding of A3 + B6** in number-theoretic language. Its *value* is diagnostic: it shows that **A3 (Fredholm determinant zeros) is the attractor** — the unique atom-seeing structure that arithmetic (Weil/Connes), analytic topology (Ray–Singer torsion denominator), and the project's own invention all independently converge to. Its *cost* is the same as A3's: the wall survives as "locate the determinant zero in the continuum." **No new material; no new lower bound; not a route.** [O][J]

**Precise missing lemma (so this is falsifiable, per charter §5):**
> **(L-W)**: For the gauge-invariant transfer operator `T_O` on `H_phys=L²(𝒜/𝒢,μ̄)`, prove that `D_O(z)` has a zero at `|z|=e^{am₁}>1` with `m₁` bounded below `(L,a)`-uniformly, **without** first assuming the gap. Equivalently: a **Weil-type positivity on the flux-tube (periodic-orbit) side** that is *provable* rather than equivalent to the spectral statement.
> Acid-test: (L-W) must fail for Maxwell (`f₁²=0`, no zero inside any `|z|>1`). It does. [J]
> **Status:** (L-W) over the full test space ⇔ the gap (circular, §3.3). No restricted-support version is known to close it. This *is* H1/A3-rest. Not solved here.

---

## 4. Citation verification log (charter §2 — checked this session, [V])

- **Douglas**, "The Yang–Mills Millennium problem," *Nature Rev. Phys.* **8**, 86–97 (2026), DOI 10.1038/s42254-025-00909-2 — **verified**: advocates rigorous stochastic quantization + strong-coupling expansion; states *"no clear path is known"* to the gap; **no arithmetic content**. (Refutes the task's Douglas-arithmetic premise.)
- **Connes**, "Trace formula in noncommutative geometry and the zeros of the Riemann zeta function," *Selecta Math. (N.S.)* **5**(1) 29–106 (1999), DOI 10.1007/s000290050042 — **verified**; zeros as absorption spectrum, trace formula on adele classes.
- **Burnol**, "The Explicit Formula and a Propagator," arXiv:math/9809119 (1998); "The Explicit Formula in simple terms," arXiv:math/9810169 (1998) — **verified**; Weil positivity as covariance `≥0`; `−log|x|` = 2D free boson propagator; author judges temporal interpretation "misguided."
- **Gaiotto–Kapustin–Komargodski–Seiberg**, "Theta, Time Reversal, and Temperature," *JHEP* **05** (2017) 091, arXiv:1703.00501 — **verified**; θ=π discrete 't Hooft anomaly mixing time reversal with 1-form center symmetry; θ=π vacuum cannot be trivially gapped.
- **Cheeger–Müller theorem** (analytic torsion = Reidemeister torsion), Müller *JAMS* **6** (1993) 721 — **verified**; normalized torsion metric-independent ⇒ rigid.
- **Minhyong Kim**, "Arithmetic gauge theory: a brief introduction," *Mod. Phys. Lett. A* **33**(29) (2018) 1830012, arXiv:1712.07602 — **verified**; arithmetic Chern–Simons/BF, primes-as-knots, no mass-gap result.
- Atiyah–Patodi–Singer (1975/1976); Ray–Singer (1971); Chamseddine–Connes (*CMP* **186** 1997); Connes–Marcolli (2008); Gurau (*CMP* **304** 2011); Costello–Gwilliam (CUP 2017/2021); Barrett–Crane (*JMP* **39** 1998); Nekrasov (*ATMP* **7** 2003); FHT (*Ann. Math.* **172** 2010) — standard, cited from primary venues [R].
- **Rejected as non-peer-reviewed / crank-adjacent** (found in search, not used as evidence): several Zenodo/ResearchGate/preprints.org "mass gap proofs" invoking eta invariants, Selberg–Hilbert–Pólya, "Langlands," "Kaleidoscopic Partition Filter," arithmetic Laplacian, etc. None verified, several with implausible claims. Not cited. [J]

---

## 5. Honest verdict

- **Clay still 0%.** [J] This sweep adds **no** new lower bound, no new non-positivity material, no acid-test-passing gap tool.
- **True trophies of the sweep** (organizational, not theorem-strength): [J]
  1. The **death map is closed** under all five wildcard domains.
  2. **A3 (Fredholm determinant / `D_O(z)` zeros) is the unique attractor** of every "deep" atom-seeing wildcard (Weil/Connes-sorption/Ray–Singer/analytic-torsion/Selberg) — confirming it, not the positivity or index faces, as the frontier object.
  3. The **Douglas-arithmetic premise is false**; his actual lead (stochastic quantization) is the existing G1/L3 line.
- **The wall is unmoved:** still `(H*)` = edge atom `f₁²>0` + isolation at `α_s~O(1)`, equivalently `(L-W)` / A3-zero-in-continuum / H1. No wildcard domain touches it. [O]
- **No QED; no false hope.** If a future handle exists in these domains, it must crack **(L-W)** — a *provable* flux-tube-side positivity — which none of the surveyed mathematics provides.

*Refs.* This file; registry `docs/new-tools-registry.md` §A3 (A3 identity), §B3 (rigidity), §B6 (complete monotonicity); `docs/BRIEF-FOR-MODEL.md` §3–4; `research/unconditional/TRANSPORT-FIVE-WAY-VERDICT.md`. External: citations in §4.
