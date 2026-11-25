# Z2-TQFT — the TQFT / geometric-Langlands / cohomology paradigm: a **no-go** for the 4D mass gap

> **2026-11-24 · Wave Z2 — the topological/cohomological paradigm.** The task: can a
> topological or cohomological invariant of the **moduli space of connections** *force* a positive
> 4D pure-SU(N) mass gap `m₁ > 0`? Test on 2D/3D TQFT (Chern–Simons, BF, 2D Yang–Mills). Acid
> test: free massless must fail. Kill criteria: (i) invariant **metric-independent** (cannot see a
> mass); (ii) it only gives **topological order**; (iii) it **reduces to H1**. Then the coordinator's
> sharpening (§3): does a **moduli-space cohomology class** whose nonvanishing in the continuum limit
> imply a positive **Hamiltonian spectral gap**, *without* probability/variance estimates?
> **Witness:** `research/v3/w7/z2_tqft.py` (`py -3.12 -I`, numpy 2.4.6, deterministic; **14/14 PASS**,
> exit 0; results in `z2_tqft_results.json`). Panels P1–P5 below.
> **Reads.** `research/lit-v3/LIT-ERA-03.md`, `LIT-ERA-12.md`, `docs/paper-v2/ch07-nogos.tex`
> (§rigidity), `Q6-H1-STATUS.md`.
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED)

The paradigm is a **no-go** for the 4D mass gap, and it fails **all three** stated kill criteria
simultaneously. There is no order of argument in which a topological/cohomological invariant of
`𝒜/𝒢` produces `m₁ > 0`.

1. **[NG] Metric-independence (kill i).** A cohomology class of a topological space is a
   **metric- and coupling-independent** object: it is invariant under the connected `β`-family that
   runs interacting YM into the free Gaussian at `β = ∞` (ch07 §rigidity). A mass is a dimensionful,
   metric- and coupling-*dependent* spectral threshold. The invariant is the **same** at `β = ∞`
   (where `m = 0`) as at `β = 2.5`; it cannot see a mass. Witness **P1**: the Witten index is a
   constant `1` across a deformation that changes the spectral gap by `×16`.

2. **[NG] The operator content of a class is a zero mode, not a gap (kill ii).** The only rigorous
   bridge from a cohomology class to a spectrum is Hodge/Morse/Witten: **class ↔ harmonic form ↔
   `ker(Δ)` ↔ zero mode**. A nonvanishing class therefore certifies a **ground state** (a vacuum,
   a topological sector, a ground-state degeneracy) — the **opposite** of a gap *above* the vacuum.
   "Nonzero cohomology ⇒ positive gap" is a category error. Witness **P4**: a nonzero generator of
   `H¹(T²)` (`b₁ = 2`, constant) coexists with a spectral gap that a metric deformation drives to `0`.

3. **[NG] Reduces to H1 (kill iii).** The two realizations in which "cohomology gives an exact mass"
   actually works — **Seiberg–Witten** BPS masses `= |∫_{H₁} λ_SW|` and **geometric Langlands**
   (local systems `= Hom(π₁,G)/G`, i.e. `H¹` data) — are (a) supersymmetric, (b) **gapless at the
   symmetric point** (monopoles/dyons become massless at `u = ±Λ²`), and (c) **literally `H₁`/`H¹`**.
   As a pure-YM gap tool they are the campaign's **H1 cell**, unsolved. The dichotomy of the
   coordinator's question (§3) is a false one: the invariant is metric-independent (kill i) **and**,
   if forced to speak about the gap, it is the transfer-matrix gap (kill iii). No third option.

**Acid test.** The invariant is **constant** across free↔interacting (P1, P2), so it does **not**
"fail on the free field" in the discriminating sense — it is blind at *both* ends. The only object
that even resembles a topological mass is the 3D Deser–Jackiw–Templeton mass `m = μe²` (§7, P3): it
**requires the metric Maxwell coupling `e²`**, vanishes as `e → 0` (correctly), is **tree-level**, and
exists **only in odd dimensions** — no 4D analogue for pure YM.

**Answer.** No topological/cohomological invariant of the moduli space of connections forces a
positive 4D mass gap. **Kill criteria (i), (ii), (iii) all triggered.**

---

## 1 · What the paradigm actually names

Three inequivalent objects are conflated under "TQFT / geometric Langlands / cohomology":

| Object | What it is | Metric/coupling dependence |
|---|---|---|
| **Characteristic / equivariant class** on `𝒜/𝒢` | e.g. `∫tr F∧F` (instanton number), the Atiyah–Bott equivariant classes, Chern classes `∈ H^•(BG)` | **none** (topological) |
| **3D TQFT** (Chern–Simons, BF) | exactly solvable; no local dof; WRT/Jones invariants; ground-state degeneracies | **none** (topological) |
| **Geometric Langlands** (Kapustin–Witten 2007) | topologically twisted `𝒩=4` on `Σ×C` → 2D sigma model into the **Hitchin moduli**; S-duality ↔ GL; local systems | **none** *after the twist* (topological) |
| **Seiberg–Witten** BPS spectrum | exact masses `= |Z|`, `Z` = period over an `H₁` cycle of the SW curve | BPS (protected); SUSY |

Every entry shares the defining property: **the object is a class/invariant of a topological space or
a topological field theory, hence independent of the metric and the coupling.** That single property
is what makes the whole paradigm blind to a mass — and it is a *theorem-level* statement, not a
failure of cleverness. The rest of this note makes it precise in each of the three guises and adds
the coordinator's Hamiltonian-gap question.

A note on the *moduli space of connections* specifically. In 2D it is `ℳ_flat(Σ,G) = Hom(π₁Σ,G)/G`
(moduli of flat connections = local systems); in 4D the analogous objects are the instanton moduli
`{F⁺ = 0}/𝒢` and the equivariant cohomology `H^*_𝒢(𝒜)` computed by Atiyah–Bott via the
Yang–Mills Morse function. Atiyah–Bott show the YM functional over a Riemann surface is **perfect**:
its Morse theory yields the **Betti numbers** of `𝒜/𝒢` — pure topology, metric-independent
(`P_t(N₀(2,1)) = ((1+t³)^{2g} − t^{2g}(1+t)^{2g})/((1−t²)(1−t⁴))`). Perfectness is exactly the
statement that all the Morse data collapses to topological invariants; nothing spectral survives.

---

## 2 · Kill (i): metric-independence, and why indices cannot carry a mass

**Theorem (rigidity; ch07 §rigidity).** *Any index `I(β) ∈ ℤ` that is locally constant on `(0,∞)`
and defined by a homotopy invariant of the gauge-field configuration space satisfies `I(β) = I(∞)`,
the free-field value, for every `β` along the Wilson family.* `[T]`

The Wilson family `β ↦ μ_β` is continuous in `β` and connects every finite `β` to the Gaussian limit
`β → ∞`. A homotopy invariant is unchanged by continuous deformation; in the Gaussian limit the only
index available is the trivial one. Hence **no discrete invariant can be positive for YM and zero for
free YM** — it is the *same* at both. This is the precise reason "the gap is not topological": the gap
must be a **`β`-non-analytic jump**, and a jump is the one thing a deformation-invariant cannot encode.

**Witness P1 (Witten index vs gap) `[T]+[V]`.** `𝒩=2` SUSY QM with `H_B = ½p² + ½(W')² − ½W''`,
`H_F = ½p² + ½(W')² + ½W''`, `W = ½λx²`. Discretized (`N = 800`, Dirichlet) and diagonalized:

| `λ` | `n_B(E<λ/2)` | `n_F(E<λ/2)` | **index** | **gap** `= E_F0` |
|---|---|---|---|---|
| 0.15 | 1 | 0 | **1** | 0.152 |
| 0.30 | 1 | 0 | **1** | 0.300 |
| 0.60 | 1 | 0 | **1** | 0.600 |
| 1.20 | 1 | 0 | **1** | 1.200 |
| 2.40 | 1 | 0 | **1** | 2.400 |

The **index is the constant integer `1`** (a genuine cohomology/deformation invariant — the dimension
of the `Q`-cohomology), while the **spectral gap equals `λ`** and spans a factor `15.8` across the
same deformation. The Witten index (Witten 1982, *Supersymmetry and Morse theory*: `index =
Σ (−1)^p b_p` = Euler characteristic, a topological invariant) is **exactly the quantity that cannot
change**, and the gap is **exactly the quantity that does**. A symmetric double well is the other
topological value (`index = 0`); the gap is never the index.

> **The mechanism.** Cohomology is a deformation invariant; the gap is a deformation-sensitive
> threshold. The Witten/Morse Laplacian `Δ_t = e^{−tf} d e^{tf} d^* + …` has `ker Δ_t ≅ H^*` for every
> `t` (the index), but its **first nonzero eigenvalue** `λ₁(t)` — the spectral gap — depends on `t`
> (the metric). The class controls `ker`; the gap is *outside* `ker`.

---

## 3 · [Coordinator] Moduli-space cohomology class vs Hamiltonian spectral gap

> **Question.** Is there a **topological invariant** (a characteristic class, or an **equivariant
> cohomology class** on `𝒜/𝒢`) whose **nonvanishing in the continuum limit implies a positive mass
> gap** — *without* any probability/variance estimate? Or is such a class metric-independent (cannot
> see a mass), or does it reduce to the transfer-matrix gap?

**Answer: no such invariant exists, and both horns of the dichotomy bite at once `[NG]`.**

### 3.1 The dichotomy is not exclusive — it is exhaustive-and-both

Let `ω ∈ H^k(𝒜/𝒢)` be a characteristic or equivariant class with `[ω] ≠ 0` in the continuum limit.
Two — and only two — things can be meant by "`[ω] ≠ 0` forces `m₁ > 0`":

- **(A) `ω` is used as a topological invariant.** Then `ω` is invariant under every metric deformation
  and under the connected `β`-family. It is the same class at `β = ∞` (`m = 0`). It **cannot** be the
  metric-dependent number `m₁ > 0`. Kill (i).
- **(B) `ω` is used as an operator (Hodge/Morse) datum.** Then `ω` corresponds to a **harmonic form**,
  i.e. an element of `ker(Δ)`, i.e. a **zero mode**. A zero mode is a **ground state / vacuum
  degeneracy / topological sector** — the **opposite** of a gap *above* the vacuum. To turn the class
  into a *gap* one must use the **spectral** (Hodge/transfer) Laplacian, whose gap `λ₁` is
  **metric-dependent**; at that point the statement is **no longer topological** and is precisely the
  **transfer-matrix gap** (= campaign H1). Kill (iii).

So the class is metric-independent (A ⇒ cannot see a mass) **and**, if made to speak about the gap,
collapses to the transfer-matrix gap (B ⇒ reduces to H1). The question's dichotomy is a **false
dichotomy**: there is no third option in which a *topological* class yields a *metric* gap.

### 3.2 Witness P4: a nonvanishing class coexisting with a vanishing gap `[T]+[V]`

The finite-dimensional stand-in for "cohomology of `𝒜/𝒢` vs the YM Hamiltonian gap" is the **squashed
torus** `T²(s)` (circumferences `1, s`). Its cohomology is a genuine nonvanishing class, metric-
independent; its spectral gap is metric-dependent:

| `s` | `b₀` | `b₁` | `b₂` | `λ₁(s) = 4π² min(1, 1/s²)` | `λ₁·s²` |
|---|---|---|---|---|---|
| 1 | 1 | **2** | 1 | 39.478 | 39.478 |
| 2 | 1 | **2** | 1 | 9.870 | 39.478 |
| 4 | 1 | **2** | 1 | 2.467 | 39.478 |
| 8 | 1 | **2** | 1 | 0.617 | 39.478 |
| 16 | 1 | **2** | 1 | 0.154 | 39.478 |

The generator of `H¹(T²)` is **nonzero and constant** (`b₁ = 2`) for every `s`, while the spectral gap
`λ₁(s) → 0` under the metric deformation (`λ₁·s² = 4π²` exactly). **A nonvanishing cohomology class
does not force a positive gap.** The class is a *vacuum* datum; the gap is the first *non*-zero
eigenvalue.

The YM reading is exact. On `T²` the **flat-connection moduli** `Hom(ℤ²,G)/G` is the **zero-energy
(vacuum) sector**; its cohomology labels the **ground states**; the gap to the first *non-flat*
connection is set by the Yang–Mills action `∝ g² C₂` — the coupling-dependent scale, not the class
(§5, P2).

### 3.3 Why "without a variance estimate" cannot rescue it

The coordinator's phrasing ("without probability/variance estimates") correctly identifies that the
campaign's H1 is a *variance* bound (`χ = ∫s⁻¹dν ≤ 1/λ`). The hope is that a **cohomological** route
bypasses variance. It cannot, because the missing information is not of the variance type at all: it
is the **metric-dependent spectral threshold** itself. Any route that extracts `m₁` from a class must
introduce the metric — and the cheapest metric functional that carries the threshold is exactly the
transfer-matrix/H1 object. There is no metric-free carrier of a dimensionful mass. **Kill (iii).**

---

## 4 · Kill (ii): topological order ≠ a mass gap; anomalies forbid the wrong thing

### 4.1 The cohomological *classification* of gapped phases **presupposes** the gap

The modern statement that "gapped phases are classified by generalized cohomology / cobordism"
**Freed–Hopkins 2021** (`arXiv:1604.06527`, *Geom. Topol.* 25, 1165) is a classification of
**invertible topological field theories** under the standing principle **"gapped ⇒ TQFT."** The input
is *spacetime dimension and symmetry group*; the output is the group of deformation classes
(Thom/Madsen–Tillmann bordism, Anderson duals). It tells you, **given that a phase is gapped**, which
SPT/topological-order class it is. It is **silent on whether the microscopic theory is gapped at all**.
Reading it backwards — "a nonzero class ⇒ a gap" — is a category error: the classifies-a-gap arrow
points the other way. `[NG]` (kill ii).

### 4.2 't Hooft anomalies force the *complement* of a trivially-gapped phase

An 't Hooft anomaly (generalized symmetries: **Gaiotto–Kapustin–Seiberg–Willett 2015**,
`arXiv:1412.5148`, *JHEP* 02 (2015) 172) forbids a **trivially gapped** IR — a unique, symmetric,
non-degenerate vacuum. Its allowed alternatives are: **spontaneous breaking, gaplessness, or
topological order.** An anomaly can force "topologically ordered" (which is *gapped but not trivial*),
but it **never forces a positive gap** — its other branch is *gapless*. So the anomaly logic has the
**wrong sign** for the gap question. `[NG]`

For **4D pure YM at `θ = 0`** there is, in fact, **no 0-form 't Hooft anomaly** to invoke: pure YM has
no fermions, and the archetypal nonperturbative global anomaly — **Witten, *An SU(2) anomaly*, Phys.
Lett. B117 (1982) 324** — requires an **odd number of Weyl doublets** and `π₄(SU(2)) = ℤ₂`; for
`N > 2`, `π₄(SU(N)) = 0`, and pure YM has no fermions at all. (ch07 already records: "no 0-form 't
Hooft anomaly at `θ = 0`; blind.") So even the anomaly route has nothing to stand on here. `[T]`

### 4.3 3D Chern–Simons / BF: a topological number with no mass parameter

**Chern–Simons** (Witten 1989, *Quantum field theory and the Jones polynomial*, *Comm. Math. Phys.*
121, 351) is an exactly solvable TQFT: **no local degrees of freedom**, so there is no local field to
carry a mass. Its invariants (WRT/Jones; ground-state degeneracy on `T² = |k|`) depend on the level
`k` only — there is **no mass parameter** in the theory. Asking "is the CS mass gap positive" is
**ill-posed**: a TQFT has no propagating local mode whose gap could be measured. **BF theory** is the
same (partition function = Ray–Singer torsion / flat-connection count; metric-independent). These are
2D/3D **specimens of the failure**, not sources of a 4D gap. `[NG]` (kill ii)

**Witness P2 (2D Yang–Mills, exact) `[T]+[V]`.** 2D YM is exactly soluble (Migdal 1975; Witten 1991,
*Comm. Math. Phys.* 141, 153) and makes the split explicit. For SU(2), `E_j = g² C₂(j)`,
`C₂(j) = j(j+1)`:

| `g²` | gap `= E_{1/2}` | shape `E_j/g²` |
|---|---|---|
| 0.05 | 0.0375 | 0, 0.75, 2.00, 3.75, 6.00 |
| 0.20 | 0.1500 | 0, 0.75, 2.00, 3.75, 6.00 |
| 1.00 | 0.7500 | 0, 0.75, 2.00, 3.75, 6.00 |

The **topological data** — which representations exist, the Casimir *shape* `E_j/g² = C₂(j)` — is
**coupling-independent**. The **gap** is the coupling-**dependent** scale `g² C₂(½) = ¾g²`, linear in
`g²`, and `→ 0` as `g → 0`. **Witten 1991 makes this sharp:** the topological (weak-coupling /
zero-area) limit of 2D YM **is** the Reidemeister–Ray–Singer torsion — i.e. the topological invariant
**lives at the free, gapless end**. The gap is the *deviation* from the TQFT at finite coupling.

---

## 5 · Kill (iii): the exact "cohomology gives a mass" instances are all `H₁` and all SUSY

There *are* settings where cohomology gives an **exact** mass. Each is instructive because each
**triggers kill (iii)** and is **SUSY / protected**, hence inapplicable to pure non-SUSY YM.

- **Seiberg–Witten `[T]` (Nucl. Phys. B426 (1994) 19, `hep-th/9407087`).** The exact BPS mass is
  `M = √2 |Z|`, `Z = n_e a + n_m a_D = ∮_{γ ∈ H₁(Σ_SW)} λ_SW`: the central charge is the **period of
  the Seiberg–Witten differential over a one-cycle** — literally `H₁` data. But: (i) it is
  `𝒩=2`-SUSY; (ii) the lightest BPS states become **massless** at the `SU(2)`-symmetric points
  `u = ±Λ²` (monopole/dyon condensation) — the theory is **gapless there**; (iii) the masses are
  **charged, Coulomb-branch** states, not the `0⁺⁺` glueball. Pure YM is non-SUSY and has no such
  moduli. **All three kill criteria at once.** `[NG]`

- **Geometric Langlands `[T]` (Kapustin–Witten 2007, `hep-th/0604151`, *Commun. Num. Theor. Phys.*
  1, 1).** The GL correspondence relates Hecke eigensheaves on `Bun_G(Σ)` to **local systems** — flat
  `G`-connections, i.e. points of `Hom(π₁Σ, G)/G`, the `H¹` (nonabelian `H¹`) data. The physical
  statement is an **S-duality of the topologically twisted `𝒩=4` theory**, whose twisted observables
  are **metric-independent**. `𝒩=4` is **conformal** (no mass gap) and the twist makes the theory a
  TQFT. So GL is about a *gapless* SUSY theory and about *topological* (`H¹`/local-system) data. It
  computes things like the Jones polynomial / the GL equivalence — `β`-independent — not a 4D
  pure-YM mass. `[NG]`

- **The campaign's `H1`.** `Q6-H1-STATUS.md`: `H1: χ := Var_μ(ξ_{k→0}) = ∫s⁻¹dν ≤ 1/λ`, a
  variance/susceptibility bound, *strictly weaker* than the gap and **satisfied by gapless models**.
  The cohomological route does not prove this inequality and does not even reach its level: it gives
  either a `β`-independent class (kill i) or the transfer gap itself (kill iii). **No independent
  content.** `[NG]`

The pun is exact and load-bearing: the TQFT/GL/cohomology data is **`H¹`/`H₁`** (local systems, cycle
periods) in the literal topological sense, and as a 4D gap tool it is the campaign's **H1 cell** in the
dynamical sense. Either way, **"reduces to H1."**

---

## 6 · The 2D/3D TQFT acid-test table

| Theory | Dim | Topological invariant | Gap? | Verdict for the paradigm |
|---|---|---|---|---|
| Chern–Simons `U(1)_k` | 3 | WRT/Jones; deg. on `T² = \|k\|` | **no local dof** (ill-posed) | specimen of failure; invariant is `k`, not a mass `[NG]` |
| BF | 3,4 | Ray–Singer torsion; flat-conn count | **no local dof** | specimen; metric-independent `[NG]` |
| 2D Yang–Mills | 2 | rep labels / Casimir shape (P2) | `gap ∝ g² C₂`, `→0` at free | topological data at the **free** end `[T]+[V]` |
| Maxwell–CS (DGT) | 3 | CS level `μ` | `m = μe²` **needs `e²`** | only "topological mass"; §7 `[T]+[V]` |
| 4D pure YM | 4 | instanton `ℤ`, `H^*(𝒜/𝒢)` | **target** `m₁ > 0` | invariant is `β`-indep.; blind `[NG]` |

Every row fails the acid test in the same way: the invariant is **constant** across free↔interacting,
or the "mass" is a **tree-level parameter** requiring the metric coupling.

---

## 7 · The strongest steelman: the 3D Deser–Jackiw–Templeton "topological mass"

The best case for "a topological term generates a mass" is 3D **Maxwell–Chern–Simons**, whose photon
is genuinely massive (Deser–Jackiw–Templeton 1982, *Phys. Rev. Lett.* 48, 975). The transverse
inverse-propagator block is `M(ω) = a·I + b·J`, `a = −(ω²−k²)/e²`, `b = iμω`, `J = [[0,1],[−1,0]]`;
`det M = a² + b² = 0` gives the physical pole `ω² − μe²ω − k² = 0`, so at `k = 0`, **`m = μe²`**.

**Witness P3 `[T]+[V]`** (nine `(μ, e²)` points): `m = μe²` to machine precision; `m/(μe²) = 1.000`
throughout; `m = 0` if `μ = 0` (no CS — gapless free Maxwell, *acid test fails correctly*) and `m = 0`
if `e² = 0` (**the topological term alone gives no mass**); `m` linear in `μ`.

**Why this is not a counterexample to the paradigm.** `m = μe²`:
- **carries the metric coupling `e²` as a factor** — it is *not* a function of a topological invariant
  alone; it needs the metric Maxwell term. Remove the metric (`e → 0`) and the mass goes to zero.
- is a **tree-level parameter** (the CS coefficient `μ`), not a dynamically generated spectral gap; the
  DGT literature stresses the *quantization* `μg²/4π ∈ ℤ` of the dimensionless ratio — still a
  parameter, not a threshold above a vacuum.
- exists **only in odd dimensions** (`ε^{μνρ}` needs `d` odd; the 4D analogue `F∧F` is a total
  derivative and gives **no** gauge-boson mass). **There is no 4D "topological mass" for pure YM.**
  `[NG]`

This is the honest best case, and it fails the test on **four** counts (metric coupling, tree-level,
odd dimension, parameter-not-threshold). It does not transfer.

---

## 8 · The constructive no-go: the only salvage leaves the paradigm (the exact missing step)

Per charter §1 we tried to *build* a topological/cohomological tool that sees the gap. The attempt
pins the exact missing step `[O]`:

- **Try: spectral (`L²`-Hodge) cohomology of the instanton moduli with the induced metric.** The
  lowest Hodge-Laplacian eigenvalue on `ℳ_inst` could be the gap. But `ℳ_inst` is **noncompact**
  (instanton size `ρ ∈ (0,∞)`), so the Laplacian has **continuous spectrum down to `0`** — no gap. On
  any compact piece, `L²` cohomology `≅` ordinary Betti numbers (Hodge) — **topological again**.
  Either fails noncompactness or reverts to topology. `[NG]`
- **Try: a metric-weighted cohomology.** To break the deformation-invariance you must weight the class
  by a **positive functional of the metric** `w(g)`. The moment you do, the object is **no longer
  topological**: it is a metric functional whose positivity is exactly the **W / H1 / atom** object
  (`w₁ = f₁²/M₀ ≥ w₀`, or `χ ≤ 1/λ`). This is the **exact missing step made explicit**: the only way
  a cohomological bread carries the gap is to add a metric weight, and the resulting object *is* the
  campaign's open H1/atom cell, not a topological invariant.

So the paradigm is **binary**: either the object is topological (**blind**, kill i) or it is a metric
functional (**H1/atom**, kill iii). There is no intermediate object that is *both* topological *and*
gap-carrying. This is a **constructive** no-go, not merely an adjudication.

---

## 9 · Relation to the existing no-go ledger (ch07)

This wave **re-derives and sharpens** two ch07 entries from the topological side:

- **Rigidity no-go (ch07 §rigidity).** Panel P1 gives it a number: index `= 1` constant, gap `×16`.
  The rigidity theorem is the reason "the gap is not topological."
- **"Topology. Dead by rigidity."** §5 adds the *positive* diagnosis: the cohomology route does not
  merely fail — it *would* work only in the SUSY/`H₁` corner (Seiberg–Witten), which is the wrong
  theory (SUSY, gapless at the symmetric point, wrong channel). This is the exact reason the corner is
  empty for pure YM.
- New: **the moduli-space cohomology-class vs Hamiltonian-gap question (§3)** is closed with the
  both-horns argument (metric-independent **and** reduces to the transfer gap) plus the `T²(s)` witness.

---

## 10 · Ledger and references

| # | Item | Claim used | Status |
|---|---|---|---|
| W-Morse | Witten 1982, *Supersymmetry and Morse theory*, *J. Diff. Geom.* 17, 661 (DOI 10.4310/jdg/1214437492) | index = Euler char. = topological; `Δ_t` gap is metric | **VERIFIED** |
| AB | Atiyah–Bott 1983, *The Yang–Mills equations over Riemann surfaces*, *Phil. Trans. R. Soc. A* 308(1505), 523 (DOI 10.1098/rsta.1983.0017) | YM Morse is perfect ⇒ **Betti numbers** (topological), not a spectrum | **VERIFIED** |
| KW | Kapustin–Witten 2007, *Electric-magnetic duality and the geometric Langlands program*, *Commun. Num. Theor. Phys.* 1, 1 (`hep-th/0604151`, DOI 10.4310/CNTP.2007.v1.n1.a1) | twisted `𝒩=4` → Hitchin sigma model; local systems = `H¹`; topological | **VERIFIED** |
| SW | Seiberg–Witten 1994, *Electric-magnetic duality, monopole condensation…*, *Nucl. Phys.* B426, 19 (`hep-th/9407087`) | BPS mass `= \|∮_{H₁} λ\|`; **gapless** at `u=±Λ²`; SUSY, charged | **VERIFIED** |
| FH | Freed–Hopkins 2021, *Reflection positivity and invertible topological phases*, `arXiv:1604.06527`, *Geom. Topol.* 25, 1165 | classifies gapped/invertible phases under **"gapped ⇒ TQFT"**; *presupposes* gap | **VERIFIED** |
| GKSW | Gaiotto–Kapustin–Seiberg–Willett 2015, *Generalized global symmetries*, *JHEP* 02 (2015) 172 (`arXiv:1412.5148`) — **four** authors (no Komargodski) | 't Hooft anomaly forbids trivially gapped ⇒ gap **or** topological order | **VERIFIED** |
| W-SU2 | Witten 1982, *An SU(2) anomaly*, *Phys. Lett.* B117, 324 (DOI 10.1016/0370-2693(82)90728-6) | `π₄(SU(2))=ℤ₂` anomaly needs **fermions**; `π₄(SU(N))=0` `N>2`; pure YM none | **VERIFIED** |
| W-CS | Witten 1989, *Quantum field theory and the Jones polynomial*, *Comm. Math. Phys.* 121, 351 (DOI 10.1007/BF01217730) | CS exactly solvable TQFT; no local dof | **VERIFIED** |
| Migdal | Migdal 1975, *Recursion equations in gauge field theories*, *Sov. Phys. JETP* 42, 413 | 2D YM exact | **VERIFIED** |
| W-2D | Witten 1991, *On quantum gauge theories in two dimensions*, *Comm. Math. Phys.* 141, 153 (DOI 10.1007/BF02100009) | weak-coupling limit = **Ray–Singer torsion** (topological = free end) | **VERIFIED** |
| DGT | Deser–Jackiw–Templeton 1982, *Three-dimensional massive gauge theories*, *Phys. Rev. Lett.* 48, 975 (DOI 10.1103/PhysRevLett.48.975) | `m = μe²`; tree-level; 3D (odd) only | **VERIFIED** |

**Correction caught:** the generalized-symmetry paper has **four** authors — Gaiotto, Kapustin,
Seiberg, Willett — **not** "Gaiotto–Kapustin–Komargodski–Seiberg–Willett"; there is no Komargodski on
this paper. Do not cite a fifth author.

Witness: `research/v3/w7/z2_tqft.py` — **14/14 PASS**, exit 0; `z2_tqft_results.json`.
Panels: **P1** SUSY-QM index `=1` constant vs gap `×16`; **P2** 2D-YM rep data constant vs gap `∝g²`;
**P3** DGT `m=μe²` (needs `μ` **and** `e²`); **P4** `T²(s)` cohomology `b₁=2` constant vs `λ₁→0`;
**P5** ledger.

---

## 11 · One-paragraph handoff

The topological / geometric-Langlands / cohomology paradigm is a **no-go** for the 4D mass gap,
because a cohomology class `is` a metric- and coupling-independent object, and a mass is a metric- and
coupling-dependent threshold. It fails the three stated kill criteria simultaneously: it is
**metric-independent** (blind — constant across the free↔interacting `β`-family; P1, P2); its only
operator content is a **zero mode / vacuum / topological sector**, the opposite of a gap *above* the
vacuum (P4; the coordinator's nonvanishing-class question is closed — a nonzero `H¹` class coexists
with a tunable-to-zero gap, and turning the class into a gap uses the metric Hodge/transfer Laplacian,
i.e. it *reduces to the transfer gap*); and the only exact "cohomology gives a mass" instances —
Seiberg–Witten `= H₁` periods, geometric Langlands `= H¹` local systems — are SUSY, gapless at the
symmetric point, and *literally* `H₁`, i.e. the campaign's **H1 cell**, unsolved. The single
"topological mass" (3D DGT `m = μe²`) needs the metric coupling `e²`, is tree-level, is odd-
dimensional, and has no 4D analogue. Constructively, the only way to make a cohomological object carry
the gap is to weight it by the metric — at which point it is the W/H1/atom cell, not topology. **The
paradigm has no independent content: it is either blind or it *is* H1.** No QED. Clay 0%.
