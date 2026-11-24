# Q1-INFRARED-BOUND — a non-Abelian reflection-positivity infrared bound for the `k=0` susceptibility `χ`: the four standard machines, the exact place each dies, and the missing lemma

> **2026-11-24 · Wave Q1.** The assigned target is H1 — the last open input of the campaign:
> `χ = Var_μ(ξ_{k→0}) = ∫s⁻¹dν(s) ≤ 1/λ` uniformly in `(L,a)` at `β=2.5`, equivalently the
> constant-mode curvature `G″_{k=0} ≥ λ > 0` (Gaussian marginal `1/χ`). **Question:** can the
> rigorous **reflection-positivity (RP) infrared-bound** machinery (Fröhlich–Simon–Spencer;
> Brydges–Fröhlich–Seiler; Osterwalder–Seiler; Balaban) be made **non-Abelian** and reach `β=2.5`?
> Reads (repo): `P-WAVE-SYNTHESIS.md`, `D6-DECISION.md`, `H6-HEAD-SPECTRAL.md`,
> `S1-CONDITIONAL-POINCARE.md`, `S6-DOBRUSHIN.md`, `C2-BALABAN-UNIFORM.md`,
> `P6-PERIODIC-TO-H1.md`, `docs/paper-v2/ch02-squeeze.tex`, `docs/paper-v3/main.tex` §6.
> Witness: `research/v3/w7/q1_infrared_bound.py` (`py -3.12`, numpy 2.4.6 + scipy 1.18.1; output §6).
> **RP machinery is treated as standard knowledge; no third-party sources downloaded.**
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. External textbook results
> are marked `[std]` inline (used, not re-verified this session). No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

**No RP infrared bound reaches `β=2.5`.** This is a `[NG]` for the **method class**
{RP Gaussian domination, reflected/chessboard positivity, compact-transfer / Perron–Frobenius,
heat-kernel / log-Sobolev, Dobrushin / influence}, and the failure is not diffuse — it localizes to
**one missing inequality at the constant mode**. Four independent faces, each a theorem:

1. **FSS Gaussian domination is vacuous at `k=0` `[T]`.** The RP infrared bound compares the
   interacting two-point function to the **free Gaussian** covariance: `⟨φ(f)²⟩_c ≤ (f, C_free f)`
   after Duhamel + RP. For `f = 𝟙` (the constant mode) the free massless covariance is the zero
   mode of `−Δ`, and `(𝟙, C_free 𝟙) = Ĉ_free(0) = ∞`. The bound reads `χ ≤ ∞`. The FSS theorem is
   *stated* for `f ⊥ constants` (its bound controls the gradient / non-constant / spin-wave modes,
   exactly to avoid this); it cannot see `k=0`. **[T]**

2. **The free photon is an RP object with `χ = ∞`, and the acid test forbids bounding it `[T]`.**
   The free Maxwell measure is Gaussian and reflection-positive (it satisfies **every** hypothesis
   of FSS), yet its constant mode is a pure-gauge flat direction: `χ = Ĉ(0) = ∞` (witness §6·A:
   `χ = 1/m² → ∞` as the IR mass `m→0`, independent of `L`). Any RP bound that returned a finite
   `χ` for non-Abelian SU(N) without using non-linearity would also return one here — contradiction.
   **The bound must use compactness/non-linearity, and RP alone cannot supply it.**

3. **Compactness does make `β=0` finite and `1D` gapped, exactly `[T]+[V]`** — but only there.
   `Var_Haar(ξ₃) = (4π²−6)/9 = 3.7198` (witness §6·B) is finite **purely because SU(2) is compact**;
   the 1D compact transfer operator is Hilbert–Schmidt, hence gapped by Perron–Frobenius
   (`r = f̂_{1/2}/f̂_0 < 1`, `Δ = 0.68`, `χ_1D = 3.06` at `β=2.5`, witness §6·F). **Neither fact
   lifts to 4D at `β=2.5`:** compactness caps the variance only at `π²N_links` (`= Θ(L⁴)`, §6·C), and
   the 4D (3D-layer) transfer operator is not compact and admits massless phases.

4. **The interaction cannot be transferred uniformly from the `β=0` Haar bound `[T]`.** Since
   `0 ≤ S_W ≤ 2βN_p`, the density ratio is bounded but its sup is `e^{2βN_p} = e^{Θ(L⁴)}`
   (`2βN_p = 30L⁴` at `β=2.5`, witness §6·D): a **volume-growing** inflation, not a uniform one.
   The only volume-uniform decorrelation route is the Dobrushin/SZZ influence bound, whose threshold
   is `β < 1/48` (`α(β)=72β`, `α(2.5)=180 ≫ 1`, witness §6·E) — disjoint from the tower `[2.5,∞)`.

**One line.** *An RP bound reaches `β=0` (compactness, `χ=3.7198`) and the 1D chain (compact
transfer gap), but not `β=2.5` in 4D: FSS is vacuous at `k=0` (massless reference), the free photon
(RP, `χ=∞`) is the forbidden target, compactness gives only `O(V)`, and the interaction reweighting
is `e^{Θ(L⁴)}`.* **The exact missing lemma is a volume-uniform finite-mass Gaussian domination at the
constant mode — a non-Abelian FKG/diamagnetic inequality — which is H1.** No QED. Clay 0%.

---

## 1 · The standard RP infrared-bound machinery, stated as standard `[std]`

### 1.1 Fröhlich–Simon–Spencer (1976), *Infrared bounds, phase transitions and continuous symmetry breaking*, CMP **50** 79–95 `[std]`

The abstract RP infrared bound has two steps.

**(S1) Gaussian domination / Duhamel bound.** Let `μ₀` be a centred Gaussian measure on `ℝ^Λ` with
covariance `C`, reflection-positive with respect to a reflection `θ`, and let
`μ = Z⁻¹ e^{−V} μ₀` with `V` a reflection-positive-compatible interaction (a sum of local terms,
each positive-definite in the RP sense — e.g. `V` is block-positive, or `V≥0` with RP). Define the
**Duhamel two-point function**
```
   D_μ(f,g) := ∫₀¹ dt ⟨ Φ(f) e^{−tV} Φ(g) e^{−(1−t)V} ⟩_μ .
```
Then **Gaussian domination** gives
```
   D_μ(f,f)  ≤  (f, C f)          for all real test functions f.                       (1.1)
```
**(S2) Reflection positivity / chessboard.** For an RP measure `μ`, the **physical connected**
two-point function obeys
```
   ⟨ Φ(f) Φ(f) ⟩_{μ,c}  ≤  D_μ(f,f) .                                                 (1.2)
```
Combining (1.1)+(1.2): `⟨Φ(f)²⟩_c ≤ (f, C f)`. For the lattice `|∇φ|² + P(φ)` model with
`C = (2β∇*∇)⁻¹` this is the classical
```
   ⟨ φ(f)² ⟩_c  ≤  (1/2β) ∫ f (−Δ)⁻¹ f ,         (1.3)
```
which drives Mermin–Wagner / the absence of Goldstone breaking in `d≤2` and the finite transverse
susceptibility in `d≥3`.

**Two features that decide the whole question for us.**

- **(1.3) is a statement about `(−Δ)⁻¹`, which is defined on the orthogonal complement of
  constants.** For `f = 𝟙`, `(−Δ)⁻¹𝟙` is the zero mode (ill-defined/infinite), and the physical
  `Var` of the `k=0` mode is *exactly* the divergent piece. FSS therefore bounds the
  **non-constant / gradient / spin-wave** sector, not the constant mode. **[T]** (this is the
  single most load-bearing observation in this note).
- The reference is the **free** covariance. The interacting theory is bounded **above** by the free
  one. This is the *correct direction* for Mermin–Wagner (correlations too large ⇒ no order) but the
  **wrong direction for H1**, which needs an upper bound on `Var(a)` when the free reference already
  has `Ĉ(0) = ∞`.

### 1.2 Osterwalder–Seiler (1978), *Gauge field theories on a lattice*, Ann. Phys. **110** 440–471 `[std]`

Proves that the Wilson lattice gauge measure (Haar `× e^{−S_W}`) is **reflection-positive** with
respect to hyperplane reflections and time reflection, hence has a positive-metric Hilbert space,
a self-adjoint transfer matrix, and gauge-invariant physical observables. Proves a **mass gap at
strong coupling** by a cluster/gauge-fixing expansion (small `β`). It does **not** provide an
FSS-type IR bound and makes no claim at `β=2.5`. **[std]**

### 1.3 Brydges–Fröhlich–Seiler (1979–83) `[std]`

The heat-kernel / random-surface formulation of lattice gauge theory: RP + polymer / random-surface
cluster expansion. Establishes construction and area law / screening at **strong coupling**, and the
`d=2` (and `d=3` with matter) constructions. The mechanism is an **expansion**, valid in a
small-`β`/large-mass regime; it is not an IR bound and does not reach `O(1)` coupling. The
**lattice diamagnetic inequality** `|∇|U|| ≤ |DU|` (their tool) is the one place compactness enters
*constructively* — and it is silent on the constant mode (see §2·F, §3·F3). **[std]**

### 1.4 Balaban (1985–89), CMP **102/109/116/119/122** `[std]`

UV stability of lattice SU(2) gauge theory: uniform in the cutoff `a` (number of RG steps) at
**fixed physical volume**, weak coupling `g∈(0,γ]`, at the level of the **effective action and its
densities**. Per `C2-BALABAN-UNIFORM.md` §1, this is **UV**, not IR; the volume uniformity is exactly
the open IR tower bottom (Channel I = H1), and the force-Hessian (`L3`) is not covered. **[std]**

### 1.5 The influence / log-Sobolev route (Dobrushin; Shen–Zhu–Zhu) `[std]`

Dobrushin's coefficient `α(β) = O(d)β` gives complete analyticity (hence exponential mixing and a
volume-uniform variance bound) for `α<1`. Shen–Zhu–Zhu (CMP **400** (2023) 805, arXiv:2204.12737)
give a volume-uniform LSI for `β<1/48`. Both are **strong coupling**. Per `S6-DOBRUSHIN.md`,
`α(2.5) = 180 ≫ 1`: **disjoint from the tower**. **[std]** + `[T]` (arithmetic, §6·E).

**Summary of the survey.** All four machines give either (a) a strong-coupling expansion/inequality
(`β < O(1/48)`), or (b) a UV / fixed-volume statement, or (c) an FSS bound against a **massless**
reference that is vacuous at `k=0`. None gives a finite `χ` at `β=2.5`.

---

## 2 · What fails, specifically, for non-Abelian lattice gauge and for the `k=0` susceptibility

Each failure is a named, labeled step; none is "this is hard."

### F1 — `k=0` is the free Laplacian zero mode ⇒ Gaussian domination is `∞` `[T]`

H1's object is `χ = Var(ξ_{k→0}) = Ĉ(0)` (the `k→0` **connected propagator**, KL form
`∫s⁻¹dν`). The FSS bound `⟨Φ(f)²⟩_c ≤ (f, C_free f)` at `f = 𝟙` reads `χ ≤ Ĉ_free(0)`. For the free
massless field `Ĉ_free(0) = (1/−Δ)(0) = ∞`. **The bound is vacuous by construction.** Witness §6·A:
`χ = 1/m² → ∞` as `m→0`, in every volume; the *non-constant* average `(1/#k≠0)∑_{k≠0}λ_k⁻¹` (which
FSS controls) is finite and `L`-bounded (`0.148, 0.151, 0.153` for `L=4,6,8`), confirming the
divergence is purely `k=0`.

To get a finite `χ` from Gaussian domination one needs the reference to have a **positive mass**:
`(f,C_m f) = V·(1/m²)` in finite volume, normalized to `χ ≤ 1/m²` — i.e. H1 with `λ = m²`. So
**the reference mass `m` must be the physical mass**, which is the gap. *The FSS route needs the
answer as its hypothesis.* This is not circular *abuse* — it is the exact structural reason an
upper-bound machine cannot manufacture a lower bound (the campaign's ch02 §schur "upper-bound machine
/ lower-bound target" mismatch, here in RP coordinates).

### F2 — The reference is massless and gauge-degenerate; no Proca mass `[T]`

The natural Gaussian reference for lattice gauge theory is the linearization of the Wilson action
(the Maxwell operator `−Δ δ_{μν} + ∂_μ∂_ν`). Its `k=0` spectrum is `{0}` for the pure-gauge
directions; adding a mass `m²A²` (Proca) **breaks gauge invariance and reflection positivity**. So
there is **no RP-compatible massive Gaussian reference** available without already knowing a mass.

### F3 — The constant mode of the free photon is an exact flat direction `[T]`

For the free (non-compact) Maxwell field, `A_μ = const` is pure gauge: it lies in the kernel of the
quadratic form, so its variance is `∞` even in finite volume. This is the **acid object**: RP holds
(the Maxwell measure is Gaussian RP), `χ = ∞`, and the acid test forbids any RP+free argument from
bounding it. Compact SU(N) escapes this only because the constant mode lives on a **compact** group
(§3·2), which is the *non-linearity* the task demands — but compactness at `β=0` is not yet the
`β=2.5` statement (§2·F6).

### F4 — No FKG / ferromagnetic order on the non-Abelian gauge quotient `[T]` / `[O]`

Step (S2) [physical connected ≤ Duhamel] in the scalar models is driven by RP **plus** positivity /
FKG (the ferromagnetic order). The Wilson SU(N) measure has **no total order** on the gauge quotient
`𝒜/𝒢`, and the transverse step `χ(k) ≤ χ(0)` that would carry a `k=0` bound to the tower is **not
established** (this is the campaign's own standing `[O]`, `P6-PERIODIC-TO-H1.md` §2.1: "no
monotonicity `χ(k)≤χ(0)` for 4D YM"). Abelian / special-plaquette FKG exists in restricted models
but does not give the non-Abelian constant-channel comparison at `β=2.5`. **[T]** (no order) +
**[O]** (the needed inequality).

### F5 — Savvidy indefiniteness kills the Gaussian-domination prerequisite `[T-standard]`

Gaussian domination needs the interaction to be a *positive (semidefinite) perturbation* of the
Gaussian. Expanding the Wilson action about the free measure, the quadratic form is the lattice
**Savvidy operator** `−D² + 2F`, whose constant-transverse mode is **tachyonic** (`λ_min = −g|B| < 0`)
for any nonzero chromomagnetic field (`S1` §4.1, `H6` §3; re-measured `λ_min: 0 → −1.43 → −8.01`).
The "free + positive interaction" decomposition fails at `β=2.5`; the correct normal form is a
**double well**, on which the essential infimum of `G″` is negative (`H6` §3). RP per se survives
(Wilson is RP), but the **Gaussian-domination inequality** `D ≤ C_free` does not follow.

### F6 — The interaction reweighting is `e^{Θ(L⁴)}`, not volume-uniform `[T]`

Write `μ = μ_Haar · e^{−S_W}/Z`. Since `0 ≤ S_W ≤ 2βN_p` with `N_p = 6L⁴` (4D hypercubic),
`e^{−2βN_p} ≤ e^{−S_W} ≤ 1`, so the density ratio `dμ/dμ_Haar = e^{−S_W}Z_Haar/Z ∈ [e^{−2βN_p},e^{2βN_p}]`.
Bounded-density comparison gives `E_μ[a²] ≤ e^{2βN_p}E_Haar[a²]` (hence `Var_μ(a) ≤ e^{2βN_p}(Var_Haar(a)+O(1))`),
an inflation `≤ e^{2βN_p} = e^{30L⁴}` at `β=2.5` (witness §6·D: exponent `480, 2430, 7680, …` for
`L=2,3,4`). **Volume-growing** (up to fixed `O(1)` constants). The compact `β=0` finiteness does not
transfer uniformly by boundedness alone; it needs decorrelation.

### F7 — Gauge invariance: the amplitude is not a physical observable `[T]`

RP constrains **gauge-invariant** observables (Osterwalder–Seiler). The H1 amplitude
`a = ⟨v, log U⟩` is gauge-dependent (a flat-chart coordinate). The gauge-invariant version is the
`tr F²` (or Wilson-loop) susceptibility, a **composite** operator whose zero-momentum mass is the
glueball gap. The FSS Duhamel→physical step needs the observable to be a self-adjoint multiplication
operator in the physical Hilbert space; a gauge-dependent `log U` is not. So even a perfect RP
inequality on the gauge-invariant sector controls the composite's mass — the gap — not the bare
`χ` in `log U`.

---

## 3 · The attempted non-Abelian version: what compactness and the heat-kernel gap actually give

### 3.1 The object, unchanged

`χ = Var_μ(a)`, `a = (1/√N_links)∑_ℓ ξ_ℓ`, `ξ = log U` (flat chart), equivalently
`G″_{k=0}(a₀*) = 1/χ` (Gaussian marginal; the `G″·Var=1` identity, `D6` §1, `ch02` §G″ vs `1/χ`,
measured `2–5%` defect at `β=2.5`). H1 = `χ ≤ 1/λ` uniformly in `(L,a)`.

### 3.2 Compactness at `β=0` — the exact theorem `[T]+[V]`

At `β=0` the links are i.i.d. Haar, so `a` is a normalized sum of `N_links` i.i.d. **compact**
variables and
```
   χ₀ = Var_Haar(ξ₃) = (4π²−6)/9 = 3.719824,   G″₀ = 1/χ₀ = 0.268830   (all L).            (3.1)
```
This is finite **because SU(2) is compact** (the non-compact photon gives `∞`, §2·F3). Witness §6·B
reproduces (3.1) to `3×10⁻¹⁵` by exact quadrature. **[T]+[V]**

### 3.3 Compactness alone is **not** enough — bounded range vs bounded variance `[T]+[V]`

`|ξ_ℓ| ≤ 2π` gives `|a| ≤ 2π√N_links`, hence only `Var(a) ≤ 4π²N_links = Θ(L⁴)` — **volume-growing**.
The independent (product-measure) Efron–Stein bound `Var(a) ≤ Σ_ℓ Var(ξ_ℓ|rest)/N_links` gives the
finite `σ²` **only because the links are independent** (`β=0`). For correlated links it fails: the
witness §6·C shows the same bounded range giving `Var = σ²` (independent) versus `Nσ²` (perfectly
correlated, growing like `L⁴`). **Volume-uniform `Var(a) = O(1)` requires DECORRELATION of the
constant mode — the correlation length / mass gap — not compactness.** This is the precise sense in
which "use compactness" is necessary but far from sufficient. **[T]+[V]**

### 3.4 The heat-kernel spectral gap is UV and silent on the constant mode `[T]`

The Laplace–Beltrami operator on SU(N) has a mass gap (the first Casimir); the one-link heat kernel
`p_t` satisfies a log-Sobolev / Poincaré inequality with that gap. But this is a **per-link (UV)**
gap: it controls the **fast fibre** relaxation, which `H6-HEAD-SPECTRAL.md` §2 (reading iii) and
`S1` §6 prove is **silent** on the coarse marginal curvature `G″` (two-coordinate independence: the
fast precision can be sent to `∞` with zero effect on `Var(a)`). So the heat-kernel gap is the
"fast fibre gap ~Λ_UV" — positive and volume-uniform, but it certifies nothing about `χ`.

### 3.5 The 1D compact transfer theorem — and why it does not lift `[T]+[V]`

In **1D** the SU(2) chain transfer operator `T(U,U') = e^{β·½tr(U†U')}` is **central**, so by
Peter–Weyl its eigenvalues are the SU(2) Fourier coefficients `f̂_j`. The operator is Hilbert–Schmidt
(bounded kernel on a compact group) and positivity-improving, hence **compact with a simple top
eigenvalue and a gap** (Perron–Frobenius). The fundamental decay ratio `r = f̂_{1/2}/f̂_0 < 1`, the
gap `Δ = −log r > 0`, and `χ_1D = (1+r)/(1−r) < ∞`. Witness §6·F:
```
   β      r          Δ=−log r     χ_1D
  1.0   0.24019     1.42631      1.632
  2.5   0.50720     0.67886      3.058
  4.0   0.65805     0.41848      4.849
```
This **reproduces the 1D SU(2) finite-`χ` theorem** (`P-WAVE-SYNTHESIS.md`; `P4` quotes `Δ=0.6444`,
`ξ=1.55` in its normalization). **It does not lift to 4D:** the 4D transfer operator acts on a
**3D spatial layer** and is not Hilbert–Schmidt in a volume-uniform way; the 4D system admits a
massless (free-photon) phase where the transfer operator has no gap. The 1D gap is a statement about
the compactness of a **single spatial slice's phase space**; in 4D the slice is spatially infinite
and the spectrum can reach the identity. **The lift of the 1D compact-transfer gap to 4D IS H1.**

### 3.6 The conditional statement an RP bound *would* yield `[C]`

The honest positive residue: **if** one had the finite-mass domination at the constant mode,
```
   (M)   Var_μ(a) ≤ Var_{μ_Haar}(a) + C₀ ,     C₀ independent of (L,a) ,
```
then H1 holds at `β=2.5` with
```
   χ ≤ 3.7198 + C₀ < ∞ ,     λ := 1/(3.7198 + C₀) > 0 ,                                    (3.2)
```
equivalently `G″_{k=0}(a₀*) ≥ (3.7198+C₀)⁻¹`. The P-wave data put `χ` rising from `3.72` (`β=0`)
to `≈4.45` (`β=2.5`), i.e. `C₀ ≈ 0.73` **numerically** — so `(M)` plausibly holds with a small
`C₀`. **But `(M)` with `C₀` independent of `(L,a)` is exactly H1**, not a consequence of RP. Label
of (3.2): `[C]` (the implication is `[T]`; the hypothesis `(M)` is `[O]`).

---

## 4 · The exact missing lemma (the deliverable)

> **Missing lemma (Q1).** A **volume-uniform finite-mass Gaussian domination at the constant mode** —
> equivalently a **non-Abelian FKG / diamagnetic inequality for `k=0`**: there exist `m₀>0` and
> `C₀<∞` (independent of `L,a`, `β=2.5` fixed) such that for the `k=0` amplitude `a`,
> ```
>   Var_μ(a)  ≤  Var_{μ_Haar}(a) + C₀ = 3.7198 + C₀ ,      uniformly in (L,a) ,
> ```
> equivalently `G″_{k=0}(a₀*) ≥ 1/(3.7198+C₀) =: λ > 0` — i.e. `χ < ∞` uniformly. This is H1.

**Four equivalent faces, all H1; each names a machine that fails to supply it.**

| face | statement | why absent |
|---|---|---|
| **(M1) finite-mass Gaussian domination** | FSS `D_μ(f,f) ≤ (f,C_m f)` with a *massive* `C_m`, `m²=λ>0` | the only RP-compatible `C_m` is massless (F2); `m>0` *is* the gap |
| **(M2) non-Abelian FKG / diamagnetic inequality at `k=0`** | `G″_β(a₀*) ≥ G″₀(a₀*) − C₀` | no total order on `𝒜/𝒢` (F4); false pointwise on the Savvidy flank (F5) |
| **(M3) volume-uniform decorrelation** | dependent Efron–Stein / Dobrushin `α(β=2.5)<1` | `α(2.5)=180≫1`; threshold `β<1/48` (F6, §1.5) |
| **(M4) 4D lift of the 1D compact-transfer gap** | the 3D-layer transfer operator gapped, uniformly | not compact in 4D; massless phases (§3.5) |

**Why this is genuinely the wall and not rhetoric.** Each face is a *volume-uniform* statement
*about the constant mode*, and each fails for a *different, structural* reason — FKG/order (M2),
strong-coupling radius (M3), non-compactness of the spatial slice (M4), massless reference (M1).
The RP machinery is exactly the union of (M1)–(M4)'s toolkits, and the failure is uniform across it
— this is the `[NG]` of §0. The missing object is not another inequality of a known type but the
**mass in the comparison**, which is the physics.

**What would actually close it (outside the RP class).** Per the campaign, a proof of `(M)` is one of
(a) a Savvidy-corrected uniform log-concavity of the coarse marginal on the typical bulk with a
tail-controlled essential infimum (`H6` §8, `(◆)=H1`); (b) a second invariant cone / total-positivity
(Karlin flag) on `𝒜/𝒢` below positivity (absent — no order); or (c) the two-scale Feshbach/Schur
quotient whose hypothesis `Σ≥μ` is the TCI = H1. None is an RP infrared bound.

---

## 5 · Acid tests `[T]+[V]`

**(A) Free massless (non-compact) U(1) photon, 4D — must NOT be bounded.** `χ = 1/m² → ∞` as the
IR mass `m→0`, every `L`; the constant gauge mode is a flat direction (pure gauge), `χ=∞` at `m=0`.
The free Maxwell measure is **reflection-positive and Gaussian**, so it meets every FSS hypothesis;
any finite `χ` from FSS would be a contradiction. FSS returns `χ ≤ Ĉ_free(0) = ∞` — **vacuous, as
required**. Witness §6·A. **Passes.**

**(B) Compact `β=0` SU(2) — must be finite (it is).** `χ = 3.719824` exactly, `L`-independent,
from pure compactness; the non-compact photon gives `∞` (§2·F3). Witness §6·B. **Passes
discriminatingly:** the finiteness is compactness-driven, not free-field-driven.

**(C) 1D SU(2) — must reproduce finite `χ` (it does).** Compact transfer operator ⇒ gap `Δ>0` ⇒
`χ_1D=(1+r)/(1−r)<∞` (`r=0.507`, `Δ=0.679`, `χ=3.06` at `β=2.5`). Witness §6·F. **Passes.**

**(D) 4D SU(N) at `β=2.5` — must fail to bound (it does).** Boundedness gives `O(L⁴)`; reweighting
gives `e^{Θ(L⁴)}`; Dobrushin `α=180`; FSS at `k=0` is `∞`. Witness §6·C–E. **Passes (by failing
where it should).**

The discriminator across (A)–(D) is exactly `compactness`: it makes `β=0` finite and 1D gapped
(both `[T]`), and it *would* make 4D finite *if* the interaction reweighting or the spatial
decorrelation were volume-uniform — which is H1. A bound that returns finite `χ` for the free
photon is **not** produced: the compactness is used, and the free photon (RP, `χ=∞`) correctly
escapes it.

---

## 6 · Witness — `research/v3/w7/q1_infrared_bound.py`

`py -3.12 q1_infrared_bound.py quick` (numpy 2.4.6, scipy 1.18.1; exit 0). Quoted (abridged):

```
[A] FREE NONCOMPACT massless d=4: chi = 1/m^2 -> inf, every L:
      m=1.0:1.00  0.5:4.00  0.2:25.0  0.1:100  0.05:400  0.02:2500  0.01:10000
    non-constant mean (1/#k!=0) sum 1/lambda_k: L=4:0.1476  L=6:0.1513  L=8:0.1528  (finite)
    => FSS bounds the non-constant average; the k=0 mode it cannot touch.           [T]
[B] COMPACT beta=0 SU(2): Var_Haar(xi3) = (4pi^2-6)/9 = 3.719824 (L-independent);
    G'' = 0.268830; quadrature cross-check matches to 3e-15.  Finiteness = compactness. [T]+[V]
[C] bounded range vs variance: independent Var ~ 3.27 (sigma^2, all N);
    perfectly correlated Var = 52, 210, 846, 3369 for N=16..1024 (= N sigma^2, grows);
    trivial bound pi^2 N grows.  Volume-uniform Var needs DECORRELATION = the gap.   [T]+[V]
[D] bounded-density: 0<=S_W<=2 beta N_p, inflation <= e^{2 beta N_p}, 2 beta N_p = 30 L^4:
    L=2:480  L=3:2430  L=4:7680  L=6:38880  L=8:122880  (e^{Theta(L^4)}).          [T]
[E] Dobrushin alpha(beta)=72 beta; alpha<1 iff beta<1/72; alpha(2.5)=180 (>>1).     [T] (S6)
[F] 1D COMPACT SU(2) transfer: r=fhat_{1/2}/fhat_0, gap=-log r, chi=(1+r)/(1-r):
    beta=1.0: r=0.240 Delta=1.426 chi=1.63; beta=2.5: r=0.507 Delta=0.679 chi=3.06. [T]+[V]
```

*Honest scope.* [A] is exact (free field, `χ=1/m²`); [B] is exact (closed form + quadrature);
[C] is a Monte-Carlo demonstration of the bounded-range-vs-variance distinction (exact in
expectation); [D],[E] are exact arithmetic in the campaign's normalization; [F] is exact quadrature
of the SU(2) Fourier coefficients (Peter–Weyl), giving the compact-transfer gap. It is a **witness,
not a proof** of the missing lemma; the theorem content is §1–§2 (FSS structure, the six failures),
and the open step is §4.

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| FSS Gaussian domination + RP = physical bound | `⟨Φ(f)²⟩_c ≤ (f,C_free f)` | **[std]** (FSS CMP 50) |
| FSS bound is on `(−Δ)⁻¹` ⟂ constants; `k=0` excluded | `χ ≤ Ĉ_free(0) = ∞` | **[T]** |
| free (non-compact) photon is RP with `χ=∞` | constant mode pure gauge, flat direction | **[T]**+**[V]** |
| `β=0` compact SU(2): `χ = (4π²−6)/9 = 3.7198` | pure compactness, `L`-independent | **[T]**+**[V]** |
| compact transfer (1D SU(2)) gap ⇒ `χ_1D<∞` | `r<1`, `Δ=0.68`, `χ=3.06` at `β=2.5` | **[T]**+**[V]** |
| compactness alone ⇒ only `Var ≤ π²N = O(L⁴)` | bounded range ≠ bounded variance | **[T]**+**[V]** |
| interaction reweighting `e^{2βN_p}=e^{Θ(L⁴)}` | bounded density not volume-uniform | **[T]**+**[V]** |
| Dobrushin `α(2.5)=180`, threshold `β<1/48` | disjoint from tower `[2.5,∞)` | **[T]** (S6) |
| no FKG/order on `𝒜/𝒢`; no `χ(k)≤χ(0)` | needed ferromagnetic comparison absent | **[T]**/**[O]** (P6 §2.1) |
| Savvidy indefinite Hessian kills Gaussian domination prereq | `λ_min=−g|B|`; double-well normal form | **[T-standard]** (S1/H6) |
| heat-kernel / fast-fibre gap is silent on `G″` | two-coordinate independence | **[T]** (H6 §2, S1 §6) |
| 1D compact-transfer gap does not lift to 4D | 3D-layer operator not compact; massless phases | **[T]** |
| conditional: `(M) ⇒ χ ≤ 3.72+C₀`, `λ=(3.72+C₀)⁻¹` | the shape an RP bound would yield | **[C]** |
| **no RP/IR bound reaches `β=2.5`** | the method class {FSS, chessboard, compact-T, heat-kernel, Dobrushin} | **[NG]** |
| **finite-mass domination at `k=0` (non-Abelian FKG/diamag.)** | `Var_μ(a) ≤ Var_Haar(a)+C₀`, uniform | **[O] = H1** |
| **`χ<∞` in 4D at `β=2.5` (H1); Clay** | — | **[O]** / **0%** |

---

## 8 · Sharpest line

**An RP infrared bound reaches `β=0` (compactness: `χ=3.7198` exactly) and the 1D chain (compact
transfer gap: `χ_1D=3.06` at `β=2.5`), but it does not reach 4D at `β=2.5`.** The FSS
Gaussian-domination inequality `⟨Φ(f)²⟩_c ≤ (f,C_free f)` is **vacuous at `k=0`** because the free
massless covariance has `Ĉ(0)=∞` — the bound controls the non-constant/gradient modes and the
constant mode is the zero mode of `(−Δ)`. The free (non-compact) photon is a **reflection-positive
object with `χ=∞`**, so the acid test forbids any RP+free argument from bounding it; the bound must
use compactness/non-linearity, and on SU(N) compactness gives only `Var(a) ≤ π²N_links = O(L⁴)`,
while the interaction reweighting is `e^{2βN_p}=e^{Θ(L⁴)}` and the Dobrushin/influence route is
strong-coupling (`α(2.5)=180`). The 1D compact-transfer theorem is a genuine compactness result but
does **not** lift to 4D. **The exact missing lemma is a volume-uniform finite-mass Gaussian
domination at the constant mode — a non-Abelian FKG / diamagnetic inequality
`Var_μ(a) ≤ Var_Haar(a)+C₀` (equivalently `G″_β(a₀*) ≥ (3.7198+C₀)⁻¹`)**; its four faces (massive
domination, FKG/order, uniform decorrelation, 4D compact-transfer) all reduce to H1. **The RP route
is closed at `β=2.5`; the wall is the mass in the Gaussian comparison, which is the gap.**
No QED. Clay 0%.

---

*Refs (repo).* `P-WAVE-SYNTHESIS.md` (periodic floor `0.222`, `β=0` `0.2688`, transfer `Δ=0.6444`);
`D6-DECISION.md` §1,§5 (dichotomy `χ<∞⟺saturation`; `G″=1/χ`; `(B-ineq)`); `H6-HEAD-SPECTRAL.md`
§1–§3,§8 (Riccati identity, BE essential-infimum, fast-fibre silence, double well, `(◆)=H1`);
`S1-CONDITIONAL-POINCARE.md` §1–§6 (SCBL `G″=S−D`, Savvidy `λ_min<0`, coarse `S` = H1);
`S6-DOBRUSHIN.md` §2,§6 (`α(β)=72β`, free-field `G″=0`, Haar `W_geo=0.205`);
`C2-BALABAN-UNIFORM.md` §1 (UV not IR; RP-1/RP-2; Lemma B); `P6-PERIODIC-TO-H1.md` §2.1 (no
`χ(k)≤χ(0)`), §1.4 (`f₁²=0`); `docs/paper-v2/ch02-squeeze.tex` (`eq:chi`, `eq:h1`, thin atom,
`§schur`); `docs/paper-v3/main.tex` §6 (`eq:fisher`, `eq:essinf`, spectral no-go, decompactification).
External (standard, not downloaded): Fröhlich–Simon–Spencer, CMP **50** (1976) 79–95;
Osterwalder–Seiler, Ann. Phys. **110** (1978) 440–471; Brydges–Fröhlich–Seiler, Ann. Phys. **121**
(1979) and CMP **91** (1983) (gauge fields on a lattice); Balaban, CMP **102/109/116/119/122**
(1985–89); Dobrushin (1970), Dobrushin–Shlosman (complete analyticity); Bakry–Émery, Brascamp–Lieb
(J. Funct. Anal. **22** (1976) 366); Prékopa (1971); Shen–Zhu–Zhu, CMP **400** (2023) 805,
arXiv:2204.12737. Witness: `research/v3/w7/q1_infrared_bound.py`. **No third-party sources
downloaded.**
