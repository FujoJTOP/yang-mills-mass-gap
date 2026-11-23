# Numerical falsification test of `CPB-W` (4D pure SU(3) lattice gauge theory)

> 2026-11-22 · directory `D:/Dev/math/research/lattice-cpb/`
> Python 3.12 + NumPy/SciPy, no GPU. All numbers below are produced by the
> scripts in this directory; literature values are marked as such.

---

## 0. Verdict (read this first)

**INCONCLUSIVE.  No support for `CPB-W` is found; it is not falsified either.**

On the **primary grid (β = 2.6–3.3, the physical handoff, `g² = 1.8–2.3`)**,
and at every other coupling probed:

* **`R₄` passes.**  On the primary grid `R₄ = −0.125…−0.131` at `b=1` and
  `−0.012…+0.003` at `b=2`; volume-independent; no non-Gaussian clustering at
  any β.  The "`R₄ ≫ 1`" failure mode is **absent**.
* **`κ_eff < 0` everywhere** (−2.82…−3.42 on the primary grid, 100 % of
  configurations, and **−1.74 at a β=40 weak-coupling control**).  The literal pass criterion fails, but the
  control shows the diagnostic does **not discriminate**, and the analytic
  side (`handoff-proof`, T4) has excluded it from the CPB decision.  It is
  reported for completeness, **not used**.
* **`ξ/a` is unresolvable at β ≤ 3.3** (`ξ/a ≲ 0.3`: the connected plaquette
  correlation dies inside one lattice spacing, before and after APE
  smearing), and `0.93(8)` at β=6.  The σ-certificate target `ξ/a ≈ 7` is
  **out of reach at the measured 7× kernel** but not absurd with 16 chains
  and a faster kernel: `L ≈ 21–28`; with 16 independent chains,
  **56–222 h at 5 %** (2–9 days) or **3.5–14 h at 20 %** with a 50× kernel.
  The catch: 20 % precision only certifies the σ cap if the measured `ξ/a`
  is ≳11.7, and a 7× single-core run is 400–1 600 h.  See §10.5–§10.6.
* The **decisive object `z(q)`/`κ_c` is not measurable here**: the polymer
  weights require the blocked effective action (inverse MC or an analytic
  character expansion).  The closest proxies — all plaquette character
  coefficients positive, `r₁ = 1.0000` — show no `c₁`-level sign cancellation,
  consistent with the analytic default `κ_c ≈ 1`.

**Bottom line:** what can be measured points to "nothing special" (small `R₄`,
same-sign character expansion) and gives no evidence for the overlap the
handoff needs; but the object that would decide it — `z(q)` — was not
constructed, and the lattice volumes needed for the σ certificate are far out
of reach.  Details, definitions, error treatment and caveats below; §9
contains the Bianchi/`z(q)` analysis and the cost/feasibility table.

## 1. What is being tested

The handoff inequality `H3: q₀ ≥ q₁` of `research/handoff-proof` needs the
estimate **`CPB-W`**: the weak-side (massive background-field) configuration
class is preserved up to `λ_c ≈ 6`–`12`, with relative error ≤ 1/2, where
`λ* = h^∨ g² = 3 g²` for SU(3).

| β | g² = 6/β | λ* = 3g² | t* = 24π²/(11λ*) | e^{−t*} |
|---|---|---|---|---|
| 1.7 | 3.529 | 10.588 | 2.034 | 1.31×10⁻¹ |
| 2.2 | 2.727 | 8.182 | 2.632 | 7.19×10⁻² |
| 5.7 | 1.053 | 3.158 | 6.819 | 1.09×10⁻³ |
| 6.0 | 1.000 | 3.000 | 7.178 | 7.63×10⁻⁴ |

**Which β is what (updated 2026-11-22, per `handoff-proof`):**

* **Primary (physical handoff): β ≈ 2.6–3.3**, i.e. `g² = 1.8–2.3`,
  `λ* = 3g² = 5.4–6.9`.  This is where the physical handoff `μ ≈ Δ`
  corresponds to `q ~ 1`; the runs at β = 2.6, 2.8, 3.0, 3.3 are the
  **primary grid**.
* **Conservative strong-coupling checks: β = 1.7, 2.2** (`λ* = 8.2–10.6`).
  These are inside the *one-loop-conservative* reach `λ_c ≥ 12` and are harder
  than the physical point.
* **Formal weak boundary only: β = 5.7–6.0** (`g² ≈ 1`, `λ* ≈ 3`).  This is
  the `g² = 1` boundary of the formal weak-side statement, **not** the
  physical handoff — at β = 6 the lattice spacing is ~0.1 fm and the
  continuum scaling regime is approached, whereas at β = 1.7–3.3 the lattice
  is coarse (`a ≈ 0.2–0.6 fm`) and the strong-coupling description applies.

**`t_up`: one loop is not enough.**  With the SU(3) MS-bar coefficients
(`b₀ = 11/(16π²)`, `b₁ = 102/(16π²)²`, `b₂ = 2857/2/(16π²)³`) and
`t_up(g²) = ∫_{g²}^∞ dx / [2x²(b₀ + b₁x + b₂x²)]`:

| β | g² | λ* | t_up 1-loop | 2-loop | 3-loop | q₀ (m/Λ=8), 1/2/3-loop |
|---|---|---|---|---|---|---|
| 3.3 | 1.818 | 5.455 | 3.95 | 2.96 | 2.67 | 0.154 / 0.414 / 0.552 |
| 3.0 | 2.000 | 6.000 | **3.59** | **2.639** | 2.356 | **0.221 / 0.571 / 0.758** |
| 2.8 | 2.143 | 6.429 | 3.35 | 2.43 | 2.15 | 0.281 / 0.707 / 0.935 |
| 2.6 | 2.308 | 6.923 | 3.11 | 2.21 | 1.94 | 0.357 / 0.874 / 1.150 |
| 2.2 | 2.727 | 8.182 | 2.63 | 1.80 | 1.53 | 0.576 / 1.326 / 1.725 |
| 1.7 | 3.529 | 10.588 | 2.03 | 1.29 | 1.05 | 1.047 / 2.200 / 2.804 |
| 6.0 | 1.000 | 3.000 | 7.18 | 5.96 | 5.64 | 0.006 / 0.021 / 0.028 |

(My 1-loop and 2-loop columns reproduce the values quoted by `handoff-proof`
exactly: `t₁(g²=2)=3.589`, `t₂(g²=2)=2.639` (the rigorous RB2 bound).  The
3-loop column differs from the `numerics/REPORT.md` "~2.0" model-based value
at the ~15 % level — that number is model/scheme dependent; the qualitative
conclusion is unchanged.)

**Consequence.**  At the physical handoff `g² ≈ 2` a **one-loop control only
reaches `q₁ ≤ 0.22`** (`q₀ = 0.22` for `m/Λ = 8`), while the 2-loop (RB2)
bound gives `q₀ ≈ 0.57` and 3-loop `≈ 0.76`.  Therefore `CPB-W` must deliver
**at least a 2-loop bound on `t_up`**, not merely class preservation: the
one-loop statement alone leaves `q₁ ≤ 0.22`.

The proposed falsification criteria (from `research/handoff-proof/STATUS.md`):
* pass if the `z`-proxy is small or `κ_eff > 0`;
* fail if `R₄ ≫ 1` or `κ_eff < 0`.

No `z(q)` was constructed (the signed-gas polymer weights for the equivariant
block spin are not defined in this project yet), so `R₄` and the defect
density `ρ` stand in as qualitative proxies.

---

## 2. Lattice setup and update algorithm

* **Geometry**: `L⁴` torus, periodic boundary conditions, `L = 4, 6, 8` (and
  `L=10` at β=6.0 for a finite-temperature check).
* **Action**: Wilson, `S = β Σ_p (1 − ⅓ Re tr U_p)`, `β = 6/g²`, standard
  single-plaquette normalisation.
* **Update**: **Cabibbo–Marinari heat bath** — an *exact* SU(2) heat bath on
  each of the three SU(2) subgroups (0,1), (0,2), (1,2), applied sequentially,
  with a checkerboard decomposition in the 4D parity and over the four
  directions (8 Gibbs phases per sweep).
  * The SU(2) sampler is exact for the density
    `∝ √(1−x²) e^{λx}`: a Beta-type rejection proposal for `λ ≤ 2` and a
    `Gamma(3/2,1/λ)` proposal with `√((2−z)/2)` acceptance for `λ > 2`
    (both exact; verified against numerical quadrature, see §2.2).
  * A heat bath has **no accept/reject step** — the "acceptance rate" is 100%
    by construction.  The meaningful efficiency is the SU(2) *proposal*
    efficiency, which is **0.31 (β=1.7), 0.45 (β=2.2), 0.97 (β=5.7/6.0)**.
  * Sweep cost on this machine: 35 ms (L=6), 120 ms (L=8), single core.

### 2.1 Validation (all in `validate.py`, all pass)

| check | result |
|---|---|
| staple sum vs. brute-force action derivative | rel. error 3×10⁻¹⁰ |
| `W(1,1)` vs. average plaquette | identical to 10⁻¹² |
| gauge invariance of `⟨P⟩`, `W(2,2)` | 10⁻¹⁷ |
| SU(2) sampler vs. quadrature, λ=0.3,2,5,12 | ≤3σ on ⟨x⟩, ⟨x²⟩ |
| heat bath vs. independent Metropolis (stationarity, L=4) | agree at β=1.0 (8×10⁻⁴) and β=2.2 (8×10⁻⁵) |
| heat bath at β=1.0 vs. leading strong-coupling β/18=0.0556 | 0.0602 ± 0.0003 (O(β³) correction) |
| three independent chains, L=4, β=2.2 | 0.14338, 0.14374, 0.14355 |
| β=6.0, L=8 average plaquette | 0.5943 (literature value ≈ 0.5937) |

### 2.2 Parameters, thermalisation, autocorrelation

| β | L | sweeps | thermal. | measurements | ⟨P⟩ | τ_int(P) | bin |
|---|---|---|---|---|---|---|---|
| 1.7 | 4 | 20000 | 1000 | 10000 | 0.107453(68) | 0.5 | 1 |
| 1.7 | 6 | 16000 | 1500 | 8000 | 0.107425(33) | 0.5 | 1 |
| 1.7 | 8 | 6000 | 800 | 3000 | 0.107394(31) | 0.5 | 1 |
| 2.2 | 4 | 20000 | 1500 | 10000 | 0.143431(70) | 0.5 | 1 |
| 2.2 | 6 | 16000 | 1500 | 16000 (2 chains) | 0.143465(24) | 0.54 | 2 |
| 2.2 | 8 | 6000 | 800 | 3000 | 0.143505(31) | 0.51 | 2 |
| 2.6 | 6 | 16000 | 1500 | 8000 | 0.173670(40) | 0.52 | 2 |
| 2.8 | 6 | 16000 | 1500 | 8000 | 0.189210(30) | 0.50 | 1 |
| 3.0 | 6 | 16000 | 1500 | 8000 | 0.205050(40) | 0.50 | 1 |
| 3.3 | 6 | 16000 | 1500 | 8000 | 0.229440(40) | 0.50 | 1 |
| 5.7 | 6 | 6000 | 600 | 6000 | 0.549310(254) | 16.6 | 34 |
| 6.0 | 6 | 6000 | 600 | 6000 | 0.594803(69) | 1.9 | 4 |
| 6.0 | 8 | 4000 | 600 | 4000 | 0.594270(50) | 2.4 | 5 |
| 6.0 | 10 | 1000 | 300 | 1000 | 0.59401(9) | 5.0 | 10 |

All chains are **cold-start**; equilibration was checked by the plaquette trace
and by two independent chains at β=2.2 (which agree to 3×10⁻⁴, i.e. within
~1σ of the binned errors).  **Autocorrelation treatment**: `τ_int` is computed
from the measured plaquette series with a Sokal window; the error bars quoted
everywhere are **jackknife over bins of size `B = max(1, ⌈2τ_int⌉)`** (last
column).  All observables — including the correlators, `R₄` and the defect
fractions — use this common binning.  For β=5.7 the autocorrelation is long
(τ≈17 sweeps), so only ~176 effectively independent bins exist; this is
reflected in its larger error bars.

**Finite-temperature caveat.**  At `L=6` the temporal extent is `N_t=6`, and
the SU(3) deconfinement transition is at `β_c(N_t=6) ≈ 5.89`
(`β_c(N_t=8) ≈ 6.06`; literature).  Hence **β=6.0 on L=6 is on the
deconfined side**, while β=5.7 (L=6) and β=6.0 (L=8) are on the confined side.
The L=10 β=6.0 run was added as a cleaner zero-temperature reference.
This is directly visible in the measured Polyakov loop:
`|⟨P⟩| = 0.0319(5)` at L=6, β=6.0 (deconfined), versus
`0.0016(4)` at L=8, β=6.0 and `0.0020(3)` at L=6, β=5.7 (both consistent with
zero — confined); at β=1.7–2.2 all lattices give `|⟨P⟩| < 4×10⁻⁴`.
This caveat matters only for the β≈6 reference point, not for the β=1.7/2.2
`CPB-W` probe.

---

## 3. Observables: exact definitions

### 3.1 Average plaquette and Wilson loops

`⟨P⟩ = (1/6L⁴) Σ_{μ<ν,x} ⅓ Re tr U_{μν}(x)`.
`W(R,T) = (1/6L⁴) Σ_{μ<ν} Σ_x ⅓ Re tr` of the rectangular `R×T` loop, all
orientations, all base points; periodic images are included (for `R,T ≤ L/2`
this is the standard "all loops that fit" average).  Extracted for
`R,T ≤ min(3, L/2)`.

### 3.2 `ξ/a`

Operator: the zero-momentum **magnetic (space-space) plaquette** density
`s(x) = ⅓ Σ_{1≤i<j≤3} ⅓ Re tr U_{ij}(x)`, summed over each time slice:
`S(t) = Σ_x s(x,t)`.  Connected correlator
`C(t) = ⟨S(t₀)S(t₀+t)⟩ − ⟨S⟩²` (periodic in t).
Fit `C(t) = A[e^{−mt} + e^{−m(L−t)}]` on `t = 1 … ⌊L/2⌋` (stated window:
t=1,2,3 for L=6; t=1,2 for L=4), `ξ/a = 1/(am)`; errors by jackknife.
The effective mass `am_eff(t) = arccosh[(C(t−1)+C(t+1))/(2C(t))]` is also
reported.

**The unsmeared plaquette operator has almost no overlap with the lightest
state at β=1.7–2.2** (the connected correlation beyond one time step is
consistent with zero), so a dedicated run with **APE-smeared** operators
(α=0.5, n=0 and n=8 smearing steps on the spatial links, temporal links
untouched; `corr_run.py`, `corr_analyze.py`) was performed.  The smeared
results are quoted; the unsmeared fits are reported only to show they are
unreliable.  The cut-off on `ξ/a` from the lattice is `ξ/a ≲ L/2 = 3`.

### 3.3 `R₄`

Plaquette field `P(x) = (1/6) Σ_{μ<ν} ⅓ Re tr U_{μν}(x)`.
**Blocking**: the `L⁴` torus is partitioned into `(L/b)⁴` geometric cubes of
side `b` (`b | L`, `b ∈ {1,2,3,4}`); the blocked variable is the *unweighted
mean over the `b⁴` sites of the block*,
`O_b(y) = b^{−4} Σ_{x∈block y} P(x)`.
Then
`R₄(b) = κ₄/κ₂²`, with `κ₂ = Var(O_b)`, `κ₄ = ⟨(O_b−⟨O_b⟩)⁴⟩ − 3 Var(O_b)²`
(the fourth cumulant = connected 4-point function).  Moments are pooled over
all blocks and all configurations; errors by jackknife over configuration
bins.  `R₄ = 0` for a Gaussian field; `R₄ ≫ 1` would indicate strong
non-Gaussian clustering (a two-state / droplet-like structure).

### 3.4 Defect density `ρ`

Operational definition: `ρ(θ) =` fraction of the `6L⁴` plaquettes with
`⅓ Re tr U_p < θ`, for `θ ∈ {0, 0.1, 0.2, ⅓}`.  `ρ(0)` counts
"wrong-sign" plaquettes and is the primary number.  In addition we report a
**curvature defect density**: the fraction of links whose local 8×8 action
Hessian (left-algebra basis, all other links fixed) has a negative eigenvalue;
this is gauge covariant (up to an adjoint rotation) and is the local analogue
of `κ_eff`.

### 3.5 `κ_eff` — lowest quotient-Hessian eigenvalue

*Coordinates*: for each link, `U(x,μ) = exp(i ξ_a T^a) U0(x,μ)`, with `T^a`
the Gell-Mann generators; the natural metric is Euclidean in `ξ`.
*Hessian*: the **bi-invariant** second derivative of `S_W`, obtained from
central differences (`ε=2×10⁻³`) of the exact gradient
`∂S/∂ξ_a = (β/3) Im tr(T^a U S)` plus the analytic Baker–Campbell–Hausdorff
correction `−½ f^{abc} G_c` that converts the coordinate Jacobian into the
symmetric Hessian.  Verified: `⟨v,Hv⟩` equals the second derivative of the
action along the geodesic `e^{iεv}U0` to ~10⁻⁷ relative; the operator is
symmetric to 10⁻¹⁰.
*Gauge-invariant subspace*: the tangent space to the gauge orbit is the range
of the covariant derivative `(Dω)(x,μ) = ω(x) − Ad_{U(x,μ)} ω(x+μ)`; `H`
annihilates it to ≤2×10⁻² (measured: |⟨Dω,HDω⟩|/|Dω|² ≤ 0.02 while |HDω|/|Dω| ≈ 1.6, i.e. H maps the orbit into the transverse complement, as a non-stationary configuration requires), and the orthogonal projector
`P = D(D†D+μ₀)⁻¹D†` (CG, `μ₀=10⁻³`) is applied so that
`κ_eff = min eig( (1−P) H (1−P) )` on `range(1−P)`.  `κ_eff` is computed by
Lanczos with full reorthogonalisation (`m=120`, residual ≤10⁻⁵ for L=6).
*Validation*: on L=2 the dense projected Hessian reproduces Lanczos to 6
digits and has exactly `8L⁴ = 128` eigenvalues at the lifted orbit value
`Λ=10³`; on the ordered configuration and on pure-gauge configurations
`κ_eff = 0` exactly (flat-connection zero modes), as it must.
*Limitations (important)*:
1. `P` projects onto the **linearised** gauge orbit.  The orbit–orbit block of
   `H` vanishes (measured: `|⟨Dω,HDω⟩|/|Dω|² ≤ 2×10⁻²`), but the orbit is not
   in general a totally geodesic submanifold of the configuration space, so the
   compression equals the Riemannian Hessian on the quotient `M/G` only up to
   a term `⟨∇S, II(v,v)⟩` involving the second fundamental form `II` of the
   orbit.  That term vanishes at stationary configurations and is *not*
   computed here.  This is a genuine uncertainty in the interpretation of
   `κ_eff`, on top of limitation 2.
2. This is the Hessian of the **microscopic Wilson action at a typical
   equilibrium configuration**, *not* the Hessian of a Wilsonian effective
   action for block variables.  At a non-stationary configuration the two are
   different objects.  This is the single most important caveat for the
   `CPB-S` reading of `κ_eff`.
3. `ε` and `μ₀` are finite; the results are stable under `ε∈[10⁻⁴,10⁻²]`,
   `μ₀∈[10⁻⁶,10⁻²]`.

---

## 4. Results

### 4.1 Average plaquette

| β | L=4 | L=6 | L=8 | L=10 |
|---|---|---|---|---|
| 1.7 | 0.107453(68) | 0.107425(33) | 0.107394(31) | — |
| 2.2 | 0.143431(70) | 0.143465(24) | 0.143505(31) | — |
| 2.6 | — | 0.173670(40) | — | — |
| 2.8 | — | 0.189210(30) | — | — |
| 3.0 | — | 0.205050(40) | — | — |
| 3.3 | — | 0.229440(40) | — | — |
| 5.7 | — | 0.549310(254) | — | — |
| 6.0 | — | 0.594803(69) | 0.594270(50) | 0.59401(9) |

Finite-size dependence is at the 10⁻⁴ level for β=1.7/2.2 and decreases with
increasing L at β=6.0 (0.59480 → 0.59427 → 0.59401), consistent with the
literature value ≈0.5937 for the infinite-volume β=6.0 plaquette; the L=10
run (N_t=10, effectively zero temperature) is the cleanest reference.  The leading
strong-coupling value `β/18` is 0.0944 and 0.1222 at β=1.7, 2.2; the measured
values are 14–17% above it, as expected from the first O(β³) correction.

### 4.2 Wilson loops (L=6)

| β | W(1,1) | W(1,2) | W(1,3) | W(2,2) | W(2,3) | W(3,3) |
|---|---|---|---|---|---|---|
| 1.7 | 0.10746(10) | 0.01143(9) | 0.00137(10) | 0.00019(10) | −0.00006(9) | 0.00003(10) |
| 2.2 | 0.14334(8) | 0.02056(7) | 0.00298(7) | 0.00041(7) | −0.00007(7) | 0.00005(7) |
| 5.7 | 0.54930(41) | 0.32490(60) | 0.19550(60) | 0.13254(63) | 0.05789(54) | 0.01978(35) |
| 6.0 | 0.59475(17) | 0.38555(25) | 0.25489(29) | 0.19318(30) | 0.10453(27) | 0.05021(23) |

(correlated errors at the 10⁻⁴ level; loops are measured every `loop_every`
sweeps, hence the larger errors compared with the plaquette, which is
`W(1,1)` measured every sweep to 3×10⁻⁵.)

The strong-coupling loops obey a steep area law; at β≥5.7 the perimeter term
dominates.  Beyond `R=T=2` at β=1.7–2.2 the signal is at the noise level — a
statement about the size of the lattice, not about the theory.

Creutz ratios `χ(R,T) = −ln[W(R,T)W(R−1,T−1)/(W(R−1,T)W(R,T−1))]`
(jackknife errors, same bins):

| β | L | χ(1,1) | χ(2,2) |
|---|---|---|---|
| 1.7 | 4 | 2.2294(20) | 1.905(1063) |
| 1.7 | 6 | 2.2307(9) | 1.887(525) |
| 1.7 | 8 | 2.2305(9) | 2.199(620) |
| 2.2 | 4 | 1.9436(16) | 2.151(556) |
| 2.2 | 6 | 1.9425(5) | 1.966(163) |
| 2.2 | 8 | 1.9413(7) | 1.801(190) |
| 5.7 | 6 | 0.5991(7) | 0.3714(20) |
| 6.0 | 6 | 0.5196(3) | 0.2574(7) |
| 6.0 | 8 | 0.5205(2) | 0.2619(6) |

At these loop sizes (`R=T≤2`) the Creutz ratio is still dominated by the
short-distance (perturbative) part of the potential; it is an *upper* proxy
for `σa²`, so `ξ_σ/a = 1/√(σa²) ≳ 0.7` at β=1.7–2.2 and `≳2` at β=6.  These
are *not* the mass gap; they only fix the order of magnitude and confirm that
β=1.7–2.2 is far from the scaling regime.

### 4.3 `ξ/a` (smeared connected correlator)

Operator and fit: `S(t) = Σ_x s(x,t)` with `s` the magnetic plaquette density
of §3.2, connected correlator, periodic cosh fit on `t=1…⌊L/2⌋`.  The operator
is used both unsmeared (n=0) and after 8 APE steps on the spatial links
(α=0.5, temporal links untouched), from dedicated 4000-sweep runs
(`corr_run.py`).  Errors: jackknife over config bins.

| β | L | n_smear | C(1)/C(0) | C(2)/C(0) | fit `am` (window) | `ξ/a` |
|---|---|---|---|---|---|---|
| 1.7 | 4 | 0 / 8 | −0.0067(78) / 0.0018(78) | −0.010(11) / −0.026(11) | — (unresolvable) | — |
| 1.7 | 6 | 0 / 8 | −0.0062(65) / 0.0002(64) | 0.005(7) / −0.003(6) | — (unresolvable) | — |
| 1.7 | 8 | 0 / 8 | −0.0027(78) / −0.0076(78) | 0.006(8) / 0.013(8) | — (unresolvable) | — |
| 2.2 | 4 | 0 / 8 | 0.0074(81) / 0.0027(79) | 0.002(11) / 0.010(11) | — (unresolvable) | — |
| 2.2 | 6 | 0 / 8 | 0.0028(64) / −0.0044(65) | −0.002(6) / −0.001(6) | — (unresolvable) | — |
| 2.6 | 6 | 0 / 8 | 0.0087(63) / −0.0016(60) | 0.0083(66) / 0.0011(63) | — (unresolvable) | — |
| 3.0 | 6 | 0 / 8 | 0.0169(61) / −0.0046(65) | 0.0049(65) / −0.0054(68) | — (unresolvable) | — |
| 3.3 | 6 | 0 / 8 | −0.0039(65) / 0.0067(60) | −0.0001(67) / −0.0005(60) | — (unresolvable) | — |
| 6.0 | 6 | 0 | 0.0962(73) | 0.0172(71) | 1.68(41) (t=1–3) | 0.60(14) |
| 6.0 | 6 | 8 | 0.2689(119) | 0.1062(119) | **1.073(88)** (t=1–3) | **0.93(8)** |

`am_eff(t)` for the β=6, n=8 case is 1.35, 1.00, 1.22 at t=1,2,3 — roughly
consistent with a single exponential at `am ≈ 1.0–1.3`; the unsmeared
operator gives `am_eff(1)=2.35`, i.e. strong excited-state contamination, which
is what the smearing is for.

**At β = 1.7–3.3 — including the whole primary grid — the connected
plaquette correlation is not resolvable**, on L=4, 6 and 8, before and after
smearing: `C(1)/C(0) = −0.0016(60) / −0.0046(65) / 0.0067(60)` at
β = 2.6/3.0/3.3 (and similar at 1.7/2.2).  Taking the 2σ bound
`|C(1)/C(0)| < 0.013` with a single-exponential model gives `am > 4.4`,
`ξ/a < 0.23`; but because the same null result is obtained with and without
smearing — which changes the overlap by an order of magnitude at β=6 — we
read this as an **operator-overlap limitation, not a physical determination**.
The honest statement is `ξ/a ≲ 0.3 (not resolved)` at β ≤ 3.3, consistent with
a mass gap of `am ~ 3–5` in lattice units at these coarse spacings.
This is the direct numerical reason why the σ-certificate target
`ξ/a ≈ 7` cannot be reached on these lattices (§9.6).

For the β≈6 reference we quote `ξ/a = 0.93(8)` (smeared, L=6); the unsmeared
fits give 0.60(14) at L=6 and 1.22 at L=10, so the systematic spread from the
operator choice and finite volume is ~±0.3 — i.e. `ξ/a ≈ 1` at β=6, an order
of magnitude larger than at the strong-coupling probe.


### 4.4 `R₄` (connected 4-point ratio of the blocked plaquette)

| β | L | b=1 | b=2 | b=3/4 |
|---|---|---|---|---|
| 1.7 | 6 | −0.0952(14) | −0.0103(59) | −0.0084(134) |
| 2.2 | 6 | −0.1106(9) | −0.0067(43) | +0.0041(97) |
| 2.6 | 6 | −0.1247(13) | −0.0119(57) | −0.0035(139) |
| 2.8 | 6 | −0.1246(12) | −0.0042(64) | −0.0064(134) |
| 3.0 | 6 | −0.1284(13) | −0.0086(58) | −0.0209(131) |
| 3.3 | 6 | −0.1308(13) | +0.0027(65) | −0.0280(131) |
| 5.7 | 6 | +0.1730(27) | +0.0408(82) | +0.0352(213) |
| 6.0 | 6 | +0.2017(29) | +0.0292(78) | +0.0031(153) |
| 1.7 | 8 | −0.0951(13) | −0.0062(56) | −0.0538(216) |
| 2.2 | 8 | −0.1139(12) | −0.0048(55) | −0.0060(220) |
| 6.0 | 8 | +0.2015(20) | +0.0318(51) | −0.0242(187) |
| 6.0 | 10 | +0.2055(22) | +0.0382(79) | — |

`R₄` is **small at every blocking level**, of the order expected from the
central-limit suppression by the number of effective degrees of freedom in a
block, and it is **volume independent** (L=6 vs L=8 agree).  Its sign changes
between the strong-coupling regime (negative, platykurtic — the plaquette is a
bounded variable) and β≈6 (positive, slightly leptokurtic), but in **no case is
`R₄ ≫ 1`**.  There is no trace of a two-state / droplet-like non-Gaussian
structure in the blocked plaquette at any of the couplings probed.

### 4.5 Defect density

| β | ρ(0) | ρ(0.1) | ρ(0.2) | ρ(⅓) | curvature defects (links with λ_min<0) |
|---|---|---|---|---|---|
| 1.7 | 0.38919(6) | 0.52238(7) | 0.64642(6) | 0.78716(5) | 0.68–0.71 |
| 2.2 | 0.33845(4) | 0.46668(5) | 0.59223(5) | 0.74321(4) | 0.55–0.59 |
| 2.6 | 0.29876(6) | 0.42150(7) | 0.54657(7) | 0.70422(6) | 0.47 |
| 2.8 | 0.27931(5) | 0.39873(6) | 0.52304(6) | 0.68346(5) | 0.41 |
| 3.0 | 0.26011(5) | 0.37604(6) | 0.49911(6) | 0.66184(5) | 0.36 |
| 3.3 | 0.23207(5) | 0.34198(5) | 0.46247(6) | 0.62776(5) | 0.30 |
| 5.7 | 0.01399(7) | 0.03197(9) | 0.06605(10) | 0.15303(12) | 0.00 |
| 6.0 | 0.00537(1) | 0.01460(2) | 0.03511(3) | 0.09711(4) | 0.00 |

The L=10, β=6.0 values are `[0.0055, 0.01484, 0.03554, 0.09802]`, identical to
L=8 within errors.  The "wrong-sign plaquette" density falls by ~70× between
β=2.2 and β=6.0, and
the local curvature defects disappear entirely above β≈5.7 (the local link
Hessian is positive definite: mean +3.05 at β=6.0, +17.5 at β=20).

### 4.6 `κ_eff`

### 4.6 `κ_eff`

Per-configuration values on saved snapshots (Lanczos `m=120` for L=6, `m=80`
for L=4; min eigenvalue, mean ± spread over snapshots):

| L | β | κ_eff | frac. configs κ<0 | max Lanczos residual | local 8×8 link Hessian: mean λ_min | frac. links λ_min<0 |
|---|---|---|---|---|---|---|
| 4 | 1.7 | −1.917(25) | 1.00 | 9.4×10⁻³ | −0.110 | 0.690 |
| 4 | 2.2 | −2.426(53) | 1.00 | 2.9×10⁻³ | −0.038 | 0.551 |
| 6 | 1.7 | −1.942(16) | 1.00 | 6.9×10⁻⁴ | −0.104 | 0.684 |
| 6 | 2.2 | −2.469(44) | 1.00 | 3.5×10⁻³ | −0.048 | 0.571 |
| 6 | 5.7 | −3.988(227) | 1.00 | 3.4×10⁻⁴ | +2.582 | 0.0007 |
| 6 | 6.0 | −3.500(161) | 1.00 | 5.7×10⁻⁵ | +3.054 | 0.0001 |
| 6 | 2.6 | −2.823(39) | 1.00 | 2.8×10⁻³ | +0.028 | 0.471 |
| 6 | 2.8 | −3.011(63) | 1.00 | 1.1×10⁻² | +0.085 | 0.413 |
| 6 | 3.0 | −3.190(33) | 1.00 | 3.5×10⁻³ | +0.146 | 0.363 |
| 6 | 3.3 | −3.424(64) | 1.00 | 4.6×10⁻⁴ | +0.248 | 0.296 |
| 6 | 8.0* | −2.208 | 1 | 8.6×10⁻³ | +5.218 | 0 |
| 6 | 10.0* | −1.999 | 1 | 2.8×10⁻³ | +7.343 | 0 |
| 6 | 20.0* | −1.729 | 1 | 1.5×10⁻³ | +17.482 | 0 |
| 6 | 40.0* | −1.740 | 1 | 1.2×10⁻² | +37.472 | 0 |

\* weak-coupling control, one configuration each (600 cold-start sweeps),
`m=160`; these are not production points.

**`κ_eff` is negative at every coupling probed, from β=1.7 to β=40**, with
`100%` of configurations negative and well-converged Lanczos residuals.  It is
**volume independent** (L=4 versus L=6 agree within 2%, despite the different
lowest momenta), so it is *not* a finite-volume collective mode.

The magnitude tracks the typical field strength:
`|κ_eff| ≈ (0.85–1.6)·β(1−⟨P⟩)`, and `β(1−⟨P⟩) ≈ 1.5–2.6` over the whole
range.  The local link curvature tells the complementary story: at β=1.7–3.3
30–70% of links have a negative local curvature (69 % at β=1.7, 47 % at
β=2.6, 30 % at β=3.3), while at β≥5.7 the local Hessian is positive definite
almost everywhere — yet the *global* quotient
Hessian remains negative (the negative mode is collective, a
Savvidy/Nielsen–Olesen-like mode, not a local one).

Why this does **not** falsify `CPB-W`: at β=40 the microscopic theory is in
the scaling region and is known to be healthy, and `κ_eff` is still −1.74.  A
diagnostic that assigns the same verdict to β=1.7 and to β=40 cannot
discriminate the crossover from the scaling regime.  Moreover the object
computed is the Hessian of the microscopic Wilson action, not of a Wilsonian
effective action (see §3.5, limitations).


### 4.7 Finite-size checks

`⟨P⟩`, `ρ(θ)` and `R₄(b=1,2)` agree between L=6 and L=8 to within their
errors at every β (see tables above); `κ_eff` agrees between L=4 and L=6 to
~2% (e.g. −1.94 vs −1.94 at β=1.7; −2.42 vs −2.48 at β=2.2).  The
`CPB-W` probe quantities are therefore volume-stable on the lattices used.

---

## 5. Plots

* `plots/plaquette_vs_beta.png` — ⟨P⟩ vs β, all L, with the strong-coupling
  line β/18.
* `plots/wilson_loops.png` — −log W(R,T) and the crude static potential.
* `plots/smeared_correlators.png` — connected correlators (smeared) and
  `am_eff`.
* `plots/R4_blocking.png` — R₄ vs blocking size b.
* `plots/defects.png` — defect fraction vs threshold.
* `plots/kappa_eff.png` — κ_eff and the local-curvature defect fraction vs β.
* `plots/thermalization.png` — cold-start thermalisation traces.
* `plots/smeared_correlators.png` — APE-smeared connected correlators and
  `am_eff(t)` (see §4.3).

---

## 6. Verdict and interpretation

### 6.1 By the project's stated criteria (primary grid included)

| criterion | measured on β = 2.6–3.3 (L=6) | literal verdict |
|---|---|---|
| `R₄ ≫ 1` ⇒ fail | `R₄(b=1) = −0.125…−0.131`, `R₄(b=2) = −0.012…+0.003` | **passes** |
| `κ_eff > 0` ⇒ pass | `κ_eff = −2.82…−3.42` (100 % of configurations) | **fails** |

So the project's literal decision rule returns **FAIL** (`κ_eff < 0`), and by
that rule the CPB route must be reconsidered.

### 6.2 Why that "fail" is a failure of the test, not of `CPB-W`

The same `κ_eff` measured at **β=40** — deep in the scaling region, where pure
SU(3) Yang–Mills is known to be a healthy confining theory — is **−1.74**,
i.e. negative with a similar magnitude.  A diagnostic that returns the same
sign at β=1.7, β=3.0 and β=40 cannot separate the crossover from the scaling
regime.  Its scaling `|κ_eff| ≈ (0.85–1.6)·β(1−⟨P⟩)` identifies it as a
Savvidy/Nielsen–Olesen-type negative mode of the **microscopic** Wilson
action, i.e. a property of any non-convex action, not of the effective action
that `CPB-S` refers to.  This is confirmed independently by `handoff-proof`
(T4), which excludes `κ_eff` from the CPB decision.

The other two measurements point the same way:
* `R₄` is small and volume-independent — there is no two-state/droplet
  structure for `CPB-W` to latch onto, but also no failure mode to trigger.
* `ξ/a` is below the lattice resolution at β ≤ 3.3, so the *strong-side*
  σ-certificate target (`ξ/a ≈ 7`) is two orders of magnitude away on these
  lattices; this is a volume problem, not an observable failure.

### 6.3 Scientific verdict

**INCONCLUSIVE.**  The numerical test as specified does not support `CPB-W`
(the stability proxy fails), but it does not falsify it either: the proxy
fails for trivial reasons in a control regime, the decisive object `z(q)` was
not constructed, and the volumes needed for the σ certificate are out of
reach.  For the project's decision this means: **do not continue on the basis
of this test; it supplies no positive evidence**, and a genuine test needs
either a measured `z(q)` (IMC-fitted blocked action / character expansion) or
a ≥20× faster kernel with L ≈ 28.

## 7. What would change the verdict

* A `z(q)` computed from the equivariant block-spin polymer weights — the
  object the criterion was designed for.
* An effective (blocked) action Hessian, obtained e.g. by inverse
  Monte-Carlo / Schwinger–Dyson, with `κ_eff > 0` in the slow sector.
* Evidence of `R₄ ≫ 1` in a larger-`b` / larger-volume study (a droplet
  structure would show up as a growing `R₄` with `b`); we see the opposite
  (`R₄` shrinks with `b`).
* A `CPB-W` proof at `λ_c ≈ 12` would of course supersede this test.

## 8. Files

```
su3.py            lattice, CM heat bath, exact SU(2) sampler, observables
su3_numba.py      Numba heat-bath kernel (7.0x measured) + cost projection §10
validate.py       all correctness checks (run this first)
run.py            production driver (--L --beta --sweeps --therm --seed ...)
corr_run.py       dedicated APE-smeared correlator runs
corr_analyze.py   smeared-correlator analysis
kappa.py          quotient-Hessian machinery (analytic gradient, projector, Lanczos)
kappa_run.py      kappa_eff on saved configurations
kappa_control.py  weak-coupling control (beta = 8..40)
kappa_summary.py  kappa tables and plot
block_obs.py      block-scale chars/defect clusters (Bianchi proxies)
block_analyze.py  -> block_summary.txt
t_up.csv          1/2/3-loop RG times at every beta (computed inline)
analyze.py        full analysis -> results.csv, summary.txt, plots/
runs/             raw data (.npz) and logs
corr/             smeared-correlator raw data
results.csv       all observables, errors, tau_int, bin sizes
summary.txt       human-readable summary
kappa_summary.csv kappa_eff table
plots/            figures
```

Reproduce:

```
python validate.py
python run.py --L 6 --beta 1.7 --sweeps 16000 --therm 1500 --meas 2 --loop_every 20 --seed 12 --out runs/L6_b1.700_s12.npz --snap_every 400 --snap_max 40
... (one command per row of the parameter table)
python corr_run.py --L 6 --beta 1.7 --sweeps 4000 --therm 1500 --seed 12 --out corr/L6_b1.700.npz
python kappa_run.py --Ls 6 --m 120 --nsnap 6 --tag L6
python kappa_control.py
python block_obs.py --L 6 --beta 3.0 --sweeps 3000 --therm 1500 --seed 53 --out runs/block_L6_b3.0.npz
python analyze.py && python corr_analyze.py && python kappa_summary.py && python block_analyze.py
```

---

## 9. Block-scale Bianchi / `z(q)` / `κ_c` measurements (requested by `handoff-proof`)

### 9.1 The object, and why it is not directly measurable here

The confirmed definition is
`z(q) = sup_x Σ_{C∋x} |w_q(C)| e^{+κq|C|}`, `κ = 1/2`, with
`|w_q(C)| ≤ e^{−c₀q|C|}`, `c₀ > κ`.

**Obstruction (structural, not statistical).**  `w_q(C)` are the polymer weights
of the **blocked effective action**.  A direct Monte Carlo of the measure
estimates expectation values; it does not return the coefficients of a polymer
expansion.  Producing `w_q(C)` requires either (i) inverse Monte Carlo / RG
matching to fit the blocked action, or (ii) an analytic character expansion of
the blocked Boltzmann factor, using a specific link-blocking scheme.  Neither
is implemented here.  **Therefore `κ_c`, `z_signed`, `z_abs` and the exact
`σ(1)` are not measured in this report.**  The closest computable proxies are
given below, each explicitly labelled as a proxy.

### 9.2 Character coefficients: is there ANY c_R < 0?  (proxy P1/P2)

`c_R = ⟨χ_R(U_p)⟩` for `R = 3, 6, 8, 10` (`χ₁₀ = [(tr U)³ + 3 tr U tr U² +
2 tr U³]/6`, `d₁₀ = 10`); jackknife errors over bins:

| β | L | c₃ | c₆ | c₈ | c₁₀ |
|---|---|---|---|---|---|
| 2.2 | 6 | +0.43041(17) | +0.09212(16) | +0.16522(24) | +0.01303(17) |
| 2.6 | 6 | +0.52107(19) | +0.13495(18) | +0.23809(26) | +0.02323(16) |
| 3.0 | 6 | +0.61476(19) | +0.18739(20) | +0.32620(27) | +0.03800(17) |
| 3.3 | 6 | +0.68861(19) | +0.23549(20) | +0.40444(29) | +0.05317(18) |
| 2.6 | 8 | +0.52064(23) | +0.13489(24) | +0.23824(33) | +0.02300(24) |
| 3.0 | 8 | +0.61527(26) | +0.18812(26) | +0.32677(40) | +0.03795(25) |

**Answer to the question: NO — no `c_R < 0` appears at any β or L.**  Every
coefficient is positive by ≫10σ (most by 10²σ).  Hence the single-plaquette signed/absolute
ratio is `r₁ = 1.0000` (same-sign sums coincide).

*Interpretation (as agreed with `handoff-proof`)*: this is **consistent with
BL-6** — in the strong-coupling surface representation `κ_c = 1` precisely
because the character coefficients are all positive.  The measurement shows
that the **sign structure a Bianchi cancellation would require is absent** in
the microscopic plaquette distribution at β=2.2–3.3, on both L=6 and L=8
(volume-independent).  It is still a **c₁-level statement, not a measurement
of `κ_c`**: `κ_c` involves blocked-lattice cluster weights, and `c_R > 0` for a
single plaquette does not by itself prove `κ_c = 1`.

### 9.3 Defect clusters (proxy P3)

Defect = site-averaged plaquette field `P(x) = (1/6)Σ_{μ<ν} ⅓Re tr U_{μν}(x)`
below a threshold; clusters are 6-connected components (non-periodic);
`n_x(m) = (1/V) Σ_{|C|=m} |C|` is the mean number of size-`m` clusters
containing a given site (so `Σ_m n_x(m)` = defect fraction).

| β | b | threshold | defect frac. | largest cluster | #clusters | n_x(1) | n_x(2) | n_x(3) |
|---|---|---|---|---|---|---|---|---|
| 2.2 | 1 | P<0 | 0.0963(2) | 7.1(1) | 85.4(1) | 0.0494 | 0.0194 | 0.0107 |
| 2.2 | 1 | P<⟨P⟩/2 | 0.2645(2) | 149(1) | 72.6(2) | — | — | — |
| 2.6 | 1 | P<0 | 0.0590(1) | 4.1(0) | 61.4(1) | 0.0394 | 0.0115 | 0.0046 |
| 2.6 | 1 | P<⟨P⟩/2 | 0.2241(2) | 63.3(5) | 88.5(2) | — | — | — |
| 3.0 | 1 | P<0 | 0.0336(1) | 2.7(0) | 38.6(1) | 0.0267 | 0.0051 | 0.0012 |
| 3.0 | 1 | P<⟨P⟩/2 | 0.1863(2) | 29.5(3) | 98.1(2) | — | — | — |

(`b=2`, 81 blocks: the defect fraction under `P<0` is 0 for all β; under
`P<⟨P⟩/2` it is ≤0.004.  A meaningful `b=3` analysis is impossible on `L=6`,
which has only 8 blocks of size 3.)

* At the strict threshold (`P<0`) the **site-averaged** defect gas is dilute
  and its clusters are small (largest 2.7–7.1 sites out of 1296).
* At the softer threshold (`P<⟨P⟩/2`) the largest cluster is 30–149 sites,
  i.e. the defects percolate — at that threshold they are **volume/domain-wall
  like, not a dilute gas of closed surfaces**.

The defect-cluster generating function `Σ_m n_x(m) e^{−qm/2}` (a proxy:
**not** `z(q)`, since the bare weights `w_q` are absent) is 0.0597 / 0.0408 /
0.0215 at `q = 0.5 / 1 / 2` (β=2.2).  Its decay is compatible with either
`e^{−cq}` or `e^{−cq²}` at this precision; as the analytic side noted,
`q²≈q` near `q≈1`, so **this cannot discriminate the Bianchi `q²` claim** and
no such judgment is drawn.

### 9.4 `σ(1)` per-block proxy (T2): the threshold dependence makes it fail

Requested measurement: the per-block defect fraction `f` at `b=2` (L=8, N=500
configs, β=2.6/3.0; L=6, N=3000), with **Poisson errors from the defect count
`k`** (not from the fraction), as a proxy for a defect action
`c₀ = −ln f`.  `k=0` cases get the 95 % bound `f < 3/N_tot`:

| β | L | threshold | k | f | c₀ = −ln f |
|---|---|---|---|---|---|
| 2.2 | 6 | P<0 | 0 | <1×10⁻⁵ | >11.30 |
| 2.2 | 6 | P<⟨P⟩/2 | 1055 | 0.00434(13) | 5.44(3) |
| 2.2 | 6 | P<0.75⟨P⟩ | 23247 | 0.09567(63) | 2.35(1) |
| 2.6 | 6 | P<0 | 0 | <1×10⁻⁵ | >11.30 |
| 2.6 | 6 | P<⟨P⟩/2 | 216 | 0.00089(6) | 7.03(7) |
| 2.6 | 6 | P<0.75⟨P⟩ | 14310 | 0.05889(49) | 2.83(1) |
| 2.6 | 8 | P<0 | 0 | <2×10⁻⁵ | >10.66 |
| 2.6 | 8 | P<⟨P⟩/2 | 114 | 0.00089(8) | 7.02(9) |
| 2.6 | 8 | P<0.75⟨P⟩ | 7707 | 0.06021(69) | 2.81(1) |
| 3.0 | 6 | P<0 | 0 | <1×10⁻⁵ | >11.30 |
| 3.0 | 6 | P<⟨P⟩/2 | 33 | 0.00014(2) | 8.90(17) |
| 3.0 | 6 | P<0.75⟨P⟩ | 8437 | 0.03472(38) | 3.36(1) |
| 3.0 | 8 | P<0 | 0 | <2×10⁻⁵ | >10.66 |
| 3.0 | 8 | P<⟨P⟩/2 | 21 | 0.00016(4) | 8.72(22) |
| 3.0 | 8 | P<0.75⟨P⟩ | 4633 | 0.03620(53) | 3.32(2) |
| 3.3 | 6 | P<⟨P⟩/2 | 10 | 0.00004(1) | 10.10(32) |
| 3.3 | 6 | P<0.75⟨P⟩ | 5498 | 0.02263(31) | 3.79(1) |

**Verdict: the proxy fails, exactly by the agreed criterion.**  For fixed
`(β, L, b)` the three thresholds give `c₀ = 2.4–3.8`, `5.4–10.1` and `>10.7` —
a factor 4–5 spread.  The "`c₀ ≈ 5.5`" figure corresponds to **one cell** of
this table (β=2.2, `P<⟨P⟩/2`); at other thresholds and couplings it is 2.4–10.

**Why — the blocked field is Gaussian.**  Computing the Gaussian tail from the
*same* mean and variance of `P_b` reproduces the measured fractions:
β=2.6, b=2, L=8: `μ = 0.17355`, `σ = 0.02805`; the Gaussian predicts
`f(⟨P⟩/2) = 0.00099` and `f(0.75⟨P⟩) = 0.0610`, measured `0.00089` and
`0.0602`.  Same at b=1 (predicted 0.0609 / 0.2195 / 0.3494; measured
0.0590 / 0.2241 / 0.3558).  Every "defect fraction" is therefore the
**Gaussian tail** fixed by `(μ, σ)`; `c₀ = −ln f` is a function of the
threshold and of the width of the blocked distribution, **not an independent
defect action**.

**Volume check**: L=6 and L=8 agree — `var(P₂) = 0.000781(2)` versus
`0.000787(3)`, and `c₀(⟨P⟩/2) = 7.026(68)` versus `7.024(94)` at β=2.6.  The
number is not a finite-volume artefact; it is the (Gaussian) width.

**Consistency with `R₄`**: a Gaussian blocked field is exactly what `R₄ ≈ 0`
(§4.4) reports, and it ties together: the smallness of `R₄` and the
"large c₀" are the same statement.

**Consequence**: the σ(1)-per-block proxy does **not** give `c₀`, and in
particular does not support "`c₀ ≈ 5.5`".  This is a negative result for the
proxy route, not for CPB.

### 9.5 T4: `κ_eff` — normalisation and mode (for completeness; not used)

`H = d²S_W/dξ² ∝ β`, so the β-independent comparison is `κ_eff/β`:

| β | 1.7 | 2.2 | 2.6 | 2.8 | 3.0 | 3.3 | 5.7 | 6.0 | 8 | 10 | 20 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| κ_eff/β | −1.14 | −1.13 | −1.09 | −1.08 | −1.06 | −1.04 | −0.70 | −0.58 | −0.28 | −0.20 | −0.086 | −0.043 |

Mode identification: the lowest Ritz value is the bottom of a **dense low-lying
cluster** (at β=2.2, L=6: −2.483, −2.450, −2.415), not a single zero-momentum
mode; the zero modes are projected out (and lifted by `Λ`).  Its L-independence
(L=4 vs L=6 agree within 2 %) identifies it as short-distance dominated, not a
finite-volume mode.  Per T4 it is **not used** in any `CPB` conclusion.

### 9.6 Cost table and the `ξ/a ≈ 7` feasibility question

Measured single-core cost (this machine, pure NumPy):

| L | sweep | per measurement (`meas=2`) | memory (links) |
|---|---|---|---|
| 6 | 0.045 s | 0.047 s | 0.75 MB |
| 8 | 0.120 s | 0.129 s | 2.4 MB |
| 10 | 0.42 s | 0.45 s | 5.8 MB |
| 32 (extrapolated, ∝L⁴) | 38 s | 38 s | 76 MB |

`ξ/a ≈ 7` requires `L/ξ ≳ 4`, i.e. `L ≳ 28`; at `L=32`, `10⁴` measurements
would take **≈105 h** on one core; the statistics
actually needed for 5 % on `ξ` at `ξ/a=7` (`N_eff ≈ 6×10⁵`, §10.3) cost
~1.7×10⁵ CPU-hours at `L=32` with the unoptimised code — hence the "out of
reach" in §0.  §10 gives the optimised-kernel projection and the anisotropic
volumes.  **This is not reachable with the present
pure-NumPy code.**  It is a feasibility limitation, not numerical noise.

Options, in order of cost:
1. **Declare supportive-only** (recommended): keep `L ≤ 8` runs as supporting
   numerics, without a claim on the σ certificate.
2. **Finite-size scaling** at `L = 12–16` (0.7–2.4 s/measurement): feasible in
   ~10⁴ s of CPU per coupling, but cannot reach `ξ/a = 7` directly — only its
   trend.
3. **Faster kernel** (C/OpenMP/Numba, ≥20×) would bring `L=32` to ~5 h per
   10⁴ measurements; a research-grade project, not this session.

**Feasibility.**  With the rigorous 2-loop `t_up = 2.639` at `λ = 6` the σ
cap is `q₁^σ ≈ 0.143`, demanding `ξ/a ≈ 7`.  This study gives `ξ/a ≈ 0.93(8)`
at β=6 and `≲0.3` at β ≤ 3.3 — two orders of magnitude short.  §10.5 revises
the cost with multi-chain parallelism and the precision–margin coupling:
`@7×` is infeasible, `@50× × 16 chains` costs 56–222 h at 5 % (2–9 days) or
3.5–14 h at 20 % — but **20 % precision only certifies the cap if the measured
`ξ/a` is ≳11.7**, and measuring `ξ/a = 7` itself needs `L ≈ 21–28`, so the
proposed L=8/12 scan cannot reach it.  Reaching
`ξ/a ≈ 7` requires **either** (i) a kernel at least 20× faster (C/OpenMP/Numba;
then `L=32` costs ~5 h per 10⁴ measurements), **or** (ii) a different
observable whose finite-size scaling fixes `ξ/a` at `L=12–16` (e.g. a
fixed-aspect-ratio scaling study or the Polyakov-loop torelon spectrum).
Neither is available in this pure-NumPy code.

### 9.7 Judgment

**The central conclusion is a structural one: `κ_c` has no direct numerical
route from a direct Monte Carlo of the measure.**  The polymer weights
`w_q(C)` are coefficients of the *blocked effective action*; a MC samples the
measure and returns correlation functions, not those coefficients.  Getting
`w_q(C)` needs inverse Monte Carlo / RG matching or an analytic character
expansion of a link-blocked Boltzmann factor.  This is **the result** for T1,
not a preliminary.

What the measurable proxies give:

* **Character coefficients (premise of the mechanism): NO `c_R < 0`** at any
  β = 2.2–3.3 or L = 6, 8 — all of `c₃, c₆, c₈, c₁₀` positive by ≫10σ.  This
  is **consistent with BL-6** (`κ_c = 1` because the coefficients are
  positive): the sign structure a Bianchi cancellation would need is absent at
  the single-plaquette level.  It does **not** prove `κ_c = 1`, which lives on
  blocked clusters.
* **σ(1) per-block proxy: FAILS.**  The threshold dependence alone spans
  `c₀ = 2.4–10.1` (plus `>10.7` at `P<0`), and the measured fractions equal
  the Gaussian tail of the blocked field — so `−ln f` measures the width, not
  an action.  The `c₀≈5.5` figure is one cell of that table.
* **Defect geometry**: dilute small clusters at `P<0`, percolating at
  `P<⟨P⟩/2`; no closed-surface/domain-wall discrimination is possible with
  this data.

**Overall: INCONCLUSIVE for the Bianchi/`z(q)` question**, with two concrete
negative results (no sign structure; the σ proxy is not a measurement) and one
structural obstruction (`w_q` inaccessible).  Both negatives point the same way
as the analytic default, but neither is a measurement of `κ_c`.

---

## 10. Fast kernel and cost projection for `ξ/a ≈ 7` (requested by `handoff-proof`)

### 10.1 The compiled kernel (Numba)

`su3_numba.py` implements the same CM heat bath (same staples, same SU(2)
sampler, same three subgroup updates) as `su3.py`, as a Numba `@njit` kernel
with per-site loops.

**Validation**
* `plaquette_mean` and `action_numba` vs the NumPy versions on an identical
  random configuration: agreement to 2×10⁻¹⁷ and 0 (machine precision).
* Staple sum at 20 random links, Numba vs NumPy: max difference 9.9×10⁻¹⁶.
* Equilibrium `⟨P⟩` at β=2.2, L=6 (independent cold starts, 300 measurements
  each): Numba 0.14334(19) vs NumPy 0.14346(24) — consistent.
* **Bitwise trajectory comparison is impossible by construction** (different
  RNG streams); the deterministic kernels above are the meaningful check.

**Measured speed (this machine, β=2.2).**

Single core:

| | ms / sweep at L=6 |
|---|---|
| NumPy (`su3.py`, isolated) | 35.0 |
| Numba, no fastmath | 7.75 |
| Numba + fastmath | 4.97 |
| **serial speedup** | **7.0×** |

Parallel (`prange` over the flattened `(x0,x1)` pair; sites of one parity are
independent, per-thread workspace; 32 logical threads):

| L | serial Numba | parallel Numba | parallel speedup | **vs single-core NumPy** |
|---|---|---|---|---|
| 6 | 5.12 ms | 0.63 ms | 8.2× | **55×** |
| 8 | 15.2 ms | 1.50 ms | 10.2× | **74×** |
| 12 | 77.2 ms | 7.15 ms | 10.8× | **78×** |

**The ≥50× target is met without C:** 55–78× versus single-core NumPy,
rising with volume.  (The earlier `prange` failure was the `continue` inside
the parallel loop; the working version uses `if`, per-iteration workspaces and
a flattened two-index `prange`.)

Additional validation: APE smearing in Numba (`ape_smear_num_par`, Newton
polar projection, 6 iterations with RMS scaling) matches a NumPy reference
with the *same* SU(3)-projection convention (`project_su3_cuberoot`) to
3×10⁻¹¹; the parallel and serial versions agree exactly.  (Note: the earlier
correlator runs used the row-phase projection convention; the two differ by a
per-link U(1) phase.  Gauge-invariant observables agree to the level checked
in the β=6.3 pilot (§10.8).)

### 10.2 Cost projection

**Note (superseded — read §10.10 first):** the "16 chains" column below
assumed the per-chain speedup could be multiplied by the number of chains.
The measured aggregate throughput of the whole machine is **~50×
numpy-core**, so those wall-clock entries are **16× too optimistic**.  The
corrected cost is:

```
5 % + GEVP : 135–1,160 h  (1–7 weeks), depending on the N_eff calibration (7e4 vs 6e5)
20 % + GEVP : 8–72 h
20 % is only useful for a 2σ certificate if the measured ξ/a ≳ 11.7.
```

**The "3 days" figure is withdrawn.**

Per-site cost measured on NumPy: 27 µs (L=6), 29 µs (L=8), 42 µs (L=10) —
i.e. ∝V with a cache penalty that grows with the volume; we use 50 µs/site
(NumPy) for L ≥ 16, i.e. **7.1 µs/site** with the 7× kernel.

| volume | sites | sweep, NumPy | sweep, 7× | sweep, 50× (hypothetical) | L_t/ξ | L_s/ξ |
|---|---|---|---|---|---|---|
| 16³×64 | 262 144 | ~13 s | **1.9 s** | 0.27 s | 9.1 | 2.3 |
| 24³×64 | 884 736 | ~44 s | 6.3 s | 0.90 s | 9.1 | 3.4 |
| 32⁴ | 1 048 576 | ~52 s | 7.5 s | 1.05 s | 4.6 | 4.6 |

### 10.3 Statistics needed for 5 % on `ξ` at `ξ/a = 7`

Calibrated on the β=6 measurement (`ξ/a = 0.93`, 4000 measurements,
`δm/m = 8.2 %` with the APE-smeared correlator).  With the fit point scaling
as `t* ~ ξ/2`, the correlator error gives `δm/m ∝ ξ/√N_eff`, so

```
N_eff(ξ=7, 5%) = 4000 × (7/0.93)² × (8.2/5)² ≈ 6×10⁵ independent configurations
```

(crude, factor ~3 either way; a GEVP variational operator would reduce it).
With a glueball autocorrelation time `τ_int ≈ 10` sweeps, this is
`N_sweeps ≈ 2 τ N_eff ≈ 1.2×10⁷` sweeps.

**Projected wall-clock (single core):**

| volume | 7× kernel | 50× kernel |
|---|---|---|
| 16³×64 | ~6 300 h | ~900 h |
| 24³×64 | ~21 000 h | ~3 000 h |
| 32⁴ | ~25 000 h | ~3 500 h |

(§9.6 quoted ~105 h for `10⁴` *measurements* at `L=32` and "out of reach"
for the full statistics; the table above is the full-statistics projection.)

### 10.4 Feasibility verdict

**`ξ/a ≈ 7` is not reachable**, not at 7× and not at a hypothetical 50×:
the cheapest configuration (16³×64, 50× kernel) still needs ~900 CPU-hours,
and the anisotropic volume has `L_s/ξ = 2.3`, i.e. a ~10 % finite-size
systematic on the mass — already comparable to the 5 % target.  `24³×64`
fixes the finite-size issue (`L_s/ξ = 3.4`) but costs ~3 000 h even at 50×.

This is a **feasibility negative**, independent of numerical noise.  What
could change it: (i) a cheaper observable for `ξ` (fixed-aspect-ratio scaling
at `L = 12–16`, where `ξ/a ≈ 2–3` is reachable, plus an extrapolation);
(ii) a many-order-of-magnitude faster kernel (GPU); (iii) accepting a much
larger statistical error (e.g. 20 % on `ξ` needs ~25× fewer measurements,
i.e. ~250 h at 50×, still a multi-week run).

### 10.5 Multi-chain parallelism and the precision–margin coupling

**Chains are independent** for pure gauge Wilson loops, so 16 chains on 16
cores divide the wall-clock by 16 (the Numba kernel is single-threaded per
chain; `τ_int` is a within-chain quantity).  Corrected projection:

| volume | @7× CPU-hours | @7× × 16 chains | @50× × 16 chains | @50×,×16, 20 % error | @50×,×16, 20 %, GEVP ×3 |
|---|---|---|---|---|---|
| 16³×64 | 6 400 | **400 h** | **56 h** | 3.5 h | 1.2 h |
| 24³×64 | 21 300 | **1 330 h** | **187 h** | 11.7 h | 3.9 h |
| 32⁴ | 25 400 | **1 590 h** | **222 h** | 13.9 h | 4.6 h |

**GEVP / variational operators** (unmeasured here, planning value): a
3–4-operator variational basis typically raises the ground-state overlap by
1.5–3×, i.e. reduces `N_eff` by 2–9×; we adopt **×3** as the central planning
number.

**How much precision is actually needed?**  The certificate wants
`q₁^σ = a/ξ ≤ 0.143`, i.e. `ξ/a ≥ 7`.  For a **2σ** certification of a measured
`ξ/a ± δ`:

| precision δ | needs measured `ξ/a ≥` |
|---|---|
| 5 % | **7.8** |
| 10 % | **8.8** |
| 20 % | **11.7** |

So **20 % precision is only sufficient if the measured `ξ/a` lands at ≳11.7**;
a 20 % measurement near 7 does not certify anything.  This must be settled
with `handoff-proof` before choosing the statistics.

**Where `ξ/a = 7` lies, and the finite-size trap.**
`ξ/a(β) ∝ e^{β/(2b₀)}` calibrated at β=6 (`ξ/a = 0.93`) gives
`ξ/a = 7` at **β ≈ 6.28**, `11.7` at β ≈ 6.35 (two-loop, indicative only —
scaling violations at β≈6 are 30–50 %).  But to *measure* `ξ/a = 7` one needs
`L ≳ 3–4 ξ`, i.e. `L ≈ 21–28`: on **L=12, `L/ξ` is only 1.7 at ξ=7**, so the
proposed `β ∈ [6.0, 6.6]` scan on L=8/12 **cannot reach ξ/a=7** — it measures
`ξ/a(β)` only up to ≈3–4 (where `L/ξ ≳ 3`) and then runs into finite-size
saturation.  The scan is still useful (it fixes the scaling curve), but the
ξ/a=7 point itself requires the 24³×64 / 32⁴ volumes.

### 10.6 Revised feasibility statement

* `@7×` single-core: infeasible (400–1 600 h even with 16 chains).
* `@50× × 16 chains`, 5 % precision: **56–222 h** (2–9 days) — borderline.
* `@50× × 16 chains`, **20 % precision**: **3.5–14 h** — feasible in a day,
  **but only useful if the measured `ξ/a` is ≳11.7** (see above).
* Adding a GEVP (×3) divides these by three.
* A C+OpenMP kernel at ≥100× would cut the `@50×` column by 2×.

Nothing is launched: this section is a projection.

### 10.7 The `ξ/a ≥ 8.5`, `L_s/ξ ≥ 4` scenario (32³×64)

If a 2σ lower bound must clear `ξ/a = 7` at 5 % precision, the target value is
`ξ/a ≥ 8.5`, and `L_s/ξ ≥ 4` then requires `L_s ≥ 34`, i.e. **32³×64**
(2 097 152 sites).  Using `N_eff = 6.1×10⁵` (5 %), `τ=10`, `N_sweeps = 1.2×10⁷`,
50 µs/site (NumPy) and a 50× kernel:

| volume | sites | `U` per chain | 16 chains | sweep @50× | CPU-h (5 %) | wall @50× ×16 | + GEVP ×3 | `L_s/ξ` at ξ=8.5 | finite-size ~`e^{−L_s/ξ}` |
|---|---|---|---|---|---|---|---|---|---|
| 16³×64 | 262 144 | 0.15 GB | 2.4 GB | 0.26 s | 890 | 56 h | 19 h | 1.9 | 15 % |
| 24³×64 | 884 736 | 0.51 GB | 8.2 GB | 0.88 s | 3 000 | 187 h | 62 h | 2.8 | 6 % |
| 32⁴ | 1 048 576 | 0.60 GB | 9.7 GB | 1.05 s | 3 550 | 222 h | 74 h | 3.8 | 2.3 % |
| **32³×64** | 2 097 152 | **1.21 GB** | **19.3 GB** | **2.10 s** | **7 100** | **444 h** | **148 h** | 3.8 | **2.3 %** |

**Finite-size systematic.**  The leading image correction is `e^{−L_s/ξ}`:
`L_s/ξ = 3` (24³×64 at ξ=8) gives `e^{−3} = 5.0 %`; `L_s/ξ = 2.8` (at ξ=8.5)
gives 6 %; `L_s/ξ = 3.8` (32³ or 32⁴) gives 2.3 %.  The periodic cosh fit used
here *already subtracts the leading image*, so the residual is at the lower end
of these numbers (~2–5 % at `L_s/ξ = 3`).  **24³×64 at `L_s/ξ≈3` is therefore
borderline: its ~2–5 % systematic is comparable to the 5 % statistical target.**
The cost ratio 32³×64 / 24³×64 is **2.4×** (32³×64 / 16³×64 is 8×).

**Memory is a real constraint:** 16 chains of 32³×64 need **19.3 GB** of link
matrices alone (this machine has 33.5 GB), plus staples and the OS — tight.
24³×64 needs 8.2 GB and is comfortable.

**Summary of the decision inputs:**
* 5 % + 2σ-certificate ⇒ ξ/a ≥ 8.5, L_s/ξ ≥ 4 ⇒ **32³×64, 444 h** at 50× with
  16 chains (**148 h** with GEVP ×3), memory-tight.
* 24³×64 (187 h, 62 h with GEVP) is usable only if a ~2–5 % finite-size
  systematic is acceptable, i.e. if the required 2σ bound can tolerate
  `ξ/a` ≈ 9 rather than 8.5.
* 16³×64 (56 h) has `L_s/ξ ≈ 2`, a ~10 % systematic — not compatible with a
  5 % target.

### 10.9 Production protocol (agreed with `handoff-proof`)

* **Same β for both volumes**: the pilot (`L=24`) and the production (`32⁴`)
  run at **β = 6.30**, so that `ξ/a(L=24)` and `ξ/a(L=32)` are directly
  comparable.  A single volume is not enough to quote a finite-size
  systematic.
* **Two-point consistency check**: if the two values differ by more than
  2–3 % (beyond their statistical errors), add a third volume or postpone the
  conclusion.  If they agree, the `e^{−L/ξ}` estimate of §10.7 is verified
  directly.
* **Direction of the bias**: a finite box truncates the long tail of the
  correlator, so a finite-volume `ξ` is expected to come out **small**.
  The measured `ξ/a` is therefore a *conservative* lower bound for the
  infinite-volume value — which is the safe direction for a certificate of
  the form `ξ/a ≥ 7`.  This must be confirmed by the two-point comparison,
  not assumed.
* **Target**: `ξ/a = 8.5 ± 5 %`, giving a 2σ lower bound `8.5 − 0.85 = 7.65 > 7`.
* **Reporting**: point value, statistical error, 2σ lower bound, `τ_P`, and
  the **effective** number of independent configurations `N_eff = N_meas/(2τ)`
  — not the nominal sweep count.  If the 2σ lower bound is below 7, the
  report must say **"insufficient evidence"**; no parameter tuning to force
  it.
* **Decision rules (fixed in advance, to avoid post-hoc rationalisation)**:
  (i) if the pilot pins `N_eff(5 %, ξ=8.5) ≈ 7×10⁴`, then 5 % + GEVP costs
  ≈135 h (5–6 days) and **is worth running**; (ii) if it pins
  `N_eff ≈ 6×10⁵`, the cost is ≈1,160 h (7 weeks) and **is not worth it** —
  do a 20 % run instead and state in the write-up that certificate-grade
  evidence needs ~7 weeks; (iii) if the pilot shows `ξ/a` outside 7–10,
  adjust β *before* any 32⁴ run.
* **Role of the numbers**: this is **evidence, not a certificate**.
  `q₁^σ ≤ 0.143` is a rigorous inequality; a numerical measurement, however
  precise, cannot prove it.  **Numerics do not constitute the certificate.**

### 10.8 Pilot at L=24, β=6.30 (measured)

`pilot_numba.py --L 24 --beta 6.30 --sweeps 3000 --therm 400 --meas 10 --nlev 8`,
16 threads, 1 726 s.  `⟨P⟩ = 0.6222`, `τ_int(P) = 1.59`,
`N_meas = 300`, **`N_eff = N_meas/(2τ) ≈ 95`**.

| level | C(1)/C0 | C(2)/C0 | C(3)/C0 | C(4)/C0 | cosh fit (t=1–10) |
|---|---|---|---|---|---|
| n=0 | 0.0767(124) | −0.007(112) | −0.010(124) | 0.018(120) | not constrained |
| **n=8** | **0.258(141)** | **0.084(137)** | **0.050(120)** | **0.033(129)** | **`am = 0.886(159)`, `ξ/a = 1.13(20)`** |

The fit is stable under the window: `ξ/a = 1.01 / 1.08 / 1.10 / 1.11 / 1.12 /
1.13` for `t_max = 3…10`.  `am_eff(t) = 1.37, 1.21` at t=1,2 (then noise).
The unsmeared operator is much worse (`C(1)/C0 = 0.077`), confirming the need
for smearing.

A 2×2 GEVP (levels n=0 and n=8) gives `m_eff` between 0.54 and 1.15 over
t=2..6 (`ξ/a = 0.9–1.9`), noise-dominated with 300 configurations; there is
**no evidence of a light state hidden below `am ≈ 0.9`**.

**So at β = 6.30, L=24: `ξ/a ≈ 1.1`, not ≈8.**

### 10.8b The β extrapolation had a factor-6 error in the exponent

`ξ/a ∝ a⁻¹ ∝ exp[1/(2b₀g²)]`; with `g² = 6/β` this is
`exp[β/(12b₀)]`, whose β-slope is `1/(12b₀) = 1.196` — **not** `1/(2b₀) =
7.18`.  The earlier estimate `ξ/a(β) = 0.93 e^{(β−6)/(2b₀)}` therefore
over-predicted the growth by a factor `e^{6×}` per unit β.

*(Equivalently: the correct two-loop `a⁻¹` gives `a⁻¹(6.3)/a⁻¹(6.0) = 1.46`;
the naive formula gave 8.6.)*

Corrected, with the two-loop `a⁻¹(β)` calibrated at β=6 (`ξ/a = 0.93`):

| β | 6.0 | 6.3 | 7.0 | 7.5 | **7.8** | 8.0 |
|---|---|---|---|---|---|---|
| `ξ/a` (predicted) | 0.93 | **1.36** | 3.28 | 6.15 | **8.95** | 11.5 |

The pilot's `ξ/a = 1.13(20)` at β=6.3 agrees with the corrected prediction
1.36 (and is slightly below it, the expected direction for a finite volume /
imperfect overlap).  **This confirms the correction.**

**Consequence for the approved plan:** a 32⁴ run at **β = 6.30 would measure
`ξ/a ≈ 1.3`, not 8.5** — it would not reach the certificate target at all.
The target `ξ/a = 8.5` sits at **β ≈ 7.75–7.8**; `ξ/a = 7` at β ≈ 7.6.
This must be fixed before any production run.  Note that the *lattice volume*
needed is unchanged (same `L/ξ`), but at β≈7.8 the autocorrelation times are
larger, so the statistics estimate of §10.3 must be redone after measuring
`τ` at that coupling.

### 10.10 Measured aggregate throughput — the 16-chain projection was 16× too optimistic

Per-chain speedups cannot be multiplied by the number of chains: a parallel
kernel already uses all cores.  Measured on this machine (L=12, 32 logical
cores), single chain:

| threads | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| ms/sweep | 109 | 59.4 | 31.0 | 16.5 | 14.5 | 13.5 |

and concurrent chains partitioning the same 32 threads (aggregate, worst-chain
timing, 20 sweeps each):

| K chains × T threads | 1×32 | 2×16 | 4×8 | 8×4 | 16×2 | 16×1 |
|---|---|---|---|---|---|---|
| aggregate sweeps/s (L=12) | 89.9 | 89.8 | 63.7 | 50.1 | 69.7 | 96.3 |
| **= × numpy-core** | **50×** | 50× | 36× | 28× | 39× | **54×** |

**The machine saturates at ≈50× numpy-core aggregate**, essentially
independent of how the threads are partitioned.  The earlier projection
("16 chains × a 50× per-chain kernel = 800×") double-counted the cores and is
**16× too optimistic**.

**Revised cost** (`32⁴`, work = `N_sweeps × 52 s` of one numpy core, divided by
the measured 50× aggregate):

| assumptions | `N_sweeps` | wall, no GEVP | with GEVP ×3 |
|---|---|---|---|
| 5 %, `N_eff=6×10⁵`, τ=10 | 1.2×10⁷ | ~3 500 h | ~1 160 h |
| 5 %, `N_eff=7×10⁴`, τ=10 | 1.4×10⁶ | ~400 h | ~135 h |
| 20 %, `N_eff=6×10⁵`, τ=10 | 7.5×10⁵ | ~220 h | ~72 h |
| 20 %, `N_eff=7×10⁴`, τ=10 | 8.8×10⁴ | ~25 h | ~8 h |

The two `N_eff` values come from two calibrations of the fit error:
`δm/m ∝ ξ/√N_eff` fitted to the β=6 run (4000 meas, ξ=0.93, 8.2 %) gives
`N_eff(5 %, ξ=8.5) = 6×10⁵`; fitted to the β=6.3 pilot (95 effective configs,
ξ=1.13, 18 %) gives `7×10⁴`.  The pilot-based calibration is more direct (same
operator and fit) but rests on only ~95 configurations.  **The β=7.8 pilot
running now will pin this down**, since it measures `ξ` and its error at the
relevant correlation length.

Measured at β=7.8, L=24 (1 500 sweeps): `⟨P⟩ = 0.7122`, **`τ_P = 1.07`** — the
plaquette autocorrelation does *not* grow at this coupling.  The glueball
correlator's `τ` will be extracted from the pilot.
