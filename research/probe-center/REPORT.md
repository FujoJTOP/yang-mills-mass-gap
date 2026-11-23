# REPORT — Center-vortex probe (NP-7): can Z(3) vortices supply `κ_c < 1`?

> 2026-11-22 · session `probe-center`, task from the main Yang–Mills session.
> Scope: NP-7 in `research/new-paths/SURVEY.md`. Does **not** change `gluing`
> R-numbers, `bianchi` BL-numbers, or `handoff` Prop B.
> Deliverable numbers: `probe_*.npz`, `analyze.py`, `plots/`.

---

## 0. Verdict (read first)

**`κ_c < 1` cannot come from center vortices for the object `CPB-W` needs.
NP-7 is closed as a route to the mass gap.**

Three independent reasons, in order of strength:

1. **Exact (algebra, not numerics).** A center vortex multiplies a Wilson loop
   in representation `R` by `χ_R(z)/d_R`, `z ∈ Z(3)`. For every representation
   with **N-ality 0** — the adjoint, and *every* colour-singlet (gauge-invariant
   local) operator — this factor is exactly `1`. So the center-vortex phase
   carries **zero** information into the neutral sector. The mass gap lives on
   neutral local operators; a mechanism that is exactly neutral-blind cannot
   supply their sign cancellation. Numerically confirmed: the centre-projected
   adjoint loop is `1` to machine zero (§4, `max|W_A−1| = 0`).

2. **Numerical (sector separation).** Maximal-centre-gauge (MCG) projection
   reproduces the **fundamental** (charged, N-ality 1) string tension
   (centre dominance, `χ_proj/χ_full ≈ 1.5–1.7` at (2,2), β=6) and is
   **identically trivial** in the adjoint. The adjoint's short-range area law
   (Casimir-scaled, `χ_A(2,2) ≈ (9/4)·χ_F(2,2)`) is therefore **entirely
   outside** the vortex sector. This is the "signal only in the charged sector"
   branch of the task's closing criterion, demonstrated directly.

3. **Measured (signed/absolute activity).** Under the physical (real-phase)
   convention the ratio is `−1/2` for **every** vortex plaquette and **every**
   cluster, independent of size; under the orientation convention it is `≈ 0`
   (CP symmetry), also size-independent. Neither is a size-dependent dynamical
   cancellation. At the project's couplings (β = 2.6–3.3) the projected ensemble
   is **fully percolating** (one spanning cluster per configuration), so there
   is no dilute cluster-size distribution for a "decay with cluster size" to
   live in at all (§3).

**Discipline — "vortices exist" ≠ "vortices give a sign cancellation".**
Vortices provably *exist* here (density, percolation, centre dominance; §2) and
are the confinement mechanism of the **charged** sector. That is not the same
as providing the neutral-sector sign cancellation `CPB-W` requires, and the
data show they do not.

**Caveat on the word "closed".** Following `SURVEY.md` §1.2 this is a
*mechanism-insufficient / no-known-route* verdict, **not** "proved false":
centre vortices as a sign carrier are excluded for the neutral sector by
group theory; whether the *vortex ensemble* dynamically contributes to a mass
gap by some other route is not addressed here.

---

## 1. Method

### 1.1 Configurations

Reused `research/lattice-cpb/runs/` (4D pure SU(3) Wilson, periodic, L⁴):

| β range | role | L | source |
|---|---|---|---|
| 2.6, 2.8, 3.0, 3.3 | **primary grid** (project handoff) | 6 | `L6_b*_s5*.npz` (40 snapshots) |
| 1.7, 2.2 | conservative strong coupling | 6 | `L6_b*_s1*/s2*.npz` (40) |
| 5.7, 6.0 | scaling-region control | 6 | `L6_b*_s3*.npz` (30) |
| 1.7, 2.2, 6.0 | finite-volume | 8 | `L8_b*_s*.npz` (20) |
| 6.0 | finite-volume | 10 | `L10_b6.0_s41.npz` (5) |
| **2.6, 3.0** | **finite-volume, generated here** | **8** | `runs_L8_b*.npz` (40) |

The L=8 primary-grid couplings were **not** in `runs/`; they were generated
here with the same Cabibbo–Marinori heat bath (`su3_numba.sweep_par`, κ = β/3),
warm-started from the equilibrated L8 β=2.2 configuration, 2500 thermalisation
sweeps + 40 snapshots spaced 20 sweeps. Reproduced ⟨P⟩ ≈ 0.17 at β=2.6 (matches
L6) — a κ=β/6 slip was caught and corrected.

### 1.2 Maximal centre gauge (MCG)

Maximise `R[g] = Σ_{x,μ} |Tr_F U_μ(x)|²`. This functional **equals** the adjoint
Landau functional, since `Tr_adj(M) = |Tr_F M|² − 1` for `M ∈ SU(3)`
(so MCG = adjoint Landau = "indirect MCG"). Per-site maximisation by Riemannian
ascent: for `f(g)=Σ_a|Tr(g C_a)|²` the Wirtinger gradient is
`∂f/∂g^† = Σ_a Tr(g C_a) C_a^†`; iterate `g ← polar(gradient)` with
monotone acceptance. Two restarts (identity, random), keep larger `R`.
Convergence: `R` rises monotonically and saturates (e.g. β=6, L6: `R: 0.99→4.6`
over ~40 sweeps). Validation: `W_F(1,1) = ⟨P⟩ = 0.174` at β=2.6 (table), and the
raw-vs-fixed vortex density drops `0.64 → 0.27` at β=6.

### 1.3 Z(3) projection and vortices

`U_μ(x) → z_μ(x) I`, `z = ω^k`, `k = round(3θ/2π) mod 3`, `θ = arg Tr U_μ(x)`
(ambiguous if `|Tr U| < 0.05`; fraction `≤ 2×10⁻⁴`). Plaquette flux
`k_{μν}(x) = k_μ(x)+k_ν(x+μ)−k_μ(x+ν)−k_ν(x) (mod 3)`; nonzero flux on a lattice
plaquette = one vortex plaquette (a dual plaquette). Vortex worldsheets are
closed 2-chains (`∂ = 0` by the Bianchi identity on the Z(3) links).
**Clusters** = connected components of vortex plaquettes sharing a lattice link
(edge-adjacency of surface tiles).

### 1.4 Observables

* Fundamental Wilson loop `W_F(R,T) = ⟨Re Tr/3⟩`; adjoint
  `W_A = (9·⟨|Tr/3|²⟩ − 1)/8` (using `Tr_adj = |Tr_F|²−1`), from full and
  projected links.
* Fundamental (charged) Polyakov correlator `⟨P(x)P*(x+r)⟩/9`; adjoint (neutral)
  connected `⟨a(x)a(x+r)⟩ − ⟨a⟩²`, `a = (|Tr P|²−1)/8`.
* Vortex density, cluster sizes, per-cluster `k=1` vs `k=2` counts.

---

## 2. Vortex existence and structure

Measured on all snapshots (mean ± s.e.m. of the density; counts are means):

| β | L | nsnap | `R` | density | `n_vortex` | `n_clusters` | k=1 | k=2 | D2 | max cluster |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.7 | 6 | 40 | 3.82 | 0.618±0.0010 | 4806 | 1.00 | 2393 | 2413 | −0.0042 | 4955 |
| 2.2 | 6 | 40 | 3.84 | 0.602±0.0010 | 4685 | 1.00 | 2344 | 2341 | +0.0008 | 4762 |
| **2.6** | **6** | 40 | 3.85 | **0.589±0.0008** | 4579 | 1.00 | 2283 | 2296 | −0.0027 | 4657 |
| **2.8** | **6** | 40 | 3.85 | **0.582±0.0007** | 4526 | 1.00 | 2259 | 2267 | −0.0018 | 4594 |
| **3.0** | **6** | 40 | 3.86 | **0.573±0.0010** | 4454 | 1.00 | 2238 | 2217 | +0.0046 | 4577 |
| **3.3** | **6** | 40 | 3.88 | **0.564±0.0009** | 4384 | 1.00 | 2194 | 2189 | +0.0012 | 4495 |
| 5.7 | 6 | 30 | 4.33 | 0.355±0.0018 | 2762 | 1.10 | 1378 | 1384 | −0.0022 | 2923 |
| 6.0 | 6 | 30 | 4.52 | 0.301±0.0028 | 2342 | 1.00 | 1172 | 1170 | +0.0007 | 2590 |
| 1.7 | 8 | 20 | 3.82 | 0.618±0.0007 | 15194 | 1.00 | 7613 | 7581 | +0.0021 | 15377 |
| 2.2 | 8 | 20 | 3.83 | 0.601±0.0008 | 14760 | 1.00 | 7387 | 7373 | +0.0009 | 14951 |
| **2.6** | **8** | 40 | 3.84 | **0.589±0.0006** | 14473 | 1.00 | 7251 | 7222 | +0.0020 | 14630 |
| **3.0** | **8** | 40 | 3.85 | **0.574±0.0004** | 14103 | 1.00 | 7061 | 7042 | +0.0013 | 14296 |
| 6.0 | 8 | 20 | 4.48 | 0.309±0.0015 | 7600 | 1.15 | 3796 | 3804 | −0.0010 | 7883 |
| 6.0 | 10 | 5 | 4.47 | 0.311±0.0004 | 18674 | 1.40 | 9309 | 9366 | −0.0031 | 18749 |

**L=6 vs L=8 agree to <1 %** at every common β (finite-volume control OK);
the independently generated L8 β=2.6/3.0 reproduce the L6 density to
0.589 vs 0.589 and 0.574 vs 0.573.

Interpretation:
* Vortices **exist** and are abundant at every β. Density falls with β
  (strong coupling → dense), and the ensemble **percolates**: on every
  configuration analysed there is **one spanning cluster** carrying essentially
  all vortex plaquettes. This is the expected confined-phase vortex condensate.
* Consequently NP-7's "cluster-size distribution with a decaying activity"
  has no regime to act in: there is no dilute population of finite clusters on
  L = 6–10 at β = 1.7–6.0. (`plots/cluster_sizes.png`.)

---

## 3. Signed / absolute activity ratio `z_signed/z_abs`

Three explicit definitions, because the phrase is ambiguous and only one can be
a `κ_c` proxy:

**(D1) Physical real-phase sign.** A vortex plaquette carries centre phase
`ω^k`; its contribution to a real (e.g. `Re Tr` fundamental) observable is
`Re ω^k = −1/2` for both `k = 1, 2`. Thus
`z_signed/z_abs = Σ Re ω^{k_p} / Σ |ω^{k_p}| = −1/2`, **exactly, for every
cluster and every size.** This is the group-theoretic number `cos(120°)`; it is
not a cancellation (all vortex weights add coherently in Re).

**(D2) Orientation / vacuum asymmetry.** `sign(k=1)=+1`, `sign(k=2)=−1`;
per cluster `(N₁−N₂)/(N₁+N₂)`. Measured ensemble value **≈ 0** (equal
vortex/antivortex populations, CP symmetry), size-independent.

**(D3) Effective-action polymer `κ_c`.** The object `CPB-W` needs — signed /
absolute ratio of blocked-action polymer weights `w_q(C)` — is **not
measurable by direct MC** (STATUS §9, and `lattice-cpb/REPORT.md` §9.1). We do
not claim to measure it. We only note the *only channel* by which a vortex can
change a real sign is the Z(3) phase (D1), which is neutral-blind (§4).

Measured (ensemble over snapshots; D2 per cluster, per configuration):

| β | L | D1 (real phase) | D2 (orientation) | cluster size bins populated |
|---|---|---|---|---|
| 2.6 | 6 | −0.5000 | −0.0027 | only `s>256` (one percolating cluster) |
| 3.0 | 6 | −0.5000 | +0.0046 | only `s>256` |
| 3.3 | 6 | −0.5000 | +0.0012 | only `s>256` |
| 2.6 | 8 | −0.5000 | +0.0020 | only `s>256` |
| 3.0 | 8 | −0.5000 | +0.0013 | only `s>256` |
| 6.0 | 8 | −0.5000 | −0.0010 | `s≈256–10⁴`; a few extra finite clusters (n_cl≈1.15) |

For completeness the *per-plaquette* orientation variance is `1` (k=±1) and
the *per-cluster* `|N₁−N₂|/(N₁+N₂)` is `≈ 0.007–0.012` for the percolating
cluster — the `1/√s`-type averaging of independent plaquette signs, not a
polymer-weight cancellation. **No bin with a systematic size trend exists: every
cluster is the spanning one.**

**Why `|D1| = 0.5 < 0.8` is *not* a `CPB-W` pass.** It is tempting to read
`|z_signed/z_abs| = 1/2 ≤ 0.8` as satisfying `CPB-W`'s `κ_c ≲ 0.80`. It is not:
(i) it is `|cos 120°|`, present identically for every cluster of every size, so
it carries no scale-dependent suppression; (ii) it is the phase a vortex imparts
to a **fundamental loop**, i.e. it acts in the charged (N-ality-1) sector; a
`CPB-W` block weight is a neutral object and sees factor `1`, not `−1/2`. So D1
is a charged-sector constant, not the neutral-sector `κ_c`.

**Conclusion for §3.** No size-dependent dynamical cancellation exists in either
D1 or D2 — D1 is a size-independent constant `−1/2`, D2 is a size-independent
constant `≈0`. The one quantity that *does* decay with size is the
**charged-sector loop expectation** `Re⟨ω^{N(γ)}⟩` — i.e. the fundamental area
law itself (§4). That is the task's "signal only in the charged sector" closing
branch.

---

## 4. Sector test — charged vs neutral

Averaged over snapshots, β = 6, L = 6 (30 snapshots) and L = 8 (20):

| observable | full | centre-projected | verdict |
|---|---|---|---|
| `W_F(1,1)` fundamental | 0.594 | 0.548 | reproduced (centre dominance) |
| `W_F(2,2)` fundamental | 0.191 | 0.136 | reproduced |
| `W_F(3,3)` fundamental | 0.049 | 0.018 | reproduced |
| `W_A(1,1)` adjoint | 0.316 | **1** | **not** reproduced |
| `W_A(3,3)` adjoint | 0.0008 | **1** | **not** reproduced |
| `max\|W_A(proj) − 1\|` (all R,T,β,L) | — | **0.0** | exactly trivial |
| fundamental Polyakov `C_F(r)` | 0.114, 0.0127, 0.0076, 0.0065 | 1, 0.051, 0.012, 0.005 | charged disorder (decays) |
| adjoint Polyakov `C_A` (r≥1) | 1.5e-5, −1.5e-4, −1.9e-4 | **~1e-16 (≡0)** | neutral: no signal |

Creutz ratios `a²σ` (β=6, L=6; L=8 agrees): fundamental full `χ_F(2,2)=0.261`,
projected `0.381` (centre dominance, over-estimates string tension by ~1.5×);
adjoint full `χ_A(2,2)=0.589 ≈ (9/4)·0.261` — **exact Casimir scaling**, i.e. the
adjoint area law at short range is carried by **non-centre** configurations,
which projection discards entirely (projected adjoint = 1). At β = 2.6–3.3 the
lattice is too coarse for a stable Creutz ratio (`W(3,3)` unresolvable /
negative), so the clean sector separation is quoted at β = 6, where projection
is a meaningful operation; the exact adjoint-blindness statement holds at every
β regardless.

**This is the decisive test.** The vortex sector reproduces the
N-ality-1 (charged) area law and is identically trivial for N-ality 0
(neutral). Every gauge-invariant local operator is a colour singlet, hence
N-ality 0, hence sees factor `1` from any centre vortex.

---

## 5. Verdict against the task's criteria

| criterion | result |
|---|---|
| **hopeful**: `z_signed/z_abs` significantly `<1` **and** decays with cluster size | **no** — D1 is the constant `−1/2`, D2 the constant `≈0`; nothing decays with cluster size; ensemble percolates (single cluster) |
| **closed**: ratio `≈1` | **yes in the neutral sector** — projected adjoint = 1 exactly; no neutral sign channel exists |
| **closed**: signal only in the charged sector | **yes** — centre dominance in fundamental, exact triviality in adjoint |
| `κ_c ≤ 0.8` for `CPB-W` (neutral) | **impossible** from centre vortices (exact neutral blindness) |
| `κ_c ≲ 0.2` for `CPB-S` (charged) | only as the **area-law decay** of charged loops — a restatement of confinement, not an independent additive-defect charge; no new mechanism |

**Bottom line.** NP-7 does not supply a `κ_c < 1` usable by `CPB-W`. The
vortex is the confinement mechanism of the charged sector — already known — and
is exactly blind to the neutral sector where the mass gap lives.

---

## 6. What this proves, and what it does not

**Proves / demonstrates**
* Exact: centre projection is trivial in the adjoint (numerically `1 ± 0`),
  and in any N-ality-0 representation — so no centre vortex can sign-flip a
  neutral operator.
* Numerical: MCG centre projection reproduces the fundamental area law
  (centre dominance) and fails entirely in the adjoint at β=6.
* Measured: vortex ensemble percolates at these couplings; signed/absolute
  activity has no size-dependent dynamical cancellation.

**Does not prove**
* A direct measurement of the effective-action `κ_c` (no numerical route;
  STATUS §9). We exhibit the *only* sign channel (Z(3) phase) and show it is
  neutral-blind; we do not claim a numeric `κ_c`.
* Any statement about a purely dynamical vortex contribution to the neutral
  mass gap outside the `κ_c` mechanism.
* The strong-coupling (β=2.6–3.3) numbers are coarse/asymptotic-only; the
  clean sector separation is at β=6. The exact adjoint-blindness argument does
  not depend on lattice spacing.

**Gaps / honesty flags**
* Snapshot errors are statistical over snapshots; autocorrelation not
  corrected (snapshots spaced 20 sweeps) — a lower bound on the true error.
* L=8 β=2.6/3.0 configurations generated in this probe (not in `runs/`); L=8
  β=6.0 has only 20 snapshots and L=10 only 5.
* Gauge fixing is a local maximum (two restarts), not a global optimum.

---

## 7. Files

`mcg.py` (MCG, projection, clusters, loops, correlators) ·
`run_probe.py` (driver) · `analyze.py` (tables/plots) ·
`gen_L8_warm.py` (L8 β=2.6/3.0 configs) · `probe_*.npz` · `plots/`.

Reproduce (Python 3.12 with numpy/scipy/numba/networkx/matplotlib):
```
python run_probe.py --npz <cfg.npz> --label <label> --rmax 3
python analyze.py
```
