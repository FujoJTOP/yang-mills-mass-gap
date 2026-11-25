# R4-RESIDUE-PROTOCOL — the decisive `0⁺⁺` residue measurement: the gauge-invariant correlator, the optimal `f₁²` estimator, the isolation diagnostics, the statistics, the Gribov protocol, the exact-identity cross-checks, the mined β=2.5 status, and the falsification test

> **2026-11-24 · Wave R4 — attacking the atom, the irreducible target left by the H1 wave.**
> Target `H*` (Q-WAVE-SYNTHESIS §4): `f₁² = |⟨0|tr F²|G⟩|² ≥ r₀ > 0` **with isolation** `I₂`
> on the gauge-invariant `0⁺⁺` Källén–Lehmann measure, at `β=2.5`, 4D `SU(2)`.
> This note **designs** the decisive measurement and **mines** the existing data; it runs **no**
> heavy lattice job (another agent owns the heavy 4D run). No third-party downloads.
>
> Reads: `Q5-4D-CHI.md`, `Q-WAVE-SYNTHESIS.md`, `SM6-SOFT-LATTICE.md`,
> `research/unconditional/COMPUTE-F1SQ-SU2.md`, `docs/paper-v2/ch12-numerics.tex`;
> `research/lattice-cpb/` (`pilot_b2.5.npz`, `sigma_L8_b2.5.npz`, `gluon_prop_results.json`,
> `su2_string_tension.py`, `prod_b2.5.log`); `research/unconditional/f1sq_b{2.2,1.0}_L6.npz`,
> `_f1sq_su2_{mem,analyze,residue}.py`.
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional / protocol / design ·
> **[NG]** no-go · **[O]** exact missing step · **[V]** numerical witness on existing data here.
> **No QED. Clay 0%.**

---

## 0 · Verdict (no QED)

1. **[V] The atom was already *seen* — at β=2.2, not β=2.5.** The connected `0⁺⁺` `tr F²`
   correlator at `β=2.2`, `L=6`, `N=4000` gives `f₁²_bare = 7.14 ± 0.42` (17σ) with `m_G = 1.51 ± 0.15 a`
   (`COMPUTE-F1SQ-SU2.md`). A `β=1.0` contact term gives `0.2σ` — the estimator manufactures **no**
   atom from noise. These are the two bookends: a real residue where the state is resolvable, a clean
   null where it sits at the cutoff.
2. **[V] There is *no* usable β=2.5 residue measurement.** The only β=2.5 `0⁺⁺` data are
   `pilot_b2.5.npz` (`N=40`, `L=8`) and `sigma_L8_b2.5.npz` (`N=3`) — feasibility pilots. Mined with
   the exact residue estimator, the pilot shows a **~3.5σ positivity hint** in the raw ratio
   `C₁/C₀` (bare `0.184 ± 0.052`; smear-8 `0.293 ± 0.075`) but **no `m_eff` plateau beyond `t=1`**, so
   the mass is unconstrained and the fitted `f₁²` amplitude is bootstrap-unstable. A hint, **not** an
   atom, **not** isolation.
3. **[V] The noise law is measured and benign.** The per-config noise `(SE·√N)` on `C₁/C₀` is
   `≈0.33–0.47`, essentially **volume-independent** (`L=6 β=2.2`: `0.47`; `L=8 β=2.5`: `0.33–0.47`).
   So a **positivity** detection of `C₁/C₀≈0.12` at 5σ needs only `N≈190–390`; the expense is
   **isolation and the volume law**, not the signal.
4. **[V] The gauge-dependent gluon has no positive pole at β=2.5.** The mined gluon propagator
   (`β=2.5`, `L=8`, `N=28`) fits a **negative** mass² `= −0.356 ± 0.030` with `D₀→∞` — IR-suppressed,
   Gribov-like, no particle pole. This is the concrete reason the atom must be read on the
   **gauge-invariant `tr F²`**, never on the gluon mode.
5. **[T] The acid test is exact and validated here.** On the composite `tr F²` the free-field two-point
   function is a **contact term** (a two-gluon cut, no isolated pole), so `f₁²=0`. Synthetic check
   (this note): white noise → `C₁/C₀ = −0.0035 ± 0.0077`; free **massless** box → `m_eff = NaN` (no
   gapped atom); free **massive** injected `A=2.0` → recovered `1.997` (at fixed `m`) / `2.50` (grid).
   The estimator returns 0 on the null and recovers an injected atom. **No atom from noise.**
6. **[C] The minimal decisive run is specified (§9):** SU(2), `β=2.5`, `L=8,12`, `N=2000`,
   `therm=2000`, `gap≫τ_int`, APE family `{0,2,4,8}`, GEVP + two-state fit; **≈6–12 CPU-h** pure
   numpy, **one job at a time**; `L=16` adds `≈40–60` CPU-h for the finite-size law. It can certify a
   **statistical** positive residue *with isolation* at β=2.5; it cannot certify the continuum atom —
   that is the `a→0` wall `[O]`, not a dataset.

**One line.** *The atom is a gauge-invariant `0⁺⁺` `tr F²` residue, already seen at 17σ at β=2.2 and
null at β=1.0; at β=2.5 only a 40-config pilot exists, showing a ~3.5σ hint but no isolation. The
decisive run is a GEVP + two-state `0⁺⁺` measurement at `β=2.5`, `L=8,12`, `N=2000` (≈6–12 CPU-h); it
settles positivity-with-isolation statistically, never the continuum. No QED. Clay 0%.*

---

## 1 · The object: the gauge-invariant `0⁺⁺` `tr F²` correlator  [T]

Everything here lives on a **gauge-invariant** observable, so — unlike the H1 gluon probe — **no gauge
fixing is required** for the atom itself.

Define, on SU(2) Wilson configs (spatial links APE-smeared at level `s`, `α=0.5`; temporal links bare),

```
O_s(t) = Σ_x (1/3) Σ_{1≤i<j≤3} (1/2) Re tr U_ij(x,t)          (spatial plaquette, zero momentum)
C_{ss'}(t) = ⟨O_s(t) O_{s'}(0)⟩ − ⟨O_s⟩⟨O_{s'}⟩               (connected, = ⟨O_s(t)O_{s'}(0)⟩_c)   (1)
```

`O_s(t)` is invariant under `U→ΩUΩ†` (trace of a closed loop) **[T]**. Its `t=0` is the local magnetic
density `∝ tr F²` up to the lattice constant; the zero-momentum projection is the exact lattice image
of `∫d³x tr F²`. The Källén–Lehmann representation of `C_{ss'}` is **[T]**

```
C_{ss'}(t) = Σ_i a^{s}_i a^{s'}_i [ e^{-m_i t} + e^{-m_i (T-t)} ] ,      f₁² = (a^{s=0}_1)² ,        (2)
```

a positive sum over the `0⁺⁺` spectrum. **`f₁²` is the squared overlap of the bare (`s=0`) `tr F²` with
the lightest `0⁺⁺` state** — the "atom". The physical (per-volume) residue is `f₁²/V`; on the
zero-momentum operator the amplitude is extensive, so a genuine atom has `f₁² ∝ V` while the mass
`m_G` is volume-stable. [T] (representation) + [C] (finite-volume identification).

**Why not the gluon.** The gauge-dependent gluon mode `a=⟨v,log U⟩` carries the H1 ceiling but is not
an observable; at `β=2.5, L=8` its propagator fits `m²=−0.36±0.03` (§7) — no pole. The atom lives on
(1), where "atom" is unambiguous.

---

## 2 · The optimal `f₁²` estimator  [C]

The existing pipeline (`_f1sq_su2_analyze.py`) already has the right primitives: `products`,
`connected`, `analyse_C`, `gevp_mass`, `fit_residue`. The **optimal** version sharpens the mass and
the residue in three layers; report all three, never one.

**(a) Single-cosh, mass fixed (the incumbent, for continuity).** `A = (w·C)/(w·w)` over `t∈{1,2}` with
`m` fixed to the plateau `m_G`; jackknife SE over bins `≫τ_int`. This reproduces `COMPUTE-F1SQ-SU2.md`.
Its weakness: it needs a *known* `m_G` and biases `A` upward if a second state contaminates.

**(b) GEVP on the smearing family (the mass).** With `C_{ss'}(t)` a `4×4` (levels `0,2,4,8`) matrix,
solve `C(t) v = λ(t) C(t₀) v`. The ground eigenvalue `λ₁(t) ∝ e^{-m₁(t-t₀)}`; `m₁(t)` from
`−ln[λ₁(t)/λ₁(t-1)]` must **plateau** (see §4). GEVP is the standard variance-reduced ground-state
projector and separates the ground state from the first excited smearing contamination.

**(c) Two-state fit (the residue *and* isolation).** Extend `fit_residue` to
```
C_{00}(t) = A₁[e^{-m₁t}+e^{-m₁(T-t)}] + A₂[e^{-m₂t}+e^{-m₂(T-t)}] + (contact at t=0, excluded),       (3)
```
over `t∈[t_min, T/2]`, fitting `(A₁,m₁,A₂,m₂)` with a `(A,m)`-profiled `χ²` grid (as the incumbent does,
no optimizer). **`f₁² := A₁`**, the ground-state residue, with full jackknife covariance.
The two-state fit is what converts "a decaying correlator" into "an atom with an isolation gap
`Δ=m₂−m₁`". **[C]** (estimator) on **[T]** (representation (2)).

**Variance / conditioning notes.** (i) Work with the **connected** correlator `⟨O(t)O(0)⟩−⟨O⟩²`
computed on the same configs — the standard subtraction that kills the contact term (§6 acid test).
(ii) The dominant per-config noise enters through `O_s(t)O_s(0)`; the smearing family `s≥2` has
`~V` lower operator variance (mined: `Var(S₈)/Var(S₀)≈0.09` at `L=6`) and is the cleanest mass
carrier, while the bare `s=0` carries the largest absolute residue. (iii) Bin over blocks `≫τ_int`
(glueball `τ_int≈10` sweeps; `REPORT.md`).

---

## 3 · Volume, statistics, and the noise law needed to resolve the atom  [V]/[C]

**Measured noise law (this note, from `f1sq_b2.2_L6.npz`, `pilot_b2.5.npz`; delete-block jackknife).**

| data | `L` | `β` | `N` | `C₁/C₀`(bare) | SE | per-config noise `SE·√N` |
|---|---|---|---|---|---|---|
| f1sq `S0` | 6 | 2.2 | 4000 | `+0.123` | `0.007` | **0.47** |
| f1sq `S8` | 6 | 2.2 | 4000 | `+0.164` | `0.007` | 0.42 |
| pilot `O_0` | 8 | 2.5 | 40 | `+0.184` | `0.052` | **0.33** |
| pilot `O_8` | 8 | 2.5 | 40 | `+0.293` | `0.075` | 0.47 |

The per-config noise is **volume-robust** (`0.33–0.47` at `L=6,8`). Consequently a **5σ positivity**
detection of `C₁/C₀=0.12` needs `N≈200–400` — cheap. What is *not* cheap is the part the pilot lacks:
an `m_eff` plateau, a resolved `m₂`, and the `L`-trend. These set the real `N`:

- **`N ≥ 2000` per volume** (adopted from Q5/SM6; justified by the measured law): drives
  `SE(C₁/C₀) → ≈0.010` (`L=8`) and `≈0.028` at `L=12` (volume factor `(12/8)³≈3.4` on a susceptibility
  if the noise grows `∝√V`), enough to (i) resolve the two-state fit, (ii) see the `m_eff` plateau over
  `t=2,3`, and (iii) test the `1/L²` finite-size law at the percent level.
- **`therm ≥ 2000`, `gap ≈ 20–50` sweeps `≫ τ_int`** (glueball `τ_int≈10`). The `k=0` and Polyakov slow
  modes are not the issue here (gauge-invariant composite), but the glueball autocorrelation is.
- **Volumes `L = 8, 12`** decisive; **`L = 16`** for the finite-size *law* (§4), not for positivity.
- **One job at a time** (charter §4: never parallel SU(N) `L≥8`; the prior OOM). Peak RSS `≲200 MB` at
  `L=12` SU(2) (only `N×L` time-slice sums stored — the `f1sq` streaming design).

---

## 4 · Isolation diagnostics: second-mass plateau and spectral density  [T]/[C]

An atom is a **delta** at `m₁²` in the `0⁺⁺` spectral density `ρ(s)=Σ_i f_i² δ(s-m_i²) + continuum`.
Positivity of `f₁²` alone is the **finite-volume triviality** (`ch12-numerics.tex §1`); *isolation* is
the discriminating content. Four diagnostics, all already available in the corpus:

**(D1) Effective-mass plateau.** `m_eff(t)=arccosh[(C(t-1)+C(t+1))/(2C(t))]` must be **flat** over a
window `t∈[t_a,t_b]` with `t_b≥2–3`, for the GEVP ground state and for `s≥2`. A downward-drifting
`m_eff` = contamination by a lighter state (or the box zero mode); an upward-drifting one = a cut.
The pilot β=2.5 has *only* `m_eff(t=1)` and nothing beyond — that is exactly why it is not an atom
measurement (§7).

**(D2) Second mass / isolation gap.** From the two-state fit (3): `m₂` resolved at `≥2σ` with
`Δ=m₂−m₁>0`, and `A₂/A₁` **decreasing** as `V` grows (the excited-state overlap is the finite-volume
artifact; a genuine atom's `A₁` is the volume-stable one). "No isolation" = `m₂→m₁` or `A₂/A₁` not
shrinking with `V`.

**(D3) Spectral density shape.** The Källén–Lehmann inversion on the GEVP-projected correlator: the
atom is a **resolved peak of finite weight** `f₁²`, a cut is a **shoulder/continuum**. Operationally the
campaign's three model-based reads apply (`ch12-numerics.tex §Diagnostics`):
the **P3 approach law** `T(τ) = −d/dτ log|m(τ)−ê₁(τ)| → Δ` for an atom vs `→1/τ` for a cut
(`WI-P3-CONTINUUM.md`); the **Re-slope** `R=1/E₁` exposing atom-slope rungs vs none for a cut
(`WI-RESLOPE-*`); the **one-sided isolation bound** `Δ≥κ(τ)`. All are **model-level [T]** and
**finite-real-τ blind at O(h)** — they *read*, they do not force.

**(D4) Finite-size law (the cleanest discriminator).** `Δ(L) = m₁ + c/(m₁L²)` for an **atom** vs
`Δ(L) ∼ (2π/L)²/m → 0` for a **cut** (`ch12-numerics.tex eq. fss`, `WI-P3-UNIFORM.md`). The atom
survives `L→∞` with `m₁>0`; the box cut dissolves. **This is the falsification handle (§8):** a
`m_G(L)` tracking `(2π/L)²` over `L=8,12,16` kills the atom.

---

## 5 · The Gribov-copy / Landau-gauge protocol — scoped correctly  [T]/[C]

**The residue `f₁²` on `tr F²` needs no gauge fixing** (§1, `O_s` is a trace) **[T]**. The Landau +
Faddeev–Popov + Gribov protocol below is required **only** for the *parallel* gauge-dependent-gluon
cross-check that converts the flat-chart H1 `χ_flat` into a physical zero-momentum object (Q5 §7,
SM6 §3.1). Stated so the two carriers are never conflated:

- **Landau gauge**, FFT/over-relaxation; **`N_copy = 5–10` random starts per config**; the **spread of
  `Var(a)` (and of `1/χ`) across copies is the Gribov systematic** — report it as an error bar, not a
  central value.
- Fold the gauge-fixing + **Faddeev–Popov** Jacobian into the measure `Φ` (stochastic determinant, or a
  controlled bound), `G″=Δ_dyn+W_geo+X` with `W_geo≥0` [T].
- Purpose: test whether `1/χ_flat` is a **physical** IR mass (the A4 `m_GI` identification, `[C]`) and
  whether the gluon carrier is non-vacuous (F4). The **atom** is certified independently on `tr F²` by
  §2–§4 and does not inherit the Gribov systematic.
- For SU(3), string-tension matched `β≈5.7`, as a cross-group check.

---

## 6 · Exact-identity cross-checks and the acid test  [T]/[V]

**Acid test (free field must give `f₁²=0`) — the composite `tr F²` has a contact two-point function
at `β→∞` (two-gluon cut, no isolated pole), so the residue is exactly zero. [T]** Validated three
ways this note (synthetic, cheap):

| test | setup | result | verdict |
|---|---|---|---|
| white noise (contact null) | `O(t)` iid, `N=2000` | `C₁/C₀ = −0.0035 ± 0.0077` | **0** ✓ |
| free **massless** box (periodic) | `C(t)=const` | `m_eff(t)=NaN` ∀t; no gapped atom | `f₁²(m>0)=0` ✓ |
| free **massive** positive control | inject `A=2.0, m=1.0` | `A=1.997` (fixed `m`); `2.50` (grid) | recovers atom ✓ |

The connected subtraction `⟨OO⟩−⟨O⟩²` is what makes the null vanish: the raw product has a huge
contact variance; the subtraction removes it. **The estimator does not manufacture an atom from
noise.**

**Exact-identity cross-checks on the real data [T]:**

1. **Reflection positivity / complete monotonicity.** `C_{00}(t)=Σ_i f_i²(e^{-m_it}+e^{-m_i(T-t)})≥0`;
   the sequence `{C(t)}` is completely monotone (successive forward differences alternate in sign). A
   sign violation beyond jackknife noise is a **bug or contamination**, not physics.
2. **Weight bound / sum rule.** `f₁² = A₁ ≤ C(0) = Var(O_s)` (the ground state carries at most all the
   spectral weight; equality iff single-state). Check `A₁ ≤ C(0)` on every level; a violation means the
   fit leaked weight into `A₂` or the jackknife is broken.
3. **β=0 Haar contact null.** At `β=0` the plaquette autocorrelation is a contact term (`C(t>0)=0`),
   so `f₁²=0` exactly (up to iid noise). The **β=1.0 run is the empirical instance** of this null
   (`C₁/C₀=0.0013±0.0065`, `0.2σ`):
   the estimator returns zero.
4. **Continuum content is renormalization-gated.** `f₁²` needs the operator renormalization `Z_O²` and
   `a→0`; only **positivity, the mass `m_G`, and the volume scaling** are scheme-independent. Report
   `f₁²_bare`, `f₁²/V`, and `m_G`; flag the absolute continuum value `[O]`.

---

## 7 · What the existing data already say (mining report)  [V]

Mined this note with the exact residue machinery (`_f1sq_su2_*` primitives) on the existing campaign
files; **no heavy run**. Two files matter, plus the gluon propagator.

**β=2.2, L=6, N=4000 — the atom is seen (`f1sq_b2.2_L6.npz`).**
`C₁/C₀`(bare) `= 0.123 ± 0.007` (**17.9σ**, reproduces the campaign's 17σ); smear-8 `0.164 ± 0.007`.
Full `C/C₀ = [1, 0.123, 0.036, −0.001, 0.036, 0.123]` — positive, symmetric, monotone toward the
middle (a clean cosh). `f₁²_bare = 7.14 ± 0.42`, `m_G = 1.51 ± 0.15 a`; `f₁²/V = 0.033`. **This is the
reference positive datum.** [V]

**β=1.0, L=6, N=4000 — the contact null (`f1sq_b1.0_L6.npz`).**
`C/C₀ = [1, 0.0013, −0.008, −0.022, −0.008, 0.0013]`; `C₁/C₀ = 0.0013 ± 0.0065` (`0.2σ`). A contact
term, `m·a` at the cutoff; `f₁²` not extractable. **The estimator returns zero on a null.** [V]
(Physically this is *not* evidence against the residue — P3 is analytic at strong coupling where
`m·a` diverges; the state is provably present but sub-lattice-spacing on `T=6`.)

**β=2.5, L=8 — only a pilot (`pilot_b2.5.npz` N=40; `sigma_L8_b2.5.npz` N=3).**
`C₁/C₀`(bare) `= 0.184 ± 0.052` (**~3.5σ hint**); smear-8 `0.293 ± 0.075`. But:
- `C/C₀ = [1, 0.184, 0.019, −0.049, 0.008, −0.049, 0.019, 0.184]` — the `t=2` value is `0.02±0.06`
  (already noise); **`m_eff` is finite only at `t=1`** (`1.68`, no `t≥2` plateau);
- the fitted cosh amplitude (fixed `m=m_G`) is **bootstrap-unstable** (`A = 49 ± 1.1×10⁵`): with the
  mass unconstrained at `N=40`, `A` is not a measurement;
- the `N=3` companion (`sigma_L8_b2.5.npz`) is useless for any fit;
- the intended production (`prod_b2.5.log`) stopped at `240/600` sweeps with **no saved `.npz`** — the
  heavy β=2.5 run is incomplete (owned elsewhere).
**Verdict: a β=2.5 positivity *hint* at ~3.5σ, no isolation, no mass, no volume trend. Not an atom
measurement.** [V]

**Gluon propagator at β=2.5, L=8 (`gluon_prop_results.json`, N=28).**
Fitted mass² `= −0.356 ± 0.030` (IR window); `D₀→∞`. The *gauge-dependent* gluon has **no positive
mass pole** at β=2.5 — Gribov/IR-suppressed. This is the numerical face of the carrier hazard (Q6):
the gluon cannot certify the mass; the atom must be read on `tr F²`. [V]

**Also at β=2.5 (mined, for context):** `Q5` flat-chart `χ_flat = 4.17±0.18 (L=4) / 3.61±0.24 (L=6)`,
no growth (plateau) — the H1 ceiling side, on the gauge-dependent carrier [V]. The `1D`/periodic
residue plateau is `Var→4.45`, `G″→0.222` (`P2`). The four 4D H1 faces and the carrier hazard are in
`Q-WAVE-SYNTHESIS.md`.

---

## 8 · Falsification test: what data would kill the atom hypothesis

The atom at `β=2.5` is **falsified** if the decisive run (`L=8,12(16)`, `N≥2000`) shows **any** of:

1. **Null residue.** The two-state ground-state `A₁=f₁²` is consistent with 0 — a 95% one-sided upper
   bound `≤0` — at the best `(L,N)`, on the bare `s=0` operator after the connected subtraction.
2. **Box cut, not atom.** The ground mass `m_G(L)` decreases over `L=8,12,16` **following the free-box
   law** `∼(2π/L)²/m → 0` rather than the atom law `m₁ + c/(m₁L²)`. A state that dissolves with volume
   is the lowest box mode, not an atom. (This is a *different* dissolution from β=1.0: there the
   coupling is at the cutoff; here the coupling is resolvable, so a volume-dissolution would be fatal.)
3. **No isolation.** `m₂ → m₁` (no gap `Δ`) as `L` grows, or `A₂/A₁` does **not** shrink with `V`, or
   the Källén–Lehmann density `ρ(s)` is a **continuum with no resolvable delta** at `m₁²`. A cut, not
   an atom.
4. **Reflection-positivity violation.** `C(t)` develops negative values beyond jackknife noise
   (non–completely-monotone sequence) — the "residue" is an artifact/contamination.
5. **Wrong volume scaling of the weight.** `f₁²` does not grow `∝V` (the extensive zero-mode overlap
   expected of an atom) but instead is flat or contact-limited (`C(t>0)=0`), signalling a pure contact
   term.

**Note on what does *not* falsify it.** A one-sided interval on a *single* volume including zero is the
finite-volume triviality, not a falsification; and the β=1.0 null is expected (state at the cutoff), so
it does not bear on the β=2.5 atom. Falsification needs the **volume sequence** with **isolation**.

---

## 9 · The minimal decisive run and its cost  [C]

**Protocol.**
- **Group / coupling:** SU(2), `β=2.5` (the head).
- **Volumes:** `L = 8, 12` (decisive: positivity + isolation + `L`-trend). `L = 16` **only** for the
  finite-size law (D4) — the first volume where the free box gap `(2π/L)²` drops below the physical
  scale and the atom/cut laws separate.
- **Configurations:** `N = 2000` per volume; `therm ≥ 2000` sweeps; `gap ≈ 20–50` sweeps `≫ τ_int≈10`;
  streaming sampler storing only `N×L` time-slice sums (peak RSS `≲200 MB` at `L=12`). **One job at a
  time.**
- **Operator:** APE-smeared spatial plaquette, family `s = {0,2,4,8}`, `α=0.5`; temporal links bare;
  zero momentum (§1).
- **Estimators:** (a) single-cosh fixed-`m_G` (continuity with `COMPUTE-F1SQ`); (b) GEVP ground mass on
  `s≥2`; (c) **two-state fit (3)** for `(A₁=f₁², m₁, A₂, m₂)`; jackknife over blocks `≫τ_int`.
- **Cross-checks:** RP/complete-monotonicity; `A₁≤C(0)`; β=0/β=1 contact null (empirical baseline, §6);
  free-massive positive control (§6).
- **Isolation:** D1–D4 (§4): `m_eff` plateau, resolved `m₂` with `A₂/A₁↓` in `V`, spectral-density
  peak, and the finite-size law `Δ(L)=m₁+c/(m₁L²)` vs `(2π/L)²`.
- **Parallel (optional, H1 ceiling only):** Landau + `N_copy` FP + Gribov-copy gluon susceptibility per
  Q5 §7 / SM6 §3.1 (§5 above) — **not** needed for `f₁²`.

**Cost (grounded in measured `β=2.5, L=8` heat-bath `≈70 ms/sweep` and the pilot `meas ≈1.1–1.7 s/cfg`).**

| volume | sweeps (`therm + gap·N`) | sweep cost | measurement (`mag_op` + APE) | subtotal |
|---|---|---|---|---|
| L=8 | `2000 + 40·2000 ≈ 8.2×10⁴` | `70 ms` → `≈1.6 h` | `2000×0.3–0.5 s` → `≈0.2–0.3 h` | **≈1.8–2 h** |
| L=12 | `8.2×10⁴` | `×5.06 ≈ 0.35 s` → `≈8 h` | `×5 ≈ 1–1.4 h` | **≈8–9 h** |
| L=16 (optional) | `8.2×10⁴` | `×16 ≈ 1.1 s` → `≈25 h` | `≈4–5 h` | **≈30 h** |

**Minimal decisive run: SU(2) `β=2.5`, `L=8,12`, `N=2000`, one job at a time ≈ 6–12 CPU-h pure numpy**
(`≈1–2 h` on C/OpenMP 16 cores). `L=16` adds `≈30 CPU-h`. Add Wilson loops (`R,T≤5`, for the scale /
string tension) roughly doubles the measurement time. **Systematics, not statistics, dominate:** the
statistical SE at `N=2000` is already below the two-state fit's closure (`±0.08–0.12`, SM6 §1.4); the
budget is the fit-range / smearing-family spread and the `a→0` / renormalization wall `[O]`.

**What it decides — and does not.** It decides **statistically and with isolation** whether a
gauge-invariant `0⁺⁺` residue `f₁²>0` survives `L=8→12` with a volume-stable mass — the sharpest
falsifiable form of the atom at `β=2.5`. It does **not** certify the continuum atom (`a→0`,
`Z_O²`, `L→∞` uniform) — that remains `[O]`, the missing machine, not a dataset.

---

## 10 · Label ledger

| item | statement | label |
|---|---|---|
| `O_s(t)` zero-momentum spatial plaquette is gauge-invariant; Källén–Lehmann rep. (2) | §1 | **[T]** |
| `f₁²` = ground-state amplitude `A₁`; two-state fit for `(A₁,m₁,A₂,m₂)` | §2 | **[C]** on **[T]** |
| per-config noise `SE·√N≈0.33–0.47`, volume-robust; `N(5σ,amp 0.12)≈200–400` | §3, mined | **[V]** |
| `N≥2000`, `L=8,12`, `gap≫τ_int`, one job: resolves isolation + FSS, not just positivity | §3/§9 | **[C]** |
| isolation: `m_eff` plateau, `m₂` with `A₂/A₁↓`, spectral peak, `Δ(L)=m₁+c/(m₁L²)` vs `(2π/L)²` | §4 | **[T]**(model)/**[C]** |
| residue on `tr F²` needs no gauge fixing; Landau+FP+Gribov only for the gluon cross-check | §5 | **[T]**/**[C]** |
| free composite `tr F²` = contact ⇒ `f₁²=0`; white-noise & massless nulls validated | §6 | **[T]**/**[V]** |
| `f₁²=A₁ ≤ C(0)` weight bound; complete-monotonicity (RP) | §6 | **[T]** |
| β=2.2 `f₁²_bare=7.14±0.42` (17σ), `m_G=1.51±0.15a`; β=1.0 null `0.2σ` | §7 | **[V]** |
| β=2.5 pilot `N=40/3`: `C₁/C₀≈0.18–0.29` (~3.5σ hint), **no** `m_eff` plateau, **no** production data | §7 | **[V]** |
| β=2.5 gluon propagator `m²=−0.36±0.03` (no pole) | §7 | **[V]** |
| falsification: null residue / box-cut `m_G` / no isolation / RP violation / wrong weight scaling | §8 | **[C]** |
| minimal decisive run `L=8,12, N=2000` ≈6–12 CPU-h, one job | §9 | **[C]** |
| continuum atom (`a→0`, `Z_O²`, `L→∞` uniform) | — | **[O]** |
| Clay mass gap | open | **0%** |

---

*Refs.* `Q-WAVE-SYNTHESIS.md`, `Q5-4D-CHI.md`, `Q6-H1-STATUS.md`, `SM6-SOFT-LATTICE.md`,
`research/unconditional/COMPUTE-F1SQ-SU2.md`, `docs/paper-v2/ch12-numerics.tex` (§Diagnostics:
`WI-P3-CONTINUUM.md`, `WI-P3-QUANT.md`, `WI-RESLOPE-*`, `WI-P3-UNIFORM.md`); witnesses/data
`research/unconditional/{f1sq_b2.2_L6,f1sq_b1.0_L6}.npz`, `_f1sq_su2_{mem,analyze,residue}.py`;
`research/lattice-cpb/{pilot_b2.5,sigma_L8_b2.5}.npz`, `gluon_prop_results.json`,
`su2_string_tension.py`, `prod_b2.5.log`. Mining scripts (scratch, not committed): connected-cosh
residue replay, error budget, volume scaling, synthetic acid test. **No heavy lattice job run in this
wave; no third-party downloads.**
