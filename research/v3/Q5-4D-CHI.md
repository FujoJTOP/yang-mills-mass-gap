# Q5-4D-CHI — a controlled 4D SU(2) test of the H1 plateau: the k=0 susceptibility **does not grow** (a free field would rise `≍L²`), it is **consistent with a plateau at the Haar value** `χ≃3.6–4.2`; but the object is a flat-chart coordinate (gauge-invariant only infinitesimally, `≤4%` under moderate local transforms), the exact force-moment `G″` is unresolved at pilot statistics (`+0.03±0.08`), and the rigorous `L`-uniform bound is still `H1`

> **2026-11-24 · Wave Q5 — the last open input, H1, attacked numerically.**
> H1: the `k=0` susceptibility `χ = Var(a)` is bounded uniformly in `L` at `β=2.5`, equivalently
> `G″_{k=0} ≥ λ > 0`. In 1D periodic SU(2) this is a plateau (`Var(a)→4.45`, `G″→0.222`,
> transfer gap `Δ=0.6444`). In 4D it is open. This note designs a **gauge-controlled** estimator of
> the `k=0` amplitude, runs a **pilot** at `L=4,6`, and states the minimal decisive run.
> Witness `research/v3/w7/q5_4d_chi.py` (`py -3.12 -I …`, **one job**, SU(2) heat bath; result JSON
> `q5_4d_chi_results.json`). **No third-party sources.**
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional / protocol · **[NG]** no-go
> for a named class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

1. **[T] The free-field acid test passes.** An exact free **massive** field gives `Var(a)=1/m²`
   (plateau, `L`-flat) and `G″=m²→0`; a free **massless** field (box-regulated) gives
   `Var/L² = 0.026` const (`∝L²`) and `G″·L² = 38.7` const (`G″∝1/L²→0`). The tool manufactures **no
   floor** for the free field.
2. **[V] The pilot data do *not* grow.** At `β=2.5`: `χ_flat = 4.17±0.18 (L=4)`, `3.61±0.24 (L=6)`,
   ratio `0.87` — the free-field growth ratio would be `≈2.0`. **Growth is excluded; the data are
   consistent with a plateau** at the Haar value `s²=(4π²−6)/9=3.72` (β=0 check: `3.52±0.16`).
3. **[V] The gauge-invariant Polyakov spatial correlator `C_P(r)` decays within `ξ≲1–2` sites** at
   both `L` → bounded (plateau) susceptibility in the N-ality-1 sector. (The `V_s·Var(P̄)` estimator is
   slow-mode-contaminated at pilot `N` and is *not* used.)
4. **[V]/[C] The exact force-moment `G″` at the `k=0` mode is *not* resolved:** KDE `G″=+0.03±0.08`,
   `1/Var=0.216`, `W_geo=+0.24` (the stable carrier), `bracket=−0.39`. The marginal is broad
   (`sd(a)=2.15`), so the curvature needs far more statistics; the Gaussianized floor `1/χ≈0.2>0`
   is real, the exact `G″≥λ` is the open half.
5. **[NG] No strictly gauge-invariant estimator with `Var↛0` exists linearly** (H1-NUMERIC-RIGHT):
   composites give `Var∝1/V` (vacuous); the fundamental-field mode — the only `Var↛0` carrier — is a
   flat-chart coordinate, gauge-invariant only infinitesimally (measured `≤4%` drift under moderate
   finite local transforms; `95%` under order-π ones). The rigorous `L`-uniform bound is **`H1` [O]**;
   no finite dataset supplies it (§7).

**One line.** *The 4D `k=0` susceptibility does not grow like a free field; it is consistent with a
plateau near the Haar value. But the plateau lives on a gauge-dependent coordinate, the exact
curvature is statistically unresolved, and the uniform-in-`L` bound remains `H1`. No QED, Clay 0%.*

---

## 1 · The object, and where gauge invariance actually lives  [T]/[NG]

H1 names `a_{k=0}` in mode normalization (`‖v‖=1`), `a = N^{-1/2}Σ_l ξ_l`, `ξ_l = log U_l` in the flat
chart, `N=4L⁴`. Two facts fix what can and cannot be measured:

- **The composite branch is vacuous [NG/T].** Any local gauge-invariant `O(x)` has `Var(Õ(k))∝1/V`
  (extensivity/CLT): the `tr F²` zero mode has `Var(q̄)=χ_q/V`, so `G″=1/Var∝V→∞` — H1's
  `Var≤1/λ` is satisfied for every `λ` at large `V`, i.e. **vacuously** (H1-NUMERIC-RIGHT §2–3;
  the susceptibility `χ_q≈0.04` is finite *and gap-blind* — finite for a free dim-4 operator too).
- **Only the fundamental field mode has `Var↛0` — and it is gauge-dependent.** `a` transforms under
  `U_l→Ω(x)U_lΩ(x+μ)†`. The **k=0** mode is invariant to **leading order** because the gauge
  variation is a gradient whose zero-momentum part vanishes (in 1D the sum telescopes exactly — which
  is why the 1D P2 result is cleanly gauge-invariant). In 4D non-abelian the `log` is not additive, so
  `a` is analytic-only: **[V] measured** — exactly invariant under global colour rotations (`rel
  4.6×10⁻¹⁵`), drifting `0.05%→5.8%` for per-site rotation size `ε=10⁻³→10⁻¹` on one config, and on
  the **ensemble mean** `≤4%` for moderate transforms (`ε≤0.4`, §5). Order-π transforms change it by
  `O(1)`.

So a genuinely gauge-invariant statement needs either a **complete gauge fixing** (Landau + FP +
Gribov control) or a **gauge-invariant two-point object** (the glueball mass `m_GI`, a decay rate).
The pilot reports the flat-chart `χ_flat` with its measured gauge systematic, plus the gauge-invariant
Polyakov correlator as a cross-check — the honest split of H1-NUMERIC-RIGHT.

---

## 2 · The estimators  [C]

| tag | definition | role |
|---|---|---|
| **χ_flat** | `χ = (1/3)⟨|A|²⟩`, `A=N^{-1/2}Σ_l ξ_l` (flat-chart `k=0`, one colour-averaged vector). Globally colour-invariant; the corpus's `Var(a)` (D2/SM6 normalization). | primary `Var(a)`; gauge-fixed proxy |
| **G″ force-moment** | exact identity `G″(v;a₀)=E[∂_v²V|a₀]−Var(∂_vV|a₀)=bracket+W_geo+X`, `V=S_W+Φ`, kernel-conditioned at `a₀≈0` (the mode; `ρ` even). `1/Var` is its Gaussian limit. | the "exact-identity / force-moment" estimator the task names |
| **C_P(r)** | gauge-invariant spatial correlator `⟨P(x)P(x+r)⟩−⟨P⟩²`, `P=½tr∏_tU_0`. | genuinely gauge-invariant plateau / correlation-length cross-check |
| **free acid** | exact Gaussian | `Var=1/m²`, `G″=m²`; massless `Var∝L²`, `G″→0` |

The `k=0` mode is **not** a plane wave (`k=2πn/L`, `n≥1`, the stored-data basis H1-NUMERIC-RIGHT flags
as gauge-dependent and blind to the soft mode); it is the uniform `n=0` mode.

---

## 3 · Acid test — free field must fail  [T]

Exact, deterministic (witness `[FREE]`):

```
 massive m      Var=1/m^2   G''=m^2    G''*Var
   1.00           1.000      1.000      1.000
   0.30          11.111      0.090      1.000
   0.03        1111.11       0.0009     1.000     <- plateau in L, G''->0 as m->0
 massless (box, k_min=2pi/L):
   L    Var        Var/L^2     G''        G''*L^2
   4    0.500      0.0312      2.000       32.0
   8    1.707      0.0267      0.586       37.5
  16    6.569      0.0257      0.152       39.0    <- Var ~ L^2 , G'' ~ 1/L^2 -> 0
```

The free massless field grows `∝L²` with `G″→0`; the free massive field plateaus with `G″=m²`. **The
tool discriminates and manufactures no floor.** [T]

---

## 4 · `β=0` Haar validation  [V]

Pure-Haar sampler check (`[HAAR]`): `χ_flat=3.52±0.16` vs exact `s²=3.7198` (≈1.25σ; the estimator is
unbiased but has per-config `sd≈2`, §5); the heat bath at `κ=0` reproduces `⟨|ξ|²⟩/3=3.7182±0.0077`
**exactly**, so the sampler is correct. `C_P(r)=[0.248, 0.005, 0.001, 0.005]` — Polyakov correlations
vanish for `r>0` as they must. `χ_P=0.281±0.035` (`≈C_P(0)`, `L`-flat). Validates sampler + estimators.

---

## 5 · Pilot data, `β=2.5`, periodic SU(2)  [V]

SU(2) Haar heat bath, cold start, `therm=800`, `gap=8`, one job at a time (charter §4). Delete-block
jackknife SEs; `τ_int` from the measured series.

| quantity | `L=4` (`N=320`) | `L=6` (`N=200`) |
|---|---|---|
| `χ_flat = Var(a)` | **`4.17 ± 0.18`** | **`3.61 ± 0.24`** |
| `τ_int(χ_flat)`, `N_eff` | `1.02`, `314` | `1.26`, `159` |
| `⟨plaq⟩` | `0.6548` | `0.6522` |
| `C_P(0),C_P(1),C_P(2),…` | `0.274, 0.108, 0.092, 0.108` | `0.239, 0.028, 0.011, 0.006, 0.011, 0.028` |
| `G″(k=0)` KDE | `+0.034 ± 0.084` | — |
| `bracket = E[H]−Var(G)` | `−0.393` | — |
| `W_geo` (proven `≥0`) | `+0.236 ± 0.009` | — |
| `X` | `+0.067 ± 0.073` | — |
| `1/Var(a)` | `0.216` | — |
| `sd(a)` | `2.15` | — |

**Independent runs** (different seeds/thermalisation) give `χ_flat(L=4)=3.8–4.2`: the per-config
`χ_flat` has `sd≈2` (heavy-tailed, near-π links), so single-volume means carry `±5%` at `N≈300`.

**Gauge systematic** (24-config subsample, base `5.10`, so the *relative* drift is the content):
`ε=0.05→+1.6%`, `0.1→+1.3%`, `0.2→0%`, `0.4→−3.7%`. The ensemble mean is stable to `≤4%` under
moderate finite local gauge transforms; exactly stable under global colour rotations. (A genuine
gauge-invariant statement still needs Landau+FP, §7.)

**Force-moment reading.** `W_geo=+0.236` is the stable, proven-`≥0` carrier (matches SM6's `≈0.25`);
`bracket` is soft-negative; `G″` itself is `+0.03±0.08` — statistically consistent with zero because the
marginal is broad (`sd(a)=2.15`). The Gaussianized `1/χ=0.216>0`; the exact `G″≥λ` is unresolved.

---

## 6 · Plateau or growth?  [V]→[C]

- `χ_flat`: `4.17 → 3.61` (`L=4→6`), **ratio `0.87`**. A free massless field would rise by
  `(6/4)²≈2.25` (equivalently `Var/L²` const, `[FREE]`). **Growth is excluded**; the ratio is
  consistent with a plateau (and with a slight finite-size decrease toward the Haar value `3.72`).
- `C_P(r)` decays within `ξ≲1–2`: `C_P(1)/C_P(0)=0.12` at `L=6` (`0.39` at `L=4`, small-volume
  enhancement). Bounded correlations ⇒ bounded susceptibility in the N-ality-1 sector.

**Verdict:** the pilot data **plateau** (do not grow like the free field), at `χ≈3.6–4.2 ≈` the Haar
floor. [V] for plateau-consistency, [C] beyond (the object is gauge-dependent and the curvature
is statistically unresolved). This is the 4D analogue of the 1D P2/D2 plateau and is fully consistent
with — but does not prove — a mass gap.

---

## 7 · The minimal decisive run  [C]

Goal: decide **statistically, with gauge control**, whether the **gauge-fixed zero-momentum
susceptibility** `χ(k=0)` plateaus or grows as `L→∞` at `β=2.5`, and resolve the exact `G″`.

- **Volumes / coupling.** SU(2), `β=2.5`, `L=8,12` (decisive) `+ L=16` (trend). One job at a time;
  never parallel SU(N) `L≥8` (charter §4 OOM history). Peak RSS `≲1 GB/job` at `L=12` SU(2).
- **Statistics.** `N≥2000` configs, `therm≥2000`, `gap≫τ` of the slow modes (the `k=0`/Polyakov modes
  have `τ` of tens of sweeps — use `gap≈50–100`, or a global heat bath on the constant mode). This
  drives the `±5%` `χ_flat` SE at `N≈300` below `±1–2%`, enough to separate plateau (`ratio≈1`) from
  `L²`-growth (`ratio≈2.25` over `8→12`) by an order of magnitude.
- **Gauge fixing (the new ingredient).** Landau gauge, over-relaxation, `N_copy=5–10` random starts per
  config; the **spread of `χ(k=0)` across copies is the Gribov systematic**; fold the Faddeev–Popov
  determinant into `Φ` (stochastic estimate or controlled bound). This turns the flat-chart proxy into
  the standard zero-momentum gluon susceptibility.
- **Estimator.** `χ(k=0)` (Landau `Var(a)`), the force-moment `G″=bracket+W_geo+X` (kernel-conditioned),
  and `1/Var`; **plus** the two genuinely gauge-invariant checks: `C_P(r)` (Polyakov, `ξ`) and the
  `tr F²` glueball mass `m_GI` (decay rate, A4) — the only objects that carry an unconditional gap
  statement.
- **What it decides / does not.** It decides the **statistical** plateau-vs-`L²`-growth of the
  gauge-fixed `k→0` susceptibility with a controlled Gribov systematic, and resolves `G″` at the mode.
  It does **not** give a rigorous `L`-uniform bound — that is `H1` [O]; a finite dataset cannot certify
  it. Cost: pure-numpy SU(2) `L=8,12`, `≈` a few CPU-hours each, one at a time; `L=16` `≈5×` `L=12`.

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| free-field tool: massive `Var=1/m²` flat, massless `Var∝L²`, `G″→0` | §3 `[FREE]` | **[T]** |
| composite `0⁺⁺` zero mode is a susceptibility (`Var∝1/V`, `G″∝V`), vacuous & gap-blind | §1; H1-NUMERIC-RIGHT | **[NG]**/**[T]** |
| `k=0` flat-chart `a` globally colour-invariant; gauge-invariant only infinitesimally | §1; measured §5 | **[V]** |
| `β=0` Haar check `χ_flat→3.72`, `C_P(r>0)→0` | §4 | **[V]** |
| pilot `χ_flat`: `4.17±0.18 (L4)`, `3.61±0.24 (L6)`, ratio `0.87` — no growth | §5–6 | **[V]** |
| gauge-invariant `C_P(r)` decays, `ξ≲1–2` → bounded (plateau) N-ality-1 sector | §5–6 | **[V]** |
| exact force-moment `G″(k=0)` unresolved (`+0.03±0.08`); `1/χ=0.22>0`, `W_geo=+0.24` | §5 | **[V]** |
| `1/χ` as physical `k→0` gluon mass / `m_GI` identification | — | **[C]** |
| **uniform-in-`L` `Var(a)≤1/λ`** (= `H1`) | — | **[O]** |
| decisive run: Landau+FP+Gribov, `L=8,12(16)`, `N≥2000`, one job at a time | §7 | **[C]** |
| Clay mass gap | still open | **0%** |

---

## 9 · Sharpest line and reproduce

> **The 4D `k=0` susceptibility does not grow — a free field would rise `∝L²` and `G″→0`, but the
> measured `χ_flat` is `L`-flat (`4.17→3.61`, ratio `0.87`, versus the free `2.25`) and the
> gauge-invariant Polyakov correlator decays within `ξ≲1–2`. So the pilot data are consistent with a
> 4D plateau at the Haar value `≈3.7` [V]. But the plateau lives on a flat-chart coordinate
> (gauge-invariant only infinitesimally, `≤4%` drift under moderate local transforms, `O(1)` under
> large ones), the exact force-moment `G″` is unresolved at pilot statistics (`+0.03±0.08`, broad
> marginal `sd=2.15`), and the rigorous `L`-uniform bound is `H1` [O]. The decisive run is a
> Landau+FP+Gribov-controlled `L=8,12(16)`, `N≥2000` zero-momentum susceptibility with the
> `tr F²` glueball mass as the gauge-invariant anchor — a statistical plateau statement, never a
> proof. No QED. Clay 0%.**

```
py -3.12 -I research/v3/w7/q5_4d_chi.py          # full pilot (~1 min, one job)
py -3.12 -I research/v3/w7/q5_4d_chi.py quick    # reduced N
```

*Refs.* `P2-PERIODIC-K0.md` (1D periodic plateau `Var→4.45`, `G″→0.222`), `D2-ANALYTIC-LIMIT.md`
(constant-mode `Var` plateau, free acid test), `P-WAVE-SYNTHESIS.md` (periodic floor `≈0.222`,
4D still `H1`), `SM6-SOFT-LATTICE.md` (decisive-run protocol, `Δ+W_geo+X`, gauge systematics),
`research/unconditional/H1-NUMERIC-RIGHT.md` (composite vacuous, field mode gauge-dependent — the
no-go this pilot respects), `docs/paper-v3/main.tex` §6.8 (decompactification);
tools `INVENT-GAUGE-INVARIANT-MASS.md` (A4 `m_GI`), `INVENT-VARIANCE-INEQ.md` (A2 `G″=S−D`).
Witness `research/v3/w7/q5_4d_chi.py`, data `research/v3/w7/q5_4d_chi_results.json`.
*No third-party sources downloaded.*
