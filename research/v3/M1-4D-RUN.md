# M1-4D-RUN — the decisive 4D SU(2) measurement at β=2.5: the flat-chart `k=0` susceptibility **χ=Var(a) plateaus** (3.81, 3.76, 3.82 at L=4,6,8 at N=2000, 3.61 at L=12; L=6→8 ratio `1.02±0.03` vs the free-massless `1.78`), the gauge-invariant Polyakov correlator **is gapped** (`ξ≲1`), but the **same coordinate in Landau gauge grows** (37→59→115) — the plateau is carrier-limited, the `0⁺⁺` glueball is **unresolved**, and the uniform-in-`L` bound is still `H1`

> **2026-11-24 · Wave M1 — the decisive 4D run specified by Q5-4D-CHI §7, executed.**
> β=2.5, SU(2), periodic 4D, L=4,6,8 each at **N=2000** stored configs (therm=800, gap=10,
> `τ_int≈1.0–1.1`, `N_eff≈1840–1950`) + an **L=12** N=300 feasibility point, **one sequential job**
> (charter §4; no OOM). Measured:
> (A) `χ = Var(a)`, `a = N^{-1/2}Σ_l log U_l`, k=0 mode — **flat chart** (continues the Q5 pilot)
> and **Landau gauge with Gribov copies** (5 random starts, 48 configs); (B) the gauge-invariant
> `0⁺⁺` correlator `⟨O(t)O(0)⟩_c`, `O=`space-summed plaquette (unsmeared + APE₂) for `m₁, f₁²`;
> (C) the gauge-invariant Polyakov spatial correlator `C_P(r)` → `ξ`.
> Witness `research/v3/w7/m1_4d_run.py`, post-analysis `research/v3/w7/m1_4d_analyze.py`,
> results `research/v3/w7/m1_4d_run_results.json`. **No third-party sources.**
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[V]** numerical witness (this run) ·
> **[C]** conditional / protocol / interpretation · **[NG]** no-go for a named class · **[O]** exact
> missing step. **No QED. Clay 0%.**

---

## 0 · Verdict (no QED; Clay 0%)

1. **[V] The `k=0` flat-chart susceptibility plateaus — from `L=4` to `L=12`.** `χ_flat = 3.807±0.074
   (L4)`, `3.760±0.072 (L6)`, `3.823±0.069 (L8)` at N=2000 each, plus a feasibility point
   `3.615±0.198 (L12, N=300)`; ratios `L4→6 = 0.987±0.027`, `L6→8 = 1.017±0.027`,
   `L8→12 = 0.945±0.058`. A free **massless** field would rise by `1.78` over `6→8` and `2.25` over
   `8→12` (`Var∝L²`, Q5 §3 acid test); it rises by `1.02` and `0.95`. **Growth is excluded at ≈28σ
   (6→8) and ≈22σ (8→12).** The plateau value `≈3.6–3.8` sits at the Haar floor
   `s²=(4π²−6)/9=3.7198` (β=0 exact) and below the 1D periodic plateau `4.45` (P2). The pilot's
   `L=4,6` plateau (`4.17, 3.61`) is confirmed and sharpened at N=2000 and extended to `L=12`.
2. **[V] The mean-subtracted variance equals the raw value** (mean-contribution `0.002–0.006`):
   `χ_flat` is a genuine variance, `⟨A⟩≈0` by symmetry — not a mean artifact.
3. **[V] The gauge-invariant Polyakov correlator is gapped.** `C_P(r)` decays to `C_P(1)/C_P(0)=0.073`
   by `r=1` at `L=8`; periodic cosh fit `ξ = 0.99, 0.61, 0.39, 0.25` (L=4,6,8,12), `ξ_{nn}=1.02,0.55,0.38,0.25`.
   `ξ≲1` lattice spacing and shrinking with `L` — a bounded (gapped) N-ality-1 sector, **consistent
   with** the plateau and a mass gap.
4. **[V]/[NG] The plateau is carrier-limited — the carrier hazard is confirmed numerically.**
   The **same log-sum coordinate**, measured in **Landau gauge** (Gribov-copy controlled), gives
   `χ_Landau = 36.98±1.86 (L4)`, `58.61±3.81 (L6)`, `115.30±10.42 (L8)` — it **grows** (ratios
   `1.59, 1.97`), `10–30×` the flat value. The growth is genuine variance (mean-contribution
   `0.07–0.14`; A-vectors point in random colour directions, `dir_align≈0.1–0.36`), not a coherent
   mean shift. So the flat-chart plateau is **not gauge-robust**: it is a property of the flat-chart
   coordinate, not a gauge-invariant statement. This is exactly Q6's carrier hazard.
5. **[V] The `0⁺⁺` glueball mass is not resolved at `L=8`.** The connected plaquette-plaquette
   correlator is flat in Euclidean time (`m_eff(t)=ln C(t)/C(t+1) ≈ 0.000–0.004` for all t, both
   unsmeared and APE₂, all L), and the exponential fit runs to the lower bound (`m₁→0.03`). Hence
   `m₁` and the residue `f₁²` are **not extractable** at L≤8 with this operator/statistics — matching
   `docs/paper-v2/ch12-numerics.tex` ("glueball two-point mass not resolved at L=8, fit m₀≳π").
6. **[O] The uniform-in-`L` bound (`H1`), and the atom, remain open.** A finite dataset at fixed `β`
   cannot certify `χ(L)≤1/λ` for all `L`; and the gauge-invariant non-self-averaging carrier (the
   `0⁺⁺` Källén–Lehmann atom with `f₁²`) is not resolved here. `H1` is a statistical plateau
   statement, never a proof.

**One line.** *At β=2.5 the flat-chart `k=0` susceptibility does not grow — it is flat to `±2%` across
`L=4,6,8` at N=2000 per volume (and `L=12` at N=300), growth excluded at ≈28σ (6→8) and ≈22σ (8→12) —
and the gauge-invariant Polyakov correlator is gapped (`ξ≲1`); **but** the same coordinate grows 10–30×
under Landau gauge, so the plateau is carrier-limited, the `0⁺⁺` glueball is unresolved, and the
rigorous uniform bound is still `H1`. No QED, Clay 0%.*

---

## 1 · The object and the estimators  [T]/[C]

H1's object is the `k=0` mode of the flat-chart log field: `ξ_l = log U_l` (the exact `SU(2)` log,
`x = θ n ∈ ℝ³`), `a = N^{-1/2}Σ_l ξ_l`, `N=4L⁴`, and

```
χ ≡ Var(a) = (1/3)⟨|A|²⟩ ,     A = N^{-1/2}Σ_l ξ_l ∈ ℝ³      (colour-averaged k=0 amplitude).
```

`⟨A⟩=0` under the Haar measure (verified: mean-contribution `≤0.006`), so `(1/3)⟨|A|²⟩` **is** the
variance. This is the Q5-pilot object (D2/SM6 normalization) and the periodic `k=0` Fourier mode of the
P-wave (β=0 exact `Var→3.7198`; 1D periodic plateau `4.45`, `G″→0.222`).

| tag | definition | role |
|---|---|---|
| **χ_flat** | `(1/3)⟨|A|²⟩`, `A=N^{-1/2}Σ_l log U_l` (flat chart) | primary `Var(a)`; continues Q5 (L=4,6) |
| **χ_Landau** | same, on **Landau-gauge-fixed** configs, `N_copy=5` random starts; Gribov spread | gauge-fixed control (Q5 §7) |
| **C_P(r)** | gauge-invariant `⟨P(x)P(x+r)⟩−⟨P⟩²`, `P=½tr∏_tU_0` → `ξ` | genuinely gauge-invariant cross-check |
| **0⁺⁺ C(t)** | `⟨O(t)O(0)⟩_c`, `O(t)=Σ_space Σ_{i<j} plaquette` (unsmeared, APE₂) → `m₁, f₁²` | scalar glueball attempt |
| **free acid** | exact Gaussian: massive `Var=1/m²`, `G″=m²`; massless `Var∝L²`, `G″→0` | discriminator (must fail for free) |

### 1.1 Acid test — the tool manufactures no floor  [T]

Exact (Q5 §3, reproduced): a free **massive** field gives `Var=1/m²` (`L`-flat) and `G″=m²`; a free
**massless** field (box-regulated, `k_min=2π/L`) gives `Var/L²≈0.026` const (`Var∝L²`) and `G″∝1/L²→0`:

```
massive m:   Var=1/m^2      1.000(m=1)   11.11(0.3)   1111(0.03)   (plateau, G''=m^2)
massless L:  Var   Var/L^2  G''          L=4: 0.500 0.0312 2.000 | L=8: 1.707 0.0267 0.586 | L=16: 6.569 0.0257 0.152
```

The measured flat-chart `χ/L² = 0.238 (L4), 0.104 (L6), 0.060 (L8)` **decreases** (not const) and
`χ` itself is flat — the data are *massive-like*, not massless-free. **The tool discriminates.** `[T]`

### 1.2 Landau gauge fixing + Gribov control  [C]

Red-black over-relaxation-free Landau gauge fixing (maximise `F=Σ Re tr V_μ`), exact `SU(2)` local
optimizer `g=proj_{SU(2)}(g_opt)` with `g_opt` the SU(2) element maximising `Re tr(g·Δ)`,
`Δ=Σ_μ(V_μ(x)+V_μ(x-μ)†)`; converges to `⟨(1−g₀)⟩<10⁻⁷` (L=8 in ≈166 iterations, `F` monotone up).
Per config, `N_copy=5` random gauge starts; the **spread of `χ` across copies is the Gribov systematic**
(48 configs). The log coordinate is computed on the fixed links; the ensemble variance is
mean-subtracted (`Var=|⟨A|²⟩−|⟨A⟩|²|/3`) so a coherent drift cannot masquerade as variance.

---

## 2 · The decisive run — results at β=2.5, N=2000 per volume  [V]

Single sequential process, one job (`charter §4`). Therm=800 sweeps, gap=10 sweeps between stored
configs. `τ_int` and `N_eff` from the measured series; delete-block jackknife SEs (20 blocks).

| quantity | **L=4** | **L=6** | **L=8** |
|---|---|---|---|
| **χ_flat = Var(a)** | **3.807 ± 0.074** | **3.760 ± 0.072** | **3.823 ± 0.069** |
| χ_flat (mean-subtracted) | 3.806 ± 0.075 | 3.759 ± 0.071 | 3.817 ± 0.068 |
| mean-contribution `|⟨A⟩|²/3` | 0.002 | 0.001 | 0.006 |
| **τ_int(χ_flat)**, **N_eff** | 1.03, 1936 | 1.03, 1947 | 1.09, 1837 |
| sd(χ_flat) per config | 3.06 | 3.19 | 3.25 |
| `⟨plaq⟩` | 0.6548 | 0.6528 | 0.6524 |
| **χ_Landau (Var)** | 36.98 ± 1.86 | 58.61 ± 3.81 | 115.30 ± 10.42 |
| τ_int(χ_Landau) | 1.86 | 2.84 | 6.58 |
| Gribov spread rel (mean/max) | 1.3% / 21% | 11.2% / 121% | 13.0% / 66% |
| **C_P(0),C_P(1),C_P(2),…** | `.270,.101,.085,.101` | `.254,.041,.025,.022,…` | `.250,.018,.007,.006,.005,…` |
| **ξ (periodic cosh)**, ξ_nn | 0.99, 1.02 | 0.61, 0.55 | 0.39, 0.38 |
| 0⁺⁺ `m_eff(t)` (unsp./APE₂) | ≈0 | ≈0 | ≈0 (unresolved) |

> **L=12 feasibility point** (N=300, gap=6, no Landau/Gribov; one job): `χ_flat = 3.615 ± 0.198`
> (`τ_int=1.31, N_eff=229`), `⟨plaq⟩≈0.652`, `C_P(1)/C_P(0)=0.018` (`ξ≲0.2`). Used only to extend the
> plateau trend (no gauge-fixed/glueball data at this volume).

### 2.1 Plateau, not growth  [V]

```
        χ_flat :  3.807 → 3.760 → 3.823 → 3.615     (L4, L6, L8, L12*)
        ratio  :    0.987±0.027 (4→6)   1.017±0.027 (6→8)   0.945±0.058 (8→12*)
   free-massless:    2.25               1.78                2.25            (excluded ≈28σ, ≈22σ)
                        * L12 = N=300 feasibility point (no Landau/Gribov)
```

Across the full range `L4→8` the ratio is `3.823/3.807 = 1.004`, and `L4→12` is `3.615/3.807 = 0.950`,
against free-massless `(12/4)²=9`. **The flat-chart `k=0` susceptibility is flat.** The plateau value
`3.6–3.8` is `1.0–1.02×` the exact Haar floor `3.7198` and below the 1D periodic `4.45`. This is the
4D analogue of the P2/D2 plateau, now at `L=8` with `±2%` errors and extended to `L=12`. `[V]`

### 2.2 Gauge-invariant cross-check: the Polyakov sector is gapped  [V]

`C_P(r)` (N-ality-1, fundamental flux) decays within one lattice spacing, faster as `L` grows
(`C_P(1)/C_P(0)=0.375, 0.163, 0.073` at L=4,6,8), with `ξ_{cosh}=0.99, 0.61, 0.39`. A correlation
length `≲1` and **decreasing** with volume is the signature of a gapped, non-critical sector —
**consistent with** the plateau and with a mass gap. It is a *different* (string/gluon) channel from
the scalar susceptibility `χ`, so it corroborates but does not prove the scalar plateau.

### 2.3 The `0⁺⁺` glueball is unresolved  [V]

The connected plaquette-plaquette correlator is flat in time: `m_eff(t)∈[−0.005,0.004]` for every
`t`, both smearing levels, all `L`; the exponential fit is driven to `m₁→0.03` (the lower bound). The
`0⁺⁺` mass and residue `f₁²` are **not extractable** at `L≤8` with this operator — the known
signal-to-noise wall (ch12-numerics). The residue `f₁²` is therefore **not measured**; the atom is
untouched. `[V]`→`[O]`

---

## 3 · The carrier hazard, made explicit  [V]/[NG]

The decisive new negative: **the same log-sum coordinate grows under gauge fixing.**

```
   χ_flat    :  3.81 → 3.76 → 3.82        (flat, ratio ≈ 1)
   χ_Landau  : 36.98 → 58.61 → 115.30     (grows, ratios 1.59, 1.97; 10–30× the flat value)
```

The growth is a **variance**, not a drift: mean-contribution `0.07–0.14` (tiny); the per-config `A`
vectors point in random colour directions (`|mean of unit A| ≈ 0.14–0.36`), so `⟨A⟩≈0` and
`Var≈⟨|A|²⟩`. Landau gauge aligns the link logs coherently (smoothing), so the zero-mode sum
`Σ_l log U_l` acquires fluctuations that grow with volume (χ ratio `1.59` and `1.97` across
`L=4→6→8`, i.e. an effective growth between `∝L²` and steeper), while the flat chart's logs cancel.
The flat-chart value is thus **not gauge-robust** — exactly the Q6 carrier hazard: on the
gauge-dependent gluon carrier H1 is gauge-dependent; on the gauge-invariant composite `tr F²` it is
vacuous (`Var∝1/V`, Q6/H1-NUMERIC-RIGHT); the gauge-invariant non-self-averaging carrier is not
constructed. The Landau value is also **ill-conditioned** (heavy tails from near-π links, `τ_int` up
to 6.6, Gribov spread to 13%–66%), so it is a *control*, not a physical susceptibility. `[V]`/`[NG]`

**Reading.** The plateau claim of §2.1 is a statement about the **flat-chart `k=0` coordinate**, the
literal H1 object, and it is decisive there. It is **not** a gauge-invariant mass-gap measurement; the
gauge-invariant anchor (Polyakov) is consistent (gapped) but lives in another channel. The honest
conclusion is the Q6 one: H1 is a necessary ceiling on a definite carrier, not the finish line, and the
carrier that would make it non-vacuous is the missing structural object `[O]`.

---

## 4 · L=12 — measured feasibility point + the full protocol  [V]/[C]

**Measured (feasibility point).** `χ_flat = 3.615 ± 0.198` (N=300, gap=6, one job, `τ_int=1.31`,
`N_eff=229`), `⟨plaq⟩≈0.652`, `C_P(1)/C_P(0)=0.018` (`ξ≲0.2`). This **extends the plateau to `L=12`**
(`L8→12` ratio `0.945±0.058`, vs free-massless `2.25`) — growth excluded at ≈22σ. Stored as `L12` in
`m1_4d_run_results.json`.

**Full `L=12` (specified, not run).** A full decisive `L=12` (**N=2000, Landau+Gribov, glueball**)
would add the gauge-fixed control and the scalar correlator at the largest volume. Measured cost
extrapolation (this machine, pure numpy, one job): sweeps `(12/8)⁴≈5.1×` (`≈78 min by measured
L=8 rate`), observables `≈15 min`, Landau fixing `700×(≈7–9 s)≈90–105 min`, Gribov `48×4×7 s≈22 min`
→ **≈3.5–4 h**, one job, **never parallel** (charter §4). It is **specified, not run**: the plateau is
already decisive (`L=4,6,8` at N=2000, `L=12` feasibility), and the machine limit forbids parallel big
jobs. `[C]`

---

## 5 · Status ledger

| item | statement | label |
|---|---|---|
| free-field tool: massive `Var=1/m²` flat, massless `Var∝L²`, `G″→0` | §1.1 | **[T]** |
| **flat-chart `χ=Var(a)` plateaus: 3.81/3.76/3.82 (L4/6/8, N=2000) + 3.61 (L12, N=300); ratios 0.99/1.02/0.95, growth excluded ≈28σ/≈22σ** | §2.1 | **[V]** |
| `χ_flat` is a genuine variance (mean-contribution `≤0.006`) | §2 | **[V]** |
| gauge-invariant Polyakov `ξ≲1`, decreasing with `L` → gapped sector, consistent | §2.2 | **[V]** |
| `0⁺⁺` glueball unresolved at `L=8` (`m_eff≈0`); `m₁, f₁²` not extractable | §2.3 | **[V]** |
| **Landau-gauge same coordinate grows 37→59→115 (10–30× flat), genuine variance** | §3 | **[V]** |
| plateau is carrier-limited: gauge-fixed coordinate not gauge-robust | §3; Q6 | **[NG]** |
| `0⁺⁺` Källén–Lehmann atom (`f₁²≥r₀`), gauge-invariant non-self-averaging carrier | — | **[O]** |
| uniform-in-`L` `Var(a)≤1/λ` (= `H1`), `G″≥λ` | — | **[O]** |
| L=12 `χ_flat` feasibility point: `3.615±0.198` (N=300), extends plateau | §4 | **[V]** |
| L=12 full run (N=2000, Landau+Gribov, glueball) | specified §4 (`≈3.5–4 h`) | **[C]** |
| Clay mass gap | open | **0%** |

---

## 6 · Sharpest line and reproduce

> **At β=2.5, SU(2), periodic 4D, N=2000 per volume, the flat-chart `k=0` susceptibility
> `χ=Var(a)=3.81, 3.76, 3.82, 3.61` at `L=4,6,8,12` — flat to `±2%` (`L12` is an N=300 feasibility
> point), the `L=6→8` ratio `1.02±0.03` and `L=8→12` ratio `0.95±0.06` against free-massless `1.78`/`2.25`
> (growth excluded at ≈28σ/≈22σ). The gauge-invariant Polyakov correlator is gapped (`ξ≲1`, decreasing
> with `L`), consistent with the plateau and a mass gap. But the same log-sum coordinate grows `10–30×`
> under Landau gauge (`χ=37→59→115`), a genuine variance with Gribov systematics `1–13%`: **the plateau
> is carrier-limited, not a gauge-invariant mass-gap statement.** The `0⁺⁺` glueball is unresolved at
> `L=8` (`m_eff≈0`), so `m₁` and `f₁²` are unmeasured and the atom is untouched. The uniform-in-`L`
> bound `H1` remains `[O]`; the run is a statistical plateau on a gauge-dependent coordinate, never a
> proof. No QED. Clay 0%.**

```
py -3.12 -I research/v3/w7/m1_4d_run.py 4 6 8      # decisive run, N=2000/vol, ~65 min, ONE job
py -3.12 -I research/v3/w7/m1_4d_run.py quick 4    # smoke test (~10 s)
py -3.12 -I research/v3/w7/m1_4d_analyze.py        # clean summary from the JSON
```

*Refs.* `Q5-4D-CHI.md` (pilot L=4,6 + decisive-run spec), `Q-WAVE-SYNTHESIS.md` (four faces of H1,
carrier hazard), `Q6-H1-STATUS.md` (carrier hazard, H1 ⇏ atom), `P-WAVE-SYNTHESIS.md` (periodic
`k=0`, Haar `3.7198`, 1D plateau `4.45`), `SM6-SOFT-LATTICE.md` (Landau+Gribov protocol),
`docs/paper-v2/ch12-numerics.tex` (composite vacuous, glueball unresolved at L=8, decisive-test
requirements), `H1-NUMERIC-RIGHT.md`; tools `INVENT-GAUGE-INVARIANT-MASS.md`,
`INVENT-VARIANCE-INEQ.md`. Witness `research/v3/w7/m1_4d_run.py`, analysis `m1_4d_analyze.py`, data
`m1_4d_run_results.json`. *No third-party sources downloaded.*
