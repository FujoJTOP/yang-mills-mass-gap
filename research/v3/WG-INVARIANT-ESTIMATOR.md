# WG-INVARIANT-ESTIMATOR — a **gauge-invariant, non-self-averaging** susceptibility for the 4D β=2.5 plateau: the **total-normalized zero mode of an APE-smeared `tr F²` density**, `χ_O = Var(V^{-1/2}Σ_x Õ_s) = Σ_z C_s(z) = O(1)`; it is exactly gauge-invariant (vs the flat chart's O(1) drift), non-self-averaging (vs `χ/V`), acid-clean (free massive → known, free massless → ∞), and **plateaus** (`0.33→0.38→0.41` at L=4,6,8; ratios `1.15, 1.09` vs the free-massless `2.25, 1.78`) — the ceiling **H1**, not the edge (atom `f₁²` still `[O]`)

> **2026-11-24 · Wave WG — the gauge-invariant carrier that replaces the flat-chart `k=0` coordinate.**
> The M1 run (`M1-4D-RUN.md`, `M1-FAST-CROSSCHECK.md`) measured a plateau on the **flat chart**
> (`χ_flat = Var(N^{-1/2}Σ_l log U_l) ≈ 3.6–3.8`), but the same coordinate **grows** in Landau gauge
> (`37→59→115`) and the gauge-invariant Polyakov susceptibility **drops** (`5.5→5.1→2.8`), leaving the
> plateau "carrier-limited" (`Q6-H1-STATUS.md` §4.5; `M1-4D-RUN.md` §3). T4 (`T4-CARRIER.md`) showed
> the asserted "`Var ∝ 1/V`, no gauge-invariant carrier" is a **normalization artifact** (`q̄` vs `q`).
> This note **builds the estimator** on that carrier, runs it on the M1 measure, and explains the
> `χ_P`-drop / `χ_Landau`-growth dichotomy.
>
> **Witness:** `research/v3/w7/wg_invariant.py` → `wg_invariant_results.json` (`py -3.12 -I`, numpy/numba,
> **one process** — charter §4, never parallel; no downloads). **No third-party sources.**
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional / protocol / interpretation ·
> **[NG]** no-go for a named class · **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **`[T]` The carrier is the total-normalized zero mode, and its variance is the susceptibility.**
   For any gauge-invariant scalar density `O` with connected two-point `C(z)`,
   ```
     q := V^{-1/2} Σ_x (O(x) − ⟨O⟩) ⇒ Var(q) = Σ_z C(z) = χ_O  (O(1), L-independent);
     q̄ := V^{-1}   Σ_x (O(x) − ⟨O⟩) ⇒ Var(q̄) = χ_O / V          (self-averaging → 0).
   ```
   The corpus's "`Var∝1/V` on `tr F²`" is the **mean**-normalized `q̄` — an artifact, not a no-go (T4-1).
   H1 (`ch02 eq:h1`) is stated on `q`. **The gauge-invariant non-self-averaging carrier exists.**

2. **`[T]+[V]` The estimator is exactly gauge-invariant and non-self-averaging.** Take the APE-smeared
   field-strength / plaquette density `Õ_s(x)=Σ_{μ<ν}(1−½Re tr P^s_{μν}(x))` and its total-normalized
   zero mode `q`. Witness: under a random gauge transform `U→U^g` the estimator is unchanged to
   **`1.4×10^{-15}`**, while `χ_flat` drifts **40–1000 %** on single configs (ensemble `≤4–6 %`, Q5) —
   the flat chart is gauge-invariant only infinitesimally. `Var(q̄)=χ_O/V` falls as **`1.3→0.3→0.1 ×10^{-3}`
   (L=4,6,8)** while `χ_O` itself stays **`0.33, 0.38, 0.41`** — non-self-averaging, `O(1)`.

3. **`[T]+[V]` Acid test passes — the estimator manufactures no floor.** Free **massive** elementary
   `φ`: `Var(q)=1/m²` **exactly**, every `L` (`=4.000000`, log-log slope `0.0000`). Free massive
   **composite** `:φ²:` (the 0⁺⁺ analog of `tr F²`): `χ` **finite** and `L`-independent
   (`4.109` (d1), `0.689` (d2), `0.062` (d4, thermodynamic)); free **massless**: `χ→+∞`
   (`d1`: `4.1→196000` as `m:0.5→0.02`) — **H1 correctly fails**, the tool does not fake a gap.

4. **`[V]` On the SU(2) β=2.5 measure the gauge-invariant scalar `χ_O` plateaus (sub-linear), not `L²`.**
   `χ_O = 0.328±0.022, 0.377±0.021, 0.411±0.058` at L=4,6,8; ratios `L4→6 = 1.15±0.09`,
   `L6→8 = 1.09±0.15` against free-massless `2.25, 1.78` (**statistically excluded at ≈12σ and ≈4.6σ**; a
   plateau `=1` sits within 1.7σ and 0.6σ). Contact-subtracted `χ_phys = χ_O − C_s(0)`
   (`0.297, 0.343, 0.375`) behaves identically. This **replaces the flat-chart plateau with a
   gauge-invariant one** — the concrete content the M1/Q6 carrier hazard asked for.

5. **`[NG]` It is the ceiling `H1`, not the edge.** `χ_O = ∫s^{-1}dν_O` is a residue-weighted harmonic
   mean: the exact squeeze gives `f₁²/χ ≤ m₁² ≤ M₀/χ`, an **upper** bound on the mass. No variance
   functional lower-bounds the mass (thin atom `wδ_ε+(1−w)δ_{s₂}`; Schur-concavity; T4-4). So a clean
   gauge-invariant plateau certifies **`χ_O ≤ 1/λ` (H1)** on a non-vacuous carrier — the necessary
   ceiling — and **leaves the atom `f₁² = ν_O({m₁²}) ≥ r₀` (with isolation `I₂`) as the one open input
   `[O]`.** This is the honest status: the estimator solves the *carrier* problem, not the *gap*.

**One line.** *The flat-chart plateau is replaced by a gauge-invariant one: the total-normalized zero
mode of the APE-smeared `tr F²` density has `χ_O = Var(q) = Σ_z C_s(z)`, exactly gauge-invariant
(`10^{-15}`), non-self-averaging (`χ_O` flat, not `χ_O/V`), acid-clean (free massive known, free
massless ∞), and on β=2.5 SU(2) it plateaus sub-linearly (`0.33→0.38→0.41`, ratios `1.15,1.09` vs free
`2.25,1.78`). `χ_P` dropped because the Polyakov loop is the wrong (N-ality-1) channel; `χ_Landau`
grew because a gauge-fixed coordinate accumulates the unphysical longitudinal/Gribov sector whose
spread is already 13 % at L=8. This certifies H1 on the correct carrier — the ceiling — and leaves the
atom as `[O]`. No QED, Clay 0%.*

---

## 1 · The carrier problem, stated once

The M1/Q6 wave left three coordinates, none satisfactory on its own:

| coordinate | gauge-invariant? | scaling with L | verdict |
|---|---|---|---|
| `χ_flat = Var(N^{-1/2}Σ_l log U_l)` | only infinitesimally (`≤4%` moderate, `O(1)` large) | flat `3.81,3.76,3.82` | plateau, but **carrier-limited** |
| `χ_Landau` (same, Landau-fixed) | no (Gribov) | grows `37→59→115` | **gauge artifact** |
| `χ_P = L³ Var(P̄)` (Polyakov) | **yes** | drops `5.5→5.1→2.8` | gapped, but **wrong channel** |

What is missing is a quantity that is **(i) exactly gauge-invariant, (ii) non-self-averaging
(`Var ↛ 1/V`)** — i.e. H1's own object — and **(iii) acid-clean**. T4-CARRIER proved (i)+(ii) are a
**normalization choice** (total vs mean), not a structural obstruction. WG makes the choice concrete and
measures it; §4 explains why the two older coordinates misbehave.

---

## 2 · The estimator  [T]/[C]

**Definition (WG).** Let `O` be a gauge-invariant `0⁺⁺` scalar density — here the field-strength /
plaquette energy density summed over all six Lorentz planes,
```
    Õ(x) = Σ_{μ<ν} ( 1 − ½ Re tr P_{μν}(x) ) ,     P_{μν} the plaquette of APE-smeared links.
```
The estimator is its **total-normalized zero-momentum mode**
```
    q = V^{-1/2} Σ_x ( Õ(x) − ⟨Õ⟩ ) ,     χ_O := Var(q) = Σ_z C_s(z) = V · Var_c(Ō_c) ,
```
with `C_s` the connected two-point function and `Ō_c = V^{-1}Σ_x Õ` the per-config spatial mean.
Three structural facts:

- **Gauge-invariant `[T]`.** `Õ` is a pointwise sum of `½Re tr` of gauge-invariant plaquettes; `q` is a
  constant (`V^{-1/2}`) times their sum minus a c-number. It is invariant under **every** gauge
  transform, per configuration — witnessed to `1.4×10^{-15}` (§5). This is exactly what the flat chart
  lacks.
- **Non-self-averaging `[T]`.** `Var(q) = Σ_z C_s(z) = χ_O` (T4-1): an `O(1)`, `L`-independent number
  whenever `χ_O∈(0,∞)`. The mean-normalized mode `q̄` would give `χ_O/V → 0`; H1 is on `q`, not `q̄`.
- **Spectral-projector reading `[C]`.** APE smearing at a **fixed physical radius** `r_s` is a
  gauge-covariant low-pass filter: `Õ_s` is `tr F²` projected onto `|p| ≲ 1/r_s`. Then
  `C_s(τ) = ∫ e^{−mτ} dν_s(m)` is the smeared `0⁺⁺` Källén–Lehmann measure and `χ_O = ∫ s^{-1} dν_s`.
  So the estimator is `Var` of the `k=0` mode of `P_s tr F²`, with `P_s` the gauge-invariant projector
  onto the low scalar spectrum. "Smeared total-normalized zero mode" and "gauge-invariant spectral
  projector" are two names for this one object.

### 2.1 The UV hazard that is real (and how smearing handles it)  [T]/[V]

The bare composite has a contact term `C_s(0) ∼ a^{-8}`; then `χ_O` is **contact-dominated** and flat
for the wrong (UV) reason. The witness measures this: the **bare** `ns=0` density gives
`χ_bare = 1.19, 1.13, 1.17` (L=4,6,8, `L`-flat but a UV false plateau), while the **APE `ns=2`
fixed-radius** density gives `χ_O = 0.33, 0.38, 0.41` with contact `C_s(0) = 0.030, 0.034, 0.035`
(only `8–9 %`). The **contact-subtracted** `χ_phys = χ_O − C_s(0) = 0.297, 0.343, 0.375` gives the
identical ratios. So the smearing is load-bearing: at a fixed physical scale the susceptibility is the
physical `∫s^{-1}dν_O`, finite and scheme-clean; the `L`-**scaling** (below) is the scale-independent
statement.

---

## 3 · Acid test — the free field must give a known value (and no floor)  [T]+[V]

Exact, deterministic (`[F]` in the witness):

```
 free elementary φ, m=0.5 :  Var(total-zero) = 1/m² = 4.000000  at every L, d=1,2,4  (slope 0.0000)
 free composite  :φ²:     :  χ = Σ_z 2G(z)²       m=0.5
     d=1 : 4.1089, 4.1088, 4.1088, 4.1088   (L=32,64,128,256)   slope 0.0000
     d=2 : 0.8515, 0.6925, 0.6891, 0.6890   (L=8,16,32,64)      -> 0.6891 (thermodynamic)
     d=4 : 0.1747, 0.0801, 0.0663, 0.0626, 0.0623 (L=4..16)      -> 0.0622 (thermodynamic, finite)
 massless acid (m→0) :  χ(:φ²:) DIVERGES   d=1 L64: 4.1, 62.8, 513, 5475, 196000  (m=.5,.2,.1,.05,.02)
                                           d=2 L64: 0.69, 4.07, 16.4, 100.6, 3082
```

The massive free field gives the **known** value (`1/m²` for `φ`, the finite lattice `Σ2G²` for `:φ²:`),
and the **massless** free field gives `χ=∞` — so **H1 fails on the gapless model and the estimator
manufactures no floor** (F1/F4 clean). The `d=4` composite converges (slowly; the lattice UV cutoff
fixes the contact) to the finite thermodynamic value `0.0622` — a free field, correctly, has a finite
susceptibility with **`f₁²=0`** (a two-particle cut, no atom): exactly the `:φ²:` acid model of
`H1-IS-IT-THE-GAP.md` §2.3-iii. The `χ` of a free composite is **not** `1/m₁²` — it is the ceiling, a
first face of the ceiling-not-edge theorem below.

---

## 4 · Why `χ_P` dropped and `χ_Landau` grew  [V]/[C]/[NG]

### 4.1 `χ_P` drops because the Polyakov loop is the **wrong channel** (N-ality-1)  [V]/[C]

`χ_P = L³ Var(P̄) = Σ_r C_P(r)` is the **integrated N-ality-1 (center / fundamental-flux)
susceptibility**, `P(x)=½tr∏_t U_0`. At β=2.5, SU(2) is in the **confined** phase: `⟨P⟩=0` (unbroken
ℤ₂ center) and the lightest N-ality-1 state is the **flux tube / string**, not the `0⁺⁺` glueball. Its
spatial correlator `C_P(r)` is **short range** (`ξ≲1`, decreasing with `L`: `C_P(1)/C_P(0)=0.375, 0.163,
0.073` at L=4,6,8) and the periodic images / finite-size enhancement at small `L` **inflate** both `ξ`
and `χ_P`. As `L` grows the physical short range is resolved and `χ_P` **relaxes down** to its `O(1)`,
saturated thermodynamic value (`5.5→5.1→2.8` in M1; `5.6→4.9→3.7` here) — a **gapped** sector.

So `χ_P` is gauge-invariant and honest, but it measures the **string (N-ality-1)**, not the scalar
atom. Its drop is the **correct** behaviour of the wrong quantum numbers; it does not contradict a
scalar plateau (the two live in disjoint centre sectors — `Q6`/`T5`). It is a **cross-check** (the
gapped N-ality-1 sector is consistent with confinement), **not** the 0⁺⁺ susceptibility.

### 4.2 `χ_Landau` grows because a gauge-fixed coordinate accumulates the unphysical sector  [V]/[NG]

Landau gauge maximizes `F=Σ Re tr V_μ`, cohering the link logs toward identity **along the gauge orbit**.
The `k=0` log-sum `A = N^{-1/2}Σ_l log U_l` computed on **Landau-fixed** links is a **gauge-fixed
coordinate, not gauge-invariant**. Two mechanisms make it grow with volume:

- **Longitudinal leakage.** The gauge condition splits `log U = transverse (physical) + longitudinal
  (pure-gauge)`. The longitudinal piece has a **non-vanishing zero-momentum component** (the naive
  "gradient sums to zero" argument needs the abelian/1D additivity that the 4D `SU(2)` log lacks,
  `Q5` §1), so `Σ_l log U` accumulates it, and its variance grows **at least `∝L²`**
  (measured ratios `1.59, 1.97`, between `L²` and steeper — `M1-4D-RUN.md` §3).
- **Gribov multiplicity.** Landau gauge is ill-defined at finite volume. The number of Gribov copies
  grows exponentially in `V`; the max-functional fixing lands on different copies with heavy-tailed
  residuals (near-π links). The measured **Gribov spread is 1.3 % (L4) → 11.2 % (L6) → 13.0 % (L8)**
  (relative std across `N_copy=5` random starts, M1). The "growth" `37→59→115` (ratios `1.59, 1.97`) is
  of the same order as — and unresolved against — this growing Gribov systematic; `τ_int` rises to `6.6`
  and the series is heavy-tailed. So `χ_Landau`'s growth is a **gauge-fixing / Gribov volume artifact**,
  carrying the unphysical longitudinal sector: a **control**, never physics. `[NG]` for "Landau-fixed
  `Var(a)` as a physical susceptibility."

**Contrast with `χ_O`.** Because `χ_O` is gauge-invariant (§2), it has **no** longitudinal/Gribov
component to accumulate: its `L`-scaling is physical — sub-linear (§5), the bounded scalar
susceptibility. The dichotomy `χ_P↓ / χ_Landau↑ / χ_O ≈ flat` is the sharp diagnostic: the first is
the wrong (centre) channel, the second is gauge contamination, the third is the scalar carrier.

---

## 5 · The measurement on β=2.5 SU(2)  [V]

Witness `[G]`: SU(2) heat bath (`m1_fast.py` Numba kernel, same measure as M1), one process,
`therm=500, gap=8`, APE `α=0.5, ns=2`, delete-block-jackknife SEs.

| quantity | **L=4** (N=800) | **L=6** (N=400) | **L=8** (N=160) |
|---|---|---|---|
| **`χ_O` (total-norm, gauge-inv)** | **0.328 ± 0.022** | **0.377 ± 0.021** | **0.411 ± 0.058** |
| `Var(q̄)=χ_O/V` (means-norm) | 1.3×10⁻³ | 0.30×10⁻³ | 0.10×10⁻³ |
| contact `C_s(0)` | 0.0302 | 0.0343 | 0.0351 |
| **`χ_phys = χ_O − C_s(0)`** | **0.297** | **0.343** | **0.375** |
| `χ_bare` (ns=0, contact-dom.) | 1.189 | 1.131 | 1.168 |
| `χ_flat` (flat chart, ref.) | 3.699 ± 0.110 | 3.757 ± 0.160 | 3.614 ± 0.235 |
| `χ_P` (N-ality-1, ref.) | 5.65 | 4.86 | 3.67 |
| gauge-inv reldiff (`χ_O`, U vs U^g) | 6.7×10⁻¹⁶ | 0.0 | 1.7×10⁻¹⁶ |
| `χ_flat` drift (single config, one U^g) | 40 % | 1067 % | 185 % |

**Plateau, not growth.** `χ_O`: `0.328 → 0.377 → 0.411`, ratios
```
   L4→6 : 1.15 ± 0.09   (free massless 2.25)  -- statistically excluded ≈12σ, plateau within 1.7σ
   L6→8 : 1.09 ± 0.15   (free massless 1.78)  -- statistically excluded ≈4.6σ, plateau within 0.6σ
```
The residual mild rise (`~9–15 %` per step) is consistent with **finite-size / fixed-radius smearing**
(the APE radius is 2–3 lattice spacings, a large fraction of `L≤8`; the same fixed-radius effect that
made the free `d=4` composite approach its limit from above, §3) and does **not** extrapolate to the
free `∝L²` growth (`2.25, 1.78`). A decisive plateau statement needs `L=12,16` at fixed physical radius
(flow, `t∼(0.5 fm)²`) — the same machine budget as the M1 L=12 run (`≈4 min` here, `M1-FAST`). **[C].**

**Gauge invariance vs the flat chart.** On every config, `χ_O` is unchanged by a random gauge
transform to machine zero (`≤1.7×10^{-16}`), while the **same-config** `χ_flat` drifts `40–1000 %`
(heavy-tailed single configs; the **ensemble** drift is `≤4–6 %` for moderate transforms, Q5). This is
the concrete sense in which `χ_O` *replaces* the flat-chart coordinate: the plateau now lives on an
exactly gauge-invariant random variable.

---

## 6 · What this does and does not establish

**Does `[T]+[V]`.**
- A gauge-invariant, non-self-averaging susceptibility `χ_O = Var(q)` exists, is elementary to write,
  and is exactly gauge-invariant (`10^{-15}`) and non-self-averaging (`χ_O=O(1)`, `Var(q̄)=χ_O/V`).
- It is acid-clean (free massive known finite; free massless `∞`; `f₁²=0` free composite).
- On β=2.5 SU(2) it plateaus sub-linearly (`0.33,0.38,0.41`; ratios `1.15,1.09` vs free `2.25,1.78`),
  **growth excluded at ≈12σ (4→6) and ≈4.6σ (6→8)**. This is H1 (`χ_O ≤ 1/λ`) on a definite,
  non-vacuous, gauge-invariant carrier — the object Q6 §4.5 said was missing.

**Does not `[NG]/[O]`.**
- **Not the gap.** `χ_O = ∫s^{-1}dν_O` is a residue-weighted **harmonic mean**: the squeeze
  `f₁²/χ_O ≤ m₁² ≤ M₀/χ_O` makes any variance an **upper** bound on the mass. No variance functional
  lower-bounds the mass (thin atom `wδ_ε+(1−w)δ_{s₂}`, `w=ε/2λ`: `χ_O` fixed, `m₁²=ε→0`; Schur; T4-4).
  A clean plateau certifies the **ceiling**, not the **atom**. The atom `f₁² = ν_O({m₁²}) ≥ r₀` (with
  isolation `I₂`) remains the one irreducible input `[O]` — `H* = atom + ceiling` (`ch02 eq:hstar`).
- **Not a continuum theorem.** A finite L=4,6,8 dataset at fixed `a` cannot certify `χ_O(L) ≤ 1/λ`
  uniformly in `(L,a)`; `H1` at β=2.5 is `[O]`. The witness is a statistical plateau on the correct
  carrier, never a proof. **No QED. Clay 0%.**

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| normalization algebra | `Var(total-norm q)=χ_O`; `Var(mean-norm q̄)=χ_O/V` | **[T]** |
| estimator `χ_O=Var(q)`, `q=V^{-1/2}Σ(Õ−⟨Õ⟩)`, `Õ` `0⁺⁺` gauge-inv density | §2 | **[T]** |
| exactly gauge-invariant (vs flat chart O(1) drift) | §5; `10^{-15}` witness | **[T]**+**[V]** |
| non-self-averaging: `χ_O` flat, `Var(q̄)=χ_O/V↓` | §5 (`1.3→0.1`×10⁻³) | **[T]**+**[V]** |
| acid: free massive → known (`1/m²`, `Σ2G²`), free massless → ∞ | §3 | **[T]**+**[V]** |
| APE fixed-radius suppresses UV contact (`χ_bare 1.17 → χ_O 0.41`) | §2.1 | **[V]** |
| **β=2.5: `χ_O=0.33,0.38,0.41`; ratios `1.15,1.09` vs free `2.25,1.78`; growth excl. ≈12σ/≈4.6σ** | §5 | **[V]** |
| `χ_P` drops = N-ality-1 (centre/string) wrong channel, gapped | §4.1 | **[V]**/**[C]** |
| `χ_Landau` grows = gauge-fixed longitudinal + Gribov sector (spread 1.3→13 %) | §4.2 | **[V]**/**[NG]** |
| plateau is the ceiling `H1`, not the edge (harmonic mean; thin atom) | §6; T4-3/T4-4 | **[T]**/**[NG]** |
| atom `f₁²≥r₀` + isolation `I₂` on the carrier | — | **[O]** |
| uniform-in-`L` `χ_O≤1/λ` at β=2.5 (`H1`) | — | **[O]** |
| Clay mass gap | unbuilt | **0 %** |

---

## 8 · Sharpest line and reproduce

> **The flat-chart plateau is replaced by a gauge-invariant one: the total-normalized zero mode
> `q=V^{-1/2}Σ_x(Õ−⟨Õ⟩)` of an APE-smeared `tr F²` density has `χ_O=Var(q)=Σ_z C_s(z)`, which is
> exactly gauge-invariant (flat chart drifts `40–1000 %` on a config; `χ_O` is stable to `10^{-15}`),
> non-self-averaging (`χ_O=O(1)`; `Var(q̄)=χ_O/V↓ 1.3→0.1×10⁻³`), and acid-clean (free massive `1/m²`
> known, free massless `∞`). On β=2.5 SU(2) it plateaus sub-linearly — `0.328,0.377,0.411` at L=4,6,8,
> ratios `1.15±0.09, 1.09±0.15` against the free-massless `2.25,1.78` (growth excluded ≈12σ, ≈4.6σ).
> `χ_P` dropped because the Polyakov loop is the wrong (N-ality-1, centre-string) channel; `χ_Landau`
> grew because a gauge-fixed coordinate accumulates the unphysical longitudinal/Gribov sector (spread
> `1.3→13 %`). This is H1 on the correct carrier — the **ceiling**; the atom `f₁²≥r₀` remains `[O]`.
> No QED. Clay 0%.**

```
py -3.12 -I research/v3/w7/wg_invariant.py         # L=4,6,8  (~95 s, one process)
py -3.12 -I research/v3/w7/wg_invariant.py quick    # L=4,6 tiny N (~15 s)
# -> research/v3/w7/wg_invariant_results.json
```

*Refs.* **Carrier / normalization:** `T4-CARRIER.md` (T4-1 `Var(q)=χ_O` vs `Var(q̄)=χ_O/V`, T4-3 squeeze,
T4-4 no-go, UV contact), `Q6-H1-STATUS.md` §1 (normalization load-bearing), §4.5 (carrier hazard),
`H1-IS-IT-THE-GAP.md` (`:φ²:` atomless acid, squeeze `f₁²/χ≤m₁²≤M₀/χ`, inverse-susceptibility),
`docs/paper-v2/ch02-squeeze.tex` (`eq:h1`, `eq:squeeze`, `eq:hstar`, thin-atom Prop., `§schur`),
`docs/paper-v2/app.tex` (`app:postaudit`), `research/unconditional/H1-NUMERIC-RIGHT.md` (composite
`q̄` vacuous — the artifact T4 repairs). **The run it replaces:** `M1-4D-RUN.md` (flat plateau
`3.81/3.76/3.82`, `χ_Landau 37→59→115`, Gribov `1.3/11.2/13.0 %`, glueball unresolved), 
`M1-FAST-CROSSCHECK.md`, `Q5-4D-CHI.md` (gauge systematics, Polyakov contaminated at pilot N),
`Q6-H1-STATUS.md` (carrier, H1 not sufficient). **Polyakov / channel:** `P3-POLYAKOV-CENTER.md`,
`T5-CENTER-VORTEX.md` (N-ality-1 vs scalar, centre-blindness). **Estimator / probe:** `PROOF-GAUSS-INVARIANT-PF.md`
(gauge-invariant transfer), `INVENT-GAUGE-INVARIANT-MASS.md` (spectral edge of `tr F²`, acid: free
Maxwell `m_GI=0`), `W-WAVE-CHARTER.md` (W4 source optimization / smearing). **Standard (textbook, not
re-fetched):** Lüscher gradient-flow renormalization; APE/HYP smearing; Källén–Lehmann; Osterwalder–Seiler
Ann. Phys. **110** (1978) 440; Schur-concavity / convex order. **Witness** `research/v3/w7/wg_invariant.py`,
data `wg_invariant_results.json` (numpy/numba, one process, no downloads). No `[C]` upgraded to `[T]`;
no QED; **Clay 0%.**
