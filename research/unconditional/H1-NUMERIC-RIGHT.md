# H1-NUMERIC-RIGHT — measuring `G″` on the CORRECT object: the gauge-invariant `0⁺⁺` soft mode is a **susceptibility** (`Var→1/V`), and the only `Var`-nonvanishing carrier is the gauge-dependent gluon, which fails the acid test at every `β`

> 2026-11-23 · SU(2) lattice, `L≤8`, streaming single process, one job at a time (charter §4). This is the correction of the `COMPUTE-GPP-K0-CONTINUUM` measurement that `K0-OVERLAP-AUDIT` proved was the **wrong object** (a gauge-dependent transverse gluon, colour-adjoint `1⁻`, not the gauge-invariant `0⁺⁺`). It measures the object the audit named — the gauge-invariant `tr F²` density — and a `β`-scan of the gluon with the **acid test done side-by-side**. Code: `research/unconditional/h1_numeric_right.py`; data: `research/unconditional/h1_right_L8_N200.json` (+ `h1_right_L8_N200.log`).

**Verdict (three-line).**
1. **The gauge-invariant `0⁺⁺` soft mode — the "zero-momentum mode of `tr F²`" the request names — is NOT an H1 carrier: it is a susceptibility.** Its amplitude variance scales as `Var(q̄) = χ/V` (measured exponent `V^{−1.01}` at β=2.4, `V^{−0.96…−1.09}` at β=1.5 across `L=4,6,8`), so `G″_dens = 1/Var = V/χ → ∞`. H1's `Var(a_{k→0}) ≤ 1/λ` is then satisfied for **every** `λ` at large `V` — **vacuously**. (A positivity of `G″` in the direction that *grows* with volume is the exact opposite of a gap.)
2. **The only gauge-consistent object with `Var ↛ 0` is the fundamental field mode (the gluon)** — but there is *no linear gauge-invariant function of the connection* (every gauge-invariant is `≥` quadratic `⇒` composite `⇒` susceptibility, leg 1). So H1 on a lattice splits into a **vacuous** gauge-invariant branch and a **gauge-dependent** field branch.
3. **The gluon branch does give `G″>0` statistically** (`+0.37, +0.11, +0.32, +0.22, +0.30` at β=`1.0…2.4`; 4/5 at `≥2.8σ`), but it is **`W_geo`-carried** (`bracket≈0`), **below the free value at every `β`** (`G″/free = 0.63, 0.18, 0.55, 0.38, 0.51 < 1`), and its **rigorous-style one-sided lower bound is negative at every `β`** (`−0.02…−0.21`). So the `+0.3` is confirmed — now across the whole strong-to-intermediate coupling range — as a **Haar/Jacobian (measure) number that the free field satisfies and exceeds**, **not** an IR-mass excess.

**Net:** no lattice `G″` of the correct (`0⁺⁺`, gauge-invariant) object certifies the gap — the object is provably wrong-type. The exact missing step is stated in §5.

---

## 0. What was measured — three objects, labelled

| tag | object | QN / gauge | role in this report |
|---|---|---|---|
| **GLUON** | `a=⟨v, log U⟩`, `v` = transverse plane wave `(ax,pol,color,n=1)`, 6 modes | `1⁻` vector, adjoint, **gauge-dependent** | the **H1-relevant measure-level** object (`K0-OVERLAP-AUDIT (c)(2)` — H1 is an LSI that *includes* the cone); a **proxy**, clearly flagged. Gives the literal `G″/W_geo/D/S` table + acid test. |
| **DENSITY** | `q(x) = 1 − (1/6)Σ_{μ<ν} ½Re tr U_p(x)` `≈ tr F²`; its `k=0` mode `q̄` and softest Fourier amplitude `q̃(k₁)` | `0⁺⁺` scalar, **gauge-invariant** | **the object the request names** ("zero-momentum mode of `tr F²`"). Measured directly — **not** a proxy. Proves it is the wrong type. |
| **FREE** | `4 sin²(k/2)` at `k=π/4` (transverse mode) | free massless | the **acid test** (`NUMERICS-AUDIT (a)`): interacting must *exceed* free to be a gap signal. |

**Object caveat, stated up-front (the request's CRITICAL item).** The measured gluon soft mode is **the gauge-DEPENDENT transverse gluon, NOT the gauge-invariant `0⁺⁺`** — it is the same object the audit flagged, and it is *not* silently re-used as "correct": it appears only as the H1-relevant **contrast**, with the acid test table. The **gauge-invariant `0⁺⁺`** is measured separately (§2) and is *exactly the suggested object* — and it is shown to be structurally wrong (a susceptibility). **Neither branch yields a certifying `G″`.**

---

## 1. GLUON `β`-scan — the requested table (`L=8`, `k=π/4`, 6 transverse modes `n=1`, `N=200`)

The exact identity (`flat chart`, Gaussian-conditioned at `a=0`; `COMPUTE-GPP-K0-CONTINUUM §1`):
`G″ = bracket + W`, `bracket = E[H|0] − Var(G|0)`, `W_geo = E[Φ″|0] − Var(Φ′|0)` (proven `≥0`), `X = −2Cov(G,Φ′|0)`, `W = W_geo+X`. Gold standard `G″` = curvature of `−log P(a)` (KDE). `S = 1/Var(a)`; **corrected defect** `D/S = 1 − G″·Var = 1 − G″/S` (W7 convention `G″=S−D`). Delete-block jackknife (20 blocks) SEs. Free value `4sin²(π/8)=0.586`.

| β | ⟨u□⟩ | **G″** | S/N | bracket | **W_geo** | X | S=1/Var | D/S | free | **G″/free** | LB_stat | **LB_rig** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.242 | **+0.372 ± 0.070** | 5.3 | +0.004 ± 0.014 | +0.254 ± 0.005 | +0.004 ± 0.009 | 0.256 | −0.45 | 0.586 | **0.63** | +0.256 | **−0.021** |
| 1.5 | 0.366 | **+0.106 ± 0.066** | 1.6 | +0.038 ± 0.029 | +0.254 ± 0.005 | +0.004 ± 0.019 | 0.261 | +0.60 | 0.586 | **0.18** | −0.002 | **−0.017** |
| 2.0 | 0.500 | **+0.320 ± 0.080** | 4.0 | +0.017 ± 0.048 | +0.242 ± 0.004 | +0.033 ± 0.028 | 0.266 | −0.20 | 0.586 | **0.55** | +0.188 | **−0.037** |
| 2.2 | 0.568 | **+0.223 ± 0.078** | 2.8 | −0.036 ± 0.076 | +0.249 ± 0.007 | +0.026 ± 0.035 | 0.281 | +0.21 | 0.586 | **0.38** | +0.094 | **−0.145** |
| 2.4 | 0.629 | **+0.298 ± 0.067** | 4.4 | −0.096 ± 0.092 | +0.256 ± 0.004 | +0.023 ± 0.026 | 0.287 | −0.04 | 0.586 | **0.51** | +0.188 | **−0.213** |

`LB_stat = G″ − 1.645·SE` (statistical). `LB_rig = (bracket+X) − 1.645·SE(bracket+X)` — keeping **only the proven `W_geo≥0`** and treating `bracket+X` statistically (`COMPUTE-GPP-K0-CONTINUUM §4`).

**Reading.**
- **`G″>0` is statistically real** at 4/5 β (`+2.8σ…+5.3σ`); inverse-variance mean `+0.258 ± 0.032` (`8.1σ` **stat**). But the identity-closure systematic is `±0.11` (§1.1), which dominates — an honest combined figure is `≈0.26 ± 0.11`, not an 8σ certificate.
- **`W_geo` is the carrier and is remarkably `β`-stable**: `+0.254, +0.254, +0.242, +0.249, +0.256` (range `0.014`). Consistent with the proof that `W_geo≥0` is a measure (Prékopa) fact; here it is the *entire* `G″` at strong coupling. `bracket ≈ 0` (`+0.04…−0.10`) over the whole range — equipartition, `β`-independent.
- **Acid test fails everywhere.** `G″/free = 0.63, 0.18, 0.55, 0.38, 0.51` — **always < 1**. The interacting coarse `G″` is *smaller* than the free massless value at the same `k` at every coupling. A gap would make it *larger*. Still no acid-test pass — now shown not to be a single-`β` accident but a `β`-uniform fact.
- **`D/S` is NOT a clean `2%`.** Corrected `D/S = −0.45, +0.60, −0.20, +0.21, −0.04`: at strong coupling the marginal is strongly non-Gaussian (the `±60%` swings); only near `β=2.4` does it approach the W7 anchor (`−0.04`, i.e. `G″` ~4% above `1/Var`). The measured `2%` is a *weak-coupling artifact*, not a `β`-independent constant.
- **`LB_rig < 0` at every β** (`−0.017…−0.213`). The proven half (`W_geo≥0`) is value-free; the sign-carrying magnitude is `W_geo≈0.25`, uncertified. 0 is not excluded — the same obstruction as before, now `β`-uniform.

### 1.1 Identity closure (honest systematic)
`bracket + W_geo + X − G″_KDE` at β=2.4: `−0.096+0.256+0.023−0.298 = −0.115`; at β=1.0: `−0.110`. A **`≈−0.11` closure defect** (the Gaussian-conditional-vs-KDE chart systematic). So `W_geo=0.25` does *not* fully explain `G″≈0.30`; `≈0.11` is unexplained by the decomposition. This bounds any "W carries X% of G″" claim to ~±0.11, confirming `NUMERICS-AUDIT (b)` ("per-mode closure misses `G″` by ~its own size").

---

## 2. DENSITY `0⁺⁺` — the suggested gauge-invariant object is a **susceptibility**

`q(x) = 1 − (1/6)Σ ½Re tr U_p ≈ tr F²/4` (gauge-invariant, `0⁺⁺`). `q̄ = (1/V)Σ_x q(x)` is **exactly the "zero-momentum mode of `tr F²`"**; `q̃(k₁) = (1/V)Σ_x q(x)e^{ik₁x}` is the softest nonzero amplitude (`k₁=2π/L`). `N=200`, `L=8`.

| β | Var(q̄) | χ=V·Var(q̄) | Var(q̃(k₁)) | **G″_dens (KDE)** | D/S_dens |
|---|---|---|---|---|---|
| 1.0 | 1.05×10⁻⁵ | 0.043 | 5.19×10⁻⁶ | +2.56×10⁵ | −0.33 |
| 1.5 | 9.29×10⁻⁶ | 0.038 | 5.18×10⁻⁶ | +1.85×10⁵ | +0.04 |
| 2.0 | 1.30×10⁻⁵ | 0.053 | 6.62×10⁻⁶ | +1.20×10⁵ | +0.21 |
| 2.2 | 1.35×10⁻⁵ | 0.055 | 5.50×10⁻⁶ | +1.48×10⁵ | +0.19 |
| 2.4 | 1.06×10⁻⁵ | 0.043 | 4.31×10⁻⁶ | +1.81×10⁵ | +0.22 |

`G″_dens` (KDE curvature of `−log P(q̃)`) is **`~10⁵`** — positive at overwhelming significance, and **meaningless for a gap**: it sits at `1/Var` of a quantity whose `Var→0` (§3); `D/S_dens` small-to-moderate (`−33%…+22%`, mostly `|D/S|≲22%`) ⇒ the marginal is **near-Gaussian** (CLT over `V` sites), i.e. a pure susceptibility fluctuation. `χ = V·Var(q̄) ≈ 0.04–0.05` is **finite and `β`-insensitive** — as expected for a dim-4 operator's susceptibility *both* free and interacting (`∫d⁴x ⟨qq⟩_c` IR-convergent for a dim-4 bilinear), so `χ` is **gap-blind** — exactly the "free-field-satisfied, value-free" failure mode of `W_geo`.

### 2.1 The `1/V` volume test (`L=4,6,8`; density-only, `N=200`)

| β | L (V) | Var(q̄) | χ=V·Var(q̄) | G″_dens=1/Var(q̃) |
|---|---|---|---|---|
| 2.4 | 4 (256) | 1.258×10⁻⁴ | 0.0322 | 1.85×10⁴ |
| 2.4 | 6 (1296) | 2.440×10⁻⁵ | 0.0316 | 1.06×10⁵ |
| 2.4 | 8 (4096) | 1.057×10⁻⁵ | 0.0433 | 1.81×10⁵ |
| 1.5 | 4 (256) | 1.542×10⁻⁴ | 0.0395 | 1.21×10⁴ |
| 1.5 | 6 (1296) | 3.259×10⁻⁵ | 0.0422 | 5.19×10⁴ |
| 1.5 | 8 (4096) | 9.289×10⁻⁶ | 0.0380 | 1.85×10⁵ |

**Measured scaling exponents** (`L4→L6`, `L6→L8`; `Var ∝ V^p`):

| β | Var(q̄) exponent | G″_dens exponent |
|---|---|---|
| 2.4 | V^−1.01 , V^−0.73 | V^+1.08 , V^+0.47 |
| 1.5 | V^−0.96 , V^−1.09 | V^+0.90 , V^+1.11 |

**β=1.5 is textbook:** `Var(q̄) ∝ V^{−1.0}`, `χ` flat (`0.039→0.042→0.038`, ±10%), `G″_dens ∝ V^{+1.0}`. The gauge-invariant `0⁺⁺` soft mode is a **susceptibility**: its variance vanishes like `1/V`, so `Var(a_{k→0}) ≤ 1/λ` holds for every `λ` at large `V` — **H1 is vacuous on this object** — and `G″_dens` diverges.

---

## 3. The structural reason (why the correct object cannot carry the gap)

For **any local operator** `O(x)`, its Fourier amplitude at fixed `k` has
`Var(Õ(k)) = C̃_O(k)/V → 0` (CLT / extensivity; `C̃_O` the connected `2`-point FT). Measured here: `Var(q̄) ∝ V^{−1}`. Consequently H1's `Var_μ(a_{k→0}) ≤ 1/λ` with `a` a **composite** is **trivially true for every `λ`** at large `V` — it carries **no gap content**. Equivalently `G″ = 1/Var → ∞` (measured `+1/V`, §2.1), the *opposite* of the bounded-below inverse propagator `k̂²+m²` H1 needs.

The only amplitudes with `Var ↛ 0` are **fundamental field modes** (`Var = 1/(k̂²+m²)`). But on the lattice the gauge field is `A_μ^a` in `su(2)`; **no linear gauge-invariant function of `A` exists** (every gauge-invariant is a trace / Wilson loop / `≥` quadratic). Hence:

> **H1's `a_{k→0}` is necessarily either (i) a composite `⇒` `Var→1/V` `⇒` vacuous, or (ii) the fundamental gauge field `⇒` gauge-dependent `⇒` not invariant.** There is no third option.

This is why the earlier measurement landed on the gluon: it is the only type with the right `Var` behaviour, and it is exactly the one the audit (correctly) refused to call `0⁺⁺`. The gauge-invariant `0⁺⁺` content lives entirely in the **two-point function** (glueball mass), a *decay-rate* statement, not a `G″`-curvature sign — `K0-OVERLAP-AUDIT (c)(3)`, now with both branches measured.

---

## 4. The `0⁺⁺` glueball mass — measured, but **not resolvable at `L=8`** here

From the connected correlator `C(r)=⟨q(x)q(x+r)⟩_c` (axis 0): the exponential/cosh fit returns `m₀ = 3.53, 2.91, 3.25` (β=2.0, 2.2, 2.4) and `nan` (β=1.0, 1.5; fit unstable). `m₀ ≳ π` is **at the lattice cutoff** — i.e. `ξ ≲ 1–2a`: on `L=8`, `β≤2.4`, the `0⁺⁺` correlator is essentially contact-like and the mass is **not resolved**. A genuine `0⁺⁺` mass (SU(2) `~1.7–1.9 GeV`, `ξ/a ≈ 7`) needs `L≈16–24` in the scaling window (`β≈2.4–2.6`) and an `a→0` extrapolation — the standard, expensive computation (`NUMERICS-AUDIT (d)`), **out of scope here and not attempted**. So even the *legitimate* `0⁺⁺` route is only reachable as standard numerical evidence, never as a curvature certificate — and is not what any campaign run measured.

---

## 5. The exact missing step (honest verdict, charter §2)

- **Gauge-invariant `0⁺⁺` branch — VACUOUS, not a gap.** `Var(q̄) ∝ 1/V` (`V^{−1.0}`), `G″_dens = 1/Var ∝ V → ∞`. H1's `Var(a_{k→0})≤1/λ` is trivially true for all `λ` at large `V`; the positivity of `G″_dens` is the *wrong direction* of scaling. **Numerically demonstrated, `L=4,6,8`, two `β`.**
- **Field (gluon) branch — right type, but gauge-dependent and acid-test-failing.** `G″>0` (`4/5` at `≥2.8σ`, mean `+0.26±0.11`), carried by `W_geo` (`bracket≈0`, `β`-stable `W_geo≈0.25`), **below free at every `β`** (`G″/free=0.18–0.63`), rigorous-style `LB_rig<0` at every `β` (`−0.02…−0.21`). Verified: `W_geo(β=2.4)=+0.256` reproduces `COMPUTE-GPP-K0-CONTINUUM`'s SU(2) `+0.254` to `0.002`.
- **The gap is not a `G″` of either object.** Closing H1 requires an inequality that makes a **`Var`-nonvanishing, gauge-invariant, field-mode** amplitude's curvature *exceed its free value* at `k→0`, uniform in `V`. No such linear gauge-invariant object exists (§3). The only gauge-invariant formulations that survive are:
  1. a **gauge-fixed + orbit-quotiented** (Gribov-aware) field-mode analysis with a *proved* quotient lower bound — not achieved (the Gribov copy problem is the wall); or
  2. the **glueball `2`-point mass** `m₀ ≥ λ` — true/positive but *evidence not proof* and not a `G″` statement.
- **Sharpest line.** The audit was right and the `β`-scan makes it uniform: the `+0.3` is `W_geo` (a Haar-measure number the free field beats at every accessible `k`); the *correct* gauge-invariant `0⁺⁺` soft mode is a susceptibility whose `G″=1/Var` **diverges** with volume (`Var∝V^{−1}`), making H1 on it **vacuous**. `G″(k→0)≥λ` is **not** a statement any of these objects can carry: composite `⇒` vacuous, field `⇒` gauge-dependent. The missing step is a **gauge-invariant field-mode with `Var↛0`** — which the lattice says does not exist linearly. That is the wall, stated honestly.

---

## 6. Caveats (all load-bearing)

- **Object.** §2 is the **true gauge-invariant `0⁺⁺`** (`tr F²` density) — *not* a proxy. §1 is the **gauge-dependent transverse gluon** (`1⁻`, adjoint) — a **proxy** for the H1-relevant measure-level object, clearly flagged, **not** the `0⁺⁺`. The gluon is never presented as "the correct object."
- **Size.** `L≤8`, SU(2), single volume per `β` (a **`β`-scan**, not an `a→0` continuum limit; `a(β)≈0.1–0.18 fm`, box `≈1.4 fm`); softest `k=π/4`, never `k→0`.
- **Closure.** The gluon decomposition has a `≈−0.11` KDE-vs-conditional systematic (§1.1); "W carries all of `G″`" is `±0.11`-limited.
- **D/S.** Not a `β`-invariant `2%`: `−0.45…+0.60`; the `2%` anchor is a weak-coupling (`β≈2.4`) artifact.
- **Glueball mass.** Unresolved at `L=8` (`m₀≳π`, `ξ≲a`); not a measurement (§4).
- **Statistics.** `N=200` configs/`β`; delete-block jackknife (20 blocks) SEs; SU(2) sampling `τ≈1.5–2.6` stored configs (block ≫ `τ`). `G″` statistical errors `±0.07–0.08`; `W_geo` `±0.004–0.007` (very well measured).
- **Memory.** One process, one `L=8` SU(2) lattice; configs never stored (in-place streaming heatbath, single `U` holder ≈0.5 MB); only `N×6` scalars + an `L`-vector correlator accumulated. **Sampler working set 61 MB** (β=1.0, matches the `~60 MB` target); the process peak reached **110 MB** on later jobs because numpy/scipy (KDE) allocations are not returned to the OS between sequential `β` jobs. Safe (single job, no parallelism), but honestly **above** the `60 MB` target at the process level — the target is met by the streaming sampler itself, not by the scipy-KDE analysis stage.

---

## 7. Reproduce

```
py -3.12 research/unconditional/h1_numeric_right.py \
   --betas 1.0 1.5 2.0 2.2 2.4 --L 8 --N 200 --therm 220 --gap 4 \
   --vol_betas 2.4 1.5 --Nvol 200 \
   --out research/unconditional/h1_right_L8_N200.json
```
Runtime 731 s total (sequential; ~100–134 s/`β`, ~8–22 s/volume-scan). Per-`β` peak RSS: β=1.0 `70 MB`, β≥1.5 process `110 MB` (see §6).

---

*Refs.* `K0-OVERLAP-AUDIT.md` (the gluon is `1⁻`/adjoint, not `0⁺⁺`; (c)(2) H1 includes the cone, (c)(3) `0⁺⁺` is a `2`-pt statement); `NUMERICS-AUDIT.md` ((a) interacting `G″`<free at finite `k`, (b) LB<0, (d) rigorous gap needs volume-uniform lower bound); `COMPUTE-GPP-K0-CONTINUUM.md` (`G″=bracket+W`, `W_geo=E[Φ″]−Var(Φ′)`, SU(2) `W_geo=+0.254`, LB `−0.035`); `CONTINUUM-VARIANCE-DEFECT.md` (`G″=S−D`, `S=1/Var`, `D/S` RG-covariant, denominator `S` = H1); `B4-VARIANCE-MAG.md` (CS sign-blind, `VarG≈2`); `PROOF-JACOBIAN-SAVES` (`W_geo≥0`, Prékopa; `Φ″≥N/12`); `PROOF-IR-H1` (measure-level LSI on `∇_∥`, includes the cone); free value `4sin²(k/2)`. Scripts `research/unconditional/h1_numeric_right.py`; data `h1_right_L8_N200.json`, `h1_right_L8_N200.log`.
