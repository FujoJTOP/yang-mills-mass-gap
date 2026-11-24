# SM6 — the lattice soft-mode scan: softest-mode conditioned curvature at β=2.5, its reanalysis on existing SU(2) data, and the minimal decisive run

> 2026-11-24 · SM6 attack on the **softest mode** of the Yang–Mills UV face — the direction that
> minimizes the conditioned curvature, not a fixed probe direction.
> Reads: `H5-HEAD-NUMERICS.md`, `H4-HEAD-COUNTEREXAMPLE.md`, `H6-HEAD-SPECTRAL.md`;
> `COMPUTE-GPP-K0-CONTINUUM.md`, `COMPUTE-GPP-SU3.md`, `H1-NUMERIC-RIGHT.md`,
> `SAVVIDY-LATTICE-CELL.md`, `K0-OVERLAP-AUDIT.md`, `NUMERICS-AUDIT.md`;
> `docs/paper-v2/ch12-numerics.tex`.
> Witness: `research/v3/w7/sm6_soft_scan.py` (+ `sm6_soft_scan_results.json`; `py -3.12 -I`, numpy/scipy,
> reads only local campaign `.npz` — `partA_su2_N{4000,800}.npz`, `k0_su2_b2.4_L8_n2.npz`; **no downloads**).
>
> **Labels (charter §2 / task instruction).** **[V]** numerical witness on existing data ·
> **[C]** conditional / protocol / design · **[O]** exact missing step. **[T]**/**[NG]** appear **only**
> where inherited verbatim from the corpus. **No QED. Clay 0%.**

---

## 0 · Verdict (no QED)

**A protocol, a reanalysis, and a qualified negative — the soft-mode margin is *not* certified positive.**

1. **[V] The existing SU(2) data hold only fixed transverse plane-wave probes, not a Hessian eigenbasis,
   and no raw links.** The stored per-mode scalars `(a,H,G[,P1,P2])` let us evaluate the conditioned
   curvature `G″(v;a₀*)` on **each** stored direction `v`, but the block Hessian's cross-couplings
   `⟨v_i, Hess v_j⟩` are not stored, so we **cannot diagonalize** to the true softest eigenvector. The best
   available proxy is the **minimum over the sampled plane waves**. By the Rayleigh inequality that
   minimum is an **upper bound on the true softest eigenvalue** (§1.5): a positive sampled minimum does
   *not* certify a positive softest margin; a negative one *does* refute it.
2. **[V] The sampled soft-mode curvature is positive but thin, and softens with volume.**
   SU(2), KDE gold-standard `G″(v;a₀*)` minimum over modes: **+0.200 ± 0.019 (L=4, β=2.5)**,
   **+0.164 ± 0.030 (L=6, β=2.5)**, **+0.142 ± 0.056 (L=8, β=2.4)** — monotone down,
   `L=8` only `2.5σ` above zero. The mean is flat `≈ +0.26…+0.33`; the *softest* direction is the one
   that drifts down. On the `L=8` full-decomposition data the exact-identity softest mode is
   `Gid = +0.080 ± 0.189` — **statistically zero**.
3. **[V]/[O] Keeping only the proven `W_geo ≥ 0`, the soft-mode lower bound is negative.**
   At `L=8` the proven carrier is stably `W_geo = +0.254` (min `+0.243`), but the action bracket + cross
   term is soft/negative (`bracket` min `−0.054`, `X` min `−0.127`); at the identity-softest mode
   `bracket+X = −0.181`, so the rigorous-style bound there is `≈ −0.48` (aggregate `−0.008`). **Zero is
   not excluded** on the softest direction.
4. **[C] All measured numbers are gauge-dependent** — a fixed coordinate plane wave in the **unfixed**
   flat chart, with no Landau fixing, no Faddeev–Popov measure, no Gribov-copy control. The gauge-invariant
   skeleton is only `W_geo ≥ 0` `[T]` and the exact identity `G″ = Δ_dyn + W_geo + X` `[T]`; their
   *measured magnitudes* are chart-dependent.
5. **[NG] Fixed plane waves are structurally blind to the Savvidy soft mode.** The head's soft mode is the
   **localized** negative (magnetic) mode of the background; the stored basis is delocalized `k=2π/L`
   plane waves, nearly orthogonal to it. So the sampled `+0.14…+0.20` is on the **stiff** sector: it is
   not evidence *for* the soft-mode margin — it is evidence that the probes miss the wall.
6. **[C] The minimal decisive run is specified** (§3): gauge-fixed (Landau) + Gribov-copy-averaged,
   **block-Hessian diagonalization** (not fixed plane waves), at `β=2.5`, `L=8,12`, `N≈2000`, one job at a
   time; `L=16` for the acid-test/volume trend. It decides the **statistical** sign of the *true* softest
   eigenmode's conditioned curvature with a gauge systematic; it does **not** turn MC into a uniform
   bound — that is H1 `[O]`.

---

## 1 · The estimator — softest mode, conditioning, error budget `[C]`

### 1.1 The object: the softest mode of the block Hessian

Fix the head coupling (`β=2.5`, SU(2)-native) and a block `B` at scale `ℓ_b≈2a` with boundary `U′`. Write
the block measure as `∏_{ℓ∈B}dU_ℓ = e^{−V(X)}dX`, `V = S_W + Φ`, `Φ=−log J_Haar` (in a gauge-fixed
realization `Φ` absorbs the gauge-fixing + Faddeev–Popov Jacobian; §1.3).

For a unit **coarse direction** `v` (linear / divergence-free, `∇·v=0` inside the block) and the coarse
coordinate `a := ⟨v, X⟩`, the **gauge-invariant conditioned coarse curvature** at the coarse minimiser
`a₀*` is the exact identity (H5 §1.1; `[T]` in corpus):

```
G″(v;a₀*) := −d²/da² log ρ_v(a)|_{a₀*} = E_{μ(·|U′)}[∂_v²V|a₀*] − Var_{μ(·|U′)}(∂_vV|a₀*)
           = Δ_dyn + W_geo + X ,      (1.1)
Δ_dyn = E[H|a₀*] − Var(G|a₀*) ,  W_geo = E[Φ″|a₀*] − Var(Φ′|a₀*) ≥ 0 [T] ,  X = −2Cov(G,Φ′|a₀*) ,
```
with `H = ∂_v²S_W`, `G = ∂_vS_W`, `Φ′=∂_vΦ`, `Φ″=∂_v²Φ`, all conditioned at `a=a₀*`.

**The softest mode** is the minimizer of the Rayleigh quotient of the block Hessian of `V` over the
admissible (linear divergence-free) directions,

```
v* := argmin_{|v|=1, ∇·v=0}  G″(v;a₀*) ,     margin := G″(v*;a₀*) .                (1.2)
```

Equivalently the smallest nonzero eigenvalue direction of the block Hessian, with the constant (gauge)
mode projected out. The *head margin* is whether `margin ≥ δ_W > 0` uniformly `[O]` (H1); the *soft-mode
scan* asks only for the sign and trend of `margin` at `β=2.5` at fixed sampled basis.

**Why the constant mode is not the soft mode.** On a small torus the `k=0` mode of the **gluon** is pure
gauge (projected out), and the `k=0` mode of the **composite** `tr F²` is a *susceptibility*
(`Var ∝ 1/V`, H1-NUMERIC-RIGHT §2–3) — `G″=1/Var → ∞`, vacuous. So the physical softest mode is the
lowest **nonzero** transverse block mode, which is what §2 samples and §3 diagonalizes.

### 1.2 Conditioning and evaluation point

- `a₀* = argmax_a ρ_v(a)` — the mode of the coarse marginal (= coarse minimiser of `S_eff=−log ρ_v`).
  On the symmetric (Haar/Wilson) measure `ρ_v` is even, so the true `a₀*=0`; the witness recovers
  `a₀* ∈ [−0.12, +0.02]` (KDE-mode noise on a broad `sd≈1.9` even distribution). We evaluate `G″` **at
  `a₀*`, not uniformly** — the strong form uniform in `a₀` is false `[NG]` (the `B=−I`, `a₀=0` barrier),
  while `G″(a₀*)≥0` is `[T]` and its *uniform lower bound* is the head.
- Conditional moments by **Nadaraya–Watson kernel at `a₀*`** (`bw=0.30·sd`) for `E[H],Var(G),E[Φ″],Var(Φ′),
  Cov(G,Φ′)`; gold-standard `G″` by the curvature of `−log P(a)` (Gaussian KDE `bw=0.20`, quadratic fit
  over `|a−a₀*|<0.8·sd`).

### 1.3 Gauge fixing and Gribov-copy control (the genuinely new ingredient)

The flat-chart `a=⟨v,log U⟩` with a **fixed coordinate** `v` is **gauge-dependent** (H5 §1.4;
K0-OVERLAP-AUDIT; H1-NUMERIC-RIGHT). A gauge-invariant realization requires:

- **Landau gauge**, FFT/over-relaxation, with **`N_copy = 5–10` random starts per config**; the **spread of
  `G″(v*;a₀*)` across copies** is the operational gauge-invariance error bar (Gribov-copy systematic).
- Fold the gauge-fixing + **Faddeev–Popov** Jacobian into `Φ` (stochastic FP determinant, or a controlled
  bound on its effect). This is the one genuinely new numerical object relative to the existing pipeline.
- The linear/divergence-free projection `P` (`v=Pv`, `P` constant in chart) makes the projection defect
  **exactly zero** (`R_proj=0`, H5 §1.3); the covariant projection would add a sign-indefinite `O(W_geo)`
  term `[NG]` and must be avoided at the head.

### 1.4 Error budget

| source | estimate | how controlled |
|---|---|---|
| statistical (`τ_int≈2`) | delete-block jackknife, 20 blocks (block ≫ `τ_int`); `±0.007–0.06` on `G″` | §2 witness + protocol |
| object/closure (KDE vs `E[H]−Var(G)+W` identity) | **`±0.08–0.12` — dominant** (H1-NUMERIC-RIGHT §1.1) | report both estimators, never one |
| gauge / Gribov-copy | **new run** — spread over `N_copy` copies | §3 |
| KDE bandwidth | `G″` stable to `±0.02` over `bw=0.12–0.35` | corpus |
| min (non-smooth) | bootstrap spread of the sampled min, `B=80–120` | witness (`se_boot_min`) |

**The statistical SE is not the bottleneck; the object/closure systemic and the (absent) gauge systematics
are.** A positive central value with `±0.1` closure and no Gribov control is a *detection on a proxy*, not
a bound.

### 1.5 Two rigorous one-sided facts about sampling `[C]`

- **Rayleigh upper bound.** For the symmetric block Hessian `K`, `λ_min(K)=min_{|v|=1}vᵀKv ≤ vᵀKv` for
  **every** sampled `v`. Hence `min_{sampled} G″ ≥ λ_min`: the sampled minimum is an **upper** bound on
  the true softest conditioned curvature. *A positive sampled minimum cannot certify positivity of the
  soft one; a negative sampled value certainly refutes it.* (And `min` over more directions only lowers
  the bound.)
- **Plane-wave blindness `[NG]`.** The stored basis is delocalized `k=2π/L` plane waves; the Savvidy soft
  mode is the **localized** lowest-mode of the background. The two have small overlap. So the stored
  minimum probes the **stiff** directions, not the wall — a *no-go for the fixed-plane-wave basis*, and the
  deepest reason the existing data cannot answer the soft-mode question.

---

## 2 · Reanalysis of existing SU(2) data at β=2.4–2.5 `[V]`

### 2.1 Data inventory (what actually exists)

| file | group/β | L | modes | stored | decomposition |
|---|---|---|---|---|---|
| `partA_su2_N4000.npz` | SU(2) **2.5** | 4 | 48 | `a,H,G` | no `Φ′/Φ″` (no `W_geo/X` split) |
| `partA_su2_N800.npz` | SU(2) **2.5** | 6 | 48 | `a,H,G` | no split |
| `k0_su2_b2.4_L8_n2.npz` | SU(2) **2.4** | 8 | 12 (6@`π/4`,6@`π/2`) | `a,H,G,P1,P2` | **full** (`Δ+W_geo+X`) |

The 48 modes are the off-diagonal transverse plane waves at `k=2π/L` (`μ≠ν`, `cos/sin`, two colour
components). **There are no raw SU(2) β≈2.5 link configurations in the repository** (the only raw configs,
`probe-center/runs_L8_b2.6*`, are 3×3 `SU(3)` at strong coupling), so **no gauge fixing can be done on
existing data** — the softest true eigenmode is unreachable without a new run.

### 2.2 The soft-mode table — minimum over sampled directions, and the `L`-trend

Witness `sm6_soft_scan.py` (`B=80` config bootstrap for the min; `N` configs as stated). `κ ≡ G″(a₀*)` KDE.

| L | β | modes | `⟨a₀*⟩` | **`κ` mean** | **`κ` min** | bootSE(min) | `κ` median | `1/Var` | identity `Gid` min |
|---|---|---|---|---|---|---|---|---|---|
| **4** | 2.5 | 48 | +0.021 | +0.268 | **+0.200** | 0.019 | +0.262 | 0.269 | — (no `Φ`) |
| **6** | 2.5 | 48 | −0.088 | +0.327 | **+0.164** | 0.030 | +0.312 | 0.272 | — (no `Φ`) |
| **8** | 2.4 | 12 | −0.117 | +0.264 | **+0.142** | 0.056 | +0.267 | 0.276 | **+0.080** (0.189) |

**Trend.** The *mean* soft-mode curvature is flat, `≈ +0.26…+0.33`; the **minimum softens with volume**:
`+0.200 → +0.164 → +0.142` (L=4,6,8), with the error on the min growing (fewer modes, more volume). At
`L=8` the KDE min is `2.5σ` from zero; the exact-identity min is `+0.080 ± 0.189` — **consistent with
zero**. The `1/Var` proxy (`≈0.27`, `L`-stable) matches the mean, not the min — the Gaussian-proxy
identity `G″≈1/Var` holds for the *typical* direction, not the soft tail.

### 2.3 The full decomposition at `L=8` — which term carries the soft mode

Per-mode (witness; sorted softest first), `b2.4, L=8`:

| mode | k | `κ`(KDE) | `bracket` | `W_geo` | `X` | `Gid` | `bracket+X` |
|---|---|---|---|---|---|---|---|
| `ax3_p1_c2_n1` | π/4 | **+0.142** | +0.182 | +0.247 | +0.049 | +0.477 | +0.230 |
| `ax2_p3_c2_n1` | π/4 | +0.173 | +0.049 | +0.257 | +0.049 | +0.355 | +0.098 |
| `ax1_p3_c2_n2` | π/2 | +0.204 | −0.002 | +0.256 | +0.062 | +0.317 | +0.061 |
| … | | | | | | | |
| `ax2_p3_c2_n2` | π/2 | +0.265 | −0.037 | +0.243 | +0.036 | +0.241 | −0.001 |
| `ax1_p2_c2_n2` | π/2 | +0.391 | **−0.054** | **+0.261** | **−0.127** | **+0.080** | **−0.181** |

Aggregates: `bracket +0.089` (min `−0.054`), `W_geo +0.254` (min `+0.243`, max `+0.264` — **the stable,
proven-`≥0` carrier**), `X −0.001` (min `−0.127`). The two "softest" notions disagree by the known
`±0.1` closure:
- **softest by KDE `κ`** (`ax3_p1_c2_n1`, `+0.142`) has `bracket+X=+0.230` — *stiff* in the identity;
- **softest by identity `Gid`** (`ax1_p2_c2_n2`, `+0.080`) has `bracket+X=−0.181` — its small positive
  `Gid` is **entirely** `W_geo`, the only proven-`≥0` term.

**Rigorous-style lower bounds (keep only `W_geo≥0`, treat `bracket+X` statistically):**
`LB_rig(aggregate) = +0.087 − 1.645·0.058 = −0.008`;
`LB_rig(softest identity mode) = −0.181 − 1.645·0.183 ≈ −0.48`.
On the softest direction the proven term does not carry the sign; **0 is not excluded**, exactly H5's
aggregate `−0.042`, sharpened here per direction.

### 2.4 What is gauge-dependent, and what is not

- **Gauge-DEPENDENT (all the measured numbers):** `a`, `κ`, `bracket`, `X`, the *value* of `W_geo`, `Gid`.
  A fixed coordinate plane wave `v` in the **unfixed** flat chart changes under `U→ΩUΩ†`; no FP measure,
  no orbit/Gribov quotient is present. These are numbers in *one implicit gauge*.
- **NOT gauge-dependent (the skeleton):**
  - the **identity** `G″ = Δ_dyn+W_geo+X` (pure differentiation, exact in any chart) `[T]`;
  - **`W_geo ≥ 0` pointwise** for every mode, `k`, `L` (Φ convex ⇒ Prékopa; `Φ″≥N/12`) `[T]` — the only
    *certified* positivity, value-free;
  - the **free-field value** `4sin²(k/2)` at a given `k` (acid test) `[T]`.
  The *magnitude* `W_geo≈+0.25` beyond 0 is **not** certified — that gap is H1 `[O]`.

### 2.5 Acid test on the softest sampled mode

Free massless value at the same `k`: `4sin²(π/L/…)` → `2.0 (L4)`, `1.0 (L6)`, `0.586 (L8)`.
Measured soft-mode minimum / free = **`0.10` (L4)**, **`0.16` (L6)**, **`0.24` (L8)** — all **`< 1`**.
The interacting softest curvature is *below* the free value at every volume: the soft-mode margin **does
not pass the finite-`k` acid test**. The ratio climbs with `L` only because the free box gap `∝k²` shrinks;
whether it reaches 1 at `L≈16` (where `4sin²(π/16)≈0.15 < W_geo`) is **not** established — the measured
`G″` also drifts down as `k` decreases (COMPUTE-GPP-K0-CONTINUUM: `G″(π/4)<G″(π/2)`). `[V]`

### 2.6 Honest limits of the reanalysis

- **Not diagonalized** — sampled minimum only; cross-couplings absent; true softest unknown (§1.5).
- **Not gauge-invariant** — unfixed flat chart (§2.4).
- **Blind to the localized soft mode** — plane waves vs the Savvidy lowest mode (§1.5).
- **`L=8` is `β=2.4`, not `2.5`**, and only 12 modes / `N=400`; the `L`-series mixes a small `β` change.
- **Not a continuum or `k→0` statement** — `ξ/a≲0.3`, single `β` per group, `k≥π/4`.
- **Dominant systematic** is the object/closure `±0.1`, larger than the statistical SE.

---

## 3 · The minimal decisive run `[C]`

Goal: decide, **statistically and with gauge control**, whether the **true softest block-Hessian
eigenmode** has `G″(v*;a₀*)>0` at `β=2.5` — the object the existing data cannot form.

### 3.1 Protocol

- **Group / coupling.** SU(2) primary (`β=2.5`, head-native). SU(3) `β≈5.7` (string-tension-matched) as a
  cross-group check.
- **Volumes.** `L = 8, 12` for the decisive sign + Gribov control; add `L=16` only for the acid-test/volume
  trend (the first volume where `4sin²(k_min)` drops below the `W_geo≈0.25` floor). One job at a time
  (charter §4: **never parallel** SU(N) `L≥8`; the prior OOM).
- **Configurations.** `N = 2000–4000` after `500` thermalisation sweeps, stored every `gap=5–10` sweeps
  (`τ_int≈2–3` stored configs ⇒ delete-block bins ≫ `τ_int`). Streaming sampler, configs discarded,
  peak RSS ≲ 200 MB at `L=12` SU(2).
- **Gauge fixing.** Landau (FFT/over-relaxation), **`N_copy=5–10` random starts per config**; the **spread
  of `G″(v*;a₀*)` across copies is the gauge systematic**. FP determinant folded into `Φ` (stochastic
  estimate or controlled bound).
- **The soft mode (the new ingredient).** Per config, on a `b=2` block, form the **block Hessian**
  `H̄ = E[Pᵀ Hess(S_W+Φ) P]` on the **linear divergence-free** subspace, **diagonalize** (dense `eigh`,
  `D≈` tens), take the lowest eigenvector `v*`; measure the exact marginal of `V_full=P v*` (KDE + kernel
  moments) on the **gauge-fixed** ensemble, at the block scale `k_b=2π/ℓ_b` (not `k→0`). This is the step
  that *finds the soft direction* the plane-wave data miss (§1.5).
- **Estimator.** (1.1) at `a₀*`; store `H,G,Φ′,Φ″` per mode so both bounds are formable:
  `LB_stat = G″ − 1.645·SE`, `LB_rig = (Δ_dyn+X) − 1.645·SE` (keep only `W_geo≥0`). Combine Gribov spread
  in quadrature with the delete-block jackknife SE and the KDE-bandwidth systematic.
- **Target.** Report `G″(v*;a₀*)`, `LB_stat`, `LB_rig`, and the **Gribov spread**, per volume, and the
  `L=8→12` trend. Tested floor `δ_W=0.1` (a conservative fraction of the measured `≈0.25`). **Honest
  expectation:** `LB_stat` may be `>0` but thin (the sampled min is already only `2.5σ` at `L=8`);
  `LB_rig` will likely stay `<0`, because `W_geo`'s magnitude beyond 0 is uncertified. The run **measures
  the correct object with gauge control**; it does not become a proof.

### 3.2 Cost estimate (grounded in measured runtimes)

Anchor: `k0` SU(2) `L=8`, `N=400`, 12 modes = **364 s** ⇒ `0.076 s` per `(mode·config)` at `L=8`
(H5 §3.2). The soft-mode run evaluates `≈D≈24–48` directions/config (to fill and apply `H̄`), then
diagonalizes a `D×D` matrix (negligible).

| volume | sites | measurement (`≈36` dirs, `N=2000`) | gauge-fix ×8 copies | subtotal (py) |
|---|---|---|---|---|
| L=8 | 4096 | `0.076·36·2000 ≈ 5.5 ks ≈ 1.5 h` | `≈0.5 h` | **≈2 h** |
| L=12 | 20736 | `×5 ≈ 7.7 h` | `≈2.5 h` | **≈10 h** |

**SU(2), `L=8,12`, one block mode, Gribov-averaged: ≈12 CPU-hours pure-Python/numpy; ≈1–2 h on C/OpenMP
16 cores; one job at a time.** `L=16` adds `≈30 h` py (`≈4–6 h` C). SU(3) cross-check `≈2–3×`. The
**decisive cost is systematics (Gribov + FP), not statistics**: `N=2000` already drives the statistical
SE below the `±0.1` object floor. A cheap **feasibility pilot** (`L=8`, `N=500`, gauge-fixed,
diagonalize once) at `≈0.5 h` py tells whether any negative eigenvector appears and bounds the Gribov
spread before committing.

### 3.3 What this run decides — and does not

- **Decides (statistically):** the sign and volume-trend of the **true softest** conditioned curvature at
  `β=2.5`, with a controlled Gribov/gauge systematic and the full `Δ_dyn+W_geo+X` decomposition. It
  upgrades "evidence on fixed stiff probes in an unfixed gauge" to "a controlled statistical statement on
  the soft mode in Landau gauge". `[C]`
- **Does not decide (rigorously):** a **uniform** `δ_W>0` independent of `(L,U′,v)`. Keeping only the
  proven `W_geo≥0` leaves `LB_rig<0`; closing it is `(◆)=Var(G|a₀*)−E[H|a₀*]<W_full` = H1 `[O]`. No
  finite lattice dataset supplies a rigorous enclosure (PROOF-VALIDATED-GPP-SIGN §3). A **protocol** is
  the correct deliverable; a **proof** is the wall.
- **Structural risk.** If the soft direction is the localized Savvidy mode, even this run reports its
  `G″` as a *typical-point* number; whether the **essential infimum** of `V_eff''` is positive (the
  spectral/H1 content, H6) is a **separate** question this run informs but does not settle.

### 3.4 Cheapest gauge-invariant cross-check (no Gribov) `[C]`

In parallel, measure the **gauge-invariant per-block effective potential** (H5 §1.4(B)): a blocked
plaquette / Polyakov-loop phase in a representation, pooled over blocks×configs, curvature of
`−log P_block(a)` at its mode. A per-block coarse spin has `O(1)` variance ⇒ **non-vacuous** (unlike the
volume-averaged `0⁺⁺` susceptibility, H1-NUMERIC-RIGHT). This needs **no gauge fixing** and directly tests
a gauge-invariant soft-mode margin; cost `≈0.3–0.6 h` py at `L=8`. It is the natural low-cost companion to
the Landau+Gribov run.

---

## 4 · Acid test, label ledger, sharpest line

### 4.1 Acid test (free field must fail)

Free massless (`Φ=0`, Gaussian `V`): `W_geo=0`, `Δ_dyn=0`, `G″=0` (no floor manufactured). At finite `k`
the free value is `4sin²(k/2)>0`; the measured soft-mode `G″` is **below** it at every `L`
(`G″min/free = 0.10, 0.16, 0.24`) — the probe is *not* a gap signal at finite `k`. The `12 %` `W_geo`
floor `[T]` is a *measure* fact the free field beats, not an IR-mass excess. A positive `δ_W` would require
the interacting `G″` to *exceed* free at `k→0` — not seen. `[T]+[V]`

### 4.2 Label ledger

| item | statement | label |
|---|---|---|
| softest mode = `argmin_{|v|=1,∇·v=0} G″(v;a₀*)`; conditioned curvature identity incl. `X` | (1.1)/(1.2) | **[T]** (identity) + **[C]** (estimator) |
| evaluation at coarse minimiser `a₀*=argmax ρ_v`; symmetric ⇒ 0; witness `a₀*∈[−0.12,+0.02]` | H4/V3 | **[T]**/**[V]** (inherited) |
| sampled minimum ≥ true softest eigenvalue (Rayleigh) | §1.5 | **[C]** |
| fixed plane waves blind to the localized Savvidy soft mode | §1.5 | **[NG]** |
| soft-mode KDE `G″` minimum: `+0.200/+0.164/+0.142` at L=4/6/8, bootSE `0.019/0.030/0.056` | §2.2 witness | **[V]** |
| `L=8` identity softest mode `Gid=+0.080±0.189` (zero-consistent) | §2.3 | **[V]** |
| proven carrier `W_geo=+0.254` (min `+0.243`) stable; `bracket+X` soft (`−0.181` at softest) | §2.3 | **[V]** |
| rigorous-style soft-mode lower bound `LB_rig(aggregate)=−0.008`, `LB_rig(softest)≈−0.48 < 0` | §2.3 | **[V]** |
| measured numbers gauge-dependent (unfixed flat chart, no FP, no Gribov) | §2.4 | **[C]** (cannot certify) |
| `W_geo≥0` pointwise, value-free; identity exact; free value `4sin²(k/2)` | §2.4 | **[T]** |
| soft-mode `G″`/free `<1` at L=4,6,8 (acid-test fail at finite k) | §2.5 | **[V]** |
| minimal decisive run: Landau + Gribov + block-Hessian diagonalization, β=2.5, L=8,12, N≈2000 | §3; ≈12 CPU-h py / 1–2 h C | **[C]** |
| gauge-invariant per-block cross-check (no Gribov), `O(1)` variance non-vacuous | §3.4 | **[C]** |
| uniform `δ_W>0` (soft-mode margin), `(★_RG)` = H1 | §3.3 | **[O]** |
| a lattice run proves the head | **no** — statistical only; rigorous is H1 | **[C]**/**[O]** |
| Clay mass gap | open | **0%** |

### 4.3 Sharpest line

**The existing SU(2) data can only sample the soft mode along fixed transverse plane waves in an unfixed
gauge — the wrong basis (blind to the localized Savvidy mode) on the wrong object (gauge-dependent). On
that sampled set the soft-mode conditioned curvature is positive but thin and **softens with volume**
(`+0.20 → +0.16 → +0.14` at `L=4,6,8`; `L=8` only `2.5σ`, identity `+0.08±0.19` zero-consistent), its only
proven carrier is the value-free `W_geo≥0`, and the rigorous-style soft-mode lower bound is negative
(`−0.48` at the identity-softest direction): **the soft-mode margin is not certified positive.** The
minimal decisive run is a gauge-fixed (Landau) + Gribov-copy-averaged **block-Hessian diagonalization** —
which finds the true soft eigenvector the plane waves miss — at `β=2.5`, `L=8,12`, `N≈2000` (≈12 CPU-h py,
one job at a time), with a no-gauge-fixing per-block effective-potential cross-check; it can turn the
sampled stiff-sector evidence into a controlled **statistical** statement on the soft mode, and nothing
more — the uniform floor `δ_W>0` is H1, not a dataset.** **No QED. Clay 0%.**

---

*Refs.* `H5-HEAD-NUMERICS.md` §1–3 (estimator, `Δ+W_geo+X`, `LB_rig`); `H4-HEAD-COUNTEREXAMPLE.md` §1
(`a₀*`, `G″(a₀*)=−ρ″/ρ`); `H6-HEAD-SPECTRAL.md` (essential infimum vs typical point, Kramers);
`COMPUTE-GPP-K0-CONTINUUM.md` §1–5 (`G″=bracket+W`, `W_geo`, `LB_rig`, volume confound);
`COMPUTE-GPP-SU3.md` (W-dominance, `1/Var`); `H1-NUMERIC-RIGHT.md` §2–3 (composite susceptibility,
`Var∝1/V`, vacuous; gluon acid-test fail); `SAVVIDY-LATTICE-CELL.md` (per-mode `G″` spread, negatives);
`K0-OVERLAP-AUDIT.md` (gluon vs `0⁺⁺`); `NUMERICS-AUDIT.md` (finite-k acid test); `docs/paper-v2/ch12-numerics.tex`
§decisive-test. Witness `research/v3/w7/sm6_soft_scan.py` (+ `sm6_soft_scan_results.json`). No third-party
sources downloaded this session.
