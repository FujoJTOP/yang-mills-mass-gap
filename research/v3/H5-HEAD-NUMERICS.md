# H5-HEAD-NUMERICS — the head inequality `G″(v;a₀*)≥δ_W` at `β=2.5`: the estimator (gauge-invariant conditioned coarse curvature, `X` included), what the existing SU(2)/SU(3) data can and cannot certify, and the minimal decisive run

> 2026-11-24 · H5 attack on the last open piece of (★): the **head** `G″(v;a₀*)≥δ_W` at `β=2.5`
> (STAR2-CLASSIFICATION §3/§5: the head is H1 `[O]`).
> Reads: `STAR2-CLASSIFICATION.md`, `V2-TRUE-BOUND.md`, `V3-COUNTEREXAMPLE-SCOPE.md`,
> `RECONCILIATION-STAR.md`, `G1-BLOCK-IDENTITY.md`, `G3-PROJECTION-TERM.md`, `S1-CONDITIONAL-POINCARE.md`;
> `COMPUTE-GPP-K0-CONTINUUM.md`, `COMPUTE-GPP-SU3.md`, `H1-NUMERIC-RIGHT.md`, `NUMERICS-AUDIT.md`,
> `OTHER-GPP-KZERO.md`, `W7-Gpp0.md`, `K0-OVERLAP-AUDIT.md`, `PROOF-VALIDATED-GPP-SIGN.md`,
> `CONTINUUM-VARIANCE-DEFECT.md`; `docs/paper-v2/ch12-numerics.tex`.
> Witness: `research/v3/w7/h5_head_reanalysis.py` (+ `h5_head_reanalysis_results.json`, `py -3.12`, reads
> only local corpus data — `partA_su2_N{4000,800}.npz`, `k0_su{2,3}_*_L8_n2.npz`; no downloads).
>
> **Labels (charter §2 / task instruction).** **[V]** numerical witness on existing data ·
> **[C]** conditional / protocol / design · **[O]** exact missing step. We add **[T]**/**[NG]** only where
> inherited verbatim from the corpus. **No QED. Clay 0%.**

---

## 0 · Verdict (no QED)

**A protocol, a reanalysis, and a clean negative — not a proof.**

1. **[V] The existing data are uniformly positive on the WRONG object.** At the exact head coupling
   `β=2.5` the (gauge-dependent, flat-chart, unfixed-gauge) transverse gluon amplitude has
   `G″(a₀*)=+0.268±0.007` (SU(2) `L=4`, 48/48 modes, `N=4000`), `+0.327±0.016` (`L=6`, `N=800`), and the
   SU(3) soft-blocked mode `+0.51±0.08`; the full-decomposition points at `β=2.4` (SU2, `L=8`) and
   `β=5.7` (SU3, `L=8`) give `W_geo=+0.254`/`+0.318`, `bracket≈0`, `X≈0`. The coarse minimiser is
   `a₀*≈0` (symmetric measure), so these already sit at the correct *evaluation point* `(★_loc)`.
2. **[V]/[O] But the one-sided bound that keeps only the proven `W_geo≥0` is NEGATIVE.** From the only
   datasets that store the `Φ`-derivative columns (`k0`, `β=2.4/5.7`, `L=8`): `LB_rig=−0.042` (SU2) /
   `−0.095` (SU3) — **0 is not excluded**. At `β=2.5` itself the stored data (`partA`) have no `Φ′/Φ″`,
   so the `W_geo/X` split cannot be formed and **no rigorous-style bound exists at the head coupling at
   all**. The positivity is carried by the *measured magnitude* `W_geo≈+0.25`, whose proven floor is `0`.
3. **[C] No existing number is the gauge-invariant conditioned coarse curvature.** All measured
   amplitudes are `a=⟨v,log U⟩` with a fixed coordinate plane wave in the **unfixed** flat chart: they
   are gauge-dependent and have no orbit/Gribov control, no Faddeev–Popov measure, and no block-scale
   (`ℓ_b`) mode (they are global `k=2π/L` Fourier modes). They cannot certify a gauge-invariant `G″`.
4. **[C] The minimal decisive run is specified** (§3): gauge-fixed (Landau) + Gribov-copy-averaged,
   block-scale softest blocked mode, full `S_W+Φ(FP)` decomposition, at `β=2.5`, `L=8,12,16`,
   `N≈2000–4000`, one job at a time. It converts "evidence on the wrong object" into "a controlled,
   gauge-invariant **statistical** one-sided bound at the head." It does **not** upgrade to a theorem:
   a rigorous uniform bound is H1 `[O]`, not a better dataset.

---

## 1 · The estimator — operator, coarse field, projection, `G″` with `X`

### 1.1 The object and the exact identity

Fix the head shell `j` (bare coupling `β_j=2N/g_j²`, here `β=2.5`) and a **block** `B` at scale
`ℓ_b≈2a` with boundary `U′`. Write the block links in the flat chart `U_ℓ=e^{X_ℓ}`,
`∏_{ℓ∈B}dU_ℓ=e^{−Φ(X)}dX` with `Φ=−log J_Haar` (and, in a gauge-fixed realization, the full
gauge-fixing + Faddeev–Popov Jacobian folded into `Φ`; §1.4). Put `V=S_W+Φ`, conditional measure
`μ_{U′}(dX)=Z^{-1}e^{−V}dX`.

For a unit **coarse direction** `v` and the coarse coordinate `a:=⟨v,X⟩`, the **gauge-invariant
conditioned coarse curvature** at the coarse field `a=a₀` is

```
G″(v;a₀) := −d²/da² log ρ_v(a)|_{a₀} = E_{μ(·|U′)}[∂_v²V | a=a₀] − Var_{μ(·|U′)}(∂_vV | a=a₀) ,   (1.1)
```

with marginal `ρ_v(a)=∫δ(⟨v,X⟩−a)e^{−V}dσ`. This is the object (★_C) of RECONCILIATION / V2 / V3
(exact, pure differentiation; `[T]` in the corpus). Splitting `V=S_W+Φ` and conditioning:

```
G″ = Δ_dyn + W_geo + X ,                                                                    (1.2)
Δ_dyn = E[H_vv|a₀] − Var(G_v|a₀)            (the Wilson "bracket"; H_vv=vᵀHess S_W v, G_v=∂_v S_W)
W_geo = E[Φ″_vv|a₀] − Var(Φ′_v|a₀)         (Haar/Jacobian; Φ″⪰(N/12)I  ⇒  W_geo≥0 pointwise [T])
X     = −2 Cov(G_v, Φ′_v|a₀)                (the cross term; sign-indefinite — the "X term")
```

**`X` must be included and is estimated, not dropped** (RECONCILIATION §1: `W_j=W_geo+X`).

### 1.2 The coarse field `a₀*` — definition and minimiser

`a₀*(U′,v) := argmax_a ρ_v(a)` — the **mode of the coarse marginal**, i.e. the coarse minimiser of
`S_eff(a)=−log ρ_v(a)` (V3-COUNTEREXAMPLE-SCOPE §3–4; RECONCILIATION §3). Evaluate `G″` **at `a₀*`,
not uniformly in `a₀`** — the strong form uniform in `a₀` is false `[NG]` (the `B→−I`, `a₀=0` barrier),
while `G″(v;a₀*)≥0` is trivial `[T]` (the marginal has a maximum there) and its *uniform lower bound*
`δ_W>0` is the head/H1. On a symmetric measure (`U↔U^{-1}`) `ρ_v` is even, so the true `a₀*=0`; the
witness recovers `a₀*∈[−0.12,+0.02]` on the stored modes (§2), consistent with 0 at KDE resolution.

### 1.3 The projection — linear/divergence-free `[T]`

The admissible class is the **linear / divergence-free** (transverse) one (STAR2 §4; G1/G3): `v` is a
transverse plane-wave direction, `div v=0`, `v=Pv` with `P` a fixed linear (constant-in-chart) block
average. For this class the flat-chart Stein identity is exact and the whole object collapses to
`E[H]−Var(G)+W_geo` (constant-`P`), and the covariant (parallel-transport) projection adds exactly one
transport defect `R_proj` (G1 §4, G3 §4),

```
G″_P = 2W_geo(w) + R_proj ,   R_proj = E[⟨(Dw)w,∇Φ⟩−⟨(Dw)w,∇S⟩] + Cov(⟨w,∇Φ⟩,div w) − Cov(⟨w,∇S⟩,div w) ,
```

which is `O(g_j)`-small in the UV tower `[C]` but **`O(W_geo)`, sign-indefinite at the head** `[NG]` for
"controlled by `W_geo` alone." A head protocol must therefore use the **linear/divergence-free** mode
(so `R_proj=0` exactly) *or* carry `R_proj` explicitly. The minimal run uses the linear class.

### 1.4 Gauge invariance — how the object is made invariant (and why the existing one is not)

The flat-chart `a=⟨v,X⟩` with a *fixed coordinate* plane wave `v` is **gauge-dependent**: under
`U→ΩUΩ†`, `X` transforms non-covariantly and `a` changes. The existing measurements therefore depend on
the (implicit, unfixed) gauge. Two admissible gauge-invariant realizations of (1.1) exist:

- **(A) Orbit-space / gauge-fixed realization (primary, practical).** Fix Landau gauge; take the
  transverse projected field `A^⊥`; the block-scale amplitude `a=⟨v^⊥,A^⊥⟩` with `v^⊥` at the *block*
  momentum `k_b=2π/ℓ_b` (not `k→0`). The measure is the gauge-fixed one, so `Φ` includes the
  gauge-fixing Jacobian **and the Faddeev–Popov determinant**; the orbit ambiguity (Gribov copies) is a
  *systematic*, controlled by averaging/reporting the spread over gauge-fixing starts (§3).
- **(B) Gauge-invariant order-parameter effective potential (cross-check).** Use a per-block
  gauge-invariant coarse variable (blocked-plaquette / Polyakov-loop phase in a representation), pool
  the values over blocks×configs, and take the curvature of `−log P_block(a)` at its mode. O(1)
  per-block variance ⇒ **non-vacuous**.

**Why (B) is not the vacuity trap.** H1-NUMERIC-RIGHT proved the *global* `k=0` Fourier mode `q̄` of the
composite `tr F²` is a susceptibility (`Var(q̄)∝V^{-1}`; `G″=1/Var∝V→∞`, H1 vacuous). That is a statement
about a **volume-averaged** mode; a **per-block** coarse spin at scale `ℓ_b` has O(1) width and O(1)
curvature. The head object is the latter. The existing runs measure a third thing — the **global
`k=2π/L` gluon amplitude in the unfixed chart** — which is neither (A) nor (B): it is gauge-dependent and
not block-local.

### 1.5 The estimator (what is actually computed)

For each config and mode `v`, the gold-standard `G″` is the curvature of `−log P(a)` at `a₀*` (Gaussian
KDE, quadratic fit over `|a−a₀*|<0.8 sd`). The decomposed estimator uses **Nadaraya–Watson kernel
conditional moments at `a₀*`** (`bw=0.30·sd`) for `E[H|a₀*], Var(G|a₀*), E[Φ″|a₀*], Var(Φ′|a₀*),
Cov(G,Φ′|a₀*)`, then (1.2). Errors: **delete-block jackknife** (20 contiguous blocks, block ≫ `τ_int`).
The one-sided bounds:

```
LB_rig  = (Δ_dyn + X) − 1.645·SE(Δ_dyn+X)        (keep ONLY the proven W_geo≥0)
LB_stat = G″ − 1.645·SE(G″)                        (trusts the full measured G″)
```

The whole content of the head is that **`LB_rig<0` while `LB_stat>0`** — i.e. the sign lives in
`W_geo`'s *magnitude*, which the theorem bounds below by 0 only.

---

## 2 · Reanalysis of the existing SU(2)/SU(3) data

Witness `h5_head_reanalysis.py` (deterministic; reads local `.npz`; delete-block jackknife `nblk=20`).
Datasets: `partA_su2_N4000` (SU2 `β=2.5, L=4`, 48 modes), `partA_su2_N800` (SU2 `β=2.5, L=6`, 48 modes)
— **exact head coupling**, store only `a,H,G`; `k0_su2_b2.4_L8_n2` (SU2 `β=2.4, L=8`, 12 modes) and
`k0_su3_b5.7_L8_n2` (SU3 `β=5.7, L=8`, 12 modes) — **full decomposition** (`a,H,G,P1,P2`).

### 2.1 The result table (evaluate at the coarse minimiser `a₀*`)

| dataset | β | L (modes) | `a₀*` | `Δ_dyn` | `W_geo` | `X` | `G″` | `1/Var` | `LB_rig` | `LB_stat` |
|---|---|---|---|---|---|---|---|---|---|---|
| partA SU2 | **2.5** | 4 (48) | +0.02 | +0.001±0.014 | — not stored — | — | **+0.268±0.007** | 0.269 | **n/a** | +0.257 |
| partA SU2 | **2.5** | 6 (48) | −0.09 | −0.003±0.040 | — not stored — | — | **+0.327±0.016** | 0.272 | **n/a** | +0.301 |
| partB SU3 | **2.5** | soft-block (W7-Gpp0) | ≈0 | E[H]=0.22,Var=0.24 ⇒ br≈−0.02 (flat −0.03) | — | — | **+0.51±0.08** | — | **n/a** | +0.38 |
| k0 SU2 | 2.4 | 8 (12) | −0.12 | +0.089±0.071 | **+0.254±0.005** | −0.001±0.036 | **+0.341±0.077** | 0.276 | **−0.042** | +0.195 |
| k0 SU3 | 5.7 | 8 (12) | −0.035 | +0.040±0.088 | **+0.318±0.006** | +0.023±0.034 | **+0.382±0.098** | 0.397 | **−0.095** | +0.242 |

(For the `k0` rows, `G″` here is the identity `Δ_dyn+W_geo+X` evaluated at `a₀*`; the KDE gold standard
at `a₀*` gives +0.264±0.041 (SU2) / +0.398±0.095 (SU3), differing from the identity by the known
closure/systematic ±0.08–0.11. The `1/Var` column is the Gaussian proxy; `G″≈1/Var` to 2–5 %. The small
`a₀*≠0` values are KDE-mode noise on a broad (`sd≈1.9`) even distribution; symmetry fixes the true mode
at 0, so these are already at `(★_loc)`'s evaluation point.)

### 2.2 What the existing data CAN say

- **`G″(a₀*)>0` statistically, at the exact head coupling.** SU2 `β=2.5`: `+0.268±0.007` (`38σ` stat,
  `L=4`) and `+0.327±0.016` (`L=6`), 48/48 modes positive; SU3 soft-block `+0.51±0.08`. The sign is not
  a single-mode or single-config fluctuation. `[V]`
- **The carrier is `W_geo`, stably.** `W_geo=+0.254` (β=2.4) / `+0.318` (β=5.7) with SE ≈0.005;
  `Δ_dyn≈0` (equipartition `E[H]≈Var(G)`); `X≈0` (`|X|≲0.04`). `1/Var` matches `G″`. `[V]`
- **The evaluation point is the minimiser.** `a₀*≈0` (symmetric measure); no barrier/shoulder ambiguity
  on the actual lattice ensemble (unlike the one-plaquette `θ=π` toy, where the mode sits at `a*≈−2`). `[V]`

### 2.3 What the existing data CANNOT say

- **Not the gauge-invariant object.** Every stored amplitude is the unfixed flat-chart gluon
  `a=⟨v,log U⟩`; no gauge fixing, no Gribov/orbit control, no FP measure. The measured `G″` is a
  gauge-dependent number (K0-OVERLAP-AUDIT; H1-NUMERIC-RIGHT). It cannot certify the gauge-invariant `G″`.
- **No rigorous-style bound at the head coupling.** `partA` (β=2.5) stores only `a,H,G` — no `Φ′/Φ″` —
  so `W_geo` and `X` cannot be separated and `LB_rig=(Δ_dyn+X)−…` cannot be formed. The only two datasets
  with the split are at `β=2.4` and `β=5.7`, and there `LB_rig=−0.042/−0.095<0`: **0 not excluded.**
- **Not a continuum or `k→0` statement.** Single `β` per group; `L≤8`; softest `k=π/4` (L8), `π/3`
  (L6), `π/2` (L4); `a(β≈2.5)≈0.15–0.2 fm`, `ξ/a≲0.3` — outside the scaling window, no `a→0`, no
  thermodynamic limit, never `k→0`. The gluon `G″` at finite `k` is **below** the free value
  (`G″/free=0.18–0.63`), so it does not even pass the acid test at finite `k` (NUMERICS-AUDIT (a)).
- **Statistics are not the bottleneck.** `N=400–4000`, `τ≈1.5–2.6`, delete-block SE on `G″` `±0.007–0.1`.
  The dominant error is the **object/closure systematic** (`±0.1`), not the statistical SE.

---

## 3 · The minimal decisive computation `[C]`

Goal: a controlled, **gauge-invariant** one-sided bound on `G″(v;a₀*)` at `β=2.5`, with the full
`Δ_dyn+W_geo+X` decomposition and a Gribov-copy systematic, and a volume scan.

### 3.1 Protocol

- **Group / coupling.** SU(2) primary (cheapest; head `β=2.5` is SU(2)-native), SU(3) as a
  string-tension-matched cross-check (SU3 `β≈5.7`, or SU3 `β=2.5` deep-SC as the W7-Gpp0 point).
- **Volumes.** `L = 8, 12, 16` (three volumes ⇒ the `1/L²` finite-size law of ch12 §diagnostics; `L=16`
  is the first volume where the softest block mode is inside the scaling window). One job at a time
  (charter §4: **never parallel** SU(N) `L≥8` — the prior OOM).
- **Configurations.** `N = 2000–4000` after `500` thermalisation sweeps, stored every `gap=5–10` sweeps
  (`τ_int≈2–3` stored configs ⇒ delete-block bins ≫ `τ_int`). Streaming sampler (configs not stored),
  peak RSS ≲ 200 MB at `L=16` SU(2).
- **Gauge fixing.** Landau (FFT/over-relaxation), **`N_copy=5–10` random starts per config**; report
  the spread of `G″` across copies as the **gauge-invariance error bar**. Fold the gauge-fixing +
  Faddeev–Popov Jacobian into `Φ` (stochastic FP estimate, or a controlled bound on its effect);
  this is the one genuinely new numerical ingredient relative to the existing pipeline.
- **Mode.** The **block-scale softest blocked mode**: build the `b=2` block lift `P`, form
  `H̄=E[Pᵀ Hess V P]`, take the softest eigenvector `v*`, measure the exact marginal of
  `V_full=P v*` (`partB`'s construction) on the **gauge-fixed** ensemble — at the block scale
  `k_b=2π/ℓ_b`, **not** `k→0`. Use the linear/divergence-free projection (§1.3) so `R_proj=0` exactly.
- **Estimator.** (1.2) at `a₀*` with the kernel-conditioned moments (§1.5); store `H,G,Φ′,Φ″` per config
  per mode so both `LB_rig` and `LB_stat` are formable. Add Gribov-copy spread in quadrature with the
  delete-block jackknife SE and the KDE-bandwidth systematic.
- **Target.** Report the achieved `LB_rig` and `LB_stat`. Set the tested floor `δ_W=0.1` (a conservative
  fraction of the measured `≈0.25–0.5`) and ask whether `LB_rig≥δ_W`. **Honest expectation:** `LB_stat`
  will be `>0` (as at β=2.4), but `LB_rig` may remain `<0` because `W_geo`'s magnitude stays uncertified
  — the run *measures the correct object*, it does not turn MC into a proof.

### 3.2 Cost estimate (grounded in measured runtimes)

Anchor: existing `k0` SU(2) `L=8`, `N=400`, 12 modes = **364 s** (≈0.076 s per mode·config at `L=8`);
`h1_right` SU(2) `L=8`, `N=200`/β ≈ 100–134 s. Scaling by `L⁴`, one mode, `N` configs:

| volume | sites | measurement (1 mode, N=2000) | sampler | + gauge-fix ×5 copies | subtotal (py) |
|---|---|---|---|---|---|
| L=8 | 4096 | ≈0.08·2000 = 160 s | ≈0.2 h | ≈0.3 h | ≈0.6 h |
| L=12 | 20736 | ≈0.38·2000 = 770 s | ≈0.6 h | ≈0.9 h | ≈1.9 h |
| L=16 | 65536 | ≈1.2·2000 = 2400 s | ≈1.8 h | ≈2.8 h | ≈5.3 h |

**SU(2), three volumes, one block mode: ≈8–10 CPU-hours pure-Python/numpy; ≈1–2 h on C/OpenMP 16
cores; safe memory if run one job at a time.** SU(3) cross-check ≈2–3× (≈1–2 days C/OpenMP). This is
*smaller* than ch12's "SU(3) L=16–24 × three couplings over days," because SU(2) at one coupling, one
block mode, is the minimal case. The decisive cost is **systematics (gauge-fixing + Gribov control +
FP), not statistics** — `N=2000` already drives the statistical SE below the ±0.1 object/systematic
floor. An `a→0` scaling window (ch12 §decisive-test) remains the far, expensive target
(`ξ/a≈7` needs `β≈7.6`, `~10⁴–10⁵` CPU-h; STATUS) and is **not** needed to test the head *floor* at
`β=2.5` — only to claim a continuum gap.

### 3.3 What this run decides — and does not

- **Decides (statistically):** whether the *gauge-invariant* conditioned coarse curvature at `β=2.5`
  has `G″(a₀*)>0` with a controlled Gribov/systematic error, with the correct `Δ_dyn+W_geo+X`
  decomposition and the volume law. This upgrades the existing evidence from the *wrong object* to the
  *right object*. `[C]`
- **Does not decide (rigorously):** a *uniform* lower bound `δ_W>0` independent of `(L,U′,v)`. Keeping
  only the proven `W_geo≥0` leaves `LB_rig<0`; closing it is exactly `(◆)` = `Var(G|a₀*)−E[H|a₀*] <
  W_full` = H1 `[O]`. No lattice dataset can supply a rigorous finite-`(β,L)` enclosure (PROOF-VALIDATED-
  GPP-SIGN §3: the 8192-dim integral / certified transfer matrix is infeasible). A *protocol* is the
  correct deliverable for the head; a *proof* is the wall.

---

## 4 · Acid test and label ledger

- **Acid test (free field must fail).** Free massless `Φ=0`, Gaussian `V`: `W_geo=0`, `Δ_dyn=0`, `G″=0`
  (no floor manufactured). At finite `k` the free value is `4sin²(k/2)>0` and the interacting `G″` is
  *below* it (`G″/free=0.18–0.63`, NUMERICS-AUDIT (a)) — so the head `G″` is **not** a gap signal at
  finite `k`; it is a curvature floor (necessary-shaped). The head inequality `δ_W>0` must be read as a
  *coarse curvature* statement, not an atom statement. `[T]`+`[V]`

| item | statement | label |
|---|---|---|
| true object = conditioned coarse curvature `G″(v;a₀)=E[∂²V\|a₀]−Var(∂V\|a₀)=Δ_dyn+W_geo+X` | (1.1)/(1.2); includes `X`; `\|a₀` load-bearing | **[T]** (identity, from corpus) |
| coarse field = minimiser `a₀*=argmax ρ_v`; symmetric ⇒ `a₀*=0` | V3 §4; witness `a₀*∈[−0.12,+0.02]` | **[V]** |
| projection = linear/divergence-free; covariant adds `R_proj=O(W_geo)` at head | G1/G3 | **[T]**/**[NG]** (inherited), **[C]** protocol |
| `G″(a₀*)>0` at `β=2.5`, gauge-dependent gluon (existing) | +0.268/+0.327/+0.51 | **[V]** |
| `W_geo>0` measured, `Δ_dyn≈0`, `X≈0` at β=2.4/5.7 | +0.254/+0.318, SE 0.005 | **[V]** |
| rigorous-style one-sided bound (only proven `W_geo≥0`) at β=2.4/5.7 | **−0.042 / −0.095 < 0** | **[V]** |
| rigorous-style bound at the head coupling β=2.5 | **not formable** (no `Φ′/Φ″` stored) | **[O]** / **[V]** (absence) |
| existing object = gauge-invariant conditioned coarse curvature | **no** (unfixed flat-chart global gluon) | **[C]** (cannot say) |
| `Var(q̄)∝1/V` vacuity of the *global k=0 composite* mode | per-block coarse spin is O(1), non-vacuous | **[T]**/**[NG]** (inherited) / **[C]** |
| uniform `δ_W>0` at the head, `(★_RG)` | = `E[Σ\|a₀]≥δ_W` = H1 | **[O]** |
| minimal decisive run: gauge-fixed + Gribov-averaged block mode, β=2.5, L=8/12/16, N≈2000 | §3; ≈1–2 h (C/OpenMP) SU(2) | **[C]** |
| a lattice test proves the head | **no** — statistical only; rigorous is H1 | **[C]**/**[O]** |
| Clay mass gap | open | **0%** |

---

## 5 · Sharpest line

**The existing corpus measures the head coupling's `G″` on the wrong object — the unfixed flat-chart
global gluon — and gets a real, statistically positive `+0.27…+0.51` whose carrier `W_geo` the theorem
bounds below by `0` only; keeping just that proven half, the one-sided lower bound is `−0.04/−0.10`, so
nothing is certified, and at `β=2.5` itself the `Φ`-split is not even stored.** The correct object is the
**gauge-invariant conditioned coarse curvature at the coarse minimiser `a₀*≈0`** (§1), with the linear
divergence-free projection and `X` included. The minimal decisive run is a gauge-fixed, Gribov-copy-
averaged, block-scale softest-mode measurement of `Δ_dyn+W_geo+X` at `β=2.5` on `L=8,12,16`,
`N≈2000–4000` (≈1–2 CPU-hours SU(2) in C/OpenMP, one job at a time) — which can convert evidence on the
wrong object into a controlled **statistical** bound on the right one, and nothing more: the uniform
floor `δ_W>0` is H1, not a dataset. **No QED. Clay 0%.**

---

*Refs.* `STAR2-CLASSIFICATION.md` §0–§7 (head = H1); `V2-TRUE-BOUND.md` Lemma V2.1 / §1; 
`V3-COUNTEREXAMPLE-SCOPE.md` §0,§3–4 (`a₀*`, `(★_loc)`); `RECONCILIATION-STAR.md` §0–4 (`G″=Δ_dyn+W_geo+X`,
`\|a₀` load-bearing); `G1-BLOCK-IDENTITY.md` §3–4 (`G″=2W_geo+Δ`); `G3-PROJECTION-TERM.md` §4–6 (`R_proj`);
`S1-CONDITIONAL-POINCARE.md` §0 (`S=1/Var`, H1); `H1-NUMERIC-RIGHT.md` §0–5 (composite vacuity; gauge-fixed
route); `K0-OVERLAP-AUDIT.md` (gluon vs `0⁺⁺`); `COMPUTE-GPP-K0-CONTINUUM.md` §1–5 (`G″=bracket+W`,
`LB_rig`); `COMPUTE-GPP-SU3.md`; `NUMERICS-AUDIT.md` (acid test, finite-k); `OTHER-GPP-KZERO.md`
(`(◆)=H1`); `W7-Gpp0.md` (β=2.5 exact marginal +0.27/+0.51); `PROOF-VALIDATED-GPP-SIGN.md` §3 (no
enclosure from MC); `CONTINUUM-VARIANCE-DEFECT.md`; `docs/paper-v2/ch12-numerics.tex` §decisive-test.
Witness `research/v3/w7/h5_head_reanalysis.py` (+ `h5_head_reanalysis_results.json`). No third-party
sources downloaded this session.
