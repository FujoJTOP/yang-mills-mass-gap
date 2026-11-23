# STATUS — lattice test of CPB-W

> 2026-11-22 · 4D pure SU(3) Wilson lattice · full report: `REPORT.md`

## Verdict

**Scientific: INCONCLUSIVE.  By the project's literal rule: FAIL (`κ_eff<0`).**

| criterion | primary grid β=2.6–3.3 (L=6) | literal |
|---|---|---|
| `R₄ ≫ 1` ⇒ fail | `R₄(b=1)=−0.125…−0.131`, `R₄(b=2)=−0.012…+0.003` | **passes** |
| `κ_eff > 0` ⇒ pass | `κ_eff = −2.82…−3.42` (100 % of configs) | **fails** |

The `κ_eff` "fail" **must not be read as a failure of `CPB-W`**: the same
`κ_eff` is **−1.74 at β=40**, where the theory is known to be healthy.  It is a
Savvidy/Nielsen–Olesen negative mode of the *microscopic* action
(`|κ_eff| ≈ (0.85–1.6)·β(1−⟨P⟩)`), and `handoff-proof` (T4) has excluded it
from the CPB decision.  It is reported, not used.

**Decision-relevant bottom line:** the test supplies **no positive evidence**
for `CPB-W`, and does not falsify it either.  `ξ/a` is unresolvable at
β ≤ 3.3 (`≲0.3`), and the σ-certificate target `ξ/a≈7` needs `L≈28–32`
(~10⁵ CPU-h in pure NumPy) — a feasibility gap.  The decisive `z(q)`/`κ_c` is
not measurable from a direct MC; the closest proxies (all character
coefficients positive, `r₁=1.0000`) are consistent with the analytic default
`κ_c≈1`.

## What was computed

* CM heat bath for SU(3) (exact SU(2) subgroup sampler), validated against
  Metropolis, strong-coupling series, and known β=6 plaquette.
* `⟨P⟩`, `W(R,T)`, `R,T ≤ 3`, `ξ/a` (APE-smeared plaquette correlator),
  `R₄(b)` for `b=1..4`, `ρ(θ)`, `κ_eff`, on `L=4,6,8` (and `L=10` at β=6).
  **ξ**: β=6 → 0.93(8); β=1.7/2.2 → **unresolvable** (ξ/a ≲ 0.3, the
  correlation dies inside one lattice spacing).
* **Primary grid (physical handoff): β = 2.6, 2.8, 3.0, 3.3** (`g²=1.8–2.3`,
  `λ*=5.4–6.9`).  **Conservative strong-coupling checks: β = 1.7, 2.2**
  (`λ*=8.2–10.6`).  **Formal weak boundary only: β = 5.7, 6.0** (`g²≈1`,
  `λ*≈3`).  The earlier "β≈6 is the handoff" reading is superseded.
* `t_up` at `g²=2`: 1-loop 3.59, 2-loop 2.639, 3-loop 2.36
  (`q₀ = 0.221 / 0.571 / 0.758` for `m/Λ=8`).  **CPB-W needs a ≥2-loop bound**;
  a one-loop control only reaches `q₁ ≤ 0.22`.  (`t_up.csv`.)

## Bianchi / z(q) / κ_c (REPORT §9)

* **κ_c has no direct numerical route** (structural conclusion, not a
  preliminary): the polymer weights `w_q(C)` belong to the blocked effective
  action; a direct MC gives correlation functions only.  Needs IMC/RG matching
  or an analytic character expansion of a link-blocked Boltzmann factor.
* **Character coefficients: NO `c_R < 0`** at any β=2.2–3.3, L=6 or 8.
  All of c₃,c₆,c₈,c₁₀ are positive by ≫10σ (e.g. β=2.6: +0.5211, +0.1349,
  +0.2381, +0.0232).  Single-plaquette signed/absolute ratio r₁=1.0000.
  → no c₁-level sign structure; **consistent with BL-6** (κ_c=1), not a
  measurement of κ_c.
* **σ(1) per-block proxy FAILS**: for fixed (β,L,b) the thresholds give
  c₀=−ln f = 2.4–3.8 / 5.4–10.1 / >10.7 (P<0.75⟨P⟩ / P<⟨P⟩/2 / P<0).
  The measured fractions equal the **Gaussian tail** computed from the same
  mean and variance (β=2.6,b=2: predicted 0.00099/0.0610 vs measured
  0.00089/0.0602) → c₀ measures the width of the blocked field, not an action.
  "c₀≈5.5" is one cell of that table.  L=6 vs L=8 agree to <1%.
* Defect clusters: dilute at P<0, percolating at P<⟨P⟩/2; no closed-surface
  vs domain-wall discrimination possible here.
* **Pilot (L=24, β=6.30, 3000 sweeps, Numba parallel)**: `ξ/a = 1.13(20)`,
  `τ_P=1.59`, `N_eff≈95`; fit stable over windows t=1..3…10 (1.01–1.13);
  smeared (n=8) far better than n=0.
* **β extrapolation corrected — factor-6 error in the exponent**: with
  `g²=6/β`, `ξ/a ∝ e^{β/(12b₀)}`, slope 1.196 not 7.18.  Two-loop:
  `ξ/a(6.3)=1.36` (pilot 1.13 ✓), **ξ/a=7 at β≈7.6, 8.5 at β≈7.8**.
  A 32⁴ run at β=6.30 would measure ξ/a≈1.3 — **do not run it at 6.30**.
* **Aggregate throughput measured**: the machine saturates at **~50×
  numpy-core total** (K×T splits 1×32…16×1 all give 28–54×).  Per-chain
  speedups cannot be multiplied by the chain count — the earlier
  "16 chains × 50×" projection was 16× too optimistic.
* **ξ/a≈7 (RB2/σ cap 0.143)**: Numba kernel **7.0× serial, 55–78× parallel**
  (vs single-core NumPy: L=6 55×, L=8 74×, L=12 78×; 32 threads, validated:
  plaquette/action 2e-17, staples 1e-15, ⟨P⟩ agrees, smearing 3e-11).  Projection (REPORT §10):
  N_eff ≈ 6×10⁵ independent configs for 5% on ξ at ξ/a=7, i.e. ~1.2×10⁷
  sweeps.  **16 independent chains / 16 cores**: 16³×64 = 400 h @7×,
  56 h @50×; 24³×64 = 1 330 / 187 h; 32⁴ = 1 590 / 222 h.
  At **20 %** precision (N_eff/16): 3.5–14 h @50××16 chains, but 20 % only
  certifies `q₁^σ≤0.143` if the measured `ξ/a ≳ 11.7` (5 % needs ≳7.8).
  Measuring `ξ/a=7` needs `L≈21–28`; a scan on L=8/12 reaches only `ξ/a≈3–4`.
  For the stricter target (5 % + 2σ ⇒ `ξ/a≥8.5`, `L_s/ξ≥4`): **32³×64 =
  444 h** @50××16 chains (148 h with GEVP×3), but 16 chains need 19.3 GB
  (machine has 33.5 GB).  24³×64 (187 h) has a 2–5 % finite-size systematic
  at `L_s/ξ≈3` — borderline for a 5 % target; cost ratio only 2.4×.

## Files

`su3.py` `validate.py` `run.py` `corr_run.py` `corr_analyze.py` `kappa.py`
`kappa_run.py` `kappa_control.py` `kappa_summary.py` `analyze.py`
`runs/` `corr/` `results.csv` `kappa_summary.csv` `summary.txt` `plots/`

## Next (for the main session)

1. Decide whether the `κ_eff` criterion should be restated for an
   **effective** (blocked) action; in its present form it is not a test.
2. Define and measure `z(q)` for the equivariant block spin, if `CPB-S`
   is to be tested numerically.
3. If a `CPB-W` proof is attempted, note that this data supports
   *smooth, Gaussian-like* blocked plaquettes (`R₄` small) but says nothing
   about the massive-background-field class itself.

## NA-3-T (Dobrushin coefficient) — 2026-11-23

`dobi.py` `dobi_run.py` `dobi_analyze.py`; results in `results_eta.csv`,
`na3t_table.txt`, `REPORT-NA3T.md`, witnesses in `runs/witness_L4_b*.npz`.

Every grid point is labelled **cert** or **est**.  A `cert` row means a
recorded witness pair (config U, link e, 18 neighbour replacements g) gives
`Σ_e' TV_lb ≥ 1`, hence `η ≥ 1` — a strict kill of uniform contraction.
`TV_lb` is the witness-function bound `½|E_p f − E_q f|` (f = Re tr(U ΔS)/‖ΔS‖_*)
minus 4σ MC error.

| β | η̂ | witness lb | cert |
|---|---|---|---|
| 2.2 | 5.74(0.23) | 1.18 | cert |
| 3.0 | 8.38(0.36) | **1.75** | cert |
| 6.0 | 17.19(0.10) | 1.82 | cert |

β=3.0 (the decisive point) is a **certified** `η ≥ 1`.  Caveats: the bound is
statistical (4σ), not interval arithmetic; it certifies kill only — `η̂ < 1`
would never have certified survival.  `m_eff` at L=4 is unresolvable, so D4 is
not assessable.  The `b=2` block kernel is not implemented.
