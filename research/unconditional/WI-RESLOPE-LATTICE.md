# WI-RESLOPE-LATTICE — lattice measurement protocol for the Re-slope signature (WI line 2, real-τ); **more robust to a continuum background than WI-SU3's raw nearest-zero/splatter, at the cost of resolving the zero ladder**

> **VERDICT (no QED; Clay 0%).** A concrete, runnable lattice protocol for the **Re-slope** atom signature
> (§3 of `WI-MODIFIED-FUNCTIONAL`): reweight the `0⁺⁺` operator (smearing family), recover the **zero ladder**
> of the full correlator's analytic continuation, find **persistent `|Im τ|` rungs**, and test whether `Re τ`
> along them is **log-linear in the atom weight `c₁` with slope `1/Ê₁`** (and the rung sits at
> `π/|Im τ| = Ê₁`). **Witness (`w7/wi_reslope_lattice.py`, exact Newton zeros): on an ATOM riding a gapped
> OR at-threshold continuum, the WI-SU3 raw nearest-zero `spread` is spoiled (`≈1.8` / `≈0.8`) — it FAILS —
> while the Re-slope gives `R=5` atom-slope rungs; a pure CUT gives `R=0`.** So the Re-slope **is** the
> continuum-robust member of the WI pair. **But** it reads the **pole log-slope = A3 (read-not-force)**,
> is **data-hungry** (needs the ladder, not one zero), and its uniform-`(L,a)` statement at `β≈2.5` is `H1`.

---

## §0 · Labels

**[T]** exact identity · **[C]** protocol / conditional diagnostic · **[O]** exact missing step (uniform `(L,a)` = `H1`) · **[V-num]** synthetic witness · **[NG]** no-go.

---

## §1 · The measured object (same as WI-SU3)

For a `0⁺⁺` interpolating operator `O_r` (zero spatial momentum, parity-even), the **full** correlator
`Z_r(τ)=w₀^r+Σ_{n≥1} c_n^r e^{−E_n τ}`, `w₀^r=|⟨O_r⟩|²>0`. **Reweighting** = same `β`, `L`, configs, vary
only the operator (smearing level `r`): it changes the overlaps `{c_n^r, w₀^r}`, **not** the levels `{E_n}`.

**Operator family `r`** (all `0⁺⁺`, zero momentum, same `β,L`): `0` = 1×1 plaquette `tr F²`; `1` = 2×2 clover;
`2` = APE/HYP 1-level; `3` = APE/HYP 3-level; `4` = Gaussian-smeared (different radius). (Same as `WI-SU3 §2`;
fix the smearing **levels** a priori.)

---

## §2 · Real-τ extraction: GEVP levels + the zero ladder

Full complex-`τ` continuation is ill-posed, so extract on **real** τ:

1. **Levels via GEVP `[C]`.** On the correlation matrix of the operator family, solve the generalized
   eigenvalue problem at a reference time-slice (with the usual `t` and `t+1` matrices). This returns
   `Ê_n` **operator-independent** (the whole point of reweighting) and, per `r`, the overlaps `ĉ_n^r, ŵ₀^r`
   (linear solves). Require the GEVP plateau first (`N_t=32–48`, `L=16–24`).
2. **Zero ladder.** With `(ĉ_n^r, Ê_n, ŵ₀^r)` fixed, find **all** zeros of `Z_r(ζ)=ŵ₀^r+Σ_n ĉ_n^r e^{−Ê_n ζ}`
   in `ζ∈ℂ` (Newton over a grid; seeds = the atom ladder `ζ₀=(ln(ĉ₁/ŵ₀)−i(2k+1)π)/Ê₁`, `k=1..K`, plus a
   coarse grid to catch spurious continuum zeros). Collect their `|Im ζ|` **rungs**.
   *(Equivalently, ESPRIT/Prony on the real-τ Hankel gives the same `(c,E)` modes — `WI-REALTIME`; then the
   zeros are Newton on the recovered sum. Use ESPRIT if the GEVP family is small.)*

**Why deeper rungs, not the nearest zero:** a continuum creates **spurious zeros at small `|Im|`** (the
`WI-ADMISSIBLE`/`WI-EMBEDDED` failure that kills the raw nearest-zero `𝒜`). The atom's **deep** rungs
(`|Im|=(2n+1)π/Ê₁`, `n≥2`) sit far from those, so they are the clean signal — hence a **ladder**, not one zero.

---

## §3 · The Re-slope decision statistic

For each rung index `n` (target `|Im|=(2n+1)π/Ê₁`):

1. For each reweighting `r`, match the zero whose `|Im|` is nearest the target (within `tol·π/Ê₁`); record
   its `Re τ`.
2. Keep rungs **persistent** across `r` (matched for every `r`).
3. Along each persistent rung, fit `Re τ` vs `ln ĉ₁^r`; slope `s_n`.

**Decision statistic:** `R = #{ n : |s_n − 1/Ê₁| < tol_slope  AND  |π/|Im τ_n| − Ê₁| < tol_rung }`.

- **Isolated atom** (simple pole of `Π` at `s=−Ê₁`): every rung has `|Im τ_n|=(2n+1)π/Ê₁` (invariant under
  `c₁`) and `Re τ_n=ln(c₁/w₀)/Ê₁` (slope `1/Ê₁`) ⇒ `R ≥ 3–5`.
- **Cut / two-particle threshold** (no pole): zeros move in both `Re` and `Im` under reweighting; no
  `c₁`-stable rung with slope `1/Ê₁` ⇒ `R = 0`.

**Witness `[V-num]`** (`w7/wi_reslope_lattice.py`, `w₀=1`, `Ê₁=1`, reweight `c₁∈{0.4,1,2.5}`, 40-level
continuum, 1% weight):

| model | RAW (WI-SU3) `𝒜_r` per `c₁` → **spread** | Re-slope **`R`** |
|---|---|---|
| **ATOM** (isolated) | `1,1,1` → **0.00** | **5** (slopes `1.000`) |
| **ATOM+GAP** (edge `2>E₁`) | `2.82, 2.69, 1.00` → **1.82** | **5** (slopes `1.003`) |
| **ATOM+THR** (edge `=E₁`) | `1.79, 1.77, 1.00` → **0.79** | **5** (slopes `0.998`) |
| **CUT** (no atom) | `2.90, 2.72, 2.51` → 0.39 | **0** |

⇒ the recurrent failure of the raw `𝒜` spread on a physical atom+continuum (`ATOM+GAP` spread `1.82`) does
**not** reach the Re-slope: `R=5` survives the continuum, `R=0` for the pure cut. **This is the protocol's claim.**

---

## §4 · Parameters (β / L / statistics)

- **Gauge group / action:** SU(3) Wilson, `β=6/g²`. (The project's SU(2) `β≈2.5` with `β=4/g²`, `g²≈1.6`,
   maps to SU(3) `β≈3.75`.)
- **Primary run:** `β ∈ {5.7, 5.85, 6.0}` (glueball window), `L³×N_t = 16³–24³ × 32–48` — the project's
   `L=4–8` only for feasibility; at small `L` the cut blurs into a dense ladder and the deep rungs are not
   resolvable.
- **Crossover link:** one stronger point scanning toward `β≈3.75` (the `α_s~O(1)` region).
- **Statistics:** `O(10³–10⁴)` configs per `β`. The Re-slope is **more data-hungry** than WI-SU3 (it must
   resolve `≥3–4` deep rungs, not one zero — Hankel conditioning is the cost, `WI-REALTIME`). Budget for it;
   a bootstrap over the GEVP fits gives the `R`-statistic error.

---

## §5 · Expected SU(3) signal

- **Physical lightest `0⁺⁺` glueball** (isolated atom, `m/√σ≈3.9`): the rungs sit at `(2n+1)π/m₁` with `Re`-slope
  `1/m₁` ⇒ **`R ≥ 3`**, `R` stable across `β` and smearing.
- **A two-glueball-cut operator** (deliberately two-particle-like): no persistent atom-slope rung ⇒ **`R=0`**
  (negative control).
- **Cross-checks:** a `2⁺⁺` operator (a clean atom in another channel) should also show `R ≥ 3`; a mixed
  operator with large continuum overlap should show degraded `R` — the contrast *is* the test.

---

## §6 · Comparison to WI-SU3 (is the Re-slope more robust?)

| | WI-SU3 (raw nearest-zero) | **WI-RESLOPE (this)** |
|---|---|---|
| statistic | spread of `𝒜=π/|Im α_{nearest}|` | `R` = #atom-slope persistent rungs |
| **continuum robustness** | **FAILS** for a physical atom+continuum (a 1% continuum moves the nearest zero; witness spread `1.82`) | **ROBUST** — deep rungs avoid the spurious nearest zeros (witness `R=5` on `ATOM+GAP`) |
| data cost | one nearest zero per `r` | the **ladder** (≥3–4 deep rungs) per `r` — **more data-hungry** |
| real-τ extraction | fit `M=2–3` exponentials → nearest zero | GEVP levels + all zeros (or ESPRIT modes → Newton on all zeros) |
| content | single-exponential-dominance (`WI-ADMISSIBLE`) | **A3** (pole log-slope `f₁²=−Res Π`) |

**Answer:** **yes, more robust to a continuum background** (the decisive practical failure of WI-SU3 is fixed:
the deep rungs are unaffected by the continuum-created small-`|Im|` zeros), **but at a higher data cost** (the
ladder must be resolved). Both are **A3 read-not-force**; neither is a gap certificate.

---

## §7 · Honest limits `[O]/[NG]`

- **A3 (read-not-force).** `(rung invariant, Re-slope=1/Ê₁)` *is* "a simple pole of `Π` at `s=−Ê₁`" = the
  Stieltjes residue `f₁²=−Res Π`. It changes the **representation**, not the content = the `A3` coordinate of `H1`.
- **Data-hungry.** Resolving `≥3–4` deep rungs from noisy real-τ data is ill-conditioned (Hankel
  condition number grows with mode count, `WI-REALTIME`); whether realistic SU(3) data reach `R≥3` decisively
  is **open** — a pilot at `β=6.0`, `L=16`, `N_t=48` is the first check.
- **Toy calibration only.** Synthetic exact exponential sums; a delta truly merged into the cut (`f₁²→0` at the
  edge) would not register (as in `WI-MODIFIED-FUNCTIONAL §6`).
- **Uniformity is `H1`.** Proving the Re-slope signature **uniform in `(L,a)` at `β≈2.5`** (the pole present,
  the deep rungs readable, through the continuum limit) is exactly `H1`. **Clay 0%, no QED.**

*Files:* `research/unconditional/WI-RESLOPE-LATTICE.md`; witness `research/unconditional/w7/wi_reslope_lattice.py`
(`py -3.12 -I`, numpy only, prints ATOM/ATOM+GAP/ATOM+THR/CUT: RAW-spread vs Re-slope `R`).
*Prior:* `WI-MODIFIED-FUNCTIONAL.md` (§3 the Re-slope), `WI-SU3.md` (the raw protocol compared here),
`WI-REALTIME.md` (ESPRIT real-τ estimator), `WI-ADMISSIBLE.md`/`WI-EMBEDDED.md` (why the raw nearest zero fails).
