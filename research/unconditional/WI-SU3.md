# WI-SU3 — a concrete SU(3) lattice protocol for the weight-invariance (WI) atom test

> **What this is.** A real, runnable measurement protocol to test the finite-`L` **atom signature** of
> `CONSTRUCT-IR-ANALYTIC-NUM §4`: for the `0⁺⁺` full correlator `Z_O(τ)=w₀+Σ_i c_i e^{−E_i τ}`, the IR-analytic
> functional `𝒜=π/Im(τ_zero)` is **weight-invariant for an isolated atom** and **weight-dependent for a cut**.
> "Weight-invariance" = `𝒜` does not move when the `0⁺⁺` interpolating operator is **reweighted** (different
> smearing), which changes the overlaps `c_i, w₀` but **not** the energies `E_i`. This is a **diagnostic**
> `[C]`, not a certificate: turning it into a uniform-`(L,a)` statement at the crossover is `H1` `[O]`.
> Synthetic witness verified; **no SU(3) data taken here.** Clay 0%; no QED.

---

## §0 · Labels

**[T]** exact identity (atom ⇒ `𝒜=m₁`, weight-independent) · **[C]** protocol / conditional diagnostic ·
**[O]** exact missing step (uniform `(L,a)` proof = `H1`) · **[V-num]** synthetic witness demonstration · **[NG]** no-go (not claimed here).

---

## §1 · The object measured

For a `0⁺⁺` interpolating operator `O_r` (zero spatial momentum, parity-even, time-reflection-even), measure the
**full** two-point function (not just connected):

```
Z_r(τ) = ⟨O_r(τ) O_r(0)⟩_full = w₀^r + Σ_{n≥1} c_n^r e^{−E_n τ},     w₀^r = |⟨O_r⟩|²,  c_n^r = |⟨0|O_r|n⟩|² ≥ 0.
```

The VEV² term `w₀ = |⟨tr F²⟩|² > 0` is **essential** (it is the constant that lets `Z` cross zero in complex `τ`;
`⟨tr F²⟩ ≠ 0` for the gluon condensate). Energies `{E_n}` are operator-independent; **reweighting** `r↦r'`
changes only `{c_n^r, w₀^r}`.

The functional `𝒜_r := π / Im(τ_zero^r)`, `τ_zero^r` the zero of the analytic continuation `Z_r(ζ)`, `ζ∈ℂ`, nearest
to the real axis. **Exact fact `[T]`:** if the lightest level `E₁` is isolated, `𝒜_r = E₁` for **every** `r`
(the zero sits at `Im τ = π/E₁`, the overlap `c₁/w₀` only shifts `Re τ`).

---

## §2 · The operator family (the "reweighting")

Use the **same** quantum numbers `0⁺⁺`, the **same** `β`, `L`, and gauge configurations — vary only the
interpolating operator. This keeps `{E_n}` fixed and reweights `{c_n, w₀}`:

| `r` | operator `O_r` (all `0⁺⁺`, zero momentum) |
|---|---|
| 0 | elementary plaquette `tr F²` (1×1) |
| 1 | `2×2` clover / blocked `tr F²` |
| 2 | APE/HYP-smeared `tr F²` (1 level) |
| 3 | APE/HYP-smeared `tr F²` (3 levels) |
| 4 | `tr F²` on a Gaussian-smeared link set (different radius) |

(Smearing changes the overlap with each state — precisely the reweighting. Do **not** change `β` or `L` across
`r`, or `{E_n}` changes and the test is void.)

---

## §3 · Parameters

- **Gauge group / action:** SU(3) Wilson, `β = 6/g²` (standard convention; note the project's SU(2) `β≈2.5` uses
  `β=4/g²`, i.e. `g²≈1.6`; the SU(3) point with the same `g²` is `β≈3.75`).
- **Primary run:** `β ∈ {5.7, 5.85, 6.0}` — the standard SU(3) glueball window where the `0⁺⁺` (and the
  `0⁺⁺*` at `≈1.6 m₁`) are resolvable; `L³×N_t` with `N_t = 32–48` for a clean exponential fit and `L = 16–24`
  (the project's `L=4–8` is only a feasibility check — at small `L` the "cut" becomes a dense level set and the
  atom/cut contrast blurs).
- **Crossover link:** one stronger-coupling point (scan toward `β≈3.75`) to connect to the project's `α_s~O(1)`.
- **Statistics:** `O(10³–10⁴)` configurations per `β` (the disconnected `w₀` is the noisy part; but see §5, `𝒜`
  is largely insensitive to `w₀` errors because the mass rides on `Im τ`).

---

## §4 · Extraction of `𝒜_r` (real-`τ` proxy, no full continuation)

Full complex-`τ` continuation from discrete real-`τ` data is ill-posed. The **practical proxy**:

1. Estimate `w₀^r = ⟨O_r⟩²` from the VEV.
2. Fit `Z_r(τ) − w₀^r = Σ_{i=1}^{M} c_i^r e^{−E_i^r τ}` with `M=2–3` (constrained multi-exponential fit, or a
   GEVP to fix `E_i`; use the same `E_i` across `r` if the fit is stable — a consistency check).
3. Solve the fitted sum `w₀^r + Σ_i ĉ_i^r e^{−Ê_i ζ}=0` in `ζ∈ℂ` (Newton from the atom guess
   `ζ₀=(ln(ĉ₁/w₀)+iπ)/Ê₁`); take the nearest zero, `𝒜_r = π/Im ζ_zero`.

**WI statistic:** `spread_r(𝒜_r) = max_r 𝒜_r − min_r 𝒜_r`. **Prediction:** atom ⇒ `spread ≈ 0`; cut / threshold
⇒ `spread > 0`.

**Synthetic validation `[V-num]`** (`w7/wi_su3.py`, numpy, machine-precision zeros):

| case | reweighting | `𝒜` | spread |
|---|---|---|---|
| **atom** `m₁=1` | overlaps `0.3 / 1 / 5` | `1.000000` (×3) | **`0.0`** |
| **atom + 2-glueball cut** | continuum overlap `0 / 0.3 / 0.8 / 2.0` | `1.000 / 0.991 / 0.970 / 0.912` | **`0.088`** |
| **cut (two comparable levels 1.00,1.35)** | weights `(1,0)/(1,.3)/(.3,1)/(1,1)` | `1.000 / 1.058 / 1.257 / 1.139` | **`0.257`** |

⇒ the atom is `𝒜`-invariant; any cut contaminates and **moves** `𝒜`. This is the quantity to measure.

---

## §5 · Prediction (what SU(3) should show)

- **Physical lightest `0⁺⁺` glueball** (believed an isolated atom at `m₁≈1.7 GeV`, `m/√σ≈3.9`): `𝒜_r ≈ m₁` for all
  smearing levels, `spread_r(𝒜_r)` consistent with 0 within errors.
- **If the lightest `0⁺⁺` region is a continuum / two-glueball threshold** (`E≈2m₁`): `𝒜_r` **moves** with
  smearing, `spread` grows — especially for operators with large overlap to the two-glueball channel.
- **Cross-checks:** (i) a `2⁺⁺` operator (a different channel, generally a clean atom) should also be
  `𝒜`-invariant; (ii) a deliberately two-particle-like operator should show large spread — a negative control.

---

## §6 · Confounders (honest)

1. **Excited atoms, not a cut.** The `0⁺⁺*` at `≈1.6 m₁` is a second *atom*; if it is not separated from `m₁` in
   the fit, `𝒜` moves for the wrong reason. Use enough `N_t` and `M` to isolate the lightest.
2. **Finite-`V` blur.** At finite `L` the continuum is a dense discrete set spaced `~1/L`; the atom/cut contrast
   requires the atom to be isolated by more than the level spacing — hence `L=16–24`, not `L=4–8`.
3. **Continuation ill-posedness.** The `M=2–3` exponential proxy has systematic error; the atom-vs-cut contrast
   must exceed the proxy's fit uncertainty (bootstrap the fits).
4. **`w₀` noise.** The condensate `⟨tr F²⟩` is noisy and subtraction-sensitive; but `𝒜`'s mass rides on `Im τ`,
   which is `w₀`-insensitive to leading order (robustness, to be quantified).
5. **Reweighting must not change `{E_n}`.** Same `β`, `L`, configs; only the operator. Over-smearing that
   effectively projects onto the ground state is fine (it is the intended reweighting) — but re-optimizing the
   smearing per operator can hide a cut; fix the smearing *levels* a priori.

---

## §7 · Honest limits

- The WI test is a **diagnostic `[C]`**; it reads the atom **at finite `L` from real-`τ` correlators** — a genuine
  new handle that does **not** need `L→∞` or a full complex continuation for the *invariance* check.
- Proving that the invariance holds **uniformly in `(L,a)` at `β≈2.5`** (atom's zero stays weight-invariant
  through the crossover and the continuum) is **`H1` `[O]`** — the same wall, in the WI coordinate.
- **No QED; Clay 0%.** This protocol tests a reduction-side signature; it does not prove the gap.

*Witness:* `w7/wi_su3.py` (numpy; three cases: atom invariant, atom+cut moves, cut moves; run `py -3.12 -I`).
*Related:* `CONSTRUCT-IR-ANALYTIC-NUM.md` (§4 signature), `CONSTRUCT-WAVEFUNCTIONAL.md`, `PNT-IR-ANALYTIC.md`,
`INVENT-LEEYANG.md`.
