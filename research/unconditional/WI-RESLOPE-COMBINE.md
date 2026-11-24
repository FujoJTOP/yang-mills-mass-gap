# WI-RESLOPE-COMBINE — Re-slope (gives `E₁`) + approach-law (gives the gap `Δ`) = a cross-validated (mass, gap) DIAGNOSTIC; but `Δ` at finite `L` is the level spacing, and the isolation (`Δ` bounded as `L→∞`) is `H1`

> **VERDICT (no QED, Clay 0%).** The two robust WI signatures combine into a **stronger read** than either
> alone: the Re-slope reads the level `E₁` and the approach-law reads the gap `Δ` to the next state, and they
> **cross-check** (both must agree on `E₁`). This is a genuine gain over each separately — a cross-validated
> `(mass, gap)` estimator, robust to a background continuum. **But it is still two reads (`E₁` = pole location =
> A3; `Δ` = a spectral read), not a certifier:** at finite `L` the spectrum is discrete so the approach is
> *always* exponential with `Δ` = the level spacing, and the **isolation** (`Δ` stays bounded away from 0 as
> `L→∞`, i.e. the power-law never kicks in) is exactly `I2 = H1`.

---

## §0 · Setup

Connected correlator `C_F(τ)=Σᵢ cᵢ e^{−Eᵢτ}` (`w₀` the VEV term); "operator" `F` = the weight vector `{cᵢ}`
(levels `{Eᵢ}` fixed under reweighting). Witness: `E=[1.0,2.5,3.0,4.0]`, atom at `E₁=1`, `Δ=E₂−E₁=1.5`.

**LINE 2 — Re-slope.** `𝒜`'s complex-τ zeros sit on `Re τ = ln(c₁/w₀)/E₁`; along the rungs
`|Im τ|=(2n+1)π/E₁` the real part slides at `dRe τ/d ln c₁ = 1/E₁` ⇒ reads `E₁`.

**LINE 1 — approach-law (P3).** `m_F(τ)=−∂_τ log C_F → e₁`; the *approach law* separates isolated-atom
(exponential `e₁+C e^{−Δτ}`, `Δ=e₂−e₁`) from a continuum edge (power-law `e₁+1/τ`).

---

## §1 · Findings (witness `w7/wi_reslope_combine.py`, exact/quadrature)

- **Re-slope recovers `E₁` exactly and robustly**: pure-atom `c₁`-sweep `Re τ = {−0.693, 0, 0.693, 1.386}`,
  fitted slope `dRe/dln c₁ = 1.000000` (= `1/E₁`); the value survives atom+continuum operators (`O_a,O_b,O_c`).
- **Approach-law separates the two regimes**:
  - atom + discrete ladder (`Δ=1.5`): `(m−E₁)e^{Δτ} → 0.7503 = (c₂/c₁)Δ` (**constant ⇒ exponential**) — the rate is the **finite-`L` spacing**.
  - true continuum edge `[E₁,E₁+W]`: `(m−E₁)τ → 0.9999` (**constant ⇒ power law `1/τ`**).
- **Cross-check passes**: `E₁` from Re-slope `= E₁` from the approach limit `= 1.000000`.

---

## §2 · Does the combination achieve anything neither does alone?

**Yes, one thing: a cross-validated `(E₁, Δ)`.** The Re-slope alone can in principle report a spurious level
if a continuum artifact mimics a rung; the approach-law confirms the **same** `E₁` is the true approach limit
(rejecting the spurious reading), and jointly reads the gap `Δ`. So the combination is a **self-consistent
two-parameter estimator**, strictly stronger than either single read (it pins the atom's existence-at-`E₁`).

**But it does not certify, and the reason is sharp:**
- `E₁` is the **pole location** = A3 (read-not-force).
- **At finite `L`** the transfer matrix is discrete ⇒ the approach is **always exponential** with `Δ` = the
  **level spacing** (positive for any finite discrete spectrum). So the "exponential vs power-law" distinction
  is **not** visible at finite `L`; it is an `L→∞` statement. The **power-law** (= no isolation) only emerges
  when the ladder spacing `Δ→0` as `L→∞`.
- Therefore the **isolation** — "does `Δ` stay bounded away from 0 as `L→∞`?" — is precisely `I2 = H1`. The
  combination estimates the *finite-`L`* spacing; whether it survives the continuum limit is the wall.

---

## §3 · Honest bottom line

- **Genuine gain:** the two lines combine into a **cross-validated (mass, gap) diagnostic**, robust to a
  background continuum — stronger than either read alone (the cross-check pins the atom's existence at `E₁`).
- **NOT a certifier:** `E₁`=A3, `Δ`=spectral read; at finite `L` everything is discrete/exponential, and the
  isolation (`Δ` bounded as `L→∞`) = `I2 = H1`. So the combination reports `(E₁, Δ_finite-L)` and leaves the
  continuum-limit survival untouched.
- **Useful as a lattice protocol:** measure Re-slope `E₁` and the approach-law `Δ` from real-τ data; the
  cross-check rejects spurious levels; but a decisive statement still needs the uniform-`(L,a)` survival = `H1`.
- **Clay = 0%.** No QED.

*Refs:* `WI-MODIFIED-FUNCTIONAL.md` (Re-slope), `WI-PLATEAU-RIGOR.md` (P3 approach-law), `WI-EMBEDDED.md`,
`FALSIFY-*`/`TOOL-*`; witness `w7/wi_reslope_combine.py`.
