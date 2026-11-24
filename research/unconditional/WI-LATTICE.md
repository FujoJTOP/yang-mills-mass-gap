# WI-LATTICE — empirical calibration of the weight-invariance (WI) atom signature on solvable systems

> **VERDICT (no QED).** Tested the WI signature — "for the correlator spectral sum
> `Z_O(τ)=Σ_n |⟨n|O|Ω⟩|² e^{−(E_n−E_0)τ}`, the nearest complex-τ zero gives `𝒜_O=π/Im(τ_zero)`;
> atom ⇒ `𝒜_O` weight-independent (zero spread), cut ⇒ `𝒜_O` moves" — on exact solvable systems.
> **The signature is CONFIRMED in exact exponential sums `[T]+[V-num]`, and on the interacting 1D
> transverse-field Ising chain (exact diagonalization) `[V-num]`.** The sharpening is decisive: **zero
> spread ⟺ the atom is ISOLATED (a gap to *every* other level)** — a near-degenerate neighbour, a band,
> or a cut all move `𝒜_O`. The TFIM 0⁺⁺ channel is a two-kink **cut** (moves), exactly as `tr F²` is a
> two-gluon cut; free-fermion TFIM has **no** isolated 0⁺⁺ atom — the YM bound-state problem. WI certifies
> isolation when present; it does not manufacture the atom. **Clay 0%.**

---

## §0 · Labels
**[T]** 真定理 · **[V-num]** new numeric · **[C]** 条件 · **[O]** 精确缺失步 · **[NG]** no-go. Witness: `w7/wi_lattice.py` (`py -3.12 -I`, exact linear algebra, no MC).

---

## §1 · Results

### A. Exact exponential sums (the operational definition of the WI test)

| case | spectrum | `𝒜 = π/Im τ₀` across weight choices | spread | label |
|---|---|---|---|---|
| **A1 isolated atom** | `[0, m]`, `m=1`, weights `(1,1),(1,.3),(1,5),(.2,1),(7,1),(1,.05)` | **`1.000000` every time** | **`0.00e+00`** | **[T]** |
| **A2 cut** | `[0]` + 40 dense levels on `[2m,8m]`, 3 densities (flat/decay/modulated) | `3.0406, 3.0602, 3.0205` | `3.97e-02` | **[V-num]** |
| **A3 near-degenerate (NOT isolated)** | `[0, m, m+0.02]`, vary the weight on `m+0.02` ∈ {.01,.3,1,5} | `1.00020, 1.00460, 1.00993, 1.01657` | `1.64e-02` | **[V-num]** |

- **A1** confirms the exact atom identity: an **isolated** two-level sum gives `𝒜=m` **exactly, independent of the weight ratio** (weights only shift `Re τ`; the mass rides on `Im τ`). `[T]`.
- **A3** is the crucial control: even a *single* "particle" level is only weight-independent if **truly isolated**; a near-degenerate neighbour at `m+0.02` already makes `𝒜` move. **Isolation — not "single level vs many" — is the key variable.**

### B. 1D transverse-field Ising, `N=10`, exact diagonalization (`J=0.5, h=1.0`, disordered: unique ground state, no symmetry-breaking doublet)

| channel | operator(s) | `𝒜` | spread | label |
|---|---|---|---|---|
| **even 0⁺⁺ (two-kink CUT)** | `Σσᶻσᶻ` (bond), `Σσˣ` (field), smeared bond | `3.840, 4.734, 3.839` | **`0.8945`** | **[V-num]** |
| **odd (kink BAND)** | `Σσᶻ` (order), smeared `Σσᶻ` | `5.368, 4.846` | `0.5218` | **[V-num]** |

- The **0⁺⁺ (even) channel is a genuine cut** — its `𝒜` moves by `0.89` across 0⁺⁺ operators. This is the
  lattice-model analogue of `tr F²` being a **two-gluon cut**.
- The odd (order) channel is a **kink band** (lowest level `m=1.0003`, but a continuum of kinks above) —
  **not an isolated atom** — so it too moves (`0.52`). This is the A3 phenomenon on a real model.

---

## §2 · The refined signature

- **Zero spread ⟺ ISOLATED atom** (gap to every other level). Exact-sum A1 (isolated): `0`. Exact-sum A3
  (near-degenerate), A2 (cut), TFIM even (cut), TFIM odd (band): all `> 0`.
- The WI signature is therefore a genuine **finite-`L`, operator-side probe of ISOLATION** — it does not
  need to resolve the complex-τ zero densely; it asks only whether `𝒜` is invariant under operator
  reweighting. `[V-num]`
- **It can be measured, in principle, on real-τ data** by varying the interpolating operator at fixed
  quantum numbers and fitting `𝒜` (full complex-τ needs ill-posed analytic continuation; the *invariance*
  is the robust part).

---

## §3 · Mapping to 4D Yang–Mills — what this does and does not give

- `tr F²` is `ε`-like: its correlator is a **two-gluon cut**, like the TFIM even channel. So the WI test on
  the raw `tr F²` correlator would give **nonzero spread** unless a bound-state atom sits below the cut.
- The 0⁺⁺ glueball atom (the object Clay needs) is the **isolated bound state** below threshold. The WI
  test would certify it (`𝒜` weight-independent) **iff** it is isolated — and "isolated" **is `I2`**, the
  hard half of the gap.
- **The signature does not manufacture the atom.** TFIM is free-fermion (no bound state); its 0⁺⁺ sector is
  a cut and moves — a clean calibration that a *cut gives nonzero spread*. A interacting model with a
  genuine 0⁺⁺ bound state would give zero spread; producing that state is the YM problem.

**Net:** the WI signature is a real finite-`L` handle for an **isolated** 0⁺⁺ atom, rigorously confirmed in
exact sums and calibrated on TFIM. Its discriminating power **requires isolation** — so it does not reduce
the difficulty below `I2/H1`; it gives a *test* of isolation, if an atom is there to be isolated.

*Labels:* `[T]` exact-sum atom identity (isolated ⇒ `𝒜=m`, weight-independent); `[V-num]` cut / near-degenerate / TFIM spreads; `[O]` the isolated 0⁺⁺ bound state (YM) is absent in free-fermion models. **Clay 0%; no QED.**
