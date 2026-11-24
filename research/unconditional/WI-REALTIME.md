# WI-REALTIME — recovering the IR-analytic atom signature from real-τ data: the robust discriminator is the *reweighting-spread*, not the raw `𝒜`

> **VERDICT (no QED).** On a lattice one measures `Z(τ)` at **real** `τ∈[0,β]`, not complex τ. We build a
> practical estimator — matrix-pencil/ESPRIT (Prony) on a Hankel matrix → modes `(c,E)` → Newton zeros of
> the recovered exponential sum in complex τ → nearest `|Im τ|` → `𝒜=π/|Im τ|` — and test it on mock
> multi-exponential data with noise. **Result:** `𝒜` is recovered **robustly for an atom** (few isolated
> levels: 0.07% error at 1% noise, ~9% at 5%) with a **near-zero reweighting-spread (~0.01)**; but for a
> **dense continuum** the raw value is **ill-conditioned** (Hankel condition number explodes `~3↦6·10⁵↦6·10¹³↦3·10¹⁶`
> from atom to `p=16`) and **unreliable**, **while the reweighting-spread stays unambiguously large**. ⇒ the
> **robust finite-L atom signature is the SPREAD (invariance of `𝒜` under operator reweighting), not the raw
> `𝒜` value** — this survives the real-τ reconstruction and is the practical handle. `[V-num]` on synthetic
> data, not a lattice measurement. Clay 0%.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The estimator

`Z(τ)=w₀ + Σᵢ fᵢ² e^{−Eᵢ τ}` (real τ). Procedure:
1. **Hankel / ESPRIT**: build `H[i,j]=y[i+j]`, SVD, shift-invariance `U₁Ψ=U₂`, eigenvalues `zⱼ`; `Eⱼ=−log zⱼ/Δτ`; solve the Vandermonde `[zⱼⁿ]c=y` for the weights `cⱼ`.
2. **Zero-finder**: Newton on the recovered exponential sum in complex τ, over a grid; take the nearest `|Im τ|>0`; `𝒜=π/|Im τ|`.
   *(Fix applied: snap numerically-nonzero "constant" modes `|E|<1e−9` to exactly 0 and reject zeros with `|Re τ|>10` — otherwise `Z` underflows to a spurious zero at `Re τ∼10¹⁶`.)*
3. **Weight-invariance**: reweight the operator (change `cⱼ`, keep `Eⱼ`); recompute `𝒜`; the **spread** is the signature.

Mock: `Z=w₀+Σ fᵢ²e^{−Eᵢτ}`, `β=12`, `N=120` samples, Gaussian noise `∝ mean|Z|`.

---

## §2 · Recovery vs noise (ATOM: `w₀=1, f₁²=0.5, m₁=1`, true `𝒜=m₁=1`)

| noise | `𝒜` recovered | error |
|---|---|---|
| 0% | `1.0000 ± 0.0000` | 0.0000 |
| 0.1% | `0.9996 ± 0.0010` | 0.0004 |
| 1% | `1.0073 ± 0.0220` | 0.0073 (0.7%) |
| 5% | `1.0933 ± 0.2360` | 0.0933 (9%) |

**The atom is recoverable from real-τ data** to sub-percent at 1% noise. `[V-num]`

---

## §3 · Recovery vs #exponentials (CUT: `w₀=1`, continuum `[1,4]`, true `𝒜`=2.76/2.59/2.54 for p=4/8/16)

| p | noise | `𝒜` recovered | note |
|---|---|---|---|
| 4 | 0% | 2.7598 | correct |
| 4 | 1% | 21.6 ± 34 | **fails** |
| 8 | 0% | 2.5905 | correct |
| 8 | 1% | 12.3 ± 1.5 | **fails** |
| 16 | 0% | 13.39 | **fails even noise-free** |

**The raw `𝒜` is unreliable for a dense continuum.** Reason (section §4): recovering many closely-spaced exponentials from finite real-τ data is exponentially ill-conditioned.

---

## §4 · Conditioning (Hankel `S₀/S_{p−1}`, noise-free)

| spectrum | `S₀/S_{p−1}` |
|---|---|
| atom (2 modes) | `3.4×10¹` |
| cut p=4 | `6.3×10⁵` |
| cut p=8 | `6.0×10¹³` |
| cut p=16 | `3.4×10¹⁶` |

The condition number grows exponentially with the number of near-degenerate modes ⇒ the raw `𝒜` (nearest-zero distance) inherits this ill-conditioning for a continuum.

---

## §5 · The robust discriminator: reweighting-spread (ATOM vs CUT)

| object | reweight | `𝒜` values | **spread** |
|---|---|---|---|
| **atom** `f₁²∈{0.2,0.5,1.5,5}` | 1% noise | `0.990, 1.000, 1.001, 1.000` | **0.011** |
| **cut p=8** random `fᵢ²` | 1% noise | `6.7, 10.8, 138, 1933` | **1927** |

Even when the raw `𝒜` is ill-conditioned, the **spread under operator reweighting separates atom from cut by ~5 orders of magnitude**. This is the practical handle: **an isolated atom gives an `𝒜` invariant under reweighting; a continuum gives an `𝒜` that moves.**

---

## §6 · Honest bottom line

- **`[T]`/`[V-num]`**: for an isolated atom, `𝒜=π/|Im τ_zero|=m₁` exactly, and real-τ data recovers it to sub-percent at 1% noise; the reweighting-spread is ~0.
- **`[NG]`**: the raw `𝒜` is **not** a reliable continuum/atom discriminator from finite real-τ data (ill-conditioned: Hankel condition `≥10⁵` for p≥4).
- **`[V-num]` (the practical handle)**: the **reweighting-spread** is the robust finite-L atom signature (atom ~0.01, cut ~10³), surviving real-τ reconstruction and moderate noise.
- **Not a proof, not a lattice measurement** — synthetic multi-exponential mock data. The real test is a lattice `0⁺⁺` correlator with several reweighted operators; proving the spread is uniformly small at `β≈2.5` (with reweighting that preserves coupling content) is `H1`.
- **Clay 0%; no QED.**

*Refs:* `CONSTRUCT-IR-ANALYTIC-NUM.md` (the weight-invariance signature found there), `PNT-IR-ANALYTIC.md`, `PNT-OPE-RHO.md` (finite-data ill-posedness). Witness: `w7/wi_realtime.py`.
