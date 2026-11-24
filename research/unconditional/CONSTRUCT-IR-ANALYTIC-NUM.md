# CONSTRUCT-IR-ANALYTIC-NUM — the IR-analytic functional `𝒜=π/Im(τ_zero)` made concrete: the atom identity is exact and weight-independent, **weight-invariance is a real finite-`L` atom signature**, but the raw value is a scale detector, and its uniform proof is `H1`

> **VERDICT (no QED).** `PNT-IR-ANALYTIC`'s functional is now *computed*, not just schematized. The clean core
> holds and is sharper than claimed; one part of the prior claim is **corrected**; and a genuine **new
> positive** appears. **(1) [T]** For an isolated `0⁺⁺` atom, the nearest complex-`τ` zero of
> `Z(τ)=w₀+Σfᵢ²e^{−Eᵢτ}` gives `𝒜=π/Im(τ_zero)=m₁` **exactly, independent of the weight ratio** (weights
> only shift `Re τ`); verified to `|Z|~10⁻¹⁶`. **(2) [V-num] CORRECTION:** at **finite `L` the raw value `𝒜`
> is a *scale detector*, not an atom certifier** — it fires on gapped-atomless `:φ²:` (`𝒜≈2.9`), on the free
> massless box (`𝒜≈3.2`, does **not** →0), and on a zero-VEV cut (`𝒜≈1.8`). The prior claim "fires iff
> `f₁²>0`, acid-clean" holds **only in infinite volume**. **(3) [V-num] NEW POSITIVE:** a **weight-invariance
> atom signature** — the atom gives `𝒜` with **zero spread** under operator reweighting, the cut gives `𝒜`
> that **moves** (3.53/2.94/3.85). This distinguishes atom from cut **at fixed `L`** — a concrete finite-`L`
> handle. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric · **[V-here]** new this round. Concrete witnesses in `w7/construct_ir_analytic_num.py` (runs `py -3.12 -I`, numpy only).

---

## §1 · What was built and computed

On a finite lattice `H` has a **discrete** spectrum, so the faithful finite-`L` object is the exponential sum
$$Z(τ)=\underbrace{w_0}_{=|\langle0|O|0\rangle|^2}+\sum_{i\ge1} w_i\,e^{-E_iτ},\qquad w_i=f_i^2\ge0,\ E_i=\text{0}^{++}\text{ excitation energies}.$$
We compute its complex-`τ` zeros by Newton from a grid (`numpy`), dedup, sort by `|Im τ|`, and evaluate
`𝒜=π/|Im(τ_nearest)|`. All results below are at machine precision (`|Z(τ_zero)|~10⁻¹⁶`).

The **`0⁺⁺` projection** enters through `Z(τ)=Tr[P_{0⁺⁺}\,O\,e^{−τH}]`-style overlaps; concretely one uses the
zero-momentum, even, connected correlator with its **VEV term `w₀`** — the zero is the **vacuum-vs-lightest-
excitation interference**, and it needs `w₀>0` (which `tr F²` has: the gluon condensate). *(This is the `(a)`
assessment: the `0⁺⁺` projection is the `w₀` term + 0-momentum even operator; without `w₀` the "vacuum vs gap"
zero is absent.)*

---

## §2 · Result `[T]` — the exact atom identity

| model | `w=(w₀,w₁)` | nearest zero `τ` | `𝒜=π/|Im τ|` | true `m₁` |
|---|---|---|---|---|
| atom `m₁=1` | `(1, 0.5)` | `−0.693 + 3.1416i` | **1.00000** | 1.0 |
| atom `m₁=1` | `(1, 3.0)` | `+1.099 + 3.1416i` | **1.00000** | 1.0 |
| atom `m₁=2` | `(1, 0.7)` | `−0.178 + 1.5708i` | **2.00000** | 2.0 |

**`Im(τ)=π/m₁` is exact and INDEPENDENT of the weight ratio** — the weights only move `Re(τ)=(1/m₁)\ln(f_1^2/w_0)`.
The atom-mass (times `π`) is carried by the **imaginary part**; the atom-**weight** by the **real part**. This is
the clean, renormalon-free content `PNT-IR-ANALYTIC` claimed — confirmed exactly.

---

## §3 · Result `[V-num]` — CORRECTION: at finite `L`, `𝒜` is a scale detector, not an atom certifier

A faithful **free massive scalar on an `N`-site ring**, `O=:φ²:` (whose `0⁺⁺` excitations are 2-particle,
`E=2ω_k`, and whose `1`-particle coupling `f₁²=0`), fed into `Z`:

| case | `𝒜` | should be (acid) | verdict |
|---|---|---|---|
| `:φ²:` cut, `m=1`, `N=8…32` | `2.90…2.96` | vanish (no atom) | **FIRES** ✗ |
| massless box, `N=8…64` | `3.41…3.27` | vanish (`G″=0`) | **FIRES, not →0** ✗ |
| zero-VEV pure cut, `N=8…32` | `1.84…1.70` | vanish | **FIRES** ✗ |

⇒ **At finite `L`, `𝒜` fires on atomless `:φ²:`, on the free massless field, and on pure cuts** — the
discriminating acid test **fails**. The zero is set by the *scale of the excitation content*, not the mass
gap. **This corrects `PNT-IR-ANALYTIC`**: "fires iff `f₁²>0`, acid-clean" is an **infinite-volume** statement;
at finite `L` the raw `𝒜` is one more gap/scale detector (same family as `SLOT-*`).

---

## §4 · Result `[V-num]` — NEW POSITIVE: weight-invariance is a finite-`L` atom signature

The **single most important computation.** Reweight the `0⁺⁺` operator (changes `{w_i}`, not `{E_i}`) and track `𝒜`:

| model | reweighting | `𝒜` | spread |
|---|---|---|---|
| **atom** `m₁=1` | `(1,0.3)`, `(1,1)`, `(1,5)` | `1.000000`, `1.000000`, `1.000000` | **`0.00×10⁰`** |
| **cut** (`:φ²:`) | uniform, `1/(ωω)`, `exp-tilt` | `3.526`, `2.938`, `3.853` | **`0.92`** |

⇒ **`𝒜` is weight-INVARIANT for an atom, weight-DEPENDENT for a cut.** This is a **concrete atom-vs-cut
discriminator that works at FIXED finite `L`** — not requiring the infinite-volume limit. The atom's zero is a
*stable, weight-independent* feature; a cut's zeros are an *unstable interference pattern*. The zero **SET**
confirms it: the atom gives a **sparse, periodic** set (`Im=π/m₁`, count `=5`), the cut a **dense, crowding**
set (count `=21`, growing with `N`).

**This is the one genuinely new handle from the construction:** *"`𝒜` invariant under operator-reweighting"*
reads the atom at finite `L`. It is the finite-`L` shadow of the atom's isolated pole.

---

## §5 · Assessment of the three questions

- **(a) `0⁺⁺` projection.** Use the zero-momentum even connected correlator with its VEV term `w₀=|⟨tr F²⟩|²`;
  the zero is vacuum-vs-lightest-excitation; `w₀>0` is essential (it holds for `tr F²`).
- **(b) `L→∞`.** The atom's zero is **fixed** (`Im=π/m₁`, sparse set, weight-invariant); the cut's zeros
  **crowd/migrate** (dense set, weight-dependent). The `L`-stable invariant is **weight-invariance**, not the
  raw value.
- **(c) Lattice measurability.** `𝒜` needs `Z(τ)` at **complex `τ`** = analytic continuation of the Euclidean
  correlator = reconstructing the exponential sum from real-`τ` data — an **ill-posed spectral reconstruction**
  (one more face of `H1`). **But** the weight-invariance test is partly measurable: **vary the `0⁺⁺` operator**
  (e.g. `tr F²` vs a smeared variant), fit `𝒜` at each, and test invariance — this reads the atom from real-`τ`
  correlators without full continuation.

---

## §6 · Honest bottom line

- **The exact identity `𝒜=m₁` for an isolated atom is `[T]`** and weight-independent (verified to `10⁻¹⁶`).
- **`PNT-IR-ANALYTIC` is corrected:** the raw `𝒜` is **acid-dirty at finite `L`** (fires on `:φ²:`, massless,
  cuts) — "fires iff `f₁²>0`" is infinite-volume only.
- **Genuine new positive `[V-num]`:** **weight-invariance of `𝒜` is a finite-`L` atom signature** — a concrete
  atom-vs-cut discriminator that does not need `L→∞`. This is real, and it was not previously identified.
- **What is still `H1`:** (i) reweighting must preserve *which states couple* (a smeared/alternate `0⁺⁺`
  operator that keeps the same spectrum content); (ii) proving the invariance *uniformly in `(L,a)`* at
  `β≈2.5` — i.e. that the atom's zero stays weight-invariant through the crossover and the continuum; (iii)
  the full complex-`τ` continuation on a real lattice.
- **Net:** `𝒜` is **more** than a mere re-encoding — its **weight-invariance** is a real finite-`L` atom
  signature — but turning that signature into a uniform certificate is `H1`. **Clay 0%; no QED.**

*Witness:* `w7/construct_ir_analytic_num.py` (sections A–E, machine-precision zeros). Not committed (parent stages). Related: `PNT-IR-ANALYTIC.md`, `INVENT-LEEYANG.md`, `SLOT-SYNTHESIS.md`, `PROMISING-SYNTHESIS.md`.
