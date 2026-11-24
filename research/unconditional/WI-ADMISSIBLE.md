# WI-ADMISSIBLE — the reweighting family, and the honest verdict: the signature detects a *single level*, not an *atom above a continuum*

> **VERDICT (no QED; this DOWNGRADES the lead).** `CONSTRUCT-IR-ANALYTIC-NUM` reported "𝒜 invariant under
> reweighting ⇔ atom". Made precise, the invariance is **exact only for a pure single level + VEV** (T1).
> **Any additional spectral level, however small, moves the nearest complex zero substantially** (T2, verified
> numerically: a **1%-weight** second level gives `𝒜=π/|Im|≈3.156` versus `𝒜=E₁=1` for the pure atom,
> because at the atom's cancellation point the second term is the *leading nonvanishing* one and *creates a
> new zero* at `|Im|≈0.996<π`). The physical KL measure is **atom + continuum**, so 𝒜 will **drift** even
> for a true atom. Hence the test is a **"single-exponential-dominance" detector, NOT an atom-vs-continuum
> certificate**: it cannot certify the Clay atom sitting above a continuum. Uniformity in `(L,a)` at `β≈2.5`
> remains `H1` regardless. **Clay 0%.**

---

## §0 · Labels and the functional

**[T]** theorem · **[C]** conditional · **[O]** exact missing step · **[NG]** no-go · **[V-num]** numeric.

Finite-`L` 0⁺⁺ correlator of `O`: `C_O(τ) = w₀(O) + Σ_i f_i²(O) e^{−E_i τ}`, with `w₀=|⟨0|O|0⟩|²`,
`f_i²=|⟨0|O|i⟩|²`, `E_i` the levels (masses). IR-analytic functional `𝒜(C) := π/|Im τ_zero|`, `τ_zero` the
zero of `C` in `ℂ` nearest the real axis. (`CONSTRUCT-IR-ANALYTIC-NUM`: `𝒜=m₁` for a pure atom; raw `𝒜` at
finite `L` is a *scale detector* that fires on cuts too — the atom signal was claimed to be **invariance
under operator reweighting**.)

---

## §1 · The three exact propositions

### (T1) Pure single level + VEV ⇒ `𝒜` exactly invariant
`C(τ)=w₀+f²e^{−Eτ}`, `w₀>0, f²>0, E>0`. Then `C=0 ⇔ e^{−Eτ}=−w₀/f²`, i.e.
```
τ_k = (1/E)[ ln(f²/w₀) − i(2k+1)π ],      |Im τ_k| = (2k+1)π/E,      𝒜 = π/|Im τ₀| = E.
```
`|Im τ|` and hence `𝒜` depend only on `E`, **not on `(w₀,f²)`**; the weights move `Re τ_zero` only.
*(One-line proof; witness §A: `𝒜=1.000000` for `f²∈{0.3,1,5}`; `Re τ_zero=ln(f²/w₀)/E` moves, `𝒜` fixed.)*

### (T2) Any second level moves `𝒜` — even a tiny one  ⟵ the decisive correction
For `C(τ)=w₀ + f₁²e^{−E₁τ} + f₂²e^{−E₂τ}` with `E₂>E₁>0`, the second exponential **oscillates faster in
`Im τ` (period `2π/E₂ < 2π/E₁`)** and, at the point where the atom cancels the VEV (`w₀+f₁²e^{−E₁τ}=0`),
the second term is the **leading nonvanishing** one — so it shifts the zero non-perturbatively and **carries
its own zeros** at `|Im|≲π/E₂ < π/E₁`. *Numerically verified* (§B/§C): `w₀=1`, atom `(f₁²=1,E₁=1)`, a
**1%-weight** level at `E₂=4` already yields `|Im τ_zero|=0.9956 < π` ⇒ `𝒜=3.156`, **not** `E₁=1`; the
nearest zero is a genuine new root (`u=2.63+4.05i` of the quartic). **So `𝒜` is NOT robust to a small higher
weight.**

### (T3) Source reweighting changes weights, never levels
Fix `H` (hence `{E_i}`). For `O,O′∈0⁺⁺`, `C_{O+cO′}(τ)=w₀(c)+Σ_i f_i²(c)e^{−E_i τ}` with
`f_i²(c)=|⟨0|O|i⟩+c⟨0|O′|i⟩|²`, **same `{E_i}`**. *(Linearity of overlaps; `[T]` at fixed `H`.)* This is an
exact weights-only dial — but by (T2) the *nearest zero* is not a level, so the dial does **not** leave `𝒜`
fixed once `≥2` levels carry weight.

---

## §2 · Consequence: what the test really detects

| spectrum | `𝒜` under reweighting |
|---|---|
| **pure single level + VEV** | **invariant** `=E₁` (T1) |
| **atom + any higher level/continuum** (physical case) | **drifts** — the higher component's faster oscillation creates closer zeros (T2) |
| **pure continuum** | drifts (density-dependent) |

⇒ **`𝒜`-invariance ⇔ "the low spectrum is effectively a single exponential"; drift ⇔ "≥2 components
carry weight".** The physical 0⁺⁺ is **atom + continuum**, hence in the **drift** class. **The test therefore
does not certify the Clay atom**; it certifies the idealized single-level limit only. (This corrects the
optimistic reading of `CONSTRUCT-IR-ANALYTIC-NUM`'s weight-invariance finding: the atom they tested was a
single level, not an atom above a continuum.)

**Is there a rescue?** Only a reweighting that changes the *atom weight* while **fixing the continuum
coupling** — i.e. `⟨0|δO|continuum⟩=0` for all continuum states, `⟨0|δO|atom⟩≠0`. Such `δO` is fine-tuned,
not generic; no operational family supplies it. So the invariance signature is **not** a usable atom
certificer for YM without that fine-tuning.

---

## §3 · The admissible family (still well-defined, for the narrower question)

Which reweightings change **weights** but not **levels** (T3), used to test *single-exponential dominance*:

- **(A1) same channel** `O,O′∈0⁺⁺` (`J^{PC}=0⁺⁺`, `p=0`, center-neutral) — else mixed basis (C4).
- **(A2) VEV no-zero** `⟨O+cO′⟩≠0` uniform sign (`w₀(c)>0`) (C1).
- **(A3) no node** `f₁(c)≠0`; detect jumps (C2).
- **(A4) O′ does not lower the continuum** below the candidate lightest level (C3) — but note: by (T2)
  *any* continuum, even above, spoils invariance, so (A4) is necessary-not-sufficient.
- **(A5) power** `f₁²(c)` varies.

**Useful reweightings:** `tr F²` at different smearing radius; different lattice discretisation of `F`;
an independent `0⁺⁺` interpolator. **Unusable:** anything changing the theory (different `β`/action —
changes levels, breaks T3); anything outside `0⁺⁺`.

**Confounds:** C1 VEV sign/zero; C2 nodes; **C3 continuum (now central: any continuum ⇒ false drift)**;
C4 channel mixing; C5 finite-`a` artifacts.

---

## §4 · Operational recipe (narrowed null vs alternative)

1. Fix `(β,a,L)`, fix `H`; vary only the source `O_c=O+cO′`.
2. Measure `C_{O_c}(τ)`; extract `𝒜(c)` (complex-τ zero, or a shared-level multi-exp fit across `c`).
3. **Decide:**
   - **"Single-exponential"**: `𝒜(c)` constant.
   - **"Multi-component"**: `𝒜(c)` drifts.
4. **Interpretation caveat (decisive):** a TRUE atom above a continuum is in the **multi-component** class,
   so a drift is **not** evidence against an atom, and invariance is **not** evidence for the physical atom
   — only for the idealized single level. The test localizes "is the low spectrum one level?", nothing more.
5. **Do NOT use real-τ effective mass as a discriminator** (a cut's edge is shared, so it looks invariant).

---

## §5 · Honest boundary (the lead is downgraded)

- **(T1) invariance is exact only for a pure single level.** For the physical atom + continuum, **(T2)**
  makes `𝒜` drift even for a genuine atom; the 1%-weight demonstration (`𝒜:1→3.156`) shows the effect is
  large, not a small correction.
- **Therefore `𝒜`-invariance is NOT an atom certificer for Yang–Mills.** It is a diagnostic of
  single-exponential dominance. The earlier "weight-invariance atom signature" is **retained only in the
  idealized sense**; as an operational atom test for the real theory it **fails** (it cannot separate a true
  atom from a single strong low resonance, and it misfires as "multi-component" on any atom + continuum).
- **What would be needed:** an extraction where the continuum is held fixed while the atom weight varies
  (requires a fine-tuned reweighting, absent generically), or a non-local functional that reads the atom
  (all such reduce to `H1`, per `PNT-OPE-RHO`/`PROMISING-SYNTHESIS`). Uniformity in `(L,a)` at `β≈2.5` is
  `H1` regardless.
- **Not a proof of anything about the gap.** This file's honest yield is negative: it removes the
  weight-invariance lead as a stand-alone atom signature.

*Refs:* `CONSTRUCT-IR-ANALYTIC-NUM.md` (empirical atom-invariance, single-level), `PNT-IR-ANALYTIC.md`
(`𝒜=π/Im τ`), `PNT-OPE-RHO.md` (finite integrals blind), `PROMISING-SYNTHESIS.md`. Witness `w7/wi_admissible.py`
(§A pure-atom invariance; §B/§C any-higher-level drift, verified to `|f|∼1e-15`).
