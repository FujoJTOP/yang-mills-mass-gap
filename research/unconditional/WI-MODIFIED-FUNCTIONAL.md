# WI-MODIFIED-FUNCTIONAL — a modified analytic functional that IS robust to a background continuum: the **log-slope of Re τ along a reweighting-persistent `|Im τ|` rung** separates a true atom (gapped *or* at threshold) from a pure cut, where the raw nearest-zero `𝒜` fails

> **VERDICT (no QED; Clay 0%).** The raw IR-analytic functional `𝒜=π/|Im τ_zero|` fails once a continuum is
> present (WI-ADMISSIBLE T2: a 1% continuum moves the *nearest* zero). **A modified functional fixes this:**
> an isolated atom puts its zeros on one vertical line `Re τ = ln(c₁/w₀)/E1`, so along each `|Im τ|` **rung**
> the real part moves **log-linearly** in the atom weight `c₁` with slope `dRe τ/d ln c₁ = 1/E1`. Under
> operator reweighting (levels fixed, weights changed) the rung **position** `|Im τ|=(2n+1)π/E1` is
> *invariant* while `Re τ` slides at that fixed slope — a **two-point signature** `(rung invariant, Re slope
> = 1/E1)` that a continuum does **not** mimic (its zeros have slope ≪ 1/E1 in both Re and Im). **Witness
> (exact Newton zeros): the atom signature survives a gapped continuum AND an at-threshold continuum —
> ATOM+GAP has 5 atom-slope rungs, ATOM+THR has 6 — while a pure CUT has 0 (slopes only 0.08–0.31).** This
> is a genuine hardening of the diagnostic over the raw `𝒜`. **But** it reads the atom's **pole log-slope =
> A3 (read-not-force)**; it needs a reweighting family + multi-zero resolution (data-hungry); and proving it
> **uniform in `(L,a)` at `β≈2.5` is `H1`.** A better *diagnostic*, not an escape from the wall.

---

## §0 · Labels

**[T]** theorem · **[C]** conditional · **[O]** exact missing step · **[NG]** no-go · **[V-num]** numeric (witness `w7/wi_modified_functional.py`).

---

## §1 · The failure being fixed

Finite-`L` 0⁺⁺ correlator: `Z(τ)=w₀+Σ_k c_k e^{−E_k τ}`, `w₀=|⟨O⟩|²>0`, zeros in `ℂ`.

- **Raw functional** `𝒜=π/|Im τ|` of the **nearest** zero. For a pure atom `w₀+c₁e^{−E1τ}`: `𝒜=E1`
  exactly, weight-invariant (`|Im τₙ|=(2n+1)π/E1`; weight only moves `Re τ`). **But** the physical KL
  measure is **atom + continuum**, and `WI-ADMISSIBLE` (T2) proved **any** second level — even a 1%-weight
  one — creates a **new zero near the atom's cancellation point** at small `|Im|`, so the *nearest* zero is
  a **spurious continuum artifact**, not the atom's. Raw `𝒜` therefore measures a continuum scale, not the
  atom. (`CONSTRUCT-IR-ANALYTIC-NUM` atom `𝒜=1.000000` vs `:φ²:` misfire `𝒜≈2.9` = this failure.)

---

## §2 · The modified functional: the reweighting-persistent rung's log-slope

**Reweighting** = changing the interpolating operator `O`; this changes the **overlaps** `{c_k}` (and `w₀`)
but **not the levels** `{E_k}` (`H` is fixed). So reweighting is a family `Z_c(τ)=w₀+Σ_k c_k e^{−E_k τ}`
parametrized by `c` (in practice, a smearing family).

**The functional** `ℒ` (atom log-slope along a rung), computed as:
1. For each reweighting `c`, find the zeros; collect the distinct **rungs** `|Im τ|`.
2. Keep rungs **persistent** across the reweighting family (same `|Im|` within tolerance, points 1 & 3).
3. Along a persistent rung, fit `Re τ` vs `ln c₁`; the **slope** is the functional.

**`ℒ` is atom-like** when there exists a persistent rung with `dRe τ/d ln c₁ = 1/E1` (and `π/|Im| = E1`).

---

## §3 · Why it is robust (the analytic reason) `[T]`

A **simple pole** of the Stieltjes transform `Π(z)=∫ρ(s)/(z+s)ds` at `s=−E1` (an isolated atom
`f₁²δ(s−E1)`) forces the atom's zeros to `τₙ=(1/E1)[ln(c₁/w₀)−i(2n+1)π]`:

- **rung position** `|Im τₙ|=(2n+1)π/E1` — depends only on `E1`, **invariant under `c₁`**;
- **real part** `Re τₙ=ln(c₁/w₀)/E1` — **log-linear in `c₁` with slope `1/E1`**.

A **branch cut** (continuum) contributes zeros whose `Re τ` moves with the density but **not** at the
atom slope `1/E1` (empirically 0.08–0.31 in the witness). So the **pair** `(rung invariant, Re-slope=1/E1)`
is a signature of the **pole**, which a cut cannot counterfeit. Because the atom's *deep* rungs are far from
the continuum-created spurious zeros (which sit at small `|Im|`), the signature survives a gapped **and** an
at-threshold continuum.

*(This is the A3 pole structure read off the zero ladder — see §6 for the honest consequence.)*

---

## §4 · Witness results `[V-num]` (exact Newton zeros; `w₀=1`, `E1=1`, reweight `c₁∈{0.4,1,2.5}`)

| model | persistent rungs | `dRe τ/d ln c₁` range | **atom-slope rungs** (`|s−1|<0.12`) |
|---|---|---|---|
| **ATOM** (isolated) | 6 | [1.000, 1.000] | **6** |
| **ATOM+GAP** (continuum edge `m_c=2>E1`, `ε=1%`) | 22 | [−0.323, 1.001] | **5** (at `|Im|=3π,5π,7π,9π,11π`) |
| **ATOM+THR** (continuum edge `m_c=E1`, `ε=1%`) | 19 | [−0.473, 1.001] | **6** |
| **CUT** (no atom) | 9 | [0.076, 0.305] | **0** |

- The **raw nearest-zero** `𝒜` degrades badly: ATOM `1.000`; ATOM+GAP `≈2.9`; ATOM+THR `≈2.04`; CUT `≈2.9`
  (indistinguishable from the atom cases) — the failure `WI-ADMISSIBLE` flagged.
- The **modified functional** `ℒ` **separates cleanly**: a true atom riding on a continuum (gapped or
  at-threshold) still exposes **5–6** atom-slope rungs; a pure cut exposes **0**.

The atom's *first* rung (`|Im|≈π`) is occasionally perturbed (ATOM+GAP: rung `3.112`, slope `0.324`), but
the **deep** rungs (`3π, 5π, …`) are clean — so `ℒ` is read off the **ladder**, not the single first zero.

---

## §5 · What it enables

A concrete lattice diagnostic: **reweight the 0⁺⁺ operator (smearing family), resolve the zero ladder,
and test for a persistent rung with log-slope `1/E1` and rung position `π/|Im|=E1`.** Positive ⟺ a genuine
atom (pole) is present, with its mass `E1` read off the rung — **robust to a continuum background**, which
the raw `𝒜` is not. This supersedes the `𝒜`-invariance test of `WI-SU3` (which detected only single-level
dominance).

---

## §6 · Honest limits `[O]/[NG]`

- **It is A3 (read-not-force).** The signature `(rung invariant, Re-slope=1/E1)` *is* "a simple pole of `Π`
  at `s=−E1`" — the Stieltjes pole `f₁²=−Res Π`. Reading the pole location/residue is `A3`; the modified
  functional changes the **representation** (zero ladder instead of one zero), **not the content**. It is
  the **A3 coordinate** of `H1`.
- **It is data-hungry.** It needs (a) a reweighting family (smearing), (b) resolving **multiple** complex
  zeros (a ladder), not one nearest zero — strictly harder than `WI-SU3`'s spread (which needed
  `10⁴–10⁶` configs and was gated by the unknown `Δ_cut`). Whether realistic SU(3) data resolve ≥3–4 deep
  rungs is **open**.
- **Toy calibration only.** Synthetic exact exponential sums; the at-threshold case is a *separate delta* at
  the cut edge (not a delta truly merged into the cut — a delta with `f₁²→0` would not register).
- **Uniformity is `H1`.** Proving this rung-log-slope signature **uniform in `(L,a)` at `β≈2.5`** (the pole
  present and readable in the continuum limit) is exactly `H1`. **Clay 0%, no QED.**

*Files:* `research/unconditional/WI-MODIFIED-FUNCTIONAL.md`; witness `research/unconditional/w7/wi_modified_functional.py`
(runs `py -3.12 -I`, numpy only). Prior: `WI-ADMISSIBLE.md`, `WI-ESCAPE.md`, `WI-THEOREM.md`, `WI-GENERAL.md`,
`WI-SU3.md`, `WI-REALTIME.md`, `CONSTRUCT-IR-ANALYTIC-NUM.md`.
