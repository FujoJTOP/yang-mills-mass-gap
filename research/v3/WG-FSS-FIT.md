# WG-FSS-FIT — falsifiable finite-size law for the gauge-invariant `χ_O`

> **2026-11-24 · Wave WG, continuation.** Fits the five gauge-invariant points
> `χ_O = 0.328, 0.377, 0.411, 0.469, 0.587` at `L = 4,6,8,12,16` (same estimator / measure /
> APE `ns=2`, `α=0.5`, `β=2.5`, SU(2)) to `χ_O(L) = χ_∞ + c/L^p`, and compares `p=2`, `p=1`,
> `p=0` (plateau), the free-`p` fit, and two log laws. **Witness:** `research/v3/w7/wg_fss_fit.py`
> (`py -3.12 -I`, numpy only, one process, no downloads). Data source: `WG-L12-L16.md` §2 /
> `w7/wg_large_results.json`.
> **Labels (charter §2).** **[V]** numerical witness · **[C]** conditional reading · **[O]** exact
> missing step. **No QED. Clay 0%.** This file is **updatable**: when L=24/32 land, replace the
> five-point arrays in the witness and re-run — no other change needed (§6).

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **`[V]` The pure plateau is excluded.** `χ_O = const` gives `χ²/dof = 29.1/4 = 7.3`
   (p-value ~ 7×10⁻⁶). A flat `χ_O` does not describe the five points.

2. **`[V]` A saturating law with `p ≤ 1` fits, and in that fit `χ_∞ > 0` at high significance —
   but the value of `χ_∞` is model-dependent.** `p=1` fixed: `χ_∞ = 0.556 ± 0.035`
   (`χ²/dof = 4.03/3 = 1.34`, acceptable). `p=2` (the `1/L²` FSS correction) is marginal-poor:
   `χ_∞ = 0.482 ± 0.023`, `χ²/dof = 2.35`. **Given** the saturating ansatz, `χ_∞/se = 16 (`p=1`)
   to 21 (`p=2`) — i.e. `χ_∞>0` *within the model*.

3. **`[V]/[C]` `χ_∞>0` is NOT robust across models — the five points do not pin it.** Models with
   **`χ_∞ = 0` fit as well or better**: a pure power `χ_O = A L^q` gives `q = 0.371 ± 0.066`,
   `χ²/dof = 0.39`; a log law `A + B ln L` gives `B = +0.148 ± 0.028`, `χ²/dof = 0.56`. The
   free-`p` 3-parameter fit is **degenerate** (minimum runs to the `p→0` boundary; `χ_∞`
   unidentifiable). Honest statement: **the data are consistent with a finite thermodynamic
   `χ_∞>0` *or* with unbounded slow growth (`L^0.37` / log); they exclude only the plateau and
   (marginally) the steep `1/L²` approach.**

**One line.** *Five points of `χ_O(L)` reject a flat plateau (`χ²/dof=7.3`); a saturating
`χ_∞=0.56±0.04 / L` law fits at `χ²/dof=1.3`, but a pure power `L^0.37` and a log law (both
`χ_∞=0`) fit equally well, so `χ_∞>0` is supported only conditionally — the data cannot yet
separate "saturates to a finite value" from "grows forever, slowly". L=24/32 will decide. No QED.
Clay 0%.*

---

## 1 · Data (identical protocol to `WG-L12-L16.md`)  [V]

| L | `χ_O` | SE | pull from weighted mean |
|---|---|---|---|
| 4  | 0.3277 | 0.0223 | −4.8 |
| 6  | 0.3774 | 0.0213 | −2.7 |
| 8  | 0.4105 | 0.0576 | −0.4 |
| 12 | 0.4688 | 0.0279 | +1.2 |
| 16 | 0.5865 | 0.0547 | +2.8 |

Weighted mean `0.393 ± 0.013`, reduced `χ² = 7.3` — the large `L` points sit above, the small
below: a monotone rise, not a plateau.

---

## 2 · Fits  [V]

Weighted LS, `w = 1/SE²`; ansatz `χ_O(L) = χ_∞ + c/L^p` unless noted. Errors from the inverse
Hessian (formal; `χ²/dof > 1` means the SEs are underestimated — see §4).

| model | `χ_∞` | `p` | `c` | `χ²/dof` | p-value |
|---|---|---|---|---|---|
| **pure plateau** `χ=const` | `0.393 ± 0.013` | 0 | — | **7.28** | 7×10⁻⁶ |
| **`p=2`** (`1/L²` FSS corr.) | `0.482 ± 0.023` | 2 (fix) | `−2.69 ± 0.57` | 2.35 | 0.07 |
| **`p=1`** (`1/L` corr.) | `0.556 ± 0.035` | 1 (fix) | `−0.963 ± 0.192` | 1.34 | 0.26 |
| **free `p`** (grid ≥0.05) | `3.36 ± 0.57` ⚠ | `0.05 ± 1.17` ⚠ | `−3.26 ± 0.62` | 0.58 | 0.63 |
| **free `p`** (grid ≥0.40) | `0.777 ± 0.075` | `0.40 ± 1.26` | `−0.80 ± 0.16` | 0.82 | 0.48 |
| **saturating log** `χ_∞+c/lnL` | `0.676 ± 0.058` | — | `−0.50 ± 0.10` | 1.38 | 0.25 |
| **pure power** `A L^q` (`χ_∞=0`) | 0 (exact) | `q=0.371 ± 0.066` | `A=0.194` | **0.39** | 0.76 |
| **log law** `A + B ln L` | undefined | `B=+0.148 ± 0.028` | `A=0.117 ± 0.054` | 0.56 | 0.64 |
| *reference* `χ_∞ + cL²` (`p=−2`, free growth) | `0.326 ± 0.018` | −2 (fix) | `+0.0010 ± 0.0002` | 0.35 | 0.79 |

⚠ **Degenerate:** the free-`p` minimum sits on the `p→0` grid boundary. As `p→0` the shift
`c/L^p → c` becomes `L`-independent, `χ_∞` and `c` are perfectly anti-correlated, and `χ_∞`
diverges. **Do not read `χ_∞ = 3.36` as physical.** The free-`p` fit therefore does **not**
determine `p` (`p = 0.05 ± 1.17`, i.e. anything from 0 up); it only shows the data prefer the
*flattest* allowed correction.

**Two clarifications on "p=2 (free)".**
- `p = +2` (correction `c/L²`) is the standard free-field / bulk FSS correction to a bounded
  susceptibility. It is **marginal**: `χ²/dof = 2.35` (p-value 0.07) — it undershoots the rise
  from `L=12→16` (it predicts `χ→0.482`, data give `0.587`).
- Genuine **free-massless growth** `χ_O ∝ L²` is the *different* object `χ_∞=0, q=2` in the pure-
  power family — already excluded by the log-log slope `0.37 ± 0.07` (`24σ` from 2, `WG-L12-L16.md`).
  The `p=−2` row is included only to expose a degeneracy: with a free additive constant, `χ_∞+cL²`
  can also mimic the five points (`χ²/dof=0.35`) because `c` is tiny and the constant carries the
  bulk. It is **not** evidence for `L²` growth (the through-origin `L²` law is dead).

---

## 3 · Reading: is `χ_∞ > 0` supported?  [V]/[C]/[O]

| claim | status | evidence |
|---|---|---|
| flat plateau (`χ_∞` const) | **excluded [V]** | `χ²/dof = 7.28`, p ~ 7×10⁻⁶ |
| `χ_∞>0` *given* `p=1` or `p=2` ansatz | supported **within model** [C] | `χ_∞/se = 16`–`21` |
| `p=2` (`1/L²`) approach | **marginal / disfavoured** [V] | `χ²/dof = 2.35` |
| `p=1` (`1/L`) approach | acceptable [V] | `χ²/dof = 1.34` |
| `χ_∞ = 0` (pure power `L^0.37`) | **equally consistent [V]** | `χ²/dof = 0.39`, `q=0.371±0.066` |
| `χ_∞ = 0` (log `A+B lnL`) | **equally consistent [V]** | `χ²/dof = 0.56`, `B>0` (diverges) |
| `χ_∞>0` **robustly** (model-independent) | **NOT established** [O] | 5 points admit `χ_∞=0` fits at equal/better `χ²` |

**Honest bottom line.** The only unconditional `[V]` statements are: (i) *not a plateau*; (ii) *not
free `L²` through the origin*. Everything about `χ_∞` is **conditional on the functional form**.
Imposing a saturating `1/L^p` (`p≤1`) yields `χ_∞ ≈ 0.56–0.78 > 0` at >10σ; but a pure power
`L^0.37` or a log law — both with `χ_∞=0`, i.e. **no finite thermodynamic limit** — fit at least as
well. With five points and `2–3` free parameters (and errors this large) the saturation regime is
**not reached**, so `χ_∞>0` cannot be claimed as a model-independent result. This is the **exact
missing step**: more volumes to break the degeneracy — `[O]`.

*Discriminator for the next round.* The viable families spread apart only at large `L`. Predicted
`χ_O` at `L=24,32` from the fits above:

| model | `χ_O(24)` | `χ_O(32)` |
|---|---|---|
| `p=1` saturating | 0.516 | 0.526 |
| `p=2` saturating | 0.477 | 0.479 |
| sat-log `χ_∞+c/lnL` | 0.519 | 0.532 |
| pure power `L^0.37` | 0.630 | 0.701 |
| log `A+B lnL` | 0.586 | 0.629 |

The saturating family clusters at `≈0.48–0.53`; the `χ_∞=0` family (power/log) at `≈0.59–0.63`
(L=24) rising to `≈0.63–0.70` (L=32). At the present SE (`~0.05`) a single `L=32` point separates
the two families at `≈2–3σ`. That is the falsifiable prediction.

---

## 4 · Caveats (formal errors, dof)  [O]

- **Formal errors.** With `χ²/dof` as low as 0.35–0.8 the inverse-Hessian SEs are optimistic (the
  model is over-fitting 5 points); with `χ²/dof>1` (p=2) they are conservative at the fitted level
  but do not capture the systematic (wrong functional form). Treat all `χ_∞±err` as *conditional*.
- **Degeneracy.** `(χ_∞, c, p)` are strongly correlated; `p` is bounded only from above by the
  rise. The free-`p` fit is unusable for `χ_∞` (§2).
- **Fixed `a`, fixed APE radius.** Same lattice spacing / smearing radius as the data; the `L`-law
  is at one `a`, not a continuum extrapolation (`WG-L12-L16.md` §4).
- **Not H1 / not the gap.** `χ_O = ∫s⁻¹dν_O` is the squeeze ceiling, never the atom; `f₁²≥r₀`,
  isolation `I₂`, and uniform-in-`L` `χ_O≤1/λ` remain `[O]`. **No QED. Clay 0%.**

---

## 5 · Status ledger

| item | statement | label |
|---|---|---|
| plateau excluded | `χ_O=const`: `χ²/dof=7.28`, p~7e-6 | **[V]** |
| `p=1` saturating fit | `χ_∞=0.556±0.035`, `χ²/dof=1.34` | **[V]** |
| `p=2` (`1/L²`) fit | `χ_∞=0.482±0.023`, `χ²/dof=2.35` (marginal) | **[V]** |
| free-`p` fit | degenerate (`p→0` boundary), `χ_∞` unidentifiable | **[V]** |
| pure power `L^0.37`, `χ_∞=0` | `q=0.371±0.066`, `χ²/dof=0.39` | **[V]** |
| log law `A+B lnL`, `χ_∞=0` | `B=+0.148±0.028`, `χ²/dof=0.56` | **[V]** |
| `χ_∞>0` model-independent | **not established**; `χ_∞=0` fits equally | **[O]** |
| next discriminator L=24/32 | sat vs power separated at ≳3σ | **[O]** |
| atom `f₁²≥r₀`; H1; gap | — | **[O]** |
| Clay mass gap | unbuilt | **0 %** |

---

## 6 · Reproduce / update

```
py -3.12 -I research/v3/w7/wg_fss_fit.py
```

To update when **L=24/32** land: append the two `(L, χ_O, SE)` triples to the arrays at the top of
`wg_fss_fit.py` (and to `wg_large_results.json` `trend.sequence`); re-run. The script refits every
model and reprints §2–§3. No other change. If `p` moves off the `p→0` boundary and `χ_∞` stabilises,
`χ_∞>0` upgrades from `[C]` to `[V]`; if the points keep climbing on the `L^0.37` line, the honest
reading flips to "no finite thermodynamic `χ_O` in this window" (a `[V]` statement, still not the gap).

*Refs.* Data: `WG-L12-L16.md`, `WG-INVARIANT-ESTIMATOR.md`; witness `w7/wg_fss_fit.py`, data
`w7/wg_large_results.json`. No `[C]` upgraded to `[T]`; no QED; **Clay 0%.**
