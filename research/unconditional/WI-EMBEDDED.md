# WI-EMBEDDED — atom embedded in a continuum: nearest-zero signature is destroyed by ANY continuum (ε\*=0), the tracking version is a weak ~6% effect — WI-ADMISSIBLE is right for the observable that can be extracted

> **VERDICT (no QED).** The WI-LATTICE ("isolation suffices") vs WI-ADMISSIBLE ("1% continuum moves 𝒜")
> split resolves by separating two observables. **(1) Nearest-zero 𝒜** (what ESPRIT / standard
> extraction returns): destroyed by **any** continuum — even `ε=10⁻⁴` (0.01% weight) sends `𝒜` from `E₁=1`
> to `≈14`, indistinguishable from a pure cut (`≈14–15`). Threshold `ε* = 0`. This **decisively confirms
> WI-ADMISSIBLE**. **(2) Atom-tracking 𝒜** (the zero whose arg is closest to `π`): retains a root near `E₁`,
> but only **~6.5%** separated from a pure cut's value, shape-dependent, and comparable to its own
> reweighting wobble (`~6%`) — **not a robust discriminator**. The physical case (atom + two-particle
> continuum, `ε` O(1)-coupled at the crossover) therefore **fails both**: the atom is not cleanly
> certifiable. **Clay 0%.**

---

## §0 · Setup and metrics

`u=e^{−τ}`; `Z=w₀+Σc_n uⁿ` (polynomial, `c_n≥0`); `|Im τ|=|arg u|/ω`; `𝒜=π/|arg u|`. Atom at level 1
(`E₁=1`), weight `f²`; continuum at `n=2..16` (two-particle threshold), total weight `ε`; VEV `w₀=1`.
Pure atom: zero at `arg=π` ⇒ `𝒜=E₁=1` for **all** `f²` (spread 0) — both metrics confirm.

Metrics: **(N) nearest-zero** = smallest `|arg|` (largest `𝒜`); **(T) atom-tracking** = arg closest to `π`.

## §1 · Results (witness `w7/wi_embedded.py`)

| `ε` (continuum wt) | 𝒜_N (nearest) | spread_N | 𝒜_T (track, f²=1) | spread_T | pure-cut 𝒜_N |
|---|---|---|---|---|---|
| 0 | 1.000 | 0 | 1.000 | 0 | — |
| 1e-4 | **14.56** | 0.41 | 1.000 | 0.061 | 15.13 |
| 1e-2 | 14.00 | 0.26 | 1.000 | 0.063 | 14.39 |
| 0.1 | 13.39 | 0.15 | 1.000 | 0.063 | 13.61 |
| 1.0 | 12.26 | 0.59 | 1.028 | 0.063 | 12.18 |

**(N)** jumps to `≈12–15` for **any** `ε>0` (a spurious near-real continuum zero becomes the nearest), and
`atom+continuum ≈ pure-cut` ⇒ **atom invisible, `ε*=0`** (true for shape rise/flat, threshold Eth=2/3).

**(T)** keeps a root near `E₁=1.000` (spread over `f²` ≈ 0.06). But a **pure cut** also has a stable
arg-near-`π` zero at **`𝒜_T≈1.065`** (a zero of the shape polynomial). So atom (`1.000`) vs cut (`1.065`):
**separation ≈ 6.5% ≈ the atom's own 6% wobble** ⇒ marginal, shape-dependent, not robust.

## §2 · Resolution and physical plausibility

- **WI-ADMISSIBLE is correct for the observable that can actually be extracted**: the nearest-zero 𝒜 (what
  ESPRIT gives, WI-REALTIME) is destroyed by any continuum — a 1% level already gives `𝒜≈3.16` (their
  number), a band gives `≈14`. Both are cut-like.
- **WI-LATTICE's "isolation suffices" is true only for a pure single level** (zero spread). The physical
  atom carries continuum weight above it, so its nearest-zero signature is cut-like.
- **The tracking workaround is weak**: 6.5% atom-vs-cut separation, comparable to the wobble, and requires
  resolving the atom root against near-real continuum zeros — exponentially ill-conditioned from `O(10)`
  real-τ points (WI-REALTIME's Hankel blow-up). At the crossover the two-glueball continuum is `O(1)`-coupled
  (`ε` not `≪1`), so this is exactly the hard regime.

**Conclusion:** the physical atom+continuum is not cleanly certifiable by WI. Consistent with WI-ADMISSIBLE,
WI-REALTIME, CONSTRUCT-IR-ANALYTIC-NUM. `[T]` pure-atom invariance; `[NG]` nearest-zero for physical case;
`[V-num]` the 6.5% tracking separation; `[O]` uniform extraction at `β≈2.5` = H1.

*Files:* `w7/wi_embedded.py` (runs clean, `py -3.12 -I`). Not committed.
