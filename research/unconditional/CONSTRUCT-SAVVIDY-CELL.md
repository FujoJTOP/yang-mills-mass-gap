# CONSTRUCT-SAVVIDY-CELL — the Savvidy negative mode **softens** the χ floor, it does not destroy it: the direct channel vanishes by scale separation, the induced softening is bounded by the quartic ⇒ the floor becomes a **softened** `H1`

> **VERDICT (no QED).** `CONSTRUCT-RG` pinned its exact missing step to a sign bound
> `Var_μ(a_{k→0}) ≤ 1/λ` on the typical-field cell where the Savvidy Hessian is **indefinite**
> (`λ_min = −g|B|`), and the BBS large-field convexity fails. This report attacks that cell directly
> and constructs the best available estimate of `Var(a_{k→0})` in it. **The negative mode does NOT
> categorically destroy the floor.** It reaches the `k→0` variance through exactly two channels:
> **(1) direct** — the negative mode is a **lowest-Landau-level helicity state** with Fourier support
> `|p| ~ √(gB)`, orthogonal to the `k→0` mode, so its direct contribution is `O((k/√(gB))²) → 0`
> `[T]`; **(2) RG-induced softening** — the cubic coupling `g a φ²` softens the coarse curvature to
> `Γ″(0) = κ − g² Var(φ²)`, with `Var(φ²)` **bounded by the quartic** (not by the negative curvature)
> `[V-num]`. The floor survives iff `κ_bare > g² Var(φ²)` — a **finite threshold shift**, not a blow-up.
> What remains open is the **uniform (across cells/β) bound on the accumulated softening**, which is
> `H1` (softened). **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件/模型 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** 本轮数值.

---

## §1 · The cell and the two channels

In the Savvidy cell the quadratic (fluctuation) operator around the coarse field is
`Hess = H₊ − |λ₋| P₋`, with `λ₋ = −g|B|` a **rank-one** negative eigenvalue. The naive
Brascamp–Lieb / Gaussian bound `Var(a) ≤ 1/(positive curvature)` fails because one curvature is
negative. We separate two ways the negative mode can reach `Var(a_{k→0})`.

## §2 · Channel 1 — direct: vanishes by scale separation `[T]`

The Savvidy negative mode is the lowest-Landau-level gluon helicity state: its fluctuation eigenvalue
is `(2n+1)gB + 2gB·s` at `n=0, s=−1`, i.e. `−gB`, and its wavefunction is a lowest-Landau-level mode
localized on the magnetic length `1/√(gB)`. **Its Fourier support is `|p| ~ √(gB)`, not `p=0`.**

The `k→0` coarse mode has support `|p| < k`. For `k < √(gB)` the supports are **disjoint**, so the
Fourier overlap is **identically zero**; with a smeared cutoff it is `O((k/√(gB))²)`. Hence the negative
mode's **own amplitude variance does not contribute to `Var(a_{k→0})`** — it contributes `→ 0` as
`k → 0`. (Witness W1: overlap² `= 0.0025, 0.01, 0.0625, 0.25, 1` at `k = 0.05…2`.)

This is the analytic face of the corpus finding (`W6-BE`) that **90–97% of the negative weight is
high-momentum**: the negative mode lives at the magnetic scale, orthogonal to the IR variance.

## §3 · Channel 2 — RG-induced softening: bounded `[V-num]/[C]`

The negative mode also reaches the coarse mode **indirectly**, through the cubic vertex `g a φ²`. Integrating
out the negative mode (toy reduction, exact quadrature) gives

> `Γ″(0) = κ − g² Var(φ²)`,  `Var(a) = 1/Γ″(0)`,

where `φ` is the negative-mode amplitude under the effective potential
`V(φ) = −½ c φ² + ¼ λ₄ φ⁴` (`c = |λ₋| > 0`). The witness computes `Var(φ²)`:

| `c=|λ₋|` | `λ₄` | `Var(φ²)` | `g²Var` (g=1) |
|---|---|---|---|
| 1.0 | 1.0 | 0.956 | 0.956 |
| 1.0 | 0.5 | 2.381 | 2.381 |
| 1.0 | 2.0 | 0.406 | 0.406 |
| 2.0 | 1.0 | 1.557 | 1.557 |
| 4.0 | 1.0 | 2.209 | 2.209 |
| 1.0 | 0.2 | 8.513 | 8.513 |

**`Var(φ²)` is finite and set by the quartic `λ₄`, NOT divergent as `c` grows.** The softening is a
**finite threshold shift**: `Var(a)` stays finite while `κ > g²Var(φ²)`, and blows up only when `κ` drops
**below** the softening (witness: `Var(a) = 0.52, 2.09, 20.9, ∞, ∞` as `κ/soft = 3, 1.5, 1.05, 1.0, 0.9`).
So the negative mode **shifts the threshold, it does not open an unbounded hole**.

## §4 · Verdict — controlled, not destructive; a *softened* `H1`

- **IR (`k ≪ √(gB)`):** channel 1 protected (overlap `→0`), channel 2 a bounded softening ⇒ the χ floor
  **survives**. The Savvidy cell is **harmless in the deep IR**.
- **Crossover (`k ~ √(gB)`):** both channels are `O(1)` ⇒ a genuine `O(1)` competition. **This is exactly
  where the Savvidy cell can bite** — and it is the crossover, not the IR, that is hard.
- The floor thus becomes a **softened inequality**

> `κ_bare  >  (accumulated Savvidy softening over scales)   ⟺   H1 (softened)`,

  i.e. the bare coarse curvature must exceed `∫ (dk/k) g²(k) Var(φ²)|_k · (RG propagation)`. At the
  crossover `g²/λ₄ = O(1)`, so this is a genuine `O(1)` competition — not a divergence.

This is consistent with the corpus variance identity `G″ = S − D` (`SLOT-SD-FIXPOINT`): the Savvidy
softening `−g²Var(φ²)` is the `D`-term (the Brascamp–Lieb defect from the non-log-concave fiber), and
the statement "`D`-term < `S`-term" is exactly `H1`.

## §5 · Honest bottom line

- **The negative mode does not destroy the χ floor** — its direct channel vanishes by scale separation
  `[T]`, and its induced softening is bounded by the quartic `[V-num]`.
- **What remains open `[O]`:** a **uniform** bound on the accumulated softening
  `∫(dk/k) g²Var(φ²)` across cells and `β`, staying below `κ_bare`. This is `H1` in the Savvidy-cell
  coordinate — sharpened from "does the negative mode kill the floor?" (answer: no) to "is the accumulated
  softening uniformly below the bare curvature?" (answer: open, `H1`).
- **Clay 0%; no QED.** This is a constructive advance on the RG cell, not a proof.

*Witness:* `w7/construct_savvidy_cell.py` (numpy, runs clean). *Refs:* `CONSTRUCT-RG.md` (the missing step),
`SLOT-SD-FIXPOINT.md` (`G″=S−D`), Savvidy/Nielsen–Olesen instability, `W6-BE` (high-momentum negative weight).
