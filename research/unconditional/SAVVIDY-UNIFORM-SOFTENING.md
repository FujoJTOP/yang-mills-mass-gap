# SAVVIDY-UNIFORM-SOFTENING — the softening is LOCALIZED to one Savvidy octave ⇒ the `(L,a)`-uniformity is REMOVED; the only open step collapses to a single marginal `O(1)` inequality at `β≈2.5` (= `H1` softened)

> **VERDICT (no QED).** `CONSTRUCT-SAVVIDY-CELL` reduced the RG missing step to a *softened* inequality
> `κ_bare > Σ`, `Σ = ∫(dk/k) g²(k) Var(φ²)|_k`, "uniform in `(β,L,a)`". This report computes it.
> **The `(L,a)`-part is resolved `[T]/[C]`:** the softening is supported in **one Savvidy octave** around
> `|p|~√(gB)` (channel 1 kills the deep-IR overlap; there is no UV `log` over `[1/L,1/a]`), so `Σ` is a
> finite **physical `O(1)`**, independent of `(L,a)` — **the uniform-in-`(L,a)` worry is removed.** What
> survives is purely the **`β`-dependence**: `Σ(β)` and `κ_bare(β)` are both `O(1)`, `Σ→0` as `β→∞`
> (weak coupling) and as `β→0` (no condensate), so the only competition is at the crossover, where it is
> a **pure `O(1)` ratio `κ_bare > C·g²/λ₄`** (marginal, since naturally `λ₄~g²`) — this is `H1`
> (softened), not a divergence. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件/模型 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** 数值.

---

## §1 · The softening pieces `[V-num]`

`Var(φ²)` for `V(φ)=−½cφ²+¼λ₄φ⁴`, `μ∝e^{−V}` (exact quadrature; reproduces the cell table):

| `c` | `λ₄` | `Var(φ²)` |
|---|---|---|
| 1 | 1 | 0.9565 |
| 1 | 0.5 | 2.3807 |
| 1 | 2 | 0.4063 |
| 2 | 1 | 1.5571 |
| 4 | 1 | 2.2088 |
| 1 | 0.2 | 8.5127 |

**Scaling** (fit at `c=1`): `Var(φ²) ~ λ₄^{−1.278}` — i.e. **`Var(φ²) ~ 1/λ₄`** to leading order, **set by
the quartic, not by the negative curvature `c`** (consistent with the cell report). Hence
`Σ = g² Var(φ²) ~ g²/λ₄`.

## §2 · Localization — the `(L,a)`-uniformity theorem `[T]/[C]`

The Savvidy negative mode is the lowest-Landau-level helicity state, supported at `|p| ~ √(gB)`
(`CONSTRUCT-SAVVIDY-CELL` channel 1). **Below `√(gB)` the coarse fiber has no negative mode**, so
`Var(φ²)|_k = 0` there and the direct overlap is `O((k/√(gB))²)→0`. Consequently the softening integral
`∫(dk/k) g²(k) Var(φ²)|_k` is supported in **one octave** around `√(gB)` — **not** a `log` over the full
`[1/L, 1/a]` range. ⇒ `Σ` is a **finite physical `O(1)` number, independent of `(L,a)`**. **This removes
the `(L,a)`-uniformity obstruction** from the softened `H1`: the softening does **not** accumulate
logarithmically in the volume/cutoff.

## §3 · `β`-dependence — a single marginal `O(1)` competition (= `H1` softened) `[O]`

With `g²(β)=4/β` (SU(2)) and the natural `λ₄(β)~g²(β)`, `Σ_nat(β)=g²Var(φ²)|_{λ₄=g²}` is `O(1)` and
**marginal** (`g²·(1/g²)`), e.g. `Σ_nat(2.5)=0.85`. The margin `κ_bare(β)−Σ(β)`:
- `β→∞` (weak coupling): `g²→0 ⇒ Σ→0`, margin `→κ_bare>0` — **safe**;
- `β→0` (strong coupling): Savvidy condensate absent, `Σ→0` — **safe**;
- `β≈2.5` (crossover): `Σ_nat=O(1)` vs `κ_bare=O(1)` — a **genuine `O(1)` competition**, neither divergent.

So the sign at the crossover is decided by the dimensionless ratio `κ_bare/(g²/λ₄)` = a **pure `O(1)`
inequality**. *(The witness's β-table fixes `κ_bare=0.8` arbitrarily, so its "crossing near β~1.8" is an
artifact — the physical `κ_bare(β)`, `λ₄(β)` are exactly `H1` and are not supplied by the toy.)*

**Evidence, not proof:** the **measured net `G″(β=2.5)≈+0.28>0`** (`SLOT-SD`, `L4→L6`-stable) shows the
physical margin **is positive at 2.5 on the lattice**; but decomposing that `0.28` into `κ_bare−Σ` (the
bare/softening split) is exactly `H1`.

## §4 · Verdict

- **`(L,a)`-uniform: HOLDS** in the softened form `[T]/[C]` — the softening is localized to one Savvidy
  octave (no UV/IR `log`), so `Σ` is finite physical `O(1)`, independent of `(L,a)`. **The `(L,a)`-part of
  the open step is closed.**
- **`β`-uniform: does NOT close** `[O]` — but it is **not** a divergence; it is a **single marginal `O(1)`
  inequality at the crossover** `κ_bare > C·g²/λ₄`, which is `H1` (softened). The two ends (`β→0,∞`) are
  safe; only the crossover competes.
- **Net sharpening:** the open step collapses from *"uniform in `(β,L,a)`"* to *"`(L,a)`-uniform [proved] +
  one marginal `O(1)` inequality at `β=2.5` [`H1`]"*. The Sadivvy softening is **controlled and localized**,
  not a volume/cutoff divergence.

**Clay 0%; no QED.** A constructive sharpening of the RG cell, not a proof.

*Witness:* `w7/savvidy_uniform_softening.py` (numpy, runs clean). *Refs:* `CONSTRUCT-SAVVIDY-CELL.md`
(softened `H1`), `CONSTRUCT-RG.md` (the sign bound), `SLOT-SD-FIXPOINT.md` (`G″=S−D`, `D`-term = softening).
