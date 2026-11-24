# IMPOSSIBLE-ANALYTIC — analytic continuation from a solvable YM limit: DEAD at an ESSENTIAL SINGULARITY (dimensional transmutation); the continuation *is* the wall, but the framing (resurgence-type non-analyticity) is a genuinely distinct, sharper obstruction

> **DIRECTIVE.** Jump out of the fixed bases (KL measure / spectral atom / χ / weak-* / H1). Impossible base: compute the gap **exactly** in a solvable YM limit, then **transfer to ℝ⁴ by an analyticity/rigidity theorem**.
>
> **VERDICT (no QED, Clay 0%).** The base **fails**, and it fails for a sharp, *new* reason — **not** the weak-* measure argument: **the gap is an essential singularity in the coupling (`m = μ·e^{−1/(2β₀g²)}·F(g²)`), so every analytic-in-`g²` transfer from a solvable point is provably blind.** At weak coupling the gap is `0` to **all orders** (all Taylor coefficients of `e^{−1/g²}` vanish); at the strong-coupling lattice end it is the **wrong (infinite-`a`) limit**; on the torus it needs a **uniform lower bound** (= the transport = H1). The continuation is the wall, but it is a **resurgence/non-analyticity** statement, distinct from the earlier measure-theoretic ones.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The three obstructions (each a [T]/[O])

**(1) Essential singularity = dimensional transmutation `[T]`.** In pure YM the only scale is `Λ = μ·e^{−1/(2β₀g²)}·…` (one-loop); every mass `m = Λ·F(g²)`. `e^{−1/g²}` has **all derivatives zero at `g=0`** (witness: `e^{−1/g²}/g^{2k}→0` for every `k`; numerically `m/μ=0` for all `g²≲10⁻²`). ⇒ the **perturbative (analytic-in-`g²`) expansion of the gap is identically 0**. **Analytic continuation in `g²` from weak coupling gives `m=0`, never the physical `g²~1` value.** The gap is *maximally* nonperturbative (zero Taylor series), so no analytic transfer reaches it — a resurgence wall (`ROOF-RESURGENCE`: the `u=2` transseries coefficient is free).

**(2) Strong-coupling is the wrong limit `[T]/[V-num]`.** The exact lattice strong-coupling gap `m·a = 4 ln(1/u(β))` **grows** as `β→0` (`2.7→14.8` for `β=2.5→0.1`), but `a→∞`: the physical `m = m·a/a → 0`. The strong-coupling expansion gives a gap in *lattice units* that vanishes in the continuum. Wrong limit.

**(3) Torus `T³×S¹`: continuity ≠ lower bound `[T]/[O]`.** `m(L)` is continuous in `L` (single phase, no transition), so `m(∞)=lim_{L→∞}m(L)`. But continuity alone gives nothing if `m(L)` may dip to `0`; a small-`L` "solvable" value is the **3D effective-theory gap** (itself unproven), and the crossover is uncontrolled. Extending to `L→∞` needs a **uniform lower bound** over all `L` — **the transport = H1**.

---

## §2 · Why no solvable limit is both exact and physical

The two **exact** ends of YM — weak coupling / UV (`g²→0`, where `m→0` by (1)) and strong coupling / IR (`β→0`, an infinite-`a` lattice by (2)) — **never meet at the physical scale**. The gap exists only at the **crossover** (`g²~1`, `α_s~O(1)`), controllable by **neither** (cf. `ROUTE-NEW-RG`: "the two rigorous ends do not meet"). The solvable limits are either **different theories** (2D: kinematically degenerate, `D−2=0`; SUSY: not pure YM; Banks–Zaks: gapless CFT; large-`N`: not exact) or **gapless/wrong-scale** ends. There is no exact gap at the physical point to continue *from*.

---

## §3 · The new content (what this base revealed that the others did not)

The earlier obstructions were **measure-theoretic** (atom = weak-* discontinuity; finite data weak-* continuous). This base gives a **distinct, sharper** statement:

> **The gap is an essential singularity in `g²`. Any procedure analytic in the coupling — perturbation theory, Borel/analytic continuation, `1/N` truncation, small-torus expansion — is *provably* blind to it, because `e^{−1/g²}` has zero Taylor series. Analytic continuation cannot manufacture a nonperturbative scale.**

This is the **resurgence coordinate** of the wall, and it is stronger than "non-analytic": it says the gap is **invisible to every analytic method, not merely hard**. It also explains *why* the solvable-limit programs (Eguchi–Kawai, small torus, deformed) all "reduce to the physical crossover": the deformation path crosses the essential singularity, and the exact value on the far side is `Λ`, not analytic data.

---

## §4 · Honest bottom line

- **The impossible base fails**, but it fails at a **new, sharper obstruction**: the gap is an **essential singularity** in `g²` (dimensional transmutation), so **no analytic-in-coupling continuation from any solvable limit reaches it**. The continuation *is* the wall.
- **Not a different exit** — it reduces (honestly) to H1 in the **resurgence/non-analyticity** coordinate, which is **distinct from** and **stronger than** the weak-* measure coordinate: it rules out the entire class of *analytic transfer* methods at once, by the structure of `e^{−1/g²}`.
- **Yield:** a clean statement — *YM's gap is non-perturbatively invisible (zero Taylor series in `g²`); any "compute exactly in a solvable limit, then continue" program is dead by construction.* The gap must be obtained **non-analytically** (a non-perturbative construction), which is `H1`.
- **Clay = 0%.** No QED.

*Witness:* `w7/impossible_analytic.py` (essential singularity `e^{−1/g²}/g^{2k}→0` ∀k; strong-coupling `m·a→` grows while `a→∞`; torus continuity-without-lower-bound).
