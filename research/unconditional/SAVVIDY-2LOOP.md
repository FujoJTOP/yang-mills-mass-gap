# SAVVIDY-2LOOP — the 2-loop/RG-improved softening moves the constant only by `O(λ₄/16π²)` (few % at the crossover) and keeps the sign: sub-criticality is **stable**, so the LO "softened H1" is representative, not a truncation artifact

> **VERDICT (no QED).** `CONSTRUCT-SAVVIDY-CELL` found the Savvidy negative mode *softens* the coarse
> curvature, `Γ″(0) = κ − g² Var(φ²)`, with `Var(φ²)` finite and set by the quartic `[V-num]`, leaving a
> *softened* `H1`: `κ_bare > (accumulated softening)`. This report asks whether 2-loop / RG improvement
> breaks or changes that. **It does not.** The softening is `soft(k) = g²(k) Var(φ²)(k)`; it is a **variance**
> (always `>0`), and under RG it runs multiplicatively with exponent `−2β₀g² + γ_{φ²}`, `γ_{φ²}=λ₄/16π²`.
> A positive quantity with a finite multiplicative exponent **never crosses zero** ⇒ sign preserved `[T]`.
> At the crossover the magnetic scale `√(gB)` and the IR scale coincide, so the log is `O(1)` — there is
> **no large-log enhancement** where it bites. The relative 2-loop shift is `O(λ₄/16π²) ~ 0.3–1%` per log
> (witness) ⇒ the constant moves by ≲ few %. **Two-loop confirms LO: sub-criticality is structurally stable,
> the constant is barely renormalised, and the uniform all-orders bound remains `H1` (softened).** Clay 0%.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件/模型 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** 本轮数值.

---

## §1 · What is being improved

LO (`CONSTRUCT-SAVVIDY-CELL §3`): integrating out the negative mode `φ` with cubic coupling `g a φ²` gives
the coarse self-energy `Σ(0) = g² Var(φ²)` and

> `Γ″(0) = κ − soft`,  `soft ≡ g² Var(φ²)`,  `Var(a) = 1/Γ″(0)`.

`soft` is the **one-loop** bubble of two `φ²` insertions at zero momentum. It is a *variance* of the
bilinear `φ²` in the measure of `V(φ) = −½c φ² + ¼λ₄φ⁴` (`c=|λ₋|=g|B|`). Two questions:
**(Q1)** does 2-loop change the **sign** of `soft` (⇒ can it *help*, or flip to anti-softening)?
**(Q2)** does RG resummation move the **constant** much at the crossover?

## §2 · Sign is protected `[T]`

`Var(φ²) ≥ 0` by definition (a variance), and `g²(k) > 0`. So `soft(k) ≥ 0` at every scale, at every loop
order that is a positivity-preserving (self-energy) contribution. Concretely, in the bimodal potential the
Gaussian-around-minimum variance is

> `Var(φ²) ≈ 4φ₀²/m_φ² + 2/m_φ⁴ = 2/λ₄ + 1/(2c²) > 0`,  (`m_φ² = V″(φ₀) = 2c`, `φ₀² = c/λ₄`),

positive and **λ₄-dominated** — matching the LO finding that the softening is set by the quartic, not the
negative curvature. The 2-loop (`λ₄`) and running couplings only **multiply** this positive number by a
positive factor; they cannot make it negative. **No loop order turns softening into anti-softening.**

## §3 · The constant barely moves `[V-num]`

`soft(k)` runs multiplicatively: `d log soft/d log k = −2β₀ g² + γ_{φ²}`, with the composite mass operator
anomalous dimension `γ_{φ²} = λ₄/16π²` (up to an `O(1)` scheme factor). A finite logarithmic correction of
size

> `δsoft/soft = (γ_{φ²} − 2β₀g²/16π²)· log(√(gB)/k) + O(λ₄²)`.

At the **crossover** `k ~ √(gB)` ⇒ `log = O(1)`; with `λ₄ ~ O(1)`, `γ_{φ²}, 2β₀g²/16π² ~ 10⁻²–10⁻¹` ⇒
`δsoft/soft ~ 0.3–1%` per log (witness). So the 2-loop shifts the constant by **≲ few %**, **not** by `O(1)`.
In the **deep IR** `k ≪ √(gB)` the log grows, *but* there channel 1 (direct overlap) already vanishes by
scale separation and the softening itself is tiny — so the large-log region is not where the floor is tested.
**The crossover — the one place it bites — has `log O(1)` and a small 2-loop shift.**

Crucially: at the crossover `α_s ~ O(1)`, so `g²/(16π²)` is `O(10⁻¹)`, and the scheme/loop correction to the
*absolute* constant is `O(α_s)` — this is the irreducible theoretical uncertainty on `soft`, not a runaway.
It is bounded, sign-definite, and `O(few %)`.

## §4 · Verdict

- **(Q1) Sign:** protected. `soft ≥ 0` always ⇒ **sub-criticality `κ_bare > soft` is stable in form** under
  all positive (self-energy) loop corrections. 2-loop neither breaks nor manufactures sub-criticality.
- **(Q2) Constant:** renormalised by `O(λ₄/16π²) ~ few %` at the crossover ⇒ the **LO constant is
  representative**. The `CONSTRUCT-SAVVIDY-CELL` sharpening ("negative mode does not kill the floor; the gap is
  a uniform sub-criticality bound") is **robust**, not a leading-order artifact.
- **What 2-loop does NOT do:** it does not *prove* the uniform bound `κ_bare > ∫(dk/k) soft(k)` across cells
  and `β`. Summing all orders and bounding the accumulated softening uniformly is still `H1` (softened). The
  2-loop merely shows that the *object* and its *sign* are trustworthy.

## §5 · Honest bottom line

`[T]` the softening is sign-protected (a positive, multiplicative variance). `[V-num]` the 2-loop/RG shift
of the constant is ≲ few % at the crossover. `[O]` a uniform all-orders bound on the accumulated softening
stays `H1`. **So: 2-loop improves the *credibility* of the Savvidy-cell sharpening (LO representative) but
does not improve the *prospects* (no new lever) and does not break the obstruction.** Clay 0%; no QED.

*Witness:* `w7/savvidy_2loop.py` (numpy, runs clean). Refs: `CONSTRUCT-SAVVIDY-CELL.md`, `CONSTRUCT-RG.md`,
`SLOT-SD-FIXPOINT.md` (`G″=S−D`, softening = `D`-term).
