# FALSIFY-FLOW — spectral flow of a deformed operator does NOT force the atom; the "only-exit" claim stands (spectral flow is a spectrum/gap object, blind to the `tr F²` overlap)

> **VERDICT (no QED).** The attempt to falsify the only-exit claim ("forcing the atom needs YM-structure or
> `J→∞`=H1") via **spectral flow of a deformed `H(λ)`** **fails** — for two independent reasons, both
> `[T]`. So the only-exit claim is **NOT falsified**; spectral flow is a third *reading / index-type* object,
> blind to the atom. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.

---

## §1 · The route and why it fails

**Route.** Deform the transfer operator `H(λ)` by a one-parameter auxiliary (mass / twist); the **spectral
flow** = net number of eigenvalues crossing a reference as `λ` varies, a **topological integer**. Claim to
test: a nonzero flow FORCES the `tr F²` atom `f₁² = |⟨0|tr F²|G⟩|² > 0` at the floor, without YM-structure=H1
or `J→∞`.

**Obstruction 1 `[T]` — spectral vs overlap.** The flow is a property of the **spectrum** (`H(λ)`); the atom
`f₁²` is a property of the **operator overlap** `⟨0|O|G⟩`. These are **independent data**: the same `H(λ)` with
the same flow supports different `f₁²`. **Witness** (`w7/falsify_flow.py`, exact): with `H(λ)=diag(λ−1,λ+1)`
the flow is `1` for *every* operator `O`, yet
- `O = σ_x` ⇒ `f₁² = 1.0` (atom present);
- `O = 𝟙` ⇒ `f₁² = 0.0` (NO atom, **yet flow = 1**);
- `O = σ_z` ⇒ `f₁² = 0.0` (NO atom, **yet flow = 1**).

⇒ **a nonzero spectral flow does NOT force `f₁² > 0`** — the crossing eigenvalue can be `tr F²`-orthogonal.
Spectral flow is a "gap detector" (crossing ⇔ gapless), not an atom certifier — the same gap-vs-atom
dichotomy as `BASE-INDEX` (I2/I3): a flow certifies `gap=0` (a crossing exists), never `f₁²>0`.

**Obstruction 2 `[T]` — rigidity.** Spectral flow is an **integer, constant on connected components** of the
deformation parameter `λ` (it can change only at degeneracies). A `β`-dependent atom floor is therefore
**not** constrainable by a flow that is topological across the physical regime — the same rigidity no-go that
killed `BASE-INDEX` (I1). Witness: a `β`-deformation `H_s=diag(1+s,2+s)` with no crossing has flow `0`,
constant over `s` — it carries no `s`-dependent atom information.

---

## §2 · What spectral flow *is* (and is not)

- **Is:** a genuine global, non-integral-in-`ν` object (eigenvalue crossings / a topological count of states).
  It does detect the *gap closing* (a crossing).
- **Is not:** a determinant of the matrix element `⟨0|tr F²|G⟩`. The atom weight is a *residue / overlap*, not
  a spectral count. A flow cannot manufacture an overlap.

---

## §3 · Exact missing step `[O]`

The atom needs the **matrix element** `⟨0|tr F²|G⟩`. The spectrum `H(λ)` and its flow do **not** determine
it — the matrix element is precisely the "YM-structure" (which `O=tr F²` couples to which states) or the
A3/H1 residue. So forcing the atom remains where the campaign put it: **nonlinear YM-structure, or
`J→∞` = H1.** Spectral flow is a third reading-type coordinate, not an escape.

---

## §4 · Honest bottom line

- **The only-exit claim is NOT falsified.** Spectral flow of a deformed operator is a **reading / gap-detector**
  object (topological, rigid), blind to the `tr F²` overlap; a nonzero flow coexists with `f₁²=0` (witness).
- **What it adds:** one more clean illustration of the campaign's central dichotomy (spectrum ≠ overlap;
  gap ≠ atom), and a concrete `[T]` counterexample to "flow forces the atom."
- **Clay = 0%.** The atom still needs YM-structure or `J→∞` = H1.

*Files:* `research/unconditional/FALSIFY-FLOW.md`, `research/unconditional/w7/falsify_flow.py`.
*Refs:* `BASE-INDEX.md` (rigidity, gap-detector), `TOOL-SYNTHESIS.md` (values=integral, zeros/poles=A3),
`lean/YMStatement/H1Sharpening.lean` (H1 ⊊ gap).
