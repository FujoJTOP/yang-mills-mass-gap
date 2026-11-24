# FALSIFY-TOPOLOGICAL — topological/homotopy invariants do NOT falsify "only exit = YM structure or H1": the winding invariant is a genuine atom-discriminator only as a weak-* DISCONTINUOUS functional, hence not computable from finite data ⇒ it reduces to A3/infinite

> **FALSIFICATION ATTEMPT: FAILED, for a precise reason.** I tried to force an atom from finite data via a
> **topological/homotopy invariant** — the winding number of the Herglotz function `Π(z)=∫dν/(s−z)` on a
> contour enclosing the floor (an argument-principle / degree-of-map count of poles, stronger than the
> zero-count that `TOOL-ZEROCOUNT` refuted). It does **not** falsify the claim. The winding **is** a genuine
> atom-vs-continuum discriminator (atom → closed loop winding `−1`; continuum → a branch cut, no closed loop),
> **but it is a weak-* DISCONTINUOUS functional**: the atom is exactly a `δ` (weak-* discontinuity) in the
> boundary density, while every finite datum is weak-* continuous. So the invariant is **not computable from
> finite data** — to evaluate it you must go to the real axis / boundary values = **A3 / infinite data**. The
> topological route therefore reduces to the already-named exit, and **the only-exit claim STANDS.**

---

## §0 · Labels

**[T]** 真定理 · **[NG]** no-go · **[O]** 精确缺失步 · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The candidate and the witness (`w7/falsify_topological.py`, runs clean)

`Π(z)=∫dν(s)/(s−z)` is a Herglotz (positive-imaginary) function. Atom = isolated pole at `s₁` (residue `−f`); continuum = a **branch cut** on the real axis.

**(1) Weak-* density `[T]/[V-num]`.** An **atomless** narrow bump (total mass `=f`, width `h=0.03`) matches the atom at off-axis points to `|Π_atom−Π_bump| = 3.6×10⁻³, 7.2×10⁻⁴, 3.6×10⁻⁴, 1.9×10⁻⁴`. So **any finite set of point-evaluations is fooled** (this is OPE-RHO in the Herglotz language).

**(2) Off-axis contour is BLIND `[T]/[V-num]`.** On a circle with `Im z ≥ 0.1` (atom `s₁∈ℝ` outside), `winding(Π_atom) = 0.000000` and `winding(Π_bump) = 0.000000` — **the topological invariant is weak-* continuous off-axis and identical.**

**(3) The atomless measure is feasible `[T]`.** `dν=ρ(s)ds`, `ρ≥0` on `[s₁,∞)`, has **no atom at the floor** yet is a valid `ν`: isolation gap can be `0`. So **no positivity/topological count can exclude it** — the same skeleton `TOOL-ZEROCOUNT` hit.

**(4) To see the atom you must cross the real axis `[T]/[V-num]`.** A contour that **encloses** the atom (at `s₁∈ℝ`) must cross ℝ. There `winding(Π_atom) = −1.0000` (pole enclosed), but `Π_bump` has a **branch cut**: its jump across the axis at `x=s₁+0.05` is `|Π(+i0)−Π(−i0)| = 34.2` (finite, `=2πρ`) — a cut, not a closed loop. So the discriminator lives **on the real axis** (the boundary density: a `δ` for the atom, smooth `ρ` for the continuum), which is **A3 / infinite data**.

---

## §2 · Why it fails (the precise reason)

The atom is **exactly** the weak-* discontinuity of `Π`'s boundary value (a `δ` in the density `ρ`). Every
**finite** datum (moments, point values, Taylor coefficients) is a **weak-*-continuous** integral, and atomless
measures are weak-* **dense** (`[T]`). Hence **any functional that distinguishes the atom from a matching bump
must be weak-* discontinuous.** The winding number *is* such a functional (it jumps `0→−1`), which is why it
discriminates — **and exactly why it is not finitely computable**: a finite sample of `Π` is matched by an
atomless bump, so the winding cannot be extracted from finite data. Evaluating it needs the boundary value
across the real axis = analytic continuation = **A3 / infinite data**.

This is the same obstruction as `TOOL-ZEROCOUNT` (positivity admits `s₂=s₁`), now made precise via weak-*
density: **homotopy/degree invariants cannot escape it, because they are weak-* discontinuous and the finite
data is not.**

---

## §3 · Verdict and yield

- **`[NG]`** topological/homotopy **forcing** of the atom: fails; reduces to A3 / infinite data.
- **`[T]`** the sharp characterization: **the atom is exactly the weak-* discontinuity (a `δ`) of the boundary
  density; every atom-forcing invariant must be weak-* discontinuous; topological invariants are, hence are
  non-finitely-computable = A3/infinite.**
- **The only-exit claim STANDS** for the topological route (confirmed, not falsified). The only routes that
  could still force are the two already named: a **nonlinear YM-structure constraint** (not "`ν` a positive
  measure") or **infinite data / uniform isolation `J→∞` = H1**.

*(Genuine small yield: the exact weak-* statement.) No QED; Clay 0%.*
