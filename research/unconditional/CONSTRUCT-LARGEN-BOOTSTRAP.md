# CONSTRUCT-LARGEN-BOOTSTRAP — the large-N loop-space bootstrap does not close: the closure deficit (non-rectangular loop topologies vs rectangles) grows exponentially, so no finite truncation pins σ

> **VERDICT (no QED).** The Makeenko–Migdal loop equation is **closed in `W`** at `N=∞` (a product of the
> same functional, no `n`-point tower — a real structural win over the DSE). But a finite truncation to
> simple/rectangular loops is **not closed**: in `D=4` the loop equation sources **exponentially many
> non-rectangular loop topologies** (polyominoes `~4.06^n` vs rectangles `d(n)=O(n^ε)`), so the closure
> deficit grows without bound and **no finite truncation determines `σ`**; moreover `D=4` forces **two
> functionals** (area + perimeter), so a one-scalar bootstrap is underdetermined by construction. `D=2`
> closes *exactly* (the area law is multiplicative at a self-intersection), which is why 2D is solvable.
> **The bootstrap yields neither a converged `σ` nor `m_{0⁺⁺}/√σ`.** `[NG]` for this construction.
> Clay 0%.

---

## §0 · Setup and labels

Large-N loop equation (`N→∞`, the `1/N²` term drops):
```
(∂/∂x_μ)(δ/δσ_{μν}(x)) W(C)  =  λ ∮ dy_ν δ(C_{xy}) W(C_{xy}) W(C_{yx})        (MM)
```
closed in `W` (RHS = product of the same functional). **Labels:** `[T]` exact · `[V-num]` numeric · `[NG]` no-go for the named class · `[O]` exact missing step. Acid: vanish on free massless *and* `:φ²:`.

---

## §1 · The construction and the three checks

The "bootstrap the master field" program would truncate loop space to simple loops `W(m,n)` (rectangles
`m×n`), impose the lattice MM equation + reflection positivity + large-N factorization, and solve for `W`
and hence `σ`. Witness `w7/construct_largen_bootstrap.py` (numpy + `math`, `<1 s`) tests the three
structural facts that decide whether this can work.

### (1) `D=2` closes exactly `[T]`
The exact 2D master field is `W(C)=exp(−λA/2)`. Because the signed area is additive at a self-intersection
(Stokes), a figure-eight of two plaquettes (area 2) factors *exactly*:
```
W(fig8) = exp(−λ·2/2) = exp(−λ/2)² = W(plaq)²   (|diff| = 0 to machine zero)
```
⇒ the "extra" loop the MM equation generates is **already determined** by the area law; **no deficit**.
This is precisely why 2D large-N YM is a solved master field.

### (2) `D=4` closure deficit grows exponentially `[T]/[V-num]`
Adding a plaquette at a link of a rectangle produces a **connected square-footprint** = a polyomino. The
number of polyomino topologies (OEIS A000105) vastly exceeds the number of rectangles of the same area:

| area `n` | polyominoes | rectangles | deficit | deficit/total |
|---|---|---|---|---|
| 1 | 1 | 1 | 0 | 0.000 |
| 2 | 1 | 1 | 0 | 0.000 |
| 3 | 2 | 1 | 1 | 0.500 |
| 4 | 5 | 2 | 3 | 0.600 |
| 5 | 12 | 1 | 11 | 0.917 |
| 6 | 35 | 2 | 33 | 0.943 |
| 7 | 108 | 1 | 107 | 0.991 |
| 8 | 369 | 2 | 367 | 0.995 |
| 9 | 1285 | 2 | 1283 | 0.998 |

Polyominoes grow `~4.06^n`; rectangles grow `d(n)=O(n^ε)`. **⇒ the number of loop topologies the MM
equation *needs* but a rectangle-truncation *omits* grows exponentially.** A truncation to simple loops is
therefore never closed: at every order the omitted unknowns overwhelm the retained equations.

### (3) `D=4` forces two functionals `[V-num]`
A realistic 4D master field carries area **and** perimeter (Coulomb) dependence, `W(m,n)=exp(−σ·mn − c·(m+n))`.
Reading a *pure-area* `σ_eff = −ln W/(mn)` gives a value that **drifts with loop size** (witness: `σ_eff`
runs `1.10 → 0.65` over `1×1 … 4×4`, spread `0.45`) — one scalar cannot absorb the drift. A two-functional
fit `Y = σ·A + c·P` (`A=mn, P=2(m+n)`) recovers the coefficients **exactly** (`σ=0.500`, `c=0.150` in the
`2(m+n)` convention, i.e. `c=0.300` in the `(m+n)` convention; residual `1.8×10⁻¹⁵`, rank 2). ⇒ **two
independent functionals are forced** (the "surface sum" of `PNT-LARGEN-LOOP`), so a one-unknown-`σ`
bootstrap is **underdetermined by construction**.

---

## §2 · Does the bootstrap converge? Does it yield `m_{0⁺⁺}/√σ`?

- **Converge:** No. The closure deficit grows exponentially (fact 2) and the closure is underdetermined
  in any one-scalar scheme (fact 3); only `D=2` closes (fact 1). A finite rectangle-truncation is a
  *truncation with uncontrolled error*, not a convergent bootstrap.
- **Yield `m_{0⁺⁺}/√σ`:** No. `σ` is not pinned. Worse, the `0⁺⁺` mass is the lightest eigenvalue of the
  **linearized fluctuation around the master field** `W` (`PNT-LARGEN-LOOP`); a bootstrap that does not
  produce `W` (only `σ`-level data at best) cannot linearize, so it gives **no** `m_{0⁺⁺}` even if `σ` were
  pinned. And the loop equation lives in the **1-form channel** (`σ`, confinement), while `m_{0⁺⁺}` is
  `0`-form — the channel wall.

---

## §3 · What is genuinely true (the structural wins)

- The MM equation **is** closed in `W` with **no `n`-point tower** — a genuine advantage over the DSE
  (whose 1-contour equation sources 2, 3, …). This is why the large-N route is the natural "closed" object.
- `D=2` is **exactly** a closed master field (`W=exp(−λA/2)`) — the only dimension where the truncation
  closes, and the calibration for facts (2)–(3).
- The obstruction is **structural and combinatorial** (polyomino growth) + **kinematic** (two functionals),
  not a tuning problem: it does not shrink with a better truncation scheme.

---

## §4 · Honest bottom line

- **`[NG]` for the large-N loop-space bootstrap** as a `σ`/`m_{0⁺⁺}` machine. The deficit grows
  exponentially; the one-scalar closure is underdetermined; the 0⁺⁺ mass needs a linearization the
  bootstrap does not supply; and the loop equation sits in the 1-form channel.
- **Exact missing step `[O]`:** a *closed* 4D master field `W(C)` — equivalently a summable surface
  functional pinning both functionals — from which one can linearize and read the `0⁺⁺` eigenvalue. That is
  the same object as `ROOF-LARGEN-MASTER`'s "construct the 4D master field", here shown to be *not*
  reachable by loop-space truncation + positivity alone.
- **Genuine residue:** the closed-in-`W` property and the 2D calibration; they confirm the loop equation is
  the right *kind* of closed object but sit on the far side of the 0-form/1-form channel wall.
- **Clay 0%.** No QED.

---

*Files:* `research/unconditional/CONSTRUCT-LARGEN-BOOTSTRAP.md` (this file), `research/unconditional/w7/construct_largen_bootstrap.py`.
*Refs:* `PNT-LARGEN-LOOP.md`, `PNT-LARGEN-DUAL-COMBO.md`, `ROOF-LARGEN-MASTER.md`, `ROOF-CONFINEMENT-GAP.md`.
