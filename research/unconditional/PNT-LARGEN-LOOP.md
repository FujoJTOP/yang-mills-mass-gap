# PNT-LARGEN-LOOP — closing the 4D Makeenko–Migdal loop equation: **structurally better than the DSE tower (closed in one functional `W(C)`, no n-point tower), but not a closed scalar gap equation — the master field is unconstructed, and the loop equation carries σ (1-form), not `m_{0⁺⁺}` (0-form)**

> **VERDICT (no QED).** The large-N Makeenko–Migdal (MM) loop equation is **closed in the single
> functional `W(C)`** — its RHS is a *product* `W(C₁)W(C₂)` of the same functional, so there is **no
> infinite n-point tower** (this is the precise sense in which it is **structurally better than the
> SLOT-DSE tower**). But it does **not** yield a closed scalar gap equation for `m_{0⁺⁺}(λ)`:
> **(i)** in D=4 the spanning surface of a loop is **not** boundary-determined (unlike the D=2 signed
> area), so `W(C)` is a **sum over surfaces** (a string functional), not `exp(−σA)` — two independent
> functionals (Coulomb + area) are forced by the same equation; **(ii)** the D=4 **master field (the
> vacuum `W(C)`) is not constructed** (Singer did D=2 only); **(iii)** the loop equation controls
> **σ** (exponential decay of `W` = 1-form / confinement), while `m_{0⁺⁺}` (0-form) needs the
> **linearized fluctuation spectrum around the unknown master field**. **Structurally better,
> comparable difficulty — and the same 0-form/1-form channel wall.** Clay 0%.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** numeric witness. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · Why the loop equation is *closed* (no tower) `[T]`

The Makeenko–Migdal equation for the Wilson-loop functional `W(C) = ⟨(1/N)Tr P exp ∮_C A⟩` at large N
(`λ = g²N` fixed):

```
∂_μ (δ/δσ_{μν}(x)) W(C)  =  λ ∮_C dy_ν δ^{(4)}(x−y) · W(C_{xy}) W(C_{yx})          (MM)
```

At a self-intersection `x = y`, the loop `C` splits into two sub-loops `C₁`, `C₂`; the RHS is
`W(C₁)·W(C₂)` — a **product of the same functional** `W`. So **(MM) is one closed equation for `W(C)`**.
There is **no sequence of higher n-point functions**: contrast the Dyson–Schwinger tower (`SLOT-DSE-FIXPOINT`),
where the 1-contour equation sources the 2-contour, which sources the 3-contour, … — an infinite
hierarchy. **In this precise, structural sense the large-N loop equation is better than the DSE tower.**
`[T]`

The catch is the *domain*: `W(C)` is a functional on the **infinite-dimensional space of closed curves**,
and (MM) is a nonlinear, non-local functional differential equation there. "No tower" ≠ "solvable."

---

## §2 · D=2: the exact closure, and *why* it works `[T]`

In D=2 the signed area of a closed curve is a **boundary functional** (Stokes):
`A(C) = ½∮_C (x dy − y dx)`. At a self-intersection the signed areas are **additive**, `A(C) = A(C₁)+A(C₂)`,
so `W(C) = exp(−λA(C)/2)` is **multiplicative**, `W(C) = W(C₁)W(C₂)`, and (MM) closes.

Witness (`w7/pnt_largen_loop_witness.py`, exact): figure-eight of two unit squares sharing a vertex —
`A(fig8)=2=A₁+A₂`, `W(fig8)=0.496585=W(sq1)W(sq2)`. This is the mechanism of **Singer's 2D master field**
(and Lévy's rigorous planar result): the *signed area being boundary-determined* is what makes the loop
equation solvable in D=2. `[T]`

---

## §3 · D=4: the obstruction — surface non-uniqueness ⇒ string functional `[NG]`

The D=2 closure uses a fact that **fails in D=4**: in D=4, a 1D boundary does **not** determine the area of
its spanning surface. Witness: the unit-square boundary in the `z=0` plane is spanned by the flat square
(area `1`) **and** by a tent of apex height `h` (area `2√(h²+¼)`: `1.414`, `2.236`, `4.123` at `h=0.5,1,2`).
Same boundary, many surfaces, different areas. `[V-num]`

Consequences `[NG]`:

- **`W(C)` is a sum over spanning surfaces** (a string functional), **not** a function of a single boundary
  functional `A(C)`. So `W ≠ exp(−σA)`; the pure area law is **not** the D=4 MM solution.
- The MM equation forces **two independent functionals**: the short-distance **Coulomb/perimeter** behavior
  `W ∼ exp(−c·L(C))` (coefficient `c`) and the long-distance **area law** `W ∼ exp(−σ·A(C))`
  (string tension `σ`). **`c` and `σ` are physically independent**, both forced. A **single scalar `σ`
  cannot** satisfy the shape-dependent (MM) equation for all loops.
- **No closed scalar gap equation** follows. The exact D=4 solution needs a *shape functional*, and the
  **D=4 master field `W(C)` (the vacuum solution) is not constructed** (Singer: D=2 only).

---

## §4 · The channel wall: the loop equation carries σ (1-form), not `m_{0⁺⁺}` (0-form) `[T]`

Even granting a solution `W(C)`, the loop equation's *large-loop* content is the **exponential decay rate
= the string tension σ** — a statement about the **1-form / center / confinement** sector. The `0⁺⁺` glueball
mass is the **0-form** object (the correlator of `tr F²`, or the lightest **fluctuation** of the master field).
`[T]`

To extract `m_{0⁺⁺}` one must **linearize (MM) around the vacuum solution** and read the eigenvalue
spectrum of the loop-space fluctuation operator (the Migdal–Makeenko "quantum string"). This requires the
**unknown** master field first. So the loop equation gives `σ` directly and `m_{0⁺⁺}` only after solving the
hard part — and it sits on the **1-form side of the same 0-form/1-form channel wall** that `ROOF-CONFINEMENT-GAP`,
`BASE-ORDER-PARAMETER`, and `ROOF-0FORM-DUAL` already mapped. `[T]`

---

## §5 · Structural comparison to SLOT-DSE `[T]/[O]`

| | DSE (`SLOT-DSE`) | large-N loop (MM) |
|---|---|---|
| unknowns | infinite tower of n-point functions | **one** functional `W(C)` |
| closure | needs infinitely many equations | **closed** (RHS = `W·W`) |
| difficulty | infinite tower | **one** infinite-dim nonlinear functional equation |
| what it gives | `m²(0)` schematically | **σ** (1-form) directly |

⇒ **Structurally better** (no tower) `[T]`, **but not easier** `[O]`: the single equation lives on
infinite-dimensional loop space, and the extraction of `m_{0⁺⁺}` needs linearization around the
**unconstructed** vacuum. The obstruction is **not** the DSE tower; it is *surface/string non-uniqueness +
unconstructed master field + the 0-form/1-form channel mismatch*.

---

## §6 · Controlled approximations (what exists, and what it does/doesn't give)

- **Area + Coulomb two-parameter ansatz**: forced by (MM), but it is a *two*-function fit, not a closed
  scalar equation; it parametrizes `σ` and `c`, does not predict `m_{0⁺⁺}`.
- **Eguchi–Kawai reduction** (PRL 48 (1982) 1114): at `N=∞` the theory reduces to a matrix model, but the
  reduction **loses spacetime separation** ⇒ it carries the loop/`σ` content and **not** the `0⁺⁺` correlator
  (confirms §4's channel point).
- **`1/N` + loop equation**: gives the leading area law; the glueball spectrum still needs linearization.
- **Loop equation + large-N bootstrap** (recent): constrained numerics on `W(C)`, no closed gap formula.

No controlled approximation produces a closed scalar equation whose solution is `m_{0⁺⁺}(λ)`. `[O]`

---

## §7 · Does it predict `m_{0⁺⁺}/√σ`? No.

The loop equation alone predicts neither an absolute `m_{0⁺⁺}` nor the ratio. The lattice ratio
`m_{0⁺⁺}/√σ ≈ 3.5–3.9` is an **input**, not an output of (MM). `[NG]`

---

## §8 · Honest bottom line

- **Structurally better than DSE** `[T]`: the large-N MM equation is **closed in one functional `W(C)`**, with
  no infinite n-point tower.
- **But not a closed scalar gap equation** `[O]`: D=4 surface non-uniqueness makes `W(C)` a string functional
  (two independent functionals `c, σ`, no single scalar); the D=4 master field is unconstructed; the loop
  equation carries `σ` (1-form), and `m_{0⁺⁺}` (0-form) needs linearization around the unknown vacuum.
- **Exact missing step:** construct the D=4 master field `W(C)` (the vacuum solution), then linearize the
  fluctuation operator and read its lightest `0⁺⁺` eigenvalue. Equivalently, the `0⁺⁺` mass is on the far
  side of the same 0-form/1-form channel wall — the loop equation is the right *kind* of closed object,
  but it lives in the wrong channel for `m_{0⁺⁺}`.
- **Clay 0%; no QED.** Witness `w7/pnt_largen_loop_witness.py` runs clean.

*Refs:* Makeenko–Migdal, Phys. Lett. B88 (1979) 135; Nucl. Phys. B188 (1981) 269. Witten, master field,
in *Recent Developments in Gauge Theories* (1980). Singer, 2D master field (1990s). Lévy, *The master field
on the plane* (Astérisque 2017). Eguchi–Kawai, Phys. Rev. Lett. 48 (1982) 1114. Migdal, loop-space /
quantum-string linearization. Prior: `ROOF-LARGEN-MASTER.md`, `ROOF-CONFINEMENT-GAP.md`,
`SLOT-DSE-FIXPOINT.md`, `INVENT-SYNTHESIS.md`.
