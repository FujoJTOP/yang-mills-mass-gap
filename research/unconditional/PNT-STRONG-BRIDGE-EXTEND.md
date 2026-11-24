# PNT-STRONG-BRIDGE-EXTEND — `m_{0⁺⁺} ≥ 4√σ` is physically plausible to `β=2.5` but PROVEN only at strong coupling; `β*` cannot be pushed rigorously

> **VERDICT (no QED).** The strong-coupling bridge `m_{0⁺⁺} ≥ 4√σ` (ROOF-CONFINEMENT-GAP) is a *local*
> statement: both the 0⁺⁺ mass and the string tension are computed in the **same convergent polymer
> expansion** in the weight `u(β)=I₂(β)/I₁(β)`, so their ratio `4√ln(1/u)` is a controlled lower bound —
> but only where that expansion converges. Four extension routes were tried. **(i) Polymer/KP resummation**
> is pinned at the KP/SZZ radius (`β*≈8.9×10⁻³`, benchmark `1/48`); constant optimization buys `<0.5`
> orders. **(ii) Padé/conformal interpolation** gives a *positive heuristic* interpolant (tree formula
> `m/√σ=3.30` at `β=2.5`, lattice `≈3.7`) but **no rigorous error bound** — a numerical estimate, not a
> certificate. **(iii) Alternative polymers/contours** hit a structural wall: the activity at `β=2.5` is
> `O(0.5)`, so no small-activity polymer exists there. **(iv) A monotone quantity** has no proven
> `β`-monotonicity. **`m_{0⁺⁺} ≥ 4√σ` is therefore proven at strong coupling, physically plausible to `β=2.5`,
> but not proven there. `β*` stays `≤8.9×10⁻³`. Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.
Acid: `≤0`/`=0` on free massless *and* discriminatingly on `:φ²:`.

---

## §1 · The bridge and why it is strong-coupling-local

From `ROOF-CONFINEMENT-GAP`: area law `σa² = −ln u + O(u)` (Osterwalder–Seiler) and the 0⁺⁺ gap
`m₁ a = 4 ln(1/u)` (P3 four-wall tube, Gupta hep-lat/9807028; Münster) share the **same polymer weight**
`u = I₂(β)/I₁(β)`. Hence

```
m_{0⁺⁺} / √σ = 4 √( ln(1/u) )  ≥ 4   ⟺   u ≤ e^{-1} ≈ 0.368.
```

This is a genuine `[C]` inequality **where the polymer expansion converges** (two convergent series sharing
`u`). It is *not* transportable: at large `β` the confinement/area law is non-perturbative (only the
asymptotic large-`β` area law of Cao–Nissim–Sheffield exists, which gives no transfer-matrix gap), and the
polymer expansion in `u` diverges. So the bridge's *proof* is local, even though its *formula* extrapolates.

---

## §2 · Four extension attempts (the numbers)

`u(β) = I₂(β)/I₁(β)`, SU(2), from `w7/pnt_strong_bridge_extend.py`:

| β | u(β) | `4√ln(1/u)` | note |
|---|---|---|---|
| `8.9×10⁻³` (residue β*) | 0.0022 | 9.89 | rigorous regime |
| `1/48` (SZZ benchmark) | 0.0052 | 9.17 | rigorous regime |
| 0.5 | 0.124 | 5.78 | |
| 1.0 | 0.240 | 4.78 | |
| 1.625 | 0.368 | 4.00 | `u=1/e` boundary |
| **2.5** (crossover) | **0.507** | **3.30** | lattice `≈3.7` |

`u(2.5)/u(1/48) ≈ 97` — the expansion parameter is ~2 orders of magnitude outside the convergence radius.

**(i) Polymer / Kotecký–Preiss resummation `[NG]`.** The rigorous convergence thresholds are `β<1/48`
(Shen–Zhu–Zhu `[T]`) and the residue-specific `β*≈8.9×10⁻³` (BASE-INFINITE-VOL). Constant optimization
buys `<0.5` orders (BASE-INFINITE-VOL); it cannot close a ~2-order gap to 2.5.

**(ii) Padé / conformal interpolation `[C-heuristic]`.** The tree formula stays **positive** for all `u<1`
(`3.30` at `β=2.5`), so a positive interpolant exists. A Padé test on a mock series with the KP-scale radius
(`f=(1+x/u*)^{-1/2}`, `u*=0.1`) gives diagonal `[M/M]` accurate to ~1% and off-diagonal `[L/M]` spread ~14%
at `u=0.5`: Padé can be numerically good but has **no rigorous error bound**. So interpolation yields a
*plausibility* argument (`m/√σ ≈ 3.3`, consistent with lattice `3.7`), **not a proof**.

**(iii) Alternative polymer/contour `[NG]`.** At `β=2.5` the fundamental plaquette/character weight is
`O(0.5)` (character ratio `a₁/a₀≈0.57`, ROOF-PEIERLS-0FORM). Any polymer or contour built on the character
expansion has activity `O(1)`; there is no polynomial-complexity polymer with small activity, so no
Kotecký–Preiss radius reaches it. This is structural, not a constant to tune.

**(iv) Monotone quantity `[NG]`.** No proven `β`-monotonicity of the residue/atom weight (ROOF-TRANSPORT-SCALAR;
empirically the lattice residue even *dissolves* at `β≈1` as a resolution artifact). So no "monotone
continuation" argument survives.

---

## §3 · What is true, exactly

- `[T]` **Rigorous:** `m_{0⁺⁺} ≥ 4√σ` for `β < 1/48` (and `β*≈8.9×10⁻³` for the residue-specific version).
- `[C-plausible]` **Not rigorous:** the same ratio formula gives `≈3.3` at `β=2.5`, consistent with lattice
  `≈3.7`; a Padé interpolant stays positive. This is evidence, not proof.
- `[NG]` **Rigorously out of reach by these four methods:** no polymer/contour/monotone/Padé argument
  extends the *proof* to `β=2.5`. The exact step is the same UV/IR-separation / KP-radius wall
  (`χ_lattice ≲ a²`) in the confinement-mass coordinate.

---

## §4 · Honest bottom line

- **The bridge cannot be extended rigorously to `β=2.5`.** `β*` stays `≤8.9×10⁻³`; the gap to 2.5 is ~2
  orders and irreducible by polymer resummation, alternative contours, or monotone arguments.
- **But the bridge is not vacuous:** it is a genuine rigorous strong-coupling theorem (`m≥4√σ`), and its
  ratio formula is *physically* plausible all the way to `β=2.5` (`3.3` vs lattice `3.7`) — so it is a useful
  *semiclassical/strong-coupling anchor*, not a proof of the crossover gap.
- **The exact missing step** is unchanged: a rigorous, scheme-independent lower bound on the 0⁺⁺ mass at
  `β≈2.5`, uniform in `(L,a)` — equivalently `χ_lattice ≲ a²`, the same strut as the rest of the campaign.
- **Clay = 0%.** No new species closes the gap; the strong-coupling bridge is confirmed local.

*Witness:* `research/unconditional/w7/pnt_strong_bridge_extend.py` (py -3.12, numpy+scipy, runs clean).
*Refs:* ROOF-CONFINEMENT-GAP, BASE-INFINITE-VOL, ROOF-PEIERLS-0FORM, ROOF-TRANSPORT-SCALAR, P3
(PROOF-P3-UNCONDITIONAL), Osterwalder–Seiler 1978, Shen–Zhu–Zhu CMP 400 (2023), Gupta hep-lat/9807028.
