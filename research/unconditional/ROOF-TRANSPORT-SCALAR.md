# ROOF-TRANSPORT-SCALAR — the scalar transport reduces to ONE monotonicity statement (`R_vac` numerically non-decreasing in `β`); the naive scalar `1/√χ` is the wrong variable, `β`-continuity is acid-dirty, and the weak-coupling IR anchor does not exist

> **VERDICT (no QED).** The base reduced the roof to one positive scalar `R_vac ∈ (0,1]` (equivalently the
> residue `f₁² = W_tube·R_vac`). We test the *transport* of this scalar across the crossover in its minimal
> form. Three structural facts first: **(i)** the lattice-unit variable `1/√χ = m·a` is the **wrong** transport
> variable — it `→0` as `β→∞` at fixed *physical* gap (a scaling fact, not a gap loss); **(ii)** the transport
> is **one-sided** — the strong side (`β<1/48`, polymer) has a rigorous *IR* floor, but the weak side
> (Balaban) is *UV* control, **not** an IR lower bound, so there is **no right anchor** to match to;
> **(iii)** `β`-continuity (compactness on a finite lattice) gives only the **free-field box gap**
> (`4sin²(π/L)`), which is **acid-dirty** and `→0` as `L→∞`. With `1/√χ` excluded and continuity vacuous,
> **the entire transport collapses to a single monotonicity statement**: *if `R_vac` (equivalently `f₁²`) is
> non-decreasing in `β`, the strong-coupling floor extends to all `β ≥ 1/48` and the roof closes* (modulo
> isolation I2 and `a→0`). **That monotonicity is unproven, empirically obscured by a lattice-resolution
> artifact (`β≈1` "dissolution"), and no correlation inequality supplies it.** The scalar transport is the
> wall in its sharpest *monotone* coordinates — pinned, not closed. **Clay 0%.**

---

## §0 · Labels, iron laws, acid test (carried)

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.
**IL-i** non-positivity (fail on free massless `f₁²=0`). **IL-ii** nonabelian. **IL-iii** no small parameter.
**Acid**: `≤0`/`=0` on free massless *and discriminatingly* on gapped-atomless `:φ²:`.

---

## §1 · Which scalar is transportable? `R_vac`, not `1/√χ` `[T]`

The base offers two "single scalars"; only one can be transported:

- **`1/√χ = G″ = m·a` (lattice units)** — **WRONG.** Hold the *physical* gap `m` fixed and let `a(β)` carry
  the asymptotic-freedom scaling (`a ~ e^{−cβ}`); then `m·a → 0` as `β→∞` (witness §A). So "`1/√χ` stays
  `>0` across `[1/48,∞)`" is unachievable **for a scaling reason** (the lattice unit shrinks), *not* because
  the gap is lost. Transporting `1/√χ` would either be trivially true or trivially false; it carries no
  `a→0`-invariant information.
- **`R_vac ∈ (0,1]`** — **correct.** It is a *pure* partition-function ratio `Z_tube/Z_0`, hence
  `a→0`-invariant (both `Z`'s share the UV divergence; the defect is local). And `f₁² = W_tube·R_vac` with
  `W_tube = a₀u^{4/a}` explicit and `β`-monotone (`BASE-INFINITE-VOL §1`). **All the β-dependence that
  decides survival is in `R_vac`; all the β-dependence of magnitude is in `W_tube`.**

⇒ **The transport variable is `R_vac` (equivalently `f₁²`).** The naive `1/√χ` must be dropped.

---

## §2 · The transport is one-sided: no rigorous weak-coupling IR anchor `[T]`

`BASE-PHASE-RIGOROUS` split the base into two solid ends — but they are **not two IR anchors**:

| side | what is rigorous | IR scalar lower bound? |
|---|---|---|
| strong `β<1/48` (polymer/KP, SZZ) | `f₁²_L ≥ r₀(β)>0` `[T]`, i.e. `R_vac ≥ 1−O(u)>0` | **YES** (the floor) |
| weak `β→∞` (Balaban) | UV stability / renormalization `[T]` | **NO** — Balaban controls UV divergences, not the IR mass |

So the "two-sided matching" has **only one anchor** (the strong-coupling floor). The weak side gives *no* IR
lower bound on `R_vac`. The transport is therefore **one-sided**: extend the left floor rightward. There is
no right-hand value to match to, and the middle `[1/48, 2.5]` is uncovered by any rigorous IR bound. **This
already rules out "two-sided matching" as a route** — there is no second side.

---

## §3 · `β`-continuity is free but acid-dirty `[T]`

`R_vac(L,β)` is continuous in `β` (partition functions are continuous for finite `L`). On a *finite* lattice
and the compact interval `[1/48, 2.5]`, a continuous strictly-positive function attains `min > 0` — so
"positivity across the interval" is **trivially true at fixed `(L,a)`**. But:

- the guaranteed positive value is the **free-field box gap** `4sin²(π/L)` (witness §B: `2.0→0.59→…→2.4×10⁻³`
  as `L=4→…→128`), which is **exactly the free field's gap** — **acid-dirty** (it "proves" a gap for the
  massless free field);
- and it `→0` as `L→∞`.

So `β`-continuity is the same vacuity as the PF single-cone (`BASE-PF` acid failure): it certifies a box gap,
never the thermodynamic gap. **Transport by continuity is vacuous; it cannot be used.**

*(This is the scalar-coordinate form of `BASE-PF`'s `[T]` sign-blindness and `BASE-NONLOCAL-BOUND`'s
"variational = upper only.")*

---

## §4 · The reduction: the transport is ONE monotonicity statement `[C]`

Excluding `1/√χ` (§1), the vacant right anchor (§2), and `β`-continuity (§3), the only remaining structure
is **monotonicity**. Write `R_vac = e^{−ΔF(β)}` (`ΔF≥0` the *tube defect free energy*). Then on
`[1/48, 2.5]`:

- **If `R_vac` is NON-DECREASING in `β`** (`ΔF` non-increasing): `inf_{β≥1/48} R_vac = R_vac(1/48)` = the
  **strong-coupling value ≥ 1−O(u) > 0** (rigorous). Then `f₁²(β) = W_tube(β)·R_vac(β) ≥ W_tube(β)·R_vac(1/48) > 0`
  for **all** `β ≥ 1/48` (witness §C, case 1). **The floor transports; the roof closes** (mod I2, `a→0`).
- **If `R_vac` is NON-INCREASING** (`ΔF` non-decreasing): `inf = R_vac(2.5)`, which needs a **right anchor**
  that does not rigorously exist (§2) — transport fails.

> **Theorem T1 (the transport in scalar form).** *A non-decreasing `R_vac` in `β` is SUFFICIENT to extend
> the strong-coupling residue floor across the crossover.* `[C]` — the implication is exact; the hypothesis
> is the open step.

**This is the sharpest possible reduction of the "transport machine":** not "build a rigorous 4D gauge RG,"
but **prove one monotonicity** — `R_vac` (equivalently `f₁²`) does not decrease in `β`.

---

## §5 · The exact failure `[NG]/[O]`

The monotonicity of T1 fails to be established for three independent reasons:

1. **Unproven `[NG]`.** No monotonicity of a Källén–Lehmann atom weight (`f₁²`) or of a defect free-energy
   ratio in `β` is proved for nonabelian lattice gauge. The multiplicative residue flow
   `∂_k f₁² = 2γ_O f₁²` (`SLOT-RICCATI`) is **linear, with `0` absorbing** — it *preserves* a sign but
   **cannot create** monotonicity in `β`.
2. **Empirically obscured `[V-num]`.** The *lattice-measured* `f₁²` is **non-monotone** — it rises then
   "dissolves into the cutoff" at `β≈1` (witness §D: `0.06→0.42→0.08`). This is a **lattice-resolution
   artifact** (the atom leaks into the cutoff continuum at coarse spacing), so it does *not* disprove
   monotonicity of the *physical* `R_vac` — but it **blocks any monotonicity proof anchored on the measured
   data**, and it is the only data we have.
3. **No inequality supplies it `[NG]/[O]`.** The natural tools all fail for `R_vac`:
   FKG/Griffiths need *ferromagnetic* positivity, absent for the full-irrep Wilson action (`BASE-INFINITE-VOL §2(D)`);
   Gaussian domination is an *upper* bound; and — decisively — `ΔF` is a **defect free energy** that tracks the
   **surface tension / confinement**, so **`ΔF` monotone in `β` is essentially the confinement statement = the
   target** (circular). Witness §C shows the two monotone cases cleanly; neither is provable.

**Exact failing step.** The hypothesis of T1 — *"the tube defect free energy `ΔF(β)` is non-increasing in
`β`"* — is an unproven, data-obscured, and (in the only available coordinates) **circular** statement. That
single monotonicity **is** the transport wall.

---

## §6 · Honest verdict

- **The scalar transport reduces to ONE monotonicity** (Theorem T1): `R_vac` non-decreasing in `β` ⟹ the
  strong-coupling residue floor extends across the crossover. This is the minimal form of the charter's
  "transport machine," and it is a genuine reduction: it replaces "build a 4D gauge RG" with "prove one
  monotonicity."
- **But the monotonicity does not close**: it is unproven, the measured residue is non-monotone by a
  resolution artifact, and the only structural handle (`ΔF` vs `β`) is the confinement statement itself.
- **Two scalars excluded along the way**: `1/√χ` (wrong — scaling `→0`), and `β`-continuity (acid-dirty —
  box gap). The transport variable must be the `a→0`-invariant `R_vac`; and the transport is **one-sided**
  (no rigorous weak-coupling IR anchor).
- **Scorecard:** (1 nonperturbative ✓) (2 sign-correct lower bound: T1 *shape* ✓ / monotonicity ✗) (3
  non-circular: implication ✓ / hypothesis ✗ `=confinement`) (4 fails free field ✓ — T1 uses the strong floor,
  not the box gap) (5 uniform `(L,a)` ✗ `=I2`). **2/5**, now in *monotone* coordinates.
- **Clay 0%; no QED.** The roof is not laid; but the transport is now **one named monotonicity** — the
  sharpest, most attackable form it has taken.

*Refs:* `BASE-INFINITE-VOL.md` §1–§2, `BASE-NONLOCAL-BOUND.md` §1–§4, `BASE-PHASE-RIGOROUS.md`,
`BASE-PF-FIXPOINT.md`, `SLOT-RICCATI-FLOW.md`, `SLOT-SYNTHESIS.md`, `BASE-SYNTHESIS.md`; witness
`w7/base_transport_scalar.py`.
*Status: **transport = one monotonicity [C]; monotonicity [NG]/[O]; Clay 0%.***
