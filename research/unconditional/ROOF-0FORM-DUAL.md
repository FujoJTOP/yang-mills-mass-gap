# ROOF-0FORM-DUAL — a literal 0-form 't Hooft algebra does **not** exist; the genuine 0-form "dual" is the **trace-anomaly Ward identity** (`tr F² ↔ dilatation current`), which degenerates to the gluon condensate `= H1`

> **VERDICT (no QED).** The roof program asks for a 0-form (0⁺⁺) analog of the Wilson–'t Hooft algebra
> `W(C)B(C′)=B(C′)W(C)e^{2πi·lk(C,C′)/N}` — a local `0⁺⁺` operator paired with a dual defect operator by a
> nontrivial phase. **It does not exist.** Three exact, independent obstructions (form-degree/linking,
> center neutrality, absence of a 0-form symmetry) make it empty. What **is** constructible is the *other*
> kind of 0-form dual: the `0⁺⁺` operator `tr F²` is the **source of the trace anomaly**, tied to the
> **dilatation current** `J_D` by an exact **Ward identity** `θ^μ_μ = (β/2g) tr F²` — a **linear anomaly
> equation with a c-number coefficient**, no phase, no linking. This degenerates to the **gluon condensate
> `⟨tr F²⟩ = H1`** (matching `BASE-FIELDMEASURE`). So the roof is built on the **anomaly**, not on an
> algebra; and the 1-form→0-form transfer is blocked **exactly at the linking/center step** — the same
> "non-transferable step" `BASE-MINIMAL-MODEL` flagged. **Clay 0%.**

---

## §0 · Labels, iron laws, acid test (carried)

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.
**IL-i** non-positivity · **IL-ii** nonabelian · **IL-iii** no small parameter. **Acid**: `=0` on free massless *and discriminatingly* on gapped-atomless `:φ²:`.
The gap object is the lightest `0⁺⁺` **KL atom weight** `f₁²=ν({m₁²})=|⟨0|tr F²|G⟩|²` of the local **0-form** operator `tr F²`.

---

## §1 · The question, stated precisely

In 4D pure `SU(N)`, the classical global structure is the **1-form center symmetry** `Z_N`. Its electric
charged operator is the Wilson loop `W(C)`; its magnetic dual is the 't Hooft loop `B(C′)` (a twist defect on
a surface bounded by `C′`); together they obey `[W(C),B(C′)] = W(C)B(C′)(e^{2πi·lk(C,C′)/N}−1)` — a **phase
= exponential of the linking number** `lk` of the two **loops**.

**Question.** Is there a 0-form analog — a **local `0⁺⁺` operator** (`tr F²`) paired with a **dual defect**
operator by a nontrivial `e^{2πi·…}` phase? If yes, the roof could be attacked with a 't Hooft-type order
parameter in the *right* channel; if no, we learn *why* the 3D Z₂ magnetic mechanism cannot transfer.

---

## §2 · Three exact obstructions (the 't Hooft analogue is empty)

### Obstruction 1 — form-degree / linking `[T]`
The 't Hooft phase is `exp(2πi·lk/N)`, and `lk` is a topological invariant of **two closed cycles**
(`W1`: Gauss integral gives `lk(Hopf)=−1.0000`, `lk(unlinked)=0`). A **local (0-form) operator has point
support** — a 0-cycle whose boundary is the **empty set**. No 2-surface has a point as its boundary, so
`lk(point, loop) ≡ 0`, and the phase is `1`. **A 0-form operator is 't Hooft-neutral.** Only extended
(1-form loops, 2-form surfaces) operators can carry a nonzero 1-form charge.

### Obstruction 2 — center neutrality `[T]` / [V-num]
`tr F² = Σ_a Σ_{μν}(F^a_{μν})²` is a **color singlet** (adjoint-valued, no free color index), hence has
**zero N-ality**. The center `z = e^{2πi/N}` acts on the algebra as `Ad_{zI}`; since `zI` is central,
`Ad_{zI}(X) = zIz̄I·X = X` (`W2`: `|Ad_{zI}(X)−X|_max = 2.2×10⁻¹⁶`, `tr F²` invariant to 6 decimals). By
contrast a **fundamental** Wilson loop `tr U → tr(zU) = z·tr U` carries the charge (`W2`). **So a 't Hooft
phase can appear only with a charged 1-form loop — never with the singlet `tr F²`.**

### Obstruction 3 — absence of a 0-form symmetry `[T]`
A 't Hooft-type algebra is the statement that two operators are charged under **dual** generalized
symmetries. The Hodge partner of a `p`-form symmetry in `D=4` is a `(D−p−2)`-form: `1↔1` (the center,
exists), but `0↔2`. A 0-form "dual" would need a **0-form symmetry** whose charged local operator pairs with
a **2-form** (surface) disorder defect. Pure YM has **no internal 0-form global symmetry** in the `0⁺⁺`
sector: (a) the gauge group is gauged — no flavor symmetry; (b) the θ-angle is `2π`-periodic — no continuous
0-form symmetry; (c) the dilatation "symmetry" is explicitly broken by the trace anomaly (`β(g)≠0`) and is
*anomalous*, not a symmetry. **The `0↔2` pairing is empty.** (`tr F²` is an operator, not a symmetry order
parameter — there is no symmetry for it to be the order parameter *of*.)

**Conclusion `[T]`.** No literal 0-form 't Hooft algebra exists. This is not a dynamical failure but a
structural one — it is the precise content of the "channel mismatch" `BASE-ORDER-PARAMETER` measured, now
proved at the level of form degree + center + symmetry.

---

## §3 · What IS constructible: the anomaly-dual `[T]`

The `0⁺⁺` operator is not symmetry-charged, but it is the **source of the trace (dilatation) anomaly** —
an exact operator Ward identity:
```
θ^μ_μ = (β(g)/2g) · tr F² ,      θ^μ_μ = ∂_μ J_D^μ  (J_D^μ = x_ν θ^{μν} the dilatation current) .
```
(convention: `F² ≡ F^a_{μν}F^{aμν}`; overall sign of `β/2g` aside). So in the `0⁺⁺` sector the **"dual
partner" of `tr F²` is the dilatation current `J_D`**, and the **"algebra" is the anomaly equation** — a
**linear** Ward identity with a **c-number coefficient** `β/2g`, *no phase, no linking*. This is a genuine
0-form dual structure, but of **anomaly type**, not 't Hooft type.

The `0⁻⁺` companion is analogous: `tr F∧F = ∂_μ K^μ` (Chern–Simons), the **chiral/topological anomaly**,
whose dual is the topological current `K^μ`; its vacuum curvature is the topological susceptibility
`χ_top` — the **Witten–Veneziano `0⁻⁺`** object, which `BASE-FIELDMEASURE` already showed is a **channel
mismatch** (wrong channel for the `0⁺⁺` gap).

**So the two 0-form operators split by anomaly type:**
| operator | channel | anomaly | dual partner | what it floors |
|---|---|---|---|---|
| `tr F²` | `0⁺⁺` | trace (dilatation) | `J_D` | **gluon condensate `⟨tr F²⟩ = H1`** |
| `tr F∧F` | `0⁻⁺` | chiral/topological | `K^μ` (Chern–Simons) | `χ_top` (wrong channel) |

**The `0⁺⁺` anomaly-dual degenerates to the gluon condensate** — exactly `BASE-FIELDMEASURE`'s "the source
can only be a `0⁺⁺` anomaly, and the unique one is the gluon condensate `= H1`." The roof round thus
**re-derives the same object from the operator side**.

---

## §4 · When does the 1-form mechanism transfer? — never, and here is exactly where it stops `[T]`

`BASE-MINIMAL-MODEL` showed the 3D Z₂ discriminator is the 't Hooft loop = **dual magnetization**, a
**0-form** object *there* because in 3D the dual of a 0-form symmetry is a 1-form and the dual Ising
magnetization is local. **The transfer to 4D blocks at exactly Obstruction 1+3:** in 4D the magnetic sector
is **1-form** (a loop), and a 0-form operator cannot link it (§2). So the 3D "magnetic order parameter" has
no 4D 0-form image — the "non-transferable step" is now proved, not just observed.

What *would* a genuine 0-form 't Hooft partner require? A **0-form symmetry** in pure YM (absent) whose
**2-form** disorder defect is a surface. No such symmetry exists at `θ=0` in the confined vacuum (`BASE-INDEX`
recorded the same: the 1-form center is a confinement order parameter, not an anomaly forbidding a gap). So
the roof cannot be a 't Hooft order parameter; it must be the **anomaly**.

---

## §5 · Acid test `[V-num]`

| adversary | value | passes? |
|---|---|---|
| free **massless** gauge (`F=0`) | `tr F² = 0` (`W3`) | ✔ returns 0 |
| gapped-atomless **`:φ²:`** (no gauge field; `β_g=0`) | gauge-anomaly source `= 0` (`W3`) | ✔ returns 0 |
| **interacting** gauge (`β(g)≠0`, `F≠0`) | `(|β|/2g)tr F² ≈ 0.117 ≠ 0` (`W3`) | ✔ fires |

**Acid-clean** (vanishes on both adversaries, fires only on the interacting gauge theory). **Honest caveat
`[O]`:** a nonzero anomaly source does **not** by itself force a gap — a *gapless* interacting theory can
carry a condensate. It passes the acid test but does not close the roof; its **floored magnitude is `H1`**.

---

## §6 · Verdict and consequence for the roof

- **The literal 0-form 't Hooft algebra does not exist** — three exact obstructions (§2), witnessed
  numerically (`W1`: linking is 1-cycle-only; `W2`: `tr F²` center-neutral). `[T]/[NG]`.
- **The genuine 0-form dual is the trace-anomaly Ward identity** `tr F² ↔ J_D` — a linear anomaly equation
  with c-number coefficient, **no phase**. `[T]`.
- **It degenerates to the gluon condensate `⟨tr F²⟩ = H1`** — the same roof object, re-derived from the
  operator side, consistent with `BASE-FIELDMEASURE`. `[C]`.
- **The 3D→4D transfer is blocked at the linking/center step** (a point links nothing; `tr F²` carries no
  N-ality) — the precise "non-transferable step." `[T]`.
- **Consequence:** the roof must be built on the **anomaly / gluon-condensate floor** (a measure statement),
  **not** on a 't Hooft-type order parameter. No new species; the target is again `H1`, now shown to be
  *forced* to be anomaly-type — the only 0-form object that is acid-clean, non-local, and in-channel.

**Clay 0%; no QED.**

*Refs:* `BASE-ORDER-PARAMETER.md` (channel mismatch), `BASE-MINIMAL-MODEL.md` (3D Z₂ mechanism, transfer
block), `BASE-FIELDMEASURE.md` (source = `0⁺⁺` anomaly = gluon condensate), `BASE-INDEX.md`; 't Hooft 1978
(Wilson–'t Hooft algebra); trace anomaly `θ^μ_μ=(β/2g)tr F²`. Witness: `w7/roof_0form_dual.py`.
