# ROOF-TRANSFER — the 3D Z₂ 0-form ridge is EXACT; the 4D span is blocked at TWO steps (form-degree jump `0→2` **and** *there is no transition to latch onto*)

> **VERDICT (no QED).** The 3D Z₂ specimen's magnetic-condensation mechanism **does** transfer as an
> **exact 0-form criterion** — but only into (i) genuinely 3D systems and (ii) *exact dimensional
> reductions* of 4D (small circle / single time slice), i.e. **deformed** theories. It does **not**
> transfer to genuine zero-temperature 4D pure YM, and the failure has **two** independent, exact steps,
> only the first of which was known: **(B1) form-degree / channel mismatch** — 3D's 't Hooft line dualizes
> to a **0-form spin**, 4D's dualizes to a **non-local 2-surface** (1-form symmetry); and **(B2) — new,
> and deeper — an order-parameter criterion needs a *transition*, but zero-T 4D pure YM is a *single* phase
> with only a *crossover* (no order parameter to condense).** Sharper than the base round: 3D Z₂ grows the
> *confinement order-parameter structure*, **not** the single-phase `0⁺⁺` mass gap. **Clay 0%.**

---

## §0 · Labels, iron laws, acid

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** this round.
**IL-i** non-positivity (fail on free massless) · **IL-ii** nonabelian · **IL-iii** no small parameter.
**Acid**: at the gapless point (β_c) the criterion must read "not gapped" (`m=0`, gap `=0`); and the pure-massless free field must give `0`.

---

## §1 · The construction: the exact 0-form criterion in the solvable specimen (3D Z₂ ≡ dual 3D Ising)

**The exact duality `[T]` (Wegner/Kramers–Wannier).** 3D Z₂ lattice gauge ↔ 3D Ising spin model. The
magnetic object — the **'t Hooft loop** `B(C)` (a line in 3D) — is **dual to the local Ising spin `σ`**,
a genuine **0-form** operator. The confining phase (`area law`, small gauge `β_g`) ↔ the **ordered**
Ising phase (`m ≠ 0`); deconfined ↔ disordered (`m = 0`). The 0-form "glueball/mass" channel is the
**energy operator `ε`**, whose correlation **gap** `Δ_ε` is the finite-size-safe 0-form gap.

**Exact 0-form criterion (the ridge).**
> *Gapped `⟺` `Δ_ε > 0`; and the two phases are distinguished by the 0-form order parameter `m`. The
> gapless point is the unique critical point `β_c`; **both phases are gapped.***

**Witness (exact, no MC):** 3D Ising transfer matrix on an `L×L=3×3` slice (`roof_transfer.py`), eigenvalues
classified by the global Z₂ parity:

| `β` | `Δ_ε` (even, = 0-form gap) | `Δ_σ` (odd) |
|---|---|---|
| 0.12 | 3.359 | 1.580 |
| 0.22 | **1.645** | 0.465 |
| 0.26 | **1.477** ← minimum (finite-`L` critical region) | 0.192 |
| 0.32 | 1.901 | 0.036 |
| 0.42 | 3.123 | 0.002 |

`Δ_ε` has a **clear minimum** near `β≈0.26` (thermodynamic `β_c≈0.2217`, finite-`L` shifted) and **rises
on both sides** ⇒ gapped away from criticality, gapless at it `[T/V-num]`. (`Δ_σ → 0` at large `β` is
**vacuum tunnelling** between the two ordered vacua on the finite torus — SSB, *not* a physical gapless
mode `[T]`.) The exact 2D-Ising magnetization curve (`m→0` at `T_c`, Onsager) calibrates the order-parameter
half. **The ridge is built and exact — but it is a criterion for a TRANSITION.**

---

## §2 · The exact reduction of 4D to the ridge (single time slice / small circle)

- **`L₄ = 1` `[T]`:** 4D Z₂ gauge on `T³×S¹` with a *single* time slice **is** 3D Z₂ gauge — the exact
  0-form criterion holds verbatim on that slice.
- **Small circle `L₄ ↓` `[C]`:** integrating out the short temporal direction gives an effective 3D Z₂
  gauge theory with `β_{3D}^{eff}(β, L₄)`; the 0-form criterion holds **in the reduced (deformed) theory**.

So the ridge **does** cover all exact reductions to 3D. The obstruction is *only* in leaving them.

---

## §3 · Where it breaks — two exact steps

**(B1) Form-degree / channel mismatch `[NG]`.** A 3D 't Hooft *line* bounds a surface whose dual is a
**point** (0-form spin). In 4D the 't Hooft *loop* dualizes to a **non-local 2-surface**, and the center
symmetry it diagnoses is a **1-form** symmetry. **There is no local 0-form dual in 4D** (consistent with
`BASE-ORDER-PARAMETER`: the 4D magnetic objects are 1-form, Schur-orthogonal to the neutral `0⁺⁺`).
Elitzur bars any local gauge-invariant order parameter; the only non-local `0⁺⁺` object (A3 Fredholm det)
is READ-not-FORCE. So B1 is the *known* channel wall.

**(B2) No transition to latch onto `[T]` — NEW.** The ridge criterion is an **order-parameter** criterion:
it certifies a *phase transition*. But **zero-temperature 4D pure SU(N) is a single confining phase**: as
the lattice coupling varies there is **no thermodynamic transition**, only the strong-to-weak **crossover**
at `β≈2.5`; the continuum theory is the `β→∞` limit inside that one phase. **An order parameter has
nothing to condense at a crossover.** The finite-`T` deconfinement transition (which *does* have an order
parameter — the Polyakov loop, a 1-form/center object) is a statement about **deconfinement**, not the
`0⁺⁺` mass gap; and decompactifying `T→0` **removes that transition while keeping the `0⁺⁺` gap**. So even
if the form degree were right, the criterion would be aimed at the wrong kind of object.

> **The base round's implicit claim is sharpened here:** 3D Z₂ grows the **confinement order-parameter
> structure** (a transition, a non-local order parameter, a gap-except-at-criticality), **not** the
> single-phase `0⁺⁺` mass gap. The analogy is exact for *confinement*, not for the zero-T gap. B2 is the
> precise reason the order-parameter route cannot reach Clay.

---

## §4 · The other candidate reductions — each reproduces the 1-form side, none the 0-form

| reduction | what it gives | where it breaks |
|---|---|---|
| **4D Z_N gauge** | self-dual `gauge↔gauge` (Wegner); a confinement transition exists | **no spin dual** ⇒ no 0-form order parameter; only 1-form (`'t Hooft` loop ↔ `Wilson` loop) |
| **Eguchi–Kawai (large-N, T⁴→point)** | reproduces **Wilson loops** = the **string tension** (1-form) | *spacetime separation is lost* (one site) ⇒ **no `0⁺⁺` mass**; solves the area-law, not the gap |
| **4D compact U(1)** | Coulomb phase | monopole **instantons are irrelevant in 4D** ⇒ massless photon; **acid-fails (massless)**; no gap |
| **small-circle `T³×S¹`** | exact 0-form criterion (via §2) | **deformed** theory; decompactification = the wall (`= I2 = H1`) |

**Every solvable 4D reduction delivers the 1-form/confinement datum; the 0-form `0⁺⁺` gap is the one thing
none of them carries.** This is the roof's honest coverage: it spans the reductions and the 3D ridge, and
**the 4D span is exactly the channel-`0` / no-transition gap.**

---

## §5 · What the roof does and does not cover

**Covers (exact):** the 3D Z₂ ridge (0-form criterion + 0-form gap `[T/V-num]`); every exact reduction of
4D to 3D (single slice, small circle `[T/C]`); a calibrated, transfer-matrix-checked instrument.

**Does not cover (the two broken steps):** (B1) the local/0-form dual does not exist in 4D (channel wall);
(B2) zero-T 4D has no transition for an order parameter to condense (crossover) — the deepest new point.

**What would actually complete the roof:** a criterion that lives in a **single phase** (no transition,
acid-clean at `β_c`) and in the **0-form** channel — i.e. a **field-measure functional** on `μ`
(`PROFILE` P1) / equivalently the gluon-condensate floor, **not an order parameter**. The base round named
exactly this (`BASE-SYNTHESIS §7`); ROOF-TRANSFER now proves the *order-parameter* route is doubly
off-target, so the only roof type left is the field-measure one.

---

## §6 · Honest bottom line

- **Constructed:** the exact 0-form condensation criterion and 3D gap (3D Z₂ ⊂ exact reductions), with an
  exact transfer-matrix witness. **This is a real, exact positive object.**
- **Blocked (two exact steps):** (B1) form-degree `0→2` at decompactification; (B2) **no transition** in
  zero-T 4D for an order parameter. B2 is new and sharper than the base round.
- **Corrects the base round:** 3D Z₂ models *confinement*, not the single-phase mass gap; its criterion is
  the wrong *kind* of tool for the zero-T `0⁺⁺` gap, independently of the channel wall.
- **Clay = 0%.** The ridge is exact; the 4D span is the channel-`0` / no-transition gap. **No QED.**

*Refs:* `BASE-SYNTHESIS.md`, `BASE-MINIMAL-MODEL.md`, `BASE-ORDER-PARAMETER.md`, `BASE-DUALITY.md`,
`PROFILE.md`; witness `w7/roof_transfer.py` (+ `.log`).
