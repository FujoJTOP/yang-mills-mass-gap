# PNT-PERCOLATION-CHI — the "level-set percolation ⇒ floored χ" inequality is **FALSE as stated**, and the true percolation inequality runs the **WRONG direction** for `H1`; percolation is not the strut

> **VERDICT (no QED).** The proposed strut — a percolation event (a macroscopic cluster of high-`tr F²`
> plaquettes) forces `χ = Var(tr F²) ≥ c > 0` — **fails**, for three independent reasons. **(A) `[T]`
> Category error:** percolation is a property of a *single configuration*, `χ` is an *ensemble* variance; a
> deterministic percolating field has `Var = 0` (witness). **(B) `[T]` Wrong direction:** the true
> percolation inequality (FK/GKS) is `χ ≥ E|C₀|` (mean cluster size) — percolation drives `χ` **large**,
> which is the **gapless** signal; `H1` needs the **opposite**, an *upper* bound `χ ≤ 1/λ`. **(C) `[NG]`
> FKG absent** for full-rep `SU(N)` (character weights are not ferromagnetic), so even the sign-only
> percolation step is unavailable. The correct-shaped lemma — **"no percolation ⇒ `χ ≤ c`"** — is true in
> principle but its hypothesis ("no percolation at `β=2.5`") **is the gap**. Percolation is **not the
> strut**; it is at best a detector of the gapless side. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.
Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The two levels: configuration vs ensemble (the bare implication is a category error)

`χ := Var(tr F²)` (or `Var` of the coarse mode) is an **ensemble** quantity — a second moment of the
measure `μ`. "Level-set percolation" `{A_θ has a spanning cluster}`, `A_θ := {p : tr F²(p) > θ}`, is a
property of **one** field configuration. A single configuration carries **no variance**.

**Counterexample `[T]` (witness §A).** Take the degenerate measure `δ_{φ≡1}`: the level set is the whole
torus, so it **percolates with probability 1**, yet `Var(Φ) = 0` exactly. More generally: any measure
supported on percolating configurations has `P(percolation)=1` and `Var` set by the *spread among those
configurations*, not by the percolation event. So

> **percolation (a config-level {0,1} event) does NOT imply `χ = Var ≥ c` (an ensemble second moment).**

The implication is not merely unproven — it is **ill-typed**. (This is the `PNT` analogue of the corpus's
recurring "gap-detector ≠ atom-certifier": a detection event is not a certificate.)

---

## §2 · The true percolation inequality runs the wrong way for `H1`

Where percolation *does* constrain a susceptibility, it is via **FK/GKS** on a ferromagnetic (spin)
representation: `⟨σ₀σ_x⟩ ≥ P(0↔x)`, hence

> **`χ = Σ_x ⟨σ₀σ_x⟩ ≥ Σ_x P(0↔x) = E|C₀|`   (mean FK cluster size).**

So **percolation ⇒ mean cluster size diverges ⇒ `χ` = ∞** — i.e. percolation is the **gapless** signal.
The gap (`H1`) requires the **opposite**: an **upper** bound `χ ≤ 1/λ` (then `m₁² ≥ f₁²/χ ≥ f₁²λ`). Thus:

> **Even if the percolation inequality were true, it would prove `χ` LARGE = gapless — the negation of
> what the gap needs.**

**Witness §B (exact, 1×4 Ising ring, `q=2` FK, `p=1−e^{−2K}`).** `χ = Σ_x⟨σ₀σ_x⟩` and `E|C₀|` coincide and
**grow together** with `K`: `(χ, E|C₀|) = (1.487, 1.487) → (2.135, 2.135) → (2.810, 2.810) → (3.669, 3.669)`
at `K = 0.2, 0.4, 0.6, 1.0`, with connectivity `P(0↔2)` rising `0.078 → 0.868`. Percolation and `χ` rise
**together** — a lower bound, the gapless direction.

---

## §3 · "Percolation" is θ-dependent and not a single event

The level-set event depends on the threshold `θ`. At generic/low `θ` the level set percolates **trivially**
even for a massive (gapped) field (a random level set occupies `~50%` of plaquettes at `θ≈median` and
percolates in `d≥2`). At **high** `θ` the question "does `{tr F²>θ}` percolate?" **is** the subcriticality
of the high-amplitude clusters — which is exactly the statement that the field's large fluctuations are
**finite-range**, i.e. that `χ` is finite, i.e. that the theory **gaps**. So the only *non-trivial* form of
the percolation question is already the gap, and the inequality as stated (which does not fix/quantify `θ`)
is either trivially true or the gap itself.

---

## §4 · FKG is unavailable for full-rep `SU(N)`

The one route that would give *any* quantitative percolation→`χ` link is FK/FKG, which requires a
**ferromagnetic / totally-ordered spin** representation. The full-rep `SU(N)` Wilson weight is a **character
sum** `Σ_r a_r(β) χ_r(U)` — no `q`-state, no spin representation, no FKG (`BASE-FIELDMEASURE` R2;
`INVENT-PERCOLATION` (ii)). So even the sign-only step `Cov(1_P, Φ) ≥ 0` is **unavailable** for the
physical gauge measure.

---

## §5 · The correct-shaped lemma, and why it re-encodes `H1`

The **right direction** is the contrapositive, and it is a standard fact:

> **"the level set `{tr F²>θ}` is *subcritical* (no spanning cluster, finite mean cluster size) ⇒ `χ ≤ c`."**

This is a genuine **upper** bound (what `H1` wants). But its hypothesis — **"no percolation at `β=2.5`"** —
is *exactly* the statement that the high-amplitude `0⁺⁺` regions do not percolate macroscopically, which is
**equivalent to the gap**. Proving subcriticality of the level set at the crossover is not a free lemma; it
*is* `H1` (in the percolation coordinate).

**Exact failing step `[O]`:** a **uniform subcriticality bound** for the high-`tr F²` level-set clusters at
`β=2.5` (equivalently a uniform upper bound on the mean cluster size `E|C₀| ≤ c`). This is the strut, and it
is not supplied by percolation theory for the physical (non-ferromagnetic, full-rep) gauge measure.

---

## §6 · Verdict

- **The proposed inequality `percolation ⇒ χ ≥ c` is FALSE** (config vs ensemble, §1) **and points the
  wrong way** for `H1` (§2); FKG is absent (§4).
- **The correct-shaped lemma `no-percolation ⇒ χ ≤ c` is real but its hypothesis is the gap (§5).**
- **Percolation is NOT the strut.** It is, at best, a **detector of the gapless side** (percolation ⇒ `χ`
  large), the mirror image of what the gap needs.
- **What survives:** the percolation *picture* localizes `H1` to a single, clean analytic question —
  **uniform subcriticality of the high-`tr F²` clusters at `β=2.5`** — but does not reduce its difficulty.

*Labels:* §1 category error `[T]`; §2 FK/GKS direction `[T]`/`[V-num]` (witness §B exact); §4 `[NG]`; §5
`[O]`. Acid passes trivially (percolation of a free Gaussian level set is generic ⇒ non-discriminating).
**No QED; Clay 0%.**

*Witness:* `w7/pnt_percolation_chi.py` (runs clean, `py -3.12 -I`, exact 1×4 Ising/FK enumeration).
*Refs:* `INVENT-PERCOLATION.md`, `BASE-FIELDMEASURE.md`, `ROOF-UNIFORM-CHI.md`, `INVENT-SYNTHESIS.md`.
