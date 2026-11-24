# BASE-MINIMAL-MODEL — the smallest rigorously controlled model with the creature's anatomy: **3D Z₂ lattice gauge ≡ dual 3D Ising**

> **VERDICT.** The minimal model that has **all four** features — (a) a gapless fixed point, (b) a genuinely gapped
> phase, (c) a **provably non-local** discriminating order parameter that finite local data cannot read, (d) a
> minimal exclusion mechanism — is **3D Z₂ lattice gauge theory**, rigorously controlled through its **duality to
> the 3D Ising model**. The anatomy it exposes is exactly `PROFILE`'s: the discriminator is **magnetic / non-local**
> (the `'t Hooft` loop), and the mechanism is **defect (vortex/monopole) condensation**. This *explains* why all
> eight `SLOT-*` local classes were blind — none of them sampled the magnetic sector — and it names the 4D target
> as **magnetic-defect condensation at `α_s~O(1)`** = `H1 = I2`. **No QED; Clay 0%.**

---

## §0 · Why a *base* and not another slot attack

Eight `SLOT-*` rounds each built a `k→0` relation and died at the same step. The base round asks a different
question: **what is the smallest rigorous environment in which the creature actually stands?** If we can exhibit a
solvable model where a *non-local* order parameter separates a gapped phase from a gapless point, we obtain the
*transferable anatomy* (the shape of the discriminator and the form of the mechanism) — not a 4D proof, but the
face's skeleton on a stand.

---

## §1 · Candidate table (feature × model) — honest

Legend: ✔ has it (rigorous) · ✘ lacks it · ~ only formally / ill-defined.

| Model | (a) gapless fixed pt | (b) gapped phase (spectral) | (c) NON-LOCAL discriminator, local data blind | (d) mechanism | transfer to 4D? |
|---|---|---|---|---|---|
| **2D pure YM** (exact, Migdal/Witten) | ✔ free limit | ✔ exact (group-Laplacian gap) | **✘** — one coupling: local `⟨tr F²⟩` and string tension `σ` are both functions of `e²`; **no genuine local/non-local split** | kinematic: compactness (no propagating d.o.f.) | **✘** wrong mechanism (kinematic, not dynamical) |
| **Large-N single plaquette / Gross–Witten** | ✔ (3rd-order pt) | ~ no Hamiltonian ⇒ "mass gap" ill-defined | ✔ eigenvalue density is non-local in entries | emergent eigenvalue support | ✘ not a QFT spectrum |
| **0+1d gauge QM** (particle on `G`) | ✘ (never gapless) | ✔ discrete (compact) | — | compactness | ✘ no gapless/gapped dichotomy |
| **3D Z₂ gauge ≡ 3D Ising** | ✔ critical point `β_c` | ✔ confined (=Ising ordered) phase, exponential clustering `[T]` | **✔** `'t Hooft` loop; **Elitzur** ⇒ no local gauge-invariant order parameter `[T]` | **magnetic vortex condensation** (dual magnetization) | **partially** — anatomy ✔, mechanism conjectural in 4D |

**⇒ 3D Z₂ gauge is the unique candidate with all four.** (2D YM is the cleanest *exact* solvable but fails (c)
genuinely: with a single coupling there is no real local-vs-non-local separation; its gap is *kinematic*.)

---

## §2 · The anatomy in 3D Z₂ gauge (with the duality dictionary)

3D Z₂ lattice gauge is dual to the 3D Ising model (Wegner 1971). The dictionary is the whole point:

| gauge theory (electric) | dual Ising (magnetic) |
|---|---|
| gauge field `U_ℓ ∈ {±1}` on links | Ising spin `s` |
| **local** gauge-invariant: plaquette `U_p`, small Wilson loop | local energy `s_i s_j` |
| **Wilson loop** `W(C)` (electric) | dual high-T (domain-wall) expansion |
| **`'t Hooft loop** `T(C)` (magnetic) | **magnetization `m`** ⟵ the order parameter |
| confined phase (gapped) | **ordered** phase (`m≠0`) |
| deconfined (gapless) | **disordered** phase (`m=0`) |
| gap ⇔ exponential clustering | gap ⇔ exponential clustering `[T]` (Glimm–Jaffe–Spencer) |

**Feature (c), rigorously.** Elitzur's theorem `[T]`: a local gauge symmetry cannot break spontaneously, so **no
local gauge-invariant observable can be an order parameter**. The confinement order parameter must be **non-local**
— it is the `'t Hooft` loop (or the vortex free energy), the dual image of the Ising magnetization. This is a
**rigorous realization of `PROFILE` P1**: the discriminator lives on the field measure in a way **no local /
spectral datum can read**.

**Witness (illustration, not proof).** `w7/base_minimal_model.py` (3D Ising, `L=8`, Metropolis, `<2 s`):

```
  beta     |m| (dual 't Hooft)      E (local / plaquette)   phase
  0.1500   0.0681                   -0.4951                DISORDERED (deconfined/gapless)
  0.1900   0.1265                   -0.6949                DISORDERED
  0.2100   0.2196                   -0.8787                critical
  0.2217   0.3710                   -1.1011                critical
  0.2350   0.6362                   -1.5739                ORDERED (confined/gapped)
  0.2600   0.7996                   -2.0743                ORDERED
  0.3000   0.9065                   -2.5158                ORDERED
```

The **order parameter** `|m|` (the non-local `'t Hooft` / dual-magnetization object) turns on sharply below
`β_c≈0.2217`, while the **local** energy `E` is smooth through `β_c`. The rigorous "local data cannot read it"
statement is Elitzur `[T]`; the numerics only illustrate. **Not a proof about 4D.**

---

## §3 · The minimal exclusion mechanism, and exactly what fails to transfer

**Mechanism `[T]` in 3D Z₂.** The gapless (critical) point is excluded on the ordered side because the **dual
magnetization `m≠0`** — equivalently, **magnetic vortices condense** — which gaps the electric spectrum by
exponential clustering. The schematic exclusion relation is

```
κ = 0  (gapless / critical)  is a solution  of  the symmetric relation   ∀β,
but for β > β_c   the dual order parameter m ≠ 0  INHOMOGENEOUSLY gaps the spectrum:
        m ≠ 0   ⇒   ξ < ∞   ⇒   Δ > 0 .
```

This is **exactly the geometry the eight `SLOT-*` routes kept demanding and never had**: an **inhomogeneous,
sign-definite source at the zero**, supplied here by a **non-local magnetic order parameter**.

**What transfers to 4D pure `SU(N)`.** The **anatomy** transfers:
- the discriminator must be the **`'t Hooft` loop / center vortex** (magnetic, non-local) — Elitzur holds in 4D too `[T]`,
  so **no local electric datum can be the order parameter**;
- the mechanism must be **magnetic-defect condensation** (the dual-superconductor picture), **not** a positivity
  inequality and **not** a local Ward identity.

**What does NOT transfer — the exact missing step.** In 3D Z₂ the defect condensation is **proven** (Peierls /
duality to Ising + cluster expansion). In 4D `SU(N)` at `α_s~O(1)` there is **no proof** that the magnetic
monopoles/vortices condense, equivalently that the `'t Hooft` loop obeys an area law. **That unproven step is
`H1 = I2`.** The base therefore **isolates the wall** to one named statement rather than "the whole theory":

> **4D target = "the `'t Hooft` loop / magnetic-defect condensate is nontrivial at the crossover."**

---

## §4 · What the base buys us (the actual progress)

1. **A [T]-backed explanation of the eight `SLOT-*` no-gos.** The discriminator is **magnetic / non-local**
   (Elitzur `[T]`). Every `SLOT-*` local class — moments, single PF cone, scheme-clean DSE, variance `Hess`,
   quadratic minorant, bootstrap moments — is **electric / local**, so it **cannot** read it. Their blindness is
   not coincidence; it is the theorem.
2. **The correct *type* of missing object, concretely.** Not a scalar function of local correlations, but a
   **non-local magnetic order parameter** (`'t Hooft` loop / vortex free energy). This is `PROFILE` P1 made
   explicit: a **field-measure (non-local) functional**, not a spectral statistic.
3. **The correct *form* of the mechanism.** **Defect condensation** (dual-superconductor), i.e. an **inhomogeneous**
   gapping of the electric spectrum — matching `SLOT-SYNTHESIS` §4's "(S3) inhomogeneous source at the zero".
4. **A sharp 4D target.** "Prove magnetic-defect condensation / `'t Hooft` area law at `α_s~O(1)`" — a **known,
   named** open problem, the standard dual-superconductor route to confinement.

---

## §5 · Is this progress or a relabel? — honest

It is **structural progress of a specific kind**: a rigorous lower-dimensional model **derives** (not asserts) two
things the campaign only conjectured — (i) *why* the discriminator is non-local (Elitzur + duality), and (ii) *why*
every local positivity/scheme tool was blind (they are all electric). It converts "the whole theory is
intractable" into **one named, structured target** (magnetic-defect condensation). It does **not** prove 4D: the
condensation step is itself open, and this base cannot supply it. **A foundation, not a roof.**

*Labels:* duality & phase structure & gap in 3D Z₂/Ising `[T]` (Wegner 1971; Elitzur 1975; Glimm–Jaffe–Spencer
1976); 2D YM exactness `[T]` (Migdal 1975; Witten 1991); Gross–Witten `[T]` (Gross–Witten 1980; Wadia 1980);
the 4D transfer & defect condensation `[O]` (open) = `H1=I2`. Witness = illustration only `[V-num]`.
**No QED. Clay 0%.**

*Refs:* `PROFILE.md`; `SLOT-SYNTHESIS.md`; `KEY-DIRECT-DERIVATION.md`; witness `w7/base_minimal_model.py`.
