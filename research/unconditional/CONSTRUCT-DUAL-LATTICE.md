# CONSTRUCT-DUAL-LATTICE — the vortex-ring 0⁺⁺ is a **continuum edge** (`f₁²=0`), not a KL atom: the dual superconductor delivers a renormalon-free **mass** but no **atom**

> **VERDICT (no QED).** Constructing the dual-superconductor `0⁺⁺` concretely shows it splits into two
> candidates, and **neither delivers the Clay atom** `f₁²>0` at the physical ratio `3.9√σ`. The **dual Higgs**
> `m_H` is discrete (atom-*like*) but its ratio is short (`1.36–3.72√σ`) and its `tr F²` overlap is **free** in
> the minimal dual Abelian–Higgs (AHM). The **vortex ring**, which does reach `≈3.9√σ`, has a **continuous
> radius** ⇒ continuous spectrum ⇒ smooth spectral density ⇒ **`f₁² = 0` (continuum edge, not an atom)**.
> A KL atom would require a **stable discrete ring** (an energy minimum at finite `R*>ξ`), which the minimal
> dual GL does not supply. So the dual picture's "renormalon-free" advantage survives for the **mass only**;
> the **atom is `H1` in the dual coordinate** — confirming and *explaining* `PNT-DUAL-CONDENSATE` structurally.

---

## §0 · Labels

**[V-num]** new numeric · **[C]** conditional · **[O]** exact missing step · **[NG]** no-go. Witness:
`w7/construct_dual_lattice.py` (runs clean, numpy only, `py -3.12 -I`).

---

## §1 · The dual-superconductor construction

Dual Abelian–Higgs (London limit), from `ROOF-DUALSUPERCONDUCTOR` / `PNT-DUAL-CONDENSATE`:

```
m_V = g v            (dual photon ≡ 1⁻⁻),
m_H = √(2λ) v        (dual Higgs),
κ_GL = m_H / m_V,
σ_ring = π v² ln κ_GL    (vortex line tension),
ξ = 1/m_H                 (core radius).
```

The two `0⁺⁺`-channel objects in the dual sector are **(a) the dual Higgs** (a genuine scalar particle of
mass `m_H`, discrete) and **(b) the vortex ring** (a closed dual flux loop = closed chromo-/dual-electric
tube). The lattice evidence that the physical `0⁺⁺` sits at `m_{0⁺⁺}/√σ ≈ 3.9` points at **(b)** (the ring),
because the dual-Higgs ratio tops out at `≈3.72√σ` (at `κ_GL=10`, witness §5) and the ring is the natural
heavy `0⁺⁺`.

---

## §2 · The ring energy is monotone ⇒ continuum (not an atom)

Standard vortex-ring energy (thin core, superfluid/dual-superconductor form):

```
E(R) = A·R·( ln(8R/ξ) − 7/4 ),      dE/dR = A·( ln(8R/ξ) − 3/4 ).
```

For every physical radius `R ≥ ξ`, `8R/ξ ≥ 8 ⇒ ln ≥ 2.079 > 0.75 ⇒ dE/dR > 0` (witness §1:
`min dE/dR = 1.33 > 0` on `[ξ, 40ξ]`). So:

- **`E(R)` is strictly increasing** on the physical domain;
- the spectrum `{E(R) : R > ξ}` is a **continuum** `[E_min, ∞)`, `E_min = E(ξ) > 0`;
- **there is no isolated level** — no radius is dynamically preferred.

A **continuum** spectral measure has a smooth density and **no δ-atom**. The `tr F²` overlap
`w(R) = ⟨0|tr F²|ring(R)⟩` (supported on the flux core, largest for small rings) gives the continuum density

```
ρ_trF²(s) = w(R(s))² / E′(R(s)),        R(s) = E⁻¹(s),
```

which is finite and smooth (witness §2: density range `[4.5×10⁻⁷, 0.75]`, no δ). The KL atom is the **delta
weight** at the lightest mass; the spectral mass in a **shrinking** window around `E_min` **→ 0**
(witness §3: `6.8×10⁻² → 0 → 0 → 0` as `ε = 10⁻¹…10⁻⁴`). Hence

> **`f₁² = ν({m₁²}) = 0`.** The vortex ring is a **continuum edge, not a KL atom.**

---

## §3 · Why there is no atom — and what *would* give one

An atom requires an **isolated** energy level. A ring can be an atom **only** if its energy has a
**minimum at finite `R*>ξ`** (a stable discrete ring). The minimal dual GL has **no trapping term** that
does this: line tension alone gives `E ∝ R ln R`, monotone. Adding a schematic curvature/trapping term
`+Bξ²/R` *does* create a minimum (`R*≈1.36ξ`, witness §4, slope ≈ 0) — **that** would be a δ-atom — but
**that term is not in the minimal dual action**.

A finite-volume torus quantizes `R` (spacing `~1/V`), giving discrete levels that **accumulate into the
continuum as `V→∞`** — a finite-volume artifact, not a genuine infinite-volume atom. (This is the same
`Δ_min( finite-V ) → 0` phenomenon as `SLOT-BOOTSTRAP` / `PNT-OPE-RHO`.)

⇒ To make the ring carry a real atom, one needs either a **stabilizing mechanism** (a discrete ring
`R*` selected by physics the minimal AHM lacks) or a **non-minimal `tr F²` coupling** that fixes the
overlap — both are extra input, i.e. `H1`.

---

## §4 · Dual-GL scaling of the continuum edge

Leading `2πRσ` at `R=ξ` (witness §5):

| `κ_GL` | `σ_ring` | `E_min/√σ` |
|---|---|---|
| 2 | 2.18 | 4.64 |
| 5 | 5.06 | 2.83 |
| 10 | 7.23 | 1.69 |

So the ring continuum **starts** at `O(1–2)√σ` — a few `√σ` up. The heavy `≈3.9√σ` `0⁺⁺` (if a ring) is
simply **further up the same continuum**, not a separate isolated state. The continuum edge is the light
`0⁺⁺` scale; there is no atom at `3.9`.

---

## §5 · The two dual candidates, and why neither is the Clay atom

| dual object | discrete? | ratio `m/√σ` | `tr F²` overlap `f₁²` |
|---|---|---|---|
| **dual Higgs** `m_H` | **yes** (scalar particle) | `1.36–3.72` (**short** of 3.9) | **free** in minimal AHM (`0≤f₁²<∞`) |
| **vortex ring** | **no** (continuum) | reaches `≈3.9` | **`f₁²=0`** (continuum edge) |

The dual Higgs is the only **discrete** state, but (i) its ratio is short and (ii) the minimal AHM fixes
*masses*, not the `tr F²` overlap — the pole `f₁² δ(s−m_H²)` needs a **non-minimal** `tr F²↔|φ|²` coupling
left free by the minimal action (`PNT-DUAL-CONDENSATE`). The vortex ring reaches the ratio but is a
continuum. **So the dual superconductor gives a renormalon-free MASS but no ATOM** — the atom is `H1`.

---

## §6 · Honest bottom line

- **What was constructed:** the vortex-ring `0⁺⁺` concretely — its energy `E(R)=AR(ln(8R/ξ)−7/4)`, its
  monotone slope, its `tr F²` continuum spectral density, and its continuum edge at `O(1–2)√σ`.
- **The answer to the key question:** the vortex ring carries **`f₁²=0` — a continuum edge, NOT a KL atom**.
  The continuous radius is the obstruction; an atom needs a stable discrete ring the minimal dual GL lacks.
- **What this explains structurally:** the dual picture's renormalon-free advantage is **mass-only**; the
  Clay atom `f₁²>0` (with isolation) is **not delivered** — it is `H1` in the dual coordinate, now shown to
  be a *continuum-vs-atom* obstruction rather than a renormalization one.
- **Not a proof; Clay 0%; no QED.** Labels [V-num]/[C]/[O]/[NG].

*Refs:* `PNT-DUAL-CONDENSATE.md`, `ROOF-DUALSUPERCONDUCTOR.md`, `PNT-LARGEN-DUAL-COMBO.md`,
`SLOT-BOOTSTRAP.md`, `PNT-OPE-RHO.md` (usc/lsc floor-blindness = the same finite-data topology).
Witness: `w7/construct_dual_lattice.py`.
