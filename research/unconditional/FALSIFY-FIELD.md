# FALSIFY-FIELD — the field-level route does NOT falsify "only exit is YM structure or H1": functionals of the gauge field reduce to `ν` (blind, OPE-RHO) or to a different channel (⇒ dynamical bridge = H1); "YM structure" and "H1" are the SAME exit

> **VERDICT: the falsification FAILS — the only-exit claim survives, sharpened.** A functional of the gauge
> measure `μ` that uses data *beyond* the 0⁺⁺ spectral measure `ν` does not give a third exit. Either it is a
> functional of the 0⁺⁺ correlator (= a functional of `ν`) → **blind (OPE-RHO, upper bounds only)**; or it
> lives in a *different* channel (1-form / topological) → by the gauge-invariant Hilbert-space decomposition it
> is **Schur-orthogonal to the 0⁺⁺ atom** ⇒ to connect it to `f₁²` needs a **dynamical bridge** = `H1`. The
> only genuine non-integral field object (`TOOL-HARDY`'s `N(ε)=∫|G|²`, quadratic/singular kernel, escapes
> OPE-RHO) detects atom *presence* but its weight needs the pole location (=gap) and `ε→0`/`L→∞` (=H1).
> Conclusion: **"use YM structure" is not a loophole around `H1` — proving the needed field→atom implication
> IS `H1`.** They are one exit, not two. `Clay 0%.`

---

## §0 · Setup and labels

`f₁² = |⟨0|tr F²|G⟩|²`, `G` = lightest `0⁺⁺` eigenstate; `ν` = the 0⁺⁺ KL spectral measure. `μ` = the gauge
field measure. **[T]/[C]/[O]/[NG]**. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The two cases — why no field functional forces `f₁²` without H1

A field functional `F(μ)` forces `f₁²>0` only if `F ∈ S ⟹ f₁²>0` for some checkable set `S`. Two cases:

**(A) `F` is a functional of the 0⁺⁺ correlator = of `ν`.** Then `F` is a finite integral/Laplace/moment of
`ν` (e.g. the Lüscher-flow energy `E(t)=∫ s e^{−2st}dν`, moments, OPE coefficients). These are bounded-kernel
integrals ⇒ **OPE-RHO applies verbatim**: the atom is usc, the data lsc ⇒ blind; at best an **upper** bound.
Witness §2: an atomless bump matches every finite moment (and `E(t)`) → no lower bound.

**(B) `F` uses field data not determined by `ν`** (Wilson loops, `⟨tr Fⁿ⟩` `n≠2`, instanton/vortex density,
topology). By gauge-invariant Hilbert-space decomposition these live in **other sectors** (1-form center /
topological), **Schur-orthogonal** to the 0⁺⁺ scalar `tr F²` (`[T]`, `ROOF-CONFINEMENT-GAP` kinematic
orthogonality). So `F` carries **no kinematic information** about `f₁²` (witness §3: block-diagonal toy —
varying the 1-form sector leaves the 0⁺⁺ atom spectrum unchanged). To connect `F` to `f₁²` requires a
**dynamical** statement about the specific YM measure (confinement ⇒ 0⁺⁺ atom; `n`-point structure ⇒ atom).
Any such statement is a theorem about `μ` — **that is `H1`** (`m_{0⁺⁺}≥4√σ` is its only proven instance,
`β<1/48`, `[C]`; its extension is `H1`).

⇒ **In both cases, forcing `f₁²` from field data either fails (blind) or reduces to a dynamical YM theorem =
`H1`. There is no third exit.** "YM structure" and "H1" are the same exit.

---

## §2 · Finite-L triviality vs the continuum object

At finite `L` the transfer matrix is compact ⇒ discrete spectrum ⇒ the lightest 0⁺⁺ level **is** an atom and
`f₁²>0` automatically (witness §1: any nonzero overlap gives `f₁²>0`). The Clay content is the **`a→0`,
`L→∞` *survival* of an isolated atom** (a uniformity) — the atom is observed to *dissolve into the cutoff*
as `β` decreases. That survival is `H1`, not a finite-L field functional.

## §3 · The one genuine non-integral field object (HARDY) is still H1 for the weight

`N(ε)=∫_{ε<|x−m₁²|<δ}|G(x+i0)|²dx` (`G` the Stieltjes transform) is **quadratic with a singular kernel** ⇒
outside the bounded-kernel integral class ⇒ OPE-RHO does not apply. Witness §4 (exact Hilbert transforms):
**atom → `N·ε → const>0`; smooth continuum → `N·ε → 0`** — a robust presence detector. **But** the weight
needs the pole **location** `m₁²` (locating = the gap), the `ε→0` limit (=A3), and the `L→∞` survival of the
isolated pole (=isolation = `H1`). And at finite `L` every level is a pole ⇒ the divergence is uniform/trivial.

---

## §4 · Honest status

- **Falsification FAILS.** No field-functional (Wilson, `⟨tr Fⁿ⟩`, instanton, gradient flow, entanglement)
  forces `f₁²` from finite data without `H1`. The only-exit claim is confirmed.
- **Sharpening:** "YM structure" is *not* a separate exit — it is the concrete realization of `H1`; proving the
  field→atom implication *is* `H1`. One exit.
- **What the field route *does* give:** (i) a real **presence** detector that escapes OPE-RHO (`N(ε)`, robust
  to continuum, atom vs smooth cut); (ii) the finite-L atom is trivially `>0`, so the whole problem is the
  `a→0,L→∞` survival + isolation = `H1`.
- **Clay 0%; no QED.**

*Witness:* `w7/falsify_field.py` (numpy, runs clean). *Refs:* `TOOL-SYNTHESIS.md` (HARDY, OPE-RHO class),
`ROOF-CONFINEMENT-GAP.md` (kinematic orthogonality), `PNT-OPE-RHO.md` (atom usc / data lsc),
`SAVVIDY-WI-SYNTHESIS.md`.
