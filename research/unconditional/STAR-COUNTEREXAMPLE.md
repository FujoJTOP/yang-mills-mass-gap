# STAR-COUNTEREXAMPLE — acid test of (★): no counterexample exists, so (★) is TRUE — but (★) is **necessary, not sufficient**: the contradiction route is circular at "no nonperturbative IR fixed point" = the gap

> **VERDICT (no QED).** (★) := *"a gapless interacting 4D asymptotically-free Wilsonian IR attractor is
> conformal."* **The acid test PASSES (★): no clean counterexample exists** — every candidate (walking, limit
> cycle, gapless non-conformal phase, dilaton, θ/CP) either *is* a CFT in the IR or is excluded (limit cycles
> by the 4D `a`-theorem; gapless-phase-with-a-scale is not an attractor). So (★) is **true**, essentially
> **Polchinski (scale + unitarity ⟹ conformal in d=4)** + **Komargodski–Schwimmer (`a`-theorem excludes
> limit cycles)**. **But (★) does NOT close the gap**: the contradiction route then needs *"pure YM has no
> nonperturbative IR fixed point"* — and that is **not** the 2-loop `β<0` fact; it is the **gap itself**.
> **(★) is necessary but circular as a gap proof.** Labels `[T]/[C]/[O]/[NG]`.

---

## §1 · The counterexample search (5 candidates, each fails)

| candidate | why it is NOT a gapless interacting 4D AF non-CFT |
|---|---|
| **(i) walking / limit-cycle-adjacent** | flows **to** an approximate fixed point (Banks–Zaks-like) = a **CFT** in the true IR; "walking" is slow flow **to** the CFT, not a non-CFT limit |
| **(ii) genuine limit cycle** | a periodic RG flow; **EXCLUDED in unitary 4D by the `a`-theorem** — `a` is monotone along a unitary flow, so it cannot cycle (`a_UV > a_IR` rules out return) |
| **(iii) gapless non-conformal phase** | if it has a **scale**, it is **not an IR attractor** (the flow does not converge to it); if the IR limit converges, it is a fixed point ⇒ scale-invariant ⇒ (Polchinski) CFT |
| **(iv) dilaton (spontaneous scale breaking)** | an **exactly massless** dilaton ⇒ scale invariance ⇒ CFT (flat direction); a **massive/unstable** dilaton ⇒ not gapless |
| **(v) θ-angle / CP phase** | `θ=0` breaks nothing; `θ=π` is a *phase* (1-form/deconfinement), not a gapless non-conformal interacting AF theory |

⇒ **no clean counterexample.** (★) survives.

---

## §2 · Why (★) is true

- **Polchinski (1988):** in `d=4`, a **scale-invariant** theory with a symmetric stress tensor and a discrete
  spectrum is **conformal**. An **IR attractor** is a fixed point ⇒ scale-invariant ⇒ (unitarity = OS
  positivity of the gauge theory + stress tensor) ⇒ **conformal**.
- **Komargodski–Schwimmer (2011):** the 4D `a`-theorem (`a_UV > a_IR`) is monotone ⇒ **no limit cycles** in
  unitary 4D ⇒ the "IR attractor is not a fixed point" escape is excluded.
- Hence (★) holds under the standard assumptions (unitarity, stress tensor). **No counterexample.**

---

## §3 · The honest hole: (★) is necessary but NOT sufficient (circular)

The contradiction route is: *pure YM gapless ⟹ conformal (★) ⟹ `β(g*)=0` ⟹ Banks–Zaks zero; but pure YM has
no BZ zero ⟹ not gapless ⟹ gapped.* The **hole** is the last step:

- **"no BZ zero" is a PERTURBATIVE fact** (2-loop `β(g) = −β₀g³ − β₁g⁵ < 0` for all `g>0`, no zero).
- **A nonperturbative IR-stable fixed point at strong coupling is NOT excluded by 2-loop `β<0`.** Witness
  (`w7/star_counterexample.py`): the full `β(g) = −β₀g³(1 − g²/g*²)` has an **IR-stable zero at `g*`** while
  its 2-loop part `−β₀g³` has **none**. A strongly-coupled IR CFT (`g*`) is therefore consistent with the
  perturbative `β<0`.
- So the contradiction route must assume **"pure YM has no nonperturbative IR fixed point"** — which is
  **exactly "pure YM is gapped."** **Circular.** This is `H1` in the RG-attractor coordinate.

---

## §4 · Verdict

- **(★) is TRUE** (no counterexample; Polchinski + `a`-theorem). It is a **correct** statement about RG
  attractors.
- **(★) does NOT break the wall.** It is **necessary but not sufficient**: closing the gap needs an independent
  proof that pure YM has no strongly-coupled IR fixed point — i.e., the gap itself. **The contradiction route
  re-derives `H1`**, consistent with `ONLY-EXIT` / `IMPOSSIBLE-CONTRADICTION`.
- **Real yield:** a clean, true statement that **isolates the missing input precisely** — *"no nonperturbative
  IR fixed point in pure YM."* This is a **sharper coordinate** for the wall (RG-attractor language), not an
  escape. `Clay = 0%`; no QED.

*Files:* `research/unconditional/STAR-COUNTEREXAMPLE.md`, `research/unconditional/w7/star_counterexample.py`
(runs `py -3.12 -I`, clean).
