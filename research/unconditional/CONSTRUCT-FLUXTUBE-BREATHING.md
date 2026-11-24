# CONSTRUCT-FLUXTUBE-BREATHING — the 0⁺⁺ as a closed flux loop's radial (breathing) mode: closed-form `m_breath = 2√(σ/κ)`, an *excited* 0⁺⁺ that overshoots the lightest glueball; atom character hinges on `κ = H1`

> **VERDICT.** Model the neutral 0⁺⁺ glueball as a closed flux loop (rigid string, `E(R)=2πσR+πκ/R`,
> MASTERFIELD). Quantize the transverse radial coordinate `R(t)=R*+δR cos ωt`. This is a harmonic
> oscillator with inertia `M=2πσR*` (tube mass) and stiffness `E''(R*)=2πκ/R*³`, giving the **closed form**
> `m_breath = √(E''/M) = 2√(σ/κ)`. At the MASTERFIELD fit `κ=0.193`: `m_ground=2π√(2κ)=3.90√σ` (=lattice)
> and `m_breath = 4.55√σ` — the breathing mode is an **excited 0⁺⁺**, overshooting the lightest glueball by
> **17%**. It sits between the one-phonon (`1/R*=3.2`) and two-phonon (`2/R*=6.4`) thresholds. **Atom vs
> continuum is decided by the rigid-string stabilization itself (`κ>0` = `H1`), not by the breathing
> computation.** Labels [C]/[V-num]/[O]; no QED.

---

## §1 · Setup and the closed form

Rigid-string closed loop (MASTERFIELD): `E(R)=2πσR+πκ/R`, interior minimum `R*=√(κ/2σ)`,
ground-state mass `m_ground = E(R*) = 2π√(2σκ)`. The **breathing mode** is the radial collective
coordinate `R(t)`:
- inertia (tube mass, `E=mc²`): `M_eff = σ·(2πR*)`;
- stiffness: `E''(R*) = d²/dR²(2πσR+πκ/R)|_{R*} = 2πκ/R*³`;
- frequency: `ω² = E''/M_eff = κ/(σR*⁴) = 4σ/κ` ⇒ **`m_breath = ω = 2√(σ/κ)`** (σ=1 in the witness).

Ratio: `m_breath/m_ground = 1/(π√(2κ))`.

## §2 · Numbers ([V-num], witness table)

| κ | R* | m_ground/√σ | m_breath/√σ | ratio | 1/R* | 2/R* |
|---|---|---|---|---|---|---|
| 0.138 | 0.263 | 3.30 | 5.38 | 1.63 | 3.81 | 7.61 |
| 0.193 | 0.311 | **3.90** | **4.55** | 1.17 | 3.22 | 6.44 |
| 0.250 | 0.354 | 4.44 | 4.00 | 0.90 | 2.83 | 5.66 |

At `κ=0.193` (MASTERFIELD → `m_ground=3.90`=lattice): `m_breath=4.55√σ`, **overshoots the lightest 0⁺⁺
by 17%** — consistent with an *excited* 0⁺⁺ (e.g. `0⁺⁺*`/`0⁺⁺`≈1.6 on the lattice; the ratio here spans
1.17–1.63 over the κ band). So the breathing-mode identification does **not** give the *lightest* 0⁺⁺.

## §3 · Atom or continuum? (the key question)

On the loop's circle, transverse phonons are quantized at `n/R*`. `m_breath` lies in `(1/R*, 2/R*)`:
- it **cannot** decay into **two** phonons (`m_breath < 2/R*`);
- single-phonon decay (`1/R* < m_breath`) requires a `J^{PC}=0^{++}` phonon combination; the leading
  model does not fix these selection rules.

Within the rigid-string model **both** `m_ground` and `m_breath` are **discrete** levels (the interior
minimum makes the loop a stabilized bound state) ⇒ model-level atoms (`f₁²>0`). **But** this is exactly the
hypothesis that is unproven: in the **minimal dual GL** (no rigidity) `E(R)` is monotonic and the ring is a
**continuum edge**, `f₁²=0` (CONSTRUCT-DUAL-LATTICE). The interior minimum — hence the atom — exists **iff**
the rigidity `κ>0` is physical.

## §4 · Exact failing step [O]

The breathing mode is a **closed-form discrete level in the rigid-string model**, but delivering a physical
atom `f₁²>0` requires (i) the rigid-string stabilization `κ>0` (= `H1` in string coordinates), and (ii) the
nonperturbative `tr F²` overlap + `J^{PC}` selection rules (not fixed by the effective string). Without
them it dissolves into the transverse continuum (`f₁²=0`). Same invariant: **`κ = H1`**.

*Files:* `w7/construct_fluxtube_breathing.py` (runs clean). Labels [C]/[V-num]/[O]; Clay 0%; no QED.
