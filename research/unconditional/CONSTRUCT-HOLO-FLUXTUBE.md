# CONSTRUCT-HOLO-FLUXTUBE — the holographic flux tube gives a 0⁺⁺ **atom** (discrete KK level) with the right O(1) ratio, but by *assuming* the confining bulk geometry (the atom is the dual of `H1`)

> **VERDICT (no QED).** The AdS/CFT Wilson-loop construction (Maldacena 1998; Rey–Yee 1998) in a
> confining bulk does give a genuine **0⁺⁺ atom**: the glueball is a *discrete* normalizable bulk KK level
> (hard-wall zero `j_{2,1}/z0`, or a Witten-cigar KK mode), **not a continuum**. In the **Witten D4 cigar**
> the ratio is a pure O(1) number `m_{0⁺⁺}/√σ ≈ 3.5–4` (Brower–Mathur–Tan, hep-th/0003115), matching the
> lattice `≈3.9`. **But** the atom's discreteness and the confining geometry are *one and the same input*:
> the IR wall / cigar **is** the dual of confinement. Holography therefore **constructs an atom with the
> right shape by postulating the wall** — deriving that wall from Yang–Mills *is* `H1`. And the naive
> hard-wall ratio `j_{2,1}/√λ'` is O(1) only at `λ'≈1.7`, **outside** the controlled SUGRA limit. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 (within the stated model) · **[C]** 条件定理 (model/conjecture) · **[O]** 精确缺失步 ·
**[NG]** no-go · **[V-num]** new numeric. Since AdS/CFT is a conjecture, every holographic statement here
is **[C]** for real Yang–Mills, **[T]** only *within* the bulk model.

---

## §1 · The construction (what the flux tube is)

**String / flux tube.** In a confining holographic background a Wilson loop `W(C)` is computed by a
fundamental string ending on `C` at the boundary (`z→0`) (Maldacena, hep-th/9803002; Rey–Yee,
hep-th/9803001). For a `T×R` rectangle with `T→∞`, the string dips to the IR region; in a **confining**
background (hard wall at `z=z0`, or the Witten D4 cigar at `u=u_Λ`) it sits at the IR tip and gives a
**linear potential** `E(R)=σR`. The flux tube is that string; its **transverse/breathing fluctuations are
the glueball modes**.

**0⁺⁺ glueball.** The `0⁺⁺` is the lightest **bulk scalar** fluctuation (the dilaton in the simplest
identification). Its mass is a KK eigenvalue of the bulk mode equation with an IR boundary condition.

---

## §2 · Hard-wall spectrum (computable, [T] within the model)

Hard-wall AdS₅, dilaton `0⁺⁺` (`ν=2`): the mode equation is the Bessel equation, and the spectrum is
`m_n = j_{2,n}/z0`. The witness gives

```
j_{2,n}   = [5.13562, 8.41724, 11.61984, 14.79595]
m_n/m_1   = [1.000, 1.639, 2.263, 2.881]
```

The tower `1 : 1.64 : 2.26 : 2.88` reproduces the lattice glueball ratios `1 : ~1.6 : ~2.3 : ~2.9`
(Morningstar–Peardon). **The spectrum is discrete (discrete Bessel zeros) ⇒ the `0⁺⁺` is an isolated
level (an ATOM), a normalizable bulk mode ⇔ a CFT state.**

---

## §3 · The ratio `m_{0⁺⁺}/√σ` — hard wall vs Witten cigar

**String tension at the IR wall.** A Nambu–Goto string at `z=z0` has induced metric `(L²/z0²)·diag(−1,1)`, so

```
σ = L² / (2πα' · z0²) ,      λ' := L²/(2πα')
√σ = √λ' / z0 ,              m_{0⁺⁺}/√σ = j_{2,1}/√λ' .
```

The ratio is a pure number set by the geometry. The witness scan (`[V-num]`):

| `λ' = L²/(2πα')` | `m_{0⁺⁺}/√σ = j_{2,1}/√λ'` |
|---|---|
| 1.000 | 5.136 |
| **1.734** | **3.900** |
| 3.000 | 2.965 |
| 10.00 | 1.624 |
| 30.00 | 0.938 |
| 100.0 | 0.514 |

To match the lattice `3.9` one needs `λ' ≈ 1.73` — **small**, *outside* the supergravity regime
(`λ'≫1`). In the controlled regime the naive hard-wall ratio `→0`. **So the hard-wall `m/√σ` is a FIT at
`λ'≈1.7`, not a prediction.** `[O]`

**Witten cigar — the genuine O(1) result `[C]`.** In the Witten D4 model (hep-th/9803131) the geometry is
genuinely confining (a cigar cap at `u=u_Λ`). Both the glueball masses and the string tension are set by
the **same** IR scale `M_KK = (2/3)(u_Λ/R³)^{1/2}`, hence `m_{0⁺⁺} ~ c·M_KK`, `√σ ~ c'·M_KK`, and

```
m_{0⁺⁺}/√σ = c/c' = O(1)   (independent of λ').
```

Brower–Mathur–Tan (hep-th/0003115) find the `0⁺⁺` lightest with `m_{0⁺⁺}/√σ ≈ 3.5–4`, and the ordering
`m(0⁺⁺) < m(2⁺⁺) < m(1⁺⁺)`, **matching the lattice `≈3.9`**. *(⚠ the precise number is often quoted, not
verbatim-verified here — treat `≈3.5–4` as the model's O(1) output.)*

---

## §4 · Atom or continuum? — **ATOM**, but the atom is the dual of `H1`

- **Discrete ⇒ atom.** The hard-wall (Bessel) and cigar (KK) spectra are **discrete**; the `0⁺⁺` is an
  isolated level, a normalizable bulk mode dual to a CFT state with `f₁²>0`. So holography **does produce
  an atom**, unlike the dual-superconductor vortex-ring (which was a continuum edge, `f₁²=0`,
  `CONSTRUCT-DUAL-LATTICE`).
- **But the discreteness is the wall.** The IR wall / cigar — what makes the spectrum discrete and
  confining — **is the dual of confinement**. The atom exists *because* the bulk has an IR tip. Deriving
  that tip from Yang–Mills (i.e. proving the confining geometry is the correct dual of pure YM) **is
  `H1`**. So the holographic atom is **the dual of the gap, imported as a boundary condition**, not derived
  from the microscopic theory.
- **Channel note.** The `0⁺⁺` here is a *bulk scalar* (dilaton), the right channel — consistent with
  `BASE-ORDER-PARAMETER`'s rule that only a bulk, 0-form object can sit in `0⁺⁺`.

---

## §5 · Verdict and what it teaches

- **Positive `[C]`:** the holographic flux tube **constructs a 0⁺⁺ atom** with the right O(1) ratio
  (`≈3.5–4`, lattice `3.9`) and the right tower (`1 : 1.6 : 2.3 : 2.9`). This is the cleanest example so
  far of a *genuine* atom in a controlled model — **and it re-shows that the atom is the right target**
  (mass ≠ atom: the dual-superconductor gave only a mass/continuum; holography gives a genuine discrete
  atom).
- **Exact failing step `[O]`:** the confining IR geometry (wall / cigar) is **assumed**, not derived — it
  is `H1` in the bulk-geometry coordinate. The naive hard-wall ratio is O(1) only at `λ'≈1.7`
  (uncontrolled). Holography is a **fittable model with one input (the IR scale / geometry) = `H1`**, not
  a derivation.
- **Contrast with the flux-tube rigidity route (`CONSTRUCT-MASTERFIELD`):** that route gave
  `m/√σ = 2π√(2κ)` with `κ` (rigidity) missing; this route gives `m/√σ = j_{2,1}/√λ'` (hard wall) or
  O(1) (cigar) with the **IR geometry / `λ'`** missing. **Same structural position: one model input = `H1`,
  just a different coordinate (bulk geometry vs flux-tube rigidity).**

*Labels:* hard-wall Bessel spectrum & `0⁺⁺` atom: `[T]` *within the hard-wall model*; Witten O(1) ratio:
`[C]` (cited, not recomputed); the atom = dual of `H1`: `[O]`; the `λ'` scan: `[V-num]`. **No QED.**

*Refs:* `w7/construct_holo_fluxtube.py`; `ROOF-HOLOGRAPHY.md`, `CONSTRUCT-MASTERFIELD.md`,
`CONSTRUCT-DUAL-LATTICE.md`, `BASE-ORDER-PARAMETER.md`.
*Citations:* Maldacena hep-th/9803002; Rey–Yee hep-th/9803001; Witten hep-th/9803131;
Brower–Mathur–Tan hep-th/0003115; Morningstar–Peardon hep-lat/9901004.
