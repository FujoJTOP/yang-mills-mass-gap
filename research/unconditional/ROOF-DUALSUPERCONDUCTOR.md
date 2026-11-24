# ROOF-DUALSUPERCONDUCTOR — dual Ginzburg–Landau gives a **renormalon-free** 0⁺⁺ mass, but it certifies a *condensate gap*, not the KL *atom*, and it *assumes* monopole condensation = the wall

> **VERDICT (no QED).** Dual superconductivity / dual Ginzburg–Landau ('t Hooft–Mandelstam) **does** give
> a 0⁺⁺ mass that is **structurally renormalon-free**: `m_{0⁺⁺}=m_H=√(2λ)·v` is a finite **algebraic
> mean-field gap** (no loops, no factorial growth), and it is built from the **dimension-2** monopole
> condensate `v` (natural Borel pole at `u=1`, UV-finite on the dual lattice), **not** the dimension-4
> `⟨tr F²⟩` with its `u=2` ambiguity. So it **sidesteps** the exact wall that killed `ROOF-CONDENSATE-FLOOR`.
> **But** (i) it *assumes* monopole condensation (deriving it from QCD **is** confinement = `H1`); (ii) it
> certifies a **condensate gap**, not the **KL atom** (`f₁²>0`) — the same atom/gap distinction as
> everywhere; (iii) the naive "0⁺⁺ = dual Higgs" identification gives `m_{0⁺⁺}/√σ ≲ 3.7` at `κ_GL≲10`,
> short of the lattice `≈3.9`, i.e. the true 0⁺⁺ is *heavier than the plain dual Higgs* (composite/ring).
> **A renormalon-free BASE, not a key. Clay 0%.**

---

## §0 · Labels, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** witness.
**Acid**: the dual-GL gap returns `0` for the free massless gauge field (`v=0 ⇒ m_H=m_V=0`) — *acid-clean
there*. **It does NOT pass the discriminating `:φ²:` test**: dual GL is a gauge-only construction (a
`!φ²:` composite has no gauge field, so dual GL does not apply), and it predicts a mass for *any*
condensing gauge theory without knowing the atom — see §4.

---

## §1 · The construction (existing theory)

Dual superconductor picture of confinement ('t Hooft–Mandelstam): the QCD vacuum is a **dual
superconductor**; magnetic monopoles condense; the electric flux of a `q\bar q` pair is squeezed into a
**dual Abrikosov vortex** (the flux tube) with tension `σ`. Effective **dual Ginzburg–Landau** model

```
L = −¼ F̃_{μν}F̃^{μν} + |D_μ φ|² + (λ/2)(|φ|² − v²)² ,   D = ∂ − i g̃ Ã ,
```
with `φ` the dual (magnetic) Higgs = the monopole condensate, `Ã` the dual photon. Standard spectrum:

| mode | mass | QCD reading |
|---|---|---|
| dual photon | `m_V = g̃ v` | `1⁻⁻` (vector glueball) |
| dual Higgs | `m_H = √(2λ) v` | `0⁺⁺` (scalar glueball) |
| vortex (flux tube) | tension `σ` | string tension |

GL parameter `κ_GL := m_H/m_V = √(2λ)/g̃`; vortex tension (London limit, type-II, large `κ_GL`)
`σ ≃ π v² ln κ_GL` (BPS anchor at `κ_GL=1/√2`: `σ = π v²` exactly).

---

## §2 · The 0⁺⁺ mass is a finite algebraic gap

The radial mode of the dual condensate is the `0⁺⁺` glueball candidate. Its mass is the curvature of the
GL potential at the minimum: `V(v)= (λ/2)(v²−v₀²)²`, `dV/d|φ|=0 ⇒ |φ|=v₀`, so

```
m_{0⁺⁺} = m_H = √(V''(v₀)) = √(2λ)·v₀ .                                              [T]
```
This is a **tree/mean-field, algebraic, definite** number — **no loop expansion, no factorial growth**
(witness §2). Contrast with the dimension-4 gluon condensate `⟨tr F²⟩`, whose perturbative coefficients
`a_n ~ 2ⁿ Γ(n+1+2b)` grow **factorially** (`a_n/a_{n−1} ~ 2(n+2b) → ∞`), giving a Borel pole at **`u=2`**
and an ambiguity `~Λ⁴` (witness §3; `ROOF-CONDENSATE-FLOOR`).

---

## §3 · The key question: renormalon-free? **Yes, structurally.**

- The dual-GL mass is **algebraic** (§2) ⇒ **no factorial series** ⇒ **no `u=2` IR renormalon**.
- It is built from `v`, the **monopole condensate**, a **dimension-2** object (natural Borel pole at
  `u = 2/2 = 1`), and on the **dual lattice** `v` is a **UV-finite order parameter** (monopole density) —
  a *different* operator from the dimension-4 `⟨tr F²⟩`.
- Hence the dual GL mass **does not reproduce** the `u=2` ambiguity that blocks the gluon-condensate floor.
  **This is a genuine structural advantage**: the strut, expressed in dual variables, is a *renormalon-free*
  quantity.

---

## §4 · What it costs (why it is a base, not a key)

1. **It assumes monopole condensation.** Dual GL is *effective*; deriving it from the QCD measure by
   monopole condensation **is** the confinement mechanism = the wall (`H1`). The dual parameters `(v, g̃, λ)`
   are **fitted to lattice data**, not computed from first principles. So the renormalon-freedom is bought
   by relocating the problem to "does the dual condensate exist", not by solving it. `[O]`
2. **Condensate gap ≠ KL atom.** `m_H` is the mass of a *background-condensate* excitation; it does **not**
   certify the Källén–Lehmann atom weight `f₁² = |⟨0|tr F²|G⟩|² > 0`. A gauge theory can have a dual
   condensate gap while the `tr F²` spectral atom is thin or absent — the same **atom/gap** distinction the
   whole campaign keeps hitting. So dual GL gives a **gap** but not the **atom** (`[NG]` for the atom). This
   is why it fails the *discriminating* `:φ²:` acid: it is blind to the atom.
3. **The ratio is short at moderate `κ_GL`.** With `g̃=1, v=1`: `m_H/√σ = κ_GL/√(π ln κ_GL)` (London),
   giving `1.36 (κ=2), 2.22 (κ=5), 3.72 (κ=10)` (witness §1). The lattice `m_{0⁺⁺}/√σ ≈ 3.9` needs
   `κ_GL ≈ 11⁺`, whereas dual-GL fits give `κ_GL ~ O(1)`. Reading: the **plain dual Higgs is lighter than the
   true 0⁺⁺**; the 0⁺⁺ glueball is a more **composite** object (vortex-ring / two-vortex excitation), whose
   mass is not just `m_H`. `[V-num]`

---

## §5 · Relation to `⟨tr F²⟩` and `σ`

- `σ` (vortex tension) and `m_H` share the **same condensate `v`**; in London, `m_H/√σ = κ_GL/√(π ln κ_GL)`
  — an O(1)–O(4) function of `κ_GL`, no renormalon.
- `⟨tr F²⟩` (the u=2-ambiguous object) and the dual condensate `v` are **different operators** (dim-4 vs
  dim-2, electric vs dual). Dual GL does **not** give a scheme-independent floor of `⟨tr F²⟩`; it replaces
  that demand by **"`v>0`"** — a renormalon-free statement but a *different* (and unproven) one.

---

## §6 · Honest bottom line

- **Renormalon-free: YES** — dual GL's `0⁺⁺` mass is a finite algebraic gap from a dim-2 condensate, free of
  the `u=2` ambiguity that blocks `⟨tr F²⟩`.
- **A key: NO** — it *assumes* monopole condensation (the wall), certifies a *condensate gap* not the KL
  *atom*, and its plain dual-Higgs mass is short of the lattice 0⁺⁺ ratio unless `κ_GL ≫ 1`.
- **Net:** dual GL is the **cleanest renormalon-free *base*** for the strut found so far — it shows the strut
  *can* be expressed without the `u=2` renormalon, at the price of the (unproven) monopole condensate. The
  remaining work is unchanged in name (`H1`/monopole condensation) but **restated in a renormalon-free
  dialect**: "does the dual condensate `v` exist, provably, in the continuum?" — a cleaner statement of the
  same wall.
- **Clay 0%; no QED.**

*Refs:* `ROOF-CONDENSATE-FLOOR.md`, `ROOF-SYNTHESIS.md`, `BASE-FIELDMEASURE.md`; witness
`w7/roof_dualsuperconductor.py` (runs clean, numpy+math, `py -3.12 -I`). 't Hooft 1978; Mandelstam 1976;
Nielsen–Olesen 1973; dual-GL lattice fits (Suzuki et al.). Labels `[T]/[C]/[O]/[NG]/[V-num]`.
