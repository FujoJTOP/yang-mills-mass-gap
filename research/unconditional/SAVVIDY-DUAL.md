# SAVVIDY-DUAL — the Savvidy negative mode in the dual (monopole) language: the condensate gives a PHYSICAL IR cutoff on the softening, but sub-criticality is `H1` in dual coordinates, **not automatic**

> **VERDICT (no QED).** In the dual ('t Hooft–Mandelstam) language the Savvidy negative mode is the
> **would-be longitudinal mode of the dual photon**, eaten by the monopole condensate `v` (dual-photon mass
> `m_V = g̃ v`). The condensate therefore **does** regulate the softening `Var(φ²)` that the electric
> calculation left cutoff only by the quartic `λ₄`: below `m_V` the dual photon is gapped (Meissner), so the
> accumulated softening is **IR-finite** and scales as `1/m_V² ∝ 1/v²` (witness), not the massless
> `1/k²` runaway. **This is a genuine structural improvement.** But the cutoff scale **is** the condensate:
> sub-criticality `S_dual < κ_bare ⟺ m_V > m_V* ⟺ v > v*`, and deriving `v` (= confinement) is exactly
> `H1`. **The dual picture relocates the sub-criticality into the condensate; it does not make it automatic.**
> Same wall, dual coordinates. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The dual dictionary (electric ↔ dual)

The electric statement (`CONSTRUCT-SAVVIDY-CELL`): on a constant chromomagnetic background `B`, the gluon's
lowest Landau level has eigenvalue `−gB` (Savvidy 1977, *Phys. Lett. B* **71**, 133) — an instability at
`|p| ~ √(gB)`. Integrating it out softens the coarse-mode curvature by `Γ″(0) = κ_bare − g² Var(φ²)`.

The dual statement: in the 't Hooft–Mandelstam dual superconductor (`'t Hooft 1978`; Mandelstam,
*Phys. Rep.* **23** (1976) 245), the chromomagnetic field is confined into **Abrikosov vortices**
(Nielsen–Olesen, *Nucl. Phys. B* **61** (1973) 45) — the **dual Meissner effect**. Monopoles condense,
`⟨φ⟩ = v`, giving

  dual photon mass `m_V = g̃ v`,  dual Higgs mass `m_H = √(2λ) v`,  GL parameter `κ_GL = m_H/m_V`.

**Dictionary `[T]/[C]` (standard dual-superconductor lore, model-level):**

| electric | dual |
|---|---|
| Savvidy negative mode (gluon LLL, `−gB`, `|p|~√(gB)`) | would-be **longitudinal mode of the dual photon** |
| the instability (magnetic vacuum is unstable) | the **dual-Meissner** tendency to form Abrikosov vortices |
| the soft direction `φ` (negative curvature) | the mode **eaten by the Higgs mechanism**; its mass is `m_V = g̃v` |

So the "soft mode" whose fluctuations the electric calculation integrates out is, in the dual language, the
mode that the monopole condensate **Higgses**: it acquires the mass `m_V = g̃v`.

---

## §2 · The softening is IR-finite in the dual language — the condensate bounds `Var(φ²)` `[V-num]`

The accumulated softening is `S = ∫ dln k · g² Var(φ²) · prop(k)`. Its `prop(k)` is the two-point function
of the soft mode:
- **electric (massless)**: `prop(k) = 1/k²` ⇒ `S` is **IR-divergent** (`∫ dln k /k² = ∫ dk/k³ → ∞`);
- **dual (massive `m_V`)**: `prop(k) = 1/(k² + m_V²)` ⇒ `S` **saturates**.

Witness `w7/savvidy_dual_witness.py` (exact quadrature):

```
P1 ELECTRIC   kIR=1e-2,1e-3,1e-4,1e-5 -> S = 5e3, 5e5, 5e7, 5e9   (IR-divergent)
P2 DUAL       m_V=0.05..5.0           -> S*m_V^2 -> 9.19 ≈ ln(kUV/kIR)=9.21
P3 threshold  kappa_bare=0.5          -> m_V* = 4.286 ; sub-critical iff m_V > m_V*
```

So `S_dual ∝ 1/m_V² ∝ 1/v²`: **the monopole condensate provides a physical IR cutoff on the softening.**
This is a **genuine improvement** over the electric language, where the only cutoff was the UV quartic `λ₄`
(a microscopic, uncontrolled parameter). In the dual language the cutoff is the confinement scale `m_V` —
physical and IR.

---

## §3 · But sub-criticality is NOT automatic — it is `H1` in dual coordinates `[O]`

Sub-criticality is `S_dual < κ_bare`, i.e. (witness P3)

  `m_V > m_V* = 4.29` (for `κ_bare = 0.5`)  ⟺  `v > v* = m_V*/g̃`.

The cutoff scale `m_V = g̃ v` **is** the monopole condensate `v` — whose derivation is exactly the
dual-superconductor statement of confinement, `H1` (`PNT-DUAL-CONDENSATE`, `KAPPA-DUAL*`,
`KAPPA-HOLO-COMPUTE`: the dual route is *mass-only*, and `v` / `κ_GL` = the confinement scale = `H1`).
Thus the dual picture converts

  *(need uniform sub-criticality of the softening)*  ⟶  *(is the monopole condensate `v` large enough?)*

which is the **same `H1`**. The dual condensate **bounds** the softening (finite, `1/v²`) but does **not
guarantee** it beats `κ_bare`: that needs `v > v*`, i.e. the condensate computed — `H1`.

**Not automatic. `[NG]` for "the dual condensate makes sub-criticality automatic."**

---

## §4 · Verdict

- **`[T]/[C]`** The dual dictionary is standard: Savvidy negative mode ↔ Higgsed dual-photon longitudinal
  mode; monopole condensate `v` gives `m_V = g̃v`; the instability ↔ dual Meissner / vortex formation.
- **`[V-num]`** The condensate **regulates** the softening: `S_dual ∝ 1/m_V²`, IR-finite, versus the
  massless `1/k²` runaway. A physical IR cutoff replaces the uncontrolled UV quartic. **Real structural
  improvement** — the softening is tied to the confinement scale, not to a UV parameter.
- **`[O]`** Sub-criticality `S_dual < κ_bare ⟺ v > v*` = `H1` in dual coordinates. **Not automatic.**
- **`[NG]`** "The dual condensate makes sub-criticality automatic" is false.

**Bottom line.** The dual language does **not** break the wall; it **makes the wall physical**: the
softening's cutoff is `v`, the condensate, and asking whether it is sub-critical is asking for `v` — the
same `H1` as `PNT-DUAL-CONDENSATE` / `KAPPA-DUAL` / `KAPPA-HOLO`. It is the *cleanest* form of those prior
results: **the dual condensate bounds the softening (`[T]/[V-num]`), and sub-criticality is `v > v*` (`[O]`
= `H1`).** Clay 0%; no QED.

*Refs:* `CONSTRUCT-SAVVIDY-CELL.md` (the electric softened-inequality result this dualises);
`PNT-DUAL-CONDENSATE.md` (dual GL: mass-only, atom = H1); `KAPPA-DUAL*`, `KAPPA-HOLO-COMPUTE.md`
(dual/holographic: the IR scale = H1). Witness: `w7/savvidy_dual_witness.py`.
