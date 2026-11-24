# ROOF-HOLOGRAPHY — holographic QCD predicts a definite 0⁺⁺ formula and ratios, but RELOCATES the ambiguity into the IR boundary condition; it does not derive a renormalon-free `m_{0⁺⁺}`

> **VERDICT (no QED).** Holographic-QCD models — hard wall (Polchinski–Strassler), soft wall
> (Karch–Katz–Son–Stephanov), Witten's D4-brane, Sakai–Sugimoto — give **definite formulas** for the
> `0⁺⁺` glueball mass and reproduce the lattice ordering/ratios qualitatively. **But in every such model
> the absolute mass is a dimensionless number times an INPUT IR scale** (`z₀` for the hard wall, the
> dilaton slope `c` for the soft wall, the compactification radius `R` for the D4-brane). So the
> dimensional-transmutation / renormalon ambiguity is **RELOCATED into the IR boundary condition, not
> removed**: holography predicts `m_{0⁺⁺}/√σ` (scale-free ratio) but not `m_{0⁺⁺}` itself. And it is **not
> rigorous**: there is no proven duality between pure 4D YM and these geometries. **It is a model of the
> target (`H1`), not a derivation of it. Clay 0%.**

---

## §0 · Labels, acid test

**[T]** rigorous (within the stated model) · **[C]** conditional / conjectural-duality · **[O]** exact
missing step · **[NG]** no-go for the class · **[V-num]** numeric. Acid: a faithful detector must vanish on
the free massless field; holography is a strong-coupling model (no free-field limit in the same sense), so
the acid test is **structural**, not a numeric free-field check — noted honestly.

---

## §1 · What holography actually gives (definite formulas)

**Hard wall** (Polchinski–Strassler 2001, hep-th/0101126; Boschi-Filho–Braga): AdS₅ cut at `z=z₀`; the
`0⁺⁺` scalar glueball solves `d_z(z^{-3}d_z φ) + m² z^{-5}φ = 0`, `φ = z² J₂(mz)`, IR BC `J₂(mz₀)=0`:

```
m_{0⁺⁺} z₀ = j_{2,1} = 5.1356 ,   m_n/m_1 = 1 : 1.639 : 2.263 : 2.881   [V-num, this round]
```

The ratio `m_{0⁺⁺*}/m_{0⁺⁺} = 1.639` matches the lattice `0⁺⁺*/0⁺⁺ ≈ 1.6` — a genuine success.

**Soft wall** (Karch–Katz–Son–Stephanov 2006, hep-th/0602224): dilaton `Φ = c z²`, linear Regge
`m_n² = 4c(n+1)` ⇒ `m_n/√(4c) = 1 : 1.414 : 1.732 : 2.0` (Regge slope fixed by `c`).

**Witten D4-brane / BMT** (Witten 1998 hep-th/9803131; Brower–Mathur–Tan, *Nucl. Phys.* B587 (2000) 249,
hep-th/0003115): strong-coupling spectrum from the AdS₇×S⁴ black-hole background; the **ordering
`m(0⁺⁺) < m(2⁺⁺) < m(1⁺⁺)`** (dilaton `0⁺⁺` degenerate with tensor `2⁺⁺`), reproducing the lattice
lightest-scalar feature. The lightest `0⁺⁺` is often quoted near `4.07 M_KK` (M_KK = inverse
compactification radius); ⚠ this specific figure was **not verbatim-verified** from the source this round
(the abstract/ordering were), so it is stated as "commonly quoted", per the corpus honesty rule.

---

## §2 · The ratio `m_{0⁺⁺}/√σ`, and consistency with the strong-coupling bridge `m ≥ 4√σ`

Holographic models are **scale-free**: `m_{0⁺⁺} = k/z₀` and `√σ = k′/z₀` for the same `z₀`, so they
predict the **ratio** `m_{0⁺⁺}/√σ = k/k′` (model-dependent; holography lands at `O(3–4)`, lattice
`m(0⁺⁺)/√σ ≈ 3.5–4`, Michael–Teper 1993; Morningstar–Peardon 1999 hep-lat/9901004). **No tension with the
strong-coupling bridge `m_{0⁺⁺} ≥ 4√σ`**: that bound is a `u ≤ e⁻¹` (deep strong-coupling) statement where
`m/√σ = 4√ln(1/u) > 4`; at the crossover the ratio drops toward the physical `≈3.7`, exactly as lattice
and holography both show. Holography is **consistent** with the bridge and with lattice, in the ratio.

---

## §3 · The bulk "strut": the dilaton / IR geometry (an INPUT)

The operator `tr F²` is dual to the **bulk dilaton** (since `1/g² ~ e^{-Φ}`, the trace anomaly lives in
the dilaton sector); the `0⁺⁺` glueball is a dilaton–metric fluctuation. Its mass is set by the **IR
geometry**: the hard-wall position `z₀`, the soft-wall dilaton slope `c`, or the D4 compactification
radius `R`. **The strut, in the holographic picture, is precisely the IR boundary condition / dilaton
profile — and that is what is put in by hand.**

---

## §4 · Renormalon-free? No — relocated, not removed

Because the strut is an **input IR scale**, holography gives:

```
m_{0⁺⁺} = (dimensionless number) × (input IR scale z₀ / √c / 1/R).
```

The **absolute** gap is therefore *not* computed; only **ratios** are predicted. The
dimensional-transmutation / `u=2`-renormalon ambiguity of `ROOF-CONDENSATE-FLOOR` is **relocated into the
choice of `z₀` / `c` / `R`**, not eliminated. Holography is a **strong-coupling model**; its own
systematic errors (`α′`, supergravity truncation, the model-to-YM map) replace the renormalon with a
different uncontrolled input. In the language of the roof: holography **parametrizes** the missing strut
`χ_lattice ≲ a²` by an IR scale and **assumes** the answer, rather than deriving it.

---

## §5 · Conjectural status (why this is not a proof)

AdS/CFT is a **conjecture** even for `N=4` SYM; Witten's D4-brane model and Sakai–Sugimoto are **further
deformations** that share qualitative features of pure YM (confinement, a glueball tower, a mass gap) but
are **not proven dual to 4D pure YM**. No model here is a rigorous statement about the Clay theory. So the
labels are `[C]`/conjectural throughout; nothing here is `[T]` for YM itself.

---

## §6 · Honest bottom line

- **Holography predicts definite formulas and ratios** (`m_{0⁺⁺}z₀ = 5.136`; `m_{0⁺⁺*}/m_{0⁺⁺}=1.64`;
  `m_n² ∝ n`), reproduces the lattice ordering and `m/√σ ≈ 3.5–4`, and is **consistent** with the
  strong-coupling bridge `m ≥ 4√σ`.
- **It does NOT give a renormalon-free absolute `m_{0⁺⁺}`.** The strut is the **dilaton / IR scale
  `z₀,c,R`**, an input; the ambiguity is **relocated** into the IR boundary condition, not removed.
- **Not rigorous:** no proven YM↔geometry duality; the model-to-YM map is uncontrolled.
- **Net for the roof:** holography is a **model of the missing strut** — it tells you the strut is an IR
  scale that must be *generated*, which is exactly `H1` / dimensional transmutation. It **confirms the
  shape** of the gap (a single IR scale) and **confirms the consistency** of `m/√σ` with the strong-coupling
  bridge, but it does **not** supply the uniform `χ_lattice ≲ a²` floor. **H1 remains open.**

*Refs:* Polchinski–Strassler hep-th/0101126; Karch–Katz–Son–Stephanov hep-th/0602224; Witten hep-th/9803131;
Brower–Mathur–Tan hep-th/0003115 (*Nucl. Phys.* B587 (2000) 249); Sakai–Sugimoto hep-th/0412141;
lattice Morningstar–Peardon hep-lat/9901004, Michael–Teper PLB 311 (1993) 408. Witness:
`w7/roof_holography.py`.
