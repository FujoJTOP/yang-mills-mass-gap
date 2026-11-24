# KAPPA-HOLO-COMPUTE — the dual string is EXACTLY Nambu–Goto at leading α′ (κ = 0); the atom is GEOMETRIC, not rigid-string; the α′σ ≈ 0.16 coincidence is not a derivation

> **VERDICT (no QED).** The dual string in a confining holographic bulk (Witten D4 cigar / hard wall /
> soft wall) is a **pure Nambu–Goto string at leading order in α′** — it carries **no worldsheet rigidity**,
> so `κ_holographic = 0` in the MASTERFIELD sense. The 0⁺⁺ gap/atom in these models is supplied by the
> **IR geometry** (hard-wall Bessel zeros; Witten-cigar KK tower), *not* by a rigid-string `κ`. Thus
> holography does **not** compute the flux-tube rigidity `κ ∈ [0.138,0.193]`; `κ` (flux-tube) and the
> holographic string are **different objects**. The natural α′-correction rigidity is `O(α′σ)`; numerically
> `α′σ = 1/(2π) ≈ 0.159` happens to sit inside the band — a **normalization coincidence, not a derivation**
> (the `O(1)` coefficient is model-dependent and the object differs).

Labels: **[T]** proved/standard · **[C]** conditional/model · **[O]** exact missing step · **[NG]** no-go for the named class · **[V-num]** numeric.

---

## §1 · The dual string is Nambu–Goto at leading α′ — κ = 0

The holographic Wilson-loop string has worldsheet action
```
S = (1/2πα′) ∫ d²σ √(−det G_ab),   G_ab = G_MN(X) ∂_a X^M ∂_b X^N ,
```
where `G_MN` is the bulk (string-frame) metric `[T]` (Maldacena; Rey–Yee). At **leading order in α′** this
is *exactly* the Nambu–Goto action — there is **no `(∂²X)²` rigidity term**. The rigidity would be a
**four-derivative** worldsheet operator, which enters only at higher order in the `α′` (worldsheet-derivative)
expansion: schematically
```
(κ/2)(∂²X)²  ∼  O(α′) × (two-derivative NG term) ,
```
so `κ_leading = 0` and `κ` is suppressed by `α′/ℓ²` at worldsheet length scale `ℓ` (witness: `1/ξ²` scaling).
**The dual string is a thin NG string.** `[T]`

## §2 · The gap/atom is GEOMETRIC, not rigid-string

The 0⁺⁺ mass in holographic models is a **bulk KK mode's** energy, fixed by the IR geometry:

- **Hard wall** (Polchinski–Strassler / Boschi–Filho–Braga): transverse bulk scalar with `ν = 2` gives
  `m_n z₀ = j_{2,n}` — the **Bessel zeros** `j_{2,1}=5.1356, j_{2,2}=8.4172, …`, a discrete tower
  `1 : 1.64 : 2.26 : 2.88` (matches lattice `1 : ~1.6 : ~2.3 : ~2.9`). The discreteness — hence the **atom**
  `f₁²>0` — comes from the **wall** (a boundary condition on the bulk mode), not from any worldsheet rigidity.
- **Witten D4 cigar**: both `m` and `√σ` are set by the same `M_KK`; the ratio `m/√σ = O(1) ≈ 3.5–4`
  (Brower–Mathur–Tan), with the correct ordering `m(0⁺⁺)<m(2⁺⁺)<m(1⁺⁺)` — again a **geometric** statement.

So `[T]`: in the holographic construction the atom exists and is gapped, **but its gap is the IR geometry's
KK scale**, not `2π√(2σκ)` from a rigid worldsheet. **The two mechanisms are distinct.** `[T]/[C]`

## §3 · Why κ (MASTERFIELD) ≠ holographic string rigidity

The MASTERFIELD `κ` is the **rigidity of the Yang–Mills flux tube** — a 4-derivative term
`W(C)=Σexp(−σ·Area − κ∮K²ds)` stabilized around an interior minimum `ρ*=√(κ/2σ)`. The holographic string is a
**fundamental** string in a confining *bulk*; its 0⁺⁺ is a **bulk KK quantum**, stabilized by the warp
factor / wall, not by a flux-tube rigidity. Conflating them is a type error (cf. `CONSTRUCT-WIDTH-KAPPA-LATTICE`:
PK `α` ≠ MASTERFIELD `κ`). **No controlled map gives `κ_holographic = κ_MASTERFIELD`.** `[NG]`

## §4 · The α′σ ≈ 0.16 coincidence (flagged, not derived)

If one nonetheless asks "what is the natural α′-rigidity in units of the tension?", the only dimensionless
combination is `α′σ`. With the standard normalization `σ = 1/(2πα′)` this is
```
α′σ = 1/(2π) ≈ 0.15915  ∈ [0.138, 0.193]   [V-num]
```
— sitting **inside** the MASTERFIELD band. This is a **coincidence of normalization**, not a derivation:
(i) the `O(1)` coefficient of the induced rigidity is model-dependent (it depends on `λ = g²N`, the warp
factor, the dilaton); (ii) the object (fundamental-string α′ correction) is not the flux-tube rigidity that
`m/√σ = 2π√(2κ)` uses; (iii) the holographic atom is geometric and does not need this rigidity at all.
Recording it `[V-num]` only so it is not later mistaken for a derivation.

## §5 · Verdict

- **`κ_holographic = 0`** at leading order — the dual string is exactly NG `[T]`.
- **The atom is geometric** (hard-wall Bessel / cigar KK), not rigid-string `[T]/[C]`.
- **Holography does not compute the flux-tube `κ`** — different object `[NG]`.
- **`α′σ ≈ 0.159 ∈ [0.138,0.193]` is a coincidence**, not a prediction `[V-num]/[NG]`.
- **Exact missing step `[O]`:** the flux-tube rigidity `κ` (equivalently the core thickness `~1/√σ`), which
  holography replaces by *assuming* an IR wall / cigar — the same `H1` in the bulk-geometry coordinate
  (`CONSTRUCT-HOLO-FLUXTUBE`).

**Clay 0%; no QED.** This computes the dual worldsheet rigidity, honestly: it is zero, and the holographic
atom rides on the geometry, not on `κ`.

*Refs:* `CONSTRUCT-MASTERFIELD.md`, `CONSTRUCT-HOLO-FLUXTUBE.md`, `CONSTRUCT-KAPPA-DERIVE.md`,
`CONSTRUCT-WIDTH-KAPPA-LATTICE.md`, `CONSTRUCT-WORLDSHEET-CFT.md`; Maldacena hep-th/9803002,
Rey–Yee hep-th/9803001, Witten hep-th/9803131, Brower–Mathur–Tan hep-th/0003115.
