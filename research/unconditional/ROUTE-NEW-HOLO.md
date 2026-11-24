# ROUTE-NEW-HOLO — can a bulk-energy / NEC / c-theorem argument bound χ ≤ 1/λ (or m₀⁺⁺) from the dual geometry, without assuming the IR wall? **No — it reduces to the IR wall scale = H1.**

> **VERDICT `[NG]`.** The holographic route gives an **exact** relation `m₀⁺⁺ z₀ = j_{2,1} = 5.1356` and places the *whole* spectrum and χ on the single scale `z₀` (the IR wall) — so bounding the gap or χ **⟺ bounding `z₀`**. But bulk positivity (NEC), the holographic c-theorem, and bulk-energy are **scale-invariant** and cannot fix `z₀`: pure AdS is NEC-**saturated** (`T_{MN}k^Mk^N = −Λ·0 = 0` for every null `k`) for *any* `z₀`. A dynamical wall sets `z₀` from a potential scale that is itself the free confinement scale. **The holographic bound = assuming the IR wall at a finite scale `z₀` = the confinement scale = H1 in the bulk-geometry coordinate.** This is *why* holography predicts the scale-free **ratio** `m/√σ = O(1)` (which is how it "matched" 3.9) but never the absolute scale. Same wall, graviton/dilaton coordinate. Clay 0%; no QED.

---

## §0 · Labels

`[T]` 真定理 · `[C]` 条件定理 · `[O]` 精确缺失步 · `[NG]` no-go · `[V-num]` new numeric.

## §1 · The exact in-model relation (real, `[T]`)

Hard-wall holography (Polchinski–Strassler; Boschi–Filho–Braga) gives the 0⁺⁺ glueball spectrum **exactly**: `m_n z₀ = j_{2,n}` (Bessel zeros of `J₂`), so

  `m₀⁺⁺ = 5.1356 / z₀`,  `m₀⁺⁺*/m₀⁺⁺ = j_{2,2}/j_{2,1} = 1.639` (lattice ~1.6–1.64).

Every level, and hence `χ = Σ f_n²/m_n²`, is a function of `z₀` alone. Consequence (`[T]`):
- `m₀⁺⁺ ≥ δ  ⟺  z₀ ≤ j_{2,1}/δ`  (bound the gap from below ⟺ bound the IR wall from **above**);
- `χ ≤ 1/λ  ⟺  z₀ ≤ (bound)/√λ`.

Both are statements about one number, `z₀`.

## §2 · Why bulk positivity cannot fix `z₀` (`[T]`)

- **NEC is saturated in pure AdS** (`[T]`): `T_{MN}k^Mk^N = −Λ g_{MN}k^Mk^N = 0` for every null `k`, for **any** `Λ` and any `z₀` (witness: `6.0 × 0 = 0`). NEC constrains the geometry's **shape**, never the **scale**.
- **Holographic c-theorem** (Freedman–Gubser–Klebanov–Rocha):
  fixes the number of IR degrees of freedom `a(IR)`, not the scale at which the flow ends.
- **Dynamical wall / dilaton flow / cigar:** sets `z₀` from a superpotential / potential scale `Λ`, which is itself the **free confinement scale** (`[C]`, model input).

So no bulk-energy / NEC / c-theorem inequality pins `z₀`.

## §3 · The consequence (`[NG]`)

Holography's success at `m/√σ ≈ 3.9` is **kinematic, not a gap**: with `σ ~ c/z₀²`, `m/√σ = j_{2,1}/c` is a **pure number** (witness: `5.1356/c`), independent of `z₀` — it fixes the geometry ratio `c`, not the absolute scale. Bounding `m₀⁺⁺` or `χ` needs `z₀` bounded, i.e. the IR wall assumed at a finite scale = **the confinement scale = H1**.

## §4 · Honest bottom line

- The holographic/gravity route **does not** bound `χ ≤ 1/λ` or the 0⁺⁺ mass independently of the IR wall.
- It **reduces to** `z₀` (the wall/cigar scale) = the confinement scale = **H1 in the bulk-geometry coordinate**.
- Genuine yield: the exact `m₀⁺⁺ z₀ = j_{2,1}` relation and a clean explanation of why holographic *ratios* succeed while the absolute scale is never predicted — a structural identification, not a lever.
- **Clay 0%.** No QED.

*Refs:* `KAPPA-HOLO-COMPUTE.md`, `ROOF-HOLOGRAPHY.md`, `CONSTRUCT-HOLO-FLUXTUBE.md`, `ONLY-EXIT-SYNTHESIS.md`; witness `w7/route_new_holo.py`.
