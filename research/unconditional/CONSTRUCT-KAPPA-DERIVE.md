# CONSTRUCT-KAPPA-DERIVE — the rigid-string rigidity κ is **not** derivable from σ or the universal Lüscher term; it is `m/√σ` re-expressed (= the gap), and its number is the microscopic core datum = `H1`

> **VERDICT (no QED).** Three routes were built to derive `κ` (the dimensionless rigidity in
> `m_{0⁺⁺}/√σ = 2π√(2κ)`, fitting `κ∈[0.138,0.193]`). **None derives the number.** **(1) [NG] the
> Lüscher/width route:** the flux-tube width's `ln R` slope is **universal** (`D_⊥/(2πσ)`, `κ`-independent
> to `5×10⁻⁶`); `κ` only shifts the *additive constant* (the effective core radius) — a UV/microscopic
> datum the width does not determine. **(2) [T] circularity:** inverting the relation gives
> `κ = m²/(8π²σ)`, so *deriving κ ≡ deriving `m/√σ` = the gap* — κ is not a softer handle. **(3) [V-num]
> dimensional analysis:** `κ = c·σ ξ²` (dimensionless); at the glueball scale `ξ~1/√σ`, `κ~c ∈ O(1)`,
> consistent with the fitted band but `c` is the microscopic coefficient. **Exact missing input:** the
> microscopic flux-tube rigidity (core radius / profile) = `H1` in string coordinates.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric witness.
Units `σ=1`; `κ` = coefficient of `∫K²ds` in the rigid surface-sum of `CONSTRUCT-MASTERFIELD`.

---

## §1 · The Lüscher/width route cannot fix κ `[NG]`

Gaussian rigid worldsheet, transverse modes `x_⊥`: propagator `G(k)=1/(σk²+κk⁴)`. The flux-tube
(coincident transverse) width is

`w²(R) = D_⊥ ∫ d²k/(2π)² · 1/(σk²+κk⁴) = D_⊥/(2π) ∫_{1/R}^{1/a} dk /(σk+κk³)`.

Closed form: `∫dk/(k(1+βk²)) = ½ ln(k²/(1+βk²))`, `β=κ/σ`. As `R→∞` the lower limit gives `+ln R`, so

`w²(R) = [D_⊥/(2πσ)]·ln R + const(κ)`,   **slope `D_⊥/(2πσ)` — independent of κ.**

Witness (`w7/construct_kappa_derive.py`): slope `0.318310` at `κ=0`, `0.318305` at `κ=0.193` (spread
`4.6×10⁻⁶`), matching the predicted `D_⊥/(2πσ)=1/π`. `κ` moves only the **intercept** (shift `−0.764`
from `κ=0` to `0.193`) — i.e. the effective core radius `r₀`, a short-distance/microscopic datum.

⇒ The universal Lüscher coefficient (Aharony–Komargodski, universal through `1/L⁵`) is a statement about
the *conformal* transverse modes and is **κ-blind**; the rigidity is precisely the *non-universal*
`1/L⁷`-and-beyond data (`PROOF-EFFECTIVE-STRING`). **The width cannot derive κ.**

---

## §2 · κ is `m/√σ` re-expressed — not an independent handle `[T]`

The master-field relation is `m_{0⁺⁺} = 2π√(2σκ)`, hence

`m² = 8π²σκ  ⇒  κ = m²/(8π²σ) = (m/√σ)²/(8π²)`.

Witness: `m/√σ=3.90 ⇒ κ=0.1926`; `m/√σ=3.30 ⇒ κ=0.1379`. **κ and `m/√σ` are the same unknown.** So
"derive κ from first principles" is **logically identical** to "derive the `0⁺⁺` mass given σ" — the gap
itself, not a softer datum. This is the cleanest sense in which the rigid-string model, though a genuine
advance (`CONSTRUCT-MASTERFIELD`), does not create a new lever: its one input *is* the answer.

---

## §3 · Dimensional analysis pins the scaling, not the number `[V-num]`

`κ` is dimensionless (`∫K²ds`: `K~1/L`, `K²~1/L²`, `ds~L` ⇒ dimensionless). The only scale is the core
radius `ξ` and the tension `σ`, so `κ = c·σξ²` with `c` an O(1) coefficient. At the glueball scale
`ξ ~ 1/√σ` (the tube is ~as wide as it is long, `PROOF-EFFECTIVE-STRING`), so

`κ = c·σ·(1/σ) = c`.

The fitted band `κ∈[0.138,0.193]` ⇒ `c∈[0.138,0.193]`: an **O(1) number**, exactly the expected size —
but `c` is set by the **microscopic flux-tube profile** (the Nielsen–Olesen / dual-Abrikosov core
structure at `β≈2.5`), which is itself the nonperturbative confinement-scale datum.

---

## §4 · PK-EFT RG: κ is a boundary condition, not a flow output `[C]/[O]`

In the Polyakov–Kleinert effective string the rigidity is a 4-derivative worldsheet coupling (marginal-
ish), renormalized along the worldsheet RG. Its **IR value is a boundary condition supplied by the UV
completion** (the microscopic gauge theory at the confinement scale `Λ`); the EFT flow does not select it.
So the EFT can propagate `κ` from one scale to another but **cannot manufacture its value** — the same
structural fact as `ROOF-CONDENSATE-FLOOR`/`TRANSMUTE-STRING-SPECTRUM` for the core thickness: it is a
`Λ`-datum.

---

## §5 · Bottom line

- **Lüscher/width route: `[NG]`** — slope universal, κ pushes only the additive constant.
- **Circularity: `[T]`** — `κ = m²/(8π²σ)`; deriving κ ≡ deriving the gap.
- **Dimensional analysis: `[V-num]`** — `κ~c∈O(1)`, right magnitude, but `c` is microscopic.
- **Exact missing input `[O]`:** the microscopic flux-tube rigidity / core radius (the `ξ` at which the
  tube forms) — equivalently `Λ` — = `H1` in string coordinates, matching `PROOF-EFFECTIVE-STRING`
  ("missing = the core thickness ~1/√σ") and `TRANSMUTE-STRING-SPECTRUM` ("a non-EFT string-core
  coercivity"). No new lever; the rigid-string model's single parameter is the gap itself.

*Refs.* Lüscher–Symanzik–Weisz, NPB **173** (1980) 365; Lüscher–Weisz, JHEP **07** (2002) 049
(hep-lat/0207003); Aharony–Komargodski, JHEP **05** (2013) 118 (arXiv:1302.6257, long-string universality
through `1/L⁵`); Polyakov, *Gauge Fields and Strings* (rigid string, rigidity RG); Kleinert (extrinsic-
curvature rigid string). Internal: `CONSTRUCT-MASTERFIELD.md` (the relation), `PROOF-EFFECTIVE-STRING.md`,
`TRANSMUTE-STRING-SPECTRUM.md`. Witness: `w7/construct_kappa_derive.py`. **Clay 0%; no QED.**
