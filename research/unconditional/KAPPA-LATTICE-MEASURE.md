# KAPPA-LATTICE-MEASURE — the sub-leading `1/R²` width coefficient IS a clean, `r₀`-free rigidity observable (`c₂=(D−2)κ_PK/4π`), but it is <0.1% of `w²` at large `R` (below current precision); and `κ_PK` ≠ MASTERFIELD `κ` (no dictionary)

> **VERDICT (no QED).** The rigid-string rigidity does have a *clean* lattice observable — the
> **sub-leading `1/R²` coefficient** of the flux-tube width, `c₂=(D−2)/(4π)·κ_PK` (σ=1 units), which is
> **positive, proportional to the rigidity, and carries NO `r₀` (UV-core) degeneracy** (unlike the constant
> term `w₀²` of `WIDTH-KAPPA`, which is `r₀`-degenerate). This is strictly better than the previously
> considered channels. **But** (i) `c₂/R²` is only ~3.8% of `w²` at `R=2`, falling to **<0.1% at `R≥8`** —
> 10–50× below current lattice width precision, so **not extractable from existing data**; (ii) it yields the
> **Polyakov–Kleinert rigidity `κ_PK`**, *not* the MASTERFIELD `κ` (no established dictionary, `WIDTH-KAPPA`);
> (iii) NG's own nonlinear corrections contaminate at the same order (fixed, calculable, but a systematic).
> The MASTERFIELD `κ∈[0.138,0.193]` therefore still comes **only** from `m_{0⁺⁺}/√σ` = the gap (circular).
> **Clay 0%.**

---

## §0 · Labels
**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric witness.
Units `σ=1`, UV cutoff `a=1` (`KAPPA-DERIVE` convention); `κ_PK` = coefficient in `G(k)=1/(σk²+κ_PK k⁴)`.

---

## §1 · The observable and its derivation `[T]`

Transverse propagator `G(k)=1/(σk²+κ_PK k⁴)`. The (coincident-transverse) flux-tube width is

`w²(R) = (D−2)/(2π) ∫_{1/R}^{1/a} dk / (k(σ + κ_PK k²))`.

Closed form (`σ=1, a=1`, verified against `scipy.quad` to `10⁻¹⁷`):

`w²(R) = (D−2)/(2π)·[ ln R − ½ ln(1+κ_PK) + ½ ln(1+κ_PK/R²) ]`
      `= (D−2)/(2π)·ln R  − (D−2)/(4π)·ln(1+κ_PK)  + (D−2)/(4π)·κ_PK/R² + O(R⁻⁴)`.

- **ln R slope** `(D−2)/(2π)`: universal, **κ-blind** (`KAPPA-DERIVE §1`).
- **constant** `−(D−2)/(4π)ln(1+κ_PK)`: a rigidity signal, but shifted identically by `w₀²` (UV core) ⇒
  **`r₀`-degenerate** (`WIDTH-KAPPA §2–3`).
- **sub-leading `1/R²` coefficient**  **`c₂ = (D−2)/(4π)·κ_PK`**  ⇒ **positive, `∝ κ_PK`, and `r₀`-FREE.**

`c₂` is the *only* width channel where the rigidity appears without the UV-core degeneracy. NG (`κ_PK=0`)
gives `c₂=0` at the quadratic level (pure log), so `c₂>0` is a **direct PK-rigidity signature**.

**Protocol.** Measure `w²(R)` on `R×T` Wilson/Polyakov loops; fit `w² = A ln R + B + C/R²` with `A` fixed
to the universal `(D−2)/(2π)`; the free `C=c₂` gives `κ_PK = 4π c₂/(D−2)`. Systematics: subtract NG's
nonlinear (determinant) corrections to `w²` at the same order (fixed coefficients, in principle known); use
large `R` (e.g. `R/√σ ≥ 4`, `R ≥ 8/√σ`) where the effective string is reliable; sum transverse directions.

---

## §2 · Resolvability `[V-num]`

From `w7/kappa_lattice_measure.py` (`κ_PK=0.193`, `c₂=0.0307`):

| `R` | `c₂/R²` | `w²(R)` | `(c₂/R²)/w²` |
|---|---|---|---|
| 2 | 0.00768 | 0.200 | **3.84 %** |
| 4 | 0.00192 | 0.415 | 0.46 % |
| 6 | 0.00085 | 0.543 | 0.16 % |
| 8 | 0.00048 | 0.634 | **0.076 %** |
| 12 | 0.00021 | 0.763 | 0.028 % |

To resolve `c₂` at `R≈8–10` needs **absolute** precision `δ(w²) ≲ 3×10⁻⁴` (rel `≲0.05–0.08%`). The width is
reliable only at large `R`, where the effect is *smallest*; at small `R` (where it reaches ~4%) the effective
string has large corrections. **Current lattice width precision is a few %** ⇒ roughly **10–50× short**.
A noiseless 3-parameter fit recovers `c₂` (`c₂_fit=0.0298` vs `0.0307`), confirming the protocol is correct,
but with realistic noise it is unresolved.

Similarly, the **energy/potential** rigidity enters only at `O(1/L⁷)` beyond the universal
(Aharony–Komargodski) `1/L⁵` — an even smaller channel.

---

## §3 · Published flux-tube data: NG-consistent, `c₂` unresolved

The lattice flux-tube width literature (Lüscher 1981; Caselle–Gliozzi et al.; Bali–Schlichter–Schilling;
Caselle–Necco) and effective-string tests (Aharony–Komargodski; Athenodorou–Bringoltz–Teper) report:
the `ln R` slope matches `(D−2)/(2πσ)` (universal), the intrinsic width `w₀ ≈ 0.2 fm`, and **overall
Nambu–Goto consistency**; the sub-leading rigidity is consistent with zero at current precision. Hence the
data do **not** resolve `c₂` (let alone to a `κ_PK` number at the few-% precision the sign/dictionary needs).
⚠ The specific `w₀²` and precision figures here are representative of that literature — verify against the
source tables (as `WIDTH-KAPPA` also flagged).

---

## §4 · The two remaining obstructions even if `c₂` were resolved `[O]`

1. **Normalization mismatch:** `c₂` measures the **dimensionful PK rigidity `α`/`κ_PK`** (worldsheet EFT
   coefficient). The MASTERFIELD `κ` (dimensionless, `m_{0⁺⁺}/√σ=2π√(2κ)`) is a **different parameter**;
   `WIDTH-KAPPA §4–5` found **no clean dictionary**. So a measured `c₂` gives a *new number*, not the
   MASTERFIELD `κ` of the flux-tube-glueball model.
2. **Tautology of the fitted band:** the MASTERFIELD `κ∈[0.138,0.193]` is *defined* by `κ=(m_{0⁺⁺}/√σ)²/(8π²)`
   (`KAPPA-DERIVE §2`). The lattice "gives" this band **exactly to the extent it measures the glueball mass
   `m_{0⁺⁺}/√σ`** = the gap itself. So `κ` is not an independent handle; measuring it *is* measuring the gap.

---

## §5 · Verdict

- The **sub-leading `1/R²` width coefficient is the correct, clean, `r₀`-free** lattice rigidity observable
  (`c₂=(D−2)κ_PK/4π > 0`), a genuine improvement over the constant-term channel — **worth measuring with
  future higher-statistics width data** (needs ~10–50× better precision at large `R`).
- **It is not extractable from existing data**, and even then gives `κ_PK`, not MASTERFIELD `κ`.
- The policy-level conclusion is unchanged: **`κ` (MASTERFIELD) is not an independent lattice handle** —
  it equals `m_{0⁺⁺}/√σ` in re-expression, so "measure κ" = "measure the gap." `κ=H1` stands.

*Refs:* `CONSTRUCT-KAPPA-DERIVE.md`, `CONSTRUCT-WIDTH-KAPPA-LATTICE.md`, `CONSTRUCT-MASTERFIELD.md`;
lattice width: Lüscher 1981, Münster–Weisz, Caselle–Gliozzi, Bali–Schlichter–Schilling, Caselle–Necco;
effective string: Aharony–Komargodski 2010, Athenodorou–Bringoltz–Teper 2011. Witness:
`w7/kappa_lattice_measure.py`.
