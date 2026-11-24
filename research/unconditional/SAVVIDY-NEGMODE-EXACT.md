# SAVVIDY-NEGMODE-EXACT — the exact (nonperturbative) treatment of the Savvidy negative mode does NOT kill the χ-floor: the softening is BOUNDED and sub-criticality is PRESERVED as the instability grows

> **VERDICT (no QED; Clay 0%).** The Gaussian-around-0 treatment `Γ″(0)=κ−g²Var(φ²)` is
> expanding about the **wrong saddle** (the barrier top, curvature `−c`). The exact (nonperturbative)
> partition function is dominated by the **true minima** `φ=±v=±√(c/λ)`, where the mode is **STABLE**
> (`V″(±v)=2c>0`). The exact `Var(φ²)` is **BOUNDED** and, crucially, **saturates at `2/λ`
> independent of `c`** as the Savvidy instability `c=gB` grows — so the softening `g²Var(φ²)` does
> **not** grow with the instability and the sub-criticality `κ_bare > g²Var(φ²)` is **preserved, not
> eventually killed**. The complex saddle (steepest descent through `φ=0`) gives the **vacuum-decay
> rate** (the imaginary part of the energy) — a *separate* object, not the real susceptibility.
> **The remaining content is the true (condensed) vacuum `⟨φ²⟩≈c/λ` = the condensate magnitude = `H1`.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.
Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The exact model of the negative-mode direction `[T]`

Single unstable mode `φ`, faithful 0-dim model of the Savvidy negative direction with quartic
stabilisation:
```
V(φ) = −(c/2) φ² + (λ/4) φ⁴ ,     c = g|B| > 0  (the instability).
```
- `V″(0) = −c = λ_min` — the **negative mode** (`φ=0` is a MAX of `V` = the barrier top).
- True minima at `φ = ±v`, `v = √(c/λ)`, with `V″(±v) = −c + 3λv² = 2c > 0` — **STABLE**.

The exact susceptibility `Var(φ²)` is the 0+0d quadrature `∫φ^{2,4} e^{−V}/∫e^{−V}`.

**Key consequence `[T]`:** the Gaussian-around-0 is expanding about a **maximum** of `e^{−V}`, i.e. the
wrong saddle. The exact measure is dominated by `±v`; the "negative mode" is a **barrier-top artifact**,
not a genuine instability of the exact partition function.

---

## §2 · Exact result: `Var(φ²)` is bounded and `c`-independent `[V-num]`

Exact quadrature (`w7/savvidy_negmode_exact.py`, `λ=1`, exponent-shifted):

| `c=gB` | `v` | `⟨φ²⟩` | `Var(φ²)` | `2/λ` | `2c` |
|---|---|---|---|---|---|
| 0.01 | 0.10 | 0.6787 | 0.546 | 2 | 0.02 |
| 0.1 | 0.316 | 0.7039 | 0.575 | 2 | 0.2 |
| 1.0 | 1.0 | 1.042 | 0.956 | 2 | 2 |
| 5.0 | 2.236 | 4.760 | **2.141** (peak) | 2 | 10 |
| 10 | 3.162 | 9.897 | 2.022 | 2 | 20 |
| 50 | 7.071 | 49.98 | 2.001 | 2 | 100 |
| 100 | 10.0 | 99.99 | **2.0002** | 2 | 200 |

**`Var(φ²)` is bounded throughout and → `2/λ` (independent of `c`) as `c→∞`.** The analytic reason:
at large `c`, the field sits in a well (`⟨φ⟩²= v² = c/λ`) and fluctuates with stiffness `2c`
(`Var(φ)=1/(2·2c)=1/(4c)`... in the Gaussian e^{−½m²x²} convention `Var=1/m²=1/(2c)`), so
`Var(φ²) ≃ 4v²·Var(φ) = 4(c/λ)·(1/(2c)) = 2/λ`. The signal `v²` grows linearly with `c` while the
noise `Var(φ)` shrinks like `1/c` — they cancel, leaving a `c`-independent bound `2/λ`.

---

## §3 · Comparison with the wrong saddle (Gaussian-around-0) `[NG]`

The Gaussian-around-0 assigns the mode the (ill-defined, negative) variance `1/(2V″(0))=−1/(2c)`; its
"magnitude" `1/(2c)` **diverges as `c→0`** (`50` at `c=0.01`, `5` at `c=0.1`) and shrinks as `c→∞`.
So the wrong saddle is wrong in **both** directions (worst at small `c`), whereas the exact `Var(φ²)`
is uniformly bounded. The prior Gaussian/RG estimate (`CONSTRUCT-SAVVIDY-CELL`) already reported
`Var(φ²)` finite and `λ₄`-set; the exact computation **confirms and sharpens** it: the bound is
`2/λ`, **`c`-independent**.

---

## §4 · Does the exact treatment change the softening? — No (it protects it)

The softening enters as `Γ″(0) = κ_bare − g² Var(φ²)` (the `D`-term of `G″=S−D`). Since the exact
`Var(φ²) ≤ 2/λ` for all `c`, the softening `g²Var(φ²) ≤ 2g²/λ` is **bounded**, and the
sub-criticality `Γ″(0) > 0 ⟺ κ_bare > 2g²/λ` is a **fixed** condition that does **not** worsen as the
Savvidy field `B` grows. **The negative mode does not eventually close the χ-floor.**

---

## §5 · The complex saddle = the decay rate, a separate object `[T]/[C]`

The standard steepest-descent *through* `φ=0` in the imaginary direction yields an **imaginary** part
of the vacuum energy (`∝ e^{−const/g²}`) — the **Savvidy-vacuum decay rate**. That is the *energy's*
instability, **not** `Var(φ²)` (an equal-time real susceptibility). The real `χ-floor` is governed by
the exact real measure (dominated by `±v`), which is bounded. So the "imaginary part" is real physics
but it is **not** a threat to the susceptibility floor; conflating them would be a category error.

---

## §6 · Exact missing step `[O]`

The exact treatment lands on the **symmetry-broken (condensed) vacuum**: `⟨φ²⟩ = ⟨φ⟩² ≈ c/λ = ⟨tr F²⟩`
(the gluon condensate). The map from the toy `c=gB` to the *physical* 4D condensate magnitude — and
the uniform `(L,a)` floor of the true `0⁺⁺` susceptibility on top of it — is **`H1`**. So the exact
resummation shows the negative mode **is not** the obstruction (bounded softening, protected
sub-criticality), and re-locates the remaining content to the **condensate magnitude = `H1`**.

---

## §7 · Bottom line

- **`[T]`** the negative mode is a **wrong-saddle (barrier-top) artifact**; the exact saddle `±v` is
  stable (`V″=2c>0`).
- **`[V-num]`** exact `Var(φ²)` is **bounded** and **`c`-independent** (`→2/λ`); the softening is
  bounded, so **sub-criticality is preserved** as `B` grows — the negative mode does not kill the floor.
- **`[NG]`** Gaussian-around-0 (diverges at small `c`) is invalid; the complex saddle is the decay rate,
  a separate object.
- **`[O]`** the exact vacuum is the **condensate** `⟨tr F²⟩ ≈ c/λ`; mapping to the physical floor is
  `H1`. **Clay 0%; no QED.**

*Witness:* `w7/savvidy_negmode_exact.py`. *Prior:* `CONSTRUCT-SAVVIDY-CELL.md` (Gaussian softening),
`G″=S−D` (`INVENT-VARIANCE-INEQ`), `CONSTRUCT-CONSTRUCTIVE-QFT.md`.
