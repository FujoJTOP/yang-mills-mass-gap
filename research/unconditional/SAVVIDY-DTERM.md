# SAVVIDY-DTERM — the Savvidy softening in the exact variance-identity form: `D = D₊ + E|⟨φ₋|c⟩|²/(g|B|)`, `D₊ ≤ 0` (BL) ⇒ `D < S ⟺ E|⟨φ₋|c⟩|²/(g|B|) < S_bare`; uniform `D<S` at the crossover is `H1`

> **VERDICT (no QED).** The Savvidy softening is the `D`-term of the established identity
> `G″ = S − D`. Writing `K⁻¹ = K₊⁻¹ − P₋/|λ₋|` (the Savvidy Hessian is `H₊` minus a **rank-one**
> negative mode), the defect splits exactly into a Brascamp–Lieb part `D₊ ≤ 0` (positive fiber) plus a
> positive Savvidy term `E[|⟨φ₋|c⟩|²]/(g|B|)`. Hence the two-sided statement:
> **`D < S ⟺ G″ > 0`** is implied by **`E[|⟨φ₋|c⟩|²]/(g|B|) < S_bare`**. This is exact. Deep IR the overlap
> vanishes by scale separation ⇒ `D/S → 0` (uniform). At the crossover `k ~ √(gB)` the ratio is `O(1)` and
> **`D<S` is NOT uniformly established — that open inequality is `H1` (softened).** **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件/模型 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** 本轮数值.

**Definitions (exact, SLOT-SD-FIXPOINT / INVENT-VARIANCE-INEQ).** Flat chart, `μ ∝ e^{−V}dX`,
coarse mode `v`, `a=⟨v,X⟩`, condition on `a=0`; `h_n = vᵀ Hess V v`, `f=∂_nV`, `c=Π₀ Hess V v`,
`K=Π₀ Hess V Π₀`. Then
```
G″(0) = S − D,   S = E_μ[h_n − ⟨K⁻¹c,c⟩|0],   D = Var_μ(f|0) − E_μ[⟨K⁻¹c,c⟩|0].
```

## §1 · Exact decomposition in the Savvidy cell `[T]`

Savvidy Hessian: `Hess = H₊ − |λ₋| P₋`, `λ₋ = −g|B|`, `P₋` **rank one** (lowest-Landau-level mode `φ₋`).
The fiber inverse is `K⁻¹ = K₊⁻¹ − P₋/|λ₋|`, so
```
⟨K⁻¹ c, c⟩ = ⟨K₊⁻¹ c, c⟩ − |⟨φ₋|c⟩|² / (g|B|).
```
Substituting into `D`:
```
D = [Var(f|0) − E⟨K₊⁻¹c,c⟩]  +  E[ |⟨φ₋|c⟩|² ] / (g|B|)
  =:        D₊                 +      D_Sav
```
- **`D₊ ≤ 0` `[T]`**: on the positive (log-concave) fiber `K₊ ⪰ 0`, Brascamp–Lieb gives
  `Var(f|0) ≤ E⟨K₊⁻¹c,c⟩`. Witness (1): `D₊ = −0.72, −0.03, −0.02, …` (always `≤ 0`).
- **`D_Sav = E[|⟨φ₋|c⟩|²]/(g|B|) ≥ 0`** — the negative-mode enhancement.

Hence the **exact two-sided bound**
```
D  ≤  D_Sav  =  E[ |⟨φ₋|c⟩|² ] / (g|B|),
```
and therefore the **conditional exact theorem**
```
(★)   E[ |⟨φ₋|c⟩|² ] / (g|B|)  <  S_bare       ⟹   D < S   ⟺   G″ > 0.
```
`(★)` is `D<S` written purely in the Savvidy cell — the variance-identity form of the softened sub-criticality.

## §2 · Toy RG closure and `D/S` `[C]/[V-num]`

Project to the negative-mode amplitude with `V(φ) = −½cφ² + ¼λ₄φ⁴` (`c=|λ₋|>0`); integrating out `φ`
gives `Γ″(0) = κ − g²Var(φ²)`, `Var(a)=1/Γ″`. Identify `S = κ` (bare) and `D = g²Var(φ²)` (the softening),
so `G″ = κ − g²Var(φ²)` and **`D/S = g²Var(φ²)/κ`**. Witness (2)/(3):

| `c` | `λ₄` | `Var(φ²)` | `D/S` @ `κ=1.5` | `@2.0` | `@3.0` |
|---|---|---|---|---|---|
| 1 | 1 | 0.956 | 0.64 | 0.48 | 0.32 |
| 1 | 0.5 | 2.381 | **1.59** | 1.19 | 0.79 |
| 1 | 2 | 0.406 | 0.27 | 0.20 | 0.14 |
| 2 | 1 | 1.557 | **1.04** | 0.78 | 0.52 |
| 4 | 1 | 2.209 | **1.47** | 1.10 | 0.74 |
| 1 | 0.2 | 8.513 | **5.68** | 4.26 | 2.84 |

`D/S = soft/κ` is **`O(1)`** and **can exceed 1**: e.g. `λ₄=0.5, κ=1.5 ⇒ D/S=1.59 ⇒ D>S ⇒ G″<0`.
The sign flip is exactly at `κ = soft` (`D/S=1, G″=0`); below it, `Var(a)=∞` (floor destroyed). So `D<S`
holds **iff `κ > g²Var(φ²)`** — an `O(1)` competition, not an automatic inequality.

## §3 · Uniformity: deep IR yes, crossover no `[T]/[O]`

- **Deep IR (`k ≪ √(gB)`) `[T]`**: the negative mode `φ₋` is a lowest-Landau-level state with Fourier
  support `|p| ~ √(gB)`; the `k→0` mode has support `|p|<k`, so the overlap `⟨φ₋|c⟩` is suppressed by
  `(k/√(gB))² → 0`. Witness (4): `overlap² = 0.010, 0.083, 0.200, 0.500, 0.800` at
  `k/√(gB)=0.1,0.3,0.5,1,2`. Hence `D_Sav → 0` and **`D/S → 0` uniformly in the deep IR** — the floor
  survives there with margin `S_bare`.
- **Crossover (`k ~ √(gB)`) `[O]`**: both `D₊` and `D_Sav` are `O(1)` ⇒ `D/S = O(1)`; whether
  `(★)` holds uniformly across cells/β **is open**. This is `H1` in the Savvidy-cell coordinate (softened).

## §4 · Verdict

- **Exact bound constructed `[T]`:** `D ≤ E[|⟨φ₋|c⟩|²]/(g|B|)`, and `(★)` `⟹ D<S ⟺ G″>0`. The Savvidy
  softening is the `D`-term, and `D<S` is the statement `E[|⟨φ₋|c⟩|²]/(g|B|) < S_bare`.
- **`D/S` at β=2.5 `[V-num]`:** `O(1)` — `0.3–0.6` for `κ` comfortably above `soft`, but **`>1`** (floor
  lost) when `g²Var(φ²)` approaches or exceeds `κ`. **NOT uniformly `<1`.**
- **Uniform `D<S`: NOT established `[O]/[NG]`.** Deep IR holds (`D/S→0`); the crossover `O(1)` inequality
  is `H1`. The Savvidy cell thus does not destroy the floor categorically — it converts `H1` into the
  precise inequality `(★)`, whose uniform form remains open.

*Witness:* `w7/savvidy_dterm.py` (numpy, `py -3.12 -I`, runs clean). *Refs:* `CONSTRUCT-SAVVIDY-CELL.md`
(softening), `SLOT-SD-FIXPOINT.md` (`G″=S−D`), `INVENT-VARIANCE-INEQ.md` (BL sign `D≤0 ⇔ K⪰0`),
`W6-BE` (high-momentum negative weight). **Clay 0%; no QED.**
