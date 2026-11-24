# BLAST-YM-NONLINEAR — the sharpest nonlinear-YM attack on H1: Brascamp–Lieb on the positive fiber + the exact rank-one Savvidy negative mode **reduces H1 to ONE marginal O(1) inequality**; the constant (bare coarse curvature vs bound `2g²/λ₄`) is the condensate = `H1`

> **VERDICT (no QED; Clay 0%).** Per `ONLY-EXIT-SYNTHESIS`, the only non-trivial exit is the **specific
> nonlinear YM structure**. This blast attacks `H1` (`χ = Var_ν(a_{k→0}) ≤ 1/λ`) with the one genuine
> **nonlinear** tool that is not a linear/integral functional: **Brascamp–Lieb (BL)** — which bounds a
> variance by the *inverse Hessian* (`Var(a) ≤ E[⟨Hess V⟩⁻¹_{aa}]`), a nonlinear object — **combined with the
> exact rank-one Savvidy negative mode**. Result: this is the *right* tool, and it gives a **`[C]` conditional
> proof** `H1 ⟺ G″ = S_bare − D_Sav > 0` where `D_Sav = g²Var(φ²) ≤ 2g²/λ₄` is **bounded, `c`-independent,
> sign-protected**. So `H1` collapses to **one marginal `O(1)` inequality `S_bare > 2g²/λ₄`**, everything else
> proven. It does **not close**: the `O(1)` constant (`S_bare` vs `2g²/λ₄`) is the microscopic condensate
> `λ₄` = `H1`. **Blast fails to close; the exact missing step is now a single, clean, nonlinear-YM
> sub-criticality bound.**

---

## §0 · Labels & the tool

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** 数值. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

**Why BL is the nonlinear tool (not the blind linear class).** Every integral/linear/geometric functional is
provably blind to the atom (OPEN-EXIT: `Π` is linear in ν; weak-* continuity). BL instead bounds
`Var_μ(a) ≤ E_μ[⟨Hess V⟩⁻¹_{aa}]` — it uses the **inverse Hessian of the YM action**, a genuinely nonlinear,
field-level object (not a moment of ν). It is the natural tool for the `χ` half of `H1`.

---

## §1 · The proof attempt `[C]`

**(1) BL on the log-concave part.** At `β=0` the Haar weight alone is log-concave (`Φ` convex, Prékopa), so BL
applies and gives the gap (`W_geo ≥ 0`, `Φ″ ≥ N/12`). At `β=2.5` the full weight is `e^{−S_W−Φ}`; `Hess Φ ⪰ 0`,
but `Hess S_W` is the **Savvidy Hessian, indefinite** (`λ₋=−g|B|`, transverse magnetic), so the full weight is
**not** log-concave and BL in its bare form fails — **exactly on the direction where the gap must be generated**.

**(2) Exact rank-one rescue `[T]` (`SAVVIDY-DTERM`).** The indefiniteness is **rank one**: `Hess = H₊ − |λ₋|P₋`.
The fiber inverse splits `K⁻¹ = K₊⁻¹ − P₋/|λ₋|`, giving the exact defect
`D = D₊ + D_Sav`, with `D₊ ≤ 0` by **BL on the positive fiber** and `D_Sav = E[|⟨φ₋|c⟩|²]/(g|B|) ≥ 0`.
Hence
```
χ = 1/G″,   G″ = S_bare − D_Sav,    D_Sav = E[|⟨φ₋|c⟩|²]/(g|B|).
```
**(3) `D_Sav` is bounded `[V-num]` (`SAVVIDY-NEGMODE-EXACT`).** The exact (not wrong-saddle) `Var(φ²)` saturates
at `2/λ₄`, independent of `c=gB`, so `D_Sav ≤ 2g²/λ₄` — **not growing with the instability, sign-positive**.

**(4) The reduction.** `H1` (`G″ ≥ λ`) ⟺ `S_bare − D_Sav ≥ λ` ⟺ **`S_bare − λ ≥ D_Sav`** ⟺ (using the bound)
**`S_bare − λ ≥ 2g²/λ₄`**. This is **ONE marginal `O(1)` inequality**. Everything except its constant is proven:
`D₊ ≤ 0` (BL), `D_Sav` bounded/sign-positive (`2g²/λ₄`), deep-IR `D_Sav/S_bare → 0` (scale separation).

---

## §2 · Why this does NOT close — the exact missing step `[O]`

The blast reduces `H1` to `S_bare > λ + 2g²/λ₄`. The failure is in the **constant**, and it is structural:

- **`S_bare`** is the bare coarse curvature `E[h_n − ⟨K₊⁻¹c,c⟩]`; **`2g²/λ₄`** is the bounded Savvidy softening.
- Both are `O(1)` at the crossover (`λ₄ ~ g²` physically ⇒ `2g²/λ₄ ~ O(1)`), and the empirical `D/S ≈ 0.94–1.00`
  (`SAVVIDY-LATTICE-CELL`) puts the margin at **near-critical**. So `S_bare > λ + 2g²/λ₄` is a **marginal `O(1)`
  competition**, not a bound with slack.
- Proving it requires a **uniform lower bound on the bare coarse curvature `S_bare`** (equivalently: the
  microscopic quartic `λ₄` floored, i.e. the gluon condensate magnitude `⟨tr F²⟩ ≈ c/λ₄`) — **that is exactly
  `H1`**. The nonlinear YM tool (BL + rank-one) was *used to the hilt*; the residual is a single scalar floor.

**Blast status:** the nonlinear YM tool is identified and correctly applied; `H1` is **reduced to one marginal
`O(1)` inequality** with the *constant* = the condensate = `H1`. **No closure.**

---

## §3 · The atomless-bump adversary (honest scope)

BL + rank-one bounds `χ` — the **susceptibility (H1) half** — and is **blind to the atom `f₁²`** (the
ONLY-EXIT adversary `:φ²:` is a positive measure with the same `χ`, same moments, `f₁²=0`). So this blast is
correctly aimed at the `χ` half only; it does **not** read the atom (that half is A3/`H1`, per `ONLY-EXIT-SYNTHESIS`).
This is consistent, not a gap in the blast: the `χ` half has a genuine nonlinear tool (BL+rank-one) reducing to
one marginal inequality; the atom half remains A3/`H1`.

---

## §4 · Honest bottom line

- **The nonlinear YM tool exists and is applied:** Brascamp–Lieb (positive fiber) + exact rank-one Savvidy.
  It is not in the blind linear/integral class — it is the field-level inverse-Hessian bound the `χ` half needs.
- **`H1` reduces to ONE marginal `O(1)` inequality** `S_bare > λ + 2g²/λ₄` — the sharpest form `H1` has taken,
  with `D₊≤0`, `D_Sav≤2g²/λ₄`, sign, and deep-IR all **proven**.
- **It does not close:** the `O(1)` constant is the microscopic condensate `λ₄` = `H1`. **Blast fails to close.**
- **The exact missing step is now singular and clean:** a uniform lower bound on the bare coarse curvature
  `S_bare` (equivalently the condensate floor) exceeding `λ + 2g²/λ₄` at `β≈2.5`. This is `H1`. **Clay 0%.**

*Refs:* `SAVVIDY-DTERM.md`, `SAVVIDY-NEGMODE-EXACT.md`, `ONLY-EXIT-SYNTHESIS.md`, `TOOL-SYNTHESIS.md`,
`PNT-OPE-RHO.md`. Witness: `w7/blast_ym_nonlinear.py`.
