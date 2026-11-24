# SAVVIDY-MULTICELL-RG — the accumulated softening is IR-MARGINAL (log-log divergent, IR-dominated); sub-criticality is decided at the Landau pole, not perturbatively = `H1`

> **VERDICT (no QED).** `CONSTRUCT-SAVVIDY-CELL` showed the Savvidy negative mode does *not* kill the χ
> floor locally: it softens it by a bounded amount (`Γ″(0)=κ_bare − g²Var(φ²)`), leaving a **softened
> inequality** `κ_bare > ∫(dk/k)g²Var(φ²)` to be checked *uniformly across RG cells*. This file computes
> that accumulated softener with asymptotic freedom and finds it is **not** a benign finite correction:
> `S = C_RG·Var(φ²)·ln( ln(k_UV/Λ)/ln(k_IR/Λ) )` — a **log-log** object that is **finite for `k_IR>Λ` but
> grows without bound as `k_IR→Λ⁺`**, and is **IR-dominated** (almost all of it accrues in the last few
> cells near `Λ`). So `κ_bare > S` is **not generic**: for a natural `O(1)` `κ_bare`, `S = O(10–100)·Var`
> already overshoots. Sub-criticality is decided exactly at `k_IR~Λ` where `g²→∞` (Landau pole,
> perturbation theory fails) — i.e. it **is** `H1` at the IR endpoint. This *sharpens* SAVVIDY-CELL: "the
> floor survives, softened" is **not** a free pass; the softening is an IR-marginal bound whose convergence
> is the nonperturbative IR.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`. Companion: `CONSTRUCT-SAVVIDY-CELL.md`.

---

## §1 · The setup and the closed form `[T]`

One-loop asymptotic freedom in pure SU(N) (`g²=λ`, 't Hooft): `g²(k) = C_RG/ln(k/Λ)` with `C_RG=24π²/11`
for SU(2) (matches `lean/YMStatement/Conditional.lean`'s `C_RG`). The accumulated Savvidy softener over the
flow from UV to IR is

```
S = ∫_{k_IR}^{k_UV} (dk/k) g²(k) Var(φ²) = C_RG · Var(φ²) · ∫ du/u = C_RG · Var · ln( ln(k_UV/Λ)/ln(k_IR/Λ) ),   u=ln(k/Λ).
```

Witness confirms closed form vs numeric to `~10⁻⁷` (single process, numpy). The one-loop `Var(φ²)` is
finite (set by the quartic; `0.41…8.51` in the SAVVIDY-CELL witness, not divergent as the negative curvature
grows) — so the *entire* accumulation structure is in the `ln ln` factor.

---

## §2 · It is IR-MARGINAL and IR-DOMINATED `[T]/[V-num]`

Fixed `k_UV=100Λ`, `Var=1`, `C_RG≈21.53`:

| `k_IR/Λ` | `ln(k_IR/Λ)` | `S` |
|---|---|---|
| 10 | 2.303 | 14.93 |
| 5 | 1.609 | 22.64 |
| 3 | 1.099 | 30.86 |
| 2 | 0.693 | 40.78 |
| 1.5 | 0.406 | 52.33 |
| 1.2 | 0.182 | 69.54 |
| 1.1 | 0.095 | 83.50 |
| 1.01 | 0.010 | 132.16 |

`S` is finite for every `k_IR>Λ`, but **grows without bound as `k_IR→Λ⁺`** (like `ln(1/ln k_IR)`), i.e. the
softening is **IR-marginal**. That it is **IR-dominated** follows directly from `dS/d ln k = g²(k)Var ∝
1/ln(k/Λ)`: per-decade it is `5.5` at `50Λ`, `9.4` at `10Λ`, `19.6` at `3Λ`, `53.1` at `1.5Λ` — the last
few cells near `Λ` carry almost all of `S`. (This is the integral face of W6-BE's "90–97% of the negative
weight is high-momentum" in the *coupling*, inverted: the softening weight is low-momentum.)

---

## §3 · Sub-criticality is NOT generic `[V-num]`

With `Var` in the SAVVIDY-CELL range and a natural `κ_bare ~ O(1)` (IR units), the margin `κ_bare − S`:

| `Var` | `S(k_IR=2Λ)` | `S(k_IR=1.1Λ)` | margin at `κ_bare=1` |
|---|---|---|---|
| 0.41 | 16.7 | 34.2 | −15.7 / −33.2 |
| 1.00 | 40.8 | 83.5 | −39.8 / −82.5 |
| 8.51 | 347.0 | 710.6 | −346 / −710 |

So `S = O(10–100)·Var`, **far above** a natural `O(1)` `κ_bare`; the margin is **negative**. Sub-criticality
`κ_bare > S` therefore needs `κ_bare` to be large enough to beat the accumulated `ln ln` — it is not a
generic perturbative fact.

---

## §4 · The exact step `[O]`

`κ_bare > C_RG · Var · ln( ln(k_UV/Λ)/ln(k_IR/Λ) )` must hold **as `k_IR→Λ`**, but there `g²(k)→∞` (the
Landau pole) and one-loop perturbation theory is invalid. So the convergence of the accumulated softener
is decided **exactly at the nonperturbative IR endpoint `k_IR~Λ`** — there is no perturbative cell-by-cell
argument. This is `H1` in the RG-flow coordinate (the IR endpoint).

**Caveat (normalization):** the comparison `κ_bare` vs `S` is in the same RG-normalized curvature units;
the *structural* result (log-log, IR-marginal, IR-dominated, decided at `k_IR~Λ`) is normalization-independent.

---

## §5 · Honest bottom line

- **The accumulation does NOT stay sub-critical automatically.** It is a **log-log**, **IR-marginal**,
  **IR-dominated** softener that diverges at the Landau pole; for natural `κ_bare~O(1)` it already overshoots.
- **This sharpens SAVVIDY-CELL, not reverses it:** the negative mode does not *destroy* the floor by scale
  separation (that part is `[T]`), but the surviving inequality is an **IR-marginal bound decided at
  `k_IR~Λ`** = `H1` — not a benign finite correction.
- **No new lever.** The multi-cell accumulation reduces, once more, to a nonperturbative IR statement =
  `H1`. **Clay 0%; no QED.**

*Companion:* `CONSTRUCT-SAVVIDY-CELL.md` (local softening, scale separation), `CONSTRUCT-RG.md` (the flow),
`CONSTRUCT-SYNTHESIS.md`. Witness: `w7/savvidy_multicell_rg.py`.
