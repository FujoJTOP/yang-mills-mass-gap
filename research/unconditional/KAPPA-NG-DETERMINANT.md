# KAPPA-NG-DETERMINANT — the exact Nambu–Goto determinant does **not** yield the stabilizing rigidity: `κ_NG = −1/(6π) ≈ −0.053` (wrong sign); the needed `κ≈0.19` is microscopic/UV, so `κ` stays `H1`

> **VERDICT (no QED).** The construction flood left one input `κ` in the rigid-string mass
> `m_{0⁺⁺}/√σ = 2π√(2κ)`, fitted to `[0.138,0.193]`, with `κ` = the coefficient of the extrinsic-curvature
> term `(κ/2)∮K² ds` in `E(ρ)=2πσρ + πκ/ρ`. **Question: does the exact NG / Lüscher–Weisz determinant
> supply this κ?** Computed: **No.** Tree-level NG is the free 2D CFT (no `k⁴`, `κ_tree=0`); the exact
> 4-derivative NG quartic is `V₄ = −(1/8)(A−B)² − (1/2)C²` (light-cone `−2R² − ½(P−Q)²`), and its 1-loop
> (Lüscher) effect on a closed ring contributes a **negative** `1/ρ` term `−(D−2)/(12ρ) = −1/(6ρ)`.
> In the MASTERFIELD normalization this is an effective `κ_NG = −1/(6π) ≈ −0.053` — **wrong sign**,
> **not in `[0.14,0.19]`**. The NG correction *destabilises*; the `+κ` that stabilises the neutral ring is
> an **extra, microscopic (flux-tube-core) irreducible coupling** — an **irrelevant dim-4 operator**
> (`CONSTRUCT-WORLDSHEET-CFT` `[T]`), whose value is set by the UV completion, not the NG IR expansion.
> **`κ` is not derived; it stays `= H1` (or needs a microscopic core computation).**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** numeric witness. Units `σ=1`, `√σ=1`. Witness: `w7/kappa_ng_determinant.py` (runs `py -3.12 -I`, numpy+sympy).

---

## §1 · (A) Tree-level NG is the free 2D CFT ⇒ `κ_tree = 0` `[T]`

In static gauge the NG action `σ∫√det(η_ab + ∂_aX^i∂_bX^i)` has the quadratic part `(σ/2)∫(∂_aX^i)²`:
a free 2D CFT of `D−2` bosons. The quadratic kernel is `ω² + k²` — **linear in `k²`, coefficient of `k⁴`
identically 0**. So the rigidity is **not a tree-level NG effect**: `κ_tree = 0`.

## §2 · (B) The exact 4-derivative NG quartic `[T]`

Expanding `√det` to quartic order (sympy, exact): with `A=∂₀X·∂₀X`, `B=∂₁X·∂₁X`, `C=∂₀X·∂₁X`,
```
V₄ = (1/2)(AB − C²) − (1/8)(A+B)² = −(1/8)(A−B)² − (1/2)C²
   = −2R² − (1/2)(P−Q)²      (light-cone: R=∂₊X·∂₋X, P=∂₊X², Q=∂₋X²)
```
(verified: expansion − claim = 0). This is the standard **Lüscher–Weisz quartic**; it is a *negative-definite*
combination, **not** the positive `(∂²X)²` rigidity that `E(ρ)` needs.

## §3 · (C) The 1-loop (Lüscher) `1/ρ` coefficient is NEGATIVE `[T]/[V-num]`

The closed-string Lüscher term is `E_L = −π(D−2)/(6L)`; for a ring `L = 2πρ` ⇒ `E_L = −(D−2)/(12ρ)`.
For `D=4`: coefficient of `1/ρ` is **`−1/6 = −0.1667`** (the `π`s cancel). In the MASTERFIELD
normalization (stabiliser `+πκ/ρ`), the NG "effective rigidity" is
```
κ_NG = (−1/6)/π = −1/(6π) = −0.0531  < 0   (wrong sign; in [0.138,0.193]? NO).
```

## §4 · (D) Minimum existence ⇒ NG alone cannot stabilise `[T]/[V-num]`

`E(ρ)=2πσρ + c/ρ` with `c = πκ − 1/6` has an interior minimum iff `c>0`. With `σ=1`:

| κ | c = πκ − 1/6 | interior min? | m/√σ |
|---|---|---|---|
| 0 (pure NG) | −0.167 | **no** | — |
| −0.053 (κ_NG) | −0.333 | **no** | — |
| 0.138 | +0.267 | yes | 2.590 |
| 0.193 | +0.440 | yes | 3.324 |

Pure NG+Lüscher `E(ρ)` is **strictly increasing** (witness) ⇒ **no interior minimum** — reproducing the
`[T]` "pure NG has no neutral closed-string mass" (`PROOF-EFFECTIVE-STRING`). A positive `1/ρ` coefficient
(a genuine `κ>1/(6π)`) is required; NG supplies the **opposite** sign.

---

## §5 · Conclusion: `κ` is microscopic/UV, not NG

- **NG does not generate κ.** It generates `κ_tree=0` at tree level and a **negative** `1/ρ` (destabilising)
  at 1-loop; its "effective κ" is `−1/(6π)`.
- **The `+κ≈0.19` is extra.** It is the flux-tube **core** extrinsic-curvature coupling — an **irrelevant
  dim-4 operator** (`(∂²X)²`, worldsheet scaling dim 4 > 2, `CONSTRUCT-WORLDSHEET-CFT` `[T]`) whose
  coefficient is set by the **UV completion / tube thickness**, not by the NG IR determinant.
- **Sharpening:** not only does NG fail to give κ — it gives the **wrong sign**. The microscopic rigidity
  must both *supply* `+κ` and *overcome* the NG-induced `−1/6` (so the required micro-κ is even a bit
  larger than the pure fit; e.g. with Lüscher included, `κ=0.193` gives `m/√σ=3.32`, not 3.90 — the fit
  shifts **up**, not down).
- **Honest caveat:** "κ as the induced `(∂²X)²` Wilson coefficient at a cutoff" is scheme-dependent
  (irrelevant operator); any such value is UV-set and non-universal. This *reinforces* the conclusion.

**Exact missing input** `[O]`: the microscopic flux-tube-core rigidity (the transverse profile / thickness
of the confining tube) — equivalently `Λ`/the confinement scale at the core — `= H1` in string coordinates.
**No NG derivation of κ exists; Clay 0%; no QED.**

*Refs:* `CONSTRUCT-MASTERFIELD.md`, `CONSTRUCT-WORLDSHEET-CFT.md`, `CONSTRUCT-KAPPA-DERIVE.md`,
`CONSTRUCT-WIDTH-KAPPA-LATTICE.md`, `CONSTRUCT-LUSCHER.md`, `PROOF-EFFECTIVE-STRING`,
`TRANSMUTE-STRING-SPECTRUM`. Witness: `w7/kappa_ng_determinant.py`.
