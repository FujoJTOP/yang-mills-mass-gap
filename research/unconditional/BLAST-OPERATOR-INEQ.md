# BLAST-OPERATOR-INEQ — the direct operator inequality EXISTS (Sherman–Morrison + Cauchy–Schwarz) but is SOFTENED by the rank-one negative mode: `Var(a_{k→0}) ≤ 1/(μ(1−g|B|α))`; the sharp `K⁻¹ ≥ λ` fails exactly as `g|B|α → 1`

> **VERDICT.** A direct operator inequality on the coarse-mode variance **does exist** and is clean:
> for the Savvidy/D-term Hessian `K = H₊ − c φ̂φ̂ᵀ` (`c = g|B| > 0`, `H₊ ≻ 0`, rank-one **negative** mode
> along `φ̂`), Sherman–Morrison + Cauchy–Schwarz give exactly
> `⟨v, K⁻¹v⟩ ≤ ⟨v, H₊⁻¹v⟩/(1 − cα) ≤ 1/(μ(1 − cα))`, `α = φ̂ᵀH₊⁻¹φ̂`, `μ = λ_min(H₊)`, valid iff `K ≻ 0`
> (⟺ `cα < 1`). **But it is SOFTENED** by the negative-mode factor `(1 − cα)⁻¹`, and it is **sharp in the
> negative-mode direction** `v ∥ φ̂` (Cauchy–Schwarz equality) — so the softening is *real*, not a slack loss.
> **The rank-one negative mode does NOT allow a sharp `K⁻¹ ≥ λ`**: the bound diverges as `g|B|α → 1`
> (negative-mode instability). H1 becomes the **softened** condition `μ(1 − g|B|α) ≥ λ` — the negative-mode
> sub-criticality, = `SAVVIDY-DTERM`'s `(★)`. No sharp blast; the exact obstruction is pinned.

---

## §0 · Setup

Coarse mode `v` (`|v|=1`), `a = ⟨v, X⟩`, marginal variance `Var(a_{k→0}) = ⟨v, K⁻¹v⟩` for a Gaussian with
precision `K`. Savvidy/D-term Hessian: `K = H₊ − c φ̂φ̂ᵀ`, `c = g|B| > 0`, `H₊ ≻ 0` (positive part),
`φ̂` the rank-one **negative** direction (`|φ̂|=1`). `α := φ̂ᵀH₊⁻¹φ̂ > 0`. `K ≻ 0 ⟺ cα < 1` (normalizable).

## §1 · The inequality `[T]`

**Sherman–Morrison** (exact): `K⁻¹ = H₊⁻¹ + c (H₊⁻¹φ̂)(H₊⁻¹φ̂)ᵀ/(1 − cα)`. Hence for every `v`:
```
⟨v, K⁻¹v⟩ = ⟨v, H₊⁻¹v⟩ + c·⟨v, H₊⁻¹φ̂⟩² / (1 − cα).          (exact)
```
**Cauchy–Schwarz** in the `H₊⁻¹` inner product: `⟨v,H₊⁻¹φ̂⟩² ≤ ⟨v,H₊⁻¹v⟩·α`. Therefore
```
⟨v, K⁻¹v⟩ ≤ ⟨v, H₊⁻¹v⟩ · (1 + cα/(1−cα)) = ⟨v, H₊⁻¹v⟩ / (1 − cα) ≤ 1 / (μ(1 − cα)).     [T]
```
The Cauchy–Schwarz step is an **equality iff `v ∥ φ̂`** — i.e. the bound is **attained in the negative-mode
direction**, which is the *worst* direction. So the inequality is **sharp where it matters**; the softening
`(1−cα)⁻¹` is not slack.

## §2 · Corollary (H1 in the softened form)

`Var(a_{k→0}) ≤ 1/(μ(1 − g|B|α))`, so
```
H1 :  Var ≤ 1/λ   ⟺   μ·(1 − g|B|α) ≥ λ.                                              [C]
```
This is exactly `SAVVIDY-DTERM`'s softened sub-criticality `S_bare > g²Var(φ²)` — the rank-one negative-mode
term `g|B|α` is the softening.

## §3 · The exact obstruction `[O]`

- The rank-one negative mode **increases** `Var(a)` (the `+c⟨·⟩²/(1−cα)` term, `c=g|B|>0`).
- The upper bound is therefore **softened** by `(1 − g|B|α)⁻¹`, which **diverges as `g|B|α → 1⁺`**
  (the negative mode approaching instability: `K` loses positivity and `Var → ∞`). Witness: at fixed `H₊,φ̂`,
  the bound goes `1.28 → 3.00 → 8.99 → 89.9 → 899` as `cα = 0.3 → 0.7 → 0.9 → 0.99 → 0.999`.
- **Hence no sharp `K⁻¹ ≥ λ` follows from the positive part alone.** Any proof of `Var ≤ 1/λ` must supply the
  uniform bound `μ(1 − g|B|α) ≥ λ` across the crossover — i.e. the negative-mode overlap `g|B|α` staying
  uniformly sub-critical. That is precisely `H1` (softened) / the `D<S` rank-one condition. **The "blast"
  does not bypass the wall; it reproduces it as a single, sharp, softened inequality.**

**Sharpness note:** the bound is exact (not loose) in the dangerous direction `v ∥ φ̂`; the numerical
"tightness" ratio `<v,K⁻¹v>/((vᵀH₊⁻¹v)/(1−cα))` reaches `1` for `v ∥ φ̂` and is `<1` otherwise (witness:
`0.979` for `v ∥ H₊⁻¹φ̂`, the wrong vector for equality).

## §4 · Verdict

- **Direct operator inequality: EXISTS `[T]`** — Sherman–Morrison + Cauchy–Schwarz, `Var ≤ 1/(μ(1−g|B|α))`.
- **Sharp `K⁻¹ ≥ λ`: FAILS `[NG]`** — the rank-one negative mode softens it, with the bound sharp (and
  divergent) in the negative-mode direction.
- **Exact missing input `[O]`** — uniform `μ(1 − g|B|α) ≥ λ` at the crossover = negative-mode sub-criticality
  = `D<S` (`SAVVIDY-DTERM`) = `H1` in the Savvidy-cell coordinate.

**No QED. Clay 0%.** The blast yields a clean theorem + a sharp softened inequality; the wall stands as the
negative-mode sub-criticality.

*Witness:* `research/unconditional/w7/blast_operator_ineq.py` (`py -3.12 -I`, clean). *Refs:* `SAVVIDY-DTERM.md`,
`CONSTRUCT-SAVVIDY-CELL.md`, `SAVVIDY-VAR-PHI2.md`, `ONLY-EXIT-SYNTHESIS.md`.
