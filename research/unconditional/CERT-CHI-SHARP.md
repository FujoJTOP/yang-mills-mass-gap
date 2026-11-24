# CERT-CHI-SHARP — the χ-half certificate `μ(1−gBα) ≥ λ` evaluated: LHS is `O(gB)`, ranging over `[0, (2/3)gB]` with overlap strength; `λ=1/χ` is the same `O(1)` scale at β=2.5 ⇒ **exactly marginal, not closeable**

> **VERDICT (no QED, Clay 0%).** The χ-half of H1 is `μ(1−gBα) ≥ λ` (BLAST-OPERATOR-INEQ /
> BLAST-YM-NONLINEAR), where `K=H₊−cφ̂φ̂ᵀ`, `c=gB`, `μ=λ_min(H₊)`, `α=φ̂ᵀH₊⁻¹φ̂`, `λ=1/χ`. **[V-num]**
> the LHS is `O(gB)` (ranging over `[0,(2/3)gB]` as the negative-mode overlap grows) and the RHS `λ=1/χ`
> is the **same `O(1)` scale** at β=2.5 — **no hierarchy, no free pass.** The two sides **cross to 0** in the
> fully-aligned limit. **The certificate is exactly marginal; it cannot be closed by any `O(1)` number I can
> compute** — `μ, α, gB, λ` are all nonperturbative at β=2.5. `[O]` fixed step = the nonperturbative
> Savvidy-Hessian overlap at the crossover.

---

## §0 · Labels

**[T]** 真定理 · **[V-num]** new numeric · **[O]** 精确缺失步 · **[NG]** no-go. Witness
`w7/cert_chi_sharp.py` (runs clean, `py -3.12 -I`).

---

## §1 · The exact algebra

`Var(a_{k→0}) = ⟨v,K⁻¹v⟩ = ⟨v,H₊⁻¹v⟩ + c⟨v,H₊⁻¹φ̂⟩²/(1−cα) ≤ ⟨v,H₊⁻¹v⟩/(1−cα) ≤ 1/(μ(1−cα))`,
valid iff `cα<1`. So `H1 (Var ≤ 1/λ) ⟸ μ(1−gBα) ≥ λ`. The certificate is the **product** `μ·(1−gBα)`.

---

## §2 · The exact numbers (Savvidy Landau model, `[T-model]`)

One-loop Landau levels `λ_{n,s}=(2n+1)gB+2gB·s`. Sorted (×`gB`): `−1, +1, +3, +3, +5, +5, …`.
So `c=gB` (the negative mode), `μ=gB` (the `(1,−1)` mode). The negative mode (`n=0`) overlaps a positive
mode **only** at `n=0`, i.e. the `(0,+1)` mode at `3gB`:

| overlap model | `α` | `cα` | **LHS `μ(1−cα)`** |
|---|---|---|---|
| `n=0`-overlap (physical) | `1/(3gB)` | `1/3` | **`(2/3) gB`** |
| fully degenerate `H₊=μ·I` | `1/gB` | `1` | **`0`** |

As the overlap strength `t:=⟨neg|μ-mode⟩∈[0,1]` grows, **LHS `= gB·(1−t²)` sweeps `gB → 0`.** At `t→1`
(negative mode aligns with the lightest positive mode) the LHS is **exactly 0** — the negative mode eats the
*entire* floor. `SU(3)` has the same per-mode `(2n+1+2s)` structure (color only adds degeneracy) ⇒ same range.

---

## §3 · The `λ` side and the crossing

`λ = 1/χ = G″` measured on the lattice (corpus, **finite volume**): `G″ ≈ 0.28` (L=4), `0.04` (L=8). In
lattice units both LHS `= O(gB)` and `λ = 1/χ` are `O(1)` — **the same order**. There is **no scale
hierarchy**. The sides **cross**: LHS `→0<λ` at full overlap (`t→1`); for partial overlap LHS `=O(gB)` and
must **beat** `λ=1/χ` — a genuine `O(1)` comparison decided by the **nonperturbative overlap `t`** and the
**Savvidy scale `gB`**, neither computable exactly at β=2.5.

---

## §4 · Verdict

- **`[NG]` the certificate is NOT closeable by an exact number.** `μ, α, gB` (the Savvidy Hessian) and
  `λ=1/χ` are all nonperturbative at β=2.5; the cleanest model gives LHS `∈[0,(2/3)gB]`, **the same `O(1)`
  scale as `1/χ`**.
- **`[O]` exact step:** the **nonperturbative overlap `t`** of the Savvidy negative mode with the lightest
  positive coarse mode, at the crossover — equivalently `D<S` (`SAVVIDY-DTERM` `(★)`) = H1.
- **`[V-num]` sharp finding:** at full alignment the floor is **exactly eaten to 0** — the χ-half is *exactly
  marginal*, which is the structural reason the empirical `D/S ≈ 0.94–1.00` (SAVVIDY-LATTICE-CELL) sits at
  the edge. The positive `G″=+0.28` (L=4) is finite-volume.

**No certificate. `Clay = 0%`; no QED.**
