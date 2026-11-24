# CERT-CHI-VARIATIONAL — a genuinely action-specific variational bound on `Var(a_{k→0})` exists (diagonal = single-link variance, computable, O(1)) but its off-diagonal part is *exactly* the gap = H1; no route reaches `Var ≤ 1/λ` at β≈2.5

> **VERDICT (no QED, Clay 0%).** Constructed the sharpest **action-specific** trial-function bound on the
> coarse-mode variance `Var(a_{k→0}) = χ` (the χ half of H1). It splits **exactly** into a **diagonal**
> (single-link / computable) term and an **off-diagonal** (connected two-point) term. The diagonal is a genuine
> action-specific `O(1)` number — computed here for SU(2), `σ²_link(β): 0.322(Haar) → 0.102(β=2.5)` — but it is
> **insufficient**. The off-diagonal is the *entire* difficulty: its sharpest bound is the spectral (KL) one
> `χ ≤ Var(X)/Δ`, which **needs `Δ = E₁−E₀ > 0` = the gap = H1**. Every available action-specific inequality
> that could bound the off-diagonal fails for a structural reason (Gaussian domination → Wilson loops not the
> log-link; FKG → absent for the full-rep character sum; Brascamp–Lieb → softened by the Savvidy rank-one
> negative mode). **The bound reduces to H1.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The object and the exact split `[T]`

For `a = Σ_x v_x X_x` (`X` = su(2) log-link, `Σv_x²=1`):

```
Var(a) = Σ_x v_x² · σ²_link                        (DIAGONAL: single-link, computable)
       + Σ_{x≠y} v_x v_y · Cov(X_x, X_y)           (OFF-DIAGONAL: connected)
```

This is an identity. The diagonal is a genuine **action-specific** quantity: the single-link marginal under the
Wilson weight `e^{β Re tr U}`. The off-diagonal is the connected two-point sum.

## §2 · The diagonal `[V-num]` (real, computable, insufficient)

SU(2), marginal `∝ sin²α · e^{2β cos α}`, `α = class angle ∈[0,π]` (witness `w7/cert_chi_variational.py`):

| β | 0 (Haar) | 0.5 | 1.0 | 1.5 | 2.5 |
|---|---|---|---|---|---|
| `Var(α) = σ²_link` | 0.322467 | 0.292384 | 0.229033 | 0.170698 | 0.102208 |

Haar check exact: `σ²_link(0) = π²/12 − 1/2 = 0.322467`. The diagonal is a real, falsifiable, action-specific
`O(1)` term that **drops** as β grows (ordering/confinement). But it is a *diagonal* `O(1)/V_eff` piece: it does
not grow toward the weak-coupling (gapless) direction, so it cannot by itself carry the crossover.

## §3 · The off-diagonal is H1 `[T]/[O]`

The sharpest bound on the off-diagonal is the spectral/KL representation `χ = Σ_n |⟨0|X|n⟩|²/(E_n−E₀)`, giving

> **`χ ≤ Var(X)/Δ`, `Δ = E₁−E₀`** — valid **iff `Δ > 0`**; at `Δ=0` the bound is `+∞`. (Toy KL check exact for
> `Δ = 0.1,0.5,1`: `χ = 8,1.6,0.8 ≤ Var(X)/Δ = 10,2,1`.)

So the off-diagonal bound **is** the gap. This is the same squeeze-with-gap as `BLAST-TRANSFER` / `PNT-CHI-CERTIFICATE`.

## §4 · No action-specific inequality closes the off-diagonal `[NG]`

- **Gaussian domination** (Fröhlich–Simon–Spencer) bounds `E[W(C)W(C′)]` for **Wilson-loop characters**, not for the
  log-link `X`; it does not control `Cov(X_x,X_y)`.
- **FKG / positive association** needs a ferromagnetic (increasing) coupling; the full-representation SU(N) Wilson
  weight is an alternating **character sum** ⇒ **no FKG** ⇒ no sign control on the off-diagonal.
- **Brascamp–Lieb** needs log-concavity; the Savvidy Hessian is **indefinite** (rank-one negative) at β≈2.5, so BL
  survives only on the positive fiber, **softened** by `(1−g|B|α)^{-1}` (`BLAST-OPERATOR-INEQ`).

## §5 · Verdict

The sharpest **action-specific variable/action-specific** bound = [diagonal `σ²_link(β)` (real, computable, O(1))] +
[off-diagonal `≤ Var(X)/Δ` (KL, needs the gap)]. The diagonal is genuine but insufficient; the off-diagonal is
**exactly H1**. **No route reaches `Var ≤ 1/λ` at β≈2.5 by this variational decomposition** — it reduces to H1.
This is consistent with `ONLY-EXIT` (the χ half collapses to one marginal inequality whose constant is the gap).

*Refs:* `BLAST-TRANSFER.md`, `BLAST-OPERATOR-INEQ.md`, `BLAST-YM-NONLINEAR.md`, `PNT-CHI-CERTIFICATE.md`,
`ONLY-EXIT-SYNTHESIS.md`; witness `w7/cert_chi_variational.py`.
