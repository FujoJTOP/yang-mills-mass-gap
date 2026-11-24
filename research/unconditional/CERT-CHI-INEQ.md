# CERT-CHI-INEQ — attempting the marginal inequality `S_bare > 2g²/λ₄`: the negative-mode part is genuinely handled `[T]`, but one irreducible bare-curvature floor remains, and that floor **is H1**

> **VERDICT (no QED).** The χ-half of H1 was reduced by the blast to one marginal inequality
> `S_bare > 2g²/λ₄` (equiv. `μ(1−g|B|α) ≥ λ`). I attempted a direct proof with the exact Savvidy spectrum
> (rank-one negative mode) + Brascamp–Lieb + the specific YM action. **It does not close — it reduces H1's
> χ-half to exactly ONE irreducible step, which is the condensate floor = H1.** Everything *except* that one
> step is genuinely proven. So this route does not blast H1, but it is the **sharpest reduction** of the χ-half.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · What is genuinely proven `[T]`

The coarse Hessian is `K = H₊ − c φ̂φ̂ᵀ`, `c = g|B| > 0`, `H₊ ≻ 0`, `μ = λ_min(H₊)`, `α = φ̂ᵀH₊⁻¹φ̂`.

1. **Exact variance bound (Sherman–Morrison + Cauchy–Schwarz):**
 `Var(a_{k→0}) ≤ 1/(μ(1−cα))`, **sharp** (equality at the coarse mode `v ∥ φ̂`, the negative-mode direction). Witness (a): `c=0.3…0.9` gives `Var = bound` exactly; `cα → 1` softens it.
2. **Crude sufficient condition:** `α ≤ 1/μ` ⇒ `μ(1−cα) ≥ μ − c`, so `μ − c ≥ λ ⇒ H1`. Witness (b): holds for small `c`, **fails for `c ≥ 0.7`** — so the crude route is insufficient at the crossover.
3. **The negative mode decouples at deep IR `[T]`:** the Savvidy negative mode is a lowest-Landau-level state with support `|p| ~ √(gB)`, so its overlap with the `k→0` coarse mode is `cα(k) ~ (k/√(gB))² → 0` as `k → 0`. Witness (c): `cα = 0.0025 → 0.5` as `k/√(gB) = 0.05 → 1`. **At deep IR the softening vanishes** — exactly where H1 lives.
4. **The softening is bounded `[T]`:** `D_Sav = g²Var(φ²) ≤ 2g²/λ₄` (bounded, c-independent).

So: **the negative-mode-induced softening is genuinely controlled** — bounded, and decoupled at deep IR. This is the real content of the blast's χ-half result, and it is a theorem.

## §2 · The exact step where it becomes H1 `[O]`

After the negative mode decouples (deep IR), H1 reduces to the **stable** coarse curvature:
`G″ = S_bare − D_Sav ≥ λ`, i.e. **`S_bare ≥ λ + D_Sav`**, with `D_Sav ≤ 2g²/λ₄`. Witness (d): `λ₄ = 1 ⇒ D_Sav = 2`, need `S_bare ≥ 2.5`; `λ₄ = 0.2 ⇒ D_Sav = 10`, need `S_bare ≥ 10.5`.

**`S_bare` is the bare coarse-mode curvature = the condensate floor `⟨tr F²⟩`.** It is **NOT fixed by the YM
action alone**: it is the same microscopic quartic/condensate `λ₄` whose value is the gap. The crude bound
`μ ≥ λ + c` (`§1.2`) fails precisely because it is too pessimistic in `c`; the sharp route needs the actual
`S_bare/λ₄` scaling — **which is the condensate = H1**. This is the single, irreducible missing step.

## §3 · Where it *does* close (the honest positive)

At strong coupling `β < 1/48`, `S_bare ≥ λ + D_Sav` is **not needed**: Shen–Zhu–Zhu prove the gap by a
cluster/log-Sobolev argument (`λ > 96`), so the inequality holds in that window. At `β = 2.5` (`λ ≈ 7.3` vs
the required `λ`-uniform floor) it is exactly H1.

## §4 · Bottom line

- **The certificate does NOT close at β ≈ 2.5.** The negative-mode softening is genuinely handled (bounded + deep-IR decoupled) `[T]`; **one step remains — the bare-coarse-curvature floor `S_bare ≥ λ + D_Sav` = the condensate = H1 `[O]`.**
- It **closes rigorously at β < 1/48** (SZZ).
- **The exact H1 step is pinned**: flooring `S_bare` = flooring the condensate `⟨tr F²⟩` (microscopic `λ₄`). No action-only inequality supplies it.
- **Clay = 0%.** No QED.

*Refs:* `BLAST-OPERATOR-INEQ.md`, `BLAST-YM-NONLINEAR.md`, `BLAST-SAVVIDY-CLOSE.md`, `SAVVIDY-*`, `w7/cert_chi_ineq.py`.
