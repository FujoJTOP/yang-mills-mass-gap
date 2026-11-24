# IMPOSSIBLE-ENTANGLEMENT — the entanglement-gap base does **not** force the 0⁺⁺ gap without ν: the entanglement spectrum is an **EDGE / 1-form** object (Li–Haldane), Schur-orthogonal to the 0⁺⁺ bulk scalar, and a functional of the vacuum measure

> **VERDICT (no QED, Clay 0%).** The "impossible base" — *the gap = the entanglement gap of the YM vacuum
> (modular Hamiltonian / entanglement spectrum); prove it forces the 0⁺⁺ gap without going through the measure ν* —
> **does not force it.** Three exact obstructions, and a decisive witness. **(1) Wedge/Bisognano–Wichmann:**
> the modular Hamiltonian of a wedge is `2π·K_boost` — **continuous spectrum, no gap** `[T]`. **(2) Area law ⟹
> gap is FALSE:** the free massless chain is gapless (`phys gap→0`) with **log** entanglement entropy; Hastings'
> theorem is the **converse** (gap ⟹ area law) `[T]/[NG]`. **(3) The entanglement spectrum is the EDGE
> spectrum, not the bulk gap** (Li–Haldane): two SSH phases with the **same bulk gap** have **opposite**
> entanglement gap (witness: trivial `ent_gap=4.02` vs topological `ent_gap=0.00`) — it reads the boundary /
> topological data, i.e. the **1-form/center sector**, Schur-orthogonal to the 0⁺⁺ bulk scalar. And the
> entanglement spectrum is a **functional of the vacuum state = the measure μ**, so "without going through ν"
> is impossible — reading its gap *is* reading μ. **The base reduces to H1.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The three sub-bases, each fails

**(i) "The vacuum's entanglement Hamiltonian has a spectral gap."** For a wedge region the modular Hamiltonian
is the **boost generator** (Bisognano–Wichmann: `K = 2π K_boost`); the boost has **continuous spectrum** (the
modular operator `Δ = e^{−K}` has spectrum `[0,∞)`), hence **no gap** `[T]`. For a general region `K_A = −log ρ_A`
is **non-local** and is not the boost; its gap is not controlled by any general theorem, and for a free massless
field the modular spectrum is **gapless** `[T]`.

**(ii) "The area law forces a correlation gap."** FALSE `[NG]`. Hastings' theorem is the **converse** — a gap
implies an area law `[T]`; the reverse does **not** hold. The free massless field has an area law (with a log
coefficient) and is **gapless** (witness: 1D massless chain, `phys gap ≈ 0.026 → 0`, entropy `S` grows like `log l`).
A gapless system with an area law is a direct counterexample to "area law ⟹ gap."

**(iii) "The entanglement spectrum of the `tr F²` region is gapped ⟹ the 0⁺⁺ gap."** FALSE `[NG]`, two reasons:
- **Edge, not bulk (Li–Haldane) `[T]`.** The low-lying entanglement spectrum of a region tracks the **edge /
  boundary** spectrum, not the bulk gap. The witness makes this decisive: two SSH phases with the **identical
  bulk gap** `2|t₁−t₂| = 1.0` — the **trivial** phase (no edge state) has `ent_gap = 4.02` (gapped), the
  **topological** phase (two zero edge modes) has `ent_gap = 0.00` (gapless). **Same bulk gap, opposite
  entanglement gap** ⇒ the entanglement gap is *not* the bulk gap; it encodes boundary/topological data.
- **Functional of μ `[O]`.** The reduced density matrix `ρ_A` is the partial trace of the vacuum state, i.e. of
  the **YM measure μ**. The entanglement spectrum is a **functional of μ** — reading its gap is reading μ. So
  "without going through ν" is impossible; the entanglement route is a change of *representation* of the same
  object, not a new handle.

---

## §2 · The witness (`w7/impossible_entanglement.py`, exact, free-fermion correlation-matrix)

| phase | bulk gap `2|t₁−t₂|` | #edge modes | entanglement gap `min|ξ|` | entropy `S` |
|---|---|---|---|---|
| trivial `t₁=1.0,t₂=0.5` | 1.0 | 0 | **4.019** | 0.178 |
| **topological** `t₁=0.5,t₂=1.0` | **1.0** | **2** | **0.000** | 0.699 |
| massless `t₁=t₂=1.0` | `≈0.026→0` | — | 1.30 (→0 as `l→∞`) | 1.14→1.20 (**log**) |

- **Same bulk gap (1.0), opposite entanglement gap** (4.02 vs 0.00): the entanglement spectrum reads the
  **edge**, not the bulk — the precise statement of the Li–Haldane obstruction.
- **Massless:** `phys gap→0` *and* entropy grows like `log l` (not a strict area law) — and it is **gapless**:
  area law/log does not imply a gap, and the entanglement gap `1.30 ≠` the physical gap `0.026`.
- **`ent_gap ≠ phys_gap`** everywhere (4.02 vs 1.0; 1.30 vs 0.026): the entanglement gap is a *different*
  quantity.

---

## §3 · The two structural obstructions (why it must fail)

1. **Channel mismatch (edge = 1-form).** The entanglement spectrum is a **boundary / edge** object. In gauge
   theory the boundary carries the **1-form / center** sector, which is **Schur-orthogonal** to the bulk 0-form
   **0⁺⁺** scalar `tr F²` (the same mismatch that killed ŧ Hooft / center-vortex / defect routes,
   `BASE-ORDER-PARAMETER`, `ROOF-0FORM-DUAL`). So the entanglement spectrum targets the **wrong channel** for
   the 0⁺⁺ gap.
2. **It is a functional of the measure.** `ρ_A = Tr_{A^c} μ`; the entanglement spectrum is determined by μ. It
   does not provide an independent handle — it is μ in entanglement clothing.

---

## §4 · Acid test + verdict

- **Acid test:** the free massless field is **gapless** (`phys gap→0`). An "entanglement gap ⟹ physical gap"
  argument **misfires** on the topological SSH phase (entanglement-gapless, bulk-gapped) — it is **false as an
  implication**, hence acid-dirty. An argument that *does* respect the acid test must use extra structure that
  distinguishes YM from the free field — which is the measure μ = `H1`.
- **Verdict `[NG]`:** the entanglement-gap base does **not** force the 0⁺⁺ gap without ν. The exact failing
  steps: (i) wedge modular spectrum is continuous (no gap); (ii) area law does not imply a gap (Hastings
  converse); (iii) the entanglement spectrum is an edge/1-form object (Li–Haldane), Schur-orthogonal to the
  0⁺⁺ bulk scalar, and is a functional of μ — reading its gap is reading μ = `H1`.
- **Genuine yield (a teaching, not a lever):** the entanglement spectrum of the YM vacuum is a well-defined,
  interesting **new object** (modular spectrum), but it is an **edge/1-form** object and a **functional of the
  measure**; it cannot be the bulk 0⁺⁺ gap certifier. Consistent with `ONLY-EXIT`.

*Files:* `w7/impossible_entanglement.py` (runs clean, `py -3.12 -I`). Labels `[T]/[C]/[O]/[NG]/[V-num]`; NO QED;
**Clay 0%**.
