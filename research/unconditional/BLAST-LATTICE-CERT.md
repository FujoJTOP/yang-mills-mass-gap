# BLAST-LATTICE-CERT — a certified L-UNIFORM compactness bound `0 ≤ Var(a_{k→0}) ≤ 4π²`; it does NOT close H1, and it localizes the wall to the `a→0` (renormalization) step, NOT `L→∞`

> **VERDICT (no QED, Clay 0%).** A **rigorous** bound on the lattice susceptibility is obtained:
> `0 ≤ χ_lattice := Var(a_{k→0}) ≤ (2π)² ≈ 39.5`, from compactness + the bounded principal log
> (`|a| ≤ π`), **uniformly in `L`** `[T]`. It **certifies** that the finite-`L` variance is bounded by a
> **constant independent of `L`** — hence the **`L→∞` step is NOT the obstruction** for this route (a genuine
> sharpening). It does **not** close H1: the ceiling is ~**9.5×** above the measured value (`3.22–4.17` at
> `β=2.4, L=8`), far above the `G″`-scale; and critically it is a **lattice-units** bound — the **physical**
> `χ ≤ 1/λ` needs the lattice→continuum map `Z(a)`, which compactness does **not** control. **The wall is
> exactly the `a→0` (dimensional-transmutation) step.** Tightening it (to the needed scale, or to a positivity
> certificate of the *exact* `G″`, not the Gaussian proxy) requires the specific YM measure — which by the
> campaign's weak-\* discontinuity result no continuous/bounded-kernel functional provides = `H1`.

---

## §0 · Labels
**[T]** 真定理 · **[V-num]** measured · **[O]** exact missing step · **[NG]** no-go. Witness: `w7/blast_lattice_cert.py`.

## §1 · The certified bound `[T]`
`a_{k→0} = ⟨v, X⟩`, `X` the branch-safe principal log (`SU(2)` angles `θ∈[0,π]`, centre-twist corrected),
`v` a unit mode ⇒ `|X|_op ≤ π ⇒ |a| ≤ π`. Hence, for the same measure,
`Var(a) ≤ (range a)² ≤ (2π)² = 4π² ≈ 39.5`. **This holds for every `L`** (the compact-group diameter does
not grow with volume). With `Var ≥ 0`: the **certified interval is `[0, 39.5]`, `L`-uniform.**

## §2 · Measured vs ceiling `[V-num]`
`SU(2)`, `β=2.4`, `L=8`, `N=400`, 12 transverse modes: `Var(a) ∈ [3.22, 4.17]` (mean `3.64`),
`1/Var(a) ∈ [0.240, 0.310]` (matches the corpus `1/Var≈0.27`). The ceiling exceeds the max by **9.47×**.
The rigorous-style **positivity** interval from these samples is `[−0.035, …]` (`G″` one-sided), **does not
exclude 0** — so no certificate of `G″>0`. The Gaussian-proxy floor `G″ ≥ 1/(4π²) ≈ 0.025` is **weak and
model-dependent** (`G″=1/Var` only in the Gaussian proxy; the exact `G″` is the curvature of `−log P(a)`, not
`1/Var`).

## §3 · Limit-localization — the round's sharpening
- **`L→∞` is NOT the wall** `[T]`: the compactness ceiling `4π²` is independent of `L`.
- **`a→0` IS the wall** `[O]`: H1 is the *physical* `χ ≤ 1/λ`; the ceiling is in **lattice units**, and the
  lattice→physical conversion `Z(a)` is uncontrolled by compactness. As `a→0`, `χ_phys` is set by `Z(a)`, not
  by the lattice constant. So the certificate is `L`-uniform but `a`-useless.

## §4 · Honest bottom line
- **A rigorous L-uniform lattice bound `χ_lattice ≤ 4π²` is certified** — and it proves **`L→∞` is not the
  obstruction**, only the **`a→0`/renormalization** step is. That is a real (if modest) sharpening: it tells us
  *where* a certificate must attack (the continuum map, not the volume).
- **H1 is not certified**: the bound is ~1 order too loose at the lattice scale, and a-uncontrolled; tightening
  it to a positivity certificate of the exact physical `G″` needs the action-specific YM structure = `H1`.
- **Clay = 0%.** No QED. This is an honest certification attempt that pins the limit: **the wall is `a→0`
  (dimensional transmutation), not `L→∞`.**
