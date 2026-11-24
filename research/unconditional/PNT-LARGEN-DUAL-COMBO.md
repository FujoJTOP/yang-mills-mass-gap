# PNT-LARGEN-DUAL-COMBO — no dual master field exists; the bridge breaks at large-N Abelian projection; the combined object's structure and predictions are pinned

> **VERDICT (no QED).** The two renormalon-free routes — (A) the large-N master field (spectral pole, `1/N`)
> and (B) dual superconductivity (dual-Higgs mass from a dim-2 condensate) — **do not meet in any
> established "dual master field."** They live on structurally different objects: the large-N master field
> is the **gauge-loop** master field `W(C)` of the Makeenko–Migdal equation; the monopole condensate `v` is
> an **effective dual-sector** object that requires 't Hooft Abelian projection, which at large `N` is
> gauge-dependent and uncontrolled (the monopole dominance behind dual superconductivity is a lattice
> observation, not a theorem). The bridge would need a **large-N controlled Abelian projection** plus a
> proof that `v` is an order parameter of the master field — neither exists. What *is* true and pinned:
> `v` (VEV, dim 1) `~ Λ`, the condensate `v²` (dim 2) `~ Λ²` by dimensional transmutation, both `N`-independent
> under 't Hooft scaling; the dim-2 Borel pole sits at `u=1` (vs `u=2` for `tr F²`) — structurally cleaner;
> and a combined object would land in the consistency band `m_{0⁺⁺}/√σ ≈ 3.3–3.9`. But the "prediction" is a
> **one-parameter fit** (`κ_GL` free), because the combined object's job would be exactly to *fix* `κ_GL` from
> the master field — which is the missing bridge. **Clay 0%.**

---

## §1 · The bridge question

- **Route A** gives `m_{0⁺⁺}(λ)` as a physical spectral pole with a `1/N` handle; its obstruction is
  *4D master-field existence / closing the Makeenko–Migdal hierarchy*.
- **Route B** gives a renormalon-free 0⁺⁺ mass from the dual Abelian–Higgs curvature, `m_H=√(2λ)v`,
  because the dim-2 condensate has a Borel pole at `u=1` (not `u=2`); its obstruction is *monopole
  condensation = confinement = `H1`*.

A combined "dual master field" would inherit **both** the `1/N` handle and the `u=1` Borel structure.
Question: does such an object exist, and does `v` appear as the master-field gap?

---

## §2 · What is true `[T]/[V-num]`

- **Dimensional transmutation.** The dual VEV `v` (dim 1) `~ Λ`; the condensate `v²` (dim 2) `~ Λ²`. Under
  't Hooft scaling (`λ=g²N` fixed, `Λ` `N`-independent) both are `N`-independent. So `v ~ Λ` (VEV) /
  `v² ~ Λ²` (condensate) — but this is *dimensional analysis*, not a master-field identification.
- **Borel structure differs.** `tr F²` (dim 4) ⇒ leading IR renormalon at `u=2` (witness: `a_n ~ Γ(n+1+2)`,
  ratio `a_n/a_{n-1} ~ n+2`). The dim-2 condensate ⇒ `u=1` (`a_n ~ Γ(n+1+1)`, ratio `~ n+1`). Different,
  structurally cleaner — the earlier claim (sidesteps the `u=2` ambiguity) is correct.
- **Consistency band.** Dual GL predicts `m_H/√σ = κ_GL/√(π ln κ_GL)` (witness: 1.36 / 2.22 / 3.72 at
  `κ_GL=2/5/10`). Large-N extrapolation gives `m_{0⁺⁺}/√σ → 3.30` (`1/N²` fit from `m(SU2)=3.79`,
  `m(SU8)=3.33`). Lattice `≈3.9`. **All lie in `3.3–3.9`.**

---

## §3 · The obstruction `[NG]/[O]`

- **No controlled Abelian projection at large N.** Dual superconductivity rests on 't Hooft's Abelian
  projection, which is *gauge-dependent*; monopole dominance is a lattice finding, not a theorem, and the
  Abelian-projection map is not known to commute with `N→∞` or to produce a master field for the monopole
  sector. There is no known "dual master field" reducing to `W(C)`.
- **`v` is not the master-field gap.** The master field's gap is `m_{0⁺⁺}(λ) ~ Λ`; `v ~ Λ` numerically, but
  the *identification* `v = ` (master-field order parameter) is unproven — a coincidence of scale, not an
  identity.

---

## §4 · The hypothetical combined object, and the honesty check

If a dual master field existed, its payoff would be to **fix `κ_GL`** (= `m_H/m_V`) from the large-N saddle,
converting the one-parameter dual-GL **fit** into a **prediction** of `m_{0⁺⁺}/√σ`. Witness: dual GL reaches
the lattice value `3.9` at `κ_GL ≈ 10.6`; large-N gives `3.30`. **Without the dual master field, `κ_GL` is
free** — so "`m_{0⁺⁺}/√σ≈3.3–3.9`" is a *consistency band*, not a derivation. This is the crux: the combined
object's whole content is the missing bridge, so its "prediction" cannot be counted until the bridge exists.

---

## §5 · Exact failing step `[O]`

A **large-N controlled Abelian projection** `Π`: (gauge master field `W(C)`) → (dual monopole sector), with
`v` shown to be an order parameter of `Π`. Neither `Π` nor the identification exists; Abelian projection is
gauge-dependent at finite `N` and uncontrolled at large `N`. This is the same wall as `H1`, in the
projection coordinate.

---

## §6 · Bottom line

- **No dual master field.** The two renormalon-free routes do not meet; the bridge breaks at large-N
  Abelian projection `[NG]/[O]`.
- **What is pinned:** `v ~ Λ` / `v² ~ Λ²` (dim transmutation, `N`-independent); dim-2 Borel pole `u=1` vs
  `u=2`; consistency band `m_{0⁺⁺}/√σ ≈ 3.3–3.9` is a **fit, not a prediction** (`κ_GL` free).
- **No QED; Clay 0%.** The combined object is a real target (it would fix `κ_GL`), but its existence is the
  obstruction, not a solution.

*Witness:* `w7/pnt_largen_dual_combo.py` (math only, `py -3.12 -I`).
*Refs:* `ROOF-LARGEN-MASTER.md`, `ROOF-DUALSUPERCONDUCTOR.md`, `INVENT-SYNTHESIS.md §2,§5`.
