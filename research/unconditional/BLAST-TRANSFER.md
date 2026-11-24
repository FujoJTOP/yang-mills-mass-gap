# BLAST-TRANSFER — the transfer-matrix/variational route to H1: the sharpest variational upper bound on `Var(a_{k→0})` is the operator-range bound `(range/2)²` (vacuous); tightening it to `1/λ` needs the vacuum spectral gap — circular (= H1)

> **VERDICT (no QED).** Direct attack on H1 (`χ = Var_ν(a_{k→0}) ≤ 1/λ`) via the transfer matrix: build a
> **variational upper bound** on the coarse-mode variance (maximise `Var_ψ(a)` over trial states). The
> sharpest such bound is **exact and available** — `Var_Ω(a) ≤ ((a_max − a_min)/2)²`, the operator-range
> bound, attained by the equal superposition of the max/min eigenvectors — **but for the coarse mode it is the
> field-range bound, vastly larger than `1/λ` ⇒ vacuous `[T]/[V-num]`.** Tightening it to `1/λ` requires the
> vacuum spectral gap: `χ = Σ|c_n|²/m_n² ≤ C(0)/Δ²` where `Δ` is the lightest coupled mass. So the only
> non-vacuous variational bound is the **squeeze-with-gap**, and it needs `Δ² ≥ λ C(0)` — a **gap lower bound**
> `[T]`. At finite `L` the compact transfer matrix supplies the **box gap** `4 sin²(π/L) → 0`, giving H1
> **vacuously**; the `a→0, L→∞` uniform version is the **gap itself** `[O]`. **The blast is circular in the
> transfer-matrix coordinate** — it reduces H1 to the gap. Consistent with the weak-* `ONLY-EXIT` picture.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.
Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · Setup (transfer matrix / Hamiltonian)

On the cylinder `T³_L×ℝ`, let the coarse mode be the zero-momentum smeared `a_{k→0}` (a gauge-invariant
functional of `tr F²` at `k=0`), `ξ` its centred version. With transfer matrix `T=e^{-aH}`, `H≥0`, vacuum `Ω`
and masses `m_n = E_n − E_0`, the KL measure is `ν = Σ_n |c_n|² δ(s−m_n²)`, `c_n=⟨n|ξ|Ω⟩`, and

```
χ  =  Var_ν(a_{k→0})  =  ∫ s⁻¹ dν  =  Σ_n |c_n|² / m_n²
   =  ⟨ξ | (H−E₀)⁻² P_⊥ | ξ⟩  =  ∫₀^∞ t · C(t) dt ,   C(t)=Σ_n|c_n|²e^{−m_n t}.
```

**H1:** `χ ≤ 1/λ`. (The campaign's target: this is *necessary* for the gap, the squeeze below.)

---

## §2 · The variational upper bound (the directive): operator-range bound `[T]` — and why it is vacuous `[V-num]`

Maximising the variance of `a` over **all** states `|ψ⟩` gives the **sharp variational upper bound**

```
Var_Ω(a) ≤ max_{|ψ|=1} Var_ψ(a) = ((a_max − a_min)/2)² ,
```

attained at `ψ = (|a_max⟩ + |a_min⟩)/√2`. This is a genuine `[T]` upper bound on the coarse-mode variance
(`Ω` is one state). **But for `a_{k→0}` it is the operator-range bound** `~ (range of tr F² on the torus)²`
— an `O(Λ⁴ · V²)`-type object, **vastly larger than `1/λ`** ⇒ **vacuous for H1**. Any "maximise the test
function" is this bound and no better (the max-variance state is the two-extremal-eigenvector superposition,
witness Part 1). So the *unconstrained* variational upper bound exists but is `≫ 1/λ`.

---

## §3 · Tightening to `1/λ` needs the vacuum gap — the squeeze `[T]`

The vacuum variance is `χ = Σ_n |c_n|²/m_n²`. Writing `C(0)=Σ|c_n|² = Var(ξ)` (equal-time) and `Δ² = min_n m_n²`
over the states `ξ` couples to,

```
χ = Σ_n |c_n|²/m_n²  ≤  (Σ_n |c_n|²)/Δ²  =  C(0)/Δ²      (the squeeze's upper half).
```

⇒ `χ ≤ 1/λ` **iff** `C(0)/Δ² ≤ 1/λ`, i.e. **`Δ² ≥ λ·C(0)`** — a **lower bound on the lightest coupled mass
relative to the equal-time variance**. That is a **gap lower bound**, not an H1-independent datum. **This is
the `[T]` squeeze**, and it is the *only* variational bound sharper than the range bound — because it uses the
vacuum spectral decomposition (the gap).

---

## §4 · Direction obstruction: variational/Dirichlet/BL give *lower* bounds on variances `[T]`

The natural variational tools (Rayleigh–Ritz, the Dirichlet principle for the resolvent, Brascamp–Lieb) yield
compute `χ = ||(H−E₀)⁻¹ξ||²` from **below** (restricting the trial space/subspace ⇒ `χ_lower ↑ χ`) and bound
energies from above — the **opposite direction** to H1's *upper* bound on `χ`. Witness Part 3: subspace
Dirichlet bounds ascend to `χ` from below, never above. Thus the transfer-matrix machinery the campaign has
used (Temple/Weinstein on `E₁`, BL on the fiber) is aimed at the **wrong inequality** for H1: those supply
**lower** bounds on variances / **upper** bounds on energies, whereas H1 is an **upper** bound on a variance.

---

## §5 · Finite `L`: the compact transfer matrix gives H1 *vacuously* via the box gap `[V-num]`

On a finite `L` lattice the transfer matrix is compact ⇒ discrete spectrum, finite gap. But the only gap
available without input is the **kinematic box gap** `Δ²_box = 4 sin²(π/L) → 0` (free-scalar check, witness
Part 2: zero-mode `Var = 1/(L·m²)` for mass `m`, `∞` for massless; the massless coarse mode's `χ ~ L² → ∞`).
So at finite `L`, `χ ≤ C(0)/Δ²_box` is finite — H1 holds **vacuously**, with the box gap, and **fails uniformity
as `L→∞`**. This is exactly the `FINITE-VOL-CERTIFY` / `CONSTRUCT-SAVVIDY-CELL` observation (`G″=+0.28` at `L=4`
→ `+0.04` at `L=8`: the small-volume box signature, not an IR floor).

---

## §6 · 3+1D, `L→∞`: exact missing input `[O]`

The transfer-matrix route reduces H1 to:

```
a uniform-in-(a,L)  lower bound  Δ²_phys ≥ λ·C(0)   (the 0⁺⁺ lightest mass vs the equal-time variance).
```

That *is* the gap (`m₁² ≥ λ C(0)` is the `PNT-CHI-CERTIFICATE` Lemma-G sufficiency `m₁² ≥ λM`). So the blast is
**circular in the transfer-matrix coordinate**: it needs the gap to prove H1.

---

## §7 · Verdict

- `[T]` A genuine variational upper bound on `Var(a_{k→0})` exists — the operator-range bound `((a_max−a_min)/2)²`
  — but it is **vacuous** (`≫ 1/λ`).
- `[T]` The only bound sharper than it is the **squeeze with the gap** (`χ ≤ C(0)/Δ²`), which requires
  `Δ² ≥ λ C(0)` — a gap lower bound.
- `[T]` Variational/Dirichlet/BL principles run the **wrong direction** for H1 (lower, not upper, bounds on
  variances).
- `[V-num]` Finite `L`: H1 holds vacuously via the box gap `4 sin²(π/L) → 0`.
- `[O]` `3+1D, L→∞`: the missing input is the uniform gap `Δ²_phys ≥ λ C(0)` — the gap = H1.
- `[NG]` No transfer-matrix variational upper bound on `χ` below the squeeze-with-gap exists in this route.

**BLAST-TRANSFER: circular; reduces H1 to the gap.** No new handle; consistent with `ONLY-EXIT-SYNTHESIS`
(the atom/gap is the weak-* discontinuity; the transfer-matrix route supplies only the vacuum continuum data
plus the box gap). **Clay 0%; no QED.**

*Refs:* `PNT-TEMPLE-LOWER.md`, `PNT-CHI-CERTIFICATE.md`, `FINITE-VOL-CERTIFY.md`, `CONSTRUCT-SAVVIDY-CELL.md`,
`ONLY-EXIT-SYNTHESIS.md`, `PROFILE.md`; witness `w7/blast_transfer.py`.
