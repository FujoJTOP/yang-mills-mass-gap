# WI-CONTINUUM-SUBTRACT — no constructible local continuum-cancelling subtraction: densities are quadratic (a parallel pair cancels the atom too), and an exact cancellation needs a continuum-null local `0⁺⁺` operator that the `:φ²:` KL cut forbids

> **VERDICT (no QED).** The continuum-cancelling subtraction is **definable but not constructible** for
> local `0⁺⁺` operators. Three exact obstructions. **(1) `[T]` densities are quadratic in the operator**:
> a linear reweighting `O−cO′` gives `ρ_{O−cO′}=|g_O−c g_O′|²` (cross terms); it cancels the continuum
> only if `g_O∥g_O′` (Cauchy–Schwarz), which cancels the **atom too** — so a **generic pair cannot cancel
> the cut while keeping the atom**. **(2) `[T]` an exact cancellation** (a designed pair with identical
> continuum vector, different atom) **works — but requires a *continuum-null* operator** `O−O′`
> (`⟨0|O−O′|n⟩=0` for all continuum `n`). **(3) `[T]` such a local operator does not exist**: a local
> `0⁺⁺` operator has a Källén–Lehmann **cut** from multi-glueball states (exactly the `:φ²:` acid test);
> a pure-single-exponential `2`-point function is the acid-violating free-composite object. The only
> continuum-null object is the **non-local glueball eigenvector** — constructing it **is the gap = `H1`**.
> **Clay 0%.**

---

## §0 · Setup and labels

Spectral function of a `0⁺⁺` operator `O`:
```
Z_O(τ) = w₀ + f₁² e^{−E₁ τ} + ∫_θ^∞ ρ_O(E) e^{−E τ} dE,     ρ_O(E)=Σ_n |⟨0|O|n⟩|² δ(E−E_n),
```
atom at `E₁=m₁²`, weight `f₁²`; continuum threshold `θ>E₁`.
Labels: **[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go.

---

## §1 · Densities are quadratic; a generic pair cannot cancel the cut and keep the atom `[T]`

For a linear reweighting `O_c = O − c O′` the continuum matrix elements are `⟨0|O_c|n⟩ = ⟨0|O|n⟩ − c⟨0|O′|n⟩`,
so the continuum **density** is
```
ρ_{O_c}(E) = ρ_O(E) + c² ρ_{O′}(E) − 2c Re ω_{O,O′}(E),
```
with `ω` the cross-density. This is **quadratic** in `c` — the operator does **not** subtract densities
linearly. Minimising `∫ρ_{O_c}` over `c` leaves residual `> 0` unless the matrix-element vectors
`g_O=(⟨0|O|n⟩)`, `g_{O′}` are **parallel** (Cauchy–Schwarz equality), i.e. `⟨0|O′|n⟩ = λ_n⟨0|O|n⟩`
with `|λ_n|=const` on the continuum.

**Witness (A)** — generic (non-parallel) pair: min continuum residual / total = **0.24**, and at the
minimising `c*` the **atom amplitude is 0.57 ≠ 0**. You cannot simultaneously cancel the cut and keep the
atom: the continuum vector that must be annihilated by `g_O − c g_O′` is not the atom.

---

## §2 · An exact cancellation exists — iff a continuum-null operator exists `[T]`

Take a **designed** pair with the **same** continuum vector (`g_{O′}=g_O`) but a **different** atom
amplitude (`a_{O′}≠a_O`). Then `O − O′` has continuum amplitude `g_O − g_O = 0` and atom amplitude
`a_O − a_{O′} ≠ 0`:

**Witness (B)** — continuum residual at `c=1`: **0.0000e+00**; atom amplitude `−2.0` (pure atom recovered).

So the subtraction device **exists mathematically**, conditioned on an operator pair whose continuum
matrix-element vectors are **proportional**. Equivalently, it needs a **continuum-null local `0⁺⁺`
operator** `O_c` with `⟨0|O_c|n⟩=0` for every continuum state `n` and `⟨0|O_c|atom⟩≠0`.

---

## §3 · The obstruction: no local continuum-null `0⁺⁺` operator `[T]`

A local `0⁺⁺` operator's `2`-point function is a Källén–Lehmann integral with spectral density
`ρ_O(E)=Σ_n|⟨0|O|n⟩|²δ(E−E_n)`. In **interacting** Yang–Mills the multi-glueball states contribute
`⟨0|O|GG…⟩≠0` (the operator is local ⇒ it couples to all states of its quantum numbers), so `ρ_O` has a
**branch cut** above the two-glueball threshold. A density that is a **single δ** (a pure single
exponential `2`-point function) is precisely the object that **fails the discriminating `:φ²:` acid test**
— a composite whose cut is absent behaves like a free Wick square. Hence:

> **A local `0⁺⁺` operator with `⟨0|O_c|n⟩=0` on the whole continuum does not exist in interacting YM.**
> The only continuum-null object is the **non-local glueball creation operator** (the exact eigenvector of
> the transfer matrix in the `0⁺⁺` channel). Writing it down **is** the `0⁺⁺` eigenvalue problem = `H1`.

This is the structural reason behind the WI-ADMISSIBLE downgrade: the quantity one would want to subtract
(the continuum) is glued to the operator by locality; separating it requires the non-local eigenvector.

---

## §4 · `𝒜(O) − 𝒜(O′)` does not cancel `[NG]`

`𝒜 = π/|Im τ_zero|` is the **zero location** of `Z_O` — a **nonlinear** functional of the whole spectral
function (the cut enters through the branch structure, not additively). Therefore `𝒜(O) − 𝒜(O′)` is the
difference of two independently **drifted** values; it cancels nothing. The nearest-zero contamination is
also threshold/weight dependent (witness C: a **spread** low-density continuum shifts `𝒜` by only
`~10⁻³`; a **concentrated** extra level creates a new zero and a large jump, per WI-ADMISSIBLE). So no
generic `𝒜`-difference is a continuum-canceller.

---

## §5 · Verdict

- **Definable:** given a pair with proportional continuum matrix elements, `O − cO′` cancels the cut exactly
  and leaves a pure atom (witness B).
- **Not constructible:** such a pair needs a **continuum-null local `0⁺⁺` operator**, which the
  Källén–Lehmann cut (`:φ²:` acid test) forbids. The only continuum-null object is the **non-local
  glueball eigenvector = `H1`**.
- **No generic device:** the zero-location `𝒜` is nonlinear; `𝒜(O)−𝒜(O′)` cancels nothing.
- **So** the continuum-cancelling subtraction is a correct idealisation but is **exactly as hard as the
  gap** — it requires knowing the continuum eigenstates. It does not reduce `H1`; it re-expresses it in
  the "continuum-null operator" coordinate. **Clay 0%; no QED.**

*Refs:* `WI-ADMISSIBLE.md` (the downgrade), `WI-THEOREM.md`, `WI-GENERAL.md`, `CONSTRUCT-IR-ANALYTIC-NUM.md`;
witness `w7/wi_continuum_subtract.py`.
