# CERT-CHI-LATTICE — the `a→0` wall is exactly `χ_lattice ≲ a²` = dimensional transmutation = `H1`; no renormalization-condition / step-scaling / monotone argument supplies an `a→0` certificate

> **VERDICT (no QED).** BLAST-LATTICE-CERT left the wall at `a→0` (the lattice→physical map `Z(a)` is
> uncontrolled by compactness). I tried the three named routes for a rigorous `a→0` certificate of
> `χ_phys ≤ 1/λ`. **None works, and the reason is a clean identity:** compactness gives a **flat** ceiling
> `χ_lattice ≤ (2π)²`, while the physical demand `χ_phys ≤ 1/λ` is `χ_lattice ≤ a²/λ` — an **`O(a²)`**
> ceiling. The gap between the two is the `a²` coefficient, which is `χ_phys` itself (the inverse mass scale).
> **`χ_lattice ≲ a² ⟺ χ_phys < ∞ ⟺ the gap = `H1`.** The `a→0` extrapolation is **not a separate tool — it
> IS dimensional transmutation.** **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** 数值. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The exact scaling that makes the wall

Two ceilings on the **dimensionless** lattice susceptibility `χ_lattice`:

- **Compactness `[T]`** (BLAST-LATTICE-CERT): `0 ≤ χ_lattice ≤ (2π)² ≈ 39.5`, **a-independent**. It does **not** shrink as `a→0` — it is a flat cap on the group manifold diameter.
- **Physical finiteness `[T]`** (dimensional analysis): with lightest mass `m`, `χ_lattice(a) = c·(am)²` (`c` an O(1) number) ⇒ `χ_lattice → 0` like `a²` as `a→0`. Free massive field exactly: `χ_phys = 1/m² ⇒ χ_lattice = a²/m²`.
- **The demand `[T]`** (ROOF-UNIFORM-CHI): `χ_phys ≤ 1/λ ⟺ χ_lattice ≤ a²/λ`, an **`O(a²)`** ceiling.

⇒ The gap between "compactness `O(1)`" and "demand `O(a²)`" **is the `a²` coefficient**, i.e. `χ_phys` = the inverse mass scale = **the gap**. Witness §3: `demand/compactness = (a²/λ)/4π²` shrinks like `a²` (6.3e-3 → 1.0e-5 over `a=0.5→0.02`), while compactness stays flat at 39.5.

**So `χ_lattice ≲ a²` is not a technical extrapolation step — it is the statement that a mass scale exists.**

---

## §2 · The three named routes all fail

- **(i) Renormalization condition `[NG]`.** A renorm condition fixes a **UV** coupling at a scale `μ`; `χ_phys` depends on the **IR** mass scale, which no UV condition supplies. The condition sets the `a`-dependence of the coupling (asymptotic freedom — perturbative, UV), not the IR constant.
- **(ii) Step-scaling `[NG]`.** The RG step gives the **anomalous dimension** of `tr F²` (perturbative, UV). It fixes the **exponent** `p` in `χ_lattice ~ a^p` in the controlled UV regime, but the **constant** (the IR mass scale) is the nonperturbative input. Same as Balaban: UV control, not the IR gap.
- **(iii) Monotone quantity in `a` `[NG]`.** The only monotone candidates — the effective curvature `κ(a)` or `R_vac(a)` — **are** the existing H1 attempts (`TRANSPORT-SCALAR`: `R_vac` monotone in β is circular; `UNIFORM-CHI`: `χ_lattice ≲ a²` = dimensional transmutation). A monotone-in-`a` quantity that lower-bounds the `a²` coefficient reads the gap; it does not force it.

---

## §3 · The decisive identification

- At **strong coupling** (`β<1/48`, the rigorously controlled region) `χ_lattice = O(1)`, **a-independent** — there is **no small-`a` regime** in the proven region.
- Any `a→0` extrapolation **must cross the crossover**, where the `O(a²)` coefficient (the IR mass scale) is exactly the nonperturbative input.
- **Exponent underdetermination `[V-num]`**: `χ_lattice ~ a^p` with `p=2` (mass gap) and `p=0` (strong coupling) are **both** compatible with compactness `0 ≤ χ_lattice ≤ 4π²`. The exponent `p` **is** the gap; compactness cannot fix it. (Witness §6.)

⇒ **The `a→0` map is `H1` in the `a`-scaling coordinate.** No `a→0` certificate exists without the crossover control — which is `H1`.

---

## §4 · Honest bottom line

- **No `a→0` certificate.** The three named routes (renorm-cond, step-scaling, monotone) each fail; each needs the nonperturbative IR mass scale.
- **The wall is identified exactly:** `χ_lattice ≲ a²` = `χ_phys < ∞` = dimensional transmutation = `H1`. Compactness is `O(1)` and `a`-flat; the demand is `O(a²)`; the gap is the gap.
- **Consistent with ONLY-EXIT and BLAST-LATTICE-CERT:** the `a→0` continuum step is not a separate obstacle to be "certified"; it is the mass scale itself.
- **Clay 0%; no QED.**

*Witness:* `w7/cert_chi_lattice.py` (runs clean, `py -3.12 -I`, numpy). *Refs:* `BLAST-LATTICE-CERT.md`,
`ROOF-UNIFORM-CHI.md`, `TRANSPORT-SCALAR.md`, `ONLY-EXIT-SYNTHESIS.md`.
