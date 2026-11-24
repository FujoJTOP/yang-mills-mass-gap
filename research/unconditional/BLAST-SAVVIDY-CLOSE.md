# BLAST-SAVVIDY-CLOSE — the Savvidy H1 attack assembled: bounded Var(φ²) + (L,a)-uniformity + D-split are real, but the multi-cell accumulation does NOT close (overshoots toward the Landau pole); exact step = the IR endpoint

> **VERDICT (no QED, Clay 0%).** This is the sharpest concrete attack on `H1` assembled from the Savvidy
> results. The **positive pieces are real theorems**: the exact defect split `D = D₊(≤0) +
> E[|⟨φ₋|c⟩|²]/(g|B|)`, the bounded softening `Var(φ²)→2/λ₄` (c-independent), sign-protection, and
> `(L,a)`-uniformity by one-octave support. **They do NOT close `κ_bare > Σ` uniformly.** Under the multi-cell
> RG accumulation the softening is **IR-marginal** and **overshoots** for `O(1) κ_bare`: sub-criticality holds
> only for `k_IR ≳ 90Λ`, and fails toward the Landau pole. The **exact step** is the IR endpoint `k_IR ~ Λ`,
> where the perturbative `g²(k)=C_RG/ln(k/Λ)` diverges and the true sub-criticality is a **nonperturbative**
> (confinement-scale) statement = `H1`. The one-octave vs multi-cell tension is itself the `H1` content.

---

## §0 · Labels
**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** 数值.

## §1 · The assembled pieces (all real)

- **Exact split `[T]`** (`SAVVIDY-DTERM`): the Savvidy Hessian `= H₊ − P₋/|λ₋|` (rank-one negative mode,
  `λ₋=−g|B|`) gives `D = D₊ + E[|⟨φ₋|c⟩|²]/(g|B|)` with `D₊ ≤ 0` (Brascamp–Lieb on the positive fiber). So
  `D ≤ E[|⟨φ₋|c⟩|²]/(g|B|)` and `D<S ⟺` the rank-one term `< S_bare`.
- **Bounded softening `[T]/[V-num]`** (`SAVVIDY-VAR-PHI2`, `NEGMODE-EXACT`): `Var(φ²)=2/λ₄+2/(gB)²+…→2/λ₄`,
  **c-independent** (witness: `c=10,20,40,80,160 → 2.020, 2.005, 2.001, 2.0003, 2.0001`). The softening
  `g²Var(φ²)` is **O(1) and g-independent** at fixed `a=c/√λ₄` (witness: `g=1.0…0.3 → 4.000` constant).
- **`(L,a)`-uniform `[T]`** (`SAVVIDY-UNIFORM-SOFTENING`): the negative mode is LLL at `|p|~√(gB)`, so the
  softening is supported in **one octave**, not a log over `[1/L,1/a]` — no volume/cutoff divergence.

## §2 · The two accumulation pictures (they disagree)

**(A) One-octave** (`Σ_oct = g²Var·ln2`): bounded, `Σ_oct = 0.35, 0.66, 1.40` for `g²Var = 0.5, 0.95, 2.02`.
This *would* close for small `κ_bare` — **but** it assumes the negative mode contributes in exactly one octave
and that `κ_bare ≳ 1.4`, neither established.

**(B) Multi-cell RG** (`Σ = (2C_RG/λ₄)·ln( ln(k_UV/Λ)/ln(k_IR/Λ) )`, `C_RG=24π²/11=21.53`): **IR-marginal and
divergent** toward the Landau pole — `Σ = 29.9, 45.3, 61.7, 81.6, 104.6, 139.1, 167.0, 195.8, 264.3` for
`k_IR/Λ = 10, 5, 3, 2, 1.5, 1.2, 1.1, 1.05, 1.01`. For `κ_bare ~ 1–5` the margin `κ_bare−Σ` is
**−80 … −166** (massive overshoot). Sub-criticality `κ_bare>Σ` holds **only for `k_IR ≳ 90Λ`** (crossover
`k_IR* = 89.97Λ` at `κ_bare=1`).

**The two pictures differ by ~2 orders** because `λ₄(k)`'s scaling is not pinned. Whether the softening is
one-octave-bounded (A) or full-flow IR-marginal (B) **is itself the H1 content**: the RG accumulation is
controlled only if the negative-mode softening is genuinely localized to one octave at every scale, and
proving that localization is a nonperturbative statement about the effective action.

## §3 · The exact step (pinned)

`[O]` **The IR endpoint `k_IR ~ Λ`.** The multi-cell `Σ` diverges there because the one-loop AF coupling
`g²(k)=C_RG/ln(k/Λ)→∞` at the Landau pole; perturbative cell-by-cell accumulation is invalid in the last few
cells, which carry almost all of `Σ`. A uniform `κ_bare > Σ` therefore requires a **nonperturbative bound on
the accumulated softening at `k~Λ`** — i.e. at the confinement scale — which is exactly `H1` (the crossover),
not a perturbative RG statement. The physical cutoff (the negative mode dies below the gluon mass `~Λ`) gives
only a **marginal O(1)** with an uncontrolled constant, not a closed inequality.

## §4 · Honest bottom line
- **The Savvidy route is assembled and its positive pieces are real `[T]`** (exact D-split, bounded `Var(φ²)→2/λ₄`,
  sign-protection, `(L,a)`-uniformity).
- **It does NOT close `H1`.** Under multi-cell accumulation the softening overshoots for `O(1) κ_bare`;
  sub-criticality holds only in the UV (`k_IR ≳ 90Λ`), failing toward the Landau pole.
- **Exact step `[O]`:** the IR endpoint `k_IR~Λ`, where the perturbative softening diverges and the true
  sub-criticality is nonperturbative (confinement scale) = `H1`.
- **Clay = 0%.** The blast failed, and it pinned the failure to the single nonperturbative IR bound it cannot supply.

*Refs:* `SAVVIDY-DTERM.md`, `SAVVIDY-VAR-PHI2.md`, `SAVVIDY-NEGMODE-EXACT.md`, `SAVVIDY-UNIFORM-SOFTENING.md`,
`SAVVIDY-MULTICELL-RG.md`, `w7/blast_savvidy_close.py`.
