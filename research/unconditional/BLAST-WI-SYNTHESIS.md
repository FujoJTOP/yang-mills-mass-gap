# BLAST-WI-SYNTHESIS — the final H1 blast fails 0/8 (every route reduces H1 to one marginal inequality / a→0 / IR endpoint); returning to the two WI lines yields genuine `[C]` diagnostics, but they are A3/I2 *reads* gated by `L→∞ = H1`

> **VERDICT (no QED, Clay 0%).** Two final waves. **Blast H1 (8 routes): 0/8 — H1 cannot be blasted.**
> Every route collapses H1 onto the *same* irreducible form in a different coordinate — one marginal `O(1)`
> inequality, the `a→0` lattice→physical map, the nonperturbative IR/Landau endpoint, a circular squeeze, or the
> renormalon σ-ambiguity — and **none touches the atom half** (`f₁²>0`). **WI continuation (8 routes): genuine
> `[C]` diagnostics** (P3 approach-law `T(τ)`, Re-slope `R` statistic, quantitative `Δ ≥ κ(τ)`, cross-validated
> `(E₁,Δ)`), acid-clean and continuum-robust *in the ideal*, **but** (i) a single finite `L` cannot decide
> (compactness ⇒ every discrete spectrum looks exponential), (ii) the discriminating content is the `L→∞`/`t→∞`
> scaling of `δ_L` = **isolation `I2` = `H1`**, (iii) on **finite real-τ data the Re-slope is blind** (`O(h)`)
> and collapses to P3. Both waves confirm `ONLY-EXIT`: the only exit is `H1` (the atom = the weak-*
> discontinuity of the boundary density).

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · Blast H1 (8 routes): 0/8

| route (file) | coordinate in which it *is* H1 |
|---|---|
| **Operator inequality** `BLAST-OPERATOR-INEQ` | `[T]` Sherman–Morrison `Var(a) ≤ 1/(μ(1−g|B|α))`, **sharply softened by the rank-one negative mode** (equality at `v∥φ̂`) ⇒ H1 ⟺ `μ(1−g|B|α) ≥ λ` = `D<S` |
| **YM nonlinear (Brascamp–Lieb)** `BLAST-YM-NONLINEAR` | `[T]` the genuine nonlinear field-level tool reduces the **χ half** to ONE marginal `O(1)` inequality `S_bare > 2g²/λ₄` (constant = condensate = H1); **blind to the atom** (`:φ²:` same χ, `f₁²=0`) |
| **Exclude the atomless bump** `BLAST-EXCLUDE-BUMP` | `[T]` no RP/OPE/Ward/nonlinear property separates the flat bump (matches `M₀,χ,A` to `O(h)`, `f₁²=0`); excluding it = isolation `I2` = H1 |
| **Certified lattice** `BLAST-LATTICE-CERT` | `[T]` `0 ≤ χ_lattice ≤ 4π²` L-uniform, but the wall is **`a→0`** (lattice→physical map uncontrolled by compactness), not `L→∞`; positivity interval `−0.035` excludes nothing |
| **Transfer/variational** `BLAST-TRANSFER` | `[T]` circular — sharpest variational upper bound vacuous; the only sharper one is squeeze-with-gap `χ≤C(0)/Δ²` (needs the gap); Dirichlet/Rayleigh/BL give the **wrong direction** (lower bounds on χ) |
| **Close Savvidy** `BLAST-SAVVIDY-CLOSE` | `[T]` pinned to the **IR/Landau endpoint**: one-octave `Σ≈0.35–1.40` vs multi-cell IR-marginal `29.9→264` (margin −80…−166); whether the softening is one-octave-bounded or full-flow marginal is **itself H1** |
| **String tension** `BLAST-STRING-TENSION` | `[C]` `H1 ⟸ (β<1/48 bridge) ∧ (λ ≤ 16σ/M)`; fails at 2.5; **proves only χ, never `f₁²>0`** (tight at the atomless bump) |
| **Resummation** `BLAST-RESUMMATION` | `[V-num]` diagnostic only: Borel–Padé recovers PV but the physical `χ = PV + σe^{−t_c/u}` has a **free σ** (20–45% ambiguity) = the condensate = H1 |

**0 / 8 blast. The χ half collapses to one marginal `O(1)` inequality (or a→0 / IR endpoint / σ); the atom half is untouched.**

---

## §2 · WI continuation (8 routes): genuine `[C]` diagnostics, all gated by `H1`

| route (file) | result |
|---|---|
| **P3 continuum** `WI-P3-CONTINUUM` | `[T]` the approach law **survives a continuum and is acid-clean** (Watson: atom below threshold ⇒ exponential rate `Δ`; pure continuum ⇒ power-law), **operator-choice-robust** — the robust isolator the complex-zero `𝒜` is not |
| **P3 rigor** `WI-P3-RIGOR` | `[T]` infinite-volume discriminator `T(τ):=−d/dτ log|m−ê₁|` (atom ⇒ `T→Δ`; cut ⇒ `T→1/τ→0`); **finite-L: every discrete spectrum looks exponential ⇒ a single L cannot decide** |
| **P3 quant** `WI-P3-QUANT` | `[T]` quantitative one-sided bound `Δ ≥ κ(τ)→Δ` (sharper than qualitative); **but read-not-force (needs the mass `e₁`), window-limited** (thin `Δ` under-certified) |
| **P3 uniform** `WI-P3-UNIFORM` | `[T]` FSS law `Δ(L)=m₁+c/(m₁L²)` (atom) vs `(2π/L)²/m→0` (cut): **uniform `Δ(L,a)≥δ` ⟺ `m₁>0` ⟺ the gap** |
| **Re-slope lattice** `WI-RESLOPE-LATTICE` | `[V-num]` with exact complex-τ zeros the Re-slope `R` statistic is robust (ATOM+GAP `R=5` where the raw nearest-zero spread 1.82 fails); **but exact complex-τ = A3**, more data-hungry |
| **Re-slope limit** `WI-RESLOPE-LIMIT` | `[T]` the Re-slope is **A3 by necessity** (reads the simple pole); finite real-τ is **blind** (`atom=bump` `O(h)`); **fragile** (atom+weak continuum breaks rung-invariance) |
| **Re-slope finite** `WI-RESLOPE-FINITE` | `[T]` no independent finite-real-data Re-slope — it **collapses to P3** (plateau + exponential-vs-power-law), reads the gap, needs `t→∞` = H1 |
| **Combine** `WI-RESLOPE-COMBINE` | `[C]` a **cross-validated `(E₁,Δ)`** estimator (Re-slope gives `E₁`, approach-law gives `Δ`), robust to continuum in the ideal — a real diagnostic, but `E₁`=A3, uniformity=`H1` |

**8 / 8: genuine tools (acid-clean, continuum-robust in the ideal), but A3/I2 reads, not certificates.**

---

## §3 · The unified conclusion (both waves confirm ONLY-EXIT)

- **The χ half** (susceptibility) reduces — via a genuine nonlinear tool (Brascamp–Lieb + the rank-one Savvidy split) — to **one marginal `O(1)` inequality** whose constant is the condensate `⟨tr F²⟩` = `H1`. No blast.
- **The atom half** (`f₁²>0` with isolation) cannot be reached by any finite-data tool: the atom is the **weak-* discontinuity** of the boundary density; every finite datum is weak-* continuous and atomless measures are weak-* dense. The WI tools *read* it only via analytic continuation (= A3) or infinite data (= H1).
- **A single lattice cannot decide**: compactness ⇒ every discrete spectrum looks exponential; the verdict is the `L→∞`/`a→0` scaling of the isolation `δ_L` = `I2` = `H1`.

⇒ **The only exit is `H1`** — the `a→0`, `L→∞` survival + isolation of the `0⁺⁺` atom. Both the blast and the WI continuation re-derive this from opposite directions.

---

## §4 · The genuine tools yielded (the campaign's best)

All `[C]`/`[V-num]`, acid-clean, lattice-testable — **diagnostics, not certificates**:
- **P3 approach-law** `T(τ)` (atom ⇒ exponential/`T→Δ`; cut ⇒ power-law/`T→1/τ`), robust to a continuum, operator-robust.
- **Quantitative isolation** `Δ ≥ κ(τ)→Δ` (one-sided, window-limited).
- **Re-slope `R` statistic** (robust with exact complex-τ zeros; A3 in practice).
- **Cross-validated `(E₁,Δ)`** estimator.
- **FSS law** `Δ(L)=m₁+c/(m₁L²)`.
- **Rigorous `m_{0⁺⁺} ≥ 4√σ`** (strong-coupling anchor); **`[C]` `H1 ⟸ (β<1/48) ∧ (λ≤16σ/M)`**.

---

## §5 · Honest bottom line

- **Blast H1: 0/8.** H1 cannot be blasted; it is the irreducible wall, reached from 8 coordinates.
- **WI continuation: genuine `[C]` diagnostics**, but all are A3/I2 *reads* gated by `L→∞ = H1`; finite real data is blind.
- **Both waves confirm ONLY-EXIT**: the only exit is `H1`.
- **Clay = 0%.** No QED. The campaign now has (1) a rigorous conditional reduction `H1 ⇒ gap` (machine-checked in Lean), (2) a proof `H1` is weaker than the gap (the atom is blind to H1 alone), (3) an exhaustive map that the only exit is `H1` itself, and (4) a set of the best available lattice diagnostics for the atom/isolation — but **no proof of H1 and no gap**.

*Refs:* all `BLAST-*.md` (8), all `WI-*.md` (~13), `ONLY-EXIT-SYNTHESIS.md`, `TOOL-SYNTHESIS.md`,
`SAVVIDY-WI-SYNTHESIS.md`, `lean/YMStatement/H1Sharpening.lean`.
