# G2C-MINWEIGHT-RIGIDITY — restoring the atom to the susceptibility inequality: the reduction is exact, the floor *is* the gap

> 2026-11-23 · Test of the atom-rigidity / minimum-weight proposal: prove every isolated low-lying `0⁺⁺` atom has `f₁²≥r₀>0`, then H1 (`χ≤1/λ`) plus min-weight forces a gap. **判词 (verdict):** the CONDITIONAL is an exact, one-line theorem — and it is the cleanest statement of why the bounded-`χ` counterexample is the *only* obstruction. But the ANTECEDENT `f₁²≥r₀` is *not* a consequence of reflection positivity, positivity of norms, or the `f⊗f̄→1` recoupling. Those supply only `f₁²≥0` (a *sign*), which the two-mass generalized free field — an admissible OS theory with an arbitrarily small lightest residue — saturates. Min-weight is therefore a *dynamical* (confining) statement, verbatim P3's I1, proven only at `β<1/48`, fixed `a`; uniform in `a` it is the open Clay core. **Score 1/5; not achieved.**

## 1 · The reduction is exact — and it is literally "put the atom back"

Let `ν` be the KL measure of the `0⁺⁺` interpolator (`ν≥0`, `supp ν⊂[m₁²,∞)`, `f₁²=ν({m₁²})`). The susceptibility is the `M₋₁` moment, `χ=∫s⁻¹dν=M₋₁`. The lightest atom contributes its weight over its mass; all other mass contributes non-negatively, so

$$\chi=\int s^{-1}d\nu\;\ge\;\frac{\nu(\{m_1^2\})}{m_1^2}=\frac{f_1^2}{m_1^2}\quad\Longrightarrow\quad \boxed{\,m_1^2\;\ge\;\frac{f_1^2}{\chi}\,.}$$ 

Grant H1 (`χ≤1/λ`) **and** min-weight (`f₁²≥r₀`): then `m₁²≥r₀/χ≥r₀λ>0` — a gap `m₁≥√(r₀λ)`, uniform whenever `r₀,λ` are. On the counterexample `ν=wδ_ε+(1−w)δ_S` (`w/ε≤1/λ`): min-weight forces `w≥r₀`, hence `χ≥r₀/ε→∞`, contradicting `χ≤1/λ`; so `ε≥r₀λ>0`. The logic is impeccable. It needs **only** the `M₋₁` moment plus the atom — **no** second-level threshold `Λ`, unlike the residue–level-gap lemma (`GAP-G2-edge`) or the two-scale-Poincaré chain of `RESIDUE-CONDITIONAL`. That is the genuine (small) nugget: H1 + atom-floor ⟹ gap with no isolation input, because `χ` absorbs the atom and a heavy tail only *lowers* the admissible `χ`.

## 2 · Why `χ` alone cannot do it — the moment is blind to the edge

This is *conditional*, not a proof, by the **direction obstruction** (`BREAK-IN-11`, `GAP-G2-edge`): `χ=M₋₁` is a *total*, and the moment-ratio lemma gives `m₁²≤M₀/χ` — an **upper** bound on the edge, never lower. Turning `χ` into a *lower* bound on `m₁²` requires inserting the *atom* `f₁²` (not a moment) into the numerator, exactly as in the boxed identity. `f₁²` is the one spectral datum positivity cannot resolve: the route correctly isolates the only missing object, but does not supply it.

## 3 · Positivity, norms, recoupling give the sign, not the floor

- **Positivity of norms.** `f₁²=|⟨0|O|G⟩|²` is a modulus-square; `≥0` is all it certifies, saturated at `0` by `O=:φ²:` in the free scalar. The distance from `≥0` to `≥r₀` is dynamics.
- **`f⊗f̄→1` recoupling.** Its coefficient `1>0` (`PROOF-P3 §3`) makes the leading strong-coupling *polymer* positive (`a₀>0`), fixing a **sign/normalization**, not a weight floor. The residue `f₁²=a₀(1−O(u))` carries `a₀=c₁(β)u^p`, `β`-dependent and `→0` as `u→0`; the floor `r₀` came from the **KP/OS confinement control**, not the recoupling.
- **No quantization.** The residue is a continuous modulus-square — `Z_O²` running in YM, arbitrary `w∈(0,1)` in a GFF. Nothing forces a discrete multiple; no rigidity/quantization principle exists.

**Decisive acid-level counterexample.** The two-mass generalized free field, `ν=wδ_{m₁²}+(1−w)δ_{m₂²}`, `m₁<m₂`, is reflection-positive and OS-admissible (Källén–Lehmann; Streater–Wightman). Its lightest atom is *isolated* below any `Δ>m₁²`, below the gap to `m₂²`, and has weight `w` **arbitrarily small**. So min-weight is **false as a structural theorem** — it is not implied by RP/positivity/recoupling; it is a confinement statement (area law ⇒ single flux tube is the unique leading weight, `PROOF-P3`).

## 4 · What min-weight *is*: I1, not a new theorem

Min-weight at strong coupling **is** `PROOF-P3`'s I1: `f₁²≥r₀(β)=½c₁u^p>0`, `V,L`-uniform, `β<1/48`, fixed `a`. It fails to transport across the Savvidy crossover (`INVENT-TRANSPORT-INEQ` needs `(CP)`); uniform-in-`a` min-weight is `H*(r₀)` (`RESIDUE-CONDITIONAL`), which is the gap. Nothing in §1–§3 crosses that wall.

## 5 · Five-point scorecard (target: uniform min-weight / atom rigidity)

1. **Nonperturbative — PARTIAL.** The reduction `m₁²≥f₁²/χ` is exact, nonperturbative; the antecedent is known only at `β<1/48`.
2. **Sign-controlling `m_eff²(k→0)≥λ` — HOLLOW.** Never touches `m_eff²`; demands a *lower* bound on an atom — the precise direction `ν≥0` cannot supply.
3. **Non-circular — FAIL.** Uniform min-weight `≡` I1 `≡` `H*(r₀)` `=` the concentration face of the gap.
4. **Non-SUSY + fails free field — FAIL as a structural tool.** The two-mass GFF is an admissible OS theory with an arbitrarily small lightest atom: min-weight is falsified outside confinement, i.e. it is dynamical, not axiomatic.
5. **Uniform `≥λ` beyond `1/48` — FAIL.** Same analyticity/transport wall (`a→0`, `β→2.5`).

**Total: 1/5. `5/5` NOT achieved.**

## 6 · Exact missing step

`f₁²≥r₀` with `r₀` **uniform in `(V,a)`** — equivalently the atom `` = the unprotected Schwinger order parameter. Given it, `m₁≥√(r₀λ)` follows in one line from `χ≤1/λ` (and *only* from it — no isolation). Positivity/recoupling fix its sign, not its size; the GFF shows no axiom can.

**Sharpest line.** The route is correct and clarifying: bounded `χ` coexists with a light edge *iff* the edge weight is tiny, so H1 + min-weight is genuinely sufficient and `m₁²≥f₁²/χ` is its sharpest single packaging. But it merely *names* the obstruction — an atomic floor that positivity provably cannot give and the GFF sets to zero. The minimum weight is not rigid; it is the dynamics.

*Refs.* `RESIDUE-CONDITIONAL` (H*(r₀)); `GAP-G2-edge` (residue–level-gap lemma; moment-ratio = upper; "RP theories with arbitrarily small gaps exist"); `BREAK-IN-11` (direction obstruction); `PROOF-P3-UNCONDITIONAL` (I1, `f⊗f̄→1` recoupling); `PROOF-ISOLATION` (I2 open); `PROOF-F1SQ` (free `:φ²:` has `f₁²=0`); `INVENT-TRANSPORT-INEQ` (`(CP)`); `D3-SUSCEPTIBILITY-VARIATIONAL` (`χ=M₋₁`, `dρ≥0`); `REDUCTION-AUDIT`. *External.* Osterwalder–Schrader CMP **31** (1973) 83, **42** (1975) 281; Källén (1952), Lehmann (1954); Streater–Wightman (1964, GFF is OS); Osterwalder–Seiler Ann. Phys. **110** (1978) 440; Shen–Zhu–Zhu CMP **400** (2023) 805, arXiv:2204.12737.
