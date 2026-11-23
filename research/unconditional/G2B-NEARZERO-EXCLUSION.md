# G2B-NEARZERO-EXCLUSION — can any structure forbid a tiny-weight near-zero atom?

> 2026-11-23 · The exact missing sufficient step. Target: prove `ν` has no weight in `(0,Δ)`, or that any sub-`Δ` weight is bounded *below*, to contradict `χ≤1/λ`. **判词 (verdict): NO-GO, 1/5. The exclusion is FALSE as a universal statement, and the one YM structure the counterexample lacks — the trace anomaly — is provably BLIND to it. The susceptibility route is not merely insufficient; for the anomaly operator it is *structurally incapable* of seeing the light state.** No fabricated `5/5`.

## 1 · The counterexample every structure must face — and it is legitimate

`ν = wδ_ε + (1−w)δ_S`, `w/ε ≤ 1/λ` ⇒ `χ = ∫s^{-1}dν ≤ 1/λ`, yet `inf supp ν = ε → 0`. This is not pathological: a **generalized (Gaussian) free field** with spectral measure `ν` is OS/Wightman-positive **iff `ν ≥ 0`** (Jaffe). So `ν` *is* the KL measure of a genuine OS theory — bounded susceptibility, gap `0`. Any exclusion must kill **this** object. It satisfies RP and KL exactly; the only listed datum it does not carry is the trace anomaly (§3).

## 2 · The four named mechanisms — each dies

- **Reflection positivity** ⇒ `ν ≥ 0` only. Every positivity functional is Schur-concave, hence an **upper** bound on totals (GAP-G2-edge; INEQ-T3). The counterexample *is* RP. Cannot lower-bound an atom.
- **Clustering** ⇒ the cluster rate is `√(inf supp ν) = √ε`, *i.e. the gap itself*. Circular: assuming clustering at `m>0` assumes the answer.
- **Exact residue / `:`-normal-ordering** ⇒ the free `:φ²:` template gives `f₁²(m²)=0`, but that is an **isolated zero at a fixed pole** (an *equality*); `(0,Δ)` is an **open band**. Analyticity lets `f₁²` cross `(0,Δ)` without touching a discrete zero — a normal-ordering identity pins one point, never forbids a band. And in YM no symmetry forces even the zero (PROOF-TRANSPORT-ANALYTIC §4; PROOF-DISSOLUTION-SPECTRAL-DENSITY §4).
- **KL representation** ⇒ any positive measure; admits the counterexample verbatim. No further constraint.

## 3 · The trace anomaly — the sharpest shot, and why it is BLIND

At the pole the trace Ward is **homogeneous in `f₁²`**: `f₁²⟨G|T^μ_μ|G⟩ = (β/2g)f₁² g_GGG ⇒ f₁² cancels` (INEQ-T4). It is an *identity* relating `g_GGG`, **not** an inequality on the edge — *zero* constraint on `f₁²`. Worse, the dilaton soft-scaling `⟨0|T^μ_μ|G⟩ ∝ m_G²` gives `f₁² ∝ f_G² m⁴`, so the light atom's χ-contribution is `f₁²/m² ∝ m² → 0`: **light states are χ-invisible for the anomaly operator.** No `λ` in `χ≤1/λ` can ever detect them. (Imposing even this scaling, `w=cε²`, keeps `χ→1/S` bounded while the atom sits at `ε→0`; the residue–level-gap lemma reports `f₁² ≥ cε² → 0`, exactly tight — the weight *is* `<Δ` and *shrinks*.)

| `ε=m²` | naive `w=ε/λ`: `χ` | dilaton `w=cε²`: `χ` | gap `√ε` |
|---|---|---|---|
| 1e-4 | 1.010 | 0.0101 | 1e-2 |
| 1e-8 | 1.010 | 0.0100 | 1e-4 |
| 1e-10 | 1.010 | 0.0100 | 1e-5 |

## 4 · Acid test — free massless must survive

Free massless, `O=φ`: `ν=Zδ_0`, `χ=∫dν/s=∞`. (a) `χ≤1/λ` is **vacuous** there (hypothesis fails) — correct behaviour. (b) it **owns** an atom at `0`, so any *unconditional* "no near-zero atom" is **FALSE** on it. (c) The conditional version (`χ<∞ ⇒ no sub-Δ atom`) is refuted by the finite-χ generalized-free measure of §1. So the exclusion is either vacuous (the open band `(0,Δ)` lets the exact atom at `0` slip through — a gapless theory) or false. The acid test kills it: it cannot distinguish gapped from gapless.

## 5 · Exact missing step + flip condition

The forbiddance must be a *lower* bound surviving `m→0`: uniformly in `(β∈[β_SC,2.5], a→0)`,
  **`Λ_β − R̄_β ≥ κ > 0`**, `R̄:=M₁/M₀`, `Λ:=inf supp(ν−f₁²δ_{m₁²})` —
equivalently a uniform-in-`m` **lower** bound on the trace decay constant `f_G` (one-state dominance / spectral saturation of the dilaton). This is verbatim `H*` = `I2` (GAP-G2-edge residue–level-gap lemma; PROOF-F1SQ "exact missing inequality"; G2-ASSEMBLE-100 §5). The anomaly's scaling is an *upper* bound `f₁²≤Cm⁴`; the missing step asks for the opposite sign. **Same wall.**

## 5-point verdict

1. **Nonperturbative** — the candidates are; but none *excludes* (as an exclusion it is false). ✗
2. **Controls SIGN / lower-bounds edge** — **FAIL**: positivity Schur-concave; anomaly gives upper `f₁²≤Cm⁴`.
3. **Non-circular** — **FAIL**: only confinement/gap forbids it.
4. **Fails free field honestly** — **FAIL**: vacuous or false (§4).
5. **Uniform `m_eff²(k→0)≥λ`** — **FAIL**: edge/atom coordinate ≠ curvature; no uniform `λ`.

**1/5. NOT `5/5`.**

**Sharpest line.** The tiny-weight near-zero atom is *explicitly realizable* as a generalized free field (`ν≥0` is all OS needs), and the trace anomaly does not touch it: the Ward cancels `f₁²` at the pole (identity, not inequality) and the dilaton scaling makes it `∝m⁴` — i.e. the anomaly does not *forbid* the light state, it makes it *invisible to the one functional (`χ`) the program measures*. Excluding it requires a level-`(n+1)` isolation lower bound `Λ>R̄` — a two-scale separation = `H*`, the upper-gap rung — obtainable only from the binding dynamics, not from RP, KL, OPE, or the anomaly.

*Refs.* REDUCTION-AUDIT (§b: `χ` is a harmonic mean, blind to `inf supp`); GAP-G2-edge (residue–level-gap lemma; RP upper-only); PROOF-F1SQ ("`Λ−R̄≥κ`, strictly above the gap"); PROOF-DISSOLUTION-SPECTRAL-DENSITY (`:φ²:` template, no YM selection rule); G2-ASSEMBLE-100 (§5 `I2=H*`); INEQ-T4 (trace Ward homogeneous in `f₁²`); GOLDS-POS-MASS (positivity pins only `m=0`); FORGE5-BOOTSTRAP-EXCLUDE (5-pt rubric, Schur-concavity). External: Jaffe, *Phys. Rev.* **158** (1967) 1454 (generalized free fields, `ν≥0 ⇒ OS`); Källén 1952 / Lehmann 1954; Osterwalder–Schrader, CMP **31** (1973) 83; Caron-Huot–Pokraka–Zahraee, JHEP **01** (2024) 195 [arXiv:2309.04472] (positivity forces minimal residue → 0); Goldstone–Salam–Weinberg 1962; dilaton soft theorem (anomaly ⇒ coupling `∝m²`).
