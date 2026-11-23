# Direct glueball 4-point bootstrap — assessment

> 2026-11-23 · Probe of one specific untried route: SDPB on the actual
> `⟨tr F²·tr F²·tr F²·tr F²⟩` (or `⟨T T T T⟩`) crossing, **not** the defect bootstrap.
> **Verdict: ill-posed as stated — pure YM is not a CFT — and even its well-posed shadow lands in
> the *same* anomaly-protected 0⁺⁺/2⁺⁺ sector, so it does not lift IN-11 and cannot supply the
> lower bound the gap needs. Not worth running in this form.**

## 0. The fatal premise: no conformal crossing exists
Conformal crossing, conformal blocks, and SDPB all presuppose conformal symmetry. Pure 4D YM has
none in the IR: it confines, and `T^μ_μ = (β/2g) tr F² ≠ 0` explicitly breaks scale invariance. Its
only fixed point is the UV Gaussian one (asymptotically free; no interacting CFT anywhere in the
flow). So `⟨tr F²⟩⁴` and `⟨TTTT⟩` are **not** conformally covariant and there is **no crossing
equation to solve**. Imposing conformal crossing on real YM silently assumes symmetry the theory
breaks; whatever SDPB returns describes a fictitious CFT, not the glueball.

## (a) Crossing + OPE — meaningful only in a hypothetical CFT
For `O = tr F²` (0⁺⁺) the OPE `O×O` runs over the identity, `T_μν` (Δ=4, conserved), and a tower of
**0⁺⁺, 2⁺⁺, 4⁺⁺ …** (even spin, C=+ P=+). The identity and T blocks are "protected" — coefficients
fixed by the Ward identity / central charge — the rest unknown. The selection rule is the crux:
`O×O` is **even parity**, so **0⁻⁺ never appears.**

## (b) Can SDPB bound the lightest 0⁺⁺ dimension?
Only a weak bound, and the wrong object. In a genuine CFT `tr F²` at a fixed point is exactly
marginal (Δ=4) and `T` is conserved (Δ=4): **both external dimensions are pinned**, so the usual
`(Δ_σ, Δ_ε)` island scan has no free parameter. What remains is a scan over the *assumed exchanged
gap*, whose feasible edge is a **one-sided bound / island on a fictitious CFT operator** — not the
YM glueball mass (in a gapped theory "dimension" and "mass" are different objects). Freeing `Δ_O`
instead just bootstraps an abstract scalar CFT, not YM.

## (c) Does the direct 4-point lift IN-11? — NO; it sees only the protected sector
IN-11's content is that the lightest bulk state sits in the C-even P-even **0⁺⁺/2⁺⁺** sector (weak
form b′: no lighter 0⁻⁺). The direct `⟨tr F²⟩⁴` 4-point is exactly `0⁺⁺×0⁺⁺` → **the same
even/even-parity sector the defect bootstrap already covers**, and it is equally **blind to 0⁻⁺**.
`tr F²` and `T_μν` are precisely the anomaly-sourced operators; the 4-point sees *all* even-parity
channels (more than a single protected one) but **only anomaly-protected ones**. The unprotected
0⁻⁺ — the channel IN-11 turns on — is invisible; reaching it needs `⟨tr F F̃⟩` correlators that are
total-divergence / contact-term ambiguous (exactly the obstruction in `GAP-IN11.md`). So the direct
route **consumes** IN-11 rather than proving it.

## (d) Genuine bound, or only island/upper?
Island/upper only — and for actual YM not even that. This is `III-UNIQUENESS` restated: SDPB gives
pointwise bounds and islands, never global exclusion or existence; "more points" never reaches the
uniqueness theorem that is absent for 4D CFTs. The gap is a **lower** bound (lightest mass ≥ Δ > 0);
crossing bounds run the other way. Even the well-posed cousin — the 4D **S-matrix** bootstrap of
glueball scattering — must *assume* the gap to define asymptotic states and returns allowed
regions / upper bounds on couplings, never a proof.

## Worth running? What would it give?
**Not as stated.** A conformal SDPB crossing on YM is the wrong EFT. Forced, it yields at most a
one-sided CFT bound/island, inside the **same** 0⁺⁺/2⁺⁺ sector the defect bootstrap already reaches —
no new channel, no lift of IN-11, wrong direction, and not a genuine YM observable. The gap stays
framework-invariant here too: G2, G3 and U are untouched.
