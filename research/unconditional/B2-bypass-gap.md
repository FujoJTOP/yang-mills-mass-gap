# B2 — BYPASS: any route to the Clay statement that need not prove the gap?

> 2026-11-23 · **Verdict: no bypass.** Clay is `AX ∧ NT ∧ GAP` with `GAP` a conjunct, so soundness forces every proof to derive `Δ>0` for the constructed `μ*`. The only honest form of the question — "can `GAP` be reached *without a lower bound on the spectral edge* `f₁²`/H1/`λ`?" — is answered NO on all four routes, each collapsing to a named dependency. The pure-YM-specific residue is problem (iii), already ruled unclosable (`III-UNIQUENESS`). No fabricated proof.

## Logical skeleton (and the conjunction lemma)
`Clay(μ*) = AX(OS0–4) ∧ NT(Λ>0, non-Gaussian) ∧ GAP(Δ>0) ∧ [Wilson construction]`. `GAP` is literally a conjunct of the conclusion, so "a reduction avoiding the gap as an intermediate" is impossible at the level of *derivation*: any sound proof entails `GAP(μ*)`. (Trivial, but it pins the real question: **is every route to `Δ>0` a lower bound on `inf supp ν`?**)

## (a) A different axiom "acceptable as existence of YM"
Every *construction-level* axiom is satisfied by a gapless witness:
- `AX` + clustering + unique vacuum ⇒ the **massless free scalar** — OS0–4, power-law `⟨φ(x)φ(0)⟩=1/(4π²x²)→0`, unique vacuum, `Δ=0` (Osterwalder–Schrader). So even `AX ∧ OS4 ∧ unique vacuum ⇏ GAP`.
- `AX` + NT + `Λ>0` (trace anomaly) ⇒ a **Banks–Zaks** flow: asymptotically free, `Λ>0` by transmutation, non-trivial, IR-interacting CFT, `Δ=0`.

Hence a construction-level `A` does **not** imply GAP. Strengthen `A` until it does, and the strengthening *is* "no massless transverse edge" = `f₁²>0` = GAP. "Acceptable as existence of YM" equivocates between OS-only (⇏gap) and OS+gap (assumes the conclusion); W1/W3 leave no middle. For a **pure-YM-specific** `A` the gapless alternative is the conformal/Coulomb IR = problem (iii), open.

## (c) Does the construction force the gap?
`NT ∧ OS ∧ clustering ⇒ GAP`? **False.** W1 (free) and W2 (interacting CFT: rigorous in 2D minimal models; N=4/3D-WF for the gauge type) are non-trivial, OS, power-law-clustering, gapless. The pure-YM instance — "the measure is not an interacting 4D CFT" — is exactly (iii): `III-UNIQUENESS` shows no 4D CFT-uniqueness theorem exists and SDPB yields pointwise islands, not global exclusion. So (c) = (iii) = the gap. `TWO-ATOMS` states the same: existence ⊥ gap for fixed-volume OS; entanglement enters only through the construction residue.

## (d) Clay is `Λ>0` existence, not a quantitative bound
Two real facts, both insufficient. (i) The trace anomaly fixes `Λ>0`, never `Δ>0`: it forces states to exist, not to be heavy (`GAP-IN11` §4); BZ has `Λ>0 ∧ Δ=0`. (ii) Bare `Δ>0` is **strictly weaker** than the project's uniform `f₁²≥r₀` — the edge may be an *atomless* threshold (`f₁²=0`, `Δ>0`) — so `RESIDUE-CONDITIONAL` is sufficient, not necessary (`W4`). But the obstruction to `Δ>0` is a lower bound on `inf supp ν`, the **same** sign-definite wall as the residue: every linear/positivity functional gives *upper* bounds; forbidding the massless edge needs a Schur-convex coercivity (total ≤ atom), absent from RP/Ward/anomaly/GZ (`EDGE`, `LEVER`). Uniformity finishes it: `Δ_∞=liminf Δ_{L,a}>0` ≡ a uniform `δ>0` along the RT (`IN-9`, `P1-L0`). (d) is the right *target*, not a method.

## (b) Spectral dichotomy
"`GAP ∨ X`, then prove `¬X`" helps only if `X` is decidable without edge control **and** `¬X` excludes a massless mode. But `¬GAP` is *exactly* `0∈supp ν` (massless atom or continuum); the dichotomy is the excluded middle on `0∈supp ν` — a tautology. To be useful, `X` must be a massless phase that an anomaly *forbids*, and anomalies constrain phases without producing a gap (`T2`(iv), `GAP-IN11`). So (b) either decides the edge (= gap) or names the forbidden structure "massless transverse mode" (= gap).

## Bypass-impossibility lemma
**Lemma.** Let `A` be a hypothesis having a model in the construction class with `Δ=0`. Then `A ⇏ GAP`; no sound derivation of Clay uses `A` alone, and any *sufficient* `A'` must already exclude the massless edge, i.e. assert `f₁²>0` — gap-strength.

*Proof.* `GAP` is conjunctive in Clay (soundness). The witnesses W1 (rigorous), W2 (rigorous in 2D), W3 (AF-with-IR-fixed-point) satisfy every construction-level axiom (OS, RP, power-law clustering, anomaly/`Λ>0`, unique vacuum). An `A` they satisfy cannot entail `GAP`. ∎

*Honest dependency.* At the construction level the lemma is a **theorem** (W1 suffices). At the **pure-YM-specific** level the witness is not a free model but the conformal/Coulomb IR alternative; excluding it is problem (iii), itself open. So the lemma reads: at the first level the bypass is **false**, at the second it **is** the gap.

**Bottom line.** `λ>2C_UV` (`SCORE-HYPOTHESIS`), `f₁²≥r₀` (`RESIDUE-CONDITIONAL`), the "one door" (`LEVER`), and "construction ⊥ gap" (`TWO-ATOMS`) are one object: the `k→0` transverse edge. `Clay ⟺ [IN-0…IN-9] + (f₁²>0)`, and the parenthetical is `IN-10/IN-11` ≡ the gap.
