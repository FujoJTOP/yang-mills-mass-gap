# INEQ-T4 — protected norm (a-anomaly / trace anomaly) vs. the individual 0⁺⁺ residue

**Verdict.** The 0⁺⁺ channel *is* the one anomaly-sourced scalar channel — so the intuition is sharp —
but the anomaly protects the **wrong object**. A Ward/STI identity at the pole **does construct** a
protected same-channel quantity: the cubic glueball coupling `g_GGG` (equivalently the forward trace
charge). It is **not** the one-particle residue `f₁²=⟨0|¼trF²|G⟩²`, which **factors out** of every such
identity. The edge weight stays unprotected *at the pole* — for a structural reason, not a missing computation.

**Protected vs. residue.** In pure YM the exact operator identity
`T^μ_μ=(β/2g)·¼trF²≡(β/2g)O_S` is the unique bridge into 0⁺⁺. Its protected contents are (i) the
coefficient `β/2g` (RG-invariant) and (ii) the **forward** matrix element `⟨G|T_μν|G⟩=2m_G²(…)` — the
equivalence principle, a **diagonal** bilinear. The KL residue `f₁²=|⟨0|O_S|G⟩|²` is instead
**off-diagonal**: a vacuum↔state decay constant.

**Construction — STI at the pole.** Factorize `Γ_μν=⟨O_S T_μν O_S⟩` and `⟨O_S O_S O_S⟩` at `p²=m_G²`:
each residue is `f₁·(vertex)·f₁`. The trace Ward identity `g^{μν}Γ_μν=(β/2g)⟨O_SO_SO_S⟩` gives
`f₁²⟨G|T^μ_μ|G⟩=(β/2g)·f₁²·g_GGG`, hence
`⟨G|T^μ_μ|G⟩=(β/2g)·g_GGG`.
So the anomaly **does** pin a same-channel object from the conserved-current normalization — but it pins
`g_GGG` (the cubic coupling / forward trace charge), and **`f₁² cancels identically**. Every Ward/STI at
the pole is homogeneous in `f₁²`: the residue multiplies both sides and factors out. An *identity*, not
an *inequality*.

**Why the residue escapes — three legs.**
1. *Diagonal ≠ off-diagonal.* Charges/anomalies fix same-state bilinears; a decay constant is
   vacuum↔state. Conservation `p^μ⟨0|T_μν|G⟩=0` fixes only the transverse tensor structure and leaves the
   scalar coefficient — hence `f₁` — free. A state can carry its protected charge while projecting
   arbitrarily weakly onto any fixed interpolator.
2. *Aggregate ≠ edge.* The UV object the anomaly controls — the `⟨T^μ_μ T^ν_ν⟩∼β²` short-distance
   coefficient, i.e. the `x→0` OPE coefficient — is an *integral/moment* of the whole measure `ν_S`.
   Positivity at fixed total cannot bound the lightest atom (moment lemma, localized to one channel).
3. *Homogeneity.* `f₁²` cancels at the pole (above): the identities relate products of couplings, never
   solve for the edge.

**a-anomaly.** It is the `E₄` coefficient in `⟨T^μ_μ⟩_g`, entering `⟨TTT⟩`/UV contact data — again
aggregate/short-distance, untouched by the IR pole. In the asymptotically free UV the `⟨TT⟩` norm is the
free central charge; protection is UV, the residue is IR.

**Honest verdict.** The 0⁺⁺ residue is the protected-norm attack's *best shot* — same channel as the
trace anomaly — and it still fails: the anomaly protects the normalization, the diagonal charge, and
`g_GGG`, **never** the edge overlap `f₁²`. **Flip condition:** a theorem forcing the trace decay
constant `⟨0|T^μ_μ|G⟩` (equivalently `f₁`) to saturate the protected forward charge — a
one-state-dominance / spectral-saturation statement. None known; it is tantamount to the IR solution.
