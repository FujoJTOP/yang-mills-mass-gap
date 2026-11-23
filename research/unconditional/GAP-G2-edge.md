# GAP-G2-edge — can ANY structure give a deterministic LOWER bound on the lightest mass or its residue?

> 2026-11-23 · Dedicated edge attack. **Verdict: no known structure. The direction wall is structural, and I can state the exact reduction that every candidate hits.**

**Setup.** A gauge-invariant interpolator in channel `c` gives, by reflection positivity, `ν ≥ 0`, `supp ν ⊂ [m_c²,∞)`, `f₁²=ν({m_c²})`, `M_k=∫sᵏdν`. The moment-ratio lemma (`M_a/M_b ≤ (m_c²)^{a−b}` for `a<b`, any finite moments, negative `k` included) gives `m_c² ≤ (M_a/M_b)^{1/(a−b)}` — always **upper**. A lower bound on `m_c` is `ν([0,M²))=0`; a lower bound on `f₁²` is the edge *weight*.

**Moment data cannot pin the edge (rigorous, elementary).** `M_0=M_1=2` is realised by both `δ_0+δ_2` (edge 0) and `2δ_1` (edge 1). Finitely many positive moments do not determine `inf supp`, nor the atom there. Every lower bound on the edge is therefore an independent, non-moment datum.

**Lemma (rigorous): residue ⟸ level gap.** Let the lightest state be an isolated atom at `m₁²` and `Λ=inf supp(ν−f₁²δ_{m₁²})>m₁²`. Then
`f₁² ≥ (ΛM_0 − M_1)/(Λ − m₁²)`, equality iff the rest is one atom at `Λ`.
*Proof.* `M_0=f₁²+ρ_0`, `M_1=f₁²m₁²+ρ_1`, `ρ_1≥Λρ_0`; eliminate `ρ_0`. ∎
So **a lower bound on the lightest residue is equivalent to a lower bound on the second level**, given the total moments. Edge control at level `n` needs level `n+1` — an infinite regress, or a proven continuum threshold. In YM the next level is the two-glueball threshold, a dynamical binding input. **Edge control ≡ a two-scale separation theorem, one rung above the gap.**

**The four candidates.**
- *Trace anomaly / stress-tensor norm.* `T^μ_μ=(β/2g)trF²` makes `0⁺⁺` non-empty; the `1/x⁸` coefficient of `⟨TT⟩` is `c>0`. But these are **total/UV** data — the lemma turns totals into *upper* bounds. For dim-4 operators the total moment diverges and only the anomaly (log) coefficient is scheme-independent; no finite protected moment exists. Even a finite `M_{−1}=∫dν/s` obeys `m₁²≤M_0/M_{−1}` — again **upper**.
- *Determinant / reflection positivity.* RP is a PSD condition ⇒ Cauchy–Schwarz/det bound matrix elements **from above**; PSD cannot manufacture an individual `≥`. RP is scale-covariant — RP theories with arbitrarily small gaps exist — so RP alone cannot force `m>0`.
- *Individual vs averaged OPE.* Protected individual coefficients are only identity, `T_{μν}`, conserved currents. The glueball residue is an individual coefficient to an **unprotected** operator; the averaged total is protected. The distinction sharpens the wall; it does not cross it.
- *Lattice–continuum.* The only rigorous gap is strong coupling, `β<β₀` (Osterwalder–Seiler gap; SZZ radius 1/48), `~10²×` below `β_c≈2.47`. A continuum residue bound needs uniform-in-`a` control **plus** the same `Λ`; neither exists.

**What edge control would require (exact).** A functional `Φ` on positive measures, built from protected data + provable bounds, with `Φ ≤ inf supp ν`, *not* implied by finitely many positive moments. A sufficient package: (a) `M_0,M_1` (protected norm + positivity ratio) **and** (b) a lower bound `Λ` on the excited threshold with `ΛM_0>M_1`. Equivalently: prove the lightest glueball is isolated with a gap to the next level — the mass gap shifted one rung up. No anomaly, RP, OPE or lattice structure supplies it unconditionally.

**Only theorems:** the moment-ratio lemma (upper), the residue–level-gap Lemma (relative, given `Λ`), anomaly channel-normalisation, and the strong-coupling lattice gap. Nothing gives a deterministic continuum lower bound.
