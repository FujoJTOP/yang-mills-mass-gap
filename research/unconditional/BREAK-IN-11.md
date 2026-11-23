# BREAK-IN-11 — can `m_{0⁻⁺} ≥ m_{0⁺⁺}` (a pseudoscalar lower bound) be proven?

**Verdict. No.** Every route gives an *upper* bound on `m_{0⁻⁺}` or an approximate value; none is a theorem.
`m_{0⁻⁺} ≥ m_{0⁺⁺}` stays physical input.

**Direction obstruction (rigorous).** Let `ν` be the Källén–Lehmann measure of a gauge-invariant interpolator in
channel `c`, support ⊂ `[m_c²,∞)`. Positivity gives every moment `M_k=∫sᵏdν ≥ 0`, so every positive moment
*ratio* is a mean ≥ `m_c²`: `m_c² ≤ (M_a/M_b)^{1/(a−b)}`. Any positivity-only structure — dispersion, KL,
susceptibility, sum rules — bounds the lightest mass from **above**. A *lower* bound concerns the spectral *edge*
(`f₁²=ν({m_c²})` bounded below); positivity bounds the *total* measure and cannot resolve the edge. The same wall
blocks the weak form (a sub-`δ` bound in every non-`0⁺⁺/2⁺⁺` channel): the edge problem replicated per channel.

**Routes.**
- *Variational / min–max.* `E₁ ≤ ⟨φ|H|φ⟩/⟨φ|φ⟩`: a genuine **upper** bound on `m_{0⁻⁺}` for any trial operator.
  Rigorous lower bounds (Temple; Kato–Weinstein) need a lower bound on the rest of the spectrum (`E₂`, the gap
  below) — the target itself. Circular. Blocked.
- *QCD sum rules, "controlled errors."* Not theorems: they assume resonance dominance, continuum threshold `s₀`,
  OPE truncation, Borel stability; "errors" are estimates, 4D duality unproven. Exact content = moment lemma ⇒
  **upper** bounds only.
- *Topology.* `q=tr F F̃=∂·K` and KL: `χ_t=Σf_n²m_n² ≥ f₁²m_{0⁻⁺}²`, i.e. `m_{0⁻⁺}² ≤ χ_t/f₁²` — **upper** bound
  (same lemma; cf. Veneziano `m²∼χ_t/f²`). Anomaly fixes the pseudoscalar *to* `χ_t`, never below.
- *Large-N.* Single-trace masses `∼N⁰`; `m_{0⁻⁺}/m_{0⁺⁺}→const` is `N`-independent but the constant is dynamical
  (master field/string unsolved). No inequality, no ordering.
- *Bootstrap on `0⁻⁺`.* The defect bootstrap is fed by `tr F²`, `T_{μν}`, whose normalizations are protected.
  `tr F F̃` is a total divergence: contact-term/scheme ambiguous, coupling `f₁` unprotected (may decouple), not
  stress-tensor-sourced. The bootstrap cannot see `0⁻⁺`; a lower bound would need `f₁` bounded *below* — the
  unprotected quantity (`TWENTY-PATHS` E3: bootstrap not a proof).

**Genuine theorems only:** Rayleigh–Ritz upper bound; the positivity moment-ratio lemma (⇒ upper bound from any
KL/dispersion input); the `χ_t` identity; trace anomaly ⇒ channel non-empty (existence, no ordering).

**Sharpest line.** `m_{0⁻⁺} ≥ m_{0⁺⁺}` needs a *lower* bound on a spectral edge. Every machine either bounds the
total measure (⇒ upper), lacks a protected edge coupling, or presupposes the gap below. None is a theorem; IN-b
remains physical.
