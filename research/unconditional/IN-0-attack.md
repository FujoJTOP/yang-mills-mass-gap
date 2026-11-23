# IN-0 attack

**Restatement.** IN-0 = the 2026 construction's 16 inputs (hexagonmath `hexagon:2611.00492v1`, not an arXiv id): RG1–RG5, UV1, PH1–PH3, AN1, Y1–Y4, B1, C1 — minus the assumed `sym` field (`submission/YM/Phase/Mixture.lean:142`; consumed at `checks/SourceScalar.lean:34`), reassigned to IN-2. Jointly they give OS E0–E4, a symmetric conserved stress tensor with Ward identity, the trace anomaly, one-loop `b₀=11h^∨/(48π²)`, subsequential fixed-volume UV and infinite-volume limits, OPE bounds, and the lattice bridge. They contain no gap, uniqueness, or confinement (`16-input-audit` §2).

**True status: proven-for-existence / negative-for-gap — not "certified external inputs."** The committee certified the OpenAI existence theorem, not hexagonmath's 16 labels (`16-input-audit` §1,§6). Tally:

- Discharged in the unrefereed, library-pinned 9,470-module build or by published imports: RG1–RG3, RG5, UV1, PH1–PH3, AN1 (E0–E3 theorems, E4 via `CompleteOS`), Y1 (§7.1), Y2 (eq.13), Y4 (§7.4), B1, C1 — **14**.
- Derived: RG4 (RG2+RG3) — not independent.
- **Open residue: Y3**, limit clause only. "The submission does not prove `g_W(r)` is independent of the subsequential infinite-volume limit" (§7.4); OPE coefficients' limit-independence likewise unproved (§7.5); "not all-orders perturbative matching" (`checks.md` A4).
- Route-negative: an explicit OS model satisfies all 16 with clustering yet has no spectral gap (prior-work §3.4). So IN-0 ⇒ OS theory only.

**Attack on Y3.** Rigorous piece: short-distance data are functionals of the limit state μ. Linear data (the trace amplitude) are affine, hence phase-mixing-invariant if every pure phase agrees; the nonlinear `g_W` (inverse √ of a blocked coefficient) needs genuine phase uniqueness. μ=μ′ ⇒ all data agree; a singleton limit point makes Y3 automatic. Exact missing step = **uniqueness of the subsequential limit on the gauge-invariant local algebra** (IN-6; assessment Open Q1). AN1's qualitative E4 cannot supply it (Ising phases cluster yet differ locally); a uniform rate would.

**Flip / plausibility.** Discharge = unique infinite-volume vacuum (or a uniform clustering rate). Likely true — strictly weaker than a gap, universally expected — but no known 4D route; IN-6 is itself conditional (W1₀/W1₁). Y3 is genuinely open but weakly load-bearing; its flip is the existing IN-6 wall, not a new reduction.
