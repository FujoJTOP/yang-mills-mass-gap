# IN-7 — Non-triviality (non-Gaussianity)

**Verdict.** IN-7 is **not closed**; `conditional-theorem.md`'s "trace anomaly + IN-0" overstates it. Proven is only the reduction `IN-7 ⟺ N1 ∧ N2`; the load-bearing unproven step is **N2 (`Λ>0`)**.

**(1) Precise statement.** Three candidate meanings; only the third is IN-7:
- *Not concentrated on one trivial vacuum* (`μ ≠ δ_orbit`) is **non-degeneracy**, supplied by MG/IN-6 — **not** IN-7's content.
- *Wilson area-law `σ>0`* is an **alternative** route (E2-σ), gated by the same IR wall (I-CONF) — not its meaning.
- *Non-Gaussian*: `¬IsGaussianMeasure μ` (`¬IsGeneralizedFreeField`). By Lemma C this is "some `n≥3` cumulant ≠0"; the free-vanishing gauge-invariant witness is `g_ren≠0` (Lemma G), equivalently `Λ>0` (Lemma DT). "`g_ren/Λ>0`" is shorthand for `g_ren≠0 ∧ Λ>0`.
- Connected **4-point ≠0 is vacuous** (free `tr F²` is a Wick square); even for even measures `κ₄=0 ⇏ Gaussian` — `κ₆≠0` (Lean `four_point_insufficient`, no axioms). Need every `n≥3`.

**(2) True status.** Unconditional (core Lean): Lemma C, the `κ₄=0/κ₆≠0` counterexample, Rademacher `κ₄≠0`; paper-level Lemmas G, S (scheme-independence), DT. **Theorem N** (`N1∧N2 ⟹ non-Gaussian`) is complete. **N1** (short-distance `g_ren` converges, matches one loop) is nearly in IN-0: `Y3`,`Y4`, and assessment §7.4 give `g_W(r)>0`, `log(1/r)g_W(r)²→1/(2b₀)`, though limit-independence is unstated. **N2 is absent**: IN-0 is short-distance only ("no intrinsic scale or spectral gap asserted") and never excludes `g₀≡0`.

**(3) Advance — exact missing step.** The quoted IN-0 content (`T^μ_μ=(β/2g)trF²`, one-loop `b₀`) holds for the Gaussian theory with both sides zero (`β=0`); it **cannot** entail `g_ren≠0`. So IN-7 reduces to one limit-selection statement: **the subsequential limit is not the free trajectory.** That is N2 — an IR/scale-fixing input, not UV.

**(4) Flip condition.** IN-7 → theorem iff one proves `Λ>0` for the limit (excludes free collapse; physically, the construction fixes a nonzero renormalized coupling at a finite scale). IN-7 → **false/vacuous** iff the limit is Gaussian (`Λ=0`, `g_ren≡0`, all `n≥3` cumulants vanish). The `σ>0` route reaches the same conclusion only through the separate IR wall — not a flip.
