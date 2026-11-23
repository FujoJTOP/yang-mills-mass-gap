# IN-5 — overlap + two certificates (H3 × 2)

**Precise statement.** Fix the RT (`Λ, C=C_RG=24π²/11`, Wilson, `b=2`, lattice, `a→0` then `L→∞`). With `λ_G(μ)=h^∨g²(μ)`, `1/g²=B ln(μ/Λ)` gives `λ_G(μ)≤λ_w ⟺ μ≥μ_up:=Λe^{C_RG/λ_w}`; and `q_i(μ)=inv_i/μ` (`inv_m=m`, `inv_σ=√σ`), antitone. IN-5 = for each `i∈{m,σ}`, ∃ explicit reach `(λ_w^i,q_c^i,δ_W^i,δ_S^i)`:
- weak (IN-3): `μ≥μ_up^i ⟹ cert_i≥δ_W^i`;
- strong (IN-4): `0<μ≤inv_i/q_c^i ⟹ cert_i≥δ_S^i`;
- **overlap:** `q_c^i ≤ q₀^i := inv_i/μ_up^i ⟺ inv_i/Λ ≥ q_c^i e^{C_RG/λ_w^i}`.

Two separate H3's (`μ*^i=inv_i/q_c^i` may differ; collapse only if `μ_up` is shared).

**True status.** NOT implied by IN-3∧IN-4. The machine-checked core `handoff_gap` (`Handoff.lean`, `GL7Dual.lean`) *consumes* `overlap`; it never produces it. IN-3/IN-4 (existence of positive constants) relate neither `q_c^i,λ_w^i` to the invariant `inv_i/Λ`. Using the cited `q=m/μ` with `m` the physical RG-invariant gap, the irreducible content is literally `Δ/Λ = m/Λ ≥ q_c^m e^{C_RG/λ_w^m} > 0` (σ analog). So **IN-5 is equivalent to a quantitative lower bound on the very gap the theorem concludes** — extra content, not a gluing technicality.

**Advance — pinning the ask.** Reduction (GL7 Prop 1): shared `μ_up` ⟹ `m/Λ ≥ e^{t_up}·max(q_c^m,(m/√σ)q_c^σ)`; generally two inequalities. Natural point `λ_w≈7.28`, strict 2-loop `t_up=2.08` (`e^t=8.0`), `m/Λ=8`, `√σ/Λ=2`: m-side met (`q_c^m≤1`=OS78); **σ-side needs `q_c^σ≤0.25` (ξ/a≈4; `0.143`/ξ/a≈7 at λ=6).** Proven area law gives only `q_c^σ≈1`. Even granting Prop B (`z≤e^{-c₀q²}`, `c₀=1→q_c=0.83`) is short; `q_c^σ=0.25` needs `c₀≈11`. So IN-5 asks of IN-4: extend the *same* signed-gas reach, string sector, `ξ/a≈1→≈4`.

**Flip + plausibility.** `IN-5 ⟺ m/Λ≥q_c^m e^{t_up} ∧ √σ/Λ≥q_c^σ e^{t_up}`. Flips: (i) larger `m/Λ,√σ/Λ` (lattice ~4–8, need ~8+); (ii) `q_c^σ→0.14–0.25`; (iii) `κ_c<1` enabling `CPB-S^σ`. Verdict: with proven reaches IN-5 fails for σ by ≈4×; it is borderline-feasible only if the strong side reaches `ξ/a≈4`, and rests on the circular `m/Λ` input. IN-5 is not provable without a quantitative gap bound — it is where the number lives.
