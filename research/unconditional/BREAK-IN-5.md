# BREAK-IN-5 — σ-channel budget: the exact wall

**Exact reduction.** GL7's collapse (`GL7Dual.dual_overlap_iff`) gives, per channel,
`max(q_c^m, c·q_c^σ) ≤ R_m e^{−t_up}`, with `R_m=m/Λ`, `R_σ=√σ/Λ`, `c=m/√σ`, `t_up=C_RG/λ_w`.
Natural point: `λ_w=7.28`, strict 2-loop `t_up=2.081`, `e^{t_up}=8.009`; `m/Λ=8`, `√σ/Λ=2`.
Then `R_m e^{−t_up}=0.999` (m-threshold) and `q₀^σ=R_σ e^{−t_up}=0.2497`, `c=4`.

- **m-side:** `q_c^m ≤ 0.999`. Met by OS78 (`q_c^m=1`) and by Prop B at `c₀=1` (`0.833≤1`).
- **σ-side:** `c·q_c^σ ≤ 0.999`, i.e. `4·q_c^σ ≤ 1`. **This is the failing inequality.**

Prop B (`z(q)≤Ce^{−c₀q²}`, crossing `z=½` at `q_c^σ=√(ln2/c₀)`) makes the whole budget one number:

> **c₀ = T₂/m² ≥ ln2/(q₀^σ)² = ln2·(e^{t_up}Λ/√σ)² = 11.11.**

**Route A (reach ξ/a≈4, q_c^σ≤0.25 → c₀≥11.1): fails.** On the q²-scaling route this is dead
twice: BL-5 (`q²−q|_{q=1}=0`) gives zero gain at the crossover, so only absolute `σ(1)` enters;
and the sole theorem-bearing convergence (Kotecký–Preiss) holds for `β<β₀≤1/48` (SZZ 2204.12737)
while `ξ/a=4` sits at `β≈2.47`, ~100× outside. The q-linear (domain-wall) branch gives `z~1/q` —
worse, not better.

**Route B (bound c₀ from above → unsatisfiability): impossible.** BL-3 proves `σ(q)≥c₀q²` is
*strictly dynamical*: the `x₄=0` slab is a closed 2-chain with q-independent per-block cost, so
`c₀` may be 0; no `c₀` lower bound and no `T₂≤C·m²` upper bound exist. Reflection positivity gives
only lower bounds on `T₂`. So neither `c₀≥11` nor `c₀≤11` is provable. (Rescue via `κ_c<1` needs
`κ_c≲0.2`, BL-6 — no mechanism; default `κ_c≈1`.)

**The exact wall — which fails and by how much.**
`c·q_c^σ = 4·√(ln2/c₀)` must be `≤ 1`. With `c₀=1`: `4·0.833=3.33 > 1` — **fails by 3.33× in the
LHS** (`= (0.833/0.25)`), i.e. **11.1× in c₀** (`=3.33²`), i.e. reach must go `ξ/a:1→4`. The
`c=4=m/√σ` factor is the entire margin. The best physically expected `c₀` is O(1) (`T₂~m²`) or
`(Λ/m)²≈0.016` (`T₂~Λ²`) — both ≪ 11.

**Verdict.** No piece closes. Within all current inputs the σ-budget is **unsatisfiable**: it needs
`c₀=T₂/m²≥11.1`, and no mechanism supplies more than O(1); the two would-be sources (`σ(1)` large,
`κ_c<1`) are closed by BL-5/BL-6. It is *not* a theorem that `T₂≤11m²` for the true theory (`c₀` is
an unconstrained dynamical ratio), so "unsatisfiable" means *unsupplied*, not *refuted*. The wall
is the single dimensionless number `T₂/m²`, short ≈11×.
