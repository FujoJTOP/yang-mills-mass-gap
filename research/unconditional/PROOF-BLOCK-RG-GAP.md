# PROOF-BLOCK-RG-GAP — block-level RG: the right variable, the wrong regime at β=2.5

> **Verdict.** The idea fixes exactly what IR-H1 diagnosed. The block link after decimation, `U_b = P_b U`, **is** the coarse mode (the `∇∥` content), so a block-level gap feeds the **coarse slot `λ`** of the exact split — it would *supply H1*, not rename it: the first scheme whose gap lives on the right side. But it **does not close at β=2.5**. One block-decimation does **not** reach the strong-coupling cluster regime `u_eff·K<1`. Asymptotic freedom runs slowly: β=2.5 ⟺ g²=1.6, `u=I₂/I₁(2.5)≈0.51`, and one factor-2 IR step gives `g²_eff≈1.76`, `β_eff≈2.3` — four times above the heuristic threshold `β_eff≲0.6` (where `u=I₂/I₁=1/7`) and two orders above the rigorous `β<1/48`. The missing estimate is precisely **an RG contraction from the scaling window into the P3 domain**, which is the deep-IR flow = the gap (IN-9). Net: a RENAME on the **correct** variable — strictly better than PROOF-H2-GLOBAL — not a closure.

## The map (Migdal–Kadanoff / Fisher)
Use a covariant block kernel `T_b(U_b|U)` with `post(U)=U`, `U_b=P_b U` a block-averaged Wilson line, block factor `b=2` (MK bond-moving then one-plaquette decimation; or Fisher's controlled real-space step). Push forward to the block measure `μ_b ∝ e^{−S_eff(U_b)}`. The rigorous object is the RG map `R_b: β ↦ β_eff` measured in the strong-coupling polymer norm `‖·‖_SC` whose unit is the plaquette activity `u_eff·K`. Convergence of the P3 (polymer) expansion needs `‖R_b(β)‖_SC < 1`.

## Does one step reach `u_eff·K<1`? **NO.**
At β=2.5, `u≈0.51`, while `u·K<1` with `K≳2d−1=7` needs `u≲1/7`. One step moves `g²` by `O(g⁴)=O(0.16)` (asymptotic freedom: `dg²/dlnμ = −(11N/3)g⁴/(8π²)`), giving `β_eff≈2.3`. The threshold `u=1/7` sits at `β_eff≈0.58`; reaching it takes `≈9–11` factor-2 steps (MK does not beat this order) — i.e. coarse-graining to the confinement scale `a_b∼1/m`. **The threshold is crossed only in the deep IR, and proving the flow reaches it *is* proving confinement** (IN-9). One step is two orders short of the rigorous domain.

## The RG inequality exists — and it captures `∇∥` (this is the real gain)
Disintegration is exact: `Var_μ(F)=E_{μ_b}Var_{fib}(F|U_b)+Var_{μ_b}(E[F|U_b])`, with `|∇F|²=|∇∥F|²+|∇⊥F|²` (PROVE-H3). So
`Δ_fine ≥ min(Λ_UV, ε_b)/(1+2S/ε_b)`, `ε_b =` Poincaré constant of `μ_b`,
the REPAIR-2SCALE `m=min(λ/2, Λ_UV/(1+2S/λ))` with **`λ=ε_b` supplied by the block gap, not assumed**. `∇⊥` is the *proved* `Λ_UV`; `∇∥` is captured by `ε_b` because `U_b` is the coarse mode. **A rigorous `ε_b>0` at β=2.5 would prove the gap.** The construct is sound; only its hypothesis fails.

## Why the block variable *is* the obstacle
Coarse-graining preserves the correlation length `ξ=1/m`; `ξ/a_b` shrinks by `b` each step. Cluster convergence needs `ξ/a_b≲1`, reached only when `a_b∼ξ`, i.e. at the mass scale. The low-pass keeps exactly the low-`k` modes — the soft mode `G″(0)≈−0.165` (VAR-FINAL) lives **in `U_b`**. The block variable captures `∇∥` for the identical reason it inherits the soft mode that defeats the expansion.

**Sharpest line.** The block-RG idea answers IR-H1 at the structural level (right variable, and the fine⟸block inequality `Δ_fine≥min(Λ_UV,ε_b)/(1+2S/ε_b)` is exact). But the only missing input is `ε_b>0` at β=2.5, and that requires a one-step contraction `‖R_b(2.5)‖_SC<1` into the P3 domain that asymptotic freedom forbids (`β_eff≈2.3≫0.6`). The remaining estimate — an RG contraction to strong coupling — is the deep-IR flow: the gap, relocated one scale out on the correct mode. A RENAME, not a closure; gap #7 stands.

*Refs.* Migdal 1975; Kadanoff; Fisher, JSP 1995 (real-space RG); Osterwalder–Seiler 1978; Shen–Zhu–Zhu arXiv:2204.12737 (β<1/48); Balaban CMP 109/116/122; PROOF-IR-H1; PROOF-H2-GLOBAL; H2-globalization; REPAIR-2SCALE; PROVE-H2-uvgap; PROVE-H3; VAR-FINAL; IN-9; W2-regularity-RG.
