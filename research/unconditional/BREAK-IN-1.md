# BREAK-IN-1 — marginal gauge-covariant BPHZ, sharpened route: the scale sum closes; Savvidy does not

**Exact target.** J=⌊log(1/aΛ)⌋, g_j²=1/(2β₀j)+O(j⁻²), β₀=11N/(48π²). UVS(L) needs the factorial bound |S^{(n)}_a(f)| ≤ K^n n!‖f‖^n uniformly in a (Minlos). Sharpened route: track one coupling g_j Tr F² to β₀ exact, keep R_j = dim>4 + higher-loop residual, demand the W4 §1 trapping lemma
(T1) 0<g_j≤g*, ‖R_j‖_X≤r j^{-θ};  (T2) ‖R_{j+1}‖_X ≤ e^{-c}‖R_j‖_X + Cg_j⁴;  (T3) g_{j+1}=g_j−½β₀g_j³+O(g_j⁵), uniform in cutoff and volume.

**Exact scale-sum estimate (the thing to prove).** Uniformity of (T2),(T3) over J reduces to two harmonic sums:
(S1) Σ_{j≤J} g_j⁴ = (2β₀)⁻⁴ Σ j⁻⁴ ↑ ζ(4)/(16β₀⁴) < ∞  [computed: → 2.33×10⁵];
(S2) coupling error δg_J ≤ Σ_{j≤J}(j/J)^{3/2} g_j⁵ = O(J^{-3/2} log J) → 0.
The naive wall Σg_j² ~ (1/β₀)log J **never appears**: β₀ tracks it exactly, and every higher power g_j^p (p≥4) is absolutely summable. That is the entire content of log-resummation — the marginal log log is *absorbed into g_j²*, not summed (W4 line 12). No 4-loop coefficients needed; β₀ suffices for a-uniform boundedness.

**Block-by-block.** Block b = scales [bΔ,(b+1)Δ): marginal residual Σ_block g_j⁴ = O(b⁻⁴) (summable); irrelevant dim>4 contract by e^{-cΔ}; coupling tracking O(b^{-5/2}) (summable). Σ_b b⁻⁴<∞, Σ_b b^{-5/2}<∞ ⟹ sup_J‖R_J‖_X<∞. **The scale sum closes. No log power and no tree breaks it.**

**Exact break — Savvidy, not combinatorics.** (T2)'s factor e^{-c} is derived only from the free Gaussian. At d=4 the marginal Tr F² yields no amplitude decay, and on the Savvidy negative-mode cone the form TrF²+R_j is not positive: the lowest fluctuation eigenvalue is ∝ −g|F| < 0 (tachyonic), at *every* scale j, not one. There e^{-c}→e^{+|c|}, so a single step's large-field uniformity fails. The obstruction is a measure-convexity / large-field bound at the one step, certified by that negative Hessian eigenvalue — **not a log power, not a BKAR tree**. Hairer's degree condition (needs d<4) and Balaban's uniformity both stop at exactly this step.

**Verdict.** The sharpened route kills the log-sum fork W3 flagged — Σg⁴ converges, the log log is not a wall. The honest residual is the Savvidy large-field step: analytic, all-orders, present at every scale. Not closed; and specifically not a log.
