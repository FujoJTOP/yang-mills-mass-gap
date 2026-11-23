# HARDFORGE-CONSTRUCTIVE-FULL — forcing `G″≥δ` on the Savvidy cone from norm + tail + Jacobian; the one step that dies

> **Verdict (判词).** A complete step-by-step attempt: the exact SCBL identity and the geometric positivity are theorems, the coercivity-weighted jet norm and the large-field value-tail are banked constructions, and the *shallow*-cone fiber sign closes. Together they pin the value, the β=0 sign, and the shallow sign — but the passage to `G″=S−D≥δ` **on the whole cone** is the cone variance bound **(★)=H1**, and it fails there. The attempt is written to the wall, not past it.

**Step 0 — setting.** Fixed volume, RT line, blocks `j≤J`, `μ_j:=Λ_{j−1}`, shell `[Λ_{j−1},Λ_j]`; `K_j=(1/2g_j²)∫Tr F²+R_j`; trap `𝔗={0<g_j≤g*, ‖R_j‖_X≤rg_j², 2g_j|F|<μ_j²}`. Coercivity `Q[F]=∫Tr F²≥0`. Measure `μ∝e^{−V}dX`, `V=S_W+Φ`, `Φ=−log J_Haar`, `Φ″⪰c_Φ=N/12` (SU2 1/6, SU3 1/4).

**Step 1 (exact identity — theorem).** Mode `v`, `a=⟨v,X⟩`, conditioned `a=0`:
`G″_j(v)=E_μ[H_vv+Φ″_vv|0]−Var_μ(G_v+Φ′_v|0)=S(v)−D(v)`,
`S=E_μ[h_n−⟨K^{−1}c,c⟩|0]`, `D=Var(f|0)−E_μ⟨K^{−1}c,c⟩|0`, `K=Π₀ Hess V Π₀`, `c=Π₀ Hess V v`. Exact, no log-concavity assumed; `S` = Gaussian–Schur curvature, `D` = BL defect: `D≤0⇔K⪰0`, `D>0⇔K` indefinite (Savvidy).

**Step 2 (geometric positivity — theorem).** `Φ` convex, `g″≥1/6` ⇒ Prékopa ⇒ `W_geo≥0` on every mode incl. `k→0`; pointwise `Φ″⪰N/12`. The Haar fiber thus adds a strictly positive `N/12` to both `h_n` and `K`.

**Step 3 (coercivity-weighted jet norm — construction).** `‖R‖_X:=sup_{|v|=1}sup_A e^{Q_{μ_j}[F]/4g_j²}|δ²_{vv}R(A)|ℓ(A)`. It dominates the shell second variation `δ²K_j=g_j^{−2}(−D²+2adF_B)+δ²R_j`, so the trap's contraction/amplification is measured in the **Hessian**, not the value. Forced (TRAPPING §1).

**Step 4 (deep-cone value-tail — explicit, banked).** On `{|F|>μ²/2g}` the peak of `e^{−c|F|²/g²}e^{t(2g|F|−μ²)}` at `|F|=μ²/2g` gives `sup=e^{−cμ_j⁴/4g_j⁴}`. Along the tower `g_j→0`, so the deep cone carries mass `≤e^{−cμ⁴/g⁴}` uniformly; at the box edge `μ_j=2π/L` it is an `a`-fixed, summable factor.

**Step 5 (shallow cone — the only positive half).** Split `C_j={2g_j|F|>μ_j²}` into shallow `S={2g_j|F|≤c_Φ}` and deep. On `S`: `K^V=K^S+Φ″⪰−2g|F|+N/12≥0` (using `c_Φ=N/12`), so BL fires, `D≤0`, and `G″≥S`. So `G″_j(v)≥S(v)−D_deep`, and `D_deep≤e^{−cμ⁴/g⁴}·(…)`. **This is the entire positive content: shallow fiber ⇒ `D≤0`; deep ⇒ tail.**

**Step 6 (the conversion — FAILS HERE).** The needed inequality to force `G″≥δ` is
`(†)  D(v) ≤ S(v)−δ`, i.e. `Var_μ(f|0) ≤ E_μ[h_n−⟨K^{−1}c,c⟩|0]−δ`.
Step 4 bounds a **value** `∫e^{−V}`; `D` is the **variance** of the force `f=G_v+Φ′`, a second moment. Coercivity gives the wrong edge: `G″=E[H]−Var≤E[H]`, so `E[H]≥0` signs nothing, and `Var≥0` subtracts. The tail controls `P(|F|>μ²/2g)`, but on the band `N/12<2g|F|≤μ²/2g` (and the typical `|F|~g` at the box edge) the field is **not** tail-suppressed, `K` is indefinite, and the variance is `O(β_j)`, unbounded a priori. Equivalently `(†)` is
`(★)  Var_{μ_j}(Pᵀ∇S|0) ≤ E_{μ_j}[Pᵀ Hess S P|0]+W−δ` on `{2g_j|F|>μ_j²}`.
The three ingredients together deliver value (Step 4), geometric sign (Step 2), and shallow-fiber sign (Step 5) — but the cone variance **cannot** be lowered; **this is the exact missing step**.

**Step 7 (`k→0` — kills non-circularity).** At `μ_j=2π/L→0` the cone never empties (`2g|F|<μ²` fails for `|F|~g`), `Var(G)→0`, `G″→E[H]=m_eff²(k→0)`. So the needed `S≥δ` there **is** `m_eff²≥δ` = constraint (5). Using it to prove `G″≥δ` is circular.

**Step 8 (assembly — conditional).** *Were* `(†)` uniform, Steps 1–5 give `G″_j(v)≥δ_W` in `j,a,L` on `𝔗`; preservation (`g_{j+1}²<g_j²`, `μ_{j+1}=eμ_j`) keeps `𝔗` invariant; the multiscale transport closes the marginal log (`Σg_j⁴<∞`); and `lim_{k→0}G″≥δ` is the fixed-volume gap. Everything but `(†)` is banked.

## 5-point verdict
1. **Nonperturbative — PASSES as framework** (BBS trap, exact identity, no order-by-order); the theorem it yields is **conditional** on `(†)`.
2. **Controls SIGN — FAILS.** That *is* `(★)`: the lower bound on the cone variance/Hessian sign is unproven. Norm+tail+Jacobian pin the **value** (Step 4), the **β=0 sign** (Step 2), the **shallow sign** (Step 5) — never the cone variance.
3. **Non-circular — FAILS**, tied to (2): at `k→0` the `S≥δ` input is `m_eff²≥δ`=(5).
4. **Non-SUSY + fails for free field — PASSES.** SCBL gives `G″=0` **exactly** for `f₁²=0` (`S=D=0`); no SUSY used.
5. **`m_eff²(k→0)≥λ` uniform — FAILS.** Identical to `(★)=H1`.

**Exact missing step: Step 6, inequality `(†)=(★)`, the cone variance bound.** A value bound, a geometric sign, and a shallow-fiber sign do not compose into a cone-variance sign. Not closed — and the failure is named, not papered.

*Refs.* PROOF-TRAPPING-REGION §1–4; PROOF-MULTISCALE-4D §2,§4; PROOF-JACOBIAN-SAVES §1–2; PROOF-SAVVIDY-EXACT-SIGN §1,§3; PROOF-G1-BPHZ §3; INVENT-VARIANCE-INEQ §1–5; GAP-G1-largefield §1–3.
