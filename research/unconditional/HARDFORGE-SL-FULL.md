# HARDFORGE-SL-FULL — multi-scale stochastic localization: a complete attempted proof, and the exact step where it dies

> **Verdict (判词).** The construction closes four of the five hard constraints and **fails exactly one: the per-scale backward cost.** Stochastic localization's tilt is a *quadratic form* `½⟨x,Tx⟩` with arbitrary positive `T`; choosing `T=Σ_b t_b Π_b` (per-RG-shell tilts) is the requested multi-scale localization, and it makes `Hess S+T≻0` on *every* shell at once (Step 2) — each scale's relevant mode is convexified. Brascamp–Lieb on the tilted, log-concave `μ_T` then gives `Var_{μ_T}(a_v) ≤ ⟨(Hess S+T)^{-1}v,v⟩`, whose coarse-mode value at `T→0` is exactly `1/G″_∞` (Step 3). The wall is Step 4: un-tilting costs `e^{osc}`, `osc=osc(½⟨x,Tx⟩) ≲ ½π²Σ_b t_b N_b`, and at β=2.5 the Savvidy concavity is **bulk** — 22–49% of Hessian directions negative (W6-BE-at-2.5), i.e. `Θ(N_□)` concave modes — so `Σ_b t_b N_b = Θ(N_□)` and the transport is `e^{cN_□}=e^{ct/a²}`. Telescoping cannot help: the *count* of concave modes is volume-extensive, not the magnitude of `t`. The exact failing inequality is `(BT)`. The free-field acid test passes (`→λ=0`); nonperturbativity, sign control, non-circularity are intact. **Only constraint (5) fails, at Step 4 = H1.**

## Step 0 — Objects and goal
Flat chart `X=log U∈g^E`, `S=S_W+Φ` (`Φ=−log J` convex, `Φ″⪰N/12`, PROOF-JACOBIAN-SAVES). `μ∝e^{−S}dX`. Shells: `Π_b` projects onto `|k|∈[1/((b+1)a),1/(ba)]`; `N_b=`#modes, `Σ_b N_b=E\,dim G=Θ(N_□)`. **Goal (H1):** the coarse mode `v=v_∞` (shell `b=∞`, `k→0`) obeys `G″_∞:=E_μ[P_∞HP_∞|U′]−Var_μ(⟨P_∞∇S,v⟩|U′)≥λ>0`, uniform in `(L,a)`. `σ_b:=−λ_min(Hess S|_{shell b})` is the *bare* Savvidy deficit (not the gap); at β=2.5, `σ_b≤1.12β≈2.8` (W6-BE)/PROOF-SAVVIDY-EXACT-SIGN.

## Step 1 — Single-scale SL (exact)
Observe `Z=√t\,X+g`. Then `p_t∝exp(−S+⟨c_t,x⟩−½t|x|²)` and `Hess(−log p_t)=Hess S+tI` **exactly** (`λ↦λ+t` — the one flow that convexifies, PROOF-STOCHASTIC-LOCALIZATION §1; contrast `λ↦λe^{2tλ}`, PROOF-FLOW-CONVEXITY-TRANSFER). Log-concavity once `t>max_b σ_b≈2.8`. BL: `Var_{μ_t}(a_v)≤⟨(Hess S+tI)^{-1}v,v⟩_{μ_t}`.

## Step 2 — Multi-scale tilt (the construction)
Standard SL uses `T=tI` and convexifies only at `t>2.8` — the wrong end. Declare `T:=Σ_b t_bΠ_b`, `t_b:=σ_b+ε`. Then
`Hess S+T` has eigenvalues `λ_b+t_b=λ_b+σ_b+ε≥ε`, so `Hess S+T⪰εI≻0`: `μ_T∝exp(−S−½⟨x,Tx⟩+⟨c,x⟩)` is log-concave **and each shell is lifted by its own `t_b`**. Non-SUSY; a genuine probability tilt.

## Step 3 — BL energy on the coarse mode
`Var_{μ_T}(a_v)≤⟨(Hess S+T)^{-1}v,v⟩=Σ_b(λ_b+t_b)^{-1}|v_b|²`. On `v=v_∞`, reducing the fiber, the `T→0` right side is `S(v)−D(v)=G″_∞`, reciprocal `1/G″_∞` (PROOF-STOCHASTIC-LOCALIZATION §2). The SL energy *is* the coarse reciprocal — an identity, saturated for Gaussian `μ_T`.

## Step 4 — the backward transfer: WHERE IT DIES
Since `dμ=w\,dμ_T` with `w=e^{−ψ}`, `ψ=½⟨x,Tx⟩−⟨c,x⟩`, Holley–Stroock gives `Var_μ(a_v)≤e^{osc}Var_{μ_T}(a_v)`, with `osc ≤ ½π²Σ_b t_b N_b` (compactness `|x_e|≤π`). At β=2.5 the Savvidy-negative directions are **bulk**: 22–49% of directions (W6-BE-at-2.5), 90–97% of them in the fiber (PROOF-COARSE-POINCARE (a)(iii)), so `Σ_b t_b N_b = Σ_b(σ_b+ε)N_b = Θ(N_□)`, with `osc=½π²Θ(N_□)=ct/a²`, `c\gtrsim0.2`. **Telescoping fails two ways.** (a) A fine schedule is no cheaper than one jump: the per-mode Gaussian charge is `∫_0^{t_b}½\,dτ/(λ_b+τ)`, and for a Savvidy mode `λ_b=−σ_b` this integrand has a pole at `τ=σ_b∈(0,t_b]`, so each concave mode contributes an `O(1)` (compactness-regulated) charge whatever the step size. (b) The *number* of such modes is `Θ(N_□)`, so the total is `e^{Θ(N_□)}` regardless of `t_b` size. Hence
> **(BT)**  `Var_μ(a_v) ≤ e^{cN_□}\big/(G″_∞+t_∞)`, and `Var_μ(a_v)≤1/λ` requires `G″_∞ ≥ λ\,e^{cN_□}−t_∞ \gg λ`.

The backward cost does not merely degrade the bound — it demands a *strictly stronger* positivity than H1.

## Step 5 — Assembly and the three roads
Road (i) full-measure transfer: obstructed by (BT). Road (ii) convexify the coarse marginal alone: the object to lift is `G″_∞`, whose semi-concavity *is* `−G″_∞` — circular (assumes `G″≥λ`). Road (iii) P3 anchor (β<1/48: `T=0` free, `λ_μ>0`, SZZ; PROOF-P3-UNCONDITIONAL) transported to 2.5 by `Z_O` (INVENT-TRANSPORT-INEQ): only through `(CP)`, the single-phase/continuum wall = H1-continuum. All three terminate at H1.

## Free-field acid test (mandatory)
Massless free field: `Hess S=A⪰0` with a `k=0` zero mode; `σ_b≡0`, so `t_b=ε→0` may be taken *freely*, `osc=½ε{π²N_□}→0`. Then `Var_{μ_T}(zero mode)=1/ε→∞`, giving `λ=0`, `G″=0`, `f₁²=0`. **PASS** — no fabricated positivity, and the vanishing `osc` (vs YM's `osc=Θ(N_□)`) is precisely what discriminates the two: the cost is *caused by* the Savvidy concavity.

## Five-constraint verdict (honest)
1. **Nonperturbative — PASS.** Exact `Hess+tI` identity, BL, Holley–Stroock; no `g`/`1/N` expansion.
2. **Sign control (not value) — PASS on the tilted cone.** `Hess S+T⪰ε` is a true, value-free sign statement on the Savvidy cone; it is sign-definite and requires no magnitude estimate.
3. **Non-circular — PASS (structure).** `σ_b` is *bare* curvature, not the gap; no gap-as-input. The *closing* step is nonetheless circular (Road ii) or H1 (Road i).
4. **Non-SUSY + fails for free field — PASS.** Probabilistic tilt, not localization-cohomology (contrast H1-T7); the free field returns `λ=0` (acid test above).
5. **Uniform lower bound `m_eff²(k→0)≥λ` — FAIL.** The only non-circular output is `≥λe^{−cN_□}→0` (Road i) or circular (Road ii). **Exact failing step = Step 4: the per-scale backward cost `e^{osc}=e^{Θ(N_□)}` — i.e. H1.**

*Refs.* R. Eldan, ECP **18** (2013) no. 63; Chen, GAFA **31** (2021) 34; Klartag, Mem. AMS **249** (2017); Klartag–Lehec (2022); Brascamp–Lieb, JFA **22** (1976) 366; Holley–Stroock, Ann. Probab. **15** (1987) 1541. Internal: PROOF-STOCHASTIC-LOCALIZATION (`Hess+tI`, energy `=G″`, `1/G″`); PROOF-FLOW-CONVEXITY-TRANSFER (`λe^{2tλ}`); INVENT-CONVEXIFICATION ((a),(b) transfer `b/a=e^{osc}`); PROOF-JACOBIAN-SAVES (`Φ″⪰N/12`); PROOF-COARSE-POINCARE (H1, fiber-indefinite); PROOF-P3-UNCONDITIONAL (`r₀`, β<1/48); INVENT-TRANSPORT-INEQ (`(CP)`); W6-BE-at-2.5 (`−1.12β`, 22–49% neg); PROOF-SAVVIDY-EXACT-SIGN; H1-T7-localization (localization ≠ SL).
