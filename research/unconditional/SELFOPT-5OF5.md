# SELFOPT-5OF5 — the self-optimizing trajectory of the SL-multiscale tool: 4/5, and the exact final residual

> **Verdict (判词).** After optimization rounds over **every** named refinement of the backward step, **no construction reaches 5/5.** The multi-scale stochastic-localization (SL) tool is maximal at **4/5** (#1,2,3,4 PASS; #5 FAILS). The loop's *one genuine gain* is the **Doob–Meyer / quadratic-variation (QV) martingale**, which **removes** the `e^{Θ(N_□)}` volume-exponential the Holley–Stroock backward incurred (`O(1)` per slice, not a product over modes — Gaussian check `∫₀^∞dτ/(h+τ)²=h⁻¹`). But an equality is not a bound: it is *tight for Gaussian μ*, and closing it lands on the **exact single-mode inequality** `G″_∞≥λ`, which **is H1 = the mass gap**. The obstruction is relocated and shrunk, not removed.

## The tool, and the four genuine passes
Chart `X=log U`, `S=S_W+Φ` (`Φ=−log J`, `Φ″⪰N/12`). Multi-scale tilt `T=Σ_b t_bΠ_b`, `t_b=σ_b+ε` (`σ_b=` bare Savvidy deficit, not the gap). Exact `Hess(−log p_T)=Hess S+T⪰εI≻0` (**#2**, a value-free sign flip on every shell). Brascamp–Lieb on the log-concave `μ_T`: `Var_{μ_T}(a_v)≤⟨(Hess S+T)⁻¹v,v⟩`. Free field `σ_b≡0⇒t_b=ε→0`, `Var→1/ε→∞`, `λ=0`,`f₁²=0` (**#4** genuinely). No `g`/`1/N` expansion (**#1**); only bare curvature fed in (**#3** structurally).

## The internal loop — what each refinement achieved
1. **Per-scale cone-emptying** (SL-REFINED, SL-UV-CONE-EMPTY): *rate* improves, exponential survives — the cone is **per-configuration** but SL's tilt is a **fixed operator**, so `t_b=sup_xσ_b(x)` (a supremum is measure-blind; the cone's thinness cannot lower it), and β=2.5 negatives are **90–97% high-momentum**, `Σσ_bN_b^{cone}=Θ(N_□)` (`count×depth=∫|F|³` is volume-extensive).
2. **UV-convergence `Σosc_b`** (SL-SCALE-SUM): the sum is `O(1)` **iff the cutoff coupling empties the cone**, i.e. `g_J<0.26`, `β≳59` — *not* β=2.5, and even there only via **typical-field** scaling, which the sup (K1) does not honour; at β=2.5 `osc=Θ(N_□)`.
3. **W_geo-on-noncone** (SL-UV-CONE-EMPTY §b): `W_geo≥0` is a real drift but touches **no** `t_b`; it changes `G″`'s value, never `osc`.
4. **Doob–Meyer/QV martingale** (SL-MARTINGALE): `Var_μ(a)=∫₀^∞E[vᵀA_t²v]dt` — **additive, `O(1)`/slice, the `e^{Θ(N_□)}` genuinely dies**; the mass sits only in the un-convexified window `[0,t*]`, whose strict excess over Gaussian is the fibre-Savvidy defect `D(v)`.
5. **Self-consistent IR gap equation** (SL-IR-SINGLE): `m_eff²=M(t*)−t*` is the identity `m=m` in the Gaussian/RPA closure, or H1 in the exact one — **empty**; `A_t` is a martingale with `A_0=Var_μ(a)` **deterministic and equal to the unknown**, so stopping/isoperimetric arguments bound the decay *rate*, never `A_0`.

## The exact final obstruction
The QV martingale collapses the entire residual to one scalar inequality on the single `k→0` mode:
> **(BT′)**  `Var_μ(a_{k→0}) ≤ t*‖Cov_μ‖·Var_μ(a) + |v|²/ε`,  closing **iff** `t*‖Cov_μ‖<1`,
equivalently `G″_∞≥λ`, i.e. **H1**:
```
∃λ=λ(β,N)>0: ∀U′,∀coarse v,  E_μ[PᵀHP|U′][v,v] − Var_μ(⟨Pᵀ∇S,v⟩|U′) ≥ λ‖v‖²,  uniform in (L,a).
```
The fiber `H_ff` is indefinite (90–97% of Savvidy negatives lie there), so the Gaussian–Schur `S(v)` is sign-**undefined** and the defect `D(v)>0`; the `k→0` mode is itself concave (`PᵀHP` min ≈−0.3, `σ_∞≈0.3`), and its return ratio `Var_μ(a)/Var_{μ_T}(a)=1+t_∞/m_eff²` **is** `m_eff²`. Loosening the scalar to the compactness bound `Var≤π²` yields only the **false positive** a massless *compact* theory shares ([#4] fails off-`ℝ^V`), and on `ℝ^V` the free zero mode gives `∫₀dt/t=∞` (`λ=0`).

## Sharpest insight
The SL flow has **two independent walls**, and the round shows both terminate at H1: the **backward** cost `e^{Θ(N_□)}` (extensive, sup-forced because the tilt must log-concavify the *bulk fibre* Savvidy negativity at its worst configuration) and the **forward** energy `1/G″_∞`. The QV martingale *provably removes* the first (product→sum, `O(1)`/slice) — the literal task goal of #5's backward step is met — only to **expose the second exactly**: one magnitude `G″_∞≥λ` on the one `k→0` mode, whose martingale `A_0` *is* the answer. So #5 is not slack; in the uniform regime it is H1 **twice**. SL is **KLS-orthogonal**: it bounds the isoperimetric *decay rate*, never the absolute *initial datum* (KLS is scale-free; the gap is scale-covariant).

## Grade
**4/5** (#1,2,3,4 PASS, #5 FAILS at (BT′)=H1). The QV martingale is a real 4→4.5 advance — the exponential bookkeeping is gone — but the residual is the gap, so no honest 5/5. **Final obstruction = H1 = `Var_μ(a_{k→0})≤1/λ` uniformly in `(L,a)`,** reached from both the backward (`e^{Θ(N_□)}`) and, once optimized, the forward (`1/G″_∞`) side; closure needs a non-positivity, scale-covariant clustering bound that does not yet exist (REDTEAM-SPEC).

*Refs (internal).* HARDFORGE-SL-FULL (BT); FORGE5-SL-BACKWARD/REFINED/SCALE-SUM/UV-CONE-EMPTY/MARTINGALE/IR-SINGLE (the six refinements); INVENT-VARIANCE-INEQ (`G″=S−D`); PROOF-COARSE-POINCARE, PROOF-IR-H1 (H1=gap); PROOF-SAVVIDY-EXACT-SIGN; W6-BE-at-2.5; PROVE-Gpp; PROOF-JACOBIAN-SAVES; REDTEAM-SPEC. Eldan ECP 18 (2013) 63; Chen GAFA 31 (2021) 34; Brascamp–Lieb JFA 22 (1976) 366; Holley–Stroock Ann. Probab. 15 (1987) 1541; Doob h-transform; Itô isometry.
