# FORGE5-SL-MARTINGALE — the Doob–Meyer / importance-sampling backward for the coarse mode: the volume-exponential dies, the residual is exactly `D(v)`

> **2026-11-23 · Direct attack on constraint #5.** Replace the Holley–Stroock backward reweight `e^{(t*/2)|x|²−t*⟨x,Y⟩}` (cost `e^{Θ(N_□)}`, FORGE5-SL-BACKWARD §5) by the **exact variance-transport (Doob–Meyer quadratic-variation) identity** of the SL martingale, which transmits `Var_μ(a)` *additively* — one scalar per time slice, no product over modes — so **the literal task goal is met: the per-mode cost is `O(1)` and the volume-exponential is genuinely removed** (verified against the Gaussian `1/h`, no `e^{t*V}`). But an equality carries no bound. Splitting at the coarse convexification threshold `t*=σ_∞+ε`: the log-concave tail `t≥t*` is volume-free `O(1/ε)`, while the **mass sits in the non-convex window `[0,t*]`**, where `vᵀA_t²v=Var_{μ_t}(a)²+|A_ty|²` and the identity is **saturated for Gaussian μ**. The only strict excess it admits is the non-Gaussian fibre coupling, whose stationary form is `D(v)` (INVENT-VARIANCE `G″=S−D`) — i.e. **the same H1**. Important sampling pays `χ²=e^{Θ(N_□)}`; a Doob `h`-transform targeting the gap core is circular (its Radon–Nikodym derivative *integrates the Savvidy volume*). **Verdict: #5 fails again, but the obstruction has moved and shrunk** — from a volume-exponential `e^{Θ(N_□)}` to one un-convexified-window contraction whose residual is exactly the fibre-Savvidy excess `D(v)`. 1 pass; 2,3 fail at closing; 4 letter-only; 5 = H1.

## 1. The martingale and the quadratic-variation identity (exact, no reweight)
Continuous SL `dZ_t=Xdt+dB_t`, `μ_t=law(X|Z_t)`, `p_t∝exp(−S+⟨Z_t,x⟩−(t/2)|x|²)`, `Hess(−log p_t)=Hess S+tI` (PROOF-STOCHASTIC-LOCALIZATION §1). Innovation `dW_t=dZ_t−E_{μ_t}[X]dt` is Brownian; the **conditional mean** `M_t^f:=E_{μ_t}[f]` is a martingale (`M_0=E_μf`, `M_∞=f(X)`) with `dM_t^f=Cov_{μ_t}(f,X)dW_t`. For `f=a=⟨v,·⟩`, `Cov_{μ_t}(a,X)=A_t v` (`A_t:=Cov_{μ_t}(X)` — a *supermartingale*, not the martingale), so Itô isometry gives
> **(QV)**  `Var_μ(a)=E(M_∞−M_0)²=∫₀^∞ E|A_t v|² dt=∫₀^∞ E[vᵀA_t²v]dt.`

No `dμ/dμ_T`, no `e^{osc}`: additive, `O(1)` per slice. Gaussian check `A_t=(H+tI)^{−1}`: `∫₀^∞vᵀ(H+t)^{−2}v\,dt=vᵀH^{−1}v=Var(a)`, mass at `t∈[0,h]`, `h=G″_∞`.

## 2. Tail free; window carries the mass
The coarse shell is lifted at `t*=σ_∞+ε`; for `t≥t*` BL on the log-concave `μ_t` gives `vᵀA_t²v≤|v|²/(t−σ_∞)²`, so `I_>:=\int_{t*}^\infty...\,≤|v|²/ε` (volume-free). By the Gaussian check the *whole* variance lives at `t∈[0,h]⊂[0,t*]` (since `h=G″_∞≪t*`); `I_<:=\int₀^{t*}...` carries it. **The mass sits precisely where log-concavity fails.**

## 3. Window decomposition — where `D` re-enters
Split `v`/fibre `y`: `vᵀA_t²v=A_{vv}²+|A_{vy}|²`, `A_{vv}=Var_{μ_t}(a)`, `A_{vy}=Cov_{μ_t}(a,y)`. The first term is the target's own variance (circular); the second is the fibre coupling — for Gaussian it supplies the Schur completion and (QV) is tight, so the **strict excess over Gaussian is the non-Gaussian fibre response**, measured `≈2%` (PROVE-Gpp) and equal to the BL defect `D(v)` (INVENT-VARIANCE), positive iff the fibre is Savvidy-indefinite. Monotonicity `A_t⪯A_0` gives `vᵀA_t²v≤‖A_t‖·vᵀA_t v`, hence the **scalar** contraction
> `I_< ≤ t*·‖Cov_μ‖·Var_μ(a).`
This is the real gain: the `e^{Θ(N_□)}` **product over modes collapses to one scalar** `‖Cov_μ‖` (max fibre variance), because (QV) never multiplies modes.

## 4. Exact failure
(QV) is an equality, *tight for Gaussian μ*; YM's coarse marginal is Gaussian to `G″≈1/Var` within `2%`, so any bound from (QV) alone is `≈98%` the Gaussian/Schur answer and cannot beat it. The strict content is that `2%` non-Gaussian excess `D(v)>0`. With §2–3:
> **(BT′)**  `Var_μ(a) ≤ t*‖Cov_μ‖·Var_μ(a) + |v|²/ε`, closing only if `‖Cov_μ‖<1/t*`, i.e. fibre curvature `≳t*=σ_∞+ε`.

But the target is `λ≪t*`: `‖Cov_μ‖≈1/λ` is the object H1 must bound (`G″≥λ`). On the compact torus the sup cutoff even yields the false positive `Var≤π²`. Equivalently, the marginal covariance must obey the Fisher/Riccati inequality `g_t'≤−c g_t²` (`c`=λ-scale) on `[0,t*]`, i.e. its LSI constant `≥λ` — `G″_∞≥λ`. **Importance sampling** (`dμ/dμ_T=e^{−ψ}`) reproduces `χ²(μ‖μ_T)=e^{Θ(N_□)}−1`; the **Doob `h`-transform** setting `Var(a)=σ²` imposes the answer by fiat (its `dμ/dh` integrates the Savvidy volume). None evades (BT′).

## 5. Five-point verdict (honest)
1. **Nonperturbative — PASS.** (QV) exact, Itô isometry; no expansion, no reweight.
2. **Sign control — FAIL.** (QV) is sign-neutral; positivity needs `∫₀^{t*}G″` positive — no value-free sign isolated from H1.
3. **Non-circular — FAIL at closing.** The window integral's Gaussian part *is* the target; the strict part is `D(v)`; (BT′) assumes `G″≥λ`.
4. **Non-SUSY / free field — PASS only off-compact.** Massless `ℝ^V`: `∫₀ dt/t=∞` at the free zero mode, `λ=0` ✓; but on the compact torus compactness gives `Var≤π²`, a false positive a massless *compact* theory shares (same caveat as FORGE5-SL-BACKWARD).
5. **Uniform `Var_μ(a)≤1/λ` — FAIL.** Exact failure `(BT′)` ⟺ `D(v)<S(v)` = `G″≥λ` = H1.

**Sharpest line.** The quadratic-variation identity *is* the right un-tilting — it costs a sum, not a product, so the `e^{Θ(N_□)}` volume-exponential genuinely dies and "`O(1)` per mode" is literally achieved. But an equality is not a bound: it is tight for Gaussian μ, and the only strict gap it admits is the non-Gaussian excess of the one un-convexified window, which is `D(v)` — the same fibre-Savvidy object the Holley–Stroock backward, INVENT-VARIANCE and the RG transport already isolated. The martingale buys **bookkeeping** (product → one scalar `‖Cov_μ‖`), not closure: **#5 = H1, relocated to `[0,t*]` and shrunk to the single contraction `t*‖Cov_μ‖<1`.** Doob-`h` to the core is circular; IS pays `χ²=e^{cN}`; no re-arrangement of the flow beats the saturation.

*Refs.* Eldan, ECP **18** (2013) 63; Chen, GAFA **31** (2021) 34; Klartag, Mem. AMS **249** (2017); Brascamp–Lieb, JFA **22** (1976) 366; Holley–Stroock, Ann. Probab. **15** (1987) 1541; Doob (h-transform); Efron–Stein (1981). Internal: HARDFORGE-SL-FULL (Step 4 = H1, (BT)); FORGE5-SL-BACKWARD (`e^{t*V}`, acid-test compactness caveat); FORGE5-SL-IR-SINGLE (`A_0`=unknown; here corrected — the martingale is `M_t=E_{μ_t}[a]`, (QV) is its exact transport); PROOF-STOCHASTIC-LOCALIZATION (`A_t=(Hess S+tI)^{-1}`, free `∫dt/t`); INVENT-VARIANCE-INEQ (`G″=S−D`, `D`=fibre BL defect, `D/S≈2%`); PROOF-COARSE-POINCARE (H1=`m_eff²`, fibre `H_ff` indefinite); PROOF-TRANSPORT-RG (transport = signed magnitude, no lower bound); PROOF-JACOBIAN-SAVES (`G″≈1/Var`); PROVE-Gpp.
