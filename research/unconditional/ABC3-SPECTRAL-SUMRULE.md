# ABC3-SPECTRAL-SUMRULE — SHAPE A round 3: the spectral-moment self-consistency `x²−E[H]x+E[H]≥0` is **vacuous** at the measured per-mode `E[H]≈2`; it bites only above the threshold `E[H]≥4` (or `E[H]+W≥2√E[H]`), which no `N=2,3` measurement reaches

> 2026-11-23 · 判词 (verdict). Decided rigorously: **no.** A2's spectral rep `S=Σ_n|⟨0|(∇S)_v|n⟩|²/m_n²`, capped by `S≤E[H]/m₁²` (SD sum rule `Σ_n|c_n|²=E[H]`) and coupled to `G″=m₁²=E[H]−S`, does **not** give `x=m₁²≥λ` non-circularly at β=2.5 — it is **empty** there. The closure is `x²−E[H]x+E[H]≥0`, discriminant `Δ=E[H](E[H]−4)`. At the measured **per-mode `E[H]≈2.0`** (COMPUTE-GPP-SU3: 2.03/1.99 at SU3 β=5.7; W7 SU2: 2.46) `Δ=−4<0` ⇒ the quadratic holds for **every** `x`, pinning nothing. It bites only at `E[H]≥4` (with the Haar term, `E[H]+W≥2√E[H]`) — a large-effective-d.o.f. threshold no SU(2)/SU(3) datum reaches. The full Hessian cannot raise `E` (SD is *per unit direction*), higher moments and their ratios are **upper** bounds on `m₁²` (wrong sign), and UV/asymptotic-freedom pinning touches only total weight, not IR-dominated `S`. The unconditioned SD closure is exact but gap-blind. **1 nominal pass (vacuous for YM); fails sign, non-circularity, uniformity — one object: the IR weight = H1.**

## (1) Honest check first — the discriminant is negative (verified)
`x≥E(1−1/x) ⟺ x²−Ex+E≥0` (multiply by `x>0`). Roots `(E±√(E(E−4)))/2` exist only for `E≥4`. Measured: `E=2.0⇒Δ=−4.0`; `E=2.46⇒−3.79`; `E=3.18⇒−2.61`. All `<0`. The "lower bound" `x≥x_+` is therefore `x≥−∞` — satisfied by `x=0`, `x=1`, every value; it constrains the gap to nothing. Physically `G″>0` is real, but `E[H]−Var≈0` (equipartition) makes the bracket `E[H]−S≈0`; the positivity lives in `W`, not in this `x`.

## (2) Full Hessian does not raise `E`; higher moments are upper-signed
The SD identity `E[H_n]=E[f²_n]` is **per direction**. Summing over an orthonormal basis gives `E[tr H]=E‖f‖²` — **extensive** in `V` and colors — but the per-mode gap `G″(0)[v]=E[H_vv]−S_v+W_v` uses a *unit* `v`: `H_vv` is `O(1)≈2` whatever `N_c`, and `Σ_n|c_n|²=E[H_vv]`, not `N_c·E[H_vv]`. The trace is the wrong object; pairing it with one mode's `x` is `V`-dependent and meaningless — "use all modes/colors" cannot lift `E` above 4. Higher moments `S_k=Σ|c_n|²/m_n^{2k}` satisfy `m_n²≥m₁²⇒S_k≤S_{k-1}/m₁²⇒m₁²≤S_{k-1}/S_k` — an **upper** bound on `m₁²` (toy m₁=1,m₂=5: ratios 1.923→1.038→1.001, all `≥` true 1). A lower bound needs the opposite sign. Moments *pin the distribution* but never *lower-bound its infimum* — that is the gap.

## (3) UV finiteness / asymptotic freedom cannot constrain the low-q weight
`S_1=Σ|c_n|²/m_n²` is **IR-dominated** (`1/m_n²` up to `m_n→m₁`); the high-`m_n` tail contributes negligibly. AF pins `Σ|c_n|²=E[H]` and the UV tail — neither enters the IR weight that sets `S`. UV-weighted sums `Σ|c_n|²m_n^{2k}` are UV-dominated, diverge in the continuum, and do not feed `G″=E−S`. ABC2-DUALITY §1 verbatim: the sum rule fixes *total* weight, not the low-q distribution — a low-q bound must carry the gap.

## (4) Adding `W` — the exact missing step
Physical curvature is `G″=m₁²=E−S+W` (k0-CONSTRAINT-EFFPOT), so `x²−(E+W)x+E≥0`, `Δ=(E+W)²−4E`; threshold `E+W≥2√E`. At `E=2` this needs `W≥2√2−2=0.828`; the **proven** floor `W_geo≥N/12` gives only `0.25` (SU3), measured `0.32` — `2.6×` short, all still `Δ<0`. The missing `0.58` is the exact step, and it is the Haar term the *spectral rep `S` does not control*: `W` supplies the margin independently (a rename of `W_geo≥0`, not a moment-method gain).

## (5) Gap equation with the unconditioned SD identity
`Var(f|0)≤E[f²|0]` is free, but `E[f²]=E[H]` is **unconditioned**; conditioning breaks it (`E[f²|0]≠E[H|0]`), and that break *is* the sign-indefinite part = H1. Unconditioned it is equipartition `Var(f)=E[H]`, gap-blind. Free massless: `E[H]=0`, `S=∞` (massless δ at q=0) ⇒ cap `∞≤0` **unsatisfiable** ⇒ no `x`, `f₁²=0` — the acid test correctly rejects the free field.

## Five-point scorecard
1. **Nonperturbative — PASS.** Exact spectral rep + SD sum rule, all β.
2. **Sign-controlling lower bound — FAIL.** `Δ<0` at `E≈2` ⇒ vacuously true; moments upper-signed.
3. **Non-circular — FAIL at closing.** `S≤E/m₁²` *is* the planar gap (ABC2-DUALITY §1); conditioning = gap.
4. **Non-SUSY / fails free field — PASS (free field).** `E[H]=0,S=∞` ⇒ no `x` ⇒ λ=0; but the YM certificate is `x≥0` (vacuous).
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** Threshold `E≥4` / `E+W≥2√E` unreached; `1/N²` non-uniform.

**Sharpest line.** The verdict is the discriminant `E[H](E[H]−4)`: at `E[H]≈2` it is `−4`, so `x²−E[H]x+E[H]≥0` holds for every `x` and the "self-consistency" is the tautology `x≥0`. It bites only once `E[H]+W` reaches `2√E[H]` — needing `W≥0.83` where the proven floor is `0.25`. Below that, the closure is the cleanest algebraic statement of *why* the IR weight is H1: sign-blind, UV-blind, upper-signed in every moment — true for all `x`, hence about nothing.

*Refs.* ABC2-DUALITY (§1 spectral rep = planar gap, IR weight); ABC2-SYNTHESIS (§4 conditioning = gap; §5 scorecard); SUB1-COV-BOUND (SD unconditioned, gap-blind); INVENT-VARIANCE-INEQ (`G″=S−D`, `S=1/Var`); CONTINUUM-VARIANCE-DEFECT (denominator `S` is IR/H1); FORGE5-CONE-VARIANCE (`E[H]=E[f²]` equipartition); k0-CONSTRAINT-EFFPOT (`G″=Δ_dyn+W_geo+X`); COMPUTE-GPP-SU3 (per-mode `E[H]≈2.0`, `W_geo≈0.32`); W7-Gpp0 (SU2 `E[H]=2.46`); PROVE-Gpp-positive (`G″≈1/Var`); PROOF-JACOBIAN-SAVES (`W_geo≥N/12`); PROOF-IR-H1; PROOF-COARSE-POINCARE; Källén–Lehmann / Schwinger–Dyson; Brascamp–Lieb 1976; Prékopa 1971.
