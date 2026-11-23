# FORGE5-SL-BACKWARD — multi-scale / terminal-slice stochastic localization against the `t→0` gap

> **2026-11-23 · Fifth-hard-constraint forge: multi-scale SL or a direct `t→0` slice without the exponential backward cost?** **Verdict: cannot satisfy all 5.** SL genuinely flips the *tilted* sign (`Hess S_t = Hess S + tI`) and literally passes the free-field acid test — but the *mechanism* it uses to escape the backward cost is **compactness of the target, not the interaction**, and the cost it cannot escape is **volume-exponential**, not the UV `e^{ct/a²}` (a red herring at `k→0`). **Passes (1); fails (2), (3), (5); passes (4) in letter only.** The exact failure is a cheap-or-vacuous dichotomy on the `k→0` mode driven by its range `√V`.

## (1) Nonperturbative — PASS
No loop expansion; SL/LSI and Brascamp–Lieb are exact. Fine.

## The multi-scale attempt, and why `t` cannot be per-scale
The tilt `(t/2)|x|² = Σ_k (t/2)|x̂_k|²` shifts **every** Fourier eigenvalue by exactly `+t`. There is no "`t`-per-scale": a frequency-selective tilt `⟨x,Tx⟩`, `T=Σ_j t_jP_j`, either localizes **UV only** — `k→0` untouched, gap unchanged — or **IR only**, which is the single-mode problem below at the same cost. So "multi-scale SL" is either *useless* or *identical* to one scale. The convex window `t ≳ 1.12β ≈ 2.8` is an **absolute** `t` set by the Savvidy magnitude, not a scale.

## (3) Non-circular — FAIL: the backward ratio *is* the gap
For the `k→0` mode `a` (Gaussian-exact, verified numerically), the tilt-and-return variance ratio is `Var_μ(a)/Var_{μ_{t*}}(a) = 1 + t*·Var_μ(a) = 1 + t*/m_eff²`: `≈10` for gapped YM (`t*=2.8`, `m_eff²≈0.3`) but `→10⁶` for a near-massless mode. **Certifying "ratio `=O(1)`" is certifying `Var_μ(a)=O(1)` — the gap itself.** The only *a-priori* bound is Holley–Stroock, `Var_μ(a) ≤ e^{osc(p_0/p_{t*})}Var_{μ_{t*}}(a)`.

## (5) Uniform lower bound — FAIL: the a-priori backward factor is `e^{cV}`
`p_0/p_{t*} = e^{+(t*/2)|x|² - t*⟨x,Y⟩}` depends on `|x|² = Σ_ℓ|X_ℓ|²` **globally**; over `X ∈ torus^V` its oscillation is `(t*/2)·(2π)²N·V ∼ t*V`. So the certified ratio is `e^{t*V}` — **volume-exponential** (already `>e^{200}` at `V=16`) — and cannot be localized to the target block, since the ratio is a product over *all* links. The dichotomy is exact: the **block** mode has range `O(1)` (cheap backward, `osc~t*R²=O(1)`) but does **not probe `k→0`** — its "gap" is trivial boundedness a massless *compact* theory shares (false positive); the **`k→0`** mode has range `√V` (`osc~t*V`) and is where the gap lives. No scale is both cheap and gap-relevant. Hence no uniform `λ`.

## (2) Sign at the cone — FAIL; (4) acid test — pass in letter only
SL flips `Hess S` on the Savvidy cone pointwise for `t>1.12β`, but that certifies the **tilted** sign; the measured un-tilted `G″=+0.5` is the `+W` SL's `+t` replaces. At `t=0` there is no control — the indefinite fiber — so the cone sign is unproven. The acid test passes (massless `ℝ` free field: `a` unbounded ⟹ `osc=∞` ⟹ `λ=0`) but only because the free target is **non-compact**; the positivity is compactness-driven, not `f₁²`-driven.

## Sharpest line
The gradient flow leaves the `k→0` mode *invariant* (multiplier `→1`); SL *shifts* it (`m_eff²↦m_eff²+t`). SL therefore genuinely *changes* the gap it must prove and must **un-shift** it — and the un-shift is exactly the gap (`ratio=1+t/m_eff²`), whose only a-priori bound is `e^{t*V}`. Both flows die on the `k→0` mode, from opposite sides: the return undoes the `+tI` that did the work. **Constraints 2, 3, 5 fail; 1 passes; 4 passes only by non-compactness of the free field, not by the interaction.**

*Refs.* Eldan, ECP **18** (2013) 63; Chen, GAFA **31** (2021) 34; Klartag, Mem. AMS **249** (2017); Brascamp–Lieb 1976; Holley–Stroock, Ann. Probab. **15** (1987). Internal: PROOF-STOCHASTIC-LOCALIZATION (same wall, terminal slice); INVENT-CONVEXIFICATION (periodic no-go, `e^{2βN_□}` ratio); PROOF-FLOW-CONVEXITY-TRANSFER (`λ↦λe^{2tλ}`, `k→0` invariant); PROOF-COARSE-POINCARE (H1=`m_eff²`); PROOF-JACOBIAN-SAVES (`W_geo≥0`, `+0.5` vs Savvidy `−1.12β`); PROOF-SAVVIDY-EXACT-SIGN (`−g|B|`); W6-BE-at-2.5.
