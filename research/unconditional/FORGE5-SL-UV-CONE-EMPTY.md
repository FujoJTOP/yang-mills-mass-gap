# FORGE5-SL-UV-CONE-EMPTY — split the SL backward cost at the cone boundary: the *sup* kills the split, and the β=2.5 cone is high-momentum, not UV-empty

> 2026-11-23 · FORGE5 attack on constraint #5, the `e^{cN_□}` backward of HARDFORGE-SL Step 4. Proposed split: `osc=½π²Σ_b t_bN_b` comes only from the concave (cone) region — the UV trap empties the cone (a), and outside it `Hess S⪰0` so `t_b=ε` with `W_geo≥0` supplying the drift (b). **Verdict: 4/5; #5 fails, and the split is the wrong shape to fix it.** **(K1 — fatal):** the tilt is a *single global operator* log-concavifying the *whole* support, so `t_b=sup_xσ_b(x)`; a supremum is measure-blind — the cone's *thinness* cannot lower it. **(K2 — empirical):** at β=2.5 the measured negatives are **high-momentum** (90–97% fluctuation block, W6-BE §3) — the cone *lives in the UV*, so the trap's emptying (the deep-UV `g→0` statement) is *inverted* at the physical coupling; and the `k→0` goal mode is itself concave (`PᵀHP` min ≈−0.3). Hence `osc=Θ(N_□)` at β=2.5.

## K1 — the split is measure-blind because the transfer is a sup
`supp μ_T=supp μ` = the whole torus, so Brascamp–Lieb needs `Hess S(x)+T⪰0` at **every** `x`, and Holley–Stroock charges `osc=ess sup ψ−ess inf ψ` over the same set. Hence `t_b≥σ_b^{max}:=sup_x(−λ_min(Hess S(x)|_{Π_b}))`. The cone `{x: λ_min<0}` is a set of *configurations*; the tilt is one operator indexed by *shells*. No `T` is `ε` "outside the cone" and large "inside": `T` must dominate the worst case — the deeply concave large-field tail, positive (if tiny) measure. **Thin ≠ absent ≠ free** — for a supremum, even `e^{−c/a⁴}`-rare suffices. The `e^{−S}` weight suppresses the cone's *contribution to `G″`*, never its *cost in `osc`*.

## K2 — at β=2.5 the cone is the high-momentum block, the opposite of "UV-empty"
Shell `b` is concave for some config ⇔ `k_b²<2g\,F_max`; only `k_b²>2gF_max` is cone-free (`t_b=ε`). At β=2.5 (`g²=2N/β≈2.4`, `g≈1.55`) `k_*:=√(2gF_max)` reaches the cutoff: the measured negatives are 90–95% *in the fluctuation (UV) block*, only 3–10% on the low-momentum slice — evidence `k_*≳k_max`, so essentially **no shell is cone-free**; `osc=Θ(N_□)` verbatim (the 22–49%). The trap's emptying (`2g_j|F|<μ_j²`) is the `g_j→0` AF statement of the *UV tower* — exactly where the gap is not; at `g≈1.5` it does not apply.

## The `k→0` goal mode is *in* the cone
The blocked coarse mode — the H1 mode — has `PᵀHP` min ≈−0.3 at β=2.5 (a slow color-adjoint fluctuation sees `−2g|F|` with no `k²` to rescue it). The split **cannot exclude the goal mode**: it sits on the expensive side.

## (b) `W_geo≥0 ⇒ G″≥S` neither holds pointwise nor lowers `osc`
`W_geo≥0` is a true floor (`≥N/12`), but `G″=E[H]−Var(G)+W_geo+X` is **not** `≥S`: `Var(G)~0.2` and sign-indefinite `X` are not absorbed by an `N/12` floor, and `W_geo` is anti-correlated with `|F|` (`ρ≳3–16`). Decisively, `W_geo` never touches `t_b`: it is a *drift*, not a sign-of-`Hess S` statement — it changes `G″`'s value, never `osc`.

## Five-point verdict
1. **Nonperturbative — PASS.** Exact cone split, Savvidy spectrum, HS/BL; no loop expansion.
2. **Sign control — PASS (on `μ_T`).** `Hess S+T⪰0` value-free; but `σ_b^{max}` is bare curvature, not the gap.
3. **Non-circular — PASS (structure), still vacuous.** `σ_b`, `k_*` bare; residual `Var≤e^{osc}/G″_∞` again demands `G″_∞≥λe^{osc}`.
4. **Non-SUSY + fails free field — PASS.** Free field `F≡0⇒k_*=0⇒σ_b≡0⇒osc=0` (cheap backward) yet `λ=0` from the zero mode: cheap backward is necessary, never sufficient; the *nonempty* cone is what distinguishes YM.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** `osc=Θ(N_□)` at β=2.5; no `λ`.

## Exact failure
`osc=½π²Σ_{b:\,k_b<k_*}(σ_b^{max}+ε)N_b`, with `σ_b^{max}=sup_x`-**forced** by log-concavity on all of `supp μ` (K1) and `k_*≳k_max` at β=2.5 (K2). No rescue from the cone's measure or from `W_geo`: a supremum ignores measure, a drift ignores Hessian sign. And even granting `O(1)` backward, Step 3's forward energy is `1/G″_∞`, so #5 would only *relax* to `G″_∞≥λ` (**= H1 exactly**), not close. The split neither reaches `O(1)` backward at β=2.5 nor, if it did, would it be 5/5.

**Sharpest line.** The split's *premise* is right — only concave modes cost — but the cost is a **supremum**, which sees the cone at full strength however thin it is; and the β=2.5 cone is the *high-momentum* block (where W6 finds the negatives), so the intended saving (UV emptying) acts precisely in the AF UV where the gap is *not*, and is void at `g≈1.5` where it is. The `k→0` goal mode sits on the expensive side. The exponential is the **count of high-momentum + coarse cone modes** — the Savvidy phase-space volume, i.e. the gap in mode-counting clothes.

*Refs.* HARDFORGE-SL-FULL (Step 4, (BT)); PROOF-TRAPPING-REGION (§2 cone-empty on `𝔗`, `g_j→0`); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`, `N_neg∝∫|F|²`); W6-BE-at-2.5 (§3 negatives high-momentum, 90–95% fluctuation block, 22–49%); PROOF-COARSE-POINCARE (a)(i) `PᵀHP` min `−0.3` at `k→0`; PROOF-SAVVIDY-WGEO-DOMINATION (`ρ≳3–16`, `W_geo` anti-correlated); PROOF-JACOBIAN-SAVES (`W_geo≥N/12`, `G″=Δ+W+X`); FORGE5-SL-BACKWARD (compactness backward); FORGE5-CONE-VARIANCE; Holley–Stroock, Ann. Probab. **15** (1987) 1541; Brascamp–Lieb, JFA **22** (1976) 366.
