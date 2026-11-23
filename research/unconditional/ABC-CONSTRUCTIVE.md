# ABC-CONSTRUCTIVE — SHAPE B: the coercivity norm, value-tail, `W_geo`, and exact cone count share ONE exponential; it purchases the AF-UV cone and nothing where the gap lives; the sign is (★)=H1

> **Verdict (判词).** A completed nonperturbative 4D YM construction (BBS trap + AF log-resummation + cone-emptying + the exact identity `G″=E[H]−Var+W`) is non-circular on its bulk, but making its **last** step — the Savvidy-cone Hessian — sign-controlling via the **combined structure** does *not* close the marginal sign. The combination is not four handles: on the cone they share **one exponential**, `Q=∫TrF² ≥ Vμ_j⁴/4g_j²`, which a supremum and conditioning divide out. It genuinely signs the shallow band and the AF-UV deep cone; it fails at β=2.5, where the cone is the *high-momentum* block and the `k→0` mode sits inside it. The required inequality there is exactly **(★)=H1**. **Honest output: exact failing estimate, not a completion.**

**§1 The combined object is one exponential.** Trap `𝔗={0<g_j≤g*, ‖R_j‖_X≤rg_j², 2g_j|F|>μ_j²}` with jet norm `‖R‖_X=e^{Q/4g_j²}|δ²R|ℓ` (TRAPPING §1). On the cone `|F|>μ_j²/2g_j`, so `Q≥Vμ_j⁴/4g_j²`. Hence the **coercivity multiplier** `e^{Q/4g_j²}` and the **large-field value-tail** `e^{−cμ_j⁴/g_j⁴}` are the *same* exponential — one applied to the remainder's second variation, one to the cone's probability. Holding `‖R‖≤rg²` forces `|δ²R|≤rg²e^{−cμ_j⁴/g_j⁴}`: the large field makes the remainder *exponentially smaller*, never larger — it cannot reach the cone depth `−2g_j|F|≤−μ_j²`. Coercivity bounds the action's *value*, not its *Hessian* (G1-fixed-volume (c)).

**§2 What the combination genuinely closes.** Split `C_j`.
- **Shallow** `{2g_j|F|≤c_Φ=N/12}`: `Φ″⪰N/12` ⇒ `K^V=K^S+Φ″⪰−2g|F|+N/12≥0`, BL fires, `D≤0`, `G″≥S`. Closed (HARDFORGE Step 5).
- **AF-UV deep cone** (`g_j→0`, `μ_j=Λe^j`): destabilizing contribution `≤ depth×measure = 2g|F|·e^{−cμ_j⁴/g_j⁴} ≤ μ_j²e^{−cμ_j⁴/g_j⁴}→0`. Closed.

Both are real advances over a pure norm. Both live in the asymptotically-free UV.

**§3 Where it fails — a supremum is measure-blind.** At β=2.5 (`g≈1.5`) the measured negatives are **90–97% in the high-momentum fluctuation block**; the cone is *not* UV-thin, and the `k→0` goal mode is itself concave (`PᵀHP` min≈−0.3). The trap's emptying `2g_j|F|<μ_j²` is the `g_j→0` statement — exactly where the gap is *not*. The trapping norm is a **supremum**: thinness cannot lower it (even `e^{−c/a⁴}`-rare configs are in the sup), so the coercivity weight suppresses the cone's *contribution to `G″`*, never its *cost*. On the returned measure the required inequality is
```
(★)  Var_μ(⟨v,∇S⟩|0) ≤ E_μ[H_vv|0]+W−δ ⟺ Var(a)≤1/δ ⟺ m_eff²(k→0)≥δ = H1,
```
and no exponential is available (`μ~g` ⇒ weight `O(1)`). `W_geo≥0` is a **drift, not a Hessian sign** — anti-correlated with `|F|` (ρ≳3–16), and an `N/12` floor cannot absorb `Var(G)~0.2` plus sign-indefinite `X`. The **exact count/depth** makes the residue terminal, not tameable: `N_neg∝(gB)²`, `λ_min=−g|B|`, `Tr(H_-)≤8C_N∫|F|³`, all monotone increasing in `|F|`; the only true positivity, `Tr H=4Tr(−D²)≥0`, yields `Tr(H_-)≤Tr(H_+)` and saturates at `0.94`. **No positivity is extractable.**

**§4 Five-point verdict.**
1. **Nonperturbative — PASS (framework).** Trap + AF resummation `Σ_j g_j⁴=ζ(4)/(16β₀⁴)<∞`; the sign is BPHZ-invisible.
2. **Controls SIGN — FAIL.** Shallow band + AF-UV deep cone are signed (new); the β=2.5 high-momentum cone and the `k→0` mode are indefinite — that is (★).
3. **Non-circular — PASS on the bulk** ((★) strictly stronger than H1, `⊥` gap); **FAIL at the terminal `k→0` step = H1.** As a path to the gap: a rename, zero slack.
4. **Non-SUSY + fails free field — PASS.** `W_geo≥0` (nonabelian `Φ″≥N/12`) is the only provable positivity and it vanishes for the Gaussian free field (`F≡0 ⇒ σ≡0`, `λ=0`).
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** Exactly H1.

**Sharpest line.** Coercivity + value-tail + Jacobian + exact count are not four handles: on the cone they share **one** exponential `Q≥Vμ⁴/4g²`, which a supremum and conditioning divide out. The combination buys the shallow band and the whole AF-UV deep cone — and is void at β=2.5, where the cone is high-momentum, the goal mode is *inside* it, and the weight is `O(1)`. The marginal sign closes everywhere the cone is exponentially improbable and nowhere it is probable; the gap lives exactly where it is probable. **Not a completion — the exact failing estimate is (★)=H1 on the physical-coupling cone.**

*Refs.* PROOF-TRAPPING-REGION §1–4; PROOF-MULTISCALE-4D §2,§4; FORGE5-CONSTRUCTIVE; PROOF-SAVVIDY-EXACT-SIGN §1–3; PROOF-JACOBIAN-SAVES §1–2; HARDFORGE-CONSTRUCTIVE-FULL §3–6; FORGE5-CONE-VARIANCE; FORGE5-SL-UV-CONE-EMPTY (K1 sup measure-blind, K2 high-momentum cone, `ρ≳3–16`); PROVE-Gpp-positive §1,§4–5; GAP-G1-largefield §1,§3; G1-fixed-volume (c); PROOF-IR-H1; Bauerschmidt–Brydges–Slade LNM 2242 (2019); Balaban CMP 109/116/122; Savvidy PLB 71 (1977) 133.
