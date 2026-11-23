# G1-L3-TOWER-POINCARE — is the block Poincaré `Λ_UV` uniform over the `J∼log(1/a)` tower, or does it degrade per shell?

> **判词 (verdict).** It **degrades, and the degradation *is* the obstruction** — with one refinement that decides the file. The **only proven** per-shell constant is the block-local Holley–Stroock gap `Λ_UV,j=e^{−2β_jN_p}/C_G`. Since `β_j=2N/g_j²∼4Nβ₀\,j` is **monotone increasing into the UV**, this constant falls shell-by-shell; the tower minimum `inf_jΛ_UV,j` sits at the **deepest** shell and equals `(aΛ)^{c}`, `c=2·(4Nβ₀N_p)/ln2>0` — a **positive power of `a`**, not `a`-uniform. This is **not** the coupling scale-sum (`Σ_jg_j⁴<∞`, which closed, G1.1): that bounds a *sum*, the tower is a *min over shells* — different objects, and the closed one does not touch the open one. Crucially the collapse is a **bound artifact**, not geometric: the physical gap (`c_L=4π²/L²` box scale, `1/a_j²` high-pass band in the UV) is `a`-uniform and *grows* toward the UV, anti-correlated with the collapse. The one inequality freezing `Λ_UV,j` at its physical floor is exactly `(★)`, the one-sided variance upper bound. So the tower *should* be `a`-uniform; the degradation is precisely the missing `(★)∧(★dec)`.

## 1 · The object
Tower of conditionals `ν_j=P_{j*}\mu_{j+1}` down the `J∼⌊log₂(1/aΛ)⌋` shells, each block map `P_j` with `Lip(P_j)=C_blk=O(1)`, `a`-independent. L3's block-Poincaré leg (G1-UNIFORM-ZA §1 **(iii)**) is: `inf_{j≤J}Λ_UV,j ≥ δ>0` **uniformly in `a`**. Just existence per shell is `(ii)`, proven; the whole content is uniformity.

## 2 · The one proven constant collapses monotonically
PROVE-H2: fiber = bounded perturbation of Haar, `osc_ξ(β_jH_fib)≤2β_jN_p(B)`, so Holley–Stroock gives
```
Λ_UV,j = e^{−2β_jN_p(B)}/C_G ,   N_p(B)=O(1), uniform in U′,L, boundary.
```
Because `β_j∼4Nβ₀\,j` strictly increases with shell index, `Λ_UV,j` strictly **decreases**: the deepest shell `j=J` is the worst. With `J∝ln(1/aΛ)`, `e^{−2β_JN_p}=(aΛ)^{2·4Nβ₀N_p/ln2}=(aΛ)^{c}`. Over `J` shells this is not a benign product but a `min` collapse to `(aΛ)^c→0`. **This leg alone returns L3 false.** Note the scale-sum `Σ_jg_j⁴=(2β₀)⁻⁴ζ(4)<∞` (BREAK-IN-1) is irrelevant here — it bounds a *sum of coupling powers*, never the *worst-shell constant*.

## 3 · SZZ is the genuine `a`-uniform LSI — but off the tower's path
Lemma SC (PROOF-COARSE-POINCARE (c)): at `β<1/48` SZZ give a volume-, hence `a`-uniform, LSI and `Var_ν(φ)≤(C_blk²/λ_μ)∫|∇φ|²dν` — a true `a`-uniform anchor. But SZZ is a **high-temperature** (strong-coupling) theorem, `β<1/48≈0.02`, while the renormalised tower occupies `β_j∈[2.5,∞)`: the box scale sits at `β≈2.5≫1/48`, and asymptotic freedom drives `β_j→∞` into the UV. The `β`-windows **do not overlap**. "An SZZ-type bound at each effective coupling" is *directionally wrong*: SZZ's regime is where the flow starts (strong coupling), the tower flows the other way (weak coupling). The only tool living at every `β` is Holley–Stroock (§2), and *its* constant is the collapsing one.

## 4 · The physical constant is `a`-uniform — and `(★)` certifies it
The gap the tower *should* have is the free one: `c_L=4π²/L²` at the box scale (`a`-uniform at fixed `L`, P1-L0), and the high-pass band gap `≈2/a_j²≥c_L` deeper in the UV; both sit above the **`a`-independent** floor `c_L`. The retargeted single-block identity `G″_j=E[H_j]−Var_j+W_j≥c_L+δ_W−Var_j` (Lemma B `W≥δ_W` scale-independent) shows the physical floor is `c_L+δ_W>0` **provided `Var_j` is controlled** — from above, which is exactly
```
(★)  Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E⟨PᵀH_jP v,v⟩+W_j(v)−δ_W ,
```
uniformly in `j,U′,a`. So `(★)` is precisely the inequality that **replaces the collapsing HS constant by the freezing physical floor**, for every shell. HS falls with `β_j` though the true gap rises, because it overcounts the interaction as oscillation `2β_jN_p` in the deep UV — where the theory is asymptotically free and the gap is nearly the free one. That overcount is the entire degradation; `(★)` is its cancellation.

## 5 · Transport across shells is a second, separate leg
Even granting per-shell `(★)`, the tower constant is the `min` over shells, and `δ_W` must survive the off-diagonal score covariance: Gershgorin gives `inf_j[c_L+δ_W−Σ_{B′≠B}|Cov(s_B,s_{B′})|]`, finite only with band decay
```
(★dec)  |Cov_{μ(·|U′)}(s_B,s_{B′})| ≤ C e^{−κ\,dist(B,B′)/a},  κ>0.
```
This is open for the same sign reason (L3-band-decay): the band fixes the range `O(a)`, the rate is the Savvidy-indefinite shell Hessian — a *separate* leg from the per-shell degradation of §2–4.

## 6 · Ledger / honest boundary
| leg | status |
|---|---|
| per-shell existence `Λ_UV,j>0` (HS) | **proven** (collapses `∼(aΛ)^c`) |
| SZZ `a`-uniform LSI at `β<1/48` | **proven**, **off the RT path** |
| tower `inf_jΛ_UV,j≥δ` `a`-uniform | **open** = `(★)` |
| cross-shell band decay | **open** = `(★dec)` |
| physical floor `c_L+δ_W>0` | **plausible**, certified by `(★)` |

**Sharpest line.** The tower does not stay `a`-uniform with proven tools — the block Poincaré constant collapses monotonically to `(aΛ)^c` at the deepest shell, and that degradation *is* the `(iii)` obstruction of L3, not a geometric invariant: the free gap is `a`-uniform and rising, so the collapse is the Holley–Stroock overcount of a theory that is asymptotically free exactly where the collapse is worst. One-sided `(★)` plus `(★dec)` removes the degradation and makes the block-Poincaré leg of L3 `a`-uniform — still fixed-volume existence only, never the gap, never `L→∞`.

*Refs.* PROOF-COARSE-POINCARE (c); G1-L3 §2; G1-UNIFORM-ZA §1–3; PROOF-L3 §1–4; L3-scale-uniform §2; L3-band-decay; PROVE-H2-uvgap (Λ_UV(b)); G1-EXISTENCE-THEOREM (ledger); G1-fixed-volume; GAP-G1-largefield §1(i),(ii); BREAK-IN-1 (S1 scale-sum); P1-L0; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737 (β<1/48).
