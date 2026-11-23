# PROOF-L3 — attempting the one open UV input directly

> **Verdict (判词). NOT closed.** L3 is the **a-uniform per-shell marginal log-concavity** `G″_j≥δ` across the `J≈log(1/a)` tower. All three routes assemble the *lower*-bound terms (`E[H]≥free gap`, `W≥δ_W` scale-independent, `Λ_UV` provable only as `e^{−2β_jN_p}→0`) and all three die on the **same** missing inequality: an a-uniform **upper** bound on the force variance, `Var(G_v)≤E[H_vv]+W(v)−δ`. That is the one-sidedness the Savvidy-indefinite shell Hessian forbids. Fixed volume removes the IR/zero-mode residue, so L3 is genuinely UV and `⊥` gap — the only good news. G1 remains conditional on L3 alone; I removed nothing.

## 1. L3, precisely
Fix `L`. Along the RT (`g_j²∼1/(2b₀log(1/a_jΛ))`, `β_j=2N/g_j²`, shells `j=j*,…,J`, `J∼⌊log(1/aΛ)⌋`). With the exact identity `G″_j=E_{μ(·|U′)}[PᵀH_jP]−Var_{μ(·|U′)}(Pᵀ∇S_j)+W_j` (PROVE-Gpp), L3 is: the normalized generator is bounded, `sup_a Z_a^{-1}‖H_a‖_D<∞` (SPEC-RES′), **and** the renormalized moments are a-uniform (UVS(L)). The second half is the substance: `UVS(L)` moments scale with `N_□∼L⁴/a⁴→∞`, so operator-norm *finiteness* (compactness gives `O(N_□)`) is not a-uniform — a-uniform tails require log-concavity. Hence L3 ⟺ the **missing inequality**
```
(★)  Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E_{μ(·|U′)}⟨PᵀH_jP v,v⟩ + W_j(v) − δ   ∀|v|=1, ∀j, uniformly in U′,a,
```
plus, for the fixed-`L`, `a→0` count `N_{E′}∼1/a⁴`, the decay `|Cov(s_B,s_{B′})|≤Ce^{−κ|B−B′|/a}` (`★dec`).

## 2. Route (i) — one-sided Savvidy sign
Cone-empty (GAP-G1(i)) gives `H_vv⪰−2g_j|F_j|≥−O(g_j²)`, i.e. `H⪰0` on the typical region: a **lower** edge. But `G″≤E[H]` because `Var≥0` *subtracts* (MARGINAL, EDGE-f₁²). So `H⪰0 ⇒ E[H]≥0` says **nothing** about `E[H]−Var`; the variance can drive it arbitrarily negative, and the sign we get is the opposite of what `G″≥δ` needs. It also cannot bound the upper edge (`H_vv≤β_jC_N#□∼β_j`), which the norm requires. **Dies: wrong edge.**

## 3. Route (ii) — Lemma B + free-box-gap (GL3/GL4)
`W≥δ_W` (SU2 ½g, SU3 ¾g, *scale- and β-independent*) and the free box gap `c_L=4π²/L²` (a-uniform at fixed `L`) give the two **positive** terms. `G″≥c_L+δ_W−Var`. Missing: an upper bound on `Var(G_v)`. The only generic bound is `Var(G_v)≤E[G_v²]∼β_j²` (Cauchy–Schwarz), off from the needed `E[H]∼β_j` by a **factor `β_j`**; the bridge is Gaussian equipartition, exact only as `g_j→0`. The rigorous `Λ_UV` (Holley–Stroock, PROVE-H2) is `e^{−2β_jN_p}`: at `β_j∼2.6` (`a=10⁻³`) this is `10^{−107}` — catastrophically **non-uniform**, not the physical free gap `c/a_j²`. **Dies: no `Var` upper bound; HS gap `→0`.**

## 4. Route (iii) — band-limited Stein/score
The identity `Cov(F,s_B)=−E⟨∇_⊥F,∇_⊥g_B⟩`, `‖∇_⊥g_B‖²≤Var(s_B)/Λ_UV`, mass `m=min(λ/2,Λ_UV/(1+2S/λ))` is **correct** (REPAIR-2SCALE, VERIFY-FINAL(2)). But `m` is a-uniform ⟺ `Λ_UV` is a-uniform (since `Var(s_B)∼β_j²→∞`), i.e. **exactly (★)**. With the only proven `Λ_UV=e^{−2β_jN_p}`, `m→0` doubly-exponentially. The "band limit `|k|≥π/2a`" fixes the *range* `O(a)` and hence the volume count `N_{E′}∼1/a⁴` — it does **not** touch the *scale* factor `e^{2β_jN_p}`. **Dies: scale factor, not volume factor.**

## 5. Sharpest line / what closing buys
The deep-UV shell is ordered (`β_j→∞`), so `|F_j|=O(g_j)`, the Savvidy cone `√(2g_j|F_j|)∼g_j` lies below `Λ_{j+1}`, and the physical `G″_j≈W≥δ_W` — L3 is **true**. But the bridge from `Var(G_v)≤E[G_v²]∼β_j²` to `Var(G_v)≲E[H]∼β_j` is strict one-loop exactness, *asymptotic* (renormalon, BREAK-IN-7); Balaban small-field/large-field RG is the right frame yet reduces to this bootstrap (gauge-fixed, incomplete). No method signs the non-Gaussian remainder `Var−E[H]` because the shell measure is not log-concave. **If (★)∧(★dec) held, G1's ledger loses its sole open line and fixed-volume existence (OS0/1/3, subsequential) becomes unconditional.** It does not: `G1 ⟸ L3` still.

## Refs
A2 §L3; G1-EXISTENCE; G1-fixed-volume §2–3; L3-scale-uniform; L3-band-decay; PROVE-Gpp §4; A1-twosided-UV §4; GAP-G1 §1,§3; PROVE-H2; H2-b-uniform §1–3; REPAIR-2SCALE; VERIFY-FINAL(2); W6-BE; W7-RERUN; P1-L0; BREAK-IN-7; EDGE-f₁².
