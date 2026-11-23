# D-BBS-LOWRANK-KATO — low-rank negative subspace + Kato/Weyl transfer + momentum separation

> 2026-11-23 · forging the *scale-uniform* `G″_j≥δ` along the proposed split `H=H_++H_−`. **判词 (verdict): the split is *exact* and the *UV tower closes scale-uniformly*, but the transfer theorem is **circular exactly at `k→0`** — the one shell where the gap must be printed.** (a) `H_−` is the `n=0` Landau level of the charged transverse pair: rank `r_−∝g²|F|²𝓥` (∝`∫|F|²`), depth `g|F|`, all inside `(−g|B|,0)`; every shell with `λ_j²>c_F‖F‖_∞` is **cone-free with margin `λ_j²−c_F‖F‖_∞`** (Weyl). (b) Kato's relative constant is `a_j=c_F‖F‖_∞/λ_j²`, `<1` **iff** `λ_j²>c_F‖F‖_∞` — so positivity transfers cleanly on the whole UV tower, rate-uniform. (c) **Momentum separation fails**: the Savvidy band *bottoms at `k_∥=0`* (`λ_{0,−}=k_∥²−gB`), so the negative subspace is **softest at zero momentum** — the same regime as the coarse mode — and the measured `PᵀHP` min is `≈−0.3<0` (FORGE5 K2). Overlap `ε=‖P_−v‖²` is `O(1)`, not small. The residual is `G″`'s variance term at `k→0` — **H1**.

## (a) Rank/depth per shell — exact, scale-uniform
On a constant abelian `F₁₂=B`, the transverse Hessian is `H=−D²+2adF`, ladder `λ_{n,s}=k_∥²+(2n+1)g|B|+2s\,gB` (PROOF-SAVVIDY-EXACT-SIGN). Only the `s=−` `n=0` rung dips: `λ_{0,−}=k_∥²−g|B|`, negative on `k_∥²<g|B|`, bottom `−g|B|`. Rank (`d=4`, `V=2adF`): per charged pair `r_−=g²B²𝓥/8π²` (∝`∫|F|²`); weight `Tr H_−≤8C_N∫|F|³` (Lieb–Thirring). Since `H⪰−D²−c_F‖F‖_∞` (`c_F:=‖2adF‖/‖F‖_∞=O(1)`), on any Littlewood–Paley shell with `P_j(−D²)P_j⪰λ_j²`,
```
P_j H P_j ⪰ (λ_j² − c_F‖F‖_∞) P_j .            (Weyl)
```
Thus **every shell above `k_⋆:=√(c_F‖F‖_∞)` is positive-definite with a margin that *grows* as `λ_j²`** — the UV tower is cone-free, uniformly, and the negative modes are confined to the `O(1)` shells with `λ_j<k_⋆`. `r_−/dim≈10⁻³` at β=2.5 (count-small), but *rank-smallness is a statement about counts*, whereas the gap is the **min** — where one negative eigenvalue is decisive. Count `E[neg]/E[pos]≈0.94` (trace-half-positivity `Tr H=4Tr(−D²)≥0`) already records this: the negative *fraction* is small, the *bottom* is not.

## (b) Kato/Weyl transfer — clean on the UV tower
`B=2adF` is `A`-bounded (`A:=−D²`): `‖Bψ‖≤c_F‖F‖_∞‖ψ‖≤a_j‖Aψ‖` with `a_j=c_F‖F‖_∞/λ_j²<1` iff `λ_j²>k_⋆²`. Then `A+B⪰(1−a_j)A⪰(1−a_j)λ_j²>0` (Reed–Simon IV, Thm X.12; Kato 1966). This is the intended transfer and it is **correct and scale-uniform** on the UV half: relative bound `a_j→0` toward the UV, margin `λ_j²−k_⋆²` strictly positive and rising. No small-rank hypothesis is even needed here — Weyl suffices, and the low-rank structure only sharpens the *count* of bad directions.

## (c) Momentum separation fails — the band is soft at `k=0`
The premise "coarse low-momentum ⊥ negatives high-momentum" is **false at the operator bottom**: `λ_{0,−}=k_∥²−gB` is *minimized at* `k_∥=0`. The negatives have 3-4 momentum `k_∥≤√(gB)` but *include* `k=0`; in 1-2 they are Landau orbitals of width `ℓ=1/√(gB)`. The `k→0` coarse mode and the deepest negative share `k_∥=0`, so overlap is not suppressed. The spectral split gives, with `ε=‖P_−v‖²`,
```
⟨v,Hv⟩ ≥ h_+‖v_+‖² − g|B|‖v_−‖² ,  h_+=bottom of positive spectrum ,
```
need `h_+(1−ε) − g|B|ε ≥ δ`. Two facts kill it at `k→0`: (i) `h_+` is the **free zero-mode value**, `→0` (no `λ_j²` floor at `k=0`), so the LHS `→−g|B|ε<0` for **any** `ε>0`; (ii) measured `PᵀHP` min `≈−0.3`, i.e. `ε=O(1)`. With the *shell* floor instead, the requirement is `ε≤(λ_j²−k_⋆²−δ)/(λ_j²−k_⋆²+g|B|)`, whose numerator is negative once `λ_j²<k_⋆²+δ` — no `ε∈[0,1]` helps. **The coarse mode sits on the expensive side** (FORGE5 K1–K2; `k_⋆≳k_max` at β=2.5).

## What survives, acid test, and the one claw
- **Survives (true, scale-uniform):** the UV-tower theorem `λ_min(P_jHP_j)≥λ_j²−k_⋆²>0` for `λ_j>k_⋆`; the rank bound `r_−≤C_N∫|F|²`; trace-half-positivity `Tr H_−≤Tr H_+≤8C_N∫|F|³`. These are the *UV half* of BBS, done.
- **Acid test:** free field `F≡0⇒k_⋆=0⇒` every shell trivially positive ⇒ the theorem holds on the massless field. So as a **gap** tool it is **vacuous** (it would "prove" a gap for a free massless field). It is a *coercivity* statement, not a sign-of-gap statement — it fails the free-field discriminator that a gap-closer must pass.
- **One claw (det leg, not gap):** because the BBS marginal-log resummation is a **Fredholm determinant** (C-BBS-MARGINALLOG (A)–(B)), `log det H = Tr log H` is *exactly* where low-rank pays: `Σ_{i∈neg}log|λ_i|` has at most `r_−` terms each `≤log(g|B|)`, so `Tr log H_−` is computable and trace-bounded. This closes the **one-loop/marginal-log leg**, not the gap.

## Exact residual
`G″=E[H]−Var(G)+W`. Kato+Weyl+low-rank bound **only** `E[H]` and only where `λ_j²>k_⋆²`; `−Var` (the `e^{−S}`-weighted fluctuation, `≈0.99` of the bare margin, PROOF-COARSE-POINCARE) and the `k→0` cone are untouched. The transfer is **circular at `k→0`**: it needs `h_+>0` there, i.e. an induced mass — which *is* the gap. Precisely missing: `∃ε<δ/(δ+g|B|)` with `PᵀHP|_{\rm coarse}≥h_+>0`, refuted at β=2.5 by `−0.3`.

**Sharpest line.** Low-rank+Kato genuinely *empties the UV cone* at rate `λ_j²−k_⋆²` — a real, scale-uniform half-theorem — but the Savvidy band is *soft at `k=0`*, so "high-momentum negatives" is a misread: the negatives and the coarse mode share the soft bottom, `ε=O(1)`, and at `λ_j→0` the floor `h_+→0` makes the transfer sign-blind (and free-field-positive). Low rank buys a determinant leg, not the min. The residual is `Var−W` at the cone — **H1**.

*Refs.* PROOF-SAVVIDY-EXACT-SIGN (spectrum, `r_−∝g²B²𝓥`, depth `−g|B|`, `Tr H_−≤8C_N∫|F|³`); FORGE5-SL-UV-CONE-EMPTY (K1 sup, K2 `k_⋆≳k_max`, `PᵀHP` min `−0.3`); W6-BE-at-2.5 (min `−1.1β`, negatives 90–97% fluctuation block); PROOF-COARSE-POINCARE (`G″=E[H]−Var`, H1, margin 99% cancelled); C-BBS-MARGINALLOG (A/B port, Fredholm det); P2-BBS-to-gauge; CONSTRUCTION-METHODS §verdict; Kato, *Perturbation Theory for Linear Operators* (1966); Reed–Simon IV, Thm X.12; Weyl (1912); Lieb–Thirring (1976); Ky Fan/Cauchy interlacing.
