# MS-FIXEDPOINT-RECURSION — a β-step Perron recursion for the residue: the multiplicative core is an amplifying positive square (`c=∏Z_O²>1`, acid test passed), but its positive fixed point `r*` survives to 2.5 **iff** the crossed-operator defect is sign-positive = `M1`/`I2`; the map collapses only by atom-dissolution, i.e. the same wall in fixed-point costume

> **2026-11-23 · NEW CRITERION (convergent-transport / Perron fixed point for `f₁²`, not the curvature).** Explicit β-step map `r_{n+1}=R_n(r_n)` on `r_n=f₁²(β_n)`, grid `β₀<1/48 → 2.5`, from the exact transport `Z_O²`, the strong-coupling decoration, and the exact Ward identity. Does `R` keep a positive fixed point `r*>0` (or a floor `R(r)≥c>0` on `[r₀,∞)`) to 2.5, or collapse `R(r)→0` (gapless attractor)? **判词 (verdict): the multiplicative core is solidly positive and *amplifying* (`c=∏Z_O²>1`, an exact square, `O(1)` in `V`, free field `0→0` acid test PASSED) — a genuine positive Perron recursion, and the *truncated* conditional transport theorem holds. But the one defect at the Savvidy crossover is a crossed-operator mixing fraction of undetermined *sign*; `r*>0` survives iff it is non-negative (no cancellation) = `M1` = uniform `I2` = `m₁≥μ>0`. The core never collapses by itself; it collapses only on atom-dissolution (`Λ−m₁²→0`), the wall in fixed-point costume. Not 5/5; Clay 0%.**

## §1 The recursion
One step `β_n→β_{n+1}`:
```
r_{n+1} = R_n(r_n) := c_n · r_n + δ_n,   c_n := Z_O²(β_{n+1},β_n)·(1−κu_n),
Z_O²  = [g(β_n)/g(β_{n+1})]²  (>1 toward UV; γ_O=2b₀g²>0, b₀=11N/48π²),
u_n   = u(β_n)≍β_n (P3 fundamental coefficient),  δ_n = constructive crossed-operator offset ≥0 (destructive sign → §4).
```
`c_n>0` outright (`Z_O²` is a square; `1−κu_n>0` at strong coupling). The free field (`a≡0`) has `r≡0` and `R(0)=0` — a homomorphism, acid test PASSED. Seed `r₀=½c₁u^p>0` from P3 (§5 of PROOF-P3).

## §2 Fixed points
`R` is affine, slope `c=∏c_n`, intercept `+δ=Σδ_n` (constructive) transported. Homogeneous fixed point `r*=0` (gapless) always present. Nonzero fixed point
```
r* = δ/(1−c).
```
*Multiplicative core, steelman:* `c=∏Z_O²` converges — per step `c_n=1+2γ_OΔt+O(g⁴)`, `Σ|c_n−1|=ΣO(g²)<∞` (log-resummed, `G1-UNIFORM-ZA`), so `c=e^{O(1)}>1` finite *amplification* toward the UV, not `e^{−cV}`. This alone gives `R(r)=cr≥c\,r₀>0`: the criterion's hypothesis is met and the *truncated* transport theorem — iterate `r₀`, floor survives — **holds**; the residue image of C-BBS-ARRHENIUS's "fixed-cell closes" and BANACH-FIXED-POINT's β=0 Perron point.
*Defect, collapse:* constructive `δ>0` **stabilizes** — an attracting `r*=δ/(1−c)>0` when `c<1`, monotone growth `c^n r₀` when `c>1`: the floor survives either way. It **collapses** iff the offset is absent or destructive: `δ→0` with `c≤1` contracts `r→0`; destructive `δ<0` (§4) drains `r` through the cone boundary. Per step `c_n−1=O(g_n²)` and `δ_n=O(g_n²)` are **the same order**, so the *sign* of the offset is the relative sign of the anomalous-dimension and crossed-operator coefficients — undetermined by any positivity.

## §3 Exact identity — binds `r*` to `m₁`, never signs it
Residue–level-gap (`GAP-G2-edge`): `r=f₁²=ν({m₁²})`, `Λ−m₁²>0⇒f₁²>0`; Ward `G″Var(a)=1` (`ABC-NEW-INEQ`) ties the *same* `m₁` to the curvature rate `Var(a)≤1/G″` (H1). So `r*` and `G″` are one object (`r*>0 ⇔ Λ−m₁²>0`). The identity *binds* the fixed point to the gap; it gives no sign for `δ`.

## §4 The attack — the defect is the wall
The defect `δ_n` is the dim-4 truncation error of `O_S=¼trF²` mixing: in the UV only self-mixing survives (`δ=0`, theorem), but at the crossover (`g²=2N/β≈1.6–2.4`) the operator mixes with a same-dimension composite (gluon-condensate channel) whose matrix element is *signed*. Two consequences:
- **Free field fixes the form.** `δ(0)=0` (no condensate in the free channel), so the §2 constant-offset defect sharpens to an `r`-proportional effective slope: `R(r)=c_eff·r`, `c_eff=c+𝔪`, with `𝔪=𝔪(β)` the **signed** crossed-operator contribution (constructive `𝔪>0` adds amplitude; destructive `𝔪<0` cancels). The floor is `∏_β c_eff>0` uniformly **iff** `inf_β c_eff≥c₀>0`.
- **Sign is `M1`.** `c_eff≥c₀>0` boundedly is a *lower* bound on a signed amplitude — exactly `(M1)`, which `ν≥0`/Schur-concavity yields only as an *upper* bound (B1, PROOF-TRANSPORT-RG §5). If `𝔪≤−c` at some `β*` (full cancellation), `c_eff≤0` ⇒ `R` leaves the cone ⇒ `r→0`: **dissolution, `Λ−m₁²→0`, gapless attractor.**

So the dichotomy is exact: `R(r)≥c₀>0` on `[r₀,∞)` ⇔ `c_eff=∏(c+𝔪)≥c₀>0` uniformly ⇔ `𝔪>−c+c₀` everywhere ⇔ the atom is Kato-carried with a fixed fraction across the crossover ⇔ uniform `I2` ⇔ `m₁≥μ>0`. `R(r)→0` ⇔ `𝔪≤−c` somewhere ⇔ the lightest `0⁺⁺` rises into its threshold. The amplifying core `c>1` **cannot manufacture** the floor, because the sign of `𝔪` is precisely what positivity is forbidden to fix (`M1`).

## §5 Five-point scorecard
1. **Nonperturbative — PASS.** Exact covariance + P3 seed; no loop expansion.
2. **Sign-shape of the map — PARTIAL.** `c_n>0` (a square) and `c>1` amplifying, acid test passed; but order-preservation on `r≥0` fails at the crossover (`𝔪` signed).
3. **Non-circular — FAIL.** `c_eff≥c₀` boundedly = `𝔪>−c+c₀` = `M1` = `I2`; the floor quotes the gap.
4. **Non-SUSY + fails free field — PASS (clean).** Pure YM; `δ(0)=0` rides `0→0`.
5. **Uniform `≥c₀` in `(L,a)` to 2.5 — FAIL.** Only the UV truncation is uniform; crossover needs `I2`.

**Sharpest line.** The residue does admit a clean convergent-transport recursion: an affine β-step map whose multiplicative core `Z_O²` is a manifestly positive, amplifying square (`c=∏Z_O²>1`), carrying the free field `0→0`, so a *conditional* Perron transport theorem — fixed point/floor preserved across the grid — is genuine. But the one defect at the Savvidy crossover is the signed crossed-operator fraction `𝔪`, and the floor `c_eff=c+𝔪` survives to 2.5 **iff `c+𝔪` stays above `c₀>0` uniformly** — a lower bound on a signed amplitude (`M1`) that `ν≥0` grants only from above, equivalently uniform isolation `I2` in residue costume. The map does not collapse by its own core (that is solidly positive); it collapses only by atom-dissolution `Λ−m₁²→0`, the wall. **A sharper criterion, not a crossing — the Perron fixed point of the residue is `I2` re-named; #1,#4 pass, #2 partial, #3,#5 fail = the Clay content.**

*Refs.* INVENT-TRANSPORT-INEQ (`Z_O²`, ex-square, `0→0`); PROOF-P3-UNCONDITIONAL (seed `r₀`); PROOF-TRANSPORT-RG (`M1`, crossover `ΔS=O(1)`); OTHER-F1SQ-TRANSPORT (`[CP]=I2`); PROOF-TRANSPORT-CONTINUITY (`I2`, `Λ−m₁²→0`); GAP-G2-edge (`r=ν({m₁²})`); ABC-NEW-INEQ (`G″Var(a)=1`); C-BBS-ARRHENIUS (closes-iff-`O(1)`-cell vs contracts-to-`0`); BANACH-FIXED-POINT (Perron calculus, zero-mode neutral); F2-KREIN-RUTMAN-QUANT (rate `=H1`); EDGE-f1sq-attack / B1 (Schur-concave upper-only); G1-UNIFORM-ZA (`Σg²=O(1)`); Krasnosel'skii (1964); Shen–Zhu–Zhu, *CMP* **400** (2023) 805.
