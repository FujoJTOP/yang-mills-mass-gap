# POTENTIAL-THEORY-CAPACITY — species CAP: Wiener criterion / equilibrium capacity on the coarse zero-mode marginal. The sign is *free* (capacity is a min of a manifestly positive Dirichlet energy — no indefinite Hessian), and the free-field acid test is genuinely non-vacuous (flat zero mode ⇒ infinite 2nd moment ⇒ zero capacity). But finite volume makes every sign bound vacuous (compact support), and the non-vacuous *uniform* content — a capacitary conductance floor — is H1's tail-coercivity in flow costume.

> **判词 (verdict).** Forge: treat the coarse `k→0` mode `a` as the 1D diffusion `da=−V̄′(a)dt+√2dB`, marginal `ν(da)=e^{−V̄}da/Z`, and bound `Var_ν(a)` by **potential theory** — the Wiener 2nd-moment test, the equilibrium capacity `Cap_ν(A,Aᶜ)=inf{∫|∇f|²dν:f|_A=1,f|_{Aᶜ}=0}`, and the capacitary conductance `Φ_eq=inf_A Cap/ν(A)`. Target `G″(k→0)≥λ ⟺ Var(a)≤1/λ`. **Verdict: a new *geometry* and the cleanest non-vacuous acid test in the family — but 2/5, not 5/5.** The exact criterion `Var(a)<∞ ⟺ ∫(1+a²)e^{−V̄}da<∞` is a **tail-coercivity** condition, **strictly weaker than pointwise `G″≥λ`** (a flat bulk with a `>2log` tail passes). It fails the free field *for real* (free massless `V̄≡0` ⇒ divergent series ⇒ `Var=∞`), unlike B3, whose fiber LSI was blind to the free transverse mode. But finite volume (compact SU(N)) makes `Var<∞,λ_1>0,Cap>0` **trivial**, so #2 is vacuous; and the uniform floor `Φ_eq≥c` — the Savvidy dip as an Arrhenius barrier — is H1. **#1,#4 pass; #2(vacuous),#3,#5 fail.**

## 1 Forge and dictionary
Transform to the `V̄`-Langevin form: scale `s′=e^{V̄}`, speed `m∝e^{−V̄}dx`, `V̄=−log(pushforward of μ)`; on the zero mode `G″≈V̄″` (COMPUTE-GPP; AQFT-CAUSAL). The lever is **capacity**, not curvature: `Cap` is `inf` of the **positive** form `∫|∇f|²dν⩾0`, a **convex** Dirichlet-min — whereas `G″⪰λ` is a **min-eigenvalue of the indefinite** `PᵀHP` (min ≈−0.3, MARGINAL). **This trades "sign an indefinite operator" for "min a positive energy."**

## 2 Sign is free — and vacuous in finite volume
`λ_1≥Φ_eq²/4` (conductance ⇒ gap; Sinclair 1992, Jerrum–Sinclair 1989). The RHS is manifestly `⩾0`: the sign is **automatic**, no Savvidy contest. But compact SU(N) (finite `L`) ⇒ `ν` compactly supported ⇒ `Var(a)⩽diam²<∞`, `λ_1>0`, `Cap>0` **for any measure** — so the sign bound holds **vacuously**. All content is `L`-uniformity (#5), exactly as B3's `#4` was vacuous.

## 3 The wall, as an Arrhenius barrier — mild at the target scale
The max-flow dual `Cap≥1/E(θ)`, `E(θ)=∫θ²/ν`, makes `Cap` **constructive**; but high `V̄` (small `e^{−V̄}`) is an **expensive** flow — the Savvidy dip (`V̄″<0`, saddle) is a barrier, `λ_1≳e^{−2H}` (Eyring–Kramers). **But** at the target thermal shell `|a|⩽λ^{−1/2}` the measured dip gives `H∼|V̄″|·λ^{−1}≈0.3/λ`, which is `O(1)` for `λ=O(1)` — **no exponential suppression at the target scale.** The obstruction is thus a **bootstrap on λ** needing an *a priori* seed `λ⩾c`, not a wall.

## 4 Acid test — non-vacuous (the star)
Free massless (scaling limit): the zero mode is **flat/non-compact**, `V̄≡0`, `∫(1+a²)e^{−V̄}=∞`, the Wiener series over `{|a|≍2ᵏ}` **diverges**, `Cap=0` (thin/polar, Kakutani), `Var=∞`. YM confined: finite. **The tool fires on the free field (returns no bound) and only on it** — a **geometric** discriminator (flat-direction moment), distinct from OTHER-H1's Haar term and B3's blindness. No SUSY.

## Exact criterion
```
Var_ν(a)<∞  ⟺  ∫(1+a²)e^{−V̄(a)}da < ∞  ⟺  Σ_k 2ᵏ·(dyadic |a|≍2ᵏ mass) converges
            ⟺  2-capacity of {|a|>R} summable in R   (Wiener 2nd-moment test)
uniform Var(a)≤1/λ  ⟺  V̄ coerces ½λa² at the tail (asymptotic, NOT pointwise V̄″≥λ)
λ_1 ≥ Φ_eq²/4,  Φ_eq = inf_A Cap_ν(A,Aᶜ)/ν(A),  Cap ≥ 1/E(θ)   [max-flow; convex, positive]
Savvidy dip  =  high-V̄ saddle  =  small e^{−V̄}  =  costly flow  =  H∼0.3/λ  =  O(1) at target
free massless: V̄≡0 ⇒ series diverges ⇒ Cap=0 ⇒ Var=∞   [acid test met, non-vacuously]
```

## 5 Five-point scorecard
1. **Nonperturbative — PASS.** Wiener/Feller test, Dirichlet-principle capacity, max-flow: exact, all β, no loops.
2. **Sign-controlling LOWER bound — FAIL (vacuous).** Sign is free (`λ_1≥Φ_eq²/4⩾0`) but trivial on compact finite volume; useful sign is #5.
3. **Non-circular — FAIL at closing (relaxed).** Functional of the marginal `ν`; B3's two-coordinate obstruction blocks fiber supply. The uniform input `Φ_eq≥c` **is** H1. But it demands only **asymptotic tail-coercivity** (a flow/occupation bound), **strictly weaker** than pointwise `G″≥λ` — a relaxation, not a rename.
4. **Non-SUSY + fails free field — PASS (cleanest, non-vacuous).** Flat free mode ⇒ zero capacity; confined YM ⇒ positive. Genuinely discriminating geometry.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** Only a bootstrap on `λ` (barrier `O(1)` at target scale) needing an a-priori seed; no `L`-uniform proof.

**Sharpest line.** Capacity replaces the wall's hard part — *signing an indefinite Hessian* — with a **convex, positive, constructive** Dirichlet min, and earns the only **non-vacuous** free-field discriminator in the species family (the flat zero mode has literally zero capacity, not a blind fiber term). What it does **not** buy: on compact finite volume positivity is trivial, and the `L`-uniform floor `Φ_eq≥c` is exactly H1 — in the costume of "no large Arrhenius barrier on the thermal shell," a **milder-looking** (barrier `O(1)`, not exponential) yet still **unclosed** estimate. New shape of the key; same lock.

*Refs.* B3-ISOPERIMETRIC (Cheeger, two-coordinate obstruction, acid-test-blindness); PROOF-COARSE-POINCARE (Lemma SC, `PᵀHP` min, `G″` exact); OTHER-H1-POINCARE (Haar-term discriminator); C8-SPECIES-NOGO (empty-referent pattern); MARGINAL (`PᵀHP` −0.3); B4-VARIANCE-MAG (`D/S≈2%`); PROOF-JACOBIAN-SAVES (`G″=+0.51`). External: Wiener, *J. Math. Phys.* **3** (1924) 127; Feller, *Ann. of Math.* **55** (1952) 468; Kakutani, *Proc. Imp. Acad. Tokyo* **20** (1944) 706; Doyle–Snell, *Random Walks and Electric Networks*, MAA Carus **22** (1984); Sinclair, *Combin. Probab. Comput.* **1** (1992) 351; Jerrum–Sinclair, *SIAM J. Comput.* **18** (1989) 1149; Khas'minskii, *Theory Probab. Appl.* **5** (1960) 179; Eyring, *J. Chem. Phys.* **3** (1935) 107; Kramers, *Physica* **7** (1940) 284; Cheeger (1970); Bakry–Émery (1985); Evans–Gariepy, *Measure Theory and Fine Properties of Functions* (capacity, Wiener criterion).
