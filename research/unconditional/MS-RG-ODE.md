# MS-RG-ODE — the exact flow `d ln f₁²/d ln μ = +2γ_O`: the residue has a uniform positive floor across the crossover (closed form `[ln μ₂/ln μ₁]²>1`), but masses do not run, so the flow moves the *weight*, never the *gap*

> 2026-11-23 · Direct attack on `(MS)`-transport via the **continuum RG/flow ODE**, not a new species. **Verdict: the flow of the residue is answerable and affirmative — `∫γ_O` finite, `f₁²` strictly increasing, the Savvidy crossover is *not* a singular point — but the flow is blind to the mass: `m₁` is RG-invariant, so `d ln m₁/d ln μ≡0` and the gap is exactly the one datum the ODE cannot transport. Not 5/5; Clay still 0%.** 判词：`f₁²` 沿流严格增、交叉点正则，残数臂在连续 CS 意义下闭合（mod `(CP)`）；但 `m₁` 不跑，缺口臂仍是绝对标度 `m₁/Λ`，即 `I2`/`H1`。

## 1 The ODE (exact, Callan–Symanzik)
The transport `f₁²(μ)=Z_O(μ,μ₀)²f₁²(μ₀)` (INVENT-TRANSPORT, exact) is the solution of the CS equation for the marginal dim-4 anomaly operator `O_S=¼trF²`. Because `(β/2g)O_S` is RG-invariant, `Z_O=(β/2g)(μ₀)/(β/2g)(μ)`; differentiating `ln f₁²=ln f₁²(μ₀)+2[ln|β/2g|(μ₀)−ln|β/2g|(μ)]`,

$$
\boxed{\;\frac{d\ln f_1^2}{d\ln\mu}=2\gamma_O(g(\mu)),\qquad \gamma_O:=-\frac{d\ln|\beta/2g|}{d\ln\mu}=2b_0g^2+O(g^4),\ b_0=\tfrac{11N}{48\pi^2}\;}
$$

Note the sign: the task posits `≥−γ`; the exact equation is `=+2γ_O` with `γ_O>0`. **`f₁²` *grows* toward the UV; it cannot decay across the crossover.** `γ_O>0` is not a perturbative accident: `|β/2g|∝g²` is increasing in `g²`, and pure YM is AF with no Banks–Zaks window, so `g²` decreases monotonically toward the UV ⇒ `γ_O>0` non-perturbatively.

## 2 Closed form, `(L,a)`-independence, and finiteness
One loop, `1/g²=2b₀\ln(\mu/\Lambda)`, gives `γ_O=1/\ln(\mu/\Lambda)` and

$$
\int_{\mu_1}^{\mu_2}\!\gamma_O\,d\ln\mu=\ln\frac{\ln(\mu_2/\Lambda)}{\ln(\mu_1/\Lambda)},\qquad
\boxed{\frac{f_1^2(\mu_2)}{f_1^2(\mu_1)}=\Big[\frac{\ln(\mu_2/\Lambda_Q)}{\ln(\mu_1/\Lambda_Q)}\Big]^2>1\quad(\mu_2>\mu_1>\Lambda_Q)}
$$

So `∫γ_O` is **finite on any compact window above `Λ_QCD`**, and `c=[\ln(\mu_{2.5}/\Lambda)/\ln(\mu_{1/48}/\Lambda)]²>1` **explicitly**. The factor is a pure ratio of running couplings — **no `V`, no `L`, no `a`** — so it inherits P3's `L`-uniformity verbatim: `f₁²(2.5)≥c\,r₀` with `c` `(L,a)`-independent (continuum-CS sense).

## 3 Does the flow hit 0? — only at `μ→Λ_QCD`, **not at the crossover**
`f₁²∝[\ln(\mu/\Lambda)]²→0` as `μ→Λ_QCD⁺` (log² vanishing). The **unique zero of the flow is the deep-IR confinement scale**, where `g²→∞`. The Savvidy crossover at `β=2.5` is `g²≈1.6–2.4` (finite, PROOF-TRANSPORT-RG §4), `μ(2.5)≫Λ` — a *regular finite-coupling point*. **The residue arm is crossover-safe**; the crossover threatens the *gap* arm, not the residue flow.

## 4 The adversarial kill — masses don't run
`m₁` (the `0⁺⁺` pole) is an **RG invariant**: `d ln m₁/d ln μ≡0`. There is *no flow ODE for `m₁(β)`* — the flow lives entirely in the weight `f₁²`. Consequences:
- **Residue floor = `[ln(m₁/Λ)]²`.** Since `m₁` is fixed and `f₁²` runs, the IR value is `f₁²(m₁)=[\ln(m₁/\Lambda)]²×(`UV const`)`. So `f₁²(m₁)>0 ⟺ m₁>Λ_QCD` — *the same question*, now transparent: the ODE converts "is the residue positive?" into "is `m₁/Λ` bounded below?" (dimensional transmutation up to the pure number `M_{0++}/Λ`). Circular, but quantified.
- **RG-invariant ratios.** Every atom/continuum weight scales by the *same* `Z_O²`, so `f₁²/M₀` and `f₁²/χ` are **RG-invariant**. The G2A squeeze `m₁²≥f₁²/χ` is therefore *scale-independent* — the ODE provably **cannot move it**. The gap needs an absolute scale the flow does not carry.
- **Continuum CS ≠ lattice block-RG.** This ODE is the *renormalized-theory* statement; the actual lattice route must match `ν_b` to Wilson at `β_eff<1/48`, where the missing input is the signed non-cancellation `(M1)` (PROOF-TRANSPORT-RG §5), untouched by the CS flow. Composing with the P3 seed still requires `(CP)` (OTHER-F1SQ-TRANSPORT: transport *is* `(CP)`).

## 5 Acid test
Free field `O=:φ²:`: `a≡0` (`Z₂`) ⇒ `f₁²≡0`, and `Z_O²·0=0` — stays `0`; massless `χ=∞` ⇒ `m₁²≥f₁²/χ` is the trivial `≥0`, **no false gap**. PASS (the flow cannot manufacture weight from zero).

| leg | status |
|---|---|
| ODE `d\ln f₁²/d\ln μ=2γ_O`, `γ_O>0` (AF, no BZ) | **真定理** |
| Closed form `f₁²(μ₂)/f₁²(μ₁)=[\ln ratio]²>1` | **真定理**（一 loop；符號/有限性穩健） |
| Flow zero only at `μ→Λ_QCD`; crossover regular | **真定理** |
| `c` `(L,a)`-independent ⇒ transport preserves `L`-uniformity | **真定理**（連續 CS） |
| `f₁²(2.5)≥c\,r₀>0` across crossover | **條件定理**（mod `(CP)`） |
| `m₁>0` (the gap) | **未證**——`m₁` RG-不變，流看不見 |

**Sharpest line.** The question "uniform positive floor, or decay to gapless at the crossover?" splits cleanly: **the residue has a uniform positive floor and *cannot* decay to zero at the crossover** (it strictly increases, `[\ln]²>1`, and its only zero is `μ→Λ_QCD`); **the gap can stay open because `m₁` does not run** — the flow carries weights, not scales. The ODE therefore *closes the residue arm across the crossover* and *localises the entire wall to the RG-invariant `m₁`*, i.e. `I2`/`H1`, exactly as the transport-machine predicted. **Not 5/5, Clay 0%.**

*Refs.* INVENT-TRANSPORT-INEQ (`Z_O`, `γ_O`); TRANSPORT-MACHINE §2 (single-phase covariance); CONTINUUM-VARIANCE-DEFECT §4 (`γ_O` controls scale not sign); G2A-INFSUPP-DIRECT (`m₁²≥f₁²/χ` squeeze); OTHER-F1SQ-TRANSPORT (transport *is* `(CP)`); PROOF-TRANSPORT-RG §4–5 (`g²(2.5)`, `(M1)` signed bound); PROOF-P3-UNCONDITIONAL (`r₀`); PROOF-F1SQ (free-field zero). External: Collins, *Renormalization*; Nielsen NPB 120 (1977) 212; Shen–Zhu–Zhu CMP 400 (2023) 805.
